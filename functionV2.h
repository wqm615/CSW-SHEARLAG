//变高度梁相关函数文件
#include <set>
#include <sstream>
#include "StructDef.h"

extern vector<CContLoad> ContLoad;
extern vector<CBd> Bd;
extern vector<CUnifLoad> UnifLoad;
extern vector<double> Spans;
extern vector<double> Diaphragm;
extern vector<int> MidPiers,SidePiers;
extern double L;
extern double E;
extern VarSecBeam var_sec;
extern vector<vector<CPlateGeom> >plate3;
extern vector<size_t> SecDataCount,Mzero;
extern vector<double> deltas;
extern vector<CSecPropty2> xIdI;
extern Cmatrix<double> rst;
extern Ccsw csw;

#include "output.h"

void CalcElemq(CFem &fem,const CUnifLoad &q);
inline double CalcRou(const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv=50);

inline void ifs_ignore(istream & is,size_t LineCount)
{
	for(size_t i=1;i<=LineCount;++i)
		is.ignore(1000,'\n');
}

inline bool InputUnifSec(ifstream & ifs)  //变高度梁输入数据
{
	//读入跨径分布
	ifs_ignore(ifs,3);
	size_t tmp;
	ifs>>tmp;
	Spans.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>Spans[i];
		L+=Spans[i];
	}
	//读入中墩固结情况
	ifs_ignore(ifs,3);
	MidPiers.resize(tmp-1);
	for(size_t i=0;i<MidPiers.size();++i)
		ifs>>MidPiers[i];
	//读入边墩支座情况
	ifs_ignore(ifs,3);
	SidePiers.resize(2);
	for(size_t i=0;i<2;++i)
		ifs>>SidePiers[i];
	
	//读入隔板情况
	ifs_ignore(ifs,3);
	ifs>>tmp;
	Diaphragm.resize(tmp);
	for(size_t i=0;i<tmp;++i)
		ifs>>Diaphragm[i];

	std::sort(Diaphragm.begin(),Diaphragm.end());

	//读入固定参数
	ifs_ignore(ifs,4);
	ifs>>E;
	//读入3块翼缘板参数
	ifs_ignore(ifs,3);
	ifs>>tmp;
	var_sec.plate1.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>var_sec.plate1[i].x>>var_sec.plate1[i].thickness;
	}
	ifs_ignore(ifs,3);
	ifs>>tmp;
	var_sec.plate2.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>var_sec.plate2[i].x>>var_sec.plate2[i].thickness;
	}
	ifs_ignore(ifs,3);
	ifs>>tmp;
	var_sec.plate3_zd.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>var_sec.plate3_zd[i].x>>var_sec.plate3_zd[i].thickness;
	}
	//读入顶底板腹板位置
	ifs_ignore(ifs,3);
	ifs>>tmp;
	var_sec.web_top.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>var_sec.web_top[i].x>>var_sec.web_top[i].thickness;
	}
	ifs_ignore(ifs,3);
	ifs>>tmp;
	var_sec.web_bot_zd.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>var_sec.web_bot_zd[i].x>>var_sec.web_bot_zd[i].thickness;
	}
	//读入截面变化参数
	ifs_ignore(ifs,4);
	ifs>>var_sec.H_delta>>var_sec.Hz>>var_sec.Hc>>var_sec.curve.pala;

	//读入线形
	ifs_ignore(ifs,3);
	ifs>>tmp;
	var_sec.curve.keypoint.resize(tmp+1);
	for(size_t i=0;i<tmp;++i)
	{
		ifs>>var_sec.curve.keypoint[i].first>>var_sec.curve.keypoint[i].second;
		if(var_sec.curve.keypoint[i].first<0||var_sec.curve.keypoint[i].first>L+1e-8) 
		{
			cout<<"线形输入有误!\n";
			return false;
		}
	}
	var_sec.curve.keypoint.back().first=L+1;
	var_sec.curve.keypoint.back().second=var_sec.curve.keypoint[tmp-1].second;
	var_sec.PreCalculate();

	//读入均布荷载
	ifs_ignore(ifs,3);
	ifs>>tmp;
	UnifLoad.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		if(UnifLoad[i].x1<0||UnifLoad[i].x1>L+1e-8) return false;
		if(UnifLoad[i].x2<0||UnifLoad[i].x2>L+1e-8) return false;
		ifs>>UnifLoad[i].x1>>UnifLoad[i].x2>>UnifLoad[i].q1>>UnifLoad[i].q2;
		if(abs(UnifLoad[i].x2-L)<=1e-8) UnifLoad[i].x2=L-1e-8;
	}
	//读入集中荷载
	ifs_ignore(ifs,3);
	ifs>>tmp;
	ContLoad.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		if(ContLoad[i].x<0||ContLoad[i].x>L+1e-8) return false;
		ifs>>ContLoad[i].x>>ContLoad[i].f;
		if(abs(ContLoad[i].x-L)<=1e-8) ContLoad[i].x=L-1e-8;
	}
	//读入波形钢腹板厚度数据
	ifs_ignore(ifs,3);
	ifs>>tmp;
	csw._data.resize(tmp,2);
	for(size_t i=1;i<=tmp;++i)
	{
		ifs>>csw._data(i,1)>>csw._data(i,2);
	}
	ifs.close();
	return true;
}

inline void GenFemMod(CFem &fem,const size_t ndiv)
{
	fem.AddMat(CMaterial(1,"Concrete",E,E/2.0/(1+0.2)));
	set<Ckeypos> keypos;  //均布荷载、集中荷载、边界
	for(size_t i=0;i<var_sec.curve.keypoint.size()-2;++i)
	{
		keypos.insert(Ckeypos(var_sec.curve.keypoint[i].first,1));
	}
	keypos.insert(Ckeypos(0,0));
	{
		double tmp=0;
		for(size_t i=0;i<Spans.size();++i)
		{
			tmp+=Spans[i];
			keypos.insert(Ckeypos(tmp,0));
		}
	}
	for(size_t i=0;i<ContLoad.size();++i)
	{
		keypos.insert(Ckeypos(ContLoad[i].x,0));
	}
	//keypos.insert(Ckeypos(L,0));
	//计算delta_min用于控制I'计算的步长
	const double delta=L/ndiv;
	double delta_min=delta;
	set<Ckeypos>::iterator it1=keypos.begin(),it2=++keypos.begin(),it_end=--keypos.end();
	{
		double tmp;
		for(it1=keypos.begin();it1!=it_end;++it1)
		{
			tmp=it2->x-it1->x;
			if(tmp<delta_min) delta_min=tmp;
			++it2;
		}
	}
	delta_min/=max(2*ndiv+6,size_t(100));//该公式必定使得delta_min小于最小单元的长度的1/2
	//结点
	{
		double tmp;
		CSecPropty2 buffer;
		SecDataCount.clear();
		deltas.clear();
		size_t num=0;
		it1=keypos.begin(),it2=++keypos.begin();

		fem.AddNode(++num,it1->x);
		size_t ndiv2=(it2->x-it1->x)/delta;
		if (ndiv2==0) ++ndiv2;
		SecDataCount.push_back(ndiv2+1);
		deltas.push_back((it2->x-it1->x)/ndiv2);

		buffer.x=it1->x;	
		var_sec.Calculate(buffer.x+min(1e-6,delta_min),buffer);
		plate3.push_back(var_sec.plate3);
		buffer.dI=var_sec.Calculate_dI(buffer.x+min(1e-6,delta_min),min(1e-6,delta_min)*0.9);
		xIdI.push_back(buffer);
		for(size_t i=1;i<ndiv2;++i)
		{
			tmp=it1->x+i*deltas.back();
			fem.AddNode(++num,tmp,0);
			buffer.x=tmp;
			var_sec.Calculate(buffer.x,buffer);
			plate3.push_back(var_sec.plate3);
			buffer.dI=var_sec.Calculate_dI(buffer.x,delta_min);
			xIdI.push_back(buffer);
		}
		for(++it1;it1!=it_end;++it1)
		{
			++it2;
			fem.AddNode(++num,it1->x);
			size_t ndiv2=(it2->x-it1->x)/delta;
			if (ndiv2==0) ++ndiv2;
			SecDataCount.push_back(ndiv2+1);
			deltas.push_back((it2->x-it1->x)/ndiv2);

			buffer.x=it1->x;
			var_sec.Calculate(buffer.x,buffer);
			plate3.push_back(var_sec.plate3);
			plate3.push_back(var_sec.plate3);
			if(it1->flag)
			{
				buffer.dI=var_sec.Calculate_dI(buffer.x-delta_min,delta_min*0.9);
				xIdI.push_back(buffer);
				buffer.dI=var_sec.Calculate_dI(buffer.x+delta_min,delta_min*0.9);
				xIdI.push_back(buffer);
			}
			else 
			{
				buffer.dI=var_sec.Calculate_dI(buffer.x,delta_min);
				xIdI.push_back(buffer);
				xIdI.push_back(buffer);
			}

			for(size_t i=1;i<ndiv2;++i)
			{
				tmp=it1->x+i*deltas.back();
				fem.AddNode(++num,tmp,0);
				buffer.x=tmp;
				var_sec.Calculate(buffer.x,buffer);
				plate3.push_back(var_sec.plate3);
				buffer.dI=var_sec.Calculate_dI(buffer.x,delta_min);
				xIdI.push_back(buffer);
			}
		}
		fem.AddNode(++num,L,0);
		buffer.x=L;
		var_sec.Calculate(buffer.x-min(1e-6,delta_min),buffer);
		plate3.push_back(var_sec.plate3);
		buffer.dI=var_sec.Calculate_dI(buffer.x-min(1e-6,delta_min),min(1e-6,delta_min)*0.9);
		xIdI.push_back(buffer);
	}
	//单元
	{
		pair<double,double> AandI;
		std::stringstream ostr;
		size_t i=0;
		CFem::iterator_node it2=fem.NodeHead();
		for(CFem::iterator_node it=fem.NodeHead();it!=--fem.NodeEnd();++it)
		{
			ostr.str("");
			ostr<<"实常数"<<++i;
			++it2;
			var_sec.Calculate((it->x)/2+(it2->x)/2);
			fem.AddReal(CReal3d(i,ostr.str(),var_sec.sec_propty.A,var_sec.sec_propty.I,var_sec.sec_propty.I,
				var_sec.sec_propty.I,Cmatrix<long double>(1,2)));
			fem.AddElem(CBeam3d(i,0,"空间梁单元"),i,i+1,1,i);
		}
	}
	//边界
	/*if(SidePiers[0])
	{
		fem.AddBond(1,0,fem_aufunc::FZ);
		fem.AddBond(1,0,fem_aufunc::FY);
		fem.AddBond(1,0,fem_aufunc::MX);
		if(SidePiers[0]==2) fem.AddBond(1,0,fem_aufunc::MY);
	}
	if(SidePiers[1])
	{
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::FZ);
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::FY);
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::MX);
		if(SidePiers[1]==2) fem.AddBond(fem.NodeCount(),0,fem_aufunc::MY);
	}
	{
		double pos=0;
		size_t BdNodeNum;
		pos+=Spans[0];
		BdNodeNum=fem.NodeNearby(pos,0,0);
		fem.AddBond(BdNodeNum,0,fem_aufunc::FX);
		fem.AddBond(BdNodeNum,0,fem_aufunc::MX);
		fem.AddBond(BdNodeNum,0,fem_aufunc::FZ);
		fem.AddBond(BdNodeNum,0,fem_aufunc::FY);
		if(MidPiers[0]==2) fem.AddBond(BdNodeNum,0,fem_aufunc::MY);
		for(size_t i=1;i<MidPiers.size();++i)
		{
			pos+=Spans[i];
			BdNodeNum=fem.NodeNearby(pos,0,0);
			fem.AddBond(BdNodeNum,0,fem_aufunc::MX);
			fem.AddBond(BdNodeNum,0,fem_aufunc::FZ);
			fem.AddBond(BdNodeNum,0,fem_aufunc::FY);
			if(MidPiers[i]==2) fem.AddBond(BdNodeNum,0,fem_aufunc::MY);
		}
	}*/
	if(!MidPiers.empty())
	{
		fem.AddBond(fem.NodeNearby(Spans[0],0,0),0,fem_aufunc::FX);
	}
	else
	{
		if(SidePiers[0]) fem.AddBond(1,0,fem_aufunc::FX);
		else fem.AddBond(fem.NodeCount(),0,fem_aufunc::FX);
	}
	if(SidePiers[0])
	{
		fem.AddBond(1,0,fem_aufunc::FZ);
		fem.AddBond(1,0,fem_aufunc::FY);
		fem.AddBond(1,0,fem_aufunc::MX);
		if(SidePiers[0]==2) fem.AddBond(1,0,fem_aufunc::MY);
	}
	if(SidePiers[1])
	{
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::FZ);
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::FY);
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::MX);
		if(SidePiers[1]==2) fem.AddBond(fem.NodeCount(),0,fem_aufunc::MY);
	}
	{
		double pos=0;
		size_t BdNodeNum;
		for(size_t i=0;i<MidPiers.size();++i)
		{
			pos+=Spans[i];
			BdNodeNum=fem.NodeNearby(pos,0,0);
			fem.AddBond(BdNodeNum,0,fem_aufunc::MX);
			fem.AddBond(BdNodeNum,0,fem_aufunc::FZ);
			fem.AddBond(BdNodeNum,0,fem_aufunc::FY);
			if(MidPiers[i]==2) fem.AddBond(BdNodeNum,0,fem_aufunc::MY);
		}
	}
	//均布荷载
	for(size_t i=0;i<UnifLoad.size();++i)
	{
		CalcElemq(fem,UnifLoad[i]);
	}
	//集中荷载
	for(size_t i=0;i<ContLoad.size();++i)
	{
		fem.AddLoadn(fem.NodeNearby(ContLoad[i].x,0,0),ContLoad[i].f);
	}
	fem.Calculate(sparse_matrix());
}

inline void CalcElemq(CFem &fem,const CUnifLoad &q)
{
	vector<size_t> Nnode;
	CFem::const_iterator_node it1=fem.NodeHead();
	for(;it1!=fem.NodeEnd();++it1)
		if(it1->x>q.x1) break;
	--it1;
	CFem::const_iterator_node it2=fem.NodeEnd();
	--it2;
	for(;it2!=fem.NodeHead();--it2)
		if(it2->x<q.x2) break;
	CFem::const_iterator_node it(it1);
	++it;
	if(it1==it2)
	{
		fem.AddLoade(it1->Num(),q.q1,q.q2,(q.x1-it1->x)/(it->x-it1->x),(q.x2-it1->x)/(it->x-it1->x));
		return;
	}
	fem.AddLoade(it1->Num(),q.q1,q.Calcq(it->x),(q.x1-it1->x)/(it->x-it1->x),1);
	for(it1=it;it1!=it2;++it1)
	{
		++it;
		fem.AddLoade(it1->Num(),q.Calcq(it1->x),q.Calcq(it->x),0,1);
	}
	++it;
	fem.AddLoade(it1->Num(),q.Calcq(it1->x),q.q2,0,(q.x2-it1->x)/(it->x-it1->x));
}

inline void GenEqSet(CeqSetD & eqsd,const CFem &fem)  //计算公式为《单个U最终版》公式16
{
	const double G=E/2.0/(1+0.2);
	const double coef=-2.8*G/E;
	const double rval=-7.0/6/E;
	//(1)方程组
	Ceq eq;
	eq.AddItem(Citem(1,2));
	eq.AddItem(Citem(1,1));
	eq.AddItem(Citem(1,0));
	eqsd.AddEq(eq);
	//(2)节段数据
    eqsd.SetSecNum(SecDataCount);
	//(3)步长
	eqsd.SetDelta(deltas);
	//(4)mat的生成
	Cmatrix<double> mat;
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	eqsd.ReSizeCoef();
	{
		size_t index=0;
		double n;
		const double b1=var_sec.plate1.back().x;
		const double b2=var_sec.plate2.back().x;
		const double b3=var_sec.plate3_zd.back().x;
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],4);
			n=1.0/(1.0-7/8.0*xIdI[index].Is/xIdI[index].I);
			mat(1,1)=1;
			mat(1,2)=xIdI[index].dI/xIdI[index].Is;
			mat(1,3)=coef*n/xIdI[index].Is*(xIdI[index].I1/pow(b1,2)+xIdI[index].I2/pow(b2,2)+xIdI[index].I3/pow(b3,2));
			mat(1,4)=rval*n/xIdI[index].I*(*it_rst)(1,3);
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				n=1.0/(1.0-7/8.0*xIdI[index].Is/xIdI[index].I);
				mat(j,1)=1;
				mat(j,2)=xIdI[index].dI/xIdI[index].Is;
				mat(j,3)=coef*n/xIdI[index].Is*(xIdI[index].I1/pow(b1,2)+xIdI[index].I2/pow(b2,2)+xIdI[index].I3/pow(b3,2));
				mat(j,4)=rval*n/xIdI[index].I*(*it_rst)(2,3);
				++it_rst;
			}
			++index;
			eqsd.SetCoef(mat,1,i);
		}
	}
	//边界条件
	CeqBD eqbd;
	eqbd.AddItem(CitemBD(1,1,1,1,1));
	eqbd.SetValue(0);
	eqsd.AddEqb(eqbd);
	for(size_t i=1;i<SecDataCount.size();++i)
	{
		eqbd.initilize();
		eqbd.AddItem(CitemBD(i,1,1,SecDataCount[i-1],1));
		eqbd.AddItem(CitemBD(i+1,1,1,1,-1));
		eqbd.SetValue(0);
		eqsd.AddEqb(eqbd);
		eqbd.initilize();
		eqbd.AddItem(CitemBD(i,1,0,SecDataCount[i-1],1));
		eqbd.AddItem(CitemBD(i+1,1,0,1,-1));
		eqbd.SetValue(0);
		eqsd.AddEqb(eqbd);
	}
	eqbd.initilize();
	eqbd.AddItem(CitemBD(SecDataCount.size(),1,1,SecDataCount.back(),1));
	eqbd.SetValue(0);
	eqsd.AddEqb(eqbd);
	eqsd.Calculate();
}

void CalcShearLagAdv(const CFem& fem,const CeqSetD &eqsd)
{
	CeqSetD eqsd_rev;
	eqsd.Reverse(eqsd_rev);
	eqsd_rev.Calculate();
	ofstream ofs("rev2.dat");
	eqsd_rev.OutputAll(ofs);
	ofs.close();
	Mzero.clear();
	double MaxM=0;
	for(CFem::const_iterator_rst it_rst=fem.EForceHead();it_rst!=fem.EForceEnd();++it_rst)
	{
		if(MaxM<abs((*it_rst)(1,5))) MaxM=abs((*it_rst)(1,5));
		if(MaxM<abs((*it_rst)(2,5))) MaxM=abs((*it_rst)(2,5));
	}
	Cmatrix<double> u1,u1_rev,rst_rev;
	eqsd.CalcVarRst(1,1,u1);
	eqsd_rev.CalcVarRst(1,1,u1_rev);
	u1_rev.revered();
	u1_rev*=-1;

	rst.resize(u1.row(),8);
	rst_rev.resize(u1_rev.row(),1);
	eqsd.GetRst(rst,2);
	eqsd_rev.GetRst(rst_rev,1);
	rst_rev.revered();
	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,2)=rst(i,2)/2+rst_rev(i,1)/2;
	}

	CFem::const_iterator_node it_node=fem.NodeHead();
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	size_t pos=1;
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		for(size_t j=1;j<SecDataCount[i];++j)
		{
			rst(pos,1)=it_node->x;
			rst(pos,5)=u1(pos,1)/2+u1_rev(pos,1)/2;
			if(abs((*it_rst)(1,5))/MaxM>1e-6)
			{
				rst(pos,3)=1.0-3.0/4.0*E*xIdI[pos-1].Is/(-(*it_rst)(1,5))*rst(pos,5);
				rst(pos,4)=(*it_rst)(1,5);
			}
			else
			{
				rst(pos,3)=rst(pos,4)=0;
				Mzero.push_back(pos);
			}
			++it_node;
			++pos;
			++it_rst;
		}
		--it_rst;
		rst(pos,1)=it_node->x;
		rst(pos,5)=u1(pos,1)/2+u1_rev(pos,1)/2;
		if(abs((*it_rst)(2,5))/MaxM>1e-6)
		{
			rst(pos,3)=1.0-3.0/4.0*E*xIdI[pos-1].Is/(-(*it_rst)(2,5))*rst(pos,5);
			rst(pos,4)=(*it_rst)(2,5);
		}
		else
		{
			rst(pos,3)=rst(pos,4)=0;
			Mzero.push_back(pos);
		}
		++pos;
		++it_rst;
	}
	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,6)=CalcRou(i,1);
		rst(i,7)=CalcRou(i,2);
		rst(i,8)=CalcRou(i,3);
	}
}

inline double CalcSecStress(const size_t index,const double y,const double z,const double bi)
{
	return E*z*(rst(index,4)/E/xIdI[index-1].I-(1.0-pow(y/bi,3)-3.0/4.0*xIdI[index-1].Is/xIdI[index-1].I)*rst(index,5));//注意这里U'和理论推导是反方向的
}

inline double LineInt(const vector<double> &data,const double delta)
{
	double rtn=(data[0]+data.back())/2;
	for(size_t i=1;i<data.size()-1;++i)
	{
		rtn+=data[i];
	}
	return rtn*=delta;
}

inline double CalcStressLineInt(const size_t pos/*截面索引*/,const double y,const size_t index/*翼缘板标志1-3*/,
	const double L,const size_t ndiv=50) //竖向线积分
{
	/*const double L=CalcLength(y,index);*/
	const double delta=L/ndiv;
	double Zbeg=-xIdI[pos-1].C;
	if(index==1||index==2) Zbeg+=xIdI[pos-1].H;
	else Zbeg+=L;
	double bf;
	if(index==1) bf=var_sec.plate1.back().x;
	else if(index==2) bf=var_sec.plate2.back().x;
	else bf=var_sec.plate3_zd.back().x;
	vector<double> StressRst(ndiv+1);
	for(size_t i=0;i<=ndiv;++i)
	{
		StressRst[i]=CalcSecStress(pos,y,Zbeg,bf);
		Zbeg-=delta;
	}
	return LineInt(StressRst,delta);
}

inline double CalcLength(const double y,const vector<CPlateGeom> &plate)
{
	if (y<0||y>plate.back().x) return -1;
	if (y==0) return plate[0].thickness;
	size_t pos=0;
	for(;pos<plate.size();++pos)
	{
		if(plate[pos].x>=y) break;
	}
	return plate[pos-1].thickness+(plate[pos].thickness-plate[pos-1].thickness)
		/(plate[pos].x-plate[pos-1].x)*(y-plate[pos-1].x);
}

inline double CalcRou(const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv)
{
	//计算翼缘面积
	vector<CPlateGeom> plate;
	if(index==1) plate=var_sec.plate1;
	else if(index==2)  plate=var_sec.plate2;
	else  plate=plate3[pos-1];
	double area=0;
	for(size_t i=0;i<plate.size()-1;++i)
	{
		area+=(plate[i+1].x-plate[i].x)/2*(plate[i+1].thickness+plate[i].thickness);
	}
	//计算翼缘平均厚度
	double L=area/plate.back().x;
	//面积分
	double bf=plate.back().x;
	const double delta=bf/ndiv;
	vector<double> StrLineInt(ndiv+1);
	double y=0;
	for(size_t i=0;i<ndiv;++i)
	{
		StrLineInt[i]=CalcStressLineInt(pos,y,index,CalcLength(y,plate),1);
		y+=delta;
	}
	StrLineInt[ndiv]=CalcStressLineInt(pos,bf,index,CalcLength(bf,plate),1);
	double AreaInt=LineInt(StrLineInt,delta);
	double LineInt=CalcStressLineInt(pos,bf,index,L,1);
	return AreaInt/LineInt/bf;
}

inline void OutputOrig(ostream & os)
{
	os<<"跨径分布:\n";
	for(size_t i=0;i<Spans.size();++i)
	{
		os<<Spans[i]<<ends;
	}
	os<<"\n\n中墩固结情况:\n";
	for(size_t i=0;i<MidPiers.size();++i)
	{
		os<<MidPiers[i]<<ends;
	}
	os<<"\n\n边墩固结情况:\n";
	for(size_t i=0;i<SidePiers.size();++i)
	{
		os<<SidePiers[i]<<ends;
	}
	os<<"\n\n固定参数:\n";
	os<<"E:"<<E<<endl<<endl;

	var_sec.output(os);

	os<<"\n均布荷载:\n";
	for(size_t i=0;i<UnifLoad.size();++i)
	{
		os<<UnifLoad[i].x1<<ends<<UnifLoad[i].x2<<ends<<UnifLoad[i].q1<<ends<<UnifLoad[i].q2<<endl;
	}
	
	os<<"\n集中荷载:\n";
	for(size_t i=0;i<ContLoad.size();++i)
	{
		os<<ContLoad[i].x<<ends<<ContLoad[i].f<<endl;
	}	
}

