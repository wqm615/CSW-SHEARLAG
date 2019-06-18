//变高度梁相关函数文件
#include <set>
#include <sstream>
#include "StructDef.h"

extern vector<CContLoad> ContLoad;
extern vector<CBd> Bd;
extern vector<CUnifLoad> UnifLoad;
extern vector<double> Spans;
extern vector<int> MidPiers,SidePiers;
extern double L;
extern double E,b1,b2,b3;
extern CTopPlate top_plate1,top_plate2;  //顶板参数
extern CBotPlate bot_plate_s;            //支点底板参数
extern double B,h_delta,Hs,Hc;  //B底板全宽,h_delta底板厚度变化量,Hs支点梁高,Hc跨中梁高,抛物线次数
extern Ccurve curve;
extern vector<size_t> SecDataCount,Mzero;
extern vector<double> deltas;
extern vector<Cx_I_dI> xIdI;
extern Cmatrix<double> rst;

#include "output.h"

void CalcElemq(CFem &fem,const CUnifLoad &q);

inline void ifs_ignore(istream & is,size_t LineCount)
{
	for(size_t i=1;i<=LineCount;++i)
		is.ignore(1000,'\n');
}

inline pair<double,double> CalcAandI(double x)
{
	CBotPlate bot_plate=bot_plate_s.GetBotPlate(B,curve(-h_delta,0,x));
	return CSection(bot_plate,curve(Hc,Hs,x)).CalcPropety(top_plate1,top_plate2);
}

inline double CalcIdiff(double x,const double delta=0.01)
{
	double I1=CalcAandI(x-delta).second,I2=CalcAandI(x+delta).second;
	return (I2-I1)/delta/2;
}

inline void CalcI(double x,double&I1,double&I2,double&I3)
{
	CBotPlate bot_plate=bot_plate_s.GetBotPlate(B,curve(-h_delta,0,x));
	CSection(bot_plate,curve(Hc,Hs,x)).CalcI(top_plate1,top_plate2,I1,I2,I3);
}

inline bool InputUnifSec(ifstream & ifs)  //等高度梁输入数据
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
	//读入固定参数
	ifs_ignore(ifs,4);
	ifs>>E>>b1>>b2>>b3;
	//读入顶板参数
	ifs_ignore(ifs,4);
	ifs>>top_plate1.A>>top_plate1.I>>top_plate1.C;
	ifs>>top_plate2.A>>top_plate2.I>>top_plate2.C;
	//读入支点底板参数
	ifs_ignore(ifs,4);
	ifs>>bot_plate_s.A>>bot_plate_s.I>>bot_plate_s.C;
	//读入截面变化参数
	ifs_ignore(ifs,4);
	ifs>>B>>h_delta>>Hs>>Hc>>curve.pala;
	//读入线形
	ifs_ignore(ifs,3);
	ifs>>tmp;
	curve.keypoint.resize(tmp+1);
	for(size_t i=0;i<tmp;++i)
	{
		if(curve.keypoint[i].first<0||curve.keypoint[i].first>L) return false;
		ifs>>curve.keypoint[i].first>>curve.keypoint[i].second;
	}
	curve.keypoint.back().first=L+1;
	curve.keypoint.back().second=curve.keypoint[tmp-1].second;
	//读入均布荷载
	ifs_ignore(ifs,3);
	ifs>>tmp;
	UnifLoad.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		if(UnifLoad[i].x1<0||UnifLoad[i].x1>L) return false;
		if(UnifLoad[i].x2<0||UnifLoad[i].x2>L) return false;
		ifs>>UnifLoad[i].x1>>UnifLoad[i].x2>>UnifLoad[i].q1>>UnifLoad[i].q2;
	}
	//读入集中荷载
	ifs_ignore(ifs,3);
	ifs>>tmp;
	ContLoad.resize(tmp);
	for(size_t i=0;i<tmp;++i)
	{
		if(ContLoad[i].x<0||ContLoad[i].x>L) return false;
		ifs>>ContLoad[i].x>>ContLoad[i].f;
	}
	ifs.close();
	return true;
}

inline void GenFemMod(CFem &fem,const size_t ndiv)
{
	fem.AddMat(CMaterial(1,"Concrete",E,E/2.0/(1+0.2)));
	set<Ckeypos> keypos;  //均布荷载、集中荷载、边界
	for(size_t i=0;i<curve.keypoint.size()-1;++i)
	{
		keypos.insert(Ckeypos(curve.keypoint[i].first,1));
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
	keypos.insert(Ckeypos(L,0));
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
		Cx_I_dI buffer;
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
		CalcI(buffer.x+min(1e-6,delta_min),buffer.I1,buffer.I2,buffer.I3);
		buffer.dI=CalcIdiff(buffer.x+1.1*delta_min,delta_min);
		xIdI.push_back(buffer);
		for(size_t i=1;i<ndiv2;++i)
		{
			tmp=it1->x+i*deltas.back();
			fem.AddNode(++num,tmp,0);
			buffer.x=tmp;
			CalcI(buffer.x,buffer.I1,buffer.I2,buffer.I3);
			buffer.dI=CalcIdiff(buffer.x,delta_min);
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
			CalcI(buffer.x,buffer.I1,buffer.I2,buffer.I3);
			if(it1->flag)
			{
				buffer.dI=CalcIdiff(buffer.x-1.1*delta_min,delta_min);
				xIdI.push_back(buffer);
				buffer.dI=CalcIdiff(buffer.x+1.1*delta_min,delta_min);
				xIdI.push_back(buffer);
			}
			else 
			{
				buffer.dI=CalcIdiff(buffer.x,delta_min);
				xIdI.push_back(buffer);
				xIdI.push_back(buffer);
			}

			for(size_t i=1;i<ndiv2;++i)
			{
				tmp=it1->x+i*deltas.back();
				fem.AddNode(++num,tmp,0);
				buffer.x=tmp;
				CalcI(buffer.x,buffer.I1,buffer.I2,buffer.I3);
				buffer.dI=CalcIdiff(buffer.x,delta_min);
				xIdI.push_back(buffer);
			}
		}
		fem.AddNode(++num,L,0);
		buffer.x=L;
		CalcI(buffer.x-min(1e-6,delta_min),buffer.I1,buffer.I2,buffer.I3);
		buffer.dI=CalcIdiff(buffer.x-1.1*delta_min,delta_min);
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
			AandI=CalcAandI((it->x)/2+(it2->x)/2);
			fem.AddReal(CReal3d(i,ostr.str(),AandI.first,AandI.second,AandI.second,AandI.second,Cmatrix<long double>(1,2)));
			fem.AddElem(CBeam3d(i,0,"空间梁单元"),i,i+1,1,i);
		}
	}
	//边界
	fem.AddBond(1,0,fem_aufunc::FX);
	if(SidePiers[0])
	{
		fem.AddBond(1,0,fem_aufunc::FZ);
		fem.AddBond(1,0,fem_aufunc::MX);
		if(SidePiers[0]==2) fem.AddBond(1,0,fem_aufunc::MY);
	}
	if(SidePiers[1])
	{
		fem.AddBond(fem.NodeCount(),0,fem_aufunc::FZ);
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

inline void GenEqSet(CeqSetD & eqsd,const CFem &fem)
{
	const double G=E/2.0/(1+0.2);
	const double coef=-22.4*G/E;
	const double rval=-28.0/3/E;
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
		double Is;
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],4);
			Is=xIdI[index].I1+xIdI[index].I2+xIdI[index].I3;
			mat(1,1)=1;
			mat(1,2)=xIdI[index].dI/Is;
			mat(1,3)=coef/Is*(xIdI[index].I1/pow(b1,2)+xIdI[index].I2/pow(b2,2)+xIdI[index].I3/pow(b3,2));
			mat(1,4)=rval/Is*(*it_rst)(1,3);
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				Is=xIdI[index].I1+xIdI[index].I2+xIdI[index].I3;
				mat(j,1)=1;
				mat(j,2)=xIdI[index].dI/Is;
				mat(j,3)=coef/Is*(xIdI[index].I1/pow(b1,2)+xIdI[index].I2/pow(b2,2)+xIdI[index].I3/pow(b3,2));
				mat(j,4)=rval/Is*(*it_rst)(2,3);
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

void CalcShearLag(const CFem& fem,const CeqSetD &eqsd)
{
	Mzero.clear();
	double MaxM=0;
	for(CFem::const_iterator_rst it_rst=fem.EForceHead();it_rst!=fem.EForceEnd();++it_rst)
	{
		if(MaxM<abs((*it_rst)(1,5))) MaxM=abs((*it_rst)(1,5));
		if(MaxM<abs((*it_rst)(2,5))) MaxM=abs((*it_rst)(2,5));
	}
	Cmatrix<double> u1;
	eqsd.CalcVarRst(1,1,u1);
	rst.resize(u1.row(),8);
	eqsd.GetRst(rst,2);
	CFem::const_iterator_node it_node=fem.NodeHead();
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	size_t pos=1;
	double Is;
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		Is=xIdI[pos-1].I1+xIdI[pos-1].I2+xIdI[pos-1].I3;
		rst(pos,1)=(it_node)->x;
		if(abs((*it_rst)(1,5))/MaxM>1e-6)
		{
			rst(pos,3)=1.0-3.0/4.0*E*Is/(-(*it_rst)(1,5))*u1(pos,1);
			rst(pos,4)=(*it_rst)(1,5);
		}
		else 
		{
			rst(pos,3)=rst(pos,4)=0;
			Mzero.push_back(pos);
		}
		++pos;
		++it_node;
		for(size_t j=2;j<SecDataCount[i];++j)
		{
			++it_rst;
			Is=xIdI[pos-1].I1+xIdI[pos-1].I2+xIdI[pos-1].I3;
			rst(pos,1)=it_node->x;
			if(abs((*it_rst)(1,5))/MaxM>1e-6)
			{
				rst(pos,3)=1.0-3.0/4.0*E*Is/(-(*it_rst)(1,5))*u1(pos,1);
				rst(pos,4)=(*it_rst)(1,5);
			}
			else
			{
				rst(pos,3)=rst(pos,4)=0;
				Mzero.push_back(pos);
			}
			++it_node;
			++pos;
		}
		Is=xIdI[pos-1].I1+xIdI[pos-1].I2+xIdI[pos-1].I3;
		rst(pos,1)=it_node->x;
		if(abs((*it_rst)(2,5))/MaxM>1e-6)
		{
			rst(pos,3)=1.0-3.0/4.0*E*Is/(-(*it_rst)(2,5))*u1(pos,1);
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
	double Is;
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		for(size_t j=1;j<SecDataCount[i];++j)
		{
			Is=xIdI[pos-1].I1+xIdI[pos-1].I2+xIdI[pos-1].I3;
			rst(pos,1)=it_node->x;
			rst(pos,5)=u1(pos,1)/2+u1_rev(pos,1)/2;
			if(abs((*it_rst)(1,5))/MaxM>1e-6)
			{
				rst(pos,3)=1.0-3.0/4.0*E*Is/(-(*it_rst)(1,5))*rst(pos,5);
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
		Is=xIdI[pos-1].I1+xIdI[pos-1].I2+xIdI[pos-1].I3;
		rst(pos,1)=it_node->x;
		rst(pos,5)=u1(pos,1)/2+u1_rev(pos,1)/2;
		if(abs((*it_rst)(2,5))/MaxM>1e-6)
		{
			rst(pos,3)=1.0-3.0/4.0*E*Is/(-(*it_rst)(2,5))*rst(pos,5);
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
	os<<"E\tb1\tb2\tb3\n";
	os<<E<<ends<<b1<<ends<<b2<<ends<<b3<<endl;

	os<<"\n顶板参数:\n";
	os<<"A1\tI1\tC1\tA2\tI2\tC2\n";
	os<<top_plate1.A<<ends<<top_plate1.I<<ends<<top_plate1.C<<ends;
	os<<top_plate2.A<<ends<<top_plate2.I<<ends<<top_plate2.C<<endl;

	os<<"\n顶板参数:\n";
	os<<"A3\tI3\tC3\n";
	os<<bot_plate_s.A<<ends<<bot_plate_s.I<<ends<<bot_plate_s.C<<endl;

	os<<"\n截面变化参数:\n";
	os<<"B\tΔh\tH支\tH中\t曲线次数\n";
	os<<B<<ends<<h_delta<<ends<<Hs<<ends<<Hc<<ends<<curve.pala<<endl;

	os<<"\n线形:\n";
	for(size_t i=0;i<curve.keypoint.size()-1;++i)
		os<<curve.keypoint[i].first<<ends<<curve.keypoint[i].second<<endl;

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
