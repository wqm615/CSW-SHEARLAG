//等高度梁相关函数文件
#include <set>
#include "StructDef.h"

void CalcElemq(CFem &fem,const CUnifLoad &q);

extern vector<CContLoad> ContLoad;
extern vector<CBd> Bd;
extern vector<CUnifLoad> UnifLoad;
extern double Area,L,I1,I2,I3,I,b1,b2,b3,E,Is,Iw,SecDepth/*梁高*/,SecZc/*形心轴高度*/;
extern vector<size_t> SecDataCount,Mzero;;
extern vector<double> deltas;
extern Cmatrix<double> rst;
extern vector<CPlateGeom> plate[3];
#include "output.h"

inline void InputUnifSec(ifstream & ifs)  //等高度梁输入数据
{
	size_t BdCount,UnfLoadCount,ContLoadCount;
	string buffer;
	ifs>>buffer>>BdCount;
	ifs>>buffer>>buffer;
	Bd.resize(BdCount);
	for(size_t i=0;i<BdCount;++i)
	{
		ifs>>Bd[i].x>>Bd[i].Fd;
	}
	ifs>>buffer>>UnfLoadCount;
	UnifLoad.resize(UnfLoadCount);
	for(size_t i=0;i<UnfLoadCount;++i)
	{
		ifs>>UnifLoad[i].x1>>UnifLoad[i].x2>>UnifLoad[i].q1>>UnifLoad[i].q2;
	}
	ifs>>buffer>>ContLoadCount;
	ContLoad.resize(ContLoadCount);
	for(size_t i=0;i<ContLoadCount;++i)
	{
		ifs>>ContLoad[i].x>>ContLoad[i].f;
	}
}

inline void GenFemMod(CFem &fem,const size_t ndiv)
{
	//材料几何参数
	fem.AddMat(CMaterial(1,"Concrete",E,E/2.0/(1+0.2)));
	fem.AddReal(CReal3d(1,"Null",Area,I,I,I,Cmatrix<long double>(1,2)));
	set<double> keypos;  //均布荷载、集中荷载、边界
	for(size_t i=0;i<Bd.size();++i)
	{
		if(Bd[i].Fd==2)  //20160517新加入的条件，之前发布的程序版本，没有该条件！
		keypos.insert(Bd[i].x);
	}
	keypos.insert(0);
	for(size_t i=0;i<ContLoad.size();++i)
	{
		keypos.insert(ContLoad[i].x);
	}
	keypos.insert(L);
	SecDataCount.clear();
	deltas.clear();
	const double delta=L/ndiv;
	set<double>::iterator it1=keypos.begin(),it2=keypos.begin(),it_end=keypos.end();
	--it_end;
	size_t num=0;
	//结点
	for(it1=keypos.begin();it1!=it_end;++it1)
	{
		++it2;
		fem.AddNode(++num,*it1);
		size_t ndiv2=(*it2-*it1)/delta;
		if (ndiv2==0) ++ndiv2;
		SecDataCount.push_back(ndiv2+1);
		deltas.push_back((*it2-*it1)/ndiv2);
		for(size_t i=1;i<ndiv2;++i)
			fem.AddNode(++num,*it1+i*deltas.back(),0);
	}
	fem.AddNode(++num,L,0);
	//单元
	for(size_t i=1;i<fem.NodeCount();++i)
		fem.AddElem(CBeam3d(i,0,"空间梁单元"),i,i+1,1,1);
	
	//边界
	for(size_t i=0;i<Bd.size();++i)
	{
		size_t BdNodeNum=fem.NodeNearby(Bd[i].x,0,0);
		fem.AddBond(BdNodeNum,0,static_cast<fem_aufunc::INDEX>(Bd[i].Fd));
		fem.AddBond(BdNodeNum,0,fem_aufunc::FY);
		fem.AddBond(BdNodeNum,0,fem_aufunc::MX);
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
	const double n=1.0/(1-7.0/8*Is/I);
	const double coef=-14*n*G/5/E/Is*(I1/b1/b1+I2/b2/b2+I3/b3/b3);
	const double rval=-7*n/6/E/I;
	//(1)方程组
	Ceq eq;
	eq.AddItem(Citem(1,2));
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
	for(size_t i=1;i<=SecDataCount.size();++i)
	{
		mat.resize(SecDataCount[i-1],3);
		mat(1,1)=1;
		mat(1,2)=coef;
		mat(1,3)=rval*(*it_rst)(1,3);
		for(size_t j=2;j<=SecDataCount[i-1];++j)
		{
			mat(j,1)=1;
			mat(j,2)=coef;
			mat(j,3)=rval*(*it_rst)(2,3);
			++it_rst;
		}
		eqsd.SetCoef(mat,1,i);
	}
	//(5)边界条件
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


inline void CalcShearLag(const CFem& fem,const CeqSetD &eqsd)
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
	rst.resize(u1.row(),4);
	eqsd.GetRst(rst,2);
	CFem::const_iterator_node it_node=fem.NodeHead();
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	size_t pos=1;
	for(size_t i=0;i<SecDataCount.size();++i)
	{
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

inline void CalcShearLagAdv(const CFem& fem,const CeqSetD &eqsd)  //取向前、向后（通过逆序计算）的平均值
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
		
	rst.resize(u1.row(),5);//第5列为u的一阶导数
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
	os<<"L,E,Area,I1,I2,I3,Iw,b1,b2,b3\n";
	os<<L<<ends<<E<<ends<<Area<<ends<<I1<<ends<<I2<<ends<<I3<<ends<<Iw<<ends<<b1<<ends<<b2<<ends<<b3<<endl;
	os<<"\n边界条件:\n";
	for(size_t i=0;i<Bd.size();++i)
	{
		os<<Bd[i].x<<ends<<fem_aufunc::MOMENT[Bd[i].Fd]<<endl;
	}
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

inline double CalcSecStress(const size_t index,const double y,const double z,const double bi)
{
	return E*z*(rst(index,4)/E/I-(1.0-pow(y/bi,3)-3.0/4.0*Is/I)*rst(index,5));
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
	double Zbeg=-SecZc;
	if(index==1||index==2) Zbeg+=SecDepth;
	else Zbeg+=L;
	double bf;
	if(index==1) bf=b1;
	else if(index==2) bf=b2;
	else bf=b3;
	vector<double> StressRst(ndiv+1);
	for(size_t i=0;i<=ndiv;++i)
	{
		StressRst[i]=CalcSecStress(pos,y,Zbeg,bf);
		Zbeg-=delta;
	}
	return LineInt(StressRst,delta);
}

inline double CalcLength(const double y,const size_t index) //i=1 2 3分别表示翼缘板1 2 3
{
	if (y<0||y>plate[index-1].back().x) return -1;
	if (y==0) return plate[index-1][0].thickness;
	size_t pos=0;
	for(;pos<plate[index-1].size();++pos)
	{
		if(plate[index-1][pos].x>=y) break;
	}
	return plate[index-1][pos-1].thickness+(plate[index-1][pos].thickness-plate[index-1][pos-1].thickness)
		/(plate[index-1][pos].x-plate[index-1][pos-1].x)*(y-plate[index-1][pos-1].x);
}

inline double CalcRou(const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv=50)
{
	//计算翼缘面积
	double area=0;
	for(size_t i=0;i<plate[index-1].size()-1;++i)
	{
		area+=(plate[index-1][i+1].x-plate[index-1][i].x)/2*(plate[index-1][i+1].thickness+plate[index-1][i].thickness);
	}
	//计算翼缘平均厚度
	double L=area/plate[index-1].back().x;
	//面积分
	double bf;
	if(index==1) bf=b1;
	else if(index==2) bf=b2;
	else bf=b3;
	const double delta=bf/ndiv;
	vector<double> StrLineInt(ndiv+1);
	double y=0;
	for(size_t i=0;i<ndiv;++i)
	{
		StrLineInt[i]=CalcStressLineInt(pos,y,index,CalcLength(y,index),2);
		y+=delta;
	}
	StrLineInt[ndiv]=CalcStressLineInt(pos,bf,index,CalcLength(bf,index),2);
	double AreaInt=LineInt(StrLineInt,delta);
	double LineInt=CalcStressLineInt(pos,bf,index,L,2);
	return AreaInt/LineInt/bf;
}