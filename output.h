#include "CeqSetD.h"
#include <fem\fem.h>

inline void OutputOrig(ostream & os);
inline void OutputFem(const CFem& fem,ostream &ofs_fem)
{
	fem.OutputData(ofs_fem);
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	ofs_fem<<left<<"\n单元内力结果如下:\n"<<setw(4)<<"单元"<<setw(12)<<"Qi"<<setw(12)<<"Qj"<<setw(12)<<"Mi"<<setw(12)<<"Mj"<<endl;
	for(size_t i=1;it_rst!=fem.EForceEnd();++it_rst)
	{
		ofs_fem<<left<<setw(4)<<i<<setw(12)<<(*it_rst)(1,3)<<setw(12)<<(*it_rst)(2,3)<<
			setw(12)<<-(*it_rst)(1,5)<<setw(12)<<-(*it_rst)(2,5)<<endl;
		++i;
	}
	ofs_fem<<left<<"\n节点位移结果如下:\n"<<setw(4)<<"节点"<<setw(12)<<"Uz"<<setw(12)<<"roty"<<endl;
	it_rst=fem.NDispHead();
	for(size_t i=1;it_rst!=fem.NDispEnd();++it_rst)
	{
		ofs_fem<<left<<setw(6)<<i<<setw(16)<<(*it_rst)(1,3)<<setw(16)<<(*it_rst)(1,5)<<endl;
		++i;
	}
}

inline void OutputRst(ostream&ofs_rst)
{
	Cmatrix<double> MaxAndMin(2,2);
	if(Mzero.empty()) 
	{
		MaxAndMin(1,1)=MaxAndMin(2,1)=rst(1,1);
		MaxAndMin(1,2)=MaxAndMin(2,2)=rst(1,3);
	}
	else
	{
		for(size_t i=1;i<=rst.row();++i)
		{
			if(i!=Mzero[0])
			{
				MaxAndMin(1,1)=MaxAndMin(2,1)=rst(i,1);
				MaxAndMin(1,2)=MaxAndMin(2,2)=rst(i,3);
				break;
			}
		}
	}
	size_t pos_index=1;
	
	ofs_rst<<"坐标,U,λ,M,U',ρ1,ρ2,ρ3\n";
	for(size_t j=0;j<Mzero.size();++j)
	{
		for(size_t i=pos_index;i<Mzero[j];++i)
		{
			ofs_rst<<rst(i,1)<<","<<rst(i,2)<<","<<rst(i,3)<<","<<rst(i,4)<<","<<rst(i,5)<<","<<rst(i,6);
			ofs_rst<<","<<rst(i,7)<<","<<rst(i,8)<<endl;
			if(MaxAndMin(1,2)<rst(i,3))
			{
				MaxAndMin(1,2)=rst(i,3);
				MaxAndMin(1,1)=rst(i,1);
			}
			if(MaxAndMin(2,2)>rst(i,3))
			{
				MaxAndMin(2,2)=rst(i,3);
				MaxAndMin(2,1)=rst(i,1);
			}
		}
		ofs_rst<<rst(Mzero[j],1)<<","<<rst(Mzero[j],2)<<","<<"NULL,"<<rst(Mzero[j],4)<<","<<rst(Mzero[j],5);
		ofs_rst<<",NULL,NULL,NULL\n";
		pos_index=Mzero[j]+1;
	}
	for(size_t i=pos_index;i<=rst.row();++i)
	{
		ofs_rst<<rst(i,1)<<","<<rst(i,2)<<","<<rst(i,3)<<","<<rst(i,4)<<","<<rst(i,5)<<","<<rst(i,6);
		ofs_rst<<","<<rst(i,7)<<","<<rst(i,8)<<endl;
		if(MaxAndMin(1,2)<rst(i,3))
		{
			MaxAndMin(1,2)=rst(i,3);
			MaxAndMin(1,1)=rst(i,1);
		}
		if(MaxAndMin(2,2)>rst(i,3))
		{
			MaxAndMin(2,2)=rst(i,3);
			MaxAndMin(2,1)=rst(i,1);
		}
	}
	ofs_rst<<endl<<"最大值: ,坐标, "<<MaxAndMin(1,1)<<",λ,"<<MaxAndMin(1,2)<<endl;
	ofs_rst<<"最小值: ,坐标,"<<MaxAndMin(2,1)<<",λ,"<<MaxAndMin(2,2)<<endl;
}

inline void OutputRst3U(ostream&ofs_rst)
{
	Cmatrix<double> MaxAndMin(2,2);
	if(Mzero.empty()) 
	{
		MaxAndMin(1,1)=MaxAndMin(2,1)=rst(1,1);
		MaxAndMin(1,2)=MaxAndMin(2,2)=rst(1,9);
	}
	else
	{
		for(size_t i=1;i<=rst.row();++i)
		{
			if(i!=Mzero[0])
			{
				MaxAndMin(1,1)=MaxAndMin(2,1)=rst(i,1);
				MaxAndMin(1,2)=MaxAndMin(2,2)=rst(i,9);
				break;
			}
		}
	}
	size_t pos_index=1;
	
	ofs_rst<<"坐标,U1,U2,U3,U1',U2',U3',M,λ,ρ1,ρ2,ρ3\n";
	for(size_t j=0;j<Mzero.size();++j)
	{
		for(size_t i=pos_index;i<Mzero[j];++i)
		{
			ofs_rst<<rst(i,1)<<","<<rst(i,2)<<","<<rst(i,3)<<","<<rst(i,4)<<","<<rst(i,5)<<","<<rst(i,6);
			ofs_rst<<","<<rst(i,7)<<","<<rst(i,8)<<","<<rst(i,9)<<","<<rst(i,10)<<","<<rst(i,11)<<","<<rst(i,12)<<endl;
			if(MaxAndMin(1,2)<rst(i,9))
			{
				MaxAndMin(1,2)=rst(i,9);
				MaxAndMin(1,1)=rst(i,1);
			}
			if(MaxAndMin(2,2)>rst(i,9))
			{
				MaxAndMin(2,2)=rst(i,9);
				MaxAndMin(2,1)=rst(i,1);
			}
		}
		ofs_rst<<rst(Mzero[j],1)<<","<<rst(Mzero[j],2)<<","<<rst(Mzero[j],3)<<","<<rst(Mzero[j],4)<<","<<rst(Mzero[j],5)
			<<","<<rst(Mzero[j],6)<<","<<rst(Mzero[j],7)<<","<<rst(Mzero[j],8)<<",NULL,NULL,NULL,NULL\n";
		pos_index=Mzero[j]+1;
	}
	for(size_t i=pos_index;i<=rst.row();++i)
	{
		ofs_rst<<rst(i,1)<<","<<rst(i,2)<<","<<rst(i,3)<<","<<rst(i,4)<<","<<rst(i,5)<<","<<rst(i,6);
			ofs_rst<<","<<rst(i,7)<<","<<rst(i,8)<<","<<rst(i,9)<<","<<rst(i,10)<<","<<rst(i,11)<<","<<rst(i,12)<<endl;
		if(MaxAndMin(1,2)<rst(i,9))
		{
			MaxAndMin(1,2)=rst(i,9);
			MaxAndMin(1,1)=rst(i,1);
		}
		if(MaxAndMin(2,2)>rst(i,9))
		{
			MaxAndMin(2,2)=rst(i,9);
			MaxAndMin(2,1)=rst(i,1);
		}
	}
	ofs_rst<<endl<<"最大值: ,坐标, "<<MaxAndMin(1,1)<<",λ,"<<MaxAndMin(1,2)<<endl;
	ofs_rst<<"最小值: ,坐标,"<<MaxAndMin(2,1)<<",λ,"<<MaxAndMin(2,2)<<endl;
}

inline void Output(string FileName,const CFem &fem, const CeqSetD& eqsd)
{
	ofstream ofs_fem(FileName+"fem.out");
	ofstream ofs_cf(FileName+"cf.dat");
	ofstream ofs_rst(FileName+"rst.csv");
	OutputOrig(ofs_cf);
	eqsd.OutputAll(ofs_cf);
	OutputFem(fem,ofs_fem);
	OutputRst(ofs_rst);
	ofs_fem.close();
	ofs_cf.close();
	ofs_rst.close();
}

inline void Output3U(string FileName,const CFem &fem, const CeqSetD& eqsd)
{
	ofstream ofs_fem(FileName+"fem.out");
	ofstream ofs_cf(FileName+"cf.dat");
	ofstream ofs_rst(FileName+"rst.csv");
	ofstream ofs_orig(FileName+"orig.dat");
	OutputOrig(ofs_orig);
	eqsd.OutputAll(ofs_cf);
	OutputFem(fem,ofs_fem);
	OutputRst3U(ofs_rst);
	ofs_fem.close();
	ofs_cf.close();
	ofs_rst.close();
	ofs_orig.close();
}
