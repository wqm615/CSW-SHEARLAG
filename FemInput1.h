#include <fem\fem.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace fem_aufunc;

inline void InputData2d(CFem2d & fem,istream & ifs)
{
	size_t num_node,num_mat,num_real,num_elem,num_fe,num_qe,num_fn,num_bd;
	/*string buffer;
	ifs>>buffer;*/
	ifs>>num_node;
	for(size_t i=1;i<=num_node;i++)
	{
		size_t num;
		long double x,y;
		ifs>>num>>x>>y;
		fem.AddNode(num,x,y);
	}
	
	ifs>>num_mat;
	for(size_t i=1;i<=num_mat;i++)
	{
		size_t num;
		string name;
		long double e;
		ifs>>num>>name>>e;
		fem.AddMat(CMaterial(num,name,e));
	}

	ifs>>num_real;
	for(size_t i=1;i<=num_real;i++)
	{
		size_t num;
		string name;
		long double area,I;
		ifs>>num>>name>>area>>I;
		fem.AddReal(CReal2d(num,area,I,Cmatrix<long double>(),name));
	}

	ifs>>num_elem;
	for(size_t i=1;i<=num_elem;i++)
	{
		size_t num;
		string name;
		size_t ni,nj,nmat,nreal;
		ifs>>num>>name>>ni>>nj>>nmat>>nreal;
		fem.AddElem(num,ni,nj,nmat,nreal,name);
	}

	ifs>>num_fe;
	for(size_t i=1;i<=num_fe;++i)
	{
		long double value;
		long double dist;
		size_t num,index;
		ifs>>num>>value>>dist>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoade(num,value,dist,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_qe;
	for(size_t i=1;i<=num_qe;++i)
	{
		long double q1,q2;
		long double d1,d2;
		size_t num,index;
		ifs>>num>>q1>>q2>>d1>>d2>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoade(num,q1,q2,d1,d2,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_fn;
	for(size_t i=1;i<=num_fn;++i)
	{
		long double value;
		size_t num,index;
		ifs>>num>>value>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoadn(num,value,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_bd;
	for(size_t i=1;i<=num_bd;++i)
	{
		long double value;
		size_t node_num,num;
		ifs>>node_num>>value>>num;
		num=num+fem_aufunc::FX-1;
		fem.AddBond(node_num,value,static_cast<fem_aufunc::INDEX>(num));
	}
}

inline void OutputData2d(CFem2d &fem,ostream&os)
{
	os<<"������Ϣ����:\n";
	os<<fem.NodeCount()<<ends<<fem.ElemCount()<<ends<<fem.MatCount()<<ends<<fem.RealCount()<<endl;
	os<<endl<<"�����Ϣ����:\n";
	for(CFem2d::iterator_node it=fem.NodeHead();it!=fem.NodeEnd();++it)
	{
		os<<it->Num()<<ends<<it->x<<ends<<it->y<<endl;
	}

	os<<endl<<"������Ϣ����:\n";
	for(CFem2d::iterator_mat it=fem.MatHead();it!=fem.MatEnd();++it)
	{
		os<<it->Num()<<ends<<it->name<<ends<<it->E<<ends<<it->G<<endl;
	}

	os<<endl<<"ʵ������Ϣ����:\n";
	for(CFem2d::iterator_real it=fem.RealHead();it!=fem.RealEnd();++it)
	{
		os<<it->Num()<<ends<<it->name<<ends<<it->area<<ends<<it->iyy<<endl;
	}

	cout<<endl<<"��Ԫ��Ϣ����:\n";
	for(CFem2d::const_iterator_elem it=fem.ElemHead();it!=fem.ElemEnd();++it)
	{
		os<<it->Num()<<ends<<it->name<<ends<<it->NumNodei()<<ends<<it->NumNodej()<<ends<<it->NumMat()<<ends<<it->NumReal()<<endl;
	}

	cout<<"\n��Ԫ���к�������:\n";
	for(CFem2d::const_iterator_Fe it=fem.FeHead();it!=fem.FeEnd();++it)
	{os<<it->GetNumElem()<<ends<<it->value<<ends<<it->dist<<ends<<fem_aufunc::FORCE[it->type]<<endl;}

	os<<"\n��Ԫ�ߺ�������:\n";
	for(CFem2d::const_iterator_qe it=fem.QeHead();it!=fem.QeEnd();++it)
	{os<<it->GetNumElem()<<ends<<it->q1<<ends<<it->q2<<ends<<it->disti<<ends<<it->distj<<ends<<fem_aufunc::FORCE[it->type]<<endl;}

	os<<"\n����������:\n";
	for(CFem2d::const_iterator_Fn it=fem.FnHead();it!=fem.FnEnd();++it)
	{os<<it->GetNumNode()<<ends<<it->value<<ends<<fem_aufunc::FORCE[it->type]<<endl;}

	os<<"\nԼ����Ϣ����:\n";
	for(CFem2d::const_iterator_bd it=fem.BondHead();it!=fem.BondEnd();++it)
	{
		CBond::bd_type rst;
		os<<it->Num()<<ends;
		it->GetBond(rst);
		for(CBond::bd_type::iterator _it=rst.begin();_it!=rst.end();++_it)
		{os<<_it->second<<ends<<fem_aufunc::FORCE[_it->first]<<",";}
		os<<endl;
	}
}

inline void InputDataCsw(CFem & fem,ifstream & ifs)
{
	if(!ifs) return;
	size_t num_node,num_mat/*,num_real*/,num_elem,num_fe,num_qe,num_fn,num_bd;

	ifs>>num_mat;
	for(size_t i=1;i<=num_mat;i++)
	{
		size_t num;
		string name;
		long double e,g;
		ifs>>num>>name>>e>>g;
		g=e/2/(1+g);
		fem.AddMat(CMaterial(num,name,e,g));
	}

	ifs>>num_node;
	for(size_t i=1;i<=num_node;i++)
	{
		size_t num;
		long double x,y,z;
		ifs>>num>>x>>y>>z;
		fem.AddNode(num,x,y,z);
	}
	
	ifs>>num_elem;
	for(size_t i=1;i<=num_elem;i++)
	{
		size_t num;
		string name;
		size_t ni,nj,nmat,nreal;
		long double theta;
		ifs>>num>>name>>ni>>nj>>nmat>>nreal>>theta;
		fem.AddElem(CBeamcsw(num,theta,name),ni,nj,nmat,nreal);
	}

	ifs>>num_fe;
	for(size_t i=1;i<=num_fe;++i)
	{
		long double value;
		long double dist;
		size_t num,index;
		ifs>>num>>value>>dist>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoade(num,value,dist,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_qe;
	for(size_t i=1;i<=num_qe;++i)
	{
		long double q1,q2;
		long double d1,d2;
		size_t num,index;
		ifs>>num>>q1>>q2>>d1>>d2>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoade(num,q1,q2,d1,d2,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_fn;
	for(size_t i=1;i<=num_fn;++i)
	{
		long double value;
		size_t num,index;
		ifs>>num>>value>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoadn(num,value,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_bd;
	for(size_t i=1;i<=num_bd;++i)
	{
		long double value;
		size_t node_num,num;
		ifs>>node_num>>value>>num;
		num=num+fem_aufunc::FX-1;
		fem.AddBond(node_num,value,static_cast<fem_aufunc::INDEX>(num));
	}
}

inline void InputData(CFem & fem,ifstream & ifs)
{
	if(!ifs) return;
	size_t num_node,num_mat,num_real,num_elem,num_fe,num_qe,num_fn,num_bd;
	Cmatrix<long double> point(4,2);
	ifs>>num_node;
	for(size_t i=1;i<=num_node;i++)
	{
		size_t num;
		long double x,y,z;
		ifs>>num>>x>>y>>z;
		fem.AddNode(num,x,y,z);
	}
	ifs>>num_mat;
	for(size_t i=1;i<=num_mat;i++)
	{
		size_t num;
		string name;
		long double e,g;
		ifs>>num>>name>>e>>g;
		g=e/2/(1+g);
		fem.AddMat(CMaterial(num,name,e,g));
	}
	ifs>>num_real;
	for(size_t i=1;i<=num_real;i++)
	{
		size_t num;
		string name;
		long double area,iyy,izz,j;
		ifs>>num>>name>>area>>j>>iyy>>izz>>point(1,1)>>point(1,2)>>point(2,1)
			>>point(2,2)>>point(3,1)>>point(3,2)>>point(4,1)>>point(4,2);
		fem.AddReal(CReal3d(num,name,area,j,iyy,izz,point));
	}
	ifs>>num_elem;
	for(size_t i=1;i<=num_elem;i++)
	{
		size_t num;
		string name;
		size_t ni,nj,nmat,nreal;
		long double theta;
		ifs>>num>>name>>ni>>nj>>nmat>>nreal>>theta;
		fem.AddElem(CBeam3d(num,theta,name),ni,nj,nmat,nreal);
	}

	ifs>>num_fe;
	for(size_t i=1;i<=num_fe;++i)
	{
		long double value;
		long double dist;
		size_t num,index;
		ifs>>num>>value>>dist>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoade(num,value,dist,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_qe;
	for(size_t i=1;i<=num_qe;++i)
	{
		long double q1,q2;
		long double d1,d2;
		size_t num,index;
		ifs>>num>>q1>>q2>>d1>>d2>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoade(num,q1,q2,d1,d2,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_fn;
	for(size_t i=1;i<=num_fn;++i)
	{
		long double value;
		size_t num,index;
		ifs>>num>>value>>index;
		index=index+fem_aufunc::FX-1;
		fem.AddLoadn(num,value,static_cast<fem_aufunc::INDEX>(index));
	}

	ifs>>num_bd;
	for(size_t i=1;i<=num_bd;++i)
	{
		long double value;
		size_t node_num,num;
		ifs>>node_num>>value>>num;
		num=num+fem_aufunc::FX-1;
		fem.AddBond(node_num,value,static_cast<fem_aufunc::INDEX>(num));
	}
}

inline void OutputData(CFem &fem,ostream&os)
{
	os<<"������Ϣ����:\n";
	os<<fem.NodeCount()<<ends<<fem.ElemCount()<<ends<<fem.MatCount()<<ends<<fem.RealCount()<<endl;
	os<<endl<<"�����Ϣ����:\n";
	for(CFem::iterator_node it=fem.NodeHead();it!=fem.NodeEnd();++it)
	{
		os<<it->Num()<<ends<<it->x<<ends<<it->y<<ends<<it->z<<endl;
	}

	os<<endl<<"������Ϣ����:\n";
	for(CFem::const_iterator_mat it=fem.MatHead();it!=fem.MatEnd();++it)
	{
		Cmatrix<long double> info;
		(*it)->GetInfo(info);
		os<<(*it)->Num()<<ends<<(*it)->name<<ends<<info<<endl;
	}

	os<<endl<<"ʵ������Ϣ����:\n";
	for(CFem::const_iterator_real it=fem.RealHead();it!=fem.RealEnd();++it)
	{
		vector<pair<string,long double> > info;
		(*it)->GetInfo(info);
		os<<"���"<<ends<<(*it)->Num()<<endl<<"����"<<ends<<(*it)->name<<endl<<info<<endl<<endl;
	}

	os<<"��Ԫ��Ϣ����:\n";
	for(CFem::const_iterator_elem it=fem.ElemHead();it!=fem.ElemEnd();++it)
	{
		os<<(*it)->Num()<<ends<<(*it)->name<<ends<<(*it)->NumNodei()<<ends<<(*it)->NumNodej()<<ends<<(*it)->NumMat()<<ends<<(*it)->NumReal()<<endl;
	}

	os<<"\n��Ԫ���к�������:\n";
	for(CFem::const_iterator_Fe it=fem.FeHead();it!=fem.FeEnd();++it)
	{os<<it->GetNumElem()<<ends<<it->value<<ends<<it->dist<<ends<<fem_aufunc::FORCE[it->type]<<endl;}

	os<<"\n��Ԫ�ߺ�������:\n";
	for(CFem::const_iterator_qe it=fem.QeHead();it!=fem.QeEnd();++it)
	{os<<it->GetNumElem()<<ends<<it->q1<<ends<<it->q2<<ends<<it->disti<<ends<<it->distj<<ends<<fem_aufunc::FORCE[it->type]<<endl;}

	os<<"\n����������:\n";
	for(CFem::const_iterator_Fn it=fem.FnHead();it!=fem.FnEnd();++it)
	{os<<it->GetNumNode()<<ends<<it->value<<ends<<fem_aufunc::FORCE[it->type]<<endl;}

	os<<"\nԼ����Ϣ����:\n";
	for(CFem::const_iterator_bd it=fem.BondHead();it!=fem.BondEnd();++it)
	{
		CBond::bd_type rst;
		os<<it->Num()<<ends;
		it->GetBond(rst);
		for(CBond::bd_type::iterator _it=rst.begin();_it!=rst.end();++_it)
		{os<<_it->second<<ends<<fem_aufunc::FORCE[_it->first]<<",";}
		os<<endl;
	}
}

inline void OutPutRst(CFem &fem,ostream&os)
{
	os<<"���λ������:"<<endl;
	CFem::const_iterator_dofs it_dof=fem.NDofsHead();
	CFem::const_iterator_node it_n=fem.NodeHead();
	for(CFem::const_iterator_rst it=fem.NDispHead();it!=fem.NDispEnd();++it)
	{
		os<<"�ڵ�:"<<it_n->Num()<<endl;
		for(set<fem_aufunc::INDEX>::const_iterator it_d=it_dof->begin();it_d!=it_dof->end();++it_d)
			os<<fem_aufunc::MOMENT[*it_d]<<" ";
		os<<endl;
		for(size_t i=1;i<=it->line();++i)
			os<<(*it)(1,i)<<" ";
		os<<endl<<endl;
		++it_dof;
		++it_n;
	}

	os<<endl<<"��Ԫ��������:"<<endl;
	CFem::const_iterator_elem it_e=fem.ElemHead();
	vector<fem_aufunc::INDEX> rst;
	for(CFem::const_iterator_rst it=fem.EForceHead();it!=fem.EForceEnd();++it)
	{
		(*it_e)->RegDofs(rst);
		os<<"��Ԫ:"<<(*it_e)->Num()<<endl;
		for(size_t i=0;i<rst.size();++i)
			os<<fem_aufunc::FORCE[rst[i]]<<" ";
		os<<endl;
		os<<*it<<endl<<endl;
		++it_e;
	}

	os<<"\n֧�㷴������:"<<endl;
	os<<"�ڵ� ";
	for(size_t i=1;i<=fem.GetRForce()->line();++i)
		os<<fem_aufunc::MOMENT[i-1]<<" ";
	os<<endl;
	size_t i=1;
	for(CFem::const_iterator_bd it=fem.BondHead();it!=fem.BondEnd();++it)
	{
		os<<it->Num()<<ends;
		for(size_t j=1;j<=fem.GetRForce()->line();++j)
			os<<(*(fem.GetRForce()))(i,j)<<ends;
		os<<endl;
		++i;
		
	}
	
	it_e=fem.ElemHead();
	os<<"\n��ԪӦ������:"<<endl;
	for(CFem::const_iterator_rst it=fem.EStressHead();it!=fem.EStressEnd();++it)
	{
		os<<"��Ԫ:"<<(*it_e)->Num()<<endl;
		os<<*it<<endl<<endl;
		++it_e;
	}
}

inline void OutPutRst(CFem &fem,ostream &os,ostream &os2)  //���θָ��嵥Ԫ��������
{
	using namespace std;
	os<<"���λ������:"<<endl;
	CFem::const_iterator_dofs it_dof=fem.NDofsHead();
	CFem::const_iterator_node it_n=fem.NodeHead();
	os.precision(5);
	os.setf(ios::left);
	const size_t width=13;
	os<<"�ڵ�,UX,UY,UZ,ROTX,ROTY,ROTZ,WRAP"<<endl;
	for(CFem::const_iterator_rst it=fem.NDispHead();it!=fem.NDispEnd();++it)
	{
		os<<it_n->Num();
		for(size_t i=1;i<=it->line();++i)
			os<<","<<(*it)(1,i);
		os<<endl;
		++it_dof;
		++it_n;
	}

	os<<endl<<"��Ԫ��������:"<<endl;
	os<<"��Ԫ,i&j,FX,FY,FZ,MX,MY,MZ,B"<<endl;
	CFem::const_iterator_elem it_e=fem.ElemHead();
	for(CFem::const_iterator_rst it=fem.EForceHead();it!=fem.EForceEnd();++it)
	{
		os<<(*it_e)->Num()<<",i";
		for(size_t i=1;i<=it->line();++i)
			os<<","<<(*it)(1,i);
		os<<endl;
		os<<(*it_e)->Num()<<",j";
		for(size_t i=1;i<=it->line();++i)
			os<<","<<(*it)(2,i);
		os<<endl;
		++it_e;
	}

	os<<"\n֧�㷴������:"<<endl;
	os<<"�ڵ�";
	for(size_t i=1;i<=fem.GetRForce()->line();++i)
		os<<","<<fem_aufunc::MOMENT[i-1];
	os<<endl;
	size_t i=1;
	for(CFem::const_iterator_bd it=fem.BondHead();it!=fem.BondEnd();++it)
	{
		os<<it->Num();
		for(size_t j=1;j<=fem.GetRForce()->line();++j)
			os<<","<<(*(fem.GetRForce()))(i,j);
		os<<endl;
		++i;
	}
	
	it_e=fem.ElemHead();
	os2<<"��ԪӦ������:"<<endl;
	for(CFem::const_iterator_rst it=fem.EStressHead();it!=fem.EStressEnd();++it)
	{
		os2<<"��Ԫ:"<<(*it_e)->Num()<<endl;
		os2<<"���,y,z,��mi,��wi,��mi+��wi,��mj,��wj,��mj+��wj,��mi,��mj,��t,��tw,��tk+��tw,��i,��j"<<endl;
		for(size_t i=1;i<=it->row();++i)
		{
			os2<<(*it)(i,1)<<","<<(*it)(i,2)<<","<<(*it)(i,3);
			for(size_t j=4;j<=it->line();++j)
			{
				os2<<","<<(*it)(i,j);
			}
			os2<<endl;
		}
		++it_e;
	}
}

