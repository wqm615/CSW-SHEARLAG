#ifndef _CEQSET_
#define _CEQSET_
#include <string>
#include <vector>
#include "matrix\matrix.h"
#include "Mathfunc.h"
using namespace std;

struct Citem    //΢�ַ��̵���
{
	size_t var_index;  //�������,��1��ʼ��
	size_t order;      //��������
	Citem(size_t _index=0,size_t _ord=0):var_index(_index),order(_ord){};
};

struct CitemB:public Citem  //�߽緽�̵���
{
	size_t _index;  //��ɢ��ţ���1��ʼ
	double _coef;  //ϵ��
	CitemB(size_t VarIndex,size_t _ord,size_t index,double coef):Citem(VarIndex,_ord),_index(index),_coef(coef){}
	CitemB():_index(0),_coef(0){}
};

class Ceq  //΢�ַ���
{
public:
	vector<Citem> _equ;
private:
	Cmatrix<double> _coef;   //ϵ���Լ�ֵ�������һ��Ϊֵ
public:
	Ceq(){};
	Ceq(const vector<Citem>& equ,const Cmatrix<double> &coef):_equ(equ),_coef(coef){};
	void initilize(){_equ.clear();_coef.clear();}
	void AddItem(const Citem& _it){_equ.push_back(_it);}
	void SetCoef(const Cmatrix<double> & coef){_coef=coef;}
	void output(ostream & os) const
	{
		for(size_t i=0;i<_equ.size();++i)
			os<<_equ[i].var_index<<ends<<_equ[i].order<<endl;
		os<<_coef<<endl;
	}
	friend class CeqSet;
	friend class CeqSetD;
};

class CeqB   //�߽緽��
{
	vector<CitemB> _eqb;
	double _rvalue;
public:
	CeqB():_rvalue(0){}
	CeqB(const vector<CitemB>& eqb,double rvalue=0):_eqb(eqb),_rvalue(rvalue){};
	void initilize(){_eqb.clear(); _rvalue=0;}
	void AddItem(const CitemB& item){_eqb.push_back(item);}
	void SetValue(double value){_rvalue=value;}
	void output(ostream &os)
	{
		for(size_t i=0;i<_eqb.size();++i)
			os<<_eqb[i].var_index<<ends<<_eqb[i].order<<ends<<_eqb[i]._index<<ends<<_eqb[i]._coef<<",";
		os<<_rvalue<<endl;
	}
	friend class CeqSet;
};

class CeqSet
{
	double _delta;
	vector<CeqB> _eqbs;
protected:
	vector<Ceq> _eqs;
	Cmatrix<double> _mat;
	Cmatrix<double> _rval;
	Cmatrix<double> _rst;  //���
	Cmatrix<double> _rst_orig;
	vector<vector<double> >_rst_expand; //ԭʼ�⣬�������������ݵ㣬������������׵���
protected:
	vector<size_t> _VarDimon;  //�洢��������߽���
private:
	vector<size_t> _index;//������ʼ��ָʾ���飬����Ϊ������+1�����1��Ϊ�ܸյ�����
public:
	CeqSet(double delta=1.0):_delta(delta){};
	void initilize(){_eqs.clear();_VarDimon.clear();_mat.clear();_rst_expand.clear();
	_delta=1;_index.clear();_eqbs.clear();_rval.clear();_rst.clear();_rst_orig.clear();}
	inline void AddEq(const Ceq& _eq){_eqs.push_back(_eq);}
	inline void AddEqb(const CeqB& _eqb){_eqbs.push_back(_eqb);}
	inline void SetDelta(double delta){_delta=delta;}
	inline void OutputAll(ostream &);
	inline void output(ostream &);
	inline void GetRst(Cmatrix<double>& rst,size_t pos_dest,size_t pos_orig=1) const
	{
		rst.assign(1,pos_dest,_rst.row(),pos_dest+_VarDimon.size()-1,_rst,1,pos_orig);
	}
private:
	inline void CalcVarDimon(); //���������������
private:
	inline void CalcMat();
public:
	inline void Calculate();
public:
	inline bool CalcVarRst(size_t var_index,size_t ord,Cmatrix<double>&rst);  // ����ĳ��������׵���
	inline double CalcVarRst(size_t var_index,size_t ord,size_t index);
	inline double check();
};

void CeqSet::output(ostream & os)
{os<<_rst<<endl;}

void CeqSet::OutputAll(ostream& os)
{
	for(size_t i=0;i<_eqs.size();++i)
	{
		os<<"����"<<i+1<<":\n";
		_eqs[i].output(os);
		os<<endl;
	}
	os<<"�߽�����:\n";
	for(size_t i=0;i<_eqbs.size();++i) _eqbs[i].output(os);
	os<<"\n��������:\n";
	for(size_t i=0;i<_VarDimon.size();++i)
		os<<_VarDimon[i]<<ends;
	os<<endl;
	os<<"\n�����ھ����е�λ��:\n";
	for(size_t i=0;i<_index.size();++i)
		os<<_index[i]<<ends;
	os<<endl<<"\n����:"<<endl;
	os<<_mat;
	os<<endl<<"\nֵ����:"<<endl;
	os<<_rval;
	os<<endl<<"\nԭʼ�⣺"<<endl<<_rst_orig;
	os<<endl<<"\n��չ��:"<<endl;
	for(size_t i=0;i<_rst_expand.size();++i)
	{
		for(size_t j=0;j<_rst_expand[i].size();++j)
		{
			os<<_rst_expand[i][j]<<ends;
		}
		os<<endl;
	}
	os<<endl<<"�⣺"<<endl;
	os<<_rst;
	os<<endl<<"У��ֵ��"<<check();
}

void CeqSet::CalcVarDimon()
{
	size_t max_index=0;   //��������
	for(size_t i=0;i<_eqs.size();++i)
	{
		for(size_t j=0;j<_eqs[i]._equ.size();++j)
		{
			max_index>_eqs[i]._equ[j].var_index?max_index:max_index=_eqs[i]._equ[j].var_index;
			
		}
	}
	this->_VarDimon.resize(max_index);
	for(size_t i=0;i<_eqs.size();++i)
	{
		for(size_t j=0;j<_eqs[i]._equ.size();++j)
		{
			if(_VarDimon[_eqs[i]._equ[j].var_index-1]<_eqs[i]._equ[j].order) _VarDimon[_eqs[i]._equ[j].var_index-1]=_eqs[i]._equ[j].order;		
		}
	}
	_index.resize(_VarDimon.size()+1);
	_index[0]=1;
	const size_t VarNum=_eqs[0]._coef.row();  //ע��ÿ�����̵���ɢ����������ͬ
	for(size_t i=0;i<_VarDimon.size();++i)
	{
		_index[i+1]=_VarDimon[i]+VarNum+_index[i];
	}
}

void CeqSet::CalcMat()
{
	_mat.clear();
	_mat.resize(_index.back()-1,_index.back()-1,0);
	vector<double> buffer;  //���ڴ洢��ֵ�ϵ��
	const size_t VarNum=_eqs[0]._coef.row();
	//΢�ַ���
	for(size_t i=0;i<_eqs.size();++i)
	{   //����
		for(size_t j=0;j<_eqs[i]._equ.size();++j)
		{   //��
			Diff(_eqs[i]._equ[j].order,buffer);
			for(size_t k=0;k<VarNum;++k)
			{   //��ɢ��
				for(size_t kk=0;kk<buffer.size();++kk)
					_mat(i*VarNum+k+1,_index[_eqs[i]._equ[j].var_index-1]+k+kk)+=buffer[kk]/pow(_delta,int(_eqs[i]._equ[j].order))*_eqs[i]._coef(k+1,j+1);
			}
		}
	}
	//�߽�����
	size_t pos=VarNum*_eqs.size()+1;
	for(size_t i=0;i<_eqbs.size();++i)
	{
		//����
		for(size_t j=0;j<_eqbs[i]._eqb.size();++j)
		{
			//��
			Diff(_eqbs[i]._eqb[j].order,buffer);
			for(size_t k=0;k<buffer.size();++k)
			{
				_mat(pos+i,_index[_eqbs[i]._eqb[j].var_index-1]+_eqbs[i]._eqb[j]._index-1+k)+=buffer[k]/pow(_delta,int(_eqbs[i]._eqb[j].order))*_eqbs[i]._eqb[j]._coef;
			}
		}
	}
	//�����ұ�
	_rval.resize(_mat.row(),1);
	pos=0;
	for(size_t i=0;i<_eqs.size();++i)
	{
		for(size_t j=0;j<VarNum;++j)
		{
			++pos;
			_rval(pos,1)=_eqs[i]._coef(j+1,_eqs[i]._coef.line());	
		}
	}
	for(size_t i=0;i<_eqbs.size();++i)
	{
		_rval(++pos,1)=_eqbs[i]._rvalue;
	}
}

void CeqSet::Calculate()
{
	CalcVarDimon();
	CalcMat();
	Cmatrix<double> tmp(_mat);
	tmp.solve_dolt(_rval,_rst_orig);
	_rst.resize(_eqs[0]._coef.row(),_VarDimon.size());
	_rst_expand.resize(_VarDimon.size());
	for(size_t j=1;j<_index.size();++j)
	{
		for(size_t i=1;i<=_rst.row();++i)
		{
			_rst(i,j)=_rst_orig(_index[j-1]+i-1,1);
		}
	}
	size_t pos=0;
	for(size_t i=0;i<_rst_expand.size();++i)
	{
		_rst_expand[i].resize(_eqs[0]._coef.row()+_VarDimon[i]);
		for(size_t j=0;j<_rst_expand[i].size();++j)
		{
			_rst_expand[i][j]=_rst_orig(++pos,1);
		}
	}
}

double CeqSet::CalcVarRst(size_t var_index,size_t ord,size_t index)
{
	vector<double> buffer;
	double rtn=0;
	Diff(ord,buffer);
	for(size_t j=0;j<buffer.size();++j)
	{
			rtn+=buffer[j]*_rst_expand[var_index-1][index+j-1]/pow(_delta,int(ord));
	}
	return rtn;
}

bool CeqSet::CalcVarRst(size_t var_index,size_t ord,Cmatrix<double>&rst)
{
	rst.clear();
	if(ord>_VarDimon[var_index-1]) 
	{
		return false;
	}
	vector<double> buffer;
	Diff(ord,buffer);
	rst.resize(_eqs[0]._coef.row(),1,0);
	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,1)=CalcVarRst(var_index,ord,i);
	}
	return true;
}

double CeqSet::check()
{
	Cmatrix<double> buffer1; 
	double rtn=0,rtn_d=0;
	for(size_t i=0;i<_eqs.size();++i)
	{
		buffer1.resize(_eqs[i]._coef.row(),_eqs[i]._coef.line()-1);
		Cmatrix<double> buffer2(_eqs[i]._coef.row(),_eqs[i]._coef.line());
		for(size_t j=0;j<_eqs[i]._equ.size();++j)
		{
			CalcVarRst(_eqs[i]._equ[j].var_index,_eqs[i]._equ[j].order,buffer2);
			buffer1.assign(1,j+1,buffer1.row(),j+1,buffer2);
		}
		for(size_t j=1;j<=buffer1.row();++j)
		{
			for(size_t k=1;k<=buffer1.line();++k)
				rtn_d+=buffer1(j,k)*_eqs[i]._coef(j,k);
			rtn_d-=_eqs[i]._coef(j,_eqs[i]._coef.line());
			rtn+=abs(rtn_d);
			rtn_d=0;
		}
	}
	for(size_t i=0;i<_eqbs.size();++i)
	{
		for(size_t j=0;j<_eqbs[i]._eqb.size();++j)
		{
			rtn_d+=CalcVarRst(_eqbs[i]._eqb[j].var_index,_eqbs[i]._eqb[j].order,_eqbs[i]._eqb[j]._index)*_eqbs[i]._eqb[j]._coef;
		}
		/*cout<<"\nrtn_d="<<rtn_d;*/
		rtn_d-=_eqbs[i]._rvalue;
		rtn+=abs(rtn_d);
		rtn_d=0;
	}
	return rtn;
}

#endif