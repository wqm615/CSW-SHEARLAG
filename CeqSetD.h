#ifndef _CEQSETD_
#define _CEQSETD_
#include "DIM.h"
#include <fem\fem.h>
struct CitemBD:public CitemB
{
	size_t SecNum;    //从1开始
	CitemBD(size_t sec_num,size_t VarIndex,size_t _ord,size_t index,double coef)
		:CitemB(VarIndex,_ord,index,coef),SecNum(sec_num){};
	CitemBD():SecNum(0){};
};

class CeqBD
{
	vector<CitemBD> _eqbd;
	double _rvalue;
public:
	CeqBD():_rvalue(0){}
	CeqBD(const vector<CitemBD>& eqbd,double rvalue=0):_eqbd(eqbd),_rvalue(rvalue){};
	void initilize(){_eqbd.clear(); _rvalue=0;}
	void AddItem(const CitemBD& item){_eqbd.push_back(item);}
	void SetValue(double value){_rvalue=value;}
	void output(ostream &os) const
	{
		for(size_t i=0;i<_eqbd.size();++i)
			os<<_eqbd[i].SecNum<<ends<<_eqbd[i].var_index<<ends<<_eqbd[i].order
			<<ends<<_eqbd[i]._index<<ends<<_eqbd[i]._coef<<",";
		os<<_rvalue<<endl;
	}
	friend class CeqSetD;
};

class CeqSetD: public CeqSet
{
public:
	vector<CeqBD> _eqbds;

public:
	Cmatrix <Cmatrix<double> > _coefs;    //行号为方程id，列号为分段id
	vector<double> _deltas;    //存储每个分段的步长
	vector<size_t> SecDataCount;   //每一分段的离散数据点数
	Cmatrix<size_t> Lindex;        //变量在总刚中列位置的起始点指示矩阵，从1开始计算，行号为变量id，列号为分段id,
	                               //依分段号、变量号、离散点号排序
	Cmatrix<size_t> Rindex;        //变量在总刚中行位置的起始点指示矩阵，从1开始计算，行号为方程id，列号为分段id
public:
	CeqSetD(){};
	void initilize(){SecDataCount.clear();_coefs.clear();_deltas.clear();
	CeqSet::initilize();Lindex.clear();Rindex.clear();_eqbds.clear();}
public:
	inline void SetDelta(const vector<double> &delta)
	{
		_deltas=delta;
	}
	inline void SetSecNum(const vector<size_t>& sec_data_count){SecDataCount=sec_data_count;}
	inline void ReSizeCoef(){_coefs.resize(_eqs.size(),SecDataCount.size());}
	inline void InputCoef(istream& is);
	inline void SetCoef(const Cmatrix<double>& coef,size_t eqs_id,size_t sec_id);
	inline void SetCoefRvalue(vector<double>& _rvalue,size_t eqs_id,size_t sec_id);
	inline Cmatrix <Cmatrix<double> >* GetCoefs(){return &_coefs;}
	inline void AddEqb(const CeqBD& eqbd){_eqbds.push_back(eqbd);}
public:
	inline void CalcVarDimon(); //计算变量的最大阶数
	inline void CalcMat();
	inline void Calculate();     
public:
	inline void OutputAll(ostream &) const;
public:
	inline bool CheckData();   //计算_delta、_coefs与SecNum、eqs.size()是否匹配
public:
	inline bool CalcVarRst(size_t var_index,size_t ord,Cmatrix<double>&rst) const;  // 计算某变量任意阶导数
	inline bool CalcVarRst(size_t ord,Cmatrix<double>&rst) const;  // 计算所有变量任意阶导数
	inline bool CalcVarRst(size_t var_index,size_t ord,Cmatrix<double>&rst,size_t secnum) const;  // 计算某变量任意阶导数
	inline double CalcVarRst(size_t var_index,size_t ord,size_t index,size_t secnum) const;
	inline double check() const;
public:
	inline void Reverse(CeqSetD &eqsd) const;   //将差分数据倒序
};

bool CeqSetD::CheckData()
{
	if(_deltas.size()!=SecDataCount.size())
	{
		cout<<"_deltas分段数目与SecNum不相等!\n";
		return false;
	}
	if(_coefs.line()!=SecDataCount.size())
	{
		cout<<"_coefs列数与SecNum不相等!\n";
		return false;
	}
	if(_coefs.row()!=_eqs.size())
	{
		cout<<"_coefs行数与方程数量不相等!\n";
		return false;
	}
	for(size_t i=1;i<=_coefs.row();++i)  //方程号
	{
		for(size_t j=1;j<=_coefs.line();++j)  // 分段号
		{
			if (_coefs(i,j).line()!=_eqs[i-1]._equ.size()+1)
			{
				cout<<"第"<<i<<"方程第"<<j<<"分段的系数矩阵列数不对!\n";
				return false;
			}
			if (_coefs(i,j).row()!=SecDataCount[j-1])
			{
				cout<<"第"<<i<<"方程第"<<j<<"分段的数量不对!\n";
				return false;
			}
		}
	}
	return true;
}

void CeqSetD::InputCoef(istream& is)
{
	_coefs.resize(_eqs.size(),SecDataCount.size());
	for(size_t i=1;i<=_coefs.row();++i)
	{
		for(size_t j=1;j<=_coefs.line();++j)
		{
			_coefs(i,j).resize(SecDataCount[j-1],_eqs[i-1]._equ.size()+1);
			is>>_coefs(i,j);
		}
	}
}

void CeqSetD::SetCoef(const Cmatrix<double>& coef,size_t eqs_id,size_t sec_id)
{
	_coefs(eqs_id,sec_id)=coef;
}

void CeqSetD::SetCoefRvalue(vector<double> &_rvalue,size_t eqs_id,size_t sec_id)
{
	for(size_t i=1;i<=_coefs(eqs_id,sec_id).row();++i)
		_coefs(eqs_id,sec_id)(i,8)=_rvalue[i-1];
}

void CeqSetD::OutputAll(ostream & os) const
{
	for(size_t i=0;i<_eqs.size();++i)
	{
		os<<"方程"<<i+1<<":\n";
		_eqs[i].output(os);
	}
	os<<"方程系数:\n";
	for(size_t i=1;i<=SecDataCount.size();++i)
	{
		os<<"分段"<<i<<":\n";
		for(size_t j=1;j<=_eqs.size();++j)
		{
			os<<"方程"<<j<<":\n";
			os<<_coefs(j,i)<<endl;
		}
		os<<endl;
	}
	os<<"边界条件:\n";
	for(size_t i=0;i<_eqbds.size();++i)
	{
		os<<"边界方程"<<i+1<<endl;
		_eqbds[i].output(os);
		os<<endl;
	}
	os<<"\n总刚：\n"<<_mat<<endl;
	os<<"\n原始解向量:\n"<<_rst_orig<<endl;
	os<<"\n解:\n"<<_rst<<endl;
	os<<"\n扩展解\n";
	for(size_t i=0;i<_VarDimon.size();++i)
	{
		for(size_t j=0;j<_rst_expand[i].size();++j)
			os<<_rst_expand[i][j]<<ends;
		os<<endl;
	}
	os<<"\n校验值:"<<check()<<endl;
}

void CeqSetD::CalcVarDimon()
{
	size_t max_index=0;   //变量个数
	for(size_t i=0;i<_eqs.size();++i)
	{
		for(size_t j=0;j<_eqs[i]._equ.size();++j)
		{
			max_index>_eqs[i]._equ[j].var_index?max_index:max_index=_eqs[i]._equ[j].var_index;
			
		}
	}
	_VarDimon.resize(max_index);
	for(size_t i=0;i<_eqs.size();++i)
	{
		for(size_t j=0;j<_eqs[i]._equ.size();++j)
		{
			if(_VarDimon[_eqs[i]._equ[j].var_index-1]<_eqs[i]._equ[j].order) _VarDimon[_eqs[i]._equ[j].var_index-1]=_eqs[i]._equ[j].order;		
		}
	}
	size_t tmp=1;
	Lindex.resize(_VarDimon.size(),SecDataCount.size());
	for(size_t j=1;j<=SecDataCount.size();++j)
	{
		for(size_t i=1;i<=_VarDimon.size();++i)
		{
			Lindex(i,j)=tmp;
			tmp=Lindex(i,j)+SecDataCount[j-1]+_VarDimon[i-1];
		}
	}
	Rindex.resize(_eqs.size(),SecDataCount.size());
	tmp=1;
	for(size_t j=1;j<=SecDataCount.size();++j)
	{
		for(size_t i=1;i<=_eqs.size();++i)
		{
			Rindex(i,j)=tmp;
			tmp+=SecDataCount[j-1];
		}
	}
}

void CeqSetD::CalcMat()
{
	_mat.clear();
	_mat.resize(Lindex(Lindex.row(),Lindex.line())+_VarDimon.back()+SecDataCount.back()-1
		,Lindex(Lindex.row(),Lindex.line())+_VarDimon.back()+SecDataCount.back()-1,0);
	vector<double> buffer;  //用于存储差分的系数
	//微分方程
	for(size_t m=1;m<=SecDataCount.size();++m)
	{
		const size_t VarNum=SecDataCount[m-1];
		for(size_t i=0;i<_eqs.size();++i)
		{   //方程
			for(size_t j=0;j<_eqs[i]._equ.size();++j)
			{   //项
				Diff(_eqs[i]._equ[j].order,buffer);
				for(size_t k=0;k<VarNum;++k)
				{   //离散点
					for(size_t kk=0;kk<buffer.size();++kk)
					{
						_mat(Rindex(i+1,m)+k,Lindex(_eqs[i]._equ[j].var_index,m)+k+kk)
							+=buffer[kk]/pow(_deltas[m-1],int(_eqs[i]._equ[j].order))*_coefs(i+1,m)(k+1,j+1);
					}
				}
			}
		}
	}
	
	//边界条件
	size_t pos=Rindex(_eqs.size(),SecDataCount.size())+SecDataCount.back();
	for(size_t i=0;i<_eqbds.size();++i)
	{
		//方程
		for(size_t j=0;j<_eqbds[i]._eqbd.size();++j)
		{
			//项
			Diff(_eqbds[i]._eqbd[j].order,buffer);
			for(size_t k=0;k<buffer.size();++k)
			{
				_mat(pos+i,Lindex(_eqbds[i]._eqbd[j].var_index,_eqbds[i]._eqbd[j].SecNum)+_eqbds[i]._eqbd[j]._index+k-1)
					+=buffer[k]/pow(_deltas[_eqbds[i]._eqbd[j].SecNum-1],int(_eqbds[i]._eqbd[j].order))*_eqbds[i]._eqbd[j]._coef;
			}
		}
	}
	
	//方程右边
	_rval.resize(_mat.row(),1);
	pos=0;
	for(size_t k=0;k<SecDataCount.size();++k)
	{
		for(size_t i=0;i<_eqs.size();++i)
		{
			for(size_t j=0;j<SecDataCount[k];++j)
			{
				++pos;
				_rval(pos,1)=_coefs(i+1,k+1)(j+1,_coefs(i+1,k+1).line());	
			}
		}
	}
	for(size_t i=0;i<_eqbds.size();++i)
	{
		_rval(++pos,1)=_eqbds[i]._rvalue;
	}

}

void CeqSetD::Calculate()
{
	if(!CheckData()) throw;
	CalcVarDimon();
	CalcMat();
	Cmatrix<double> tmp(_mat);
	tmp.solve_dolt(_rval,_rst_orig);
	size_t DataCount=0;
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		DataCount+=SecDataCount[i];
	}
	_rst.resize(DataCount,_VarDimon.size());
	for(size_t i=0;i<_VarDimon.size();++i)
	{
		size_t pos=0;
		for(size_t j=0;j<SecDataCount.size();++j)
		{
			for(size_t k=1;k<=SecDataCount[j];++k)
				_rst(++pos,i+1)=_rst_orig(Lindex(i+1,j+1)+k-1,1);
		}
		
	}
	_rst_expand.resize(_VarDimon.size());
	for(size_t i=0;i<_rst_expand.size();++i)
	{
		size_t pos=0;
		_rst_expand[i].resize(_rst.row()+_VarDimon[i]*SecDataCount.size());
		for(size_t j=0;j<SecDataCount.size();++j)
		{
			for(size_t k=1;k<=SecDataCount[j]+_VarDimon[i];++k)
				_rst_expand[i][pos++]=_rst_orig(Lindex(i+1,j+1)+k-1,1);
		}
	}
}

double CeqSetD::CalcVarRst(size_t var_index,size_t ord,size_t index,size_t secnum) const
{
	if(index>SecDataCount[secnum-1])
	{
		cout<<"CalcVarRst:index越界!\n";
		throw;
	}
	vector<double> buffer;
	double rtn=0;
	Diff(ord,buffer);
	for(size_t j=0;j<buffer.size();++j)
	{
			rtn+=buffer[j]*_rst_orig(Lindex(var_index,secnum)+index-1+j,1)/pow(_deltas[secnum-1],int(ord));
	}
	return rtn;
}

bool CeqSetD::CalcVarRst(size_t var_index,size_t ord,Cmatrix<double>&rst) const
{
	rst.clear();
	if(ord>_VarDimon[var_index-1]) 
	{
		return false;
	}
	vector<double> buffer;
	Diff(ord,buffer);
	rst.resize(_rst.row(),1,0);
	size_t pos=0;
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		for(size_t j=1;j<=SecDataCount[i];++j)
		rst(++pos,1)=CalcVarRst(var_index,ord,j,i+1);
	}
	return true;
}

bool CeqSetD::CalcVarRst(size_t ord,Cmatrix<double>&rst) const
{
	rst.resize(_rst.row(),_VarDimon.size(),0);
	vector<double> buffer;
	for(size_t var_index=1;var_index<=_VarDimon.size();++var_index)
	{
		if(ord>_VarDimon[var_index-1]) 
		{
			return false;
		}
		Diff(ord,buffer);
		size_t pos=0;
		for(size_t i=0;i<SecDataCount.size();++i)
		{
			for(size_t j=1;j<=SecDataCount[i];++j)
			rst(++pos,var_index)=CalcVarRst(var_index,ord,j,i+1);
		}
	}
	return true;
}

bool CeqSetD::CalcVarRst(size_t var_index,size_t ord,Cmatrix<double>&rst,size_t secnum) const
{
	rst.clear();
	if(ord>_VarDimon[var_index-1]) 
	{
		return false;
	}
	vector<double> buffer;
	Diff(ord,buffer);
	rst.resize(SecDataCount[secnum-1],1,0);
	size_t pos=0;
	for(size_t j=1;j<=SecDataCount[secnum-1];++j)
		rst(++pos,1)=CalcVarRst(var_index,ord,j,secnum);
	return true;
}

double CeqSetD::check() const
{
	Cmatrix<double> buffer1; 
	double rtn=0,rtn_d=0;
	for(size_t m=0;m<SecDataCount.size();++m)
	{
		for(size_t i=0;i<_eqs.size();++i)
		{
			buffer1.resize(_coefs(i+1,m+1).row(),_coefs(i+1,m+1).line()-1);
			Cmatrix<double> buffer2;
			for(size_t j=0;j<_eqs[i]._equ.size();++j)
			{
				CalcVarRst(_eqs[i]._equ[j].var_index,_eqs[i]._equ[j].order,buffer2,m+1);
				buffer1.assign(1,j+1,buffer1.row(),j+1,buffer2);
			}
			for(size_t j=1;j<=buffer1.row();++j)
			{
				for(size_t k=1;k<=buffer1.line();++k)
					rtn_d+=buffer1(j,k)*_coefs(i+1,m+1)(j,k);
				rtn_d-=_coefs(i+1,m+1)(j,_coefs(i+1,m+1).line());
				rtn+=abs(rtn_d);
				rtn_d=0;
			}
		}
	}
	for(size_t i=0;i<_eqbds.size();++i)
	{
		for(size_t j=0;j<_eqbds[i]._eqbd.size();++j)
		{
			rtn_d+=CalcVarRst(_eqbds[i]._eqbd[j].var_index,_eqbds[i]._eqbd[j].order,_eqbds[i]._eqbd[j]._index,
				_eqbds[i]._eqbd[j].SecNum)*_eqbds[i]._eqbd[j]._coef;
		}
		rtn_d-=_eqbds[i]._rvalue;
		rtn+=abs(rtn_d);
		rtn_d=0;
	}
	return rtn;
}

void CeqSetD::Reverse(CeqSetD &eqsd) const
{
	eqsd.initilize();
	eqsd._eqs=_eqs;
	eqsd.SecDataCount.assign(SecDataCount.rbegin(),SecDataCount.rend());
	eqsd._deltas.assign(_deltas.rbegin(),_deltas.rend());
	eqsd._coefs=_coefs;
	eqsd._coefs.revered(ROW);
	for(size_t i=1;i<=eqsd._coefs.row();++i)
	{
		for(size_t j=1;j<=eqsd._coefs.line();++j)
		{
			eqsd._coefs(i,j).revered();
		}
	}
	//将奇数阶系数乘-1
	for(size_t i=0;i<eqsd._eqs.size();++i)
	{
		for(size_t j=0;j<eqsd._eqs[i]._equ.size();++j)
		{
			if(eqsd._eqs[i]._equ[j].order%2==1)
			{
				for(size_t k=0;k<eqsd.SecDataCount.size();++k)
				{
					for(size_t kk=1;kk<=eqsd.SecDataCount[k];++kk)
					{
						eqsd._coefs(i+1,k+1)(kk,j+1)*=-1;
					}
				}
			}
		}
	}
	eqsd._eqbds.resize(_eqbds.size());
	{
		size_t j=_eqbds.size()-1;
		for(size_t i=0;i<_eqbds.size();++i)
		{
			eqsd._eqbds[j]._rvalue=_eqbds[i]._rvalue;
			eqsd._eqbds[j]._eqbd.resize(_eqbds[i]._eqbd.size());
			for(size_t k=0;k<_eqbds[i]._eqbd.size();++k)
			{
				eqsd._eqbds[j]._eqbd[k].SecNum=SecDataCount.size()+1-_eqbds[i]._eqbd[k].SecNum;
				eqsd._eqbds[j]._eqbd[k].order=_eqbds[i]._eqbd[k].order;
				eqsd._eqbds[j]._eqbd[k].var_index=_eqbds[i]._eqbd[k].var_index;
				if(eqsd._eqbds[j]._eqbd[k].order%2==0)	eqsd._eqbds[j]._eqbd[k]._coef=_eqbds[i]._eqbd[k]._coef;
				else eqsd._eqbds[j]._eqbd[k]._coef=-_eqbds[i]._eqbd[k]._coef;
				eqsd._eqbds[j]._eqbd[k]._index=SecDataCount[_eqbds[i]._eqbd[k].SecNum-1]+1-_eqbds[i]._eqbd[k]._index;
			}
			--j;	
		}
	}
};

#endif