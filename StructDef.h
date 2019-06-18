#ifndef _STRUCTDEF_
#define _STRUCTDEF_
struct CUnifLoad
{
	double x1,x2,q1,q2;  //单位m,kN/m
	double Calcq(const double x) const
	{
		if(x<x1||x>x2) return 0;
		return q1+(q2-q1)/(x2-x1)*(x-x1);
	}
};

struct CContLoad
{
	double x,f;   //单位m,kN/m
	CContLoad(double _x=0,double _f=0):x(_x),f(_f){};
};

struct CBd
{
	double x;
	size_t Fd;
};

struct CPlateProp
{
	double A,I,C;
	CPlateProp(double _A=0,double _I=0,double _C=0):A(_A),I(_I),C(_C){};
	void output(ostream &os)
	{
		os<<"A:"<<A<<",I:"<<I<<",C:"<<C<<endl;
	}
};

struct CTopPlate
{
	double A,I,C;
	CTopPlate(double _A=0,double _I=0,double _C=0):A(_A),I(_I),C(_C){};
};

struct CBotPlate:public CTopPlate
{
	CBotPlate(double _A=0,double _I=0,double _C=0):CTopPlate(_A,_I,_C){};
	CBotPlate(const CTopPlate&_TopPlate):CTopPlate(_TopPlate){};
	CBotPlate GetBotPlate(const double _B,const double _hd)
	{
		CBotPlate rtn;
		rtn.A=A+_B*_hd;
		rtn.C=(A*(C+_hd)+_B*_hd*_hd/2)/rtn.A;
		rtn.I=I+A*pow(rtn.C-C-_hd,2)+_B*pow(_hd,3)/12+_B*_hd*pow(rtn.C-_hd/2,2);
		return rtn;
	}
};

struct CSection
{
	double H;
	CBotPlate botplate;
	CSection(const CBotPlate & _botplate,double _H):botplate(_botplate),H(_H){};
	CSection(){H=0;}
	pair<double,double> CalcPropety(const CTopPlate& TopPlate1,const CTopPlate &TopPlate2)
	{
		pair<double,double> rtn;
		rtn.first=TopPlate1.A+TopPlate2.A+botplate.A;
		double C=(TopPlate1.A*(H-TopPlate1.C)+TopPlate2.A*(H-TopPlate2.C)+botplate.A*botplate.C)/rtn.first;
		rtn.second=TopPlate1.I+TopPlate1.A*pow(H-TopPlate1.C-C,2)+TopPlate2.I+TopPlate2.A*pow(H-TopPlate2.C-C,2)
			+botplate.I+botplate.A*pow(C-botplate.C,2);
		return rtn;
	}
	void CalcI(const CTopPlate& TopPlate1,const CTopPlate &TopPlate2,double &I1,double &I2,double &I3)
	{
		double A=TopPlate1.A+TopPlate2.A+botplate.A;
		double C=(TopPlate1.A*(H-TopPlate1.C)+TopPlate2.A*(H-TopPlate2.C)+botplate.A*botplate.C)/A;
		I1=TopPlate1.I+TopPlate1.A*pow(H-TopPlate1.C-C,2);
		I2=TopPlate2.I+TopPlate2.A*pow(H-TopPlate2.C-C,2);
		I3=botplate.I+botplate.A*pow(C-botplate.C,2);
	}
};

class Ccurve
{
	//double vals,valc;  //支点、跨中数值
public:
	vector<pair<double,size_t> > keypoint;
	double pala;  //次数
public:
	Ccurve(double _pala,const vector<pair<double,size_t> >& _keypoint):
	pala(_pala),keypoint(_keypoint){};
	Ccurve(double _pala=0):pala(_pala){};
	void clear(){pala=0;keypoint.clear();}
	void setpala(double _pala)
	{pala=_pala;}
	void setkeypoint(const vector<pair<double,size_t> >& _keypoint){keypoint=_keypoint;}
	void addkeypoint(const pair<double,size_t>& _keypoint){keypoint.push_back(_keypoint);}
public:
	void output(ostream & os)
	{
		os<<"次数:"<<ends<<pala<<endl<<"曲线"<<endl;
		for(size_t i=0;i<keypoint.size();++i)
			os<<keypoint[i].first<<","<<keypoint[i].second<<endl;
	}
	double operator()(double valc,double vals,double x)const
	{
		vector<double> index(2);
		index[0]=valc;
		index[1]=vals;
		if(x<keypoint.front().first||x>keypoint.back().first) 
		{	
			cout<<"Ccurve::operator()参数x有误!";
			return -1;
		}
		size_t pos=0;
		for(;pos<keypoint.size();++pos)
		{
			if(x==keypoint[pos].first) 
			{
				return index[keypoint[pos].second];
			}
			if(x<keypoint[pos].first) break;
		}
		--pos;
		if(keypoint[pos].second==keypoint[pos+1].second) return index[keypoint[pos].second];
		const double delta=vals-valc;
		const double L=keypoint[pos+1].first-keypoint[pos].first;
		double dist;
		if(!keypoint[pos].second) dist=x-keypoint[pos].first;
		else dist=keypoint[pos+1].first-x;
		return delta/pow(L,pala)*pow(dist,pala)+valc;
	}
};

struct I_dI  
{
	double I1,I2,I3,Iw,I,Is,C,H;
	double dI;
};

struct Cx_I_dI:public I_dI   //仅functionV中使用该类
{
	double x;
	Cx_I_dI(double _x=0):x(_x){};
};

struct Ckeypos
{
	double x;
	int flag;  //1表示变高度梁分段位置，2表示隔板位置，0为其它情况
	Ckeypos(double _x=0,int _flag=0):x(_x),flag(_flag){};
};

bool operator<(const Ckeypos&v1,const Ckeypos&v2)
{
	if(abs(v1.x-v2.x)<=1e-6) return false;	
	else return v1.x<v2.x;
}

struct CPlateGeom
{
	double x;
	double thickness;
	CPlateGeom(double _x=0,double _thickness=0):x(_x),thickness(_thickness){};
};

struct CSecPropty
{
	double H,C,A,I1,I2,I3,Iw,I;  //C为截面形心至截面底缘的距离
	CSecPropty():H(0),C(0),A(0),I1(0),I2(0),I3(0),Iw(0),I(0){};
};

struct CSecPropty2:public CSecPropty
{
	double x,dI,Is;
	double dI1,dI2,dI3,dIw;
	CSecPropty2():CSecPropty(),x(0),dI(0),Is(0){};
	void Calculate() {Is=I1+I2+I3;I=Is+Iw;dI=dI1+dI2+dI3+dIw;}
	void assign(const CSecPropty & sec_proty)
	{
		H=sec_proty.H;
		A=sec_proty.A;
		C=sec_proty.C;
		I1=sec_proty.I1;
		I2=sec_proty.I2;
		I3=sec_proty.I3;
		Iw=sec_proty.Iw;
		Calculate();	
	}
};

class VarSecBeam
{
public:
	double H_delta/*底板厚度变化*/,Hz,Hc;
public:
	vector<CPlateGeom> plate1,plate2,web_top; 
	vector<CPlateGeom> plate3,web_bot;
	vector<CPlateGeom> plate3_zd,web_bot_zd;  //支点处的板3、底板腹板
public:
	Ccurve curve;
public:
	CPlateProp propety1,propety2,propety_web_top;  //保存顶板的特性
	CPlateProp propety3,propety_web_bot;  //保存底板的特性
public:
	CSecPropty sec_propty;
private:
	double H_top_max;  //顶板腹板位置的厚度,用于计算波形钢腹板高度
public:
	VarSecBeam(){};
	VarSecBeam(double h_delta,vector<CPlateGeom>&_plate1,vector<CPlateGeom>&_plate2,
		vector<CPlateGeom>&_plate3_zd,vector<CPlateGeom>&_web_top,vector<CPlateGeom>&_web_bot)
		:H_delta(h_delta),plate1(_plate1),plate2(_plate2),plate3_zd(_plate3_zd),web_top(_web_top)
	,web_bot(_web_bot){PreCalculate();};
public:
	inline void PreCalculate();  //计算两块顶板的实常数以及顶板腹板位置的厚度
	inline bool Calculate(const double x,CSecPropty2& SecProp);
	inline bool Calculate(const double x);
	inline double Calculate_dI(const double x,const double dx=1e-5);
	inline void Calculate_dI(const double x,CSecPropty2& SecProp,const double dx=1e-5);
	inline double CalculateHw(const double x);  //计算腹板高度
private:
	inline void CalcProp(const vector<CPlateGeom>& plate,CPlateProp& propety);
public:
	void output(ostream &os)
	{
		os<<"翼缘参数1:\n";
		for(size_t i=0;i<plate1.size();++i)
			os<<plate1[i].x<<","<<plate1[i].thickness<<endl;
		os<<"\n翼缘参数2:\n";
		for(size_t i=0;i<plate2.size();++i)
			os<<plate2[i].x<<","<<plate2[i].thickness<<endl;
		os<<"\n支点翼缘参数3:\n";
		for(size_t i=0;i<plate3_zd.size();++i)
			os<<plate3_zd[i].x<<","<<plate3_zd[i].thickness<<endl;
		os<<"\n顶板腹板参数:\n";
		for(size_t i=0;i<web_top.size();++i)
			os<<web_top[i].x<<","<<web_top[i].thickness<<endl;
		os<<"\n支点底板腹板参数:\n";
		for(size_t i=0;i<web_bot_zd.size();++i)
			os<<web_bot_zd[i].x<<","<<web_bot_zd[i].thickness<<endl;
		os<<"\n截面变化参数:\nΔh(底板厚度变化值) H支(支点梁高)  H中(跨中梁高)\n";
		os<<H_delta<<ends<<Hz<<ends<<Hc<<endl;
		os<<"\n梁高曲线(0为跨中截面,1为支点截面,起点坐标为距墩的距离):\n";
		curve.output(os);
	}
	void outputprop(ostream &os)
	{
		os<<"plate1:\n";
		propety1.output(os);
		os<<"\nplate2:\n";
		propety2.output(os);
		os<<"\nplate3:\n";
		propety3.output(os);
		os<<"\nweb_top:\n";
		propety_web_top.output(os);
		os<<"\nweb_bot:\n";
		propety_web_bot.output(os);
	}
};

void VarSecBeam::CalcProp(const vector<CPlateGeom>& plate,CPlateProp& propety)
{
	propety.A=0;
	propety.C=0;
	propety.I=0;
	double tmp=0;
	for(size_t i=0;i<plate.size()-1;++i)
	{
		tmp=(plate[i+1].x-plate[i].x)*(plate[i+1].thickness+plate[i].thickness)/2.0;
		propety.C+=tmp*(pow(plate[i].thickness,2)+pow(plate[i+1].thickness,2)+plate[i+1].thickness*plate[i].thickness)
			/(plate[i].thickness+plate[i+1].thickness)/3.0;
		propety.A+=tmp;
		propety.I+=(plate[i+1].x-plate[i].x)/12*(pow(plate[i+1].thickness,2)+pow(plate[i].thickness,2))*(plate[i].thickness+plate[i+1].thickness);
	}
	propety.C/=propety.A;
	propety.I-=propety.A*pow(propety.C,2);
}

void VarSecBeam::PreCalculate()
{
	//计算propety1
	CalcProp(plate1,propety1);
	CalcProp(plate2,propety2);
	CalcProp(web_top,propety_web_top);

	plate3.resize(plate3_zd.size());
	web_bot.resize(web_bot_zd.size());
	for(size_t i=0;i<plate3_zd.size();++i)
	{
		plate3[i].x=plate3_zd[i].x;
	}
	for(size_t i=0;i<web_bot_zd.size();++i)
	{
		web_bot[i].x=web_bot_zd[i].x;
	}
	H_top_max=0;
	for(size_t i=1;i<web_top.size();++i)
	{
		H_top_max=max(H_top_max,web_top[i].thickness);
	}
}

bool VarSecBeam::Calculate(const double x)
{
	if(x<0||x>curve.keypoint.back().first) return false;
	double h_delta=curve(-H_delta,0,x);
	sec_propty.H=curve(Hc,Hz,x);

	for(size_t i=0;i<plate3.size();++i)
	{
		plate3[i].thickness=plate3_zd[i].thickness+h_delta;
	}
	for(size_t i=0;i<web_bot.size();++i)
	{
		web_bot[i].thickness=web_bot_zd[i].thickness+h_delta;
	}
	CalcProp(plate3,propety3);
	CalcProp(web_bot,propety_web_bot);
	sec_propty.A=propety1.A+propety2.A+propety3.A+propety_web_top.A+propety_web_bot.A;

	sec_propty.C=propety3.A*propety3.C+propety_web_bot.A*propety_web_bot.C;
	sec_propty.C+=propety1.A*(sec_propty.H-propety1.C)+propety2.A*(sec_propty.H-propety2.C);
	sec_propty.C+=propety_web_top.A*(sec_propty.H-propety_web_top.C);
	sec_propty.C/=sec_propty.A;
	sec_propty.A*=2;
	sec_propty.I1=(propety1.I+propety1.A*pow(sec_propty.H-sec_propty.C-propety1.C,2))*2;
	sec_propty.I2=(propety2.I+propety2.A*pow(sec_propty.H-sec_propty.C-propety2.C,2))*2;
	sec_propty.Iw=(propety_web_top.I+propety_web_top.A*pow(sec_propty.H-sec_propty.C-propety_web_top.C,2))*2;
	sec_propty.I3=(propety3.I+propety3.A*pow(sec_propty.C-propety3.C,2))*2;
	sec_propty.Iw+=(propety_web_bot.I+propety_web_bot.A*pow(sec_propty.C-propety_web_bot.C,2))*2;
	sec_propty.I=sec_propty.I1+sec_propty.I2+sec_propty.I3+sec_propty.Iw;
	return true;
}

bool VarSecBeam::Calculate(const double x,CSecPropty2 &sec_prop)
{
	Calculate(x);
	sec_prop.assign(sec_propty);
	return true;
}

double VarSecBeam::Calculate_dI(const double x,const double dx)
{
	Calculate(x-dx);
	double If=sec_propty.I1+sec_propty.I2+sec_propty.I3;
	Calculate(x+dx);
	double Ib=sec_propty.I1+sec_propty.I2+sec_propty.I3;
	return (Ib-If)/dx/2;
}


void VarSecBeam::Calculate_dI(const double x,CSecPropty2& SecProp,const double dx)
{
	Calculate(x-dx);
	double If1=sec_propty.I1;
	double If2=sec_propty.I2;
	double If3=sec_propty.I3;
	double Iw=sec_propty.Iw;
	Calculate(x+dx);
	SecProp.dI1=(sec_propty.I1-If1)/dx/2;
	SecProp.dI2=(sec_propty.I2-If2)/dx/2;
	SecProp.dI3=(sec_propty.I3-If3)/dx/2;
	SecProp.dIw=(sec_propty.Iw-Iw)/dx/2;
	SecProp.dI=SecProp.dI1+SecProp.dI2+SecProp.dI3+SecProp.dIw;
}

double VarSecBeam::CalculateHw(const double x)
{
	double H=curve(Hc,Hz,x);
	//往下开始写
	double h_delta=curve(-H_delta,0,x);
	double h_bot=web_bot_zd[0].thickness;
	for(size_t i=1;i<web_bot_zd.size();++i)
	{
		h_bot=max(h_bot,web_bot_zd[i].thickness);
	}
	h_bot+=h_delta;
	return H-h_bot-H_top_max;
};

struct Ccsw   //波形钢腹板
{
	Cmatrix<double> _data;   //每列依次为x坐标、厚度,单位m
inline double tw(double x);
};

double Ccsw::tw(double x)
{
	size_t index=1;
	if (_data.row()<=1) 
	{
		cout<<"波形钢腹板厚度数据为空!";
		throw;
	}
	if(x==_data(_data.row(),1)) 
		return _data(_data.row(),2);
	if(x>_data(_data.row(),1)||x<_data(1,1))
	{
		cout<<"纵坐标越界!";
		throw;
	}
	for(;index<=_data.row();++index)
	{
		if(x<_data(index,1)) {--index;break;}
	}
	return _data(index,2)+( _data(index+1,2)- _data(index,2))/( _data(index+1,1)- _data(index,1))*(x-_data(index,1));
};

struct Csection
{
	pair<size_t,double> beg,end;  //首单元号、起点位置距首单元i节点的距离比
	                                //末单元号，终点位置距末单元j节点的距离比
};

struct CMvldPos
{
	vector<Csection> pos_uniload_max,pos_uniload_min;  //列依次为首单元号、起点位置距首单元i节点的距离比,末单元号，终点位置距末单元j节点的距离比
	size_t pos_contload_max,pos_contload_min;                        //集中荷载作用的节点号
	double q;                                   //均布荷载，应为负值
	double F;                                   //集中荷载，应为负值
};

struct CMvldRst    //移动荷载加载位置
{
	vector<CUnifLoad>uniload_mx;
	vector<CUnifLoad>uniload_mn;
	pair<size_t,size_t>pos_ctload;
	inline void clear()
	{
		uniload_mx.clear();
		uniload_mn.clear();
		pos_ctload.first=pos_ctload.second=0;
	}
	inline void output(ostream &os)
	{
		os<<"正值区间:";
		for(size_t i=0;i<uniload_mx.size();++i)
		{
			os<<" ["<<uniload_mx[i].x1<<","<<uniload_mx[i].x2<<"]";
		}
		os<<"\n负值区间:";
		for(size_t i=0;i<uniload_mn.size();++i)
		{
			os<<" ["<<uniload_mn[i].x1<<","<<uniload_mn[i].x2<<"]";
		}
		os<<"\n最大、最小值位置:"<<pos_ctload.first<<","<<pos_ctload.second<<endl;
	}
};

struct CFluLine
{
	size_t NodeNum;
	double x;  //节点坐标
	Cmatrix<double> flu_line;   //坐标、值
	CMvldRst MvldRst;
	inline CFluLine():NodeNum(0),x(0){};
	
	inline void Calculate(const double q=-10.5,const double F=-360); 
	                                                           //影响线正区间的坐标，负区间的坐标，最大最小值的节点位置,pos_ctload编号为0表示最大或最小值为0
															   //返回值最大最小值,根据pos_ctload来判断影响线是否存在正值或负值
	inline void output(ostream& os);
	inline void output_rst(ostream&os);
};

inline double findzero(double y1,double y2)  //距起点的距离比
{
	return abs(y1)/(abs(y1)+abs(y2));
}
inline void CalcElemq(CFem &fem,const CUnifLoad &q);

void CFluLine::output(ostream &os)
{
	os<<"结点号码："<<NodeNum<<endl;
	os<<"x,M"<<endl;
	for(size_t i=1;i<=flu_line.row();++i)
		os<<flu_line(i,1)<<","<<flu_line(i,2)<<endl;
}

void CFluLine::output_rst(ostream&os)
{
	os<<"节点号:"<<NodeNum<<endl;
	MvldRst.output(os);
	os<<endl;
}

////已测试该函数准确
void CFluLine::Calculate(const double q,const double F)
{
	const double presion=1e-6;
	MvldRst.clear();
	size_t index=1;
	CUnifLoad unifload;
	unifload.q1=unifload.q2=q;
	for(;index<=flu_line.row();++index)
	{
		if(abs(flu_line(index,2))>presion) break;
	}
	if(index==flu_line.row()+1) return;
	if(index==flu_line.row())
	{
		unifload.x1=flu_line(index-1,1);
		unifload.x2=flu_line(index,1);
		if(flu_line(index,2)>0)
		{
			MvldRst.uniload_mx.push_back(unifload);
			MvldRst.pos_ctload.first=index;
			return;
		}
		else 
		{
			MvldRst.uniload_mn.push_back(unifload);
			MvldRst.pos_ctload.second=index;
			return;
		}
	}
	//unifload.x1=flu_line(index,1);
	index>=2?unifload.x1=flu_line(index-1,1):unifload.x1=flu_line(index,1);
	for(size_t i=index+1;i<=flu_line.row();++i)
	{
		if(abs(flu_line(i,2))<=presion) 
		{/*
			if(abs(flu_line(i-1,2))<=presion)
				pos_zero.back()=flu_line(i,1);
			else*/
			unifload.x2=flu_line(i,1);
			flu_line(i-1,2)>0?MvldRst.uniload_mx.push_back(unifload):MvldRst.uniload_mn.push_back(unifload);
			unifload.x1=unifload.x2;
		}
		else if(flu_line(i-1,2)*flu_line(i,2)<0&&abs(flu_line(i-1,2))>presion)
		{
			unifload.x2=findzero(flu_line(i-1,2),flu_line(i,2))*(flu_line(i,1)-flu_line(i-1,1))+flu_line(i-1,1);
			flu_line(i-1,2)>0?MvldRst.uniload_mx.push_back(unifload):MvldRst.uniload_mn.push_back(unifload);
			unifload.x1=unifload.x2;
		}
	}
	if(abs(flu_line(flu_line.row(),2))>presion)
	{
		unifload.x2=flu_line(flu_line.row(),1);
		flu_line(flu_line.row(),2)>0?MvldRst.uniload_mx.push_back(unifload):MvldRst.uniload_mn.push_back(unifload);
	}
	double max_value=flu_line(index,2),min_value=flu_line(index,2);
	MvldRst.pos_ctload.first=index;
	MvldRst.pos_ctload.second=index;
	for(size_t i=index+1;i<=flu_line.row();++i)
	{
		if(flu_line(i,2)>max_value)
		{
			max_value=flu_line(i,2);
			MvldRst.pos_ctload.first=i;
		}
		if(flu_line(i,2)<min_value)
		{
			min_value=flu_line(i,2);
			MvldRst.pos_ctload.second=i;
		}
	}
	if(max_value<0) MvldRst.pos_ctload.first=0;
	if(min_value>0) MvldRst.pos_ctload.second=0;
}

#endif