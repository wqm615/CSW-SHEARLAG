extern vector<size_t> SecDataCountBak;
extern vector<double> deltasBak;
extern vector<size_t> DiagIndexBak;
extern vector<bool> ConsolidNode;
inline size_t IndexAndSecToPos(const size_t index,const size_t sec);
inline size_t IndexAndSecToNodeNum(const size_t index,const size_t sec);
inline void NodeNumToIndexAndSec(const size_t NodeNum,vector< pair<size_t,size_t> >&IndexAndSec);
void CalcShearLagAdv3U(Cmatrix<double>& _rst,const CFem& fem,const CeqSetD &eqsd,const CeqSetD&eqsd_rev,const size_t index,const size_t secnum,const size_t NodeNum,const double MaxM);

inline bool InputUnifSec2(ifstream & ifs)  //变高度梁输入数据,不包括集中荷载和均布荷载
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

	//读入波形钢腹板厚度数据
	ifs_ignore(ifs,3);
	ifs>>tmp;
	csw._data.resize(tmp,2);
	for(size_t i=1;i<=tmp;++i)
	{
		ifs>>csw._data(i,1)>>csw._data(i,2);
	}
	ifs.close();

	//计算翼缘板平均厚度
	for(size_t i=0;i<var_sec.plate1.size()-1;++i)
	{
		La1+=(var_sec.plate1[i+1].x-var_sec.plate1[i].x)/2*(var_sec.plate1[i+1].thickness+var_sec.plate1[i].thickness);
	}
	La1/=var_sec.plate1.back().x;
	for(size_t i=0;i<var_sec.plate2.size()-1;++i)
	{
		La2+=(var_sec.plate2[i+1].x-var_sec.plate2[i].x)/2*(var_sec.plate2[i+1].thickness+var_sec.plate2[i].thickness);
	}
	La2/=var_sec.plate2.back().x;

	return true;
}

inline void GenFemMod3UShear2(CFem &fem,const size_t ndiv)  //此函数用于影响线的计算,无需荷载数据
{
	fem.AddMat(CMaterial(1,"Concrete",E,E/2.0/(1+0.2)));
	set<Ckeypos> keypos;  //均布荷载、集中荷载、边界
	pair<set<Ckeypos>::iterator,bool > status;
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

	vector<size_t> buffer;
	for(size_t i=0;i<Diaphragm.size();++i)
	{
		status=keypos.insert(Ckeypos(Diaphragm[i],0));
		buffer.push_back(distance(keypos.begin(),status.first));
	}

	DiagIndex.resize(keypos.size(),0);
	for(size_t i=0;i<buffer.size();++i)
		DiagIndex[buffer[i]]=1;
	//size_t tmp=0;
	//for(set<Ckeypos>::iterator it_tmp=keypos.begin();it_tmp!=keypos.end();++it_tmp)
	//{
	//	cout<<it_tmp->x<<ends<<DiagIndex[tmp]<<endl;
	//	++tmp;
	//}
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
		var_sec.Calculate_dI(buffer.x+min(1e-6,delta_min),buffer,min(1e-6,delta_min)*0.9);
		xIdI.push_back(buffer);
		for(size_t i=1;i<ndiv2;++i)
		{
			tmp=it1->x+i*deltas.back();
			fem.AddNode(++num,tmp,0);
			buffer.x=tmp;
			var_sec.Calculate(buffer.x,buffer);
			plate3.push_back(var_sec.plate3);
			var_sec.Calculate_dI(buffer.x,buffer,delta_min);
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
				var_sec.Calculate_dI(buffer.x-delta_min,buffer,delta_min*0.9);
				xIdI.push_back(buffer);
				var_sec.Calculate_dI(buffer.x+delta_min,buffer,delta_min*0.9);
				xIdI.push_back(buffer);
			}
			else 
			{
				var_sec.Calculate_dI(buffer.x,buffer,delta_min);
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
				var_sec.Calculate_dI(buffer.x,buffer,delta_min);
				xIdI.push_back(buffer);
			}
		}
		fem.AddNode(++num,L,0);
		buffer.x=L;
		var_sec.Calculate(buffer.x-min(1e-6,delta_min),buffer);
		plate3.push_back(var_sec.plate3);
		var_sec.Calculate_dI(buffer.x-min(1e-6,delta_min),buffer,min(1e-6,delta_min)*0.9);
		xIdI.push_back(buffer);
	}
	//单元
	{
		pair<double,double> AandI;
		std::stringstream ostr;
		size_t i=0;
		CFem::iterator_node it2=fem.NodeHead();
		double hw,tw;
		double fai,el,As;   //剪切变形影响系数,单元长度，波希钢腹板面积
		double G=7.9e7;
		for(CFem::iterator_node it=fem.NodeHead();it!=--fem.NodeEnd();++it)
		{
			ostr.str("");
			ostr<<"实常数"<<++i;
			++it2;
			var_sec.Calculate((it->x)/2+(it2->x)/2);
			hw=var_sec.CalculateHw((it->x)/2+(it2->x)/2);
			tw=csw.tw((it->x)/2+(it2->x)/2)/1000;
			As=hw*tw;
			el=it2->x-it->x;
			fai=12*E*var_sec.sec_propty.I/G/As/pow(el,2);
			//cout<<i<<ends<<el<<ends<<E/G<<ends<<var_sec.sec_propty.I<<ends<<fai<<endl;
			fem.AddReal(CReal3d(i,ostr.str(),var_sec.sec_propty.A,var_sec.sec_propty.I,var_sec.sec_propty.I,
				var_sec.sec_propty.I,Cmatrix<long double>(1,2),0,fai));
			fem.AddElem(CBeam3d(i,0,"空间梁单元"),i,i+1,1,i);
		}
	}
	//边界
	ConsolidNode.resize(fem.NodeCount(),false);
	//添加纵向约束
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
			if(MidPiers[i]==2) 
			{
				fem.AddBond(BdNodeNum,0,fem_aufunc::MY);
				ConsolidNode[BdNodeNum-1]=true;
			}
		}
	}
	SecDataCountBak=SecDataCount;
	deltasBak=deltas;
	DiagIndexBak=DiagIndex;
	ofstream ofs("femMvLd.dat");
	fem.OutputData(ofs);
	ofs.close();
}

inline void GenEqSet3U2(CeqSetD & eqsd,const CFem &fem)
{
	eqsd.initilize();
	const double G=E/2.0/(1+0.2);
	//(1)方程组
	Ceq eq;
	eq.AddItem(Citem(1,2));
	eq.AddItem(Citem(2,2));
	eq.AddItem(Citem(3,2));
	eq.AddItem(Citem(1,1));
	eq.AddItem(Citem(2,1));
	eq.AddItem(Citem(3,1));
	eq.AddItem(Citem(1,0));
	eqsd.AddEq(eq);
	eq._equ.back().var_index=2;
	eqsd.AddEq(eq);
	eq._equ.back().var_index=3;
	eqsd.AddEq(eq);

	//(2)节段数据
    eqsd.SetSecNum(SecDataCount);
	//(3)步长
	eqsd.SetDelta(deltas);
	//(4)mat的生成
	Cmatrix<double> mat;
	CFem::const_iterator_rst it_rst1=fem.EForceHead(),it_rst2=fem.EForceHead(),it_rst3=fem.EForceHead();
	/*ofstream ofs_I("I2.csv");
	ofs_I<<"I1,I2,I3,Iw,I,dI1,dI2,dI3,dIw,I\n";
	for(size_t i=0;i<xIdI.size();++i)
	{
		ofs_I<<xIdI[i].I1<<","<<xIdI[i].I2<<","<<xIdI[i].I3<<","<<xIdI[i].Iw<<","<<xIdI[i].I<<",";
		ofs_I<<xIdI[i].dI1<<","<<xIdI[i].dI2<<","<<xIdI[i].dI3<<","<<xIdI[i].dIw<<","<<xIdI[i].dI<<endl;
	}
	ofs_I.close();*/
	eqsd.ReSizeCoef();
	{
		size_t index1=0,index2=0,index3=0;
		const double b1=var_sec.plate1.back().x;
		const double b2=var_sec.plate2.back().x;
		const double b3=var_sec.plate3_zd.back().x;
		double dIs;
		size_t i=1;
		for(;i<=SecDataCount.size();++i)
		{
			//方程1
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index1].dI-xIdI[index1].dIw;
			mat(1,1)=9.0/16/xIdI[index1].I*pow(xIdI[index1].I1,2)-9.0/14*xIdI[index1].I1;
			mat(1,2)=9.0/16*xIdI[index1].I1*xIdI[index1].I2/xIdI[index1].I;
			mat(1,3)=9.0/16*xIdI[index1].I1*xIdI[index1].I3/xIdI[index1].I;
			mat(1,4)=9.0/8*xIdI[index1].I1*xIdI[index1].dI1/xIdI[index1].I-9.0/14*xIdI[index1].dI1;
			mat(1,4)-=9.0/16*pow(xIdI[index1].I1,2)*dIs/pow(xIdI[index1].I,2);
			mat(1,5)=xIdI[index1].I1*xIdI[index1].dI2/xIdI[index1].I+xIdI[index1].I2*xIdI[index1].dI1/xIdI[index1].I;
			mat(1,5)-=xIdI[index1].I2*xIdI[index1].I1*dIs/pow(xIdI[index1].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=xIdI[index1].I1*xIdI[index1].dI3/xIdI[index1].I+xIdI[index1].I3*xIdI[index1].dI1/xIdI[index1].I;
			mat(1,6)-=xIdI[index1].I1*xIdI[index1].I3*dIs/pow(xIdI[index1].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index1].I1/E/pow(b1,2);
			mat(1,8)=-3.0/4*xIdI[index1].I1/E/xIdI[index1].I*(*it_rst1)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index1].dI1-xIdI[index1].I1*xIdI[index1].dI/xIdI[index1].I)*(*it_rst1)(1,5)/E/xIdI[index1].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index1;
				dIs=xIdI[index1].dI-xIdI[index1].dIw;
				mat(j,1)=9.0/16/xIdI[index1].I*pow(xIdI[index1].I1,2)-9.0/14*xIdI[index1].I1;
				mat(j,2)=9.0/16*xIdI[index1].I1*xIdI[index1].I2/xIdI[index1].I;
				mat(j,3)=9.0/16*xIdI[index1].I1*xIdI[index1].I3/xIdI[index1].I;
				mat(j,4)=9.0/8*xIdI[index1].I1*xIdI[index1].dI1/xIdI[index1].I-9.0/14*xIdI[index1].dI1;
				mat(j,4)-=9.0/16*pow(xIdI[index1].I1,2)*dIs/pow(xIdI[index1].I,2);
				mat(j,5)=xIdI[index1].I1*xIdI[index1].dI2/xIdI[index1].I+xIdI[index1].I2*xIdI[index1].dI1/xIdI[index1].I;
				mat(j,5)-=xIdI[index1].I2*xIdI[index1].I1*dIs/pow(xIdI[index1].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=xIdI[index1].I1*xIdI[index1].dI3/xIdI[index1].I+xIdI[index1].I3*xIdI[index1].dI1/xIdI[index1].I;
				mat(j,6)-=xIdI[index1].I1*xIdI[index1].I3*dIs/pow(xIdI[index1].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index1].I1/E/pow(b1,2);
				mat(j,8)=-3.0/4*xIdI[index1].I1/E/xIdI[index1].I*(*it_rst1)(2,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index1].dI1-xIdI[index1].I1*xIdI[index1].dI/xIdI[index1].I)*(*it_rst1)(2,5)/E/xIdI[index1].I;
				++it_rst1;
			}
			eqsd.SetCoef(mat,1,i);
			//方程2
			dIs=xIdI[index2].dI-xIdI[index2].dIw;
			mat(1,1)=9.0/16*xIdI[index2].I1*xIdI[index2].I2/xIdI[index2].I;
			mat(1,2)=9.0/16*pow(xIdI[index2].I2,2)/xIdI[index2].I-9.0/14*xIdI[index2].I2;
			mat(1,3)=9.0/16*xIdI[index2].I2*xIdI[index2].I3/xIdI[index2].I;
			mat(1,4)=xIdI[index2].I2*xIdI[index2].dI1/xIdI[index2].I+xIdI[index2].I1*xIdI[index2].dI2/xIdI[index2].I;
			mat(1,4)-=xIdI[index2].I1*xIdI[index2].I2*dIs/pow(xIdI[index2].I,2);
			mat(1,4)*=9.0/16;
			mat(1,5)=9.0/8*xIdI[index2].I2*xIdI[index2].dI2/xIdI[index2].I-9.0/14*xIdI[index2].dI2;
			mat(1,5)-=9.0/16*pow(xIdI[index2].I2,2)*dIs/pow(xIdI[index2].I,2);
			mat(1,6)=xIdI[index2].I2*xIdI[index2].dI3/xIdI[index2].I+xIdI[index2].I3*xIdI[index2].dI2/xIdI[index2].I;
			mat(1,6)-=xIdI[index2].I2*xIdI[index2].I3*dIs/pow(xIdI[index2].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index2].I2/E/pow(b2,2);
			mat(1,8)=-3.0/4*xIdI[index2].I2/E/xIdI[index2].I*(*it_rst2)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index2].dI2-xIdI[index2].I2*xIdI[index2].dI/xIdI[index2].I)*(*it_rst2)(1,5)/E/xIdI[index2].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index2;
				dIs=xIdI[index2].dI-xIdI[index2].dIw;
				mat(j,1)=9.0/16*xIdI[index2].I1*xIdI[index2].I2/xIdI[index2].I;
				mat(j,2)=9.0/16*pow(xIdI[index2].I2,2)/xIdI[index2].I-9.0/14*xIdI[index2].I2;
				mat(j,3)=9.0/16*xIdI[index2].I2*xIdI[index2].I3/xIdI[index2].I;
				mat(j,4)=xIdI[index2].I2*xIdI[index2].dI1/xIdI[index2].I+xIdI[index2].I1*xIdI[index2].dI2/xIdI[index2].I;
				mat(j,4)-=xIdI[index2].I1*xIdI[index2].I2*dIs/pow(xIdI[index2].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=9.0/8*xIdI[index2].I2*xIdI[index2].dI2/xIdI[index2].I-9.0/14*xIdI[index2].dI2;
				mat(j,5)-=9.0/16*pow(xIdI[index2].I2,2)*dIs/pow(xIdI[index2].I,2);
				mat(j,6)=xIdI[index2].I2*xIdI[index2].dI3/xIdI[index2].I+xIdI[index2].I3*xIdI[index2].dI2/xIdI[index2].I;
				mat(j,6)-=xIdI[index2].I2*xIdI[index2].I3*dIs/pow(xIdI[index2].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index2].I2/E/pow(b2,2);
				mat(j,8)=-3.0/4*xIdI[index2].I2/E/xIdI[index2].I*(*it_rst2)(2,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index2].dI2-xIdI[index2].I2*xIdI[index2].dI/xIdI[index2].I)*(*it_rst2)(2,5)/E/xIdI[index2].I;
				++it_rst2;
			}
			eqsd.SetCoef(mat,2,i);
			//方程3
			dIs=xIdI[index3].dI-xIdI[index3].dIw;
			mat(1,1)=9.0/16*xIdI[index3].I1*xIdI[index3].I3/xIdI[index3].I;
			mat(1,2)=9.0/16*xIdI[index3].I2*xIdI[index3].I3/xIdI[index3].I;
			mat(1,3)=9.0/16*pow(xIdI[index3].I3,2)/xIdI[index3].I-9.0/14*xIdI[index3].I3;	
			mat(1,4)=xIdI[index3].I3*xIdI[index3].dI1/xIdI[index3].I+xIdI[index3].I1*xIdI[index3].dI3/xIdI[index3].I;
			mat(1,4)-=xIdI[index3].I1*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
			mat(1,4)*=9.0/16;
			mat(1,5)=xIdI[index3].I2*xIdI[index3].dI3/xIdI[index3].I+xIdI[index3].I3*xIdI[index3].dI2/xIdI[index3].I;
			mat(1,5)-=xIdI[index3].I2*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=9.0/8*xIdI[index3].I3*xIdI[index3].dI3/xIdI[index3].I-9.0/14*xIdI[index3].dI3;
			mat(1,6)-=9.0/16*pow(xIdI[index3].I3,2)*dIs/pow(xIdI[index3].I,2);
			mat(1,7)=9.0/5*G*xIdI[index3].I3/E/pow(b3,2);
			mat(1,8)=-3.0/4*xIdI[index3].I3/E/xIdI[index3].I*(*it_rst3)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index3].dI3-xIdI[index3].I3*xIdI[index3].dI/xIdI[index3].I)*(*it_rst3)(1,5)/E/xIdI[index3].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index3;
				dIs=xIdI[index3].dI-xIdI[index3].dIw;
				mat(j,1)=9.0/16*xIdI[index3].I1*xIdI[index3].I3/xIdI[index3].I;
				mat(j,2)=9.0/16*xIdI[index3].I2*xIdI[index3].I3/xIdI[index3].I;
				mat(j,3)=9.0/16*pow(xIdI[index3].I3,2)/xIdI[index3].I-9.0/14*xIdI[index3].I3;	
				mat(j,4)=xIdI[index3].I3*xIdI[index3].dI1/xIdI[index3].I+xIdI[index3].I1*xIdI[index3].dI3/xIdI[index3].I;
				mat(j,4)-=xIdI[index3].I1*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=xIdI[index3].I2*xIdI[index3].dI3/xIdI[index3].I+xIdI[index3].I3*xIdI[index3].dI2/xIdI[index3].I;
				mat(j,5)-=xIdI[index3].I2*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=9.0/8*xIdI[index3].I3*xIdI[index3].dI3/xIdI[index3].I-9.0/14*xIdI[index3].dI3;
				mat(j,6)-=9.0/16*pow(xIdI[index3].I3,2)*dIs/pow(xIdI[index3].I,2);
				mat(j,7)=9.0/5*G*xIdI[index3].I3/E/pow(b3,2);
				mat(j,8)=-3.0/4*xIdI[index3].I3/E/xIdI[index3].I*(*it_rst3)(2,3);  //Q=(*it_rst)(2,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index3].dI3-xIdI[index3].I3*xIdI[index3].dI/xIdI[index3].I)*(*it_rst3)(2,5)/E/xIdI[index3].I;
				++it_rst3;
			}
			eqsd.SetCoef(mat,3,i);
			if(SecDataCount[i-1]==SecDataCountBak[i-1])
			{
				++index1;
				++index2;
				++index3;
			}
			else 
			{
				++i;
				break;
			}
		}
		for(;i<=SecDataCount.size();++i)
		{
			//方程1
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index1].dI-xIdI[index1].dIw;
			mat(1,1)=9.0/16/xIdI[index1].I*pow(xIdI[index1].I1,2)-9.0/14*xIdI[index1].I1;
			mat(1,2)=9.0/16*xIdI[index1].I1*xIdI[index1].I2/xIdI[index1].I;
			mat(1,3)=9.0/16*xIdI[index1].I1*xIdI[index1].I3/xIdI[index1].I;
			mat(1,4)=9.0/8*xIdI[index1].I1*xIdI[index1].dI1/xIdI[index1].I-9.0/14*xIdI[index1].dI1;
			mat(1,4)-=9.0/16*pow(xIdI[index1].I1,2)*dIs/pow(xIdI[index1].I,2);
			mat(1,5)=xIdI[index1].I1*xIdI[index1].dI2/xIdI[index1].I+xIdI[index1].I2*xIdI[index1].dI1/xIdI[index1].I;
			mat(1,5)-=xIdI[index1].I2*xIdI[index1].I1*dIs/pow(xIdI[index1].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=xIdI[index1].I1*xIdI[index1].dI3/xIdI[index1].I+xIdI[index1].I3*xIdI[index1].dI1/xIdI[index1].I;
			mat(1,6)-=xIdI[index1].I1*xIdI[index1].I3*dIs/pow(xIdI[index1].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index1].I1/E/pow(b1,2);
			mat(1,8)=-3.0/4*xIdI[index1].I1/E/xIdI[index1].I*(*it_rst1)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index1].dI1-xIdI[index1].I1*xIdI[index1].dI/xIdI[index1].I)*(*it_rst1)(1,5)/E/xIdI[index1].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index1;
				dIs=xIdI[index1].dI-xIdI[index1].dIw;
				mat(j,1)=9.0/16/xIdI[index1].I*pow(xIdI[index1].I1,2)-9.0/14*xIdI[index1].I1;
				mat(j,2)=9.0/16*xIdI[index1].I1*xIdI[index1].I2/xIdI[index1].I;
				mat(j,3)=9.0/16*xIdI[index1].I1*xIdI[index1].I3/xIdI[index1].I;
				mat(j,4)=9.0/8*xIdI[index1].I1*xIdI[index1].dI1/xIdI[index1].I-9.0/14*xIdI[index1].dI1;
				mat(j,4)-=9.0/16*pow(xIdI[index1].I1,2)*dIs/pow(xIdI[index1].I,2);
				mat(j,5)=xIdI[index1].I1*xIdI[index1].dI2/xIdI[index1].I+xIdI[index1].I2*xIdI[index1].dI1/xIdI[index1].I;
				mat(j,5)-=xIdI[index1].I2*xIdI[index1].I1*dIs/pow(xIdI[index1].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=xIdI[index1].I1*xIdI[index1].dI3/xIdI[index1].I+xIdI[index1].I3*xIdI[index1].dI1/xIdI[index1].I;
				mat(j,6)-=xIdI[index1].I1*xIdI[index1].I3*dIs/pow(xIdI[index1].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index1].I1/E/pow(b1,2);
				mat(j,8)=-3.0/4*xIdI[index1].I1/E/xIdI[index1].I*(*it_rst1)(2,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index1].dI1-xIdI[index1].I1*xIdI[index1].dI/xIdI[index1].I)*(*it_rst1)(2,5)/E/xIdI[index1].I;
				++it_rst1;
			}
			eqsd.SetCoef(mat,1,i);
			//方程2
			dIs=xIdI[index2].dI-xIdI[index2].dIw;
			mat(1,1)=9.0/16*xIdI[index2].I1*xIdI[index2].I2/xIdI[index2].I;
			mat(1,2)=9.0/16*pow(xIdI[index2].I2,2)/xIdI[index2].I-9.0/14*xIdI[index2].I2;
			mat(1,3)=9.0/16*xIdI[index2].I2*xIdI[index2].I3/xIdI[index2].I;
			mat(1,4)=xIdI[index2].I2*xIdI[index2].dI1/xIdI[index2].I+xIdI[index2].I1*xIdI[index2].dI2/xIdI[index2].I;
			mat(1,4)-=xIdI[index2].I1*xIdI[index2].I2*dIs/pow(xIdI[index2].I,2);
			mat(1,4)*=9.0/16;
			mat(1,5)=9.0/8*xIdI[index2].I2*xIdI[index2].dI2/xIdI[index2].I-9.0/14*xIdI[index2].dI2;
			mat(1,5)-=9.0/16*pow(xIdI[index2].I2,2)*dIs/pow(xIdI[index2].I,2);
			mat(1,6)=xIdI[index2].I2*xIdI[index2].dI3/xIdI[index2].I+xIdI[index2].I3*xIdI[index2].dI2/xIdI[index2].I;
			mat(1,6)-=xIdI[index2].I2*xIdI[index2].I3*dIs/pow(xIdI[index2].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index2].I2/E/pow(b2,2);
			mat(1,8)=-3.0/4*xIdI[index2].I2/E/xIdI[index2].I*(*it_rst2)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index2].dI2-xIdI[index2].I2*xIdI[index2].dI/xIdI[index2].I)*(*it_rst2)(1,5)/E/xIdI[index2].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index2;
				dIs=xIdI[index2].dI-xIdI[index2].dIw;
				mat(j,1)=9.0/16*xIdI[index2].I1*xIdI[index2].I2/xIdI[index2].I;
				mat(j,2)=9.0/16*pow(xIdI[index2].I2,2)/xIdI[index2].I-9.0/14*xIdI[index2].I2;
				mat(j,3)=9.0/16*xIdI[index2].I2*xIdI[index2].I3/xIdI[index2].I;
				mat(j,4)=xIdI[index2].I2*xIdI[index2].dI1/xIdI[index2].I+xIdI[index2].I1*xIdI[index2].dI2/xIdI[index2].I;
				mat(j,4)-=xIdI[index2].I1*xIdI[index2].I2*dIs/pow(xIdI[index2].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=9.0/8*xIdI[index2].I2*xIdI[index2].dI2/xIdI[index2].I-9.0/14*xIdI[index2].dI2;
				mat(j,5)-=9.0/16*pow(xIdI[index2].I2,2)*dIs/pow(xIdI[index2].I,2);
				mat(j,6)=xIdI[index2].I2*xIdI[index2].dI3/xIdI[index2].I+xIdI[index2].I3*xIdI[index2].dI2/xIdI[index2].I;
				mat(j,6)-=xIdI[index2].I2*xIdI[index2].I3*dIs/pow(xIdI[index2].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index2].I2/E/pow(b2,2);
				mat(j,8)=-3.0/4*xIdI[index2].I2/E/xIdI[index2].I*(*it_rst2)(2,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index2].dI2-xIdI[index2].I2*xIdI[index2].dI/xIdI[index2].I)*(*it_rst2)(2,5)/E/xIdI[index2].I;
				++it_rst2;
			}
			eqsd.SetCoef(mat,2,i);
			//方程3
			dIs=xIdI[index3].dI-xIdI[index3].dIw;
			mat(1,1)=9.0/16*xIdI[index3].I1*xIdI[index3].I3/xIdI[index3].I;
			mat(1,2)=9.0/16*xIdI[index3].I2*xIdI[index3].I3/xIdI[index3].I;
			mat(1,3)=9.0/16*pow(xIdI[index3].I3,2)/xIdI[index3].I-9.0/14*xIdI[index3].I3;	
			mat(1,4)=xIdI[index3].I3*xIdI[index3].dI1/xIdI[index3].I+xIdI[index3].I1*xIdI[index3].dI3/xIdI[index3].I;
			mat(1,4)-=xIdI[index3].I1*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
			mat(1,4)*=9.0/16;
			mat(1,5)=xIdI[index3].I2*xIdI[index3].dI3/xIdI[index3].I+xIdI[index3].I3*xIdI[index3].dI2/xIdI[index3].I;
			mat(1,5)-=xIdI[index3].I2*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=9.0/8*xIdI[index3].I3*xIdI[index3].dI3/xIdI[index3].I-9.0/14*xIdI[index3].dI3;
			mat(1,6)-=9.0/16*pow(xIdI[index3].I3,2)*dIs/pow(xIdI[index3].I,2);
			mat(1,7)=9.0/5*G*xIdI[index3].I3/E/pow(b3,2);
			mat(1,8)=-3.0/4*xIdI[index3].I3/E/xIdI[index3].I*(*it_rst3)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index3].dI3-xIdI[index3].I3*xIdI[index3].dI/xIdI[index3].I)*(*it_rst3)(1,5)/E/xIdI[index3].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index3;
				dIs=xIdI[index3].dI-xIdI[index3].dIw;
				mat(j,1)=9.0/16*xIdI[index3].I1*xIdI[index3].I3/xIdI[index3].I;
				mat(j,2)=9.0/16*xIdI[index3].I2*xIdI[index3].I3/xIdI[index3].I;
				mat(j,3)=9.0/16*pow(xIdI[index3].I3,2)/xIdI[index3].I-9.0/14*xIdI[index3].I3;	
				mat(j,4)=xIdI[index3].I3*xIdI[index3].dI1/xIdI[index3].I+xIdI[index3].I1*xIdI[index3].dI3/xIdI[index3].I;
				mat(j,4)-=xIdI[index3].I1*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=xIdI[index3].I2*xIdI[index3].dI3/xIdI[index3].I+xIdI[index3].I3*xIdI[index3].dI2/xIdI[index3].I;
				mat(j,5)-=xIdI[index3].I2*xIdI[index3].I3*dIs/pow(xIdI[index3].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=9.0/8*xIdI[index3].I3*xIdI[index3].dI3/xIdI[index3].I-9.0/14*xIdI[index3].dI3;
				mat(j,6)-=9.0/16*pow(xIdI[index3].I3,2)*dIs/pow(xIdI[index3].I,2);
				mat(j,7)=9.0/5*G*xIdI[index3].I3/E/pow(b3,2);
				mat(j,8)=-3.0/4*xIdI[index3].I3/E/xIdI[index3].I*(*it_rst3)(2,3);  //Q=(*it_rst)(2,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index3].dI3-xIdI[index3].I3*xIdI[index3].dI/xIdI[index3].I)*(*it_rst3)(2,5)/E/xIdI[index3].I;
				++it_rst3;
			}
			eqsd.SetCoef(mat,3,i);
			++index1;
			++index2;
			++index3;
		}
	}
	//边界条件
	CeqBD eqbd;
	{
		size_t var_index=1;
		eqbd.AddItem(CitemBD(1,var_index,1,1,1));//分段号 变量号 阶数 离散点号 系数
		eqbd.SetValue(0);
		eqsd.AddEqb(eqbd);
		for(size_t i=1;i<SecDataCount.size();++i)
		{
			eqbd.initilize();
			eqbd.AddItem(CitemBD(i,var_index,1,SecDataCount[i-1],1));
			eqbd.AddItem(CitemBD(i+1,var_index,1,1,-1));
			eqbd.SetValue(0);
			eqsd.AddEqb(eqbd);

			eqbd.initilize();
			eqbd.AddItem(CitemBD(i,var_index,0,SecDataCount[i-1],1));
			eqbd.AddItem(CitemBD(i+1,var_index,0,1,-1));
			eqbd.SetValue(0);
			eqsd.AddEqb(eqbd);
		}
		eqbd.initilize();
		eqbd.AddItem(CitemBD(SecDataCount.size(),var_index,1,SecDataCount.back(),1));
		eqbd.SetValue(0);
		eqsd.AddEqb(eqbd);
		eqbd.initilize();
	}
	for(size_t var_index=2;var_index<=3;++var_index)
	{
		if(DiagIndex[0]==0)
		{
			eqbd.AddItem(CitemBD(1,var_index,1,1,1));//分段号 变量号 阶数 离散点号 系数
			eqbd.SetValue(0);
			eqsd.AddEqb(eqbd);
		}
		else
		{
			eqbd.AddItem(CitemBD(1,var_index,0,1,1));//分段号 变量号 阶数 离散点号 系数
			eqbd.SetValue(0);
			eqsd.AddEqb(eqbd);
		}
		for(size_t i=1;i<SecDataCount.size();++i)
		{
			eqbd.initilize();
			if(DiagIndex[i]==0)
			{
				eqbd.AddItem(CitemBD(i,var_index,1,SecDataCount[i-1],1));
				eqbd.AddItem(CitemBD(i+1,var_index,1,1,-1));
				eqbd.SetValue(0);
				eqsd.AddEqb(eqbd);

				eqbd.initilize();
				eqbd.AddItem(CitemBD(i,var_index,0,SecDataCount[i-1],1));
				eqbd.AddItem(CitemBD(i+1,var_index,0,1,-1));
				eqbd.SetValue(0);
				eqsd.AddEqb(eqbd);
			}
			else
			{
				eqbd.AddItem(CitemBD(i,var_index,0,SecDataCount[i-1],1));
				eqbd.SetValue(0);
				eqsd.AddEqb(eqbd);

				eqbd.initilize();
				eqbd.AddItem(CitemBD(i+1,var_index,0,1,1));
				eqbd.SetValue(0);
				eqsd.AddEqb(eqbd);
			}
		}
		eqbd.initilize();
		if(DiagIndex.back()==0)
		{
			eqbd.AddItem(CitemBD(SecDataCount.size(),var_index,1,SecDataCount.back(),1));
			eqbd.SetValue(0);
			eqsd.AddEqb(eqbd);
		}
		else
		{
			eqbd.AddItem(CitemBD(SecDataCount.size(),var_index,0,SecDataCount.back(),1));
			eqbd.SetValue(0);
			eqsd.AddEqb(eqbd);
		}
		eqbd.initilize();
	}

	eqsd.Calculate();
}

//已测试该函数准确
inline void ModifySecDelDia(const size_t pos)  //根据集中荷载位置x,修改SecDataCount,deltas,DiagIndex
{
	if(pos==0||pos==1) return;
	SecDataCount=SecDataCountBak;
	deltas=deltasBak;
	DiagIndex=DiagIndexBak;
	size_t index=0;
	vector<size_t>::iterator it=SecDataCount.begin();
	vector<double>::iterator it2=deltas.begin();
	size_t pos2=1;
	for(;index<SecDataCountBak.size();++index)
	{
		pos2+=SecDataCountBak[index]-1;
		if(pos2==pos) return;
		else if(pos2>pos) break;
	}
	pos2=pos-(pos2-SecDataCountBak[index]+1);
	advance(it,index);
	advance(it2,index);
	it=SecDataCount.insert(it,pos2+1);
	++it;
	(*it)-=pos2;
	deltas.insert(it2,deltasBak[index]);
	it=DiagIndex.begin();
	advance(it,index+1);
	DiagIndex.insert(it,0);
}

//应先运行 GenFemMod3UShear2
inline void CalcMvLoad(Cmatrix<double>&rst_max,Cmatrix<double>&rst_min,CFem& _fem,CFluLine & _flu,const double MaxM,double q=-10.5,double F=-360)  //
{																		        
	_flu.Calculate(q,F);
	CeqSetD eqsd;
	vector< pair<size_t,size_t> > IndexAndSec;
	Cmatrix<double> _rst;
	rst_max.resize(1,5); //依次为M λ rou1 rou2 rou3
	rst_min.resize(1,5);
	CeqSetD eqsd_rev;
	time_t timep1,timep2;   //time_t为int型
	tm *tim_local;        //一个time结构的结构体
//最大值加载
	cout<<"   最大值加载:\n";
	if(_flu.MvldRst.pos_ctload.first)
	{
		_fem.ClearLoadAll();
		cout<<"\t01 有限元加载……";
		for(size_t i=0;i<_flu.MvldRst.uniload_mx.size();++i)
			CalcElemq(_fem,_flu.MvldRst.uniload_mx[i]);
		_fem.AddLoadn(_flu.MvldRst.pos_ctload.first,F);
		cout<<"OK!\n";
		cout<<"\t02 有限元计算……";
		_fem.Calculate(sparse_matrix());
		cout<<"OK!\n";
		cout<<"\t03 ModifySecDelDia……";
		ModifySecDelDia(_flu.MvldRst.pos_ctload.first);
		cout<<"OK!\n";
		cout<<"\t04 GenEqSet3U2……";
		time (&timep1);	//获得从1900年1月1日到现在所经过GM时间的秒数
		tim_local=localtime(&timep1);
		cout<<tim_local->tm_mon+1<<"-"<<tim_local->tm_mday<<ends<<tim_local->tm_hour<<":"<<tim_local->tm_min<<":"<<tim_local->tm_sec<<"……";
		GenEqSet3U2(eqsd,_fem);
		time(&timep2);
		cout<<"OK!  (Time consuming:"<<timep2-timep1<<"s)\n";
		cout<<"\t05 eqsd.Reverse……";
		eqsd.Reverse(eqsd_rev);
		cout<<"OK!\n";
		cout<<"\t06 eqsd_rev.Calculate……";
		time (&timep1);	//获得从1900年1月1日到现在所经过GM时间的秒数
		tim_local=localtime(&timep1);
		cout<<tim_local->tm_mon+1<<"-"<<tim_local->tm_mday<<ends<<tim_local->tm_hour<<":"<<tim_local->tm_min<<":"<<tim_local->tm_sec<<"……";
		eqsd_rev.Calculate();
		time(&timep2);
		cout<<"OK!  (Time consuming:"<<timep2-timep1<<"s)\n";
		cout<<"\t07 NodeNumToIndexAndSec……";
		NodeNumToIndexAndSec(_flu.NodeNum,IndexAndSec);
		cout<<"OK!\n";
		if(IndexAndSec.size()==2) rst_max.resize(2,5);
		for(size_t i=0;i<IndexAndSec.size();++i)
		{
			cout<<"\t08 CalcShearLagAdv3U……";
			CalcShearLagAdv3U(_rst,_fem,eqsd,eqsd_rev,IndexAndSec[i].first,IndexAndSec[i].second,_flu.NodeNum,MaxM);
			cout<<"OK!\n";
			rst_max(i+1,1)=_rst(1,7);
			rst_max(i+1,2)=_rst(1,8);
			rst_max(i+1,3)=_rst(1,9);
			rst_max(i+1,4)=_rst(1,10);
			rst_max(i+1,5)=_rst(1,11);
		}
	}
	else rst_max.clear();
//最小值加载
	cout<<"   最小值加载:\n";
	if(_flu.MvldRst.pos_ctload.second)
	{
		_fem.ClearLoadAll();
		cout<<"\t01 有限元加载……";
		for(size_t i=0;i<_flu.MvldRst.uniload_mn.size();++i)
			CalcElemq(_fem,_flu.MvldRst.uniload_mn[i]);
		_fem.AddLoadn(_flu.MvldRst.pos_ctload.second,F);
		cout<<"OK!\n";
		cout<<"\t02 有限元计算……";
		_fem.Calculate(sparse_matrix());
		cout<<"OK!\n";
		cout<<"\t03 ModifySecDelDia……";
		ModifySecDelDia(_flu.MvldRst.pos_ctload.second);
		cout<<"OK!\n";
		cout<<"\t04 GenEqSet3U2……";
		time (&timep1);	//获得从1900年1月1日到现在所经过GM时间的秒数
		tim_local=localtime(&timep1);
		cout<<tim_local->tm_mon+1<<"-"<<tim_local->tm_mday<<ends<<tim_local->tm_hour<<":"<<tim_local->tm_min<<":"<<tim_local->tm_sec<<"……";
		GenEqSet3U2(eqsd,_fem);
		time(&timep2);
		cout<<"OK!  (Time consuming:"<<timep2-timep1<<"s)\n";
		cout<<"\t05 eqsd.Reverse……";
		eqsd.Reverse(eqsd_rev);
		cout<<"OK!\n";
		cout<<"\t06 eqsd_rev.Calculate……";
		time (&timep1);	//获得从1900年1月1日到现在所经过GM时间的秒数
		tim_local=localtime(&timep1);
		cout<<tim_local->tm_mon+1<<"-"<<tim_local->tm_mday<<ends<<tim_local->tm_hour<<":"<<tim_local->tm_min<<":"<<tim_local->tm_sec<<"……";
		eqsd_rev.Calculate();
		time(&timep2);
		cout<<"OK!  (Time consuming:"<<timep2-timep1<<"s)\n";
		cout<<"\t07 NodeNumToIndexAndSec……";
		NodeNumToIndexAndSec(_flu.NodeNum,IndexAndSec);
		cout<<"OK!\n";
		if(IndexAndSec.size()==2) rst_min.resize(2,5);
		for(size_t i=0;i<IndexAndSec.size();++i)
		{
			cout<<"\t08 CalcShearLagAdv3U……";
			CalcShearLagAdv3U(_rst,_fem,eqsd,eqsd_rev,IndexAndSec[i].first,IndexAndSec[i].second,_flu.NodeNum,MaxM);
			cout<<"OK!\n";
			rst_min(i+1,1)=_rst(1,7);	
			rst_min(i+1,2)=_rst(1,8);
			rst_min(i+1,3)=_rst(1,9);
			rst_min(i+1,4)=_rst(1,10);
			rst_min(i+1,5)=_rst(1,11);
		}
	}
	else rst_min.clear();
}

//已通过测试
inline void NodeNumToIndexAndSec(const size_t NodeNum,vector< pair<size_t,size_t> >&IndexAndSec)  //根据节点号计算SecDataCount中的index和Sec
{
	IndexAndSec.resize(1);
	int count=int(NodeNum)-1;
	size_t i=0;
	for(;i<SecDataCount.size();++i)
	{
		count-=SecDataCount[i]-1;
		if(count<=0) break;
	}
	IndexAndSec[0].second=i+1;
	IndexAndSec[0].first=SecDataCount[i]+count;
	if(count==0&&i!=SecDataCount.size()-1)
	{
		IndexAndSec.resize(2);
		IndexAndSec[1].first=1;
		IndexAndSec[1].second=IndexAndSec[0].second+1;
	}
}

inline size_t IndexAndSecToNodeNum(const size_t index,const size_t sec)  //根据SecDataCount中的index和Sec,返回节点号
{
	size_t rtn=0;
	for(size_t i=1;i<sec;++i)
	{
		rtn+=SecDataCount[i-1]-1;
	}
	return rtn+index;
}

//已通过测试
inline size_t IndexAndSecToPos(const size_t index,const size_t sec)  //根据SecDataCount中的index和Sec,返回SecDataCountBak中的pos,该pos为plate3、xIdI的索引
{
	if(sec>SecDataCount.size())
	{
		cout<<"NodeNumToPos():sec越界\n";
		throw;
	}
	if(index>SecDataCount[sec-1])
	{
		cout<<"NodeNumToPos():index越界\n";
		throw;
	}
	size_t i=1,pos=0;
	for(;i<sec;++i)
	{
		if(SecDataCount[i-1]==SecDataCountBak[i-1])
			pos+=SecDataCount[i-1];
		else 
		{
			if(sec-i==1) 
			{
					pos+=SecDataCount[i-1]+index-1;
					return pos;
			}
			else
			{
				pos+=SecDataCountBak[i-1];
				for(size_t j=i+1;j<sec-1;++j)
				{
					pos+=SecDataCountBak[j-1];
				}
				pos+=index;
				return pos;
			}
		}
	}
	pos+=index;
	return pos;
}

//已通过测试
double CalcFluLine(CFem& _fem,vector<CFluLine> &_flu)   //返回最大值
{
	{
		size_t consoild_count=0;
		for(size_t i=0;i<ConsolidNode.size();++i)
			if(ConsolidNode[i]) ++consoild_count;
		_flu.resize(_fem.NodeCount()+consoild_count);
	}
	for(size_t i=0;i<_flu.size();++i) 
		_flu[i].flu_line.resize(_fem.NodeCount(),2);
	for(size_t i=0;i<_flu.size();++i)
	{
		_flu[i].flu_line.resize(_fem.NodeCount(),2);
	}
	CFem::const_iterator_node it_node=_fem.NodeHead();
	for(size_t i=1;i<=_fem.NodeCount();++i)
	{
		_fem.ClearLoadAll();
		_fem.AddLoadn(i,-1);
		_fem.Calculate(sparse_matrix());
		CFem::const_iterator_rst it_rst=_fem.EForceHead();
		//_flu[0].NodeNum=1;
		_flu[0].flu_line(i,1)=it_node->x;
		_flu[0].flu_line(i,2)=(*it_rst)(1,5);
		size_t index=1;
		for(size_t j=2;j<=_fem.NodeCount();++j)
		{
			//_flu[index].NodeNum=j;
			_flu[index].flu_line(i,1)=it_node->x;
			_flu[index].flu_line(i,2)=(*it_rst)(2,5);
			if(ConsolidNode[j-1])
			{
				++it_rst;
				++index;
				//_flu[index].NodeNum=j;
				_flu[index].flu_line(i,1)=it_node->x;
				_flu[index].flu_line(i,2)=(*it_rst)(1,5);
				--it_rst;
			}
			++it_rst;
			++index;
		}
		++it_node;
	}
	it_node=_fem.NodeHead();
	size_t index=0;
	for(size_t i=1;i<=_fem.NodeCount();++i)
	{
		_flu[index].NodeNum=i;
		_flu[index].x=it_node->x;
		if(ConsolidNode[i-1])
		{
			++index;
			_flu[index].NodeNum=i;
			_flu[index].x=it_node->x;
		}
		++it_node;
		++index;
	}
	double rtn=0;
	for(size_t i=0;i<_flu.size();++i)
	{
		for(size_t j=1;j<=_flu[i].flu_line.row();++j)
			rtn=max(rtn,abs(_flu[i].flu_line(j,2)));
	}
	return rtn*360;
}

void CalcShearLagAdv3U(Cmatrix<double>& _rst,const CFem& fem,const CeqSetD &eqsd,const CeqSetD&eqsd_rev,const size_t index,const size_t secnum,const size_t NodeNum,const double MaxM) //secmum,index,pos从1开始算,
{
	const size_t pos=IndexAndSecToPos(index,secnum);
	//const size_t NodeNum=IndexAndSecToNodeNum(index,secnum);
	_rst.resize(1,11);
	Cmatrix<double> rst_rev(1,6);//rst、rst_rev 每列依次为U1 U2 U3 U1' U2' U3' M λ  rou1 rou2 rou3
	size_t index_rev,secnum_rev;
	//CFem::const_iterator_node it=fem.NodeHead();
	//advance(it,NodeNum-1);
	//_rst(1,1)=it->x;
	_rst(1,1)=eqsd.CalcVarRst(1,0,index,secnum);
	_rst(1,2)=eqsd.CalcVarRst(2,0,index,secnum);
	_rst(1,3)=eqsd.CalcVarRst(3,0,index,secnum);
	_rst(1,4)=eqsd.CalcVarRst(1,1,index,secnum);
	_rst(1,5)=eqsd.CalcVarRst(2,1,index,secnum);
	_rst(1,6)=eqsd.CalcVarRst(3,1,index,secnum);

	secnum_rev=eqsd.SecDataCount.size()-secnum+1;
	index_rev=eqsd.SecDataCount[secnum-1]-index+1;
	
	rst_rev(1,1)=eqsd_rev.CalcVarRst(1,0,index_rev,secnum_rev);
	rst_rev(1,2)=eqsd_rev.CalcVarRst(2,0,index_rev,secnum_rev);
	rst_rev(1,3)=eqsd_rev.CalcVarRst(3,0,index_rev,secnum_rev);
	rst_rev(1,4)=-eqsd_rev.CalcVarRst(1,1,index_rev,secnum_rev);
	rst_rev(1,5)=-eqsd_rev.CalcVarRst(2,1,index_rev,secnum_rev);
	rst_rev(1,6)=-eqsd_rev.CalcVarRst(3,1,index_rev,secnum_rev);
	for(size_t i=1;i<=6;++i)
	{
		_rst(1,i)=_rst(1,i)/2+rst_rev(1,i)/2;
	}
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	if(index==1)
	{
		advance(it_rst,NodeNum-1);
		_rst(1,7)=(*it_rst)(1,5);
	}
	else
	{
		advance(it_rst,NodeNum-2);
		_rst(1,7)=(*it_rst)(2,5);
	}
	if(abs(_rst(1,7))/MaxM>1e-6)
	{
		_rst(1,8)=1.0-3.0/4.0*E/_rst(1,7)*(xIdI[pos-1].I1*_rst(1,4)+xIdI[pos-1].I2*_rst(1,5)+xIdI[pos-1].I3*_rst(1,6));
	}
	else
	{
		_rst(1,8)=_rst(1,7)=0;
		//Mzero.push_back(pos);
	}
	_rst(1,9)=CalcRou3U(_rst,pos,1);
	_rst(1,10)=CalcRou3U(_rst,pos,2);
	_rst(1,11)=CalcRou3U(_rst,pos,3);
}