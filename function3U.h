extern double La1,La2;
extern vector<double> Diaphragm;
extern vector<size_t> DiagIndex;//隔板在keypos中的位置，从0开始算
inline double CalcRou3U(const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv=50);
inline double CalcRou3U(const Cmatrix<double>&_rst,const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv=50);

inline void GenFemMod3U(CFem &fem,const size_t ndiv)
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
	ofstream ofs("4.0fem.dat");
	fem.OutputData(ofs);
	ofs.close();
	fem.Calculate(sparse_matrix());
}

inline void GenFemMod3UShear(CFem &fem,const size_t ndiv,bool bcalc=true)
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
	if(!bcalc)	return;
	//边界
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

inline void GenEqSet3U(CeqSetD & eqsd,const CFem &fem)
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
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	ofstream ofs_I("I.csv");
	ofs_I<<"I1,I2,I3,Iw,I,dI1,dI2,dI3,dIw,I\n";
	for(size_t i=0;i<xIdI.size();++i)
	{
		ofs_I<<xIdI[i].I1<<","<<xIdI[i].I2<<","<<xIdI[i].I3<<","<<xIdI[i].Iw<<","<<xIdI[i].I<<",";
		ofs_I<<xIdI[i].dI1<<","<<xIdI[i].dI2<<","<<xIdI[i].dI3<<","<<xIdI[i].dIw<<","<<xIdI[i].dI<<endl;
	}
	ofs_I.close();
	eqsd.ReSizeCoef();
	{
		size_t index=0;
		const double b1=var_sec.plate1.back().x;
		const double b2=var_sec.plate2.back().x;
		const double b3=var_sec.plate3_zd.back().x;
		double dIs;
		//方程1
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index].dI-xIdI[index].dIw;
			mat(1,1)=9.0/16/xIdI[index].I*pow(xIdI[index].I1,2)-9.0/14*xIdI[index].I1;
			mat(1,2)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
			mat(1,3)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
			mat(1,4)=9.0/8*xIdI[index].I1*xIdI[index].dI1/xIdI[index].I-9.0/14*xIdI[index].dI1;
			mat(1,4)-=9.0/16*pow(xIdI[index].I1,2)*dIs/pow(xIdI[index].I,2);
			mat(1,5)=xIdI[index].I1*xIdI[index].dI2/xIdI[index].I+xIdI[index].I2*xIdI[index].dI1/xIdI[index].I;
			mat(1,5)-=xIdI[index].I2*xIdI[index].I1*dIs/pow(xIdI[index].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=xIdI[index].I1*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI1/xIdI[index].I;
			mat(1,6)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index].I1/E/pow(b1,2);
			mat(1,8)=-3.0/4*xIdI[index].I1/E/xIdI[index].I*(*it_rst)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index].dI1-xIdI[index].I1*xIdI[index].dI/xIdI[index].I)*(*it_rst)(1,5)/E/xIdI[index].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				dIs=xIdI[index].dI-xIdI[index].dIw;
				mat(j,1)=9.0/16/xIdI[index].I*pow(xIdI[index].I1,2)-9.0/14*xIdI[index].I1;
				mat(j,2)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
				mat(j,3)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
				mat(j,4)=9.0/8*xIdI[index].I1*xIdI[index].dI1/xIdI[index].I-9.0/14*xIdI[index].dI1;
				mat(j,4)-=9.0/16*pow(xIdI[index].I1,2)*dIs/pow(xIdI[index].I,2);
				mat(j,5)=xIdI[index].I1*xIdI[index].dI2/xIdI[index].I+xIdI[index].I2*xIdI[index].dI1/xIdI[index].I;
				mat(j,5)-=xIdI[index].I2*xIdI[index].I1*dIs/pow(xIdI[index].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=xIdI[index].I1*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI1/xIdI[index].I;
				mat(j,6)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index].I1/E/pow(b1,2);
				mat(j,8)=-3.0/4*xIdI[index].I1/E/xIdI[index].I*(*it_rst)(2,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index].dI1-xIdI[index].I1*xIdI[index].dI/xIdI[index].I)*(*it_rst)(2,5)/E/xIdI[index].I;
				++it_rst;
			}
			++index;
			eqsd.SetCoef(mat,1,i);
		}
		
		//方程2
		index=0;
		it_rst=fem.EForceHead();
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index].dI-xIdI[index].dIw;
			mat(1,1)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
			mat(1,2)=9.0/16*pow(xIdI[index].I2,2)/xIdI[index].I-9.0/14*xIdI[index].I2;
			mat(1,3)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
			
			mat(1,4)=xIdI[index].I2*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI2/xIdI[index].I;
			mat(1,4)-=xIdI[index].I1*xIdI[index].I2*dIs/pow(xIdI[index].I,2);
			mat(1,4)*=9.0/16;

			mat(1,5)=9.0/8*xIdI[index].I2*xIdI[index].dI2/xIdI[index].I-9.0/14*xIdI[index].dI2;
			mat(1,5)-=9.0/16*pow(xIdI[index].I2,2)*dIs/pow(xIdI[index].I,2);
			mat(1,6)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
			mat(1,6)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index].I2/E/pow(b2,2);
			mat(1,8)=-3.0/4*xIdI[index].I2/E/xIdI[index].I*(*it_rst)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index].dI2-xIdI[index].I2*xIdI[index].dI/xIdI[index].I)*(*it_rst)(1,5)/E/xIdI[index].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				dIs=xIdI[index].dI-xIdI[index].dIw;
				mat(j,1)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
				mat(j,2)=9.0/16*pow(xIdI[index].I2,2)/xIdI[index].I-9.0/14*xIdI[index].I2;
				mat(j,3)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
				mat(j,4)=xIdI[index].I2*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI2/xIdI[index].I;
				mat(j,4)-=xIdI[index].I1*xIdI[index].I2*dIs/pow(xIdI[index].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=9.0/8*xIdI[index].I2*xIdI[index].dI2/xIdI[index].I-9.0/14*xIdI[index].dI2;
				mat(j,5)-=9.0/16*pow(xIdI[index].I2,2)*dIs/pow(xIdI[index].I,2);
				mat(j,6)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
				mat(j,6)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index].I2/E/pow(b2,2);
				mat(j,8)=-3.0/4*xIdI[index].I2/E/xIdI[index].I*(*it_rst)(2,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index].dI2-xIdI[index].I2*xIdI[index].dI/xIdI[index].I)*(*it_rst)(2,5)/E/xIdI[index].I;
				++it_rst;
			}
			++index;
			eqsd.SetCoef(mat,2,i);
		}

		//方程3
		index=0;
		it_rst=fem.EForceHead();
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index].dI-xIdI[index].dIw;
			mat(1,1)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
			mat(1,2)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
			mat(1,3)=9.0/16*pow(xIdI[index].I3,2)/xIdI[index].I-9.0/14*xIdI[index].I3;	
			mat(1,4)=xIdI[index].I3*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI3/xIdI[index].I;
			mat(1,4)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,4)*=9.0/16;
			mat(1,5)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
			mat(1,5)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=9.0/8*xIdI[index].I3*xIdI[index].dI3/xIdI[index].I-9.0/14*xIdI[index].dI3;
			mat(1,6)-=9.0/16*pow(xIdI[index].I3,2)*dIs/pow(xIdI[index].I,2);
			mat(1,7)=9.0/5*G*xIdI[index].I3/E/pow(b3,2);
			mat(1,8)=-3.0/4*xIdI[index].I3/E/xIdI[index].I*(*it_rst)(1,3);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index].dI3-xIdI[index].I3*xIdI[index].dI/xIdI[index].I)*(*it_rst)(1,5)/E/xIdI[index].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				dIs=xIdI[index].dI-xIdI[index].dIw;
				mat(j,1)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
				mat(j,2)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
				mat(j,3)=9.0/16*pow(xIdI[index].I3,2)/xIdI[index].I-9.0/14*xIdI[index].I3;	
				mat(j,4)=xIdI[index].I3*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI3/xIdI[index].I;
				mat(j,4)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
				mat(j,5)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=9.0/8*xIdI[index].I3*xIdI[index].dI3/xIdI[index].I-9.0/14*xIdI[index].dI3;
				mat(j,6)-=9.0/16*pow(xIdI[index].I3,2)*dIs/pow(xIdI[index].I,2);
				mat(j,7)=9.0/5*G*xIdI[index].I3/E/pow(b3,2);
				mat(j,8)=-3.0/4*xIdI[index].I3/E/xIdI[index].I*(*it_rst)(2,3);  //Q=(*it_rst)(2,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index].dI3-xIdI[index].I3*xIdI[index].dI/xIdI[index].I)*(*it_rst)(2,5)/E/xIdI[index].I;
				++it_rst;
			}
			++index;
			eqsd.SetCoef(mat,3,i);
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

inline void GenEqSet3U(CeqSetD & eqsd,const CFem &fem,const Cmatrix<double>& elem_force)
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
	ofstream ofs_I("I.csv");
	ofs_I<<"I1,I2,I3,Iw,I,dI1,dI2,dI3,dIw,I\n";
	for(size_t i=0;i<xIdI.size();++i)
	{
		ofs_I<<xIdI[i].I1<<","<<xIdI[i].I2<<","<<xIdI[i].I3<<","<<xIdI[i].Iw<<","<<xIdI[i].I<<",";
		ofs_I<<xIdI[i].dI1<<","<<xIdI[i].dI2<<","<<xIdI[i].dI3<<","<<xIdI[i].dIw<<","<<xIdI[i].dI<<endl;
	}
	ofs_I.close();
	eqsd.ReSizeCoef();
	{
		size_t index=0;
		const double b1=var_sec.plate1.back().x;
		const double b2=var_sec.plate2.back().x;
		const double b3=var_sec.plate3_zd.back().x;
		double dIs;
		size_t elem_index=1;
		//方程1
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index].dI-xIdI[index].dIw;
			mat(1,1)=9.0/16/xIdI[index].I*pow(xIdI[index].I1,2)-9.0/14*xIdI[index].I1;
			mat(1,2)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
			mat(1,3)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
			mat(1,4)=9.0/8*xIdI[index].I1*xIdI[index].dI1/xIdI[index].I-9.0/14*xIdI[index].dI1;
			mat(1,4)-=9.0/16*pow(xIdI[index].I1,2)*dIs/pow(xIdI[index].I,2);
			mat(1,5)=xIdI[index].I1*xIdI[index].dI2/xIdI[index].I+xIdI[index].I2*xIdI[index].dI1/xIdI[index].I;
			mat(1,5)-=xIdI[index].I2*xIdI[index].I1*dIs/pow(xIdI[index].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=xIdI[index].I1*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI1/xIdI[index].I;
			mat(1,6)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index].I1/E/pow(b1,2);
			mat(1,8)=-3.0/4*xIdI[index].I1/E/xIdI[index].I*elem_force(elem_index,2);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index].dI1-xIdI[index].I1*xIdI[index].dI/xIdI[index].I)*elem_force(elem_index,1)/E/xIdI[index].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				dIs=xIdI[index].dI-xIdI[index].dIw;
				mat(j,1)=9.0/16/xIdI[index].I*pow(xIdI[index].I1,2)-9.0/14*xIdI[index].I1;
				mat(j,2)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
				mat(j,3)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
				mat(j,4)=9.0/8*xIdI[index].I1*xIdI[index].dI1/xIdI[index].I-9.0/14*xIdI[index].dI1;
				mat(j,4)-=9.0/16*pow(xIdI[index].I1,2)*dIs/pow(xIdI[index].I,2);
				mat(j,5)=xIdI[index].I1*xIdI[index].dI2/xIdI[index].I+xIdI[index].I2*xIdI[index].dI1/xIdI[index].I;
				mat(j,5)-=xIdI[index].I2*xIdI[index].I1*dIs/pow(xIdI[index].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=xIdI[index].I1*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI1/xIdI[index].I;
				mat(j,6)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index].I1/E/pow(b1,2);
				mat(j,8)=-3.0/4*xIdI[index].I1/E/xIdI[index].I*elem_force(elem_index,4);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index].dI1-xIdI[index].I1*xIdI[index].dI/xIdI[index].I)*elem_force(elem_index,3)/E/xIdI[index].I;
				++elem_index;
			}
			++index;
			eqsd.SetCoef(mat,1,i);
		}
		
		//方程2
		index=0;
		elem_index=1;
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index].dI-xIdI[index].dIw;
			mat(1,1)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
			mat(1,2)=9.0/16*pow(xIdI[index].I2,2)/xIdI[index].I-9.0/14*xIdI[index].I2;
			mat(1,3)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
			
			mat(1,4)=xIdI[index].I2*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI2/xIdI[index].I;
			mat(1,4)-=xIdI[index].I1*xIdI[index].I2*dIs/pow(xIdI[index].I,2);
			mat(1,4)*=9.0/16;

			mat(1,5)=9.0/8*xIdI[index].I2*xIdI[index].dI2/xIdI[index].I-9.0/14*xIdI[index].dI2;
			mat(1,5)-=9.0/16*pow(xIdI[index].I2,2)*dIs/pow(xIdI[index].I,2);
			mat(1,6)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
			mat(1,6)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,6)*=9.0/16;
			mat(1,7)=9.0/5*G*xIdI[index].I2/E/pow(b2,2);
			mat(1,8)=-3.0/4*xIdI[index].I2/E/xIdI[index].I*elem_force(elem_index,2);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index].dI2-xIdI[index].I2*xIdI[index].dI/xIdI[index].I)*elem_force(elem_index,1)/E/xIdI[index].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				dIs=xIdI[index].dI-xIdI[index].dIw;
				mat(j,1)=9.0/16*xIdI[index].I1*xIdI[index].I2/xIdI[index].I;
				mat(j,2)=9.0/16*pow(xIdI[index].I2,2)/xIdI[index].I-9.0/14*xIdI[index].I2;
				mat(j,3)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
				mat(j,4)=xIdI[index].I2*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI2/xIdI[index].I;
				mat(j,4)-=xIdI[index].I1*xIdI[index].I2*dIs/pow(xIdI[index].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=9.0/8*xIdI[index].I2*xIdI[index].dI2/xIdI[index].I-9.0/14*xIdI[index].dI2;
				mat(j,5)-=9.0/16*pow(xIdI[index].I2,2)*dIs/pow(xIdI[index].I,2);
				mat(j,6)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
				mat(j,6)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,6)*=9.0/16;
				mat(j,7)=9.0/5*G*xIdI[index].I2/E/pow(b2,2);
				mat(j,8)=-3.0/4*xIdI[index].I2/E/xIdI[index].I*elem_force(elem_index,4);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index].dI2-xIdI[index].I2*xIdI[index].dI/xIdI[index].I)*elem_force(elem_index,3)/E/xIdI[index].I;
				++elem_index;
			}
			++index;
			eqsd.SetCoef(mat,2,i);
		}

		//方程3
		index=0;
		elem_index=1;
		for(size_t i=1;i<=SecDataCount.size();++i)
		{
			mat.resize(SecDataCount[i-1],8);
			dIs=xIdI[index].dI-xIdI[index].dIw;
			mat(1,1)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
			mat(1,2)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
			mat(1,3)=9.0/16*pow(xIdI[index].I3,2)/xIdI[index].I-9.0/14*xIdI[index].I3;	
			mat(1,4)=xIdI[index].I3*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI3/xIdI[index].I;
			mat(1,4)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,4)*=9.0/16;
			mat(1,5)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
			mat(1,5)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
			mat(1,5)*=9.0/16;
			mat(1,6)=9.0/8*xIdI[index].I3*xIdI[index].dI3/xIdI[index].I-9.0/14*xIdI[index].dI3;
			mat(1,6)-=9.0/16*pow(xIdI[index].I3,2)*dIs/pow(xIdI[index].I,2);
			mat(1,7)=9.0/5*G*xIdI[index].I3/E/pow(b3,2);
			mat(1,8)=-3.0/4*xIdI[index].I3/E/xIdI[index].I*elem_force(elem_index,2);  //Q=(*it_rst)(1,3),但正方向和理论公式的剪力方向相反
			mat(1,8)+=3.0/4*(xIdI[index].dI3-xIdI[index].I3*xIdI[index].dI/xIdI[index].I)*elem_force(elem_index,1)/E/xIdI[index].I;
			
			for(size_t j=2;j<=SecDataCount[i-1];++j)
			{
				++index;
				dIs=xIdI[index].dI-xIdI[index].dIw;
				mat(j,1)=9.0/16*xIdI[index].I1*xIdI[index].I3/xIdI[index].I;
				mat(j,2)=9.0/16*xIdI[index].I2*xIdI[index].I3/xIdI[index].I;
				mat(j,3)=9.0/16*pow(xIdI[index].I3,2)/xIdI[index].I-9.0/14*xIdI[index].I3;	
				mat(j,4)=xIdI[index].I3*xIdI[index].dI1/xIdI[index].I+xIdI[index].I1*xIdI[index].dI3/xIdI[index].I;
				mat(j,4)-=xIdI[index].I1*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,4)*=9.0/16;
				mat(j,5)=xIdI[index].I2*xIdI[index].dI3/xIdI[index].I+xIdI[index].I3*xIdI[index].dI2/xIdI[index].I;
				mat(j,5)-=xIdI[index].I2*xIdI[index].I3*dIs/pow(xIdI[index].I,2);
				mat(j,5)*=9.0/16;
				mat(j,6)=9.0/8*xIdI[index].I3*xIdI[index].dI3/xIdI[index].I-9.0/14*xIdI[index].dI3;
				mat(j,6)-=9.0/16*pow(xIdI[index].I3,2)*dIs/pow(xIdI[index].I,2);
				mat(j,7)=9.0/5*G*xIdI[index].I3/E/pow(b3,2);
				mat(j,8)=-3.0/4*xIdI[index].I3/E/xIdI[index].I*elem_force(elem_index,4);  //Q=(*it_rst)(2,3),但正方向和理论公式的剪力方向相反
				mat(j,8)+=3.0/4*(xIdI[index].dI3-xIdI[index].I3*xIdI[index].dI/xIdI[index].I)*elem_force(elem_index,3)/E/xIdI[index].I;
				++elem_index;
			}
			++index;
			eqsd.SetCoef(mat,3,i);
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

void CalcShearLagAdv3U(const CFem& fem,const CeqSetD &eqsd)
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
	Cmatrix<double> du,du_rev,rst_rev;
	eqsd.CalcVarRst(1,du);
	eqsd_rev.CalcVarRst(1,du_rev);
	du_rev.revered();
	du_rev*=-1;

	rst.resize(du.row(),12);
	rst_rev.resize(du_rev.row(),3);
	eqsd.GetRst(rst,2);
	eqsd_rev.GetRst(rst_rev,1);
	rst_rev.revered();
	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,2)=rst(i,2)/2+rst_rev(i,1)/2;
		rst(i,3)=rst(i,3)/2+rst_rev(i,2)/2;
		rst(i,4)=rst(i,4)/2+rst_rev(i,3)/2;
	}

	CFem::const_iterator_node it_node=fem.NodeHead();
	CFem::const_iterator_rst it_rst=fem.EForceHead();
	size_t pos=1;
	//rst 每列依次为x U1 U2 U3 U1' U2' U3' M λ  rou1 rou2 rou3
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		for(size_t j=1;j<SecDataCount[i];++j)
		{
			rst(pos,1)=it_node->x;
			rst(pos,5)=du(pos,1)/2+du_rev(pos,1)/2;
			rst(pos,6)=du(pos,2)/2+du_rev(pos,2)/2;
			rst(pos,7)=du(pos,3)/2+du_rev(pos,3)/2;
			if(abs((*it_rst)(1,5))/MaxM>1e-6)
			{
				rst(pos,9)=1.0-3.0/4.0*E/(*it_rst)(1,5)*(xIdI[pos-1].I1*rst(pos,5)+xIdI[pos-1].I2*rst(pos,6)+xIdI[pos-1].I3*rst(pos,7));
				rst(pos,8)=(*it_rst)(1,5);
			}
			else
			{
				rst(pos,9)=rst(pos,8)=0;
				Mzero.push_back(pos);
			}
			++it_node;
			++pos;
			++it_rst;
		}
		--it_rst;
		rst(pos,1)=it_node->x;
		rst(pos,5)=du(pos,1)/2+du_rev(pos,1)/2;
		rst(pos,6)=du(pos,2)/2+du_rev(pos,2)/2;
		rst(pos,7)=du(pos,3)/2+du_rev(pos,3)/2;
		if(abs((*it_rst)(2,5))/MaxM>1e-6)
		{
			rst(pos,9)=1.0-3.0/4.0*E/(*it_rst)(2,5)*(xIdI[pos-1].I1*rst(pos,5)+xIdI[pos-1].I2*rst(pos,6)+xIdI[pos-1].I3*rst(pos,7));
			rst(pos,8)=(*it_rst)(2,5);
		}
		else
		{
			rst(pos,9)=rst(pos,8)=0;
			Mzero.push_back(pos);
		}
		++pos;
		++it_rst;
	}
	//计算翼缘平均厚度
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

	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,10)=CalcRou3U(i,1);
		rst(i,11)=CalcRou3U(i,2);
		rst(i,12)=CalcRou3U(i,3);
	}
}

void CalcShearLagAdv3U(const CFem& fem,const CeqSetD &eqsd,const Cmatrix<double>& elem_force)  //elem_force 1-4列分别为Mi Qi Mj Qj
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
	Cmatrix<double> du,du_rev,rst_rev;
	eqsd.CalcVarRst(1,du);
	eqsd_rev.CalcVarRst(1,du_rev);
	du_rev.revered();
	du_rev*=-1;

	rst.resize(du.row(),12);
	rst_rev.resize(du_rev.row(),3);
	eqsd.GetRst(rst,2);
	eqsd_rev.GetRst(rst_rev,1);
	rst_rev.revered();
	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,2)=rst(i,2)/2+rst_rev(i,1)/2;
		rst(i,3)=rst(i,3)/2+rst_rev(i,2)/2;
		rst(i,4)=rst(i,4)/2+rst_rev(i,3)/2;
	}

	CFem::const_iterator_node it_node=fem.NodeHead();
	size_t elem_index=1;
	size_t pos=1;
	//rst 每列依次为x U1 U2 U3 U1' U2' U3' M λ  rou1 rou2 rou3
	for(size_t i=0;i<SecDataCount.size();++i)
	{
		for(size_t j=1;j<SecDataCount[i];++j)
		{
			rst(pos,1)=it_node->x;
			rst(pos,5)=du(pos,1)/2+du_rev(pos,1)/2;
			rst(pos,6)=du(pos,2)/2+du_rev(pos,2)/2;
			rst(pos,7)=du(pos,3)/2+du_rev(pos,3)/2;
			if(abs(elem_force(elem_index,1))/MaxM>1e-6)
			{
				rst(pos,9)=1.0-3.0/4.0*E/elem_force(elem_index,1)*(xIdI[pos-1].I1*rst(pos,5)+xIdI[pos-1].I2*rst(pos,6)+xIdI[pos-1].I3*rst(pos,7));
				rst(pos,8)=elem_force(elem_index,1);
			}
			else
			{
				rst(pos,9)=rst(pos,8)=0;
				Mzero.push_back(pos);
			}
			++it_node;
			++pos;
			++elem_index;
		}
		--elem_index;
		rst(pos,1)=it_node->x;
		rst(pos,5)=du(pos,1)/2+du_rev(pos,1)/2;
		rst(pos,6)=du(pos,2)/2+du_rev(pos,2)/2;
		rst(pos,7)=du(pos,3)/2+du_rev(pos,3)/2;
		if(abs(elem_force(elem_index,3))/MaxM>1e-6)
		{
			rst(pos,9)=1.0-3.0/4.0*E/elem_force(elem_index,3)*(xIdI[pos-1].I1*rst(pos,5)+xIdI[pos-1].I2*rst(pos,6)+xIdI[pos-1].I3*rst(pos,7));
			rst(pos,8)=elem_force(elem_index,3);
		}
		else
		{
			rst(pos,9)=rst(pos,8)=0;
			Mzero.push_back(pos);
		}
		++pos;
		++elem_index;
	}
	//计算翼缘平均厚度
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

	for(size_t i=1;i<=rst.row();++i)
	{
		rst(i,10)=CalcRou3U(i,1);
		rst(i,11)=CalcRou3U(i,2);
		rst(i,12)=CalcRou3U(i,3);
	}
}

inline double CalcSecStress3U(const size_t index,const double y,const double z,const size_t plate_index)
{
	double bi;
	if(plate_index==1) bi=var_sec.plate1.back().x;
	else if(plate_index==2) bi=var_sec.plate2.back().x;
	else bi=var_sec.plate3_zd.back().x;
	return E*z*(rst(index,8)/E/xIdI[index-1].I+(1.0-pow(y/bi,3))*rst(index,4+plate_index)-3.0/4.0/xIdI[index-1].I
		*(xIdI[index-1].I1*rst(index,5)+xIdI[index-1].I2*rst(index,6)+xIdI[index-1].I3*rst(index,7)));
}

inline double CalcSecStress3U(const Cmatrix<double>&_rst,const size_t index,const double y,const double z,const size_t plate_index)
{
	double bi;
	if(plate_index==1) bi=var_sec.plate1.back().x;
	else if(plate_index==2) bi=var_sec.plate2.back().x;
	else bi=var_sec.plate3_zd.back().x;
	return E*z*(_rst(1,7)/E/xIdI[index-1].I+(1.0-pow(y/bi,3))*_rst(1,3+plate_index)-3.0/4.0/xIdI[index-1].I
		*(xIdI[index-1].I1*_rst(1,4)+xIdI[index-1].I2*_rst(1,5)+xIdI[index-1].I3*_rst(1,6)));
}

inline double CalcStressLineInt3U(const size_t pos/*截面索引*/,const double y,const size_t index/*翼缘板标志1-3*/,
	const double L,const size_t ndiv=50) //竖向线积分
{
	/*const double L=CalcLength(y,index);*/
	const double delta=L/ndiv;
	double Zbeg=xIdI[pos-1].C;
	if(index==1||index==2) Zbeg-=xIdI[pos-1].H;
	else Zbeg-=L;
	vector<double> StressRst(ndiv+1);
	for(size_t i=0;i<=ndiv;++i)
	{
		StressRst[i]=CalcSecStress3U(pos,y,Zbeg,index);
		Zbeg+=delta;
	}
	return LineInt(StressRst,delta);
}

inline double CalcStressLineInt3U(const Cmatrix<double>& _rst,const size_t pos/*截面索引*/,const double y,const size_t index/*翼缘板标志1-3*/,
	const double L,const size_t ndiv=50) //竖向线积分
{
	/*const double L=CalcLength(y,index);*/
	const double delta=L/ndiv;
	double Zbeg=xIdI[pos-1].C;
	if(index==1||index==2) Zbeg-=xIdI[pos-1].H;
	else Zbeg-=L;
	vector<double> StressRst(ndiv+1);
	for(size_t i=0;i<=ndiv;++i)
	{
		StressRst[i]=CalcSecStress3U(_rst,pos,y,Zbeg,index);
		Zbeg+=delta;
	}
	return LineInt(StressRst,delta);
}

inline double CalcRou3U(const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv)
{
	//计算翼缘面积
	double L=0,bf;
	vector<CPlateGeom> plate;
	if(index==1) {L=La1;plate=var_sec.plate1;}
	else if(index==2)  {L=La2;plate=var_sec.plate2;}
	else  
	{
		for(size_t i=0;i<plate3[pos-1].size()-1;++i)
		{
			L+=(plate3[pos-1][i+1].x-plate3[pos-1][i].x)/2*(plate3[pos-1][i+1].thickness+plate3[pos-1][i].thickness);
		}
		L/=plate3[pos-1].back().x;
		plate=plate3[pos-1];
	}
	bf=plate.back().x;
	//计算翼缘平均厚度
	//面积分
	const double delta=bf/ndiv;
	vector<double> StrLineInt(ndiv+1);
	double y=0;
	for(size_t i=0;i<ndiv;++i)
	{
		StrLineInt[i]=CalcStressLineInt3U(pos,y,index,CalcLength(y,plate),1);
		y+=delta;
	}
	StrLineInt[ndiv]=CalcStressLineInt3U(pos,bf,index,CalcLength(bf,plate),1);
	double AreaInt=LineInt(StrLineInt,delta);
	double LineInt=CalcStressLineInt3U(pos,bf,index,L,1);
	return AreaInt/LineInt/bf;
}

inline double CalcRou3U(const Cmatrix<double>&_rst,const size_t pos/*截面索引*/,const size_t index/*翼缘板标志1-3*/,const size_t ndiv)
{
	//计算翼缘面积
	double L=0,bf;
	vector<CPlateGeom> plate;
	if(index==1) {L=La1;plate=var_sec.plate1;}
	else if(index==2)  {L=La2;plate=var_sec.plate2;}
	else  
	{
		for(size_t i=0;i<plate3[pos-1].size()-1;++i)
		{
			L+=(plate3[pos-1][i+1].x-plate3[pos-1][i].x)/2*(plate3[pos-1][i+1].thickness+plate3[pos-1][i].thickness);
		}
		L/=plate3[pos-1].back().x;
		plate=plate3[pos-1];
	}
	bf=plate.back().x;
	//计算翼缘平均厚度
	//面积分
	const double delta=bf/ndiv;
	vector<double> StrLineInt(ndiv+1);
	double y=0;
	for(size_t i=0;i<ndiv;++i)
	{
		StrLineInt[i]=CalcStressLineInt3U(_rst,pos,y,index,CalcLength(y,plate),1);
		y+=delta;
	}
	StrLineInt[ndiv]=CalcStressLineInt3U(_rst,pos,bf,index,CalcLength(bf,plate),1);
	double AreaInt=LineInt(StrLineInt,delta);
	double LineInt=CalcStressLineInt3U(_rst,pos,bf,index,L,1);
	return AreaInt/LineInt/bf;
}