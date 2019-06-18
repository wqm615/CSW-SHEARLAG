//适用于变高度梁,不考虑腹板处刚度,分布竖向荷载、集中竖向荷载作用下的计算
//考虑隔板
#include <iostream>
#include <fstream>
#include "DIM.h"
#include "CeqSetD.h"
#include <string>
#include <fem\fem.h>
#include <iomanip>
#include "StructDef.h"
#include "functionV2.h"
#include "function3U.h"

vector<CBd > Bd;  //全局边界条件
vector<CContLoad > ContLoad;  //全局集中荷载
vector<CUnifLoad> UnifLoad;
vector<double> Spans;
vector<double> Diaphragm;
vector<int> MidPiers,SidePiers(2);
double L=0;
double E;
VarSecBeam var_sec;
vector<CSecPropty2> xIdI;
vector<vector<CPlateGeom> >plate3;
vector<size_t> SecDataCount,Mzero;
vector<double> deltas;
Cmatrix<double> rst;
double La1=0,La2=0;  //翼缘板1、2的平均厚度
vector<size_t> DiagIndex;//隔板在keypos中的位置，从0开始算
Ccsw csw;

int main()
{
	cout<<"请输入数据文件(扩展名为txt,勿输入扩展名):";
	char name[1000]="测试4.2";
	cin.getline(name,1000);
	string FileName(name);;
	ifstream ifs(FileName+".txt");
	if(!ifs)
	{
		cout<<"文件打开失败!\n";
		system("pause");
		return 0;
	}
	size_t ndiv;   //有限元网格及差分法步长尺寸
	ifs_ignore(ifs,3);
	ifs>>ndiv;
	//有限元模型的输入、建立与计算
	if(	!InputUnifSec(ifs)) return 0;
	ofstream ofs(FileName+"_orig.dat");
	OutputOrig(ofs);
	ofs.close();
	//return 0;
	CFem fem;
	GenFemMod3U(fem,ndiv);
	cout<<"有限元计算完成!\n";
	//差分方程的建立
	CeqSetD eqsd;
	GenEqSet3U(eqsd,fem);
	cout<<"差分计算完成!\n";
	//剪力滞系数相关结果的整理
	CalcShearLagAdv3U(fem,eqsd);
	//计算输出
	Output3U(FileName,fem,eqsd);
	system("pause");
}

