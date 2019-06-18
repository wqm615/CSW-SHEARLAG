//�����ڱ�߶���,�����Ǹ��崦�ն�,�ֲ�������ء�����������������µļ���
//���Ǹ���
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

vector<CBd > Bd;  //ȫ�ֱ߽�����
vector<CContLoad > ContLoad;  //ȫ�ּ��к���
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
double La1=0,La2=0;  //��Ե��1��2��ƽ�����
vector<size_t> DiagIndex;//������keypos�е�λ�ã���0��ʼ��
Ccsw csw;

int main()
{
	cout<<"�����������ļ�(��չ��Ϊtxt,��������չ��):";
	char name[1000]="����4.2";
	cin.getline(name,1000);
	string FileName(name);;
	ifstream ifs(FileName+".txt");
	if(!ifs)
	{
		cout<<"�ļ���ʧ��!\n";
		system("pause");
		return 0;
	}
	size_t ndiv;   //����Ԫ���񼰲�ַ������ߴ�
	ifs_ignore(ifs,3);
	ifs>>ndiv;
	//����Ԫģ�͵����롢���������
	if(	!InputUnifSec(ifs)) return 0;
	ofstream ofs(FileName+"_orig.dat");
	OutputOrig(ofs);
	ofs.close();
	//return 0;
	CFem fem;
	GenFemMod3U(fem,ndiv);
	cout<<"����Ԫ�������!\n";
	//��ַ��̵Ľ���
	CeqSetD eqsd;
	GenEqSet3U(eqsd,fem);
	cout<<"��ּ������!\n";
	//������ϵ����ؽ��������
	CalcShearLagAdv3U(fem,eqsd);
	//�������
	Output3U(FileName,fem,eqsd);
	system("pause");
}

