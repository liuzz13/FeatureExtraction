/*

*/

#include<iostream>
#include"CreateFeature.h"
using namespace std;


int fLow = 60;	//��ֹƵ������
int fHigh = 3400;	//��ֹƵ������



/*
��������������ļ���  ����ļ���  �ļ��б�  �ܴ���  ��ֹƵ��ѡ���ֹƵ�ʿɲ��
FeatureExtraction.exe inDir outDir wavlist.txt 15 60 3400

*/
int main(int argc, char* argv[])
{
	////�������淶
	//if (argc == 6)
	//{
	//	string inputDir = argv[0];
	//	string outputDir = argv[1];
	//	string wavlistFile = argv[2];
	//	int energyBandNum = atoi(argv[3]);
	//	int fLow = atoi(argv[4]);
	//	int fHigh = atoi(argv[5]);
	//	if (energyBandNum == 0)
	//	{
	//		cout << "�ܴ����������";
	//		return 0;
	//	}
	//	if (fLow == 0 || fHigh == 0 || fLow>=fHigh)
	//	{
	//		cout << "��ֹƵ��ѡ���������";
	//		return 0;
	//	}
	//}
	//else if (argc == 4)
	//{
	//	string inputDir = argv[0];
	//	string outputDir = argv[1];
	//	string wavlistFile = argv[2];
	//	int energyBandNum = atoi(argv[3]);
	//	if (energyBandNum == 0)
	//	{
	//		cout << "�ܴ����������";
	//		return 0;
	//	}
	//}
	//else
	//{
	//	cout << "���������������";
	//	return 0;
	//}
	int fHigh = 3400;
	int fLow = 60;
	int energyBandNum = 15;

	CCreateFeature CreateFeature(fHigh, fLow, energyBandNum);

	CreateFeature.FBankExtraction("E:/��ѵ/��ע��ѵ/�㶫���ݱ�ע/eng167.wav");





	return 0;
}