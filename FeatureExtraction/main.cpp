/*

*/

#include<iostream>
#include"CreateFeature.h"
using namespace std;


int fLow = 60;	//截止频率下限
int fHigh = 3400;	//截止频率上限



/*
输入参数：输入文件夹  输出文件夹  文件列表  能带数  截止频率选项（截止频率可不填）
FeatureExtraction.exe inDir outDir wavlist.txt 15 60 3400

*/
int main(int argc, char* argv[])
{
	////检查输入规范
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
	//		cout << "能带数输入错误";
	//		return 0;
	//	}
	//	if (fLow == 0 || fHigh == 0 || fLow>=fHigh)
	//	{
	//		cout << "截止频率选项输出错误";
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
	//		cout << "能带数输入错误";
	//		return 0;
	//	}
	//}
	//else
	//{
	//	cout << "输入参数个数有误";
	//	return 0;
	//}
	int fHigh = 3400;
	int fLow = 60;
	int energyBandNum = 15;

	CCreateFeature CreateFeature(fHigh, fLow, energyBandNum);

	CreateFeature.FBankExtraction("E:/培训/标注培训/广东数据标注/eng167.wav");





	return 0;
}