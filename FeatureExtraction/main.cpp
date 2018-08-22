/*

*/

#include <iostream>
#include <fstream>
#include"CreateFeature.h"
using namespace std;


/*
输入参数：输入文件夹  输出文件夹  文件列表  能带数  截止频率选项（截止频率可不填）
FeatureExtraction.exe inDir outDir wavlist.txt 15 60 3400

*/
int main(int argc, char* argv[])
{
	//检查输入规范
	string inputDir;
	string outputDir;
	string wavlistFile;
	int energyBandNum;
	int fLow;
	int fHigh;
	if (argc == 7)
	{
		inputDir = argv[1];
		outputDir = argv[2];
		wavlistFile = argv[3];
		energyBandNum = atoi(argv[4]);
		fLow = atoi(argv[5]);
		fHigh = atoi(argv[6]);
		if (energyBandNum == 0)
		{
			cout << "能带数输入错误";
			system("pause");
			return 0;
		}
		if (fLow == 0 || fHigh == 0 || fLow>=fHigh)
		{
			cout << "截止频率选项输出错误";
			system("pause");
			return 0;
		}
	}
	else if (argc == 5)
	{
		inputDir = argv[1];
		outputDir = argv[2];
		wavlistFile = argv[3];
		energyBandNum = atoi(argv[4]);
		fLow = -1;
		fHigh = -1;
		if (energyBandNum == 0)
		{
			cout << "能带数输入错误";
			system("pause");
			return 0;
		}
	}
	else
	{
		cout << "输入参数个数有误";
		system("pause");
		return 0;
	}

	//进行特征提取
	CCreateFeature CreateFeature(fHigh, fLow, energyBandNum);

	ifstream file(wavlistFile);
	string line; 
	while (getline(file, line))
	{
		ostringstream wavfile;
		wavfile << inputDir << '\\' << line;
		bool temp = CreateFeature.FBankExtraction(wavfile.str());
		if (temp == 0)
		{
			system("pause");
			return 0;
		}

		//写入文件
		ostringstream featurefile;
		featurefile << outputDir << '\\' << line << ".txt";
		FILE *outputfile = fopen(featurefile.str().c_str(), "wb");
		//写头文件
		int frameNum = CreateFeature.frameNum;
		int frameShift = 100000;
		short filterSize = energyBandNum * 4;
		short flag = 73;
		fwrite(&frameNum, sizeof(int), 1, outputfile);
		fwrite(&frameShift, sizeof(int), 1, outputfile);
		fwrite(&filterSize, sizeof(short), 1, outputfile);
		fwrite(&flag, sizeof(short), 1, outputfile);
			
	}

 	system("pause");
	return 0;
}