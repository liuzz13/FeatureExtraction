/*

*/

#include <iostream>
#include <fstream>
#include <io.h>
#include <windows.h>
#include"CreateFeature.h"
using namespace std;


/*
输入参数：输入文件夹  输出文件夹  文件列表  能带数  截止频率选项（截止频率可不填）
FeatureExtraction.exe inDir outDir wavlist.txt 15 60 3400
*/

bool CreatDir(string inputDir, string outputDir);

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

	CreatDir(inputDir, outputDir);	//输出文件下建立对应的文件夹
	ifstream file(wavlistFile);
	string line; 
	while (getline(file, line))
	{
		//进行特征提取
		CCreateFeature CreateFeature(fHigh, fLow, energyBandNum);
		ostringstream wavfile;
		wavfile << inputDir << '\\' << line;
		bool flag = CreateFeature.FBankExtraction(wavfile.str());
		if (flag == 0)
		{
			system("pause");
			return 0;
		}

		//写入文件
		ostringstream featurefile;
		featurefile << outputDir << '\\' << line << ".fbank";
		FILE *outputfile = fopen(featurefile.str().c_str(), "wb");
		//写头文件
		int frameNum = CreateFeature.frameNum;
		int frameShift = 100000;
		short filterSize = energyBandNum * 4;
		short mark = 73;
		fwrite(&frameNum, sizeof(int), 1, outputfile);
		fwrite(&frameShift, sizeof(int), 1, outputfile);
		fwrite(&filterSize, sizeof(short), 1, outputfile);
		fwrite(&mark, sizeof(short), 1, outputfile);
		cout << CreateFeature.Fbank.size();
		for (int i = 0; i < CreateFeature.Fbank.size(); i++)
		{
			float temp = CreateFeature.Fbank[i];
			fwrite(&temp, sizeof(float), 1, outputfile);
		}
			
		//fwrite(&CreateFeature.Fbank, sizeof(float), CreateFeature.Fbank.size(), outputfile);
		fclose(outputfile);
	}

 	system("pause");
	return 0;
}


bool CreatDir(string inputDir, string outputDir)
{
	//string path = inputDir;
	struct _finddata_t fileinfo;
	long hFile = 0;
	string p;
	vector<string> name;
	vector<string> dir;
	if ((hFile = _findfirst(p.assign(inputDir).append("\\*").c_str(), &fileinfo)) != -1)
	{
		//获取指定文件夹下所有的文件名
		while (_findnext(hFile, &fileinfo) == 0)
		{
			//获取所有的文件夹名称
			if (fileinfo.attrib == _A_SUBDIR)
				dir.push_back(fileinfo.name);
		}
	}
	_findclose(hFile);

	//在输出文件夹中创建对应的子文件夹
	for (int i = 1; i < dir.size(); i++)
	{
		string dirName = outputDir + "\\" + dir[i];
		if (CreateDirectory(dirName.c_str(), NULL) == 0)
		{
			cout << "创建文件夹失败";
			return 0;
		}
	}

	return 1;
}