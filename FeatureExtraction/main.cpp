/*

*/

#include <iostream>
#include <fstream>
#include"CreateFeature.h"
using namespace std;


/*
��������������ļ���  ����ļ���  �ļ��б�  �ܴ���  ��ֹƵ��ѡ���ֹƵ�ʿɲ��
FeatureExtraction.exe inDir outDir wavlist.txt 15 60 3400

*/
int main(int argc, char* argv[])
{
	//�������淶
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
			cout << "�ܴ����������";
			system("pause");
			return 0;
		}
		if (fLow == 0 || fHigh == 0 || fLow>=fHigh)
		{
			cout << "��ֹƵ��ѡ���������";
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
			cout << "�ܴ����������";
			system("pause");
			return 0;
		}
	}
	else
	{
		cout << "���������������";
		system("pause");
		return 0;
	}

	//����������ȡ
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

		//д���ļ�
		ostringstream featurefile;
		featurefile << outputDir << '\\' << line << ".txt";
		FILE *outputfile = fopen(featurefile.str().c_str(), "wb");
		//дͷ�ļ�
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