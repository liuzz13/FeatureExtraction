/*

*/

#include <iostream>
#include <fstream>
#include <io.h>
#include <windows.h>
#include"CreateFeature.h"
using namespace std;


/*
��������������ļ���  ����ļ���  �ļ��б�  �ܴ���  ��ֹƵ��ѡ���ֹƵ�ʿɲ��
FeatureExtraction.exe inDir outDir wavlist.txt 15 60 3400
*/

bool CreatDir(string inputDir, string outputDir);

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

	CreatDir(inputDir, outputDir);	//����ļ��½�����Ӧ���ļ���
	ifstream file(wavlistFile);
	string line; 
	while (getline(file, line))
	{
		//����������ȡ
		CCreateFeature CreateFeature(fHigh, fLow, energyBandNum);
		ostringstream wavfile;
		wavfile << inputDir << '\\' << line;
		bool flag = CreateFeature.FBankExtraction(wavfile.str());
		if (flag == 0)
		{
			system("pause");
			return 0;
		}

		//д���ļ�
		ostringstream featurefile;
		featurefile << outputDir << '\\' << line << ".fbank";
		FILE *outputfile = fopen(featurefile.str().c_str(), "wb");
		//дͷ�ļ�
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
		//��ȡָ���ļ��������е��ļ���
		while (_findnext(hFile, &fileinfo) == 0)
		{
			//��ȡ���е��ļ�������
			if (fileinfo.attrib == _A_SUBDIR)
				dir.push_back(fileinfo.name);
		}
	}
	_findclose(hFile);

	//������ļ����д�����Ӧ�����ļ���
	for (int i = 1; i < dir.size(); i++)
	{
		string dirName = outputDir + "\\" + dir[i];
		if (CreateDirectory(dirName.c_str(), NULL) == 0)
		{
			cout << "�����ļ���ʧ��";
			return 0;
		}
	}

	return 1;
}