#include<iostream>
#include<cstring>
#include<math.h>
#include<vector>
#include <complex>
#include <bitset>
using namespace std;

class CCreateFeature
{
public:
	bool FBankExtraction(string fList);
	CCreateFeature(int fHigh, int fLow, int fliterNum);
	//bool ReadWav(string fList, int fHigh);	//传入文件名，以及高频截止频率，用于检测高频截止频率是否合法

private:
	short *sample;	//语音数据
	struct Swavhead;
	int filterNum=15;
	int fHigh;
	int fLow;
	int FFTLen=256;
	
	void PreEmphasise(const short *data, int len, float *out);	//对语音进行预加重
	void InitHamming(float *hamWin, int hamWinSize);	//初始化汉明窗
	void HammingWindow(float *buf, float *win, int frmLen, float* data);	//对每帧语音进行加窗
	void ComputeFFT(float *data, int frameSampleLen, vector<complex<float> > &vecList);	//做FFT
	void FFT(const unsigned long &ulN, vector<complex<float>> &vecList);
	void CalPowSpectrum(vector<complex<float>> &vecList, vector<float> &powSpectrum); 
	void InitFilt(float **FiltWeight, int num_filt);
	void CreateFilt(float **FiltWeight, int num_filt, int Fs, int high, int low);
};