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
	CCreateFeature(int fHigh, int fLow, int fliterNum);
	~CCreateFeature();
	bool FBankExtraction(string fList);
	vector<float> Fbank;	//存储Fbank特征
	int frameNum;
	
private:
	short *sample;	//语音数据
	struct Swavhead;	//语音数据头
	int filterNum;	//梅尔滤波器三角滤波器个数，等于输入的能带数
	int fHigh;	//截止频率上限
	int fLow;	//截止频率下限
	int FFTLen=256;	//FFT点数
	
	void PreEmphasise(const short *data, int len, float *out);	//对语音进行预加重
	void InitHamming(float *hamWin, int hamWinSize);	//初始化汉明窗
	void HammingWindow(float *buf, float *win, int frmLen, float* data);	//对每帧语音进行加窗
	void ComputeFFT(float *data, int frameSampleLen, vector<complex<float> > &vecList);	//语音信号做FFT
	void FFT(const unsigned long &ulN, vector<complex<float>> &vecList);	//FFT计算函数
	void CalPowSpectrum(vector<complex<float>> &vecList, vector<float> &powSpectrum); //计算功率谱
	void InitFilt(float **FiltWeight, int num_filt);	//初始化滤波器
	void CreateFilt(float **FiltWeight, int num_filt, int Fs, int high, int low);	//计算滤波器在每一点上的效果
	void Mel_EN(float **w, int num_filt, vector<float> &powSpectrum);	//信号做梅尔滤波函数
};