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
	//bool ReadWav(string fList, int fHigh);	//�����ļ������Լ���Ƶ��ֹƵ�ʣ����ڼ���Ƶ��ֹƵ���Ƿ�Ϸ�

private:
	short *sample;	//��������
	struct Swavhead;
	int filterNum=15;
	int fHigh;
	int fLow;
	int FFTLen=256;
	
	void PreEmphasise(const short *data, int len, float *out);	//����������Ԥ����
	void InitHamming(float *hamWin, int hamWinSize);	//��ʼ��������
	void HammingWindow(float *buf, float *win, int frmLen, float* data);	//��ÿ֡�������мӴ�
	void ComputeFFT(float *data, int frameSampleLen, vector<complex<float> > &vecList);	//��FFT
	void FFT(const unsigned long &ulN, vector<complex<float>> &vecList);
	void CalPowSpectrum(vector<complex<float>> &vecList, vector<float> &powSpectrum); 
	void InitFilt(float **FiltWeight, int num_filt);
	void CreateFilt(float **FiltWeight, int num_filt, int Fs, int high, int low);
};