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
	vector<float> Fbank;	//�洢Fbank����
	int frameNum;
	
private:
	short *sample;	//��������
	struct Swavhead;	//��������ͷ
	int filterNum;	//÷���˲��������˲�������������������ܴ���
	int fHigh;	//��ֹƵ������
	int fLow;	//��ֹƵ������
	int FFTLen=256;	//FFT����
	
	void PreEmphasise(const short *data, int len, float *out);	//����������Ԥ����
	void InitHamming(float *hamWin, int hamWinSize);	//��ʼ��������
	void HammingWindow(float *buf, float *win, int frmLen, float* data);	//��ÿ֡�������мӴ�
	void ComputeFFT(float *data, int frameSampleLen, vector<complex<float> > &vecList);	//�����ź���FFT
	void FFT(const unsigned long &ulN, vector<complex<float>> &vecList);	//FFT���㺯��
	void CalPowSpectrum(vector<complex<float>> &vecList, vector<float> &powSpectrum); //���㹦����
	void InitFilt(float **FiltWeight, int num_filt);	//��ʼ���˲���
	void CreateFilt(float **FiltWeight, int num_filt, int Fs, int high, int low);	//�����˲�����ÿһ���ϵ�Ч��
	void Mel_EN(float **w, int num_filt, vector<float> &powSpectrum);	//�ź���÷���˲�����
};