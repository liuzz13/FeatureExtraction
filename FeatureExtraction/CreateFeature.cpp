#include"CreateFeature.h"

extern double Pi = 3.1415926536;
extern int frameLen = 25;	//֡��
extern int frameSpace = 10; //֡��


struct CCreateFeature::Swavhead
{
	char ChunkID[4];                    // "RIFF"��־  
	unsigned int ChunkSize;     // �ļ�����(WAVE�ļ��Ĵ�С, ����ǰ8���ֽ�)  
	char Format[4];                     // "WAVE"��־  
	char SubChunk1ID[4];                // "fmt "��־  
	unsigned int SubChunk1Size; // �����ֽ�(����)  
	unsigned short int AudioFormat;     // ��ʽ���(10HΪPCM��ʽ����������)  
	unsigned short int NumChannels;     // ͨ����(������Ϊ1, ˫����Ϊ2)  
	unsigned short int SampleRate;      // ������(ÿ��������), ��ʾÿ��ͨ���Ĳ����ٶ�  
	unsigned int ByteRate;          // ������Ƶ���ݴ�������, ��ֵΪ:ͨ����*ÿ������λ��*ÿ����������λ��/8  
	unsigned short int BlockAlign;      // ÿ����������λ��(���ֽ���), ��ֵΪ:ͨ����*ÿ����������λֵ/8  
	unsigned short int BitsPerSample;   // ÿ����������λ��, ��ʾÿ�������и�������������λ��.  
	char SubChunk2ID[4];                // ���ݱ��"data"  
	unsigned int SubChunk2Size; // �������ݵĳ��� 
};

CCreateFeature::CCreateFeature(int fHighIn, int fLowIn, int fliterNumIn)
{
	fHigh = fHighIn;
	fLow = fLowIn;
	filterNum = fliterNumIn;
	Fbank.clear();
}

CCreateFeature::~CCreateFeature()
{
	Fbank.clear();
}

bool CCreateFeature::FBankExtraction(string fList)
{
	//���ļ�
	FILE *wavfile = fopen(fList.c_str(), "rb");
	if (wavfile == NULL)
	{
		cout << "�ļ�" << fList << "�򿪴���";
		return 0;
	}


	//���ļ�ͷ��Ϣ
	Swavhead wav;
	fread(&wav, sizeof(struct Swavhead), 1, wavfile);

	if (fHigh > wav.SampleRate / 2)
	{
		cout << "��Ƶ��ֹƵ�����ô���";
		return 0;
	}
	
	//����Ƿ�Ϊwav�ļ�
	if (wav.Format == "WAVE")
	{
		cout << "�򿪲�Ϊwav�ļ�";
		return 0;
	}

	//��û��ָ����ֹƵ�����趨ΪȫƵ��
	if (fHigh == -1)
		fHigh = wav.SampleRate / 2;
	if (fLow == -1)
		fLow = 0;

	int sampleSize = wav.SubChunk2Size / wav.BlockAlign;	//����������
	frameNum = sampleSize / (wav.SampleRate / (1000 / frameSpace));	//������֡��
	int frameSampleLen = wav.SampleRate / (1000 / frameLen);	//ÿ֡�����Ĳ�������
	int frameShift = wav.SampleRate / (1000 / frameSpace);	//֡�Ʋ�������

	//������������ĩβ����һ֡�Ĳ���
	int frameActualLen = (frameNum - 1)*frameShift + frameSampleLen;
	sample = new short[frameActualLen];
	fread(sample, sizeof(short), sampleSize, wavfile);

	for (int i = sampleSize; i < frameActualLen; i++)
		sample[i] = 0;

	fclose(wavfile);

	float *sampleAfterPreEmphasise = new float[frameActualLen];
	PreEmphasise(sample, frameActualLen, sampleAfterPreEmphasise);	//����������Ԥ���أ��������sampleAfterPreEmphasise
	
	//��ʼ��������
	float *hamWin = new float[frameSampleLen];
	InitHamming(hamWin, frameSampleLen);

	//��ʼ��÷���˲���
	float **filtWeight = new float *[filterNum];
	for (int i = 0; i < filterNum; i++)
		filtWeight[i] = new float[FFTLen / 2 + 1];

	InitFilt(filtWeight, filterNum);
	CreateFilt(filtWeight, filterNum, wav.SampleRate, fHigh, fLow);

	//��ÿ֡�������д���
	vector<complex<float> > vecList;
	vector<float> powSpectrum;
	Fbank.clear();

	for (int i = 0; i < frameNum; i++)
	{
		float *data = new float[frameSampleLen];
		memcpy(data, &sampleAfterPreEmphasise[frameShift*i], frameSampleLen * 4);	//��ȡÿһ֡

		float *dataAfterWin = new float[frameSampleLen];
		HammingWindow(data, hamWin, frameSampleLen, dataAfterWin);	//�Ӵ�
		ComputeFFT(dataAfterWin, frameSampleLen, vecList);	//����FFT
		CalPowSpectrum(vecList, powSpectrum);	//����Ƶ�׼��㹦����
		Mel_EN(filtWeight, filterNum, powSpectrum);

		delete[]data;
		delete[]dataAfterWin;
		vecList.clear();
		powSpectrum.clear();
	}
		
	delete[] sample;
	delete[] sampleAfterPreEmphasise;
	delete[] hamWin;
	for (int i = 0; i < filterNum; i++)
		delete[] filtWeight[i];
	delete[] filtWeight;

	return 1;
}


/*
���������ݽ���Ԥ����
���룺
const short *data��Ԥ����ǰ����������
int len���洢�������ݵ������С
float *out�����ڷ���Ԥ���غ���
�������
*/
void CCreateFeature::PreEmphasise(const short *data, int len, float *out)//Ԥ����
{
	out[0] = data[0];
	for (int i = 1; i < len; i++)
		out[i] = float(data[i]) - 0.97 * float(data[i - 1]);
}

/*
��ʼ��������
���룺
double *hamWin�����غ���������
int hamWinSize����������С
�������
*/
void CCreateFeature::InitHamming(float *hamWin, int hamWinSize)
{
	for (int i = 0; i < hamWinSize; i++)
	{
		hamWin[i] = (float)(0.54 - 0.46 * cos(2 * Pi * (float)i / ((float)hamWinSize - 1)));
	}
}

/*
��ĳ֡�������мӴ�
���룺
double *dataIn��������������
double *win����
int frmLen��֡��
double* dataOut���Ӵ�����
�������
*/
void CCreateFeature::HammingWindow(float *dataIn, float *win, int frmLen, float* dataOut)
{
	for (int i = 0; i<frmLen; i++)
	{
		dataOut[i] = dataIn[i] * win[i];
	}
}


/*
���㸵��Ҷ���������������ݽ��в��������֮�����FFT��������ص�vecList��
���룺
float *data��һ֡����������
int frameSampleLen���������ݳ���
vector<complex<float>> &vecList������Ҷ����
*/
void CCreateFeature::ComputeFFT(float *data, int frameSampleLen, vector<complex<float> > &vecList)
{
	for (int i = 0; i<FFTLen ; ++i)
	{
		if (i<frameSampleLen)
		{
			complex<float> temp(data[i]);
			vecList.push_back(temp);
		}
		else
		{
			complex<float> temp(0);
			vecList.push_back(temp);
		}
	}
	FFT(frameSampleLen, vecList);
}
/*
����Ҷ�任
�ο���https://github.com/Linzecong/MFCC-DTW/blob/master/WaveFunction.cpp
*/
void CCreateFeature::FFT(const unsigned long &ulN, vector<complex<float> > &vecList)
{
	//�õ�����

	unsigned long ulPower = 0; //����
	unsigned long ulN1 = ulN - 1;
	while (ulN1 > 0)
	{
		ulPower++;
		ulN1 /= 2;
	}
	//����

	bitset<sizeof(unsigned long)* 8> bsIndex; //����������
	unsigned long ulIndex; //��ת������
	unsigned long ulK;
	for (unsigned long p = 0; p < ulN; p++)
	{
		ulIndex = 0;
		ulK = 1;
		bsIndex = bitset<sizeof(unsigned long)* 8>(p);
		for (unsigned long j = 0; j < ulPower; j++)
		{
			ulIndex += bsIndex.test(ulPower - j - 1) ? ulK : 0;
			ulK *= 2;
		}

		if (ulIndex > p)
		{
			complex<float> c = vecList[p];
			vecList[p] = vecList[ulIndex];
			vecList[ulIndex] = c;
		}
	}

	//������ת����

	vector<complex<float> > vecW;
	for (unsigned long i = 0; i < ulN / 2; i++)
	{
		vecW.push_back(complex<float>(cos(2 * i * Pi / ulN), -1 * sin(2 * i * Pi / ulN)));
	}

	//����FFT

	unsigned long ulGroupLength = 1; //�εĳ���
	unsigned long ulHalfLength = 0; //�γ��ȵ�һ��
	unsigned long ulGroupCount = 0; //�ε�����
	complex<float> cw; //WH(x)
	complex<float> c1; //G(x) + WH(x)
	complex<float> c2; //G(x) - WH(x)
	for (unsigned long b = 0; b < ulPower; b++)
	{
		ulHalfLength = ulGroupLength;
		ulGroupLength *= 2;
		for (unsigned long j = 0; j < ulN; j += ulGroupLength)
		{
			for (unsigned long k = 0; k < ulHalfLength; k++)
			{
				cw = vecW[k * ulN / ulGroupLength] * vecList[j + k + ulHalfLength];
				c1 = vecList[j + k] + cw;
				c2 = vecList[j + k] - cw;
				vecList[j + k] = c1;
				vecList[j + k + ulHalfLength] = c2;
			}
		}
	}
}

/*
��ʼ���˲���
���룺
float **FiltWeight��������ϵ������
int num_filt�������˲�������
*/

void CCreateFeature::InitFilt(float **FiltWeight, int num_filt)
{
	int i, j;
	for (i = 0; i<num_filt; i++)
		for (j = 0; j < FFTLen / 2 + 1; j++)
			FiltWeight[i][j] = 0.0;
}


/*
����÷���˲���ÿ�����ϵ��˲�Ч��
�ο���https://github.com/Linzecong/MFCC-DTW/blob/master/WaveFunction.cpp
���룺
float **w�����ڷ��صõ������˲���Ч��
int num_filt�������˲�������
int Fs������������
int high����ֹƵ������
int low����ֹƵ������
*/
void CCreateFeature::CreateFilt(float **w, int num_filt, int Fs, int high, int low)
{
	float df = (float)Fs / (float)FFTLen;    // FFT interval
	int indexlow = round((float)FFTLen*(float)low / (float)Fs); // FFT index of low freq limit
	int indexhigh = round((float)FFTLen*(float)high / (float)Fs); // FFT index of high freq limit

	float melmax = 2595.0*log10(1.0 + (float)high / 700.0); // mel high frequency
	float melmin = 2595.0*log10(1.0 + (float)low / 700.0);  // mel low frequency
	float melinc = (melmax - melmin) / (float)(num_filt + 1); //mel half bandwidth
	//float melcenters[num_filt];        // mel center frequencies
	//float fcenters[num_filt];          // Hertz center frequencies
	//int indexcenter[num_filt];         // FFT index for Hertz centers
	//int indexstart[num_filt];   //FFT index for the first sample of each filter
	//int indexstop[num_filt];    //FFT index for the last sample of each filter
	float *melcenters = new float[num_filt];
	float *fcenters = new float[num_filt];
	int *indexcenter = new int[num_filt];
	int *indexstart = new int[num_filt];
	int *indexstop = new int[num_filt];

	float increment, decrement; // increment and decrement of the left and right ramp
	float sum = 0.0;
	int i, j;
	for (i = 1; i <= num_filt; i++)
	{
		melcenters[i - 1] = (float)i*melinc + melmin;   // compute mel center frequencies
		fcenters[i - 1] = 700.0*(pow(10.0, melcenters[i - 1] / 2595.0) - 1.0); // compute Hertz center frequencies
		indexcenter[i - 1] = round(fcenters[i - 1] / df); // compute fft index for Hertz centers		 
	}
	for (i = 1; i <= num_filt - 1; i++)  // Compute the start and end FFT index of each channel
	{
		indexstart[i] = indexcenter[i - 1];
		indexstop[i - 1] = indexcenter[i];
	}
	indexstart[0] = indexlow;
	indexstop[num_filt - 1] = indexhigh;
	for (i = 1; i <= num_filt; i++)
	{
		increment = 1.0 / ((float)indexcenter[i - 1] - (float)indexstart[i - 1]); // left ramp
		for (j = indexstart[i - 1]; j <= indexcenter[i - 1]; j++)
			w[i - 1][j] = ((float)j - (float)indexstart[i - 1])*increment;
		decrement = 1.0 / ((float)indexstop[i - 1] - (float)indexcenter[i - 1]);    // right ramp
		for (j = indexcenter[i - 1]; j <= indexstop[i - 1]; j++)
			w[i - 1][j] = 1.0 - ((float)j - (float)indexcenter[i - 1])*decrement;
	}

	for (i = 1; i <= num_filt; i++)     // Normalize filter weights by sum
	{
		for (j = 1; j <= FFTLen / 2 + 1; j++)
			sum = sum + w[i - 1][j - 1];
		for (j = 1; j <= FFTLen / 2 + 1; j++)
			w[i - 1][j - 1] = w[i - 1][j - 1] / sum;
		sum = 0.0;
	}

	delete[]melcenters;
	delete[]fcenters;
	delete[]indexcenter;
	delete[]indexstart;
	delete[]indexstop;

}

/*
���㹦����
���룺
vector<complex<float>> &vecList��FFT��Ƶ��
vector<float> &powSpectrum�����ڷ��ؼ���õ��Ĺ�����
*/
void CCreateFeature::CalPowSpectrum(vector<complex<float>> &vecList, vector<float> &powSpectrum)
{
	int i;
	float temp;
	for (i = 0; i < FFTLen / 2 + 1; i++)
	{
		temp = vecList[i].real()*vecList[i].real() + vecList[i].imag()*vecList[i].imag();
		powSpectrum.push_back(temp);
	}
}

/*
����÷���˲�����
���룺
float **w��÷���˲����˲�Ȩ��
int num_filt�������˲�������
vector<float>& vec_mag��������
*/
void CCreateFeature::Mel_EN(float **w, int num_filt, vector<float>& vec_mag) // computes log energy of each channel
{
	int i, j;
	float *M_energy = new float[num_filt];
	for (i = 0; i < num_filt; i++)    // set initial energy value to 0
		M_energy[i] = 0.0F;

	for (i = 0; i < num_filt; i++)
	{
		for (j = 0; j < FFTLen / 2 + 1; j++)
			M_energy[i] = M_energy[i] + w[i][j] * vec_mag[j];
		M_energy[i] = (float)(log(M_energy[i]));
		Fbank.push_back(M_energy[i]);
	}

}