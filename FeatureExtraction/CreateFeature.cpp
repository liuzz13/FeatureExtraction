
#include"CreateFeature.h"

extern double Pi = 3.1415926536;
extern int frameLen = 25;	//帧长
extern int frameSpace = 10; //帧移


struct CCreateFeature::Swavhead
{
	char ChunkID[4];                    // "RIFF"标志  
	unsigned int ChunkSize;     // 文件长度(WAVE文件的大小, 不含前8个字节)  
	char Format[4];                     // "WAVE"标志  
	char SubChunk1ID[4];                // "fmt "标志  
	unsigned int SubChunk1Size; // 过渡字节(不定)  
	unsigned short int AudioFormat;     // 格式类别(10H为PCM格式的声音数据)  
	unsigned short int NumChannels;     // 通道数(单声道为1, 双声道为2)  
	unsigned short int SampleRate;      // 采样率(每秒样本数), 表示每个通道的播放速度  
	unsigned int ByteRate;          // 波形音频数据传输速率, 其值为:通道数*每秒数据位数*每样本的数据位数/8  
	unsigned short int BlockAlign;      // 每样本的数据位数(按字节算), 其值为:通道数*每样本的数据位值/8  
	unsigned short int BitsPerSample;   // 每样本的数据位数, 表示每个声道中各个样本的数据位数.  
	char SubChunk2ID[4];                // 数据标记"data"  
	unsigned int SubChunk2Size; // 语音数据的长度 
};

CCreateFeature::CCreateFeature(int fHigh, int fLow, int fliterNum)
{
	fHigh = fHigh;
	fLow = fLow;
	fliterNum = fliterNum;
}


bool CCreateFeature::FBankExtraction(string fList)
{
	//读文件
	//string fList = "E:/培训/标注培训/广东数据标注/eng167.wav";
	FILE *wavfile = fopen(fList.c_str(), "rb");

	//读文件头信息
	Swavhead wav;
	fread(&wav, sizeof(struct Swavhead), 1, wavfile);

	if (fHigh > wav.SampleRate / 2)
	{
		cout << "高频截止频率设置错误";
		return 0;
	}

	int sampleSize = wav.SubChunk2Size / wav.BlockAlign;	//采样点总数
	int frameNum = sampleSize / (wav.SampleRate / (1000 / frameSpace));	//语音总帧数
	int frameSampleLen = wav.SampleRate / (1000 / frameLen);	//每帧语音的采样点数

	sample = new short[sampleSize];
	fread(sample, sizeof(short), sampleSize, wavfile);	//读入整个语音数据

	fclose(wavfile);


	float *sampleAfterPreEmphasise = new float[sampleSize];
	PreEmphasise(sample, sampleSize, sampleAfterPreEmphasise);	//对语音进行预加重，结果存入sampleAfterPreEmphasise

	//初始化汉明窗
	float *hamWin = new float[frameSampleLen];
	InitHamming(hamWin, frameSampleLen);

	//初始化梅尔滤波器
	float **filtWeight = new float *[filterNum];
	for (int i = 0; i < filterNum; i++)
		filtWeight[i] = new float[FFTLen / 2 + 1];

	InitFilt(filtWeight, filterNum);
	CreateFilt(filtWeight, filterNum, wav.SampleRate, fHigh, fLow);

	//对每帧语音进行处理
	vector<complex<float> > vecList;
	vector<float> powSpectrum;
	for (int i = 0; i < frameNum; i++)
	{
		float *data = new float[frameSampleLen];
		//if(memcpy(data, &sampleAfterPreEmphasise[frameSampleLen*i], frameSampleLen * wav.BlockAlign)==NULL);	//提取一帧的信息
		//	break;
		memcpy(data, &sampleAfterPreEmphasise[frameSampleLen*i], frameSampleLen * wav.BlockAlign);

		float *dataAfterWin = new float[frameSampleLen];
		HammingWindow(data, hamWin, frameSampleLen, dataAfterWin);	//加窗
		ComputeFFT(dataAfterWin, frameSampleLen, vecList);	//计算FFT
		CCreateFeature::CalPowSpectrum(vecList, powSpectrum);	//根据频谱计算功率谱


		delete[]data;
		delete[]dataAfterWin;
	}

	vecList.clear();
	delete []sample;
	delete []sampleAfterPreEmphasise;

	return 1;
}


/*
对语音数据进行预加重
输入：const short *data：预加重前语音数据流
	  int len：存储语音数据的数组大小
	  float *out：用于返回预加重后结果
输出：无
*/
void CCreateFeature::PreEmphasise(const short *data, int len, float *out)//预加重
{
	for (int i = len - 1; i <= 1; i--)
	{
		out[i] = float(data[i]) - 0.97 * float(data[i - 1]);
	}
	out[0] = data[0];
}

/*
初始化汉明窗
输入：double *hamWin：返回汉明窗数据
	  int hamWinSize：汉明窗大小
输出：无
*/
void CCreateFeature::InitHamming(float *hamWin, int hamWinSize)
{
	for (int i = 0; i < hamWinSize; i++)
	{
		hamWin[i] = (float)(0.54 - 0.46 * cos(2 * Pi * (float)i / ((float)hamWinSize - 1)));
	}
}

/*
对某帧语音进行加窗
输入：double *dataIn：输入语音数据
	  double *win：窗
	  int frmLen：帧长
	  double* dataOut：加窗后结果
输出：无
*/
void CCreateFeature::HammingWindow(float *dataIn, float *win, int frmLen, float* dataOut)
{
	for (int i = 0; i<frmLen; i++)
	{
		dataOut[i] = dataIn[i] * win[i];
	}
}


/*
计算傅里叶参数，主要是对语音数据进行补零操作，之后进行FFT，结果返回到vecList中
输入：float *data：一帧的语音数据
	  int frameSampleLen：语音数据长度
	  vector<complex<float>> &vecList：傅里叶参数
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
傅里叶变换
参考自https://github.com/Linzecong/MFCC-DTW/blob/master/WaveFunction.cpp
*/
void CCreateFeature::FFT(const unsigned long &ulN, vector<complex<float> > &vecList)
{
	//得到幂数

	unsigned long ulPower = 0; //幂数
	unsigned long ulN1 = ulN - 1;
	while (ulN1 > 0)
	{
		ulPower++;
		ulN1 /= 2;
	}
	//反序

	bitset<sizeof(unsigned long)* 8> bsIndex; //二进制容器
	unsigned long ulIndex; //反转后的序号
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

	//计算旋转因子

	vector<complex<float> > vecW;
	for (unsigned long i = 0; i < ulN / 2; i++)
	{
		vecW.push_back(complex<float>(cos(2 * i * Pi / ulN), -1 * sin(2 * i * Pi / ulN)));
	}

	//计算FFT

	unsigned long ulGroupLength = 1; //段的长度
	unsigned long ulHalfLength = 0; //段长度的一半
	unsigned long ulGroupCount = 0; //段的数量
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


void CCreateFeature::CalPowSpectrum(vector<complex<float>> &vecList, vector<float> &powSpectrum)
{
	int i;
	float temp;
	for (i = 1; i <= FFTLen / 2 + 1; i++)
	{
		temp = vecList[i - 1].real()*vecList[i - 1].real() + vecList[i - 1].imag()*vecList[i - 1].imag();
		powSpectrum.push_back(temp);
	}
}


void CCreateFeature::InitFilt(float **FiltWeight, int num_filt)
{
	int i, j;
	for (i = 0; i<num_filt; i++)
	for (j = 0; j<FFTLen / 2 + 1; j++)
		*(*(FiltWeight + i) + j) = 0.0;
}


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