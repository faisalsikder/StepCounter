#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

using namespace std;

class detect{
private:
public:
	detect(){};
	void removeInitialBias(float *data, const int points);
	void four1(float *data, const int points, const int isign);
	int readDataSize(string filename);
	int readData(string filename, float *data);
	void normalizeData(float *data1, int size);
	void doFftAndInvFft(float *data1, float *inputData, int samplingRate, int numData, int fc, int frequency);
	void detectMidSwing(int startPos, int endPos, float *data1, float *data2, int size);
	void detectStartAndEnd(float *data1, float *data2, int startpos,int endpos);
	void detectTurn(int startPos, int endPos, float *data1, float *data2, int size);
	void dc_shift(float *data, const int points,int factor);
	void detectTurnHillView(int startPost, int endPost, float *data1, float *data2, int size, int minTh, int maxTh, float minVal, float maxVal);
	void detectTurnPitView(int startPost, int endPost, float *data1, float *data2, int size, int minTh, int maxTh, float minVal, float maxVal);
};