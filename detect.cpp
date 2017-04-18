/**
*  Created by Faisal Sikder
*  University of Miami
*  This code is part of the step counter
*  Created 2013
*/

#include "detect.h"


void detect::removeInitialBias(float *data, const int points){

	int i, n;

	float  average = 0, sum = 0;

	n = points;

	for (i = 0; i < 5; i++) {

		sum = sum + data[i];

	}

	average = sum / 5;

	for (i = 0; i < n; i++){

		data[i] = data[i] - average;

	}
}

void detect::dc_shift(float *data, const int points,int factor){

	int i, n;

	float  average = 0, sum = 0;

	n = points;

	for (i = 0; i < n; i++) {

		sum = sum + data[factor*i];

	}

	average = sum / n;

	for (i = 0; i < n; i++){

		data[factor*i] = data[factor*i] - average;

	}
}


void detect::four1(float *data, const int points, const int isign) {

#define SWAP(a,b)tempr=(a);(a)=(b);(b)=tempr

	//the complex array is real+complex so the array 

	//as a size n = 2* number of complex samples

	// real part is the data[index] and 

	//the complex part is the data[index+1]

	int n = 2 * points, i, j, m, mmax, istep;

	float theta, wtemp, wpr, wpi, wr, wi, tempr, tempi;

	float pi = 6.28318530717959 / 2;



	//binary inversion (note that the indexes 

	//start from 0 witch means that the

	//real part of the complex is on the even-indexes 

	//and the complex part is on the odd-indexes

	j = 0;

	for (i = 0; i<n / 2; i += 2) {

		if (j > i) {

			//swap the real part

			SWAP(data[j], data[i]);

			//swap the complex part

			SWAP(data[j + 1], data[i + 1]);

			// checks if the changes occurs in the first half

			// and use the mirrored effect on the second half

			if ((j / 2)<(n / 4)){

				//swap the real part

				SWAP(data[(n - (i + 2))], data[(n - (j + 2))]);

				//swap the complex part

				SWAP(data[(n - (i + 2)) + 1], data[(n - (j + 2)) + 1]);

			}

		}

		m = n / 2;

		while (m >= 2 && j >= m) {

			j -= m;

			m = m / 2;

		}

		j += m;

	}

	//Danielson-Lanzcos routine 

	mmax = 2;

	//external loop

	while (n > mmax)

	{

		istep = mmax << 1;

		theta = isign*(2 * pi / mmax);

		wtemp = sin(0.5*theta);

		wpr = -2.0*wtemp*wtemp;

		wpi = sin(theta);

		wr = 1.0;

		wi = 0.0;

		//internal loops

		for (m = 1; m<mmax; m += 2) {

			for (i = m; i <= n; i += istep) {

				j = i + mmax;

				tempr = wr*data[j - 1] - wi*data[j];

				tempi = wr*data[j] + wi*data[j - 1];

				data[j - 1] = data[i - 1] - tempr;

				data[j] = data[i] - tempi;

				data[i - 1] += tempr;

				data[i] += tempi;

			}

			wr = (wtemp = wr)*wpr - wi*wpi + wr;

			wi = wi*wpr + wtemp*wpi + wi;

		}

		mmax = istep;

	}

}
int detect::readDataSize(string filename){
	vector <float> waistAccXaxis;

	int fc = 1, fStop = 1, samplingRate = 50;

	ifstream myfile;

	myfile.open(filename.c_str());

	if (!myfile.is_open())  // check file is open, quit if not

	{

		std::cerr << "failed to open file\n";

	}

	else {

		float number = 0;

		while (myfile >> number){    //

			waistAccXaxis.push_back(number);

		}

	}

	myfile.close();
	return  waistAccXaxis.size();
}

int detect::readData(string filename, float *inputData){
	int i = 0;

	ifstream myfile;

	myfile.open(filename.c_str());

	if (!myfile.is_open())  // check file is open, quit if not

	{
		std::cerr << "failed to open file\n";
	}

	else {

		float number = 0;
		
		while (myfile >> number){    //
			inputData[i] = number;
			i++;

		}

	}

	myfile.close();
	
	return i;
}


void detect::normalizeData(float *data1, int numData){
	float maxr, maxi, mini, minr,min;

	maxr = -1000, minr = 1000;

	maxi = -1000, mini = 1000;

	

	for (int i = 0; i < numData; i++){

		if (data1[2 * i] > maxr) maxr = data1[2 * i];

		if (data1[2 * i] < minr) minr = data1[2 * i];

		if (data1[2 * i + 1] > maxi) maxi = data1[2 * i + 1];

		if (data1[2 * i + 1] < mini) mini = data1[2 * i + 1];



	}

	min = 1000;

	for (int i = 0; i < numData; i++){

		if (data1[2 * i]> 0) data1[2 * i] = 2 * data1[2 * i] / ((maxr - minr));

		if (data1[2 * i] < 0) {

			data1[2 * i] = 2 * data1[2 * i] / ((maxr - minr));

			if (min > data1[2 * i]) min = data1[2 * i];

		}

		// cout << "min =" << min << endl;

		if (data1[2 * i + 1] > 0) data1[2 * i + 1] = 2 * data1[2 * i + 1] / ((maxi - mini));

		if (data1[2 * i + 1] < 0) data1[2 * i + 1] = 2 * data1[2 * i + 1] / ((maxi - mini));

	}
	//return 0;
}

void detect::doFftAndInvFft(float *data1, float *inputData, int samplingRate, int numData, int fc, int frequency){


	for (int i = 0; i < numData; i++){

		data1[2 * i] = inputData[i];

		data1[2 * i + 1] = 0;

	}



	for (int i = numData; i < fc; i++){

		data1[2 * i] = 0;

		data1[2 * i + 1] = 0;

	}



	four1(data1, fc, 1);

	int fStop = (fc / samplingRate) * frequency;

	for (int i = fStop; i < (fc - fStop); i++){
		data1[2 * i] = 0;
		data1[2 * i + 1] = 0;
	}

	four1(data1, fc, -1);
}



void detect::detectMidSwing(int startPos, int endPos, float *data1, float *data2, int numData){
	int i = 0, k;
	float minVal = 1000, maxVal = -1000, minTh, maxTh, tempPick = 0, tempSwigPick;

	for (i = 0; i < numData; i++){

		if (data1[2 * i] < minVal){
			minVal = data1[2 * i];
		}
		if (data1[2 * i] > maxVal){
			maxVal = data1[2 * i];
		}
	}
	minTh = minVal*.4;
	maxTh = maxVal*.5;
	i = 0;
	while (i<numData)
	{
		if (data1[2 * i] >= maxTh){
			tempPick = data1[2 * i];
			while (tempPick <= data1[2 * i])
			{
				tempPick = data1[2 * i];
				i++;
			}
			data2[2 * (i - 1)] = -2;
			while (data1[2 * i] >= 0)
			{
				i++;
			}
		}

		if (data1[2 * i] <= minTh){
			tempPick = data1[2 * i];
			while (tempPick >= data1[2 * i])
			{
				tempPick = data1[2 * i];
				i++;
			}
			data2[2 * (i - 1)] = -1;

			while (data1[2 * i] <= minTh){
				i++;
			}
		}
		i++;
	}

}

void detect::detectStartAndEnd(float *data1, float *data2, int startPos ,int endPos){
	int i = 0, k, zeroPos, minTh = startPos, maxTh = 0;
	float minVal = 1000, maxVal = -1000;
	//detect start position
	bool found = false;

	//get the min position
	for (i = startPos; i <= endPos; i++){
		if (data1[2 * i] < minVal){
			minVal = data1[2 * i];
			minTh = i;
		}
	}
	//cout << minVal << endl;
	//get max position on the right
	for (i = minTh; i>startPos; i--){
		if (data1[2 * i] > maxVal){
			maxVal = data1[2 * i];
			maxTh = i;
		}
	}

	//now get the zero
	i = minTh;
	
	while (i > maxTh)
	{
		if (data1[2 * i] >= 0){
			zeroPos = i;
			break;
		}
		i--;
	}
	
	//now detect the left starting position
	
	i = minTh-2;
	while (i >= zeroPos){
		if (((data1[2 * (i - 1)] - data1[2 * i]) <= .001)&&(data1[i * 2] >= minVal*.3)){
			data2[2 * i] = 2;
			found = true;
			break;
		}
		i--;
	}
	if (!found){
		data2[2 * zeroPos] = 2;
	}
	maxVal = -1000;
	maxTh = 0;
	//get the right max postion
	for (i = minTh; i < endPos; i++){
		if (data1[2 * i] > maxVal){
			maxVal = data1[2 * i];
			maxTh = i;
		}
	}

	//now get the right zero position
	i = minTh;
	zeroPos = data1[2 * (endPos - 2)];
	while (i < maxTh)
	{
		if (data1[2 * i] >= 0){
			zeroPos = i;
			break;
		}
		i++;
	}

	//cout << maxTh << " " << minTh << " " << endPos << endl;
	//right position
	found = false;
	i = minTh+2;
	while (i <= zeroPos){
		if (((data1[2 * (i + 1)] - data1[2 * i]) <= .001)&&(data1[i*2]>=minVal*.3)){
			data2[2 * i] = 2;
			found = true;
			break;
		}
		i++;
	}
	if (!found){
		data2[2 * zeroPos] = 2;
	}
}

void detect::detectTurn(int startPost,int endPost, float *data1, float *data2, int numData){
	int minTh=0, maxTh=0,i;
	float minVal = 1000, maxVal = -1000;

	int loopEnd = ((endPost - startPost) / 2 + startPost);
	//detect 1st turn
	for (i = startPost; i <= loopEnd; i++){
		
		if (data1[2 * i] < minVal){
			minVal = data1[2 * i];
			minTh = i;
		}
		if (data1[2 * i] > maxVal){
			maxVal = data1[2 * i];
			maxTh = i;
		}
	}

	if (maxTh > minTh){
		detectTurnHillView(startPost,endPost,data1, data2, numData, minTh, maxTh, minVal, maxVal);
	}
	else{
		detectTurnPitView(startPost, endPost, data1, data2, numData, minTh, maxTh, minVal, maxVal);
	}

}

void detect::detectTurnPitView(int startPost, int endPost, float *data1, float *data2, int numData, int minTh, int maxTh, float minVal, float maxVal){
	int i = 0, k;
	float tempPick = 0,zeroPos=0;

	i = maxTh ;
	tempPick = 1000;
	while (i<minTh)
	{
		if (data1[2 * i] <= 0){
			zeroPos = i;
			break;
		}
		i++;
	}
	//data2[(int)zeroPos * 2] = 3;

	i = zeroPos-1;
	while (i >= maxTh){
		if ((data1[2 * (i - 1)] - data1[2 * i]) <= .001){
			data2[2 * i] = 2;
			break;
		}
		i--;
	}

	i = zeroPos+1;
	while (i <= minTh){
		if ((data1[2 * i] - data1[2 * (i + 1)]) <= .001){
			data2[2 * i] = 2;
			break;
		}
		i++;
	}


	minVal = 100;
	maxVal = -100;
	maxTh = 0;
	minTh = 0;
	int loopStart = ((endPost - startPost) / 2 + startPost)-20;
	//detect 2nd turn
	for (i = loopStart; i < endPost; i++){
		
		if (data1[2 * i] < minVal){
			minVal = data1[2 * i];
			minTh = i;
		}
		if (data1[2 * i] > maxVal){
			maxVal = data1[2 * i];
			maxTh = i;
		}
	}
	
	
	i = minTh;
	tempPick = 1000;
	while (i<maxTh)
	{
		if (data1[2 * i] >= 0){
			zeroPos = i;
			break;
		}
		i++;
	}
	//data2[(int)zeroPos * 2] = 3;
	i = zeroPos - 1;
	while (i >= minTh){
		if ((data1[2 * i] - data1[2 * (i - 1)]) < .001){
			data2[2 * i] = 2;
			break;
		}
		i--;
	}

	i = zeroPos + 1;
	while (i <= maxTh){
		if ((data1[2 * (i + 1)] - data1[2 * i]) <= .001){
			data2[2 * i] = 2;
			break;
		}
		i++;
	}
	
}

void detect::detectTurnHillView(int startPost, int endPost, float *data1, float *data2, int numData, int minTh, int maxTh, float minVal, float maxVal){
	int i = 0, k;
	float tempPick = 0, zeroPos = 0;


	i = minTh;
	tempPick = 1000;
	while (i<maxTh)
	{
		if (data1[2 * i] >= 0){
			zeroPos = i;
			break;
		}
		i++;
	}
	//data2[(int)zeroPos * 2] = 3;

	i = zeroPos - 1;
	while (i >= minTh){
		if ((data1[2 * i] - data1[2 * (i - 1)]) <= .001){
			data2[2 * i] = 2;
			break;
		}
		i--;
	}

	i = zeroPos + 1;
	while (i <= maxTh){
		if ((data1[2 * (i+1)] - data1[2 * i]) <= .001){
			data2[2 * i] = 2;
			break;
		}
		i++;
	}


	minVal = 100;
	maxVal = -100;
	maxTh = 0;
	minTh = 0;
	int loopStart = ((endPost - startPost) / 2 + startPost) - 20;
	//detect 2nd turn
	for (i = loopStart / 2; i < endPost; i++){
		
		if (data1[2 * i] < minVal){
			minVal = data1[2 * i];
			minTh = i;
		}
		if (data1[2 * i] > maxVal){
			maxVal = data1[2 * i];
			maxTh = i;
		}
	}


	i = maxTh;
	while (i<minTh)
	{
		if (data1[2 * i] <= 0){
			zeroPos = i;
			break;
		}
		i++;
	}
	//data2[(int)zeroPos * 2] = 3;
	i = zeroPos - 1;
	while (i >= maxTh){
		if ((data1[2 * (i-1)] - data1[2 * i]) < .001){
			data2[2 * i] = 2;
			break;
		}
		i--;
	}

	i = zeroPos + 1;
	while (i <= minTh){
		if ((data1[2 * i] - data1[2 * (i+1)]) <= .001){
			data2[2 * i] = 2;
			break;
		}
		i++;
	}

}