/**
*  Created by Faisal Sikder
*  University of Miami
*  This code is part of the step counter
*  Created 2013
*/
#include "detect.h"
#include "positiondetect.h"
#include <algorithm>    // std::sort

vector<positiondetect> getStartEndPosition(string file,string outfile);
vector<positiondetect> getTurnPosition(int startPos, int endPos, string file, string outfile);
vector<positiondetect> getSteps(int startPos, int endPos, string file, string outfile);


bool sort_function(positiondetect i, positiondetect j) { return (i.position<j.position); }

int main() {

	detect obj;
	int fc = 1, fStop = 1, samplingRate = 50, frequency = 4, numData=0;
	float *data1, *data2, *inputData,*inputData2;
	string outFile;
	ofstream fout;
	
	vector<positiondetect> startEnd = getStartEndPosition("WaistAccXaxis.txt","start-end-Waist-Acc-X.txt");
	vector<positiondetect> startEndY = getStartEndPosition("WaistAccYaxis.txt", "start-end-Waist-Acc-Y.txt");
	

	int startPosition = startEnd[0].position;
	int endPosition = startEnd[startEnd.size() - 1].position;

	vector<positiondetect> steps = getSteps(startPosition, endPosition,"AnkelGyroZaxis.txt", "steps-Ankel-Gyro-Z.txt");

	vector<positiondetect> turn = getTurnPosition(startPosition+100,endPosition,"WaistMagYaxis.txt","turn-Waist-Mag-Y.txt");
	vector<positiondetect> turnX = getTurnPosition(startPosition+100,endPosition,"WaistMagXaxis.txt", "turn-Waist-Mag-X.txt");
	//determind turn and turnX
	int turnEnd = 0;
	if ((turn.size() == turnX.size()) && turn.size() == 4){
		for (int i = 0; i<turn.size(); i++){
			if (abs(turn[i].position - turnX[i].position) >= 4){
				turn[i].position = (turn[i].position + turnX[i].position) / 2 - 5;
			}
		}
		turnEnd = turn[3].position;
	}
	//start-end-position
	if ((startEnd.size() == startEndY.size()) && startEnd.size() == 4){
		for (int i = 0; i<startEnd.size(); i++){
			if (abs(startEnd[i].position - startEndY[i].position) >= 8){
				startEnd[i].position = (startEnd[i].position + startEndY[i].position) / 2;
			}
		}
		if (startEnd[2].position < turnEnd){
			startEnd[2].position = turnEnd;
		}
	}
	
	startPosition = startEnd[0].position;
	endPosition = startEnd[startEnd.size() - 1].position;
	/*
	outFile = "allvalue-together.txt";
	vector<positiondetect> allValue;

	fout.open(outFile.c_str(), ios::out);

	for (int i = 0; i < startEnd.size(); i++){
		positiondetect obj = startEnd[i];
		fout << obj.position << " " << obj.value<<endl;
		allValue.push_back(startEnd[i]);
	}

	for (int i = 0; i < steps.size(); i++){
		positiondetect obj = steps[i];
		fout << obj.position << " " << obj.value << endl;
		allValue.push_back(steps[i]);
	}
	for (int i = 0; i < turn.size(); i++){
		positiondetect obj = turn[i];
		fout << obj.position << " " << obj.value << endl;
		allValue.push_back(turn[i]);
	}
	fout.close();
	*/
	outFile = "allvalue-together-cont.txt";

	//std::sort(allValue.begin(), allValue.end(), sort_function);

	fout.open(outFile.c_str(), ios::out);
	int j = 0;
	bool found = false;
	for (int i = 0; i <=(endPosition+10); i++){
		//fout << i<<" ";
		if (i == startPosition){
			fout << "2 ";
		}
		else if (i == endPosition){
			fout << "2 ";
		}
		else{
			fout << "0 ";
		}
		
		found = false;
		for (j = 0; j < startEnd.size(); j++){
			if (i == startEnd[j].position){
				fout << "5 ";
				found = true;
			}
		}
		if (!found){
			fout << "3 ";
		}

		found = false;
		for (j = 0; j < turn.size(); j++){
			if (i == turn[j].position){
				fout << "8 ";
				found = true;
			}
		}
		if (!found){
			fout << "6 ";
		}

		found = false;
		for (j = 0; j < steps.size(); j++){
			if (i == steps[j].position){
				fout << "11";
				found = true;
			}
		}
		if (!found){
			fout << "9";
		}
		fout << endl;
	}
	
	fout.close();

	return 0;
}

vector<positiondetect> getStartEndPosition(string infile,string outfile){

	detect obj;
	int fc = 1, fStop = 1, samplingRate = 50, frequency = 4, numData = 0;
	float *data1, *data2, *inputData, *inputData2;
	string outFile;
	ofstream fout;

	//first get starting point using Waist ACC X and Y
	numData = obj.readDataSize(infile);
	inputData = new float[numData];
	numData = obj.readData(infile, inputData);


	obj.removeInitialBias(inputData, numData);


	while (fc < numData){

		fc = fc * 2;

	}


	data1 = new float[2 * fc];
	frequency = 2;
	//do FFT and then Inv FFT
	obj.doFftAndInvFft(data1, inputData, samplingRate, numData, fc, frequency);

	obj.dc_shift(data1, numData, 2);
	//normalize fourier data
	obj.normalizeData(data1, numData);

	data2 = new float[2 * fc];
	for (int i = 0; i <= numData; i++){
		data2[i * 2] = 0;
	}
	obj.detectStartAndEnd(data1, data2, 0, numData / 2);
	obj.detectStartAndEnd(data1, data2, numData / 2, numData - 1);

	


	fout.open(outfile.c_str(), ios::out);
	vector<positiondetect> startEnd;
	for (int i = 0; i < numData; i++){

		fout << data1[2 * i] << " " << data2[2 * i] << endl;
		if (data2[2 * i] == 2){
			positiondetect temp(i, 5);
			startEnd.push_back(temp);
		}
	}

	fout.close();
	return startEnd;
}

vector<positiondetect> getTurnPosition(int startPos, int endPos, string infile, string outfile){
	detect obj;
	int fc = 1, fStop = 1, samplingRate = 50, frequency = 4, numData = 0;
	float *data1, *data2, *inputData, *inputData2;
	string outFile;
	ofstream fout;
	//detect turn
	frequency = 1;
	
	fc = 1;
	numData = 0;
	numData = obj.readDataSize(infile);
	inputData = new float[numData];
	numData = obj.readData(infile, inputData);

	obj.removeInitialBias(inputData, numData);


	while (fc < numData){

		fc = fc * 2;

	}


	data1 = new float[2 * fc];

	//do FFT and then Inv FFT
	obj.doFftAndInvFft(data1, inputData, samplingRate, numData, fc, frequency);

	//do DC shift again
	obj.dc_shift(data1, numData, 2);

	//normalize fourier data
	obj.normalizeData(data1, numData);


	data2 = new float[2 * fc];
	for (int i = 0; i < numData; i++){
		data2[i * 2] = 0;
	}
	obj.detectTurn(startPos,endPos, data1, data2, numData);

	

	fout.open(outfile.c_str(), ios::out);
	vector<positiondetect> turn;

	for (int i = 0; i < numData; i++){

		fout << data1[2 * i] << " " << data2[2 * i] << endl;
		if (data2[2 * i] != 0){
			positiondetect temp(i, 8);
			turn.push_back(temp);
		}

	}

	fout.close();
	return turn;
}

vector<positiondetect> getSteps(int startPos, int endPos, string infile, string outfile){
	detect obj;
	int fc = 1, fStop = 1, samplingRate = 50, frequency = 4, numData = 0;
	float *data1, *data2, *inputData, *inputData2;
	string outFile;
	ofstream fout;
	frequency = 8;
	fc = 1;
	numData = 0;
	numData = obj.readDataSize(infile);
	inputData = new float[numData];
	numData = obj.readData(infile, inputData);

	obj.removeInitialBias(inputData, numData);


	while (fc < numData){

		fc = fc * 2;

	}


	data1 = new float[2 * fc];

	//do FFT and then Inv FFT
	obj.doFftAndInvFft(data1, inputData, samplingRate, numData, fc, frequency);

	//normalize fourier data
	obj.normalizeData(data1, numData);

	data2 = new float[2 * fc];
	for (int i = 0; i <= numData; i++){
		data2[i * 2] = 0;
	}

	obj.detectMidSwing(startPos,endPos, data1, data2, numData);

	//outFile = "fft-post-data-a-g-z.txt";

	fout.open(outfile.c_str(), ios::out);
	
	vector<positiondetect> steps;
	int gotPick = 0;
	int firstFound = 0;
	for (int i = 0; i < numData; i++){

		fout << data1[2 * i] << " " << data2[2 * i] << endl;
		if (data2[2 * i] != 0){
			if (data2[2 * i] == -2){
				gotPick = 2;
				firstFound = 0;
			}

			if (data2[2 * i] == -1 && gotPick == 2 && firstFound==0){
				firstFound = i;
				gotPick--;
			}else if (data2[2 * i] == -1 && gotPick==1){
				//toe off
				positiondetect temp2(firstFound, 11);
				steps.push_back(temp2);
				//heal strike 
				positiondetect temp(i, 11);
				steps.push_back(temp);
				gotPick--;
				firstFound = 0;
			}
		}
	}

	fout.close();
	return steps;
}