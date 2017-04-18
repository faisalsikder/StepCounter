////////////////////////////////////
/**
*  Created by Faisal Sikder
*  University of Miami
*  This code is part of the step counter
*  Created 2013
*/


#include <iostream>
#include <fstream>

//using namespace std;

//1
int readData(float a[],int size){
	int i = 0;
	ifstream fin;
	string file = "data.txt";
	//Open the files
	fin.open(file.c_str(), ios::in);
	
	if (!fin) {
		cout << "Could Not Open Input File" << endl;
		return -1;
	}
	
	i = 0;
	while (!fin.eof() && i<size){
		fin >> a[i];
		i++;
	}
	
	fin.close();
	
	return i;
}

void writeintofile(float b[], int size,string output){
	//string output = "datahaar.txt";
	ofstream fout;

	fout.open(output.c_str(), ios::out);

	for (int i = 0; i<size; i++){
		fout << b[i] << endl;
	}
	fout.close();
}


int main(){
	int size = 512, i, j, half = size / 2;
	int itr = log2(size);
	float *a = new float[size];
	float *b = new float[size];
	float twosqrt = sqrt(2.0);
	i=readData(a,size);
	while (itr>0){
		i = 0;
		j = 0;
		
		while (i < size){
			b[j] = (a[i] + a[i + 1]) / twosqrt;
			b[half + j] = (a[i] - a[i + 1]) / twosqrt;
			i += 2;
			j++;
		}
		for (i = 0; i < size; i++){
			a[i] = b[i];
		}
		itr--;
	}
	writeintofile(b, size, "datahaar-mid.txt");
	for (i = 0; i<half; i++){
		a[i] = b[i];
	}
	for (i = half; i<size; i++){
		a[i] = 0;
	}
	itr = log2(size);
	while (itr>0){
		i = 0;
		j = 0;

		while (i < size){
			b[j] = (a[i] + a[i + 1]) / twosqrt;
			b[half + j] = (a[i] - a[i + 1]) / twosqrt;
			i += 2;
			j++;
		}
		for (i = 0; i < size; i++){
			a[i] = b[i];
		}
		itr--;
	}
	writeintofile(b, size,"datahaar.txt");
}