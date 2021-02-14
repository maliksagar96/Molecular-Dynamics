#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "molcDyn.h"

using namespace std;

/*** To write initial configuration into a structure.dat file ***/
void structureFile(double x_[],double y_[], double z_[]) {
	fstream positions;
	positions.open("structure.dat", ios::out);
	positions<<"#Initializing at default lattice positions\n\n";
	positions<<particles<<" atoms\n\n";
	positions<<"1 atom types\n\n";
	positions<<0.0<<"\t"<<setprecision(10)<<boxLength<<" xlo xhi\n";
	positions<<0.0<<"\t"<<boxLength<<" ylo yhi\n";
	positions<<0.0<<"\t"<<boxLength<<" zlo zhi\n\n";
	positions<<"Atoms\n\n";

	for(int i = 0;i< particles;i++) {
		positions<<setw(3)<<i+1<<setw(3)<<1<<setw(21)<<x_[i]<<setw(21)<<y_[i]<<setw(21)<<z_[i]<<endl;
	}	
}

void writeData(){
	int a = 10, b = 15, c = 120;
	ofstream dataFile;
	dataFile.open ("example.dat");
	dataFile << "Writing this to a file.\n";
	dataFile << a <<" "<<b<<" "<<c;
	dataFile.close(); 
}

void printConstants() {
	cout<<"Sigma = "<<sigma<<endl;
	cout<<"Rc = "<<rc<<endl;
	cout<<"Urc = "<<urc<<endl;
	cout<<"Frc = "<<frc<<endl;
}

void equilVector(vector<double>& x1, vector<double>& x2) {
	
	for(int i = 0;i<x1.size();i++) {
		x2[i] = x1[i];
	}
}

void printVelocity(double vx_[]) {
	int i;
	double avg = 0;
	for(i=0;i<particles;i++) {
		avg = avg + vx_[i];
	}
	cout<<"Net velocity = "<<avg<<endl;
}
