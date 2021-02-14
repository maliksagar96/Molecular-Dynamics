#include <iostream>
#include <time.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "molcDyn.h"

using namespace std;

long int  particles = 50;
const double bondLength = 1.0, FEND = 0.0;		

const double sigma = 0.8 * bondLength, epsilon = 1.0, mass = 1.0;		//Effectively Epsilon is 1.
const double rc = sigma * pow(2.0,0.16666666);
double dt = 0.005;
double boxLength = 40.10 * bondLength;

//FEND = in +ve x-direction for particle at (0,0,0) and in -ve X direction for last particle in chain

double dt2by2 = (dt*dt)/(2.0*mass);
double frc = epsilon * (12.0*(pow(sigma, 12)/pow(rc,13))-(6.0*(pow(sigma,6)/pow(rc,7))));
double sig6 = pow(sigma,6)/pow(rc,6);
//double urc = epsilon*sig6*(sig6-1);
double urc = epsilon*((pow(sigma,12.0)/pow(rc,12.0)) - (pow(sigma,6.0)/pow(rc,6.0)));
//double urc = epsilon*((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/(rc*rc*rc*rc*rc*rc*rc*rc*rc*rc*rc*rc)) - (1/(rc*rc*rc*rc*rc*rc)));
double stiffness = 100.0;

double potentialEnergy = 0.0, kineticEnergy = 0.0, totalEnergy = 0.0;

double x[1000];double y[1000];double z[1000];
double vx[1000];double vy[1000];double vz[1000];
double fx[1000];double fy[1000];double fz[1000];
double newForceX[1000];double newForceY[1000];double newForceZ[1000];

double randInRange(double min, double max) {
	return min + (max - min)*(rand()/(double)(RAND_MAX));;
}

void posInitLJ() {
	int i, maxCount, countx = 1, county = 1, countz = 1;
	double x_init = 1.0, y_init = 1, z_init = 1.0; 
	double sepFactor = sigma * 1.2;
	
	maxCount = (int)(boxLength/sepFactor);
	
	x[0] = x_init;
	y[0] = y_init;
	z[0] = z_init;


	
	for(i = 1;i<particles;i++) {
		x[i+1] = 0;
		y[i+1] = 0;
		z[i+1] = z_init;
	}
	
	for(i = 1;i<particles;i++) {		
		
		if(countx == maxCount) {
			x[i] = x_init;	
			countx = 1;
		}
		
		else {
			x[i] = x[i-1] + sepFactor;
			countx++;
		}
	}
	
	countx=1;
	
	for(i=1;i<particles;i++) {
		
		
		
		if(countx == maxCount) {
			y[i] = y[i-1] + sepFactor;
			countx = 1;
			county++;
		}
	
	
		else {
			y[i] = y[i-1];
			countx++;
		}
		

		if(county > maxCount){
			y[i] = y_init;
			countx = 1;
			county = 1;
		}		
		
	}
	
	for(i=1;i<particles;i++) {
		
		if(countz == (maxCount*maxCount)) {
			z[i] = z[i-1] + sepFactor;
			countz = 1;
		}
		
		else{
			z[i] = z[i-1];
			countz++;
		}
	}

	structureFile(x,y,z);	
	
}

void velocityInit(double vmin, double vmax, int seed) {
	srand(seed);
	int i;
	double avgVx = 0.0, avgVy = 0.0, avgVz = 0.0;
	vx[0] = 0;
	vy[0] = 0;
	vz[0] = 0;
	
//Initializing velocity compoenents from vmin to vmax.	
	for(i = 1;i<particles;i++) {
		vx[i] = randInRange(vmin, vmax);
		vy[i] = randInRange(vmin, vmax);
		vz[i] = randInRange(vmin, vmax);
	}
	
	for(i=0;i<particles;i++) {
			avgVx = avgVx + vx[i];
			avgVy = avgVy + vy[i];
			avgVz = avgVz + vz[i];
	}
	
	avgVx = avgVx/particles;
	avgVy = avgVy/particles;
	avgVz = avgVz/particles;
	
	for(i = 0;i<particles;i++) {
		vx[i] = vx[i] - avgVx;
		vy[i] = vy[i] - avgVy;
		vz[i] = vz[i] - avgVz;
	}
}


void forceInit() {
	int i;
	
//Initializing Forces to 0.
	for(i = 0;i<particles;i++) {
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;
		
	}
	
	for(i=0;i<particles;i++) {
		newForceX[i] = 0;
		newForceY[i] = 0;
		newForceZ[i] = 0;
	}		
}


void posLine() {
	int i;
	boxLength = 2*particles*bondLength;
	double x_init = boxLength/4, y_init = boxLength/2, z_init = boxLength/2; 

	double l = bondLength;
	
	x[0] = x_init;
	y[0] = y_init;
	z[0] = z_init;
	
	for(i = 1;i<particles;i++) {
		
		x[i] = x[i-1] + l;
		y[i] = y_init;
		z[i] = z_init;
	}
	
	structureFile(x,y,z);
}



