#ifndef MOLCDYN_H
#define MOLCDYN_H

#include <vector>

using namespace std;

/***INIT.CPP***/
/***Defining Parameter which won't change throughout the program***/
extern long int particles;
extern double boxLength;
extern const double sigma, epsilon, mass,bondLength ,relaxLen, FEND, rc;
extern double dt;
extern double dt2by2, stiffness;
extern double frc, urc;
extern double potentialEnergy, kineticEnergy, totalEnergy;

/***Defining all the dynamic quantites***/
extern double x[1000], y[1000], z[1000];
extern double vx[1000], vy[1000], vz[1000];
extern double fx[1000], fy[1000], fz[1000];
extern double newForceX[1000], newForceY[1000], newForceZ[1000];

/***Defining functions to initialise position,force and velocity***/
/**********************************************************************************
	**posInit2D and velocityInit2D are for 2D initialization and initates all z = 0.
	**Test contains code for testing of 2 particles only. Initialises pos, velocity and force all equal to 0.
	  Position of 1st particle is 1.5 units away.
*/

void posInitLJ();
void posLine();
void velocityInit(double, double, int);
void forceInit();
void velocityInit2D(double, double, int);
void test();

double randInRange(double, double );
void verlet(int, int);

/***MOLCDYN.CPP***/
/***Defining updating functions***/
void updatePos();
void updateVelocity();
void updateForce();
void updateForceLJ();
void velocityVerlet(int, int);
void verlet(int, int);
void energyPlot(int, int);
void dataCollection(int, int);
void verlet(int, int);
void energyEquil1(int, int);
void energyEquil2(int, int);

/*--------------------------------------------------------------------------------------------*/

/***CALC.CPP***/
/***Defining functions to calculate stuff***/ 
double calcAvg(vector<double>&);
double calcRg(vector<double>&, vector<double>&,vector<double>&, double, double, double);
double endDistance(double, double, double, double, double, double);
double springForce(double, double, double);
double ljForce(double, double, double);
double minDistance(double, double, double);			//For periodic Boundary Conditions
double ppp(double, double);
double calcKE(double vx_[], double vy_[], double vz_[]);
double ljPE(double);
double springPE(double);
/*--------------------------------------------------------------------------------------------*/

/***PUBLISH.CPP***/
void writeData();
void trajectoryFile();
void structureFile(double x_[],double y_[], double z_[]);
void trajectoryFile(vector<double>& ,vector<double>& , vector<double>& ,
					vector<double>& ,vector<double>& , vector<double>& ,int );			//To print out any vector quantity at a time
void printVector(vector<double>&);
void printVelocity(double vx_[]);
void equilVector(vector<double>&, vector<double>&);
void printConstants();
void finalData();



#endif