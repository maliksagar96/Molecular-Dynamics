#include <iostream>
#include "molcDyn.h"
#include <cmath>
#include <time.h>
#include <vector>
#include <fstream>

using namespace std;

int main () {
	
	long int i,j;

	const int snap = 250;
	const int cycles = pow(10,5);
	
	const double vmin = -0.50;
	const double vmax = 0.50;
	
	posInitLJ();
	/*forceInit();
	velocityInit(vmin, vmax, 212);
	
	energyPlot(cycles, snap);
	*/
	return 0;
}

