#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>

#include "CDT2d.h"

/**
STIT Tessellation is simulated until time t = time_stop, 
restricted to a window W of side length a (integer greater than 0).
**/
/*
argc: amount of arguments.
argv[]: array with the parameters.
argv[0]: ./main
argv[1]: lifetime option (0: perimeter, 1: area)
argv[2]: STIT option (0: ISO, 1: DRECT, 2: ELLIP)

All STIT:
	argv[3]: timeStop
	argv[4]: window side length
	argv[5]: omega
	argv[6]: seed

STIT DRECT / ELLIP:
	argv[7]: nDir (number of directions)
	argv[8 + 2*i]: angle dir[i]
	argv[8 + 2*i + 1]: probability dir[i]
	for each i = 0,...,nDir - 1
		
DRECT:
	argv[8 + 2*nDir]: ellipse half axis b

ELLIP:
	argv[7]: ellipse half axis b

*/
int main(int argc, char *argv[]){
	
	unsigned long lifeOption = atoi(argv[1]);
	unsigned long option = atoi(argv[2]);
	double timeStop = atof(argv[3]);
	unsigned long a = atoi(argv[4]);
	double omg = atof(argv[5]);
	int seed = atoi(argv[6]);
	unsigned int feasible = 1;
	
	CDT2d cdt;
	
	EmptyCDT2d(&cdt,a);
	InitCDT2d(&cdt,seed);

	if ((omg < 0.0) || (omg >= 0.5)){
		fprintf(stderr,"Invalid omega parameter value (0 <= omg < 0.5).\n");
		exit(EXIT_FAILURE);
	}
	
	// STIT ISO
	if (option == 0){
		feasible = STIT2dIso(&cdt, timeStop, lifeOption, omg);
	}	
	// STIT Anisotropic
	else if (option == 1){
		unsigned long nDir = atoi(argv[7]);
		double bEllip = atof(argv[8 + 2*nDir]);
		double *angleDir, *probDir;
		if (nDir > MAX_NUMBER_OF_DIR){
			fprintf(stderr,"The number of directions must be less or equal than %d.\n",MAX_NUMBER_OF_DIR);
			exit(EXIT_FAILURE);		
		}
		angleDir = Malloc(MAX_NUMBER_OF_DIR, double);
		probDir = Malloc(MAX_NUMBER_OF_DIR, double);
		int iDir = 0;
		while (iDir < nDir){
			angleDir[iDir] = atof(argv[8 + 2*iDir]);
			probDir[iDir] = atof(argv[8 + 2*iDir + 1]);
			iDir++;
		}
		feasible = STIT2dAnisoDisturbed(&cdt, timeStop, angleDir, probDir, nDir, bEllip, lifeOption, omg);
		Free(angleDir);
		Free(probDir);
	}
	//STIT Ellip
	else{
		double bEllip = atof(argv[7]);
		feasible = STIT2dAnisoEllip(&cdt, timeStop, bEllip, lifeOption, omg);
	}
	if(feasible){
		NoBoundary(&cdt);
		PlotCDT2d(&cdt);
	}
	StatSTIT(&cdt,feasible);
	FreeCDT2d(&cdt,feasible);
	
	return 0;
}
