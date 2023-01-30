#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
argv[2]: STIT option (0: Iso, 1: Aniso, 2: AnisoDisturbed, 3: AnisoEllip)

All STIT:
	argv[3]: timeStop
	argv[4]: timeGaussStop
	argv[5]: sigma
	argv[6]: window side length
	argv[7]: omega

STIT Anisotropic / Anisotropic Disturbed:
	argv[8]: nDir (number of directions)
	argv[9 + 2*i]: angle dir[i]
	argv[9 + 2*i + 1]: probability dir[i]
	for each i = 0,...,nDir - 1
	
Anisotropic Disturbed:
	argv[9 + 2*nDir]: ellipse half axis b

Anisotropic Elliptic:
	argv[8]: ellipse half axis b
*/
int main(int argc, char *argv[]){
	
	unsigned long lifeOption = atoi(argv[1]);
	unsigned long option = atoi(argv[2]);
	double timeStop = atof(argv[3]);
	double timeGaussStop = atof(argv[4]);
	double sigma = atof(argv[5]);
	unsigned long a = atoi(argv[6]);
	double omg = atof(argv[7]);
	
	CDT2d cdt;
	
	EmptyCDT2d(&cdt,a);
	InitCDT2d(&cdt);

	// STIT Isotropic
	if (option == 0){
		STIT2dIso(&cdt, timeStop, lifeOption, omg);
		//Gauss modification
		NoBoundary(&cdt);
		STIT2dGauss(&cdt, timeGaussStop, sigma, lifeOption, omg);
	}	
	// STIT Anisotropic
	else if ((option == 1) || (option == 2)){
		unsigned long nDir = atoi(argv[8]);
		double *angleDir, *probDir;
		if (nDir > MAX_NUMBER_OF_DIR){
			fprintf(stderr,"The number of directions must be less or equal than %d.\n",MAX_NUMBER_OF_DIR);
			exit(EXIT_FAILURE);		
		}
		angleDir = Malloc(MAX_NUMBER_OF_DIR, double);
		probDir = Malloc(MAX_NUMBER_OF_DIR, double);
		int iDir = 0;
		while (iDir < nDir){
			angleDir[iDir] = atof(argv[9 + 2*iDir]);
			probDir[iDir] = atof(argv[9 + 2*iDir + 1]);
			iDir++;
		}
		if (option == 1){
			STIT2dAniso(&cdt, timeStop, angleDir, probDir, nDir, lifeOption, omg);
			//Gauss modification
			NoBoundary(&cdt);
			STIT2dGauss(&cdt, timeGaussStop, sigma, lifeOption, omg);
		}
		else{
			double bEllip = atof(argv[9 + 2*nDir]);
			STIT2dAnisoDisturbed(&cdt, timeStop, angleDir, probDir, nDir, bEllip, lifeOption, omg);		
			//Gauss modification
			NoBoundary(&cdt);
			STIT2dGaussDisturbed(&cdt, timeGaussStop, sigma, bEllip, lifeOption, omg);
		}

		Free(angleDir);
		Free(probDir);
	}
	else{
		double bEllip = atof(argv[8]);
		STIT2dAnisoEllip(&cdt, timeStop, bEllip, lifeOption, omg);
		//Gauss modification
		NoBoundary(&cdt);
		STIT2dGaussDisturbed(&cdt, timeGaussStop, sigma, bEllip, lifeOption, omg);
	}
	
	//ImageCDT2d(&cdt, &image, "CDT2d");
	StatSTIT(&cdt);
	PlotCDT2d(&cdt);
	FreeCDT2d(&cdt);
	
	return 0;
}
