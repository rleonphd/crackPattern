#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include <limits.h>

#include "CDT2d.h"

/**
	@function: EmptyCDT2d
	---------------------
**/
int EmptyCDT2d(CDT2d *cdt, unsigned long size){

	cdt->size = size;
	
	cdt->nSize = Malloc(2, unsigned long);
	cdt->nSize[0] = 1024; cdt->nSize[1] = 1024;
	
	cdt->Polygon = Malloc(MAX_NUMBER_OF_POLYGONS, POLYGON);
	cdt->I_Segment = Malloc(MAX_NUMBER_OF_I_SEGMENTS, I_SEGMENT);
	
	cdt->numberOfPolygons = 0;
	cdt->numberOfI_Segments = 0;
	
	cdt->Phi = Calloc(NUMBER_OF_INT, double);
	cdt->PhiCum = Malloc(NUMBER_OF_INT, double);
	
	return 0;
}

/** 
	@function: InitCDT2d
	--------------------
**/
void InitCDT2d(CDT2d *cdt, int seed){
	unsigned long is;
	double winSize = (double) (cdt->size);
	
	// random seed is initialized.
 	//srandom(time(NULL));
 	srandom(seed);
 	
	cdt->Polygon[0].numberOfVertices = 4;
	cdt->Polygon[0].V = (double **)calloc2d(4, 2, sizeof(double));

	cdt->Polygon[0].V[0][0] = 0.0; cdt->Polygon[0].V[0][1] = 0.0;
	cdt->Polygon[0].V[1][0] = winSize; cdt->Polygon[0].V[1][1] = 0.0;
	cdt->Polygon[0].V[2][0] = winSize; cdt->Polygon[0].V[2][1] = winSize;
	cdt->Polygon[0].V[3][0] = 0.0; cdt->Polygon[0].V[3][1] = winSize;
	
	WidthFunction(cdt, cdt->Polygon);
	PerimeterPolygon(cdt->Polygon);
	AreaPolygon(cdt->Polygon);
	RoundnessPolygon(cdt->Polygon);
	
	(cdt->numberOfPolygons)++;
	
	// window edge (0,0) --> (size,0)
	cdt->I_Segment[0].x0 = 0.0, cdt->I_Segment[0].y0 = 0.0; 
	cdt->I_Segment[0].x1 = winSize, cdt->I_Segment[0].y1 = 0.0;
	
	// window edge (size,0) --> (size,size)
	cdt->I_Segment[1].x0 = winSize, cdt->I_Segment[1].y0 = 0.0; 
	cdt->I_Segment[1].x1 = winSize, cdt->I_Segment[1].y1 = winSize;
	
	// window edge (size,size) --> (0, size)
	cdt->I_Segment[2].x0 = winSize, cdt->I_Segment[2].y0 = winSize; 
	cdt->I_Segment[2].x1 = 0.0; cdt->I_Segment[2].y1 = winSize;
	
	// window edge (0, size) --> (0,0)
	cdt->I_Segment[3].x0 = 0.0, cdt->I_Segment[3].y0 = winSize; 
	cdt->I_Segment[3].x1 = 0.0; cdt->I_Segment[3].y1 = 0.0;
	
	cdt->numberOfI_Segments += 4;

	for(is = 0; is < 4; is++){
		cdt->I_Segment[is].beta = 0.0;
		cdt->I_Segment[is].length = winSize;
	} 	
}

/**
	@function: FreeCDT2d
	--------------------
**/
int FreeCDT2d(CDT2d *cdt, unsigned int feasible){
	unsigned long i;
	
	/*Free Polygons*/
	for(i = 0; i < cdt->numberOfPolygons; i++){
		Free2d(((cdt->Polygon) + i)->V);
		Free2d(((cdt->Polygon) + i)->H);
		Free(((cdt->Polygon) + i)->B);
	}
	Free(cdt->Polygon);
	
	if(feasible){
		/*Free No Boundary Polygons*/
		for(i = 0; i < cdt->numberOfNBPolygons; i++){
			Free2d(((cdt->NoBoundaryPolygon) + i)->V);
			Free2d(((cdt->NoBoundaryPolygon) + i)->H);
			Free(((cdt->NoBoundaryPolygon) + i)->B);
		}
		Free(cdt->NoBoundaryPolygon);
	}
	
	/*Free CDT2d*/
	Free(cdt->I_Segment);
	Free(cdt->nSize);
	Free(cdt->Phi);
	Free(cdt->PhiCum);
	
	return 0;
}

/**
	@function SetPhiIso
	-------------------
	Sets the phi distribution function and phi cumulative 
	distribution function for isotropic STIT.
**/
int SetPhiIso(CDT2d *cdt){
	unsigned long j;
	
	for(j = 0; j < NUMBER_OF_INT; j++){
		cdt->Phi[j] = 1.0 / NUMBER_OF_INT;
		cdt->PhiCum[j] = (double)(j + 1) / NUMBER_OF_INT;		
	}
	
	return 0;
}

/**
	@function SetPhiAniso
	-------------------
	Sets the phi distribution function and phi cumulative 
	distribution function for anisotropic STIT.
**/
int SetPhiAniso(CDT2d *cdt, double *angleDir, double *probDir, unsigned long nDir){
	unsigned long j, iDir;
	double sumProb = 0.;
	
	for(iDir = 0; iDir < nDir; iDir++){
		// checking if the angle is negative or greater than 1.
		if  ((angleDir[iDir] < 0.0) || (angleDir[iDir] >= 1.0)){
			fprintf(stderr,"The angle of the normal direction must be in the interval [0,pi[\n");
			exit(EXIT_FAILURE);
		}
		// checking if the probability is less or equal than 0 or greater than 1.
		if  ((probDir[iDir] <= 0.0) || (probDir[iDir] >= 1.0)){
			fprintf(stderr,"The probability of the normal direction must be in the interval ]0,1[\n");
			exit(EXIT_FAILURE);
		}
	}
	#define EPS 1.0E-16
	// checking if the sum of the probabilities is 1.
	for(iDir = 0; iDir < nDir; iDir++)
		sumProb += probDir[iDir];
	if (fabs(1 - sumProb) > EPS){
		fprintf(stderr,"The sum of probabilities is not 1.\n");
		exit(EXIT_FAILURE);
	}
	#undef EPS
	
	//computing phi discrete
	for(iDir = 0; iDir < nDir; iDir++){
		// the index j of the normal direction is computed ==> j = floor(angle*NUM_OF_INT) 
		j = floor(angleDir[iDir]*NUMBER_OF_INT);
		//the probability in the position j is updated.
		cdt->Phi[j] += probDir[iDir];
	}
	/* 
	the cummulative distribution function is computed, where 
	cumm_phi[0] = phi[0], and
	cumm_phi[j] = cumm_phi[j - 1] + phi[j] for all j = 1,...,N_OF_INT 
	*/
	cdt->PhiCum[0] = cdt->Phi[0];
	for(j = 1; j < NUMBER_OF_INT; j++){
		cdt->PhiCum[j] = cdt->PhiCum[j - 1] + cdt->Phi[j];
	}
	
	return 0;
}

/**
@function SetPhiAnisoEllip

Sets the phi distribution function and phi cumulative 
distribution function for anisotropic elliptic STIT.
**/
int SetPhiAnisoEllip(CDT2d *cdt, double bEllip){
	unsigned long j;
	int N2;
	double angle, area;
	
	N2 = NUMBER_OF_INT / 2;
	for(j = 0; j < N2 - 1; j++){
		// angle for the index j.
		angle = (j + 1) * M_PI / NUMBER_OF_INT;
		/* computing the area of the ellipse sector [0,angle_j] 
		   divided by the half area of the ellipse (cummulative distribution function).*/
		area = (1.0 / M_PI) * atan(tan(angle) / bEllip);
		/* assign to position j of the discretisation.*/
		cdt->PhiCum[j] = area;
	}
	
	cdt->PhiCum[N2 - 1] = 0.5;
	for(j = 0; j < N2 - 1; j++){
		/* assign to position j of the discretisation. */
		cdt->PhiCum[N2 + j] = 1.0 - cdt->PhiCum[N2 - j - 2];
	}
	cdt->PhiCum[NUMBER_OF_INT - 1] = 1.0;
	
	cdt->Phi[0] = cdt->PhiCum[0];
	for(j = 1; j < NUMBER_OF_INT; j++){
		cdt->Phi[j] = cdt->PhiCum[j] - cdt->PhiCum[j - 1];
	}
	
	return 0;
}

/**
	@function: InitTime
	-------------------
**/
void InitTime(CDT2d *cdt, unsigned long lOption){
	unsigned long j;
	double bPhi = 0.0;
	
	// The lifetime of the polygon is determined. The lifetime
	// is exponentially distribuited with parameter bPhi, which is
	// the weighted sum of the width and the weights are the probabilities
	// Phi[j].
	
	// Perimeter
	if(lOption == 0){
		for(j = 0; j < NUMBER_OF_INT; j++){
			bPhi += cdt->Polygon[0].B[j]*cdt->Phi[j];
		}
		cdt->Polygon[0].tau = RandomExponential(bPhi);
	}
	// Area
	else{
		cdt->Polygon[0].tau = RandomExponential(cdt->Polygon[0].area);
	}
	cdt->Polygon[0].beta = 0.0;
	cdt->Polygon[0].time = cdt->Polygon[0].beta + cdt->Polygon[0].tau;
}

/**
	@function RandomDirection
	-------------------------
	Computes a random direction with the distribution phi_discrete.
	
	[IN]
	@param cdt: CDT2d.
	@param randNum: random number uniformly distributed.
	
	[OUT]
	@param j: index of the normal direction.
**/
int RandomDirection(CDT2d *cdt, double randNum){
	unsigned long j;
	
	for(j = 0; j < NUMBER_OF_INT; j++){
		if (randNum < cdt->PhiCum[j])
			break;
	}
	return j;
}

/**
	@function RandomIntersectionPolygon
	-----------------------------------
	Computes the division line that improves the roundness of the new polygons.
	
	[IN]
	@param cdt: CDT2d.
	@param i: polygon number which is divided.
	@param ip: polygon number of the new cell.
	@param lOption: 0 (LSTIT lifetime) - 1 (LAREA lifetime).
	@param omg: angle for checking the avoid small angle (ASA) condition.
**/
#define EPS 1.0E-12
unsigned int RandomIntersectionPolygon(CDT2d *cdt, unsigned long i, unsigned long ip, unsigned long lOption, double omg){
	unsigned long j, k, k1, k2;
	double A[2], B[2], C[2], magA, magB, magC, cosAB, cosCD, omg1, omg2, bPhi;
	double deltaX, deltaY, det, t;
	
	/*variables for computing the mesh*/
	double randNum, alpha, cosAlpha, sinAlpha, H0, H1, deltaH, d;
	double a1, a2, p1, p2, r1, r2, rd, maxRD;
	unsigned long ind, nMesh = 25, s, maxSize = 50, maxRej = 100, ir = 0;
	
	POLYGON P;
	/*Polygons for computing the mesh*/
	POLYGON P1,P2;
	/*Line which maximize the minimum roundness*/
	double bestLine[2][2] = {{0,0},{0,0}};
	/*Best Polygons*/
	POLYGON MP1,MP2;
	/* k0[1] - k0[0] is the number of vertices between the first and the second intersection point. */
	/*auxiliar variable for computing vertices of each line*/
	double line[2][2] = {{0,0},{0,0}};
	/*auxiliar variable for the intersection points*/
	unsigned long pint[2];

	/* auxiliary variable for checking if a new random line must be generated.*/
	int DETcheck = 0;
	
	// copy of the i-th Polygon into polygon P.
	P.numberOfVertices = cdt->Polygon[i].numberOfVertices;
	P.V = (double **)calloc2d(P.numberOfVertices, 2, sizeof(double));
	for(k = 0; k < P.numberOfVertices; k++){
		P.V[k][0] = cdt->Polygon[i].V[k][0];
		P.V[k][1] = cdt->Polygon[i].V[k][1];
	}
	P.time = cdt->Polygon[i].time;
	
	/*initialize best polygons*/
	MP1.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	MP2.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	MP1.numberOfVertices = 0;
	MP2.numberOfVertices = 0;
	
	/*initialize auxiliar polygons*/
	P1.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	P2.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	
	omg = M_PI*omg;
			
	do{
		/*
		procedure for obtaining a random line that maximizes
		the minimum roundness between the new cells.
		*/
		
		/*random direction*/
		randNum = RandomUniform();
		j = RandomDirection(cdt, randNum);

		/*direction angle*/
		alpha = M_PI * (double)j / (double)NUMBER_OF_INT;
		cosAlpha = cos(alpha);
		sinAlpha = sin(alpha);
		
		/*computing deltaH of the polygon in the direction j*/
		H0 = cdt->Polygon[i].H[j][0];
		H1 = cdt->Polygon[i].H[j][1];
		/*printf("[H0,H1] = [%.9f,%.9f]\n",H0,H1);*/
		deltaH = H1 - H0;
		
		maxRD = 0;
		for(s = 1; s < nMesh; s++){
			/*computing d for each interval in the mesh*/
			d = H0 + (deltaH)*((double)s / (double)nMesh);
			/*computing intersection points*/
			ind = 0;
			for(k = 0; k < P.numberOfVertices; k++){
				k1 = (k + 1) % P.numberOfVertices;
				deltaX = P.V[k1][0] - P.V[k][0];
				deltaY = P.V[k1][1] - P.V[k][1];
				det = deltaX * cosAlpha + deltaY * sinAlpha;
				if(fabs(det) > EPS){
					// the segment between vertices k and (k + 1) is parametrized 
					// with parameter t in [0,1]. The parameter value of the 
					// intersection point is computed.
					t = (d - P.V[k][0] * cosAlpha - P.V[k][1] * sinAlpha) / det;
					if(t > EPS && t < (1.0 - EPS)){
						//printf("IND: %ld - t: %.20f\n",ind,t);
						line[ind][0] = P.V[k][0] + t * deltaX;
						line[ind][1] = P.V[k][1] + t * deltaY;
						pint[ind] = k;
						ind++;
					}
				}
			}
			if(ind < 2){
				//No intersection!
				//printf("No intersection\n");
				continue;
			}
			
			P1.V[0][0] = line[1][0];
			P1.V[0][1] = line[1][1];
			P1.V[1][0] = line[0][0];
			P1.V[1][1] = line[0][1];
			P1.numberOfVertices = 2;
			
			P2.V[0][0] = line[0][0];
			P2.V[0][1] = line[0][1];
			P2.V[1][0] = line[1][0];
			P2.V[1][1] = line[1][1];
			P2.numberOfVertices = 2;

			for(k = 0; k < P.numberOfVertices; k++) {
				k1 = (k + pint[0] + 1) % P.numberOfVertices;
				k2 = (k + pint[1] + 1) % P.numberOfVertices;
				if(P1.numberOfVertices < (pint[1] - pint[0] + 2)){
					P1.V[P1.numberOfVertices][0] = P.V[k1][0];
					P1.V[P1.numberOfVertices][1] = P.V[k1][1];
					P1.numberOfVertices++;
				}
				if(P2.numberOfVertices < (P.numberOfVertices - (pint[1] - pint[0]) + 2)){
					P2.V[P2.numberOfVertices][0] = P.V[k2][0];
					P2.V[P2.numberOfVertices][1] = P.V[k2][1];
					P2.numberOfVertices++; 
				}
			}
			
			/*computing area and perimeter of Polygons P1 and P2*/
			p1 = PerimeterV(P1.V,P1.numberOfVertices);
			a1 = AreaV(P1.V,P1.numberOfVertices);
			r1 = 4*M_PI*a1 / (p1*p1);
			
			p2 = PerimeterV(P2.V,P2.numberOfVertices);
			a2 = AreaV(P2.V,P2.numberOfVertices);
			r2 = 4*M_PI*a2 / (p2*p2);
			
			/*computing sum of squares of roundness between Polygons P1 and P2*/
			rd = r1*r1 + r2*r2;
			if (rd > maxRD){
				maxRD = rd;
				/*current best line*/
				bestLine[0][0] = line[0][0];
				bestLine[0][1] = line[0][1];
				bestLine[1][0] = line[1][0];
				bestLine[1][1] = line[1][1];
				
				/*current best polygons*/	
				for(k = 0; k < P1.numberOfVertices; k++){
					MP1.V[k][0] = P1.V[k][0];
					MP1.V[k][1] = P1.V[k][1];
				}
				MP1.numberOfVertices = P1.numberOfVertices;
				for(k = 0; k < P2.numberOfVertices; k++){
					MP2.V[k][0] = P2.V[k][0];
					MP2.V[k][1] = P2.V[k][1];
				}
				MP2.numberOfVertices = P2.numberOfVertices;
			}
		}
		DETcheck = 1;
		
		// checking small angle 
		A[0] = P.V[pint[0]][0] - bestLine[0][0];
		A[1] = P.V[pint[0]][1] - bestLine[0][1];
		B[0] = bestLine[1][0] - bestLine[0][0];
		B[1] = bestLine[1][1] - bestLine[0][1];
		C[0] = P.V[pint[1]][0] - bestLine[1][0];
		C[1] = P.V[pint[1]][1] - bestLine[1][1];
		magA = sqrt(A[0]*A[0] + A[1]*A[1]);
		magB = sqrt(B[0]*B[0] + B[1]*B[1]);
		magC = sqrt(C[0]*C[0] + C[1]*C[1]);
		cosAB = (A[0]*B[0] + A[1]*B[1]) / (magA*magB);
		cosCD = -(B[0]*C[0] + B[1]*C[1]) / (magB*magC);
		// computing angle between A and B vector
		omg1 = acos(cosAB);
		// computing angle between C and D vector
		omg2 = acos(cosCD);
		// checking if omg1 < omg and PI - omg1 < omg
		if ((omg1 < omg) || (M_PI - omg1) < omg){
			DETcheck = 0;
		}
		// checking if omg2 < omg and PI - omg2 < omg
		if ((omg2 < omg) || (M_PI - omg2) < omg){
			DETcheck = 0;
		}
		
		if (DETcheck == 1) 
			break;
		else
			ir++;
		
	}while(ir < maxRej);
	
	if (DETcheck == 0){
		Free2d(P.V);
		Free2d(P1.V);
		Free2d(P2.V);
		Free2d(MP1.V);
		Free2d(MP2.V);
		return 0;
	}
	else{
		Free2d(cdt->Polygon[i].V);
		Free2d(cdt->Polygon[i].H);
		Free(cdt->Polygon[i].B);

		//printf("New Polygons: %ld -- %ld\n",i,ip);
		// vertices of the new two polygons are determined.
		cdt->Polygon[i].V = (double **)calloc2d(MP1.numberOfVertices, 2, sizeof(double));
		cdt->Polygon[i].numberOfVertices = MP1.numberOfVertices;
		for(k = 0; k < MP1.numberOfVertices; k++){
			cdt->Polygon[i].V[k][0] = MP1.V[k][0];
			cdt->Polygon[i].V[k][1] = MP1.V[k][1];
		}

		cdt->Polygon[ip].V = (double **)calloc2d(MP2.numberOfVertices, 2, sizeof(double));
		cdt->Polygon[ip].numberOfVertices = MP2.numberOfVertices;
		for(k = 0; k < MP2.numberOfVertices; k++){
			cdt->Polygon[ip].V[k][0] = MP2.V[k][0];
			cdt->Polygon[ip].V[k][1] = MP2.V[k][1];
		}
		
		// The width function of Polygon i is computed.
		WidthFunction(cdt, (cdt->Polygon) + i);
		PerimeterPolygon((cdt->Polygon) + i);
		AreaPolygon((cdt->Polygon) + i);
		RoundnessPolygon((cdt->Polygon) + i);
		  
		// The width function of Polygon ip is computed.
		WidthFunction(cdt, (cdt->Polygon) + ip);
		PerimeterPolygon((cdt->Polygon) + ip);
		AreaPolygon((cdt->Polygon) + ip);
		RoundnessPolygon((cdt->Polygon) + ip);
	
		// The lifetime of the polygon is determined. The lifetime
		// is exponentially distribuited with parameter bPhi, which is
		// the weighted sum of the width and the weights are the probabilities
		// Phi[j].
	
		// Lifetime of Polygon i.
		// Perimeter
		if(lOption == 0){
			bPhi = 0.;
			for(j = 0; j < NUMBER_OF_INT; j++){
				bPhi += (cdt->Polygon[i].B[j]) * (cdt->Phi[j]);
			}
			cdt->Polygon[i].tau = RandomExponential(bPhi);
		}
		// Area
		else{
			cdt->Polygon[i].tau = RandomExponential(cdt->Polygon[i].area);
		}
		cdt->Polygon[i].beta = P.time;
		cdt->Polygon[i].time = cdt->Polygon[i].beta + cdt->Polygon[i].tau;
	
		// Lifetime of Polygon ip.
		// Perimeter
		if(lOption == 0){
			bPhi = 0.;
			for(j = 0; j < NUMBER_OF_INT; j++){
				bPhi += (cdt->Polygon[ip].B[j]) * (cdt->Phi[j]);
			}
	
			cdt->Polygon[ip].tau = RandomExponential(bPhi);
		}
		// Area
		else{
			cdt->Polygon[ip].tau = RandomExponential(cdt->Polygon[ip].area);
		}
		cdt->Polygon[ip].beta = P.time;
		cdt->Polygon[ip].time = cdt->Polygon[ip].beta + cdt->Polygon[ip].tau;
	
		// The two vertices of the I_Segment are assgined.
		cdt->I_Segment[cdt->numberOfI_Segments].x0 = bestLine[0][0];
		cdt->I_Segment[cdt->numberOfI_Segments].y0 = bestLine[0][1]; 
		cdt->I_Segment[cdt->numberOfI_Segments].x1 = bestLine[1][0]; 
		cdt->I_Segment[cdt->numberOfI_Segments].y1 = bestLine[1][1];
		// The length of the I_Segment is computed.
		deltaX = bestLine[1][0] - bestLine[0][0];
		deltaY = bestLine[1][1] - bestLine[0][1];
		cdt->I_Segment[cdt->numberOfI_Segments].length = sqrt(deltaX*deltaX + deltaY*deltaY);
		// The birthtime of the I_Segment is assigned.
		cdt->I_Segment[cdt->numberOfI_Segments].beta = P.time;
	
		Free2d(P.V);
		Free2d(P1.V);
		Free2d(P2.V);
		Free2d(MP1.V);
		Free2d(MP2.V);
		return 1;
	}
}

/**
	@function RandomIntersectionPolygonDisturbed
	--------------------------------------------
	Computes the division line that improves the roundness of the new polygons.
	The function is for the anisotropic discrete disturbed and anisotropic ellipse model.
	
	[IN]
	@param cdt: CDT2d.
	@param i: polygon number which is divided.
	@param ip: polygon number of the new cell.
	@param lOption: 0 (LSTIT lifetime) - 1 (LAREA lifetime).
	@param omg: angle for checking the avoid small angle (ASA) condition.
**/

unsigned int RandomIntersectionPolygonDisturbed(CDT2d *cdt, unsigned long i, unsigned long ip, double bEllip, unsigned long lOption, double omg){
	unsigned long j, k, k1, k2, jGamma;
	double A[2], B[2], C[2], magA, magB, magC, cosAB, cosCD, omg1, omg2, bPhi, gamma;

	double deltaX, deltaY, det, t;
	
	/*variables for computing the mesh*/
	double randNum, alpha, cosAlpha, sinAlpha, H0, H1, deltaH, d;
	double a1, a2, p1, p2, r1, r2, rd, maxRD;
	unsigned long ind, nMesh = 25, s, maxSize = 50, maxRej = 100, ir = 0;
	
	POLYGON P;
	/*Polygons for computing the mesh*/
	POLYGON P1,P2;
	/*Line which maximize the minimum roundness*/
	double bestLine[2][2] = {{0,0},{0,0}};
	/*Best Polygons*/
	POLYGON MP1,MP2;
	/* k0[1] - k0[0] is the number of vertices between the first and the second intersection point. */
	/*auxiliar variable for computing vertices of each line*/
	double line[2][2] = {{0,0},{0,0}};
	/*auxiliar variable for the intersection points*/
	unsigned long pint[2];

	/* auxiliary variable for checking if a new random line must be generated.*/
	int DETcheck = 0;
	
	// copy of the i-th Polygon into polygon P.
	P.numberOfVertices = cdt->Polygon[i].numberOfVertices;
	P.V = (double **)calloc2d(P.numberOfVertices, 2, sizeof(double));
	for(k = 0; k < P.numberOfVertices; k++){
		P.V[k][0] = cdt->Polygon[i].V[k][0];
		P.V[k][1] = cdt->Polygon[i].V[k][1];
	}
	P.time = cdt->Polygon[i].time;
	
	/*initialize best polygons*/
	MP1.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	MP2.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	MP1.numberOfVertices = 0;
	MP2.numberOfVertices = 0;
	
	/*initialize auxiliar polygons*/
	P1.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	P2.V = (double **)calloc2d(maxSize, 2, sizeof(double));
	
	omg = M_PI*omg;
			
	do{
		/*
		procedure for obtaining a random line that maximizes
		the minimum roundness between the new cells.
		*/
		
		/*random direction*/
		randNum = RandomUniform();
		j = RandomDirection(cdt, randNum);
		gamma = RandomElliptic(bEllip);
		jGamma = floor(gamma * NUMBER_OF_INT / M_PI);
		j = (j + jGamma) % (NUMBER_OF_INT);

		/*direction angle*/
		alpha = M_PI * (double)j / (double)NUMBER_OF_INT;
		cosAlpha = cos(alpha);
		sinAlpha = sin(alpha);
		
		/*computing deltaH of the polygon in the direction j*/
		H0 = cdt->Polygon[i].H[j][0];
		H1 = cdt->Polygon[i].H[j][1];
		/*printf("[H0,H1] = [%.9f,%.9f]\n",H0,H1);*/
		deltaH = H1 - H0;
		
		maxRD = 0;
		for(s = 1; s < nMesh; s++){
			/*computing d for each interval in the mesh*/
			d = H0 + (deltaH)*((double)s / (double)nMesh);
			/*computing intersection points*/
			ind = 0;
			for(k = 0; k < P.numberOfVertices; k++){
				k1 = (k + 1) % P.numberOfVertices;
				deltaX = P.V[k1][0] - P.V[k][0];
				deltaY = P.V[k1][1] - P.V[k][1];
				det = deltaX * cosAlpha + deltaY * sinAlpha;
				if(fabs(det) > EPS){
					// the segment between vertices k and (k + 1) is parametrized 
					// with parameter t in [0,1]. The parameter value of the 
					// intersection point is computed.
					t = (d - P.V[k][0] * cosAlpha - P.V[k][1] * sinAlpha) / det;
					if(t > EPS && t < (1.0 - EPS)){
						//printf("IND: %ld - t: %.20f\n",ind,t);
						line[ind][0] = P.V[k][0] + t * deltaX;
						line[ind][1] = P.V[k][1] + t * deltaY;
						pint[ind] = k;
						ind++;
					}
				}
			}
			if(ind < 2){
				//No intersection!
				//printf("No intersection\n");
				continue;
			}
			
			P1.V[0][0] = line[1][0];
			P1.V[0][1] = line[1][1];
			P1.V[1][0] = line[0][0];
			P1.V[1][1] = line[0][1];
			P1.numberOfVertices = 2;
			
			P2.V[0][0] = line[0][0];
			P2.V[0][1] = line[0][1];
			P2.V[1][0] = line[1][0];
			P2.V[1][1] = line[1][1];
			P2.numberOfVertices = 2;

			for(k = 0; k < P.numberOfVertices; k++){
				//printf("V[%ld] = (%.6f,%.6f)\n",k,P.V[k][0],P.V[k][1]);
				k1 = (k + pint[0] + 1) % P.numberOfVertices;
				k2 = (k + pint[1] + 1) % P.numberOfVertices;
				if(P1.numberOfVertices < (pint[1] - pint[0] + 2)){
					P1.V[P1.numberOfVertices][0] = P.V[k1][0];
					P1.V[P1.numberOfVertices][1] = P.V[k1][1];
					P1.numberOfVertices++;
				}
				if(P2.numberOfVertices < (P.numberOfVertices - (pint[1] - pint[0]) + 2)){
					P2.V[P2.numberOfVertices][0] = P.V[k2][0];
					P2.V[P2.numberOfVertices][1] = P.V[k2][1];
					P2.numberOfVertices++; 
				}
			}
			
			/*computing area and perimeter of Polygons P1 and P2*/
			p1 = PerimeterV(P1.V,P1.numberOfVertices);
			a1 = AreaV(P1.V,P1.numberOfVertices);
			r1 = 4*M_PI*a1 / (p1*p1);
			
			p2 = PerimeterV(P2.V,P2.numberOfVertices);
			a2 = AreaV(P2.V,P2.numberOfVertices);
			r2 = 4*M_PI*a2 / (p2*p2);
			
			/*computing sum of squares of roundness*/
			rd = r1*r1 + r2*r2;
			if (rd > maxRD){
				maxRD = rd;
				/*current best line*/
				bestLine[0][0] = line[0][0];
				bestLine[0][1] = line[0][1];
				bestLine[1][0] = line[1][0];
				bestLine[1][1] = line[1][1];
				
				/*current best polygons*/	
				for(k = 0; k < P1.numberOfVertices; k++){
					MP1.V[k][0] = P1.V[k][0];
					MP1.V[k][1] = P1.V[k][1];
				}
				MP1.numberOfVertices = P1.numberOfVertices;
				for(k = 0; k < P2.numberOfVertices; k++){
					MP2.V[k][0] = P2.V[k][0];
					MP2.V[k][1] = P2.V[k][1];
				}
				MP2.numberOfVertices = P2.numberOfVertices;
			}
		}
		DETcheck = 1;
		
		// checking small angle 
		A[0] = P.V[pint[0]][0] - bestLine[0][0];
		A[1] = P.V[pint[0]][1] - bestLine[0][1];
		B[0] = bestLine[1][0] - bestLine[0][0];
		B[1] = bestLine[1][1] - bestLine[0][1];
		C[0] = P.V[pint[1]][0] - bestLine[1][0];
		C[1] = P.V[pint[1]][1] - bestLine[1][1];
		magA = sqrt(A[0]*A[0] + A[1]*A[1]);
		magB = sqrt(B[0]*B[0] + B[1]*B[1]);
		magC = sqrt(C[0]*C[0] + C[1]*C[1]);
		cosAB = (A[0]*B[0] + A[1]*B[1]) / (magA*magB);
		cosCD = -(B[0]*C[0] + B[1]*C[1]) / (magB*magC);
		// computing angle between A and B vector
		omg1 = acos(cosAB);
		// computing angle between C and D vector
		omg2 = acos(cosCD);
		//printf("OMG1: %.6f OMG2: %.6f PI-OMG1: %.6f PI-OMG2: %.6f OMG: %.6f\n",omg1,omg2,(M_PI - omg1),(M_PI - omg2),omg);
		// checking if omg1 < omg and PI - omg1 < omg
		if ((omg1 < omg) || (M_PI - omg1) < omg){
			DETcheck = 0;
		}
		// checking if omg2 < omg and PI - omg2 < omg
		if ((omg2 < omg) || (M_PI - omg2) < omg){
			DETcheck = 0;
		}
		
		if (DETcheck == 1) 
			break;
		else
			ir++;
		
	}while(ir < maxRej);

	if (DETcheck == 0){
		Free2d(P.V);
		Free2d(P1.V);
		Free2d(P2.V);
		Free2d(MP1.V);
		Free2d(MP2.V);
		return 0;
	}
	else{
		Free2d(cdt->Polygon[i].V);
		Free2d(cdt->Polygon[i].H);
		Free(cdt->Polygon[i].B);

		//printf("New Polygons: %ld -- %ld\n",i,ip);
		// vertices of the new two polygons are determined.
		cdt->Polygon[i].V = (double **)calloc2d(MP1.numberOfVertices, 2, sizeof(double));
		cdt->Polygon[i].numberOfVertices = MP1.numberOfVertices;
		for(k = 0; k < MP1.numberOfVertices; k++){
			cdt->Polygon[i].V[k][0] = MP1.V[k][0];
			cdt->Polygon[i].V[k][1] = MP1.V[k][1];
		}

		cdt->Polygon[ip].V = (double **)calloc2d(MP2.numberOfVertices, 2, sizeof(double));
		cdt->Polygon[ip].numberOfVertices = MP2.numberOfVertices;
		for(k = 0; k < MP2.numberOfVertices; k++){
			cdt->Polygon[ip].V[k][0] = MP2.V[k][0];
			cdt->Polygon[ip].V[k][1] = MP2.V[k][1];
		}
		
		// The width function of Polygon i is computed.
		WidthFunction(cdt, (cdt->Polygon) + i);
		PerimeterPolygon((cdt->Polygon) + i);
		AreaPolygon((cdt->Polygon) + i);
		RoundnessPolygon((cdt->Polygon) + i);
		  
		// The width function of Polygon ip is computed.
		WidthFunction(cdt, (cdt->Polygon) + ip);
		PerimeterPolygon((cdt->Polygon) + ip);
		AreaPolygon((cdt->Polygon) + ip);
		RoundnessPolygon((cdt->Polygon) + ip);
	
		// The lifetime of the polygon is determined. The lifetime
		// is exponentially distribuited with parameter bPhi, which is
		// the weighted sum of the width and the weights are the probabilities
		// Phi[j].
	
		// Lifetime of Polygon i.
		// Perimeter
		if(lOption == 0){
			bPhi = 0.;
			for(j = 0; j < NUMBER_OF_INT; j++){
				bPhi += (cdt->Polygon[i].B[j]) * (cdt->Phi[j]);
			}
			cdt->Polygon[i].tau = RandomExponential(bPhi);
		}
		// Area
		else{
			cdt->Polygon[i].tau = RandomExponential(cdt->Polygon[i].area);
		}
		cdt->Polygon[i].beta = P.time;
		cdt->Polygon[i].time = cdt->Polygon[i].beta + cdt->Polygon[i].tau;
	
		// Lifetime of Polygon ip.
		// Perimeter
		if(lOption == 0){
			bPhi = 0.;
			for(j = 0; j < NUMBER_OF_INT; j++){
				bPhi += (cdt->Polygon[ip].B[j]) * (cdt->Phi[j]);
			}
	
			cdt->Polygon[ip].tau = RandomExponential(bPhi);
		}
		// Area
		else{
			cdt->Polygon[ip].tau = RandomExponential(cdt->Polygon[ip].area);
		}
		cdt->Polygon[ip].beta = P.time;
		cdt->Polygon[ip].time = cdt->Polygon[ip].beta + cdt->Polygon[ip].tau;
	
		// The two vertices of the I_Segment are assgined.
		cdt->I_Segment[cdt->numberOfI_Segments].x0 = bestLine[0][0];
		cdt->I_Segment[cdt->numberOfI_Segments].y0 = bestLine[0][1]; 
		cdt->I_Segment[cdt->numberOfI_Segments].x1 = bestLine[1][0]; 
		cdt->I_Segment[cdt->numberOfI_Segments].y1 = bestLine[1][1];
		// The length of the I_Segment is computed.
		deltaX = bestLine[1][0] - bestLine[0][0];
		deltaY = bestLine[1][1] - bestLine[0][1];
		cdt->I_Segment[cdt->numberOfI_Segments].length = sqrt(deltaX*deltaX + deltaY*deltaY);
		// The birthtime of the I_Segment is assigned.
		cdt->I_Segment[cdt->numberOfI_Segments].beta = P.time;
	
		Free2d(P.V);
		Free2d(P1.V);
		Free2d(P2.V);
		Free2d(MP1.V);
		Free2d(MP2.V);
		return 1;
	}
}

#undef EPS

/**
	@function STIT2dIso
	-------------------
	Isotropic STIT2d
**/
unsigned int STIT2dIso(CDT2d *cdt, double timeStop, unsigned long lOption, double omg){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	unsigned int feasible = 1;
	
	// phi and cumulative distribution function are computed.
	SetPhiIso(cdt);
	
	//lifetime of the first polygon (window) is computed.
	InitTime(cdt, lOption);
	
	// it assumes that the process continues.
	stitStop = 1;
	do{
		numPol = cdt->numberOfPolygons;
  	q = 0;
  	for(i = 0; i < numPol; i++){
			// checking if the polygon i must be divided.
			if(cdt->Polygon[i].time < timeStop){
				if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
					printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
					return 1;
				}
				if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
					printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
					return 1;
				}
				feasible = RandomIntersectionPolygon(cdt, i, numPol + q, lOption, omg);
				if(feasible == 0){
					q = 0;
					break;
				}
				q++;
				(cdt->numberOfPolygons)++;
				(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return feasible;
}

/**
	@function STIT2dAnisoDisturbed
	-------------------

	Anisotropic Disturbed STIT2d
**/
unsigned int STIT2dAnisoDisturbed(CDT2d *cdt, double timeStop, double *angleDir, double *probDir, unsigned long nDir, double bEllip, unsigned long lOption, double omg){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	unsigned int feasible = 1;
	
	// phi and cumulative distribution function are computed.
	SetPhiAniso(cdt, angleDir, probDir, nDir);
	
	//lifetime of the first polygon (window) is computed.
	InitTime(cdt, lOption);
	
	// it assumes that the process continues.
	stitStop = 1;
	do{
		numPol = cdt->numberOfPolygons;
  	q = 0;
  	for(i = 0; i < numPol; i++){
  		// checking if the polygon i must be divided.
  		if(cdt->Polygon[i].time < timeStop){				
				if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
					printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
					return 1;
				}
				if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
					printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
					return 1;
				}
				feasible = RandomIntersectionPolygonDisturbed(cdt, i, numPol + q, bEllip, lOption, omg);
				if(feasible == 0){
					q = 0;
					break;
				}
				q++;
				(cdt->numberOfPolygons)++;
				(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return feasible;
}

/**
@function STIT2dAnisoEllip

Anisotropic Elliptic STIT2d
**/
unsigned int STIT2dAnisoEllip(CDT2d *cdt, double timeStop, double bEllip, unsigned long lOption, double omg){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	unsigned int feasible = 1;
	
	// phi and cumulative distribution function are computed.
	SetPhiAnisoEllip(cdt, bEllip);
	
	//lifetime of the first polygon (window) is computed.
	InitTime(cdt, lOption);
	
	// it assumes that the process continues.
	stitStop = 1;
	do{
		numPol = cdt->numberOfPolygons;
    	q = 0;
    	for(i = 0; i < numPol; i++){
    		// checking if the polygon i must be divided.
    		if(cdt->Polygon[i].time < timeStop){
					if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
						printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
						return 1;
					}
					if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
						printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
						return 1;
					}
					feasible = RandomIntersectionPolygonDisturbed(cdt, i, numPol + q, bEllip, lOption, omg);
					if(feasible == 0){
						q = 0;
						break;
					}
					q++;
					(cdt->numberOfPolygons)++;
					(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    	if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return feasible;
}

/**	
	@function: WidthFunction 
	-------------------------
	Computes the width function of a polygon.
	
	[IN]
    @param Polygon: polygon for which the width function is computed.
    
    [OUT]
    @param Polygon->H: support function (fixed direction)
    @param Polygon->B: width of each fixed direction.
    @param Polygon->maxWidth: maximal width of the polygon = max(b)
    @param Polygon->minWidth: minimal width of the polygon = min(b)
**/
int WidthFunction(CDT2d *cdt, POLYGON *Polygon){
	unsigned long j, k;
	double cosDir, sinDir, s, sMin, sMax, bMax, bMin, angDir;
  
	Polygon->H = (double **)calloc2d(NUMBER_OF_INT, 2, sizeof(double));
	Polygon->B = (double *)calloc(NUMBER_OF_INT,sizeof(double));
  
	for(j = 0; j < NUMBER_OF_INT; j++){
		angDir = (double)j * M_PI / (double)NUMBER_OF_INT;
		cosDir = cos(angDir);
		sinDir = sin(angDir);
		sMin = DOUBLE_MAX;
		sMax = DOUBLE_MIN;

		for(k = 0; k < Polygon->numberOfVertices; k++){
			s = Polygon->V[k][0] * cosDir + Polygon->V[k][1] * sinDir;
			if(s < sMin) sMin = s;
			if(s > sMax) sMax = s;
		}
		Polygon->H[j][0] = sMin;
		Polygon->H[j][1] = sMax;
	}

	bMax = DOUBLE_MIN;
	bMin = DOUBLE_MAX;
	for(j = 0; j < NUMBER_OF_INT; j++) {
		Polygon->B[j] = Polygon->H[j][1] - Polygon->H[j][0];
		if(Polygon->B[j] > bMax) bMax = Polygon->B[j];
		if(Polygon->B[j] < bMin) bMin = Polygon->B[j];
	}
	Polygon->maxWidth = bMax;
	Polygon->minWidth = bMin;
	Polygon->aspectRatio = bMin / bMax;
	
	Polygon->wML = 1. / ((cdt->size - Polygon->B[0]) * (cdt->size - Polygon->B[NUMBER_OF_INT / 2]));
  
  return 0;
}

/**
	@function: Perimeter
	----------------------
**/
double PerimeterV(double **V, unsigned long nov){
	unsigned long k, k1;
	double deltaX, deltaY, p = 0.;
	
	for(k = 0; k < nov; k++){
		k1 = (k + 1) % nov;
		deltaX = V[k1][0] - V[k][0];
		deltaY = V[k1][1] - V[k][1];
		p += sqrt(deltaX*deltaX + deltaY*deltaY); 
	}
	
	return p;
}

/**
	@function: Area
	----------------------
**/
double AreaV(double **V, unsigned long nov){
	unsigned long k, k1;
	double a = 0.;
	
	for(k = 0; k < nov; k++){
		k1 = (k + 1) % nov;
		a += (V[k][0]*V[k1][1]) - (V[k][1]*V[k1][0]); 
	}
	a = 0.5*a;
	return a;
}

/**
	@function: PerimeterPolygon
	----------------------
**/
int PerimeterPolygon(POLYGON *Polygon){
	unsigned long k, k1;
	double deltaX, deltaY;
	
	Polygon->perimeter = 0.0;
	for(k = 0; k < Polygon->numberOfVertices; k++){
		k1 = (k + 1) % Polygon->numberOfVertices;
		deltaX = Polygon->V[k1][0] - Polygon->V[k][0];
		deltaY = Polygon->V[k1][1] - Polygon->V[k][1];
		Polygon->perimeter += sqrt(deltaX*deltaX + deltaY*deltaY); 
	}
	
	return 0;
}

/**
	@function: AreaPolygon
	----------------------
**/
int AreaPolygon(POLYGON *Polygon){
	unsigned long k, k1;
	
	Polygon->area = 0.0;
	for(k = 0; k < Polygon->numberOfVertices; k++){
		k1 = (k + 1) % Polygon->numberOfVertices;
		Polygon->area += (Polygon->V[k][0]*Polygon->V[k1][1]) - (Polygon->V[k][1]*Polygon->V[k1][0]); 
	}
	Polygon->area = 0.5*(Polygon->area);
	return 0;
}

/**
	@function: RoundnessPolygon
	----------------------
**/
int RoundnessPolygon(POLYGON *Polygon){
	
	Polygon->roundness = 4*M_PI*Polygon->area / (Polygon->perimeter * Polygon->perimeter);
	return 0;
}

/**
	@function RandomExponential
	---------------------------
	Simulates an exponentially distributed random number.
	
	[IN]
	@param lambda: the parameter of the exponential distribution.
	
	[OUT]
	@return a real-valued pseudorandom number.
	
	@warning should be initialized using InitRandomizer() or srandom()  
**/
double RandomExponential(double lambda){
	return -log(RandomUniform()) / lambda;
}

/**
	@function RandomElliptic
	---------------------------
	Generates a random number for the elliptic distribution
	with the cumulative distribution function 
	phi_cum(alpha) = area of ellipse sector [0,alpha] divided by half 
	area of the ellipse with half axes 1 and bEllip (less than 1)
	and alpha in [0,pi[.
	
	[IN]
	@param bEllip: minor half axis (vertical and less than 1).
	[OUT]
	@return random number.
**/
double RandomElliptic(double bEllip){
	double as, c, r, s, tas, x, phi;

	as = asin(bEllip);
	tas = tan(as/2.0);

	do{
		x = (as - bEllip*log(tas)) * RandomUniform();
		r = RandomUniform();
		phi = x;
		if(phi > as) {
			phi = 2.0 * (atan(tas * exp((x-as) / bEllip)));
			r = bEllip * ((1.0/sin(phi) -1.0) * r + 1.0);
		}
		else 
			r = (1.0 - bEllip) * r + bEllip;
	
		c = cos(phi);
		s = sin(phi);
	}while(r <= bEllip / sqrt(bEllip*bEllip*c*c + s*s));

	if(random() % 2)
		phi = M_PI - phi;

	return phi;
}

/**
	@function: StatSTIT
	--------------------
	Computes statistics of STIT simulation.
**/
void StatSTIT(CDT2d *cdt, unsigned int feasible){
	unsigned long i;
	double sumWeights = 0., sumValues, mean, var, std, cv;
	FILE *statFile;
		
	//Writing stats in a file.
	statFile = fopen("statSTIT.txt","w");
	
	if(feasible){
		fprintf(statFile,"%ld\n",cdt->numberOfPolygons);
	
		//Computing sum(wML)
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumWeights += cdt->Polygon[i].wML;
	
		//Area
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * cdt->Polygon[i].area;
		mean = sumValues / sumWeights;
	
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * (cdt->Polygon[i].area - mean) * (cdt->Polygon[i].area - mean);
		
		var = sumValues / sumWeights;
		std = sqrt(var);
		cv = std / mean;
		
		fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);
	
		//Perimeter
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * cdt->Polygon[i].perimeter;
		mean = sumValues / sumWeights;
	
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * (cdt->Polygon[i].perimeter - mean) * (cdt->Polygon[i].perimeter - mean);
		
		var = sumValues / sumWeights;
		std = sqrt(var);
		cv = std / mean;
		
		fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);
	
		//Roundness
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * cdt->Polygon[i].roundness;
		mean = sumValues / sumWeights;
	
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * (cdt->Polygon[i].roundness - mean) * (cdt->Polygon[i].roundness - mean);
		
		var = sumValues / sumWeights;
		std = sqrt(var);
		cv = std / mean;
		
		fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

		//Max width
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * cdt->Polygon[i].maxWidth;
		mean = sumValues / sumWeights;
	
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * (cdt->Polygon[i].maxWidth - mean) * (cdt->Polygon[i].maxWidth - mean);
		
		var = sumValues / sumWeights;
		std = sqrt(var);
		cv = std / mean;
		
		fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

		//Min width
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * cdt->Polygon[i].minWidth;
		mean = sumValues / sumWeights;
	
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * (cdt->Polygon[i].minWidth - mean) * (cdt->Polygon[i].minWidth - mean);
		
		var = sumValues / sumWeights;
		std = sqrt(var);
		cv = std / mean;
		
		fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

		//Aspect ratio
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * cdt->Polygon[i].aspectRatio;
		mean = sumValues / sumWeights;
	
		sumValues = 0.;
		for(i = 0; i < cdt->numberOfPolygons; i++)
			sumValues += cdt->Polygon[i].wML * (cdt->Polygon[i].aspectRatio - mean) * (cdt->Polygon[i].aspectRatio - mean);
		
		var = sumValues / sumWeights;
		std = sqrt(var);
		cv = std / mean;
		
		fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);
	}
	else{
		fprintf(statFile,"%d\n",-1);
		for(i = 0; i < 6; i++)
			fprintf(statFile,"%.12f,%.12f,%.12f\n",0.,0.,0.);
	}

	fclose(statFile);
}

/**
@function NoBoundary
**/
void NoBoundary(CDT2d *cdt){
	unsigned int i, k, state;
	
	cdt->NoBoundaryPolygon = Malloc(MAX_NUMBER_OF_POLYGONS, POLYGON);
	cdt->numberOfNBPolygons = 0;
	
	for(i = 0; i < cdt->numberOfPolygons; i++){
		state = 1;
		for(k = 0; k < cdt->Polygon[i].numberOfVertices; k++){
if ((cdt->Polygon[i].V[k][0] == 0.0 || cdt->Polygon[i].V[k][0] == 1.0) || (cdt->Polygon[i].V[k][1] == 0.0 || cdt->Polygon[i].V[k][1] == 1.0)){
				state = 0;
				break;
			}
		}
		if (state){
			CopyPolygon(cdt->NoBoundaryPolygon + cdt->numberOfNBPolygons, cdt->Polygon + i);
			(cdt->numberOfNBPolygons)++;
		}
	}
}

/**
@function CopyPolygon
**/
void CopyPolygon(POLYGON *NewPolygon, POLYGON *Polygon){
	unsigned int j, k;

	NewPolygon->numberOfVertices = Polygon->numberOfVertices;
	NewPolygon->V = (double **)calloc2d(NewPolygon->numberOfVertices, 2, sizeof(double));
	
	for(k = 0; k < NewPolygon->numberOfVertices; k++){
		NewPolygon->V[k][0] = Polygon->V[k][0];
		NewPolygon->V[k][1] = Polygon->V[k][1];
	}
	
	NewPolygon->H = (double **)calloc2d(NUMBER_OF_INT, 2, sizeof(double));
	for(j = 0; j < NUMBER_OF_INT; j++){
		NewPolygon->H[j][0] = Polygon->H[j][0];
		NewPolygon->H[j][1] = Polygon->H[j][1];
	}
	
	NewPolygon->B = (double *)calloc(NUMBER_OF_INT,sizeof(double));
	for(j = 0; j < NUMBER_OF_INT; j++){
		NewPolygon->B[j] = Polygon->B[j];
	}
	
	NewPolygon->aspectRatio = Polygon->aspectRatio;
	NewPolygon->perimeter = Polygon->perimeter;
	NewPolygon->area = Polygon->area;
	NewPolygon->roundness = Polygon->roundness;
	NewPolygon->beta = Polygon->beta;
	NewPolygon->tau = Polygon->tau;
	NewPolygon->time = Polygon->time;
	NewPolygon->maxWidth = Polygon->maxWidth;
	NewPolygon->minWidth = Polygon->minWidth;
	NewPolygon->wML = Polygon->wML;
}

/**
@function PlotCDT2d
**/
void PlotCDT2d(CDT2d *cdt){
	unsigned long i;
	char buffer[15];
	FILE *file;
	
	sprintf(buffer,"CDT2d.tex");
	file = fopen(buffer, "w");
	
	fprintf(file, "\\documentclass[10pt]{article}\n");
	fprintf(file, "\\usepackage{tikz}\n");
	fprintf(file,"\\usepackage[active,pdftex,tightpage]{preview}\n");
	fprintf(file,"\\PreviewEnvironment[{[]}]{tikzpicture}\n");
	fprintf(file,"\\pagestyle{empty}\n");
	fprintf(file, "\\begin{document}\n");
	fprintf(file, "\\begin\{tikzpicture}[x = 10cm,y = 10cm]\n");
   
	for(i = 0; i < cdt->numberOfI_Segments; i++)
		fprintf(file, "\\draw[thick](%lf,%lf) -- (%lf,%lf);\n", cdt->I_Segment[i].x0, cdt->I_Segment[i].y0, cdt->I_Segment[i].x1, cdt->I_Segment[i].y1);

	fprintf(file, "\\end\{tikzpicture}\n");
	fprintf(file, "\\end\{document}\n");

	fclose(file);
}

/** Allocation of a 2d array. 
    Allocates a continuous area in memory for a 2d array x[0..m-1][0..n-1] 
    of items of size ny. The area is not initialized.   
    @param m [IN] number of rows  
    @param n [IN] number of collumns  
    @param ny [IN] size of an item (element)  
    @return either a null pointer if unable to allocate, or a void pointer
            to the allocated area
***/
void **malloc2d(unsigned long m, unsigned long n, size_t ny) {
  unsigned long i;
  unsigned char *x;
  void* *y;
   
  /* sanity check */
  if( m==0 || n==0 || ny==0) {
    assert( !(m==0 || n==0 || ny==0) );
    return NULL;
  }
  if(m>ULONG_MAX/n) { 
    assert(m<=ULONG_MAX/n); 
    return NULL; 
  }
  if(m*n>ULONG_MAX/ny) { 
    assert(m*n<=ULONG_MAX/ny); 
    return NULL; 
  }

  /* work part */
  if((x=(unsigned char *)malloc(m*n*ny))==NULL) {
    assert(x!=NULL); 
    return NULL; 
  }
  if((y=(void *)malloc(m*sizeof(void *)))==NULL) {
    assert(y!=NULL); 
    free((char *)x); 
    return NULL; 
  }

  for(i=0;i<m;i++) 
    y[i]=&x[i*n*ny];
  
  return (void **)y;
}

/** Allocation and initialization of a 2d array. 
    Allocates a continuous area in memory for a 2d array x[0..m-1][0..n-1] 
    of items of size ny. The area is initialized to all bits zero.   
    @param m [IN] number of rows  
    @param n [IN] number of collumns  
    @param ny [IN] size of an item (element)  
    @return either a null pointer if unable to allocate, or a void pointer to the allocated area
**/
void **calloc2d(unsigned long m, unsigned long n, size_t ny) {
  unsigned long i;
  unsigned char *x;
  void* *y;

  /* sanity check */
  if( m==0 || n==0 || ny==0) {
    assert( !(m==0 || n==0 || ny==0) );
    return NULL;
  }
  if(m>ULONG_MAX/n) { 
    assert(m<=ULONG_MAX/n); 
    return NULL; 
  }
  if(m*n>ULONG_MAX/ny) { 
    assert(m*n<=ULONG_MAX/ny); 
    return NULL; 
  }

  /* work part */
  if((x=(unsigned char *)calloc(m*n,ny))==NULL) {
    assert(x!=NULL); 
    return NULL;
  }
  if((y=(void *)malloc(m*sizeof(void *)))==NULL) {
    assert(y!=NULL); 
    free((char *)x); 
    return NULL; 
  }
  for(i=0;i<m;i++) 
     y[i]=&x[i*n*ny];

  return (void **)y;
}

/** Delallocation of a 2d array.
    Deallocates the memory area pointed to by x, that was allocated
    by a previous malloc2d or calloc2d. If x is null, no action occurs.
    @param x [IN] the 2d array to be allocated    
**/
void free2d(
  char **x
  )
{
  /* sanity check */
  if(x==NULL) return;

  free((char *)&(x[0][0])); 
  free((char *)x); 
  x = NULL;
} 

double RandomUniform(void)
{  
  return (double)random()/(double)RAND_MAX; 
}
