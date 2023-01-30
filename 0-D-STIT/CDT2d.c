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
**/
void InitCDT2d(CDT2d *cdt){
	unsigned long is;
	double winSize = (double) (cdt->size);
	
	// random seed is initialized.
 	srandom(time(NULL));
 	
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
**/
int FreeCDT2d(CDT2d *cdt){
	unsigned long i;
	
	Free(cdt->nSize);
	
	for(i = 0; i < cdt->numberOfPolygons; i++){
		Free2d(((cdt->Polygon) + i)->V);
		Free2d(((cdt->Polygon) + i)->H);
		Free(((cdt->Polygon) + i)->B);
	}
	Free(cdt->Polygon);
	
	Free(cdt->I_Segment);
	Free(cdt->Phi);
	Free(cdt->PhiCum);
	
	return 0;
}

/**
@function SetPhiIso

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
	#define EPS 1.0E-12
	// checking if the sum of the probabilities is 1.
	for(iDir = 0; iDir < nDir; iDir++){
		sumProb += probDir[iDir];
	}
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
**/
void InitTime(CDT2d *cdt, unsigned long lOption){
	unsigned long j;
	double bPhi = 0.0;
	
	// The lifetime of the polygon is determined. The lifetime
	// is exponentially distribuited with parameter bPhi, which is
	// the weighted sum of the width and the weights are the probabilities
	// Phi[j].
	
	//Perimeter
	if(lOption == 0){
		for(j = 0; j < NUMBER_OF_INT; j++){
			bPhi += cdt->Polygon[0].B[j]*cdt->Phi[j];
		}
		cdt->Polygon[0].tau = RandomExponential(bPhi);
	}
	//Area
	else{
		cdt->Polygon[0].tau = RandomExponential(cdt->Polygon[0].area);
	}
	cdt->Polygon[0].beta = 0.0;
	cdt->Polygon[0].time = cdt->Polygon[0].beta + cdt->Polygon[0].tau;
}

/**
@function RandomDirection

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
@function RandomLinePolygon

Isotropic, Anisotropic discrete
**/
LINE RandomLinePolygon(CDT2d *cdt, unsigned long i){
	unsigned long j;
	double d, randNum;
	LINE dirLine;
	
	do{
		randNum = RandomUniform();
		// random direction index.
		j = RandomDirection(cdt, randNum);
		// printf("j: %ld\n",j);
		// the signed distance dist is computed.
		d = cdt->Polygon[i].maxWidth * RandomUniform() + cdt->Polygon[i].H[j][0];
	// checking if the line is accepted.		
	}while(d > cdt->Polygon[i].H[j][1]); 

	// if the line is accepted, angle alpha is computed.
	dirLine.alpha = M_PI * (double)j / (double)NUMBER_OF_INT;
	dirLine.dist = d;
	dirLine.index = j;

	return dirLine;
}

/**
@function RandomLinePolygonDisturbed

Anisotropic discrete disturbed, Anisotropic Elliptic
**/
LINE RandomLinePolygonDisturbed(CDT2d *cdt, unsigned long i, double bEllip){
	unsigned long j, jGamma;
	double d, randNum, gamma;
	LINE dirLine;
	
	do{
		randNum = RandomUniform();
		// random direction index.
		j = RandomDirection(cdt, randNum);
		// printf("j: %ld\n",j);
		gamma = RandomElliptic(bEllip);
		jGamma = floor(gamma * NUMBER_OF_INT / M_PI);
		j = (j + jGamma) % (NUMBER_OF_INT);
		// the signed distance dist is computed.
		d = cdt->Polygon[i].maxWidth * RandomUniform() + cdt->Polygon[i].H[j][0];
	// checking if the line is accepted.		
	}while(d > cdt->Polygon[i].H[j][1]); 

	// if the line is accepted, angle alpha is computed.
	dirLine.alpha = M_PI * (double)j / (double)NUMBER_OF_INT;
	dirLine.dist = d;
	dirLine.index = j;

	return dirLine;
}

/**
@function RandomIntersectionPolygon
**/
#define EPS 1.0E-32
void RandomIntersectionPolygon(CDT2d *cdt, unsigned long i, unsigned long ip, LINE dirLine, unsigned long lOption){
	unsigned long j, k, k1, k2;
	// k0[1] - k0[0] is the number of vertices between the first and the second intersection point.
	unsigned long k0[2]; 
	// ind means index of intersection points (0 and 1).
	unsigned long ind = 0;
	double cosAlpha, det, deltaX, deltaY, sinAlpha, t, alpha, d, bPhi;
	double newVertex[2][2];

	// auxiliary variable for checking if the angle of the intersection between I_Segment and
	// the side of the polygon where the intersection point lies, is not close to zero. 
	// DETcheck = 1 (angle is ok) and DETcheck = 0 (angle is close to zero or equal to zero for
	// all the side of the polygons where an intersection point could be). Then, if DETcheck = 0
	// a new random line must be generated.
	int DETcheck = 1;

	POLYGON P;
  
	// copy of Polygon in to polygon P.
	P.numberOfVertices = cdt->Polygon[i].numberOfVertices;
	P.V = (double **)calloc2d(P.numberOfVertices, 2, sizeof(double));
	for(k = 0; k < P.numberOfVertices; k++) {
		P.V[k][0] = cdt->Polygon[i].V[k][0];
		P.V[k][1] = cdt->Polygon[i].V[k][1];
	}
	
	P.time = cdt->Polygon[i].time;
	
	alpha = dirLine.alpha;
	d = dirLine.dist;
  
	do{
		if (DETcheck == 0){	
			dirLine = RandomLinePolygon(cdt, i);
			alpha = dirLine.alpha;
			d = dirLine.dist;
			ind = 0;
			DETcheck = 1;
		}
		cosAlpha = cos(alpha);
		sinAlpha = sin(alpha);
	
		for(k = 0; k < P.numberOfVertices; k++) {
			k1 = (k + 1) % P.numberOfVertices;
			deltaX = P.V[k1][0] - P.V[k][0];
			deltaY = P.V[k1][1] - P.V[k][1];
			det = deltaX * cosAlpha + deltaY * sinAlpha;
			if(fabs(det) > EPS) {
				// the segment between vertices k and (k + 1) is parametrized 
				// with parameter t in [0,1]. The parameter value of the 
				// intersection point is computed.
				t = (d - P.V[k][0] * cosAlpha - P.V[k][1] * sinAlpha) / det;
			  	if(t >= 0.0 && t < 1.0){
					newVertex[ind][0] = P.V[k][0] + t * deltaX;
					newVertex[ind][1] = P.V[k][1] + t * deltaY;
					k0[ind] = k;
					ind++; 
				}				
			}
		}
		if (ind < 2)
			// the line is rejected, because there are less than 2 intersection points. 
			DETcheck = 0;
			
		if (DETcheck == 1) break;
		
	}while(1);

	Free2d(cdt->Polygon[i].V);
	Free2d(cdt->Polygon[i].H);
	Free(cdt->Polygon[i].B);

	// vertices of the new two polygons are determined.
	cdt->Polygon[i].V = (double **)calloc2d(k0[1] - k0[0] + 2, 2, sizeof(double));
	cdt->Polygon[i].V[0][0] = newVertex[1][0];
	cdt->Polygon[i].V[0][1] = newVertex[1][1];
	cdt->Polygon[i].V[1][0] = newVertex[0][0];
	cdt->Polygon[i].V[1][1] = newVertex[0][1];
	cdt->Polygon[i].numberOfVertices = 2;

	cdt->Polygon[ip].V = (double **)calloc2d(P.numberOfVertices - (k0[1] - k0[0]) + 2, 2, sizeof(double));
	cdt->Polygon[ip].V[0][0] = newVertex[0][0];
	cdt->Polygon[ip].V[0][1] = newVertex[0][1];
	cdt->Polygon[ip].V[1][0] = newVertex[1][0];
	cdt->Polygon[ip].V[1][1] = newVertex[1][1];
	cdt->Polygon[ip].numberOfVertices = 2;
  
	for(k = 0; k < P.numberOfVertices; k++) {
		k1 = (k + k0[0] + 1) % P.numberOfVertices;
		k2 = (k + k0[1] + 1) % P.numberOfVertices;
		if(cdt->Polygon[i].numberOfVertices < k0[1] - k0[0] + 2){
			cdt->Polygon[i].V[cdt->Polygon[i].numberOfVertices][0] = P.V[k1][0];
			cdt->Polygon[i].V[cdt->Polygon[i].numberOfVertices][1] = P.V[k1][1];
			cdt->Polygon[i].numberOfVertices++;
		}
		if(cdt->Polygon[ip].numberOfVertices < P.numberOfVertices - (k0[1] - k0[0]) + 2){
			cdt->Polygon[ip].V[cdt->Polygon[ip].numberOfVertices][0] = P.V[k2][0];
			cdt->Polygon[ip].V[cdt->Polygon[ip].numberOfVertices][1] = P.V[k2][1];
			cdt->Polygon[ip].numberOfVertices++; 
		}
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
	cdt->I_Segment[cdt->numberOfI_Segments].x0 = newVertex[0][0];
	cdt->I_Segment[cdt->numberOfI_Segments].y0 = newVertex[0][1]; 
	cdt->I_Segment[cdt->numberOfI_Segments].x1 = newVertex[1][0]; 
	cdt->I_Segment[cdt->numberOfI_Segments].y1 = newVertex[1][1];
	// The length of the I_Segment is computed.
	deltaX = newVertex[1][0] - newVertex[0][0];
	deltaY = newVertex[1][1] - newVertex[0][1];
	cdt->I_Segment[cdt->numberOfI_Segments].length = sqrt(deltaX*deltaX + deltaY*deltaY);
	// The birthtime of the I_Segment is assigned.
	cdt->I_Segment[cdt->numberOfI_Segments].beta = P.time;

	Free2d(P.V);
}
#undef EPS

/*
@function RandomIntersectionPolygonDisturbed
*/
#define EPS 1.0E-32
void RandomIntersectionPolygonDisturbed(CDT2d *cdt,unsigned long i,unsigned long ip, LINE dirLine, double bEllip, unsigned long lOption){
	unsigned long j, k, k1, k2;
	// k0[1] - k0[0] is the number of vertices between the first and the second intersection point.
	unsigned long k0[2]; 
	// ind means index of intersection points (0 and 1).
	unsigned long ind = 0;
	double cosAlpha, det, deltaX, deltaY, sinAlpha, t, alpha, d, bPhi;
	double newVertex[2][2];

	// auxiliary variable for checking if the angle of the intersection between I_Segment and
	// the side of the polygon where the intersection point lies, is not close to zero. 
	// DETcheck = 1 (angle is ok) and DETcheck = 0 (angle is close to zero or equal to zero for
	// all the side of the polygons where an intersection point could be). Then, if DETcheck = 0
	// a new random line must be generated.
	int DETcheck = 1;

	POLYGON P;
  
	// copy of Polygon in to polygon P.
	P.numberOfVertices = cdt->Polygon[i].numberOfVertices;
	P.V = (double **)calloc2d(P.numberOfVertices, 2, sizeof(double));
	for(k = 0; k < P.numberOfVertices; k++) {
		P.V[k][0] = cdt->Polygon[i].V[k][0];
		P.V[k][1] = cdt->Polygon[i].V[k][1];
	}
	
	P.time = cdt->Polygon[i].time;
	
	alpha = dirLine.alpha;
	d = dirLine.dist;
  
	do{
		if (DETcheck == 0){	
			dirLine = RandomLinePolygonDisturbed(cdt, i, bEllip);
			alpha = dirLine.alpha;
			d = dirLine.dist;
			ind = 0;
			DETcheck = 1;
		}
		cosAlpha = cos(alpha);
		sinAlpha = sin(alpha);
	
		for(k = 0; k < P.numberOfVertices; k++) {
			k1 = (k + 1) % P.numberOfVertices;
			deltaX = P.V[k1][0] - P.V[k][0];
			deltaY = P.V[k1][1] - P.V[k][1];
			det = deltaX * cosAlpha + deltaY * sinAlpha;
			if(fabs(det) > EPS) {
				// the segment between vertices k and (k + 1) is parametrized 
				// with parameter t in [0,1]. The parameter value of the 
				// intersection point is computed.
				t = (d - P.V[k][0] * cosAlpha - P.V[k][1] * sinAlpha) / det;
			  	if(t >= 0.0 && t < 1.0){
					newVertex[ind][0] = P.V[k][0] + t * deltaX;
					newVertex[ind][1] = P.V[k][1] + t * deltaY;
					k0[ind] = k;
					ind++; 
				}				
			}
		}
		if (ind < 2)
			// the line is rejected, because there are less than 2 intersection points. 
			DETcheck = 0;
			
		if (DETcheck == 1) break;
		
	}while(1);

	Free2d(cdt->Polygon[i].V);
	Free2d(cdt->Polygon[i].H);
	Free(cdt->Polygon[i].B);

	// vertices of the new two polygons are determined.
	cdt->Polygon[i].V = (double **)calloc2d(k0[1] - k0[0] + 2, 2, sizeof(double));
	cdt->Polygon[i].V[0][0] = newVertex[1][0];
	cdt->Polygon[i].V[0][1] = newVertex[1][1];
	cdt->Polygon[i].V[1][0] = newVertex[0][0];
	cdt->Polygon[i].V[1][1] = newVertex[0][1];
	cdt->Polygon[i].numberOfVertices = 2;

	cdt->Polygon[ip].V = (double **)calloc2d(P.numberOfVertices - (k0[1] - k0[0]) + 2, 2, sizeof(double));
	cdt->Polygon[ip].V[0][0] = newVertex[0][0];
	cdt->Polygon[ip].V[0][1] = newVertex[0][1];
	cdt->Polygon[ip].V[1][0] = newVertex[1][0];
	cdt->Polygon[ip].V[1][1] = newVertex[1][1];
	cdt->Polygon[ip].numberOfVertices = 2;
  
	for(k = 0; k < P.numberOfVertices; k++) {
		k1 = (k + k0[0] + 1) % P.numberOfVertices;
		k2 = (k + k0[1] + 1) % P.numberOfVertices;
		if(cdt->Polygon[i].numberOfVertices < k0[1] - k0[0] + 2){
			cdt->Polygon[i].V[cdt->Polygon[i].numberOfVertices][0] = P.V[k1][0];
			cdt->Polygon[i].V[cdt->Polygon[i].numberOfVertices][1] = P.V[k1][1];
			cdt->Polygon[i].numberOfVertices++;
		}
		if(cdt->Polygon[ip].numberOfVertices < P.numberOfVertices - (k0[1] - k0[0]) + 2){
			cdt->Polygon[ip].V[cdt->Polygon[ip].numberOfVertices][0] = P.V[k2][0];
			cdt->Polygon[ip].V[cdt->Polygon[ip].numberOfVertices][1] = P.V[k2][1];
			cdt->Polygon[ip].numberOfVertices++; 
		}
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
	cdt->I_Segment[cdt->numberOfI_Segments].x0 = newVertex[0][0];
	cdt->I_Segment[cdt->numberOfI_Segments].y0 = newVertex[0][1]; 
	cdt->I_Segment[cdt->numberOfI_Segments].x1 = newVertex[1][0]; 
	cdt->I_Segment[cdt->numberOfI_Segments].y1 = newVertex[1][1];
	// The length of the I_Segment is computed.
	deltaX = newVertex[1][0] - newVertex[0][0];
	deltaY = newVertex[1][1] - newVertex[0][1];
	cdt->I_Segment[cdt->numberOfI_Segments].length = sqrt(deltaX*deltaX + deltaY*deltaY);
	// The birthtime of the I_Segment is assigned.
	cdt->I_Segment[cdt->numberOfI_Segments].beta = P.time;

	Free2d(P.V);
}
#undef EPS


/**
@function STIT2dIso

Isotropic STIT2d
**/
int STIT2dIso(CDT2d *cdt, double timeStop, unsigned long lOption){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	// current line candidate in the rejection process.	
	LINE dirLine;
	
	// phi and cumulative distribution function are computed.
	SetPhiIso(cdt);
	
	//lifetime of the first polygon (window) is computed.
	InitTime(cdt,lOption);
	
	// it assumes that the process continues.
	stitStop = 1;
	do{
		numPol = cdt->numberOfPolygons;
    	q = 0;
    	for(i = 0; i < numPol; i++){
    		// checking if the polygon i must be divided.
      		if(cdt->Polygon[i].time < timeStop){
				dirLine = RandomLinePolygon(cdt, i);
				
				if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
					printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
					return 1;
				}
				if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
					printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
					return 1;
				}
				
				RandomIntersectionPolygon(cdt, i, numPol + q, dirLine, lOption);
				q++;
				(cdt->numberOfPolygons)++;
				(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    	if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return 0;
}

/**
@function STIT2dAniso

Anisotropic STIT2d
**/
int STIT2dAniso(CDT2d *cdt, double timeStop, double *angleDir, double *probDir, unsigned long nDir, unsigned long lOption){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	// current line candidate in the rejection process.	
	LINE dirLine;
	
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
				dirLine = RandomLinePolygon(cdt, i);
				
				if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
					printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
					return 1;
				}
				if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
					printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
					return 1;
				}
				
				RandomIntersectionPolygon(cdt, i, numPol + q, dirLine, lOption);
				q++;
				(cdt->numberOfPolygons)++;
				(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    	if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return 0;
}

/**
@function STIT2dAnisoDisturbed

Anisotropic Disturbed STIT2d
**/
int STIT2dAnisoDisturbed(CDT2d *cdt, double timeStop, double *angleDir, double *probDir, unsigned long nDir, double bEllip, unsigned long lOption){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	// current line candidate in the rejection process.	
	LINE dirLine;
	
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
				dirLine = RandomLinePolygonDisturbed(cdt, i, bEllip);
				
				if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
					printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
					return 1;
				}
				if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
					printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
					return 1;
				}
				
				RandomIntersectionPolygonDisturbed(cdt, i, numPol + q, dirLine, bEllip, lOption);
				q++;
				(cdt->numberOfPolygons)++;
				(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    	if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return 0;
}

/**
@function STIT2dAnisoEllip

Anisotropic Elliptic STIT2d
**/
int STIT2dAnisoEllip(CDT2d *cdt, double timeStop, double bEllip, unsigned long lOption){
	unsigned long i;
	
	// state variable, 0: stit process stops, 1: stit process continues.
	unsigned int stitStop;
	
	// number of polygons in the current level of the branching tree, which are 
	// divided before time_stop. The branching tree describes the cell division 
	// process. Also q is an auxiliary variable for storing new polygons.
	unsigned int q;
	
	// current number of polygons.
	unsigned int numPol;
	
	// current line candidate in the rejection process.	
	LINE dirLine;
	
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
				dirLine = RandomLinePolygonDisturbed(cdt, i, bEllip);
				
				if(cdt->numberOfPolygons == MAX_NUMBER_OF_POLYGONS){
					printf("Number of Polygons exceeds %dld\n", MAX_NUMBER_OF_POLYGONS);
					return 1;
				}
				if(cdt->numberOfI_Segments == MAX_NUMBER_OF_I_SEGMENTS) {
					printf("Number of I_Segments exceeds %dld\n", MAX_NUMBER_OF_I_SEGMENTS);
					return 1;
				}
				
				RandomIntersectionPolygonDisturbed(cdt, i, numPol + q, dirLine, bEllip, lOption);
				q++;
				(cdt->numberOfPolygons)++;
				(cdt->numberOfI_Segments)++;
			}
		}
		// the proccess stops.
    	if(q == 0) stitStop = 0;
  	}while(stitStop);
  	
	return 0;
}

/**	
@function: WidthFunction 

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
	
	Polygon->wML = 1. / ((cdt->size - Polygon->B[0]) * (cdt->size - Polygon->B[64]));

  return 0;
}

/**
@function: PerimeterPolygon
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
**/
int RoundnessPolygon(POLYGON *Polygon){
	
	Polygon->roundness = 4*M_PI*Polygon->area / (Polygon->perimeter * Polygon->perimeter);
	return 0;
}

/**
@function RandomExponential

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
Computes statistics of a STIT.
**/
void StatSTIT(CDT2d *cdt){
	unsigned long i;
	double sumWeights = 0., sumValues, mean, var, std, cv;
	FILE *statFile;
	
	//Writing stats in a file.
	statFile = fopen("statSTIT.txt","w");
	fprintf(statFile,"%ld\n",cdt->numberOfNBPolygons);
	
	//Computing sum(wML)
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumWeights += cdt->NoBoundaryPolygon[i].wML;
	
	//Area
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * cdt->NoBoundaryPolygon[i].area;
	mean = sumValues / sumWeights;
	
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * (cdt->NoBoundaryPolygon[i].area - mean) * (cdt->NoBoundaryPolygon[i].area - mean);
		
	var = sumValues / sumWeights;
	std = sqrt(var);
	cv = std / mean;
		
	fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);
	
	//Perimeter
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * cdt->NoBoundaryPolygon[i].perimeter;
	mean = sumValues / sumWeights;
	
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * (cdt->NoBoundaryPolygon[i].perimeter - mean) * (cdt->NoBoundaryPolygon[i].perimeter - mean);
		
	var = sumValues / sumWeights;
	std = sqrt(var);
	cv = std / mean;
		
	fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

	//Roundness
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * cdt->NoBoundaryPolygon[i].roundness;
	mean = sumValues / sumWeights;
	
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * (cdt->NoBoundaryPolygon[i].roundness - mean) * (cdt->NoBoundaryPolygon[i].roundness - mean);
		
	var = sumValues / sumWeights;
	std = sqrt(var);
	cv = std / mean;

	fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

	//Max width
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * cdt->NoBoundaryPolygon[i].maxWidth;
	mean = sumValues / sumWeights;
	
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * (cdt->NoBoundaryPolygon[i].maxWidth - mean) * (cdt->NoBoundaryPolygon[i].maxWidth - mean);
		
	var = sumValues / sumWeights;
	std = sqrt(var);
	cv = std / mean;
		
	fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

	//Min width
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * cdt->NoBoundaryPolygon[i].minWidth;
	mean = sumValues / sumWeights;
	
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * (cdt->NoBoundaryPolygon[i].minWidth - mean) * (cdt->NoBoundaryPolygon[i].minWidth - mean);
		
	var = sumValues / sumWeights;
	std = sqrt(var);
	cv = std / mean;
		
	fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

	//Aspect ratio
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * cdt->NoBoundaryPolygon[i].aspectRatio;
	mean = sumValues / sumWeights;
	
	sumValues = 0.;
	for(i = 0; i < cdt->numberOfNBPolygons; i++)
		sumValues += cdt->NoBoundaryPolygon[i].wML * (cdt->NoBoundaryPolygon[i].aspectRatio - mean) * (cdt->NoBoundaryPolygon[i].aspectRatio - mean);
		
	var = sumValues / sumWeights;
	std = sqrt(var); 
	cv = std / mean;
		
	fprintf(statFile,"%.12f,%.12f,%.12f\n",mean,std,cv);

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
