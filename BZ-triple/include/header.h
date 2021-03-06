// Header file for BZ-couple
#ifndef BZ_TRIPLE_HEADER
#define BZ_TRIPLE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double af(double,double,double,char);
void af2(double* ,double* ,double**, double **, double**, double**, int, char );
void cf2(double* ,double* ,double**, double **, double**, int);
double max(double **);
void randSpiral(double **,double **,double **,double **);
double cf(double,double);
void spiralGen(double **,double **);


//#define DATASET "1305191S"
#define SKIP 2000 //at SKIP steps an image is saved
#define ASTART 0.5
#define LIMIT 1e-3
#define TSTART 0.0
#define TEND 500
#define NX 580
#define NY 580
#define PI 3.14159

#define k_1 3.0
#define k_12 -1.0

#define k_2 3.0
#define k_21 -1.0


#define EPSILON 0.08 
#define Q 0.0015
#define F1 1.5
#define F2 1.5
#define F3 1.5

#define DA 1.0 // diffusion of activator
#define DC 0.5 //diffusion of inhibitor

// steps in space and time
#define DT 0.001
#define dx 0.2

// corecting factor to safe tiff image with a nearly max intensity
#define corfac 0.5

// which startconditions
#define	StartCondition 2// 1- cross, 2- random, 3- homogeneous,default is taking the latest file


// coupling strenght
//#define I0 0.01


#define DATASET "1402240S"
#define I0 0.01
#endif
