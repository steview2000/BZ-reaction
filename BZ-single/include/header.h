#include <stdio.h>
#include <stdlib.h>

//#define DATASET "1301170S"
#define SKIP 400 //at SKIP steps an image is saved
#define ASTART 0.5
#define LIMIT 1e-3
#define TSTART 0.0
#define TEND 100
#define NX 128//1024//256
#define NY 128//1024//256
#define PI 3.14159

#define k_1 3.0
#define k_12 -1.0

#define k_2 3.0
#define k_21 -1.0

#define EPSILON 0.08
#define Q 0.0015
#define F 1.5

#define DA 1.0 // diffusion of activator
#define DC 0.5 //diffusion of inhibitor

// steps in space and time
#define DT 0.001
#define dx 0.2

// corecting factor to safe tiff image with a nearly max intensity

// coupling strenght
#define I0 0.00

// instant output
#define GnuOut 0

// which startconditions
#define	StartCondition 1// 1- cross, 2- random, 3- homogeneous,default is taking the latest file

#define f_f 0.0 //forcing frequency

#define SPEED 1.0
#define DATASET "16081821" 
