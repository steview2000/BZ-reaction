#include <time.h>
#include "../include/header.h"


// calculates the derivative of a matrix with respect to x
// periodic boundary condictions are assumed
// N is the length of z

double max(double **x){
	double mx=0;
	int i,j;
	for(i=0;i<NY;i++){
		for(j=0;j<NX;j++){
			if(mx<x[i][j]) mx=x[i][j];
		}
	}

	return mx;
}

void diffxy2_p(double *x, double *y, double **z,double *ddx[], double *ddy[],int N){
	int i,j;
	double h;
	h=dx;
	// calculate first in x-direction (j)
/*	if(h != (x[N-1]-x[0])/(N-1)){
		printf("Problem with \"h\"!\n");
		exit(1);
	};
*/
	for(i=0;i<N;i++){
		for(j=1;j<(N-1);j++){
			ddx[i][j] = (z[i][j+1]-2*z[i][j]+z[i][j-1])/(h*h);	
			ddy[j][i] = (z[j+1][i]-2*z[j][i]+z[j-1][i])/(h*h);	
		}
		ddx[i][0] = (z[i][1]-2*z[i][0]+z[i][N-1])/(h*h);	
		ddx[i][N-1] =(z[i][0]-2*z[i][N-1]+z[i][N-2])/(h*h);	

		ddy[0][i] = (z[1][i]-2*z[0][i]+z[N-1][i])/(h*h);	
		ddy[N-1][i] =(z[0][i]-2*z[N-1][i]+z[N-2][i])/(h*h);	
	};
};

	
// calculates the derivative of a matrix with respect to x
// zero boundary condictions are assumed
// N is the length of z
	
void diffxy2_cOld(double *x, double *y, double **z,double *ddx[], double *ddy[],int N){
	int i,j;
	double h;
	h=dx;
	for(i=0;i<N;i++){
		for(j=1;j<(N-1);j++){
			ddx[i][j] = (z[i][j+1]-2*z[i][j]+z[i][j-1])/(h*h);	
			ddy[j][i] = (z[j+1][i]-2*z[j][i]+z[j-1][i])/(h*h);	
		}
		ddx[i][0] = (z[i][1]-z[i][0])/(h*h);	
		ddx[i][N-1] = (-z[i][N-1]+z[i][N-2])/(h*h);	


		ddy[0][i] =	 (z[1][i]-z[0][i])/(h*h);	
		ddy[N-1][i] = (-z[N-1][i]+z[N-2][i])/(h*h);	
	};
};



void diffxy2_c(double *x, double *y, double **z,double *ddxy2[],int N){
	int i,j;
	double h;
	h=dx;
	// calculate first in x-direction (j)
	// this is the same as diffxy2_Old but it is almost twice as fast
	for(i=1;i<(N-1);i++){
		for(j=1;j<(N-1);j++){
			ddxy2[i][j] = (z[i][j+1]-4*z[i][j]+z[i][j-1]+z[i+1][j]+z[i-1][j] )/(h*h);	
		}
		ddxy2[i][0] = (z[i][1]-3*z[i][0]+z[i+1][0]+z[i-1][0] )/(h*h);	
		ddxy2[i][N-1] = (-3*z[i][N-1]+z[i][N-2]+z[i+1][N-1]+z[i-1][N-1] )/(h*h);	

	};
	for(j=1;j<(N-1);j++){
		ddxy2[0][j] = (z[0][j+1]-3*z[0][j]+z[0][j-1]+z[1][j])/(h*h);	
	}
	ddxy2[0][0] = (z[0][1]-3*z[0][0]+z[1][0]+z[0][0] )/(h*h);	
	ddxy2[0][N-1] = (-3*z[0][N-1]+z[0][N-2]+z[1][N-1]+z[0][N-1] )/(h*h);	

	for(j=1;j<(N-1);j++){
		ddxy2[N-1][j] = (z[N-1][j+1]-4*z[N-1][j]+z[N-1][j-1]+z[N-1][j]+z[N-2][j] )/(h*h);	
	}
	ddxy2[N-1][0] = (z[N-1][1]-3*z[N-1][0]+z[N-1][0]+z[N-2-1][0] )/(h*h);	
	ddxy2[N-1][N-1] = (-3*z[N-1][N-1]+z[N-1][N-2]+z[N-1][N-1]+z[N-2][N-1] )/(h*h);	

};

// Inhibitor
double cf(double A, double C){
	double d;
//	printf("e in function %lf\n",e);
//	d = A*C*C-k*C; // for Turing
	d = A-C;
	return d;
};


// Activator
double af(double A, double C,double IF, char flag, double Cmax){
	double d;
	double f;
	if(flag=='1') f=F1;
	else if(flag=='2') f=F2;
	else {
		printf("Error with f!!\n\n");
		exit(1);
	};

//	printf("e in function %lf\n",e);
//	d = -A*C*C+f*(1-A); // for Turing
	d = 1.0/EPSILON*(A-A*A-(f*C+(Cmax-C)*IF)*(A-Q)/(A+Q));
	return d;
};

void af2(double *x, double *y, double **a0,double **c0,double **IF, double **a, int N,char flag,double Cmax){
	int i,j;
	double k1,k2,k3,k4;
	double alpha,sp;
	double aM[N],bM[N],cM[N],r[N],u[N];
	volatile double dt,temp,h;
	
	temp = dx;
	h=temp;
	temp = DT;
	dt = temp;	

	if(flag=='1') sp=1.;
	else if(flag=='2') sp = SPEED2;
	else {
		printf("Error with f!!\n\n");
		exit(1);
	};
	alpha = sp*DA*dt/(h*h);
	// Crank-Nicholson to solve the diffusion problem (based on eq. 20.3.16 of numerical recipes)
	// Defining the 3 vectors of the triang matric:
	
		
	
	// calculate first in x-direction (j 461.890,  75.3157)
	// this is the same as diffxy2_Old but it is almost twice as fast
	
	// 
	for(i=1;i<(N-1);i++){
		for(j=1;j<(N-1);j++){
			a[i][j] = a0[i][j]+ (alpha/6)* (a0[i-1][j+1]+a0[i+1][j+1]+a0[i-1][j-1]+a0[i+1][j-1]-20*a0[i][j]+4*a0[i][j-1]+4*a0[i+1][j]+4*a0[i-1][j]+4*a0[i][j+1] );
		}
		a[i][0] = a0[i][0]+ (alpha/6)* (a0[i-1][1]+a0[i+1][1]+a0[i-1][0]+a0[i+1][0]-20*a0[i][0]+4*a0[i][0]+4*a0[i+1][0]+4*a0[i-1][0] + 4*a0[i][1]);
		a[i][N-1] = a0[i][N-1]+ (alpha/6)* (a0[i-1][N-1]+a0[i+1][N-1]+a0[i-1][N-2]+a0[i+1][N-2]-20*a0[i][N-1]+4*a0[i][N-2]+4*a0[i+1][N-1]+4*a0[i-1][N-1]+4*a0[i][N-1] );
	};
	for(j=1;j<(N-1);j++){
		a[0][j] = a0[0][j]+ alpha/6* (a0[0][j+1]+a0[1][j+1]+a0[0][j-1]+a0[1][j-1]-20*a0[0][j]+4*a0[0][j-1]+4*a0[1][j]+4*a0[0][j] +4*a0[0][j+1]);
		a[N-1][j] = a0[N-1][j]+ alpha/6* (a0[N-2][j+1]+a0[N-1][j+1]+a0[N-2][j-1]+a0[N-1][j-1]-20*a0[N-1][j]+4*a0[N-1][j-1]+4*a0[N-1][j]+4*a0[N-2][j] +4*a0[N-1][j+1]);
	}
	a[0][0] = a0[0][0]+ alpha/6* (a0[0][1]+a0[1][1]+a0[0][0]+a0[1][0]-20*a0[0][0]+4*a0[0][0]+4*a0[1][0]+4*a0[0][0] +4*a0[0][1])/(h*h);
	a[0][N-1] = a0[0][N-1]+ dt* (a0[0][N-1]+a0[1][N-1]+a0[0][N-2]+a0[1][N-2]-20*a0[0][N-1]+4*a0[0][N-2]+4*a0[1][N-1]+4*a0[0][N-1]+4*a0[0][N-1] );

	a[N-1][0] = a0[N-1][0]+ alpha/6* (a0[N-2][1]+a0[N-1][1]+a0[N-2][0]+a0[N-1][0]-20*a0[N-1][0]+4*a0[N-1][0]+4*a0[N-1][0]+4*a0[N-2][0] +4*a0[N-1][1]);
	a[N-1][N-1] = a0[N-1][N-1]+ alpha/6* (a0[N-2][N-1]+a0[N-1][N-1]+a0[N-2][N-2]+a0[N-1][N-2]-20*a0[N-1][N-1]+4*a0[N-1][N-2]+4*a0[N-1][N-1]+4*a0[N-2][N-1] + 4*a0[N-1][N-1] );

	
	// the following assumes boundaries with a=0 and c=0 at the boundaries:
//	for(i=1;i<(N-1);i++){
//		for(j=1;j<(N-1);j++){
//			a[i][j] = a0[i][j]+ (alpha/6)* (a0[i-1][j+1]+a0[i+1][j+1]+a0[i-1][j-1]+a0[i+1][j-1]-20*a0[i][j]+4*a0[i][j-1]+4*a0[i+1][j]+4*a0[i-1][j]+4*a0[i][j+1] );	
//		}
//		a[i][0] = a0[i][0]+ (alpha/6)* (a0[i-1][1]+a0[i+1][1]-20*a0[i][0]+4*a0[i+1][0]+4*a0[i-1][0]+4*a0[i][1]);
//		a[i][N-1] = a0[i][N-1]+ (alpha/6)* (a0[i-1][N-2]+a0[i+1][N-2]-20*a0[i][N-1]+4*a0[i][N-2]+4*a0[i+1][N-1]+4*a0[i-1][N-1]);	
//	};
//	for(j=1;j<(N-1);j++){
//		a[0][j] = a0[0][j]+ (alpha/6)* (a0[1][j+1]+a0[1][j-1]-20*a0[0][j]+4*a0[0][j-1]+4*a0[1][j]+4*a0[0][j+1]);	
//		a[N-1][j] = a0[N-1][j]+ (alpha/6)* (a0[N-2][j+1]+a0[N-2][j-1]-20*a0[N-1][j]+4*a0[N-1][j-1]+4*a0[N-2][j]+4*a0[N-1][j+1] );	
//	}
//	a[0][0] = a0[0][0]+ (alpha/6)* (a0[1][1]-20*a0[0][0]+4*a0[1][0]+4*a0[0][1] );	
//	a[0][N-1] = a0[0][N-1]+ (alpha/6)* (a0[1][N-2]-20*a0[0][N-1]+4*a0[0][N-2]+4*a0[1][N-1]);	
//	
//	a[N-1][0] = a0[N-1][0]+ (alpha/6)* (a0[N-2][1]-20*a0[N-1][0]+4*a0[N-2][0]+4*a0[N-1][1] );	
//	a[N-1][N-1] = a0[N-1][N-1]+ (alpha/6)* (a0[N-2][N-2]-20*a0[N-1][N-1]+4*a0[N-1][N-2]+4*a0[N-2][N-1]);	

	// 4th order Runge-Kutta for the reaction term:
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			k1 = sp*dt*af(a0[i][j],c0[i][j],IF[i][j],flag, Cmax);
		//	k2 = dt*af(a0[i][j]+0.5*k1,c0[i][j],IF[i][j]);
		//	k3 = dt*af(a0[i][j]+0.5*k2,c0[i][j],IF);
		//	k4 = dt*af(a0[i][j]+k3,c0[i][j],IF);

			a[i][j] = a[i][j]+k1;//1/6.0*k1+1/3.0*k2+1/3.0*k3+1/6.0*k4;
		//	Imean=Imean+IF[i][j];
		//	if(IF[i][j]>0) printf("i: %i\t j: %lf IF: %lf\t\n",i,j,IF[i][j]);
		}
	}
//	printf("Imean: %lf\t",Imean);

}

void cf2(double *x, double *y, double **a0,double **c0, double **c, int N,char flag){
	int i,j;
	double k1,k2,k3,k4;
	double alpha,sp;
	double aM[N],bM[N],cM[N],r[N],u[N];
	volatile double dt,temp,h;
	
	temp = dx;
	h=temp;
	temp = DT;
	dt = temp;	

	if(flag=='1') sp=1.;
	else if(flag=='2') sp = SPEED2;
	else {
		printf("Error with f!!\n\n");
		exit(1);
	};

	alpha = sp*DC*dt/(h*h);
	// Crank-Nicholson to solve the diffusion problem (based on eq. 20.3.16 of numerical recipes)
	// calculate first in x-direction (j 461.890,  75.3157)
	// this is the same as diffxy2_Old but it is almost twice as fast
	
	for(i=1;i<(N-1);i++){
		for(j=1;j<(N-1);j++){
			c[i][j] = c0[i][j]+ (alpha/6)* (c0[i-1][j+1]+c0[i+1][j+1]+c0[i-1][j-1]+c0[i+1][j-1]-20*c0[i][j]+4*c0[i][j-1]+4*c0[i+1][j]+4*c0[i-1][j]+4*c0[i][j+1] );
		}

		c[i][0] = c0[i][0]+ alpha/6* (c0[i-1][1]+c0[i+1][1]+c0[i-1][0]+c0[i+1][0]-20*c0[i][0]+4*c0[i][0]+4*c0[i+1][0]+4*c0[i-1][0] + 4*c0[i][1]);
		c[i][N-1] = c0[i][N-1]+ alpha/6* (c0[i-1][N-1]+c0[i+1][N-1]+c0[i-1][N-2]+c0[i+1][N-2]-20*c0[i][N-1]+4*c0[i][N-2]+4*c0[i+1][N-1]+4*c0[i-1][N-1]+4*c0[i][N-1] );
	};
	
	for(j=1;j<(N-1);j++){
		c[0][j] = c0[0][j]+ (alpha/6)* (c0[0][j+1]+c0[1][j+1]+c0[0][j-1]+c0[1][j-1]-20*c0[0][j]+4*c0[0][j-1]+4*c0[1][j]+4*c0[0][j] +4*c0[0][j+1]);	
		c[N-1][j] = c0[N-1][j]+ alpha/6* (c0[N-2][j+1]+c0[N-1][j+1]+c0[N-2][j-1]+c0[N-1][j-1]-20*c0[N-1][j]+4*c0[N-1][j-1]+4*c0[N-1][j]+4*c0[N-2][j] +4*c0[N-1][j+1]);	

	}
		c[0][0] = c0[0][0]+ alpha/6* (c0[0][1]+c0[1][1]+c0[0][0]+c0[1][0]-20*c0[0][0]+4*c0[0][0]+4*c0[1][0]+4*c0[0][0] +4*c0[0][1])/(h*h);	
		c[0][N-1] = c0[0][N-1]+ dt* (c0[0][N-1]+c0[1][N-1]+c0[0][N-2]+c0[1][N-2]-20*c0[0][N-1]+4*c0[0][N-2]+4*c0[1][N-1]+4*c0[0][N-1]+4*c0[0][N-1] );	

		c[N-1][0] = c0[N-1][0]+ alpha/6* (c0[N-2][1]+c0[N-1][1]+c0[N-2][0]+c0[N-1][0]-20*c0[N-1][0]+4*c0[N-1][0]+4*c0[N-1][0]+4*c0[N-2][0] +4*c0[N-1][1]);	
		c[N-1][N-1] = c0[N-1][N-1]+ alpha/6* (c0[N-2][N-1]+c0[N-1][N-1]+c0[N-2][N-2]+c0[N-1][N-2]-20*c0[N-1][N-1]+4*c0[N-1][N-2]+4*c0[N-1][N-1]+4*c0[N-2][N-1] + 4*c0[N-1][N-1] );	

	// the following assumes boundaries with a=0 and c=0 at the boundaries:
//	for(i=1;i<(N-1);i++){
//		for(j=1;j<(N-1);j++){
//			c[i][j] = c0[i][j]+ alpha/6* (c0[i-1][j+1]+c0[i+1][j+1]+c0[i-1][j-1]+c0[i+1][j-1]-20*c0[i][j]+4*c0[i][j-1]+4*c0[i+1][j]+4*c0[i-1][j]+4*c0[i][j+1] );	
//		}
//		c[i][0] = c0[i][0]+ alpha/6* (c0[i-1][1]+c0[i+1][1]+c0[i-1][0]+c0[i+1][0]-20*c0[i][0]+4*c0[i][0]+4*c0[i+1][0]+4*c0[i-1][0] + 4*c0[i][1]);	
//		c[i][N-1] = c0[i][N-1]+ alpha/6* (c0[i-1][N-1]+c0[i+1][N-1]+c0[i-1][N-2]+c0[i+1][N-2]-20*c0[i][N-1]+4*c0[i][N-2]+4*c0[i+1][N-1]+4*c0[i-1][N-1]+4*c0[i][N-1] );	
//	};
//	for(j=1;j<(N-1);j++){
//		c[0][j] = c0[0][j]+ alpha/6* (c0[0][j+1]+c0[1][j+1]+c0[0][j-1]+c0[1][j-1]-20*c0[0][j]+4*c0[0][j-1]+4*c0[1][j]+4*c0[0][j] +4*c0[0][j+1]);	
//	}
//		c[0][0] = c0[0][0]+ alpha/6* (c0[0][1]+c0[1][1]+c0[0][0]+c0[1][0]-20*c0[0][0]+4*c0[0][0]+4*c0[1][0]+4*c0[0][0] +4*c0[0][1])/(h*h);	
//		c[0][N-1] = c0[0][N-1]+ dt* (c0[0][N-1]+c0[1][N-1]+c0[0][N-2]+c0[1][N-2]-20*c0[0][N-1]+4*c0[0][N-2]+4*c0[1][N-1]+4*c0[0][N-1]+4*c0[0][N-1] );	
//
//	for(j=1;j<(N-1);j++){
//		c[N-1][j] = c0[N-1][j]+ alpha/6* (c0[N-2][j+1]+c0[N-1][j+1]+c0[N-2][j-1]+c0[N-1][j-1]-20*c0[N-1][j]+4*c0[N-1][j-1]+4*c0[N-1][j]+4*c0[N-2][j] +4*c0[N-1][j+1]);	
//
//	}
//		c[N-1][0] = c0[N-1][0]+ alpha/6* (c0[N-2][1]+c0[N-1][1]+c0[N-2][0]+c0[N-1][0]-20*c0[N-1][0]+4*c0[N-1][0]+4*c0[N-1][0]+4*c0[N-2][0] +4*c0[N-1][1]);	
//		c[N-1][N-1] = c0[N-1][N-1]+ alpha/6* (c0[N-2][N-1]+c0[N-1][N-1]+c0[N-2][N-2]+c0[N-1][N-2]-20*c0[N-1][N-1]+4*c0[N-1][N-2]+4*c0[N-1][N-1]+4*c0[N-2][N-1] + 4*c0[N-1][N-1] );	
//

	// 4th order Runge-Kutta for the reaction term:
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			k1 = sp*dt*cf(a0[i][j],c0[i][j]);
		//	k2 = dt*af(a0[i][j]+0.5*k1,c0[i][j],IF[i][j]);
		//	k3 = dt*af(a0[i][j]+0.5*k2,c0[i][j],IF);
		//	k4 = dt*af(a0[i][j]+k3,c0[i][j],IF);

			c[i][j] = c[i][j]+k1;//1/6.0*k1+1/3.0*k2+1/3.0*k3+1/6.0*k4;
		//	Imean=Imean+IF[i][j];
		//	if(IF[i][j]>0) printf("i: %i\t j: %lf IF: %lf\t\n",i,j,IF[i][j]);
		}
	}
//	printf("Imean: %lf\t",Imean);

}



void randSpiral(double **a0, double **c0,double **a, double **c){
	// this function takes the previous image and randomizes it
 	// create 8x8 windows
	int NumBox = NX/32;
		
 	int pos[NumBox*NumBox],temp,lmax,kmax,posx,posy;
	int i,j;
	int k,l;
	time_t seconds;
	time(&seconds);
	sleep(1);
	srand((unsigned int) seconds);
	// create set of shuffeld numbers from 0 to 15
	for(i=0;i<(NumBox*NumBox);i++){
		pos[i] = i;
//		posy[i] = i;
	}
		
	for(i=0;i<(NumBox*NumBox);i++){
		j = (rand()%(NumBox*NumBox));
//		k = (rand()%16);
		temp = pos[i]; 
		pos[i] = pos[j];
		pos[j] = temp;
	} 

	for(i=0;i<(NumBox*NumBox);i++){
		j = (rand()%NumBox);
//		k = (rand()%16);
		temp = pos[i]; 
		pos[i] = pos[j];
		pos[j] = temp;
	} 
	
//	for(i=0;i<(16*16);i++){
//		printf("%d\t",pos[i]);
//	} 
	printf("\n");

	kmax = NY/NumBox;
	lmax = NX/NumBox;

	for(i=0;i<NumBox;i++){
		for(j=0;j<NumBox;j++){
			for(k=0;k<kmax;k++){
				for(l=0;l<lmax;l++){
					posy = (pos[j*NumBox+i]/NumBox);
					posx = (pos[j*NumBox+i]%NumBox);
					a[i*kmax+k][j*lmax+l] = a0[posx*kmax+k][posy*lmax+l];	
					c[i*kmax+k][j*lmax+l] = c0[posx*kmax+k][posy*lmax+l];	
				}
			}
		}
	}
 
}





//routine to solve tridiagonal system of equation (from Numerical recipes, $2.4)
void tridag(double *a, double *b, double *c, double *r, double *u,int n){
	int j;
	double bet;
	double gam[n];
	
	if(b[0] == 0.0){
	 	printf("Error 1 in tridag!\n");
		exit(1);	
	}
	u[0]=r[0]/(bet=b[0]);
	for(j=1;j<n;j++){
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if(bet==0.0) printf("Error 2 in tridag");
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for(j=(n-2);j>=0;j--)
		u[j]-=gam[j+1]*u[j+1];

}

void spiralGen(double **a0, double **c0,double ss,double umax,double vmax, double umin, double vmin){
	int	i,j,M;
	M = NX/2;
	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			a0[i][j] = ss;
			c0[i][j] = ss;
		};
	};
	for(j=0;j<0.8*M;j++){
		for(i=M;i<(M+40);i++){
			a0[i][j] = umax;
			c0[i][j] = vmax;
		}
		for (i=(M-80);i<M;i++){
			a0[i][j] = umin;
			c0[i][j] = vmax;//(1+ss)/2.;
		}
	};
}
