#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gnuplot_i.h"

double diffx_p(double *, double*, double*, int);

int main(){
	double x[100];
	double z[100][100],ddz[100][100];
	int i;
	
	gnuplot_ctrl *h1;

	for(i=0;i<100;i++){
			x[i] = i/10.;
		for(j=0;j<100;j++){
			y[j] = j/10.;
			z[i][j] = sin(x[i])+cos(y[j]);
		}
	}

	diffx_p(x,z,dz,100);
//	diffx_p(x,dz,ddz,100);
	h1=gnuplot_init();
	gnuplot_resetplot(h1);
	gnuplot_plot_xy(h1,x,z,100,"sin(x)");
	gnuplot_plot_xy(h1,x,dz,100,"d sin(x)/dx");
//	gnuplot_plot_xy(h1,x,ddz,100,"dd sin(x)/ddx");

	getchar();
	gnuplot_close(h1);	

	return 0;
}

//calculates the 2nd derivative in x and y direction of a two dim matrix.
//
double diffxy2_p(double *x, double *y, double *z,double *ddx, double *ddy,int N){
	double h;
	int i;
	// calculate first in x-direction (j)
	h = (x[N-1]-x[0])/N;

	for(i=0;i<N;i++){
		for(j=1;j<(N-1);j++){
			ddx[i][j] = (z[i][j+1]-2*z[i][j]+z[i][j-1])/(h*h);	
			ddy[j][i] = (z[j+1][i]-2*z[j][i]+z[j-1][i])/(h*h);	
		}
		ddx[i][0] = (z[i][1]-2*z[i][0]+z[i][N-1])/(h*h);	
		ddx[i][N-1] =(z[i][0]-2*z[i][N-1]+z[i][N-2])/(h*h);	

		ddy[0][i] = (z[1][i]-2*z[0][i]+z[N-1][i])/(h*h);	
		ddy[N-1][i] =(z[0][i]-2*z[N-1][i]+z[N-2][i])/(h*h);	
	}


}
