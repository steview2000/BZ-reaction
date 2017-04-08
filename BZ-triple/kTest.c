#include <stdio.h>

#include "include/header.h"

double tridag(double*,double*,double*,double*,double *,int N);

int main(){
	int i;
	double a[4],b[4],c[4],r[4],u[4];
	
	a[1]=5.0;a[2]=2.0;a[3]=-2.0;
	b[0]=1.0;b[1]=3.0;b[2]=7;b[3]=8.0;
	c[0]=2.0;c[1]=8.0;c[2]=1.0;

	r[0]=-1.0;r[1]=60.0;r[2]=69.0;r[3]=14.0;
	tridag(a,b,c,r,u,4);
	for(i=0;i<4;i++) printf("%lf\n",u[i]);
	

	return 0;
}
