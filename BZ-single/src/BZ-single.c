#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sys/types.h>
#include <signal.h>
#include <time.h>

#include "../../analysisPic/include/tiffFunc.h"

#include "../include/gnuplot_i.h"
#include "../include/header.h"

// Here is where to store the data 
#define DataPath "./Data/"

double af(double,double,double,char);
double af2(double* ,double* ,double**, double **, double**, double**, int, char );
double cf(double,double);
double cf2(double* ,double* ,double**, double **, double**,  int );
double randSpiral(double **, double **, double **,double **);
void sigfun(int);
 
double err;

int main(){
	double tstart,tend,t,lambda,phi,maxInt,minInt;
	double M,ss,vmin,vmax,umin,umax,u,v,vOld;
	volatile double temp,h,dt,tempOld;
	double x[NX],y[NY],intensity;
	double r[NX*NY],angle[NX*NY];
	uint8 image[NX*NY],imageF[NX*NY];
	int startcond,i,j,count,nostop;
	time_t starttime;	

	// activator
	double **a = (double **)malloc(NX*sizeof(double *));
	double **a0 = (double **)malloc(NX*sizeof(double *));

	//inhibitor
	double **c = (double **)malloc(NX*sizeof(double *));
	double **c0 = (double **)malloc(NX*sizeof(double *));
	
	// Forcing
	double **IF = (double **)malloc(NX*sizeof(double *));

	char ac,minIntF,maxIntF;

	char command[100];
	char filenameAmp[100],filenamediffa[100],filenamediffc[100],filenamediffaTiff[100],filenamediffFTiff[100],filenameRecord[100],filenameInfo[100];	
;	
	FILE *fpa,*fpc,*fpstart, *fpA, *fpRec, *fpI;

	//gnuplot variables
	//gnuplot_ctrl *h1;
	//gnuplot_ctrl *h2;

	// If you want to save as TIFF
	uint32 ImSize[2];
	ImSize[0]=NY;
	ImSize[1]=NX;
	
	// catch Ctrl+C
	(void) signal(SIGINT,sigfun);
	
	// Create the Data path in case it doesn't exist yet
	sprintf(command,"mkdir %s",DataPath);
	system(command);

	//create folder for dataset
	sprintf(command,"mkdir %s/%s",DataPath,DATASET);
	system(command);
	sprintf(command,"mkdir %s/%s/Movie",DataPath,DATASET);
	system(command);
	sprintf(filenameInfo,"%s/%s/info.txt",DataPath,DATASET);

	if((fpI=fopen(filenameInfo,"r"))){
		fclose(fpI);
		printf("DataSet %s exists already!!\n",DATASET);	
		printf("Do you want to override? (Y/N)");
		ac=getchar();
		printf("%c",ac);
		if((ac=='Y')||(ac=='y')){
			printf("Understood");	
			sprintf(command,"rm -r %s/%s/Movie/*",DataPath,DATASET);
			system(command);
			sprintf(command,"rm %s/%s/point1.dat",DataPath,DATASET);
			system(command);
		}
		else{
			printf("Change the Dataset in header.h and try again!\n\n");
			exit(-1);
		}
	};

	fpI = fopen(filenameInfo,"w");
	fprintf(fpI,"DataSet: %s\nNX: %i\tNY: %i\ndt: %lf\t h: %lf\nepsilon: %lf\nq: %lf\nf: %lf\nI0: %lf\nSKIP: %i\tforcing freq: %lf (Homogeneous Forcing)\nSPEED: %f\n",DATASET,NX,NY,DT,dx,EPSILON,Q,F,I0,SKIP,f_f,SPEED);


	sprintf(filenamediffa,"a-1.dat");
	sprintf(filenamediffc,"c-1.dat");
	sprintf(filenameAmp,"%s/%s/point1.dat",DataPath,DATASET);
	sprintf(filenameRecord,"%s/%s/%srecordinfo.txt",DataPath,DATASET,DATASET);

	// this should bring h and dt to numbers that are acurrate for the computer (see Numerical recipe pg. 230)
	temp=dx;
	h=temp;
	temp=DT;
	dt=temp;

	// starting conditions
	startcond= StartCondition;// 1- cross, 2- random, 3- homogeneous,default is taking the latest file
	for(i=0;i<NX;i++){
		c[i] = (double *)malloc(NY*sizeof(double));
		c0[i] = (double *)malloc(NY*sizeof(double));

		a[i] = (double *)malloc(NY*sizeof(double));
		a0[i] = (double *)malloc(NY*sizeof(double));

		IF[i] = (double *)malloc(NY*sizeof(double));

	};

	//h1=gnuplot_init();
//	h2=gnuplot_init();
	
	fpstart = fopen("diffStart1.txt","w");
	
	printf("dt/(h*h): %lf\n",dt/(h*h));

	//gnuplot_cmd(h1,"set pm3d map;set term x11");
	//gnuplot_cmd(h1,"set size square");
	//gnuplot_cmd(h1,"set title \"Cell 1\"");
//	gnuplot_cmd(h1,"set xrange [0:%i];set yrange [0:%i]",NX,NY);	
	
//	gnuplot_cmd(h2,"set pm3d map;set term x11");
//	gnuplot_cmd(h2,"set size square");
//	gnuplot_cmd(h2,"set title \"a\"");
	
	
	
	// start value
	//calculate the steady state:
	ss = -0.5*(Q-1+F)+sqrt((Q-1+F)*(Q-1+F)*0.25+Q*(1+F));
	// calculate the local minimum and maximum
	vmin = 1000;
	umin = 0;
	vmax = 0;
	umax = 0;
	vOld = 1000;
	u = Q*1.1;
	v = u*(u-1)/F*(Q+u)/(Q-u);
	nostop=1;
	while (nostop){
		vOld = v;
		v = u*(u-1)/F*(Q+u)/(Q-u);
		if (u<ss){
			if (v<vmin){
				vmin=v;
				umin=u;
			}
		}
		if (u>ss){
			if (v>vmax){
				vmax=v;
				umax=u;
			}
			if (vOld>v){
				nostop=0;
			}
		}
		u = 1.1*u;
	}	
	
	fprintf(fpI,"ss: %f\numin:%f\nvmin: %f\numax: %f\n vmax: %f\n",ss,umin,vmin,umax,vmax);
	fclose(fpI);

	printf("ss: %f umin:%f vmin: %f umax: %f vmax: %f\n",ss,umin,vmin,umax,vmax);
	M=NX/2;	
	switch(startcond){
		case 1:
	// just use the points of the minimum and the maximum, this gives a nice spiral (does not work for the excitable state)
			printf("Case 1\n");
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					if (i<(M)){
					 	a0[i][j] = umin;//2.0e-2;
					}
					else{
						 a0[i][j]= umax;//1;//5e-1;
					}
					if(j>M){
					 	c0[i][j] = vmin;//3e-2;
					}
					else{
						c0[i][j] = vmax;//1.;//0.1;
					}
				};
			};

		break;
		case 2:
			printf("Case 2\n");
		
			fpa = fopen(filenamediffa,"r");
			fpc = fopen(filenamediffc,"r");
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					fscanf(fpa,"%lf\n",&a0[i][j]);
					fscanf(fpc,"%lf\n",&c0[i][j]);
				}
					fscanf(fpa,"\n");
					fscanf(fpc,"\n");
				
			}
			fclose(fpa);
			fclose(fpc);
			randSpiral(a0,c0,a,c);
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					a0[i][j] = a[i][j];
					c0[i][j] = c[i][j];
				}
			}	
		
		break;
			//a0[100][100]=10;
		case 3:
			printf("Case 3\n");
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					a0[i][j] = 0.05;
					c0[i][j] = 0.05;
				};
			};
			//a0[100][100]=10;
			
		case 4:
			printf("Case 4\n");
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
		break;

		default:// take the last file
	
			fpa = fopen(filenamediffa,"r");
			fpc = fopen(filenamediffc,"r");
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					fscanf(fpa,"%lf\n",&a0[i][j]);
					fscanf(fpc,"%lf\n",&c0[i][j]);
				}
					fscanf(fpa,"\n");
					fscanf(fpc,"\n");
				
			}
			fclose(fpa);
			fclose(fpc);
		break;
	}

	// saving the start condition and calculate arrays for x, y and I
	for(i=0;i<NX;i++){
		if (i<NX/2){
			x[i]=(i-NX/4)*h;
		}
		else x[i]=(i-3*NX/4)*h; 

		for(j=0;j<NY;j++) {
			if (j<NY/2){
				y[j]=(j-NY/4)*h;
			}
			else y[j]=(j-3*NY/4)*h;

			fprintf(fpstart,"%lf\n",a0[i][j]);
			IF[i][j]=I0;
			r[j+i*NY] = sqrt(x[i]*x[i]+y[j]*y[j]);	
			phi = atan(y[j]/(x[i]+0.00000000001));
			if(x[i]<0){
				if(y[j] >=0) phi = phi+PI;
				if(y[j] < 0) phi = phi-PI; 
				}
			angle[j+i*NY] = phi;
		}
	  	fprintf(fpstart,"\n");
	};
	fclose(fpstart);


	//gnuplot_resetplot(h1);
	//gnuplot_cmd(h1,"splot \"diffStart1.txt\"");
	//getchar();
	count=0;	
	err = 100;
	t=0;
	fpA = fopen(filenameAmp,"w");
	fpRec = fopen(filenameRecord,"w");

	fprintf(fpA,"#Time:\ta[100,100]\tc[100,100]\tmin\tmax\n");
	fclose(fpA);
	starttime = time(NULL);
	while((err>LIMIT)&&(t<TEND)){

		count++;
		//IF=0.01;
		af2(x,y,a0,c0,IF,a,NX,'1');
		cf2(x,y,a0,c0,c,NX);
		intensity=0;
		for(i=0;i<NX;i++){
		   for(j=0;j<NY;j++){
		//	c[i][j] = c0[i][j]+dt*cf(a0[i][j],c0[i][j]);
			c0[i][j] = c[i][j];
			a0[i][j] = a[i][j];
			if(a0[i][j]<0) a0[i][j]=0;
			if(c0[i][j]<0) c0[i][j]=0;
			if(a0[i][j]>1) a0[i][j]=1;
			if(c0[i][j]>1) c0[i][j]=1;
		//	intensity=intensity+a0[i][j];
   		   };//for loop for j
		};
		intensity = intensity/(NX*NY*1.0);
		t=t+dt;

		if(count%SKIP==0){
			sprintf(filenamediffaTiff,"%s/%s/Movie/Cell-%07.2f.tif",DataPath,DATASET,t);
			if (I0 !=0) sprintf(filenamediffFTiff,"%s/%s/Movie/Forcing-%07.2f.tif",DataPath,DATASET,t);

			printf("t: %f count: %i\n",t,count/SKIP);
			fpa = fopen(filenamediffa,"w");
			fpc = fopen(filenamediffc,"w");
			fpA = fopen(filenameAmp,"a");

			maxInt = 0;
			minInt = 1;
			maxIntF = 0;
			minIntF = 1;

			for(i=0;i<NX;i++){
			for(j=0;j<NY;j++){
				fprintf(fpa,"%lf\n",a[i][j]);///(c[i][j]+a[i][j]));
				fprintf(fpc,"%lf\n",c[i][j]);///(c[i][j]+a[i][j]));
				if (c[i][j]<minInt) minInt=c[i][j];	
				if (c[i][j]>maxInt) maxInt=c[i][j];	

		//		if (IF[i][j]<minIntF) minIntF = IF[i][j];	
		//		if (IF[i][j]>maxIntF) maxIntF = IF[i][j];	
			  }	
				fprintf(fpc,"\n");
				fprintf(fpa,"\n");
			};
			printf("max: %lf\tmin: %lf\t",maxInt,minInt);
			//printf("maxF: %lf\tminF: %lf\t",maxIntF,minIntF);

			fclose(fpa);
			fclose(fpc);
				
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					image[j+i*NY] = (c[i][j]-minInt)*255/(maxInt-minInt);
					imageF[j+i*NY] = (IF[i][j])*128/I0;
				}
			}
			WriteTiffgray8(filenamediffaTiff,ImSize,image);
//			WriteTiffgray8(filenamediffFTiff,ImSize,imageF);
			fprintf(fpRec,"%lf\t%s\n",t,filenamediffaTiff);
			fprintf(fpA,"%lf\t%lf\t%lf\t%lf\t%lf\n",t,a[100][100],c[100][100],maxInt,minInt);
			fclose(fpA);
		//if(GnuOut==1){
		//	gnuplot_resetplot(h1);
			//if (NX>256) gnuplot_cmd(h1,"splot \"%s\" every 2:2 ",filenamediffc);
			//else  gnuplot_cmd(h1,"splot \"%s\"",filenamediffc);
		//	gnuplot_resetplot(h2);
		//	gnuplot_cmd(h2,"splot \"%s\" ",filenamediffc);
		//	 usleep(500000);
		//};
		};
	};//end while loop for time
		printf("time: %d\n",time(NULL)-starttime);
//		gnuplot_resetplot(h2);
//		gnuplot_cmd(h2,"splot \"%s\" ",filenamediffc);
	//	getchar();
//	}; //end for loop for e
	fclose(fpRec);
	//gnuplot_close(h1);
//	gnuplot_close(h2);

	return 0;

}//end main

void sigfun(int sig)
{
	printf("Save the state and close the programm!\n\n");
	err=0;
//	printf("Error: %lf\n",err);
}

