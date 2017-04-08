#include <complex.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <signal.h>
#include <time.h>

#include "/home/sweiss/analysis/include/tiffFunc.h"
//#include "/home/stevie/Work/analysisPic/include/tiffFunc.h"
#include "../include/header.h"

#define DataPath "/data/sweiss/BZ"

#define DATASETCORE "999Core1"

void sigfun(int);
 
double err;

int main(int argc, char* argv[]){
	double r, tstart,tend,t,lambda,phi,maxInt,minInt;
	double M,ss1,ss2,vmin1,vmin2,vmax1,vmax2,umin1,umin2,umax1,umax2,u1,v1,u2,v2,vOld1,vOld2;
	double csx1, csx2;
	volatile double temp,h,dt,tempOld;
	double x[NX],y[NY],intensity,max1,max2,vt1,vt2,phasediff;
	uint8 image[NX*NY];
	int i,j,count,nostop;
	time_t starttime;	

	// activator
	double **a01 = (double **)malloc(NX*sizeof(double *));
	double **a02 = (double **)malloc(NX*sizeof(double *));
	double **a1 = (double **)malloc(NX*sizeof(double *));
	double **a2 = (double **)malloc(NX*sizeof(double *));

	//inhibitor
	double **c01 = (double **)malloc(NX*sizeof(double *));
	double **c02 = (double **)malloc(NX*sizeof(double *));
	double **c1 = (double **)malloc(NX*sizeof(double *));
	double **c2 = (double **)malloc(NX*sizeof(double *));

	double **IF = (double **)malloc(NX*sizeof(double *));

	char ac;

	char command[100];
	char filenameAmp[100],filenamediffa1[100],filenamediffc1[100],filenamediffaTiff1[100];	
	char filenamediffa2[100],filenamediffc2[100],filenamediffaTiff2[100];	

	FILE *fpa1,*fpc1,*fpstart1, *fpA;
	FILE *fpa2,*fpc2,*fpstart2;


	// If you want to save as TIFF
	uint32 ImSize[2];
	ImSize[0]=NY;
	ImSize[1]=NX;
	
	// catch Ctrl+C
	(void) signal(SIGINT,sigfun);

	//create folder for dataset
	sprintf(command,"mkdir %s/%s",DataPath,DATASETCORE);
	system(command);
	sprintf(command,"mkdir %s/%s/Movie",DataPath,DATASETCORE);
	system(command);
	sprintf(filenameAmp,"%s/%s/info.txt",DataPath,DATASETCORE);

	fpA=fopen(filenameAmp,"w");
	fprintf(fpA,"DataSet: %s\nNX: %i\tNY: %i\ndt: %lf\t h: %lf\n epsilon: %lf\n q: %lf\n f1: %lf\n f2: %lf\n I0: %lf\n",DATASETCORE,NX,NY,DT,dx,EPSILON,Q,F1,F2,I0);

	fclose(fpA);
	
	sprintf(filenamediffa1,"a-1.dat");
	sprintf(filenamediffc1,"c-1.dat");
	sprintf(filenameAmp,"%s/%s/point1.dat",DataPath,DATASETCORE);
	
	sprintf(filenamediffa2,"a-2.dat");
	sprintf(filenamediffc2,"c-2.dat");

	// Shifts of both centers in x directions
	csx1 = atof(argv[1]);
	csx2 = atof(argv[2]);
	phasediff = atof(argv[3]);
	printf("csx1: %lf, csx2: %lf phasediff: %lf\n",csx1,csx2,phasediff);
	phasediff = phasediff/1.5;
	// this should bring h and dt to numbers that are acurrate for the computer (see Numerical recipe pg. 230)
	temp=dx;
	h=temp;
	temp=DT;
	dt=temp;

	// starting conditions
	for(i=0;i<NX;i++){
		c01[i] = (double *)malloc(NY*sizeof(double));
		c02[i] = (double *)malloc(NY*sizeof(double));
		c1[i] = (double *)malloc(NY*sizeof(double));
		c2[i] = (double *)malloc(NY*sizeof(double));

		a01[i] = (double *)malloc(NY*sizeof(double));
		a02[i] = (double *)malloc(NY*sizeof(double));
		a1[i] = (double *)malloc(NY*sizeof(double));
		a2[i] = (double *)malloc(NY*sizeof(double));

		IF[i] = (double *)malloc(NY*sizeof(double));

	};


//	h1=gnuplot_init();
//	h2=gnuplot_init();
	
	fpstart1 = fopen("diffStart1.txt","w");
	fpstart2 = fopen("diffStart2.txt","w");
	
	printf("dt/(h*h): %lf\n",dt/(h*h));

	M=NX/2;

	//calculate the steady state:
	printf("F1: %lf\tF2:%lf\tQ: %lf\n",F1,F2,Q);
	ss1 = -0.5*(Q-1+F1)+sqrt((Q-1+F1)*(Q-1+F1)*0.25+Q*(1+F1));
	ss2 = -0.5*(Q-1+F2)+sqrt((Q-1+F2)*(Q-1+F2)*0.25+Q*(1+F2));

	// calculate the local minimum and maximum
	vmin1 = 1000;umin1 = 0;vmax1 = 0; umax1 = 0;vOld1 = 1000;
	vmin2 = 1000;umin2 = 0;vmax2 = 0; umax2 = 0;vOld2 = 1000;
	u1 = Q*1.001;
	v1 = u1*(u1-1)/F1*(Q+u1)/(Q-u1);
	u2 = Q*1.001;
	v2 = u2*(u2-1)/F2*(Q+u2)/(Q-u2);


	nostop=1;
	while (nostop){
		vOld1 = v1;
		v1 = u1*(u1-1)/F1*(Q+u1)/(Q-u1);
		vOld2 = v2;
		v2 = u2*(u2-1)/F2*(Q+u2)/(Q-u2);
		if (u1<ss1){
			if (v1<vmin1){
				vmin1=v1;
				umin1=u1;
			}
		}
		if (u1>ss1){
			if (v1>vmax1){
				vmax1=v1;
				umax1=u1;
			}
		}

		if (u2<ss2){
			if (v2<vmin2){
				vmin2=v2;
				umin2=u2;
			}
		}
		if (u2>ss2){
			if (v2>vmax2){
				vmax2=v2;
				umax2=u2;
			}
			if (u1>ss1){
				if ((vOld2>v2)&&(vOld2>v2)) nostop=0;
			}
		}
		u1 = 1.01*u1;
		u2 = 1.01*u2;
	}	

	vt1 = (vmin1+vmax1)/2;
	vt2 = (vmin2+vmax2)/2;


	printf("ss1: %f umin1:%f vmin1: %f umax1: %f vmax1: %f vt1: %f\n",ss1,umin1,vmin1,umax1,vmax1,vt1);
	printf("ss2: %f umin2:%f vmin2: %f umax2: %f vmax2: %f vt2: %f\n",ss2,umin2,vmin2,umax2,vmax2,vt2);
	
	// start value
	printf("Calculating Start values\n");
	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			if (i<M){
			 	a01[NX-i-1][j]=5e-3;
			}
			else{
				 a01[NX-i-1][j]=5e-1;
			}
			if(j<(M+csx1)){
			 	c01[NX-i-1][j] = 1e-2;
			}
			else{
				c01[NX-i-1][j] = 25.0e-2;
			}
		};
	};
	for(i=0;i<NX;i++){
		for(j=0;j<NY;j++){
			if (i<M){
				a02[NX-i-1][j]=5e-3;
			}
			else{
				 a02[NX-i-1][j]=5e-1;
			}
			if(j<(M+csx2)){
				c02[NX-i-1][j] = 1e-2;
			}
			else{
				c02[NX-i-1][j] = 25.0e-2;
			}
		};
	};

	// saving the start condition and calculate arrays for x, y and I
	for(i=0;i<NX;i++){
		x[i]=i*h;
		y[i]=i*h;
		for(j=0;j<NY;j++) {
			fprintf(fpstart1,"%lf\n",a01[i][j]);
			fprintf(fpstart2,"%lf\n",a02[i][j]);
			IF[i][j]=0;
		}
	  	fprintf(fpstart1,"\n");
	  	fprintf(fpstart2,"\n");

	};
	fclose(fpstart1);
	fclose(fpstart2);


	count=0;	
	err = 100;
	t=0;
	fpA = fopen(filenameAmp,"w");
	fprintf(fpA,"#time\t a1\t c1\t intFac1\t a2\t c2\t intFac2\n");
	fclose(fpA);


	starttime = time(NULL);
	while((err>LIMIT)&&(t<80)){
		count++;
		max1=max(c01);
		max2=max(c02);

		//IF=0.01;
		af2(x,y,a01,c01,IF,a1,NX,'1',vmax1);
		cf2(x,y,a01,c01,c1,NX);
		intensity=0;
		for(i=0;i<NX;i++){
		   for(j=0;j<NY;j++){
			c01[i][j] = c1[i][j];
			a01[i][j] = a1[i][j];
			if(a01[i][j]<0) a01[i][j]=0;
			if(c01[i][j]<0) c01[i][j]=0;
			if(a01[i][j]>1) a01[i][j]=1;
			if(c01[i][j]>1) c01[i][j]=1;
		//	intensity=intensity+a0[i][j];
   		   };//for loop for j
		};
		//intensity = intensity/(NX*NY*1.0);
// and now the same thing for the second cell:
	
		// I need a 2d-array for the forcing to submit to the function

		af2(x,y,a02,c02,IF,a2,NX,'2',vmax2);
		cf2(x,y,a02,c02,c2,NX);
		intensity=0;
		for(i=0;i<NX;i++){
		   for(j=0;j<NY;j++){
			c02[i][j] = c2[i][j];
			a02[i][j] = a2[i][j];

			if(a02[i][j]<0) a02[i][j]=0;
			if(c02[i][j]<0) c02[i][j]=0;
			if(a02[i][j]>1) a02[i][j]=1;
			if(c02[i][j]>1) c02[i][j]=1;
		//	intensity=intensity+a0[i][j];
   		   };//for loop for j
		};
	
		t=t+dt;

		//shm[NX*NY]=0;	
		if(count%SKIPF==0){
			sprintf(filenamediffaTiff1,"%s/%s/Movie/Cell1-%07.2f.tif",DataPath,DATASETCORE,t);
			sprintf(filenamediffaTiff2,"%s/%s/Movie/Cell2-%07.2f.tif",DataPath,DATASETCORE,t);
			//printf("%s\n",filenamediffaTiff1);
			printf("t: %f count: %i\n",t,count/SKIPF);
			fpa1 = fopen(filenamediffa1,"w");
			fpc1 = fopen(filenamediffc1,"w");
			fpa2 = fopen(filenamediffa2,"w");
			fpc2 = fopen(filenamediffc2,"w");

			fpA = fopen(filenameAmp,"a");

			maxInt=0;
			minInt=1;
			for(i=0;i<NX;i++){
			for(j=0;j<NY;j++){
				fprintf(fpa1,"%lf\n",a1[i][j]);///(c[i][j]+a[i][j]));
				fprintf(fpc1,"%lf\n",c1[i][j]);///(c[i][j]+a[i][j]));
				if (c1[i][j]>maxInt) maxInt=c1[i][j];
				if (c1[i][j]<minInt) minInt=c1[i][j];

			  }	
				fprintf(fpc1,"\n");
				fprintf(fpa1,"\n");
			};

			for(i=0;i<NX;i++){
			for(j=0;j<NY;j++){
				image[j+i*NY] = (c1[i][j]-minInt)*255/(maxInt-minInt);
			};
			};
		//	corfac = maxInt;
			WriteTiffgray8(filenamediffaTiff1,ImSize,image);	
			fprintf(fpA,"%lf\t %lf\t %lf\t %lf\t",t,a1[100][100],c1[100][100],255/(maxInt-minInt));

			maxInt = 0;
			minInt = 1;
			for(i=0;i<NX;i++){
			for(j=0;j<NY;j++){
				fprintf(fpa2,"%lf\n",a2[i][j]);///(c[i][j]+a[i][j]));
				fprintf(fpc2,"%lf\n",c2[i][j]);///(c[i][j]+a[i][j]));
				if (c2[i][j]>maxInt) maxInt=c2[i][j];	
				if (c2[i][j]<minInt) minInt=c2[i][j];	

			  }	
				fprintf(fpc2,"\n");
				fprintf(fpa2,"\n");
			};

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					image[j+i*NY] = (c2[i][j]-minInt)*255/(maxInt-minInt);
			  }	
			};
			fprintf(fpA,"%lf\t %lf\t %lf\n",a2[100][100],c2[100][100],255/(maxInt-minInt));


		//	corfac = maxInt;
			WriteTiffgray8(filenamediffaTiff2,ImSize,image);	

			fclose(fpa1);
			fclose(fpc1);
			fclose(fpa2);
			fclose(fpc2);


			fclose(fpA);
		//	gnuplot_resetplot(h1);
		//	gnuplot_cmd(h1,"splot \"%s\" ",filenamediffa1);
		//	gnuplot_resetplot(h2);
		//	gnuplot_cmd(h2,"splot \"%s\" ",filenamediffc1);
		//	usleep(200000);

		};

		if ((t>=(80-phasediff-dt/2)) && (t<(80-phasediff+dt/2))){
			sprintf(command,"cp ?-1.dat %s/%s/",DataPath,DATASETCORE);
			system(command);
			printf("1st data save!\n");
			printf("t: %lf\t 3-phase.: %lf\n",t,3-phasediff-dt/2,3-phasediff+dt/2);
		}

	};//end while loop for time
	// Save the last state into the DataSet directory
	sprintf(command,"cp ?-2.dat %s/%s/",DataPath,DATASETCORE);
	system(command);
	
	printf("time: %d\n",time(NULL)-starttime);
	//gnuplot_resetplot(h2);
	//gnuplot_cmd(h2,"splot \"%s\" ",filenamediffc1);
		//sleep(1);
//	getchar();
//	}; //end for loop for e

	//gnuplot_close(h1);
	//gnuplot_close(h2);

	return 0;

}//end main

void sigfun(int sig)
{
	printf("Save the state and close the programm!\n\n");
	err=0;
//	printf("Error: %lf\n",err);
}

