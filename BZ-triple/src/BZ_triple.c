#include <complex.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <signal.h>
#include <time.h>

#include "../../analysisPic/include/tiffFunc.h"
#include "../include/header.h"

// Here is where to store the data 
#define DataPath "./Data/"

void sigfun(int);
 
double err;

int main(){
	double r, tstart,tend,t,lambda,phi,maxInt1v,minInt1v,maxInt2v,minInt2v,maxInt3v,minInt3v,maxInt1u,minInt1u,maxInt2u,minInt2u,maxInt3u,minInt3u;
	double M,ss1,ss2,vmin1,vmin2,vmax1,vmax2,umin1,umin2,umax1,umax2,u1,v1,u2,v2,vOld1,vOld2;
	double ss3,vmin3,vmax3,umin3,umax3,u3,v3,vOld3;
	volatile double temp,h,dt,tempOld;
	double x[NX],y[NY],max1,max2,max3,vt1,vt2,vt3;
	double corefac1u,corefac1v,corefac2u,corefac2v,corefac3u,corefac3v;
	uint8 image[NX*NY];
	int i,j,count,nostop,search1,search2,search3;
	time_t starttime;	

	// activator
	double **a01 = (double **)malloc(NX*sizeof(double *));
	double **a02 = (double **)malloc(NX*sizeof(double *));
	double **a03 = (double **)malloc(NX*sizeof(double *));

	double **a1 = (double **)malloc(NX*sizeof(double *));
	double **a2 = (double **)malloc(NX*sizeof(double *));
	double **a3 = (double **)malloc(NX*sizeof(double *));

	//inhibitor
	double **c01 = (double **)malloc(NX*sizeof(double *));
	double **c02 = (double **)malloc(NX*sizeof(double *));
	double **c03 = (double **)malloc(NX*sizeof(double *));
	double **c1 = (double **)malloc(NX*sizeof(double *));
	double **c2 = (double **)malloc(NX*sizeof(double *));
	double **c3 = (double **)malloc(NX*sizeof(double *));

//	printf("Size of c: %i\n\n",sizeof(c));

// forcing
	double **IF1 = (double **)malloc(NX*sizeof(double *));
	double **IF2 = (double **)malloc(NX*sizeof(double *));
	double **IF3 = (double **)malloc(NX*sizeof(double *));

	char ac,cont;

	char command[100];
	char filenameAmp[100],filenamediffa1[100],filenamediffc1[100],filenamediffaTiff1[100],filenamediffcTiff1[100];	
	char filenamediffa2[100],filenamediffc2[100],filenamediffaTiff2[100],filenamediffcTiff2[100];	
	char filenamediffa3[100],filenamediffc3[100],filenamediffaTiff3[100],filenamediffcTiff3[100];	

	FILE *fpa1,*fpc1,*fpstart1, *fpA;
	FILE *fpa2,*fpc2,*fpstart2;
	FILE *fpa3,*fpc3,*fpstart3;

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
	sprintf(filenameAmp,"%s/%s/info.txt",DataPath,DATASET);

	if((fpA=fopen(filenameAmp,"r"))!=0){
		printf("DataSet %s exists already!!\n",DATASET);	
		printf("Do you want to override? (Y/N)");
		ac=getchar();
		getchar();
		if((ac=='Y')||(ac=='y')){
			sprintf(command,"rm -r %s/%s/Movie/*",DataPath,DATASET);
			system(command);
			sprintf(command,"rm %s/%s/point1.dat",DataPath,DATASET);
			system(command);
		}
		else{
			printf("Change the Dataset in header.h and try again!\n\n");
			exit(-1);
		}
		fclose(fpA);	
	}

	fpA=fopen(filenameAmp,"w");
	fprintf(fpA,"DataSet: %s\nNX: %i\tNY: %i\ndt: %lf\t h: %lf\nepsilon: %lf\nq: %lf\nf1: %lf\nf2: %lf\nf3: %lf\nI0: %lf\n",DATASET,NX,NY,DT,dx,EPSILON,Q,F1,F2,F3,I0);

		// this should bring h and dt to numbers that are acurrate for the computer (see Numerical recipe pg. 230)
	temp=dx;
	h=temp;
	temp=DT;
	dt=temp;

	// starting conditions
//	StartCondition=3;// 1- cross, 2- random, 3- homogeneous, 4- target,default is taking the latest file
	for(i=0;i<NX;i++){
		c01[i] = (double *)malloc(NY*sizeof(double));
		c1[i] = (double *)malloc(NY*sizeof(double));
		a01[i] = (double *)malloc(NY*sizeof(double));
		a1[i] = (double *)malloc(NY*sizeof(double));

		c02[i] = (double *)malloc(NY*sizeof(double));
		c2[i] = (double *)malloc(NY*sizeof(double));
		a02[i] = (double *)malloc(NY*sizeof(double));
		a2[i] = (double *)malloc(NY*sizeof(double));
		
		c03[i] = (double *)malloc(NY*sizeof(double));
		c3[i] = (double *)malloc(NY*sizeof(double));
		a03[i] = (double *)malloc(NY*sizeof(double));
		a3[i] = (double *)malloc(NY*sizeof(double));
		
		IF1[i] = (double *)malloc(NY*sizeof(double));
		IF2[i] = (double *)malloc(NY*sizeof(double));
		IF3[i] = (double *)malloc(NY*sizeof(double));

	};

	fpstart1 = fopen("diffStart1.txt","w");
	fpstart2 = fopen("diffStart2.txt","w");
	fpstart3 = fopen("diffStart3.txt","w");

	printf("dt/(h*h): %lf\n",dt/(h*h));

	M=NX/2;

	//calculate the steady state:
	printf("F1: %lf\tF2: %lf\tF3: %lfQ: %lf\n",F1,F2,F3,Q);
	ss1 = -0.5*(Q-1+F1)+sqrt((Q-1+F1)*(Q-1+F1)*0.25+Q*(1+F1));
	ss2 = -0.5*(Q-1+F2)+sqrt((Q-1+F2)*(Q-1+F2)*0.25+Q*(1+F2));
	ss3 = -0.5*(Q-1+F3)+sqrt((Q-1+F3)*(Q-1+F3)*0.25+Q*(1+F3));

	// calculate the local minimum and maximum
	// Note: The maximum is without forcing, the minimum is with forcing
	vmin1 = 1000;umin1 = 0;vmax1 = 0; umax1 = 0;vOld1 = 1000;
	vmin2 = 1000;umin2 = 0;vmax2 = 0; umax2 = 0;vOld2 = 1000;
	vmin3 = 1000;umin3 = 0;vmax3 = 0; umax3 = 0;vOld3 = 1000;

	u1 = Q*1.0001;
	v1 = u1*(u1-1)/F1*(Q+u1)/(Q-u1);

	u2 = Q*1.0001;
	v2 = u2*(u2-1)/F2*(Q+u2)/(Q-u2);

	u3 = Q*1.0001;
	v3 = u3*(u3-1)/F3*(Q+u3)/(Q-u3);

	search1 = 1;
	search2 = 1;
	search3 = 1;

	nostop=1;
	while (nostop){
		vOld1 = v1;
		v1 = u1*(u1-1)/F1*(Q+u1)/(Q-u1);
		vOld2 = v2;
		v2 = u2*(u2-1)/F2*(Q+u2)/(Q-u2);
		vOld3 = v3;
		v3 = u3*(u3-1)/F3*(Q+u3)/(Q-u3);
		// Cell1
		if (vmin1>1){
			if (v1>vOld1){
				vmin1=v1;
				umin1=u1;
			}
		}
		else{
			if ((vmax1==0)&&(v1<vOld1)){
				vmax1=v1;
				umax1=u1;
				search1 = 0;
			}
		}
		//Cell2
		if (vmin2>1){
			if (v2>vOld2){
				vmin2=v2;
				umin2=u2;
			}
		}
		else{
			if ((vmax2==0)&&(v2<vOld2)){
				vmax2=v2;
				umax2=u2;
				search2 = 0;
			}
		}
		//Cell3
		if (vmin3>1){
			if (v3>vOld3){
				vmin3=v3;
				umin3=u3;
			}
		}
		else{
			if ((vmax3==0)&&(v3<vOld3)){
				vmax3=v3;
				umax3=u3;
				search3 = 0;
			}
		}

		u1 = 1.01*u1;
		u2 = 1.01*u2;
		u3 = 1.01*u3;
		nostop = (search1 || search2 || search3);
//		printf("search1: %i\t search2: %i\t search3: %i\tnostop: %i\n",search1,search2,search3,nostop);
	}

	fprintf(fpA,"\nvmin1: \tvmax1: \tumin1: \tumax1: \tvmin2: \tvmax2: \tumin2: \tumax2: \tvmin3: \tvmax3: \tumin3: \tumax3: \n");
	fprintf(fpA,"\n %lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf\n",vmin1,vmax1,umin1,umax1,vmin2,vmax2,umin2,umax2,vmin3,vmax3,umin3,umax3);


	fclose(fpA);
	
	sprintf(filenamediffa1,"a-1.dat");
	sprintf(filenamediffc1,"c-1.dat");
	sprintf(filenameAmp,"/home/stwe/Data/%s/point1.dat",DATASET);
	
	sprintf(filenamediffa2,"a-2.dat");
	sprintf(filenamediffc2,"c-2.dat");

	sprintf(filenamediffa3,"a-3.dat");
	sprintf(filenamediffc3,"c-3.dat");


	vt1 = (vmin1+vmax1)/2;
	vt2 = (vmin2+vmax2)/2;
	vt3 = (vmin3+vmax3)/2;

	printf("ss1: %f umin1:%f vmin1: %f umax1: %f vmax1: %f vt1: %f\n",ss1,umin1,vmin1,umax1,vmax1,vt1);
	printf("ss2: %f umin2:%f vmin2: %f umax2: %f vmax2: %f vt2: %f\n",ss2,umin2,vmin2,umax2,vmax2,vt2);
	printf("ss3: %f umin3:%f vmin3: %f umax3: %f vmax3: %f vt3: %f\n",ss3,umin3,vmin3,umax3,vmax3,vt3);

	// start value
	switch(StartCondition){
		case 1:
			printf("Case 1\n");
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					if (i<M){
					 	a01[i][j]=5e-3;
						a02[i][j]=5e-3;
						a03[i][j]=5e-3;
					}
					else{
						 a01[i][j]=5e-1;
						 a02[i][j]=5e-1;
						 a03[i][j]=5e-1;
					}
					if(j<M){
					 	c01[i][j] = 1e-2;
	 					c02[i][j] = 1e-2;
						c03[i][j] = 1e-2;
					}
					else{
						c01[i][j] = 25.0e-2;
						c02[i][j] = 25.0e-2;
						c03[i][j] = 25.0e-2;
					}
				};
			};
		break;
		case 2:
			printf("Case 2\n");
	
			fpa1 = fopen(filenamediffa1,"r");
			fpc1 = fopen(filenamediffc1,"r");
			fpa2 = fopen(filenamediffa2,"r");
			fpc2 = fopen(filenamediffc2,"r");
			fpa3 = fopen(filenamediffa3,"r");
			fpc3 = fopen(filenamediffc3,"r");


			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					fscanf(fpa1,"%lf\n",&a01[i][j]);
					fscanf(fpc1,"%lf\n",&c01[i][j]);
					fscanf(fpa2,"%lf\n",&a02[i][j]);
					fscanf(fpc2,"%lf\n",&c02[i][j]);
					fscanf(fpa3,"%lf\n",&a03[i][j]);
					fscanf(fpc3,"%lf\n",&c03[i][j]);
				}
					fscanf(fpa1,"\n");
					fscanf(fpc1,"\n");
					fscanf(fpa2,"\n");
					fscanf(fpc2,"\n");
					fscanf(fpa3,"\n");
					fscanf(fpc3,"\n");
			}
			fclose(fpa1);
			fclose(fpc1);
			fclose(fpa2);
			fclose(fpc2);
			fclose(fpa3);
			fclose(fpc3);

			randSpiral(a01,c01,a1,c1);
			sleep(1);
			randSpiral(a02,c02,a2,c2);
			sleep(1);
			randSpiral(a03,c03,a3,c3);

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					a01[i][j] = a1[i][j];
					c01[i][j] = c1[i][j];
					a02[i][j] = a2[i][j];
					c02[i][j] = c2[i][j];
					a03[i][j] = a3[i][j];
					c03[i][j] = c3[i][j];
				}
			}	
		break;
			//a0[100][100]=10;
		case 3:
			printf("Case 2\n");
			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					a01[i][j] = 0.05;
					c01[i][j] = 0.05;
					a02[i][j] = 0.05;
					c02[i][j] = 0.05;
					a03[i][j] = 0.05;
					c03[i][j] = 0.05;
				};
			};
		break;
		case 4:
			printf("No Case 4\n");
			exit(1);
		//	spiralGen(a01,c01,a02,c02,a03,c03);
		break;
		default:// take the last file
	
			fpa1 = fopen(filenamediffa1,"r");
			fpc1 = fopen(filenamediffc1,"r");
			fpa2 = fopen(filenamediffa2,"r");
			fpc2 = fopen(filenamediffc2,"r");
			fpa3 = fopen(filenamediffa3,"r");
			fpc3 = fopen(filenamediffc3,"r");

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					fscanf(fpa1,"%lf\n",&a01[i][j]);
					fscanf(fpc1,"%lf\n",&c01[i][j]);
					fscanf(fpa2,"%lf\n",&a02[i][j]);
					fscanf(fpc2,"%lf\n",&c02[i][j]);
					fscanf(fpa3,"%lf\n",&a03[i][j]);
					fscanf(fpc3,"%lf\n",&c03[i][j]);
				}
					fscanf(fpa1,"\n");
					fscanf(fpc1,"\n");
					fscanf(fpa2,"\n");
					fscanf(fpc2,"\n");
					fscanf(fpa3,"\n");
					fscanf(fpc3,"\n");
			}
			fclose(fpa1);
			fclose(fpc1);
			fclose(fpa2);
			fclose(fpc2);
			fclose(fpa3);
			fclose(fpc3);
		break;
	}
	// saving the start condition and calculate arrays for x, y and I
	for(i=0;i<NX;i++){
		x[i]=i*h;
		y[i]=i*h;
		for(j=0;j<NY;j++) {
			fprintf(fpstart1,"%lf\n",a01[i][j]);
			fprintf(fpstart2,"%lf\n",a02[i][j]);
			fprintf(fpstart3,"%lf\n",a03[i][j]);

			IF1[i][j]=I0;
			IF2[i][j]=I0;
			IF3[i][j]=I0;
		}
	  	fprintf(fpstart1,"\n");
	  	fprintf(fpstart2,"\n");
	  	fprintf(fpstart3,"\n");
	};
	fclose(fpstart1);
	fclose(fpstart2);
	fclose(fpstart3);

	count=0;	
	err = 100;
	t=0;
	fpA = fopen(filenameAmp,"w");
	fprintf(fpA,"#time\t a1\t c1\t coreFac1u\t coreFac1v\t a2\t c2\t a3\t c3\t coreFac2u\t coreFac2v\t coreFac3u\t coreFac3v\n");
	fclose(fpA);


	starttime = time(NULL);
	while(err>LIMIT){
		if ((t>TEND) && (cont !='y')) {
			err =0;
			//printf("TEND has been reached!\n\nContinue? (y/n):");
			//cont = getchar();
			//getchar();
			//if (cont == 'y') printf("Simulations continues .....\n\n");
			//else {
			//	printf("cont: %c",cont);
			//	putchar(cont);
			//	printf("Simulations will be stopped!\n\n");
			//	err = 0;
			//	}
		}
		count++;
		// Wait until the other cell has read the data
//		while (shm[NX*NY]==0){
		//	printf("Not yet! %lf\n",shm[NX*NY]);
//			usleep(1);
//		}
	// we start with the first cell
		max1=max(c01);
		max2=max(c02);
		max3=max(c03);

		// Cell1
		for(i=0;i<NX;i++){
		   for(j=0;j<NY;j++){
				IF1[i][j] = I0*(( c02[i][j] >vt2 )*0.5 + (c03[i][j] >vt3 )*0.5); 
				IF2[i][j] = I0*(( c01[i][j] >vt1 )*0.5 + (c03[i][j] >vt3 )*0.5); 
				IF3[i][j] = I0*(( c02[i][j] >vt2 )*0.5 + (c01[i][j] >vt1 )*0.5); 
			}
		}

		af2(x,y,a01,c01,IF1,a1,NX,'1');
		cf2(x,y,a01,c01,c1,NX);
		af2(x,y,a02,c02,IF2,a2,NX,'2');
		cf2(x,y,a02,c02,c2,NX);
		af2(x,y,a03,c03,IF3,a3,NX,'3');
		cf2(x,y,a03,c03,c3,NX);
for(i=0;i<NX;i++){
		   for(j=0;j<NY;j++){
			c01[i][j] = c1[i][j];
			a01[i][j] = a1[i][j];
			if(a01[i][j]<0) a01[i][j]=0;
			if(c01[i][j]<0) c01[i][j]=0;
			if(a01[i][j]>1) a01[i][j]=1;
			if(c01[i][j]>1) c01[i][j]=1;

			c02[i][j] = c2[i][j];
			a02[i][j] = a2[i][j];
			if(a02[i][j]<0) a02[i][j]=0;
			if(c02[i][j]<0) c02[i][j]=0;
			if(a02[i][j]>1) a02[i][j]=1;
			if(c02[i][j]>1) c02[i][j]=1;

			c03[i][j] = c3[i][j];
			a03[i][j] = a3[i][j];
			if(a03[i][j]<0) a03[i][j]=0;
			if(c03[i][j]<0) c03[i][j]=0;
			if(a03[i][j]>1) a03[i][j]=1;
			if(c03[i][j]>1) c03[i][j]=1;


   		   };//for loop for j
		};
		
		t=t+dt;

		//shm[NX*NY]=0;	
		if(count%SKIP==0){
			sprintf(filenamediffaTiff1,"/home/stwe/Data/%s/Movie/Cell-u-1-%07.2f.tif",DATASET,t);
			sprintf(filenamediffaTiff2,"/home/stwe/Data/%s/Movie/Cell-u-2-%07.2f.tif",DATASET,t);
			sprintf(filenamediffaTiff3,"/home/stwe/Data/%s/Movie/Cell-u-3-%07.2f.tif",DATASET,t);

			sprintf(filenamediffcTiff1,"/home/stwe/Data/%s/Movie/Cell-v-1-%07.2f.tif",DATASET,t);
			sprintf(filenamediffcTiff2,"/home/stwe/Data/%s/Movie/Cell-v-2-%07.2f.tif",DATASET,t);
			sprintf(filenamediffcTiff3,"/home/stwe/Data/%s/Movie/Cell-v-3-%07.2f.tif",DATASET,t);

			//printf("%s\n",filenamediffaTiff1);
			printf("t: %f count: %i\n",t,count/SKIP);
			fpa1 = fopen(filenamediffa1,"w");
			fpc1 = fopen(filenamediffc1,"w");
			fpa2 = fopen(filenamediffa2,"w");
			fpc2 = fopen(filenamediffc2,"w");
			fpa3 = fopen(filenamediffa3,"w");
			fpc3 = fopen(filenamediffc3,"w");

			fpA = fopen(filenameAmp,"a");

			maxInt1v=0;minInt1v=1;
			maxInt2v=-100;minInt2v=1;
			maxInt3v=-100;minInt3v=1;

			maxInt1u=-100;minInt1u=1;
			maxInt2u=-100;minInt2u=1;
			maxInt3u=-100;minInt3u=1;

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					fprintf(fpa1,"%lf\n",a1[i][j]);
					fprintf(fpc1,"%lf\n",c1[i][j]);
					fprintf(fpa2,"%lf\n",a2[i][j]);
					fprintf(fpc2,"%lf\n",c2[i][j]);
					fprintf(fpa3,"%lf\n",a3[i][j]);
					fprintf(fpc3,"%lf\n",c3[i][j]);

					if (c1[i][j]>maxInt1v) maxInt1v=c1[i][j];
					if (c1[i][j]<minInt1v) minInt1v=c1[i][j];
					if (c2[i][j]>maxInt2v) maxInt2v=c2[i][j];
					if (c2[i][j]<minInt2v) minInt2v=c2[i][j];
					if (c3[i][j]>maxInt3v) maxInt3v=c3[i][j];
					if (c3[i][j]<minInt3v) minInt3v=c3[i][j];
                                                            
					if (a1[i][j]>maxInt1u) maxInt1u=a1[i][j];
					if (a1[i][j]<minInt1u) minInt1u=a1[i][j];
					if (a2[i][j]>maxInt2u) maxInt2u=a2[i][j];
					if (a2[i][j]<minInt2u) minInt2u=a2[i][j];
					if (a3[i][j]>maxInt3u) maxInt3u=a3[i][j];
					if (a3[i][j]<minInt3u) minInt3u=a3[i][j];

				}	
				fprintf(fpc1,"\n");
				fprintf(fpa1,"\n");	
				fprintf(fpc2,"\n");
				fprintf(fpa2,"\n");
				fprintf(fpc3,"\n");
				fprintf(fpa3,"\n");
			};
			
			corefac1v =255./maxInt1v;
			corefac2v =255./maxInt2v;
			corefac3v =255./maxInt3v;

			corefac1u =255./maxInt1u;
			corefac2u =255./maxInt2u;
			corefac3u =255./maxInt3u;

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					image[j+i*NY] = c1[i][j]*corefac1v;
				};
			};
			WriteTiffgray8(filenamediffcTiff1,ImSize,image);	
			//for(i=0;i<NX;i++){
			//	for(j=0;j<NY;j++){
			//		image[j+i*NY] = a1[i][j]*corefac1u;
			//	};
			//};
			//WriteTiffgray8(filenamediffaTiff1,ImSize,image);	

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					image[j+i*NY] = c2[i][j]*corefac2v;
			  }	
			};
			WriteTiffgray8(filenamediffcTiff2,ImSize,image);

			//for(i=0;i<NX;i++){
			//	for(j=0;j<NY;j++){
			//		image[j+i*NY] = a2[i][j]*corefac2u;
			//  }	
			//};
			//WriteTiffgray8(filenamediffaTiff2,ImSize,image);	

			for(i=0;i<NX;i++){
				for(j=0;j<NY;j++){
					image[j+i*NY] = c3[i][j]*corefac3v;
			  }	
			};
			WriteTiffgray8(filenamediffcTiff3,ImSize,image);	

			//for(i=0;i<NX;i++){
			//	for(j=0;j<NY;j++){
			//		image[j+i*NY] = a3[i][j]*corefac3u;
			//  }	
			//};
			//WriteTiffgray8(filenamediffaTiff3,ImSize,image);	
			fprintf(fpA,"%lf\t %lf\t %lf\t %lf\t %lf\t",t,a1[250][250],c1[250][250],corefac1u,corefac1v);
			fprintf(fpA,"%lf\t %lf\t %lf\t %lf\t",a2[250][250],c2[250][250],corefac2u,corefac2v);
			fprintf(fpA,"%lf\t %lf\t %lf\t %lf\n",a3[250][250],c3[250][250],corefac3u,corefac3v);


			fclose(fpa1); fclose(fpc1);
			fclose(fpa2); fclose(fpc2);
			fclose(fpa3); fclose(fpc3);


			fclose(fpA);

		};
	};//end while loop for time
	// Save the last state into the DataSet directory
	sprintf(command,"cp a-?.dat /home/stwe/Data/%s/",DATASET);
	system(command);
	sprintf(command,"cp c-?.dat /home/stwe/Data/%s/",DATASET);
	system(command);

	printf("time: %d\n",time(NULL)-starttime);

	return 0;

}//end main

void sigfun(int sig)
{
	printf("Save the state and close the programm!\n\n");
	err=0;
//	printf("Error: %lf\n",err);
}

