// removes minimal gray image and divides by the background
// image is rescaled afterward

// syntax: "./flatten <file1.tif> <file2.tif> ... 

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "dirent.h"
#include "sys/types.h"
#include "../include/tiffFunc.h"

#define BackgroundFilename "background.tif"
//#define MinGrayFilename "mingray.tif"
#define ysec 172         // that is for the relevant section of the image
#define xsec 161
#define secWidth 410
#define secHeight 410

int main(int argc, char* argv[]){
  int i,j,k;
  uint32 ImSize[2];
  uint32 ImHeight,ImWidth;
  char inname[100];
  double avb,avm;
  CheckTiffSize(BackgroundFilename, ImSize);

  ImHeight = ImSize[0];
  ImWidth = ImSize[1];

  if(xsec+secWidth>ImWidth){
  	printf("ROI to wide!\n Check the values!!\n");
	exit(-1);
  }
  if(ysec+secHeight>ImHeight){
  	printf("ROI to high!\n Check the values!!\n");
	exit(-1);
  }

  printf("Imsize: %d %d\n",ImSize[0],ImSize[1]);
  uint8 *bg = (uint8*) malloc(sizeof(uint8) *ImHeight*ImWidth);
  uint8 fimage[ImHeight*ImWidth];
  uint8 *oimage = (uint8*) malloc(sizeof(uint8) *ImHeight*ImWidth);
  float image[ImHeight*ImWidth];


  ReadTiffgray8(BackgroundFilename, ImSize, bg);
    if(argc<2) {
    printf("\n\nUsage: ./flatten <file1> <file2> ... \n");
    exit(-1);
  }

 //flat everything
 for(i=1;i<argc;i++){
    char outname[] = "F-";
    printf("Lese: %s\n",argv[i]);
    if(-1 == ReadTiffgray8(argv[i],ImSize,oimage)) {printf("Images have different sizes !!\n");exit(-1);}
    avb=avTIFF(bg, ImSize, ysec, xsec, secHeight, secWidth);
    avm=avTIFF(oimage, ImSize,ysec,xsec,secHeight,secWidth);
    printf("avb: %f\navm: %f\n",avb,avm);
    	
    for(j=0;j<(ImHeight*ImWidth);j++){
	   image[j]=((float)oimage[j])/((float)bg[j]);
	 //  if(image[j]<0) printf("problem!!\n");
	  // image[j]=image[j]/avm*60;
	   image[j]=image[j]*avb;
	   image[j] = (image[j] < 0.0) ? 0 : image[j];
  	   image[j] = (image[j] > 255.0) ? 255.0 : image[j];
	   fimage[j]=(unsigned int)image[j]; 
    }
    strcpy(inname,argv[i]);
    strcat(outname,inname); 
    printf("outname: %s\n",outname);
    WriteTiffgray8(outname,ImSize,fimage);
 }

  return 0;
}
