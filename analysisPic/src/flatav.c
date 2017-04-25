// averages over some images gray image
// therefore it uses all available images in a folder 


#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "dirent.h"
#include "sys/types.h"
#include "tiffFunc.h"
//#include "sys/dir.h"

#define BackgroundFilename "avImage.tif"

int main(){
   int filenumber,tifnumber,l,i,j;
   char end[3];
   char name[100];
   char *dir = "./";
   struct dirent *dp;
   DIR *dfd;
   char *filename [1000];   
   uint32 ImHeight,ImWidth;
   uint32 ImSize[2];
   char backname[100];// = (char*) malloc(30+strlen(dir));
   char minname[100];// = (char*) malloc(30+strlen(dir));   

// read directory 
   if((dfd = opendir(dir)) == NULL){
	fprintf(stderr, "dirwalk: can't open %s\n",dir);
	return;
   }    
   filenumber=0;
   tifnumber=0;
   while ((dp = readdir(dfd)) != NULL) {
	if (strcmp(dp -> d_name, ".") ==0 || strcmp(dp -> d_name, "..") == 0)
 		continue;
	if (strlen(dir)+strlen(dp->d_name)+2 > sizeof(name))
	   fprintf(stderr, "dirwalk: name %s/%s too long\n",dir,dp->d_name);
	else{
	   l=strlen(dp -> d_name);
	   sprintf(name,"%s\n", dp->d_name);
	   for(i=0;i<3;i++) end[2-i]= name[l-i-1];
//	   printf("%s\n",end);
           if((name[0]=='l')&&(end[0] == 't') && (end[1] == 'i') && (end[2] == 'f')){
 		i=0;
		filename[tifnumber] = (char*) malloc(l);
	     	strcpy(filename[tifnumber],dp->d_name);			
		tifnumber++;
	   }
	}
	filenumber++;
  }
// Check the Image Size for the first image
   CheckTiffSize(filename[0],ImSize);
 
   ImHeight = ImSize[0];
   ImWidth = ImSize[1];  
   printf("Height: %i\nWidth: %i\n",ImHeight,ImWidth);  
  
   uint8 *image = malloc(ImHeight*ImWidth*sizeof(uint8));   
   uint8 *imageBG = malloc(ImHeight*ImWidth*sizeof(uint8));   

   double bg[ImHeight*ImWidth];
 
   closedir(dfd);

  for(i = 0;i<(ImHeight*ImWidth);i++){
	 bg[i] = 0;
	 image[i] = 0;
  }
  

    
 for(l=0;l<tifnumber;l++){
    CheckTiffSize(filename[l],ImSize);
    printf("%s\n",filename[l]);
    if(ImSize[0]*ImSize[1] != ImHeight*ImWidth){
	printf("Images have different sizes !!\n");
	exit(-1);
     }
    if(ReadTiffgray8(filename[l],ImSize, image) == -1){ printf("Images have different sizes!!\n");exit(-1);}
    for(i=0;i<(ImHeight*ImWidth);i++){
	   bg[i]=bg[i] + ((double)image[i])/tifnumber;
    }

 }

  ReadTiffgray8("background.tif" ,ImSize, imageBG); //Reads background image
  
  double bgi[ImHeight*ImWidth];

  // flatten:
  for(i=0;i<(ImHeight*ImWidth);i++) bgi[i] = bg[i]/((double) imageBG[i]); 

  FILE *fout;
  fout=fopen("F-avImage.txt","w");
  for(i=0;i<ImWidth;i++){
	for(j=0;j<ImHeight;j++){
		fprintf(fout,"%lf ",bgi[i*ImHeight+j]);
	}
	fprintf(fout,"\n");
  }	  
  fclose(fout);
  
  
  return 0;
}

