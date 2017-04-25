/* tiffFunc.
 * several functions for reading and writing TIFF-files (based on the libTiff - Library)
 *
 *
 * Stephan Weiss (2007)
 */

#include "tiffio.h"
#include "stdio.h"
#include "../include/tiffFunc.h"


// Returns the size of an tiff-Image
void CheckTiffSize(char *arg,uint32 *ImSize){
	TIFF* tif;
	tif = TIFFOpen(arg, "r");
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &ImSize[1]);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ImSize[0]);

	TIFFClose(tif);
	if (ImSize[1] > ImMaxWidth){
		printf("Width of the image to large !!\n");
		exit(-1);
	}
	if (ImSize[0] > ImMaxHeight){
		printf("Height of the image to large !!\n");
		exit(-1);
	}
}

// Reads TIFF images and gives them back as matrix
int ReadTiffgray8(char *arg, uint32 *ImSize, uint8 *image){
	int i,j;
	uint32 w, h;
	TIFF* tif;
	tif=TIFFOpen(arg,"r");
	if (tif) {
		size_t npixels;
		uint32 *raster;

		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		npixels = w*h;
		raster = (uint32 *) _TIFFmalloc(npixels * sizeof (uint32));
		if (raster != NULL){
	  	  if (TIFFReadRGBAImage(tif, w, h, raster, 0)){
			for(i=0;i<h;i++){
			  for(j=0; j<w; j++){
		    		image[(h-1-i)*w+j] = (uint8)raster[i*w+j];
			  }
			}
	          }
		}
	        _TIFFfree(raster);
	}
	else {
	  printf("File %s not available !\n",arg);
	  exit(-1);
	}	
	TIFFClose(tif);	
	return 1;
}

void WriteTiffgray8(char *arg, uint32 *ImSize,uint8 *image){
	int i,j;
	unsigned long		stripBytes;
	unsigned long		numStrips;
	unsigned short		bitsPerSample;
	unsigned short		samplesPerPixel;
	unsigned short		photometric;
	uint32 ImHeight,ImWidth;
	ImHeight = ImSize[0];
	ImWidth = ImSize[1];
	numStrips = ImHeight;
	stripBytes = ImWidth;

	//image=malloc(sizeof(uint32)*ImSize[0]*ImSize[1]);
	TIFF* tiffout;
	//write lines to out.tif
	tiffout=TIFFOpen(arg,"w");
	
	TIFFSetField( tiffout, TIFFTAG_ROWSPERSTRIP,  (unsigned short)1 );
	TIFFSetField( tiffout, TIFFTAG_BITSPERSAMPLE,  (unsigned short)8 );
	TIFFSetField( tiffout, TIFFTAG_SAMPLESPERPIXEL,  (unsigned short)1);

	TIFFSetField(tiffout, TIFFTAG_IMAGEWIDTH,ImWidth);
	TIFFSetField(tiffout, TIFFTAG_IMAGELENGTH,ImHeight);
	TIFFSetField(tiffout, TIFFTAG_PHOTOMETRIC, (unsigned short) 1 );
	TIFFSetField(tiffout, TIFFTAG_PLANARCONFIG, (unsigned short) 1);
	TIFFSetField( tiffout, TIFFTAG_COMPRESSION, ( unsigned short ) 1 );	// None
	
	unsigned char ima[ImHeight*ImWidth];
	for (i=0;i<(ImWidth*ImHeight);i++) ima[i] = (unsigned char)image[i];

	char * pointim;
	unsigned char* srcPtr = (unsigned char*)( ima );
	unsigned long strip;
	for ( strip = 0; strip < numStrips; strip++ )
	{
		TIFFWriteEncodedStrip( tiffout, strip, srcPtr, stripBytes );
		srcPtr += ImWidth;
	}
	TIFFClose(tiffout);

}

double avTIFF(uint8 *image,uint32 *ImSize, int y, int x, int secHeight, int secWidth){
	double av;
	int i,j;
	uint8 im[ImSize[0]][ImSize[1]];
	av=0;
	for(i=y;i <= (y+secHeight);i++){
	  for(j=x;j <= (x+secWidth);j++){
	  	av=av+(double)image[i*ImSize[0]+j]/(secWidth*secHeight);
	  }
	}
	return av;
}

