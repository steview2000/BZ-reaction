#include "tiffio.h"
#define ImMaxWidth 10000 // Maximal Width of the Image
#define ImMaxHeight 10000 // Maximal Height of the Image

void CheckTiffSize(char *arg, uint32 *ImSize);
int ReadTiffgray8(char *arg, uint32 *ImSize, uint8* image);
void WriteTiffgray8(char *arg, uint32 *ImSize, uint8 *image);
double avTIFF( uint8 *image,uint32 *ImSize, int y, int x, int secHeight, int secWidth);
 
