#ifndef LUV_H
#define LUV_H

#include <cv.h>
#include <highgui.h>

/*-------Range dos valores de LUV---------*/
/*      0 < L < 100 -> range: 101         */
/*    -84 < u < 175 -> range: 260         */
/*   -125 < v < 87  -> range: 213         */
/*----------------------------------------*/

int RGBToLUV(IplImage* bitmap, int step, int nChannels);
int LUVToRGB(IplImage* bitmap, int step, int nChannels);
int RGBToXYZ(IplImage* bitmap, int step, int nChannels);
int XYZToRGB(IplImage* bitmap, int step, int nChannels);

void rgb2luv(int R, int G, int B,  int *luv);
void luv2rgb(int* vec, int* rgb);
void convertion(IplImage* bitmap, char* conv, int step, int nChannels);

#endif
