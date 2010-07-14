// static char rcsid[] = "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/patch/apps/interpolate.c,v 1.1 2010/07/14 09:20:07 xudong Exp $";
/*
 *  interpolate.c                                 ~soi/(version)/src/libM.d
 *
 *  Interpolation functions.
 *
 *  Contents:
 *	linint2, linint2d	bilinear
 *      ccint2, ccint2d		cubic convolution
 *
 *    The following function is declared local only:
 *      ccker
 *
 *  Responsible:  Rick Bogart			RBogart@solar.Stanford.EDU
 *
 *  Usage:
 *	See descriptions of functions
 *
 *  Bugs:
 *	See descriptions of functions
 *
 *  Planned updates:
 *	Optimization
 *
 *  Revision history is at end of file.
 */

// #include <soi_util.h>
#include <math.h>
#include "jsoc_main.h"

/*
 *  Cubic convolution interpolation, based on GONG software, received Apr-88
 *  Interpolation by cubic convolution in 2 dimensions with 32-bit float data
 *
 *  Reference:
 *    R.G. Keys in IEEE Trans. Acoustics, Speech, and Signal Processing, 
 *    Vol. 29, 1981, pp. 1153-1160.
 *
 *  Usage:
 *    float ccint2 (float *f, int nx, int ny, double x, double y)
 *    double ccint2d (double *f, int nx, int ny, double x, double y)
 *
 *  Bugs:
 *    Currently interpolation within one pixel of edge is not
 *      implemented.  If x or y is in this range or off the picture, the
 *      function value returned is NaN.
 */

				/*  Calculate the interpolation kernel.  */
static void ccker (double *u, double s) {
  double s2, s3;

  s2= s * s;
  s3= s2 * s;
  u[0] = s2 - 0.5 * (s3 + s);
  u[1] = 1.5*s3 - 2.5*s2 + 1.0;
  u[2] = -1.5*s3 + 2.0*s2 + 0.5*s;
  u[3] = 0.5 * (s3 - s2);
}

					/*  Cubic convolution interpolation  */
float ccint2 (float *f, int nx, int ny, double x, double y) {
  double  ux[4], uy[4];
  double sum;
  int ix, iy, ix1, iy1, i, j;

  if (x < 1. || x >= (float)(nx-2) || y < 1. || y >= (float)(ny-2))
    return DRMS_MISSING_FLOAT;

  ix = (int)x;
  ccker (ux,  x - (double)ix);
  iy = (int)y;
  ccker (uy,  y - (double)iy);

  ix1 = ix - 1;
  iy1 = iy - 1;
  sum = 0.;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      sum = sum + f[(iy1+i) * nx + ix1 + j] * uy[i] * ux[j];
  return (float)sum;
}

double ccint2d (double *f, int nx, int ny, double x, double y) {
  double  ux[4], uy[4];
  double sum;
  int ix, iy, ix1, iy1, i, j;

  if (x < 1. || x >= (float)(nx-2) || y < 1. || y >= (float)(ny-2))
    return DRMS_MISSING_FLOAT;

  ix = (int)x;
  ccker (ux,  x - (double)ix);
  iy = (int)y;
  ccker (uy,  y - (double)iy);

  ix1 = ix - 1;
  iy1 = iy - 1;
  sum = 0.;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      sum = sum + f[(iy1+i) * nx + ix1 + j] * uy[i] * ux[j];
  return sum;
}

/*
 *  Bilinear interpolation, based on pipeLNU remap.
 *
 *  Reference: 
 *     Abramowitz, Milton, and Stegun, Irene, "Handbook of Mathematical
 *     Functions," Dover Publications, New York, 1964.
 *
 *  Usage:
 *    float linint2 (float *f, int nx, int ny, double x, double y)
 *    double linint2d (double *f, int nx, int ny, double x, double y)
 *
 */

float linint2 (float *f, int nx, int ny, double x, double y) {

   double p0, p1, q0, q1;			/* weights		*/
   int ilow, jlow;				/* selected CCD pixels	*/
   float *fptr;					/* input array temp	*/
   double u;                        

   ilow = floor (x);
   jlow = floor (y);
   if (ilow < 0) ilow = 0;
   if (ilow+1 >= nx) ilow -= 1;
   if (jlow < 0) jlow = 0;
   if (jlow+1 >= ny) jlow -= 1;
   p1 = x - ilow;
   p0 = 1.0 - p1;
   q1 = y - jlow;
   q0 = 1.0 - q1;
			/*  Bilinear interpolation from four data points. */
   u = 0.0;
   fptr = f + jlow*nx + ilow;
   u += p0 * q0 * *fptr++;
   u += p1 * q0 * *fptr;
   fptr += nx - 1;
   u += p0 * q1 * *fptr++;
   u += p1 * q1 * *fptr;
   return (float)u;
}

double linint2d (double *f, int nx, int ny, double x, double y) {

   double p0, p1, q0, q1;			/* weights		*/
   int ilow, jlow;				/* selected CCD pixels	*/
   double *fptr;				/* input array temp	*/
   double u;                        

   ilow = floor (x);
   jlow = floor (y);
   if (ilow < 0) ilow = 0;
   if (ilow+1 >= nx) ilow -= 1;
   if (jlow < 0) jlow = 0;
   if (jlow+1 >= ny) jlow -= 1;
   p1 = x - ilow;
   p0 = 1.0 - p1;
   q1 = y - jlow;
   q0 = 1.0 - q1;
			/*  Bilinear interpolation from four data points. */
   u = 0.0;
   fptr = f + jlow*nx + ilow;
   u += p0 * q0 * *fptr++;
   u += p1 * q0 * *fptr;
   fptr += nx - 1;
   u += p0 * q1 * *fptr++;
   u += p1 * q1 * *fptr;
   return u;
}

/*
$Id: interpolate.c,v 1.1 2010/07/14 09:20:07 xudong Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/patch/apps/interpolate.c,v $
$Author: xudong $
 * $Log: interpolate.c,v $
 * Revision 1.1  2010/07/14 09:20:07  xudong
 * for patch handeling
 *
 * Revision 1.5  2003/11/06  01:14:37  rick
 * added functions ccint2d and linint2d
 *
 *
 * Revision 1.4  1998/07/07 16:56:56  giles
 * Changed return value of ccint2 from (-1) to (NaN) in the case where
 * the target location is within one pixel of the edge of the edge of
 * the array.
 *
 * Revision 1.3  1996/07/16  22:50:01  schou
 * Changed sum from float to double to speed up things.
 *
 * Revision 1.2  1995/03/14  23:40:45  kay
 * initial version - an interpolation package
 * */
