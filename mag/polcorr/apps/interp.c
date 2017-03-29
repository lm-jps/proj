/*
 * Module name:		interp.c
 *
 * Description:
 *
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 */

#include <math.h>

/* Get interpolation index of x0 in vector x */
double get_index(int sinlat, int n, double *x, double x0)
{
    int i, i0, i1;
    double ind;
    
    if (x0 <= x[0]) {
        i0 = 0; i1 = 1;
    } else if (x0 > x[n - 1]) {
        i0 = n - 2; i1 = n - 1;
    } else {
        for (int i = 0; i < n - 1; i++) {
            if (x0 > x[i] && x0 <= x[i + 1]) {
                i0 = i; i1 = i + 1;
                break;
            }
        }
    }
    /* Sine Latitude? */
    if (sinlat) {
        ind = i0 + (sin(x0) - sin(x[i0])) / (sin(x[i1]) - sin(x[i0]));
    } else {
        ind = i0 + (x0 - x[i0]) / (x[i1] - x[i0]);
    }
    return ind;
}

/* Bilinear interpolation (from obs2helio).
 * If any neighboring point is NaN the result is NaN.
 */
double linintd(double *f, int nx, int ny, double x, double y)
{
    double p0, p1, q0, q1;      // weights
    int ilow, jlow;             // selected pixels
    double *fptr;               // input array temp
    double u;
    
    ilow = (int)floor(x);
    jlow = (int)floor(y);
    if (ilow < 0) ilow = 0;
    if (ilow + 1 >= nx) ilow -= 1;
    if (jlow < 0) jlow = 0;
    if (jlow + 1 >= ny) jlow -= 1;
    p1 = x - ilow;
    p0 = 1.0 - p1;
    q1 = y - jlow;
    q0 = 1.0 - q1;
    
    /* Bilinear interpolation from four data points. */
    u = 0.0;
    fptr = f + jlow * nx + ilow;
    u += p0 * q0 * (*fptr++);
    u += p1 * q0 * (*fptr);
    fptr += nx - 1;
    u += p0 * q1 * (*fptr++);
    u += p1 * q1 * (*fptr);
    return u;
}
