/*
 * Green function potential field solver
 * Adapted from Wiegelmann's bfield.c
 * Takes in bz0 magnetogram
 * Computes and outputs Bpf
 * Adapted by X. Sun
 * Version 1.0 Feb 18 2010
 * Version 1.1 Mar 16 2010
 * Changed zoff Sep 25 2020
 *    Moved from relax.c
 */

#include <math.h>


/* ========================================== */
/* ========= Potential Field Solver ========= */
/* ========================================== */

void green(double *bz0, 
           double *Bx, double *By, double *Bz, 
           int nx, int ny, int nz, 
           int verb)
{
    double *x, *y, *z, *Pot0;
    double zoff, dummy1, h, r, rx, ry, rz;
    int nynz, nxnynz;
    int i, i2, ix, iy, iz, ix1, iy1;

    if (verb) {printf("Vector magnetogram loaded\n"); fflush(stdout);}

    h = 1.0;
    nynz = ny * nz;
    nxnynz = nx * ny * nz;
    zoff = 1.0 / sqrt(2 * M_PI);        // Sakurai (1982) (Eq 2.14)

    x = (double *) calloc(nx, sizeof(double));
    y = (double *) calloc(ny, sizeof(double));
    z = (double *) calloc(nz, sizeof(double));
    Pot0 = (double *) calloc(nxnynz, sizeof(double));

    for (ix = 0; ix < nx; ix++) {x[ix] = ix * 1.0;}
    for (iy = 0; iy < ny; iy++) {y[iy] = iy * 1.0;}
    for (iz = 0; iz < nz; iz++) {z[iz] = zoff + iz * 1.0;}
    
    /* Calculate Potential */

    double c, yt, t;    // kahan summation
#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,iz,dummy1,ix1,iy1,rx,ry,rz,r,i2,c,yt,t)
#endif
    for (ix = 0; ix < nx; ix++) {
        if (verb) {printf("percent finished = %lf\n", 100.0 * ix / nx); fflush(stdout);}
        for (iy = 0; iy < ny; iy++)
        for (iz = 0; iz < nz; iz++) {
            i = ix * nynz + iy*nz + iz;
            dummy1 = 0.0;
            c = 0.0;
            for (ix1 = 0; ix1 < nx; ix1++)
            for (iy1 = 0; iy1 < ny; iy1++) {
                i2 = ny * ix1 + iy1;
                rx = x[ix] - x[ix1];
                ry = y[iy] - y[iy1];
                rz = z[iz];
                r = sqrt(rx * rx + ry * ry + rz * rz);
                yt = - bz0[i2] / r - c;
                t = dummy1 + yt;
                c = (t - dummy1) - yt;
                dummy1 = t;
            }	// ix1, iy1
            Pot0[i] = dummy1 / 2.0 / M_PI;
        }	// iy, iz
    }	// ix
    
    /* B_pot */

#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,iz)
#endif
    for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++)
    for (iz = 0; iz < nz; iz++)
    {
        i = ix * nynz + iy * nz + iz;
        Bx[i] = GRADX(Pot0,i);
        By[i] = GRADY(Pot0,i);
        Bz[i] = GRADZ(Pot0,i);
    }
    
    /* Clean up */
    free(x); free(y); free(z); free(Pot0);
    
    if (verb) printf("Potential field done\n");
}
