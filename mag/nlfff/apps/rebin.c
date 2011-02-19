/*
 * Rebin magnetogram and field cube
 * Adapted from Wiegelmann's rebin.c & rebin2d.c
 * Adapted by X. Sun
 * Version 1.0 Mar 22 2010
 */

/* Rebin magnetogram to 2 ^ (level - 1) smaller
 * Include level = 1 case
 * Adapted from rebin2d.c
 */
void rebin2d(double *Bx, double *By, double *Bz, 
             double *Bx2, double *By2, double *Bz2, 
             int nx, int ny, int level)
{
    int i, ix, iy;
    int i2, ix2, iy2, nx2, ny2;
    int k, kx, ky;
    int d;
    double dd;
    
    d = 1;
    for (k = level; k > 1; k--) { d *= 2; }
    dd = 1.0 * d * d;
    
    nx2 = nx / d; ny2 = ny / d;
 printf("rebin2d"); fflush(stdout);
    for (iy2 = 0; iy2 < ny2; iy2++)
    for (ix2 = 0; ix2 < nx2; ix2++)
    {
        i2 = ny2 * ix2 + iy2;
        ix = ix2 * d;
        iy = iy2 * d;
        i = ny * ix + iy;
        Bx2[i2] = 0.0; By2[i2] = 0.0; Bz2[i2] = 0.0;
        for (ky = 0; ky < d; ky++)
        for (kx = 0; kx < d; kx++)
        {
            Bx2[i2] += Bx[i + kx * ny + ky];
            By2[i2] += By[i + kx * ny + ky];
            Bz2[i2] += Bz[i + kx * ny + ky];
        }
        Bx2[i2] /= dd; By2[i2] /= dd; Bz2[i2] /= dd;
    }
     printf("rebin2d_done"); fflush(stdout);
}

/* Rebin 3d field cube to twice big in every dimension
 * Adapted from rebin.c
 */
void rebin3d(double *Bx, double *By, double *Bz, 
             double *Bx2, double *By2, double *Bz2, 
             int nx, int ny, int nz)
{
    int nynz = ny * nz, nxnynz = nx * ny * nz;
    int i, ix, iy, iz;
    int i2, ix2, iy2, iz2;
    int i3, ix3, iy3, iz3;
    
    for(ix2 = 0; ix2 < 2 * nx; ix2++)
    for(iy2 = 0; iy2 < 2 * ny; iy2++)
    for(iz2 = 0; iz2 < 2 * nz; iz2++)
    {
        ix = ix2 / 2;
        iy = iy2 / 2;
        iz = iz2 / 2;
        i = ix * nynz + iy * nz + iz;
        i2 = ix2 * 4 * nynz + iy2 * 2 * nz + iz2;

        ix3 = ix2 % 2;	/* For Trilinear interpolation */
        iy3 = iy2 % 2;
        iz3 = iz2 % 2;
        i3 = i;
        if (ix < nx - 1) { i3 = i3 + ix3 * nynz; }
        if (iy < ny - 1) { i3 = i3 + iy3 * nz; }
        if (iz < nz - 1) { i3 = i3 + iz3; }
        Bx2[i2] = 0.5 * (Bx[i] + Bx[i3]);
        By2[i2] = 0.5 * (By[i] + By[i3]);
        Bz2[i2] = 0.5 * (Bz[i] + Bz[i3]);
    }
}
