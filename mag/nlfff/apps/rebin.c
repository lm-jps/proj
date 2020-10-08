/*
 * Rebin magnetogram and field cube
 * Adapted from Wiegelmann's rebin.c & rebin2d.c
 * Adapted by X. Sun
 * V1 Mar 22 2010
 * V2 Feb 23 2019: Added mask array to rebin2d => rebin2d4
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
    //printf("rebin2d"); fflush(stdout);
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
    //printf("rebin2d_done"); fflush(stdout);
}


/*
 * Same as rebin 2d but with 4 arrays
 */
void rebin2d4(double *Bx, double *By, double *Bz, double *mask,
             double *Bx2, double *By2, double *Bz2, double *mask2,
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
    printf("nx2=%d, ny2=%d, level=%d, d=%d\n", nx2, ny2, level, d); fflush(stdout);
    for (iy2 = 0; iy2 < ny2; iy2++)
        for (ix2 = 0; ix2 < nx2; ix2++)
        {
            i2 = ny2 * ix2 + iy2;
            ix = ix2 * d;
            iy = iy2 * d;
            i = ny * ix + iy;
            Bx2[i2] = 0.0; By2[i2] = 0.0; Bz2[i2] = 0.0;
            mask2[i2] = 0.0;
            for (ky = 0; ky < d; ky++)
                for (kx = 0; kx < d; kx++)
                {
                    Bx2[i2] += Bx[i + kx * ny + ky];
                    By2[i2] += By[i + kx * ny + ky];
                    Bz2[i2] += Bz[i + kx * ny + ky];
                    mask2[i2] += mask[i + kx * ny + ky];
                }
            Bx2[i2] /= dd; By2[i2] /= dd; Bz2[i2] /= dd;
            mask2[i2] /= dd;
        }
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

/* New version from Wiegelmann's code rebin3d */
void rebin3d_v2(double *Bx, double *By, double *Bz,
                double *Bx2, double *By2, double *Bz2,
                int nx, int ny, int nz)
{
    int nynz = ny * nz, nxnynz = nx * ny * nz;
    int i, ix, iy, iz;
    int i1, i2, i3, i4, i5;
    int ix2, iy2, iz2;
    
    // Write values of input field to corresponding places in coarser output field
    for (ix2 = 0; ix2 < 2 * nx; ix2++)
    for (iy2 = 0; iy2 < 2 * ny; iy2++)
    for (iz2 = 0; iz2 < 2 * nz; iz2++) {
        ix = ix2 / 2;
        iy = iy2 / 2;
        iz = iz2 / 2;
        i = ix * nynz + iy * nz + iz;
        i1 = ix2 * 4 * nynz + iy2 * 2 * nz + iz2;
        if (ix2 % 2 == 0 && iy2 % 2 == 0 && iz2 % 2 == 0){
            Bx2[i1] = Bx[i];
            By2[i1] = By[i];
            Bz2[i1] = Bz[i];
        }
    }
    
    // From now on we interpolate only on the rebinned grid (smaller field is not used anymore)
    // Interpolation in x- and y-direction for each height z
    for (ix2 = 0; ix2 < (2 * nx) - 1; ix2++)
    for (iy2 = 0; iy2 < (2 * ny) - 1; iy2++)
    for (iz2 = 0; iz2 < (2 * nz) - 1; iz2++) {
        i1 = ix2 * 4 * nynz + iy2 * 2 * nz + iz2;
        // linear interpolation in x-direction for each height
        if (iy2 % 2 == 0 && ix2 % 2 != 0){
            i2 = (ix2 + 1) * 4 * nynz + iy2 * 2 * nz + iz2;
            i3 = (ix2 - 1) * 4 * nynz + iy2 * 2 * nz + iz2;
            Bx2[i1] = 0.5 * (Bx2[i2] + Bx2[i3]);
            By2[i1] = 0.5 * (By2[i2] + By2[i3]);
            Bz2[i1] = 0.5 * (Bz2[i2] + Bz2[i3]);
        }
        // linear interpolation in y-direction for each height
        else if (iy2 % 2 != 0 && ix2 % 2 == 0){
            i2 = ix2 * 4 * nynz + (iy2 + 1) * 2 * nz + iz2;
            i3 = ix2 * 4 * nynz + (iy2 - 1) * 2 * nz + iz2;
            Bx2[i1] = 0.5 * (Bx2[i2] + Bx2[i3]);
            By2[i1] = 0.5 * (By2[i2] + By2[i3]);
            Bz2[i1] = 0.5 * (Bz2[i2] + Bz2[i3]);
        }
        // bilinear interpolation in x- and y-direction wherever linear interpolation is not possible
        // due to not yet calculated values
        else if (iy2 % 2 != 0 && ix2 % 2 != 0){
            i2 = (ix2 + 1) * 4 * nynz + (iy2 - 1) * 2 * nz + iz2;
            i3 = (ix2 - 1) * 4 * nynz + (iy2 - 1) * 2 * nz + iz2;
            i4 = (ix2 + 1) * 4 * nynz + (iy2 + 1) * 2 * nz + iz2;
            i5 = (ix2 - 1) * 4 * nynz + (iy2 + 1) * 2 * nz + iz2;
            Bx2[i1] = 0.25*(Bx2[i2] + Bx2[i3] + Bx2[i4] + Bx2[i5]);
            By2[i1] = 0.25*(By2[i2] + By2[i3] + By2[i4] + By2[i5]);
            Bz2[i1] = 0.25*(Bz2[i2] + Bz2[i3] + Bz2[i4] + Bz2[i5]);
        }
    }
    // Linear interpolation in z-direction based on previously calculated (x,y)-values
    for (iz2 = 0; iz2 < (2 * nz) - 1; iz2++){
        if (iz2 % 2 != 0){
            for (iy2 = 0; iy2 < (2 * ny) - 1; iy2++)
            for (ix2 = 0; ix2 < (2 * nx) - 1; ix2++){
                i1 = ix2 * 4 * nynz + iy2 * 2 * nz + iz2;
                i2 = ix2 * 4 * nynz + iy2 * 2 * nz + (iz2 + 1);
                i3 = ix2 * 4 * nynz + iy2 * 2 * nz + (iz2 - 1);
                Bx2[i1] = 0.5 * (Bx2[i2] + Bx2[i3]);
                By2[i1] = 0.5 * (By2[i2] + By2[i3]);
                Bz2[i1] = 0.5 * (Bz2[i2] + Bz2[i3]);
            }
        }
    }
    // Dublicate the n-1 elements of the input array to get the nth element of the output
    // array because their interpolates would lie outside the input array
    for (ix2 = 0; ix2 < 2 * nx; ix2++)
    for (iy2 = 0; iy2 < 2 * ny; iy2++)
    for (iz2 = 0; iz2 < 2 * nz; iz2++){
        i1 = ix2 * 4 * nynz + iy2 * 2 * nz + iz2;
        if (ix2 == (2 * nx) - 1){
            i2 = (ix2 - 1) * 4 * nynz + iy2 * 2 * nz + iz2;
            Bx2[i1] = Bx2[i2];
            By2[i1] = By2[i2];
            Bz2[i1] = Bz2[i2];
        }
        else if (iy2 == (2 * ny) - 1){
            i2 = ix2 * 4 * nynz + (iy2 - 1) * 2 * nz + iz2;
            Bx2[i1] = Bx2[i2];
            By2[i1] = By2[i2];
            Bz2[i1] = Bz2[i2];
        }
        else if (iz2 == (2 * nz) - 1){
            i2 = ix2 * 4 * nynz + iy2 * 2 * nz + (iz2 - 1);
            Bx2[i1] = Bx2[i2];
            By2[i1] = By2[i2];
            Bz2[i1] = Bz2[i2];
        }
    }
}
