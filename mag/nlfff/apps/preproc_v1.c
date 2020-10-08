/*
 * Magnetogram preprocessing
 * Adapted from Wiegelmann's prepro.c
 * Takes in bz0 magnetogram
 * Computes and outputs processed mag
 * Adapted by X. Sun
 * Version 1.0 Mar 21 2010
 *	Parameters needs to be tweaked
 */

// mu1 = 1.0; mu2 = 1.0; mu3 = 0.001; mu4 = 0.01; mu5 = 0.0;
// maxit = 10000;

#ifndef KAHAN_SUM
#define KAHAN_SUM
#include "kahan_sum.c"
#endif

void preproc(double *Bx, double *By, double *Bz, 
             int nx, int ny, 
             double mu1, double mu2, double mu3, double mu4, double mu5, 
             int maxit, int verb)
{
    double *x, *y, *r;
    double *Bxo, *Byo, *Bzo;
    double *Bx1, *By1, *Bz1;
    double *LapBx, *LapBy, *LapBz;
    double *term4a, *term4b, *term4c;
    double Bave, dt;
    double dx, dy, dxdy;
    double Emag, ihelp;
    double dL, L, L1, L2, L12, L3, L4;
    double oldL12, oldL3, oldL4, oldL;
    double dL12, dL3, dL4;
    double term1a, term1b, term1c, term1c1, term1c2; 
    double term2a, term2b, term2c;
    double term2a1, term2a2, term2b1, term2b2, term2c1, term2c2;
    double term5;
    int nxny;
    int i, ix, iy, it;
    int i1, i2, ix1, iy1;
    int i3, i4, ix2, iy2;

    nxny = nx * ny;

    LapBx = (double *) calloc(nxny, sizeof(double));
    LapBy = (double *) calloc(nxny, sizeof(double));
    LapBz = (double *) calloc(nxny, sizeof(double));
    Bx1 = (double *) calloc(nxny, sizeof(double));
    By1 = (double *) calloc(nxny, sizeof(double));
    Bz1 = (double *) calloc(nxny, sizeof(double));
    Bxo = (double *) calloc(nxny, sizeof(double));
    Byo = (double *) calloc(nxny, sizeof(double));
    Bzo = (double *) calloc(nxny, sizeof(double));
    x = (double *) calloc(nx, sizeof(double));
    y = (double *) calloc(ny, sizeof(double));
    r = (double *) calloc(nxny, sizeof(double));
    term4a = (double *) calloc(nxny, sizeof(double));
    term4b = (double *) calloc(nxny, sizeof(double));
    term4c = (double *) calloc(nxny, sizeof(double));

    // For Kahan summation
    double *Bave_arr = (double *) calloc(nxny, sizeof(double));
    double *Emag_arr = (double *) calloc(nxny, sizeof(double));
    double *ihelp_arr = (double *) calloc(nxny, sizeof(double));
    double *L4_arr = (double *) calloc(nxny, sizeof(double));
    
    double *term1a_arr = (double *) calloc(nxny, sizeof(double));
    double *term1b_arr = (double *) calloc(nxny, sizeof(double));
    double * term1c1_arr = (double *) calloc(nxny, sizeof(double));
    double *term1c2_arr = (double *) calloc(nxny, sizeof(double));
    double *term2a1_arr = (double *) calloc(nxny, sizeof(double));
    double *term2a2_arr = (double *) calloc(nxny, sizeof(double));
    double *term2b1_arr = (double *) calloc(nxny, sizeof(double));
    double *term2b2_arr = (double *) calloc(nxny, sizeof(double));
    double *term2c1_arr = (double *) calloc(nxny, sizeof(double));
    double *term2c2_arr = (double *) calloc(nxny, sizeof(double));
    double *L3_arr = (double *) calloc(nxny, sizeof(double));

    if (verb) printf("Raw vectormagnetogram loaded \n");
    Bave = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy)
#endif
    for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
        i = ny * ix + iy;
        Bave_arr[i] = sqrt(Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i]);
    }
    kahan_sum_sgl(Bave_arr, nxny, &Bave);   // single thread to keep consistent
    Bave = Bave / nx / ny;
    printf("\n Bave= %lf",Bave);

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxny; i++) {
        Bx[i] = Bx1[i] = Bxo[i] = Bx[i] / Bave;
        By[i] = By1[i] = Byo[i] = By[i] / Bave;
        Bz[i] = Bz1[i] = Bzo[i] = Bz[i] / Bave;
    }
    dx = 1.0 / (nx - 1); dy = 1.0 / (ny - 1);
    dxdy = dx * dy;
    for (ix = 0; ix < nx; ix++) x[ix] = ix * dx;
    for (iy = 0; iy < ny; iy++) y[iy] = iy * dy;

#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy)
#endif
    for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
        i = ny * ix + iy;
        r[i] = sqrt(x[ix] * x[ix] + y[iy] * y[iy]);
    }

    /* Do the Preprocessing */
    
    dt = 0.01;
    Emag = ihelp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxny; i++) {
        Emag_arr[i] = Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i];
        ihelp_arr[i] = r[i] * (Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i]);
    }
    kahan_sum_sgl(Emag_arr, nxny, &Emag);
    kahan_sum_sgl(ihelp_arr, nxny, &ihelp);
    printf("\n Emag= %lf, ihelp= %lf",Emag,ihelp);

    /* Start Iteration */
    L = 10.0; dL=1.0; oldL=10.0;
    it = -1;
    printf("\n L12, L3, L4 are multiplied with 1.0e6");
    while(it < maxit && dL > 0.0001 && dt > 1.0e-6)
    {
        it++;
        term1a = term1b = term1c = term1c1 = term1c2 = 0.0;
        term2a = term2b = term2c = 0.0;
        term2a1 = term2b1 = term2c1 = 0.0;
        term2a2 = term2b2 = term2c2 = 0.0;
        L1 = L2 = L3 = L4 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy)
#endif
        for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < ny; iy++)
        {
            i = ny * ix + iy;
            /* Term 1 Force */
            term1a_arr[i] = Bx[i] * Bz[i];
            term1b_arr[i] = By[i] * Bz[i];
            term1c1_arr[i] = Bz[i] * Bz[i];
            term1c2_arr[i] = Bx[i] * Bx[i] + By[i] * By[i];
            /* Term 2 Torque */
            term2a1_arr[i] = x[ix] * Bz[i] * Bz[i];
            term2a2_arr[i] = x[ix] * (Bx[i] * Bx[i] + By[i] * By[i]);
            term2b1_arr[i] = y[iy] * Bz[i] * Bz[i];
            term2b2_arr[i] = y[iy] * (Bx[i] * Bx[i] + By[i] * By[i]);
            term2c1_arr[i] = y[iy] * Bx[i] * Bz[i];
            term2c2_arr[i] = x[ix] * By[i] * Bz[i];
            /* Term 3 Data */
            double term3a = (Bx[i] - Bxo[i]) * (Bx[i] - Bxo[i]) / Emag;
            double term3b = (By[i] - Byo[i]) * (By[i] - Byo[i]) / Emag;
            double term3c = (Bz[i] - Bzo[i]) * (Bz[i] - Bzo[i]) / Emag;
            L3_arr[i] = term3a + term3b + term3c;
        }
        
        kahan_sum_sgl(term1a_arr, nxny, &term1a);
        kahan_sum_sgl(term1b_arr, nxny, &term1b);
        kahan_sum_sgl(term1c1_arr, nxny, &term1c1);
        kahan_sum_sgl(term1c2_arr, nxny, &term1c2);
        kahan_sum_sgl(term2a1_arr, nxny, &term2a1);
        kahan_sum_sgl(term2a2_arr, nxny, &term2a2);
        kahan_sum_sgl(term2b1_arr, nxny, &term2b1);
        kahan_sum_sgl(term2b2_arr, nxny, &term2b2);
        kahan_sum_sgl(term2c1_arr, nxny, &term2c1);
        kahan_sum_sgl(term2c2_arr, nxny, &term2c2);
        kahan_sum_sgl(L3_arr, nxny, &L3);
        
        term1a = term1a / Emag;
        term1b = term1b / Emag;
        term1c = (term1c1 - term1c2) / Emag;
        term2a = (term2a1 - term2a2) / ihelp;
        term2b = (term2b1 - term2b2) / ihelp;
        term2c = (term2c1 - term2c2) / ihelp;
        L1 = term1a * term1a + term1b * term1b + term1c * term1c;
        L2 = term2a * term2a + term2b * term2b + term2c * term2c;
        L12 = L1 + L2;
        L3 = L3 / Emag;
        
        /* Laplace Bx,By,Bz */
        
#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,ix1,ix2,iy1,iy2,i1,i2,i3,i4)
#endif
        for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < ny; iy++)
        {
            i = ny * ix + iy;
            ix1 = ix + 1; if (ix1 == nx) ix1 = 0;
            ix2 = ix - 1; if (ix2 == -1) ix2 = nx - 1;
            iy1 = iy + 1; if (iy1 == ny) iy1 = 0;
            iy2 = iy - 1; if (iy2 == -1) iy2 = ny - 1;
            i1 = ny * ix1 + iy;
            i2 = ny * ix2 + iy;
            i3 = ny * ix + iy1;
            i4 = ny * ix + iy2;
            LapBx[i] = -4 * Bx[i] + Bx[i1] + Bx[i2] + Bx[i3] + Bx[i4];
            LapBy[i] = -4 * By[i] + By[i1] + By[i2] + By[i3] + By[i4];
            LapBz[i] = -4 * Bz[i] + Bz[i1] + Bz[i2] + Bz[i3] + Bz[i4];
            L4_arr[i] = LapBx[i] * LapBx[i] + LapBy[i] * LapBy[i] + LapBz[i] * LapBz[i];
        }
        
        kahan_sum_sgl(L4_arr, nxny, &L4);
        L4 = dxdy * L4 / Emag;
        L = mu1 * L1 + mu2 * L2 + mu3 * L3 + mu4 * L4;
        
        /* Term 4 Smoothing */
#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,ix1,ix2,iy1,iy2,i1,i2,i3,i4)
#endif
        for(ix = 0; ix < nx; ix++)
        for(iy = 0; iy < ny; iy++)
        {
            i = ny * ix + iy;
            ix1 = ix + 1; if (ix1 == nx) ix1 = 0;
            ix2 = ix - 1; if (ix2 == -1) ix2 = nx - 1;
            iy1 = iy + 1; if (iy1 == ny) iy1 = 0;
            iy2 = iy - 1; if (iy2 == -1) iy2 = ny - 1;
            i1 = ny * ix1 + iy;
            i2 = ny * ix2 + iy;
            i3 = ny * ix + iy1;
            i4 = ny * ix + iy2;
            term4a[i] = - 8 * LapBx[i] + 2 * LapBx[i1] 
                        + 2 * LapBx[i2] + 2 * LapBx[i3] + 2 * LapBx[i4];
            term4b[i] = - 8 * LapBy[i] + 2 * LapBy[i1] 
                        + 2 * LapBy[i2] + 2 * LapBy[i3] + 2 * LapBy[i4];
            term4c[i] = - 8 * LapBz[i] + 2 * LapBz[i1] 
                        + 2 * LapBz[i2] + 2 * LapBz[i3] + 2 * LapBz[i4];
        }
        if (it > 0) {
            dL12 = dL3 = dL4 = 0.0;
            if (L12 != 0.0) dL12 = fabs(L12 - oldL12) / L12;
            if (L3 != 0.0) dL3 = fabs(L3 - oldL3) / L3;
            if (L4 != 0.0) dL4 = fabs(L4 - oldL4) / L4;
            if (L != 0.0) dL = (oldL - L) / L;
        }
        if (it % 100 == 0 && verb == 2) {
            printf("%i L12= %lf, L3= %lf, L4= %lf, L=%lf, dL= %lf\n", 
                   it, 1.0e6*L12, 1.0e6*L3, 1.0e6*L4, 1.0e6*L, dL);
        }
        if (oldL < L && it > 0) {
            dt = dt / 2;
            if (verb == 2) printf("dt reduced: %i, dt=%lf, dL=%lf, oldL=%lf, L=%lf\n",
                                  it, dt, dL, 1.0e6*oldL, 1.0e6*L);
            dL = 1.0;
            oldL = L;
        } else {
            if (it > 0) dt = dt * 1.001;
            oldL12 = L12; oldL3 = L3;
            oldL4 = L4;
            oldL = L;
        }
        
        /* Iterate for new B */
#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy)
#endif
        for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < ny; iy++)
        {
            i = ny * ix + iy;
            Bx1[i] = Bx[i] + dt * 
                     (mu1 * (-2 * term1a * Bz[i] + 4 * term1c * Bx[i])
                      - mu3 * 2 * (Bx[i] - Bxo[i])
                      - mu4 * term4a[i]
                      + mu2 * 
                          (4 * term2a * x[ix] * Bx[i] + 4 * term2b * y[iy] * Bx[i] 
                           - 2 * term2c * y[iy] * Bz[i]));
            By1[i] = By[i] + dt * 
                     (mu1 * (-2 * term1b * Bz[i] + 4 * term1c * By[i])
                      - mu3 * 2 * (By[i] - Byo[i])
                      - mu4 * term4b[i]
                      + mu2 * 
                          (4 * term2a * x[ix] * By[i] + 4 * term2b * y[iy] * By[i]
                           + 2 * term2c * x[ix] * Bz[i]));
            Bz1[i] = Bz[i] + dt * (- mu3 * 2 * (Bz[i] - Bzo[i]) - mu4 * term4c[i]);
        }
        
        /* Update B */
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
        for (i = 0; i < nxny; i++) {
            Bx[i] = Bx1[i];
            By[i] = By1[i];
            Bz[i] = Bz1[i];
        }

    } /*End of Iteration*/
    
    if (verb) printf("End of Iteration\n");
    if (verb == 2) {
        printf("mu3= %lf \t mu4= %lf \t mu5= %lf\n", mu3, mu4, mu5);
        printf("%i L12= %lf, L3= %lf, L4= %lf, L=%lf, dL= %lf\n", 
        it, 1.0e6*L12, 1.0e6*L3, 1.0e6*L4, 1.0e6*L,dL);
    }
    
    /* Output data */
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxny; i++) {
        Bx[i] = Bx[i] * Bave;
        By[i] = By[i] * Bave;
        Bz[i] = Bz[i] * Bave;
    }
    
    /* Clean up */
    free(LapBx); free(LapBy); free(LapBz);
    free(Bx1); free(By1); free(Bz1);
    free(Bxo); free(Byo); free(Bzo);
    free(x); free(y); free(r);
    free(term4a); free(term4b); free(term4c);
    free(term1a_arr); free(term1b_arr);
    free(term1c1_arr); free(term1c2_arr);
    free(term2a1_arr); free(term2a2_arr);
    free(term2b1_arr); free(term2b2_arr);
    free(term2c1_arr); free(term2c2_arr);
    free(L3_arr);
    free(Bave_arr); free(Emag_arr); free(ihelp_arr); free(L4_arr);
}
