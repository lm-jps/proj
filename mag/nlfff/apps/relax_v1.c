/*
 * Relaxation step for NLFFF model code
 * Adapted from Wiegelmann's relax1.c
 * Takes in initial guess and model parameters
 * Computes B
 * Adapted by X. Sun
 * Version 1.0 Feb 18 2010
 * Version 1.1 Mar 17 2010
 *	Took out green(), removed B_PF
 * Version 2.0 Mar 26 2010
 *	Move boundary substitution out to main
 *	Move macros to main
 *	Work with multigrid now
 */

#ifndef PI
#include <math.h>
#define PI (M_PI)
#endif 

#ifndef KAHAN_SUM
#define KAHAN_SUM
#include "kahan_sum.c"
#endif

/* Macros */
#ifndef Macro
#define Macro
#define MCENTERGRAD(f,id) ((f[i+id]-f[i-id])/(2*h))
#define MLEFTGRAD(f,id)   ((-3*f[i]+4*f[i+id]-f[i+2*id])/(2*h))
#define MRIGHTGRAD(f,id)  ((+3*f[i]-4*f[i-id]+f[i-2*id])/(2*h))
#define GRADX(f,i) ((ix>0 && ix<nx-1) ? (MCENTERGRAD(f,nynz)) : ((ix==0) ? (MLEFTGRAD(f,nynz)) : ((ix==nx-1) ? (MRIGHTGRAD(f,nynz)) : (0.0))))
#define GRADY(f,i) ((iy>0 && iy<ny-1) ? (MCENTERGRAD(f,nz)) : ((iy==0) ? (MLEFTGRAD(f,nz)) : ((iy==ny-1) ? (MRIGHTGRAD(f,nz)) : (0.0))))
#define GRADZ(f,i) ((iz>0 && iz<nz-1) ? (MCENTERGRAD(f,1)) : ((iz==0) ? (MLEFTGRAD(f,1)) : ((iz==nz-1) ? (MRIGHTGRAD(f,1)) : (0.0))))
#endif

/*===========================================*/
/*=============== Prototypes ================*/
/*===========================================*/


/* Originally in relax.c */
double zeit();

/* Originally in optimization.c */
double calculateL(double *Bx, double *By, double *Bz, 
                  double *DivB,
                  double *oxa, double *oya, double *oza,
                  double *oxb, double *oyb, double *ozb,
                  double *oxbx, double *oxby, double *oxbz,
                  double *oxjx, double *oxjy, double *oxjz,
                  double *odotb,
                  double *wa0, double *wb0, 
                  double *wax, double *way, double *waz,
                  double *wbx, double *wby, double *wbz,
                  double *Fx, double *Fy, double *Fz,
                  int nx, int ny, int nz,
                  double dx, double dy, double dz, double *helpL_arr);
                  
double calculateLi(double *Bx, double *By, double *Bz, 
                   double *oxa, double *oya, double *oza,
                   double *oxb, double *oyb, double *ozb,
                   int nx, int ny, int nz, int nd,
                   double dx, double dy, double dz, double *helpL_arr);


/*=======================================================*/
/*===================== NLFFF Main ======================*/
/*=======================================================*/

void relax(double *Bx, double *By, double *Bz, 
           int nx, int ny, int nz, 
           int nd, double wf, double wd, int maxit, int verb)
{
    double Lx, Ly, Lz;
    double *Bx1, *By1, *Bz1, *Bx2, *By2, *Bz2;
    double *ox, *oy, *oz, *oxbx, *oxby, *oxbz, *odotb, *DivB;
    double *oxjx, *oxjy, *oxjz, *Fx, *Fy, *Fz;
    double *oxa, *oya, *oza, *oxb, *oyb, *ozb;
    double *wa0, *wb0;
    double *wax, *way, *waz, *wbx, *wby, *wbz;
    double dx, dy, dz, h;
    double mue;
    int nynz, nxnynz, nxmax, nymax, nzmax;
    int diagstep;
    
    int i, it, ix, iy, iz;
    double L, Li;
    double time1;
    double nave, Bave;
    double oldL;
    double prevL, newL, gradL;
    double *prof, *xprof, *yprof, *zprof;
    int statcount;
    
    L = oldL = 0.0;
    statcount = 0;

    nynz = ny * nz;
    nxnynz = nx * ny * nz;
    
    nxmax = nx; nymax = ny; nzmax = nz;
    nave = sqrt(1.0 * (nx - 2 * nd - 1) * (ny - 2 * nd - 1));

    Lx = 1.0 * (nx - 1) / nave;
    Ly = 1.0 * (ny - 1) / nave;
    Lz = 1.0 * (nz - 1) / nave;
    dx = dy = dz = Lx / (nx - 1);

    Bx1 = (double *) calloc(nxnynz, sizeof(double));
    By1 = (double *) calloc(nxnynz, sizeof(double));
    Bz1 = (double *) calloc(nxnynz, sizeof(double));
    Bx2 = (double *) calloc(nxnynz, sizeof(double));
    By2 = (double *) calloc(nxnynz, sizeof(double));
    Bz2 = (double *) calloc(nxnynz, sizeof(double));
    ox = (double *) calloc(nxnynz, sizeof(double));
    oy = (double *) calloc(nxnynz, sizeof(double));
    oz = (double *) calloc(nxnynz, sizeof(double));
    oxbx = (double *) calloc(nxnynz, sizeof(double));
    oxby = (double *) calloc(nxnynz, sizeof(double));
    oxbz = (double *) calloc(nxnynz, sizeof(double));
    odotb = (double *) calloc(nxnynz, sizeof(double));
    DivB = (double *) calloc(nxnynz, sizeof(double));
    oxjx = (double *) calloc(nxnynz, sizeof(double));
    oxjy = (double *) calloc(nxnynz, sizeof(double));
    oxjz = (double *) calloc(nxnynz, sizeof(double));
    Fx = (double *) calloc(nxnynz, sizeof(double));
    Fy = (double *) calloc(nxnynz, sizeof(double));
    Fz = (double *) calloc(nxnynz, sizeof(double));
    oxa = (double *) calloc(nxnynz, sizeof(double));
    oya = (double *) calloc(nxnynz, sizeof(double));
    oza = (double *) calloc(nxnynz, sizeof(double));
    oxb = (double *) calloc(nxnynz, sizeof(double));
    oyb = (double *) calloc(nxnynz, sizeof(double));
    ozb = (double *) calloc(nxnynz, sizeof(double));
    wax = (double *) calloc(nxnynz, sizeof(double));
    way = (double *) calloc(nxnynz, sizeof(double));
    waz = (double *) calloc(nxnynz, sizeof(double));
    wa0 = (double *) calloc(nxnynz, sizeof(double));
    wbx = (double *) calloc(nxnynz, sizeof(double));
    wby = (double *) calloc(nxnynz, sizeof(double));
    wbz = (double *) calloc(nxnynz, sizeof(double));
    wb0 = (double *) calloc(nxnynz, sizeof(double));
    
    double *helpL_arr = (double *) calloc(nxnynz, sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxnynz; i++) {
        ox[i] = oy[i] = oz[i] = oxbx[i] = oxby[i] = oxbz[i] = odotb[i] = DivB[i] = 0.0;
        oxjx[i] = oxjy[i] = oxjz[i] = Fx[i] = Fy[i] = Fz[i] = 0.0;
    }

    iz = 0; Bave = 0.0;
    for (ix = nd; ix < nx - nd; ix++)
    for (iy = nd; iy < ny - nd; iy++) {
        i = nynz * ix + nz * iy + iz;
        Bave = Bave + sqrt(Bx[i] * Bx[i] + By[i] * By[i] + Bz[i]*Bz[i]);
    }
    Bave = Bave / (nx - 2 * nd) / (ny - 2 * nd);

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxnynz; i++) {
        Bx[i] = Bx[i] / Bave;
        By[i] = By[i] / Bave;
        Bz[i] = Bz[i] / Bave;
    }

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxnynz; i++) {
        Bx1[i] = Bx2[i] = Bx[i];
        By1[i] = By2[i] = By[i];
        Bz1[i] = Bz2[i] = Bz[i];

    }

    /* ================= Weighting Function ================= */

    h = dx;  // FOR GRAD!

    if (verb) printf("Calculate wa and wb\n");

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxnynz; i++) {
        wa0[i] = wf; wb0[i] = wd;   // Oct 4 2020
    }

    if (nd > 1) {
        prof = (double*) calloc(nd, sizeof(double));
        xprof = (double *) calloc(nx, sizeof(double));
        yprof = (double *) calloc(ny, sizeof(double));
        zprof = (double *) calloc(nz, sizeof(double));
        for (i = 0; i < nd; i++) prof[i] = 0.5 + 0.5 * cos(1.0 * i / (nd - 1) * PI);
        for (ix = 0; ix < nx; ix++) xprof[ix] = 1.0;
        for (iy = 0; iy < ny; iy++) yprof[iy] = 1.0;
        for (iz = 0; iz < nz; iz++) zprof[iz] = 1.0;
        for (i = 0; i < nd; i++) {
            xprof[nx - nd + i] = prof[i];
            xprof[nd - i - 1] = prof[i]; 
            yprof[ny - nd + i] = prof[i];
            yprof[nd - i - 1] = prof[i];
            zprof[nz - nd + i] = prof[i];
        }

#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,iz)
#endif
        for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < ny; iy++)
        for (iz = 0; iz < nz; iz++) {
            i = ix * nynz + iy * nz + iz;
            wa0[i] = xprof[ix] * yprof[iy] * zprof[iz];
            wb0[i] = xprof[ix] * yprof[iy] * zprof[iz];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxnynz; i++) {
        wax[i] = GRADX(wa0,i);
        way[i] = GRADY(wa0,i);
        waz[i] = GRADZ(wa0,i);
        wbx[i] = GRADX(wb0,i);
        wby[i] = GRADY(wb0,i);
        wbz[i] = GRADZ(wb0,i);
    }

    /* ================= Relaxation Loop ================= */

    it = -1;
    mue = 0.1 * dx * dx;

    // Added by xudong mar 18 2010
    // When L > oldL we need to back up one step, B, B1, B2 all get restored
    // and L gets recomputed in the next step, using the restored B
    // Ideally we'll get L = oldL, but rounding errors may cause this to fail
    // Added a restore flag
    int restore = 0;

    while (it < maxit && statcount < 10 && mue > (1.0e-7 * dx * dx)) 
    {
        it++;
        memset(helpL_arr, 0.0, nxnynz * sizeof(double));
        L = calculateL(Bx, By, Bz, 
                       DivB, 
                       oxa, oya, oza, 
                       oxb, oyb, ozb, 
                       oxbx, oxby, oxbz, 
                       oxjx, oxjy, oxjz, 
                       odotb, 
                       wa0, wb0, 
                       wax, way, waz, 
                       wbx, wby, wbz, 
                       Fx, Fy, Fz, 
                       nx, ny, nz,
                       dx, dy, dz, helpL_arr);

        if (it == 0)
            oldL = L;

        if (restore) L = oldL;

        if (it > 0 && L > oldL) {
            restore = 1;
            mue = mue / 2.0;
            if (verb) {
                printf("mue reduced, mue= %lf \t mue/dx^2 = %lf\n", mue, mue / (dx * dx));
                printf("oldL= %lf \t L=%lf\n", oldL, L);
            }
            it--;
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
            for (i = 0; i < nxnynz; i++) {	// step too large, refused
                Bx[i] = Bx1[i] = Bx2[i];
                By[i] = By1[i] = By2[i];
                Bz[i] = Bz1[i] = Bz2[i];
            }
        } else {
            restore = 0;
            mue = mue * 1.01;
            oldL = L;
        }
        
        if (restore == 0) {
            
        /* L minimize */

#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,iz)
#endif
            for (ix = 1; ix < nx - 1; ix++)
            for (iy = 1; iy < ny - 1; iy++)
            for (iz = 1; iz < nz - 1; iz++) {
                i = ix * nynz + iy * nz + iz;
                Bx1[i] = Bx[i] + mue * Fx[i];
                By1[i] = By[i] + mue * Fy[i];
                Bz1[i] = Bz[i] + mue * Fz[i];
            }
        }
        
        /* Some tests with (Anti) Symmetric Boundary-conditions */
        diagstep = 10;
        if (it % diagstep == 0) /* diagnostic */
        {
            time1 = zeit();
            memset(helpL_arr, 0.0, nxnynz * sizeof(double));
            Li = calculateLi(Bx, By, Bz, 
                             oxa, oya, oza, 
                             oxb, oyb, ozb, 
                             nx, ny, nz, nd,
                             dx, dy, dz, helpL_arr);
            if (verb) printf("%i L= %.4f, L_i= %.4f,  cpu %.1f s", it, L, Li, time1);

            /* NEW: calc gradient (dL/dt)L for stopping rule */
            if (it==0) {
                prevL = 2.0 * L;
                newL = L;
            } else {
                prevL = newL;
                newL = L;
            }
            gradL = fabs((newL - prevL) / newL);
            if (verb) printf(",  gradL/L= %lf\n", gradL);
            if (gradL < 0.0001) {
                statcount++;
                if (verb) printf("*** STATIONARY STATE count: %i *** grad L/L= %lf \n", statcount, gradL);
            }
            if (gradL > 0.0001)
                statcount = 0;

        }	/* End Diagnostic */

#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
        for (i = 0; i < nxnynz; i++) {
            Bx2[i] = Bx[i];
            By2[i] = By[i];
            Bz2[i] = Bz[i];
        }
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
        for (i = 0; i < nxnynz; i++) {
            Bx[i] = Bx1[i];
            By[i] = By1[i];
            Bz[i] = Bz1[i];
        }
    }	// while

    /* Normalize B */
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
    for (i = 0; i < nxnynz; i++) {
        Bx[i] = Bx[i] * Bave;
        By[i] = By[i] * Bave;
        Bz[i] = Bz[i] * Bave;
    }

    /* Clean up */
    free(Bx1); free(By1); free(Bz1);
    free(Bx2); free(By2); free(Bz2);
    free(ox); free(oy); free(oz);
    free(oxbx); free(oxby); free(oxbz);
    free(odotb);
    free(oxjx); free(oxjy); free(oxjz);
    free(Fx); free(Fy); free(Fz);
    free(oxa); free(oya); free(oza);
    free(oxb); free(oyb); free(ozb);
    free(wa0); free(wb0);
    free(wax); free(way); free(waz);
    free(wbx); free(wby); free(wbz);
    if (nd > 1) {
        free(prof); free(xprof); free(yprof); free(zprof);
    }
    free(helpL_arr);
}

/* ========================================== */
/* ====== Routine fuer CPU Zeitmessung ====== */
/* ========================================== */

double zeit()
{
    static double tima = 0.;
    double tim, t;
    tim = (double)clock() / CLOCKS_PER_SEC;
    t = tim - tima;
    tima = tim;
    return t;
}

/* ========================================== */
/* ============== Optimization ============== */
/* ========================================== */

double calculateL(double *Bx, double *By, double *Bz, 
                  double *DivB,
                  double *oxa, double *oya, double *oza,
                  double *oxb, double *oyb, double *ozb,
                  double *oxbx, double *oxby, double *oxbz,
                  double *oxjx, double *oxjy, double *oxjz,
                  double *odotb,
                  double *wa0, double *wb0, 
                  double *wax, double *way, double *waz,
                  double *wbx, double *wby, double *wbz,
                  double *Fx, double *Fy, double *Fz,
                  int nx, int ny, int nz,
                  double dx, double dy, double dz, double *helpL_arr)
{
    double h;
    double bx, by, bz, cbx, cby, cbz, fx, fy, fz;
    double divB, b2, helpL;
    double o2a, o2b, term1x, term2x, term3x, term4x, term5ax, term5bx;
    double term1y, term2y, term3y, term4y, term5ay, term5by;
    double term1z, term2z, term3z, term4z, term5az, term5bz;
    double term6x, term6y, term6z, term7x, term7y, term7z;
    int ix, iy, iz, i;
    int nynz = ny * nz, nxnynz = nx * ny * nz;

    h = dx;
    helpL = 0.0;

#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,iz,bx,by,bz,cbx,cby,cbz,fx,fy,fz,divB,b2,o2a,o2b)
#endif
    for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++)
    for (iz = 0; iz < nz; iz++)
    {
        i = ix * nynz + iy * nz + iz;
        bx = Bx[i];
        by = By[i];
        bz = Bz[i];
        b2 = (bx * bx + by * by + bz * bz);

        cbx = GRADY(Bz,i) - GRADZ(By,i);
        cby = GRADZ(Bx,i) - GRADX(Bz,i);
        cbz = GRADX(By,i) - GRADY(Bx,i);

        fx = cby * bz - cbz * by;
        fy = cbz * bx - cbx * bz;
        fz = cbx * by - cby * bx;
        divB = GRADX(Bx,i) + GRADY(By,i) + GRADZ(Bz,i);
        DivB[i] = divB;

        if (b2 > 0.0) {
            oxa[i] = (1.0 / b2) * fx;
            oya[i] = (1.0 / b2) * fy;
            oza[i] = (1.0 / b2) * fz;
            oxb[i] = (1.0 / b2) * (divB * bx);
            oyb[i] = (1.0 / b2) * (divB * by);
            ozb[i] = (1.0 / b2) * (divB * bz);
            o2a = oxa[i] * oxa[i] + oya[i] * oya[i] + oza[i] * oza[i];
            o2b = oxb[i] * oxb[i] + oyb[i] * oyb[i] + ozb[i] * ozb[i];
            helpL_arr[i] = wa0[i] * b2 * o2a + wb0[i] * b2 * o2b;
//            helpL = helpL + wa0[i] * b2 * o2a + wb0[i] * b2 * o2b;
        }
        
        oxbx[i] = oya[i] * bz - oza[i] * by;
        oxby[i] = oza[i] * bx - oxa[i] * bz;
        oxbz[i] = oxa[i] * by - oya[i] * bx;
        odotb[i] = oxb[i] * bx + oyb[i] * by + ozb[i] * bz;
        oxjx[i] = oya[i] * cbz - oza[i] * cby;
        oxjy[i] = oza[i] * cbx - oxa[i] * cbz;
        oxjz[i] = oxa[i] * cby - oya[i] * cbx;
    }	// ix, iy, iz

    // Kahan summation
    kahan_sum(helpL_arr, nxnynz, &helpL);
    
    helpL = helpL * (dx * dy * dz);
    
    /*10.07.02 diagnostik */
#ifdef _OPENMP
#pragma omp parallel for private (i, ix, iy, iz, o2a, o2b, term1x, term2x, term3x, term4x, term5ax, term5bx, term6x, term7x, term1y, term2y, term3y, term4y, term5ay, term5by, term6y, term7y, term1z, term2z, term3z, term4z, term5az, term5bz, term6z, term7z)
#endif
    for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++)
    for (iz = 0; iz < nz; iz++)
    {
        i = ix * nynz + iy * nz + iz;
        term1x = GRADY(oxbz,i) - GRADZ(oxby,i);
        term1y = GRADZ(oxbx,i) - GRADX(oxbz,i);
        term1z = GRADX(oxby,i) - GRADY(oxbx,i);

        term2x = oxjx[i];
        term2y = oxjy[i];
        term2z = oxjz[i];

        term3x = GRADX(odotb,i);
        term3y = GRADY(odotb,i);
        term3z = GRADZ(odotb,i);

        term4x = oxb[i] * DivB[i];
        term4y = oyb[i] * DivB[i];
        term4z = ozb[i] * DivB[i];

        o2a = oxa[i] * oxa[i] + oya[i] * oya[i] + oza[i] * oza[i];
        o2b = oxb[i] * oxb[i] + oyb[i] * oyb[i] + ozb[i] * ozb[i];
        term5ax = Bx[i] * o2a;
        term5ay = By[i] * o2a;
        term5az = Bz[i] * o2a;
        term5bx = Bx[i] * o2b;
        term5by = By[i] * o2b;
        term5bz = Bz[i] * o2b;

        /* Terms regarding weighting function */
        term6x = oxby[i] * waz[i] - oxbz[i] * way[i];
        term6y = oxbz[i] * wax[i] - oxbx[i] * waz[i];
        term6z = oxbx[i] * way[i] - oxby[i] * wax[i];
        term7x = odotb[i] * wbx[i];
        term7y = odotb[i] * wby[i];
        term7z = odotb[i] * wbz[i];

        Fx[i] = wa0[i] * (term1x - term2x + term5ax) + wb0[i] * (term3x - term4x + term5bx) + term6x + term7x;
        Fy[i] = wa0[i] * (term1y - term2y + term5ay) + wb0[i] * (term3y - term4y + term5by) + term6y + term7y;
        Fz[i] = wa0[i] * (term1z - term2z + term5az) + wb0[i] * (term3z - term4z + term5bz) + term6z + term7z;
    }

    return helpL; 
}


/* ========================================== */
/* ======== Optimization (inner box) ======== */
/* ========================================== */

double calculateLi(double *Bx, double *By, double *Bz, 
                   double *oxa, double *oya, double *oza,
                   double *oxb, double *oyb, double *ozb,
                   int nx, int ny, int nz, int nd,
                   double dx, double dy, double dz, double *helpL_arr)
{
    double h;
    double bx, by, bz;
    double b2, helpL;
    double o2a, o2b, o2;
    int ix, iy, iz, i;
    int nynz = ny * nz, nxnynz = nx * ny * nz;

    h = dx;
    helpL = 0.0;

#ifdef _OPENMP
#pragma omp parallel for private (i,ix,iy,iz,bx,by,bz,b2,o2a,o2b,o2)
#endif
    for (ix = nd; ix < nx - nd; ix++)
    for (iy = nd; iy < ny - nd; iy++)
    for (iz = 0; iz < nz - nd; iz++)
    {
        i = ix * nynz + iy * nz + iz;
        bx = Bx[i];
        by = By[i];
        bz = Bz[i];
        b2 = (bx * bx + by * by + bz * bz);
        o2a = oxa[i] * oxa[i] + oya[i] * oya[i] + oza[i] * oza[i];
        o2b = oxb[i] * oxb[i] + oyb[i] * oyb[i] + ozb[i] * ozb[i];
        o2 = o2a + o2b;
        helpL_arr[i] = b2 * o2;
//        helpL = helpL + b2*o2;
    }
    
    // Kahan summation
    kahan_sum(helpL_arr, nxnynz, &helpL);
    helpL = helpL * (dx * dy * dz);
    return helpL;
}
