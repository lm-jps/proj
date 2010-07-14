/*
 * Name:		hccsss_pkg.c
 *
 * Description:		This is the C package of HCCSSS model.
 *			Including harmonic expansion, and field vector computation
 *                      Original code of the main algorithm is in IDL by Xuepu Zhao.
 *			Restructured and optimized by Xudong Sun.
 *			Refer to pfssnotes.pdf and Zhao et al (1995) for details.
 *
 * Function List:
 *			void hchf(float r, double *rrr, double *D_rrr, int lmax, float apar)
 *			void cshf(float r, double *rrr, double *D_rrr, int lc, float apar)
 *			void gh_hccsss(float *g, float *h, struct Grid *grid, int lmax, 
 *				float *gc, float *hc, int lc, float apar)
 *			void hccsss0_hc(float *g, float *h, struct Point *pt, 
 *				float *Bvec, int lmax, float apar)
 *			void hccsss0_cs(float *g, float *h, struct Point *pt, 
 *				float *Bvec, int lmax, float apar)
 *			void hccsss(float *g, float *h, float *gc, float *hc, struct Grid *grid, 
 *				float *Br, float *Bt, float *Bp, int lmax, int lc, int apar)
 *
 * Calling:		pfss_pkg.c, fieldline_pkg.c
 * Called by:		jgh2hccsss.c
 *
 * Original Source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu) 
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Nov 05 2009
 *			v1.0a		Nov 10 2009
 *			v1.0b		Dec 07 2009
 *
 * Issues:
 *			v1.0
 *			Added radial component function for HCCSSS model
 *			Tested and Added coefficient computation above cusp surface, gh_hccsss()
 *			Used LAPACK for matrix inversition instead of old self-written code as 
 *			there was some strange memory bug that wouldn't work on 64bit machine
 *			Added hccsss0_hc() and hccsss0_cs()
 *			Tested hccsss() and worked fine if placed after trace() in fieldline_pkg.c
 *			v1.0a
 *			free components of gridt
 *			v1.0b
 *			fixed a bug in hccsss where gridt->rr was not defined (!)
 *
 */
 

#include "hccsss_pkg.h"


/* ################ Radial functions for pfss_pkg.h ################ */


/* Radial part of HCCSSS below cusp surface */
void hchf(float r, double *rrr, double *D_rrr, int lmax, float apar)
{
    int Dptr;
    double R02, r2;
    double aR0, aR0l;
    double ar, arl, arl1;

    R02 = R0 * R0; r2 = r * r;
    aR0 = R0 + apar; aR0l = 1.0;
    ar = r + apar; arl = 1.0; arl1 = ar;
    for (int l = 0; l <= lmax; l++) {
        Dptr = l;
        // Need to check with Xuepu
        rrr[Dptr] = R02 * aR0l / ((l + 1.0) * arl1);
        // eta included
        D_rrr[Dptr] = - R02 * aR0l / (r2 * arl);
        aR0l *= aR0; arl = arl1; arl1 *= ar;
    }
}


/* Radial part of HCCSSS above cusp surface */
void cshf(float r, double *rrr, double *D_rrr, int lc, float apar)
{
    int Dptr;
    double RC2, r2;
    double aRC, aRC2, aRCl, aRCl1, aRC2l1;
    double aRS, aRS2, aRS2l1;
    double ar, ar2, arl, arl1, ar2l1;

    RC2 = RCP * RCP; r2 = r * r;
    aRC = RCP + apar; aRCl = 1.0; aRCl1 = aRC;
    aRS = RSS + apar; aRS2 = aRS * aRS; aRS2l1 = aRS;
    ar = r + apar; arl = 1.0; arl1 = ar;
    for (int l = 0; l <= lc; l++) {
        Dptr = l;
        double denom = ((l + 1.0) * aRS2l1 + l * aRC2l1) / (aRCl * aRS2l1);
        // Need to check with Xuepu
        rrr[Dptr] = RC2 * (1.0 / arl1 - arl / aRS2l1) / denom;
        // eta included
        D_rrr[Dptr] = - (RC2 / r2) * 
                ((l + 1.0) / arl + l * arl1 / aRS2l1) / denom;
        aRCl = aRCl1; aRCl1 *= aRC;
        aRS2l1 *= aRS2;
        arl = arl1; arl1 *= ar;
    }
}



/* ################ Harmonic coefficients above RCP ################ */


/* Computing harmonic coefficients above cusp surface
 * Using Schatten's least square technique, see X. Zhao
 * et al. (1995). Original code in IDL by X. Zhao,
 * Adapted into C with matrix inversion from LAPACK library
 */

void gh_hccsss(float *g, float *h, struct Grid *grid, int lmax, 
        float *gc, float *hc, int lc, float apar)
{
    /* Comput field values on cusp surface */
    struct Grid *gridt;
    int nt, np, npnt;
    float *Brcp, *Btcp, *Bpcp;
    void (*rfunc)(float, double *, double *, int, float);

    nt = grid->nt; np = grid->np; npnt = np * nt;
    gridt = (struct Grid *)(malloc(sizeof(struct Grid)));
    gridt->nr = 1; gridt->nt = nt; gridt->np = np;
    gridt->rr = (float *)(malloc(gridt->nr * sizeof(float)));
    gridt->th = (float *)(malloc(gridt->nt * sizeof(float)));
    gridt->ph = (float *)(malloc(gridt->np * sizeof(float)));
    gridt->rr[0] = RCP;
    for (int j = 0; j < nt; j++) gridt->th[j] = grid->th[j];
    for (int i = 0; i < np; i++) gridt->ph[i] = grid->ph[i];

    Brcp = (float *)(malloc(npnt * sizeof(float)));
    Btcp = (float *)(malloc(npnt * sizeof(float)));
    Bpcp = (float *)(malloc(npnt * sizeof(float)));

    rfunc = hchf;

    Bcube(g, h, gridt, Brcp, Btcp, Bpcp, lmax, apar, rfunc);
//    printf("\n%f, %f, %f\n", Brcp[0], Brcp[1], Brcp[2]);

    /* Functions needed for matrix ab */
    int lc1 = lc + 1;
    int lc12 = lc1 * (lc1 + 1) / 2;
    int nplc1 = np * lc1;
    int ntlc12 = nt * lc12;

    // Legendre
    double *leg = (double *)(calloc(ntlc12, sizeof(double)));
    double *D_leg = (double *)(calloc(ntlc12, sizeof(double)));
    for (int j = 0; j < nt; j++)
        pdpth(gridt->th[j], (leg + j * lc12), (D_leg + j * lc12), lc);

    // Cos(m*phi), Sin(m*phi)
    double *cmph = (double *)(malloc(nplc1 * sizeof(double)));
    double *smph = (double *)(malloc(nplc1 * sizeof(double)));
    for (int i = 0; i < np; i++) {
        for (int m = 0; m <= lc; m++) {
            cmph[i * lc1 + m] = cos(m * gridt->ph[i]);
            smph[i * lc1 + m] = sin(m * gridt->ph[i]);
        }
    }

    // Sin(theta)
    double *sth = (double *)(malloc(nt * sizeof(double)));
    for (int j = 0; j < nt; j++)
        sth[j] = sin(gridt->th[j]);

    // ar1, ar2, ar3
    double *ar1 = (double *)(malloc(lc1 * sizeof(double)));
    double *ar2 = (double *)(malloc(lc1 * sizeof(double)));
    double *ar3 = (double *)(malloc(lc1 * lc1 * sizeof(double)));

    double ars = apar + RSS, arc = apar + RCP;
    double ars2 = ars * ars, arc2 = arc * arc;
    double ars2l1 = ars, arc2l1 = arc;
    double Kl;

    for (int l = 0; l <= lc; l++) {
        Kl = (RCP / arc) * (ars2l1 - arc2l1) / 
                ((l + 1) * ars2l1 + l * arc2l1);
        ar1[l] = 1.0;
        ar2[l] = -1.0 * Kl;
        for (int m = 0; m <= l; m++)
            ar3[l * lc1 + m] = m * Kl;
        ars2l1 *= ars2; arc2l1 *= arc2;
    }

    /* Generate matrix ab */
    int ixx = lc1 * lc1, iyy = 3 * npnt;
    double *ab = (double *)(malloc(ixx * iyy * sizeof(double)));

    int lm, ptr0, ptr, ptrt0, ptrt, ptrp;
    // This part is consistent with Xuepu's IDL code,
    // not with his paper, but essentially the same
    ptr = 0;
    for (int l = 0; l <= lc; l++) {
    for (int m = 0; m <= l; m++) {
        lm = l * lc1 + m;
        ptrt0 = (l + 1) * l / 2 + m;
        for (int j = 0; j < nt; j++) {
        ptrt = j * lc12 + ptrt0;
        for (int i = 0; i < np; i++) {
            ptrp = i * lc1 + m;
            ab[ptr++] = ar1[l] * cmph[ptrp] * leg[ptrt];
            ab[ptr++] = ar2[l] * cmph[ptrp] * D_leg[ptrt];
            ab[ptr++] = ar3[lm] * smph[ptrp] * leg[ptrt] / sth[j];
        }}
    }}

    for (int l = 1; l <= lc; l++) {
    for (int m = 1; m <= l; m++) {
        lm = l * lc1 + m;
        ptrt0 = (l + 1) * l / 2 + m;
        for (int j = 0; j < nt; j++) {
        ptrt = j * lc12 + ptrt0;
        for (int i = 0; i < np; i++) {
            ptrp = i * lc1 + m;
            ab[ptr++] = ar1[l] * smph[ptrp] * leg[ptrt];
            ab[ptr++] = ar2[l] * smph[ptrp] * D_leg[ptrt];
            ab[ptr++] = - ar3[lm] * cmph[ptrp] * leg[ptrt] / sth[j];
        }}
    }}

    /* Generate ab * T(ab) */
    double *abtab = (double *)(calloc(ixx * ixx, sizeof(double)));
    for (int ia = 0; ia < ixx; ia++) {
    for (int ib = 0; ib < ixx; ib++) {
        for (int ix = 0; ix < iyy; ix++)
            abtab[ia * ixx + ib] += ab[ia * iyy + ix] * ab[ib * iyy + ix];
    }}

    /* Invert ab * T(ab) */
    double *abtab1 = (double *)(calloc(ixx * ixx, sizeof(double)));
    for (int i = 0; i < ixx * ixx; i++) abtab1[i] = abtab[i];

    int info;
    int *ipvt = (int *)(malloc(ixx * sizeof(int)));
    dgetrf_(&ixx, &ixx, abtab1, &ixx, ipvt, &info);

/*
    for (int i = 0; i < ixx; i++) {
      for (int j = 0; j < ixx; j++) {
        printf("%10.4f", abtab1[i * ixx + j]);
      }
      printf("\n");
    }
    printf("\n");
*/

    int lwork = 200 * ixx;
    double *work = (double *)(malloc(lwork * ixx * sizeof(double)));
    dgetri_(&ixx, abtab1, &ixx, ipvt, work, &lwork, &info);

/*
    for (int i = 0; i < ixx; i++) {
      for (int j = 0; j < ixx; j++) {
        printf("%10.4f", abtab1[i * ixx + j]);
      }
      printf("\n");
    }
*/

    /* Create B */
    double *Bout = (double *)(malloc(iyy * sizeof(double)));
    double sign;
    ptr = 0;
    for (int j = 0; j < nt; j++) {
    for (int i = 0; i < np; i++) {
        ptr0 = j * np + i;
        sign = (Brcp[ptr0] < 0) ? -1.0 : 1.0;
        Bout[ptr++] = Brcp[ptr0] * sign;
        Bout[ptr++] = Btcp[ptr0] * sign;
        Bout[ptr++] = Bpcp[ptr0] * sign;
    }}

    /* Multiply ab and B */
    double *abB = (double *)(calloc(ixx, sizeof(double)));
    for (int ix = 0; ix < ixx; ix++) {
    for (int iy = 0; iy < iyy; iy++) {
        abB[ix] += ab[ix * iyy + iy] * Bout[iy];
    }}

    /* Multiply abtab and abB for ghout */
    double *ghout = (double *)(calloc(ixx, sizeof(double)));
    for (int i = 0; i < ixx; i++) {
    for (int ix = 0; ix < ixx; ix++) {
        ghout[i] += abtab1[i * ixx + ix] * abB[ix];
    }}

    /* Output */
    ptr = 0;
    for (int l = 0; l <= lc; l++) {
    for (int m = 0; m <= l; m++) {
        gc[l * lc1 + m] = ghout[ptr++];
    }}
    for (int l = 1; l <= lc; l++) {
    for (int m = 1; m <= l; m++) {
        hc[l * lc1 + m] = ghout[ptr++];
    }}
    for (int l = 0; l <= lc; l++) {
        for (int m = l + 1; m <= lc; m++) {
            gc[l * lc1 + m] = NAN;
            hc[l * lc1 + m] = NAN;
        }
        hc[l * lc1] = 0.0;
    }

    // Clean up
    free(gridt->rr); free(gridt->th); free(gridt->ph); free(gridt);
    free(Brcp); free(Btcp); free(Bpcp);
    free(leg); free(D_leg);
    free(cmph); free(smph); free(sth);
    free(ar1); free(ar2); free(ar3);
    free(ab); free(abtab); free(abtab1);
    free(ipvt); free(work);
    free(Bout); free(abB); free(ghout);
}



/* ################ HCCSSS Field Vector ################ */


/* Assemble gh into B vectors on a single point, below cusp surface */
void hccsss0_hc(float *g, float *h, struct Point *pt, 
        float *Bvec, int lmax, float apar)
{
    void (*rfunc)(float, double *, double *, int, float);
    rfunc = hchf;
    Bpoint(g, h, pt, Bvec, lmax, apar, rfunc);
}


/* Assemble gh into B vectors on a single point, above cusp surface */
void hccsss0_cs(float *g, float *h, struct Point *pt, 
        float *Bvec, int lmax, float apar)
{
    void (*rfunc)(float, double *, double *, int, float);
    rfunc = cshf;
    Bpoint(g, h, pt, Bvec, lmax, apar, rfunc);
}


/* Assemble gh into B vectors for HCCSSS on a 3d grid */
void hccsss(float *g, float *h, float *gc, float *hc, struct Grid *grid, 
        float *Br, float *Bt, float *Bp, int lmax, int lc, float apar)
{
    int nt = grid->nt, np = grid->np;
    int npnt = np * nt;

    int nr = grid->nr;
    int nr0, nr1;	// Number of points below/above cusp surface
    for (nr0 = 0; nr0 < nr; nr0++)
        if (grid->rr[nr0] > RCP) break;
    nr1 = nr - nr0;

    void (*rfunc)(float, double *, double *, int, float);
    struct Grid *gridt = (struct Grid *)(malloc(sizeof(struct Grid)));
    gridt->nt = nt; gridt->np = np;
    gridt->th = (float *)(malloc(nt * sizeof(float)));
    gridt->ph = (float *)(malloc(np * sizeof(float)));
    for (int j = 0; j < nt; j++) gridt->th[j] = grid->th[j];
    for (int i = 0; i < np; i++) gridt->ph[i] = grid->ph[i];

    // Below cusp surface
    if (nr0 > 0) {
        rfunc = hchf;
        gridt->nr = nr0;
        gridt->rr = (float *)(malloc(nr0 * sizeof(float)));
        for (int n = 0; n < nr0; n++) gridt->rr[n] = grid->rr[n];
        Bcube(g, h, gridt, Br, Bt, Bp, lmax, apar, rfunc);
        free(gridt->rr);
    }

    // Above cusp surface
    if (nr1 > 0) {
        rfunc = cshf;
        gridt->nr = nr1;
        gridt->rr = (float *)(malloc(nr1 * sizeof(float)));
        for (int n = 0; n < nr1; n++) gridt->rr[n] = grid->rr[n + nr0];
        int Dpt = nr0 * npnt;
        float *Brt = Br + Dpt, *Btt = Bt + Dpt, *Bpt = Bp + Dpt;
        Bcube(gc, hc, gridt, Brt, Btt, Bpt, lc, apar, rfunc);
        // Trace down to Rcp for each point to determine the sign
        void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
							void (*)(float *, float *, struct Point *, float *, int, float));
        void (*Bfunc)(float *, float *, struct Point *, float *, int, float);
		integrator = Euler;
        Bfunc = hccsss0_cs;
        struct Point *pt0 = (struct Point *)(malloc(sizeof(struct Point)));
        struct Point *pt = (struct Point *)(malloc(sizeof(struct Point)));
        pt->r = RCP;
        struct FldLn *fl = (struct FldLn *)(malloc(sizeof(struct FldLn)));
        float Bvec[3];
        for (int n = 0; n < nr1; n++) {
            pt0->r = gridt->rr[n];
        for (int j = 0; j < nt; j++) {
            pt0->t = gridt->th[j];
        for (int i = 0; i < np; i++) {
            pt0->p = gridt->ph[i];
            trace(gc, hc, pt0, fl, lc, integrator, Bfunc, gridt->rr[n], RCP, apar);
            if (fl->num) {
                pt->t = fl->tt[fl->num];
                pt->p = fl->pp[fl->num];
                hccsss0_hc(g, h, pt, Bvec, lmax, apar);
                if (Bvec[0] < 0) {	// Reverse
                    Dpt = n * npnt + j * np + i;
                    Brt[Dpt] *= -1; Btt[Dpt] *= -1; Bpt[Dpt] *= -1;
                }
            }
        }}}
        free(pt); free(pt0); free(fl);
        free(gridt->rr);
    }

    free(gridt->th); free(gridt->ph); free(gridt);
}
