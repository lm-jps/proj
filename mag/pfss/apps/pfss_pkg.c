/*
 * Name:		pfss_pkg.c
 *
 * Description:		This is the C package of PFSS model.
 *			Including harmonic expansion, and field vector computation
 *                      Original code of the main algorithm is in IDL by Xuepu Zhao.
 *			Restructured and optimized by Xudong Sun.
 *			Refer to pfssnotes.pdf for details.
 *
 * Function List:
 *			void rdr(float r, double *rrr, double *D_rrr, int lmax, float apar)
 *			void pdpth(float th, double *leg, double *D_leg, int lmax)
 *			void csmph(float ph, double *mgh, double *D_mgh, 
 *				float *g, float *h, int lmax)
 *			void gh_pfss(float *map, float *g, float *h, struct Grid *grid, 
 *				int lmax, int sinlat)
 *			void Brtp(float r, float sint, float *Br0, float *Bt0, float *Bp0, 
 *				struct CombinedCoeff *CmbCoeff, int lmax)
 *			void Bcube(float *g, float *h, struct Grid *grid, float *Br,
 *                              float *Bt, float *Bp, int lmax, float apar, void (*rfunc)()) 
 *                      void Bpoint(float *g, float *h, struct Point *pt, 
 *                              float *Bvec, int lmax, float apar, void (*rfunc)())
 *			void pfss(float *g, float *h, struct Grid *grid, 
 *				float *Br, float *Bt, float *Bp, int lmax, float apar)
 *			void pfss0(float *g, float *h, struct Point *pt, 
 *				float *Bvec, int lmax, float apar)
 *
 * Called by:		jsynop2gh.c, jpfssfield.c, jgh2wsa.c, jgh2hccsss.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Jul 26 2009
 *			v1.1		Aug 06 2009
 *			v1.2		Aug 08 2009
 *			v1.3		Oct 13 2009
 *			v1.4		Oct 29 2009
 *
 * Issues:
 *			v1.1
 *			Added gh computation on latitude grid
 *			Moved constants and data structures to pfss.h, called from pfss.h
 *			v1.2
 *			Fixed a bug in pfss0 (free memories)
 *			v1.3
 *			Changed gh_pfss(), use grid structure instead of np, nt
 *			v1.4
 *			Added Bcube() and Bpoint()
 *                      Take the major part of pfss and pfss0 out and put them into Bcube
 *			and Bpoint so they could be shared by hccsss etc as well
 *			Added apar to rdr(), pfss() and pfss0() so it comforms with the
 *			prototype used Bcube and Bpoint
 *
 */


#include "pfss_pkg.h"


/* ################ Coefficients in three direction ################ */


/* Radial part for PFSS */
void rdr(float r, double *rrr, double *D_rrr, int lmax, float apar)
{
    int Dptr;
    double R0r, rRS, rRS2;
    double R0rl1, rRS2l1, R0RS2l1;
    R0r = R0/r; rRS = r/RS; rRS2 = rRS*rRS;
    R0rl1 = R0r; rRS2l1 = rRS; R0RS2l1 = R0RS;

    for (int l = 0; l <= lmax; l++) {
        Dptr = l;
        double denom = l + 1 + l * R0RS2l1;
        rrr[Dptr] = R0 * R0rl1 * (1 - rRS2l1) / denom;
        D_rrr[Dptr] = -R0r * R0rl1 * (l + 1 + l * rRS2l1) / denom;
        R0rl1 *= R0r; rRS2l1 *= rRS2; R0RS2l1 *= R0RS2;
    }
}


/* Theta part, Legendre function
 * Schmidt normalized to (2-dm0)*2/(2l-1)
 * P_m^m = (2-dm0)^0.5 ((2m-1)!!/2m!!)^0.5 (1-x^2)^{m/2}, where dm0=1, m=0; dm0=0, m>0
 * P_{l+1}^m = ((l+1)^2-m^2)^{-0.5}[(2l+1) x P_l^m - (l^2-m^2)^0.5 P_{l-1}^m]
 * See pfssnotes.pdf for details
 */
void pdpth(float th, double *leg, double *D_leg, int lmax)
{
    int Dptt;
    double cost, sint;
    double t0, t1, t2, dt0, dt1, dt2, expmi;
    double m2, l1, l12, s1, s2, f1;

    cost = cos(th); sint = sin(th);
    for (int m = 0; m <= lmax; m++) {
        t1 = (m == 0) ? 1.0 : sqrt(2.0);
        dt1 = t1 * m * cost / sint;
        for (int mi = 0; mi < m; mi++) {
            expmi = sint * sqrt((2.0 * mi + 1.0) / (2.0 * mi + 2.0));
            t1 *= expmi; dt1 *= expmi;
        }
        Dptt = (m + 1) * m / 2 + m;
        leg[Dptt] = t1; D_leg[Dptt] = dt1;
        m2 = (double)(m) * m; t0 = 0.0; dt0 = 0.0; s1 = 0.0;
        for (int l = m; l < lmax; l++) {
            l1 = l + 1.0; l12 = l1 * l1;
            s2 = sqrt(l12 - m2);
            f1 = 2.0 * l + 1.0;
            t2 = (f1 * cost * t1 - s1 * t0) / s2;
            dt2 = (f1 * (cost * dt1 - sint * t1) - s1 * dt0) / s2;
            Dptt = (l + 2) * (l + 1) / 2 + m;
            leg[Dptt] = t2; D_leg[Dptt] = dt2;  
            s1 = s2; t0 = t1; t1 = t2; dt0 = dt1; dt1 = dt2;                  
        }
    }
}


/* Phi part, combine g and h with phi */
void csmph(float ph, double *mgh, double *D_mgh, float *g, float *h, int lmax)
{
    int Dptp, Dpt;
    double sinmp, cosmp;

    for (int m = 0; m <= lmax; m++) {
    cosmp = cos(m * ph); sinmp = sin(m * ph);
        for (int l = m; l <= lmax; l++) {
            Dptp = (l + 1) * l / 2 + m;
            Dpt = l * (lmax + 1) + m;
            mgh[Dptp] = g[Dpt] * cosmp + h[Dpt] * sinmp;
            D_mgh[Dptp] = -m * g[Dpt] * sinmp + m * h[Dpt] * cosmp;
        }
    }
}



/* ################ GH Coefficient ################ */



/* Use numerical integration to get gh */
void gh_pfss(float *map, float *g, float *h, struct Grid *grid, 
	int lmax, int sinlat)
{
    double *cmph, *smph;	// cos(m*phi), sin(m*phi)
    double *leg, *D_leg;	// P_l^m(cos(theta)), dP/dcos
    double glm, hlm;

    int np = grid->np, nt = grid->nt;
    float ph, th, kl;
    float kth, sinthdth, mdcosth[nt];
    int Dpt, Dptp, Dptt;
    int i, j, l, m;
    int lmax1 = lmax + 1;
    int lmax12 = lmax1 * (lmax1 + 1) / 2;
    int nplmax1 = np * lmax1;
    int ntlmax12 = nt * lmax12;

    cmph = (double *)(calloc(nplmax1, sizeof(double)));
    smph = (double *)(calloc(nplmax1, sizeof(double)));
    leg = (double *)(calloc(ntlmax12, sizeof(double)));
    D_leg = (double *)(calloc(ntlmax12, sizeof(double)));

    for (i = 0; i < np; i++) {
        ph = grid->ph[i];
//        ph = i * 1.0 / np * 2 * M_PI;	// FFT grid
//        ph = (i + 1.0) / np * 2 * M_PI;	// For SSW
        for (m = 0; m <= lmax; m++) {
            Dptp = i * lmax1 + m;
            cmph[Dptp] = cos(m * ph); smph[Dptp] = sin(m * ph);
        }
    }

    for (j = 0; j < nt; j++) {
        th = grid->th[j];
        mdcosth[j] = cos((nt - j - 1) * M_PI / nt) - 
            cos((nt - j) * M_PI / nt);	// For lat grid
        Dptt = j * lmax12;
        pdpth(th, (leg + Dptt), (D_leg + Dptt), lmax);
    }

    /* Integrate */
    sinthdth = 2.0 / nt;
    for (l = 0; l <= lmax; l++) {
    kl = (2.0 * l + 1.0) / (2.0 * np);
    for (m = 0; m <= l; m++) {
        glm = 0.0; hlm = 0.0;
        for (j = 0; j < nt; j++) {
        kth = sinlat ? sinthdth : mdcosth[j];
        for (i = 0; i < np; i++) {
            Dpt = j * np + i;
            Dptp = j * lmax12 + l * (l + 1) / 2 + m;
            Dptt = i * lmax1 + m;
            glm += map[Dpt] * leg[Dptp] * cmph[Dptt] * kth;
            hlm += map[Dpt] * leg[Dptp] * smph[Dptt] * kth;
        }	// i
        }	// j
        g[l * lmax1 + m] = glm * kl;
        h[l * lmax1 + m] = hlm * kl;
    }	// m
    h[l * lmax1] = 0.0;
    }	// l

    free(cmph); free(smph);
    free(leg); free(D_leg);
}



/* ################ Field Strength ################ */


/* Compute Bvec at a point, geometries in CmbCoeff */
void Brtp(float r, float sint, float *Br0, float *Bt0, float *Bp0, 
	struct CombinedCoeff *CmbCoeff, int lmax)
{
    int Dptr, Dptt, Dptp;
    double fr, dfr, ft, dft, fp, dfp;
    double x1 = 0.0, x2 = 0.0, x3 = 0.0;
    for (int l = 0; l <= lmax; l++) {
        Dptr = l;
        Dptt = (l + 1) * l / 2;
        Dptp = (l + 1) * l / 2;
        for (int m = 0; m <= l; m++) {
            fr = CmbCoeff->rrr[Dptr];
            dfr = CmbCoeff->D_rrr[Dptr];
            ft = CmbCoeff->leg[Dptt + m];
            dft = CmbCoeff->D_leg[Dptt + m];
            fp = CmbCoeff->mgh[Dptp + m];
            dfp = CmbCoeff->D_mgh[Dptp + m];
            x1 += (dfr * ft * fp);
            x2 += (fr * dft * fp);
            x3 += (fr * ft * dfp);
        }
    }
    *Br0 = -x1;
    *Bt0 = -x2 / r;
    *Bp0 = -x3 / r / sint;
}


/* Common field function for cube used by all PFSS-like models,
 * Coefficient pre-computed and stored so more efficient
 */
void Bcube(float *g, float *h, struct Grid *grid, float *Br,
	float *Bt, float *Bp, int lmax, float apar, 
	void (*rfunc)(float, double *, double *, int, float))
{
    int np = grid->np, nt = grid->nt, nr = grid->nr;
    int npnt = np * nt;
    int lmax1 = lmax + 1;
    int lmax12 = lmax1 * (lmax1 + 1) / 2;
    int nrlmax1 = nr * lmax1;
    int ntlmax12 = nt * lmax12;
    int nplmax12 = np * lmax12;
	
    // Collection of g and h with grid parameters
    int sz = sizeof(double);
    double *rrr = (double *)(calloc(nrlmax1, sz));
    double *D_rrr = (double *)(calloc(nrlmax1, sz));
    double *leg = (double *)(calloc(ntlmax12, sz));
    double *D_leg = (double *)(calloc(ntlmax12, sz));
    double *mgh = (double *)(calloc(nplmax12, sz));
    double *D_mgh = (double *)(calloc(nplmax12, sz));

    // Combine coeffs
    for (int n = 0; n < nr; n++)
        rfunc(grid->rr[n], (rrr + n * lmax1), (D_rrr + n * lmax1), lmax, apar);
    for (int j = 0; j < nt; j++)
        pdpth(grid->th[j], (leg + j * lmax12), (D_leg + j * lmax12), lmax);
    for (int i = 0; i < np; i++)
        csmph(grid->ph[i], (mgh + i * lmax12), (D_mgh + i * lmax12), g, h, lmax);

    // Coeff for each point
    struct CombinedCoeff *CmbCoeff = 
        (struct CombinedCoeff *)(malloc(sizeof(struct CombinedCoeff)));

    for (int n = 0; n < nr; n++) {
        float r = grid->rr[n];
        CmbCoeff->rrr = rrr + n * lmax1;
        CmbCoeff->D_rrr = D_rrr + n * lmax1;
        for (int j = 0; j < nt; j++) {
            float sint = sin(grid->th[j]);
            CmbCoeff->leg = leg + j * lmax12;
            CmbCoeff->D_leg = D_leg + j * lmax12;
            for (int i = 0; i < np; i++) {
                int Dpt = n * npnt + j * np + i;
                float *B1 = Br + Dpt;
                float *B2 = Bt + Dpt;
                float *B3 = Bp + Dpt;
                CmbCoeff->mgh = mgh + i * lmax12;
                CmbCoeff->D_mgh = D_mgh + i * lmax12;
                // Point by point
                Brtp(r, sint, B1, B2, B3, CmbCoeff, lmax);
            }
        }
    }

    free(rrr); free(D_rrr);
    free(leg); free(D_leg);
    free(mgh); free(D_mgh);
    free(CmbCoeff);
}

/* Common field function for point used by all PFSS-like models */
void Bpoint(float *g, float *h, struct Point *pt, 
	float *Bvec, int lmax, float apar, 
	void (*rfunc)(float, double *, double *, int, float))
{
    int lmax1 = lmax + 1;
    int lmax12 = lmax1 * (lmax1 + 1) / 2;
    int sz = sizeof(double);

    struct CombinedCoeff *CmbCoeff = 
        (struct CombinedCoeff *)(malloc(sizeof(struct CombinedCoeff)));
    CmbCoeff->rrr = (double *)(calloc(lmax1, sz));
    CmbCoeff->D_rrr = (double *)(calloc(lmax1, sz));
    CmbCoeff->leg = (double *)(calloc(lmax12, sz));
    CmbCoeff->D_leg = (double *)(calloc(lmax12, sz));
    CmbCoeff->mgh = (double *)(calloc(lmax12, sz));
    CmbCoeff->D_mgh = (double *)(calloc(lmax12, sz));

    rfunc(pt->r, CmbCoeff->rrr, CmbCoeff->D_rrr, lmax, apar);
    pdpth(pt->t, CmbCoeff->leg, CmbCoeff->D_leg, lmax);
    csmph(pt->p, CmbCoeff->mgh, CmbCoeff->D_mgh, g, h, lmax);

    float B1, B2, B3;
    Brtp(pt->r, sin(pt->t), &B1, &B2, &B3, CmbCoeff, lmax);
    Bvec[0] = B1; Bvec[1] = B2; Bvec[2] = B3;

    free(CmbCoeff->rrr); free(CmbCoeff->D_rrr);
    free(CmbCoeff->leg); free(CmbCoeff->D_leg);
    free(CmbCoeff->mgh); free(CmbCoeff->D_mgh);
    free(CmbCoeff);
}



/* ############## Calls for PFSS ############## */


/* Assemble gh into B vectors on a single point */
void pfss0(float *g, float *h, struct Point *pt, 
	float *Bvec, int lmax, float apar)
{
    void (*rfunc)(float, double *, double *, int, float);
    rfunc = rdr;
    Bpoint(g, h, pt, Bvec, lmax, 0.0, rfunc);
}


/* Assemble gh into B vectors on the grid */
void pfss(float *g, float *h, struct Grid *grid, 
	float *Br, float *Bt, float *Bp, int lmax, float apar)
{
    void (*rfunc)(float, double *, double *, int, float);
    rfunc = rdr;
    Bcube(g, h, grid, Br, Bt, Bp, lmax, 0.0, rfunc);
}
