/*
 * Name:		pfss_pkg.c
 *
 * Description:	
 *  This is the C package of PFSS model.
 *	Including harmonic expansion, and field vector computation
 *  Original code of the main algorithm is in IDL by Xuepu Zhao.
 *	Restructured and optimized by Xudong Sun.
 *	Refer to XXX for details.
 *
 * Version:		v1.0		Apr 26 2017
 *
 */


#define R0      (1.0)		// Photosphere
#define RS      (2.5)		// Source surface
#define R0RS    (R0 / RS)
#define R0RS2   (R0RS * R0RS)


/* 
 * Radial part for PFSS
 *
 */

void rdr (double r, int lmax, double *rrr, double *D_rrr)
{
    
    double R0r = R0/r, rRS = r/RS, rRS2 = rRS*rRS;
    double R0rl1 = R0r, rRS2l1 = rRS, R0RS2l1 = R0RS;
    
    for (int l = 0; l <= lmax; l++) {
        double denom = l + 1 + l * R0RS2l1;
        rrr[l] = R0 * R0rl1 * (1 - rRS2l1) / denom;
        D_rrr[l] = -R0r * R0rl1 * (l + 1 + l * rRS2l1) / denom;
        R0rl1 *= R0r; rRS2l1 *= rRS2; R0RS2l1 *= R0RS2;
    }
    
}


// =================================================


/*
 * Phi part, combine g and h with phi
 *
 */

void csmph (double ph, int lmax, double *g, double *h, double *mgh, double *D_mgh)
{
    
    for (int m = 0; m <= lmax; m++) {
        double cosmp = cos(m * ph), sinmp = sin(m * ph);
        for (int l = m; l <= lmax; l++) {
            int Dptp = (l + 1) * l / 2 + m;
            int Dpt = l * (lmax + 1) + m;
            mgh[Dptp] = g[Dpt] * cosmp + h[Dpt] * sinmp;
            D_mgh[Dptp] = - g[Dpt] * sinmp + h[Dpt] * cosmp;        // m absorbed to mlegs
        }
    }
    
}


// =================================================


/* 
 * Theta part, Legendre function
 * Schmidt normalized to (2-dm0)*2/(2l-1)
 * P_m^m = (2-dm0)^0.5 ((2m-1)!!/2m!!)^0.5 (1-x^2)^{m/2}, where dm0=1, m=0; dm0=0, m>0
 * P_{l+1}^m = ((l+1)^2-m^2)^{-0.5}[(2l+1) x P_l^m - (l^2-m^2)^0.5 P_{l-1}^m]
 * Input: colat th, lmax
 * Output: Plm(costh), (for Bt) dPlm(costh)/dth,
 *      and (for Bp) m*Plm(costh)/sinth: array (lmax+1)*(lmax+2)/2
 *
 */

void pdpth (double th, int lmax, double *leg, double *D_leg, double *mlegs)
{
    int Dptt;
    double cost, sint;
    double t0, t1, t2, dt0, dt1, dt2, t0s, t1s, t2s, expmi;
    double m2, l1, l12, s1, s2, f1;
    
    cost = cos(th); sint = sin(th);
    for (int m = 0; m <= lmax; m++) {
        t1 = (m == 0) ? 1.0 : sqrt(2.0);
        t1s = t1;
        dt1 = t1 * m * cost;
        for (int mi = 0; mi < m; mi++) {
            expmi = sqrt((2.0 * mi + 1.0) / (2.0 * mi + 2.0));
            t1 *= (sint * expmi);
            if (mi != 0) {       // take care of the case t=0 or pi
                dt1 *= sint;    // one less sinth
                t1s *= sint;
            }
            dt1 *= expmi;
            t1s *= expmi;
        }
        Dptt = (m + 1) * m / 2 + m;
        leg[Dptt] = t1; D_leg[Dptt] = dt1; mlegs[Dptt] = t1s * m;
        
        m2 = (double)(m) * m; t0 = 0.0; t0s = 0.0; dt0 = 0.0; s1 = 0.0;
        for (int l = m; l < lmax; l++) {        // skip if m = lmax
            l1 = l + 1.0; l12 = l1 * l1;
            s2 = sqrt(l12 - m2);
            f1 = 2.0 * l + 1.0;
            t2 = (f1 * cost * t1 - s1 * t0) / s2;
            t2s = (f1 * cost * t1s - s1 * t0s) / s2;        // one less sinth
            dt2 = (f1 * (cost * dt1 - sint * t1) - s1 * dt0) / s2;
            Dptt = (l + 2) * (l + 1) / 2 + m;
            leg[Dptt] = t2; D_leg[Dptt] = dt2; mlegs[Dptt] = t2s * m;   // P(l+1,m)
            s1 = s2; t0 = t1; t1 = t2; t0s = t1s; t1s = t2s;
            dt0 = dt1; dt1 = dt2;
        }
    }
}


// =================================================


/*
 * Use numerical integration to get gh
 * p & t are grid locations, kp & kt are area weighting
 * for Q-map grids that include 0/2pi in longitude
 * and 0/pi in colatitude, the edges are weighted by 0.5
 * Under our normalization scheme, sum(kp#kt)=1.
 * this must be specified in the main code
 *
 */

void gh_pfss (double *br, int np, int nt, int lmax,
              double *p, double *t, double *kp, double *kt,
              double *g, double *h)
{
    
    int Dpt, Dptp, Dptt;        // pointer offset
    int lmax1 = lmax + 1;
    int lmax12 = lmax1 * (lmax1 + 1) / 2;
    int nplmax1 = np * lmax1;
    int ntlmax12 = nt * lmax12;
    
    // // cos(m*phi), sin(m*phi)
    
    double *cmph = (double *)(calloc(nplmax1, sizeof(double)));
    double *smph = (double *)(calloc(nplmax1, sizeof(double)));
    
    for (int i = 0; i < np; i++) {
        for (int m = 0; m <= lmax; m++) {
            Dptp = i * lmax1 + m;
            cmph[Dptp] = cos(m * p[i]);
            smph[Dptp] = sin(m * p[i]);
        }
    }
    
    // Legendre, P(costh), dP/dcos, m*P/sinth
    
    double *leg = (double *)(calloc(ntlmax12, sizeof(double)));
    double *D_leg = (double *)(calloc(ntlmax12, sizeof(double)));
    double *mlegs = (double *)(calloc(ntlmax12, sizeof(double)));
    
    for (int j = 0; j < nt; j++) {
        Dptt = j * lmax12;
        pdpth(t[j], lmax,
              (leg + Dptt), (D_leg + Dptt), (mlegs + Dptt));
    }
    
    // Integrate
    
    double w;
    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            double glm = 0.0, hlm = 0.0;
            for (int j = 0; j < nt; j++) {
                for (int i = 0; i < np; i++) {
                    w = kp[i] * kt[j];
                    Dpt = j * np + i;
                    Dptp = j * lmax12 + l * (l + 1) / 2 + m;
                    Dptt = i * lmax1 + m;
                    glm += (br[Dpt] * leg[Dptp] * cmph[Dptt] * w);
                    hlm += (br[Dpt] * leg[Dptp] * smph[Dptt] * w);
                }	// i
            }	// j
            g[l * lmax1 + m] = glm * (2 * l + 1.);
            h[l * lmax1 + m] = hlm * (2 * l + 1.);
        }	// m
        h[l * lmax1] = 0.0;
    }   // l
    g[0] = 0.0;     // no monopole
    
    free(cmph); free(smph);
    free(leg); free(D_leg); free(mlegs);
    
}


// =================================================



/*
 * Compute PFSS field vectors Bp, Bt, Br
 * on a regular grid p, t, r
 *
 */

void pfss_cube (double *g, double *h, int lmax,
                double *p, double *t, double *r, int np, int nt, int nr,
                double *Bp, double *Bt, double *Br)
{
    
    int npnt = np * nt;
    int lmax1 = lmax + 1;
    int lmax12 = lmax1 * (lmax1 + 1) / 2;
    int nrlmax1 = nr * lmax1;
    int ntlmax12 = nt * lmax12;
    int nplmax12 = np * lmax12;
    
    // Collection of g and h with grid parameters
    
    double *rrr = (double *)(calloc(nrlmax1, sizeof(double)));
    double *D_rrr = (double *)(calloc(nrlmax1, sizeof(double)));
    double *leg = (double *)(calloc(ntlmax12, sizeof(double)));
    double *D_leg = (double *)(calloc(ntlmax12, sizeof(double)));
    double *mlegs = (double *)(calloc(ntlmax12, sizeof(double)));
    double *mgh = (double *)(calloc(nplmax12, sizeof(double)));
    double *D_mgh = (double *)(calloc(nplmax12, sizeof(double)));
    
    // Combine coeffs
    
    for (int n = 0; n < nr; n++) {      // radial
        int Dpt = n * lmax1;      // offset
        rdr(r[n], lmax, (rrr + Dpt), (D_rrr + Dpt));
    }
    
    for (int j = 0; j < nt; j++) {      // theta
        int Dpt = j * lmax12;
        pdpth(t[j], lmax, (leg + Dpt), (D_leg + Dpt), (mlegs + Dpt));
    }
    
    for (int i = 0; i < np; i++) {      // phi
        int Dpt = i * lmax12;
        csmph(p[i], lmax, g, h, (mgh + Dpt), (D_mgh + Dpt));
    }
    
    // Compute each point
    
    for (int n = 0; n < nr; n++) {
        
        double *rrr_t = rrr + n * lmax1;
        double *D_rrr_t = D_rrr + n * lmax1;
        
        for (int j = 0; j < nt; j++) {
            
            double *leg_t = leg + j * lmax12;
            double *D_leg_t = D_leg + j * lmax12;
            double *mlegs_t = mlegs + j * lmax12;
            
            for (int i = 0; i < np; i++) {
                
                double *mgh_t = mgh + i * lmax12;
                double *D_mgh_t = D_mgh + i * lmax12;
                
                double x1 = 0.0, x2 = 0.0, x3 = 0.0;
                
                for (int l = 0; l <= lmax; l++) {
                    
                    int Dptr = l;
                    int Dptt = (l + 1) * l / 2;
                    int Dptp = (l + 1) * l / 2;
                    
                    for (int m = 0; m <= l; m++) {
                        
                        x1 += (D_rrr_t[Dptr] * leg_t[Dptt + m] * mgh_t[Dptp + m]);
                        x2 += (rrr_t[Dptr] * D_leg_t[Dptt + m] * mgh_t[Dptp + m]);
                        x3 += (rrr_t[Dptr] * mlegs_t[Dptt + m] * D_mgh_t[Dptp + m]);
                        
                    }   // m
                }   // l
                
                int Dpt = n * npnt + j * np + i;
                Br[Dpt] = - x1;
                Bt[Dpt] = - x2 / r[n];
                Bp[Dpt] = - x3 / r[n];
                
            }   // i
        }   // j
    }   // n
    
    // Clean up
    
    free(rrr); free(D_rrr);
    free(leg); free(D_leg); free(mlegs);
    free(mgh); free(D_mgh);
    
}

