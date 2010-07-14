/*
 * Name:		wsa_pkg.c
 *
 * Description:		This is the C package of WSA model.
 *			Including harmonic expansion, and field vector computation
 *                      Original code of the main algorithm is in IDL mainly by C. N. Arge,
 *			partly written for MDI data in IDL by Xudong Sun.
 *			Now restructured and optimized by Xudong Sun.
 *
 * Function List:
 *			float gcd(float lon1, float lat1, float lon2, float lat2)
 *			void getbd(int *fpmap, float *fte, float *fpph, float *fpth, 
 *				int np, int nt)
 *			float ang_sep(float ph, float th, float *fpmap, int np, int nt)
 *			void fte_global(float *g, float *h, struct Grid *grid, float *fte, 
 *				float *fpph, float *fpth, float *brss, int lmax, void (*integrator)())
 *			void fte_subearth(float *g, float *h, struct Grid *gridsub, float *fte_sub, 
 *				float *fpph_sub, float *fpth_sub, float *brss_sub, int lmax, 
 *				void (*integrator)())
 *			void propagate(double *t, float *v, float *b, double *t1AU, float *v1AU, 
 *				float *b1AU, int np)
 *			void interpol(float *y0, double *x0, int n0, float *y, double *x, int n)
 *			void fixts1au(double *t, float *v, float *b, int n, 
 *				double *date, float *speed, float *imf, int num)
 *			void wsa_ss(float *g, float *h, int lmax, int rk4, 
 *				struct Grid *grid, struct Gridsub *gridsub, 
 *				float *v_g, float *brss_g, float *fte_g, float *fpph_g, float *fpth_g, 
 *				float *v_s, float *brss_s, float *fte_s, float *fpph_s, float *fpth_s)
 *			void wsa_1AU(struct Gridsub *gridsub, float *v_s, float *brss_s, float *date, 
 *				float *speed, float *imf, int tslen)
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *			IDL code developed by C. N. Arge (Nick.Arge@Kirtland.af.mil)
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Jul 29 2009
 *			v1.1		Jul 31 2009
 *			v1.2		Aug 03 2009
 *			v1.2a		Aug 06 2009
 *			v1.2b		Aug 07 2009
 *			v1.3		Oct 29 2009
 *			v1.4		Apr 19 2010
 *
 * Issues:
 *			v1.0
 *			Fix uncomputed points?
 *			Time and Solar coordinate conversion code ?????
 *			Not sure about the crude angular separation algorithm :-(
 *			v1.1
 *			Added suntime.c using Arge's algorithm (need ctime later)
 *			Added wsa_ss function for external use
 *			v1.2
 *			Added propagate function
 *			Added outside control of integrator
 *			v1.2a
 *			Changed Gridsub->jd, from float to double
 *			Added interpolation code to get equi-spaced temporal grid
 *			Added wsa_1AU()
 *			Moved constants and data structures to pfss.h, called from there
 *			v1.2b
 *			Changed imf polarity definition (to geocentric)
 *			v1.3
 *			Added Bfunc to fte_global() and fte_subearth()
 *			v1.4
 *			Changed time to sec, changed imf to int, date to double
 *
 */

#include "wsa_pkg.h"


/* ################ Empirical Function ################ */


/* WSA empirical function for speed (km/s) at 2.5Rsun (Arge)
 * This version is for MDI, 2.5d resolution, latitude format
 * Refer to X. Sun's PFSS paper for details.
 */
float wsa_f(float f, float th)
{
     float fs_exp, bb, kk, th0, th_exp, all_exp;
     float fac0, expo, fac1, spd;
     fs_exp = 0.22; bb = 6.08; kk = 3.30;
     th0 = 2.50; th_exp = 1.50; all_exp = 3.40;
     fac0 = pow((1. + fabs(f)), fs_exp);
     expo =  1 - pow((th / th0), th_exp);
     fac1 = bb - kk * exp(expo) / exp(1.0);
     spd = 265. + 1.5 * pow(fac1, all_exp) / fac0;
     return spd;
}



/* ################ Angular separation ################ */


/* Great circle distance between pt1 and pt2 */
float gcd(float lon1, float lat1, float lon2, float lat2)
{
    double dlon = lon1 - lon2;
    double dlat = lat1 - lat2;
    if ((fabs(dlon) < 1.E-5) && (fabs(dlat) < 1.E-5))
        return 0.;
    double num1 = cos(lat2) * sin(dlon);
    double num2 = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon);
    double num = sqrt(num1 * num1 + num2 * num2);
    double den = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon);
    float result = fabs(atan2(num, den));
    return (result);
}


/* Smooth foot point map to avoid small closed structure */
void smoothfp(int *fpmap, int np, int nt)
{
    int Dpt, Dpt0;
    for (int j = 0; j < nt; j++) {
        for (int i = 0; i < np; i++) {
            Dpt = j * np + i;
            if (fpmap[Dpt]) continue;	// Already open
            int count = 0;
            for (int ii = -1; ii <= 1; ii++) {
                for (int jj = -1; jj <= 1; jj++) {
                    if (ii == 0 && jj == 0) continue;
                    int lon = i + ii;
                    int lat = j + jj;
                    if (lon < 0) lon += np;	// Avoid % operator
                    if (lon >= np) lon -= np;
                    if (lat < 0) {
                        lat = abs(lat);
                        lon = lon + np / 2;
                        if (lon >= np) lon -= np;
                    }
                    if (lat >= nt) {
                        lat = 2 * nt - lat - 1;
                        lon = lon + np / 2;
                        if (lon >= np) lon -= np;
                    }
                    Dpt0 = lat * np + lon;
                    if (fpmap[Dpt0]) count++;
                }
            }
            if (count >= 5) fpmap[Dpt] = 1;
        }
    }
}


/* At photosphere, note open/close at each grid point
 * based on downward tracing result, then do smoothing
 */
void getbd(int *fpmap, float *fte, float *fpph, float *fpth, 
	int np, int nt)
{
    float dph = 360. / np, dth = 180. / nt;
    // Make map
    int Dpt, Dpt0;
    for (int j = 0; j < nt; j++) {
        for (int i = 0; i < np; i++) {
            Dpt = j * np + i;
            if (fabs(fte[Dpt]) > 1.E-5) {
                int ii = (int)(fpph[Dpt] / DTOR / dph);
                int jj = nt - 1 - (int)(fpth[Dpt] / DTOR / dth);
                Dpt0 = jj * np + ii;
                fpmap[Dpt0] = 1;		// For open foot point
            }
        }
    }
    // "Smoothing"
    for (int i = 0; i < 15; i++)
        smoothfp(fpmap, np, nt);
}


/* Find minimum angular separation from one foot point
 * to the edge of the coronal hole, start from 40d, 
 * gradually narrow down around in the closed area. (Arge)
 */
float ang_sep(float ph, float th, int *fpmap, int np, int nt)
{
    float dph = 360. / np, dth = 180. / nt;
    float xx = ph, yy = M_PI / 2. - th;
    float dist = 181.0, dist_t = 0.0;

    int i, ii, jj;
    int Dpt;
    float dx, dy;
    int x_grid, y_grid;
    int box_x, box_y;
    int left, right, up, down;
    float lon, lat;

    dx = dph * DTOR;
    dy = dth * DTOR;
    x_grid = (int)(ph / dx);
    y_grid = (int)((M_PI - th) / dy);
    box_x = (int)(BOX / dph / cos(yy));
    box_y = (int)(BOX / dth);
    if (box_x >= (np / 2)) {
        left = 0;
        right = np - 1;
    } else {
        left = x_grid - box_x;
        right = x_grid + box_x;
    }
    up = y_grid + box_y;
    down = y_grid - box_y;
    for (int aa = left; aa <= right; aa++) {
        if (aa < 0)
            i = aa + np;
        else if (aa >= np)
            i = aa - np;
        else
            i = aa;
        for (int bb = down; bb <= up; bb++) {
            if (bb < 0) {
                ii = i + np / 2;
                if (ii >= np) ii -= np;
                jj = abs(bb);
            } else if (bb >= nt) {
                ii = i + np / 2;
                if (ii >= np) ii -= np;
                jj = 2 * nt - bb - 1;
            } else {
                ii = i;
                jj = bb;
            }
            Dpt = jj * np + ii;
            if ((ii != x_grid) && (jj != y_grid) && !fpmap[Dpt]) {
                lon = (ii + 0.5) / np * 2 * M_PI;
                lat = (jj + 0.5) / nt * M_PI -  M_PI / 2.;
                dist_t = gcd(xx, yy, lon, lat) / DTOR;
                if (dist_t < dist) dist = dist_t;
            }
        }
    }
    return dist;
}



/* ################ Getting the physical variables ################ */


/* Global field line tracing for flux tube expansion
 * factor and foot point locations. Brss as a byproduct.
 */
void fte_global(float *g, float *h, struct Grid *grid, float *fte, 
	float *fpph, float *fpth, float *brss, int lmax, 
	void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float)))
{
    int np = grid->np, nt = grid->nt;
    int failed = 0;	// Number of points where tracing fails
    // Starting point
    struct Point *pt0 = (struct Point *)(malloc(sizeof(struct Point)));
    pt0->r = RS;
    struct FldLn *fl = (struct FldLn *)(malloc(sizeof(struct FldLn)));
    fl->num = 0;
    void (*Bfunc)(float *, float *, struct Point *, float *, int, float);
    Bfunc = pfss0;

    for (int i = 0; i < np; i++) {
        for (int j = 0; j < nt; j++) {
            pt0->t = grid->th[j];
            pt0->p = grid->ph[i];
            // Trace field line from source surface to photosphere
            trace(g, h, pt0, fl, lmax, integrator, Bfunc, RS, R0 + 0.001, 0.0);
            // Result
            int Dpt = j * np + i;
            brss[Dpt] = fl->brss;
            if (fl->num) {		// Field tracing success
                fte[Dpt] = fl->fte;
                fpth[Dpt] = fl->tt[fl->num];
                fpph[Dpt] = fl->pp[fl->num];
            } else {
                fte[Dpt] = 0.;	// Field tracing failed
                failed++;
            }
//          printf("i=%d, j=%d, fte=%f\n", i, j, fte[Dpt]);
        }
    }
    free(fl); free(pt0);
//    Fix failed points??
//    if (failed) fte_global_fix(fte, fpph, fpth);
}


/* Subearth field line tracing, need global result first
 * Results are matrices fte, fpph, fpth, brss (NPHx3)
 * Row 1 will be subearth result, row 0 is DTH south, 
 * row 2 DTH north, following the regular map convention.
 */
void fte_subearth(float *g, float *h, struct Gridsub *gridsub, float *fte_sub, 
	float *fpph_sub, float *fpth_sub, float *brss_sub, int lmax, 
	void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float)))
{
    int np = gridsub->np, nt = gridsub->nt;
    int failed = 0;	// Number of points where tracing fails
    // Starting point
    struct Point *pt0 = (struct Point *)(malloc(sizeof(struct Point)));
    pt0->r = RS;
    struct FldLn *fl = (struct FldLn *)(malloc(sizeof(struct FldLn)));
    fl->num = 0;
    void (*Bfunc)(float *, float *, struct Point *, float *, int, float);
    Bfunc = pfss0;

    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < np; i++) {
            pt0->p = gridsub->ph[i];
            int Dpt = j * np + i;
            pt0->t = gridsub->th[Dpt];
            // Trace field line from source surface to photosphere
            trace(g, h, pt0, fl, lmax, integrator, Bfunc, RS, R0 + 0.001, 0.0);
            // Result
            brss_sub[Dpt] = fl->brss;
            if (fl->num) {		// Field tracing success
                fte_sub[Dpt] = fl->fte;
                fpth_sub[Dpt] = fl->tt[fl->num];
                fpph_sub[Dpt] = fl->pp[fl->num];
            } else {
                fte_sub[Dpt] = 0.;	// Field tracing failed
                failed++;
            }
        }
    }
    free(fl); free(pt0);
//    Fix failed points??
//    if (failed) fte_subearth_fix(fte_sub, fpph_sub, fpth_sub);
}



/* ################ Propagating to 1AU ################ */


/* v, b: Series at source surface; v1AU, b1AU: Series at 1AU */
void propagate(double *t, float *v, float *b, double *t1AU, float *v1AU, 
	float *b1AU, int np)
{
    double step = PROPDIST * RSUN / NSTEP;
    double t1AU_old[np];
    float v1AU_old[np];
    float b1AU_old[np];
    for (int i = 0; i < np; i++) {
        t1AU[i] = t[i];
        v1AU[i] = v[i];
        b1AU[i] = b[i];
    }
    for (int rrr = 0; rrr < NSTEP; rrr++) {
        for (int i = 0; i < np; i++) {
            t1AU_old[i] = t1AU[i];
            v1AU_old[i] = v1AU[i];
            b1AU_old[i] = b1AU[i];
        }
        // Trial step
        for (int i = 0; i < np; i++) {
            t1AU[i] = t1AU_old[i] + (step / v1AU_old[i]);	// Changed t1AU to sec Apr 19 2010
        }
        // Adjustment
        v1AU[0] = v1AU_old[0];
        b1AU[0] = b1AU_old[0];
        for (int qq = 1; qq < np; qq++) {
            // Stream parcels interacting
            if (t1AU[qq] < t1AU[qq - 1]) {
                t1AU[qq] = t1AU[qq - 1] + 0.0001;
                double dtime = t1AU[qq] - t1AU_old[qq];
                float v0 = step / dtime;	// Changed dtime to sec Apr 19 2010
                float v1 = v1AU_old[qq];
                v1AU[qq] = sqrt(2.0 / (1.0 / (v0 * v0) + 1.0 / (v1 * v1)));
                if (v1AU[qq] < v1AU[qq - 1])
                    v1AU[qq] = v1AU[qq - 1];
                b1AU[qq] = (b1AU_old[qq] + b1AU_old[qq - 1]) / 2.0;    
            } else {
                v1AU[qq] = v1AU_old[qq];
                b1AU[qq] = b1AU_old[qq];
            }
        }
    }
    float rar2 = (PROPDIST / RS) * (PROPDIST / RS);
    for (int i = 0; i < np; i++) b1AU[i] /= rar2;
}


/* Linear interpolation function */
void interpol(float *y0, double *x0, int n0, float *y, double *x, int n)
{
     if (n0 < 2) {
         printf("Not enough points in given series.\n");
         return;
     }
     int i1, i2;
     for (int i = 0; i < n; i++) {
         if (x[i] <= x0[0]) {
             i1 = 1; i2 = 0;
         } else if (x[i] > x0[n0 - 1]) {
             i1 = n0 - 2; i2 = n0 - 1;
         } else {
             for (int j = 0; j < n0 - 2; j++)
                 if (x0[j] < x[i] && x0[j + 1] >= x[i]) {
                     i1 = j; i2 = j + 1;
                     break;
                 }
         }
         y[i] = y0[i1] + (y0[i2] - y0[i1]) /
                 (x0[i2] - x0[i1]) * (x[i] - x0[i1]);
     }
}


/* This function takes the predicted series at 1AU,
 * average and interpolate them onto a uniform grid.
 */
void fixts1au(double *t, float *v, float *b, int n, 
	double *date, float *speed, int *imf, int num)
{
     // Temporary arrays
     double t0[n];
     float v0[n], b0[n];
     for (int i = 0; i < n; i++) {
         t0[i] = 0; v0[i] = 0; b0[i] = 0;
     }
     // Average parcels that arrive at same time
     int n1 = 0, n2 = 1, n0 = 0;
     while (n1 < n) {
         while (n2 < n && fabs(t[n2] - t[n1]) < 90.)	// 15min, changed to secs Apr 19
             n2++;
         for (int i = n1; i < n2; i++) {
             t0[n0] += t[i]; v0[n0] += v[i]; b0[n0] += b[i];
         }
         int cc = n2 - n1;
         t0[n0] /= cc; v0[n0] /= cc; b0[n0] /= cc;
         n1 = n2; n2 = n1 + 1; n0++;
     }
     // Interpolate to date, clean up
     float imf_t[num];
     interpol(v0, t0, n0, speed, date, num);
     interpol(b0, t0, n0, imf_t, date, num);
     // For Bx at 1AU, from (+Br) is -, away (-Br) is +
     for (int i = 0; i < num; i++) {
         if (imf_t[i] > 0)
             imf[i] = -1;
         else if (imf_t[i] < 0)
             imf[i] = 1;
         else
             imf[i] = 0;
     }
}



/* ################ Main WSA ################ */


/* This is the main WSA code, computes parameters
 * at the source surface. Note there are requirements
 * to the dimensions of the array. (No checking here...)
 */
void wsa_ss(float *g, float *h, int lmax, int rk4, 
	struct Grid *grid, struct Gridsub *gridsub, 
	float *v_g, float *brss_g, float *fte_g, float *fpph_g, float *fpth_g, 
	float *v_s, float *brss_s, float *fte_s, float *fpph_s, float *fpth_s)
{
    int np = grid->np, nt = grid->nt;
    int npnt = np * nt;
    int i, j, Dpt;
    int *fpmap;
    float ang;

    void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float));
    if (rk4) integrator = RK4; else integrator = Euler;	// In fieldline_pkg.c

    // Get global parameters
    fte_global(g, h, grid, fte_g, fpph_g, fpth_g, brss_g, lmax, integrator);
    // Get subearth parameters
    fte_subearth(g, h, gridsub, fte_s, fpph_s, fpth_s, brss_s, lmax, integrator);
    // Get boundary of CH
    fpmap = (int *)(calloc(npnt, sizeof(int)));
    getbd(fpmap, fte_g, fpph_g, fpth_g, np, nt);
    // Get angular separation, then speed (global)
    for (j = 0; j < nt; j++) {
        for (i = 0; i < np; i++) {
            Dpt = j * np + i;
            ang = ang_sep(fpph_g[Dpt], fpth_g[Dpt], fpmap, np, nt);
            v_g[Dpt] = wsa_f(fte_g[Dpt], ang);
        }
    }
    // Get angular separation, then speed (subearth)
    for (j = 0; j < 3; j++) {
        for (i = 0; i < np; i++) {
            Dpt = j * np + i;
            ang = ang_sep(fpph_s[Dpt], fpth_s[Dpt], fpmap, np, nt);
            v_s[Dpt] = wsa_f(fte_s[Dpt], ang);
        }
    }
    free(fpmap);
}


/* This is the main WSA code for 1AU prediction
 * v_s, brss_s from subearth computation, date is predefined
 * 1AU time grid. We propagate SW to 1AU ballistically and include
 * interaction between fast/slow parcels. Results are then
 * interpolated to the predefined grid.
 */
void wsa_1AU(struct Gridsub *gridsub, float *v_s, float *brss_s, double *date, 
	float *speed, int *imf, int tslen)
{
    int i, j;
    int np = gridsub->np;

    // Temporary arrays
    double t0[np], t1au[np];
    float v0[np], b0[np], v1au[np], b1au[np];
    for (i = 0; i < np; i++)
        t0[i] = gridsub->time[np - 1 - i];

    // Make date double
//    double date_tmp[tslen];
//    for (i = 0; i < tslen; i++) date_tmp[i] = date[i];
    float *speed_tmp;
    int *imf_tmp;

    for (j = 0; j < 3; j++) {
        for (i = 0; i < np; i++) {
            v0[i] = v_s[j * np + np - 1 - i];
            b0[i] = brss_s[j * np + np - 1 - i];
        }
        propagate(t0, v0, b0, t1au, v1au, b1au, np);
        speed_tmp = speed + j * tslen; imf_tmp = imf + j * tslen;
        fixts1au(t1au, v1au, b1au, np, date, speed_tmp, imf_tmp, tslen);
    }
}

