/*
 *  rdfitc.c						~rick/src/ringfit/rdfitc
 *
 *  Responsible:  Rick Bogart                   	RBogart@spd.aas.org
 *     
 *  "Ring" fitting based on the code developed and used by Sarbani Basu
 *    and H. M. Antia for fitting 3-d solar acoustic power spectra; tbased
 *    on SOI module ringfit_v0.86
 *
 *  Main program begins around line 800
 *
 *  Bugs:
 *    Timing is not re-zeroed for each I/O record pair processed
 *    default values of lmin and lmax are appropriate for MDI full-disc
 *	resolution and areas of size 15 deg
 *    non-default value of nmax seems to force premature termination; in
 *	particular, nmax <= 3 causes crash due to going beyond array bounds
 *	of freq, wid
 *
 *  Planned development:
 *    1.0 parametrize test criteria for bad fits
 *    Use drms_keyword_lookup()
 *    Remove unnecessary defined constants, predefined size declarations
 *    Add parameters for fitting coefficients
 *    Remove unused variables, provide more descriptive names
 *    Generalize bfgs() to use function call in parameter
 *
 *  Revision history is at end of file.
 */

#include <jsoc_main.h>
						      /*  module identifier  */
char *module_name = "ringfit_bba";
char *module_desc = "ring fitting using method of Basu and Antia";
char *version_id = "1.0";

ModuleArgs_t module_args[] = {
  {ARG_DATASET,	"in", "", "input dataset"}, 
  {ARG_STRING, "guessfile", "",
      "name of file containing initial frequency guesses"},
  {ARG_STRING, "out", "fit.out",
      "name of output data series (or file) containing fit coefficients"},
  {ARG_INT, "nmin", "0", "lowest radial order to fit"},
  {ARG_INT, "nmax", "10", "highest radial order to fit"},
  {ARG_INT, "lmin", "80", "lowest degree to fit"},
  {ARG_INT, "lmax", "3000", "highest degree to fit"},
  {ARG_INT, "fmin", "900", "lowest frequency to fit [uHz]"},
  {ARG_INT, "fmax", "7000", "highest frequency index to fit"},
  {ARG_INT, "bfgsct", "125", "number of iterations for bfgs"},
  {ARG_INT, "linminct", "15", "number of iterations for linmin"},
  {ARG_FLOAT, "ux", "0.0", "initial guess for zonal flow speed [m/s]"},
  {ARG_FLOAT, "uy", "0.0", "initial guess for merridional flow speed [m/s]"},
  {ARG_FLAG, "n", "0", "no fitting (diagnostics only)"},      
  {ARG_FLAG, "v", "0", "verbose output"},      
  {ARG_FLAG, "x", "0", "extra reporting of number of iterations"},      
  {ARG_FLAG, "2", "0", "calculate second derivatives explicitly"},      
  {}
};
       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "Cadence", "LonSpan",
    "T_START", "T_STOP", "Coverage", "Size", "Width", "Height",
    "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel",
    "MapScale", "MapProj", "Map_PA", "RSunRef", "PosAng", "MAI",
    "Apode_f", "Apode_k_min", "Apode_k_max", "APODIZNG", "APOD_MIN", "APOD_MAX"};

#include <unistd.h>
#include <math.h>
#include <time.h>
#include "keystuff.c"

#define VSCALE	(-0.4372)
#define RSUN_MM	(696.0)

static char *time_elapsed (int cct, int *sum) {
/*
 *  print a string containing utime reported to millisec plus wall-clock
 *    time since last zeroing call
 *  allows for one wrap-around in clock, which occurs every 36 minutes,
 *    but results will still be garbage if consecutive calls are spaced
 *    more than 72 minutes apart.
 */
  static clock_t t_last, t_curr;
  static time_t wt_last, wt_curr;
  static double cloadavg[3], loadavg[3];
  static double t_count, scale = 1.0 / CLOCKS_PER_SEC;
  static int loadct, first = 1;
  static char str[64];
  int n;

  if (first) {
    t_last = clock ();
    wt_last = time (NULL);
    first = 0;
    *sum = 0;
    for (n = 0; n < 3; n++) loadavg[n] = 0.0;
    loadct = 0;
  }
  t_curr = clock ();
  t_count += scale * (t_curr - t_last);
/*
  if (t_curr < t_last) t_count += 4294.9673;
*/
  wt_curr = time (NULL);
  getloadavg (cloadavg, 3);
  for (n = 0; n < 3; n++) loadavg[n] += cloadavg[n];
  loadct++;
  
  sprintf (str, "%.3f (%d) sec (avg load: %.2f %.2f %.2f)", t_count,
      (int)(wt_curr - wt_last), loadavg[0]/loadct, loadavg[1]/loadct,
      loadavg[2]/loadct);
  t_last = t_curr;
  if (!cct) {
    *sum += wt_curr - wt_last;
    t_last = clock ();
    wt_last = time (NULL);
    t_count = 0.0;
    for (n = 0; n < 3; n++) loadavg[n] = 0.0;
    loadct = 0;
  }
  return str;
}

static int gauelm (int n, int num, double *a, double *x, double *det,
    int *indx, int lj, int flag) {

  double r1, t1;
  int i, j, k, km, l;

  if (n <= 0 || n > lj) return 101;

  if (flag < 2) {
						    /*  perform elimination  */
    *det = 1.0;
    for (k = 0; k < n-1; k++) {
				       /*  find maximum element in column k  */
      r1 = 0.0;
      km = k;
      for (l = k; l < n; l++)
	if (fabs (a[l*lj + k]) > r1) {
	  r1 = fabs (a[l*lj + k]);
	  km = l;
	}
      indx[k] = km;
      if (km != k) {
				      /*  interchange the rows if necessary  */
	for (l = k; l < n; l++) {
	  t1 = a[k*lj + l];
	  a[k*lj + l] = a[km*lj + l];
	  a[km*lj + l] = t1;
	}
	*det *= -1;
      }
      *det *= a[k*lj + k];
      if (a[k*lj + k] == 0.0) return 121;
      for (l = k+1; l < n; l++) {
	a[l*lj + k] /= a[k*lj + k];
	for (i = k+1; i < n; i++)
	  a[l*lj + i] -= a[l*lj + k] * a[k*lj + i];
      }
    }
    *det *= a[(n-1)*lj + n-1];
    indx[n-1] = n - 1;
    if (a[(n-1)*lj + n-1] == 0) return 121;
    if (flag == 1) return 0;
  }
			/*  solution for the num different right-hand sides  */
  for (j = 0; j < num; j++) {
						   /*  forward-substitution  */
    int ljj = lj * j;
    for (k = 0; k < n-1; k++) {
      if (k != indx[k]) {
	t1 = x[ljj + k];
	x[ljj + k] = x[ljj + indx[k]];
	x[ljj + indx[k]] = t1;
      }
      for (l = k+1; l < n; l++)
	x[ljj + l] -= a[l*lj + k] * x[ljj + k];
    }
						      /*  back-substitution  */
    x[ljj + n-1] /= a[(n-1)*lj + n-1];
    for (k = n-2; k >= 0; k--) {
      for (l = n-1; l >= k+1; l--)
	x[ljj + k] -= x[ljj + l] * a[k*lj + l];
      x[ljj + k] /= a[k*lj + k];
    }
  }
  return 0;
}

static void fcn1 (int n, double *x, double *f, double *g,
    double *kx, double *ky, double *k2, double *omeg, double *ppow,
    double pk, double *fm, double ki, int npix) {
/*
 *  Given values for the set of fitting parameters x[13], form the sums:
 *    log (likelihood function) = f = sum of log M_i + O_i / M_i
 *	(M_i are fits for the n data values kx, ky, omeg;
 *	 O_i are the n values of power, ppow[i]);
 *    merit function = fm = sum of (O_i / M_i - 1)^2;
 *    g[13] = sum of ? (partial derivatives of fm w.r.t parameters?)
 */
  static double vf = 1.0e-4;
  double pk1, k, kx2, ky2, bk, dk, om1, wd, o1, d1, p1;
  double logk, kpk1, bko1p1, bktrm;
  double p1x0, p1x1, p1x2, p1x3, p1x4, p1x5, p1x8, p1x9,
      p1x10, p1x11, p1x12, p1x13, pix11, pix12, pi, dpix;
  int i;
							/*  initialize sums  */
  *f = 0.0;
  *fm = 0.0;
  for (i = 0; i < n; i++) g[i] = 0.0;
				     /*  limit values of fitting parameters  */
/*
  if (x[0] > 100.0) x[0] = 100.0;
  if (x[0] < 1.0) x[0] = 1.0;
  if (x[5] > 100.0) x[5] = 100.0;
  if (x[5] < 1.0) x[5] = 1.0;
  if (x[6] > 100.0) x[6] = 100.0;
  if (x[6] < 1.0) x[6] = 1.0;
  if (x[6] < 20.0 && x[5] < 20.0) x[6] = 20.0;
*/
  if (x[7] > 1000.0) x[7] = 1000.0;
  if (x[7] < -1000.0) x[7] = -1000.0;
  if (x[8] > 8.0e4) x[8] = 8.0e4;
  if (x[8] < -8.0e4) x[8] = -8.0e4;
  if (x[9] > 5.0e4) x[9] = 5.0e4;
  if (x[9] < -5.0e4) x[9] = -5.0e4;
  if (x[10] > 4000.0) x[10] = 4000.0;
  if (x[10] < -2000.0) x[10] = -2000.0;
  if (x[11] > 2.0e4) x[11] = 2.0e4;
  if (x[11] < -2.0e4) x[11] = -2.0e4;
  if (x[12] > 5000.0) x[12] = 5000.0;
  if (x[12] < -5000.0) x[12] = -5000.0;

  pk1 = pk + x[10] * vf;
  bk = x[12] * vf;

  for (i = 0; i < npix; i++) {
    kx2 = kx[i] * kx[i];
/*
    k = hypot (kx[i], ky[i]);
    k2 = k * k;
    k2 = kx2 + ky[i] * ky[i];
*/
    k = sqrt (k2[i]);
    kx2 /= k2[i];
    ky2 = kx[i] * ky[i] / k2[i];
    dk = k - ki;
/*
    kpk1 = pow (k, pk1);
*/
    logk = log (k);
    kpk1 = exp (pk1 * logk);

    om1 = x[1] * kpk1 * logk * vf;
    wd = x[4] + x[11] * dk * vf;
    o1 = (omeg[i] - x[1] * kpk1 - x[2] * vf * kx[i]
        - x[3] * vf * ky[i]) / wd;
    bko1p1 = 1.0 + bk * o1;
    bktrm = bk * bk + bko1p1 * bko1p1;
    d1 = o1 * o1 + 1.0;
    p1 = exp (x[0] + vf * (x[7] * dk + x[8] * kx2
        + x[9] * ky2)) * bktrm / d1;
    p1x1 = p1;
    p1x0 = p1 * 2 * bk * bko1p1 / bktrm - 2 * p1 * o1 / d1;
    p1x0 /= wd;
    p1x2 = -p1x0 * kpk1;
    p1x3 = -p1x0 * kx[i] * vf;
    p1x4 = -p1x0 * ky[i] * vf;
    p1x5 = -p1x0 * o1;
    p1x8 = p1 * dk * vf;
    p1x9 = p1 * kx2 * vf;
    p1x10 = p1 * ky2 * vf;
    p1x11 = -p1x0 * om1;
    p1x12 = -p1x0 * o1 * dk * vf;
    p1x13 = p1 * 2 * (bk + bko1p1 * o1) * vf / bktrm;
    pix11 = exp (x[5]) / (k2[i] * k);
    pix12 = exp (x[6]) / (k2[i] * k2[i]);

    pi = p1 + pix11 + pix12;
    *f += log (pi) + ppow[i] / pi;
    dpix = 1.0 - ppow[i] / pi;
    *fm += dpix * dpix;
    dpix /= pi;
    g[0] += p1x1 * dpix;
    g[1] += p1x2 * dpix;
    g[2] += p1x3 * dpix;
    g[3] += p1x4 * dpix;
    g[5] += pix11 * dpix;      
    g[6] += pix12 * dpix;      
    g[4] += p1x5 * dpix;
    g[7] += p1x8 * dpix;
    g[8] += p1x9 * dpix;
    g[9] += p1x10 * dpix;
    g[10] += p1x11 * dpix;
    g[11] += p1x12 * dpix;
    g[12] += p1x13 * dpix;
  }
}

static void d2fcn (int n, double *x, double *f, double *g, double *dd,
    double *k,  double *kx, double *ky, double *k2, double *omeg,
    double *ppow, double pk, double ki, int npix) {
  static double vf = 1.0e-4;
  double dp1[20*20];
  double bk, bkf, d1, om1, o1, pi, pix11, pix12, pk1, wd;
  double p1, p1x0, p1x1, p1x2, p1x3, p1x4, p1x5, p1x8, p1x9,
      p1x10, p1x11, p1x12, p1x13;
  double pi2, spow, pow2m1, scale;
  double dk, dp1xx, dp1x6;
  double bko1p1, dptrm;
  double k3, k4, kx2, ky2;
  double kxk, kpk1, logk;
  int i, j, p, col, row;

  *f = 0.0;
  for (j = 0; j < n; j++) {
    row = n * j;
    g[j] = 0.0;
    for (i = 0; i < n; i++)
      dd[i + row] = 0.0;
  }
  pk1 = x[10] * vf + pk;

  for (i = 0; i < npix; i++) {
/*
    k2 = kx[i] * kx[i] + ky[i] * ky[i];
    k = sqrt (k2[i]);
*/
    k3 = k2[i] * k[i];
    k4 = k2[i] * k2[i];
    kxk = kx[i] / k[i];
    kx2 = vf * kxk * kxk;
    ky2 = vf * kxk * ky[i] / k[i];
    dk = k[i] - ki;
    logk = log (k[i]);
    kpk1 = exp (pk1 * logk);

    bk = x[12] * vf;
    om1 = x[1] * kpk1 * logk * vf;
    wd = x[4] + x[11] * dk * vf;
    o1 = (omeg[i] - x[1] * kpk1 - x[2] * vf * kx[i] -
        x[3] * vf * ky[i]) / wd;
    d1 = o1 * o1 + 1.0;
    bko1p1 = 1.0 + bk * o1;
    bkf = bk * bk + bko1p1 * bko1p1;
    p1 = exp (x[0] + vf * x[7] * dk + x[8] * kx2 +
        x[9] * ky2) * bkf / d1;
    p1x0 = 2 * ((p1 * bk * bko1p1 / bkf) - p1 * o1 / d1);
    p1x0 /= wd;
    p1x1 = p1;
    p1x2 = -p1x0 * kpk1;
    p1x3 = -p1x0 * (kx[i] * vf);
    p1x4 = -p1x0 * (ky[i] * vf);
    p1x5 = -p1x0 * o1;
    p1x8 = p1 * dk * vf;
    p1x9 = p1 * kx2;
    p1x10 = p1 * ky2;
    p1x11 = -p1x0 * om1;
    p1x12 = -p1x0 * o1 * dk * vf;
    p1x13 = p1 * 2 * (bk + bko1p1 * o1) * vf / bkf;
    dp1xx = (2 * p1) * (bk * bk / bkf -
        4 * o1 * bk * bko1p1 / (d1 * bkf) +
	4 * o1 * o1 / d1 / d1 - 1 / d1);
    dp1x6 = (4*o1*bk + 2 - 2*o1*(2*bk + 2*o1*bko1p1) / d1) * vf * p1 / bkf;
    dp1[0 + 20*0] = p1x1;
    dp1[0 + 20*1] = p1x2;
    dp1[0 + 20*2] = p1x3;
    dp1[0 + 20*3] = p1x4;
    dp1[0 + 20*4] = p1x5;
    dp1[0 + 20*7] = p1x8;
    dp1[0 + 20*8] = p1x9;
    dp1[0 + 20*9] = p1x10;
    dp1[0 + 20*10] = p1x11;
    dp1[0 + 20*11] = p1x12;
    dp1[0 + 20*12] = p1x13;
    dp1[1 + 20*(0)] = p1x2;
    dp1[2 + 20*(0)] = p1x3;
    dp1[3 + 20*(0)] = p1x4;
    dp1[4 + 20*(0)] = p1x5;
    dp1[7 + 20*(0)] = p1x8;
    dp1[8 + 20*(0)] = p1x9;
    dp1[9 + 20*(0)] = p1x10;
    dp1[10 + 20*(0)] = p1x11;
    dp1[11 + 20*(0)] = p1x12;
    dp1[12 + 20*(0)] = p1x13;
    dptrm = dp1xx * o1 / wd + p1x0;
    dp1[1 + 20*(1)] = dp1xx * (kpk1 / wd) * (kpk1 / wd);
    dp1[1 + 20*(2)] = dp1xx * kpk1 * kx[i] * vf / wd / wd;
    dp1[1 + 20*(3)] = dp1xx * kpk1 * ky[i] * vf / wd / wd;
    dp1[1 + 20*(4)] = dptrm * kpk1 / wd;
    dp1[1 + 20*(7)] = p1x2 * dk * vf;
    dp1[1 + 20*(8)] = p1x2 * kx2;
    dp1[1 + 20*(9)] = p1x2 * ky2;
    dp1[1 + 20*(10)] = dp1xx * om1 * kpk1 / wd / wd +
        p1x0 * om1 / x[1];
    dp1[1 + 20*(11)] = dptrm * dk * vf * kpk1 / wd;
    dp1[1 + 20*(12)] = -dp1x6 * kpk1 / wd;
    dp1[2 + 20*(1)] = dp1[1 + 20*(2)];
    dp1[3 + 20*(1)] = dp1[1 + 20*(3)];
    dp1[4 + 20*(1)] = dp1[1 + 20*(4)];
    dp1[7 + 20*(1)] = dp1[1 + 20*(7)];
    dp1[8 + 20*(1)] = dp1[1 + 20*(8)];
    dp1[9 + 20*(1)] = dp1[1 + 20*(9)];
    dp1[10 + 20*(1)] = dp1[1 + 20*(10)];
    dp1[11 + 20*(1)] = dp1[1 + 20*(11)];
    dp1[12 + 20*(1)] = dp1[1 + 20*(12)];
    dp1[2 + 20*(2)] = dp1xx * (kx[i] * vf / wd) * (kx[i] * vf / wd);
    dp1[2 + 20*(3)] = dp1xx * kx[i] * ky[i] * (vf / wd) * (vf / wd);
		       /*  changed from presumed erroneous code in original  */
    dp1[2 + 20*(4)] = dptrm * kx[i] * vf / wd;
    dp1[2 + 20*(7)] = p1x3 * dk * vf;
    dp1[2 + 20*(8)] = p1x3 * kx2;
    dp1[2 + 20*(9)] = p1x3 * ky2;
    dp1[2 + 20*(10)] = dp1xx * kx[i] * vf * om1 / wd / wd;
    dp1[2 + 20*(11)] = dptrm * dk * kx[i] * vf *vf / wd;
    dp1[2 + 20*(12)] = -dp1x6 * kx[i] * vf / wd;
    dp1[3 + 20*(2)] = dp1[2 + 20*(3)];
    dp1[4 + 20*(2)] = dp1[2 + 20*(4)];
    dp1[7 + 20*(2)] = dp1[2 + 20*(7)];
    dp1[8 + 20*(2)] = dp1[2 + 20*(8)];
    dp1[9 + 20*(2)] = dp1[2 + 20*(9)];
    dp1[10 + 20*(2)] = dp1[2 + 20*(10)];
    dp1[11 + 20*(2)] = dp1[2 + 20*(11)];
    dp1[12 + 20*(2)] = dp1[2 + 20*(12)];
    dp1[3 + 20*(3)] = dp1xx * (ky[i] * vf / wd) * (ky[i] * vf / wd);
    dp1[3 + 20*(4)] = dptrm * ky[i] * vf / wd;
    dp1[3 + 20*(7)] = p1x4 * dk * vf;
    dp1[3 + 20*(8)] = p1x4 * kx2;
    dp1[3 + 20*(9)] = p1x4 * ky2;
    dp1[3 + 20*(10)] = dp1xx * ky[i] * vf * om1 / wd / wd;
    dp1[3 + 20*(11)] = dptrm * dk * ky[i] * vf * vf / wd;
    dp1[3 + 20*(12)] = -dp1x6 * ky[i] * vf / wd;
    dp1[4 + 20*(3)] = dp1[3 + 20*(4)];
    dp1[7 + 20*(3)] = dp1[3 + 20*(7)];
    dp1[8 + 20*(3)] = dp1[3 + 20*(8)];
    dp1[9 + 20*(3)] = dp1[3 + 20*(9)];
    dp1[10 + 20*(3)] = dp1[3 + 20*(10)];
    dp1[11 + 20*(3)] = dp1[3 + 20*(11)];
    dp1[12 + 20*(3)] = dp1[3 + 20*(12)];
    dp1[4 + 20*(4)] = (dptrm + p1x0) * o1 / wd;
    dp1[4 + 20*(7)] = p1x5 * dk * vf;
    dp1[4 + 20*(8)] = p1x5 * kx2;
    dp1[4 + 20*(9)] = p1x5 * ky2;
    dp1[4 + 20*(10)] = dptrm * om1 / wd;
    dp1[4 + 20*(11)] = (dptrm + p1x0) * dk * o1 * vf / wd;
    dp1[4 + 20*(12)] = -dp1x6 * o1 / wd;
    dp1[7 + 20*(4)] = dp1[4 + 20*(7)];
    dp1[8 + 20*(4)] = dp1[4 + 20*(8)];
    dp1[9 + 20*(4)] = dp1[4 + 20*(9)];
    dp1[10 + 20*(4)] = dp1[4 + 20*(10)];
    dp1[11 + 20*(4)] = dp1[4 + 20*(11)];
    dp1[12 + 20*(4)] = dp1[4 + 20*(12)];
    dp1[7 + 20*(7)] = p1 * dk * dk * vf * vf;
    dp1[7 + 20*(8)] = p1 * dk * kx2;
    dp1[7 + 20*(9)] = p1 * dk * ky2;
    dp1[7 + 20*(10)] = p1x11 * dk * vf;
    dp1[7 + 20*(11)] = p1x12 * dk * vf;
    dp1[7 + 20*(12)] = p1x13 * dk * vf;
    dp1[8 + 20*(7)] = dp1[7 + 20*(8)];
    dp1[9 + 20*(7)] = dp1[7 + 20*(9)];
    dp1[10 + 20*(7)] = dp1[7 + 20*(10)];
    dp1[11 + 20*(7)] = dp1[7 + 20*(11)];
    dp1[12 + 20*(7)] = dp1[7 + 20*(12)];
    dp1[8 + 20*(8)] = p1 * kx2 * kx2;
    dp1[8 + 20*(9)] = p1 * kx2 * ky2;
    dp1[8 + 20*(10)] = p1x11 * kx2;
    dp1[8 + 20*(11)] = p1x12 * kx2;
    dp1[8 + 20*(12)] = p1x13 * kx2;
    dp1[9 + 20*(8)] = dp1[8 + 20*(9)];
    dp1[10 + 20*(8)] = dp1[8 + 20*(10)];
    dp1[11 + 20*(8)] = dp1[8 + 20*(11)];
    dp1[12 + 20*(8)] = dp1[8 + 20*(12)];
    dp1[9 + 20*(9)] = p1 * ky2 * ky2;
    dp1[9 + 20*(10)] = p1x11 * ky2;
    dp1[9 + 20*(11)] = p1x12 * ky2;
    dp1[9 + 20*(12)] = p1x13 * ky2;
    dp1[10 + 20*(9)] = dp1[9 + 20*(10)];
    dp1[11 + 20*(9)] = dp1[9 + 20*(11)];
    dp1[12 + 20*(9)] = dp1[9 + 20*(12)];
    dp1[10 + 20*(10)] = dp1xx * (om1 / wd) * (om1 / wd) -
        p1x0 * om1 * logk * vf;
    dp1[10 + 20*(11)] = dptrm * om1 * dk * vf / wd;
    dp1[10 + 20*(12)] = -dp1x6 * om1 / wd;
    dp1[11 + 20*(10)] = dp1[10 + 20*(11)];
    dp1[12 + 20*(10)] = dp1[10 + 20*(12)];
    dp1[11 + 20*(11)] = (dptrm + p1x0) * o1 *
        (dk * vf / wd) * dk * vf;
    dp1[11 + 20*(12)] = -dp1x6 * o1 * dk * vf / wd;
    dp1[12 + 20*(11)] = dp1[11 + 20*(12)];
    dp1[12 + 20*(12)] = 2 * p1 * vf * vf * d1 / bkf;

    pix11 = exp(x[5]) / k3;
    pix12 = exp(x[6]) / k4;
    pi = 1.0 / (p1 + pix11 + pix12);
    pi2 = pi * pi;
    spow = pi * ppow[i];
    *f += -log (pi) + spow;

    p1x1 *= pi;
    p1x2 *= pi;
    p1x3 *= pi;
    p1x4 *= pi;
    p1x5 *= pi;
    p1x8 *= pi;
    p1x9 *= pi;
    p1x10 *= pi;
    p1x11 *= pi;
    p1x12 *= pi;
    pix11 *= pi;
    pix12 *= pi;
    pow2m1 = spow + spow - 1.0;
    g[0] += p1x1 - spow * p1x1;
    g[1] += p1x2 - spow * p1x2;
    g[2] += p1x3 - spow * p1x3;
    g[3] += p1x4 - spow * p1x4;
    g[4] += p1x5 - spow * p1x5;
    g[5] += pix11 - spow * pix11;
    g[6] += pix12 - spow * pix12;
    g[7] += p1x8 - spow * p1x8;
    g[8] += p1x9 - spow * p1x9;
    g[9] += p1x10 - spow * p1x10;
    g[10] += p1x11 - spow * p1x11;
    g[11] += p1x12 - spow * p1x12;
    g[12] += p1x13 - spow * p1x13;

    for (j = 0; j < n; j++) {
      int rows = 20 * j;
      for (col = 0; col < n; col++) {
        int cols = 20 * col;
	if ((j != 6) && (j != 5) && (col != 6) && (col != 5))
	  dd[j + n*col] += dp1[j + cols] * pi - dp1[rows] * dp1[cols] * pi2 -
	      spow * dp1[j + cols] * pi +
	      2 * spow * dp1[rows] * dp1[cols] * pi2;
      }
    }

    dd[5 + n*5] += pix11 * (1 + spow) + pix11 * pix11 * (2 * spow - 1);
    dd[6 + n*(6)] += pix12 * (1 + spow) + pix12 * pix12 * (2 * spow - 1);
    scale = pix11 * pow2m1;
/*
    for (p = 0; p < 13; p++) {
      if (p == 5) continue;
      dd[5 + n*p] += scale * dpdx[p];
    }
*/
    dd[5 + n*(0)] += scale * p1x1;
    dd[5 + n*(1)] += scale * p1x2;
    dd[5 + n*(2)] += scale * p1x3;
    dd[5 + n*(3)] += scale * p1x4;
    dd[5 + n*(4)] += scale * p1x5;
    dd[5 + n*(6)] += scale * pix12;
    dd[5 + n*(7)] += scale * p1x8;
    dd[5 + n*(8)] += scale * p1x9;
    dd[5 + n*(9)] += scale * p1x10;
    dd[5 + n*(10)] += scale * p1x11;
    dd[5 + n*(11)] += scale * p1x12;
    dd[5 + n*(12)] += scale * p1x13;
    row = 5 * n;
    for (p = 0; p < 13; p++) {
      if (p == 5) continue;
      dd[p + row] = dd[5 + p*n];
    }
    scale = pix12 * pow2m1;
/*
    for (p = 0; p < 13; p++) {
      if (p == 5) continue;
      dd[6 + n*p] += scale * dpdx[p];
    }
*/
    dd[6 + n*(0)] += scale * p1x1;
    dd[6 + n*(1)] += scale * p1x2;
    dd[6 + n*(2)] += scale * p1x3;
    dd[6 + n*(3)] += scale * p1x4;
    dd[6 + n*(4)] += scale * p1x5;
    dd[6 + n*(7)] += scale * p1x8;
    dd[6 + n*(8)] += scale * p1x9;
    dd[6 + n*(9)] += scale * p1x10;
    dd[6 + n*(10)] += scale * p1x11;
    dd[6 + n*(11)] += scale * p1x12;
    dd[6 + n*(12)] += scale * p1x13;
    row = 6 * n;
    for (p = 0; p < 13; p++) {
      if (p == 5) continue;
      dd[p + row] = dd[6 +p*n];
    }
  }
}

/*
static double flnm (void (*fcn)(int, double *, double *, double *), 
    double x, double *df, double *v, double *x0, int n) {
*/
static double flnm (double x, double *df, double *v, double *x0, int n,
    double *rkx, double *rky, double *k2, double *omeg, double *ppow,
    double pk, double *fm, double rki, int npix, int *flnm_ct) {

  int i, n2;
  double flnm_ret;

  (*flnm_ct)++;
  n2 = n + n;
  for (i = 0; i < n; i++) v[n+i] = x0[i] + x * v[i];
/*
  (*fcn) (n, &v[n], &flnm_ret, &v[n2]);
*/
  fcn1 (n, &v[n], &flnm_ret, &v[n2], rkx, rky, k2, omeg, ppow, pk, fm, rki, npix);
  *df = 0.0;
  for (i = 0; i < n; i++) *df += v[i] * v[n2+i];
  return flnm_ret;
}

static int linmin (double *x1, double *x2, double *f1, double *df1,
    double reps, double aeps,
/*
    void (*f)(int, double *, double *, double *),
*/
    double *v, double *xi, int n, double *rkx, double *rky,
    double *k2, double *omeg, double *ppow, double pk, double *fm, double rki,
    int npix, int nit, int *flnm_ct, int *linmin_iter) {
  double dfa, df0, df12, df2, dx, dx1, fa, fc, f0, f2, r, r1, xa, x0;
  int i, qb, status;

  double rho = 0.01, sigma = 0.1;
  double t1 = 9.0, t2 = 0.9, t3 = 0.5;
/*
  int nit = 15;
*/

  status = 0;
  qb = 0;

/*
  f2 = flnm ((*f), *x2, &df2, v, xi, n);
*/
  f2 = flnm (*x2, &df2, v, xi, n, rkx, rky, k2, omeg, ppow, pk, fm, rki, npix,
      flnm_ct);

  dx1 = *x2 - *x1;
  f0 = *f1;
  df0 = *df1;
  x0 = *x1;

  for (i = 0; i < nit; i++) {
    *linmin_iter = i + 1;
    fc = f0 + df0 * rho * (*x2 - x0);
    if ((fabs (df2) <= -sigma * df0) && (f2 <= fc)) {
      *x1 = *x2;
      *f1 = f2;
      *df1 = df2;
      return status;
    }
    if (!qb)
      if ((f2 > fc) || (f2 > *f1) || (df2 >= 0.0)) qb = 1;
    df12 = (f2 - *f1) / dx1;
    r = 2 * df2 + *df1 - 3 * df12;
    r1 = 3 * df12 - df2 - *df1;
    r1 *= r1;
    r1 -= *df1 * df2;
    dx = 0.0;
    if (r1 > 0.0) {
      r1 = sqrt (r1);
      if (r < 0.0) r1 *= -1;
      r += r1;
      dx = -df2 * dx1 / r;
    } else {
      double r2 = 2 * (df12 - df2);
      if (r2 != 0.0) dx = dx1 * df2 / r2;
    }

    if (qb) {
      if (dx < (-t2 * dx1)) dx = -t2 * dx1;
      if (dx > (-t3 * dx1)) dx = -t3 * dx1;
      xa = *x2 + dx;
/*
      fa = flnm ((*f), xa, &dfa, v, xi, n);
*/
      fa = flnm (xa, &dfa, v, xi, n, rkx, rky, k2, omeg, ppow, pk, fm, rki,
          npix, flnm_ct);
      fc = f0 + df0 * rho * (xa - x0);
      if ((fabs (dfa) <= (-sigma * df0)) && (fa <= fc)) {
        *x1 = xa;
	*f1 = fa;
	*df1 = dfa;
	return status;
      } else if ((fa > fc) || (fa >= *f1) || (dfa > 0.0)) {
        *x2 = xa;
	f2 = fa;
	df2 = dfa;
      } else {
        *x1 = xa;
	*f1 = fa;
	*df1 = dfa;
      }
      dx1 = *x2 - *x1;
      if ((fabs (dx1) < (reps * fabs (*x2))) || (fabs (dx1) < aeps)) {
	status = 135;
	if (f2 <= f0) {
	  *x1 = *x2;
	  *f1 = f2;
	  *df1 = df2;
	  status = 33;
	}
	return status;
      }
    } else {
      if ((dx < (*x2 - *x1)) && (dx > 0.0)) dx = *x2 - *x1;
      if ((dx > t1 * (*x2 - *x1)) || (dx <= 0.0)) dx = t1 * (*x2 - *x1);
      *x1 = *x2;
      *x2 += dx;
      dx1 = dx;
      *f1 = f2;
      *df1 = df2;
/*
      f2 = flnm (f, *x2, &df2, v, xi, n);
*/
      f2 = flnm (*x2, &df2, v, xi, n, rkx, rky, k2, omeg, ppow, pk, fm, rki,
          npix, flnm_ct);
    }
  }
  if (f2 <= f0) {
    *f1 = f2;
    *x1 = *x2;
    *df1 = df2;
    return 22;
  } else if ((*f1 <= f0) && (*x1 != x0)) {
    return 22;
  }
  return 133;
}

static int bfgs (int n, double *x, double *f, double *g, double *h,
    double reps, double aeps, double *wk,
/*
    void (*fcn)(int, double *, double *, double *)
*/
    double *rkx, double *rky, double *k2, double *omeg, double *ppow, double pk,
    double *fm, double rki, int npix, int nit, int linmin_nit, int *flnm_ct,
    int *bfgs_iter, int *linmin_iter) {
/*
 *  Implements (?) the BFGS (Broyden, Fletcher, Goldfarb, Shanno)
 *    algorithm for updating the Hessian matrix
 *  The function fcn() is invoked to form the log of the likelihood function
 *    f and fill the array g[i] for the parameter set x[n] and the data.
 *  The work array wk is internal and could be mallocd to 3n.
 */

  double df, dg, f1, ghg, gi, h1, h2, r1, x1, x2;
  int i, ier, it, j, n2, qc, status = 0;
/*
  int nit = 200;
*/

  if (n < 1)
    return 111;
/*
  (*fcn) (n, x, f, g);
*/
  fcn1 (n, x, f, g, rkx, rky, k2, omeg, ppow, pk, fm, rki, npix);

  df = fabs (*f);
  if (df == 0.0) df = 1.0;
  n2 = n + n;
  h2 = 1.0;
  for (it = 0; it < nit; it++) {
    double df1 = 0.0;
    *bfgs_iter = it + 1;
    h1 = h2;
    h2 = 0.0;
    for (i = 0; i < n; i++) {
      wk[i] = 0.0;
      for (j = 0; j < n; j++) {
	h2 += fabs (h[i + n*j]);
	wk[i] -= h[i + n*j] * g[j];
      }
      df1 += wk[i] * g[i];
    }
    if (df1 == 0.0) {
      if (fabs (h2 / h1) > 1.3) status = 44;
      return status;
    }

    x1 = 0.0;
    x2 = -2 * df / df1;
    if (x2 > 1.0) x2 = 1.0;
    f1 = *f;
    if (x2 <= 0.0) x2 = 1.0;
/*
    linmin (&x1, &x2, &f1, &df1, reps, aeps, &ier, (*fcn),
        wk, x, n, linmin_nit);
*/
    ier = linmin (&x1, &x2, &f1, &df1, reps, aeps, wk, x, n, rkx, rky,
        k2, omeg, ppow, pk, fm, rki, npix, linmin_nit, flnm_ct, linmin_iter);
    if (ier > 0) status = ier;
    if (ier > 100) return status;

    qc = 1;
    for (i = 0; i < n; i++) {
      x[i] = wk[n+i];
      wk[n+i] = x1 * wk[i];
      if ((fabs (wk[n+i]) > aeps) &&
          (fabs (wk[n+i]) > reps * fabs (x[i]))) qc = 0;
      gi = wk[n2+i];
      wk[n2+i] -= g[i];
      g[i] = gi;
    }
    df = *f - f1;
    *f = f1;
    if (qc) {
      if (fabs (h2 / h1) > 1.3) status = 44;
      return status;
    }

    for (i = 0; i < n; i++) {
      wk[i] = 0.0;
      for (j = 0; j < n; j++) {
	wk[i] += h[i + n*j] * wk[n2+j];
      }
    }
    ghg = dg = 0.0;
    for (i = 0; i < n; i++) {
      dg += wk[n2+i] * wk[n+i];
      ghg += wk[i] * wk[n2+i];
    }
    r1 = (1.0 + (ghg / dg)) / dg;
    for (j = 0; j < n; j++)
      for (i = 0; i < n; i++)
        h[i + n*j] += r1 * wk[n+i] * wk[n+j] -
	    (wk[n+i] * wk[j] + wk[i] * wk[n+j]) / dg;
  }
  return 122;
}

#define RSUNMM		(696.00)

static int get_scaling_values (DRMS_Record_t *rec, double *dnu, double *dkx,
    double *dky) {
  double rsun;
  double key_dbl;
  int status;
  char *key_str;
  int need_dkx = 0, need_dky = 0, need_dnu = 0;
					    /*  do a little sanity checking  */
  key_str = drms_getkey_string (rec, "WCSNAME", &status);
  if (key_str) {
    if (strcasecmp (key_str, "k_x/k_y/omega")) {
      fprintf (stderr, "Warning: unknown WCS Type %s\n", key_str);
    } else {
      key_str = drms_getkey_string (rec, "CTYPE1", &status);
      if (strcasecmp (key_str, "WAVE-NUM"))
	fprintf (stderr, "Warning: inconsistent value of CTYPE1: %s\n", key_str);
      key_str = drms_getkey_string (rec, "CTYPE2", &status);
      if (strcasecmp (key_str, "WAVE-NUM"))
	fprintf (stderr, "Warning: inconsistent value of CTYPE2: %s\n", key_str);
      key_str = drms_getkey_string (rec, "CTYPE3", &status);
      if (strcasecmp (key_str, "FREQ-ANG"))
	fprintf (stderr, "Warning: inconsistent value of CTYPE3: %s\n", key_str);
    }
  } else {
    fprintf (stderr, "Warning: no WCSTYPE; k_x/k_y/omega assumed\n");
  }
/*  should really be able to parse CUNIT and convert accordingly; for now..  */
  key_str = drms_getkey_string (rec, "CUNIT1", &status);
  if (key_str) {
    if (strcasecmp (key_str, "Mm-1")) {
      fprintf (stderr, "Warning: unparsed value of CUNIT1: %s\n", key_str);
      fprintf (stderr, "  using DELTA_K instead\n");
      need_dkx = 1;
    } else {
      *dkx = drms_getkey_double (rec, "CDELT1", &status);
      if (status || isnan (*dkx)) {
	fprintf (stderr, "Warning: no valid value for CDELT1\n");
	fprintf (stderr, "  using DELTA_K instead\n");
	need_dkx = 1;
      }
    }
  }
  key_str = drms_getkey_string (rec, "CUNIT2", &status);
  if (key_str) {
    if (strcasecmp (key_str, "Mm-1")) {
      fprintf (stderr, "Warning: unparsed value of CUNIT2: %s\n", key_str);
      fprintf (stderr, "  using DELTA_K instead\n");
      need_dky = 1;
    } else {
      *dky = drms_getkey_double (rec, "CDELT2", &status);
      if (status || isnan (*dky)) {
	fprintf (stderr, "Warning: no valid value for CDELT2\n");
	fprintf (stderr, "  using DELTA_K instead\n");
	need_dky = 1;
      }
    }
  }
  key_str = drms_getkey_string (rec, "CUNIT3", &status);
  if (key_str) {
    if (strcasecmp (key_str, "s-1")) {
      fprintf (stderr, "Warning: unparsed value of CUNIT3: %s\n", key_str);
      fprintf (stderr, "  using D_OMEGA or DELTA_NU instead\n");
      need_dnu = 1;
    } else {
      *dnu = drms_getkey_double (rec, "CDELT3", &status);
      if (status || isnan (*dnu)) {
	fprintf (stderr, "Warning: no valid value for CDELT3\n");
	fprintf (stderr, "  using D_OMEGA or DELTA_NU instead\n");
	need_dnu = 1;
      }
    }
  }
  
  if (need_dkx || need_dky) {
    key_dbl = drms_getkey_double (rec, "DELTA_K", &status);
    if (status || isnan (key_dbl)) {
      fprintf (stderr, "Error: no valid value for delta-k\n");
      return 1;
    }
    if (need_dkx) *dkx = key_dbl;
    if (need_dky) *dky = key_dbl;
  }
  rsun = drms_getkey_double (rec, "RSunRef", &status);
  if (status || isnan (rsun)) {
    rsun = RSUNMM;
    fprintf (stderr, "Warning: no valid value for RSunRef; ");
    fprintf (stderr, "using %f\n", rsun);
  }
  *dkx *= rsun;
  *dky *= rsun;

  if (need_dnu) {
    *dnu = drms_getkey_double (rec, "DELTA_NU", &status);
    if (status || isnan (*dnu)) {
      *dnu = drms_getkey_double (rec, "D_OMEGA", &status);
      if (status || isnan (*dnu)) {
	fprintf (stderr, "Error: no valid value for delta-nu\n");
	return 1;
      } else *dnu *= 0.5e6 / M_PI;
    }
  } else *dnu *= 0.5e6 / M_PI;

  return 0;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *pspec;
  FILE *unit11, *unit22, *runlog;
  double *rkl, *rk0, *rku;
  double *freq, *wid;
  double *frl, *rk;
  double *rkk, *rkx, *rky, *k2, *omeg, *ppow;
  double *d2, *h, *hs;
  double par[20], wk[900], par0[20];
  double g[20];
  double aeps, amp, back, bk1, bk2, det, dkx, dky, dnu;
  double f, fr, f0, f00;
  double pk1, reps, kmax, rkw, rl0, sn, wd;
  double pk, fm, rki, xf;
  double dval;
  float *spec;
  int intg[20];
  int c;
  int i, j, jj, k, kx, k0, l;
  int bfgs_iter, bfgs_status, flnm_ct, linmin_iter;
  int iers, fp0, fplo, fpup, fp;
  int nkx, nky, nnu, nfr, nstp;
  int n, nk, nt, nx2, ny2;
  int fit_ct, cvg_ct, cvg_minf, cvg_maxf;
  int npx, npxmax, status, total_clock;
  int rgn, rgnct, segct, isegnum, osegnum, logsegnum, drms_output, dispose;
  char logfile[DRMS_MAXPATHLEN], outfile[DRMS_MAXPATHLEN];
  char source[DRMS_MAXQUERYLEN], recid[DRMS_MAXQUERYLEN];
  char line[1024], module_ident[64];

  double scale[] = {1.0, 1.0, VSCALE, VSCALE, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0e-4};
  double rscale  = 1.0 / RSUN_MM;
  int npar = 13;
  char *hdrstr[] = {
    "          A0", "           c", "          Ux", "          Uy",
    "          w0", "          B1", "          B2", "          A1",
    "          A2", "          A3", "           p", "          w1",
    "           S"};
  char *hdrstr2[] = {
    "        D_A0", "         D_c", "        D_Ux", "        D_Uy",
    "        D_w0", "        D_B1", "        D_B2", "        D_A1",
    "        D_A2", "        D_A3", "         D_p", "        D_w1",
    "         D_S"};

  int nmin = params_get_int (params, "nmin");
  int nmax = params_get_int (params, "nmax");
				     /*  should actually be nmax - nmin + 1  */
  int nct = nmax + 2;
  int lmin = params_get_int (params, "lmin");
  int lmax = params_get_int (params, "lmax");
				     /*  should actually be lmax - lmin + 1  */
  int lct = lmax + 1;
  int fmin = params_get_int (params, "fmin");
  int fmax = params_get_int (params, "fmax");
  int bfgs_nit = params_get_int (params, "bfgsct");
  int linmin_nit = params_get_int (params, "linminct");

  double ux_guess = params_get_double (params, "ux") / VSCALE;
  double uy_guess = params_get_double (params, "uy") / VSCALE;

  char *inds = strdup (params_get_str (params, "in"));
  char *guessfile =  strdup (params_get_str (params, "guessfile"));
  char *outser =  strdup (params_get_str (params, "out"));

  int verbose = params_isflagset (params, "v");
  int extra_reporting = params_isflagset (params, "x");
  int calc_seconds = params_isflagset (params, "2");
  int no_fits = params_isflagset (params, "n");
 
  double p = 0.0;			      /*  uninitialized in original  */

/*
  time_elapsed (0, &total_clock);
*/
  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if ((unit11 = fopen (guessfile, "r")) == NULL) {
    fprintf (stderr, "Error: unable to read frequency guess table \"%s\"\n",
	guessfile);
    return 1;
  }
  d2 = (double *)malloc (npar * npar * sizeof (double));
  h = (double *)malloc (npar * npar * sizeof (double));
  hs = (double *)malloc (npar * npar * sizeof (double));
  freq = (double *)calloc (nct * lct, sizeof (double));
  wid = (double *)calloc (nct * lct, sizeof (double));
		 /*  Read in the initial guess for frequency and line-width  */
  i = j = 0;
	 /*  Read frequencies and widths from table into [nmax*lmax] arrays  */
  while ((c = fgetc (unit11)) != EOF) {
    if (c == '\n') {
      line[j] = '\0';
      if (sscanf (line, "%d %d %lf %lf", &n, &l, &fr, &wd) == 4) {
        if ((n <= nmax + 1) && (n >= nmin) && (l <= lmax) && (l >= lmin)) {
	  if (wd > 600.0) wd = 600.0;
/* this should probably be removed  */
	  if (fr <= freq[n + nct*(l-1)])
	    fr = freq[n + nct*(l-1)] + 0.5;

	  freq[n + l*nct] = fr;
	  wid[n + l*nct] = 0.5 * wd;
	}
      }
      j = 0;
      i++;
    } else line[j++] = (char)c;
  }
  fclose (unit11);
  if (verbose) printf ("%d entries read from guess table\n", i);

  ids = drms_open_records (drms_env, inds, &status);
  if (status) {
    fprintf (stderr, "Error: drms_open_records returned %d for dataset:\n", status);
    fprintf (stderr, "  %s\n", inds);
    return 0;
  }
		     /*  code originally written to process single spectrum 
		      processing of multiple regions should be parallelized  */
  rgnct = ids->n;
  if (rgnct < 1) {
    fprintf (stderr, "Error: no records found in requested dataset:\n");
    fprintf (stderr, "  %s\n", inds);
    return 0;
  }
  irec = ids->records[0];
  segct = irec->segments.num_total;
  if (segct != 1) {
    for (n = 0; n < segct; n++) {
      iseg = drms_segment_lookupnum (irec, n);
      if (iseg->info->naxis == 3) break;
      isegnum = n;
    }
    if (n >= segct) {
      fprintf (stderr, "found no segmemt of rank 3 in input dataset\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
  } else {
    isegnum = 0;
    iseg = drms_segment_lookupnum (irec, isegnum);
  }
  if (!iseg) {
    fprintf (stderr, "Error, could not open data segment\n");
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
					/*  check output data series struct  */
  orec = drms_create_record (drms_env, outser, DRMS_TRANSIENT, &status);
  if (orec) {
    segct = drms_record_numsegments (orec);
    if (segct < 1) {
      fprintf (stderr,
	  "Error: no data segment in output series %s\n", outser);
      drms_close_records (ids, DRMS_FREE_RECORD);
      drms_close_record (orec, DRMS_FREE_RECORD);
      return 1;
    }
    oseg = drms_segment_lookup (orec, "fit.out");
    if (oseg && oseg->info->protocol == DRMS_GENERIC)
      osegnum = oseg->info->segnum;
    else {
    /*  either segment fit.out does not exist, or it has the wrong protocol  */
      for (n = 0; n < segct; n++) {
	oseg = drms_segment_lookupnum (orec, n);
	if (oseg->info->protocol != DRMS_GENERIC) continue;
	osegnum = n;
   	break;
      }
      if (n == segct) {
	fprintf (stderr,
	    "Error: no generic data segment in output series %s\n", outser);
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
      fprintf (stderr, "         writing to segment %s\n", oseg->info->name);
    }
    logsegnum = -1;
    oseg = drms_segment_lookup (orec, "Log");
    if (oseg && oseg->info->protocol == DRMS_GENERIC)
      logsegnum = oseg->info->segnum;
    drms_close_record (orec, DRMS_FREE_RECORD);
    drms_output = 1;
  } else {
	   /*  Can't create in named series, assume it's a filename instead  */
    if (rgnct > 1) {
      fprintf (stderr,
	  "Query produced %d matching records with output to single file;",
	  rgnct);
      fprintf (stderr, " must limit to 1\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    strcpy (outfile, outser);
    drms_output = 0;
  }
  dispose = (no_fits) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
	/*  loop over all input records, creating an output record for each  */
  for (rgn = 0; rgn < rgnct; rgn++) {
    time_elapsed (0, &total_clock);
    irec = ids->records[rgn];
    drms_sprint_rec_query (source, irec);
    iseg = drms_segment_lookupnum (irec, isegnum);
    pspec = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
    if (status) {
      fprintf (stderr, "Warning: could not read segment %d\n", isegnum);
      fprintf (stderr, "       in %s; skipped\n", source);
      continue;
    }
    if (drms_output) {
      char *suffix;
      orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
      oseg = drms_segment_lookupnum (orec, osegnum);
      drms_segment_filename (oseg, outfile);
				  /*  strip off the .generic tag if present  */
      suffix = strstr (outfile, ".generic");
      if (suffix) *suffix = '\0';
      if (logsegnum >= 0) {
        oseg = drms_segment_lookupnum (orec, logsegnum);
        drms_segment_filename (oseg, logfile);
				  /*  strip off the .generic tag if present  */
        suffix = strstr (logfile, ".generic");
        if (suffix) *suffix = '\0';
	runlog = fopen (logfile, "w");
      }
    } else runlog = stderr;
    unit22 = fopen (outfile, "w");
    if (!unit22) {
      fprintf (stderr, "Error: could not open output file %s\n", outfile);
      return 1;
    }

    nkx = pspec->axis[0];
    nky = pspec->axis[1];
    nnu = pspec->axis[2];
		  /*  Warning: this "22" must be at least 2 larger than nct  */
    rkl = (double *)malloc (22 * nnu * sizeof (double));
    rk0 = (double *)malloc (22 * nnu * sizeof (double));
    rku = (double *)malloc (22 * nnu * sizeof (double));

    if (get_scaling_values (irec, &dnu, &dkx, &dky)) {
      fprintf (stderr, "Error: region %d lacks appropriate scaling keywords\n",
	  rgn);
      continue;
    }
    nstp = (int)(rint (8.0 / dnu + 0.5) + 0.1);
    if (nstp < 1) nstp = 1;
    nfr = (int)(rint (100.0 / dnu + 0.5) + 0.1);
    if (nfr < 10) nfr = 10;
/*
    kscale = drms_getkey_double (irec, "MAPSCALE", &status);
    if (status || isnan (kscale)) {
      fprintf (stderr, "Warning: required keyword MAPSCALE missing: assumed 0.125\n");
      kscale = 0.125;
    }
    dlon = kscale * M_PI / 180.0;
    if (verbose) printf ("original image scale = %.8f rad/pixel\n", dlon); 
*/
					       /*  Assumes azimuthal mapping  */
/*
    keyval = drms_getkey_string (irec, "MapProj", &status);
    if (keyval) {
      char *projection = (char *)keyval;
      if (strncmp (projection, "Postel", 6)) {
	fprintf (stderr,
	    "Warning: unknown projection %s; azimuthal symmetry assumed\n",
	    projection);
      }
    dlat = dlon;
    } else {
      fprintf (stderr,
	  "Warning: unknown projection; azimuthal symmetry assumed\n");
      dlat = dlon;
    }
*/
    spec = (float *)pspec->data;
    dval = drms_getkey_double (irec, "LOG_BASE", &status);
    if (!status && isfinite (dval)) {
      double scaled_log_base = dval / exp (1.0);
      int ntot = nkx * nky * nnu;
      for (n = 0; n < ntot; n++) spec[n] = (float)exp (scaled_log_base * spec[n]);
    }
    nt = nnu;
/*
    dkx = 2.0 * M_PI / (dlon * nkx);
    dky = 2.0 * M_PI / (dlat * nky);
*/
			/*  dk is /rad (360 deg / map width), effectively dl  */
    kmax = dkx * (nkx/2 - 3);
    if (verbose) printf ("resolution: dkx: %.2f, dky: %.2f, dnu: %f, kmax = %.2f\n",
	dkx, dky, dnu, kmax);
    fprintf (unit22, "# fit of %s\n", inds);
    fprintf (unit22, "# %s version %s\n", module_name, version_id);
    fprintf (unit22, "# nmin=%d nmax=%d lmin=%d lmax=%d fmin=%d fmax=%d\n",
	nmin, nmax, lmin, lmax, fmin, fmax);
    fprintf (unit22, "# bfgsct=%d linminct=%d\n", bfgs_nit, linmin_nit);
    fprintf (unit22, "#n        l        k        nu      w0");
    for (n = 2; n < 4; n++) {
      fprintf (unit22, "%s%s", hdrstr[n], hdrstr2[n]);
      if (calc_seconds) fprintf (unit22, "            ");
    }
    fprintf (unit22, " qual");
    fprintf (unit22, "  l_guess bfgs           f          fm");
    for (n = 0; n < 2; n++) {
      fprintf (unit22, "%s%s", hdrstr[n], hdrstr2[n]);
      if (calc_seconds) fprintf (unit22, "            ");
    }
    for (n = 4; n < npar; n++) {
      fprintf (unit22, "%s%s", hdrstr[n], hdrstr2[n]);
      if (calc_seconds) fprintf (unit22, "             ");
    }
/*
    fprintf (unit22, "#n        nu  l_guess        l bfgs");
    fprintf (unit22, "           f          fm");
    for (n = 0; n < npar; n++) fprintf (unit22, "%s", hdrstr[n]);
    if (calc_seconds)
      for (n = 0; n < npar; n++) fprintf (unit22, "%s            ", hdrstr2[n]);
    else
      for (n = 0; n < npar; n++) fprintf (unit22, "%s", hdrstr2[n]);
*/
/*
    fprintf (unit22, "         pk          sn");
    if (extra_reporting) fprintf (unit22, " flnm bfgs linm");
*/
    fprintf (unit22, "\n");

    fplo = (int)(fmin / dnu);
    if (fplo <= nfr) fplo = nfr + 1;
    fpup = (int)(fmax / dnu);
    if (fpup > (nt - nfr)) fpup = nt - nfr;
    if (verbose) printf ("fplo = %d fpup = %d\n", fplo, fpup);

    reps = 1.0e-6;
    frl = (double *)malloc (lct * sizeof (double));
    rk = (double *)malloc (lct * sizeof (double));
/*
 *  Interpolate l-values from guess table to target frequencies
 *    for each n, set up [lct] arays of l and f values
 *    interpolated values stored in [22 * nnu] array rk0
 *    based on ??, in arrays rkl and rku
 *
 *  Frequencies are assumed monotone increasing with l for given n
 */
    for (fp = fplo - nfr; fp <= fpup + nfr; fp++) {
      f0 = fp * dnu;
      for (n = 0; n <= nmax; n++) {
	nk = 0;
	for (l = 0; l < lct; l++) {
	  if (freq[n + nct*l] > 0.0) {
	    frl[nk] = freq[n + nct*l];
	    rk[nk] = l;
	    nk++;
	  }
	}
	for (i = 0; i < nk; i++) if (frl[i] > f0) break;
	if (i >= nk) rk0[fp + nnu*n] = lct + n - 1;
	else {
	  xf = (frl[i] - f0) / (frl[i] - frl[i-1]);
	  rk0[fp + nnu*n] = rk[i] - xf * (rk[i] - rk[i-1]);
	}
	if (f0 > frl[nk-1]) {
	/* for frequencies beyond tabulated range in guess table, set to max  */
	  rk0[fp + nnu*n] = lct + n - 1;
	}
      }
/*
 *  Set target bounds in [22 * nnu] arrays rkl and rku
 *  For n > 0 ridges, boundaries are 1/3 of distance to adjoining ridges
 *  For n = 0, lower bound is as above, upper bound is equidistant from target
 */
      rkl[fp] = (2 * rk0[fp] + rk0[fp + nnu]) / 3.0;
      rku[fp] = rk0[fp] + (rk0[fp] - rkl[fp]);
      if (rku[fp] > kmax) rku[fp] = kmax;
      for (n = 1; n < 18; n++) {
	rku[fp + nnu*n] = (2 * rk0[fp + nnu*n] + rk0[fp + nnu*(n-1)]) / 3.0;
	if (rku[fp + nnu*n] > kmax) rku[fp + nnu*n] = kmax;
	rkl[fp + nnu*n] = (2 * rk0[fp + nnu*n] + rk0[fp + nnu*(n+1)]) / 3.0;
      }
      if (verbose) {
	printf ("%8.3f", f0);
	for (j = nmin; j <= nmax; j++) printf ("%11.4e", rk0[fp + nnu*j]);
	printf ("\n");
      }
    }
    free (frl);
    free (rk);

    nx2 = nkx / 2 + 1;
    ny2 = nky / 2 + 1;
    npxmax = nkx * nky * (2 * nfr + 1);
    rkk = (double *)malloc (npxmax * sizeof (double));
    rkx = (double *)malloc (npxmax * sizeof (double));
    rky = (double *)malloc (npxmax * sizeof (double));
    k2 = (double *)malloc (npxmax * sizeof (double));
    omeg = (double *)malloc (npxmax * sizeof (double));
    ppow = (double *)malloc (npxmax * sizeof (double));
    for (n = nmin; n <= nmax; n++) {
      int row = nnu * n;
      int reject_ux = 0, reject_bfgs =0, reject_k = 0, reject_fm = 0,
          reject_sn = 0, reject_l = 0, reject_tot = 0;

      pk = (n) ? 0.75 : 0.5;
      iers = 200;
      f00 = 6000.0;
      cvg_ct = 0;
      cvg_minf = fpup;
      cvg_maxf = fplo;

      if (no_fits) {
	int expfct = fpup - fplo + 1;
	printf
	    ("n = %2d, row = %d, flo = %.2f fhi = %.2f fct = %d, k0 = %.2f, ku = %.2f\n",
	    n, row, fplo * dnu, fpup * dnu, expfct,  rk0[fp0 + row],
	    rku[fp0 + row]);
      }
      for (fp0 = fplo, fit_ct = 0; fp0 <= fpup; fp0 += nstp, fit_ct++) {
	double powmax = 0.0;

	if (no_fits) continue;
	if (iers == 300) continue;
					              /*  Changed 2004.03.30  */
	f0 = fp0 * dnu;
	k0 = (int)rk0[fp0 + row];
	rki = rk0[fp0 + row];
	kx = (int)(rk0[fp0 + row] / dkx) + nx2;
	rkw = rk0[fp0 + row] - rkl[fp0 + row];
	if (rkw > (rku[fp0 + row] - rk0[fp0 + row]))
          rkw = rku[fp0 + row] - rk0[fp0 + row];
	if (verbose) {
          printf ("region %.3f %.3f %.3f %.3f\n", rkl[fp0 + row],
	      rku[fp0 + row], f0, rk0[fp0 + row]);
          printf ("%d %d\n", kx, fp0);
          printf ("%.4f %.4f %.4f %.4f\n",
              rkl[fp0 + row] / dkx + nx2, rku[fp0 + row] / dkx + nx2,
	      nx2 - rkl[fp0 + row] / dkx, nx2 - rku[fp0 + row] / dkx);
	}
	if (kx > (nkx - 3)) {
	  if (verbose) printf ("kx = %d > %d; skipped\n", kx, nkx - 3);
	  continue;
	}
	if ((k0 < lmin) || (k0 > lmax)) {
	  if (verbose) printf ("k0 = %d outside [%d, %d]; skipped\n", k0,
	      lmin, lmax);
	  continue;
	}
	if (freq[n + nct*(k0-1)] > 0.0)
          p = (freq[n + nct*k0] - freq[n + nct*(k0-1)]) * k0 /
	      freq[n + nct*(k0-1)];
	if (verbose) printf ("p = %f\n", p);
	if ((p > 0.3) && (p < 0.5) && (k0 < 495)) pk = p;
	npx = 0;
/*
      oml = f0 - 3.6 * dnu;
      omu = f0 + 3.6 * dnu;
*/
	for (i = fp0 - nfr; i <= fp0 + nfr; i++) {
					              /*  Changed 2007.08.10  */
	  double om = i * dnu;
	  double klow = rkl[i + row];
	  double khigh = rku[i + row];
	  for (j = 0; j < nkx; j++) {
	    double rlx = (j + 1 - nx2) * dkx;
	    double rlx2 = rlx * rlx;
					              /*  Changed 2007.08.10  */
	    int powrow = j + nkx*nky*i;
	    for (k = 0; k < nky; k++) {
	      double rly = (k + 1 - ny2) * dky;
	      double kk = sqrt (rlx2 + rly * rly);
	      if ((kk > klow) && (kk < khigh)) {
		rkk[npx] = kk;
		rkx[npx] = rlx;
		rky[npx] = rly;
		k2[npx] = rlx2 + rly * rly;
		omeg[npx] = om;
		ppow[npx] = spec[powrow + nky*k];
		if (ppow[npx] > powmax) powmax = ppow[npx];
		npx++;
	      }
	    }
	  }
	}

	if ((npx < 20) || (powmax <= 0.0)) {
	  if (verbose)
	    printf ("npx = %d < 20 || powmax = %.1f <= 0.0; skipped\n", npx,
		powmax);
	  continue;
	}
						    /*  The initial guesses  */
/*  Parameter correspondences (notation of Basu et al. ApJ 512, 458 and
 *    Basu & Antia, Solar Phys 192, 469):
 *    par[0]	A_0		7	A_1
 *	1	c		8	A_2
 *	2	U_x		9	A_3
 *	3	U_y		10	p
 *	4	w_0		11	w_1
 *	5	B_1		12	S
 *	6	B_2
 */
	if ((fp0 == fplo) || (iers > 100)) {
	  par[1] = f0 / exp (pk * log (rk0[fp0 + row]));
	  par[2] = ux_guess;
          par[3] = uy_guess;
/*
 * change to NaN?
 */
	  par[0] = log (powmax);
	  par[4] = wid[n + nct*k0];
	  if (par[4] < dnu) par[4] = dnu;
	  if (par[4] > 100.0) par[4] = 100.0;
/*
	par[6] = par[1] * pow (rk0 / rkp, pk);
*/
	  par[7] = par[10] = par[11] = 0.0;
	  par[8] = 10.0 * (rki - 200.0);
	  if (par[8] < 0.0) par[8] = 0.0;
	  par[9] = 6.0 * (rki - 300.0);
	  if (par[9] < 0.0) par[9] = 0.0;
				      /*  par[12] = 0 implies symmetric fit  */
	  par[12] = -200.0;
	  par[5] = log (powmax) + 3 * log (rk0[fp0 + row]) - 5.0;
	  par[6] = par[5] + log (rk0[fp0 + row]);
	} else {
	  for (jj = 0; jj < npar; jj++) par[jj] = par0[jj];
	  back = par[6] - log (rki);
	  if (back < par[5]) back = par[5];
	  par[5] = back;
	  par[6] = back + log (rki);
/*
	par[12] = 500.0;
*/
	}
	reps = 1.e-5;
	aeps = 1.e-7;

	for (i = 0; i < npar; i++) {
	  for (j = 0; j < npar; j++) h[j + npar*i] = 0.0;
	  h[i + npar*i] = 1.0;
	}
	if (verbose)
          printf ("%d %f %f %f %f\n", npx, par[1], par[4], par[5], par[6]);
	flnm_ct = 0;
	bfgs_status = bfgs (npar, par, &f, g, h, reps, aeps, wk,
/*
          (*fcn1),
*/
            rkx, rky, k2, omeg, ppow, pk, &fm, rki, npx, bfgs_nit, linmin_nit,
	    &flnm_ct, &bfgs_iter, &linmin_iter);
	if (verbose) {
          printf ("%5d %5d %12.4e\n", bfgs_status, flnm_ct, f);
          for (i = 0; i < npar; i++) {
            printf ("%12.4e", par[i]);
	    if (i%5 == 4) printf ("\n");
          }
          printf ("\n");
          for (i = 0; i < npar; i++) printf
	      ("%12.4e", sqrt (scale[i] * scale[i] * fabs (h[i + npar*i])));
          printf ("\n");
	}
	fcn1 (npar, par, &f, g, rkx, rky, k2, omeg, ppow, pk, &fm, rki, npx);
	pk1 = pk + 1.0e-4 * par[10];
	fm /= npx - npar;
	rl0 = f0 / par[1];
	if (rl0 > 0.0) rl0 = exp (log (rl0) / pk1);
	amp = par[0];
	bk1 = par[5] - 3 * log (fabs (rl0));
	bk2 = par[6] - 4 * log (fabs (rl0));
	sn = amp - ((bk1 > bk2) ? bk1 : bk2);
	if (verbose) printf ("merit %12.4e %.4f\n", fm, sn);

	if (calc_seconds) {
          d2fcn (npar, par, &f, g, d2, rkk, rkx, rky, k2, omeg, ppow, pk, rki,
	      npx);
          for (i = 0; i < npar; i++) {
	    for (j = 0; j < npar; j++) hs[j + npar*i] = 0.0;
	    hs[i + npar*i] = 1.0;
          }
	  if (gauelm (npar, npar, d2, hs, &det, intg, npar, 0) > 0) return 111;
          if (verbose) {
            for (i = 0; i < npar; i++) printf ("%12.5e",
	        sqrt (scale[i] * scale[i] * fabs (hs[i + npar*i])));
            printf ("\n");
          }
	}

	if ((f0 - f00) > 120) iers = 200;
	if (((f0 - f00) > 220) && (f0 > 3500) && (n <= 1)) iers = 300;
	if (((f0 - f00) > 300) && (f0 > 4000) && (n >= 2)) iers = 300;
/*
	if ((fabs (par[2]) < 500.0) && (bfgs_status < 100) &&
*/
        if ((bfgs_status < 100) &&
            (fabs (rl0 - rki) < 50.0) && (fm < 26.5) && (sn > -2.0) &&
	    (sn < 45.0) && (rl0 > 40.0)) {
	  if (!cvg_ct) cvg_minf = fp0;
	  cvg_maxf = fp0;
	  cvg_ct++;
/*
    fprintf (unit22, "#n         l      k       nu      w0");
    fprintf (unit22, "           f          fm");
    for (n = 2; n < 4; n++)
      fprintf (unit22, "%s  D_%s", hdrstr[n], hdrstr[n]);
    fprintf (unit22, " l_guess bfgs");
    for (n = 0; n < 2; n++)
      fprintf (unit22, "%s  D_%s", hdrstr[n], hdrstr[n]);
    for (n = 4; n < npar; n++) {
      fprintf (unit22, "%s  D_%s", hdrstr[n], hdrstr[n]);
      if (calc_seconds) fprintf (unit22, "            ");
    }
*/
	  fprintf (unit22, "%2d %8.3f %8.5f %9.3f %7.3f", n, rl0, rl0 * rscale,
	      f0, par[4]*scale[4]);
	  for (i=2; i<4; i++) {
	    fprintf (unit22, "%12.4e%12.4e", par[i] * scale[i],
		sqrt (scale[i] * scale[i] * fabs (h[i + npar*i])));
	    if (calc_seconds) fprintf (unit22, "%12.4e",
	        sqrt (scale[i] * scale[i] * fabs (hs[i + npar*i])));
	  }
/*
   Actual quality flag needs to go here
*/
	  fprintf (unit22, "%5i", -1);
	  fprintf (unit22, " %8.3f %4d%12.4e%12.4e", rki, bfgs_status, f, fm);
	  for (i=0; i<2; i++) {
	    fprintf (unit22, "%12.4e%12.4e", par[i] * scale[i],
		sqrt (scale[i] * scale[i] * fabs (h[i + npar*i])));
	    if (calc_seconds) fprintf (unit22, "%12.4e",
	        sqrt (scale[i] * scale[i] * fabs (hs[i + npar*i])));
	  }
	  for (i=4; i<npar; i++) {
	    fprintf (unit22, "%12.4e%12.4e", par[i] * scale[i],
		sqrt (scale[i] * scale[i] * fabs (h[i + npar*i])));
	    if (calc_seconds) fprintf (unit22, "%12.4e",
	        sqrt (scale[i] * scale[i] * fabs (hs[i + npar*i])));
	  }
/*
	  fprintf (unit22, "%2d%10.3f%9.3f%9.3f%5d", n, f0, rki, rl0,
		bfgs_status);
	  fprintf (unit22, "%12.4e%12.4e", f, fm);
*/
/*
 *  scaling of asymmetry parameter, calibration of u_x, u_y parameters to m/s
 *    added 01.11.05
 */
/*
	  for (i=0; i<npar; i++)
	    fprintf (unit22, "%12.4e", par[i] * scale[i]);
	  for (i=0; i<npar; i++) {
	    fprintf (unit22, "%12.4e",
		sqrt (scale[i] * scale[i] * fabs (h[i + npar*i])));
	    if (calc_seconds) fprintf (unit22, "%12.4e",
	      sqrt (scale[i] * scale[i] * fabs (hs[i + npar*i])));
	  }
*/
/*
	  fprintf (unit22, "%12.4e%12.4e\n", pk, sn);
	  if (extra_reporting)
	    fprintf (unit22, "%5d%5d%5d", flnm_ct, bfgs_iter, linmin_iter);
*/
	  fprintf (unit22, "\n");
	  iers = 0;
	  for (j = 0; j < npar; j++) par0[j] = par[j];
	  f00 = f0;
	} else {
	  if (fabs (par[2] - ux_guess) >= 500.0) reject_ux++;
	  if (bfgs_status >= 100) reject_bfgs++;
	  if (fabs (rl0 - rki) >= 50.0) reject_k++;
	  if (fm >= 26.5) reject_fm++;
	  if (sn <= -2.0 || sn >= 45.0) reject_sn++;
	  if (rl0 <= 40.0) reject_l++;
	  reject_tot++;
	}
	if ((fit_ct % 10) == 9) time_elapsed (fit_ct, &total_clock); 
      }
      fprintf (runlog, "End iteration n = %d; time = %s\n", n,
          time_elapsed (0, &total_clock));
      fprintf (runlog, "  %d converged", cvg_ct);
      if (cvg_ct > 1)
	fprintf (runlog, ", first @ %.0f, last @ %.0f", cvg_minf * dnu, cvg_maxf * dnu);
      else if (cvg_ct) fprintf (runlog, ", @ %.0f", cvg_minf * dnu);
      fprintf (runlog, "\n");
      fprintf (runlog,
          "  rejected: %d (ux); %d (bfgs); %d (k); %d (fm); %d (sn); %d (l); %d (total)\n",
          reject_ux, reject_bfgs, reject_k, reject_fm, reject_sn, reject_l,
	  reject_tot);
    }
    fclose (unit22);

    if (drms_output) {
      int kstat = 0;
      int keyct = sizeof (propagate) / sizeof (char *);
      char *key_str;
      double key_dbl;
      int key_int, crcl_known = 1;
      TIME key_tim;
				    /*  copy designated keywords from input  */
      kstat += propagate_keys (orec, irec, propagate, keyct);
	     /*  if necessary, construct CarrTime from CR:CL, or vice-versa  */
      key_str = drms_getkey_string (irec, "CarrTime", &status);
      if (status) {
        key_int = drms_getkey_int (irec, "CarrRot", &status);
	if (status) crcl_known = 0;
	else {
	  key_dbl = drms_getkey_double (irec, "CMLon", &status);
	  if (status) crcl_known = 0;
	}
	if (crcl_known) {
	  char CarrTime[9];
	  snprintf (CarrTime, 9, "%d:%03.0f", key_int, key_dbl);
	  kstat += check_and_set_key_str (orec, "CarrTime", CarrTime);
	}
      } else {
	key_int = drms_getkey_int (irec, "CarrRot", &status);
	if (status) {
	  sscanf (key_str, "%d:%lf", &key_int, &key_dbl);
	  kstat += check_and_set_key_int (orec, "CarrRot", key_int);
	  kstat += check_and_set_key_double (orec, "CMLon", key_dbl);
	}
      }

      kstat += check_and_set_key_str (orec, "Module", module_ident);
      kstat += check_and_set_key_str (orec, "Source", source);
      kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
      kstat += check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
      if (kstat) {
	drms_sprint_rec_query (recid, orec);
	fprintf (stderr, "Error writing key values(s) to record %s\n", recid);
	fprintf (stderr, "      series may not have appropriate structure\n");
     	drms_close_record (orec, DRMS_FREE_RECORD);
        continue;
      }
      gethostname (line, 1024);
      fprintf (runlog, "Total wall clock time = %d sec on %s\n", total_clock, line);
      if (runlog != stderr) fclose (runlog);
      drms_close_record (orec, dispose);
    }
  }
  drms_close_records (ids, DRMS_FREE_RECORD);
  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  v 0.7  09.07.14		ported from SOI ringfit v 0.86
 *    0.8  09.08.25		fixed scaling of ux_guess and uy_guess and test
 *			on par[2]; trimmed output format slightly; fixed bug in
 *			setting of Dnu from Domega
 *    0.9  09.09.04		output to DRMS as possible, with processing of
 *			multiple records; copying and setting of a few keys;
 *			modified output format to conform to standard
 *	09.11.04		fixed key setting; modified to use and prefer
 *			WCS keywords for frequency scaling, and to no longer
 *			use propagated mapping keys from tracking
 *	09.11.16		corrected rscale
 *    1.0  09.12.09		added numerous keywords to propagation list
 *      10.01.07		removed constraints on values of A0, B1, and
 *      		B2, so that fits are insensitive to power spectrum
 *	      		normalization
 *	10.02.01		fixed treatment of params_get_str()    
 *				fixed setting of Source key and diagnostics and
 *			error recovery for multiple region processing
 *	10.02.21		minor stderr format mods; added Created key
 *	10.03.05		added BLD_VERS key
 *	10.04.22		corrected rscale (again! maybe it's right this
 *			time); changed defaults for nmax, lmin, and lmax to be
 *			more appropriate for HMI data; added quality flag column
 *			to output; removed default values for input dataset and
 *			guessfile; added three more apodization keywords to
 *			propagation list
 *  v 1.0 frozen 2010.04.23
 */	
