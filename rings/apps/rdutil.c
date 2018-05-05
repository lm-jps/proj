/*
 *  rdutil.c
 *
 *  library of utility functions for dealing with ring-diagram fit files
 *
 *    autoweed_vel	automatically weed velocities for inversions
 *    read_fit		reads n, l, f, and d_f from a fit file
 *    read_fit_v	reads n, l, f, d_f, u_i and d_ui from a fit file
 *
 *  The functions in this file are included in module versions as indicated:
 *	function	first included in
 *
 *	autoweed_vel	rdvinv v 0.92
 *	interp_vel	rdvinv v 0.92
 *	read_fit_v	rdvinv v 0.92
 *
 *  Corresponding functions in older unversioned source files are included
 *    in the following modules:
 *
 *	autoweed_vel	rdfitc v 1.2
 *
 *  Additional functions in older unversioned source files to be included
 *    in future are:
 *
 *	autoweed	rdsinv v 0.6
 *	freqdif		rdsinv v 0.6
 *	interp_freq	rdsinv v 0.6
 *	read_fit	rdsinv v 0.6
 *
 *  Functions referenced internally are (or will be):
 *	function	called in
 *
 *	divdif		interp_vel, (freqdif, interp_freq)
 *	line		(autoweed)
 *	llsq		(autoweed)
 *	nearst		divdif
 *	svd		(svdevl)
 *	svdevl		(llsq)
 *
 *  Bugs:
 *    read_fit() and read_fit_v() assume that the frequency error estimate
 *	is in column 5 of the fit table; for fit files produced with versions
 *	of rdfitc below 1.2, this is not the case
 *    The functions read_fit() and read_fit_v() could be merged (using varargs?)
 *    In autoweed_vel(), maxn and tol should be calling arguments, not hard-
 *	coded. Similarly for nmax in divdif(), the array lengths of
 *	n_num and fb in interp_vel(), and possibly the value of aeps in
 *	interp_vel()
 *    range in autoweed_vel() should perhaps be parametrized
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int autoweed_vel (int *n, double *l, double *ux, double *uy, int *mask,
    int npts) {
/*
 *  Returns 'mask' of length npts: mask[i] is 0 for rejected modes.
 *    Ux and Uy are weeded together. Values are rejected when their difference
 *    from the mean exceeds 5 (tol?) std. dev., both the mean and std dev
 *    being calculated over the central 3/5 of the l range for the order
 */
  const int maxn = 8;
  const double tol = 5.0;

  int i, j, offset;
  int n_num[maxn];
  double num;
  double sumx, sumxx, sumy, sumyy, meanx, meany, stdx, stdy;
  double range, uxmin, uxmax, uymin, uymax;
  double lmin, lmax;
	  
  for (i = 0; i < maxn; i++) n_num[i] = 0;
  for (i = 0; i < npts; i++) {
    if (n[i] < maxn) n_num[n[i]]++;
    mask[i] = 0;
  }

  offset = 0;
  for (i = 0; i < maxn; i++) {
    num = 0.0;
    sumx = 0.0;
    sumy = 0.0;
    sumxx = 0.0;
    sumyy = 0.0;
							    /*  get l limits  */
    lmin = lmax = l[offset];
    for (j = offset; j < n_num[i] + offset; j++) {
      if (l[j] < lmin) lmin = l[j];
      if (l[j] > lmax) lmax = l[j];
    }
    range = 0.2 * (lmax - lmin);
    lmin += 2.0 * range;
    lmax -= 2.0 * range;
    for (j = offset; j < n_num[i] + offset; j++) {
      if ((l[j] <= lmax) && (l[j] >= lmin)) {
	sumx += ux[j];
	sumxx += ux[j]*ux[j];
	sumy += uy[j];
	sumyy += uy[j]*uy[j];
	num++;
      }
    }
    meanx = sumx / num;
    meany = sumy / num;
    stdx = num*sumxx - sumx*sumx;
    stdy = num*sumyy - sumy*sumy;
    stdx /= (num * (num - 1));
    stdy /= (num * (num - 1));
    stdx = sqrt (stdx);
    stdy = sqrt (stdy);
    uxmin = meanx - tol * stdx;
    uxmax = meanx + tol * stdx;
    uymin = meany - tol * stdy;
    uymax = meany + tol * stdy;
    
    for (j = offset; j < n_num[i] + offset; j++) {
      if (ux[j] <= uxmax && ux[j] >= uxmin &&
	  uy[j] <= uymax && uy[j] >= uymin) mask[j] = 1;
     }
    offset += n_num[i];
  }
  return 0;
}

int nearst (double xb, double x[], int ntab) {
  int low, igh, mid;

  low=0; igh=ntab-1;
  if ((xb < x[low]) != (xb < x[igh]) ) {
      /*  If the point is within the range of table, locate it by bisection  */
    while (igh-low > 1) {
      mid = (low + igh) / 2;
      if ((xb < x[mid]) == (xb < x[low])) low = mid;
      else igh = mid;
    }
  }
  if (fabs (xb - x[low]) < fabs (xb - x[igh])) return low;
  else return igh;
}

int divdif (double xb, double x[], double f[], int *nuse, int ntab,
    double fb[], double aeps, double *dfb, double *ddfb) {
  int i, j, k, next, in, ip, nit, ier, nmax=10;
  double err, px, dpx, ddpx, xn[11], xd[11];
						  /*  Find the nearest point  */
  next = nearst (xb, x, ntab);
  fb[1] = f[next];
  xd[1] = f[next];
  xn[1] = x[next];
  ier = 0;
  px = 1.0;
				     /*  Initialisation for the derivatives   */
  *dfb = *ddfb = 0.0;
  dpx = ddpx = 0.0;
		    /*  points between in and ip are used for interpolation   */
  ip = in = next;
		   /*  maximum number of points to be used for interpolation  */
  nit = *nuse;
  if (nmax < nit) nit = nmax;
  if (ntab < nit) nit = ntab;
  if (*nuse > nmax || *nuse > ntab) ier = 22;
  if (*nuse < 1) {
    ier = 21;
    nit = 6;
    if (nmax < nit) nit = nmax;
    if (ntab < nit) nit = ntab;
  }
  *nuse = 1;
	  		  /*  calculate successive interpolation polynomials  */
  for(j = 2; j <= nit; ++j) {
				     /*  choose the next nearest point to xb  */
    if (in <= 0 ) {
      ip = ip + 1;
      next = ip;
    } else if (ip >= ntab-1) {
      in = in - 1;
      next = in;
    } else if (fabs (xb - x[ip+1]) < fabs (xb - x[in-1]) ) {
      ip = ip + 1;
      next = ip;
    } else {
      in = in - 1;
      next = in;
    }
				     /*  calculating the divided differences  */
    xd[j] = f[next];
    xn[j] = x[next];
    for(k = j - 1; k >= 1; --k)
      xd[k] = (xd[k+1] - xd[k]) / (xn[j] - xn[k]);
					     /*  calculating the derivatives  */
    ddpx = ddpx * (xb - xn[j-1]) + 2.0 * dpx;
    dpx = dpx * (xb - xn[j-1]) + px;
    *dfb = *dfb + dpx * xd[1];
    *ddfb = *ddfb + ddpx * xd[1];

    px = px * (xb - xn[j-1]);
    err = xd[1] * px;
    fb[j] = fb[j-1] + err;
    *nuse = j;

    if (fabs (err) < aeps) return ier;
  }
  return 23;
}

int interp_vel (int *n, double *l, double *f, double *ef, double *ux, 
      double *eux, double *uy, double *euy, int npts) {
/*
 *   Takes power spectra fit parameters, and interpolates them to integer l.
 *   Data are returned in the input arrays - **data are overwritten!**
 */
  int i, j, nuse, ierr, offset=0;
  int n_num[13];
  double fb[10], ll;
  double dfb, ddfb;
  double *inp_l, *inp_f, *inp_ef, *inp_ux, *inp_uy, *inp_eux, *inp_euy;
  
  double aeps = 1e-7;
				  /*  count number of frequencies for each n  */
		/*  note that it is assumed that all frequencies for each n 
					  are grouped together in the arrays  */
  for (i = 0; i < 13; i++) n_num[i] = 0;
  for (i = 0; i < npts; i++) (n_num[n[i]])++;
  
  inp_l = (double *)malloc (npts * sizeof(double));
  inp_f = (double *)malloc (npts * sizeof(double));
  inp_ef = (double *)malloc (npts * sizeof(double));
  inp_ux = (double *)malloc (npts * sizeof(double));
  inp_uy = (double *)malloc (npts * sizeof(double));
  inp_eux = (double *)malloc (npts * sizeof(double));
  inp_euy = (double *)malloc (npts * sizeof(double));
							     /*  loop over n  */
  for (i=0; i < 10; i++) {
    for (j=0; j<npts; j++) {
      inp_l[j] = 0.0;
      inp_f[j] = 0.0;
      inp_ef[j] = 0.0;
      inp_ux[j] = 0.0;
      inp_uy[j] = 0.0;
      inp_eux[j] = 0.0;
      inp_euy[j] = 0.0;
    }
    for(j=0; j<n_num[i]; j++) {
      inp_l[j] = l[j+offset];
      inp_f[j] = f[j+offset];
      inp_ef[j] = ef[j+offset];
      inp_ux[j] = ux[j+offset];
      inp_uy[j] = uy[j+offset];
      inp_eux[j] = eux[j+offset];
      inp_euy[j] = euy[j+offset];
    }
    for (j = 0; j<n_num[i]; j++) {
      ll = (int)(inp_l[j]+0.5);
      nuse = 3;
      ierr = divdif (ll, inp_l, inp_f, &nuse, n_num[i], fb, aeps, &dfb, &ddfb);
      f[j+offset] = fb[nuse];
      nuse = 3;
      ierr = divdif (ll, inp_l, inp_ef, &nuse, n_num[i], fb, aeps, &dfb, &ddfb);
      ef[j+offset] = fb[nuse];
      nuse = 3;
      ierr = divdif (ll, inp_l, inp_ux, &nuse, n_num[i], fb, aeps, &dfb, &ddfb);
      ux[j+offset] = fb[nuse];
      nuse = 3;
      ierr = divdif (ll, inp_l, inp_uy, &nuse, n_num[i], fb, aeps, &dfb, &ddfb);
      uy[j+offset] = fb[nuse];
      nuse = 3;
      ierr = divdif (ll, inp_l, inp_eux, &nuse, n_num[i], fb, aeps, &dfb, &ddfb);
      eux[j+offset] = fb[nuse];
      nuse = 3;
      ierr = divdif (ll, inp_l, inp_euy, &nuse, n_num[i], fb, aeps, &dfb, &ddfb);
      euy[j+offset] = fb[nuse];
      l[j+offset] = ll;
    }
    offset += n_num[i];
  }
  free (inp_l);
  free (inp_f);
  free (inp_ef);
  free (inp_ux);
  free (inp_uy);
  free (inp_eux);
  free (inp_euy);

  return 0;
}

int read_fit_v (FILE *fpt, int *npts, int **n, double **l, double **f,
    double **ef, double **ux, double **eux, double **uy, double **euy) {
/*
 *  Given an open file pointer, and reads n, l, frequencies, velocity
 *    parameters, and velocity parameter error estimates, and returns
 *    arrays for same of length corresponding to the total number of
 *    tabulated mode fits
 */
  int i, nlines;
  char buffer[8192];
  if (ferror (fpt)) return -1;
  if (feof (fpt)) rewind (fpt);
  
  nlines = 0;
  while (!feof (fpt))     {
    fgets (buffer, 8192, fpt);
    if (buffer[0] != '#' && !feof (fpt)) nlines++;
  }
  if (nlines == 0) return 1;
  *n = (int *)malloc (nlines * sizeof(int));
  *l = (double *)malloc (nlines * sizeof(double));
  *f = (double *)malloc (nlines * sizeof(double));
  *ef = (double *)malloc (nlines * sizeof(double));
  *ux = (double *)malloc (nlines * sizeof(double));
  *eux = (double *)malloc (nlines * sizeof(double));
  *uy = (double *)malloc (nlines * sizeof(double));
  *euy = (double *)malloc (nlines * sizeof(double));
  rewind (fpt);
  
  for (i=0; i<nlines; i++) {
    fgets (buffer, 8192, fpt);
    while (buffer[0] == '#') fgets (buffer, 8192, fpt);
    sscanf (buffer, "%i %lf %*f %lf %lf %lf %lf %lf %lf", &(*n)[i], &(*l)[i],
	&(*f)[i], &(*ef)[i], &(*ux)[i], &(*eux)[i], &(*uy)[i], &(*euy)[i]);
  } 
  *npts = nlines;
  
  return 0;
}

/*
 *  Revision History
 *
 *  v 0.9	First numbered version created 2018.03.12 from old_rdutil.c
 *	Modified autoweed_vel() to initialize all mask values to 0
 *	Modified read_fit_v() to malloc arrays to exact number of mode fits
 *	  (and to not calloc *l)
 *	Modified interp_vel() to free internally malloc'd inp_* arrays
 *	Removed functions autoweed(), freqdif(), interp_freq(), and read_fit()
 *	  and supporting functions referenced only in rdsinv v 0.6
 *	Replaced hard-coded tolerance values in autoweed_vel() with variable
 *	  (which is set to a hard-coded value!)
 *	Removed unused file pointer in interp_vel
 *  v 0.9 frozen 2018.05.02
 */
