/*
 *  C version of subroutine ola_ in ola_xy_v13.f
 *
 *  Bugs
 *    Relies on FORTRAN lapack subs dsytrf, dsyytrs
 *    Derivative estimates from divdif are not used
 *
 *  Revision history is at end of file
 */

/*
 *  Global declaration of arrays whose dimeansions are dependent on only
 *    the (constant) numbers of target radii and/or radial samples in
 *    the kernel file
 */

double *rtrg, *rcgx, *rcgy, *quartx, *quarty, *solnx, *solny, *errx, *erry;
double *aa1x, *aa1y, *fakc;
double *avcx, *avcy, *domx, *domy, *errorx, *errory;
double *avccx, *avccy, *ctil1x, *ctil1y, *ctil2x, *ctil2y;
/*
  wmoden = (int *)malloc (kernrvals * sizeof (int));
  wmodel = (int *)malloc (kernrvals * sizeof (int));
*/

extern void dsytrf_ (char *, int *, double *, int*, int *, double *, int *,
    int *);

extern void dsytrs_ (char *, int *, int *, double *, int *, int *, double *,
    int *, int *);

void fwidth (int np, int num, double *rad, double *avc, double *widthc) {
  double *dx, *y;
  double fb[10];
  double anor, dfb, ddfb, sker, sum;
  double q1 = 0.25, q2 = 0.5, q3 = 0.75;
  double reps = 1.0e-5;
  int i, j, nuse;

  dx = (double *)malloc (np * sizeof (double));
  y = (double *)malloc (np * sizeof (double));

  for (j = 0; j < np - 1; j++) dx[j] = rad[j] - rad[j + 1];
  y[0] = sum = 0.0;
  for (j = 0; j < np - 1; j++) {
    sker = 0.5 * (avc[j] + avc[j+1]);
    sum += sker * dx[j];
    y[j+1] = sum;
  }
  anor = y[np-1];
  for (i = 0; i < np; i++) y[i] /= anor;
  nuse = 4;
  divdif (q1, y, rad, &nuse, np, fb, reps, &dfb, &ddfb);
  widthc[0] = fb[nuse];
  nuse = 4;
  divdif (q2, y, rad, &nuse, np, fb, reps, &dfb, &ddfb);
  widthc[1] = fb[nuse];
  nuse = 4;
  divdif (q3, y, rad, &nuse, np, fb, reps, &dfb, &ddfb);
  widthc[2] = fb[nuse];
  free (y);
  free (dx);
}

int solve (int la, int norder, int cols, double *a, double *b) {
  double *work;
  static int *ipiv = NULL;
  int n, lda, lwork, info, nrhs, ldb;
  static int last = 0, optsize;
  char uplo;

  uplo = 'U';
  n = norder;
  lda = la;
  if (norder != last) {
    ipiv = (int *)realloc (ipiv, norder * sizeof (int));
    lwork = -1;
    work = (double *)malloc (sizeof (double));
    dsytrf_ (&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    if (info) {
      fprintf (stderr, "Error in solve():\n");
      fprintf (stderr, "  initialization call to dsytrf_ returned %d\n", info);
      return info;
    }
    optsize = work[0];
    free (work);
    uplo = 'U';
    n = norder;
    lda = la;
    last = norder;
  }
  lwork = optsize;
  work = (double *)malloc (lwork * sizeof (double));
  dsytrf_ (&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
  if (info) {
    fprintf (stderr, "Error in solve():\n");
    fprintf (stderr, "  call to dsytrf_ returned %d\n", info);
    free (work);
    return info;
  }

  uplo = 'U';
  n = norder;
  nrhs = cols;
  ldb = lda = la;
  dsytrs_ (&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info) {
    fprintf (stderr, "Error in solve():\n");
    fprintf (stderr, "  call to dsytrs_ returned %d\n", info);
  }

  free (work);
  return info;
}

int ola_xy (double *dl_inp, int *n_inp, double *f_inp, double *ux_inp,
    double *euxinp, double *uy_inp, double *euyinp, int nmodes,
    double *rad, int *wmoden, int *wmodel,
    double *weight, int nrpts, int nct, int lct,
    int *moden, int *model, double *krnamp, int krnmodect,
    int *nummodes, double ob, double oe, int numtrg, double mu) {
  double *ov1, *covx, *covy, *suml, *vv1x, *vv1y;
  double freq[lct][nct][2], err[lct][nct][2], fdif[lct][nct];
  double dthcx[3], dthcy[3];
  double sume1, sume2, sumcx, sumcy;
  double sum1, sum2, sum3, sa, sb, sc, term, sumccx, sumccy;
  double sumc1x, sumc2x, sumc1y, sumc2y;
  double sum1x, er1x, s1x, sum1y, er1y, s1y;
  float *ak;
  static float *rd = NULL;
  float rs, rms;
  int krnmode, nmode, np, /* nsp, */ numker, nktmp, /* ntab, */ nmdtmp;
  int nrhs, ierr;
  int i, j, jj, k, l, ni, nik, target;
  int status;
						     /*  initialize matrices  */
  memset (fdif, 0, lct * nct * sizeof (double));
  memset (freq, 0, lct * nct * 2 * sizeof (double));
  memset (err, 0, lct * nct * 2 * sizeof (double));

  ov1 = (double *)calloc (nrpts * nrpts, sizeof (double));
  covx = (double *)calloc (nrpts * nrpts, sizeof (double));
  covy = (double *)calloc (nrpts * nrpts, sizeof (double));
  suml = (double *)calloc (3 * nrpts * nrpts, sizeof (double));

  vv1x = (double *)calloc (numtrg * nrpts, sizeof (double));
  vv1y = (double *)calloc (numtrg * nrpts, sizeof (double));
  memset (ctil1x, 0, nrpts * numtrg * sizeof (double));
  memset (ctil1y, 0, nrpts * numtrg * sizeof (double));
  memset (ctil2x, 0, nrpts * numtrg * sizeof (double));
  memset (ctil2y, 0, nrpts * numtrg * sizeof (double));
						 /*  read in observed errors  */
  for (i = 0; i < nmodes; i++) {
    int ni = n_inp[i];
    int l = dl_inp[i];
    double f = 0.001 * f_inp[i];
							    /*  special hack  */
if (ni >= nct || l >= lct) continue;
    if (f > ob && f < oe) {
      freq[l][ni][0] = ux_inp[i];
      freq[l][ni][1] = uy_inp[i];
      err[l][ni][0] = euxinp[i];
      err[l][ni][1] = euyinp[i];
      fdif[l][ni] = f_inp[i];
    }
  }

  for (krnmode = 0, nmode = 0; krnmode < krnmodect; krnmode++) {
    int ll = model[krnmode];
    int nn = moden[krnmode];
    if (ll >= lct || nn >= nct) continue;
    double dif1 = freq[ll][nn][0];
    double dif2 = freq[ll][nn][1];
    double erri1 = err[ll][nn][0];
    double erri2 = err[ll][nn][1];
    if (erri1 <= 1.0e-11) continue;
    if (erri2 <= 1.0e-11) continue;
    for (j = 0; j < nrpts; j++)
      fakc[j + nrpts * nmode] = krnamp[j + nrpts * krnmode];
    wmoden[nmode] = nn;
    wmodel[nmode] = ll;
    errorx[nmode] = erri1;
    errory[nmode] = erri2;
    domx[nmode] = dif1;
    domy[nmode] = dif2;
    nmode++;
  }
  sume1 = sume2 = 0.0;
  for (i = 0; i < nmode; i++) {
    sume1 += 1.0 / (errorx[i] * errorx[i]);
    sume2 += 1.0 / (errory[i] * errory[i]);
  }
						 /*  set covariance matrices  */
  sumcx = sumcy = 0.0;
  for (j = 0; j < nmode; j++) {
    sumcx += errorx[j] * errorx[j];
    covx[nrpts*j + j] = errorx[j] * errorx[j];
    sumcy += errory[j] * errory[j];
    covy[nrpts*j + j] = errory[j] * errory[j];
  }
  sumcx /= nmode;
  sumcy /= nmode;
					       /*  Now the OLA related stuff  */
						  /*  set the overlap matrix  */
  int nrpsq = nrpts * nrpts;
  for (i = 0; i < nmode; i++) {
    int inrp = i * nrpts;
    for (j = 0; j < nmode; j++) {
      int jnrp = j * nrpts;
      sum2 = sum3 = sa = sb = sc = 0.0;
      for (k = 0; k < nrpts; k++) {
	term = weight[k] * fakc[inrp + k] * fakc[jnrp + k];
	sa += term * rad[k] * rad[k];
	sb += term * rad[k];
	sc += term;
      }
      suml[jnrp + i] = suml[inrp + j] = sa;
      suml[nrpsq + jnrp + i] = suml[nrpsq + inrp + j] = sb;
      suml[2*nrpsq + jnrp + i] = suml[2*nrpsq + inrp + j] = sc;
    }
  }
						 /*  the Lagrange multiplier  */
  int nmnrp = nmode * nrpts;
  for (i = 0, nik = 0; i < nmode; i++) {
    int inrp = i * nrpts;
    sum1 = 0.0;
    for (k = 0; k < nrpts; k++, nik++) sum1 += weight[k] * fakc[nik];
    ov1[i + nmnrp] = ov1[nmode + inrp] = 0.5 * sum1;
    suml[i + nmnrp] = suml[nmode + inrp] = 0.5 * sum1;
  }

  numker = nmode + 1;
  *nummodes = nmode;
					     /*  need loop over X and Y here  */
  for (target = 0; target < numtrg; target++) {
						  /*  set the overlap matrix  */
    double rtrgsq = rtrg[target] * rtrg[target];
    int nrpt2 = nrpts * nrpts;
    for (i = 0; i < nmode; i++) {
      int inrp = i * nrpts;
      for (j = i; j < nmode; j++) {
	sum1 = suml[nrpts*j + i] -
	    2 * rtrg[target] * suml[nrpt2 + nrpts*j + i] +
	    rtrgsq * suml[2*nrpt2 + nrpts*j + i];
	ov1[i + j*nrpts] = ov1[j + inrp] = sum1;
      }
    }
					   /*  set inhomogeneous part of rhs  */
    for (j = 0; j < numtrg; j++)
      vv1x[nmode + j*nrpts] = vv1y[nmode + j*nrpts] = 0.5;
						      /*  set matrix for lhs  */
    for (j = 0; j < nmode; j++) {
      int jnrp = j * nrpts;
      for (i = 0; i < nmode; i++) {
	aa1x[i + jnrp] = ov1[i + jnrp] + mu * covx[jnrp + i] / sumcx;
	aa1y[i + j*nrpts] = ov1[i + jnrp] + mu * covy[jnrp + i] / sumcy;
     }
    }
    for (j = 0; j < numker; j++) {
      int jnrp = j * nrpts;
      for (i = nmode; i < numker; i++) {
	aa1x[i + jnrp] = aa1y[i + jnrp] = ov1[i + jnrp];
	aa1x[j + i*nrpts] = aa1y[j + i*nrpts] = ov1[j + i*nrpts];
      }
      for (k = 0; k < numtrg; k++) {
	ctil1x[j + k*nrpts] = vv1x[j + k*nrpts];
	ctil1y[j + k*nrpts] = vv1y[j + k*nrpts];
      }
    }
						  /*  solving the equations  */
    status = solve (nrpts, numker, 1, aa1x, ctil1x);
    status = solve (nrpts, numker, 1, aa1y, ctil1y);
    jj = 0;
    for (i = 0; i < numker; i++) {
      ctil2x[i + target*nrpts] = ctil1x[i + jj*nrpts];
      ctil2y[i + target*nrpts] = ctil1y[i + jj*nrpts];
    }
					   /*  Calculating averaging kernel  */
    for (j = 0; j < nrpts; j++) {
      sumccx = sumccy = 0.0;
      for (i = 0; i < nmode; i++) {
	sumccx += ctil1x[i + jj*nrpts] * fakc[j + i*nrpts];
	sumccy += ctil1y[i + jj*nrpts] * fakc[j + i*nrpts];
      }
      avcx[j] = sumccx;
      avccx[numtrg*j + target] = sumccx;
      avcy[j] = sumccy;
      avccy[numtrg*j + target] = sumccy;
    }
					 /*  Finding cg of averaging kernel  */
    sumc1x = sumc2x = sumc1y = sumc2y = 0.0;
     for (j = 0; j < nrpts; j++) {
      sumc1x += weight[j] * rad[j] * avcx[j];
      sumc2x += weight[j] * avcx[j];
      sumc1y += weight[j] * rad[j] * avcy[j];
      sumc2y += weight[j] * avcy[j];
    }
    fwidth (nrpts, numtrg, rad, avcx, dthcx);
    fwidth (nrpts, numtrg, rad, avcy, dthcy);
							    /*  the solution  */
    sum1x = er1x = sum1y = er1y = 0.0;
    for (i = 0; i < nmode; i++) {
      sum1x += ctil1x[i] * domx[i];
      er1x += ctil1x[i] * ctil1x[i] * covx[nrpts*i + i];
      sum1y += ctil1y[i] * domy[i];
      er1y += ctil1y[i] * ctil1y[i] * covy[nrpts*i + i];
    }
    s1x = fabs (dthcx[0] - dthcx[2]);
    s1y = fabs (dthcy[0] - dthcy[2]);
    er1x = sqrt (er1x);
    er1y = sqrt (er1y);
    rcgx[target] = sumc1x / sumc2x;
    rcgy[target] = sumc1y / sumc2y;
    for (i = 0; i < 3; i++) {
      quartx[i + 3*target] = dthcx[i];
      quarty[i + 3*target] = dthcy[i];
    }
    solnx[target] = sum1x;
    solny[target] = sum1y;
    errx[target] = er1x;
    erry[target] = er1y;
  }

  free (ov1);
  free (covx);
  free (covy);
  free (suml);
  free (vv1x);
  free (vv1y);
  return 0;
}

/*
 *  Revision History
 *  v 1.2 based on ola_xy_v12.f with modifications:
 *	Removed commented code, unused variables GG, anorml, omax, omin, ndim,
 *	  w, delta, and fb (and assignments of same)
 *	Removed FORTRAN read variables ii, ji
 *	Removed no longer needed variables inum, kk
 *	Changed variable qverb to verbose
 *	Declared initialized matrix arrays (ov1, cov, suml, vv1x, vv1y, ctil1x,
 *	  ctil1y, ctil2x, ctil2y) dynamic with calloc and free
 *	Declared initialized matrix arrays (fdif, freq, err) dynamic with static
 *	  dimensions
 *	Declared all other arrays static, NULL initialized and realloc'd, so
 *	  their contents are preserved up to maximum previous length
 *	Eliminated trap on nmodes > 90000
 *	Removed unused setting of term in setting of
 *	Removed pointless loop over j=1,1
 *	Added some debugging diagnostics
 *  v 1.2 frozen 2012.07.31
 *  v 1.3
 *	Replaced references to FORTRAN code in ola_subs.f to locally defined
 *	  functions (or those defined in rdutil.c)
 *	Removed debugging diagnostics, some unused and redundant variables;
 *	  addition of some minor efficiency variables
 *	Removed fixed size assignments in fwidth in favor of allocated memory
 *	Removed writing of averaging kernels and inversion coefficients, just
 *	  pass back needed arguments
 *	Accept kernel modes up to n=nct, l=lct rather than hardcoded limits of
 *	  50, 2000
 *	All arrays with sizes dependent on the values of numtrgs and nrpts
 *	  (except those which must be zeroed on each call) are passed in as
 *	  arguments or declared global; use independent arrays for x and y
 *	  solutions
 *  v 1.3 frozen 2019.01.31
 */
