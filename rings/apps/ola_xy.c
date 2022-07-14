/*
 *  C version of subroutine ola_ in ola_xy_v12.f
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
							      /*  prototypes  */
extern void dsytrf_ (char *, int *, double *, int*, int *, double *, int *,
    int *);

extern void dsytrs_ (char *, int *, int *, double *, int *, int *, double *,
    int *, int *);

extern void solve_ (int *, int *, int *, int *, double *, double *,
    int *, int *, int *);

extern void fwidth_ (int *, int*, int *, double *, double *, double *, int *);

int solve (int la, int lb, int norder, int cols, double *a, double *b) {
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
fprintf (stderr, "  initialization call to dsytrf_  for order %d returned %f\n",
norder, work[0]);
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
fprintf (stderr, "LDA = %d (%d), N = %d (%d)\n", lda, la, n, norder);
    free (work);
    return info;
  }

  uplo = 'U';
  n = norder;
  nrhs = cols;
  lda = la;
  ldb = lb;
  dsytrs_ (&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  if (info) {
    fprintf (stderr, "Error in solve():\n");
    fprintf (stderr, "  call to dsytrs_ returned %d\n", info);
  }

  free (work);
  return info;
}

int ola_xy (double *dl_inp, int *n_inp, double *f_inp, double *ux_inp,
    double *euxinp, double *uy_inp, double *euyinp, int nmodes, char *filek1,
    double *x0, double *rcgx, double *rcgy, double *quartx, double *quarty, 
    double *solnx, double *solny, double *errx, double *erry, char *filker,
    char *filcoef, int qave, int qcoef, int verbose, double ob, double oe,
    int num, double beg, double endd, double amu) {
  FILE *lun19, *lun21, *lun27, *lun29, *lun37;
  double *ov1, *cov, *suml, *vv1x, *vv1y;
  static double *aa1x = NULL, *aa1y = NULL;
  static double *avccx = NULL, *avccy= NULL;
  static double *avcx = NULL, *avcy = NULL, *cencx = NULL, *cency = NULL;
  static double *ckintx = NULL, *ckinty = NULL;
  double *ctil1x = NULL, *ctil1y = NULL, *ctil2x = NULL, *ctil2y = NULL;
  static double *rad = NULL, *fakc = NULL, *rl = NULL;
  static double *om = NULL, *dom = NULL, *error = NULL, *weight = NULL;
/*
  double freq[2000][50][2], err[2000][50][2];
  double fdif[2000][50];
*/
  double freq[3000][50][2], err[3000][50][2];
  double fdif[3000][50];
  double dthcx[3], dthcy[3];
  double dif1, dif2, erri1, erri2, sume1, sume2, sumcx, sumcy;
  double cc, dx, ooo, ot;
  double sum1, sum2, sum3, sa, sb, sc, term, sumccx, sumccy;
  double sumc1x, sumc2x, sumc1y, sumc2y;
  double sum1x, er1x, s1x, sum1y, er1y, s1y;
  float *ak;
  static float *rd = NULL;
  float rs, rms;
  static int *n = NULL, *ipiv = NULL;
  int nmode, np, nsp, numker, nktmp, ntab, nmdtmp;
  int lll, nnn;
  int nrhs, ierr;
  int i, j, jj, k, nn, /* inum,*/ itnum;
  int status;

  int npt = 4000, nmd = 3100, nx0 = 100;

fprintf (stderr, "begin ola_xy()\n");
  aa1x = (double *)realloc (aa1x, nmd * nmd * sizeof (double));
  aa1y = (double *)realloc (aa1y, nmd * nmd * sizeof (double));
  avccx = (double *)realloc (avccx, nx0 * npt * sizeof (double));
  avccy = (double *)realloc (avccy, nx0 * npt * sizeof (double));
  avcx = (double *)realloc (avcx, npt * sizeof (double));
  avcy = (double *)realloc (avcy, npt * sizeof (double));
  cencx = (double *)realloc (cencx, nx0 * sizeof (double));
  cency = (double *)realloc (cency, nx0 * sizeof (double));
  ckintx = (double *)realloc (ckintx, nx0 * sizeof (double));
  ckinty = (double *)realloc (ckinty, nx0 * sizeof (double));
  rad = (double *)realloc (rad, npt * sizeof (double));
  fakc = (double *)realloc (fakc, npt * nmd * sizeof (double));
  rl = (double *)realloc (rl, nmd * sizeof (double));
  om = (double *)realloc (om, nmd * sizeof (double));
  dom = (double *)realloc (dom, 2 * nmd * sizeof (double));
  error = (double *)realloc (error, 2 * nmd * sizeof (double));
  weight = (double *)realloc (weight, npt * sizeof (double));
  n = (int *)realloc (n, nmd * sizeof (int));
  ipiv = (int *)realloc (ipiv, nmd * sizeof (int));

  if (num > nx0) return 1;
  if (verbose) {
    printf ("frequency limits: %.3f - %.3f mHz\n", ob, oe);
    printf ("error trade-off parameter = %.5f\n", amu);
  }
  lun21 = fopen (filek1, "r");
  if (!lun21) {
    fprintf (stderr, "Error: unable to open input file %s\n", filek1);
    return 1;
  }

  if (qave) {
				 /*  open output files for averaging kernels  */
    char *filkern = (char *)malloc (sizeof (char) * (strlen (filker) + 3));
    sprintf (filkern, "Ux_%s", filker);
    lun19 = fopen (filkern, "w");
    if (!lun19) {
      fprintf (stderr, "Error: unable to open kernel output file %s\n", filkern);
      fclose (lun21);
      return 1;
    }
    sprintf (filkern, "Uy_%s", filker);
    lun29 = fopen (filkern, "w");
    if (!lun29) {
      fprintf (stderr, "Error: unable to open kernel output file %s\n", filkern);
      fclose (lun21);
      fclose (lun19);
      return 1;
    }
    free (filkern);
  }
  if (qcoef) {
			    /*  open output files for inversion coefficients  */
    char *filkern = (char *)malloc (sizeof (char) * (strlen (filcoef) + 3));
    sprintf (filkern, "Ux_%s", filcoef);
    lun27 = fopen (filkern, "w");
    if (!lun27) {
      fprintf (stderr, "Error: unable to open coefficient output file %s\n",
	  filkern);
      fclose (lun21);
      if (qave) {
	fclose (lun19);
	fclose (lun29);
      }
      return 1;
    }
    sprintf (filkern, "Uy_%s", filcoef);
    lun37 = fopen (filkern, "w");
    if (!lun37) {
      fprintf (stderr, "Error: unable to open coefficient output file %s\n",
	  filkern);
      fclose (lun37);
      fclose (lun21);
      if (qave) {
	fclose (lun19);
	fclose (lun29);
      }
      return 1;
    }
    free (filkern);
  }
						     /*  initialize matrices  */
  bzero (fdif, sizeof (fdif));
  bzero (freq, sizeof (freq));
  bzero (err, sizeof (err));

  ov1 = (double *)calloc (nmd * nmd, sizeof (double));
  cov = (double *)calloc (2 * nmd * nmd, sizeof (double));
  suml = (double *)calloc (3 * nmd * nmd, sizeof (double));

  vv1x = (double *)calloc (num * nmd, sizeof (double));
  vv1y = (double *)calloc (num * nmd, sizeof (double));
  ctil1x = (double *)calloc (nmd * nx0, sizeof (double));
  ctil1y = (double *)calloc (nmd * nx0, sizeof (double));
  ctil2x = (double *)calloc (nmd * nx0, sizeof (double));
  ctil2y = (double *)calloc (nmd * nx0, sizeof (double));
						 /*  read in observed errors  */
  nsp = 0;
  for (i = 0; i < nmodes; i++) {
    int ni = n_inp[i];
    int l = dl_inp[i];
    double f = 0.001 * f_inp[i];
    if (f > ob && f < oe) {
      freq[l][ni][0] = ux_inp[i];
      freq[l][ni][1] = uy_inp[i];
      err[l][ni][0] = euxinp[i];
      err[l][ni][1] = euyinp[i];
      fdif[l][ni] = f_inp[i];
      nsp++;
    }
  }

  if (verbose) {
    printf ("nsp = %d\n", nsp);
    printf ("%d frequencies and velocities\n", i);
    printf ("MESSAGE: READ IN DATA\n");
  }
fprintf (stderr, "ola_xy(): reading kernel\n");
						     /*  read in the kernels  */
						 /*  read in radius and mass  */
  fscanf (lun21, "%g %g", &rs, &rms);
							    /*  read in mesh  */
  fscanf (lun21, "%d", &np);
fprintf (stderr, "ola_xy(): np = %d, npt = %d\n", np, npt);
  rd = (float *)realloc (rd, np * sizeof (float));
  for (i = 0; i < np; i++) fscanf (lun21, "%g", &(rd[i]));
  ak = (float *)malloc (np * sizeof (float));
  if (verbose) printf ("R = %.7e, M = %.3e\n", rs, rms);
  nmode = 0;
  for (i = 0; i < 20000; i++) {
    if (fscanf (lun21, "%d %d %lg", &lll, &nnn, &ooo) != 3) break;
    for (nn = 0; nn < np; nn++) fscanf (lun21, "%g", &(ak[nn]));
    if (lll >= 2000) continue;
    if (nnn >= 50) continue;
    dif1 = freq[lll][nnn][0];
    dif2 = freq[lll][nnn][1];
    erri1 = err[lll][nnn][0];
    erri2 = err[lll][nnn][1];
    if (erri1 <= 1.0e-11) continue;
    if (erri2 <= 1.0e-11) continue;
    ot = fdif[lll][nnn];
					      /* turning floats into doubles  */
/*
    for (jj = 0, inum = 0; jj < np - 1; jj++, inum++) {
      fakc[npt*nmode + inum] = ak[jj];
      rad[inum] = rd[jj];
*/
    for (jj = 0; jj < np - 1; jj++) {
      fakc[npt*nmode + jj] = ak[jj];
      rad[jj] = rd[jj];
    }
    n[nmode] = nnn;
    rl[nmode] = lll;
    om[nmode] = 0.001 * ot;
    error[2*nmode] = erri1;
    error[2*nmode + 1] = erri2;
    dom[2*nmode] = dif1;
    dom[2*nmode + 1] = dif2;
    nmode++;
  }
  free (ak);

//  np = inum;
  if (verbose) printf ("np = %d\n", np);
  sume1 = sume2 = 0.0;
  for (i = 0; i < nmode; i++) {
    sume1 += 1.0 / (error[2*i] * error[2*i]);
    sume2 += 1.0 / (error[2*i + 1] * error[2*i + 1]);
  }
  if (verbose) printf ("read kernels %d\n", nmode);
						   /*  set covariance matrix  */
fprintf (stderr, "ola_xy(): setting covariance matrix\n");
  sumcx = sumcy = 0.0;
  for (j = 0; j < nmode; j++) {
    sumcx += error[2*j] * error[2*j];
    cov[2*nmd*j + 2*j + 0] = error[2*j] * error[2*j];
    sumcy += error[2*j + 1] * error[2*j + 1];
    cov[2*nmd*j + 2*j + 1] = error[2*j + 1] * error[2*j + 1];
  }
  sumcx /= nmode;
  sumcy /= nmode;
  if (verbose) {
    printf ("sumc: %g, %g\n", sumcx, sumcy);
    printf ("calc. integration weights\n");
  }
						 /*  get integration weights  */
fprintf (stderr, "ola_xy(): calculating integration weights\n");
fprintf (stderr, "ola_xy(): np = %d, npt = %d\n", np, npt);
  weight[0] = 0.5 * (rad[0] - rad[1]);
  for (j = 0; j < np - 2; j++) weight[j+1] = 0.5 * (rad[j] - rad[j+2]);
  weight[np-1] = 0.5 * (rad[np-2] - rad[np-1]);
					       /*  Now the OLA related stuff  */
fprintf (stderr, "ola_xy(): calculating overlap matrix\n");
  if (verbose) printf ("calc overlap matrix\n");
						  /*  set the overlap matrix  */
  for (i = 0; i < nmode; i++) {
    for (j = 0; j < nmode; j++) {
      sum2 = sum3 = sa = sb = sc = 0.0;
      for (k = 0; k < np; k++) {
	term = weight[k] * fakc[npt*i + k] * fakc[npt*j + k];
	sa += term * rad[k] * rad[k];
	sb += term * rad[k];
	sc += term;
      }
      suml[nmd*nmd*0 + nmd*j + i] = suml[nmd*nmd*0 + nmd*i + j] = sa;
      suml[nmd*nmd*1 + nmd*j + i] = suml[nmd*nmd*1 + nmd*i + j] = sb;
      suml[nmd*nmd*2 + nmd*j + i] = suml[nmd*nmd*2 + nmd*i + j] = sc;
    }
  }
						 /*  the Lagrange multiplier  */
  for (i = 0; i < nmode; i++) {
    sum1 = 0.0;
    for (k = 0; k < np; k++) sum1 += weight[k] * fakc[k + npt*i];
    ov1[i + nmode*nmd] = ov1[nmode + i*nmd] = 0.5 * sum1;
    suml[i + nmode*nmd] = suml[nmode + i*nmd] = 0.5 * sum1;
  }

  numker = nmode + 1;
  if (verbose) printf ("numker = %d\n", numker);
							/*  set target radii  */
  if (verbose) printf ("calc. targets\n");
  ntab = np;
  dx = (endd - beg) / (num - 1.0);
  for (i = 0; i < num; i++) x0[i] = beg + i * dx;
					     /*  need loop over X and Y here  */
  for (itnum = 0; itnum < num; itnum++) {
						  /*  set the overlap matrix  */
    if (verbose) printf ("calc. ov1\n");
    for (i = 0; i < nmode; i++) {
      for (j = i; j < nmode; j++) {
	sum1 = suml[nmd*nmd*0 + nmd*j + i] -
	    2 * x0[itnum] * suml[nmd*nmd*1 + nmd*j + i] +
	    x0[itnum] * x0[itnum] * suml[nmd*nmd*2 + nmd*j + i];
	ov1[i + j*nmd] = ov1[j + i*nmd] = sum1;
      }
    }
					   /*  set inhomogeneous part of rhs  */
    for (j = 0; j < num; j++)
      vv1x[nmode + j*nmd] = vv1y[nmode + j*nmd] = 0.5;
						      /*  set matrix for lhs  */
    if (verbose)  printf ("calc. lhs\n");
    for (j = 0; j < nmode; j++) {
      for (i = 0; i < nmode; i++) {
	cc = cov[j*2*nmd + i*2 + 0] / sumcx;
	aa1x[i + j*nmd] = ov1[i + j*nmd] + amu * cc;
 	cc = cov[j*2*nmd + i*2 + 1] / sumcy;
	aa1y[i + j*nmd] = ov1[i + j*nmd] + amu * cc;
     }
    }
    if (verbose) printf ("numker = %d %d\n", numker, nmode);
    for (j = 0; j < numker; j++) {
      for (i = nmode; i < numker; i++) {
	aa1x[i + j*nmd] = aa1y[i + j*nmd] = ov1[i + j*nmd];
	aa1x[j + i*nmd] = aa1y[j + i*nmd] = ov1[j + i*nmd];
      }
      for (k = 0; k < num; k++) {
	ctil1x[j + k*nmd] = vv1x[j + k*nmd];
	ctil1y[j + k*nmd] = vv1y[j + k*nmd];
      }
    }
						  /*  solving the equations  */
/*
    nrhs = 1;
    nktmp = numker;
    nmdtmp = nmd;
    solve_ (&nmdtmp, &nmdtmp, &nktmp, &nrhs, aa1x, ctil1x, &ierr, ipiv, &verbose);
*/
//fprintf (stderr, "ola_xy(): calling solve() for ux\n");
    status = solve (nmd, nmd, numker, 1, aa1x, ctil1x);

    if (verbose) printf ("done Ux\n");
/*
    nrhs = 1;
    nktmp = numker;
    nmdtmp = nmd;
    solve_ (&nmdtmp, &nmdtmp, &nktmp, &nrhs, aa1y, ctil1y, &ierr, ipiv, &verbose);
*/
//fprintf (stderr, "ola_xy(): calling solve() for uy\n");
    status = solve (nmd, nmd, numker, 1, aa1y, ctil1y);

    if (verbose) printf ("done Uy\n");
    jj = 0;
    for (i = 0; i < numker; i++) {
      ctil2x[i + itnum*nmd] = ctil1x[i + jj*nmd];
      ctil2y[i + itnum*nmd] = ctil1y[i + jj*nmd];
    }
					   /*  Calculating averaging kernel  */
    for (j = 0; j < np; j++) {
      sumccx = sumccy = 0.0;
      for (i = 0; i < nmode; i++) {
	sumccx += ctil1x[i + jj*nmd] * fakc[j + i*npt];
	sumccy += ctil1y[i + jj*nmd] * fakc[j + i*npt];
      }
      avcx[j] = sumccx;
      avccx[nx0*j + itnum] = sumccx;
      avcy[j] = sumccy;
      avccy[nx0*j + itnum] = sumccy;
    }
					 /*  Finding cg of averaging kernel  */
    sumc1x = sumc2x = sumc1y = sumc2y = 0.0;
     for (j = 0; j < np; j++) {
      sumc1x += weight[j] * rad[j] * avcx[j];
      sumc2x += weight[j] * avcx[j];
      sumc1y += weight[j] * rad[j] * avcy[j];
      sumc2y += weight[j] * avcy[j];
    }
    if (verbose) printf ("int. %g\n", sumc1x);
    cencx[0] = sumc1x / sumc2x;
    ckintx[0] = sumc2x;
    cency[0] = sumc1y / sumc2y;
    ckinty[0] = sumc2y;
    if (verbose) printf ("calling width for Ux\n");
    fwidth_ (&npt, &np, &num, rad, avcx, dthcx, &verbose);
    if (verbose) printf ("width done for Ux: %g %g %g\n",
	dthcx[0], dthcx[1], dthcx[2]);
    if (verbose) printf ("calling width for Uy\n");
    fwidth_ (&npt, &np, &num, rad, avcy, dthcy, &verbose);
    if (verbose) printf ("width done for Uy: %g %g %g\n",
	dthcy[0], dthcy[1], dthcy[2]);
							    /*  the solution  */
    sum1x = er1x = sum1y = er1y = 0.0;
//fprintf (stderr, "ola_xy(): sums zeroed\n");
    for (i = 0; i < nmode; i++) {
//fprintf (stderr, "ola_xy(): incrementing sums  for %d\n", i);
      sum1x += ctil1x[i] * dom[2*i];
      er1x += ctil1x[i] * ctil1x[i] * cov[2*nmd*i + 2*i];
      sum1y += ctil1y[i] * dom[2*i + 1];
      er1y += ctil1y[i] * ctil1y[i] * cov[2*nmd*i + 2*i + 1];
//fprintf (stderr, "ola_xy(): sums incremented for %d\n", i);
    }
/*
    sum1y = er1y = 0.0;
    for (i = 0; i < nmode; i++) {
      sum1y += ctil1y[i] * dom[2*i + 1];
      er1y += ctil1y[i] * ctil1y[i] * cov[2*nmd*i + 2*i + 1];
    }
*/
//fprintf (stderr, "ola_xy(): sums calculated\n");
    s1x = fabs (dthcx[0] - dthcx[2]);
    s1y = fabs (dthcy[0] - dthcy[2]);
    er1x = sqrt (er1x);
    er1y = sqrt (er1y);
    if (verbose) printf ("%g %g %g\n", x0[itnum], sum1x, er1x);
    rcgx[itnum] = cencx[0];
    rcgy[itnum] = cency[0];
    for (i = 0; i < 3; i++) {
      quartx[i + 3*itnum] = dthcx[i];
      quarty[i + 3*itnum] = dthcy[i];
    }
    solnx[itnum] = sum1x;
    solny[itnum] = sum1y;
    errx[itnum] = er1x;
    erry[itnum] = er1y;
/*
    fprintf (lun10, "%7.5f %13.5e %6.4f %6.4f %6.4f %6.4f %6.4f %13.5e %13.5e\n",
	x0[itnum], amu, cencx[0], dthcx[2], dthcx[1], dthcx[0], s1x, sum1x,
	er1x);
    fprintf (lun20, "%7.5f %13.5e %6.4f %6.4f %6.4f %6.4f %6.4f %13.5e %13.5e\n",
	x0[itnum], amu, cency[0], dthcy[2], dthcy[1], dthcy[0], s1y, sum1y,
	er1y);
*/
  }

  if (qave) {
				         /*  write out the averaging kernels  */
    fprintf (lun19, "# radius ");
    fprintf (lun29, "# radius ");
    for (i = 0; i < num; i++) fprintf (lun19, "%13.5e", x0[i]);
    for (i = 0; i < num; i++) fprintf (lun29, "%13.5e", x0[i]);
    for (j = 0; j < np; j += 2) {
      fprintf (lun19, "%14.6e", rad[j]);
      for (k = 0; k < num; k++) fprintf (lun19, "%14.6e", avccx[nx0*j + k]);
      fprintf (lun29, "%14.6e", rad[j]);
      for (k = 0; k < num; k++) fprintf (lun29, "%14.6e", avccy[nx0*j + k]);
    }
    fclose (lun19);
    fclose (lun29);
  }
  if (qcoef) {
				        /*  write out inversion coefficients  */
    fprintf (lun27, "# ");
    for (i = 0; i < num; i++) fprintf (lun27, "%13.5e", x0[i]);
    fprintf (lun27, "\n");
    fprintf (lun37, "# ");
    for (i = 0; i < num; i++) fprintf (lun37, "%13.5e", x0[i]);
    fprintf (lun37, "\n");
    for (j = 0; j < numker - 1; j++) {
      for (i = 0; i < num; i++) fprintf (lun27, "%14.6e", ctil2x[i + j * nmd]);
      fprintf (lun27, "\n");
      for (i = 0; i < num; i++) fprintf (lun37, "%14.6e", ctil2y[i + j * nmd]);
      fprintf (lun37, "\n");
    }
    fclose (lun27);
    fclose (lun37);
  }
  fclose (lun21);

  free (ov1);
  free (cov);
  free (suml);
  free (vv1x);
  free (vv1y);
  free (ctil1x);
  free (ctil1y);
  free (ctil2x);
  free (ctil2y);
  
  return 0;
}
