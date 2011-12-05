/*
 *  C version of subroutine ola_ in ola_xf.f included in rdvinv_v04
 *
 *  N.B.: this code has a bug, and seg faults at line 230
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
							     /*  prototypes  */
extern void solve_ (int *, int *, int *, int *, double *, double *,
    int *, int *, int *);

extern void fwidth_ (int *, int*, int *, double *, double *, double *, int *);

void ola_xy (double *dl_inp, int *n_inp, double *f_inp, double *ux_inp,
    double *euxinp, double *uy_inp, double *euyinp, int nmodes,
    char *filek1, char *filsolx, char *filsoly, char *filker, char *filcoef,
    int qave, int qcoef, int verbose, double ob, double oe, int num,
    double beg, double endd, double amu) {
  FILE *lun10, *lun19, *lun20, *lun21, *lun27, *lun29, *lun37;
  double *ov1, *cov, *suml, *vv1x, *vv1y;
  double freq[2000][50][2], err[2000][50][2];
  double fdif[2000][50];
  double dthcx[3], dthcy[3];
  double dif1, dif2, erri1, erri2, sume1, sume2, sumcx, sumcy;
  double cc, dx, ooo, ot;
  double sum1, sum2, sum3, sa, sb, sc, term, sumccx, sumccy;
  double sumc1x, sumc2x, sumc1y, sumc2y;
  double sum1x, er1x, s1x, sum1y, er1y, s1y;
  float *ak;
  float rs, rms;
  int nmode, np, nsp, numker, ntab;
  int lll, nnn;
  int nrhs, ierr, ipiv;
  int i, j, jj, k, nn, inum, itnum;

  int npt = 4000, nmd = 3100, nx0 = 100;

  double *aa1x = (double *)malloc (nmd * nmd * sizeof (double));
  double *aa1y = (double *)malloc (nmd * nmd * sizeof (double));
  double *avccx = (double *)malloc (nx0 * npt * sizeof (double));
  double *avccy = (double *)malloc (nx0 * npt * sizeof (double));
  double *avcx = (double *)malloc (npt * sizeof (double));
  double *avcy = (double *)malloc (npt * sizeof (double));
  double *cencx = (double *)malloc (nx0 * sizeof (double));
  double *cency = (double *)malloc (nx0 * sizeof (double));
  double *ckintx = (double *)malloc (nx0 * sizeof (double));
  double *ckinty = (double *)malloc (nx0 * sizeof (double));
  double *ctil1x = (double *)malloc (nmd * nx0 * sizeof (double));
  double *ctil1y = (double *)malloc (nmd * nx0 * sizeof (double));
  double *ctil2x = (double *)malloc (nmd * nx0 * sizeof (double));
  double *ctil2y = (double *)malloc (nmd * nx0 * sizeof (double));
  double *x0 = (double *)malloc (nx0 * sizeof (double));
  double *rad = (double *)malloc (npt * sizeof (double));
  double *fakc = (double *)malloc (npt * nmd * sizeof (double));
  double *rl = (double *)malloc (nmd * sizeof (double));
  double *om = (double *)malloc (nmd * sizeof (double));
  double *dom = (double *)malloc (2 * nmd * sizeof (double));
  double *error = (double *)malloc (2 * nmd * sizeof (double));
  double *weight = (double *)malloc (npt * sizeof (double));
  float *rd = (float *)malloc (npt * sizeof (float));
  int *n = (int *)malloc (nmd * sizeof (int));

  if (verbose) {
    printf ("frequency limits: %.3f - %.3f mHz\n", ob, oe);
    printf ("error trade-off parameter = %.5f\n", amu);
  }
  lun21 = fopen (filek1, "r");
  if (!lun21) {
    fprintf (stderr, "Error: unable to open input file %s\n", filek1);
    return;
  }
  lun10 = fopen (filsolx, "w");
  if (!lun10) {
    fprintf (stderr, "Error: unable to open output file %s\n", filsolx);
    return;
  }
  lun20 = fopen (filsoly, "w");
  if (!lun20) {
    fprintf (stderr, "Error: unable to open output file %s\n", filsoly);
    return;
  }
  if (qave) {
				 /*  open output files for averaging kernels  */
    char *filkern = (char *)malloc (sizeof (char) * (strlen (filker) + 3));
    sprintf (filkern, "Ux_%s", filker);
    lun19 = fopen (filkern, "w");
    if (!lun19) {
      fprintf (stderr, "Error: unable to open kernel output file %s\n", filkern);
      return;
    }
    sprintf (filkern, "Uy_%s", filker);
    lun29 = fopen (filkern, "w");
    if (!lun29) {
      fprintf (stderr, "Error: unable to open kernel output file %s\n", filkern);
      return;
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
      return;
    }
    sprintf (filkern, "Uy_%s", filcoef);
    lun37 = fopen (filkern, "w");
    if (!lun37) {
      fprintf (stderr, "Error: unable to open coefficient output file %s\n",
	  filkern);
      return;
    }
    free (filkern);
  }
						     /*  initialize matrices  */
  memset (fdif, 0, sizeof (fdif));
  memset (freq, 0, sizeof (freq));
  memset (err, 0, sizeof (err));

  ov1 = (double *)calloc (nmd * nmd, sizeof (double));
  cov = (double *)calloc (2 * nmd * nmd, sizeof (double));
  suml = (double *)calloc (3 * nmd * nmd, sizeof (double));

  vv1x = (double *)calloc (num * nmd, sizeof (double));
  vv1y = (double *)calloc (num * nmd, sizeof (double));
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
						     /*  read in the kernels  */
						 /*  read in radius and mass  */
  fscanf (lun21, "%g %g", &rs, &rms);
							    /*  read in mesh  */
  fscanf (lun21, "%d", &np);
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
    if (erri1 <= 1.0e-12) continue;
    if (erri2 <= 1.0e-12) continue;
    ot = fdif[lll][nnn];
					      /* turning floats into doubles  */
    for (jj = 0, inum = 0; jj < np - 1; jj++, inum++) {
      fakc[npt*nmode + inum] = ak[jj];
      rad[inum] = rd[jj];
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

  np = inum;
  if (verbose) printf ("np = %d\n", np);
  sume1 = sume2 = 0.0;
  for (i = 0; i < nmode; i++) {
    sume1 += 1.0 / (error[2*i] * error[2*i]);
    sume2 += 1.0 / (error[2*i + 1] * error[2*i + 1]);
  }
  if (verbose) printf ("read kernels %d\n", nmode);
						   /*  set covariance matrix  */
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
  weight[0] = 0.5 * (rad[0] - rad[1]);
  for (j = 0; j < np - 2; j++) weight[j+1] = 0.5 * (rad[j] - rad[j+2]);
  weight[np-1] = 0.5 * (rad[np-2] - rad[np-1]);
					       /*  Now the OLA related stuff  */
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
    for (k = 0; k < np; k++) sum1 += weight[k] * fakc[npt*i + k];
    ov1[nmd*nmode + i] = ov1[nmd*i + nmode] = 0.5 * sum1;
    suml[nmd*nmode + i] = suml[nmd*i + nmode] = 0.5 * sum1;
  }

  numker = nmode + 1;
  if (verbose) printf ("numker = %d\n", numker);
							/*  set target radii  */
  if (verbose) printf ("calc. targets\n");
  ntab = np;
  dx = (endd - beg) / (num - 1.0);
  for (i = 0; i < num; i++) x0[i] = beg + i * dx;
fprintf (stderr, "writing to lun10\n");
  fprintf (lun10, "#target rad, err-suppr. param, CG of av. kernel, \n");
  fprintf (lun10, " 1st,2nd,3rd quart. pt, (3rd-1st)quart.pt., soln, err\n");
fprintf (stderr, "writing to lun20\n");
  fprintf (lun20, "#target rad, err-suppr. param, CG of av. kernel, \n");
  fprintf (lun20, " 1st,2nd,3rd quart. pt, (3rd-1st)quart.pt., soln, err\n");						   /* Now for the Solutions  */
					     /*  need loop over X and Y here  */
  for (itnum = 0; itnum < num; itnum++) {
						  /*  set the overlap matrix  */
    if (verbose) printf ("calc. ov1\n");
    for (i = 0; i < nmode; i++) {
      for (j = i; j < nmode; j++) {
	sum1 = suml[nmd*nmd*0 + nmd*j + i] -
	    2 * x0[itnum] * suml[nmd*nmd*1 + nmd*j + i] +
	    x0[itnum] * x0[itnum] * suml[nmd*nmd*2 + nmd*j + i];
	ov1[nmd*j + i] = ov1[nmd*i + j] = sum1;
      }
    }
fprintf (stderr, "overlap matrix set\n");
fprintf (stderr, "allocation on vv1x/y: %d\n", num * nmd);
fprintf (stderr, "max indices on vv1x/y: %d\n", num*nmd + nmode);
					   /*  set inhomogeneous part of rhs  */
    for (j = 0; j < num; j++)
      vv1x[j*nmd + nmode] = vv1y[j*nmd + nmode] = 0.5;
						      /*  set matrix for lhs  */
    if (verbose)  printf ("calc. lhs\n");
    for (j = 0; j < nmode; j++) {
      for (i = 0; i < nmode; i++) {
	cc = cov[j*2*nmd + i*2 + 0] / sumcx;
	aa1x[nmd*j + i] = ov1[nmd*j + i] + amu * cc;
 	cc = cov[j*2*nmd + i*2 + 1] / sumcy;
	aa1y[nmd*j + i] = ov1[nmd*j + i] + amu * cc;
     }
    }
    if (verbose) printf ("numker = %d %d\n", numker, nmode);
    for (j = 0; j < numker; j++) {
      for (i = nmode; i < numker; i++) {
	aa1x[nmd*j + i] = aa1y[nmd*j + i] = ov1[j*nmd + i];
	aa1x[nmd*i + j] = aa1y[nmd*i + j] = ov1[i*nmd + j];
      }
      for (k = 0; k < num; k++) {
	ctil1x[nmd*k + j] = vv1x[nmd*k + j];
	ctil1y[nmd*k + j] = vv1y[nmd*k + j];
      }
    }
						  /*  solving the equations  */
fprintf (stderr, "Before FTN solve nmd = %d, numker = %d\n", nmd, numker);
    if (verbose) printf ("solving eqn\n");
    nrhs = 1;
    solve_ (&nmd, &nmd, &numker, &nrhs, aa1x, ctil1x, &ierr, &ipiv, &verbose);
    if (verbose) printf ("done Ux\n");
    nrhs = 1;
    solve_ (&nmd, &nmd, &numker, &nrhs, aa1y, ctil1y, &ierr, &ipiv, &verbose);
    if (verbose) printf ("done Uy\n");
    for (i = 0; i < numker; i++) {
      ctil2x[itnum*nmd + i] = ctil1x[i];
      ctil2y[itnum*nmd + i] = ctil1y[i];
    }
fprintf (stderr, "After FTN solve nmd = %d, numker = %d\n", nmd, numker);
					   /*  Calculating averaging kernel  */
fprintf (stderr, "Calculating averaging kernel: np = %d\n", np);
    for (j = 0; j < np; j++) {
fprintf (stderr, "j = %d, nmode = %d\n", j, nmode);
      sumccx = sumccy = 0.0;
fprintf (stderr, "allocations of ctil1x/y, fakc = %d, %d\n", nmd*nx0, npt*nmd);
fprintf (stderr, "max indices on ctil1x/y, fakc = %d, %d\n", nmode, np + npt*nmode);
      for (i = 0; i < nmode; i++) {
if (i >= 17515 || i <= 1) fprintf (stderr, "i = %d, nmode = %d; (i < nmode) ? %d\n",
i, nmode, (i < nmode));
	sumccx += ctil1x[i] * fakc[npt*i + j];
	sumccy += ctil1y[i] * fakc[npt*i + j];
      }
fprintf (stderr, "finished i loop\n");
      avcx[j] = sumccx;
      avccx[nx0*j + itnum] = sumccx;
      avcy[j] = sumccy;
      avccy[nx0*j + itnum] = sumccy;
    }
					 /*  Finding cg of averaging kernel  */
fprintf (stderr, "Finding CG of averaging kernel\n");
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
    for (i = 0; i < nmode; i++) {
      sum1x += ctil1x[i] * dom[2*i];
      er1x += ctil1x[i] * ctil1x[i] * cov[2*nmd*i + 2*i];
      sum1y+= ctil1y[i] * dom[2*i + 1];
      er1y += ctil1y[i] * ctil1y[i] * cov[2*nmd*i + 2*i + 1];
    }
    s1x = fabs (dthcx[0] - dthcx[2]);
    s1y = fabs (dthcy[0] - dthcy[2]);
    er1x = sqrt (er1x);
    er1y = sqrt (er1y);
    if (verbose) printf ("%g %g %g\n", x0[itnum], sum1x, er1x);
    fprintf (lun10, "%7.5f %13.5e %6.4f %6.4f %6.4f %6.4f %6.4f %13.5e %13.5e\n",
	x0[itnum], amu, cencx[0], dthcx[2], dthcx[1], dthcx[0], s1x, sum1x,
	er1x);
    fprintf (lun20, "%7.5f %13.5e %6.4f %6.4f %6.4f %6.4f %6.4f %13.5e %13.5e\n",
	x0[itnum], amu, cency[0], dthcy[2], dthcy[1], dthcy[0], s1y, sum1y,
	er1y);
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
  fclose (lun10);
  fclose (lun20);
  fclose (lun21);
  return;
}
