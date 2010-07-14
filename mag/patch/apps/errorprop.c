#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
static void ccker (double *u, double s) {
  double s2, s3;

  s2= s * s;
  s3= s2 * s;
  u[0] = s2 - 0.5 * (s3 + s);
  u[1] = 1.5*s3 - 2.5*s2 + 1.0;
  u[2] = -1.5*s3 + 2.0*s2 + 0.5*s;
  u[3] = 0.5 * (s3 - s2);
}
*/

int errorprop (float *bTotal, float *bAzim, float *bInc, float *bFill, float *ebT, float *ebA, float *ebI, float *ebF,
               float *ebTbA, float *ebTbI, float *ebTbF, float *ebIbA, float *ebFbA, float *ebFbI, 
               double lon, double lat, double lonc, double latc, double pAng, int nx, int ny, double x, double y,
               double *BtVariance, double *BpVariance, double *BrVariance) {

/* Inputs:
*  bTotal, bFill, bInc, bAzim: magnetic field inverted by an inversion code. They are
*                              field strength, inclination, and azimuth.
* ebT, ebF, ebI, ebA:          variance of bTotal, bInc, bAzim.
* ebTbF, ebFbI, ebFbA:         covariance of bTotal, bFill, bInc, bAzim.
* ebTbI, ebTbA, ebIbA:         covariance of bTotal, bFill, bInc, bAzim.
* lon, lat:                    heliographic coordinates of the location of map that the data
*                              is projected.
* lonc, latc:                  heliographic coordinates of the image disk center.
* pAng:                        position angle of the heliographic north pole, measured eastward
*                              from the north.
* nx, ny:                      size of the array.
* x, y:                        location of the pixel in image coordinates.
*
* Outputs:
* BrVariance, BtVariance, BpVariance: variances of these three components of magnetic
*                                     field.
*
* Note:                 Interpolation useed here is a cubic convolution interpolation.
*/

  static double raddeg = M_PI / 180.;
  double ux[4], uy[4];
  double *b, *fill, *inc, *azim, *weigh;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double *dBrdBtotal, *dBrdFill, *dBrdInc, *dBrdAzim;
  double *dBtdBtotal, *dBtdFill, *dBtdInc, *dBtdAzim;
  double *dBpdBtotal, *dBpdFill, *dBpdInc, *dBpdAzim;
  double *errBT, *errFILL, *errINC, *errAZ;
  double *errBTFILL, *errBTINC, *errBTAZ, *errFILLINC, *errFILLAZ, *errINCAZ;
  double BrSigma2, BtSigma2, BpSigma2;
  int xbox = 4, ybox = 4, pixSize = xbox * ybox;
  int ix, iy, ix1, iy1, i, j;

  if (x < 1. || x >= (float)(nx-2) || y < 1. || y >= (float)(ny-2))
  return (1);

  b = (double *)malloc (pixSize * sizeof (double));
  fill = (double *)malloc (pixSize * sizeof (double));
  inc = (double *)malloc (pixSize * sizeof (double));
  azim = (double *)malloc (pixSize * sizeof (double));
  weigh = (double *)malloc (pixSize * sizeof (double));

  dBrdBtotal = (double *)malloc (pixSize * sizeof (double));
  dBrdFill = (double *)malloc (pixSize * sizeof (double));
  dBrdInc = (double *)malloc (pixSize * sizeof (double));
  dBrdAzim = (double *)malloc (pixSize * sizeof (double));
  dBtdBtotal = (double *)malloc (pixSize * sizeof (double));
  dBtdFill = (double *)malloc (pixSize * sizeof (double));
  dBtdInc = (double *)malloc (pixSize * sizeof (double));
  dBtdAzim = (double *)malloc (pixSize * sizeof (double));
  dBpdBtotal = (double *)malloc (pixSize * sizeof (double));
  dBpdFill = (double *)malloc (pixSize * sizeof (double));
  dBpdInc = (double *)malloc (pixSize * sizeof (double));
  dBpdAzim = (double *)malloc (pixSize * sizeof (double));

  errBT = (double *)malloc (pixSize * sizeof (double));
  errFILL = (double *)malloc (pixSize * sizeof (double));
  errINC = (double *)malloc (pixSize * sizeof (double));
  errAZ = (double *)malloc (pixSize * sizeof (double));
  errBTFILL = (double *)malloc (pixSize * sizeof (double));
  errBTINC = (double *)malloc (pixSize * sizeof (double));
  errBTAZ = (double *)malloc (pixSize * sizeof (double));
  errFILLINC = (double *)malloc (pixSize * sizeof (double));
  errFILLAZ = (double *)malloc (pixSize * sizeof (double));
  errINCAZ = (double *)malloc (pixSize * sizeof (double));

  a11 = -sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - lonc);
  a12 =  sin(latc) * cos(pAng) * sin(lon - lonc) + sin(pAng) * cos(lon - lonc);
  a13 = -cos(latc) * sin(lon - lonc);
  a21 = -sin(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc))
        - cos(lat) * cos(latc) * sin(pAng);
  a22 =  sin(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc))
        + cos(lat) * cos(latc) * cos(pAng);
  a23 = -cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat);
  a31 =  cos(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc))
        - sin(lat) * cos(latc) * sin(pAng);
  a32 = -cos(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc))
        + sin(lat) * cos(latc) * cos(pAng);
  a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc);

  ix = (int)x;
  ccker (ux,  x - (double)ix);
  iy = (int)y;
  ccker (uy,  y - (double)iy);

  ix1 = ix - 1;
  iy1 = iy - 1;
  int yOff = 0, iData = 0;
  int yAll = 0, iAll = 0;
  for (j = 0; j < ybox; j++) {
    int i = 0;
    yOff = j * xbox;
    yAll = (iy1 + j) * nx;

    for (i = 0; i < xbox; i++) {
      iData = yOff + i;
      iAll = yAll + (ix1 + i);

      weigh[iData] = uy[j] * ux[i]; /* weights from the cubic convolution interpolation */
      if (isnan(bTotal[iAll])) return(1);
      b[iData] = bTotal[iAll]; /* field strength from the pixels used for the interpolation */
      if (isnan(bFill[iAll])) return(1);
      fill[iData] = bFill[iAll]; /* fill factor */
      if (isnan(bInc[iAll])) return(1);
      inc[iData] = raddeg * bInc[iAll]; /* inclination */
      if (isnan(bAzim[iAll])) return(1);
      azim[iData] = raddeg * bAzim[iAll]; /* azimuth */
      if (isnan(ebT[iAll])) return(1);
      errBT[iData] = ebT[iAll]; /* variance of field strength */
      if (isnan(ebF[iAll])) return(1);
      errFILL[iData] = ebF[iAll]; /* variance of fill factor */
      if (isnan(ebI[iAll])) return(1);
      errINC[iData] = ebI[iAll]; /* variance of inclination */
      if (isnan(ebA[iAll])) return(1);
      errAZ[iData] = ebA[iAll]; /* variance of azimuth */
      if (isnan(ebTbF[iAll])) return(1);
      errBTFILL[iData] = ebTbF[iAll]; /* covariance of field and fill factor */
      if (isnan(ebTbI[iAll])) return(1);
      errBTINC[iData] = ebTbI[iAll]; /* covariance of field and inclination */
      if (isnan(ebTbA[iAll])) return(1);
      errBTAZ[iData] = ebTbA[iAll];  /* covariance of field and azimuth */
      if (isnan(ebFbI[iAll])) return(1);
      errFILLINC[iData] = ebFbI[iAll];  /* covariance of fill factor and inclination */
      if (isnan(ebFbA[iAll])) return(1);
      errFILLAZ[iData] = ebFbA[iAll];  /* covariance of fill factor and azimuth */
      if (isnan(ebIbA[iAll])) return(1);
      errINCAZ[iData] = ebIbA[iAll]; /* covariance of inclination and azimuth */
    }
  }

  for (j = 0; j < pixSize; j++) {

      dBtdBtotal[j] = weigh[j] * fill[j] * (a11 * sin(inc[j]) * cos(azim[j]) + 
                                            a12 * sin(inc[j]) * sin(azim[j]) + a13 * cos(inc[j]));
      dBtdFill[j] = weigh[j] * b[j] * (a11 * sin(inc[j]) * cos(azim[j]) + 
                                       a12 * sin(inc[j]) * sin(azim[j]) + a13 * cos(inc[j]));
      dBtdInc[j] = weigh[j] * b[j] * fill[j] * (a11 * cos(inc[j]) * cos(azim[j]) + 
                                                a12 * cos(inc[j]) * sin(azim[j]) - a13 * sin(inc[j]));
      dBtdAzim[j] = weigh[j] * b[j] * fill[j] * (-a11 * sin(inc[j]) * sin(azim[j]) + 
                                                 a12 * sin(inc[j]) * cos(azim[j]));
      dBpdBtotal[j] = weigh[j] * fill[j] * (a21 * sin(inc[j]) * cos(azim[j]) + 
                                            a22 * sin(inc[j]) * sin(azim[j]) + a23 * cos(inc[j]));
      dBpdFill[j] = weigh[j] * b[j] * (a21 * sin(inc[j]) * cos(azim[j]) + 
                                       a22 * sin(inc[j]) * sin(azim[j]) + a23 * cos(inc[j]));
      dBpdInc[j] = weigh[j] * b[j] * fill[j] * (a21 * cos(inc[j]) * cos(azim[j]) + 
                                                a22 * cos(inc[j]) * sin(azim[j]) - a23 * sin(inc[j]));
      dBpdAzim[j] = weigh[j] * b[j] * fill[j] * (-a21 * sin(inc[j]) * sin(azim[j]) + 
                                                 a22 * sin(inc[j]) * cos(azim[j]));
      dBrdBtotal[j] = weigh[j] * fill[j] * (a31 * sin(inc[j]) * cos(azim[j]) + 
                                            a32 * sin(inc[j]) * sin(azim[j]) + a33 * cos(inc[j]));
      dBrdFill[j] = weigh[j] * b[j] * (a31 * sin(inc[j]) * cos(azim[j]) + 
                                       a32 * sin(inc[j]) * sin(azim[j]) + a33 * cos(inc[j]));
      dBrdInc[j] = weigh[j] * b[j] * fill[j] * (a31 * cos(inc[j]) * cos(azim[j]) + 
                                                a32 * cos(inc[j]) * sin(azim[j]) - a33 * sin(inc[j]));
      dBrdAzim[j] = weigh[j] * b[j] * fill[j] * (-a31 * sin(inc[j]) * sin(azim[j]) + 
                                                 a32 * sin(inc[j]) * cos(azim[j]));
  }

  BrSigma2 = 0.0;
  BtSigma2 = 0.0;
  BpSigma2 = 0.0;

  for (j = 0; j < pixSize; j++) {

      BtSigma2 = BtSigma2 + dBtdBtotal[j] * dBtdBtotal[j] * errBT[j] + dBtdFill[j] * dBtdFill[j] * errFILL[j] + 
                            dBtdInc[j] * dBtdInc[j] * errINC[j] + dBtdAzim[j] * dBtdAzim[j] * errAZ[j] + 
                            2.0 * dBtdBtotal[j] * dBtdFill[j] * errBTFILL[j] + 2.0 * dBtdBtotal[j] * dBtdInc[j] * errBTINC[j] +
                            2.0 * dBtdBtotal[j] * dBtdAzim[j] * errBTAZ[j] + 2.0 * dBtdFill[j] * dBtdInc[j] * errFILLINC[j]  + 
                            2.0 * dBtdFill[j] * dBtdAzim[j] * errFILLAZ[j] + 2.0 * dBtdInc[j] * dBtdAzim[j] * errINCAZ[j];
      BpSigma2 = BpSigma2 + dBpdBtotal[j] * dBpdBtotal[j] * errBT[j] + dBpdFill[j] * dBpdFill[j] * errFILL[j] + 
                            dBpdInc[j] * dBpdInc[j] * errINC[j] + dBpdAzim[j] * dBpdAzim[j] * errAZ[j] + 
                            2.0 * dBpdBtotal[j] * dBpdFill[j] * errBTFILL[j] + 2.0 * dBpdBtotal[j] * dBpdInc[j] * errBTINC[j] +
                            2.0 * dBpdBtotal[j] * dBpdAzim[j] * errBTAZ[j] + 2.0 * dBpdFill[j] * dBpdInc[j] * errFILLINC[j] +
                            2.0 * dBpdFill[j] * dBpdAzim[j] * errFILLAZ[j] + 2.0 * dBpdInc[j] * dBpdAzim[j] * errINCAZ[j];
      BrSigma2 = BrSigma2 + dBrdBtotal[j] * dBrdBtotal[j] * errBT[j] + dBrdFill[j] * dBrdFill[j] * errFILL[j] + 
                            dBrdInc[j] * dBrdInc[j] * errINC[j] + dBrdAzim[j] * dBrdAzim[j] * errAZ[j] + 
                            2.0 * dBrdBtotal[j] * dBrdFill[j] * errBTFILL[j] + 2.0 * dBrdBtotal[j] * dBrdInc[j] * errBTINC[j] +
                            2.0 * dBrdBtotal[j] * dBrdAzim[j] * errBTAZ[j] + 2.0 * dBrdFill[j] * dBrdInc[j] * errFILLINC[j] +
                            2.0 * dBrdFill[j] * dBrdAzim[j] * errFILLAZ[j] + 2.0 * dBrdInc[j] * dBrdAzim[j] * errINCAZ[j];
  }
  *BrVariance = BrSigma2;
  *BtVariance = BtSigma2;
  *BpVariance = BpSigma2;
  return 0;
}

/*
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/patch/apps/errorprop.c,v $
$Author: xudong $
*/

