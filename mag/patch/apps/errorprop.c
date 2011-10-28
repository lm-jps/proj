#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int errorprop (float *bTotal, float *bAzim, float *bInc, float *ebT, float *ebA, float *ebI,
               float *ebTbA, float *ebTbI, float *ebIbA, 
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
  
  // Xudong Oct 18 2011: Fill factor variances covariances are all NaNs. removed
  // NNB interpolation, just 1 point
  // Fixed definition of azi for all derivatives

  static double raddeg = M_PI / 180.;
  double b, inc, azim;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double dBrdBtotal, dBrdInc, dBrdAzim;
  double dBtdBtotal, dBtdInc, dBtdAzim;
  double dBpdBtotal, dBpdInc, dBpdAzim;
  double errBT, errINC, errAZ;
  double errBTINC, errBTAZ, errINCAZ;
  double BrSigma2, BtSigma2, BpSigma2;
  int ix, iy;

  if (x < 1. || x >= (float)(nx-2) || y < 1. || y >= (float)(ny-2))
  return (1);
  
  a11 = - sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - lonc);
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
  iy = (int)y;

  int iAll = iy * nx + ix;

  if (isnan(bTotal[iAll])) return(1);
  b = bTotal[iAll]; /* field strength from the pixels used for the interpolation */
  if (isnan(bInc[iAll])) return(1);
  inc = raddeg * bInc[iAll]; /* inclination, in deg */
  if (isnan(bAzim[iAll])) return(1);
  azim = raddeg * bAzim[iAll]; /* azimuth, in deg */
  if (isnan(ebT[iAll])) return(1);
  errBT = ebT[iAll]; /* variance of field strength, in G^2 */
  if (isnan(ebI[iAll])) return(1);
  errINC = ebI[iAll]; /* variance of inclination, in rad^2 */
  if (isnan(ebA[iAll])) return(1);
  errAZ = ebA[iAll]; /* variance of azimuth, in rad^2 */
  if (isnan(ebTbI[iAll])) return(1);
  errBTINC = ebTbI[iAll]; /* covariance of field and inclination, in G.rad */
  if (isnan(ebTbA[iAll])) return(1);
  errBTAZ = ebTbA[iAll];  /* covariance of field and azimuth, in G.rad */
  if (isnan(ebIbA[iAll])) return(1);
  errINCAZ = ebIbA[iAll]; /* covariance of inclination and azimuth, in rad^2 */


  dBpdBtotal = (- a11 * sin(inc) * sin(azim) + a12 * sin(inc) * cos(azim) + a13 * cos(inc));
  dBpdInc = b * (- a11 * cos(inc) * sin(azim) + a12 * cos(inc) * cos(azim) - a13 * sin(inc));
  dBpdAzim = b * (- a11 * sin(inc) * cos(azim) - a12 * sin(inc) * sin(azim));
  
  dBtdBtotal = (- a21 * sin(inc) * sin(azim) + a22 * sin(inc) * cos(azim) + a23 * cos(inc));
  dBtdInc = b * (- a21 * cos(inc) * sin(azim) + a22 * cos(inc) * cos(azim) - a23 * sin(inc));
  dBtdAzim = b * (- a21 * sin(inc) * cos(azim) - a22 * sin(inc) * sin(azim));
  
  dBrdBtotal = (- a31 * sin(inc) * sin(azim) + a32 * sin(inc) * cos(azim) + a33 * cos(inc));
  dBrdInc = b * (- a31 * cos(inc) * sin(azim) + a32 * cos(inc) * cos(azim) - a33 * sin(inc));
  dBrdAzim = b * (- a31 * sin(inc) * cos(azim) - a32 * sin(inc) * sin(azim));

  BrSigma2 = 0.0;
  BtSigma2 = 0.0;
  BpSigma2 = 0.0;

  BtSigma2 = BtSigma2 + dBtdBtotal * dBtdBtotal * errBT + 
                        dBtdInc * dBtdInc * errINC + 
                        dBtdAzim * dBtdAzim * errAZ + 
                        2.0 * dBtdBtotal * dBtdInc * errBTINC +
                        2.0 * dBtdBtotal * dBtdAzim * errBTAZ +
                        2.0 * dBtdInc * dBtdAzim * errINCAZ;
  BpSigma2 = BpSigma2 + dBpdBtotal * dBpdBtotal * errBT +
                        dBpdInc * dBpdInc * errINC + 
                        dBpdAzim * dBpdAzim * errAZ + 
                        2.0 * dBpdBtotal * dBpdInc * errBTINC +
                        2.0 * dBpdBtotal * dBpdAzim * errBTAZ +
                        2.0 * dBpdInc * dBpdAzim * errINCAZ;
  BrSigma2 = BrSigma2 + dBrdBtotal * dBrdBtotal * errBT + 
                        dBrdInc * dBrdInc * errINC + 
                        dBrdAzim * dBrdAzim * errAZ + 
                        2.0 * dBrdBtotal * dBrdInc * errBTINC +
                        2.0 * dBrdBtotal * dBrdAzim * errBTAZ +
                        2.0 * dBrdInc * dBrdAzim * errINCAZ;

  *BrVariance = BrSigma2;
  *BtVariance = BtSigma2;
  *BpVariance = BpSigma2;
  
  return 0;
}



