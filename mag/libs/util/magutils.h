#ifndef __MAGUTILS_H__
#define __MAGUTILS_H__

#include "drms.h"

#ifdef __cplusplus
extern "C" {
#endif

void imagefromchebyshev(double *image, int m, int n, int order, double coef[], double *xdelta, double *ydelta);
int img2helioVector (double bxImg, double byImg, double bzImg, double *bxHelio, double *byHelio, double *bzHelio, double lon, double lat, double lonc, double latc, double pAng);
int noisemaskimag4twindow(int xDim, int yDim, float xcen, float ycen, float rsun, float vrcenter, int maskid, float *image);
int obstime2maskid(TIME tobs);

#ifdef __cplusplus
}
#endif

#endif /* __MAGUTILS_H__ */
