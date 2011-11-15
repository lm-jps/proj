#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/libs/stats/dstats2.c,v 1.3 2011/11/15 20:35:19 kehcheng Exp $"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define NBINS 65536
#define OK 0
#define TOO_FEW_GOOD_POINTS -2

int dstats2(int n, double arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood)
{
    int i, bin, nv = 0;
    double dmin = DBL_MAX, dmax = -DBL_MAX;
    double s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, avg, var;
    static int hist[NBINS];
    double delta;

    *min = *max = *medn = *mean = *sig = *skew = *kurt = __builtin_nan("");

    for (i = 0; i < n; ++i) {
	double t = arr[i];
	if (isnan(t))
	    continue;
	if (dmin > t)
	    dmin = t;
	if (dmax < t)
	    dmax = t;
	s += t;
	++nv;
    }
    *ngood = nv;

    if (nv < 2) {
	if (nv == 1)
	    *min = *max = *medn = *mean = s;
	return TOO_FEW_GOOD_POINTS;
    }

    if (dmin == dmax) {
	*min = *max = *medn = *mean = dmin;
	*sig = *skew = *kurt = 0;
	return OK;
    }

    delta = (dmax - dmin) / NBINS;
    memset(hist, 0, sizeof(int)*NBINS);
    avg = s / nv;
    for (i = 0; i < n; ++i) {
	double tmp = arr[i];
	double t = tmp - avg;
	double tt;
	if (isnan(tmp))
	    continue;
	tt = t * t;
	s2 += tt;
	tt *= t;
	s3 += tt;
	tt *= t;
	s4 += tt;
	bin = (tmp - dmin) / delta;
	if (bin >= NBINS)
	    bin = NBINS - 1;
	++hist[bin];
    }

    bin = 0;
    n = hist[bin];
    while (2*n < nv)
	n += hist[++bin];
    *medn = dmin + delta * (bin + 0.5);

    *min = dmin;
    *max = dmax;
    *mean = avg;
    var = s2 / (nv - 1);
    *sig = sqrt(var);
    *skew = s3 / (var * *sig * nv);
    *kurt = s4 / (var * var * nv) - 3.0;

    return OK;
}
