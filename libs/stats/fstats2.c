#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/libs/stats/fstats2.c,v 1.2 2010/07/21 00:47:35 phil Exp $"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define NBINS 65536
#define OK 0
#define TOO_FEW_GOOD_POINTS -2

int fstats2(int n, float arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood)
{
    int i, bin, nv = 0;
    float fmin = FLT_MAX, fmax = FLT_MIN;
    double s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, avg, var;
    static int hist[NBINS];
    double delta;

    *min = *max = *medn = *mean = *sig = *skew = *kurt = __builtin_nan("");

    for (i = 0; i < n; ++i) {
	float t = arr[i];
	if (isnanf(t))
	    continue;
	if (fmin > t)
	    fmin = t;
	if (fmax < t)
	    fmax = t;
	s += t;
	++nv;
    }
    *ngood = nv;

    if (nv < 2) {
	if (nv == 1)
	    *min = *max = *medn = *mean = s;
	return TOO_FEW_GOOD_POINTS;
    }

    if (fmin == fmax) {
	*min = *max = *medn = *mean = fmin;
	*sig = *skew = *kurt = 0;
	return OK;
    }

    delta = (fmax - fmin) / NBINS;
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
	bin = (tmp - fmin) / delta;
	if (bin >= NBINS)
	    bin = NBINS - 1;
	++hist[bin];
    }

    bin = 0;
    n = hist[bin];
    while (2*n < nv)
	n += hist[++bin];
    *medn = fmin + delta * (bin + 0.5);

    *min = fmin;
    *max = fmax;
    *mean = avg;
    var = s2 / (nv - 1);
    *sig = sqrt(var);
    *skew = s3 / (var * *sig * nv);
    *kurt = s4 / (var * var * nv) - 3.0;

    return OK;
}
