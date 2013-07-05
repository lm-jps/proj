#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/libs/stats/dstats.c,v 1.4 2013/07/05 18:35:19 phil Exp $"

#include <stdlib.h>
#include <math.h>
#include <float.h>

//
// med - returns median of an n-element double array arr[]
//       based on Numerical Recipes select()
//       Note1: NaN's must be removed first.
//       Note2: This function rearranges elements of arr[] on exit.
//       Note3: For even n, this functions returns element n/2 in the
//              sorted array.
//       Note4: +- inf must also be removed.

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

static double median(int n, double arr[])
{
    int k, i, ir, j, l, mid;
    double a, temp;

    if (n < 1)
	return __builtin_nan("");

    k = n / 2;
    l = 0;
    ir = n - 1;
    for (;;) {
	if (ir <= l + 1) {
	    if (ir == l + 1 && arr[ir] < arr[l]) {
		SWAP(arr[l], arr[ir])
	    }
	    return arr[k];
	} else {
	    mid = (l + ir) >> 1;
	    SWAP(arr[mid], arr[l + 1])
		if (arr[l] > arr[ir]) {
		SWAP(arr[l], arr[ir])
	    }
	    if (arr[l + 1] > arr[ir]) {
		SWAP(arr[l + 1], arr[ir])
	    }
	    if (arr[l] > arr[l + 1]) {
		SWAP(arr[l], arr[l + 1])
	    }
	    i = l + 1;
	    j = ir;
	    a = arr[l + 1];
	    for (;;) {
		do
		    i++;
		while (arr[i] < a);
		do
		    j--;
		while (arr[j] > a);
		if (j < i)
		    break;
		SWAP(arr[i], arr[j])
	    }
	    arr[l + 1] = arr[j];
	    arr[j] = a;
	    if (j >= k)
		ir = j - 1;
	    if (j <= k)
		l = i;
	}
    }
}

#define OK 0
#define OUT_OF_MEMORY -1
#define TOO_FEW_GOOD_POINTS -2

int dstats(int n, double arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood)
{
    int i;
    int nv = 0;
    double *dat, dmin = DBL_MAX, dmax = -DBL_MAX;
    double s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, avg, var;

    *min = *max = *medn = *mean = *sig = *skew = *kurt = __builtin_nan("");

    if (!(dat = (double *) malloc(n * sizeof(double))))
	return OUT_OF_MEMORY;

    for (i = 0; i < n; ++i) {
	double t = arr[i];
        int fpclass = fpclassify(t);
	if (fpclass == FP_NAN || fpclass == FP_INFINITE)
	    continue;
	dat[nv++] = t;
	if (dmin > t)
	    dmin = t;
	if (dmax < t)
	    dmax = t;
	s += t;
    }
    *ngood = nv;

    if (nv < 2) {
	if (nv == 1)
	    *min = *max = *medn = *mean = s;
	free(dat);
	return TOO_FEW_GOOD_POINTS;
    }

    if (dmin == dmax) {
	*min = *max = *medn = *mean = dmin;
	*sig = *skew = *kurt = 0;
	free(dat);
	return OK;
    }

    avg = s / nv;
    for (i = 0; i < nv; ++i) {
	double tt;
	double t = dat[i] - avg;
	tt = t * t;
	s2 += tt;
	tt *= t;
	s3 += tt;
	tt *= t;
	s4 += tt;
    }

    *min = dmin;
    *max = dmax;
    *medn = median(nv, dat);
    *mean = avg;
    var = s2 / (nv - 1);
    *sig = sqrt(var);
    *skew = s3 / (var * *sig * nv);
    *kurt = s4 / (var * var * nv) - 3.0;

    free(dat);
    return OK;
}
