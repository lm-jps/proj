//Place holder only to get the make build_lev1_iris working
#include <drms.h>
#include "lev0lev1.h"

#define XXX_BAD_GEOMETRY -1
#define XXX_BAD_PIXCOUNT -2

#define NBINS 1048576
#define MINOUT -256
#define MAXOUT 1048320

static int hist[NBINS];
static float dark2[2270912];	// 2072*1096 = 2270912
static float flat2[2270912];	// 2072*1096 = 2270912

////////////////////////////////
int do_flat_iris(LEV0LEV1 *info)
{
    int is_dark = info->darkflag;
    float *flat = info->adataff;
    float *dark = info->adatadark;
    short *in = info->adata0;
    int *out = info->dat1.adata1A;
    int datavals = info->datavals;
    int missvals = info->missvals;
    int totalvals = datavals + missvals;
    int nx = info->nx;
    int ny = info->ny;
    int sx = info->sumx;
    int sy = info->sumy;
    int totalpix = nx * ny;
    int ndark, nflat;
    int min, max, medn;
    short tmp;
    int itmp;
    float ftmp,flatsum,darksum;
    double dtmp, s, s2, s3, s4, ss;
    int i, j, k, l;
    int idx, IDX;

    // sanity checks
    if (sx < 1 || sy < 1 || (nx*sx != 4144 && nx*sx != 2072) || ny*sy != 1096)
	return XXX_BAD_GEOMETRY;
    if (totalpix < totalvals)
	return XXX_BAD_PIXCOUNT;

    //
    // bin down dark and flatfield
    //
    if (sx > 1 || sy > 1) {
	for (j=0; j<ny; ++j) {
	    for (i=0; i<nx; ++i) {
		idx = j*nx + i;
		nflat = ndark = sx*sy;
		flatsum = darksum = 0;
		for (l=j*sy; l<(j+1)*sy; ++l) {
		    for (k=i*sx; k<(i+1)*sx; ++k) {
			IDX = l*nx*sx + k;
			if (isnan(dark[IDX]))
			    --ndark;
			else
			    darksum += dark[IDX];
			if (isnan(flat[IDX]))
			    --nflat;
			else
			    flatsum += flat[IDX];
		    }
		}
		dark2[idx] = ndark ? darksum/ndark : DRMS_MISSING_FLOAT;
		flat2[idx] = nflat ? flatsum/nflat : DRMS_MISSING_FLOAT;
	    }
	}
    }

    //
    // apply dark and flatfield correction
    //
    memset(hist, 0, 4*NBINS);
    datavals = 0;
    for (i=0; i<totalpix; ++i) {
	tmp = in[i];
	if (DRMS_MISSING_SHORT == tmp) {
	    out[i] = DRMS_MISSING_INT;
	} else {
	    if (sx == 1 && sy == 1)
		ftmp = is_dark ? tmp : (tmp - dark[i]) / flat[i];
	    else
		ftmp = is_dark ? tmp : (tmp - dark2[i]) / flat2[i];
	    if (ftmp >= MINOUT && ftmp < MAXOUT) {
		itmp = roundf(ftmp);
		out[i] = itmp;
		++hist[itmp-MINOUT];
		++datavals;
	    } else {
		out[i] = DRMS_MISSING_INT;
	    }
	}
    }
    info->datavals = datavals;
    info->missvals = totalvals - datavals;

    //
    // do flip
    //
    if (datavals) {
	switch(info->winflip) {
	case 1:	// top-bottom flip
	    {
		int tmp[4144];
		for (j=0; j<ny/2; ++j) {
		    memcpy(tmp, &out[j*nx], 4*nx);
		    memcpy(&out[j*nx], &out[(ny-j-1)*nx], 4*nx);
		    memcpy(&out[(ny-j-1)*nx], tmp, 4*nx);
		}
	    }
	    break;
	case 2:	// left-right flip
	    {
		int tmp2;
		for (j=0; j<ny; ++j)
		    for (i=0; i<nx/2; ++i) {
			tmp2 = out[j*nx+i];
			out[j*nx+i] = out[j*nx+(nx-1-i)];
			out[j*nx+(nx-1-i)] = tmp2;
		    }
	    }
	    break;
	case 3:	// both
	    {
		int tmp[4144], tmp2;
		for (j=0; j<ny/2; ++j) {
		    memcpy(tmp, &out[j*nx], 4*nx);
		    memcpy(&out[j*nx], &out[(ny-j-1)*nx], 4*nx);
		    memcpy(&out[(ny-j-1)*nx], tmp, 4*nx);
		}
		for (j=0; j<ny; ++j)
		    for (i=0; i<nx/2; ++i) {
			tmp2 = out[j*nx+i];
			out[j*nx+i] = out[j*nx+(nx-1-i)];
			out[j*nx+(nx-1-i)] = tmp2;
		    }
	    }
	    break;
	}
    }


    //
    // do statistics
    //
    info->datamin = info->datamax = info->datamedn = DRMS_MISSING_INT;
    info->datamean = info->data_rms = info->dataskew = info->datakurt = DRMS_MISSING_DOUBLE;

    min = 0;
    max = NBINS-1;
    if (datavals) {
	while (hist[min] == 0)
	    ++min;
	info->datamin = min + MINOUT;

	while (hist[max] == 0)
	    --max;
	info->datamax = max + MINOUT;

	medn = min;
	i = hist[medn];
	while (2*i < datavals)
	    i += hist[++medn];
	info->datamedn = medn + MINOUT;

	s = s2 = s3 = s4 = 0.0;
	for (i=min; i<=max; ++i) {
	    double ii = i + MINOUT;
	    s += (dtmp = ii*hist[i]);
	    s2 += (dtmp *= ii);
	    s3 += (dtmp *= ii);
	    s4 += (dtmp *= ii);
	}
	s /= datavals;
	info->datamean = s;
	ss = s*s;
	s2 /= datavals;
	s3 /= datavals;
	s4 /= datavals;
	if (datavals > 1) {
	    dtmp = datavals * (s2-ss) / (datavals-1);
	    info->data_rms = sqrt(dtmp);
	    if (dtmp > 0.0) {
		info->dataskew = (s3 - s * (3*s2 - 2*ss)) / 
		    (dtmp*info->data_rms);
		info->datakurt = (s4 - 4*s*s3 + 3*ss*(2*s2-ss)) / 
		    (dtmp*dtmp) - 3;
	    }
	}
    }
    return 0;
}

