#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <endian.h>
#include <string.h>
#include <math.h>
#include <drms_types.h>
#include <printk.h>
#include "imgdecode_iris.h"

#define MASK(n) 			((1u << (n)) - 1u)
#define MAX(a,b)			((a) < (b) ? (b) : (a))
#define MIN(a,b)			((a) > (b) ? (b) : (a))

static int IMGX, IMGY, npixels;

/////////////////////////////////
int imgstat_iris(IMG *img, STAT *stat)
/////////////////////////////////
{
    double s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, t, ss;
    unsigned *h = img->hist;
    unsigned n = img->datavals;
    unsigned u;
    int i;

    if (n == 0) {
        stat->min = DRMS_MISSING_SHORT;
        stat->max = DRMS_MISSING_SHORT;
        stat->median = DRMS_MISSING_SHORT;
        stat->mean = DRMS_MISSING_DOUBLE;
        stat->rms = DRMS_MISSING_DOUBLE;
        stat->skew = DRMS_MISSING_DOUBLE;
        stat->kurt = DRMS_MISSING_DOUBLE;
	return 0;
    }

    if (img->reopened) {
	//
	// This image has been reopened.  Need to reconstruct histogram.
	//
	for (i = 0; i < MAXHIST; ++i)
	    h[i] = 0;
	n = 0;
	npixels = img->nx * img->ny;
	for (i = 0; i < npixels; ++i) {
	    if (img->dat[i] == BLANK)
		continue;
	    ++h[img->dat[i]];
	    ++n;
	}
	img->datavals = n;
    }

    memset(stat, 0, sizeof(STAT));

    while (h[stat->min] == 0)
	++stat->min;

    stat->max = MAXHIST - 1;
    while (h[stat->max] == 0)
	--stat->max;

    stat->median = stat->min;
    i = h[stat->median];
    while (2*i < n)
	i += h[++stat->median];

    for (i = stat->min; i <= stat->max; ++i) {
	double x = i;
	s += (t = x*h[i]);
	s2 += (t *= x);
	s3 += (t *= x);
	s4 += (t *= x);
    }

    s /= n; 
    ss = s*s;
    stat->mean = s;
    s2 /= n; 
    s3 /= n; 
    s4 /= n;
    if (n > 1) {
	t = n * (s2 - ss) / (n-1);
	stat->rms = sqrt(t);
    }
    if (stat->rms > 0.0) {
	stat->skew = (s3 - s * (3*s2 - 2*ss)) / (t*stat->rms);
	stat->kurt = (s4 - 4*s*s3 + 3*ss*(2*s2-ss)) / (t*t) - 3;
    }

    return 0;
}

//////////////////////////////////////////////////////////////////////////////
static void put(short *dat, int TAP, int r, int c, int n, unsigned short *pix)
//////////////////////////////////////////////////////////////////////////////
{
    int i;

    switch (TAP) {
    case 0:	// 4-port: E,F,G,H
	if (r < IMGY/2)
	    memcpy(dat+IMGX*r+c, pix, 2*n);
	else if (r < IMGY)
	    for (i=c; i<c+n; ++i)
		dat[IMGX*(r-IMGY/2)+IMGX-1-i] = *pix++;
	else if (r < IMGY+IMGY/2)
	    for (i=c; i<c+n; ++i)
		dat[IMGX*(2*IMGY-1-r)+IMGX-1-i] = *pix++;
	else
	    memcpy(dat+IMGX*(2*IMGY+IMGY/2-1-r)+c, pix, 2*n);
	break;
    case 1:	// 2-port: E,F
    case 9:	// like TAP 0 but only output E and F quadrants
	if (r < IMGY)
	    memcpy(dat+IMGX*r+c, pix, 2*n);
	else
	    for (i=c; i<c+n; ++i)
		dat[IMGX*(r-IMGY)+IMGX-1-i] = *pix++;
	break;
    case 2:	// 2-port: F,G
	if (r < IMGY/2)
	    for (i=c; i<c+n; ++i)
		dat[IMGX*r+IMGX-1-i] = *pix++;
	else
	    for (i=c; i<c+n; ++i)
		dat[IMGX*(IMGY+IMGY/2-1-r)+IMGX-1-i] = *pix++;
	break;
    case 3:	// 2-port: G,H
    case 10:	// like TAP 0 but only output G and H quadrants
	if (r < IMGY)
	    for (i=c; i<c+n; ++i)
		dat[IMGX*(IMGY-1-r)+IMGX-1-i] = *pix++;
	else
	    memcpy(dat+IMGX*(2*IMGY-1-r)+c, pix, 2*n);
	break;
    case 4:	// 2-port: E,H (NOT H,E!)
	if (r < IMGY/2)
	    memcpy(dat+IMGX*r+c, pix, 2*n);
	else
	    memcpy(dat+IMGX*(IMGY+IMGY/2-1-r)+c, pix, 2*n);
	break;
    case 5:	// 1-port: E
	memcpy(dat+IMGX*r+c, pix, 2*n);
	break;
    case 6:     // 1-port: F
	for (i=c; i<c+n; ++i)
	    dat[IMGX*r+IMGX-1-i] = *pix++;
	break;
    case 7:	// 1-port: G
	for (i=c; i<c+n; ++i)
	    dat[IMGX*(IMGY-1-r)+IMGX-1-i] = *pix++;

	break;
    case 8:	// 1-port: H
	memcpy(dat+IMGX*(IMGY-1-r)+c, pix, 2*n);
	break;
    }
}

/////////////////////////////////
int imgdecode_iris_init_hack(IMG *img)
/////////////////////////////////
{
    int i;

    img->datavals = 0;
    img->npackets = 0;
    img->nerrors = 0;
    img->last_pix_err = 0;
    img->first_packet_time = UINT64_MAX;
    for (i = 0; i < MAXPIXELS; ++i)
	img->dat[i] = BLANK;
    for (i = 0; i < MAXHIST; ++i)
	img->hist[i] = 0;
}

//////////////////////////////////////////////
int imgdecode_iris(unsigned short *impdu, IMG *img)
//////////////////////////////////////////////
{
    unsigned short old, final, diff, w, h, skp, tak;
    int i, N, K, bits2go, wordcnt, nzero;
    unsigned u, bitbuf, sgn, low, fs, kmask, nmask, nmk;
    unsigned offset, ndecoded;
    uint64_t uu;
    const unsigned short *p;
    static unsigned short pix[7000];	// max no. of pixels per 869-word
					// packet is (869-2) * (16/2) = 6936
    static unsigned short *lut[256];
    static CROPTABLE cropt[4096];
    char fname[64];
    FILE *fp;

    //
    // number of trailing zero bits in 0 .. 255
    //
    static char nz[256] = {
	0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0 
    };

    //
    // img->nx,ny is for a thin, horizontal image
    //  IMGX,Y is for a thin, vertical image
    // 
    IMGX = img->ny;
    IMGY = img->nx;
    npixels = IMGX * IMGY;
	
//
// swap bytes if necessary
//
#if __BYTE_ORDER == __LITTLE_ENDIAN
    for (i=0; i<PACKETWORDS; ++i)
	impdu[i] = (impdu[i]>>8) + (impdu[i]<<8);
#endif

//
// initialize img struct and read various tables.
//
    if (img->initialized) 
	goto ___DECODE_START___;

    u = impdu[4] & 0x7ffu;	// APID
    if ( u != APID_IRIS_SCIENCE) {
	return IMGDECODE_BAD_APID;
    }
    img->apid = u;
    img->telnum = impdu[11] >> 14;
    img->fsn = ((impdu[11] & 0x3fff) << 16) + impdu[12];
    img->isysn = impdu[11] >> 14;
    //img->fid = ((impdu[13] & 0xff) << 16) + impdu[14];
    img->fid = ((impdu[14] & 0xff) << 16) + impdu[13];
    img->cropid = impdu[15] >> 4;
    img->overflow = impdu[15] & 1;
    img->headerr = (impdu[15] >> 1) & 1;
    img->luid = impdu[17] >> 8;
    img->tap = impdu[16] >> 12;
    if (img->tap == 0 && img->cropid > 0)
	if (img->isysn == 2)
	    img->tap = 9;
	else if (img->isysn == 1)
	    img->tap = 10;
    u = impdu[16] & 0xff;	// compression ID
    if (u == 0 || u == 128) {
	// raw mode
	img->N = 16;
	img->K = 0;
	img->R = 0;
    } else {
	// compressed mode
	img->N = u >> 3;
	if (img->N > 14)
	    return IMGDECODE_BAD_N;
	img->K = u & 0x7;
	img->R = (impdu[16] >> 8) & 0xf;
    }

// TEMP
//Removed 30May2012 See K.C. mail 11:45
//        img->N = 16;
//        img->K = 0;
//        img->R = 0;


// The following has been moved to imgdecode_init_hack()
//#if 0
    img->datavals = 0;
    img->npackets = 0;
    img->nerrors = 0;
    img->last_pix_err = 0;
    img->first_packet_time = UINT64_MAX;

    for (i = 0; i < MAXPIXELS; ++i)
	img->dat[i] = BLANK;

    for (i = 0; i < MAXHIST; ++i)
	img->hist[i] = 0;
//#endif

    img->initialized = 1;

___DECODE_START___:

    //
    // read inverse lookup table
    //
    if (img->luid && !lut[img->luid]) {
	if (isIRIS(img->apid))
	    snprintf(fname, 64, TABLE_DIR "/lu/iris/ilu%u", img->luid);
	else
	    return IMGDECODE_BAD_APID;
	fp = fopen(fname, "r");
	if (!fp) 
	    return IMGDECODE_NO_LOOKUP_TABLE;
	fscanf(fp, "%u", &u);
	if (u != img->luid)
	    return IMGDECODE_LOOKUP_ID_MISMATCH;
	lut[img->luid] = (unsigned short *) malloc(16384*2);
	if (!lut[img->luid])
	    return IMGDECODE_OUT_OF_MEMORY;
	for (i=0; i<16384; ++i)
	    if (1 != fscanf(fp, "%hu", lut[img->luid] + i)) {
		fclose(fp);
		free(lut[img->luid]);
		lut[img->luid] = 0;
		return IMGDECODE_BAD_LOOKUP_TABLE;
	    }
	fclose(fp);
    }

    //
    // read crop table
    //
    if (cropt[img->cropid].totalpix)
	img->totalvals = cropt[img->cropid].totalpix;
    else if (img->cropid) {
	if (isIRIS(img->apid))
	    snprintf(fname, 64, TABLE_DIR "/crop/iris/crop%u", img->cropid);
	else
	    return IMGDECODE_BAD_APID;
	fp = fopen(fname, "r");
	if (!fp)
	    return IMGDECODE_NO_CROP_TABLE;
	fscanf(fp, "%u", &u);
	if (u != img->cropid)
	    return IMGDECODE_CROP_ID_MISMATCH;
	fscanf(fp, "%hu %hu", &w, &h);
	if (!(w == IMGX && h == IMGY) && !(w == IMGX/2 && h == IMGY*2)) {
printk("bad_geo: w=%u  h=%u  IMGX=%d  IMGY=%d\n", w, h, IMGX, IMGY);
printk("bad_geo: cropid=%d  fsn=%u\n", img->cropid, img->fsn);
	    return IMGDECODE_BAD_CROP_GEOMETRY;
        }
	cropt[img->cropid].width = w;
	cropt[img->cropid].height = h;
	cropt[img->cropid].skip = (unsigned short *) malloc(2*h);
	cropt[img->cropid].take = (unsigned short *) malloc(2*h);
	cropt[img->cropid].offset = (unsigned *) malloc(4*h);
	if (!cropt[img->cropid].offset)
	    return IMGDECODE_OUT_OF_MEMORY;
	offset = 0;
	for (i=0; i<h; ++i) {
	    if (2 != fscanf(fp, "%hu %hu", &skp, &tak)) {
		fclose(fp);
		free(cropt[img->cropid].skip);
		free(cropt[img->cropid].take);
		free(cropt[img->cropid].offset);
		return IMGDECODE_BAD_CROP_TABLE;
	    }
	    if (skp + tak > w)
		return IMGDECODE_BAD_CROP_SKIP_TAKE;
	    cropt[img->cropid].skip[i] = skp;
	    cropt[img->cropid].take[i] = tak;
	    cropt[img->cropid].offset[i] = offset;
	    offset += tak;
	}
	fclose(fp);
	cropt[img->cropid].totalpix = offset;
	img->totalvals = cropt[img->cropid].totalpix;
    } else
	img->totalvals = npixels;

    //
    // update first_packet_time if necessary
    //
    uu = ((uint64_t)impdu[7] << 32) + ((uint64_t)impdu[8] << 16) + impdu[9];
    if (uu < img->first_packet_time)
	img->first_packet_time = uu;

    //
    // check overflow and header error flags
    //
    if (!img->overflow)
	img->overflow = impdu[15] & 1;
    if (!img->headerr) 
	img->headerr = (impdu[15] >> 1) & 1;


    //
    // where does this packet start?
    //
    offset = ((impdu[17] & 0xff) << 16) + impdu[18];
    if (offset >= img->totalvals)
	return IMGDECODE_BAD_OFFSET;

    p = impdu + PACKETHEADERWORDS;

    //
    // raw mode
    //
    if (img->N == 16) {
	for (i = 0; i < PACKETDATAWORDS && img->totalvals > offset+i; ++i)
	    pix[i] = p[i] & 0x3fff;	// 14 bits only !!!
	ndecoded = i;
    }

    //
    // compressed mode
    //
    else {
	old = pix[0] = *p++;
	ndecoded = 1;
	bitbuf = *p++;
       	bits2go = 16;
	K = img->K;
	N = img->N;
	nmk = N - K;
	kmask = MASK(K);
	nmask = MASK(nmk);
	final = impdu[PACKETWORDS - 1];

	wordcnt = PACKETDATAWORDS - 3;
	for (i = PACKETWORDS - 2; i > PACKETHEADERWORDS + 1; --i) {
	    if (impdu[i]) break;
	    --wordcnt;
	}

	// At this point wordcnt is the number of words, excluding the
	// zero-fills, yet to be read.  
	// ===== However, when the final pixel saturates and has some 
	// ===== leading zero bits in its top N-K bits, we may have
	// ===== mistaken a non-fill word which happens to be zero for 
	// ===== a zero-fill.  This is taken care of below. 

	while (wordcnt > 0 || (bits2go && bitbuf)) {
	    if (bits2go < 2 + K) {
		bitbuf += *p++ << bits2go;
		--wordcnt;
		bits2go += 16;
	    }

	    //
	    // decode low K+1 bits
	    //
	    low = bitbuf & kmask;
	    bitbuf >>= K;
	    sgn = bitbuf & 0x1;
	    bitbuf >>= 1;
	    bits2go -= 1 + K;

	    //
	    // decode FS
	    //
	    if (bitbuf == 0) {
		if (bits2go < 9) {
		    fs = bits2go;
		    if (wordcnt > 0) {
			bitbuf = *p++;
			--wordcnt;
			bits2go = 16;
		    } else
			goto __DECOMPRESS_FAILURE__;
		    if (bitbuf & 0xff) {
			nzero = nz[bitbuf & 0xff];
			fs += nzero;
			if (fs > 8)
			    goto __DECOMPRESS_FAILURE__;
		    } else
			goto __DECOMPRESS_FAILURE__;
		} else
		    goto __DECOMPRESS_FAILURE__;
	    } else if (bitbuf & 0xff) {
		nzero = nz[bitbuf & 0xff];
		fs = nzero;
	    } else if (bitbuf & 0x1ff) {
		nzero = 8;
		fs = nzero;
	    } else
		goto __DECOMPRESS_FAILURE__;
	    bitbuf >>= (nzero + 1);
	    bits2go -= nzero + 1;

	    //
	    // SATURATION case
	    //
	    if (fs == 8) {
		if (bits2go < nmk) {
		    u = nmk - bits2go;
		    fs = bitbuf;

		    // The condition is (wordcnt >= 0), NOT (wordcnt > 0)!
		    // See comments above with ===== mark.
		    if (wordcnt >= 0) {
			bitbuf = *p++;
			--wordcnt;
			fs += (bitbuf & MASK(u)) << bits2go;
			bitbuf >>= u;
			bits2go = 16 - u;
		    } else
			goto __DECOMPRESS_FAILURE__;
		} else {
		    fs = bitbuf & nmask;
		    bitbuf >>= nmk;
		    bits2go -= nmk;
		}
	    }

	    //
	    // be done with this pixel
	    //
	    diff = (fs << K) + low;
	    if (sgn)
		old -= diff;
	    else
		old += diff;
	    pix[ndecoded++] = old;
	}

	if (old != final)
	    goto __DECOMPRESS_FAILURE__;
    }

    goto __DECOMPRESS_SUCCESS__;

__DECOMPRESS_FAILURE__:

    ++img->nerrors;
    u = offset + ndecoded;
    if (u < img->totalvals && u > img->totalvals - 4) {
	img->last_pix_err = 1;
	img->datavals += ndecoded;
	++img->npackets;
	goto __POST_PROCESSING__;
    }
    return IMGDECODE_DECOMPRESS_ERROR;

__DECOMPRESS_SUCCESS__:

    u = img->datavals + ndecoded;
    if (u > img->totalvals && !img->reopened) {
	++img->nerrors;
	return IMGDECODE_TOO_MANY_PIXELS;
    }
    if (offset + ndecoded > img->totalvals) {
	++img->nerrors; 
	return IMGDECODE_BAD_OFFSET;
    }
    img->datavals = u;
    ++img->npackets;

__POST_PROCESSING__:

    for (i = 0; i < ndecoded; ++i) {
	//
	// inverse lookup
	//
	if (img->luid)
	    pix[i] = lut[img->luid][pix[i]];
	//
	// undo right shift
	//
	if (img->R)
	    pix[i] <<= img->R;
	//
	// update histogram
	//
	if (img->hist)
	    ++img->hist[pix[i]];
    }

//
// uncrop and reconstruct image
//
    if (img->cropid) {
	unsigned *o = cropt[img->cropid].offset;
	unsigned short *s = cropt[img->cropid].skip;
	unsigned short *t = cropt[img->cropid].take;
	int done = 0;
	h = cropt[img->cropid].height;
	for (i = 0; ndecoded > 0 && i < h; ++i) {
	    if (o[i] > offset || (h > i+1 && o[i+1] <= offset))
		continue;
	    skp = offset - o[i];
	    tak = MIN(ndecoded, t[i] - skp);
	    if (tak) {
		skp += s[i];
		put(img->dat, img->tap, i, skp, tak, pix+done);
		ndecoded -= tak;
		offset += tak;
		done += tak;
	    }
	}
    } else {
	int done = 0;
	w = IMGX;
	h = IMGY;
	switch(img->tap) {
	case 0: 
	case 1: 
	case 3: 
	case 9:
	case 10:
	    w = IMGX/2; h = IMGY*2; 
	    break;
	}
	i = offset / w;
	while (ndecoded) {
	    skp = offset % w;
	    tak = MIN(ndecoded, w - skp);
	    if (tak) {
		put(img->dat, img->tap, i, skp, tak, pix+done);
		ndecoded -= tak;
		offset += tak;
		done += tak;
		++i;
	    }
	}
    }

    return 0;
}
