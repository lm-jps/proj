#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <endian.h>
#include <string.h>
#include <math.h>
#include "imgdecode.h"

#define MASK(n) 			((1u << (n)) - 1u)
#define MAX(a,b)			((a) < (b) ? (b) : (a))
#define MIN(a,b)			((a) > (b) ? (b) : (a))

/////////////////////////////////
int imgstat(IMG *img, STAT *stat)
/////////////////////////////////
{
    double s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, t, ss;
    unsigned *h = img->hist;
    unsigned n = img->datavals;
    unsigned u;
    int i;

    if (n == 0)
	return -1;

    if (img->reopened) {
	//
	// This image has been reopened.  Need to reconstruct histogram.
	//
	for (i = 0; i < MAXHIST; ++i)
	    h[i] = 0;
	n = 0;
	for (i = 0; i < MAXPIXELS; ++i) {
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
	s += (t = i*h[i]);
	s2 += (t *= i);
	s3 += (t *= i);
	s4 += (t *= i);
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
    case 0:	// 1-port: E,F,G,H
	if (r < 2048)
	    memcpy(dat+4096*r+c, pix, 2*n);
	else if (r < 4096)
	    for (i=c; i<c+n; ++i)
		dat[4096*(r-2048)+4095-i] = *pix++;
	else if (r < 6144)
	    for (i=c; i<c+n; ++i)
		dat[4096*(8191-r)+4095-i] = *pix++;
	else
	    memcpy(dat+4096*(10239-r)+c, pix, 2*n);
	break;
    case 1:	// 2-port: E,F
	if (r < 4096)
	    memcpy(dat+4096*r+c, pix, 2*n);
	else
	    for (i=c; i<c+n; ++i)
		dat[4096*(r-4096)+4095-i] = *pix++;
	break;
    case 2:	// 2-port: F,G
	if (r < 2048)
	    for (i=c; i<c+n; ++i)
		dat[4096*r+4095-i] = *pix++;
	else
	    for (i=c; i<c+n; ++i)
		dat[4096*(6143-r)+4095-i] = *pix++;
	break;
    case 3:	// 2-port: G,H
	if (r < 4096)
	    for (i=c; i<c+n; ++i)
		dat[4096*(4095-r)+4095-i] = *pix++;
	else
	    memcpy(dat+4096*(8191-r)+c, pix, 2*n);
	break;
    case 4:	// 2-port: E,H (NOT H,E!)
	if (r < 2048)
	    memcpy(dat+4096*r+c, pix, 2*n);
	else
	    memcpy(dat+4096*(6143-r)+c, pix, 2*n);
	break;
    case 5:	// 1-port: E
	memcpy(dat+4096*r+c, pix, 2*n);
	break;
    case 6:     // 1-port: F
	for (i=c; i<c+n; ++i)
	    dat[4096*r+4095-i] = *pix++;
	break;
    case 7:	// 1-port: G
	for (i=c; i<c+n; ++i)
	    dat[4096*(4095-r)+4095-i] = *pix++;

	break;
    case 8:	// 1-port: H
	memcpy(dat+4096*(4095-r)+c, pix, 2*n);
	break;
    }
}

//////////////////////////////////////////////
int imgdecode(unsigned short *impdu, IMG *img)
//////////////////////////////////////////////
{
    unsigned short old, final, diff, w, h, skp, tak;
    int i, N, K, bits2go, wordcnt, nzero;
    int headerr, overlfow;
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
    if (u != APID_HMI_SCIENCE_1 && u != APID_HMI_SCIENCE_2 &&
	u != APID_AIA_SCIENCE_1 && u != APID_AIA_SCIENCE_2) {
	return IMGDECODE_BAD_APID;
    }
    img->apid = u;
    img->telnum = impdu[11] >> 14;
    img->fsn = ((impdu[11] & 0x3fff) << 16) + impdu[12];
    img->fid = ((impdu[13] & 0xff) << 16) + impdu[14];
    img->cropid = impdu[15] >> 4;
    img->overflow = impdu[15] & 1;
    img->headerr = (impdu[15] >> 1) & 1;
    img->luid = impdu[17] >> 8;
    img->tap = impdu[16] >> 12;
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

    img->datavals = 0;
    img->npackets = 0;
    img->nerrors = 0;
    img->last_pix_err = 0;
    img->first_packet_time = UINT64_MAX;

    for (i = 0; i < MAXPIXELS; ++i)
	img->dat[i] = BLANK;

    for (i = 0; i < MAXHIST; ++i)
	img->hist[i] = 0;

    img->initialized = 1;

___DECODE_START___:

    //
    // read inverse lookup table
    //
    if (img->luid && !lut[img->luid]) {
	snprintf(fname, 64, TABLE_DIR "/lu/ilu%u", img->luid);
	fp = fopen(fname, "r");
	if (!fp) 
	    return IMGDECODE_NO_LOOKUP_TABLE;
	fscanf(fp, "%u", &u);
	if (u != img->luid)
	    return IMGDECODE_LOOKUP_ID_MISMATCH;
	lut[img->luid] = (unsigned short *) malloc(MAXHIST*2);
	if (!lut[img->luid])
	    return IMGDECODE_OUT_OF_MEMORY;
	for (i=0; i<MAXHIST; ++i)
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
	snprintf(fname, 64, TABLE_DIR "/crop/crop%u", img->cropid);
	fp = fopen(fname, "r");
	if (!fp)
	    return IMGDECODE_NO_CROP_TABLE;
	fscanf(fp, "%u", &u);
	if (u != img->cropid)
	    return IMGDECODE_CROP_ID_MISMATCH;
	fscanf(fp, "%hu %hu", &w, &h);
	if (!(w == 4096 && h == 4096) && !(w == 2048 && h == 8192))
	    return IMGDECODE_BAD_CROP_GEOMETRY;
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
	img->totalvals = MAXPIXELS;

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
    if (u < img->totalvals && u > img->totalvals - 3) {
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
	w = h = 4096;
	switch(img->tap) {
	case 0: 
	case 1: 
	case 3: 
	    w = 2048; h = 8192; 
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

#ifdef MAKE_IMGDECODE_MAIN

#include <fitsio.h>

void img2fits(IMG *i, STAT *s)
{
    int status = 0;
    char fn[64];
    long naxes[] = {4096,4096};
    long fpix[] = {1,1};
    short blank = BLANK;
    fitsfile *fp;
    int missvals;

    fn[0] = '!';
    snprintf(fn+1, 63, "%08u.fits[compress R]", i->fsn);
    fits_create_file(&fp, fn, &status);
    if (status) {
	fits_report_error(stderr, status);
	return;
    }
    fits_create_img(fp, 16, 2, naxes, &status);
    fits_write_pix(fp, TSHORT, fpix, MAXPIXELS, i->dat, &status);
    fits_write_key(fp, TSHORT, "BLANK", &blank, "", &status);
    fits_write_key(fp, TINT, "FSN", &i->fsn, "", &status);
    fits_write_key(fp, TINT, "DATAVALS", &i->datavals, "", &status);
    missvals = i->totalvals - i->datavals;
    fits_write_key(fp, TINT, "MISSVALS", &missvals, "", &status);
    fits_write_key(fp, TSHORT, "DATAMIN", &s->min, "", &status);
    fits_write_key(fp, TSHORT, "DATAMAX", &s->max, "", &status);
    fits_write_key(fp, TSHORT, "DATAMEDN", &s->median, "", &status);
    fits_write_key(fp, TDOUBLE, "DATAMEAN", &s->mean, "", &status);
    fits_write_key(fp, TDOUBLE, "DATA_RMS", &s->rms, "", &status);
    fits_write_key(fp, TDOUBLE, "DATASKEW", &s->skew, "", &status);
    fits_write_key(fp, TDOUBLE, "DATAKURT", &s->kurt, "", &status);
    fits_write_key(fp, TINT, "NUMPKTS", &i->npackets, "", &status);
    fits_write_key(fp, TINT, "NUMERRS", &i->nerrors, "", &status);
    fits_write_key(fp, TLOGICAL, "EOIERR", &i->last_pix_err, 
	    "decompression error at end of image?", &status);
    fits_write_chksum(fp, &status);
    fits_close_file(fp, &status);
    fits_report_error(stderr, status);
}

int main()
{
    IMG img;
    STAT imstat;
    unsigned char buf[1788];
    unsigned short *buf2 = (unsigned short *)buf;
    unsigned fsn, ofsn, apid;
    int npkt, ierr;
    
    img.initialized = 0;
    img.reopened = 0;
    ofsn = 0;
    npkt = 0;

    while(fread(buf,1,1788,stdin) == 1788) {
	apid = ((buf[18] << 8) + buf[19]) & 0x7ff;
	if (!(apid==400 || apid==410 || apid==500 || apid==510))
	    continue;
	fsn = (buf[32] << 24) + (buf[33] << 16) + (buf[34] << 8) + buf[35]; 
	if (fsn != ofsn && img.initialized) {
	    imgstat(&img, &imstat);
	    img2fits(&img, &imstat);
	    img.initialized = 0;
	    img.reopened = 0;
	}

	ierr = imgdecode(buf2+5,&img);
	if (ierr)
	    fprintf(stderr,"packet %d return code %d\n", npkt, ierr);
	++npkt;
	ofsn = fsn;
    }

    // write final image
    if (img.initialized) {
	imgstat(&img, &imstat);
	img2fits(&img, &imstat);
    }

    return 0;
}

#endif
