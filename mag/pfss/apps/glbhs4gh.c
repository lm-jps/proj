/*
 * Name:		glbhs4gh.c
 *
 * Description:		Special case of working part in jhelio2mlat.c and jqdotprod.c
 *
 * Function List:
 *			void helio2mlat(float *map, float *map_mlat, 
 *				int map_cols, int map_rows, int lmax, int map_lmax)
 *			void qdotprod(float *map_mlat, float *out_coeff, 
 *				int map_rows, int lmin, int lmax, float sinBdelta)
 *
 * Called by:		jsynop2gh.c
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Jul 25 2009
 *			v1.0a		Aug 06 2009
 *
 * Issues:
 *			v1.0a
 *			Added prototypes for Fortran fuctions (scopy, setplm, sgemm)
 *			
 */

#include <fftw3.h>
#include "astro.h"


#define absval(x)		(((x) < 0) ? (-(x)) : (x))
#define minval(x,y)		(((x) < (y)) ? (x) : (y))
#define maxval(x,y)		(((x) < (y)) ? (y) : (x))
#define very_small		(1.0e-30)
#define is_very_small(x)	(absval(x) < very_small)
#define is_even(x) 		(!((x)%2))
#define is_odd(x) 		((x)%2)


/* Simplified version of jhelio2mlat.c
 * No apodizing, weighted to 1, no mean substraction,
 * no normalize, assuming no bad points.
 */

void helio2mlat(float *map, float *map_mlat, int map_cols, int map_rows, int lmax, int map_lmax)
{
    int i, row, col;
    int map_cols2 = map_cols / 2;
    int nout = (lmax + 1) * 2;
    int lfft, nfft, nok;
    double norm = 1.0, normx;
    float *buf, *wbuf, *bp;
    float *ip, *inp = map;
    float *op, *outp = map_mlat;
    fftwf_plan fftwp;

    lfft = 2 * map_lmax;
    nfft = lfft + 2;

    // Get working buffer
    buf = (float *)(malloc(nfft * sizeof(float)));
    wbuf = (float *)(malloc(nfft * sizeof(float)));
    fftwp = fftwf_plan_r2r_1d(lfft, buf, wbuf, FFTW_R2HC, FFTW_ESTIMATE);

    // Do it
    for (row = 0; row < map_rows; row++) {
        // 0 at central meridian
        // First copy right side of meridian
        bp = buf;
        ip = inp + map_cols * row + map_cols2;
        for (col = 0; col <= map_cols2; col++) {
            if (!isnan(*ip)) {
                *bp++ = *ip;	// weight 1
            } else
                *bp++ = 0.0;
            ip++;
        }
        // Then do left side of meridian
        bp = buf + lfft - map_cols2;
        ip = inp + map_cols * row;
        for (col = 0; col < map_cols2; col++) {
            if (!isnan(*ip)) {
                *bp++ = *ip;
            } else
                *bp++ = 0.0;
            ip++;
        }

        // For missing (no zero_miss or normalize)
        norm = 1. / lfft;

        // Fourier transform
        fftwf_execute(fftwp);

        // Transpose, normalize
        // Real part
        for (int i = 0; i < nout / 2; i++) {
            op = outp + 2 * i * map_rows + row;
            *op = wbuf[i] * norm;
        }
        // Imaginary part, m = 0 is 0
        *(outp + map_rows + row) = 0.0;
        // Use normx to get complex conjugate
        normx = -norm;
        for (int i = 0; i < nout / 2; i++) {
            op = outp + (2 * i + 1) * map_rows + row;
            *op = wbuf[lfft - i] * normx;
        }          
    } // End of row
}







/* Prototypes of some FORTRAN function */

void scopy_(int *n, float *x, const int *incx, float *y, const int *incy);
extern int setplm_(int *lmin, int *lmax, int *m, int *nx, int *indx, 
        double *x, int *nplm, double *plm);
void sgemm_(const char *transa, const char *transb,
        int *l, int *n, int *m, float *alpha,
        const float *a, int *lda, float *b, int *ldb,
        float *beta, void *c, int *ldc, int na, int nb);



/* Simplified version of jqdotprod.c
 * No bad images, single point in time series
 * timechunk is 1, lchunk is 1
 * Need to include <astro.h> <math.h> in main
 */

void qdotprod(float *map_mlat, float *out_coeff, int map_rows, int lmin, int lmax, float sinBdelta)
{
    int i, l, m, modd, meven;

    float *oddpart, *evenpart, *inptr, *workptr, *mptr, *outptr;
    float *folded, *masks, *real4evenl, *real4oddl, *imag4evenl, *imag4oddl, *outx;
    /* used for setting up plm's */
    double *plm, *saveplm, *latgrid;
    int *indx;
    double *plmptr;
    float *maskptr;

    int nrecs = 1;	// Single point time series
    int lmax1 = lmax + 1;
    int mx, msize = lmax1, foldedsize;
    int imagesize = map_rows * 2 * msize;	// Size of map_mlat

    /* for scopy call */
    int increment1 = 1, increment2 = 2;

    /* arguments for sgemm call */
    char transpose[] = "t";
    char normal[] = "n";
    float one = 1.0;
    float zero = 0.0;
    float cnorm; /* Constant to get proper normalization. */

    int lfirst, llast, ifirst, ilast, lstart, ldone;
    int lfirsteven, llasteven, nevenl, lfirstodd, llastodd, noddl;
    int fromoffset, tooffset, imageoffset;

    int nlat = map_rows / 2, latx, ilatx, poslatx, neglatx, moffset;
    /* SGI's like odd leading dimensions of the first array in sgemm */
    int nlatx = 2 * (nlat / 2) + 1;
    /* make nlatx divisible by 4 on linux systems */
    #ifdef __linux__
        if (nlat % 4) nlatx = 4 * (nlat / 4 + 1);
        else nlatx = nlat;
    #endif

    int snx;
    int maxnsn = nrecs;	// Equals to 1
    int nsn = maxnsn, fournsn = 4 * nsn;
    int lchunksize = msize;	// lmax + 1

    real4evenl = (float *)(malloc(nlatx * maxnsn * sizeof(float)));
    real4oddl = (float *)(malloc(nlatx * maxnsn * sizeof(float)));
    imag4evenl = (float *)(malloc(nlatx * maxnsn * sizeof(float)));
    imag4oddl = (float *)(malloc(nlatx * maxnsn * sizeof(float)));
    outx = (float *)(malloc(maxnsn * 2 * lchunksize * sizeof(float)));

    plm = (double *)(malloc(lmax1 * nlat * sizeof(double))); 
    saveplm = (double *)(malloc(lmax1 * nlat * 2 * sizeof(double)));
    latgrid = (double *)(malloc(nlat * sizeof(double)));
    for (i = 0; i < nlat; i++) latgrid[i] = (i + 0.5) * sinBdelta;

    indx = (int *)(malloc(lmax1 * sizeof(int)));
    for (l = 0; l <= lmax; l++) indx[l] = l;

    masks = (float *)(malloc(nlat * lchunksize * sizeof(float)));

    foldedsize = 4 * nlat * lmax1 * maxnsn;
    folded = (float *)(malloc(foldedsize * sizeof(float)));

    oddpart = (float *)(malloc(nlat * sizeof(float)));
    evenpart = (float *)(malloc(nlat * sizeof(float)));

    inptr = map_mlat;
    imageoffset = 0;		// irec = 0
    // For each m, re and im
    for (mx = 0; mx < 2 * msize; mx++) {
        moffset = mx * map_rows;
        mptr = inptr + moffset;
        for (latx = 0; latx < nlat; latx++) {
            poslatx = nlat + latx; neglatx = nlat - 1 - latx;
            evenpart[latx] = mptr[poslatx] + mptr[neglatx];
            oddpart[latx] = mptr[poslatx] - mptr[neglatx]; 
        }
        workptr = folded + imageoffset + moffset;
        scopy_ (&nlat, evenpart, &increment1, workptr, &increment1);
        workptr += nlat;
        scopy_ (&nlat, oddpart, &increment1, workptr, &increment1);
    }

    /* We now have folded data for a chunk of sn's */
    /* Now do Jesper's tricks */
    /* ldone is the last l for which plm's have been set up */
    ldone = -1;

    /* loop on each chunk of l's */
    //lchunkfirst = lmin / lchunksize;	// 0
    //lchunklast = lmax / lchunksize;	// 0

    /* assumes the output datasets are in order with time increasing most rapidly */

    cnorm = sqrt(2.) * sinBdelta;	// norm = 1

    lfirst = 0;	// lchunk = 0
    llast = lmax;

    /* get the first and last indeces into the l-m array */
    ifirst = lfirst * (lfirst + 1) / 2;		// 0
    ilast = (llast + 1) * (llast + 2) / 2 - 1;	// (lmax+1)(lmax+2)/2-1

    outptr = out_coeff;	// (lmax+1)(lmax+2)

    /* loop on each m */
    for (m = 0; m <= llast; m++) {
        modd = is_odd(m);
        meven = !modd;
        lstart = maxval(lfirst, m); /* no l can be smaller than this m */

        /* set up masks (plm's) for this m and chunk in l */
        if ((lstart - 1) == ldone) {
            /* get saved plms if any */
            if ((lstart - 2) >= m)
                for (latx = 0; latx < nlat; latx++)
                    plm[(lstart - 2) * nlat + latx] = saveplm[m * nlat + latx];
            if ((lstart - 1) >= m)
                for (latx = 0; latx < nlat; latx++)
                    plm[(lstart - 1) * nlat + latx] = saveplm[msize * nlat + m * nlat + latx];
            /* then set up the current chunk */
            setplm_ (&lstart, &llast, &m, &nlat, indx, latgrid, &nlat, plm); 
        } else {
            /* This fixes the lmin != 0 problem */
            setplm_ (&m, &llast, &m, &nlat, indx, latgrid, &nlat, plm); 
        }

        /* save plm's for next chunk */
        if ((llast-1) >= m)
            for (latx = 0; latx < nlat; latx++)
               saveplm[m * nlat + latx] = plm[(llast - 1) * nlat + latx];
        for (latx = 0; latx < nlat; latx++)
            saveplm[msize * nlat + m * nlat + latx] = plm[llast * nlat + latx];
        ldone = llast;

        /* copy plm's into masks */
        /* note that this converts from double to single precision */
        /* the test prevents underflows which gobble CPU time */

        plmptr = plm + nlat * lstart;
        maskptr = masks;
        latx = nlat * (llast - lstart + 1);
        for (ilatx = 0; ilatx < latx; ilatx++)
            maskptr[ilatx] = plmptr[ilatx];

        snx = 0;	// Single chunk
        /* select folded data for real/imag l's and this m 
           into temporay arrays for matrix multiply */
        /* TO DO - pull offset calculations out of loop */
        /* New code with odd leading dimension */
        scopy_ (&nlat,
                folded + nlat * (4 * m + modd) + snx * imagesize,
                &increment1,
                real4evenl + snx * nlatx,
                &increment1);
        scopy_ (&nlat,
                folded + nlat * (4 * m + meven) + snx * imagesize,
                &increment1,
                real4oddl + snx * nlatx,
                &increment1);
        scopy_ (&nlat,
                folded + nlat * (4 * m + 2 + modd) + snx * imagesize,
                &increment1,
                imag4evenl + snx * nlatx,
                &increment1);
        scopy_ (&nlat,
                folded + nlat * (4 * m + 2 + meven) + snx * imagesize,
                &increment1,
                imag4oddl + snx * nlatx,
                &increment1);

        /* do even l's */
        lfirsteven = is_even(lstart) ? lstart : lstart + 1;
        llasteven = is_even(llast) ? llast : llast - 1;
        nevenl = (llasteven - lfirsteven) / 2 + 1; /* number of even l's */
        /* do real part */
        /* All parts used to have alpha=&one, now have alpha=&cnorm */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &nevenl,   /* number of even l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,    /* scalar multiplier of op(A) */
                real4evenl,  /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat * (lfirsteven - lstart), /* matrix B */
                &map_rows,  /* 2*nlat, use every other row (nlat long) of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn * 2 * (lfirsteven - lstart), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */
        /* do imag part */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &nevenl,   /* number of even l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,    /* scalar multiplier of op(A) */
                imag4evenl,  /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat * (lfirsteven - lstart), /* matrix B */
                &map_rows,  /* 2*nlat, use every other nlat-long row of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn * ( 2 * (lfirsteven - lstart) + 1), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */

        /* do odd l's */
        lfirstodd = is_odd(lstart) ? lstart : lstart + 1;
        llastodd = is_odd(llast) ? llast : llast - 1; 
        noddl = (llastodd - lfirstodd) / 2 + 1; /* number of odd l's */
        /* do real part */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &noddl,    /* number of odd l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,   /* scalar multiplier of op(A) */
                real4oddl,   /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat * (lfirstodd - lstart), /* matrix B */
                &map_rows,  /* 2*nlat, use every other nlat-long row of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn * 2 * (lfirstodd - lstart), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */
        /* do imag part */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &noddl,    /* number of odd l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,    /* scalar multiplier of op(A) */
                imag4oddl,  /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat * (lfirstodd - lstart), /* matrix B */
                &map_rows,  /* 2*nlat, use every other nlat-long row of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn * (2 * (lfirstodd - lstart) + 1), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */

        /* copy outx into out sds */
        /* alternate real and imaginary values in out - as in pipeLNU */
        for (l = lstart; l <= llast; l++) {	// 0
             fromoffset = 2 * nsn * (l - lstart);		// 0
             tooffset = 2 * nsn * (l * (l + 1) / 2 + m - ifirst);
             scopy_ (&nsn,
                     outx + fromoffset,
                     &increment1,
                     outptr + tooffset,
                     &increment2);
             scopy_ (&nsn,
                     outx + fromoffset + nsn,
                     &increment1,
                     outptr + tooffset + 1,
                     &increment2);
        } /* end loop through l's for this m */

    } /* end loop on m */
}
