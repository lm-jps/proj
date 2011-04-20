/*
 * Module name:		synop2gh.c
 *
 * Description:		Compute harmonic expansion coefficient for given synoptic map
 *			Two modes are provided: (0) Numerical integration (1) FFT
 *			For FFT method, look in jhelio2mlat.c and jqdotprod.c
 *			Description of Numerical integration method and normalization
 *			problem of FFT can be found in pfssnotes.pdf
 *
 * Calling:		glbhs4gh.c, pfss_pkg.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *			jv2ts.c module in JSOC pipeline
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Jul 25 2009
 *			v1.0a		Aug 06 2009
 *			v1.0b		Sep 17 2009
 *			v1.1		Oct 12 2009
 *			v1.2		Mar 31 2009
 *			v1.3		Apr 02 2009
 *          v1.4        Apr 05 2009
 *
 * Issues:
 *			v1.0
 *			Input synoptic must be in sinlat! (How to check keywords for it?)
 *			Set CAR_ROT and LMAX as prime key for output series
 *			Output arrays have variable dimensions
 *			Need to confirm keywords to accommodate different synoptic maps
 *			v1.0a
 *			Added some keyword check (sinlat: CTYPE2, lon0, etc.)
 *			Added computation on latitude maps
 *			Calling pfss.h instead
 *			v1.0b
 *			Set default method to mag (FFT=0) as FFT method takes data at
 *			0, dphi, 2dphi... rather than 0.5dphi, 1.5dphi, 2.5dphi
 *			v1.1
 *			For mag method, shift the coordinate left by half width of the orignal
 *			resolution (so the first pixel is centered at 0).
 *			Added a grid variable to send the grid coords into gh_pfss()
 *			v1.2
 *			Moved keywords checking into the loop
 *			v1.3
 *			Added check for cell center and fft
 *			Added link to synop in compliance with su_xudong.hmi_gh
 *          v1.4
 *          Splitted g and h
 *
 * Example:
 *			jsynop2gh "in=su_xudong.hmi_synop_polfil[1928]" out=su_xudong.hmi_gh
 *			jsynop2gh "in=su_xudong.wsosynop[1988]" out=su_xudong.hmi_gh DPH0=0 LMAX=9 FFT=0
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define	DTOR	(M_PI / 180.)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#include "pfss.h"
#include "pfss_pkg.h"
#include "glbhs4gh.c"
#include "timing.c"

char *module_name = "jsynop2gh";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_INT, "LMAX", "200", "Maximum l number desired."},
    {ARG_INT, "FFT", "1", "Set to one for FFT, numerical integration by default"},
    {ARG_FLOAT, "DPH0", "0.1", "Original longitude resolution."},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=all messages"},
    {ARG_END}
};


/* ############## Fixing Normalization ############## */

/* This function is called after FFT, in order to convert the
 * complex y_l^m into real gh that fits the magnetic convention.
 * We use (ee pfssnotes.pdf for details):
 * g_lm = sqrt[(2l+1)(2-d_0)] / 2 * Re(y_l^m)
 * h_lm = sqrt[(2l+1)(2-d_0)] / 2 * Im(y_l^m)
 * where d_0 = 1 (m=0), 0 (m!=0) 
 */

void ylm2gh(ylm, g, h, lmax)
    float *ylm, *g, *h;
    int lmax;
{
    int l, m;
    int lmax1 = lmax + 1;
    float *gl, *hl, *yl, *ym;
    float kl[lmax1], km[lmax1];

    for (l = 0; l <= lmax; l++) kl[l] = sqrt(2.0 * l + 1.0);
    km[0] = 0.5;
    for (m = 1; m <= lmax; m++) km[m] = sqrt(0.5);

    for (l = 0; l <= lmax; l++) {
        gl = g + l * lmax1;
        hl = h + l * lmax1;
        yl = ylm + l * (l + 1);
        for (m = 0; m <= l; m++) {
            ym = yl + m * 2;
            gl[m] = ym[0] * kl[l] * km[m];
            hl[m] = ym[1] * kl[l] * km[m];
        }
        hl[0] = 0.0;
    }
}


/* ################## Main Module ################## */

int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS, *outRS;
    int irec, nrecs;
    DRMS_Record_t *inRec, *outRec;
    DRMS_Segment_t *inSeg, *outSeg_g, *outSeg_h;
    DRMS_Array_t *inArray, *outArray_g, *outArray_h;
    DRMS_Link_t *srcLink;
    float *inData, *g, *h;

    char *inchecklist[] = {"CAR_ROT", "CDELT1", "CDELT2", "CTYPE2", 
                           "LON_FRST"};		// For now
    char *outchecklist[] = {"CAR_ROT", "LON0", "LMAX"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int verbflag;
    int outDims[2];
    int lmax, lmax1, lmax2;
    int fft, sinlat;
    int map_lmax;
    int car_rot;
    float lon_frst, lon0, dph0;
    double cdelt1, cdelt2;
    float sinBdelta;
    char *ctype2;

    int i, l, m, itest;
    int np, nt, npnt;

    struct Grid *grid;

    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);

    /* Temporary arrays for FFT */
    float *map_mlat, *ylm;

    /* Get parameters */
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    lmax = params_get_int(&cmdparams, "LMAX");
    fft = params_get_int(&cmdparams, "FFT");
    dph0 = params_get_float(&cmdparams, "DPH0");
    verbflag = params_get_int(&cmdparams, "VERB");

    lmax1 = lmax + 1; lmax2 = lmax + 2;
    outDims[0] = outDims[1] = lmax1;	// Output array dimensions

    /* Open input */
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;

    /* Check for essential input keywords */
    inRec_tmp = inRS->records[0];
    for (itest = 0; itest < ARRLENGTH(inchecklist); itest++) {
        inkeytest = hcon_lookup_lower(&inRec_tmp->keywords, inchecklist[itest]);
        if (!inkeytest) {
            fprintf(stderr, "ERROR: required input keyword %s is missing\n", inchecklist[itest]);
            return 0;
        }
    }

    /* Check for essential output keywords */
    outRec_tmp = drms_create_record(drms_env, outQuery, DRMS_TRANSIENT, &status);
    if (status) DIE("Unable to open record in output series");
    for (itest = 0; itest < ARRLENGTH(outchecklist); itest++) {
        outkeytest = hcon_lookup_lower(&outRec_tmp->keywords, outchecklist[itest]);
        if (!outkeytest || outkeytest->info->islink || outkeytest->info->recscope == 1) {
            fprintf(stderr, "Output keyword %s is missing/constant/a link.\n", outchecklist[itest]);
            return 0;
        }
    }
    drms_close_record(outRec_tmp, DRMS_FREE_RECORD);

    /* Create output */
    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output recordset not created");

    /* Do this for each record */
    for (irec = 0; irec < nrecs; irec++)
    {
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }

        /* Input record and data */
        inRec = inRS->records[irec];
        inSeg = drms_segment_lookupnum(inRec, 0);	// Single segment
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No data file found. \n");
            drms_free_array(inArray);
            continue;
        }
        inData = (float *)inArray->data;
        np = inArray->axis[0]; nt = inArray->axis[1];

        /* Checking keywords */
        ctype2 = drms_getkey_string(inRec_tmp, "CTYPE2", &status);
        sinlat = (strcmp(ctype2, "Sine Latitude") == 0 || strcmp(ctype2, "CRLT-CEA") == 0) ? 1 : 0;
        if (!sinlat && fft)
            DIE("Input not in sine latitude, can't use FFT");
        cdelt1 = fabs(drms_getkey_double(inRec_tmp, "CDELT1", &status));
        cdelt2 = fabs(drms_getkey_double(inRec_tmp, "CDELT2", &status));
        map_lmax = (int)(180. / cdelt1);
        if (lmax > map_lmax) DIE("Request lmax greater than input lmax");
        sinBdelta = cdelt2;
        // FFT requires grid at 0, dph, 2dph... Rebinning shifts cell center
        if ((fabs(cdelt1 - dph0) >= 0.1 * dph0) && fft)
        	DIE("Not original data resolution, can't use FFT");
        
        /* Make sure there's no missing points */
        npnt = np * nt;
        for (i = 0; i < npnt; i++) if (isnan(inData[i])) break;
        if (i < npnt && isnan(inData[i])) {
            printf("Missing points in data. \n");
            drms_free_array(inArray);
            continue;
        }

        /* Output data */
        g = (float *)calloc(lmax1 * lmax1, sizeof(float));
        h = (float *)calloc(lmax1 * lmax1, sizeof(float));
        outArray_g = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, g, &status);
        outArray_h = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, h, &status);

        /* This is the working part. Two modes for choice. */
        if (fft) {

            if (verbflag) SHOW("starting FFT method... ");

            printf("%d, %d\n", nt, lmax1);

            /* FFT, in 3 steps, see glbhs4gh.c */
            /* Step 1: Adapted from jhelio2mlat.c */
            map_mlat = (float *)(malloc(nt * 2 * lmax1 * sizeof(float)));
SHOW("0");
            helio2mlat(inData, map_mlat, np, nt, lmax, map_lmax);
SHOW("1");
            /* Step 2: Adapted from jqdotprod.c */
            ylm = (float *)(malloc(lmax1 * lmax2 * sizeof(float)));
            qdotprod(map_mlat, ylm, nt, 0, lmax, sinBdelta);
SHOW("2");
            /* Step 3: Re-normalization for mag convention */
            ylm2gh(ylm, g, h, lmax);
            free(map_mlat); free(ylm);
            if (verbflag) SHOW("done.\n");
            
        } else {
        
            /* Specify grid */
            grid = (struct Grid *)(malloc(sizeof(struct Grid)));
            grid->np = np; grid->nt = nt;
            grid->ph = (float *)(calloc(np, sizeof(float)));
            grid->th = (float *)(calloc(nt, sizeof(float)));
            for (int i = 0; i < np; i++)
                grid->ph[i] = (i + 0.5) / np * 2 * M_PI - dph0 / 2.0 * DTOR;
            for (int j = 0; j < nt; j++)
                grid->th[j] = sinlat ? acos(1.0 - 2.0 * (nt - j - 0.5) / nt) : 
                    ((nt - j - 0.5) / nt * M_PI);
            /* Traditional numerical integration, see pfss_pkg.c */
            if (verbflag) SHOW("starting Mag method... ");
            gh_pfss(inData, g, h, grid, lmax, sinlat);
            if (verbflag) SHOW("done.\n");
            free(grid);
        
        }

        // Filling in nan
        for (l = 0; l <= lmax; l++) {
            for (m = l + 1; m <= lmax; m++) {
                g[l * lmax1 + m] = DRMS_MISSING_FLOAT;
                h[l * lmax1 + m] = DRMS_MISSING_FLOAT;
            }
        }

        drms_free_array(inArray);

        /* Outpur record */
        outRec = outRS->records[irec];
        outSeg_g = drms_segment_lookup(outRec, "g");
        outSeg_h = drms_segment_lookup(outRec, "h");
        for (i = 0; i < 2; i++) {
            outSeg_g->axis[i] = outArray_g->axis[i];	// For variable dimensions
            outSeg_h->axis[i] = outArray_h->axis[i];
        }
        outArray_g->parent_segment = outSeg_g;
        outArray_h->parent_segment = outSeg_h;

        /* Result writing */
        status = drms_segment_write(outSeg_g, outArray_g, 0);
        if (status) DIE("Problem writing g file");
        //drms_free_array(outArray_g);
        status = drms_segment_write(outSeg_h, outArray_h, 0);
        if (status) DIE("Problem writing h file");
        //drms_free_array(outArray_h);
        free(g); free(h);

        /* Set keywords */
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        drms_copykey(outRec, inRec, "T_REC");
        drms_copykey(outRec, inRec, "T_OBS");
        car_rot = drms_getkey_int(inRec, "CAR_ROT", NULL);
        lon_frst = drms_getkey_float(inRec, "LON_FRST", NULL);
        lon0 = car_rot * 360. - (lon_frst + 360.);	// This is for now
        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_setkey_float(outRec, "LON0", lon0);
        drms_setkey_int(outRec, "LMAX", lmax);
        drms_setkey_int(outRec, "FFT", fft);
        drms_setkey_int(outRec, "SINLAT", sinlat);
        drms_setkey_int(outRec, "NP", np);
        drms_setkey_int(outRec, "NT", nt);

        /* Set links */
        srcLink = hcon_lookup_lower(&outRec->links, "SYNOP");
        if (srcLink) {
            drms_setlink_static(outRec, "SYNOP", inRec->recnum);
        }

        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                     irec, wt - wt1, ct - ct1);
        }
    }

    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);

    if (verbflag) {
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
                wt - wt0, ct - ct0);
    }

    return(DRMS_SUCCESS);
}
