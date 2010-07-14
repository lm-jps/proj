/*
 * Module name:		jgh2pfss.c
 *
 * Description:		Compute PFSS field on a give grid. Read in harmonic coefficients
 *			first. Refer to pfssnotes.pdf for details.
 *
 * Calling:		pfss_pkg.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v0.0 (original)	Jun 14 2000
 *			v1.0		Jul 25 2009
 *			v1.0a		Aug 06 2009
 *			v1.0b		Nov 10 2009
 *			v1.1		Apr 03 2010
 *          v1.2        Apr 06 2010
 *
 * Issues:
 *			v1.0
 *			Seriously, keywords? Output specifying SINLAT
 *			Maybe incorporating irregular grids from data series?
 *			v1.0a
 *			Added several keywords
 *			Calling pfss.h instead
 *			Need to develop an external data series for grid spec
 *			v1.0b
 *			free(grid) moved out of the main loop, freed components of grid
 *			v1.1
 *			Changes to comply with su_xudong.hmi_gh and su_xudong.hmi_pfss
 *          v1.2
 *          Split g and h
 *
 * Example:
 *			jgh2pfss "in=su_xudong.hmi_gh[1928][0]" out=su_xudong.hmi_pfss
 *			jgh2pfss "in=su_xudong.hmi_gh[1988][0]" out=su_xudong.hmi_pfss LMAX=9 NP=72 NT=30 NR=2 SINLAT=1
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
#include "timing.c"

char *module_name = "jgh2pfss";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_INT, "LMAX", "72", "Maximum l number desired."},
    {ARG_INT, "SINLAT", "0", "Set to zero for latitude grid; lat by default."},
    {ARG_INT, "NP", "360", "Number of Grid points in longitude."},
    {ARG_INT, "NT", "180", "Number of Grid points in latitude."},
    {ARG_INT, "NR", "16", "Number of Grid points in radial direction."},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=all messages"},
    {ARG_END}
};


/* ################## Main Module ################## */

int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS, *outRS;
    int irec, nrecs;
    DRMS_Record_t *inRec, *outRec, *linkedRec;
    DRMS_Segment_t *inSeg_g, *inSeg_h;
    DRMS_Segment_t *outSegBr, *outSegBt, *outSegBp;
    DRMS_Array_t *inArray_g, *inArray_h;
    DRMS_Array_t *outArrayBr, *outArrayBt, *outArrayBp;
    DRMS_Link_t *synopLink, *ghLink;
    float *inData_g, *inData_h;
    float *g, *h;
    float *Br, *Bt, *Bp;

    char *inchecklist[] = {"CAR_ROT", "LON0", "LMAX"};
    char *outchecklist[] = {"CAR_ROT", "LON0", "LMAX", "CTYPE2"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int verbflag;
    int outDims[3];
    int lmax, lmax1;
    int want_lmax, have_lmax;
    int car_rot, sinlat;
    int np, nt, nr;

    int i, j, n;
    int l, m, itest;
    int n0, n1, n2;

    struct Grid *grid;

    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);

    /* Get parameters */
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    want_lmax = params_get_int(&cmdparams, "LMAX");
    sinlat = params_get_int(&cmdparams, "SINLAT");
    np = params_get_int(&cmdparams, "NP");
    nt = params_get_int(&cmdparams, "NT");
    nr = params_get_int(&cmdparams, "NR");
    verbflag = params_get_int(&cmdparams, "VERB");

    outDims[0] = np; outDims[1] = nt; outDims[2] = nr;

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
    have_lmax = drms_getkey_int(inRec_tmp, "LMAX", &status);
    if (want_lmax > have_lmax) DIE("Request lmax greater than input lmax");

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

    /* Specify grid */
    grid = (struct Grid *)(malloc(sizeof(struct Grid)));
    grid->np = np; grid->nt = nt; grid->nr = nr;
    grid->ph = (float *)(calloc(np, sizeof(float)));
    grid->th = (float *)(calloc(nt, sizeof(float)));
    grid->rr = (float *)(calloc(nr, sizeof(float)));
    for (i = 0; i < np; i++)
        grid->ph[i] = (i + 0.5) / np * 2 * M_PI;
    for (j = 0; j < nt; j++)
        grid->th[j] = sinlat ? acos(1.0 - 2.0 * (nt - j - 0.5) / nt) : ((nt - j - 0.5) / nt * M_PI);
    if (nr != 1) {
        for (n = 0; n < nr; n++)
            grid->rr[n] = (RS - R0) / (nr - 1) * n + R0;
    } else {
        grid->rr[0] = RS;
        printf("Only one point specified in radial direction, compute source surface.\n");
    }
    
    /*
    for (i = 0; i < np; i++) printf("%8.1f", grid->ph[i] / DTOR);
    printf("\n");
    for (j = 0; j < nt; j++) printf("%8.1f", cos(grid->th[j]) * 15.);
    printf("\n");
    for (n = 0; n < nr; n++) printf("%8.2f", grid->rr[n]);
    printf("\n");
    */

    // Do this for each record set
    for (irec = 0; irec < nrecs; irec++) {

        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }
        
        /* Input record and data */
        inRec = inRS->records[irec];
        inSeg_g = drms_segment_lookup(inRec, "g");
        inArray_g = drms_segment_read(inSeg_g, DRMS_TYPE_FLOAT, &status);
        if (status) DIE("No g data file found. \n");
        inSeg_h = drms_segment_lookup(inRec, "h");
        inArray_h = drms_segment_read(inSeg_h, DRMS_TYPE_FLOAT, &status);
        if (status) DIE("No h data file found. \n");
        inData_g = (float *)inArray_g->data; inData_h = (float *)inArray_h->data;

        // Check for lmax, truncate g, h
        car_rot = drms_getkey_int(inRec, "CAR_ROT", &status);
        n0 = inArray_g->axis[0];		// have_lmax + 1
        n1 = inArray_g->axis[1];		// have_lmax + 1
        if (inArray_g->axis[0] != inArray_g->axis[1] || 
            inArray_h->axis[0] != inArray_h->axis[1] ||
            (inArray_g->axis[1] - 1) != have_lmax || 
            (inArray_h->axis[1] - 1) != have_lmax)
            DIE("Wrong input data dimension");
        if (want_lmax > have_lmax)
            DIE("Request lmax greater than available lmax");
        lmax = want_lmax;
        lmax1 = lmax + 1;
        g = (float *)(calloc(lmax1 * lmax1, sizeof(float)));
        h = (float *)(calloc(lmax1 * lmax1, sizeof(float)));
        for (int l = 0; l <= lmax; l++) {
            for (int m = 0; m <= l; m++) {
                g[l * lmax1 + m] = inData_g[l * n0 + m];
                h[l * lmax1 + m] = inData_h[l * n0 + m];
            }
        }
        g[0] = 0.0; h[0] = 0.0;	// No monopole
        drms_free_array(inArray_g); drms_free_array(inArray_h);

/*        
        for (l = 0; l <= lmax; l++) {
            for (m = 0; m <= l; m++) {
                printf("%10.5f", g[l * lmax1 + m]);
            }
            printf("\n");
        }      
        for (l = 0; l <= lmax; l++) {
            for (m = 0; m <= l; m++) {
                printf("%10.5f", h[l * lmax1 + m]);
            }
            printf("\n");
        } 
*/

        /* Output data */
        outArrayBr = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayBt = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayBp = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        Br = (float *)outArrayBr->data;
        Bt = (float *)outArrayBt->data;
        Bp = (float *)outArrayBp->data;

        /* This is the working part */
        if (verbflag) SHOW("start pfss computation... ");
        pfss(g, h, grid, Br, Bt, Bp, lmax, 0.0);
        if (verbflag) SHOW("done.\n");
        free(g); free(h);

        /* Output record (3 segs) */
        outRec = outRS->records[irec];
        outSegBr = drms_segment_lookupnum(outRec, 0);
        outSegBt = drms_segment_lookupnum(outRec, 1);
        outSegBp = drms_segment_lookupnum(outRec, 2);
        for (i = 0; i < 3; i++) {
            outSegBr->axis[i] = outArrayBr->axis[i];
            outSegBt->axis[i] = outArrayBt->axis[i];
            outSegBp->axis[i] = outArrayBp->axis[i];
        }
        outArrayBr->parent_segment = outSegBr;
        outArrayBt->parent_segment = outSegBt;
        outArrayBp->parent_segment = outSegBp;

        /* Result writing */
        status = drms_segment_write(outSegBr, outArrayBr, 0);
        if (status) DIE("problem writing Br");
        drms_free_array(outArrayBr);
        status = drms_segment_write(outSegBt, outArrayBt, 0);
        if (status) DIE("problem writing Bt");
        drms_free_array(outArrayBt);
        status = drms_segment_write(outSegBp, outArrayBp, 0);
        if (status) DIE("problem writing Bp");
        drms_free_array(outArrayBp);

        /* Set keywords */
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_copykey(outRec, inRec, "LON0");
        drms_setkey_int(outRec, "LMAX", lmax);
        drms_setkey_float(outRec, "RSS", RS);
        drms_setkey_float(outRec, "CDELT1", - 360. / np);
        drms_setkey_float(outRec, "CRPIX1", 1.0);
        drms_setkey_float(outRec, "CRVAL1", car_rot  * 360. - 180. / np);
        drms_setkey_float(outRec, "CRPIX2", nt / 2.0 + 0.5);
        drms_setkey_float(outRec, "CRVAL2", 0.0);
        if (sinlat) {
            drms_setkey_string(outRec, "CTYPE2", "CRLT-CEA");
            drms_setkey_float(outRec, "CDELT2", 2. / nt);
            drms_setkey_string(outRec, "CUNIT2", "Sine latitude");
        } else {
            drms_setkey_string(outRec, "CTYPE2", "CRLT-CAR");
            drms_setkey_float(outRec, "CDELT2", 180. / nt);
            drms_setkey_string(outRec, "CUNIT2", "degree");
        }
        drms_setkey_float(outRec, "CDELT3", (RS - R0) / (nr - 1));
        drms_setkey_float(outRec, "CRPIX3", 1.0);
        drms_setkey_float(outRec, "CRVAL3", 1.0);
        
        /* Set links */
        ghLink = hcon_lookup_lower(&outRec->links, "HARMONICS");
        if (ghLink) {
            drms_setlink_static(outRec, "HARMONICS", inRec->recnum);
            linkedRec = drms_link_follow(outRec, "HARMONICS", &status);
        }
        if (!status) {
	        synopLink = hcon_lookup_lower(&linkedRec->links, "SYNOP");
    	    if (synopLink) {
            	drms_setlink_static(outRec, "SYNOP", linkedRec->recnum);
        	}
        }

        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                     irec, wt - wt1, ct - ct1);
        }
    }

    free(grid->rr); free(grid->th); free(grid->ph); free(grid);

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
