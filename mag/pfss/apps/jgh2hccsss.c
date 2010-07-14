/*
 * Module name:		jgh2hccsss.c
 *
 * Description:		Compute HCCSSS field on a give grid. Read in harmonic coefficients
 *			below cusp surface. Calculate coefficient above the cusp surface, then
 *			field vectors. Refer to pfssnotes.pdf and Zhao et al (1995) for details.
 *
 * Calling:		pfss_pkg.c
 *
 * Original Source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v0.0 (orininal)
 *			v1.0		Oct 20 2009
 *			v1.0a		Nov 10 2009
 *			v1.1		Apr 13 2009
 *
 * Issues:
 *			v1.0
 *			Need to develop a data series for external grid spec
 *			v1.0a
 *			free(grid) moved out of the main loop, freed components of grid
 *			v1.1
 *			adapted for new jsds, new headers, new keywords
 *
 * Example:		jgh2hccsss "in=su_xudong.hmi_gh[1928][0]" out=su_xudong.hmi_hccsss
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
#include "hccsss_pkg.h"

#include "timing.c"

char *module_name = "jgh2hccsss";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_INT, "LMAX", "72", "Maximum l number used below rcp."},
    {ARG_INT, "LC", "15", "Maximum l number desired above rcp."},
    {ARG_INT, "SINLAT", "0", "Set to zero for latitude grid; lat by default."},
    {ARG_INT, "NP", "360", "Number of Grid points in longitude."},
    {ARG_INT, "NT", "180", "Number of Grid points in latitude."},
    {ARG_INT, "NR", "31", "Number of Grid points in radial direction."},
    {ARG_FLOAT, "DR", "0.1", "Step size in radial direction."},
    {ARG_FLOAT, "APAR", "0.2", "Parameterized length scale of horizontal current"},
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
    DRMS_Segment_t *inSegG, *inSegH;
    DRMS_Segment_t *outSegBr, *outSegBt, *outSegBp, *outSegGc, *outSegHc;
    DRMS_Array_t *inArrayG, *inArrayH;
    DRMS_Array_t *outArrayBr, *outArrayBt, *outArrayBp, *outArrayGc, *outArrayHc;
    DRMS_Link_t *synopLink, *ghLink;
    float *inDataG, *inDataH;
    float *g, *h, *gc, *hc;
    float *Br, *Bt, *Bp;

    char *inchecklist[] = {"CAR_ROT", "LON0", "LMAX"};
    char *outchecklist[] = {"CAR_ROT", "LON0", "LMAX", "LC", "CTYPE2"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int outDims[3], outDimsGHc[2];
    int lmax, lmax1, lc, lc1;
    int want_lmax, have_lmax;
    int sinlat, car_rot, verbflag;
    int np, nt, nr;
    float apar, dr;

    int i, j, n;
    int l, m, itest;
    int n0, n1;

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
    lc = params_get_int(&cmdparams, "LC");
    sinlat = params_get_int(&cmdparams, "SINLAT");
    np = params_get_int(&cmdparams, "NP");
    nt = params_get_int(&cmdparams, "NT");
    nr = params_get_int(&cmdparams, "NR");
	dr = params_get_float(&cmdparams, "DR");
    apar = params_get_float(&cmdparams, "APAR");
    verbflag = params_get_int(&cmdparams, "VERB");

    outDims[0] = np; outDims[1] = nt; outDims[2] = nr;
    lc1 = lc + 1;
    outDimsGHc[0] = lc1; outDimsGHc[1] = lc1;

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
/*
    grid->rr[0] = R0; grid->rr[1] = RCP; grid->rr[2] = 3.0;
    grid->rr[3] = 5.0; grid->rr[4] = 10.0; grid->rr[5] = RSS;	// for now
 */
    if (nr != 1) {
        for (n = 0; n < nr; n++)
            grid->rr[n] = R0 + dr * n;
    } else {
        grid->rr[0] = RS;
        printf("Only one point specified in radial direction, compute source surface.\n");
    }
    
    // Do this for each record set
    for (irec = 0; irec < nrecs; irec++) {

        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }

        /* Input record and data */
        inRec = inRS->records[irec];
        
        inSegG = drms_segment_lookup(inRec, "g");
        inArrayG = drms_segment_read(inSegG, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No G data file found. \n");
            drms_free_array(inArrayG);
            continue;
        }
        inDataG = (float *)inArrayG->data;
        
        inSegH = drms_segment_lookup(inRec, "h");
        inArrayH = drms_segment_read(inSegH, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No H data file found. \n");
            drms_free_array(inArrayH);
            continue;
        }
        inDataH = (float *)inArrayH->data;
        
        car_rot = drms_getkey_int(inRec, "CAR_ROT", &status);

        // Check for lmax, truncate g, h
        n0 = inArrayG->axis[0];		// have_lmax + 1
        n1 = inArrayG->axis[1];		// have_lmax + 1
        if (n0 != n1 || (n0 - 1) != have_lmax)
            DIE("Wrong input data dimension");
        if (want_lmax > have_lmax)
            DIE("Request lmax greater than available lmax");
        lmax = want_lmax;
        lmax1 = lmax + 1;
        g = (float *)(calloc(lmax1 * lmax1, sizeof(float)));
        h = (float *)(calloc(lmax1 * lmax1, sizeof(float)));
        for (int l = 0; l <= lmax; l++) {
            for (int m = 0; m <= l; m++) {
                g[l * lmax1 + m] = inDataG[l * n0 + m];
                h[l * lmax1 + m] = inDataH[l * n0 + m];
         //       printf("g[%d,%d]=%f, h[%d,%d]=%f\n", l, m, g[l * lmax1 + m], l, m, h[l * lmax1 + m]);
            }
        }
        g[0] = 0.0; h[0] = 0.0;	// No monopole
        drms_free_array(inArrayG); drms_free_array(inArrayH);

        /* Output data */
        outArrayBr = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayBt = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayBp = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayGc = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGHc, NULL, &status);
        outArrayHc = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGHc, NULL, &status);
        Br = (float *)outArrayBr->data;
        Bt = (float *)outArrayBt->data;
        Bp = (float *)outArrayBp->data;
        gc = (float *)outArrayGc->data;
        hc = (float *)outArrayHc->data;

        /* This is the working part */
        if (verbflag) SHOW("start computing coefficient... ");
        gh_hccsss(g, h, grid, lmax, gc, hc, lc, apar);
        if (verbflag) SHOW("done.\n");
        if (verbflag) SHOW("start computing field... ");
        hccsss(g, h, gc, hc, grid, Br, Bt, Bp, lmax, lc, apar);
        if (verbflag) SHOW("done.\n");
/*
        for (int ii = 0; ii <= lc; ii++) {
          for (int jj = 0; jj <= lc; jj++) {
           printf("%12.4f", gc[ii*lc1+jj]);
          }
          printf("\n");
        }
        for (int ii = 0; ii <= lc; ii++) {
          for (int jj = 0; jj <= lc; jj++) {
           printf("%12.4f", hc[ii*lc1+jj]);
          }
          printf("\n");
        }
 */
        free(g); free(h);

        /* Output record */
        outRec = outRS->records[irec];
        outSegBr = drms_segment_lookup(outRec, "Br");
        outSegBt = drms_segment_lookup(outRec, "Bt");
        outSegBp = drms_segment_lookup(outRec, "Bp");
        outSegGc = drms_segment_lookup(outRec, "gc");
        outSegHc = drms_segment_lookup(outRec, "hc");

        for (i = 0; i < 3; i++) {
            outSegBr->axis[i] = outArrayBr->axis[i];
            outSegBt->axis[i] = outArrayBt->axis[i];
            outSegBp->axis[i] = outArrayBp->axis[i];
        }
        for (i = 0; i < 2; i++) {
            outSegGc->axis[i] = outArrayGc->axis[i];
            outSegHc->axis[i] = outArrayHc->axis[i];
        }
        outArrayBr->parent_segment = outSegBr;
        outArrayBt->parent_segment = outSegBt;
        outArrayBp->parent_segment = outSegBp;
        outArrayGc->parent_segment = outSegGc;
        outArrayHc->parent_segment = outSegHc;

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
        status = drms_segment_write(outSegGc, outArrayGc, 0);
        if (status) DIE("problem writing Gc");
        drms_free_array(outArrayGc);
        status = drms_segment_write(outSegHc, outArrayHc, 0);
        if (status) DIE("problem writing Hc");
        drms_free_array(outArrayHc);

        /* Set keywords */
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_copykey(outRec, inRec, "LON0");
        drms_setkey_int(outRec, "LMAX", lmax);
        drms_setkey_int(outRec, "LC", lc);
        drms_setkey_float(outRec, "RCP", RCP);
        drms_setkey_float(outRec, "RSS", RSS);
        drms_setkey_float(outRec, "APAR", apar);
        drms_setkey_float(outRec, "CDELT1", -360. / np);
        drms_setkey_float(outRec, "CRPIX1", 1.0);
        drms_setkey_float(outRec, "CRVAL1", car_rot  * 360. - 180. / np);
        drms_setkey_float(outRec, "CRPIX2", nt / 2.0 + 0.5);
        drms_setkey_float(outRec, "CRVAL2", 0.0);
        if (sinlat) {
            drms_setkey_string(outRec, "CTYPE2", "CRLT-CEA");
            drms_setkey_float(outRec, "CDELT2", 2. / nt);
            drms_setkey_string(outRec, "CUNIT2", "Sine Latitude");
        } else {
            drms_setkey_string(outRec, "CTYPE2", "CRLT-CAR");
            drms_setkey_float(outRec, "CDELT2", 180. / nt);
            drms_setkey_string(outRec, "CUNIT2", "degree");
        }
        drms_setkey_float(outRec, "CDELT3", dr);
        drms_setkey_float(outRec, "CRPIX3", 1.0);
        drms_setkey_float(outRec, "CRVAL3", 1.0);
        
        /* Set links */
        ghLink = hcon_lookup_lower(&outRec->links, "HARMONICS");
        if (synopLink) {
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
