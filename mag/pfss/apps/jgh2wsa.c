/*
 * Module name:		jgh2wsa.c
 *
 * Description:		Code for WSA mode using standard PFSS extrapolation.
 *			Take harmonic coefficients and trace open field line
 *			from the source surface, globally and on subearth lines.
 *			Outputs are: global maps at 1R (v, Brss, fte, FP);
 *			subearth series at 1R (v, Brss, fte, FP); 1AU series 
 *			(v, IMF polarity, time). Note that subearth series include
 *			results on grid point north/south to the subearth point as
 *			an estimation of error (Arge 2000), and 1AU series are 
 *			interpolated to a uniform grid, with an error estimation.
 *			Note that results are computed on a uniform latitude grid.
 *			For detailed method, refer to X. Sun's paper on PFSS.
 *
 * Calling:		pfss_pkg.c, fieldline_pkg.c, wsa_pkg.c, suntime.c, timing.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *			IDL code developed by C. N. Arge (Nick.Arge@Kirtland.af.mil)
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v0.0		Jun 14 2000
 *			v1.0		Jul 31 2009
 *			v1.1		Aug 05 2009
 *			v1.1a		Aug 07 2009
 *			v1.1b		Nov 10 2009
 *			v2.0		Apr 19 2009
 *
 * Issues:
 *			v1.0
 *			Keywords (lon0: left most longitude???)
 *			gh are from sinlat maps: might be normalizaition issues (rings)?
 *			No 1AU time series yet
 *			v1.1
 *			Added choice of integrator, RK4 or Euler
 *			Added propagation to 1AU, need to decide on time grids somehow?
 *			Time changed from float to double
 *			Fixed the bug: car_rot, lon0 should be specified by individual rec
 *			For now the 1AU series has only (jd-2440000) as time (what else???)
 *			v1.1a
 *			Added keywords (lon0, etc.), time format?
 *			Calling pfss.h instead
 *			imf changed to geocentric (away +; from -)
 *			v1.1b
 *			free(grid) moved outside the loop, freed components of the grid
 *			v2.0
 *			New I/O, changed time functions using HMI lib
 *
 * Example:
 *			jgh2wsa "in=su_xudong.hmi_gh[1928][0]" out=su_xudong.hmi_wsa LMAX=9
 *
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <astro.h>
#include <sys/time.h>
#include <sys/resource.h>

#define	DTOR	(M_PI / 180.)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#include "pfss.h"
#include "wsa_pkg.h"
//#include "suntime.c"
#include "timing.c"


/* ################## Ephemeris ################## */

// Converts TIME to the solar b angle.
float time2solarb(TIME t)
{
    int crot;
    double L, B;
    HeliographicLocation(t, &crot, &L, &B);
    return (float)B;
}

// Converts Carrington rotation and longitude to TIME
TIME carlon2time(int car, float lon)
{
    return HeliographicTime(car, (double)lon);
}

/* ############################################### */


char *module_name = "jgh2wsa";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_INT, "LMAX", "9", "Maximum l number desired."},
    {ARG_INT, "NP", "144", "Number of Grid points in longitude."},
    {ARG_INT, "NT", "72", "Number of Grid points in latitude."},
    {ARG_INT, "RK4", "0", "Set non zero to use RK4 in field line tracing; default is 0."},
    {ARG_INT, "TSLEN", "120", "Number of points in 1AU series, default is 120."},
    {ARG_DOUBLE, "TSTEP", "14400.0", "Step of 1AU time series in seconds"},
    {ARG_FLOAT, "LON_CMP", "60.0", "CMP lon of the last mag wrt left edge, default is 60."},
    {ARG_DOUBLE, "DAYSINADV", "5.0", "Dt of last point in TS (1AU) wrt CMP lon (source surface)."},
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
    DRMS_Segment_t *outSegGlb_v, *outSegGlb_brss, *outSegGlb_f, *outSegGlb_ph, *outSegGlb_th;
    DRMS_Segment_t *outSegSub_v, *outSegSub_brss, *outSegSub_f, *outSegSub_ph, *outSegSub_th;
    DRMS_Segment_t *outSegTS1AU_vx, *outSegTS1AU_imf;
    DRMS_Array_t *inArrayG, *inArrayH;
    DRMS_Array_t *outArrayGlb_v, *outArrayGlb_brss, *outArrayGlb_f, *outArrayGlb_ph, *outArrayGlb_th;
    DRMS_Array_t *outArraySub_v, *outArraySub_brss, *outArraySub_f, *outArraySub_ph, *outArraySub_th;
    DRMS_Array_t *outArrayTS1AU_vx, *outArrayTS1AU_imf;
    DRMS_Link_t *synopLink, *ghLink;
    float *inDataG, *inDataH;
    float *g, *h;
    float *v_g, *brss_g, *fte_g, *fpph_g, *fpth_g;
    float *v_s, *brss_s, *fte_s, *fpph_s, *fpth_s;
    double *date;
    float *speed;
    int *imf;

    char *inchecklist[] = {"CAR_ROT", "LON0", "LMAX"};
    char *outchecklist[] = {"CAR_ROT", "LMAX"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int verbflag;
    int outDimsGlb[2], outDimsSub[2], outDimsTS1AU[2];
    int lmax, lmax1;
    int want_lmax, have_lmax;
    int np, nt, npnt, np3;
    float dph, dth;
    int rk4;
    int tslen;
    double tstep, t_frst, t_last, t_cmp;
    float lon_cmp;
    double daysinadv;
    int car_rot;
    float lon0;

    int i, j, n, itest;
    int l, m;
    int n0, n1;

    struct Grid *grid;
    struct Gridsub *gridsub;

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
    np = params_get_int(&cmdparams, "NP");
    nt = params_get_int(&cmdparams, "NT");
    rk4 = params_get_int(&cmdparams, "RK4");
    tslen = params_get_int(&cmdparams, "TSLEN");
    tstep = params_get_double(&cmdparams, "TSTEP");
    lon_cmp = params_get_float(&cmdparams, "LON_CMP");
    daysinadv = params_get_double(&cmdparams, "DAYSINADV");
    verbflag = params_get_int(&cmdparams, "VERB");

    npnt = np * nt; np3 = np * 3;
    dph = 360. / np;
    dth = 180. / nt;

    /* Set up dimesions for out put */
    // Global[np*nt*5]: v, Brss, fte, ph_fp, th_fp
    outDimsGlb[0] = np; outDimsGlb[1] = nt;
    // Sub[np*3*5]: v, Brss, fte, ph_fp, th_fp
    outDimsSub[0] = np; outDimsSub[1] = 3;
    // TS1AU[tslen*3*3]: time, v, IMF polarity
    outDimsTS1AU[0] = tslen; outDimsTS1AU[1] = 3;

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

    /* Specify global grid */
    grid = (struct Grid *)(malloc(sizeof(struct Grid)));
    grid->np = np; grid->nt = nt;
    grid->ph = (float *)(calloc(np, sizeof(float)));
    grid->th = (float *)(calloc(nt, sizeof(float)));
    for (i = 0; i < np; i++)
        grid->ph[i] = (i + 0.5) / np * 2 * M_PI;
    for (j = 0; j < nt; j++)
        grid->th[j] = (nt - j - 0.5) / nt * M_PI;

    /* Specify subearth grid: renewed for each record */
    gridsub = (struct Gridsub *)(malloc(sizeof(struct Gridsub)));
    gridsub->np = np; gridsub->nt = nt;
    gridsub->time = (double *)(calloc(np, sizeof(double)));
    gridsub->b0 = (float *)(calloc(np, sizeof(float)));
    gridsub->ph = (float *)(calloc(np, sizeof(float)));
    gridsub->th = (float *)(calloc(np3, sizeof(float)));

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
        inSegG = drms_segment_lookup(inRec, "g");
        inArrayG = drms_segment_read(inSegG, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No data file found. \n");
            drms_free_array(inArrayG);
            continue;
        }
        inDataG = (float *)inArrayG->data;
        inSegH = drms_segment_lookup(inRec, "h");
        inArrayH = drms_segment_read(inSegH, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No data file found. \n");
            drms_free_array(inArrayH);
            continue;
        }
        inDataH = (float *)inArrayH->data;

        /* Parameters */
        car_rot = drms_getkey_int(inRec, "CAR_ROT", &status);
        lon0 = drms_getkey_float(inRec, "LON0", &status);
        for (i = 0; i < np; i++) {
            gridsub->ph[i] = (i + 0.5) / np * 2 * M_PI;
            gridsub->time[i] = carlon2time(car_rot, lon0 + gridsub->ph[i] / DTOR);	// Time
            gridsub->b0[i] = time2solarb(gridsub->time[i]);			// Get b angle
            gridsub->th[np + i] = M_PI / 2. - gridsub->b0[i] * DTOR;		// subearth
            gridsub->th[i] = gridsub->th[np + i] + dth * DTOR;		// dth south
            gridsub->th[2 * np + i] = gridsub->th[np + i] - dth * DTOR;	// dth north
        }

        // Check for lmax, truncate g, h
        n0 = inArrayG->axis[0];		// have_lmax + 1
        n1 = inArrayG->axis[1];		// have_lmax + 1
        if (n0 != n1 || (n0 - 1) != have_lmax || inArrayH->axis[0] != inArrayH->axis[1] || inArrayH->axis[0] != n0)
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
            }
        }
        g[0] = 0.0; h[0] = 0.0;	// No monopole
        drms_free_array(inArrayG); drms_free_array(inArrayH);

        /* Output data */
        outArrayGlb_v = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGlb, NULL, &status);
        outArrayGlb_brss = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGlb, NULL, &status);
        outArrayGlb_f = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGlb, NULL, &status);
        outArrayGlb_ph = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGlb, NULL, &status);
        outArrayGlb_th = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsGlb, NULL, &status);

        outArraySub_v = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsSub, NULL, &status);
        outArraySub_brss = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsSub, NULL, &status);
        outArraySub_f = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsSub, NULL, &status);
        outArraySub_ph = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsSub, NULL, &status);
        outArraySub_th = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsSub, NULL, &status);

        outArrayTS1AU_vx = drms_array_create(DRMS_TYPE_FLOAT, 2, outDimsTS1AU, NULL, &status);
        outArrayTS1AU_imf = drms_array_create(DRMS_TYPE_INT, 2, outDimsTS1AU, NULL, &status);

        v_g = (float *)outArrayGlb_v->data;			// Predicted speed (global)
        brss_g = (float *)outArrayGlb_brss->data;	// Source surface field (global)
        fte_g = (float *)outArrayGlb_f->data;		// Expansion factor (global)
        fpph_g = (float *)outArrayGlb_ph->data;		// Foot point phi (global)
        fpth_g = (float *)outArrayGlb_th->data;		// Foot point theta (global)

        v_s = (float *)outArraySub_v->data;			// Predicted speed (subearth)
        brss_s = (float *)outArraySub_brss->data;	// Source surface field (subearth)
        fte_s = (float *)outArraySub_f->data;		// Expansion factor (subearth)
        fpph_s = (float *)outArraySub_ph->data;		// Foot point phi (subearth)
        fpth_s = (float *)outArraySub_th->data;		// Foot point theta (subearth)

        speed = (float *)outArrayTS1AU_vx->data;		// Speed at 1AU
        imf = (int *)outArrayTS1AU_imf->data;			// Imf polarity at 1AU

        date = (double *)calloc(tslen, sizeof(double));			// Time series at 1AU

        /* Specify time series */
        t_cmp = carlon2time(car_rot, lon0 + lon_cmp);
        t_last = t_cmp + (daysinadv + 0.5) * SECSINDAY;
        t_last = (int)(t_last / tstep) * tstep;
        t_frst = t_last - (tslen - 1) * tstep;
        for (i = 0; i < tslen; i++) {
            date[i] = t_last - (tslen - i - 1) * tstep;
        }


        /* This is the working part */
        
        // Get parameters at source surface
        if (verbflag) SHOW("tracing from source surface... ");
        wsa_ss(g, h, lmax, rk4, grid, gridsub, 
                v_g, brss_g, fte_g, fpph_g, fpth_g, 
                v_s, brss_s, fte_s, fpph_s, fpth_s);
        if (verbflag) SHOW("done.\n");
        
        // Get time series at 1AU
        // Note v_s and brss_s are in reverse order, temporally
        if (verbflag) SHOW("propagating to 1AU... ");
        wsa_1AU(gridsub, v_s, brss_s, date, speed, imf, tslen);
        if (verbflag) SHOW("done.\n");

        /* Outpur record */
        outRec = outRS->records[irec];
        outSegGlb_v = drms_segment_lookup(outRec, "Glb_v");
        outSegGlb_brss = drms_segment_lookup(outRec, "Glb_brss");
        outSegGlb_f = drms_segment_lookup(outRec, "Glb_f");
        outSegGlb_ph = drms_segment_lookup(outRec, "Glb_ph");
        outSegGlb_th = drms_segment_lookup(outRec, "Glb_th");
        outSegSub_v = drms_segment_lookup(outRec, "Sub_v");
        outSegSub_brss = drms_segment_lookup(outRec, "Sub_brss");
        outSegSub_f = drms_segment_lookup(outRec, "Sub_f");
        outSegSub_ph = drms_segment_lookup(outRec, "Sub_ph");
        outSegSub_th = drms_segment_lookup(outRec, "Sub_th");
        outSegTS1AU_vx = drms_segment_lookup(outRec, "TS1AU_vx");
        outSegTS1AU_imf = drms_segment_lookup(outRec, "TS1AU_imf");
        for (i = 0; i < 2; i++) {
            outSegGlb_v->axis[i] = outArrayGlb_v->axis[i];
            outSegGlb_brss->axis[i] = outArrayGlb_brss->axis[i];
            outSegGlb_f->axis[i] = outArrayGlb_f->axis[i];
            outSegGlb_ph->axis[i] = outArrayGlb_ph->axis[i];
            outSegGlb_th->axis[i] = outArrayGlb_th->axis[i];
            outSegSub_v->axis[i] = outArraySub_v->axis[i];
            outSegSub_brss->axis[i] = outArraySub_brss->axis[i];
            outSegSub_f->axis[i] = outArraySub_f->axis[i];
            outSegSub_ph->axis[i] = outArraySub_ph->axis[i];
            outSegSub_th->axis[i] = outArraySub_th->axis[i];
            outSegTS1AU_vx->axis[i] = outArrayTS1AU_vx->axis[i];
            outSegTS1AU_imf->axis[i] = outArrayTS1AU_imf->axis[i];
        }
        outArrayGlb_v->parent_segment = outSegGlb_v;
        outArrayGlb_brss->parent_segment = outSegGlb_brss;
        outArrayGlb_f->parent_segment = outSegGlb_f;
        outArrayGlb_ph->parent_segment = outSegGlb_ph;
        outArrayGlb_th->parent_segment = outSegGlb_th;
        outArraySub_v->parent_segment = outSegSub_v;
        outArraySub_brss->parent_segment = outSegSub_brss;
        outArraySub_f->parent_segment = outSegSub_f;
        outArraySub_ph->parent_segment = outSegSub_ph;
        outArraySub_th->parent_segment = outSegSub_th;
        outArrayTS1AU_vx->parent_segment = outSegTS1AU_vx;
        outArrayTS1AU_imf->parent_segment = outSegTS1AU_imf;

        /* Result writing */
        status = drms_segment_write(outSegGlb_v, outArrayGlb_v, 0);
        if (status) DIE("problem writing Glb_v");
        drms_free_array(outArrayGlb_v);
        status = drms_segment_write(outSegGlb_brss, outArrayGlb_brss, 0);
        if (status) DIE("problem writing Glb_brss");
        drms_free_array(outArrayGlb_brss);
        status = drms_segment_write(outSegGlb_f, outArrayGlb_f, 0);
        if (status) DIE("problem writing Glb_f");
        drms_free_array(outArrayGlb_f);
        status = drms_segment_write(outSegGlb_ph, outArrayGlb_ph, 0);
        if (status) DIE("problem writing Glb_ph");
        drms_free_array(outArrayGlb_ph);
        status = drms_segment_write(outSegGlb_th, outArrayGlb_th, 0);
        if (status) DIE("problem writing Glb_th");
        drms_free_array(outArrayGlb_th);
        status = drms_segment_write(outSegSub_v, outArraySub_v, 0);
        if (status) DIE("problem writing Sub_v");
        drms_free_array(outArraySub_v);
        status = drms_segment_write(outSegSub_brss, outArraySub_brss, 0);
        if (status) DIE("problem writing Sub_brss");
        drms_free_array(outArraySub_brss);
        status = drms_segment_write(outSegSub_f, outArraySub_f, 0);
        if (status) DIE("problem writing Sub_f");
        drms_free_array(outArraySub_f);
        status = drms_segment_write(outSegSub_ph, outArraySub_ph, 0);
        if (status) DIE("problem writing Sub_ph");
        drms_free_array(outArraySub_ph);
        status = drms_segment_write(outSegSub_th, outArraySub_th, 0);
        if (status) DIE("problem writing Sub_th");
        drms_free_array(outArraySub_th);
        status = drms_segment_write(outSegTS1AU_vx, outArrayTS1AU_vx, 0);
        if (status) DIE("problem writing TS1AU_vx");
        drms_free_array(outArrayTS1AU_vx);
        status = drms_segment_write(outSegTS1AU_imf, outArrayTS1AU_imf, 0);
        if (status) DIE("problem writing TS1AU_imf");
        drms_free_array(outArrayTS1AU_imf);

        /* Set keywords */
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);

        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_copykey(outRec, inRec, "LON0");
        drms_setkey_int(outRec, "LMAX", lmax);
        drms_setkey_float(outRec, "RSS", RS);
        drms_setkey_int(outRec, "RK4", rk4);
        drms_setkey_float(outRec, "FL_STEP", STEP);
        drms_setkey_float(outRec, "LON_CMP", lon_cmp);

        drms_setkey_time(outRec, "T_CMP", t_cmp);
        drms_setkey_double(outRec, "DAYSINADV", daysinadv);
        drms_setkey_time(outRec, "T_START", t_frst);
        drms_setkey_time(outRec, "T_END", t_last);
        drms_setkey_double(outRec, "T_STEP", tstep);

        drms_setkey_float(outRec, "CRPIX1_Glb", 1.0);
        drms_setkey_float(outRec, "CRVAL1_Glb", car_rot  * 360. - 180. / np);
        drms_setkey_float(outRec, "CDELT1_Glb", - 360. / np);
        drms_setkey_float(outRec, "CRPIX2_Glb", nt / 2.0 + 0.5);
        drms_setkey_float(outRec, "CRVAL2_Glb", 0.0);       
        drms_setkey_float(outRec, "CDELT2_Glb", 180. / nt);

        drms_setkey_float(outRec, "CRPIX1_Sub", 1.0);
        drms_setkey_float(outRec, "CRVAL1_Sub", car_rot  * 360. - 180. / np);
        drms_setkey_float(outRec, "CDELT1_Sub", - 360. / np);
        drms_setkey_float(outRec, "CRPIX2_Sub", 2.0);
        drms_setkey_float(outRec, "CRVAL2_FRST_Sub", gridsub->b0[0]); 
        drms_setkey_float(outRec, "CRVAL2_LAST_Sub", gridsub->b0[nt - 1]);     
        drms_setkey_float(outRec, "CDELT2_Glb", 180. / nt);

        drms_setkey_float(outRec, "CRPIX1_TS1AU", 1.0);
        drms_setkey_time(outRec, "CRVAL1_TS1AU", t_frst);
        drms_setkey_double(outRec, "CDELT1_TS1AU", tstep);
        drms_setkey_float(outRec, "CRPIX2_TS1AU", 2.0);
        drms_setkey_float(outRec, "CRVAL2_TS1AU", 0.0); 
        drms_setkey_float(outRec, "CDELT2_TS1AU", 180. / nt);

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

    // Clean up
    free(grid->th); free(grid->ph);
    free(grid);
    free(gridsub->time); free(gridsub->b0); free(gridsub->th); free(gridsub->ph);
    free(gridsub);

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
