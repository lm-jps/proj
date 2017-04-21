/*
 * Module name:		test_d4vm.c
 *
 * Description:
 *			DAVE4VM C version
 *
 * Calling:
 *					
 *
 * Original by:		Peter Schuck, Jacob Hageman
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Mar 01 2010
 *			v1.1		Mar 15 2010
 *

 * Issues:
 *			v1.0
 *			1. Designed to deal with a continuous time series, computing result
 *			with each pair of neighbouring records.
 *			2. Added "cadence" keyword to stride over a minimum fixed period between
 *			two records, so the result can be computed at a lower cadence.
 *			3. Commented out #infdef USE_INTEL_LIB and DEBUG, seems not properly set up
 *			4. Will have to make sure the stack size is big enough (at least ~70MB)
 *			run "limit s u" to set unlimited
 *			v1.1
 *			Updated according to Jacob's new code
 *
 *
 * Example:
 *			test_d4vm_old 'in=su_xudong.d4vmin[2004.08.31_23:55:00_TAI/10m]'
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fstats.h>

#define PI	(M_PI)
#define	DTOR	(PI / 180.)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#include "timing.c"



/* ################## Prototypes ################## */

extern void d4vm_preproc_(float *bx0, float *bx1, 
                          float *by0, float *by1, 
                          float *bz0, float *bz1, 
                          double *dx, double *dy, double *dt, 
                          float *bx, float *bxx, float *bxy, 
                          float *by, float *byx, float *byy, 
                          float *bz, float *bzx, float *bzy, float *bzt,
                          long *xsize, long *ysize);

void d4vm_kernel_(double *dx, double *dy, 
                  long *ksize, float *k_th,
                  float *k_x, float *k_y, 
                  float *k_xx, float *k_yy, float *k_xy);

void d4vm_matrix_(float *a, 
                  float *bx, float *bxx, float *bxy, 
                  float *by, float *byx, float *byy, 
                  float *bz, float *bzx, float *bzy, float *bzt, 
                  float *k_th, float *k_x, float *k_y,
                  float *k_xx, float *k_yy, float *k_xy,
                  long *xsize, long *ysize, long *ksize);

void d4vm_solver_(float *a, 
                  float *u0, float *v0, float *w0, 
                  float *ux, float *vx, float *wx, 
                  float *uy, float *vy, float *wy,
                  long *xsize, long *ysize, double *threshold);


char *module_name = "dave4vm4velocity";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", "hmi.sharp_cea_720s", "Input data series."},
    {ARG_STRING, "out", "su_yang.hmi_dave4vm4velo",  "Output data series."},
    {ARG_INT, "ksz0", "19", "Width of window for convolution, must be odd."},
    {ARG_INT, "ksz1", "19", "Height of window for convolution, must be odd."},
    {ARG_DOUBLE, "threshold", "0.", "Threshold for resolving aperture problem."},
    {ARG_DOUBLE, "cadence", "240.0", "Minimum stride between two records."},
    {ARG_INT, "VERB", "2", "Level of verbosity: 0=errors/warnings; 1=status; 2=all"},
    {ARG_END}
};


/* ################## Main Module ################## */

int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS;
    int irec, nrecs;
    DRMS_Record_t *inRec0, *inRec1, *outRec;
    DRMS_Segment_t *inSeg0, *inSeg1;
    DRMS_Segment_t *outSegV; 
    DRMS_Array_t *inArrayX0, *inArrayY0, *inArrayZ0, *inArrayX1, *inArrayY1, *inArrayZ1;
    DRMS_Array_t *outArrayV, *outArrayVy; 
    float *outDataVx, *outDataVy, *outDataVz; 
    int harpnum;
    
    float *Bx0, *Bx1, *By0, *By1, *Bz0, *Bz1;				// ancillary
    float *Bx, *Bxx, *Bxy, *By, *Byx, *Byy, *Bz, *Bzx, *Bzy, *Bzt;	// processed
    float *k_th, *k_x, *k_y, *k_xx, *k_yy, *k_xy;			// convolution kernels
    float *a;								// matrix
    float *Vx, *Vy, *Vz;						// outputs
    float *DVxDx, *DVyDx, *DVzDx;
    float *DVxDy, *DVyDy, *DVzDy;
    float onedegree4L = 6.955e5 * DTOR; // 1 degree corresponse to km

    char *inchecklist[] = {"T_REC","T_OBS","HARPNUM"};		// For now
    char *outchecklist[] = {"T_REC", "K0", "K1", "THRESH","HARPNUM"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int outDims[2] = {0, 0};
    long nx, ny, nxny, ksize[2];	// Seems compatible with long in Fortran?
    int verbflag;
    double threshold, cadence;
    TIME t0, t1;	// Need to check format
    double dx, dy, dt;

    int i, l, m, itest;
    int nsuccess = 0;
    
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
    ksize[0] = params_get_int(&cmdparams, "ksz0");
    ksize[1] = params_get_int(&cmdparams, "ksz1");
    threshold = params_get_double(&cmdparams, "threshold");
    cadence = params_get_double(&cmdparams, "cadence");
    verbflag = params_get_int(&cmdparams, "VERB");
    for (i = 0; i < 2; i++) {
        if (ksize[i] % 2 != 1) {
            DIE("Kernel size must be odd");
        }
    }

    /* Open input */
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    if (inRS->n < 2) DIE("At least 2 records are needed");
    nrecs = inRS->n;

//printf("nrecs=%d\n", nrecs);

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

    /* Read in the first record and save it */
    inRec0 = inRS->records[0];
    t0 = drms_getkey_time(inRec0, "T_REC", NULL);	// NEED TO CHECK!!!!!
    harpnum = drms_getkey_int(inRec0, "HARPNUM", NULL);
    
    inSeg0 = drms_segment_lookup(inRec0, "Bp");
    inArrayX0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    Bx0 = (float *)inArrayX0->data;
    nx = inArrayX0->axis[0]; ny = inArrayX0->axis[1]; nxny = nx * ny;
    outDims[0] = nx; outDims[1] = ny;
    inSeg0 = drms_segment_lookup(inRec0, "Bt");
    inArrayY0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    By0 = (float *)inArrayY0->data;
    for (int ii = 0; ii < nxny; ii++) By0[ii] *= -1.;
    inSeg0 = drms_segment_lookup(inRec0, "Br");
    inArrayZ0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    Bz0 = (float *)inArrayZ0->data;

    /* Do this for record #1 to #nrecs-1 */
    for (irec = 1; irec < nrecs; irec++)
    {
        /* Input record and data */
        inRec1 = inRS->records[irec];
        t1 = drms_getkey_time(inRec1, "T_REC", NULL);	// NEED TO CHECK!!!!!
        dt = t1 - t0;
        if (dt < cadence) {
            printf("Cadence higher than limit, skipping current record. \n");
            drms_free_array(inArrayX1);
            continue;
        }
        if (harpnum != drms_getkey_time(inRec1, "HARPNUM", NULL)) {
            printf("HARP number incorrect, skipping current record. \n");
            drms_free_array(inArrayX1);
            continue;
        }

        inSeg1 = drms_segment_lookup(inRec1, "Bp");	// Single segment
        inArrayX1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
//printf("t0=%f, t1=%f\n", t0, t1);
//printf("dt=%f\n",t1-t0);
        Bx1 = (float *)inArrayX1->data;
        if (inArrayX1->axis[0] != nx || inArrayX1->axis[1] != ny)
            DIE("Wrong data dimension for current record, stop. \n"); 
        inSeg1 = drms_segment_lookup(inRec1, "Bt");
        inArrayY1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
        By1 = (float *)inArrayY1->data;
        for (int ii = 0; ii < nxny; ii++) By1[ii] *= -1.;
        inSeg1 = drms_segment_lookup(inRec1, "Br");
        inArrayZ1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
        Bz1 = (float *)inArrayZ1->data;
        
        /* Time measure */
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", nsuccess);
        }

        /* ======================== */
        /* This is the working part.*/
        /* ======================== */

        // Preproc
        Bx = (float *)calloc(sizeof(float), nxny);
        Bxx = (float *)calloc(sizeof(float), nxny);
        Bxy = (float *)calloc(sizeof(float), nxny);
        By = (float *)calloc(sizeof(float), nxny);
        Byx = (float *)calloc(sizeof(float), nxny);
        Byy = (float *)calloc(sizeof(float), nxny);
        Bz = (float *)calloc(sizeof(float), nxny);
        Bzx = (float *)calloc(sizeof(float), nxny);
        Bzy = (float *)calloc(sizeof(float), nxny);
        Bzt = (float *)calloc(sizeof(float), nxny);

        DVxDx = (float *)calloc(sizeof(float), nxny);
        DVyDx = (float *)calloc(sizeof(float), nxny);
        DVzDx = (float *)calloc(sizeof(float), nxny);
        DVxDy = (float *)calloc(sizeof(float), nxny);
        DVyDy = (float *)calloc(sizeof(float), nxny);
        DVzDy = (float *)calloc(sizeof(float), nxny);

        // Kernel
        k_th = (float *)calloc(sizeof(float), ksize[0] * ksize[1]);
        k_x = (float *)calloc(sizeof(float), ksize[0] * ksize[1]);
        k_y = (float *)calloc(sizeof(float), ksize[0] * ksize[1]);
        k_xx = (float *)calloc(sizeof(float), ksize[0] * ksize[1]);
        k_yy = (float *)calloc(sizeof(float), ksize[0] * ksize[1]);
        k_xy = (float *)calloc(sizeof(float), ksize[0] * ksize[1]);

        // Matrix
        a = (float *)calloc(sizeof(float), 10 * 10 * nxny);

        // Output
        Vx = (float *)calloc(sizeof(float), nxny);
        Vy = (float *)calloc(sizeof(float), nxny);
        Vz = (float *)calloc(sizeof(float), nxny);

        dx = onedegree4L * drms_getkey_float(inRec0, "CDELT1", &status);
        dy = onedegree4L * drms_getkey_float(inRec0, "CDELT2", &status);

//printf("dx = %f, dy = %f\n", dx, dy);
//SHOW("111\n");
        // Preprocess data
        d4vm_preproc_(Bx0, Bx1, By0, By1, Bz0, Bz1, 
                      &dx, &dy, &dt, 
                      Bx, Bxx, Bxy, By, Byx, Byy, Bz, Bzx, Bzy, Bzt,
                      &nx, &ny);
//SHOW("112\n");
        // Generate kernel
        d4vm_kernel_(&dx, &dy, ksize,
                     k_th, k_x, k_y, k_xx, k_yy, k_xy);
//SHOW("113\n");

        // Process data
        d4vm_matrix_(a, 
                     Bx, Bxx, Bxy, By, Byx, Byy, Bz, Bzx, Bzy, Bzt,
                     k_th, k_x, k_y, k_xx, k_yy, k_xy,
                     &nx, &ny, ksize);
//SHOW("114\n");
        // Solve matrix
        d4vm_solver_(a, 
                     Vx, Vy, Vz, 
                     DVxDx, DVyDx, DVzDx, 
                     DVxDy, DVyDy, DVzDy, 
                     &nx, &ny, &threshold);
//SHOW("115\n");

        /* Output record, create one at a time */
        outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
        if (status) DIE("Output record not created");

        outSegV = drms_segment_lookup(outRec, "Vp");
        outArrayV = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, Vx, &status);
        outArrayV->israw = 0;
        outArrayV->bzero = 0.0;
        outArrayV->bscale = 0.0001;
        /* Result writing */
        status = drms_segment_write(outSegV, outArrayV, 0);
        if (status) DIE("Problem writing file V");
        free(outArrayV);

        outSegV = drms_segment_lookup(outRec, "Vt");
        int ii;
        for (ii = 0; ii < nxny; ii++) Vy[ii] *= -1.0;
        outArrayV = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, Vy, &status);
        outArrayV->israw = 0;
        outArrayV->bzero = 0.0;
        outArrayV->bscale = 0.0001;
        /* Result writing */
        status = drms_segment_write(outSegV, outArrayV, 0);
        if (status) DIE("Problem writing file V");
        free(outArrayV);

        outSegV = drms_segment_lookup(outRec, "Vr");
        outArrayV = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, Vz, &status);
        outArrayV->israw = 0;
        outArrayV->bzero = 0.0;
        outArrayV->bscale = 0.0001;
        /* Result writing */
        status = drms_segment_write(outSegV, outArrayV, 0);
        if (status) DIE("Problem writing file V");
        free(outArrayV);

        /* Set keywords */
        drms_copykey(outRec, inRec0, "HARPNUM");
        drms_copykey(outRec, inRec0, "T_REC");
        drms_copykey(outRec, inRec0, "T_OBS");
// add aditional keywords
        drms_copykey(outRec, inRec0, "DATE__OBS");
        drms_copykey(outRec, inRec0, "CRPIX1");
        drms_copykey(outRec, inRec0, "CRPIX2");
        drms_copykey(outRec, inRec0, "CRVAL1");
        drms_copykey(outRec, inRec0, "CRVAL2");
        drms_copykey(outRec, inRec0, "CDELT1");
        drms_copykey(outRec, inRec0, "CDELT2");
        drms_copykey(outRec, inRec0, "IMCRPIX1");
        drms_copykey(outRec, inRec0, "IMCRPIX2");
        drms_copykey(outRec, inRec0, "IMCRVAL1");
        drms_copykey(outRec, inRec0, "IMCRVAL2");
        drms_copykey(outRec, inRec0, "CROTA2");
        drms_copykey(outRec, inRec0, "CRDER1");
        drms_copykey(outRec, inRec0, "CRDER2");
        drms_copykey(outRec, inRec0, "CSYSER1");
        drms_copykey(outRec, inRec0, "CSYSER2");
        drms_copykey(outRec, inRec0, "DSUN_OBS");
        drms_copykey(outRec, inRec0, "CRLN_OBS");
        drms_copykey(outRec, inRec0, "CRLT_OBS");
        drms_copykey(outRec, inRec0, "CAR_ROT");
        drms_copykey(outRec, inRec0, "OBS_VR");
        drms_copykey(outRec, inRec0, "OBS_VW");
        drms_copykey(outRec, inRec0, "OBS_VN");
        drms_copykey(outRec, inRec0, "RSUN_OBS");
        drms_copykey(outRec, inRec0, "INSTRUME");
        drms_copykey(outRec, inRec0, "CAMERA");
        drms_copykey(outRec, inRec0, "LATDTMIN");
        drms_copykey(outRec, inRec0, "LONDTMIN");
        drms_copykey(outRec, inRec0, "LATDTMAX");
        drms_copykey(outRec, inRec0, "LONDTMAX");
        drms_copykey(outRec, inRec0, "OMEGA_DT");
        drms_copykey(outRec, inRec0, "LAT_MIN");
        drms_copykey(outRec, inRec0, "LON_MIN");
        drms_copykey(outRec, inRec0, "LAT_MAX");
        drms_copykey(outRec, inRec0, "LON_MAX");
        drms_copykey(outRec, inRec0, "NOAA_AR");
        drms_copykey(outRec, inRec0, "NOAA_NUM");
        drms_copykey(outRec, inRec0, "NOAA_ARS");

        drms_setkey_time(outRec, "T_REC_STAR", drms_getkey_time(inRec0, "T_REC", NULL));
        drms_setkey_time(outRec, "T_OBS_STAR", drms_getkey_time(inRec0, "T_OBS", NULL));
        drms_setkey_time(outRec, "T_REC_STOP", drms_getkey_time(inRec1, "T_REC", NULL));
        drms_setkey_time(outRec, "T_OBS_STOP", drms_getkey_time(inRec1, "T_OBS", NULL));

// add statistical keywords

        double statMin, statMax, statMedn, statMean, statSig, statSkew, statKurt;
        int statNgood, ipixels;
         if (fstats(outDims[0] * outDims[1], Vx, &statMin, &statMax, &statMedn, &statMean, &statSig,
            &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");

        ipixels = outDims[0] * outDims[1];
        drms_setkey_int(outRec, "TOTVALS_1", ipixels);
        drms_setkey_int(outRec, "DATAVALS_1", statNgood);
        ipixels = outDims[0] * outDims[1]-statNgood;
        drms_setkey_int(outRec, "MISSVALS_1", ipixels);
        drms_setkey_double(outRec, "DATAMIN_1", statMin);
        drms_setkey_double(outRec, "DATAMAX_1", statMax);
        drms_setkey_double(outRec, "DATAMEDN_1", statMedn);
        drms_setkey_double(outRec, "DATAMEAN_1", statMean);
        drms_setkey_double(outRec, "DATARMS_1", statSig);
        drms_setkey_double(outRec, "DATASKEW_1", statSkew);
        drms_setkey_double(outRec, "DATAKURT_1", statKurt);

         if (fstats(outDims[0] * outDims[1], Vy, &statMin, &statMax, &statMedn, &statMean, &statSig,
            &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");

        ipixels = outDims[0] * outDims[1];
        drms_setkey_int(outRec, "TOTVALS_2", ipixels);
        drms_setkey_int(outRec, "DATAVALS_2", statNgood);
        ipixels = outDims[0] * outDims[1]-statNgood;
        drms_setkey_int(outRec, "MISSVALS_2", ipixels);
        drms_setkey_double(outRec, "DATAMIN_2", statMin);
        drms_setkey_double(outRec, "DATAMAX_2", statMax);
        drms_setkey_double(outRec, "DATAMEDN_2", statMedn);
        drms_setkey_double(outRec, "DATAMEAN_2", statMean);
        drms_setkey_double(outRec, "DATARMS_2", statSig);
        drms_setkey_double(outRec, "DATASKEW_2", statSkew);
        drms_setkey_double(outRec, "DATAKURT_2", statKurt);

         if (fstats(outDims[0] * outDims[1], Vz, &statMin, &statMax, &statMedn, &statMean, &statSig,
            &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");

        ipixels = outDims[0] * outDims[1];
        drms_setkey_int(outRec, "TOTVALS_3", ipixels);
        drms_setkey_int(outRec, "DATAVALS_3", statNgood);
        ipixels = outDims[0] * outDims[1]-statNgood;
        drms_setkey_int(outRec, "MISSVALS_3", ipixels);
        drms_setkey_double(outRec, "DATAMIN_3", statMin);
        drms_setkey_double(outRec, "DATAMAX_3", statMax);
        drms_setkey_double(outRec, "DATAMEDN_3", statMedn);
        drms_setkey_double(outRec, "DATAMEAN_3", statMean);
        drms_setkey_double(outRec, "DATARMS_3", statSig);
        drms_setkey_double(outRec, "DATASKEW_3", statSkew);
        drms_setkey_double(outRec, "DATAKURT_3", statKurt);

// end of statistical keywords


        drms_setkey_int(outRec, "K0", ksize[0]);
        drms_setkey_int(outRec, "K1", ksize[1]);
        drms_setkey_double(outRec, "THRESH", threshold);
        drms_close_record(outRec, DRMS_INSERT_RECORD);

        /* Clean up, move to next record */
        drms_free_array(inArrayX0);
        drms_free_array(inArrayY0);
        drms_free_array(inArrayZ0);
        free(Bx); free(Bxx); free(Bxy);
        free(By); free(Byx); free(Byy);
        free(Bz); free(Bzx); free(Bzy); free(Bzt);
        free(DVxDx); free(DVyDx); free(DVzDx);
        free(DVxDy); free(DVyDy); free(DVzDy);
        free(a);
        free(k_th); free(k_x); free(k_y);
        free(k_xx); free(k_yy); free(k_xy);
        free(Vy); free(Vz); free(Vx);

// set the second data to the first for next calculation

        inRec0 = inRec1;
        inSeg0 = inSeg1;
        inArrayX0 = inArrayX1;
        inArrayY0 = inArrayY1;
        inArrayZ0 = inArrayZ1;
        Bx0 = Bx1;
        By0 = By1;
        Bz0 = Bz1;
        t0 = t1;
        
        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                     nsuccess, wt - wt1, ct - ct1);
        }

        nsuccess++;
    }

    drms_free_array(inArrayX1);
    drms_free_array(inArrayY1);
    drms_free_array(inArrayZ1);
    drms_close_records(inRS, DRMS_FREE_RECORD);
    return(DRMS_SUCCESS);
}
