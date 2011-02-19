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
 *			Undated according to Jacob's new code
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


char *module_name = "test_d4vm";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", "su_xudong.d4vmin", "Input data series."},
    {ARG_STRING, "out", "su_xudong.d4vmout",  "Output data series."},
    {ARG_INT, "ksz0", "23", "Width of window for convolution, must be odd."},
    {ARG_INT, "ksz1", "23", "Height of window for convolution, must be odd."},
    {ARG_DOUBLE, "dx", "348.36", "X scaling."},
    {ARG_DOUBLE, "dy", "348.36", "Y scaling."},
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
    DRMS_Segment_t *outSegV, *outSegDVDx, *outSegDVDy;
    DRMS_Array_t *inArrayX0, *inArrayY0, *inArrayZ0, *inArrayX1, *inArrayY1, *inArrayZ1;
    DRMS_Array_t *outArrayV, *outArrayDVDx, *outArrayDVDy;
    float *inData0, *inData1, *outDataV, *outDataDVDx, *outDataDVDy;

    float *Bx0, *Bx1, *By0, *By1, *Bz0, *Bz1;				// ancillary
    float *Bx, *Bxx, *Bxy, *By, *Byx, *Byy, *Bz, *Bzx, *Bzy, *Bzt;	// processed
    float *k_th, *k_x, *k_y, *k_xx, *k_yy, *k_xy;			// convolution kernels
    float *a;								// matrix
    float *Vx, *Vy, *Vz;						// outputs
    float *DVxDx, *DVyDx, *DVzDx;
    float *DVxDy, *DVyDy, *DVzDy;

    char *inchecklist[] = {"T_OBS"};		// For now
    char *outchecklist[] = {"T_OBS", "T_NEXT", "K0", "K1", "THRESH"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int outDims[3] = {0, 0, 3};
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
    dx = params_get_double(&cmdparams, "dx");
    dy = params_get_double(&cmdparams, "dy");
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
    inSeg0 = drms_segment_lookupnum(inRec0, 0);
    inArrayX0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    if (status) {
        drms_free_array(inArrayX0);
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("No data file found for record 0, stop");
    }
    t0 = drms_getkey_time(inRec0, "T_OBS", NULL);	// NEED TO CHECK!!!!!
//printf("t0=%f\n", t0);
//printf("status=%d\n", status);
    Bx0 = (float *)inArrayX0->data;
    nx = inArrayX0->axis[0]; ny = inArrayX0->axis[1]; nxny = nx * ny;
    outDims[0] = nx; outDims[1] = ny;
    inSeg0 = drms_segment_lookupnum(inRec0, 1);
    inArrayY0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    By0 = (float *)inArrayY0->data;
    inSeg0 = drms_segment_lookupnum(inRec0, 2);
    inArrayZ0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    Bz0 = (float *)inArrayZ0->data;

    /* Do this for record #1 to #nrecs-1 */
    for (irec = 1; irec < nrecs; irec++)
    {
        /* Input record and data */
        inRec1 = inRS->records[irec];
        inSeg1 = drms_segment_lookupnum(inRec1, 0);	// Single segment
        inArrayX1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No data file found for record %d, skipping it. \n", irec);
            drms_free_array(inArrayX1);
            continue;
        }
        t1 = drms_getkey_time(inRec1, "T_OBS", NULL);	// NEED TO CHECK!!!!!
//printf("t=%f\n", t1);
//printf("dt=%f\n",t1-t0);
        dt = t1 - t0;
        if (dt < cadence) {
            printf("Cadence higher then limit, skipping current record. \n");
            drms_free_array(inArrayX1);
            continue;
        }
        Bx1 = (float *)inArrayX1->data;
        if (inArrayX1->axis[0] != nx || inArrayX1->axis[1] != ny)
            DIE("Wrong data dimension for current record, stop. \n"); 
        inSeg1 = drms_segment_lookupnum(inRec1, 1);
        inArrayY1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
        By1 = (float *)inArrayY1->data;
        inSeg1 = drms_segment_lookupnum(inRec1, 2);
        inArrayZ1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
        Bz1 = (float *)inArrayZ1->data;
        
        /* Time measure */
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", nsuccess);
        }

        /* Output data */
        outArrayV = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayDVDx = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outArrayDVDy = drms_array_create(DRMS_TYPE_FLOAT, 3, outDims, NULL, &status);
        outDataV = (float *)outArrayV->data;
        outDataDVDx = (float *)outArrayDVDx->data;
        outDataDVDy = (float *)outArrayDVDy->data;


        /* ======================== */
        /* This is the working part.*/
        /* ======================== */

        // Input
//        Bx0 = inData0; By0 = inData0 + nxny; Bz0 = inData0 + nxny * 2;
//        Bx1 = inData1; By1 = inData1 + nxny; Bz1 = inData1 + nxny * 2;
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
        Vx = outDataV; Vy = outDataV + nxny; Vz = outDataV + nxny * 2;
        DVxDx = outDataDVDx; DVyDx = outDataDVDx + nxny; DVzDx = outDataDVDx + nxny * 2;
        DVxDy = outDataDVDy; DVyDy = outDataDVDy + nxny; DVzDy = outDataDVDy + nxny * 2;
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


        // Process data
        d4vm_matrix_(a, 
                     Bx, Bxx, Bxy, By, Byx, Byy, Bz, Bzx, Bzy, Bzt,
                     k_th, k_x, k_y, k_xx, k_yy, k_xy,
                     &nx, &ny, ksize);

        // Solve matrix
        d4vm_solver_(a, 
                     Vx, Vy, Vz, 
                     DVxDx, DVyDx, DVzDx, 
                     DVxDy, DVyDy, DVzDy, 
                     &nx, &ny, &threshold);


        /* Output record, create one at a time */
        outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
        if (status) DIE("Output record not created");
        outSegV = drms_segment_lookupnum(outRec, 0);
        outSegDVDx = drms_segment_lookupnum(outRec, 1);
        outSegDVDy = drms_segment_lookupnum(outRec, 2);
        for (i = 0; i < 3; i++) {
            outSegV->axis[i] = outArrayV->axis[i];	// For variable dimensions
            outSegDVDx->axis[i] = outArrayDVDx->axis[i];
            outSegDVDy->axis[i] = outArrayDVDy->axis[i];
        }
        outArrayV->parent_segment = outSegV;
        outArrayDVDx->parent_segment = outSegDVDx;
        outArrayDVDy->parent_segment = outSegDVDy;

        /* Set keywords */
        drms_copykey(outRec, inRec0, "T_OBS");
        drms_setkey_time(outRec, "T_NEXT", t1);
        drms_setkey_int(outRec, "K0", ksize[0]);
        drms_setkey_int(outRec, "K1", ksize[1]);
        drms_setkey_double(outRec, "THRESH", threshold);

        /* Result writing */
        status = drms_segment_write(outSegV, outArrayV, 0);
        if (status) DIE("Problem writing file V");
        status = drms_segment_write(outSegDVDx, outArrayDVDx, 0);
        if (status) DIE("Problem writing file DVDx"); 
        status = drms_segment_write(outSegDVDy, outArrayDVDy, 0);
        if (status) DIE("Problem writing file DVDy"); 
        drms_close_record(outRec, DRMS_INSERT_RECORD);
        drms_free_array(outArrayV);
        drms_free_array(outArrayDVDx);
        drms_free_array(outArrayDVDy);

        /* Clean up, move to next record */
        drms_free_array(inArrayX0);
        drms_free_array(inArrayY0);
        drms_free_array(inArrayZ0);
        free(Bx); free(Bxx); free(Bxy);
        free(By); free(Byx); free(Byy);
        free(Bz); free(Bzx); free(Bzy); free(Bzt);
        free(a);
        free(k_th); free(k_x); free(k_y);
        free(k_xx); free(k_yy); free(k_xy);
        inRec0 = inRec1;
        inSeg0 = inSeg1;
        inArrayX0 = inArrayX1;
        inArrayY0 = inArrayY1;
        inArrayZ0 = inArrayZ1;
        inData0 = inData1;
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
