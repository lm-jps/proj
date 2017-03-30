/*
 *  cgem_prep.c
 *
 *  Example:
 *	limit s u
 *  cgem_pdfi "in=hmi_test.cgem_pdfi_in[11158][2011.02.15_01:24/30m]" "out=hmi_test.cgem_pdfi_out"
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "jsoc_main.h"
#include "astro.h"
#include "wcs4pdfi.c"

// Legacy macros

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"

#define TESTMODE    (1)

// Constants
// Blon equivalent to Bx, Blat equivalent to By

#ifdef TESTMODE
char *inSegNames[] = {"Bloncoe", "Blatcoe", "Brllcoe",
                      "Vloncoe", "Vlatcoe", "Vlosllcoe",
                      "Lloncoe", "Llatcoe", "Lrllcoe"};
char *outSegNames[] = {"Elonpdfi", "Elatpdfi", "Erllpdfi"};
char *outSegBNames[] = {"Blon0", "Blat0", "Brll0",
                        "Blon1", "Blat1", "Brll1"};
#else
char *inSegNames[] = {"Bx", "By", "Bz",
                      "Vx", "Vy", "Vz",
                      "Lx", "Ly", "Lz"};
char *outSegNames[] = {"Expdfi", "Eypdfi", "Ezpdfi"};
char *outSegBNames[] = {"Bx0", "By0", "Bz0",
                        "Bx1", "By1", "Bz1"};
#endif

// Parameters

struct reqInfo {
    int n, m;       // input already has padding, n: cols; m: rows
    int npad, mpad; // padding on each side
    int n_o, m_o;   // n = npad * 2 + n_o; m = mpad * 2 + m_o
    double a, b, c, d;
    TIME t, trec;   // use T_OBS for time, trec for prime key
};

// FORTRAN function

extern void pdfi_wrapper4jsoc_ss_(int *m, int *n, double *rsun, double *a, double *b, double *c, double *d,
                                  double *bloncoe0, double *blatcoe0, double *brllcoe0,
                                  double *bloncoe1, double *blatcoe1, double *brllcoe1,
                                  double *vloncoe0, double *vlatcoe0, double *vlosllcoe0,
                                  double *vloncoe1, double *vlatcoe1, double *vlosllcoe1,
                                  double *lloncoe0, double *llatcoe0, double *lrllcoe0,
                                  double *lloncoe1, double *llatcoe1, double *lrllcoe1,
                                  double *tjul0, double *tjul1,
                                  double *blon0, double *blat0, double *brll0,
                                  double *blon1, double *blat1, double *brll1,
                                  double *elonpdfi, double *elatpdfi, double *erllpdfi,
                                  double *srll, double *hmll,
                                  double *tjulhalf);

// C functions

// Reading 9 input arrays
int getInputArr(DRMS_Record_t *inRec, struct reqInfo *rInfo,
                double **bloncoe_ptr, double **blatcoe_ptr, double **brllcoe_ptr,
                double **vloncoe_ptr, double **vlatcoe_ptr, double **vloscoe_ptr,
                double **lloncoe_ptr, double **llatcoe_ptr, double **lrllcoe_ptr);

// Writing 3 output arrays
int writeOutputArr(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, DRMS_Record_t *outRec,
                   struct reqInfo *rInfo0, struct reqInfo *rInfo1,
                   double *elonpdfi, double *elatpdfi, double *erllpdfi,
                   double *blon0, double *blat0, double *brll0,
                   double *blon1, double *blat1, double *brll1);

// =====================================

char *module_name = "cgem_pdfi";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input data series."},
    {ARG_STRING, "out", kNotSpecified, "Output data series."},
    {ARG_INT, "npad", "50", "Padding in column."},
    {ARG_INT, "mpad", "50", "Padding in row."},
//    {ARG_DOUBLE, "bthr", "200.", "Threshold for masking."},
//    {ARG_DOUBLE, "bthr_r", "0.", "Threshold for relaxation."},
    {ARG_FLAG, "n", "", "No padding, overrides npad & mpad."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    char *inQuery, *outQuery;    // Data series query
    
    inQuery = (char *) params_get_str(&cmdparams, "in");
    outQuery = (char *) params_get_str(&cmdparams, "out");
    int noPadding = params_isflagset(&cmdparams, "n");
    int mpad = params_get_int(&cmdparams, "mpad");
    int npad = params_get_int(&cmdparams, "npad");
    if (noPadding) {
        mpad = npad = 0;
    }
//    double bthr = params_get_double(&cmdparams, "bthr");
//    double bthr_relax = params_get_double(&cmdparams, "bthr_r");
    
    /* Input Data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs <= 1) {
        DIE("Input data series error");
    }
    int cgemnum = drms_getkey_int(inRS->records[0], "CGEMNUM", &status);
    for (int irec = 1; irec < nrecs; irec++) {
        if (drms_getkey_int(inRS->records[irec], "CGEMNUM", &status) != cgemnum) {
            DIE("Input data series error");     // Check cgemnum
        }
    }
    
    /* Output */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs - 1, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        DIE("Output data series error");
    }
    
    /* Loop */
    
    // Useful variables
    
    double *bloncoe0, *blatcoe0, *brllcoe0;     // B vectors
    double *bloncoe1, *blatcoe1, *brllcoe1;
    
    double *vloncoe0, *vlatcoe0, *vlosllcoe0;     // V vectors
    double *vloncoe1, *vlatcoe1, *vlosllcoe1;
    
    double *lloncoe0, *llatcoe0, *lrllcoe0;     // los vectors
    double *lloncoe1, *llatcoe1, *lrllcoe1;

    double *blon0, *blat0, *brll0;      // staggered output B vector
    double *blon1, *blat1, *brll1;
    
    double *elonpdfi, *elatpdfi, *erllpdfi;     // staggered output E vector
    double *srll, *hmll;                // Poynting flux, helicity flux density
    
    int n, m;
    double a, b, c, d;
    TIME t0, t1, thalf;
    TIME tjul0, tjul1, tjulhalf;      // Julian day required by FORTRAN (with offset)
    double rsun_ref = 6.96e5;       // fixed for HMI
    
    /* Loop */
    
    for (int irec = 0; irec < nrecs - 1; irec++) {
        
        // Input, n x m
        
        DRMS_Record_t *inRec0 = inRS->records[irec];
        DRMS_Record_t *inRec1 = inRS->records[irec+1];
        
        struct reqInfo rInfo0, rInfo1;
        rInfo0.npad = rInfo1.npad = npad;
        rInfo0.mpad = rInfo1.mpad = mpad;
        
        if (getInputArr(inRec0, &rInfo0,
                        &bloncoe0, &blatcoe0, &brllcoe0,
                        &vloncoe0, &vlatcoe0, &vlosllcoe0,
                        &lloncoe0, &llatcoe0, &lrllcoe0)) {
            printf("Input array error, record #%d skipped", irec);
            continue;
        }
        
        if (getInputArr(inRec1, &rInfo1,
                        &bloncoe1, &blatcoe1, &brllcoe1,
                        &vloncoe1, &vlatcoe1, &vlosllcoe1,
                        &lloncoe1, &llatcoe1, &lrllcoe1)
            || rInfo0.m != rInfo1.m || rInfo0.n != rInfo1.n) {
            printf("Input array error, record #%d skipped", irec);
            continue;
        }
        
        // Needed parameters
        
        n = rInfo0.n; m = rInfo0.m;
        a = rInfo0.a; b = rInfo0.b; c = rInfo0.c; d = rInfo0.d;
        t0 = rInfo0.t; t1 = rInfo1.t;
        tjul0 = t0 / SECINDAY; tjul1 = t1 / SECINDAY;
        
        printf("dt=%36.30f\n", t1 - t0);
        
        /*
         a = 1.444096088409423828125000000000;
         b = 1.720951676368713378906250000000;
         c = 0.000000000000000000000000000000;
         d = 0.288267821073532104492187500000;
         t0 = 0.0;
         t1 = t0 + 1439.999994635581970214843750000000;     // for now
         */
        
        // Output arrays, n x m; rec0 as default
        
        blon0 = (double *) (malloc((n + 1) * m * sizeof(double)));
        blat0 = (double *) (malloc(n * (m + 1) * sizeof(double)));
        brll0 = (double *) (malloc(n * m * sizeof(double)));
        blon1 = (double *) (malloc((n + 1) * m * sizeof(double)));
        blat1 = (double *) (malloc(n * (m + 1) * sizeof(double)));
        brll1 = (double *) (malloc(n * m * sizeof(double)));
        
        elonpdfi = (double *) (malloc(n * (m + 1) * sizeof(double)));
        elatpdfi = (double *) (malloc((n + 1) * m * sizeof(double)));
        erllpdfi = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
        
        srll = (double *) (malloc(n * m * sizeof(double)));
        hmll = (double *) (malloc(n * m * sizeof(double)));
        
        // Run it in FORTRAN
        // New wrapper makes sure that a, b, c, d are within bounds (Oct 21 2016)
        
        pdfi_wrapper4jsoc_ss_(&m, &n, &rsun_ref, &a, &b, &c, &d,
                              bloncoe0, blatcoe0, brllcoe0, bloncoe1, blatcoe1, brllcoe1,
                              vloncoe0, vlatcoe0, vlosllcoe0, vloncoe1, vlatcoe1, vlosllcoe1,
                              lloncoe0, llatcoe0, lrllcoe0, lloncoe1, llatcoe1, lrllcoe1,
                              &tjul0, &tjul1,
                              blon0, blat0, brll0, blon1, blat1, brll1,
                              elonpdfi, elatpdfi, erllpdfi, srll, hmll,
                              &tjulhalf);
        thalf = tjulhalf * SECINDAY;
        SHOW("Fortran done\n");
        
        // Output, finally save as n x m
        
        DRMS_Record_t *outRec = outRS->records[irec];
        
        if (writeOutputArr(inRec0, inRec1, outRec,
                           &rInfo0, &rInfo1,
                           elonpdfi, elatpdfi, erllpdfi,
                           blon0, blat0, brll0,
                           blon1, blat1, brll1)) {
            printf("Output array error, record #%d skipped", irec);
            free(elonpdfi); free(elatpdfi); free(erllpdfi);     // otherwise freed in function
            free(blon0); free(blat0); free(brll0);
            free(blon1); free(blat1); free(brll1);
        }
        SHOW("here1\n");
        
        // Clean up
        
        free(bloncoe0); free(blatcoe0); free(brllcoe0);
        free(bloncoe1); free(blatcoe1); free(brllcoe1);
        free(vloncoe0); free(vlatcoe0); free(vlosllcoe0);
        free(vloncoe1); free(vlatcoe1); free(vlosllcoe1);
        free(lloncoe0); free(llatcoe0); free(lrllcoe0);
        free(lloncoe1); free(llatcoe1); free(lrllcoe1);
        
        // Not written out now
        free(srll); free(hmll);
        
    }
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    drms_close_records(inRS, DRMS_FREE_RECORD);

    //
    
    return DRMS_SUCCESS;
    
}   // DoIt




/*
 * Get input data
 * To be consistent, copy over instead of carrying DRMS_Arrays over
 *
 */

int getInputArr(DRMS_Record_t *inRec, struct reqInfo *rInfo,
                double **bloncoe_ptr, double **blatcoe_ptr, double **brllcoe_ptr,
                double **vloncoe_ptr, double **vlatcoe_ptr, double **vlosllcoe_ptr,
                double **lloncoe_ptr, double **llatcoe_ptr, double **lrllcoe_ptr)
{
    int status = 0;
    
    // Need to compute a/b/c/d from WCS
    /*
    rInfo->a = drms_getkey_double(inRec, "MINCOLAT", &status); if (status) return 1;
    rInfo->b = drms_getkey_double(inRec, "MAXCOLAT", &status); if (status) return 1;
    rInfo->c = drms_getkey_double(inRec, "MINLON", &status); if (status) return 1;
    rInfo->d = drms_getkey_double(inRec, "MAXLON", &status); if (status) return 1;
     */
    
    rInfo->t = drms_getkey_time(inRec, "T_OBS", &status); if (status) return 1;
    rInfo->trec = drms_getkey_time(inRec, "T_REC", &status); if (status) return 1;
    
    // Segments & check
    
    double *data_ptr[9];        // holder for pointers
    DRMS_Segment_t *inSeg[9];   // 9 segments
    
    inSeg[0] = drms_segment_lookup(inRec, inSegNames[0]);
    rInfo->n_o = inSeg[0]->axis[0] - 1; rInfo->m_o = inSeg[0]->axis[1] - 1;     // n + 1 cols, m + 1 rows
    
    rInfo->n = rInfo->n_o + 2 * rInfo->npad;
    rInfo->m = rInfo->m_o + 2 * rInfo->mpad;
    
    if (rInfo->n <=0 || rInfo->m <=0) {return 1;}
    for (int i = 1; i < 9; i++) {
        inSeg[i] = drms_segment_lookup(inRec, inSegNames[i]);
        if (inSeg[i]->axis[0] != (rInfo->n_o + 1) || inSeg[i]->axis[1] != (rInfo->m_o + 1)) {
            return 1;
        }
    }
    
    // Compute a/b/c/d
    // Assuming CRLN-CAR or CRLN-HG
    
    double crpix1 = drms_getkey_double(inRec, "CRPIX1", &status); if (status) return 1;
    double crpix2 = drms_getkey_double(inRec, "CRPIX2", &status); if (status) return 1;
    double crval1 = drms_getkey_double(inRec, "CRVAL1", &status); if (status) return 1;
    double crval2 = drms_getkey_double(inRec, "CRVAL2", &status); if (status) return 1;
    double cdelt1 = drms_getkey_double(inRec, "CDELT1", &status); if (status) return 1;
    double cdelt2 = drms_getkey_double(inRec, "CDELT2", &status); if (status) return 1;
    
    rInfo->c = (crval1 + (1. - rInfo->npad - crpix1) * cdelt1) * RADSINDEG;
    rInfo->d = (crval1 + (rInfo->n_o + rInfo->npad - crpix1) * cdelt1) * RADSINDEG;
    rInfo->a = (90. - (crval2 + (rInfo->m_o + rInfo->mpad - crpix2) * cdelt2)) * RADSINDEG;
    rInfo->b = (90. - (crval2 + (1. - rInfo->mpad - crpix2) * cdelt2)) * RADSINDEG;
    
    printf("n=%d, m=%d\n", rInfo->n, rInfo->m);
    printf("n_o=%d, m_o=%d\n", rInfo->n_o, rInfo->m_o);
    printf("a=%36.30f, b=%36.30f\nc=%36.30f, d=%36.30f\n", rInfo->a, rInfo->b, rInfo->c, rInfo->d);
    
    // Read segments
    
    DRMS_Array_t *inArray[9];
    double *inData[9];
    
    for (int i = 0; i < 9; i++) {
        inArray[i] = drms_segment_read(inSeg[i], DRMS_TYPE_DOUBLE, &status);
        if (status) {
            for (int j = i - 1; j >=0; j--) {
                drms_free_array(inArray[j]);
            }
            return 1;
        }
        inData[i] = (double *) inArray[i]->data;
    }
    
    // Copy over
    
    int n1m1 = (rInfo->n + 1) * (rInfo->m + 1);
    
    *bloncoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    *blatcoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    *brllcoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    
    *vloncoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    *vlatcoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    *vlosllcoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    
    *lloncoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    *llatcoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    *lrllcoe_ptr = (double *) (malloc(n1m1 * sizeof(double)));
    
    // Copy & padding
    
    int col, col_o, row, row_o, lead, lead_o, idx, idx_o;
    for (row_o = 0; row_o < rInfo->m_o + 1; row_o++) {
        row = row_o + rInfo->mpad;
        lead = row * (rInfo->n + 1);
        lead_o = row_o * (rInfo->n_o + 1);
        for (col_o = 0; col_o < rInfo->n_o + 1; col_o++) {
            col = col_o + rInfo->npad;
            idx = lead + col;
            idx_o = lead_o + col_o;
            (*bloncoe_ptr)[idx] = (inData[0])[idx_o];
            (*blatcoe_ptr)[idx] = (inData[1])[idx_o];
            (*brllcoe_ptr)[idx] = (inData[2])[idx_o];
            (*vloncoe_ptr)[idx] = (inData[3])[idx_o];
            (*vlatcoe_ptr)[idx] = (inData[4])[idx_o];
            (*vlosllcoe_ptr)[idx] = (inData[5])[idx_o];
            (*lloncoe_ptr)[idx] = (inData[6])[idx_o];
            (*llatcoe_ptr)[idx] = (inData[7])[idx_o];
            (*lrllcoe_ptr)[idx] = (inData[8])[idx_o];
        }
    }
    
    /*
    memcpy(*bloncoe_ptr, inData[0], n1m1 * sizeof(double));
    memcpy(*blatcoe_ptr, inData[1], n1m1 * sizeof(double));
    memcpy(*brllcoe_ptr, inData[2], n1m1 * sizeof(double));
    memcpy(*vloncoe_ptr, inData[3], n1m1 * sizeof(double));
    memcpy(*vlatcoe_ptr, inData[4], n1m1 * sizeof(double));
    memcpy(*vlosllcoe_ptr, inData[5], n1m1 * sizeof(double));
    memcpy(*lloncoe_ptr, inData[6], n1m1 * sizeof(double));
    memcpy(*llatcoe_ptr, inData[7], n1m1 * sizeof(double));
    memcpy(*lrllcoe_ptr, inData[8], n1m1 * sizeof(double));
     */
    
    // Clean up
    
    for (int i = 0; i < 9; i++) {
        drms_free_array(inArray[i]);
    }
        
    return 0;

}



/*
 * Write output
 * No compression for now
 * Writing staggered E
 *
 */

int writeOutputArr(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, DRMS_Record_t *outRec,
                   struct reqInfo *rInfo0, struct reqInfo *rInfo1,
                   double *elonpdfi, double *elatpdfi, double *erllpdfi,
                   double *blon0, double *blat0, double *brll0,
                   double *blon1, double *blat1, double *brll1)
{
    int status = 0;
    
    // Main output
    
    DRMS_Segment_t *outSeg[3];
    for (int i = 0; i < 3; i++) {
         outSeg[i] = drms_segment_lookup(outRec, outSegNames[i]);
    }
    
    int n = rInfo0->n, m = rInfo0->m;
    
    double *outData[3] = {elonpdfi, elatpdfi, erllpdfi};
    DRMS_Array_t *outArray[3];
    int dims_arr[3][2] = {{n, m+1}, {n+1, m}, {n+1, m+1}};
    
    for (int i = 0; i < 3; i++) {
        int dims[2] = {dims_arr[i][0], dims_arr[i][1]};
        outArray[i] = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims, outData[i], &status);
        if (status) {
            for (int j = i - 1; j >=0; j--) {
                drms_free_array(outArray[j]);
            }
            return 1;
        }
        outSeg[i]->axis[0] = outArray[i]->axis[0]; outSeg[i]->axis[1] = outArray[i]->axis[1];
        status = drms_segment_write(outSeg[i], outArray[i], 0);
    }
    
    for (int i = 0; i < 3; i++) {
        drms_free_array(outArray[i]);
    }
    
    // Keywords, T_REC/T_REC1 are both prime keys
    
    drms_copykeys(outRec, inRec0, 0, 0);     // copy all keys, for now use t0
    
    drms_setkey_int(outRec, "NPAD", rInfo0->npad);
    drms_setkey_int(outRec, "MPAD", rInfo0->mpad);
    drms_setkey_time(outRec, "T_REC0", rInfo0->trec);
    drms_setkey_time(outRec, "T_REC1", rInfo1->trec);
    drms_setkey_time(outRec, "T_REC", (rInfo0->trec + rInfo1->trec) / 2.);
    drms_setkey_time(outRec, "T_OBS0", rInfo0->t);
    drms_setkey_time(outRec, "T_OBS1", rInfo1->t);
    drms_setkey_time(outRec, "T_OBS", (rInfo0->t + rInfo1->t) / 2.);
    
    // Set a/b/c/d as mean of two input records
    
    drms_setkey_double(outRec, "MINCOLAT", (rInfo0->a + rInfo1->a) / 2.);
    drms_setkey_double(outRec, "MAXCOLAT", (rInfo0->b + rInfo1->b) / 2.);
    drms_setkey_double(outRec, "MINLON", (rInfo0->c + rInfo1->c) / 2.);
    drms_setkey_double(outRec, "MAXLON", (rInfo0->d + rInfo1->d) / 2.);
    
    // Update geometry keywords, CRPIXn are per record
    
    double crval1_0 = drms_getkey_double(inRec0, "CRVAL1", &status);
    double crval1_1 = drms_getkey_double(inRec1, "CRVAL1", &status);
    double crval2_0 = drms_getkey_double(inRec0, "CRVAL2", &status);
    double crval2_1 = drms_getkey_double(inRec1, "CRVAL2", &status);
    
    // Try wcs4pdfi functions
    double cdelt1 = drms_getkey_double(inRec0, "CDELT1", &status);
    double cdelt2 = drms_getkey_double(inRec0, "CDELT2", &status);
    double crpix1 = drms_getkey_double(inRec0, "CRPIX1", &status);
    double crpix2 = drms_getkey_double(inRec0, "CRPIX2", &status);
    double a, b, c, d;
    double crval1 = crval1_0, crval2 = crval2_0;
    int n_o = n - rInfo0->npad * 2, m_o = m - rInfo0->mpad * 2;
    printf("crval1=%f, crval2=%f, crpix1=%f, crpix2=%f, cdelt1=%f, cdelt2=%f\n",
           crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
    wcs2pdfi(n_o, m_o, COE, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2,
             &a, &b, &c, &d);   // TE defined in wcs2pdfi.c
    printf("a=%f, b=%f, c=%f, d=%f\n", a, b, c, d);
    pdfi2wcs(n_o, m_o, COE, a, b, c, d,
             &crval1, &crval2, &crpix1, &crpix2, &cdelt1, &cdelt2);
    printf("crval1=%f, crval2=%f, crpix1=%f, crpix2=%f, cdelt1=%f, cdelt2=%f\n",
           crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
    
    drms_setkey_double(outRec, "CRVAL1", (crval1_0 + crval1_1) / 2.);
    drms_setkey_double(outRec, "CRVAL2", (crval2_0 + crval2_1) / 2.);
    // Ex, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_000", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_000", (2. + m) / 2.);
    // Ey, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_001", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_001", (1. + m) / 2.);
    // Ez, COE grid, (n+1)*(m+1)
    drms_setkey_double(outRec, "CRPIX1_002", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_002", (2. + m) / 2.);
    // Bx0, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_003", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_003", (1. + m) / 2.);
    // By0, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_004", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_004", (2. + m) / 2.);
    // Bz0, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_005", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_005", (1. + m) / 2.);
    // Bx1, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_006", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_006", (1. + m) / 2.);
    // By1, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_007", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_007", (2. + m) / 2.);
    // Bz1, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_008", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_008", (1. + m) / 2.);
    
    // Supplementary output
    
    DRMS_Segment_t *outSegB[6];
    for (int i = 0; i < 6; i++) {
        outSegB[i] = drms_segment_lookup(outRec, outSegBNames[i]);
    }
    
    double *outDataB[6] = {blon0, blat0, brll0, blon1, blat1, brll1};
    DRMS_Array_t *outArrayB[6];
    int dimsB_arr[6][2] = {{n+1, m}, {n, m+1}, {n, m},
                           {n+1, m}, {n, m+1}, {n, m}};
    
    for (int i = 0; i < 6; i++) {
        int dimsB[2] = {dimsB_arr[i][0], dimsB_arr[i][1]};
        outArrayB[i] = drms_array_create(DRMS_TYPE_DOUBLE, 2, dimsB, outDataB[i], &status);
        if (status) {
            for (int j = i - 1; j >=0; j--) {
                drms_free_array(outArrayB[j]);
            }
            return 1;
        }
        outSegB[i]->axis[0] = outArrayB[i]->axis[0]; outSegB[i]->axis[1] = outArrayB[i]->axis[1];
        status = drms_segment_write(outSegB[i], outArrayB[i], 0);
    }
    
    for (int i = 0; i < 6; i++) {
        drms_free_array(outArrayB[i]);
    }

    //
    
    return 0;
}

