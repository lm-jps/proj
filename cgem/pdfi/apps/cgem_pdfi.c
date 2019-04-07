/*
 *  cgem_prep.c
 *
 *  Example:
 *	limit s u
 *  cgem_pdfi "in=hmi_test.cgem_pdfi_in[11158][2011.02.15_01:24/24m]" "out=hmi_test.cgem_pdfi_out"
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
#include "pdfi_vers.h"

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

//#define TESTMODE    (1)

// Adjust padding sizes such that m, n are dividable by UPADDING, must be even
#define UPADDING	(12)

#define MAXTDIFF    0.01    // relative allowed difference of cadence for pairing input frames

// Constants
// Blon equivalent to Bx, Blat equivalent to By

#ifdef TESTMODE
char *inSegNames[] = {"Bx", "By", "Bz",
                      "Vx", "Vy", "Vz",
                      "Lx", "Ly", "Lz"};
char *outSegNames[] = {"Expdfi", "Eypdfi", "Ezpdfi",
                       "dExdz", "dEydz", "Srz", "Hmz"};
char *outSegBNames[] = {"Bx0", "By0", "Bz0",
                        "Bx1", "By1", "Bz1"};
#else
char *inSegNames[] = {"Bloncoe", "Blatcoe", "Brllcoe",
                      "Vloncoe", "Vlatcoe", "Vlosllcoe",
                      "Lloncoe", "Llatcoe", "Lrllcoe"};
char *outSegNames[] = {"Elonpdfi", "Elatpdfi", "Erllpdfi", "Erllind",
                       "dElondr", "dElatdr", "Srll", "Hmll"};
char *outSegBNames[] = {"Blon0", "Blat0", "Brll0",
                        "Blon1", "Blat1", "Brll1"};
char *outSegMNames[] = {"Mcoell", "Mcoll", "Mcell", "Mtell", "Mpell"};
#endif

// Parameters

struct reqInfo {
    int n, m;       // input already has padding, n: cols; m: rows
    int npad0, mpad0;		// Initial guess of padding on each side
    int npadl, npadr, mpadb, mpadt;		// padding on left/right, bottom/top
    int n_o, m_o;   // n = npadl + n_o + npadr; m = mpadb + m_o + npadt
    double a, b, c, d;
    double bmin;
    TIME t, trec;   // use T_OBS for time, trec for prime key
};

struct pdfi_result {        // results from PDFI calculation
    double srtot, hmtot;
};

// FORTRAN function

extern void pdfi_wrapper4jsoc_ss_(int *m, int *n, double *rsun,
                                  double *a, double *b, double *c, double *d, double *bmin,
                                  double *bloncoe0, double *blatcoe0, double *brllcoe0,
                                  double *bloncoe1, double *blatcoe1, double *brllcoe1,
                                  double *vloncoe0, double *vlatcoe0, double *vlosllcoe0,
                                  double *vloncoe1, double *vlatcoe1, double *vlosllcoe1,
                                  double *lloncoe0, double *llatcoe0, double *lrllcoe0,
                                  double *lloncoe1, double *llatcoe1, double *lrllcoe1,
                                  double *tjul0, double *tjul1,
                                  double *blon0, double *blat0, double *brll0,
                                  double *blon1, double *blat1, double *brll1,
                                  double *elonpdfi, double *elatpdfi, double *erllpdfi, double *erllind,
                                  double *delondr, double *delatdr,
                                  double *srll, double *srtot, double *hmll, double *hmtot,
                                  double *mcoell, double *mcoll, double *mcell, double *mtell, double *mpell,
                                  double *tjulhalf);

// C functions

/* Find pairs of images according to T_REC and delt */
int pairTRecs(DRMS_RecordSet_t *inRS, double delt, int *pairedRecNums, int *npairs);

/* Get patch and padding dimensions */
void getRInfo(DRMS_Record_t *inRec, struct reqInfo *rInfo);

// Reading 9 input arrays
int getInputArr(DRMS_Record_t *inRec, struct reqInfo *rInfo,
                double *bloncoe, double *blatcoe, double *brllcoe,
                double *vloncoe, double *vlatcoe, double *vlosllcoe,
                double *lloncoe, double *llatcoe, double *lrllcoe);

// Writing 3 output arrays
int writeOutputArr(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, DRMS_Record_t *outRec,
                   struct reqInfo *rInfo0, struct reqInfo *rInfo1,
                   double *elonpdfi, double *elatpdfi, double *erllpdfi, double *erllind,
                   double *delondr, double *delatdr, double *srll, double *hmll,
                   double *blon0, double *blat0, double *brll0,
                   double *blon1, double *blat1, double *brll1,
                   double *mcoell, double *mcoll, double *mcell, double *mtell, double *mpell,
                   struct pdfi_result *res, int wrt_dbl);

// Figure out padding sizes
void pad_int_gen_ss(struct reqInfo *rInfo);

// =====================================

char *module_name = "cgem_pdfi";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input data series."},
    {ARG_STRING, "out", kNotSpecified, "Output data series."},
    {ARG_DOUBLE, "delt", "720.", "Temporal spacing between input frames"},
    {ARG_DOUBLE, "bmin", "250.", "Threshold of field value for masking PDFI results."},
    {ARG_INT, "npad", "50", "Initial guess of padding in column."},
    {ARG_INT, "mpad", "50", "Initial guess of padding in row."},
    {ARG_FLAG, "d", "", "Output in double instead of scaled int"},
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
    int mpad = params_get_int(&cmdparams, "mpad"); if (mpad < 0) mpad = 0;
    int npad = params_get_int(&cmdparams, "npad"); if (npad < 0) npad = 0;
    int wrt_dbl = params_isflagset(&cmdparams, "d");
    
    double delt = fabs(params_get_double(&cmdparams, "delt"));
    double bmin = fabs(params_get_double(&cmdparams, "bmin"));
    
    double rsun_ref = 6.96e5;       // fixed for HMI
    
    /* Input Data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs <= 1) {
        DIE("Input data series error");
    }
    
    // Gather all T_RECs in input, determine available pairs of output
    
    int *pairedRecNums = (int *) (malloc(nrecs * sizeof(int)));
    int npairs;
    if (pairTRecs(inRS, delt, pairedRecNums, &npairs)) {
        free(pairedRecNums);
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("No pairs of input images satisfy requested cadence.");
    }
    printf("Found %d pairs in %d frames.\n", npairs, nrecs); fflush(stdout);
    
    /* Output */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, npairs, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        free(pairedRecNums);
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("Output data series error");
    }
    
    /* Loop */
    
    int ipair = 0;
    for (int irec = 0; irec < nrecs; irec++) {
        
        if (pairedRecNums[irec] >= nrecs) continue;     // no pairing, skip frame
        
        printf("==============\nProcessing rec pairs #%d of %d\n", ipair + 1, npairs); fflush(stdout);
        
        // Input, n x m
        
        DRMS_Record_t *inRec0 = inRS->records[irec];
        DRMS_Record_t *inRec1 = inRS->records[pairedRecNums[irec]];
        
        struct reqInfo rInfo0, rInfo1;
        rInfo0.npad0 = rInfo1.npad0 = npad;
        rInfo0.mpad0 = rInfo1.mpad0 = mpad;
        rInfo0.bmin = rInfo1.bmin = bmin;
        
        // Patch & padding info
        
        getRInfo(inRec0, &rInfo0);
        getRInfo(inRec1, &rInfo1);
        if (rInfo0.n != rInfo1.n || rInfo0.m != rInfo1.m) {
            printf("Input array error, record #%d skipped", ipair + 1);
            continue;
        }
        
        // Input data arrays with padding, need calloc for zeros
        
        int n1m1 = (rInfo0.n + 1) * (rInfo0.m + 1);
        
        double *bloncoe0 = (double *) (calloc(n1m1, sizeof(double)));      // B vectors
        double *blatcoe0 = (double *) (calloc(n1m1, sizeof(double)));
        double *brllcoe0 = (double *) (calloc(n1m1, sizeof(double)));
        
        double *vloncoe0 = (double *) (calloc(n1m1, sizeof(double)));      // V vectors
        double *vlatcoe0 = (double *) (calloc(n1m1, sizeof(double)));
        double *vlosllcoe0 = (double *) (calloc(n1m1, sizeof(double)));
        
        double *lloncoe0 = (double *) (calloc(n1m1, sizeof(double)));
        double *llatcoe0 = (double *) (calloc(n1m1, sizeof(double)));      // los vectors
        double *lrllcoe0 = (double *) (calloc(n1m1, sizeof(double)));
        
        double *bloncoe1 = (double *) (calloc(n1m1, sizeof(double)));
        double *blatcoe1 = (double *) (calloc(n1m1, sizeof(double)));
        double *brllcoe1 = (double *) (calloc(n1m1, sizeof(double)));
        
        double *vloncoe1 = (double *) (calloc(n1m1, sizeof(double)));
        double *vlatcoe1 = (double *) (calloc(n1m1, sizeof(double)));
        double *vlosllcoe1 = (double *) (calloc(n1m1, sizeof(double)));
        
        double *lloncoe1 = (double *) (calloc(n1m1, sizeof(double)));
        double *llatcoe1 = (double *) (calloc(n1m1, sizeof(double)));
        double *lrllcoe1 = (double *) (calloc(n1m1, sizeof(double)));

        
        // Get data
        
        if (getInputArr(inRec0, &rInfo0,
                        bloncoe0, blatcoe0, brllcoe0,
                        vloncoe0, vlatcoe0, vlosllcoe0,
                        lloncoe0, llatcoe0, lrllcoe0)) {
            printf("Input array error, record #%d skipped", ipair);
            continue;
        }
        
        if (getInputArr(inRec1, &rInfo1,
                        bloncoe1, blatcoe1, brllcoe1,
                        vloncoe1, vlatcoe1, vlosllcoe1,
                        lloncoe1, llatcoe1, lrllcoe1)) {
            printf("Input array error, record #%d skipped", ipair);
            continue;
        }
        
        // Needed parameters
        
        int n = rInfo0.n, m = rInfo0.m;
        double a = (rInfo0.a + rInfo1.a) / 2.;      // changed to mean of two frames Jan 19 2019 XS
        double b = (rInfo0.b + rInfo1.b) / 2.;
        double c = (rInfo0.c + rInfo1.c) / 2.;
        double d = (rInfo0.d + rInfo1.d) / 2.;
        TIME t0 = rInfo0.t, t1 = rInfo1.t;
        double tjul0 = t0 / SECINDAY, tjul1 = t1 / SECINDAY, tjulhalf;
        TIME trec = (rInfo0.trec + rInfo1.trec) / 2.;
        
        char trec_str[100];
        sprint_time(trec_str, trec, "TAI", 0);
        printf("TREC=[%s], dt=%lf, djul=%lf...\n", trec_str, t1 - t0, tjul1 - tjul0);
        
        printf("n=%d, m=%d, rsun_ref=%f\n", n, m, rsun_ref);
        printf("a=%f, b=%f, c=%f, d=%f\n", a, b, c, d);
        printf("tjul0=%f, tjul1=%f\n", tjul0, tjul1);
        
        /*
        int idx = 48890;
        printf("blon0=%f, blat0=%f, brll0=%f\n", bloncoe0[idx], blatcoe0[idx], brllcoe0[idx]);
        printf("vlon0=%f, vlat0=%f, vlosll0=%f\n", vloncoe0[idx], vlatcoe0[idx], vlosllcoe0[idx]);
        printf("llon0=%f, llat0=%f, lrll0=%f\n", lloncoe0[idx], llatcoe0[idx], lrllcoe0[idx]);
        printf("blon1=%f, blat1=%f, brll1=%f\n", bloncoe1[idx], blatcoe1[idx], brllcoe1[idx]);
        printf("vlon1=%f, vlat1=%f, vlosll1=%f\n", vloncoe1[idx], vlatcoe1[idx], vlosllcoe1[idx]);
        printf("llon1=%f, llat1=%f, lrll1=%f\n", lloncoe1[idx], llatcoe1[idx], lrllcoe1[idx]);
         */
        
        // Output arrays, rec0 as default, freed in writeOutputArr
        
        double *blon0 = (double *) (calloc((n + 1) * m, sizeof(double)));
        double *blat0 = (double *) (calloc(n * (m + 1), sizeof(double)));
        double *brll0 = (double *) (calloc(n * m, sizeof(double)));
        double *blon1 = (double *) (calloc((n + 1) * m, sizeof(double)));
        double *blat1 = (double *) (calloc(n * (m + 1), sizeof(double)));
        double *brll1 = (double *) (calloc(n * m, sizeof(double)));
        
        double *elonpdfi = (double *) (calloc(n * (m + 1), sizeof(double)));
        double *elatpdfi = (double *) (calloc((n + 1) * m, sizeof(double)));
        double *erllpdfi = (double *) (calloc((n + 1) * (m + 1), sizeof(double)));
        double *erllind = (double *) (calloc((n + 1) * (m + 1), sizeof(double)));
        
        double *delondr = (double *) (calloc(n * (m + 1), sizeof(double)));
        double *delatdr = (double *) (calloc((n + 1) * m, sizeof(double)));
        
        double *srll = (double *) (calloc(n * m, sizeof(double)));
        double *hmll = (double *) (calloc(n * m, sizeof(double)));
        
        double *mcoell = (double *) (calloc((n + 1) * (m + 1), sizeof(double)));
        double *mcoll = (double *) (calloc((n - 1) * (m - 1), sizeof(double)));
        double *mcell = (double *) (calloc(n * m, sizeof(double)));
        double *mtell = (double *) (calloc(n * (m + 1), sizeof(double)));
        double *mpell = (double *) (calloc((n + 1) * m, sizeof(double)));
        
        struct pdfi_result res;
        double srtot = 0., hmtot = 0.;
        
        // Run it in FORTRAN
        // New wrapper makes sure that a, b, c, d are within bounds (Oct 21 2016)
        
        pdfi_wrapper4jsoc_ss_(&m, &n, &rsun_ref, &a, &b, &c, &d, &bmin,
                              bloncoe0, blatcoe0, brllcoe0, bloncoe1, blatcoe1, brllcoe1,
                              vloncoe0, vlatcoe0, vlosllcoe0, vloncoe1, vlatcoe1, vlosllcoe1,
                              lloncoe0, llatcoe0, lrllcoe0, lloncoe1, llatcoe1, lrllcoe1,
                              &tjul0, &tjul1,
                              blon0, blat0, brll0, blon1, blat1, brll1,
                              elonpdfi, elatpdfi, erllpdfi, erllind, delondr, delatdr,
                              srll, &srtot, hmll, &hmtot,
                              mcoell, mcoll, mcell, mtell, mpell,
                              &tjulhalf);
        res.srtot = srtot; res.hmtot = hmtot;
        
        // Output, finally save as n x m
        
        DRMS_Record_t *outRec = outRS->records[ipair];
        
        if (writeOutputArr(inRec0, inRec1, outRec,
                           &rInfo0, &rInfo1,
                           elonpdfi, elatpdfi, erllpdfi, erllind,
                           delondr, delatdr, srll, hmll,
                           blon0, blat0, brll0,
                           blon1, blat1, brll1,
                           mcoell, mcoll, mcell, mtell, mpell,
                           &res, wrt_dbl)) {
            printf("Output array error, record #%d skipped", ipair);
            free(elonpdfi); free(elatpdfi); free(erllpdfi);     // otherwise freed in function
            free(delondr); free(delatdr);
            free(srll); free(hmll);
            free(blon0); free(blat0); free(brll0);
            free(blon1); free(blat1); free(brll1);
            free(mcoell); free(mcoll); free(mcell); free(mtell); free(mpell);
        }

        // Clean up
        
        free(bloncoe0); free(blatcoe0); free(brllcoe0);
        free(bloncoe1); free(blatcoe1); free(brllcoe1);
        free(vloncoe0); free(vlatcoe0); free(vlosllcoe0);
        free(vloncoe1); free(vlatcoe1); free(vlosllcoe1);
        free(lloncoe0); free(llatcoe0); free(lrllcoe0);
        free(lloncoe1); free(llatcoe1); free(lrllcoe1);
        
        SHOW("done.\n")
        
        // NEXT!!!
        
        ipair++;
        
    }
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    drms_close_records(inRS, DRMS_FREE_RECORD);

    free(pairedRecNums);
    
    //
    
    return DRMS_SUCCESS;
    
}   // DoIt


// ========================================

/*
 * Find pairs of images according to T_REC and delt
 * For each record, find the record that closest to T_REC+delt
 * Store the record numbers in paiedRecNums[], total pairs in npairs
 *
 */

int pairTRecs(DRMS_RecordSet_t *inRS, double delt, int *pairedRecNums, int *npairs)
{
    
    int status = 0;
    
    int nrecs = inRS->n;
    int *cgemnums = (int *) (malloc(nrecs * sizeof(int)));       // CGEMNUM
    TIME *trecs = (TIME *) (malloc(nrecs * sizeof(TIME)));       // T_REC
    
    for (int irec = 0; irec < nrecs; irec++) {
        cgemnums[irec] = drms_getkey_int(inRS->records[irec], "CGEMNUM", &status);
        trecs[irec] = drms_getkey_time(inRS->records[irec], "T_REC", &status);
    }
    
    *npairs = 0;
    for (int i = 0; i < nrecs; i++) {
        pairedRecNums[i] = nrecs;       // no findings
        for (int j = 0; j < nrecs; j++) {
            if (j == i || cgemnums[j] != cgemnums[i] || trecs[j] <= trecs[i]) continue;
            if (fabs(trecs[j] - trecs[i] - delt) < (MAXTDIFF * delt)) {
                pairedRecNums[i] = j;
                (*npairs)++;
                break;
            }
        }   // j
    }   // i
    
    free(cgemnums); free(trecs);
    
    if (*npairs == 0) return 1;
    
    return 0;
    
}

// ========================================

/*
 * Get patch & padding dimensions
 *
 */
 
void getRInfo(DRMS_Record_t *inRec, struct reqInfo *rInfo)
{
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, inSegNames[0]);
    rInfo->n_o = inSeg->axis[0] - 1; rInfo->m_o = inSeg->axis[1] - 1;     // n + 1 cols, m + 1 rows
    
    pad_int_gen_ss(rInfo);        // find out the right padding size, adapted from George's code
    
    rInfo->n = rInfo->n_o + rInfo->npadl + rInfo->npadr;
    rInfo->m = rInfo->m_o + rInfo->mpadb + rInfo->mpadt;
    
//    printf("m_o=%d, m=%d, mpadb=%d, mpadt=%d\n", rInfo->m_o, rInfo->m, rInfo->mpadb, rInfo->mpadt);
//    printf("n_o=%d, n=%d, npadl=%d, npadr=%d\n", rInfo->n_o, rInfo->n, rInfo->npadl, rInfo->npadr);
}

// ========================================

/*
 * Get input data
 * To be consistent, copy over instead of carrying DRMS_Arrays over
 *
 */

int getInputArr(DRMS_Record_t *inRec, struct reqInfo *rInfo,
                double *bloncoe, double *blatcoe, double *brllcoe,
                double *vloncoe, double *vlatcoe, double *vlosllcoe,
                double *lloncoe, double *llatcoe, double *lrllcoe)
{
    int status = 0;
    
    rInfo->t = drms_getkey_time(inRec, "T_OBS", &status); if (status) return 1;
    rInfo->trec = drms_getkey_time(inRec, "T_REC", &status); if (status) return 1;
    
    // Segments
    
    DRMS_Segment_t *inSeg[9];   // 9 segments
    
    if (rInfo->n <= 0 || rInfo->m <= 0) {return 1;}
    for (int i = 0; i < 9; i++) {
        inSeg[i] = drms_segment_lookup(inRec, inSegNames[i]);
        if (inSeg[i]->axis[0] != (rInfo->n_o + 1) || inSeg[i]->axis[1] != (rInfo->m_o + 1)) {
            return 1;
        }
    }
    
    // Compute a/b/c/d
    // Assuming CRLN-CAR or HGLN-CAR
    
    double crpix1 = drms_getkey_double(inRec, "CRPIX1", &status); if (status) return 1;
    double crpix2 = drms_getkey_double(inRec, "CRPIX2", &status); if (status) return 1;
    double crval1 = drms_getkey_double(inRec, "CRVAL1", &status); if (status) return 1;
    double crval2 = drms_getkey_double(inRec, "CRVAL2", &status); if (status) return 1;
    double cdelt1 = drms_getkey_double(inRec, "CDELT1", &status); if (status) return 1;
    double cdelt2 = drms_getkey_double(inRec, "CDELT2", &status); if (status) return 1;
    
    // Fixed Apr 6 2019
    rInfo->c = (crval1 + (1.0 - rInfo->npadl - crpix1) * cdelt1) * RADSINDEG;   					// min lon
    rInfo->d = (crval1 + (rInfo->n_o + 1.0 + rInfo->npadr - crpix1) * cdelt1) * RADSINDEG;			// max lon
    rInfo->a = (90. - (crval2 + (rInfo->m_o + 1.0 + rInfo->mpadt - crpix2) * cdelt2)) * RADSINDEG;	// min co-lat
    rInfo->b = (90. - (crval2 + (1.0 - rInfo->mpadb - crpix2) * cdelt2)) * RADSINDEG;			    // max co-lat
    
//    printf("n=%d, m=%d\n", rInfo->n, rInfo->m);
//    printf("n_o=%d, m_o=%d\n", rInfo->n_o, rInfo->m_o);
//    printf("a=%36.30f, b=%36.30f\nc=%36.30f, d=%36.30f\n", rInfo->a, rInfo->b, rInfo->c, rInfo->d);
    
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
    
    // Copy & padding
    
    int col, col_o, row, row_o, lead, lead_o, idx, idx_o;
    for (row_o = 0; row_o < rInfo->m_o + 1; row_o++) {
        row = row_o + rInfo->mpadb;
        lead = row * (rInfo->n + 1);
        lead_o = row_o * (rInfo->n_o + 1);
        for (col_o = 0; col_o < rInfo->n_o + 1; col_o++) {
            col = col_o + rInfo->npadl;
            idx = lead + col;
            idx_o = lead_o + col_o;
            bloncoe[idx] = (inData[0])[idx_o];
            blatcoe[idx] = (inData[1])[idx_o];
            brllcoe[idx] = (inData[2])[idx_o];
            vloncoe[idx] = (inData[3])[idx_o];
            vlatcoe[idx] = (inData[4])[idx_o];
            vlosllcoe[idx] = (inData[5])[idx_o];
            lloncoe[idx] = (inData[6])[idx_o];
            llatcoe[idx] = (inData[7])[idx_o];
            lrllcoe[idx] = (inData[8])[idx_o];
        }
    }
    
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
                   double *elonpdfi, double *elatpdfi, double *erllpdfi, double *erllind,
                   double *delondr, double *delatdr, double *srll, double *hmll,
                   double *blon0, double *blat0, double *brll0,
                   double *blon1, double *blat1, double *brll1,
                   double *mcoell, double *mcoll, double *mcell, double *mtell, double *mpell,
                   struct pdfi_result *res, int wrt_dbl)
{
    int status = 0;
    
    // Main output
    
    DRMS_Segment_t *outSeg[8];      // added erllind Jan 28 2019
    for (int i = 0; i < 8; i++) {
         outSeg[i] = drms_segment_lookup(outRec, outSegNames[i]);
    }
    
    int n = rInfo0->n, m = rInfo0->m;
    
    double *outData[8] = {elonpdfi, elatpdfi, erllpdfi, erllind, delondr, delatdr, srll, hmll};
    DRMS_Array_t *outArray[8];
    int dims_arr[8][2] = {{n, m+1}, {n+1, m}, {n+1, m+1}, {n+1, m+1}, {n, m+1}, {n+1, m}, {n, m}, {n,m}};
    
    for (int i = 0; i < 8; i++) {
        int dims[2] = {dims_arr[i][0], dims_arr[i][1]};
        outArray[i] = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims, outData[i], &status);
        if (status) {
            for (int j = i - 1; j >=0; j--) {
                drms_free_array(outArray[j]);
            }
            return 1;
        }
        outSeg[i]->axis[0] = outArray[i]->axis[0]; outSeg[i]->axis[1] = outArray[i]->axis[1];
//        if ((i < 3) && wrt_dbl == 0) {              // Added Jan 18 2019
        if (wrt_dbl == 0) {              // Modified Jan 21 2019
            outArray[i]->israw = 0;        // compressed
            outArray[i]->bzero = outSeg[i]->bzero;
            outArray[i]->bscale = outSeg[i]->bscale;
        } else {
            outArray[i]->israw = 1;
        }
        status = drms_segment_write(outSeg[i], outArray[i], 0);
        if (status) return 1;
    }
    
    for (int i = 0; i < 8; i++) {
        drms_free_array(outArray[i]);
    }

    // Field output
    
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
        if (wrt_dbl == 0) {
            outArrayB[i]->israw = 0;        // compressed
            outArrayB[i]->bzero = outSegB[i]->bzero;
            outArrayB[i]->bscale = outSegB[i]->bscale;
        } else {
            outArrayB[i]->israw = 1;
        }
        status = drms_segment_write(outSegB[i], outArrayB[i], 0);
        if (status) return 1;
    }
    
    for (int i = 0; i < 6; i++) {
        drms_free_array(outArrayB[i]);
    }
    
    // Mask output
    
    DRMS_Segment_t *outSegM[5];
    for (int i = 0; i < 5; i++) {
        outSegM[i] = drms_segment_lookup(outRec, outSegMNames[i]);
    }
    
    double *outDataM[5] = {mcoell, mcoll, mcell, mtell, mpell};
    DRMS_Array_t *outArrayM[5];
    int dimsM_arr[5][2] = {{n+1, m+1}, {n-1, m-1}, {n, m}, {n, m+1}, {n+1, m}};
    
    for (int i = 0; i < 5; i++) {
        int dimsM[2] = {dimsM_arr[i][0], dimsM_arr[i][1]};
        outArrayM[i] = drms_array_create(DRMS_TYPE_CHAR, 2, dimsM, NULL, &status);
        if (status) {
            for (int j = i - 1; j >=0; j--) {
                drms_free_array(outArrayM[j]);
            }
            return 1;
        }
        // Copy double to char
        char *mask = (char *) outArrayM[i]->data;
        int npix = dimsM[0] * dimsM[1];
        for (int j = 0; j < npix; j++) {
            mask[j] = (outDataM[i])[j];
        }
        outSegM[i]->axis[0] = outArrayM[i]->axis[0]; outSegM[i]->axis[1] = outArrayM[i]->axis[1];
        outArrayM[i]->israw = 0;        // always compressed
        outArrayM[i]->bzero = outSegM[i]->bzero;
        outArrayM[i]->bscale = outSegM[i]->bscale;
        status = drms_segment_write(outSegM[i], outArrayM[i], 0);
        if (status) return 1;
    }
    
    for (int i = 0; i < 5; i++) {
        drms_free_array(outArrayM[i]);
        free(outDataM[i]);
    }
    
    // Keywords
    
    drms_copykeys(outRec, inRec0, 0, kDRMS_KeyClass_Explicit);     // copy all keys, for now use t0
    
    drms_setkey_string(outRec, "PDFIVERS", pdfi_vers);        // pdfi_vers.h
    
    drms_setkey_double(outRec, "SRTOT", res->srtot);            // Integrated Poynting flux
    drms_setkey_double(outRec, "HMTOT", res->hmtot);            // Helicity flux
    printf("srtot=%lf, hmtot=%lf\n", res->srtot, res->hmtot);
    
    drms_setkey_double(outRec, "BMIN", rInfo0->bmin);       // Bthreshold
    drms_setkey_int(outRec, "NPADL", rInfo0->npadl);        // padding sizes
    drms_setkey_int(outRec, "NPADR", rInfo0->npadr);
    drms_setkey_int(outRec, "MPADB", rInfo0->mpadb);
    drms_setkey_int(outRec, "MPADT", rInfo0->mpadt);
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
    
    double cdelt1 = drms_getkey_double(inRec0, "CDELT1", &status);        // assuming two records have same pixel size
    double cdelt2 = drms_getkey_double(inRec0, "CDELT2", &status);
    double crpix1 = drms_getkey_double(inRec0, "CRPIX1", &status);
    double crpix2 = drms_getkey_double(inRec0, "CRPIX2", &status);
    
    crval1_0 += ((rInfo0->npadr - rInfo0->npadl) / 2. * cdelt1);        // modify if asymmetric
    crval1_1 += ((rInfo1->npadr - rInfo1->npadl) / 2. * cdelt1);
    crval2_0 += ((rInfo0->mpadt - rInfo0->mpadb) / 2. * cdelt2);
    crval2_1 += ((rInfo1->mpadt - rInfo1->mpadb) / 2. * cdelt2);
    
    // Try wcs4pdfi functions
    /*
     double a, b, c, d;
     double crval1 = crval1_0, crval2 = crval2_0;
     int n_o = n - rInfo0->npadl - rInfo0->npadr, m_o = m - rInfo0->mpadb -rInfo0->mpadt;
     printf("crval1=%f, crval2=%f, crpix1=%f, crpix2=%f, cdelt1=%f, cdelt2=%f\n",
     crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
     wcs2pdfi(n_o, m_o, COE, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2,
     &a, &b, &c, &d);   // TE defined in wcs2pdfi.c
     printf("a=%f, b=%f, c=%f, d=%f\n", a, b, c, d);
     pdfi2wcs(n_o, m_o, COE, a, b, c, d,
     &crval1, &crval2, &crpix1, &crpix2, &cdelt1, &cdelt2);
     printf("crval1=%f, crval2=%f, crpix1=%f, crpix2=%f, cdelt1=%f, cdelt2=%f\n",
     crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
     */
    
    drms_setkey_double(outRec, "CRVAL1", (crval1_0 + crval1_1) / 2.);
    drms_setkey_double(outRec, "CRVAL2", (crval2_0 + crval2_1) / 2.);
    // Elonpdfi, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_000", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_000", (2. + m) / 2.);
    // Elatpdfi, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_001", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_001", (1. + m) / 2.);
    // Erllpdfi, COE grid, (n+1)*(m+1)
    drms_setkey_double(outRec, "CRPIX1_002", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_002", (2. + m) / 2.);
    // Erllind, COE grid, (n+1)*(m+1)
    drms_setkey_double(outRec, "CRPIX1_003", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_003", (2. + m) / 2.);
    // dElondr, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_004", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_004", (2. + m) / 2.);
    // dElatdr, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_005", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_005", (1. + m) / 2.);
    // Srll, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_006", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_006", (1. + m) / 2.);
    // Hmll, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_007", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_007", (1. + m) / 2.);
    // Blon0, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_008", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_008", (1. + m) / 2.);
    // Blat0, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_009", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_009", (2. + m) / 2.);
    // Brll0, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_010", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_010", (1. + m) / 2.);
    // Blon1, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_011", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_011", (1. + m) / 2.);
    // Blat1, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_012", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_012", (2. + m) / 2.);
    // Brll1, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_013", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_013", (1. + m) / 2.);
    // Mcoell, COE grid, (n+1)*(m+1)
    drms_setkey_double(outRec, "CRPIX1_014", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_014", (2. + m) / 2.);
    // Mcoll, CO grid, (n-1)*(m-1)
    drms_setkey_double(outRec, "CRPIX1_015", n / 2.);
    drms_setkey_double(outRec, "CRPIX2_015", m / 2.);
    // Mcell, CE grid, n*m
    drms_setkey_double(outRec, "CRPIX1_016", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_016", (1. + m) / 2.);
    // Mtell, TE grid, n*(m+1)
    drms_setkey_double(outRec, "CRPIX1_017", (1. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_017", (2. + m) / 2.);
    // Mpell, PE grid, (n+1)*m
    drms_setkey_double(outRec, "CRPIX1_018", (2. + n) / 2.);
    drms_setkey_double(outRec, "CRPIX2_018", (1. + m) / 2.);
    
    drms_setkey_string(outRec, "BUNIT_000", "V cm^(-1)");
    drms_setkey_string(outRec, "BUNIT_001", "V cm^(-1)");
    drms_setkey_string(outRec, "BUNIT_002", "V cm^(-1)");
    drms_setkey_string(outRec, "BUNIT_003", "V cm^(-1)");
    drms_setkey_string(outRec, "BUNIT_004", "V cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_005", "V cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_006", "erg cm^(-2) sec^(-1)");
    drms_setkey_string(outRec, "BUNIT_007", "Mx^2 cm^(-2) sec^(-1) ");
    drms_setkey_string(outRec, "BUNIT_008", "Mx cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_009", "Mx cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_010", "Mx cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_011", "Mx cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_012", "Mx cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_013", "Mx cm^(-2)");
    drms_setkey_string(outRec, "BUNIT_014", " ");
    drms_setkey_string(outRec, "BUNIT_015", " ");
    drms_setkey_string(outRec, "BUNIT_016", " ");
    drms_setkey_string(outRec, "BUNIT_017", " ");
    drms_setkey_string(outRec, "BUNIT_018", " ");
    
    //
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    return 0;
}

// Compute padding sizes
// Adapted from George's code

void pad_int_gen_ss(struct reqInfo *rInfo)

{
    
    // UPADDING = 12
    
	rInfo->mpadb = (((rInfo->m_o + 2 * rInfo->mpad0 + (rInfo->m_o % 2)) / UPADDING) * UPADDING - rInfo->m_o) / 2;
	rInfo->mpadt = rInfo->mpadb + (rInfo->m_o % 2);
	
	rInfo->npadl = (((rInfo->n_o + 2 * rInfo->npad0 + (rInfo->n_o % 2)) / UPADDING) * UPADDING - rInfo->n_o) / 2;
	rInfo->npadr = rInfo->npadl + (rInfo->n_o % 2);
	
}
