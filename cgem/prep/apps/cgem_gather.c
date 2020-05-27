/*
 *  cgem_gather.c
 *
 *  Example:
 %  cgem_gather "prep=hmi_test.cgem_prep[11158][2011.02.15_12:00/252m]" "vel=hmi_test.cgem_vel[11158][2011.02.15_12:00/252m]" "out=hmi_test.cgem_pdfi_in"
 %  cgem_gather "prep=hmi_test.cgem_prep[12673][2017.09.05_05:24/252m]" "vel=hmi_test.cgem_vel[12673][2017.09.05_05:24/252m]" "out=hmi_test.cgem_pdfi_in"
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "astro.h"

// Legacy macros

#define PI              (M_PI)
#define RADSINDEG        (PI/180.)
#define RAD2ARCSEC        (648000./M_PI)
#define SECINDAY        (86400.)
#define FOURK            (4096)
#define FOURK2          (16777216)
#define DTTHRESH        (518400.)                    // 6 day
#define INTERPOPT       (0)         // Interpolation

#define MAXTDIFF        (10.)       // Max allowed T_REC diff in sec

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"
#define dpath    "/home/jsoc/cvs/Development/JSOC"

// =====================================

// Input segment names in prep series and corresponding output segment names
char *prepSegNames[] = {"Bx", "By", "Bz", "Lx", "Ly", "Lz", "vlos_mag"};
char *outPrepSegNames[] = {"Bloncoe", "Blatcoe", "Brllcoe", "Lloncoe", "Llatcoe", "Lrllcoe", "Vlosllcoe"};
double scalePrepSegs[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// Input segment names in vecl series and corresponding output segment names
char *velSegNames[] = {"Vx", "Vy"};
char *outVelSegNames[] = {"Vloncoe", "Vlatcoe"};
double scaleVelSegs[] = {0.001, 0.001};

// Units for segments
char *bunits[] = {"Mx/cm^2", "Mx/cm^2", "Mx/cm^2", "km/s", "km/s", "m/s", " ", " ", " "};

// =====================================

// Pair records in two input series
int pairInputs(DRMS_RecordSet_t *prepRS, DRMS_RecordSet_t *velRS, int *pairedRecNums, int *npairs);

// Check and get CGEMNUM, T_REC, dimensions
int checkInputs(DRMS_RecordSet_t *inRS, int *cgemnums, TIME *trecs, int *nxs, int *nys);

// Copy over segment with scaling
int copyData(DRMS_Record_t *inRec, DRMS_Record_t *outRec,
             char *inSegName, char *outSegName, double scale);

// =====================================

char *module_name = "cgem_gather";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "prep", kNotSpecified, "Input data series containing field, LOS vectors and Dopplergram."},
    {ARG_STRING, "vel", kNotSpecified, "Input data series containing FLCT result."},
    {ARG_STRING, "out", kNotSpecified, "Output data series."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    char *prepQuery = (char *) params_get_str(&cmdparams, "prep");
    char *velQuery = (char *) params_get_str(&cmdparams, "vel");
    char *outQuery = (char *) params_get_str(&cmdparams, "out");
    
    /* Input data */
    
    DRMS_RecordSet_t *prepRS = drms_open_records(drms_env, prepQuery, &status);
    int nrecs = prepRS->n;      // fiducial rec number
    if (status || nrecs == 0 || !prepRS) DIE("Input prep data series error");
    
    DRMS_RecordSet_t *velRS = drms_open_records(drms_env, velQuery, &status);
    if (status || velRS->n == 0 || !velRS) DIE("Input vel data series error");
    
    int nPrepSegs = ARRLENGTH(prepSegNames);
    int nVelSegs = ARRLENGTH(velSegNames);
    
    /* Pair two input data series through CGEMNUM, T_REC, and array sizes */
    
    int npairs = 0;      // Final pairs of records
    int *pairedRecNums = (int *) (malloc(nrecs * sizeof(int))); // paired record number in velRS
    if (pairInputs(prepRS, velRS, pairedRecNums, &npairs)) {
        free(pairedRecNums);
        DIE("Pairing two input data series error or no pairs exists.");
    }
    printf("Found %d pairs for %d frames of prep data.\n", npairs, nrecs); fflush(stdout);
    
    /* Output */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, npairs, outQuery, DRMS_PERMANENT, &status);
    if (status || !outRS) {
        free(pairedRecNums);
        DIE("Error creating output series\n");
    }
    
    /* Loop */
    
    int ipair = 0;
    for (int irec = 0; irec < nrecs; irec++) {
        
        status = 0;     // reset
        if (pairedRecNums[irec] >= nrecs) continue;     // no pairing, skip frame
        
        printf("==============\nProcessing rec pairs #%d of %d, ", ipair + 1, npairs); fflush(stdout);
        
        DRMS_Record_t *prepRec = prepRS->records[irec];
        DRMS_Record_t *velRec = velRS->records[pairedRecNums[irec]];
        
        DRMS_Record_t *outRec = outRS->records[ipair];
        
        TIME trec = drms_getkey_time(prepRec, "T_REC", &status);
        char trec_str[100];
        sprint_time(trec_str, trec, "TAI", 0);
        printf("TREC=[%s].... ", trec_str);
        
        // Copy over segments in prepRS
        
        for (int iSeg = 0; iSeg < nPrepSegs; iSeg++) {
            if (copyData(prepRec, outRec,
                         prepSegNames[iSeg], outPrepSegNames[iSeg], scalePrepSegs[iSeg])) {
                printf("\nCopy segment %s error, segment skipped.\n", prepSegNames[iSeg]);
                status = 1;
                break;
            }
        }
        if (status) continue;       // skip record, T_REC=NAN
        
        // Copy over segments in velRS, velocity from m/s to km/s
        
        for (int iSeg = 0; iSeg < nVelSegs; iSeg++) {
            if (copyData(velRec, outRec,
                         velSegNames[iSeg], outVelSegNames[iSeg], scaleVelSegs[iSeg])) {
                printf("\nCopy segment %s error, segment skipped.\n", velSegNames[iSeg]);
                status = 1;
                break;
            }
        }
        if (status) continue;       // skip record, T_REC=NAN
        
        // Keys
        
        drms_copykeys(outRec, velRec, 0, kDRMS_KeyClass_Explicit);  // first, copy keywords from velRec
        drms_copykeys(outRec, prepRec, 0, kDRMS_KeyClass_Explicit); // then, overwrite most from prepRec
        
        drms_setkey_time(outRec, "DATE_PREP", drms_getkey_time(prepRec, "DATE", &status));
        drms_setkey_time(outRec, "DATE_VEL", drms_getkey_time(velRec, "DATE", &status));
        
        double crpix1 = drms_getkey_double(prepRec, "CRPIX1", &status);
        double crpix2 = drms_getkey_double(prepRec, "CRPIX2", &status);
        double crval1 = drms_getkey_double(prepRec, "CRVAL1", &status);
        double crval2 = drms_getkey_double(prepRec, "CRVAL2", &status);
        double cdelt1 = drms_getkey_double(prepRec, "CDELT1", &status);
        double cdelt2 = drms_getkey_double(prepRec, "CDELT2", &status);
        
        DRMS_Segment_t *outSeg = drms_segment_lookupnum(outRec, 0);
        int rows = outSeg->axis[1], cols = outSeg->axis[0];
        
        double minlat = crval2 + (0.5 - crpix2) * cdelt2;           // Edge is at 0.5 and rows+0.5
        double maxlat = crval2 + (rows + 0.5 - crpix2) * cdelt2;
        double minlon = crval1 + (0.5 - crpix1) * cdelt1;
        double maxlon = crval1 + (cols + 0.5 - crpix1) * cdelt1;
        
        drms_setkey_double(outRec, "MINCOLAT", (90. - maxlat) * RADSINDEG);         // a
        drms_setkey_double(outRec, "MAXCOLAT", (90. - minlat) * RADSINDEG);         // b
        drms_setkey_double(outRec, "MINLON", minlon * RADSINDEG);       // c
        drms_setkey_double(outRec, "MAXLON", maxlon * RADSINDEG);       // d
        
        int nSegs = ARRLENGTH(bunits);      // BUNIT, no checking
        char key[64];
        for (int iSeg = 0; iSeg < nSegs; iSeg++) {          // Segment units
            snprintf(key, 64, "BUNIT_%03d", iSeg);
            drms_setkey_string(outRec, key, bunits[iSeg]);
        }
        
        TIME val, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(outRec, "DATE", tnow);
        
        //
        
        SHOW("done.\n");
        
        // NEXT!!!
        
        ipair++;
        
    }   // irec
    
    /* Clean up */
    
    free(pairedRecNums);
    drms_close_records(prepRS, DRMS_FREE_RECORD);
    drms_close_records(velRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);

    return 0;
    
}   // Do it


// =====================================

/*
 * Pair records in two input series
 * Check cgemnum, t_rec, array dims
 *
 */

int pairInputs(DRMS_RecordSet_t *prepRS, DRMS_RecordSet_t *velRS, int *pairedRecNums, int *npairs)
{
    int status = 0;
    
    // Check prepRS
    
    int nrecs_p = prepRS->n;
    int *cgemnums_p = (int *) (malloc(nrecs_p * sizeof(int)));       // CGEMNUM
    TIME *trecs_p = (TIME *) (malloc(nrecs_p * sizeof(TIME)));       // T_REC
    int *nxs_p = (int *) (malloc(nrecs_p * sizeof(int)));
    int *nys_p = (int *) (malloc(nrecs_p * sizeof(int)));
    
    if (checkInputs(prepRS, cgemnums_p, trecs_p, nxs_p, nys_p)) {
        free(cgemnums_p); free(trecs_p); free(nxs_p); free(nys_p);
        printf("P\n");
        return 1;
    }
    
    // Check veclRS
    
    int nrecs_v = velRS->n;
    int *cgemnums_v = (int *) (malloc(nrecs_v * sizeof(int)));       // CGEMNUM
    TIME *trecs_v = (TIME *) (malloc(nrecs_v * sizeof(TIME)));       // T_REC
    int *nxs_v = (int *) (malloc(nrecs_v * sizeof(int)));
    int *nys_v = (int *) (malloc(nrecs_v * sizeof(int)));
    
    if (checkInputs(velRS, cgemnums_v, trecs_v, nxs_v, nys_v)) {
        free(cgemnums_p); free(trecs_p); free(nxs_p); free(nys_p);
        free(cgemnums_v); free(trecs_v); free(nxs_v); free(nys_v);
        printf("V\n");
        return 1;
    }
    
    // Pair with velRS
    
    *npairs = 0;
    for (int i = 0; i < nrecs_p; i++) {
        pairedRecNums[i] = nrecs_p;       // no findings
        for (int j = 0; j < nrecs_v; j++) {
            if (cgemnums_v[j] == cgemnums_p[i] &&
                nxs_v[j] == nxs_p[i] && nys_v[j] == nys_p[i] &&
                fabs(trecs_v[j] - trecs_p[i]) < MAXTDIFF) {
                pairedRecNums[i] = j;
                (*npairs)++;
                break;
            }
        }   // j
    }   // i
    
    // Clean up
    
    free(cgemnums_p); free(trecs_p); free(nxs_p); free(nys_p);
    free(cgemnums_v); free(trecs_v); free(nxs_v); free(nys_v);
    
    if (*npairs == 0) return 1;     // failed
    
    return 0;
}


// =====================================

/*
 * Check and get CGEMNUM, T_REC, dimensions
 *
 */


int checkInputs(DRMS_RecordSet_t *inRS, int *cgemnums, TIME *trecs, int *nxs, int *nys)
{
    int status = 0;
    
    int nrecs = inRS->n;
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        DRMS_Record_t *inRec = inRS->records[irec];
        
        cgemnums[irec] = drms_getkey_int(inRec, "CGEMNUM", &status);
        if (status) return status;
        trecs[irec] = drms_getkey_time(inRec, "T_REC", &status);
        if (status) return status;
        
        DRMS_Segment_t *inSeg = drms_segment_lookupnum(inRec, 0);
        if (status) return 1;
        nxs[irec] = inSeg->axis[0]; nys[irec] = inSeg->axis[1];
        
        int nSegs = drms_record_numsegments(inRec);
        for (int iSeg = 1; iSeg < nSegs; iSeg++) {
            inSeg = drms_segment_lookupnum(inRec, iSeg);
            if (status) return 1;
            if (nxs[irec] != inSeg->axis[0] ||
                nys[irec] != inSeg->axis[1]) {
                printf("nx=%d, ny=%d\n", inSeg->axis[0], inSeg->axis[1]);
                return 1;
            }
        }
        
    }
    
    return 0;
}


// =====================================

/*
 * Copy over segment with scaling
 *
 */

int copyData(DRMS_Record_t *inRec, DRMS_Record_t *outRec,
             char *inSegName, char *outSegName, double scale)
{
    int status = 0;
    
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, inSegName);
    DRMS_Segment_t *outSeg = drms_segment_lookup(outRec, outSegName);
    if (!inSeg || !outSeg) return 1;
    
    DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
    if (status || !inArray) return 1;
    double *inData = (double *) inArray->data;
    
    // Copy over & scale
    
    int dims[2] = {inArray->axis[0], inSeg->axis[1]};
    int nsize = dims[0] * dims[1];
    double *outData = (double *) (malloc(nsize * sizeof(double)));
    memcpy(outData, inData, nsize * sizeof(double));
    for (int i = 0; i < nsize; i++) {
        outData[i] *= scale;
    }
    drms_free_array(inArray);
    
    // Output
    
    DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims, outData, &status);
    if (status) {
        free(outData); return 1;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    drms_free_array(outArray);
    if (status) return 1;
    
    return 0;
}
