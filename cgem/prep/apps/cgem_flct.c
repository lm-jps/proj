/*
 *  cgem_flct.c
 *
 *	This module maps Plate-Carree Bz into Mercator map
 *  Performs FLCT, scale the velocity (how?)
 *  then interpolate vx, vy back to P-C
 *  for time step i and i+2, save to time (t(i)+t(i+2))/2.
 *
 *	Author:
 *		Xudong Sun
 *
 *	Version:
 *              v0.0 Oct 12 2016
 *	Notes:
 *		v0.0
 *
 *	Example Calls:
 *      > cgem_flct "in=hmi_test.cgem_prep[11158][2011.02.15_12:00-2011.02.15_12:24]" "out=hmi_test.cgem_vel"
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "astro.h"
#include "cartography_cgem.c"
#include "finterpolate.h"
#include <complex.h>
#include <fftw3.h>
//#include "flct4jsoc.c"
#include "flctsubs.h"           // impletmented by GHF version 1.6
#include "flct_vers.h"

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)

// 0 for Wiener, 1 for Cubic convolution, 2 for Bilinear, 3 for nnb
#define INTERPOPT       (0)

#define ARRLENGTH(ARR)  (sizeof(ARR) / sizeof(ARR[0]))

//#define TEST        1
#define NATMAP      1       // Use George's mapping code
#define MAXTDIFF    0.01    // relative allowed difference of cadence for pairing input frames

// Some other things
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"
#define dpath    "/home/jsoc/cvs/Development/JSOC"

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crval1 and crval2 are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crval1)*cosa + (wy-crval1)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crval2)*cosa - (wx-crval2)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crval1)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crval2)

// Ephemeris information
struct ephemeris {
    double crln_obs, crlt_obs;
    double rSun, asd, pa;
    double crval1, crval2, crpix1, crpix2, cdelt1, cdelt2;
    TIME t_rec, t_obs;
    char *ctype1, *ctype2;
};

// Mapping information
struct reqInfo {
    double xc, yc;
    int ncol, nrow;
    double dx, dy;
    int proj;
};

// FLCT options
struct flctOpt {
    double thresh, sigma, kr;
    int transp, absflag, doFilter, interpolate, biascor, verbose;
    int skip, poffset, qoffset, natmap;
    char vers[100];
};

/* ========================================================================================================== */

char *mapName[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
    "LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
    "stereographic", "orthographic", "LambertAzimuthal"};
char *wcsCode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
    "SIN", "ZEA"};

/* ========================================================================================================== */

/* Find pairs of images according to T_REC and delt */
int pairTRecs(DRMS_RecordSet_t *inRS, double delt, int *pairedRecNums, int *npairs);

/* Map Plate-Carree to Mercator */
int PC2Mercator(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, char *segQuery,
                int *dims_mer, double **bz0_mer_ptr, double **bz1_mer_ptr);

void pc2m(int *dims, float *bz, struct ephemeris *ephem,
          int *dims_mer, float *bz_mer);

/* Scale veclocity output */
void scaleflct(DRMS_Record_t *inRec, int *dims_mer, double *v);

/* Map Mercator back to Plate-Carree */
int Mercator2PC(DRMS_Record_t *inRec, int *dims_mer, double *vx_mer, double *vy_mer, double *vm_mer,
                int *dims, float *vx, float *vy, float *vm);

void m2pc(int *dims_mer, float *map_mer, struct ephemeris *ephem,
          int *dims, float *map);

/* Output */
int writeV(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, DRMS_Record_t *outRec,
           int *dims, struct flctOpt *fOpt, float *vx, float *vy, float *vm);

int writeSupp(DRMS_Record_t *inRec, DRMS_Record_t *outRec, int *dims_mer,
              double *vx_mer, double *vy_mer, double *vm_mer, double *bz0_mer, double *bz1_mer);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Mapping function */
void performSampling(float *inData, float *outData, float *u, float *v, int *dims_in, int *dims_out);

/* Nearest neighbor interpolation */
float nnb (float *f, int nx, int ny, float x, float y);


/* ========================================================================================================== */

char *module_name = "cgem_flct";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input series."},
    {ARG_STRING, "out", kNotSpecified, "Output series."},
    {ARG_FLOAT, "delt", "1440.", "Temporal spacing between FLCT input frames"},
    {ARG_STRING, "seg", "Bz", "Iutput segment."},
    {ARG_DOUBLE, "sigma", "5.", "For FLCT. Images modulated by gaussian of width sigma pixels."},
    {ARG_FLAG, "k", "", "For FLCT, evoke filtering below. Default no filtering"},
    {ARG_DOUBLE, "thresh", "200.", "For FLCT. Pixels below threshold ignored. In Gauss if greater than 1, relative value of max value if between 0 and 1. Default 200G"},
    {ARG_DOUBLE, "kr", "0.25", "For FLCT. Filter subimages w/gaussian w/ roll-off wavenumber kr."},
    {ARG_FLAG, "a", "", "Force threshold to be absolute value"},
    {ARG_FLAG, "q", "", "Quiet mode for FLCT"},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series and keywords */
    
    char *inQuery = NULL, *segQuery = NULL;      // Input query
    char *outQuery = NULL;                      // Output query
    
    inQuery = (char *) params_get_str(&cmdparams, "in");
    segQuery = (char *) params_get_str(&cmdparams, "seg");
    outQuery = (char *) params_get_str(&cmdparams, "out");
    
    if (!strcmp(inQuery, kNotSpecified) ||
        !strcmp(segQuery, kNotSpecified) ||
        !strcmp(outQuery, kNotSpecified)) {
        DIE("Input/output not available");
    }
    
    double delt = params_get_double(&cmdparams, "delt");
    
    // FLCT options
    struct flctOpt fOpt;
    fOpt.sigma = params_get_double(&cmdparams, "sigma");
    fOpt.doFilter = params_isflagset(&cmdparams, "k");
    fOpt.thresh = params_get_double(&cmdparams, "thresh");
    fOpt.kr = params_get_double(&cmdparams, "kr");
    fOpt.verbose = !(params_isflagset(&cmdparams, "q"));
    fOpt.absflag = params_isflagset(&cmdparams, "a");
    // Some additional default options for FLCT
    fOpt.transp = 1;
    fOpt.skip = 0;
    fOpt.poffset = 0;
    fOpt.qoffset = 0;
    fOpt.interpolate = 0;
    fOpt.biascor = 0;
    strcpy(fOpt.vers, flct_vers);       // version
#if NATMAP
    fOpt.natmap = 1;
#else
    fOpt.natmap = 0;
#endif
    
    /* Input Data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs <= 2 || !inRS) DIE("Input data series error");
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRS->records[0], segQuery);        // check segment
    if (!inSeg) {
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("Requested segment not available");
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

    /*
    int ii = 0;
    for (int irec = 0; irec < nrecs; irec++) {
        if (pairedRecNums[irec] == nrecs) continue;
        printf("pair %d: #%d & #%d, ", ii + 1, irec, pairedRecNums[irec]);
        TIME tobs0 = drms_getkey_time(inRS->records[irec], "T_OBS", &status);
        TIME tobs1 = drms_getkey_time(inRS->records[pairedRecNums[irec]], "T_OBS", &status);
        TIME trec0 = drms_getkey_time(inRS->records[irec], "T_REC", &status);
        TIME trec1 = drms_getkey_time(inRS->records[pairedRecNums[irec]], "T_REC", &status);
        char trec_str[100];
        sprint_time(trec_str, (trec0 + trec1) / 2., "TAI", 0);
        printf("T_REC=[%s], dt=%f\n", trec_str, tobs1 - tobs0);
        ii++;
    }
     */
    
    /* Output Data */

    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, npairs, outQuery, DRMS_PERMANENT, &status);
    if (status || !outRS) {
        free(pairedRecNums);
        DIE("Error creating output series\n");
    }
    
    /* Loop */
    
    int ipair = 0;
    for (int irec = 0; irec < nrecs; irec++) {
        
        if (pairedRecNums[irec] >= nrecs) continue;     // no pairing, skip frame
        
        DRMS_Record_t *inRec0 = inRS->records[irec];
        DRMS_Record_t *inRec1 = inRS->records[pairedRecNums[irec]];
        
        DRMS_Record_t *outRec = outRS->records[ipair];
        printf("==============\nProcessing rec pairs #%d of %d, ", ipair + 1, npairs); fflush(stdout);
        
        // Time info and pixel size
        
        TIME tobs0 = drms_getkey_time(inRec0, "T_OBS", &status);
        TIME tobs1 = drms_getkey_time(inRec1, "T_OBS", &status);
        TIME tobs = (tobs0 + tobs1) / 2.;
        
        TIME trec0 = drms_getkey_time(inRec0, "T_REC", &status);
        TIME trec1 = drms_getkey_time(inRec1, "T_REC", &status);
        TIME trec = (trec0 + trec1) / 2.;
        
        double dt = tobs1 - tobs0;
        char trec_str[100];
        sprint_time(trec_str, trec, "TAI", 0);
        printf("TREC=[%s], dt=%f... \n", trec_str, dt);
        
        double rSun_ref = drms_getkey_double(inRec0, "RSUN_REF", &status);        // meter, assuming two frames identical
        double cdelt = drms_getkey_double(inRec0, "CDELT1", &status);
        double ds = cdelt * RADSINDEG * rSun_ref;
//        printf("cd=%f, ds=%f, dt=%f\n", cdelt, ds, dt);
        
        // Check bz size
        
        DRMS_Segment_t *inSeg0 = NULL, *inSeg1 = NULL;
        inSeg0 = drms_segment_lookup(inRec0, segQuery);
        inSeg1 = drms_segment_lookup(inRec1, segQuery);
        if (!inSeg0 || !inSeg1) {
            SHOW("\nSegment check error, frames skipped.\n");
            continue;
        }
        if ((inSeg0->axis[0] != inSeg1->axis[0]) ||
            (inSeg0->axis[1] != inSeg1->axis[1])) {
            SHOW("\nImage dimensions do not match, frames skipped.\n");
            continue;
        }
        int nx = inSeg0->axis[0], ny = inSeg0->axis[1], nxny = nx * ny;
        int dims[2] = {nx, ny};
        
        // Output data
        
        float *vx = (float *) (malloc(nxny * sizeof(float)));
        float *vy = (float *) (malloc(nxny * sizeof(float)));
        float *vm = (float *) (malloc(nxny * sizeof(float)));

#if NATMAP          // Use George's mapping
        
        // No error checking
        
        DRMS_Array_t *inArray0 = drms_segment_read(inSeg0, DRMS_TYPE_DOUBLE, &status);
        if (status) {
            SHOW("\nArray read error... ");
            return 1;
        }
        double *bz0 = (double *)inArray0->data;

        DRMS_Array_t *inArray1 = drms_segment_read(inSeg1, DRMS_TYPE_DOUBLE, &status);
        if (status) {
            SHOW("\nArray read error... ");
            return 1;
        }
        double *bz1 = (double *)inArray1->data;
        
        struct ephemeris ephem0, ephem1;
        getEphemeris(inRec0, &ephem0);
        getEphemeris(inRec1, &ephem1);
        
        double latmin0 = (ephem0.crval2 - ny / 2. * ephem0.cdelt2) * RADSINDEG;
        double latmin1 = (ephem1.crval2 - ny / 2. * ephem1.cdelt2) * RADSINDEG;
        double latmax0 = (ephem0.crval2 + ny / 2. * ephem0.cdelt2) * RADSINDEG;
        double latmax1 = (ephem1.crval2 + ny / 2. * ephem1.cdelt2) * RADSINDEG;
        double latmin = (latmin0 + latmin1) / 2.;
        double latmax = (latmax0 + latmax1) / 2.;
        
//        printf("crlva2_0=%lf, cdelt2_0=%lf\n", ephem0.crval2, ephem0.cdelt2);
//        printf("latmin=%lf, latmax=%lf\n", latmin, latmax);
        
        double *vx_d = (double *) (malloc(nxny * sizeof(double)));      // need double
        double *vy_d = (double *) (malloc(nxny * sizeof(double)));
        double *vm_d = (double *) (malloc(nxny * sizeof(double)));
        
        if (flct_pc(fOpt.transp, bz0, bz1, dims[0], dims[1], dt, ds, fOpt.sigma,
                    vx_d, vy_d, vm_d, fOpt.thresh, fOpt.absflag, fOpt.doFilter, fOpt.kr, fOpt.skip, fOpt.poffset, fOpt.qoffset,
                    fOpt.interpolate, latmin, latmax, fOpt.biascor, fOpt.verbose)) {
            drms_free_array(inArray0); drms_free_array(inArray1);
            SHOW("FLCT error, frames skipped.\n");
            continue;
        }
        
        for (int i = 0; i < nxny; i++) {
            vx[i] = vx_d[i];
            vy[i] = vy_d[i];
            vm[i] = vm_d[i];
        }
        
        free(vx_d); free(vy_d); free(vm_d);
        drms_free_array(inArray0); drms_free_array(inArray1);
        
#else               // NATMAP
        
        // Map to Mercator
        
        int dims_mer[2];       // P-C map sizes (read); mercator map sizes (TBD)
        double *bz0_mer = NULL, *bz1_mer = NULL;
        
        if (PC2Mercator(inRec0, inRec1, segQuery, dims_mer, &bz0_mer, &bz1_mer)) {
            if (bz0_mer) free(bz0_mer); if (bz1_mer) free(bz1_mer);
            SHOW("Carree to Mercator mapping error, frames skipped.\n");
            continue;
        }
        
        // FLCT
        
        int nx_mer = dims_mer[0], ny_mer = dims_mer[1], nxny_mer = nx_mer * ny_mer;
        double *vx_mer = (double *) (malloc(nxny_mer * sizeof(double)));
        double *vy_mer = (double *) (malloc(nxny_mer * sizeof(double)));
        double *vm_mer = (double *) (malloc(nxny_mer * sizeof(double)));
        
        // Threshold performed on Mercator map
        // when mapped back to P-C, zeros (masked) reduces velocity
        // This is taken care in George's native mapping method but not here
        
        if (flct(fOpt.transp, bz0_mer, bz1_mer, dims_mer[0], dims_mer[1], dt, ds, fOpt.sigma,
                 vx_mer, vy_mer, vm_mer, fOpt.thresh, fOpt.absflag, fOpt.doFilter, fOpt.kr, fOpt.skip, fOpt.poffset, fOpt.qoffset,
                 fOpt.interpolate, fOpt.biascor, fOpt.verbose)) {
            if (bz0_mer) free(bz0_mer); if (bz1_mer) free(bz1_mer);
            free(vx_mer); free(vy_mer); free(vm_mer);
            SHOW("FLCT error, frames skipped.\n");
            continue;
        }
        
        // Map back to P-C
        // Need to do scaling: values scaled by cos(lat) * pix / sec

        scaleflct(inRec_target, dims_mer, vx_mer);
        scaleflct(inRec_target, dims_mer, vy_mer);

        if (Mercator2PC(inRec1, dims_mer, vx_mer, vy_mer, vm_mer, dims, vx, vy, vm)) {
            if (bz0_mer) free(bz0_mer); if (bz1_mer) free(bz1_mer);
            free(vx_mer); free(vy_mer); free(vm_mer); free(vx); free(vy); free(vm);
            SHOW("Mercator to Carree mapping error, frames skipped.\n");
            continue;
        }
        
        // Rewrite vmask. If interpolated set to 0.5
        /*
        for (int i = 0; i < nxny; i++) {
            if (vm[i] > 1.e-6) {
                vm[i] = 0.5;
            } else {
                vm[i] = 0;
            }
        }
         */

#ifdef TEST
        if (writeSupp(inRec1, outRec, dims_mer,
                      vx_mer, vy_mer, vm_mer, bz0_mer, bz1_mer)) {
            SHOW("Output intermediate data error, carrying on.\n");
        }
#else
        if (bz0_mer) free(bz0_mer); if (bz1_mer) free(bz1_mer);
        if (vx_mer) free(vx_mer); if (vy_mer) free(vy_mer); if (vz_mer) free(vm_mer);
#endif
        
#endif          // NATMAP
        
        // Output, vx, vy, vm freed inside writeV

        if (writeV(inRec0, inRec1, outRec, dims, &fOpt, vx, vy, vm)) {
            SHOW("Output error, frames skipped.\n");
            continue;
        }
 
        SHOW("done.\n")
        
        // NEXT!!!
        
        ipair++;

    }       // irec
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    free(pairedRecNums);
    
    return DRMS_SUCCESS;
    
}	// DoIt


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


// =============================================

/*
 * Plate-Carree to Mercator
 *
 */

int PC2Mercator(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, char *segQuery,
                int *dims_mer, double **bz0_mer_ptr, double **bz1_mer_ptr)
{
    int status = 0;
    
    // Ephemeris
    
    struct ephemeris ephem0, ephem1;
    if (getEphemeris(inRec0, &ephem0) || getEphemeris(inRec1, &ephem1)) {
        SHOW("\nEphemeris error... ");
        return 1;
    }
    if (strstr(ephem0.ctype1, "CAR") == NULL ||
        strstr(ephem0.ctype2, "CAR") == NULL ||
        strstr(ephem1.ctype1, "CAR") == NULL ||
        strstr(ephem1.ctype2, "CAR") == NULL) {
        SHOW("\nMust be in Plate Carree ");
        return 1;
    }
    
    // Carree map size
    
    DRMS_Segment_t *inSeg0 = NULL, *inSeg1 = NULL;
    inSeg0 = drms_segment_lookup(inRec0, segQuery);
    inSeg1 = drms_segment_lookup(inRec1, segQuery);
    int nx = inSeg0->axis[0], ny = inSeg0->axis[1];
    int nxny = nx * ny, dims[2] = {nx, ny};

    // Determin geometry using frame 0
    // Mercater map size using Brian's formula
    
    double yc_mer = ephem0.crval2 * RADSINDEG;
    double dy_mer = ephem0.cdelt2 * RADSINDEG;
    int nx_mer = nx;
    double th_i = yc_mer - dy_mer * ny / 2.;
    double th_f = yc_mer + dy_mer * ny / 2.;
    int ny_mer = ceil(log((tan(th_f) + 1. / cos(th_f)) / (tan(th_i) + 1. / cos(th_i))) / dy_mer);
    int nxny_mer = nx_mer * ny_mer;
    dims_mer[0] = nx_mer; dims_mer[1] = ny_mer;
//    printf("nx=%d, ny=%d, nx_mer=%d, ny_mer=%d\n", nx, ny, nx_mer, ny_mer);
    
    // Frame 0
    
    DRMS_Array_t *inArray0 = drms_segment_read(inSeg0, DRMS_TYPE_FLOAT, &status);
    if (status) {
        SHOW("\nArray read error... ");
        return 1;
    }
    float *bz0 = (float *)inArray0->data;
    float *bz0_mer = (float *) (malloc(nxny_mer * sizeof(float)));
    
    pc2m(dims, bz0, &ephem0, dims_mer, bz0_mer);
    
    *bz0_mer_ptr = (double *) (malloc(nxny_mer * sizeof(double)));      // copy float to double
    for (int i = 0; i < nxny_mer; i++) {
        (*bz0_mer_ptr)[i] = bz0_mer[i];
    }
    drms_free_array(inArray0); free(bz0_mer);

    // Frame 1
    
    DRMS_Array_t *inArray1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
    if (status) {
        SHOW("\nArray read error... ");
        free(*bz0_mer_ptr);
        return 1;
    }
    float *bz1 = (float *)inArray1->data;
    float *bz1_mer = (float *) (malloc(nxny_mer * sizeof(float)));
    
    pc2m(dims, bz1, &ephem1, dims_mer, bz1_mer);
    
    *bz1_mer_ptr = (double *) (malloc(nxny_mer * sizeof(double)));      // copy float to double
    for (int i = 0; i < nxny_mer; i++) {
        (*bz1_mer_ptr)[i] = bz1_mer[i];
    }
    drms_free_array(inArray1); free(bz1_mer);
    
    //
    
    return 0;
}

// =============================================

/*
 * Perform mapping from Plate Carre to Mercator
 *
 */

void pc2m(int *dims, float *bz, struct ephemeris *ephem,
          int *dims_mer, float *bz_mer)
{

    int nx = dims[0], ny = dims[1], nxny = nx * ny;
    int nx_mer = dims_mer[0], ny_mer = dims_mer[1], nxny_mer = nx_mer * ny_mer;
    
    int isCarr = (strstr(ephem->ctype1, "CRLN") != NULL) ? 1 : 0;
    //    printf("isCarr=%d\n", isCarr);
    double xc = ((isCarr) ? (ephem->crval1 - ephem->crln_obs) : (ephem->crval1)) * RADSINDEG;
    double yc = ephem->crval2 * RADSINDEG;
    double dx = ephem->cdelt1 * RADSINDEG;
    double dy = ephem->cdelt2 * RADSINDEG;
    double colc = ephem->crpix1 - 1;
    double rowc = ephem->crpix1 - 1;
    //  printf("xc=%f, yc=%f, dx=%f, dy=%f\n", xc, yc, dx, dy);
    
    double xc_mer = xc, yc_mer = yc, dx_mer = dx, dy_mer = dy;      // all identical
    double rowc_mer = (ny_mer - 1.) / 2., colc_mer = (nx_mer - 1.) / 2.;
    
    // Coordinates
    
    float *col_work = (float *) (malloc(nxny_mer * sizeof(float)));
    float *row_work = (float *) (malloc(nxny_mer * sizeof(float)));
    
    // Mapping
    
    double x, y, x_mer, y_mer, lat, lon;
    int idx;
    for (int row_mer = 0; row_mer < ny_mer; row_mer++) {
        y_mer = (row_mer - rowc_mer) * dy_mer;
        for (int col_mer = 0; col_mer < nx_mer; col_mer++) {
            x_mer = (col_mer - colc_mer) * dx_mer;
            idx = row_mer * nx_mer + col_mer;
            if (plane2sphere(x_mer, y_mer, yc_mer, xc_mer, &lat, &lon, MERCATOR)) {     // requested
                col_work[idx] = DRMS_MISSING_FLOAT;
                row_work[idx] = DRMS_MISSING_FLOAT;
                continue;
            }
            if (sphere2plane(lat, lon, yc, xc, &x, &y, RECTANGULAR)) {       // orig
                col_work[idx] = DRMS_MISSING_FLOAT;
                row_work[idx] = DRMS_MISSING_FLOAT;
                continue;
            }
            col_work[idx] = x / dx + colc;      // index in input images
            row_work[idx] = y / dy + rowc;
        }
    }
    
    // Interpolate
    
    performSampling(bz, bz_mer, col_work, row_work, dims, dims_mer);
    
    // Change NaNs back to 0
    
    for (idx = 0; idx < nxny_mer; idx++) {
        if (isnan(bz_mer[idx])) bz_mer[idx] = 0.;
    }

    free(col_work); free(row_work);
    
}


// =============================================

/*
 * Scale veclocity output by cos(lat) * dx / ds
 * in m/s
 *
 */

void scaleflct(DRMS_Record_t *inRec, int *dims_mer, double *v)
{
 
    int status = 0;
    
    struct ephemeris ephem;
    status = getEphemeris(inRec, &ephem);
    
    // Determine lat, then scale
    
    int nx_mer = dims_mer[0], ny_mer = dims_mer[1], nxny_mer = nx_mer * ny_mer;
    double rowc_mer = (ny_mer - 1.) / 2., colc_mer = (nx_mer - 1.) / 2.;
    double dx_mer = ephem.cdelt1 * RADSINDEG;
    double dy_mer = ephem.cdelt2 * RADSINDEG;
    int isCarr = (strstr(ephem.ctype1, "CRLN") != NULL) ? 1 : 0;
    double xc_mer = ((isCarr) ? (ephem.crval1 - ephem.crln_obs) : (ephem.crval1)) * RADSINDEG;
    double yc_mer = ephem.crval2 * RADSINDEG;
    
    double x_mer, y_mer, lat, lon;
    int idx;
    for (int row_mer = 0; row_mer < ny_mer; row_mer++) {
        y_mer = (row_mer - rowc_mer) * dy_mer;
        for (int col_mer = 0; col_mer < nx_mer; col_mer++) {
            x_mer = (col_mer - colc_mer) * dx_mer;
            idx = row_mer * nx_mer + col_mer;
            if (plane2sphere(x_mer, y_mer, yc_mer, xc_mer, &lat, &lon, MERCATOR)) {     // requested
                v[idx] = DRMS_MISSING_DOUBLE;
                continue;
            }
            v[idx] *= cos(lat);
        }
    }

}


// =============================================

/*
 * Map Mercator back to Plate Carree
 *
 */

int Mercator2PC(DRMS_Record_t *inRec, int *dims_mer, double *vx_mer, double *vy_mer, double *vm_mer,
                int *dims, float *vx, float *vy, float *vm)
{
    int status = 0;
    
    // Ephem
    
    struct ephemeris ephem;
    if (getEphemeris(inRec, &ephem)) {
        SHOW("\nEphemeris error... ");
        return 1;
    }
    
    // Map
    
    int nx_mer = dims_mer[0], ny_mer = dims_mer[1], nxny_mer = nx_mer * ny_mer;
    
    float *vx_mer_t = (float *) (malloc(nxny_mer * sizeof(float)));      // copy double to float
    float *vy_mer_t = (float *) (malloc(nxny_mer * sizeof(float)));
    float *vm_mer_t = (float *) (malloc(nxny_mer * sizeof(float)));
    for (int i = 0; i < nxny_mer; i++) {
        vx_mer_t[i] = vx_mer[i];
        vy_mer_t[i] = vy_mer[i];
        vm_mer_t[i] = vm_mer[i];
    }
    
    m2pc(dims_mer, vx_mer_t, &ephem, dims, vx);
    m2pc(dims_mer, vy_mer_t, &ephem, dims, vy);
    m2pc(dims_mer, vm_mer_t, &ephem, dims, vm);
    
    free(vx_mer_t); free(vy_mer_t); free(vm_mer_t);
    
    return 0;
}

// =============================================

/*
 * Perform mapping from Mercator to Plate Carre
 *
 */

void m2pc(int *dims_mer, float *map_mer, struct ephemeris *ephem,
          int *dims, float *map)
{
    
    int nx = dims[0], ny = dims[1], nxny = nx * ny;
    int nx_mer = dims_mer[0], ny_mer = dims_mer[1], nxny_mer = nx_mer * ny_mer;

    int isCarr = (strstr(ephem->ctype1, "CRLN") != NULL) ? 1 : 0;
    //    printf("isCarr=%d\n", isCarr);
    double xc = ((isCarr) ? (ephem->crval1 - ephem->crln_obs) : (ephem->crval1)) * RADSINDEG;
    double yc = ephem->crval2 * RADSINDEG;
    double dx = ephem->cdelt1 * RADSINDEG;
    double dy = ephem->cdelt2 * RADSINDEG;
    double colc = ephem->crpix1 - 1;
    double rowc = ephem->crpix2 - 1;
    
    double xc_mer = xc, yc_mer = yc, dx_mer = dx, dy_mer = dy;      // all identical
    double rowc_mer = (ny_mer - 1.) / 2., colc_mer = (nx_mer - 1.) / 2.;
    
    // Coordinates
    
    float *col_work = (float *) (malloc(nxny * sizeof(float)));
    float *row_work = (float *) (malloc(nxny * sizeof(float)));
    
    // Mapping
    
    double x, y, x_mer, y_mer, lat, lon;
    int idx;
    for (int row = 0; row < ny; row++) {
        y = (row - rowc) * dy;
        for (int col = 0; col < nx; col++) {
            x = (col - colc) * dx;
            idx = row * nx + col;
            if (plane2sphere(x, y, yc, xc, &lat, &lon, RECTANGULAR)) {     // requested
                col_work[idx] = DRMS_MISSING_FLOAT;
                row_work[idx] = DRMS_MISSING_FLOAT;
                continue;
            }
            if (sphere2plane(lat, lon, yc_mer, xc_mer, &x_mer, &y_mer, MERCATOR)) {       // orig
                col_work[idx] = DRMS_MISSING_FLOAT;
                row_work[idx] = DRMS_MISSING_FLOAT;
                continue;
            }
            col_work[idx] = x_mer / dx_mer + colc_mer;      // index in input images
            row_work[idx] = y_mer / dy_mer + rowc_mer;
        }
    }
    
    // Interpolate
    
    performSampling(map_mer, map, col_work, row_work, dims_mer, dims);
    
    // Change NaNs back to 0
    
    for (idx = 0; idx < nxny; idx++) {
        if (isnan(map[idx])) map[idx] = 0.;
    }
    
    free(col_work); free(row_work);

}


// =============================================

/*
 * Output
 *
 */

int writeV(DRMS_Record_t *inRec0, DRMS_Record_t *inRec1, DRMS_Record_t *outRec,
           int *dims, struct flctOpt *fOpt, float *vx, float *vy, float *vm)
{
    int status = 0;
    
    DRMS_Segment_t *outSeg;
    DRMS_Array_t *outArray;
    
    int nx = dims[0], ny = dims[1], nxny = nx * ny;
    
    // Vx
    
    outSeg = drms_segment_lookup(outRec, "Vx");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, vx, &status);
    if (status) {
        SHOW("\nOutput vx error... ");
        free(vx); free(vy); free(vm);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput vx write error... ");
        free(vx); free(vy); free(vm);
        return status;
    }
    drms_free_array(outArray);
    
    // Vy

    outSeg = drms_segment_lookup(outRec, "Vy");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, vy, &status);
    if (status) {
        SHOW("\nOutput vy error... ");
        free(vy); free(vm);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput vy write error... ");
        free(vy); free(vm);
        return status;
    }
    drms_free_array(outArray);
    
    // Vm
    
    outSeg = drms_segment_lookup(outRec, "Vm");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, vm, &status);
    if (status) {
        SHOW("\nOutput vm error... ");
        free(vm);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput vm write error... ");
        free(vm);
        return status;
    }
    drms_free_array(outArray);
    
    // Keywords
    
    drms_copykeys(outRec, inRec0, 0, kDRMS_KeyClass_Explicit);     // most keywords copied from inRec0
    
    TIME tobs0 = drms_getkey_time(inRec0, "T_OBS", &status);
    TIME tobs1 = drms_getkey_time(inRec1, "T_OBS", &status);
    TIME tobs = (tobs0 + tobs1) / 2.;
    drms_setkey_time(outRec, "T_OBS", tobs);
    
    TIME trec0 = drms_getkey_time(inRec0, "T_REC", &status);
    TIME trec1 = drms_getkey_time(inRec1, "T_REC", &status);
    TIME trec = (trec0 + trec1) / 2.;
    drms_setkey_time(outRec, "T_REC", trec);
    drms_setkey_time(outRec, "T_REC0", trec0);
    drms_setkey_time(outRec, "T_REC1", trec1);
    
    TIME val, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    // FLCT options
    
    drms_setkey_double(outRec, "SIGMA", fOpt->sigma);
    drms_setkey_double(outRec, "BTHRESH", fOpt->thresh);
    if (fOpt->absflag || fOpt->thresh > 1.0) {
        drms_setkey_int(outRec, "ABSTHR", 1);
    } else {
        drms_setkey_int(outRec, "ABSTHR", 0);
    }
    if (fOpt->doFilter) {
        drms_setkey_double(outRec, "KR", fOpt->kr);
    } else {
        drms_setkey_double(outRec, "KR", 0.0);
    }
    drms_setkey_int(outRec, "BIASCORR", fOpt->biascor);
    drms_setkey_int(outRec, "NATMAP", fOpt->natmap);
    drms_setkey_string(outRec, "FLCTVERS", fOpt->vers);
    
    drms_setkey_string(outRec, "BUNIT_000", "m/s");
    drms_setkey_string(outRec, "BUNIT_001", "m/s");
    drms_setkey_string(outRec, "BUNIT_002", " ");
    
    return 0;
}

// =============================================

/*
 * Output supplementary data
 *
 */

int writeSupp(DRMS_Record_t *inRec, DRMS_Record_t *outRec, int *dims_mer,
              double *vx_mer, double *vy_mer, double *vm_mer, double *bz0_mer, double *bz1_mer)
{
    int status = 0;
    
    DRMS_Segment_t *outSeg;
    DRMS_Array_t *outArray;
    
    int nx_mer = dims_mer[0], ny_mer = dims_mer[1], nxny_mer = nx_mer * ny_mer;
    
    // Vx_mer
    
    outSeg = drms_segment_lookup(outRec, "Vx_mer");
    outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims_mer, vx_mer, &status);
    if (status) {
        SHOW("\nOutput vx_mer error... ");
        free(vx_mer); free(vy_mer); free(vm_mer); free(bz0_mer); free(bz1_mer);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput vx_mer write error... ");
        free(vx_mer); free(vy_mer); free(vm_mer); free(bz0_mer); free(bz1_mer);
        return status;
    }
    drms_free_array(outArray);
    
    // Vy_mer
    
    outSeg = drms_segment_lookup(outRec, "Vy_mer");
    outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims_mer, vy_mer, &status);
    if (status) {
        SHOW("\nOutput vy_mer error... ");
        free(vy_mer); free(vm_mer); free(bz0_mer); free(bz1_mer);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput vy_mer write error... ");
        free(vy_mer); free(vm_mer); free(bz0_mer); free(bz1_mer);
        return status;
    }
    drms_free_array(outArray);
    
    // Vm_mer
    
    outSeg = drms_segment_lookup(outRec, "Vm_mer");
    outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims_mer, vm_mer, &status);
    if (status) {
        SHOW("\nOutput vm_mer error... ");
        free(vm_mer); free(bz0_mer); free(bz1_mer);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput vm_mer write error... ");
        free(vm_mer); free(bz0_mer); free(bz1_mer);
        return status;
    }
    drms_free_array(outArray);

    // Bz0_mer
    
    outSeg = drms_segment_lookup(outRec, "Bz0_mer");
    outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims_mer, bz0_mer, &status);
    if (status) {
        SHOW("\nOutput bz0_mer error... ");
        free(bz0_mer); free(bz1_mer);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput bz0_mer write error... ");
        free(bz0_mer); free(bz1_mer);
        return status;
    }
    drms_free_array(outArray);

    // Bz1_mer
    
    outSeg = drms_segment_lookup(outRec, "Bz1_mer");
    outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims_mer, bz1_mer, &status);
    if (status) {
        SHOW("\nOutput bz1_mer error... ");
        free(bz1_mer);
        return status;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("\nOutput bz1_mer write error... ");
        free(bz1_mer);
        return status;
    }
    drms_free_array(outArray);
    
    return 0;
}

// =============================================

/*
 * Fetch ephemeris info from a DRMS record
 * No error checking for now
 *
 */

int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem)
{
    
    int status = 0;
    
    double crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation
    double sina = sin(crota2 * RADSINDEG);
    double cosa = cos(crota2 * RADSINDEG);
    
    ephem->pa = - crota2 * RADSINDEG;
    ephem->crlt_obs = drms_getkey_float(inRec, "CRLT_OBS", &status);
    ephem->crln_obs = drms_getkey_float(inRec, "CRLN_OBS", &status);
    
    ephem->crval1 = drms_getkey_double(inRec, "CRVAL1", &status);
    ephem->crval2 = drms_getkey_double(inRec, "CRVAl2", &status);
    ephem->crpix1 = drms_getkey_double(inRec, "CRPIX1", &status);
    ephem->crpix2 = drms_getkey_double(inRec, "CRPIX2", &status);
    ephem->cdelt1 = drms_getkey_double(inRec, "CDELT1", &status);
    ephem->cdelt2 = drms_getkey_double(inRec, "CDELT2", &status);
    
    ephem->asd = drms_getkey_double(inRec, "RSUN_OBS", &status);
    if (status) {
        double dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
        double rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
        if (status) rSun_ref = 6.96e8;
        ephem->asd = asin(rSun_ref / dSun);
    }
    ephem->rSun = ephem->asd * RAD2ARCSEC / ephem->cdelt1;
    
    ephem->t_rec = drms_getkey_time(inRec, "T_REC", &status);
    ephem->t_obs = drms_getkey_time(inRec, "T_OBS", &status);
    
    ephem->ctype1 = drms_getkey_string(inRec, "CTYPE1", &status);
    ephem->ctype2 = drms_getkey_string(inRec, "CTYPE2", &status);
    
    return 0;
    
}



// =============================================

/*
 * Perform interpolation using lib function
 *
 */

void performSampling(float *inData, float *outData, float *u, float *v,
                     int *dims_in, int *dims_out)
{
    // Use lib interpolation method
    struct fint_struct pars;
    
    switch (INTERPOPT) {
        case 0:			// Wiener, 6 order, 1 constraint
            init_finterpolate_wiener(&pars, 6, 1, 6, 2, 1, 1, NULL, dpath);
            break;
        case 1:			// Cubic convolution
            init_finterpolate_cubic_conv(&pars, 1., 3.);
            break;
        case 2:			// Bilinear
            init_finterpolate_linear(&pars, 1.);
            break;
        case 3:			// Near neighbor
        default:
            break;
    }
    
    int ind;
    int ncol = dims_out[0], nrow = dims_out[1];
    if (INTERPOPT == 3) {			// Aug 6 2013, Xudong
        for (int row = 0; row < nrow; row++) {
            for (int col = 0; col < ncol; col++) {
                ind = row * ncol + col;
                outData[ind] = nnb(inData, dims_in[0], dims_in[1], u[ind], v[ind]);
            }
        }
    } else {
        // lib function
        finterpolate(&pars, inData, u, v, outData,
                     dims_in[0], dims_in[1], dims_in[0], ncol, nrow, ncol, DRMS_MISSING_FLOAT);
    }
    
}


// =============================================

/*
 * Near neighbor interpolation
 *
 */

float nnb (float *f, int nx, int ny, float x, float y)
{
    
    if (x <= -0.5 || y <= -0.5 || x > nx - 0.5 || y > ny - 0.5)
        return DRMS_MISSING_FLOAT;
    int ilow = floor (x);
    int jlow = floor (y);
    int i = ((x - ilow) > 0.5) ? ilow + 1 : ilow;
    int j = ((y - jlow) > 0.5) ? jlow + 1 : jlow;
    return f[j * nx + i];
    
}
