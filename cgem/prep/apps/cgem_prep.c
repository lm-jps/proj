/*
 *  cgem_prep.c
 *
 *  Example:
 %  cgem_prep "in=hmi_test.bvcutout_720s[11158][2011.02.15_12:00/252m]" -s "tref=2011.02.15_12:00:00_TAI" "cols=600" "rows=600" "lonref=18.6" "latref=-20.4" "dx=0.03" "dy=0.03"
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "astro.h"
#include "cartography.c"
#include "img2helioVector.c"
#include "finterpolate.h"

// Legacy macros

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)
#define DTTHRESH        (518400.)					// 6 day
#define INTERPOPT       (0)         // Interpolation

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"
#define dpath    "/home/jsoc/cvs/Development/JSOC"

#define NT          (5)             // chunk size
#define FLIP_THR    (120)           // azimuth flipping threshold

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

// WSC code
char *wcsCode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
    "SIN", "ZEA"};


// Ephemeris information
struct ephemeris {
	double disk_lonc, disk_latc;
	double disk_xc, disk_yc;
    double rSun, asd, pa;
    TIME t_rec;
};

// Requested mapping information
struct reqInfo {
    double lonref, latref;
    int ncol, nrow;
    double dx, dy;
    TIME tref;
    int proj;
};


// =====================================

// Doppler correction subroutine, return 0 if success
int correct_dop(DRMS_Record_t *inRec, DRMS_Record_t *outRec);

// Azimuth correction subroutine, return 0 if success
int correct_azi(DRMS_Record_t *inRec, DRMS_Record_t *outRec, float *azi_work,
                int edgeRec, int nx, int ny, int nt);

// Mapping subroutine, return 0 if success
int cgem_mapping(DRMS_Record_t *inRec, DRMS_Record_t *outRec, struct reqInfo *req);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Compute the pixel coordinates at which the map is sampled */
void findCoord(struct reqInfo *req, struct ephemeris *ephem,
			   float *u_out, float *v_out, float *mapCenter);

/* Performing local vector transformation */
void vectorTransform(float *fld, float *inc, float *az, float *bx, float *by, float *bz,
                     int *dims, struct ephemeris *ephem);

/* Convert LOS vectors to local xyz coordinate */
void getLOSVector(struct reqInfo *req, struct ephemeris *ephem, float *mapCenter,
                  float *lx, float *ly, float *lz);

/* Mapping function */
void performSampling(float *inData, float *outData, float *u, float *v, int *dims_in, int *dims_out);

/* Nearest neighbor interpolation */
float nnb (float *f, int nx, int ny, float x, float y);

/* Azimuth correction */
extern void azim_mindiff_jsoc_(float *azi_work, float *azi_out, int *nx, int *ny, int *nt2, int *flip_thr);

// =====================================

char *module_name = "cgem_prep";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input data series."},
    {ARG_STRING, "medv", "hmi_test.medv", "Intermediate series for Doppler correction."},
    {ARG_STRING, "meda", "hmi_test.meda", "Intermediate series for azimuth correction."},
    {ARG_STRING, "out", "hmi_test.cgem_prep", "Output data series."},
    {ARG_FLAG, "a", "", "No azimuth correction."},
    {ARG_FLAG, "s", "", "Save intermediate steps."},
    {ARG_NUME, "map", "carree", "Projetion method, carree by default.",
        "carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
    {ARG_FLOAT, "lonref", "0", "Reference patch center Stonyhurst lon, in deg."},
	{ARG_FLOAT, "latref", "0", "Reference patch center Stonyhurst lat, in deg."},
	{ARG_INT, "cols", "500", "Columns of output cutout."},
	{ARG_INT, "rows", "500", "Rows of output cutout."},
	{ARG_FLOAT, "dx", "0.03", "X pixel size, unit depending on projection (default deg)."},
	{ARG_FLOAT, "dy", "0.03", "Y pixel size, unit depending on projection (default deg)."},
	{ARG_STRING, "tref", kNotSpecified, "Reference time."},

    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;

    /* Get data series */
    
    char *inQuery, *outQuery, *medvQuery, *medaQuery;    // Data series query
    
    inQuery = (char *) params_get_str(&cmdparams, "in");
    outQuery = (char *) params_get_str(&cmdparams, "out");
    medvQuery = (char *) params_get_str(&cmdparams, "medv");
    medaQuery = (char *) params_get_str(&cmdparams, "meda");
    
    /* Arguments */
    
    int do_azi = !(params_isflagset(&cmdparams, "a"));      // skip azi crrection if flag is set
    int saveInterm = params_isflagset(&cmdparams, "s");     // save intermediate steps, no saving by default (no flag)
    int nt = round(fabs(NT));
    
    /* Mapping arguments */
    
    struct reqInfo req;
    req.lonref = params_get_float(&cmdparams, "lonref") * RADSINDEG;
    req.latref = params_get_float(&cmdparams, "latref") * RADSINDEG;
    req.ncol = params_get_int(&cmdparams, "cols");
    req.nrow = params_get_int(&cmdparams, "rows");
    req.tref = params_get_time(&cmdparams, "tref");
    req.dx = params_get_float(&cmdparams, "dx") * RADSINDEG;		// deg to rad, for now
    req.dy = params_get_float(&cmdparams, "dy") * RADSINDEG;
    req.proj = params_get_int(&cmdparams, "map");
    
    // Some checking
    DRMS_RecLifetime_t intLife = saveInterm ? DRMS_PERMANENT : DRMS_TRANSIENT;      // type of intermediate series
    int closeAction = saveInterm ? DRMS_INSERT_RECORD : DRMS_FREE_RECORD;        // type of closing method
    nt = (nt % 2) ? nt : nt + 1;        // odd
    int nt2 = nt / 2;
//    printf("map=%d, nt=%d, nt2=%d\n", req.proj, nt, nt2);
    
    /* Input Data */
    // Assumes records are sorted by T_REC (for single CGEMNUM it should be so)
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs == 0 || !inRS) DIE("Input data series error");
    if (do_azi && nrecs < nt) DIE("Input data series not long enough for azi correction");
    
    // Check reference time, if more than 6 day from first record, use first record
    
    TIME t0 = drms_getkey_time(inRS->records[0], "T_REC", &status);
    if (fabs(req.tref - t0) > DTTHRESH) {
        req.tref = t0;
    }
    
    // All cutouts must have same size
    
    int nx, ny;
    for (int irec = 0; irec < nrecs; irec++) {
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, "azimuth");
        if (!inSeg) DIE("Can't check array dimension!");
        if (irec == 0) {
            nx = inSeg->axis[0]; ny = inSeg->axis[1];
            continue;
        }
        if ((inSeg->axis[0] != nx) || (inSeg->axis[1] != ny)) {
            DIE("Array sizes do not match!");
        }
    }
    int nxny = nx * ny;
    
    // =============================================
    
    /* Correcting Doppler */
    
    printf("==============\nStart correcting Doppler. %d image(s) in total.\n", nrecs);
    
    // Temporary dataset
    
    DRMS_RecordSet_t *medvRS = drms_create_records(drms_env, nrecs, medvQuery, intLife, &status);
    if (status || !medvRS) {DIE("Error in Doppler correciton series");}
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *medvRec = medvRS->records[irec];
        
        printf("Frame %d of %d... ", irec + 1, nrecs);
        
        if (correct_dop(inRec, medvRec)) {          // error
            printf("Doppler correction error, frame %s skipped.\n", irec);
            continue;
        }
        
        printf("done.\n");
        
    }
    
    drms_close_records(inRS, DRMS_FREE_RECORD);     // original data series closed
    
    // =============================================
    
    /* Correcting Azimuth */
    
    DRMS_RecordSet_t *medaRS = NULL;
    
    if (do_azi) {
        
        printf("==============\nStart correcting Azimuth.\n");
        
        medaRS = drms_create_records(drms_env, nrecs, medaQuery, intLife, &status);
        if (status || !medaRS) {DIE("Error in azimuth correciton series");}
        
        // Read in all azimuth
        
        printf("Reading all azimuths... ");
        
        float *azi_cube = (float *) (malloc(nxny * nrecs * sizeof(float)));
        for (int irec = 0; irec < nrecs; irec++) {
            
            DRMS_Record_t *medvRec = medvRS->records[irec];
            DRMS_Segment_t *inSeg = drms_segment_lookup(medvRec, "azimuth");
            if (!inSeg) DIE("Can't retrieve azimtuh!");
            
            DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status || !inArray) return 1;
            
            float *inData = (float *)inArray->data;
            float *azi_cube_lead = azi_cube + nxny * irec;      // first loc
            memcpy(azi_cube_lead, inData, nxny * sizeof(float));        // copy over
            
            drms_free_array(inArray);
            
        }
        
        printf("done.\n");
        
        // Construct in/out array, adapted from Welsch's azim_main.f90
        
        float *azi_work = (float *) (malloc(nxny * nt * sizeof(float)));    // working array
        
        for (int irec = 0; irec < nrecs; irec++) {
            
            DRMS_Record_t *medvRec = medvRS->records[irec];
            DRMS_Record_t *medaRec = medaRS->records[irec];
            
            printf("Frame %d of %d... ", irec + 1, nrecs);
            
            int edgeRec = ((irec < nt2) || (irec >= (nrecs - nt2)));	// first/last NT/2 recs
            float *azi_work_lead, *azi_cube_lead;
            
            // Place desired frame at #nt2
            
            if (irec <= nt2) {      // first nt2 + 1
                
                azi_work_lead = azi_work + (nt2 - irec) * nxny;
                azi_cube_lead = azi_cube;
                memcpy(azi_work_lead, azi_cube_lead, nxny * (nt - nt2 + irec) * sizeof(float));
                
            } else if (irec >= (nrecs - nt2)) {     // last nt2
                
                azi_work_lead = azi_work;
                azi_cube_lead = azi_cube + (irec - nt2) * nxny;
                memcpy(azi_work_lead, azi_cube_lead, nxny * (nrecs - irec + nt2) * sizeof(float));
                
            } else {
                
                // frame #0-#nt is corrected, so effectively propagated backward
                
                for (int i = 0; i < nt - 1; i++) {
                    memcpy(azi_work + i * nxny, azi_work + (i + 1) * nxny, nxny * sizeof(float));
                }
                azi_work_lead = azi_work + (nt - 1) * nxny;
                azi_cube_lead = azi_cube + (irec + nt2) * nxny;
                memcpy(azi_work_lead, azi_cube_lead, nxny * sizeof(float));
                
            }
            
            // Correct azimuth, calling Fortran code
            // For first or last nt2 frames, simply copy over frame nt+1
            // For other frames, results written out directly, but also saved at
            // frame #nt2 in azi_work
            
            if (correct_azi(medvRec, medaRec, azi_work, edgeRec, nx, ny, nt)) {
                printf("Azimuth correction error, frame %s skipped.\n", irec);
                continue;
            }
            
            printf("done.\n");
            
        }
        
        // Clean up
        
        free(azi_work);
        free(azi_cube);
        
        
    } else {
        
        printf("Azimuth correction skipped by request.\n");
        
    }
    
    // =============================================
    
    /* Mapping and compute LOS vectors */
    
    // Do nothing for now
    
    if (do_azi) {
        inRS = medaRS;              // use corrected azimuth
    } else {
        inRS = medvRS;              // use uncorrected azimuth
    }

    printf("==============\nStart mapping.\n");
        
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status || !outRS) {DIE("Error in output series");}
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *outRec = outRS->records[irec];
        
        printf("Frame %d of %d... ", irec + 1, nrecs);
        
        if (cgem_mapping(inRec, outRec, &req)) {          // error
            printf("Mapping error, frame %d skipped.\n", irec);
            continue;
        }
        
        printf("done.\n");
        
    }

    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    // =============================================
    
    /* Clean up */
    
    drms_close_records(medvRS, closeAction);
    if (do_azi) drms_close_records(medaRS, closeAction);
    
    return 0;

}   // DoIt

// =============================================

/* 
 * Top level Doppler correction subroutine
 */

int correct_dop(DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
    int status = 0;
    
    // Read azimuth
    
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, "Dopplergram");
    if (!inSeg) return 1;
    
    DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status || !inArray) return 1;
    
    float *dop_in = (float *) inArray->data;          // actual array
    int nx = inArray->axis[0], ny = inArray->axis[1], nxny = nx * ny;
    int dims[2] = {nx, ny};
    
    // Correction
    
    float *dop_out = (float *) (malloc(nxny * sizeof(float)));      // output array
    double vr = drms_getkey_double(inRec, "OBS_VR", &status);    // Keywords
    
    // ===================
    // DO CORRECTIONS HERE
    // ===================
    
    // copy over & subtract vr as placeholder during testing
    for (int i = 0; i < nxny; i++) {
    	dop_out[i] = dop_in[i] - vr;
    }
    float ddop = 0;			// correction
    drms_free_array(inArray);
    
    // Output
    
    DRMS_Segment_t *outSeg = drms_segment_lookup(outRec, "Dopplergram");
    if (!outSeg) return 1;
    
    DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, dop_out, &status);
    if (status) {free(dop_out); return 1;}
    
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) return status;
    
    drms_free_array(outArray);
    
    // Links and keywords
    
    drms_copykeys(outRec, inRec, 0, 0);     // copy all keys
    drms_setkey_float(outRec, "DDOP", ddop);		// correction
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    DRMS_Link_t *link = hcon_lookup_lower(&outRec->links, "DATA");
    if (link) drms_link_set("DATA", outRec, inRec);
    
    return 0;
}

// =============================================

/*
 * Top level azimuth correction subroutine
 */

int correct_azi(DRMS_Record_t *inRec, DRMS_Record_t *outRec, float *azi_work,
                int edgeRec, int nx, int ny, int nt)
{
    int status = 0;
    
    int nt2 = nt / 2;
    int flip_thr = FLIP_THR;
    int nxny = nx * ny;
    
    float *azi_out = (float *) (malloc(nxny * sizeof(float)));      // output array
    
    // For first nt2 or last nt2, simply copy over
    // for the rest, do correction, then copy back to frame #nt2 in azi_work
    
    if (edgeRec) {
        
        memcpy(azi_out, azi_work + nxny * nt2, nxny * sizeof(float));
        
    } else {
        
        azim_mindiff_jsoc_(azi_work, azi_out, &nx, &ny, &nt2, &flip_thr);
        
    }
    
    // copy back
    
    memcpy(azi_work + nxny * nt2, azi_out, nxny * sizeof(float));
    
    // Output
    
    DRMS_Segment_t *outSeg = drms_segment_lookup(outRec, "azimuth");
    if (!outSeg) return 1;
    
    int dims[2] = {nx, ny};
    DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, azi_out, &status);
    if (status) {free(azi_out); return 1;}
    
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) return status;
    
    drms_free_array(outArray);
    
    // Links and keywords
    
    drms_copykeys(outRec, inRec, 0, 0);     // copy all keys
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    DRMS_Link_t *link = hcon_lookup_lower(&outRec->links, "DATA");
    if (link) drms_link_set("DATA", outRec, inRec);
    
    return 0;
    
}

// =============================================

/*
 * Top level mapping subroutine
 * Also obtain LOS vectors
 */

int cgem_mapping(DRMS_Record_t *inRec, DRMS_Record_t *outRec, struct reqInfo *req)
{
    int status = 0;
    
    // Map info
    
    struct ephemeris ephem;
    
    if (getEphemeris(inRec, &ephem)) {
        SHOW("Mapping: get ephemeris error\n");
        return 1;
    }

    // Read data
    
    DRMS_Segment_t *inSeg;
    DRMS_Array_t *inArray_fld, *inArray_inc, *inArray_azi, *inArray_dop;
    
    float *fld, *inc, *azi, *dop;
    
    inSeg = drms_segment_lookup(inRec, "field");
    inArray_fld = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    fld = (float *)inArray_fld->data;
    
    inSeg = drms_segment_lookup(inRec, "inclination");
    inArray_inc = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    inc = (float *)inArray_inc->data;
    
    inSeg = drms_segment_lookup(inRec, "azimuth");
    inArray_azi = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    azi = (float *)inArray_azi->data;
    
    inSeg = drms_segment_lookup(inRec, "Dopplergram");
    inArray_dop = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    dop = (float *)inArray_dop->data;
    
    int dims_in[2] = {inArray_fld->axis[0], inArray_fld->axis[1]};
    
    // Vector transform
    
    float *bx = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    float *by = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    float *bz = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    
    vectorTransform(fld, inc, azi, bx, by, bz, dims_in, &ephem);
    
    // Find position
    
    int dims_out[2] = {req->ncol, req->nrow};
    float *u_out = (float *) (malloc(dims_out[0] * dims_out[1] * sizeof(float)));
    float *v_out = (float *) (malloc(dims_out[0] * dims_out[1] * sizeof(float)));
    
    float mapCenter[2];		// map center coordinate
    findCoord(req, &ephem, u_out, v_out, mapCenter);
    
    // Sampling
    
    float *bx_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));
    float *by_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));
    float *bz_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));
    float *dop_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));
    
    performSampling(bx, bx_map, u_out, v_out, dims_in, dims_out);
    performSampling(by, by_map, u_out, v_out, dims_in, dims_out);
    performSampling(bz, bz_map, u_out, v_out, dims_in, dims_out);
    performSampling(dop, dop_map, u_out, v_out, dims_in, dims_out);
    
    // Obtain LOS vectors lx, ly, lz
    
    float *lx_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));
    float *ly_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));
    float *lz_map = (float *) (malloc(req->ncol * req->nrow * sizeof(float)));

    getLOSVector(req, &ephem, mapCenter, lx_map, ly_map, lz_map);
    
    // Clean up
    
    drms_free_array(inArray_fld);
    drms_free_array(inArray_inc);
    drms_free_array(inArray_azi);
    drms_free_array(inArray_dop);
    free(u_out);
    free(v_out);
    
    // Output
    
    DRMS_Segment_t *outSeg;
    DRMS_Array_t *outArray;
    
    outSeg = drms_segment_lookup(outRec, "Bx");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, bx_map, &status);
    if (status) {free(bx_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }
    
    outSeg = drms_segment_lookup(outRec, "By");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, by_map, &status);
    if (status) {free(by_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }
    
    outSeg = drms_segment_lookup(outRec, "Bz");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, bz_map, &status);
    if (status) {free(bz_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }
    
    outSeg = drms_segment_lookup(outRec, "Dopplergram");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, dop_map, &status);
    if (status) {free(dop_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }
    
    outSeg = drms_segment_lookup(outRec, "Lx");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, lx_map, &status);
    if (status) {free(lx_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }
    
    outSeg = drms_segment_lookup(outRec, "Ly");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, ly_map, &status);
    if (status) {free(ly_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }
    
    outSeg = drms_segment_lookup(outRec, "Lz");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, lz_map, &status);
    if (status) {free(lz_map); return 1;}
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;		// always compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) { return status; } else { drms_free_array(outArray); }

	// keys

    drms_copykeys(outRec, inRec, 0, 0);     // copy all keys
    
    drms_setkey_float(outRec, "CRPIX1", req->ncol/2. + 0.5);
	drms_setkey_float(outRec, "CRPIX2", req->nrow/2. + 0.5);
	drms_setkey_float(outRec, "CRVAL1", mapCenter[0]);
	drms_setkey_float(outRec, "CRVAL2", mapCenter[1]);
	drms_setkey_float(outRec, "CDELT1", req->dx / RADSINDEG);
	drms_setkey_float(outRec, "CDELT2", req->dy / RADSINDEG);
	drms_setkey_string(outRec, "CUNIT1", "degree");
	drms_setkey_string(outRec, "CUNIT2", "degree");
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
		
	char key[64];
	snprintf (key, 64, "CRLN-%s", wcsCode[req->proj]);
	drms_setkey_string(outRec, "CTYPE1", key);
	snprintf (key, 64, "CRLT-%s", wcsCode[req->proj]);
	drms_setkey_string(outRec, "CTYPE2", key);
	drms_setkey_float(outRec, "CROTA2", 0.0);
    
    return 0;
}


// =============================================

/*
 * Convert LOS vectors to local xyz coordinate
 *
 */

void getLOSVector(struct reqInfo *req, struct ephemeris *ephem, float *mapCenter,
                  float *lx, float *ly, float *lz)
{
    
    double lonc = mapCenter[0] * RADSINDEG;
    double latc = mapCenter[1] * RADSINDEG;
    
    double x, y;        // coordinate in final map
    double lat, lon;    // lon, lat of each pix
    double lx0, ly0, lz0;
    
    int ind_map;
    for (int row = 0; row < req->nrow; row++) {
        for (int col = 0; col < req->ncol; col++) {
            
            ind_map = row * req->ncol + col;
            x = (col + 0.5 - req->ncol/2.) * req->dx;
            y = (row + 0.5 - req->nrow/2.) * req->dy;
            
            if (plane2sphere (x, y, latc, lonc, &lat, &lon, req->proj)) {
                lx[ind_map] = DRMS_MISSING_FLOAT;
                ly[ind_map] = DRMS_MISSING_FLOAT;
                lz[ind_map] = DRMS_MISSING_FLOAT;
                continue;
            }
            
            // convert (disk center lon = 0)
            
            img2helioVector (0.0, 0.0, 1.0, &lx0, &ly0, &lz0,
                             lon, lat, 0.0, ephem->disk_latc, ephem->pa);
            
            lx[ind_map] = lx0;
            ly[ind_map] = ly0;
            lz[ind_map] = lz0;
            
        }
    }
    
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
    
    float crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation
    double sina = sin(crota2 * RADSINDEG);
    double cosa = cos(crota2 * RADSINDEG);
    
    ephem->pa = - crota2 * RADSINDEG;
    ephem->disk_latc = drms_getkey_float(inRec, "CRLT_OBS", &status) * RADSINDEG;
    ephem->disk_lonc = drms_getkey_float(inRec, "CRLN_OBS", &status) * RADSINDEG;
    
    float crvalx = drms_getkey_float(inRec, "CRVAL1", &status);
    float crvaly = drms_getkey_float(inRec, "CRVAl2", &status);
    float crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);
    float crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
    float cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
    ephem->disk_xc = PIX_X(0.0,0.0) - 1.0;		// Center of disk in pixel, starting at 0
    ephem->disk_yc = PIX_Y(0.0,0.0) - 1.0;
    
    float dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
    float rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
    if (status) rSun_ref = 6.96e8;
    
    ephem->asd = asin(rSun_ref/dSun);
    ephem->rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
    
    ephem->t_rec = drms_getkey_time(inRec, "T_REC", &status);
    
    return 0;
    
}


// =============================================


/*
 * Compute the coordinates to be sampled on cutout image
 *
 */

void findCoord(struct reqInfo *req, struct ephemeris *ephem, 
			   float *u_out, float *v_out, float *mapCenter)
{
 
    // Rotate
    // Differential rotation rate in urad/s
    // proj/lev1.5_hmi/libs/lev15/rotcoef_file.txt
    
    double latc = req->latref;
    double difr = 2.7139 - 0.405 * pow(sin(latc),2) - 0.422 * pow(sin(latc),4);
    double dt = ephem->t_rec - req->tref;
    double lonc = req->lonref + dt * difr * 1.0e-6;      // urad/s to rad
    mapCenter[0] = lonc / RADSINDEG;
    mapCenter[1] = latc / RADSINDEG;
    
    // print some info
    char tstr[100];
    sprint_time(tstr, ephem->t_rec, "TAI", 0);
//    printf("t_rec=%s, lonc=%f\n", tstr, lonc/RADSINDEG);
    
    double x, y;        // coordinate in final map
    double lat, lon;    // lon, lat of each pix
    double u, v;    // pixel address
    
    int ind_map;
    for (int row = 0; row < req->nrow; row++) {
        for (int col = 0; col < req->ncol; col++) {
            
            ind_map = row * req->ncol + col;
            x = (col + 0.5 - req->ncol/2.) * req->dx;
            y = (row + 0.5 - req->nrow/2.) * req->dy;
            
            //
            
            if (plane2sphere (x, y, latc, lonc, &lat, &lon, req->proj)) {
                u_out[ind_map] = -1;
                v_out[ind_map] = -1;
                continue;
            }
            
            if (sphere2img (lat, lon, ephem->disk_latc, 0.0, &u, &v,
                            ephem->disk_xc, ephem->disk_yc, ephem->rSun, ephem->pa,
                            0., 0., 0., 0.)) {
                u_out[ind_map] = -1;
                v_out[ind_map] = -1;
                continue;
            }
            
            u_out[ind_map] = u;
            v_out[ind_map] = v;

        }
    }

}


// =============================================

/*
 * Performing local vector transformation
 *  xyz: z refers to vertical (radial) component, x EW (phi), y NS
 *
 */

void vectorTransform(float *fld, float *inc, float *azi, float *bx, float *by, float *bz,
                     int *dims, struct ephemeris *ephem)
{
    
    int ncol = dims[0], nrow = dims[1];
    
    double b_xi, b_eta, b_zeta;
    double bx0, by0, bz0;
    
    double x, y;        // solar x, y
    double lat, lon;	// lat / lon for current point
    double rho, sinlat, coslat, sig, mu, chi;
    
    int ind;
    for (int row = 0; row < nrow; row++) {
        for (int col = 0; col < ncol; col++) {
            
            ind = row * ncol + col;
            
            x = (col - ephem->disk_xc) / ephem->rSun;      // normalized by disk radius
            y = (row - ephem->disk_yc) / ephem->rSun;
            
            // find lat/lon (disk center lon = 0)
            if (img2sphere(x, y, ephem->asd, ephem->disk_latc, 0.0, ephem->pa,
                           &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi)) {
                bx[ind] = by[ind] = bz[ind] = DRMS_MISSING_FLOAT;
                continue;
            }
            
            // field vector
            
            b_xi = - fld[ind] * sin(inc[ind] * RADSINDEG) * sin(azi[ind] * RADSINDEG);
            b_eta = fld[ind] * sin(inc[ind] * RADSINDEG) * cos(azi[ind] * RADSINDEG);
            b_zeta = fld[ind] * cos(inc[ind] * RADSINDEG);
            
            // convert (disk center lon = 0)
            
            img2helioVector (b_xi, b_eta, b_zeta, &bx0, &by0, &bz0,
                             lon, lat, 0.0, ephem->disk_latc, ephem->pa);
            
            bx[ind] = bx0;
            by[ind] = by0;
            bz[ind] = bz0;
            
        }
    }

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
