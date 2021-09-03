/*
 *  hmib2ptr.c
 *
 *	This module takes hmi.B_720s and convert into three field components (Bp, Bt, Br)
 *  and optionally Stony Hurst longitude, latitude
 *  and three error images (Bp_err, Bt_err, Br_err)
 *  Potentially to be used in export system
 *
 *	Author:
 *		Xudong Sun
 *
 *	Version:
 *              v0.0 Aug 28 2015
 *
 *	Notes:
 *		v0.0
 *      Need to check externally that the input is hmi.B_720s
 *
 *	Example Calls:
 *      hmib2ptr "in=hmi.B_720s[2011.02.15_00:00]" "out=hmi_test.Bptr_720s" "requestid=test" -e -l
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
#include "fstats.h"
//#include "cartography.c"
//#include "img2helioVector.c"
//#include "vecErrProp.c"      // new version of errorprop

#define PI              (M_PI)
#define TWOPI           (2*M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2    (16777216)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

// Some other things
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define FREE_ARR(arr) {if (arr) {free(arr);}}
#define DRMS_FREE_ARR(arr) {if (arr && (arr->data)) {drms_free_array(arr);}}

#define kNotSpecified "Not Specified"

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

// Ephemeris information
struct ephemeris {
    double disk_lonc, disk_latc;
    double disk_xc, disk_yc;
    double rSun, asd, pa;
};

//=================================================================================================

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Convert field vector, calls img2helioVector */
void get_bptr(struct ephemeris *ephem, int *dims, float *fld, float *inc, float *azi,
              float *bp, float *bt, float *br, float *lon, float *lat);

/* Compute error, calls errprop_rad */
void get_bptr_err(struct ephemeris *ephem, int *dims, float *fld, float *inc, float *azi,
                  float *lon, float *lat,
                  float *err_fld, float *err_inc, float *err_azi, float *cc_fi, float *cc_fa, float *cc_ia,
                  float *err_bp, float *err_bt, float *err_br);

/* Vector transformation */
int img2helioVector (double bxImg, double byImg, double bzImg, double *bxHelio, 
                     double *byHelio, double *bzHelio, double lon, double lat,
                     double lonc, double latc, double pAng);

/* Error propagation */
void vecErrProp(float fld, float inc, float azi,
                float var_fld, float var_inc, float var_azi, float cov_fi, float cov_fa, float cov_ia,
                double lon, double lat, double lonc, double latc, double pa,
                double *var_bp, double *var_bt, double *var_br);

/* From Cartography.c */
int img2sphere (double x, double y, double ang_r, double latc, double lonc,
    double pa, double *rho, double *lat, double *lon, double *sinlat,
    double *coslat, double *sig, double *mu, double *chi);

//=================================================================================================

char *module_name = "sharp";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input B series."},
    {ARG_STRING, "out", kNotSpecified, "Input Bptr series."},
    {ARG_STRING, "requestid", kNotSpecified, "Request ID."},
    {ARG_INT, "ambweak", "2", "Disambiguation method. 0: potential acute, 1: random, 2: radial acute"},
    {ARG_FLAG, "l", "0", "Flag for lat/lon output."},
    {ARG_FLAG, "e", "0", "Flag for error output."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get parameters */
    
    char *inQuery = (char *) params_get_str(&cmdparams, "in");
    char *outQuery = (char *) params_get_str(&cmdparams, "out");
    char *requestid = (char *) params_get_str(&cmdparams, "requestid");
    int ambweak = params_get_int(&cmdparams, "ambweak");
    int do_lonlat = params_isflagset(&cmdparams, "l");
    int do_error = params_isflagset(&cmdparams, "e");
    
    if (ambweak < 0 || ambweak > 2) ambweak = 2;
    
    /* Input and output data */

    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    
    if (!inRS)
    {
        DIE("No records specified by record-set specification.\n");
    }
    
    int nrecs = inRS->n;
    if (status || nrecs == 0) DIE("Input records error.");
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("Output records not created.");
    }
    
    /* Work */
    
    printf("Processing %d records:\n", nrecs);
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        // Input record
        DRMS_Record_t *inRec = inRS->records[irec];
        
        if (!inRec)
        {
            fprintf(stderr, "No record with index %d in record set.\n", irec);
            continue;
        }
        
        TIME t_rec = drms_getkey_time(inRec, "T_REC", &status);
        
        if (status != DRMS_SUCCESS)
        {
            fprintf(stderr, "Unable to obtain keyword value for T_REC.\n");
            continue;            
        }
        
        char t_rec_str[100];
        sprint_time(t_rec_str, t_rec, "TAI", 0);
        printf("Processing record #%d, T_REC=%s\n", irec, t_rec_str);

        // Ephemeris
        
        struct ephemeris ephem;
        status = getEphemeris(inRec, &ephem); // can never return anything but 'success'.
        
        // Input images
        
        DRMS_Segment_t *inSeg_fld = drms_segment_lookup(inRec, "field");
        if (!inSeg_fld)
        {
            fprintf(stderr, "Missing required segment 'field'.\n");
            continue;
        }
        
        DRMS_Segment_t *inSeg_inc = drms_segment_lookup(inRec, "inclination");
        if (!inSeg_inc)
        {
            fprintf(stderr, "Missing required segment 'inclination'.\n");
            continue;
        }
        
        DRMS_Segment_t *inSeg_azi = drms_segment_lookup(inRec, "azimuth");
        if (!inSeg_azi)
        {
            fprintf(stderr, "Missing required segment 'azimuth'.\n");
            continue;
        }
        
        DRMS_Segment_t *inSeg_amb = drms_segment_lookup(inRec, "disambig");
        if (!inSeg_amb)
        {
            fprintf(stderr, "Missing required segment 'disambig'.\n");
            continue;
        }
        
        int status_fld = 0, status_inc = 0, status_azi = 0, status_amb = 0;
        DRMS_Array_t *inArray_fld = NULL, *inArray_inc = NULL, *inArray_azi = NULL, *inArray_amb = NULL;
        inArray_fld = drms_segment_read(inSeg_fld, DRMS_TYPE_FLOAT, &status_fld);
        inArray_inc = drms_segment_read(inSeg_inc, DRMS_TYPE_FLOAT, &status_inc);
        inArray_azi = drms_segment_read(inSeg_azi, DRMS_TYPE_FLOAT, &status_azi);
        inArray_amb = drms_segment_read(inSeg_amb, DRMS_TYPE_CHAR, &status_amb);
        if (status_fld || status_inc || status_azi || status_amb) {
            printf("Reading input arrays error, record skipped\n");
            DRMS_FREE_ARR(inArray_fld); DRMS_FREE_ARR(inArray_inc);
            DRMS_FREE_ARR(inArray_azi); DRMS_FREE_ARR(inArray_amb);
            continue;
        }
        
        // Input error images if needed
        
        DRMS_Array_t *inArray_err_fld = NULL, *inArray_err_inc = NULL, *inArray_err_azi = NULL;
        DRMS_Array_t *inArray_cc_fi = NULL, *inArray_cc_fa = NULL, *inArray_cc_ia = NULL;
        if (do_error) {
            DRMS_Segment_t *inSeg_err_fld = drms_segment_lookup(inRec, "field_err");
            DRMS_Segment_t *inSeg_err_inc = drms_segment_lookup(inRec, "inclination_err");
            DRMS_Segment_t *inSeg_err_azi = drms_segment_lookup(inRec, "azimuth_err");
            DRMS_Segment_t *inSeg_cc_fi = drms_segment_lookup(inRec, "field_inclination_err");
            DRMS_Segment_t *inSeg_cc_fa = drms_segment_lookup(inRec, "field_az_err");
            DRMS_Segment_t *inSeg_cc_ia = drms_segment_lookup(inRec, "inclin_azimuth_err");
            
            if (!inSeg_err_fld || !inSeg_err_inc || !inSeg_err_azi || !inSeg_cc_fi || !inSeg_cc_fa || !inSeg_cc_ia) 
            {
                fprintf(stderr, "Missing one or more required error segments.\n");
                continue;
            }            
            
            int status_err_fld = 0, status_err_inc = 0, status_err_azi = 0;
            int status_cc_fi = 0, status_cc_fa = 0, status_cc_ia = 0;

            inArray_err_fld = drms_segment_read(inSeg_err_fld, DRMS_TYPE_FLOAT, &status_err_fld);
            inArray_err_inc = drms_segment_read(inSeg_err_inc, DRMS_TYPE_FLOAT, &status_err_inc);
            inArray_err_azi = drms_segment_read(inSeg_err_azi, DRMS_TYPE_FLOAT, &status_err_azi);
            inArray_cc_fi = drms_segment_read(inSeg_cc_fi, DRMS_TYPE_FLOAT, &status_cc_fi);
            inArray_cc_fa = drms_segment_read(inSeg_cc_fa, DRMS_TYPE_FLOAT, &status_cc_fa);
            inArray_cc_ia = drms_segment_read(inSeg_cc_ia, DRMS_TYPE_FLOAT, &status_cc_ia);
            if (status_err_fld || status_err_inc || status_err_azi ||
                status_cc_fi || status_cc_fa || status_cc_ia) {
                printf("Reading input uncertainty arrays error, record skipped\n");
                DRMS_FREE_ARR(inArray_err_fld); DRMS_FREE_ARR(inArray_err_inc); DRMS_FREE_ARR(inArray_err_azi);
                DRMS_FREE_ARR(inArray_cc_fi); DRMS_FREE_ARR(inArray_cc_fa); DRMS_FREE_ARR(inArray_cc_ia);
                DRMS_FREE_ARR(inArray_fld); DRMS_FREE_ARR(inArray_inc);
                DRMS_FREE_ARR(inArray_azi); DRMS_FREE_ARR(inArray_amb);
                continue;
            }
        }
        
        // Some pre-processing
        
        int ncols = inArray_fld->axis[0], nrows = inArray_fld->axis[1];     // FOURK
        int npix = ncols * nrows;
        int dims[2] = {ncols, nrows};
        
        float *fld = (float *) (inArray_fld->data);
        float *inc = (float *) (inArray_inc->data);
        float *azi = (float *) (inArray_azi->data);
        char *amb = (char *) (inArray_amb->data);
        for (int ipix = 0; ipix < npix; ipix++) {
            inc[ipix] *= RADSINDEG;
            // Bug fixed Sep 2 2021 XS
            if ((amb[ipix] >> ambweak) % 2) azi[ipix] += 180.;     // bitwise operator
            azi[ipix] *= RADSINDEG;
        }
        
        float *err_fld = NULL, *err_inc = NULL, *err_azi = NULL;
        float *cc_fi = NULL, *cc_fa = NULL, *cc_ia = NULL;
        if (do_error) {
            err_fld = (float *) (inArray_err_fld->data);
            err_inc = (float *) (inArray_err_inc->data);
            err_azi = (float *) (inArray_err_azi->data);
            cc_fi = (float *) (inArray_cc_fi->data);
            cc_fa = (float *) (inArray_cc_fa->data);
            cc_ia = (float *) (inArray_cc_ia->data);
            for (int ipix = 0; ipix < npix; ipix++) {
                err_inc[ipix] *= RADSINDEG;
                err_azi[ipix] *= RADSINDEG;         // conver to radian
            }
        }

        // Bptr, also calculate lon/lat
        
        float *bp = (float *) (malloc(npix * sizeof(float)));
        float *bt = (float *) (malloc(npix * sizeof(float)));
        float *br = (float *) (malloc(npix * sizeof(float)));
        float *lon = (float *) (malloc(npix * sizeof(float)));      // lon/lat in rad
        float *lat = (float *) (malloc(npix * sizeof(float)));
        get_bptr(&ephem, dims, fld, inc, azi, bp, bt, br, lon, lat);    // working function
        
        // Error if requested
        
        float *err_bp = NULL, *err_bt = NULL, *err_br = NULL;
        if (do_error) {
            err_bp = (float *) (malloc(npix * sizeof(float)));
            err_bt = (float *) (malloc(npix * sizeof(float)));
            err_br = (float *) (malloc(npix * sizeof(float)));
            get_bptr_err(&ephem, dims, fld, inc, azi, lon, lat,
                         err_fld, err_inc, err_azi, cc_fi, cc_fa, cc_ia,
                         err_bp, err_bt, err_br);    // working function, all input in rad
        }
        
        // Convert lon/lat to deg if needed
        
        if (do_lonlat) {
            for (int ipix = 0; ipix < npix; ipix++) {
                lon[ipix] /= RADSINDEG;
                lat[ipix] /= RADSINDEG;
            }
        }
        
        // Some clean up
        
        DRMS_FREE_ARR(inArray_fld); DRMS_FREE_ARR(inArray_inc);
        DRMS_FREE_ARR(inArray_azi); DRMS_FREE_ARR(inArray_amb);
        if (do_error) {
            DRMS_FREE_ARR(inArray_err_fld); DRMS_FREE_ARR(inArray_err_inc); DRMS_FREE_ARR(inArray_err_azi);
            DRMS_FREE_ARR(inArray_cc_fi); DRMS_FREE_ARR(inArray_cc_fa); DRMS_FREE_ARR(inArray_cc_ia);
        }
        if (!do_lonlat) {
            FREE_ARR(lon); FREE_ARR(lat);
        }
        
        // Output
        
        DRMS_Record_t *outRec = outRS->records[irec];
        if (!outRec)
        {
            fprintf(stderr, "No output record with index %d in record set.\n", irec);
            continue;
        }
        
        DRMS_Segment_t *outSeg_bp = drms_segment_lookup(outRec, "Bp");
        if (!outSeg_bp)
        {
            fprintf(stderr, "Missing required output segment 'Bp'.\n");
            continue;
        }
        
        DRMS_Segment_t *outSeg_bt = drms_segment_lookup(outRec, "Bt");
        if (!outSeg_bt)
        {
            fprintf(stderr, "Missing required output segment 'Bt'.\n");
            continue;
        }

        DRMS_Segment_t *outSeg_br = drms_segment_lookup(outRec, "Br");
        if (!outSeg_br)
        {
            fprintf(stderr, "Missing required output segment 'Br'.\n");
            continue;
        }
        
        int status_bp = 0, status_bt = 0, status_br = 0;
        DRMS_Array_t *outArray_bp = NULL, *outArray_bt = NULL, *outArray_br = NULL;
        outArray_bp = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, bp, &status_bp);
        outArray_bt = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, bt, &status_bt);
        outArray_br = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, br, &status_br);
        if (status_bp || status_bt || status_br) {
            printf("Writing output arrays error, record skipped\n");
            DRMS_FREE_ARR(outArray_bp); DRMS_FREE_ARR(outArray_bt); DRMS_FREE_ARR(outArray_br);
            continue;
        }
        outArray_bp->israw = outArray_bt->israw = outArray_br->israw = 0;		// always compressed
        outArray_bp->bzero = outSeg_bp->bzero; outArray_bp->bscale = outSeg_bp->bscale;
        outArray_bt->bzero = outSeg_bt->bzero; outArray_bt->bscale = outSeg_bt->bscale;
        outArray_br->bzero = outSeg_br->bzero; outArray_br->bscale = outSeg_br->bscale;

        status_bp = drms_segment_write(outSeg_bp, outArray_bp, 0);
        status_bt = drms_segment_write(outSeg_bt, outArray_bt, 0);
        status_br = drms_segment_write(outSeg_br, outArray_br, 0);
        DRMS_FREE_ARR(outArray_bp); DRMS_FREE_ARR(outArray_bt); DRMS_FREE_ARR(outArray_br);
        
        if (status_bp || status_bt || status_br) {
            printf("Writing output arrays error, record skipped\n");
            continue;
        }
        
        // Optional output
        
        if (do_error) {
            DRMS_Segment_t *outSeg_err_bp = drms_segment_lookup(outRec, "Bp_err");
            if (!outSeg_err_bp)
            {
                fprintf(stderr, "Missing required output segment 'Bp_err'.\n");
                continue;
            }
            
            DRMS_Segment_t *outSeg_err_bt = drms_segment_lookup(outRec, "Bt_err");
            if (!outSeg_err_bt)
            {
                fprintf(stderr, "Missing required output segment 'Bt_err'.\n");
                continue;
            }

            DRMS_Segment_t *outSeg_err_br = drms_segment_lookup(outRec, "Br_err");
            if (!outSeg_err_br)
            {
                fprintf(stderr, "Missing required output segment 'Br_err'.\n");
                continue;
            }
            
            int status_err_bp = 0, status_err_bt = 0, status_err_br = 0;
            DRMS_Array_t *outArray_err_bp = NULL, *outArray_err_bt = NULL, *outArray_err_br = NULL;
            outArray_err_bp = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, err_bp, &status_err_bp);
            outArray_err_bt = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, err_bt, &status_err_bt);
            outArray_err_br = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, err_br, &status_err_br);
            if (status_err_bp || status_err_bt || status_err_br) {
                printf("Writing error output arrays error, record skipped\n");
                DRMS_FREE_ARR(outArray_err_bp); DRMS_FREE_ARR(outArray_err_bt); DRMS_FREE_ARR(outArray_err_br);
                continue;
            }
            outArray_err_bp->israw = outArray_err_bt->israw = outArray_err_br->israw = 0;		// always compressed, fixed Dec 23 2015
            outArray_err_bp->bzero = outSeg_err_bp->bzero; outArray_err_bp->bscale = outSeg_err_bp->bscale;
            outArray_err_bt->bzero = outSeg_err_bt->bzero; outArray_err_bt->bscale = outSeg_err_bt->bscale;
            outArray_err_br->bzero = outSeg_err_br->bzero; outArray_err_br->bscale = outSeg_err_br->bscale;
            
            status_err_bp = drms_segment_write(outSeg_err_bp, outArray_err_bp, 0);
            status_err_bt = drms_segment_write(outSeg_err_bt, outArray_err_bt, 0);
            status_err_br = drms_segment_write(outSeg_err_br, outArray_err_br, 0);
            DRMS_FREE_ARR(outArray_err_bp); DRMS_FREE_ARR(outArray_err_bt); DRMS_FREE_ARR(outArray_err_br);
            
            if (status_err_bp || status_err_bt || status_err_br) {
                printf("Writing error output arrays error, record skipped\n");
                continue;
            }
        }
        
        if (do_lonlat) {
            DRMS_Segment_t *outSeg_lon = drms_segment_lookup(outRec, "lon");
            if (!outSeg_lon)
            {
                fprintf(stderr, "Missing required output segment 'lon'.\n");
                continue;
            }
            
            DRMS_Segment_t *outSeg_lat = drms_segment_lookup(outRec, "lat");
            if (!outSeg_lon)
            {
                fprintf(stderr, "Missing required output segment 'lat'.\n");
                continue;
            }

            int status_lon = 0, status_lat = 0;
            DRMS_Array_t *outArray_lon = NULL, *outArray_lat = NULL;
            outArray_lon = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, lon, &status_lon);
            outArray_lat = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, lat, &status_lat);
            if (status_lon || status_lat) {
                printf("Writing lon/lat arrays error, record skipped\n");
                DRMS_FREE_ARR(outArray_lon); DRMS_FREE_ARR(outArray_lat);
                continue;
            }
            outArray_lon->israw = outArray_lat->israw = 0;		// always compressed
            outArray_lon->bzero = outSeg_lon->bzero; outArray_lon->bscale = outSeg_lon->bscale;
            outArray_lat->bzero = outSeg_lat->bzero; outArray_lat->bscale = outSeg_lat->bscale;
            
            status_lon = drms_segment_write(outSeg_lon, outArray_lon, 0);
            status_lat = drms_segment_write(outSeg_lat, outArray_lat, 0);
            DRMS_FREE_ARR(outArray_lon); DRMS_FREE_ARR(outArray_lat);
            
            if (status_lon || status_lat) {
                printf("Writing lon/lat arrays error, record skipped\n");
                continue;
            }
        }
        
        // Keywords

        drms_copykeys(outRec, inRec, 0, 0);     // copy all keys
        
        TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(outRec, "DATE", tnow);
        drms_setkey_int(outRec, "AMBWEAK", ambweak);
        drms_setkey_string(outRec, "RequestID", requestid);
        drms_setkey_string(outRec, "BUNIT_000", "Mx/cm^2");
        drms_setkey_string(outRec, "BUNIT_001", "Mx/cm^2");
        drms_setkey_string(outRec, "BUNIT_002", "Mx/cm^2");
        drms_setkey_string(outRec, "BUNIT_003", "Mx/cm^2");
        drms_setkey_string(outRec, "BUNIT_004", "Mx/cm^2");
        drms_setkey_string(outRec, "BUNIT_005", "Mx/cm^2");
        drms_setkey_string(outRec, "BUNIT_006", "degree");
        drms_setkey_string(outRec, "BUNIT_007", "degree");
        
    } // irec
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    return 0;
    
}	// DoIt


//=================================================================================================

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
    
    float crvalx = 0.0;
    float crvaly = 0.0;
    float crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);
    float crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
    float cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
    ephem->disk_xc = PIX_X(0.0,0.0) - 1.0;		// Center of disk in pixel, starting at 0
    ephem->disk_yc = PIX_Y(0.0,0.0) - 1.0;
    
    float dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
    float rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
    if (status) rSun_ref = 6.96e8;
    
    ephem->asd = asin(rSun_ref/dSun);
    ephem->rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;       // in pixel
    
    return 0;
    
}

/* Convert field vector, calls img2helioVector */
void get_bptr(struct ephemeris *ephem, int *dims, float *fld, float *inc, float *azi,
              float *bp, float *bt, float *br, float *lon, float *lat)
{
    int ncols = dims[0], nrows = dims[1];
    int ipix = 0;
    
    double lonc = 0., latc = ephem->disk_latc;        // iput
    double xc = ephem->disk_xc, yc = ephem->disk_yc;
    double asd = ephem->asd, pa = ephem->pa, rSun = ephem->rSun;
    double x, y, b_xi, b_eta, b_zeta;
    double rho, lat_t, lon_t, sinlat_t, coslat_t, sig, mu, chi;     // output
    double bx_t, by_t, bz_t;
//    printf("xc=%f, yc=%f, rSun=%f\n", xc, yc, rSun);
    
    for (int row = 0; row < nrows; row++) {
        for (int col = 0; col < ncols; col++) {
            ipix = row * ncols + col;
            
            // lon/lat
            x = (col - xc) / rSun; y = (row - yc) / rSun;
            if (img2sphere(x, y, asd, latc, lonc, pa,
                           &rho, &lat_t, &lon_t, &sinlat_t, &coslat_t,
                           &sig, &mu, &chi)) {                      // fail
                lon[ipix] = lat[ipix] = DRMS_MISSING_FLOAT;
                bp[ipix] = bt[ipix] = br[ipix] = DRMS_MISSING_FLOAT;
                continue;
            }
            if (lon_t > PI) lon_t -= TWOPI;
            lon[ipix] = lon_t; lat[ipix] = lat_t;
            
            // Bvec
            b_xi = - fld[ipix] * sin(inc[ipix]) * sin(azi[ipix]);
            b_eta = fld[ipix] * sin(inc[ipix]) * cos(azi[ipix]);
            b_zeta = fld[ipix] * cos(inc[ipix]);
            img2helioVector(b_xi, b_eta, b_zeta, &bx_t, &by_t, &bz_t,
                            lon_t, lat_t, lonc, latc, pa);
            bp[ipix] = bx_t;
            bt[ipix] = by_t * (-1);
            br[ipix] = bz_t;
        }
    }
}

/* Compute error, calls errprop_rad */
void get_bptr_err(struct ephemeris *ephem, int *dims, float *fld, float *inc, float *azi,
                  float *lon, float *lat,
                  float *err_fld, float *err_inc, float *err_azi, float *cc_fi, float *cc_fa, float *cc_ia,
                  float *err_bp, float *err_bt, float *err_br)
{
    int ncols = dims[0], nrows = dims[1];
    int npix = ncols * nrows;
    
    double lonc = 0., latc = ephem->disk_latc, pa = ephem->pa;		// Dec 23 2015, changed lonc from disc_lonc to 0.
    
    float var_fld, var_inc, var_azi;            // variance
    float cov_fi, cov_fa, cov_ia;               // covariance
    
    double var_bp, var_bt, var_br;
    
    for (int ipix = 0; ipix < npix; ipix++) {
        if (isnan(lon[ipix])) {
        	err_bp[ipix] = err_bt[ipix] = err_br[ipix] = DRMS_MISSING_FLOAT;
        	continue;
        }
        
        // variance & covariance
        var_fld = err_fld[ipix] * err_fld[ipix];
        var_inc = err_inc[ipix] * err_inc[ipix];
        var_azi = err_azi[ipix] * err_azi[ipix];
        cov_fi = cc_fi[ipix] * err_fld[ipix] * err_inc[ipix];
        cov_fa = cc_fa[ipix] * err_fld[ipix] * err_azi[ipix];
        cov_ia = cc_ia[ipix] * err_inc[ipix] * err_azi[ipix];
        
        // Do it
        vecErrProp(fld[ipix], inc[ipix], azi[ipix], var_fld, var_inc, var_azi, cov_fi, cov_fa, cov_ia,
                   lon[ipix], lat[ipix], lonc, latc, pa, &var_bp, &var_bt, &var_br);
        err_bp[ipix] = sqrt(var_bp);
        err_bt[ipix] = sqrt(var_bt);
        err_br[ipix] = sqrt(var_br);
    }
}

// From img2helioVector.c
int img2helioVector (double bxImg, double byImg, double bzImg, double *bxHelio, 
                     double *byHelio, double *bzHelio, double lon, double lat,
                     double lonc, double latc, double pAng) {
/*
*     perform tramsformation of a vector from image location (lon, lat) to
*     heliographic center. The formula is from Hagyard (1987), and further
*     developed by Gary & Hagyard (1990).
*
* Arguments:
*
*    bxImg, byImg, bzImg: three components of vector magnetic field on image
*                         coordinates.
*    lon, lat:            heliographic coordinates of the location where the vector field
*                         measured. They are in radians.
*    lonc, latc:          heliographic coordinates of the image disk center. They are in
*                         radians.
*    pAng:                position angle of the heliographic north pole, measured eastward
*                         from the north. It's in radians.   
*/

    static double raddeg = M_PI / 180.;
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;

    a11 = -sin(latc) * sin(pAng) * sin(lon - lonc) + cos(pAng) * cos(lon - lonc);
    a12 =  sin(latc) * cos(pAng) * sin(lon - lonc) + sin(pAng) * cos(lon - lonc);
    a13 = -cos(latc) * sin(lon - lonc);
    a21 = -sin(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc)) 
          - cos(lat) * cos(latc) * sin(pAng);
    a22 =  sin(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc))
          + cos(lat) * cos(latc) * cos(pAng);
    a23 = -cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat);
    a31 =  cos(lat) * (sin(latc) * sin(pAng) * cos(lon - lonc) + cos(pAng) * sin(lon - lonc))
          - sin(lat) * cos(latc) * sin(pAng);
    a32 = -cos(lat) * (sin(latc) * cos(pAng) * cos(lon - lonc) - sin(pAng) * sin(lon - lonc))
          + sin(lat) * cos(latc) * cos(pAng);
    a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc);

    *bxHelio = a11 * bxImg + a12 * byImg + a13 * bzImg;
    *byHelio = a21 * bxImg + a22 * byImg + a23 * bzImg;
    *bzHelio = a31 * bxImg + a32 * byImg + a33 * bzImg;

return 0;
}


// Error propagation for vector field
// Adapted from Yang's errorprop.c
// All angles in radian
// By or Bt=-By doesn't matter (all terms are squared)

void vecErrProp(float fld, float inc, float azi,
                float var_fld, float var_inc, float var_azi, float cov_fi, float cov_fa, float cov_ia,
                double lon, double lat, double lonc, double latc, double pa,
                double *var_bp, double *var_bt, double *var_br)

{
    
    // Coefficients
    
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
    a11 = - sin(latc) * sin(pa) * sin(lon - lonc) + cos(pa) * cos(lon - lonc);
    a12 =  sin(latc) * cos(pa) * sin(lon - lonc) + sin(pa) * cos(lon - lonc);
    a13 = - cos(latc) * sin(lon - lonc);
    a21 = - sin(lat) * (sin(latc) * sin(pa) * cos(lon - lonc) + cos(pa) * sin(lon - lonc))
            - cos(lat) * cos(latc) * sin(pa);
    a22 =  sin(lat) * (sin(latc) * cos(pa) * cos(lon - lonc) - sin(pa) * sin(lon - lonc))
            + cos(lat) * cos(latc) * cos(pa);
    a23 = - cos(latc) * sin(lat) * cos(lon - lonc) + sin(latc) * cos(lat);
    a31 =  cos(lat) * (sin(latc) * sin(pa) * cos(lon - lonc) + cos(pa) * sin(lon - lonc))
            - sin(lat) * cos(latc) * sin(pa);
    a32 = - cos(lat) * (sin(latc) * cos(pa) * cos(lon - lonc) - sin(pa) * sin(lon - lonc))
            + sin(lat) * cos(latc) * cos(pa);
    a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc);
    
    // Derivative
    
    double dBpdfld, dBpdinc, dBpdazi;
    double dBtdfld, dBtdinc, dBtdazi;
    double dBrdfld, dBrdinc, dBrdazi;
    
    dBpdfld = (- a11 * sin(inc) * sin(azi) + a12 * sin(inc) * cos(azi) + a13 * cos(inc));
    dBpdinc = fld * (- a11 * cos(inc) * sin(azi) + a12 * cos(inc) * cos(azi) - a13 * sin(inc));
    dBpdazi = fld * (- a11 * sin(inc) * cos(azi) - a12 * sin(inc) * sin(azi));
    
    dBtdfld = (- a21 * sin(inc) * sin(azi) + a22 * sin(inc) * cos(azi) + a23 * cos(inc));
    dBtdinc = fld * (- a21 * cos(inc) * sin(azi) + a22 * cos(inc) * cos(azi) - a23 * sin(inc));
    dBtdazi = fld * (- a21 * sin(inc) * cos(azi) - a22 * sin(inc) * sin(azi));
    
    dBrdfld = (- a31 * sin(inc) * sin(azi) + a32 * sin(inc) * cos(azi) + a33 * cos(inc));
    dBrdinc = fld * (- a31 * cos(inc) * sin(azi) + a32 * cos(inc) * cos(azi) - a33 * sin(inc));
    dBrdazi = fld * (- a31 * sin(inc) * cos(azi) - a32 * sin(inc) * sin(azi));
    
    // Sum
    
    *var_bp = dBpdfld * dBpdfld * var_fld + dBpdinc * dBpdinc * var_inc + dBpdazi * dBpdazi * var_azi +
                2.0 * dBpdfld * dBpdinc * cov_fi + 2.0 * dBpdfld * dBpdazi * cov_fa + 2.0 * dBpdinc * dBpdazi * cov_ia;
    
    *var_bt = dBtdfld * dBtdfld * var_fld + dBtdinc * dBtdinc * var_inc + dBtdazi * dBtdazi * var_azi +
                2.0 * dBtdfld * dBtdinc * cov_fi + 2.0 * dBtdfld * dBtdazi * cov_fa + 2.0 * dBtdinc * dBtdazi * cov_ia;
    
    *var_br = dBrdfld * dBrdfld * var_fld + dBrdinc * dBrdinc * var_inc + dBrdazi * dBrdazi * var_azi +
                2.0 * dBrdfld * dBrdinc * cov_fi + 2.0 * dBrdfld * dBrdazi * cov_fa + 2.0 * dBrdinc * dBrdazi * cov_ia;
    
}

// from Cartography.c
int img2sphere (double x, double y, double ang_r, double latc, double lonc,
    double pa, double *rho, double *lat, double *lon, double *sinlat,
    double *coslat, double *sig, double *mu, double *chi) {
/*
 *  Map projected coordinates (x, y) to (lon, lat) and (rho | sig, chi)
 *  
 *  Arguments:
 *    x }	    Plate locations, in units of the image radius, relative
 *    y }		to the image center
 *    ang_r	    Apparent semi-diameter of sun (angular radius of sun at
 *			the observer's tangent line)
 *    latc	    Latitude of disc center, uncorrected for light travel time
 *    lonc	    Longitude of disc center
 *    pa	    Position angle of solar north on image, measured eastward
 *			from north (sky coordinates)
 *  Return values:
 *    rho	    Angle point:sun_center:observer
 *    lon	    Heliographic longitude
 *    lat	    Heliographic latitude
 *    sinlat	    sine of heliographic latitude
 *    coslat	    cosine of heliographic latitude
 *    sig	    Angle point:observer:sun_center
 *    mu	    cosine of angle between the point:observer line and the
 *			local normal
 *    chi	    Position angle on image measured westward from solar
 *			north
 *
 *  All angles are in radians.
 *  Return value is 1 if point is outside solar radius (in which case the
 *    heliographic coordinates and mu are meaningless), 0 otherwise.
 *  It is assumed that the image is direct; the x or y coordinates require a
 *    sign change if the image is inverted.
 *
 */
  static double ang_r0 = 0.0, sinang_r = 0.0, tanang_r = 0.0;
  static double latc0 = 0.0, coslatc = 1.0, sinlatc = 0.0;
  double cosr, sinr, sinlon, sinsig;

  if (ang_r != ang_r0) {
    sinang_r = sin (ang_r);
    tanang_r = tan (ang_r);
    ang_r0 = ang_r;
  }
  if (latc != latc0) {
    sinlatc = sin (latc);
    coslatc = cos (latc);
    latc0 = latc;
  }
  *chi = atan2 (x, y) + pa;
  while (*chi > 2 * M_PI) *chi -= 2 * M_PI;
  while (*chi < 0.0) *chi += 2 * M_PI;
             /*  Camera curvature correction, no small angle approximations  */
  *sig = atan (hypot (x, y) * tanang_r);
  sinsig = sin (*sig);
  *rho = asin (sinsig / sinang_r) - *sig;
  if (*sig > ang_r) return (-1);
  *mu = cos (*rho + *sig);
  sinr = sin (*rho);
  cosr = cos (*rho);

  *sinlat = sinlatc * cos (*rho) + coslatc * sinr * cos (*chi);
  *coslat = sqrt (1.0 - *sinlat * *sinlat);
  *lat = asin (*sinlat);
  sinlon = sinr * sin (*chi) / *coslat;
  *lon = asin (sinlon);
  if (cosr < (*sinlat * sinlatc)) *lon = M_PI - *lon;
  *lon += lonc;
  while (*lon < 0.0) *lon += 2 * M_PI;
  while (*lon >= 2 * M_PI) *lon -= 2 * M_PI;
  return (0);
}
