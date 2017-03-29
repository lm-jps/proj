/*
 * Module name:		update_poldb.c
 *
 * Description:
 *          Collect 360 deg of polar field data from Synoptic maps
 *          centered at Mar 06/Sep 08, save orig and smoothed version
 *          maybe overridden by T_CEN if not available
 *          Fit down to 60 latitude, using Stereographic and bilinear regardless
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *                  Image fitting code by Rasmus Larsen (rmlarsen@gmail.com)
 *
 * External function:
 *			heliographic_coord.c    polar_remap.c   interp.c    fitimage.c
 *
 * Version:
 *			v1.0		Sep 10 2013
 *
 * Issues:
 *			v1.0
 *
 *
 * Example:
 *  update_poldb year=2012 ns=0 in=hmi.synoptic_mr_720s out=su_xudong.hmi_poldb scale=1 -o
 *
 */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

// Heliographic coord, lib function
#include "heliographic_coords.c"

// Polar projection
#include "polar_remap.c"

// Image fitting
#include "fitimage.c"

/* ################## Preset ################## */

// Image size
#define IMGDIM	(600)
// Lower limit of latitude for fitting
#define LAT0    (60.)

#define LAT1    (62.)
#define LAT2    (75.)

#define EPSL (1.e-5)
#define SECSINDAY (86400.0)
#define	DTOR	(M_PI / 180.)

#define S_DATE ".03.06_02:00:00_UTC"
#define N_DATE ".09.08_02:00:00_UTC"

#define EPOCH_DIFF ((double)220924792.000) /* 1970.01.01_00:00:00_UTC */

#define INSEGNAME	"synopMr"

/* ############# Macros ############# */

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}
#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s (status=%d)\n", msg, status); return(status);}

#define FREE_ARR(arr) {if (arr) {free(arr);}}

/* ############# Keywords ############# */

char *synop_keys[] = {
    "CTYPE1", "CTYPE2", "CDELT1", "CDELT2", "CUNIT1", "CUNIT2", "WCSNAME"};

/* ################## Aux function ################## */

/* Get current year */
int currentYear(void)
{
    time_t timeval;
    struct tm *tp;
    time (&timeval);
    tp = gmtime(&timeval);
    return (tp->tm_year + 1900);
}

/* Convert TIME to current year */
int time2Year(TIME t)
{
    time_t timeval = (time_t)(t + EPOCH_DIFF);
    struct tm *tp = gmtime(&timeval);
    return (tp->tm_year + 1900);
}

/* ################## Checking functions ################## */

/* Check parameters, set center time */
int checkParams(int year, int ns, int nmax, float scale, TIME *t_cen_ptr);

/* Check output availability */
int checkOutput(char *outRSName, int year, int ns, int ovwt);

/* Check input availability, determine crot, etc. */
int checkInput(char *inRSName, char *insegName, int ns, TIME t_cen,
                 int *crot, int *crot1, double *lonc);

/* Check input record, keys, called by checkInput */
int checkInputRec(char *inRSName, char *insegName, int carr,
                  double *cdelt1, double *cdelt2, char **ctype1, char **ctype2, int *nx, int *ny);

/* ################## Working functions ################## */

/* Top level working wrapper
 * Get polar data, remap, smooth, remap back
 * Image allocated inside the function
 */
int createPole(float **data_raw, float **data, float **img_raw, float **img,
               int *dataDims, int *imgDims,
               char *inRSName, char *insegName, int ns, int nmax, int crot, int crot1, double *lonc);

/* Get maps from two CRs, obtain lat/lon and
 * dimensions of the polar data
 */
int patchMaps(double **map, int *mapDims, double **lons, double **lats, int *dataDims,
              char *inRSName, char *insegName, int crot, int crot1, double *lonc);

/* ################## Output functions ################## */

/* Output data, free arrays inside and reset to NULL */
int writeOutput(float *data_raw, float *data, float *img_raw, float *img, int *dataDims, int *imgDims,
                char *inRSName, char *outRSName, char *insegName, int year, int ns, int nmax, int crot, double lonc, TIME t_cen, float scale);

/* Compute mean field between LAT1 and LAT2 and above LAT2 */
void getMean(float *mean1, float *mean2, float *data, int *dims,
             int ns, double crpix2, double crval2, double cdelt2);

/* ################################################# */
/* ################## Main Module ################## */
/* ################################################# */

// Manually correct keywords
//#define CORRECT_EPHEM 1

char *module_name = "update_poldb";

ModuleArgs_t module_args[] =
{
    {ARG_INT, "year", NULL, "Year of the data"},
    {ARG_INT, "ns", NULL, "North or South, 0 for N, 1 for S"},
    {ARG_INT, "nmax", "5", "Max order of Chebychev for smoothing"},
    {ARG_STRING, "in", NULL, "Input data series, name only."},
    {ARG_STRING, "out", NULL, "Output data series, name only."},
    {ARG_FLOAT, "scale", "1.", "Multiplier to scale the polar field"},
    {ARG_TIME, "t_cen", "1980.01.01_UTC", "Manual-set center time"},
    {ARG_STRING, "inseg", INSEGNAME, "segment name"},
    {ARG_FLAG, "o", "", "Flag to overwrite existing data records"},
    {ARG_END}
};

int DoIt(void)
{

    int status = DRMS_SUCCESS;
    
    /* ============== */
    /* Get parameters */
    /* ============== */
    
    char *inRSName = (char *)params_get_str(&cmdparams, "in");
    char *outRSName = (char *)params_get_str(&cmdparams, "out");
    char *insegName = (char *)params_get_str(&cmdparams, "inseg");
    int year = params_get_int(&cmdparams, "year");
    int ns = params_get_int(&cmdparams, "ns");
    int nmax = params_get_int(&cmdparams, "nmax");
    float scale = params_get_float(&cmdparams, "scale");
    TIME t_cen = params_get_time(&cmdparams, "t_cen");
    int ovwt = params_isflagset(&cmdparams, "o");       // overwrite existing
    
    /* ================ */
    /* Check parameters */
    /* ================ */
    
    SHOW("Check parameters...\n");
    
    // t_cen will be updated if use default
    status = checkParams(year, ns, nmax, scale, &t_cen);
    
    if (status) DIE("Stop.");
    
    /* ================== */
    /* Open output record */
    /* ================== */
    
    SHOW("Check output data series...\n");
    
    status = checkOutput(outRSName, year, ns, ovwt);
    
    if (status) DIE("Stop.");
    
    /* =================== */
    /* Check input records */
    /* =================== */
    
    int crot, crot1;     // map center carr
    double lonc;      // center pixel longitude
    int dataDims[2], imgDims[2] = {IMGDIM, IMGDIM};
    
    SHOW("Check input data series...\n");
    
    // Determine maps to use, and array dimensions
    status = checkInput(inRSName, insegName, ns, t_cen,
                        &crot, &crot1, &lonc);
    
    if (status) DIE("Stop.");

    /* ===================== */
    /* Get data, then smooth */
    /* ===================== */
    
    // data_raw: cutout of polar data
    // data: smoothed polar data, used for filling
    // img: intermediate, projected polar view
    float *data_raw = NULL, *data = NULL, *img_raw = NULL, *img = NULL;
    
    SHOW("Processing...\n");
    
    // All the work!
    status = createPole(&data_raw, &data, &img_raw, &img,
                        dataDims, imgDims,
                        inRSName, insegName, ns, nmax, crot, crot1, &lonc);
    
    if (status) {
        FREE_ARR(data_raw); FREE_ARR(data); FREE_ARR(img); FREE_ARR(img_raw);
        DIE("Stop\n");
    }
    
    /* ============= */
    /* Output record */
    /* ============= */
    
    SHOW("Writing output...\n");
    
    // Arrays freed inside if success
    status = writeOutput(data_raw, data, img_raw, img, dataDims, imgDims,
                         inRSName, outRSName, insegName, year, ns, nmax, crot, lonc, t_cen, scale);
    
    if (status) {
        FREE_ARR(data_raw); FREE_ARR(data);
        FREE_ARR(img_raw); FREE_ARR(img);
        DIE("Stop\n");
    }
    
    /* ============= */
    
    return DRMS_SUCCESS;
    
}   // DoIt


// ========================================

/* Check parameters, set center time */

int checkParams(int year, int ns, int nmax, float scale, TIME *t_cen_ptr)
{
    
    if ((year > currentYear()) || (year < 1996)) {
        SHOW("Year is out of range. ");  // HMI/MDI only
        return 1;
    }
    if ((ns < 0) || (ns > 1)) {
        SHOW("Please set ns to 0 for north, 1 for south. ");
        return 1;
    }
    if (nmax <= 0) {
        SHOW("Please set nmax to a positive integer. ");
        return 1;
    }
    if (scale <= 0) {
        SHOW("Please set scale to be positive. ");
        return 1;
    }
    printf("  Update %s polar field database for year %4d.\n",
           ((ns) ? "southern" : "northern"), year);
    
    // Check time, reject time that differs from the default by 30 days
    
    char t_cen_str[100];
    TIME t_cen_tmp;

    snprintf(t_cen_str, 100, "%4d%s", year, ((ns) ? S_DATE : N_DATE));
    t_cen_tmp = sscan_time(t_cen_str);      // Default, 09.08 for N, 03.06 for S
    
    int useTCen = (year == time2Year(*t_cen_ptr)) &&
                  ((t_cen_tmp - *t_cen_ptr) < (30. * SECSINDAY));   // Compare with default
    
    // Final t_cen
    
    if (!useTCen) {                 // Default center time
        *t_cen_ptr = t_cen_tmp;     // Update t_cen in main()
        printf("  Use default center time: %s\n", t_cen_str); fflush(stdout);
    } else {                        // User specified center time
        sprint_time(t_cen_str, *t_cen_ptr, "UTC", 0);
        printf("  Use user-specified center time: %s\n", t_cen_str); fflush(stdout);
    }
    
    return 0;
    
}


// ========================================

/* Check output availability */

int checkOutput(char *outRSName, int year, int ns, int ovwt)
{

    int status = 0;
    
    char outQuery[100];
    snprintf(outQuery, 100, "%s[%4d][%1d]", outRSName, year, ns);
    
    DRMS_RecordSet_t *outRS = NULL;
    outRS = drms_open_records(drms_env, outQuery, &status);
    if (status) {
        if (outRS) drms_close_records(outRS, DRMS_FREE_RECORD);
        SHOW("Output record set error. ");
        return status;
    }
     
    if (outRS->n != 0 && !ovwt) {
        if (outRS) drms_close_records(outRS, DRMS_FREE_RECORD);
        SHOW("Output record exists, overwriting not requested. ");
        return 1;
    }
    
    return 0;
    
}


// ========================================

/* Check input availability, determine crot, etc */

int checkInput(char *inRSName, char *insegName, int ns, TIME t_cen,
                  int *crot, int *crot1, double *lonc)
{
 
    int status = 0;
    int carr, carr1;
    
    // Center lon, carr_num for t_cen
    
    double lon_cen, lat_cen;
    HeliographicLocation(t_cen, &carr, &lon_cen, &lat_cen);
    printf("  Center pixel: crot=%4d, lon_cen=%f\n", carr, lon_cen); fflush(stdout);
    
    // Input query for crot, check keywords

    int nx = 0, ny = 0;     // output array size
    double cdelt1 = DRMS_MISSING_DOUBLE, cdelt2 = DRMS_MISSING_DOUBLE;
    char *ctype1 = NULL, *ctype2 = NULL;
    
    status = checkInputRec(inRSName, insegName, carr, &cdelt1, &cdelt2, &ctype1, &ctype2, &nx, &ny);
    if (status) return status;
    
    // Generally requires two frames, since lon_cen is not at map center
    // One full map only when the center longitude is within half of cdelt1 of 180-cdelt1/2.
    
    int x_off = round((lon_cen - (180. - fabs(cdelt1) / 2.)) / fabs(cdelt1));
    carr1 = (x_off == 0) ? carr : ((x_off < 0) ? (carr + 1) : (carr - 1));
    
    // Input query for a second map, crot+/-1
    
    int nx_1 = 0, ny_1 = 0;     // output array size
    double cdelt1_1 = DRMS_MISSING_DOUBLE, cdelt2_1 = DRMS_MISSING_DOUBLE;
    char *ctype1_1 = NULL, *ctype2_1 = NULL;
    
    status = checkInputRec(inRSName, insegName, carr1, &cdelt1_1, &cdelt2_1, &ctype1_1, &ctype2_1, &nx_1, &ny_1);
    if (status) return status;
    
    // Map dimensions and types must match
    
    if (carr != carr1) {
        if (isnan(cdelt1) || isnan (cdelt2) || isnan(cdelt1_1) || isnan(cdelt2_1) ||
            !ctype1 || !ctype2 || !ctype1_1 || !ctype2_1 ||
            fabs(cdelt1 - cdelt1_1) > EPSL ||
            fabs(cdelt2 - cdelt2_1) > EPSL ||
            strcmp(ctype1, ctype1_1) != 0 ||
            strcmp(ctype2, ctype2_1) != 0 ||
            (strcmp(ctype2, "CRLT-CEA") != 0 && strcmp(ctype2, "Sine Latitude") != 0) ||  // has to be sine-lat
            !nx || !ny || !nx_1 || !ny_1 ||
            nx != nx_1 || ny != ny_1)
        {
            SHOW("  Keywords of two input records don't match.\n");
            FREE_ARR(ctype1); FREE_ARR(ctype2);
            FREE_ARR(ctype1_1); FREE_ARR(ctype2_1);
            return 1;
        }
    }
    
    // Output
    
    *crot = carr;
    *crot1 = carr1;
    *lonc = lon_cen;
    
    // Clean up
    
    FREE_ARR(ctype1); FREE_ARR(ctype2);
    FREE_ARR(ctype1_1); FREE_ARR(ctype2_1);
    
    return 0;
}


// ========================================

/* Check input record, keys, called by checkInput */

int checkInputRec(char *inRSName, char *insegName, int carr,
                  double *cdelt1, double *cdelt2, char **ctype1, char **ctype2, int *nx, int *ny)
{
    
    int status = 0;
    
    char inQuery[100];
    snprintf(inQuery, 100, "%s[%4d]", inRSName, carr);
    printf("  Check %s\n", inQuery);
    
    DRMS_RecordSet_t *inRS = NULL;
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) {
        if (inRS) drms_close_records(inRS, DRMS_FREE_RECORD);
        SHOW("  Input record error.\n");
        return ((status) ? status : 1);
    }
    
    DRMS_Record_t *inRec = inRS->records[0];
    *cdelt1 = drms_getkey_double(inRec, "CDELT1", &status);
    *cdelt2 = drms_getkey_double(inRec, "CDELT2", &status);
    *ctype1 = drms_getkey_string(inRec, "CTYPE1", &status);
    *ctype2 = drms_getkey_string(inRec, "CTYPE2", &status);
    
    // Check segment & array
    
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, insegName);
    DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) {
        drms_close_records(inRS, DRMS_FREE_RECORD);
        SHOW(" Array read error.\n");
        return status;
    }
    
    // Determine map size
    
    *nx = inArray->axis[0]; *ny = inArray->axis[1];
    
    // Clean up
    
    drms_free_array(inArray);
    drms_close_records(inRS, DRMS_FREE_RECORD);
    
    printf("    CDELT1=%f, CDELT2=%f, CTYPE1=%s, CTYPE2=%s\n",
           *cdelt1, *cdelt2, *ctype1, *ctype2);
    
    return 0;
    
}


// ========================================

/* Top level working wrapper
 * Get polar data, remap, smooth, remap back
 * Image allocated inside the function
 */

int createPole(float **data_raw, float **data, float **img_raw, float **img,
               int *dataDims, int *imgDims,
               char *inRSName, char *insegName, int ns, int nmax, int crot, int crot1, double *lonc)
{
 
    int status = 0;
    
    // Open maps
    
    double *map = NULL;      // full patched map
    double *lons = NULL, *lats = NULL;     // longitude, latitude
    int mapDims[2];         // full map size
    
    status = patchMaps(&map, mapDims, &lons, &lats, dataDims,       // All these are updated
                       inRSName, insegName, crot, crot1, lonc);
    if (status) return status;
    
    // Allocate arrays
    
    int nx = dataDims[0], ny = dataDims[1], nxny = nx * ny;
    int imgDim2 = imgDims[0] * imgDims[1];
    
    printf("  nx=%d, ny=%d\n", nx, ny);
    
    *data_raw = (float *) (malloc(nxny * sizeof(float)));
    *data = (float *) (malloc(nxny * sizeof(float)));
    *img_raw = (float *) (malloc(imgDim2 * sizeof(float)));
    *img = (float *) (malloc(imgDim2 * sizeof(float)));
    
    double *image = (double *) (malloc(imgDim2 * sizeof(double)));  // operate in double
    
    // data_raw: direct cutout of the polar region
    
    int offset = (ns) ? 0 : (mapDims[0] * (mapDims[1] - ny));      // south start from 0
    
    for (int count = 0; count < nxny; count++) {
        (*data_raw)[count] = map[count+offset];
    }

    // Stereographic projection to polar view
    
    synop2pole(ns, 1, image, imgDims, map, mapDims, lons, lats, (LAT0*DTOR));   // in polar_remap.c
    
    // Smoothing
    
    for (int count = 0; count < imgDim2; count++) (*img_raw)[count] = image[count];     // copy out raw
    fitimage(image, imgDims[0], imgDims[1], nmax, 1);
    for (int count = 0; count < imgDim2; count++) (*img)[count] = image[count];     // copy out
    
    // Stereographic projection back to Carrington map, saved back to map
    
    pole2synop(ns, 1, image, imgDims, map, mapDims, lons, lats, (LAT0*DTOR));   // in polar_remap.c
    
    // data: smoothed polar region
    
    for (int count = 0; count < nxny; count++) {
        (*data)[count] = map[count+offset];
    }
    
    // Clean up
    
    FREE_ARR(image);
    FREE_ARR(map);
    FREE_ARR(lons); FREE_ARR(lats);
    
    return 0;
    
}

// ========================================

/* Get maps from two CRs, obtain lat/lon and 
 * dimensions of the polar data
 * lonc will be updated according to resolution
 */

int patchMaps(double **map, int *mapDims, double **lons, double **lats, int *dataDims,
              char *inRSName, char *insegName, int crot, int crot1, double *lonc)
{
    
    int status = 0;
    char inQuery[100], inQuery1[100];
    DRMS_RecordSet_t *inRS = NULL, *inRS1 = NULL;
    DRMS_Record_t *inRec = NULL, *inRec1 = NULL;
    DRMS_Segment_t *inSeg = NULL, *inSeg1 = NULL;
    DRMS_Array_t *inArray = NULL, *inArray1 = NULL;
    
    // Open CROT

    snprintf(inQuery, 100, "%s[%4d]", inRSName, crot);
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status) {
        SHOW(" Getting data error.\n");
        return status;
    }
    inRec = inRS->records[0];
    inSeg = drms_segment_lookup(inRec, insegName);
    inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) {
        SHOW(" Getting data error.\n");
        if (inArray) drms_free_array(inArray);
        drms_close_records(inRS, DRMS_FREE_RECORD);
        return status;
    }
    
    // Open CROT1 uf necessary
    
    if (crot != crot1) {
        snprintf(inQuery1, 100, "%s[%4d]", inRSName, crot1);
        inRS1 = drms_open_records(drms_env, inQuery1, &status);
        if (status) {
            SHOW(" Getting data error.\n");
            drms_free_array(inArray);
            drms_close_records(inRS, DRMS_FREE_RECORD);
            return status;
        }
        inRec1 = inRS1->records[0];
        inSeg1 = drms_segment_lookup(inRec1, insegName);
        inArray1 = drms_segment_read(inSeg1, DRMS_TYPE_FLOAT, &status);
        if (status) {
            SHOW(" Getting data error.\n");
            drms_free_array(inArray);
            if (inArray1) drms_free_array(inArray1);
            drms_close_records(inRS, DRMS_FREE_RECORD);
            drms_close_records(inRS1, DRMS_FREE_RECORD);
            return status;
        }
    }
    
    // Lon/lat of full map
    
    int xsz = inArray->axis[0], ysz = inArray->axis[1];
    mapDims[0] = xsz; mapDims[1] = ysz;
    
    *lons = (double *) (malloc(xsz * sizeof(double)));
    *lats = (double *) (malloc(ysz * sizeof(double)));
    
    double cdelt1 = drms_getkey_double(inRec, "CDELT1", &status);
    double cdelt2 = drms_getkey_double(inRec, "CDELT2", &status);
    double crpix1 = drms_getkey_double(inRec, "CRPIX1", &status);
    double crpix2 = drms_getkey_double(inRec, "CRPIX2", &status);
    double crval1 = drms_getkey_double(inRec, "CRVAL1", &status);   // Carrington time
#ifdef CORRECT_EPHEM
    crpix1 += 1;        // crpix2 += 0.5; crval1 += fabs(cdelt1) / 2.0;     // no need now, longitude start from 0.1
#endif
    double crval2 = drms_getkey_double(inRec, "CRVAL2", &status);
    printf("cdelt1=%f, crpix1=%f, crval1=%f\n", cdelt1, crpix1, crval1);
    
    // Outmap relative coords
    
    for (int col = 0; col < xsz; col++) {
        (*lons)[col] = (col + 1) * fabs(cdelt1) * DTOR - M_PI;      // -179.9, 180; changed to multiple of 0.1 on Apr 30 XS
    }
    for (int row = 0; row < ysz; row++) {
        (*lats)[row] = asin((row - (crpix2 - 1.)) * cdelt2 + crval2);
    }
    
    // Size of polar data
    
    int row = 0;
    while (row < (ysz / 2)) {
        if (fabs((*lats)[row++] / DTOR) < LAT0) break;
    }
    
    dataDims[0] = xsz; dataDims[1] = row;
    
    printf("  Orig map size: %d x %d\n", xsz, ysz);
    printf("  Polar data size: %d x %d\n", dataDims[0], dataDims[1]);
    
    // Map
    
    int xysz = xsz * ysz;
    *map = (double *) (malloc(xysz * sizeof(double)));      // 3600x1440
    
    float *synop = (float *) inArray->data;
    
    if (crot1 == crot) {

        // Need only one map
        for (int count = 0; count < xysz; count++) {(*map)[count] = synop[count];}      // copy everything

    } else {
        
        float *synop1 = (float *) inArray1->data;
        
        // longitude of map center
        
        double crlon = crot * (double)360. - crval1;            // generally 180
        float colc = (*lonc - crlon) / fabs(cdelt1) + crpix1 - 1;   // column for requested lonc, start from 0
        colc = (xsz % 2) ? round(colc) : (floor(colc) + 0.5);
        *lonc = (colc + 1 - crpix1) * fabs(cdelt1) + crlon;     // actual center longitude, modified
        printf("  colc=%f, lonc=%f\n", colc, *lonc);
        
        // Copy & patch data
        
        int col0;
        if (crot1 > crot) {     // crot on the right
            col0 = round(colc + (xsz - 1) / 2.);        // right edge, col0 included
            for (int row = 0; row < ysz; row++) {
                int offset_s = row * xsz;
                int offset = row * xsz + xsz - 1 - col0;
                for (int col = 0; col <= col0; col++) {(*map)[offset + col] = synop[offset_s + col];}        // 0, col0
                int offset1 = row * xsz - col0 - 1;
                for (int col = col0 + 1; col < xsz; col++) {(*map)[offset1 + col] = synop1[offset_s + col];}    // col0 + 1, xsz-1
            }
        } else {                // crot on the left
            col0 = round(colc - (xsz - 1) / 2.);        // left edge, col0 included
            for (int row = 0; row < ysz; row++) {
                int offset_s = row * xsz;
                int offset = row * xsz - col0;
                for (int col = col0; col < xsz; col++) {(*map)[offset + col] = synop[offset_s + col];}      // col0, xsz-1
                int offset1 = row * xsz + xsz - col0;
                for (int col = 0; col < col0; col++) {(*map)[offset1 + col] = synop1[offset_s + col];}      // 0, col0-1
            }
        }
        printf("  col0=%d, lonc=%f, map[0,100]=%f, map[0,101]=%f\n", col0, *lonc, (*map)[100l*xsz], (*map)[101l*xsz]);
    }
    
    // Clean up
    
    drms_free_array(inArray);
    drms_close_records(inRS, DRMS_FREE_RECORD);
    if (crot1 != crot) {
        drms_free_array(inArray1);
        drms_close_records(inRS1, DRMS_FREE_RECORD);
    }
    
    return 0;
    
}

// ========================================

// ========================================

/* Output data, free arrays inside and reset to NULL */

int writeOutput(float *data_raw, float *data, float *img_raw, float *img, int *dataDims, int *imgDims,
                char *inRSName, char *outRSName, char *insegName, int year, int ns, int nmax, int crot, double lonc, TIME t_cen, float scale)
{
    
    int status = 0;
    
    // Create output
    
    DRMS_Record_t *outRec = drms_create_record(drms_env, outRSName, DRMS_PERMANENT, &status);
    if (status) {
        SHOW("  Creating output record failed.\n");
        return status;
    }
    
    // Segment
    
    DRMS_Segment_t *outSeg_data = NULL, *outSeg_data_raw = NULL;
    DRMS_Segment_t *outSeg_img = NULL, *outSeg_img_raw = NULL;
    outSeg_data = drms_segment_lookup(outRec, "data");          // smoothed
    outSeg_data_raw = drms_segment_lookup(outRec, "data_raw");      // original
    outSeg_img = drms_segment_lookup(outRec, "image");          // polar view
    outSeg_img_raw = drms_segment_lookup(outRec, "image_raw");          // polar view original
    if (!outSeg_data || !outSeg_data_raw || !outSeg_img || !outSeg_img_raw) {
        drms_close_record(outRec, DRMS_FREE_RECORD);
        SHOW("  Creating output segment failed.\n");
        return 1;
    }
    
    // Arrays
    
    DRMS_Array_t *outArray_data = drms_array_create(DRMS_TYPE_FLOAT, 2, dataDims, data, &status);
	if (status) {SHOW("  Creating data array failed.\n"); return 1;}
    DRMS_Array_t *outArray_data_raw = drms_array_create(DRMS_TYPE_FLOAT, 2, dataDims, data_raw, &status);
	if (status) {SHOW("  Creating data_raw array failed.\n"); return 1;}
    DRMS_Array_t *outArray_img = drms_array_create(DRMS_TYPE_FLOAT, 2, imgDims, img, &status);
	if (status) {SHOW("  Creating img array failed.\n"); return 1;}
    DRMS_Array_t *outArray_img_raw = drms_array_create(DRMS_TYPE_FLOAT, 2, imgDims, img_raw, &status);
	if (status) {SHOW("  Creating img array failed.\n"); return 1;}
    
    for (int i = 0; i < 2; i++) {
        outSeg_data->axis[i] = outArray_data->axis[i];
        outSeg_data_raw->axis[i] = outArray_data_raw->axis[i];
        outSeg_img->axis[i] = outArray_img->axis[i];
        outSeg_img_raw->axis[i] = outArray_img_raw->axis[i];
    }
    
    outArray_data->israw = 0;		// always compressed
    outArray_data_raw->israw = 0;
    outArray_img->israw = 0;
    outArray_img_raw->israw = 0;
    
	outArray_data->bzero = outSeg_data->bzero;
    outArray_data_raw->bzero = outSeg_data_raw->bzero;
    outArray_img->bzero = outSeg_img->bzero;
    outArray_img_raw->bzero = outSeg_img_raw->bzero;
    
	outArray_data->bscale = outSeg_data->bscale;
    outArray_data_raw->bscale = outSeg_data_raw->bscale;
    outArray_img->bscale = outSeg_img->bscale;
    outArray_img_raw->bscale = outSeg_img_raw->bscale;
    
    // Write
    
    status = drms_segment_write(outSeg_data, outArray_data, 0);
	if (status) {SHOW("  Writing data array failed.\n"); return 1;}
    status = drms_segment_write(outSeg_data_raw, outArray_data_raw, 0);
	if (status) {SHOW("  Writing data_raw array failed.\n"); return 1;}
    status = drms_segment_write(outSeg_img, outArray_img, 0);
	if (status) {SHOW("  Writing img array failed.\n"); return 1;}
    status = drms_segment_write(outSeg_img_raw, outArray_img_raw, 0);
	if (status) {SHOW("  Writing img_raw array failed.\n"); return 1;}
    
    // Keywords
    
    char inQuery[100];
    snprintf(inQuery, 100, "%s[%4d]", inRSName, crot);
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    DRMS_Record_t *inRec = inRS->records[0];
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, insegName);
    int nKeys = ARRLENGTH(synop_keys);
    for (int iKey = 0; iKey < nKeys; iKey++) {
        drms_copykey(outRec, inRec, synop_keys[iKey]);
    }
    int ysz = inSeg->axis[1];
    double yc = drms_getkey_double(inRec, "CRPIX2", &status);
#ifdef CORRECT_EPHEM
    // yc += 0.5;       // Apr 30
#endif
    double yref = drms_getkey_double(inRec, "CRVAL2", &status);
    
    drms_setkey_int(outRec, "YEAR", year);
    drms_setkey_int(outRec, "NS", ns);
    drms_setkey_int(outRec, "NMAX", nmax);
    drms_setkey_int(outRec, "CAR_ROT", crot);
    drms_setkey_double(outRec, "CARRTIME", crot * (double)360. - lonc);
    drms_setkey_time(outRec, "T_OBS", t_cen);
    drms_setkey_float(outRec, "SCALE", scale);
    
    drms_setkey_double(outRec, "CRPIX1", (1. + dataDims[0]) / 2.);
    double crpix2 = (ns) ? yc : (yc - ysz + dataDims[1]);
    drms_setkey_double(outRec, "CRPIX2", crpix2);
    drms_setkey_double(outRec, "CRVAL1", crot * (double)360. - lonc);
    drms_setkey_double(outRec, "CRVAL2", yref);
    
    char timebuf[1024];
	double UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
	sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
	drms_setkey_string(outRec, "DATE", timebuf);
    
    // Mean field between LAT1 and LAT2 and above LAT2
    
    double cdelt2 = drms_getkey_double(inRec, "CDELT2", &status);
    float mean1, mean2;
    getMean(&mean1, &mean2, data_raw, dataDims, ns, crpix2, yref, cdelt2);
    drms_setkey_float(outRec, "MEAN62", mean1);
    drms_setkey_float(outRec, "MEAN75", mean2);
    
    //
    
    drms_free_array(outArray_data); data = NULL;
    drms_free_array(outArray_data_raw); data_raw = NULL;
    drms_free_array(outArray_img); img = NULL;
    drms_free_array(outArray_img_raw); img_raw = NULL;
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_record(outRec, DRMS_INSERT_RECORD);

    return 0;
    
}

// =====================

/* Compute mean field between LAT1 and LAT2 and above LAT2 */

void getMean(float *mean1, float *mean2, float *data, int *dims,
             int ns, double crpix2, double crval2, double cdelt2)
{
    
    double lat;
    int row1, row2;

    for (int row = 0; row < dims[1]; row++) {
        lat = fabs(asin((row - crpix2) * cdelt2 + crval2) / DTOR);
        if (ns) {       // south
            if (lat >= LAT2) row2 = row;
            if (lat >= LAT1) row1 = row;
        } else {
            if (lat <= LAT2) row2 = row;
            if (lat <= LAT1) row1 = row;
        }
    }
    printf("  row1=%d, row2=%d\n", row1, row2);
    
    //
    
    double sum1 = 0.0, sum2 = 0.0;
    int count1 = 0, count2 = 0;
    double val;
    if (ns) {
        for (int row = 0; row < row2; row++) {
            for (int col = 0; col < dims[0]; col++) {
                val = data[row * dims[0] + col];
                if (isnan(val)) continue;
                sum2 += val; count2++;
            }
        }
        for (int row = row2; row < row1; row++) {
            for (int col = 0; col < dims[0]; col++) {
                val = data[row * dims[0] + col];
                if (isnan(val)) continue;
                sum1 += val; count1++;
            }
        }
    } else {
        for (int row = row1; row < row2; row++) {
            for (int col = 0; col < dims[0]; col++) {
                val = data[row * dims[0] + col];
                if (isnan(val)) continue;
                sum1 += val; count1++;
            }
        }
        for (int row = row2; row < dims[1] - 1; row++) {
            for (int col = 0; col < dims[0]; col++) {
                val = data[row * dims[0] + col];
                if (isnan(val)) continue;
                sum2 += val; count2++;
            }
        }
    }
    
    //
    
    *mean1 = (float) (sum1 / count1);
    *mean2 = (float) (sum2 / count2);
    
    printf("  mean1=%f, mean2=%f\n", *mean1, *mean2);
    
}
