/*
 * Module name:		polcorr.c
 *
 * Description:
 *          
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * External function:
 *			polar_remap.c   interp.c    fitimage.c
 *
 * Version:
 *			v1.0		Sep 10 2013
 *
 * Issues:
 *			v1.0
 *
 *
 * Example:
 *  polcorr in=hmi.synoptic_mr_720s'[2114]' out=su_xudong.synop_polfil method=TEMP_SPAT
 *
 */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "set_statistics.c"

// Heliographic coord, lib function
#include "heliographic_coords.c"

// Polar projection
#include "polar_remap.c"

// Image fitting
#include "fitimage.c"

/* ################## Preset ################## */

// Image size
#define IMGDIM	(600)

#define K_TRANS (1.21)       // Scaling factor to be multiplied on low-lat field

#define SECSINDAY ((double)86400.0)
#define SECSINYEAR    ((double)31536000.)
#define	DTOR	(M_PI / 180.)
#define TWOPI   (2 * M_PI)
#define CARRLONINYEAR (4817.5)           // 365./27.2753*360.

#define EPOCH_DIFF ((double)220924792.000) /* 1970.01.01_00:00:00_UTC */

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

/* ############# Data structures ############# */

/* Option */
struct option {
    int method;
    double lat0, latfil, latfil1;
    int nflag, sflag;
    int order;
    int interpOnly;
    int amb;
};

/* Map attributes */
struct attrib {
    // prime key (year) for polar data base series
    int n_year[2], s_year[2];
    int n_extrap, s_extrap;     // 0 for interp, 1 for forward extrap, -1 for backward extrap
    int sinlat;     // whether the map is equally spaced in sinlat or lat
    int carr;
    double carr_time;          // Carrington time of map center
    double *lons, *lats;
    int mapDims[2];
};

// Polar data base series
// static const char poldb_name[] = "su_xudong.hmi_poldb";
char poldb_name[100];

/* ############# Function templates ############# */

/* Check necessary input keywords, get geometry, check database if needed */
int getInputInfo(DRMS_Record_t *inRec, struct attrib *att, struct option *opt);

/* Check availability of polar data base series */
int checkPoldb(struct attrib *att, struct option *opt);

/* Find year in polar data base */
int findYearInPoldb(struct attrib *att, int y, TIME t, int ns);

/* Default fill method, following Sun et al. (2011) */
int temp_spat(float *outMap, struct attrib *att, struct option *opt);

/* Spatial smoothing only */
int spat_2d(float *outMap, struct attrib *att, struct option *opt);

/* Interpolate polar field using time series of polar field */
int interpPoldb(int ns, int rowfil, double *map, struct attrib *att, struct option *opt);

/* Write output */
int writeOutput(char *outQuery, DRMS_Record_t *inRec, float **outMap,
                struct attrib *att, struct option *opt);


/* ################################################# */
/* ################## Main Module ################## */
/* ################################################# */

//#define DEVELOP 1

char *module_name = "polcorr";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL, "Output data series, name only."},
    {ARG_STRING, "poldb", "hmi.polar_db", "Polar field database, name only."},
    {ARG_NUME, "method", "TEMP_SPAT", "Correction scheme.", "TEMP_SPAT, SPAT_2D"},
    {ARG_DOUBLE, "lat0", "60.", "Start latitude used for fitting."},
    {ARG_DOUBLE, "latfil", "75.", "Start latitude for filling in, no lower than 75."},
    {ARG_DOUBLE, "latfil1", "62.", "Start latitude for reduced noise."},
    {ARG_INT, "nflag", "1", "Filling for north pole: 1 for yes, 0 for no."},
    {ARG_INT, "sflag", "1", "Filling for south pole: 1 for yes, 0 for no."},
    {ARG_INT, "order", "5", "Highest power of polynomials."},
    {ARG_INT, "amb", "2", "disambiguation method, 0 for potential, 1 for random, 2 for radial acute"},
    {ARG_FLAG, "i",	"", "Force to run in interpolation only mode."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    struct option opt;      // argument options
    struct attrib att;      // attribute of each record
    
    /* Get parameters */
    
    char *inQuery = (char *)params_get_str(&cmdparams, "in");
    char *outQuery = (char *)params_get_str(&cmdparams, "out");
    strcpy(poldb_name, (char *)params_get_str(&cmdparams, "poldb"));
    
    opt.method = params_get_int(&cmdparams, "method");
    opt.lat0 = fabs(params_get_double(&cmdparams, "lat0"));
    opt.latfil = fabs(params_get_double(&cmdparams, "latfil"));
    opt.latfil1 = fabs(params_get_double(&cmdparams, "latfil1"));
    opt.nflag = params_get_int(&cmdparams, "nflag");
    opt.sflag = params_get_int(&cmdparams, "sflag");
    opt.order = params_get_int(&cmdparams, "order");
    opt.interpOnly = params_isflagset(&cmdparams, "i");
    opt.amb = params_get_int(&cmdparams, "amb");
    
    if (opt.latfil < 75.) opt.latfil = 75.;
    if (opt.latfil1 < opt.lat0) opt.latfil1 = opt.lat0;
    if (opt.latfil1 > opt.latfil) opt.latfil1 = opt.latfil;

    // printf("method=%d\n", opt.method);
    
    /* Check input records */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) {
        DIE("No input data found.\n");
    }
    int nrecs = inRS->n;
    
    /* Do this for each record */
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        printf("Record %d of %d\n", irec + 1, nrecs);
        
        // Input record
        
        DRMS_Record_t *inRec = inRS->records[irec];
        
        // Check necessary keywords, availability of database
        // Returned lats, lons, carr, carr_time, year of database, all stored in "att"
        
        att.lons = NULL; att.lats = NULL;
        status = getInputInfo(inRec, &att, &opt);       // here
        if (status) {
            SHOW("record skipped.\n");
            FREE_ARR(att.lons); FREE_ARR(att.lats);
            continue;
        }
        
        printf("  N_Y0=%d, N_Y1=%d\n  S_Y0=%d, S_Y1=%d\n",
               att.n_year[0], att.n_year[1], att.s_year[0], att.s_year[1]);
        printf("  sinlat=%d, carr=%d, carr_time=%f\n", att.sinlat, att.carr, att.carr_time);
        
        // Read array
        
        DRMS_Segment_t *inSeg = drms_segment_lookupnum(inRec, 0);
        DRMS_Array_t *inArray = NULL;
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status) {
            SHOW("  Input data not found, record skipped.\n");
            if (!inArray) drms_free_array(inArray);
            continue;
        }
        
        float *map = (float *) inArray->data;
        int nx = att.mapDims[0], ny = att.mapDims[1], nxny = nx * ny;
        
        // Output map
        
        float *outMap = (float *) (malloc(nxny * sizeof(float)));
        for (int count = 0; count < nxny; count++) {
            outMap[count] = map[count];                 // copy over
        }
        drms_free_array(inArray);
        
        // Do the correction

        switch (opt.method) {
            default:
            case (0):
                status = temp_spat(outMap, &att, &opt);
                break;
            case (1):
                status = spat_2d(outMap, &att, &opt);
                break;
        }
       
        FREE_ARR(att.lons); FREE_ARR(att.lats);
        if (status) {
            SHOW("  Polar filling error, record skipped.\n");
        }
        
        // Output

//        FREE_ARR(outMap);

        status = writeOutput(outQuery, inRec, &outMap, &att, &opt);
        if (status) {
            SHOW("record skipped.\n");
            FREE_ARR(outMap);
            continue;
        }
        
    }
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    
    //
    
    return DRMS_SUCCESS;
    
}   // DoIt


// =============================================

/* Check necessary input keywords, get geometry,
 * i.e. lat/lon, sinlat, check database if needed
 */

int getInputInfo(DRMS_Record_t *inRec, struct attrib *att, struct option *opt)
{
    int status = 0;
    
    // Keywords
    
    char *ctype2 = drms_getkey_string(inRec, "CTYPE2", &status);
    if (strcmp(ctype2, "CRLT-CEA") == 0) {
        att->sinlat = 1;
    } else if (strcmp(ctype2, "CRLT-CAR") == 0) {
        att->sinlat = 0;
    } else {
        SHOW("  CTYPE2 unknown, ");
        return 1;
    }
    
    att->carr = drms_getkey_int(inRec, "CAR_ROT", &status);
    att->carr_time = drms_getkey_double(inRec, "CRVAL1", &status);
    if (status) {
        SHOW("  CRVAL1 unknown, ");
        return status;
    }
    
    // Lat, lon
    
    DRMS_Segment_t *inSeg = drms_segment_lookupnum(inRec, 0);       // changed to 
    int xsz = inSeg->axis[0], ysz = inSeg->axis[1];
    att->mapDims[0] = xsz;
    att->mapDims[1] = ysz;
    att->lons = (double *) (malloc(xsz * sizeof(double)));
    att->lats = (double *) (malloc(ysz * sizeof(double)));
    
    double cdelt1 = drms_getkey_double(inRec, "CDELT1", &status);
    double cdelt2 = drms_getkey_double(inRec, "CDELT2", &status);
    double xc = drms_getkey_double(inRec, "CRPIX1", &status);
    double yc = drms_getkey_double(inRec, "CRPIX2", &status);
    double xref = att->carr_time;   // Carrington time
    double yref = drms_getkey_double(inRec, "CRVAL2", &status);
    
    for (int col = 0; col < xsz; col++) {
        att->lons[col] = (col + 1) * fabs(cdelt1) * DTOR - M_PI;      // -179.9, 180
    }
    double y;
    for (int row = 0; row < ysz; row++) {
        y = (row - (yc - 1.)) * (cdelt2) + yref;
        att->lats[row] = (att->sinlat) ? asin(y) : (y * DTOR);
    }
    
    // Determine opt parameter for specific fill methods
    
    switch (opt->method) {
        default:
        case (0):       // TEMP_SPAT
            // Get the records that bracket carr_time (opt->n_year, opt->s_year)
            status = checkPoldb(att, opt);
            if (status) {
                SHOW("  Polar data not available, ");
                return status;
            }
            break;
        case (1):       // SPAT_2D
            break;
    }
    
    // Check interpolation-only mode for TEMP_SPAT
    
    if (opt->interpOnly) {
        // Skip when either is extrapolation
        int n = (opt->nflag && att->n_extrap);
        int s = (opt->sflag && att->s_extrap);
        if (n || s) {
            SHOW("  In interpolation only mode but data base not ready, ");
            return 1;
        }
    }
    
    //
    
    return 0;
}

// =============================================

/* Check availability of polar data base series
 * Store record prime keys (year) in opt->n_year, opt->s_year
 */

int checkPoldb(struct attrib *att, struct option *opt)
{
    
    int status = 0;
    att->n_year[0] = att->n_year[1] = 0;
    att->s_year[0] = att->s_year[1] = 0;
    att->n_extrap = att->s_extrap = 0;
    
    int carr = att->carr;
    double lonc = carr * 360. - att->carr_time;
    TIME t = HeliographicTime(carr, lonc);
    int y = time2Year(t);
    
    printf("y=%d\n", y);
    
    // check north
    
    if (opt->nflag) {
        if (findYearInPoldb(att, y, t, 0)) return 1;
    }
    
    // check south
    
    if (opt->sflag) {
        if (findYearInPoldb(att, y, t, 1)) return 1;
    }
    
    //

    return 0;
    
}

// =============================================

/* Find correponding years in the data base
 * if extrapolation, the element is 0. if exceeds the last
 * record by one year, return error
 */

int findYearInPoldb(struct attrib *att, int y, TIME t, int ns)
{
    
    int status = 0;
    if (ns != 0 && ns != 1) return 1;
    
    // place pointer
    
    int *year = (ns) ? att->s_year : att->n_year;
    int *extrap = (ns) ? &(att->s_extrap) : &(att->n_extrap);       // need to modify
    
    //
    
    char inQuery[100], inQuery_tmp[100];
    TIME t0, t1, t2, tt;
    int y0, y1, y2, yy;
    
    snprintf(inQuery, 100, "%s[%4d-%4d][%d]", poldb_name, y - 1, y + 1, ns);     // 3 years
    DRMS_RecordSet_t *inRS = NULL;
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) {        // 1-3 records
        if (inRS) drms_close_records(inRS, DRMS_FREE_RECORD);
        return 1;
    }

    // Treat extrapolation as interpolation
    // For now, backward extrapolation uses the first two years
    // and a linear extrapolation
    // Forward extrapolation uses the last year, and scales
    // the mean to match the lower latitudes
    
    t0 = drms_getkey_time(inRS->records[0], "T_OBS", &status);
    y0 = drms_getkey_int(inRS->records[0], "YEAR", &status);
    printf("inQuery=%s, n=%d\n", inQuery, inRS->n);
    switch (inRS->n) {
        case 1:
            if (t <= t0) {
                // need to have more than 2 years
                snprintf(inQuery_tmp, 100, "%s[%4d][%d]", poldb_name, y0 + 1, ns);
                DRMS_RecordSet_t *inRS_tmp = NULL;
                inRS_tmp = drms_open_records(drms_env, inQuery_tmp, &status);
                if (status || inRS->n != 1) {
                    if (inRS_tmp) drms_close_records(inRS_tmp, DRMS_FREE_RECORD);
                    return 1;
                }
                drms_close_records(inRS_tmp, DRMS_FREE_RECORD);
                //
                year[0] = y0; year[1] = y0 + 1; *extrap = -1;
            }
            if (t > t0) {year[0] = y0; year[1] = y0; *extrap = 1;}
            break;
        case 2:
            t1 = drms_getkey_time(inRS->records[1], "T_OBS", &status);
            y1 = drms_getkey_int(inRS->records[1], "YEAR", &status);
            if (t1 < t0) {tt = t1; t1 = t0; t0 = tt; yy = y1; y1 = y0; y0 = yy;}
            if (fabs(t1 - t0) >= (2.1 * SECSINYEAR)) { status = 1; break; }
            if (t <= t0) {      // backward extrap
                year[0] = y0; year[1] = y1; *extrap = -1;
            } else if (t >= t1) {       // forward extrap
                year[0] = y1; year[1] = y1; *extrap = 1;
            } else {        // interp
                year[0] = y0; year[1] = y1; *extrap = 0;
            }
            break;
        case 3:
            t1 = drms_getkey_time(inRS->records[1], "T_OBS", &status);
            y1 = drms_getkey_int(inRS->records[1], "YEAR", &status);
            t2 = drms_getkey_time(inRS->records[2], "T_OBS", &status);
            y2 = drms_getkey_int(inRS->records[2], "YEAR", &status);
            if (t1 < t0) {tt = t1; t1 = t0; t0 = tt; yy = y1; y1 = y0; y0 = yy;}
            if (t2 < t1) {tt = t2; t2 = t1; t1 = tt; yy = y2; y2 = y1; y1 = yy;}
            if (t1 < t0) {tt = t1; t1 = t0; t0 = tt; yy = y1; y1 = y0; y0 = yy;}
            *extrap = 0;
            if (t < t1) {
                year[0] = y0; year[1] = y1;
            } else {
                year[0] = y1; year[1] = y2;
            }
            break;
        default:
            break;
    }
    drms_close_records(inRS, DRMS_FREE_RECORD);

    //
    
    return status;
}

// =============================================

/* Default fill method, following Sun et al. (2011) */

int temp_spat(float *outMap, struct attrib *att, struct option *opt)
{
    
    // Copy to double
    
    int xsz = att->mapDims[0], ysz = att->mapDims[1], xysz = xsz * ysz;
    double *map = (double *) (calloc(xysz, sizeof(double)));
    for (int count = 0; count < xysz; count++) {
        map[count] = outMap[count];
    }
    
    // Projected image
    
    int imgDims[2] = {IMGDIM, IMGDIM}, imgDim2 = imgDims[0] * imgDims[1];
    double *img_n = (double *) (calloc(imgDim2, sizeof(double)));
    double *img_s = (double *) (calloc(imgDim2, sizeof(double)));
    
    // Create mask to transition from smoothing filling (above latfil)
    // to data (above latfil)
    
    int rowfil = 0, rowfil1 = 0;
    while (rowfil < (ysz / 2)) {
        if (fabs(att->lats[rowfil++] / DTOR) < opt->latfil) break;
    }
    while (rowfil1 < (ysz / 2)) {
        if (fabs(att->lats[rowfil1++] / DTOR) < opt->latfil1) break;
    }
    
    printf("rowfil=%d, rowfil1=%d\n", rowfil, rowfil1);
    
    int nxny1 = xsz * rowfil1, nxny = xsz * rowfil;
    float *mask_n = (float *) (calloc(nxny1, sizeof(float)));
    float *mask_s = (float *) (calloc(nxny1, sizeof(float)));
    
    // Above rowfil, all interpolated data
    
    for (int row = 0; row < rowfil; row++) {
        for (int col = 0; col < xsz; col++) {
            mask_s[row * xsz + col] = 1.0;
            mask_n[(rowfil1 - 1 - row) * xsz + col] = 1.0;
        }
    }
    
    // Between rowfil and rowfil1, a linear scale
    
    float d_row = 1.0 / (rowfil1 + 1.0 - rowfil);
    for (int row = rowfil; row < rowfil1; row++) {
        for (int col = 0; col < xsz; col++) {
            mask_s[row * xsz + col] = d_row * (rowfil1 - row);
            mask_n[(rowfil1 - 1 - row) * xsz + col] = d_row * (rowfil1 - row);
        }
    }
    
//   for (int row = 0; row < rowfil1; row++) printf("%d, %f, %f\n", row, mask_n[row*xsz+xsz-1], mask_s[row*xsz+xsz-1]);
    
    // Do north

#ifndef DEVELOP
    if (opt->nflag) {
        
        // Fill in north
 
        
        if (interpPoldb(0, rowfil, map, att, opt)) {    // map now filled down to rowfil
            FREE_ARR(map);
            FREE_ARR(img_n); FREE_ARR(img_s);
            FREE_ARR(mask_n); FREE_ARR(mask_s);
            return 1;
        }
        
        // Smooth north, smoothed filling included (same as spat_2d)
        
        synop2pole(0, att->sinlat, img_n, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        fitimage(img_n, imgDims[0], imgDims[1], opt->order, 1);
        pole2synop(0, att->sinlat, img_n, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        
        // Recombine with data
        
        int count0 = xysz - nxny1, count1 = xysz;
        for (int count = count0; count < count1; count++) {
            if (isnan(outMap[count])) outMap[count] = 0;        // otherwise nan*0=nan
            outMap[count] = map[count] * mask_n[count - count0] +
                            outMap[count] * (1. - mask_n[count - count0]);
        }
        
    }
    
    // Do south
    
    if (opt->sflag) {
        
        // Fill in south
        
        if (interpPoldb(1, rowfil, map, att, opt)) {        // map now filled down to rowfil
            FREE_ARR(map);
            FREE_ARR(img_n); FREE_ARR(img_s);
            FREE_ARR(mask_n); FREE_ARR(mask_s);
            return 1;
        }
        
        // Smooth south, smoothed filling included (same as spat_2d)
        
        synop2pole(1, att->sinlat, img_s, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        fitimage(img_s, imgDims[0], imgDims[1], opt->order, 1);
        pole2synop(1, att->sinlat, img_s, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        
        // Recombine with data, map is smoothed, outmap is original
        
        int count0 = 0, count1 = nxny1;
        for (int count = count0; count < count1; count++) {
            if (isnan(outMap[count])) outMap[count] = 0;        // otherwise nan*0=nan
            outMap[count] = map[count] * mask_s[count - count0] +
                            outMap[count] * (1. - mask_s[count - count0]);
        }
    
    }
#endif

    //
    
    FREE_ARR(map);
    FREE_ARR(img_n); FREE_ARR(img_s);
    FREE_ARR(mask_n); FREE_ARR(mask_s);
    
    return 0;
    
}

// =============================================

/* Method 2: 2D spatial smoothing only */

int spat_2d(float *outMap, struct attrib *att, struct option *opt)
{
    
    // Copy to double
    
    int xsz = att->mapDims[0], ysz = att->mapDims[1], xysz = xsz * ysz;
    double *map = (double *) (calloc(xysz, sizeof(double)));
    for (int count = 0; count < xysz; count++) {
        map[count] = outMap[count];
    }
    
    // Projected image
    
    int imgDims[2] = {IMGDIM, IMGDIM}, imgDim2 = imgDims[0] * imgDims[1];
    double *img_n = (double *) (calloc(imgDim2, sizeof(double)));
    double *img_s = (double *) (calloc(imgDim2, sizeof(double)));
    
    // Create mask to transition from smoothing filling (above latfil)
    // to data (above latfil)
    
    int rowfil = 0, rowfil1 = 0;
    while (rowfil < (ysz / 2)) {
        if (fabs(att->lats[rowfil++] / DTOR) < opt->latfil) break;
    }
    while (rowfil1 < (ysz / 2)) {
        if (fabs(att->lats[rowfil1++] / DTOR) < opt->latfil1) break;
    }
    
    printf("rowfil=%d, rowfil1=%d\n", rowfil, rowfil1);
    
    int nxny1 = xsz * rowfil1, nxny = xsz * rowfil;
    float *mask_n = (float *) (calloc(nxny1, sizeof(float)));
    float *mask_s = (float *) (calloc(nxny1, sizeof(float)));
    
    // Above rowfil, all interpolated data
    
    for (int row = 0; row < rowfil; row++) {
        for (int col = 0; col < xsz; col++) {
            mask_s[row * xsz + col] = 1.0;
            mask_n[(rowfil1 - 1 - row) * xsz + col] = 1.0;
        }
    }
    
    // Between rowfil and rowfil1, a linear scale
    
    float d_row = 1.0 / (rowfil1 + 1.0 - rowfil);
    for (int row = rowfil; row < rowfil1; row++) {
 //       printf("z=%f\n",d_row * (rowfil1 - row));
        for (int col = 0; col < xsz; col++) {
            mask_s[row * xsz + col] = d_row * (rowfil1 - row);
            mask_n[(rowfil1 - 1 - row) * xsz + col] = d_row * (rowfil1 - row);
        }
    }

    // North
    
#ifndef DEVELOP
    if (opt->nflag) {
        
        // In polar_remap.c, mapping map to image
        synop2pole(0, att->sinlat, img_n, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        fitimage(img_n, imgDims[0], imgDims[1], opt->order, 1);
        pole2synop(0, att->sinlat, img_n, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        
        // Recombine with data
        
        int count0 = xysz - nxny1, count1 = xysz;
        for (int count = count0; count < count1; count++) {
            if (isnan(outMap[count])) outMap[count] = 0;        // otherwise nan*0=nan
            outMap[count] = map[count] * mask_n[count - count0] +
                            outMap[count] * (1. - mask_n[count - count0]);
        }
        
    }
    
    // South
    
    if (opt->sflag) {

        // In polar_remap.c, mapping map to image
        synop2pole(1, att->sinlat, img_s, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        fitimage(img_s, imgDims[0], imgDims[1], opt->order, 1);
        pole2synop(1, att->sinlat, img_s, imgDims, map, att->mapDims,
                   att->lons, att->lats, (opt->lat0 * DTOR));
        
//        printf("map[0]=%f, outmap[0]=%f\n", map[0], outMap[0]);
        
        // Recombine with data, map is smoothed, outmap is original
        
        int count0 = 0, count1 = nxny1;
        for (int count = count0; count < count1; count++) {
            if (isnan(outMap[count])) outMap[count] = 0;        // otherwise nan*0=nan
            outMap[count] = map[count] * mask_s[count - count0] +
                            outMap[count] * (1. - mask_s[count - count0]);
        }
        
//        printf("map[0]=%f, outmap[0]=%f\n", map[0], outMap[0]);
        
    }
#endif
    
    //
    
    FREE_ARR(map);
    FREE_ARR(img_n); FREE_ARR(img_s);
    
    return 0;
    
}

// =============================================

/* Interpolate polar field using time series of polar field */

int interpPoldb(int ns, int rowfil, double *map, struct attrib *att, struct option *opt)
{
    
    int status = 0;
    
    // interpolation or extrapolation
    
    // Treat extrapolation as interpolation
    // For now, backward extrapolation uses the first two years
    // and a linear extrapolation
    // Forward extrapolation uses the last year, and scales
    // the mean to match the lower latitudes
    
    int *year = (ns) ? att->s_year : att->n_year;
    int extrap = (ns) ? att->s_extrap : att->n_extrap;
    
    // Read in records, already checked
    // Get necessary keywords, arrays, etc.
    
    char inQuery[100];
    DRMS_RecordSet_t *inRS0 = NULL, *inRS1 = NULL;
    DRMS_Record_t *inRec0 = NULL, *inRec1 = NULL;
    DRMS_Segment_t *inSeg0 = NULL, *inSeg1 = NULL;
    DRMS_Array_t *inArray0 = NULL, *inArray1 = NULL;
    
    double crpix1, crval1, cdelt1, crpix2, crval2, cdelt2;
    int xsz0, ysz0, xsz1, ysz1;
    double x, y;
    
    double carr_time0, carr_time1;
    double lonc0, lonc1;
    float scale0, scale1;
    double *data0, *data1;
    double *lons_p0 = NULL, *lats_p0 = NULL;        // lat/lon of polar data
    double *lons_p1 = NULL, *lats_p1 = NULL;
    
    snprintf(inQuery, 100, "%s[%4d][%d]", poldb_name, year[0], ns);     // already checked
    inRS0 = drms_open_records(drms_env, inQuery, &status);
//    printf("%s, %d, %d\n", inQuery, status, inRS0->n);
    if (status || inRS0->n < 1) {
        SHOW("  Open polar record error, ");
        return 1;
    }
    inRec0 = inRS0->records[0];
    //
    crpix1 = drms_getkey_double(inRec0, "CRPIX1", &status);
    crval1 = drms_getkey_double(inRec0, "CRVAL1", &status);
    cdelt1 = drms_getkey_double(inRec0, "CDELT1", &status);
    crpix2 = drms_getkey_double(inRec0, "CRPIX2", &status);
    crval2 = drms_getkey_double(inRec0, "CRVAL2", &status);
    cdelt2 = drms_getkey_double(inRec0, "CDELT2", &status);
    carr_time0 = drms_getkey_double(inRec0, "CARRTIME", &status);       // should be same as crpix1
    //
    scale0 = drms_getkey_double(inRec0, "SCALE", &status);
    inSeg0 = drms_segment_lookup(inRec0, "data");
    xsz0 = inSeg0->axis[0]; ysz0 = inSeg0->axis[1];
    inArray0 = drms_segment_read(inSeg0, DRMS_TYPE_DOUBLE, &status);    // read as double
    if (status) {
        SHOW("  Read polar data array error, ");
        return 1;
    }
    data0 = (double *) inArray0->data;
    //
    lons_p0 = (double *) (calloc(xsz0, sizeof(double)));
    lats_p0 = (double *) (calloc(ysz0, sizeof(double)));
    lonc0 = (ceil(carr_time0 / 360.) * 360. - carr_time0) * DTOR;
    for (int col = 0; col < xsz0; col++) {
        lons_p0[col] = (col + 1 - crpix1) * fabs(cdelt1) * DTOR + lonc0;
    }
    for (int row = 0; row < ysz0; row++) {
        y = (row - crpix2 + 1) * cdelt2 + crval2;
        lats_p0[row] = (att->sinlat) ? asin(y) : (y * DTOR);
    }
    
    if (extrap != 1) {      // not forward extrapolation
        snprintf(inQuery, 100, "%s[%4d][%d]", poldb_name, year[1], ns);     // already checked
        inRS1 = drms_open_records(drms_env, inQuery, &status);
        if (status) {
            SHOW("  Open polar record error, ");
            return 1;
        }
        inRec1 = inRS1->records[0];
        //
        crpix1 = drms_getkey_double(inRec1, "CRPIX1", &status);
        crval1 = drms_getkey_double(inRec1, "CRVAL1", &status);
        cdelt1 = drms_getkey_double(inRec1, "CDELT1", &status);
        crpix2 = drms_getkey_double(inRec1, "CRPIX2", &status);
        crval2 = drms_getkey_double(inRec1, "CRVAL2", &status);
        cdelt2 = drms_getkey_double(inRec1, "CDELT2", &status);
        carr_time1 = drms_getkey_double(inRec1, "CARRTIME", &status);       // should be same as crpix1
        //
        scale1 = drms_getkey_double(inRec1, "SCALE", &status);
        inSeg1 = drms_segment_lookup(inRec1, "data");
        xsz1 = inSeg1->axis[0]; ysz1 = inSeg1->axis[1];
        inArray1 = drms_segment_read(inSeg1, DRMS_TYPE_DOUBLE, &status);    // read as double 
        if (status) {
            SHOW("  Read polar data array error, ");
            if (inArray0) drms_free_array(inArray0);
            return 1;
        }
        data1 = (double *) inArray1->data;
        //
        lonc1 = (ceil(carr_time1 / 360.) * 360. - carr_time1) * DTOR;
    } else {        // forward extrapolation
        float mean62 = drms_getkey_float(inRec0, "MEAN62", &status);     // Mean field 62 & above
        float mean75 = drms_getkey_float(inRec0, "MEAN75", &status);     // Mean field 75 & above
        scale1 = scale0 * (mean62 / mean75 * K_TRANS);          // Scale field above 75
        data1 = data0;
        xsz1 = xsz0; ysz1 = ysz0;
        carr_time1 = carr_time0 + CARRLONINYEAR;     // assuming one year later
        lonc1 = lonc0;      // assuming same longitude (do not agree with carr_time1)
    }
    
    //
    lons_p1 = (double *) (calloc(xsz1, sizeof(double)));
    lats_p1 = (double *) (calloc(ysz1, sizeof(double)));
    for (int col = 0; col < xsz1; col++) {
        lons_p1[col] = (col + 1 - crpix1) * fabs(cdelt1) * DTOR + lonc1;
    }
    for (int row = 0; row < ysz1; row++) {
        y = (row - crpix2 + 1) * cdelt2 + crval2;
        lats_p1[row] = (att->sinlat) ? asin(y) : (y * DTOR);
    }
    
    // Interpolation, directly put values back to map
           
    double b0, b1;      // field
    printf("t=%f, t0=%f, t1=%f\n", att->carr_time, carr_time0, carr_time1);
    printf("scale0=%f, scale1=%f\n", scale0, scale1);
    
    double kt = (att->carr_time - carr_time0) / (carr_time1 - carr_time0);
    
    int row0, row1;     // rows to be filled
    if (ns) {
        row0 = 0; row1 = rowfil;
    } else {
        row0 = att->mapDims[1] - rowfil; row1 = att->mapDims[1];        // used as < row1
    }
    
    printf("kt=%f\n", kt);
    
    double x0, x1, y0, y1;
    double lon;
    for (int row = row0; row < row1; row++) {
        y0 = get_index(att->sinlat, ysz0, lats_p0, att->lats[row]);
        y1 = get_index(att->sinlat, ysz1, lats_p1, att->lats[row]);
        for (int col = 0; col < att->mapDims[0]; col++) {
            //
            lon = att->lons[col];
            while (lon < lons_p0[0]) {lon += TWOPI;}
            while (lon >= (lons_p0[0] + TWOPI)) {lon -= TWOPI;}
            x0 = get_index(0, xsz0, lons_p0, lon);
            //
            lon = att->lons[col];
            while (lon < lons_p1[0]) {lon += TWOPI;}
            while (lon >= (lons_p1[0] + TWOPI)) {lon -= TWOPI;}
            x1 = get_index(0, xsz1, lons_p1, lon);
            //
//            if (row == row0 && col % 100 == 0) {printf("x0=%f, x1=%f\n", x0, x1);}
            //
            b0 = linintd(data0, xsz0, ysz0, x0, y0) * scale0;
            b1 = linintd(data1, xsz1, ysz1, x1, y1) * scale1;
            //
            map[row * att->mapDims[0] + col] = b0 * (1. - kt) + b1 * kt;
        }
    }
    
    //
    
    drms_free_array(inArray0);
    drms_close_records(inRS0, DRMS_FREE_RECORD);
    if (extrap != 1) {
        drms_free_array(inArray1);
        drms_close_records(inRS1, DRMS_FREE_RECORD);
    }
    
    FREE_ARR(lons_p0); FREE_ARR(lats_p0);
    FREE_ARR(lons_p1); FREE_ARR(lats_p1);
    
    return 0;
    
}

// =============================================

/* Write output */

int writeOutput(char *outQuery, DRMS_Record_t *inRec, float **outMap,
                struct attrib *att, struct option *opt)
{
    int status = 0;
    
    DRMS_Record_t *outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        SHOW("  Output record not created, record skipped.\n");
        return 1;
    }
    
    DRMS_Segment_t *outSeg = drms_segment_lookupnum(outRec, 0);
    DRMS_Array_t *outArray = NULL;
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, att->mapDims, (*outMap), &status);
    if (status) {
        SHOW("  Creating output array failed, ");
        if (!outArray) {
            drms_free_array(outArray);
            *outMap = NULL;
        }
        drms_close_record(outRec, DRMS_FREE_RECORD);
        return 1;
    }
    
    outSeg->axis[0] = outArray->axis[0];
    outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;                // compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) {
        SHOW("  Writing output array failed, ");
        if (!outArray) drms_free_array(outArray);
        *outMap = NULL;
        drms_close_record(outRec, DRMS_FREE_RECORD);
        return 1;
    }
    
    // Keywords
    
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);    // copy all relevant keys?
    
    drms_copykey(outRec, inRec, "CAR_ROT");
    drms_copykey(outRec, inRec, "T_OBS");
    drms_copykey(outRec, inRec, "T_REC");
    drms_setkey_int(outRec, "AMB", opt->amb);
    
    set_statistics(outSeg, outArray, 1);        // Jan 18 2018
    
    switch (opt->method) {
        default:
        case (0):           // keys for TEMP_SPAT
            drms_setkey_string(outRec, "METHOD", "TEMP_SPAT");
            drms_setkey_int(outRec, "N_EXTRAP", att->n_extrap);
            drms_setkey_int(outRec, "S_EXTRAP", att->s_extrap);
            drms_setkey_double(outRec, "LAT0", opt->lat0);
            drms_setkey_double(outRec, "LATFIL", opt->latfil);
            drms_setkey_double(outRec, "LATFIL1", opt->latfil1);
            drms_setkey_int(outRec, "NMAX", opt->order);
            break;
        case (1):           // keys for SPAT_2D
            drms_setkey_string(outRec, "METHOD", "SPAT_2D");
            drms_setkey_int(outRec, "NMAX", opt->order);
            break;
    }
    
    char timebuf[1024];
    double UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
    drms_setkey_string(outRec, "DATE", timebuf);
    
    //
    
    drms_free_array(outArray);
    drms_close_record(outRec, DRMS_INSERT_RECORD);
    
    return 0;
    
}
