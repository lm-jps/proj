/*
 *  multi_cut.c
 *
 *  This module creates cutout from several full disk series per Yang's request:
 *      hmi.S_720s:     I0-5, Q0-5, U0-5, V0-5
 *      hmi.Ic_720s:    continuum
 *      hmi.Ic_noLimbDark_720s:     continuum
 *      hmi.V_720s:     Dopplergram
 *      Derived maps:   lon/lat/mu
 *  Track with differential rotation rate.
 *
 *  Input:
 *      tstr:   Time range, to be combined with data series names for query
 *      arid:   AR identifier, can be any integer number
 *      tref:   Reference time
 *
 *  Output:
 *      out:    Output data series name
 *
 *  Optional input:
 *      s/ic/icFlat/v:  Stokes, continnum, flattened continuum, Doppler series
 *      lonref/latref:  Reference Stonyhurst lon/lat at tref (default 0/0)
 *      cols/rows:  Cutout image size (default 500/500)
 *      c:      Set flag to track with Carrington rate (default no)
 *      noaa_ar/harpnum:    NOAA AR number/HARP number if known
 *
 *  Adapted from cgem_cutout.c
 *
 *  Author:
 *      Xudong Sun
 *
 *  Version
 *      v0.0    Apr 29 2020
 *
 *  Example calls:
 *      > multi_cut tstr="2011.02.15_12:00_TAI/24m" arid="11158" "out=hmi_test.multi_cut_720s" "tref=2011.02.14_12:00:00_TAI" "cols=960" "rows=960" "lonref=5.6" "latref=-20.4" "noaa_ar=11158" "harpnum=377"
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
#include "cartography.c"

// /cvs/JSOC/proj/lev1.5_hmi/libs/lev15/rotcoef_file.txt
double diffrot[3] = {2.71390, -0.4050, -0.4220};

#define PI              (M_PI)
#define RADSINDEG       (PI/180.)
#define RAD2ARCSEC      (648000./M_PI)
#define SECINDAY        (86400.)
#define FOURK           (4096)
#define FOURK2          (16777216)
#define DTTHRESH        (SECINDAY*7.)

#define MAXTDIFF        (10.)       // Max allowed T_REC diff in sec
#define MAXSTRLEN       (256)

#define ARRLENGTH(ARR)  (sizeof(ARR) / sizeof(ARR[0]))

// Macros
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified   "Not Specified"
#define dpath           "/home/jsoc/cvs/Development/JSOC"

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
    TIME t_rec;
};

// Mapping information
struct patchInfo {
    double lonref, latref;
    int cols, rows;
    TIME tref;
    int arid, noaa_ar, harpnum;
};

/* Cutout segment names, input identical to output */
char *sSegs[] = {"I0", "Q0", "U0", "V0", "I1", "Q1", "U1", "V1",
                 "I2", "Q2", "U2", "V2", "I3", "Q3", "U3", "V3",
                 "I4", "Q4", "U4", "V4", "I5", "Q5", "U5", "V5"};
char icSegIn[] = "continuum", icSegOut[] = "continuum";             // continuum segment name
char icFlatSegIn[] = "continuum", icFlatSegOut[] = "continuum0";  // flattened continuum segment name
char vSegIn[] = "Dopplergram", vSegOut[] = "Dopplergram";           // Doppler segment name
char lonSegOut[] = "lon", latSegOut[] = "lat", muSegOut[] = "mu";
char *segUnits[] = {"DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S",
                    "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S",
                    "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S", "DN/S",     // s
                    "DN/S",     // ic
                    " ",        // icFlat, normalized
                    "m/s",      // v
                    "degree", "degree", " "     // lon/lat/mu
                   };

/* ====================================================== */

/* Check inputs */
int checkInputs(DRMS_RecordSet_t *sRS, DRMS_RecordSet_t *icRS, DRMS_RecordSet_t *icFlatRS, DRMS_RecordSet_t *vRS);

/* Create Cutout record */
int createCutRecord(DRMS_Record_t *sRec, DRMS_Record_t *icRec, DRMS_Record_t *icFlatRec, DRMS_Record_t *vRec,
                    DRMS_Record_t *outRec, struct patchInfo *pInfo);

/* Determine location of cutout */
int findCoord(DRMS_Record_t *inRec, struct patchInfo *pInfo, double *diff_a, int *ll, int *ur);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Get cutout and write segment */
int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec,
                int *ll, int *ur, char *outSegName, char *inSegName);

/* Compute and write lon/lat/mu for cutout */
int writellm(DRMS_Record_t *outRec, DRMS_Record_t *inRec, int *ll, int *ur);

/* Set all keywords, no error checking for now */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct patchInfo *pInfo, int *ll);

/* ====================================================== */

char *module_name = "multi_cut";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "tstr", kNotSpecified, "CGEM dataset ID."},
    {ARG_INT, "arid", kNotSpecified, "AR ID."},
    {ARG_STRING, "out", kNotSpecified, "Output cutout series."},
    {ARG_STRING, "s", "hmi.S_720s", "Input Stokes series, name only."},
    {ARG_STRING, "ic", "hmi.Ic_720s", "Input continuum series, name only."},
    {ARG_STRING, "icFlat", "hmi.Ic_NoLimbDark_720s", "Input flattened continuum serie, name only."},
    {ARG_STRING, "v", "hmi.V_720s", "Input Doppler series, name only."},
    {ARG_INT, "noaa_ar", "-1", "Set NOAA number if positive, otherwise equivalent to CGEMNUM."},
    {ARG_INT, "harpnum", "-1", "Set HARP number if positive."},
    {ARG_STRING, "tref", kNotSpecified, "Reference time."},
    {ARG_FLOAT, "lonref", "0", "Reference patch center Stonyhurst lon, in deg."},
    {ARG_FLOAT, "latref", "0", "Reference patch center Stonyhurst lat, in deg."},
    {ARG_INT, "cols", "500", "Columns of output cutout."},
    {ARG_INT, "rows", "500", "Rows of output cutout."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    char *sDSName = (char *) params_get_str(&cmdparams, "s");
    char *icDSName = (char *) params_get_str(&cmdparams, "ic");
    char *icFlatDSName = (char *) params_get_str(&cmdparams, "icFlat");
    char *vDSName = (char *) params_get_str(&cmdparams, "v");
    char *outQuery = (char *) params_get_str(&cmdparams, "out");
    
    /* Get arguments */
    
    char *tstr = (char *) params_get_str(&cmdparams, "tstr");
    struct patchInfo pInfo;
    pInfo.lonref = params_get_float(&cmdparams, "lonref") * RADSINDEG;
    pInfo.latref = params_get_float(&cmdparams, "latref") * RADSINDEG;
    pInfo.cols = params_get_int(&cmdparams, "cols");
    pInfo.rows = params_get_int(&cmdparams, "rows");
    pInfo.tref = params_get_time(&cmdparams, "tref");
    pInfo.arid = params_get_int(&cmdparams, "arid");
    pInfo.noaa_ar = params_get_int(&cmdparams, "noaa_ar");
    pInfo.harpnum = params_get_int(&cmdparams, "harpnum");
    
    /* Obtain query strings */
    
    char sQuery[MAXSTRLEN], icQuery[MAXSTRLEN], icFlatQuery[MAXSTRLEN], vQuery[MAXSTRLEN];
//    sprint_time(tstr, trec, "TAI", 0);
    snprintf(sQuery, MAXSTRLEN, "%s[%s]\n", sDSName, tstr);
    snprintf(icQuery, MAXSTRLEN, "%s[%s]\n", icDSName, tstr);
    snprintf(icFlatQuery, MAXSTRLEN, "%s[%s]\n", icFlatDSName, tstr);
    snprintf(vQuery, MAXSTRLEN, "%s[%s]\n", vDSName, tstr);
    SHOW(sQuery); SHOW(icQuery); SHOW(icFlatQuery); SHOW(vQuery);
    
    /* Input data */
    
    DRMS_RecordSet_t *sRS = drms_open_records(drms_env, sQuery, &status);        // Input data series
    int nrecs_s = sRS->n;
    if (status || nrecs_s == 0 || !sRS) DIE("Input Stokes data series error.");
    
    DRMS_RecordSet_t *icRS = drms_open_records(drms_env, icQuery, &status);        // Input data series
    int nrecs_ic = icRS->n;
    if (status || nrecs_ic == 0 || !icRS) DIE("Input continuum data series error.");
    
    DRMS_RecordSet_t *icFlatRS = drms_open_records(drms_env, icFlatQuery, &status);        // Input data series
    int nrecs_icFlat = icFlatRS->n;
    if (status || nrecs_icFlat == 0 || !icFlatRS) DIE("Input flattened continuum data series error.");
    
    DRMS_RecordSet_t *vRS = drms_open_records(drms_env, vQuery, &status);        // Input data series
    int nrecs_v = vRS->n;
    if (status || nrecs_v == 0 || !vRS) DIE("Input Doppler data series error.");
    
    /* Check if all T_RECs are in order and match */
    
    if (checkInputs(sRS, icRS, icFlatRS, vRS)) {
        drms_close_records(sRS, DRMS_FREE_RECORD);
        drms_close_records(icRS, DRMS_FREE_RECORD);
        drms_close_records(icFlatRS, DRMS_FREE_RECORD);
        drms_close_records(vRS, DRMS_FREE_RECORD);
        DIE("Input records do not match.");
    }
    
    printf("Found %d groups of frames for [%s].\n", nrecs_s, tstr); fflush(stdout);
    
    /* Start */
    
    for (int irec = 0; irec < nrecs_s; irec++) {
        
        status = 0;
        
        printf("==============\nProcessing frame #%d of %d:\n", irec + 1, nrecs_s); fflush(stdout);
        
        /* Records at work */
        
        DRMS_Record_t *sRec = sRS->records[irec];
        DRMS_Record_t *icRec = icRS->records[irec];
        DRMS_Record_t *icFlatRec = icFlatRS->records[irec];
        DRMS_Record_t *vRec = vRS->records[irec];
        
        /* Create Cutout record */
        
        DRMS_Record_t *outRec = NULL;
        outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
        if (status) {        // if failed
            printf("Creating cutout failed, frame #%d skipped.\n", irec + 1);
            if (!outRec) {
                drms_close_record(outRec, DRMS_FREE_RECORD);
            }
            continue;
        }

        if (createCutRecord(sRec, icRec, icFlatRec, vRec,       // DO THE WORK!!!
                            outRec, &pInfo)) {
            printf("Creating cutout failed, frame #%d skipped.\n", irec + 1);
            drms_close_record(outRec, DRMS_FREE_RECORD);
            continue;
        }
        
        /* Done */
        
        drms_close_record(outRec, DRMS_INSERT_RECORD);
        
        printf("Frame #%d done.\n", irec + 1);
        
    }
    
    /* Clean up */
    
    drms_close_records(sRS, DRMS_FREE_RECORD);
    drms_close_records(icRS, DRMS_FREE_RECORD);
    drms_close_records(icFlatRS, DRMS_FREE_RECORD);
    drms_close_records(vRS, DRMS_FREE_RECORD);
    
    return 0;
    
}   // DoIt


// ===================================================================
// ===================================================================
// ===================================================================

/*
 * Pair records in four input series
 * T_REC must match in order
 *
 */

int checkInputs(DRMS_RecordSet_t *sRS, DRMS_RecordSet_t *icRS, DRMS_RecordSet_t *icFlatRS, DRMS_RecordSet_t *vRS)
{
    
    int status = 0;
    
    int nrecs_s = sRS->n, nrecs_ic = icRS->n, nrecs_icFlat = icFlatRS->n, nrecs_v = vRS->n;
    
    if (!((nrecs_s == nrecs_ic) && (nrecs_ic == nrecs_icFlat) && (nrecs_icFlat == nrecs_v))) {
        return 1;
    }
    
    for (int irec = 0; irec < nrecs_s; irec++) {
        TIME trec_s = drms_getkey_time(sRS->records[irec], "T_REC", &status);
        if (status) return 1;
        TIME trec_ic = drms_getkey_time(icRS->records[irec], "T_REC", &status);
        if (status || fabs(trec_s - trec_ic) > MAXTDIFF) return 1;
        TIME trec_icFlat = drms_getkey_time(icFlatRS->records[irec], "T_REC", &status);
        if (status || fabs(trec_s - trec_icFlat) > MAXTDIFF) return 1;
        TIME trec_v = drms_getkey_time(vRS->records[irec], "T_REC", &status);
        if (status || fabs(trec_s - trec_v) > MAXTDIFF) return 1;
    }
    
    return 0;
    
}


/*
 * Create Cutout record: top level subroutine
 * Do the loops on segments and set the keywords here
 * Work is done in writeCutout routine below
 *
 */

int createCutRecord(DRMS_Record_t *sRec, DRMS_Record_t *icRec, DRMS_Record_t *icFlatRec, DRMS_Record_t *vRec,
                    DRMS_Record_t *outRec, struct patchInfo *pInfo)
{
    
    int status = 0;
    
    // Obtain dimension of cutout in full disk
    
    int ll[2], ur[2];       // lower left and upper right coord
    if (findCoord(sRec, pInfo, diffrot, ll, ur)) {
        printf("Coordinate unreasonable\n");
        return 1;
    }
    printf("ll=%d,%d\n", ll[0], ll[1]);
    
    // Coutout Stokes
    
    int iSeg, nSegs = ARRLENGTH(sSegs);
    for (iSeg = 0; iSeg < nSegs; iSeg++) {
        if (writeCutout(outRec, sRec, ll, ur, sSegs[iSeg], sSegs[iSeg])) {
            printf("Stokes cutout fails for %s.\n", sSegs[iSeg]);
            return 1;   // failed
        }
    }
    SHOW("Stokes cutout done.\n");
    
    // Cutout Ic
    
    if (writeCutout(outRec, icRec, ll, ur, icSegOut, icSegIn)) {
        printf("Continuum cutout fails.\n");
        return 1;   // failed
    }
    SHOW("Continuum cutout done.\n");
    
    // Cutout flattened Ic
    
    if (writeCutout(outRec, icFlatRec, ll, ur, icFlatSegOut, icFlatSegIn)) {
        printf("Flattened continuum cutout fails.\n");
        return 1;   // failed
    }
    SHOW("Flattened continuum cutout done.\n");
    
    // Cutout flattened Ic
    
    if (writeCutout(outRec, vRec, ll, ur, vSegOut, vSegIn)) {
        printf("Doppler cutout fails.\n");
        return 1;   // failed
    }
    SHOW("Doppler cutout done.\n");
    
    // Compute & write out lon/lat/mu
    
    if (writellm(outRec, sRec, ll, ur)) {
        printf("lon/lat/mu output fails.\n");
        return 1;
    }
    SHOW("lon/lat/mu output done.\n");
    
    // Keywords & Links
    
    setKeys(outRec, sRec, pInfo, ll);      // set keywords
    SHOW("Keyword set done.\n");
    
    return 0;
    
}

/*
 * Determine the location of the patch
 *
 */

int findCoord(DRMS_Record_t *inRec, struct patchInfo *pInfo, double *diff_a, int *ll, int *ur)
{
    
    int status = 0;
    
    // Get ephemeris
    
    struct ephemeris ephem;
    if (getEphemeris(inRec, &ephem)) {
        SHOW("Get ephemeris error\n");
        return 1;
    }
    
    // Rotate
    // Differential rotation rate in urad/s
    
    double lat_c = pInfo->latref;
    double difr = diff_a[0] + diff_a[1] * pow(sin(lat_c),2) + diff_a[2] * pow(sin(lat_c),4);
    double dt = drms_getkey_time(inRec, "T_REC", &status) - pInfo->tref;
    double lon_c = pInfo->lonref + dt * difr * 1.0e-6;      // urad/s to rad
    printf("lon_c=%f, lat_c=%f\n", lon_c/RADSINDEG, lat_c/RADSINDEG);
    
    double xc, yc;        // patch center pixel in 4K image
    if (sphere2img (lat_c, lon_c, ephem.disk_latc, 0., &xc, &yc, // disk_lonc=0 (Stonyhurst)
                    ephem.disk_xc, ephem.disk_yc, ephem.rSun,      // in pixel
                    ephem.pa, 0., 0., 0., 0.)) {
        SHOW("Patch center on far side\n");
        return 1;
    }
    
    // Location
    
    ll[0] = round(xc - pInfo->cols / 2.0);
    ll[1] = round(yc - pInfo->rows / 2.0);
    ur[0] = ll[0] + pInfo->cols - 1;
    ur[1] = ll[1] + pInfo->rows - 1;
    
    return 0;
    
}


/*
 * Fetch ephemeris info from a DRMS record
 * No error checking for now
 *
 */

int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem)
{
    
    int status = 0;
    
    float crota2 = drms_getkey_float(inRec, "CROTA2", &status);    // rotation
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
    ephem->disk_xc = PIX_X(0.0,0.0) - 1.0;        // Center of disk in pixel, starting at 0
    ephem->disk_yc = PIX_Y(0.0,0.0) - 1.0;
    
    float dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
    float rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
    if (status) rSun_ref = 6.96e8;
    
    ephem->asd = asin(rSun_ref/dSun);
    ephem->rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
    
    ephem->t_rec = drms_getkey_time(inRec, "T_REC", &status);
    
    return 0;
    
}

/*
 * Get cutout and write segment
 *
 */

int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec,
                int *ll, int *ur, char *outSegName, char *inSegName)
{
    
    int status = 0;
    
    DRMS_Segment_t *inSeg = NULL, *outSeg = NULL;
    DRMS_Array_t *cutoutArray = NULL;
    
    inSeg = drms_segment_lookup(inRec, inSegName);
    if (!inSeg) return 1;
    
    /* Geometry info */
    
    int nx, ny, nxny;        // lower-left and upper right coords
    nx = ur[0] - ll[0] + 1;
    ny = ur[1] - ll[1] + 1;
    nxny = nx * ny;
    
    /* Read */
    
    if (inSeg->axis[0] == FOURK && inSeg->axis[1] == FOURK) {        // for full disk ones
        cutoutArray = drms_segment_readslice(inSeg, DRMS_TYPE_DOUBLE, ll, ur, &status);
        if (status) return 1;
    } else {
        return 1;
    }
    
    /* Write out */
    
    outSeg = drms_segment_lookup(outRec, outSegName);
    if (!outSeg) return 1;
    outSeg->axis[0] = cutoutArray->axis[0];
    outSeg->axis[1] = cutoutArray->axis[1];
    cutoutArray->israw = 0;        // always compressed
    cutoutArray->bzero = outSeg->bzero;
    cutoutArray->bscale = outSeg->bscale;        // Same as inArray's
    status = drms_segment_write(outSeg, cutoutArray, 0);
    drms_free_array(cutoutArray);
    if (status) return 1;
    
    //    printf("%s done\n", outSegName);
    
    return 0;
    
}

/*
 * Compute and write lon/lat/mu for cutout
 * Use cartography, all info from inRec, ll, ur
 *
 */
int writellm(DRMS_Record_t *outRec, DRMS_Record_t *inRec, int *ll, int *ur)
{

    int status = 0;
    
    // Check segments first
    
    DRMS_Segment_t *outSeg_lon = NULL, *outSeg_lat = NULL, *outSeg_mu = NULL;
    outSeg_lon = drms_segment_lookup(outRec, lonSegOut);
    outSeg_lat = drms_segment_lookup(outRec, latSegOut);
    outSeg_mu = drms_segment_lookup(outRec, muSegOut);
    if (outSeg_lon == NULL || outSeg_lat == NULL || outSeg_mu == NULL) {
        return 1;
    }
    
    // Geometry
    
    struct ephemeris ephem;
    getEphemeris(inRec, &ephem);
    
    // Arrays
    
    int nrows = ur[1] - ll[1] + 1;
    int ncols = ur[0] - ll[0] + 1;
    int npix = nrows * ncols;
    
    double *lon = (double *) (malloc(npix * sizeof(double)));
    double *lat = (double *) (malloc(npix * sizeof(double)));
    double *mu = (double *) (malloc(npix * sizeof(double)));
    
    // Compute
    
    int idx;    // index
    double x, y;
    double rho_t, lat_t, lon_t, sinlat_t, coslat_t, sig_t, mu_t, chi_t;
    
    for (int row = 0; row < nrows; row++) {
        for (int col = 0; col < ncols; col++) {
            idx = row * ncols + col;
            x = (col + ll[0] - ephem.disk_xc) / ephem.rSun;
            y = (row + ll[1] - ephem.disk_yc) / ephem.rSun;
            if (img2sphere(x, y, ephem.asd, ephem.disk_latc, 0.0, ephem.pa,
                           &rho_t, &lat_t, &lon_t, &sinlat_t, &coslat_t, &sig_t, &mu_t, &chi_t)) {
                lon[idx] = DRMS_MISSING_DOUBLE;
                lat[idx] = DRMS_MISSING_DOUBLE;
                mu[idx] = DRMS_MISSING_DOUBLE;
                continue;
            }
            lon[idx] = lon_t / RADSINDEG;
            lat[idx] = lat_t / RADSINDEG;
            mu[idx] = mu_t;
        }
    }
    
    // Write out
    
    int dims[2] = {ncols, nrows};
    double *outArrs[3] = {lon, lat, mu};
    DRMS_Segment_t *outSegs[3] = {outSeg_lon, outSeg_lat, outSeg_mu};
    for (int iSeg = 0; iSeg < 3; iSeg++) {
        DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims, outArrs[iSeg], &status);
        if (status) {
            for (int i = iSeg; i < 3; i++) { free(outArrs[i]); }
            return 1;
        }
        outSegs[iSeg]->axis[0] = outArray->axis[0];
        outSegs[iSeg]->axis[1] = outArray->axis[1];
        outArray->israw = 0;        // always compressed
        outArray->bzero = outSegs[iSeg]->bzero;
        outArray->bscale = outSegs[iSeg]->bscale;
        status = drms_segment_write(outSegs[iSeg], outArray, 0);
        if (status) {
            for (int i = iSeg; i < 3; i++) { free(outArrs[i]); }
            return 1;
        }
        drms_free_array(outArray);      // success
    }
    
    return 0;
}



/* Set keywords */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct patchInfo *pInfo, int *ll)
{
    
    int status = 0;
    
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
    drms_setkey_int(outRec, "ARID",  pInfo->arid);
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    drms_setkey_int(outRec, "NOAA_AR",  pInfo->noaa_ar);
    drms_setkey_int(outRec, "HARPNUM",  pInfo->harpnum);
    
    drms_setkey_double(outRec, "DIFF_A0", diffrot[0]);
    drms_setkey_double(outRec, "DIFF_A2", diffrot[1]);
    drms_setkey_double(outRec, "DIFF_A4", diffrot[2]);
    
    // Geometry
    
    struct ephemeris ephem;
    getEphemeris(inRec, &ephem);
    // Defined as disk center's pixel address wrt lower-left of cutout
    drms_setkey_float(outRec, "CRPIX1", ephem.disk_xc - ll[0] + 1.);
    drms_setkey_float(outRec, "CRPIX2", ephem.disk_yc - ll[1] + 1.);
    drms_setkey_float(outRec, "IMCRPIX1", ephem.disk_xc + 1.);
    drms_setkey_float(outRec, "IMCRPIX2", ephem.disk_yc + 1.);
    // Always 0.
    drms_setkey_float(outRec, "CRVAL1", 0);
    drms_setkey_float(outRec, "CRVAL2", 0);
    drms_setkey_float(outRec, "IMCRVAL1", 0);
    drms_setkey_float(outRec, "IMCRVAL2", 0);
    
    // Unit
    int nSegs = ARRLENGTH(segUnits);
    for (int iSeg = 0; iSeg < nSegs; iSeg++) {
        DRMS_Segment_t *outSeg = NULL;
        outSeg = drms_segment_lookupnum(outRec, iSeg);      // ordered
        if (!outSeg) continue;
        // Set Bunit
        char bunit_xxx[MAXSTRLEN];
        sprintf(bunit_xxx, "BUNIT_%03d", iSeg);
        //printf("%s, %s\n", bunit_xxx, CutBunits[iSeg]);
        drms_setkey_string(outRec, bunit_xxx, segUnits[iSeg]);
    }
    
}
