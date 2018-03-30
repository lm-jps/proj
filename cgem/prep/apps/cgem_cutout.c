/*
 *  cgem_cutout.c   
 *
 *	This module creates cutout of b vectors and v
 *	User needs to specify the detailed geometry of the cutout
 *	including the center of the reference time in lat/lon and frame size
 *	The module takes full disk b (disambuated)
 *	and differentially track for the specified time interval
 *	Adapted from sharp.c
 *
 *	Author:
 *		Xudong Sun
 *
 *	Version:
 *              v0.0 Jan 07 2016
 *	Notes:
 *		v0.0
 *
 *	Example Calls:
 *      > cgem_cutout "b=hmi.B_720s[2011.02.15_12:00/2h]" "dop=hmi_test.doppcal_720s[2011.02.15_12:00/2h]" "out=hmi_test.bvcutout_720s" "tref=2011.02.14_12:00:00_TAI" "cols=960" "rows=960" "lonref=5.6" "latref=-20.4" "cgemnum=11158"
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

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)
#define DTTHRESH        (518400.)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

#define DISAMB_AZI		1
#define DOPPCAL         0       // Calibrate Doppler here

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
    int cgemnum;
    int proj;
};

/* Cutout segment names, input identical to output */
char *BSegs[] = {"vlos_mag", "inclination", "azimuth", "field", "disambig", "conf_disambig"};
char *CutSegs[] = {"vlos_mag", "inclination", "azimuth", "field", "disambig", "conf_disambig"};
char *CutBunits[] = {"cm/s", "degree", "degree", "Mx/cm^2", " ", " "};

/* ========================================================================================================== */

/* Get all input data series */
int getInputRS(DRMS_RecordSet_t **bRS_ptr, DRMS_RecordSet_t **dopRS_ptr,
			   char *bQuery, char *dopQuery, struct patchInfo *pInfo);

/* Create Cutout record */
int createCutRecord(DRMS_Record_t *bRec, DRMS_Record_t *dopRec,
					DRMS_Record_t *outRec, struct patchInfo *pInfo);

/* Determine location of cutout */
int findCoord(DRMS_Record_t *inRec, struct patchInfo *pInfo, int *ll, int *ur);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Get cutout and write segment */
int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec,
                int *ll, int *ur, double doppbias, char *SegName);

/* Set all keywords, no error checking for now */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct patchInfo *pInfo, int *ll);

/* Compute Dopper correction due to spacecraft and rotation; same as cgem_prep */
int getDoppCorr_patch(DRMS_Record_t *inRec, int *ll, int *ur, float *doppCorr);

/* ========================================================================================================== */

char *module_name = "cgem_cutout";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "b", kNotSpecified, "Input B series."},
//	{ARG_STRING, "dop", kNotSpecified, "Input Doppler series."},
	{ARG_STRING, "dop", kNotSpecified, "Data sereis containing Doppler bias correction"},
	{ARG_STRING, "out", kNotSpecified, "Output Sharp cutout series."},
	{ARG_INT, "cgemnum", "-1", "CGEM dataset ID."},
	{ARG_FLOAT, "lonref", "0", "Reference patch center Stonyhurst lon, in deg."},
	{ARG_FLOAT, "latref", "0", "Reference patch center Stonyhurst lat, in deg."},
	{ARG_INT, "cols", "500", "Columns of output cutout."},
	{ARG_INT, "rows", "500", "Rows of output cutout."},
	{ARG_STRING, "tref", kNotSpecified, "Reference time."},
	{ARG_END}
};

int DoIt(void)
{

    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    char *bQuery = NULL, *dopQuery = NULL;      // Input query
    char *outQuery = NULL;                      // Output query
    
    bQuery = (char *) params_get_str(&cmdparams, "b");
    dopQuery = (char *) params_get_str(&cmdparams, "dop");
    outQuery = (char *) params_get_str(&cmdparams, "out");
    
    /* Get arguments */
    
    struct patchInfo pInfo;
    pInfo.lonref = params_get_float(&cmdparams, "lonref") * RADSINDEG;
    pInfo.latref = params_get_float(&cmdparams, "latref") * RADSINDEG;
    pInfo.cols = params_get_int(&cmdparams, "cols");
    pInfo.rows = params_get_int(&cmdparams, "rows");
    pInfo.tref = params_get_time(&cmdparams, "tref");
    pInfo.cgemnum = params_get_int(&cmdparams, "cgemnum");
    
    /* Read input data, check everything */

    DRMS_RecordSet_t *bRS = NULL, *dopRS = NULL;        // Input data series
    // tref set to the first image t_rec if greater than 6 days ahead
    if (getInputRS(&bRS, &dopRS, bQuery, dopQuery, &pInfo))
        DIE("Input data error.");
    int nrecs = bRS->n;

    /* Start */
    
    printf("==============\nStart. %d image(s) in total.\n", nrecs);
    
    for (int irec = 0; irec < nrecs; irec++) {

        /* Records at work */
        
        DRMS_Record_t *bRec = NULL, *dopRec = NULL;
        
        bRec = bRS->records[irec];
		dopRec = dopRS->records[irec];         // already checked in getInputRS
        
        TIME trec = drms_getkey_time(bRec, "T_REC", &status);
        
        /*
        char tstr[100], queryStr[100];
        sprint_time(tstr, trec, "TAI", 0);
        snprintf(queryStr, 100, "%d%s[%s]\n", irec, bRec->seriesinfo->seriesname, tstr);
        SHOW(queryStr);
         */
        
        /* Create Cutout record */
        
        DRMS_Record_t *outRec = NULL;
        outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
        if (status) {		// if failed
            printf("Creating cutout failed, image #%d skipped.\n", irec);
            continue;
        }
        
        if (createCutRecord(bRec, dopRec, outRec, &pInfo)) {		// do the work
            printf("Creating cutout failed, image #%d skipped.\n", irec);
            drms_close_record(outRec, DRMS_FREE_RECORD);
            continue;
        }
        
        /* Done */
        
        drms_close_record(outRec, DRMS_INSERT_RECORD);
        
        printf("Image #%d done.\n", irec);
        
    }
    
    /* Clean up */
    
    drms_close_records(bRS, DRMS_FREE_RECORD);
//    drms_close_records(dopRS, DRMS_FREE_RECORD);
    
	return 0;
	
}	// DoIt


// ===================================================================
// ===================================================================
// ===================================================================


/*
 * Get input data series, including b and dop
 * Need all records to match, otherwise quit
 *
 */

int getInputRS(DRMS_RecordSet_t **bRS_ptr, DRMS_RecordSet_t **dopRS_ptr, 
			   char *bQuery, char *dopQuery, struct patchInfo *pInfo)
{
    
    int status = 0;
    
    *bRS_ptr = drms_open_records(drms_env, bQuery, &status);
    if (status || (*bRS_ptr)->n == 0) return 1;
    
    *dopRS_ptr = drms_open_records(drms_env, dopQuery, &status);
    if (status || (*dopRS_ptr)->n == 0) return 1;
    
    // Compare T_REC for b and dop, must be identical

    
    if ((*bRS_ptr)->n != (*dopRS_ptr)->n) return 1;
    int nrecs = (*bRS_ptr)->n;
    
    DRMS_Record_t *bRec_t = NULL, *dopRec_t = NULL;		// temporary recs for utility
    
    for (int i = 0; i < nrecs; i++) {
        bRec_t = (*bRS_ptr)->records[i];
        dopRec_t = (*dopRS_ptr)->records[i];
        if (drms_getkey_time(bRec_t, "T_REC", &status) !=
             drms_getkey_time(dopRec_t, "T_REC", &status))
            return 1;
    }
    
    // Correct t_ref
    
    TIME t0 = drms_getkey_time((*bRS_ptr)->records[0], "T_REC", &status);
    if (fabs(pInfo->tref - t0) > DTTHRESH)
        pInfo->tref = t0;
    
    return 0;
    
}


/*
 * Create Cutout record: top level subroutine
 * Do the loops on segments and set the keywords here
 * Work is done in writeCutout routine below
 *
 */

int createCutRecord(DRMS_Record_t *bRec, DRMS_Record_t *dopRec,
					DRMS_Record_t *outRec, struct patchInfo *pInfo)
{
    
    int status = 0;
    
    // Obtain dimension of cutout in full disk
    
    int ll[2], ur[2];
    if (findCoord(bRec, pInfo, ll, ur)) {
        printf("Coordinate unreasonable\n");
        return 1;
    }
    printf("ll=%d,%d\n", ll[0], ll[1]);
    
    // Cutout Doppler
    
    /*
    if (writeCutout(outRec, dopRec, ll, ur, "Dopplergram")) {
        printf("Doppler cutout failed\n");
        return 1;
    }
    printf("Dopplergram cutout done.\n");
    */
    
    // Coutout B
    
    int iSeg, nSegs = ARRLENGTH(BSegs);
    double doppbias = drms_getkey_double(dopRec, "DOPPBIAS", &status);
    
    for (iSeg = 0; iSeg < nSegs; iSeg++) {
        if (writeCutout(outRec, bRec, ll, ur, doppbias, BSegs[iSeg])) {
            printf("B cutout fails for %s\n", BSegs[iSeg]);
            break;
        }
    }
    if (iSeg != nSegs) return 1;		// failed
    printf("Vector B cutout done.\n");
    
    // Keywords & Links
    
    setKeys(outRec, bRec, pInfo, ll);      // set keywords
	drms_copykey(outRec, dopRec, "DOPPBIAS");     // doppler correction
    
    return 0;
    
}


/*
 * Determine the location of the patch
 *
 */

int findCoord(DRMS_Record_t *inRec, struct patchInfo *pInfo, int *ll, int *ur)
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
    // proj/lev1.5_hmi/libs/lev15/rotcoef_file.txt
    
    double lat_c = pInfo->latref;
    double difr = 2.7139 - 0.405 * pow(sin(lat_c),2) - 0.422 * pow(sin(lat_c),4);
    double dt = drms_getkey_time(inRec, "T_REC", &status) - pInfo->tref;
    double lon_c = pInfo->lonref + dt * difr * 1.0e-6;      // urad/s to rad
    printf("lon_c=%f\n", lon_c/RADSINDEG);
    
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
    ephem->rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
    
    ephem->t_rec = drms_getkey_time(inRec, "T_REC", &status);
    
    return 0;
    
}


/*
 * Get cutout and write segment
 * DISAMB_AZI == 1 to apply disambiguation to azimuth
 *
 */

int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec,
                int *ll, int *ur, double doppbias, char *SegName)
{
    
    int status = 0;
    
    DRMS_Segment_t *inSeg = NULL, *outSeg = NULL;
    DRMS_Array_t *cutoutArray = NULL;
    
    inSeg = drms_segment_lookup(inRec, SegName);
    if (!inSeg) return 1;
    
    /* Geometry info */
    
    int nx, ny, nxny;		// lower-left and upper right coords
    nx = ur[0] - ll[0] + 1;
    ny = ur[1] - ll[1] + 1;
    nxny = nx * ny;
    
    /* Read */
    
    if (inSeg->axis[0] == FOURK && inSeg->axis[1] == FOURK) {		// for full disk ones
        cutoutArray = drms_segment_readslice(inSeg, DRMS_TYPE_DOUBLE, ll, ur, &status);
        if (status) return 1;
    } else {
        return 1;
    }
    
    /* Disambig */
    
#if DISAMB_AZI
    if (!strcmp(SegName, "azimuth")) {
        DRMS_Segment_t *disambSeg = NULL;
        disambSeg = drms_segment_lookup(inRec, "disambig");
        if (!disambSeg) {drms_free_array(cutoutArray); return 1;}
        DRMS_Array_t *disambArray = drms_segment_readslice(disambSeg, DRMS_TYPE_CHAR, ll, ur, &status);
        if (status) return 1;
        double *azimuth = (double *) cutoutArray->data;
        char *disamb = (char *) disambArray->data;
        for (int n = 0; n < nxny; n++) {
            // Feb 12 2014, use bit #2 for full disk, lowest bit for patch
            if ((disamb[n] / 4) % 2) azimuth[n] += 180.;
        }
        drms_free_array(disambArray);
    }
#endif
    
#if DOPPCAL
    if (!strcmp(SegName, "vlos_mag")) {
        // Doppler correction due to spacecraft v and diff rot
        printf("doppbias=%f\n", doppbias);
        float *doppCorr = (float *) (malloc(nxny * sizeof(float)));
        getDoppCorr_patch(inRec, ll, ur, doppCorr);
        double *dop = (double *) cutoutArray->data;
        for (int i = 0; i < nxny; i++) {
            dop[i] = dop[i] - doppCorr[i] - doppbias;                // cm/s
        }
        free(doppCorr);
    }
#endif
    
    /* Write out */
    
    outSeg = drms_segment_lookup(outRec, SegName);
    if (!outSeg) return 1;
    outSeg->axis[0] = cutoutArray->axis[0];
    outSeg->axis[1] = cutoutArray->axis[1];
    cutoutArray->israw = 0;		// always compressed
    cutoutArray->bzero = outSeg->bzero;
    cutoutArray->bscale = outSeg->bscale;		// Same as inArray's
    status = drms_segment_write(outSeg, cutoutArray, 0);
    drms_free_array(cutoutArray);
    if (status) return 1;
    
//    printf("%s done\n", SegName);
    
    return 0;

}


/* Set keywords */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct patchInfo *pInfo, int *ll)
{
    
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
    drms_setkey_int(outRec, "CGEMNUM",  pInfo->cgemnum);
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);

    // Geometry
    
    int status = 0;
    
    float disk_xc = drms_getkey_float(outRec, "CRPIX1", &status);
    float disk_yc = drms_getkey_float(outRec, "CRPIX2", &status);
    // Defined as disk center's pixel address wrt lower-left of cutout
    drms_setkey_float(outRec, "CRPIX1", disk_xc - ll[0] + 1.);
    drms_setkey_float(outRec, "CRPIX2", disk_yc - ll[1] + 1.);
    drms_setkey_float(outRec, "IMCRPIX1", disk_xc + 1.);
    drms_setkey_float(outRec, "IMCRPIX2", disk_yc + 1.);
    // Always 0.
    drms_setkey_float(outRec, "CRVAL1", 0);
    drms_setkey_float(outRec, "CRVAL2", 0);
    drms_setkey_float(outRec, "IMCRVAL1", 0);
    drms_setkey_float(outRec, "IMCRVAL2", 0);
    
    // Jan 2 2014 XS
    int nSeg = ARRLENGTH(CutSegs);
    for (int iSeg = 0; iSeg < nSeg; iSeg++) {
        DRMS_Segment_t *outSeg = NULL;
        outSeg = drms_segment_lookup(outRec, CutSegs[iSeg]);
        if (!outSeg) continue;
        // Set Bunit
        char bunit_xxx[20];
        sprintf(bunit_xxx, "BUNIT_%03d", iSeg);
        //printf("%s, %s\n", bunit_xxx, CutBunits[iSeg]);
        drms_setkey_string(outRec, bunit_xxx, CutBunits[iSeg]);
    }
    
}


// =============================================

/*
 * Compute Dopper correction due to spacecraft and rotation
 * Adapted from Brian's code doppcal_estimate.f90, in cm/s
 * Based on Phil's email: need to be careful about the sign of V
 * Doppler velocity has redshift as pos (toward Sun center)
 * OBS_VR/VW/VN is in heliocetric coordinate, so away from Sun center is pos
 * To perform correction, one needs to project OBS_VX onto the LOS vector
 * in *HELIOCENTRIC* coordniate (away from Sun is positive)
 * Then simply subtract this number from the Dopplergram
 *
 */

int getDoppCorr_patch(DRMS_Record_t *inRec, int *ll, int *ur, float *doppCorr)
{
    
    int status = 0;
    
    // Info
    
    int nx = ur[0] - ll[0] + 1, ny = ur[1] - ll[1] + 1, nxny = nx * ny;
    
    struct ephemeris ephem;
    if (getEphemeris(inRec, &ephem)) {
        SHOW("Mapping: get ephemeris error\n");
        return 1;
    }
    
    // Lon/lat
    
    double *lon_nobp = (double *) (malloc(nxny * sizeof(double)));      // assuming b=0, p=0, in rad
    double *lat_nobp = (double *) (malloc(nxny * sizeof(double)));
    double *lon = (double *) (malloc(nxny * sizeof(double)));
    double *lat = (double *) (malloc(nxny * sizeof(double)));
    double x, y, lon_t, lat_t;
    double rho, sinlat, coslat, sig, mu, chi;
    int ind;
    
    printf("xc=%f, yc=%f, lonc=%f\nnx=%d, ny=%d\n", ephem.disk_xc, ephem.disk_yc, ephem.disk_lonc, nx, ny);
    printf("ll[0]=%d, ll[1]=%d\nasd=%f, pa=%f, rSun=%f\n", ll[0], ll[1], ephem.asd, ephem.pa, ephem.rSun);
    
    for (int row = 0; row < ny; row++) {
        for (int col = 0; col < nx; col++) {
            
            ind = row * nx + col;
            x = (ll[0] + col - ephem.disk_xc) / ephem.rSun;       // normalized pixel address
            y = (ll[1] + row - ephem.disk_yc) / ephem.rSun;       // normalized pixel address
            
            img2sphere (x, y, ephem.asd, 0.0, 0.0, 0.0,
                        &rho, &lat_t, &lon_t, &sinlat, &coslat, &sig, &mu, &chi);
            lon_nobp[ind] = lon_t; lat_nobp[ind] = lat_t;
            
            img2sphere (x, y, ephem.asd, ephem.disk_latc, 0.0, ephem.pa,
                        &rho, &lat_t, &lon_t, &sinlat, &coslat, &sig, &mu, &chi);
            lon[ind] = lon_t; lat[ind] = lat_t;
            
        }
    }
    
    // Correction, following Brian's method
    
    double vr = drms_getkey_double(inRec, "OBS_VR", &status) * 100.;    // m/s => cm/s
    double vw = drms_getkey_double(inRec, "OBS_VW", &status) * 100.;
    double vn = drms_getkey_double(inRec, "OBS_VN", &status) * 100.;
    double k = sin(ephem.asd);      // ~(1/215.)
    
    double crota2 = ephem.pa * (-1.);        // In rad already
    double obsv_x = vw * cos(crota2) + vn * sin(crota2);
    double obsv_y = - vw * sin(crota2) + vn * cos(crota2);
    
    double rsun_ref = 6.96e10;                // cm
    double coslatc = cos(ephem.disk_latc);
    double sinlon;
    double sinlat_nobp, sinlon_nobp, coslat_nobp;
    
    for (int row = 0; row < ny; row++) {
        for (int col = 0; col < nx; col++) {
            
            ind = row * nx + col;
            
            sinlat = sin(lat[ind]);
            sinlon = sin(lon[ind]);
            sinlat_nobp = sin(lat_nobp[ind]);
            sinlon_nobp = sin(lon_nobp[ind]);
            coslat_nobp = cos(lat_nobp[ind]);
            
            doppCorr[ind] = vr - obsv_y * sinlat_nobp * k
                            - obsv_x * sinlon_nobp * coslat_nobp * k;        // cm/s
            
            doppCorr[ind] += (rsun_ref * sinlon * coslatc * 1.e-6 *
                              (2.71390 - 0.405000 * pow(sinlat, 2) - 0.422000 * pow(sinlat, 4)));        // cm/s
            
        }
    }
    
    //
    
    free(lon_nobp); free(lat_nobp);     // clean up
    free(lon); free(lat);
    
    return 0;
    
}
