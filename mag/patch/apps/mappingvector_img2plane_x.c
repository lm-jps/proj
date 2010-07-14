/*
 * This jsoc module maps vector magnetic field from image coordinates to 
 * heliographic planer. The code uses Rick's cartography package for projection.
 * The code finds corresponding image coordinates for sample grid point on the plane
 * and interpolates for values on the image. Vector transformation is performed
 * after the interpolation. To supress aliasing we oversample the image before binning
 * it down to the final resolution, using Jesper's frebin code. One can also set the patch
 * center via providing the center Lon/Lat, or a reference data record which contains all
 * info to compute center Lon/Lat. Assuming RoI rotates at a Carrington rate this functionality
 * essentialy can track the RoI over disk passage.
 */

// Author:  Y. Liu, Apr 7 2010
// Adated:  X. Sun, Apr 8 2010, for patches with tracking functionality
// Modified:
//          X. Sun, Jun 30 2010, added prototypes for oversampling/rebin

// mappingvector_img2plane_x "in=su_xudong.test_patch_vec_a[2010.03.29_12:00:00_TAI][2]" "out=su_xudong.test_mercator_vec" "LONWIDTH=15.0" "LATWIDTH=9.0"
// mappingvector_img2plane_x "in=su_xudong.test_patch_vec_a[2010.03.29_12:12:00_TAI/1h][2]" "out=su_xudong.test_mercator_vec" "LONWIDTH=15.0" "LATWIDTH=9.0" "ref=[2010.03.29_12:12:00_TAI][2]"



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "astro.h"
#include "fstats.h"
#include "cartography.c"
#include "interpolate.c"
#include "img2helioVector.c"
#include "errorprop.c"

#define PI              (M_PI)
#define RADSINDEG       (PI/180)
#define RAD2ARCSEC      (648000. / M_PI)
#define SECINDAY		(86400.)

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)


/* ========================================================================================================== */



char *module_name = "mappingvector_img2plane_x";

ModuleArgs_t module_args[] =
{
   {ARG_STRING, "in", NULL, "Input data records."},
   {ARG_STRING, "out", NULL, "Output data series."},
   {ARG_STRING, "ref", "norec", "Reference record for alignment"},
   {ARG_STRING, "PROJECTION", "MERCATOR", ""},
   {ARG_FLOAT, 	"XUNITMAP", "0.03", "X resolution in mapped image"}, 	// pixel resolution in the mapping image, in degree.
   {ARG_FLOAT, 	"YUNITMAP", "0.03", "Y resolution in mapped image"}, 	// correspond to a resolution of 0.5 arc-sec/pixel.
   {ARG_FLOAT, 	"SAMPLE", "1", "Ratio of over-sampling"}, 				// the ratio of pixel sampling rate and final resolution, oversampling for anti-aliasing
   {ARG_FLOAT, 	"LONWIDTH", "0.", "X demension"}, 						// dimension, in degree
   {ARG_FLOAT, 	"LATWIDTH", "0.", "Y demension"},
   {ARG_FLOAT, 	"REFLON", "370.", "Reference lon for alignment"},		// reference lon/lat of patch center
   {ARG_FLOAT, 	"REFLAT", "100.", "Reference lat for alignment"},
   {ARG_INT, 	"COVAR", "0", ""},
   {ARG_FLAG, 	"e", "0", ""}, 											// set this flag to compute errors.
   {ARG_END}
};


int DoIt(void)
{
    int status = DRMS_SUCCESS;
    int nrecs, irec;
    
    DRMS_RecordSet_t *inRS, *outRS, *refRS = NULL;
    DRMS_Record_t *inRec, *outRec, *refRec = NULL;
    
    DRMS_Segment_t *inSeg;
    DRMS_Array_t *inArray_ambig;
    DRMS_Array_t *inArray_bTotal, *inArray_bAzim, *inArray_bIncl, *inArray_bFill;
    DRMS_Array_t *inArray_errbT, *inArray_errbAz, *inArray_errbIn, *inArray_errbF;
    DRMS_Array_t *inArray_bTbI, *inArray_bTbA, *inArray_bTbF, *inArray_bAbI, *inArray_bAbF, *inArray_bIbF;
    
    DRMS_Segment_t *outSeg;
    DRMS_Array_t *outArray;
    
    DRMS_Link_t *patchLink;
    
    char *inQuery, *outQuery, *refQuery;
    char *mapping;
    int projection, covar, doerror;
    float xunitMap, yunitMap, sample;
    float LonWidth, LatWidth;
    float ref_lon, ref_lat;
    
    char *ctype1, *ctype2;
    int useOrigWidth;
    TIME ref_trec, trec, dt;

    int xDim = 4096, yDim = 4096;
    float rSun, xCenter, yCenter;
    float lonc, latc, pa;
    float cdelt, rSun_ref, dSun;
    double maplatc, maplonc, maplatl, maplonl, maplatr, maplonr;
    double x, y, lon, lat, coslat, sinlat;
    double mu, rho, sig, chi, xi, zeta;
    double bx_tmp, by_tmp, bz_tmp;
    double bx_helio, by_helio, bz_helio;
    double asd;
    float xcoor, ycoor, xcen, ycen, x_halfwidth, y_halfwidth;
    float crvalx, crvaly, crpix1, crpix2, crota2, sina, cosa, pcx, pcy;      
    
    int xMap, yMap, xMap0, yMap0;
    float maplat_size, maplon_size;	// from HWIDTHn
    
    int llx, lly, bmx, bmy;		// lower left coordinate and dimension of bitmap
    
    float *bxHel, *byHel, *bzHel;		// remapped, oversampled
    float *bxHel0, *byHel0, *bzHel0;	// remapped, final
    float *bx, *by, *bz;		// full disk
    char *ambig;		// bitmap from disambiguation
    float *bTotal, *bIncl, *bAzim, *bFill;		// full disk
    float *errBx, *errBy, *errBz;		// remapped, oversampled
    float *errBx0, *errBy0, *errBz0;		// remapped, final
    float *errbT, *errbAz, *errbIn, *errbF;
    float *errbTbI, *errbTbA, *errbTbF, *errbAbI, *errbAbF, *errbIbF;
    float *zeroArray;
    double bxSigma2, bySigma2, bzSigma2;
    
    int outDims[2];
    
    /* Get parameters */
    
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    refQuery = (char *)params_get_str(&cmdparams, "ref");
    mapping = (char *)params_get_str(&cmdparams, "PROJECTION");

    xunitMap = cmdparams_get_float(&cmdparams, "XUNITMAP", &status);
    yunitMap = cmdparams_get_float(&cmdparams, "YUNITMAP", &status);
    sample = cmdparams_get_float(&cmdparams, "SAMPLE", &status);   
    if (sample < 1.) sample = 1.;	// Sampling rate

    LonWidth = cmdparams_get_float(&cmdparams, "LONWIDTH", &status);
    LatWidth = cmdparams_get_float(&cmdparams, "LATWIDTH", &status);
    ref_lon = cmdparams_get_float(&cmdparams, "REFLON", &status);
    ref_lat = cmdparams_get_float(&cmdparams, "REFLAT", &status);

    covar = cmdparams_get_int(&cmdparams, "COVAR", &status);
    doerror = (cmdparams_isflagset(&cmdparams, "e") != 0);

    
    /* Set projection method */
    
    if (strcmp(mapping, "RECTANGULAR") == 0 || strcmp(mapping, "Rectangular") == 0 ||
        strcmp(mapping, "rectangular") == 0) {
                        projection = 0; ctype1 = "HGLN-CAR"; ctype2 = "HGLT-CAR";}
    if (strcmp(mapping, "CASSINI") == 0 || strcmp(mapping, "Cassini") == 0 || 
        strcmp(mapping, "cassini") == 0) {
                        projection = 1; ctype1 = "HGLN-CAS"; ctype2 = "HGLT-CAS";}
    if (strcmp(mapping, "MERCATOR") == 0 || strcmp(mapping, "Mercator") == 0 ||
        strcmp(mapping, "mercator") == 0) {
                        projection = 2; ctype1 = "HGLN-MER"; ctype2 = "HGLT-MER";}
    if (strcmp(mapping, "CYLEQA") == 0 || strcmp(mapping, "Cyleqa") == 0 ||
        strcmp(mapping, "cyleqa") == 0) {
                        projection = 3; ctype1 = "HGLN-CYL"; ctype2 = "HGLT-CYL";}
    if (strcmp(mapping, "SINEQA") == 0 || strcmp(mapping, "Sineqa") == 0 ||
        strcmp(mapping, "sineqa") == 0) {
                        projection = 4; ctype1 = "HGLN-SIE"; ctype2 = "HGLT-SIE";}
    if (strcmp(mapping, "GNOMONIC") == 0 || strcmp(mapping, "Gnomonic") == 0 ||
        strcmp(mapping, "gnomonic") == 0) {
                        projection = 5; ctype1 = "HGLN-TAN"; ctype2 = "HGLT-TAN";}
    if (strcmp(mapping, "POSTEL") == 0 || strcmp(mapping, "Postel") == 0 ||
        strcmp(mapping, "postel") == 0) {
                        projection = 6; ctype1 = "HGLN-POS"; ctype2 = "HGLT-POS";}
    if (strcmp(mapping, "STEREOGRAPHIC") == 0 || strcmp(mapping, "Stereographic") == 0 ||
        strcmp(mapping, "stereographic") == 0) {
                        projection = 7; ctype1 = "HGLN-STE"; ctype2 = "HGLT-STE";}
    if (strcmp(mapping, "ORTHOGRAPHIC") == 0 || strcmp(mapping, "Orthographic") == 0 ||
        strcmp(mapping, "orthographic") == 0) {
                        projection = 8; ctype1 = "HGLN-SIN"; ctype2 = "HGLT-SIN";}
    if (strcmp(mapping, "LAMBERT") == 0 || strcmp(mapping, "Lambert") == 0 ||
        strcmp(mapping, "lambert") == 0) {
                        projection = 9; ctype1 = "HGLN-LAM"; ctype2 = "HGLT-LAM";}
                        
    /* Find the reference record for alignment, if not found, refRec = NULL */
     
    refRS = drms_open_records(drms_env, refQuery, &status);
    if (!status) {
        refRec = refRS->records[0]; // open patch data record
        ref_trec = drms_getkey_time(refRec, "T_REC", &status);
        
        dSun = drms_getkey_float(refRec, "DSUN_OBS", &status);
        rSun_ref = drms_getkey_float(refRec, "RSUN_REF", &status);
        if (status) rSun_ref = 6.96e8;
        cdelt = drms_getkey_float(refRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
        asd = asin(rSun_ref/dSun);
        rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;

        crvalx = drms_getkey_float(refRec, "IMCRVAL1", &status);	// center of solar disc in arcsec
        crvaly = drms_getkey_float(refRec, "IMCRVAL2", &status);
        crpix1 = drms_getkey_float(refRec, "IMCRPIX1", &status);	// disk center in ccd, original
        crpix2 = drms_getkey_float(refRec, "IMCRPIX2", &status);
        crota2 = drms_getkey_float(refRec, "CROTA2", &status);	// rotation
        sina = sin(crota2 * RADSINDEG); 
        cosa = cos(crota2 * RADSINDEG);

        xCenter = PIX_X(0.0,0.0) - 1.0;		// Center of disk, starting at 0
        yCenter = PIX_Y(0.0,0.0) - 1.0;

        latc = drms_getkey_float(refRec, "CRLT_OBS", &status) * RADSINDEG;
        lonc = drms_getkey_float(refRec, "CRLN_OBS", &status) * RADSINDEG;
        pa = - crota2 * RADSINDEG;

        pcx = drms_getkey_float(refRec, "CRPIX1", &status);
        pcy = drms_getkey_float(refRec, "CRPIX2", &status);
        xcoor = pcx - 1.;					// Center of patch, starting at 0
        ycoor = pcy - 1.;

        x = (xcoor - xCenter) / rSun;
        y = (ycoor - yCenter) / rSun;
        if (img2sphere (x, y, asd, latc, lonc, pa, &rho, &maplatc, &maplonc,
                        &sinlat, &coslat, &sig, &mu, &chi)) {
            refRec = NULL;
        } else {
            ref_lat = maplatc / RADSINDEG;
            ref_lon = maplonc / RADSINDEG;
        }

    }
    status = 0;
    if (refRS != NULL) drms_close_records(refRS, DRMS_FREE_RECORD);
                  
    /* Input data */
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;

    /* Output data */
    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output recordset not created");

    // For each record
    
    for (irec = 0; irec < nrecs; irec++)
    {
    
        /* Input data keywords */

        inRec = inRS->records[irec]; // open patch data record
        trec = drms_getkey_time(inRec, "T_REC", &status);
                
        dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
        rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
        if (status) rSun_ref = 6.96e8;
        cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
        asd = asin(rSun_ref/dSun);
        rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;

        crvalx = drms_getkey_float(inRec, "IMCRVAL1", &status);	// center of solar disc in arcsec
        crvaly = drms_getkey_float(inRec, "IMCRVAL2", &status);
        crpix1 = drms_getkey_float(inRec, "IMCRPIX1", &status);	// disk center in ccd, original
        crpix2 = drms_getkey_float(inRec, "IMCRPIX2", &status);
        crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation
        sina = sin(crota2 * RADSINDEG); 
        cosa = cos(crota2 * RADSINDEG);

        xCenter = PIX_X(0.0,0.0) - 1.0;		// Center of disk, starting at 0
        yCenter = PIX_Y(0.0,0.0) - 1.0;

        latc = drms_getkey_float(inRec, "CRLT_OBS", &status) * RADSINDEG;
        lonc = drms_getkey_float(inRec, "CRLN_OBS", &status) * RADSINDEG;
        pa = - crota2 * RADSINDEG;

        pcx = drms_getkey_float(inRec, "CRPIX1", &status);
        pcy = drms_getkey_float(inRec, "CRPIX2", &status);
        xcoor = pcx - 1.;					// Center of patch, starting at 0
        ycoor = pcy - 1.;

        x_halfwidth = drms_getkey_int(inRec, "HWIDTH1", &status);
        y_halfwidth = drms_getkey_int(inRec, "HWIDTH2", &status);
        
        /* Dimensions of bitmap */
        
        llx = (int)xcoor - (int)x_halfwidth;
        lly = (int)ycoor - (int)y_halfwidth;
        bmx = 2 * (int)x_halfwidth + 1;
        bmy = 2 * (int)y_halfwidth + 1;
        
        /* Find binding lon, lat (maplatc/maplonc, maplatl/maplonl, maplatr/maplonr)  */

        useOrigWidth = 1;		// Use patch width derived from HWIDTH until overridden by input parameters
        
        x = (xcoor - xCenter) / rSun;
        y = (ycoor - yCenter) / rSun;
        if (img2sphere (x, y, asd, latc, lonc, pa, &rho, &maplatc, &maplonc,
                        &sinlat, &coslat, &sig, &mu, &chi))
            DIE("problem finding map center");

        x = ((xcoor - x_halfwidth) - xCenter) / rSun;
        y = ((ycoor - y_halfwidth) - yCenter) / rSun;
        if (img2sphere (x, y, asd, latc, lonc, pa, &rho, &maplatl, &maplonl,
                    &sinlat, &coslat, &sig, &mu, &chi)) {
            printf("problem finding left boundary, use command line parameter\n");
            useOrigWidth = 0;
        }

        x = ((xcoor + x_halfwidth) - xCenter) / rSun;
        y = ((ycoor + y_halfwidth) - yCenter) / rSun;
        if (img2sphere (x, y, asd, latc, lonc, pa, &rho, &maplatr, &maplonr,
                    &sinlat, &coslat, &sig, &mu, &chi)) {
            printf("problem finding right boundary, use command line parameter\n");
            useOrigWidth = 0;
        }
        
        /* Find geometry of the remapped image
         * try HWIDTH first, overrides with input paramters if necessary
         * xMap/yMap is dimension; MapLat/MapLon is center lon/lat
         */
        
        if (useOrigWidth) {
            if (fabs(maplonr - maplonl) < M_PI) {
                maplon_size = fabs(maplonr - maplonl) / (xunitMap * RADSINDEG);
            } else {
                maplon_size = (2 * M_PI - fabs(maplonr - maplonl)) / (xunitMap * RADSINDEG);
            }
            maplat_size = fabs(maplatr - maplatl) / (yunitMap * RADSINDEG);
        }
        
        if (LonWidth > 1.E-5) {
            xMap0 = rint (LonWidth / xunitMap);
            xMap = rint (LonWidth / xunitMap * sample);
        } else {
            if (!useOrigWidth) DIE("no x width available\n");
            xMap0 = rint (maplon_size);
            xMap = rint (maplon_size * sample);
        }
        
        if (LatWidth > 1.E-5) {
            yMap0 = rint (LatWidth / yunitMap);
            yMap = rint (LatWidth / yunitMap * sample);
        } else {
            if (!useOrigWidth) DIE("no y width available\n");
            yMap0 = rint (maplat_size);
            yMap = rint (maplat_size * sample);
        }
        
        printf("xMap=%d, yMap=%d\n", xMap, yMap);
        
        /* Align, using reference record, or reference lat/lon */
        
        dt = trec - ref_trec;
        if (fabs(dt) < 7. * SECINDAY) {		// No more than 7 days allowed
            maplatc = ref_lat * RADSINDEG;
            maplonc = ref_lon * RADSINDEG;
            while (maplonc < 0.) maplonc += (2. * PI);
            while (maplonc >= 2. * PI) maplonc -= (2. * PI);
            printf("Use ref lon/lat: ");
        } else {
            printf("Use orig lon/lat: ");
        }
        
        printf("MapLon = %f, MapLat = %f\n", maplonc / RADSINDEG, maplatc / RADSINDEG);

        /* Allocate arrays */

        bx = (float *) malloc(xDim * yDim * sizeof(float));
        by = (float *) malloc(xDim * yDim * sizeof(float));
        bz = (float *) malloc(xDim * yDim * sizeof(float));

        bxHel = (float *) malloc(xMap * yMap * sizeof(float));
        byHel = (float *) malloc(xMap * yMap * sizeof(float));
        bzHel = (float *) malloc(xMap * yMap * sizeof(float));

        if (doerror) 
        {
        errBx = (float *) malloc(xMap * yMap * sizeof(float));
        errBy = (float *) malloc(xMap * yMap * sizeof(float));
        errBz = (float *) malloc(xMap * yMap * sizeof(float));
        }

        /* Read in segments */

        inSeg = drms_segment_lookup(inRec, "ambig_flag");
        inArray_ambig = drms_segment_read(inSeg, DRMS_TYPE_CHAR, &status);
        if (status)
          DIE("problem reading file");
        ambig = (char *)inArray_ambig->data;
        if (inArray_ambig->axis[0] != bmx || inArray_ambig->axis[1] != bmy) {
            printf("%d, %d, %d, %d\n", 
                inArray_ambig->axis[0], bmx, inArray_ambig->axis[1], bmy);
            fflush(stdout);
            DIE("wrong ambig dimension");
        }

        inSeg = drms_segment_lookup(inRec, "field");
        inArray_bTotal = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file");
        bTotal = (float *)inArray_bTotal->data;

        inSeg = drms_segment_lookup(inRec, "azimuth");
        inArray_bAzim = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file");
        bAzim = (float *)inArray_bAzim->data;

        inSeg = drms_segment_lookup(inRec, "inclination");
        inArray_bIncl = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file");
        bIncl = (float *)inArray_bIncl->data;

        inSeg = drms_segment_lookup(inRec, "alpha_mag");
        inArray_bFill = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file");
        bFill = (float *)inArray_bFill->data;

        if (doerror)
        {
        inSeg = drms_segment_lookup(inRec, "field_err");
        inArray_errbT = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file errbT");
        errbT = (float *)inArray_errbT->data;

        inSeg = drms_segment_lookup(inRec, "azimuth_err");
        inArray_errbAz = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file errbAz");
        errbAz = (float *)inArray_errbAz->data;

        inSeg = drms_segment_lookup(inRec, "inclination_err");
        inArray_errbIn = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file errbIn");
        errbIn = (float *)inArray_errbIn->data;

        inSeg = drms_segment_lookup(inRec, "alpha_err");
        inArray_errbF = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
          DIE("problem reading file errbF");
        errbF = (float *)inArray_errbF->data;

        /* Covariances generally not in use */
        if (!covar) {
            zeroArray = (float *)malloc(xDim * yDim * sizeof(float));
            errbTbI = errbTbA = errbTbF = errbAbI = errbAbF = errbIbF = zeroArray;
            for (int ii = 0; ii < xDim * yDim; ii++) zeroArray[ii] = 0.0;
        } else {
            inSeg = drms_segment_lookup(inRec, "field_Inclination_err");
            inArray_bTbI = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
              DIE("problem reading file errbTbI");
            errbTbI = (float *)inArray_bTbI->data;

            inSeg = drms_segment_lookup(inRec, "field_az_err");
            inArray_bTbA = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
              DIE("problem reading file errbTbA");
            errbTbA = (float *)inArray_bTbA->data;

            inSeg = drms_segment_lookup(inRec, "field_alpha_err");
            inArray_bTbF = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
              DIE("problem reading file errbTbF");
            errbTbF = (float *)inArray_bTbF->data;

            inSeg = drms_segment_lookup(inRec, "inclin_azimuth_err");
            inArray_bAbI = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
              DIE("problem reading file errbAbI");
            errbAbI = (float *)inArray_bAbI->data;

            inSeg = drms_segment_lookup(inRec, "azimuth_alpha_err");
            inArray_bAbF = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
              DIE("problem reading file errbAbF");
            errbAbF = (float *)inArray_bAbF->data;

            inSeg = drms_segment_lookup(inRec, "inclination_alpha_err");
            inArray_bIbF = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
              DIE("problem reading file errbIbF");
            errbIbF = (float *)inArray_bIbF->data;
        }	// covar
        }	// doerror
        
        /* Convert data to Bxyz */
        
        int kx, ky, kOff;
        char flag;		// bitmap flag
        int jy = 0, yOff = 0, iData = 0;
        int ix = 0;
        for (jy = 0; jy < yDim; jy++)
        {
            ix = 0;
            yOff = jy * xDim;
            ky = jy - lly;
            for (ix = 0; ix < xDim; ix++)
            {
                iData = yOff + ix;
                kx = ix - llx;
                
                /* zero azi pointing North, zero incl pointing out from sun */
                bx[iData] = - bTotal[iData] * bFill[iData] * sin(bIncl[iData] * RADSINDEG) * sin(bAzim[iData] * RADSINDEG);
                by[iData] = bTotal[iData] * bFill[iData] * sin(bIncl[iData] * RADSINDEG) * cos(bAzim[iData] * RADSINDEG);
                bz[iData] = bTotal[iData] * bFill[iData] * cos(bIncl[iData] * RADSINDEG);

                /* zero azi pointing North, zero incl pointing into sun */
/*
                bx[iData] = - bTotal[iData] * bFill[iData] * sin(bIncl[iData] * RADSINDEG) * sin(bAzim[iData] * RADSINDEG);
                by[iData] = bTotal[iData] * bFill[iData] * sin(bIncl[iData] * RADSINDEG) * cos(bAzim[iData] * RADSINDEG);
                bz[iData] = - bTotal[iData] * bFill[iData] * cos(bIncl[iData] * RADSINDEG);
*/
            
                /* Disambiguation */
                
                if (kx < 0 || kx >= bmx || ky < 0 || ky > bmy) {
                	continue;
                } else {
                	kOff = ky * bmx + kx;
                	if (ambig[kOff]) {		// 180
                		bx[iData] *= -1.; by[iData] *= -1.;
                	} 
                }
            }
        }
        
        /* Start interpolation */
        
        jy = 0; yOff = 0; iData = 0;
        for (jy = 0; jy < yMap; jy++)
        {
            ix = 0;
            yOff = jy * xMap;
            for (ix = 0; ix < xMap; ix++)
            {
                iData = yOff + ix;
                x = (ix + 0.5 - xMap/2.) * xunitMap;
                y = (jy + 0.5 - yMap/2.) * yunitMap;

                /* map grid [x, y] corresponds to the point [lon, lat] in the heliographic coordinates. 
                 * the [x, y] are in radians with respect of the center of the map [xcMap, ycMap].
                 * projection methods could be Mercator, Lambert, and many others. [maplonc, mapLatc]
                 * is the heliographic longitude and latitude of the map center. Both are in degree.    
                 */

                if (plane2sphere (x * RADSINDEG, y * RADSINDEG, maplatc, maplonc, 
                                  &lat, &lon, projection)) {
                    bxHel[iData] = DRMS_MISSING_FLOAT;
                    byHel[iData] = DRMS_MISSING_FLOAT;
                    bzHel[iData] = DRMS_MISSING_FLOAT;
                    if (doerror) 
                    {
                    errBx[iData] = bxSigma2;
                    errBy[iData] = bySigma2;
                    errBz[iData] = bzSigma2;
                    }
                    continue;
                }

                /* map the grid [lon, lat] in the heliographic coordinates to [xi, zeta], a point in the
                 * image coordinates. The image properties, xCenter, yCenter, rSun, pa, ecc and chi are given.
                 */

                if (sphere2img (lat, lon, latc, lonc, &xi, &zeta, xCenter/rSun, yCenter/rSun, 1.0,
                                pa, 1., 0., 0., 0.)) {
                    bxHel[iData] = DRMS_MISSING_FLOAT;
                    byHel[iData] = DRMS_MISSING_FLOAT;
                    bzHel[iData] = DRMS_MISSING_FLOAT;
                    if (doerror)
                    {
                    errBx[iData] = bxSigma2;
                    errBy[iData] = bySigma2;
                    errBz[iData] = bzSigma2;
                    }
                    continue;
                }
                xi *= rSun;
                zeta *= rSun;

                /* Interpolation is carried out here. The method is the cubic convolution. */

                bx_tmp = ccint2 (bx, xDim, yDim, xi, zeta);
                by_tmp = ccint2 (by, xDim, yDim, xi, zeta);
                bz_tmp = ccint2 (bz, xDim, yDim, xi, zeta);
                
                /* Projection */
                
                bx_helio = 0; by_helio = 0; bz_helio = 0;
                
                img2helioVector (bx_tmp, by_tmp, bz_tmp,
                                 &bx_helio, &by_helio, &bz_helio,
                                 lon, lat, lonc, latc, pa);

                bxHel[iData] = bx_helio;
                byHel[iData] = by_helio;
                bzHel[iData] = bz_helio;
                
                /* If an estimate of the error is needed, call the subroutine errorprop.c.
                 * The image location is xi (x-axis) and zeta (y-axis).
                 */

                if (doerror)
                {
                    if (errorprop (bTotal, bAzim, bIncl, bFill, 
                                   errbT, errbAz, errbIn, errbF, errbTbA, errbTbI, errbTbF, errbAbI, errbAbF, errbIbF, 
                                   lon, lat, lonc, latc, pa, xDim, yDim, xi, zeta, 
                                   &bxSigma2, &bySigma2, &bzSigma2)) 
                    {
                    errBx[iData] = DRMS_MISSING_FLOAT;
                    errBy[iData] = DRMS_MISSING_FLOAT;
                    errBz[iData] = DRMS_MISSING_FLOAT;
                    continue;
                    }
                    if (!(bxSigma2 >= 0.0 || bySigma2 >= 0.0 || bzSigma2 >= 0.0))    
                    {
                    errBx[iData] = DRMS_MISSING_FLOAT;
                    errBy[iData] = DRMS_MISSING_FLOAT;
                    errBz[iData] = DRMS_MISSING_FLOAT;
                    continue;
                    } 
                    errBx[iData] = sqrt(bxSigma2);
                    errBy[iData] = sqrt(bySigma2);
                    errBz[iData] = sqrt(bzSigma2);
                }

            }	// ix
        }	//jy
        
        /* Rebinning */
        
        bxHel0 = (float *) malloc(xMap0 * yMap0 * sizeof(float));
        byHel0 = (float *) malloc(xMap0 * yMap0 * sizeof(float));
        bzHel0 = (float *) malloc(xMap0 * yMap0 * sizeof(float));
        if (doerror) 
        {
        errBx0 = (float *) malloc(xMap0 * yMap0 * sizeof(float));
        errBy0 = (float *) malloc(xMap0 * yMap0 * sizeof(float));
        errBz0 = (float *) malloc(xMap0 * yMap0 * sizeof(float));
        }
        
        // Wrapper for Jesper's code
        /*
        frebin(bxHel, xMap, yMap, bxHel0, xMap0, yMap0);
        frebin(byHel, xMap, yMap, byHel0, xMap0, yMap0);
        frebin(bzHel, xMap, yMap, bzHel0, xMap0, yMap0);
        if (doerror)
        {
        frebin(errBx, xMap, yMap, errBx0, xMap0, yMap0);
        frebin(errBy, xMap, yMap, errBy0, xMap0, yMap0);
        frebin(errBz, xMap, yMap, errBz0, xMap0, yMap0);
        }
        */
        
        // This is for now
        for (int ii = 0; ii < xMap0 * yMap0; ii++) {
            bxHel0[ii] = bxHel[ii];
            byHel0[ii] = byHel[ii];
            bzHel0[ii] = bzHel[ii];
        }
        if (doerror)
        {
        for (int ii = 0; ii < xMap0 * yMap0; ii++) {
            errBx0[ii] = errBx[ii];
            errBy0[ii] = errBy[ii];
            errBz0[ii] = errBz[ii];
        }
        }

        /* Stats */

      
        /* Output */

        outDims[0] = xMap0; outDims[1] = yMap0;
        outRec = outRS->records[irec];     
        
        // write Bx as the first segment
        outSeg = drms_segment_lookup(outRec, "Bx");
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, bxHel0, &status);
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file");
        drms_free_array(outArray);

        // write By as the second segment
        outSeg = drms_segment_lookup(outRec, "By");
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, byHel0, &status);
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file");
        drms_free_array(outArray);

        // write Bz as the third segment
        outSeg = drms_segment_lookup(outRec, "Bz");
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, bzHel0, &status);
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file");
        drms_free_array(outArray);
        
        if (doerror)
        {
        // write errBx as the fourth segment
        outSeg = drms_segment_lookup(outRec, "errBx");
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, errBx0, &status);
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file errBx");
        drms_free_array(outArray);

        // write errBy as the fourth segment
        outSeg = drms_segment_lookup(outRec, "errBy");
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, errBy0, &status);
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file errBy");
        drms_free_array(outArray);

        // write errBz as the fourth segment
        outSeg = drms_segment_lookup(outRec, "errBz");
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, errBz0, &status);
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file errBz");
        drms_free_array(outArray);
        }  
        
        /* Clean up */
        
        drms_free_array(inArray_ambig);
        drms_free_array(inArray_bTotal); drms_free_array(inArray_bAzim);
        drms_free_array(inArray_bIncl); drms_free_array(inArray_bFill);
        free(bx); free(by); free(bz);
        free(bxHel); free(byHel); free(bzHel);
                
        if (doerror) {
            drms_free_array(inArray_errbT); drms_free_array(inArray_errbAz);
            drms_free_array(inArray_errbIn); drms_free_array(inArray_errbF);
            free(errBx); free(errBy); free(errBz);
            if (!covar) {
                free(zeroArray);
            } else {
                drms_free_array(inArray_bTbI); drms_free_array(inArray_bTbA);
                drms_free_array(inArray_bTbF); drms_free_array(inArray_bAbI);
                drms_free_array(inArray_bAbF); drms_free_array(inArray_bIbF);
            }
        }

        /* Keywords */
        
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
		drms_copykey(outRec, inRec, "T_REC");
        drms_copykey(outRec, inRec, "PNUM");
        drms_setkey_string(outRec, "BUNIT", "Gauss");
        drms_setkey_string(outRec, "CTYPE1", ctype1);
        drms_setkey_string(outRec, "CTYPE2", ctype2);
        drms_setkey_float(outRec, "CRPIX1", xMap / 2. + 0.5);  
        drms_setkey_float(outRec, "CRPIX2", yMap / 2. + 0.5);
        drms_setkey_float(outRec, "CDELT1", xunitMap);
        drms_setkey_float(outRec, "CDELT2", yunitMap);
        drms_setkey_float(outRec, "CRVAL1", maplonc / RADSINDEG);
        drms_setkey_float(outRec, "CRVAL2", maplatc / RADSINDEG);
        drms_setkey_string(outRec, "CUNIT1", "degree");
        drms_setkey_string(outRec, "CUNIT2", "degree");
        drms_setkey_string(outRec, "PROJECT", mapping);
        
        /* Link */
        
        patchLink = hcon_lookup_lower(&outRec->links, "PATCH");
        if (patchLink) {
            drms_setlink_static(outRec, "PATCH", inRec->recnum);
        }
        
    }		// irec
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);

    return 0;

}
