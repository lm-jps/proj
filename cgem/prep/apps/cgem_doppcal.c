/*
 *  cgem_doppcal.c
 *
 *  This module calculates covective blueshift using the method in
 *  Welsch, Fisher, & Sun (2013, ApJ, 765, 98)
 *  Fortran functions written by B. Welsch
 *  Wrapper written by X. Sun
 *
 *  Input:
 *      Vector field record containting:
 *          1. azimuth (not disambiguated)
 *          2. conf_disambig
 *          3. confid
 *          4. disambig
 *          5. field
 *          6. inclination
 *          7. vlos (uncorrected for S/C, diff. rot., or convective blueshift)
 *
 *  Output:
 *      Doppler velocity bias
 *      Dopplergram corrected for spacecraft motion and rotation
 *      Map of LOS PILs
 *      Map of radial field PILs
 *
 *  Example:
 *      % cgem_doppcal "in=hmi.B_720s[2011.02.15_12:00/1h]" "out=hmi_test.doppcal_720s"
 *
 *  Note:
 *      Probably useful to save vlos corrected for diff.rot and cov. blueshift
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "diffrot.h"        // double diff_rot[3];
#include "jsoc_main.h"

// Macros

#define PI              (M_PI)
#define FOURK           (4096)
#define FOURK2          (16777216)
#define RADSINDEG        (PI/180.)
#define RAD2ARCSEC        (648000./M_PI)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

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

// Segments

const char *inSegNames[] = {"azimuth", "conf_disambig", "confid_map", "disambig",
                            "field", "inclination", "vlos_mag"};
const char *outSegNames[] = {"vlos_mag_corr", "pil_bl", "pil_br"};

const int nInSegs = (ARRLENGTH(inSegNames));
const int nOutSegs = (ARRLENGTH(outSegNames));

// Ephemeris information

struct ephemeris {
    double crlt_obs, crln_obs;
    double obs_vr, obs_vn, obs_vw;
    double rSun_obs, rSun_ref, crota2;
    double crpix1, crpix2, cdelt;
    TIME t_rec;
};

// =====================================

/* Read input */
int getInput(DRMS_Record_t *inRec, double *data_in, struct ephemeris *ephem);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Write output */
int writeOutput(DRMS_Record_t *inRec, DRMS_Record_t *outRec, int writeImg, 
                float *data_out, float doppcal_bias, double a0, double a2, double a4);

/* Doppler correction */
extern void doppcal_estimate_ (int *naxis1, int *naxis2, int *naxis3i, int *naxis3o,
                               double *data_in, float *data_out,
                               double *crpix, double *obs_v, double *rsun_obs, double *rsun_ref, double *crlt_obs, double *crota2,
                               float *doppcal_bias,
                               double *max_ang_in, double *thresh_blos_in, double *pix_in,
                               double *thresh_bmag_in, int *rad_los_sep_in,
                               double *diff_a0_in, double *diff_a2_in, double *diff_a4_in);

// =====================================

char *module_name = "cgem_doppcal";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input data series."},
    {ARG_STRING, "out", kNotSpecified, "Output data series."},
    {ARG_FLAG, "w", "", "Set flag to write out images."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    char *inQuery = (char *) params_get_str(&cmdparams, "in");
    char *outQuery = (char *) params_get_str(&cmdparams, "out");
    int writeImg = params_isflagset(&cmdparams, "w");
    
    double a0 = diffrot[0], a2 = diffrot[1], a4 = diffrot[2];
//    printf("a0=%lf, a2=%lf, a4=%lf\n", a0, a2, a4);
    
    /* Input data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs == 0 || !inRS) {
        DIE("Input data series error");
    }
    
    /* Output data */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status || !outRS) {
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("Error in output series");
    }
    
    printf("==============\nStart processing. %d image(s) in total.\n", nrecs);
    
    /* Main loop */
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        printf("==============\nProcessing frame %d of %d...\n", irec + 1, nrecs);
        
        status = 0;
        struct ephemeris ephem;
        
        /* Get input data */
        
        printf("Reading input... ");
        DRMS_Record_t *inRec = inRS->records[irec];
        
        double *data_in = (double *) (malloc(FOURK2 * nInSegs * sizeof(double)));       // FD only
        
        status = getInput(inRec, data_in, &ephem);
        if (status) {
            SHOW("input read error, record skipped.\n");
            free(data_in);
            continue;
        }
        
        printf("done.\n");
        
        /* Perform correction */
        
        printf("Perform correction... "); fflush(stdout);
        
        // Input
        int naxis1 = FOURK, naxis2 = FOURK, naxis3i = nInSegs, naxis3o = nOutSegs;
        double crpix[2] = {ephem.crpix1, ephem.crpix2};
        double obs_v[3] = {ephem.obs_vr, ephem.obs_vn, ephem.obs_vw};
        double rsun_obs = ephem.rSun_obs, rsun_ref = ephem.rSun_ref;
        double crlt_obs = ephem.crlt_obs, crota2 = ephem.crota2;
        
        // Optional input, specify them for ease of interfacing
        double max_ang_in = 60., thresh_blos_in = 60., thresh_bmag_in = 250., pix_in = ephem.cdelt;
        int rad_los_sep_in = 2;
        
        // Output
        float doppcal_bias = 0.;     // output bias of Doppler vel. zero pt., subtract for true
        float *data_out = (float *) (malloc(FOURK2 * nOutSegs * sizeof(float)));
        
        // Fortran function
        
        doppcal_estimate_ (&naxis1, &naxis2, &naxis3i, &naxis3o,
                           data_in, data_out,
                           crpix, obs_v, &rsun_obs, &rsun_ref, &crlt_obs, &crota2,
                           &doppcal_bias,
                           &max_ang_in, &thresh_blos_in, &pix_in,
                           &thresh_bmag_in, &rad_los_sep_in,
                           &a0, &a2, &a4);
        
        printf("done.\n");
        printf("doppcal_bias = %f\n", doppcal_bias);
        
        /* Output */
        
        printf("Writing output... ");
        
        DRMS_Record_t *outRec = outRS->records[irec];
        
        status = writeOutput(inRec, outRec, writeImg, data_out, doppcal_bias, a0, a2, a4);
        
        free(data_in); free(data_out);      // Clean up regardless of status
        
        if (status) {
            SHOW("output error, record skipped.\n");
            continue;
        }
        
        printf("done.\n");
        
    }   // irec
    
	/* Clean up */
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    drms_close_records(inRS, DRMS_FREE_RECORD);
    
    return DRMS_SUCCESS;
    
}   // DoIt



// =============================================

/*
 * Read input segments into a data cube
 * Get ephemeris; if not FD then fail
 */

int getInput(DRMS_Record_t *inRec, double *data_in, struct ephemeris *ephem)
{
    
    int status = 0;
    
    // Ephemeris, no error checking
    
    getEphemeris(inRec, ephem);
    
    // Read in segments in a loop, nInSeg as global const
    
    for (int iSeg = 0; iSeg < nInSegs; iSeg++) {
        
        DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, inSegNames[iSeg]);
        if (!inSeg) {
            printf("segment %s error... ", inSegNames[iSeg]);
            return -1;
        }
        
        DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);    // everything in double
        if (status) {
            printf("array %s erro... ", inSegNames[iSeg]);
            return 1;
        }
        
        int nx = inArray->axis[0], ny = inArray->axis[1];
        int nxny = nx * ny;     // FOURK2
        if (nx != FOURK || ny != FOURK) {
            printf("array %s is not full disk... ", inSegNames[iSeg]);
            drms_free_array(inArray);
            return 2;
        }
        
        double *inData = (double *) (inArray->data);
        double *data_in_lead = data_in + iSeg * nxny;
        memcpy(data_in_lead, inData, nxny * sizeof(double));        // destination<-source copy over
        
        drms_free_array(inArray);
        
    }   // iSeg
    
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
    
    ephem->crota2 = drms_getkey_double(inRec, "CROTA2", &status);    // rotation
    double sina = sin(ephem->crota2 * RADSINDEG);
    double cosa = cos(ephem->crota2 * RADSINDEG);
    
    ephem->crlt_obs = drms_getkey_double(inRec, "CRLT_OBS", &status);
    ephem->crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status);
    
    double crvalx = drms_getkey_double(inRec, "CRVAL1", &status);
    double crvaly = drms_getkey_double(inRec, "CRVAl2", &status);
    double crpix1 = drms_getkey_double(inRec, "CRPIX1", &status);
    double crpix2 = drms_getkey_double(inRec, "CRPIX2", &status);
    double cdelt = drms_getkey_double(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
    ephem->crpix1 = PIX_X(0.0,0.0);        // Real CRPIX1, starting at 1
    ephem->crpix2 = PIX_Y(0.0,0.0);
    ephem->cdelt = cdelt;
    
    ephem->rSun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
    if (status) {
        ephem->rSun_ref = 6.96e8;
    }
    ephem->rSun_obs = drms_getkey_double(inRec, "RSUN_OBS", &status);
    if (status) {
        double dSun = drms_getkey_double(inRec, "DSUN_OBS", &status);
        ephem->rSun_obs = asin(ephem->rSun_ref / dSun) * RAD2ARCSEC;       // in arcsec
    }
    
    ephem->obs_vr = drms_getkey_double(inRec, "OBS_VR", &status);
    ephem->obs_vn = drms_getkey_double(inRec, "OBS_VN", &status);
    ephem->obs_vw = drms_getkey_double(inRec, "OBS_VW", &status);
    
    ephem->t_rec = drms_getkey_time(inRec, "T_REC", &status);
    
    return 0;
    
}


// =============================================

/*
 * Write output, loop over segments
 * Copy over keywords, set key doppcal_bias
 */

int writeOutput(DRMS_Record_t *inRec, DRMS_Record_t *outRec, int writeImg,
                float *data_out, float doppcal_bias, double a0, double a2, double a4)
{
    
    int status = 0;
    
    // Write segments in a loop, nOutSeg as global const
    
    if (writeImg) {
    
	    for (int iSeg = 0; iSeg < nOutSegs; iSeg++) {
    
 	       DRMS_Segment_t *outSeg = drms_segment_lookup(outRec, outSegNames[iSeg]);
	        if (!outSeg) {
				printf("segment %s error... ", outSegNames[iSeg]);
	            return -1;
	        }
        
	        int nx = FOURK, ny = FOURK, nxny = nx * ny;
	        int dims[2] = {nx, ny};
        
	        // Write as float, internal coversion to appropriate type
        
	        float *outData = (float *) (malloc(nxny * sizeof(float)));
			float *data_out_lead = data_out + iSeg * nxny;
	        memcpy(outData, data_out_lead, nxny * sizeof(float));        // destination<-source copy over
        
	        DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, outData, &status);
	        if (status) {
	            printf("array %s create error... ", outSegNames[iSeg]);
	            free(outData);
	            return 1;
	        }
        
	        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
			outArray->israw = 0;        // always compressed
	        outArray->bzero = outSeg->bzero;
	        outArray->bscale = outSeg->bscale;
        
	        status = drms_segment_write(outSeg, outArray, 0);
	        if (status) {
	            printf("array %s write error... ", outSegNames[iSeg]);
	            free(outData);
	            return 2;
	        }
        
	        drms_free_array(outArray);      // free outData along
    
	    }       // iSeg
	    
	}
    
    // Keywords
    
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
    
    drms_setkey_float(outRec, "DOPPBIAS", doppcal_bias);
    drms_setkey_string(outRec, "BUNIT_000", "cm/s");
    drms_setkey_string(outRec, "BUNIT_001", " ");
    drms_setkey_string(outRec, "BUNIT_002", " ");
    
    drms_setkey_double(outRec, "DIFFF_A0", a0);         // Jan 28 2019
    drms_setkey_double(outRec, "DIFFF_A2", a2);
    drms_setkey_double(outRec, "DIFFF_A4", a4);
    
    TIME val, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    // Link
    
    DRMS_Link_t *link = hcon_lookup_lower(&outRec->links, "BDATA");
    if (link) drms_link_set("BDATA", outRec, inRec);
    
    //
    
    return 0;
}

