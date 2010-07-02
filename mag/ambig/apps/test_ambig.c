/*
 * Module name:		test_ambig.c
 *
 * Description:		Testing wrapper for Graham Barnes' disambiguation module
 *
 * Original source:	Fortran wrapper by Graham and Igor (~graham/ME0/HMI/ambig.c)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * External function:
 *			ambig.f
 *
 * Version:		v1.0		Jan 25 2009
 *
 * Issues:
 *			v1.0
 *
 *
 * Example:
 * limit s u
 * test_ambig "in=su_keiji.vmagf_2d_720s_3_fltprf[2010.03.29_12:00:00_TAI]" "out=su_xudong.test_vmagf_a" "mask=su_xudong.test_armask_720s" "geometry=2"
 * test_ambig "in=su_xudong.test_patch_vec[2010.03.29_12:00:00_TAI][2]" "out=su_xudong.test_patch_vec_a" "geometry=1"
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define PI	(M_PI)
#define	DTOR	(PI / 180.)
#define Rad2arcsec	(3600. * 180. / PI) 


#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#include "timing.c"

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

char *module_name = "test_ambig";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_STRING, "mask", " ", "Bitmap series name, as AR mask, for full disk"},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=all messages"},
    {ARG_INT, "geometry", "1", "1 for patch; 2 for full disk."},
    {ARG_INT, "npad", "200", "Pixel number to pad with zeros for potential field."},
    {ARG_INT, "nap", "10", "Pixel number to apodize for potential field."},
    {ARG_INT, "ntx", "20", "Tile number in x (lon) direction for pf on a sphere."},
    {ARG_INT, "nty", "20", "Tile number in y (lat) direction for pf on a sphere."},
    {ARG_FLOAT, "bthresh", "3.e2", "Threshold field strength."},
    {ARG_INT, "seed", "1", "Input random number seed (seed>0)."},
    {ARG_INT, "neq", "10", "Num of reconfigurations attempted at each temperature setting."},
    {ARG_FLOAT, "lambda", "1.", "Weighting factor between div. and vert. current density."},
    {ARG_FLOAT, "tfac0", "2.", "Input factor to scale initial temperature (tfac0>0)."},
    {ARG_FLOAT, "tfactr", "0.990", "Input factor to reduce temperature (0<tfactr<1)."},
    {ARG_END}
};

/* ##### Prototypes for external Fortran functions ##### */
extern void ambig_(int *geometry,
       int *nx, int *ny,
       int *npad,
       float *xcen, float *ycen,
       int *verb,
       float *lambda,
       int *neq, float *tfactr, float *tfac0,
       int *seed,
       int *ntx, int *nty,
       float *Bx, float *By, float *Bz, float *dBt, float *probBa,
       int *bitmap,
       float *radius,
       int *nap, float *bthresh);

/* ################## Main Module ################## */

int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS, *outRS;
    int irec, nrecs;
    DRMS_Record_t *inRec, *outRec;

    int nseg = 10;
    char *segName[10] = {"field", "inclination", "azimuth", "alpha_mag", 
                         "field_err", "inclination_err", "alpha_err", 
                         "field_inclination_err", "field_alpha_err", "inclination_alpha_err"};
    
    DRMS_Segment_t *inSeg[nseg], *outSeg_flag, *outSeg_prob;
    DRMS_Array_t *inArray[nseg], *outArray_flag, *outArray_prob;
    float *inData[nseg];

    // Link
    DRMS_Link_t *magLink;

    // Bitmap mask, comes with patch, from mask data series for full disk
    // Read in as char *mask, converted to int *bitmap later
    // Query for mask data created from series name and time
    char *maskSeries, *maskQuery = NULL, *trec_str = NULL;
    DRMS_RecordSet_t *maskRS;
    DRMS_Record_t *maskRec;
    DRMS_Segment_t *maskSeg;
    DRMS_Array_t *maskArray;
    char *mask;
    int mask_id;	// id of main feature identified in patch

    int geometry;
    int npad, nap, ntx, nty, seed, neq;
    float bthresh, lambda, tfac0, tfactr;
    int verbflag;
    int outDims[2];

    int i, j, l, m;
    int nx, ny, nxny, nxnye, nxnyg, nerode, ngrow;

    TIME t_rec;

    float radius, xcen, ycen;  // Pointing information to be computed and passed to ambig
    float crpix1, crpix2, cdelt1, cdelt2;  // Keywords
    float im_scale, rsun_ref, dsun_obs;  // More keywords
    float crvalx, crvaly, cdelt, crota2, sina, cosa;
    // for patch
    float center_x, center_y;
    int hw_x, hw_y;
    int ll[2], ur[2];

    float *Bx, *By, *Bz, *dBt;	// Arrays to be computed + passed to ambig
    float *probBa;  // Array to be computed in ambig and returned

    // Inputs to be read from DRMS:  Use arrays where they need to be
    // saved and passed to ambig, otherwise just use as floats temporarily,
    // to compute Bx,By,Bz. - ELW-20100226  (was B_los; B_transvers; Azimuth)
    float Bazm, Bmag, Bfil, Binc;  	// Instances used to compute Bxyz(i,j)
    float BmagBfil, SinBinc, BmagBfilSinBinc, CosBinc;
    float varBfil, varBinc, varBmag, varBfilBmag, varBfilBinc, varBmagBinc;
    float xx, yy, r2, rad2;
    int *bitmap = NULL;
    int *erodemap, *growmap;
    char *ambig_flag;

    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);


/* cmdparams defined in jsoc_main.h:  extern CmdParams_t cmdparams; */

    /* Get parameters */
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    maskSeries = (char *)params_get_str(&cmdparams, "mask");
    verbflag = params_get_int(&cmdparams, "VERB");
    geometry = params_get_int(&cmdparams, "geometry");
    npad = params_get_int(&cmdparams, "npad");
    nap = params_get_int(&cmdparams, "nap");
    ntx = params_get_int(&cmdparams, "ntx");
    nty = params_get_int(&cmdparams, "nty");
    seed = params_get_int(&cmdparams, "seed");
    neq = params_get_int(&cmdparams, "neq");
    bthresh = params_get_float(&cmdparams, "bthresh");
    lambda = params_get_float(&cmdparams, "lambda");
    tfac0 = params_get_float(&cmdparams, "tfac0");
    tfactr = params_get_float(&cmdparams, "tfactr");

    /* Check parameters */
    if (nap > npad) {nap = npad; SHOW("nap set to npad\n");}
    if (seed <= 0) {seed = 1; SHOW("seed set to 1\n");}
    if (lambda < 0) {lambda = 1.0; SHOW("lambda set to 1.0\n");}
    if (tfac0 < 0) {tfac0 = 2.; SHOW("tfac0 set to 2.0\n");}
    if (neq <= 0) {neq = 20; SHOW("neq set to 20\n");}
    if ((tfactr < 0) || (tfactr > 1)) {tfactr = 0.995; SHOW("tfactr set to 0.995\n");}
    if (geometry < 1 || geometry > 2) {DIE("Unknown geometry flag");}

    /* Open input */

    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;

    /* Create output */
    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output recordset not created");

    /* Do this for each record */
    for (irec = 0; irec < nrecs; irec++) {
    
        /* Measure time */
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }
        
        /* Input record and data */
        inRec = inRS->records[irec];
        
        /* Processing keywords */

        t_rec = drms_getkey_time(inRec, "T_REC", &status);

        crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation
        sina = sin(crota2 * DTOR); 
        cosa = cos(crota2 * DTOR);
        cdelt = drms_getkey_float(inRec, "CDELT1", &status);	// in arcsec, assuming dx=dy

        rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
        if (status) rsun_ref = 6.96e8;
        dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status);
   	    radius = asin(rsun_ref / dsun_obs) * Rad2arcsec / cdelt;

        if (geometry == 1) {
            crvalx = drms_getkey_float(inRec, "IMCRVAL1", &status);	// center of solar disc in arcsec
            crvaly = drms_getkey_float(inRec, "IMCRVAL2", &status);
            crpix1 = drms_getkey_float(inRec, "IMCRPIX1", &status);	// disk center in ccd, original
            crpix2 = drms_getkey_float(inRec, "IMCRPIX2", &status);

            center_x = drms_getkey_float(inRec, "CRPIX1", &status) - 1; // Center
            center_y = drms_getkey_float(inRec, "CRPIX2", &status) - 1;
            hw_x = drms_getkey_int(inRec, "HWIDTH1", &status);      // Half width
            hw_y = drms_getkey_int(inRec, "HWIDTH2", &status);

            // center of patch relative to disk center in arcsec
            xcen = center_x - (PIX_X(0.0,0.0) - 1.);
            ycen = center_y - (PIX_Y(0.0,0.0) - 1.);
printf("xcen=%f, ycen=%f\n", xcen, ycen); fflush(stdout);

            // Coordinates of cut out rectangle
            ll[0] = (int)center_x - hw_x; ll[1] = (int)center_y - hw_y;
            ur[0] = (int)center_x + hw_x; ur[1] = (int)center_y + hw_y;         
        } else {
            crvalx = drms_getkey_float(inRec, "CRVAL1", &status);	// center of solar disc in arcsec
            crvaly = drms_getkey_float(inRec, "CRVAL2", &status);
            crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);	// disk center in ccd, original
            crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);

            xcen = PIX_X(0.0,0.0);
            ycen = PIX_Y(0.0,0.0);
        }

        /* Find segments and read in data */
        for (i = 0; i < nseg; i++) {
            inSeg[i] = drms_segment_lookup(inRec, segName[i]);
            
            if (geometry == 1) {
                inArray[i] = drms_segment_readslice(inSeg[i], DRMS_TYPE_FLOAT, ll, ur, &status);
            } else {
                inArray[i] = drms_segment_read(inSeg[i], DRMS_TYPE_FLOAT, &status);
            }
            
            if (status) {
                for (j = 0; j <= i; j++)
  	            drms_free_array(inArray[i]);
                DIE("Segment reading error \n");
            }
            inData[i] = (float *)inArray[i]->data;
        }
	
        /* Dimensions */
        nx = inArray[0]->axis[0]; ny = inArray[0]->axis[1];
        if (geometry == 1) {
            if (nx != 2 * hw_x + 1 || ny != 2 * hw_y + 1) DIE("wrong cut out dimension");
        }
        nxny = nx * ny;
        printf("nx=%d, ny=%d\n", nx, ny);
        
        /* Get mask, if any */
        if (geometry == 1) {
            maskSeg = drms_segment_lookup(inRec, "bitmap");		// from segment, for patch
            mask_id = drms_getkey_int(inRec, "MASK", &status);		// main AR in patch
        } else {
            // Create query string for mask, from series name and time
            // e.g. "su_turmon.armask_45s" and "[2010.03.29_12:00:00_TAI]"
            maskQuery = (char *)malloc(100 * sizeof(char));
            trec_str = (char *)malloc(30 * sizeof(char));
            sprint_time(trec_str, t_rec, "TAI", 0);
            sprintf(maskQuery, "%s[%s]", maskSeries, trec_str);
            // Get mask record, continue if not found
            maskRS = drms_open_records(drms_env, maskQuery, &status);
            if (status || maskRS->n != 1) {
                maskSeg = NULL;				// no mask used;
                printf("no mask used\n");
            } else {
                maskRec = maskRS->records[0];
                maskSeg = drms_segment_lookup(maskRec, "bitmap");
                mask_id = 1;
                printf("%s, %s\n", trec_str, maskSeries);
                printf("#%d: t_inrec = %lf; t_mask = %lf\n", 
                       irec, t_rec, drms_getkey_time(maskRec, "T_REC", &status));
            }
        }

        // Allocate arrays to pass to ambig()
        // Re-created for each record as patch could have different dimensions
        Bx = (float *) malloc(nxny * sizeof(float));
        By = (float *) malloc(nxny * sizeof(float));
        Bz = (float *) malloc(nxny * sizeof(float));
        dBt = (float *) malloc(nxny * sizeof(float));
    
        //
        // Parse inData, and compute Bx,By,Bz along the way
        //

        for (i = 0; i < nxny; i++) {

            // Field
            Bmag = inData[0][i];		// "field"
            Binc = inData[1][i];		// "inclination"
            Bazm = inData[2][i];		// "azimuth"
            Bfil = inData[3][i];		// "alpha_mag"

            // Variances
            varBmag = inData[4][i];	    	// "field_err"
            varBinc = inData[5][i];		// "inclination_err"
            //varBfil = inData[6][i];		// "alpha_err"
// XXXX zero this until it's clear the fill factor is being computed
            varBfil = 0.;

            // Covariances
            varBmagBinc = inData[7][i];		// "field_inclination_err"
            varBfilBmag = inData[8][i];		// "field_alpha_err"
            varBfilBinc = inData[9][i];		// "inclination_alpha_err"

            // Compute and set arrays for ambig.  First check for NANs in input
            // data, and set array elements accordingly.
            if (isnan(Bmag) || isnan(Binc) || isnan(Bazm) || isnan(Bfil)) {

                dBt[i] = 1.E9;
                Bx[i] = 0.; By[i] = 0.; Bz[i] = 0.;

            } else {

                Binc *= DTOR, Bazm *= DTOR; 		// Convert to radians

                BmagBfil = Bmag * Bfil;				// Vars for computation
                CosBinc = cos(Binc); SinBinc=sin(Binc);
                BmagBfilSinBinc = BmagBfil * SinBinc;

		// Arrays for ambig
        
                // zero azi pointing West, zero incl pointing into sun
                /*
                Bx[i] = BmagBfilSinBinc * cos(Bazm);
                By[i] = BmagBfilSinBinc * sin(Bazm);
                Bz[i] = - BmagBfil * cos(Binc);
                */

                // zero azi pointing North, zero incl pointing out from sun
                Bx[i] =-BmagBfilSinBinc * sin(Bazm);
                By[i] = BmagBfilSinBinc * cos(Bazm);
                Bz[i] = BmagBfil * cos(Binc);
                               
                // zero azi pointing North, zero incl pointing into sun
                /*
                Bx[i] =-BmagBfilSinBinc * sin(Bazm);
                By[i] = BmagBfilSinBinc * cos(Bazm);
                Bz[i] =-BmagBfil * cos(Binc);
                */
                               
                // If Bmag=0 or Bfil=0, set dBt to something large, else compute:
                dBt[i] = (BmagBfil==0) ? 1.E9 : 
                                         (BmagBfil * 
                                          sqrt (SinBinc * SinBinc * (varBfil / (Bfil * Bfil) + varBmag / (Bmag * Bmag)
      	                                        + 2.*varBfilBmag/BmagBfil)
                                                + varBinc * CosBinc * CosBinc
                                                + (varBfilBinc / Bfil + varBmagBinc / Bmag) * SinBinc * CosBinc));

            }
        }
        printf("Arrays read. \n");
 
        /* Free input arrays */
        for (i = 0; i < nseg; i++) {
           drms_free_array(inArray[i]);
        }

        if (geometry == 1 && maskSeg == NULL) DIE("Error for patch mask");	// mandatory for patch
        if (maskSeg != NULL) {
            maskArray = drms_segment_read(maskSeg, DRMS_TYPE_CHAR, &status);
            if (status) DIE("Error reading mask");
            mask = (char *)maskArray->data;
            if (maskArray->axis[0] != nx || maskArray->axis[1] != ny)
                DIE("Mask has wrong dimension");
            // Create bitmap
            bitmap = (int *)calloc(nxny, sizeof(int));
            for (i = 0; i < nxny; i++) {
               bitmap[i] = (mask[i] == mask_id) ? 1 : 0;
            }
            drms_free_array(maskArray);
        } else {
            // Construct bitmap based on transverse field strength and distance from disk center.
            bitmap = (int *)calloc(nxny, sizeof(int));
            rad2 = radius * radius;
            for (i = 0; i < nxny; i++) {
               j = i / nx;
               xx = i - j * nx - xcen + 1;
               yy = j - ycen + 1;
               r2 = xx * xx + yy * yy;
               if (r2 < rad2) {
                  bitmap[i] = (sqrt(Bx[i] * Bx[i] + By[i] * By[i]) > bthresh * (2. - sqrt(1. - r2 / rad2))) ? 1 : 0;
                  //bitmap[i] = (sqrt(Bx[i] * Bx[i] + By[i] * By[i]) > bthresh) ? 1 : 0;
               } else {
                  bitmap[i] = 0;
               }
            }
            // Erode this bitmap to remove isolated above threshold pixels. 
            nerode = 1;
            nxnye = (nx + 2 * nerode) * (ny + 2 * nerode);
            erodemap = (int *)calloc(nxnye, sizeof(int));
            for (i = 0; i < nxnye; i++) {
               erodemap[i] = 1;
            }
            for (i = 0; i < nxny; i++) {
               if (bitmap[i] == 0) {
                  j = i + nerode * (nx + 1 + 2 * (i / nx + nerode));

                  erodemap[j - 2 * nerode - nx - 1] = 0;
                  erodemap[j - 2 * nerode - nx] = 0;
                  erodemap[j - 2 * nerode - nx + 1] = 0;
                  erodemap[j - 1] = 0;
                  erodemap[j] = 0;
                  erodemap[j + 1] = 0;
                  erodemap[j + 2 * nerode + nx - 1] = 0;
                  erodemap[j + 2 * nerode + nx] = 0;
                  erodemap[j + 2 * nerode + nx + 1] = 0;
               }
            }
            // Now grow the bitmap to provide some padding around above threshold pixels. 
            ngrow = 5;
            nxnyg = (nx + 2 * ngrow) * (ny + 2 * ngrow);
            growmap = (int *)calloc(nxnyg, sizeof(int));
            for (i = 0; i < nxnyg; i++) {
               growmap[i] = 0;
            }
            for (i = 0; i < nxnye; i++) {
               if (erodemap[i] == 1) {
                  j = i + ngrow * 2 * (i / (nx + 2 * nerode));
                  for (l = 0; l <= 2 * ngrow; l++) {
                     for (m = 0; m <= 2 * ngrow; m++) {
                        growmap[j + l + m * (nx + 2 * (nerode + ngrow))] = 1;
                     }
                  }
               }
            }
            // Extract the final bitmap from the padded/eroded/grown bitmap.
            for (i = 0; i < nxny; i++) {
               bitmap[i] = growmap[i + (nerode + ngrow) * (nx + 1 + 2 * (i / nx + nerode + ngrow))];
            }
        }

        /* Output data */
        outDims[0] = nx; outDims[1] = ny;	// Output array dimensions
        probBa = (float *) malloc(nxny * sizeof(float));
        ambig_flag = (char *) calloc(nxny, sizeof(char));

        /* This is the working part */
        ambig_(&geometry,
               &nx, &ny,
               &npad,
               &xcen, &ycen,
               &verbflag,
               &lambda,
               &neq, &tfactr, &tfac0,
               &seed,
               &ntx, &nty,
               Bx, By, Bz, probBa,
               dBt, bitmap,
               &radius,
               &nap, &bthresh);
        /* Flag pixels for which the azimuth angle needs to be changed */
        // Updated by xudong jun 30: Bx[i] > 0 instead of By[i] < 0
        for (i = 0; i < nxny; i++) {
            ambig_flag[i] = (Bx[i] > 0.) ? 1 : 0;
        }

        /* Stats */
        int nancount = 0, outsz = outDims[0] * outDims[1]; 
        for (i = 0; i < outsz; i++) {
            if (isnan(probBa[i])) nancount++;
        }

 
        /* Output segment */
        outArray_flag = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, ambig_flag, &status);
        if (status) DIE("Error creating flag array");
        outArray_prob = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, probBa, &status);
        if (status) DIE("Error creating prob array");

        /* Output record */
        outRec = outRS->records[irec];
        outSeg_flag = drms_segment_lookup(outRec, "ambig_flag");
        outSeg_prob = drms_segment_lookup(outRec, "confidence");
        for (i = 0; i < 2; i++) {
          outSeg_flag->axis[i] = outArray_flag->axis[i];	// For variable dimensions
          outSeg_prob->axis[i] = outArray_prob->axis[i];
        }
        outArray_flag->parent_segment = outSeg_flag;
        outArray_prob->parent_segment = outSeg_prob;

        /* Result writing */
        status = drms_segment_write(outSeg_flag, outArray_flag, 0);
        if (status) DIE("Problem writing flag file");
        drms_free_array(outArray_flag);
        status = drms_segment_write(outSeg_prob, outArray_prob, 0);
        if (status) DIE("Problem writing prob file");
        drms_free_array(outArray_prob);

        /* Set keywords */
        // Prime key
    	drms_copykey(outRec, inRec, "T_REC");
    	if (geometry == 1) {
    	    drms_copykey(outRec, inRec, "PNUM");
    	}
        // Info
    	drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        // Parameters
// XXXX update these to reflect geometry flag?
        drms_setkey_int(outRec, "npad", npad);
        drms_setkey_int(outRec, "nap", nap);
        drms_setkey_int(outRec, "ntx", ntx);
        drms_setkey_int(outRec, "nty", nty);
        drms_setkey_int(outRec, "seed", seed);
        drms_setkey_int(outRec, "neq", neq);
        drms_setkey_float(outRec, "bthresh", bthresh);
        drms_setkey_float(outRec, "lambda", lambda);
        drms_setkey_float(outRec, "tfac0", tfac0);
        drms_setkey_float(outRec, "tfactr", tfactr);
        // stats
        drms_setkey_int(outRec, "TOTVALS", outsz);
        drms_setkey_int(outRec, "DATAVALS", outsz - nancount);
        drms_setkey_int(outRec, "MISSVALS", nancount);
    
        /* Set link */
        if (geometry == 1) {
            magLink = hcon_lookup_lower(&outRec->links, "PATCH");
            if (magLink) {
                drms_setlink_static(outRec, "PATCH", inRec->recnum);
            }
        } else {
            magLink = hcon_lookup_lower(&outRec->links, "MDATA");
            if (magLink) {
                drms_setlink_static(outRec, "MDATA", inRec->recnum);
            }
        }

        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                   irec, wt - wt1, ct - ct1);
        }

        /* Clean up */
        free(Bx); free(By); free(Bz); free(dBt);  // ELW20100226 - mainly for sanity
        if (bitmap != NULL) { free(bitmap); bitmap = NULL; }
        if (maskQuery != NULL) { free(maskQuery); maskQuery = NULL; }
        if (trec_str != NULL) { free(trec_str); trec_str = NULL; }

    }

    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);

    if (verbflag) {
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
                wt - wt0, ct - ct0);
    }


    
    return(DRMS_SUCCESS);
}
