/*
 * Module name:		disambig.c
 *
 * Description:		Wrapper for Graham Barnes' disambiguation module
 *
 * Original source:	Fortran wrapper by Graham and Igor (~graham/ME0/HMI/ambig.c)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * External function:
 *			ambig.f
 *
 * Version:
 *			v1.0		Jan 25 2009
 *			v1.1		Mar 21 2012
 *			v1.2		Mar 27 2012
 *			v1.3		Apr 09 2012
 *			v2.0		Mon 04 2012
 *
 * Issues:
 *			v1.0
 *			v1.1
 *			Added 3 solutions for weak field
 *			Re-defined confidence array (0-100, 100 is good)
 *			v1.2
 *			Added preliminary full disk noise mask. the I/O needs to be updated once
 *			the format of the mask is determined
 *			v1.3
 *			Added query for speed
 *			v2.0
 *			Record all three weak field solutions
 *			Implemented noise mask, added -l flag for control
 *			Clean-up
 *
 *
 * Example:
 * limit s u
 *
 *
 */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "copy_me_keys.c"
// Added Mar 21 for geometry
#include "cartography.c"
#include "noisemask.c"

#define PI	(M_PI)
#define	DTOR	(PI / 180.)
#define Rad2arcsec	(3600. * 180. / PI)
#define FOURK (4096)
#define FOURK2	(16777216)

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

// Quality bit for disambiguation
// Now occupy lower 8 bits
// lowest bit disambiguated or not
// 2nd last bit annealed or not
// 3rd last bit strong or not
#define ambCode0 0x0000 // not disambiguated, turned
#define ambCode1 0x0001 // weak not annealed (pf/radial acute/random)
#define ambCode2 0x0003 // weak annealed
#define ambCode3 0x0007 // strong annealed

// switches
#define QUALMAP 1		// QUALMAP not updated if 0

int getAmbCode (char prob)
{
	if (prob < 11)
		return ambCode0;
	else if (prob < 51)
		return ambCode1;
	else if (prob < 61)
		return ambCode2;
	else
		return ambCode3;
}

// Wrapper for Yang's mask function

int createMask(float obs_vr, float radius, float x0, float y0, double *noiseMask)
{return noisemask(4096, 4096, (int)x0, (int)y0, (int)radius, obs_vr, noiseMask);}

//====================

char *module_name = "disambig_new";	/* Module name */
char *version_id = "2012 Jun 04";  /* Version number */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=minimal messages; 2=all messages"},
		{ARG_FLAG, "l", "", "Flag to use linear noise threshold (AMBTHRn)"},
    {ARG_INT, "ksig", "1", "Multiplier of the noise mask"},
    {ARG_INT, "offset", "50", "Constant added to the noise mask"},
    {ARG_INT, "AMBGMTRY", "0", "0 for automatic selection; 1 for planar; 2 for spherical."},
//	{ARG_INT, "AMBWEAK", "1", "0 for random; 1 for potential field; 2 for most radial."},
    {ARG_INT, "AMBNEROD", "1", "Number of pixels by which to erode map of above threshold pixels."},
    {ARG_INT, "AMBNGROW", "5", "Number of pixels by which to grow eroded map."},
    {ARG_INT, "AMBNPAD", "20", "Pixel number to pad with zeros for potential field."},
    {ARG_INT, "AMBNAP", "10", "Pixel number to apodize for potential field."},
    {ARG_INT, "AMBNTX", "20", "Tile number in x (lon) direction for pf on a sphere."},
    {ARG_INT, "AMBNTY", "20", "Tile number in y (lat) direction for pf on a sphere."},
    {ARG_FLOAT, "AMBBTHR0", "2.0e2", "Threshold field strength."},
    {ARG_FLOAT, "AMBBTHR1", "3.0e2", "Threshold field strength."},
    {ARG_INT, "AMBSEED", "1", "Input random number seed (seed>0)."},
    {ARG_INT, "AMBNEQ", "10", "Num of reconfigurations attempted at each temperature setting."},
    {ARG_FLOAT, "AMBLMBDA", "1.", "Weighting factor between div. and vert. current density."},
    {ARG_FLOAT, "AMBTFCT0", "2.", "Input factor to scale initial temperature (tfac0>0)."},
    {ARG_FLOAT, "AMBTFCTR", "0.990", "Input factor to reduce temperature (0<tfactr<1)."},
    {ARG_END}
};


/* ##### Prototypes for external Fortran functions ##### */
extern void ambig_(int *geometry,
				   int *weak,
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
				   int *nap, float *bthresh0, float *bthresh1);

/* ################## Main Module ################## */

int DoIt(void)
{

	char ambcodev[50];
	sprintf(ambcodev,"%s %s", module_name, version_id);

	int status = DRMS_SUCCESS;
	char *inQuery, *outQuery;
	DRMS_RecordSet_t *inRS, *outRS;
 	int irec, nrecs;
	DRMS_Record_t *inRec, *outRec;
	
	int nseg = 10;
    char *segName[10] = {"field", "inclination", "azimuth", "alpha_mag", 
		"field_err", "inclination_err", "alpha_err", 
		"field_inclination_err", "field_alpha_err", "inclination_alpha_err"};
    
    DRMS_Segment_t *inSeg[nseg], *inSeg_conf, *inSeg_qual, *outSeg_flag, *outSeg_prob, *outSeg_conf, *outSeg_qual;
    DRMS_Array_t *inArray[nseg], *inArray_conf, *inArray_qual, *outArray_flag, *outArray_prob, *outArray_conf, *outArray_qual;
    float *inData[nseg];
	
    // Link
    DRMS_Link_t *magLink;
	
	int useMask, useMask0;		// user specified & for each record
	double *noiseMask = NULL;		// full disk mask; real size mask
	
	DRMS_Segment_t *maskSeg;		// this is actually bitmap
    DRMS_Array_t *maskArray;
    char *mask;
	char harpflag;
    int mask_id, harpnum;	// id of main feature identified in patch
	
    int geometry, weak;
    int npad, nap, ntx, nty, seed, neq;
    float bthresh0, bthresh1, lambda, tfac0, tfactr;
    int verbflag;
    int outDims[2];
    int ksig, offset;
	
    int i, j, l, m, i0, j0;
    int nx, ny, nxny, nx0, ny0, nxny0, nxnye, nxnyg, nerode, ngrow;
	
    TIME t_rec;
	float obs_vr;		// radial speed
//	int vr_min, vr_max;
	float x0, y0;		// disk center, plate scale
	
    float radius, xcen, ycen;  // Pointing information to be computed and passed to ambig
    float crpix1, crpix2, cdelt1, cdelt2;  // Keywords
    float im_scale, rsun_ref, dsun_obs;  // More keywords
    float crvalx, crvaly, cdelt, crota2, sina, cosa;
    // for patch
    int ll0[2], ur0[2];
    int ll[2], ur[2];
    float minlon, minlat, maxlon, maxlat;
	
    float *Bx, *By, *Bz, *dBt;	// Arrays to be computed + passed to ambig
    float *probBa;  // Array to be computed in ambig and returned
	
    // Inputs to be read from DRMS:  Use arrays where they need to be
    // saved and passed to ambig, otherwise just use as floats temporarily,
    // to compute Bx,By,Bz. - ELW-20100226  (was B_los; B_transvers; Azimuth)
    float Bazm, Bmag, Bfil, Binc;  	// Instances used to compute Bxyz(i,j)
    float BmagBfil, SinBinc, BmagBfilSinBinc, CosBinc;
    float varBfil, varBinc, varBmag, varBfilBmag, varBfilBinc, varBmagBinc;
    float xx, yy, r2, rad2;
    char *confidence;			// changed Oct 13 2011
    int *bitmap = NULL;
    int *erodemap, *growmap;
    char *ambig_flag;
    int *qual_map, *inData_qual, qual;						// Jun 22 Xudong
    char *confid_map, *inData_conf, *inData_conf_fd;
	
    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);

	/* Get parameters */
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    verbflag = params_get_int(&cmdparams, "VERB");
    geometry = params_get_int(&cmdparams, "AMBGMTRY");
	useMask = !(params_isflagset(&cmdparams, "l"));
	ksig = params_get_int(&cmdparams, "ksig");
    offset = params_get_int(&cmdparams, "offset");
	
//	weak = params_get_int(&cmdparams, "AMBWEAK");
	weak = 1;		// Mar 21 2012 set to be constant as now we compute all thress options
    nerode = params_get_int(&cmdparams, "AMBNEROD");
    ngrow = params_get_int(&cmdparams, "AMBNGROW");
    npad = params_get_int(&cmdparams, "AMBNPAD");
    nap = params_get_int(&cmdparams, "AMBNAP");
    ntx = params_get_int(&cmdparams, "AMBNTX");
    nty = params_get_int(&cmdparams, "AMBNTY");
    seed = params_get_int(&cmdparams, "AMBSEED");
    neq = params_get_int(&cmdparams, "AMBNEQ");
    bthresh0 = params_get_float(&cmdparams, "AMBBTHR0");
    bthresh1 = params_get_float(&cmdparams, "AMBBTHR1");
    lambda = params_get_float(&cmdparams, "AMBLMBDA");
    tfac0 = params_get_float(&cmdparams, "AMBTFCT0");
    tfactr = params_get_float(&cmdparams, "AMBTFCTR");
	
    /* Check parameters */
    if (nap > npad) {nap = npad; SHOW("nap set to npad\n");}
    if (seed <= 0) {seed = 1; SHOW("seed set to 1\n");}
    if (lambda < 0) {lambda = 1.0; SHOW("lambda set to 1.0\n");}
    if (tfac0 < 0) {tfac0 = 2.; SHOW("tfac0 set to 2.0\n");}
    if (neq <= 0) {neq = 20; SHOW("neq set to 20\n");}
    if ((tfactr < 0) || (tfactr > 1)) {tfactr = 0.995; SHOW("tfactr set to 0.995\n");}
    if (geometry < 0 || geometry > 2) {DIE("Unknown geometry flag");}
    if (weak < 0 || weak > 2) {DIE("Unknown weak flag");}
    if (nerode <= 0) {nerode = 1; SHOW("nerode set to 1\n");}
    if (ngrow <= 0) {ngrow = 5; SHOW("ngrow set to 5\n");}
	
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

		/* Set harp flag based on whether a HARPNUM is found */
        harpnum = drms_getkey_int(inRec, "HARPNUM", &status);	// HARP number
		// Changed Jun 28 by Xudong
		if (harpnum == DRMS_MISSING_INT) {
			harpflag = 0;
		}
		else {
			harpflag = 1;
		}
		status = 0;
		printf("harpflag=%d\n",harpflag);
		
		/* If specified, obtain full-disk noise mask */
		
		useMask0 = useMask && (!harpflag);		// no noise mask for harp now
		if (useMask0) {
			obs_vr = drms_getkey_float(inRec, "OBS_VR", &status);
			printf("status=%d, obs_vr=%f\n",status, obs_vr);
			if (harpflag) x0 = drms_getkey_double(inRec, "IMCRPIX1", &status)-1.;
				else x0 = drms_getkey_double(inRec, "CRPIX1", &status)-1.;
			if (harpflag) y0 = drms_getkey_double(inRec, "IMCRPIX2", &status)-1.;
				else y0 = drms_getkey_double(inRec, "CRPIX2", &status)-1.;
			noiseMask = (double *) (malloc(FOURK2 * sizeof(double)));
			/* Create the full disk mask using Yang's code, TBD */
			/****************************************************/
			if (createMask(obs_vr, radius, x0, y0, noiseMask)) {		// error
				if (noiseMask) free(noiseMask); noiseMask = NULL;
				useMask0 = 0; continue;
			}
			/****************************************************/
		}
		
printf("here\n");

        /* Depending on whether it's a patch or not, read in the
		 appropriate keywords and data segments. */
		
        if (harpflag) {
            // If unspecified, set geometry based on patch size
            if (geometry == 0) {
                minlon = drms_getkey_float(inRec, "LON_MIN", &status);
                minlat = drms_getkey_float(inRec, "LAT_MIN", &status);
                maxlon = drms_getkey_float(inRec, "LON_MAX", &status);
                maxlat = drms_getkey_float(inRec, "LAT_MAX", &status);
                if (maxlon - minlon < 20. && maxlat - minlat < 20.) {
                    geometry = 1;
                    printf("using planar geometry\n");
                } else {
                    geometry = 2;
                    printf("using spherical geometry\n");
                }
            }
            // Read the mask array; no longer used
            maskSeg = drms_segment_lookup(inRec, "bitmap");		// from segment, for patch
            nx0 = maskSeg->axis[0]; ny0 = maskSeg->axis[1];		// Get array dimensions from mask array.
            nxny0 = nx0 * ny0;
            // Pointing information.
            crvalx = drms_getkey_float(inRec, "IMCRVAL1", &status);	// center of solar disc in arcsec
            crvaly = drms_getkey_float(inRec, "IMCRVAL2", &status);
            crpix1 = drms_getkey_float(inRec, "IMCRPIX1", &status);	// disk center in ccd, original
            crpix2 = drms_getkey_float(inRec, "IMCRPIX2", &status);
            // Coordinates of cut out rectangle
            ll0[0] = drms_getkey_float(inRec, "CRPIX1", &status) - 1;	// lower left corner
            ll0[1] = drms_getkey_float(inRec, "CRPIX2", &status) - 1;
            ur0[0] = ll0[0] + nx0 - 1;					// upper right corner
            ur0[1] = ll0[1] + ny0 - 1;
            // Pad the actual patch
            nx = nx0 + 2*npad;
            ny = ny0 + 2*npad;
            ll[0] = ll0[0] - npad;
            ll[1] = ll0[1] - npad;
            ur[0] = ur0[0] + npad;
            ur[1] = ur0[1] + npad;
            // Ensure that this stays within the image
            if (ll[0] < 0) {
                nx = nx + ll[0];
                ll[0] = 0;
            }
            if (ll[1] < 0) {
                ny = ny + ll[1];
                ll[1] = 0;
            }
            if (ur[0] > 4095) {
                nx = nx + ur[0] - 4095;
                ur[0] = 4095;
            }
            if (ur[1] > 4095) {
                ny = ny + ur[1] - 4095;
                ur[1] = 4095;
            }
            nxny = nx * ny;
//            printf("nx0=%d, ny0=%d\n", nx0, ny0);
            // Different pointing information depending on geometry
            // Added by G. Barnes Aug 9 2011
            if (geometry == 1) {
				// center of patch relative to disk center
				xcen = (float)ll0[0] + 0.5*((float)nx0 - 1.) - (PIX_X(0.0,0.0) - 1.);
				ycen = (float)ll0[1] + 0.5*((float)ny0 - 1.) - (PIX_Y(0.0,0.0) - 1.);
            }
            if (geometry == 2) {
				// disk center relative to lower left corner of padded patch
				xcen = PIX_X(0.0,0.0) - (float)ll[0];
				ycen = PIX_Y(0.0,0.0) - (float)ll[1];
            }
            // Find segments and read in data
            for (i = 0; i < nseg; i++) {
                inSeg[i] = drms_segment_lookup(inRec, segName[i]);
                inArray[i] = drms_segment_readslice(inSeg[i], DRMS_TYPE_FLOAT, ll, ur, &status);
                if (status) {
                    for (j = 0; j <= i; j++)
						drms_free_array(inArray[i]);
                    DIE("Segment reading error \n");
                }
                inData[i] = (float *)inArray[i]->data;
            }
            // Jun 22 Xudong
#if QUALMAP==1            
            inSeg_qual = drms_segment_lookup(inRec, "qual_map");
            inArray_qual = drms_segment_readslice(inSeg_qual, DRMS_TYPE_INT, ll0, ur0, &status);
            if (status) DIE("Segment reading error \n");
            inData_qual = (int *) inArray_qual->data;
            
            inSeg_conf = drms_segment_lookup(inRec, "confid_map");
			//            inArray_conf = drms_segment_readslice(inSeg_conf, DRMS_TYPE_CHAR, ll0, ur0, &status);
            inArray_conf = drms_segment_read(inSeg_conf, DRMS_TYPE_CHAR, &status);
            if (status) DIE("Segment reading error \n");
            inData_conf_fd = (char *) inArray_conf->data;
            inData_conf = (char *) calloc(nxny0, sizeof(char));
            for (i = 0; i < nx0; i++) {
                for (j = 0; j < ny0; j++) {
                    inData_conf[i + j * nx0] = inData_conf_fd[ll0[0] + i + (ll0[1] + j) * 4096];
                }
            }
#endif
        }		// harp case done
		else {
            /* Full disk case */
            if (geometry != 2) geometry = 2;	// force spherical geometry for full disk
            // added by X. Sun, to be consistent with harp case
			i0 = 0; j0 = 0;
			ll0[0] = ll0[1] = ll[0] = ll[1] = 0;
			ur0[0] = ur0[1] = ur[0] = ur[1] = 0;
            crvalx = drms_getkey_float(inRec, "CRVAL1", &status);	// center of solar disc in arcsec
            crvaly = drms_getkey_float(inRec, "CRVAL2", &status);
            crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);	// disk center in ccd, original
            crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
            xcen = PIX_X(0.0,0.0);
            ycen = PIX_Y(0.0,0.0);
            // Find segments and read in data
            for (i = 0; i < nseg; i++) {
                inSeg[i] = drms_segment_lookup(inRec, segName[i]);
                inArray[i] = drms_segment_read(inSeg[i], DRMS_TYPE_FLOAT, &status);
                if (status) {
                    for (j = 0; j <= i; j++)
						drms_free_array(inArray[i]);
                    DIE("Segment reading error \n");
                }
                inData[i] = (float *) inArray[i]->data;
            }
            // Jun 22 Xudong
#if QUALMAP==1            
            inSeg_qual = drms_segment_lookup(inRec, "qual_map");
            inArray_qual = drms_segment_read(inSeg_qual, DRMS_TYPE_INT, &status);
            if (status) DIE("Segment reading error \n");
            inData_qual = (int *) inArray_qual->data;
            inSeg_conf = drms_segment_lookup(inRec, "confid_map");
            inArray_conf = drms_segment_read(inSeg_conf, DRMS_TYPE_CHAR, &status);
            if (status) DIE("Segment reading error \n");
            inData_conf = (char *) inArray_conf->data;
#endif	
            // Dimensions
            nx = inArray[0]->axis[0]; ny = inArray[0]->axis[1];
            nxny = nx * ny;
            nx0 = nx; ny0 = ny; nxny0 = nxny;
        }		// full disk case done
//        printf("nx=%d, ny=%d\n", nx, ny);
//        printf("xcen=%f, ycen=%f\n", xcen, ycen); fflush(stdout);

		printf("ll=[%d,%d], ur=[%d,%d], nx=%d, ny=%d\n", ll[0], ll[1], ur[0], ur[1], nx, ny);
		printf("ll0=[%d,%d], ur0=[%d,%d], nx0=%d, ny0=%d\n", ll0[0], ll0[1], ur0[0], ur0[1], nx0, ny0);
		printf("radius=%f, xcen=%f, ycen=%f\n", radius, xcen, ycen);
		
		
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
            varBfil = 0.;			// set to 0. since fill factor not computed
			
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
				
                // zero azm pointing North, zero inc pointing out from sun
                Bx[i] =-BmagBfilSinBinc * sin(Bazm);
                By[i] = BmagBfilSinBinc * cos(Bazm);
                Bz[i] = BmagBfil * cos(Binc);
				
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
		
		/* Create bitmap, with different methods depending on whether */
        /* it's a patch or full disk */
		
		if (harpflag) {		// harp
            bitmap = (int *)calloc(nxny, sizeof(int));
            /* bitmap no longer being used; all pixels treated the same */
            // only pixels in unpadded rectangle are disambiguated
            // XXXX add one extra pixel to account for one-sided differences?
            i0 = ll0[0] - ll[0];// i < ur[0] - ur0[0]; i++ 
            j0 = ll0[1] - ll[1];// j < ur[1] - ur0[1]; j++
            for (i = i0; i < nx - ur[0] + ur0[0]; i++) {
                for (j = j0; j < ny - ur[1] + ur0[1]; j++) {
					bitmap[i + j*nx] = 1;						// Question: what to do when we have noise mask?
                }
            }
        } else {			// full disk
            bitmap = (int *)calloc(nxny, sizeof(int));
            rad2 = radius * radius;
			if (useMask0) {						// Added Mar 27 2012
				for (i = 0; i < nxny; i++) {
					bitmap[i] = (sqrt(Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i]) > (noiseMask[i] * ksig + offset)) ? 1 : 0;
				}
			} else {  // Construct bitmap based on transverse field strength and distance from disk center.
				for (i = 0; i < nxny; i++) {
					j = i / nx;
					xx = i - j * nx - xcen + 1;
					yy = j - ycen + 1;
					r2 = xx * xx + yy * yy;
					if (r2 < rad2) {
						bitmap[i] = (sqrt(Bx[i] * Bx[i] + By[i] * By[i]) > bthresh1 + (bthresh0 - bthresh1) * sqrt(1. - r2 / rad2)) ? 1 : 0;
					} else {
						bitmap[i] = 0;
					}
				}
			}
            // Erode this bitmap to remove isolated above threshold pixels. 
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
        }		// end creating noise bitmap
		
		free(erodemap); free(growmap);
		
		probBa = (float *) calloc(nxny, sizeof(float));

		/* Before Bx, By, Bz change, find the radial acute solution */
		// Added Mar 21
		int index, index1;
		char *ambig_radial = (char *) calloc(nxny0, sizeof(char));		// solution for radial acute
		double bz0, bz1;		// two versions of radial, zonal and meridional field
		double pa = -1. * crota2 * DTOR;
		double lonc = drms_getkey_float(inRec, "CRLN_OBS", &status) * DTOR;
		double latc = drms_getkey_float(inRec, "CRLT_OBS", &status) * DTOR;
		double xx, yy;
		double lat, lon;
		double asd = asin(rsun_ref / dsun_obs);
		double rho, sinlat, coslat, sig, mu, chi;
		double a31, a32, a33;
		for (i = 0; i < nx0; i++) {
			for (j = 0; j < ny0; j++) {
				index = i + j * nx0;
				index1 = i + i0 + (j + j0) * nx;
				// find lat/lon
				xx = (ll0[0] + i - (crpix1 - 1)) / radius;		// radius computed earlier
				yy = (ll0[1] + j - (crpix2 - 1)) / radius;
				if (img2sphere(xx, yy, asd, latc, lonc, pa,
							   &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi))
					continue;		// off disk or something
				// from img2helioVector.c
				a31 =  cos(lat) * (sin(latc) * sin(pa) * cos(lon - lonc) + cos(pa) * sin(lon - lonc))
				- sin(lat) * cos(latc) * sin(pa);
				a32 = - cos(lat) * (sin(latc) * cos(pa) * cos(lon - lonc) - sin(pa) * sin(lon - lonc))
				+ sin(lat) * cos(latc) * cos(pa);
				a33 =  cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc);
				bz0 = a31 * Bx[index1] + a32 * By[index1] + a33 * Bz[index1];
				bz1 = - a31 * Bx[index1] - a32 * By[index1] + a33 * Bz[index1];
				ambig_radial[index] = (fabs(bz0) < fabs(bz1)) ? 1 : 0;
            }
        }

		
        /* This is the working part */
		
		printf("start\n"); fflush(stdout);
		
        ambig_(&geometry,
               &weak,
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
               &nap, &bthresh0, &bthresh1);
		
		
        /* Output data */
        outDims[0] = nx0; outDims[1] = ny0;	// Output array dimensions
		
		printf("nx=%d, ny=%d\n", nx0, ny0);
		printf("i0=%d, j0=%d\n", i0, j0);

		
		// Updated Oct 12, confidence is now char
        confidence = (char *) calloc(nxny0, sizeof(char));
        ambig_flag = (char *) calloc(nxny0, sizeof(char));
		
        /* Flag pixels for which the azimuth angle needs to be changed */
        for (i = 0; i < nx0; i++) {
            for (j = 0; j < ny0; j++) {
				index = i + j * nx0;
				index1 = i + i0 + (j + j0) * nx;
				// Updated Oct 12, now a char from 0 to 100
				// Updated Mar 21, still from 0 to 100, but 100 is good.
				confidence[index] = (char) (probBa[index1] * 100);
				ambig_flag[index] = (Bx[index1] > 0.) ? 1 : 0;
            }
        }
		
		/* The above solution is for potential only, now update the random solution as second bit */
		char pot;		// result from potential
		char r;
		srand((unsigned)time(NULL));	// might not be the best random number generator
		for (i = 0; i < nx0; i++) {
            for (j = 0; j < ny0; j++) {
				index = i + j * nx0;
				pot = (ambig_flag[index] % 2);
				if (confidence[index] > 51) {
					ambig_flag[index] += pot * 2;
				} else {
					r = rand() % 2;
					ambig_flag[index] += r * 2;
				}
			}
		}
		
		/* Update the radial acute solution as third bit */
		for (i = 0; i < nx0; i++) {
            for (j = 0; j < ny0; j++) {
				index = i + j * nx0;
				pot = (ambig_flag[index] % 2);
				if (confidence[index] > 51) {
					ambig_flag[index] += pot * 4;
				} else {
					ambig_flag[index] += ambig_radial[index] * 4;
				}
			}
		}
		free(ambig_radial);
		
		
        /* Bit maps */
        // Jun 22 Xudong
#if QUALMAP==1 
        qual_map = (int *) calloc(nxny0, sizeof(int));
        confid_map = (char *) calloc(nxny0, sizeof(char));
        for (i = 0; i < nx0; i++) {
            for (j = 0; j < ny0; j++) {
                // Jun 22 Xudong: update highest bit for now, should be lowest bit in final version
                qual_map[i + j * nx0] = inData_qual[i + j * nx0] | getAmbCode(confidence[i + j * nx0]);
                // Jun 22 Xudong: not updated for now, waiting for final scheme
                confid_map[i + j * nx0] = inData_conf[i + j * nx0];
            }
        }
        drms_free_array(inArray_qual); drms_free_array(inArray_conf);
#endif
		
        /* Stats */
        int nancount = 0, outsz = outDims[0] * outDims[1]; 
        for (i = 0; i < outsz; i++) {
            if (isnan(confidence[i])) nancount++;
        }
		
        /* Output segment */
        outArray_flag = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, ambig_flag, &status);
        if (status) DIE("Error creating flag array");
        outArray_prob = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, confidence, &status);
        if (status) DIE("Error creating prob array");
        outArray_qual = drms_array_create(DRMS_TYPE_INT, 2, outDims, qual_map, &status);
        if (status) DIE("Error creating qual array");
        outArray_conf = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, confid_map, &status);
        if (status) DIE("Error creating conf array");
		
        /* Output record */
        outRec = outRS->records[irec];
        outSeg_flag = drms_segment_lookup(outRec, "disambig");
        outSeg_prob = drms_segment_lookup(outRec, "conf_disambig");
#if QUALMAP==1 
        outSeg_qual = drms_segment_lookup(outRec, "info_map");
        outSeg_conf = drms_segment_lookup(outRec, "confid_map");
#endif
        for (i = 0; i < 2; i++) {
			outSeg_flag->axis[i] = outArray_flag->axis[i];	// For variable dimensions
			outSeg_prob->axis[i] = outArray_prob->axis[i];
#if QUALMAP==1
			outSeg_qual->axis[i] = outArray_qual->axis[i];
			outSeg_conf->axis[i] = outArray_conf->axis[i];
#endif
        }
        outArray_flag->parent_segment = outSeg_flag;
        outArray_prob->parent_segment = outSeg_prob;
#if QUALMAP==1
        outArray_qual->parent_segment = outSeg_qual;
        outArray_conf->parent_segment = outSeg_conf;
#endif
        // Compressed
		//        outArray_prob->bscale = 0.01;
		//        outArray_prob->bzero = 0.0;
		
        /* Result writing */
        status = drms_segment_write(outSeg_flag, outArray_flag, 0);
        if (status) DIE("Problem writing flag file");
        drms_free_array(outArray_flag);
        status = drms_segment_write(outSeg_prob, outArray_prob, 0);
        if (status) DIE("Problem writing prob file");
        drms_free_array(outArray_prob);
#if QUALMAP==1
        status = drms_segment_write(outSeg_qual, outArray_qual, 0);
        if (status) DIE("Problem writing qual file");
        drms_free_array(outArray_qual);
        status = drms_segment_write(outSeg_conf, outArray_conf, 0);
        if (status) DIE("Problem writing conf file");
        drms_free_array(outArray_conf);
#endif
        /* Set keywords */
        // Prime key
    	drms_copykey(outRec, inRec, "T_REC");
    	if (harpflag) {
    	    drms_copykey(outRec, inRec, "HARPNUM");
    	}
    	drms_copykey(outRec, inRec, "RUNNUM");
        // Info
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        // Parameters
        copy_me_keys(inRec, outRec);
        copy_geo_keys(inRec, outRec);
        if (harpflag) copy_patch_keys(inRec, outRec);
        drms_setkey_string(outRec, "BUNIT_026", " ");
        drms_setkey_string(outRec, "BUNIT_027", " ");
        // Disambiguation
        drms_setkey_int(outRec, "AMBPATCH", harpflag);	// Jun 22 Xudong
        drms_setkey_int(outRec, "AMBGMTRY", geometry);
        drms_setkey_int(outRec, "AMBWEAK", weak);
        drms_setkey_int(outRec, "AMBNEROD", nerode);
        drms_setkey_int(outRec, "AMBNGROW", ngrow);
        drms_setkey_int(outRec, "AMBNPAD", npad);
        drms_setkey_int(outRec, "AMBNAP", nap);
        drms_setkey_int(outRec, "AMBNTX", ntx);
        drms_setkey_int(outRec, "AMBNTY", nty);
        drms_setkey_int(outRec, "AMBSEED", seed);
        drms_setkey_int(outRec, "AMBNEQ", neq);
        drms_setkey_float(outRec, "AMBBTHR0", bthresh0);
        drms_setkey_float(outRec, "AMBBTHR1", bthresh1);
        drms_setkey_float(outRec, "AMBLMBDA", lambda);
        drms_setkey_float(outRec, "AMBTFCT0", tfac0);
        drms_setkey_float(outRec, "AMBTFCTR", tfactr);
        // stats
        drms_setkey_int(outRec, "TOTVALS", outsz);
        drms_setkey_int(outRec, "DATAVALS", outsz - nancount);
        drms_setkey_int(outRec, "MISSVALS", nancount);
        // Code version
		drms_setkey_string(outRec, "CODEVER5", "$Id: disambig.c,v 1.13 2012/09/07 03:25:13 xudong Exp $");
		drms_setkey_string(outRec, "AMBCODEV", ambcodev);
		// Maskinfo
		if (useMask) {
			drms_setkey_int(outRec, "USEMASK", 1);
//			drms_setkey_string(outRec, "MASKINFO", maskQuery);
			drms_setkey_int(outRec, "KSIG", ksig);
			drms_setkey_int(outRec, "C_NOISE", offset);
		} else {
			drms_setkey_int(outRec, "USEMASK", 0);
			drms_setkey_int(outRec, "KSIG", 0);
			drms_setkey_int(outRec, "C_NOISE", offset);
		}
		
        /* Set link */
        if (harpflag) {
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
		
		// Free mask
		if (noiseMask) free(noiseMask); noiseMask = NULL;
		
		/* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                   irec, wt - wt1, ct - ct1);
        }
		
        /* Clean up */
        free(Bx); free(By); free(Bz); free(dBt);  // ELW20100226 - mainly for sanity
        if (harpflag) free(inData_conf);
        if (bitmap != NULL) { free(bitmap); bitmap = NULL; }
		free(probBa);
		

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
