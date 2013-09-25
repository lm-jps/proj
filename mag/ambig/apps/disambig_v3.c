/*
 * Module name:		disambig.c
 *
 * Description:		Wrapper for the disambiguation module
 *
 * Original source:	Fortran wrapper by Graham and Igor (~graham/ME0/HMI/ambig.c)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * External function:
 *			ambig.f90
 *
 * Version:
 *			v1.0		Jan 25 2009
 *			v1.1		Mar 21 2012
 *			v1.2		Mar 27 2012
 *			v1.3		Apr 09 2012
 *			v2.0		Jun 04 2012
 *			v3.0		Sep 17 2012
 *			v3.1		Jun 04 2013
 *			v3.2		Jun 11 2013
 *      v3.3    Jun 26 2013
 *			v3.4		Jul 26 2013
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
 *			v3.0
 *			Rewritten by Xudong
 *			Changed bitmap and input for ambig.f90
 *			removed dBt, bthresh0,1, bitmap now generated in wrapper
 *			reinstate ambweak for possible NRT usage
 *			Q: fixed geometry argumenet & probBa freeing
 *			Q: xcen start from 0 for planar geo, 1 for spherical geo?
 *			v3.1
 *			Fixed bug in erode/grow
 *			Updated noise mask series name
 *			Updated geometry keyword
 *			v3.2
 *			Assigned offdisk values in conf_disambig to 0, using bitmap
 *      v3.3
 *      Fixed a but in getNoiseMask() and noiseMask.c, now return non-zero value if fails
 *			v3.4
 *			Take out radial acute and random solution for HARPs, fixed DATE and other keywords in an earlier checkin
 *
 *
 * Example:
 *  disambig_v3 "in=hmi_test.MEharp_720s_fd10[1420][2012.02.20_17:00:00_TAI]" "out=hmi_test.Bharp_720s_fd10"
 *  disambig_v3 "in=hmi_test.MEharp_720s_fd10[1420][2012.02.20_23:00:00_TAI]" "out=hmi_test.Bharp_720s_fd10" \
		"AMBGMTRY=2" "AMBBTHR0=200." "AMBBTHR1=300." "AMBTFCTR=0.99" "AMBNEQ=10" -l
 *
 */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "copy_me_keys.c"
#include "cartography.c"
#include "noisemask.c"
#include "timing.c"

#define PI	(M_PI)
#define	DTOR	(PI / 180.)
#define Rad2arcsec	(3600. * 180. / PI)
#define FOURK (4096)
#define FOURK2	(16777216)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define BB(bx,by,bz) (sqrt(bx*bx+by*by+bz*bz))
#define BH(bx,by) (sqrt(bx*bx+by*by))

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
#define ambCode0 0x0000 // 000, not disambiguated
#define ambCode1 0x0001 // 001, weak not annealed (pf/radial acute/random)
#define ambCode2 0x0003 // 011, weak annealed
#define ambCode3 0x0007 // 111, strong annealed

// switches
#define QUALMAP 1		// QUALMAP not updated if 0

// Switches
//#define OLD 1			// Old version of ambig_()
//#define OLDBITMAP 1	// Old bitmap
//#define WRITEBITMAP 1		// Output bitmap, etc. for test purpose
#define NO_DISAMBIG 1

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

// ============================================

// Get ephemeris info, note only three params are used: radius, center location

int getEphemeris(DRMS_Record_t *inRec, float *radius, float *pix_x, float *pix_y);

// Obtain size, keywords, etc.

int getGeometry(DRMS_Record_t *inRec, int npad, int *geometry,
				int *ll, int *ur, int *ll0, int *ur0,
				float *radius, float *xcen, float *ycen);

// Obtain data arrays

int getData(DRMS_Record_t *inRec, float *Bx, float *By, float *Bz, float *dBt,
			int *ll, int *ur, int *ll0, int *ur0);

// Create bitmap

int getBitmap(DRMS_Record_t *inRec, int *bitmap, float *Bx, float *By, float *Bz,
			  int geometry, int *useMask, int nrt, int *offset, int nerode, int ngrow,
			  int *ll, int *ur, int *ll0, int *ur0, float bthresh0, float bthresh1, char *maskQuery);		// offset changed to pointer Jun 27 Xudong

// Get noise esimation from Yang's code

int getNoiseMask(DRMS_Record_t *inRec, float *noiseMask, int *ll, int *ur, char *maskQuery);

// Get noise esimation using a linear threshold

int getLinearMask(DRMS_Record_t *inRec, float *noiseMask, int *ll, int *ur, 
				  float bthresh0, float bthresh1);

// Get radial acute solution

void getRadialAcute(DRMS_Record_t *inRec, char *ambig_r, float *Bx, float *By, float *Bz,
					int *ll, int *ur, int *ll0, int *ur0);

// Erode and grow the bitmap to clean up isolated pixels

void erodeAndGrow(DRMS_Record_t *inRec, int *bitmap, int nerode, int ngrow,
				  int *ll, int *ur, int *ll0, int *ur0);

// Gather and write out arrays

int writeData(DRMS_Record_t *outRec, DRMS_Record_t *inRec,
			  float *Bx, float *probBa, char *ambig_r,
			  int *ll, int *ur, int *ll0, int *ur0, int harpflag);		// Jul 26 2013 XS, added harpflag

// Output auxiliary arrays

int writeAuxData(DRMS_Record_t *outRec, DRMS_Record_t *inRec, 
				 float *Bx, float *By, float *Bz, float *probBa, float *dBt, int *bitmap,
				 int *ll, int *ur);

// ============================================

/* ##### Prototypes for external Fortran functions ##### */

#ifdef OLD
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
#else
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
				   float *Bx, float *By, float *Bz, float *probBa,
				   int *bitmap,
				   float *radius,
				   int *nap);
#endif

//====================

char *module_name = "disambig_v3";	/* Module name */
char *version_id = "2013 Jul 26";  /* Version number */

#ifdef OLD
char *segName[] = {"field", "inclination", "azimuth", "alpha_mag", 
					"field_err", "inclination_err", "alpha_err", 
					"field_inclination_err", "field_alpha_err", "inclination_alpha_err"};
#else
char *segName[] = {"field", "inclination", "azimuth", "alpha_mag"};
#endif

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=minimal messages; 2=all messages"},
	{ARG_FLAG, "l", "", "Flag to use linear noise threshold (AMBTHRn)"},
	{ARG_FLAG, "n", "", "Flag for nrt series"},
    {ARG_INT, "OFFSET", "20", "Constant added to the noise mask"},
    {ARG_INT, "AMBGMTRY", "0", "0 for automatic selection; 1 for planar; 2 for spherical."},
	{ARG_INT, "AMBWEAK", "1", "0 for random; 1 for potential field; 2 for most radial."},	// now defunct
    {ARG_INT, "AMBNEROD", "1", "Number of pixels by which to erode map of above threshold pixels."},
    {ARG_INT, "AMBNGROW", "5", "Number of pixels by which to grow eroded map."},
    {ARG_INT, "AMBNPAD", "20", "Pixel number to pad with zeros for potential field."},
    {ARG_INT, "AMBNAP", "10", "Pixel number to apodize for potential field."},
    {ARG_INT, "AMBNTX", "20", "Tile number in x (lon) direction for pf on a sphere."},
    {ARG_INT, "AMBNTY", "20", "Tile number in y (lat) direction for pf on a sphere."},
    {ARG_FLOAT, "AMBBTHR0", "2.0e2", "Threshold field strength."},
    {ARG_FLOAT, "AMBBTHR1", "4.0e2", "Threshold field strength."},
    {ARG_INT, "AMBSEED", "1", "Input random number seed (seed>0)."},
    {ARG_INT, "AMBNEQ", "20", "Num of reconfigurations attempted at each temperature setting."},
    {ARG_FLOAT, "AMBLMBDA", "1.", "Weighting factor between div. and vert. current density."},
    {ARG_FLOAT, "AMBTFCT0", "2.", "Input factor to scale initial temperature (tfac0>0)."},
    {ARG_FLOAT, "AMBTFCTR", "0.990", "Input factor to reduce temperature (0<tfactr<1)."},
    {ARG_END}
};

/* ################## Main Module ################## */

int DoIt(void)
{
	
	int status = DRMS_SUCCESS;
	
	// Time measuring
	
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);
	
	// Versioning
	
	char ambcodev[50];
	sprintf(ambcodev,"%s %s", module_name, version_id);
	
	// Get parameters
	
    char *inQuery = (char *)params_get_str(&cmdparams, "in");
    char *outQuery = (char *)params_get_str(&cmdparams, "out");
    int verbflag = params_get_int(&cmdparams, "VERB");
    int geometry = params_get_int(&cmdparams, "AMBGMTRY");
	int useMask = !(params_isflagset(&cmdparams, "l"));		// if not set, use default noise mask
	int nrt = params_isflagset(&cmdparams, "n");			// nrt series
    int offset = params_get_int(&cmdparams, "OFFSET");
	
	int weak = params_get_int(&cmdparams, "AMBWEAK");		// for possible nrt usage
    int nerode = params_get_int(&cmdparams, "AMBNEROD");
    int ngrow = params_get_int(&cmdparams, "AMBNGROW");
    int npad = params_get_int(&cmdparams, "AMBNPAD");
    int nap = params_get_int(&cmdparams, "AMBNAP");
    int ntx = params_get_int(&cmdparams, "AMBNTX");
    int nty = params_get_int(&cmdparams, "AMBNTY");
    int seed = params_get_int(&cmdparams, "AMBSEED");
    int neq = params_get_int(&cmdparams, "AMBNEQ");
    float bthresh0 = params_get_float(&cmdparams, "AMBBTHR0");
    float bthresh1 = params_get_float(&cmdparams, "AMBBTHR1");
    float lambda = params_get_float(&cmdparams, "AMBLMBDA");
    float tfac0 = params_get_float(&cmdparams, "AMBTFCT0");
    float tfactr = params_get_float(&cmdparams, "AMBTFCTR");
	
	// Check parameters
	
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
	
	// Open input
	
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    int nrecs = inRS->n;
	
	// Do this for each record */
    for (int irec = 0; irec < nrecs; irec++) {
		
		// Measure time
		
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }
		
		// Input record, prime keys, code modes
		
        DRMS_Record_t *inRec = inRS->records[irec];
		TIME t_rec = drms_getkey_time(inRec, "T_REC", &status);
		int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);	// HARP number
		int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
		int useMask_t = useMask;		// using reconstructed mask image
		int geometry_t = geometry;
		
		// Get geometry
				
		int ll[2], ur[2];		// lower-left and upper-right corner after padding	
		int ll0[2], ur0[2];		
		float radius, xcen, ycen;
		
		if (getGeometry(inRec, npad, &geometry_t, 
						ll, ur, ll0, ur0,
						&radius, &xcen, &ycen)) {		// if failed
			printf("Getting geometry failed, image #%d skipped.\n", irec);
			continue;
		}
		
		int nx = ur[0] - ll[0] + 1;			// size with padding
		int ny = ur[1] - ll[1] + 1;
		int nx0 = ur0[0] - ll0[0] + 1;		// without padding
		int ny0 = ur0[1] - ll0[1] + 1;
		int nxny = nx * ny, nxny0 = nx0 * ny0;
		
		printf("ll=[%d,%d], ur=[%d,%d], nx=%d, ny=%d\n", ll[0], ll[1], ur[0], ur[1], nx, ny);
		printf("ll0=[%d,%d], ur0=[%d,%d], nx0=%d, ny0=%d\n", ll0[0], ll0[1], ur0[0], ur0[1], nx0, ny0);
		printf("radius=%f, xcen=%f, ycen=%f\n", radius, xcen, ycen);
		
		printf("Geometry determined. \n");
		
		// Get data
		
		float *Bx = (float *) (malloc(nxny * sizeof(float)));
        float *By = (float *) (malloc(nxny * sizeof(float)));
        float *Bz = (float *) (malloc(nxny * sizeof(float)));
		float *dBt = (float *) (malloc(nxny * sizeof(float)));
		
		if (getData(inRec, Bx, By, Bz, dBt, ll, ur, ll0, ur0)) {	// if failed
			free(Bx); free(By); free(Bz); free(dBt);
			printf("Getting data failed, image #%d skipped.\n", irec);
			continue;
		}
		
		printf("Arrays read. \n");
 
		// Get bitmap, noise mask used inside getBitmap()
		
		int *bitmap = (int *) (calloc(nxny, sizeof(int)));
		
		int offset_t = offset;
		char maskQuery[100];        // query name for noise mask, sep 23
        
		if (getBitmap(inRec, bitmap, Bx, By, Bz, 
					  geometry, &useMask_t, nrt, &offset_t, nerode, ngrow,
					  ll, ur, ll0, ur0, bthresh0, bthresh1, maskQuery)) {			// Offset changed to pointer Jun 27 2013 Xudong
			free(bitmap);
			free(Bx); free(By); free(Bz);
			printf("Creating bitmap failed, image #%d skipped.\n", irec);
			continue;
		}
 
		printf("Bitmap created. \n");
		
		// Obtain radial acute solution before Bx changes
		
		char *ambig_r = NULL;
		if (!harpflag) {		// Skipped for HARPs, Jul 26
		  ambig_r = (char *) calloc(nxny0, sizeof(char));
  		getRadialAcute(inRec, ambig_r, Bx, By, Bz, ll, ur, ll0, ur0);
  		printf("Radial acute solution determined\n");
  	}
		
		// This is the working part
		
		printf("start...\n");
		
		float *probBa = (float *) (malloc(nxny * sizeof(float)));
		
		printf("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
			   geometry_t, nx, ny, weak, npad, seed, neq, ntx, nty, nap);
		printf("%f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
			   xcen, ycen, radius, lambda, tfactr, tfac0, radius, bthresh0, bthresh1);

#ifndef NO_DISAMBIG
#ifdef OLD
		ambig_(&geometry_t,
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
#else
        ambig_(&geometry_t,
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
               bitmap,
               &radius,
               &nap);
#endif
#endif
		
		printf("done\n");
 
		
		// Output, one record each time
		
		DRMS_Record_t *outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
		if (status) {		// if failed
			printf("Creating output record failed, image #%d skipped.\n", irec);
			continue;
		}

#ifdef WRITEBITMAP				// output bitmap, etc. for testing purpose
		
		if (writeAuxData(outRec, inRec, Bx, By, Bz, probBa, dBt, bitmap, ll, ur)) {
			printf("Output test arrays error, record skipped\n");
			continue;
		}
		
#else
		
		if (writeData(outRec, inRec, Bx, probBa, ambig_r, ll, ur, ll0, ur0, harpflag)) {		// Jul 26, added harpflag
			printf("Output error, record skipped\n");
			free(Bx); free(By); free(Bz);
			free(bitmap); free(probBa); free(ambig_r);
			continue;
		}
		
#endif
		
		// Keywords
		
		// =======================================
		drms_copykey(outRec, inRec, "T_REC");
    	if (harpflag) drms_copykey(outRec, inRec, "HARPNUM");
        // Parameters
        copy_me_keys(inRec, outRec);
        copy_geo_keys(inRec, outRec);
        if (harpflag) copy_patch_keys(inRec, outRec);
        drms_setkey_string(outRec, "BUNIT_026", " ");
        drms_setkey_string(outRec, "BUNIT_027", " ");
        drms_setkey_int(outRec, "DOFFSET", offset_t);		// Jun 27 Xudong
		// Info
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        char timebuf[1024];
        double UNIX_epoch = -220924792.0;	/* 1970.01.01_00:00:00_UTC */
        sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
        drms_setkey_string(outRec, "DATE", timebuf);
        // Disambiguation
        drms_setkey_int(outRec, "AMBPATCH", harpflag);	// Jun 22 Xudong
        drms_setkey_int(outRec, "AMBGMTRY", geometry_t);	// Jun 04 2013
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
        // Code version
		drms_setkey_string(outRec, "CODEVER5", "$Id: disambig_v3.c,v 1.9 2013/09/25 20:56:42 xudong Exp $");
		drms_setkey_string(outRec, "AMBCODEV", ambcodev);
		// Maskinfo
		if (useMask_t) {            // Sep 25, changed to useMask_t, NOISEMASK
			drms_setkey_int(outRec, "NOISEMASK", 1);
            drms_setkey_string(outRec, "MASKINFO", maskQuery);
		} else {
			drms_setkey_int(outRec, "NOISEMASK", 0);
            drms_setkey_string(outRec, "MASKINFO", "Linear");
		}
		// For debug
		drms_setkey_string(outRec, "CODENAME", module_name);

		// =======================================
		
		// Link
		
		DRMS_Link_t *link = NULL;
		if (harpflag) {
            link = hcon_lookup_lower(&outRec->links, "PATCH");
            if (link) drms_link_set("PATCH", outRec, inRec);
        } else {
            link = hcon_lookup_lower(&outRec->links, "MDATA");
            if (link) drms_link_set("MDATA", outRec, inRec);
        }
		
		// Clean up
		
		drms_close_record(outRec, DRMS_INSERT_RECORD);

#ifndef WRITEBITMAP
		free(Bx); free(By); free(Bz); free(dBt);
		free(probBa); free(bitmap);
#endif
		
		if (harpflag) free(ambig_r);		// Jul 26
		
		// Time measure
		
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                   irec, wt - wt1, ct - ct1);
        }
		
		
	} // irec
	
	drms_close_records(inRS, DRMS_FREE_RECORD);
	
	if (verbflag) {
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
			   wt - wt0, ct - ct0);
    }
	
	return 0;
	
}	// main


// ======================================================================

// Get ephemeris info, note only three params are used: radius, center location

int getEphemeris(DRMS_Record_t *inRec, float *radius, float *pix_x, float *pix_y)
{
	
	int status = 0;
	
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
	
	float crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation
	float sina = sin(crota2 * DTOR); 
	float cosa = cos(crota2 * DTOR);
	float cdelt = drms_getkey_float(inRec, "CDELT1", &status);	// in arcsec, assuming dx=dy
	float rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
	if (status) rsun_ref = 6.96e8;
	float dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status);
	
	*radius = asin(rsun_ref / dsun_obs) * Rad2arcsec / cdelt;		// in pixel
	
	float crvalx, crvaly, crpix1, crpix2;
	if (harpflag) {
		crvalx = drms_getkey_float(inRec, "IMCRVAL1", &status);	// center of solar disc in arcsec
		crvaly = drms_getkey_float(inRec, "IMCRVAL2", &status);
		crpix1 = drms_getkey_float(inRec, "IMCRPIX1", &status);	// disk center in ccd, original
		crpix2 = drms_getkey_float(inRec, "IMCRPIX2", &status);
	} else {
		crvalx = drms_getkey_float(inRec, "CRVAL1", &status);	// center of solar disc in arcsec
		crvaly = drms_getkey_float(inRec, "CRVAL2", &status);
		crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);	// disk center in ccd, original
		crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
	}
	
	*pix_x = PIX_X(0.0,0.0);
	*pix_y = PIX_Y(0.0,0.0);
	
	return 0;
	
}

// ======================================================================

// Obtain size, keywords, etc.

int getGeometry(DRMS_Record_t *inRec, int npad, int *geometry,
				int *ll, int *ur, int *ll0, int *ur0,
				float *radius, float *xcen, float *ycen)
{
	
	int status = 0;
	int nx, ny, nx0, ny0;
	
	// Ephemeris
	
	float pix_x, pix_y;		// PIX_X, PIX_Y
	if (getEphemeris(inRec, radius, &pix_x, &pix_y)) {		// no error checking for now
		SHOW("Getting ephemeris error\n");
		return 1;
	}
	
	// Start
	
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
	
	if (harpflag) {		// harp
		
		// If unspecified, set geometry based on patch size
		
		if (*geometry == 0) {
			float minlon = drms_getkey_float(inRec, "LON_MIN", &status);
			float minlat = drms_getkey_float(inRec, "LAT_MIN", &status);
			float maxlon = drms_getkey_float(inRec, "LON_MAX", &status);
			float maxlat = drms_getkey_float(inRec, "LAT_MAX", &status);
			if (fabs(maxlon - minlon) < 20. && fabs(maxlat - minlat) < 20.) {
				*geometry = 1;
				printf("Geometry selection: using planar geometry\n");
			} else {
				*geometry = 2;
				printf("Geometry selection: using spherical geometry\n");
			}
		}
		
		// Mask array for size
		
		DRMS_Segment_t *maskSeg = drms_segment_lookup(inRec, "bitmap");		// from segment, for patch
		nx0 = maskSeg->axis[0]; ny0 = maskSeg->axis[1];		// Get array dimensions from mask array.
			
		// Coordinates of cut out rectangle
		
		ll0[0] = drms_getkey_int(inRec, "CRPIX1", &status) - 1;	// lower left corner
		ll0[1] = drms_getkey_int(inRec, "CRPIX2", &status) - 1;
		ur0[0] = ll0[0] + nx0 - 1;					// upper right corner
		ur0[1] = ll0[1] + ny0 - 1;
		
		// Pad the actual patch
		
		nx = nx0 + 2 * npad;
		ny = ny0 + 2 * npad;
		
		ll[0] = ll0[0] - npad;
		ll[1] = ll0[1] - npad;
		ur[0] = ur0[0] + npad;
		ur[1] = ur0[1] + npad;
		
		// Ensure that this stays within the image
		// Assuming original harp within image
		
		if (ll[0] < 0) {
			nx += ll[0];
			ll[0] = 0;
		}
		if (ll[1] < 0) {
			ny += ll[1];
			ll[1] = 0;
		}
		if (ur[0] > 4095) {
			nx -= (ur[0] - 4095);		// fixed a bug (originially +=)
			ur[0] = 4095;
		}
		if (ur[1] > 4095) {
			ny -= (ur[1] - 4095);
			ur[1] = 4095;
		}
		
		// xcen, ycen
		// center of patch relative to disk center (start from 0)
		if (*geometry == 1) {	
			*xcen = (float)ll0[0] + 0.5 * ((float)(nx0) - 1.) - (pix_x - 1.);
			*ycen = (float)ll0[1] + 0.5 * ((float)(ny0) - 1.) - (pix_y - 1.);
		}
		
		// disk center relative to lower left corner of padded patch (start from 1)
		if (*geometry == 2) {	
			*xcen = pix_x - (float)ll[0];
			*ycen = pix_y - (float)ll[1];
		}
		
		
	} else {		// full disk
		
		*geometry = 2;
		
		*xcen = pix_x;	// start from 1
		*ycen = pix_y;
		
		nx = ny = nx0 = ny0 = FOURK;
		ll[0] = ll[1] = 0;
		ur[0] = ur[1] = FOURK - 1;
		ll0[0] = ll0[1] = 0;
		ur0[0] = ur0[1] = FOURK - 1;
		
	}	// harpflag
	
//	printf("nx=%d, ny=%d, nx0=%d, ny0=%d\n", nx, ny, nx0, ny0);
	
	return 0;
	
}


// ======================================================================

// Obtain data arrays

int getData(DRMS_Record_t *inRec, float *Bx, float *By, float *Bz, float *dBt,
			int *ll, int *ur, int *ll0, int *ur0)
{
	
	int status = 0;
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	int nx0 = ur0[0] - ll0[0] + 1;
	int ny0 = ur0[1] - ll0[1] + 1;
	int nxny = nx * ny, nxny0 = nx0 * ny0;
	
	// Read out from input record
	
	int nseg = ARRLENGTH(segName);
	DRMS_Segment_t *inSeg;
	DRMS_Array_t *inArray[nseg], *inArray_qual, *inArray_conf;
	float *inData[nseg];
	int *inData_qual;
	char *inData_conf;
	
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
	
	if (harpflag) {		// harp
		
		for (int i = 0; i < nseg; i++) {
			DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, segName[i]);
			inArray[i] = drms_segment_readslice(inSeg, DRMS_TYPE_FLOAT, ll, ur, &status);
			if (status) {
				for (int j = 0; j < i; j++) drms_free_array(inArray[i]);
				printf("Segment reading error\n");
				return 1;
			}
			inData[i] = (float *) inArray[i]->data;
		}
		
	} else {		// full disk
		
		for (int i = 0; i < nseg; i++) {
			inSeg = drms_segment_lookup(inRec, segName[i]);
			inArray[i] = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
			if (status) {
				for (int j = 0; j <= i; j++) drms_free_array(inArray[i]);
				printf("Segment reading error\n");
				return 1;
			}
			inData[i] = (float *) inArray[i]->data;
		}
		
	}
	
	// Convert to Bx, By, Bz, copy out qual & conf
	// Changed to double for precision
	
	double Bmag, Binc, Bazm, Bfil;
	double BmagBfil, CosBinc, SinBinc;
	
	double varBmag, varBinc, varBfil;
	double varBmagBinc, varBfilBmag, varBfilBinc, BmagBfilSinBinc;
	
	for (int i = 0; i < nxny; i++) {
		
		Bmag = inData[0][i];			// "field"
		Binc = inData[1][i] * DTOR;		// "inclination"
		Bazm = inData[2][i] * DTOR;		// "azimuth"
		Bfil = inData[3][i];			// "alpha_mag"
		
		BmagBfil = Bmag * Bfil;
		CosBinc = cos(Binc);
		SinBinc = sin(Binc);
		BmagBfilSinBinc = BmagBfil * SinBinc;
		
		if (isnan(Bmag) || isnan(Binc) || isnan(Bazm) || isnan(Bfil)) {
			
			Bx[i] = 0.; By[i] = 0.; Bz[i] = 0.;
			
		} else {
			
			Bx[i] = - BmagBfilSinBinc * sin(Bazm);
			By[i] = BmagBfilSinBinc * cos(Bazm);
			Bz[i] = BmagBfil * CosBinc;
			
		}
		
#ifdef OLD
		
		// Variances
		varBmag = inData[4][i];	    	// "field_err"
		varBinc = inData[5][i];		// "inclination_err"
		//varBfil = inData[6][i];		// "alpha_err"
		varBfil = 0.;			// set to 0. since fill factor not computed
		
		// Covariances
		varBmagBinc = inData[7][i];		// "field_inclination_err"
		varBfilBmag = inData[8][i];		// "field_alpha_err"
		varBfilBinc = inData[9][i];		// "inclination_alpha_err"

		if (isnan(Bmag) || isnan(Binc) || isnan(Bazm) || isnan(Bfil)) {
			dBt[i] = 1.E9;
		} else {
			dBt[i] = (BmagBfil==0) ? 1.E9 : 
			(BmagBfil * 
			 sqrt (SinBinc * SinBinc * (varBfil / (Bfil * Bfil) + varBmag / (Bmag * Bmag)
										+ 2.*varBfilBmag/BmagBfil)
				   + varBinc * CosBinc * CosBinc
				   + (varBfilBinc / Bfil + varBmagBinc / Bmag) * SinBinc * CosBinc));
		}
		
//		dBt[i] = 0.;

#endif
		
	}
	
	for (int i = 0; i < nseg; i++) drms_free_array(inArray[i]);
	
	//
	
	return 0;
	
}

// ======================================================================

// Create bitmap

int getBitmap(DRMS_Record_t *inRec, int *bitmap, float *Bx, float *By, float *Bz,
			  int geometry, int *useMask, int nrt, int *offset, int nerode, int ngrow,
			  int *ll, int *ur, int *ll0, int *ur0, float bthresh0, float bthresh1, char *maskQuery)
{
	
	int status = 0;
	
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	int nx0 = ur0[0] - ll0[0] + 1;
	int ny0 = ur0[1] - ll0[1] + 1;
	int nxny = nx * ny, nxny0 = nx0 * ny0;
	
	int i0 = ll0[0] - ll[0];	// offset
	int j0 = ll0[1] - ll[1];
	
	float *noiseMask = (float *) (malloc(nxny * sizeof(float)));
	
	float offset_t;		// constant to be added on noise estimation

#ifdef OLDBITMAP
	
	printf("%d, %d, %d, %d\n", i0, j0, i0 + nx0, j0 + ny0);
	
	int idx;
	if (harpflag) {
		for (int i = i0; i < (i0 + nx0); i++) {
			for (int j = j0; j < (j0 + ny0); j++) {
				idx = j * nx + i;
				bitmap[idx] = 1;
			}
		}
	} else {
		// Test only
	}
	
#else
	
	// Get noise estimation and thresholding

	if (*useMask &&
		!(getNoiseMask(inRec, noiseMask, ll, ur, maskQuery))) {	// get Yang's noise estimate
        
		offset_t = *offset;		// specified by argument
		SHOW("Using noise mask estimate\n");
		
	} else {
		
		*useMask = 0;	// originally 0, or failed to fetch Yang's noise estimate
		
		if (getLinearMask(inRec, noiseMask, ll, ur, 
						  bthresh0, bthresh1)) {	// get linear noise estimate
			SHOW("Error getting error estimate\n");
			return 1;
		}
		
		offset_t = 0;
		SHOW("Using linear noise estimate\n");
		
	}
	
	// Create bitmap
	
	printf("offset=%f\n", offset_t);

	for (int i = 0; i < nxny; i++) {
		if (isnan(noiseMask[i])) {
			bitmap[i] = 0;
		} else {
			bitmap[i] = ((BH(Bx[i],By[i])) > (noiseMask[i] + offset_t)) ? 2 : 0;	// new def
		}
	}
	
	// Erode and grow
	
	erodeAndGrow(inRec, bitmap, nerode, ngrow, ll, ur, ll0, ur0);
	
#endif
	
	*offset = offset_t;		// Jun 27 Xudong
	free(noiseMask);
	
	//
	
	return 0;

}

// ======================================================================

// Get noise esimation from Yang's code

int getNoiseMask(DRMS_Record_t *inRec, float *noise, int *ll, int *ur, char *maskQuery)
{
	
	int status = 0;
	
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	
	//
	
	TIME tobs = drms_getkey_time(inRec, "T_OBS", &status);
	float vr = drms_getkey_float(inRec, "OBS_VR", &status);
	
	double *noise_fd = (double *) (malloc(FOURK2 * sizeof(double)));	// full disk
	
	int xDim = FOURK, yDim = FOURK;
	float r, x0, y0;
	if (getEphemeris(inRec, &r, &x0, &y0)) {
		SHOW("Getting ephemeris error\n");
		return 1;
	}
	x0 -= 1; y0 -= 1;
	
	status = noisemask(tobs, xDim, yDim, x0, y0, r, vr, noise_fd, maskQuery);
    
    if (status) {               // Added Jun 26 2013, Xudong
        free(noise_fd);
        return status;
    }
	
	// Cutout
	
	int idx, idx_fd;
	int full_disk = ((nx == FOURK) && (ny == FOURK));
	float r_thresh = 0.995 * r;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			idx = j * nx + i;
			idx_fd = (j + ll[1]) * FOURK + (i + ll[0]);
			noise[idx] = noise_fd[idx_fd];
			if (full_disk && (hypot(i * 1.0 + ll[0] - x0, j * 1.0 + ll[1] - y0) > r_thresh))
				noise[idx] = 1.e4;
		}
	}
	
	free(noise_fd);
	
	return 0;
	
}

// ======================================================================

// Get noise esimation using a linear threshold
// More like cosine...

int getLinearMask(DRMS_Record_t *inRec, float *noiseMask, int *ll, int *ur,
				  float bthresh0, float bthresh1)
{
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	int nxny = nx * ny;
	
	// Get radius (in pixel), xcen, ycen (start from 1)
	
	float radius, xcen, ycen;
	if (getEphemeris(inRec, &radius, &xcen, &ycen)) {
		SHOW("Getting ephemeris error\n");
		return 1;
	}
	xcen -= 1; ycen -= 1;		// start from 0
	float rad2 = radius * radius;
	
	// Noise estimation
	
	int idx;
	float xx, yy, mu2;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			idx = j * nx + i;
			xx = i + ll[0] - xcen;		// full disk has ll=[0,0]
			yy = j + ll[1] - ycen;
			mu2 = (xx * xx + yy * yy) / rad2;
			if (mu2 < 1) {
				noiseMask[idx] = bthresh1 + (bthresh0 - bthresh1) * sqrt(1. - mu2);
			} else {
				noiseMask[idx] = DRMS_MISSING_FLOAT;
			}
		}
	}
	
	return 0;
}

// ======================================================================

// Get radial acute solution

void getRadialAcute(DRMS_Record_t *inRec, char *ambig_r, float *Bx, float *By, float *Bz,
					int *ll, int *ur, int *ll0, int *ur0)
{
	
	int status = 0;
	
	int nx = ur[0] - ll[0] + 1;
	int nx0 = ur0[0] - ll0[0] + 1;
	int ny0 = ur0[1] - ll0[1] + 1;
	
	int i0 = ll0[0] - ll[0];	// offset
	int j0 = ll0[1] - ll[1];
	
	float r, x0, y0;
	if (getEphemeris(inRec, &r, &x0, &y0)) {
		SHOW("Getting ephemeris error\n");
		return;
	}
	x0 -= 1; y0 -= 1;
	
	float crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation	
	double lonc = drms_getkey_float(inRec, "CRLN_OBS", &status) * DTOR;
	double latc = drms_getkey_float(inRec, "CRLT_OBS", &status) * DTOR;
	double rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
	if (status) rsun_ref = 6.96e8;
	double dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status);
	
	double pa = -1. * crota2 * DTOR;
	double asd = asin(rsun_ref / dsun_obs);
	
	double xx, yy;
	double lat, lon;
	double rho, sinlat, coslat, sig, mu, chi;
	double a31, a32, a33;
	
	double bz0, bz1;		// two versions of radial, zonal and meridional field
	int idx, idx1;
	
	for (int i = 0; i < nx0; i++) {
		for (int j = 0; j < ny0; j++) {
			idx = i + j * nx0;
			idx1 = i + i0 + (j + j0) * nx;
			// find lat/lon
			xx = (ll0[0] + i - x0) / r;		// radius computed earlier
			yy = (ll0[1] + j - y0) / r;
			if (img2sphere(xx, yy, asd, latc, lonc, pa,
						   &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi))
				continue;		// off disk or something
			// from img2helioVector.c
			a31 = cos(lat) * (sin(latc) * sin(pa) * cos(lon - lonc) + cos(pa) * sin(lon - lonc))
					- sin(lat) * cos(latc) * sin(pa);
			a32 = - cos(lat) * (sin(latc) * cos(pa) * cos(lon - lonc) - sin(pa) * sin(lon - lonc))
					+ sin(lat) * cos(latc) * cos(pa);
			a33 = cos(lat) * cos(latc) * cos(lon - lonc) + sin(lat) * sin(latc);
			//
			bz0 = a31 * Bx[idx1] + a32 * By[idx1] + a33 * Bz[idx1];
			bz1 = - a31 * Bx[idx1] - a32 * By[idx1] + a33 * Bz[idx1];
			//
			ambig_r[idx] = (fabs(bz0) < fabs(bz1)) ? 1 : 0;
		}
	}
	
	//
	
	return;
		
}


// ======================================================================

// Erode and grow the bitmap to clean up isolated pixels
// In harp mode, the bitmap is not grown back

void erodeAndGrow(DRMS_Record_t *inRec, int *bitmap, int nerode, int ngrow,
				  int *ll, int *ur, int *ll0, int *ur0)
{
	
	int status = 0;
	
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	int harpflag = (harpnum == DRMS_MISSING_INT) ? 0 : 1;		// full disk vs harp
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	int nxny = nx * ny;
	int nx0 = ur0[0] - ll0[0] + 1;
	int ny0 = ur0[1] - ll0[1] + 1;
	int nxny0 = nx0 * ny0;
	
	int i, j, l, m;
	
	printf("nerode=%d, ngrow=%d\n", nerode, ngrow);
	
	// Erode this bitmap to remove isolated above threshold pixels.
	
	int nxe = nx + 2 * nerode;
	int nye = ny + 2 * nerode;
	int nxnye = nxe * nye;
	
	int *erodemap = (int *)calloc(nxnye, sizeof(int));
	for (i = 0; i < nxnye; i++) erodemap[i] = 1;
	
	for (i = 0; i < nxny; i++) {
		if (bitmap[i] == 0) {
			j = i + nerode * (nx + 1 + 2 * (i / nx + nerode));
			for (l = -nerode; l <= nerode; l++) {		// 2013.05.16
				for (m = -nerode; m <= nerode; m++) {
					erodemap[j + l + m * (nx + 2 * nerode)] = 0;
				}
			}
		}
	}
	
	// For HARP, set all remaining pixels in the bounding box to 1, paddings to 0, then return
	// For full disk, set all eroded pixels to 0
	
	int idx, idx1;
	
	if (harpflag) {
		
		int i0 = ll0[0] - ll[0], i1 = ur0[0] - ll[0];
		int j0 = ll0[1] - ll[1], j1 = ur0[1] - ll[1];
		printf("i0=%d,j0=%d,i1=%d,j1=%d\n", i0, j0, i1, j1);
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				idx = j * nx + i;
				idx1 = (j + nerode) * nxe + (i + nerode);
				if (erodemap[idx1] == 0) { bitmap[idx] = 1; }
				if (i < i0 || i > i1 || j < j0 || j > j1) { bitmap[idx] = 0; }
			}
		}
		printf("here\n");
		free(erodemap);
		return;
		
	} else {
		
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				idx = j * nx + i;
				idx1 = (j + nerode) * nxe + (i + nerode);
				if (erodemap[idx1] == 0) { bitmap[idx] = 0; }
			}
		}
		
	}

	// Grow for full disk mode
	
	int nxg = nxe + 2 * ngrow;
	int nyg = nye + 2 * ngrow;
	int nxnyg = nxg * nyg;
	int *growmap = (int *)calloc(nxnyg, sizeof(int));
	
	for (i = 0; i < nxnye; i++) {
		if (erodemap[i]) {
			j = i + ngrow * (nxe + 1 + 2 * (i / nxe + ngrow));
			for (l = -ngrow; l <= ngrow; l++) {		// 2013.05.16
				for (m = -ngrow; m <= ngrow; m++) {
					growmap[j + l + m * (nxe + 2 * ngrow)] = 1;
				}
			}
		}
	}
	
	// Extract value (tentative)
	
	int tmp_bit;
	for (i = 0; i < nxny; i++) {
		tmp_bit = growmap[i + (nerode + ngrow) * (nx + 1 + 2 * (i / nx + nerode + ngrow))];
		if (tmp_bit && !bitmap[i]) { bitmap[i] = 1; }
	}

	free(erodemap); free(growmap);
	
	//
	
	return;
	
}

// ======================================================================

// Gather and write out arrays

int writeData(DRMS_Record_t *outRec, DRMS_Record_t *inRec,
			  float *Bx, float *probBa, char *ambig_r,
			  int *ll, int *ur, int *ll0, int *ur0, int harpflag)
{
	
	int status = 0;
	
	// Dimension
	
	int i0 = ll0[0] - ll[0];		// padding offset
	int j0 = ll0[1] - ll[1];
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	int nxny = nx * ny;
	int nx0 = ur0[0] - ll0[0] + 1;
	int ny0 = ur0[1] - ll0[1] + 1;
	int nxny0 = nx0 * ny0;
	
	int outDims[2] = {nx0, ny0};
	
	// Flagging pixels
	
	char *confidence = (char *) calloc(nxny0, sizeof(char));
	char *ambig_flag = (char *) calloc(nxny0, sizeof(char));
	
	int idx, idx1;
	for (int i = 0; i < nx0; i++) {
		for (int j = 0; j < ny0; j++) {
			idx = i + j * nx0;
			idx1 = (i + i0) + (j + j0) * nx;
			// char from 0 to 100, but 100 is good.
			confidence[idx] = (char) (probBa[idx1] * 100);
			ambig_flag[idx] = (Bx[idx1] > 0.) ? 1 : 0;
		}
	}
	
	// For HARP, either 7 (111) or 0 (000)
	// Jul 26 2013, XS
	
	if (harpflag) {
	
	  for (int i = 0; i < nxny0; i++) {
	    ambig_flag[i] = (ambig_flag[i] > 0) ? 7 : 0;
	  }
	
	} else {
	
  	// Add random solution as second bit
	
	  char pot, r;
	  srand((unsigned)time(NULL));	// random number generator
	  for (int i = 0; i < nx0; i++) {
		  for (int j = 0; j < ny0; j++) {
			  idx = i + j * nx0;
			  pot = (ambig_flag[idx] % 2);
			  if (confidence[idx] > 61) {		// for weak field only
			  	ambig_flag[idx] += pot * 2;
			  } else {
				  r = rand() % 2;
				  ambig_flag[idx] += r * 2;
		  	}
		  }
	  }
	
	  // Add radial acute solution as third bit
	
  	for (int i = 0; i < nx0; i++) {
		  for (int j = 0; j < ny0; j++) {
			  idx = i + j * nx0;
			  pot = (ambig_flag[idx] % 2);
			  if (confidence[idx] > 61) {
				  ambig_flag[idx] += pot * 4;
			  } else {
				  ambig_flag[idx] += ambig_r[idx] * 4;
			  }
		  }
	  }
	
	}		// end full disk
	
	// Read in bitmap, zero out off-disk pixels in conf_disambig
	// Jun 11 2013
	// Also zero out ambig_flag for off-disk pixels, Jul 26
	// Sep 23: Full disk has no bitmap, fixed
	
	if (harpflag) {
  	DRMS_Segment_t *inSeg_bitmap = drms_segment_lookup(inRec, "bitmap");
		DRMS_Array_t *inArray_bitmap = drms_segment_read(inSeg_bitmap, DRMS_TYPE_CHAR, &status);
		if (status) {
			printf("Bitmap reading error \n");
			return 1;
		}
		char *bitmap = (char *) inArray_bitmap->data;
  	for (int i = 0; i < nx0; i++) {
			for (int j = 0; j < ny0; j++) {
		  	idx = i + j * nx0;
		  	if (!bitmap[idx]) {
		   		confidence[idx] = 0;
		   		ambig_flag[idx] = 0;		// Jul 26
		  	}
			}
		}
		drms_free_array(inArray_bitmap);
	} else {						// Sep 23
		float x0 = drms_getkey_float(inRec, "CRPIX1", &status) - 1.;
		float y0 = drms_getkey_float(inRec, "CRPIX2", &status) - 1.;
		float rsun_obs = drms_getkey_float(inRec, "RSUN_OBS", &status);
		float dx = drms_getkey_float(inRec, "CDELT1", &status);
		float r = rsun_obs / dx;
		printf("x0=%f, y0=%f, r=%f, nx0=%d, ny0=%d\n", x0, y0, r, nx0, ny0);
		for (int i = 0; i < nx0; i++) {
			for (int j = 0; j < ny0; j++) {
		  	idx = i + j * nx0;
		  	if (hypot((i - x0), (j - y0)) > r) {
		  	  confidence[idx] = 0;
		   		ambig_flag[idx] = 0;
		  	}
			}
		}
	}

	// Copy over quality maps, update if necessary
	
	DRMS_Segment_t *inSeg_qual = drms_segment_lookup(inRec, "qual_map");
	DRMS_Array_t *inArray_qual = drms_segment_readslice(inSeg_qual, DRMS_TYPE_INT, ll0, ur0, &status);
	if (status) {
		printf("Qual map reading error \n");
		return 1;
	}
	int *inData_qual = (int *) inArray_qual->data;
	
	DRMS_Segment_t *inSeg_conf = drms_segment_lookup(inRec, "confid_map");
	DRMS_Array_t *inArray_conf = drms_segment_readslice(inSeg_conf, DRMS_TYPE_CHAR, ll0, ur0, &status);
	if (status) {
		printf("Confid map reading error \n");
		return 1;
	}
	char *inData_conf = (char *) inArray_conf->data;
	
	int *qual_map = (int *) calloc(nxny0, sizeof(int));
	char *confid_map = (char *) calloc(nxny0, sizeof(char));
	
	for (int i = 0; i < nx0; i++) {
		for (int j = 0; j < ny0; j++) {
			idx = i + j * nx0;
			// Jun 22 Xudong: not updated for now, waiting for final scheme
			confid_map[idx] = inData_conf[idx];
			qual_map[idx] = inData_qual[idx];
#ifdef QUALMAP
			qual_map[idx] = qual_map[idx] | getAmbCode(confidence[idx]);
#endif
		}
	}
	
	// Creating output arrays
	
	DRMS_Array_t *outArray_flag = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, ambig_flag, &status);
	if (status) { printf("Error creating flag array\n"); return 1; }
	DRMS_Array_t *outArray_prob = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, confidence, &status);
	if (status) { printf("Error creating prob array\n"); return 1; }
	DRMS_Array_t *outArray_qual = drms_array_create(DRMS_TYPE_INT, 2, outDims, qual_map, &status);
	if (status) { printf("Error creating qual array\n"); return 1; }
	DRMS_Array_t *outArray_conf = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, confid_map, &status);
	if (status) { printf("Error creating conf array\n"); return 1; }
	
	// Output
	
	DRMS_Segment_t *outSeg_flag = drms_segment_lookup(outRec, "disambig");
	DRMS_Segment_t *outSeg_prob = drms_segment_lookup(outRec, "conf_disambig");
	DRMS_Segment_t *outSeg_qual = drms_segment_lookup(outRec, "info_map");
	DRMS_Segment_t *outSeg_conf = drms_segment_lookup(outRec, "confid_map");
	
	for (int i = 0; i < 2; i++) {
		outSeg_flag->axis[i] = outArray_flag->axis[i];	// For variable dimensions
		outSeg_prob->axis[i] = outArray_prob->axis[i];
		outSeg_qual->axis[i] = outArray_qual->axis[i];
		outSeg_conf->axis[i] = outArray_conf->axis[i];
	}
	
	outArray_flag->parent_segment = outSeg_flag;
	outArray_prob->parent_segment = outSeg_prob;
	outArray_qual->parent_segment = outSeg_qual;
	outArray_conf->parent_segment = outSeg_conf;
	
	status = drms_segment_write(outSeg_flag, outArray_flag, 0);
	if (status) { printf("Problem writing flag file\n"); return 1; }
	status = drms_segment_write(outSeg_prob, outArray_prob, 0);
	if (status) { printf("Problem writing prob file\n"); return 1; }
	status = drms_segment_write(outSeg_qual, outArray_qual, 0);
	if (status) { printf("Problem writing qual file\n"); return 1; }
	status = drms_segment_write(outSeg_conf, outArray_conf, 0);
	if (status) { printf("Problem writing conf file\n"); return 1; }
	
	// Clean up
	
	drms_free_array(inArray_qual);
	drms_free_array(inArray_conf);
	
	drms_free_array(outArray_flag);
	drms_free_array(outArray_prob);
	drms_free_array(outArray_qual);
	drms_free_array(outArray_conf);
	
	//
	
	return 0;
	
}


// ======================================================================

// Write out aux arrays, no error checking for now

int writeAuxData(DRMS_Record_t *outRec, DRMS_Record_t *inRec, 
				 float *Bx, float *By, float *Bz, float *probBa, float *dBt, int *bitmap,
				 int *ll, int *ur)
{
	
	int status = 0;
	
	int nx = ur[0] - ll[0] + 1;
	int ny = ur[1] - ll[1] + 1;
	
	int outDims[2] = {nx, ny};		// Output array dimensions
	
	DRMS_Array_t *outArray_Bx = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, Bx, &status);
	DRMS_Array_t *outArray_By = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, By, &status);
	DRMS_Array_t *outArray_Bz = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, Bz, &status);
	DRMS_Array_t *outArray_probBa = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, probBa, &status);
	DRMS_Array_t *outArray_dBt = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, dBt, &status);
	DRMS_Array_t *outArray_bitmap = drms_array_create(DRMS_TYPE_INT, 2, outDims, bitmap, &status);
	
	DRMS_Segment_t *outSeg_Bx = drms_segment_lookup(outRec, "Bx");
	DRMS_Segment_t *outSeg_By = drms_segment_lookup(outRec, "By");
	DRMS_Segment_t *outSeg_Bz = drms_segment_lookup(outRec, "Bz");
	DRMS_Segment_t *outSeg_probBa = drms_segment_lookup(outRec, "probBa");
	DRMS_Segment_t *outSeg_dBt = drms_segment_lookup(outRec, "dBt");
	DRMS_Segment_t *outSeg_bitmap = drms_segment_lookup(outRec, "bitmap");
	for (int i = 0; i < 2; i++) {
		outSeg_Bx->axis[i] = outArray_Bx->axis[i];	// For variable dimensions
		outSeg_By->axis[i] = outArray_By->axis[i];
		outSeg_Bz->axis[i] = outArray_Bz->axis[i];
		outSeg_probBa->axis[i] = outArray_probBa->axis[i];
		outSeg_dBt->axis[i] = outArray_dBt->axis[i];
		outSeg_bitmap->axis[i] = outArray_bitmap->axis[i];
	}
	outArray_Bx->parent_segment = outSeg_Bx;
	outArray_By->parent_segment = outSeg_By;
	outArray_Bz->parent_segment = outSeg_Bz;
	outArray_probBa->parent_segment = outSeg_probBa;
	outArray_dBt->parent_segment = outSeg_dBt;
	outArray_bitmap->parent_segment = outSeg_bitmap;
	
	status = drms_segment_write(outSeg_Bx, outArray_Bx, 0);
	status = drms_segment_write(outSeg_By, outArray_By, 0);
	status = drms_segment_write(outSeg_Bz, outArray_Bz, 0);
	status = drms_segment_write(outSeg_probBa, outArray_probBa, 0);
	status = drms_segment_write(outSeg_dBt, outArray_dBt, 0);
	status = drms_segment_write(outSeg_bitmap, outArray_bitmap, 0);
	
	drms_free_array(outArray_Bx);
	drms_free_array(outArray_By);
	drms_free_array(outArray_Bz);
	drms_free_array(outArray_probBa);
	drms_free_array(outArray_dBt);
	drms_free_array(outArray_bitmap);
	
	return 0;
	
}
