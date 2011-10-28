/*
 *  bmap.c
 *
 *	This module maps internal B vector cutout to a specified coordinate
 *	and vector decomposition, with proper WCS coordinate info
 *
 *	Coordinates include:
 *	(1) CCD cutout (original Helioprojective-Cartesian, gnomonic)
 *	(2) Various Helioprojective coordinates (Lambert, Postel, etc.)
 *
 *	Vector decomposition include:
 *	(1)	local x, y, z (with z being radial, x-y plane tanget to solar surface)
 *	(2) field, inclination, azimuth (original)
 *	(3) LoS, transverse, azimuth (wrt plane of sky)
 *
 *	Example:
 *
 *	bmap "in=su_xudong.Bharp_720s[12][2011.02.15_00:00:00_TAI]" "out=su_xudong.Btest_jul5" -v -a -c "rep=fia"
 *	bmap "in=su_xudong.Bharp_720s[12][2011.02.15_00:00:00_TAI]" "out=su_xudong.Btest_jul5" -v "cols=480" "rows=400" "lat=-20.4" "lon=35.0"
 *
 *	Author:
 *			X. Sun (xudong@sun.stanford.edu)
 *
 *	Adapted from:
 *			R. Bogart (mtrack.c)
 *			Y. Liu, Apr 07 2010 (mappingvector_img2plane.c)
 *			X. Sun, Apr 08 2010 (mappingvector_img2plane_x.c)
 *
 *	Version:
 *			v0.0	Jul 07 2011
 *			v0.1	Sep 22 2011
 *			v0.2  Oct 07 2011
 *      v0.3  Oct 17 2011
 *
 *
 *
 *
 *	Notes:
 *			v0.0
 *			v0.1
 *			Added -x flag for local Cartesian vector transformation
 *			The original "xyz" representation does r, theta, phi decomposition
 *			New -x flag uses patch center as the only direction reference
 *			v0.2
 *			Add spherical flag -s, so y is flipped, pointing south as theta component (use with care)
 *      v0.3
 *      Add error propagation option doerr
 *      Now only support field_err, inclination_err and azimuth_err in cutout mode
 *      or Br_err, Bt_err and Bp_err in mapping mode
 *
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
#include "cartography.c"
#include "fresize.h"
#include "finterpolate.h"
#include "img2helioVector.c"
#include "copy_me_keys.c"
#include "errorprop.c"

#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include <omp.h>

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2    (16777216)

// Nyqvist rate at disk center is 0.03 degree. Oversample above 0.015 degree
#define NYQVIST		(0.015)

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

/* ========================================================================================================== */

/* Data structure for mapping info */

struct mapInfo {
	// fixed
	int cutOut, localOpt, globalOpt, fullDisk;
	int localCart, spheric;
	int latlon, noDisamb, fillNan;
	int autoTrack, autoSize;
	int verbose;
  int doerr;
	int nbin, gauss;
	int mapOpt, repOpt, interpOpt;
	float ref_lon, ref_lat;
	float xscale, yscale;
	int ref_nrow, ref_ncol;
	// variable
	float xc, yc;		// reference point
	int nrow, ncol;		// size
};

/* ========================================================================================================== */

char *mapName[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
	"LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
	"stereographic", "orthographic", "LambertAzimuthal"};
char *wcsCode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
	"SIN", "ZEA"};
char *repName[] = {"Local Bx, By, Bz", "Field, inclination, azimuth", 
	"Blos, Btrans, azimuth"};
char *segName[] = { "bx", "by", "bz",
	"field", "inclination", "azimuth",
	"blon", "btrans", "azimuth"};
char *interpName[] = {"Wiener", "Cubic Convolution", "Bilinear", "Nearest Neighbor"};

/* ========================================================================================================== */

/* Determine reference point coordinate and patch size according to input */

int findPosition(DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Read full disk images, assuming always 4096x4096, field in CCD xyz representation */

int readFullDisk(DRMS_Record_t *inRec, struct mapInfo *mInfo, float *bx_img, float *by_img, float *bz_img);

/* Mapping function */

int performMapping(DRMS_Record_t *inRec, struct mapInfo *mInfo,
				   float *bx_img, float *by_img, float *bz_img,
				   float *bx_map, float *by_map, float *bz_map);

/* Get error estimation */

int getError(DRMS_Record_t *inRec, struct mapInfo *mInfo,
             float *b1_err, float *b2_err, float *b3_err,
             int *qual_map, char *confid_map);

/* Performing local vector transformation */

int localVectorTransform(DRMS_Record_t *inRec, struct mapInfo *mInfo, 
						 float *bx_map, float *by_map, float *bz_map);


/* Convert result to requested vector presentation */

void repVectorTransform(float *bx_map, float *by_map, float *bz_map, 
				  float *b1, float *b2, float *b3,
				  struct mapInfo *mInfo);

/* Write to output series */

int writeMap(DRMS_Record_t *outRec, struct mapInfo *mInfo,
			 float *b1, float *b2, float *b3);

/* Write error estimation to output series */

int writeError(DRMS_Record_t *outRec, struct mapInfo *mInfo,
               float *b1_err, float *b2_err, float *b3_err,
               int *qual_map, char *confid_map);

/* For testing, output Lat/Lon */

void outputLatLon(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Set all keywords */

void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct mapInfo *mInfo);


/* Get relevant ephemeris */

int getEphemeris(DRMS_Record_t *inRec, double *disk_lonc, double *disk_latc, 
				 double *disk_xc, double *disk_yc, double *rSun, double *asd, double *pa);

/* Nearest neighbor interpolation */

float nnb (float *f, int nx, int ny, double x, double y);

/* Wrapper for Jesper's rebin code */

void frebin (float *image_in, float *image_out, int nx, int ny, int nbin, int gauss);


/* ========================================================================================================== */

char *module_name = "bmap";

ModuleArgs_t module_args[] =
{
	{ARG_STRING, "in", "Not Specified", "Input vector image series."},
	{ARG_STRING, "out", "Not Specified", "Output map series."},
	{ARG_NUME, "map", "carree", "Projetion method, carree by default.",
		"carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
	{ARG_NUME, "rep", "xyz", "Vector representation, xyz by default.",
		"xyz, fia, lta"},
	{ARG_NUME, "interp", "wiener", "Interpolation method, higher order polynomial by default.",
		"wiener, cubic, bilinear, nearest"},
	{ARG_FLOAT, "lon", "0.", "Reference longitude for output map center."},
	{ARG_FLOAT, "lat", "0.", "Reference latitude for output map center."},
	{ARG_FLOAT, "xscale", "0.03", "Scale in x direction at map center, in degrees."},
	{ARG_FLOAT, "yscale", "0.03", "Scale in x direction at map center, in degrees."},
	{ARG_INT, "cols", "0", "Columns in output map"},
	{ARG_INT, "rows", "0", "Rows in output map"},
	{ARG_INT, "nbin", "3", "Oversampling rate for anti-aliasing. See documentation."},
  {ARG_INT, "doerr", "1", "Perform error propagation, use near neighbor so far."},
	{ARG_FLAG, "c",	"", "Simple full pixel cutout without interpolation, override all mapping options."},
	{ARG_FLAG, "l",	"", "Force vector representation in local tangent coordinates, for cutout mode only."},
	{ARG_FLAG, "g",	"", "Force vector representation in plane of sky coordinates, for remap mode only."},
	{ARG_FLAG, "a",	"", "Automatic tracking, override lon/lat."},
	{ARG_FLAG, "b", "", "Boxcar smoothing instead of Gaussian."},
	{ARG_FLAG, "f", "", "Full-disk mode, no mapping allowed."},
	{ARG_FLAG, "r", "", "Use original resolution, suppressing oversampling."},
	{ARG_FLAG, "e", "", "Output latitude, longitude."},
	{ARG_FLAG, "z", "", "Skip disambiguation solution."},
	{ARG_FLAG, "v", "", "Verbose mode."},
	{ARG_FLAG, "x", "", "Uses patch center direction for vector decomposition"},
	{ARG_FLAG, "q", "", "Fill regions outside HARP defined area as NaN"},
	{ARG_FLAG, "s", "", "Spherical option, flip By as Btheta"},
	{ARG_END}
};

int DoIt(void)
{
	
	int status = DRMS_SUCCESS;
	int nrecs, irec;
	
	struct mapInfo mInfo;
	
	char *inQuery, *outQuery;
	
	DRMS_RecordSet_t *inRS = NULL, *outRS = NULL;
	DRMS_Record_t *inRec = NULL, *outRec = NULL;
	
	/* Get parameters */
    
  inQuery = (char *) params_get_str(&cmdparams, "in");
  outQuery = (char *) params_get_str(&cmdparams, "out");
	
	mInfo.mapOpt = params_get_int(&cmdparams, "map");
	mInfo.repOpt = params_get_int(&cmdparams, "rep");
	mInfo.interpOpt = params_get_int(&cmdparams, "interp");
	
	mInfo.ref_lon = params_get_float(&cmdparams, "lon");
	mInfo.ref_lat = params_get_float(&cmdparams, "lat");
	
	mInfo.xscale = params_get_float(&cmdparams, "xscale");
	mInfo.yscale = params_get_float(&cmdparams, "yscale");
	mInfo.nbin = params_get_int(&cmdparams, "nbin");
  
  mInfo.doerr = params_get_int(&cmdparams, "doerr");
	
	if ((mInfo.xscale / mInfo.nbin) > NYQVIST || (mInfo.yscale / mInfo.nbin) > NYQVIST)
		mInfo.nbin = MAX((round(mInfo.xscale / NYQVIST)),(round(mInfo.yscale / NYQVIST)));
	if (mInfo.xscale <= NYQVIST && mInfo.yscale <= NYQVIST)
		mInfo.nbin = 1;
	mInfo.nbin = mInfo.nbin / 2 * 2 + 1;
	if (params_isflagset(&cmdparams, "r")) mInfo.nbin = 1;
	
//	printf("nbin=%d\n", mInfo.nbin);
	
	mInfo.ref_ncol = params_get_int(&cmdparams, "cols");
	mInfo.ref_nrow = params_get_int(&cmdparams, "rows");
	if (mInfo.ref_ncol <= 0 || mInfo.ref_nrow <= 0)
		mInfo.autoSize = 1;
	else
		mInfo.autoSize = 0;
	
	mInfo.localCart = params_isflagset(&cmdparams, "x");
	mInfo.cutOut = params_isflagset(&cmdparams, "c");
	mInfo.localOpt = params_isflagset(&cmdparams, "l");
	mInfo.globalOpt = params_isflagset(&cmdparams, "g");
	
	mInfo.autoTrack = params_isflagset(&cmdparams, "a");
	mInfo.gauss = !(params_isflagset(&cmdparams, "b"));
	
	mInfo.fullDisk = params_isflagset(&cmdparams, "f");
	if (mInfo.fullDisk) {
		mInfo.cutOut = 1;
		mInfo.ref_ncol = mInfo.ref_nrow = mInfo.nrow = mInfo.ncol = FOURK;
	}
	
	mInfo.latlon = params_isflagset(&cmdparams, "e");
	mInfo.noDisamb = params_isflagset(&cmdparams, "z");
	mInfo.verbose = params_isflagset(&cmdparams, "v");
	mInfo.fillNan = params_isflagset(&cmdparams, "q");
	mInfo.spheric = params_isflagset(&cmdparams, "s");
	
	/* Check and print parameters */
	
	if (mInfo.verbose) {
		printf("==============\nMapping %s:\n", inQuery);
		if (mInfo.cutOut) {
			printf("Full pixel cutout\n");
			if (mInfo.localOpt) printf("Use local tangent plane as reference\n");
		} else {
			printf("Mapping: %s\n", mapName[mInfo.mapOpt]);
			printf("Interpolation: %s\n", interpName[mInfo.interpOpt]);
			printf("Scale: %f, %f deg\n", mInfo.xscale, mInfo.yscale);
			if (mInfo.globalOpt) printf("Use plane of sky as reference\n");
		}
		printf("Representation: %s\n", repName[mInfo.repOpt]);
		if (mInfo.autoTrack) {
			printf("Automatic centered\n");
		} else {
			printf("Use reference lon/lat: %f, %f\n", mInfo.ref_lon, mInfo.ref_lat);
		}
		if (mInfo.autoSize) {
			printf("Automatic sized\n");
		} else {
			printf("Map size: %d x %d\n", mInfo.ref_ncol, mInfo.ref_nrow);
		}
		if (mInfo.fullDisk) printf("Full disk mode\n");
		if (mInfo.latlon) printf("Output latitude/longitude\n");
		if (mInfo.noDisamb) printf("No disambiguation\n");
	}
	
	/* Input data */
	
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;
	
    /* Output data */

    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output recordset not created");
	
	/* Start */
	
	printf("==============\nStart mapping, %d image(s) in total.\n", nrecs);

//	nrecs = 2;
	for (irec = 0; irec < nrecs; irec++) {
		
		inRec = inRS->records[irec];
        TIME trec = drms_getkey_time(inRec, "T_REC", &status);
		
		if (mInfo.verbose) {
			printf("==============\nImage #%d: ", irec);
			int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
			if (!mInfo.fullDisk && !status) printf("[%d]", harpnum);
			char trec_str[40];
			sprint_time(trec_str, trec, "TAI", 0);
			printf("[%s]\n", trec_str);
		}
		
		// Check and determine all relevant geometry
		
		if (mInfo.verbose) SHOW("Determine geometry...");
		
		if (!mInfo.fullDisk &&
			findPosition(inRec, &mInfo)) {
			printf(" erroneous geometry, image skipped.\n");
			continue;
		}
		
		if (mInfo.verbose) printf(" done.\n");
		
		// Read all relevant segments, convert field to image bxyz
		
		float *bx_img = NULL, *by_img = NULL, *bz_img = NULL;		// in CCD bxyz first
		
		bx_img = (float *) malloc(FOURK * FOURK * sizeof(float));
		by_img = (float *) malloc(FOURK * FOURK * sizeof(float));
		bz_img = (float *) malloc(FOURK * FOURK * sizeof(float));
		
		if (mInfo.verbose) SHOW("Read input data...");
		
		if (readFullDisk(inRec, &mInfo, bx_img, by_img, bz_img)) {
			if (mInfo.verbose) printf(" segment read error, image skipped.\n");
			if (bx_img) free(bx_img); if (by_img) free(by_img); if (bz_img) free(bz_img);
			continue;
		}
		
		if (mInfo.verbose) printf(" done.\n");

		// Start mapping, oversampling/smoothing if neeed
		
		float *bx_map = NULL, *by_map = NULL, *bz_map = NULL;	// intermediate maps, in CCD bxyz representation
		int mapsize = mInfo.ncol * mInfo.nrow;
		
		if (!mInfo.fullDisk) {
			bx_map = (float *) malloc(mapsize * sizeof(float));
			by_map = (float *) malloc(mapsize * sizeof(float));
			bz_map = (float *) malloc(mapsize * sizeof(float));
		} else {
			bx_map = bx_img;
			by_map = by_img;
			bz_map = bz_img;
		}

		if (mInfo.verbose) SHOW("Mapping...");
		
		if (!mInfo.fullDisk &&
			performMapping(inRec, &mInfo,
						   bx_img, by_img, bz_img,
						   bx_map, by_map, bz_map)) {
			printf("Image #%d skipped.", irec);
			if (mInfo.verbose) printf(": mapping error\n"); else printf("\n");
			if (bx_img) free(bx_img); if (by_img) free(by_img); if (bz_img) free(bz_img);
			if (bx_map) free(bx_map); if (by_map) free(by_map); if (bz_map) free(bz_map);
			continue;
		}
		
		if (mInfo.verbose) printf(" done.\n");
      
    // Added Oct 17 2011: Do error
    
    float *b1_err = NULL, *b2_err = NULL, *b3_err = NULL;   // final error arrays
    int *qual_map = NULL; char *confid_map = NULL;
    
    if (mInfo.doerr) {
      if (mInfo.verbose) SHOW("Getting error estimation...");
      
      if (!mInfo.fullDisk) {
        b1_err = (float *) malloc(mapsize * sizeof(float));
        b2_err = (float *) malloc(mapsize * sizeof(float));
        b3_err = (float *) malloc(mapsize * sizeof(float));
        qual_map = (int *) malloc(mapsize * sizeof(int));
        confid_map = (char *) malloc(mapsize * sizeof(char));
      } else {
        b1_err = (float *) malloc(FOURK * FOURK * sizeof(float));
        b2_err = (float *) malloc(FOURK * FOURK * sizeof(float));
        b3_err = (float *) malloc(FOURK * FOURK * sizeof(float));
        qual_map = (int *) malloc(FOURK * FOURK * sizeof(int));
        confid_map = (char *) malloc(FOURK * FOURK * sizeof(char));
      }

if (1) {     
      if (getError(inRec, &mInfo, b1_err, b2_err, b3_err, qual_map, confid_map)) {
        printf("Image #%d skipped.", irec);
        if (mInfo.verbose) printf(": error estimation error\n"); else printf("\n");
        continue;
      }
}      
      if (mInfo.verbose) printf(" done.\n");
    }
		
		// Representation
		// Note that the definition of representation depends on the coordinate
		
		// In cutout format (wrt plane of sky (CCD)):
		// xyz: z refers to LoS, outward positive; x to the right; y up
		// fia: inclination within [0,180], 0 pointing +z; azimuth 0 pointing +y, CCW.
		// lta: longitudinal along +z, transverse in plane of sky
		
		// In remap format (local):
		// xyz: z refers to vertical (radial) component, x largely to west (phi), y north (-theta)
		// if -x is specified, change to local xyz
		// if -s is specified, change to spherical setup, y is flipped
		// fia: inclination 0 pointing +z, azimuth 0 pointing +y, CCW
		// lta: vertical/horizontal/azimuth

		if (mInfo.verbose) SHOW("Vector transformation...");
		
		if (localVectorTransform(inRec, &mInfo, bx_map, by_map, bz_map)) {
			printf("Image #%d skipped.", irec);
			if (mInfo.verbose) printf(": vector transformation error\n"); else printf("\n");
			if (bx_img) free(bx_img); if (by_img) free(by_img); if (bz_img) free(bz_img);
			if (bx_map) free(bx_map); if (by_map) free(by_map); if (bz_map) free(bz_map);
			continue;
		}

		// Final output arrays, convert to output format
		
		float *b1 = NULL, *b2 = NULL, *b3 = NULL;
		
		b1 = (float *) malloc(mapsize * sizeof(float));
		b2 = (float *) malloc(mapsize * sizeof(float));
		b3 = (float *) malloc(mapsize * sizeof(float));
		
		repVectorTransform(bx_map, by_map, bz_map, b1, b2, b3, &mInfo);
		
		/*
		 for (int row = 0; row < mInfo.nrow; row++) {
		 for (int col = 0; col < mInfo.ncol; col++) {
		 int ind_map = row * mInfo.ncol + col;
		 printf("%f, %f, %f\n", b1[ind_map], b2[ind_map], b3[ind_map]);				
		 }
		 }
		 */
		
		if (mInfo.verbose) printf(" done.\n");
		
		// Clean up
		
		if (bx_img) free(bx_img); if (by_img) free(by_img); if (bz_img) free(bz_img);
		if (!mInfo.fullDisk) {
			if (bx_map) free(bx_map); if (by_map) free(by_map); if (bz_map) free(bz_map);
		}
		
		// Output
		// note that b1, b2, b3 are freed in writeMap 
		
    outRec = outRS->records[irec];  
		
		if (mInfo.verbose) SHOW("Output...");
		
		if (writeMap(outRec, &mInfo, b1, b2, b3)) {
			printf("Image #%d skipped.", irec);
			if (mInfo.verbose) printf(": output error\n"); else printf("\n");
			if (b1) free(b1); if (b2) free(b2); if (b3) free(b3);
			continue;
		}
		
		if (mInfo.verbose) printf(" done.\n");
    
    // Write error if necessary
    
    if (mInfo.doerr) {
      if (mInfo.verbose) SHOW("Output error estimation...");
      if (writeError(outRec, &mInfo, b1_err, b2_err, b3_err, qual_map, confid_map)) {
        printf("Image #%d skipped.", irec);
        if (mInfo.verbose) printf(": error estimation output error\n"); else printf("\n");
        if (b1_err) free(b1_err); if (b2_err) free(b2_err); if (b3_err) free(b3_err);
        continue;
      }
      if (mInfo.verbose) printf(" done.\n");     
    }
		
		// For testing, output lat/lon
		
		if (mInfo.latlon) {
			if (mInfo.verbose) SHOW("Output lat/lon...");
			outputLatLon(outRec, inRec, &mInfo);
			if (mInfo.verbose) printf(" done.\n");
		}
		
		// Keywords
		
		drms_copykey(outRec, inRec, "T_REC");
    drms_copykey(outRec, inRec, "HARPNUM");
		drms_copykey(outRec, inRec, "RUNNUM");		// for Monte-Carlo
		
		setKeys(outRec, inRec, &mInfo);
		
		//
		
		printf("Image #%d done.\n", irec);
		
	}	// irec
	
	drms_close_records(inRS, DRMS_FREE_RECORD);
	
  drms_close_records(outRS, DRMS_INSERT_RECORD);
	
	return 0;
	
}	// end DoIt



/* ======================================= */


/* 
 * Determine reference point coordinate and patch size according to input
 * xc, yc are the coordinate of patch center, in pixel (full or pm 0.5) if cutout, in degrees if mapping
 * ncol and nrow are the final size
 *
 */

int findPosition(DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	int status = 0;
	
	// Get all ephemeris
	
	double disk_latc, disk_lonc, disk_xc, disk_yc, rSun, asd, pa;
	
	if (getEphemeris(inRec, &disk_lonc, &disk_latc, &disk_xc, &disk_yc, &rSun, &asd, &pa))
		return 1;
	
	/*
	printf("disk_latc=%f, disk_lonc=%f, disk_xc=%f, disk_yc=%f, rSun=%f, asd=%f, pa=%f,\n", 
		   disk_latc/RADSINDEG, disk_lonc/RADSINDEG, disk_xc, disk_yc, rSun, asd, pa/RADSINDEG);*/
	
	/* ==== Look for center ==== */
	
	double lonc, latc;
	
	float minlon, maxlon, minlat, maxlat;
	minlon = drms_getkey_float(inRec, "MINLON", &status); if (status) return 1;		// Stonyhurst lon
	maxlon = drms_getkey_float(inRec, "MAXLON", &status); if (status) return 1;
	minlat = drms_getkey_float(inRec, "MINLAT", &status); if (status) return 1;
	maxlat = drms_getkey_float(inRec, "MAXLAT", &status); if (status) return 1;
	
	if (mInfo->autoTrack) {
		
		// use HARP defined lat/lon center
		
		latc = (maxlat + minlat) / 2.;
		lonc = (maxlon + minlon) / 2. + disk_lonc / RADSINDEG;
		
	} else {
		
		// use specified
		
		latc = mInfo->ref_lat;
		lonc = mInfo->ref_lon;
		
		// It's possible that reference pixel and disk center belong to two Carrington rotations
		// We make sure the difference is less than 180 deg, i.e. always using the current image's CAR_ROT
		
		if ((lonc * RADSINDEG - disk_lonc) >= PI) lonc -= 360.;		// ref point to the east of disk center
		if ((lonc * RADSINDEG - disk_lonc) <= -PI) lonc += 360.;	// ref point to the west of disk center
		
	}
	
	double xi, zeta;		// reference point in CCD pixel
	
	if (mInfo->cutOut) {
		
		if (sphere2img (latc * RADSINDEG, lonc * RADSINDEG, disk_latc, disk_lonc, &xi, &zeta, 
						 disk_xc/rSun, disk_yc/rSun, 1.0, pa, 1., 0., 0., 0.))
			return 1;			// off disk
		
		xi *= rSun;
		zeta *= rSun;
		
		if (!mInfo->autoSize) {
			if (mInfo->ncol % 2 == 0) mInfo->xc = (int) xi + 0.5; else mInfo->xc = round(xi);
			if (mInfo->nrow % 2 == 0) mInfo->yc = (int) zeta + 0.5; else mInfo->yc = round(zeta);
		} else {
			mInfo->xc = round(xi);
			mInfo->yc = round(zeta);
		}
		
	} else {
		
		mInfo->xc = lonc;
		mInfo->yc = latc;
		
	}
	
	/* ==== Determine size ==== */
	
	if (mInfo->autoSize) {		// based on HARP size
		
		if (mInfo->cutOut) {
			
			DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, "disambig");
			mInfo->ncol = inSeg->axis[0] / 2 * 2 + 1;
			mInfo->nrow = inSeg->axis[1] / 2 * 2 + 1;		// full pixel center requires odd dim
			
		} else {
			
			mInfo->ncol = round((maxlon - minlon) / mInfo->xscale);
			mInfo->nrow = round((maxlat - minlat) / mInfo->yscale);
			
		}
		
	} else {			// whatever that is requested
		
		mInfo->ncol = mInfo->ref_ncol;
		mInfo->nrow = mInfo->ref_nrow;
		
	}
	
//	printf("\n%f, %f, %d, %d\n", mInfo->xc, mInfo->yc, mInfo->ncol, mInfo->nrow);
	
	return 0;
	
}


/* ======================================= */


/* 
 * Read full disk images, assuming always 4096x4096
 * Field in CCD xyz representation, with bz as line-of-sight, bx along CCDx, by along CCDy
 *
 */

int readFullDisk(DRMS_Record_t *inRec, struct mapInfo *mInfo, float *bx_img, float *by_img, float *bz_img)
{
	
	int status = 0;
	
	DRMS_Segment_t *inSeg;
	DRMS_Array_t *inArray_ambig;
    DRMS_Array_t *inArray_bTotal, *inArray_bAzim, *inArray_bIncl, *inArray_bFill;
	
	char *ambig;
	float *bTotal, *bAzim, *bIncl, *bFill;
	
	/* Read in segments */
	
	inSeg = drms_segment_lookup(inRec, "disambig");
	inArray_ambig = drms_segment_read(inSeg, DRMS_TYPE_CHAR, &status);
	if (status) return 1;
	ambig = (char *)inArray_ambig->data;
	
	inSeg = drms_segment_lookup(inRec, "field");
	inArray_bTotal = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bTotal = (float *)inArray_bTotal->data;
	
	inSeg = drms_segment_lookup(inRec, "azimuth");
	inArray_bAzim = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bAzim = (float *)inArray_bAzim->data;
	
	inSeg = drms_segment_lookup(inRec, "inclination");
	inArray_bIncl = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bIncl = (float *)inArray_bIncl->data;
	
	inSeg = drms_segment_lookup(inRec, "alpha_mag");
	inArray_bFill = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bFill = (float *)inArray_bFill->data;

	/* Convert CCD xyz */
	
	int llx, lly;		// lower-left corner
	int bmx, bmy;		// bitmap size
	
	int fullDisk = 0;		// if bitmap is full disk
	if (inArray_ambig->axis[0] == FOURK && inArray_ambig->axis[1] == FOURK)
		fullDisk = 1;
	
	llx = (int)(drms_getkey_float(inRec, "CRPIX1", &status)) - 1;
	lly = (int)(drms_getkey_float(inRec, "CRPIX2", &status)) - 1;
	
	bmx = inArray_ambig->axis[0];
	bmy = inArray_ambig->axis[1];
	
	int kx, ky, kOff;
	int jy = 0, yOff = 0, iData = 0;
	int ix = 0;
	int xDim = FOURK, yDim = FOURK;
	
	for (jy = 0; jy < yDim; jy++)
	{
		ix = 0;
		yOff = jy * xDim;
		ky = jy - lly;
		for (ix = 0; ix < xDim; ix++)
		{
			iData = yOff + ix;
			kx = ix - llx;
			
			/* zero azi pointing up, zero incl pointing out from sun */
			bx_img[iData] = - bTotal[iData] * bFill[iData] * sin(bIncl[iData] * RADSINDEG) * sin(bAzim[iData] * RADSINDEG);
			by_img[iData] = bTotal[iData] * bFill[iData] * sin(bIncl[iData] * RADSINDEG) * cos(bAzim[iData] * RADSINDEG);
			bz_img[iData] = bTotal[iData] * bFill[iData] * cos(bIncl[iData] * RADSINDEG);
            
			/* Disambiguation */

			if (mInfo->noDisamb) continue;
			
			if (mInfo->fullDisk || fullDisk) {
				if (ambig[iData]) {
					bx_img[iData] *= -1.; by_img[iData] *= -1.;
				}
				continue;
			}
			
			if (kx < 0 || kx >= bmx || ky < 0 || ky >= bmy) {
				continue;
			} else {
				kOff = ky * bmx + kx;
				if (ambig[kOff]) {		// 180
					bx_img[iData] *= -1.; by_img[iData] *= -1.;
				} 
			}

		}
	}
	
	/* Clean up */
	
	drms_free_array(inArray_ambig);
	drms_free_array(inArray_bTotal); drms_free_array(inArray_bAzim);
	drms_free_array(inArray_bIncl); drms_free_array(inArray_bFill);
	
	return 0;
	
}


/* ======================================= */


/* 
 * Mapping function
 * A simple cutout if requested
 * Otherwise oversampling by nbin, then smoothing using a Gaussian
 * Requested coordinates in final map is transformed back to CCD coordinate,
 * where interpolation is performed on CCD xyz components
 *
 */

int performMapping(DRMS_Record_t *inRec, struct mapInfo *mInfo,
				   float *bx_img, float *by_img, float *bz_img,
				   float *bx_map, float *by_map, float *bz_map)
{
	
	int status = 0;
	
	if (mInfo->cutOut) {		// for cutout, simply copy over
		
		int nrow = mInfo->nrow, ncol = mInfo->ncol;
		float xc = mInfo->xc, yc = mInfo->yc;
		
		// printf("ncol=%d, nrow=%d, xc=%f, yc=%f\n", ncol, nrow, xc, yc);
		
		int ind_map, ind_img;
		int x, y;
		
		for (int row = 0; row < nrow; row++) {
			for (int col = 0; col < ncol; col++) {
				
				ind_map = row * ncol + col;
				x = round(col + xc + 0.5 - ncol / 2.0);
				y = round(row + yc + 0.5 - nrow / 2.0);
				
				if (x < 0 || x >= FOURK || y < 0 || y >= FOURK) {
					bx_map[ind_map] = DRMS_MISSING_FLOAT;
					by_map[ind_map] = DRMS_MISSING_FLOAT;
					bz_map[ind_map] = DRMS_MISSING_FLOAT;
				} else {
					ind_img = y * FOURK + x;
					bx_map[ind_map] = bx_img[ind_img];
					by_map[ind_map] = by_img[ind_img];
					bz_map[ind_map] = bz_img[ind_img];
				//	printf("%f, %f, %f\n", bx_img[ind_img], by_img[ind_img], bz_img[ind_img]);
				}
				
			}
		}
		
		// perform NaN filling if necessary
		if (mInfo->fillNan) {
		  int col0 = drms_getkey_float(inRec, "CRPIX1", &status) - 1.;
		  int row0 = drms_getkey_float(inRec, "CRPIX2", &status) - 1.;
		  DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, "disambig");
		  int col1 = col0 + inSeg->axis[0];
		  int row1 = row0 + inSeg->axis[1];
  		for (int row = 0; row < nrow; row++) {
  			for (int col = 0; col < ncol; col++) {
    			ind_map = row * ncol + col;
  		  	x = round(col + xc + 0.5 - ncol / 2.0);
			  	y = round(row + yc + 0.5 - nrow / 2.0);
			  	if (x < col0 || x >= col1 || y < row0 || y >= row1) {
					  bx_map[ind_map] = DRMS_MISSING_FLOAT;
					  by_map[ind_map] = DRMS_MISSING_FLOAT;
					  bz_map[ind_map] = DRMS_MISSING_FLOAT;
				  }
  			}
  		}
		
		}
		
	} else {				// mapping
		
		// Now get all ephemeris
		
		double disk_latc, disk_lonc, disk_xc, disk_yc, rSun, asd, pa;
		
		if (getEphemeris(inRec, &disk_lonc, &disk_latc, &disk_xc, &disk_yc, &rSun, &asd, &pa))
			return 1;
		
		// Oversampling output
		
		int ncol0, nrow0;		// oversampled map size, same as ncol/nrow if cutout
		float xscale0, yscale0;	// oversampling resolution
		
		ncol0 = mInfo->ncol * mInfo->nbin + (mInfo->nbin / 2) * 2;	// pad with nbin/2 on edge to avoid NAN
		nrow0 = mInfo->nrow * mInfo->nbin + (mInfo->nbin / 2) * 2;
		xscale0 = mInfo->xscale / mInfo->nbin;
		yscale0 = mInfo->yscale / mInfo->nbin;
		
//		printf("ncol0=%d, nrow0=%d, xscale0=%f, yscale0=%f\n", ncol0, nrow0, xscale0, yscale0);
		
		float *bx_map0 = NULL, *by_map0 = NULL, *bz_map0 = NULL;	// intermediate maps, in CCD xyz representation
		
		bx_map0 = (float *) malloc(ncol0 * nrow0 * sizeof(float));
    by_map0 = (float *) malloc(ncol0 * nrow0 * sizeof(float));
    bz_map0 = (float *) malloc(ncol0 * nrow0 * sizeof(float));
		
		double lonc = mInfo->xc, latc = mInfo->yc;
		
		float *xi_out, *zeta_out;		// coordinates of sampled point
		
		xi_out = (float *) malloc(ncol0 * nrow0 * sizeof(float));
		zeta_out = (float *) malloc(ncol0 * nrow0 * sizeof(float));
		
		double x, y;			// map coord
		double lat, lon;		// helio coord
		double xi, zeta;		// image coord (for one point)
		
		int ind_map;
		
//		printf("xc=%f, yc=%f\n", lonc, latc);
		
		// Determine sampling coord
		
		
		for (int row0 = 0; row0 < nrow0; row0++) {
			for (int col0 = 0; col0 < ncol0; col0++) {
				
				ind_map = row0 * ncol0 + col0;
				
				x = (col0 + 0.5 - ncol0/2.) * xscale0;
        y = (row0 + 0.5 - nrow0/2.) * yscale0;
				
				/* map grid [x, y] corresponds to the point [lon, lat] in the heliographic coordinates. 
         * the [x, y] are in radians with respect of the center of the map [xcMap, ycMap].
         * projection methods could be Mercator, Lambert, and many others. [maplonc, mapLatc]
         * is the heliographic longitude and latitude of the map center. Both are in degree.    
         */
				
        if (plane2sphere (x * RADSINDEG, y * RADSINDEG, latc * RADSINDEG, lonc * RADSINDEG, 
                          &lat, &lon, mInfo->mapOpt)) {
          xi_out[ind_map] = -1;
          zeta_out[ind_map] = -1;
          continue;
        }
				
        /* map the grid [lon, lat] in the heliographic coordinates to [xi, zeta], a point in the
         * image coordinates. The image properties, xCenter, yCenter, rSun, pa, ecc and chi are given.
         */
				
        if (sphere2img (lat, lon, disk_latc, disk_lonc, &xi, &zeta, 
                        disk_xc/rSun, disk_yc/rSun, 1.0, pa, 1., 0., 0., 0.)) {
          xi_out[ind_map] = -1;
          zeta_out[ind_map] = -1;
          continue;
        }
				
        xi_out[ind_map] = xi * rSun;
        zeta_out[ind_map] = zeta * rSun;
				
			}
		}
		
		/* Interpolation */
		
		struct fint_struct pars;
		
		switch (mInfo->interpOpt) {
			default:
			case 0:			// Wiener, 6 order, 1 constraint
				init_finterpolate_wiener(&pars, 6, 1, 6, 2, 1, 1, NULL);
				break;
			case 1:			// Cubic convolution
				init_finterpolate_cubic_conv(&pars, 1., 3.);
				break;
			case 2:			// Bilinear
				init_finterpolate_linear(&pars, 1.);
				break;
			case 3:			// Nearest neighbor
				break;
		}
		
		if (mInfo->interpOpt == 3) {
			for (int row0 = 0; row0 < nrow0; row0++) {
				for (int col0 = 0; col0 < ncol0; col0++) {
					ind_map = row0 * ncol0 + col0;
					bx_map0[ind_map] = nnb(bx_img, FOURK, FOURK, xi_out[ind_map], zeta_out[ind_map]);
					by_map0[ind_map] = nnb(by_img, FOURK, FOURK, xi_out[ind_map], zeta_out[ind_map]);
					bz_map0[ind_map] = nnb(bz_img, FOURK, FOURK, xi_out[ind_map], zeta_out[ind_map]);
				}
			}
		} else {
			finterpolate(&pars, bx_img, xi_out, zeta_out, bx_map0, 
						 FOURK, FOURK, FOURK, ncol0, nrow0, ncol0, DRMS_MISSING_FLOAT);
			finterpolate(&pars, by_img, xi_out, zeta_out, by_map0, 
						 FOURK, FOURK, FOURK, ncol0, nrow0, ncol0, DRMS_MISSING_FLOAT);
			finterpolate(&pars, bz_img, xi_out, zeta_out, bz_map0, 
						 FOURK, FOURK, FOURK, ncol0, nrow0, ncol0, DRMS_MISSING_FLOAT);
		}
		
		/* Rebinning/smoothing */
		
		frebin(bx_map0, bx_map, ncol0, nrow0, mInfo->nbin, mInfo->gauss);
		frebin(by_map0, by_map, ncol0, nrow0, mInfo->nbin, mInfo->gauss);
		frebin(bz_map0, bz_map, ncol0, nrow0, mInfo->nbin, mInfo->gauss);
		
		/* Clean up */
		
		free(bx_map0); free(by_map0); free(bz_map0);
		free(xi_out); free(zeta_out);
		
	}
	
	return 0;

}

/* ======================================= */

/* Get error estimation */

int getError(DRMS_Record_t *inRec, struct mapInfo *mInfo,
             float *b1_err, float *b2_err, float *b3_err,
             int *qual_map, char *confid_map)
{
  int status = 0;
  
  int prop = (mInfo->cutOut && mInfo->localOpt == 1) || 
              (!mInfo->cutOut && mInfo->globalOpt == 0);
  
  DRMS_Segment_t *inSeg;
  
  DRMS_Array_t *inArray_bTotal, *inArray_bAzim, *inArray_bIncl;
  DRMS_Array_t *inArray_qual, *inArray_confid;
  DRMS_Array_t *inArray_errbT, *inArray_errbAz, *inArray_errbIn;
  DRMS_Array_t *inArray_bTbI, *inArray_bTbA, *inArray_bAbI;
  
  float *bTotal, *bAzim, *bIncl;
  int *qual; char *confid;
  float *errbT0, *errbAz0, *errbIn0;
  float *errbT, *errbAz, *errbIn;
  float *errbTbI0, *errbTbA0, *errbAbI0;
  float *errbTbI, *errbTbA, *errbAbI;
  
  // Read full disk images
  
	inSeg = drms_segment_lookup(inRec, "field");
	inArray_bTotal = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bTotal = (float *)inArray_bTotal->data;
	
	inSeg = drms_segment_lookup(inRec, "azimuth");
	inArray_bAzim = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bAzim = (float *)inArray_bAzim->data;
	
	inSeg = drms_segment_lookup(inRec, "inclination");
	inArray_bIncl = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (status) return 1;
	bIncl = (float *)inArray_bIncl->data;
	  
  // Bitmaps
  
  inSeg = drms_segment_lookup(inRec, "qual_map");
  inArray_qual = drms_segment_read(inSeg, DRMS_TYPE_INT, &status);
  if (status) return 1;
  qual = (int *)inArray_qual->data;
  
  inSeg = drms_segment_lookup(inRec, "confid_map");
  inArray_confid = drms_segment_read(inSeg, DRMS_TYPE_CHAR, &status);
  if (status) return 1;
  confid = (char *)inArray_confid->data;
  
  // Errors, need to convert to variances
  
  inSeg = drms_segment_lookup(inRec, "field_err");
  inArray_errbT = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
  if (status) return 1;
  errbT0 = (float *)inArray_errbT->data;
  
  inSeg = drms_segment_lookup(inRec, "azimuth_err");
  inArray_errbAz = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
  if (status) return 1;
  errbAz0 = (float *)inArray_errbAz->data;
  
  inSeg = drms_segment_lookup(inRec, "inclination_err");
  inArray_errbIn = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
  if (status) return 1;
  errbIn0 = (float *)inArray_errbIn->data;
  
  errbT = (float *) (malloc(FOURK2 * sizeof(float)));
  errbAz = (float *) (malloc(FOURK2 * sizeof(float)));
  errbIn = (float *) (malloc(FOURK2 * sizeof(float)));
  
  for (int i = 0; i < FOURK2; i++) {
    errbT[i] = errbT0[i] * errbT0[i];
    if (fabs(errbAz0[i]) > 180.) errbAz0[i] = 180.;
    if (fabs(errbIn0[i]) > 180.) errbIn0[i] = 180.;
    errbAz[i] = errbAz0[i] * errbAz0[i] * RADSINDEG * RADSINDEG;
    errbIn[i] = errbIn0[i] * errbIn0[i] * RADSINDEG * RADSINDEG;
  }
  
  // Correlation coefficients, need to convert to covariances
  
  if (prop) {
    
    inSeg = drms_segment_lookup(inRec, "field_inclination_err");
    inArray_bTbI = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    errbTbI0 = (float *)inArray_bTbI->data;
    
    inSeg = drms_segment_lookup(inRec, "field_az_err");
    inArray_bTbA = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    errbTbA0 = (float *)inArray_bTbA->data;
    
    inSeg = drms_segment_lookup(inRec, "inclin_azimuth_err");
    inArray_bAbI = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status) return 1;
    errbAbI0 = (float *)inArray_bAbI->data;
    
    errbTbI = (float *) (malloc(FOURK2 * sizeof(float)));
    errbTbA = (float *) (malloc(FOURK2 * sizeof(float)));
    errbAbI = (float *) (malloc(FOURK2 * sizeof(float)));
    
    for (int i = 0; i < FOURK2; i++) {
      errbTbI[i] = errbTbI0[i] * errbT0[i] * errbIn0[i] * RADSINDEG;
      errbTbA[i] = errbTbA0[i] * errbT0[i] * errbAz0[i] * RADSINDEG;
      errbAbI[i] = errbAbI0[i] * errbAz0[i] * errbIn0[i] * RADSINDEG * RADSINDEG;
    }
  }
  
  // Ephemeris
	
	double disk_latc, disk_lonc, disk_xc, disk_yc, rSun, asd, pa;
	
	if (getEphemeris(inRec, &disk_lonc, &disk_latc, &disk_xc, &disk_yc, &rSun, &asd, &pa))
		return 1;
  
  // Get error arrays
  
  double x, y;          // map coord
  double lat, lon;      // spherical coord
  double xi, zeta;      // image coord, round to full pixel
  
  int nrow = mInfo->nrow, ncol = mInfo->ncol;
  float xc = mInfo->xc, yc = mInfo->yc;
  float xscale = mInfo->xscale, yscale = mInfo->yscale;
  double lonc = mInfo->xc, latc = mInfo->yc;
  
  int ind_map, ind_img;

  int dx0, dy0, ind_img0;		// address for cutout maps
  inSeg = drms_segment_lookup(inRec, "disambig");
  int col0 = drms_getkey_float(inRec, "CRPIX1", &status) - 1.;		// offset for qual_map and confid_map
  int row0 = drms_getkey_float(inRec, "CRPIX2", &status) - 1.;
  int nx0 = inSeg->axis[0];
  int ny0 = inSeg->axis[1];
  int col1 = col0 + nx0;
  int row1 = row0 + ny0;
  
  double btSigma2, bpSigma2, brSigma2;
  
  if (mInfo->verbose) {
    if (prop) SHOW("\nError propagation...") else SHOW("\nCopy variance...");
  }

  for (int row = 0; row < nrow; row++) {
    for (int col = 0; col < ncol; col++) {
      
      ind_map = row * ncol + col;
      
      // First, get image coord
      
      if (mInfo->cutOut) {

				xi = round(col + xc + 0.5 - ncol / 2.0);
				zeta = round(row + yc + 0.5 - nrow / 2.0);
        
      } else {
        
        x = (col + 0.5 - ncol/2.) * xscale;
        y = (row + 0.5 - nrow/2.) * yscale;
        
        if (plane2sphere (x * RADSINDEG, y * RADSINDEG, latc * RADSINDEG, lonc * RADSINDEG, 
                          &lat, &lon, mInfo->mapOpt)) {
          qual_map[ind_map] = DRMS_MISSING_CHAR;
          confid_map[ind_map] = DRMS_MISSING_CHAR;
          b1_err[ind_map] = DRMS_MISSING_FLOAT;
          b2_err[ind_map] = DRMS_MISSING_FLOAT;
          b3_err[ind_map] = DRMS_MISSING_FLOAT;
          continue;
        }

        if (sphere2img (lat, lon, disk_latc, disk_lonc, &xi, &zeta, 
                        disk_xc/rSun, disk_yc/rSun, 1.0, pa, 1., 0., 0., 0.)) {
          qual_map[ind_map] = DRMS_MISSING_CHAR;
          confid_map[ind_map] = DRMS_MISSING_CHAR;
          b1_err[ind_map] = DRMS_MISSING_FLOAT;
          b2_err[ind_map] = DRMS_MISSING_FLOAT;
          b3_err[ind_map] = DRMS_MISSING_FLOAT;
          continue;
        }
				
        xi *= rSun; xi = round(xi);
        zeta *= rSun; zeta = round(zeta);     // nearest neighbor
        
      }

      // Now let's get confid_map and qual_map

      dx0 = xi - col0; dy0 = zeta - row0;
      if (dx0 < 0 || dy0 < 0 || dx0 >= nx0 || dy0 >= ny0) {
        qual_map[ind_map] = DRMS_MISSING_INT;
        confid_map[ind_map] = DRMS_MISSING_CHAR;
      } else {
        ind_img0 = dy0 * nx0 + dx0;
        qual_map[ind_map] = qual[ind_img0];
        confid_map[ind_map] = confid[ind_img0];
      }
            
      // Error, note we need to use variances
      // Currently only support fia and local xyz (brt) output
      
      ind_img = round(zeta * FOURK + xi);
      
//      if (ind_map < 10) printf("ind_map=%d, xi=%f, zeta=%f\n", ind_map, xi, zeta);
      
      if (prop) {
        
        if (errorprop (bTotal, bAzim, bIncl, 
                       errbT, errbAz, errbIn, errbTbA, errbTbI, errbAbI, 
                       lon, lat, disk_lonc, disk_latc, pa, FOURK, FOURK, xi, zeta, 
                       &btSigma2, &bpSigma2, &brSigma2)) {
          b1_err[ind_map] = DRMS_MISSING_FLOAT;
          b2_err[ind_map] = DRMS_MISSING_FLOAT;
          b3_err[ind_map] = DRMS_MISSING_FLOAT;
//          SHOW("dod");
          continue;
        }
//        
        b1_err[ind_map] = sqrt(bpSigma2);
        b2_err[ind_map] = sqrt(btSigma2);
        b3_err[ind_map] = sqrt(brSigma2);

      } else {

        b1_err[ind_map] = sqrt(errbT[ind_img]);
        b2_err[ind_map] = sqrt(errbIn[ind_img])/RADSINDEG;
        b3_err[ind_map] = sqrt(errbAz[ind_img])/RADSINDEG;
        
      }
      
    }
  }
  
  // perform NaN filling if necessary
  if (mInfo->cutOut && !mInfo->localOpt && mInfo->fillNan) {
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        ind_map = row * ncol + col;
        int x1 = round(col + xc + 0.5 - ncol / 2.0);
        int y1 = round(row + yc + 0.5 - nrow / 2.0);
        if (x1 < col0 || x1 >= col1 || y1 < row0 || y1 >= row1) {
          b1_err[ind_map] = DRMS_MISSING_FLOAT;
          b2_err[ind_map] = DRMS_MISSING_FLOAT;
          b3_err[ind_map] = DRMS_MISSING_FLOAT;
        }
      }
    }
  }
  
  // Clean up
  
  free(errbT); free(errbAz); free(errbIn);
  drms_free_array(inArray_bTotal); drms_free_array(inArray_bAzim);
  drms_free_array(inArray_bIncl);
  drms_free_array(inArray_qual); drms_free_array(inArray_confid);
  drms_free_array(inArray_errbT); drms_free_array(inArray_errbAz);
  drms_free_array(inArray_errbIn);
  if (prop) {
    free(errbTbI); free(errbTbA); free(errbAbI);
    drms_free_array(inArray_bTbI); drms_free_array(inArray_bTbA);
    drms_free_array(inArray_bAbI);
  }

  return 0;
  
}



/* ======================================= */


/* Performing local vector transformation */
								 
int localVectorTransform(DRMS_Record_t *inRec, struct mapInfo *mInfo, 
						 float *bx_map, float *by_map, float *bz_map)
{

	if (mInfo->cutOut && (! mInfo->localOpt))
		return 0;

	if ((! mInfo->cutOut) && mInfo->globalOpt)
		return 0;
	
	// Ephemeris
	
	double disk_latc, disk_lonc, disk_xc, disk_yc, rSun, asd, pa;
	
	if (getEphemeris(inRec, &disk_lonc, &disk_latc, &disk_xc, &disk_yc, &rSun, &asd, &pa))
		return 1;
	
	// Find lat/lon for each point
	// Then perform vector transformation
	
	int mapsize = mInfo->ncol * mInfo->nrow;
	int ind_map;
	
	double x, y;		// point corrd
	double lat, lon;	// lat / lon for current point
	double latc = mInfo->yc, lonc = mInfo->xc;	// patch center
	
	double bx_tmp, by_tmp, bz_tmp;
	double rho, sinlat, coslat, sig, mu, chi;
	
//	printf("xc=%f, yc=%f\n", mInfo->xc, mInfo->yc);
	
	// For local Cartesian
	if (mInfo->repOpt == 0 && mInfo->localCart == 1) {
		
		if (mInfo->cutOut) {
			x = (int) mInfo->xc;
			y = (int) mInfo->yc;
			x = (x - disk_xc) / rSun;
			y = (y - disk_yc) / rSun;
			img2sphere(x, y, asd, disk_latc, disk_lonc, pa,
					   &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi);
		} else {
			lat = latc * RADSINDEG;
			lon = lonc * RADSINDEG;
		}
	}
	
	
	for (int row = 0; row < mInfo->nrow; row++) {
		for (int col = 0; col < mInfo->ncol; col++) {
			
			ind_map = row * mInfo->ncol + col;
			
			// Find lat/lon for each point
				
			if (mInfo->repOpt != 0 || mInfo->localCart == 0) {
				
				if (mInfo->cutOut) {
					
					x = (int)(col + mInfo->xc + 0.5 - mInfo->ncol / 2.0);
					y = (int)(row + mInfo->yc + 0.5 - mInfo->nrow / 2.0);
					x = (x - disk_xc) / rSun;
					y = (y - disk_yc) / rSun;
					
					if (img2sphere(x, y, asd, disk_latc, disk_lonc, pa,
								   &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi)) {
						bx_map[ind_map] = DRMS_MISSING_FLOAT;
						by_map[ind_map] = DRMS_MISSING_FLOAT;
						bz_map[ind_map] = DRMS_MISSING_FLOAT;
						continue;
					}
					
				} else {
					
					x = (col + 0.5 - mInfo->ncol / 2.) * mInfo->xscale;
					y = (row + 0.5 - mInfo->nrow / 2.) * mInfo->yscale;
					
					if (plane2sphere (x * RADSINDEG, y * RADSINDEG, latc * RADSINDEG, lonc * RADSINDEG, 
									  &lat, &lon, mInfo->mapOpt)) {
						bx_map[ind_map] = DRMS_MISSING_FLOAT;
						by_map[ind_map] = DRMS_MISSING_FLOAT;
						bz_map[ind_map] = DRMS_MISSING_FLOAT;
						continue;
					}
				}
				
			}
			
			// printf("%f, %f\n", lat/RADSINDEG, lon/RADSINDEG);
			
			// Vector transformation
			
			bx_tmp = by_tmp = bz_tmp = 0;
			
			img2helioVector (bx_map[ind_map], by_map[ind_map], bz_map[ind_map],
							 &bx_tmp, &by_tmp, &bz_tmp,
							 lon, lat, disk_lonc, disk_latc, pa);
			
			bx_map[ind_map] = bx_tmp;
			by_map[ind_map] = by_tmp;
			bz_map[ind_map] = bz_tmp;
			
		}
	}

	return 0;
	
}

/* ======================================= */

/* Convert result to requested vector presentation */

void repVectorTransform(float *bx_map, float *by_map, float *bz_map, 
						float *b1, float *b2, float *b3,
						struct mapInfo *mInfo)
{
	
	int mapsize = mInfo->ncol * mInfo->nrow;
	float k = 1;
	
	switch (mInfo->repOpt) {
		default:
		case 0:				// xyz
		  if ((!mInfo->cutOut && !mInfo->localCart && mInfo->spheric) ||
		      (mInfo->cutOut && mInfo->localOpt && !mInfo->localCart && mInfo->spheric))
		    k = -1;
			for (int i = 0; i < mapsize; i++) {
				b1[i] = bx_map[i];
				b2[i] = by_map[i] * k;		// Btheta pointing south
				b3[i] = bz_map[i];
			}
			break;
		case 1:				// field, inclination, azimuth
			for (int i = 0; i < mapsize; i++) {
				b1[i] = sqrt(bx_map[i] * bx_map[i] + by_map[i] * by_map[i] + bz_map[i] * bz_map[i]);
				b2[i] = acos(bz_map[i] / b1[i]) / RADSINDEG;
				b3[i] = atan2(-bx_map[i], by_map[i]) / RADSINDEG;
				if (b3[i] < 0.) b3[i] += 360.;
			}
			break;
		case 2:				// longitudinal/transverse/azimuth or vertical/horizontal/azimuth
			for (int i = 0; i < mapsize; i++) {
				b1[i] = bz_map[i];
				b2[i] = sqrt(bx_map[i] * bx_map[i] + by_map[i] * by_map[i]);
				b3[i] = atan2(-bx_map[i], by_map[i]) / RADSINDEG;
				if (b3[i] < 0.) b3[i] += 360.;
			}
			break;
	}
	
}


/* ======================================= */

/* Write to output series */

int writeMap(DRMS_Record_t *outRec, struct mapInfo *mInfo,
             float *b1, float *b2, float *b3)
{
	
	int status = 0;
	
	int outDims[2];
	outDims[0] = mInfo->ncol; outDims[1] = mInfo->nrow;
	
	DRMS_Segment_t *outSeg;
	DRMS_Array_t *outArray;
	
	// write b1 as the 1st segment
	outSeg = drms_segment_lookupnum(outRec, 0);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, b1, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
	
	// write b2 as the 2nd segment
	outSeg = drms_segment_lookupnum(outRec, 1);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, b2, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
	
	// write b3 as the 3rd segment
	outSeg = drms_segment_lookupnum(outRec, 2);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, b3, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
	
	return 0;
	
}

/* ======================================= */

/* Write error estimation to output series */

int writeError(DRMS_Record_t *outRec, struct mapInfo *mInfo,
               float *b1_err, float *b2_err, float *b3_err,
               int *qual_map, char *confid_map)
{
  
  int status = 0;
	
	int outDims[2];
	outDims[0] = mInfo->ncol; outDims[1] = mInfo->nrow;
	
	DRMS_Segment_t *outSeg;
	DRMS_Array_t *outArray;
	
	// write b1_err as the 4th segment
	outSeg = drms_segment_lookupnum(outRec, 3);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, b1_err, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
	
	// write b2_err as the 5th segment
	outSeg = drms_segment_lookupnum(outRec, 4);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, b2_err, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
	
	// write b3_err as the 6th segment
	outSeg = drms_segment_lookupnum(outRec, 5);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, b3_err, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
  
  // write qual_map
  if (1) {
	outSeg = drms_segment_lookup(outRec, "qual_map");
	outArray = drms_array_create(DRMS_TYPE_INT, 2, outDims, qual_map, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);

	
	// write confid_map
	outSeg = drms_segment_lookup(outRec, "confid_map");
	outArray = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, confid_map, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 1;
	drms_free_array(outArray);
	} 
  return 0;
  
}


/* ======================================= */

/* For testing, output Lat/Lon */

void outputLatLon(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	
	if (!mInfo->latlon) return;
	
	// Ephemeris
	
	double disk_latc, disk_lonc, disk_xc, disk_yc, rSun, asd, pa;
	
	if (getEphemeris(inRec, &disk_lonc, &disk_latc, &disk_xc, &disk_yc, &rSun, &asd, &pa))
		return;
	
	// Arrays
	
	float *carrLat = NULL, *carrLon = NULL;
	int mapsize = mInfo->nrow * mInfo->ncol;
	
	carrLat = (float *) malloc(mapsize * sizeof(float));
	carrLon = (float *) malloc(mapsize * sizeof(float));
	
	// Below is similar to localVectorTransform
	
	// Find lat/lon for each point
	// Then perform vector transformation
	
	int ind_map;
	
	double x, y;		// point corrd
	double lat, lon;	// lat / lon for current point
	double latc = mInfo->yc, lonc = mInfo->xc;	// patch center
	
	double bx_tmp, by_tmp, bz_tmp;
	double rho, sinlat, coslat, sig, mu, chi;
	
/*	printf("nrow=%d, ncol=%d, disk_xc=%f, disk_yc=%f, rSun=%f\n", 
		   mInfo->nrow, mInfo->ncol, disk_xc, disk_yc, rSun);
*/
	
	//	printf("xc=%f, yc=%f\n", mInfo->xc, mInfo->yc);
	
	for (int row = 0; row < mInfo->nrow; row++) {
		for (int col = 0; col < mInfo->ncol; col++) {
			
			ind_map = row * mInfo->ncol + col;
			
			// Find lat/lon for each point
			
			if (mInfo->cutOut) {
				
				if (!mInfo->fullDisk) {
					x = (int)(col + mInfo->xc + 0.5 - mInfo->ncol / 2.0);
					y = (int)(row + mInfo->yc + 0.5 - mInfo->nrow / 2.0);
				} else {
					x = col;
					y = row;
				}

				x = (x - disk_xc) / rSun;
				y = (y - disk_yc) / rSun;
				
				if (img2sphere(x, y, asd, disk_latc, disk_lonc, pa,
							   &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi)) {
					carrLat[ind_map] = DRMS_MISSING_FLOAT;
					carrLon[ind_map] = DRMS_MISSING_FLOAT;
					continue;
				}
				
			} else {
				
				x = (col + 0.5 - mInfo->ncol / 2.) * mInfo->xscale;
				y = (row + 0.5 - mInfo->nrow / 2.) * mInfo->yscale;
				
				if (plane2sphere (x * RADSINDEG, y * RADSINDEG, latc * RADSINDEG, lonc * RADSINDEG, 
								  &lat, &lon, mInfo->mapOpt)) {
					carrLat[ind_map] = DRMS_MISSING_FLOAT;
					carrLon[ind_map] = DRMS_MISSING_FLOAT;
					continue;
				}
				
			}
			
			carrLat[ind_map] = lat / RADSINDEG;
			carrLon[ind_map] = lon / RADSINDEG;
			
		}
	}
	
	// Output, no error checking
	
	int status = 0;
	
	int outDims[2];
	outDims[0] = mInfo->ncol; outDims[1] = mInfo->nrow;
	
	DRMS_Segment_t *outSeg = NULL;
	DRMS_Array_t *outArray;
	
	// write lat as the first segment
	outSeg = drms_segment_lookup(outRec, "lat");
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, carrLat, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	drms_free_array(outArray);
	
	// write lon as the second segment
	outSeg = drms_segment_lookup(outRec, "lon");
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, carrLon, &status);
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
	outArray->parent_segment = outSeg;
	status = drms_segment_write(outSeg, outArray, 0);
	drms_free_array(outArray);
	
}


/* ======================================= */

/* Set all keywords, no error checking for now */

void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	int status = 0;
	char key[64];
	
	// Info
	drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
	drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
	
	// All other keywords
	copy_me_keys(inRec, outRec);
	copy_patch_keys(inRec, outRec);
	copy_ambig_keys(inRec, outRec);
	
	// Geometric
	
	if (mInfo->fullDisk) {
		copy_geo_keys(inRec, outRec);
		return;
	}
	
	double disk_latc, disk_lonc, disk_xc, disk_yc, rSun, asd, pa;
	
	if (getEphemeris(inRec, &disk_lonc, &disk_latc, &disk_xc, &disk_yc, &rSun, &asd, &pa))
		return;
	
	drms_copykey(outRec, inRec, "CAR_ROT");
	
	switch (mInfo->repOpt) {
		default:
		case 0:
			drms_setkey_string(outRec, "BUNIT_000", "Gauss");
			drms_setkey_string(outRec, "BUNIT_001", "Gauss");
			drms_setkey_string(outRec, "BUNIT_002", "Gauss");
			break;
		case 1:
			drms_setkey_string(outRec, "BUNIT_000", "Gauss");
			drms_setkey_string(outRec, "BUNIT_001", "degree");
			drms_setkey_string(outRec, "BUNIT_002", "degree");
			break;
		case 2:
			drms_setkey_string(outRec, "BUNIT_000", "Gauss");
			drms_setkey_string(outRec, "BUNIT_001", "Gauss");
			drms_setkey_string(outRec, "BUNIT_002", "degree");
			break;
	}
	
	if (mInfo->cutOut) {
		
		drms_setkey_float(outRec, "CRPIX1", disk_xc - mInfo->xc + (mInfo->ncol - 1.) / 2. + 1.);
		drms_setkey_float(outRec, "CRPIX2", disk_yc - mInfo->yc + (mInfo->nrow - 1.) / 2. + 1.);
		
		drms_setkey_float(outRec, "CRVAL1", 0);
		drms_setkey_float(outRec, "CRVAL2", 0);
		
		drms_copykey(outRec, inRec, "CDELT1");
		drms_copykey(outRec, inRec, "CDELT2");
		drms_copykey(outRec, inRec, "CUNIT1");
		drms_copykey(outRec, inRec, "CUNIT2");
		drms_copykey(outRec, inRec, "CTYPE1");
		drms_copykey(outRec, inRec, "CTYPE2");
		drms_copykey(outRec, inRec, "CROTA2");
		
		drms_setkey_string(outRec, "PROJECTION", "Cutout");
		
	} else {
		
		drms_setkey_float(outRec, "CRPIX1", mInfo->ncol/2. + 0.5);
		drms_setkey_float(outRec, "CRPIX2", mInfo->nrow/2. + 0.5);
		
		drms_setkey_float(outRec, "CRVAL1", mInfo->xc);
		drms_setkey_float(outRec, "CRVAL2", mInfo->yc);
		drms_setkey_float(outRec, "CDELT1", mInfo->xscale);
		drms_setkey_float(outRec, "CDELT2", mInfo->yscale);
		drms_setkey_string(outRec, "CUNIT1", "degree");
		drms_setkey_string(outRec, "CUNIT2", "degree");
		
		snprintf (key, 64, "CRLN-%s", wcsCode[mInfo->mapOpt]);
		drms_setkey_string(outRec, "CTYPE1", key);
		snprintf (key, 64, "CRLT-%s", wcsCode[mInfo->mapOpt]);
		drms_setkey_string(outRec, "CTYPE2", key);
		drms_setkey_float(outRec, "CROTA2", 0.0);
		
		drms_setkey_string(outRec, "PROJECTION", mapName[mInfo->mapOpt]);
		
	}

}


/* ======================================= */


/* Get relevant ephemeris */

int getEphemeris(DRMS_Record_t *inRec, double *disk_lonc, double *disk_latc, 
				 double *disk_xc, double *disk_yc, double *rSun, double *asd, double *pa)
{
	
	int status = 0;
	
	float dSun, rSun_ref, cdelt;
	float crvalx, crvaly, crpix1, crpix2, crota2;
	double sina, cosa;
	
	dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
	rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
	if (status) rSun_ref = 6.96e8;
	cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
	*asd = asin(rSun_ref/dSun);
	
	*rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
	
	float testFulldisk = drms_getkey_float(inRec, "IMCRPIX1", &status);		// Full disk does not have this key
	
	if (status) {		// full disk
		crvalx = drms_getkey_float(inRec, "CRVAL1", &status);	// center of solar disc in arcsec
		crvaly = drms_getkey_float(inRec, "CRVAL2", &status);
		crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);	// disk center in ccd, original
		crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
	} else {
		crvalx = drms_getkey_float(inRec, "IMCRVAL1", &status);	// center of solar disc in arcsec
		crvaly = drms_getkey_float(inRec, "IMCRVAL2", &status);
		crpix1 = drms_getkey_float(inRec, "IMCRPIX1", &status);	// disk center in ccd, original
		crpix2 = drms_getkey_float(inRec, "IMCRPIX2", &status);
	}
	
	crota2 = drms_getkey_float(inRec, "CROTA2", &status);	// rotation
	sina = sin(crota2 * RADSINDEG); 
	cosa = cos(crota2 * RADSINDEG);
	
	*disk_xc = PIX_X(0.0,0.0) - 1.0;		// Center of disk, starting at 0
	*disk_yc = PIX_Y(0.0,0.0) - 1.0;
	
	*disk_latc = drms_getkey_float(inRec, "CRLT_OBS", &status) * RADSINDEG;
	*disk_lonc = drms_getkey_float(inRec, "CRLN_OBS", &status) * RADSINDEG;
	*pa = - crota2 * RADSINDEG;
	
	return 0;
	
}

/* ############# Nearest neighbour interpolation ############### */

float nnb (float *f, int nx, int ny, double x, double y)
{
	
	if (x <= -0.5 || y <= -0.5 || x > nx - 0.5 || y > ny - 0.5)
		return DRMS_MISSING_FLOAT;
	int ilow = floor (x);
	int jlow = floor (y);
	int i = ((x - ilow) > 0.5) ? ilow + 1 : ilow;
	int j = ((y - jlow) > 0.5) ? jlow + 1 : jlow;
	return f[j * nx + i];

}

/* ################## Wrapper for Jesper's rebin code ################## */

void frebin (float *image_in, float *image_out, int nx, int ny, int nbin, int gauss)
{
	
	struct fresize_struct fresizes;
	int nxout, nyout, xoff, yoff;
	int nlead = nx;
	
	nxout = nx / nbin; nyout = ny / nbin;
	if (gauss && nbin != 1)
		init_fresize_gaussian(&fresizes, (nbin / 2), (nbin / 2 * 2), nbin);		// for nbin=3, sigma=1, half truncate width=2
	else
		init_fresize_bin(&fresizes, nbin);
	xoff = nbin / 2 + nbin / 2;
	yoff = nbin / 2 + nbin / 2;
	fresize(&fresizes, image_in, image_out, nx, ny, nlead, nxout, nyout, nxout, xoff, yoff, DRMS_MISSING_FLOAT);

}
