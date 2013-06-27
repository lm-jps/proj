/*
 *  sharp.c
 *
 *	This module creates the pipeline space weather harps
 *	It is a hard-coded strip-down version of bmap.c
 *	It takes the Mharp and Bharp series and crete the following quantities
 *  Series 1: Sharp_CEA
 *	  CEA remapped magnetogram, bitmap, continuum, doppler (same size in map coordinate, need manual spec?)
 *	  CEA remapped vector field (Br, Bt, Bp) (same as above)
 *    Space weather indices based on vector cutouts (step 2)
 *  Series 2: Sharp_cutout:
 *	  cutouts of magnetogram, bitmap, continuum, doppler (HARP defined, various sizes in CCD pixels)
 *	  cutouts of all vector data segments (same as above)
 *	Series 3: Other remaps
 *
 *	Author:
 *		Xudong Sun; Monica Bobra
 *
 *	Version:
 *		v0.0	Jul 02 2012
 *		v0.1	Jul 23 2012
 *		v0.2	Sep 04 2012
 *		v0.3	Dec 18 2012
 *		v0.4	Jan 02 2013
 *      v0.5    Jan 23 2012
 *
 *	Notes:
 *		v0.0
 *		Mharp & Bharp must be fully specified; other input are series names only
 *		All input records need to match, otherwise quit
 *		Mapping parameters depend on keywords of each record only, not necessarily consistent for now
 *		Cutout doesn't work for char segments yet (drms_segment_readslice bug)
 *		SW indices require ephemeris info which is not passed properly as of now
 *		v0.1
 *		Fixed char I/O thanks to Art
 *		SW indices fixed
 *		Added doppler and continuum
 *              Added other keywords: HEADER (populated by cvs build version), DATE_B
 *		v0.3
 *		Fixed memory leakage of 0.15G per rec; denoted with "Dec 18"
 *		v0.4
 *		Took out convert_inplace(). Was causing all the images to be int
 *      v0.5
 *      Corrected ephemeris keywords, added argument mInfo for setKeys()
 *
 *	Example:
 *	sharp "mharp=hmi.Mharp_720s[1404][2012.02.20_10:00]" \
 "bharp=hmi_test.Bharp_720s_fd10[1404][2012.02.20_10:00]" \
 "dop=hmi.V_720s[2012.02.20_10:00]" \
 "cont=hmi.Ic_720s[2012.02.20_10:00]" \
 "sharp_cea=su_xudong.Sharp_CEA" "sharp_cut=su_xudong.Sharp_Cut"
 *	For comparison:
 *	bmap "in=hmi_test.Bharp_720s_fd10[1404][2012.02.20_10:00]" \
 "out=hmi_test.B_720s_CEA" -s -a "map=cyleqa"
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
#include "sw_functions.c"

//#include <mkl.h> // Comment out mkl.h, which can only run on solar3
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

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

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

#define kNotSpecified "Not Specified"

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

#define DISAMB_AZI		1
#define XSCALE			0.03
#define YSCALE			0.03
#define NBIN			3
#define INTERP			0
#define dpath    "/home/jsoc/cvs/Development/JSOC"

/* ========================================================================================================== */

// Space weather keywords
struct swIndex {
	float mean_vf;
    float count_mask;
	float absFlux;
	float mean_hf;
	float mean_gamma;
	float mean_derivative_btotal;
	float mean_derivative_bh;
	float mean_derivative_bz;
	float mean_jz;
	float us_i;
	float mean_alpha;
	float mean_ih;
	float total_us_ih;
	float total_abs_ih;
	float totaljz;
	float totpot;
	float meanpot;
	float area_w_shear_gt_45;
	float meanshear_angle;
	float area_w_shear_gt_45h;
	float meanshear_angleh;
    float mean_derivative_btotal_err;
    float mean_vf_err;
    float mean_gamma_err;
    float mean_derivative_bh_err;
    float mean_derivative_bz_err;
    float mean_jz_err;
    float us_i_err;
    float mean_alpha_err;
    float mean_ih_err;
    float total_us_ih_err;
    float total_abs_ih_err;
    float totaljz_err;
    float meanpot_err;
    float totpot_err;
    float meanshear_angle_err;
};

// Mapping method
enum projection {
	carree,
	cassini,
	mercator,
	cyleqa,
	sineqa,
	gnomonic,
	postel,
	stereographic,
	orthographic,
	lambert
};

// WSC code
char *wcsCode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
	"SIN", "ZEA"};

// Ephemeris information
struct ephemeris {
	double disk_lonc, disk_latc;
	double disk_xc, disk_yc;
	double rSun, asd, pa;
};

// Mapping information
struct mapInfo {
	float xc, yc;		// reference point: center
	int nrow, ncol;		// size
	float xscale, yscale;	// scale
	int nbin;
	enum projection proj;	// projection method
	struct ephemeris ephem;		// ephemeris info
	float *xi_out, *zeta_out;	// coordinate on full disk image to sample at
};


/* ========================================================================================================== */

/* Get all input data series */
int getInputRS(DRMS_RecordSet_t **mharpRS_ptr, DRMS_RecordSet_t **bharpRS_ptr,
			   char *mharpQuery, char *bharpQuery);

/* Check if Mharp and Bharp match */
int compareHarp(DRMS_RecordSet_t *mharpRS, DRMS_RecordSet_t *bharpRS);

/* Get other data series */
int getInputRS_aux(DRMS_RecordSet_t **inRS_ptr, char *inQuery, DRMS_RecordSet_t *harpRS);

/* Find record from record set with given T_rec */
int getInputRec_aux(DRMS_Record_t **inRec_ptr, DRMS_RecordSet_t *inRS, TIME trec);

/* Create CEA record */
int createCeaRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *bharpRec,
					DRMS_Record_t *dopRec, DRMS_Record_t *contRec,
					DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr);

/* Mapping single segment, wrapper */
int mapScaler(DRMS_Record_t *sharpRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec,
			  struct mapInfo *mInfo, char *segName);

/* Mapping vector magnetogram, wrapper */
int mapVectorB(DRMS_Record_t *sharpRec, DRMS_Record_t *bharpRec, struct mapInfo *mInfo);

/* Mapping errors of vector magnetogram, wraper */
int mapVectorBErr(DRMS_Record_t *sharpRec, DRMS_Record_t *bharpRec, struct mapInfo *mInfo);

/* Determine reference point coordinate and patch size according to input */
int findPosition(DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Compute the coordinates at which the full disk image is sampled */
void findCoord(struct mapInfo *mInfo);

/* Mapping function */
int performSampling(float *outData, float *inData, struct mapInfo *mInfo);

/* Performing local vector transformation */
void vectorTransform(float *bx_map, float *by_map, float *bz_map, struct mapInfo *mInfo);

/* Map and propogate errors */
int getBErr(float *bx_err, float *by_err, float *bz_err,
			DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Read full disk vector magnetogram */
int readVectorB(DRMS_Record_t *inRec, float *bx_img, float *by_img, float *bz_img);

/* Read variances and covariances of vector magnetograms */
int readVectorBErr(DRMS_Record_t *bharpRec,
				   float *bT, float *bI, float *bA,
				   float *errbT, float *errbI, float *errbA,
				   float *errbTbI, float *errbTbA, float *errbIbA);

// ===================

/* Create Cutout record */
int createCutRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *bharpRec,
					DRMS_Record_t *dopRec, DRMS_Record_t *contRec,
					DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr);

/* Get cutout and write segment */
int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec, char *SegName);

// ===================

/* Compute space weather indices, no error checking for now */
void computeSWIndex(struct swIndex *swKeys_ptr, DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Set space weather indices, no error checking for now */
void setSWIndex(DRMS_Record_t *outRec, struct swIndex *swKeys_ptr);

/* Set all keywords, no error checking for now */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct mapInfo *mInfo);

// ===================

/* Nearest neighbor interpolation */
float nnb (float *f, int nx, int ny, double x, double y);

/* Wrapper for Jesper's rebin code */
void frebin (float *image_in, float *image_out, int nx, int ny, int nbin, int gauss);

/* ========================================================================================================== */

/* Remap segment names */
#define BR_SEG_CEA		"Br"
#define BT_SEG_CEA		"Bt"
#define BP_SEG_CEA		"Bp"
#define BR_ERR_SEG_CEA		"Br_err"
#define BT_ERR_SEG_CEA		"Bt_err"
#define BP_ERR_SEG_CEA		"Bp_err"

/* Cutout segment names, input identical to output */
char *MharpSegs[] = {"magnetogram", "bitmap"};
char *BharpSegs[] = {"inclination", "azimuth", "field", "vlos_mag", "dop_width", "eta_0",
	"damping", "src_continuum", "src_grad", "alpha_mag", "chisq",
	"conv_flag", // fixed
	"info_map", "confid_map",
	"inclination_err", "azimuth_err", "field_err", "vlos_err", "alpha_err",
	"field_inclination_err", "field_az_err", "inclin_azimuth_err",
	"field_alpha_err","inclination_alpha_err", "azimuth_alpha_err",
	"disambig", "conf_disambig"};
// For stats
char *CutSegs[] = {"magnetogram", "bitmap", "Dopplergram", "continuum",
	"inclination", "azimuth", "field", "vlos_mag", "dop_width", "eta_0",
	"damping", "src_continuum", "src_grad", "alpha_mag", "chisq",
	"conv_flag", // fixed
	"info_map", "confid_map",
	"inclination_err", "azimuth_err", "field_err", "vlos_err", "alpha_err",
	"field_inclination_err", "field_az_err", "inclin_azimuth_err",
	"field_alpha_err","inclination_alpha_err", "azimuth_alpha_err",
	"disambig", "conf_disambig"};
char *CEASegs[] = {"magnetogram", "bitmap", "Dopplergram", "continuum", "disambig",
	BR_SEG_CEA, BT_SEG_CEA, BP_SEG_CEA, BR_ERR_SEG_CEA, BT_ERR_SEG_CEA, BP_ERR_SEG_CEA};

/* ========================================================================================================== */

char *module_name = "sharp";
char *version_id = "2013 Jun 26";  /* Version number */
int seed;

ModuleArgs_t module_args[] =
{
	{ARG_STRING, "mharp", kNotSpecified, "Input Mharp series."},
	{ARG_STRING, "bharp", kNotSpecified, "Input Bharp series."},
	{ARG_STRING, "dop", kNotSpecified, "Input Doppler series."},
	{ARG_STRING, "cont", kNotSpecified, "Input Continuum series."},
	{ARG_STRING, "sharp_cea", kNotSpecified, "Output Sharp CEA series."},
	{ARG_STRING, "sharp_cut", kNotSpecified, "Output Sharp cutout series."},
    {ARG_INT,    "seed", "987654", "Seed for the random number generator."},
	{ARG_END}
};

int DoIt(void)
{
    int errbufstat=setvbuf(stderr, NULL, _IONBF, BUFSIZ);
    int outbufstat=setvbuf(stdout, NULL, _IONBF, BUFSIZ);
    
	int status = DRMS_SUCCESS;
	int nrecs, irec;
	
	char *mharpQuery, *bharpQuery;
	char *dopQuery, *contQuery;
	char *sharpCeaQuery, *sharpCutQuery;
	
	DRMS_RecordSet_t *mharpRS = NULL, *bharpRS = NULL;
	DRMS_RecordSet_t *dopRS = NULL, *contRS = NULL;
	
	/* Get parameters */
    
	mharpQuery = (char *) params_get_str(&cmdparams, "mharp");
	bharpQuery = (char *) params_get_str(&cmdparams, "bharp");
	dopQuery = (char *) params_get_str(&cmdparams, "dop");
	contQuery = (char *) params_get_str(&cmdparams, "cont");
	sharpCeaQuery = (char *) params_get_str(&cmdparams, "sharp_cea");
	sharpCutQuery = (char *) params_get_str(&cmdparams, "sharp_cut");
    
    seed = params_get_int(&cmdparams, "seed");
	
	/* Get input data, check everything */
	
	if (getInputRS(&mharpRS, &bharpRS, mharpQuery, bharpQuery))
		DIE("Input harp data error.");
	nrecs = mharpRS->n;
	
	if (getInputRS_aux(&dopRS, dopQuery, mharpRS))
		DIE("Input doppler data error.");
	
	if (getInputRS_aux(&contRS, contQuery, mharpRS))
		DIE("Input continuum data error.");
	
	/* Start */
	
	printf("==============\nStart. %d image(s) in total.\n", nrecs);
    
	for (irec = 0; irec < nrecs; irec++) {
		
		/* Records in work */
		
		DRMS_Record_t *mharpRec = NULL, *bharpRec = NULL;
		mharpRec = mharpRS->records[irec];
		bharpRec = bharpRS->records[irec];
		
		TIME trec = drms_getkey_time(mharpRec, "T_REC", &status);
		
		struct swIndex swKeys;
		
		DRMS_Record_t *dopRec = NULL, *contRec = NULL;
		if (getInputRec_aux(&dopRec, dopRS, trec)) {
			printf("Fetching Doppler failed, image #%d skipped.\n", irec);
			continue;
		}
		if (getInputRec_aux(&contRec, contRS, trec)) {
			printf("Fetching continuum failed, image #%d skipped.\n", irec);
			continue;
		}
        
		/* Create CEA record */
		
		DRMS_Record_t *sharpCeaRec = drms_create_record(drms_env, sharpCeaQuery, DRMS_PERMANENT, &status);
		if (status) {		// if failed
			printf("Creating CEA failed, image #%d skipped.\n", irec);
			continue;
		}
		
		if (createCeaRecord(mharpRec, bharpRec, dopRec, contRec, sharpCeaRec, &swKeys)) {		// do the work
			printf("Creating CEA failed, image #%d skipped.\n", irec);
			drms_close_record(sharpCeaRec, DRMS_FREE_RECORD);
			continue;
		}		// swKeys updated here
		
		drms_close_record(sharpCeaRec, DRMS_INSERT_RECORD);
		
		/* Create Cutout record */
		
		DRMS_Record_t *sharpCutRec = drms_create_record(drms_env, sharpCutQuery, DRMS_PERMANENT, &status);
		if (status) {		// if failed
			printf("Creating cutout failed, image #%d skipped.\n", irec);
			continue;
		}
		
		if (createCutRecord(mharpRec, bharpRec, dopRec, contRec, sharpCutRec, &swKeys)) {		// do the work
			printf("Creating cutout failed, image #%d skipped.\n", irec);
			drms_close_record(sharpCutRec, DRMS_FREE_RECORD);
			continue;
		}		// swKeys used here
		
		drms_close_record(sharpCutRec, DRMS_INSERT_RECORD);
		
		/* Done */
		
		printf("Image #%d done.\n", irec);
		
	} // irec
    
	
	drms_close_records(mharpRS, DRMS_FREE_RECORD);
	drms_close_records(bharpRS, DRMS_FREE_RECORD);
	drms_close_records(dopRS, DRMS_FREE_RECORD);				// Dec 18 2012
	drms_close_records(contRS, DRMS_FREE_RECORD);				// Dec 18 2012
	
	return 0;
	
}	// DoIt


// ===================================================================
// ===================================================================
// ===================================================================


/*
 * Get input data series, including mHarp and bharp
 * Need all records to match, otherwise quit
 *
 */

int getInputRS(DRMS_RecordSet_t **mharpRS_ptr, DRMS_RecordSet_t **bharpRS_ptr,
			   char *mharpQuery, char *bharpQuery)
{
	
	int status = 0;
	
	*mharpRS_ptr = drms_open_records(drms_env, mharpQuery, &status);
    if (status || (*mharpRS_ptr)->n == 0) return 1;
	
	*bharpRS_ptr = drms_open_records(drms_env, bharpQuery, &status);
    if (status || (*bharpRS_ptr)->n == 0) return 1;
	
	if (compareHarp((*mharpRS_ptr), (*bharpRS_ptr))) return 1;
	
	return 0;
	
}

/*
 * Check if Mharp and Bharp match
 *
 */

int compareHarp(DRMS_RecordSet_t *mharpRS, DRMS_RecordSet_t *bharpRS)
{
	
	int status = 0;
	int nrecs = mharpRS->n;
	
	DRMS_Record_t *mharpRec_t = NULL, *bharpRec_t = NULL;		// temporary recs for utility
	
    if (bharpRS->n != nrecs) {
		return 1;		// return 1 if different
	}
	
	for (int i = 0; i < nrecs; i++) {
		mharpRec_t = mharpRS->records[i];
		bharpRec_t = bharpRS->records[i];
		if ((drms_getkey_int(mharpRec_t, "HARPNUM", &status) !=
			 drms_getkey_int(bharpRec_t, "HARPNUM", &status)) ||
			(drms_getkey_time(mharpRec_t, "T_REC", &status) !=
			 drms_getkey_time(bharpRec_t, "T_REC", &status)))
		{
			return 1;
		}
	}
	
	return 0;
	
}

/*
 * Get other data series, check all T_REC are available
 *
 */

int getInputRS_aux(DRMS_RecordSet_t **inRS_ptr, char *inQuery, DRMS_RecordSet_t *harpRS)
{
	
	int status = 0;
	
	*inRS_ptr = drms_open_records(drms_env, inQuery, &status);
	if (status || (*inRS_ptr)->n == 0) return status;
	
	// Check if all T_rec are available, need to match both ways
	int n = harpRS->n, n0 = (*inRS_ptr)->n;
	
	for (int i0 = 0; i0 < n0; i0++) {
		DRMS_Record_t *inRec = (*inRS_ptr)->records[i0];
		TIME trec0 = drms_getkey_time(inRec, "T_REC", &status);
		TIME trec = 0;
		for (int i = 0; i < n; i++) {
			DRMS_Record_t *harpRec = harpRS->records[i];
			trec = drms_getkey_time(harpRec, "T_REC", &status);
			if (fabs(trec0 - trec) < 10) break;
		}
		if (fabs(trec0 - trec) >= 10) return 1;
	}
	
	for (int i = 0; i < n; i++) {
		DRMS_Record_t *harpRec = harpRS->records[i];
		TIME trec = drms_getkey_time(harpRec, "T_REC", &status);
		TIME trec0 = 0;
		for (int i0 = 0; i0 < n0; i0++) {
			DRMS_Record_t *inRec = (*inRS_ptr)->records[i0];
			trec0 = drms_getkey_time(inRec, "T_REC", &status);
			if (fabs(trec0 - trec) < 10) break;
		}
		if (fabs(trec0 - trec) >= 10) return 1;
	}
	
	return 0;
	
}

/*
 * Find record from record set with given T_rec
 *
 */

int getInputRec_aux(DRMS_Record_t **inRec_ptr, DRMS_RecordSet_t *inRS, TIME trec)
{
	
	int status = 0;
	
	int n = inRS->n;
	for (int i = 0; i < n; i++) {
		*inRec_ptr = inRS->records[i];
		TIME trec0 = drms_getkey_time((*inRec_ptr), "T_REC", &status);
		if (fabs(trec0 - trec) < 10) return 0;
	}
	
	return 1;
	
}




/*
 * Create CEA record: top level subroutine
 * Also compute all the space weather keywords here
 *
 */

int createCeaRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *bharpRec,
					DRMS_Record_t *dopRec, DRMS_Record_t *contRec,
					DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr)
{
	
	int status = 0;
	DRMS_Segment_t *inSeg;
	DRMS_Array_t *inArray;
	
	struct mapInfo mInfo;
	mInfo.proj = (enum projection) cyleqa;		// projection method
	mInfo.xscale = XSCALE;
	mInfo.yscale = YSCALE;
	mInfo.nbin = NBIN;
	
	// Get ephemeris
	
	if (getEphemeris(mharpRec, &(mInfo.ephem))) {
		SHOW("CEA: get ephemeris error\n");
		return 1;
	}
	
	// Find position
	
	if (findPosition(mharpRec, &mInfo)) {
		SHOW("CEA: find position error\n");
		return 1;
	}
	
	// Create xi_out, zeta_out array in mInfo:
	// Coordinates to sample in original full disk image
	
	int ncol0, nrow0;		// oversampled map size
	ncol0 = mInfo.ncol * mInfo.nbin + (mInfo.nbin / 2) * 2;	// pad with nbin/2 on edge to avoid NAN
	nrow0 = mInfo.nrow * mInfo.nbin + (mInfo.nbin / 2) * 2;
	
	mInfo.xi_out = (float *) (malloc(ncol0 * nrow0 * sizeof(float)));
	mInfo.zeta_out = (float *) (malloc(ncol0 * nrow0 * sizeof(float)));
	
	findCoord(&mInfo);		// compute it here so it could be shared by the following 4 functions
	
	// Mapping single segment: Mharp, etc.
    
	if (mapScaler(sharpRec, mharpRec, mharpRec, &mInfo, "magnetogram")) {
		SHOW("CEA: mapping magnetogram error\n");
		return 1;
	}
	if (mapScaler(sharpRec, mharpRec, mharpRec, &mInfo, "bitmap")) {
		SHOW("CEA: mapping bitmap error\n");
		return 1;
	}
	printf("Magnetogram mapping done.\n");
	if (mapScaler(sharpRec, dopRec, mharpRec, &mInfo, "Dopplergram")) {
		SHOW("CEA: mapping dopplergram error\n");
		return 1;
	}
	printf("Dopplergram mapping done.\n");
	
	if (mapScaler(sharpRec, contRec, mharpRec, &mInfo, "continuum")) {
		SHOW("CEA: mapping continuum error\n");
		return 1;
	}
	printf("Intensitygram mapping done.\n");
	
	if (mapScaler(sharpRec, bharpRec, mharpRec, &mInfo, "conf_disambig")) {
		SHOW("CEA: mapping conf_disambig error\n");
		return 1;
	}
	printf("Conf disambig mapping done.\n");
    
	// Mapping vector B
	
	if (mapVectorB(sharpRec, bharpRec, &mInfo)) {
		SHOW("CEA: mapping vector B error\n");
		return 1;
	}
	printf("Vector B mapping done.\n");
	
	// Mapping vector B errors
	
	if (mapVectorBErr(sharpRec, bharpRec, &mInfo)) {
		SHOW("CEA: mapping vector B uncertainty error\n");
		return 1;
	}
	printf("Vector B error done.\n");
	
	// Keywords & Links
	
	drms_copykey(sharpRec, mharpRec, "T_REC");
	drms_copykey(sharpRec, mharpRec, "HARPNUM");
	
	DRMS_Link_t *mHarpLink = hcon_lookup_lower(&sharpRec->links, "MHARP");
	if (mHarpLink) drms_link_set("MHARP", sharpRec, mharpRec);
	DRMS_Link_t *bHarpLink = hcon_lookup_lower(&sharpRec->links, "BHARP");
	if (bHarpLink) drms_link_set("BHARP", sharpRec, bharpRec);
	
	 setKeys(sharpRec, bharpRec, &mInfo);            // Set all other keywords
	drms_copykey(sharpRec, mharpRec, "QUALITY");		// copied from los records
	
	// Space weather
	
	computeSWIndex(swKeys_ptr, sharpRec, &mInfo);		// compute it!
	printf("Space weather indices done.\n");
	
	setSWIndex(sharpRec, swKeys_ptr);	// Set space weather indices
	
	// Stats
	
	int nCEASegs = ARRLENGTH(CEASegs);
	for (int iSeg = 0; iSeg < nCEASegs; iSeg++) {
		DRMS_Segment_t *outSeg = drms_segment_lookupnum(sharpRec, iSeg);
		DRMS_Array_t *outArray = drms_segment_read(outSeg, DRMS_TYPE_FLOAT, &status);
		int stat = set_statistics(outSeg, outArray, 1);
		//		printf("%d => %d\n", iSeg, stat);
		drms_free_array(outArray);
	}
	
	free(mInfo.xi_out);
	free(mInfo.zeta_out);
	return 0;
	
}


/*
 * Mapping a single segment
 * Read in full disk image, utilize mapImage for mapping
 * then write the segment out, segName same in in/out Rec
 *
 */

int mapScaler(DRMS_Record_t *sharpRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec,
			  struct mapInfo *mInfo, char *segName)
{
	
	int status = 0;
	int nx = mInfo->ncol, ny = mInfo->nrow, nxny = nx * ny;
	int dims[2] = {nx, ny};
	
	// Input full disk array
	
	DRMS_Segment_t *inSeg = NULL;
	inSeg = drms_segment_lookup(inRec, segName);
	if (!inSeg) return 1;
	
	DRMS_Array_t *inArray = NULL;
	inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (!inArray) return 1;
	
	float *inData;
	int xsz = inArray->axis[0], ysz = inArray->axis[1];
	if ((xsz != FOURK) || (ysz != FOURK)) {		// for bitmap, make tmp full disk
		float *inData0 = (float *) inArray->data;
		inData = (float *) (calloc(FOURK2, sizeof(float)));
		int x0 = (int) drms_getkey_float(harpRec, "CRPIX1", &status) - 1;
		int y0 = (int) drms_getkey_float(harpRec, "CRPIX2", &status) - 1;
		int ind_map;
		for (int row = 0; row < ysz; row++) {
			for (int col = 0; col < xsz; col++) {
				ind_map = (row + y0) * FOURK + (col + x0);
				inData[ind_map] = inData0[row * xsz + col];
			}
		}
		drms_free_array(inArray); inArray = NULL;
	} else {
		inData = (float *) inArray->data;
	}
	
	// Mapping
	
	float *map = (float *) (malloc(nxny * sizeof(float)));
	if (performSampling(map, inData, mInfo))
	{if (inArray) drms_free_array(inArray); free(map); return 1;}
	
	// Write out
	
	DRMS_Segment_t *outSeg = NULL;
	outSeg = drms_segment_lookup(sharpRec, segName);
	if (!outSeg) return 1;
	
    //	DRMS_Type_t arrayType = outSeg->info->type;
	DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, map, &status);
	if (status) {if (inArray) drms_free_array(inArray); free(map); return 1;}
	
	// convert to needed data type
	
    //	drms_array_convert_inplace(outSeg->info->type, 0, 1, outArray);		// Jan 02 2013
	
	outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    //	outArray->parent_segment = outSeg;
	outArray->israw = 0;		// always compressed
	outArray->bzero = outSeg->bzero;
	outArray->bscale = outSeg->bscale;
	
	status = drms_segment_write(outSeg, outArray, 0);
	if (status) return 0;
	
	if (inArray) drms_free_array(inArray);
	if ((xsz != FOURK) || (ysz != FOURK)) free(inData);			// Dec 18 2012
	if (outArray) drms_free_array(outArray);
	return 0;
	
}


/*
 * Mapping vector magnetogram
 *
 */

int mapVectorB(DRMS_Record_t *sharpRec, DRMS_Record_t *bharpRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	int nx = mInfo->ncol, ny = mInfo->nrow, nxny = nx * ny;
	int dims[2] = {nx, ny};
	
	// Read in segments, filling factor assume to be 1
	
	float *bx_img = (float *) (malloc(FOURK2 * sizeof(float)));
	float *by_img = (float *) (malloc(FOURK2 * sizeof(float)));
	float *bz_img = (float *) (malloc(FOURK2 * sizeof(float)));
	
	if (readVectorB(bharpRec, bx_img, by_img, bz_img)) {
		printf("Read full disk image error\n");
		free(bx_img); free(by_img); free(bz_img);
		return 1;
	}
	
	// Mapping
	
	float *bx_map = NULL, *by_map = NULL, *bz_map = NULL;	// intermediate maps, in CCD bxyz representation
	
	bx_map = (float *) (malloc(nxny * sizeof(float)));
	if (performSampling(bx_map, bx_img, mInfo))
	{free(bx_img); free(by_img); free(bz_img); free(bx_map); return 1;}
	
	by_map = (float *) (malloc(nxny * sizeof(float)));
	if (performSampling(by_map, by_img, mInfo))
	{free(bx_img); free(by_img); free(bz_img); free(bz_map); return 1;}
	
	bz_map = (float *) (malloc(nxny * sizeof(float)));
	if (performSampling(bz_map, bz_img, mInfo))
	{free(bx_img); free(by_img); free(bz_img); free(bz_map); return 1;}
	
	free(bx_img); free(by_img); free(bz_img);
	
	// Vector transform
	
	vectorTransform(bx_map, by_map, bz_map, mInfo);
	
	for (int i = 0; i < nxny; i++) by_map[i] *= -1;		// positive theta pointing south
	
	// Write out
	
	DRMS_Segment_t *outSeg;
	DRMS_Array_t *outArray;
	
	float *data_prt[3] = {bx_map, by_map, bz_map};
	char *segName[3] = {BP_SEG_CEA, BT_SEG_CEA, BR_SEG_CEA};
	
	for (int iSeg = 0; iSeg < 3; iSeg++) {
		outSeg = drms_segment_lookup(sharpRec, segName[iSeg]);
		outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, data_prt[iSeg], &status);
		outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        //		outArray->parent_segment = outSeg;
		outArray->israw = 0;
		outArray->bzero = outSeg->bzero;
		outArray->bscale = outSeg->bscale;
		status = drms_segment_write(outSeg, outArray, 0);
		if (status) return 1;
		drms_free_array(outArray);
	}
	
	//
	
	return 0;
	
}


/*
 * Mapping vector magnetogram errors
 *
 */

int mapVectorBErr(DRMS_Record_t *sharpRec, DRMS_Record_t *bharpRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	
	int nx = mInfo->ncol, ny = mInfo->nrow, nxny = nx * ny;
	int dims[2] = {nx, ny};
	
	// Compute propogated errors, using nearest neighbour interpolation
	
	float *bx_err = (float *) (malloc(nxny * sizeof(float)));
	float *by_err = (float *) (malloc(nxny * sizeof(float)));
	float *bz_err = (float *) (malloc(nxny * sizeof(float)));
	
	if (getBErr(bx_err, by_err, bz_err, bharpRec, mInfo)) {
		free(bx_err); free(by_err); free(bz_err);
		return 1;
	}
	
	// Write out
	
	DRMS_Segment_t *outSeg;
	DRMS_Array_t *outArray;
	
	float *data_prt[3] = {bx_err, by_err, bz_err};
	char *segName[3] = {BP_ERR_SEG_CEA, BT_ERR_SEG_CEA, BR_ERR_SEG_CEA};
	
	for (int iSeg = 0; iSeg < 3; iSeg++) {
		outSeg = drms_segment_lookup(sharpRec, segName[iSeg]);
		outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, data_prt[iSeg], &status);
		outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        //		outArray->parent_segment = outSeg;
		outArray->israw = 0;
		outArray->bzero = outSeg->bzero;
		outArray->bscale = outSeg->bscale;
		status = drms_segment_write(outSeg, outArray, 0);
		if (status) return 1;
		drms_free_array(outArray);
	}
	
	//
	
	return 0;
	
}



/*
 * Determine reference point coordinate and patch size according to keywords
 * xc, yc are the coordinate of patch center, in degrees
 * ncol and nrow are the final size
 *
 */

int findPosition(DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	int harpnum = drms_getkey_int(inRec, "HARPNUM", &status);
	TIME trec = drms_getkey_time(inRec, "T_REC", &status);
	float disk_lonc = drms_getkey_float(inRec, "CRLN_OBS", &status);
	
	/* Center coord */
	
	float minlon = drms_getkey_float(inRec, "LONDTMIN", &status); if (status) return 1;		// Stonyhurst lon
	float maxlon = drms_getkey_float(inRec, "LONDTMAX", &status); if (status) return 1;
	float minlat = drms_getkey_float(inRec, "LATDTMIN", &status); if (status) return 1;
	float maxlat = drms_getkey_float(inRec, "LATDTMAX", &status); if (status) return 1;
	
	// A bug fixer for HARP (per M. Turmon)
	// When AR is below threshold, "LONDTMIN", "LONDTMAX" will be wrong
	// Also keywords such as "SIZE" will be NaN
	// We compute minlon & minlat then by
	// LONDTMIN(t) = LONDTMIN(t0) + (t - t0) * OMEGA_DT
	
    //	float psize = drms_getkey_float(inRec, "SIZE", &status);
    //	if (psize != psize) {
    
    if (minlon != minlon || maxlon != maxlon) {		// check lons instead of SIZE
		TIME t0 = drms_getkey_time(inRec, "T_FRST1", &status); if (status) return 1;			// changed from T_FRST to T_FRST1, T_FRST may not exist
		double omega = drms_getkey_double(inRec, "OMEGA_DT", &status); if (status) return 1;
		char firstRecQuery[100], t0_str[100];
		sprint_time(t0_str, t0, "TAI", 0);
		snprintf(firstRecQuery, 100, "%s[%d][%s]", inRec->seriesinfo->seriesname, harpnum, t0_str);
		DRMS_RecordSet_t *tmpRS = drms_open_records(drms_env, firstRecQuery, &status);
		if (status || tmpRS->n != 1) return 1;
		DRMS_Record_t *tmpRec = tmpRS->records[0];
		double minlon0 = drms_getkey_double(tmpRec, "LONDTMIN", &status); if (status) return 1;
		double maxlon0 = drms_getkey_double(tmpRec, "LONDTMAX", &status); if (status) return 1;
		minlon = minlon0 + (trec - t0) * omega / SECINDAY;
		maxlon = maxlon0 + (trec - t0) * omega / SECINDAY;
		printf("%s, %f, %f\n", firstRecQuery, minlon, maxlon);
	}
	
	mInfo->xc = (maxlon + minlon) / 2. + disk_lonc;
	mInfo->yc = (maxlat + minlat) / 2.;
	
	/* Size */
	
	mInfo->ncol = round((maxlon - minlon) / mInfo->xscale);
	mInfo->nrow = round((maxlat - minlat) / mInfo->yscale);
	
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
	float crpix1 = drms_getkey_float(inRec, "IMCRPIX1", &status);
	float crpix2 = drms_getkey_float(inRec, "IMCRPIX2", &status);
	float cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
	printf("cdelt=%f\n",cdelt);
	ephem->disk_xc = PIX_X(0.0,0.0) - 1.0;		// Center of disk in pixel, starting at 0
	ephem->disk_yc = PIX_Y(0.0,0.0) - 1.0;
	
	float dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
	float rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
	if (status) rSun_ref = 6.96e8;
	
	ephem->asd = asin(rSun_ref/dSun);
	ephem->rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
	
	return 0;
	
}


/*
 * Compute the coordinates to be sampled on full disk image
 * mInfo->xi_out & mInfo->zeta_out
 * This is oversampled, its size is ncol0 & nrow0 as shown below
 *
 *
 */

void findCoord(struct mapInfo *mInfo)
{
	
	int ncol0 = mInfo->ncol * mInfo->nbin + (mInfo->nbin / 2) * 2;	// pad with nbin/2 on edge to avoid NAN
	int nrow0 = mInfo->nrow * mInfo->nbin + (mInfo->nbin / 2) * 2;
	
	float xscale0 = mInfo->xscale / mInfo->nbin * RADSINDEG;		// oversampling resolution
	float yscale0 = mInfo->yscale / mInfo->nbin * RADSINDEG;		// in rad
	
	double lonc = mInfo->xc * RADSINDEG;	// in rad
	double latc = mInfo->yc * RADSINDEG;
	
	double disk_lonc = (mInfo->ephem).disk_lonc;
	double disk_latc = (mInfo->ephem).disk_latc;
	
	double rSun = (mInfo->ephem).rSun;
	double disk_xc = (mInfo->ephem).disk_xc / rSun;
	double disk_yc = (mInfo->ephem).disk_yc / rSun;
	double pa = (mInfo->ephem).pa;
	
	// Temp pointers
	
	float *xi_out = mInfo->xi_out;
	float *zeta_out = mInfo->zeta_out;
	
	// start
	
	double x, y;		// map coord
	double lat, lon;	// helio coord
	double xi, zeta;	// image coord (for one point)
	
	int ind_map;
	
	for (int row0 = 0; row0 < nrow0; row0++) {
		for (int col0 = 0; col0 < ncol0; col0++) {
			
			ind_map = row0 * ncol0 + col0;
			
			x = (col0 + 0.5 - ncol0/2.) * xscale0;		// in rad
			y = (row0 + 0.5 - nrow0/2.) * yscale0;
			
			/* map grid [x, y] corresponds to the point [lon, lat] in the heliographic coordinates.
			 * the [x, y] are in radians with respect of the center of the map [xcMap, ycMap].
			 * projection methods could be Mercator, Lambert, and many others. [maplonc, mapLatc]
			 * is the heliographic longitude and latitude of the map center. Both are in degree.
			 */
			
			if (plane2sphere (x, y, latc, lonc, &lat, &lon, (int) mInfo->proj)) {
				xi_out[ind_map] = -1;
				zeta_out[ind_map] = -1;
				continue;
			}
			
			/* map the grid [lon, lat] in the heliographic coordinates to [xi, zeta], a point in the
			 * image coordinates. The image properties, xCenter, yCenter, rSun, pa, ecc and chi are given.
			 */
			
			if (sphere2img (lat, lon, disk_latc, disk_lonc, &xi, &zeta,
							disk_xc, disk_yc, 1.0, pa, 0., 0., 0., 0.)) {
				xi_out[ind_map] = -1;
				zeta_out[ind_map] = -1;
				continue;
			}
			
			xi_out[ind_map] = xi * rSun;
			zeta_out[ind_map] = zeta * rSun;
			
		}
	}
	
}


/*
 * Sampling function
 * oversampling by nbin, then binning using a Gaussian
 * save results in outData, always of float type
 *
 */

int performSampling(float *outData, float *inData, struct mapInfo *mInfo)
{
	
	int status = 0;
	
	int ncol0 = mInfo->ncol * mInfo->nbin + (mInfo->nbin / 2) * 2;	// pad with nbin/2 on edge to avoid NAN
	int nrow0 = mInfo->nrow * mInfo->nbin + (mInfo->nbin / 2) * 2;
	
	float *outData0 = (float *) (malloc(ncol0 * nrow0 * sizeof(float)));
	
	float *xi_out = mInfo->xi_out;
	float *zeta_out = mInfo->zeta_out;
	
	// Interpolation
	
	struct fint_struct pars;
	int interpOpt = INTERP;		// Use Wiener by default, 6 order, 1 constraint
	
	switch (interpOpt) {
		case 0:			// Wiener, 6 order, 1 constraint
			init_finterpolate_wiener(&pars, 6, 1, 6, 2, 1, 1, NULL, dpath);
			break;
		case 1:			// Cubic convolution
			init_finterpolate_cubic_conv(&pars, 1., 3.);
			break;
		case 2:			// Bilinear
			init_finterpolate_linear(&pars, 1.);
			break;
		default:
			return 1;
	}
	
	finterpolate(&pars, inData, xi_out, zeta_out, outData0,
				 FOURK, FOURK, FOURK, ncol0, nrow0, ncol0, DRMS_MISSING_FLOAT);
	
	// Rebinning, smoothing
	
	frebin(outData0, outData, ncol0, nrow0, mInfo->nbin, 1);		// Gaussian
	free(outData0);		// Dec 18 2012
	
	//
	
	return 0;
	
}


/*
 * Performing local vector transformation
 *  xyz: z refers to vertical (radial) component, x EW (phi), y NS
 *
 */

void vectorTransform(float *bx_map, float *by_map, float *bz_map, struct mapInfo *mInfo)
{
	
	int ncol = mInfo->ncol;
	int nrow = mInfo->nrow;
	
	float xscale = mInfo->xscale * RADSINDEG;		// in rad
	float yscale = mInfo->yscale * RADSINDEG;
	
	double lonc = mInfo->xc * RADSINDEG;	// in rad
	double latc = mInfo->yc * RADSINDEG;
	
	double disk_lonc = (mInfo->ephem).disk_lonc;
	double disk_latc = (mInfo->ephem).disk_latc;
	
	double rSun = (mInfo->ephem).rSun;
	double disk_xc = (mInfo->ephem).disk_xc / rSun;
	double disk_yc = (mInfo->ephem).disk_yc / rSun;
	double pa = (mInfo->ephem).pa;
	
	int ind_map;
	double x, y;
	double lat, lon;	// lat / lon for current point
	
	double bx_tmp, by_tmp, bz_tmp;
	
	//
	
	for (int row = 0; row < mInfo->nrow; row++) {
		for (int col = 0; col < mInfo->ncol; col++) {
			
			ind_map = row * mInfo->ncol + col;
			
			x = (col + 0.5 - ncol / 2.) * xscale;
			y = (row + 0.5 - nrow / 2.) * yscale;
			
			if (plane2sphere (x, y, latc, lonc, &lat, &lon, (int) mInfo->proj)) {
				bx_map[ind_map] = DRMS_MISSING_FLOAT;
				by_map[ind_map] = DRMS_MISSING_FLOAT;
				bz_map[ind_map] = DRMS_MISSING_FLOAT;
				continue;
			}
			
			bx_tmp = by_tmp = bz_tmp = 0;
			
			img2helioVector (bx_map[ind_map], by_map[ind_map], bz_map[ind_map],
							 &bx_tmp, &by_tmp, &bz_tmp,
							 lon, lat, disk_lonc, disk_latc, pa);
			
			bx_map[ind_map] = bx_tmp;
			by_map[ind_map] = by_tmp;
			bz_map[ind_map] = bz_tmp;
			
		}
	}
	
}



/*
 * Map and propogate vector field errors
 *
 */

int getBErr(float *bx_err, float *by_err, float *bz_err,
			DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	
	// Get variances and covariances, filling factor assume to be 1
	
	float *bT = (float *) (malloc(FOURK2 * sizeof(float)));	// field
	float *bI = (float *) (malloc(FOURK2 * sizeof(float)));	// inclination
	float *bA = (float *) (malloc(FOURK2 * sizeof(float)));	// azimuth
	
	float *errbT = (float *) (malloc(FOURK2 * sizeof(float)));
	float *errbI = (float *) (malloc(FOURK2 * sizeof(float)));
	float *errbA = (float *) (malloc(FOURK2 * sizeof(float)));
	
	float *errbTbI = (float *) (malloc(FOURK2 * sizeof(float)));
	float *errbTbA = (float *) (malloc(FOURK2 * sizeof(float)));
	float *errbIbA = (float *) (malloc(FOURK2 * sizeof(float)));
	
	if (readVectorBErr(inRec,
					   bT, bI, bA,
					   errbT, errbI, errbA,
					   errbTbI, errbTbA, errbIbA)) {
		printf("Read full disk variances & covariances error\n");
		free(bT); free(bI); free(bA);
		free(errbT); free(errbI); free(errbA);
		free(errbTbI); free(errbTbA); free(errbIbA);
		return 1;
	}
	
	// Size
	
	int ncol = mInfo->ncol;
	int nrow = mInfo->nrow;
	
	float xscale = mInfo->xscale * RADSINDEG;		// in rad
	float yscale = mInfo->yscale * RADSINDEG;
	
	double lonc = mInfo->xc * RADSINDEG;	// in rad
	double latc = mInfo->yc * RADSINDEG;
	
	double disk_lonc = (mInfo->ephem).disk_lonc;
	double disk_latc = (mInfo->ephem).disk_latc;
	
	double rSun = (mInfo->ephem).rSun;
	double disk_xc = (mInfo->ephem).disk_xc / rSun;
	double disk_yc = (mInfo->ephem).disk_yc / rSun;
	double pa = (mInfo->ephem).pa;
	
	// Start
	
	double x, y;          // map coord
	double lat, lon;      // spherical coord
	double xi, zeta;      // image coord, round to full pixel
	
	int ind_map, ind_img;
	
	double bpSigma2, btSigma2, brSigma2;		// variances after prop
	
	for (int row = 0; row < mInfo->nrow; row++) {
		for (int col = 0; col < mInfo->ncol; col++) {
			
			ind_map = row * mInfo->ncol + col;
			
			x = (col + 0.5 - ncol / 2.) * xscale;
			y = (row + 0.5 - nrow / 2.) * yscale;
			
			if (plane2sphere (x, y, latc, lonc, &lat, &lon, (int) mInfo->proj)) {
				bx_err[ind_map] = DRMS_MISSING_FLOAT;
				by_err[ind_map] = DRMS_MISSING_FLOAT;
				bz_err[ind_map] = DRMS_MISSING_FLOAT;
				continue;
			}
			
			if (sphere2img (lat, lon, disk_latc, disk_lonc, &xi, &zeta,
							disk_xc, disk_yc, 1.0, pa, 0., 0., 0., 0.)) {
				bx_err[ind_map] = DRMS_MISSING_FLOAT;
				bx_err[ind_map] = DRMS_MISSING_FLOAT;
				bx_err[ind_map] = DRMS_MISSING_FLOAT;
				continue;
			}
			
			xi *= rSun; xi = round(xi);
			zeta *= rSun; zeta = round(zeta);     // nearest neighbor
			
			ind_img = round(zeta * FOURK + xi);
			
			if (errorprop(bT, bA, bI,
						  errbT, errbA, errbI, errbTbA, errbTbI, errbIbA,
						  lon, lat, disk_lonc, disk_latc, pa, FOURK, FOURK, xi, zeta,
						  &btSigma2, &bpSigma2, &brSigma2)) {
				bx_err[ind_map] = DRMS_MISSING_FLOAT;
				by_err[ind_map] = DRMS_MISSING_FLOAT;
				bz_err[ind_map] = DRMS_MISSING_FLOAT;
				continue;
			}
			
			bx_err[ind_map] = sqrt(bpSigma2);
			by_err[ind_map] = sqrt(btSigma2);
			bz_err[ind_map] = sqrt(brSigma2);
			
		}
	}
	
	//
	
	free(bT); free(bI); free(bA);
	free(errbT); free(errbI); free(errbA);
	free(errbTbI); free(errbTbA); free(errbIbA);
	return 0;
	
}



/*
 * Read full disk vector magnetograms
 * Fill factor is 1, use default disambiguity resolution
 *
 */

int readVectorB(DRMS_Record_t *inRec, float *bx_img, float *by_img, float *bz_img)
{
	
	int status = 0;
	
	DRMS_Segment_t *inSeg;
	DRMS_Array_t *inArray_ambig;
	DRMS_Array_t *inArray_bTotal, *inArray_bAzim, *inArray_bIncl;
	
	char *ambig;
	float *bTotal, *bAzim, *bIncl;
	
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
	
	// Convert CCD xyz
	
	int llx, lly;		// lower-left corner
	int bmx, bmy;		// bitmap size
	
	llx = (int)(drms_getkey_float(inRec, "CRPIX1", &status)) - 1;
	lly = (int)(drms_getkey_float(inRec, "CRPIX2", &status)) - 1;
	
	bmx = inArray_ambig->axis[0];
	bmy = inArray_ambig->axis[1];
	
	int kx, ky, kOff;
	int ix = 0, jy = 0, yOff = 0, iData = 0;
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
			
			// zero azi pointing up, zero incl pointing out from sun
			bx_img[iData] = - bTotal[iData] * sin(bIncl[iData] * RADSINDEG) * sin(bAzim[iData] * RADSINDEG);
			by_img[iData] = bTotal[iData] * sin(bIncl[iData] * RADSINDEG) * cos(bAzim[iData] * RADSINDEG);
			bz_img[iData] = bTotal[iData] * cos(bIncl[iData] * RADSINDEG);
            
			// Disambiguation
			
			if (kx < 0 || kx >= bmx || ky < 0 || ky >= bmy) {
				continue;
			} else {
				kOff = ky * bmx + kx;
				if (ambig[kOff] % 2) {		// 180
					bx_img[iData] *= -1.; by_img[iData] *= -1.;
				}
			}
		}
	}
	
	// Clean up
	
	drms_free_array(inArray_ambig);
	drms_free_array(inArray_bTotal);
	drms_free_array(inArray_bAzim);
	drms_free_array(inArray_bIncl);
	
	return 0;
	
}


/*
 * Read variances and covariances of vector magnetograms
 *
 */

int readVectorBErr(DRMS_Record_t *inRec,
				   float *bT, float *bI, float *bA,
				   float *errbT, float *errbI, float *errbA,
				   float *errbTbI, float *errbTbA, float *errbIbA)
{
	
	int status = 0;
	
	float *data_ptr[9];
	char *segName[9] = {"field", "inclination", "azimuth",
		"field_err", "inclination_err", "azimuth_err",
		"field_inclination_err", "field_az_err", "inclin_azimuth_err"};
	DRMS_Segment_t *inSegs[9];
	DRMS_Array_t *inArrays[9];
	
	// Read full disk images
	
	for (int iSeg = 0; iSeg < 9; iSeg++) {
		
		inSegs[iSeg] = drms_segment_lookup(inRec, segName[iSeg]);
		inArrays[iSeg] = drms_segment_read(inSegs[iSeg], DRMS_TYPE_FLOAT, &status);
		data_ptr[iSeg] = (float *) inArrays[iSeg]->data;
		
	}
	
	float *bT0 = data_ptr[0], *bI0 = data_ptr[1], *bA0 = data_ptr[2];
	float *errbT0 = data_ptr[3], *errbI0 = data_ptr[4], *errbA0 = data_ptr[5];
	float *errbTbI0 = data_ptr[6], *errbTbA0 = data_ptr[7], *errbIbA0 = data_ptr[8];
	
	// Convert errors to variances, correlation coefficients to covariances
	
	for (int i = 0; i < FOURK2; i++) {
		
		if (fabs(errbI0[i]) > 180.) errbI0[i] = 180.;
		if (fabs(errbA0[i]) > 180.) errbA0[i] = 180.;
		
		bT[i] = bT0[i];
		bI[i] = bI0[i];
		bA[i] = bA0[i];
		
		errbT[i] = errbT0[i] * errbT0[i];
		errbI[i] = errbI0[i] * errbI0[i] * RADSINDEG * RADSINDEG;
		errbA[i] = errbA0[i] * errbA0[i] * RADSINDEG * RADSINDEG;
		
		errbTbI[i] = errbTbI0[i] * errbT0[i] * errbI0[i] * RADSINDEG;
        errbTbA[i] = errbTbA0[i] * errbT0[i] * errbA0[i] * RADSINDEG;
        errbIbA[i] = errbIbA0[i] * errbI0[i] * errbA0[i] * RADSINDEG * RADSINDEG;
		
	}
	
	//
	
	for (int iSeg = 0; iSeg < 9; iSeg++) drms_free_array(inArrays[iSeg]);
	
	return 0;
	
}


/*
 * Create Cutout record: top level subroutine
 * Do the loops on segments and set the keywords here
 * Work is done in writeCutout routine below
 *
 */

int createCutRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *bharpRec,
					DRMS_Record_t *dopRec, DRMS_Record_t *contRec,
					DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr)
{
	
	int status = 0;
	
	int iHarpSeg;
	int nMharpSegs = ARRLENGTH(MharpSegs), nBharpSegs = ARRLENGTH(BharpSegs);
	
	// Cutout Mharp
	
	for (iHarpSeg = 0; iHarpSeg < nMharpSegs; iHarpSeg++) {
		if (writeCutout(sharpRec, mharpRec, mharpRec, MharpSegs[iHarpSeg])) {
			printf("Mharp cutout fails for %s\n", MharpSegs[iHarpSeg]);
			break;
		}
	}
	if (iHarpSeg != nMharpSegs) {
		SHOW("Cutout: segment number unmatch\n");
		return 1;		// if failed
	}
	printf("Magnetogram cutout done.\n");
	
	// Cutout Doppler
	
	if (writeCutout(sharpRec, dopRec, mharpRec, "Dopplergram")) {
		printf("Doppler cutout failed\n");
		return 1;
	}
	printf("Dopplergram cutout done.\n");
	
	// Cutout Continuum
	
	if (writeCutout(sharpRec, contRec, mharpRec, "continuum")) {
		printf("Continuum cutout failed\n");
		return 1;
	}
	printf("Intensitygram cutout done.\n");
	
	// Coutout Bharp
	
	for (iHarpSeg = 0; iHarpSeg < nBharpSegs; iHarpSeg++) {
		if (writeCutout(sharpRec, bharpRec, mharpRec, BharpSegs[iHarpSeg])) {
			printf("Bharp cutout fails for %s\n", BharpSegs[iHarpSeg]);
			break;
		}
	}
	if (iHarpSeg != nBharpSegs) return 1;		// if failed
	printf("Vector B cutout done.\n");
	
	// Keywords & Links
	
	drms_copykey(sharpRec, mharpRec, "T_REC");
	drms_copykey(sharpRec, mharpRec, "HARPNUM");
	
	DRMS_Link_t *mHarpLink = hcon_lookup_lower(&sharpRec->links, "MHARP");
	if (mHarpLink) drms_link_set("MHARP", sharpRec, mharpRec);
	DRMS_Link_t *bHarpLink = hcon_lookup_lower(&sharpRec->links, "BHARP");
	if (bHarpLink) drms_link_set("BHARP", sharpRec, bharpRec);
	
	setSWIndex(sharpRec, swKeys_ptr);	// Set space weather indices
	setKeys(sharpRec, bharpRec, NULL);              // Set all other keywords, NULL specifies cutout
	
	// Stats
	
	int nCutSegs = ARRLENGTH(CutSegs);
	for (int iSeg = 0; iSeg < nCutSegs; iSeg++) {
		DRMS_Segment_t *outSeg = drms_segment_lookupnum(sharpRec, iSeg);
		DRMS_Array_t *outArray = drms_segment_read(outSeg, DRMS_TYPE_FLOAT, &status);
		set_statistics(outSeg, outArray, 1);
		drms_free_array(outArray);
	}
	
	return 0;
	
}


/*
 * Get cutout and write segment
 * Change DISAMB_AZI to apply disambiguation to azimuth
 *
 */

int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec, char *SegName)
{
	
	int status = 0;
	
	DRMS_Segment_t *inSeg = NULL, *outSeg = NULL;
	DRMS_Array_t *cutoutArray = NULL;
	//	DRMS_Type_t arrayType;
	
	int ll[2], ur[2], nx, ny, nxny;		// lower-left and upper right coords
	
	/* Info */
	
	inSeg = drms_segment_lookup(inRec, SegName);
	if (!inSeg) return 1;
	
	nx = (int) drms_getkey_float(harpRec, "CRSIZE1", &status);
	ny = (int) drms_getkey_float(harpRec, "CRSIZE2", &status);
	nxny = nx * ny;
	ll[0] = (int) drms_getkey_float(harpRec, "CRPIX1", &status) - 1; if (status) return 1;
	ll[1] = (int) drms_getkey_float(harpRec, "CRPIX2", &status) - 1; if (status) return 1;
	ur[0] = ll[0] + nx - 1; if (status) return 1;
	ur[1] = ll[1] + ny - 1; if (status) return 1;
	
	if (inSeg->axis[0] == nx && inSeg->axis[1] == ny) {			// for bitmaps, infomaps, etc.
		cutoutArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
		if (status) return 1;
	} else if (inSeg->axis[0] == FOURK && inSeg->axis[1] == FOURK) {		// for full disk ones
		cutoutArray = drms_segment_readslice(inSeg, DRMS_TYPE_DOUBLE, ll, ur, &status);
		if (status) return 1;
	} else {
		return 1;
	}
	
	/* Adding disambiguation resolution to cutout azimuth? */
	
#if DISAMB_AZI
	if (!strcmp(SegName, "azimuth")) {
		DRMS_Segment_t *disambSeg = NULL;
		disambSeg = drms_segment_lookup(inRec, "disambig");
		if (!disambSeg) {drms_free_array(cutoutArray); return 1;}
		DRMS_Array_t *disambArray;
		if (disambSeg->axis[0] == nx && disambSeg->axis[1] == ny) {
			disambArray = drms_segment_read(disambSeg, DRMS_TYPE_CHAR, &status);
			if (status) {drms_free_array(cutoutArray); return 1;}
		} else {
            drms_free_array(cutoutArray);
			return 1;
		}
		double *azimuth = (double *) cutoutArray->data;
		char *disamb = (char *) disambArray->data;
		for (int n = 0; n < nxny; n++) {
			if (disamb[n]) azimuth[n] += 180.;
		}
		drms_free_array(disambArray);
	}
#endif
	
	/* Write out */
	
	outSeg = drms_segment_lookup(outRec, SegName);
	if (!outSeg) return 1;
    //	drms_array_convert_inplace(outSeg->info->type, 0, 1, cutoutArray);	// Jan 02 2013
	outSeg->axis[0] = cutoutArray->axis[0];
	outSeg->axis[1] = cutoutArray->axis[1];
    //	cutoutArray->parent_segment = outSeg;
	cutoutArray->israw = 0;		// always compressed
    cutoutArray->bzero = outSeg->bzero;
    cutoutArray->bscale = outSeg->bscale;		// Same as inArray's
	status = drms_segment_write(outSeg, cutoutArray, 0);
	drms_free_array(cutoutArray);
	if (status) return 1;
	
	return 0;
	
}


/*
 * Compute space weather indices, no error checking for now
 * Based on M. Bobra's swharp_vectorB.c
 * No error checking for now
 *
 */

void computeSWIndex(struct swIndex *swKeys_ptr, DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	int nx = mInfo->ncol, ny = mInfo->nrow;
	int nxny = nx * ny;
	int dims[2] = {nx, ny};
    
	// Get bx, by, bz, mask
	
	// Use HARP (Turmon) bitmap as a threshold on spaceweather quantities
	DRMS_Segment_t *bitmaskSeg = drms_segment_lookup(inRec, "bitmap");
	DRMS_Array_t *bitmaskArray = drms_segment_read(bitmaskSeg, DRMS_TYPE_INT, &status);
	int *bitmask = (int *) bitmaskArray->data;		// get the previously made mask array
	
	//Use conf_disambig map as a threshold on spaceweather quantities
	DRMS_Segment_t *maskSeg = drms_segment_lookup(inRec, "conf_disambig");
	DRMS_Array_t *maskArray = drms_segment_read(maskSeg, DRMS_TYPE_INT, &status);
	int *mask = (int *) maskArray->data;		// get the previously made mask array
	
	DRMS_Segment_t *bxSeg = drms_segment_lookup(inRec, BP_SEG_CEA);
	DRMS_Array_t *bxArray = drms_segment_read(bxSeg, DRMS_TYPE_FLOAT, &status);
	float *bx = (float *) bxArray->data;		// bx
	
	DRMS_Segment_t *bySeg = drms_segment_lookup(inRec, BT_SEG_CEA);
	DRMS_Array_t *byArray = drms_segment_read(bySeg, DRMS_TYPE_FLOAT, &status);
	float *by = (float *) byArray->data;		// by
	for (int i = 0; i < nxny; i++) by[i] *= -1;
	
	DRMS_Segment_t *bzSeg = drms_segment_lookup(inRec, BR_SEG_CEA);
	DRMS_Array_t *bzArray = drms_segment_read(bzSeg, DRMS_TYPE_FLOAT, &status);
	float *bz = (float *) bzArray->data;		// bz
    
	DRMS_Segment_t *bz_errSeg = drms_segment_lookup(inRec, BR_ERR_SEG_CEA);
	DRMS_Array_t *bz_errArray = drms_segment_read(bz_errSeg, DRMS_TYPE_FLOAT, &status);
	float *bz_err = (float *) bz_errArray->data;		// bz_err
    
	DRMS_Segment_t *by_errSeg = drms_segment_lookup(inRec, BT_ERR_SEG_CEA);
	DRMS_Array_t *by_errArray = drms_segment_read(by_errSeg, DRMS_TYPE_FLOAT, &status);
	float *by_err = (float *) by_errArray->data;		// by_err
	//for (int i = 0; i < nxny; i++) by_err[i] *= -1;
    
	DRMS_Segment_t *bx_errSeg = drms_segment_lookup(inRec, BP_ERR_SEG_CEA);
	DRMS_Array_t *bx_errArray = drms_segment_read(bx_errSeg, DRMS_TYPE_FLOAT, &status);
	float *bx_err = (float *) bx_errArray->data;		// bx_err
	
	// Get emphemeris
	
	//float cdelt1_orig = drms_getkey_float(inRec, "CDELT1",   &status);
	float  cdelt1      = drms_getkey_float(inRec, "CDELT1",   &status);
	float  dsun_obs    = drms_getkey_float(inRec, "DSUN_OBS",   &status);
	double rsun_ref   = drms_getkey_double(inRec, "RSUN_REF", &status);
	double rsun_obs   = drms_getkey_double(inRec, "RSUN_OBS", &status);
	float imcrpix1    = drms_getkey_float(inRec, "IMCRPIX1", &status);
	float imcrpix2    = drms_getkey_float(inRec, "IMCRPIX2", &status);
	float crpix1      = drms_getkey_float(inRec, "CRPIX1", &status);
	float crpix2      = drms_getkey_float(inRec, "CRPIX2", &status);
	
	// Temp arrays
	
	float *bh      = (float *) (malloc(nxny * sizeof(float)));
	float *bt      = (float *) (malloc(nxny * sizeof(float)));
	float *jz      = (float *) (malloc(nxny * sizeof(float)));
	float *jz_smooth = (float *) (malloc(nxny * sizeof(float)));
	float *bpx     = (float *) (malloc(nxny * sizeof(float)));
	float *bpy     = (float *) (malloc(nxny * sizeof(float)));
	float *bpz     = (float *) (malloc(nxny * sizeof(float)));
	float *derx    = (float *) (malloc(nxny * sizeof(float)));
	float *dery    = (float *) (malloc(nxny * sizeof(float)));
	float *derx_bt = (float *) (malloc(nxny * sizeof(float)));
	float *dery_bt = (float *) (malloc(nxny * sizeof(float)));
	float *derx_bh = (float *) (malloc(nxny * sizeof(float)));
	float *dery_bh = (float *) (malloc(nxny * sizeof(float)));
	float *derx_bz = (float *) (malloc(nxny * sizeof(float)));
	float *dery_bz = (float *) (malloc(nxny * sizeof(float)));
	float *bt_err  = (float *) (malloc(nxny * sizeof(float)));
	float *bh_err  = (float *) (malloc(nxny * sizeof(float)));
    float *jz_err  = (float *) (malloc(nxny * sizeof(float)));
    float *jz_err_squared = (float *) (malloc(nxny * sizeof(float)));
    float *jz_err_squared_smooth = (float *) (malloc(nxny * sizeof(float)));
    float *jz_rms_err = (float *) (malloc(nxny * sizeof(float)));
	//spaceweather quantities computed
    
    
	if (computeAbsFlux(bz_err, bz , dims, &(swKeys_ptr->absFlux), &(swKeys_ptr->mean_vf),  &(swKeys_ptr->mean_vf_err),
                       &(swKeys_ptr->count_mask), mask, bitmask, cdelt1, rsun_ref, rsun_obs))
    {
		swKeys_ptr->absFlux = DRMS_MISSING_FLOAT;		// If fail, fill in NaN
		swKeys_ptr->mean_vf = DRMS_MISSING_FLOAT;
        swKeys_ptr->mean_vf_err = DRMS_MISSING_FLOAT;
        swKeys_ptr->count_mask  = DRMS_MISSING_INT;
	}
    
	for (int i = 0; i < nxny; i++) bpz[i] = bz[i];
	greenpot(bpx, bpy, bpz, nx, ny);
	
	computeBh(bx_err, by_err, bh_err, bx, by, bz, bh, dims, &(swKeys_ptr->mean_hf), mask, bitmask);
    
	if (computeGamma(bz_err, bh_err, bx, by, bz, bh, dims, &(swKeys_ptr->mean_gamma), &(swKeys_ptr->mean_gamma_err),mask, bitmask))
	{
        swKeys_ptr->mean_gamma     =  DRMS_MISSING_FLOAT;
        swKeys_ptr->mean_gamma_err =  DRMS_MISSING_FLOAT;
    }
	
	computeB_total(bx_err, by_err, bz_err, bt_err, bx, by, bz, bt, dims, mask, bitmask);
	
	if (computeBtotalderivative(bt, dims, &(swKeys_ptr->mean_derivative_btotal), mask, bitmask, derx_bt, dery_bt, bt_err, &(swKeys_ptr->mean_derivative_btotal_err)))
    {
		swKeys_ptr->mean_derivative_btotal = DRMS_MISSING_FLOAT;
		swKeys_ptr->mean_derivative_btotal_err = DRMS_MISSING_FLOAT;
    }
	
	if (computeBhderivative(bh, bh_err, dims, &(swKeys_ptr->mean_derivative_bh), &(swKeys_ptr->mean_derivative_bh_err), mask, bitmask, derx_bh, dery_bh))
    {
		swKeys_ptr->mean_derivative_bh = DRMS_MISSING_FLOAT;
        swKeys_ptr->mean_derivative_bh_err = DRMS_MISSING_FLOAT;
	}
    
	if (computeBzderivative(bz, bz_err, dims, &(swKeys_ptr->mean_derivative_bz), &(swKeys_ptr->mean_derivative_bz_err), mask, bitmask, derx_bz, dery_bz))
    {
		swKeys_ptr->mean_derivative_bz = DRMS_MISSING_FLOAT; // If fail, fill in NaN
        swKeys_ptr->mean_derivative_bz_err = DRMS_MISSING_FLOAT;
    }
	
	computeJz(bx_err, by_err, bx, by, dims, jz, jz_err, jz_err_squared, mask, bitmask, cdelt1, rsun_ref, rsun_obs,
              derx, dery);
    
    
    if(computeJzsmooth(bx, by, dims, jz, jz_smooth, jz_err, jz_rms_err, jz_err_squared_smooth, &(swKeys_ptr->mean_jz),
                       &(swKeys_ptr->mean_jz_err), &(swKeys_ptr->us_i), &(swKeys_ptr->us_i_err), mask, bitmask, cdelt1,
                       rsun_ref, rsun_obs, derx, dery))
    {
        swKeys_ptr->mean_jz            = DRMS_MISSING_FLOAT;
		swKeys_ptr->us_i               = DRMS_MISSING_FLOAT;
        swKeys_ptr->mean_jz_err        = DRMS_MISSING_FLOAT;
        swKeys_ptr->us_i_err           = DRMS_MISSING_FLOAT;
	}
    
	if (computeAlpha(jz_err, bz_err, bz, dims, jz, jz_smooth, &(swKeys_ptr->mean_alpha), &(swKeys_ptr->mean_alpha_err), mask, bitmask, cdelt1, rsun_ref, rsun_obs))
    {
		swKeys_ptr->mean_alpha         = DRMS_MISSING_FLOAT;
        swKeys_ptr->mean_alpha_err     = DRMS_MISSING_FLOAT;
    }
	
	if (computeHelicity(jz_err, jz_rms_err, bz_err, bz, dims, jz, &(swKeys_ptr->mean_ih), &(swKeys_ptr->mean_ih_err), &(swKeys_ptr->total_us_ih), &(swKeys_ptr->total_abs_ih),
                        &(swKeys_ptr->total_us_ih_err), &(swKeys_ptr->total_abs_ih_err), mask, bitmask, cdelt1, rsun_ref, rsun_obs))
    {
		swKeys_ptr->mean_ih            = DRMS_MISSING_FLOAT;
		swKeys_ptr->total_us_ih        = DRMS_MISSING_FLOAT;
  		swKeys_ptr->total_abs_ih       = DRMS_MISSING_FLOAT;
        swKeys_ptr->mean_ih_err        = DRMS_MISSING_FLOAT;
        swKeys_ptr->total_us_ih_err    = DRMS_MISSING_FLOAT;
        swKeys_ptr->total_abs_ih_err   = DRMS_MISSING_FLOAT;
	}
    
	if (computeSumAbsPerPolarity(jz_err, bz_err, bz, jz, dims, &(swKeys_ptr->totaljz), &(swKeys_ptr->totaljz_err),
								 mask, bitmask, cdelt1, rsun_ref, rsun_obs))
    {
		swKeys_ptr->totaljz            = DRMS_MISSING_FLOAT;
        swKeys_ptr->totaljz_err        = DRMS_MISSING_FLOAT;
	}
	
	if (computeFreeEnergy(bx_err, by_err, bx, by, bpx, bpy, dims,
						  &(swKeys_ptr->meanpot), &(swKeys_ptr->meanpot_err), &(swKeys_ptr->totpot), &(swKeys_ptr->totpot_err),
						  mask, bitmask, cdelt1, rsun_ref, rsun_obs))
    {
		swKeys_ptr->meanpot            = DRMS_MISSING_FLOAT; // If fail, fill in NaN
		swKeys_ptr->totpot             = DRMS_MISSING_FLOAT;
        swKeys_ptr->meanpot_err        = DRMS_MISSING_FLOAT;
        swKeys_ptr->totpot_err         = DRMS_MISSING_FLOAT;
	}
    
	if (computeShearAngle(bx_err, by_err, bh_err, bx, by, bz, bpx, bpy, bpz, dims,
						  &(swKeys_ptr->meanshear_angle), &(swKeys_ptr->meanshear_angle_err), &(swKeys_ptr->area_w_shear_gt_45),
						  mask, bitmask)) {
		swKeys_ptr->meanshear_angle    = DRMS_MISSING_FLOAT; // If fail, fill in NaN
		swKeys_ptr->area_w_shear_gt_45 = DRMS_MISSING_FLOAT;
        swKeys_ptr->meanshear_angle_err= DRMS_MISSING_FLOAT;
	}
	
	// Clean up the arrays
	
	drms_free_array(bitmaskArray);		// Dec 18 2012 Xudong
	drms_free_array(maskArray);
	drms_free_array(bxArray);
	drms_free_array(byArray);
	drms_free_array(bzArray);
	
	free(bh); free(bt); free(jz); free(jz_smooth);
	free(bpx); free(bpy); free(bpz);
	free(derx); free(dery);
	free(derx_bt); free(dery_bt);
	free(derx_bz); free(dery_bz);
	free(derx_bh); free(dery_bh);
	free(bt_err); free(bh_err);  free(jz_err);
    free(jz_err_squared); free(jz_rms_err);
    free(jz_err_squared_smooth);
}

/*
 * Set space weather indices, no error checking for now
 *
 */

void setSWIndex(DRMS_Record_t *outRec, struct swIndex *swKeys_ptr)
{
	drms_setkey_float(outRec, "USFLUX",  swKeys_ptr->mean_vf);
	drms_setkey_float(outRec, "MEANGAM", swKeys_ptr->mean_gamma);
	drms_setkey_float(outRec, "MEANGBT", swKeys_ptr->mean_derivative_btotal);
	drms_setkey_float(outRec, "MEANGBH", swKeys_ptr->mean_derivative_bh);
	drms_setkey_float(outRec, "MEANGBZ", swKeys_ptr->mean_derivative_bz);
	drms_setkey_float(outRec, "MEANJZD", swKeys_ptr->mean_jz);
	drms_setkey_float(outRec, "TOTUSJZ", swKeys_ptr->us_i);
	drms_setkey_float(outRec, "MEANALP", swKeys_ptr->mean_alpha);
	drms_setkey_float(outRec, "MEANJZH", swKeys_ptr->mean_ih);
	drms_setkey_float(outRec, "TOTUSJH", swKeys_ptr->total_us_ih);
	drms_setkey_float(outRec, "ABSNJZH", swKeys_ptr->total_abs_ih);
	drms_setkey_float(outRec, "SAVNCPP", swKeys_ptr->totaljz);
	drms_setkey_float(outRec, "MEANPOT", swKeys_ptr->meanpot);
	drms_setkey_float(outRec, "TOTPOT",  swKeys_ptr->totpot);
	drms_setkey_float(outRec, "MEANSHR", swKeys_ptr->meanshear_angle);
	drms_setkey_float(outRec, "SHRGT45", swKeys_ptr->area_w_shear_gt_45);
    drms_setkey_float(outRec, "CMASK",   swKeys_ptr->count_mask);
    drms_setkey_float(outRec, "ERRBT",   swKeys_ptr->mean_derivative_btotal_err);
    drms_setkey_float(outRec, "ERRVF",   swKeys_ptr->mean_vf_err);
    drms_setkey_float(outRec, "ERRGAM",  swKeys_ptr->mean_gamma_err);
    drms_setkey_float(outRec, "ERRBH",   swKeys_ptr->mean_derivative_bh_err);
    drms_setkey_float(outRec, "ERRBZ",   swKeys_ptr->mean_derivative_bz_err);
    drms_setkey_float(outRec, "ERRJZ",   swKeys_ptr->mean_jz_err);
    drms_setkey_float(outRec, "ERRUSI",  swKeys_ptr->us_i_err);
    drms_setkey_float(outRec, "ERRALP",  swKeys_ptr->mean_alpha_err);
    drms_setkey_float(outRec, "ERRMIH",  swKeys_ptr->mean_ih_err);
    drms_setkey_float(outRec, "ERRTUI",  swKeys_ptr->total_us_ih_err);
    drms_setkey_float(outRec, "ERRTAI",  swKeys_ptr->total_abs_ih_err);
    drms_setkey_float(outRec, "ERRJHT",  swKeys_ptr->totaljz_err);
    drms_setkey_float(outRec, "ERRMPOT", swKeys_ptr->meanpot_err);
    drms_setkey_float(outRec, "ERRTPOT", swKeys_ptr->totpot_err);
    drms_setkey_float(outRec, "ERRMSHA", swKeys_ptr->meanshear_angle_err);
};

/*
 * Set all keywords, no error checking for now
 *
 */

void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	copy_me_keys(inRec, outRec);
	copy_patch_keys(inRec, outRec);
	copy_geo_keys(inRec, outRec);
	copy_ambig_keys(inRec, outRec);

    int status = 0;
	
	// Change a few geometry keywords for CEA records
	if (mInfo != NULL) {
        
        drms_setkey_float(outRec, "CRPIX1", mInfo->ncol/2. + 0.5);
		drms_setkey_float(outRec, "CRPIX2", mInfo->nrow/2. + 0.5);
		
		drms_setkey_float(outRec, "CRVAL1", mInfo->xc);
		drms_setkey_float(outRec, "CRVAL2", mInfo->yc);
		drms_setkey_float(outRec, "CDELT1", mInfo->xscale);
		drms_setkey_float(outRec, "CDELT2", mInfo->yscale);
		drms_setkey_string(outRec, "CUNIT1", "degree");
		drms_setkey_string(outRec, "CUNIT2", "degree");
		
		char key[64];
		snprintf (key, 64, "CRLN-%s", wcsCode[(int) mInfo->proj]);
		drms_setkey_string(outRec, "CTYPE1", key);
		snprintf (key, 64, "CRLT-%s", wcsCode[(int) mInfo->proj]);
		drms_setkey_string(outRec, "CTYPE2", key);
		drms_setkey_float(outRec, "CROTA2", 0.0);
		
	} else {
        
        float disk_xc = drms_getkey_float(inRec, "IMCRPIX1", &status);
        float disk_yc = drms_getkey_float(inRec, "IMCRPIX2", &status);
        float x_ll = drms_getkey_float(inRec, "CRPIX1", &status);
        float y_ll = drms_getkey_float(inRec, "CRPIX2", &status);
        // Defined as disk center's pixel address wrt lower-left of cutout
        drms_setkey_float(outRec, "CRPIX1", disk_xc - x_ll + 1.);
		drms_setkey_float(outRec, "CRPIX2", disk_yc - y_ll + 1.);
		// Always 0.
		drms_setkey_float(outRec, "CRVAL1", 0);
		drms_setkey_float(outRec, "CRVAL2", 0);
		
	}
	
	char timebuf[1024];
	float UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
	double val;
	
	val = drms_getkey_double(inRec, "DATE",&status);
	drms_setkey_double(outRec, "DATE_B", val);
	sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
	drms_setkey_string(outRec, "DATE", timebuf);
	
	// set cvs commit version into keyword HEADER
	char *cvsinfo = strdup("$Id: sharp.c,v 1.15 2013/06/27 01:58:39 xudong Exp $");
	//   status = drms_setkey_string(outRec, "HEADER", cvsinfo);
	status = drms_setkey_string(outRec, "CODEVER7", cvsinfo);
	
};

//
//

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

