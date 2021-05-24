/*
 *  smarp.c   
 *
 *	This module creates the pipeline for Space Weather MDI Active Region Patches (SMARPs).
 *	It is a modified version of sharp.c, created by Xudong Sun and Monica Bobra.
 *	It takes the mdi.mtarp series to create the following:
 *
 *      Series 1: mdi.smarp_cea_96m
 *	          CEA remapped magnetogram, bitmap, continuum (same size in map coordinate)
 *                Space weather indices based on line-of-sight magnetogram in the cutout series
 *
 *      Series 2: mdi.smarp_96m
 *	          cutouts of magnetogram, bitmap, continuum, (TARP defined, various sizes in CCD pixels)
 *                Space weather indices based on line-of-sight magnetogram in the cutout series
 *           
 *	Author:
 *		Monica Bobra; Xudong Sun
 *
 *	Version:
 *              v0.0 9 February 2018 Monica Bobra
 *
 *	Notes:
 *		v0.0 No explicit notes
 *		v0.1 Adjusted to take all bitmap values > 36
 *
 *	Example Call:
 *      > smarp "mharp=mdi.Mtarp[13643][2010.10.14_20:48:00_TAI]" "sharp_cea=su_mbobra.smarp_cea_96m" cont="mdi.fd_Ic_interp[2010.10.14_20:48:00_TAI]" "sharp_cut=su_mbobra.smarp_96m"
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
#include "smarp_functions.c"

//#include <mkl.h> // Comment out mkl.h, which can only run on solar3
#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include <omp.h>

#define PI              (M_PI)
#define RADSINDEG	(PI/180.)
#define RAD2ARCSEC	(648000./M_PI)
#define SECINDAY	(86400.)
#define FOURK		(1024)
#define FOURK2          (1048576)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

// FOR HMI: Nyqvist rate at disk center is 0.03 degree. Oversample above 0.015 degree
// FOR HMI: Nyqvist rate at disk center is 0.12 degree. Oversample above 0.06 degree
#define NYQVIST		(0.06)

// Maximum variation of LONDTMAX-LONDTMIN
#define MAXLONDIFF	(1.2e-4)

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

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)
#define XSCALE			0.12
#define YSCALE			0.12
#define NBIN			3
#define INTERP			0
#define dpath    "/home/jsoc/cvs/Development/JSOC"

/* ========================================================================================================== */

// Space weather keywords
struct swIndex {
    float mean_vf;
    float count_mask;
    float absFlux;
    float mean_derivative_los;
    float Rparam;
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
int getInputRS(DRMS_RecordSet_t **mharpRS_ptr, char *mharpQuery);

/* Get other data series */
int getInputRS_aux(DRMS_RecordSet_t **inRS_ptr, char *inQuery, DRMS_RecordSet_t *harpRS);

/* Find record from record set with given T_rec */
int getInputRec_aux(DRMS_Record_t **inRec_ptr, DRMS_RecordSet_t *inRS, TIME trec);

/* Create CEA record */
int createCeaRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *contRec, DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr);

/* Mapping single segment, wrapper */
int mapScaler(DRMS_Record_t *sharpRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec,
			  struct mapInfo *mInfo, char *segName);

/* Determine reference point coordinate and patch size according to input */
int findPosition(DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Compute the coordinates at which the full disk image is sampled */
void findCoord(struct mapInfo *mInfo);

/* Mapping function */
int performSampling(float *outData, float *inData, struct mapInfo *mInfo, int interpOpt);

// ===================

/* Create Cutout record */
int createCutRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *contRec, DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr);

/* Get cutout and write segment */
int writeCutout(DRMS_Record_t *outRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec, char *SegName);

// ===================

/* Compute space weather indices */
void computeSWIndex(struct swIndex *swKeys_ptr, DRMS_Record_t *inRec, struct mapInfo *mInfo);

/* Set space weather indices */
void setSWIndex(DRMS_Record_t *outRec, struct swIndex *swKeys_ptr);

/* Set all keywords */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *mharpRec, struct mapInfo *mInfo);

// ===================

/* Nearest neighbor interpolation */
float nnb (float *f, int nx, int ny, double x, double y);

/* Wrapper for Jesper's rebin code */
void frebin (float *image_in, float *image_out, int nx, int ny, int nbin, int gauss);

/* ========================================================================================================== */

/* Cutout segment names, input identical to output */
char *MharpSegs[] = {"magnetogram", "bitmap"};
char *CutSegs[] = {"magnetogram", "bitmap", "continuum"};
char *CEASegs[] = {"magnetogram", "bitmap", "continuum"};
// For BUNIT
char *CutBunits[] = {"Mx/cm^2", " ", "DN/s"};
char *CEABunits[] = {"Mx/cm^2", " ", "DN/s"};
/* ========================================================================================================== */

char *module_name = "smarp";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "mharp", kNotSpecified, "Input Mharp series."},
    {ARG_STRING, "cont", kNotSpecified, "Input Continuum series."},
    {ARG_STRING, "sharp_cea", kNotSpecified, "Output Sharp CEA series."},
    {ARG_STRING, "sharp_cut", kNotSpecified, "Output Sharp cutout series."},
    {ARG_END}
};

int DoIt(void)
{
        int errbufstat = setvbuf(stderr, NULL, _IONBF, BUFSIZ);
        int outbufstat = setvbuf(stdout, NULL, _IONBF, BUFSIZ);
	int status = DRMS_SUCCESS;
	int nrecs, irec;
	char *mharpQuery; 
        char *contQuery;
	char *sharpCeaQuery, *sharpCutQuery;
	DRMS_RecordSet_t *mharpRS = NULL;
	DRMS_RecordSet_t *contRS = NULL;

	/* Get parameters */
    
	mharpQuery = (char *) params_get_str(&cmdparams, "mharp");
	sharpCeaQuery = (char *) params_get_str(&cmdparams, "sharp_cea");
	sharpCutQuery = (char *) params_get_str(&cmdparams, "sharp_cut");
        contQuery = (char *) params_get_str(&cmdparams, "cont");
	
	/* Get input data, check everything */
        if (getInputRS(&mharpRS, mharpQuery))
            DIE("Input harp data error.");
	    nrecs = mharpRS->n;

	if (getInputRS_aux(&contRS, contQuery, mharpRS))
	    DIE("Input continuum data error.");	

	/* Start */
	
	printf("==============\nStart. %d image(s) in total.\n", nrecs);
    
	for (irec = 0; irec < nrecs; irec++) {
		
		/* Records in work */
		
		DRMS_Record_t *mharpRec = NULL;
		DRMS_Record_t *contRec = NULL;

		mharpRec = mharpRS->records[irec];
                TIME trec = drms_getkey_time(mharpRec, "T_REC", &status);
        
		struct swIndex swKeys;

		if (getInputRec_aux(&contRec, contRS, trec)) {
			printf("Fetching Continuum failed, image #%d skipped.\n", irec);
			continue;
		}

	        printf("Obtained all the data \n");
        
		/* Create CEA record */

		DRMS_Record_t *sharpCeaRec = drms_create_record(drms_env, sharpCeaQuery, DRMS_PERMANENT, &status);
		if (status) {		// if failed
			printf("Creating CEA failed, image #%d skipped.\n", irec);
			continue;
		}
		if (createCeaRecord(mharpRec, contRec, sharpCeaRec, &swKeys)) {		// do the work
			printf("Creating CEA failed, image #%d skipped.\n", irec);
			drms_close_record(sharpCeaRec, DRMS_FREE_RECORD);
			continue;
		}		// swKeys updated here
		
		drms_close_record(sharpCeaRec, DRMS_INSERT_RECORD);

	        printf("Created CEA record \n");
				
		/* Create Cutout record */
		
		DRMS_Record_t *sharpCutRec = drms_create_record(drms_env, sharpCutQuery, DRMS_PERMANENT, &status);
		if (status) {		// if failed
			printf("Creating cutout failed, image #%d skipped.\n", irec);
			continue;
		}
		
		if (createCutRecord(mharpRec, contRec, sharpCutRec, &swKeys)) {		// do the work
			printf("Creating cutout failed, image #%d skipped.\n", irec);
			drms_close_record(sharpCutRec, DRMS_FREE_RECORD);
			continue;
		}		// swKeys used here
		drms_close_record(sharpCutRec, DRMS_INSERT_RECORD);
	        printf("Created CUT record \n");
		/* Done */
		
		printf("Image #%d done.\n", irec);
		
	} // irec
    
	
	drms_close_records(mharpRS, DRMS_FREE_RECORD);
	drms_close_records(contRS, DRMS_FREE_RECORD);
	
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

int getInputRS(DRMS_RecordSet_t **mharpRS_ptr, char *mharpQuery)
{
	int status = 0;
	*mharpRS_ptr = drms_open_records(drms_env, mharpQuery, &status);
        if (status || (*mharpRS_ptr)->n == 0) return 1;      
	return 0;	
}

/*
 * Get other data series, check all T_REC are available
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
 */

int createCeaRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *contRec, DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr)
{
	
	int status = 0;
	DRMS_Segment_t *inSeg;
	DRMS_Array_t *inArray;
	int val;

	struct mapInfo mInfo;
	mInfo.proj = (enum projection) cyleqa;		// projection method
	mInfo.xscale = XSCALE;
	mInfo.yscale = YSCALE;
	
    int ncol0, nrow0;		// oversampled map size
	
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
	
	// ========================================
	// Do this for all bitmaps, Aug 12 2013 XS
	// ========================================
	
        mInfo.nbin = 1;			// for bitmaps. suppress anti-aliasing
	ncol0 = mInfo.ncol;
	nrow0 = mInfo.nrow;
	
	mInfo.xi_out = (float *) (malloc(ncol0 * nrow0 * sizeof(float)));
	mInfo.zeta_out = (float *) (malloc(ncol0 * nrow0 * sizeof(float)));
	
	findCoord(&mInfo);		// compute it here so it could be shared by the following 4 functions
	
	if (mapScaler(sharpRec, mharpRec, mharpRec, &mInfo, "bitmap")) {
		SHOW("CEA: mapping bitmap error\n");
		return 1;
	}
	printf("Bitmap mapping done.\n");
	
        free(mInfo.xi_out);
	free(mInfo.zeta_out);

	// ========================================
	// Do this again for floats, Aug 12 2013 XS
	// ========================================
	// Create xi_out, zeta_out array in mInfo:
	// Coordinates to sample in original full disk image
	
	mInfo.nbin = NBIN;
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
	printf("Magnetogram mapping done.\n");
	 
	if (mapScaler(sharpRec, contRec, mharpRec, &mInfo, "continuum")) {
		SHOW("CEA: mapping continuum error\n");
		return 1;
	}
	printf("Intensitygram mapping done.\n");

	// Keywords & Links
	copy_patch_keys(mharpRec, sharpRec);
	copy_geo_keys(mharpRec, sharpRec);
        // rename HARPNUM to TARPNUM
	val = drms_getkey_double(mharpRec, "HARPNUM", &status);
        drms_setkey_double(sharpRec, "TARPNUM", val);	
	// copy everything else 
	drms_copykey(sharpRec, mharpRec, "T_REC");
	drms_copykey(sharpRec, mharpRec, "CDELT1");
	drms_copykey(sharpRec, mharpRec, "RSUN_OBS");
	drms_copykey(sharpRec, mharpRec, "DSUN_OBS");
	drms_copykey(sharpRec, mharpRec, "OBS_VR");
	drms_copykey(sharpRec, mharpRec, "OBS_VW");
	drms_copykey(sharpRec, mharpRec, "OBS_VN");
        drms_copykey(sharpRec, mharpRec, "CRLN_OBS");
        drms_copykey(sharpRec, mharpRec, "CRLT_OBS");
	drms_copykey(sharpRec, mharpRec, "CAR_ROT");
	drms_copykey(sharpRec, mharpRec, "SIZE_SPT");
	drms_copykey(sharpRec, mharpRec, "AREA_SPT");
        drms_copykey(sharpRec, mharpRec, "DATE__OBS");
        drms_copykey(sharpRec, mharpRec, "T_OBS");
        drms_copykey(sharpRec, mharpRec, "T_MAXPIX");
	drms_copykey(sharpRec, mharpRec, "QUALITY");
	drms_copykey(sharpRec, mharpRec, "NPIX_SPT");
	drms_copykey(sharpRec, mharpRec, "ARS_NCLN");
	drms_copykey(sharpRec, mharpRec, "ARS_MODL");
	drms_copykey(sharpRec, mharpRec, "ARS_EDGE");
	drms_copykey(sharpRec, mharpRec, "ARS_BETA");
	drms_copykey(sharpRec, mharpRec, "T_MID1");
	drms_copykey(sharpRec, mharpRec, "T_CMPASS");

	DRMS_Link_t *mHarpLink = hcon_lookup_lower(&sharpRec->links, "MTARP");
	if (mHarpLink) {
            drms_link_set("MTARP", sharpRec, mharpRec);
        }

        // set other keywords
        setKeys(sharpRec, mharpRec, &mInfo);

	// set space weather keywords
        computeSWIndex(swKeys_ptr, sharpRec, &mInfo);	 
        printf("Space weather indices done.\n");
	    setSWIndex(sharpRec, swKeys_ptr);
	         
	// set statistical keywords (e.g. DATAMIN, DATAMAX, etc.)	
	//int nCEASegs = ARRLENGTH(CEASegs);
	int nCEASegs = 3;	
        for (int iSeg = 0; iSeg < 3; iSeg++) {
		DRMS_Segment_t *outSeg = drms_segment_lookupnum(sharpRec, iSeg);
		DRMS_Array_t *outArray = drms_segment_read(outSeg, DRMS_TYPE_FLOAT, &status);
		int stat = set_statistics(outSeg, outArray, 1);
		//printf("%d => %d\n", iSeg, stat);
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
 */

int mapScaler(DRMS_Record_t *sharpRec, DRMS_Record_t *inRec, DRMS_Record_t *harpRec,
			  struct mapInfo *mInfo, char *segName)
{
	
	int status = 0;
	int nx = mInfo->ncol, ny = mInfo->nrow, nxny = nx * ny;
	int dims[2] = {nx, ny};
	int interpOpt = INTERP;		// Aug 12 XS, default, overridden below for bitmaps and conf_disambig
	
	// Input full disk array
	
	DRMS_Segment_t *inSeg = NULL;
	inSeg = drms_segment_lookup(inRec, segName);
	if (!inSeg) return 1;
	
	DRMS_Array_t *inArray = NULL;
	inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
	if (!inArray) return 1;

    if (!strcmp(segName, "conf_disambig") || !strcmp(segName, "bitmap")) {
        // Moved out so it works for FD conf_disambig as well
        // Jan 2 2014 XS
        interpOpt = 3;		// Aug 12 XS, near neighbor
    }
	    
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
	if (performSampling(map, inData, mInfo, interpOpt))		// Add interpOpt for different types, Aug 12 XS
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
 * Determine reference point coordinate and patch size according to keywords
 * xc, yc are the coordinate of patch center, in degrees
 * ncol and nrow are the final size
 */

int findPosition(DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	int harpnum = drms_getkey_int(inRec, "TARPNUM", &status);
	TIME trec = drms_getkey_time(inRec, "T_REC", &status);
	float disk_lonc = drms_getkey_float(inRec, "CRLN_OBS", &status);
	
	/* Center coord */
    // Changed into double Jun 16 2014 XS
	
	double minlon = drms_getkey_double(inRec, "LONDTMIN", &status); if (status) return 1;		// Stonyhurst lon
	double maxlon = drms_getkey_double(inRec, "LONDTMAX", &status); if (status) return 1;
	double minlat = drms_getkey_double(inRec, "LATDTMIN", &status); if (status) return 1;
	double maxlat = drms_getkey_double(inRec, "LATDTMAX", &status); if (status) return 1;
	
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
    // Rounded to 1.d3 precision first. Jun 16 2014 XS
    // The previous fix does not work. LONDTMAX-LONDTMIN varies from frame to frame
    // Need to find out the maximum possible difference, MAXLONDIFF (1.2e-4)
    // Now, ncol = (maxlon-minlon)/xscale, if the decimal part is outside 0.5 \pm (MAXLONDIFF/xscale)
    // proceed as it is. else, all use floor on ncol
	
	float dpix = (MAXLONDIFF / mInfo->xscale) * 1.5;		// "danger zone"
	float ncol = (maxlon - minlon) / mInfo->xscale;
	float d_ncol = fabs(ncol - floor(ncol) - 0.5);			// distance to 0.5
	if (d_ncol < dpix) {
		mInfo->ncol = floor(ncol);
	} else {
		mInfo->ncol = round(ncol);
	}

	mInfo->nrow = round((maxlat - minlat) / mInfo->yscale);
	
	return 0;
	
}


/*
 * Fetch ephemeris info from a DRMS record
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
 */

int performSampling(float *outData, float *inData, struct mapInfo *mInfo, int interpOpt)
{
	
	int status = 0;
	int ind_map;
	
	int ncol0 = mInfo->ncol * mInfo->nbin + (mInfo->nbin / 2) * 2;	// pad with nbin/2 on edge to avoid NAN
	int nrow0 = mInfo->nrow * mInfo->nbin + (mInfo->nbin / 2) * 2;
	
	// Changed Aug 12 2013, XS, for bitmaps
	float *outData0;
	if (interpOpt == 3 && mInfo->nbin == 1) {
        outData0 = outData;
	} else {
        outData0 = (float *) (malloc(ncol0 * nrow0 * sizeof(float)));
	}
	
	float *xi_out = mInfo->xi_out;
	float *zeta_out = mInfo->zeta_out;
	
	// Interpolation
	
	struct fint_struct pars;
	// Aug 12 2013, passed in as argument now
	
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
		case 3:			// Near neighbor
            break;
		default:
			return 1;
	}
	
	//printf("interpOpt = %d, nbin = %d ", interpOpt, mInfo->nbin);
	if (interpOpt == 3) {			// Aug 6 2013, Xudong
	  	for (int row0 = 0; row0 < nrow0; row0++) {
            for (int col0 = 0; col0 < ncol0; col0++) {
                ind_map = row0 * ncol0 + col0;
                outData0[ind_map] = nnb(inData, FOURK, FOURK, xi_out[ind_map], zeta_out[ind_map]);
            }
        }
	} else {
        finterpolate(&pars, inData, xi_out, zeta_out, outData0,
                     FOURK, FOURK, FOURK, ncol0, nrow0, ncol0, DRMS_MISSING_FLOAT);
	}
	
	// Rebinning, smoothing
	
	if (interpOpt == 3 && mInfo->nbin == 1) {
        return 0;
	} else {
        frebin(outData0, outData, ncol0, nrow0, mInfo->nbin, 1);		// Gaussian
        free(outData0);		
	}

	return 0;
	
}

/*
 * Create Cutout record: top level subroutine
 * Do the loops on segments and set the keywords here
 * Work is done in writeCutout routine below
 */

int createCutRecord(DRMS_Record_t *mharpRec, DRMS_Record_t *contRec, DRMS_Record_t *sharpRec, struct swIndex *swKeys_ptr)
{
	
	int status = 0;
	int val;	
	int iHarpSeg;
	int nMharpSegs = ARRLENGTH(MharpSegs);
	
	// Cutout Mharp
	
	for (iHarpSeg = 0; iHarpSeg < nMharpSegs; iHarpSeg++) {
		if (writeCutout(sharpRec, mharpRec, mharpRec, MharpSegs[iHarpSeg])) {
			printf("Mharp cutout fails for %s\n", MharpSegs[iHarpSeg]);
			printf("iHarpSeg nMharpSegs %d %d \n",iHarpSeg,nMharpSegs);
			break;
		}
	}
	if (iHarpSeg != nMharpSegs) {
		SHOW("Cutout: segment number mismatch\n");
		return 1;		// if failed
	}
	printf("Magnetogram cutout done.\n");

	// Cutout Continuum
	
	if (writeCutout(sharpRec, contRec, mharpRec, "continuum")) {
		printf("Continuum cutout failed\n");
		return 1;
	}
	printf("Intensitygram cutout done.\n");
			
	// Keywords & Links
	copy_patch_keys(mharpRec, sharpRec);
	copy_geo_keys(mharpRec, sharpRec);
        // rename HARPNUM to TARPNUM
	val = drms_getkey_double(mharpRec, "HARPNUM", &status);
        drms_setkey_double(sharpRec, "TARPNUM", val);	
	// copy everything else 
	drms_copykey(sharpRec, mharpRec, "T_REC");
	drms_copykey(sharpRec, mharpRec, "CDELT1");
	drms_copykey(sharpRec, mharpRec, "RSUN_OBS");
	drms_copykey(sharpRec, mharpRec, "DSUN_OBS");
	drms_copykey(sharpRec, mharpRec, "OBS_VR");
	drms_copykey(sharpRec, mharpRec, "OBS_VW");
	drms_copykey(sharpRec, mharpRec, "OBS_VN");
        drms_copykey(sharpRec, mharpRec, "CRLN_OBS");
        drms_copykey(sharpRec, mharpRec, "CRLT_OBS");
	drms_copykey(sharpRec, mharpRec, "CAR_ROT");
	drms_copykey(sharpRec, mharpRec, "SIZE_SPT");
	drms_copykey(sharpRec, mharpRec, "AREA_SPT");
        drms_copykey(sharpRec, mharpRec, "DATE__OBS");
        drms_copykey(sharpRec, mharpRec, "T_OBS");
        drms_copykey(sharpRec, mharpRec, "T_MAXPIX");
	drms_copykey(sharpRec, mharpRec, "QUALITY");
	drms_copykey(sharpRec, mharpRec, "NPIX_SPT");
	drms_copykey(sharpRec, mharpRec, "ARS_NCLN");
	drms_copykey(sharpRec, mharpRec, "ARS_MODL");
	drms_copykey(sharpRec, mharpRec, "ARS_EDGE");
	drms_copykey(sharpRec, mharpRec, "ARS_BETA");
	drms_copykey(sharpRec, mharpRec, "T_MID1");
	drms_copykey(sharpRec, mharpRec, "T_CMPASS");

	DRMS_Link_t *mHarpLink = hcon_lookup_lower(&sharpRec->links, "MTARP");
	if (mHarpLink) {
            drms_link_set("MTARP", sharpRec, mharpRec);
        }
	
	setSWIndex(sharpRec, swKeys_ptr);	// Set space weather indices
	setKeys(sharpRec, mharpRec, NULL);      // Set all other keywords, NULL specifies cutout
	
	// Stats
       	
	int nCutSegs = 3; 
	for (int iSeg = 0; iSeg < 3; iSeg++) {
		DRMS_Segment_t *outSeg = drms_segment_lookupnum(sharpRec, iSeg);
		DRMS_Array_t *outArray = drms_segment_read(outSeg, DRMS_TYPE_FLOAT, &status);
		set_statistics(outSeg, outArray, 1);
		drms_free_array(outArray);
	}
	
	return 0;
	
}

/*
 * Get cutout and write segment
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
	//printf("SegName=%s\n",SegName); fflush(stdout);
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
	//printf("line1068\n"); fflush(stdout);
	return 0;
	
}


/*
 * Compute space weather indices
 */

void computeSWIndex(struct swIndex *swKeys_ptr, DRMS_Record_t *inRec, struct mapInfo *mInfo)
{
	
	int status = 0;
	int nx = mInfo->ncol, ny = mInfo->nrow;
	int nxny = nx * ny;
	int dims[2] = {nx, ny};
	
	// Get mask + use TARP bitmap as a threshold on spaceweather quantities
	DRMS_Segment_t *bitmaskSeg = drms_segment_lookup(inRec, "bitmap");
	DRMS_Array_t *bitmaskArray = drms_segment_read(bitmaskSeg, DRMS_TYPE_INT, &status);
	int *bitmask = (int *) bitmaskArray->data;		// bitmap
	    
    // Get line-of-sight magnetogram
    DRMS_Segment_t *losSeg = drms_segment_lookup(inRec, "magnetogram");
    DRMS_Array_t *losArray = drms_segment_read(losSeg, DRMS_TYPE_FLOAT, &status);
    float *los = (float *) losArray->data;          // los
    	 
	// Get emphemeris
	float  cdelt1_orig = drms_getkey_float(inRec, "CDELT1",   &status);
	float  dsun_obs    = drms_getkey_float(inRec, "DSUN_OBS",   &status);
	double rsun_ref    = drms_getkey_double(inRec, "RSUN_REF", &status);
	double rsun_obs    = drms_getkey_double(inRec, "RSUN_OBS", &status);
	float imcrpix1     = drms_getkey_float(inRec, "IMCRPIX1", &status);
	float imcrpix2     = drms_getkey_float(inRec, "IMCRPIX2", &status);
	float crpix1       = drms_getkey_float(inRec, "CRPIX1", &status);
	float crpix2       = drms_getkey_float(inRec, "CRPIX2", &status);
    
    // convert cdelt1_orig from degrees to arcsec
    float cdelt1       = (atan((rsun_ref*cdelt1_orig*RADSINDEG)/(dsun_obs)))*(1/RADSINDEG)*(3600.);

	// Temp arrays
	float *derx_los = (float *) (malloc(nxny * sizeof(float)));
	float *dery_los = (float *) (malloc(nxny * sizeof(float)));
     
    // define some values for the R calculation
    int scale = round(2.0/cdelt1);
    int nx1 = nx/scale;
    int ny1 = ny/scale;
    int nxp = nx1+40; // same comment as above
    int nyp = ny1+40; // why is this a +40 pixel size? is this an MDI pixel?
    float *rim     = (float *)malloc(nx1*ny1*sizeof(float));
    float *p1p0    = (float *)malloc(nx1*ny1*sizeof(float));
    float *p1n0    = (float *)malloc(nx1*ny1*sizeof(float));
    float *p1p     = (float *)malloc(nx1*ny1*sizeof(float));
    float *p1n     = (float *)malloc(nx1*ny1*sizeof(float));
    float *p1      = (float *)malloc(nx1*ny1*sizeof(float));
    float *pmap    = (float *)malloc(nxp*nyp*sizeof(float));
    float *p1pad   = (float *)malloc(nxp*nyp*sizeof(float));
    float *pmapn   = (float *)malloc(nx1*ny1*sizeof(float));
    
	// Compute three spaceweather quantities, USFLUX, MEANGBZ, R_VALUE, on LOS data
	
	if (computeAbsFlux_los(los, dims, &(swKeys_ptr->absFlux), &(swKeys_ptr->mean_vf),
                           &(swKeys_ptr->count_mask), bitmask, cdelt1, rsun_ref, rsun_obs))
        {
		swKeys_ptr->absFlux = DRMS_MISSING_FLOAT;		// If fail, fill in NaN
		swKeys_ptr->mean_vf = DRMS_MISSING_FLOAT;
        swKeys_ptr->count_mask  = DRMS_MISSING_INT;
	    }
    
	if (computeLOSderivative(los, dims, &(swKeys_ptr->mean_derivative_los), 
                                bitmask, derx_los, dery_los))
        {
		swKeys_ptr->mean_derivative_los = DRMS_MISSING_FLOAT; // If fail, fill in NaN
        }
	
	if (computeR_los(los, dims, &(swKeys_ptr->Rparam), cdelt1, rim, p1p0, p1n0,
                     p1p, p1n, p1, pmap, nx1, ny1, scale, p1pad, nxp, nyp, pmapn))
        {
		swKeys_ptr->Rparam = DRMS_MISSING_FLOAT;		// If fail, fill in NaN
        }
    	
	// Clean up the arrays
      	drms_free_array(bitmaskArray);
        drms_free_array(losArray);
        // free arrays related to LOS derivative
	    free(derx_los); free(dery_los);
        // free the arrays that are related to the r calculation     
        free(rim);
        free(p1p0);
        free(p1n0);
        free(p1p);
        free(p1n);
        free(p1);
        free(pmap);
        free(p1pad);
        free(pmapn);
}

/*
 * Set space weather indices
 */

void setSWIndex(DRMS_Record_t *outRec, struct swIndex *swKeys_ptr)
{
    drms_setkey_float(outRec, "USFLUXL",  swKeys_ptr->mean_vf);
    drms_setkey_float(outRec, "MEANGBL", swKeys_ptr->mean_derivative_los);
    drms_setkey_float(outRec, "R_VALUE", swKeys_ptr->Rparam);
    drms_setkey_float(outRec, "CMASKL", swKeys_ptr->count_mask);
};

/*
 * Set all keywords, no error checking for now
 */

void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *mharpRec, struct mapInfo *mInfo)
{
    
        int status = 0;
	
	// Change a few geometry keywords for CEA & cutout records
	if (mInfo != NULL) {        // CEA
	  printf("Calculating CEA keys\n");
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
            	// Set BUNIT for each segment        
        	int nSeg = 3;
        	for (int iSeg = 0; iSeg < nSeg; iSeg++) {
            	DRMS_Segment_t *outSeg = NULL;
            	outSeg = drms_segment_lookup(outRec, CEASegs[iSeg]);
            	if (!outSeg) continue;
            	char bunit_xxx[20];
            	sprintf(bunit_xxx, "BUNIT_%03d", iSeg);
            	//printf("%s, %s\n", bunit_xxx, CEABunits[iSeg]);
            	drms_setkey_string(outRec, bunit_xxx, CEABunits[iSeg]);
        		}
		
	} else {        // Cutout
        
        	float disk_xc, disk_yc;
        	disk_xc = drms_getkey_float(mharpRec, "IMCRPIX1", &status);
       		disk_yc = drms_getkey_float(mharpRec, "IMCRPIX2", &status);
        	float x_ll = drms_getkey_float(mharpRec, "CRPIX1", &status);
        	float y_ll = drms_getkey_float(mharpRec, "CRPIX2", &status);
        	// Defined as disk center's pixel address wrt lower-left of cutout
        	drms_setkey_float(outRec, "CRPIX1", disk_xc - x_ll + 1.);
		drms_setkey_float(outRec, "CRPIX2", disk_yc - y_ll + 1.);
		// Always 0.
		drms_setkey_float(outRec, "CRVAL1", 0);
		drms_setkey_float(outRec, "CRVAL2", 0);
        
        	// Jan 2 2014 XS
        	int nSeg = ARRLENGTH(CutSegs);
        	for (int iSeg = 0; iSeg < nSeg; iSeg++) 
                       {
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
	
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
	
    val = drms_getkey_time(mharpRec, "DATE", &status);
    drms_setkey_time(outRec, "DATE", tnow);
	
    // set cvs commit version into keyword CODEVER7
    char *cvsinfo  = strdup("$Id");
    char *cvsinfo2 = smarp_functions_version();
    char cvsinfoall[2048];
    strcat(cvsinfoall,cvsinfo);
    strcat(cvsinfoall,"\n");
    strcat(cvsinfoall,cvsinfo2);
    status = drms_setkey_string(outRec, "CODEVER7", cvsinfoall);

};

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
