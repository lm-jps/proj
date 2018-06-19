/*
 *  cgem_prep.c
 *
 *  Example:
 %  cgem_prep "in=hmi_test.cgem_cutout_720s[11158][2011.02.15_12:00/252m]" "out=hmi_test.cgem_prep" "tref=2011.02.15_12:00:00_TAI" "cols=600" "rows=600" "lonref=18.6" "latref=-20.4" "dx=0.03" "dy=0.03"
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
#include "img2helioVector.c"
#include "finterpolate.h"
#include "diffrot.h"            // double diffrot[3]

// Legacy macros

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)
#define DTTHRESH        (SECINDAY*7.)					// 7 day
#define INTERPOPT       (0)         // Interpolation

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"
#define dpath    "/home/jsoc/cvs/Development/JSOC"

#define vSegName	"vlos_mag"

#define NT          (5)             // chunk size
#define FLIP_THR    (120)           // azimuth flipping threshold

//#define TEST    (1)                 // Test output corrected azi/dopp in cutout
//#define NOAZI   (1)                 // No azimuth correction
//#define NODOP   (1)                 // No Doppler correction

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

// WSC code
char *wcsCode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
    "SIN", "ZEA"};

// Output segments
char *outSegNames[] = {"Bx", "By", "Bz", "Lx", "Ly", "Lz", vSegName};

// Ephemeris information
struct ephemeris {
	double disk_lonc, disk_latc;
	double disk_xc, disk_yc;
    double rSun, asd, pa;
    TIME t_rec;
};

// Requested mapping information
struct reqInfo {
    double lonref, latref;
    int ncol, nrow;
    double dx, dy;
    TIME tref;
    int proj;
};

// =====================================

/* Read named array from record */
int readDRMSArrayF(DRMS_Record_t *inRec, float *data, const char *segName);     // float version
int readDRMSArrayD(DRMS_Record_t *inRec, double *data, const char *segName);    // double version

/* Correct azimuth */
int correct_azi(DRMS_RecordSet_t *inRS, int irec, int nt, int *dims,
                float *azi_work, float *azi_corr);

/* Correct Dopplergram */
int correct_dop(DRMS_Record_t *inRec, int *dims, float *dop_corr);

/* Compute Dopper correction due to spacecraft and rotation */
int getDoppCorr(DRMS_Record_t *inRec, int *dims, double *doppCorr);

/* Mapping field and LOS vectors and Doppler */
int cgem_mapping(DRMS_Record_t *inRec, struct reqInfo *req, int *dims_in, float *mapCenter,
                 float *azi_corr, float *dop_corr,
                 float *bx_map, float *by_map, float *bz_map,
                 float *lx_map, float *ly_map, float *lz_map, float *dop_map);

/* Output */
int writeOutput(DRMS_Record_t *inRec, DRMS_Record_t *outRec, struct reqInfo *req, float *mapCenter,
                float *bx_map, float *by_map, float *bz_map,
                float *lx_map, float *ly_map, float *lz_map, float *dop_map);

/* Output cutout for test */
int writeOutputTest(DRMS_Record_t *inRec, DRMS_Record_t *outRec, int *dims,
                    float *azi_corr, float *dop_corr);

/* Get ephemeris information */
int getEphemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

/* Compute the pixel coordinates at which the map is sampled */
void findCoord(struct reqInfo *req, struct ephemeris *ephem,
               float *u_out, float *v_out, float *mapCenter);

/* Performing local vector transformation */
void vectorTransform(float *fld, float *inc, float *az, float *bx, float *by, float *bz,
                     int *dims, struct ephemeris *ephem);

/* Convert LOS vectors to local xyz coordinate */
void getLOSVector(struct reqInfo *req, struct ephemeris *ephem, float *mapCenter,
                  float *lx, float *ly, float *lz);

/* Mapping function */
void performSampling(float *inData, float *outData, float *u, float *v, int *dims_in, int *dims_out);

/* Nearest neighbor interpolation */
float nnb (float *f, int nx, int ny, float x, float y);

/* Azimuth correction */
extern void azim_mindiff_jsoc_(float *azi_work, float *azi_out, int *nx, int *ny, int *nt2, int *flip_thr);

// =====================================

char *module_name = "cgem_prep";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input data series."},
    {ARG_STRING, "out", "hmi_test.cgem_prep", "Output data series."},
    {ARG_NUME, "map", "carree", "Projetion method, carree by default.",
        "carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
    {ARG_DOUBLE, "lonref", "0", "Reference patch center Stonyhurst lon, in deg."},
	{ARG_DOUBLE, "latref", "0", "Reference patch center Stonyhurst lat, in deg."},
	{ARG_INT, "cols", "480", "Columns of output cutout."},
	{ARG_INT, "rows", "480", "Rows of output cutout."},
	{ARG_DOUBLE, "dx", "0.03", "X pixel size, unit depending on projection (default deg)."},
	{ARG_FLOAT, "dy", "0.03", "Y pixel size, unit depending on projection (default deg)."},
	{ARG_STRING, "tref", kNotSpecified, "Reference time."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;

    /* Get data series */
    
    char *inQuery, *outQuery;    // Data series query
    
    inQuery = (char *) params_get_str(&cmdparams, "in");
    outQuery = (char *) params_get_str(&cmdparams, "out");
    
    /* Arguments */
    
    int nt = round(fabs(NT));           // Number of records for batch azimuth correction
    nt = (nt % 2) ? nt : nt + 1;        // odd
    int nt2 = nt / 2;
    
    struct reqInfo req;
    req.lonref = params_get_double(&cmdparams, "lonref") * RADSINDEG;
    req.latref = params_get_double(&cmdparams, "latref") * RADSINDEG;
    req.ncol = params_get_int(&cmdparams, "cols");
    req.nrow = params_get_int(&cmdparams, "rows");
    req.tref = params_get_time(&cmdparams, "tref");
    req.dx = params_get_double(&cmdparams, "dx") * RADSINDEG;		// deg to rad, for now
    req.dy = params_get_double(&cmdparams, "dy") * RADSINDEG;
    req.proj = params_get_int(&cmdparams, "map");
    
    /* Input data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs == 0 || !inRS) DIE("Input data series error");
#ifndef NOAZI
    if (nrecs < nt) DIE("Input data series not long enough for azi correction");
#endif
    
    // Check reference time, if more than 6 days from first record, use first record
    
    TIME t0 = drms_getkey_time(inRS->records[0], "T_REC", &status);
    if (fabs(req.tref - t0) > DTTHRESH) {
        req.tref = t0;
    }
    
    // All cutouts must have same size
    
    int nx, ny;
    for (int irec = 0; irec < nrecs; irec++) {
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, "azimuth");
        if (!inSeg) DIE("Can't check array dimension!");
        if (irec == 0) {
            nx = inSeg->axis[0]; ny = inSeg->axis[1];
            continue;
        }
        if ((inSeg->axis[0] != nx) || (inSeg->axis[1] != ny)) {
            DIE("Array sizes do not match!");
        }
    }
    int nxny = nx * ny, dims[2] = {nx, ny};
    
    printf("\n=============\nStart processing, %d in total.\n", nrecs);
    
    // Pre-load (nt - 1) azimuth into the working array for correction
    // Correction starts at frame irec = nt2, frames #0, #1, ..., #nt-2 already loaded
    // No move all frames in working array forward by 1
    // read in frame irec + nt2 to the end of the work array
    
    printf("=============\nPre-loading %d azimuth maps: ", nt - 1);
    
    float *azi_work = (float *) (malloc(nxny * nt * sizeof(float)));    // working array for azimuth correction
    float *azi_work_lead;           // location in azi_work
    
    for (int irec = 0; irec < nt - 1; irec++) {
        
        DRMS_Record_t *inRec = inRS->records[irec];
        azi_work_lead = azi_work + nxny * (irec + 1);
        if (readDRMSArrayF(inRec, azi_work_lead, "azimuth")) {
            free(azi_work);
            DIE("Can't retrieve azimtuh!");
        }
        printf("%d ", irec + 1);
        
    }
    
    printf("done.\n");
    
    /* Output data */
    
     DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
     if (status || !outRS) {DIE("Error in output series");}
    
    /* Main loop */
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        printf("=============\nProcessing frame %d of %d:", irec + 1, nrecs);
        
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *outRec = outRS->records[irec];
        
        int cgemnum = drms_getkey_int(inRec, "CGEMNUM", &status);
        TIME t_rec = drms_getkey_time(inRec, "T_REC", &status);
        char t_rec_str[100];
        sprint_time(t_rec_str, t_rec, "TAI", 0);
        
        printf(" [%d][%s]\n", cgemnum, t_rec_str);
        
        /* ================== */
        /* Correcting azimuth */
        /* ================== */
        
        printf("correcting azimuth: ");
        
        float *azi_corr = (float *) (malloc(nxny * sizeof(float)));     // Target cutout array
        
        if (correct_azi(inRS, irec, nt, dims, azi_work, azi_corr)) {    // work
            free(azi_corr);
            printf("error, frame skipped.\n");
            continue;
        }
        
        printf("done.\n");
        
        /* ================== */
        /* Correcting Doppler */
        /* ================== */
        
        printf("correcting Doppler: ");
        
        float *dop_corr = (float *) (malloc(nxny * sizeof(float)));  // Target cutout array
        
        if (correct_dop(inRec, dims, dop_corr)) {                       // work
            free(azi_corr); free(dop_corr);
            printf("error, frame skipped.\n");
            continue;
        }
        
        printf("done.\n");
        
#ifdef TEST
        
        /* ========================= */
        /* Output cutout for testing */
        /* ========================= */
        
        printf("output cutout for testing: ");
        
        if (writeOutputTest(inRec, outRec, dims, azi_corr, dop_corr)) {   // Clean up done inside subroutine
            printf("error, frame skipped.\n");
            continue;
        }
        
        printf("done.\n");
        
#else
        
        /* ============================= */
        /* Mapping field and LOS vectors */
        /* ============================= */
        
        printf("mapping: ");
        
        float *bx_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float *by_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float *bz_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float *lx_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float *ly_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float *lz_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float *dop_map = (float *) (malloc(req.ncol * req.nrow * sizeof(float)));
        float mapCenter[2];        // map center coordinate
        
        if (cgem_mapping(inRec, &req, dims, mapCenter,
                         azi_corr, dop_corr,
                         bx_map, by_map, bz_map, lx_map, ly_map, lz_map, dop_map)) {
            free(azi_corr); free(dop_corr);
            free(bx_map); free(by_map); free(bz_map);
            free(lx_map); free(ly_map); free(lz_map); free(dop_map);
            printf("error, frame skipped.\n");
            continue;
        }
        
        free(azi_corr); free(dop_corr);
        printf("done.\n");
        
        /* ====== */
        /* Output */
        /* ====== */
        
        printf("output: ");
        
        if (writeOutput(inRec, outRec, &req, mapCenter,           // Clean up done inside subroutine
                        bx_map, by_map, bz_map, lx_map, ly_map, lz_map, dop_map)) {
            printf("error, frame skipped.\n");
            continue;
        }
        
        printf("done.\n");
        
#endif      // TEST
        
    }   // irec
    
    /* Clean up */
    
    free(azi_work);
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    return 0;
    
}       // DoIt


// =============================================

/*
 * Read named array from record
 *
 */

int readDRMSArrayF(DRMS_Record_t *inRec, float *data, const char *segName)
{
 
    int status = 0;
    
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, segName);
    if (!inSeg) return 1;
    
    DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status || !inArray) return 1;
    
    int nx = inArray->axis[0], ny = inSeg->axis[1];
    int nxny = nx * ny;     // already checked
    
    float *inData = (float *)inArray->data;
    memcpy(data, inData, nxny * sizeof(float));        // copy over
    drms_free_array(inArray);
    
    return 0;
    
}

int readDRMSArrayD(DRMS_Record_t *inRec, double *data, const char *segName)
{
    
    int status = 0;
    
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRec, segName);
    if (!inSeg) return 1;
    
    DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
    if (status || !inArray) return 1;
    
    int nx = inArray->axis[0], ny = inSeg->axis[1];
    int nxny = nx * ny;     // already checked
    
    double *inData = (double *)inArray->data;
    memcpy(data, inData, nxny * sizeof(double));        // copy over
    drms_free_array(inArray);
    
    return 0;
    
}

// =============================================

/*
 * Correct azimuth
 * Input record set, currect record number, chunk size, chunks of azimuth to be corrected, dims
 * Output corrected frame
 *
 */

int correct_azi(DRMS_RecordSet_t *inRS, int irec, int nt, int *dims,
                float *azi_work, float *azi_corr)
{
    
    int status = 0;
    
    int nt2 = nt / 2;
    int nx = dims[0], ny = dims[1], nxny = nx * ny;
    int flip_thr = FLIP_THR;
    
#ifdef NOAZI
    
    printf("no correction by request, ");
    
    DRMS_Record_t *inRec = inRS->records[irec];
    if (readDRMSArrayF(inRec, azi_corr, "azimuth")) return 1;
    
#else       // default
    
    int edgeRec = ((irec < nt2) || (irec >= (inRS->n - nt2)));    // first/last NT/2 recs
    
    if (edgeRec) {                      // direct copy
        
        DRMS_Record_t *inRec = inRS->records[irec];
        if (readDRMSArrayF(inRec, azi_corr, "azimuth")) return 1;
        
    } else {                            // insert record #irec+nt2 to last frame of azi_work, then run correction
        
        DRMS_Record_t *inRec = inRS->records[irec + nt2];
        if (readDRMSArrayF(inRec, azi_corr, "azimuth")) return 1;
        
        for (int i = 0; i < nt - 1; i++) {      // Move every frame in azi_work forward
            memcpy(azi_work + i * nxny, azi_work + (i + 1) * nxny, nxny * sizeof(float));
        }
        memcpy(azi_work + (nt - 1) * nxny, azi_corr, nxny * sizeof(float));  // Insert new frame at the end
        
        // Correction using Brian's module
        
        azim_mindiff_jsoc_(azi_work, azi_corr, &nx, &ny, &nt2, &flip_thr);
        
    }
    
#endif
    
    return 0;
    
}

// =============================================

/*
 * Correct Doppler correction
 *
 */

int correct_dop(DRMS_Record_t *inRec, int *dims, float *dop_corr)
{
    
    int status = 0;
    
    int nx = dims[0], ny = dims[1], nxny = nx * ny;
    
    // Read
    
    double *dop_in = (double *) (malloc(nxny * sizeof(double)));
    if (readDRMSArrayD(inRec, dop_in, vSegName)) return 1;

#ifdef NODOP
    
    printf("no correction by request, ");
    
    for (int i = 0; i < nxny; i++) {
        dop_corr[i] = dop_in[i] / 100.;                // cm/s => m/s
    }

#else
    
    // Correction
    
    double *dv = (double *) (malloc(nxny * sizeof(double)));      // correct at eacu pixel in cm/s
    double doppbias = drms_getkey_double(inRec, "DOPPBIAS", &status);    // Keywords
    
    if (getDoppCorr(inRec, dims, dv)) {
        free(dv);
        return 1;
    }
    
    // copy over & subtract vr as placeholder during testing
    
    for (int i = 0; i < nxny; i++) {
        dop_corr[i] = (dop_in[i] - dv[i] - doppbias) / 100.;                // cm/s => m/s
    }
    
    // Clean up
    
    free(dv);
    
#endif
    
    free(dop_in);
    
    return 0;
    
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

int getDoppCorr(DRMS_Record_t *inRec, int *dims, double *doppCorr)
{
    
    int status = 0;
    
    int nx = dims[0], ny = dims[1], nxny = nx * ny;
    
    struct ephemeris ephem;
    if (getEphemeris(inRec, &ephem)) {
        return 1;
    }
    
 //   printf("xc=%f, yc=%f, lonc=%f\nnx=%d, ny=%d\n", ephem.disk_xc, ephem.disk_yc, ephem.disk_lonc, nx, ny);
 //   printf("asd=%f, pa=%f, rSun=%f\n", ephem.asd, ephem.pa, ephem.rSun);
    
    // Lon/lat
    
    double *lon_nobp = (double *) (malloc(nxny * sizeof(double)));      // assuming b=0, p=0, in rad
    double *lat_nobp = (double *) (malloc(nxny * sizeof(double)));
    double *lon = (double *) (malloc(nxny * sizeof(double)));
    double *lat = (double *) (malloc(nxny * sizeof(double)));
    double x, y, lon_t, lat_t;
    double rho, sinlat, coslat, sig, mu, chi;
    int ind;
    
    for (int row = 0; row < ny; row++) {
        for (int col = 0; col < nx; col++) {
            
            ind = row * nx + col;
            x = (col - ephem.disk_xc) / ephem.rSun;       // normalized pixel address
            y = (row - ephem.disk_yc) / ephem.rSun;       // normalized pixel address
            
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
    double k = sin(ephem.asd);      // arcsin(1/215.)
    
    double crota2 = ephem.pa * (-1.);        // In rad already
    double obsv_x = vw * cos(crota2) + vn * sin(crota2);
    double obsv_y = - vw * sin(crota2) + vn * cos(crota2);
    
    double rsun_ref = 6.96e10;                // cm
    double coslatc = cos(ephem.disk_latc);
    double sinlon;
    double sinlat_nobp, sinlon_nobp, coslat_nobp;
    
//    printf("diffrot=%lf, %lf, %lf\n", diffrot[0], diffrot[1], diffrot[2]);
    
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
                              (diffrot[0] + diffrot[1] * pow(sinlat, 2) + diffrot[2] * pow(sinlat, 4)));        // cm/s
            
        }
    }
    
    //
    
    free(lon_nobp); free(lat_nobp);     // clean up
    free(lon); free(lat);
    
    return 0;
    
}

// =============================================

/*
 * Perform mapping on field, los vector, and Doppler
 *
 */

int cgem_mapping(DRMS_Record_t *inRec, struct reqInfo *req, int *dims_in, float *mapCenter,
                 float *azi_corr, float *dop_corr,
                 float *bx_map, float *by_map, float *bz_map,
                 float *lx_map, float *ly_map, float *lz_map, float *dop_map)
{
    
    int status = 0;
    
    struct ephemeris ephem;
    if (getEphemeris(inRec, &ephem)) {
        return 1;
    }
    
    // Read data
    
    float *fld = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    if (readDRMSArrayF(inRec, fld, "field")) {
        free(fld);
        return 1;
    }
    
    float *inc = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    if (readDRMSArrayF(inRec, inc, "inclination")) {
        free(fld); free(inc);
        return 1;
    }
    
    float *azi = azi_corr, *dop = dop_corr;     // corrected
    
    // Vector transform
    
    float *bx = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    float *by = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    float *bz = (float *) (malloc(dims_in[0] * dims_in[1] * sizeof(float)));
    
    vectorTransform(fld, inc, azi, bx, by, bz, dims_in, &ephem);        // no error checking
    
    // Find position
    
    int dims_out[2] = {req->ncol, req->nrow};
    float *u_out = (float *) (malloc(dims_out[0] * dims_out[1] * sizeof(float)));
    float *v_out = (float *) (malloc(dims_out[0] * dims_out[1] * sizeof(float)));
    
    findCoord(req, &ephem, u_out, v_out, mapCenter);                    // no error checking; Stonyhurst
    
    // Sampling
    
    performSampling(bx, bx_map, u_out, v_out, dims_in, dims_out);
    performSampling(by, by_map, u_out, v_out, dims_in, dims_out);
    performSampling(bz, bz_map, u_out, v_out, dims_in, dims_out);
    performSampling(dop, dop_map, u_out, v_out, dims_in, dims_out);
    
    // Obtain LOS vectors
    
    getLOSVector(req, &ephem, mapCenter, lx_map, ly_map, lz_map);       // no error checking
    
    // Clean up
    
    free(fld); free(inc);
    free(bx); free(by); free(bz);
    free(u_out); free(v_out);
    
    return 0;
    
}

// =============================================

/*
 * Output all arrays, clean up
 *
 */

int writeOutput(DRMS_Record_t *inRec, DRMS_Record_t *outRec, struct reqInfo *req, float *mapCenter,
                float *bx_map, float *by_map, float *bz_map,
                float *lx_map, float *ly_map, float *lz_map, float *dop_map)
{
    
    int status = 0;
    
    int dims_out[2] = {req->ncol, req->nrow};
    int segNum = 0;
    
    DRMS_Segment_t *outSeg;
    DRMS_Array_t *outArray;
    
    float *outArrays[] = {bx_map, by_map, bz_map, lx_map, ly_map, lz_map, dop_map};
    int nSegs = ARRLENGTH(outSegNames), nArrays = ARRLENGTH(outArrays);
    if (nSegs != nArrays) {
        for (int i = 0; i < nArrays; i++) free(outArrays[i]);
        return 1;
    }
        
    // Write
    
    for (int iSeg = 0; iSeg < nSegs; iSeg++) {
        outSeg = drms_segment_lookup(outRec, outSegNames[iSeg]);
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims_out, outArrays[iSeg], &status);
        if (status) {
            for (int i = iSeg; i < nSegs; i++) free(outArrays[i]);
            return 1;
        }
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->israw = 0;
        outArray->bzero = outSeg->bzero;
        outArray->bscale = outSeg->bscale;
        status = drms_segment_write(outSeg, outArray, 0);
        drms_free_array(outArray);
        if (status) {
            for (int i = iSeg + 1; i < nSegs; i++) free(outArrays[i]);
            return 1;
        }
    }
    
    // Keys
    
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
    
    drms_setkey_float(outRec, "CRPIX1", req->ncol/2. + 0.5);
    drms_setkey_float(outRec, "CRPIX2", req->nrow/2. + 0.5);
    drms_setkey_float(outRec, "CRVAL1", mapCenter[0]);      // Stonyhurst
    drms_setkey_float(outRec, "CRVAL2", mapCenter[1]);
    drms_setkey_float(outRec, "CDELT1", req->dx / RADSINDEG);
    drms_setkey_float(outRec, "CDELT2", req->dy / RADSINDEG);
    drms_setkey_string(outRec, "CUNIT1", "degree");
    drms_setkey_string(outRec, "CUNIT2", "degree");
    
    TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
    tnow = (double)time(NULL);
    tnow += UNIX_epoch;
    drms_setkey_time(outRec, "DATE", tnow);
    
    char key[64];
    snprintf (key, 64, "CRLN-%s", wcsCode[req->proj]);
    drms_setkey_string(outRec, "CTYPE1", key);
    snprintf (key, 64, "CRLT-%s", wcsCode[req->proj]);
    drms_setkey_string(outRec, "CTYPE2", key);
    drms_setkey_float(outRec, "CROTA2", 0.0);
    
    drms_setkey_string(outRec, "BUNIT_000", "Mx/cm^2");
    drms_setkey_string(outRec, "BUNIT_001", "Mx/cm^2");
    drms_setkey_string(outRec, "BUNIT_002", "Mx/cm^2");
    drms_setkey_string(outRec, "BUNIT_003", " ");
    drms_setkey_string(outRec, "BUNIT_004", " ");
    drms_setkey_string(outRec, "BUNIT_005", " ");
    drms_setkey_string(outRec, "BUNIT_006", "m/s");
 
    return 0;
    
}

// =============================================

/*
 * Output cutout arrays for testing
 *
 */

int writeOutputTest(DRMS_Record_t *inRec, DRMS_Record_t *outRec, int *dims,
                    float *azi_corr, float *dop_corr)
{
    
    int status = 0;
    
    DRMS_Segment_t *outSeg;
    DRMS_Array_t *outArray;
    
    // Write
    
    outSeg = drms_segment_lookup(outRec, "azi_corr");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, azi_corr, &status);
    if (status) {
        free(azi_corr); free(dop_corr);
        return 1;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    drms_free_array(outArray);
    if (status) {
        free(dop_corr);
        return 1;
    }
    
    outSeg = drms_segment_lookup(outRec, "dop_corr");
    outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, dop_corr, &status);
    if (status) {
        free(dop_corr);
        return 1;
    }
    outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
    outArray->israw = 0;
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    drms_free_array(outArray);
    if (status) {
        return 1;
    }
    
    // Keyw
    
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
    
    int noazi = 0, nodop = 0;
#ifdef NOAZI
    noazi = 1;
#endif
#ifdef NODOP
    nodop = 1;
#endif
    
    drms_setkey_int(outRec, "NOAZI", noazi);
    drms_setkey_int(outRec, "NODOP", nodop);
    
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
    
    float crota2 = drms_getkey_float(inRec, "CROTA2", &status);    // rotation
    double sina = sin(crota2 * RADSINDEG);
    double cosa = cos(crota2 * RADSINDEG);
    
    ephem->pa = - crota2 * RADSINDEG;
    ephem->disk_latc = drms_getkey_float(inRec, "CRLT_OBS", &status) * RADSINDEG;
    ephem->disk_lonc = drms_getkey_float(inRec, "CRLN_OBS", &status) * RADSINDEG;
    
    float crvalx = drms_getkey_float(inRec, "CRVAL1", &status);
    float crvaly = drms_getkey_float(inRec, "CRVAl2", &status);
    float crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);
    float crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
    float cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
    ephem->disk_xc = PIX_X(0.0,0.0) - 1.0;        // Center of disk in pixel, starting at 0
    ephem->disk_yc = PIX_Y(0.0,0.0) - 1.0;
    
    float dSun = drms_getkey_float(inRec, "DSUN_OBS", &status);
    float rSun_ref = drms_getkey_float(inRec, "RSUN_REF", &status);
    if (status) rSun_ref = 6.96e8;
    
    ephem->asd = asin(rSun_ref/dSun);                               // in rad
    ephem->rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;       // in pixel
    
    ephem->t_rec = drms_getkey_time(inRec, "T_REC", &status);
    
    return 0;
    
}



// =============================================

/*
 * Performing local vector transformation
 *  xyz: z refers to vertical (radial) component, x EW (phi), y NS
 *
 */

void vectorTransform(float *fld, float *inc, float *azi, float *bx, float *by, float *bz,
                     int *dims, struct ephemeris *ephem)
{
    
    int ncol = dims[0], nrow = dims[1];
    
    double b_xi, b_eta, b_zeta;
    double bx0, by0, bz0;
    
    double x, y;        // solar x, y
    double lat, lon;    // lat / lon for current point
    double rho, sinlat, coslat, sig, mu, chi;
    
    int ind;
    for (int row = 0; row < nrow; row++) {
        for (int col = 0; col < ncol; col++) {
            
            ind = row * ncol + col;
            
            x = (col - ephem->disk_xc) / ephem->rSun;      // normalized by disk radius
            y = (row - ephem->disk_yc) / ephem->rSun;
            
            // find lat/lon (disk center lon = 0)
            if (img2sphere(x, y, ephem->asd, ephem->disk_latc, 0.0, ephem->pa,
                           &rho, &lat, &lon, &sinlat, &coslat, &sig, &mu, &chi)) {
                bx[ind] = by[ind] = bz[ind] = DRMS_MISSING_FLOAT;
                continue;
            }
            
            // field vector
            
            b_xi = - fld[ind] * sin(inc[ind] * RADSINDEG) * sin(azi[ind] * RADSINDEG);
            b_eta = fld[ind] * sin(inc[ind] * RADSINDEG) * cos(azi[ind] * RADSINDEG);
            b_zeta = fld[ind] * cos(inc[ind] * RADSINDEG);
            
            // convert (disk center lon = 0)
            
            img2helioVector (b_xi, b_eta, b_zeta, &bx0, &by0, &bz0,
                             lon, lat, 0.0, ephem->disk_latc, ephem->pa);
            
            bx[ind] = bx0;
            by[ind] = by0;
            bz[ind] = bz0;
            
        }
    }
    
}

// =============================================

/*
 * Compute the coordinates to be sampled on cutout image
 *
 */

void findCoord(struct reqInfo *req, struct ephemeris *ephem,
               float *u_out, float *v_out, float *mapCenter)
{
    
    // Rotate
    // Differential rotation rate in urad/s
    // proj/lev1.5_hmi/libs/lev15/rotcoef_file.txt
    
    double latc = req->latref;
    double difr = 2.7139 - 0.405 * pow(sin(latc),2) - 0.422 * pow(sin(latc),4);
    double dt = ephem->t_rec - req->tref;
    double lonc = req->lonref + dt * difr * 1.0e-6;      // urad/s to rad
    mapCenter[0] = lonc / RADSINDEG;
    mapCenter[1] = latc / RADSINDEG;
    
    // print some info
    char tstr[100];
    sprint_time(tstr, ephem->t_rec, "TAI", 0);
    //    printf("t_rec=%s, lonc=%f\n", tstr, lonc/RADSINDEG);
    
    double x, y;        // coordinate in final map
    double lat, lon;    // lon, lat of each pix
    double u, v;    // pixel address
    
    int ind_map;
    for (int row = 0; row < req->nrow; row++) {
        for (int col = 0; col < req->ncol; col++) {
            
            ind_map = row * req->ncol + col;
            x = (col + 0.5 - req->ncol/2.) * req->dx;
            y = (row + 0.5 - req->nrow/2.) * req->dy;
            
            //
            
            if (plane2sphere (x, y, latc, lonc, &lat, &lon, req->proj)) {
                u_out[ind_map] = -1;
                v_out[ind_map] = -1;
                continue;
            }
            
            if (sphere2img (lat, lon, ephem->disk_latc, 0.0, &u, &v,
                            ephem->disk_xc, ephem->disk_yc, ephem->rSun, ephem->pa,
                            0., 0., 0., 0.)) {
                u_out[ind_map] = -1;
                v_out[ind_map] = -1;
                continue;
            }
            
            u_out[ind_map] = u;
            v_out[ind_map] = v;
            
        }
    }
    
}


// =============================================

/*
 * Convert LOS vectors to local xyz coordinate
 *
 */

void getLOSVector(struct reqInfo *req, struct ephemeris *ephem, float *mapCenter,
                  float *lx, float *ly, float *lz)
{
    
    double lonc = mapCenter[0] * RADSINDEG;
    double latc = mapCenter[1] * RADSINDEG;
    
    double x, y;        // coordinate in final map
    double lat, lon;    // lon, lat of each pix
    double lx0, ly0, lz0;
    
    int ind_map;
    for (int row = 0; row < req->nrow; row++) {
        for (int col = 0; col < req->ncol; col++) {
            
            ind_map = row * req->ncol + col;
            x = (col + 0.5 - req->ncol/2.) * req->dx;
            y = (row + 0.5 - req->nrow/2.) * req->dy;
            
            if (plane2sphere (x, y, latc, lonc, &lat, &lon, req->proj)) {
                lx[ind_map] = DRMS_MISSING_FLOAT;
                ly[ind_map] = DRMS_MISSING_FLOAT;
                lz[ind_map] = DRMS_MISSING_FLOAT;
                continue;
            }
            
            // convert (disk center lon = 0)
            
            img2helioVector (0.0, 0.0, 1.0, &lx0, &ly0, &lz0,
                             lon, lat, 0.0, ephem->disk_latc, ephem->pa);
            
            lx[ind_map] = lx0;
            ly[ind_map] = ly0;
            lz[ind_map] = lz0;
            
        }
    }
    
}

// =============================================

/*
 * Perform interpolation using lib function
 *
 */


void performSampling(float *inData, float *outData, float *u, float *v,
                     int *dims_in, int *dims_out)
{
    // Use lib interpolation method
    struct fint_struct pars;
    
    switch (INTERPOPT) {
        case 0:            // Wiener, 6 order, 1 constraint
            init_finterpolate_wiener(&pars, 6, 1, 6, 2, 1, 1, NULL, dpath);
            break;
        case 1:            // Cubic convolution
            init_finterpolate_cubic_conv(&pars, 1., 3.);
            break;
        case 2:            // Bilinear
            init_finterpolate_linear(&pars, 1.);
            break;
        case 3:            // Near neighbor
        default:
            break;
    }
    
    int ind;
    int ncol = dims_out[0], nrow = dims_out[1];
    if (INTERPOPT == 3) {            // Aug 6 2013, Xudong
        for (int row = 0; row < nrow; row++) {
            for (int col = 0; col < ncol; col++) {
                ind = row * ncol + col;
                outData[ind] = nnb(inData, dims_in[0], dims_in[1], u[ind], v[ind]);
            }
        }
    } else {
        // lib function
        finterpolate(&pars, inData, u, v, outData,
                     dims_in[0], dims_in[1], dims_in[0], ncol, nrow, ncol, DRMS_MISSING_FLOAT);
    }
    
}


// =============================================

/*
 * Near neighbor interpolation
 *
 */

float nnb (float *f, int nx, int ny, float x, float y)
{
    
    if (x <= -0.5 || y <= -0.5 || x > nx - 0.5 || y > ny - 0.5)
        return DRMS_MISSING_FLOAT;
    int ilow = floor (x);
    int jlow = floor (y);
    int i = ((x - ilow) > 0.5) ? ilow + 1 : ilow;
    int j = ((y - jlow) > 0.5) ? jlow + 1 : jlow;
    return f[j * nx + i];
    
}


