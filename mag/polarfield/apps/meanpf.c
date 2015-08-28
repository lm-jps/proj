/*
 * Module name:		meanpf.c
 *
 * Description:
 *  Taking ab HMI magnetogram and compute mean field
 *  in various lat/lon windows, adapted from polmf.c
 *  The original, double 180x4 segment is splitted into
 *  four 180-element, Rice compressed segments
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *  v1.0		Aug 27 2015
 *
 * Issues:
 *  v1.0    No error checking as of now
 *
 *
 * Example:
 *  meanpf "in=hmi.M_720s[2011.11.01]"
 *
 */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cartography.c"

/* ############# Macros ############# */

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define SHOW(msg) {printf("%s", msg); fflush(stdout);}
#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s (status=%d)\n", msg, status); return(status);}

#define FREE_ARR(arr) {if (arr) {free(arr);}}
#define DRMS_FREE_ARR(arr) {if (arr) {drms_free_array(arr);}}

#define EPSL (1.0E-5)
#define RAD2DEG       (180. / M_PI)
#define RAD2ARCSEC      (648000. / M_PI)

// apodize radius
#define APOD    (0.995)

// Br threshold for FWTLAT
#define BRTHRESH    (120.)

// Test
//#define TESTOUT     1

/* ############# Structures ############# */

struct ephemeris {
    double rsun_obs;    // solar disk, in arcsec
    double asd;         // solar disk, in radian
    double xc, yc;      // disk center in pixel, lower-left (0,0)
    double cdelt;       // plate scale, in arcsec, assuming x & y same
    double p0, b0;      // p angle, b angle, in radian
    double rsun;        // in pixel
    double da;          // pixel area, in 1.e16 cm2 (1 Mm2)
};

struct fwt {
    double fluxn[2];       // flux below 40 lat for north [pos,neg], in 1e20 Mx
    double fluxs[2];       // flux below 40 lat for south [pos,neg], in 1e20 Mx
    double latn[2];        // flux weighted centroid for north
    double lats[2];        // flux weighted centroid for south
};

/* ############# Pre-emb ############# */

/* Timing by T. P. Larson */

double getwalltime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000.0 + tv.tv_usec/1000.0;
}

double getcputime(double *utime, double *stime)
{

  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  *utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec/1000.0;
  *stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec/1000.0;
  return *utime + *stime;
}


// Get ephemeris
int get_ephemeris(DRMS_Record_t *inRec, struct ephemeris *ephem);

// Get lon/lat for each pixel
void get_lonlat(struct ephemeris *ephem, int *dims, double *lon, double *lat, double *wght);

// Convert LOS field to radial
void los2radial(struct ephemeris *ephem, int *dims, double *bl, double *br, double *wght);

// Compute mean values in certain bins
void compute_mf(int *dims, double *br, double *bl, double *lon, double *lat,
                double *wght, int *outDims, double *mf);

// Compute flux weighted lat
void compute_fwtlat(int *dims, double *br, double *lon, double *lat, double *wght, struct fwt *fwtlat);

// Set polar field keywords
void set_mf_keys(DRMS_Record_t *outRec, int *outDims, double *mf, struct fwt *fwtlat);


/* ################################################# */
/* ################## Main Module ################## */
/* ################################################# */

char *module_name = "meanpf";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", "", "Input query"},
    {ARG_STRING, "out", "hmi.meanpf_720s", "Output query"},
    {ARG_END}
};

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
    
    /* ============== */
    /* Get parameters */
    /* ============== */
    
    char *inQuery = (char *)params_get_str(&cmdparams, "in");
    char *outQuery = (char *)params_get_str(&cmdparams, "out");
    
    /* ============= */
    /* Input records */
    /* ============= */
    
    DRMS_RecordSet_t *inRS = NULL;
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("Input records error.");
    
    int nrecs = inRS->n;
    
    /* ============== */
    /* Output records */
    /* ============== */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output records not created.");
    
    /* ============ */
    /* Loop records */
    /* ============ */
    
    printf("Start, %d records in total\n", nrecs);
    
    for (int irec = 0; irec < nrecs; irec++) {

        // Time measure
        
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("Record %d of %d done, %.2f ms wall time, %.2f ms cpu time\n",
               irec + 1, nrecs, wt - wt1, ct - ct1);
        
        // Timing
        
        wt1 = getwalltime();
        ct1 = getcputime(&ut1, &st1);
        
        // Input ecord
        
        DRMS_Record_t *inRec = inRS->records[irec];
        TIME t_rec = drms_getkey_time(inRec, "T_REC", &status);
        
        char t_rec_str[100];
        sprint_time(t_rec_str, t_rec, "TAI", 0);
        printf("==============\nRecord #%d, [%s]\n", irec, t_rec_str);
        
        struct ephemeris ephem;
        status = get_ephemeris(inRec, &ephem);
        if (status) {
            SHOW("Input record ephemeris error, record skipped.\n");
        }
        
        DRMS_Segment_t *inSeg = drms_segment_lookupnum(inRec, 0);
        DRMS_Array_t *inArray = NULL;
        inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
        if (status) {
            SHOW("Input array read error, record skipped.\n");
            if (inArray) drms_free_array(inArray);
            continue;
        }
        
        int nx = inArray->axis[0], ny = inArray->axis[1];       // size
        int nxny = nx * ny;
        int dims[2] = {nx, ny};
        
        // Lon, Lat
        
        SHOW("Get lon/lat...\n");
        
        double *lon = NULL, *lat = NULL, *wght = NULL;
        lon = (double *) (calloc(sizeof(double), nxny));
        lat = (double *) (calloc(sizeof(double), nxny));
        wght = (double *) (calloc(sizeof(double), nxny));
        
        get_lonlat(&ephem, dims, lon, lat, wght);
        
        // LOS to Radial
        
        double *bl = NULL, *br = NULL;
        bl = (double *) (inArray->data);
        br = (double *) (calloc(sizeof(double), nxny));
        
        SHOW("Convert bl to br...\n");
        los2radial(&ephem, dims, bl, br, wght);
        
        // Compute MF
        
        SHOW("Compute mf...\n");
        
        int nlat = 180;
        int outDims_mf[2] = {nlat,4};
        double *mf = (double *) (calloc(sizeof(double), outDims_mf[0]*outDims_mf[1]));
        /*
         double *mf_br = mf;
         double *mf_bl = mf + row_offset;
         double *w = mf + row_offset * 2;
         double *num = mf + row_offset * 3;
         */
        
        compute_mf(dims, br, bl, lon, lat, wght, outDims_mf, mf);
        
        // Compute flux weighted centroid
        
        struct fwt fwtlat;
        compute_fwtlat(dims, br, lon, lat, wght, &fwtlat);
        
        fwtlat.fluxn[0] *= ephem.da;            // 1e20 Mx
        fwtlat.fluxn[1] *= ephem.da;
        fwtlat.fluxs[0] *= ephem.da;
        fwtlat.fluxs[1] *= ephem.da;
        
        // Split into four 180-element arrays
        
        double *mf_br = (double *) (malloc(nlat * sizeof(double)));
        double *mf_bl = (double *) (malloc(nlat * sizeof(double)));
        double *w = (double *) (malloc(nlat * sizeof(double)));
        double *num = (double *) (malloc(nlat * sizeof(double)));
        for (int i = 0; i < nlat; i++) {
            mf_br[i] = mf[i];
            mf_bl[i] = mf[nlat + i];
            w[i] = mf[nlat * 2 + i];
            num[i] = mf[nlat * 3 + i];
        }
        
        // Clean up #1
        
        drms_free_array(inArray);       // bl
        FREE_ARR(br); FREE_ARR(wght); FREE_ARR(lon); FREE_ARR(lat);

        // Output record
        
        SHOW("Output...\n");
        
        int outDims[1] = {nlat};
        DRMS_Record_t *outRec = outRS->records[irec];
        DRMS_Segment_t *outSeg_br = drms_segment_lookup(outRec, "mf_br");
        DRMS_Array_t *outArray_br = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, mf_br, &status);
        status = drms_segment_write(outSeg_br, outArray_br, 0);
        
        DRMS_Segment_t *outSeg_bl = drms_segment_lookup(outRec, "mf_bl");
        DRMS_Array_t *outArray_bl = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, mf_bl, &status);
        status = drms_segment_write(outSeg_bl, outArray_bl, 0);
        
        DRMS_Segment_t *outSeg_w = drms_segment_lookup(outRec, "w");
        DRMS_Array_t *outArray_w = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, w, &status);
        status = drms_segment_write(outSeg_w, outArray_w, 0);
        
        DRMS_Segment_t *outSeg_num = drms_segment_lookup(outRec, "num");
        DRMS_Array_t *outArray_num = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, num, &status);
        status = drms_segment_write(outSeg_num, outArray_num, 0);
        
        // Keywords
        
        SHOW("Set keywords...\n");
        
        drms_copykey(outRec, inRec, "T_REC");
        
        char timebuf[1024];
        double UNIX_epoch = -220924792.0;	/* 1970.01.01_00:00:00_UTC */
        sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
        drms_setkey_string(outRec, "DATE", timebuf);
        
        drms_setkey_double(outRec, "THRESH", BRTHRESH);
        drms_setkey_double(outRec, "DA", ephem.da);
        
        drms_copykey(outRec, inRec, "QUALITY");
        drms_copykey(outRec, inRec, "CRLT_OBS");
        drms_copykey(outRec, inRec, "CRLN_OBS");
        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_copykey(outRec, inRec, "OBS_VR");
        drms_copykey(outRec, inRec, "INTERVAL");        // for MDI
        
        // Keywords of mf
        
        set_mf_keys(outRec, outDims_mf, mf, &fwtlat);       // do nothing for now
        
        // Clean up
        
        DRMS_FREE_ARR(outArray_br); DRMS_FREE_ARR(outArray_bl);
        DRMS_FREE_ARR(outArray_w); DRMS_FREE_ARR(outArray_num);
        FREE_ARR(mf);
        
        // Time measure
        
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("Record %d of %d done, %.2f ms wall time, %.2f ms cpu time\n",
               irec + 1, nrecs, wt - wt1, ct - ct1);

    }   // irec
    
    /* ============= */
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    wt = getwalltime();
    ct = getcputime(&ut, &st);
    printf("==============\nTotal time spent: %.2f ms wall time, %.2f ms cpu time\n",
           wt - wt0, ct - ct0);
    
    return DRMS_SUCCESS;
    
}   // DoIt

/* ############# Helper functions ############# */

// Get ephemeris

int get_ephemeris(DRMS_Record_t *inRec, struct ephemeris *ephem)
{
    
    int status = 0;
    
    ephem->rsun_obs = drms_getkey_double(inRec, "RSUN_OBS", &status);   // in arcsec
    if (status) return status;
    ephem->asd = ephem->rsun_obs / RAD2ARCSEC;      // in rad
    
    ephem->xc = drms_getkey_double(inRec, "CRPIX1", &status) - 1.;      // start from (0,0)
    if (status) return status;
    ephem->yc = drms_getkey_double(inRec, "CRPIX2", &status) - 1.;
    if (status) return status;
    
    double cdelt1 = drms_getkey_double(inRec, "CDELT1", &status);
    if (status) return status;
    double cdelt2 = drms_getkey_double(inRec, "CDELT2", &status);
    if (status) return status;
    if (fabs(cdelt1 - cdelt2) > EPSL) return 1;
    ephem->cdelt = cdelt1;
    
    ephem->rsun = ephem->rsun_obs / ephem->cdelt;     // in pixel
    
    ephem->p0 = drms_getkey_double(inRec, "CROTA2", &status) / RAD2DEG * (-1.);
    if (status) return status;
    
    ephem->b0 = drms_getkey_double(inRec, "CRLT_OBS", &status) / RAD2DEG;
    if (status) return status;
    
    double dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status) / 1.e6;   // in Mm, or 1.e8 cm
    double dx = cdelt1 * ((dsun_obs - 696.) / RAD2ARCSEC);       // pixel size in Mm, or 1.e8 cm
    ephem->da = dx * dx;        // in 1.e16 cm^2
    
    return 0;
    
}

// Get lon/lat for each pixel

void get_lonlat(struct ephemeris *ephem, int *dims, double *lon, double *lat, double *wght)
{
    
    int ncol = dims[0], nrow = dims[1];
    int idx;
    int status = 0;
    
    double x, y, r;        // dist to center, in solar radius
    double rho, sinlat, coslat, sig, mu, chi;
    double lat0, lon0;
    
    for (int row = 0; row < nrow; row++) {
        y = (row - ephem->yc) / ephem->rsun;
        for (int col = 0; col < ncol; col++) {
            x = (col - ephem->xc) / ephem->rsun;
            r = hypot(x, y);
            idx = row * ncol + col;
            if (r <= APOD) {
                status = img2sphere(x, y, ephem->asd, ephem->b0, 0.0, ephem->p0,
                                    &rho, &lat0, &lon0, &sinlat, &coslat, &sig, &mu, &chi);
                if (status) {
                    lon[idx] = lat[idx] =  wght[idx] = DRMS_MISSING_DOUBLE;
                } else {
                    lon[idx] = ((lon0 > M_PI) ? (lon0 - M_PI * 2) : lon0) * RAD2DEG;
                    lat[idx] = lat0 * RAD2DEG;
                    wght[idx] = 1.0 / mu;       // changed to mu from cos(rho) on Apr 16
                }
            } else {
                lon[idx] = lat[idx] = wght[idx] = DRMS_MISSING_DOUBLE;
            }
        }
    }
    
}

// Convert LOS field to radial

void los2radial(struct ephemeris *ephem, int *dims, double *bl, double *br, double *wght)
{
    
    int nxny = dims[0] * dims[1];
    int idx;
    
    double x, y, r;        // dist to center, in solar radius
    double crho;
    
    for (int idx = 0; idx < nxny; idx++) {
        if (isnan(wght[idx])) {
            br[idx] = bl[idx] = DRMS_MISSING_DOUBLE;
        } else {
            br[idx] = bl[idx] * wght[idx];      // wght = 1/mu
        }
    }
    
}

// Compute mean values in certain bins

void compute_mf(int *dims, double *br, double *bl, double *lon, double *lat,
                double *wght, int *outDims, double *mf)
{
    
    int row_offset = outDims[0];        // 180
    
    double *mf_br = mf;
    double *mf_bl = mf + row_offset;
    double *w = mf + row_offset * 2;
    double *num = mf + row_offset * 3;
    
    int nx = dims[0], ny = dims[1];
    int nxny = nx * ny;
    
    int ngood = 0;
    int ilat;
    double s, v1, v2;
    
    // Int Loop
    for (int idx = 0; idx < nxny; idx++) {
        if (isnan(lon[idx])) continue;
        ngood++;
        // lon
        if (fabs(lon[idx]) > 45.) continue;       // dob't use lon gt \pm 45
        // lat
        ilat = floor(lat[idx] + 90.);       // 0S, 180N
        if (ilat >= 180 || ilat < 0) continue;
        // compute
        mf_br[ilat] += (br[idx] * wght[idx]);
        mf_bl[ilat] += (bl[idx]);
        w[ilat] += wght[idx];
        num[ilat] ++;
    }
    
    // MF loop
    for (ilat = 0; ilat < 180; ilat ++) {
        num[ilat] = fabs(num[ilat]);
        if (num[ilat] < EPSL) {
            mf_br[ilat] = DRMS_MISSING_DOUBLE;
            mf_bl[ilat] = DRMS_MISSING_DOUBLE;
            continue;
        }
        mf_br[ilat] /= (w[ilat]);
        mf_bl[ilat] /= (num[ilat]);
    }
    
}

// Compute flux weighted lat

void compute_fwtlat(int *dims, double *br, double *lon, double *lat, double *wght, struct fwt *fwtlat)
{
    
    int nx = dims[0], ny = dims[1];
    int nxny = nx * ny;
    
    double fln[2] = {0, 0}, fls[2] = {0, 0};        // product of flux and lat
    double f, fl;
    
    for (int idx = 0; idx < nxny; idx++) {
        // sanction
        if (fabs(lon[idx]) > 45. || fabs(lat[idx]) > 40. || fabs(br[idx]) < BRTHRESH) continue;
        // flux, flux times lat
        f = fabs(br[idx] * wght[idx]);
        fl = f * lat[idx];
        // north pos
        if (lat[idx] >= 0 && br[idx] >= 0) {
            fwtlat->fluxn[0] += f; fln[0] += fl;
        }
        // north neg
        if (lat[idx] >= 0 && br[idx] < 0) {
            fwtlat->fluxn[1] += f; fln[1] += fl;
        }
        // south pos
        if (lat[idx] < 0 && br[idx] >= 0) {
            fwtlat->fluxs[0] += f; fls[0] += fl;
        }
        // sorth neg
        if (lat[idx] < 0 && br[idx] < 0) {
            fwtlat->fluxs[1] += f; fls[1] += fl;
        }
    }
    
    fwtlat->latn[0] = fln[0] / fwtlat->fluxn[0];
    fwtlat->latn[1] = fln[1] / fwtlat->fluxn[1];
    fwtlat->lats[0] = fls[0] / fwtlat->fluxs[0];
    fwtlat->lats[1] = fls[1] / fwtlat->fluxs[1];
    
    fwtlat->fluxn[0] *= (1.e-4);        // in 1.e4 G^2
    fwtlat->fluxn[1] *= (-1.e-4);
    fwtlat->fluxs[0] *= (1.e-4);
    fwtlat->fluxs[1] *= (-1.e-4);
    
}

// Set polar field keywords

void set_mf_keys(DRMS_Record_t *outRec, int *outDims, double *mf, struct fwt *fwtlat)
{
    
    int row_offset = outDims[0];        // 180
    
    double *mf_br = mf;
    double *mf_bl = mf + row_offset;
    double *w = mf + row_offset * 2;
    double *num = mf + row_offset * 3;
    
    int idxn, idxs;
    
    // 70-90 band
    
    double bn3 = 0, bs3 = 0, wn3 = 0, ws3 = 0, numn3 = 0, nums3 = 0;
    for (int i = 0; i < 20; i++) {        // 20 lat bins
        idxn = 179 - i;
        if (num[idxn] > EPSL) {     // skip NaN's and average available ones
            bn3 += mf_br[idxn] * w[idxn];
            wn3 += w[idxn];
            numn3 += num[idxn];
        }
        //
        idxs = i;
        if (num[idxs] > EPSL) {
            bs3 += mf_br[idxs] * w[idxs];
            ws3 += w[idxs];
            nums3 += num[idxs];
        }
    }
    
    // 50-60 and 60-70
    
    double bn2 = 0, bs2 = 0, wn2 = 0, ws2 = 0, numn2 = 0, nums2 = 0;
    double bn1 = 0, bs1 = 0, wn1 = 0, ws1 = 0, numn1 = 0, nums1 = 0;
    for (int i = 0; i < 10; i++) {
        idxn = 159 - i;
        bn2 += mf_br[idxn] * w[idxn];
        wn2 += w[idxn];
        numn2 += num[idxn];
        //
        idxn = 149 - i;
        bn1 += mf_br[idxn] * w[idxn];
        wn1 += w[idxn];
        numn1 += num[idxn];
        //
        idxs = i + 20;
        bs2 += mf_br[idxs] * w[idxs];
        ws2 += w[idxs];
        nums2 += num[idxs];
        //
        idxs = i + 30;
        bs1 += mf_br[idxs] * w[idxs];
        ws1 += w[idxs];
        nums1 += num[idxs];
    }
    
    // Set keys
    
    drms_setkey_double(outRec, "BANDN3", bn3 / wn3);
    drms_setkey_double(outRec, "BANDS3", bs3 / ws3);
    
    drms_setkey_double(outRec, "BANDN2", bn2 / wn2);
    drms_setkey_double(outRec, "BANDS2", bs2 / ws2);
    
    drms_setkey_double(outRec, "BANDN1", bn1 / wn1);
    drms_setkey_double(outRec, "BANDS1", bs1 / ws1);
    
    drms_setkey_double(outRec, "CAPN2", (bn2 + bn3) / (wn2 + wn3));
    drms_setkey_double(outRec, "CAPS2", (bs2 + bs3) / (ws2 + ws3));
    
    drms_setkey_double(outRec, "CAPN1", (bn1 + bn2 + bn3) / (wn1 + wn2 + wn3));
    drms_setkey_double(outRec, "CAPS1", (bs1 + bs2 + bs3) / (ws1 + ws2 + ws3));
    
    drms_setkey_double(outRec, "NUMN3", numn3);
    drms_setkey_double(outRec, "NUMS3", nums3);
    
    drms_setkey_double(outRec, "NUMN2", numn2);
    drms_setkey_double(outRec, "NUMS2", nums2);
    
    drms_setkey_double(outRec, "NUMN1", numn1);
    drms_setkey_double(outRec, "NUMS1", nums1);
    
    drms_setkey_double(outRec, "WN3", wn3);
    drms_setkey_double(outRec, "WS3", ws3);
    
    drms_setkey_double(outRec, "WN2", wn2);
    drms_setkey_double(outRec, "WS2", ws2);
    
    drms_setkey_double(outRec, "WN1", wn1);
    drms_setkey_double(outRec, "WS1", ws1);
    
    // Set keys
    
    drms_setkey_double(outRec, "FLUXN_P", fwtlat->fluxn[0]);
    drms_setkey_double(outRec, "FLUXN_N", fwtlat->fluxn[1]);
    drms_setkey_double(outRec, "FLUXS_P", fwtlat->fluxs[0]);
    drms_setkey_double(outRec, "FLUXS_N", fwtlat->fluxs[1]);
    
    drms_setkey_double(outRec, "FLATN_P", fwtlat->latn[0]);
    drms_setkey_double(outRec, "FLATN_N", fwtlat->latn[1]);
    drms_setkey_double(outRec, "FLATS_P", fwtlat->lats[0]);
    drms_setkey_double(outRec, "FLATS_N", fwtlat->lats[1]);
    
}
