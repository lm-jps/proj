/*
 * Module name:		usflux.c
 *
 * Description:
 *  Taking ab HMI magnetogram and compute unsigned flux
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *  v1.0		Apr 19 2016
 *
 * Issues:
 *  v1.0    No error checking as of now
 *
 *
 * Example:
 *  usflux "in=hmi.M_720s[2011.11.01]"
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
#define APOD    (0.990)
#define BRTHRESH    (0.0)
#define BLTHRESH    (0.0)

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

// Compute flux
void compute_flux(int *dims, double *br, double *bl, double *wght,
                  double *ntflux_bl, double *ntflux_br, double *usflux_bl, double *usflux_br,
                  double *w_bl, double *w_br, double *r_bl, double *r_br);



/* ################################################# */
/* ################## Main Module ################## */
/* ################################################# */

char *module_name = "usflux";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", "", "Input query"},
    {ARG_STRING, "out", "hmi_test.usflux_720s", "Output query"},
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
        
        wt1 = getwalltime();
        ct1 = getcputime(&ut, &st);
        
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
        
        // Compute flux
        
        SHOW("Compute usflux...\n");
        
        // LOS flux and radial (read) flux, for net and unsigned, not scaled by area
        double ntflux_bl = 0, ntflux_br = 0, usflux_bl = 0, usflux_br = 0;
        // sum of pixel, pixel area, percentage above threshold
        double w_bl = 0, w_br = 0, r_bl = 0, r_br = 0;
        
        compute_flux(dims, br, bl, wght,
                     &ntflux_bl, &ntflux_br, &usflux_bl, &usflux_br,
                     &w_bl, &w_br, &r_bl, &r_br);
        
        // Clean up #1
        
        drms_free_array(inArray);       // bl
        FREE_ARR(br); FREE_ARR(wght); FREE_ARR(lon); FREE_ARR(lat);
        
        // Output record
        
        SHOW("Output...\n");
        printf("mean bl: %.4f G\n", ntflux_bl / w_bl);
        
        DRMS_Record_t *outRec = outRS->records[irec];
        
        drms_copykey(outRec, inRec, "T_REC");
        
        char timebuf[1024];
        double UNIX_epoch = -220924792.0;	/* 1970.01.01_00:00:00_UTC */
        sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
        drms_setkey_string(outRec, "DATE", timebuf);
        
        double k = ephem.da * 1.e-4;            // 1d20 cm^2
        drms_setkey_double(outRec, "NTFLUX_BL", ntflux_bl * k);
        drms_setkey_double(outRec, "NTFLUX_BR", ntflux_br * k);
        drms_setkey_double(outRec, "USFLUX_BL", usflux_bl * k);
        drms_setkey_double(outRec, "USFLUX_BR", usflux_br * k);

        drms_setkey_double(outRec, "BRTHRESH", BRTHRESH);
        drms_setkey_double(outRec, "BLTHRESH", BRTHRESH);
        drms_setkey_double(outRec, "W_BL", w_bl);
        drms_setkey_double(outRec, "W_BR", w_br);
        drms_setkey_double(outRec, "R_BL", r_bl);
        drms_setkey_double(outRec, "R_BR", r_br);
        drms_setkey_double(outRec, "DA", ephem.da);
        drms_setkey_double(outRec, "APOD", APOD);
        
        drms_copykey(outRec, inRec, "QUALITY");
        drms_copykey(outRec, inRec, "CRLT_OBS");
        drms_copykey(outRec, inRec, "CRLN_OBS");
        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_copykey(outRec, inRec, "OBS_VR");
        
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

// Compute flux

void compute_flux(int *dims, double *br, double *bl, double *wght,
                  double *ntflux_bl, double *ntflux_br, double *usflux_bl, double *usflux_br,
                  double *w_bl, double *w_br, double *r_bl, double *r_br)
{
    int nx = dims[0], ny = dims[1];
    int nxny = nx * ny;
    
    // Int Loop
    for (int idx = 0; idx < nxny; idx++) {
        if (isnan(wght[idx])) continue;
        (*w_bl) += 1;
        (*w_br) += wght[idx];
        // compute
        if (fabs(bl[idx]) >= BLTHRESH) {
            (*ntflux_bl) += (bl[idx]);
            (*usflux_bl) += fabs(bl[idx]);
            (*r_bl) += 1;
        }
        if (fabs(br[idx]) >= BRTHRESH) {
            (*ntflux_br) += (br[idx] * wght[idx]);
            (*usflux_br) += fabs(br[idx] * wght[idx]);
            (*r_br) += wght[idx];
        }
    }
    
    (*r_bl) /= (*w_bl);
    (*r_br) /= (*w_br);
    
}

