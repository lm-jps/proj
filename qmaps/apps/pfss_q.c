/*
 *  pfss_q.c
 *
 *  This module performs smoothing on full-resolution synoptic maps with
 *  3 deg fwhm Gaussian and resample to 361x181, following Toth+ (2011).
 *  PFSS is computed to lmax=120, on a 361x181x50 grid.
 *  Additional Br is saved for 10 layers at 1441x721 where Q will be computed
 *
 *  Author:
 *      Xudong Sun, Based on Xuepu Zhao's Fortran code
 *
 *  Version
 *      v0.0    Apr 18 2017
 *
 *  Notes:
 *      v0.0
 *          Fixed grid, implicitly hardcoded
 *          Smoothing is ad hoc, only appropriate for significant down-sample
 *
 *
 *  Example call:
 *      pfss_q "in=hmi.synoptic_mr_polfil_720s[2099]" "out=hmi_test.pfss_synop"
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "radius.h"
#include "pfss_func.c"

#define PI              (M_PI)
#define TWOPI           (M_PI*2.)
#define HPI             (M_PI/2.)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define NR_PFSS     (ARRLENGTH(r_pfss))
#define NR_OUT      (ARRLENGTH(r_out))

// =================================

// Smoothing synoptic map
void smoothSynop (double *br_o, double *p_o, double *t_o, int np_o, int nt_o, double sig,
                  double *br, double *p, double *t, int np, int nt);

// Gaussian kernel
void createGK (double *gK, int nx, int ny, double sigx, double sigy, double yoff);

// Extract a submap for each row with proper padding
void patchMap (double *br0, int nx0, int ny0, int yc,
               double *br, int nx, int ny);

// Template for output
void writeOutput (DRMS_Record_t *outRec, double *data, int *dims, int ndims, const char *segName);


/* ====================================================================================== */

char *module_name = "pfss_q";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input synoptic map."},
    {ARG_STRING, "out", kNotSpecified, "Output PFSS series."},
    {ARG_DOUBLE, "sig", "1.5", "Sigma of smoothing kernel in degree."},
    {ARG_INT, "np", "361", "Grid points in longitude for input map and PFSS cube."},
    {ARG_INT, "nt", "181", "Grid points in colatitude for input map and PFSS cube."},
    {ARG_INT, "np_q", "1441", "Grid points in longitude for Q-map."},
    {ARG_INT, "nt_q", "721", "Grid points in colatitude for Q-map."},
    {ARG_INT, "lmax", "120", "Maximum L for spherical harmonics."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    const char *inQuery = (const char *) params_get_str(&cmdparams, "in");
    const char *outQuery = (const char *) params_get_str(&cmdparams, "out");
    const double sig = params_get_double(&cmdparams, "sig");
    
    const int np = params_get_int(&cmdparams, "np");
    const int nt = params_get_int(&cmdparams, "nt");
    const int np_q = params_get_int(&cmdparams, "np_q");
    const int nt_q = params_get_int(&cmdparams, "nt_q");
    const int lmax = params_get_int(&cmdparams, "lmax");
    
    const int nr = NR_PFSS, nr_q = NR_OUT;
    
    /* Input Data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs < 1) {
        DIE("Input data series error");
    }
    
    /* Output */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        DIE("Output data series error");
    }

    /* Loop */
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        printf("==========================\n");
        printf("Processing record #%d of %d:\n", irec + 1, nrecs);
        printf("==========================\n");
        fflush(stdout);
        
        // Input/output
        
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *outRec = outRS->records[irec];
        
        // Read input data
        
        DRMS_Segment_t *inSeg = drms_segment_lookupnum(inRec, 0);
        
        DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
        if (status) { SHOW("input error\n"); continue; }
        
        int np_o = inArray->axis[0], nt_o = inArray->axis[1];
        int npix_o = np_o * nt_o;
        double *br_o = (double *) inArray->data;
        
        // Output data
        
        double *br_sm = (double *) (malloc(np * nt * sizeof(double)));
        
        // Grids, can be computed from keywords
        
        double *p_o = (double *) (malloc(np_o * sizeof(double)));       // Input HMI/MDI map at 0.0/0.1...
        double *t_o = (double *) (malloc(nt_o * sizeof(double)));       // Input HMI/MDI map with sine lat
        double crvalx;      // longitude of first column
        for (int idx = 0; idx < np_o; idx++) { p_o[idx] = idx * TWOPI / np_o; }
        for (int idx = 0; idx < nt_o; idx++) { t_o[idx] = HPI - asin((idx + 0.5) * 2. / nt_o - 1.); }
        
        double *p = (double *) (malloc(np * sizeof(double)));           // Output including poles
        double *t = (double *) (malloc(nt * sizeof(double)));
        double *r = (double *) (malloc(nr * sizeof(double)));
        double *kp = (double *) (malloc(np * sizeof(double)));      // area, sum to 1
        double *kt = (double *) (malloc(nt * sizeof(double)));
        for (int idx = 0; idx < np; idx++) {
            p[idx] = idx * TWOPI / (np - 1);
            kp[idx] = (idx == 0 || idx == (np - 1)) ? 0.5 : 1.0;
            kp[idx] /= (np - 1);
        }
        for (int idx = 0; idx < nt; idx++) {
            t[idx] = (nt - 1 - idx) * PI / (nt - 1);
            double y0 = (idx == 0) ? -1. : cos((nt - 0.5 - idx) * PI / (nt - 1));      // south pole
            double y1 = (idx == (nt - 1)) ? 1. : cos((nt - 1.5 - idx) * PI / (nt - 1));    // north pole
            kt[idx] = (y1 - y0) / 2.;
        }
        for (int idx = 0; idx < nr; idx++) { r[idx] = r_pfss[idx]; }
        
        double *p_q = (double *) (malloc(np_q * sizeof(double)));       // Q-map at higher resolution
        double *t_q = (double *) (malloc(nt_q * sizeof(double)));
        double *r_q = (double *) (malloc(nr_q * sizeof(double)));
        for (int idx = 0; idx < np_q; idx++) { p_q[idx] = idx * TWOPI / (np_q - 1); }
        for (int idx = 0; idx < nt_q; idx++) { t_q[idx] = (nt_q - 1 - idx) * PI / (nt_q - 1); }
        for (int idx = 0; idx < nr_q; idx++) { r_q[idx] = r_out[idx]; }
        
        // =============
        // Smoothing
        // =============
        
        SHOW("Smoothing synoptic map... ");
        smoothSynop (br_o, p_o, t_o, np_o, nt_o, sig,
                     br_sm, p, t, np, nt);
        SHOW("done.\n");
        
        // =============
        // Calculate SPH coeffcients
        // =============
        
        double *g = (double *) (calloc((lmax + 1) * (lmax + 1), sizeof(double)));
        double *h = (double *) (calloc((lmax + 1) * (lmax + 1), sizeof(double)));
        
        SHOW("Calculating spherical harmonics... ");
        gh_pfss (br_sm, np, nt, lmax, p, t, kp, kt, g, h);
        SHOW("done.\n");
        
        // =============
        // PFSS cube
        // =============
        
        double *Bp = (double *) (malloc(np * nt * nr * sizeof(double)));
        double *Bt = (double *) (malloc(np * nt * nr * sizeof(double)));
        double *Br = (double *) (malloc(np * nt * nr * sizeof(double)));
        
        SHOW("Calculating PFSS cube... ");
        pfss_cube (g, h, lmax, p, t, r, np, nt, nr, Bp, Bt, Br);
        SHOW("done.\n");
        
        // =============
        // PFSS cube at height for Q
        // =============

        double *Bp_q = (double *) (malloc(np_q * nt_q * nr_q * sizeof(double)));
        double *Bt_q = (double *) (malloc(np_q * nt_q * nr_q * sizeof(double)));
        double *Br_q = (double *) (malloc(np_q * nt_q * nr_q * sizeof(double)));
        
        SHOW("Calculating PFSS cube for Q... ");
        pfss_cube (g, h, lmax, p_q, t_q, r_q, np_q, nt_q, nr_q, Bp_q, Bt_q, Br_q);
        SHOW("done.\n");
        
        // Write out
        
        /*
        int dims0[2] = {np, nt};
        writeOutput (outRec, br_sm, dims0, 2, "Br0");       // arrays freed in subroutine
         */
        free(br_sm);
        
        int dims[3] = {np, nt, nr};
        writeOutput (outRec, Br, dims, 3, "Br");
        writeOutput (outRec, Bp, dims, 3, "Bp");
        writeOutput (outRec, Bt, dims, 3, "Bt");
        
        int dims_q[3] = {np_q, nt_q, nr_q};
        writeOutput (outRec, Br_q, dims_q, 3, "Brq");
        
        // Keywords
        
        int car_rot = drms_getkey_int(inRec, "CAR_ROT", &status);
        double crval1 = drms_getkey_double(inRec, "CRVAL1", &status);
        double crpix1 = drms_getkey_double(inRec, "CRPIX1", &status);
        double cdelt1 = drms_getkey_double(inRec, "CDELT1", &status);
        // New *grid center* of 1s pixel identical with original map
        double crval1_n = crval1 + (1. - crpix1) * cdelt1;
        
        drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
        
        drms_setkey_int(outRec, "CAR_ROT", car_rot);
        
        drms_setkey_int(outRec, "LMAX", lmax);
        drms_setkey_double(outRec, "SIGMA", sig);
        
        drms_setkey_double(outRec, "CRPIX1", 1.);
        printf("crval1=%f, crpix1=%f, cdelt1=%f, crva11_n=%f\n", crval1, crpix1, cdelt1, crval1_n);
        drms_setkey_double(outRec, "CRVAL1", crval1_n);
        drms_setkey_double(outRec, "CDELT1", (-1.) * 360. / (np - 1.));
        drms_setkey_string(outRec, "CUNIT1", "degree");
        drms_setkey_string(outRec, "CTYPE1", "CRLN-CAR");
        
        drms_setkey_double(outRec, "CRPIX2", (1. + nt) / 2.);
        drms_setkey_double(outRec, "CRVAL2", 0.0);
        drms_setkey_double(outRec, "CDELT2", 180. / (nt - 1.));
        drms_setkey_string(outRec, "CUNIT2", "degree");
        drms_setkey_string(outRec, "CTYPE2", "CRLT-CAR");
        
        TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(outRec, "DATE", tnow);

        // Clean up
        
        free(p_o); free(t_o);
        free(p); free(t); free(r); free(kp); free(kt);
        free(p_q); free(t_q); free(r_q);
        free(g); free(h);
        free(Bt_q); free(Bp_q);
        drms_free_array(inArray);
        
    }
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);

    //
    
    return DRMS_SUCCESS;
    
}       // DoIt


/*
 * Smoothing Br synopitc map with specified Gaussian kernel sigma=sig
 *
 * Input longitude assumed at 0, 2pi*1/np_o, 2pi*2/np_o
 * Input latitude assumed at asin(0.5/(nt_o/2))-pi, etc. (sine lat)
 * Output longitude assumed at 0, 2pi*1/(np-1), etc., with
 *  last column at 360., repeating the first
 * Output latitude assumed at 0, 2pi*1/(nt-1)-pi/2, etc., with
 *  last row at 180., both south/north poles are covered
 *
 * Naive implementation:
 * The poles are always averaging of the last row
 * For desired output lon/lat, find out correspoding sigmax, sigmay
 * Find out the center pixel address to nearest full pixel
 * Then create a Gaussian kernel with [sigmax,sigmay], convolve row by row
 * For edge in X direction, wrap around
 * For edge in Y direction, repeat
 *
 */

void smoothSynop (double *br_o, double *p_o, double *t_o, int np_o, int nt_o, double sig,
                  double *br, double *p, double *t, int np, int nt)
{
    
    double sig_rad = sig * RADSINDEG;
    
    // Grid parameters
    
    double dx_o = TWOPI / np_o;
    double dy_o = ((double) 2.0) / nt_o;        // Input, assuming uniform lon and sinlat
    
    // Find out sigmax & sigmay, constant for same row
    
    double *sigmax = (double *) (malloc(nt * sizeof(double)));
    double *sigmay = (double *) (malloc(nt * sizeof(double)));
    double tmin, tmax;      // lower/upper colat
    for (int idx = 0; idx < nt; idx++) {            // in pixel of original map
        sigmax[idx] = sig_rad / dx_o / sin(t[idx]);
        tmin = MIN(MAX(t[idx] - sig_rad, 0.), PI);
        tmax = MAX(MIN(t[idx] + sig_rad, PI), 0.);
        sigmay[idx] = MAX((cos(tmin) - cos(t[idx])) / dy_o, (cos(t[idx]) - cos(tmax)) / dy_o);
    }
    sigmax[0] = sigmax[nt - 1] = np_o / 2.;         // poles
    
    // Find out kernel size, from filter_image.pro
    
    int *npixx = (int *) (malloc(nt * sizeof(int)));
    int *npixy = (int *) (malloc(nt * sizeof(int)));
    for (int idx = 0; idx < nt; idx++) {
        npixx[idx] = 2 * floor(MIN(MAX(3. * 2.355 * sigmax[idx], 3.), np_o) / 2.) + 1;
        npixy[idx] = 2 * floor(MIN(MAX(3. * 2.355 * sigmay[idx], 3.), nt_o) / 2.) + 1;
    }
    
    // Find out center pixel in original map, with y offset
    
    int *x_c = (int *) (malloc(np * sizeof(int)));
    int *y_c = (int *) (malloc(nt * sizeof(int)));
    double *y_off = (double *) (malloc(nt * sizeof(double)));
    for (int idx = 0; idx < np; idx++) {
        x_c[idx] = ((int) (round((p[idx] - p_o[0]) / dx_o))) % np_o;
    }
    for (int idx = 0; idx < nt; idx++) {
        y_c[idx] = round((cos(t[idx]) - cos(t_o[0])) / dy_o);
        y_off[idx] = (cos(t[idx]) - cos(t_o[0])) / dy_o - y_c[idx];
    }
    
    // Do it row by row
    
    for (int row = 0; row < nt; row++) {
        
        // Create kernel
        
        int npix = npixx[row] * npixy[row];
        double *gK = (double *) (malloc(npix * sizeof(double)));
        
        createGK (gK, npixx[row], npixy[row], sigmax[row], sigmay[row], y_off[row]);
        
        // Create a submap for each row with proper padding
        
        int x_patch = npixx[row] / 2;
        int np_t = np_o + x_patch * 2, nt_t = npixy[row], npix_t = np_t * nt_t;
        double *br_t = (double *) (malloc(npix_t * sizeof(double)));
        
        patchMap (br_o, np_o, nt_o, y_c[row],
                  br_t, np_t, nt_t);
        
        // Do it for each pixel in the row
        
        for (int col = 0; col < np; col++) {
            
            double sum = 0.;
            int x_lead = x_patch + x_c[col] - (npixx[row] / 2);      // starting X coord
            
            // Loop over each row/col in br_t
            
            int idx, idx_t;
            for (int y = 0; y < npixy[row]; y++) {
                for (int x = 0; x < npixx[row]; x++) {
                    idx = y * npixx[row] + x;       // index in kernel
                    idx_t = y * np_t + (x_lead + x);
                    sum += (br_t[idx_t] * gK[idx]);
                }
            }
            
            br[row * np + col] = sum;
            
        }   // col
        
        // Clean up
        
        // printf("row #%d done.\n", row); fflush(stdout);
        free(gK); free(br_t);
        
    }   // row

    // Test printing
    // ============================
    /*
     for (int idx = 0; idx < 10; idx++) { printf("lon_o[%d]=%f\n", idx, p_o[idx]/RADSINDEG); }
     for (int idx = 0; idx < 10; idx++) { printf("lat_o[%d]=%f\n", idx, 90. - t_o[idx]/RADSINDEG); }
     */
    /*
     for (int idx = 0; idx < 10; idx++) { printf("lon[%d]=%f\n", idx, p[idx]/RADSINDEG); }
     for (int idx = 0; idx < 10; idx++) { printf("lat[%d]=%f\n", idx, 90. - t[idx]/RADSINDEG); }
     */
    /*
     for (int idx = 0; idx < 11; idx++) { printf("sigmax[%d]=%f\n", idx, sigmax[idx]); }
     for (int idx = 0; idx < 11; idx++) { printf("sigmay[%d]=%f\n", idx, sigmay[idx]); }
     for (int idx = 80; idx < 91; idx++) { printf("sigmax[%d]=%f\n", idx, sigmax[idx]); }
     for (int idx = 80; idx < 91; idx++) { printf("sigmay[%d]=%f\n", idx, sigmay[idx]); }
     */
    /*
     for (int idx = 0; idx < 11; idx++) { printf("npixx[%d]=%d\n", idx, npixx[idx]); }
     for (int idx = 0; idx < 11; idx++) { printf("npixy[%d]=%d\n", idx, npixy[idx]); }
     for (int idx = 80; idx < 91; idx++) { printf("npixx[%d]=%d\n", idx, npixx[idx]); }
     for (int idx = 80; idx < 91; idx++) { printf("npixy[%d]=%d\n", idx, npixy[idx]); }
     */
    /*
     for (int idx = 357; idx < 361; idx++) { printf("x_c[%d]=%d\n", idx, x_c[idx]); }
     for (int idx = 1; idx < 11; idx++) { printf("y_c[%d]=%d, y_off[%d]=%f\n", idx, y_c[idx], idx, y_off[idx]); }
     for (int idx = 80; idx < 91; idx++) { printf("y_c[%d]=%d, y_off[%d]=%f\n", idx, y_c[idx], idx, y_off[idx]); }
     */
    /*
     int npix = 11 * 11;
     double *gK = (double *) (malloc(npix * sizeof(double)));
     createGK (gK, 11, 11, 1.5, 1.5, -0.3);
     int idx = 0;
     double sum = 0.0;
     for (int row = 0; row < 11; row++) {
     for (int col = 0; col < 11; col++) {
     printf("%9.6f", gK[idx]);
     sum += gK[idx++];
     }
     printf("\n");
     }
     printf("sum=%f\n", sum);
     free(gK);
     */
    // ================================
    
    // Clean up
    
    free(sigmax); free(sigmay); free(npixx); free(npixy);
    free(x_c); free(y_c); free(y_off);
    
}



/*
 * Create a normalized Gaussian kernel
 * Consider subpixel offset in y direction
 * as nrows is generally large, additional pixel on one side
 * shouldn't matter that much
 *
 */

void createGK (double *gK, int nx, int ny, double sigx, double sigy, double yoff)
{
    
    // 1D X
    
    double *psfx = (double *) (malloc(nx * sizeof(double)));
    double xc = (nx - 1.) / 2.0;
    double xfac = 1.0 / (sqrt(TWOPI) * sigx);
    double x;
    for (int idx = 0; idx < nx; idx++) {
        x = (idx - xc) / sigx;
        psfx[idx] = exp(- x * x / 2.) * xfac;
//        printf("psfx[%d]=%f\n", idx, psfx[idx]);
    }

    // 1D Y
    
    double *psfy = (double *) (malloc(ny * sizeof(double)));
    double yc = (ny - 1.) / 2.0 + yoff;
    double yfac = 1.0 / (sqrt(TWOPI) * sigy);
    double y;
    for (int idx = 0; idx < ny; idx++) {
        y = (idx - yc) / sigy;
        psfy[idx] = exp(- y * y / 2.) * yfac;
//        printf("psfy[%d]=%f\n", idx, psfy[idx]);
    }
    
    // 2D
    
    double sum = 0.;
    int idx = 0;
    for (int row = 0; row < ny; row++) {
        for (int col = 0; col < nx; col++) {
            gK[idx] = psfx[col] * psfy[row];
            sum += gK[idx++];
        }
    }
    
    // Normalize
    
    int npix = nx * ny;
    for (int idx = 0; idx < npix; idx++) {
        gK[idx] /= sum;
    }
    
    // Clean up
    
    free(psfx); free(psfy);
}


/*
 * Extract a submap of size (nx, ny), centered at
 * (xc, yc), with propper padding
 *
 */

void patchMap (double *br0, int nx0, int ny0, int yc,
               double *br, int nx, int ny)
{
    // No error checking now
    
    int x_patch = (nx - nx0) / 2;
    int ymin = yc - ny / 2;
    
    // Copy row by row, central portion
    
    for (int y = 0; y < ny; y++) {
        
        int y0 = y + ymin;      // row number in original map
        int idx0;               // starting index in original map
        
        if (y0 < 0) {                       // below south pole
            idx0 = 0;
        } else if (y0 > (ny0 - 1)) {        // above north pole
            idx0 = (ny0 - 1) * nx0;
        } else {                            // regular
            idx0 = y0 * nx0;
        }
        
        int idx = y * nx + x_patch;      // starting index in new map
        
        for (int x = 0; x < nx0; x++) {     // copy original map
            br[idx + x] = br0[idx0 + x];
        }
        
    }
    
    // Padding in x direction, periodically
    
    for (int y = 0; y < ny; y++) {
        
        int idx = y * nx;
        
        for (int x = 0; x < x_patch; x++) {
            br[idx + x] = br[idx + x + nx0];
            br[idx + x_patch + nx0 + x] = br[idx + x_patch + x];
        }
        
    }
    
}


/*
 * Template for output, no error checking
 *
 */

void writeOutput (DRMS_Record_t *outRec, double *data, int *dims, int ndims, const char *segName)
{

    int status = 0;
    DRMS_Segment_t *outSeg = drms_segment_lookup(outRec, segName);
    DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_DOUBLE, ndims, dims, data, &status);
    for (int idx = 0; idx < ndims; idx++) {
        outSeg->axis[idx] = outArray->axis[idx];
    }
    outArray->israw = 0;		// compressed
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    status = drms_segment_write(outSeg, outArray, 0);
    drms_free_array(outArray);
    
}

