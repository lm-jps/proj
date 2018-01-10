/*
 *  cgem_map.c
 *
 *	This module maps a scaler map of one projection to another
 *  Written as prototype for functions needed for cgem_prep
 *  One application is to map a Plate-Carree Br map into Mercater and back
 *
 *	Author:
 *		Xudong Sun
 *
 *	Version:
 *              v0.0 Jul 18 2016
 *	Notes:
 *		v0.0
 *
 *	Example Calls:
 *      > cgem_map "in=su_xudong.B_720s_carree[377][2011.02.15_00:00]" "seg=Bz" "out=su_xudong.B_720s_mercator" "map=Mercator" "cols=600" "rows=428" -a
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
#include "cartography_cgem.c"
#include "finterpolate.h"

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

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
#define dpath    "/home/jsoc/cvs/Development/JSOC"

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

// Ephemeris information
struct ephemeris {
    double disk_lonc, disk_latc;
    double disk_xc, disk_yc;
    double rSun, asd, pa;
};

// Mapping information
struct reqInfo {
    double xc, yc;
    int ncol, nrow;
    double dx, dy;
    int proj;
};

/* ========================================================================================================== */

char *mapName[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
    "LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
    "stereographic", "orthographic", "LambertAzimuthal"};
char *wcsCode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
    "SIN", "ZEA"};

/* ========================================================================================================== */

/* ========================================================================================================== */

char *module_name = "cgem_map";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input series."},
    {ARG_STRING, "out", kNotSpecified, "Output series."},
    {ARG_STRING, "seg", kNotSpecified, "Iutput segment."},
    {ARG_NUME, "map", "carree", "Projetion method, carree by default.",
        "carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
    {ARG_FLOAT, "xref", "0", "Reference output map x center, ususally Carrington lon."},
    {ARG_FLOAT, "yref", "0", "Reference output map y center, ususally Carrington lon."},
    {ARG_FLOAT, "dx", "0", "X pixel size, unit depending on projection (default deg)."},
    {ARG_FLOAT, "dy", "0", "Y pixel size, unit depending on projection (default deg)."},
    {ARG_FLAG, "a", "", "Automatic, i.e. adopting original center and scale, override xyref, dxy."},
    {ARG_INT, "cols", "500", "Columns of output."},
    {ARG_INT, "rows", "500", "Rows of output."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    char *inQuery = NULL, *segQuery = NULL;      // Input query
    char *outQuery = NULL;                      // Output query
    
    inQuery = (char *) params_get_str(&cmdparams, "in");
    segQuery = (char *) params_get_str(&cmdparams, "seg");
    outQuery = (char *) params_get_str(&cmdparams, "out");
    
    if (!strcmp(inQuery, kNotSpecified) ||
        !strcmp(segQuery, kNotSpecified) ||
        !strcmp(outQuery, kNotSpecified)) {
        DIE("Input/output not available");
    }
    
    /* Get arguments */
    
    int automap = params_isflagset(&cmdparams, "a");
    
    struct reqInfo req;
    req.xc = params_get_float(&cmdparams, "xref") * RADSINDEG;
    req.yc = params_get_float(&cmdparams, "yref") * RADSINDEG;
    req.ncol = params_get_int(&cmdparams, "cols");
    req.nrow = params_get_int(&cmdparams, "rows");
    req.dx = params_get_float(&cmdparams, "dx") * RADSINDEG;		// deg to rad, for now
    req.dy = params_get_float(&cmdparams, "dy") * RADSINDEG;
    req.proj = params_get_int(&cmdparams, "map");
    
    if (req.ncol <= 2 || req.ncol > 4000 ||
        req.nrow <= 2 || req.nrow > 4000) {
        DIE("requested map too small/too big");
    }

    /* Input Data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs == 0 || !inRS) DIE("Input data series error");
    DRMS_Segment_t *inSeg = drms_segment_lookup(inRS->records[0], segQuery);        // check segment
    if (!inSeg) {
        drms_close_records(inRS, DRMS_FREE_RECORD);
        DIE("Requested segment not available");
    }
    
    /* Output Data */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status || !outRS) {
        drms_close_records(outRS, DRMS_FREE_RECORD);
        DIE("Error creating output series");
    }
    
    /* Loop */
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *outRec = outRS->records[irec];
        printf("Processing rec #%d of %d:\n", irec+1, nrecs);
        
        // Determine center, scale, xysize, etc.
        
        double xc_o, yc_o, dx_o, dy_o, colc_o, rowc_o;  // original
        double crln_obs;
        int ncol_o, nrow_o, proj_o;
        char proj_str[4], *ctype1;
        double xc, yc, dx, dy;          // requested
        int ncol, nrow, proj;
        
        // input projection, no error checking for now
        xc_o = drms_getkey_double(inRec, "CRVAL1", &status) * RADSINDEG;        // lonc
        yc_o = drms_getkey_double(inRec, "CRVAL2", &status) * RADSINDEG;        // latc
        crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status) * RADSINDEG;
        xc_o -= crln_obs;       // Stonyhurst
        colc_o = drms_getkey_double(inRec, "CRPIX1", &status) - 1;
        rowc_o = drms_getkey_double(inRec, "CRPIX2", &status) - 1;
        dx_o = drms_getkey_double(inRec, "CDELT1", &status) * RADSINDEG;
        dy_o = drms_getkey_double(inRec, "CDELT2", &status) * RADSINDEG;
        ctype1 = drms_getkey_string(inRec, "CTYPE1", &status);
        strncpy(proj_str, ctype1+5, 3);
        for (proj_o = 0; proj_o++; proj_o < ARRLENGTH(wcsCode)) {
            if (strcmp(proj_str, wcsCode[proj_o])== 0) break;
        }
        if (proj_o >= ARRLENGTH(wcsCode) || status) {
            printf("Erroneous geometry, image skipped.\n");
            continue;
        } else {
        	proj_o -= 1;
        }
        if (automap) {
            xc = xc_o; yc = yc_o;
            dx = dx_o; dy = dy_o;
            // determin maximum size
            ;
        } else {
            xc = req.xc - crln_obs; yc = req.yc;
            dx = req.dx; dy = req.dy;
        }
        ncol = req.ncol; nrow = req.nrow;
        proj = req.proj;
        printf("xc_o=%f, yc_o=%f, dx_o=%f, dy_o=%f, proj_o=%s\n",
               xc_o/RADSINDEG, yc_o/RADSINDEG, dx_o/RADSINDEG, dy_o/RADSINDEG, mapName[proj_o]);
        printf("xc=%f, yc=%f, dx=%f, dy=%f, proj=%s\n",
               xc/RADSINDEG, yc/RADSINDEG, dx/RADSINDEG, dy/RADSINDEG, mapName[proj]);
        
        // Input image
        
        inSeg = drms_segment_lookup(inRec, segQuery);
        DRMS_Array_t *inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status) return 1;
        float *map_in = (float *)inArray->data;
        ncol_o = inSeg->axis[0]; nrow_o = inSeg->axis[1];
        
        // Coordinates
        
        int npix = ncol * nrow;
        float *col_o = (float *) (malloc(npix * sizeof(float)));
        float *row_o = (float *) (malloc(npix * sizeof(float)));
        double rowc = (nrow - 1.) / 2., colc = (ncol - 1.) / 2.;
        printf("colc=%f, rowc=%f, colc_o=%f, rowc_o=%f\n", colc, rowc, colc_o, rowc_o);
        
        // Convert
        
        printf("Mapping..."); fflush(stdout);
        double x, y, x_o, y_o, lon, lat;
        int idx;
        for (int row = 0; row < nrow; row++) {
            y = (row - rowc) * dy;
            for (int col = 0; col < ncol; col++) {
                x = (col - colc) * dx;
                idx = row * ncol + col;
                if (plane2sphere(x, y, yc, xc, &lat, &lon, proj)) {     // requested
                    col_o[idx] = DRMS_MISSING_FLOAT;
                    row_o[idx] = DRMS_MISSING_FLOAT;
                    continue;
                }
                if (sphere2plane(lat, lon, yc_o, xc_o, &x_o, &y_o, proj_o)) {       // orig
                    col_o[idx] = DRMS_MISSING_FLOAT;
                    row_o[idx] = DRMS_MISSING_FLOAT;
                    continue;
                }
                col_o[idx] = x_o / dx_o + colc_o;      // index in input images
                row_o[idx] = y_o / dy_o + rowc_o;
            }
        }
        printf(" done\n"); fflush(stdout);
        
        // Interpolate
        
        printf("Interpolate..."); fflush(stdout);
        float *map_out = (float *) (malloc(npix * sizeof(float)));      // output array
        struct fint_struct pars;
        init_finterpolate_linear(&pars, 0.);        // bi-linear
        
        finterpolate(&pars, map_in, col_o, row_o, map_out,
                     ncol_o, nrow_o, ncol_o, ncol, nrow, ncol, DRMS_MISSING_FLOAT);
        
        free(col_o); free(row_o);       // Clean up
        drms_free_array(inArray);
        printf(" done\n"); fflush(stdout);
        
        // Output
        
         printf("Write data..."); fflush(stdout);
        int outDims[2] = {ncol, nrow};
        DRMS_Segment_t *outSeg = drms_segment_lookupnum(outRec, 0);
        DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, map_out, &status);
        if (status) {
            SHOW("Output array error");
            free(map_out);
            continue;
        }
        outSeg->axis[0] = outArray->axis[0]; outSeg->axis[1] = outArray->axis[1];
        outArray->israw = 0;		// always compressed
        outArray->bzero = outSeg->bzero;
        outArray->bscale = outSeg->bscale;
        
        status = drms_segment_write(outSeg, outArray, 0);
        drms_free_array(outArray);
        if (status) {
            SHOW("Output array write error");
            continue;
        }
        
        // Keywords
        
        drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);     // copy all keys
        
        TIME val, trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(outRec, "DATE", tnow);

        drms_setkey_float(outRec, "CRPIX1", ncol/2. + 0.5);
        drms_setkey_float(outRec, "CRPIX2", nrow/2. + 0.5);
        
        drms_setkey_float(outRec, "CRVAL1", (xc + crln_obs) / RADSINDEG);
        drms_setkey_float(outRec, "CRVAL2", yc / RADSINDEG);
        drms_setkey_float(outRec, "CDELT1", dx / RADSINDEG);
        drms_setkey_float(outRec, "CDELT2", dy / RADSINDEG);
        drms_setkey_string(outRec, "CUNIT1", "degree");
        drms_setkey_string(outRec, "CUNIT2", "degree");
        
        char key[64];
        snprintf (key, 64, "CRLN-%s", wcsCode[proj]);
        drms_setkey_string(outRec, "CTYPE1", key);
        snprintf (key, 64, "CRLT-%s", wcsCode[proj]);
        drms_setkey_string(outRec, "CTYPE2", key);
        drms_setkey_float(outRec, "CROTA2", 0.0);
        
        printf(" done\n"); fflush(stdout);
        
    }       // irec
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    return DRMS_SUCCESS;
    
}	// DoIt

