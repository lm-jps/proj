/*
 * Module name:		jpolfil.c
 *
 * Description:		This module takes a synoptic map and fill in polar regions, using
 *			certain interpolation algorithms. By default, the polar region is
 *			remapped and a least square fitting is used to determine the missing
 *			points. See notes below for detailed options.
 *
 * Original source:	Image fitting code by Rasmus Larsen (rmlarsen@gmail.com)
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Oct 01 2009
 *			v1.1		Oct 05 2009
 *			v1.2		Mar 13 2010
 *			v1.3		Mar 16 2010
 *			v1.4		Mar 26 2010
 *			v1.5		Mar 30 2010
 *      v1.6    Nov 13 2011
 *
 * Issues:
 *			v1.0
 *			Only default option is implemented
 *			Both poles are filled above a certain latitude
 *			(could implement it so only NaNs are filled)
 *			v1.1
 *			Added pvflag & ftflag intermediate output data series
 *			Default pv series: su_xudong.polarview
 *			Default ft series: su_xudong.polarfit
 *			Stereographic projection project to z=1 rather than z=0
 *			v1.2
 *			Removed pv and ft. Original version in jpolfil0.c
 *			v1.3
 *			Moved everything from polarproj.c and fitimage.c to here
 *			Options, flags all made global so don't have to be passed around
 *			v1.4
 *			Added/modify some keywords according to su_yang.hmi_synop
 *			v1.5
 *			Moved fitimage and timing out
 *			Added a link to original synoptic map (subject to change)
 *      v1.6
 *      Updated keywords list for official hmi synoptic map
 *
 * Example:		jpolfil "in=hmi.Synoptic_Mr_720s[1928]" out=hmi_test.Synoptic_Mr_polfil_720s
 *
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define	DTOR	(M_PI / 180.)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

// For time measuring, by T. Larson
#include "timing.c"

// For 2D image fitting, by R. Larsen
#include "fitimage.c"

// For stats, by K. Chu
#include "fstats.c"


/* ############## Structures, Global Variables, Prototypes ############## */

#define PROJNM	STEREOGRAPHIC
#define INTPNM	BILINEAR
#define IMGDIM	600

/* Filling option flags */
struct fil_opt {
    double lat0;
    double latfil;
    int nflag;
    int sflag;
};

/* Method numerator */
enum methodnums {
    SPATIAL_2D
};
char *methodlist[] = {"SPATIAL_2D"};

/* Interpolation enumerators */
enum intp {
    BILINEAR
};
char *intplist[] = {"BILINEAR"};

/* Projection enumerators */
enum proj {
    STEREOGRAPHIC
};
char *projlist[] = {"STEREOGRAPHIC"};


/* Global variables */
struct fil_opt option;
enum methodnums methodnum;
enum intp intpname;
enum proj projname;
int verbflag;


/* 2D spatial fitting */
void spatial_2d(float *synop, double *lons, double *lats, 
                int nlons, int nlats, int sinlat, int order, int imgDims[]);

/* Polar fitting, from polarproj.c */
double get_index(int sinlat, int n, double x[], double x0);
double linintd(double *f, int nx, int ny, double x, double y);

void stereo_fwd(int south, double lon, double lat, double *x, double *y, 
                double xc, double yc, double rx, double ry, int nx, int ny);
void stereo_bwd(int south, double *lon, double *lat, double x, double y, 
                double xc, double yc, double rx, double ry, int nx, int ny);

void polar_proj_fwd(int south, int sinlat, double *image, int nx, int ny,
                    double *map, int nlons, int nlats, double lons[], double lats[], double lat_lim);
void polar_proj_bwd(int south, int sinlat, double *image, int nx, int ny,
                    double *map, int nlons, int nlats, double lons[], double lats[], double lat_lim);




/* ################## Main Module ################## */


char *module_name = "jpolfil";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL, "Output data series."},
    {ARG_STRING, "METHOD", "SPATIAL_2D", "Filling method."},
    {ARG_STRING, "INTERP", "BILINEAR", "Interpolation method."},
    {ARG_STRING, "PROJECT", "STEREOGRAPHIC", "Map projection method."},
    {ARG_DOUBLE, "LAT0", "65.", "Starting latitude used for fitting."},
    {ARG_DOUBLE, "LATFIL", "75.", "Starting latitude for filling in."},
    {ARG_INT, "NFLAG", "1", "Filling for north pole: 1 for yes, 0 for no."},
    {ARG_INT, "SFLAG", "1", "Filling for south pole: 1 for yes, 0 for no."},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=all messages"},
    {ARG_INT, "ORDER", "7", "Highest power of polynomials."},
    {ARG_END}
};

/* Keywords to copy over */
char *synop_keys[] = {"TELESCOP", "INSTRUME", "WAVELNTH", "BUNIT", "CONTENT", "HISTORY", "COMMENT", 
                      "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", 
                      "CUNIT1", "CUNIT2", "CRDER1", "CRDER2", "CSYSER1", "CSYSER2", "WCSNAME", 
                      "T_REC", "T_REC_epoch", "T_REC_step", "T_REC_unit", "CADENCE", "DATASIGN", "T_OBS", "T_START", "T_STOP", "T_ROT", "T_EARTH", 
                      "CARRTIME", "B0_ROT", "B0_FRST", "B0_LAST", "EARTH_B0", "LON_FRST", "LON_LAST", "LON_STEP", "W_OFFSET", "W_WEIGHT", 
                      "MAG_NUM", "MAG_FRST", "MAG_LAST", "MAG_ROT", "HWNWIDTH", "EQPOINTS", "NSIGMA", "CARSTRCH", "DIFROT_A", "DIFROT_B", "DIFROT_C",
                      "DSUN_OBS", "CRLN_OBS", "CRLT_OBS", "OBS_VR", "OBS_VN", "OBS_VW", "QUALITY", "MAPMMAX", "MAPLGMAX", "MAPLGMIN", "SINBDIVS",
                      "MAPBMAX", "MAPRMAX", "LGSHIFT", "INTERPO", "MCORLEV", "MOFFSET", "FRTIMWDN", "FRAVEPNT", "FRWINDOW", "SYNDRORA",
                      "SYNDRO_A", "SYNDRO_B", "SYNDRO_C", "OSYNDR_A", "OSYNDR_B", "OSYNDR_C"};

int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS, *outRS;
    int irec, nrecs;
    DRMS_Record_t *inRec, *outRec;
    DRMS_Segment_t *inSeg, *outSeg;
    DRMS_Array_t *inArray, *outArray;
    DRMS_Link_t *srcLink;
    float *inData, *outData;

    char *inchecklist[] = {"CAR_ROT", "CTYPE2"};
    char *outchecklist[] = {"CAR_ROT", "METHOD", "LAT0", "LATFIL", "NFLAG", "SFLAG"};
    DRMS_Record_t *inRec_tmp, *outRec_tmp;
    DRMS_Keyword_t *inkeytest, *outkeytest;

    int outDims[2];
    char *method, *interp, *project;
    int methodflag = 0, intpflag = 0, projflag = 0;
    int imgDims[3];	// Both north and south

    char *ctype2;
    int sinlat, order;
    double lon0;

    double *lons, *lats;

    int i, j, itest;
    int n0, n1, mapsz;
    int img_n, img_s, imgsz;

    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);
    
    // For stats
    double dat_min, dat_max, dat_medn, dat_mean, dat_sig, dat_skew, dat_kurt;
    int dat_ngood;

    /* Get parameters */
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    method = (char *)params_get_str(&cmdparams, "METHOD");
    for (itest = 0; itest < ARRLENGTH(methodlist) && !methodflag; itest++) {
        if (!strcmp(method, methodlist[itest])) {
            methodflag = 1;
            methodnum = (enum methodnums)(itest);
        }
    }
    if (!methodflag) DIE("Unknow method");
    
    interp = (char *)params_get_str(&cmdparams, "INTERP");
    for (itest = 0; itest < ARRLENGTH(intplist) && !intpflag; itest++) {
        if (!strcmp(interp, intplist[itest])) {
            intpflag = 1;
            intpname = (enum intp)(itest);
        }
    }
    if (!intpflag) DIE("Unknow interpolation method");
    
    project = (char *)params_get_str(&cmdparams, "PROJECT");
    for (itest = 0; itest < ARRLENGTH(projlist) && !projflag; itest++) {
        if (!strcmp(project, projlist[itest])) {
            projflag = 1;
            projname = (enum proj)(itest);
        }
    }
    if (!projflag) DIE("Unknow projection method");

    // For 2D spatial fitting
    option.lat0 = params_get_double(&cmdparams, "LAT0");
    option.latfil = params_get_double(&cmdparams, "LATFIL");
    option.nflag = params_get_int(&cmdparams, "NFLAG");
    option.sflag = params_get_int(&cmdparams, "SFLAG");
    verbflag = params_get_int(&cmdparams, "VERB");
    order = params_get_int(&cmdparams, "ORDER");

    /* Open input */
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;

    /* Check for essential input keywords */
    inRec_tmp = inRS->records[0];
    for (itest = 0; itest < ARRLENGTH(inchecklist); itest++) {
        inkeytest = hcon_lookup_lower(&inRec_tmp->keywords, inchecklist[itest]);
        if (!inkeytest) {
            fprintf(stderr, "ERROR: required input keyword %s is missing\n", inchecklist[itest]);
            return 0;
        }
    }

    /* Check for essential output keywords */
    outRec_tmp = drms_create_record(drms_env, outQuery, DRMS_TRANSIENT, &status);
    if (status) DIE("Unable to open record in output series");
    for (itest = 0; itest < ARRLENGTH(outchecklist); itest++) {
        outkeytest = hcon_lookup_lower(&outRec_tmp->keywords, outchecklist[itest]);
        if (!outkeytest || outkeytest->info->islink || outkeytest->info->recscope == 1) {
            fprintf(stderr, "Output keyword %s is missing/constant/a link.\n", outchecklist[itest]);
            return 0;
        }
    }
    drms_close_record(outRec_tmp, DRMS_FREE_RECORD);

    /* Create output */
    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output recordset not created");

    /* Do this for each record */
    for (irec = 0; irec < nrecs; irec++)
    {
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }

        /* Input record and data */
        inRec = inRS->records[irec];
        inSeg = drms_segment_lookupnum(inRec, 0);	// Single segment
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status) {
            printf("No data file found. \n"); fflush(stdout);
            drms_free_array(inArray);
            continue;
        }
        inData = (float *)inArray->data;

        /* Parameters */
        n0 = inArray->axis[0]; n1 = inArray->axis[1];
        mapsz = n0 * n1;
        ctype2 = drms_getkey_string(inRec, "CTYPE2", &status);
        sinlat = (strcmp(ctype2, "Sine Latitude") == 0 || strcmp(ctype2, "CRLT-CEA") == 0) ? 1 : 0;

        /* Output data */
        outDims[0] = n0; outDims[1] = n1;
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        outData = (float *)outArray->data;
        for (i = 0; i < mapsz; i++) {
            outData[i] = inData[i];	// Copy data (including NAN)
        }
        drms_free_array(inArray);

        // For 2D spatial fitting
        imgDims[0] = IMGDIM; imgDims[1] = IMGDIM; imgDims[2] = 2;

        /* Map geometry */
        lons = (double *)malloc(n0 * sizeof(double));
        lats = (double *)malloc(n1 * sizeof(double));
        for (i = 0; i < n0; i++) lons[i] = (i - n0 / 2.0) / n0 * 360. * DTOR;
        for (j = 0; j < n1; j++)
            lats[j] = (j + 0.5 - n1 / 2.0) / (n1 / 2.0);
        if (sinlat) {
            for (j = 0; j < n1; j++) lats[j] = asin(lats[j]);
        } else {
            for (j = 0; j < n1; j++) lats[j] = lats[j] * M_PI;
        }

        /* This is the working part */
        switch (methodnum) {
            case ((enum methodnums)(0)):
                spatial_2d(outData, lons, lats, n0, n1, sinlat, order, imgDims);
                break;
            default:
                spatial_2d(outData, lons, lats, n0, n1, sinlat, order, imgDims);
                break;
        }

        /* Outpur record */
        outRec = outRS->records[irec];
        outSeg = drms_segment_lookupnum(outRec, 0);
        outSeg->axis[0] = outArray->axis[0];
        outSeg->axis[1] = outArray->axis[1];
        outArray->parent_segment = outSeg;
        
        /* Set link */
        srcLink = hcon_lookup_lower(&outRec->links, "SYNOPMR_ORIG");
        if (srcLink) {
            drms_setlink_static(outRec, "SYNOPMR_ORIG", inRec->recnum);
        }

        /* Set keywords */
        // Prime key
        drms_copykey(outRec, inRec, "CAR_ROT");
        // Methods
        drms_setkey_string(outRec, "METHOD", method);
        drms_setkey_double(outRec, "LAT0", option.lat0);
        drms_setkey_double(outRec, "LATFIL", option.latfil);
        drms_setkey_int(outRec, "NFLAG", option.nflag);
        drms_setkey_int(outRec, "SFLAG", option.sflag);

        /* Stats */
        status = fstats(mapsz, outData, &dat_min, &dat_max, &dat_medn, &dat_mean, 
        				&dat_sig, &dat_skew, &dat_kurt, &dat_ngood);
        
        /* Changed keywords */
        drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        drms_setkey_int(outRec, "TOTVALS", mapsz);
        if (!status) {
	        drms_setkey_int(outRec, "DATAVALS", dat_ngood);
    	    drms_setkey_int(outRec, "MISSVALS", mapsz - dat_ngood);
    	    drms_setkey_double(outRec, "DATAMIN", dat_min);
    	    drms_setkey_double(outRec, "DATAMAX", dat_max);
    	    drms_setkey_double(outRec, "DATAMEDN", dat_medn);
    	    drms_setkey_double(outRec, "DATAMEAN", dat_mean);
    	    drms_setkey_double(outRec, "DATARMS", dat_sig);
    	    drms_setkey_double(outRec, "DATASKEW", dat_skew);
    	    drms_setkey_double(outRec, "DATAKURT", dat_kurt);
    	}
    	drms_copykey(outRec, inRec, "T_REC");
    	
    	/* Result writing */
        status = drms_segment_write(outSeg, outArray, 0);
        if (status) DIE("problem writing data");
        drms_free_array(outArray);
        free(lons); free(lats);
    	
        // Copy over
        for (itest = 0; itest < ARRLENGTH(synop_keys); itest++) {
            drms_copykey(outRec, inRec, synop_keys[itest]);
        }

        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                     irec, wt - wt1, ct - ct1);
        }
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

}  // End of module








/* ################## 2D Spatial Fitting ################## */


/* Retuns a filled map using 2d spatial fitting */
void spatial_2d(synop, lons, lats, nlons, nlats, sinlat, order, imgDims)
    float *synop;
    double *lons, *lats;
    int nlons, nlats, sinlat, order;
    int imgDims[];
{
    int i, j;
    int dj;	// rows of filling
    int mapsz = nlons * nlats;
    int imgsz = imgDims[0] * imgDims[1];
    double *lats_n, *lats_s; 
    double *map;
    double *image_n, *image_s;

    /* Additional geometry */
    lats_n = (double *)malloc(nlats * sizeof(double));
    lats_s = (double *)malloc(nlats * sizeof(double));
    for (j = 0; j < nlats; j++) {
        lats_n[j] = (j + 1. - nlats / 2.0) / (nlats / 2.0);	// north up to 90d
        lats_s[j] = (j - nlats / 2.0) / (nlats / 2.0);	// south up to 90d
    }
    if (sinlat) {
        for (j = 0; j < nlats; j++) {
            lats_n[j] = asin(lats_n[j]); lats_s[j] = asin(lats_s[j]);
        }
    } else {
        for (j = 0; j < nlats; j++) {
            lats_n[j] = lats_n[j] * M_PI; lats_s[j] = lats_s[j] * M_PI;
        }
    }

    /* Map and images */
    map = (double *)malloc(mapsz * sizeof(double));	// This is the working map
    for (i = 0; i < mapsz; i++) map[i] = synop[i];
    image_n = (double *)malloc(imgsz * sizeof(double));
    image_s = (double *)malloc(imgsz * sizeof(double));

    /* Remapping */
    if (verbflag) SHOW("remapping... ");
    polar_proj_fwd(0, sinlat, image_n, imgDims[0], imgDims[1], 
            map, nlons, nlats, lons, lats_n, (option.lat0) * DTOR);
    polar_proj_fwd(1, sinlat, image_s, imgDims[0], imgDims[1], 
            map, nlons, nlats, lons, lats_n, (option.lat0) * DTOR);
    if (verbflag) SHOW("done.\n");

    /* Fitting */
    if (verbflag) SHOW("fitting... ");
    fitimage(image_n, imgDims[0], imgDims[1], order, 1);
    fitimage(image_s, imgDims[0], imgDims[1], order, 1);
    if (verbflag) SHOW("done.\n")

    /* Mapping back */
    if (verbflag) SHOW("mapping back... ");
    polar_proj_bwd(0, sinlat, image_n, imgDims[0], imgDims[1], 
            map, nlons, nlats, lons, lats_n, (option.lat0) * DTOR);
    polar_proj_bwd(1, sinlat, image_s, imgDims[0], imgDims[1], 
            map, nlons, nlats, lons, lats_n, (option.lat0) * DTOR);
    if (verbflag) SHOW("done.\n");

    /* Filling */
    if (sinlat) {
        dj = floor((1. - sin(fabs(option.latfil))) / (2. / nlats));
    } else {
        dj = floor((90. - fabs(option.latfil)) / (360. / nlats));
    }
    if (option.nflag) {
        for (j = nlats - dj; j < nlats; j++)
            for (i = 0; i < nlons; i++)
                synop[j * nlons + i] = map[j * nlons + i];
    }
    if (option.sflag) {
        for (j = 0; j < dj; j++)
            for (i = 0; i < nlons; i++)
                synop[j * nlons + i] = map[j * nlons + i];
    }

    /* Clean up */
    free(lats_n); free(lats_s);
    free(map); free(image_n); free(image_s);
}



/* ################ Interpolations ################ */


/* Get interpolation index of x0 in vector x */
double get_index(sinlat, n, x, x0)
    int sinlat, n;
    double x[], x0;
{
    int i, i0, i1;
    double ind;

    if (x0 <= x[0]) {
        i0 = 0; i1 = 1;
    } else if (x0 > x[n - 1]) {
        i0 = n - 2; i1 = n - 1;
    } else {
        for (int i = 0; i < n - 1; i++) {
            if (x0 > x[i] && x0 <= x[i + 1]) {
                i0 = i; i1 = i + 1;
                break;          
            }
        }
    }
    /* Sine Latitude? */
    if (sinlat) { 
        ind = i0 + (sin(x0) - sin(x[i0])) / (sin(x[i1]) - sin(x[i0]));
    } else { 
        ind = i0 + (x0 - x[i0]) / (x[i1] - x[i0]);
    }
    return ind;
}


/* Bilinear interpolation (from obs2helio).
 * If any neighboring point is NaN the result is NaN.
 */ 
double linintd(f, nx, ny, x, y)
    double *f;
    int nx, ny;
    double x, y; 
{   
    double p0, p1, q0, q1;      // weights
    int ilow, jlow;             // selected pixels
    double *fptr;               // input array temp
    double u;

    ilow = (int)floor(x);
    jlow = (int)floor(y);
    if (ilow < 0) ilow = 0;
    if (ilow + 1 >= nx) ilow -= 1;
    if (jlow < 0) jlow = 0;
    if (jlow + 1 >= ny) jlow -= 1;
    p1 = x - ilow;
    p0 = 1.0 - p1;
    q1 = y - jlow;
    q0 = 1.0 - q1;

    /* Bilinear interpolation from four data points. */
    u = 0.0;
    fptr = f + jlow * nx + ilow;
    u += p0 * q0 * (*fptr++);
    u += p1 * q0 * (*fptr);
    fptr += nx - 1;
    u += p0 * q1 * (*fptr++);
    u += p1 * q1 * (*fptr);
    return u;
}



/* ################ Map Projection Function ################ */

/* Prototype:
 * proj_func(south, lon, lat, x, y, xc, yc, rx, ry, nx, ny)
 * south: 0 or 1
 * lon, lat: in radius
 * x, y: image coord, left-bottom as (0, 0)
 * xc, yc: center pixel coord
 * rx, ry: radius of equator, in pixel
 * nx, ny: dimension of image
 */


/* Forward stereographic projection */
void stereo_fwd(south, lon, lat, 
            x, y, xc, yc, rx, ry, nx, ny)
    int south;
    double lon, lat;
    double *x, *y, xc, yc, rx, ry;
    int nx, ny;
{
    double x0, y0;	// Normalized
    double r, th;	// Polar coord

    r = cos(lat) / (1 + sin(fabs(lat)));
    th = lon;
    x0 = r * cos(th);
    y0 = r * sin(th);
    *x = x0 * rx + xc;
    *y = y0 * ry + yc;
}


/* Backward stereographic projection */
void stereo_bwd(south, lon, lat, 
            x, y, xc, yc, rx, ry, nx, ny)
    int south;
    double *lon, *lat;
    double x, y, xc, yc, rx, ry;
    int nx, ny;
{
    double x0, y0;	// Nomalized
    double r, th;	// Polar coord

    x0 = (x - xc) / rx;
    y0 = (y - yc) / ry;
    r = sqrt(x0 * x0 + y0 * y0);
    th = atan2(y0, x0);
    *lon = th;
    *lat = M_PI / 2.0 - atan(r) * 2.0;
    if (south) *lat = -(*lat);
}



/* ################ Polar Projections ################ */


/* This is an ad hoc function that projects the poles of a map to a specified 
 * image, using certain map projection. The result is a polar view from one of 
 * the poles. Reverse projection formula is used to find the map coordinates
 * of a certain image point. A bilinear interpolation is then used to get the 
 * value. in the image. During projection, the pole is projected to the center
 * of the image, and the lower latitude limit is projected to an elipse that is 
 * tangent to the edge of the image. Note there's no error checking here. Any 
 * neighboring missing data during interpolation results in NaN. This code handles
 * double precision data. Both sinlat and lat format can be used.
 */
void polar_proj_fwd(south, sinlat, image, nx, ny,
           map, nlons, nlats, lons, lats, lat_lim)
    int south, sinlat;
    double *image;
    int nx, ny;
    double *map;
    int nlons, nlats;
    double lons[], lats[], lat_lim;
{
    int i, j;
    double lon, lat, ind_lon, ind_lat;
    double x, y, xc, yc, r_lim, rx, ry;
    void (*proj_bwd)();
    double (*interpolate)();

    switch (projname) {
        case ((enum proj)(0)):
            proj_bwd = stereo_bwd;
            break;
        default:
            proj_bwd = stereo_bwd;
            break;
    }
    switch (intpname) {
        case ((enum intp)(0)):
            interpolate = linintd;
            break;
        default:
            interpolate = linintd;
            break;
    }

    // Find center
    xc = (nx - 1) / 2.0;
    yc = (ny - 1) / 2.0;
    // Normalized radius
    r_lim = cos(lat_lim) / (1 + sin(fabs(lat_lim)));
    rx = xc / r_lim;
    ry = yc / r_lim;

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            /* Projection */
            (*proj_bwd)(south, &lon, &lat,
                    (double)i, (double)j, xc, yc, rx, ry, nx, ny);
            if (fabs(lat) < fabs(lat_lim)) {
                image[j * nx + i] = NAN;
                continue;
            }
            /* Get index */
            ind_lon = get_index(0, nlons, lons, lon);
            ind_lat = get_index(sinlat, nlats, lats, lat);
            /* Do interpolation */
            image[j * nx + i] = 
                (*interpolate)(map, nlons, nlats, ind_lon, ind_lat);
        }
    }
}


/* This is an ad hoc function that projects a polar projection image back to
 * a map. Forward projection formula is used to find corresponding image point
 * for each map point and an interpolation is done on the image. Handle doubles.
 */
void polar_proj_bwd(south, sinlat, image, nx, ny,
           map, nlons, nlats, lons, lats, lat_lim)
    int south, sinlat;
    double *image;
    int nx, ny;
    double *map;
    int nlons, nlats;
    double lons[], lats[], lat_lim;
{
    int i, j, j0, j1;
    double x, y, xc, yc, r_lim, rx, ry;
    void (*proj_fwd)();
    double (*interpolate)();

    switch (projname) {
        case ((enum proj)(0)):
            proj_fwd = stereo_fwd;
            break;
        default:
            proj_fwd = stereo_fwd;
            break;
    }
    switch (intpname) {
        case ((enum intp)(0)):
            interpolate = linintd;
            break;
        default:
            interpolate = linintd;
            break;
    }

    // Find center
    xc = (nx - 1) / 2.0;
    yc = (ny - 1) / 2.0;
    // Normalized radius
    r_lim = cos(lat_lim) / (1 + sin(fabs(lat_lim)));
    rx = xc / r_lim;
    ry = yc / r_lim;

    if (south) {
        j0 = 0;
        j1 = ceil(get_index(sinlat, nlats, lats, -fabs(lat_lim)));
    } else {
        j0 = floor(get_index(sinlat, nlats, lats, fabs(lat_lim)));
        j1 = nlats - 1;
    }

    for (j = j0; j <= j1; j++) {
        for (i = 0; i < nlons; i++) {
            /* Projection */
            (*proj_fwd)(south, lons[i], lats[j],
                    &x, &y, xc, yc, rx, ry, nx, ny);
            /* Do interpolation */
            if (x < 0 || x >= nx || y < 0 || y >= ny) {
                continue;	// Do nothing
            } else {
                map[j * nlons + i] = 
                    (*interpolate)(image, nx, ny, x, y);
            }
        }
    }
}

