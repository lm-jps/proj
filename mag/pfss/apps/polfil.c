/*
 * Module name:		polfil.c
 *
 * Description:		This module takes a synoptic map and fill in polar regions, using
 *			certain interpolation algorithms. By default, the polar region is
 *			remapped and a least square fitting is used to determine the missing
 *			points. See notes below for detailed options.
 *
 * Original source:	Image fitting code by Rasmus Larsen (rmlarsen@gmail.com)
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Oct 01 2009
 *			v1.1		Oct 05 2009
 *			v1.2		Mar 13 2010
 *			v1.3		Mar 16 2010
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
 *
 *
 * Example:		polfil "in=su_xudong.synop720[1920]" out=su_xudong.testpolfil
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

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif


/* ################## Timing from T. Larson ################## */

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
  *utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec / 1000.0;
  *stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec / 1000.0;
  return *utime + *stime;
}



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

/* Image fitting */
void dgels_(const char *, int *, int *, int *, double *, int *, double *, 
	    int *, double *, int *, int *);	// Fortran LAPACK
void cheby_basis(int n, int order, double x[], double cheb[], int ldc);

void fitimage(double image[], int m, int n, int order, int fill);




/* ################## Main Module ################## */


char *module_name = "polfil";	/* Module name */

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
    {ARG_INT, "ORDER", "3", "Highest power of polynomials."},
    {ARG_END}
};

int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS, *outRS;
    int irec, nrecs;
    DRMS_Record_t *inRec, *outRec;
    DRMS_Segment_t *inSeg, *outSeg;
    DRMS_Array_t *inArray, *outArray;
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
    int car_rot;
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

    /* Get parameters */
    inQuery = params_get_str(&cmdparams, "in");
    outQuery = params_get_str(&cmdparams, "out");
    method = params_get_str(&cmdparams, "METHOD");
    for (itest = 0; itest < ARRLENGTH(methodlist) && !methodflag; itest++) {
        if (!strcmp(method, methodlist[itest])) {
            methodflag = 1;
            methodnum = (enum methodnums)(itest);
        }
    }
    if (!methodflag) DIE("Unknow method");
    interp = params_get_str(&cmdparams, "INTERP");
    for (itest = 0; itest < ARRLENGTH(intplist) && !intpflag; itest++) {
        if (!strcmp(interp, intplist[itest])) {
            intpflag = 1;
            intpname = (enum intp)(itest);
        }
    }
    if (!intpflag) DIE("Unknow interpolation method");
    project = params_get_str(&cmdparams, "PROJECT");
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
    ctype2 = drms_getkey_string(inRec_tmp, "CTYPE2", &status);
    sinlat = strcmp(ctype2, "Sine Latitude") ? 0 : 1;

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

        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }

        /* Parameters */
        car_rot = drms_getkey_int(inRec, "CAR_ROT", &status);
        n0 = inArray->axis[0]; n1 = inArray->axis[1];
        mapsz = n0 * n1;

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

        /* Result writing */
        status = drms_segment_write(outSeg, outArray, 0);
        if (status) DIE("problem writing data");
        drms_free_array(outArray);
        free(lons); free(lats);

        /* Set keywords */
        drms_setkey_int(outRec, "CAR_ROT", car_rot);
        drms_setkey_string(outRec, "METHOD", method);
        drms_setkey_double(outRec, "LAT0", option.lat0);
        drms_setkey_double(outRec, "LATFIL", option.latfil);
        drms_setkey_int(outRec, "NFLAG", option.nflag);
        drms_setkey_int(outRec, "SFLAG", option.sflag);

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
}



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



/* ################ Chebyshev Polynomials ################ */


/* Generate matrix with Chebychev polynomials T_0,T_1,...,T_{order-1}
 * evaluated on the gridpoints in x.
 */
void cheby_basis(n, order, x, cheb, ldc)
    int n, order;
    double x[], cheb[];
    int ldc;
{
    int i, s;

    for (i = 0; i < n; i++) {
        cheb[i] = 1;
        cheb[i + ldc] = x[i];
    }  
    for (s = 2; s < order; s++) {
        for (i = 0; i < n; i++) {
            cheb[s * ldc + i] = 
                2 * x[i] * cheb[(s - 1) * ldc + i] - cheb[(s - 2) * ldc + i];
        }
    }
}  



/* ################ Image Fitting ################ */


void fitimage(image, m, n, order, fill)
    double image[];
    int m, n, order, fill;
{
    int i, j, k, l, s, t, npix, lda, mn;
    int nterms, info, lwork, *map[2], nimgs;
    double *A, *B, *x, *y, *work, *chebx, *cheby;
    double pix;

    if (m <= 0 || n <= 0 || order <= 1)
        return;
    nimgs = 1;

    /* Allocate memory for matrix, right-hand sides etc. */
    lda = m * n;
    A = (double *)malloc(m * n * order * order * sizeof(double));
    B = (double *)malloc(m * n * sizeof(double));
    map[0] = (int *)malloc(n * m * sizeof(int));
    map[1] = (int *)malloc(n * m * sizeof(int));
    x = (double *)malloc(n * sizeof(double));
    y = (double *)malloc(m * sizeof(double));
    chebx = (double *)malloc(order * n * sizeof(double));
    cheby = (double *)malloc(order * m * sizeof(double));
    nterms = order * order;

    /* Set up abscissa grid on [-1;1]x[-1;1] and Chebychev polynomials. */
    for (i = 0; i < n; i++) 
        x[i] = (double)(2.0 * i) / (n - 1) - 1;
    cheby_basis(n, order, x, chebx, n);
    for (j = 0; j < m; j++)
        y[j] = (double)(2.0 * j) / (m - 1) - 1;
    cheby_basis(m, order, y, cheby, m);

    /* Extract "good" pixels and set up matrix with bivariate polynomials
       in the corresponding points. */
    npix = 0;
    for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
        pix = image[i * m + j];
        if (!isnan(pix)) {
            B[npix] = pix;
            map[0][npix] = i;
            map[1][npix] = j;
            npix++;
        }
    }
    for (s = 0; s < order; s++)
    for (t = 0; t < order; t++)
    for (l = 0; l < npix; l++) {
        i = map[0][l];
        j = map[1][l];
        A[(m * n) * (s * order + t) + l] = chebx[n * s + i] * cheby[m * t + j];
    }

    /* Extract remaining right-hand sides */
    for (i = 0; i < npix; i++) {
        j = map[0][i] * m + map[1][i];
        B[i] = image[j];
    }

    /* Solve the least-squares problems. */
    mn = MIN(npix, nterms);
    lwork = MAX(1, mn + MAX(npix, MAX(nterms, 1)) * 32);
    work = (double *)malloc(lwork * sizeof(double));

    dgels_("n", &npix, &nterms, &nimgs, A, &lda, B, &lda, work, &lwork, &info);

    if (info < 0) {
        fprintf(stderr, "DGELS: The %dth argument had an illegal value.\n", -info);
        exit(1);
    }

    /* Overwrite image with the smoothed images corresponding to the
       least-squares solutions. */
    if (fill) {
        for (i = 0; i < n; i++)
	for (j = 0; j < m; j++) {
	    image[i * m + j] = 0;
        }
        for (s = 0; s < order; s++)
	for (t = 0; t < order; t++)
	for (i = 0; i < n; i++)
	for (j = 0; j < m; j++) {
	    image[i * m + j] +=
                B[s * order + t] * chebx[n * s + i] * cheby[m * t + j];
        }
    } else {
        for (l = 0; l < npix; l++) {
	    i = map[0][l];
	    j = map[1][l];
	    image[i * m + j] = 0;
        }
        for (s = 0; s < order; s++) 
	for (t = 0; t < order; t++)
	for (l = 0; l < npix; l++) {
	    i = map[0][l];
	    j = map[1][l];
	    image[i * m + j] +=
                B[s * order + t] * chebx[n * s + i] * cheby[m * t + j];
	}
    }

    /* Free memory */
    free(map[0]); free(map[1]);
    free(x); free(y);
    free(chebx); free(cheby);
    free(A); free(B); free(work);
}
