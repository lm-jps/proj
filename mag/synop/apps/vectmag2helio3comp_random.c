/*  this JSOC module is a combination of 3 SOI modules: v2helio, helio2mlat, qdotprod
 *  ported by Tim Larson
 *  removes check on mapbmin, mapbmax, and maprows, these keywords are no longer carried
 */

/*
 *  v2helio.c                            ~soi/(version)/src/modules/v2helio.c
 *
 *  This module interpolates CCD velocity data to estimate values that
 *    would be obtained at specified equal increments of heliographic 
 *    longitude and sine of latitude.  Apodization and corrections for
 *    solar rotation and limbshift are included. 
 *
 *  Responsible:  Kay Leibrand                   KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *    This module is under development.  Look for ??? in comments.
 *
 *  Planned updates:
 *    Decide when and where and if default values should be used.
 *
 */

/*
 *  helio2mlat.c                     ~soi/(version)/src/modules/helio2mlat.c
 *
 *  Description:
 *     Adapted by Kay Leibrand from pipeLNU shtfft
 *     fft by map_rows.  Use FORTRAN library Cmlib for float version.
 *     extends data to 360 degrees prior to transform but only saves
 *     cols up to lmax.  Performs transpose needed prior to dotprod.
 *
 *  Responsible:  Kay Leibrand                   KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *    This module is under development.  Look for ??? in comments.
 *
 *  Planned updates:
 *    Restructure to use functions?
 *    Refine parameter definitions and names to be consistent with
 *       new keywords and ancillary data flow.
 *    Fix phase
 *
 */

/*
 *  qdotprod.c                     ~soi/(version)/src/modules/qdotprod.c
 *
 *  Description:
 *    Conversion of Jesper Schou's FORTRAN q(uick)dotprod to C 
 *    from file ~schou/testdot/testdot4c.f 
 *    uses FORTRAN functions from blas library
 *    contains optimizations, calculates masks, allows "chunking" in l 
 *
 *  Responsible:  Kay Leibrand                   KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *    Inadequate checks for consistency between data specifications, 
 *      i.e. dataset names, and LMIN, LMAX, LCHUNK parameters.
 *    This module is under development.  Look for ??? and TBD in comments.
 *    Normalization in plm's is different from PHS pipeLNU
 *
 *  Planned updates:
 *    Fix known bugs. 
 *    Use a consistent pointer style, i.e. pointer arith vs [] 
 *    Restructure to use functions.
 *    Refine parameter definitions and names to be consistent with
 *       new keywords and ancillary data flow.
 *
 *  Revision history is at end of file.
 */

/* Normalization of the resulting time-series when norm != 0:

  Let the observed surface behaviour of an oscillation be given by 
  Re(exp(-i\omega t) Y_l^m (\theta,\phi)) with (Y_l^m)^2 normalized
  to have an average of 1 over the unit sphere (this is 4\pi times
  the usual definition where the integral is 1). Assume that the
  whole (4\pi) Sun is observed with no velocity projection factor.
  Then the resulting time series is given by exp(-i\omega t).

  This is equivalent to preserving the rms value of the real parts
  of the oscillations.

  And maybe I got the implementation straight.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fftw3.h>
#include "jsoc_main.h"
#include "fstats.h"
#include "drms_dsdsapi.h"
#include "errlog.h"
#include "projection.h"
#include "atoinc.h"
//#include "magutils.h"
//#include "cartography.h"
//#include "img2helioVector.c"

// include a couple of subroutines
//#include "obs2helio.c"
//#include "setplm2.c"
//#include "apodize.c"
// end

#define ARRLENGTH(ARR)  (sizeof(ARR)/sizeof(ARR[0]))
#define PI              (M_PI)

#define absval(x)		(((x) < 0) ? (-(x)) : (x))
#define minval(x,y)		(((x) < (y)) ? (x) : (y))
#define maxval(x,y)		(((x) < (y)) ? (y) : (x))
#define very_small		(1.0e-30)
#define is_very_small(x)	(absval(x) < very_small)
#define is_even(x) 		(!((x)%2))
#define is_odd(x) 		((x)%2)

#define RADSINDEG 	(PI/180)
#define ARCSECSINRAD 	(3600*180/PI)
#define DAYSINYEAR	(365.2425)
#define SECSINDAY	(86400)
#define TAU_A           (499.004783806) // light travel time in seconds, = 1 AU/c
//#define TAU_A		(499.004782)	// this value used in old v2helio
#define TCARR		(25.38)		// days
#define RTRUE		(6.96000000e8)	// meters
#define AU		(149597870691)	// meters/au
//#define AU		(1.49597870e11)	// this value used in old v2helio

#define QUAL_NODATA	 (0x80000000)
#define QUAL_MIXEDCALVER (0x00000001)
#define kNOTSPECIFIED  "not specified"

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}

// define remap-img2helioVector
extern int remap_img2helioVector (double bxImg, double byImg, double bzImg, double *bxHelio,
                     double *byHelio, double *bzHelio, double lon, double lat,
                     double lonc, double latc, double pAng); 

// Defined in synop-cartography.c:
extern int synop_plane2sphere (double x, double y, double latc, double lonc, double *lat, double *lon, int projection);
extern int synop_img2sphere(double x, double y, double ang_r, double latc, double lonc, double pa, double *rho, double *lat, double *lon, double *sinlat, double *coslat, double *sig, double *mu, double *chi);
extern int synop_sphere2img(double lat, double lon, double latc, double lonc, double *x, double *y, double xcenter, double ycenter, double rsun, double peff, double ecc, double chi, int xinvrt, int yinvrt);
extern int synop_sphere2plane(double lat, double lon, double latc, double lonc, double *x, double *y, int projection);

char *module_name = "vectmag2helio3comp_random";
//char *cvsinfo_jv2ts = "cvsinfo: $Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/vectmag2helio3comp_random.c,v 1.2 2016/05/17 16:31:31 yliu Exp $";

ModuleArgs_t module_args[] = 
{
// these inputs are common to all modules
   {ARG_STRING, "in",           NULL,           "input data records"},
   {ARG_STRING, "tsout",        kNOTSPECIFIED,  "output data series"},
   {ARG_STRING, "segin",        kNOTSPECIFIED,  "input data segment"},
   {ARG_STRING, "segout",       kNOTSPECIFIED,  "output data segment"},
   {ARG_STRING, "histlink",     "HISTORY",      "name of link to ancillary dataseries for processing metadata. specify \"none\" to suppress warning messages."},
   {ARG_STRING, "srclink",      "SRCDATA",      "name of link to source data"}, // used only for jv2helio and jhelio2mlat output
   {ARG_STRING, "v2hout",       kNOTSPECIFIED,  "output data series for jv2helio"},
   {ARG_STRING, "h2mout",       kNOTSPECIFIED,  "output data series for jhelio2mlat"},
   {ARG_INT,    "PERM",         "1",            "set to 0 for transient records, nonzero for permanent records"},
   {ARG_INT,    "VERB",         "1",            "option for level of verbosity: 0 gives only error and warning messages; >0 prints messages outside of loop; >1 prints messages inside loop; >2 for debugging output"},
   {ARG_INT,    "FORCEOUTPUT",  "0",            "option to specify behavior on missing inputs; 0 gives an error on missing or duplicate inputs, >0 makes outputs regardless"},
   {ARG_STRING, "TAG",          "none",         "this parameter sets a keyword of the same name to the same value, usually used as an additional primekey"},
   {ARG_STRING, "VERSION",      kNOTSPECIFIED,  "this parameter sets a keyword of the same name to the same value, useful for selecting obsolete records"},
   {ARG_INT,    "CALVERKEY",    "2",            "short integer bit mask determining which 4-bit fields of CALVER64 are permitted to change on input.  set the bit to disallow change of the corresponding nybble."},
// these are from jqdotprod
   {ARG_INT,    "LMIN",         "0",            "minimum l in output"},
//   {ARG_INT,    "LMAX",         "0",            "maximum l in output, take from input if 0"},   /* if 0, default is LMAX of in_sds */
   {ARG_INT,    "LCHUNK",       "0",            "increment in l on output, default is lmax+1"},   /* if 0, is LMAX+1 */ 
   {ARG_INT,    "NORM",         "1",            "set to nonzero to use cnorm=sinbdelta*sqrt(2) in sgemm call, otherwise use cnorm=1"},   /* Uses old norm if =0, see below */ 
   {ARG_TIME,   "TSTART",       NULL,           "start time of first output record"},
   {ARG_STRING, "TTOTAL",       NULL,           "total length of time in output"},
   {ARG_STRING, "TCHUNK",       kNOTSPECIFIED,  "length of output timeseries"},
// these are from jhelio2mlat
   {ARG_INT,    "LMAX",         "300",       "maximum l (maximum m) in the output, cannot be greater than MAPMMAX", NULL},
   {ARG_INT,    "SUBMEAN",      "0",         "nonzero subtracts mean of input image", NULL},
   {ARG_INT,    "NORMLIZE",     "0",         "nonzero multiplies by sqrt((fraction nonmissing)/2) for each row", NULL},
   {ARG_INT,    "CENTLONG",     "1",         "nonzero centers the longitude Fourier transform on the center of the remapped image", NULL},
   {ARG_INT,    "ZEROMISS",     "0",         "zero sets any row containing a NaN to DRMS_MISSING", NULL},
   {ARG_INT,    "LGAPOD",       "0",         "nonzero apodizes in longitude", NULL},
   {ARG_DOUBLE, "LGAPMIN",      "60.0",      "start of longitude apodization, degrees", NULL},
   {ARG_DOUBLE, "LGAPWIDT",     "10.0",      "width of longitude apodization, degrees", NULL},
// these are from jv2helio
   {ARG_INT,     "MAPMMAX",     "5402",      "determines mapcols", ""},	/* determines mapcols, default value is 3*512 */
   {ARG_INT,     "SINBDIVS",    "2160",       "number of increments in sin latitude from 0 to 1", ""},	/* # of = increments in sinB from sin(0) to sin(PI/2) */
   {ARG_FLOAT,   "MAPRMAX",     "0.998",      "maximum image radius", ""},
   {ARG_INT,     "NAN_BEYOND_RMAX",    "0",  "set to nonzero to set output to DRMS_MISSING outside MAPRMAX, otherwise uses 0.0 outside MAPRMAX"},  /* Non 0 sets data outside RMAX MISSING */   
   {ARG_FLOAT,   "MAPLGMAX",    "90.0",      "longitude maximum, degrees", ""},	/* degrees */    
   {ARG_FLOAT,   "MAPLGMIN",    "-90.0",     "longitude minimum, degrees", ""},
   {ARG_FLOAT,   "MAPBMAX",     "90.0",      "latitude maximum, degrees, also used for minimum", ""},  
   {ARG_INT,     "LGSHIFT",     "3",         "option for longitude shift: 0=none; 1=fixed rate; 2=nearest degree; 3=nearest tenth of a degree", ""}, /* 0=none; 1=fixed rate; 2=nearest Degree */
   {ARG_TIME,    "REF_T0",      "1987.01.03_17:31:12_TAI", "reference time for computing fixed rate longitude shift", ""},
   {ARG_FLOAT,   "REF_L0",      "0.0",       "reference longitude for computing fixed rate longitude shift ", ""},
   {ARG_FLOAT,   "SOLAR_P",     "999.0",     "P-angle; if unset, taken from keywords", ""},	/* can't use D_MISSING here */
   {ARG_FLOAT,   "PSIGN",       "1.0",       "sign of SOLAR_P", ""},	/* Sign of P. For MWO data. */
   {ARG_FLOAT,   "PERR",        "0.0",       "fixed P-angle error, likely -0.22", ""},	/* Fixed P-angle error. Maybe -0.22. */
   {ARG_FLOAT,   "IERR",        "0.0",       "error in Carrington inclination, likely -0.10", ""},	/* Error in Carrington inclination. Maybe -0.10. */
   {ARG_TIME,    "REF_TB0",     "2001.06.06_06:57:22_TAI", "reference time for computing correction to P and B angles, roughly when B0=0", ""},
   {ARG_INT,     "INTERPO",     "1",         "option for interpolation: 0=bilinear; 1=cubic convolution", ""},	/* 2 methods - see soi_fun.h */
   {ARG_INT,     "APODIZE",     "0",         "option for apodization: 0=none; 1=use true solar coordinates; 2=use ideal solar coordinates (b0=0)", ""},	/* see soi_fun.h or apodize.c */
   {ARG_FLOAT,   "APINNER",     "0.90",      "start of apodization in fractional image radius", ""},	/* start of apodization */
   {ARG_FLOAT,   "APWIDTH",     "0.05",      "width of apodization in fractional image radius", ""},	/* width of apodization */
   {ARG_INT,     "APEL",        "0",         "set to nonzero for elliptical apodization described by APX and APY", ""},	/* do elliptical apodization described by apx and apy */
   {ARG_FLOAT,   "APX",         "1.00",      "divide the x position by this before applying apodization", ""},	/* divide the x position by this before applying apodization */
   {ARG_FLOAT,   "APY",         "1.00",      "divide the y position by this before applying apodization", ""},	/* divide the y position by this before applying apodization */
   {ARG_INT,     "VCORLEV",     "0",         "option for velocity correction: 0=none; 1=subtract a model of differential rotation; 2=also divide by line of sight projection factor for purely radial velocities", ""}, 	/* 3 levels - see soi_fun.h*/
   {ARG_INT,     "MCORLEV",     "0",         "option for magnetic correction: 0=none; 1=line of sight; 2=radial", ""}, 	/* 2 levels - see soi_fun.h*/
   {ARG_INT,     "MOFFSETFLAG", "0",         "set to nonzero to get BFITZERO from input record and subtract from data before interpolating", ""}, 	/* 1=apply BFITZERO correction*/
   {ARG_FLOAT,   "OUTSCALE",    "1.0",       "bscale to use for output", ""},    /* scale for output */
   {ARG_FLOAT,   "OUTBIAS",     "0.0",       "bzero to use for output", ""},    /* bias for scaled output */
   {ARG_INT,     "DISTORT",     "0",         "option for distortion correction: 0=none; 1=full disk(fd) data; 2=vector-weighted(vw) data", ""}, /* 0 for none, 1 for FD, 2 for vw */
   {ARG_FLOAT,   "CUBIC",       "7.06E-9",   "cubic distortion in fd units", ""}, /* Cubic distortion in FD units */
   {ARG_FLOAT,   "TILTALPHA",   "2.59",      "tilt of CCD, degrees", ""}, /* TILT of CCD in degrees */
   {ARG_FLOAT,   "TILTBETA",    "56.0",      "direction of CCD tilt, degrees", ""}, /* Direction of TILT in degrees */
   {ARG_FLOAT,   "TILTFEFF",    "12972.629", "effective focal length", ""}, /* Effective focal length */
   {ARG_INT,     "OFLAG",       "0",         "set to nonzero to force reading orientation from keyword, otherwise \"SESW\" is assumed)", ""},  /* Non 0 skips checko (SESW assumed) */
   {ARG_INT,     "DATASIGN",    "0",         "value to multiply data; set to 0 to take DATASIGN from keyword, or 1.0 if not found", ""}, 	/* Non 0 forces datasign to value*/
   {ARG_INT,     "MAXMISSVALS", "0",         "if >0, this becomes threshold on MISSVALS from keyword", ""},  /* max. allowed MISSING pixels */
   {ARG_INT,     "CARRSTRETCH", "1",         "set to nonzero to correct for differential rotation according to DIFROT_[ABC]", ""},  /* 0 - don't correct for diff rot, 1 - correct */
   {ARG_FLOAT,   "DIFROT_A",    "13.562",    "A coefficient in differential rotation adjustment (offset)", ""}, /* A coefficient in diff rot adj (offset) */
   {ARG_FLOAT,   "DIFROT_B",    "-2.04",     "B coefficient (to sin(lat) ^ 2)", ""},  /* B coefficient (to sin(lat) ^ 2) */
   {ARG_FLOAT,   "DIFROT_C",    "-1.4875",   "C coefficient (to sin(lat) ^ 4)", ""},  /* C coefficient (to sin(lat) ^ 4) */
   {ARG_FLOAT,   "RESCALE",     "0.333333",  "Scale factor."}, // YLiu

   {ARG_END}
};

#include "saveparm.c"
#include "timing.c"
#include "set_history.c"
#include "calversfunctions.c"

// define remap-obs2helio
extern int remap_SetDistort(int dist, double cubic, double alpha, double beta, double feff, LIBPROJECTION_Dist_t *dOut);

extern int remap_obs2helio(float *V, float *U, int xpixels, int ypixels, double x0, double y0, double BZero, double P,
              double   S, double rsun, double Rmax, int	interpolation, int cols, int rows, double Lmin,
              double Ldelta, double Ladjust, double sinBdelta, double smajor, double sminor, double sangle,
              double xscale, double yscale, const char *orientation, int mag_correction, int velocity_correction,
              double obs_vr, double obs_vw, double obs_vn, double vsign, int NaN_beyond_rmax, int carrStretch,
              const LIBPROJECTION_Dist_t *distpars, float diffrotA, float diffrotB, float diffrotC, 
              LIBPROJECTION_RotRate_t *rRates, int size);
// define remap-apodize
extern int remap_apodize(float *data, double b0, int cols, int rows, double Lmin, double Ldelta, double sinBdelta, 
            int apodlevel, double apinner, double apwidth, int apel, double apx, double apy);

// define remap-setplm2
extern void remap_setplm2(int lmin,int lmax,int m,long nx,int *indx,double *x,long nplm,double *plm,double *dplm);

//char *getshtversion(void);

/* forward decls for fortran functions */
void scopy_(int *, float *, int *, float *, int *);
void setplm_(int *, int *, int *, int *, int *, double *, int *, double *);
void sgemm_(); /* give up */
void do_boxcar(float *image_in, float *image_out, int in_nx, int in_ny, float fscale);

typedef enum
{
   V2HStatus_Success,
   V2HStatus_MissingParameter,
   V2HStatus_IllegalTimeFormat,
   V2HStatus_TimeConvFailed,
   V2HStatus_Unimplemented,
   V2HStatus_IllegalOrientation
} V2HStatus_t;

static void CheckO(const char *orientation, V2HStatus_t *stat)
{
   /* check for legal orientation string */
   static char o[16];
   char *c = o;

   if(!orientation)
   {
      *stat = V2HStatus_MissingParameter;
   } 
   else if (4 != sscanf(orientation, "%[NS]%[EW]%[NS]%[EW]", c++, c++, c++, c++))
   {
      *stat = V2HStatus_IllegalOrientation;
   }
   else if ((o[0] == o[2]) && (o[1] == o[3]))
   { 
      *stat = V2HStatus_IllegalOrientation;
   }
   else if ((o[0] !=o [2]) && (o[1] != o[3])) 
   {
      *stat = V2HStatus_IllegalOrientation;
   }
   else 
   {
      *stat = V2HStatus_Success;   
   }
}

static inline void ParameterDef(int status, 
				char *pname, 
				double defaultp, 
				double *p, 
				long long recnum,
                                int verbflag)
{
   /* logs warning and sets parameter to default value */
   if (status != 0)
   {
      *p = defaultp;
      if (verbflag)
         fprintf(stderr, "WARNING: default value %g used for %s, status = %d, recnum = %lld\n", defaultp, pname, status, recnum);
   }
}

#define PARAMETER_ERROR(PNAME) \
if (status != DRMS_SUCCESS) \
{ \
  fprintf(stderr, "ERROR: problem with required keyword %s, status = %d, T_REC = %s, recnum = %lld, histrecnum = %lld\n", PNAME, status, trecstr, inrec->recnum, histrecnum); \
  drms_free_array(inarrBtotal); \
  drms_free_array(inarrBincl); \
  drms_free_array(inarrBazim); \
  drms_free_array(inarrBdisamb); \
  return 0; \
}


int DoIt(void)
{ 
  int newstat=0;
  DRMS_RecChunking_t chunkstat = kRecChunking_None;
  int fetchstat = DRMS_SUCCESS;
  int status = DRMS_SUCCESS;
  char *inrecquery = NULL;
  char *outseries = NULL;
  char *segnamein = NULL;
  char *segnameout = NULL;
  DRMS_RecordSet_t *inrecset = NULL;
  DRMS_Record_t *inrec = NULL;
  DRMS_Record_t *outrec = NULL;
  DRMS_Segment_t *segin = NULL;
  DRMS_Segment_t *segout = NULL;
  DRMS_Array_t *inarrBdisamb = NULL, *inarrBtotal = NULL, *inarrBincl = NULL, *inarrBazim = NULL;
  DRMS_Array_t *outarr = NULL;
  DRMS_Type_t usetype = DRMS_TYPE_FLOAT;
  DRMS_RecLifetime_t lifetime;
  long long histrecnum=-1;
  int length[2];

  double wt0, wt1, wt2, wt3, wt;
  double ut0, ut1, ut2, ut3, ut;
  double st0, st1, st2, st3, st;
  double ct0, ct1, ct2, ct3, ct;

  TIME trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
  char trecstr[100], tstartstr[100];

  int quality;
  float *bTotal, *bIncl, *bAzim, *disamb;
  int xDim, yDim, iData, ix, jy, yOff;
//  double bxHel, byHel, bzHel;
//  float *bRadial, *bTheta, *bPhi;
  float *h2mptr;
  float *v2hbr, *v2hbt, *v2hbp;

// from jqdotprod, norm changed to normflag
  float *oddpart, *evenpart, *inptr, *workptr, *mptr, *outptr;
  float *folded, *masks, *real4evenl, *real4oddl, *imag4evenl, *imag4oddl, *outx; 

/* used for setting up plm's */
  double *plm, *saveplm, *latgrid;
  int *indx;

  double sinBdelta;
  double *plmptr;
  float *maskptr;

  int nsn, fournsn, snx, maxnsn;
  int lchunksize, lchunkfirst, lchunklast, lchunk, l;
  int lmin, lmax;
  int msize, mx, foldedsize;
  int maprows, nlat, imagesize, latx, poslatx, neglatx, moffset, nlatx;
  int i, m, modd, meven;
  int lfirst, llast, ifirst, ilast, lstart, ldone;
  int lfirsteven, llasteven, nevenl, lfirstodd, llastodd, noddl;
  int fromoffset, tooffset, imageoffset;

  int increment1 = 1; /* for scopy call */
  int increment2 = 2; /* for scopy call */

/* arguments for sgemm call */
  char transpose[] = "t";
  char normal[] = "n";
  float one = 1.0;
  float zero = 0.0;
  int normflag;
  float cnorm; /* Constant to get proper normalization. */

  double tstart, tepoch, tstep, tround, cadence, nseconds, chunksecs, cadence0;
  char *ttotal, *tchunk;
  int nrecs, irec, trecind, bc, iset, ntimechunks;
  int *bad;
  int ndt;

// from jhelio2mlat, i, lmax duplicated.
  double mean, norm=1.0, normx;
  int subtract_mean, normalize, cent_long, zero_miss, lgapod;
  double lgapmin, lgapwidt, lgapmax, lon;

  int row, col; 
  int lfft, mapped_lmax; 
  int mapcols2;
  int nfft, nmean, nok, nout;

  float *buf, *bp, *ip, *inp, *op, *outp, *weight, *wp;
  float *wbuf;
  fftwf_plan fftwp;

// from jv2helio, row, mapped_lmax, maprows, sinBdelta  duplicated
  V2HStatus_t vstat = V2HStatus_Success;
  const char *orientationdef = "SESW    ";
  char *orientation = NULL;
  int paramsign;
  int longitude_shift, velocity_correction, interpolation, apodization;
  int mag_correction;
  int mag_offset;
  int sinb_divisions, mapcols, nmax, nmin;
  int carrStretch = 0;
  float diffrotA = 0.0;
  float diffrotB = 0.0;
  float diffrotC = 0.0;
  double tobs, tmearth, tref, trefb0;
  double smajor, sminor, sangle;
  double xscale, yscale, imagescale;
  int xpixels, ypixels, pixels;
  double obs_vr, obs_vw, obs_vn;
  double b0, bmax, bmin;
  double longmax, longmin, longmax_adjusted, longmin_adjusted, longinterval;
  double p0, p, rmax;
  double ierr, perr, psign;
  double x0, y0;
  double obsdist, longshift, obsl0, refl0, mapl0, longrate, rtrue, rsun, S;
  double rsunDef, rsunobs;
  int obsCR;
  int apel;
  double apinner, apwidth, apx, apy;
  double scale, bias;
  double colsperdeg;

  double satrot, instrot;
  double dsignout, vsign;
  int distsave;
  double cubsave, tiltasave, tiltbsave, tiltfsave;
  LIBPROJECTION_Dist_t distP;

  int errbufstat=setvbuf(stderr, NULL, _IONBF, BUFSIZ);
  int outbufstat=setvbuf(stdout, NULL, _IONBF, BUFSIZ);

  wt0=getwalltime();
  ct0=getcputime(&ut0, &st0);

  inrecquery = (char *)cmdparams_save_str(&cmdparams, "in", &newstat);
  outseries = (char *)cmdparams_save_str(&cmdparams, "tsout", &newstat);
  int tsflag = strcmp(kNOTSPECIFIED, outseries);
  segnamein = (char *)cmdparams_save_str(&cmdparams, "segin", &newstat);
  segnameout = (char *)cmdparams_save_str(&cmdparams, "segout", &newstat);
  int seginflag = strcmp(kNOTSPECIFIED, segnamein);
  int segoutflag = strcmp(kNOTSPECIFIED, segnameout);
  int verbflag = cmdparams_save_int(&cmdparams, "VERB", &newstat);
  int permflag = cmdparams_save_int(&cmdparams, "PERM", &newstat);
  if (permflag)
    lifetime = DRMS_PERMANENT;
  else
    lifetime = DRMS_TRANSIENT;
  int forceoutput = cmdparams_save_int(&cmdparams, "FORCEOUTPUT", &newstat);
  char *tag = (char *)cmdparams_save_str(&cmdparams, "TAG", &newstat);
  char *version = (char *)cmdparams_save_str(&cmdparams, "VERSION", &newstat);
  int verflag = strcmp(kNOTSPECIFIED, version);
  unsigned short calverkey = (unsigned short)cmdparams_save_int(&cmdparams, "CALVERKEY", &newstat);

  char *histlinkname = (char *)cmdparams_save_str(&cmdparams, "histlink", &newstat);
  char *srclinkname = (char *)cmdparams_save_str(&cmdparams, "srclink", &newstat);

  char *v2hout = (char *)cmdparams_save_str(&cmdparams, "v2hout", &newstat);
  char *h2mout = (char *)cmdparams_save_str(&cmdparams, "h2mout", &newstat);
  int v2hflag = strcmp(kNOTSPECIFIED, v2hout);
  int h2mflag = strcmp(kNOTSPECIFIED, h2mout);
  int histflag = strncasecmp("none", histlinkname, 4);

  if (!v2hflag && !h2mflag && !tsflag)
  {
    fprintf(stderr, "ERROR: no outputs specified.\n");
    return 1; 
  }

  lmin=cmdparams_save_int(&cmdparams, "LMIN", &newstat);
  lmax=cmdparams_save_int(&cmdparams, "LMAX", &newstat);
  lchunksize=cmdparams_save_int(&cmdparams, "LCHUNK", &newstat);
  normflag=cmdparams_save_int(&cmdparams, "NORM", &newstat);

  tstart=cmdparams_save_time(&cmdparams, "TSTART", &newstat);
  sprint_time(tstartstr, tstart, "TAI", 0);
  ttotal=(char *)cmdparams_save_str(&cmdparams, "TTOTAL", &newstat);
  nseconds=atoinc(ttotal);

/*
  if (strcmp(kNOTSPECIFIED, ttotal))
  {
    nseconds=atoinc(ttotal);
  }
  else
    nseconds=0.0;
*/
/*
  status=drms_names_parseduration(&ttotal, &nseconds, 1);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr, "ERROR: problem parsing TTOTAL, = %s\n", ttotal);
    return 1; 
  }
*/
  tchunk=(char *)cmdparams_save_str(&cmdparams, "TCHUNK", &newstat);
  if (strcmp(kNOTSPECIFIED, tchunk))
  {
    chunksecs=atoinc(tchunk);
/*
    status=drms_names_parseduration(&tchunk, &chunksecs, 1);
    if (status != DRMS_SUCCESS)
      newstat = newstat | CPSAVE_UNKNOWN_ERROR;
*/
  }
  else if (!tsflag)
  {
    fprintf(stderr, "ERROR: TCHUNK must be specified if no tsout is given.\n");
    return 1; 
  }
  else
    chunksecs=0.0;


  subtract_mean = cmdparams_save_int(&cmdparams, "SUBMEAN", &newstat);
  normalize = cmdparams_save_int(&cmdparams, "NORMLIZE", &newstat);
  /* CENTLONG=1 centers the longitude Fourier transform on the center
     of the remapped image */
  cent_long = cmdparams_save_int(&cmdparams, "CENTLONG", &newstat);
  /* ZEROMISS=1 sets missing data to 0,
     ZEROMISS=0 fills the output row with missing */
  zero_miss = cmdparams_save_int(&cmdparams, "ZEROMISS", &newstat);
  lgapod = cmdparams_save_int(&cmdparams, "LGAPOD", &newstat);
  lgapmin = cmdparams_save_double(&cmdparams, "LGAPMIN", &newstat);
  lgapwidt = cmdparams_save_double(&cmdparams, "LGAPWIDT", &newstat);
  lgapmax = lgapmin+lgapwidt;

  int checko = cmdparams_save_int(&cmdparams, "OFLAG", &newstat);
  int NaN_beyond_rmax = cmdparams_save_int(&cmdparams, "NAN_BEYOND_RMAX", &newstat);
  int maxmissvals = cmdparams_save_int(&cmdparams, "MAXMISSVALS", &newstat);
  float rescale = params_get_float(&cmdparams, "RESCALE"); // YLiu

//  float beyondrmax = cmdparams_save_float(&cmdparams, "BEYONDRMAX", &newstat);

  carrStretch = cmdparams_save_int(&cmdparams, "CARRSTRETCH", &newstat);
  diffrotA = cmdparams_save_float(&cmdparams, "DIFROT_A", &newstat);
  diffrotB = cmdparams_save_float(&cmdparams, "DIFROT_B", &newstat);
  diffrotC = cmdparams_save_float(&cmdparams, "DIFROT_C", &newstat);

  longrate = 360.0 / TCARR - 360.0 / DAYSINYEAR; // degrees per day 
  longrate /= SECSINDAY; // degrees per sec 

  apodization = cmdparams_save_int(&cmdparams, "APODIZE", &newstat);
  apinner = cmdparams_save_double(&cmdparams, "APINNER", &newstat);  
  apwidth = cmdparams_save_double(&cmdparams, "APWIDTH", &newstat);
  apel = cmdparams_save_int(&cmdparams, "APEL", &newstat);
  apx = cmdparams_save_double(&cmdparams, "APX", &newstat);
  apy = cmdparams_save_double(&cmdparams, "APY", &newstat);
  longitude_shift = cmdparams_save_int(&cmdparams, "LGSHIFT", &newstat);
  mag_correction = cmdparams_save_int(&cmdparams, "MCORLEV", &newstat);
  mag_offset = cmdparams_save_int(&cmdparams, "MOFFSETFLAG", &newstat);
  velocity_correction = cmdparams_save_int(&cmdparams, "VCORLEV", &newstat);
  interpolation = cmdparams_save_int(&cmdparams, "INTERPO", &newstat);
  paramsign = cmdparams_save_int(&cmdparams, "DATASIGN", &newstat);
  rmax = cmdparams_save_double(&cmdparams, "MAPRMAX", &newstat);
  refl0 = cmdparams_save_double(&cmdparams, "REF_L0", &newstat);

  distsave = cmdparams_save_int(&cmdparams, "DISTORT", &newstat);
  cubsave = cmdparams_save_double(&cmdparams, "CUBIC", &newstat);
  tiltasave = cmdparams_save_double(&cmdparams, "TILTALPHA", &newstat);
  tiltbsave = cmdparams_save_double(&cmdparams, "TILTBETA", &newstat);
  tiltfsave = cmdparams_save_double(&cmdparams, "TILTFEFF", &newstat);

  scale = cmdparams_save_double(&cmdparams, "OUTSCALE", &newstat);
  bias = cmdparams_save_double(&cmdparams, "OUTBIAS", &newstat);
  p0 = cmdparams_save_double(&cmdparams, "SOLAR_P", &newstat);
  psign = cmdparams_save_double(&cmdparams, "PSIGN", &newstat);
  perr = cmdparams_save_double(&cmdparams, "PERR", &newstat);
  ierr = cmdparams_save_double(&cmdparams, "IERR", &newstat);
  trefb0 = cmdparams_save_time(&cmdparams, "REF_TB0", &newstat);

  remap_SetDistort(distsave, cubsave, tiltasave, tiltbsave, tiltfsave, &distP);

  tref = cmdparams_save_time(&cmdparams, "REF_T0", &newstat);

   // determine mapcols and adjust longmin and longmax */
  mapped_lmax = cmdparams_save_int(&cmdparams, "MAPMMAX", &newstat); 
  longmax = cmdparams_save_double(&cmdparams, "MAPLGMAX", &newstat); /* degrees */
  longmin = cmdparams_save_double(&cmdparams, "MAPLGMIN", &newstat); /* degrees */
  longinterval = (180.0) / mapped_lmax;	                     /* degrees */
 
   // This does not always handle the case where 1/longinterval is an integer correctly.

   // the next two statement do nothing, right? 
   // why do nmin and max keep getting set with different RHSs? 
  nmin = (int)(longmin / longinterval); // round towards 0 
  nmax = (int)(longmax / longinterval); // round towards 0 
  colsperdeg = mapped_lmax / 180.0;
  nmin = (int)(longmin * colsperdeg); // round towards 0 
  nmax = (int)(longmax * colsperdeg); // round towards 0 
  mapcols = nmax - nmin + 1;
  longmin_adjusted = nmin * longinterval;
  longmax_adjusted = nmax * longinterval;

   // determine maprows, bmax, bmin, and sinBdelta
  sinb_divisions = cmdparams_save_int(&cmdparams, "SINBDIVS", &newstat);
  sinBdelta = 1.0/sinb_divisions;
  bmax = cmdparams_save_double(&cmdparams, "MAPBMAX", &newstat);     // degrees
  bmin = -bmax; 
  nmax = (int)(sin(RADSINDEG*bmax)*sinb_divisions); // round towards 0 
  maprows = 2*nmax;

/* define vector B 
      xDim = 4096;
      yDim = 4096;
      bRadial = (float *) malloc(xDim * yDim * sizeof(float));
      bTheta = (float *) malloc(xDim * yDim * sizeof(float));
      bPhi = (float *) malloc(xDim * yDim * sizeof(float));
*/

  if (normflag == 0)
    cnorm=1.0;
  else
    cnorm = sqrt(2.)*sinBdelta;

  if (lmax > mapped_lmax || lmin > lmax)
  {
    fprintf(stderr, "ERROR: must have MAPMMAX >= LMAX >= LMIN, MAPMMAX = %d, LMAX= %d, LMIN = %d\n", mapped_lmax, lmax, lmin);
    return 1;
  }

  if (newstat) 
  {
    fprintf(stderr, "ERROR: problem with input arguments, status = %d, diagnosis follows\n", newstat);
    cpsave_decode_error(newstat);
    return 1;
  }  
  else if (savestrlen != strlen(savestr)) 
  {
    fprintf(stderr, "ERROR: problem with savestr, savestrlen = %d, strlen(savestr) = %d\n", savestrlen, (int)strlen(savestr));
    return 1;
  }

  DRMS_Record_t *tempoutrec;
  DRMS_Link_t *histlink;
  int itest;
/*
// cvsinfo used to be passed in the call to set_history. now this information is encoded in CVSTAG, which is defined by a compiler flag in the make.
  char *cvsinfo;
  cvsinfo = (char *)malloc(1024);
  strcpy(cvsinfo,"$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/vectmag2helio3comp_random.c,v 1.2 2016/05/17 16:31:31 yliu Exp $");
  strcat(cvsinfo,"\n");
  strcat(cvsinfo,getshtversion());
*/
// assume all output dataseries link to the same dataseries for HISTORY
  if (tsflag)
  {
    tempoutrec = drms_create_record(drms_env, outseries, DRMS_TRANSIENT, &status);
    if (status != DRMS_SUCCESS) 
    {
     fprintf(stderr,"ERROR: couldn't open a record in output dataseries %s, status = %d\n", outseries, status);
     return 1;
    }

    DRMS_Keyword_t *outkeytest = hcon_lookup_lower(&tempoutrec->keywords, "MAPMMAX");
    if (outkeytest != NULL && outkeytest->info->recscope == 1)
    {
      int mapmmaxout=drms_getkey_int(tempoutrec, "MAPMMAX", &status);
      if (mapmmaxout != mapped_lmax)
      {
       fprintf(stderr,"ERROR: output MAPMMAX=%d does not match input parameter MAPMMAX=%d, status = %d\n", mapmmaxout, mapped_lmax, status);
       return 1;
      }
    }

    outkeytest = hcon_lookup_lower(&tempoutrec->keywords, "SINBDIVS");
    if (outkeytest != NULL && outkeytest->info->recscope == 1)
    {
      int sinbdivsout=drms_getkey_int(tempoutrec, "SINBDIVS", &status);
      if (sinbdivsout != sinb_divisions)
      {
       fprintf(stderr,"ERROR: output SINBDIVS=%d does not match input parameter SINBDIVS=%d, status = %d\n", sinbdivsout, sinb_divisions, status);
       return 1;
      }
    }

// set up ancillary dataseries for processing metadata
    if (histflag)
    {
      histlink = hcon_lookup_lower(&tempoutrec->links, histlinkname);
      if (histlink != NULL) 
      {
        histrecnum=set_history(histlink);
        if (histrecnum < 0)
        {
          drms_close_record(tempoutrec, DRMS_FREE_RECORD);
          return 1;
        }
      }
      else
      {
        fprintf(stderr,"WARNING: could not find history link in output dataseries\n");
      }
    }

// these must be present in the output dataseries and variable, not links or constants   
    char *outchecklist[] = {"T_START", "QUALITY", "LMIN", "LMAX", "NDT"};
    for (itest=0; itest < ARRLENGTH(outchecklist); itest++)
    {
      DRMS_Keyword_t *outkeytest = hcon_lookup_lower(&tempoutrec->keywords, outchecklist[itest]);
      if (outkeytest == NULL || outkeytest->info->islink || outkeytest->info->recscope == 1)
      {
        fprintf(stderr, "ERROR: output keyword %s is either missing, constant, or a link\n", outchecklist[itest]);
        drms_close_record(tempoutrec, DRMS_FREE_RECORD);
        return 1;
      }
    }

    cadence0=drms_getkey_float(tempoutrec, "T_STEP", &status);
    tepoch=drms_getkey_time(tempoutrec, "T_START_epoch", &status);
    tstep=drms_getkey_float(tempoutrec, "T_START_step", &status);
    tround=drms_getkey_float(tempoutrec, "T_START_round", &status);
    if (fmod(tstart-tepoch,tstep) > tround/2)
    {
      sprint_time(trecstr, tepoch, "TAI", 0);
      fprintf(stderr, "ERROR: output dataseries seems incompatible with input parameters (tstep must divide tstart-tepoch): TSTART = %s, T_START_epoch = %s, tstep = %f\n", 
                                                                                                                            tstartstr, trecstr, tstep);
      drms_close_record(tempoutrec, DRMS_FREE_RECORD);
      return 1;
    }
    if (chunksecs == 0.0)
      chunksecs = tstep;
    else if (fmod(chunksecs,tstep))
    {
      fprintf(stderr, "ERROR: output dataseries seems incompatible with input parameters (tstep must divide chunksecs): chunksecs = %f, tstep = %f\n", chunksecs, tstep);
      drms_close_record(tempoutrec, DRMS_FREE_RECORD);
      return 1;
    }

    drms_close_record(tempoutrec, DRMS_FREE_RECORD);
  }

  if (v2hflag)
  {
    tempoutrec = drms_create_record(drms_env, v2hout, DRMS_TRANSIENT, &status);
    if (status != DRMS_SUCCESS) 
    {
     fprintf(stderr,"ERROR: couldn't open a record in output dataseries %s, status = %d\n", v2hout, status);
     return 1;
    }

// set up ancillary dataseries for processing metadata
    if (histflag && histrecnum < 0)
    {
      histlink = hcon_lookup_lower(&tempoutrec->links, histlinkname);
      if (histlink != NULL) 
      {
        histrecnum=set_history(histlink);
        if (histrecnum < 0)
        {
          drms_close_record(tempoutrec, DRMS_FREE_RECORD);
          return 1;
        }
      }
      else
      {
        fprintf(stderr,"WARNING: could not find history link in output dataseries\n");
      }
    }

// these must be present in the output dataseries and variable, not links or constants   
    char *outchecklist[] = {"T_REC", "QUALITY", "CRPIX1", "CRVAL1", "CDELT1", "CRPIX2", "CROTA2", "CDELT2" };
    for (itest=0; itest < ARRLENGTH(outchecklist); itest++)
    {
      DRMS_Keyword_t *outkeytest = hcon_lookup_lower(&tempoutrec->keywords, outchecklist[itest]);
      if (outkeytest == NULL || outkeytest->info->islink || outkeytest->info->recscope == 1)
      {
        fprintf(stderr, "ERROR: output keyword %s is either missing, constant, or a link\n", outchecklist[itest]);
        drms_close_record(tempoutrec, DRMS_FREE_RECORD);
        return 1;
      }
    }

    drms_close_record(tempoutrec, DRMS_FREE_RECORD);
  }

  if (h2mflag)
  {
    tempoutrec = drms_create_record(drms_env, h2mout, DRMS_TRANSIENT, &status);
    if (status != DRMS_SUCCESS) 
    {
     fprintf(stderr,"ERROR: couldn't open a record in output dataseries %s, status = %d\n", h2mout, status);
     return 1;
    }

// set up ancillary dataseries for processing metadata
    if (histflag && histrecnum < 0)
    {
      histlink = hcon_lookup_lower(&tempoutrec->links, histlinkname);
      if (histlink != NULL) 
      {
        histrecnum=set_history(histlink);
        if (histrecnum < 0)
        {
          drms_close_record(tempoutrec, DRMS_FREE_RECORD);
          return 1;
        }
      }
      else
      {
        fprintf(stderr,"WARNING: could not find history link in output dataseries\n");
      }
    }

// these must be present in the output dataseries and variable, not links or constants   
    char *outchecklist[] = {"T_REC", "QUALITY", "CRPIX1", "CDELT1", "CRPIX2", "CDELT2" };
    for (itest=0; itest < ARRLENGTH(outchecklist); itest++)
    {
      DRMS_Keyword_t *outkeytest = hcon_lookup_lower(&tempoutrec->keywords, outchecklist[itest]);
      if (outkeytest == NULL || outkeytest->info->islink || outkeytest->info->recscope == 1)
      {
        fprintf(stderr, "ERROR: output keyword %s is either missing, constant, or a link\n", outchecklist[itest]);
        drms_close_record(tempoutrec, DRMS_FREE_RECORD);
        return 1;
      }
    }

    drms_close_record(tempoutrec, DRMS_FREE_RECORD);
  }


  if (fmod(nseconds,chunksecs) != 0.0)
  {
    fprintf(stderr, "ERROR: input parameters seem incompatible (chunksecs must divide totalsecs): totalsecs = %f, chunksecs = %f\n", nseconds, chunksecs);
    return 1;
  }
  ntimechunks=nseconds/chunksecs;

  inrecset = drms_open_recordset(drms_env, inrecquery, &status);
  if (status != DRMS_SUCCESS || inrecset == NULL)
  {
    fprintf(stderr, "ERROR: problem opening input recordset: status = %d\n", status);
    return 1;
  }

  int nrecsin = drms_count_records(drms_env, inrecquery, &status);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr, "ERROR: problem counting input records: status = %d, nrecs = %d\n", status, nrecsin);
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }
  if (nrecsin == 0)
  {
    fprintf(stderr, "ERROR: input recordset contains no records. if such was intended use jretile instead.\n");
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }

//the above replaces the following.  drms_open_recordset() no longer fills in the number of records.
/*
  if (inrecset->n == 0)
  {
    fprintf(stderr, "ERROR: input recordset contains no records. if such was intended use jretile instead.\n");
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }
*/

  if (verbflag) 
    printf("input recordset opened, nrecs = %d\n", nrecsin);

  inrec = drms_recordset_fetchnext(drms_env, inrecset, &fetchstat, &chunkstat, NULL);

// these must be present in the input dataseries
  char *inchecklist[] = {"T_REC", "QUALITY", "T_OBS", "CRLT_OBS", "CRLN_OBS", "CADENCE",
                     //     "SAT_ROT", "INST_ROT", "IM_SCALE",
                          "CDELT1", "CDELT2"};

  DRMS_Keyword_t *inkeytest;
  for (itest=0; itest < ARRLENGTH(inchecklist); itest++)
  {
    inkeytest = hcon_lookup_lower(&inrec->keywords, inchecklist[itest]);
    if (inkeytest == NULL)
    {
      fprintf(stderr, "ERROR: required input keyword %s is missing\n", inchecklist[itest]);
      drms_close_records(inrecset, DRMS_FREE_RECORD);
      return 1;
    }
  }

  int readrsunref=0;
  double rsunref;
  inkeytest = hcon_lookup_lower(&inrec->keywords, "RSUN_REF");
  if (inkeytest == NULL)
    rtrue = RTRUE/AU;
  else if (inkeytest->info->recscope == 1)
  {
    rsunref = drms_getkey_double(inrec, "RSUN_REF", &status);
    ParameterDef(status, "RSUN_REF", RTRUE, &rsunref, inrec->recnum, 1);
    rtrue=rsunref/AU;
  }
  else
    readrsunref=1;

  trec = drms_getkey_time(inrec, "T_REC", &status);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr, "ERROR: problem with required parameter T_REC: status = %d, recnum = %lld\n", status, inrec->recnum);
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }
  sprint_time(trecstr, trec, "TAI", 0);

  cadence=drms_getkey_float(inrec, "CADENCE", &status);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr, "ERROR: problem with required parameter CADENCE: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }

  if (!forceoutput)
  {
    if (nrecsin != nseconds/cadence)
    {
      fprintf(stderr, "ERROR: input recordset does not contain a record for every slot.\n");
      drms_close_records(inrecset, DRMS_FREE_RECORD);
      return 1;
    }
  }

  if (tsflag && cadence != cadence0)
  {
    fprintf(stderr, "ERROR: input CADENCE does not match output T_STEP: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }

  nrecs=chunksecs/cadence;
  maxnsn=nsn=nrecs;
  fournsn=4*nsn;

  if (verbflag) 
    printf("ntimechunks = %d, recs per timechunk = %d\n", ntimechunks, nrecs);


  if (trec >= tstart + nseconds)
  {
    fprintf(stderr, "ERROR: no records processed: first input record is after last output record: T_REC = %s\n", trecstr);
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }

  while (trec < tstart && chunkstat != kRecChunking_NoMoreRecs)
  {
    fprintf(stderr, "WARNING: input record will not be included in output: T_REC = %s, TSTART = %s \n", trecstr, tstartstr);
    inrec = drms_recordset_fetchnext(drms_env, inrecset, &fetchstat, &chunkstat, NULL);
    if (inrec != NULL)
    {
      trec = drms_getkey_time(inrec, "T_REC", &status);
      sprint_time(trecstr, trec, "TAI", 0);
    }
  }
  if (chunkstat == kRecChunking_NoMoreRecs)
  {
    fprintf(stderr,"ERROR: no records processed: last input record is before first output record: T_REC = %s\n", trecstr);
    drms_close_records(inrecset, DRMS_FREE_RECORD);
    return 1;
  }

  msize = lmax+1;
  if (lchunksize == 0) lchunksize = msize;

  nlat = maprows/2;
  imagesize = maprows*2*msize; /* out image could be smaller than in */

  nout = 2 * (lmax + 1);

  lfft = 2 * mapped_lmax;
  nfft = lfft + 2;
  mapcols2 = mapcols/2;

/* Let's try this since SGI's like odd leading dimensions of the first
   array in sgemm */
//  nlatx=2*(nlat/2)+1;

/* make nlatx divisible by 4 on linux systems */
//#ifdef __linux__
  if (nlat % 4)
    nlatx=4*(nlat/4+1);
  else
    nlatx=nlat;
//#endif

  DRMS_RecordSet_t *outrecset, *v2hrecset, *h2mrecset;
  if (tsflag)
  {
    real4evenl = (float *)(malloc (nlatx*maxnsn*sizeof(float)));
    real4oddl = (float *)(malloc (nlatx*maxnsn*sizeof(float)));
    imag4evenl = (float *)(malloc (nlatx*maxnsn*sizeof(float)));
    imag4oddl = (float *)(malloc (nlatx*maxnsn*sizeof(float)));
    outx = (float *)(malloc (maxnsn*2*lchunksize*sizeof(float)));

    plm = (double *)(malloc ((lmax+1)*nlat*sizeof(double)));  
    saveplm = (double *)(malloc ((lmax+1)*nlat*2*sizeof(double)));
    latgrid = (double *)(malloc (nlat*sizeof(double)));
    for (i=0; i<nlat; i++) latgrid[i] = (i+0.5)*sinBdelta;

    indx = (int *)(malloc ((lmax+1)*sizeof(int)));
    for (l=0; l<=lmax; l++) indx[l]=l;

    masks = (float *)(malloc (nlat*lchunksize*sizeof(float)));
    foldedsize = 4*nlat*(lmax+1)*maxnsn;
    folded = (float *)(malloc (foldedsize*sizeof(float)));
    oddpart = (float *)(malloc (nlat*sizeof(float)));
    evenpart = (float *)(malloc (nlat*sizeof(float)));

    lchunkfirst = lmin/lchunksize;
    lchunklast = lmax/lchunksize;
    int nlchunks = (lchunklast - lchunkfirst) + 1;
    int nrecsout = nlchunks*ntimechunks;
    outrecset = drms_create_records(drms_env, nrecsout, outseries, lifetime, &status);
    if (status != DRMS_SUCCESS || outrecset == NULL)
    {
      fprintf(stderr,"ERROR: unable to create records in output dataseries %s, status = %d\n", outseries, status);
      drms_close_records(inrecset, DRMS_FREE_RECORD);
      return 1;
    }

  }

  if (v2hflag)
  {
    v2hrecset = drms_create_records(drms_env, nrecsin, v2hout, lifetime, &status);
    if (status != DRMS_SUCCESS || v2hrecset == NULL)
    {
      fprintf(stderr,"ERROR: unable to create records in output dataseries %s, status = %d\n", v2hout, status);
      drms_close_records(inrecset, DRMS_FREE_RECORD);
      return 1;
    }
  }

  if (h2mflag)
  {
    h2mrecset = drms_create_records(drms_env, nrecsin, h2mout, lifetime, &status);
    if (status != DRMS_SUCCESS || h2mrecset == NULL)
    {
      fprintf(stderr,"ERROR: unable to create records in output dataseries %s, status = %d\n", h2mout, status);
      drms_close_records(inrecset, DRMS_FREE_RECORD);
      return 1;
    }
  }


  bad = (int *)(malloc (nrecs*sizeof(int)));
  v2hbr = (float *)(malloc(maprows*mapcols*sizeof(float)));
  v2hbt = (float *)(malloc(maprows*mapcols*sizeof(float)));
  v2hbp = (float *)(malloc(maprows*mapcols*sizeof(float)));
  h2mptr = (float *)(malloc(maprows*nout*sizeof(float)));

    /* get working buffer */
  buf = (float *)malloc(nfft * sizeof(float));

  wbuf = (float *)malloc(nfft * sizeof(float));
  fftwp = fftwf_plan_r2r_1d(lfft, buf, wbuf, FFTW_R2HC, FFTW_ESTIMATE);

    /* get weight array for apodizing */
  weight = (float *)malloc(nfft * sizeof(float));

  wp = weight;
  for (col=0; col<mapcols; col++) 
  {
    if (lgapod) 
    {
      lon=abs(col-mapcols2)*360.0/(2*mapped_lmax);
      if (lon < lgapmin)
        *wp++=1.0;
      else if (lon < lgapmax)
        *wp++=0.5+0.5*cos(PI*(lon-lgapmin)/lgapwidt);
      else
        *wp++=0.0;
    }
    else 
      *wp++ = 1.0;
  }

  status=drms_stage_records(inrecset, 1, 0);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr, "ERROR: drms_stage_records returned status = %d\n", status);
    return 1;
  }

  unsigned long long calversout, calvers;
  int calversunset=1;
  unsigned int nybblearrout[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int fixflagarr[16]            = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for (i=0;i<16;i++)
  {
    if (getbits(calverkey,i,1))
      fixflagarr[i]=1;
  }

  int mixflag=0;

  int nrecords=0;
  int nsegments=0;
  int error=0;
  int nodata=0;
  int irecout=0;
  int iv2hrec=0;
  int ih2mrec=0;
  for (iset=0; iset < ntimechunks && chunkstat != kRecChunking_NoMoreRecs; iset++)
  {
    sprint_time(tstartstr, tstart, "TAI", 0);
    if (verbflag)
    {
      wt1=getwalltime();
      ct1=getcputime(&ut1, &st1);
      printf("processing timechunk %d, tstart = %s\n", iset, tstartstr);
    }

    if (trec >= tstart+chunksecs)
    {
      if (forceoutput)
      {
        nodata=1;
        goto skip_norecs;
      }
      else
      {
        fprintf(stderr, "ERROR: no data for timechunk beginning at %s\n", tstartstr);
        error++;
        tstart+=chunksecs;
        continue;
      }
    }

    while (trec < tstart && chunkstat != kRecChunking_NoMoreRecs)
    {
      inrec = drms_recordset_fetchnext(drms_env, inrecset, &fetchstat, &chunkstat, NULL);
      if (inrec != NULL)
      {
        trec = drms_getkey_time(inrec, "T_REC", &status);
        sprint_time(trecstr, trec, "TAI", 0);
      }
    }

    bc=0;
    nodata=0;
    mixflag=0;
    calversunset=1;
    trecind=(trec-tstart+cadence/2)/cadence;

    for (irec=0; irec < nrecs && chunkstat != kRecChunking_NoMoreRecs; irec++)
    {

      if (trecind > irec)
//some inputs were missing
      {
        if (forceoutput)
        {
          bad[bc++]=irec;
          continue;
        }
        else
        {
          fprintf(stderr, "ERROR: some input records missing, T_START = %s, T_REC = %s, irec = %d\n", tstartstr, trecstr, irec);
          error++;
          goto continue_outer_loop;
        }
      }

      while (trecind < irec)
//T_REC is duplicated in input
      {
        if (forceoutput)
        {
          inrec = drms_recordset_fetchnext(drms_env, inrecset, &fetchstat, &chunkstat, NULL);
          if (inrec != NULL)
          {
            trec = drms_getkey_time(inrec, "T_REC", &status);
            sprint_time(trecstr, trec, "TAI", 0);
            trecind=(trec-tstart+cadence/2)/cadence;
          }
        }
        else
        {
          fprintf(stderr, "ERROR: some input records have duplicate T_REC, T_START = %s, T_REC = %s, irec = %d\n", tstartstr, trecstr, irec);
          error++;
          goto continue_outer_loop;
        }
      }

      if (verbflag > 1) 
      {
        wt2=getwalltime();
        ct2=getcputime(&ut2, &st2);
        printf("  processing record %d\n", irec);
      }

      quality=drms_getkey_int(inrec, "QUALITY", &status);
      if (status != DRMS_SUCCESS || (quality & QUAL_NODATA)) //may want stricter test on quality here
      {
        bad[bc++]=irec;
        if (verbflag > 2)
          fprintf(stderr, "SKIP: image rejected based on quality: T_REC = %s, quality = %08x\n", trecstr, quality);
        goto skip;
      }

      if (tsflag)
      {
        if (calversunset)
        {
          calversout=drms_getkey_longlong(inrec, "CALVER64", &status);
          if (status != DRMS_SUCCESS)
            calversout = 0;
          else
            calversout = fixcalver64(calversout);
          calversunset=0;

          for (i=0;i<16;i++)
            nybblearrout[i]=getbits(calversout,4*i+3,4);

        }

        calvers=drms_getkey_longlong(inrec, "CALVER64", &status);
        if (status != DRMS_SUCCESS)
          calvers = 0;
        else
          calvers = fixcalver64(calvers);

        for (i=0;i<16;i++)
        {
          int nybble=getbits(calvers,4*i+3,4);
          if (fixflagarr[i])
          {
            if (nybble != nybblearrout[i])
            {
              fprintf(stderr, "ERROR: input data has mixed values for field %d of CALVER64: %d and %d, recnum = %lld, histrecnum = %lld\n", i, nybblearrout[i], nybble, inrec->recnum, histrecnum);
              error++;
              goto continue_outer_loop;
            }
          }
          else
          {
            if (nybble < nybblearrout[i])
              nybblearrout[i]=nybble;
          }
        }

        if (!mixflag && calvers != calversout)
          mixflag=1;
      }

/*
      if (seginflag)
        segin = drms_segment_lookup(inrec, segnamein);
      else
        segin = drms_segment_lookupnum(inrec, 0);
*/

// read vector B data to *bTotal, *bIncl, *bAzim, *disamb

      segin = drms_segment_lookup(inrec, "disambig");
      if (segin != NULL)
        inarrBdisamb = drms_segment_read(segin, usetype, &status);

      if (segin == NULL || inarrBdisamb == NULL || status != DRMS_SUCCESS)
      {
        fprintf(stderr, "ERROR: problem with input segment or array: status = %d, T_REC = %s, recnum = %lld, histrecnum = %lld\n", status, trecstr, inrec->recnum, histrecnum);
        return 0;
      }
      disamb = (float *)inarrBdisamb->data;

      segin = drms_segment_lookup(inrec, "field");
      if (segin != NULL)
        inarrBtotal = drms_segment_read(segin, usetype, &status);

      if (segin == NULL || inarrBtotal == NULL || status != DRMS_SUCCESS)
      {
        fprintf(stderr, "ERROR: problem with input segment or array: status = %d, T_REC = %s, recnum = %lld, histrecnum = %lld\n", status, trecstr, inrec->recnum, histrecnum);
        drms_free_array(inarrBdisamb);
        return 0;
      }
      bTotal = (float *)inarrBtotal->data;

      segin = drms_segment_lookup(inrec, "inclination");
      if (segin != NULL)
        inarrBincl = drms_segment_read(segin, usetype, &status);

      if (segin == NULL || inarrBincl == NULL || status != DRMS_SUCCESS)
      {
        fprintf(stderr, "ERROR: problem with input segment or array: status = %d, T_REC = %s, recnum = %lld, histrecnum = %lld\n", status, trecstr, inrec->recnum, histrecnum);
        drms_free_array(inarrBdisamb);
        drms_free_array(inarrBtotal);
        return 0;
      }
      bIncl = (float *)inarrBincl->data;

      segin = drms_segment_lookup(inrec, "azimuth");
      if (segin != NULL)
        inarrBazim = drms_segment_read(segin, usetype, &status);

      if (segin == NULL || inarrBazim == NULL || status != DRMS_SUCCESS)
      {
        fprintf(stderr, "ERROR: problem with input segment or array: status = %d, T_REC = %s, recnum = %lld, histrecnum = %lld\n", status, trecstr, inrec->recnum, histrecnum);
        return 0;
        drms_free_array(inarrBdisamb);
        drms_free_array(inarrBtotal);
        drms_free_array(inarrBincl);
      }
      bAzim = (float *)inarrBazim->data;

// end of reading vector B

      if (maxmissvals > 0) 
      {
        int missvals = drms_getkey_int(inrec, "MISSVALS", &status);
        PARAMETER_ERROR("MISSVALS")
        if (missvals > maxmissvals) 
        {
          bad[bc++]=irec;
          if (verbflag > 1)
            fprintf(stderr, "SKIP: %d pixels MISSING, max allowed is %d: T_REC = %s, recnum = %lld\n", missvals, maxmissvals, trecstr, inrec->recnum);
          drms_free_array(inarrBazim);
          drms_free_array(inarrBdisamb);
          drms_free_array(inarrBtotal);
          drms_free_array(inarrBincl);
          goto skip;
        }
      }

      tobs = drms_getkey_time(inrec, "T_OBS", &status);
      PARAMETER_ERROR("T_OBS")

      // MDI keyword was OBS_B0
      b0 = drms_getkey_double(inrec, "CRLT_OBS", &status);
      PARAMETER_ERROR("CRLT_OBS")

      // MDI keyword was OBS_L0
      obsl0 = drms_getkey_double(inrec, "CRLN_OBS", &status);
      PARAMETER_ERROR("CRLN_OBS")

      if (p0 == 999.0) 
      {
         // MDI keyword was SOLAR_P = -(SAT_ROT + INST_ROT)
/*
         satrot = drms_getkey_double(inrec, "SAT_ROT", &status);
         PARAMETER_ERROR("SAT_ROT")
         instrot = drms_getkey_double(inrec, "INST_ROT", &status);
         PARAMETER_ERROR("INST_ROT")
         p=-(satrot+instrot);
*/
         double crota = drms_getkey_double(inrec, "CROTA2", &status);
         PARAMETER_ERROR("CROTA2")
         p=-crota;
      } 
      else 
      {
         p = p0;
      }

      p = psign * p ;
      b0 = b0 + ierr * sin((tobs - trefb0) / 31557600. * 2 * PI);
      p = p + perr - ierr * cos((tobs - trefb0) / 31557600. * 2 * PI);

      // S_MAJOR, S_MINOR, S_ANGLE, X_SCALE, Y_SCALE were MDI keywords, not used for HMI, but still necessary for GONG data
      smajor = drms_getkey_double(inrec, "S_MAJOR", &status);
      ParameterDef(status, "S_MAJOR", 1.0, &smajor, inrec->recnum, 0);

      sminor = drms_getkey_double(inrec, "S_MINOR", &status);
      ParameterDef(status, "S_MINOR", 1.0, &sminor, inrec->recnum, 0);

      sangle = drms_getkey_double(inrec, "S_ANGLE", &status);
      ParameterDef(status, "S_ANGLE", 0.0, &sangle, inrec->recnum, 0);

      /*
      our calculation of CDELTi does not follow WCS conventions.  it should be the plate scale
      at the reference pixel (disk center), but instead we use the average between center and limb.
      this is taken into account in the calculation of rsun below.
      */
      xscale = drms_getkey_double(inrec, "CDELT1", &status);
      PARAMETER_ERROR("CDELT1")
      yscale = drms_getkey_double(inrec, "CDELT2", &status);
      PARAMETER_ERROR("CDELT2")

      // use xscale and yscale for the following check, then set to 1.0 for the call to obs2helio
      if (xscale != yscale)
      {
        fprintf(stderr, "ERROR: CDELT1 != CDELT2 not supported, CDELT1 = %f, CDELT2 = %f: T_REC = %s, recnum = %lld \n", xscale, yscale, trecstr, inrec->recnum);
        drms_free_array(inarrBazim);
        drms_free_array(inarrBdisamb);
        drms_free_array(inarrBtotal);
        drms_free_array(inarrBincl);
        error++;
        goto continue_outer_loop;
      }
      imagescale=xscale;
      xscale=1.0;
      yscale=1.0;

/*
      imagescale = drms_getkey_double(inrec, "IM_SCALE", &status);
      PARAMETER_ERROR("IM_SCALE")
*/

      if (paramsign != 0)
      {
         vsign = paramsign;
      }
      else
      {
         vsign = drms_getkey_double(inrec, "DATASIGN", &status);
         ParameterDef(status, "DATASIGN", 1.0, &vsign, inrec->recnum, 1);
      }

      if (velocity_correction) 
      {
         obs_vr = drms_getkey_double(inrec, "OBS_VR", &status);
         ParameterDef(status, "OBS_VR", 0.0, &obs_vr, inrec->recnum, 1);

         obs_vw = drms_getkey_double(inrec, "OBS_VW", &status);
         ParameterDef(status, "OBS_VW", 0.0, &obs_vw, inrec->recnum, 1);

         obs_vn = drms_getkey_double(inrec, "OBS_VN", &status);
         ParameterDef(status, "OBS_VN", 0.0, &obs_vn, inrec->recnum, 1);
      }

      // MDI keyword was OBS_DIST, in AU
      obsdist = drms_getkey_double(inrec, "DSUN_OBS", &status) / AU;
      // note that an incorrect value of 1.49597892e11 has sometimes been used to convert between OBS_DIST and DSUN_OBS when porting data from DSDS to DRMS 
      ParameterDef(status, "OBS_DIST", 1.0, &obsdist, inrec->recnum, 1);
      if (readrsunref)
      {
        rsunref = drms_getkey_double(inrec, "RSUN_REF", &status);
        ParameterDef(status, "RSUN_REF", RTRUE, &rsunref, inrec->recnum, 1);
        rtrue=rsunref/AU;
      }
      S = rtrue / obsdist; // radians - approx. arcsin(rtrue/obsdist), but don't undo this approximation, because it is assumed in obs2helio

      rsunobs = drms_getkey_double(inrec, "RSUN_OBS", &status);
      if (status == DRMS_SUCCESS) 
        rsun = rsunobs/imagescale;  //this calculation of rsun assumes approximation of imagescale mentioned in comment above
      else
      {
        rsun = drms_getkey_double(inrec, "R_SUN", &status);
        if (status != DRMS_SUCCESS)
          rsun = ARCSECSINRAD * S / sqrt(1.0 - S * S) / imagescale;	   
      }

      if (longitude_shift == 1) 
      {
         tmearth = tobs+TAU_A*(1.0-obsdist); 
         longshift = (obsl0-refl0)+longrate*(tmearth-tref); // degrees 
         while (longshift > 180.0) longshift-=360.0; 
         while (longshift < -180.0) longshift+=360.0;
      }
      else if (longitude_shift == 2) // Shift center to nearest Carrington Degree 
      {
         longshift =  obsl0 - (int)(obsl0);
         if (longshift > 0.5) longshift -= 1.0;
      }
      else if (longitude_shift == 3)  // Shift center to nearest tenth of a degree 
      {
         longshift = (obsl0 * 10 - (int)(obsl0 * 10)) / 10;
         if (longshift > 0.5) longshift -= 1.0;
      }
      else
      {
         longshift = 0.0;
      }

      mapl0 = obsl0 - longshift;

      xpixels = inarrBtotal->axis[0];
      ypixels = inarrBtotal->axis[1];
      pixels  = xpixels * ypixels;

      // MDI keyword was X0
      x0 = drms_getkey_double(inrec, "CRPIX1", &status);
      ParameterDef(status, "CRPIX1", xpixels / 2, &x0, inrec->recnum, 1);
      x0 -= 1.0;

      // MDI keyword was Y0
      y0 = drms_getkey_double(inrec, "CRPIX2", &status);
      ParameterDef(status, "CRPIX2", ypixels / 2, &y0, inrec->recnum, 1);
      y0 -= 1.0;

      if (mag_offset) 
      {
        float *dat = (float *)inarrBtotal->data;
        double bfitzero = drms_getkey_double(inrec, "BFITZERO", &status);
        PARAMETER_ERROR("BFITZERO")
        int i;

        if (!isnan(bfitzero)) 
        {
          for (i = 0; i < pixels; ++i) 
          {
            dat[i] -= (float)bfitzero;
          }

        }
      }

      if (checko)
      {
        orientation = drms_getkey_string(inrec, "ORIENT", &status);
        PARAMETER_ERROR("ORIENT")
        CheckO(orientation, &vstat);
        if (vstat != V2HStatus_Success)
        { 
          fprintf(stderr,"ERROR: illegal ORIENT: T_REC = %s, recnum = %lld\n", trecstr, inrec->recnum);
          drms_free_array(inarrBazim);
          drms_free_array(inarrBdisamb);
          drms_free_array(inarrBtotal);
          drms_free_array(inarrBincl);
          free(orientation);
          error++;
          goto continue_outer_loop;
        }
      }
      else
      {
         orientation = strdup(orientationdef);
      }

// convert Btotal, inclination and azimuth to Br, Btheta and Bphi


      xDim = inarrBtotal->axis[0];
      yDim = inarrBtotal->axis[1];
      double bxHel, byHel, bzHel;
      float *bRadial, *bTheta, *bPhi;
      bRadial = (float *)malloc(xDim * yDim * sizeof(float));
      bTheta = (float *)malloc(xDim * yDim * sizeof(float));
      bPhi = (float *)malloc(xDim * yDim * sizeof(float));

        double xx, yy, lon, lat, coslat, sinlat;
        double mu, rho, sig, chi;
        jy = 0; iData = 0; yOff = 0;
        for (jy = 0; jy < yDim; jy++)
          {
            ix = 0;
            yy = (double)jy - y0;
            yy /= rsun;
            yOff = jy * xDim;

            for (ix = 0; ix < xDim; ix++)
              {
                iData = yOff + ix;
                xx = (double)ix - x0;
                xx /= rsun;
                if (isnan(bTotal[iData]))
                  {
                    bRadial[iData] = DRMS_MISSING_FLOAT;
                    bTheta[iData] = DRMS_MISSING_FLOAT;
                    bPhi[iData] = DRMS_MISSING_FLOAT;
                    continue;
                  }

                  synop_img2sphere (xx, yy, asin(S), b0 * RADSINDEG, obsl0 * RADSINDEG, p * RADSINDEG, &rho, &lat, &lon,
                    &sinlat, &coslat, &sig, &mu, &chi);
/*
                if (img2sphere (xx, yy, asin(S), b0 * RADSINDEG, obsl0 * RADSINDEG, p * RADSINDEG, &rho, &lat, &lon,
                    &sinlat, &coslat, &sig, &mu, &chi))
                  {
                    bRadial[iData] = DRMS_MISSING_FLOAT;
                    bTheta[iData] = DRMS_MISSING_FLOAT;
                    bPhi[iData] = DRMS_MISSING_FLOAT;
                    continue;
                  }
*/
                double bx = 0.0, by = 0.0, bz = 0.0;
//                if (disamb[iData] >= 4.0) bAzim[iData] += 180.0; //-- radial acute 
                if ((int)(disamb[iData]/2)%2 == 1) bAzim[iData] += 180.0; //-- random assumption
//                if ((int)(disamb[iData])%2 == 1) bAzim[iData] += 180.0; // -- potential field solution
                bx = -bTotal[iData] * sin(bIncl[iData] * RADSINDEG)
                     * sin(bAzim[iData] * RADSINDEG);
                by = bTotal[iData] * sin(bIncl[iData] * RADSINDEG)
                     * cos(bAzim[iData] * RADSINDEG);
                bz = bTotal[iData] * cos(bIncl[iData] * RADSINDEG);
                     // Azimuth angle is defined here to increase counter-clockwisely. The zero angle points 
                     // to the right. 
                     // transform the magnetic vector from image coordinates to the heliographic coordinates.
                remap_img2helioVector (bx, by, bz, &bxHel, &byHel, &bzHel, lon, lat, obsl0 * RADSINDEG, b0 * RADSINDEG, p * RADSINDEG);
                bPhi[iData] = bxHel;
                bTheta[iData] = -byHel;
                bRadial[iData] = bzHel;
                }
            } 

//      if (status =       obs2helio((float *)inarr->data, 
        if (status =       remap_obs2helio(bRadial,
	  			   v2hbr,
				   xpixels, 
				   ypixels, 
				   x0, 
				   y0, 
				   b0 * RADSINDEG, 
				   p * RADSINDEG, 
				   S, 
				   rsun, 
				   rmax,
				   interpolation, 
				   mapcols, 
				   maprows, 
				   longmin_adjusted * RADSINDEG, 
				   longinterval * RADSINDEG,
				   longshift * RADSINDEG, 
				   sinBdelta, 
				   smajor, 
				   sminor, 
				   sangle * RADSINDEG, 
				   xscale, 
				   yscale, 
				   orientation, 
				   mag_correction,
				   velocity_correction, 
				   obs_vr, 
				   obs_vw, 
				   obs_vn, 
				   vsign, 
				   NaN_beyond_rmax,
				   carrStretch,
				   &distP,
				   diffrotA,
				   diffrotB,
				   diffrotC,
				   NULL,
				   0))
      {
        fprintf(stderr, "ERROR: failure in obs2helio: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
        drms_free_array(inarrBazim);
        drms_free_array(inarrBdisamb);
        drms_free_array(inarrBtotal);
        drms_free_array(inarrBincl);
        free(bPhi); free(bTheta); free(bRadial);
        free(orientation);
        error++;
        goto continue_outer_loop;
      }

      if (status =      remap_apodize(v2hbr,
				 b0 * RADSINDEG, 
				 mapcols, 
				 maprows,
				 longmin_adjusted * RADSINDEG,
				 longinterval * RADSINDEG,
				 sinBdelta,
				 apodization, 
				 apinner, 
				 apwidth, 
				 apel, 
				 apx, 
				 apy)) 
      { 
        fprintf(stderr, "ERROR: failure in apodize: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
        error++;
        goto continue_outer_loop;
      }


        if (status =       remap_obs2helio(bTheta,
                                   v2hbt,
                                   xpixels,
                                   ypixels,
                                   x0,
                                   y0,
                                   b0 * RADSINDEG,
                                   p * RADSINDEG,
                                   S,
                                   rsun,
                                   rmax,
                                   interpolation,
                                   mapcols,
                                   maprows,
                                   longmin_adjusted * RADSINDEG,
                                   longinterval * RADSINDEG,
                                   longshift * RADSINDEG,
                                   sinBdelta,
                                   smajor,
                                   sminor,
                                   sangle * RADSINDEG,
                                   xscale,
                                   yscale,
                                   orientation,
                                   mag_correction,
                                   velocity_correction,
                                   obs_vr,
                                   obs_vw,
                                   obs_vn,
                                   vsign,
                                   NaN_beyond_rmax,
                                   carrStretch,
                                   &distP,
                                   diffrotA,
                                   diffrotB,
                                   diffrotC,
                                   NULL,
                                   0))
      {
        fprintf(stderr, "ERROR: failure in obs2helio: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
        drms_free_array(inarrBazim);
        drms_free_array(inarrBdisamb);
        drms_free_array(inarrBtotal);
        drms_free_array(inarrBincl);
        free(bPhi); free(bTheta); free(bRadial);
        free(orientation);
        error++;
        goto continue_outer_loop;
      }

      if (status =      remap_apodize(v2hbt,
                                 b0 * RADSINDEG,
                                 mapcols,
                                 maprows,
                                 longmin_adjusted * RADSINDEG,
                                 longinterval * RADSINDEG,
                                 sinBdelta,
                                 apodization,
                                 apinner,
                                 apwidth,
                                 apel,
                                 apx,
                                 apy))
      {
        fprintf(stderr, "ERROR: failure in apodize: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
        error++;
        goto continue_outer_loop;
      }

        if (status =       remap_obs2helio(bPhi,
                                   v2hbp,
                                   xpixels,
                                   ypixels,
                                   x0,
                                   y0,
                                   b0 * RADSINDEG,
                                   p * RADSINDEG,
                                   S,
                                   rsun,
                                   rmax,
                                   interpolation,
                                   mapcols,
                                   maprows,
                                   longmin_adjusted * RADSINDEG,
                                   longinterval * RADSINDEG,
                                   longshift * RADSINDEG,
                                   sinBdelta,
                                   smajor,
                                   sminor,
                                   sangle * RADSINDEG,
                                   xscale,
                                   yscale,
                                   orientation,
                                   mag_correction,
                                   velocity_correction,
                                   obs_vr,
                                   obs_vw,
                                   obs_vn,
                                   vsign,
                                   NaN_beyond_rmax,
                                   carrStretch,
                                   &distP,
                                   diffrotA,
                                   diffrotB,
                                   diffrotC,
                                   NULL,
                                   0))
      {
        fprintf(stderr, "ERROR: failure in obs2helio: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
        drms_free_array(inarrBazim);
        drms_free_array(inarrBdisamb);
        drms_free_array(inarrBtotal);
        drms_free_array(inarrBincl);
        free(bPhi); free(bTheta); free(bRadial);
        free(orientation);
        error++;
        goto continue_outer_loop;
      }

      drms_free_array(inarrBazim);
      drms_free_array(inarrBdisamb);
      drms_free_array(inarrBtotal);
      drms_free_array(inarrBincl);
      free(bPhi); free(bTheta); free(bRadial);
      free(orientation);

      if (status =      remap_apodize(v2hbp,
                                 b0 * RADSINDEG,
                                 mapcols,
                                 maprows,
                                 longmin_adjusted * RADSINDEG,
                                 longinterval * RADSINDEG,
                                 sinBdelta,
                                 apodization,
                                 apinner,
                                 apwidth,
                                 apel,
                                 apx,
                                 apy))
      {
        fprintf(stderr, "ERROR: failure in apodize: status = %d, T_REC = %s, recnum = %lld\n", status, trecstr, inrec->recnum);
        error++;
        goto continue_outer_loop;
      }

      if (v2hflag)
      {
//        outrec = drms_create_record(drms_env,  v2hout, lifetime,  &status);
        outrec=v2hrecset->records[iv2hrec++];
        drms_copykeys(outrec, inrec, 0, kDRMS_KeyClass_Explicit);
        DRMS_Link_t *histlink = hcon_lookup_lower(&outrec->links, histlinkname);
        DRMS_Link_t *srclink = hcon_lookup_lower(&outrec->links, srclinkname);
        if (histlink != NULL)
          drms_setlink_static(outrec, histlinkname,  histrecnum);
        if (srclink != NULL)
          drms_setlink_static(outrec, srclinkname,  inrec->recnum);
        if (segoutflag)
          segout = drms_segment_lookup(outrec, segnameout);
        else
          segout = drms_segment_lookupnum(outrec, 0);

// -- YLiu
        int mapcols_out = mapcols * rescale + 0.5;
        int maprows_out = maprows * rescale + 0.5;
        float *br_out, *bt_out, *bp_out;
        
        if (!(br_out = (float *) malloc(mapcols_out*maprows_out*4))) DIE("MALLOC_FAILED");
        if (!(bt_out = (float *) malloc(mapcols_out*maprows_out*4))) DIE("MALLOC_FAILED");
        if (!(bp_out = (float *) malloc(mapcols_out*maprows_out*4))) DIE("MALLOC_FAILED");

        do_boxcar(v2hbr, br_out, mapcols, maprows, rescale); // YLiu
        do_boxcar(v2hbt, bt_out, mapcols, maprows, rescale); // YLiu
        do_boxcar(v2hbp, bp_out, mapcols, maprows, rescale); // YLiu

        length[0]=mapcols_out;
        length[1]=maprows_out;

//        outarr = drms_array_create(usetype, 2, length, v2hbr, &status);
        outarr = drms_array_create(usetype, 2, length, br_out, &status); // YLiu
//        drms_setkey_int(outrec, "TOTVALS", maprows*mapcols); // YLiu
        drms_setkey_int(outrec, "TOTVALS", maprows_out*mapcols_out);
        set_statistics(segout, outarr, 1);
        outarr->bzero=segout->bzero;
        outarr->bscale=segout->bscale;
        status=drms_segment_write(segout, outarr, 0);
        free(outarr);

          segout = drms_segment_lookupnum(outrec, 1);
//        outarr = drms_array_create(usetype, 2, length, v2hbt, &status);
        outarr = drms_array_create(usetype, 2, length, bt_out, &status);
        outarr->bzero=segout->bzero;
        outarr->bscale=segout->bscale;
        status=drms_segment_write(segout, outarr, 0);
        free(outarr);

          segout = drms_segment_lookupnum(outrec, 2);
//        outarr = drms_array_create(usetype, 2, length, v2hbp, &status);
        outarr = drms_array_create(usetype, 2, length, bp_out, &status);
        outarr->bzero=segout->bzero;
        outarr->bscale=segout->bscale;
        status=drms_segment_write(segout, outarr, 0);
        free(outarr);

        if (status != DRMS_SUCCESS)
        {
          fprintf(stderr, "ERROR: problem writing output segment: status = %d, T_REC = %s, input recnum = %lld, histrecnum = %lld\n", status, trecstr, inrec->recnum, histrecnum);
          return 0;
        }

//        drms_copykey(outrec, inrec, "T_REC");
        drms_setkey_int(outrec, "QUALITY", quality);
        drms_setkey_int(outrec, "MAPMMAX", (mapped_lmax - 0.5) * rescale + 0.5);
        drms_setkey_int(outrec, "SINBDIVS", (sinb_divisions-0.5) * rescale + 0.5);
        drms_setkey_float(outrec, "CRPIX1", mapcols_out/2.0+ 0.5);
        drms_setkey_float(outrec, "CRVAL1", mapl0);
        drms_setkey_float(outrec, "CROTA1", 0.0);
        drms_setkey_float(outrec, "CDELT1", longinterval/rescale);
        drms_setkey_float(outrec, "CRPIX2", maprows_out/2.0 + 0.5);
        drms_setkey_float(outrec, "CRVAL2", 0.0);
        drms_setkey_float(outrec, "CROTA2", 0.0);
        drms_setkey_float(outrec, "CDELT2", sinBdelta/rescale);
        drms_setkey_string(outrec, "CTYPE1", "CRLN_CEA");
        drms_setkey_string(outrec, "CTYPE2", "CRLT_CEA");
        drms_setkey_string(outrec, "CUNIT1", "deg");
        drms_setkey_string(outrec, "CUNIT2", "sinlat");

// set keywords for magnetic pipeline
        drms_setkey_float(outrec, "MAPRMAX", rmax);
        drms_setkey_float(outrec, "MAPBMAX", bmax);
        drms_setkey_float(outrec, "MAPLGMAX", longmax_adjusted);
        drms_setkey_float(outrec, "MAPLGMIN", longmin_adjusted);
        drms_setkey_int(outrec, "INTERPO", interpolation);
        drms_setkey_int(outrec, "LGSHIFT", longitude_shift);
        drms_setkey_int(outrec, "MCORLEV", mag_correction);
        drms_setkey_int(outrec, "MOFFSET", mag_offset);
        drms_setkey_int(outrec, "CARSTRCH", carrStretch);
        drms_setkey_float(outrec, "DIFROT_A", diffrotA);
        drms_setkey_float(outrec, "DIFROT_B", diffrotB);
        drms_setkey_float(outrec, "DIFROT_C", diffrotC);

        dsignout=vsign*drms_getkey_double(inrec, "DATASIGN", &status);
        if (status != DRMS_SUCCESS) 
          dsignout=vsign;
        dsignout/=fabs(dsignout);
        drms_setkey_int(outrec, "DATASIGN", (int)dsignout);

        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(outrec, "DATE", tnow);

//        drms_close_record(outrec, DRMS_INSERT_RECORD);
      }

      if (verbflag > 1) 
      {
        wt=getwalltime();
        ct=getcputime(&ut, &st);
        fprintf(stdout, 
                "    remap done, %.2f ms wall time, %.2f ms cpu time\n", 
                wt-wt2,
                ct-ct2);
      }

      if (!h2mflag && !tsflag)
        goto skip;

      inp=v2hbr;
      outp=h2mptr;

      mean = 0.0;
      if (subtract_mean)  /* get mean of entire remapped image */
      {
        nmean = 0;
        ip = inp;
        for (row=0; row<maprows; row++)
        {
          for (col=0; col<mapcols; col++) 
          {
            if (!isnan(*ip)) 
            {
              nmean++; 
              mean += *ip;
            }
            ip++;
          }
        }
        mean /= (nmean ? nmean : 1);
      }

      for (row=0; row<maprows; row++) 
      {

        if (cent_long == 0) 
        {

/* Old code with 0 at beginning of array */

          nok = 0; 
          bp = buf; 
          wp = weight;
          ip = inp + mapcols * row;
          for (col=0; col<mapcols; col++) 
          {
            if (!isnan(*ip)) 
            {
              *bp++ = *wp * (*ip - mean);
               nok += 1;
            }
            else 
              *bp++ = 0.0;
            ip++;
            wp++;
          }

          /* zero fill */
          for (i=col; i<nfft; i++)
            *bp++ = 0.0;
        }
        else 
        {

/* New code with 0 at center meridian */
/* Assumes that input array is symmetric around meridian */

          nok = 0;

/* First copy right side of meridian */

          bp = buf;
          ip = inp + mapcols * row+mapcols2;
          wp = weight + mapcols2;
          if (subtract_mean) 
          {
            for (col=0; col<=mapcols2; col++) 
            {
              if (!isnan(*ip)) 
              {
                *bp++ = *wp * (*ip - mean);
                nok += 1;
              }
              else 
                *bp++ = 0.0;
              ip++;
              wp++;
            }
          }
          else 
          {
            for (col=0; col<=mapcols2; col++) 
            {
              if (!isnan(*ip)) 
              {
                *bp++ = *wp * *ip;
                nok += 1;
              }
              else 
                *bp++ = 0.0;
              ip++;
              wp++;
            }
          }

/* Then do left side of meridian */
          bp = buf+lfft-mapcols2;
          ip = inp + mapcols * row;
          wp = weight;
          if (subtract_mean) 
          {
            for (col=0; col<mapcols2; col++) 
            {
              if (!isnan(*ip)) 
              {
                *bp++ = *wp * (*ip - mean); 
                nok += 1;
              }
              else 
                *bp++ = 0.0;
              ip++;
              wp++;
            }
          }
          else 
          {
            for (col=0; col<mapcols2; col++) 
            {
              if (!isnan(*ip)) 
              {
                *bp++ = *wp * *ip;
                nok += 1;
              }
              else 
                *bp++ = 0.0;
              ip++;
              wp++;
            }
          }

/* Finally zero fill */
          bp = buf+mapcols2+1;
          for (i=0; i<lfft-mapcols; i++)
            *bp++ = 0.0;

        } /* End of copying if statement */

        if ((zero_miss == 0) && (nok != mapcols)) 
        {

/* Stuff with missing */
          for (i=0; i<nout; i++) 
          {  
            op = outp + i * maprows + row;
            *op = DRMS_MISSING_FLOAT;
          }
        }
        else 
        {
          if (normalize) 
            norm = 1./sqrt(2*(double)nfft/nok)/lfft;
          else
            norm = 1./lfft;

            /* Fourier transform */
          fftwf_execute(fftwp);

            /* transpose, normalize, this is where most time is spent */
            /* First do real part */
          for (i=0; i<nout/2; i++) 
          {  
            op = outp + 2*i * maprows + row;
            *op = wbuf[i]*norm;
          }
            /* Do imaginary part */
            /* Imaginary part of m=0 is 0*/
          *(outp+row+maprows)=0.0;
            /* Use normx to get the complex conjugate */
          normx=-norm;
          for (i=1; i<nout/2; i++) 
          {  
            op = outp + 2*i * maprows + row + maprows;
            *op = wbuf[lfft-i]*normx;
          }

        }
      } /* Next row */

      if (h2mflag)
      {
//        outrec = drms_create_record(drms_env, h2mout, lifetime, &status);
        outrec=h2mrecset->records[ih2mrec++];
        drms_copykeys(outrec, inrec, 0, kDRMS_KeyClass_Explicit);
        DRMS_Link_t *histlink = hcon_lookup_lower(&outrec->links, histlinkname);
        DRMS_Link_t *srclink = hcon_lookup_lower(&outrec->links, srclinkname);
        if (histlink != NULL)
          drms_setlink_static(outrec, histlinkname,  histrecnum);
        if (srclink != NULL)
          drms_setlink_static(outrec, srclinkname,  inrec->recnum);
        if (segoutflag)
          segout = drms_segment_lookup(outrec, segnameout);
        else
          segout = drms_segment_lookupnum(outrec, 0);
        length[0]=maprows;
        length[1]=nout;
        outarr = drms_array_create(usetype, 2, length, h2mptr, &status);
        drms_setkey_int(outrec, "TOTVALS", maprows*nout);
        set_statistics(segout, outarr, 1);
        outarr->bzero=segout->bzero;
        outarr->bscale=segout->bscale;
        status=drms_segment_write(segout, outarr, 0);
        free(outarr);

        if (status != DRMS_SUCCESS)
        {
          fprintf(stderr, "ERROR: problem writing output segment: status = %d, T_REC = %s, input recnum = %lld, histrecnum = %lld\n", status, trecstr, inrec->recnum, histrecnum);
          return 0;
        }

//        drms_copykey(outrec, inrec, "T_REC");
        drms_setkey_int(outrec, "QUALITY", quality);
        drms_setkey_int(outrec, "MAPMMAX", mapped_lmax);
        drms_setkey_int(outrec, "SINBDIVS", sinb_divisions);
        drms_setkey_int(outrec, "LMAX", lmax);
        drms_setkey_double(outrec, "CRPIX1", maprows/2.0 + 0.5);
        drms_setkey_double(outrec, "CRVAL1", 0.0);
        drms_setkey_double(outrec, "CROTA1", 0.0);
        drms_setkey_double(outrec, "CDELT1", sinBdelta);
        drms_setkey_double(outrec, "CRPIX2", 1.0);
        drms_setkey_double(outrec, "CRVAL2", 0.0);
        drms_setkey_double(outrec, "CROTA2", 0.0);
        drms_setkey_double(outrec, "CDELT2", 1.0);
        drms_setkey_string(outrec, "CTYPE1", "CRLT_CEA");
        drms_setkey_string(outrec, "CTYPE2", "CRLN_FFT");
        drms_setkey_string(outrec, "CUNIT1", "rad");
        drms_setkey_string(outrec, "CUNIT2", "m");

        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(outrec, "DATE", tnow);
//        drms_close_record(outrec, DRMS_INSERT_RECORD);
      }

      if (verbflag > 1) 
      {
        wt2=getwalltime();
        ct2=getcputime(&ut2, &st2);
        fprintf(stdout, 
                "    fft and transpose done, %.2f ms wall time, %.2f ms cpu time\n", 
                wt2-wt,
                ct2-ct);
      }

      if (!tsflag)
        goto skip;

      inptr=h2mptr;
      imageoffset = imagesize * irec; 
      for (mx = 0; mx < 2*msize; mx++) /* for each m, re and im */
      {
        moffset = mx * maprows;
        mptr = inptr + moffset;
        for (latx = 0; latx < nlat; latx++)
        {
          poslatx = nlat + latx; neglatx = nlat - 1 - latx;
          evenpart[latx] = mptr[poslatx] + mptr[neglatx];
          oddpart[latx] = mptr[poslatx] - mptr[neglatx];   
        }
        workptr = folded + imageoffset + moffset;
        scopy_ (&nlat, evenpart, &increment1, workptr, &increment1);
        workptr += nlat;
        scopy_ (&nlat, oddpart, &increment1, workptr, &increment1);
      }

      skip:
      inrec = drms_recordset_fetchnext(drms_env, inrecset, &fetchstat, &chunkstat, NULL);

      if (inrec != NULL) 
      {
        trec = drms_getkey_time(inrec, "T_REC", &status);
        PARAMETER_ERROR("T_REC")
        trecind=(trec-tstart+cadence/2)/cadence;
        sprint_time(trecstr, trec, "TAI", 0);
      }

    } /* end loop on each input image for this timechunk */

    if (verbflag)
      printf("  number of bad images = %d\n", bc);

    if (!tsflag)
      goto continue_outer_loop;

//needed if recordset does not extend to end of a timechunk
    while (irec < nrecs)
    {
      bad[bc++]=irec;
      irec++;
    }

    if (bc == nrecs)
    {
      nodata=1;
    }
    else
    {
      while (bc > 0)
      {
        imageoffset=imagesize*bad[--bc];
        for (i=0;i<imagesize;i++) folded[imageoffset+i]=0.0;
      }
    }

    if (verbflag)
    {
      wt2=getwalltime();
      ct2=getcputime(&ut2, &st2);
      fprintf(stdout, 
                "  images processed, %.2f ms wall time, %.2f ms cpu time\n", 
                wt2-wt1,
                ct2-ct1);
    }

// skip to here if no input records in a given timechunk
   skip_norecs:  

   /* we now have folded data for a chunk of sn's */
   /* now do Jesper's tricks */

   /* ldone is the last l for which plm's have been set up */
    ldone=-1;

   /* loop on each chunk of l's */
    for (lchunk = lchunkfirst; lchunk <= lchunklast; lchunk++)
    {
      lfirst = lchunk * lchunksize; 
      llast = lfirst + lchunksize - 1;
      lfirst = maxval(lfirst,lmin);
      llast = minval(llast,lmax);
      /* get the first and last indexes into the l-m array */
      ifirst = lfirst*(lfirst+1)/2;
      ilast = llast*(llast+1)/2+llast;

      if (verbflag > 1) 
      {
        wt3=getwalltime();
        ct3=getcputime(&ut3, &st3);
        printf("  processing lchunk %d, lmin = %d, lmax = %d\n", lchunk, lfirst, llast);
      }

//      outrec = drms_create_record(drms_env, outseries, lifetime, &status);
      outrec=outrecset->records[irecout++];
/*
      if (status != DRMS_SUCCESS || outrec == NULL)
      {
        fprintf(stderr,"ERROR: unable to open record in output dataseries %s, status = %d, histrecnum = %lld\n", outseries, status, histrecnum);
        return 0;
      }
*/
      if (histlink != NULL)
        drms_setlink_static(outrec, histlinkname,  histrecnum);

      if (nodata)
        goto skip_nodata;

      /* now the size of the output array is known */
      length[0]=2*nrecs;           /* accomodate re & im parts for each sn */
      length[1]=(ilast-ifirst+1);  /* for each l & m, lfirst <= l <= llast */

      outarr = drms_array_create(usetype, 2, length, NULL, &status);

      if (segoutflag)
        segout = drms_segment_lookup(outrec, segnameout);
      else
        segout = drms_segment_lookupnum(outrec, 0);

      if (segout == NULL || outarr == NULL || status != DRMS_SUCCESS)
      {
        fprintf(stderr,"ERROR: problem with output segment or data array: lfirst = %d, llast = %d, length = [%d, %d], status = %d, iset = %d, T_START= %s, histrecnum = %lld", 
                                                                          lfirst, llast, length[0], length[1], status, iset, tstartstr, histrecnum);
        return 0; 
      }

      outptr = (float *)(outarr->data);

      /* loop on each m */
      for (m = 0; m <= llast; m++)
      {

        modd = is_odd(m);
        meven = !modd;
        lstart = maxval(lfirst,m); /* no l can be smaller than this m */

         /* set up masks (plm's) for this m and chunk in l */

        if ((lstart-1) == ldone)
        {
           /* get saved plms if any */
          if ((lstart - 2) >= m)
            for (latx = 0; latx < nlat; latx++)
              plm[(lstart-2)*nlat + latx] = saveplm[m*nlat + latx];
          if ((lstart - 1) >= m)
            for (latx = 0; latx < nlat; latx++)
              plm[(lstart-1)*nlat + latx] = saveplm[msize*nlat + m*nlat + latx];
  
           /* then set up the current chunk */
//          setplm_ (&lstart, &llast, &m, &nlat, indx, latgrid, &nlat, plm); 
          remap_setplm2(lstart, llast, m, nlat, indx, latgrid, nlat, plm, NULL); 

        }
        else
        {
         /* This fixes the lmin != 0 problem */
//          setplm_ (&m, &llast, &m, &nlat, indx, latgrid, &nlat, plm);
          remap_setplm2(m, llast, m, nlat, indx, latgrid, nlat, plm, NULL);
        }

         /* save plm's for next chunk */
        if ((llast-1) >= m)
          for (latx = 0; latx < nlat; latx++)
            saveplm[m*nlat + latx] = plm[(llast - 1)*nlat + latx];
        for (latx = 0; latx < nlat; latx++)
          saveplm[msize*nlat + m*nlat + latx] = plm[llast*nlat + latx];
        ldone=llast;

         /* copy plm's into masks */
         /* note that this converts from double to single precision */
         /* the test prevents underflows which gobble CPU time */
/* Hmmm... looks like if statement is not needed. Weird...
         for (l = lstart; l <= llast; l++) {
            moffset = l * nlat;
            for (latx = 0; latx < nlat; latx++) {
               plmptr = plm + moffset + latx;
               if (is_very_small (*plmptr))
                  masks [(l-lstart)*nlat + latx] = 0.0;
               else
                  masks [(l-lstart)*nlat + latx] = *plmptr;
            }
            plmptr = plm + moffset;
            maskptr = masks+(l-lstart)*nlat;
            dscopy_(&nlat,plmptr,maskptr);
         } 
*/
        plmptr=plm+nlat*lstart;
        maskptr=masks;
        latx=nlat*(llast-lstart+1);
//        dscopy_(&latx,plmptr,maskptr);
        int ilatx;
        for (ilatx=0;ilatx<latx;ilatx++)
          maskptr[ilatx]=plmptr[ilatx];

         /* for each sn in snchunk */
//         for (sn = infsn; sn <= inlsn; sn++) {
        for (snx=0; snx<nrecs; snx++)
        {
            /* select folded data for real/imag l's and this m 
               into temporay arrays for matrix multiply */
//            snx = sn - infsn;
            /* TO DO - pull offset calculations out of loop */
/* New code with odd leading dimension */
          scopy_ (&nlat, 
                  folded + nlat*(4*m+modd) + snx*imagesize,
                  &increment1, 
                  real4evenl + snx*nlatx,
                  &increment1);
          scopy_ (&nlat, 
                  folded + nlat*(4*m+meven) + snx*imagesize, 
                  &increment1,
                  real4oddl + snx*nlatx,
                  &increment1);
          scopy_ (&nlat, 
                  folded + nlat*(4*m+2+modd) + snx*imagesize, 
                  &increment1,
                  imag4evenl + snx*nlatx,
                  &increment1);
          scopy_ (&nlat, 
                  folded + nlat*(4*m+2+meven) + snx*imagesize, 
                  &increment1,
                  imag4oddl + snx*nlatx,
                  &increment1);
        } /* end loop through snchunk */ 


         /* do even l's */
        lfirsteven = is_even(lstart) ? lstart : lstart+1;
        llasteven = is_even(llast) ? llast : llast-1;
        nevenl = (llasteven-lfirsteven)/2 + 1; /* number of even l's */
         /* do real part */
         /* All parts used to have alpha=&one, now have alpha=&cnorm */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &nevenl,   /* number of even l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,    /* scalar multiplier of op(A) */
                real4evenl,  /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat*(lfirsteven-lstart), /* matrix B */
                &maprows,  /* 2*nlat, use every other row (nlat long) of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn*2*(lfirsteven-lstart), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */
         /* do imag part */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &nevenl,   /* number of even l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,    /* scalar multiplier of op(A) */
                imag4evenl,  /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat*(lfirsteven-lstart), /* matrix B */
                &maprows,  /* 2*nlat, use every other nlat-long row of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn*(2*(lfirsteven-lstart)+1), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */

         /* do odd l's */
        lfirstodd = is_odd(lstart) ? lstart : lstart+1;
        llastodd = is_odd(llast) ? llast : llast-1; 
        noddl = (llastodd-lfirstodd)/2 + 1; /* number of odd l's */
         /* do real part */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &noddl,    /* number of odd l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm ,   /* scalar multiplier of op(A) */
                real4oddl,   /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat*(lfirstodd-lstart), /* matrix B */
                &maprows,  /* 2*nlat, use every other nlat-long row of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn*2*(lfirstodd-lstart), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */
         /* do imag part */
        sgemm_ (transpose, /* form of op(A) */ 
                normal,    /* form of op(B) */ 
                &nsn,      /* number of sn's */
                &noddl,    /* number of odd l's for this m */
                &nlat,     /* number of latitudes */
                &cnorm,    /* scalar multiplier of op(A) */
                imag4oddl,  /* matrix A */
                &nlatx,     /* use every nlat-long row of A */ 
                masks + nlat*(lfirstodd-lstart), /* matrix B */
                &maprows,  /* 2*nlat, use every other nlat-long row of B */ 
                &zero,     /* scalar multiplier of C */
                outx + nsn*(2*(lfirstodd-lstart)+1), /* matrix C (output) */ 
                &fournsn,  /* use every fourth nsn-long row of C */
                1,         /* length of transpose character string */
                1);        /* length of normal character string */

         /* copy outx into out sds */
         /* alternate real and imaginary values in out - as in pipeLNU */
        for (l = lstart; l <= llast; l++)
        {
          fromoffset = 2*nsn*(l-lstart);
          tooffset = 2*nsn*(l*(l+1)/2 + m -ifirst);
          scopy_ (&nsn,            
                  outx+fromoffset, 
                  &increment1, 
                  outptr+tooffset, 
                  &increment2);
          scopy_ (&nsn, 
                  outx+fromoffset+nsn, 
                  &increment1, 
                  outptr+tooffset+1, 
                  &increment2);
        } /* end loop through l's for this m */


      } /* end loop on m */


      outarr->bzero=segout->bzero;
      outarr->bscale=segout->bscale;
      status=drms_segment_write(segout, outarr, 0);
      drms_free_array(outarr);
      nsegments++;

      if (status != DRMS_SUCCESS)
      {
        fprintf(stderr, "ERROR: problem writing output segment: status = %d, T_START = %s, LMIN = %d, LMAX = %d, histrecnum = %lld\n", status, tstartstr, lfirst, llast, histrecnum);
        return 0;
      }

      skip_nodata:

      drms_setkey_int(outrec, "LMIN", lfirst);
      drms_setkey_int(outrec, "LMAX", llast);
      drms_setkey_time(outrec, "T_START", tstart);
      drms_setkey_time(outrec, "T_STOP", tstart+chunksecs);
      drms_setkey_time(outrec, "T_OBS", tstart+chunksecs/2);
      drms_setkey_time(outrec, "DATE__OBS", tstart);
      drms_setkey_string(outrec, "TAG", tag);
      if (verflag)
         drms_setkey_string(outrec, "VERSION", version);

      for (i=0;i<16;i++)
        setbits(calversout,4*i+3,4,nybblearrout[i]);
      drms_setkey_longlong(outrec, "CALVER64", calversout);

      if (nodata)
        drms_setkey_int(outrec, "QUALITY", QUAL_NODATA);
      else if (mixflag)
        drms_setkey_int(outrec, "QUALITY", QUAL_MIXEDCALVER);
      else
        drms_setkey_int(outrec, "QUALITY", 0);

      // these could be constant, but set them just in case
      drms_setkey_int(outrec, "MAPMMAX", mapped_lmax);
      drms_setkey_int(outrec, "SINBDIVS", sinb_divisions);
      drms_setkey_float(outrec, "T_STEP", cadence);

      ndt=chunksecs/cadence;
      drms_setkey_int(outrec, "NDT", ndt);

      tnow = (double)time(NULL);
      tnow += UNIX_epoch;
      drms_setkey_time(outrec, "DATE", tnow);

//      drms_close_record(outrec, DRMS_INSERT_RECORD);
      nrecords++;

      if (verbflag > 1) 
      {
        wt=getwalltime();
        ct=getcputime(&ut, &st);
        fprintf(stdout, 
                "    %.2f ms wall time, %.2f ms cpu time\n", 
                wt-wt3,
                ct-ct3);
      }


    } /* end loop on each chunk of l's */

    if (verbflag)
    {
      wt1=getwalltime();
      ct1=getcputime(&ut1, &st1);
      fprintf(stdout, "SHT of timechunk %d complete: %.2f ms wall time, %.2f ms cpu time\n", iset, 
              wt1-wt2, ct1-ct2);
    }

//    free(v2hbt); free(v2hbp);
    continue_outer_loop:
    tstart+=chunksecs;
  } /* end loop on each time chunk */


  if (chunkstat != kRecChunking_LastInRS && chunkstat != kRecChunking_NoMoreRecs)
    fprintf(stderr, "WARNING: input records remain after last output record: chunkstat = %d\n", (int)chunkstat);

  drms_close_records(inrecset, DRMS_FREE_RECORD);
  if (tsflag)
    drms_close_records(outrecset, DRMS_INSERT_RECORD);
  if (v2hflag)
    drms_close_records(v2hrecset, DRMS_INSERT_RECORD);
  if (h2mflag)
    drms_close_records(h2mrecset, DRMS_INSERT_RECORD);


  wt=getwalltime();
  ct=getcputime(&ut, &st);
  if (verbflag && tsflag) 
  {
    printf("number of records created  = %d\n", nrecords);
    printf("number of segments created = %d\n", nsegments);
  }
  if (verbflag) 
  {
    fprintf(stdout, "total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
            wt-wt0, ct-ct0);
  }

  if (!error)
    printf("module %s successful completion\n", cmdparams.argv[0]);
  else
    printf("module %s failed to produce %d timechunks: histrecnum = %lld\n", cmdparams.argv[0], error, histrecnum);

  return 0;
}

void do_boxcar(float *image_in, float *image_out, int in_nx, int in_ny, float fscale)
{
  int iscale, nvector, vec_half;
  int inx, iny, outx, outy, i, j;
  float val;

      iscale = 1.0/fscale + 0.5;
      nvector = iscale;
      vec_half = nvector/2;

      int in_go = (iscale-1)/2.0 + 0.5;
      int out_nx = in_nx * fscale + 0.5;
      int out_ny = in_ny * fscale + 0.5;

        for (outy = 0; outy < out_ny; outy += 1)
          for (outx = 0; outx < out_nx; outx += 1)
            {
            double total = 0.0;
            double weight = 0.0;
            int nn = 0;
            for (j = 0; j < nvector; j += 1)
              {
              iny = outy*iscale + in_go + j - vec_half;
              for (i = 0; i < nvector; i += 1)
                {
                inx = outx*iscale + in_go + i - vec_half;
                if (inx >= 0 && inx < in_nx && iny >=0 && iny < in_ny)
                  {
                  val = image_in[in_nx*(iny) + inx];
                  if (!drms_ismissing_float(val))
                    {
                    double w = 1.0;
                    total += w*val;
                    weight += w;
                    nn++;
                    }
                  }
                }
              }
            image_out[out_nx*outy + outx] = (nn > 0 ? total/weight : DRMS_MISSING_FLOAT);
            }

}

