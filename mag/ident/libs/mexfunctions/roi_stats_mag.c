#include "mex.h"  /* must appear first */
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "mexhead.h" /* my mex defines */
#include "Doc/roi_stats_mag_docstring.h"  /* autogenerated */


/**************************************************************

%roi_stats_mag: accumulate statistics on regions
% 
% [s,names]=roi_stats_mag(x,y,mag,geom,mode)
% * A set of per-region statistics is gathered based on regions 
% encoded in the image inputs x (containing region tags, 1..Nr, or
% 0 for no tag) and y (containing region indicators, 0/1).  The
% statistics are functions of region configuration, and of 
% line-of-sight magnetic field mag.
% * x must be Nan, or a nonnegative integer.  Pixels in the 
% range 1..Nr, where Nr is the number of regions, are treated as
% tags for being within the named region.
% Inputs x of 0 or NaN are treated as not belonging to any
% labeled region.
% * y must be NaN, 0, or 1 -- 1 indicates activity present.  Not all
% tagged pixels will be active; typically the tagged pixels are
% large blobs, and the active pixels are finer details within each
% blob.
% * For each region, a row of statistics is computed:
%    1: rgnnum = # pixels tagged (all may not be active)
%    2: rgnsize = projected (flat) area in microhemispheres (0..1e6)
%    3: rgnarea = un-projected (solid-angle) area in microhemispheres (0..1e6)
%    4: arnum = # active pixels (x == 1)
%    5: arsize = projected (flat) active area in microhemispheres (0..1e6)
%    6: ararea = unprojected (solid-angle) active area in microhemis (0..1e6)
%    7: arlat,lon = average location, weighted by un-projected area
%    9: arminlat,lon  lower corner of (lat,lon) bounding box
%   11: armaxlat,lon  upper corner of (lat,lon) bounding box
%   13: rgnbtot = sum of absolute LoS flux within the identified region
%   14: rgnbnet = net LoS flux within the identified region
%   15: rgnbpos = absolute value of total positive LoS flux
%   16: rgnbneg = absolute value of total negative LoS flux
%   17: arfwtlat,lon = flux-weighted center of active pixels
%   19: arfwtpos_lat,lon = flux-weighted center of positive flux
%   21: arfwtneg_lat,lon = flux-weighted center of negative flux
%   23: daysgone = #days until bounding box vanishes from front of disk
%   24: daysback = #days until bounding box first reappears on front
% These definitions are in a machine-readable include file,
% roi_stats_mag_def.h
% * The solar disk is at location given by center, observed at tip 
% angle beta, in degrees.
% * Optionally returns the standard text name of each of the
% computed statistics. 
% * This is implemented as a MEX file.
% 
% Inputs:
%   int  x(m,n);
%   int  y(m,n);
%   real mag(m,n);
%   real center(5);  -- [x0 y0 r_sun b p]
%   string mode;
% 
% Outputs:
%   real s(nr,ns)
%   opt string names(ns)
% 
% See Also:

% turmon june 2010, sep 2010

****************************************************************/


/* standard boilerplate */
static const char *progname = "roi_stats_mag";
#define PROGNAME roi_stats_mag
#define SHORTNAME rsm

#define NARGIN_MIN	5	   /* min number of inputs */
#define NARGIN_MAX	5	   /* max number of inputs */
#define NARGOUT_MIN	0	   /* min number of outputs */
#define NARGOUT_MAX	2	   /* max number of outputs */

#define ARG_X     0
#define ARG_Y     1
#define ARG_MAG   2
#define ARG_GEOM  3
#define ARG_MODE  4

#define ARG_S     0
#define ARG_NAMES 1

static const char *in_specs[NARGIN_MAX] = {
  "RM",
  "RM",
  "RM",
  "RV(5)",
  "SV"};
static const char *in_names[NARGIN_MAX] = {
  "x",
  "y",
  "mag",
  "geom",
  "mode"};
static const char *out_names[NARGOUT_MAX] = {
  "s",
  "names"};

#include "roi_stats_mag_defs.h"

/*
 * Solar rotation constants
 * all in deg/day
 * first has sidereal rate, less Earth's rotational rate, to make synodic rate
 * if you're sensitive to these values, you're probably doing something wrong...
 */

#define ROTATE_CONST_K0 (14.6-0.9865) 
#define ROTATE_CONST_K1 (-2.2)
#define ROTATE_CONST_K2 0.0


/*************************************************************
 *
 * HELPER ROUTINES
 *
 *************************************************************/


/*
 * extracts various (OK, 1) mode parameters from a free-form string
 *
 * returns 0 if the string was OK, an explanatory tag otherwise
 */
static
char *
extract_mode(char *mode, int *Msesw) 
{
  char *word;
  char *sep = ", \t";

  if (!mode) return "null"; /* should not happen */
  for (word = strtok(mode, sep); word; word = strtok(NULL, sep)) {
    /* sesw, or transposed */
    if      (strcasecmp(word, "sesw"  ) == 0) *Msesw = 1;
    else if (strcasecmp(word, "sene"  ) == 0) *Msesw = 0;
    /* unrecognized word: not OK */
    else return word; 
  }
  return 0; /* OK */
}

/* simple utilities */ 
#define WithinROI(y)   ((!isnan(y)) && ((y) > 0)) /* not NaN and not 0 */
#define ActiveLabel(y) ((!isnan(y)) && ((y) > 0)) /* not NaN and not 0 */

/*
 * stat_compute: calculate statistics for each region
 */
static
void
stat_compute(double *roi,   // input roi image
	     double *ar,    // input acr image
	     double *mag,   // input mgram image, or NULL
	     double **s,    // result: statistics per region (2d)
	     int Nr,        // number of regions
	     int maxx,      // image size 
	     int maxy,
	     int strx,      // stride along x, units of doubles
	     int stry,      // stride along y
	     // geom
	     double x0,     // x center, in C coords, origin at 0 
	     double y0,     // y center, in C coords, origin at 0 
	     double R,      // radius, in pixels 
	     double B,      // tilt angle, radians 
	     double P)      // twist angle, radians

{
  int x, y;               // loop counters
  int r;                  // region counter
  int offset;             // offset into arrays
  double lat, lon;        // hold current lat,lon
  double p1x, p1y, p1z;   // point P1 coords in 3D
  double p2x, p2y, p2z;   // point P2 coords in 3D
  double temp;            // workspace slot
  double area1;           // un-projected area
  double tau;             // a threshold for finding the limb
  double mag1, magM1;     // current mgram value
  const double nan = mxt_getnand(); // cache nan
  const double cosB = cos(B);
  const double sinB = sin(B);
  const double cosP = cos(P);
  const double sinP = sin(P);
  const double Rinv = 1/R;

  /* 
   * STEP 0: set initial values for accumulators 
   */
  // most stuff starts at zero, so do it all now
  bzero(&(s[0][0]), ((size_t)RS_num_stats)*Nr*sizeof(s[0][0]));
  for (r = 0; r < Nr; r++) {
    // 0: rgnnum, rgnsize, rgnarea
    // 0: arnum,  arsize,  ararea     
    // 0: arlat,  arlon
    // set greater than +2pi: arminlat, arminlon
    s[RS_arminlat][r] = s[RS_arminlon][r] =  10;
    // set lower than -2pi: armaxlat, armaxlon
    s[RS_armaxlat][r] = s[RS_armaxlon][r] = -10;
    // 0: rgnbtot, rgnbnet, rgnbpos, rgnbneg,    
    // 0: arfwtlat,    arfwtlon
    // 0: arfwtposlat, arfwtposlon
    // 0: arfwtneglat, arfwtneglon
    // -1: daysgone, daysback ("don't know")
    s[RS_daysgone][r] = s[RS_daysback][r] = -1;
  }

  /* 
   * STEP 1: loop over all pixels and accumulate region info 
   */
  for (x = 0; x < maxx; x++) 
    for (y = 0; y < maxy; y++) {
      offset = x*strx + y*stry;  // still in units of doubles
      r = roi[offset];
      if (!WithinROI(r)) continue;
      r--; /* make region number (1..Nr) into a region index (0..Nr-1) */
      mag1  = mag ? mag[offset] : nan; // abbreviate the magnetic field...
      magM1 = fabs(mag1); // ...and its magnitude
      /*
       * Find the coordinates, including lat/lon
       */
      /* find original vector coordinates P2 = (p2x, p2y, p2z) */
      p2x = x - x0;
      p2y = y - y0;
      temp = (R + p2x) * (R - p2x) - p2y*p2y; /* R*R-x*x = (R-x)*(R+x) */
      /* temp is the residual which goes into the p2z component */
      if (temp <= 0) continue; /* off-disk (should not happen) */
      p2z = sqrt(temp);   /* totally OK -- on visible side */
      area1 = R/p2z;      /* note: the original, unrotated, 1/z */
      if (area1 > 1e4) area1 = 1e4; // cap it -- could do better here
      /* rotate by p in the xy plane: P1 = rot(P2) (the p-angle) */
      p1x =  cosP * p2x + sinP * p2y;
      p1y = -sinP * p2x + cosP * p2y;
      p1z =  p2z;
      /* rotate by +beta in yz plane (x axis fixed): P2 = rot(P1) */
      p2x =  p1x;  /* no change to x */
      p2y =  cosB * p1y + sinB * p1z;
      p2z = -sinB * p1y + cosB * p1z;
      /* find (lat,lon) of the point */
      lon = atan2(p2x, p2z);  /* p2=(1,0,0), at limb, has lon=pi/2 */
      lat = asin(p2y * Rinv);

      /*
       * Update all the accumulators 
       * (listed in order)
       */
      // whole-region accumulators
      s[RS_rgnnum ][r]++;               // count pixels
      s[RS_rgnsize][r]++;               // count projected area
      s[RS_rgnarea][r] += area1;        // unprojected area, proportional to 1/z
      // ar-only accumulators: add up when this ROI pixel is also active
      if (ActiveLabel(ar[offset])) {
	s[RS_arnum ][r]++;              // count pixels
	s[RS_arsize][r]++;              // count projected area
	s[RS_ararea][r] += area1;       // unprojected area, proportional to 1/z
      }
      // area-weighted lat/lon
      s[RS_arlat][r] += lat * area1;
      s[RS_arlon][r] += lon * area1;
      // lat/lon bounding box
      if (lat < s[RS_arminlat][r]) s[RS_arminlat][r] = lat;
      if (lat > s[RS_armaxlat][r]) s[RS_armaxlat][r] = lat;
      if (lon < s[RS_arminlon][r]) s[RS_arminlon][r] = lon;
      if (lon > s[RS_armaxlon][r]) s[RS_armaxlon][r] = lon;
      // line-of-sight flux
      s[RS_rgnbtot][r] += magM1;
      s[RS_rgnbnet][r] += mag1;
      s[RS_rgnbpos][r] += (mag1 > 0) ? magM1 : 0.0;
      s[RS_rgnbneg][r] += (mag1 < 0) ? magM1 : 0.0;
      // flux-weighted lat/lon
      s[RS_arfwtlat][r] += lat * magM1;
      s[RS_arfwtlon][r] += lon * magM1;
      // flux-weighted centers of + and - flux
      s[RS_arfwtposlat][r] += (mag1 > 0) ? lat*magM1 : 0.0;
      s[RS_arfwtposlon][r] += (mag1 > 0) ? lon*magM1 : 0.0;
      s[RS_arfwtneglat][r] += (mag1 < 0) ? lat*magM1 : 0.0;
      s[RS_arfwtneglon][r] += (mag1 < 0) ? lon*magM1 : 0.0;
    } /* for (x,y) */

  /* 
   * STEP 2: normalize and do unit conversions on s array 
   */
  // 2a: convert accumulated sums into averages, in order of definition
  for (r = 0; r < Nr; r++) {
    // average location, weighted by un-projected area
    s[RS_arlat][r] /= s[RS_ararea][r];  // the sum of the weights 1/z
    s[RS_arlon][r] /= s[RS_ararea][r];
    // flux-weighted center of active pixels
    s[RS_arfwtlat   ][r] /= s[RS_rgnbtot][r];
    s[RS_arfwtlon   ][r] /= s[RS_rgnbtot][r];
    // flux-weighted center of positive flux
    s[RS_arfwtposlat][r] /= s[RS_rgnbpos][r];
    s[RS_arfwtposlon][r] /= s[RS_rgnbpos][r];
    // flux-weighted center of negative flux
    s[RS_arfwtneglat][r] /= s[RS_rgnbneg][r];
    s[RS_arfwtneglon][r] /= s[RS_rgnbneg][r];
  }
  // 2b: compute (daysgone, daysback) from lat-lon info.  
  /* Note, the visible disk is characterized by: 
   *    cos(lon) >= tan(beta) * tan(lat).  
   * Thus, find the values dt for which: 
   *    lon + (dt)*(vel) 
   * satisfies this. The above condition is independent of P.
   */
  for (r = 0; r < Nr; r++) {
    // find threshold for longitude (using average lat, still in radians)
    tau = (sinB/cosB) * tan(s[RS_arlat][r]);
    if (fabs(tau) >= 1) continue;  // beta, lat: all rotations visible!
    tau = acos(tau);     /* in the range 0..pi */
    // find angular velocity at our latitude (still in radians)
    temp = sin(s[RS_arlat][r]);      // sin of average latitude
    temp *= temp;                    // temp = sin(lat)^2
    temp =    ROTATE_CONST_K0 + 
      temp * (ROTATE_CONST_K1 + 
      temp * (ROTATE_CONST_K2));     // k0 + k1*temp + k2*temp^2, in deg/day
    temp *= (M_PI/180.0);            // convert deg/day -> radians/day
    // translate temp, an angular rate in radians/day, into a time in days
    // minlon: used for last passage behind limb
    // maxlon: used for first reappearance
    s[RS_daysgone][r] = (tau - s[RS_arminlon][r])/temp;  
    s[RS_daysback][r] = ((2*M_PI - tau) - s[RS_armaxlon][r])/temp; 
  }
  // 2c: unit conversions, in order of definition
  for (r = 0; r < Nr; r++) {
    const double rad2deg = 180.0/M_PI;
    const double size2uhemi = 1e6/(M_PI*R*R);   // 1 hemi has pi R^2 pixels
    const double area2uhemi = 1e6/(2*M_PI*R*R); // 1 hemi integrates to 2piR^2
    
    // area/sizes to microhemispheres
    // no conversion: RS_rgnnum, RS_arnum
    s[RS_rgnsize][r] *= size2uhemi;
    s[RS_rgnarea][r] *= area2uhemi;
    s[RS_arsize ][r] *= size2uhemi;
    s[RS_ararea ][r] *= area2uhemi;
    // angles to degrees
    s[RS_arlat   ][r] *= rad2deg;  s[RS_arlon   ][r] *= rad2deg;
    s[RS_arminlat][r] *= rad2deg;  s[RS_arminlon][r] *= rad2deg;
    s[RS_armaxlat][r] *= rad2deg;  s[RS_armaxlon][r] *= rad2deg;
    // no conversion: rgnbtot, rgnbnet, rgnbpos, rgnbneg
    // angles to degrees
    s[RS_arfwtlat   ][r] *= rad2deg;  s[RS_arfwtlon   ][r] *= rad2deg;
    s[RS_arfwtposlat][r] *= rad2deg;  s[RS_arfwtposlon][r] *= rad2deg;
    s[RS_arfwtneglat][r] *= rad2deg;  s[RS_arfwtneglon][r] *= rad2deg;
  }
  // 2d: plug in NaNs for empty regions
  for (r = 0; r < Nr; r++) {
    if (s[RS_arnum][r] == 0) {
      s[RS_arlat   ][r] = s[RS_arlon   ][r] = nan;
      s[RS_arminlat][r] = s[RS_arminlon][r] = nan;
      s[RS_armaxlat][r] = s[RS_armaxlon][r] = nan;
    }
  }
}
   
/*
 * Are classes all integers, 0, or NaN?
 */
static
int
count_classes(double *y,
	      int N)
{
  double Nr = 0; /* class count */
  int Nri;
  int n; 

  for (n = 0; n < N; n++)
    if (!isnan(y[n])) {
      /* below commented out for speed */
      if (y[n] < 0 /* || floor(y[n]) != y[n] */ )
	return -1;
      if (y[n] > Nr)
	Nr = y[n];
    }
  Nri = (int) Nr;
  if (Nri != Nr)
    return -1; // trouble
  else
    return Nri;
}


/*************************************************************
 *
 * DRIVER ROUTINE
 *
 *************************************************************/

/*
 * Gateway routine
 */
#ifdef StaticP  /* undefined under mex */
StaticP
#endif
void 
mexFunction(
	    int nlhs, 
	    mxArray *plhs[], 
	    int nrhs,
	    const mxArray *prhs[])
{
  int m,  n;                 // size of x
  int maxx, maxy;            // sizes, transposed if needed
  int strx, stry;            // strides in images along x and y dims
  int Nr;                    // number of regions found
  int mode_sesw;             // sesw = 1, transposed = 0
  char *mode, *word;         // label for converted mode string and substring
  double *g;                 // convenience
  double x0, y0, r, b, p;    // geom
  double *mag;
  double **stat2;
  char errstr[256];

  /* Hook for introspection (function signature, docstring) */
  if (nrhs < 0) { 
    plhs[0] = mxt_PackSignature((mxt_Signature) (-nrhs), 
				NARGIN_MIN, NARGIN_MAX, 
				NARGOUT_MIN, NARGOUT_MAX, 
				in_names, in_specs, out_names, docstring);
    return;
  }
  /*
   * check args
   */
  if ((nrhs < NARGIN_MIN) || (nrhs > NARGIN_MAX))
    mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			   "%s: Expect %d <= input args <= %d",
			   progname, NARGIN_MIN, NARGIN_MAX), errstr));
  if ((nlhs < NARGOUT_MIN) || (nlhs > NARGOUT_MAX))
    mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			   "%s: Expect %d <= output args <= %d",
			   progname, NARGOUT_MIN, NARGOUT_MAX), errstr));
  mexargparse(nrhs, prhs, in_names, in_specs, NULL, progname);

  /* get size of input image x, needed below */
  m = (int) mxGetM(prhs[ARG_X]);
  n = (int) mxGetN(prhs[ARG_X]);

  /******************************************************************
   * Read and check arguments in order supplied
   ******************************************************************/

  /*
   * image geometry
   */
  // source image (x) geometry
  g   = mxGetPr(prhs[ARG_GEOM]); // abbreviate geom
  x0 = g[0] - 1;         // x0, in C coordinates
  y0 = g[1] - 1;         // y0, in C coordinates
  r  = g[2];             // disk radius
  b  = g[3]*M_PI/180.0;  // deg -> rad
  p  = g[4]*M_PI/180.0;  // deg -> rad

  /* 
   * `mode' input
   * (required)
   */
  if ((mode = mxArrayToString(prhs[ARG_MODE])) == NULL)
    mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			   "%s: bad mode (non-string?). "
			   "Could not convert mode arg to string.", 
			   progname), errstr));
  /* no default values -- user must specify both explicitly */
  mode_sesw = -1; 
  if ((word = extract_mode(mode, &mode_sesw)) != NULL)
    mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			   "%s: bad word <%s> in mode <%.80s>", 
			   progname, word, /* need a fresh copy of the mode */
			   mxArrayToString(prhs[ARG_MODE])), errstr));
  if (mode_sesw < 0)
    mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			   "%s: need sesw/sene in mode <%.80s>", 
			   progname, /* need a fresh copy of the mode */
			   mxArrayToString(prhs[ARG_MODE])), errstr));
  mxFree(mode); // done with the mode string

  /*
   * create space for output
   */
  /* 1: find how big it must be */
  m  = mxGetM(prhs[ARG_X]);
  n  = mxGetN(prhs[ARG_X]);

  /* 2: make the space */
  /* first, check x and find #regions */
  if ((Nr = count_classes(mxGetPr(prhs[ARG_X]), m*n)) < 0)
    mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			   "%s: Regions x must be NaN or integers >= 0", 
			   progname), errstr));
  /* now, create space for the region info */
  plhs[ARG_S] = mxCreateDoubleMatrix(Nr, (mwSize) RS_num_stats, mxREAL); 

  /*
   * Depending on pixel-order mode, set up strides and bounding boxes in
   * (x,y) coordinates from their (m,n) settings.
   * Note that the (x0,y0) positions, in the geom1 and geom2 inputs, 
   * correspond to x and to y, so they do not need to be swapped.
   */
  if (mode_sesw) {
    // x (E-W axis) varies fastest -- x corresponds to matlab m coordinate
    strx = 1; stry = m; // x varies fastest, x has m points...
    maxx = m; maxy = n; // y has n points
  } else {
    // y (N-S axis) varies fastest -- x corresponds to matlab n coordinate
    strx = m; stry = 1; // y varies fastest, y has m points...
    maxx = n; maxy = m; // and x has n points
  }

  // simplify indexing into stat
  stat2 = mxt_make_matrix2(plhs[ARG_S],   -1, -1, 0.0);
  // account for empty mag case
  mag = (mxGetNumberOfElements(prhs[ARG_MAG]) > 0) ? 
    mxGetPr(prhs[ARG_MAG]) : NULL;

  /*
   * do the computation
   */
  stat_compute(mxGetPr(prhs[ARG_X]), 
	       mxGetPr(prhs[ARG_Y]), 
	       mag,
	       stat2,
	       Nr,                           // computed number of regions
	       maxx, maxy,                   // size
	       strx, stry,                   // strides in units of doubles
	       x0, y0, r, b, p);             // geom
  // for debugging 
  if (0) {
    int r, sn;
    double *stats = mxGetPr(plhs[ARG_S]);

    for (sn = 0; sn < RS_num_stats; sn++)
      printf("\t%s", RS_index2name[sn]);
    printf("\n");
    for (r = 0; r < Nr; r++) {
      printf("P%d: ", r);
      for (sn = 0; sn < RS_num_stats; sn++)
	printf("\t%.3g", stats[sn*Nr+r]);
      printf("\n");
    }
  }
  // free up this indexing tool
  mxFree(stat2);

  /*
   * plug in names if wanted
   */
  if (nlhs > ARG_NAMES) {
    plhs[ARG_NAMES] = mxCreateCharMatrixFromStrings((int) RS_num_stats, 
						    RS_index2name);
    if (plhs[ARG_NAMES] == NULL)
      mexErrMsgTxt((snprintf(errstr, sizeof(errstr),
			     "%s: Failed to create `%s' output for %d names",
			     progname, out_names[ARG_NAMES], 
			     (int) RS_num_stats), 
		    errstr));
  }
}


/* Hook for generic tail matter */
#ifdef MEX2C_TAIL_HOOK
#include "mex2c_tail.h"
#endif

