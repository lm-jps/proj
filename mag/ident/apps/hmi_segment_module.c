/*
 * Module name: hmi_segment_module.c
 * This jsoc module generates masks from LOS B and intensity
 *
 * Current (04/2011) usage:
 * hmi_segment_module 
 *      xm='hmi.M_720s[2010.07.03_12:48:00_TAI/1h]' 
 *      xp='hmi.Ic_720s[2010.07.03_12:48:00_TAI/1h]' 
 *      alpha='0,-4' beta=0.4 y=hmi.Marmask_720s
 *
 *
 * Old (6/2010) usage:
 * hmi_segment_module 
 *   "xm=su_couvidat.HMISeriesLev15245[2010.03.25_00:53:15_TAI]" 
 *   "xp=su_couvidat.HMISeriesLev15545[2010.03.25_00:53:15_TAI]" 
 *   beta=0.4 alpha="[0,-4]" T="[1,1,0.9,0]" y=su_turmon.armask2
 *
 *
 * Michael Turmon, JPL
 *  adapted from code by Yang Liu (10/2009)
 *  updated with additions by Xudong Sun (3/2010)
 *  updated 4/2011
 *
 */

#include "jsoc_main.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

// mex and jsoc inter-operability
#include "mex.h"
#include "mexhead.h"
#include "mxt_jsoc.c"
// Argument placement for the library routine we call
#include "hmi_segment.h"


// exit with error
//  (standard trick to swallow the semicolon)
//  (ensure nonzero exit even if zero status)
#define DIE(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: FATAL: %s, status=%d\n", module_name, msg, status); \
        return(status ? status : 1); \
        } while (0)
// report non-fatal error
//  (standard trick to swallow the semicolon)
//  TBD: should this affect the error status?
#define WARN(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: ERROR: %s.  Continuing.\n", module_name, msg); \
        } while (0)

// some standard JSOC stuff
#define	DTOR	(M_PI / 180.)
#define RAD2ARCSEC	(648000. / M_PI)
// Macros for WCS transformations.  
// ASSUMES: 
//   crpix1, crpix2 = CRPIX1, CRPIX2
//   sina,cosa = sin and cos of CROTA2 resp.
//   crvalx and crvaly are CRVAL1 and CRVAL2, 
//   cdelt = CDELT1 == CDELT2
// PIX_X and PIX_Y are CCD pixel addresses, 
// WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)


/* ################## Timing from T. Larson ################## */

static double getwalltime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000.0 + tv.tv_usec/1000.0;
}

static double getcputime(double *utime, double *stime)
{

  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  *utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec / 1000.0;
  *stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec / 1000.0;
  return *utime + *stime;
}


/********************************************************
 *
 * Boilerplate for input argument processing
 *
 *********************************************************/


/* 
 * command-line parameter names
 * (following the order in the driver arg list)
 */
#define kRecSetXM       "xm"
#define kRecSetXP       "xp"
#define kParIter        "iter"
#define kParT           "T"
#define kParBeta        "beta"
#define kParAlpha       "alpha"
// ctr is determined from KW's
#define kParRho         "rho"
#define kParModel1      "m1"
#define kParModel2      "m2"
// restrict to 2 models
#define kSeriesOut      "y"

char *module_name = "active region mask";
char *mex2c_progname = "hmi_segment_module";

// valid for Yang's test images
static char seg_model_qs_yang[] =
  "[0.204635504622,1.20531725063,1.00777710322,19.346252856,0.00419579448962,0.0564904666445,"
  "0.204635504622,-1.20531725063,1.00777710322,19.346252856,0.00419579448962,-0.0564904666445,"
  "0.232398997138,0,1.02955618645,2.7404836868,0.00427555589387,0,"
  "0.145253155252,0,0.992142854282,247.248235524,0.00504149434045,0,"
  "0.139157664971,0,0.950374859774,74.2478459635,0.000555807820391,0,"
  "0.0617507192555,0,0.986125619938,2854.35303874,0.00151332362565,0,"
  "0.01216845414,0,1.01070360771,89583.9686521,0.0145007628344,0]";
static char seg_model_ar_yang[] =
  "[0.192105993124,-846.854676681,0.743104620169,118603.827293,0.00499678354147,10.6891118764,"
  "0.192105993124,846.854676681,0.743104620169,118603.827293,0.00499678354147,-10.6891118764,"
  "0.126980271901,311.120029981,0.773004292296,18546.8591144,0.00277679162354,0.270699599935,"
  "0.126980271901,-311.120029981,0.773004292296,18546.8591144,0.00277679162354,-0.270699599935,"
  "0.111840866217,1725.55076542,0.432208758685,148513.00487,0.0244679802074,-43.2940104149,"
  "0.111840866217,-1725.55076542,0.432208758685,148513.00487,0.0244679802074,43.2940104149,"
  "0.0690728687575,-2565.81843681,0.13790051683,124931.169245,0.0014162131028,10.5502256214,"
  "0.0690728687575,2565.81843681,0.13790051683,124931.169245,0.0014162131028,-10.5502256214]";

// valid for Sebastien's images taken late march 2010
static char seg_model_qs_sebastien[] = 
  "[0.676700379049,0,13321.3804543,137.665818814,133580569.057,0,"
  "0.30516524961,0,13833.1964163,214.193316404,127827406.447,0,"
  "0.0181343713409,0,13300.8328273,4715.87050328,147742870.796,0]";
static char seg_model_ar_sebastien[] = 
  "[0.373249775824,-534.166827748,13470.5783991,87955.4275489,140934500.036,-21833.1675272,"
  "0.373249775824,534.166827748,13470.5783991,87955.4275489,140934500.036,21833.1675272,"
  "0.118699362858,194.068579955,13621.338073,6919.60349716,150211728.679,-33384.7350668,"
  "0.118699362858,-194.068579955,13621.338073,6919.60349716,150211728.679,33384.7350668,"
  "0.0161017226367,0,13038.102521,165305.443351,73124251.7542,0]";

ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kRecSetXM, "",           "Input mgram data records."},
   {ARG_STRING,  kRecSetXP, "",           "Input pgram data records."},
   {ARG_FLOATS,  kParIter,  "[0,0.05]",   "Iterations (RV(2))"},
   {ARG_FLOATS,  kParT,     "[1,1,0.9,0]","Temperature (RV(4))"},
   {ARG_FLOAT,   kParBeta,  "0.4",        "Smoothness"},
   {ARG_FLOATS,  kParAlpha, "[0.0,0.0]",  "Class Bias"},
   {ARG_FLOAT,   kParRho,   "100",        "Smoothness Anisotropy"},
   // assuming two classes
   // note: should these models be in DRMS instead?
   {ARG_FLOATS,  kParModel1, seg_model_qs_sebastien, "Model for Quiet Sun"},
   {ARG_FLOATS,  kParModel2, seg_model_ar_sebastien, "Model for ARs"},
   {ARG_STRING,  kSeriesOut, "",          "Output data series."},
   {ARG_INT,     "VERB",     "1",         "Verbosity: 0=errors/warnings; 1=all"},
   {ARG_END}
};


int DoIt(void)
{
  // status and input arguments
  int status = DRMS_SUCCESS;
  int verbflag;
  CmdParams_t *params = &cmdparams;
  const char *xmRecQuery = params_get_str(params, kRecSetXM);
  const char *xpRecQuery = params_get_str(params, kRecSetXP);
  const char * yRecQuery = params_get_str(params, kSeriesOut);
  // DRMS data
  DRMS_RecordSet_t *xmRD,  *xpRD,  *yRD;
  DRMS_Record_t    *xmRec, *xpRec, *yRec;
  DRMS_Segment_t   *xmSeg, *xpSeg, *ySeg;
  DRMS_Array_t   *xmArray, *xpArray, *yArray;
  DRMS_Link_t  *srcLink_M, *srcLink_P;
  // mexFunction arguments
  // (these are defined in hmi_segment.h)
  const int nrhs = MXT_Hseg_ARG_m1 + 2;   // two models
  const int nlhs = MXT_Hseg_ARG_y + 1;    // just the labeling y
  mxArray *prhs[MXT_Hseg_NARGIN_MAX];     // input args (too long is OK)
  mxArray *plhs[MXT_Hseg_NARGOUT_MAX];    // output args
  double *xm_tmp, *xp_tmp;         // holding places for array data
  // other variables
  int ds, Nds;                     // iteration
  char **keyname;                  // *keyname is a generic KW
  char *msg;
  int M, N, yDims[2];              // image sizes
  double rSun, x0, y0, p0;         // KWs
  TIME t_rec;
  char time_buf[64]; // time string for t_rec
  // wcs info
  float crvalx, crvaly, cdelt, crpix1, crpix2, crota2, sina, cosa;
  double rsun_ref, dsun_obs;
  // Time measuring -- ugly, but functional
  double wt0, wt1, wt;
  double ut0, ut1, ut;
  double st0, st1, st;
  double ct0, ct1, ct;
  // for Runtime1 keyword
  time_t now;

  // timings
  wt0 = getwalltime();
  ct0 = getcputime(&ut0, &st0);

  // verbosity level
  verbflag = params_get_int(&cmdparams, "VERB");

  // set up fixed arguments: iter, T, beta, alpha, rho, mode, m1, m2
  prhs[MXT_Hseg_ARG_iter]  = mxt_param_to_matrix(params, kParIter, 1, -1);
  prhs[MXT_Hseg_ARG_T]     = mxt_param_to_matrix(params, kParT, 1, -1);
  prhs[MXT_Hseg_ARG_beta]  = mxt_param_to_scalar(params, kParBeta);
  prhs[MXT_Hseg_ARG_alpha] = mxt_param_to_matrix(params, kParAlpha, 1, -1);
  prhs[MXT_Hseg_ARG_ctr]   = mxCreateDoubleMatrix(1, 3, mxREAL); // vals from KWs
  prhs[MXT_Hseg_ARG_rho]   = mxt_param_to_scalar(params, kParRho);
  prhs[MXT_Hseg_ARG_m1]    = mxt_param_to_matrix(params, kParModel1, 6, -1);
  prhs[MXT_Hseg_ARG_m1+1]  = mxt_param_to_matrix(params, kParModel2, 6, -1);

  // check args
  if (!prhs[MXT_Hseg_ARG_iter]  || 
      !prhs[MXT_Hseg_ARG_T]     || 
      !prhs[MXT_Hseg_ARG_beta]  || 
      !prhs[MXT_Hseg_ARG_alpha] || 
      !prhs[MXT_Hseg_ARG_ctr]   || 
      !prhs[MXT_Hseg_ARG_rho])
    DIE("failed to set up a miscellaneous parameter");
  // special-case these as they are tricky
  if (!prhs[MXT_Hseg_ARG_m1]) 
    DIE("Failed to set up model1 (quiet).  Length not /6 ?");
  if (!prhs[MXT_Hseg_ARG_m1+1]) 
    DIE("Failed to set up model2 (active).  Length not /6 ?");


  // debug printf's
  if (0) {
    int i, k;
    for (i = 0; i < 2; i++)
      for (k = 0; k < mxGetN(prhs[MXT_Hseg_ARG_m1+i]); k++)
	printf("\tmod(%d)[%d][0:5] = %.3f %g,%g %g,%g,%g\n", i, k,
	       mxGetPr(prhs[MXT_Hseg_ARG_m1+i])[6*k+0],
	       mxGetPr(prhs[MXT_Hseg_ARG_m1+i])[6*k+1],
	       mxGetPr(prhs[MXT_Hseg_ARG_m1+i])[6*k+2],
	       mxGetPr(prhs[MXT_Hseg_ARG_m1+i])[6*k+3],
	       mxGetPr(prhs[MXT_Hseg_ARG_m1+i])[6*k+4],
	       mxGetPr(prhs[MXT_Hseg_ARG_m1+i])[6*k+5]);
  }

  // open mgrams
  xmRD = drms_open_records(drms_env, (char *) xmRecQuery, &status);
  if (status || xmRD->n == 0)
    DIE("No mgram input data found");
  Nds = xmRD->n;

  // open pgrams
  xpRD = drms_open_records(drms_env, (char *) xpRecQuery, &status);
  if (status || xpRD->n == 0)
    DIE("No pgram input data found");
  if (xpRD->n != Nds)
    DIE("pgram series length does not match mgram series length");

  // open labelings
  // NOTES: Disabled on Oct 27 by Xudong Sun
  // Creating records one per time
  /*
  yRD = drms_create_records(drms_env, Nds, (char *) yRecQuery, DRMS_PERMANENT, &status);
  if (status)
    DIE("Output recordset for labelings could not be created");
  */
	
  // loop over a series of (mgram,pgram) pairs
  for (ds = 0; ds < Nds; ds++) {
    printf("begin record %d of %d\n", ds+1, Nds);fflush(stdout);
    // mark the time
    wt1 = getwalltime();
    ct1 = getcputime(&ut1, &st1);

    xmRec = xmRD->records[ds];
    xpRec = xpRD->records[ds];
	  
    // NOTES: Disabled on Oct 27 by Xudong Sun
    // Creating records one per time, moved to below
    // yRec  =  yRD->records[ds];

    t_rec = drms_getkey_time(xmRec, "T_REC", &status);
    if (status || drms_getkey_time(xpRec, "T_REC", &status) != t_rec) {
      WARN("mgram and pgram input record T_REC's don't match");
      continue;
    }
    // summary
    sprint_time(time_buf, t_rec, "TAI", 0);
    printf("\tprocessing images taken at %s\n", time_buf);
    fflush(stdout);

    // Ephemeris courtesy of Phil's macros above
    // No error checking for now
    crvalx = drms_getkey_float(xmRec, "CRVAL1", &status); // disc center in arcsec
    crvaly = drms_getkey_float(xmRec, "CRVAL2", &status);
    cdelt  = drms_getkey_float(xmRec, "CDELT1", &status); // arcsec, assumimg dx=dy
    crpix1 = drms_getkey_float(xmRec, "CRPIX1", &status); // disk center on ccd
    crpix2 = drms_getkey_float(xmRec, "CRPIX2", &status);
    crota2 = drms_getkey_float(xmRec, "CROTA2", &status); // rotation
    if (status) {
      WARN("could not get WCS from mgram");
      continue;
    }
    sina = sin(crota2 * DTOR); 
    cosa = cos(crota2 * DTOR); 
    x0 = PIX_X(0.0,0.0) - 1.0; // zero-based coordinates
    y0 = PIX_Y(0.0,0.0) - 1.0;
    rsun_ref = drms_getkey_double(xmRec, "RSUN_REF", &status);
    if (status) {
      rsun_ref = 6.96e8; // this is OK
      status = 0;
    }
    dsun_obs = drms_getkey_double(xmRec, "DSUN_OBS", &status);
    if (status) {
      WARN("could not get DSUN_OBS from mgram");
      continue;
    }
    rSun = asin(rsun_ref / dsun_obs) * RAD2ARCSEC / cdelt;
    if (verbflag)
      printf("\tX0=%.3f, Y0=%.3f, RSUN=%.3f, P0=%.2f\n", x0, y0, rSun, crota2);

    // mgram data
    xmSeg = drms_segment_lookupnum(xmRec, 0);
    if (!xmSeg) {
      WARN("problem looking up mgram");
      continue;
    }
    // load as doubles
    xmArray = drms_segment_read(xmSeg, DRMS_TYPE_DOUBLE, &status);
    if (status) {
      drms_free_array(xmArray);
      WARN("problem reading mgram");
      continue;
    }
    // input size (same for all images)
    M = xmArray->axis[0]; N = xmArray->axis[1];

    // pgram data
    xpSeg = drms_segment_lookupnum(xpRec, 0);
    if (!xpSeg) {
      drms_free_array(xmArray);
      WARN("problem looking up pgram");
      continue;
    }
    xpArray = drms_segment_read(xpSeg, DRMS_TYPE_DOUBLE, &status);
    if (status) {
      drms_free_array(xmArray);
      drms_free_array(xpArray);
      WARN("problem reading pgram");
    }
    // Check if dimensions match
    if (xpArray->axis[0] != M || xpArray->axis[1] != N) {
      drms_free_array(xmArray);
      drms_free_array(xpArray);
      WARN("pgram and mgram dimensions don't match");
      continue;
    }

    // Moved here by X. Sun Oct 27
    // Input segments exist!
    yRec = drms_create_record(drms_env, (char *) yRecQuery, DRMS_PERMANENT, &status);
    if (status) {
      WARN("Output recordset for labelings could not be created");
      continue;
    }
	
    // set up varying arguments to function: xm, xp, ctr
    if ((prhs[MXT_Hseg_ARG_xm] = mxCreateDoubleMatrix(M, N, mxREAL)) == NULL)
      DIE("failed calloc for mgram"); // make this fatal
    xm_tmp = mxGetPr(prhs[MXT_Hseg_ARG_xm]); // hold data segment
    mxSetPr(prhs[MXT_Hseg_ARG_xm], (double *) xmArray->data);
    if ((prhs[MXT_Hseg_ARG_xp] = mxCreateDoubleMatrix(M, N, mxREAL)) == NULL)
      DIE("failed calloc for pgram"); // make this fatal
    xp_tmp = mxGetPr(prhs[MXT_Hseg_ARG_xp]); // hold data segment
    mxSetPr(prhs[MXT_Hseg_ARG_xp], (double *) xpArray->data);
    // this was allocated outside the loop
    // mexfunction coordinates are one-based (origin at (1,1)), so add 1.
    // mexfunction ctr[1] direction varies fastest, ctr[0] direction
    // varies slowest, so flip around.
    // (makemrfdiscwts and clean_edge_label are affected by x/y issues)
    mxGetPr(prhs[MXT_Hseg_ARG_ctr])[0] = y0+1;  // FIXME: x, y OK?
    mxGetPr(prhs[MXT_Hseg_ARG_ctr])[1] = x0+1;
    mxGetPr(prhs[MXT_Hseg_ARG_ctr])[2] = rSun;
 
    // print some sample values
    if (verbflag)
      printf("\tM[%d,%d] = %f\n\tP[%d,%d] = %f\n", 
	     M/2, N/2, ((double *) xmArray->data)[M*N/2+M/2],
	     M/2, N/2, ((double *) xpArray->data)[M*N/2+M/2]);

    // call the function
    printf("segment_module: hmi_segment calling mexfunction.\n");fflush(stdout);
    msg = main_hmi_segment(nlhs, plhs, nrhs, prhs);
    if (msg == NULL) {
      printf("segment_module: hmi_segment returned OK from mexfunction.\n");
    } else {
      /* mexErrMsgTxt within mexFunction() will route thru this clause */
      printf("segment_module: returned via error from mexfunction.\n");
      /* fail */
      status = 1;
      DIE("mexFunction failed (mexErrMsgTxt)");
    }

    // write labeling as the first segment
    // (it is a double here, but defined in the .jsd with BITPIX=8)
    printf("segment_module: looking up segment\n");
    ySeg = drms_segment_lookupnum(yRec, 0);
    printf("segment_module: array_create\n");
    printf("segment_module: y[%d,%d] = %g\n", 
	   M/2, N/2, mxGetPr(plhs[MXT_Hseg_ARG_y])[M*N/2+M/2]);
    fflush(stdout);

    // convert the double mxArray to a char drms_array
    // (this is an independent copy of the output array)
    yArray = mxArray2drms_array(plhs[MXT_Hseg_ARG_y], DRMS_TYPE_CHAR);
    if (!yArray)
      DIE("problem creating drms array for labeling");

    // Link array to segment
    yArray->parent_segment = ySeg;
	  
    // some stats
    // updated by xudong oct 13 2010
    int missval = 0, totalval = yArray->axis[0] * yArray->axis[1];
    char *yData = (char *)yArray->data;
    for (int ii = 0; ii < totalval; ii++) {
      if (isnan(yData[ii])) missval++;
    }

    printf("segment_module: segment_write\n");fflush(stdout);
    status = drms_segment_write(ySeg, yArray, 0);
    if (status)
      DIE("problem writing drms segment for labeling");

    // free the changing part of the argument list
    // (prep: replace data segments in xm, xp)
    mxSetPr(prhs[MXT_Hseg_ARG_xm], (double *) xm_tmp);
    mxSetPr(prhs[MXT_Hseg_ARG_xp], (double *) xp_tmp);
    // now free them
    printf("segment_module: destroy array: xm, xp, y\n");fflush(stdout);
    mxDestroyArray(prhs[MXT_Hseg_ARG_xm]);
    mxDestroyArray(prhs[MXT_Hseg_ARG_xp]);
    mxDestroyArray(plhs[MXT_Hseg_ARG_y ]);

    // free input DRMS arrays
    printf("segment_module: within-loop frees\n");fflush(stdout);
    drms_free_array(xmArray); 
    drms_free_array(xpArray);
    drms_free_array(yArray);		// added by Xudong Apr 13 2011

    // Establish links to source data
    if ((srcLink_M = hcon_lookup_lower(&yRec->links, "MDATA")) != NULL)
      drms_setlink_static(yRec, "MDATA", xmRec->recnum);
    if ((srcLink_P = hcon_lookup_lower(&yRec->links, "PDATA")) != NULL)
      drms_setlink_static(yRec, "PDATA", xpRec->recnum);
            
    // Essential prime key: T_REC
    drms_copykey(yRec, xmRec, "T_REC");
	  
    // date and build version
    // updated by xudong oct 13 2010
    drms_setkey_string(yRec, "BLD_VERS", jsoc_version);
    drms_setkey_time(yRec, "DATE", CURRENT_SYSTEM_TIME);
	  
    // stats
    // updated by xudong oct 13 2010
    drms_setkey_int(yRec, "TOTVALS", totalval);
    drms_setkey_int(yRec, "DATAVALS", totalval - missval);
    drms_setkey_int(yRec, "MISSVALS", missval);

    // lets me see when the record was actually updated
    time(&now);
    drms_setkey_string(yRec, "RUNTIME1", asctime(localtime(&now)));
	  
    // Oct 27
    drms_close_record(yRec, DRMS_INSERT_RECORD);

    // Timing info
    if (verbflag) {
      wt = getwalltime();
      ct = getcputime(&ut, &st);
      printf("record %d time used: %.3f s wall, %.3f s cpu\n", 
	     ds+1, (wt - wt1)*1e-3, (ct - ct1)*1e-3);
    }
    printf("record %d done\n", ds+1); fflush(stdout);
  } 
  printf("segment_module: close record (mgram)\n");
  drms_close_records(xmRD, DRMS_FREE_RECORD);
  printf("segment_module: close record (pgram)\n");
  // note: when pgram == mgram, this segfaults
  drms_close_records(xpRD, DRMS_FREE_RECORD);
	
  // Moved into the loop Oct 27 X. Sun
  // Insert one by one
  /*
  printf("segment_module: close record (y)\n");
  drms_close_records(yRD, DRMS_INSERT_RECORD);
  */
	
  // free the unchanging arguments
  printf("segment_module: final destroy arrays\n");
  mxDestroyArray(prhs[MXT_Hseg_ARG_iter]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_T]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_beta]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_alpha]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_ctr]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_rho]);

  if (verbflag) {
    wt = getwalltime();
    ct = getcputime(&ut, &st);
    printf("total time used: %.3f s wall, %.3f s cpu\n", 
	   (wt - wt0)*1e-3, (ct - ct0)*1e-3);
  }

  // exit OK
  return DRMS_SUCCESS;
}
