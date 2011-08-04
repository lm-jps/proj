/*
 * Module name: hmi_segment_module.c
 * This jsoc module generates masks from LOS field and intensity
 *
 * We currently use hmi.M_720s and hmi.Ic_noLimbDark_720s, but this is
 * change-able by specifying a different model.  In particular,
 * to use hmi.M_720s and hmi.Ic_720s, you can use:
 *   model=/builtin/hmi.M_Ic_720s.production
 * To use hmi.Ic_noLimbDark_720s instead, use:
 *   model=/builtin/hmi.M_Ic_noLimbDark_720s.production
 * To see all builtin models, use model=/builtin/all
 *
 * You must specify M and Ic series that match up in time (T_REC) 
 * and image coordinates.  (The code checks for this.)  By coordinates,
 * we mean that XO, Y0, R_SUN, CROTA2, and CRLT_OBS, or their WCS
 * analogs, must agree.
 *
 * Accepts a verbose option (VERB=v where v is 0, 1, 2, or 3).
 *
 * Current (07/2011) usage:
 *
 * hmi_segment_module 
 *      xm='hmi.M_720s[2010.07.03_12:48:00_TAI/1h]' 
 *      xp='hmi.Ic_noLimbDark_720s[2010.07.03_12:48:00_TAI/1h]' 
 *      model=/builtin/hmi.M_Ic_noLimbDark_720s.production beta=0.4 y=hmi.Marmask_720s
 *
 * Michael Turmon, JPL
 *  adapted from code by Yang Liu (10/2009)
 *  updated with additions by Xudong Sun (3/2010)
 *  updated 4/2011
 *  updated 6/2011
 *
 */

#include "jsoc_main.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

// mex and jsoc inter-operability
#include "mex.h"
#include "mexhead.h"
#include "mxt_jsoc.c"
// Argument placement for the library routine we call
#include "hmi_segment.h"
// Layout of patch-statistics vector
#include "roi_stats_mag_defs.h"
// Classification model settings
#include "segment_modelset.h"
// Keyword propagation
#include "propagate_keys.c"


// exit with error
//  (standard trick to swallow the semicolon)
//  (ensure nonzero exit even if zero status)
#define DIE(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: FATAL: %s. (status=%d)\n", module_name, msg, status); \
        return(status ? status : 1); \
        } while (0)
// report non-fatal error
//  (standard trick to swallow the semicolon)
#define WARN(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: WARNING: %s. Continuing.\n", module_name, msg); \
        } while (0)
// facilitate verbose output
// if flag is true, print the message in the form:
// <first><module_name>: <message>"
// Usage:
// V_printf(VERB > 0, "\t", "Mask(%d) = %d\n", 2048, mask[2048]);
void
V_printf(int flag, char *first, char *format, ...) {
  va_list args;
  extern char *module_name;

  va_start(args, format);
  if (flag) {
    // first is a string, even "" -- print the module name too
    // otherwise, omit it
    if (first)
      printf("%s%s: ", first, module_name);
    vprintf(format, args);
    fflush(stdout);
  }
  va_end(args);
}


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


/************************************************************* 
 *
 * Timing
 * from T. Larson
 *
 *************************************************************
 */
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

/************************************************************* 
 *
 * Keyword maintenance
 *
 *************************************************************
 */

/*
 * printkey_vector
 * helper function: sprintf's a vector mxArray into a string
 * str of length n chars.  Uses snprintf to avoid overflow.
 */
static void
printkey_vector(char *str, int n, const mxArray *a)
{
  const int numel = mxGetNumberOfElements(a);
  const double *x = mxGetPr(a);
  int np = 0; // number of chars printed this time
  int i;

  for (i = 0; (i < numel) && (n > 0); i++, n -= np) {
    np = snprintf(str, n, "%.4g%s", x[i], (i+1 < numel) ? "," : "");
    str += np;
  }
}

static int
setkey_code_info(DRMS_Record_t *outRec)
{
  extern char *hmi_segment_version; // defined below
  char keystr[80]; // plausible FITS keyword length
  int ok;
  int not_ok = 0;

  // version number (of the present code)
  snprintf(keystr, sizeof(keystr), "%s ver: %s", module_name, hmi_segment_version);
  ok = drms_setkey_string(outRec, "ARMCODEV", keystr);
  if (ok != DRMS_SUCCESS) not_ok++;
  // url points to descriptive text in wiki
  ok = drms_setkey_string(outRec, "ARMDOCU", "http://jsoc.stanford.edu/jsocwiki/ARmaskDataSeries");
  if (ok != DRMS_SUCCESS) not_ok++;

  return not_ok;
}

static int
setkey_mask_qual(DRMS_Record_t *outRec, 
		 int nclean,
		 char *y, 
		 int M, int N)
{
  const int totalval = M*N;
  int missval, i;
  int quality;
  int ok;
  int not_ok = 0;

  // number of in-valid mask values (nan in the double version)
  for (missval = i = 0; i < totalval; i++) 
    if (y[i] == 0) 
      missval++;
  // pixel counts, updated by xudong oct 13 2010
  ok = drms_setkey_int(outRec, "TOTVALS",  totalval);
  if (ok != DRMS_SUCCESS) not_ok++;
  ok = drms_setkey_int(outRec, "DATAVALS", totalval - missval);
  if (ok != DRMS_SUCCESS) not_ok++;
  ok = drms_setkey_int(outRec, "MISSVALS", missval);
  if (ok != DRMS_SUCCESS) not_ok++;
  // nclean
  ok = drms_setkey_int(outRec, "ARM_NCLN", nclean);
  if (ok != DRMS_SUCCESS) not_ok++;
  // quality (more would be better...)
  quality = 0x0;
  quality |= (nclean > 1000) ? 1 : 0;
  ok = drms_setkey_int(outRec, "ARM_QUAL", quality);
  if (ok != DRMS_SUCCESS) not_ok++;

  return not_ok;
}

static int
setkey_mask_stats(DRMS_Record_t *outRec, 
		  double *s, 
		  seg_modelset_t *modelset)
{
  const int Nc = modelset->nclass;
  int ok;
  int not_ok = 0;
  seg_onemodel_t *m1;
  int class;
  int sn;

  // find mask value for `active' class 
  for (class = 0; class < Nc; class++) {
    m1 = &(modelset->models[class]);
    if (strstr(m1->name, "active") != 0)
      break;
  }
  // insert keys if there was a match for the active class
  // printf("mask_stats: class = %d among 0..%d\n", class, Nc-1);
  if (class < modelset->nclass) {
    // insert active-region summary keywords
    for (sn = 0; sn < RS_num_stats; sn++) {
      const char *name1 = RS_index2mask_keyname[sn];
      double stat1 = s[sn*Nc+class];
      // insert if HMI keyword name is non-null
      if (name1) {
	switch (*name1) {
	case 'i':
	  ok = drms_setkey_int(outRec, name1+1, (int) stat1);
	  break;
	case 'f':
	  ok = drms_setkey_float(outRec, name1+1, (float) stat1);
	  // printf("\tP %d:\t%s -> %.3g\n", p, name1+1, stat1);
	  break;
	default:
	  WARN("setkey_mask_stats: unknown key type (first char of HMI keyname)");
	  ok = -1;
	  break;
	} // end switch
        if (ok != DRMS_SUCCESS) not_ok++;
      }
    } // end for sn
  } else {
    not_ok = -1; // could not find `active' class
  }

  return not_ok;
}

static int
setkey_model_params(DRMS_Record_t *outRec, 
		    seg_modelset_t *modelset)
{
  int ok;
  int not_ok = 0;
  seg_onemodel_t *m1;
  int inx;
  int modelnum;
  char *key_match;
  // format: (string to find in modelset, keyword name)* (NULL, NULL)
  // currently just: quiet, active
  static char *model2key[] = {"quiet", "QUIET", "active", "ACTIVE", NULL, NULL};

  // obtain #classes from modelset
  ok = drms_setkey_int(outRec, "NCLASS", modelset->nclass);
  if (ok != DRMS_SUCCESS) not_ok++;
  // offdisk is always 0
  ok = drms_setkey_int(outRec, "OFFDISK", 0);
  if (ok != DRMS_SUCCESS) not_ok++;
  // on-patch mask is always 64
  ok = drms_setkey_int(outRec, "ON_PATCH", 64);
  if (ok != DRMS_SUCCESS) not_ok++;
  // add class numbers
  for (inx = 0; model2key[inx] != NULL; inx += 2) {
    // return value is 0..nclass-1, or -1
    modelnum = seg_modelname2modelnum(modelset, model2key[inx]);
    if (modelnum < 0) {
      not_ok++; // did not find the model
    } else {
      // add 1 because off-disk is coded as 0
      ok = drms_setkey_int(outRec, model2key[inx+1], modelnum+1);
      if (ok != DRMS_SUCCESS) not_ok++;
    }
  }

  return not_ok;
}


static int
setkey_algorithm_params(DRMS_Record_t *outRec, 
			const char *modelName, 
			mxArray *prhs[])
{
  char kstr[256];
  const int sz = sizeof(kstr);
  int ok;
  int not_ok = 0;

  // iter
  printkey_vector(kstr, sz, prhs[MXT_Hseg_ARG_iter]);
  ok = drms_setkey_string(outRec, "ARM_ITER", kstr);
  if (ok != DRMS_SUCCESS) not_ok++;
  // temperature
  printkey_vector(kstr, sz, prhs[MXT_Hseg_ARG_T]);
  ok = drms_setkey_string(outRec, "ARM_TEMP", kstr);
  if (ok != DRMS_SUCCESS) not_ok++;
  // beta (treat as string, just in case)
  printkey_vector(kstr, sz, prhs[MXT_Hseg_ARG_beta]);
  ok = drms_setkey_string(outRec, "ARM_BETA", kstr);
  if (ok != DRMS_SUCCESS) not_ok++;
  // alpha
  printkey_vector(kstr, sz, prhs[MXT_Hseg_ARG_alpha]);
  ok = drms_setkey_string(outRec, "ARM_ALFA", kstr);
  if (ok != DRMS_SUCCESS) not_ok++;
  // rho
  ok = drms_setkey_float(outRec, "ARM_RHO", 
			 (float) mxGetScalar(prhs[MXT_Hseg_ARG_rho]));
  if (ok != DRMS_SUCCESS) not_ok++;
  // patch smoothing kernel width
  // FIXME: need to have this param
  ok = drms_setkey_float(outRec, "ARM_KPAR", -1.0);
  if (ok != DRMS_SUCCESS) not_ok++;
  // FIXME: need to have this param
  // snprintf(kstr, sz, "%f,%f", 1.2345, 2.71828);
  snprintf(kstr, sz, "(unused)");
  ok = drms_setkey_string(outRec, "ARM_KWT", kstr);
  if (ok != DRMS_SUCCESS) not_ok++;
  // smoothing threshold
  // FIXME: need to have this param
  ok = drms_setkey_float(outRec, "ARM_THRS", -1.0);
  if (ok != DRMS_SUCCESS) not_ok++;
  // edge
  ok = drms_setkey_float(outRec, "ARM_EDGE", 
			 mxGetScalar(prhs[MXT_Hseg_ARG_edge]));
  if (ok != DRMS_SUCCESS) not_ok++;
  // model name
  drms_setkey_string(outRec, "ARM_MODL", modelName);
  if (ok != DRMS_SUCCESS) not_ok++;

  return not_ok;
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
#define kParEdge        "edge"
#define kParIter        "iter"
#define kParT           "T"
#define kParBeta        "beta"
// ctr is determined from KW's
#define kParRho         "rho"
#define kParModel       "model"
#define kSeriesOut      "y"
#define kVerb           "VERB"

char *module_name = "hmi_segment_module";
char *mex2c_progname = "hmi_segment_module";
char *hmi_segment_version = "2.0";  // used for keywords


ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kRecSetXM,  "",           "Input mgram data records."},
   {ARG_STRING,  kRecSetXP,  "",           "Input pgram data records."},
   {ARG_FLOAT,   kParEdge,   "2.5",        "Limb width set-to-quiet (pixels)"},
   {ARG_FLOATS,  kParIter,   "[0,0.05]",   "Iterations (RV(2))"},
   {ARG_FLOATS,  kParT,      "[1,1,0.9,0]","Temperature (RV(4))"},
   {ARG_FLOAT,   kParBeta,   "0.4",        "Smoothness"},
   {ARG_FLOAT,   kParRho,    "100",        "Smoothness Anisotropy"},
   {ARG_STRING,  kParModel,  "",           "Model set name (`/builtin/all' for list)"},
   {ARG_STRING,  kSeriesOut, "",           "Output data series."},
   {ARG_INT,     kVerb,      "1",          "Verbosity: 0=errs only; 1, 2, 3 = more"},
   {ARG_END}
};


/********************************************************
 *
 * Main routine: read arguments, loop over input images
 *
 *********************************************************/

int DoIt(void)
{
  // status and input arguments
  int status = DRMS_SUCCESS;
  char errmsg[256];              // buffer for error messages
  char *msg;                     // return code/message
  CmdParams_t *params = &cmdparams;
  const char *xmRecQuery = params_get_str(params, kRecSetXM);
  const char *xpRecQuery = params_get_str(params, kRecSetXP);
  const char  *yRecQuery = params_get_str(params, kSeriesOut);
  const double  edge_pix = params_get_float(params, kParEdge);
  const char  *modelName = params_get_str(params, kParModel);
  const int     verbflag = params_get_int(params, kVerb); // an int, not a flag
  // DRMS data
  DRMS_RecordSet_t *xmRD,  *xpRD,  *yRD;
  DRMS_Record_t    *xmRec, *xpRec, *yRec;
  DRMS_Segment_t   *xmSeg, *xpSeg, *ySeg;
  DRMS_Array_t   *xmArray, *xpArray, *yArray;
  DRMS_Link_t  *srcLink_M, *srcLink_P;
  // mexFunction arguments
  // (these are defined in hmi_segment.h)
  int nrhs;                         // size of modelset is known at runtime
  const int nlhs = MXT_Hseg_NARGOUT_MAX;  // ask for all args
  mxArray *prhs[MXT_Hseg_NARGIN_MAX];     // input args (too long is OK)
  mxArray *plhs[MXT_Hseg_NARGOUT_MAX];    // output args
  double *xm_tmp, *xp_tmp;         // holding places for array data
  // other variables
  seg_modelset_t *modelset;        // see hmi_segment_models.h
  seg_onemodel_t *m1;              // one member of the above, e.g. QS or AR
  int Nclass, class, class_quiet;  // Classes in above model set
  int ds, Nds, ds_good;            // iteration
  int M, N;                        // image sizes
  // keyword-related
  double rSun, x0, y0, p0, beta;
  TIME t_rec;
  char time_buf[64]; // time string for t_rec
  // wcs info
  float crvalx, crvaly, cdelt, crpix1, crpix2, crota2, crlt, sina, cosa;
  double rsun_ref, dsun_obs;
  // Time measuring -- ugly, but functional
  double wt0, wt1, wt;
  double ut0, ut1, ut;
  double st0, st1, st;
  double ct0, ct1, ct;

  // timings
  wt0 = getwalltime();
  ct0 = getcputime(&ut0, &st0);
  // look up model name
  if (strcmp(modelName, SEG_BUILTIN_PREFIX "/all") == 0) {
    modelset = NULL; // indicate trouble, but omit error message
  } else {
    modelset = seg_modelsetname2modelset(modelName, &msg, 1);
    if (modelset == NULL)
      fprintf(stderr, "%s: Did not recognize model <%s>\nError: %s\n", 
	      module_name, modelName, msg);
  }
  // for enlightenment, show known models
  if (modelset == NULL) {
    // explanatory messages before DIE() go to stderr
    fprintf(stderr, "Built-in models are:\n");
    seg_modelsets_fprintf(stderr, (verbflag > 0) ? 1 : 0);
    DIE("Did not find a known model.  Cannot make masks");
  }
  // print the model out
  if (verbflag > 2)
    seg_modelset_fprintf(stdout, 1, modelset);
  class_quiet = seg_modelname2modelnum(modelset, "quiet"); // quiet index in modelset
  if (class_quiet < 0)
    DIE("Quiet class not found in modelset, check modelset.models[*]->name");
  Nclass = modelset->nclass; // number of classes in the modelset
  nrhs = MXT_Hseg_ARG_m1 + Nclass; // 1 class -> pass in nrhs = ARG_m1+1 args, etc.
  if (nrhs > MXT_Hseg_NARGIN_MAX) {
    // explanatory messages before DIE() go to stderr
    fprintf(stderr, "Given modelset has %d classes, max is %d\n", 
	    Nclass, MXT_Hseg_NARGIN_MAX - MXT_Hseg_ARG_m1);
    DIE("Given modelset has too many classes, must recompile");
  }

  // set up fixed arguments: iter, T, beta, ctr, rho, mode
  prhs[MXT_Hseg_ARG_edge]  = mxCreateDoubleMatrix(1, 2, mxREAL); // set vals later
  prhs[MXT_Hseg_ARG_iter]  = mxt_param_to_matrix(params, kParIter, 1, -1);
  prhs[MXT_Hseg_ARG_T]     = mxt_param_to_matrix(params, kParT, 1, -1);
  prhs[MXT_Hseg_ARG_beta]  = mxt_param_to_scalar(params, kParBeta);
  prhs[MXT_Hseg_ARG_geom]  = mxCreateDoubleMatrix(1, 5, mxREAL); // vals from KWs
  prhs[MXT_Hseg_ARG_rho]   = mxt_param_to_scalar(params, kParRho);
  // check args and allocations
  if (!prhs[MXT_Hseg_ARG_edge]  || 
      !prhs[MXT_Hseg_ARG_iter]  || 
      !prhs[MXT_Hseg_ARG_T]     || 
      !prhs[MXT_Hseg_ARG_beta]  || 
      !prhs[MXT_Hseg_ARG_geom]  || 
      !prhs[MXT_Hseg_ARG_rho])
    DIE("failed to set up a miscellaneous parameter");
  // ok to set vals now
  mxGetPr(prhs[MXT_Hseg_ARG_edge])[0] = edge_pix; // pulled from cmd params
  mxGetPr(prhs[MXT_Hseg_ARG_edge])[1] = class_quiet+1; // in 1..Nclass (0 is offdisk)

  // another fixed argument: alpha
  prhs[MXT_Hseg_ARG_alpha] = mxCreateDoubleMatrix(1, Nclass, mxREAL);
  if (!prhs[MXT_Hseg_ARG_alpha])
    DIE("failed to set up parameter alpha");
  memcpy(mxGetPr(prhs[MXT_Hseg_ARG_alpha]), modelset->alpha, 
	 Nclass*sizeof(double));

  // remaining fixed arguments: m_1...m_Nclass
  for (class = 0; class < Nclass; class++) {
    m1 = &(modelset->models[class]); // current model, for convenience
    prhs[MXT_Hseg_ARG_m1 + class] = mxCreateDoubleMatrix(6, m1->k, mxREAL);
    if (!prhs[MXT_Hseg_ARG_m1 + class]) {
      fprintf(stderr, "Fatal: could not set up model %d (%d components, %s)\n", 
	      class, m1->k, m1->name);
      DIE("failed to set up a model parameter matrix");
    }
    memcpy(mxGetPr(prhs[MXT_Hseg_ARG_m1 + class]), 
	   m1->params, 6 * (m1->k) * sizeof(double));
  }

  // debug printf's
  if (0) {
    int k;
    for (class = 0; class < Nclass; class++) {
      double *p1 = mxGetPr(prhs[MXT_Hseg_ARG_m1+class]);
      for (k = 0; k < mxGetN(prhs[MXT_Hseg_ARG_m1+class]); k++)
	V_printf(verbflag > 2, NULL, 
		 "\tmod(%d)[%d][0:5] = %.3f %g,%g %g,%g,%g\n", class, k,
		 p1[6*k], p1[6*k+1], p1[6*k+2], p1[6*k+3], p1[6*k+4], p1[6*k+5]);
    }
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

  // loop over a series of (mgram,pgram) pairs
  for (ds_good = ds = 0; ds < Nds; ds++) {
    V_printf(verbflag, "", "begin record %d of %d\n", ds+1, Nds);

    // mark the time
    wt1 = getwalltime();
    ct1 = getcputime(&ut1, &st1);

    xmRec = xmRD->records[ds];
    xpRec = xpRD->records[ds];
	  
    /*
     * Input keys: time and WCS
     */
    t_rec = drms_getkey_time(xmRec, "T_REC", &status);
    if (status || drms_getkey_time(xpRec, "T_REC", &status) != t_rec) {
      WARN("mgram and pgram input record T_REC's don't match");
      continue;
    }
    // summary
    sprint_time(time_buf, t_rec, "TAI", 0);
    V_printf(verbflag, "\t", "processing images at T_REC = %s\n", time_buf);

    // get WCS
    crvalx = drms_getkey_float(xmRec, "CRVAL1",   &status); // disc center, arcsec
    crvaly = drms_getkey_float(xmRec, "CRVAL2",   &status);
    cdelt  = drms_getkey_float(xmRec, "CDELT1",   &status); // arcsec, assumimg dx=dy
    crpix1 = drms_getkey_float(xmRec, "CRPIX1",   &status); // disk center on ccd
    crpix2 = drms_getkey_float(xmRec, "CRPIX2",   &status);
    crota2 = drms_getkey_float(xmRec, "CROTA2",   &status); // twist
    crlt   = drms_getkey_float(xmRec, "CRLT_OBS", &status); // tilt
    if (status) {
      WARN("could not get WCS from mgram (missing image?). skipping");
      continue;
    }
    if (isnan(crota2) || isnan(crlt)  ||
	isnan(crvalx) || isnan(cdelt) || isnan(crpix1)) {
      WARN("WCS from mgram had at least one NaN (missing image?). skipping");
      continue;
    }
    // Check agreement of mgram and pgram coords; the if() is just for scope
    // (this should be in a function)
    if (1 == 1) {
      float crvalxP, crvalyP, cdeltP, crpix1P, crpix2P, crota2P, crltP;
      const double THR = 1e-4; // comparison threshold, just in case
      
      crvalxP = drms_getkey_float(xpRec, "CRVAL1",   &status); // disc center, arcsec
      crvalyP = drms_getkey_float(xpRec, "CRVAL2",   &status);
      cdeltP  = drms_getkey_float(xpRec, "CDELT1",   &status); // arcsec, assumimg dx=dy
      crpix1P = drms_getkey_float(xpRec, "CRPIX1",   &status); // disk center on ccd
      crpix2P = drms_getkey_float(xpRec, "CRPIX2",   &status);
      crota2P = drms_getkey_float(xpRec, "CROTA2",   &status); // twist
      crltP   = drms_getkey_float(xpRec, "CRLT_OBS", &status); // tilt
      if (status) {
	WARN("could not get WCS from pgram (missing image?). skipping");
	continue;
      }
      if (isnan(crota2P) || isnan(crltP)  ||
	  isnan(crvalxP) || isnan(cdeltP) || isnan(crpix1P)) {
	WARN("WCS from pgram had at least one NaN (missing image?). skipping");
	continue;
      }
      if (fabs(crvalxP - crvalx) > THR || 
	  fabs(crvalyP - crvaly) > THR || 
	  fabs(cdeltP  - cdelt)  > THR ||
	  fabs(crpix1P - crpix1) > THR || 
	  fabs(crpix2P - crpix2) > THR || 
	  fabs(crota2P - crota2) > THR || 
	  fabs(crltP   - crlt)   > THR) {
	WARN("WCS from pgram and mgram did not match. skipping");
	continue;
      }
    }
    
    // Ephemeris courtesy of Phil's macros above
    sina = sin(crota2 * DTOR); 
    cosa = cos(crota2 * DTOR); 
    p0 = -crota2;
    beta = crlt;
    x0 = PIX_X(0.0,0.0) - 1.0; // zero-based coordinates
    y0 = PIX_Y(0.0,0.0) - 1.0;
    rsun_ref = drms_getkey_double(xmRec, "RSUN_REF", &status);
    if (status) {
      rsun_ref = 6.96e8; // this is OK
      status = 0;
    }
    dsun_obs = drms_getkey_double(xmRec, "DSUN_OBS", &status);
    if (status) {
      WARN("could not get DSUN_OBS from mgram. skipping");
      continue;
    }
    rSun = asin(rsun_ref / dsun_obs) * RAD2ARCSEC / cdelt;
    V_printf(verbflag > 1, NULL, 
	     "\t\tX0=%.3f, Y0=%.3f, RSUN=%.3f, P0=%.2f, beta=%.2f\n", 
	     x0, y0, rSun, p0, beta);
    // one final check
    if (isnan(rSun) || isnan(x0) || isnan(y0) || isnan(p0) || isnan(beta)) {
      WARN("Computed disk parameters from mgram had at least one NaN. skipping");
      continue;
    }

    /*
     * Input data arrays
     */
    // mgram data
    xmSeg = drms_segment_lookupnum(xmRec, 0);
    if (!xmSeg) {
      WARN("problem looking up mgram. skipping");
      continue;
    }
    // load as doubles
    xmArray = drms_segment_read(xmSeg, DRMS_TYPE_DOUBLE, &status);
    if (status) {
      drms_free_array(xmArray);
      WARN("problem reading mgram. skipping");
      continue;
    }
    // input size (same for all images)
    M = xmArray->axis[0]; N = xmArray->axis[1];

    // pgram data
    xpSeg = drms_segment_lookupnum(xpRec, 0);
    if (!xpSeg) {
      drms_free_array(xmArray);
      WARN("problem looking up pgram. skipping");
      continue;
    }
    xpArray = drms_segment_read(xpSeg, DRMS_TYPE_DOUBLE, &status);
    if (status) {
      drms_free_array(xmArray);
      drms_free_array(xpArray);
      WARN("problem reading pgram. skipping");
      continue;
    }
    // Check if dimensions match
    if (xpArray->axis[0] != M || xpArray->axis[1] != N) {
      drms_free_array(xmArray);
      drms_free_array(xpArray);
      WARN("pgram and mgram dimensions don't match. skipping");
      continue;
    }

    // create output record, since input segments exist
    V_printf(verbflag, "\t", "create_record\n"); 
    yRec = drms_create_record(drms_env, (char *) yRecQuery, DRMS_PERMANENT, &status);
    if (status) {
      drms_free_array(xmArray);
      drms_free_array(xpArray);
      WARN("Output recordset for labeling could not be created. skipping");
      continue;
    }
	
    /*
     * Set up arguments to core function
     */
    // set up varying arguments to function: xm, xp, geom
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
    // mexfunction now (06/2011) is set up for untransposed input geometry,
    // so no need to flip x0 and y0.
    mxGetPr(prhs[MXT_Hseg_ARG_geom])[0] = x0+1;
    mxGetPr(prhs[MXT_Hseg_ARG_geom])[1] = y0+1;
    mxGetPr(prhs[MXT_Hseg_ARG_geom])[2] = rSun;
    mxGetPr(prhs[MXT_Hseg_ARG_geom])[3] = beta;
    mxGetPr(prhs[MXT_Hseg_ARG_geom])[4] = p0;
 
    // print sample input values
    V_printf(verbflag > 1, NULL, 
	     "\t\tM[%d,%d] = %f\n\t\tP[%d,%d] = %f\n", 
	     M/2, N/2, ((double *) xmArray->data)[M*N/2+M/2],
	     M/2, N/2, ((double *) xpArray->data)[M*N/2+M/2]);

    /*
     * Call the workhorse function
     */
    V_printf(verbflag, "\t", "calling mexfunction `hmi_segment'.\n");
    fflush(stdout);
    msg = main_hmi_segment(nlhs, plhs, nrhs, prhs);
    if (msg == NULL) {
      V_printf(verbflag, "\t", "returned OK from mexfunction.\n");
    } else {
      // mexErrMsgTxt within mexFunction() will route thru this clause
      V_printf(1, "\t", 
	       "returned via error from mexfunction. Message: %s\n" , msg);
      // fail
      status = 1;
      DIE("mexFunction failed (mexErrMsgTxt)");
    }

    // print sample output value
    V_printf(verbflag > 1, "\t", "y[%d,%d] = %g\n", 
	     M/2, N/2, mxGetPr(plhs[MXT_Hseg_ARG_y])[M*N/2+M/2]);
    // check stats size; it will be < #models if there are 0 active pixels
    // (TODO: add Nregions argument to stats-finder to prevent this)
    if (mxGetM(plhs[MXT_Hseg_ARG_s]) != Nclass) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Modelset has %d classes, region stats have %d, %s", 
	       Nclass,  (int) mxGetM(plhs[MXT_Hseg_ARG_s]), 
	       Nclass > mxGetM(plhs[MXT_Hseg_ARG_s]) ? 
	       "can very rarely happen" : "should NOT happen");
      WARN(errmsg);
    }

    // free the input images in the argument list
    // note: we still need the results: y, s, nclean
    // (prep: replace data segments in xm, xp)
    mxSetPr(prhs[MXT_Hseg_ARG_xm], (double *) xm_tmp);
    mxSetPr(prhs[MXT_Hseg_ARG_xp], (double *) xp_tmp);
    // free xm and xp
    V_printf(verbflag > 1, "\t", "within-loop input-image frees\n");
    mxDestroyArray(prhs[MXT_Hseg_ARG_xm]);
    mxDestroyArray(prhs[MXT_Hseg_ARG_xp]);
    // free input DRMS arrays
    drms_free_array(xmArray); 
    drms_free_array(xpArray);

    /*
     * Write labeling out
     */
    // (it is a double here, but defined in the .jsd with BITPIX=8)
    V_printf(verbflag, "\t", "looking up mask segment\n");
    ySeg = drms_segment_lookupnum(yRec, 0); // first segment
    // convert the double mxArray to a char drms_array
    // (this is an independent copy of the output array)
    yArray = mxArray2drms_array(plhs[MXT_Hseg_ARG_y], DRMS_TYPE_CHAR);
    if (!yArray)
      DIE("problem creating drms array for labeling");
    yArray->parent_segment = ySeg; // Link array to segment
    V_printf(verbflag, "\t", "writing mask segment\n");
    status = drms_segment_write(ySeg, yArray, 0);
    if (status)
      DIE("problem writing drms segment for labeling");

    /*
     * Save keywords
     */
    // mask quality keywords
    status = setkey_mask_qual(yRec, 
			      (int) mxGetScalar(plhs[MXT_Hseg_ARG_nclean]), 
			      yArray->data, M, N);
    if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Mask quality keyword insertion failed for %d keys", status);
      WARN(errmsg);
    }
    // Mask statistics keywords
    status = setkey_mask_stats(yRec, mxGetPr(plhs[MXT_Hseg_ARG_s]), modelset);
    if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Mask stats keyword insertion failed for %d keys", status);
      WARN(errmsg);
    } else if (status < 0) {
      WARN("Class with `active' not found.  Failed to insert mask stats keys");
    }
    // Magnetogram keys, including geometry, time, code versions
    status = propagate_keys_mag2mask(xmRec, yRec);
    if (status < 0) {
      WARN("Mag -> Mask key propagation failed, internal error, should not happen");
    } else if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Mag -> Mask key propagation failed for %d keys", status);
      WARN(errmsg);
    }
    // Intensitygram keys, for now just flat-field
    status = propagate_keys_intensity2mask(xpRec, yRec);
    if (status < 0) {
      WARN("Ic -> Mask key propagation failed, internal error, should not happen");
    } else if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Ic -> Mask key propagation failed for %d keys", status);
      WARN(errmsg);
    }
    // Algorithmic-parameter keywords
    status = setkey_algorithm_params(yRec, modelName, prhs);
    if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Algorithm keyword insertion failed for %d keys", status);
      WARN(errmsg);
    }
    // Model-specific keywords
    status = setkey_model_params(yRec, modelset);
    if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Class model keyword insertion failed for %d keys", status);
      WARN(errmsg);
    }
    // Code info keywords
    status = setkey_code_info(yRec);
    if (status > 0) {
      snprintf(errmsg, sizeof(errmsg), 
	       "Code info keyword insertion failed for %d keys", status);
      WARN(errmsg);
    }
    // Essential prime key: T_REC
    drms_copykey(yRec, xmRec, "T_REC");
    // date and build version
    drms_setkey_time(yRec, "DATE", CURRENT_SYSTEM_TIME);
    drms_setkey_string(yRec, "BLD_VERS", jsoc_version);
	  
    // Establish links to source data
    if ((srcLink_M = hcon_lookup_lower(&yRec->links, "MDATA")) != NULL)
      status = drms_setlink_static(yRec, "MDATA", xmRec->recnum);
    if (!srcLink_M) WARN("MDATA link failed: bad lookup");
    if (status)     WARN("MDATA link failed: bad setlink");
    if ((srcLink_P = hcon_lookup_lower(&yRec->links, "PDATA")) != NULL)
      status = drms_setlink_static(yRec, "PDATA", xpRec->recnum);
    if (!srcLink_P) WARN("PDATA link failed: bad lookup");
    if (status)     WARN("PDATA link failed: bad setlink");
            
    // close mask record
    V_printf(verbflag, "\t", "inserting mask record\n"); 
    drms_close_record(yRec, DRMS_INSERT_RECORD);
    ds_good++; // one more was OK

    // free the other results in the argument list (inputs already freed)
    V_printf(verbflag > 1, "\t", "within-loop result frees\n");
    mxDestroyArray(plhs[MXT_Hseg_ARG_y]);
    mxDestroyArray(plhs[MXT_Hseg_ARG_s]);
    mxDestroyArray(plhs[MXT_Hseg_ARG_post]);
    mxDestroyArray(plhs[MXT_Hseg_ARG_nclean]);
    // free output DRMS array
    drms_free_array(yArray);

    // Timing info
    wt = getwalltime();
    ct = getcputime(&ut, &st);
    V_printf(verbflag, "\t", "record %d time used: %.3f s wall, %.3f s cpu\n", 
	     ds+1, (wt - wt1)*1e-3, (ct - ct1)*1e-3);
    V_printf(verbflag, "\t", "record %d done.\n", ds+1);
  } 
  V_printf(verbflag, "", "close records (mgram)\n");
  drms_close_records(xmRD, DRMS_FREE_RECORD);
  V_printf(verbflag, "", "close records (pgram)\n");
  drms_close_records(xpRD, DRMS_FREE_RECORD); // if pgram == mgram, this segfaults
	
  // free the unchanging arguments
  V_printf(verbflag, "", "final destroy arrays\n");
  mxDestroyArray(prhs[MXT_Hseg_ARG_edge]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_iter]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_T]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_beta]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_alpha]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_geom]);
  mxDestroyArray(prhs[MXT_Hseg_ARG_rho]);
  for (class = 0; class < Nclass; class++)
    mxDestroyArray(prhs[MXT_Hseg_ARG_m1 + class]);

  // success record
  V_printf(verbflag, "", "successfully added %d of %d masks.\n", ds_good, Nds);
  if (ds_good < Nds)
    WARN("Some masks were not able to be computed from the input series");
  // overall timing
  wt = getwalltime();
  ct = getcputime(&ut, &st);
  V_printf(verbflag, "", "total time used: %.3f s wall, %.3f s cpu\n", 
	   (wt - wt0)*1e-3, (ct - ct0)*1e-3);

  // exit OK
  return DRMS_SUCCESS;
}

