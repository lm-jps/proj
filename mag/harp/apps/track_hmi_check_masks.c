/*
 *
 * track_hmi_check_masks: pre-run mask availability checker
 *
 * Module name: track_hmi_check_masks
 *
 * This jsoc module checks that a mask is present where-ever one could
 * be, depending on the source (mag,pgram) data sources.  It counts 
 * "unexcused absences" of masks, and if there are more than a given
 * maximum (default 0, see "thresh" argument), exits with an error code,
 * otherwise, exits cleanly.
 *
 * It is intended to be automatically run from the master HARP script, 
 * before the tracker itself is called.
 * 
 * Arguments xm (mags), xp (pgrams), and y (masks) are required.  The
 * mags and pgrams are data series names, and the masks are a recordset 
 * query giving the T_REC interval.
 *
 * For definitive:
 *   hmi.M_720s, hmi.Ic_noLimbDark_720s, and hmi.Marmask_720s
 * For NRT:
 *   use corresponding NRT series
 *
 * + Optional thresh argument, specifying the maximum number of
 *   missing masks that will trigger an "OK" status code (0) upon exit;
 *   larger than this number will trigger an error (nonzero) exit status.
 *   If thresh is a fraction strictly less than 1, it is interpreted as a 
 *   proportion of the total number of mags in the interval (rounded 
 *   down).  If 1 or larger, it is interpreted as an absolute number
 *   of masks.  Thus, if thresh=0 (default), any missing masks 
 *   at all will cause error exit status.
 * + Optional level argument, specifying whether presence of data
 *   segments is to be checked.
 * + Optional verbose argument (VERB=v where v is 0, 1, 2, or 3).
 *
 * Current typical usage:
 *
 * track_hmi_check_masks
 *      xm='hmi.M_720s'
 *      xp='hmi.Ic_noLimbDark_720s'
 *      mask='hmi.Marmask_720s[2010.07.03_TAI/1d]' 
 *
 * Michael Turmon, JPL
 *   apr 2013: created  
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

// set according to verbosity input parameter
static int verbflag;

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
	fflush(stderr); \
        } while (0)
// V_printf: facilitate verbose output
//   if flag is > 0, output is to stdout, if < 0, to stderr
//   if flag is 0, no output is made at all
// The message is printed in the form:
// <first><module_name>: <message>"
// Usage:
// V_printf(VERB > 0, "\t", "Mask(%d) = %d\n", 2048, mask[2048]);
void
V_printf(int flag, char *first, char *format, ...) {
  va_list args;
  extern char *module_name;
  FILE *fp = (flag > 0) ? stdout : stderr;

  va_start(args, format);
  if (flag != 0) {
    // first is a string, even "" -- print the module name too
    // otherwise, omit it
    if (first)
      fprintf(fp, "%s%s: ", first, module_name);
    vfprintf(fp, format, args);
    fflush(fp);
  }
  va_end(args);
}

// standard strings
#define STR_MAX 256

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

/***********************************************************************
 *
 * Declarations
 *
 ***********************************************************************/

/*
 * WCS, in the same parameterization as in HMI keywords
 */
typedef struct {
  int ok;           // 0 if not OK
  char missing[32]; // short message indicating missing data, if !ok
  // types here are double or float according to hmi.M_720s types
  double rsun_ref; // RSUN_REF || 6.96e8
  double dsun_obs; // DSUN_OBS
  float  cdelt1;   // CDELT1
  float  crval1;   // CRVAL1 -- disc center, arcsec
  float  crval2;   // CRVAL2
  float  crpix1;   // CRPIX1 -- center, ccd, orig.
  float  crpix2;   // CRPIX2
  float  crota2;   // CROTA2
  float  crlt_obs; // CRLT_OBS (converted to radians)
  float  crln_obs; // CRLN_OBS (converted to radians)
} wcs_t;


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

/***********************************************************************
 *
 * Utility routines
 *
 ***********************************************************************/


#define max(a,b) (((a) > (b)) ? (a) : (b))

static
double
wcs_diff(wcs_t *wcs1, wcs_t *wcs2)
{
  double top = 0;
  // skip dsun_obs: it's too big, and was not in original 
  // mask maker "bad-wcs" checker
  top = max(top, fabs(wcs1->cdelt1   - wcs2->cdelt1  ));
  top = max(top, fabs(wcs1->crval1   - wcs2->crval1  ));
  top = max(top, fabs(wcs1->crval2   - wcs2->crval2  ));
  top = max(top, fabs(wcs1->crpix1   - wcs2->crpix1  ));
  top = max(top, fabs(wcs1->crpix2   - wcs2->crpix2  ));
  top = max(top, fabs(wcs1->crota2   - wcs2->crota2  ));
  top = max(top, fabs(wcs1->crlt_obs - wcs2->crlt_obs));
  // also was not in original checker
  // top = max(top, fabs(wcs1->crln_obs - wcs2->crln_obs));
  return top;
}

// clean namespace
#undef max

/*
 * get ephemeris into *wcs from a DRMS record
 *   returns 0 for success, nonzero for failed getkey or NaN in result
 *   note that all fields in *wcs are named exactly as the corresponding
 *   keyword, except in lower-case
 */

static
int 
wcs_load(wcs_t *wcs, DRMS_Record_t *rec)
{
  int status, retval;
  const size_t Nmiss = sizeof(wcs->missing);

  // RSUN_REF and DSUN_OBS are doubles in the JSD, rest are floats
  status = 0;
  wcs->rsun_ref = drms_getkey_double(rec, "RSUN_REF", &status);
  if (status) wcs->rsun_ref = 6.96e8;
  status = 0;
  wcs->dsun_obs = drms_getkey_double(rec, "DSUN_OBS", &status);
  wcs->cdelt1   = drms_getkey_float (rec, "CDELT1",   &status);  // arcsec, if dx=dy
  wcs->crval1   = drms_getkey_float (rec, "CRVAL1",   &status);  // disc center, arcsec
  wcs->crval2   = drms_getkey_float (rec, "CRVAL2",   &status);
  wcs->crpix1   = drms_getkey_float (rec, "CRPIX1",   &status);  // center, ccd, orig.
  wcs->crpix2   = drms_getkey_float (rec, "CRPIX2",   &status);
  wcs->crota2   = drms_getkey_float (rec, "CROTA2",   &status);
  wcs->crlt_obs = drms_getkey_float (rec, "CRLT_OBS", &status);
  wcs->crln_obs = drms_getkey_float (rec, "CRLN_OBS", &status);
  // assume set valid status...
  strcpy(wcs->missing, "");
  retval = 0;
  // ...and alter as necessary
  if (status != 0) {
    snprintf(wcs->missing, Nmiss, "Failed getkey()");
    retval = 1;
  } else {
    if (isnan(wcs->dsun_obs)) snprintf(wcs->missing, Nmiss, "NaN %s", "DSUN_OBS");
    if (isnan(wcs->cdelt1  )) snprintf(wcs->missing, Nmiss, "NaN %s", "CDELT1");
    if (isnan(wcs->crval1  )) snprintf(wcs->missing, Nmiss, "NaN %s", "CRVAL1");
    if (isnan(wcs->crval2  )) snprintf(wcs->missing, Nmiss, "NaN %s", "CRVAL2");
    if (isnan(wcs->crpix1  )) snprintf(wcs->missing, Nmiss, "NaN %s", "CRPIX1");
    if (isnan(wcs->crpix2  )) snprintf(wcs->missing, Nmiss, "NaN %s", "CRPIX2");
    if (isnan(wcs->crota2  )) snprintf(wcs->missing, Nmiss, "NaN %s", "CROTA2");
    if (isnan(wcs->crlt_obs)) snprintf(wcs->missing, Nmiss, "NaN %s", "CRLT_OBS");
    if (isnan(wcs->crln_obs)) snprintf(wcs->missing, Nmiss, "NaN %s", "CRLN_OBS");
    // if any of the above conditions fire, wcs->missing will be non-empty
    if (strlen(wcs->missing) > 0) 
      retval = 2;
  }
  // if (retval) V_printf(verbflag>1, "wcs_load: ", "WCS fail: %s\n", wcs->missing);
  return retval;
}

/*
 * open drms record, optionally load WCS, look-up segment, load segment
 *   i.e., general-purpose drms data-getter/checker
 * Usage: Call, for example, with query = 'hmi.M_720s' and trec = '2012.05.01_TAI'.
 *
 * Always: opens the corresponding record set, makes sure it has exactly one 
 *   matching record, and checks its quality for missing-ness (0x80000000).
 * Typically returns the RecordSet in *recSetP (it is returned because it must 
 *   be freed) and the single record in *recP (which becomes invalid once *recSetP 
 *   is freed).
 * Option 0a: If recP is NULL, finds but does not return recP.  
 * Option 0b: If, further to (0a) above, recSetP is NULL, finds but does not 
 *   return recSetP.
 * In either option 0 case, the error status is returned, allowing data-present 
 *   checks.  Note, you don't use 0b without 0a.
 * Option 1: If wcs is non-null, looks up and fills in wcs_t structure.  If this
 *   fails or WCS are partly missing, it's an error.
 * Option 2a: If segname is non-null, looks up the named segment 
 *   within the record, to verify it is there.  Nothing extra is returned.
 * Option 2b: If, further to (2a) above, dataP is non-null, loads 
 *   the segment according to d_typ, storing it in *dataP.  Again, *dataP 
 *   becomes invalid once *recSetP is freed.  
 *   The data in the *dataP is in (*dataP)->data.  
 *   If unneeded, you can let d_typ = DRMS_TYPE_RAW, which is noncommittal.
 *
 * Note:
 *   If you want the data (dataP) or the record (recP), you *have* to get the 
 *   recSet, because that's what you free.  If you want the WCS, you don't
 *   need to get the recSet.
 *   
 * If success, returns NULL.  If failure, returns a descriptive string, and in this
 *   case, nothing needs to be freed.  If success, be sure to:
 *     drms_close_records(*recSetP, DRMS_FREE_RECORD);
 *   when you are done, if you requested the recSet.
 */

static
char *
load_drms_rec_seg(DRMS_RecordSet_t **recSetP, 
		  DRMS_Record_t       **recP, 
		  DRMS_Array_t       **dataP, 
		  wcs_t *wcs,
		  // remainder are inputs
		  const char *query, 
		  const char *trec,
		  const char *segname,
		  DRMS_Type_t d_typ)
{
  char rec_query[STR_MAX];
  char tag[STR_MAX];
  int  status, return_code;
  DRMS_RecordSet_t *the_recSet;
  DRMS_Record_t    *the_rec;
  DRMS_Array_t     *the_data;
  DRMS_Segment_t   *the_seg;
  
  // exclude mistakes
  if (!query || !trec)
    return "Illegal argument to load_drms_rec_seg: need query and trec";
  if (!recSetP && (recP || dataP)) 
    return "Illegal argument to load_drms_rec_seg: recSet required here";
  if (!segname && dataP)
    return "Illegal argument to load_drms_rec_seg: segname required for data";
  // null out output pointers as permitted, in case of early exit
  if (recSetP) *recSetP = NULL;
  if (recP)    *recP    = NULL;
  if (dataP)   *dataP   = NULL;
  if (wcs)     {wcs->ok = 0; strcpy(wcs->missing, ""); } // no message for now
  // descriptive tag
  snprintf(tag, sizeof(tag), "%s %s", 
	   segname ? "and segment" : "(no segment)",
	   dataP ? "and array" : "");
  // query for data
  snprintf(rec_query, sizeof(rec_query), "%s[%s]", query, trec);
  V_printf(verbflag > 2, "\t\t", "Opening record %s %s\n", rec_query, tag);
  status = 0;
  the_recSet = drms_open_records(drms_env, rec_query, &status);
  if (status != 0) {
    // (nothing to close)
    return "Could not open record at this T_REC";
  }
  // (below here, we close recSet before error out)
  if (the_recSet->n != 1) {
    drms_close_records(the_recSet, DRMS_FREE_RECORD);
    return "Could not open single record at this T_REC";
  }
  the_rec = the_recSet->records[0];
  // quality check to identify missing data 
  int quality = drms_getkey_int(the_rec, "QUALITY", &status);
  if (status != 0 || (quality & 0x80000000)) {
    drms_close_records(the_recSet, DRMS_FREE_RECORD);
    return "Missing data at this T_REC (quality)";
  }
  // option 1: load WCS
  if (wcs) {
    return_code = wcs_load(wcs, the_rec);
    if (return_code != 0) {
      drms_close_records(the_recSet, DRMS_FREE_RECORD);
      return "Missing WCS at this T_REC";
    }
  }
  // option 2a: look up segname
  if (segname) {
    the_seg = drms_segment_lookup(the_rec, segname);
    if (the_seg == NULL) {
      drms_close_records(the_recSet, DRMS_FREE_RECORD);
      return "Failed data segment lookup at this T_REC";
    }
  }
  // option 2b: data
  if (segname && dataP) {
    // (know the_seg is ok from above)
    the_data = drms_segment_read(the_seg, d_typ, &status);
    if (the_data == NULL || status != 0) {
      drms_close_records(the_recSet, DRMS_FREE_RECORD);
      return "Failed data segment read at this T_REC";
    }
  }
  // all OK: set up outputs as allowed
  if (recSetP) *recSetP = the_recSet;
  if (recP)    *recP    = the_rec;
  if (dataP)   *dataP   = the_data;
  // OTOH: close the record set if it was not returned
  if (!recSetP)
    drms_close_records(the_recSet, DRMS_FREE_RECORD);
  return NULL;
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
#define kSeriesXM    "xm"
#define kSeriesXP    "xp"
#define kQueryMask   "y"
#define kThresh      "thresh"
#define kLevel       "level"
#define kVerb        "VERB"

char *module_name = "track_hmi_check_masks";

ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kSeriesXM,   "",   "M-gram data series."},
   {ARG_STRING,  kSeriesXP,   "",   "P-gram data series."},
   {ARG_STRING,  kQueryMask,  "",   "Mask recordset query."},
   {ARG_DOUBLE,  kThresh,    "0",   "Threshold: count, or proportion in [0,1)"},
   {ARG_INT,     kLevel,     "1",   "Check level (unused at present)"},
   {ARG_INT,     kVerb,      "1",   "Verbosity: 0=errs only; 1, 2, 3 = more"},
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
  CmdParams_t    *params = &cmdparams;
  const char   *xmSeries = params_get_str(params, kSeriesXM);
  const char   *xpSeries = params_get_str(params, kSeriesXP);
  const char  *maskQuery = params_get_str(params, kQueryMask);
  const double  threshIn = params_get_double(params, kThresh);
  const int       level  = params_get_int(params, kLevel); // an int -- unused
  // related queries/series
  char *maskTrecTail;
  char maskSeries[STR_MAX];
  char magQuery[STR_MAX];
  // DRMS data
  DRMS_RecordSet_t *xmRecSet;  // for iterating over
  // counts
  int ds, Nds; // iteration
  int ds_ok_source, ds_ok_mask, ds_bad; // data record counts
  int thresh;
  // trec and wcs keys, etc.
  wcs_t mag_wcs, pgram_wcs;
  TIME trec;
  char trecStr[64];
  char trec_save[64];
  char *err;
  char err_save[STR_MAX];
  // Time measuring -- ugly, but functional
  double wt0, wt1, wt;
  double ut0, ut1, ut;
  double st0, st1, st;
  double ct0, ct1, ct;

  // static, global to this file
  verbflag = params_get_int(params, kVerb); // an int, not a flag
  
  V_printf(verbflag, "", "Beginning.\n");
  // timings
  wt0 = getwalltime();
  ct0 = getcputime(&ut0, &st0);

  // extract T_REC part of maskQuery
  if (strchr(xmSeries, '[') != NULL)
    WARN("Got T_REC qualifier [...] in mag input series, probably an error");
  if (strchr(xpSeries, '[') != NULL)
    WARN("Got T_REC qualifier [...] in pgram input series, probably an error");
  if ((maskTrecTail = strchr(maskQuery, '[')) == NULL)
    DIE("Need T_REC qualifier [...] in y (mask) input record query");
  // make mag query: hmi.M_720s[...]
  snprintf(magQuery, sizeof(magQuery), "%s%s", xmSeries, maskTrecTail);
  // re-write mask query without the T_REC part
  strncpy(maskSeries, "", sizeof(maskSeries)); // copy \0 into maskSeries
  strncpy(maskSeries, maskQuery, maskTrecTail-maskQuery); // copy up to the [

  // open mgrams, which are our base for iteration
  xmRecSet = drms_open_records(drms_env, (char *) magQuery, &status);
  if (status) {
    V_printf(-1, "", "Fatal: could not open mags for masks (%s)\n", magQuery);
    DIE("No matching mgram input data found");
  }
  Nds = xmRecSet->n; // 0 is OK here

  // set threshold for success
  if (threshIn >= 1.0) {
    // a raw count
    thresh = (int) threshIn;
  } else {
    // a proportion
    thresh = (int) (Nds * threshIn);
  }

  // loop over records
  V_printf(verbflag, "", "Examining sources and masks at %d T_RECs.\n", Nds);
  ds_ok_source = ds_ok_mask = 0;
  for (ds = 0; ds < Nds; ds++) {
    // mark the time
    wt1 = getwalltime();
    ct1 = getcputime(&ut1, &st1);
    // get t_rec
    status = 0;
    trec = drms_getkey_time(xmRecSet->records[ds], "T_REC", &status); 
    if (status) {
      V_printf(-1, "", "No T_REC at record %d of %d, skipping.\n", ds, Nds);
      continue;
    }
    sprint_time(trecStr, trec, "TAI", 0);
    V_printf(verbflag > 1, "\t", "[%d/%d] %s\n", ds, Nds, trecStr);

    // check if mag is there
    err = load_drms_rec_seg(NULL, NULL, NULL, &mag_wcs,
			    xmSeries, trecStr, "magnetogram", DRMS_TYPE_RAW);
    V_printf(err && (verbflag > 1), "\t\t", 
	     "%s[%s]: %s %s-- mag missing is OK.\n", xmSeries, trecStr, err,
	     mag_wcs.ok ? "" : mag_wcs.missing);
    // ordinary missing mag: not a problem
    if (err) continue;

    // check if mag is there
    err = load_drms_rec_seg(NULL, NULL, NULL, &pgram_wcs,
			    xpSeries, trecStr, "continuum", DRMS_TYPE_RAW);
    V_printf(err && (verbflag > 1), "\t\t", 
	     "%s[%s]: %s %s-- pgram missing is OK.\n", xpSeries, trecStr, err,
	     pgram_wcs.ok ? "" : pgram_wcs.missing);
    // ordinary missing photogram: not a problem
    if (err) continue;

    // Check agreement of mgram and pgram coords
    //   disagreement is surprising, but OK, it will prevent the mask being made
    //   the 1e-4 threshold matches the value in hmi_segment_module
    if (wcs_diff(&mag_wcs, &pgram_wcs) > 1e-4) {
      err = "WCS mismatch between pgram/mgram at same T_REC";
      V_printf(verbflag > 1, "\t", 
	       "%s[%s]: %s -- OK (but surprising).\n", xmSeries, trecStr, err);
      continue;
    }

    // check if mask is there
    err = load_drms_rec_seg(NULL, NULL, NULL, NULL, 
			    maskSeries, trecStr, "mask", DRMS_TYPE_RAW);
    // note, lower threshold for reporting than data sources
    V_printf(err && (verbflag > 0), "\t\t", 
	     "%s[%s]: %s -- Error: mask missing.\n", maskSeries, trecStr, err);

    // always bump the "OK source" count
    ds_ok_source++;
    // bump "ok mask" only if no error
    if (!err) ds_ok_mask++;
    // save one trec/error for info later
    if (err) {
      snprintf(err_save, sizeof(err_save), "%s", err);
      snprintf(trec_save, sizeof(trec_save), "%s", trecStr);
    }
      
    // Timing info
    wt = getwalltime();
    ct = getcputime(&ut, &st);
    V_printf(verbflag > 2, "\t\t", "record %d time used: %.3f s wall, %.3f s cpu\n", 
	     ds+1, (wt - wt1)*1e-3, (ct - ct1)*1e-3);
  } 
  V_printf(verbflag > 2, "", "close records (mgram)\n");
  drms_close_records(xmRecSet, DRMS_FREE_RECORD);
  
  // success record
  ds_bad = ds_ok_source - ds_ok_mask;
  V_printf(verbflag, "", "Examined NT = %d T_RECs.\n", Nds);
  V_printf(verbflag, "", "  Of NT, found %d valid (mag,pgram) pairs.\n", ds_ok_source);
  V_printf(verbflag, "", "  Of NT, found %d valid masks.\n", ds_ok_mask);
  if (ds_bad > 0) {
    // always say something if we're missing masks
    V_printf(1, "", "Some masks (%d) were missing although sources were present.\n", ds_bad);
    V_printf(1, "", "\tExample: %s[%s] -- %s.\n", maskSeries, trec_save, err_save);
  } else {
    V_printf(1, "", "No unexplained absences of masks.\n");
  }
  if (ds_bad <= thresh) {
    V_printf(1, "", "Overall mask check passed: %d missing <= %d threshold.\n",
	     ds_bad, thresh);
  } else {
    V_printf(1, "", "Mask check failed: %d missing > %d threshold\n",
	     ds_bad, thresh);
  }

  // overall timing
  wt = getwalltime();
  ct = getcputime(&ut, &st);
  V_printf(verbflag, "", "total time used: %.3f s wall, %.3f s cpu\n", 
	   (wt - wt0)*1e-3, (ct - ct0)*1e-3);

  V_printf(verbflag, "", "Done.\n");
  // exit status depends on threshold comparison
  if (ds_bad <= thresh) {
    return DRMS_SUCCESS;
  } else {
    // anything nonzero is OK, and pre-defined codes not apt
    return 1;
  }
}

