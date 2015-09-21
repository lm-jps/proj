/*
 *  datavg.c						~rick/src/drms
 *
 *  Construct mean and variance matrices of selected data segments from an
 *    input record set (same format)
 *  This is a DRMS version of SOI module doppavg generalized to a
 *    larger variety of data series and types
 *
 *  Responsible:  Rick Bogart				RBogart@spd.aas.org
 *
 *  Usage:
 *    datavg in= out= 
 *
 *  Parameters: (type   default         description)
 *	in	str	-		Input dataset or data series;
 *			if only a data series name, additional selection
 *			parameters will be required
 *      out	str	-		Output data series
 *      count	str	valid		Output series segment containing count
 *      mean	str	mean		Output series segment containing mean
 *      power	str	power		Output series segment containing variance
 *      log	str	Log		Output series segment containing run log
 *      tmid	str	unspecified	midpoint of averaging interval (CR:CL)
 *	length	float	unspecified	length of averaging interval (in deg og
 *			Carrington rotation)
 *	reject	string	Not Specified	If specified, name of a file with a
 *				list of input images to be rejected regardless
 *				of quality mask values
 *	qmask	int	0x80000000	Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	cvok	str	any		64-bit mask of acceptable values for
 *				calibration veraion of input; default -> all
 *				values acceptable
 *	cvno	str	none		64-bit mask of acceptable values for
 *				calibration veraion of input; default -> no
 *				values unacceptable
 *	copy	str	+		list of keys to be propagated as-is
 *	average	str	+		list of keys to be averaged
 *	pkey	str	T_OBS		Name of keyword to be used as index
 *			over which input records are selected for averaging
 *	qual_key str	Quality		Key name of uint type keyword describing
 *			data quality
 *	cvkey	str	CalVer64	Key name of 64-bit mask type keyword
 *			describing calibration version used
 *	roll_key str	CROTA2		Name of keyword describing roll angle
 *	pa	float	180.0		Centre of acceptable roll angles
 *	dpa	float	1.0		Maximum deviation of acceptable roll
 *			angles from nominal centre
 *	mscale	float	(SEG default)	Output BSCALE factor for mean
 *	mzero	float	(SEG default)	Output BZERO offset for mean
 *	pscale	float	(SEG default)	Output BSCALE factor for power
 *	pzero	float	(SEG default)	Output BZERO offset for power
 *	setkey	str	Not Specified	Name of special keyword to be set
 *	setval	double	Not Specified	Value for special keyword to be set, or
 *			name of key whose value is to be used
 *
 *  Flags
 *	-n	no ingestion of output; diagnostic only
 *	-o	correct individual images for orbital velocity
 *	-v	run verbose
 *
 *  Bugs:
 *    Although drms_segment_write() errors are flagged in stderr, they do
 *	not force an abort; consequently it is possible to have records
 *	created with incomplete SUMS data, or no SUMS data at all.
 *    Propagated keys are copied from the first input record, without checking
 *	for uniqueness
 *    If the input dataset is specified in the in string, the keywords for
 *	CarRot and CMLon are set to garbage, and there is no checking that
 *	the length is correct
 *    The MidTime value is geoecentric for midpoint Carrington longitude
 *    The correction for observer velocity only takes account of the radial
 *	component, and is only appropriate for Doppler data
 *    For log-base data, the LOG_BASE keyword is propagated to the output
 *	record, but it really only applies to the mean and variance; it should
 *	be a per-segment keyword
 *    The check on nominal roll angle and the averagang of CROTA* values
 *	assumes that the units are degrees
 *    The reported extremal values of mean and power are determined before
 *	possible output scaling; they may be outside the representable range
 *    Only a single value can be specified for cvok or cvno (which are treated
 *	as string arguments for local parsing, although in principle they
 *	should be arrays of ints)
 *    The protocol of the log segment, if present in the output series, is not
 *	verified to be generic
 *    The values of the primekey in the input set, as specified by pkey, are
 *	automatically averaged into a corresponding key (and D_key) in the
 *	output data set. If the default value of T_OBS is overridden, must
 *	provide appropriate keys in output series if desired; alternatively,
 *	T_OBS should be added to average list.
 *
 *  Revision history is at the end of the file.
 *
 */

#include <jsoc_main.h>
#include "keystuff.c"
#include "earth_ephem.c"

#define	DO_NOTHING	(0)
#define	DO_ABSVAL	(1)
#define	DO_SQUARE	(2)

char *module_name = "data average";
char *version_id = "1.5";

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"in", "", "input data series or dataset"}, 
  {ARG_STRING,	"out", "", "output data series"}, 
  {ARG_STRING,	"tmid", "Not Specified",
      "midpoint of averaging interval (in Carrington longitude)"}, 
  {ARG_FLOAT,	"length", "Not Specified",
      "length of averaging interval (in degrees of Carrington rotation)"}, 
  {ARG_FLOAT,	"pa", "180.0", "centre of acceptable roll angles"}, 
  {ARG_FLOAT,	"dpa", "1.0", "maximum deviation of acceptable roll angles"}, 
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_STRING,	"cvkey", "CalVer64", "keyname for Calibration Version key"}, 
  {ARG_STRING,  "cvok", "any", "Acceptable value of cvkey"},
  {ARG_STRING,  "cvno", "none", "Unacceptable value of cvkey"},
  {ARG_STRING,  "copy",  "+", "comma separated list of keys to propagate"},
  {ARG_STRING,  "average",  "+", "comma separated list of keys to average"},
  {ARG_STRING,	"pkey", "T_OBS",
      "keyname for index over which records are selected for averaging"}, 
  {ARG_STRING,	"qual_key", "Quality", "keyname for 32-bit image quality field"}, 
  {ARG_STRING,	"roll_key", "CROTA2", "keyname for input roll-angle field"}, 
  {ARG_STRING,	"count", "valid", "output data series segment containing count"}, 
  {ARG_STRING,	"mean", "mean", "output data series segment containing mean"}, 
  {ARG_STRING,	"power", "power", "output data series segment containing variance"}, 
  {ARG_STRING,	"log", "Log", "output data series segment containing run log"}, 
  {ARG_STRING,	"setkey", "Not Specified", "name of special extra key to be set"}, 
  {ARG_DOUBLE,	"setval", "Not Specified",
	"value of special extra key to be set; if invalid, name of key whose value is to be used"}, 
  {ARG_FLOAT,	"mscale", "Segment Default", "output BSCALE factor for mean"},
  {ARG_FLOAT,	"mzero", "Segment Default", "output BZERO offset for mean"},
  {ARG_FLOAT,	"pscale", "Segment Default", "output BSCALE factor for power"},
  {ARG_FLOAT,	"pzero", "Segment Default", "output BZERO offset for power"},
  {ARG_FLAG,	"n",	"", "no output records produced; diagnostics only"}, 
  {ARG_FLAG,	"o",	"", "remove effects of observer velocity"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {ARG_END}
};

	/*  list of keywords to propagate by default (if possible and unique)
						        from input to output  */
char *propagate[] = {"TELESCOP", "INSTRUME", "WCSNAME", "WCSAXES"};
		 /*  list of keywords to average by default (if possible)
						        from input to output  */
char *average[] = {"T_OBS", "OBS_VR", "OBS_VW", "OBS_VN"};
			  /*  the following belong in external utility files  */
static long long params_get_mask (CmdParams_t *params, char *arg,
    long long defval) {
/*
 *  This function parses the string associated with the command parameters
 *    argument "arg" as the hexadecimal representation of an unsigned 64-bit
 *    integer, which it returns. The string may consist of up to 16 hexadecimal
 *    characters. (They can optionally be preceded by '0x', but the string is
 *    treated as a hexadecimal representation regardless.)  If there are any
 *    extra or illegal characters in the string, the value "defval" is returned.
 */
  long long retval;
  const char *str = params_get_str (params, arg);
  char *ext;

  retval  = strtoull (str, &ext, 16);
  if (strlen (ext)) retval = defval;
  return retval;
}

static int fgetline (FILE *in, char *line, int max) {
  if (fgets (line, max, in) == NULL) return 0;
  else return (strlen (line));
}

static int read_reject_list (FILE *file, int **list) {
  int ds, sn, rec, last_rec;
  int allocd = 1024, ct = 0, gap = 0;
  char line[1024], t_str[64], estr[16];

  *list = (int *)malloc (allocd * sizeof (int));
  while (fgetline (file, line, 1024)) {
    if (strlen (line) == 1) continue;
    if (line[0] == '#') continue;
    if (sscanf (line, "%d %d %d %s", &ds, &sn, &rec, t_str) != 4) {
      if (sscanf (line, "%d %s", &rec, t_str) != 2) {
	sscanf (line, "%s", estr);
	if (strcmp (estr, "...")) continue;
	gap = 1;
	last_rec = rec;
	continue;
      }
    }
    if (gap) {
      while (rec > last_rec) {
	last_rec++;
	(*list)[ct++] = last_rec;
	if (ct >= allocd) {
	  allocd += 1024;
	  *list = (int *)realloc (*list, allocd * sizeof (int));
	}
      }
      gap = 0;
      continue;
    }
    (*list)[ct++] = rec;
    if (ct >= allocd) {
      allocd += 1024;
      *list = (int *)realloc (*list, allocd * sizeof (int));
    }
  }
  return ct;
}

static int key_params_from_dspec (const char *dspec) {
/*
 *  Establish whether input data records are determined from dataset specifier
 *  assume that if a bracket is found in the dataset specifier it is a
 *  set of well-specified records containing a properly ordered input dataset;
 *  otherwise, a dataseries to be queried
 */
  int n, nt = strlen (dspec);

  for (n = 0; n < nt; n++) if (dspec[n] == '[') return 1;
  return 0;
}

int set_stats_keys (DRMS_Record_t *rec, DRMS_Array_t *vcts, DRMS_Array_t *mean,
    DRMS_Array_t *powr, long long ntot, int mscaled, int pscaled,
    double mrepmin, double mrepmax, double prepmin, double prepmax, int verbose) {
  double *vavg = (double *)mean->data;
  double *vvar = (powr) ? (double *)powr->data : NULL;
  int *vval = (int *)vcts->data;
  double mnvmin, mnvmax, pwrmin, pwrmax, norm, norm2, scale, sig, var;
  double sumv1, sump1, sum2, sum3, sum4;
  int sumv0, sump0;
  int valmin = ntot, valmax = 0;
  int n;
  int vminloc, vmaxloc, pminloc, pmaxloc;
  int morhict, morloct, morct, porhict, porloct, porct;
  int kstat = 0;

  mnvmin = pwrmin = DBL_MAX;
  mnvmax = pwrmax = -DBL_MAX;
  sumv1 = sump1 = 0.0;
  sumv0 = sump0 = 0;
  if (mscaled || pscaled) {
    morhict = morloct = 0;
    porhict = porloct = 0;
  }

  for (n = 0; n < ntot; n++) {
    if (vval[n] < valmin) valmin = vval[n];
    if (vval[n] > valmax) valmax = vval[n];
    if (vavg[n] < mnvmin) {
      mnvmin = vavg[n];
      vminloc = n;
    }
    if (vavg[n] > mnvmax) {
      mnvmax = vavg[n];
      vmaxloc = n;
    }
    if (isfinite (vavg[n])) {
      sumv0++;
      sumv1 += vavg[n];
    }
    if (mscaled) {
      if (vavg[n] < mrepmin)  morloct++;
      if (vavg[n] > mrepmax)  morhict++;
    }
    if (vvar) {
      if (vvar[n] < pwrmin) {
	pwrmin = vvar[n];
	pminloc = n;
      }
      if (vvar[n] > pwrmax) {
	pwrmax = vvar[n];
	pmaxloc = n;
      }
      if (isfinite (vvar[n])) {
	sump0++;
	sump1 += vvar[n];
      }
      if (pscaled) {
	if (vvar[n] < prepmin)  porloct++;
	if (vvar[n] > prepmax)  porhict++;
      }
    }
  }
  
  kstat += check_and_set_key_int (rec, "CountMIN", valmin);
  kstat += check_and_set_key_int (rec, "CountMAX", valmax);
  kstat += check_and_set_key_double (rec, "MeanMIN", mnvmin);
  kstat += check_and_set_key_double (rec, "MeanMAX", mnvmax);
  if (vvar) {
    kstat += check_and_set_key_double (rec, "PowrMIN", pwrmin);
    kstat += check_and_set_key_double (rec, "PowrMAX", pwrmax);
  }

  if (sumv0) {
    scale = 1.0 / sumv0;
    sumv1 *= scale;
    sum2 = sum3 = sum4 = 0.0;
    kstat += check_and_set_key_double (rec, "MeanAVG", sumv1);
    for (n = 0; n < ntot; n++) {
      if (isfinite (vavg[n])) {
        norm = vavg[n] - sumv1;
	norm2 = norm * norm;
	sum2 += norm2;
	sum3 += norm * norm2;
	sum4 += norm2 * norm2;
      }
    }
    kstat += check_and_set_key_double (rec, "MeanRMS", sqrt (sum2 / sumv0));
    if (sumv0 > 1) {
      var = sum2 / (sumv0 - 1);
      sig = sqrt (var);
      kstat += check_and_set_key_double (rec, "MeanSKEW",
	  scale * sum3 / (var * sig));
      kstat += check_and_set_key_double (rec, "MeanKURT",
	  scale * sum4 / (var * var) - 3.0);
    }
  }

  if (sump0) {
    scale = 1.0 / sump0;
    sump1 *= scale;
    sum2 = sum3 = sum4 = 0.0;
    kstat += check_and_set_key_double (rec, "PowrAVG", sump1);
    for (n = 0; n < ntot; n++) {
      if (isfinite (vvar[n])) {
        norm = vvar[n] - sump1;
	norm2 = norm * norm;
	sum2 += norm2;
	sum3 += norm * norm2;
	sum4 += norm2 * norm2;
      }
    }
    kstat += check_and_set_key_double (rec, "PowrRMS", sqrt (sum2 / sump0));
    if (sump0 > 1) {
      var = sum2 / (sump0 - 1);
      sig = sqrt (var);
      kstat += check_and_set_key_double (rec, "PowrSKEW",
	  scale * sum3 / (var * sig));
      kstat += check_and_set_key_double (rec, "PowrKURT",
	  scale * sum4 / (var * var) - 3.0);
    }
  }

  if (verbose) {
    int rank = mean->naxis;
    int rmnd;

    printf ("Extrema: minimum value = %.4e @ %d [", mnvmin, vminloc);
    rmnd = vminloc;
    for (n = 0; n < rank; n++) {
      if (n) printf (", ");
      printf ("%d", rmnd %  mean->axis[n]);
      rmnd /= mean->axis[n];
    }
    printf ("]\n");

    printf ("         maximum value = %.4e @ %d [", mnvmax, vmaxloc);
    rmnd = vmaxloc;
    for (n = 0; n < rank; n++) {
      if (n) printf (", ");
      printf ("%d", rmnd %  mean->axis[n]);
      rmnd /= mean->axis[n];
    }
    printf ("]\n");

    if (vvar) {
      printf ("         minimum power = %.4e @ %d [", pwrmin, pminloc);
      rmnd = pminloc;
      for (n = 0; n < rank; n++) {
	if (n) printf (", ");
	printf ("%d", rmnd %  mean->axis[n]);
	rmnd /= mean->axis[n];
      }
      printf ("]\n");

      printf ("         maximum power = %.4e @ %d [", pwrmax, pmaxloc);
      rmnd = pmaxloc;
      for (n = 0; n < rank; n++) {
	if (n) printf (", ");
	printf ("%d", rmnd %  mean->axis[n]);
	rmnd /= mean->axis[n];
      }
      printf ("]\n");
    }
  }

  if (mscaled || pscaled) {
    morct = morloct + morhict;
    porct = porloct + porhict;
    if (morct) {
      fprintf (stderr,
	  "Warning: %d average values out of representable range\n", morct);
      fprintf (stderr,
	  "         for BScale = %g and BZero = %g\n", mean->bscale,
	  mean->bzero);
      fprintf (stderr,
	  "         %d values < %g, %d values > %g\n", morloct, mrepmin,
	  morhict, mrepmax);
    }
    if (porct) {
      fprintf (stderr,
	  "Warning: %d variance values out of representable range\n", porct);
      fprintf (stderr,
	  "         for PScale = %g and PZero = %g\n", powr->bscale,
	  powr->bzero);
      fprintf (stderr,
	  "         %d values < %g, %d values > %g\n", porloct, prepmin,
	  porhict, prepmax);
    }
  }

  return kstat;
}

void report_pkey_value (FILE *out, DRMS_Record_t *rec, char *key, int recnum) {
  DRMS_Keyword_t *keyword;
  char buf[128];

  fprintf (out, "record ");
  keyword = drms_keyword_lookup (rec, key, 1);
  if (!keyword) {
    fprintf (out, "%d (unknown value) [#%lld]: ", recnum, rec->recnum);
    return;
  }
  drms_keyword_snprintfval (keyword, buf, sizeof (buf));
  fprintf (out, "%s (%d) [#%lld]: ", buf, recnum, rec->recnum);
  fflush (out);
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids = NULL;
  DRMS_Record_t *irec, *orec = NULL;
  DRMS_Segment_t *iseg, *vseg = NULL, *mseg = NULL, *pseg = NULL;
  DRMS_Segment_t *logseg = NULL;
  DRMS_Array_t *data_array, *vcts, *mean, *powr;
  DRMS_Keyword_t *keywd;
  TIME tmid, tobs_rec, /* tobs, */ tfirst, tlast;
  FILE *runlog;
  double *v, *vavg, *vvar;
  double *crpix, *crval, *cdelt, *crota, *avgval;
  double *crpixv, *crvalv, *cdeltv, *crotav, *avgvalv;
  double crpix_rec, crval_rec, cdelt_rec, crota_rec, avg_rec;
  double /* tobsv, */ vobs;
  double log_base;
  double mrepmin, mrepmax, prepmin, prepmax;
  double crlnobs, crlnval, crlnvalv;
  float clstrt, clmid, clstop;
  float pa_rec, dpa;
  long long *cvgolist, *cvnolist, *cvlist;
  long long ntot, calver;
  unsigned int quality;
  int *inaxis, *vval, *reject_list, *cvfound;
  int rec, recct, segct, segnum, maxct, imgct;
  int crstrt, crmid, crstop, carrot;
  int checkseg, mscaled, pscaled;
  int kstat, status;
  int log_status;
  int writelog;
  int n, naxis, wcsaxes;
  int propct, meanct, add_defaults, cvct, cvgoct, cvnoct, cvused, cvmaxct;
  char **copykeylist, **meankeylist;
  char *source, *keystr;
  char recset_query[DRMS_MAXQUERYLEN], keyname[DRMS_MAXKEYNAMELEN];
  char calverinfo[DRMS_MAXQUERYLEN];
  char logfilename[DRMS_MAXPATHLEN];
  char module_ident[64], tbuf[64], vbuf[128], strbuf[128];

  double fp_nan = 0.0 / 0.0;
  int badqual = 0, blacklist = 0, badpa = 0, badcv = 0, rejects = 0;
  int needpowr = 1;
  int propkeyct = sizeof (propagate) / sizeof (char *);
  int meankeyct = sizeof (average) / sizeof (char *);

  char *inset = strdup (params_get_str (params, "in"));
  char *out_series = strdup (params_get_str (params, "out"));
  char *vsegname = strdup (params_get_str (params, "count"));
  char *msegname = strdup (params_get_str (params, "mean"));
  char *psegname = strdup (params_get_str (params, "power"));
  char *logsegname = strdup (params_get_str (params, "log"));
  char *tmid_str = strdup (params_get_str (params, "tmid"));
  float intrvl = params_get_float (params, "length");
  float pa_nom = params_get_float (params, "pa");
  float dpa_max = params_get_float (params, "dpa");
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  long long cvaccept = params_get_mask (params, "cvok", -1);
  long long cvreject = params_get_mask (params, "cvno", 0);
  char *calverkey = strdup (params_get_str (params, "cvkey"));
  char *rejectfile = strdup (params_get_str (params, "reject"));
  char *propagate_req = strdup (params_get_str (params, "copy"));
  char *average_req = strdup (params_get_str (params, "average"));
  char *primekey = strdup (params_get_str (params, "pkey"));
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *roll_key = strdup (params_get_str (params, "roll_key"));
  char *spec_key = strdup (params_get_str (params, "setkey"));
  double spec_val = params_get_double (params, "setval");
  char *specuse_key = strdup (params_get_str (params, "setval"));
  double mscale = params_get_double (params, "mscale");
  double mzero = params_get_double (params, "mzero");
  double pscale = params_get_double (params, "pscale");
  double pzero = params_get_double (params, "pzero");
  int no_save = params_isflagset (params, "n");
  int remove_obsvel = params_isflagset (params, "o");
  int verbose = params_isflagset (params, "v");

  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
  int check_pa = (dpa_max < 180.0);
  int set_extra_key = (strcmp (spec_key, "Not Specified")) ? 1 : 0;
  int use_other_key = isnan (spec_val);
  int filt_on_calver = (cvreject || ~cvaccept);

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);
			 /*  initialize CalVer lists (single element for now  */
  cvused = 0;
  cvmaxct = 64;
  cvfound = (int *)malloc (cvmaxct * sizeof (int));
  cvlist = (long long *)malloc (cvmaxct * sizeof (long long));
  if (filt_on_calver) {
    cvgoct = 1;
    cvgolist = (long long *)malloc (cvgoct * sizeof (long long));
    cvgolist[0] = cvaccept;
    cvnoct = 1;
    cvnolist = (long long *)malloc (cvnoct * sizeof (long long));
    cvnolist[0] = cvreject;
    if (~cvaccept) cvgoct = 0;
    if (!cvreject) cvnoct = 0;
  }
  if (!isfinite (pa_nom)) pa_nom = 180.0;
					      /*  additional initializations  */
  prepmin = mrepmin = 1.0 / 0.0;
  prepmax = mrepmax = - mrepmin;
						    /*  create output record  */
  if (no_save) {
    orec = drms_create_record (drms_env, out_series, DRMS_TRANSIENT, &status);
    if (status) {
      fprintf (stderr, "Warning: drms_create_record returned %d for data series:\n",
	  status);
      fprintf (stderr, "  %s\nKeyword setting will not be checked\n", out_series);
    }
  } else {
    orec = drms_create_record (drms_env, out_series, DRMS_PERMANENT, &status);
    if (status) {
      fprintf (stderr, "Error: drms_create_record returned %d for data series:\n",
	  status);
      fprintf (stderr, "  %s\n", out_series);
      return 1;
    }
  }
  if ((segct = orec->segments.num_total) < 1) {
    fprintf (stderr, "Error: no data segments in output data series:\n");
    fprintf (stderr, "  %s\n", out_series);
    return 1;
  }
  vseg = drms_segment_lookup (orec, vsegname);
  mseg = drms_segment_lookup (orec, msegname);
  pseg = drms_segment_lookup (orec, psegname);
  if (!vseg && !mseg && !pseg) {
    fprintf (stderr,
	"Error: output series %s does not contain any of the required segments:\n",
	out_series);
    fprintf (stderr, "       \"%s\", \"%s\", \"%s\"\n", vsegname, msegname,
	psegname);
    drms_close_record (orec, DRMS_FREE_RECORD);
    return 1;
  }
  if (vseg) {
    switch (vseg->info->type) {
      case (DRMS_TYPE_CHAR):
	maxct = SCHAR_MAX;
	break;
      case (DRMS_TYPE_SHORT):
	maxct = SHRT_MAX;
	break;
      default:
	maxct = INT_MAX;
    }
    maxct += vseg->bzero;
    maxct /= vseg->bscale;
  } else {
    fprintf (stderr,
	"Warning: output series %s does not contain the segment %s\n",
	out_series, vsegname);
    maxct = INT_MAX;
  }
  if (mseg) {
    mscaled = (mseg->info->type == DRMS_TYPE_CHAR) ||
              (mseg->info->type == DRMS_TYPE_SHORT) ||
              (mseg->info->type == DRMS_TYPE_INT) ||
              (mseg->info->type == DRMS_TYPE_LONGLONG);
  } else {
    fprintf (stderr,
	"Warning: output series %s does not contain the segment %s\n",
	out_series, msegname);
    mscaled = 0;
  }
  if (pseg) {
    pscaled = (pseg->info->type == DRMS_TYPE_CHAR) ||
              (pseg->info->type == DRMS_TYPE_SHORT) ||
              (pseg->info->type == DRMS_TYPE_INT) ||
              (pseg->info->type == DRMS_TYPE_LONGLONG);
  } else {
    fprintf (stderr,
	"Warning: output series %s does not contain the segment %s\n",
	out_series, psegname);
    needpowr = pscaled = 0;
  }
  writelog = 0;
  logseg = drms_segment_lookup (orec, logsegname);
  if (logseg && dispose == DRMS_INSERT_RECORD) {
    drms_segment_filename (logseg, logfilename);
    runlog = fopen (logfilename, "w");
    if (runlog) writelog = 1;
  }

					     /*  process input specification  */
  if (key_params_from_dspec (inset)) {
				    /*  input specified as specific data set  */
    if (!(ids = drms_open_records (drms_env, inset, &status))) {
      fprintf (stderr, "Error: (%s) unable to open input data set %s\n",
        module_ident, inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    if ((recct = ids->n) < 2) {
      fprintf (stderr, "Error: (%s) <2 records in selected input set\n",
	  module_ident);
      fprintf (stderr, "       %s\n", inset);
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    source = strdup (inset);
    for (n = 0; n < strlen (source); n++) if (source[n] == '[') {
      source[n] = '\0';
      break;
    }
    tmid = clmid = fp_nan;
    crmid = -1;
  } else {
				/*  only the input data series is named,
				    get record specifications from arguments  */
    source = strdup (inset);
    if (sscanf (tmid_str, "%d:%f", &crmid, &clmid) == 2) {
						 /*  tmid specified as CR:CL  */
      if (isnan (clmid)) {
 	fprintf (stderr,
	    "Error: target Carrington longitude = %f\n", clmid);
	return 1;
      }
      while (clmid <= 0) {
        clmid += 360;
	crmid++;
      }
      while (clmid > 360) {
        clmid -= 360;
	crmid--;
      }
      clstrt = clmid + 0.5 * intrvl;
      clstop = clmid - 0.5 * intrvl;
      crstrt = crstop = crmid;
      if (intrvl < 360.0 && clstrt <= 360.0 && clstop >= 0.0) {
	sprintf (recset_query,
	    "%s[?Car_Rot=%d and CRLN_OBS<=%.2f and CRLN_OBS>=%.2f?]", source,
	    crmid, clstrt, clstop);
      } else {
	while (clstrt > 360.0) {
	  clstrt -= 360.0;
	  crstrt--;
	}
	while (clstop < 0.0) {
	  clstop += 360.0;
	  crstop++;
	}
	if ((crstop - crstrt) > 1)
	  sprintf (recset_query,
	      "%s[?(Car_Rot=%d and CRLN_OBS<=%.2f) or (Car_Rot=%d and CRLN_OBS>=%.2f) or (CarRot>%d and CarRot<%d?]",
	      source, crstrt, clstrt, crstop, clstop, crstrt, crstop);
	else
	  sprintf (recset_query,
	      "%s[?(Car_Rot=%d and CRLN_OBS<=%.2f) or (Car_Rot=%d and CRLN_OBS>=%.2f)?]",
	      source, crstrt, clstrt, crstop, clstop);
      }
      if (!(ids = drms_open_records (drms_env, recset_query, &status))) {
	fprintf (stderr, "Error: (%s) unable to open input data set %s\n",
          module_ident, recset_query);
	fprintf (stderr, "       status = %d\n", status);
	return 1;
      }
      if ((recct = ids->n) < 2) {
	fprintf (stderr, "Error: (%s) <2 records in selected input set\n",
	    module_ident);
	fprintf (stderr, "       %s\n", recset_query);
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
      free (inset);
      inset = strdup (recset_query);
      tmid = earth_meridian_crossing (clmid, crmid);
    } else {
			      /*  tmid specified as normal date-time string  */
      tmid = sscan_time (tmid_str);
      clmid = fp_nan;
      crmid = -1;
       sprintf (recset_query, "%s[?%s between %f and %f?]", source,
	  primekey, tmid - 0.5 * intrvl, tmid + 0.5 * intrvl);
      if (!(ids = drms_open_records (drms_env, recset_query, &status))) {
	fprintf (stderr, "Error: (%s) unable to open input data set %s\n",
          module_ident, recset_query);
	fprintf (stderr, "       status = %d\n", status);
	return 1;
      }
      if ((recct = ids->n) < 2) {
	fprintf (stderr, "Error: (%s) <2 records in selected input set\n",
	    module_ident);
	fprintf (stderr, "       %s\n", recset_query);
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
      free (inset);
      inset = strdup (recset_query);
   }
  }
  if (verbose) printf ("processing %d input records\n", recct);
  propct = construct_stringlist (propagate_req, ',', &copykeylist);
					       /*  replace '+' with defaults  */
  if (propkeyct) {
    add_defaults = 0;
    for (n = 0; n < propct; n++) {
      if (!strcmp (copykeylist[n], "+")) {
	add_defaults = 1;
	strcpy (copykeylist[n], propagate[0]);
	n = propct;
      }
    }
    if (add_defaults) {
      int newct = propct + propkeyct - 1;
      copykeylist = (char **)realloc (copykeylist, newct * sizeof (char **));
      for (n = 1; n < propkeyct; n++) {
	copykeylist[n + propct - 1] = propagate[n];
      }
      propct = newct;
    }
  }
  if (verbose && propct) {
    printf ("propagating values for up to %d key(s):\n", propct);
    for (n = 0; n < propct; n++) {
      if (drms_keyword_lookup (orec, copykeylist[n], 0))
	printf ("  %s\n", copykeylist[n]);
      else
	printf ("  %s (not in output series)\n", copykeylist[n]);
    }
  }

  meanct = construct_stringlist (average_req, ',', &meankeylist);
					       /*  replace '+' with defaults  */
  if (meankeyct) {
    add_defaults = 0;
    for (n = 0; n < meanct; n++) {
      if (!strcmp (meankeylist[n], "+")) {
	add_defaults = 1;
	strcpy (meankeylist[n], average[0]);
	n = meanct;
      }
    }
    if (add_defaults) {
      int newct = meanct + meankeyct - 1;
      meankeylist = (char **)realloc (meankeylist, newct * sizeof (char **));
      for (n = 1; n < meankeyct; n++)
	meankeylist[n + meanct - 1] = average[n];
      meanct = newct;
    }
  }
  if (verbose) {
    printf ("averaging values for up to %d key(s):\n", meanct);
    for (n = 0; n < meanct; n++) {
      if (drms_keyword_lookup (orec, meankeylist[n], 0)) {
	printf ("  %s\n", meankeylist[n]);
	sprintf (keyname, "D_%s", meankeylist[n]);
	if (!drms_keyword_lookup (orec, keyname, 0))
	  printf ("    Warning: %s not in output series\n", keyname);
      } else
	printf ("  %s (not in output series)\n", meankeylist[n]);
    }
  }
  avgval = (double *)calloc (meanct, sizeof (double));
  avgvalv = (double *)calloc (meanct, sizeof (double));

  irec = ids->records[0];
						 /*  establish input segment  */
  segct = drms_record_numsegments (irec);
  segnum = 0;
  if (segct > 1) {
    fprintf (stderr,
	"input records contain multiple segments, segment must be specified\n");
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 0;
  }
  iseg = drms_segment_lookupnum (irec, segnum);
  naxis = iseg->info->naxis;
		   /*  will dimensions of input segments have to be checked?  */
  checkseg = 0;
  if (iseg->info->scope == DRMS_VARDIM) {
    if (verbose) printf ("Warning: input segment sizes may vary\n");
    checkseg = 1;
    inaxis = (int *)malloc (naxis * sizeof (int));
  }
  ntot = 1;
  for (n = 0; n < naxis; n++) {
    ntot *= iseg->axis[n];
    if (checkseg) inaxis[n] = iseg->axis[n];
  }
					     /*  check for required keywords  */
  if (remove_obsvel) {
    keywd = drms_keyword_lookup (irec, "OBS_VR", 1);
    if (!keywd) {
      fprintf (stderr, "Warning: required keyword %s not found\n", "OBS_VR");
      fprintf (stderr, "         no correction applied for observer velocity\n");
      remove_obsvel = 0;
    }
  }
  if (filt_on_calver) {
    keywd = drms_keyword_lookup (irec, calverkey, 1);
    if (!keywd) {
      fprintf (stderr, "Warning: required keyword %s not found\n", calverkey);
      fprintf (stderr, "         no calibration version filtering applied\n");
      filt_on_calver = 0;
    }
  }
					   /*  set up WCS keys for averaging  */
  wcsaxes = drms_getkey_int (irec, "WCSAXES", &status);
  if (status) wcsaxes = naxis;
  crpix = (double *)calloc (wcsaxes, sizeof (double));
  crval = (double *)calloc (wcsaxes, sizeof (double));
  cdelt = (double *)calloc (wcsaxes, sizeof (double));
  crota = (double *)calloc (wcsaxes, sizeof (double));
  crpixv = (double *)calloc (wcsaxes, sizeof (double));
  crvalv = (double *)calloc (wcsaxes, sizeof (double));
  cdeltv = (double *)calloc (wcsaxes, sizeof (double));
  crotav = (double *)calloc (wcsaxes, sizeof (double));
                       /*  create new data maps for each of three statistics  */
  vval = (int *)calloc (ntot, sizeof (int));
  vavg = (double *)calloc (ntot, sizeof (double));
  if (needpowr) {
    vvar = (double *)calloc (ntot, sizeof (double));
    if (!vvar) {
      fprintf (stderr,
	  "Error: unable to allocate %lld bytes required for output arrays\n",
	  ntot * (sizeof (int) + 2 * sizeof (double)));
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
  }
  if (!vval || !vavg) {
    fprintf (stderr,
	"Error: unable to allocate %lld bytes required for output arrays\n",
	ntot * (sizeof (int) + sizeof (double)));
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
  vcts = drms_array_create (DRMS_TYPE_INT, naxis, iseg->axis, (void *)vval,
      &status);
  mean = drms_array_create (DRMS_TYPE_DOUBLE, naxis, iseg->axis, (void *)vavg,
      &status);
  if (needpowr) powr = drms_array_create (DRMS_TYPE_DOUBLE, naxis, iseg->axis,
      (void *)vvar, &status);
		/*  set the output scaling factors to be segment defaults
				    if not overridden as argument parameters  */
  if (mseg) {
    mean->bscale = (isnan (mscale) || mscale == 0.0) ? mseg->bscale : mscale;
    mean->bzero = (isnan (mzero)) ?  mseg->bzero : mzero;
  } else {
    mean->bscale = (isnan (mscale) || mscale == 0.0) ? 1.0 : mscale;
    mean->bzero = (isnan (mzero)) ?  0.0 : mzero;
  }
  if (needpowr) {
    if (pseg) {
      powr->bscale = (isnan (pscale) || pscale == 0.0) ? pseg->bscale : pscale;
      powr->bzero = (isnan (pzero)) ?  pseg->bzero : pzero;;
    } else {
      powr->bscale = (isnan (pscale) || pscale == 0.0) ? 1.0 : pscale;
      powr->bzero = (isnan (pzero)) ?  0.0 : pzero;;
    }
  }
	       /*  determine range of representable values for scaled output  */
  if (mscaled) {
    unsigned long long maxval =
	(mseg->info->type == DRMS_TYPE_CHAR) ? SCHAR_MAX :
	(mseg->info->type == DRMS_TYPE_SHORT) ? SHRT_MAX :
	(mseg->info->type == DRMS_TYPE_INT) ? INT_MAX : LLONG_MAX;
    mrepmax = maxval;
    mrepmin = -mrepmax;
    mrepmax *= mean->bscale;
    mrepmin *= mean->bscale;
    mrepmax += mean->bzero;
    mrepmin += mean->bzero;
  }
  if (pscaled) {
    unsigned long long maxval =
	(pseg->info->type == DRMS_TYPE_CHAR) ? SCHAR_MAX :
	(pseg->info->type == DRMS_TYPE_SHORT) ? SHRT_MAX :
	(pseg->info->type == DRMS_TYPE_INT) ? INT_MAX : LLONG_MAX;
    prepmax = maxval;
    prepmin = -prepmax;
    prepmax *= powr->bscale;
    prepmin *= powr->bscale;
    prepmax += powr->bzero;
    prepmin += powr->bzero;
  }

  if (recct > maxct) {
    double bscale = 1.0;
    if (vseg) bscale = vseg->bscale;
    while (recct > maxct) {
      bscale *= 0.5;
      maxct *= 2;
    }
    vcts->bscale = bscale;
    fprintf (stderr,
        "Warning: more records in input data set than maximum valid count\n");
    if (vseg) fprintf (stderr,
        "BSCALE adjusted from %g to %g\n", vseg->bscale, bscale);
  }
	 /*  write out some keys to make sure they're okay before proceeding  */
  kstat = 0;
  kstat += check_and_set_key_int   (orec, "CarrRot", crmid);
  kstat += check_and_set_key_float (orec, "CMLon", clmid);
  kstat += check_and_set_key_time  (orec, "MidTime", tmid);
  kstat += check_and_set_key_str   (orec, "Module", module_ident);
  kstat += check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
  kstat += check_and_set_key_str   (orec, "Source", source);
  kstat += check_and_set_key_str   (orec, "Input", inset);
  kstat += check_and_set_key_time  (orec, "Created", CURRENT_SYSTEM_TIME);
  kstat += check_and_set_key_float (orec, "Interval", intrvl);
  kstat += propagate_keys (orec, irec, copykeylist, propct);
					       /*  propagate WCS string keys  */
  for (n = 0; n < wcsaxes; n++) {
    sprintf (keyname, "CTYPE%d", n + 1);
    keystr = drms_getkey_string (irec, keyname, &status);
    if (!status) kstat += check_and_set_key_str (orec, keyname, keystr);
    sprintf (keyname, "CUNIT%d", n + 1);
    keystr = drms_getkey_string (irec, keyname, &status);
    if (!status) kstat += check_and_set_key_str (orec, keyname, keystr);
  }
  if (kstat) {
    fprintf (stderr, "Error writing key value(s) to %s\n", out_series);
    fprintf (stderr, "      output series may not have appropriate structure\n");
    if (!no_save) return 1;
  }
		 	   /*  write argument list to Comment (or History)  */
/*
  if (filt_on_calver) {
    calverinfo[0] = '\0';
    if (~cvaccept) {
      sprintf (calverinfo, "%s values accepted: %016llx", calverkey, cvaccept);
      if (cvreject) strcat (calverinfo, "\n");
    }
    if (cvreject) {
      sprintf (strbuf, "%s values rejected: %016llx", calverkey, cvreject);
      strcat (calverinfo, strbuf);
    }
  } else
    sprintf (calverinfo, "%s values accepted: ANY", calverkey);
*/
  if ((keywd = drms_keyword_lookup (orec, "COMMENT", 1))) {
    append_args_tokey (orec, "COMMENT");
/*
    drms_appendhistory (orec, calverinfo, 1);
*/
  } else if ((keywd = drms_keyword_lookup (orec, "HISTORY", 1))) {
    append_args_tokey (orec, "HISTORY");
/*
    drms_appendhistory (orec, calverinfo, 1);
*/
  }
		  /*  support special hack of reading of rejection list file  */
  rejects = 0;
  if (strcmp (rejectfile, "Not Specified")) {
    FILE *rejectfp = fopen (rejectfile, "r");
    if (rejectfp) rejects = read_reject_list (rejectfp, &reject_list);
    else fprintf (stderr,
	"Warning: could not open rejection list %s; ignored\n", rejectfile);
  }
					  /*  initialize keys to be averaged  */
/*
  tobs = tobsv = 0;
*/
  imgct = 0;
  crlnval = crlnvalv = 0.0;
					      /*  loop through input records  */
  for (rec = 0; rec < recct; rec++) {
    irec = ids->records[rec];
    if (verbose) report_pkey_value (stdout, irec, primekey, rec);
    if (writelog) report_pkey_value (runlog, irec, primekey, rec);

			  /*  check for record quality, reject as applicable  */
    quality = drms_getkey_int (irec, qual_key, &status);
    if ((quality & qmask) && !status) {
      badqual++;
      if (verbose) printf ("skipped (quality)\n");
      if (writelog) fprintf (runlog, "skipped (quality)\n");
      continue;
    }
    if (rejects) {
				      /*  check against special rection list  */
      int idrec = drms_getkey_int (irec, "T_REC_index", &status);
      int match = 0;
      if (status) {
	fprintf (stderr, "Warning: \"T_REC_index\" keyword not found\n");
	fprintf (stderr, "         up to %d bad images could be processed\n",
	  rejects);
	rejects = 0;
      }
      for (n = 0; n < rejects; n++) {
	if (idrec == reject_list[n]) {
	  match = 1;
	  break;
	}
      }
      if (match) {
        blacklist++;
        if (verbose) printf ("skipped (blacklist)\n");
        if (writelog) fprintf (runlog, "skipped (blacklist)\n");
        continue;
      }
    }
	      /*  check for record calibration version; reject as applicable  */
    calver = drms_getkey_longlong (irec, calverkey, &status);
    if (filt_on_calver) {
      if (status) {
		     /*  this probably never happens if the keyword exists  */
	if (verbose) printf ("skipped\n  (no %s key value)\n", calverkey);
	if (writelog) fprintf (runlog, "skipped\n  (no %s key value)\n",
	    calverkey);
	continue;
      }
      for (cvct = 0; cvct < cvgoct; cvct++)
	if (calver == cvgolist[cvct]) break;
      if (cvct >= cvgoct) {
	badcv++;
	if (verbose)
	  printf ("skipped\n  (calibration version key %s = %016llx)\n",
	      calverkey, calver);
	if (writelog)
	  fprintf (runlog, "skipped\n  (calibration version key %s = %016llx)\n",
	      calverkey, calver);
	continue;
      }
      for (cvct = 0; cvct < cvnoct; cvct++) {
 	if (calver == cvnolist[cvct]) {
	  badcv++;
	  if (verbose)
	    printf ("skipped\n  (calibration version key %s = %016llx)\n",
	        calverkey, calver);
	  if (writelog)
	    fprintf (runlog, "skipped\n  (calibration version key %s = %016llx)\n",
	        calverkey, calver);
	  continue;
	}
      }
    }
    for (cvct = 0; cvct < cvused; cvct++) {
      if (calver == cvlist[cvct]) {
        cvfound[cvct]++;
        break;
      }
    }
    if (cvct == cvused) {
      if (cvused < cvmaxct) {
        cvlist[cvct] = calver;
	cvfound[cvct] = 1;
	cvused++;
      }
    }
			      /*  check for non-nominal image rotation angle  */
    if (check_pa) {
      pa_rec = drms_getkey_float (irec, roll_key, &status);
      if (status && check_pa) {
	fprintf (stderr, "Warning: \"%s\" keyword not found\n", roll_key);
	fprintf (stderr, "         no limits on rotation angle\n");
	check_pa = 0;
	dpa_max = 360.0;
      }
      if (isfinite (pa_rec)) {
	dpa = fabs (pa_rec - pa_nom);
	while (dpa > 180.0) dpa -= 360.0;
	while (dpa < 0.0) dpa += 360.0;
	if (dpa > dpa_max) {
	  badpa++;
	  if (verbose) printf ("skipped (rotated)\n");
	  if (writelog) fprintf (runlog, "skipped (rotated)\n");
	  continue;
	}
      }
    }
					         /*  check segment structure  */
    iseg = drms_segment_lookupnum (irec, segnum);
    if (checkseg) {
      int okay = 1;
      for (n = 0; n < naxis; n++) {
	if (inaxis[n] != iseg->axis[n]) {
	  okay = 0;
	  break;
	}
      }
      if (!okay) {
        if (verbose) printf ("skipped (dimension mismatch)\n");
        if (writelog) fprintf (runlog, "skipped (dimension mismatch)\n");
        continue;
      }
    }
    						   /*  read input data image  */
    data_array = drms_segment_read (iseg, DRMS_TYPE_DOUBLE, &status);
    if (status) {
      if (data_array) drms_free_array (data_array);
      if (verbose) printf ("skipped (segment read error)\n");
      if (writelog) fprintf (runlog, "skipped (segment read error)\n");
      else {
	fprintf (stderr, "Error reading file for record %d: ", rec);
	keywd = drms_keyword_lookup (irec, primekey, 1);
	if (!keywd) fprintf (stderr, "unknown value\n");
	else {
	  drms_keyword_snprintfval (keywd, vbuf, sizeof (vbuf));
	  fprintf (stderr, "%s\n", vbuf);
	}
      }
      continue;
    }

    tobs_rec = drms_getkey_time (irec, primekey, &status);
    if (status || time_is_invalid (tobs_rec)) {
      if (writelog) fprintf (runlog, "skipped (time invalid)\n");
      if (verbose) printf ("skipped (time invalid)\n");
      else fprintf (stderr, "error reading %s from record #%d\n", primekey, rec);
      if (data_array) drms_free_array (data_array);
      continue;
    }
    if (!imgct) tfirst = tobs_rec;
    tlast = tobs_rec;
/*
    tobs += tobs_rec;
    tobsv += tobs_rec * tobs_rec;
*/
    imgct++;

    v = (double *)data_array->data;
    vval = (int *)vcts->data;
    vavg = (double *)mean->data;
    if (needpowr) vvar = (double *)powr->data;
    vobs = (remove_obsvel) ? drms_getkey_double (irec, "OBS_VR", &status) : 0.0;
				     /*  exponentiate logarithm as necessary  */
    log_base = drms_getkey_double (irec, "LOG_BASE", &status);
    if (!status && isfinite (log_base)) {
      log_base = log (log_base);
      for (n = 0; n < ntot; n++) v[n] = exp (log_base * v[n]);
      log_status = 1;
    } else log_status = 0;
						      /*  process valid data  */
    n = ntot;
    if (needpowr) {
      while (n--) {
	if (isfinite (*v)) {
	  (*vval)++;
	  *v -= vobs;
/*
	if (function == DO_ABSVAL) *v = fabs (*v);
	if (function == DO_SQUARE) *v *= *v;
*/
          *vavg += *v;
          *vvar += *v * *v;
        }
	v++;
	vval++;
	vavg++;
	vvar++;
      }
    } else {
      while (n--) {
	if (isfinite (*v)) {
	  (*vval)++;
	  *v -= vobs;
/*
	if (function == DO_ABSVAL) *v = fabs (*v);
	if (function == DO_SQUARE) *v *= *v;
*/
          *vavg += *v;
        }
	v++;
	vval++;
	vavg++;
      }
    }
    drms_free_array (data_array);
					      /*  get WCS keys for averaging  */
    for (n = 0; n < wcsaxes; n++) {
      sprintf (keyname, "CRPIX%d", n + 1);
      crpix_rec = drms_getkey_double (irec, keyname, &status);
      if (!status && isfinite (crpix_rec)) {
	crpix[n] += crpix_rec;
	crpixv[n] += crpix_rec * crpix_rec;
      }
      sprintf (keyname, "CRVAL%d", n + 1);
      crval_rec = drms_getkey_double (irec, keyname, &status);
      if (!status && isfinite (crval_rec)) {
	crval[n] += crval_rec;
	crvalv[n] += crval_rec * crval_rec;
      }
      sprintf (keyname, "CDELT%d", n + 1);
      cdelt_rec = drms_getkey_double (irec, keyname, &status);
      if (!status && isfinite (cdelt_rec)) {
	cdelt[n] += cdelt_rec;
	cdeltv[n] += cdelt_rec * cdelt_rec;
      }
      sprintf (keyname, "CROTA%d", n + 1);
      crota_rec = drms_getkey_double (irec, keyname, &status);
      if (!status && isfinite (crota_rec)) {
        crota_rec -= pa_nom;
	while (crota_rec < -180.0) crota_rec += 360.0;
	while (crota_rec > 180.0) crota_rec -= 360.0;
	crota[n] += crota_rec;
	crotav[n] += crota_rec * crota_rec;
      }
    }
				   /*  averaged CRLN_OBS and CarrRot/CAR_ROT  */
    crlnobs = drms_getkey_double (irec, "CRLN_OBS", &status);
    if (!status && isfinite (crlnobs)) {
      carrot = drms_getkey_int (irec, "CAR_ROT", &status);
      if (status || (carrot < 1)) {
	carrot = drms_getkey_int (irec, "CarrRot", &status);
	if (status || (carrot < 1)) carrot = 0;
      }
      crlnobs = 360.0 * carrot - crlnobs;
      crlnval += crlnobs;
      crlnvalv += crlnobs * crlnobs;
    }
				        /*  get requested keys for averaging  */
    for (n = 0; n < meanct; n++) {
      avg_rec = drms_getkey_double (irec, meankeylist[n], &status);
      if (!status && isfinite (avg_rec)) {
	avgval[n] += avg_rec;
	avgvalv[n] += avg_rec * avg_rec;
      }
    }
    if (verbose) printf ("ok\n");
    if (writelog) fprintf (runlog, "ok\n");
  }
  if (checkseg) free (inaxis);
  vval = (int *)vcts->data;
  vavg = (double *)mean->data;
  if (needpowr) {
    vvar = (double *)powr->data;
    for (n = 0; n < ntot; n++) {
      if (vval[n]) {
	vavg[n] /= vval[n];
	vvar[n] /= vval[n];
	vvar[n] -= vavg[n] * vavg[n];
      } else vavg[n] = vvar[n] = fp_nan;
    }
  } else {
    for (n = 0; n < ntot; n++) {
      if (vval[n]) vavg[n] /= vval[n];
      else vavg[n] = fp_nan;
    }
  }
  if (log_status) {
    double scale = 1.0 / log_base;
    if (needpowr) {
      for (n = 0; n < ntot; n++) {
	vavg[n] = scale * log (vavg[n]);
	vvar[n] = scale * log (vvar[n]);
      }
    } else {
      for (n = 0; n < ntot; n++) vavg[n] = scale * log (vavg[n]);
    }
  }
  if (vseg)
    if (drms_segment_write (vseg, vcts, 0))
      fprintf (stderr, "Warning: unable to write to count segment\n");
  if (mseg)
    if (drms_segment_write (mseg, mean, 0))
      fprintf (stderr, "Warning: unable to write to mean segment\n");
  if (pseg)
    if (drms_segment_write (pseg, powr, 0))
      fprintf (stderr, "Warning: unable to write to variance segment\n");
						      /*  set remaining keys  */
  kstat = 0;
  kstat += check_and_set_key_int  (orec, "DataRecs", imgct);
  kstat += check_and_set_key_int  (orec, "MissRecs", recct - imgct);
							 /*  statistics keys  */
  kstat += set_stats_keys (orec, vcts, mean, powr, ntot, mscaled, pscaled,
      mrepmin, mrepmax, prepmin, prepmax, verbose);
						       /*  averaged WCS keys  */
  if (imgct) {
    kstat += check_and_set_key_time  (orec, "T_FIRST", tfirst);
    kstat += check_and_set_key_time  (orec, "T_LAST", tlast);
/*
    tobs /= imgct;
    kstat += check_and_set_key_time  (orec, primekey, tobs);
    tobsv /= imgct;
    tobsv -= tobs * tobs;
    sprintf (keyname, "D_%s", primekey);
    kstat += check_and_set_key_time  (orec, keyname, sqrt (tobsv));
*/
    for (n = 0; n < wcsaxes; n++) {
      crpix[n] /= imgct;
      sprintf (keyname, "CRPIX%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, crpix[n]);
      crpixv[n] /= imgct;
      crpixv[n] -= crpix[n] * crpix[n];
      sprintf (keyname, "D_CRPIX%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, sqrt (crpixv[n]));
      crval[n] /= imgct;
      sprintf (keyname, "CRVAL%d", n + 1);
      kstat += check_and_set_key_double  (orec, keyname, crval[n]);
      crvalv[n] /= imgct;
      crvalv[n] -= crval[n] * crval[n];
      sprintf (keyname, "D_CRVAL%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, sqrt (crvalv[n]));
      cdelt[n] /= imgct;
      sprintf (keyname, "CDELT%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, cdelt[n]);
      cdeltv[n] /= imgct;
      cdeltv[n] -= cdelt[n] * cdelt[n];
      sprintf (keyname, "D_CDELT%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, sqrt (cdeltv[n]));
      crota[n] /= imgct;
      crotav[n] /= imgct;
      crotav[n] -= crota[n] * crota[n];
      crota[n] += pa_nom;
      sprintf (keyname, "CROTA%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, crota[n]);
      sprintf (keyname, "D_CROTA%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, sqrt (crotav[n]));
    }
				 /*  averaged CRLN_OBS (and CarrRot/CAR_ROT)  */
    crlnval /= imgct;
    crlnobs = 360.0 - fmod (crlnval, 360.0);
    kstat += check_and_set_key_double (orec, "CRLN_OBS", crlnobs);
    crlnvalv /= imgct;
    crlnvalv -= crlnval * crlnval;
    kstat += check_and_set_key_double (orec, "D_CRLN_OBS", sqrt (crlnvalv));
						     /*  other averaged keys  */
    for (n = 0; n < meanct; n++) {
      avgval[n] /= imgct;
      kstat += check_and_set_key_double  (orec,  meankeylist[n], avgval[n]);
      avgvalv[n] /= imgct;
      avgvalv[n] -= avgval[n] * avgval[n];
      sprintf (keyname, "D_%s", meankeylist[n]);
      kstat += check_and_set_key_double (orec, keyname, sqrt (avgvalv[n]));
    }
    if (log_status) {
      sprintf (keyname, "LOG_BASE");
      kstat += check_and_set_key_double (orec, "LOG_BASE", exp (log_base));
    }
  }
  if (set_extra_key) {
    if (use_other_key) {
      spec_val = drms_getkey_double (orec, specuse_key, &status);
      if (status) {
        fprintf (stderr,
	    "Warning: key %s does not exist in output series or is of wront type\n",
	    specuse_key);
	fprintf (stderr, "         key %s will not be set\n", spec_key);
	set_extra_key = 0;
      }
    }
    kstat += check_and_set_key_double (orec, spec_key, spec_val);
  }
  if (kstat) {
    fprintf (stderr, "Error writing key value(s) to %s\n", out_series);
    fprintf (stderr, "      output series may not have appropriate structure\n");
    if (!no_save) return 1;
  }
  							        /*  clean up  */
  drms_close_records (ids, DRMS_FREE_RECORD);
  if (verbose) {
    printf ("record %s[:#%lld] ", out_series, orec->recnum);
    if (dispose == DRMS_FREE_RECORD) printf ("not ");
    printf ("written\n\n");
    if (badqual) printf
	("    %d input records rejected for quality matching %08x\n",
	badqual, qmask);
    if (blacklist) printf
	("    %d input records rejected from rejection list\n", blacklist);
    if (badpa) printf
	("    %d input records rejected for roll difference exceeding %.2f\n",
	badpa, dpa_max);
    if (badcv) {
      printf ("    %d input records rejected for calib version matching %016llx\n",
	  badcv, cvreject);
      printf ("                     or failing to match %016llx\n", cvaccept);
    }
    printf ("%s values used:", calverkey);
    for (cvct = 0; cvct < cvused; cvct++) {
      if (cvct) printf (",");
      printf (" %016llx (%d)", cvlist[cvct], cvfound[cvct]);
    }
    printf ("\n");
  }
  if (writelog) {
    if (badqual) fprintf (runlog,
	"    %d input records rejected for quality matching %08x\n",
	badqual, qmask);
    if (blacklist) fprintf (runlog,
	"    %d input records rejected from rejection list\n", blacklist);
    if (badpa) fprintf (runlog,
	"    %d input records rejected for roll difference exceeding %.2f\n",
	badpa, dpa_max);
    if (badcv) {
      fprintf (runlog,
	  "    %d input records rejected for calib version matching %016llx\n",
	  badcv, cvreject);
      fprintf (runlog,
	  "                     or failing to match %016llx\n", cvaccept);
    }
    fprintf (runlog, "%s values used:", calverkey);
    for (cvct = 0; cvct < cvused; cvct++) {
      if (cvct) fprintf (runlog, ",");
      fprintf (runlog, " %016llx (%d)", cvlist[cvct], cvfound[cvct]);
    }
    fprintf (runlog, "\n");
    fclose (runlog);
  }
  drms_close_record (orec, dispose);
  return 0;
}

/*
 *  Revision history
 *  (all mods by Rick Bogart unless otherwise indicated)
 *
 *  v 0.0  09.11.01	created file, based on version 5.0 of SOI module
 *		doppavg (CM version of 2004.10.25)
 *  v 0.1  09.12.01	include a certain minimal amount of input processing
 *  v 0.2  10.01.12	include minimal processing and output to DRMS
 *  v 0.3  10.01.13	checks on input size
 *  v 0.4  10.02.09	include time selection by target CR:CL and rot width;
 *		setting of a few keys
 *  v 0.5	some cleanup, added setting of a few more keys, quality
 *		checking from keyword and blacklist, stubs for propagation of
 *		fixed keys and averaging of varying keys
 *  v 0.6	Stripped commented code from doppavg; see v05 for them
 *		Added support for key propagation, output keys MidTime,
 *		  averaged, WCS keys, statistics keys
 *		Added optional correction for observer velocity (radial only)
 *		Moved construct_stringlist() function to keystuff
 *		Added trap for bad segment read
 *		Added option for no output
 *		Added diagnostics at end
 *  v 0.6 frozen 10.08.09
 *  v 0.7
 *	10.08.11	Removed default values for in/out arguments and
 *		  unspecified defaults for input set specifers
 *		Fixed bug that was preventing HMI rejection lists from being
 *		  read properly
 *	10.08.18	CRVALn, CRPIXn, CDELTn, & CROTn always averaged with
 *		  standard deviations; CTYPEn and CUNITn always propagated;
 *		  generalized to arbitrary number of WCS axes
 *		Implemented expansion of keyword requests for copied and
 *		  averaged keywords to include defaults if requested
 *		Implemented averaging of requested keywords
 *  v 0.7 frozen 10.08.19
 *  v 0.8
 *	10.09.13	S Chakraborty added code to support log-base data
 *	10.11.16	added test for acceptable roll-angles, with defaults
 *		appropriate to HMI
 *  v 0.8 frozen 10.11.16
 *  v 0.9
 *	11.01.11	added reporting of locations of extremal values, checks
 *		for out-of-range values
 *	11.02.10	added traps for inappropriate output series structure
 *  v 0.9 frozen 11.02.28
 *  v 1.0
 *	11.04.23	changed argument tobs_key to pkey, generalized treatment
 *		of prime key to different types
 *	11.06.06	added keywords for output segment names to override the
 *		default values
 *	11.06.17	allow for output series without all three segments
 *	11.06.21	added setting of T_FIRST, T_LAST keys
 *	11.06.22	added averaging of CRLN_OBS
 *  v 1.0 frozen 11.11.14
 *	12.05.24	fixed averaging of CROTA* for cases where values may
 *		wrap; SAT_ROT values no longer averaged by default; provide
 *		option for overriding roll angle keyword; added setkey, setval
 *		options
 *  v 1.1 frozen 12.08.03
 *	12.09.12	added support for acceptance and/or rejection of
 *		records with certain values of CalVer64 (or other equivalent
 *		key)
 *	12.10.04	added recording of calling params info to comment or
 *		history keyword
 *	12.10.16	added writing of verbose output to log segment if
 *		present, and recording of input recnum's; added logging of
 *		CalVer64 values used; added support for specification of
 *		tmid in date_time format
 *  v 1.2 frozen 13.04.08
 *	14.03.01	added T_OBS to default list of keys to be averaged;
 *		added verbose warnings of missing output keys for propagation
 *		or averaging; removed automatic averaging of prime key
 *  v 1.3 frozen 14.07.07
 *	14.10.09	removed functions drms_appendstr_tokey() and
 *		append_args_tokey() (now in keystuff)
 *  v 1.4 frozen 15.01.13
 *	15.02.24	added initializations of several variables; made
 *		calculation of power values conditional on presence of the
 *		appropriate segment in the output record
 *  v 1.5 frozen 15.04.08
 */
