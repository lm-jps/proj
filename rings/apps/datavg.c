/*
 *  datavg.c						~rick/src/drms
 *
 *  Construct mean and variance matrices of selected data segments from an
 *    input records set (same format)
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
 *      tmid	str	unspecified	midpoint of averaging interval (CR:CL)
 *	length	float	unspecified	length of averaging interval (in deg og
 *			Carrington rotation)
 *	reject	string	Not Specified	If specified, name of a file with a
 *				list of input images to be rejected regardless
 *				of quality mask values
 *	qmask	int	0x80000000	Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	copy	str	+		list of keys to be propagated as-is
 *	average	str	+		list of keys to be averaged
 *	tobs_key str	T_OBS		Key name of time type keyword describing
 *			observation time (midpoint) of each input image
 *	qual_key str	Quality		Key name of uint type keyword describing
 *			data quality
 *
 *  Flags
 *	-n	no ingestion of output; diagnostic only
 *	-o	correct individual images for orbital velocity
 *	-v	run verbose
 *
 *  Status:
 *    Works with limited functionality; not widely tested
 *
 *  Bugs:
 *    Propagated keys are copied from the first input record, without checking
 *	for uniqueness
 *    If the input dataset is specified in the in string, the keywords for
 *	CarRot and CMLon are set to garbage, and there is no checking that
 *	the length is correct
 *    The MidTime value is geoecentric for midpoint Carrington longitude
 *    The correction for observer velocity only takes account of the radial
 *	component, and is only appropriate for Doppler data
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
char *version_id = "0.7";

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"in", "", "input data series or dataset"}, 
  {ARG_STRING,	"out", "", "output data series"}, 
  {ARG_STRING,	"tmid", "Not Specified",
      "midpoint of averaging interval (in Carrington longitude)"}, 
  {ARG_FLOAT,	"length", "Not Specified",
      "length of averaging interval (in degrees of Carrington rotation)"}, 
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_STRING,  "copy",  "+", "comma separated list of keys to propagate"},
  {ARG_STRING,  "average",  "+", "comma separated list of keys to average"},
  {ARG_STRING,	"tobs_key", "T_OBS", "keyname for image observation time"}, 
  {ARG_STRING,	"qual_key", "Quality", "keyname for 32-bit image quality field"}, 
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
char *average[] = {"SAT_ROT", "OBS_VR", "OBS_VW", "OBS_VN"};
			/*  the following belong in external utility files  */
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
    DRMS_Array_t *powr, int ntot) {
  double *vavg = (double *)mean->data;
  double *vvar = (double *)powr->data;
  int *vval = (int *)vcts->data;
  double mnvmin, mnvmax, pwrmin, pwrmax, norm, norm2, scale, sig, var;
  double sumv1, sump1, sum2, sum3, sum4;
  int sumv0, sump0;
  int valmin = ntot, valmax = 0;
  int n;
  int kstat = 0;

  mnvmin = pwrmin = HUGE;
  mnvmax = pwrmax = -HUGE;
  sumv1 = sump1 = 0.0;
  sumv0 = sump0 = 0;
/*
  mnvavg = pwravg = mnvvar = pwrvar = mnv3 = pwr3 = mnv4 = pwr4 = 0.0;
*/
  for (n = 0; n < ntot; n++) {
    if (vval[n] < valmin) valmin = vval[n];
    if (vval[n] > valmax) valmax = vval[n];
    if (vavg[n] < mnvmin) mnvmin = vavg[n];
    if (vavg[n] > mnvmax) mnvmax = vavg[n];
    if (vvar[n] < pwrmin) pwrmin = vvar[n];
    if (vvar[n] > pwrmax) pwrmax = vvar[n];
    if (isfinite (vavg[n])) {
      sumv0++;
      sumv1 += vavg[n];
    }
    if (isfinite (vvar[n])) {
      sump0++;
      sump1 += vvar[n];
    }
  }
  
  kstat += check_and_set_key_int (rec, "CountMIN", valmin);
  kstat += check_and_set_key_int (rec, "CountMAX", valmax);
  kstat += check_and_set_key_double (rec, "MeanMIN", mnvmin);
  kstat += check_and_set_key_double (rec, "MeanMAX", mnvmax);
  kstat += check_and_set_key_double (rec, "PowrMIN", pwrmin);
  kstat += check_and_set_key_double (rec, "PowrMAX", pwrmax);

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

  return kstat;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids = NULL;
  DRMS_Record_t *irec, *orec = NULL;
  DRMS_Segment_t *iseg, *vseg = NULL, *mseg = NULL, *pseg = NULL;
  DRMS_Array_t *data_array, *vcts, *mean, *powr;
  DRMS_Keyword_t *keywd;
  TIME tmid, tobs_rec, tobs;
  double *v, *vavg, *vvar;
  double *crpix, *crval, *cdelt, *crota, *avgval;
  double *crpixv, *crvalv, *cdeltv, *crotav, *avgvalv;
  double crpix_rec, crval_rec, cdelt_rec, crota_rec, avg_rec;
  double tobsv, vobs;
/*
  double crpix1, crpix2, cdelt1, cdelt2, crota2;
  double crpix1_rec, crpix2_rec, cdelt1_rec, cdelt2_rec, crota2_rec;
*/
  float clstrt, clmid, clstop;
  long long ntot, img_size;
  unsigned int quality;
  int *inaxis, *vval, *reject_list;
  int rec, recct, segct, segnum, maxct, imgct;
  int crstrt, crmid, crstop;
  int checkseg;
  int kstat, status;
  int n, naxis, wcsaxes;
  int propct, meanct, add_defaults;
  char **copykeylist, **meankeylist;
  char *source, *keystr;
  char recset_query[DRMS_MAXQUERYLEN], keyname[DRMS_MAXKEYNAMELEN];
  char module_ident[64], tbuf[64];

  double fp_nan = 0.0 / 0.0;
  int badqual = 0, blacklist = 0, rejects = 0;
  int propkeyct = sizeof (propagate) / sizeof (char *);
  int meankeyct = sizeof (average) / sizeof (char *);

  char *inset = strdup (params_get_str (params, "in"));
  char *out_series = strdup (params_get_str (params, "out"));
  char *tmid_str = strdup (params_get_str (params, "tmid"));
  float intrvl = params_get_float (params, "length");
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  char *rejectfile = strdup (params_get_str (params, "reject"));
  char *propagate_req = strdup (params_get_str (params, "copy"));
  char *average_req = strdup (params_get_str (params, "average"));
  char *tobs_key = strdup (params_get_str (params, "tobs_key"));
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  int no_save = params_isflagset (params, "n");
  int remove_obsvel = params_isflagset (params, "o");
  int verbose = params_isflagset (params, "v");
  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);
						   /*  create output record  */
  if (no_save) {
    orec = drms_create_record (drms_env, out_series, DRMS_TRANSIENT, &status);
    if (status) {
      fprintf (stderr, "Warning: drms_create_record returned %d for data series:\n",
	  status);
      fprintf (stderr, "  %s\nKeyword setting will not be checked", out_series);
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
    if (!no_save) return 1;
  }
  vseg = drms_segment_lookup (orec, "valid");
  mseg = drms_segment_lookup (orec, "mean");
  pseg = drms_segment_lookup (orec, "power");
  switch (vseg->info->type) {
    case (DRMS_TYPE_CHAR):
      maxct = 127;
      break;
    case (DRMS_TYPE_SHORT):
      maxct = 32767;
      break;
    default:
      maxct = 2147483648;
  }
  maxct += vseg->bzero;
  maxct /= vseg->bscale;

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
    clmid = fp_nan;
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
fprintf (stderr, "specification of time target by other than CR:CL not supported\n");
return 0;
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
      for (n = 1; n < propkeyct; n++)
	copykeylist[n + propct - 1] = propagate[n];
      propct = newct;
    }
  }
  if (verbose) {
    printf ("propagating values for %d key(s):\n", propct);
    for (n = 0; n < propct; n++) printf ("  %s\n", copykeylist[n]);
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
    printf ("averaging values for %d key(s):\n", meanct);
    for (n = 0; n < meanct; n++) printf ("  %s\n", meankeylist[n]);
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
  vvar = (double *)calloc (ntot, sizeof (double));
  vcts = drms_array_create (DRMS_TYPE_INT, naxis, iseg->axis, (void *)vval,
      &status);
  mean = drms_array_create (DRMS_TYPE_DOUBLE, naxis, iseg->axis, (void *)vavg,
      &status);
  powr = drms_array_create (DRMS_TYPE_DOUBLE, naxis, iseg->axis, (void *)vvar,
      &status);
  mean->bscale = mseg->bscale;
  mean->bzero = mseg->bzero;
  powr->bscale = pseg->bscale;
  powr->bzero = pseg->bzero;
  img_size = ntot;

  if (recct > maxct) {
    double bscale = vseg->bscale;
    while (recct > maxct) {
      bscale *= 0.5;
      maxct *= 2;
    }
    vcts->bscale = bscale;
    fprintf (stderr,
        "Warning: more records in input data set than maximum valid count\n");
    fprintf (stderr,
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
		 /*  support special hack of reading of rejection list file  */
  rejects = 0;
  if (strcmp (rejectfile, "Not Specified")) {
    FILE *rejectfp = fopen (rejectfile, "r");
    if (rejectfp) rejects = read_reject_list (rejectfp, &reject_list);
    else fprintf (stderr,
	"Warning: could not open rejection list %s; ignored\n", rejectfile);
  }
					 /*  initialize keys to be averaged  */
  tobs = tobsv = 0;
  imgct = 0;
					     /*  loop through input records  */
  for (rec = 0; rec < recct; rec++) {
    irec = ids->records[rec];
    if (verbose) {
      printf ("processing record ");
      tobs_rec = drms_getkey_time (irec, tobs_key, &status);
      if (status) printf (" %d (unknown time) : ", rec);
      else {
	sprint_time (tbuf, tobs_rec, "TAI", 0);
	printf ("%s (%d): ", tbuf, rec);
      }
    }
			 /*  check for record quality, reject as applicable  */
    quality = drms_getkey_int (irec, qual_key, &status);
    if ((quality & qmask) && !status) {
      badqual++;
      if (verbose) printf ("skipped (quality)\n");
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
        continue;
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
        continue;
      }
    }
    						  /*  read input data image  */
    data_array = drms_segment_read (iseg, DRMS_TYPE_DOUBLE, &status);
    if (status) {
      if (data_array) drms_free_array (data_array);
      if (verbose) printf ("skipped (segment read error)\n");
      else {
	fprintf (stderr, "Error reading file for record %d: ", rec);
	tobs_rec = drms_getkey_time (irec, tobs_key, &status);
	if (status) fprintf (stderr, "unknown time\n");
	else {
	  sprint_time (tbuf, tobs_rec, "TAI", 0);
	  fprintf (stderr, "%s\n", tbuf);
	}
      }
      continue;
    }

    tobs_rec = drms_getkey_time (irec, tobs_key, &status);
    if (status || time_is_invalid (tobs_rec)) {
      if (verbose) printf ("skipped (time invalid)\n");
      else fprintf (stderr, "error reading %s from record #%d\n", tobs_key, rec);
      if (data_array) drms_free_array (data_array);
      continue;
    }
    tobs += tobs_rec;
    tobsv += tobs_rec * tobs_rec;
    imgct++;

    v = (double *)data_array->data;
    vval = (int *)vcts->data;
    vavg = (double *)mean->data;
    vvar = (double *)powr->data;
    vobs = (remove_obsvel) ? drms_getkey_double (irec, "OBS_VR", &status) : 0.0;
						    /*  process valid data  */
    n = ntot;
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
	crota[n] += crota_rec;
	crotav[n] += crota_rec * crota_rec;
      }
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
  }
  if (checkseg) free (inaxis);
  vval = (int *)vcts->data;
  vavg = (double *)mean->data;
  vvar = (double *)powr->data;
  for (n = 0; n < ntot; n++) {
    if (vval[n]) {
      vavg[n] /= vval[n];
      vvar[n] /= vval[n];
      vvar[n] -= vavg[n] * vavg[n];
    } else vavg[n] = vvar[n] = fp_nan;
  }
  if (drms_segment_write (vseg, vcts, 0))
    fprintf (stderr, "Warning: unable to write to count segment\n");
  if (drms_segment_write (mseg, mean, 0))
    fprintf (stderr, "Warning: unable to write to mean segment\n");
  if (drms_segment_write (pseg, powr, 0))
    fprintf (stderr, "Warning: unable to write to variance segment\n");
						     /*  set remaining keys  */
  kstat = 0;
  kstat += check_and_set_key_int  (orec, "DataRecs", imgct);
  kstat += check_and_set_key_int  (orec, "MissRecs", recct - imgct);
							/*  statistics keys  */
  kstat += set_stats_keys (orec, vcts, mean, powr, ntot);
						      /*  averaged WCS keys  */
  if (imgct) {
    tobs /= imgct;
    kstat += check_and_set_key_time  (orec, tobs_key, tobs);
    tobsv /= imgct;
    tobsv -= tobs * tobs;
    sprintf (keyname, "D_%s", tobs_key);
    kstat += check_and_set_key_time  (orec, keyname, sqrt (tobsv));
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
      sprintf (keyname, "CROTA%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, crota[n]);
      crotav[n] /= imgct;
      crotav[n] -= crota[n] * crota[n];
      sprintf (keyname, "D_CROTA%d", n + 1);
      kstat += check_and_set_key_double (orec, keyname, sqrt (crotav[n]));
    }
						    /*  other averaged keys  */
    for (n = 0; n < meanct; n++) {
      avgval[n] /= imgct;
      kstat += check_and_set_key_double  (orec,  meankeylist[n], avgval[n]);
      avgvalv[n] /= imgct;
      avgvalv[n] -= avgval[n] * avgval[n];
      sprintf (keyname, "D_%s", meankeylist[n]);
      kstat += check_and_set_key_double (orec, keyname, sqrt (avgvalv[n]));
    }
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
    printf ("written\n");
  }
  drms_close_record (orec, dispose);

  if (verbose) {
    if (badqual)
      printf ("    %d input records rejected for quality matching %08x\n",
	  badqual, qmask);
    if (blacklist)
      printf ("    %d input records rejected from rejection list\n", blacklist);
  }
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
 *
 */
