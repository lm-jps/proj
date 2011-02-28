/*
 *  rdcover.c						~rick/hmi/rings/src
 *
 *  Responsible:  Rick Bogart				RBogart@spd.aas.org
 *
 *  Report effective coverage for selected mtrack arguments
 *
 *  (Main module body begins around line 560)
 *
 *  Parameters: (type   default         description)
 *	ds	DataSet none            Target dataset
 *				A series of records such as would be used
 *		as input to mtrack. It may be specified as either a data set
 *		query or just a data series name, in which case additional
 *		selection parameters are required.
 *		The dataset is assumed to include in its ancillary data
 *		the observation time, a quality parameter, and MISSVALS
 *	reject	string	-		If specified, name of a file with a
 *				list of input images to be rejected regardless
 *				of quality mask values
 *	qmask	int	0x80000000	Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	max_miss int	0		Tolerance threshold for number of blank
 *				values in image (assumed to exclude crop)
 *	tmid	string	-		Midpoint of target interval;
 *				ignored if input specified as data set; defined
 *				as string because it can be specified in either
 *				regular time format or as CR:CL
 *	length	int	-		Length of interval, in units
 *				of data cadence; ignored if input specified as
 *				data set
 *	tstart	time	-		Start of target tracking interval;
 *				ignored if input specified as data set, or if
 *				tmid and length are specified
 *	tstop	time	-		End of target tracking interval;
 *				ignored if input specified as data set, or if
 *				tmid and length are specified
 *	trec_key string	T_REC		Keyname of time type prime keyword for
 *				input data series
 *	tobs_key string	T_OBS		Keyname of time type keyword describing
 *				observation time (midpoint) of each input image
 *	tstp_key string	Cadence		Keyname of float type keyword describing
 *				observing cadence of input data
 *	qual_key string	Quality		Keyname of int type keyword describing
 *				record quality as bitmask
 *
 *  Flags
 *	-v	run verbose
 *
 *  Notes:
 *    This module is just the front-end processing for mtrack; it should
 *	perhaps just be a flag to that module
 *
 *  Bugs:
 *    If tstep is specified and does not match the data series cadence, the
 *	reported value of coverage is erroneous
 *    Checks for validity of tstart and tstop parameters are unreliable, due
 *	to bugs in DRMS treatment of scans of invalid strings for times
 *	(ticket #177)
 *    Will not accept an input dataset specification of @*
 *    There is no verification that numerous essential key data are actually
 *	in the input dataset, in particular trec_key, tobs_key, and
 *	crot_key
 *    When a target time and length are given, the actual number of records
 *	may differ by one from the expected number; this only occurs if
 *	the target time differs by a small but non-zero amount (< 0.1 sec)
 *	from either the data record time or the midway point between data
 *	record times. This happens for example, if the length is odd and
 *	the target time differs from an actual observation time by more than
 *	0.1 usec but less than 0.05 sec, with a 1-minute cadence
 *    When the target cadence is much larger than the input cadence, many
 *	input records are read needlessly; likewise, if there are multiple
 *	records for the same reference time
 *
 *  Revision history is at end of file
 */

#include <jsoc_main.h>

						      /*  module identifier  */
char *module_name = "rdcoverage";
char *module_desc = "report input data coverage for tracking options";
char *version_id = "0.8";

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds", "", "input data series or dataset"}, 
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_INT,	"max_miss", "0", "missing values threshold for image rejection"},
  {ARG_STRING,	"tmid", "Not Specified", "midpoint of tracking interval"}, 
  {ARG_INT,	"length", "0",
      "target length of tracking interval [input cadence]"}, 
				      /*  necessitated by bug (ticket #177)  */
  {ARG_STRING,	"tstart", "Not Specified", "start of coverage interval"}, 
  {ARG_STRING,	"tstop", "Not Specified", "end of coverage interval"}, 
  {ARG_FLOAT,	"tstep", "Not specified", "temporal cadence for output"},
  {ARG_STRING,	"trec_key", "T_REC", "keyname of (slotted) prime key"}, 
  {ARG_STRING,	"tobs_key", "T_OBS", "keyname for image observation time"}, 
  {ARG_STRING,	"tstp_key", "CADENCE",  "keyname for image observation time"}, 
  {ARG_STRING,	"qual_key", "Quality",  "keyname for 32-bit image quality field"}, 
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"crot_key", "CAR_ROT", "keyname for image Carrington rotation"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {}
};

#include "keystuff.c"
#include "selstuff.c"
#include "soho_ephem.c"

/*
 *  Functions to support reading from a specified rejection list of the
 *    appropriate format
 */

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

		 /*  global declaration of missing to be initialized as NaN  */
float missing_val;

static int verify_keys (DRMS_Record_t *rec, const char *clon,
    const char *clat, double *keyscale) {
  DRMS_Keyword_t *keywd;
  
  keywd = drms_keyword_lookup (rec, clon, 1);
  if (!keywd) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for Carrington longitude of observer not found\n",
	clon);
    fprintf (stderr, "       Must supply an appropriate value for clon_key\n");
    fprintf (stderr, "       (Carrington longitude of disc center in deg)\n");
    return -1;
  }
  if (keywd->info->type == DRMS_TYPE_TIME ||
      keywd->info->type == DRMS_TYPE_STRING ||
      keywd->info->type == DRMS_TYPE_RAW) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for observer Carrington longitude is of wrong type\n",
	clon);
    fprintf (stderr, "       Must supply an appropriate value for clon_key\n");
    fprintf (stderr, "       (Carrington longitude of disc centre in deg)\n");
    return -1;
  }
  if (strncasecmp (keywd->info->unit, "deg", 3)) {
    fprintf (stderr,
	"Warning: Keyword \"%s\" for observer Carrington longitude has unit of \"%s\"\n",
	clon, keywd->info->unit);
    fprintf (stderr, "         ignored, \"deg\" assumed\n");
  }

  keywd = drms_keyword_lookup (rec, clat, 1);
  if (!keywd) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for observer heliographic latitude not found\n",
	clat);
    fprintf (stderr, "       Must supply an appropriate value for clat_key\n");
    fprintf (stderr, "       (heliographic latitude of disc centre in deg)\n");
    return -1;
  }
  if (keywd->info->type == DRMS_TYPE_TIME ||
      keywd->info->type == DRMS_TYPE_STRING ||
      keywd->info->type == DRMS_TYPE_RAW) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for observer heliographic latitude is of wrong type\n",
	clat);
    fprintf (stderr, "       Must supply an appropriate value for clat_key\n");
    fprintf (stderr, "       (heliographic latitude of disc centre in deg)\n");
    return -1;
  }
  if (strncasecmp (keywd->info->unit, "deg", 3)) {
    fprintf (stderr,
	"Warning: Keyword \"%s\" for observer heliographic latitude has unit of \"%s\"\n",
	clat, keywd->info->unit);
    fprintf (stderr, "         ignored, \"deg\" assumed\n");
  }

  return 0;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ds = NULL;
  DRMS_Record_t *irec;
  DRMS_Keyword_t *keywd;
  TIME trec, tobs, tmid, tbase, tfirst, tlast, ttrgt, tstrt, tstop;
  double carr_lon, cm_lon_start, cm_lon_stop, lon_span;
  double cm_lon_first, cm_lon_last;
  double cadence, data_cadence, coverage, t_eps, phase;
  double lat, lon, img_lat, img_lon;
  double tmid_cl, lon_cm;
  double img_xc, img_yc, img_xscl, img_yscl, img_radius, img_pa;
  double kscale;

  long long nn;
  unsigned int quality;
  int *reject_list;
  int recct, rgn, rgnct, segct, valid;
  int tmid_cr, cr_start, cr_stop;
  int cr_first, cr_last;
  int col, row, pixct, i, found, n, nr, or;
  int blankvals, no_merid_v, rejects, status;
  int need_cadence, need_ephem;
  int badpkey, badqual, badfill, badtime, blacklist;
  char rec_query[256];
  char module_ident[64], key[64], tbuf[64], ptbuf[64], ctime_str[16];
  char keyvalstr[DRMS_DEFVAL_MAXLEN];
  enum platloc {LOC_UNKNOWN, LOC_MWO, LOC_GONG_MR, LOC_GONG_LE, LOC_GONG_UD,
      LOC_GONG_TD, LOC_GONG_CT, LOC_GONG_TC, LOC_GONG_BB, LOC_GONG_ML, LOC_SOHO,
      LOC_SDO} platform = LOC_UNKNOWN;

  int need_ephem_from_time = 0;
  int need_crcl = 1;
  int check_platform = 0;
  int extrapolate = 1;
  int found_first = 0, found_last = 0;
						 /*  process command params  */
  char *inset = strdup (params_get_str (params, "ds"));
  char *rejectfile = strdup (params_get_str (params, "reject"));
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  int max_miss = params_get_int (params, "max_miss");
  char *tmid_str = strdup (params_get_str (params, "tmid"));
  char *tstrt_str = strdup (params_get_str (params, "tstart"));
  char *tstop_str = strdup (params_get_str (params, "tstop"));
  int length = params_get_int (params, "length");
  double tstep = params_get_double (params, "tstep");

  char *trec_key = strdup (params_get_str (params, "trec_key"));
  char *tobs_key = strdup (params_get_str (params, "tobs_key"));
  char *tstp_key = strdup (params_get_str (params, "tstp_key"));
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *crot_key = strdup (params_get_str (params, "crot_key"));

  int verbose = params_isflagset (params, "v");

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s: JSOC version %s\n", module_ident, jsoc_version);
	  /*  get lists of latitudes and longitudes defining region centers  */
  need_cadence = isnan (tstep);

  if (key_params_from_dspec (inset)) {
				   /*  input specified as specific data set  */
    if (!(ds = drms_open_records (drms_env, inset, &status))) {
      fprintf (stderr, "Error: (%s) unable to open input data set %s\n",
        module_ident, inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    if ((recct = ds->n) < 2) {
      printf (" <2 records in selected input set\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 0;
    }
    irec = ds->records[0];
    status = verify_keys (irec, clon_key, clat_key, &kscale);
    if (status) {
      fprintf (stderr, "Error: unable to verify keys\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    if (need_cadence) {
			 /*  output cadence not specified, use data cadence  */
      if ((keywd = drms_keyword_lookup (irec, tstp_key, 1))) {
					     /*  cadence should be constant  */
	if (keywd->info->recscope != 1) {
	  fprintf (stderr, "Warning: cadence is variable in input series %s\n",
	      inset);
	  fprintf (stderr, "         output cadence will be %f\n",
	       drms_getkey_double (irec, tstp_key, &status));
	}
      } else {
			 /*  could infer from slotting info as well, but...  */
	fprintf (stderr,
          "Error: data cadence keyword %s not in input series %s\n",
	  tstp_key, inset);
	fprintf (stderr, "       Specify desired output cadence as tstep\n");
	drms_close_records (ds, DRMS_FREE_RECORD);
	return 1;
      }
      data_cadence = drms_getkey_double (irec, tstp_key, &status);
      tstep = data_cadence;
    }
    tstrt = drms_getkey_time (irec, trec_key, &status);
    tstop = drms_getkey_time (ds->records[recct - 1], trec_key, &status);
    length = (tstop - tstrt + 1.01 * tstep) / tstep;
    tmid = 0.5 * (tstrt + tstop);
    segct = drms_record_numsegments (irec);
  } else {
				/*  only the input data series is named,
				   get record specifications from arguments  */
		   /*  get required series info from first record in series  */
					       /*  platform, cadence, phase  */
    snprintf (rec_query, 256, "%s[#^]", inset);
    if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    irec = ds->records[0];
    status = verify_keys (irec, clon_key, clat_key, &kscale);
    if (status) {
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    if ((keywd = drms_keyword_lookup (irec, tstp_key, 1))) {
					     /*  cadence should be constant  */
      if (keywd->info->recscope != 1) {
	fprintf (stderr, "Warning: cadence is variable in input series %s\n",
	    inset);
	if (need_cadence)
	  fprintf (stderr, "         output cadence will be %f\n",
	      drms_getkey_double (irec, tstp_key, &status));
      }
    } else {
			 /*  could infer from slotting info as well, but...  */
      fprintf (stderr,
          "Error: data cadence keyword %s not in input series %s\n",
	  tstp_key, inset);
      fprintf (stderr, "       Specify desired output cadence as tstep\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    data_cadence = drms_getkey_double (irec, tstp_key, &status);
    t_eps = 0.5 * data_cadence;
    if (need_cadence)
			 /*  output cadence not specified, use data cadence  */
      tstep = data_cadence;
    if ((keywd = drms_keyword_lookup (irec, "TELESCOP", 1))) {
						     /*  should be constant  */
      if (keywd->info->recscope != 1)
	fprintf (stderr, "Warning: TELESCOP is variable in input series %s\n",
	    inset);
      if (!strcmp (drms_getkey_string (irec, "TELESCOP", &status), "SDO/HMI"))
	platform = LOC_SDO;
      else if (!strcmp (drms_getkey_string (irec, "TELESCOP", &status), "SOHO"))
	platform = LOC_SOHO;
      else if (!strcmp (drms_getkey_string (irec, "TELESCOP", &status), "NSO-GONG")) {
	if ((keywd = drms_keyword_lookup (irec, "SITE", 1))) {
	  if (!strcmp (drms_getkey_string (irec, "SITE", &status), "MR"))
	    platform = LOC_GONG_MR;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "LE"))
	    platform = LOC_GONG_LE;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "UD"))
	    platform = LOC_GONG_UD;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "TD"))
	    platform = LOC_GONG_TD;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "CT"))
	    platform = LOC_GONG_CT;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "BB"))
	    platform = LOC_GONG_BB;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "ML"))
	    platform = LOC_GONG_ML;
	  else {
	    platform = LOC_GONG_MR;
	  }
	} else {
          fprintf (stderr, "Warning: unspecified GONG site: MR assumed\n");
	  platform = LOC_GONG_MR;
	}
      }
    }
    if (platform == LOC_UNKNOWN) {
      fprintf (stderr, "Warning: observing location unknown, assumed geocenter\n");
    }
    segct = drms_record_numsegments (irec);

    if (strcmp (tmid_str, "Not Specified")) {
/*  determine start and stop times from length (in units of tstep) and midtime
				   (which can be CR:CL as well as date_time) */
       if (sscanf (tmid_str, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
			  /*  tmid specified as CR:CL : need ephemeris info  */
	need_crcl = 0;
	if (platform == LOC_SDO || platform == LOC_SOHO ||
	    (platform >= LOC_GONG_MR && platform <= LOC_GONG_ML) ||
	    platform == LOC_UNKNOWN) {
	  if (platform == LOC_SOHO)
	    tmid = SOHO_meridian_crossing (tmid_cl, tmid_cr);
	  else
	    tmid = earth_meridian_crossing (tmid_cl, tmid_cr);
	  sprint_time (ptbuf, tmid, "", 0);
			  /*  adjust to phase of input, within data cadence  */
	  tbase = drms_getkey_time (irec, trec_key, &status);
	  phase = fmod ((tmid - tbase), data_cadence);
	  tmid -= phase;
	  if (phase > 0.5 * data_cadence) tmid += data_cadence;
	  if (verbose) {
	    sprint_time (tbuf, tmid, "", 0);
	    printf ("Target time %d:%05.1f = %s adjusted to %s\n",
	        tmid_cr, tmid_cl, ptbuf, tbuf);
	  }
	} else {
	  fprintf (stderr,
	      "Time specification in CR:CL for observing platform not supported\n");
	  drms_close_records (ds, DRMS_FREE_RECORD);
	  return 1;
	}
      } else {
			      /*  tmid specified as normal date-time string  */
        tmid = sscan_time (tmid_str);
      }
      tbase = drms_getkey_time (irec, trec_key, &status);
      phase = fmod ((tmid - tbase), data_cadence);
      tstrt = tmid - 0.5 * length * tstep;
      tstop = tstrt + (length - 1) * tstep;
      if (tstop <= tstrt) {
	fprintf (stderr,
	    "Error: requested end time before or at start time for ");
	fprintf (stderr, "length= %d\n", length);
	drms_close_records (ds, DRMS_FREE_RECORD);
	return 1;
      }
			  /*  adjust stop time to reflect sampling symmetry  */
      if ((fabs (phase) < 0.001 * t_eps) && length % 2)
	tstop += tstep;
      if ((fabs (phase - t_eps) < 0.001 * t_eps) && (length % 2 == 0))
	tstop += tstep;
    } else {
	       /*  tstart and tstop specified, determine midtime and length  */
      if (sscanf (tstrt_str, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
	if (platform == LOC_SOHO)
	  tstrt = SOHO_meridian_crossing (tmid_cl, tmid_cr);
	else
	  tstrt = earth_meridian_crossing (tmid_cl, tmid_cr);
      } else tstrt = sscan_time (tstrt_str);
      if (sscanf (tstop_str, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
	if (platform == LOC_SOHO)
	  tstop = SOHO_meridian_crossing (tmid_cl, tmid_cr);
	else
	  tstop = earth_meridian_crossing (tmid_cl, tmid_cr);
      } else tstop = sscan_time (tstop_str);
      tmid = 0.5 * (tstrt + tstop);
      length = (tstop - tstrt + 1.01 * tstep) / tstep;
    }
    drms_close_records (ds, DRMS_FREE_RECORD);

    snprintf (rec_query, 256, "%s[?%s > %23.16e and %s < %23.16e?]", inset,
	trec_key, tstrt - t_eps, trec_key, tstop + t_eps);
    if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
      fprintf (stderr, "       status = %d\n", status);
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    if ((recct = ds->n) < 2) {
      printf ("<2 records in selected input set\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 0;
    }
/*
    if (strcmp (tmid_str, "Not Specified")) {
      if (!(ds = select_dataset_from_time_interval (inset, tmid_str, length * tstep))){
	fprintf (stderr, "Error, unable to open dataset from series %s\n",
	  inset);
	return 1;
      }
      if ((recct = ds->n) < 2) {
	printf ("<2 records in selected input set\n");
	drms_close_records (ds, DRMS_FREE_RECORD);
	return 0;
      }
    } else {
      if (!strcmp (tstrt_str, "Not Specified") ||
	  !strcmp (tstop_str, "Not Specified")) {
	fprintf (stderr,
	    "Error: either a specific data record set must be selected as input\n");
	fprintf (stderr, "       or (tmid and length) or (tstart and tstop) must be\n");
	fprintf (stderr, "       specified\n");
	return 1;
      }
      if (!(ds = select_dataset_from_time_range (inset, tstrt_str, tstop_str))) {
	fprintf (stderr, "Error, unable to open dataset from series %s\n",
	    inset);
	return 1;
      }
      if ((recct = ds->n) < 2) {
	printf ("<2 records in selected input set\n");
	drms_close_records (ds, DRMS_FREE_RECORD);
	return 0;
      }
    }
*/
  }
			    /*  end determination of record set from params  */
  if (verbose) {
    sprint_time (ptbuf, tstrt, "", -1);
    sprint_time (tbuf, tstop, "", -1);
    printf ("checking data from %s - %s at cadence of %.1f s\n", ptbuf, tbuf,
	tstep);
    printf ("         %d data records\n", recct);
  }

  cadence = data_cadence;
  tobs = drms_getkey_time (ds->records[0], tobs_key, &status);
  if (fabs (tobs - tstrt) < cadence) {
    tfirst = tobs;
    irec = ds->records[0];
    cm_lon_first = cm_lon_start = drms_getkey_double (irec, clon_key, &status);
    cr_first = cr_start = drms_getkey_int (irec, crot_key, &status);
    found_first = 1;
  } else need_ephem_from_time = 1;

  tobs = drms_getkey_time (ds->records[recct - 1], tobs_key, &status);
  if (fabs (tobs - tstop) < cadence) {
    tlast = tobs;
    irec = ds->records[recct - 1];
    cm_lon_last = cm_lon_stop = drms_getkey_double (irec, clon_key, &status);
    cr_last = cr_stop = drms_getkey_int (irec, crot_key, &status);
    found_last = 1;
  } else need_ephem_from_time = 1;

  found = 0;
  for (nr = 0; nr < recct; nr++) {
    tobs = drms_getkey_time (ds->records[nr], tobs_key, &status);
    if (fabs (tobs - tmid) < cadence) {
      irec = ds->records[nr];
      if (need_crcl) {
        tmid_cl = drms_getkey_double (irec, clon_key, &status);
        tmid_cr = drms_getkey_int (irec, crot_key, &status);
      }
      found = 1;
    }
    if (!found_first && !time_is_invalid (tobs)) {
      irec = ds->records[nr];
      cm_lon_first = drms_getkey_double (irec, clon_key, &status);
      cr_first = drms_getkey_int (irec, crot_key, &status);
      tfirst = tobs;
      found_first = 1;
    }
  }
  if (!found_last) {
    for (nr = recct - 1; nr; nr--) {
      if (!time_is_invalid (tobs)) {
	if (tobs < tmid) {
	  cm_lon_last = tmid_cl;
	  cr_last = tmid_cr;
	} else {
	  irec = ds->records[nr];
	  cm_lon_last = drms_getkey_double (irec, clon_key, &status);
	  cr_last = drms_getkey_int (irec, crot_key, &status);
	}
	tlast = tobs;
	found_last = 1;
	break;
      }
    }
  }
  if (!found) need_ephem_from_time = 1;

  if (need_ephem_from_time) {
    double rsun, vr, vn, vw;
    TIME table_mod_time;
    if (platform == LOC_SDO) {
      earth_ephemeris (tstrt, &rsun, &img_lat, &cm_lon_start, &vr, &vn, &vw);
      cr_start = carrington_rots (tstrt, 1);
      earth_ephemeris (tmid, &rsun, &img_lat, &tmid_cl, &vr, &vn, &vw);
      tmid_cr = carrington_rots (tmid, 1);
      earth_ephemeris (tstop, &rsun, &img_lat, &cm_lon_stop, &vr, &vn, &vw);
      cr_stop = carrington_rots (tstop, 1);
    } else if (platform == LOC_SOHO) {
      soho_ephemeris (tstrt, &rsun, &img_lat, &cm_lon_start, &vr, &vn, &vw,
	  &table_mod_time);
      cr_start = carrington_rots (tstrt, 1);
      soho_ephemeris (tmid, &rsun, &img_lat, &tmid_cl, &vr, &vn, &vw,
	  &table_mod_time);
      tmid_cr = carrington_rots (tmid, 1);
      soho_ephemeris (tstop, &rsun, &img_lat, &cm_lon_stop, &vr, &vn, &vw,
	  &table_mod_time);
      cr_stop = carrington_rots (tstop, 1);
    } else {
printf ("need ephem from time\n");
      if (!found_first || !found_last) {
        fprintf (stderr, "Error: Carrington ephemeris from time not supported\n");
        fprintf (stderr, "         and no valid times in data for estimation!\n");
	return 1;
      }
		/*  estimate midpoint ephemeris by linear interpolation
					 of closest observations to midtime  */
 /*  This code has not been well tested in crossover of Carrington rotation  */
    	  /*  assume at most one rotation between first and last estimators  */
      fprintf (stderr, "Warning: Carrington ephemeris from time not supported\n");
      if (!found) {
	tmid_cr = cr_first;
	if (cr_last != cr_first) cm_lon_last -= 360.0;
	tmid_cl = cm_lon_first +
            (cm_lon_last - cm_lon_first) * (tmid - tfirst) / (tlast - tfirst);
    	      /*  assume at most one rotation between first and last estimators  */
	if (tmid_cl < 0.0) {
	  tmid_cr++;
	  tmid_cl += 360.0;
	}
      }
      fprintf (stderr, "         estimating midpoint as %d:%08.4f\n",
	  tmid_cr, tmid_cl);
	/*  long extrapolations, and lazy correction for change of rotation
					number,but only needed for lon_span  */
      cr_start = cr_first;
      cm_lon_start = cm_lon_first +
          (cm_lon_last - cm_lon_first) * (tstrt - tfirst) / (tlast - tfirst);
      cr_stop = cr_last;
      cm_lon_stop = cm_lon_first +
          (cm_lon_last - cm_lon_first) * (tstop - tfirst) / (tlast - tfirst);
      if (cm_lon_stop > cm_lon_start) cr_stop++;
    }
  }
  lon_span = cm_lon_start - cm_lon_stop;
  while (cr_start < cr_stop) {
    cr_start++;
    lon_span += 360.0;
  }
		 /*  support special hack of reading of rejection list file  */
  rejects = 0;
  if (strcmp (rejectfile, "Not Specified")) {
    FILE *rejectfp = fopen (rejectfile, "r");
    if (rejectfp) rejects = read_reject_list (rejectfp, &reject_list);
    else fprintf (stderr,
	"Warning: could not open rejection list %s; ignored\n", rejectfile);
  }

  found = 0;
					     /*  loop through input records  */
  valid = 0;
  ttrgt = tstrt;
  or = 0;
  badpkey = badqual = badfill = badtime = blacklist = 0;
  for (nr = 0; nr < recct; nr++) {
    irec = ds->records[nr];
    quality = drms_getkey_int (irec, qual_key, &status);
    if ((quality & qmask) && !status) {
						    /*  partial image, skip  */
      badqual++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s = %08x matches %08x\n", keyvalstr, qual_key, quality,
	    qmask);
      }
      continue;
    }
    tobs = drms_getkey_time (irec, tobs_key, &status);
    if (time_is_invalid (tobs)) {
							  /*  no data, skip  */
      badpkey++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s invalid\n", keyvalstr, tobs_key);
      }
      continue;
    }
    blankvals = drms_getkey_int (irec, "MISSVALS", &status);
    if (blankvals > max_miss && !status) {
						    /*  partial image, skip  */
      badfill++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s = %d > %d\n", keyvalstr, "MISSVALS", blankvals,
	    max_miss);
      }
      continue;
    }
    if (tobs == tlast) {
					 /*  same time as last record, skip  */
      badtime++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s = ", keyvalstr, tobs_key);
	drms_keyword_snprintfval (drms_keyword_lookup (irec, tobs_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	sprint_time (tbuf, tlast, "", 0);
	printf ("%s = %s\n", keyvalstr, tbuf);
      }
      continue;
    }
			     /*  replace with call to solar_ephemeris_info?  */
    img_lon = drms_getkey_double (irec, clon_key, &status);
    img_lat = drms_getkey_double (irec, clat_key, &status);
			 /*  check for record quality, reject as applicable  */
    if (rejects) {
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
        continue;
      }
    }
		      /*  loop through output records, appending time slice  */
	    /*  extrapolate first image backward to start time if necessary  */
    if (extrapolate) {
      while (ttrgt < tobs) {

	if (verbose) printf ("step %d to be extrapolated from image %s\n", or, tbuf);

	or++;
	if (or >= length) {
	  fprintf (stderr, "Error: reached output length limit\n");
	  drms_close_records (ds, DRMS_FREE_RECORD);
	  return 1;
	}
	ttrgt += tstep;
      }
      extrapolate = 0;
    }
    if (ttrgt < tobs) {
	/*  linearly interpolate individual pixels between last valid map
								and current  */
      double f, g;
      int ct, ntot = rgnct * pixct;
      float *val = (float *)malloc (ntot * sizeof (float));
      char *skip = (char *)malloc (ntot * sizeof (char));
      
      while (ttrgt < tobs) {
        f = (ttrgt - tlast) / (tobs - tlast);
	g = 1.0 - f;

	if (verbose) {
	  sprint_time (ptbuf, tlast, "", 0);
	  sprint_time (tbuf, tobs, "", 0);
	  printf ("step %d to be interpolated from images %s and %s\n",
	      or, ptbuf, tbuf);
	}

	or++;
	if (or >= length) {
	  if (nr < recct - 1)
	    fprintf (stderr,
	      "Warning: reached output limit before last input record processed\n");
	  ttrgt = tstop;
	}
	ttrgt += tstep;
      }
      free (val);
      free (skip);
    }

    if (ttrgt == tobs) {
/*
      if (verbose) printf ("step %d mapped from image %s\n", or, tbuf);
*/
      or++;
      ttrgt += tstep;
    }
    tlast = tobs;
    strcpy (ptbuf, tbuf);
    valid++;
  }

  drms_close_records (ds, DRMS_FREE_RECORD);

  coverage = (double)valid / (double)length;
  printf ("%d of %d possible input records accepted\n", valid, recct);
  printf ("    effective coverage = %.3f\n", coverage);
  if (badpkey)
    printf ("    %d input records rejected for invalid values of %s\n",
	badpkey, tobs_key);
  if (badqual)
    printf ("    %d input records rejected for quality matching %08x\n",
	badqual, qmask);
  if (badfill)
    printf ("    %d input records rejected for missing values exceeding %d\n",
	badfill, max_miss);
  if (badtime)
    printf ("    %d input records rejected for duplicate values of %s\n",
	badtime, tobs_key);
  if (blacklist)
    printf ("    %d input records rejected from rejection list\n", blacklist);
  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  10.04.17 	created this file, based on mtrack; versions < 0.7 not saved
 *  v 0.7	added some verbose reporting, extracted some utility functions
 *  v 0.7 frozen 2010.08.19
 *  10.10.28	moved quality flag check ahead of T_OBS check
 *  v 0.8 frozen 2011.02.28
 *
 */
