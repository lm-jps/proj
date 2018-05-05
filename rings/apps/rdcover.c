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
 *	max_miss int	All		Tolerance threshold for number of blank
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
 *	tstep	float	-		Temporal cadence for (track) output
 *	trec_key string	T_REC		Keyname of time type prime keyword for
 *				input data series
 *	tobs_key string	T_OBS		Keyname of time type keyword describing
 *				observation time (midpoint) of each input image
 *	tstp_key string	Cadence		Keyname of float type keyword describing
 *				observing cadence of input data
 *	qual_key string	Quality		Keyname of int type keyword describing
 *				record quality as bitmask
 *	pa	float	-		Value of position angle required for
 *				acceptance
 *	dpa	float	-		Acceptance window width for position
 *				angle
 *	pa_key string	CROTA2		Keyname of float type keyword describing
 *				position angle of (solar) north on image
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
 *    There is no verification that some essential key data are actually
 *	in the input dataset, in particular tobs_key and crot_key
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
 *    The call to drms_record_getvector uses a double rather than equivalent
 *	time type due to a bug in the function scheduled to be fixed with the
 *	DRMS 8.11 release
 *
 *  Revision history is at end of file
 */

#include <jsoc_main.h>
						       /*  module identifier  */
char *module_name = "rdcoverage";
char *module_desc = "report input data coverage for tracking options";
char *version_id = "1.2";

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds", "", "input data series or dataset"}, 
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_INT,	"max_miss", "All",
      "missing values threshold for image rejection"},
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
  {ARG_STRING,	"qual_key", "Quality",
      "keyname for 32-bit image quality field"}, 
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"crot_key", "CAR_ROT", "keyname for image Carrington rotation"},
  {ARG_FLOAT,	"pa", "Not Specified",
      "centre of acceptable roll angles [deg]"}, 
  {ARG_FLOAT,	"dpa", "Not Specified",
      "maximum deviation of acceptable roll angles [deg]"}, 
  {ARG_STRING,	"pa_key", "CROTA2", "keyname for image position angle"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {}
};

#include "earth_ephem.c"
#include "soho_ephem.c"

static int key_params_from_dspec (const char *dspec) {
/*
 *  Establish whether target times are determined from dataset specifier
 *  assume that if a bracket is found in the dataset specifier it is a
 *  set of well-specified records containing a properly ordered input dataset;
 *  otherwise, a dataseries to be queried
 */
  int n, nt = strlen (dspec);

  for (n = 0; n < nt; n++) if (dspec[n] == '[') return 1;
  return 0;
}

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

static int verify_keys (DRMS_Record_t *rec, const char *trec, const char *clon,
    const char *clat, int *slotted) {
  DRMS_Keyword_t *keywd;

  keywd = drms_keyword_lookup (rec, trec, 1);
  if (!keywd) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for slotted time not found\n", trec);
    fprintf (stderr, "       Must supply an appropriate value for trec_key\n");
    return -1;
  }
				       /*  check whether trec key is slotted  */
  *slotted = drms_keyword_isslotted (keywd);
  if (!(*slotted))
    fprintf (stderr, "Warning: Keyword \"%s\" is not slotted\n", trec);

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
  DRMS_Array_t *tvec;
  DRMS_Keyword_t *keywd;
  TIME trec, tobs, tmid, tbase, tfirst, tlast, tstrt, tstop;
  double data_cadence, coverage, t_eps, phase;
  double lat, lon, img_lat, img_lon;
  double tmid_cl;
  double img_xc, img_yc, img_xscl, img_yscl, img_radius, img_pa;
  float pa_rec, dpa;

  long long nn;
  unsigned int quality;
  int *reject_list;
  int recct, rgn, rgnct, segct, valid;
  int tmid_cr;
  int col, row, pixct, i, found, n, nr;
  int blankvals, no_merid_v, rejects, status;
  int need_cadence;
  int badpkey, badqual, badfill, badtime, badpa, blacklist;
  int pkeyslotted;
  char rec_query[256];
  char module_ident[64], key[64], tbuf[64], ptbuf[64], ctime_str[16];
  char keyvalstr[DRMS_DEFVAL_MAXLEN];
  enum platloc {LOC_UNKNOWN, LOC_MWO, LOC_GONG_MR, LOC_GONG_LE, LOC_GONG_UD,
      LOC_GONG_TD, LOC_GONG_CT, LOC_GONG_TC, LOC_GONG_BB, LOC_GONG_ML, LOC_SOHO,
      LOC_SDO} platform = LOC_UNKNOWN;

  int check_platform = 0;
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
  float pa_nom = params_get_float (params, "pa");
  float dpa_max = params_get_float (params, "dpa");
  int check_pa = (isfinite (pa_nom) && dpa_max < 180.0);

  char *trec_key = strdup (params_get_str (params, "trec_key"));
  char *tobs_key = strdup (params_get_str (params, "tobs_key"));
  char *tstp_key = strdup (params_get_str (params, "tstp_key"));
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *crot_key = strdup (params_get_str (params, "crot_key"));
  char *pang_key = strdup (params_get_str (params, "pa_key"));

  int verbose = params_isflagset (params, "v");

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s: JSOC version %s\n", module_ident, jsoc_version);
	   /*  get lists of latitudes and longitudes defining region centers  */
  need_cadence = isnan (tstep);

  if (key_params_from_dspec (inset)) {
    TIME *tvals;
				    /*  input specified as specific data set  */
    if (!(ds = drms_open_recordset (drms_env, inset, &status))) {
      fprintf (stderr, "Error: (%s) unable to open input data set %s\n",
        module_ident, inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
    if (!irec) {
      printf (" No records in selected input set\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    status = verify_keys (irec, trec_key, clon_key, clat_key, &pkeyslotted);
    if (status) {
      fprintf (stderr, "Error: unable to verify keys\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    tvec = drms_record_getvector (drms_env, inset, trec_key, DRMS_TYPE_DOUBLE,
	0, &status);
    tvals = (TIME *)tvec->data;
    recct = tvec->axis[1];
    if (recct  < 2) {
      printf ("<2 records in selected input set\n");
      drms_free_array (tvec);
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 0;
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
    tstrt = tvals[0];
    tstop = tvals[recct - 1];
    length = (tstop - tstrt + 1.01 * tstep) / tstep;
    tmid = 0.5 * (tstrt + tstop);
    segct = drms_record_numsegments (irec);
  } else {
			/*  only the input data series is named:
				    get record specifications from arguments  */
		    /*  get required series info from series template record  */
						/*  platform, cadence, phase  */
    irec = drms_template_record (drms_env, inset, &status);
    if (status) {
      fprintf (stderr, "Error: unable to create template record in %s\n", inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    status = verify_keys (irec, trec_key, clon_key, clat_key, &pkeyslotted);
    if (status) {
      fprintf (stderr, "Error: unable to verify keys from template record\n");
      drms_free_record (irec);
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
      drms_free_record (irec);
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
      else if (!strcmp (drms_getkey_string (irec, "TELESCOP", &status),
          "NSO-GONG")) {
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
      fprintf (stderr,
          "Warning: observing location unknown, assumed geocenter\n");
    }
    segct = drms_record_numsegments (irec);

    if (strcmp (tmid_str, "Not Specified")) {
	/*  determine start and stop times from length (in units of tstep)
		       and midtime (which can be CR:CL as well as date_time)  */
       if (sscanf (tmid_str, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
			   /*  tmid specified as CR:CL : need ephemeris info  */
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
	  drms_free_record (irec);
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
	drms_free_record (irec);
	return 1;
      }
			  /*  adjust stop time to reflect sampling symmetry  */
      if ((fabs (phase) < 0.001 * t_eps) && length % 2)
	tstop += tstep;
      if ((fabs (phase - t_eps) < 0.001 * t_eps) && (length % 2 == 0))
	tstop += tstep;
    } else {
			   /*  tstart and tstop specified, determine length  */
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
      if (time_is_invalid (tstrt) || time_is_invalid (tstop)) {
	fprintf (stderr,
	    "Error: if data set not explicitly selected, either tmid and length\n");
	fprintf (stderr, "       or tstart and tstop must be specified\n");
	drms_free_record (irec);
	return 1;
      }
      length = (tstop - tstrt + 1.01 * tstep) / tstep;
    }
    drms_free_record (irec);
    if (pkeyslotted) {
      TIME pbase;
      double pstep;
      int tstrt_ind, tstop_ind;
      char pkeyindx[32], pkeyepoch[32], pkeystep[32];

      sprintf (pkeyindx, "%s_index", trec_key);
      sprintf (pkeyepoch, "%s_epoch", trec_key);
      sprintf (pkeystep, "%s_step", trec_key);
      pbase = drms_getkey_time (irec, pkeyepoch, &status);
      pstep = drms_getkey_double (irec, pkeystep, &status);
      tstrt_ind = (tstrt - pbase) / pstep;
      tstop_ind = (tstop + t_eps - pbase) / pstep;
      snprintf (rec_query, 256, "%s[?%s between %d and %d?]", inset,
	  pkeyindx, tstrt_ind, tstop_ind);
    } else {
      snprintf (rec_query, 256, "%s[?%s between %23.16e and %23.16e?]", inset,
	  trec_key, tstrt - t_eps, tstop + t_eps);
    }
    if (!(ds = drms_open_recordset (drms_env, rec_query, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
      fprintf (stderr, "       status = %d\n", status);
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
    if (!irec) {
      printf (" No records in selected input set\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
    tvec = drms_record_getvector (drms_env, rec_query, trec_key, DRMS_TYPE_DOUBLE,
	0, &status);
    recct = tvec->axis[1];
    if (recct  < 2) {
      printf ("<2 records in selected input set\n");
      drms_free_array (tvec);
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 0;
    }
  }
  drms_free_array (tvec);
			     /*  end determination of record set from params  */
  if (verbose) {
    sprint_time (ptbuf, tstrt, "TAI", 0);
    sprint_time (tbuf, tstop, "TAI", 0);
    printf ("checking data from %s - %s\n  at cadence of %.1f s\n", ptbuf, tbuf,
	tstep);
    printf ("         %d data records\n", recct);
  }

  /*  At this point, we have an input record set, tstrt, tstop, and tstep,
  				all that is needed for the checking loop  */

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
  badpkey = badqual = badfill = badtime = badpa = blacklist = 0;
  tlast = tstrt - 2 * tstep;
  for (nr = 0; nr < recct; nr++) {
    quality = drms_getkey_int (irec, qual_key, &status);
    if ((quality & qmask) && !status) {
					 /*  unacceptable quality bits, skip  */
      badqual++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s = %08x matches %08x\n", keyvalstr, qual_key, quality,
	    qmask);
      }
      irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
      continue;
    }
    tobs = drms_getkey_time (irec, tobs_key, &status);
    if (time_is_invalid (tobs)) {
						   /*  invalid obstime, skip  */
      badpkey++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s invalid\n", keyvalstr, tobs_key);
      }
      irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
      continue;
    } else {
      if (tobs == tlast) {
				       /*  same obstime as last record, skip  */
	badtime++;
	if (verbose) {
	  drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	      keyvalstr, DRMS_DEFVAL_MAXLEN);
	  printf ("%s: %s = ", keyvalstr, tobs_key);
	  drms_keyword_snprintfval (drms_keyword_lookup (irec, tobs_key, 1),
	      keyvalstr, DRMS_DEFVAL_MAXLEN);
	  printf ("%s duplicated\n", keyvalstr);
	}
	irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
	continue;
      }
      tlast = tobs;
    }
    blankvals = drms_getkey_int (irec, "MISSVALS", &status);
    if ((max_miss >= 0) && (blankvals > max_miss) && !status) {
					     /*  too many blank values, skip  */
      badfill++;
      if (verbose) {
	drms_keyword_snprintfval (drms_keyword_lookup (irec, trec_key, 1),
	    keyvalstr, DRMS_DEFVAL_MAXLEN);
	printf ("%s: %s = %d > %d\n", keyvalstr, "MISSVALS", blankvals,
	    max_miss);
      }
      irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
      continue;
    }
    if (check_pa) {
			      /*  check for non-nominal image rotation angle  */
      pa_rec = drms_getkey_float (irec, pang_key, &status);
      if (status) {
	fprintf (stderr, "Warning: \"%s\" keyword not found\n", pang_key);
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
	  if (verbose) {
	    drms_keyword_snprintfval (drms_keyword_lookup (irec, pang_key, 1),
		keyvalstr, DRMS_DEFVAL_MAXLEN);
	    printf ("%s: |%s - %.2f] = %.2f > %.2f\n", keyvalstr, pang_key,
		pa_nom, dpa, dpa_max);
	  }
	  irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
	  continue;
	}
      }
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
	irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
        continue;
      }
    }

    valid++;
    irec = drms_recordset_fetchnext (drms_env, ds, &status, NULL, NULL);
  }

  drms_close_records (ds, DRMS_FREE_RECORD);

  coverage = (double)valid / (double)length;
  printf ("%d of %d possible input records accepted\n", valid, recct);
  printf ("    effective coverage = %.3f\n", coverage);
  if (badqual)
    printf ("    %d input records rejected for quality matching %08x\n",
	badqual, qmask);
  if (badpkey)
    printf ("    %d input records rejected for invalid values of %s\n",
	badpkey, tobs_key);
  if (badtime)
    printf ("    %d input records rejected for duplicate values of %s\n",
	badtime, tobs_key);
  if (badfill)
    printf ("    %d input records rejected for missing values exceeding %d\n",
	badfill, max_miss);
  if (badpa)
    printf ("    %d input records rejected for PA more than %.2f from %.2f\n",
	badpa, dpa_max, pa_nom);
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
 *  12.05.23	added pangle selection options as for datavg
 *  v 0.9 frozen 2012.08.03
 *  13.06.10	general cleanup of unnecessary and sometimes interfering
 *		code for determining ephemeris from time
 *  v 1.0 frozen 2013.11.21
 *  15.02.24	initialized tlast; removed needless (and uninitialized)
 *		interpolation code and verbose messages of same
 *  v 1.1 frozen 2015.04.08
 *  15.09.08	use template record rather than first in series for keyword
 *	verification;
 *		check for existence and slotting of trec_key;
 *		use slotted query if possible;
 *		use drms_open_recordset instead of drms_open_records to
 *	eliminate length restriction, and drms_getvector to get record count;
 *	also results in substantial speed increas (~80%);
 *		changed default value for max_miss from 0 to "All" (i.e.
 *	undefined);
 *		fixed check for multiple values of T_OBS (which is probably
 *	obsolete) to ignore missing values
 *		print verbose diagnostic times in TAI
 *		altered order of detail reports to match order of tests
 *		removed needless interpolation code and verbose messages of
 *	same (only partially removed in v 1.1)
 *  15.09.14	fixed bug in slot index calculation
 *  v 1.2 frozen 2018.03.02
 *  17.01.25	added trap for mistaken or missing specifications of time
 *	targets when required; included function key_params_from_dspec, only
 *	function from selstuff used, and included earth_ephem
 *  18.05.04	removed unused include file
 *  v 1.3 frozen 2018.05.04
 *
 */
