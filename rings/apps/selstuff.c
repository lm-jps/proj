/*
 *  selstuff.c				(to be linked from) ~rick/src/util
 *
 *  library of miscellaneous utility functions for extracting a target data
 *    set from a series by selecting on parameters, and for extracting
 *    parameter values from a dataset selection
 *
 *  key_params_from_dspec()
 *  select_dataset_from_time_interval()
 *  select_dataset_from_time_range()
 *
 *  Bugs:
 *    Lots of checking is skipped
 *    select_dataset_from_time_interval() should just call
 *	select_dataset_from_time_range(), or vice-versa
 *    The record selection based on Carrington rotation and longitude gives
 *	incomplete or screwy results sometimes, especially with HMI data
 *    Unconditionally uses geocentric ephemeris for Carrington times
 *
 *  Revision history is at end of file
 */
 
#include "earth_ephem.c"

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

DRMS_RecordSet_t *select_dataset_from_time_range (const char *series,
    char *tstrt_str, char *tstop_str) {
/*
 *  Select a dataset from a series based on specification of target times
 *    for start and stop
 *  The target times may be specified as either date_time strings or as CR:CL
 */
  DRMS_RecordSet_t *ds = NULL;
  TIME tstrt, tstop;
  double strt_cl, stop_cl;
  int strt_cr, stop_cr;
  int status;
  char rec_query[DRMS_MAXQUERYLEN];
		   /*  get required series info from first record in series  */
					       /*  platform, cadence, phase  */
/*  not currently implemented  */
/*
  snprintf (rec_query, DRMS_MAXQUERYLEN, "%s[#^]", series);
  if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
    fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
    fprintf (stderr, "       status = %d\n", status);
  }
  rec = ds->records[0];
  drms_close_records (ds, DRMS_FREE_RECORD);
*/
  if (sscanf (tstrt_str, "%d:%lf", &strt_cr, &strt_cl) == 2) {
			 /*  tstrt specified as CR:CL : need ephemeris info  */
    if (sscanf (tstop_str, "%d:%lf", &stop_cr, &stop_cl) != 2) {
      fprintf (stderr, "Error: start and stop times must be specified in same format\n");
      fprintf (stderr, "       either Date_Time or CR:CL\n");
      return ds;
    }
/*
    if (strt_cr == stop_cr) {
      snprintf (rec_query, DRMS_MAXQUERYLEN,
	  "%s[?CAR_ROT = %d and CRLN_OBS <= %f and CRLN_OBS >= %f?]", series,
	  strt_cr, strt_cl, stop_cl);
    } else if ((stop_cr - strt_cr) > 1) {
      snprintf (rec_query, DRMS_MAXQUERYLEN,
	  "%s[?(%s > %d and %s < %d) or (%s = %d and %s <= %f) or (%s = %d and %s >= %f?)]",
	  series, "CAR_ROT", strt_cr, "CAR_ROT", stop_cr,
	  "CAR_ROT", strt_cr, "CRLN_OBS", strt_cl,
	  "CAR_ROT", stop_cr, "CRLN_OBS", stop_cl);
    } else {
      snprintf (rec_query, DRMS_MAXQUERYLEN,
	  "%s[?(%s = %d and %s <= %f) or (%s = %d and %s >= %f)?]",
	  series, "CAR_ROT", strt_cr, "CRLN_OBS", strt_cl,
	  "CAR_ROT", stop_cr, "CRLN_OBS", stop_cl);
    }
*/
    tstrt = earth_meridian_crossing (strt_cl, strt_cr);
    tstop = earth_meridian_crossing (stop_cl, stop_cr);
  } else {
    tstrt = sscan_time (tstrt_str);
    tstop = sscan_time (tstop_str);
    if (time_is_invalid (tstrt) || time_is_invalid (tstop)) {
      fprintf (stderr, "Error: start and stop times must be specified in same format\n");
      fprintf (stderr, "       either Date_Time or CR:CL\n");
      return ds;
    }
  }
  snprintf (rec_query, 256, "%s[?%s >= %13.6e and %s <= %13.6e?]", series,
/*
      trec_key, tstrt - t_eps, trec_key, tstop + t_eps);
*/
      "T_REC", tstrt, "T_REC", tstop);
  if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
    fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
    fprintf (stderr, "       status = %d\n", status);
  }
  return ds;
}

DRMS_RecordSet_t *select_dataset_from_time_interval (const char *series,
    char *tmid_str, double intrvl) {
/*
 *  Select a dataset from a series based on specification of a target time
 *    and interval length
 *  The target time may be specified as either a date_time string or as CR:CL
 *    In the former case, the length specification is expected to be in seconds;
 *    in the latter case it is expected to be in degrees of Carrington rotation
 *  NO, this doesn't yet work; need to specify intrvl in seconds
 */
  DRMS_RecordSet_t *ds = NULL;
  DRMS_Record_t *rec;
  TIME tmid, tstrt, tstop;
  double tmid_cl, tstrt_cl, tstop_cl;
  int tmid_cr, tstrt_cr, tstop_cr;
  int status;
  char rec_query[DRMS_MAXQUERYLEN];
		   /*  get required series info from first record in series  */
					       /*  platform, cadence, phase  */
/*  not currently implemented  */
/*
  snprintf (rec_query, DRMS_MAXQUERYLEN, "%s[#^]", series);
  if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
    fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
    fprintf (stderr, "       status = %d\n", status);
  }
  rec = ds->records[0];
  drms_close_records (ds, DRMS_FREE_RECORD);
*/

  if (sscanf (tmid_str, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
						/*  tmid specified as CR:CL  */
    tmid = earth_meridian_crossing (tmid_cl, tmid_cr);
    tstrt_cr = tstop_cr = tmid_cr;
    tstrt_cl = tmid_cl + 0.5 * intrvl;
    while (tstrt_cl > 360.0) {
      tstrt_cl -= 360.0;
      tstrt_cr--;
    }
    tstop_cl = tmid_cl - 0.5 * intrvl;
    while (tstop_cl < 0.0) {
      tstop_cl += 360.0;
      tstop_cr++;
    }
    tstrt = earth_meridian_crossing (tstrt_cl, tstrt_cr);
    tstop = earth_meridian_crossing (tstop_cl, tstop_cr);
/*
    if (stop_cl >= 0.0 && strt_cl <= 360.0) {
      snprintf (rec_query, DRMS_MAXQUERYLEN,
	  "%s[?CAR_ROT = %d and CRLN_OBS <= %f and CRLN_OBS >= %f?]", series,
	  tmid_cr, strt_cl, stop_cl);
    } else {
      int strt_cr, stop_cr;
      strt_cr = stop_cr = tmid_cr;
    while (strt_cl > 360.0) {
      strt_cl -= 360.0;
      strt_cr--;
    }
    while (stop_cl < 0.0) {
      stop_cl += 360.0;
      stop_cr++;
    }
      if ((stop_cr - strt_cr) > 1) {
	snprintf (rec_query, DRMS_MAXQUERYLEN,
	    "%s[?(%s > %d and %s < %d) or (%s = %d and %s <= %f) or (%s = %d and %s >= %f?)]",
	    series, "CAR_ROT", strt_cr, "CAR_ROT", stop_cr,
	    "CAR_ROT", strt_cr, "CRLN_OBS", strt_cl,
	    "CAR_ROT", stop_cr, "CRLN_OBS", stop_cl);
      } else {
	snprintf (rec_query, DRMS_MAXQUERYLEN,
	    "%s[?(%s = %d and %s <= %f) or (%s = %d and %s >= %f)?]",
	    series, "CAR_ROT", strt_cr, "CRLN_OBS", strt_cl,
	    "CAR_ROT", stop_cr, "CRLN_OBS", stop_cl);
      }
    }
*/
  } else {
			      /*  tmid specified as normal date-time string  */
    tmid = sscan_time (tmid_str);
    tstrt = tmid - 0.5 * intrvl;
    tstop = tmid + 0.5 * intrvl;
  }
  snprintf (rec_query, 256, "%s[?%s >= %17.10e and %s <= %17.10e?]", series,
/*
      trec_key, tstrt - t_eps, trec_key, tstop + t_eps);
*/
      "T_REC", tstrt, "T_REC", tstop);
  if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
    fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
    fprintf (stderr, "       status = %d\n", status);
  }
  return ds;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  10.04.27		created this file from routines already in mtrack and
 *		rdcover
 *  10.05.03		added function select_dataset_from_time_range
 *  10.08.19		extended precision of time query in
 *		select_dataset_from_time_interval
 *
 */
