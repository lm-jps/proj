/*
 *  xtrackd.c						~rick/hmi/rings/src
 *
 *  Extract a subinterval from a tracked cube
 *
 *  This is a special case application of the more general module drms_rebin,
 *    with no actual rebinning, just extraction in the slow dimension only,
 *    but with special additional features to analyze and recreate the log
 *    and to set applicable keys.
 *     
 *  Responsible:  Rick Bogart			RBogart@solar.stanford.edu
 *
 *  Usage:
 *    xtrackd [[arg= ]...]
 *
 *  Parameters: (type   	default 	description)
 *	in	DataSet		-	input data set
 *	out	DataSeries	-	output data series
 *
 *  Bugs:
 *    No check for agreement of dimensions between input and output
 *    Log information not copied to output
 *    Only uses earth ephemeris
 *    Data are unconditionally treated as floats (except in set_stat_keys)
 *    Copying of LonHG and LatHG assumes Carrington tracking
 *    The following mtrack standard keys are not (yet) set:
 *	CarrTime, CarrLon, LonSpan, Input, *_bscale, *_bzero, MeanMU, MeanLat,
 *	MDI_PA_Corr, PrimeKeyString
 *
 *  Revision history is at end of file.
 */
						       /*  required include  */
#include <jsoc_main.h>
#include "keystuff.c"
#include "earth_ephem.c"
#include "soho_ephem.c"

char *module_name = "xtrackd";
char *module_desc = "extract subset from tracked region";
char *version_id = "0.7";

ModuleArgs_t module_args[] = {
  {ARG_DATASET,	"in", "", "input data set: tracked cubes"}, 
  {ARG_STRING,	"out", "", "output data series"}, 
  {ARG_INT,	"length", "", "length in timesteps", "(0,)"},
  {ARG_TIME,	"tmid", "mid-point",
      "midpoint of extraction interval, in either CR:CL or standard date_time format"},
  {ARG_FLOATS,	"mai", "NaN", "Magnetic Activity Indices"},
  {ARG_FLAG, "n", "", "do not save output (for testing)"},
  {ARG_FLAG, "v", "", "verbose mode"}, 
  {}
};

/*
 *  Find and return (the) appropriate output segment(s) in the input or
 *    output data series
 *  An appropriate series will have one 3-dimensional array segment
 *    and one generic segment for the log.
 *  If there is only one segment in the series it is assumed to hold the
 *    data cube
 *  If the series has more than one appropriate segment of each type, the
 *    "first" appropriate segment for the array will be returned, whatever
 *    that means.
 */
int check_seginfo (DRMS_Record_t *rec, int *axis, int *datasegnum,
    int *logsegnum) {
  DRMS_Segment_t *seg;
  static int need_check = 0;
  int segfound, status;
  int i, n;

  if (*datasegnum < 0 && *logsegnum < 0) {
    int segct= rec->segments.num_total;

    if (segct < 1) {
      fprintf (stderr, "Error: series %s lacks data segments\n",
	  rec->seriesinfo->seriesname);
      return -1;
    }
    segfound = 0;
    for (n = 0; n < segct; n++) {
      seg = drms_segment_lookupnum (rec, n);
      if (seg->info->protocol == DRMS_GENERIC) continue;
      if (seg->info->scope == DRMS_CONSTANT) continue;
      if (seg->info->naxis != 3) continue;
      if (!segfound) {
	*datasegnum = n;
	for (i = 0; i < 3; i++) axis[i] = seg->axis[i];
        if (seg->info->scope == DRMS_VARDIM) need_check = 1;
      }
      segfound++;
    }
    if (!segfound) {
      fprintf (stderr, "Error: series %s lacks 3-dimensional data segment\n",
	  rec->seriesinfo->seriesname);
      return -1;
    }
    if (segfound > 1) fprintf (stderr,
	"Warning: series %s contains more than one 3-d segment; using %s\n",
	rec->seriesinfo->seriesname, seg->info->name);
    segfound = 0;
    for (n = 0; n < segct; n++) {
      seg = drms_segment_lookupnum (rec, n);
      if (seg->info->scope == DRMS_CONSTANT) continue;
      if (seg->info->protocol != DRMS_GENERIC) continue;
      if (!segfound) *logsegnum = n;
      segfound++;
    }
    if (!segfound) {
      fprintf (stderr, "Warning: series %s lacks segment for log\n",
	  rec->seriesinfo->seriesname);
    }
    if (segfound > 1) fprintf (stderr,
	"Warning: series %s contains more than one generic segment; using %s\n",
	seg->record->seriesinfo->seriesname, seg->info->name);
    return 0;
  }

  if (need_check) {
    seg = drms_segment_lookupnum (rec, *datasegnum);
    for (i = 0; i < 3; i++) axis[i] = seg->axis[i];
  }
  return 0;
}
/*
double get_carrlon (TIME t) {
  TIME upd;
  double r, lat, lon, vr, vn, vw;;

  soho_ephemeris (t, &r, &lat, &lon, &vr, &vn, &vw, &upd);
  return lon;
}
*/
void analyze_log_file (char *fname, int start, int planes, float *coverage,
    int *firstvalid, int *lastvalid) {
  FILE *log = fopen (fname, "r");
  int counting = 0, ngood = 0, ntot = 0;
  int nrec, rdum, status;
  char *line = NULL;
  char dsdum[DRMS_MAXQUERYLEN];
  size_t len;

  if (!log) return;
  while (getline (&line, &len, log) > 0) {
    if (sscanf (line, "step %d ", &nrec) < 1) continue;
    if (nrec == start) counting = 1;
    if (counting) {
      status = sscanf (line, "step %d mapped from image %s", &rdum, dsdum);
      if (status == 2) {
        if (!ntot) *firstvalid = ntot;
	*lastvalid = ntot;
	ngood++;
      }
      ntot++;
      if (status == 0) {
	fprintf (stderr,
	    "Error: unexpected line in log file: coverage may be erroneous\n");
	fprintf (stderr, "  %s", line);
      }
    }
    if (ntot == planes) break;
  }
  *coverage = (double)ngood / (double)ntot;
  fclose (log);
}

int set_stat_keys (DRMS_Record_t *rec, DRMS_Array_t *data) {
  long long n, ntot, valid;
  double *dv;
  double vmn, vmx, sum, sum2, sum3, sum4, v2;
  double scal, avg, var, skew, kurt;
  float *fv;
  int i;
  int kstat = 0;

  ntot = (data->naxis) ? 1 : 0;
  for (i = 0; i < data->naxis; i++) ntot *= data->axis[i];
  if (ntot < 1) return 0;
  valid = ntot;
  vmn = 1.0 / 0.0;
  vmx = -vmn;
  sum = sum2 = sum3 = sum4 = 0.0;
  switch (data->type) {
    case DRMS_TYPE_DOUBLE:
      dv = (double *)data->data;
      for (n = 0; n < ntot; n++) {
        if (isfinite (dv[n])) {
	  if (dv[n] > vmx) vmx = dv[n];
	  if (dv[n] < vmn) vmn = dv[n];
	  sum += dv[n];
	  v2 = dv[n] * dv[n];
	  sum2 += v2;
	  sum3 += v2 * dv[n];
	  sum4 += v2 * v2;
	} else {
	  valid--;
	}
      }
      break;
    case DRMS_TYPE_FLOAT:
      fv = (float *)data->data;
      for (n = 0; n < ntot; n++) {
        if (isfinite (fv[n])) {
	  if (fv[n] > vmx) vmx = fv[n];
	  if (fv[n] < vmn) vmn = fv[n];
	  sum += fv[n];
	  v2 = fv[n] * fv[n];
	  sum2 += v2;
	  sum3 += v2 * fv[n];
	  sum4 += v2 * v2;
	} else {
	  valid--;
	}
      }
      break;
    default:
      return 0;
  }

  kstat += check_and_set_key_int (rec, "DATAVALS", valid);
  kstat += check_and_set_key_int (rec, "MISSVALS", ntot - valid);
  if (valid <= 0) return kstat;

  kstat += check_and_set_key_double (rec, "DATAMIN", vmn);
  kstat += check_and_set_key_double (rec, "DATAMAX", vmx);

  scal = 1.0 / valid;
  avg = scal * sum;
  var = scal * sum2 - avg * avg;
  skew = scal * sum3 - 3 * var * avg - avg * avg * avg;
  kurt = scal * sum4 - 4 * skew * avg - 6 * avg * avg * var -
      avg * avg * avg * avg;
  kstat += check_and_set_key_double (rec, "DATAMEAN", avg);
  if (var < 0.0) return kstat;
  kstat += check_and_set_key_double (rec, "DATARMS", sqrt (var));
  if (var == 0.0) return kstat;
  skew /= var * sqrt (var);
  kurt /= var * var;
  kurt -= 3.0;
  kstat += check_and_set_key_double (rec, "DATASKEW", skew);
  kstat += check_and_set_key_double (rec, "DATAKURT", kurt);

  return kstat;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *drs;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg, *ilog, *olog;
  DRMS_Array_t *data;
  TIME ttmid, ttstrt, ttstop;
  double tmid_cl, rsunm;
  float *mai;
  float coverage, pmid, tstep;
  float merid_v, zonal_v, lon_cm;
  int axis[3], start[3], stop[3];
  int tmid_cr;
  int pfirst, plast;
  int i, n, reci, reco, recs;
  int idatasegnum = -1, ilogsegnum = -1, odatasegnum = -1, ologsegnum = -1;
  int kstat, status;
  int need_stats, setmais;
  char pathname[2*DRMS_MAXPATHLEN+1];
  char source[DRMS_MAXQUERYLEN];
  char key[64], tbuf[64];
  char module_ident[64];

  double degrad = 180.0 / M_PI;
  int tfrominp = 0;

  char *inset = strdup (params_get_str (params, "in"));
  char *oser = strdup (params_get_str (params, "out"));
  TIME tmid = params_get_time (params, "time");
  int lngth = params_get_int (params, "length");
  int maict = params_get_int (params, "mai_nvals");
  int no_save = params_isflagset (params, "n");
  int verbose = params_isflagset (params, "v");
  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s: JSOC version %s\n", module_ident, jsoc_version);
  if (time_is_invalid (tmid)) {
		  /* requested time not in standard format YYYY.* or [M]JD_*  */
    char *tstr = strdup (params_get_str (params, "tmid"));
    if (sscanf (tstr, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
      tmid = earth_meridian_crossing (tmid_cl, tmid_cr);
    } else {
				       /* requested time not in format CR:CL  */
      fprintf (stderr,
	 "tmid not in standard format; using midtime from input record(s)\n");
      tfrominp = 1;
    }
  } else {
    double rsun, vr, vn, vw, img_lat;
    earth_ephemeris (tmid, &rsun, &img_lat, &tmid_cl, &vr, &vn, &vw);
    tmid_cr = carrington_rots (tmid, 0);
  }
					   /*  open and check input data set  */
  if (!(drs = drms_open_records (drms_env, inset, &status))) {
    fprintf (stderr, "Error (%s): unable to open input data set %s\n",
	module_ident, inset);
    fprintf (stderr, "       status = %d\n", status);
    return 1;
  }
  recs = drs->n;
  setmais = (maict == recs);
  if (setmais) {
    mai = (float *)malloc (maict * sizeof (float));
    for (i = 0; i < maict; i++) {
      snprintf (key, 64, "mai_%d_value", i);
      mai[i] = params_get_float (params, key);
    }
  }

  for (reci = 0; reci < recs; reci++) {
    irec = drs->records[reci];
    drms_sprint_rec_query (source, irec);
    if (verbose) printf ("processing input record %s\n", source);
    if (check_seginfo (irec, axis, &idatasegnum, &ilogsegnum)) {
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
    if (!(iseg = drms_segment_lookupnum (irec, idatasegnum))) {
      fprintf (stderr, "Error: could not find data segment in record:\n");
      fprintf (stderr, "       %s\n", source);
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
    for (n = 0; n < 2; n++) {
      start[n] = 0;
      stop[n] = axis[n] - 1;
    }
			    /*  the input looks okay, open the output record  */
    orec = drms_create_record (drms_env, oser, DRMS_PERMANENT, &status);
    if (status) {
      fprintf (stderr, "Error: drms_create_record returned %d for data series:\n",
	  status);
      fprintf (stderr, "       %s\n", oser);
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
    if (check_seginfo (orec, axis, &odatasegnum, &ologsegnum)) {
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
    if (!(oseg = drms_segment_lookupnum (orec, odatasegnum))) {
      fprintf (stderr, "Error: could not open output data segment\n");
      drms_close_records (drs, DRMS_FREE_RECORD);
      drms_close_record (orec, DRMS_FREE_RECORD);
      return 1;
    }
    if (reci == 0) need_stats = (drms_keyword_lookup (orec, "DATAMIN", 1) ||
	drms_keyword_lookup (orec, "DATAMAX", 1) ||
	drms_keyword_lookup (orec, "DATAMEAN", 1) ||
	drms_keyword_lookup (orec, "DATARMS", 1) ||
	drms_keyword_lookup (orec, "DATASKEW", 1) ||
	drms_keyword_lookup (orec, "DATAKURT", 1) ||
	drms_keyword_lookup (orec, "DATAVALS", 1) ||
	drms_keyword_lookup (orec, "MISSVALS", 1));

    ttmid = (TIME)drms_getkey_double (irec, "CRVAL3", &status);
    tstep = drms_getkey_float (irec, "CDELT3", &status);
    ttstrt = drms_getkey_time (irec, "T_START", &status);
    ttstop = drms_getkey_time (irec, "T_STOP", &status);
    if (time_is_invalid (ttmid)) ttmid = 0.5 * (ttstrt = ttstop);
    if (tfrominp) {
      double rsun, vr, vn, vw, img_lat;
/*
      tmid_cr = drms_getkey_int (irec, "CarrRot", &status);
      tmid_cl = drms_getkey_double (irec, "CMLon", &status);
      tmid = earth_meridian_crossing (tmid_cl, tmid_cr);
sprint_time (tbuf, tmid, "TAI", 0);
printf ("tmid = %d:%05.1f = %s\n", tmid_cr, tmid_cl, tbuf);
*/
      tmid = ttmid;
      earth_ephemeris (tmid, &rsun, &img_lat, &tmid_cl, &vr, &vn, &vw);
      tmid_cr = carrington_rots (tmid, 0);
/*
sprint_time (tbuf, tmid, "TAI", 0);
printf ("ttmid = %s\n", tbuf);
*/
    }
    pmid = drms_getkey_float (irec, "CRPIX3", &status);
    start[2] = pmid - 0.5 * lngth + 0.5;
    stop[2] = start[2] + lngth - 1;
    start[2] += (tmid - ttmid) / tstep;
    stop[2] += (tmid - ttmid) / tstep;
    if (start[2] < 0 || stop[2] >= iseg->axis[2]) {
      fprintf (stderr,
	  "Error: output interval not contained within input cube:\n");
      fprintf (stderr, "       [%d, %d] outside [0, %d]\n", start[2], stop[2],
	  iseg->axis[2] - 1);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
/*
fprintf (stderr, "[%d, %d] inside [0, %d]\n", start[2], stop[2], iseg->axis[2] - 1);
*/
    data = drms_segment_readslice (iseg, DRMS_TYPE_FLOAT, start, stop, &status);
    status = drms_segment_write (oseg, data, 0);
    if (ilogsegnum >= 0) {
      ilog = drms_segment_lookupnum (irec, ilogsegnum);
      drms_record_directory (irec, pathname, 1);
      strcat (pathname, "/");
      strncat (pathname, ilog->filename, DRMS_MAXPATHLEN);
      analyze_log_file (pathname, start[2], lngth, &coverage, &pfirst, &plast);
    }

    kstat = 0;
						       /*  likely prime keys  */
    kstat += check_and_copy_key (orec, irec, "Ident");
    kstat += check_and_copy_key (orec, irec, "LonHG");
    kstat += check_and_copy_key (orec, irec, "LatHG");
    kstat += check_and_set_key_int   (orec, "CarrRot", tmid_cr);
    kstat += check_and_set_key_float (orec, "CMLon", tmid_cl);
    lon_cm = drms_getkey_float (irec, "LonHG", &status); - tmid_cl;
    while (lon_cm < -180.0) lon_cm += 360.0;
    while (lon_cm > 180.0) lon_cm -= 360.0;
    kstat += check_and_set_key_float (orec, "LonCM", lon_cm);
    kstat += check_and_set_key_time (orec, "MidTime", tmid);
    kstat += check_and_set_key_time (orec, "T_START", tmid - 0.5 * lngth * tstep);
    kstat += check_and_set_key_time (orec, "T_STOP", tmid + 0.5 * lngth * tstep);
    kstat += check_and_set_key_float (orec, "Duration", lngth * tstep);
					      /*  miscellaneous general keys  */
    kstat += check_and_set_key_str   (orec, "Module", module_ident);
    kstat += check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str (orec, "Source", source);
								/*  WCS keys  */
    kstat += check_and_copy_key (orec, irec, "CRPIX1");
    kstat += check_and_copy_key (orec, irec, "CRPIX2");
    kstat += check_and_set_key_float (orec, "CRPIX3", pmid);
    kstat += check_and_copy_key (orec, irec, "CRVAL1");
    kstat += check_and_copy_key (orec, irec, "CRVAL2");
    kstat += check_and_set_key_double (orec, "CRVAL3", tmid);
    kstat += check_and_copy_key (orec, irec, "CDELT1");
    kstat += check_and_copy_key (orec, irec, "CDELT2");
    kstat += check_and_copy_key (orec, irec, "CDELT3");
    kstat += check_and_copy_key (orec, irec, "CUNIT1");
    kstat += check_and_copy_key (orec, irec, "CUNIT2");
    kstat += check_and_copy_key (orec, irec, "CUNIT3");
    kstat += check_and_copy_key (orec, irec, "CTYPE1");
    kstat += check_and_copy_key (orec, irec, "CTYPE2");
    kstat += check_and_copy_key (orec, irec, "CTYPE3");
    kstat += check_and_copy_key (orec, irec, "WCSNAME");
    kstat += check_and_copy_key (orec, irec, "WCSAXES");
							      /*  other keys  */
    kstat += check_and_set_key_float (orec, "Coverage", coverage);
    kstat += check_and_copy_key (orec, irec, "Width");
    kstat += check_and_copy_key (orec, irec, "Height");
    kstat += check_and_copy_key (orec, irec, "Size");
    kstat += check_and_copy_key (orec, irec, "Backgrnd");
    kstat += check_and_copy_key (orec, irec, "RejectList");
    kstat += check_and_copy_key (orec, irec, "Cadence");
    kstat += check_and_copy_key (orec, irec, "MapScale");
    kstat += check_and_copy_key (orec, irec, "Map_PA");
    kstat += check_and_copy_key (orec, irec, "MapProj");
    kstat += check_and_copy_key (orec, irec, "Interp");
    kstat += check_and_copy_key (orec, irec, "RSunRef");
    kstat += check_and_copy_key (orec, irec, "ZonalTrk");
    kstat += check_and_copy_key (orec, irec, "ZonalVel");
    kstat += check_and_copy_key (orec, irec, "MeridTrk");
    kstat += check_and_copy_key (orec, irec, "MeridVel");
    merid_v = drms_getkey_float (irec, "MeridTrk", &status);
    kstat += check_and_set_key_float (orec, "LatDrift",
	merid_v * degrad * lngth * tstep);
    zonal_v = drms_getkey_float (irec, "ZonalTrk", &status);
    kstat += check_and_set_key_float (orec, "LonDrift",
	zonal_v * degrad * lngth * tstep);
    if (setmais && isfinite (mai[reci]))
      kstat += check_and_set_key_float (orec, "MAI", mai[reci]);

    if (need_stats) kstat += set_stat_keys (orec, data);

    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s\n",
	  oser);
      fprintf (stderr, "      series may not have appropriate structure\n");
      drms_close_record (orec, DRMS_FREE_RECORD);
    }

    drms_close_record (orec, dispose);
    drms_free_array (data);
  }
/*
  drms_segment_getscaling (iseg, &bzero, &bscale);
bzero = 0.0;
bscale = 1.0;
  if (iseg->axis[0] % spbin || iseg->axis[1] % spbin) {
    fprintf (stderr,
	"Error: spatial binning parameter does not exactly divide one or more axis lengths\n");
    drms_close_records (irecs, DRMS_FREE_RECORD);
    return 0;
  }
/*
  if (time_is_invalid (tmid)) {
		 /* requested time not in standard format YYYY.* or [M]JD_*  */
/*
    char *tstr = params_get_str (params, "time");
    if (sscanf (tstr, "%d:%lf", &cr, &cl) != 2) {
				      /* requested time not in format CR:CL  */
/*
      if (sscanf (tstr, "CT%lf", &cl) != 1) {
	tstart = drms_getkey_time (irec, "start", &status);
	tstop = drms_getkey_time (irec, "stop", &status);
	tmid = ttmid = 0.5 * (tstart + tstop);
	if (verbose) {
	  sprint_time (tbuf, tmid, "TAI", 0);
	  printf ("tracking interval centered at %s\n", tbuf);
	}
	  /*  this will exactly center the extract if iplanes%2 = oplanes%2  */
/*
	pmid = iplanes / 2;
	pmidset = 1;
	pstart = pmid - (oplanes / 2);
	cl = carrington_rots (tmid, 1);
	cr = cl;
	cl = 360.0 * (1.0 + cr - cl);
      } else {
						    /*  convert CT to CR:CL  */
/*
        cr = cl;
	cl = 360.0 * (1.0 + cr - cl);
	tmid = SOHO_meridian_crossing (cl, cr);
      }
    } else tmid = SOHO_meridian_crossing (cl, cr);
  } else {
		  /*  get central meridian longitude for requested mid-time  */
/*
    cl = carrington_rots (tmid, 1);
    cr = cl;
    cl = 360.0 * (1.0 + cr - cl);
  }
  dcr = crsel - cr;
  while (dcr > 0) {
    cl += 360.0;
    dcr--;
  }
  while (dcr < 0) {
    cl -= 360.0;
    dcr++;
  }
  cmlon = lonsel - cl;
  tstart = drms_getkey_time (irec, "start", &status);
  tstop = drms_getkey_time (irec, "stop", &status);
  cadence = drms_getkey_float (irec, "Cadence", &status);
  if (!pmidset) {
    ttmid = 0.5 * (tstart + tstop);
    pmid = iplanes / 2;
    pmid += (tmid - ttmid) / cadence;
  }
  pxmid = oplanes / 2;
  pstart = pmid + 1 - (oplanes / 2) - (oplanes % 2);
  if (length % 2) {
    trstart = tmid - pxmid * cadence;
    trstop = tmid + pxmid * cadence;
  } else {
    trstart = tmid - (pxmid - 0.5) * cadence;
    trstop = tmid + (pxmid - 0.5) * cadence;
  }
  cl0 = get_carrlon (trstart);
  cl1 = get_carrlon (trstop);
  lonspan = cl0 - cl1;
  while (lonspan < 0.0) lonspan += 360.0;
  while (lonspan >= 360.00) lonspan -= 360.0;
    if (drms_segment_setscaling (oseg, 0.0, 0.5))
      fprintf (stderr, "Warning: scaling failed\n");
*/
/*
    unlink (tempname);
  }
  stat (pathname, &stat_buf);
  drms_setkey_longlong (orec, "Size", stat_buf.st_size);

  drms_close_record (orec, DRMS_INSERT_RECORD);
*/
  drms_close_records (drs, DRMS_FREE_RECORD);
  return 0;
}

/*
 *  Revision History
 *
 *  12.03.06	created this file, based on xxtrackd, stripping out CGI and
 *		export related stuff and old libraries
 */
