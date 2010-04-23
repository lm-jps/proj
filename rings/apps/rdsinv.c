/*
 * rdsinv.c				~baldner/src/rdsinv, ~rick/src/rdsinv
 *
 *  Asymptotic sound speed inversion code for ring diagrams
 *
 *  Responsible:
 * 	Charles Baldner, Rick Bogart, Sarbani Basu
 *
 *  Usage:
 *    rdsinv [-v] [in= seg= comp= out= ] ...
 *
 *  Arguments: (type	default         description)
 *    in	string	su_rsb.rdfits	Input data series name, or dataset
 *	specification, or filename; if only the series name is provided, 
 *	additional arguments are required to specify a selection range.
 *	If it cannot be opened as a data set, the argument is interpreted
 *	as the filename of a single ringfit file
 *    seg	string	fit.out
 *    cr	int	NS
 *    clon	float	NaN
 *    lon	float	NaN
 *    lat	float	NaN
 *    comp	string	yale_cb.rdVfitsc_avg
 *    out	string	/tmp		Output data series name, or directory
 *	name. If a record cannot be opened in the data series, the argument is
 *	interpreted as a directory name.
 *    no	int	6		Number of knots in frequency
 *    nw	int	10		Number of knots in w (acoustic depth)
 *    nsol	int	50		Number of solution points
 *
 *  Flags
 *	-n	do not save results (diagnostic only)
 *
 *  Bugs:
 *    The specification of any of (cr, clon, lon, lat) will force the
 *	interpretation of the in argument name as a data series or data set
 *
 *    The default values for in, seg, comp, and especially out are not
 *	appropriate for general use
 *
 *    The seg argument should not be necessary; if multiple segments are
 *	present in the input records, they should be specified as part of
 *	the data set name. In any case, the segment name for the comparison
 *	series is currently hard-coded.
 *
 *    Numerous combinations of requirements for distinctions among
 *	"prime keys" in output file names are not supported; it's not clear
 *	if it's such a good idea to try to distinguish, may be just better to
 *	overspecify.
 *
 *  
 *  Revision history is at end of file
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jsoc_main.h>
#include "keystuff.c"
#include "rdutil.c"
							     /*  prototypes  */
extern void h1h2_(int *n, int *l, double *om, double *dom,
    double *err, double *w, int *nmode, int *no, int *nw,
    int *nsol, char *filcsq, char *filsw, char *filout,
    int *qverb, int csqlen, int swlen, int outlen);
int get_comp_filepath(double lat, double lon, char *comp, char *filepath, char *compKeyWord);

char *module_name = "rdsinv";
char *version_id = "0.7";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in", "su_rsb.rdfits", "input data set (or file)"},
  {ARG_STRING, "seg", "fit.out", "data segment name"},
  {ARG_INT,    "cr", "Not Specified", "Carrington rotation for all regions"},
  {ARG_FLOAT,  "clon", "Not Specified", "Carrington longitude of CM for all regions"},
  {ARG_FLOAT,  "lon", "Not Specified", "Carrington longitude of all regions"},
  {ARG_FLOAT,  "lat", "Not Specified", "Carrington latitude of all regions"},
  {ARG_INT,    "no", "6", "number of knots in frequency"},
  {ARG_INT,    "nw", "10", "number of knots in w (accoustic depth)"},
  {ARG_INT,    "nsol", "50", "number of solution points"},
//  {ARG_STRING, "comp", "Not Specified", "Comparison frequencies"},
  {ARG_STRING, "comp", "yale_cb.rdVfitsc_avg", "Comparison frequencies"},
  {ARG_STRING, "out", "/tmp", "Output series name or directory name"},
  {ARG_FLAG,    "v",    "", "run in verbose mode"},
  {ARG_END}
};

       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "MapProj", "MapScale", "Map_PA", "Width", "Height",
    "Size", "Cadence", "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel", "MAI"};

int get_comp_filepath(double lat, double lon, char *comp, char *filepath, char *compKeyWord)	{
  int i, j, status, numrec, car, userec;
  DRMS_RecordSet_t *recordSet = NULL;
  DRMS_Record_t *record;
  DRMS_Segment_t *segment;
  double distance, latC, lonC, minDistance = 1.0e20;
  char *carTime;
		 /*  first try opening comp as a file; if it is one, return  */
  FILE *test = fopen (comp, "r");
  if (test) {
    fclose (test);
    sprintf (filepath, "%s", comp);
    return 0;
  }

  recordSet = drms_open_records (drms_env, comp, &status);
  if (status) {
    fprintf (stderr, "Error: unable to open %s\n", comp);
    fprintf (stderr, "       as either a data set or a file\n");
    return 1;
  }
  
  numrec = recordSet->n;
//    distance = (double*) malloc(numrec*sizeof(double));
  for (i=0; i<numrec; i++)	{
    record = recordSet->records[i];
    latC = drms_getkey_double(record, "LatHG", &status);
    if (status) {
	fprintf(stderr, "Keyword LatHG not found in comparison series record\n");
	return 1;
    }
      lonC = drms_getkey_double(record, "LonCM", &status);
      if (status) {
	lonC = drms_getkey_double(record, "LonHG", &status);
	if (status) {
	  fprintf (stderr, "Error: no keyword for longitude found in comparison series record\n");
	  return 1;
	}
      }

      car = drms_getkey_int(record, "CarrRot", &status);
      if (status) {
	carTime = drms_getkey_string(record, "CarrTime", &status);
	if(status)	{
	  fprintf (stderr, "Keyword CarrTime not found in comparison series record\n");
	  return 1;
	}
	sscanf(carTime, "%*i:%i", &car);
      }

      //lonC = ((double) car) - lonC;
      distance = (lon-lonC)*(lon-lonC)+(lat-latC)*(lat-latC);
      if(distance < minDistance)	{
	minDistance = distance;
	userec = i;
	sprintf(compKeyWord, "%s[%4.1lf][%4.1lf][%4i]", comp, latC, lonC, (int) car);
      }
    }
    record = recordSet->records[userec];
    if (record->sunum != -1LL && record->su == NULL)  {
      record->su = drms_getunit (record->env, record->seriesinfo->seriesname,
	  record->sunum, 1, &status);
    }
  segment = drms_segment_lookup (record, "fit.out");
  sprintf (filepath, "%s/S%05i/%s", segment->record->su->sudir, segment->record->slotnum, "fit.out");
  fprintf (stderr, "Using comparison region %s\n", compKeyWord);

  return 0;
}
		/*  these functions should probably be added to keystuff.c  */
int key_is_prime (DRMS_Record_t *rec, const char *key) {
  DRMS_Keyword_t *pkey, *keywd;
  int n, npct = rec->seriesinfo->pidx_num;

  if (!(keywd = drms_keyword_lookup (rec, key, 1))) return 0;
  if (npct) {
    for (n = 0; n < npct; n++) {
      pkey = rec->seriesinfo->pidx_keywords[n];
      if (pkey->info->recscope > 1) pkey = drms_keyword_slotfromindex (pkey);
      if (strcmp (key, pkey->info->name)) continue;
      return 1;
    }
  }
  return 0;
}

int unique_values (const char *dsqry, const char *key, DRMS_Record_t *rec) {
  int status, values;
  DRMS_Array_t *value;
  DRMS_Keyword_t *pkey, *keywd;

  if (!(keywd = drms_keyword_lookup (rec, key, 1))) return 0;
  value = drms_record_getvector (drms_env, dsqry, key,
      DRMS_TYPE_FLOAT, 1, &status);
  if (!value) return 0;
  values = value->axis[1];
  drms_free_array (value);
  return values;
}

char *printcrcl_loc (int cr, double cl, double lon, double lat) {
  static char fnam[32];
  double dlon;
						/*  construct a filename CR:CL_lonEWlatNS  */
  while (lon < 0) lon += 360;
  while (lon > 360) lon -= 360;
  if (cr > 0) {
    while (cl < 0) {
      cl += 360;
      cr++;
    } while (cl > 360) {
      cl -= 360;
      cr--;
    }
    dlon = lon - cl;
    while (dlon > 180) dlon -= 360;
    while (dlon < -180) dlon += 360;
    if (dlon >= 0) {
      if (lat >= 0) sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fN", cr, cl, dlon, lat);
      else sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fS", cr, cl, dlon, -lat);
    } else {
      if (lat >= 0) sprintf (fnam, "%04d:%05.1f_%04.1fE%04.1fN", cr, cl, -dlon, lat);
      else sprintf (fnam, "%04d:%05.1f_%04.1fE%04.1fS", cr, cl, -dlon, -lat);
    }
  } else sprintf (fnam, "%05.1f%+04.1f", lon, lat);
  return fnam;
}

int DoIt(void)	{
  CmdParams_t *params = &cmdparams;
  int status = 0;
  DRMS_RecordSet_t *recordSet = NULL;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *segment, *osegment;
  FILE *fpt;
  // these next need to be arguments
char *filcsq = "/home/baldner/src/rdsinv/fort.33";
char *filsw = "/home/baldner/src/rdsinv/fw.out";
int outlen;
  int i, j, k, rec_i, drms_input, drms_output;
  int *n_in, *n, *n_i, *l, *l_i, *mask;
  int npts_in, npts_comp, npts, nrec, seg_ct;
  double *l_in, *om_in, *err_in, *om_i, *err_i, *om, *dom, *err, *w;
  char *carstr, *suffix;
  double latc, lonc, loncCar;
  char filename[DRMS_MAXPATHLEN], compfile[DRMS_MAXPATHLEN];
  char odir[DRMS_MAXPATHLEN], ofilename[DRMS_MAXPATHLEN];
  char compKeyWord[DRMS_MAXPATHLEN];
  char rec_query[DRMS_MAXQUERYLEN], source[DRMS_MAXQUERYLEN];
  char module_ident[64];

  char *in = strdup (params_get_str (params, "in"));
  char *seg = strdup (params_get_str (params, "seg"));
  char *out = strdup (params_get_str (params, "out"));
  char *comp = strdup (params_get_str (params, "comp"));
  int cr = params_get_int (params, "cr");
  float cl = params_get_float (params, "clon");
  float lon = params_get_float (params, "lon");
  float lat = params_get_float (params, "lat");
  int no = params_get_int (params, "no");
  int nw = params_get_int (params, "nw");
  int nsol = params_get_int (params, "nsol");
  int verbose = params_isflagset (params, "v");

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);

int csqlen = strlen(filcsq);
int swlen = strlen(filsw);

  int anycr = (cr <= 0);
  int anycl = isnan (cl);
  int anylon = isnan (lon);
  int anylat = isnan (lat);

  drms_input = 1;
  if (anycr && anycl && anylon && anylat) {
		/*  no target range specified, assume implicit in input set  */
			/*  but it might be implicit in output series spec!  */
    recordSet = drms_open_records (drms_env, in, &status);
    if (status) {
      fpt = fopen (in, "r");
      if (!fpt) {
	fprintf (stderr, "Error: in specification \"%s\"\n", in);
	fprintf (stderr, "       does not specify a readable data set nor file\n");
	return 1;
      }
      drms_input = 0;
      nrec = 1;
    } else nrec = recordSet->n;
    if (nrec < 1) {
      fprintf (stderr, "No data records in set %s\n", in);
      return 1;
    }
  } else {
    strncpy (rec_query, in, DRMS_MAXQUERYLEN);
    if (!anycr) {
      snprintf (source, DRMS_MAXQUERYLEN, "[?CarrRot=%d?]", cr);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
    }
    if (!anycl) {
      float clmin = cl - 0.001;
      float clmax = cl + 0.001;
						/*  need to deal with 0/360  */
      snprintf (source, DRMS_MAXQUERYLEN, "[?CMLon>%g and CMLon<%g?]",
	clmin, clmax);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
    }
    if (!anylon) {
      float lonmin = lon - 0.01;
      float lonmax = lon + 0.01;
						/*  need to deal with 0/360  */
      snprintf (source, DRMS_MAXQUERYLEN, "[?LonHG>%g and LonHG<%g?]",
	lonmin, lonmax);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
    }
    if (!anylat) {
      float latmin = lat - 0.01;
      float latmax = lat + 0.01;
      snprintf (source, DRMS_MAXQUERYLEN, "[?LatHG>%g and LatHG<%g?]",
	  latmin, latmax);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
    }
    if (!(recordSet = drms_open_records (drms_env, rec_query, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    nrec = recordSet->n;
  }
  printf ("Processing %d records\n", nrec);

  if (drms_input) {
		/*  check that range is restricted by dataset specification  */
    irec = recordSet->records[0];
    if (anycr) anycr = unique_values (in, "CarrRot", irec) - 1;
    if (anycl) anycl = unique_values (in, "CMLon", irec) - 1;
    if (anylon) anylon = unique_values (in, "LonHG", irec) - 1;
    if (anylat) anylat = unique_values (in, "LatHG", irec) - 1;
  }

  drms_output = 1;
  orec = drms_create_record (drms_env, out, DRMS_PERMANENT, &status);
  if (status) {
    drms_output = 0;
    fprintf (stderr,
        "Warning: drms_create_record() returned %d for data series:\n", status);
    fprintf (stderr, "       %s; will be interpreted as directory name\n", out);
    strncpy (odir, out, DRMS_MAXPATHLEN);
  } else if (drms_record_directory (orec, odir, 1)) {
    fprintf (stderr,
        "Error: drms_record_directory() failure: SUMS may not be available\n");
    fprintf (stderr,
        "       or series %s may contain no segments\n", out);
    return 1;
  }

  if (drms_output) {
    int kstat = 0;
    int keyct = sizeof (propagate) / sizeof (char *);
				    /*  copy designated keywords from input  */
/*
      kstat = propagate_keys (orec, irec, propagate, keyct);
*/
    kstat += check_and_set_key_str (orec, "Module", module_ident);
/*
    kstat += check_and_set_key_str (orec, "Source", source);
    kstat += check_and_set_key_str (orec, "Compare", compKeyWord);
*/
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
					       /*  module specific keywords  */
    kstat += check_and_set_key_int (orec, "NumO", no);
    kstat += check_and_set_key_int (orec, "NumW", nw);
    kstat += check_and_set_key_int (orec, "NumSol", nsol);
				 /*  set possible prime keys as appropriate  */
    if (anycr) {
      if (key_is_prime (orec, "CarrRot")) {
	fprintf (stderr, "Error: CarrRot is prime key for output series %s\n", out);
	fprintf (stderr, "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
				     /*  this pretty much has to be true...  */
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_int (orec, "CarrRot", cr);
    if (anycl) {
      if (key_is_prime (orec, "CMLon")) {
	fprintf (stderr, "Error: CMLon is prime key for output series %s\n", out);
	fprintf (stderr, "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_double (orec, "CMLon", cl);
    if (anycr || anycl) {
      if (key_is_prime (orec, "CarrTime")) {
	fprintf (stderr, "Error: CarrTime is prime key for output series %s\n", out);
	fprintf (stderr, "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else {
      char crcl[32];
      sprintf (crcl, "%d:%03.0f", cr, cl);
      kstat += check_and_set_key_str (orec, "CarrTime", crcl);
    }
    if (anylon) {
      if (key_is_prime (orec, "LonHG")) {
	fprintf (stderr, "Error: LonHG is prime key for output series %s\n", out);
	fprintf (stderr, "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_double (orec, "LonHG", lon);
    if (anylat) {
      if (key_is_prime (orec, "LatHG")) {
	fprintf (stderr, "Error: LatHG is prime key for output series %s\n", out);
	fprintf (stderr, "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_double (orec, "LatHG", lat);

    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s:\n",
	  out);
      fprintf (stderr, "      series may not have appropriate structure;\n");
      drms_close_record (orec, DRMS_FREE_RECORD);
      if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
      return 1;
    }
  }

  for (rec_i=0; rec_i<nrec; rec_i++) {
    if (verbose) printf ("  Processing record %i\n", rec_i);
    if (drms_input) {
      irec = recordSet->records[rec_i];
      if (irec->sunum != -1LL && irec->su == NULL) {
	irec->su = drms_getunit (irec->env, irec->seriesinfo->seriesname,
	    irec->sunum, 1, &status);
      }
      segment = drms_segment_lookup(irec, seg);
      sprintf (filename, "%s/S%05i/%s", segment->record->su->sudir,
	  segment->record->slotnum, seg);
      fpt = fopen (filename, "r");
      drms_sprint_rec_query (source, irec);
    }
    else {
      fpt = fopen (in, "r");
      strcpy (source, in);
    }
    if (!fpt) {
      fprintf (stderr, "Error: could not open SUMS file for record\n");
      fprintf (stderr, "       %s\n", source);
      drms_close_records (recordSet, DRMS_FREE_RECORD);
      if (drms_output) drms_close_record (orec, DRMS_FREE_RECORD);
      return 1;
    }
    if (verbose) printf ("  %s\n", source);
    read_fit (fpt, &npts_in, &n_in, &l_in, &om_in, &err_in);
    fclose (fpt);
    n_i = (int*) malloc(npts_in*sizeof(int));
    l_i = (int*) malloc(npts_in*sizeof(int));
    om_i = (double*) malloc(npts_in*sizeof(double));
    err_i = (double*) malloc(npts_in*sizeof(double));
    n = (int*) malloc(npts_in*sizeof(int));
    l = (int*) malloc(npts_in*sizeof(int));
    om = (double*) malloc(npts_in*sizeof(double));
    dom = (double*) malloc(npts_in*sizeof(double));
    err = (double*) malloc(npts_in*sizeof(double));
    w = (double*) malloc(npts_in*sizeof(double));
    interp_freq (n_in, l_in, om_in, err_in, npts_in, n_i, l_i, om_i, err_i,
	&npts_in);
    free(n_in);
    free(l_in);
    free(om_in);
    free(err_in);
    if (drms_input) {
      latc = drms_getkey_double (irec, "LatHG", &status);
      lonc = drms_getkey_double (irec, "LonHG", &status);
      carstr = drms_getkey_string (irec, "CarrTime", &status);
      if (status) {
        cr = drms_getkey_int (irec, "CarrRot", &status);
	if (status) cr = -1;
	loncCar = drms_getkey_double (irec, "CMLon", &status);
	if (status) loncCar = 2 * lonc;
      } else sscanf (carstr, "%i:%lf", &cr, &loncCar);
      lonc = loncCar - lonc;
    } else latc = lonc = 0.0;

    status = get_comp_filepath (latc, lonc, comp, compfile, compKeyWord);
//    fprintf(stderr, "status=%i\n",status);
//    fprintf(stderr, "Opening %s\n", compfile);
    fpt = fopen(compfile, "r");
    read_fit(fpt, &npts_comp, &n_in, &l_in, &om_in, &err_in);
    fclose(fpt);
    freqdif(n_i, n_in, l_i, l_in, om_i, om_in, err_i, err_in, npts_in, npts_comp,
	n, l, om, dom, err, &npts);
    mask = (int*) malloc(npts*sizeof(int));
    for(i=0; i<npts; i++) mask[i] = 0;
    autoweed(l, n, om, dom, err, mask, npts);
    j = 0;
//    printf("npts=%i\n", npts);
    for(i=0; i<npts; i++)	{
      if(mask[i])	{
	l[j] = l[i];
	n[j] = n[i];
	om[j] = om[i];
	dom[j] = dom[i];
	err[j] = err[i];
	w[j] = log10(om[i]/(l[i]+0.5)) - 3.0;
	j++;
      }
    }
    npts = j;
//    printf("npts=%i\n", npts);
    sprintf (ofilename, "%s/%s", odir, printcrcl_loc (cr, loncCar, lonc + loncCar, latc));
    if (verbose) printf ( "  Writing solution to %s\n", ofilename);
    outlen = strlen (ofilename);

    h1h2_(n, l, om, dom, err, w, &npts, &no, &nw, &nsol, filcsq, filsw, ofilename,
	&verbose, csqlen, swlen, outlen);

    free(n_in);
    free(l_in);
    free(om_in);
    free(err_in);
    free(n_i);
    free(l_i);
    free(om_i);
    free(err_i);
    free(n);
    free(l);
    free(om);
    free(dom);
    free(err);
    free(w);
    free(mask);
    if (verbose) printf ("  Closeout record %d...\n", rec_i);
  }
  if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
  if (drms_output) drms_close_record (orec, DRMS_INSERT_RECORD);
  return 0;
}

/*
 *  Revision History 
 *
 *   November 2009 - v 0.1 module created by Charles Baldner
 *   Dec 12 2009 - v 0.5
 *     Added ability to read comparison region from drms series
 *   Jan 21 2009 - v 0.6
 *     Updated to write CompReg keyword to output series
 *  v 0.7	fixed bug in input record close; simplified keyword propagation
 *		Removed separate argument for named file input
 *		Extensive changes to driver code to support writing of multiple
 *		output files into single output directory associated with single
 *		SUMS record
 *  v 0.7 frozen 2010.04.23
 *
 */
