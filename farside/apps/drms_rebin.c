/*
 *  drms_rebin.c						~rick/src/drms
 *
 *  Rebin the values of selected segments of selected records to a desired
 *    grid (and location)
 *
 *  Responsible:  Rick Bogart				RBogart@spd.aas.org
 *
 *  Usage:
 *    drms_rebin [-nv] in= ...
 *
 *  Arguments: (type		default         description)
 *    in	DataSet		-		Input dataset
 *    out	DataSeries	in		Output dataseries
 *    seg	String		"Not Specified	Output data segment if required
 *    copy	String		"+"		List of keys to propagate from
 *			input to output records, if possible
 *
 *  Flags
 *	-n	do not save results (diagnostic only)
 *	-v	run verbose
 *
 *  Notes:
 *    The output segment structure needs to be checked for consistency with
 *	the input segment and bin and extract paremeters as follows for the
 *	different protocol cases:
 *	in variable/constant	out variable/constant	once at beginning
 *	in variable/constant	out vardim		never
 *	in vardim		out variable/constant	always
 *	in vardim		out vardim		never
 *
 *  Bugs:
 *    Recalculation of target arrays for vardim input has not been tested
 *    Does not check quality of input records, causing errors when it tries to
 *	open files for missing data
 *    Only one input segment and one output segment can be processed; there
 *	is currently no way of dealing with the situation when there are
 *	multiple possible candidate segments and the input segments must
 *	be specified in the input dataset specfiers
 *    Has not been tested with arrays of rank > 3
 *    There may be some inconsistency about use of long long's for indexing in
 *	possibly large arrays
 *
 *  Revision history is at end of file.
 */
#include <jsoc_main.h>
#include "keystuff.c"

char *module_name = "drms_rebin";
char *version_id = "0.8";
char *module_desc = "rectangular region binning";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in",   "mdi.fd_V_lev18[:#$]", "input data set"},
//  {ARG_STRING, "in",    "mdi.fd_V_lev18[2008.01.01_00:40]", "input data set"},
  {ARG_STRING, "out",   "su_rsb.tstrebin", "output data series"},
  {ARG_STRING, "seg",   "Not Specified", "output data segment"},
  {ARG_INTS,   "bin",   "{1}", "array of per-axis bin widths"},
  {ARG_INTS,   "start", "{0}", "array of input axis starting pixels"},
  {ARG_INTS,   "stop",  "{-1}", "array of input axis ending pixels"},
  {ARG_STRING, "copy",  "+",
      "comma separated list of keys to propagate forward"},
  {ARG_FLAG,   "n", "", "do not save output (for testing)"},
  {ARG_FLAG,   "v", "", "run in verbose mode"},
  {ARG_END}
};
		 /*  list of keywords to propagate by default (if possible)
						       from input to output  */
char *propagate[] = {"T_REC", "T_OBS"};

			/*  the following belong in external utility files  */
/*
 *  Parse a token separated list of character strings into an array of
 *    character strings and return the number of such strings, plus the
 *    array itself in the argument list
 *
 *  Bugs:
 *    There is no way of including the token in the parsed strings
 */
int construct_stringlist (const char *request, char token, char ***stringlist) {
  int keyct = 1;
  int m, n;
  char *req, *req0 = strdup (request);
  char c;
					/*  count the number of separators  */
  req = req0;
  while (c = *req) {
    if (c == token) {
      *req = '\0';
      keyct++;
    }
    req++;
  }
  *stringlist = (char **)malloc (keyct * sizeof (char **));
  req = req0;
  for (n = 0; n < keyct; n++) {
    (*stringlist)[n] = strdup (req);
    req += strlen (req) + 1;
  }
  for (n = 0; n < keyct; n++) {
    char *subs = (*stringlist)[n];
					     /*  remove leading white space  */
    while (isspace (c = *subs)) subs++;
    (*stringlist)[n] = subs;
					    /*  remove trailing white space  */
    if (strlen (subs)) {
      subs += strlen (subs) - 1;
      while (isspace (c = *subs)) {
	*subs = '\0';
	subs--;
      }
    }
  }
					 /*  remove empty strings from list  */
  for (n = 0; n < keyct; n++) {
    if (!strlen ((*stringlist)[n])) {
      for (m = n; m < keyct - 1; m++)
        (*stringlist)[m] = (*stringlist)[m + 1];
      keyct--;
    }
  }
  free (req0);
  return keyct;
}

							/*  main module body  */
int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  int propct;
  char **copykeylist;
  char module_ident[128];
  DRMS_RecordSet_t *drs = NULL;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *orig_array, *bind_array;
  HIterator_t *lastseg;
  double *vdat, *vbin;
  double crpix, cdelt;
  long long *nssub, *ntsub;
  long long nn, ntdat, ntbin;
  int **loc;
  int *np;
  int *in_axes, *out_axes;
  int *bin, *strt, *stop, *strt0, *stop0;
  int axis, npmax;
  int m, n, key_n, rec_n, rec_ct, seg_n, iseg_ct, oseg_ct, valid_ct;
  int isegnum, osegnum, rank, axis_check;
  int kstat, status = 0;
  int keyct = sizeof (propagate) / sizeof (char *);
  char source[DRMS_MAXQUERYLEN];
  char key[256];

  double missing_val = 0.0 / 0.0;

  char *inds = strdup (params_get_str (params, "in"));
  char *outser = strdup (params_get_str (params, "out"));
  char *out_segname = strdup (params_get_str (params, "seg"));
  int binvals = params_get_int (params, "bin_nvals");
  int strtvals = params_get_int (params, "start_nvals");
  int stopvals = params_get_int (params, "stop_nvals");
  char *propagate_req = strdup (params_get_str (params, "copy"));
  int no_save = params_isflagset (params, "n");
  int verbose = params_isflagset (params, "v");
  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;

  snprintf (module_ident, 128, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);

  drs = drms_open_records (drms_env, inds, &status);
  if (!drs) {
    fprintf (stderr, "Error: unable to open record set %s\n", inds);
    return 1;
  }
  rec_ct = drs->n;
  if (rec_ct < 1) {
    fprintf (stderr, "No records in selected set %s\n", inds);
    return 1;
  }
  if (!strcmp (outser, "in")) {
    outser = strdup (inds);
    for (n = 0; n < strlen (outser); n++) if (outser[n] == '[') {
      outser[n] = '\0';
      break;
    }
  }
 				    /*  create sample output series record  */
  orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
  if (!orec) {
    fprintf (stderr, "Error: unable to create records in series %s\n",
	outser);
    fprintf (stderr, "       drms_create_record() returned status %d\n",
	status); 
    return 1;
  }
 					 /*  check input data series struct  */
  irec = drs->records[0];
  iseg_ct = drms_record_numsegments (irec);
  valid_ct = 0;
  lastseg = NULL;
		    /*  must use iterator in case segment name is included
					   as part of record set specifier  */
		/*  must add another argument (0?) to call in next release  */
  while (iseg = drms_record_nextseg (irec, &lastseg, 1)) {
/*
printf ("\t%-10s\t%d\n", iseg->info->name, iseg->info->segnum);
*/
    if (iseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	iseg->info->protocol == DRMS_GENERIC) continue;
    if (!valid_ct) isegnum = iseg->info->segnum;
    valid_ct++;
  }
  if (valid_ct > 1) {
    fprintf (stderr, "Error: input data set %s\n", inds);
    fprintf (stderr,
	"       contains multiple segments of appropriate protocol;\n");
    fprintf (stderr, "       in segment must be specified\n");
    drms_free_record (orec);
    drms_close_records (drs, DRMS_FREE_RECORD);
    return 1;
  }
  if (valid_ct < 1) {
    fprintf (stderr, "Error: input data set %s\n", inds);
    fprintf (stderr, "       contains no segments of appropriate protocol\n");
    drms_free_record (orec);
    drms_close_records (drs, DRMS_FREE_RECORD);
    return 1;
  }
  iseg = drms_segment_lookupnum (irec, isegnum);
  rank = iseg->info->naxis;
 					/*  check output data series struct  */
  oseg_ct = drms_record_numsegments (orec);
  if (strcmp (out_segname, "Not Specified")) {
	  /*  make sure the named output segment is writeable from an array  */
    oseg = drms_segment_lookup (orec, out_segname);
    if (!oseg) {
      fprintf (stderr, "Error: no such data segment %s\n", out_segname);
      drms_free_record (orec);
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
    if (oseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	oseg->info->protocol == DRMS_GENERIC) {
      fprintf (stderr, "Error: output data segment %s\n", out_segname);
      fprintf (stderr, "       is not of appropriate protocol\n");
      drms_free_record (orec);
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
    valid_ct = 1;
    osegnum = oseg->info->segnum;
  } else {
		 /*  find the (only) output segment writeable from an array  */
    for (seg_n = 0, valid_ct = 0; seg_n < oseg_ct; seg_n++) {
      oseg = drms_segment_lookupnum (orec, seg_n);
      if (oseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	  oseg->info->protocol == DRMS_GENERIC) continue;
      if (!valid_ct) osegnum = seg_n;
      valid_ct++;
    }
    if (valid_ct > 1) {
      fprintf (stderr, "Error: output data series %s\n", outser);
      fprintf (stderr,
	  "       contains multiple segments of appropriate protocol;\n");
      fprintf (stderr, "       seg must be specified\n");
      drms_free_record (orec);
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
  }
  oseg = drms_segment_lookupnum (orec, osegnum);
  if (oseg->info->naxis != rank) {
    fprintf (stderr,
	"Error: ranks of input and output data segments do not match\n");
    drms_free_record (orec);
    drms_close_records (drs, DRMS_FREE_RECORD);
    return 1;
  }
  axis_check = (oseg->info->scope != DRMS_VARDIM);
		   /*  fill the per-axis arrays of bin, start & stop values,
				extending with the last element as necessary  */
  bin = (int *)malloc (rank * sizeof (int));
  strt = (int *)malloc (rank * sizeof (int));
  stop = (int *)malloc (rank * sizeof (int));
  strt0 = (int *)malloc (rank * sizeof (int));
  stop0 = (int *)malloc (rank * sizeof (int));
  if (binvals > rank) binvals = rank;
  if (strtvals > rank) strtvals = rank;
  if (stopvals > rank) stopvals = rank;
  for (n = 0; n < binvals; n++) {
    sprintf (key, "bin_%d_value", n);
    bin[n] = params_get_int (params, key);
  }
  while (n < rank) {
    bin[n] = bin[n-1];
    n++;
  }
  for (n = 0; n < strtvals; n++) {
    sprintf (key, "start_%d_value", n);
    strt[n] = params_get_int (params, key);
  }
  while (n < rank) {
    strt[n] = strt[n-1];
    n++;
  }
  for (n = 0; n < stopvals; n++) {
    sprintf (key, "stop_%d_value", n);
    stop[n] = params_get_int (params, key);
  }
  while (n < rank) {
    stop[n] = stop[n-1];
    n++;
  }
			  /*  save argument values for case of vardim input  */
  for (n = 0; n < rank; n++) {
    strt0[n] = strt[n];
    stop0[n] = stop[n];
  }
    /*  initialize arrays and calculate target pixel lists for first record  */
		/*  determine source locations for each binned target pixel  */
  in_axes = (int *)malloc (rank * sizeof (int));
  out_axes = (int *)malloc (rank * sizeof (int));
  nssub = (long long *)malloc (rank * sizeof (long long));
  ntsub = (long long *)malloc (rank * sizeof (long long));
  ntdat = ntbin = npmax = 1;
  nssub[0] = ntsub[0] = 1;
  for (n = 0; n < rank; n++) {
    in_axes[n] = iseg->axis[n];
					   /*  adjust start and stop values  */
    while (strt[n] < 0) strt[n] += in_axes[n];
    while (stop[n] < 0) stop[n] += in_axes[n];
    while (strt[n] >= in_axes[n]) strt[n] -= in_axes[n];
    while (stop[n] >= in_axes[n]) stop[n] -= in_axes[n];
    if (stop[n] < strt[n]) {
      int save = strt[n];
      strt[n] = stop[n];
      stop[n] = save;
    }

    axis = stop[n] - strt[n] + 1;
    out_axes[n] = axis / bin[n];
    if (axis % bin[n]) out_axes[n]++;
    ntdat *= in_axes[n];
    ntbin *= out_axes[n];
    npmax *= bin[n];
    if (n) {
      nssub[n] = nssub[n-1] * in_axes[n-1];
      ntsub[n] = ntsub[n-1] * out_axes[n-1];
    }
  }
		       /*  this check also belongs below if input is vardim  */
  if (axis_check) {
    /*  the output and the input both have fixed sizes; they'd better match  */
    for (n = 0; n < rank; n++) {
      if (out_axes[n] != oseg->axis[n]) {
	fprintf (stderr,
	    "Error: array mismatch between output segment definition\n");
	fprintf (stderr,
	    "       and desired binning on fixed input segment definition\n");
	fprintf (stderr, "axis[%d] = %d, should be %d\n", n, oseg->axis[n],
	    out_axes[n]);
	drms_free_record (orec);
	drms_close_records (drs, DRMS_FREE_RECORD);
	return 1;
      }
    }
    if (iseg->info->scope != DRMS_VARDIM) axis_check = 0;
  }
  np = (int *)calloc (ntbin, sizeof (int));
  loc = (int **)malloc (ntbin * sizeof (int *));
  for (nn = 0; nn < ntbin; nn++) loc[nn] = (int *)malloc (npmax * sizeof (int));
		/*  determine source locations for each binned target pixel  */
  for (nn = 0; nn < ntdat; nn++) {
    int trgpix;
    int intarget = 1;
    int srcpix = (nn % in_axes[0]);
    if (srcpix < strt[0] || srcpix >= stop[0]) intarget = 0;
    trgpix = (srcpix - strt[0]) / bin[0];
    for (m = 1; m < rank; m++) {
      srcpix = (nn / nssub[m]) % in_axes[m];
      if (srcpix < strt[m] || srcpix >= stop[m]) {
	intarget = 0;
	break;
      }
      trgpix += ((srcpix - strt[m]) / bin[m]) * ntsub[m];
    }
    if (intarget) {
      loc[trgpix][np[trgpix]] = nn;
      np[trgpix]++;
    }
  }

  vbin = (double *)malloc (ntbin * sizeof (double));
  int *bind_axes = (int *)malloc (rank * sizeof (int));
  memcpy (bind_axes, out_axes, rank * sizeof (int));
  bind_array = drms_array_create (DRMS_TYPE_DOUBLE, rank, out_axes,
      (void *)vbin, &status);
  bind_array->bscale = oseg->bscale;
  bind_array->bscale = oseg->bscale;
  drms_free_record (orec);

  propct = construct_stringlist (propagate_req, ',', &copykeylist);
  if (verbose) {
    printf ("propagating %d key(s):\n", propct);
    for (n = 0; n < propct; n++) printf ("  %s\n", copykeylist[n]);
    printf ("\nprocessing %d record(s) in series %s:\n", rec_ct, inds);
  }

  if (axis_check) {
    fprintf (stderr,
	"Error: either input or output segment has protocol vardim\n");
    fprintf (stderr, "       not currently supported\n");
    return 1;
  }

  for (rec_n = 0; rec_n < rec_ct; rec_n++) {
    if (verbose) printf ("record %d:\n", rec_n);
    irec = drs->records[rec_n];
    drms_sprint_rec_query (source, irec);
    if (verbose) printf ("  processing record %d: %s\n", rec_n, source);
    iseg = drms_segment_lookupnum (irec, isegnum);
    if (!iseg) {
      fprintf (stderr, "Warning: could not find segment # %d\n", isegnum);
      fprintf (stderr, "       in %s; skipped\n", source);
      continue;
    }
			/*  read data from segment into array as doubles
				       (will block until segment is staged)  */
    orig_array = drms_segment_read (iseg, DRMS_TYPE_DOUBLE, &status);
    if (!orig_array) {
      fprintf (stderr, "Warning: could not read segment # %d\n", isegnum);
      fprintf (stderr, "       in %s; skipped\n", source);
      continue;
    }
    vdat = (double *)orig_array->data;
    orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    oseg = drms_segment_lookupnum (orec, osegnum);
    if (iseg->info->scope == DRMS_VARDIM) {
	       /*  input series has variable dimensions; have they changed?  */
      int recalc = 0, recok = 1;
      for (n = 0; n < rank; n++) {
	if (in_axes[n] != iseg->axis[n]) recalc = 1;
	in_axes[n] = iseg->axis[n];
      }
      if (recalc) {
				  /*  yes, need to recalculate target array  */
	for (n = 0; n < rank; n++) {
	  strt[n] = strt0[n];
	  stop[n] = stop0[n];
	}
	ntdat = ntbin = npmax = 1;
	ntsub[0] = nssub[0] = 1;
	for (n = 0; n < rank; n++) {
					   /*  adjust start and stop values  */
	  while (strt[n] < 0) strt[n] += in_axes[n];
	  while (stop[n] < 0) stop[n] += in_axes[n];
	  while (strt[n] >= in_axes[n]) strt[n] -= in_axes[n];
	  while (stop[n] >= in_axes[n]) stop[n] -= in_axes[n];
	  if (stop[n] < strt[n]) {
	    int save = strt[n];
	    strt[n] = stop[n];
	    stop[n] = save;
	  }
	  axis = stop[n] - strt[n] + 1;
	  out_axes[n] = axis / bin[n];
	  if (axis % bin[n]) out_axes[n]++;
	  ntdat *= in_axes[n];
 	  ntbin *= out_axes[n];
	  npmax *= bin[n];
	  if (n) {
	    nssub[n] = nssub[n-1] * in_axes[n-1];
	    ntsub[n] = ntsub[n-1] * out_axes[n-1];
	  }
	}
	if (axis_check) {
				/*  the output has fixed sizes, so check it  */
	  for (n = 0; n < rank; n++) {
	    if (out_axes[n] != oseg->axis[n]) {
	      fprintf (stderr,
		  "Warning: array mismatch between output segment definition\n");
	      fprintf (stderr,
		  "       and desired binning on input segment in\n");
	      fprintf (stderr, "       %s; skipped\n", source);
	      drms_free_record (orec);
	      recok = 0;
	      break;
            }
          }
	  if (!recok) continue;
	}

	np = (int *)realloc (np, ntbin * sizeof (int));
	for (nn = 0; nn < ntbin; nn++) np[nn] = 0;
	loc = (int **)realloc (loc, ntbin * sizeof (int *));
	for (nn = 0; nn < ntbin; nn++) {
	  if (loc[nn]) free (loc[nn]);
	  loc[nn] = (int *)malloc (npmax * sizeof (int));
	}
	for (nn = 0; nn < ntdat; nn++) {
	  int trgpix;
	  int intarget = 1;
	  int srcpix = (nn % in_axes[0]);
	  if (srcpix < strt[0] || srcpix >= stop[0]) intarget = 0;
	  trgpix = (srcpix - strt[0]) / bin[0];
	  for (m = 1; m < rank; m++) {
	    srcpix = (nn / nssub[m]) % in_axes[m];
	    if (srcpix < strt[m] || srcpix >= stop[m]) {
	      intarget = 0;
	      break;
	    }
	    trgpix += ((srcpix - strt[m]) / bin[m]) * ntsub[m];
	  }
	  if (intarget) {
	    loc[trgpix][np[trgpix]] = nn;
	    np[trgpix]++;
	  }
	}
      }
    }
			      /*  for now, match scaling of output to input  */
/*
    bind_array->bscale = orig_array->bscale;
    bind_array->bzero = orig_array->bzero;
*/
    for (nn = 0; nn < ntbin; nn++) {
      vbin[nn] = 0;
      for (m = 0; m < np[nn]; m++) {
        vbin[nn] += vdat[loc[nn][m]];
      }
					 /*  adjust scaling of binned array  */
      if (np[nn]) vbin[nn] /= np[nn];
      else vbin[nn] = missing_val;
    }
    drms_free_array (orig_array);
    
    drms_segment_write (oseg, bind_array, 0);
    kstat = 0;
    for (key_n = 0; key_n < propct; key_n++) {
      if (strcmp (copykeylist[key_n], "+"))
        kstat += check_and_copy_key (orec, irec, copykeylist[key_n]);
      else kstat += propagate_keys (orec, irec, propagate, keyct);
    }
    kstat += check_and_set_key_time (orec, "Date", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_str (orec, "Source", source);
							/*  adjust WCS keys  */
    for (n = 0; n < rank; n++) {
      sprintf (key, "CTYPE%d", n + 1);
      kstat += check_and_copy_key (orec, irec, key);
      sprintf (key, "CRVAL%d", n + 1);
      kstat += check_and_copy_key (orec, irec, key);
      sprintf (key, "CRPIX%d", n + 1);
      crpix = drms_getkey_double (irec, key, &status);
      if (!status) {
        crpix -= strt[n];
        crpix /= bin[n];
	kstat += check_and_set_key_double (orec, key, crpix);
      }
      sprintf (key, "CDELT%d", n + 1);
      cdelt = drms_getkey_double (irec, key, &status);
      if (!status) {
        cdelt *= bin[n];
	kstat += check_and_set_key_double (orec, key, cdelt);
      }
    }
    drms_close_record (orec, dispose);
  }
  drms_close_records (drs, DRMS_FREE_RECORD);

  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  10.01.22	File created by R Bogart
 *  v 0.7	Binning, and setting of critical (WCS) keywords, but no
 *	image extraction
 *  v 0.7 frozen 2010.01.26
 *  v 0.8	Modified call to drms_record_nextseg for DRMS 2.0
 *		Changed indexing to quad long where appropriate
 *		Added recalculation of target arrays for vardim input
 */
