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
 *    Only one input segment and one output segment can be processed; there
 *	is currently no way of dealing with the situation when there are
 *	multiple possible candidate segments and the input segments must
 *	be specified in the input dataset specfiers
 *    The option for removal of full orbital velocity is only a stub, and
 *	unimplemented; the option for removal of radial velocity is superfluous,
 *	as it is the same as specifying offset= -OBS_VR
 *    Has not been tested with arrays of rank > 3
 *    There may be some inconsistency about use of long long's for indexing in
 *	possibly large arrays
 *
 *  Revision history is at end of file.
 */
#include <jsoc_main.h>
#include "keystuff.c"

char *module_name = "drms_rebin";
char *version_id = "0.9";
char *module_desc = "rectangular region binning";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in",    "", "input data set"},
  {ARG_STRING, "out",   "", "output data series"},
  {ARG_STRING, "seg",   "Not Specified", "output data segment"},
  {ARG_INTS,   "bin",   "{1}", "array of per-axis bin widths"},
  {ARG_INTS,   "start", "{0}", "array of input axis starting pixels"},
  {ARG_INTS,   "stop",  "{-1}", "array of input axis ending pixels"},
  {ARG_STRING, "copy",  "+",
      "comma separated list of keys to propagate forward"},
  {ARG_STRING, "scale",  "Not Specified", "scaling value (number or keyword)"},
  {ARG_STRING, "offset",  "Not Specified", "offset value (number or keyword)"},
  {ARG_INT,    "qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING, "qual_key", "Quality", "keyname for 32-bit input image quality field"}, 
  {ARG_FLAG,   "n", "", "do not save output (for testing)"},
  {ARG_FLAG,   "o", "",
      "remove orbital velocity from signal (radial component only)"},
  {ARG_FLAG,   "O", "", "remove orbital velocity from signal"},
  {ARG_FLAG,   "v", "", "run in verbose mode"},
  {ARG_END}
};
		 /*  list of keywords to propagate by default (if possible)
						       from input to output  */
char *propagate[] = {"T_REC", "T_OBS", "DATE__OBS", "TELESCOP", "INSTRUME"};

			/*  the following belong in external utility files  */
			   /*  parse a string of the form [*|/][+|-]keyname  */
char *parse_as_arith_key (const char *str, int *mult, int *add) {
  char *parsed = strdup (str);
  *mult = *add = 1;
  if (parsed[0] == '*') parsed++;
  else if (parsed[0] == '/') {
    *mult = -1;
    parsed++;
  }
  if (parsed[0] == '+') parsed++;
  else if (parsed[0] == '-') {
    *add = -1;
    parsed++;
  }
  return parsed;
}

int set_stats_keys (DRMS_Record_t *rec, DRMS_Array_t *v, long long ntot) {
  double *vavg = (double *)v->data;

  double vmin, vmax, norm, norm2, scale, sig, var;
  double sumv1, sum2, sum3, sum4;
  long long n;
  int sumv0;
  int valmin = ntot, valmax = 0;
  int kstat = 0;

  vmin = HUGE;
  vmax = -HUGE;
  sumv1 = 0.0;
  sumv0 = 0;
  for (n = 0; n < ntot; n++) {
    if (vavg[n] < vmin) vmin = vavg[n];
    if (vavg[n] > vmax) vmax = vavg[n];
    if (isfinite (vavg[n])) {
      sumv0++;
      sumv1 += vavg[n];
    }
  }
  
  kstat += check_and_set_key_double (rec, "DataMin", vmin);
  kstat += check_and_set_key_double (rec, "DataMax", vmax);

  if (sumv0) {
    scale = 1.0 / sumv0;
    sumv1 *= scale;
    sum2 = sum3 = sum4 = 0.0;
    kstat += check_and_set_key_double (rec, "DataMean", sumv1);
    for (n = 0; n < ntot; n++) {
      if (isfinite (vavg[n])) {
        norm = vavg[n] - sumv1;
	norm2 = norm * norm;
	sum2 += norm2;
	sum3 += norm * norm2;
	sum4 += norm2 * norm2;
      }
    }
    kstat += check_and_set_key_double (rec, "DataRMS", sqrt (sum2 / sumv0));
    if (sumv0 > 1) {
      var = sum2 / (sumv0 - 1);
      sig = sqrt (var);
      kstat += check_and_set_key_double (rec, "DataSkew",
	  scale * sum3 / (var * sig));
      kstat += check_and_set_key_double (rec, "DataKurt",
	  scale * sum4 / (var * var) - 3.0);
    }
  }

  return kstat;
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
  double scale, bias, vr, vw, vn;
  long long *nssub, *ntsub;
  long long nn, ntdat, ntbin;
  unsigned int qual_inp, qual_out;
  int **loc;
  int *np;
  int *in_axes, *out_axes;
  int *bin, *strt, *stop, *strt0, *stop0;
  int axis, npmax;
  int m, n, key_n, rec_n, rec_ct, seg_n, iseg_ct, oseg_ct, valid_ct;
  int isegnum, osegnum, rank, axis_check;
  int bias_mult, bias_sign, scale_mult, scale_sign;
  int kstat, status = 0;
  int keyct = sizeof (propagate) / sizeof (char *);
  char *key_scale, *key_bias;
  char source[DRMS_MAXQUERYLEN];
  char key[256], valuestr[256];

  double missing_val = 0.0 / 0.0;

  char *inds = strdup (params_get_str (params, "in"));
  char *outser = strdup (params_get_str (params, "out"));
  char *out_segname = strdup (params_get_str (params, "seg"));
  int binvals = params_get_int (params, "bin_nvals");
  int strtvals = params_get_int (params, "start_nvals");
  int stopvals = params_get_int (params, "stop_nvals");
  char *propagate_req = strdup (params_get_str (params, "copy"));
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  int no_save = params_isflagset (params, "n");
  int add_orbital_vr = params_isflagset (params, "o");
  int add_orbital_full = params_isflagset (params, "O");
  int scale_values = strcmp (params_get_str (params, "scale"), "Not Specified");
  int bias_values = strcmp (params_get_str (params, "offset"), "Not Specified");
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
 				       /*  check input data series keywords  */
  if (scale_values) {
    char *endptr;
    scale = strtod (params_get_str (params, "scale"), &endptr);
    if (strlen (endptr)) {
      key_scale = parse_as_arith_key (params_get_str (params, "scale"),
	  &scale_mult, &scale_sign);
      if (!drms_keyword_lookup (irec, key_scale, 1)) {
        fprintf (stderr, "Warning: required keyword %s not found\n", key_scale);
        fprintf (stderr, "         no scaling applied\n");
	scale_values = 0;
      }
    } else scale_mult = scale_sign = 0;
  }
  if (bias_values) {
    char *endptr;
    bias = strtod (params_get_str (params, "offset"), &endptr);
    if (strlen (endptr)) {
      key_bias = parse_as_arith_key (params_get_str (params, "offset"),
	  &bias_mult, &bias_sign);
      if (!drms_keyword_lookup (irec, key_bias, 1)) {
        fprintf (stderr, "Warning: required keyword %s not found\n", key_bias);
        fprintf (stderr, "         no offset applied\n");
	bias_values = 0;
      }
    } else bias_mult = bias_sign = 0;
  }

  if (add_orbital_vr || add_orbital_full) {
    if (!drms_keyword_lookup (irec, "OBS_VR", 1)) {
      fprintf (stderr, "Warning: required keyword %s not found\n", "OBS_VR");
      fprintf (stderr, "         no correction applied for observer velocity\n");
      add_orbital_vr = add_orbital_full = 0;
    }
    if (add_orbital_full) {
      add_orbital_vr = 0;
      if (!drms_keyword_lookup (irec, "OBS_VW", 1)) {
        fprintf (stderr, "Warning: required keyword %s not found\n", "OBS_VW");
        fprintf (stderr, "         no correction applied for observer velocity\n");
        add_orbital_full = 0;
      }
      if (!drms_keyword_lookup (irec, "OBS_VW", 1)) {
        fprintf (stderr, "Warning: required keyword %s not found\n", "OBS_VN");
        fprintf (stderr, "         no correction applied for observer velocity\n");
        add_orbital_full = 0;
      }
      if (rank != 2) {
	fprintf (stderr, "Warning: input data not image type\n");
	fprintf (stderr,
	    "         full correction of orbital velocity not supported\n");
        add_orbital_full = 0;
      }
      fprintf (stderr,
	  "Warning: full correction of orbital velocity not implemented\n");
      add_orbital_full = 0;
    }
  }
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
    if (srcpix < strt[0] || srcpix > stop[0]) intarget = 0;
    trgpix = (srcpix - strt[0]) / bin[0];
    for (m = 1; m < rank; m++) {
      srcpix = (nn / nssub[m]) % in_axes[m];
      if (srcpix < strt[m] || srcpix > stop[m]) {
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
	"Warning: either input or output segment has protocol vardim\n");
    fprintf (stderr,
	"         not currently supported; results may be garbage!\n");
/*
    return 1;
*/
  }

  for (rec_n = 0; rec_n < rec_ct; rec_n++) {
    if (verbose) printf ("record %d:\n", rec_n);
    irec = drs->records[rec_n];
    drms_sprint_rec_query (source, irec);
    orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    if (verbose) printf ("  processing record %d: %s\n", rec_n, source);
    kstat = 0;
		/*  copy in values for new record's prime keys if possible
			   (they may be recopied or even overwritten later)  */
    kstat += copy_prime_keys (orec, irec);
    for (key_n = 0; key_n < propct; key_n++) {
      if (strcmp (copykeylist[key_n], "+"))
        kstat += check_and_copy_key (orec, irec, copykeylist[key_n]);
      else kstat += propagate_keys (orec, irec, propagate, keyct);
    }
    kstat += check_and_set_key_time (orec, "Date", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_str (orec, "Input", inds);
    kstat += check_and_set_key_str (orec, "Source", source);
    kstat += check_and_set_key_str (orec, "BLD_VERS", jsoc_version);
			/*  check quality flag for presence of data segment  */
    qual_inp = drms_getkey_int (irec, qual_key, &status);
    if ((qual_inp & qmask) && !status) {
      kstat += check_and_set_key_int (orec, "Quality", 0x80000000);
      drms_close_record (orec, dispose);
      continue;
    }
    qual_out = 0x00000000;
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
    if (bias_values) {
      if (bias_mult) {
        bias = drms_getkey_double (irec, key_bias, &status);
	if (bias_mult < 0) bias = 1.0 / bias;
	bias *= bias_sign;
      }
      for (nn = 0; nn < ntdat; nn++) if (isfinite (bias)) vdat[nn] += bias;
    }
    if (scale_values) {
      if (scale_mult) {
        scale = drms_getkey_double (irec, key_scale, &status);
	if (scale_mult < 0) scale = 1.0 / scale;
	scale *= scale_sign;
      }
      for (nn = 0; nn < ntdat; nn++) if (isfinite (scale)) vdat[nn] *= scale;
    }
    if (add_orbital_vr) {
      vr = drms_getkey_double (irec, "OBS_VR", &status);
      for (nn = 0; nn < ntdat; nn++) if (isfinite (vr)) vdat[nn] -= vr;
    }
    if (add_orbital_full) {
      vr = drms_getkey_double (irec, "OBS_VR", &status);
      vw = drms_getkey_double (irec, "OBS_VW", &status);
      vn = drms_getkey_double (irec, "OBS_VN", &status);
    }
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
    kstat += check_and_set_key_int (orec, "Quality", qual_out);
							 /*  parameter keys  */
    for (n = 0; n < rank; n++) {
      sprintf (key, "binwdth_%d", n);
      kstat += check_and_set_key_int (orec, key, bin[n]);
      sprintf (key, "startpx_%d", n);
      kstat += check_and_set_key_int (orec, key, strt[n]);
      sprintf (key, "stoppx_%d", n);
      kstat += check_and_set_key_int (orec, key, stop[n]);
    }
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
      sprintf (key, "CROTA%d", n + 1);
      kstat += check_and_copy_key (orec, irec, key);
    }
    if (bias_values) {
      if (bias_mult) {
	if (bias_mult < 0) {
	  if (bias_sign < 0) snprintf (valuestr, 256, "1/(-%s)", key_bias);
	  else snprintf (valuestr, 256, "1/%s", key_scale);
	} else {
	  if (bias_sign < 0) snprintf (valuestr, 256, "-%s", key_bias);
	  else snprintf (valuestr, 256, "%s", key_bias);
	}
      } else snprintf (valuestr, 256, "%g", bias);
      kstat += check_and_set_key_str (orec, "OffsetBy", valuestr);
      if (verbose) {
        printf ("  offset by ");
	if (bias_mult) printf ("%s = ", valuestr);
	printf ("%g\n", bias);
      }
    }
    if (scale_values) {
      if (scale_mult) {
	if (scale_mult < 0) {
	  if (scale_sign < 0) snprintf (valuestr, 256, "1/(-%s)", key_scale);
	  else snprintf (valuestr, 256, "1/%s", key_scale);
	} else {
	  if (scale_sign < 0) snprintf (valuestr, 256, "-%s", key_scale);
	  else snprintf (valuestr, 256, "%s", key_scale);
	}
      } else snprintf (valuestr, 256, "%g", scale);
      kstat += check_and_set_key_str (orec, "ScaledBy", valuestr);
      if (verbose) {
        printf ("  scaled by ");
	if (scale_mult) printf ("%s = ", valuestr);
	printf ("%g\n", scale);
      }
    }
							/*  statistics keys  */
    kstat += set_stats_keys (orec, bind_array, ntbin);
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
 *  v 0.8 frozen 2010.03.08
 *  v 0.9	Put in hack to accept vardim segments on input, but not
 *	actually deal with them
 *		Added CROTAn, DATE__OBS, TELESCOP and INSTRUME to default
 *	propagated keywords
 *		Added copying of output's prime keys, setting of Input key
 *		Fixed bug that was preventing last entry of slowest-index
 *	from being properly filled
 *		Added options for removal of orbital velocity (2nd-order only
 *	a stub) and scaling and/or offset by fixed numerical value or keyword
 *	value, and setting of appropriate keywords; added (unconditional)
 *	setting of statistics keywords; added checking and setting of quality
 *	flags
 *		Moved construct_stringlist() function to keystuff
 *  v 0.9 frozen 2010.08.20
 */
