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
 *  Arguments: (type	default         description)
 *    in	DataSet	-	input dataset
 *    out	DataSeries in	output dataseries
 *    seg	str	Not Specoutput data segment if required
 *    bin	int*	{1}	array of per-axis bin widths
 *    start	int*	{0}	array of per-axis crop start pixels
 *    stop	int*	{-1}	array of per-axis crop end pixels
 *    wmin	float*	{NaN}
 *    copy	str	"+"	list of keys to propagate from input to output
 *	records, if possible
 *    scale	str	Not Spec
 *    offset	str	Not Spec
 *    qmask	Int	0x80000000
 *    qual_key	str	Quality
 *    idkey	str	Not Spec	prime key to set to pval
 *    idval	str	Not Spec	value to set for pkey
 *
 *  Flags
 *	-c	collapse output rank if possible
 *	-m	merge output into higher rank array
 *	-n	do not save results (diagnostic only)
 *	-o	remove radial component of orbital velocity from signal
 *	-O	remove radial component of vector orbital velocity from signal
 *		(unimplemented)
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
 *    Status of keyword setting not checked
 *    The merge option is a quick hack and does not work in many
 *	cases: fails with seg fault if rebinning 6 hours or more of GONG
 *	images for example
 *    The merge and collapse options cannot both be invoked; the merge option
 *	takes precedence
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
 *	large arrays
 *
 *  Revision history is at end of file.
 */

#define MODULE_VERSION_NUMBER	("1.1")
#define KEYSTUFF_VERSION "keystuff_v10.c"

#include <jsoc_main.h>
#include <fftw3.h>
#include KEYSTUFF_VERSION

char *module_name = "drms_rebin";
char *version_id = "1.1";
char *module_desc = "rectangular region binning";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in",    "", "input data set"},
  {ARG_STRING, "out",   "", "output data series"},
  {ARG_STRING, "seg",   "Not Specified", "output data segment"},
  {ARG_INTS,   "bin",   "{1}", "array of per-axis bin widths"},
  {ARG_INTS,   "start", "{0}", "array of input axis starting pixels"},
  {ARG_INTS,   "stop",  "{-1}", "array of input axis ending pixels"},
  {ARG_FLOATS, "wmin",  "{NaN}",
      "minimum wavelengths per axis for low-pass filtering"},
  {ARG_STRING, "copy",  "+",
      "comma separated list of keys to propagate forward"},
  {ARG_STRING, "scale",  "Not Specified", "scaling value (number or keyword)"},
  {ARG_STRING, "offset",  "Not Specified", "offset value (number or keyword)"},
  {ARG_INT,    "qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING, "qual_key", "Quality",
	"keyname for 32-bit input image quality field"},
  {ARG_STRING,	"idkey", "Not Specified", "name of (prime) key to set to idval"},
  {ARG_STRING,	"idval", "Not Specified",
      "value to set for idkey (:XXX for value of keyword XXX"},
  {ARG_FLAG,   "c", "", "collapse rank of output for unit axes"},
  {ARG_FLAG,   "m", "",
      "increase rank of output by merging along principal axis"},
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

static int data_filter (DRMS_Array_t *data) {
  double *wdata = (double *)data->data;
  double norm;
  float *xdata = (float *)data->data;
  long long n, ntot;
  int dp_calc;
  int *axes = data->axis;
  int rank = data->naxis;
  static char *valid = NULL;
				       /* Declarations for FFTW library calls */
  fftw_plan fplan, iplan;
  fftwf_plan rfplan, riplan;
  static fftw_complex *wform = NULL;
  static fftwf_complex *xform = NULL;

  if (data->type == DRMS_TYPE_DOUBLE) dp_calc = 1;
  else if (data->type == DRMS_TYPE_FLOAT) dp_calc = 0;
  else {
    fprintf (stderr, "Error: filtering of data of type %s not supported\n",
	drms_type_names[data->type]);
    return 1;
  }

  ntot = 1;
  for (n = 0; n < rank; n++) ntot *= axes[n];
  norm = 1.0 / ntot;
  valid = (char *)realloc (valid, ntot * sizeof (char));

  if (dp_calc) {
    wform = (fftw_complex *)realloc (wform, (ntot + 1) * sizeof (fftw_complex));
    fplan = fftw_plan_dft_r2c (rank, axes, wdata, wform, FFTW_ESTIMATE);
    iplan = fftw_plan_dft_c2r (rank, axes, wform, wdata, FFTW_ESTIMATE);
  } else {
    xform = (fftwf_complex *)realloc (xform, (ntot + 1) * sizeof (fftwf_complex));
    rfplan = fftwf_plan_dft_r2c (rank, axes, xdata, xform, FFTW_ESTIMATE);
    riplan = fftwf_plan_dft_c2r (rank, axes, xform, xdata, FFTW_ESTIMATE);
  }
						       /*  forward transform  */
  if (dp_calc) {
    for (n = 0; n < ntot; n++)
      if (isnan (wdata[n])) valid[n] = wdata[n] = 0;
      else valid[n] = 1;
    fftw_execute (fplan);
  } else {
    for (n = 0; n < ntot; n++)
      if (isnan (xdata[n])) valid[n] = xdata[n] = 0;
      else valid[n] = 1;
    fftwf_execute (rfplan);
  }
double rx2, rxfac, ry2, ryfac, rfilt;
int col, hcols, cols = axes[0];
int row, hrows, rows = axes[1];
hrows = (rows + 1) / 2;
hcols = cols / 2;		     /*  possible bug if cols odd  */
rxfac = 1.0 / cols;
ryfac = 1.0 / rows;
if (dp_calc) {
for (row = 0, n = 0; row < hrows; row++) {
  ry2 = row * ryfac;
  ry2 *= ry2;
  for (col = 0; col < cols; col++, n+=2) {
    rx2 = col * rxfac;
    rx2 *= rx2;
    rfilt = sqrt (rx2 + ry2);
if (rfilt > 0.1) wform[n] = wform[n+1] = 0.0;
  }
}
} else {
for (row = 0, n = 0; row < hrows; row++) {
  ry2 = row * ryfac;
  ry2 *= ry2;
  for (col = 0; col < cols; col++, n+=2) {
    rx2 = col * rxfac;
    rx2 *= rx2;
    rfilt = sqrt (rx2 + ry2);
if (rfilt > 0.1) xform[n] = xform[n+1] = 0.0;
  }
}
}
						       /*  inverse transform  */
  if (dp_calc) {
    fftw_execute (iplan);
    for (n = 0; n < ntot; n++) {
      if (!valid[n]) wdata[n] = NAN;
      else wdata[n] *= norm;
    }
  } else {
    fftwf_execute (riplan);
    for (n = 0; n < ntot; n++) {
      if (!valid[n]) xdata[n] = NAN;
      else xdata[n] *= norm;
    }
  }
  return 0;
}

static int check_input_series (char *inds, int *segnum) {
/*
 *  Check input data series structure for the existence of a unique segment of
 *    appropriate protocol for array segment reads, if the segment is not
 *    named as part of the dataset specification; otherwise check that the
 *    named segment has the appropriate protocol
 */
  DRMS_RecordSet_t *drs = NULL;
  DRMS_Record_t *irec;
  DRMS_Segment_t *iseg;
  HIterator_t *lastseg = NULL;
  int rec_ct;
  int status = 0, valid_ct = 0;

  *segnum = 0;
  drs = drms_open_records (drms_env, inds, &status);
  if (!drs) {
    fprintf (stderr, "Error: unable to open record set %s\n", inds);
    return 1;
  }
  rec_ct = drs->n;
  if (rec_ct < 1) {
    fprintf (stderr, "Error: no records in selected dataset %s\n", inds);
    return 1;
  }
  irec = drs->records[0];
		    /*  must use iterator in case segment name is included
					   as part of record set specifier  */
  while (iseg = drms_record_nextseg (irec, &lastseg, 1)) {
    if (iseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	iseg->info->protocol == DRMS_GENERIC) continue;
    if (!valid_ct) *segnum = iseg->info->segnum;
    valid_ct++;
  }
  if (valid_ct > 1) {
    fprintf (stderr, "Error: input data set %s\n", inds);
    fprintf (stderr,
	"       contains multiple segments of appropriate protocol;\n");
    fprintf (stderr, "       in segment must be specified\n");
    status = 1;
  }
  if (valid_ct < 1) {
    fprintf (stderr, "Error: input data set %s\n", inds);
    fprintf (stderr, "       contains no segments of appropriate protocol\n");
    status = 1;
  }
  drms_close_records (drs, DRMS_FREE_RECORD);
  return status;
}

static int check_output_series (char *series, char *segname, int *segnum,
    int *check) {
/*
 *  Check output data series structure for the existence of segments of
 *    appropriate protocol for array segment writes, if the segment is not
 *    explicitly named; otherwise check that the named segment has the
 *    appropriate protocol. If multiple acceptable output segments exist,
 *    set check argument.
 */
  DRMS_Record_t *orec;
  DRMS_Segment_t *oseg;
  int oseg_ct, seg_n;
  int status;

  orec = drms_create_record (drms_env, series, DRMS_TRANSIENT, &status);
  if (!orec) {
    fprintf (stderr, "Error: unable to create records in series %s\n", series);
    fprintf (stderr, "       drms_create_record() returned status %d\n",
	status); 
    return 1;
  }

  *check = 0;
  if (strcmp (segname, "Not Specified")) {
	  /*  make sure the named output segment is writeable from an array  */
    oseg = drms_segment_lookup (orec, segname);
    if (!oseg) {
      fprintf (stderr, "Error: no such data segment %s in series %s\n", segname,
	  series);
      drms_free_record (orec);
      return 1;
    }
    if (oseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	oseg->info->protocol == DRMS_GENERIC) {
      fprintf (stderr, "Error: output data segment %s\n", segname);
      fprintf (stderr, "       is not of appropriate protocol\n");
      drms_free_record (orec);
      return 1;
    }
    *segnum = oseg->info->segnum;
  } else {
    int valid_ct = 0;
    oseg_ct = drms_record_numsegments (orec);
		 /*  find the (only) output segment writeable from an array  */
    for (seg_n = 0; seg_n < oseg_ct; seg_n++) {
      oseg = drms_segment_lookupnum (orec, seg_n);
      if (oseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	  oseg->info->protocol == DRMS_GENERIC) continue;
      if (!valid_ct) *segnum = seg_n;
      valid_ct++;
    }
    if (valid_ct > 1) *check = 1;
    if (valid_ct < 1) {
      fprintf (stderr, "Error: output data series %s\n", series);
      fprintf (stderr, "       contains no segments of appropriate protocol\n");
      drms_free_record (orec);
      return 1;
    }
  }
  drms_free_record (orec);
  return 0;
}

static int check_output_segment (DRMS_Segment_t *seg, int rank, int *axes) {
  return 0;
}

static int find_output_segment (char *series, int rank, int *axes) {
/*
 *  Find the unique (or first) output data segment matching the requirements
 *    for appropriate protocol for array segment writes
 */
  DRMS_Record_t *orec;
  DRMS_Segment_t *oseg;
  int oseg_ct, seg_n, n;
  int status;
  int valid_ct = 0, segnum = -1;

  orec = drms_create_record (drms_env, series, DRMS_TRANSIENT, &status);
  oseg_ct = drms_record_numsegments (orec);
  for (seg_n = 0; seg_n < oseg_ct; seg_n++) {
    oseg = drms_segment_lookupnum (orec, seg_n);
    if (oseg->info->protocol == DRMS_PROTOCOL_INVALID ||
	oseg->info->protocol == DRMS_GENERIC) continue;
    if (oseg->info->naxis != rank) continue;
    if (oseg->info->scope != DRMS_VARDIM) {
      for (n = 0; n < rank; n++) {
		  /*  the output segment has fixed axis sizes, so check them  */
	if (axes[n] != oseg->axis[n]) continue;
      }
    }
    if (!valid_ct) segnum = seg_n;
    valid_ct++;
  }

  drms_free_record (orec);
  if (valid_ct > 1) {
    oseg = drms_segment_lookupnum (orec, segnum);
    fprintf (stderr, "Warning: output data series %s\n", series);
    fprintf (stderr,
	"       contains multiple segments of appropriate protocol and dimensions;\n");
    fprintf (stderr,
	"       using first matching segment: %s.\n", oseg->info->name);
    fprintf (stderr,
	"       To use another, seg must be specified\n");
  }
  drms_free_record (orec);
  return segnum;
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
  DRMS_Array_t *orig_array, *binned_array;
  double *vdat, *vbin;
  double crpix, cdelt;
  double scale, bias, vr, vw, vn;
  float *wmin;
  long long *nssub, *ntsub;
  long long nn, ntdat, ntbin;
  unsigned int qual_inp, qual_out;
  int **loc;
  int *np;
  int *in_axes, *out_axes, *oslice_strt, *oslice_stop;
  int *bin, *strt, *stop, *strt0, *stop0;
  int axis, npmax;
  int m, n, key_n, rec_n, rec_ct;
  int isegnum, osegnum, rank, crank, mrank, axis_check;
  int bias_mult, bias_sign, scale_mult, scale_sign, scaling_warned;
  int prefilter;
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
  int wminvals = params_get_int (params, "wmin_nvals");
  char *propagate_req = strdup (params_get_str (params, "copy"));
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *idkey = strdup (params_get_str (params, "idkey"));
  char *idval = strdup (params_get_str (params, "idval"));
  int collapse = params_isflagset (params, "c");
  int merge = params_isflagset (params, "m");
  int no_save = params_isflagset (params, "n");
  int add_orbital_vr = params_isflagset (params, "o");
  int add_orbital_full = params_isflagset (params, "O");
  int scale_values = strcmp (params_get_str (params, "scale"), "Not Specified");
  int bias_values = strcmp (params_get_str (params, "offset"), "Not Specified");
  int verbose = params_isflagset (params, "v");
  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;

  snprintf (module_ident, 128, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);

 				      /*  check input data series structure  */
  if (check_input_series (inds, &isegnum)) return 1;

  if (!strcmp (outser, "in")) {
    outser = strdup (inds);
    for (n = 0; n < strlen (outser); n++) if (outser[n] == '[') {
      outser[n] = '\0';
      break;
    }
  }
 				     /*  check output data series structure  */
  if (check_output_series (outser, out_segname, &osegnum, &axis_check)) {
    drms_close_records (drs, DRMS_FREE_RECORD);
    return 1;
  }
						    /*  open input data set  */
  drs = drms_open_records (drms_env, inds, &status);
  rec_ct = drs->n;
  irec = drs->records[0];
  iseg = drms_segment_lookupnum (irec, isegnum);
  crank = mrank = rank = iseg->info->naxis;
  if (merge) {
    mrank++;
    oslice_strt = (int *)malloc (mrank * sizeof (int));
    oslice_stop = (int *)malloc (mrank * sizeof (int));
		     /*  forbid both collapse and merge, too complex for now  */
    collapse = 0;
  }
	  /*  fill the per-axis arrays of bin, start, stop, and wmin values,
				extending with the last element as necessary  */
  bin = (int *)malloc (mrank * sizeof (int));
  strt = (int *)malloc (mrank * sizeof (int));
  stop = (int *)malloc (mrank * sizeof (int));
  strt0 = (int *)malloc (mrank * sizeof (int));
  stop0 = (int *)malloc (mrank * sizeof (int));
  wmin = (float *)malloc (mrank * sizeof (float));
  if (binvals > mrank) binvals = mrank;
  if (strtvals > mrank) strtvals = mrank;
  if (stopvals > mrank) stopvals = mrank;
  if (wminvals > mrank) wminvals = mrank;
  for (n = 0; n < binvals; n++) {
    sprintf (key, "bin_%d_value", n);
    bin[n] = params_get_int (params, key);
  }
  while (n < mrank) {
    bin[n] = bin[n-1];
    n++;
  }
  for (n = 0; n < strtvals; n++) {
    sprintf (key, "start_%d_value", n);
    strt[n] = params_get_int (params, key);
  }
  while (n < mrank) {
    strt[n] = strt[n-1];
    n++;
  }
  for (n = 0; n < stopvals; n++) {
    sprintf (key, "stop_%d_value", n);
    stop[n] = params_get_int (params, key);
  }
  while (n < mrank) {
    stop[n] = stop[n-1];
    n++;
  }
  for (n = 0; n < wminvals; n++) {
    sprintf (key, "wmin_%d_value", n);
    wmin[n] = params_get_float (params, key);
  }
  while (n < mrank) {
    wmin[n] = wmin[n-1];
    n++;
  }
			   /*  save argument values for case of vardim input  */
  for (n = 0; n < mrank; n++) {
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
  if (collapse) {
    int *new_axes = (int *)malloc (rank * sizeof (int));
    for (n = 0, m = 0; n < rank; n++) {
      if (out_axes[n] == 1) crank--;
      else new_axes[m++] = out_axes[n];
    }
					   /*  do not collapse to zero rank  */
    if (crank < 1) {
      crank = 1;
      new_axes[0] = 1;
    }
    if (crank == rank) collapse = 0;
    else for (n = 0; n < crank; n++) out_axes[n] = new_axes[n];
    free (new_axes);
  }
  if (merge) {
    strt[rank] = 0;
    stop[rank] = rec_ct - 1;
    out_axes[rank] = rec_ct;
  }
	/*  use first input record segment to find appropriate output segment
								   if needed  */
  if (axis_check) {
    osegnum = find_output_segment (outser, crank, out_axes);
    if (osegnum < 0) return 1;
  }
 				      /*  create sample output series record  */
  orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
  oseg = drms_segment_lookupnum (orec, osegnum);
  axis_check = (oseg->info->scope != DRMS_VARDIM);
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
  if (merge) {
    if (oseg->info->naxis != mrank) {
      fprintf (stderr,
	  "Error: rank of output data segment inconsistent with merged input\n");
      drms_free_record (orec);
      drms_close_records (drs, DRMS_FREE_RECORD);
      return 1;
    }
  } else if (!collapse && oseg->info->naxis != rank) {
    fprintf (stderr,
	"Error: ranks of input and output data segments do not match\n");
    drms_free_record (orec);
    drms_close_records (drs, DRMS_FREE_RECORD);
    return 1;
  }
  axis_check = (oseg->info->scope != DRMS_VARDIM);
	  /*  fill the per-axis arrays of bin, start, stop, and wmin values,
				extending with the last element as necessary  */
/*
  bin = (int *)malloc (rank * sizeof (int));
  strt = (int *)malloc (rank * sizeof (int));
  stop = (int *)malloc (rank * sizeof (int));
  strt0 = (int *)malloc (rank * sizeof (int));
  stop0 = (int *)malloc (rank * sizeof (int));
  wmin = (float *)malloc (rank * sizeof (float));
  if (binvals > rank) binvals = rank;
  if (strtvals > rank) strtvals = rank;
  if (stopvals > rank) stopvals = rank;
  if (wminvals > rank) wminvals = rank;
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
  for (n = 0; n < wminvals; n++) {
    sprintf (key, "wmin_%d_value", n);
    wmin[n] = params_get_float (params, key);
  }
  while (n < rank) {
    wmin[n] = wmin[n-1];
    n++;
  }
*/
			   /*  save argument values for case of vardim input  */
/*
  for (n = 0; n < rank; n++) {
    strt0[n] = strt[n];
    stop0[n] = stop[n];
  }
*/
						/*  prefilter to be applied?  */
  prefilter = isfinite (wmin[0]);
  for (n = 1; n < rank; n++) {
    if (prefilter && isnan (wmin[n])) {
      fprintf (stderr, "Warning: not all wmin values valid\n");
      fprintf (stderr, "         filtering turned off\n");
      prefilter = 0;
      break;
    }
    if (isfinite (wmin[n])) prefilter = 1;
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
	    /*  don't worry about missing axis if collapse option is invoked  */
        if (collapse && oseg->axis[n] == 0) continue;
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
  binned_array = drms_array_create (DRMS_TYPE_DOUBLE, rank, out_axes,
      (void *)vbin, &status);
/*
  binned_array->bscale = oseg->bscale;
  binned_array->bzero = oseg->bzero;
*/
  if (merge) {
    for (n = 0; n < rank; n++) {
      oslice_strt[n] = 0;
      oslice_stop[n] = binned_array->axis[n];
    }
  }
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

  if (merge) {
    orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    oseg = drms_segment_lookupnum (orec, osegnum);
    if (oseg->info->scope == DRMS_VARDIM) {
      for (n = 0; n < rank; n++) oseg->axis[n] = binned_array->axis[n];
      oseg->axis[rank] = rec_ct;
    }
    kstat = 0;
    kstat += check_and_set_key_time (orec, "Date", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_str (orec, "Input", inds);
    kstat += check_and_set_key_str (orec, "BLD_VERS", jsoc_version);
							 /*  parameter keys  */
    for (n = 0; n < rank; n++) {
      sprintf (key, "binwdth_%d", n);
      kstat += check_and_set_key_int (orec, key, bin[n]);
      sprintf (key, "startpx_%d", n);
      kstat += check_and_set_key_int (orec, key, strt[n]);
      sprintf (key, "stoppx_%d", n);
      kstat += check_and_set_key_int (orec, key, stop[n]);
    }
    irec = drs->records[0];
    drms_sprint_rec_query (source, irec);
    kstat += check_and_set_key_str (orec, "startrec", source);
    irec = drs->records[rec_ct-1];
    drms_sprint_rec_query (source, irec);
    kstat += check_and_set_key_str (orec, "stoprec", source);
  }
  scaling_warned = 0;
  for (rec_n = 0; rec_n < rec_ct; rec_n++) {
    irec = drs->records[rec_n];
    drms_sprint_rec_query (source, irec);
    if (verbose) printf ("processing record %d: %s\n", rec_n, source);
    if (merge) {
      oslice_stop[rank] = rec_n;
    } else {
      orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
			    /*  set basic key values (including optional id)  */
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
      if (strcmp (idkey, "Not Specified")) {
					    /*  set the special id key value  */
	DRMS_Type_t keytyp;
	drms_getkey (orec, idkey, &keytyp, &status);
	if (!status) {
	  TIME tval;
	  double dval;
          float fval;
	  int ival;
	  short sval;
	  char *strval;
	  switch (keytyp) {
	    case (DRMS_TYPE_SHORT):
	      sval = params_get_short (params, "idval");
	      kstat += check_and_set_key_short (orec, idkey, sval);
	      break;
	  case (DRMS_TYPE_INT):
	  ival = params_get_short (params, "idval");
	  kstat += check_and_set_key_int (orec, idkey, ival);
	  break;
	  case (DRMS_TYPE_FLOAT):
	  fval = params_get_float (params, "idval");
	  kstat += check_and_set_key_float (orec, idkey, fval);
	  break;
	  case (DRMS_TYPE_DOUBLE):
	  dval = params_get_double (params, "idval");
	  kstat += check_and_set_key_double (orec, idkey, dval);
	  break;
	  case (DRMS_TYPE_TIME):
	  tval = params_get_time (params, "idval");
	  kstat += check_and_set_key_time (orec, idkey, tval);
	  break;
	  case (DRMS_TYPE_STRING):
	  strval = strdup (params_get_str (params, "idval"));
	  kstat += check_and_set_key_str (orec, idkey, strval);
	  break;
	  default:
	  fprintf (stderr, "Warning unsupported type %s for keyword %s\n",
	      drms_type_names[keytyp], idkey);
	  }
	} else if (verbose) {
	  printf ("Warning: requested ID key %s is not in output data series\n",
	      idkey);
	}
      }
    }
			/*  check quality flag for presence of data segment  */
    qual_inp = drms_getkey_int (irec, qual_key, &status);
    if ((qual_inp & qmask) && !status) {
      if (!merge) kstat += check_and_set_key_int (orec, "Quality", 0x80000000);
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
    if (prefilter) {
					/*  apply optional prefilter to data  */
      status = data_filter (orig_array);
      fprintf (stderr, "Warning: prefiltering not implemented\n");
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
    binned_array->bscale = orig_array->bscale;
    binned_array->bzero = orig_array->bzero;
    if (verbose && !scaling_warned) {
      if (orig_array->bscale != oseg->bscale ||
	  orig_array->bzero != oseg->bzero) {
	printf
	    ("Warning: data BSCALE/BZERO differ from output series defaults\n");
	printf ("input BSCALE = %f BZERO = %f\n", orig_array->bscale,
	    orig_array->bzero);
	printf ("default output BSCALE = %f BZERO = %f\n", oseg->bscale,
	    oseg->bzero);
	scaling_warned = 1;
      }
    }
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
    
		/*  why does this code have to replicate that at line 604?  */
			/*  because out_axes was re-mallocd (without free)  */
    if (collapse) {
      DRMS_Array_t *coll_array;
      int *new_axes = (int *)malloc (rank * sizeof (int));
      int ccrank = rank;
      for (n = 0, m = 0; n < rank; n++) {
        if (out_axes[n] == 1) ccrank--;
	else new_axes[m++] = out_axes[n];
      }
					   /*  do not collapse to zero rank  */
      if (ccrank < 1) {
	ccrank = 1;
	new_axes[0] = 1;
      }
      if (ccrank == rank) collapse = 0;
      else for (n = 0; n < ccrank; n++) out_axes[n] = new_axes[n];
      free (new_axes);
      coll_array = drms_array_create (DRMS_TYPE_DOUBLE, ccrank, out_axes,
	  (void *)vbin, &status);
      coll_array->bscale = binned_array->bscale;
      coll_array->bzero = binned_array->bzero;
      drms_segment_write (oseg, coll_array, 0);
    } else if (merge) {
      drms_segment_writeslice (oseg, binned_array, oslice_strt, oslice_stop, 0);
    } else {
      kstat += check_and_set_key_int (orec, "Quality", qual_out);
      drms_segment_write (oseg, binned_array, 0);
    }
							 /*  parameter keys  */
    if (!merge) {
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
      if (!merge) kstat += check_and_set_key_str (orec, "OffsetBy", valuestr);
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
      if (!merge) kstat += check_and_set_key_str (orec, "ScaledBy", valuestr);
      if (verbose) {
        printf ("  scaled by ");
	if (scale_mult) printf ("%s = ", valuestr);
	printf ("%g\n", scale);
      }
    }
							/*  statistics keys  */
    if (!merge) {
      kstat += set_stats_keys (orec, binned_array, ntbin);
      drms_close_record (orec, dispose);
    }
  }
  drms_close_records (drs, DRMS_FREE_RECORD);
  if (merge) drms_close_record (orec, dispose);

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
 *  v 1.0	Added stubs for optional pre-filtering
 *		Added collapse option
 *		Fixed setting of array bzero
 *  v 1.0 frozen 2011.11.15
 *  v 1.1	Added option for merging output records along principal axis
 *		Skip consistency check for collapsed axes
 *		Fixed axis collapse code for case in which multiple axes are
 *	collapsed
 *		Default output FITS scaling same as input; no way to override,
 *	but warning issued if differs from output segment default
 *		Fixed check for existence of input records
 *		Added version control for included local source files
 *  v 1.1 frozen 2018.08.20
 */
