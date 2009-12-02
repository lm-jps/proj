/*
 *  pspec3.c                                         	~rick/src/pspec
 *
 *  Calculate the 3-dimensional power spectrum of an input dataset
 *  The module expects a 3-dimensional real dataset
 *
 *  Responsible:
 *      Rick Bogart                             RBogart@solar.stanford.edu
 *
 *  Usage:
 *    pspec3 [-lvx] in= pspec= ...
 *
 *  Arguments: (type   default         description)
 *    in	DataSet	-		Input dataset
 *    out	DataSeries -	Output dataseries
 *    mask_in	Float	0.9375		Inner apodization radius
 *    mask_ex	Float	1.0		Outer apodization radius
 *    apodize	Float	0.96875		Temporal apodization edge
 *    fbin	int	0		output frequency binning
 *
 *  Flags
 *	-l	output scaled log of power spectrum
 *	-v	run verbose
 *	-x	use double-precision calculation
 *
 *  Bugs:
 *    Only works with a single data record on input
 *    Requires that the input data segment be named V
 *    No checking if parameters consistent with (constant) values of output
 *	series
 *
 *  Future Updates
 *    0.8 clean up code, remove CGI-output
 *    Use drms_keyword_lookup()
 *    Create new data series as needed
 *
 *  Revision history is at end of file.
 */
#include <jsoc_main.h>
#include <fftw3.h>

char *module_name = "pspec3";
char *module_desc = "3-d power spectrum";
char *version_id = "0.7";

ModuleArgs_t module_args[] = {
/*
  {ARG_NUME,	"in",	"mdi.hmiVtrack16", "input data series: tracked cubes",
      "mdi.hmiVtrack32, mdi.hmiVtrack16, mdi.hmiVtrack5, mdi.Vtrack"}, 
  {ARG_STRING,	"record", "Not Specified", "input data record"}, 
*/
  {ARG_DATASET, "in", "mdi.rdVtrack_dp[?CarrRot=1988 and LonHG=180.0 and LatHG=0.0?]", "Input data record"},
  {ARG_DATASERIES, "pspec", "mdi.rdVpspec_dp", "Ouput data series"},
  {ARG_FLOAT, "mask_in", "0.9375", "inner radial apodization edge", "[0.0,)"},
  {ARG_FLOAT, "mask_ex", "1.0", "outer radial apodization edge", "[0.0,)"},
  {ARG_FLOAT,	"apodize", "0.96875", "temporal apodization edge",
	"[0.0,1.0]"},
  {ARG_INT,	"fbin", "0", "output frequency fbins (0-> none)"},
  {ARG_FLAG,	"l",	"",
	"output direct power spectrum rather than scaled log"},
  {ARG_FLAG,	"x", 	"", "use double-precision calculation"},
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {}
};

/*
 *  queue_stuff.c					~rick/hmi/rings/src
 *
 *  Functions to support queueing of remote requests
 *     
 *  Responsible:  Rick Bogart			RBogart@solar.stanford.edu
 *
 *  Bugs:
 *    For generating a queueing script, the module_name must be the same as
 *	name of the binary file
 *    Probably all arg values should be protected by quotes when writing to
 *	the script
 *
 *  Revision history is at end of file.
 */
#define QUEUE_DIR	("/tmp/qreq")
				      /*  needed for rsb_write_queue_script  */
static void str_compress (char *s) {
  char *source, *dest;
  source = dest = s;
  while (*source) {
    if (!isspace (*source)) {
      *dest = *source;
      dest++;
    }
    source++;
  }
  *dest = 0;
}
				      /*  needed for rsb_write_queue_script  */
static int parse_numerated (char *klist, char ***names) {
/*
 *  Parses an entry in the args list of type ARG_NUME and returns the number
 *    of possible values, determined from the range entry interpreted as a
 *    comma-separated set of strings, and a malloc'd array of strings
 *    corresponding to the enumeration choices
 */
  int found = 0, maxlen;
  char *next, *nptr, *tmp;

  if (!klist) return found;
  maxlen = strlen (klist);
  tmp = malloc (maxlen + 1);
  strcpy (tmp, klist);
  str_compress (tmp);				     /*  remove white space  */
  next = tmp;

  *names = (char **)malloc (maxlen * sizeof (char *));
  while (next) {
    nptr = next;
    if (next = (char *)strchr (nptr, ',')) {
      *next = 0;
      next++;
    }
    if (!strlen (nptr)) continue;
					      /*  nptr now points to entity  */
    (*names)[found] = (char *)malloc (strlen (nptr) + 1);
    strcpy ((*names)[found], nptr);
    found++;
  }
  free (tmp);

  return found;
}

int rsb_write_queue_script (const char *modname, int verbose) {
  ModuleArgs_t *arg = module_args;
  CmdParams_t *params = &cmdparams;
  FILE *script;
  int fd, flags = 0;
  char *notify, *scriptid;
  char template[] = "/tmp/qreq/qsub.XXXXXX";
	   /*  tempnam overrides dir if TMPDIR environment variable is set!  */
  fd = mkstemp (template);
  script = fdopen (fd, "w");
  if (!script) {
    fprintf (stderr, "Error, unable to generate queueing script %s\n",
	template);
    return 0;
  }
  
  fprintf (script, "#!/bin/csh -f\n");
  scriptid = template + strlen (template) - 6;
  fprintf (script, "set QID = %s\n", scriptid);
  fprintf (script, "set ODIR = /usr/local/www/htdocs/qres\n");
  fprintf (script, "set WWW = http://rick.stanford.edu/qres\n");
  fprintf (script, "set OUT = $ODIR/qres.%s.html\n", scriptid);

  fprintf (script, "/home/rick/bin/linux_x86_64/%s ", modname);
  while (arg->type != ARG_END) {
    if (arg->type == ARG_FLAG) {
      if (params_isflagset (params, arg->name)) {
	if (!flags) fprintf (script, "-");
	fprintf (script, "%s ", arg->name);
	flags++;
      }
    }
    arg++;
  }
  arg = module_args;
  while (arg->type != ARG_END) {
    if (arg->type == ARG_FLAG) {
      arg++;
      continue;
    }
    if (arg->type == ARG_STRING || arg->type == ARG_INTS ||
	arg->type == ARG_INTS || arg->type == ARG_DATASET)
      fprintf (script, " \\\n\t%s= \"%s\"", arg->name,
	  params_get_str (params, arg->name));
    else if (arg->type == ARG_NUME) {
      char **names;
      int nval = params_get_int (params, arg->name);

      parse_numerated (arg->range, &names);
      fprintf (script, " \\\n\t%s= %s", arg->name, names[nval]);
    } else  fprintf (script, " \\\n\t%s=  %s", arg->name,
	params_get_str (params, arg->name));
    arg++;
  }
  fprintf (script, " \\\n\tHTML_Output= 1 >& $OUT\n");
  if (notify = params_get_str (params, "CGI_email_notify")) {
    fprintf (script, "set MAIL = /bin/mail\n");
    fprintf (script, "set MSG = /tmp/qreq/qnot.%s\n", scriptid);
    fprintf (script, "echo \"Your JSOC process has completed\" > $MSG\n");
    fprintf (script, "echo \"Output can be found at $WWW/qres.%s.html\" >> $MSG\n",
	scriptid);
    fprintf (script, "$MAIL -s \"JSOC processing completed\" %s < $MSG\n",
	notify);
    fprintf (script, "rm $MSG\n");
  }

  fclose (script);
  chmod (template, 00777);
  printf ("Your module execution request has been submitted to the queue\n");
  printf ("When it is completed, the output can be found ");
  printf ("<A HREF= \"http://rick.stanford.edu/qres/qres.%s.html\"<B>here</B></A>\n",
      scriptid);
  if (verbose) printf ("script to be queued is: %s\n", template);
  return 0;
}

static void drms_segment_rename_file (DRMS_Segment_t *segment) {
  char path[DRMS_MAXPATHLEN+1], old[2*DRMS_MAXPATHLEN+3],
      new[DRMS_MAXPATHLEN+DRMS_MAXSEGNAMELEN+3];
  drms_record_directory (segment->record, path, 1);
  sprintf (old, "%s/%s", path, segment->filename);
  sprintf (new, "%s/V.fits", path);
  rename (old, new);
  strcpy (segment->filename, "V.fits");
}

static int rsb_get_max_keyvalue_int (const char *series, const char *key, int *status) {
  DB_Text_Result_t *sqres;
  int val = -1;
  char query[DRMS_MAXQUERYLEN];

  sprintf (query, "select max (%s) from %s", key, series);
  if ((sqres = drms_query_txt (drms_env->session, query)))
    val = atol (sqres->field[0][0]);
  if (status) *status = -1;
  return val;
}

static void PrintTrailer () {
  printf ("\n</PRE><DIV ALIGN=CENTER><HR>\n");
  printf ("<B>Please do not use the back button of your browser,\n");
  printf ("as it will requeue the module!</B><BR>\n");
  printf ("To return to the module executor\n");
  printf ("<A HREF= \"http://rick.stanford.edu/cgi-bin/%s\"><B>Click Here</B></A>\n",
      module_name);
  printf ("<HR>\n");
  printf ("<I>Home pages for:</I>\n");
  printf ("<A HREF = \"http://jsoc.stanford.edu/\">");
  printf ("<B>SDO-JSOC</B></A>\n");
  printf ("<B>|</B>\n");
  printf ("<A HREF = \"http://sha.stanford.edu/\">");
  printf ("<B>SHA</B></A>\n");
  printf ("<HR></DIV>\n");
  printf ("</BODY>\n");
  printf ("</HTML>\n");
  fflush (stdout);
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *irecs;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *orig, *pspec;

  double *apodization;
  double *wdata;
  double crpix, crval = 0.0;
  double dval, r, r2, ra, ra2, rx, ry, r0x, r0y, t, weight;
  float *data, *xdata;
  float ftrb, ftib, fval, vmin, vmax;
  int axes[3];
  int rank, col, cols, row, rows, plane, planes, area;
  int hcols, hrows, hplanes, opln, xcols, cols_even, rows_even;
  int segs;
  int bin, is, js, l, m, n, ntot;
  int status;
  char pathname[2*DRMS_MAXPATHLEN+1];
  char odsname[DRMS_MAXSERIESNAMELEN + 13];
/*
  char comment[64], comment0[64];
*/
  void *val;

  char *inds = params_get_str (params, "in");
  char *out_series = params_get_str (params, "pspec");
  double apode_edge = params_get_double (params, "apodize");
  double edge_inner = params_get_double (params, "mask_in");
  double edge_outer = params_get_double (params, "mask_ex");
  int fbins = params_get_int (params, "fbin");
  int log_out = params_isflagset (params, "l") ? 0 : 1;
  int dp_calc = params_isflagset (params, "x");
  int verbose = params_isflagset (params, "v");
  int forexport = (strcmp (out_series, "jsoc.exports")) ? 0 : 1;
  int from_cgi = cmdparams_exists (params, "CGI_Input");
  int html_format = cmdparams_exists (params, "HTML_Output");
/*
  int num_ds = getkey_int (params, "in_nsets");
  int new_protocol =
      (strcmp (GETKEY_str (params, "protocol"), "RDB.FITS_MERGE") == 0) ?
      VDS_FITS_MERGE : VDS_FITS;
*/
                                        /* Declaration for FFTW library call */
  fftwf_plan fplan;
  fftw_plan plan;

  if (verbose) {
    printf ("%s V %s\n", module_name, version_id);
    if (fbins > 1)
      printf ("  output binned by %d frequencies per plane\n", fbins);
  }
  irecs = drms_open_records (drms_env, inds, &status);
  if (status) {
    fprintf (stderr, "Error: drms_open_records returned %d for dataset:\n", status);
    fprintf (stderr, "  %s\n", inds);
    return 0;
  }

  if (irecs->n > 1) {
    fprintf (stderr, "Query produced %d matching records; must limit to 1\n",
	irecs->n);
    drms_close_records (irecs, DRMS_FREE_RECORD);
    return 0;
  }
  if (irecs->n > 1) {
    fprintf (stderr, "No records found in input data set %s\n", inds);
    drms_close_records (irecs, DRMS_FREE_RECORD);
    return 0;
  }
  irec = irecs->records[0];
  if (!(iseg = drms_segment_lookup (irec, "V"))) {
    fprintf (stderr, "Error, could not find segment \"V\"\n");
    drms_close_records (irecs, DRMS_FREE_RECORD);
    return 0;
  }
  if (iseg->info->naxis != 3) {
    fprintf (stderr, "Error: rank of data cube (%d) != 3\n", iseg->info->naxis);
    drms_close_records (irecs, DRMS_FREE_RECORD);
    return 0;
  }

  if (from_cgi) {
	      /*  if the job is to be queued, write out the script and exit  */
    drms_close_records (irecs, DRMS_FREE_RECORD);
    rsb_write_queue_script (module_name, verbose);
    return 0;
  }
			   /*  the input looks okay, open the output record  */
  orec = drms_create_record (drms_env, out_series, DRMS_PERMANENT, &status);
  if (status) {
    fprintf (stderr, "Error: drms_create_record returned %d for data series:\n",
	status);
    fprintf (stderr, "  %s\n", out_series);
    drms_close_records (irecs, DRMS_FREE_RECORD);
    return 0;
  }
  cols = iseg->axis[0];
  rows = iseg->axis[1];
  planes = iseg->axis[2];
  area = cols * rows;
  xcols = 2 * ((cols + 2) / 2);
  hcols = cols / 2;
  hrows = rows / 2;
  cols_even = (cols % 2) ? 0 : 1;
  rows_even = (rows % 2) ? 0 : 1;
  hplanes = planes / 2;			     /*  possible bug if planes odd  */
  if (planes % 2) hplanes++;

  apodization = (double *)malloc (area * sizeof (double));
  if (edge_outer == 0.0) {
    for (n = 0; n < area; n++)
      apodization[n] = 1.0;
  } else {
    r0x = hcols + 0.5 * cols_even;
    r0y = hrows + 0.5 * rows_even;
    for (row=0, n=0; row<rows; row++) {
      ry = (row - hrows + 0.5 * rows_even) / r0y;
      for (col=0; col<cols; col++) {
        rx = (col - hcols + 0.5 * cols_even) / r0x;
	r2 = rx * rx + ry * ry;
	r = sqrt (r2);
        if (edge_inner >= edge_outer)
	  ra = (r < edge_outer) ? 0.0 : 1.0;
	else
	  ra = (r - edge_inner) / (edge_outer - edge_inner);
	ra2 = 1.0 - ra * ra;
							/*  fixed  03.07.17  */
        apodization[n++] = (ra >= 1.0) ? 0.0 : (ra <= 0.0) ? 1.0 : ra2 * ra2;
      }
    }
  }
					   /*  Initialize FFT working space  */
  if (dp_calc) {
    wdata = (double *)malloc (planes * rows * xcols * sizeof (double));
    plan = fftw_plan_dft_r2c_3d (planes, rows, cols, wdata,
	(fftw_complex *)wdata, FFTW_ESTIMATE);
  } else {
    xdata = (float *)malloc (planes * rows * xcols * sizeof (float));
    fplan = fftwf_plan_dft_r2c_3d (planes, rows, cols, xdata,
	(fftwf_complex *)xdata, FFTW_ESTIMATE);
  }
                                                 /*  Loop on input datasets  */
/*
  for (ds = 0; ds < num_ds; ds++) {
*/
                                                   /*  open output data set  */
/*
*/
				    /*  add the argument values to header  */
/*
    sds_append_args_tohist (vout->global_attributes, arguments, params);
    if (sds_append_attrs (vout->global_attributes, vin->global_attributes))
      (*errlog) ("WARNING - problem copying overview file required info\n");
    vds_set_attribute (vout, "CONFORMS", "MISC", SDS_STRING,
        "2-d mosaic of transformed dimension 3-cubes");
    vds_set_attribute (vout, "DATAKIND", "3-d power spectrum", SDS_STRING, "");
    fsn = getkey_int (params, key_fsn);
    lsn = getkey_int (params, key_lsn);
    if (lsn == -1) lsn = vds_last_record (vin);
    dsn = getkey_int (params, key_dsn);
*/
                               /*  Loop through records within each data set */
/*
    for (sn = fsn; sn <= lsn; sn+=dsn) {
      if (verbose) (*history) ("Begin %d:%d\n", ds, sn); 
      cube = (no_vds) ?
	  sds_in (fitsname (params, ids, sn), SDS_FLOAT, SDS_FITS) :
	  VDS_select_frec (vin, 0, sn);
*/
/*
      if (!cube) {
	if (no_vds) (*errlog) ("Error reading file %s) - skipped\n",
	    fitsname (params, ids, sn));
	else (*errlog)
	    ("Error reading record #%d - skipped\n", sn);
	(*errlog) ("  soi_errno = %d\n", soi_errno);
	continue;
      }
*/
// printf ("calling drms_segment_read()\n");
    orig = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
// printf ("drms_segment_read() returned %d\n", status);
    if (status) {
      fprintf (stderr, "Error on drms_segment_read\n");
      return 0;
    }
    l = m = n = 0;
    data = (float *)orig->data;
      if (planes < 2) apode_edge = 2.0;
      for (plane = 0; plane < planes; plane++) {
        if (apode_edge < 1.0) {
	  t = (2.0 * plane) / (planes - 1.0);
	  if (t > 1.0) t = 2.0 - t;
	  t = 1.0 - t;
	  if (t > apode_edge) {
	    t = (t - apode_edge) / (1.0 - apode_edge);
	    t *= t;
	    weight = (1.0 - t);
	    weight *= weight;
	  } else weight = 1.0;
	} else weight = 1.0;
	if (dp_calc) {
	  for (row = 0; row < rows; row++) {
	    for (col = 0; col < cols; col++) {
	      dval = data[m++];
	      if (isnan (dval)) dval = 0.0;
	      wdata[n++] = weight * dval * apodization[l++];
            }
	    n += xcols - cols;
	  }
	} else {
	  for (row = 0; row < rows; row++) {
	    for (col = 0; col < cols; col++) {
	      fval = data[m++];
	      if (isnan (fval)) fval = 0.0;
	      xdata[n++] = weight * fval * apodization[l++];
            }
	    n += xcols - cols;
	  }
        }
	l = 0;
      }
                                              /*  forward Fourier transform  */
  if (dp_calc) fftw_execute (plan);
  else fftwf_execute (fplan);
						 /*  Create output data set  */
  axes[0] = cols;
  axes[1] = rows;
  opln = (fbins > 1) ? hplanes / fbins : hplanes;
  if (fbins > 1 && hplanes % fbins) opln++;
  axes[2] = opln;		             /*  possible bug if planes odd  */
  ntot = cols * rows * opln;
  data = (float *)malloc (ntot * sizeof (float));
  pspec = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, data, &status);
  if (status) {
    fprintf (stderr, "Error: couldn't create output array\n");
    return 0;
  }
						/*  copy keywords from input  */
  drms_copykey (orec, irec, "CarrTime");
  drms_copykey (orec, irec, "CarrRot");
  drms_copykey (orec, irec, "CarrLon");
  drms_copykey (orec, irec, "LonHG");
  drms_copykey (orec, irec, "LatHG");
  drms_copykey (orec, irec, "CMLon");
  drms_copykey (orec, irec, "MidTime");
  drms_copykey (orec, irec, "Duration");
  drms_copykey (orec, irec, "LonSpan");
  drms_copykey (orec, irec, "T_START");
  drms_copykey (orec, irec, "T_STOP");
  drms_copykey (orec, irec, "Coverage");
  drms_copykey (orec, irec, "ZonalTrk");
  drms_copykey (orec, irec, "MeridTrk");
  drms_copykey (orec, irec, "MapScale");
  drms_copykey (orec, irec, "MapProj");
  drms_copykey (orec, irec, "PosAng");
  drms_copykey (orec, irec, "MAI");
						       /*  and set new ones  */
  crpix = 0.5 * (cols + 1);
  drms_setkey_float (orec, "CRPIX1", crpix);
  if (cols != rows) crpix = 0.5 * (rows + 1);
  drms_setkey_float (orec, "CRPIX2", crpix);
  crpix = 1.0;
  drms_setkey_float (orec, "CRPIX3", crpix);
  drms_setkey_float (orec, "CRVAL1", crval);
  drms_setkey_float (orec, "CRVAL2", crval);
  drms_setkey_float (orec, "CRVAL3", crval);
  drms_setkey_string (orec, "CTYPE1", "inv-Mm");
  drms_setkey_string (orec, "CTYPE2", "inv-Mm");
  drms_setkey_string (orec, "CTYPE2", "inv-sec");

  dval = drms_getkey_double (irec, "MapScale", &status);
  if (!status) {
						    /*  assume cols = rows!  */
    dval = 360.0 / dval / cols;				 /*  k (inv radian)  */
    dval /= 696.0;					     /*  k (inv Mm)  */
    drms_setkey_float (orec, "CDELT1", dval);
    if (cols == rows) {
      drms_setkey_double (orec, "CDELT2", dval);
      drms_setkey_float (orec, "Delta_k", dval);
    } else {
      drms_setkey_float (orec, "Delta_kx", dval);
      dval *= (double)cols / (double)rows;
      drms_setkey_float (orec, "CDELT2", dval);
      drms_setkey_float (orec, "Delta_ky", dval);
    }
  }
  dval = (double)drms_getkey_time (irec, "T_STOP", &status);
  dval -= (double)drms_getkey_time (irec, "T_START", &status);
  fval = drms_getkey_float (irec, "CDELT3", &status);
  if (!status) dval += fval;
  if (dval <= 0.0)
    fprintf (stderr, "** Warning: t(stop) <= (tstart); ignored\n");
  else {
    dval = (1.0e6 - (1.0e6 / planes)) / dval;
    if (fbins > 1) dval *= fbins;
    drms_setkey_float (orec, "Delta_nu", dval);
    dval *= 2.0e-6 * M_PI;
    drms_setkey_float (orec, "D_OMEGA", dval);
    drms_setkey_float (orec, "CDELT3", dval);
  }
  drms_setkey_float (orec, "apode_k_min", edge_inner);
  drms_setkey_float (orec, "apode_k_max", edge_outer);
  drms_setkey_float (orec, "apode_f", apode_edge);
		/*  convert complex transform to power spectrum and reorder  */
  vmax = vmin = data[0];
  if (fbins > 1) {
    double nfac = 1.0 / fbins;
    for (plane = 0, n = 0; plane < hplanes; plane += fbins) {
      for (row = 0; row < rows; row++) {
	for (col = 0; col < cols; col++, n++) {
	      js = (row == row % hrows) ? row + hrows : row - hrows;
	      is = (col == col % hcols) ? col + hcols : col - hcols;
	      ftrb = ftib = 0.0;
	      for (bin = 0; bin < fbins; bin++) {
		m = (plane + bin) * xcols * rows + js * xcols + 2 * is;
		if (is >= hcols) {
		  m = 2 * (cols - is);
		  if (js != 0) m += (rows - js) * xcols;
		  if ((plane + bin) != 0)
		    m += (planes - plane - bin) * xcols * rows;
		}
		if (dp_calc) {
		  ftrb += wdata[m];
		  ftib += wdata[m+1];
		} else {
		  ftrb += xdata[m];
		  ftib += xdata[m+1];
	        }
	      }
	      data[n] = nfac * (ftrb*ftrb + ftib*ftib);
	      if (log_out) {
		if (data[n] > vmax) vmax = data[n];
		if (data[n] < vmin) vmin = data[n];
	      }
	}
      }
    }
    if (n > ntot) {
      fprintf (stderr, "** Error: data written beyond output array bounds\n");
      return 1;
    }
  } else {
    for (plane = 0, n = 0; plane < hplanes; plane++) {
      for (row = 0; row < rows; row++) {
	for (col = 0; col < cols; col++, n++) {
	      js = (row == row % hrows) ? row + hrows : row - hrows;
	      is = (col == col % hcols) ? col + hcols : col - hcols;
	      m = plane * xcols * rows + js * xcols + 2 * is;
	      if (is >= hcols) {
		m = 2 * (cols - is);
		if (js != 0) m += (rows - js) * xcols;
		if (plane != 0) m += (planes - plane) * xcols * rows;
	      }
	      data[n] = (dp_calc) ?
	        wdata[m]*wdata[m] + wdata[m+1]*wdata[m+1] :
	        xdata[m]*xdata[m] + xdata[m+1]*xdata[m+1];
	      if (log_out) {
		if (data[n] > vmax) vmax = data[n];
		if (data[n] < vmin) vmin = data[n];
	      }
	}
      }
    }
  }
  if (forexport) {
    int reqid = rsb_get_max_keyvalue_int (out_series, "RequestID", &status);
    char *notify;

    oseg = drms_segment_lookup (orec, "Data");
    drms_segment_rename_file (oseg);
    drms_setkey_time (orec, "ReqTime", CURRENT_SYSTEM_TIME);
    drms_setkey_string (orec, "Requestor", getlogin ());
    reqid++;
    drms_setkey_int (orec, "RequestID", reqid);
    if (notify = params_get_str (params, "CGI_email_notify"))
      drms_setkey_string (orec, "Notify", notify);
    sprintf (odsname, "%s[%d]", out_series, reqid);
  } else {
    if ((segs = orec->segments.num_total) < 1) {
      fprintf (stderr, "Error: no data segments in output data series:\n");
      fprintf (stderr, "  %s\n", out_series);
      drms_close_record (orec, DRMS_FREE_RECORD);
      drms_close_records (irecs, DRMS_FREE_RECORD);
      return 0;
    }
    for (n = 0; n < segs; n++) {
      oseg = drms_segment_lookupnum (orec, n);
      if (oseg->info->naxis == 3) break;
    }
    if (n >= segs && verbose)
      printf ("found no segmemt of rank 3, using segment %d\n", n);
    sprintf (odsname, "%s", out_series);
  }

  if (log_out) {
    double bscale, bzero;
    if (verbose) printf ("  values range from %.3e to %.3e;\n", vmin, vmax);
    if (vmin < 1.0e-10) {
      fprintf (stderr,
	  "** Warning: minimum value = %.2e; value = log (1 + data)\n", vmin);
      for (n = 0; n < ntot; n++) data[n] = log (1.0 + data[n]);
      vmin = log (1.0 + vmin);
      vmax = log (1.0 + vmax);
    } else {
      for (n = 0; n < ntot; n++) data[n] = log (data[n]);
      vmin = log (vmin);
      vmax = log (vmax);
    }
    if (verbose) printf ("  logs scaled between %f and %f\n", vmin, vmax);
/*
    drms_segment_autoscale (oseg, pspec);
*/
    bscale = (vmax - vmin) / 65000.0;
    bzero = 0.5 * (vmax + vmin);
    pspec->bscale = bscale;
    pspec->bzero = bzero;
/*
    if (drms_segment_setscaling (oseg, bzero, bscale))
      fprintf (stderr, "Warning: scaling failed, precision may be lost\n");
*/
    drms_setkey_double (orec, "Log_Base", exp (1.0));
  }
  drms_setkey_float (orec, "DataMin", vmin);
  drms_setkey_float (orec, "DataMax", vmax);

  status = drms_segment_write (oseg, pspec, 0);
  drms_segment_filename (oseg, pathname);
  if (html_format) {
    printf ("output data set is <B>%s</B>\n", odsname);
    printf ("output data cube at <A HREF= \"%s\"><B>%s</B></A>\n", pathname, pathname);
  } else if (verbose) {
    printf ("output data set is %s\n", odsname);
    printf ("output data cube at %s\n", pathname);
  }

  drms_close_record (orec, DRMS_INSERT_RECORD);
                                                      /* Clean up FFTW stuff */
  fftwf_destroy_plan (fplan);
  drms_close_records (irecs, DRMS_FREE_RECORD);
  if (html_format) PrintTrailer ();
  status = 0;
  return status;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  08.03.24	created this file, based on SOI powrspec3_v20
 *  v 0.7	fix to run with NetDRMS 2.0 (09.07.14);
 *		Modify and add to list of keywords copied from input
 *
 */
