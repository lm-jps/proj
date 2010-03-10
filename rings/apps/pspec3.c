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
 *  Arguments: (type		default         description)
 *    in	DataSet		-		Input dataset
 *    pspec	DataSeries	-		Output dataseries
 *    mask_in	Float		0.9375		Inner apodization radius
 *    mask_ex	Float		1.0		Outer apodization radius
 *    apodize	Float		0.96875		Temporal apodization edge
 *    fbin	int		0		output frequency binning
 *
 *  Flags
 *	-l	output direct tpower spectrum rather than scaled log
 *	-n	do not save results (diagnostic only)
 *	-v	run verbose
 *	-x	use double-precision calculation
 *
 *  Bugs:
 *    Requires that the input data segment be named V
 *    No checking that input data series corresponds with link in output
 *	series
 *    The value of CDELT3 in the input set overrides the value of Cadence
 *	if propagated, but is not checked for units
 *
 *  Future Updates
 *    Use drms_keyword_lookup()
 *    Create new data series as needed
 *    Parallelize over multiple records
 *    Add recreation option for target output records
 *    Add parameters for overriding essential input keywords
 *    Add parameter for overriding keyword list for propagation
 *    Modify to also work with a series of two-dimensional images - maybe
 *	not such a good idea
 *
 *  Revision history is at end of file.
 */
#include <jsoc_main.h>
#include <fftw3.h>

char *module_name = "pspec3";
char *module_desc = "3-d power spectrum";
char *version_id = "0.8";

ModuleArgs_t module_args[] = {
  {ARG_DATASET, "in", "", "Input data record"},
  {ARG_DATASERIES, "pspec", "", "Ouput data series"},
  {ARG_FLOAT, "mask_in", "0.9375", "inner radial apodization edge", "[0.0,)"},
  {ARG_FLOAT, "mask_ex", "1.0", "outer radial apodization edge", "[0.0,)"},
  {ARG_FLOAT,	"apodize", "0.96875", "temporal apodization edge",
	"[0.0,1.0]"},
  {ARG_INT,	"fbin", "0", "output frequency fbins (0-> none)"},
  {ARG_FLAG,	"l",	"",
	"output direct power spectrum rather than scaled log"},
  {ARG_FLAG,    "n", "0", "do not save output record (diagnostics only)"},      
  {ARG_FLAG,	"x", 	"", "use double-precision calculation"},
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {}
};
       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "LonSpan", "T_START", "T_STOP", "Coverage",
    "Cadence", "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel", "MapScale",
    "Size", "MapProj", "Map_PA", "PosAng", "RSunRef", "MAI"};

#include "keystuff.c"

static int cleanup (int error, DRMS_RecordSet_t *irecs, DRMS_Record_t *orec,
    DRMS_Array_t *orig, DRMS_Array_t *pspec, int dispose) {
  if (orig) drms_free_array (orig);
  if (pspec) drms_free_array (pspec);
  if (irecs) drms_close_records (irecs, DRMS_FREE_RECORD);
  if (orec) {
    if (error) drms_close_record (orec, DRMS_FREE_RECORD);
    else drms_close_record (orec, dispose);
  }
  return error;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *irecs = NULL;
  DRMS_Record_t *irec, *orec = NULL;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *orig = NULL, *pspec = NULL;

  double *apodization;
  double *wdata;
  double crpix, crval = 0.0;
  double dnu, domega, dval, r, r2, ra, ra2, rx, ry, r0x, r0y, t, tstep;
  double dataavg, datavar, powrint, normal, weight;
  double bzero, bscale, scale_range;
  float *data, *xdata;
  float ftrb, ftib, fval, vmin, vmax;
  int axes[3];
  int rank, col, cols, row, rows, plane, planes, area;
  int hcols, hrows, hplanes, opln, xcols, cols_even, rows_even;
  int rgn, rgnct, segs;
  int bin, is, js, l, m, n, ntot;
  int kstat, status;
  char *input;
  char module_ident[64];
  char pathname[2*DRMS_MAXPATHLEN+1];
  char source[DRMS_MAXQUERYLEN], recid[DRMS_MAXQUERYLEN];
  void *val;

  int keyct = sizeof (propagate) / sizeof (char *);

  char *inds = params_get_str (params, "in");
  char *out_series = params_get_str (params, "pspec");
  double apode_edge = params_get_double (params, "apodize");
  double edge_inner = params_get_double (params, "mask_in");
  double edge_outer = params_get_double (params, "mask_ex");
  int fbins = params_get_int (params, "fbin");
  int log_out = params_isflagset (params, "l") ? 0 : 1;
  int dp_calc = params_isflagset (params, "x");
  int verbose = params_isflagset (params, "v");
  int no_save = params_isflagset (params, "n");
  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
                                        /* Declaration for FFTW library call */
  fftwf_plan fplan;
  fftw_plan plan;

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) {
    printf ("%s:\n", module_ident);
    if (fbins > 1)
      printf ("  output binned by %d frequencies per plane\n", fbins);
    if (no_save)
      printf ("(diagnostic run only, no records will be written to DRMS)\n");
  }
  irecs = drms_open_records (drms_env, inds, &status);
  if (status) {
    fprintf (stderr, "Error (%s): drms_open_records() returned %d for dataset:\n",
	module_ident, status);
    fprintf (stderr, "  %s\n", inds);
    return cleanup (1, irecs, orec, orig, pspec, dispose);
  }

  rgnct = irecs->n;
  if (rgnct < 1) {
    fprintf (stderr, "No records found in input data set %s\n", inds);
    return cleanup (1, irecs, orec, orig, pspec, dispose);
  }
  input = strdup (inds);
  for (rgn = 0; rgn < rgnct; rgn++) {
    irec = irecs->records[rgn];
    drms_sprint_rec_query (source, irec);

    if (!(iseg = drms_segment_lookup (irec, "V"))) {
      fprintf (stderr, "Warning: could not find segment \"V\"\n");
      fprintf (stderr, "       in %s; skipped\n", source);
      continue;
    }
    if (iseg->info->naxis != 3) {
      fprintf (stderr, "Error: rank of data cube (%d) != 3\n", iseg->info->naxis);
      return cleanup (1, irecs, orec, orig, pspec, dispose);
    }
    if (drms_wcs_timestep (irec, 3, &tstep)) {
      dval = (double)drms_getkey_time (irec, "T_STOP", &status);
      tstep = dval - (double)drms_getkey_time (irec, "T_START", &status);
      if (tstep <= 0.0) {
	fprintf (stderr,
	    "Error: insufficient key information for frequency determination\n");
	fprintf (stderr,
	    "       in record %s; skipped\n", source);
	continue;
      }
      tstep /= planes - 1.0;
    }
			   /*  the input looks okay, open the output record  */
    orec = drms_create_record (drms_env, out_series, DRMS_PERMANENT, &status);
    if (status) {
      fprintf (stderr,
          "Error: drms_create_record returned %d for data series:\n", status);
      fprintf (stderr, "       %s\n", out_series);
      return cleanup (1, irecs, orec, orig, pspec, dispose);
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
    hplanes = planes / 2;		     /*  possible bug if planes odd  */
    if (planes % 2) hplanes++;
 
    apodization = (double *)malloc (area * sizeof (double));
    if (edge_outer == 0.0) {
      for (n = 0; n < area; n++) apodization[n] = 1.0;
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
				    /*  add the argument values to header  */
    orig = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
    if (status) {
      fprintf (stderr, "Error on drms_segment_read in\n");
      fprintf (stderr, "      %s\n", source);
      return cleanup (1, irecs, orec, orig, pspec, dispose);
    }
    l = m = n = 0;
    dataavg = datavar = 0.0;
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
	      dval *= weight * apodization[l++];
	      dataavg += dval;
	      datavar += dval * dval;
	      wdata[n++] = dval;
          }
	  n += xcols - cols;
	}
      } else {
	for (row = 0; row < rows; row++) {
	  for (col = 0; col < cols; col++) {
	      fval = data[m++];
	      if (isnan (fval)) fval = 0.0;
	      fval *= weight * apodization[l++];
	      dataavg += fval;
	      datavar += fval * fval;
	      xdata[n++] = fval;
          }
	  n += xcols - cols;
	}
      }
      l = 0;
    }
    drms_free_array (orig);
						  /*  necessary for cleanup  */
    orig = NULL;
    dataavg /= m;
    datavar /= m;
    datavar -= dataavg * dataavg;
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
				    /*  copy designated keywords from input  */
    kstat = 0;
    kstat = propagate_keys (orec, irec, propagate, keyct);
    dval = drms_getkey_double (irec, "CDELT3", &status);
    if (!status) kstat += check_and_set_key_float (orec, "Cadence", dval);
						       /*  and set new ones  */
    crpix = 0.5 * (cols + 1);
    kstat += check_and_set_key_float (orec, "CRPIX1", crpix);
    if (cols != rows) crpix = 0.5 * (rows + 1);
    kstat += check_and_set_key_float (orec, "CRPIX2", crpix);
    crpix = 1.0;
    kstat += check_and_set_key_float (orec, "CRPIX3", crpix);
    kstat += check_and_set_key_float (orec, "CRVAL1", crval);
    kstat += check_and_set_key_float (orec, "CRVAL2", crval);
    kstat += check_and_set_key_float (orec, "CRVAL3", crval);
    kstat += check_and_set_key_str (orec, "CTYPE1", "WAVE-NUM");
    kstat += check_and_set_key_str (orec, "CTYPE2", "WAVE-NUM");
    kstat += check_and_set_key_str (orec, "CTYPE3", "FREQ-ANG");
    kstat += check_and_set_key_float (orec, "apode_k_min", edge_inner);
    kstat += check_and_set_key_float (orec, "apode_k_max", edge_outer);
    kstat += check_and_set_key_float (orec, "apode_f", apode_edge);
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s:\n",
	  out_series);
      fprintf (stderr, "      series may not have appropriate structure;\n");
      fprintf (stderr, "      record %d skipped\n", rgn);
      continue;
    }

    dval = drms_getkey_double (irec, "MapScale", &status);
    if (!status) {
						    /*  assume cols = rows!  */
      dval = 360.0 / dval / cols;				 /*  k (inv radian)  */
      dval /= 696.0;					     /*  k (inv Mm)  */
      kstat += check_and_set_key_double (orec, "CDELT1", dval);
      if (cols == rows) {
	kstat += check_and_set_key_double (orec, "CDELT2", dval);
	kstat += check_and_set_key_double (orec, "Delta_k", dval);
      } else {
	kstat += check_and_set_key_double (orec, "Delta_kx", dval);
	dval *= (double)cols / (double)rows;
	kstat += check_and_set_key_double (orec, "CDELT2", dval);
	kstat += check_and_set_key_double (orec, "Delta_ky", dval);
      }
      if (kstat) {
	fprintf (stderr, "Error writing key value(s) to record in series %s\n",
	  out_series);
        fprintf (stderr, "      series may not have appropriate structure\n");
        fprintf (stderr, "      record %d skipped\n", rgn);
        continue;
      }
    }
    dnu = 1.0e6 / planes / tstep;
    if (fbins > 1) dnu *= fbins;
    kstat += check_and_set_key_double (orec, "Delta_nu", dnu);
    domega = dnu * 2.0e-6 * M_PI;
    kstat += check_and_set_key_double (orec, "D_OMEGA", domega);
    kstat += check_and_set_key_double (orec, "CDELT3", domega);
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s\n",
	  out_series);
      fprintf (stderr, "      series may not have appropriate structure\n");
      fprintf (stderr, "      record %d skipped\n", rgn);
      continue;
    }
		/*  convert complex transform to power spectrum and reorder  */
    normal = 1.0 / ntot;
    normal *= normal;
    powrint = 0.0;
    if (fbins > 1) {
      double nfac = normal / fbins / fbins;
      for (plane = 0, n = 0; plane < hplanes; plane += fbins) {
	for (row = 0; row < rows; row++) {
	  for (col = 0; col < cols; col++, n++) {
	    js = (row == row % hrows) ? row + hrows : row - hrows;
	    is = (col == col % hcols) ? col + hcols : col - hcols;
	    ftrb = ftib = 0.0;
	    for (bin = 0; bin < fbins; bin++) {
		m = (plane + bin) * xcols * rows + js * xcols + 2 * is;
		if (is >= hcols) {
		  m = 2 * (cols - is);    kstat += propagate_keys (orec, irec, propagate, keyct);

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
	    if (n == 0) vmax = vmin = data[n];
	    if (data[n] > vmax) vmax = data[n];
	    if (data[n] < vmin) vmin = data[n];
	    powrint += data[n];
	  }
	}
      }
      if (n > ntot) {
	fprintf (stderr, "Error: data written beyond output array bounds\n");
	return cleanup (1, irecs, orec, orig, pspec, dispose);
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
	    data[n] *= normal;
	    if (n == 0) vmax = vmin = data[n];
	    if (data[n] > vmax) vmax = data[n];
	    if (data[n] < vmin) vmin = data[n];
	    powrint += data[n];
	  }
	}
      }
    }
    if (dp_calc) free (wdata);
    else free (xdata);

    if (verbose) {
      printf ("  P[0] = %11.4e (cf. apodized data mean = %11.4e)\n",
	  data[cols * hrows + hcols], dataavg);
      printf ("  S(P) = %11.4e (cf. apodized data var. = %11.4e)\n", powrint, datavar);
    }
    if ((segs = orec->segments.num_total) < 1) {
      fprintf (stderr, "Error: no data segments in output data series:\n");
      fprintf (stderr, "  %s\n", out_series);
      return cleanup (1, irecs, orec, orig, pspec, dispose);
    }
    for (n = 0; n < segs; n++) {
      oseg = drms_segment_lookupnum (orec, n);
      if (oseg->info->naxis == 3) break;
    }
    if (n >= segs && verbose)
      printf ("found no segmemt of rank 3, using segment %d\n", n);
    switch (oseg->info->type) {
      case DRMS_TYPE_CHAR:
	scale_range = 250.0;
	break;
      case DRMS_TYPE_SHORT:
	scale_range = 65000.0;
	break;
      case DRMS_TYPE_INT:
	scale_range = 4.2e9;
	break;
      case DRMS_TYPE_LONGLONG:
	scale_range = 1.8e19;
	break;
      default:
	scale_range = 1.0;
    }

    if (verbose) printf ("  values range from %.3e to %.3e\n", vmin, vmax);
    if (log_out) {
      double target_min = vmax / exp (scale_range);
      if (vmin < target_min) {
	fprintf (stderr,
	    "** Warning: minimum value = %.2e; value = log (%.2e + data)\n",
	    vmin, target_min);
	for (n = 0; n < ntot; n++) data[n] = log (target_min + data[n]);
	vmin = log (target_min + vmin);
	vmax = log (target_min + vmax);
      } else {
	for (n = 0; n < ntot; n++) data[n] = log (data[n]);
	vmin = log (vmin);
	vmax = log (vmax);
      }
      if (verbose) printf ("  logs scaled between %f and %f\n", vmin, vmax);
      kstat += check_and_set_key_double (orec, "LOG_BASE", exp (1.0));
    }
  
    kstat += check_and_set_key_float (orec, "DataMin", vmin);
    kstat += check_and_set_key_float (orec, "DataMax", vmax);
    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_str (orec, "Source", source);
    kstat += check_and_set_key_str (orec, "Input", input);
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s\n",
	  out_series);
      fprintf (stderr, "      series may not have appropriate structure\n");
      fprintf (stderr, "      record %d skipped\n", rgn);
      continue;
    }
    bscale = (vmax - vmin) / scale_range;
    bzero = 0.5 * (vmax + vmin);
    pspec = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, data, &status);
    if (status) {
      fprintf (stderr, "Error: couldn't create output array\n");
      return cleanup (1, irecs, orec, orig, pspec, dispose);
    }
    pspec->bscale = bscale;
    pspec->bzero = bzero;
    status = drms_segment_write (oseg, pspec, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      return cleanup (1, irecs, orec, orig, pspec, dispose);
    }
    drms_free_array (pspec);
    pspec = NULL;

    if (verbose) {
      drms_sprint_rec_query (recid, orec);
      if (no_save)
	printf ("output data record would be %s\n", recid);
      else {
	drms_segment_filename (oseg, pathname);
	printf ("output data record is %s\n", recid);
	printf ("output data cube at %s\n", pathname);
      }
    }
    drms_close_record (orec, dispose);
  }
                                                      /* Clean up FFTW stuff */
  fftwf_destroy_plan (fplan);
/*  this is kind of silly, why not just remove error argument from cleanup?  */
  return cleanup (0, irecs, NULL, orig, pspec, dispose);
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  08.03.24	created this file, based on SOI powrspec3_v20
 *  v 0.6 frozen 2009.07.14
 *  v 0.7:	fix to run with NetDRMS 2.0 (09.07.14);
 *		Modify and add to list of keywords copied from input
 *  v 0.7 frozen 2009.09.08.05
 *  v 0.8:	general code cleanup; added copying of all standard (and
 *	additional "standard") keys;
 *		removed defaults for I/O;
 *		removed CGI output and special processing for exports;
 *		fixed bug in setting of CDELT3 and associated keywords;
 *		modified key copying to use list and to be more careful;
 *		added setting of Cadence keywords from input CDELT3;
 *		put in (approximately) correct normalization;
 *		generalized scaling of logs (and direct values) for output
 *	precision
 *		added looping over multiple input and output records
 *		added no_save option for diagnostics
 *  v 0.8 frozen 2010.01.18
 *
 */
