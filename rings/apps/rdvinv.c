/*****************************************************************
 *  rdvinv.c						~rick/src/rdvinv
 *
 *  Responsible: Charles Baldner, Rick Bogart, Sarbani Basu
 *							RBogart@spd.aas.org
 *
 *  Velocity inversion code for ring diagram power spectra fits
 *
 *  Parameters: (type   default         description)
 *	in	string	-		Input data set descriptor or filename
 *	out	string	fort.10.hmixy	Output filename
 *	seg	string	fit.out		Input data segment name
 *	cr	int	NS		Carrington rotation common to dataset
 *	clon	float	NS		CM Longitude common to dataset
 *	lon	float	NS		Carrington longitude common to dataset
 *	lat	float	NS		Heliographic latitude common to dataset
 *	amu	double	0.005		Error trade-off parameter
 *	ob	double	1.0		Lower frequency limit (mHz)
 *	oe	double	5.2		Upper frequency limit (mHz)
 *	rb	double	0.97		Lower radius limit
 *	re	double	1.00		Upper radius limit
 *	num	int	40		Number of target inversion points
 *	kernel	string	-		Data record or pathname of file
 *				containing mode kernels for inversions
 *	ave	string	Not Specified	Output file for averaging kernels
 *				(if not set, kernels are not written)
 *	coef	string	Not Specified	Output file for ?? coefficients
 *				(if not set, coefficients are not written)
 *
 *  Flags
 *	-v	run verbose
 *	
 *  Notes:
 *    Output is specified by argment 'out'.  If 'out' is a drms data series
 *	that can be opened, the inversions will be  saved to the appropriate
 *	drms records; otherwise, 'out' will be used as the base filename for
 *	the inversion outputs.
 *
 *  Bugs:
 *    If the input data is not in a drms series, an output series can't be
 *	created - rdvinv will exit with an error message.
 *    Argument 'out', when a file name, cannot contain any path information.
 *    When writing output to a file, rather than a drms record set, only one
 *	set of files can be specified; if the input contains more than one
 *	data record, each inversion will overwrite the last.
 *    There is no verification that the requested output segments have the
 *	correct protocol; probably doesn't matter anyway
 *    The constant value of amu is written needlessly for each target location
 *	for backward compatability
 *    The reading of the kernel file is still in the FORTRAN subroutine
 *
 *  Revision history is at end of file
 ******************************************************************/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jsoc_main.h>
#include "keystuff.c"
#include "old_rdutil.c"
#include "ola_xy.c"

/* prototypes */


extern void ola_(double *, int *, double *, double *, double *,
      double *, double *, int *, char *, double *, double *, double *,
      double *, double *, double *, double *, double *, double *,
      char *, char *, int *, int *, int *, double *, double *,
      int *, double *, double *, double *, int, int, int);
/*
extern void ola_(double *, int *, double *, double *, double *,
      double *, double *, int *, char *, char *, char *, char *, char *,
      int *, int *, int *, double *, double *,
      int *, double *, double *, double *,
      int, int, int, int, int);
*/
char *module_name = "rdvinv";
char *version_id = "0.9";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in", "", "Input data series or recordset"},
  {ARG_STRING, "seg", "fit.out", "Input data segment name"},
  {ARG_STRING, "uxseg", "Ux", "Output data segment name for Ux inversions"},
  {ARG_STRING, "uyseg", "Uy", "Output data segment name for Uy inversions"},
  {ARG_INT,    "cr", "Not Specified", "Carrington rotation for all regions"},
  {ARG_FLOAT,  "clon", "Not Specified",
      "Carrington longitude of CM for all regions"},
  {ARG_FLOAT,  "lon", "Not Specified", "Carrington longitude of all regions"},
  {ARG_FLOAT,  "lat", "Not Specified", "Carrington latitude of all regions"},
  {ARG_DOUBLE, "amu", "0.005", "Error trade-off parameter"},
  {ARG_DOUBLE, "ob", "1.0", "Lower frequency limit (mHz)"},
  {ARG_DOUBLE, "oe", "5.2", "Upper frequency limit (mHz)"},
  {ARG_DOUBLE, "rb", "0.97", "Lower radius limit"},
  {ARG_DOUBLE, "re", "1.00", "Upper radius limit"},
  {ARG_INT,    "num", "40", "Number of target inversion points"},
  {ARG_STRING, "out", "fort.10.hmixy", "Output filename"},
  {ARG_STRING, "kernel", "", ""},
  {ARG_STRING, "ave",  "Not Specified",
      "output file for averaging kernels (if not set, kernels are not written)"},
  {ARG_STRING, "coef", "Not Specified", ""},
  {ARG_FLAG,   "v",    "", "run in verbose mode"},
  {ARG_END}
};

       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "MapProj", "MapScale", "Map_PA", "Width", "Height",
    "Size", "Cadence", "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel", "MAI"};

int process_record (const char *filename, int drms_output, char *outfilex,
    char *outfiley, char *kernel, char *kername, char *ave, char *coef, int qave, int qcoef,
    int verbose, double ob, double oe, int num, double rb, double re,
    double amu, char *sourcerec, char *codever) {
  FILE *infile = fopen (filename, "r");
  FILE *filex = fopen (outfilex, "w");
  FILE *filey = fopen (outfiley, "w");
  double *rcgx, *rcgy, *quartx, *quarty, *solnx, *solny, *errx, *erry;
  double *l, *f, *ef, *ux, *eux, *uy, *euy, *rtrg;
  int *n, *mask;
  int npts, i, j, status;
  int lenkern, lenave, lencoef;
  status = read_fit_v (infile, &npts, &n, &l, &f, &ef, &ux, &eux, &uy, &euy);
  fclose (infile);
  if(status)	{
    fprintf(stderr, "File %s could not be read\n", filename);
    return 1;
  }

  rtrg = (double *)malloc (num * sizeof (double));
  rcgx = (double *)malloc (num * sizeof (double));
  rcgy = (double *)malloc (num * sizeof (double));
  quartx = (double *)malloc (3 * num * sizeof (double));
  quarty = (double *)malloc (3 * num * sizeof (double));
  solnx = (double *)malloc (num * sizeof (double));
  solny = (double *)malloc (num * sizeof (double));
  errx = (double *)malloc (num * sizeof (double));
  erry = (double *)malloc (num * sizeof (double));
  mask = (int *)malloc (npts * sizeof(int));
  interp_vel (n, l, f, ef, ux, eux, uy, euy, npts);
  autoweed_vel (n, l, ux, uy, mask, npts);
  j = 0;
  for (i=0; i<npts; i++) {
    if (mask[i]) {
      n[j] = n[i];
      l[j] = l[i];
      f[j] = f[i];
      ef[j] = ef[i];
      ux[j] = ux[i];
      eux[j] = eux[i];
      uy[j] = uy[i];
      euy[j] = euy[i];
      j++;
    }
  }
  npts = j;

  fprintf (filex, "# Flow inversion of %s\n", sourcerec);
  fprintf (filey, "# Flow inversion of %s\n", sourcerec);
  fprintf (filex, "# against kernel %s\n", kername);
  fprintf (filey, "# against kernel %s\n", kername);
  fprintf (filex, "# %s\n", codever);
  fprintf (filey, "# %s\n", codever);
  fprintf (filex, "# amu = %.5f ob = %.3f mHz, oe = %.3f mHz\n", amu, ob, oe);
  fprintf (filey, "# amu = %.5f ob = %.3f mHz, oe = %.3f mHz\n", amu, ob, oe);
  fprintf (filex, "# rb = %.3f, re= %.3f, num = %d\n",rb, re, num);
  fprintf (filey, "# rb = %.3f, re= %.3f, num = %d\n",rb, re, num);
  fprintf (filex,
      "# Rtrg    amu    CG of av. kernel, quartiles                soln          err\n");
  fprintf (filex, "#                           1       2       3     3-1\n");
  fprintf (filey,
      "# Rtrg    amu    CG of av. kernel, quartiles                soln          err\n");
  fprintf (filey, "#                           1       2       3     3-1\n");

  lenkern = strlen (kernel);
  lenave = (qave) ? strlen (ave) : 0;
  lencoef = (qcoef) ? strlen (coef) : 0;

  ola_ (l, n, f, ux, eux, uy, euy, &npts, kernel, rtrg, rcgx, rcgy,
      quartx, quarty, solnx, solny, errx, erry, ave, coef, &qave, &qcoef,
      &verbose, &ob, &oe, &num, &rb, &re, &amu, lenkern, lenave, lencoef);
/*
  ola_xy (l, n, f, ux, eux, uy, euy, npts, kernel, outfilex, outfiley, ave,
      coef, qave, qcoef, verbose, ob, oe, num, rb, re, amu);
*/
  for (i = 0; i < num; i++)
    fprintf (filex, "%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %12.5e %12.5e\n",
	rtrg[i], amu, rcgx[i], quartx[2 + 3*i], quartx[1 + 3*i], quartx[3*i],
	fabs (quartx[3*i] - quartx[2+ 3*i]), solnx[i], errx[i]);
  for (i = 0; i < num; i++)
    fprintf (filey, "%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %12.5e %12.5e\n",
	rtrg[i], amu, rcgy[i], quarty[2 + 3*i], quarty[1 + 3*i], quarty[3*i],
	fabs (quarty[3*i] - quarty[2+ 3*i]), solny[i], erry[i]);

  free (n);
  free (mask);
  free (l);
  free (f);
  free (ef);
  free (ux);
  free (eux);
  free (uy);
  free (euy);

  free (rtrg);
  free (rcgx);
  free (rcgy);
  free (quartx);
  free (quarty);
  free (solnx);
  free (solny);
  free (errx);
  free (erry);
  fclose (filex);
  fclose (filey);

  return 0;
}
		 /*  these functions should probably be added to keystuff.c  */
					    /*  they are also in rdsinv_v07  */
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

int unique_key_values (const char *dsqry, const char *key, DRMS_Record_t *rec) {
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

char *printcrcl_loc (int cr, double cl, double lnhg, double lncm, double lat) {
  static char fnam[32];
  double dlon = lncm;

  if (cr > 0) {
    if (isnan (cl)) {
      if (isnan (lncm)) {
				       /*  construct a filename CR_lonHG+lat  */
	sprintf (fnam, "%04d_%04.1f%+05.1f", cr, lnhg, lat);
      } else {
				      /*  construct a filename CR_lonEWlatNS  */
	if (dlon >= 0) {
	  if (lat >= 0) sprintf (fnam, "%04d_%04.1fW%04.1fN", cr, dlon, lat);
	  else if (lat > -0.05) sprintf (fnam, "%04d_%04.1fW%04.1fN",
	      cr, dlon, -lat);
          else sprintf (fnam, "%04d_%04.1fW%04.1fS", cr, dlon, -lat);
	} else if (dlon > -0.05) {
	  if (lat >= 0) sprintf (fnam, "%04d_%04.1fW%04.1fN", cr, -dlon, lat);
	  else if (lat > -0.05) sprintf (fnam, "%04d_%04.1fW%04.1fN",
	      cr, -dlon, -lat);
          else sprintf (fnam, "%04d_%04.1fW%04.1fS", cr, -dlon, -lat);
	} else {
	  if (lat >= 0) sprintf (fnam, "%04d_%04.1fE%04.1fN", cr, -dlon, lat);
	  else if (lat > -0.05) sprintf (fnam, "%04d_%04.1fE%04.1fN",
	      cr, -dlon, -lat);
          else sprintf (fnam, "%04d_%04.1fE%04.1fS", cr, -dlon, -lat);
	}
      }
      return fnam;
    }
				   /*  construct a filename CR:CL_lonEWlatNS  */
    while (cl < 0) {
      cl += 360;
      cr++;
    } while (cl > 360) {
      cl -= 360;
      cr--;
    }
    if (isnan (lncm)) {
      dlon = lnhg - cl;
      while (dlon > 180) dlon -= 360;
      while (dlon < -180) dlon += 360;
    }
    if (dlon >= 0) {
      if (lat >= 0)
        sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fN", cr, cl, dlon, lat);
      else if (lat > -0.05)
        sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fN", cr, cl, dlon, -lat);
      else sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fS", cr, cl, dlon, -lat);
    } else if (dlon > -0.05) {
      if (lat >= 0)
        sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fN", cr, cl, -dlon, lat);
      else if (lat > -0.05)
        sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fN", cr, cl, -dlon, -lat);
      else sprintf (fnam, "%04d:%05.1f_%04.1fW%04.1fS", cr, cl, -dlon, -lat);
    } else {
      if (lat >= 0)
        sprintf (fnam, "%04d:%05.1f_%04.1fE%04.1fN", cr, cl, -dlon, lat);
      else if (lat > -0.05)
        sprintf (fnam, "%04d:%05.1f_%04.1fE%04.1fN", cr, cl, -dlon, -lat);
      else sprintf (fnam, "%04d:%05.1f_%04.1fE%04.1fS", cr, cl, -dlon, -lat);
    }
  } else {
    if (isnan (dlon)) {
					  /*  construct a filename lonHG+lat  */
      sprintf (fnam, "%04.1f%+05.1f", lnhg, lat);
      return fnam;
    }
					 /*  construct a filename lonEWlatNS  */
    if (dlon >= 0) {
      if (lat >= 0) sprintf (fnam, "%04.1fW%04.1fN", dlon, lat);
      else if (lat > -0.05) sprintf (fnam, "%04.1fW%04.1fN", dlon, -lat);
      else sprintf (fnam, "%04.1fW%04.1fS", dlon, -lat);
    } else if (dlon > -0.05) {
      if (lat >= 0) sprintf (fnam, "%04.1fW%04.1fN", -dlon, lat);
      else if (lat > -0.05) sprintf (fnam, "%04.1fW%04.1fN", -dlon, -lat);
      else sprintf (fnam, "%04.1fW%04.1fS", -dlon, -lat);
    } else {
      if (lat >= 0) sprintf (fnam, "%04.1fE%04.1fN", -dlon, lat);
      else if (lat > -0.05) sprintf (fnam, "%04.1fE%04.1fN", -dlon, -lat);
      else sprintf (fnam, "%04.1fE%04.1fS", -dlon, -lat);
    }
  }
  return fnam;
}
							     /*  module body  */
int DoIt(void)	{
  CmdParams_t *params = &cmdparams;
  int status = 0;
  int drms_input, drms_output;
  int nrec, rec_i;
  DRMS_RecordSet_t *recordSet = NULL, *kern_set = NULL;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *segment, *oseg;
  FILE *fpt;
  double latc, loncm, lonhg, loncCar;
  float fval;
  int ival;
  int n;
  char odir[DRMS_MAXPATHLEN], filename[DRMS_MAXPATHLEN+5];
  char buffer[1024], outfilex[DRMS_MAXPATHLEN], outfiley[DRMS_MAXPATHLEN];
  char outdir[DRMS_MAXPATHLEN], outdirx[DRMS_MAXPATHLEN],
      outdiry[DRMS_MAXPATHLEN], kernfile[DRMS_MAXPATHLEN];
  char *carstr;
  char rec_query[DRMS_MAXQUERYLEN], source[DRMS_MAXQUERYLEN];
  char module_ident[64];

  int twosegs = 0;

  char *in = strdup (params_get_str (params, "in"));
  int cr = params_get_int (params, "cr");
  float cl = params_get_float (params, "clon");
  float lat = params_get_float (params, "lat");
  float lon = params_get_float (params, "lon");
  char *seg = strdup (params_get_str  (params, "seg"));
  char *uxseg = strdup (params_get_str  (params, "uxseg"));
  char *uyseg = strdup (params_get_str  (params, "uyseg"));
  char *kernel = strdup (params_get_str (params, "kernel"));
  char *out = strdup (params_get_str (params, "out"));
  char *ave = strdup (params_get_str (params, "ave"));
  char *coef = strdup (params_get_str (params, "coef"));
  double ob = params_get_double (params, "ob");
  double oe = params_get_double (params, "oe");
  double rb = params_get_double (params, "rb");
  double re = params_get_double (params, "re");
  double amu = params_get_double (params, "amu");
  int num = params_get_int (params, "num");
  int verbose = params_isflagset (params, "v");

  int qave = strcmp (ave, "Not Specified");
  int qcoef = strcmp (coef, "Not Specified");
  int anycr = (cr <= 0);
  int anycl = isnan (cl);
  int anylon = isnan (lon);
  int anylat = isnan (lat);

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);

  drms_input = 1;
			 /*  attempt to open the kernel as a DRMS record set  */
  kern_set = drms_open_records (drms_env, kernel, &status);
  if (status) {
			       /*  attempt to open the kernel as a file name  */
    fpt = fopen (kernel, "r");
    if (!fpt) {
      fprintf (stderr, "Error: kernel specification \"%s\"\n", kernel);
      fprintf (stderr,
          "       does not specify a readable data record nor file\n");
      return 1;
    }
    fclose (fpt);
    strcpy (kernfile, kernel);
  } else {
		       /*  check that the kernel set specifies unique record  */
    if (kern_set->n != 1) {
      fprintf (stderr, "Error: no data records in set %s\n", kernel);
      return 1;
    }
    irec = kern_set->records[0];
    drms_record_directory (irec, kernfile, 1);
    segment = drms_segment_lookup (irec, "kernel");
    if (!segment) {
			       /*  and that it contains the segment "kernel"  */
      fprintf (stderr, "Error: kernel record \"%s\"\n", kernel);
      fprintf (stderr, "       does not contain a segment \"kernel\"\n");
      drms_close_records (kern_set, DRMS_FREE_RECORD);
      return 1;
    }
    sprintf (kernfile, "%s/S%05i/kernel", segment->record->su->sudir,
	segment->record->slotnum);
    drms_close_records (kern_set, DRMS_FREE_RECORD);
  }

  if (anycr && anycl && anylon && anylat) {
		 /*  no target range specified, assume implicit in input set  */
			 /*  but it might be implicit in output series spec!  */
    recordSet = drms_open_records (drms_env, in, &status);
    if (status) {
      fpt = fopen (in, "r");
      if (!fpt) {
	fprintf (stderr, "Error: in specification \"%s\"\n", in);
	fprintf (stderr,
	    "       does not specify a readable data set nor file\n");
	return 1;
      }
      fclose (fpt);
      drms_input = 0;
      nrec = 1;
    } else nrec = recordSet->n;
    if (nrec < 1) {
      fprintf (stderr, "Error: no data records in set %s\n", in);
      return 1;
    }
  } else {
    strncpy (rec_query, in, DRMS_MAXQUERYLEN);
    if (!anycr) {
      snprintf (source, DRMS_MAXQUERYLEN, "[?CarrRot=%d?]", cr);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
    }
    if (!anycl) {
      float clmin = cl - 0.01;
      float clmax = cl + 0.01;
						/*  need to deal with 0/360  */
      if (clmin < 0.0) {
        clmin += 360.0;
	snprintf (source, DRMS_MAXQUERYLEN, "[?CMLon>%g or CMLon<%g?]",
	    clmin, clmax);
      } else if (clmax > 360.0) {
        clmax -= 360.0;
	snprintf (source, DRMS_MAXQUERYLEN, "[?CMLon>%g or CMLon<%g?]",
	    clmin, clmax);
      } else snprintf (source, DRMS_MAXQUERYLEN, "[?CMLon>%g and CMLon<%g?]",
	  clmin, clmax);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
    }
    if (!anylon) {
      float lonmin = lon - 0.01;
      float lonmax = lon + 0.01;
						/*  need to deal with 0/360  */
      if (lonmin < 0.0) {
        lonmin += 360.0;
	snprintf (source, DRMS_MAXQUERYLEN, "[?LonHG>%g or LonHG<%g?]",
	    lonmin, lonmax);
      } else if (lonmax > 360.0) {
        lonmax -= 360.0;
	snprintf (source, DRMS_MAXQUERYLEN, "[?LonHG>%g or LonHG<%g?]",
	    lonmin, lonmax);
      } else snprintf (source, DRMS_MAXQUERYLEN, "[?LonHG>%g and LonHG<%g?]",
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
  if (nrec < 1) {
    fprintf (stderr, "Error: no records found in input data set\n");
    fprintf (stderr, "       %s\n", rec_query);
  }
  if (drms_input) {
		/*  check that range is restricted by dataset specification  */
    irec = recordSet->records[0];
    if (anycr) {
      anycr = unique_key_values (in, "CarrRot", irec) - 1;
      ival = drms_getkey_int (irec, "CarrRot", &status);
      if (!status) cr = ival;
    }
    if (anycl) {
      anycl = unique_key_values (in, "CMLon", irec) - 1;
      fval = drms_getkey_float (irec, "CMLon", &status);
      if (!status) cl = fval;
    }
    if (anylon) {
      anylon = unique_key_values (in, "LonHG", irec) - 1;
      fval = drms_getkey_float (irec, "LonHG", &status);
      if (!status) lon = fval;
    }
    if (anylat) {
      anylat = unique_key_values (in, "LatHG", irec) - 1;
      fval = drms_getkey_float (irec, "LatHG", &status);
      if (!status) lat = fval;
    }
    if (verbose) printf ("Processing %d records\n", nrec);
  }
				   /*  Attempt to create output drms record  */
  drms_output = 1;
  orec = drms_create_record (drms_env, out, DRMS_PERMANENT, &status);
  if (status) {
    drms_output = 0;
    fprintf (stderr,
        "Warning: drms_create_record() returned %d for data series:\n", status);
    fprintf (stderr,
	"       %s\n       will be interpreted as directory name\n", out);
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

    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
					       /*  module specific keywords  */
    kstat += check_and_set_key_float (orec, "amu", amu);
    kstat += check_and_set_key_float (orec, "freqmin", ob);
    kstat += check_and_set_key_float (orec, "freqmax", oe);
    kstat += check_and_set_key_float (orec, "radmin", rb);
    kstat += check_and_set_key_float (orec, "radmax", re);
    kstat += check_and_set_key_str (orec, "BLD_VERS", jsoc_version);
				 /*  set possible prime keys as appropriate  */
    if (anycr) {
      if (key_is_prime (orec, "CarrRot")) {
	fprintf (stderr,
	    "Error: CarrRot is prime key for output series %s\n", out);
	fprintf (stderr,
	    "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
				     /*  this pretty much has to be true...  */
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_int (orec, "CarrRot", cr);
    if (anycl) {
      if (key_is_prime (orec, "CMLon")) {
	fprintf (stderr,
	    "Error: CMLon is prime key for output series %s\n", out);
	fprintf (stderr,
	    "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_double (orec, "CMLon", cl);
    if (anycr || anycl) {
      if (key_is_prime (orec, "CarrTime")) {
	fprintf (stderr,
	    "Error: CarrTime is prime key for output series %s\n", out);
	fprintf (stderr,
	    "       but input data set contains multiple values\n");
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
	fprintf (stderr,
	    "Error: LonHG is prime key for output series %s\n", out);
	fprintf (stderr,
	    "       but input data set contains multiple values\n");
	drms_close_record (orec, DRMS_FREE_RECORD);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	return 1;
      }
    } else kstat += check_and_set_key_double (orec, "LonHG", lon);
    if (anylat) {
      if (key_is_prime (orec, "LatHG")) {
	fprintf (stderr,
	    "Error: LatHG is prime key for output series %s\n", out);
	fprintf (stderr,
	    "       but input data set contains multiple values\n");
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
			   /*  search for separate output segments ux and uy  */
    if (drms_segment_lookup (orec, uxseg) && drms_segment_lookup (orec, uyseg))
      twosegs = 1;
    else {
			       /*  search for 1st generic segment and use it  */
      int segct = drms_record_numsegments (orec);
      if (segct < 1) {
	fprintf (stderr,
	    "Error: no data segment in output series %s\n", out);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	drms_close_record (orec, DRMS_FREE_RECORD);
	return 1;
      }
      for (n = 0; n < segct; n++) {
	oseg = drms_segment_lookupnum (orec, n);
	if (oseg->info->protocol != DRMS_GENERIC) continue;
   	break;
      }
      if (n == segct) {
	fprintf (stderr,
	    "Error: no generic data segment in output series %s\n", out);
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	drms_close_record (orec, DRMS_FREE_RECORD);
	return 1;
      }
      fprintf (stderr, "         writing to segment %s\n", oseg->info->name);
    }
  }
  if (twosegs) {
    sprintf (outdirx, "%s/%s", odir, uxseg);
    sprintf (outdiry, "%s/%s", odir, uyseg);
    mkdir (outdirx, 01755);
    mkdir (outdiry, 01755);
  } else {
    if (drms_output) {
      sprintf (outdir, "%s/%s", odir, oseg->info->name);
      mkdir (outdir, 01755);
    } else strcpy (outdir, out);
  }
						     /*  main processing loop  */
  for (rec_i=0; rec_i<nrec; rec_i++) {
    if (verbose) printf ("  Processing record %i\n", rec_i);
    if (drms_input) {
      irec = recordSet->records[rec_i];
		       /*  should use call to drms_record_directory() instead  */
      if (irec->sunum != -1LL && irec->su == NULL) {
	irec->su = drms_getunit (irec->env, irec->seriesinfo->seriesname,
	    irec->sunum, 1, &status);
      }
      segment = drms_segment_lookup(irec, seg);
      sprintf (filename, "%s/S%05i/%s", segment->record->su->sudir,
	  segment->record->slotnum, seg);
      drms_sprint_rec_query (source, irec);
      cr = drms_getkey_int (irec, "CarrRot", &status);
      latc = drms_getkey_double (irec, "LatHG", &status);
      lonhg = drms_getkey_double (irec, "LonHG", &status);
      carstr = drms_getkey_string (irec, "CarrTime", &status);
      if (status) loncCar = drms_getkey_double (irec, "CMLon", &status);
      else sscanf (carstr, "%i:%lf", &cr, &loncCar);
      loncm = drms_getkey_double (irec, "LonCM", &status);
    } else {
      strcpy (filename, in);
      strcpy (source, in);
      latc = loncm = lonhg = 0.0;
    }
    if (twosegs) {
      sprintf (outfilex, "%s/%s/%s.Ux", odir, uxseg, 
	  printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
      sprintf (outfiley, "%s/%s/%s.Uy", odir, uyseg,
	  printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
    } else {
      sprintf (outfilex, "%s/%s.Ux", outdir, 
	  printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
      sprintf (outfiley, "%s/%s.Uy", outdir,
	  printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
    }
    status = process_record (filename, drms_output, outfilex, outfiley,
	kernfile, kernel, ave, coef, qave, qcoef, verbose, ob, oe, num, rb, re, amu,
	source, module_ident);
    if (status && drms_output)	{
      fprintf (stderr, "Error processing record %d (%s); aborted\n", rec_i,
	  printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
      if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
      if (drms_output) drms_close_record (orec, DRMS_FREE_RECORD);
      return 1;
    }
  }

  if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
  if (drms_output) drms_close_record (orec, DRMS_INSERT_RECORD);

  return 0;
}

/*
 *  Revision History
 *
 *  June 2009 - v 0.1 module created, no real drms functionality
 *  09.10.2009 - v 0.2 drms functionality added - input and output from drms
 *	records now  possible, but input/output to rooted filenames also
 *	possible
 *  0.3 2009.10.22	fixed format of name of output files
 *  0.4 2010.02.11	unified processing calls for named file and drms input
 *	branches
 *  0.5 2010.02.12	Combined branches into one, added some very obvious
 *	error handling (no given input argument, input cannot be opened)
 *			removed rdutil.h inclusion
 *  v 0.5 frozen 2010.03.08
 *    0.6	Added logic for specifying input grouping of data set records
 *		by Carr Rot, CM Longitude, Longitude, and/or Latitude
 *		Added keyword propagation list (but not propagation!)
 *		Output to multiple files named by lat and lon in segment
 *		directory
 *		Fixed processing of longitude queries near 0/360
 *	2010.03.10	Fixed setting of keys for run parameters; put in trap
 *		for empty input data sets and widened clon match range
 *  v 0.6 frozen 2010.04.23
 *    0.7		Output data files are now written into per segment
 *		subdirectories; added option for naming individual segments
 *		(if there are two)
 *  v 0.7 frozen 2010.08.19
 *    0.8	Created 2011.04.22; Fixed typos in ola_xy.c (not yet used)
 *		Fixed bugs in printing of filenames (particularly the
 *		reversal of east and west), and dropped requirement for CMLon;
 *		heliographic option for file naming
 *		Removed default for kernel
 *    0.9	Created 2011.12.04; Fixed typos in ola_xy.c (not yet used)
 *  v 0.8 frozen 2012.02.01 (but included in JSOC release 6.1 2011.12.05 in
 *		identical form except for comments)
 *    0.9	added (preferred) option for specifying kernel as DRMS data
 *		record rather than filename
 */
