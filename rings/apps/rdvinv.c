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
 *	out	string	-		Output data series or directory
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
 *	nmax	int	14		Maximum acceptable order
 *	lmax	int	2000		Maximum acceptable order
 *	kernel	string	-		Data record or pathname of file
 *				containing mode kernels for inversions
 *	ave	string	Not Specified	Output file for averaging kernels
 *				(if not set, kernels are not written)
 *	coef	string	Not Specified	Output file for inversion coefficients
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
 *    nmax is set unconditionally to 8 to reflect hard-coded value in
 *	autoweed_vel
 *    If cr, clon, lon, and/or lat is specified, the input data must be in a
 *	drms series with the appropriate keywords CarrRot, CMLon, LonHG, and
 *	LatHG.
 *    It is not possible to restrict the DRMS input record set based on
 *	Stonyhurst longitude via a calling parameter, only via the recordset
 *	descriptor.
 *    There is no verification that the requested output segments have the
 *	correct protocol; probably doesn't matter anyway
 *    The constant value of amu is written needlessly for each target location
 *	for backward compatability
 *    The mode weeding for abnormal values of Ui is only applied to modes with
 *	n < 8; all higher-order modes are rejected; there are other hard-coded
 *	parameters in the autoweed function
 *    The minimum acceptable error estimate for Ui is hardcoded via a define
 *    Location specifications for the optional files with averaging kernels and
 *	or inversion coefficients are for named directories only; neither SUMS
 *	inclusion nor DRMS segment specification is supported
 *    It would be better to specify the frequency bounds in uHz rather than mHz,
 *	and as fmin and fmax rather than ob and oe
 *
 *  Revision history is at end of file
 ******************************************************************************/


#define MODULE_VERSION_NUMBER	("0.93")
#define KEYSTUFF_VERSION "keystuff_v10.c"
#define RDUTIL_VERSION "rdutil_v09.c"
#define OLAXY_VERSION "olaxy_v13.c"

char *module_name = "rdvinv";
char *version_id = MODULE_VERSION_NUMBER;

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jsoc_main.h>
#include KEYSTUFF_VERSION
#include RDUTIL_VERSION
#include OLAXY_VERSION

#define MIN_ERR	(1.0e-11)

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in", "", "Input data series or recordset"},
  {ARG_STRING, "seg", "fit.out", "Input data segment name"},
  {ARG_STRING, "out", "", "Output series or directory"},
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
  {ARG_INT,    "nmax", "14", "Maximum acceptable mode order"},
  {ARG_INT,    "lmax", "2000", "Maximum acceptable mode degree"},
  {ARG_STRING, "kernel", "", ""},
  {ARG_STRING, "ave",  "Not Specified",
      "output directory for averaging kernels (if not set, kernels are not written)"},
  {ARG_STRING, "coef", "Not Specified",
      "output directory for inversion coefficients (if not set, coefficients are not written)"},
  {ARG_FLAG,   "v",    "", "run in verbose mode"},
  {ARG_END}
};

	/*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "MapProj", "MapScale", "Map_PA", "Width", "Height",
    "Size", "Cadence", "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel", "MAI"};

int process_record (const char *filename, char *outfilex, char *outfiley,
    double *krnrval, double *weight, int nrpts, int *kmoden, int *kmodel,
    double *krnamp, int krnmodect, char *kername, int nmax, int lmax,
    int *kernlloc, int *kernlnct,
    char *avgkrnx, char *avgkrny, int wantavgkrn,
    char *invcoefx, char *invcoefy, int wantinvcoef,
    int *wmoden, int *wmodel, double ob, double oe,
    int numtrg, double rb, double re, double mu, char *sourcerec, char *codever) {
  FILE *infile = fopen (filename, "r");
  FILE *filex = fopen (outfilex, "w");
  FILE *filey = fopen (outfiley, "w");
  double *l, *f, *ef, *ux, *eux, *uy, *euy;
/*
  double fmin, fmax;
  static int *kmode = NULL, *newn = NULL, *newl = NULL;
*/
  int *moden, *mask;
  int i, j, k, krnmode, nmodes, numker;
  int allfits, weededfits, target, status;
		/*  read in all "mode" fit values of n, "l" (fractional),
						  f [uHz], Ui and D_Ui [m/s]  */
  status = read_fit_v (infile, &allfits, &moden, &l, &f, &ef,
      &ux, &eux, &uy, &euy);
  fclose (infile);
  if (status)	{
    fprintf (stderr, "File %s could not be read\n", filename);
    return 1;
  }

  mask = (int *)malloc (allfits * sizeof(int));
  interp_vel (moden, l, f, ef, ux, eux, uy, euy, allfits);
  autoweed_vel (moden, l, ux, uy, mask, allfits);
				 /*  and only use the fits surviving weeding  */
  j = 0;
  for (i = 0; i < allfits; i++) {
    if (mask[i]) {
      moden[j] = moden[i];
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
  weededfits = j;
/*
 *  The following code is intended to replace the search through the
 *    entire kernel file for each mode in ola_xy(); however, it has
 *    not been implemented; if it were, it would be formally necessary
 *    to fix the bug assuming that the kernel file contains modes starting
 *    with order 0 for each degree, which is not true for any of the kernel
 *    files available for degrees 0 and 1
 */
/*
  kmode = (int *)realloc (kmode, weededfits * sizeof (int));
  newn = (int *)realloc (newn, weededfits * sizeof (int));
  newl = (int *)realloc (newl, weededfits * sizeof (int));
  for (i = 0; i < weededfits; i++) {
    newn[i] = newl[i] = -1;
  }
  fmin = 1.0e3 * ob;
  fmax = 1.0e3 * oe;
*/
				      /*  read in observed values and errors,
					 get corresponding kernel amplitudes  */
/*
  for (i = 0, nmodes = 0; i < weededfits; i++) {
    int ni = n[i];
    int li = l[i];
    if (li >= lmax || ni >= nmax) continue;
    if (f[i] <= fmin || f[i] >= fmax) continue;
    if (eux[i] <= MIN_ERR || euy[i] <= MIN_ERR) continue;
    krnmode = kernlloc[li];
    for (j = 0; j < kernlnct[li]; j++) if (kmoden[j + krnmode] == ni) break;
    krnmode += j;
*/
		  /*  check that this is a mode for which there is a kernel!  */
/*
    if (ni != kmoden[krnmode] || li != kmodel[krnmode]) continue;
    newn[nmodes] = ni;
    newl[nmodes] = li;
    kmode[nmodes] = krnmode;
    nmodes++;
  }
  numker = nmodes + 1;
*/
		    /*  only use fits for which there are acceptable kernels  */
/*
  j = 0;
  for (i = 0, nmodes = 0; i < weededfits; i++) {
    if (n[i] == newn[nmodes] && l[i] == newl[nmodes]) {
      n[j] = n[i];
      l[j] = l[i];
      f[j] = f[i];
      ef[j] = ef[i];
      ux[j] = ux[i];
      eux[j] = eux[i];
      uy[j] = uy[i];
      euy[j] = euy[i];
      kmode[j] = kmode[i];
      nmodes++;
      j++;
    }
  }
*/
  status = ola_xy (l, moden, f, ux, eux, uy, euy, weededfits,
      krnrval, wmoden, wmodel, weight, nrpts, nmax, lmax,
      kmoden, kmodel, krnamp, krnmodect, &nmodes, ob, oe, numtrg, mu);
  if (status) return status;
							   /*  print headers  */
  fprintf (filex, "# Flow inversion of %s\n", sourcerec);
  fprintf (filey, "# Flow inversion of %s\n", sourcerec);
  fprintf (filex, "# against kernel %s\n", kername);
  fprintf (filey, "# against kernel %s\n", kername);
  fprintf (filex, "# %s\n", codever);
  fprintf (filey, "# %s\n", codever);
  fprintf (filex, "# error trade-off parameter mu = %.5f\n", mu);
  fprintf (filey, "# error trade-off parameter mu = %.5f\n", mu);
  fprintf (filex, "# mode frequency limits ob = %.3f mHz, oe = %.3f mHz\n", ob, oe);
  fprintf (filey, "# mode frequency limits ob = %.3f mHz, oe = %.3f mHz\n", ob, oe);
  fprintf (filex, "# mode limits nmax = %d, lmax = %d\n", nmax, lmax);
  fprintf (filey, "# mode limits nmax = %d, lmax = %d\n", nmax, lmax);
  fprintf (filex, "# target bounds rb = %.3f, re= %.3f, num = %d\n", rb, re, numtrg);
  fprintf (filey, "# target bounds rb = %.3f, re= %.3f, num = %d\n", rb, re, numtrg);
  fprintf (filex,
      "# Rtrg     mu    CG of av. kernel, quartiles                soln          err\n");
  fprintf (filex, "#                           1       2       3     3-1\n");
  fprintf (filey,
      "# Rtrg     mu    CG of av. kernel, quartiles                soln          err\n");
  fprintf (filey, "#                           1       2       3     3-1\n");
  for (i = 0; i < numtrg; i++)
    fprintf (filex, "%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %12.5e %12.5e\n",
	rtrg[i], mu, rcgx[i], quartx[2 + 3*i], quartx[1 + 3*i], quartx[3*i],
	fabs (quartx[3*i] - quartx[2 + 3*i]), solnx[i], errx[i]);
  for (i = 0; i < numtrg; i++)
    fprintf (filey, "%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %12.5e %12.5e\n",
	rtrg[i], mu, rcgy[i], quarty[2 + 3*i], quarty[1 + 3*i], quarty[3*i],
	fabs (quarty[3*i] - quarty[2 + 3*i]), solny[i], erry[i]);
  fclose (filex);
  fclose (filey);

  if (wantavgkrn) {
    filex = fopen (avgkrnx, "w");
    filey = fopen (avgkrny, "w");
    fprintf (filex, "# radius: ");
    fprintf (filey, "# radius: ");
    for (i = 0; i < numtrg; i++) fprintf (filex, "%13.6f", rtrg[i]);
    fprintf (filex, "\n");
    for (i = 0; i < numtrg; i++) fprintf (filey, "%13.6f", rtrg[i]);
    fprintf (filey, "\n");
    for (j = 0; j < nrpts; j++) {
      fprintf (filex, "%12.6e", krnrval[j]);
      for (k = 0; k < numtrg; k++) fprintf (filex, "%13.5e", avccx[numtrg*j + k]);
      fprintf (filex, "\n");
      fprintf (filey, "%12.6e", krnrval[j]);
      for (k = 0; k < numtrg; k++) fprintf (filey, "%13.5e", avccy[numtrg*j + k]);
      fprintf (filey, "\n");
    }
    fclose (filex);
    fclose (filey);
  }
  if (wantinvcoef) {
    filex = fopen (invcoefx, "w");
    filey = fopen (invcoefy, "w");
    fprintf (filex, "#  l  n    U_x       R: %8.5f", rtrg[0]);
    fprintf (filey, "#  l  n    U_y       R: %8.5f", rtrg[0]);
    for (i = 1; i < numtrg; i++) fprintf (filex, "%13.5f", rtrg[i]);
    fprintf (filex, "\n");
    for (i = 1; i < numtrg; i++) fprintf (filey, "%13.5f", rtrg[i]);
    fprintf (filey, "\n");
    for (j = 0; j < nmodes; j++) {
      fprintf (filex, "%4d %2d %11.4e", wmodel[j], wmoden[j], domx[j]);
      for (i = 0; i < numtrg; i++) fprintf (filex, "%13.5e", ctil2x[j + i * nrpts]);
      fprintf (filex, "\n");
      fprintf (filey, "%4d %2d %11.4e", wmodel[j], wmoden[j], domy[j]);
      for (i = 0; i < numtrg; i++) fprintf (filey, "%13.5e", ctil2y[j + i * nrpts]);
      fprintf (filey, "\n");
    }
    fclose (filex);
    fclose (filey);
  }

  free (moden);
  free (mask);
  free (l);
  free (f);
  free (ef);
  free (ux);
  free (eux);
  free (uy);
  free (euy);
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
						 /*  read in the kernel file  */
	      /*  this is somewhat time-consuming, so delay till all is okay  */
int read_kernel (char *filename, double **krnrval, int **moden, int **model,
    double **modeamp, int *kernrvals, int **krnlloc, int **krnlnct,
    int minl, int maxl) {
/*
 *  It is assumed that the modes appear in order of increasing degree, and
 *    within those for a given degree, in order of increasing mode order.
 *    This is the case for the three kernels currently in DRMS series
 *    hmi.HSKernels[BS05,cm_opal78,invbb77].
 *  krnlloc is a list of the locations within the kernel array where modes of
 *    a given degree start, of length maxl+1 (regardless of minl). krnlnct
 *    is a list of the number of modes present for the given degree.
 *  Note: it would be better if the kernel record were in DRMS with metadata
 *    values for at least lmin, lmax, and modect
 */
  double modlfreq;
  float model_radius, model_mass;
  long startmodes;
  int modldeg, modlord, n;

  int lastdeg = -1, modect = 0, maxmodect = 1000;

  FILE *kernfile = fopen (filename, "r");
  if (!kernfile) return -1;
  *krnlloc = (int *)malloc ((maxl + 1) * sizeof (int));
  *krnlnct = (int *)calloc ((maxl + 1), sizeof (int));
  fscanf (kernfile, "%g %g %d", &model_radius, &model_mass, kernrvals);
  *krnrval = (double *)malloc (*kernrvals * sizeof (double));
  for (n = 0; n < *kernrvals; n++) fscanf (kernfile, "%lg", &((*krnrval)[n]));
  *modeamp = (double *)malloc (*kernrvals * maxmodect * sizeof (double));
  *moden = (int *)malloc (maxmodect * sizeof (int));
  *model = (int *)malloc (maxmodect * sizeof (int));
  modect = 0;
  while (fscanf (kernfile, "%d %d %lg", &modldeg, &modlord, &modlfreq) == 3) {
    for (n = 0; n < *kernrvals; n++) fscanf (kernfile, "%lg",
	&((*modeamp)[n + *kernrvals * modect]));
    (*moden)[modect] = modlord;
    (*model)[modect] = modldeg;
    if (modldeg < minl) continue;
    if (modldeg > maxl) {
      fclose (kernfile);
      return modect;
    }
    if (modldeg > lastdeg) {
      lastdeg = modldeg;
      (*krnlloc)[modldeg] = modect;
    }
    (*krnlnct)[modldeg]++;
    modect++;
    if (modect >= maxmodect) {
      maxmodect += 1000;
      *modeamp = (double *)realloc (*modeamp, *kernrvals * maxmodect * sizeof (double));
      *moden = (int *)realloc (*moden, maxmodect * sizeof (int));
      *model = (int *)realloc (*model, maxmodect * sizeof (int));
    }
  }
  fclose (kernfile);
  return modect;
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
  double *krnrval, *krnamp, *weight;
  float fval;
  int *kernlloc, *kernlnct, *moden, *model, *wmoden, *wmodel;
  int ival;
  int kernrvals, modect, n;
  char odir[DRMS_MAXPATHLEN], filename[DRMS_MAXPATHLEN+5];
  char outfilex[DRMS_MAXPATHLEN], outfiley[DRMS_MAXPATHLEN];
  char avgkrnfilex[DRMS_MAXPATHLEN], avgkrnfiley[DRMS_MAXPATHLEN];
  char invcoefilex[DRMS_MAXPATHLEN], invcoefiley[DRMS_MAXPATHLEN];
  char outdir[DRMS_MAXPATHLEN], outdirx[DRMS_MAXPATHLEN],
      outdiry[DRMS_MAXPATHLEN], kernfilename[DRMS_MAXPATHLEN];
  char *carstr;
  char rec_query[DRMS_MAXQUERYLEN], source[DRMS_MAXQUERYLEN];
  char recset[DRMS_MAXQUERYLEN];
  char locname[24];
  char module_ident[64];

  int twosegs = 0;
  int maxmodect = 15000;
       /*  there are 14924 modes in the kernel file hmi.HSKernels[cm_opal78]  */

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
  char *avgkrndest = strdup (params_get_str (params, "ave"));
  char *invcoefdest = strdup (params_get_str (params, "coef"));
  double ob = params_get_double (params, "ob");
  double oe = params_get_double (params, "oe");
  double rb = params_get_double (params, "rb");
  double re = params_get_double (params, "re");
  double mu = params_get_double (params, "amu");
  int numtrg = params_get_int (params, "num");
  int lmax = params_get_int (params, "lmax");
  int nmax = params_get_int (params, "nmax");
  int verbose = params_isflagset (params, "v");

  int wantavgkrn = strcmp (avgkrndest, "Not Specified");
  int wantinvcoef = strcmp (invcoefdest, "Not Specified");
  int anycr = (cr <= 0);
  int anycl = isnan (cl);
  int anylon = isnan (lon);
  int anylat = isnan (lat);
  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);

	    /*  special hack to reflect limit in current version of autoweed  */
if (nmax > 8) {
  fprintf (stderr, "Warning: nmax reduced from %d to 8\n", nmax);
  nmax = 8;
}
  drms_input = 1;
		/*  attempt to open the kernel argument as a DRMS record set  */
  kern_set = drms_open_records (drms_env, kernel, &status);
  if (status) {
		      /*  attempt to open the kernel argument as a file name  */
    fpt = fopen (kernel, "r");
    if (!fpt) {
      fprintf (stderr, "Error: kernel specification \"%s\"\n", kernel);
      fprintf (stderr,
          "       does not specify a readable data record nor file\n");
      return 1;
    }
    fclose (fpt);
    strcpy (kernfilename, kernel);
  } else {
		       /*  check that the kernel set specifies unique record  */
    if (kern_set->n != 1) {
      fprintf (stderr, "Error: no data records in set %s\n", kernel);
      return 1;
    }
    irec = kern_set->records[0];
    drms_record_directory (irec, kernfilename, 1);
    segment = drms_segment_lookup (irec, "kernel");
    if (!segment) {
			       /*  and that it contains the segment "kernel"  */
      fprintf (stderr, "Error: kernel record \"%s\"\n", kernel);
      fprintf (stderr, "       does not contain a segment \"kernel\"\n");
      drms_close_records (kern_set, DRMS_FREE_RECORD);
      return 1;
    }
    sprintf (kernfilename, "%s/S%05i/kernel", segment->record->su->sudir,
	segment->record->slotnum);
    drms_close_records (kern_set, DRMS_FREE_RECORD);
  }

  strncpy (recset, in, DRMS_MAXQUERYLEN);
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
      snprintf (source, DRMS_MAXQUERYLEN, "[%d]", cr);
      strncat (recset, source, DRMS_MAXQUERYLEN);
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
      snprintf (source, DRMS_MAXQUERYLEN, "[%03.0f]", cl);
      strncat (recset, source, DRMS_MAXQUERYLEN);
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
      snprintf (source, DRMS_MAXQUERYLEN, "[%05.1f]", lon);
      strncat (recset, source, DRMS_MAXQUERYLEN);
    }
    if (!anylat) {
      float latmin = lat - 0.01;
      float latmax = lat + 0.01;
      snprintf (source, DRMS_MAXQUERYLEN, "[?LatHG>%g and LatHG<%g?]",
	  latmin, latmax);
      strncat (rec_query, source, DRMS_MAXQUERYLEN);
      if (anylon) strncat (recset, "[]", 2);
      snprintf (source, DRMS_MAXQUERYLEN, "[%+05.1f]", lat);
      strncat (recset, source, DRMS_MAXQUERYLEN);
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

  if (out[strlen(out) - 1] == '/') drms_output = 0;
  else drms_output = 1;
				   /*  Attempt to create output drms record  */
  if (drms_output) {
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
  } else strncpy (odir, out, DRMS_MAXPATHLEN);
  if (drms_output) {
    int kstat = 0;
    int keyct = sizeof (propagate) / sizeof (char *);

    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str   (orec, "Source", recset);
					       /*  module specific keywords  */
    kstat += check_and_set_key_float (orec, "amu", mu);
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
	    /*  allocate arrays dependent on number of targets for inversion  */
  rtrg = (double *)malloc (numtrg * sizeof (double));
  rcgx = (double *)malloc (numtrg * sizeof (double));
  rcgy = (double *)malloc (numtrg * sizeof (double));
  quartx = (double *)malloc (3 * numtrg * sizeof (double));
  quarty = (double *)malloc (3 * numtrg * sizeof (double));
  solnx = (double *)malloc (numtrg * sizeof (double));
  solny = (double *)malloc (numtrg * sizeof (double));
  errx = (double *)malloc (numtrg * sizeof (double));
  erry = (double *)malloc (numtrg * sizeof (double));
							/*  set target radii  */
  double drt = (re - rb) / (numtrg - 1);
  for (int target = 0; target < numtrg; target++)
    rtrg[target] = rb + target * drt;
						     /*  get mode amplitudes  */
  modect = read_kernel (kernfilename, &krnrval, &moden, &model, &krnamp,
      &kernrvals, &kernlloc, &kernlnct, 0, lmax);
  if (modect < 0) {
    fprintf (stderr, "Error: unable to read kernel file %s\n", kernfilename);
    drms_close_record (orec, DRMS_FREE_RECORD);
    if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
    return 1;
  }
       /*  allocate arrays dependent on number of kernel radii for inversion  */
  aa1x = (double *)malloc (kernrvals * kernrvals * sizeof (double));
  aa1y = (double *)malloc (kernrvals * kernrvals * sizeof (double));
  fakc = (double *)malloc (kernrvals * kernrvals * sizeof (double));
  avcx = (double *)malloc (kernrvals * sizeof (double));
  avcy = (double *)malloc (kernrvals * sizeof (double));
  domx = (double *)malloc (kernrvals * sizeof (double));
  domy = (double *)malloc (kernrvals * sizeof (double));
  errorx = (double *)malloc (kernrvals * sizeof (double));
  errory = (double *)malloc (kernrvals * sizeof (double));
  wmoden = (int *)malloc (kernrvals * sizeof (int));
  wmodel = (int *)malloc (kernrvals * sizeof (int));
	  /*  and those dependent on both number of kernel radii and targets  */
  avccx = (double *)malloc (numtrg * kernrvals * sizeof (double));
  avccy = (double *)malloc (numtrg * kernrvals * sizeof (double));
  ctil1x = (double *)malloc (numtrg * kernrvals * sizeof (double));
  ctil1y = (double *)malloc (numtrg * kernrvals * sizeof (double));
  ctil2x = (double *)malloc (numtrg * kernrvals * sizeof (double));
  ctil2y = (double *)malloc (numtrg * kernrvals * sizeof (double));
						 /*  set integration weights  */
  weight = (double *)malloc (modect * sizeof (double));
  weight[0] = 0.5 * (krnrval[0] - krnrval[1]);
  for (int j = 0; j < modect - 2; j++)
    weight[j+1] = 0.5 * (krnrval[j] - krnrval[j+2]);
  weight[modect-1] = 0.5 * (krnrval[modect-2] - krnrval[modect-1]);
						    /*  main processing loop  */
  for (rec_i=0; rec_i<nrec; rec_i++) {
    if (verbose) printf ("  Processing record %i : ", rec_i);
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
    if (verbose) printf ("%s\n", source);
	   /*  create appropriate output file names based on target location  */
    sprintf (locname, "%s", printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
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
    if (wantavgkrn) {
	    /*  check writeability of output directory for averaging kernels  */
      FILE *tstfile;
      sprintf (avgkrnfilex, "%s/%s.Ux", avgkrndest, locname); 
      sprintf (avgkrnfiley, "%s/%s.Uy", avgkrndest, locname);
      tstfile = fopen (avgkrnfilex, "w");
      if (!tstfile) {
	fprintf (stderr, "Warning: unable to open file %s\n", avgkrnfilex);
	fprintf (stderr, "         averaging kernels will not be written\n");
	wantavgkrn = 0;
      } else {
	fclose (tstfile);  
	tstfile = fopen (avgkrnfiley, "w");
	if (!tstfile) {
	  fprintf (stderr, "Warning: unable to open file %s\n", avgkrnfiley);
	  fprintf (stderr, "         averaging kernels will not be written\n");
	  wantavgkrn = 0;
	} else fclose (tstfile);
      }  
    }
    if (wantinvcoef) {
       /*  check writeability of output directory for inversion coefficients  */
      FILE *tstfile;
      sprintf (invcoefilex, "%s/%s.Ux", invcoefdest, locname); 
      sprintf (invcoefiley, "%s/%s.Uy", invcoefdest, locname);
      tstfile = fopen (invcoefilex, "w");
      if (!tstfile) {
	fprintf (stderr, "Warning: unable to open file %s\n", invcoefilex);
	fprintf (stderr, "         inversion coefficients will not be written\n");
	wantinvcoef = 0;
      } else {
	fclose (tstfile);  
	tstfile = fopen (invcoefiley, "w");
	if (!tstfile) {
	  fprintf (stderr, "Warning: unable to open file %s\n", invcoefiley);
	  fprintf (stderr, "         inversion coefficients will not be written\n");
	  wantinvcoef = 0;
	} else fclose (tstfile);
      }  
    }

    status = process_record (filename, outfilex, outfiley,
    	krnrval, weight, kernrvals, moden, model, krnamp,
	modect, kernel, nmax, lmax, kernlloc, kernlnct,
	avgkrnfilex, avgkrnfiley, wantavgkrn,
	invcoefilex, invcoefiley, wantinvcoef,
	wmoden, wmodel, ob, oe,
	numtrg, rb, re, mu, source, module_ident);
    if (status)	{
      fprintf (stderr, "Error processing record %d (%s); aborted\n", rec_i,
	  printcrcl_loc (cr, loncCar, lonhg, loncm, latc));
      if (drms_output) {
	if (drms_input) drms_close_records (recordSet, DRMS_FREE_RECORD);
	drms_close_record (orec, DRMS_FREE_RECORD);
	return 1;
      }
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
 *	reversal of east and west), and dropped requirement for CMLon;
 *		Heliographic option for file naming
 *		Removed default for kernel
 *    0.9	Created 2011.12.04; Fixed typos in ola_xy.c (not yet used)
 *  v 0.8 frozen 2012.02.01 (but included in JSOC release 6.1 2011.12.05 in
 *		Identical form except for comments)
 *    0.9	Added (preferred) option for specifying kernel as DRMS data
 *	record rather than filename
 *  v 0.9 frozen 2012.04.24
 *    0.91	Similar to 0.9, but includes fixes in FORTRAN code ola_xy_v12
 *		and provision for return status
 *  v 0.91 frozen 2012.04.24
 *    0.92	Created 2018.03.10; similar to 0.91 but including version
 *	control for included source files; includes rdutil_v09.c rather than
 *	old_rdutil.c; also intended to be built with ola_subs_v13.f and
 *	ola_xy_v13.f, although this can only be controlled through the
 *	Makefile
 *		Added setting of source key
 *	2018.03.24 Added verbose reporting of source record (or file) name at
 *		each step
 *	2018.04.13 Removed unused argument drms_output from call to
 *		process_record()
 *	2018.04.18 Added version control for all included C source
 *	2018.04.25 Removed default value for out argument
 *  v 0.92 frozen 2018.05.07
 *    0.93	Created 2018.11.20; similar to 0.92, but including olaxy_v13.c
 *	and eliminating functions in ola_xy*.f (Fortran compilation no longer
 *	required); removed unused variables associated with C/Fortran interface
 *		Removed verbose flag from call to ola_xy
 *		Moved reading of kernel values out of process_record() and
 *	into main module, so only done once
 *		Added trap for output name ending in '/' to not try as DRMS
 *		Moved kernel read to own function (called once from main)
 *		Added trap for array limits of lmax and nmax in ola_xy, setting of
 *	relevant limits as runtime arguments and and passing of array limits
 *	as argument; reduced default limits from 50,3000, to 14,2000
 *		Save time (~10 sec) by not reading through kernel file twice
 *	to get exact array size, but increment by 1000 (there are ~15000 mode
 *	amplitudes in the standard kernel file)
 *		Moved optional writing of averaging kernel file from ola_xy()
 *	to process_record(); changed format of averaging kernel files; changed
 *	specification for optionally written file from filename to directory
 *	name; in writing of averaging kernel, changed from writing values for
 *	every other radial target of kernel to writing every one.
 *		Moved optional writing of inversion coefficients from ola_xy()
 *	to process_record() and modified format; modified kernel file reading
 *	to keep track of locations in kernel array for each degree (and order
 *	counts for same); removed trap for ola_xy array limits no longer set;
 *	moved allocation of most arrays from ola_xy to main module; declared
 *	them global
 *		Added special hack overriding nmax>8 to 8
 *  v 0.93 frozen 2019.01.31
 */
