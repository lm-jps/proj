/*
 *  rdfitf.c					~rick/src/ringfit/rdfitf
 *
 *  The "fast" ring-fitting procedure based on fitting only a few parameters
 *    to the low-order azimuthal Fourier components in k-space of the power
 *    spectra, rather than all m values
 *
 *  DESCRIPTION NEEDS TO BE REVIEWED
 *
 *  Usage:
 *    rdfitf [-nv] [param= ...]
 *
 *  Parameters: 	(type   default         description)
 *	in	DataSet	-	Input dataset
 *	out	string	fit.out	output data series or (if series is not
 *				writeable) file name
 *	guessfile string -	name of file containing parameterizations for
 *				each fitting parameter
 *	unwrap	string	unwrap		name of the segment to which
 *		the unwrapped input power spectrum (before filtering) is to be
 *		written (if -U flag set)
 *	filter	string	filter	name of the segment to which the filter profile
 *				is to be written (if -f flag set)
 *	filtps	string	filtps	name of the segment to which the unwrapped,
 *				filtered, and subsampled input power spectrum
 *				is to be written (if -u flag set)
 *	nmin	int	0	lowest radial order to fit
 *	nmax	int	7	highest radial order to fit
 *	npar	int	6	number of parameters to fit
 *	ntht	int	400	number of azimuth values for initial unwrapping
 *		of spectrum, usually = 2*pi*nk/2; should be divisible by 4
 *	nthts	int	40	number of azimuth values after subsampling;
 *				generally a factor of 10 less than ntht,
 *				should also be divisible by 4
 *	ux	float	0.0	initial guess for zonal flow speed [m/s]
 *	uy	float	0.0	" " for meridional flow speed [m/s]
 *	smooth	int	0	amount of smoothing wanted [wavenumber bins]
 *	mdata	int	16	number of azimuthal Fourier components to keep
 *				when subsampling the data
 *	mfilt	int	8	highest azimuthal Fourier component to use in
 *				filtering the data
 *	kstart	ints	{11, 7,  6,  6,  6,  6,  7,  8}
 *				starting values for k-bins for each n
 *	kstop	ints	{80, 70, 55, 44, 37, 27, 23, 20}
 *				ending values for k-bins
 *	xtol	double	1.e-4	termination tolerance on fitting parameters
 *	ftol	double	1.e-8	termination tolerance on fitting function
 *	feps	double	1.e-16	accuracy of function values from FUNC
 *	delta0	double	0.0	initial radius of the trust region
 *	maxeval	int	750	maximum number of iterations for each fit
 *	igrad	int	1
 *	iextrap	int	1
 *	idiag	int	0
 *	iprint	int	0
 *
 *  Flags
 *	-f	write spatial wavenumber filter to appropriate segment
 *	-n	do not write out fits (diagnostic only)
 *	-u	write filtered, subsampled, unwrapped power spectrum to
 *		appropriate segment
 *	-U	write unwrapped power spectrum to appropriate segment
 *	-v	run verbose
 *	-x	run extra verbose
 *	-2	skip 2nd derivatives in parameter errors
 *
 *  Bugs:
 *    Several of the named value arguments are really just flags or enums
 *    Writing of intermediate data files (unwrapped spectra etc.) is only
 *	supported in DRMS, there is no direct FITS write to named file
 *    The test on lnbase in cart_to_polar is probably wrong, since lnbase
 *	will be returned as NaN rather than 0 if the keyword is missing
 *    Should use WCS keywords if present
 *    It is not clear what happens if npar > 6; module fails for npar < 6.
 *    Most or all of the internally written information when the verbose flag
 *	is set should be written to the Log
 *    The default values of kstart and kstop refer to the nunmber of kbins, and
 *	need to be adjusted for different sizes of tiles; they should be
 *	calculated not read in as parameters (or from a file as formerly)
 *    The decision as to whether the output specifier refers to a data series
 *	or a file name based on whether it can be referred to a series
 *	writeable by the user is not sensible, but probably not particularly
 *	dangerous.
 *    It is implictly assumed that in the input spectrum, NAXIS1 = NAXIS2,
 *	CRPIX1 = CRPIX2, CRVAL1 = CRVAL2, and CDELT1 = CDELT2; only the first
 *	and second of these are checked. Also, it is required that NAXIS1 and
 *	NAXIS2 are even
 *
 *  Future development:
 *    Process LOG_BASE key in driver, converting power spectrum to log if
 *	necessary, and remove lnbase from arguments to Fortran driver (and
 *	processing alternative in cart_to_polar)
 *    Add reprocess flag, and sufficient information in keywords from arg
 *	list to do so
 *
 *  Notes:
 *    The following parameters are passed between the driver module and the
 *    main interface Fortran function (strings [*] require corresponding string
 *    lengths):
 *	powr		array of floats, read from the input segment
 *	nkx, nky, nnu	dimensions of input power spectrum
 *      ntht, nthts, nk	dimensions of unwrapped power spectra
 *	crpix1, crpix2	reference location of center of input spectrum
 *	smooth		wavenumber smoothing parameter
 *	dk, dnu		scale factors for the axes, from the input record
 *	lnbase		logarithmic base for input power spectrum scaling
 *	nmin, nmax	limits on radial orders to fit
 *	npar		number of parameters to be fit
 *	ux, uy		initial guess values for Ux & Uy parameters
 *	inds [*]	input data record descriptor or data file name
 *      guessfile [*]	name of the file containing guess frequencies,
 *		amplitudes, widths, and background parameters
 *	unwrapped_pspec	float* of unwrapped power spectrum (returned)
 *	filter		float* of wavenumber filter (returned)
 *	filtered_pspec	double* of filtered power spectrum (returned)
 *	verbose		verbose flag value
 *	mdata		number of subsampling Fourier components
 *	mfilt		maximum Fourier component to use in filtering
 *	doptions	(see below)
 *	ioptions	(see below)
 *	iflg
 *	status		return value from call		
 *
 *    doptions(1) = xtol   :  double precision
 *      Termination tolerance on X.
 *
 *    doptions(2) = ftol   :  double precision
 *      Termination tolerance on F.
 *
 *    doptions(3) = feps   :  double precision
 *      (Approximately) The accuracy of the function values
 *      calculated in FUNC.
 *
 *    doptions(4) = delta0 :  double precision
 *      Initial radius of the trust region.
 *
 *    ioptions(1) = maxeval  :  integer
 *      On input  : Maximum number of function evaluations.
 * 
 *    ioptions(2) = igrad  :  integer
 *      If igrad is 0 then F' is calculated numerically.
 *      If igrad is 1 then the user supplied subroutine DF1 is called
 *      by subroutine dogleg from subroutine hmifits to calculate F'.
 *      If igrad is 2 then the user supplied subroutine DF1 is called
 *      to calculate F' and initially the results from GRAD are checked
 *      by comparing them to the F' obtained using finite differences.
 *
 *    ioptions(3) = iextrap  :  integer
 *      If iextrap is 1 then safeguarded quadratic/cubic interpolation
 *      and extrapolation is performed.
 *
 *    ioptions(4) = idiag  :  integer
 *      If idiag is 0 then I is used as the initial approximation to the
 *      Hessian.
 *      If idiag is 1 then the diagonal of the true Hessian is computed
 *      initially using a finite difference approximation.
 *
 *    ioptions(5) = iprint  :  integer
 *      If iprint=1 then information is printed after each iteration.
 *
 *  Revision history is at end of file.
 */

#define MODULE_VERSION_NUMBER	("0.9")
#define KEYSTUFF_VERSION "keystuff_v11.c"

char *module_name = "ringfit_ssw";
char *module_desc = "ring fitting using method of Haber and Hindman";
char *version_id = MODULE_VERSION_NUMBER;

#include <jsoc_main.h>
#include KEYSTUFF_VERSION

ModuleArgs_t module_args[] = {
  {ARG_DATASET,	"in", "", "input dataset"}, 
  {ARG_STRING, "out", "fit.out",
      "name of output data series (or file) containing fit coefficients"},
  {ARG_STRING, "guessfile", "",
      "name of file containing parameterizations for each fitting parameter"},
  {ARG_STRING, "unwrap", "unwrap",
      "segment containing optional unwrapped power spectra"},
  {ARG_STRING, "filter", "filter",
      "segment containing optional spatial wavenumber filter"},
  {ARG_STRING, "filtps", "filtps",
      "segment containing filtered, subsampled, power spectra"},
  {ARG_INT, "nmin", "0", "lowest radial order to fit", "[0,)"},
  {ARG_INT, "nmax", "7", "highest radial order to fit", "[0,)"},
  {ARG_INT, "npar", "6", "number of parameters to fit", "[6,)"},
  {ARG_INT, "ntht", "400",
      "number of azimuth values for initial unwrapping of spectrum, usually = 2*pi*nk, should be divisible by 4"},
/* should be calculated in code  */
  {ARG_INT, "nthts", "40",
      "number of azimuth values after subsampling, generally a factor of 10 less than ntht, should also be divisible by 4"},
/* should be calculated in code  */
  {ARG_INT, "nrdtot", "16",
      "number of ridges in guess file. as long as there are 16 or less the default should work"},
  {ARG_FLOAT, "ux", "0.0", "initial guess for zonal flow speed [m/s]"},
  {ARG_FLOAT, "uy", "0.0", "initial guess for meridional flow speed [m/s]"},
  {ARG_INT, "smooth", "0", "amount of smoothing wanted [wavenumber bins]"},
  {ARG_INT, "mdata", "16",
      "number of azimuthal Fourier components to keep when subsampling the data"},
  {ARG_INT, "mfilt", "8",
      "highest azimuthal Fourier component to use in filtering the data"},
  {ARG_INTS, "kstart", "{11, 7,  6,  6,  6,  6,  7,  8}",
      "starting kbin values"},
  {ARG_INTS, "kstop",  "{55, 55, 45, 44, 37, 27, 23, 20}",
      "ending kbin values"},
  {ARG_FLAG,  "f", "0", "write spatial wavenumber filter"},      
  {ARG_FLAG,  "n", "0", "do not write out fits (diagnostics only)"},      
  {ARG_FLAG,  "U", "0", "write unwrapped power spectrum"},      
  {ARG_FLAG,  "v", "0", "verbose output"}, 
  {ARG_FLAG,  "x", "0", "really verbose output"}, 
  {ARG_FLAG,  "d", "0", "do not use second derivatives in parameter errors"},      
  {ARG_FLOAT, "xtol", "1e-4", "termination tolerance on fitting parameters"}, 
  {ARG_FLOAT, "ftol", "1e-8", "termination tolerance on fitting function"}, 
  {ARG_FLOAT, "feps", "1e-16", "accuracy of function values from FUNC"}, 
  {ARG_FLOAT, "delta0", "0.0", "initial radius of the trust region"}, 
  {ARG_INT, "maxeval", "750", "maximum number of iterations for each fit"},
  {ARG_INT, "igrad", "1", "how to determine gradient, 1=user supplied DF1"},
  {ARG_INT, "iextrap", "1", "if iextrap is 1 then safeguarded quadratic/cubic interpolation and extrapolation is performed"},
  {ARG_INT, "idiag", "0", "if idiag is 0 then I is used as the initial approximation to the Hessian"},
  {ARG_INT, "iprint", "0", "if iprint=1 then information is printed after each iteration"},
  {ARG_END}
};

extern void ringanalysis_ (float *powr, int *nkx, int *nky, int *nnu, 
    int *ntht, int *nthts, int *nk, float *crpix1, float *crpix2, int *smooth,
    double *dk, double *dnu, float *lnbase, int *nmin, int *nmax, 
    int *npar, int *ntot, double *ux_guess, double *uy_guess, float *ux_fit, 
    float *d_ux, float *uy_fit, float *d_uy, float *amp, float *d_amp, float *bg, 
    float *d_bg, float *fwhm, float *d_fwhm, float *nu, float *d_nu,
    float *delnu, int *nval, float *kval, int *kbin, int *kbstrt, int *kbend,
    int *kbmin, int *kbmax, int *nrdg, float *min_func, int *n_iter, int *fit_chk,
    float *chi, float *kmin, float *kmax, int *modecount, int *good, float *unwrapped_pspec, 
    float *filter, float *filtered_pspec, double *fnu, double *width, double *ampg,
    double *bkg, int *verbose, int *mdata, int *mfilt, double *doptions, int *ioptions, 
    int *iflg, int *rgn, int *nrdtot, int *status);

       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "Cadence", "LonSpan",
    "T_START", "T_STOP", "Coverage", "Size", "Width", "Height",
    "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel",
    "MapScale", "MapProj", "Map_PA", "RSunRef", "PosAng", "MAI",
    "Apode_f", "APOD_MIN", "APOD_MAX"};

#include "read_guess.c"

#define RSUN_MM (696.0)
#define NOFITS	(0x80000000)

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *pspec;
  DRMS_Array_t *filtps_array = NULL, *filter_array = NULL, *unwrap_array = NULL; 
  FILE *unit24, *runlog;
  float *filtered_pspec;
  double doptions[4]; 
  double dk, dnu, dklast;
  double *fnu, *width, *ampg, *bkg;
  double *fnua, *fwa, *pa, *bga;
  float *spec;
  float *unwrapped_pspec, *filter;
  float *ux_fit, *d_ux, *uy_fit, *d_uy;
  float *amp, *d_amp, *bg, *d_bg, *fwhm, *d_fwhm;
  float *nu, *d_nu, *delnu, *chi, *min_func, *kval;
  float lnbase, crpix1, crpix2;
  float rsun, kmin, kmax, lmin, lmax, ell, nkc;
  unsigned int quality;
  int *nval, *n_iter, *fit_chk, *good, *nvct;
  int *kbstrt, *kbend, *kbin;
  int kbmin, kbmax;
  int ioptions[5], m[4];
  int modecount;
  int rgn, rgnct, segct, drms_output, dispose;
  int isegnum, osegnum, logsegnum, filtsegnum;
  int nmalloc, nloopmalloc;
  int unwrapsegnum, filtpowsegnum;
  int n, nkx, nky, nnu, nk, ntot, ma; 
  int nrdg, i, j;
  int status;
  char logfile[DRMS_MAXPATHLEN], outfile[DRMS_MAXPATHLEN];
  char recid[DRMS_MAXQUERYLEN];
  char module_ident[64], key[64];

  char *inds = strdup (params_get_str (params, "in"));
  char *outser =  strdup (params_get_str (params, "out"));
  char *guessfile =  strdup (params_get_str (params, "guessfile"));
  char *unwrap_segment =  strdup (params_get_str (params, "unwrap"));
  char *filter_segment =  strdup (params_get_str (params, "filter"));
  char *filtered_segment =  strdup (params_get_str (params, "filtps"));
  int nmin = params_get_int (params, "nmin");
  int nmax = params_get_int (params, "nmax");
  int npar = params_get_int (params, "npar");
  int ntht = params_get_int (params, "ntht");
  int nthts = params_get_int (params, "nthts");
  int nrdtot = params_get_int (params, "nrdtot");
  double ux_guess = params_get_double (params, "ux");
  double uy_guess = params_get_double (params, "uy");
  int smooth = params_get_int (params, "smooth");
  int kstrtct = params_get_int (params, "kstart_nvals");
  int kstopct = params_get_int (params, "kstop_nvals");

  int write_filter = params_isflagset (params, "f");
  int no_fits = params_isflagset (params, "n");
  int write_unwrap_pow = params_isflagset (params, "U");
  int write_filt_pow = params_isflagset (params, "u");
  int verbose = params_isflagset (params, "v");
  int Verbose = params_isflagset (params, "x");

  int mdata = params_get_int (params, "mdata");
  int mfilt = params_get_int (params, "mfilt");
  double xtol = params_get_double (params, "xtol");
  double ftol = params_get_double (params, "ftol");
  double feps = params_get_double (params, "feps");
  double delta0 = params_get_double (params, "delta0");
  int maxeval = params_get_int (params, "maxeval");
  int igrad = params_get_int (params, "igrad");
  int iextrap = params_get_int (params, "iextrap");
  int idiag = params_get_int (params, "idiag");
  int iprint = params_get_int (params, "iprint");
  int iflg = params_isflagset (params, "2");

  snprintf (module_ident, 64, "%s v%s", module_name, version_id);
  if (verbose) {
    printf ("%s:\n", module_ident);
    printf ("fits of %s\n", inds);
    if (no_fits)
      printf ("(diagnostic run only, no records will be written to DRMS)\n");
  }
	  /*  the following are needed strictly for the C/Fortran interface  */
  ioptions[0] = maxeval;
  ioptions[1] = igrad;
  ioptions[2] = iextrap;
  ioptions[3] = idiag;
  ioptions[4] = iprint;
  doptions[0] = xtol;
  doptions[1] = ftol;
  doptions[2] = feps;
  doptions[3] = delta0;
				    /*  sanity checks on certain parameters  */
  if (nmax + 1 > nrdtot) {
    fprintf (stderr, "Error: not enough ridges in guess file\n");
    return 0;
  }
  if (kstopct != kstrtct) {
    fprintf (stderr, "Warning: kstop vals (%d) != kstart vals (%d);", kstopct,
        kstrtct);
    if (kstrtct < kstopct) kstopct = kstrtct;
    else kstrtct = kstopct;
    fprintf (stderr, " reduced to %d\n", kstopct);
  }
  nrdg = kstopct;
  if (nrdg <= nmax) {
    fprintf (stderr,
        "Warning: insufficient kstart/kstop values for requested range\n");
    fprintf (stderr, "         nmax reduced from %d to %d\n", nmax, nrdg - 1);
    nmax = nrdg - 1;
    if (nrdg <= nmax) {
      fprintf (stderr, "Error: still insufficient kstart/kstop values \n");
      return 0;
    }
  }
  if (nmin > nmax) {
    fprintf (stderr, "Warning: nmax %d < nmin %d)\n", nmax, nmin);
    fprintf (stderr, "setting nmax=nmin");
    nmax = nmin;
  }
  if (ntht < nthts) {
    fprintf (stderr, "Error: ntht %d < nthts %d\n", ntht, nthts);
    return 1;
  }
  if (npar != 6) {
    fprintf (stderr, "Error: only 6 fitting parameters possible at this time\n");
    fprintf (stderr, "You'll need to change subroutines and guess table\n");
    fprintf (stderr, "if you want to try other parameters or functional forms\n");
    return 0;
  }

  nvct = (int *)malloc ((nmax - nmin + 1) * sizeof (int));
  kbstrt = (int *)malloc (nrdg * sizeof (int));
  kbend = (int *)malloc (nrdg * sizeof (int));
  for (n = 0; n < nrdg; n++) {
    snprintf (key, 64, "kstart_%d_value", n);
    kbstrt[n] = params_get_int (params, key);
    snprintf (key, 64, "kstop_%d_value", n);
    kbend[n] = params_get_int (params, key);
  }
						   /*  open and check input  */   
  ids = drms_open_records (drms_env, inds, &status);
  if (status) {
    fprintf (stderr, "Error: drms_open_records returned %d for dataset:\n", status);
    fprintf (stderr, "  %s\n", inds);
    return 0;
  }
  rgnct = ids->n;
  if (rgnct < 1) {
    fprintf (stderr, "dataset %s contains no records\n", inds);
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 0;
  }
		/*  this code is shared with code in rdfitc,
			      and should be turned into a separate function  */
					/*  check output data series struct  */
  orec = drms_template_record (drms_env, outser, &status);
  if (orec) {
    segct = drms_record_numsegments (orec);
    if (segct < 1) {
      fprintf (stderr,
	  "Error: no data segment in output series %s\n", outser);
      drms_close_records (ids, DRMS_FREE_RECORD);
      drms_free_record (orec);
      return 1;
    }
    oseg = drms_segment_lookup (orec, "fit.out");
    if (oseg && oseg->info->protocol == DRMS_GENERIC)
      osegnum = oseg->info->segnum;
    else {
    /*  either segment fit.out does not exist, or it has the wrong protocol  */
      for (n = 0; n < segct; n++) {
	oseg = drms_segment_lookupnum (orec, n);
	if (oseg->info->protocol != DRMS_GENERIC) continue;
	osegnum = n;
   	break;
      }
      if (n == segct) {
	fprintf (stderr,
	    "Error: no generic data segment in output series %s\n", outser);
	drms_close_records (ids, DRMS_FREE_RECORD);
	drms_free_record (orec);
	return 1;
      }
      fprintf (stderr, "         writing to segment %s.\n", oseg->info->name);
    }
    logsegnum = -1;
    oseg = drms_segment_lookup (orec, "Log");
    if (oseg && oseg->info->protocol == DRMS_GENERIC)
      logsegnum = oseg->info->segnum;
					  /*  this code is unique to rdfitf  */
    if (write_unwrap_pow) {
      oseg = drms_segment_lookup (orec, unwrap_segment);
      if (oseg && oseg->info->protocol != DRMS_GENERIC && oseg->info->naxis == 3)
	unwrapsegnum = oseg->info->segnum;
      else {
	fprintf (stderr,
	    "Warning: unwrapped spectrum requested, but no appropriate segment\n");
	fprintf (stderr, "  %s in output series %s\n", unwrap_segment, outser);
	write_unwrap_pow = 0;
      }
    }
    filtsegnum = -1;
    if (write_filter) {
      oseg = drms_segment_lookup (orec, filter_segment);
      if (oseg && oseg->info->protocol != DRMS_GENERIC && oseg->info->naxis == 2)
	filtsegnum = oseg->info->segnum;
      else {
	fprintf (stderr,
	    "Warning: filter file requested, but no appropriate segment\n");
	fprintf (stderr, "  %s in output series %s\n", filter_segment, outser);
	write_filter = 0;
      }
    }
    if (write_filt_pow) {
      oseg = drms_segment_lookup (orec, filtered_segment);
      if (oseg && oseg->info->protocol != DRMS_GENERIC && oseg->info->naxis == 3)
	filtpowsegnum = oseg->info->segnum;
      else {
	fprintf (stderr,
	    "Warning: filtered spectrum requested, but no appropriate segment\n");
	fprintf (stderr, "  %s in output series %s\n", filtered_segment, outser);
	write_filt_pow = 0;
      }
    }
    drms_free_record (orec);
    drms_output = 1;
    if (verbose)
        printf("rgnct = %d \n",rgnct);
    if (verbose && (rgnct > 1)) printf ("%d records in input dataset %s\n",
        rgnct, inds);
  } else {
	    /*  Can't create in named series, assume it's a filename instead  */
    if (rgnct > 1) {
      fprintf (stderr,
	  "Query produced %d matching records with output to single file;",
	  rgnct);
      fprintf (stderr, " must limit to 1\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    strcpy (outfile, outser);
    drms_output = 0;
  }
  dispose = (no_fits) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
					     /*  end shared code with rdfitc  */
		    /*  Read in guess table coefficients for main parameters  */
  irec = ids->records[0];
  dk = drms_getkey_double (irec, "DELTA_K", &status);
  nkc = drms_getkey_float (irec, "CRPIX1", &status);
  if (nkc != drms_getkey_float (irec, "CRPIX2", &status)) {
    fprintf (stderr,
        "Error: CRPIX1 and CRPIX2 in input spectrum must be equal\n");
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
		     /*  original code expected crpix1,2 to be 0.5 too small  */
  nkc -= 0.5;
  nk = nkc - 0.5;
  status = read_guess (guessfile, nrdtot, m, &fnua, &fwa, &pa, &bga, &ma,
      verbose);
  if (status) {
    fprintf (stderr, "Error: problem reading guessfile %s\n", guessfile);
    fprintf (stderr, "       read_guess() returned %d\n", status);
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }

  fnu = (double *)calloc (nk * nrdtot , sizeof (double));
  width = (double *)calloc (nk * nrdtot , sizeof (double));
  ampg = (double *)calloc (nk * nrdtot , sizeof (double));
  bkg = (double *)calloc (nk * nrdtot , sizeof (double));
  ntot = nmax - nmin + 1;
  kbmin = min_arr(kbstrt, ntot);
  kbmax = max_arr(kbend, ntot);
  status = make_table (fnua, fwa, pa, bga, fnu, width, ampg, bkg, m, nk, dk,
      nrdtot, kbmin, kbmax, ma, verbose);
  if (status) {
    fprintf (stderr, "Error: problem making table of guesses\n");
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
  dklast = dk;
	 /*  loop over all input spectra, creating an output record for each  */
			 /*  check to make sure you are at first good record  */
  nmalloc = 0;
  for (rgn = 0; rgn < rgnct; rgn++) {
    irec = ids->records[rgn];
    drms_sprint_rec_query (recid, irec);
    if (verbose && (rgnct > 1))
      printf ("processing record %d:\n    %s\n", rgn, recid);
    dk = drms_getkey_double (irec, "DELTA_K", &status);
    dnu = drms_getkey_double (irec, "DELTA_NU", &status);
    lnbase = drms_getkey_float (irec, "LOG_BASE", &status);
    crpix1 = drms_getkey_float (irec, "CRPIX1", &status);
    crpix2 = drms_getkey_float (irec, "CRPIX2", &status);
    if (crpix2 != crpix1) {
      fprintf (stderr, "Error: CRPIX1 != CRPIX2 in input spectrum\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
		       /*  Fortran code expects crpix1,2 to be 0.5 too small  */
    crpix1 -= 0.5;
    crpix2 -= 0.5;
    rsun = drms_getkey_float (irec, "RSunRef", &status);
    if (status || isnan (rsun)) {			     /*  from rdfitc  */
      rsun = RSUN_MM;
      fprintf (stderr, "Warning: no valid value for RSunRef; ");
      fprintf (stderr, "using %f\n", rsun);
    }
    if (dk != dklast) {
      status = make_table(fnua, fwa, pa, bga, fnu, width, 
                ampg, bkg, m, nk, dk, nrdtot, kbmin, kbmax, ma, verbose);
      if (status) {
	printf ("Problem making table of guesses \n");
	return 1;
      }
    }
    dklast = dk;
    if (verbose) printf("dk = %lf   dklast = %lf \n", dk, dklast);

    segct = irec->segments.num_total;
    if (segct != 1) {
      for (n = 0; n < segct; n++) {
	iseg = drms_segment_lookupnum (irec, n);
        if (!iseg) continue;
	if (iseg->info->naxis == 3) break;
      }
      if (n >= segct) {
	fprintf (stderr, "found no segmemt of rank 3 in input dataset\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
        return 1;
      }
      isegnum = n;
    } else {
      isegnum = 0;
      iseg = drms_segment_lookupnum (irec, isegnum);
      if (!iseg) {
	fprintf (stderr, "Error, could not open data segment\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
      if (iseg->info->naxis != 3) {
	fprintf (stderr, "Error: found no segmemt of rank 3 in input dataset\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
    }
    if (!iseg) {
      fprintf (stderr, "Error, could not open data segment\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    pspec = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
    if (status) {
      fprintf (stderr,
	  "Error, could not read data segment from input record %d\n", rgn);
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
				     /*  these are now being used as of v07.c */
    nkx = pspec->axis[0];
    nky = pspec->axis[1];
    nnu = pspec->axis[2];
    spec = (float *)pspec->data;
    if (nkx != nky) {
      fprintf (stderr, "Error: NAXIS1 = %d != NAXIS2 = %d\n", nkx, nky);
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    if (nkx % 2) {
      fprintf (stderr, "Error: NAXIS1 and NAXIS2 must be even\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    nk = nkx / 2;
    if (verbose)
      printf("nkx = %d  nky = %d  nnu = %d  nk = %d  ntot = %d\n",nkx,nky,nnu,nk,ntot);
    for (n = 0; n < nrdg; n++) {
      if (kbend[n] > nk) {
	fprintf (stderr, "Error: kstop value %d > nk = %d \n", kbend[n], nk);
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
    }
    if (nmalloc == 0) {
      unwrapped_pspec = (float *)malloc (ntht * nk * nnu * sizeof (float));
      filter = (float *)malloc (nthts * nk * sizeof (float));
      ux_fit = (float *)malloc (nk * ntot * sizeof (float));
      d_ux = (float *)malloc (nk * ntot * sizeof (float));
      uy_fit = (float *)malloc (nk * ntot * sizeof (float));
      d_uy = (float *)malloc (nk * ntot * sizeof (float));
      amp = (float *)malloc (nk * ntot * sizeof (float));
      d_amp = (float *)malloc (nk * ntot * sizeof (float));
      bg = (float *)malloc (nk * ntot * sizeof (float));
      d_bg = (float *)malloc (nk * ntot * sizeof (float));
      fwhm = (float *)malloc (nk * ntot * sizeof (float));
      d_fwhm = (float *)malloc (nk * ntot * sizeof (float));
      nu = (float *)malloc (nk * ntot * sizeof (float));
      d_nu = (float *)malloc (nk * ntot * sizeof (float));
      delnu = (float *)malloc (nk * ntot * sizeof (float));
      nval = (int *)malloc (nk * ntot * sizeof (int));
      kval = (float *)malloc (nk * ntot * sizeof (float));
      kbin = (int *)malloc (nk * ntot * sizeof (int));
      good = (int *)malloc (nk * ntot * sizeof (int));
      min_func = (float *)malloc (nk * ntot * sizeof (float));
      n_iter = (int *)malloc (nk * ntot * sizeof (int));
      fit_chk = (int *)malloc (nk * ntot * sizeof (int));
      chi = (float *)malloc (nk * ntot * sizeof (float));
      filtered_pspec = (float*)malloc (nthts * nk * nnu * sizeof (float));
      nmalloc = ntht * nk * nnu;
    } else { 
      nloopmalloc = ntht * nk * nnu;
      if (nloopmalloc != nmalloc) {
        fprintf (stderr, "Warning: bad record, go to next one");
        fprintf (stderr, "  dimensions are not the same"); 
        continue;
      }
      
    }

    if (drms_output) {
      int kstat = 0;
      int keyct = sizeof (propagate) / sizeof (char *);
      char *key_str, *suffix;
      double key_dbl;
      int key_int, crcl_known = 1;
      TIME key_tim;
      
      orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
      oseg = drms_segment_lookupnum (orec, osegnum);
      drms_segment_filename (oseg, outfile);
				  /*  strip off the .generic tag if present  */
      suffix = strstr (outfile, ".generic");
      if (suffix) *suffix = '\0';
      if (logsegnum >= 0) {
        oseg = drms_segment_lookupnum (orec, logsegnum);
        drms_segment_filename (oseg, logfile);
				  /*  strip off the .generic tag if present  */
        suffix = strstr (logfile, ".generic");
        if (suffix) *suffix = '\0';
      }
		/*  this code is or ought to be shared with code in rdfitc,
			      and should be turned into a separate function  */
      kstat += propagate_keys (orec, irec, propagate, keyct);
	     /*  if necessary, construct CarrTime from CR:CL, or vice-versa  */
      key_str = drms_getkey_string (irec, "CarrTime", &status);
      if (status) {
        key_int = drms_getkey_int (irec, "CarrRot", &status);
	if (status) crcl_known = 0;
	else {
	  key_dbl = drms_getkey_double (irec, "CMLon", &status);
	  if (status) crcl_known = 0;
	}
	if (crcl_known) {
	  char CarrTime[9];
	  snprintf (CarrTime, 9, "%d:%03.0f", key_int, key_dbl);
	  kstat += check_and_set_key_str (orec, "CarrTime", CarrTime);
	}
      } else {
	key_int = drms_getkey_int (irec, "CarrRot", &status);
	if (status) {
	  sscanf (key_str, "%d:%lf", &key_int, &key_dbl);
	  kstat += check_and_set_key_int (orec, "CarrRot", key_int);
	  kstat += check_and_set_key_double (orec, "CMLon", key_dbl);
	}
      }
      kstat += check_and_set_key_str (orec, "Module", module_ident);
      kstat += check_and_set_key_str (orec, "BLD_VERS", jsoc_version);
      kstat += check_and_set_key_str (orec, "Source", recid);
      kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
					    /*  end shared code with rdfitc  */
      kstat += check_and_set_key_int (orec, "n_min", nmin);
      kstat += check_and_set_key_int (orec, "n_max", nmax);
      kstat += check_and_set_key_int (orec, "n_param", npar);
      kstat += check_and_set_key_int (orec, "n_ang", ntht);
      kstat += check_and_set_key_int (orec, "n_ang_subsmpl", nthts);
      kstat += check_and_set_key_str (orec, "GuessTable", guessfile);
      kstat += check_and_set_key_int (orec, "GuessRidges", nrdtot);
      kstat += check_and_set_key_float (orec, "Ux_guess", ux_guess);
      kstat += check_and_set_key_float (orec, "Uy_guess", uy_guess);
      kstat += check_and_set_key_int (orec, "Smoothing", smooth);
      kstat += check_and_set_key_int (orec, "m_data", mdata);
      kstat += check_and_set_key_int (orec, "m_max_filt", mfilt);
      kstat += check_and_set_key_float (orec, "FitParmTol", xtol);
      kstat += check_and_set_key_float (orec, "FitFuncTol", ftol);
      kstat += check_and_set_key_float (orec, "FitFuncPrec", feps);
      kstat += check_and_set_key_float (orec, "trust", delta0);
      kstat += check_and_set_key_int (orec, "iter_max", maxeval);
      kstat += check_and_set_key_int (orec, "GradFlag", igrad);
      kstat += check_and_set_key_int (orec, "IntrpFlag", iextrap);
      kstat += check_and_set_key_int (orec, "HessianFlag", idiag);
      if (kstat) {
	fprintf (stderr, "Error writing key value(s) to record %d in series %s\n",
	  rgn, outser);
	fprintf (stderr, "      series may not have appropriate structure\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
     	drms_close_record (orec, DRMS_FREE_RECORD);
        return 1;
      }
    }
					    /*  end shared code with rdfitc  */

    status = 0;
				  /*  reinitialize ioptions(0) with maxeval  */
    ioptions[0] = maxeval;

    ringanalysis_ (spec, &nkx, &nky, &nnu, &ntht, &nthts, &nk, &crpix1, &crpix2, 
        &smooth, &dk, &dnu, &lnbase, &nmin, &nmax, &npar, &ntot, &ux_guess,
        &uy_guess, ux_fit, d_ux, uy_fit, d_uy, amp, d_amp, bg, d_bg, fwhm,
        d_fwhm, nu, d_nu, delnu, nval, kval, kbin, 
        kbstrt, kbend, &kbmin, &kbmax, &nrdg, min_func, n_iter, fit_chk,
	chi, &kmin, &kmax, &modecount, good, unwrapped_pspec, filter,
	filtered_pspec, fnu, width, ampg, bkg, &Verbose, &mdata, &mfilt,
	doptions, ioptions, &iflg, &rgn, &nrdtot, &status);

    drms_free_array (pspec);
    if (status) {
      quality = NOFITS;
      if (status < 10) {
	fprintf (stderr, "Unknown problem: code %d", status);
        quality |= 0x40000000;
      } else if (status < 20) {
	fprintf (stderr, "Problem unwrapping power spectrum: code %d", status);
	fprintf (stderr, "        (set in cart_to_polar subroutine)");
        quality |= 0x20000000;
      } else if (status < 30) {
	fprintf (stderr,
	    "Problem subsampling and filtering spectrum: code %d", status);
	fprintf (stderr, "        (set in fourier_filter subroutine)");
        quality |= 0x10000000;
      } else if (status == 350) {
	fprintf (stderr, "Problem code %d  - no fits for power spectrum\n", status);
	fprintf (stderr, "        %s\n", recid);
        quality |= 0x08000000;
      } else {
	fprintf (stderr, "Problem fitting power spectrum: code %d", status);
	fprintf (stderr, "        (set in hmifit or subsidiary subroutine)");
        quality |= 0x04000000;
      }
      if (drms_output) {
        drms_setkey_int (orec, "Quality", quality);
	for (n = nmin; n < nmax; n++) {
	  snprintf (key, 64, "n_%d_fits", n);
	  drms_setkey_int (orec, key, 0);
	}
        drms_close_record (orec, dispose);
      }
      continue;
    }
    quality = 0;
    if (drms_output) {
      drms_setkey_int (orec, "Quality", quality);
    }

    unit24 = fopen (outfile, "w");
    if (!unit24) {
      fprintf (stderr, "Error: could not open file %s for output\n", outfile);
      drms_close_records (ids, DRMS_FREE_RECORD);
      if (drms_output) drms_close_record (orec, DRMS_FREE_RECORD);
      return 1;
    }
    if (verbose) printf ("rgn = %d\n", rgn);

    lmin = kmin * rsun;
    lmax = kmax * rsun;
      							   /*  write header  */
    fprintf (unit24, "# input record = %s\n", recid);
    if (drms_output) {
      drms_sprint_rec_query (recid, orec);
      fprintf (unit24, "# output record = %s\n", recid);
    } else fprintf (unit24, "# output file = %s\n", outfile);
    fprintf (unit24, "# program used = %s\n", module_ident);

    fprintf (unit24, "# version %s\n", version_id);
    fprintf (unit24, "# guess file used = %s\n", guessfile);
    fprintf (unit24, "# nmin = %2d        nmax = %2d \n", nmin, nmax);
    fprintf (unit24, "# lmin = %8.4f    lmax = %8.4f   kmin = %8.4f   kmax = %8.4f \n", lmin, lmax, kmin, kmax);
    fprintf (unit24, "# delta_nu = %13.4e   delta_k = %13.4e \n", dnu, dk);
    fprintf (unit24, "# number of modes fit = %d\n", modecount);
    fprintf (unit24, "#n     l        k        nu        d_nu ");
    fprintf (unit24, "        ux         d_ux         uy         d_uy    fit    amp");
    fprintf (unit24, "        d_amp         bg         d_bg        fwhm");
    fprintf (unit24, "       d_fwhm      delnu        d_nu    k_bin  nfe");
    fprintf (unit24, "   min_func      rdchi\n");
					      /*  write out individual fits  */
    for (n = nmin; n <= nmax; n++) nvct[n] = 0;
    for (n = 0; n <  ntot*nk; n++) {  
      if (good[n] == 1) {
	ell = kval[n] * rsun;
	fprintf (unit24, "%2d %8.3f %8.5f", nval[n], ell, kval[n]);
	fprintf (unit24, "%10.3f%12.4e", nu[n], d_nu[n]);
	fprintf (unit24, "%12.4e%12.4e", ux_fit[n], d_ux[n]); 
	fprintf (unit24, "%12.4e%12.4e", uy_fit[n], d_uy[n]); 
	fprintf (unit24, "%3d%12.4e%12.4e", fit_chk[n], amp[n], d_amp[n]); 
	fprintf (unit24, "%12.4e%12.4e", bg[n], d_bg[n]); 
	fprintf (unit24, "%12.4e%12.4e", fwhm[n], d_fwhm[n]); 
	fprintf (unit24, "%12.4e%12.4e", delnu[n], d_nu[n]); 
	fprintf (unit24, "%5d%6d%12.4e", kbin[n], n_iter[n], min_func[n]);
	fprintf (unit24, "%12.3e\n", chi[n]);
	nvct[nval[n] - nmin]++;
      }
    }
    fclose (unit24);
    if (drms_output) {
      for (n = nmin; n <= nmax; n++) {
	sprintf (key, "n_%d_fits", n);
	drms_setkey_int (orec, key, nvct[n]);
      }
    }

    if (write_filter) {
      if (!filter_array) {
	int axes[] = {nthts, nk};
	filter_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, filter,
	    &status);
	if (status) {
	  fprintf (stderr, "Error: couldn't create output array for filter\n");
	  return 1;
	}
      }
        oseg = drms_segment_lookupnum (orec, filtsegnum);
        printf ("writing filter array to seg num %d\n", filtsegnum);

	if (drms_segment_write (oseg, filter_array, 0)) {
	  fprintf (stderr, "Error writing data to record in series %s\n",
	      outser);
	  fprintf (stderr,
	      "      series may not have appropriate structure\n");
	  return 1;
      }
    }
    if (write_unwrap_pow) {
        if (!unwrap_array) {
          int axes[] = {ntht, nk, nnu};
	  unwrap_array = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, unwrapped_pspec,
	    &status);
	  if (status) {
	    fprintf (stderr,
	      "Error: couldn't create output array for unwrapped power spectrum\n");
	    return 1;
	  }
	}
        oseg = drms_segment_lookupnum (orec, unwrapsegnum);
        printf ("writing unwrapped power to seg num %d\n", unwrapsegnum);
	if (drms_segment_write (oseg, unwrap_array, 0)) {
	  fprintf (stderr, "Error writing data to record in series %s\n",
	      outser);
	  fprintf (stderr,
	      "      series may not have appropriate structure\n");
	  return 1;
	}
        printf ("succeeded\n");
    }
    if (write_filt_pow) {
      if (!filtps_array) {
        int axes[] = {nthts, nk, nnu};
	filtps_array = drms_array_create (DRMS_TYPE_FLOAT, 3, axes,
	    filtered_pspec, &status);
	if (status) {
	  fprintf (stderr,
	      "Error: couldn't create output array for filtered power spectrum\n");
	  return 1;
	}
      }
      oseg = drms_segment_lookupnum (orec, filtpowsegnum);
      printf ("writing filtered unwrapped power to seg num %d\n", filtpowsegnum);
      if (drms_segment_write (oseg, filtps_array, 0)) {
	fprintf (stderr, "Error writing data to record in series %s\n", outser);
	fprintf (stderr,
	    "      series may not have appropriate structure\n");
	return 1;
      }
      printf ("succeeded\n");
    }

    if (drms_output) drms_close_record (orec, dispose);
    else if (write_filt_pow || write_unwrap_pow || write_filter) {
	fprintf (stderr,
	    "Warning: output of unwrapped spectrum and/or filter requested,\n");
	fprintf (stderr, "  but this is only supported in DRMS\n");
    }
  }
  drms_close_records (ids, DRMS_FREE_RECORD);
  return 0;
  
}

/*
 *  Revision History
 *
 *  09.09.15	created by RSB, based on rdfitc
 *  09.09.25    revised by DAH, put in all parameters necessary for rdfitf
 *                              added in use of crpix to identify center of
 *                              power spectrum
 *  09.11.03	revised by RSB, along with ringtst_v07.f, in order to compile
 *  09.11.03	":	added lnbase to arguments list; removed unwrap_powfile,
 *			filter_file, filtered_powfile
 *  10.01.13	revisions by DAH: added mdata, mfilt, and iflg to module and
 *			ringanalysis_ arguments list; added module arguments 
 *			write_filter, write_filter_pow
 *  10.01.17	revisions by RSB: added keyword propagation;
 *	transferred writing of ancillary files to main module body;
 *	implemented looping on multiple data records (in DRMS only);
 *	replaced use of limitsfile with parameters for k-bin start and stop;
 *	write first few lines of output file header from main body, removing
 *		need for passing multiple string arguments;
 *	input (and output if in DRMS) identifiers in fit file refer to records,
 *		not filenames
 *	write specific input record id to Source key
 *  10.03.10     revised by DAH to pass arrays back to c program for output to
 *               drms
 *  v 0.8	fixed bug in freeing debug arrays
 *    10.05.18	Added fclose when output record skipped
 *    10.06.01	Include new (C) read_guess function
 *		Added output keys BLD_VERS, Created
 *    10.06.03	Added output keys for parameters, results; records written
 *		without segments when no fits
 *  v 0.8 frozen 2010.08.19
 *	11.04.18	Fixed bug in determination of default input segment
 *		when multiple segments possible
 *	19.04.23	Modified treatment of CRPIX1,2 to reflect corrections
 *		in output of module pspec3; added some sanity checks; removed
 *		defaults for input power spectra and guessfile; changed flag 2
 *		to d
 */
