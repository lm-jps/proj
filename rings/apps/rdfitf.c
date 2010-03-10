/*
 *  rdfitf.c					~rick/src/ringfit/rdfitf
 *
 *  The "fast" ring-fitting procedure based on fitting only a few parameters
 *    to the low-order azimuthal Fourier components in k-space of the power
 *    spectra, rather than all m values
 *
 *  Usage:
 *    rdfitf [-nv] [param= ...]
 *
 *  Parameters: 	(type   default         description)
 *	in	DataSet	-	Input dataset
 *	out	string	fit.out	output data series or (if series is not
 *				writeable) file name
 *	guessfile string	~dhaber/ringtst/ptsa2cn44555newline.dat
 *				name of file containing parameterizations for
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
 *	nmax	int	6	highest radial order to fit
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
 *	kstart	ints	{12, 10,  8,  8,  8,  8,  9,  9}
 *				starting values for k-bins for each n
 *	kstop	ints	{50, 50, 48, 44, 37, 30, 23, 20}
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
 *    Considerable I/O is handled from the Fortran code rather than the
 *	driver, meaning that the DRMS is not used very effectively; in
 *	particular, guessfile information should be in DRMS, and read from
 *	module driver
 *    Several of the named value arguments are really just flags or enums
 *    None of the run parameters are put into the DRMS output record
 *    Writing of intermediate data files (unwrapped spectra etc.) is only
 *	supported in DRMS, there is no direct FITS write to named file
 *    The test on lnbase in cart_to_polar is probably wrong, since lnbase
 *	will be returned as NaN rather than 0 if the keyword is missing
 *    Should use WCS keywords if present
 *    It is not clear what happens if npar > 6; module fails for npar < 6.
 *    If the extra verbose flag is set, a file named calc_guess.dat is created
 *	in the current working directory (in read_guess_v3.f90)
 *    Most or all of the internally written information when the verbose flag
 *	is set should be written to the Log
 *    An RSunRef value is set as a constant in hmifits; it would be better to
 *	use the keyword value associated with the input spectrum if available
 *    The default values of kstart and kstop refer to the nunmber of kbins, and
 *	need to be adjusted for different sizes of tiles; they should be
 *	calculated not read in as parameters (or from a file as formerly)
 *    The decision as to whether the output specifier refers to a data series
 *	or a file name based on whether it can be referred to a series
 *	writeable by the user is not sensible, but probably not particularly
 *	dangerous.
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
 *	outfile	[*]	name of the output data series or file
 *      guessfile [*]	name of the file containing guess frequencies,
 *		amplitudes, widths, and background parameters
 *	unwrapped_pspec	float* of unwrapped power spectrum (returned)
 *	filter		float* of wavenumber filter (returned)
 *	filtered_pspec	float* of filtered power spectrum (returned)
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

char *module_name = "ringfit_ssw";
char *module_desc = "ring fitting using method of Haber and Hindman";
char *version_id = "0.8";

#include <jsoc_main.h>

ModuleArgs_t module_args[] = {
  {ARG_DATASET,	"in", "mdi.rdVpspec_dp[1988][180][180][0]", "input dataset"}, 
  {ARG_STRING, "out", "fit.out",
      "name of output data series (or file) containing fit coefficients"},
  {ARG_STRING, "guessfile", "/home/dhaber/ringtst/ptsa2cn44555newline.dat",
      "name of file containing parameterizations for each fitting parameter"},
  {ARG_STRING, "unwrap", "unwrap",
      "segment containing optional unwrapped power spectra"},
  {ARG_STRING, "filter", "filter",
      "segment containing optional spatial wavenumber filter"},
  {ARG_STRING, "filtps", "filtps",
      "segment containing filtered, subsampled, power spectra"},
  {ARG_INT, "nmin", "0", "lowest radial order to fit", "[0,)"},
  {ARG_INT, "nmax", "6", "highest radial order to fit", "[0,)"},
  {ARG_INT, "npar", "6", "number of parameters to fit", "[6,)"},
  {ARG_INT, "ntht", "400",
      "number of azimuth values for initial unwrapping of spectrum, usually = 2*pi*nk/2, should be divisible by 4"},
  {ARG_INT, "nthts", "40",
      "number of azimuth values after subsampling, generally a factor of 10 less than ntht, should also be divisible by 4"},
  {ARG_FLOAT, "ux", "0.0", "initial guess for zonal flow speed [m/s]"},
  {ARG_FLOAT, "uy", "0.0", "initial guess for merridional flow speed [m/s]"},
  {ARG_INT, "smooth", "0", "amount of smoothing wanted [wavenumber bins]"},
  {ARG_INT, "mdata", "16",
      "number of azimuthal Fourier components to keep when subsampling the data"},
  {ARG_INT, "mfilt", "8",
      "highest azimuthal Fourier component to use in filtering the data"},
  {ARG_INTS, "kstart", "{12, 10,  8,  8,  8,  8,  9,  9}",
      "starting kbin values"},
  {ARG_INTS, "kstop",  "{50, 50, 48, 44, 37, 30, 23, 20}",
      "ending kbin values"},
  {ARG_FLAG,  "f", "0", "write spatial wavenumber filter"},      
  {ARG_FLAG,  "n", "0", "do not write out fits (diagnostics only)"},      
  {ARG_FLAG,  "U", "0", "write unwrapped power spectrum"},      
  {ARG_FLAG,  "u", "0", "write filtered subsampled unwrapped power spectrum"},      
  {ARG_FLAG,  "v", "0", "verbose output"}, 
  {ARG_FLAG,  "x", "0", "really verbose output"}, 
  {ARG_FLAG,  "2", "0", "do not use second derivatives in parameter errors"},      
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
    double *dk, double *dnu, float *lnbase, int *nmin, int *nmax, int *npar,
    double *ux, double *uy, char *outfile, int *strln1,
    char *guessfile, int *strln2, int *kbstrt, int *kbend, int *nrdg,
    char *header, int *hdrlen, float *unwrapped_pspec, float *filter,
    float *filtered_pspec, int *verbose, int *mdata, int *mfilt, 
    double *doptions, int *ioptions, int *iflg, int *status);

       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrTime", "CarrRot", "CMLon", "LonHG", "LatHG", "LonCM",
    "MidTime", "Duration", "Cadence", "LonSpan",
    "T_START", "T_STOP", "Coverage", "Size", "Width", "Height",
    "ZonalTrk", "ZonalVel", "MeridTrk", "MeridVel",
    "MapScale", "MapProj", "Map_PA", "RSunRef", "PosAng", "MAI",
    "APODIZNG", "APOD_MIN", "APOD_MAX"};

#include "keystuff.c"

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *pspec, *array;
  float *filtered_pspec;
  double doptions[4]; 
  double dk, dnu;
  float *spec;
  float *unwrapped_pspec, *filter;
  float lnbase, crpix1, crpix2;
  int *kbstrt, *kbend;
  int ioptions[5];
  int hdrlen, strln_out, strln_guess;
  int rgn, rgnct, segct, drms_output, dispose;
  int isegnum, osegnum, logsegnum, filtsegnum;
  int unwrapsegnum, filtpowsegnum;
  int n, nkx, nky, nnu, nk;
  int nrdg;
  int status;
  char logfile[DRMS_MAXPATHLEN], outfile[DRMS_MAXPATHLEN];
  char recid[DRMS_MAXQUERYLEN];
  char header[8192], line[256], module_ident[64], key[64];

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
    printf ("%s version %s\n", module_name, version_id);
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
  }
  if (nmin > nmax) {
    fprintf (stderr, "Error: nmax (%d) < nmin (%d)\n", nmax, nmin);
    return 0;
  }
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
  orec = drms_create_record (drms_env, outser, DRMS_TRANSIENT, &status);
  if (orec) {
    segct = drms_record_numsegments (orec);
    if (segct < 1) {
      fprintf (stderr,
	  "Error: no data segment in output series %s\n", outser);
      drms_close_records (ids, DRMS_FREE_RECORD);
      drms_close_record (orec, DRMS_FREE_RECORD);
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
    drms_close_record (orec, DRMS_FREE_RECORD);
    drms_output = 1;
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
	/*  loop over all input spectra, creating an output record for each  */
  for (rgn = 0; rgn < rgnct; rgn++) {
    irec = ids->records[rgn];
    drms_sprint_rec_query (recid, irec);
    if (verbose && (rgnct > 1)) {
      printf ("processing record %d : %s\n", rgn, recid);
    }
    dk = drms_getkey_double (irec, "DELTA_K", &status);
    dnu = drms_getkey_double (irec, "DELTA_NU", &status);
    lnbase = drms_getkey_float (irec, "LOG_BASE", &status);
    crpix1 = drms_getkey_float (irec, "CRPIX1", &status);
    crpix2 = drms_getkey_float (irec, "CRPIX2", &status);
    segct = irec->segments.num_total;
    if (segct != 1) {
      for (n = 0; n < segct; n++) {
	iseg = drms_segment_lookupnum (irec, n);
	if (iseg->info->naxis == 3) break;
	isegnum = n;
      }
      if (n >= segct) {
	fprintf (stderr, "found no segmemt of rank 3 in input dataset\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
      }
    } else {
      isegnum = 0;
      iseg = drms_segment_lookupnum (irec, isegnum);
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
    nk = nkx/2;

    unwrapped_pspec = (float *)malloc (ntht * nk * nnu * sizeof (float));
    filter = (float *)malloc (nthts * nk * sizeof (float));
    filtered_pspec = (float *)malloc (nthts * nk * nnu * sizeof (float));

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
      kstat += check_and_set_key_str (orec, "Source", recid);
      if (kstat) {
	fprintf (stderr, "Error writing key value(s) to record %d in series %s\n",
	  rgn, outser);
	fprintf (stderr, "      series may not have appropriate structure\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
     	drms_close_record (orec, DRMS_FREE_RECORD);
        return 1;
      }
					    /*  end shared code with rdfitc  */
    }
    sprintf (header, "# input record = %s\n", recid);
    if (drms_output) {
      drms_sprint_rec_query (recid, orec);
      sprintf (line, "# output record = %s\n", recid);
    } else sprintf (line, "# output file = %s\n", outfile);
    strcat (header, line);
    sprintf (line, "# program used = %s\n", module_ident);
    strcat (header, line);
    sprintf (line, "# guess file used = %s", guessfile);
    strcat (header, line);
    hdrlen = strlen (header) + 1;
    strln_out = strlen (outfile) + 1;
    strln_guess = strlen (guessfile) + 1;

    status = 0;
    ringanalysis_ (spec, &nkx, &nky, &nnu, &ntht, &nthts, &nk, &crpix1, &crpix2, 
        &smooth, &dk, &dnu, &lnbase, &nmin, &nmax, &npar, &ux_guess, &uy_guess,
	outfile, &strln_out, guessfile, &strln_guess,
	kbstrt, kbend, &nrdg, header, &hdrlen,
	unwrapped_pspec, filter, filtered_pspec,
	&Verbose, &mdata, &mfilt, doptions, ioptions, &iflg, &status);
    drms_free_array (pspec);
    if (status) {
      if (status < 10)
	fprintf (stderr, "Unknown problem: code %d", status);
      else if (status < 20) {
	fprintf (stderr, "Problem unwrapping power spectrum: code %d", status);
	fprintf (stderr, "        (set in cart_to_polar subroutine)");
      } else if (status < 30) {
	fprintf (stderr,
	    "Problem subsampling and filtering spectrum: code %d", status);
	fprintf (stderr, "        (set in fourier_filter subroutine)");
      } else {
	fprintf (stderr, "Problem fitting power spectrum: code %d", status);
	fprintf (stderr, "        (set in hmifilt or subsidiary subroutine)");
      }
      if (drms_output) drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
    if (drms_output) {
      if (write_filter) {
        int axes[2];
	axes[0] = nthts;
	axes[1] = nk;
        oseg = drms_segment_lookupnum (orec, filtsegnum);
	array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, unwrapped_pspec,
	    &status);
	if (status) {
	  fprintf (stderr,
	    "Error: couldn't create output array for unwrapped power spectrum\n");
	  return 1;
	}
	if (drms_segment_write (oseg, array, 0)) {
	  fprintf (stderr, "Error writing data to record in series %s\n",
	      outser);
	  fprintf (stderr,
	      "      series may not have appropriate structure\n");
	  return 1;
	}
	drms_free_array (array);
      }
      if (write_unwrap_pow) {
        int axes[3];
	axes[0] = ntht;
	axes[1] = nk;
	axes[2] = nnu;
        oseg = drms_segment_lookupnum (orec, unwrapsegnum);
	array = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, unwrapped_pspec,
	    &status);
	if (status) {
	  fprintf (stderr,
	    "Error: couldn't create output array for unwrapped power spectrum\n");
	  return 1;
	}
	if (drms_segment_write (oseg, array, 0)) {
	  fprintf (stderr, "Error writing data to record in series %s\n",
	      outser);
	  fprintf (stderr,
	      "      series may not have appropriate structure\n");
	  return 1;
	}
	drms_free_array (array);
      }
      if (write_filt_pow) {
        int axes[3];
	axes[0] = nthts;
	axes[1] = nk;
	axes[2] = nnu;
        oseg = drms_segment_lookupnum (orec, filtpowsegnum);
	array = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, filtered_pspec,
	    &status);
	if (status) {
	  fprintf (stderr,
	    "Error: couldn't create output array for unwrapped power spectrum\n");
	  return 1;
	}
	if (drms_segment_write (oseg, array, 0)) {
	  fprintf (stderr, "Error writing data to record in series %s\n",
	      outser);
	  fprintf (stderr,
	      "      series may not have appropriate structure\n");
	  return 1;
	}
	drms_free_array (array);
      }
      drms_close_record (orec, dispose);
    } else {
      if (write_filt_pow || write_unwrap_pow || write_filter) {
	fprintf (stderr,
	    "Warning: output of unwrapped spectrum and/or filter requested,\n");
	fprintf (stderr, "  but this is only supported in DRMS\n");
      }
    }
  }

  drms_close_records (ids, DRMS_FREE_RECORD);
  return (0);
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
 */
