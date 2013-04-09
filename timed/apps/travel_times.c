/*
 *  travel_times.c						 ~rick/hmi/timed
 *
 *  Calculate travel times for selected annuli in a tracked data cube using
 *    Gabor wavelet fitting and the method of Gizon & Birch
 *  The module expects a 3-dimensional real dataset and a pixel locator data
 *    file of a special format
 *
 *  Responsible:  Junwei Zhao & Rick Bogart	      RBogart@solar.Stanford.EDU
 *
 *  Usage:
 *    travel_times [-v] in= pxloc= out= ...
 *
 *
 *  Arguments: (type		default         description)
 *	in	DataSet
 *		hmi_test.tdVtrack_synop[2010.07.08_04:00:00_TAI][][0.0][0.0]
 *			            	Input dataset
 *			A series of records containing 3-d data cubes as one
 *			of their segments. It is assumed that all records in
 *			the dataset are from the same dataseries, or at least
 *			share a common segment structure
 *	pxloc	DataSet		hmi.tdpixlist[]
 *					Annulus pixel locator data records
 *			A set of records specifying the different annuli to
 *			be used for analysis
 *	out	DataSeries	hmi_test.tdVtimes_synop
 *					Output data series
 *      copy	String		"+"	Comma separated list of keys to
 *			propagate forward; if the list includes "+", all
 *			keys in the default list of propagation keys, plus
 *			all prime keys in the input records will be copied
 *
 *  Flags
 *
 *  Notes:
 *    All the work is done by the FORTRAN subroutines mainsub1_() & mainsub2_(),
 *	whose arguments are:
 *	  (float) velo
 *	  (float) wid
 *	  (int) num_annuli
 *	  (float) time[max(num_annuli)]
 *	  (float) coef[5]
 *	  (int) n_start
 *	  (int) n_th
 *	  float *input : input data values
 *	  char *pix_locat_file (int *len_ploc_file)
 *	    name of a FOTRAN unformatted file containing the following vectors
 *	    of length n1 * n2 * num_annuli:
 *		(int) num_w
 *		(int) num_e
 *		(int) num_n
 *		(int) num_s
 *		(short) quad_w_x[num_w+1]
 *		(short) quad_w_y[num_w+1]
 *		(short) quad_e_x[num_e+1]
 *		(short) quad_e_y[num_e+1]
 *		(short) quad_n_x[num_n+1]
 *		(short) quad_n_y[num_n+1]
 *		(short) quad_s_x[num_s+1]
 *		(short) quad_s_y[num_s+1]
 *	  float *out1
 *	    array for the output of the do_fitting_() procedure
 *	  float *out2
 *	    array for the output of the gbtimes-2_() procedure
 *    Subsidiary FORTRAN code is in the following files:
 *	  subroutine	  location		  called in
 *	C_CORRELATE	sub_correlate_BLAS.f	MAINSUB1, MAINSUB2
 *	DO_FITTING	sub_do_fitting.f	MAINSUB1, MAINSUB2
 *	FILTER		sub_filter_HMI_ppline.f	MAINSUB1, MAINSUB2
 *	GBTIMES02	MAINSUB1, MAINSUB2
 *	LMFIT		sub_lmfit.f		DO_FITTING
 *	SAXPY		library			MAINSUB1, MAINSUB2, DO_FITTING
 *	SCOPY		library			MAINSUB1, MAINSUB2, DO_FITTING
 *	SHIFT		library			MAINSUB1, MAINSUB2
 *	SSCAL		library			MAINSUB1, MAINSUB2, DO_FITTING
 *    Typical processing time on a single processor for a 256*256*640 cube
 *	~
 *
 *  Bugs:
 *    Few parameter options; many of the parameters are read from pre-calculated
 *	data in the "pxloc: record. In particular, the input and output array
 *	sizes are hardcoded, but there is no checking of the actual input data.
 *      It is not clear that the module will work with any but a 512*512*640
 *	tracked cube
 *    Scaled output not properly supported; no setting of array parameter
 *    FORTRAN Code writes progress/status to stdout
 *    Missing input data are unconditionally set to 0
 *    The value written into CRLT_Mid is geocentric
 *    No documentation
 *
 *  Revision history is at the end of the file.
 *
 */

#include <jsoc_main.h>
#include "keystuff.c"
#include "earth_ephem.c"
						      /*   module identifier  */
char *module_name = "travel_times_loop_wGB";
char *version_id = "0.9.2";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in",
      "hmi_test.tdVtrack_synop[2010.07.08_04:00:00_TAI][][0.0][0.0]",
      "input record set"}, 
  {ARG_STRING, "pxloc", "hmi.tdpixlist[]", "pixel locater record set"}, 
  {ARG_STRING, "out", "hmi_test.tdVtimes_synop", "output series"}, 
  {ARG_STRING, "copy",  "+",
      "comma separated list of keys to propagate forward"},
  {ARG_FLAG, "v", "0", "verbose output"},      
  {}
};

       /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrRot", "CMLon", "LonHG", "LatHG", "LonCM", "MidTime",
    "CRVAL1", "CRVAL2", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2",
    "T_START", "T_STOP", "LonSpan", "Coverage", "Duration", "MapProj", "MapPA",
    "Zonal Trk", "ZonalVel"};

#define NUM_ANNULI_MAX	(20)

extern void mainsub1_ (float *velo, float *wid, int *num_annuli, float *time,
    float *coef, int *n_start, int *n_th, float *input,
    char *pix_locat_file, int *len_ploc_file, float *out1, float *out2);
extern void mainsub2_ (float *velo, float *wid, int *num_annuli, float *time,
    float *coef, int *n_start, int *n_th, float *input,
    char *pix_locat_file, int *len_ploc_file, float *out1, float *out2);

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *plds, *inds;
  DRMS_Record_t *ploc, *trck, *orec;
  DRMS_Segment_t *pseg, *dseg, *seg_gabor, *seg_gb;
  DRMS_Array_t *data_array, *res_gabor, *res_gb;
  FILE *timesfile, *fout;
  double mapscale_factor;
  float **ttimes;
  float *data, *input, *input2, *input3, *output, *output2, *out1, *out2;
  float *annmin, *annmax, *velo, *wid, *coef0, *coef1, *coef2, *coef3, *coef4;
  float coef[5];
  float newctr, newscale;
  long long n, ntot;
  int *num_annuli, *n_start;
  int trck_axis[2];
  int ann, i, j, k, *axes;
  int len_ploc_file, outok, ploc_version, status;
  int kstat, key_n, propct, propstd;
  int rec, rgnct, smpl, smpls;
  int need_ephem;
  char **pix_locat_file, **copykeylist;
  char *genstr, *errmsg;
  char plocfile[DRMS_MAXQUERYLEN];
  char source[DRMS_MAXQUERYLEN], recid[DRMS_MAXQUERYLEN];
  char module_ident[64], keyname[32];

  int keyct = sizeof (propagate) / sizeof (char *);

  char *plocds = strdup (params_get_str (params, "pxloc"));
  char *trckds = strdup (params_get_str (params, "in"));
  char *outser =  strdup (params_get_str (params, "out"));
  char *propagate_req = strdup (params_get_str (params, "copy"));
  int verbose = params_isflagset (params, "v");

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);
	     /*  read info from correlation pixel-locator records to be used  */
  plds = drms_open_records (drms_env, plocds, &status);
  if (status) {
    fprintf (stderr, "Error: unable to open pixel target record set %s\n",
	plocds);
    return 1;
  }
  smpls = plds->n;
  if (verbose)
    printf ("  calculating travel time maps for %d annulus sets\n", smpls);
  if (smpls != 11)
    fprintf (stderr, "Warning: untested case of %d annulus sets\n", smpls);
  num_annuli = (int *)malloc (smpls * sizeof (int));
  annmin = (float *)malloc (smpls * sizeof (float));
  annmax = (float *)malloc (smpls * sizeof (float));
  velo = (float *)malloc (smpls * sizeof (float));
  wid = (float *)malloc (smpls * sizeof (float));
  n_start = (int *)malloc (smpls * sizeof (int));
  coef0 = (float *)malloc (smpls * sizeof (float));
  coef1 = (float *)malloc (smpls * sizeof (float));
  coef2 = (float *)malloc (smpls * sizeof (float));
  coef3 = (float *)malloc (smpls * sizeof (float));
  coef4 = (float *)malloc (smpls * sizeof (float));
  ttimes = (float **)malloc (smpls * sizeof (float *));
  pix_locat_file = (char **)malloc (smpls * sizeof (char **));
  for (smpl = 0; smpl < smpls; smpl++) {
    ploc = plds->records[smpl];
    drms_sprint_rec_query (source, ploc);
				   /*  check version of pixel locator record  */
    if (smpl) {
      if (ploc_version != drms_getkey_int (ploc, "Version", &status)) {
	if (ploc_version) {
	  fprintf (stderr,
	      "Warning: inconsistent versions of pixel locator records\n");
	  fprintf (stderr,
	      "         Version number will be set to 0 in output\n");
	}
	ploc_version = 0;
      }
    } else {
      ploc_version = drms_getkey_int (ploc, "Version", &status);
      if (status) {
	fprintf (stderr,
	    "Warning: unknown version of pixel locator record(s)\n");
        ploc_version = 0;
      }
    }
				     /*  get additional key info from record  */
    num_annuli[smpl] = drms_getkey_int (ploc, "Annuli", &status);
    if (num_annuli[smpl] > NUM_ANNULI_MAX) {
      fprintf (stderr,
          "Warning: num_annuli exceeds fixed mem alloc in Fortran code\n");
      fprintf (stderr, "         reduced to %d\n", NUM_ANNULI_MAX);
      num_annuli[smpl] = NUM_ANNULI_MAX;
    }
    annmin[smpl] = drms_getkey_float (ploc, "RadInner", &status);
    annmax[smpl] = drms_getkey_float (ploc, "RadOuter", &status);
    velo[smpl] = drms_getkey_float (ploc, "PhaseVel", &status);
    wid[smpl] = drms_getkey_float (ploc, "FiltWid", &status);
    coef0[smpl] = drms_getkey_float (ploc, "AmplInit", &status);
    coef1[smpl] = drms_getkey_float (ploc, "DFreqIni", &status);
    coef2[smpl] = drms_getkey_float (ploc, "GrpTTIni", &status);
    coef3[smpl] = drms_getkey_float (ploc, "FreqInit", &status);
    coef4[smpl] = drms_getkey_float (ploc, "PhsTTIni", &status);
    n_start[smpl] = drms_getkey_float (ploc, "NStart", &status);
    ttimes[smpl] = (float *)calloc (NUM_ANNULI_MAX, sizeof (float));
    pseg = drms_segment_lookup (ploc, "times");
    if (!pseg) {
      fprintf (stderr,
          "Error: record %s\n  does not contain required segment \"times\"\n",
	  source);
      drms_close_records (plds, DRMS_FREE_RECORD);
      return 1;
    }
    drms_segment_filename (pseg, plocfile);
    timesfile = fopen (plocfile, "r");
    for (n = 0; n < num_annuli[smpl]; n++)
      fscanf (timesfile, "%f", &(ttimes[smpl][n]));
		  /*  should not be necessary, calloc should take care of it  */
    for (; n < NUM_ANNULI_MAX; n++) ttimes[smpl][n] = 0.0;
    fclose (timesfile);
							      /*  read times  */
    pseg = drms_segment_lookup (ploc, "pxloc.fortran");
    if (!pseg) {
      fprintf (stderr,
          "Error: record %s\n  does not contain required segment \"pxloc.fortran\"\n",
	  source);
      drms_close_records (plds, DRMS_FREE_RECORD);
      return 1;
    }
    drms_segment_filename (pseg, plocfile);
    pix_locat_file[smpl]= (char *)malloc (strlen (plocfile) + 1);
    strncpy (pix_locat_file[smpl], plocfile, strlen (plocfile));
    pix_locat_file[smpl][strlen (plocfile)] = '\0';
  }
  					       /*  allocate the output array  */
  axes = (int *)malloc (4 * sizeof (int));
  axes[0] = axes[1] = 256;
  axes[2] = 4;
  axes[3] = smpls;
  output = (float *)malloc (smpls * 4 * 256 * 256 * sizeof (float));
  output2 = (float *)malloc (smpls * 4 * 256 * 256 * sizeof (float));
  out1 = (float *)malloc (4 * 256 * 256 *sizeof(float));
  out2 = (float *)malloc (4 * 256 * 256 *sizeof(float));
  res_gabor = drms_array_create (DRMS_TYPE_FLOAT, 4, axes, output, &status);
  res_gb = drms_array_create (DRMS_TYPE_FLOAT, 4, axes, output2, &status);
  drms_close_records (plds, DRMS_FREE_RECORD);

  input2 = (float *)malloc (256 * 256 * 640 * sizeof(float));
  input3 = (float *)malloc (256 * 256 * 640 * sizeof(float));
						     /*  open input data set  */
  inds = drms_open_records (drms_env, trckds, &status);
  if (status) {
    fprintf (stderr, "Error: unable to open input record set %s\n", trckds);
    drms_free_array (res_gabor);
    drms_free_array (res_gb);
    return 1;
  }
  rgnct = inds->n;
  if (verbose)
    printf ("  in each of %d regions\n", rgnct);

  outok = 0;
  for (rec = 0; rec < rgnct; rec++) {
    trck = inds->records[rec];
    dseg = drms_segment_lookup (trck, "V");
    if (!dseg) {
      fprintf (stderr, "Error: no segment \"V\" in record %d; skipped\n", rec);
      continue;
    }
    for (i = 0; i < 2; i++) {
      trck_axis[i] = dseg->axis[i];
      if (i) {
        if (trck_axis[i] - trck_axis[i-1]) {
	  fprintf (stderr,
	      "Error: code as implemented requires square tracked regions\n");
	  drms_free_array (res_gabor);
	  drms_free_array (res_gb);
	  drms_close_records (inds, DRMS_FREE_RECORD);
	  return 1;
	}
      }
    }
    mapscale_factor = 256.0 / (double)trck_axis[0];
						/*  create the output record  */
    orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    if (!orec) {
      fprintf (stderr, "Error: unable to create record %d in series %s\n", rec,
	  outser);
      fprintf (stderr, "       remainder of processing skipped\n");
      drms_close_records (inds, DRMS_FREE_RECORD);
      return 0;
    }
	 /*  check the output data series struct (only needs to be done once  */
    seg_gabor = drms_segment_lookup (orec, "gabor");
    seg_gb = drms_segment_lookup (orec, "GB");
    if (!outok) {
				      /*  check for the appropriate segments  */
      if (!seg_gabor || !seg_gb) {
	fprintf (stderr, "Error: output data series %s\n", outser);
	fprintf (stderr, "       lacks segment \"GB\" and/or \"gabor\"\n");
	drms_close_records (inds, DRMS_FREE_RECORD);
	return 1;

      }
    /*  check structure of segments: in particular naxis[3] must match smpls  */
      outok = 1;
    }
    res_gabor->bscale = res_gb->bscale = 1.0;
    res_gabor->bzero = res_gb->bzero = 0.0;
				     /*  copy designated keywords from input  */
    propct = construct_stringlist (propagate_req, ',', &copykeylist);
    kstat = 0;
    propstd = 0;
					/*  copy specifically requested keys  */
    for (key_n = 0; key_n < propct; key_n++) {
      if (strcmp (copykeylist[key_n], "+"))
        kstat += check_and_copy_key (orec, trck, copykeylist[key_n]);
      else propstd = 1;
    }
    if (propstd) {
						      /*  copy standard keys  */
      kstat += propagate_keys (orec, trck, propagate, keyct);
					   /*  and any additional prime keys  */
      kstat += copy_prime_keys (orec, trck);
    }
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s:\n",
	  outser);
      fprintf (stderr, "      series may not have appropriate structure;\n");
      fprintf (stderr, "      record %d skipped\n", rec);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
    drms_sprint_rec_query (source, trck);
		    /*  to write CRLT_Mid to output, require MidTime as well  */
    need_ephem = (drms_keyword_lookup (orec, "CRLT_Mid", 1) != NULL);
    if (need_ephem)
      need_ephem = (drms_keyword_lookup (orec, "MidTime", 1) != NULL);
    if (need_ephem) {
      double rsun, clat = 0.0 / 0.0, clon, vr, vn, vw;
      TIME tmid =  drms_getkey_time (orec, "MidTime", &status);
      earth_ephemeris (tmid, &rsun, &clat, &clon, &vr, &vn, &vw);
      kstat += check_and_set_key_float (orec, "CRLT_Mid", clat);
    }
							    /*  set WCS keys  */
    kstat += check_and_set_key_int (orec, "WCSAXES", 2);
    newctr = drms_getkey_float (trck, "CRPIX1", &status);
    newctr = mapscale_factor * (newctr - 0.5) + 0.5;
    kstat += check_and_set_key_float (orec, "CRPIX1", newctr);
    newctr = drms_getkey_float (trck, "CRPIX2", &status);
    newctr = mapscale_factor * (newctr - 0.5) + 0.5;
    kstat += check_and_set_key_float (orec, "CRPIX2", newctr);
    newscale = drms_getkey_float (trck, "CDELT1", &status);
    newscale /= mapscale_factor;
    kstat += check_and_set_key_float (orec, "CDELT1", newscale);
    newscale = drms_getkey_float (trck, "CDELT2", &status);
    newscale /= mapscale_factor;
    kstat += check_and_set_key_float (orec, "CDELT2", newscale);
    genstr = drms_getkey_string (trck, "WCSNAME", &status);
    if (!status) {
	      /*  strip off the trailing '/Time' from a tracked cube WCSNAME  */
      strtok (genstr, "/");
      kstat += check_and_set_key_str (orec, "WCSNAME", genstr);
    }
					     /*  set segment descriptor keys  */
    kstat += check_and_set_key_str (orec, "TTDiff_1",
	"mean (outgoing,ingoing)");
    kstat += check_and_set_key_str (orec, "TTDiff_2", "outgoing-ingoing");
    kstat += check_and_set_key_str (orec, "TTDiff_3", "west-east");
    kstat += check_and_set_key_str (orec, "TTDiff_4", "north-south");
						 /*  update the MapScale key  */
    newscale = drms_getkey_float (trck, "MapScale", &status);
    if (!status) {
      newscale /= mapscale_factor;
      kstat += check_and_set_key_float (orec, "MapScale", newscale);
    }
					     /*  set annulus descriptor keys  */
    for (smpl = 0; smpl < smpls; smpl++) {
      sprintf (keyname, "AnnMin%02d", smpl + 1);
      kstat += check_and_set_key_float (orec, keyname, annmin[smpl]);
      sprintf (keyname, "AnnMax%02d", smpl + 1);
      kstat += check_and_set_key_float (orec, keyname, annmax[smpl]);
      sprintf (keyname, "PhsVel%02d", smpl + 1);
      kstat += check_and_set_key_float (orec, keyname, velo[smpl]);
      sprintf (keyname, "FltWid%02d", smpl + 1);
      kstat += check_and_set_key_float (orec, keyname, wid[smpl]);
    }
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s:\n",
	  outser);
      fprintf (stderr, "      series may not have appropriate structure;\n");
      fprintf (stderr, "      record %d skipped\n", rec);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
						     /*  read the input data  */
    data_array = drms_segment_read (dseg, DRMS_TYPE_FLOAT, &status);
    input = (float *)data_array->data;
    if (!input) {
      fprintf (stderr, "Error: unable to read input data record-segment %d\n",
          rec);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
    for (i=0; i<640; i++) 
      for (j=0; j<256; j++)
	for (k=0; k<256; k++)
	  input2[i*256*256 + j*256 + k] =
	      0.25*(input[i*512*512 + (2*j*512) + (2*k)] + 
	      input[i*512*512 + (2*j)*512 + (2*k+1)] +
	      input[i*512*512 + (2*j+1)*512 + (2*k)] + 
	      input[i*512*512 + (2*j+1)*512 + (2*k+1)]);
		    /*  zero out any missing data; Fourier transforms object  */
    ntot = 1;
    for (i = 0; i < data_array->naxis; i++) ntot *= data_array->axis[i];
    for (n = 0; n < ntot; n++) if (isnan (input[n])) input[n] = 0.0;

    for (ann=0; ann<3; ann++) {
      int nth = ann + 1;
      if (verbose) printf ("begin annulus %d\n", ann);
			 /*  string lengths required for C/FORTRAN interface  */
      len_ploc_file = strlen (pix_locat_file[ann]);
				   /*  Call Fortran Code subroutine MAIN_SUB  */
      coef[0] = coef0[ann];
      coef[1] = coef1[ann];
      coef[2] = coef2[ann];
      coef[3] = coef3[ann];
      coef[4] = coef4[ann];
      mainsub1_ (&velo[ann], &wid[ann], &num_annuli[ann], ttimes[ann], coef,
	  &n_start[ann], &nth, input, pix_locat_file[ann],
	  &len_ploc_file, out1, out2);

      for (j=0; j<4*256*256; j++) {
        output[ann*4*256*256+j] = out1[j]; 
        output2[ann*4*256*256+j] = out2[j];
      }
	     /*  need to reread segment because input array may have changed  */
      drms_free_array (data_array);
      data_array = drms_segment_read (dseg, DRMS_TYPE_FLOAT, &status);
      input = (float *)data_array->data;
      if (!input) {
	fprintf (stderr, "Error: unable to read input data record-segment %d\n",
            rec);
	drms_close_record (orec, DRMS_FREE_RECORD);
	ann = 11;
        continue;
      }
      for (n = 0; n < ntot; n++) if (isnan (input[n])) input[n] = 0.0;
    }

    ntot = 256 * 256 * 640;
    for (n = 0; n < ntot; n++) if (isnan (input2[n])) input2[n] = 0.0;
    for (; ann<11; ann++) {
      int nth = ann + 1;
      if (verbose) printf ("begin annulus %d\n", ann);
      for (n = 0; n < ntot; n++) input3[n] = input2[n];
			 /*  string lengths required for C/FORTRAN interface  */
      len_ploc_file = strlen (pix_locat_file[ann]);
      coef[0] = coef0[ann];
      coef[1] = coef1[ann];
      coef[2] = coef2[ann];
      coef[3] = coef3[ann];
      coef[4] = coef4[ann];
      mainsub2_ (&velo[ann], &wid[ann], &num_annuli[ann], ttimes[ann], coef,
	  &n_start[ann], &nth, input3, pix_locat_file[ann],
	  &len_ploc_file, out1, out2);
      for (j=0; j<4*256*256; j++) {
        output[ann*4*256*256 + j] = out1[j]; 
        output2[ann*4*256*256 + j] = out2[j];
      }
    }
    status = drms_segment_write (seg_gabor, res_gabor, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      drms_segment_write returned %d\n", status);
      fprintf (stderr, "      series may not have appropriate structure\n");
      fprintf (stderr, "      record %d skipped\n", rec);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
    status = drms_segment_write (seg_gb, res_gb, 0);
    if (status) {

      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      fprintf (stderr, "      record %d skipped\n", rec);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
						    /*  set a few final keys  */
    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_str (orec, "BLD_VERS", jsoc_version);
    kstat += check_and_set_key_str (orec, "Source", source);
    kstat += check_and_set_key_str (orec, "Input", trckds);
    kstat += check_and_set_key_str (orec, "PixLocate", plocds);
    kstat += check_and_set_key_int (orec, "PxLocVer", ploc_version);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s:\n",
	  outser);
      fprintf (stderr, "      series may not have appropriate structure;\n");
      fprintf (stderr, "      record %d skipped\n", rec);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }

    drms_close_record (orec, DRMS_INSERT_RECORD);
  }
  drms_close_records (inds, DRMS_FREE_RECORD);
  return 0;
}
/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  10.09.24	created this file as prototype
 *  10.10.08	copied from working and tested version file
 *	/tmp00/junwei/HMI_script/travel_time/travel_times_loops_wGB.c
 *		added comments; minor mod to dup returned string pointers;
 *		incorporated get_data function in file from original file
 *	/tmp00/junwei/HMI_script/travel_time/get_data.c
 *  v 0.5 frozen 2010.10.11
 *  10.10.11	removed extraneous (commented) code
 *		replaced hard-coded parameter file names with hard-coded
 *	parameter values
 *  10.11.15	added verbose argument and some diagnostic output; added
 *	argument for and processing of DRMS data records for correlation pixel
 *	sets; input read from DRMS data records, with option of processing
 * 	multiple records
 *  v 0.6 frozen 2011.01.06
 *  11.01.28	stripped old test code, debugging lines; changed output from
 *	named files to DRMS series; added some error checking on input and
 *	output; removed get_data, need for cfitsio; fixed possible indexing
 *	error; removed included fftw3 source from sub_filter_HMI_ppline.f;
 *	changed output to use drms utilities rather than private fits code;
 *	added several standard output keys; incorporated updated version of
 *	sub_GB_fitting_2002.f; removed debug prints from fortran routines
 *  11.02.25	fixed bug in setting of filename for fortran
 *  11.03.07	free data_array before rereading; added debug statements
 *  11.05.27	zero out any missing input data
 *  v 0.7 frozen 2011.06.08
 *  11.06.10	added propagation of version number of pixel locator records
 *  11.07.18	added setting of CRLT_Mid; added WCS keys plus numerous others
 *	to default propagation list
 *  v 0.8 frozen 2011.07,18
 *  11.12.03	removed setting of MapScale key, pending availability from
 *	pxloc data
 *  12.01.13	fixed setting of MapScale key; added setting of annuli info
 *	keys
 *  v 0.9 frozen 2012.01.14
 *  v 0.91 frozen 2012.08.04
 *  12.09.28	added option of adding to or overriding standard list of
 *	keywords from input to output
 *  13.03.12	corrected setting of  CRPIX1,2 and CDELT1,2 to values consistent
 *	with pixel locator file (supposedly); added MapProj and MapPA to the
 *	default list of propagated keywords; added modified propagation of
 *	WCSNAME value
 */
