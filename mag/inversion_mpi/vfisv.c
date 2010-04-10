/* ------------------------------------------------------------------------------------------------------
 *
 *  vfisv.c                                        ~rick/hmi/vecb/src/v09
 * 3/25/10  this is Keiji's less memory usage version with my edits for keywords * and k=0.
 * 3/23/10 THIS IS a copy of REBECCA's /v13  except with vfisv_mpi
 * this is the mpi version  vfisv.c version 10. It should have the keywords correct AND 12 err files.
 * This directory contains Rebeccas new/improved Invert code with the "scaling  factor"??
 *  C DRMS module version of Juan Borrero's driver program dmrs.f90
 *    version v2.0 modified to handle input file(s) from the DRMS containing
 *    artificial HMI data and work on the mpi.
 * Note the input parameters for the invert module have changed!!
 *
 *  Usage: vfisv [...]
 *    or (on new cluster):  mpiexec -n # vfisv [...]
 *    or (on old cluster):  ~kehcheng/mpich2/bin/mpirun -np # vfisv_mpc [...]
 *
 *  Bugs:
 *    Facilities for inverting only part of an image are for debugging,
 *        should be eliminated from final version, as they add unnecessray
 *        complexity
 *    The statistics for inverted quantities are calculated over the whole
 *        image, not just the part inverted
 *    Off-limb pixels are not handled well, the inversion results are set
 *        to 0; this means that the statistics of the inverted parameters
 *        calculated over the whole image are meaningless.
 *    The initialization of the inversion quantities is silly, an arbitrary
 *        one of the Stokes parameters; they should be calloc'd to 0 or (better)
 *        NaN
 *    1-21-10:Added a map_param which will track which of the datapts were NaN's
 *      and were then converted to Zero. We should not borther to waste computational
 *      power to invert them.
 *    There is no parameter to govern which of the quantities are inverted;
 *        it has to be modified by changing the FREE array in map_param.f90
 *        and recompiling, and output files are generated for all quantities,
 *        regardless of whether they have been calculated or not
 *    Likewise, the QUICKLOOK, STRAYLIGHT_CALCULATION, and ERRORS parameters
 *        are set in code (map_param.f90)
 *    The output FITS files are floats, not scaled shorts; an entire set of
 *        18*imgpix values have to be allocated for assembly; it might be better
 *        to write slices into the output FITS segments.
 *    Can't handle case in which some data segments (e.g. a wavlength set)
 *        are missing.
 *    Number of threads must be a divisor of number of pixels in image;
 *        2**n for full HMI images; need to generalize, and allow for unequal
 *        thread segments
 *
 * --- some comments by K.H. -----
 *
 *  I) March 26, 2010
 *     First, combining efforts done by 5:00PM, March 26, 2010
 *      (1) a lot of things such as JSOC/DRMS things done by Priya and Rick
 *            /home/priya/vecb/VFISV_new/v11/vfisv.c
 *      (2) Default parameters defined by Rebecca
 *            /home/rce/v16/vfisv.c
 *      (3) MPI things by Priya, Rick and Keiji
 *  II) March 29
 *      (1) added lines to handle pixels containing NaN or all-zero values in input
 * III) April  2
 *      (1) added if-blocks at MPI part to avoid deadlock when run with -n 1 or without specifying num. of process
 *      (2) modified DRMS part : correcting JSOC-keyword(s) including T_REC.
 *      (3) moved lines for inversion initialization to be done
 *            after malloc/calloc-ings all wrapper-arrays (and just before the inverson process)
 *      (4) modified positions and contents of printf().
 *      (5) modified how the argument at the command line will be taken into.
 *
 * ------------------------------------------------------------------------------------------------------ */

/* To speed up, only some rectangle regions or a limited-width column be processed. Set 0 to disable */
#define QUICKRUN 0

/* Only non-NAN or physically meaningful pixel will be evenly assigned to each PE. Recommend 1. */
#define EQUIAREA 1

/* 1 if the inversion array differ from those by April 5 */
#define NEWARRAY 1

/* Phil's macro */
#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

char *module_name = "vfisv";
char *version_id = "0.10";

extern void invert_ (double *, double *, double *, double *, double *);
extern void filt_init_ (int *, int *, int *, double *, double *, double *);
extern void free_init_ (int *);
extern void free_memory_ ();
extern void inv_init_ (int *, double *, double *, double *, double *, double *);
extern void line_init_ (double *, double *, double *);
extern void svd_init_ ();
extern void voigt_init_ ();
extern void wave_init_ (double *, double *, int *);

#include <jsoc_main.h>
#include <math.h>
#include <mpi.h>

#define ONECOL (spct*wlct)
#define ONEROW (spct*wlct*cols)
#define IMGCTR (cols*(rows/2) + cols/2)

ModuleArgs_t module_args[] = {
/*This is the input data in the form of 24 fits file segnments */
#if 0
  {ARG_STRING,  "in",  "su_couvidat.hmiStokes[2007.10.14_23:50:00_TAI]", "input series"},
  {ARG_STRING,  "out", "hao_rce.KD4096_Mar26", "output series"},
#else
  {ARG_STRING,  "in", "su_couvidat.HMISeriesLev1pa90[2010.03.27_06:00:00_TAI]", "input series"},  /*This is the input data in the form of 24 fits file segnments */
  {ARG_STRING,  "out", "hao_rce.juanma5", "output series"},
#endif
  {ARG_INT,  "npix", "0", "number of pixels per segment (default: all/segs)"},
  {ARG_INT,  "num_iter", "30", "number of iterations(default: 30)"},
  {ARG_INT,  "num_lambda", "49", "number of ??(default: 33)"},
  {ARG_INT,  "Num_lambda_filter", "6", "Number of filters accross the wvl (default: 6)"},
  {ARG_INT,  "Num_tunning", "6", "Number of ??(default: 6)"},
  {ARG_DOUBLE,  "svd_tolerance", "1.0e-32", "svd tolerance (default: 1.0e-32)"},
  {ARG_DOUBLE,  "chi2_stop", "1.0e-6", "chisq-stop (default: 1.0e-6)"},
  {ARG_DOUBLE,  "Polarization_threshold", "1.0e-2", "polarization threshold (default: 0.01)"},
  {ARG_DOUBLE,  "Percentage_Jump", "10.0", "Percentage Jump (default: 10%)"},
  {ARG_DOUBLE,  "Lambda_Min", "-648.0", "Intensity threshold (default: -432)"},
  {ARG_DOUBLE,  "Lambda_0", "6173.3433", "Wavelength(default:6173.3433 Angstrom )"},
  {ARG_DOUBLE,  "Lambda_B", "0.044475775", "FWHM?? (default: 0.044475775)"},
  {ARG_DOUBLE,  "Delta_Lambda", "27.0", "Delta Lambda(default: 27.0)"},
  {ARG_DOUBLE,  "Lyotfwhm", "424.0", "Lyot filter FWHM (default: 424.0)"},
  {ARG_DOUBLE,  "Wnarrow", "172.0", "W narrow (default: 172.0)"},
  {ARG_DOUBLE,  "Wspacing", "69.0", "W narrow (default: 69.0)"},
  {ARG_FLAG, "v",    "", "run verbose"},
  {ARG_FLAG, "d",    "", "turn damping off"},
#if 0
  {ARG_DOUBLE,  "Intensity_Threshold", "1e2", "Intensity threshold (default: 0.8)"},
  {ARG_DOUBLE,  "Noise_LEVEL", "5.0e3", "Intensity threshold (default: 3.0e-3)"},
  {ARG_DOUBLE,  "Continuum", "0.0", "Intensity threshold (default: 0)"},
#else
  {ARG_DOUBLE,  "Intensity_Threshold", "5e2", "Intensity threshold (default: 0.8)"},
  {ARG_DOUBLE,  "Noise_LEVEL", "0.8e1", "Intensity threshold (default: 3.0e-3)"},
  {ARG_INT,     "Continuum", "0", "Intensity threshold (default: 0)"},
#endif
  {}
};

int DoIt (void)
{
/* JSOC-DRMS variables */
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *records, *out;
  DRMS_Record_t *inRec, *outRec;
  DRMS_Segment_t *seg;
  DRMS_Array_t *stokes_array, *invrt_array,*err_array;

/* get values at commandline argument, as constant */
  const int    npixc              = params_get_int(params, "npix");
  const int    NUM_ITERATIONSc    = params_get_int(params, "num_iter");
  const int    NUM_LAMBDAc        = params_get_int(params, "num_lambda");
  const int    NUM_LAMBDA_FILTERc = params_get_int(params, "Num_lambda_filter");
  const int    NUM_TUNNINGc       = params_get_int(params, "Num_tunning");
  const int    CONTINUUMc         = params_get_int(params, "Continuum");
  const double SVD_TOLERANCEc          = params_get_double(params, "svd_tolerance");
  const double CHI2_STOPc              = params_get_double(params, "chi2_stop");
  const double POLARIZATION_THRESHOLDc = params_get_double(params, "Polarization_threshold");
  const double INTENSITY_THRESHOLDc    = params_get_double(params, "Intensity_Threshold");
  const double PERCENTAGE_JUMPc        = params_get_double(params, "Percentage_Jump");
  const double LAMBDA_MINc             = params_get_double(params, "Lambda_Min");
  const double LAMBDA_0c               = params_get_double(params, "Lambda_0");
  const double LAMBDA_Bc               = params_get_double(params, "Lambda_B");
  const double DELTA_LAMBDAc           = params_get_double(params, "Delta_Lambda");
  const double NOISE_LEVELc            = params_get_double(params, "Noise_LEVEL");
  const double LYOTFWHMc               = params_get_double(params, "Lyotfwhm");
  const double WNARROWc                = params_get_double(params, "Wnarrow");
  const double WSPACINGc               = params_get_double(params, "Wspacing");
  const char   *indsdescc = params_get_str(params, "in");
  const char   *outserc   = params_get_str(params, "out");
  const int    verbosec   = params_isflagset(params, "v");
/* then copy it to non-constants, to avoid JSOC-compiler's complain */
  int    npix;
  int    NUM_ITERATIONS, NUM_LAMBDA, NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM;
  double SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD, INTENSITY_THRESHOLD, PERCENTAGE_JUMP;
  double LAMBDA_0, LAMBDA_B, NOISE_LEVEL;
  double LAMBDA_MIN, DELTA_LAMBDA;
  double LYOTFWHM, WNARROW, WSPACING;
  char   *indsdesc, *outser;
  int    verbose;
  npix = npixc;
  NUM_ITERATIONS    = NUM_ITERATIONSc;
  NUM_LAMBDA        = NUM_LAMBDAc;
  NUM_LAMBDA_FILTER = NUM_LAMBDA_FILTERc;
  NUM_TUNNING       = NUM_TUNNINGc;
  CONTINUUM         = CONTINUUMc;
  SVD_TOLERANCE          = SVD_TOLERANCEc;
  CHI2_STOP              = CHI2_STOPc;
  POLARIZATION_THRESHOLD = POLARIZATION_THRESHOLDc;
  INTENSITY_THRESHOLD    = INTENSITY_THRESHOLDc;
  PERCENTAGE_JUMP        = PERCENTAGE_JUMPc;
  LAMBDA_0     = LAMBDA_0c;
  LAMBDA_B     = LAMBDA_Bc;
  NOISE_LEVEL  = NOISE_LEVELc;
  LAMBDA_MIN   = LAMBDA_MINc;
  DELTA_LAMBDA = DELTA_LAMBDAc;
  LYOTFWHM = LYOTFWHMc;
  WNARROW  = WNARROWc;
  WSPACING = WSPACINGc;
  verbose  = verbosec;
  indsdesc = strdup(indsdescc);
  outser   = strdup(outserc);

/* important variables for inversions */

  int list_free_params[10]={1,1,1,0,1,1,1,1,1,0};
#if 0
  double guess[10]= {15.0,90.0,45.0,0.5,40.0,50.0,0.0,0.4*1e6,0.6*1e6,1.0};
#else
  double guess[10]= {15.0,90.0,45.0,0.5,40.0,50.0,0.0,0.4*6e3,0.6*6e3,1.0};
#endif

/* inversion-related variables used by wrapper */
#if NEWARRAY == 1
  char *Resname[] = {"eta_0", "inclination", "azimuth", "damping", "dop_width", "field",
                     "vlos_mag", "src_continuum", "src_grad", "alpha_mag","chisq",
		     "field_err","inclination_err","azimuth_err", "vlos_err", "alpha_err",
		     "field_inclination_err","field_az_err","inclin_azimuth_err", "field_alpha_err", 
		     "inclination_alpha_err", "azimuth_alpha_err"};
#else
  char *Resname[] = {"eta_0", "inclination", "azimuth", "damping", "dop_width", "field",
                     "vlos_mag", "src_continuum", "src_grad", "alpha_mag",
		     "field_err","inclination_err","azimuth_err", "vlos_err", "alpha_err",
                     "field_inclination_err","field_az_err","inclin_azimuth_err", "field_alpha_err",
		     "inclination_alpha_err", "azimuth_alpha_err"};
#endif
  char segname[16];
  char *spname[] = {"I", "Q", "U", "V"};
  int spct = (sizeof (spname) / sizeof (char *));
  int wlct = 6;
  int paramct = 10;
#if NEWARRAY == 1
  int Err_ct=12;
#else
  int Err_ct=11;
#endif

/* working array etc. used by wrapper */
  double *FinalRes,*FinalErr;
  time_t startime, endtime;
  double *ddat, *res;
  double *obs, *scat,*err;
  float *data, *data0; //, *inv = NULL;
  int *nan_map;
  float *stokesi, *stokesq, *stokesu, *stokesv;
  float fval;
  int param[3];
  int rn, rec_ct, sn, seg_ct;
  int cols, rows, imgpix, imgbytes, imgloc;
  int sp, wl, nvar;
  int imgseg, segsize, ressize;
  int j,m,i,k,n ;
  int status;
/* MPI variables */
  MPI_Status mpistat;
  int mpi_tag, mpi_rank, mpi_size;
  int myrank, nprocs, numpix;
  int istart, iend;
  int *istartall, *iendall;
  void para_range(int,int,int,int *,int *);

/* time at the start */
  time(&startime);

/* Initalizing MPI */
  MPI_Status mpi_stat;
  status = 0;
  MPI_Init (&status, NULL);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);

  istartall = (int *)malloc(sizeof(int) * mpi_size);
  iendall   = (int *)malloc(sizeof(int) * mpi_size);

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {

    printf ("%s Ver. %s\n", module_name, version_id);
    if (verbose) printf ("%d CPU/core(s) be used for parallel-inversion. \n", mpi_size);
    printf("Lambda_O= %f\n",LAMBDA_0);

    records = drms_open_records (drms_env, indsdesc, &status); /*  open input record_set  */
    if (status) {DIE("drms_open_records failed.\n");}
    if ((rec_ct = records->n) == 0){DIE("No records in selected dataset.\n");}
    if ((rec_ct = records->n) >  1){fprintf (stderr, "Warning: only first record in selected set processed\n");}
    rn = 0;

    inRec = records->records[rn];
    seg = drms_segment_lookupnum (inRec, 0);
    cols = seg->axis[0];
    rows = seg->axis[1];
    imgpix = cols * rows;
    imgbytes = imgpix * sizeof (float);
    nvar = wlct * spct;

    if (npix < 1) npix = imgpix / mpi_size;
    if (imgpix % npix) {printf ("Error: npix (%d) must be a divisor of %d\n", npix, imgpix);return 0;}

    data = data0 = (float *)malloc (sizeof (float) * imgpix * nvar);
    nan_map=(int*)calloc (imgpix, sizeof (int));

    for (sp = 0; sp < spct; sp++) /* spct=4,wlct=6 */
    {
      for (wl = 0; wl < wlct; wl++)
      {
        sprintf (segname, "%s%d", spname[sp], wl);
        if ((seg = drms_segment_lookup (inRec, segname)) == NULL){
          fprintf (stderr, "Error reading segment %s of record %d\n", segname, rn);
          return 1;
        }
        /* 4 x 6 segment, 4k x 4k data points each */
        stokes_array = drms_segment_read (seg, DRMS_TYPE_FLOAT, &status);
        /* printf("segment read %s\n",segname); */
        memcpy (data, stokes_array->data, imgbytes);
        drms_free_array (stokes_array);
        data += imgpix;
      }
    }
    data = data0;
    printf("Imgpix= %d\n",imgpix);

/* Calculate "square" of solar radius in float pixel */
    float crpixx, crpixy;
    float cdeltx, cdelty;
    float rsun_ref, dsun_obs;
    float sunarc;
    crpixx   = drms_getkey_float(inRec,"CRPIX1",&status);   // center of the solar disk
    crpixy   = drms_getkey_float(inRec,"CRPIX2",&status);
    cdeltx   = drms_getkey_float(inRec,"CDELT1",&status);   // arcsec per pixel
    cdelty   = drms_getkey_float(inRec,"CDELT2",&status);
    rsun_ref = drms_getkey_float(inRec,"RSUN_REF",&status); // solar radius in meter
    dsun_obs = drms_getkey_float(inRec,"DSUN_OBS",&status); // distance from Sun to SDO in meter
    sunarc = asin(rsun_ref/dsun_obs)              // arc-sin in radian
           / 3.14159265358979e0 * 180.0 * 3600.0; // radian to arc-second
    printf("solar radius is %f in arcsec \n",sunarc);
    float fdummy;
    fdummy = sunarc / (cdeltx + cdelty) * 2.0;
    printf("solar radius is %f in CCD pix.\n",fdummy);
    if ((isnan(sunarc)) ||((isnan(crpixx)) || isnan(crpixy))) // in case something had gone wrong
    {
      sunarc = (float)(cols+rows) * 0.5;
      crpixx = (float)cols * 0.5 + 0.5;
      crpixy = (float)rows * 0.5 + 0.5;
    }
    sunarc = sunarc*sunarc;  // make it square
/* Map of invalid values (NaN or all-zero) , or off-disk */
    for (n = 0; n < imgpix; n++)
    {
      nan_map[n] = 0; // because nan_map was calloc-ed, this is not needed.
      double sumsqr;  // better be of double precision, maybe...
      sumsqr = 0.0;
      for (m = 0; m < nvar; m++){sumsqr = sumsqr + data[n + m*imgpix] * data[n + m*imgpix];}
      if (sumsqr < 1.0e-2){nan_map[n] = 1;} //     turn on flag to-be-skipped for having all-zero value
      if (isnan(sumsqr))  {nan_map[n] = 1;} //     turn on flag to-be-skipped for containing NaN
      float fpixdist, fxpix, fypix;
      int   ix, iy;
      ix = n % cols;
      iy = n / cols;
      fxpix = ((float)(ix) -(crpixx - 1.0)) * cdeltx;
      fypix = ((float)(iy) -(crpixy - 1.0)) * cdelty;
      fpixdist = fxpix * fxpix + fypix * fypix; // must be square
      if (fpixdist > sunarc) {nan_map[n] = 2;}  // turn on flag to-be-skipped for being out-of-disk
    }
    printf("data is read\n");
/* counting how many on-disk non-NaN pixel */
    int nonnan, numnan, nofdsk;
    nonnan = 0;
    numnan = 0;
    nofdsk = 0;
    for (n = 0; n < imgpix; n++)
    {
      if (nan_map[n] == 0) nonnan = nonnan + 1;
      if (nan_map[n] == 1) numnan = numnan + 1;
      if (nan_map[n] == 2) nofdsk = nofdsk + 1;
    }
    printf(" Num of pixel total (4k x 4k or 16 x 2^10)   : %8d \n", imgpix);
    printf(" Num of pixel out-of-disk                    : %8d \n", nofdsk);
    printf(" Num of pixel skipped due to all-zero or NaN : %8d \n", numnan);
    printf(" Num of pixel to be processed                : %8d \n", nonnan);

/* make equi-area (non-NaN) list */
    int irank;
    int *itmps;
    int *itmpe;
    itmps = (int *)malloc(sizeof(int) * mpi_size);
    itmpe = (int *)malloc(sizeof(int) * mpi_size);

    for(irank = 0; irank < mpi_size; irank++)
    {
      int myrank, nprocs, numpix;
      int istart, iend;
      myrank = irank;
      nprocs = mpi_size;
      numpix = nonnan;
      para_range(myrank,nprocs,numpix,&istart,&iend);
      itmps[irank] = istart;
      itmpe[irank] = iend;
    }
    int icount;
    icount = -1;
    for (n = 0; n < imgpix; n++)
    {
      if (nan_map[n] == 0)
      {
        icount = icount + 1;
        for (m = 0; m < mpi_size; m++)
        {
          if (itmps[m] == icount){istartall[m]=n;}
          if (itmpe[m] == icount){iendall[m]=n;}
        }
      }
    }
    free(itmps); // liberate
    free(itmpe);
  } // end if mpi_rank == 0
  else
  {
/* #ifndef DRMS_CLIENT */
    drms_server_end_transaction(drms_env, 1, 0);  // Kehcheng and Arts' suggestion
    db_disconnect(&drms_env->session->db_handle); // disconnect from non-primary PE to DRMS
/* #endif */
  }

  MPI_Barrier(MPI_COMM_WORLD);

/* at this moment, the PE(s) other than the primary do not know the value of imgpix etc.*/
  int ibuff[4];
  ibuff[0]=imgpix;  // packing
  ibuff[1]=nvar;
  ibuff[2]=cols;
  ibuff[3]=rows;
  MPI_Bcast(ibuff,4,MPI_INT,0,MPI_COMM_WORLD);
  imgpix = ibuff[0];
  nvar   = ibuff[1];
  cols   = ibuff[2];
  rows   = ibuff[3];
/* sending the istart-iend list array */
  MPI_Bcast(istartall,mpi_size,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(iendall,  mpi_size,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

/* large array ONLY at the primary */
  if (mpi_rank == 0)
  {
    FinalErr=(double *)malloc(sizeof(double)*imgpix*Err_ct);
    FinalRes=(double *)malloc(sizeof(double)*imgpix*paramct);
  }

/* part of input/output ; local array for ALL PE. */
  float  *dataLocal;
  double *FinalResLocal,*FinalErrLocal;
  int    *nan_mapLocal;
  myrank = mpi_rank;
  nprocs = mpi_size;
  numpix = imgpix;
#if EQUIAREA == 1
  istart=istartall[myrank];
  iend=iendall[myrank];
#else
  para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
  if (verbose)
  {
    printf("Hello, this is %2d th PE of %2d, in charge of pixels from %9d to %9d of 0 to %9d.\n",
           mpi_rank,mpi_size,istart,iend,(imgpix-1));
  }
  int imgpixlocal;
  imgpixlocal = iend - istart + 1;
  nan_mapLocal =(int    *)malloc(sizeof(int)   *imgpixlocal);
  dataLocal    =(float  *)malloc(sizeof(float) *imgpixlocal*nvar);
  FinalErrLocal=(double *)malloc(sizeof(double)*imgpixlocal*Err_ct);
  FinalResLocal=(double *)malloc(sizeof(double)*imgpixlocal*paramct);

/* tiny arrays used at processing each pixel as the invert()'s input/output arguments. */
  segsize = wlct * spct * npix;
  ressize = paramct * npix;
  obs = (double *)malloc (sizeof (double) * nvar);
  res = (double *)calloc (paramct, sizeof (double));
  scat= (double *)malloc (sizeof (double) * nvar);
  err = (double *)calloc (Err_ct,sizeof (double));

  MPI_Barrier(MPI_COMM_WORLD);

/* send partial input data to each PE */
  if (mpi_rank == 0)
  {
/* first, the primary makes copy for its own part */
    myrank = mpi_rank;
    nprocs = mpi_size;
    numpix = imgpix;
#if EQUIAREA == 1
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    for (n = istart ; n < iend+1 ; n++){for (m = 0; m < nvar; m++){dataLocal[(n-istart)*nvar+m] = data[n + m*imgpix];}}
/* then send the partials to the other PEs */
    if (mpi_size > 1)
    {
      int irank;
      for(irank = 1; irank < mpi_size; irank++)
      {
        int mpi_dest;
        int ibufsize;
        float *fbufsend;
        mpi_dest = irank;
#if EQUIAREA == 1
        istart=istartall[mpi_dest];
        iend=iendall[mpi_dest];
#else
        para_range(mpi_dest,nprocs,numpix,&istart,&iend);
#endif
        ibufsize = (iend-istart+1) * nvar;
        fbufsend= (float*)malloc(sizeof(float) * ibufsize);
        for (n = istart ; n < iend+1 ; n++){for (m = 0; m < nvar; m++){fbufsend[(n-istart)*nvar+m] = data[n + m*imgpix];}}
        mpi_tag = 1400 + irank;
        MPI_Send(fbufsend, ibufsize, MPI_REAL, mpi_dest, mpi_tag, MPI_COMM_WORLD);
        free(fbufsend);
      }
    }
  }
  else
  {
    int mpi_from = 0;
    int ibufsize;
    float *fbufrecv;
#if EQUIAREA == 1
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    ibufsize = (iend-istart+1) * nvar;
    fbufrecv = (float*)malloc(sizeof(float) * ibufsize);
    mpi_tag = 1400 + mpi_rank;
    MPI_Recv(fbufrecv, ibufsize, MPI_REAL, mpi_from, mpi_tag, MPI_COMM_WORLD, &mpistat);
    for (n = istart ; n < iend+1 ; n++){for (m = 0; m < nvar; m++){dataLocal[(n-istart)*nvar+m] = fbufrecv[(n-istart)*nvar+m];}}
    free(fbufrecv);
  } // end-if mpi_rank is 0, or not.

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) printf("input data had propagated to all PE.\n");

/* send partial non-NAN mask-map from primary PE to the others */
  if (mpi_rank == 0)
  {
/* first, the primary makes copy for its own part */
    myrank = mpi_rank;
    nprocs = mpi_size;
    numpix = imgpix;
#if EQUIAREA == 1
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    for (n = istart ; n < iend+1 ; n++){nan_mapLocal[n-istart] = nan_map[n];}
/* then send the partial to each PE */
    if (mpi_size > 1)
    {
      int irank;
      for(irank = 1; irank < mpi_size; irank++)
      {
        int mpi_dest;
        int ibufsize;
        int *ibufsend;
        mpi_dest = irank;
#if EQUIAREA == 1
        istart=istartall[mpi_dest];
        iend=iendall[mpi_dest];
#else
        para_range(mpi_dest,nprocs,numpix,&istart,&iend);
#endif
        ibufsize = (iend-istart+1) * 1;
        ibufsend= (int*)malloc(sizeof(int) * ibufsize);
        for (n = istart ; n < iend+1 ; n++){ibufsend[n-istart] = nan_map[n];}
        mpi_tag = 1500 + irank;
        MPI_Send(ibufsend, ibufsize, MPI_INTEGER, mpi_dest, mpi_tag, MPI_COMM_WORLD);
        free(ibufsend);
      }
    }
  }
  else
  {
/* non-primary PEs wait until served */
    int mpi_from = 0;
    int ibufsize;
    int *ibufrecv;
#if EQUIAREA == 1
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    ibufsize = (iend-istart+1) * 1;
    ibufrecv = (int*)malloc(sizeof(int) * ibufsize);
    mpi_tag = 1500 + mpi_rank;
    MPI_Recv(ibufrecv, ibufsize, MPI_INTEGER, mpi_from, mpi_tag, MPI_COMM_WORLD, &mpistat);
    for (n = istart ; n < iend+1 ; n++){nan_mapLocal[n-istart] = ibufrecv[n-istart];}
    free(ibufrecv);
  } // end-if mpi_rank is 0, or not.

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) printf("mask  data had propagated to all PE.\n");

/* inversion initializations : must be done by each PE. */
  if (mpi_rank == 0) printf("\n----------- inversion initializations ----------------- \n");

  line_init_(&LAMBDA_0,&LAMBDA_B,&NOISE_LEVEL);
  if (verbose){printf("done line_init for mpi_rank %d\n",mpi_rank);}
  wave_init_ (&LAMBDA_MIN,&DELTA_LAMBDA,&NUM_LAMBDA);
  if (verbose){printf("done wave_init for mpi_rank %d\n", mpi_rank );}
  filt_init_ (&NUM_LAMBDA_FILTER,&NUM_TUNNING,&CONTINUUM,&LYOTFWHM,&WNARROW,&WSPACING);
  
  /* JM Borrero & RC Elliot Apr 7, 2010 */
  /* filt_init should be replaced for a routine, written by Sebastien
     that will read from the DMRS the phases, contrasts, etcetera: ONLY READ */

  if (verbose){printf("done filt_init for mpi_rank %d\n",mpi_rank);}
  inv_init_(&NUM_ITERATIONS,&SVD_TOLERANCE,&CHI2_STOP,&POLARIZATION_THRESHOLD,&INTENSITY_THRESHOLD,&PERCENTAGE_JUMP);
  if (verbose){printf("done inv_init\n");}
  free_init_(list_free_params);
  if (verbose){printf("done list_free_params for mpi_rank %d\n", mpi_rank);}
  svd_init_();
  if (verbose){printf("done svd_init\n");}
/* By RCE & Borrero */
/* Changed index of list_free_params to 3 (instead of 4) to refer to Damping*/
  if  (list_free_params[3] == 0.0) voigt_init_();
  for (n=0; n< nvar; n++) scat[n]=0;

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) printf("\n----------- inversion initializations done ------------ \n");

/* at last, do inversion in parallel */
  if (mpi_rank == 0) printf("\n----------- finally, do inversion in parallel --------- \n");

  myrank = mpi_rank;
  nprocs = mpi_size;
  numpix = imgpix;
#if EQUIAREA == 1
  istart=istartall[myrank];
  iend  =iendall[myrank];
#else
  para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
  int pixdone;
  pixdone = 0;
  for (n = istart; n < iend + 1; n++)
  {
    if (nan_mapLocal[n-istart] == 0)
    {
#if QUICKRUN == 1
//     if ((n % rows > 2298) && (n % rows < 2400)) // intentionally off-center
     if ((n % rows > 1250) && (n % rows < 1650)) // central 100-pixel width column
     {
#endif

      for(m = 0; m < nvar; m++){obs[m]=dataLocal[(n-istart)*nvar+m];}

      /* JM Borrero & RC Eliiot Apr 7, 2010 */
      /* Right before calling invert_ there should be a call to
	 another of Sebastien's program to calculate the filter profiles
	 for the particular pixel "n". The output filters should be in the
	 form: double (Num_lambda_filter,num_lambda). Note that we want the
	 filters evaluated at the following wavelenghts:
	 FOR K=1,Num_Lambda DO WAVELENGTH[K] = Lambda_Min + K*Delta_Lambda/1000 [Units are Angstroms]
	 Lambda_Min, Delta_Lambda, Num_lambda_filter & num_lambda are all defined above */
      /* Note: I would avoid sending this program the full 512*512 matrixes (can slow down teh code)
	 but rather the 4 neighboring points to perform the bilinear interpolation */


      invert_ (obs, scat, guess, res, err);                                  // do inversion

      /* JM Borrero & RC Elliot Apr 7, 2010 */
      /* invert_ should be called in a new way:
	 invert_ (obs,scat,guess,filter,res,err)
	 where filter should be a 2D double array with
	 dimensions (Num_lambda_filter,num_lambda) eg (6,49) */

      for (j=0; j<paramct; j++){FinalResLocal[(n-istart)*paramct+j]=res[j];} // copy results to the local-small array(s)
      for (k=0; k<Err_ct;  k++){FinalErrLocal[(n-istart)*Err_ct +k]=err[k];}
#if QUICKRUN == 1
     }
     else
     {
      for (j=0; j<paramct; j++){FinalResLocal[(n-istart)*paramct+j]= 1.0+mpi_rank;} // put dummy
      for (k=0; k<Err_ct ; k++){FinalErrLocal[(n-istart)*Err_ct +k]=-1.0-mpi_rank;}
     }
#endif
      pixdone = pixdone + 1;
    }
    else
    {
      for (j=0; j<paramct; j++){FinalResLocal[(n-istart)*paramct+j]=NAN;} // put NAN
      for (k=0; k<Err_ct;  k++){FinalErrLocal[(n-istart)*Err_ct +k]=NAN;}
    } // end of if (NAN) etc.
  } // end of n-loop
  if (verbose){printf("Hello, this is %2d th PE : inversion done for %9d pixels. \n", mpi_rank, pixdone);}
  MPI_Barrier(MPI_COMM_WORLD);

/* output data are gathered to primary from each PE in charge of pixels from istart to iend */
  if (mpi_rank == 0)
  {
/* first, copy the portion the primary itself did */
    myrank = mpi_rank;
    nprocs = mpi_size;
    numpix = imgpix;
#if EQUIAREA == 1
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    for (n = istart ; n < iend+1 ; n++)
    {
      for (j=0; j<paramct; j++){FinalRes[(n*paramct)+j]=FinalResLocal[(n-istart)*paramct+j];}
      for (k=0; k<Err_ct ; k++){FinalErr[(n*Err_ct) +k]=FinalErrLocal[(n-istart)*Err_ct +k];}
    }
/* then collecting the portions done by the others */
    if (mpi_size > 1)
    {
      printf("now collecting data from PEs : ");
      int irecv;
      for (irecv = 1; irecv < mpi_size ; irecv++)
      {
        printf(" %d ",irecv);
        int mpi_from;
        nprocs = mpi_size;
        numpix = imgpix;
#if EQUIAREA == 1
        istart=istartall[irecv];
        iend=iendall[irecv];
#else
        para_range(irecv,nprocs,numpix,&istart,&iend);
#endif
        int ibufsize;
        ibufsize = (iend-istart+1) * (paramct+Err_ct);
        double *dbufrecv;
        dbufrecv = (double*)malloc(sizeof (double) * ibufsize);
        mpi_from = irecv;
        mpi_tag = 1200 + irecv;
        MPI_Recv(dbufrecv, ibufsize, MPI_DOUBLE, mpi_from, mpi_tag, MPI_COMM_WORLD, &mpistat);
        for (n = istart ; n < iend+1 ; n++)
        {
          for (j=0; j<paramct; j++){FinalRes[(n*paramct)+j]=dbufrecv[(n-istart)*(paramct+Err_ct)        +j];}
          for (k=0; k<Err_ct ; k++){FinalErr[(n*Err_ct) +k]=dbufrecv[(n-istart)*(paramct+Err_ct)+paramct+k];}
        }
        free(dbufrecv);
      }
      printf("done \n");
   }
  }
  else
  {
    int isend;
    int mpi_dest = 0;
    nprocs = mpi_size;
    numpix = imgpix;
    isend = mpi_rank;
#if EQUIAREA == 1
    istart=istartall[isend];
    iend=iendall[isend];
#else
    para_range(isend,nprocs,numpix,&istart,&iend);
#endif
    int ibufsize;
    ibufsize = (iend-istart+1) * (paramct+Err_ct);
    double *dbufsend;
    dbufsend = (double*)malloc(sizeof (double) * ibufsize);
    for (n = istart ; n < iend + 1 ; n++)
    {
      for (j=0; j<paramct; j++){dbufsend[(n-istart)*(paramct+Err_ct)        +j]=FinalResLocal[(n-istart)*paramct+j];}
      for (k=0; k<Err_ct;  k++){dbufsend[(n-istart)*(paramct+Err_ct)+paramct+k]=FinalErrLocal[(n-istart)*Err_ct +k];}
    }
    mpi_tag = 1200 + mpi_rank;
    MPI_Send(dbufsend, ibufsize, MPI_DOUBLE, mpi_dest, mpi_tag, MPI_COMM_WORLD);
    free(dbufsend);
  } // end-if mpi_rank is 0, or not.
  MPI_Barrier(MPI_COMM_WORLD); // silly..but always safe

/* the primary PE takes care of pixels none of PE had taken care of */
  if (mpi_rank == 0)
  {
    if (istartall[0] > 0)
    {
      for (n = 0 ; n < istartall[0]; n++)
      {
        for (j=0; j<paramct; j++){FinalRes[(n*paramct)+j]=NAN;}
        for (k=0; k<Err_ct ; k++){FinalErr[(n*Err_ct) +k]=NAN;}
      }
    }
    if (iendall[mpi_size-1] < imgpix - 1)
    {
      for (n = iendall[mpi_size-1] + 1; n < imgpix; n++)
      {
        for (j=0; j<paramct; j++){FinalRes[(n*paramct)+j]=NAN;}
        for (k=0; k<Err_ct ; k++){FinalErr[(n*Err_ct) +k]=NAN;}
      }
    }
  } // end-if mpi_rank is 0.

  MPI_Barrier(MPI_COMM_WORLD); // silly..but always safe
  free(FinalResLocal);         // liberate....
  free(FinalErrLocal);
  free(dataLocal);
  free(nan_mapLocal);
  MPI_Barrier(MPI_COMM_WORLD);

/* write outputs through DRMS/JSOC */
  if (mpi_rank == 0)
  {
    outRec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    if (!outRec) {fprintf (stderr, "Error creating record in series %s; abandoned\n",outser);return 1;}

/* succeed a lot of info. from the input data record */
#if 0
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_All);
#else
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit); // Phil's solution !
//    fprintf(stderr,"T_REC_step=%f\n",drms_getkey_double(outRec,"T_REC_step",0));
#endif
#if 0
    drms_copykey(outRec, inRec, "CENTER_X"); // itemized
    drms_copykey(outRec, inRec, "CENTER_Y");
    drms_copykey(outRec, inRec, "SOLAR_B0");
    drms_copykey(outRec, inRec, "SOLAR_P0");
    drms_copykey(outRec, inRec, "R_SUN");
    drms_copykey(outRec, inRec, "TELESCOP");
    drms_copykey(outRec, inRec, "INSTRUME");
    drms_copykey(outRec, inRec, "WAVELENG");
    drms_copykey(outRec, inRec, "T_REC");
    drms_copykey(outRec, inRec, "T_OBS");
#endif
#if 0
    drms_setkey_string(outRec, "T_REC", "2000.02.02_02:02:02_TAI"); // enforce T_REC dummy ones
    drms_setkey_string(outRec, "T_OBS", "2000.02.02_02:02:02_TAI"); //    also T_OBS
#endif

    char trectmp2[26];
    TIME trectmp1 = drms_getkey_time(outRec,"T_REC",&status);
    sprint_time(trectmp2,trectmp1,"TAI",0);
    printf("created record %s[%s] \n",outser,trectmp2);

/* output array will be split to individual array for each j-th parameter */
    for (j = 0; j < paramct; j++)
    {
      float *dat1;
      dat1 = (float*) calloc(imgpix, sizeof(float));
/* collect the result for each parameter over all pixels */
      for(n = 0; n < imgpix ; n++){dat1[n] = FinalRes[(n*paramct)+j];}
/* should get this info from segment lookup in output record */
      int axes[2];
      axes[0] = cols;
      axes[1] = rows;
      invrt_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, dat1, &status);
      seg = drms_segment_lookup (outRec, Resname[j]);
      if (!seg) {fprintf(stderr, "Error getting data segment %s; abandoned\n", Resname[j]);}
      if (drms_segment_write (seg, invrt_array, 0))
      {
        fprintf (stderr, "Error writing segment %d (%s); abandoned\n", j,seg->info->name);
        return 1;
      }
      else
      {
        if (verbose){printf("Results written out to %-s\n", Resname[j]);}
      }
      free(dat1);
    } // end of j-loop

    printf("printing errors now\n");
    for (k = 0; k  < Err_ct; k++)
    {
      float *dat2;
      dat2 = (float*) calloc(imgpix , sizeof(float));
      int axes[2];
      for (n=0; n < imgpix;n++){dat2[n] = FinalErr[(n*Err_ct)+k];}
      axes[0] = cols;
      axes[1] = rows;
      err_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, dat2,&status);
      seg = drms_segment_lookup (outRec, Resname[k+paramct]);
      if (!seg)
      {
        if (verbose){fprintf(stderr, "Error getting data segment %s; abandoned\n", Resname[k+paramct]);}
      }
      else
      {
        if (drms_segment_write (seg, err_array, 0))
        {
          fprintf (stderr, "Error writing to segment %d (%s); abandoned\n", k+paramct, seg->info->name);
          return 1;
        }
        else
        {
          if (verbose){printf("Errors  written out to %-s\n", Resname[k+paramct]);}
        }
      }
      free(dat2);
    } // end of k-loop
    printf("write-out done !\n");

    printf("so, close all DRMS record(s) !\n");
/* DRMS trailer and closer */
    drms_close_record (outRec, DRMS_INSERT_RECORD);
    drms_close_records (records, DRMS_FREE_RECORD);

/* how long it took */
    time(&endtime);
    printf ("%ld sec for %d profiles\n", endtime - startime, imgpix);
    printf ("%.2f profiles per second\n", (float)(npix*mpi_size) / (0.01 + (float)(endtime - startime)));

    printf("good bye !\n");
  } // end-if mpi_rank is 0.

  MPI_Barrier(MPI_COMM_WORLD);

/* say good bye to MPI things.*/
  MPI_Finalize();

  return 0;
}

/* ------------------------------- end of main-wrapper layer ------------------------------- */


/* ------------------------------- subprograms ------------------------------- */


/* --------------------------------------------------------------------------------
 *
 * Calculate the first and the last addresses of array each PE will be in charge of.
 *
 * -------------------------------------------------------------------------------- */

void para_range(int myrank, int nprocs, int numpix, int *istart, int *iend)
{
  int iwork1, iwork2, imin, idummy, n1, n2;
  n1 = 0;              // When written in Fortran, this is typically 1.
  n2 = numpix - 1;     // When written in Fortran, this is typically numpix.

  iwork1 = (n2 - n1 + 1) / nprocs;
  iwork2 = (n2 - n1 + 1) % nprocs; // mod(n2-n1+1,nprocs)
  if (iwork2 >= myrank){imin = myrank;}else{imin=iwork2;} // imin = min(myrank,iwork2)
  *istart = myrank*iwork1 + n1 + imin;
  idummy  = *istart + iwork1 - 1;
  if (iwork2 <= myrank){*iend = idummy;}else{*iend=idummy+1;}
}

/* ----------------------------- end of this file ----------------------------- */
