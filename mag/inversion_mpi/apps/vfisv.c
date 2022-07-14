 /*
 *
 *  vfisv.c                                        ~rick/hmi/vecb/src/v09
 * 3/25/10  this is Keiji's less memory usage version with my edits for keywords * and k=0.
 * 3/23/10 THIS IS a copy of REBECCA's /v13  except with vfisv_mpi
 * this is the mpi version  vfisv.c version 10. It should have the keywords correct AND 11 err files.
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
 * --- some comments added by K.H. -----
 *
 * 1) March 26, 2010
 *   First, combining efforts done by 5:00PM, March 26, 2010
 *    (1) a lot of things such as JSOC/DRMS things done by Priya and Rick
 *          /home/priya/vecb/VFISV_new/v11/vfisv.c
 *    (2) Default parameters defined by Rebecca
 *          /home/rce/v16/vfisv.c
 *    (3) MPI things by Priya, Rick and Keiji
 */

/* 
 * To speed up, only some rectangle portion be processed.
 * Please set zero (0) for formal/normal running
 */
#define KEIJISKIP 1
#define EQUIAREA  1


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
  {ARG_STRING,  "in", "su_couvidat.hmiStokes[2007.10.14_23:50:00_TAI]", "input series"},  /*This is the input data in the form of 24 fits file segnments */
  {ARG_STRING,  "out", "hao_rce.KD4096_Mar26", "output series"},
  {ARG_INT,  "npix", "0", "number of pixels per segment (default: all/segs)"},
  {ARG_INT,  "num_iter", "30", "number of iterations(default: 30)"},
  {ARG_INT,  "num_lambda", "49", "number of ??(default: 33)"},
  {ARG_INT,  "Num_lambda_filter", "6", "Number of filters accross the wvl (default: 6)"},
  {ARG_INT,  "Num_tunning", "6", "Number of ??(default: 6)"},
  {ARG_DOUBLE,  "svd_tolerance", "1.0e-32", "svd tolerance (default: 1.0e-32)"},
  {ARG_DOUBLE,  "chi2_stop", "1.0e-6", "chisq-stop (default: 1.0e-6)"},
  {ARG_DOUBLE,  "Polarization_threshold", "1.0e-2", "polarization threshold (default: 0.01)"},
  {ARG_DOUBLE,  "Intensity_Threshold", "1e2", "Intensity threshold (default: 0.8)"},
/*{ARG_DOUBLE,  "Intensity_Threshold", "0.8", "Intensity threshold (default: 0.8)"},*/
  {ARG_DOUBLE,  "Percentage_Jump", "10.0", "Percentage Jump (default: 10%)"},
  {ARG_DOUBLE,  "Lambda_Min", "-648.0", "Intensity threshold (default: -432)"},
  {ARG_DOUBLE,  "Lambda_0", "6173.3433", "Wavelength(default:6173.3433 Angstrom )"},
  {ARG_DOUBLE,  "Lambda_B", "0.044475775", "FWHM?? (default: 0.044475775)"},
  {ARG_DOUBLE,  "Delta_Lambda", "27.0", "Delta Lambda(default: 27.0)"},
  {ARG_DOUBLE,  "Noise_LEVEL", "5.0e1", "Intensity threshold (default: 3.0e-3)"},
/*{ARG_DOUBLE,  "Noise_LEVEL", "5.0e3", "Intensity threshold (default: 3.0e-3)"},*/
  {ARG_DOUBLE,  "Continuum", "0.0", "Intensity threshold (default: 0)"},
  {ARG_DOUBLE,  "Lyotfwhm", "424.0", "Lyot filter FWHM (default: 424.0)"},
  {ARG_DOUBLE,  "Wnarrow", "172.0", "W narrow (default: 172.0)"},
  {ARG_DOUBLE,  "Wspacing", "69.0", "W narrow (default: 69.0)"},
  {ARG_FLAG, "v",    "", "run verbose"},
  {ARG_FLAG, "d",    "", "turn damping off"},
  {}
};

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *records, *out;
  DRMS_Record_t *inRec, *outRec;
  DRMS_Segment_t *seg;
  DRMS_Array_t *stokes_array, *invrt_array,*err_array;

  time_t start, end;
  double *ddat, *res;
  double *obs, *scat,*err;
  double guess[10]= {15.0,90.0,45.0,0.05,40.0,50.0,0.0,0.4*1e6,0.6*1e6,1.0};
  double *FinalRes,*FinalErr;
  double SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD, INTENSITY_THRESHOLD, PERCENTAGE_JUMP;
  double LAMBDA_0, LAMBDA_B, NOISE_LEVEL;
  double LAMBDA_MIN, DELTA_LAMBDA;
  double LYOTFWHM, WNARROW, WSPACING;
  double sumchi, avgval, minval, maxval;
  float *data, *data0, *inv = NULL;
  int *nan_map;
  float *stokesi, *stokesq, *stokesu, *stokesv;
  float fval;
  int list_free_params[10]={1,1,1,1,1,1,1,1,1,1};
  int param[3];
  int rn, rec_ct, sn, seg_ct;
  int cols, rows, imgpix, imgbytes, imgloc;
  int sp, wl, nvar;
  int imgseg, segsize, ressize;
  int j,m,i,k,n ;
  int status;
  int NUM_ITERATIONS;
  int NUM_LAMBDA;
  int NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM;
  char *Resname[] = {"eta_0", "inclination", "azimuth", "damping", "dop_width",
                     "field", "vlos_mag", "src_continuum", "src_grad", "alpha_mag","field_err","inclination_err","azimuth_err", "vlos_err", "alpha_err","field_inclination_err","field_az_err","inclin_azimuth_err", "field_alpha_err", "inclination_alpha_err", "azimuth_alpha_err" };

        /*  if errors = .true. is error in filling factor propagated?  */
                                        /*  Documentation suggest not  */

  char segname[16];
                                                        /*  may eliminate  */
  char *infile;
                                                        /*  constants  */
  char *spname[] = {"I", "Q", "U", "V"};
  int spct = (sizeof (spname) / sizeof (char *));
  int wlct = 6;
  int paramct = 10; /*is it 16 or 10??*/
  int Err_ct=11;                                                      /*  parameters  */
  int npix = params_get_int (params, "npix");
  NUM_ITERATIONS = params_get_int (params, "num_iter");
  NUM_LAMBDA = params_get_int (params, "num_lambda");
  NUM_LAMBDA_FILTER = params_get_int (params, "Num_lambda_filter");
  NUM_TUNNING = params_get_int(params, "Num_tunning");
  SVD_TOLERANCE=params_get_double (params, "svd_tolerance");
  CHI2_STOP=params_get_double (params, "chi2_stop");
  POLARIZATION_THRESHOLD=params_get_double (params, "Polarization_threshold");
  INTENSITY_THRESHOLD=params_get_double (params, "Intensity_Threshold");
  PERCENTAGE_JUMP=params_get_double (params, "Percentage_Jump");
  LAMBDA_MIN=params_get_double (params, "Lambda_Min");
  LAMBDA_0=params_get_double (params, "Lambda_0");
  LAMBDA_B=params_get_double (params, "Lambda_B");
  DELTA_LAMBDA=params_get_double (params, "Delta_Lambda");
  NOISE_LEVEL=params_get_double (params, "Noise_LEVEL");
  LYOTFWHM=params_get_double (params, "Lyotfwhm");
  WNARROW=params_get_double (params, "Wnarrow");
  WSPACING=params_get_double (params, "Wspacing");
  CONTINUUM=params_get_double (params, "Continuum");
  char *indsdesc = strdup (params_get_str (params, "in"));
  int verbose = params_isflagset (params, "v");
  char *outser =strdup ( params_get_str (params, "out"));
/* MPI things .... */
  int segct;
  int mpi_rank, mpi_size;
  MPI_Status mpistat;
  int mpierr, mpilabel;
  void para_range(int,int,int,int *,int *);
  int *istartall, *iendall;

                                /*  Initialize Clock  */
  time (&start);
  printf("Lambda_O= %f\n",LAMBDA_0);

  MPI_Status mpi_stat;   /* Initalizing MPI */
  status = 0;
  MPI_Init (&status, NULL);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);

  istartall = (int *)malloc(sizeof(int) * mpi_size);
  iendall   = (int *)malloc(sizeof(int) * mpi_size);

  MPI_Barrier(MPI_COMM_WORLD);

  segct = mpi_size;
  if (mpi_rank == 0) {
    if (verbose) printf ("%s V %s\n", module_name, version_id);
                                                /*  open input record_set  */
    records = drms_open_records (drms_env, indsdesc, &status);
    if (status)
    {
      fprintf (stderr, "drms_open_records failed, in=%s, status=%d.  Aborting.\n", indsdesc, status);
      return 1;
    }
    if ((rec_ct = records->n) == 0)
    {
      fprintf (stderr, "No records in selected dataset %s  Aborting.\n", indsdesc);
      return 1;
    }
    if ((rec_ct = records->n) > 1)
    {
      fprintf (stderr, "Warning: only first record in selected set processed\n");
    }
    rn = 0;

    records = drms_open_records (drms_env, indsdesc, &status);
    if (status) {
      fprintf (stderr, "drms_open_records failed, in=%s, status=%d.  Aborting.\n", indsdesc, status);
      return 1;
    }
    if ((rec_ct = records->n) == 0)
    {
      fprintf (stderr, "No records in selected dataset %s  Aborting.\n", indsdesc);
      return 1;
    }
    if ((rec_ct = records->n) > 1)
    {
      fprintf (stderr, "Warning: only first record in selected set processed\n");
    }
    rn = 0;
    inRec = records->records[rn];
    seg = drms_segment_lookupnum (inRec, 0);
    cols = seg->axis[0];
    rows = seg->axis[1];
    imgpix = cols * rows;
    imgbytes = imgpix * sizeof (float);
    nvar = wlct * spct;

    if (npix < 1) npix = imgpix / segct;
    if (imgpix % npix) {printf ("Error: npix (%d) must be a divisor of %d\n", npix, imgpix);return 0;}

    data = data0 = (float *)malloc (sizeof (float) * imgpix * nvar);
    nan_map=(int*)calloc (imgpix, sizeof (int));

    /*spct=4,wlct=6    */
    for (sp = 0; sp < spct; sp++)
    {
      for (wl = 0; wl < wlct; wl++)
      {
        sprintf (segname, "%s%d", spname[sp], wl);
        if ((seg = drms_segment_lookup (inRec, segname)) == NULL){
          fprintf (stderr, "Error reading segment %s of record %d\n", segname, rn);
          return 1;
        }
        stokes_array = drms_segment_read (seg, DRMS_TYPE_FLOAT, &status);
        /* data is read in 1 segment at a time-so for each location;  */
        /* printf("segment read %s\n",segname); */
        /* the 24 requiired values are 4096*4096 datapoints apart */
        memcpy (data, stokes_array->data, imgbytes);
        drms_free_array (stokes_array);
        data += imgpix;
      }
    }
    data = data0;
    printf("Imgpix= %d\n",imgpix);
/* Map of pixels where input data are Nan */
    for (n = 0; n < imgpix; n++)
    {
      nan_map[n] = 0; // because nan_map was calloc-ed, this is not needed.
      double sumsqr;
      sumsqr = 0.0;
      for (m = 0; m < nvar; m++){sumsqr = sumsqr + data[n + m*imgpix] * data[n + m*imgpix];}
      if (sumsqr < 1.0e-2){nan_map[n] = 1;} // turn on flag to-be-skipped if all data is (almost) zero
      if (isnan(sumsqr))  {nan_map[n] = 1;} // turn on flag to-be-skipped if data contain NaN.
    }
    printf("data is read\n");
/* now counting how many non-NaN */ 
    int nonnan, numnan;
    nonnan = 0;
    numnan = 0;
    for (n = 0; n < imgpix; n++)
    {
      if (nan_map[n] == 0) nonnan = nonnan + 1;
      if (nan_map[n] == 1) numnan = numnan + 1;
    }
    printf(" Num of pixel total                          : %d \n", imgpix);
    printf(" Num of pixel to be processed                : %d \n", nonnan);
    printf(" Num of pixel skipped due to all-zero or NaN : %d \n", numnan);

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
  } // end if mpi_rank == 0
  else
  {
/* #ifndef DRMS_CLIENT */
    drms_server_end_transaction(drms_env, 1, 0); // Kehcheng and Arts
    db_disconnect(&drms_env->session->db_handle); 
/* #endif */
  }

/* MIND that MPI_Barrier is time-consuming and computationally-expensive,
   but makes me feel safe.... so I do use it frequently, unnecessarily */
  MPI_Barrier(MPI_COMM_WORLD);

/* at this moment, the PE(s) other than the primary do not know the value of imgpix etc.*/
  if (mpi_rank == 0)
  {
    int mpi_trgt;
    int ibufsend[4];
    ibufsend[0]=imgpix;
    ibufsend[1]=nvar;
    ibufsend[2]=cols;
    ibufsend[3]=rows;
    for (mpi_trgt=1; mpi_trgt < mpi_size; mpi_trgt++)
    {
      mpilabel = 1000 + mpi_trgt; // unique number in code.
      MPI_Send(ibufsend, 4, MPI_INT, mpi_trgt, mpilabel, MPI_COMM_WORLD);
    }
  }
  else
  {
    int mpi_from = 0;
    int ibufrecv[4];
    mpilabel = 1000 + mpi_rank;  // must have counterpart one at MPI_Send's argument.
    MPI_Recv(ibufrecv, 4, MPI_INT, mpi_from, mpilabel, MPI_COMM_WORLD, &mpistat);
    imgpix = ibufrecv[0];
    nvar   = ibufrecv[1];
    cols   = ibufrecv[2];
    rows   = ibufrecv[3];
  }
  MPI_Barrier(MPI_COMM_WORLD);

/* now sending the istart-iend list */
  if (mpi_rank == 0)
  {
    int ibufsize;
    int *ibufsend;
    int mpi_trgt;
    ibufsize = 2 * mpi_size;
    ibufsend = (int *)malloc(sizeof(int) * ibufsize);
    for (i=0; i<mpi_size;i++){ibufsend[i*2] = istartall[i];ibufsend[i*2+1]=iendall[i];}
    for (mpi_trgt=1; mpi_trgt < mpi_size; mpi_trgt++)
    {
      mpilabel = 1600 + mpi_trgt; // unique number in code.
      MPI_Send(ibufsend, ibufsize, MPI_INT, mpi_trgt, mpilabel, MPI_COMM_WORLD);
    }
  }
  else
  {
    int mpi_from = 0;
    int ibufsize;
    int *ibufrecv;
    ibufsize = 2 * mpi_size;
    ibufrecv = (int *)malloc(sizeof(int) * ibufsize);
    mpilabel = 1600 + mpi_rank;  // must have counterpart one at MPI_Send's argument.
    MPI_Recv(ibufrecv, ibufsize, MPI_INT, mpi_from, mpilabel, MPI_COMM_WORLD, &mpistat);
    for (i=0; i<mpi_size;i++){istartall[i]=ibufrecv[i*2];iendall[i]=ibufrecv[i*2+1];}
  }
  MPI_Barrier(MPI_COMM_WORLD);

/* large array ONLY at the primary */
  if (mpi_rank == 0)
  {
    FinalErr=(double *)malloc(sizeof(double)*imgpix*Err_ct);
    FinalRes=(double *)malloc(sizeof(double)*imgpix*paramct);
  }

/* small array for ALL PE. */
  float  *dataLocal;
  double *FinalResLocal,*FinalErrLocal;
  int    *nan_mapLocal;
  int istart, iend;
  int myrank, nprocs, numpix;
  myrank = mpi_rank; // copy to protect original... not needed, just for safety
  nprocs = mpi_size;
  numpix = imgpix;
#if EQUIAREA
  istart=istartall[myrank];
  iend=iendall[myrank];
#else
  para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
  printf("Hello, this is %d th PE of %d, in charge of pixel from %d to %d of 0 to %d \n", mpi_rank,mpi_size,istart,iend,(imgpix-1));
  int imgpixlocal;
  imgpixlocal = iend - istart + 1;
  nan_mapLocal =(int    *)malloc(sizeof(int)   *imgpixlocal);
  dataLocal    =(float  *)malloc(sizeof(float) *imgpixlocal*nvar);
  FinalErrLocal=(double *)malloc(sizeof(double)*imgpixlocal*Err_ct);
  FinalResLocal=(double *)malloc(sizeof(double)*imgpixlocal*paramct);

/* now declare arrays that will be used at processing each pixel as the invert() argument. */
  segsize = wlct * spct * npix;
  ressize = paramct * npix;
  obs = (double *)malloc (sizeof (double) * nvar);
  res = (double *)calloc (paramct, sizeof (double));
  scat= (double *)malloc (sizeof (double) * nvar);
  err = (double *)calloc (Err_ct,sizeof (double));
/* something inversion initializations */
  line_init_(&LAMBDA_0, &LAMBDA_B,&NOISE_LEVEL);
  printf("done line_init for  mpirank %d\n",mpi_rank);
  wave_init_ (&LAMBDA_MIN, &DELTA_LAMBDA, &NUM_LAMBDA);
  printf("done wave_init  for  mpirank %d\n", mpi_rank );
  filt_init_ (&NUM_LAMBDA_FILTER, &NUM_TUNNING, &CONTINUUM, &LYOTFWHM, &WNARROW, &WSPACING);
  printf("done filt_init for mpi_rank %d\n",mpi_rank);
  inv_init_(&NUM_ITERATIONS, &SVD_TOLERANCE, &CHI2_STOP, &POLARIZATION_THRESHOLD, &INTENSITY_THRESHOLD, &PERCENTAGE_JUMP);
  printf("done inv_init\n");
  free_init_(list_free_params);
  printf("done list_free_params for mpi_rank %d\n", mpi_rank);
  svd_init_();
  printf("done svd_init\n  ");
  if  (list_free_params[4] == 0.0) voigt_init_();
  for (n=0; n< nvar; n++) scat[n]=0;

  MPI_Barrier(MPI_COMM_WORLD);

/* now send partial input data to all PE */
  if (mpi_rank == 0)
  {
/* first, the primary makes copy for its own part */
    myrank = mpi_rank;
    nprocs = mpi_size;
    numpix = imgpix;
#if EQUIAREA
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    for (n = istart ; n < iend+1 ; n++){for (m = 0; m < nvar; m++){dataLocal[(n-istart)*nvar+m] = data[n + m*imgpix];}}
/* then send the partials to the other PEs */
    int irank;
    for(irank = 1; irank < mpi_size; irank++)
    {
      int mpi_trgt;
      int ibufsize;
      float *fbufsend;
      mpi_trgt = irank;
#if EQUIAREA
      istart=istartall[mpi_trgt];
      iend=iendall[mpi_trgt];
#else
      para_range(mpi_trgt,nprocs,numpix,&istart,&iend);
#endif
      ibufsize = (iend-istart+1) * nvar;
      fbufsend= (float*)malloc(sizeof(float) * ibufsize);
      for (n = istart ; n < iend+1 ; n++){for (m = 0; m < nvar; m++){fbufsend[(n-istart)*nvar+m] = data[n + m*imgpix];}}
      mpilabel = 1400 + irank;
      MPI_Send(fbufsend, ibufsize, MPI_REAL, mpi_trgt, mpilabel, MPI_COMM_WORLD);
      free(fbufsend);
    }
  }
  else
  {
    int mpi_from = 0;
    int ibufsize;
    float *fbufrecv;
#if EQUIAREA
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    ibufsize = (iend-istart+1) * nvar;
    fbufrecv = (float*)malloc(sizeof(float) * ibufsize);
    mpilabel = 1400 + mpi_rank;
    MPI_Recv(fbufrecv, ibufsize, MPI_REAL, mpi_from, mpilabel, MPI_COMM_WORLD, &mpistat);
    for (n = istart ; n < iend+1 ; n++){for (m = 0; m < nvar; m++){dataLocal[(n-istart)*nvar+m] = fbufrecv[(n-istart)*nvar+m];}}
    free(fbufrecv);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) printf("input data had propagated to all PE.\n");

/* now send partial mask-data to all PE */
  if (mpi_rank == 0)
  {
/* first, the primary makes copy for its own part */
    myrank = mpi_rank;
    nprocs = mpi_size;
    numpix = imgpix;
#if EQUIAREA
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    for (n = istart ; n < iend+1 ; n++){nan_mapLocal[n-istart] = nan_map[n];}
/* then send the partials to the other PEs */
    int irank;
    for(irank = 1; irank < mpi_size; irank++)
    {
      int mpi_trgt;
      int ibufsize;
      int *ibufsend;
      mpi_trgt = irank;
#if EQUIAREA
      istart=istartall[mpi_trgt];
      iend=iendall[mpi_trgt];
#else
      para_range(mpi_trgt,nprocs,numpix,&istart,&iend);
#endif
      ibufsize = (iend-istart+1) * 1;
      ibufsend= (int*)malloc(sizeof(int) * ibufsize);
      for (n = istart ; n < iend+1 ; n++){ibufsend[n-istart] = nan_map[n];}
      mpilabel = 1500 + irank;
      MPI_Send(ibufsend, ibufsize, MPI_INTEGER, mpi_trgt, mpilabel, MPI_COMM_WORLD);
      free(ibufsend);
    }
  }
  else
  {
    int mpi_from = 0;
    int ibufsize;
    int *ibufrecv;
#if EQUIAREA
    istart=istartall[myrank];
    iend=iendall[myrank];
#else
    para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
    ibufsize = (iend-istart+1) * 1;
    ibufrecv = (int*)malloc(sizeof(int) * ibufsize);
    mpilabel = 1500 + mpi_rank;
    MPI_Recv(ibufrecv, ibufsize, MPI_INTEGER, mpi_from, mpilabel, MPI_COMM_WORLD, &mpistat);
    for (n = istart ; n < iend+1 ; n++){nan_mapLocal[n-istart] = ibufrecv[n-istart];}
    free(ibufrecv);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) printf("mask data had propagated to all PE.\n");

/* now do inversion in parallel */
  if (mpi_rank == 0) printf("now doing inversion in parallel.\n");
  myrank = mpi_rank;
  nprocs = mpi_size;
  numpix = imgpix;
#if EQUIAREA
  istart=istartall[myrank];
  iend=iendall[myrank];
#else
  para_range(myrank,nprocs,numpix,&istart,&iend);
#endif
  int pixdone;
  pixdone = 0;
  for (n = istart; n < iend + 1; n++)
  {

    if (nan_mapLocal[n-istart] == 0) 
    {
#if KEIJISKIP
/*   if ((n % rows > rows * 18/32) && (n % rows < rows * 19/32)) */
     if ((n % rows > 2298) && (n % rows < 2400))
     {
#endif
      for(m = 0; m < nvar; m++){obs[m]=dataLocal[(n-istart)*nvar+m];}
      invert_ (obs, scat, guess, res, err);                                  // do inversion
      for (j=0; j<paramct; j++){FinalResLocal[(n-istart)*paramct+j]=res[j];} // copy results to the local-small array(s)
      for (k=0; k<Err_ct;  k++){FinalErrLocal[(n-istart)*Err_ct +k]=err[k];}
#if KEIJISKIP
     }
     else
     {
      for (j=0; j<paramct; j++){FinalResLocal[(n-istart)*paramct+j]=1.0+mpi_rank;} // put dummy
      for (k=0; k<Err_ct ; k++){FinalErrLocal[(n-istart)*Err_ct +k]=0.0;}
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
  printf("Hello, this is %d th PE : done inversion for %d pixels \n", mpi_rank, pixdone);
  MPI_Barrier(MPI_COMM_WORLD);

/* now output data are gathered to primary from each PE in charge of subscription range (from istart to iend) */
  if (mpi_rank == 0)
  {
/* first, copy the portion the primary did */
    myrank = mpi_rank;
    nprocs = mpi_size;
    numpix = imgpix;
#if EQUIAREA
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
/* then collecting the portions done by the other */
    printf("now collecting data from PEs : ");
    int irecv;
    for (irecv = 1; irecv < mpi_size ; irecv++)
    {
      printf(" %d ",irecv);
      int mpi_from;
      nprocs = mpi_size;
      numpix = imgpix;
#if EQUIAREA
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
      mpilabel = 1200 + irecv;
      MPI_Recv(dbufrecv, ibufsize, MPI_DOUBLE, mpi_from, mpilabel, MPI_COMM_WORLD, &mpistat);
      for (n = istart ; n < iend+1 ; n++)
      {
        for (j=0; j<paramct; j++){FinalRes[(n*paramct)+j]=dbufrecv[(n-istart)*(paramct+Err_ct)        +j];}
        for (k=0; k<Err_ct ; k++){FinalErr[(n*Err_ct) +k]=dbufrecv[(n-istart)*(paramct+Err_ct)+paramct+k];}
      }
      free(dbufrecv);
    }
  }
  else
  {
    int isend;
    int mpi_trgt = 0;
    nprocs = mpi_size;
    numpix = imgpix;
    isend = mpi_rank;
#if EQUIAREA
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
    mpilabel = 1200 + mpi_rank;
    MPI_Send(dbufsend, ibufsize, MPI_DOUBLE, mpi_trgt, mpilabel, MPI_COMM_WORLD);
    free(dbufsend);
  }
  MPI_Barrier(MPI_COMM_WORLD); // silly..but always safe

/* PE 0 takes care of pixels none of PE took care of */
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
  }

  MPI_Barrier(MPI_COMM_WORLD); // silly..but always safe
  free(FinalResLocal);         // liberate....
  free(FinalErrLocal);
  free(dataLocal);
  free(nan_mapLocal);
  MPI_Barrier(MPI_COMM_WORLD); // silly..but always safe

/* write outputs through DRMS/JSOC */

  if (mpi_rank == 0)
  {
        /*  output of results to individual FITS files for each parameter  */
    outRec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    printf("created record\n");
    if (!outRec) {fprintf (stderr, "Error creating record in series %s; abandoned\n",outser);return 1;}

    for (j = 0; j < paramct; j++)
    {
      float *dat1;
      dat1 = (float*) calloc(imgpix, sizeof(float));
      /*  should get this info from segment lookup in output record  */
      int axes[2];
      axes[0] = cols;
      axes[1] = rows;
      for(n = 0; n < imgpix ; n++){dat1[n] = FinalRes[(n*paramct)+j];}  /* collect the result for each parameter over all pixels */
      invrt_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, dat1, &status);
      seg = drms_segment_lookup (outRec, Resname[j]);
      if (!seg) {fprintf (stderr, "Error getting data segment %s; abandoned\n", Resname[j]);}
      if (drms_segment_write (seg, invrt_array, 0))
      {
        fprintf (stderr, "Error writing segment %d (%s); abandoned\n", j,seg->info->name);
        return 1;
      }
      else
      {
        printf ("Results written out to %15s\n", Resname[j]);
      }
      free(dat1);
    }

    printf("printing errors now\n");
    for (k = 0; k  < Err_ct; k++)
    {
      float *dat2;
      dat2 = (float*) calloc(imgpix , sizeof(float));
        /*  should get this info from segment lookup in output record  */
      int axes[2];
      axes[0] = cols;
      axes[1] = rows;
      for (n=0; n < imgpix;n++){dat2[n] = FinalErr[(n*Err_ct)+k];}  /* collect the result for each parameter over all pixels */
      err_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, dat2,&status);
      seg = drms_segment_lookup (outRec, Resname[k+paramct]);
      if (!seg)
      {
        fprintf (stderr, "Error getting data segment %s; abandoned\n", Resname[k+paramct]);
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
          printf ("errors written out to %15s\n", Resname[k+paramct]);
        }
      }
        drms_copykey(outRec, inRec, "CENTER_X");
        drms_copykey(outRec, inRec, "CENTER_Y");
        drms_copykey(outRec, inRec, "SOLAR_B0");
        drms_copykey(outRec, inRec, "SOLAR_P0");
        drms_copykey(outRec, inRec, "R_SUN");
        drms_copykey(outRec, inRec, "TELESCOP");
        drms_copykey(outRec, inRec, "INSTRUME");
        drms_copykey(outRec, inRec, "WAVELENG");  
   
      free(dat2);
    }
    drms_close_record (outRec, DRMS_INSERT_RECORD);

    time (&end);
    printf ("%ld sec for %d profiles\n", end-start, imgpix);
    printf ("%.2f profiles per second\n", (float) npix*segct / (float)(0.01 + end-start));

    drms_close_records (records, DRMS_FREE_RECORD);
  } // end-if mpi_rank is 0 nor not.

  MPI_Barrier(MPI_COMM_WORLD);
/* say good bye to MPI things.*/
  MPI_Finalize();

  return 0;
}

/* ------------------------------- end of main-wrapper layer ------------------------------- */


/* added by K.H. */

/* --------------------------------------------------------------------------------
 *
 * Calculate the first and the last addresses of array, each PE is in charge of. 
 *
 * -------------------------------------------------------------------------------- */

void para_range(int myrank, int nprocs, int numpix, int *istart, int *iend)
{
  int iwork1, iwork2, imin, idummy, n1, n2;
  n1 = 0;              // When written in Fortran, this is typically 1.
  n2 = numpix - 1;     // When written in Fortran, this is typically numpix.
#if 1
  iwork1 = (n2 - n1 + 1) / nprocs;
  iwork2 = (n2 - n1 + 1) % nprocs; // mod(n2-n1+1,nprocs)
  if (iwork2 >= myrank){imin = myrank;}else{imin=iwork2;} // imin = min(myrank,iwork2)
  *istart = myrank*iwork1 + n1 + imin;
  idummy  = *istart + iwork1 - 1;
  if (iwork2 <= myrank){*iend = idummy;}else{*iend=idummy+1;}
#endif
/* equi-cosine .... does not work well .... Ummmm */
#if 0
  float fwork1, fwork2;
  fwork1 = 1.0 - 2.0 * ((float) myrank      ) / ((float) nprocs);
  fwork2 = 1.0 - 2.0 * ((float) myrank + 1.0) / ((float) nprocs);
  fwork1 = acos(fwork1) / 3.14159265358979 * 2.0 * ((float) numpix);
  fwork2 = acos(fwork2) / 3.14159265358979 * 2.0 * ((float) numpix);
  iwork1 = (int) fwork1 - 1;
  iwork2 = (int) fwork2 + 1;
  if (iwork1 < n1) iwork1 = n1;
  if (iwork2 > n2) iwork2 = n2;
  *istart = iwork1;
  *iend   = iwork2;
#endif
}

/* end of a part added by K.H. */



/* ----------------------------- end of this file ----------------------------- */
