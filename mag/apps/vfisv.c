 /*
 *  vfisv.c						~rick/hmi/vecb/src/v09
 *
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
 *	should be eliminated from final version, as they add unnecessray
 *	complexity
 *    The statistics for inverted quantities are calculated over the whole
 *	image, not just the part inverted
 *    Off-limb pixels are not handled well, the inversion results are set
 *	to 0; this means that the statistics of the inverted parameters
 *	calculated over the whole image are meaningless.
 *    The initialization of the inversion quantities is silly, an arbitrary
 *	one of the Stokes parameters; they should be calloc'd to 0 or (better)
 *	NaN
 *    1-21-10:Added a map_param which will track which of the datapts were NaN's
 *      and were then converted to Zero. We should not borther to waste computational
 *      power to invert them.
 *    There is no parameter to govern which of the quantities are inverted;
 *	it has to be modified by changing the FREE array in map_param.f90
 *	and recompiling, and output files are generated for all quantities,
 *	regardless of whether they have been calculated or not
 *    Likewise, the QUICKLOOK, STRAYLIGHT_CALCULATION, and ERRORS parameters
 *	are set in code (map_param.f90)
 *    The output FITS files are floats, not scaled shorts; an entire set of
 *	18*imgpix values have to be allocated for assembly; it might be better
 *	to write slices into the output FITS segments.
 *    Can't handle case in which some data segments (e.g. a wavlength set)
 *	are missing.
 *    Number of threads must be a divisor of number of pixels in image;
 *	2**n for full HMI images; need to generalize, and allow for unequal
 *	thread segments
 *
 */

char *module_name = "vfisv";
char *version_id = "0.9";

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

#define ONECOL	(spct*wlct)
#define ONEROW	(spct*wlct*cols)
#define IMGCTR	(cols*(rows/2) + cols/2)

ModuleArgs_t module_args[] = {
  {ARG_STRING,  "in", "su_rsb.stokes_profile_yl[:#5]", "input series"},  /*This is the input data in the form of 24 fits file segnments */
  {ARG_STRING,  "out", "su_priya.Res_Err_stokes_inv", "output series"},
  {ARG_INT,  "npix", "0", "number of pixels per segment (default: all/segs)"},
  {ARG_INT,  "num_iter", "30", "number of iterations(default: 30)"},
  {ARG_INT,  "num_lambda", "33", "number of ??(default: 33)"},
  {ARG_INT,  "Num_lambda_filter", "6", "Number of filters accross the wvl (default: 6)"},
  {ARG_INT,  "Num_tunning", "6", "Number of ??(default: 6)"},
  {ARG_DOUBLE,  "svd_tolerance", "1.0e-32", "svd tolerance (default: 1.0e-32)"},
  {ARG_DOUBLE,  "chi2_stop", "1.0e-6", "chisq-stop (default: 1.0e-6)"},
  {ARG_DOUBLE,  "Polarization_threshold", "1.0e-2", "polarization threshold (default: 0.01)"},
  {ARG_DOUBLE,  "Intensity_Threshold", "0.8", "Intensity threshold (default: 0.8)"},
  {ARG_DOUBLE,  "Percentage_Jump", "10.0", "Percentage Jump (default: 10%)"},
  {ARG_DOUBLE,  "Lambda_Min", "-432.0", "Intensity threshold (default: -432)"},
  {ARG_DOUBLE,  "Lambda_0", "6173.3433", "Wavelength(default:6173.3433 Angstrom )"},
  {ARG_DOUBLE,  "Lambda_B", "0.044475775", "FWHM?? (default: 0.044475775)"},
  {ARG_DOUBLE,  "Delta_Lambda", "27.0", "Delta Lambda(default: 27.0)"},
  {ARG_DOUBLE,  "Noise_LEVEL", "3.0e-3", "Intensity threshold (default: 3.0e-3)"},
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
  DRMS_Record_t *rec;
  DRMS_Segment_t *seg;
  DRMS_Array_t *stokes_array, *invrt_array,*err_array;
  
  time_t start, end;
  double *ddat, *res;
  double *obs, *scat,*err;
  double guess[10]= {15.0,90.0,45.0,0.05,40.0,50.0,0.0,0.4,0.6,1.0};
  double *FinalRes,*FinalErr;
  double SVD_TOLERANCE, CHI2_STOP, POLARIZATION_THRESHOLD, INTENSITY_THRESHOLD, PERCENTAGE_JUMP;
  double LAMBDA_0, LAMBDA_B, NOISE_LEVEL;
  double LAMBDA_MIN, DELTA_LAMBDA;
  double LYOTFWHM, WNARROW, WSPACING;
  double sumchi, avgval, minval, maxval;
  float *data, *data0, *inv = NULL;
  int *map_nan;
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
  int segct, mpi_rank, mpi_size;
  int NUM_ITERATIONS;
  int NUM_LAMBDA;
  int NUM_LAMBDA_FILTER, NUM_TUNNING, CONTINUUM;
  char *Resname[] = {"eta_0", "inclination", "azimuth", "damping", "dop_width",
      "field", "vlos_mag", "src_continuum", "src_grad", "alpha_mag","field_err","inclination_err","azimuth_err", "vlos_err", "alpha_err","field_inclination_err","field_az_err", "inclin_azimuth_err" , "field_alpha_err", "inclination_alpha_err", "azimuth_alpha_err" };
     
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
  							/*  parameters  */
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
 

					/*  Initialize Clock  */
    time (&start);
    printf("Lambda_O= %f\n",LAMBDA_0);

  MPI_Status mpi_stat;   /* Initalizing MPI */

  status = 0;
  MPI_Init (&status, NULL);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  segct = mpi_size;
  if (mpi_rank == 0) {
    if (verbose) printf ("%s V %s\n", module_name, version_id);
						/*  open input record_set  */
    records = drms_open_records (drms_env, indsdesc, &status);
    if (status) {
      fprintf (stderr, "drms_open_records failed, in=%s, status=%d.  Aborting.\n",
	indsdesc, status);
      return 1;
    }
    if ((rec_ct = records->n) == 0) {
      fprintf (stderr, "No records in selected dataset %s  Aborting.\n",
	  indsdesc);
      return 1;
    }
    if ((rec_ct = records->n) > 1) {
      fprintf (stderr, "Warning: only first record in selected set processed\n");
    }
    rn = 0;

    records = drms_open_records (drms_env, indsdesc, &status);
    if (status) {
      fprintf (stderr, "drms_open_records failed, in=%s, status=%d.  Aborting.\n",
	  indsdesc, status);
      return 1;
    }
    if ((rec_ct = records->n) == 0) {
      fprintf (stderr, "No records in selected dataset %s  Aborting.\n",
	  indsdesc);
      return 1;
    }
    if ((rec_ct = records->n) > 1) {
      fprintf (stderr, "Warning: only first record in selected set processed\n");
    }
    rn = 0;
    rec = records->records[rn];
    seg = drms_segment_lookupnum (rec, 0);
    cols = seg->axis[0];
    rows = seg->axis[1];
    imgpix = cols * rows;
    imgbytes = imgpix * sizeof (float);
    nvar = wlct * spct;

    if (npix < 1) npix = imgpix / segct;
    if (imgpix % npix) {
      printf ("Error: npix (%d) must be a divisor of %d\n", npix, imgpix);
      return 0;
    }

    data = data0 = (float *)malloc (sizeof (float) * imgpix * nvar);
    map_nan=(int*)calloc (imgpix, sizeof (int));

    for (sp = 0; sp < spct; sp++) {      /*spct=4,wlct=6    */
      for (wl = 0; wl < wlct; wl++) {
	sprintf (segname, "%s%d", spname[sp], wl);
        if ((seg = drms_segment_lookup (rec, segname)) == NULL) {
	  fprintf (stderr, "Error reading segment %s of record %d\n",
	      segname, rn);
	  return 1;
	}
	stokes_array = drms_segment_read (seg, DRMS_TYPE_FLOAT, &status);/* data is read in 1 segment at a time-so for each location;  */
	/* printf("segment read %s\n",segname); */                	/*  the 24 requiired values are 4096*4096 datapoints apart   */
	memcpy (data, stokes_array->data, imgbytes);
	drms_free_array (stokes_array);
	data += imgpix;
      }
    }
    data = data0;
    printf("Imgpix= %d\n",imgpix);
    /*Map of pixels where input data are Nan */ for (n = 0; n < imgpix; n++) {
         for (m = 0; m < nvar; m++) {
        if (isnan (data[n + m*imgpix])) {
          for (m = 0; m < nvar; m++) data[n + m*imgpix] = 0.0;
	  map_nan[n]=1; 
          
           }
      }
 }
  }
    printf("data is read\n");
  segsize = wlct * spct * npix;
  ressize = paramct * npix;
  obs = (double *)malloc (sizeof (double) * nvar);
  res = (double *)calloc (paramct, sizeof (double));
 scat=(double *)malloc (sizeof (double) * nvar);
 err=(double *)calloc (8,sizeof (double));

 FinalErr=(double *)malloc(sizeof(double)*imgpix*10);
 FinalRes=(double *)malloc(sizeof(double)*imgpix*paramct);
							
  line_init_(&LAMBDA_0, &LAMBDA_B,&NOISE_LEVEL);
  printf("done line_init for  mpirank %d\n",mpi_rank);
  wave_init_ (&LAMBDA_MIN, &DELTA_LAMBDA, &NUM_LAMBDA);
  printf("done wave_init  for  mpirank %d\n", mpi_rank );
  filt_init_ (&NUM_LAMBDA_FILTER, &NUM_TUNNING, &CONTINUUM, &LYOTFWHM, &WNARROW, &WSPACING);
  printf("done filt_init for mpi_rank %d\n",mpi_rank);
   inv_init_(&NUM_ITERATIONS, &SVD_TOLERANCE, &CHI2_STOP, &POLARIZATION_THRESHOLD, &INTENSITY_THRESHOLD, &PERCENTAGE_JUMP);  printf("done inv_init\n");
  free_init_(list_free_params);
  printf("done list_free_params for mpi_rank %d\n", mpi_rank);
  svd_init_();
  printf("done svd_init\n  ");
  if  (list_free_params[4] == 0.0) voigt_init_();  
  for (n=0; n< nvar; n++) scat[n]=0;
  
  /*for (n=800000 ;n<800020;n++) {
    for (i=0;i<nvar;i++){
              printf("For imgpix=%d,  Obs= %9.6f",n,data[i+i*imgpix]) ;
         }
	 }*/

  for (n=0; n<16777216 ; n++) {

         for(m= 0; m<nvar; m++) {
           obs[m]=data[n+(m*imgpix)];
          }
          
           invert_ (obs, scat, guess, res, err);
	   /*for (i=0;i<nvar;i++)  printf("For imgpix=%d,  Obs= %9.6f",n,obs[i]) ;  */
            for (j=0; j<paramct; j++) FinalRes[(n*paramct)+j]=res[j];
           for (k=0; k<8; k++) FinalErr[(n*8)+k]=err[k];    
	   /*   for (j=0;j<paramct;j++)  printf("n=%d,result:%9.6f\n",n,res[j]);  */
  }

  /*
  for ( n=800000; n<800025; n++) {
       for (j=0;j<10; j++) { 
	 printf("n=%d,j=%d, Finalres=%9.6f\n",n,j,FinalRes[(n*paramct)+j]);
   }  
   } 
  */
 printf("done inverting-now to combine and write out\n");
     MPI_Finalize ();
  if (mpi_rank == 0) {
	/*  output of results to individual FITS files for each parameter  */
    rec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    printf("created record\n");
    if (!rec) {fprintf (stderr, "Error creating record in series %s; abandoned\n",
	  outser);
      return 1;
    }
    

    for (j = 0; j < paramct; j++) { 
      float *dat1;
      dat1 = (float*) calloc(imgpix, sizeof(float));
      /*  should get this info from segment lookup in output record  */ 
     int axes[2];
      axes[0] = cols;
      axes[1] = rows;
     for(n=0; n < imgpix;n++)
       {
        /* collect the result for each parameter over all pixels */
        dat1[n] = FinalRes[(n*paramct)+j];
       }
        invrt_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, dat1, &status);
       seg = drms_segment_lookup (rec, Resname[j]);
       
      if (!seg) {
	    fprintf (stderr, "Error getting data segment %s; abandoned\n", Resname[j]);	   
      }
      if (drms_segment_write (seg, invrt_array, 0)) {
        fprintf (stderr, "Error writing segment %d (%s); abandoned\n", j,
	    seg->info->name);
	    return 1;
      }else
     printf ("Results written out to %15s\n", Resname[j]);
    }
    
    printf("printing errors now\n");
  for (k = 1; k  < 8; k++) {
   float *dat2;
   dat2 = (float*) calloc(imgpix , sizeof(float));
	/*  should get this info from segment lookup in output record  */
   int axes[2];
      axes[0] = cols;
      axes[1] = rows;
	  for (n=0; n < imgpix;n++)
           {
        /* collect the result for each parameter over all pixels */
          dat2[n] = FinalErr[(n*8)+k];
            }

      err_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, dat2,&status);
      seg = drms_segment_lookup (rec, Resname[k+paramct]);
      if (!seg) {
	   fprintf (stderr, "Error getting data segment %s; abandoned\n", Resname[k+paramct]);
      } else
      if (drms_segment_write (seg, err_array, 0)) {
        fprintf (stderr, "Error writing to segment %d (%s); abandoned\n", k+paramct,
	    seg->info->name);
	    return 1;
      }else
     printf ("errors written out to %15s\n", Resname[k+paramct]);
    }

    drms_close_record (rec, DRMS_INSERT_RECORD);

    time (&end);
    printf ("%ld sec for %d profiles\n", end-start, imgpix);
    printf ("%.2f profiles per second\n",
        (float) npix*segct / (float)(0.01 + end-start));

    /*printf ("\n                   mean         min         max      valid\n");
      for (m = 0; m < paramct; m++) {
      int vct = 0;
      avgval = 0.0;
      minval = maxval = inv[m * imgpix];
      for (j = 0; j < segct; j++) {
	int soff = npix * j + m * imgpix;
	for (n = 0; n < npix; n++) {
	   if (isnan (fval = inv[soff + n])) continue;
	   vct++;
	   avgval += fval;
	   if (minval > fval) minval = fval;
           if (maxval < fval) maxval = fval;
	}
      }
      if (vct) printf ("%15s %11.4e %11.4e %11.4e %8d\n", Resname[m],
	  avgval / vct, minval, maxval, vct);
      else
 	printf ("%15s    -----       -----       -----    %8d\n", Resname[m], vct);
	}*/
    drms_close_records (records, DRMS_FREE_RECORD);
  }
  return 0;
}
