int cosmic_rays(DRMS_Record_t *record, int *image_int, int *badpix, int nbad, const float limit, unsigned char *cosmic,  int *n_cosmic, int nx, int ny)
{
  //fresize: structure obtained from init_fresize_gaussian(&fresizes,stdd,nwd,1);
  //image: input image
  //limit: detection limit in sigma (somewhere between 4.0 and 7.0)
  //cosmic: integer array of image size: will be used to write cosmic ray coordinates into [0...n_cosmic-1]
  //n_cosmic: number of cosmic ray hits
  //nx: horizontal dimension
  //ny: vertical dimension
  //

  struct fresize_struct fresizes;

  //constants
   //initialize
  static int first_time=1;

  if (first_time)
    {
      first_time=0;
      float stdd=fwhm/2.0/sqrt(2.0*log(2.0)); 
      int nwd=(int)(fwhm*2.0);  //kernel size
      init_fresize_gaussian(&fresizes,stdd,nwd,1);
    }


  int i, j;
  int count;
  float *imhp;
  imhp=(float *)(MKL_malloc(nx*ny*sizeof(float),malign));
  float sigma;
  float sum;
  int nump=nx/cent_frac*ny/cent_frac;
  int status=1;

  float *image;
  image=(float *)(MKL_malloc(nx*ny*sizeof(float),malign));

  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) if (image_int[j*nx+i] == nanval) image[j*nx+i]=NAN; else image[j*nx+i]=(float)image_int[j*nx+i];

  for (i=0; i<nbad; ++i)image[badpix[i]]=NAN; // set bad pixels to NAN

  
   fresize(&fresizes,image, imhp, nx,ny,nx,nx,ny,nx,0,0,NAN);   //convolve with gaussian
  
  #pragma omp parallel for private(i,j)
  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) imhp[j*nx+i]=image[j*nx+i]-imhp[j*nx+i]; //highpass filter

  //get stddev at center
  //*******************************************************************************************************
  sum=0.0;
  count=0;

  //calculate stddev at center 
#pragma omp parallel for reduction(+:sum,count) private(i,j)  
  for (j=ny*(cent_frac-1)/cent_frac/2; j<ny*(cent_frac+1)/cent_frac/2; ++j) for (i=nx*(cent_frac-1)/cent_frac/2; i<nx*(cent_frac+1)/cent_frac/2; ++i)
    {
      if (!isnan(imhp[j*nx+i]))
	{
	  sum += imhp[j*nx+i]*imhp[j*nx+i];
	  ++count;
	}
    }

  //detect everything that stick out more than limit*sigma
  if (count > nump/2) sigma=sqrt(sum/(float)(count)); else sigma=NAN;
  //********************************************************************************************************

  //check for outliers
  //*******************************************************************************************************
  count=0;
  if (!isnan(sigma) && sigma > sigmamin && sigma < sigmamax)
    {
      #pragma omp parallel for reduction(+:count) private(i,j)
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
	{
	  if (!isnan(imhp[j*nx+i]))
	    {
	      if (imhp[j*nx+i] > limit*sigma)
		{
		  cosmic[j*nx+i]=1;
		  ++count;
		}
	      else
		{
		  cosmic[j*nx+i]=0;
		}
	    }
	}
      status=0;
    }
  //********************************************************************************************************

  *n_cosmic=count;

    
  MKL_free(imhp);
  MKL_free(image);

  return status;
}
