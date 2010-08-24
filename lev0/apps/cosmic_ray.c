int cosmic_rays(DRMS_Record_t *record, float *image, int *badpix, int nbad, int *cosmic, int *n_cosmic, int nx, int ny)
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


  float coef[5];


      float stdd=fwhm/2.0/sqrt(2.0*log(2.0)); 
      int nwd=(int)(fwhm*kernel_size);  //kernel size
      init_fresize_gaussian(&fresizes,stdd,nwd,1);
      if (fresize == NULL){printf("can not initialize fresize\n"); exit(EXIT_FAILURE);}

      coef[0]=1.00000;
      coef[1]=0.0;
      coef[2]=0.0;
      coef[3]=0.0;
      coef[4]=0.0;

    
  int i, j;
  int count;

  float *imhp;
  imhp=(float *)(malloc(nx*ny*sizeof(float)));
  if (imhp == NULL){printf("can not allocate memory\n"); exit(EXIT_FAILURE);}

  float *img;
  img=(float *)(malloc(nx*ny*sizeof(float)));
  if (img == NULL){printf("can not allocate memory\n"); exit(EXIT_FAILURE);}

  float sigma;
  float sum;
  float sig, r, factor;
  float x0, y0, rad;
  float rsun_obs, image_scl;
  int nump=nx/cent_frac*ny/cent_frac;
  int status=1, status1, status2;
  int stat_pos=0;
  int focid, camid;


  
  x0=drms_getkey_float(record,X0_MP_key,&status);
  if (status != 0 || isnan(x0)){stat_pos=1; printf("no center\n");}

  y0=drms_getkey_float(record,Y0_MP_key,&status);
  if (status != 0 || isnan(y0)){stat_pos=1;printf("no center\n");}

  rsun_obs=drms_getkey_double(record,RSUN_OBS_key,&status1);
  image_scl=drms_getkey_float(record,IMSCL_MP_key,&status2);
  if (status1 == 0 && status2 == 0 && !isnan(image_scl) && image_scl > 0.0 && !isnan(rsun_obs) && rsun_obs > 0.0){rad=rsun_obs/image_scl;} else {stat_pos=1;printf("no radius\n"); }

  focid=drms_getkey_int(record,HCFTID_key,&status1);
  camid=drms_getkey_int(record,HCAMID_key,&status1);

  if (status1 != 0 || status2 != 0){printf("no foc id or camera\n"); exit(EXIT_FAILURE);}

  if ((camid != light_val1 && camid != light_val2) || focid < 1 || focid > 16){stat_pos=1; factor=2.0;} else factor=1;


  //copy image
  #pragma omp parallel for private(i,j)
  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) img[j*nx+i]=image[j*nx+i];

  for (i=0; i<nbad; ++i)if (badpix[i] < 0 || badpix[i] >=  16777216){printf("invalid pixel coordinate\n"); exit(EXIT_FAILURE);}
  for (i=0; i<nbad; ++i){img[badpix[i]]=NAN; cosmic[i]=badpix[i];} // set bad pixels to NAN
 
  
 


   fresize(&fresizes,img, imhp, nx,ny,nx,nx,ny,nx,0,0,NAN);   //convolve with gaussian



  #pragma omp parallel for private(i,j)
  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) imhp[j*nx+i]=img[j*nx+i]-imhp[j*nx+i]; //highpass filter

  //get stddev at center
  //*******************************************************************************************************
  sum=0.0;
  count=0;



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
  printf("sigma %f\n", sigma);
  //********************************************************************************************************

  //check for outliers
  //*******************************************************************************************************
  count=nbad;
 
 
  if (!isnan(sigma) && sigma > sigmamin && sigma < sigmamax)
    {
#pragma omp parallel for private(i,j,r)
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
	{

	  if (!isnan(imhp[j*nx+i]))
	    {
	      if (stat_pos == 0) r=sqrt(pow((float)j-y0,2)+pow((float)i-x0,2))/rad; else r=0.0;
	      if (r< 0.98){
		sig=(coef[0]+coef[1]*r+coef[2]*pow(r,2)+coef[3]*pow(r,3)+coef[4]*pow(r,4))*sigma;
		if (imhp[j*nx+i] > limit*factor*sig || image[j*nx+i] > maxval)
		{
                  #pragma omp critical(detect)
		  cosmic[count]=j*nx+i;
		  ++count;
		}
	      }
	    }
	    
	    
	}
      status=0;
     }
  //********************************************************************************************************

  *n_cosmic=count;

    

 
	
  free_fresize(&fresizes);
  free(imhp);
  free(img);




  return status;

}
