//Return codes:
//2: center fit and radius fit failed
//1: only radius fit failed
//0: success 
int limb_fit(DRMS_Record_t *record, float *image_in, double *rsun_lf, double *x0_lf, double *y0_lf, int nx, int ny, int method)
{

void cross_corr(int nx, int ny, double *a, double *b);
int is_parab(double *arr, int n, int nmax);
void free_mem(struct mempointer *memory);


 int status_res, status_gap;
  int status,status1,status2;

  //initialization
  int i, j;

  
  struct fresize_struct fresizes;

  struct fill_struct fills;

  static int config_id=0;
  static unsigned char mask[4096*4096];


 


  size_t bytes_read;
  int new_config_id;
  new_config_id=drms_getkey_int(record, "HIMGCFID", &status);

  FILE *maskfile;
 
  if (status == DRMS_SUCCESS && new_config_id >= min_imcnf && new_config_id <=max_imcnf)
    {
    if (config_id != new_config_id)
      {

	unsigned char *mk;
	mk=(unsigned char *)(malloc(sizeof(unsigned char)));

	config_id=new_config_id;

	printf("new image config\n");
	maskfile=fopen("/home/jsoc/hmi/tables/cropmask.105", "rb"); //use the same crop table for now
    if (maskfile==NULL){fputs("mask file not found", stderr); status_res=2; return status;}
    for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) 
      {
	bytes_read=fread(mk,sizeof(unsigned char),1,maskfile);
	mask[j*nx+i]=mk[0];
      }

    fclose(maskfile);
      }
    }
  else 
    {
      status_res=2; return status_res;
    }

    ////////////////////////////




  //image: input image
  //nx: horizontal dimension
  //ny: vertical dimension

  //constants 
  const int nxx=nx/4;
  const int nyy=nx/4;

  const int binx=nx/nxx;
  const int biny=nx/nyy;

  int k,l;
  
   long double akl[3][3];
   long double bl[3];
   long double cof[3];





  int rcount=0;  //counter for missing points

  double dx, ds;

  int ispar;
  int camid, focid;
  double cx, cy, rad;
  double image_scl;
  double rsun_obs;


  cx=drms_getkey_float(record,X0_MP_key,&status);
  if (status != 0 || isnan(cx)){status_res=2; return status_res;}
  *x0_lf=(double)cx;

  cy=drms_getkey_float(record,Y0_MP_key,&status);
  if (status != 0 || isnan(cx)){status_res=2; return status_res;}
  *y0_lf=(double)cy;

  rsun_obs=drms_getkey_double(record,RSUN_OBS_key,&status1);
  image_scl=drms_getkey_float(record,IMSCL_MP_key,&status2);



  if (status1 == 0 && status2 == 0 && !isnan(image_scl) && image_scl > 0.0 && !isnan(rsun_obs) && rsun_obs > 0.0){rad=(double)rsun_obs/image_scl;} else {status_res=2; return status_res;}
  *rsun_lf=(double)(rad);

 

  if (nx != ny){printf("Limb finder does not work for non-square images\n"); status_res=2; return status_res;}

  camid=drms_getkey_int(record,HCAMID_key,&status);
  focid=drms_getkey_int(record,HCFTID_key,&status);
  if ((camid != light_val1 && camid != light_val2) || focid < 1 || focid > 16){printf("not an obs image\n"); status_res=2; return status_res;}



  int nr=(int)(rad*(high-low)/2.0)*2+1; //odd number of points
  int nphi=nx/2;

  int marg_up=(int)((high-1.0)*rad+cent_err);
  int marg_lo=(int)(-(low-1.0)*rad+cent_err); // make wider margin. Important if initial center estimate is off!!
  int marg_up0=marg_up/binx-1;
  int marg_lo0=marg_lo/binx-1;
 
  //       
  struct mempointer memory;

  double *xrp, *yrp, *imrphi, *rc, *phic;

  xrp=(double *)(malloc(nr*nphi*sizeof(double)));
  yrp=(double *)(malloc(nr*nphi*sizeof(double)));
  imrphi=(double *)(calloc(nr*nphi,sizeof(double)));
  rc=(double *)(malloc(nr*sizeof(double)));
  phic=(double *)(malloc(nphi*sizeof(double)));

  double *avgphi=(double *)(calloc(nr, sizeof(double)));

 
  double *imhp=(double *)(calloc(nxx*nyy, sizeof(double)));
  double *imro=(double *)(calloc(nxx*nyy, sizeof(double)));
  float *image=(float *)(malloc(nxx*nyy*sizeof(float)));
  double *parab=(double *)(malloc(parsize*sizeof(double)));
  unsigned char *mask_p=(unsigned char *)(malloc(nx*nx*sizeof(unsigned char)));
  float *cnorm=(float *)(malloc(nx*nx*sizeof(float)));
  float *ierror=(float *)(malloc(nx*nx*sizeof(float)));

  memory.xrp=xrp;
  memory.yrp=yrp;
  memory.imrphi=imrphi;
  memory.rc=rc;
  memory.phic=phic;
  memory.avgphi=avgphi;
  memory.imhp=imhp;
  memory.imro=imro;
  memory.image=image;
  memory.parab=parab;
  memory.mask_p=mask_p;
  memory.cnorm=cnorm;
  memory.ierror=ierror;

  int imgid;

  imgid=drms_getkey_int(record,"HIMGCFID",&status);
 
  rad=rsun_obs/image_scl;
  double rad0=rad/(double)binx;
  double cx0=cx/(double)binx;
  double cy0=cy/(double)biny;
 
  rcount=0;
#pragma omp parallel  for reduction(+:rcount) private(i,j,dx)    
  for (j=0; j<ny; ++j) 
      for (i=0; i<nx; ++i)
	{
	  if (isnan(image_in[j*nx+i])){mask_p[j*nx+i]=2;} else {mask_p[j*nx+i]=0;}
	  dx=sqrt(pow((double)i-cx,2)+pow((double)j-cy,2));
	  if (dx > (rad-(marg_lo+2)) && dx < (rad+(marg_up+2)) && mask[j*nx+i] == 1)
	    {
	      if (isnan(image_in[j*nx+i])){mask_p[j*nx+i]=1; ++rcount;}
	    }
	}

 


  //fill image if there are missing points

      if (rcount > 0)
   	{
	  printf("gapfilling: pixels to fill: %d\n", rcount);	
	  status_gap=init_fill(gapfill_method, gapfill_regular, gapfill_order,gapfill_order2,gapfill_order2,&fills, NULL);
	  if (status_gap == 0){fgap_fill(&fills,image_in,nx,nx,nx,mask_p,cnorm,ierror);} else {status_res=2; free_mem(&memory); return status_res;}
	  free_fill(&fills);
	}



      init_fresize_boxcar(&fresizes,2,4);

      fresize(&fresizes,image_in, image, nx,nx,nx,nxx,nyy,nxx,0,0,NAN);   //rebin image
   
      free_fresize(&fresizes);
       
       

#pragma omp parallel for private(i,j,dx)
         for (j=1; j<(nyy-1); ++j) 
            for (i=1; i<(nxx-1); ++i)
	      {
	 	  dx=sqrt(pow((double)i-cx0,2)+pow((double)j-cy0,2));
		  if (dx > (rad0-marg_lo0) && dx < (rad0+marg_up0))
		    {
		      if(!isnan(image[j*nxx+i-1]) && !isnan(image[j*nxx+i+1]) && !isnan(image[(j+1)*nxx+i]) && !isnan(image[(j-1)*nxx+i]))
			{
			  imhp[j*nxx+i]=sqrt(pow((double)image[j*nxx+i-1]-(double)image[j*nxx+i+1],2)+pow((double)image[(j+1)*nxx+i]-(double)image[(j-1)*nxx+i],2));
			  imro[(nyy-j-1)*nxx+nxx-1-i]=imhp[j*nxx+i];
		 
			}
		    }
		 
	      }

		 


      //calculate center
      //////////////////////////////////////////////////////////////////////////
      
	 int max, xmax, ymax;
	 double fmax=0.0, favg=0.0, denom;
	 double xfra, yfra;


	
 

      cross_corr(nxx,nyy, imhp, imro);


      for (j=0; j<nyy; ++j)
	for (i=0; i<nxx; ++i)
	  {
	  if (imro[j*nxx+i] > fmax){max=j*nxx+i; fmax=imro[j*nxx+i];}
	  favg=favg+imro[j*nxx+i];
	  }


  
      if (favg/(double)(nxx*nyy)*limit_cc< fmax)
	{ 
      xmax=max % nxx;
      ymax=max/nxx;

      if (xmax*ymax > 0 && xmax < (nxx-1) && ymax < (nyy-1))
	{
	  denom = fmax*2.0 - imro[ymax*nxx+xmax-1] - imro[ymax*nxx+xmax+1];
	  xfra = ((double)xmax-0.5) + (fmax-imro[ymax*nxx+xmax-1])/denom;
	  denom = fmax*2.0 - imro[(ymax-1)*nxx+xmax] - imro[(ymax+1)*nxx+xmax];
	  yfra = ((double)ymax-0.5) + (fmax-imro[(ymax-1)*nxx+xmax])/denom;
        }

      cx0=(double)xfra/2.0+(double)nxx/4.0-0.5;
      cy0=(double)yfra/2.0+(double)nyy/4.0-0.5;

      *x0_lf=(double)cx0*(double)binx;
      *y0_lf=(double)cy0*(double)biny;
	}
      else
	{status_res=2; free_mem(&memory); return status_res;}
	 
   
      //////////////////////////////////////////////////////////

      //remap for radius determination
   

   

      for (i=0; i<nr; ++i)rc[i]=((high-low)*(double)i/(double)(nr-1)+low)*rad;
      for (j=0; j<nphi; ++j) phic[j]=(double)j/(double)nphi*2.0f*(double)M_PI;

#pragma omp parallel for private(i,j)
      for (j=0; j<nphi; ++j)
	for (i=0; i<nr; ++i)
	  {
	    xrp[j*nr+i]=rc[i]/4.0*cos(phic[j])+(*x0_lf)/4.0-0.5;
	    yrp[j*nr+i]=rc[i]/4.0*sin(phic[j])+(*y0_lf)/4.0-0.5;
	  }

   

      //remap

      int xl,xu,yl,yu;
      double rx, ry;
      double wg[2*lim+1];

#pragma omp parallel for private(i,j,xl,yl,xu,yu,rx,ry)
      for (j=0; j<nphi; ++j)
	for (i=0; i<nr; ++i)
	  {
	    xl=floor(xrp[j*nr+i]);
	    xu=xl+1;

	    yl=floor(yrp[j*nr+i]);
	    yu=yl+1;

	    rx=(double)xu-(xrp[j*nr+i]);
	    ry=(double)yu-(yrp[j*nr+i]);

	    if (xl >= 0 && xu < nxx && yl >= 0 && yu < nyy)
	    imrphi[j*nr+i]=rx*ry*imhp[yl*nxx+xl]+(1.0-rx)*ry*imhp[yl*nxx+xu]+rx*(1.0-ry)*imhp[yu*nxx+xl]+(1.0-rx)*(1.0-ry)*imhp[yu*nxx+xu];
	  }
   
 	
      fmax=0.0; max=0;
   
      for (i=0; i<(2*lim+1); ++i){wg[i]=1.0-pow(sin(M_PI/(float)(2*lim)*(float)i+M_PI/2.0),4);}

  
      if (method == 1)
	{
      for (i=0; i<nr; ++i) for (j=0; j<nphi; ++j) avgphi[i]=avgphi[i]+imrphi[j*nr+i]/(double)nphi;

      //get polynomial fit
           

      for (k=0; k<=2; ++k) for (l=0; l<=2; ++l) akl[k][l]=0.0;
      for (k=0; k<=2; ++k) bl[k]=0.0;

      for (i=0; i<nr; ++i) if (avgphi[i] > fmax){fmax=avgphi[i]; max=i;}

   
      
      if ((max-lim) >= 0 && (max+lim) < nr) for (i=(max-lim+1); i<=(max+lim-1); ++i){parab[i-max+lim]=avgphi[i]; ispar=is_parab(parab, parsize, lim);} else ispar=1;
     

	  if (ispar == 0 )
	    {
	      for (k=0; k<=2; ++k) for (l=0; l<=2; ++l) for (i=(max-lim); i<=(max+lim); ++i) akl[k][l]=akl[k][l]+wg[i-(max-lim)]*pow((long double)rc[i], k+l);
	      for (k=0; k<=2; ++k) for (i=(max-lim); i<=(max+lim); ++i) bl[k]=bl[k]+wg[i-(max-lim)]*pow((long double)rc[i], k)*(long double)avgphi[i];

      
      //invert akl
	      long double det=akl[0][0]*(akl[1][1]*akl[2][2]-akl[2][1]*akl[1][2]) - akl[1][0]*(akl[2][2]*akl[0][1]-akl[2][1]*akl[0][2]) + akl[2][0]*(akl[1][2]*akl[0][1]-akl[1][1]*akl[0][2]);
	      long double invmat[3][3];

	      invmat[0][0]=(akl[2][2]*akl[1][1]-akl[2][1]*akl[1][2])/det;
	      invmat[1][0]=-(akl[2][2]*akl[0][1]-akl[2][1]*akl[0][2])/det; invmat[0][1]=invmat[1][0];
	      invmat[2][0]=(akl[1][2]*akl[0][1]-akl[1][1]*akl[0][2])/det;  invmat[0][2]=invmat[2][0];
	      invmat[1][1]=(akl[2][2]*akl[0][0]-akl[2][0]*akl[0][2])/det;
	      invmat[2][1]=-(akl[1][2]*akl[0][0]-akl[1][0]*akl[0][2])/det; invmat[1][2]=invmat[2][1];
	      invmat[2][2]=(akl[1][1]*akl[0][0]-akl[1][0]*akl[0][1])/det;

    
	      for (k=0; k<=2; ++k) cof[k]=0.0;
	      for (k=0; k<=2; ++k) for (l=0; l<=2; ++l) cof[l]=cof[l]+invmat[k][l]*bl[k];
	      double rad_ip=(double)(-cof[1]/cof[2]/2.0);

	      *rsun_lf=(double)rad_ip;
	      status_res=0;
	    }
	  else
	    {
	  *rsun_lf=rad;
	  status_res=1;
	    }
	
	}

	 //
         
      // fit each distance seperately
      if (method == 0)
	{
	  double *rad_ipa=(double *)(malloc(nphi*sizeof(double)));

      for (j=0; j<nphi; ++j){

      for (i=0; i<nr; ++i) avgphi[i]=imrphi[j*nr+i];
      for (k=0; k<=2; ++k) for (l=0; l<=2; ++l) akl[k][l]=0.0;
      for (k=0; k<=2; ++k) bl[k]=0.0;

      fmax=0.0; max=0;
      for (i=0; i<nr; ++i) if (avgphi[i] > fmax){fmax=avgphi[i]; max=i;}
      
      if ((max-lim) >= 0 && (max+lim) < nr){for (i=(max-lim); i<=(max+lim); ++i)parab[i-max+lim]=avgphi[i]; ispar=is_parab(parab, parsize, lim);} else ispar=1;

      if (ispar == 0)
	{
	  for (k=0; k<=2; ++k) for (l=0; l<=2; ++l) for (i=(max-lim); i<=(max+lim); ++i) akl[k][l]=akl[k][l]+wg[i-(max-lim)]*pow((long double)rc[i], k+l);
	  for (k=0; k<=2; ++k) for (i=(max-lim); i<=(max+lim); ++i) bl[k]=bl[k]+wg[i-(max-lim)]*pow((long double)rc[i], k)*(long double)avgphi[i];

      
      //invert akl
	  long double det=akl[0][0]*(akl[1][1]*akl[2][2]-akl[2][1]*akl[1][2]) - akl[1][0]*(akl[2][2]*akl[0][1]-akl[2][1]*akl[0][2]) + akl[2][0]*(akl[1][2]*akl[0][1]-akl[1][1]*akl[0][2]);
	  long double invmat[3][3];
	  invmat[0][0]=(akl[2][2]*akl[1][1]-akl[2][1]*akl[1][2])/det;
	  invmat[1][0]=-(akl[2][2]*akl[0][1]-akl[2][1]*akl[0][2])/det; invmat[0][1]=invmat[1][0];
	  invmat[2][0]=(akl[1][2]*akl[0][1]-akl[1][1]*akl[0][2])/det;  invmat[0][2]=invmat[2][0];
	  invmat[1][1]=(akl[2][2]*akl[0][0]-akl[2][0]*akl[0][2])/det;
	  invmat[2][1]=-(akl[1][2]*akl[0][0]-akl[1][0]*akl[0][2])/det; invmat[1][2]=invmat[2][1];
	  invmat[2][2]=(akl[1][1]*akl[0][0]-akl[1][0]*akl[0][1])/det;
    

	      for (k=0; k<=2; ++k) cof[k]=0.0;
	      for (k=0; k<=2; ++k) for (l=0; l<=2; ++l) cof[l]=cof[l]+invmat[k][l]*bl[k];
	      rad_ipa[j]=(double)(-cof[1]/cof[2]/2.0);
      
	    }
	  else
	    {
	      rad_ipa[j]=NAN;
	    }
      }

	
     double rad_ipv=0.0;
     double rad_count=0.0;
     for (i=0; i<nphi; ++i)
       {
	 if (!isnan(rad_ipa[i])) if (fabs(rad_ipa[i]-rad) < limit_var){rad_ipv+=rad_ipa[i]; rad_count=rad_count+1.0;}
       }

     if (rad_count/(double)nphi > percent_good){*rsun_lf=(double)(rad_ipv/rad_count); status_res=0;} else {status_res=1; *rsun_lf=(double)rad;}

     free(rad_ipa);
 	}
     
       free_mem(&memory);

   return status_res;

}

void free_mem(struct mempointer *memory)
{
  free(memory->imhp);
  free(memory->imro);
  free(memory->mask_p);
  free(memory->image);
 
  free(memory->avgphi);
  free(memory->xrp);
  free(memory->yrp);
  free(memory->imrphi);
  free(memory->rc);
  free(memory->phic);
  free(memory->parab);
  free(memory->cnorm);
  free(memory->ierror);

  return;
}



int is_parab(double *arr, int n, int nmax)
{
  int i;
  int status=0;

  if (nmax < n)
    {
      for (i=(nmax+1); i<n; ++i) 
	{
	  if (arr[i] >= arr[i-1]) status=1;
	}
      for (i=0; i<nmax; ++i)
	{
	  if (arr[i+1] <= arr[i]) status=1;
	}
    }
  else
    {
      status=1;
    }
  return status;
}




void cross_corr(int nx, int ny, double *imhp, double *imro)

  {


    double a[nx][ny];
    double b[nx][ny];


    fftw_complex ac[nx][ny/2+1], bc[nx][ny/2+1];
    fftw_plan p;
    long i,j;

    double scale=1./((double)(nx*ny));


#pragma omp parallel for private(i,j)
      for (j=0; j<ny; ++j)
	for (i=0; i<nx; ++i)
	  {
	    a[i][j]=(double)imhp[j*nx+i];
	    b[i][j]=(double)imro[j*nx+i];
	  }


     //fft(a)
     p = fftw_plan_dft_r2c_2d(nx, ny, &a[0][0], &ac[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);

      
     //fft(b)
     p = fftw_plan_dft_r2c_2d(nx, ny, &b[0][0], &bc[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);  


     //fft(a)*conj(fft(b))

#pragma omp parallel for private(i,j)
       for (i = 0; i < nx; ++i){
	 for (j = 0; j < ny/2+1; ++j){
 	 ac[i][j]=ac[i][j]*conj(bc[i][j])*scale*scale;
       }
     }

     //  fft(fft(a)*conj(fft(b)),1)
     p = fftw_plan_dft_c2r_2d(nx, ny, &ac[0][0], &b[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);

     //  rearrange

#pragma omp parallel for private(i,j)
     for (i = 0; i < ny; ++i){
       for (j = 0; j < nx; ++j){
   	    imro[j*nx+i]=b[(i+nx/2) % nx][(j+ny/2) % ny];
	  }
     }
  

  }







////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
