// Based on hmi/accel/coeffj.pro
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "interpol_code.h"

#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) < (y)) ? (y) : (x))

////constants 


int initialize_interpol(struct initial *const_param, struct init_files *infiles, int nx, int ny, const char *dpath)
{

  char *namep;
  int i, j, s, t;
  float datum;
  int status;
  int nlead=nx;

 int dist_order=6;
 int maxo=dist_order+1;
 float *dist_param_front, *dist_param_side;
 float *distx, *disty;

 int gapfill_order=10;
 int gapfill_method=11;
 float gapfill_regular=0.0295;

 int interpolation_order=10;


 int gapfill_order2=gapfill_order/2;
 int nconst=1; // Number of polynomial terms exactly reproduced in temporal interpolation

 int malign=32;

 struct fill_struct fills;
 status=init_fill(gapfill_method, gapfill_regular, gapfill_order,gapfill_order2,gapfill_order2,&fills,&namep,dpath);
 printf("table for gapfill %s\n", namep);

float **dist_par_front, **dist_par_side;
 dist_par_front=(float **)(MKL_malloc(16*sizeof(float *),malign)); //pointer to distortion coefficients (for 16 focus pos.)
 dist_par_side=(float **)(MKL_malloc(16*sizeof(float *),malign));

 for (j=0; j<16; ++j)
   {
 dist_param_front=(float *)(MKL_malloc((dist_order+1)*(dist_order+1)*2*sizeof(float),malign));
 *(dist_par_front+j)=dist_param_front;

 dist_param_side=(float *)(MKL_malloc((dist_order+1)*(dist_order+1)*2*sizeof(float),malign));
 *(dist_par_side+j)=dist_param_side;
   }

 





 FILE *distfile;
 size_t bytes_read;

    float *dparam_front;
    dparam_front=(float *)(malloc(maxo*maxo*2*sizeof(float)));

    float *dparam_side;
    dparam_side=(float *)(malloc(maxo*maxo*2*sizeof(float)));

    char *filename_distfront=infiles->dist_file_front;

    distfile=fopen(filename_distfront, "rb");

    if (distfile != NULL){
      for (i=0; i<maxo*maxo*2; ++i){
	fscanf(distfile, "%f", &datum);
	dparam_front[i]=datum;
      }
    }
    fclose(distfile);

    char *filename_distside=infiles->dist_file_side;

    distfile=fopen(filename_distside, "rb");
    if (distfile != NULL){
      for (i=0; i<maxo*maxo*2; ++i){
	fscanf(distfile, "%f", &datum);
	dparam_side[i]=datum;
      }
    }
    fclose(distfile);


  


 for (j=0; j<16; ++j)
   {
     dist_param_front=*(dist_par_front+j);
     dist_param_side=*(dist_par_side+j);
     for (i=0; i<maxo*maxo*2; ++i) dist_param_front[i]=dparam_front[i];
     for (i=0; i<maxo*maxo*2; ++i) dist_param_side[i]=dparam_side[i];
   }





 float *diffrot_param;

 int diffrot_order;
 char *rotcoef_file;
 rotcoef_file=infiles->diffrot_coef;

    FILE *diffrot;
    diffrot=fopen(rotcoef_file, "r");  

    if (diffrot != NULL){
	fscanf(diffrot, "%d", &diffrot_order);
	if (diffrot_order > 0)
	  {
	    diffrot_param=(float *)(MKL_malloc((diffrot_order+1)*sizeof(float),malign));
	    for (j=0; j<=diffrot_order; ++j)
	      {
		fscanf(diffrot, "%f", &datum);
		diffrot_param[j]=datum;
	      }
	  }
    }

   

   

 for (j=0; j<16; ++j){dist_param_front=*(dist_par_front+j); const_param->dist_coef[0][j]=dist_param_front;}        //distortion parameters  //0 : front camera
 for (j=0; j<16; ++j){dist_param_side=*(dist_par_side+j); const_param->dist_coef[1][j]=dist_param_side;}


 const_param->order_dist_coef=dist_order;    // order of distortion palynominal
 const_param->diffrot_coef=diffrot_param;    // differentail rotation coefficients
 const_param->order2_rot_coef=diffrot_order; // order of differental rotation polynominal in sin(lat)^2
 const_param->order_int=interpolation_order;                 // order for spatial  interpolation
 const_param->nconst=nconst;
 const_param->fills=fills;

 const_param->code_version="INTERPOL_2.0";
 MKL_free(dist_par_front);
 MKL_free(dist_par_side);

 return status;
}





void free_interpol(struct initial *const_param)
{
  int j;
  free_fill(&const_param->fills);
  for (j=0; j<16; ++j) MKL_free(const_param->dist_coef[0][j]);
  for (j=0; j<16; ++j) MKL_free(const_param->dist_coef[1][j]);

  //  MKL_free(const_param->dist_vec[0]);
  //  MKL_free(const_param->dist_vec[1]);
  //  MKL_free(const_param->dist_vec[2]);
  //  MKL_free(const_param->dist_vec[3]);

  //  MKL_free(const_param->dist_vec);

  MKL_free(const_param->diffrot_coef);
}




int do_gapfill(float *image, unsigned char *mask, struct initial *const_param, char *ierror, int nx, int ny)
{

  int isample;
  int status;
  float t1;
  float *cnorm, *ierror_float;
  int i, j;
  int nlead=nx;
  int malign=32;

  cnorm=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
  ierror_float=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));

  struct fill_struct fills=const_param->fills;

  status=fgap_fill(&fills,image,nx,ny,nlead,mask,cnorm,ierror_float);

  for (j=0;j<ny;j++) for (i=0;i<nx;i++) ierror[j*nlead+i]=(char)(ierror_float[j*nlead+i]*64.0-64.0);
   
  MKL_free(cnorm);
  MKL_free(ierror_float);


  return status;


}


int do_interpolate(float **images, char **ierrors, float *image_out, struct keyword *key, struct keyword *key_out, struct initial *const_param, int nsample, int nx, int ny, float average, const char *dpath)
{

  void derotation(float, float, float, float, float, float, float *, int, float *, int, int);
  void derotation_full(struct fint_struct *fints, float *, float *, double, float, float, float, float, float, float, float, float, float, float, float, float, float *, int, int, int, int, float *, float *);
  long factorial(int);
  float intsincos(unsigned int n, unsigned int m);

 
    

  int i,j,isample,s,t,k;  // loop variables
  float **int_images;   // pointer to interpolated images
  unsigned char **masks;
  float *image_intout;  // pointer to single interpolated image
  double t1, t0;    // time variable (beta)
  int *status_arr;  //status of spatial interpolation
  int status;
  struct keyword keyl;  // input keyword for single image

  int nlead=nx;
  int malign=32;
  float maxext=0.0;
  float fillval=sqrt(-1.0); // NaN
  float thresold=1.0;

  //initialize
 
  int_images=(float **)(MKL_malloc(nsample*sizeof(float *), malign));
  masks=(unsigned char **)(MKL_malloc(nsample*sizeof(unsigned char *),malign));

  char *namep;
  float *par;

  float *rrsun, *bb0, *pp0, *xx0, *yy0, *ddist;
  int *ccamera, *ffocus;
  double *tsample;

  ccamera=(int *)(MKL_malloc(nsample*sizeof(int),malign));
  ffocus=(int *)(MKL_malloc(nsample*sizeof(int),malign));
  rrsun=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  bb0=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  pp0=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  xx0=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  yy0=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  ddist=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  tsample=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  
  status_arr=(int *)(MKL_malloc(nsample*sizeof(int),malign));

  float rsun, b0, p0, cent_x, cent_y, dist;
  double tint; //target time
  
 


 float *distx, *disty;
 static float *dist_vec[4];
  float *dcoef_front, *dcoef_side;
  //

  //read out keywords
  for (isample=0; isample<nsample; ++isample){
    keyl=*(key+isample);
    ccamera[isample]=keyl.camera;
    ffocus[isample]=keyl.focus;
    rrsun[isample]=keyl.rsun;
    bb0[isample]=keyl.b0;
    pp0[isample]=keyl.p0;
    xx0[isample]=keyl.xx0;
    yy0[isample]=keyl.yy0;
    ddist[isample]=keyl.dist;
    tsample[isample]=keyl.time;
     }

  
  //distx_front=const_param->dist_vec[0];
  //disty_front=const_param->dist_vec[1];
  //distx_side=const_param->dist_vec[2];
  //disty_side=const_param->dist_vec[3];



  int maxo=const_param->order_dist_coef+1;
      
  float coef_x_front[maxo][maxo], coef_y_front[maxo][maxo], coef_x_side[maxo][maxo], coef_y_side[maxo][maxo];
  float y1[maxo], x1[maxo];

  float rad=(float)nx/2.0;


  static int focus_current=0;

  float *distx_front, *disty_front, *distx_side, *disty_side;

  if (focus_current == 0)
    {
      distx_front=(float *)(malloc(nx*ny*sizeof(float)));
      disty_front=(float *)(malloc(nx*ny*sizeof(float)));
      distx_side=(float *)(malloc(nx*ny*sizeof(float)));
      disty_side=(float *)(malloc(nx*ny*sizeof(float)));
    }


  if (ffocus[0] != focus_current)
    {
      focus_current=ffocus[0];

         dcoef_front=const_param->dist_coef[0][ffocus[0]-1];
	 dcoef_side=const_param->dist_coef[1][ffocus[0]-1];


      t0=dsecnd();
 //define distortion in CCD coordinates
  
      for (s=0; s<maxo; ++s){for (t=0; t<maxo; ++t){coef_x_front[s][t]=dcoef_front[s*maxo+t]; coef_y_front[s][t]=dcoef_front[s*maxo+t+maxo*maxo];}}
      for (s=0; s<maxo; ++s){for (t=0; t<maxo; ++t){coef_x_side[s][t]=dcoef_side[s*maxo+t]; coef_y_side[s][t]=dcoef_side[s*maxo+t+maxo*maxo];}}


    

#pragma omp parallel for private(j,i,y1,x1,s,t)
      for (j=0; j<ny; ++j)
	for (i=0; i<nx; ++i)
          {
	    y1[0]=1.0;
            y1[1]=((float)j-rad)/rad; // center at CCD center, normalize to "CCD radius"
	    for (t=2; t<maxo; ++t)y1[t]=y1[t-1]*y1[1];

	    x1[0]=1.0;
	    x1[1]=((float)i-rad)/rad;
	    for (s=2; s<maxo; ++s)x1[s]=x1[s-1]*x1[1];

	   

	    distx_front[j*nx+i]=0.0;
	    disty_front[j*nx+i]=0.0;

	    distx_side[j*nx+i]=0.0;
	    disty_side[j*nx+i]=0.0;

	    for (s=0; s<maxo; ++s){
	      for (t=0; t<maxo; ++t){
		if (coef_x_front[s][t] != 0.0) distx_front[j*nx+i]=distx_front[j*nx+i]+coef_x_front[s][t]*x1[s]*y1[t]*rad;
		if (coef_y_front[s][t] != 0.0) disty_front[j*nx+i]=disty_front[j*nx+i]+coef_y_front[s][t]*x1[s]*y1[t]*rad;

		if (coef_x_side[s][t] != 0.0) distx_side[j*nx+i]=distx_side[j*nx+i]+coef_x_side[s][t]*x1[s]*y1[t]*rad;
		if (coef_y_side[s][t] != 0.0) disty_side[j*nx+i]=disty_side[j*nx+i]+coef_y_side[s][t]*x1[s]*y1[t]*rad;
	      }
	    }

	  }

      dist_vec[0]=distx_front;
      dist_vec[1]=disty_front;      	
      dist_vec[2]=distx_side;
      dist_vec[3]=disty_side;



    }


  



 
  rsun=key_out->rsun;
  b0=key_out->b0;
  p0=key_out->p0;
  cent_x=key_out->xx0;
  cent_y=key_out->yy0;
  dist=key_out->dist;
  tint=key_out->time;
  //////////////////////////

  float p0min=pp0[0],p0max=pp0[0],b0min=bb0[0],b0max=bb0[0],rsmin=rrsun[0],rsmax=rrsun[0],x0min=xx0[0],x0max=xx0[0],y0min=yy0[0],y0max=yy0[0];

  float limp0=0.1;
  float limb0=0.05;
  float limrsun=0.1;
  float limx0=1.0;
  float limy0=1.0;



  // obtain range of keywords
  for (isample=1; isample < nsample; ++isample){
    
   p0min=minval(pp0[isample], p0min);
   p0max=maxval(pp0[isample], p0max);
   b0min=minval(bb0[isample], b0min);
   b0max=maxval(bb0[isample], b0max);
   rsmin=minval(rrsun[isample], rsmin);
   rsmax=maxval(rrsun[isample], rsmax);
   x0min=minval(xx0[isample], x0min);
   x0max=maxval(xx0[isample], x0max);
   y0min=minval(yy0[isample], y0min);
   y0max=maxval(yy0[isample], y0max);
    }

  ////////////////////////


  //check for valid range of keywords

  // for (isample=0; isample < nsample; ++isample){
  //
  //  if (pp0[isample] < -2.0*M_PI || pp0[isample] > 2.0*M_PI) status_arr[isample]=1; else status_arr[isample]=0;
  //  if (bb0[isample] < -7.5/180.0*M_PI || bb0[isample] > 7.5/180.0*M_PI) status_arr[isample]=1; else status_arr[isample]=0;
  //  if (rrsun[isample] < 1600.0 || rrsun[isample] > 2048.0)  status_arr[isample]=1; else status_arr[isample]=0;
  //  if (xx0[isample] < 0 || xx0[isample] > 4095) status_arr[isample]=1; else status_arr[isample]=0;
  //  if (yy0[isample] < 0 || yy0[isample] > 4095) status_arr[isample]=1; else status_arr[isample]=0;
  //  if (ddist[isample] < 0.9 || ddist[isample] > 1.1) status_arr[isample]=1; else status_arr[isample]=0;
   
  //}

  int derot_status=0;

  if ((p0max-p0min) > limp0) derot_status=1;
  if ((b0max-b0min) > limb0) derot_status=1;
  if ((rsmax-rsmin) > limrsun) derot_status=1;
  if ((x0max-x0min) > limx0) derot_status=1;
  if ((y0max-y0min) > limy0) derot_status=1;

  printf("derot status %d\n", derot_status);
  ///////////

  for (isample=0; isample < nsample; ++isample) status_arr[isample]=0;
 

 

  t0=dsecnd();
  t1=dsecnd();

  float *xin, *yin;     //input coordinated

  xin=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
  yin=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));


  struct fint_struct fints, fints_error;
  
  
  int order_int=const_param->order_int;   //get interpolation order
  int order2_int=order_int/2;             // interpolation order /2

 
 
  
  init_finterpolate_wiener(&fints,order_int,0,1.0,2,1,1,&namep, dpath);  //initialize interpolation
  printf("table name for interpolation %s\n", namep);
  free(namep);

  init_finterpolate_linear(&fints_error, 0.0); // linear interpolation for masks

  
  double time;  // image time - target time
 
  


  float coef_x[maxo][maxo], coef_y[maxo][maxo];  //distortion coefficients in CCD coordinates
 
  float *dcoef;
 
  //calculate rotated coefficients
  ///////////////////////////////////

  // float coefx_rot[nsample][maxo][maxo], coefy_rot[nsample][maxo][maxo];
  long bin_si, bin_tj, bin_sip, bin_tjp;
  float av_offset_x=0.0;
  float av_offset_y=0.0;
  int sp, tp, ip, jp;
  float rad_term1=0.0, rad_term2=0.0, rad_term3=0.0, rad_term4=0.0;

 

  //initailize distortion coefficients in target image coordinates

  //for (isample=0; isample<nsample; ++isample) for (s=0; s<maxo; ++s) for (t=0; t<maxo; ++t){coefx_rot[isample][s][t]=0.0; coefy_rot[isample][s][t]=0.0;} //initialize
 

  for (isample=0; isample<nsample; ++isample){

      dcoef=const_param->dist_coef[ccamera[isample]][ffocus[isample]-1];
      for (s=0; s<maxo; ++s){for (t=0; t<maxo; ++t){coef_x[s][t]=dcoef[s*maxo+t]; coef_y[s][t]=dcoef[s*maxo+t+maxo*maxo];}}
  


    float xc=(xx0[isample]-rad)/rad;   //target image center (centered and normalized)
    float yc=(yy0[isample]-rad)/rad;


  //calculate rsun-change due to removing distortion, and image center change due to removing distortion:
  //only output keywords are affected: Output rsun and image center should be center of undistorted image
    //disabled for Zernike polynomial version: Too complicated to calculate
    av_offset_x=0.0;
    av_offset_y=0.0;
    rad_term1=0.0; rad_term2=0.0; rad_term3=0.0; rad_term4=0.0; 

       for (s=0; s<maxo; ++s) for (t=0; t<=(maxo-1-s); ++t) for (i=0; i<=s; ++i) for (j=0; j<=t; ++j){

	 bin_si=(float)factorial(s)/factorial(i)/factorial(s-i);
	 bin_tj=(float)factorial(t)/factorial(j)/factorial(t-j);

       av_offset_x += coef_x[s][t]*bin_si*bin_tj*pow(rrsun[isample]/rad, s+t-(i+j))*pow(xc, i)*pow(yc, j)*intsincos(s-i+2,t-j)*2.0*rad; //calculate offs offset due to distortion
       av_offset_y += coef_y[s][t]*bin_si*bin_tj*pow(rrsun[isample]/rad, s+t-(i+j))*pow(xc, i)*pow(yc, j)*intsincos(s-i,t-j+2)*2.0*rad;
    
       rad_term1 += 2.0*coef_x[s][t]*bin_si*bin_tj*pow(rrsun[isample]/rad,s+t-(i+j)+1)*pow(xc,i)*pow(yc,j)*intsincos(s-i+1,t-j);
       rad_term3 += 2.0*coef_y[s][t]*bin_si*bin_tj*pow(rrsun[isample]/rad,s+t-(i+j)+1)*pow(xc,i)*pow(yc,j)*intsincos(s-i,t-j+1);

       for (sp=0; sp<maxo; ++sp) for (tp=0; tp<=(maxo-1-s); ++tp) for (ip=0; ip<=sp; ++ip) for (jp=0; jp<=tp; ++jp){ //contour integral over (circular) limb 

	bin_sip=(float)factorial(sp)/factorial(ip)/factorial(sp-ip);
    	bin_tjp=(float)factorial(tp)/factorial(jp)/factorial(tp-jp);

    	rad_term2 += coef_x[s][t]*coef_x[sp][tp]*bin_si*bin_tj*bin_sip*bin_tjp*pow(rsun/rad,s+t+sp+tp-(i+j+ip+jp))*pow(xc, i+ip)*pow(yc, j+jp)*intsincos(s+sp-i-ip, t+tp-j-jp);
    	rad_term4 += coef_y[s][t]*coef_y[sp][tp]*bin_si*bin_tj*bin_sip*bin_tjp*pow(rsun/rad,s+t+sp+tp-(i+j+ip+jp))*pow(xc, i+ip)*pow(yc, j+jp)*intsincos(s+sp-i-ip, t+tp-j-jp);
      }
      }

   
   
   
     //apply correction due to distortion; 
       xx0[isample]=xx0[isample]-av_offset_x;  // has to be removed when nominal values are used
       yy0[isample]=yy0[isample]-av_offset_y;
       rrsun[isample]=sqrt(rrsun[isample]*rrsun[isample]-(rad_term1+rad_term2+rad_term3+rad_term4)*rad*rad); 

     
  }

     ////////////////

  t1=dsecnd()-t1;



  t1=dsecnd();

  int order2_rot_coef=const_param->order2_rot_coef;   //squared order of differential rotation coefficients

  float *rot_coef;
  rot_coef=(float *)(MKL_malloc((order2_rot_coef+1)*sizeof(float), malign)); 

  for (i=0; i<=order2_rot_coef; ++i) rot_coef[i]=const_param->diffrot_coef[i];
      


  t1=dsecnd()-t1;

  
 /////preparation for mask building
 float *ierror_out, *ierror_float;   
 unsigned char *maskp;
 ierror_float=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
 ierror_out=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
 for (j=0; j<nlead; ++j) for (i=0; i<nx; ++i) ierror_out[j*nx+i]=0.0;

 float *shift;
 shift=(float *)(malloc(nlead*ny*2*sizeof(float)));


 derotation(rsun, cent_x, cent_y, dist, p0, b0, rot_coef, order2_rot_coef, shift, nx, ny);
 ///////////////////////////////////////////////////////////////
  for (isample=0; isample<nsample; ++isample){ //loop over images
    printf("pang %f %f\n", pp0[isample], p0); 
    if (status_arr[isample] == 0){
     

    time=(double)(tsample[isample]-tint); // difference to target time
    printf("time lag: %g \n", time);

    //distx=const_param->dist_vec[ccamera[isample]*2];
    //disty=const_param->dist_vec[ccamera[isample]*2+1];

    distx=dist_vec[ccamera[isample]*2];
    disty=dist_vec[ccamera[isample]*2+1];


    
    ////calculate derotation pattern
    t1=dsecnd();
    
    if (derot_status == 1)
    derotation_full(&fints_error, distx, disty, time, rrsun[isample], rsun, xx0[isample], cent_x, yy0[isample], cent_y, ddist[isample], dist, pp0[isample], p0, bb0[isample], b0, rot_coef, order2_rot_coef,const_param->order_dist_coef, nx, ny, xin, yin);
    
       
     if (derot_status == 0)
       {
#pragma omp parallel for private(j,i)
    for (j=0; j<ny; ++j)
      for (i=0; i<nx; ++i)
	{
	  xin[j*nlead+i]=(float)i + distx[j*nlead+i]+(xx0[isample]-cent_x)+((float)i-cent_x)/rsun*(rrsun[isample]-rsun)+shift[j*nlead+i]*time;
	  yin[j*nlead+i]=(float)j + disty[j*nlead+i]+(yy0[isample]-cent_y)+((float)j-cent_y)/rsun*(rrsun[isample]-rsun)+shift[ny*nlead+j*nlead+i]*time;
	}
    }


    //////////////spatial interpoaltion
  
    image_intout=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));

    
    status_arr[isample]=finterpolate(&fints,*(images+isample), xin, yin, image_intout, nx,ny,nlead,nx,ny,nlead,fillval); //spatial interpolation
    t1=dsecnd()-t1;  
    
 
   
   *(int_images+isample)=image_intout;  //pointer to each output image
   /////////////////////////

   
   

   ///////make masks out of errors
  
   for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) ierror_float[j*nlead+i]=(float)((*(ierrors+isample))[j*nlead+i]+64)/64.0;
   finterpolate(&fints_error,ierror_float, xin, yin, ierror_out, nx,ny,nlead,nx,ny,nlead,2.0f);

  
   maskp=(unsigned char *)(MKL_malloc(nlead*ny*sizeof(unsigned char),malign));
   for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) if (ierror_out[j*nlead+i] > thresold || isnan(image_intout[j*nlead+i])) maskp[j*nlead+i]=1; else maskp[j*nlead+i]=0;
   *(masks+isample)=maskp;
 
 
    }
  
  } // end loop over images
    




  //end interpolation///////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //count valid images
  int nvalid=0;
 
  for (isample=0; isample<nsample; ++isample){nvalid=nvalid+status_arr[isample];}
  nvalid=nsample-nvalid;



 //get valid subsample

  double *tvalid;
  float **valid_images;
  unsigned char **valid_masks;

 
  tvalid=(double *)(MKL_malloc(nvalid*sizeof(double),malign));
  valid_images=(float **)(MKL_malloc(nvalid*sizeof(float *),malign));
  valid_masks=(unsigned char **)(MKL_malloc(nvalid*sizeof(unsigned char *),malign));

  int count=0;
  for (isample=0; isample<nsample; ++isample){
    if (status_arr[isample] == 0){
      tvalid[count]=tsample[isample];
      valid_images[count]=int_images[isample];
      valid_masks[count]=masks[isample];
      ++count;
    }
  }
  
   

    //temporal interpolatiom 
  //printf("nvalid: %u %u \n", nvalid, count);

  char *namep_temp;
  t1=dsecnd();
  if (average < 0.0)
    {
      status=tinterpolate(nvalid,tvalid,tint,const_param->nconst,valid_images,valid_masks,image_out,nx,ny,nlead,1,&namep_temp, fillval,dpath);  //call temporal interpolation function
    }
  else {
  
    int torder=4;
    double avgval=round(average/2.0/45.0/2.0)*2.0;
    int hwidth=(int)(avgval*1.5);
    status=taverage(nvalid,tvalid,tint,const_param->nconst,valid_images,valid_masks,image_out,nx,ny,nlead,1,&namep_temp, tavg_cosine,torder,45.0,hwidth,avgval,avgval/2.0,fillval,dpath);
  }

 
  printf("table name for temporal interpolation %s\n", namep_temp);
  free(namep_temp);

  t0=dsecnd()-t0;


  MKL_FreeBuffers();
 

 for (isample=0;isample<nsample;isample++){if (status_arr[isample] ==0){MKL_free(*(int_images+isample)); MKL_free(*(masks+isample));}}
 free_finterpolate(&fints_error);
 free_finterpolate(&fints);
 free(shift);

 MKL_free(int_images);
 
 MKL_free(masks);
 
 MKL_free(xin);
 
 MKL_free(yin);
 


 MKL_free(rot_coef);
 


 MKL_free(ccamera);
 MKL_free(ffocus);
 MKL_free(tsample);
 MKL_free(ddist); 
 MKL_free(xx0);
 MKL_free(yy0);
 MKL_free(pp0);
 MKL_free(bb0);
 MKL_free(rrsun);
 
 MKL_free(status_arr);
 
 MKL_free(ierror_out);
 MKL_free(ierror_float);
 
 if (nvalid > 0){
   MKL_free(tvalid);
   MKL_free(valid_images);
   MKL_free(valid_masks);
 }
 ///////////////////////


 return status;  //status of temporal interpolation
 
 

}


void derotation_full(struct fint_struct *fints, float *distx, float *disty, double time, float radius, float radius_t, float cent_x, float cent_x_t, float cent_y, float cent_y_t, float dist, float dist_t, float p0, float p0_t, float b0, float b0_t, float *rot_coef, int order2_rot_coef, int order_dist, int nx, int ny, float *xin, float *yin){
    

  int i,j, k, s, t; //loop variables


    
    float xyzr[3], xyzb[3];
   
    int nlead=nx;

    float xy[2], xy_t[2], xyp_t[2], xyz[3], xyz_t[3];
    float inr, ind;
    float inr_t, ind_t;
    float  xyp[2], xyo[2];
    float ikf;

    float phie, beta;
    float slat;
    float omeg_t;
    float singam, cosgam, phi, ing, inf;

    const float au=215.020*dist; // (1 astonomical unit in solar radii)
    const float au_t=215.020*dist_t;

    int maxo=order_dist+1;
 
    float maxphi=asin(1.0/au);
    float maxphi_t=asin(1.0/au_t);
    float coef_x[maxo][maxo], coef_y[maxo][maxo];

    float *xinput, *yinput;
    xinput=(float *)(malloc(nlead*ny*sizeof(float)));
    yinput=(float *)(malloc(nlead*ny*sizeof(float)));
 
  
      ///////////////////////////
   
  

  

    //  for (s=0; s<maxo; ++s){for (t=0; t<maxo; ++t){coef_x[s][t]=dparam[s*maxo+t]; coef_y[s][t]=dparam[s*maxo+t+maxo*maxo];}}

  

 
  float rad=(float)nx/2.0;
    float x1a,x2a,x3a,xx,yy, y1,y2,y3;

    float cos_p0_t=cos(p0_t);
    float sin_p0_t=sin(p0_t);
    float cos_b0_t=cos(b0_t);
    float sin_b0_t=sin(b0_t);
    float cos_b0=cos(b0);
    float sin_b0=sin(b0);
    float sin_p0=sin(p0);
    float cos_p0=cos(p0);


  
 ////////////////////////////////////////////////





    

#pragma omp parallel for private(j,i,k,xy,xy_t,xyp,inf,inr_t,xyp_t,phie,beta,ind_t,xyz_t,xyzb,slat,omeg_t,xyzr,xyz,singam,cosgam,phi,ing,ikf,xyo,x1a,x2a,x3a,y1,y2,y3,xx,yy)
   for (j=0; j<ny; j++){
       for (i=0; i< nx; i++){

       
	 
	   xy[1]=(float)j;
	   xy_t[1]=((float)j-cent_y_t)/radius_t; // normalized xy coordinate of target image in the sky (solar radius=1)
             
	   xy[0]=(float)i;
	   xy_t[0]=((float)i-cent_x_t)/radius_t;   // normalized xy coordinate of target image in the sky (solar radius=1)
	     
	   inr_t=sqrt(xy_t[0]*xy_t[0]+xy_t[1]*xy_t[1]); // radial coordinate in the sky / target
 
	     xyp_t[0]=cos_p0_t*xy_t[0]+sin_p0_t*xy_t[1];
	     xyp_t[1]=-sin_p0_t*xy_t[0]+cos_p0_t*xy_t[1];  // p_angle rotation p-angle=p0_t -> p_angle=0;


	  if (inr_t >= 1.0)
	    {
	      xy[0]=(xy[0]-cent_x_t)/inr_t+cent_x_t;  // force off-limb point to radial distance 1.0
	      xy[1]=(xy[1]-cent_y_t)/inr_t+cent_y_t;
	      xyp_t[0]=xyp_t[0]/inr_t;
	      xyp_t[1]=xyp_t[1]/inr_t;  // force off-limb point to radial distance 1.0
	      inr_t=1.0;
	      phie=maxphi_t;
	      ind_t=sqrt(1.0-1.0/au_t/au_t);
            }
	  else
            {
	      phie=inr_t*maxphi_t;
	      beta=tan(phie);
	      ind_t=(beta*au_t-beta*sqrt(1.0+beta*beta-au_t*au_t*beta*beta))/(1+beta*beta);
	    }


 
 
	  if (inr_t != 0.0) ikf=ind_t/inr_t; else ikf=(au_t-1.0)/au_t; //limit for exact center //!! infinite distance approximation
       	
	  xyz_t[0]=sqrt(1.0-ind_t*ind_t);
	  xyz_t[1]= xyp_t[0]*ikf;
	  xyz_t[2]= xyp_t[1]*ikf;  //project from target coordinates onto sphere


	  xyzb[0]=cos_b0_t*xyz_t[0]-sin_b0_t*xyz_t[2]; // rotate from b_angle=b0_t -> b_angle=0; 
          xyzb[1]=xyz_t[1];
	  xyzb[2]=sin_b0_t*xyz_t[0]+cos_b0_t*xyz_t[2];

	 

	  slat=xyzb[2];  //sine of heliographic latitude
	  
	  
	  omeg_t=1e-6*(rot_coef[1]*pow(slat,2)+rot_coef[2]*pow(slat,4)+rot_coef[0])*(float)time; 
	  
	  

	  xyzr[0]=cos(omeg_t)*xyzb[0]-sin(omeg_t)*xyzb[1];  //rotate sphere around z-axis
          xyzr[1]=sin(omeg_t)*xyzb[0]+cos(omeg_t)*xyzb[1];
	  xyzr[2]=xyzb[2];    


	  // svec[0]=-omeg*xyz[1]*cos_b0*time;
	  // svec[1]=omeg*(xyz[0]*cos_b0-xyz[2]*sin_b0)*time; // linearized rotation
	  // svec[2]=omeg*xyz[1]*sin_b0*time; //


          xyz[0]=cos_b0*xyzr[0]+sin_b0*xyzr[2]; // rotate to b-angle=b0
          xyz[1]=xyzr[1];
	  xyz[2]=-sin_b0*xyzr[0]+cos_b0*xyzr[2];


	  // if (xyz[0] >= 1.0/au)
	  //   {
	      singam=sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]);   //inverse projection into original image coordinates / calculate scale factor
	      cosgam=xyz[0];
	      phi=atan(singam/(au-cosgam));
	      ing=phi/maxphi;
	      //  }
	      // else 
	      //   {
	      //    xyz[0]=1.0/au;
	      //   inf=sqrt(1.0-1.0/au/au)/sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]);
	      //    xyz[1]=xyz[1]*inf;
	      //   xyz[2]=xyz[2]*inf;

	      //   singam=sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]);   //inverse projection into original image coordinates / calculate scale factor
	      //   cosgam=xyz[0];
	      //   phi=atan(singam/(au-cosgam));
	      //   ing=phi/maxphi;


	      //  }


	  if (singam != 0.0) ikf=ing/singam; else ikf=au_t/(au_t-1.0); //limit for exact center

	 
	  xyp[0]=xyz[1]*ikf;     // do inverse projection
	  xyp[1]=xyz[2]*ikf;
	 


	 
          xyo[0]=(cos_p0*xyp[0]-sin_p0*xyp[1]); // inverse p-angle rotation p_angle=0 -> p_angle=p0
	  xyo[1]=(sin_p0*xyp[0]+cos_p0*xyp[1]); //

	  //xin[j*nx+i]=xyo[0]*radius+cent_x-xy[0]; //scale to object / target radius / original - target
	  //yin[j*nx+i]=xyo[1]*radius+cent_y-xy[1];
	   
   


    
   
	  xx=(float)i+xyo[0]*radius+cent_x-xy[0];  //object coordinates
	  yy=(float)j+xyo[1]*radius+cent_y-xy[1];    //object coordinates

	  //  y1=(yy-rad)/rad; // center at CCD center, normalize to "CCD radius"
	  //  y2=y1*y1;
	  //  y3=y2*y1;

	  // x1a=(xx-rad)/rad;  // center at CCD center, normalize to "CCD radius"
	  // x2a=x1a*x1a;
	  // x3a=x2a*x1a;
	
  
   
	 // xin[j*nlead+i]=xx+(coef_x[2][0]*x2a+coef_x[3][0]*x3a+coef_x[1][1]*x1a*y1+coef_x[2][1]*x2a*y1+coef_x[0][2]*y2+coef_x[1][2]*x1a*y2+coef_x[0][3]*y3)*rad;
	 // yin[j*nlead+i]=yy+(coef_y[2][0]*x2a+coef_y[3][0]*x3a+coef_y[1][1]*x1a*y1+coef_y[2][1]*x2a*y1+coef_y[0][2]*y2+coef_y[1][2]*x1a*y2+coef_y[0][3]*y3)*rad;

	 xinput[j*nlead+i]=xx;
	 yinput[j*nlead+i]=yy;
     }
   }


   float *distx_out, *disty_out;
   distx_out=(float *)(malloc(nlead*ny*sizeof(float)));
   disty_out=(float *)(malloc(nlead*ny*sizeof(float)));

   finterpolate(fints,distx, xinput, yinput, distx_out, nx,ny,nlead,nx,ny,nlead,0.0);
   finterpolate(fints,disty, xinput, yinput, disty_out, nx,ny,nlead,nx,ny,nlead,0.0);
 

  
#pragma omp parallel for private(j,i)
    for (j=0; j<ny; ++j)
     for (i=0; i<nx; ++i)
       {
        xin[j*nlead+i]=xinput[j*nlead+i]+distx_out[j*nlead+i];
        yin[j*nlead+i]=yinput[j*nlead+i]+disty_out[j*nlead+i];
       }

    free(distx_out);
    free(disty_out);
    free(xinput);
    free(yinput);

}


void derotation(float radius, float cent_x, float cent_y, float dist, float p0, float b0, float *rot_coef, int order2_rot_coef, float *shift, int nx, int ny){
    

  int i, j, k;
  float xy[2], xyp[2];
  float xyz[3], sxyz[3];
  float sxyp[2], sxy[2];
  float inr, ind, ikf, phie, ing, beta;

   
  float slat, omeg;
  float singam, cosgam, phi; 
      
    
   const float au=215.020*dist; // (1 astonomical unit in solar radii) au=149597870.691 km, R_SUN=695740 km (Kuhn 2004)
   const float maxphi=asin(1.0/au);
 
   
   float sinp0=sin(p0);
   float cosp0=cos(p0);
   float sinb0=sin(b0);
   float cosb0=cos(b0);

   int nlead=nx;


#pragma omp parallel for private(i,j,xyp,inr,phie,ind,ikf,beta,xyz,slat,omeg,sxyz,sxyp,sxy)
   for (j=0; j<ny; ++j){
      for (i=0; i<nx; ++i){
   
	xy[0]=(cosp0*(i-cent_x)-sinp0*(j-cent_y))/radius;     
	xy[1]=(sinp0*(i-cent_x)+cosp0*(j-cent_y))/radius;  // normalized xy coordinate in the sky (solar radius=1)	



  // p-angle rotation

	  inr=sqrt(xy[0]*xy[0]+xy[1]*xy[1]); // radial coordinate in the sky
         
 
             
	  if (inr >= 1.0)
	    {
	      xy[0]=xy[0]/inr;
	      xy[1]=xy[1]/inr;  // force off-limb point to radial distance 1.0
	      inr=1.0;

	      phie=maxphi;
	      ind=sqrt(1.0-1.0/au/au);
            }
	  else
            {
	      phie=inr*maxphi;
	      beta=tan(phie);
              ind=(beta*au-beta*sqrt(1.0+beta*beta-au*au*beta*beta))/(1+beta*beta);
            }

	  if (inr != 0.0) ikf=ind/inr; else ikf=(au-1.0)/au; //limit for exact center


       	  xyz[0]=sqrt(1.0-ind*ind);
	  xyz[1]=xy[0]*ikf;
	  xyz[2]=xy[1]*ikf;

	  	      
	  slat=sinb0*xyz[0]+cosb0*xyz[2];  //sine of heliographic latitude (finite distance approximation)

		
	  omeg=1e-6*(rot_coef[1]*pow(slat,2)+rot_coef[2]*pow(slat,4)+rot_coef[0]);

	  sxyz[0]=-omeg*xyz[1]*cosb0;
	  sxyz[1]=omeg*(xyz[0]*cosb0-xyz[2]*sinb0); // linearized rotation
	  sxyz[2]=omeg*xyz[1]*sinb0; //

	  singam=sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]); // projection from sphere onto observing plane
	  cosgam=xyz[0];
	  phi=atan(singam/(au-cosgam));
	  ing=phi/maxphi;

	  
	  if (singam != 0.0) ikf=ing/singam; else ikf=au/(au-1.0); //limit for exact center

          sxyp[0]=sxyz[1]/ikf;
	  sxyp[1]=sxyz[2]/ikf;

	  sxy[0]=cosp0*sxyp[0]-sinp0*sxyp[1];
	  sxy[1]=sinp0*sxyp[0]+cosp0*sxyp[1];  // inverse p-angle rotation 

	  


	  shift[j*nlead+i]=sxy[0]*radius; // !! changed to sxyp
	  shift[ny*nlead+j*nlead+i]=sxy[1]*radius;
	   
     }
   }

}











long factorial(int n)
{
int i;
long res=1;
if (n >=1){
  for (i=1; i<=n; ++i) res=res*(long)i;
 }
return res;
}


float intsincos(unsigned int n, unsigned int m)
{
  int mp, np;
  float f=1.0;

  if (n % 2 == 0 && m % 2 ==0){
       
      for (np=2; np<=n; np=np+2){
	f=f*((float)np-1.0)/((float)np+(float)m);
      }
   
   
      for (mp=2; mp<=m; mp=mp+2){
	f=f*((float)mp-1.0)/(float)mp;
      }
   
    return f;
  }
  else return 0.0;

}

char *interpol_version() // Returns CVS version
{
  return strdup("$Id: interpol_code.c,v 1.2 2011/12/06 18:11:03 arta Exp $");
}


//icc -O3 -o interpol_main interpol_main.c interpol_code.c tinterpolate.c finterpolate.c cholesky_down.c gapfill.c -lmkl_em64t -lguide -lpthread -openmp
