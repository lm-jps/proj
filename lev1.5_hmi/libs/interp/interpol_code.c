// Based on hmi/accel/coeffj.pro
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "interpol_code.h"



////constants 


int initialize_interpol(struct initial *const_param)
{

  int i;
  float datum;
  int status;

 int dist_order=3;
 float *dist_param;

 int gapfill_order=5;
 int gapfill_method=11;
 float gapfill_regular=0.0025;

 int interpolation_order=10;


 int gapfill_order2=gapfill_order/2;
 int nconst=1; // Number of polynomial terms exactly reproduced in temporal interpolation

 int malign=32;

 struct fill_struct fills;
 status=init_fill(gapfill_method, gapfill_regular, gapfill_order,gapfill_order2,gapfill_order2,&fills);


 dist_param=(float *)(MKL_malloc((dist_order+1)*(dist_order+1)*2*sizeof(float),malign));
 for (i=0; i<2*(dist_order+1)*(dist_order+1); ++i) dist_param[i]=0.0;

 char coefficients[64]={"/tmp20/csoares/hmi_distortion/coef1766859/fitted-par_f09_c0.txt"};

 
 float dparam[32]={0.0,  0.0,  6.0592e-05, -8.3770e-05,  0.0, -3.4322e-04, 4.6960e-05,  0.0,  1.2345e-04, -2.0776e-04,  0.0,  0.0, -7.3936e-06,  0.0, 0.0, 0.0, 0.0, 0.0, -9.9498e-04,  2.0456e-04,  0.0,  3.1870e-05, -2.1903e-04,  0.0, -3.5076e-04,  5.2103e-05,  0.0,  0.0, -7.3271e-05, 0.0, 0.0, 0.0};
 
 for (i=0; i<32; ++i) dist_param[i]=dparam[i];
 

 int diffrot_order=2;
 float *diffrot_param;
 diffrot_param=(float *)(MKL_malloc((diffrot_order+1)*sizeof(float),malign));

 diffrot_param[0]=2.913-0.1991;  // Komm, Harvey, Howe (1994)
 diffrot_param[1]=-0.405;
 diffrot_param[2]=-0.422;

 const_param->dist_coef=dist_param;        //distortion parameters
 const_param->order_dist_coef=dist_order;    // order of distortion palynominal
 const_param->diffrot_coef=diffrot_param;    // differentail rotation coefficients
 const_param->order2_rot_coef=diffrot_order; // order of differental rotation polynominal in sin(lat)^2
 const_param->order_int=interpolation_order;                 // order for spatial  interpolation
 const_param->nconst=nconst;
 const_param->fills=fills;

 return status;
}





void free_interpol(struct initial *const_param)
{
  free_fill(&const_param->fills);
  MKL_free(const_param->dist_coef);
  MKL_free(const_param->diffrot_coef);
}


//int initialalize_gapfill(struct initial *const_param)
//{
//  int order_gapfill=const_param->gapfill_order;
//  int order_gapfill2=order_gapfill/2;
//  int status;
//  struct fill_struct fills;
//
// status=init_fill(const_param->gapfill_method, const_param->gapfill_regular, order_gapfill,order_gapfill2,order_gapfill2,&fills);
//}

int do_gapfill(float *image, unsigned char *mask, struct initial *const_param, float *ierror, int nx, int ny)
{

  int isample;
  int status;
  float t1;
  float *cnorm;

  int nlead=nx;
  int malign=32;

  cnorm=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));

  struct fill_struct fills=const_param->fills;

  status=fgap_fill(&fills,image,nx,ny,nlead,mask,cnorm,ierror);
   
  MKL_free(cnorm);


  return status;


}


int do_interpolate(float **images, float **ierrors, float *image_out, struct keyword *key, struct keyword *key_out, struct initial *const_param, int nsample, int nx, int ny)
{

  void derotation(float, float, float, float, float, float, float *, int, float *, int, int);
  long factorial(int);
  float intsincos(unsigned int n, unsigned int m);

  void derotation_grid(float *xin, float *yin, float radius, float cent_x, float cent_y, float dist, float p0, float b0, float *shift, int nx, int ny);
    

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

  //initialize
 
  int_images=(float **)(MKL_malloc(nsample*sizeof(float *), malign));
  masks=(unsigned char **)(MKL_malloc(nsample*sizeof(unsigned char *),malign));

  float *par;

  float *rrsun, *bb0, *pp0, *xx0, *yy0, *ddist;
  double *tsample;

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

  //

  //read out keywords
  for (isample=0; isample<nsample; ++isample){
    keyl=*(key+isample);
    rrsun[isample]=keyl.rsun;
    bb0[isample]=keyl.b0;
    pp0[isample]=keyl.p0;
    xx0[isample]=keyl.xx0;
    yy0[isample]=keyl.yy0;
    ddist[isample]=keyl.dist;
    tsample[isample]=keyl.time;
   }
 
  rsun=key_out->rsun;
  b0=key_out->b0;
  p0=key_out->p0;
  cent_x=key_out->xx0;
  cent_y=key_out->yy0;
  dist=key_out->dist;
  tint=key_out->time;
  //////////////////////////

  //float p0min=pp0[0],p0max=pp0[0],b0min=bb0[0],b0max=bb0[0],rsmin=rrsun[0],rsmax=rrsun[0],x0min=xx0[0],x0max=xx0[0],y0min=yy0[0],y0max=yy0[0];


  // obtain range of keywords
  //for (isample=1; isample < nsample; ++isample){
    
  //  p0min=minval(pp0[isample], p0min);
  //  p0max=maxval(pp0[isample], p0max);
  //  b0min=minval(bb0[isample], b0min);
  //  b0max=maxval(bb0[isample], b0max);
  //  rsmin=minval(rrsun[isample], rsmin);
  //  rsmax=maxval(rrsun[isample], rsmax);
  //  x0min=minval(xx0[isample], x0min);
  //  x0max=maxval(xx0[isample], x0max);
  //  y0min=minval(yy0[isample], y0min);
  //  y0max=maxval(yy0[isample], y0max);
  //  }

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

  ///////////

 
  for (isample=0; isample < nsample; ++isample) status_arr[isample]=0;
 

 

  t0=dsecnd();
  t1=dsecnd();

  float *xin, *yin, *shift;     //nput coordinated and shift

  xin=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
  yin=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
  shift=(float *)(MKL_malloc(2*nlead*ny*sizeof(float),malign));

  struct fint_struct fints, fints_error;
  
  
  int order_int=const_param->order_int;   //get interpolation order
  int order2_int=order_int/2;             // interpolation order /2

 

  init_finterpolate(&fints,fint_test);  //initialize interpolation
  init_finterpolate(&fints_error, fint_linear); // linear interpolation for masks


  float rad=(float)nx/2.0;   //normalization radius for distortion polynominal (2048)

  double time;  // image time - target time
 
  int maxo=const_param->order_dist_coef+1; //distortion polynominal order +1   



  float coef_x[maxo][maxo], coef_y[maxo][maxo];  //distortion coefficients in CCD coordinates
 
  printf("Distortion coefficient \n");

  for (s=0; s<maxo; ++s){for (t=0; t<maxo; ++t){coef_x[s][t]=const_param->dist_coef[s*maxo+t]; coef_y[s][t]=const_param->dist_coef[s*maxo+t+maxo*maxo];}}
  

  //calculate rotated coefficients
  ///////////////////////////////////

  float coefx_rot[nsample][maxo][maxo], coefy_rot[nsample][maxo][maxo];
  long bin_si, bin_tj, bin_sip, bin_tjp;
  float av_offset_x=0.0;
  float av_offset_y=0.0;
  int sp, tp, ip, jp;

  float rad_term1=0.0, rad_term2=0.0, rad_term3=0.0, rad_term4=0.0;

 

  //initailize distortion coefficients in target image coordinates

  for (isample=0; isample<nsample; ++isample) for (s=0; s<maxo; ++s) for (t=0; t<maxo; ++t){coefx_rot[isample][s][t]=0.0; coefy_rot[isample][s][t]=0.0;} //initialize
 

  float xc=(cent_x-(float)nx/2.0)/rad;   //average image center (centered and normalized)
  float yc=(cent_y-(float)ny/2.0)/rad;


  //calculate rsun-change due to removing distortion, and image center change due to removing distortion:
  //only output keywords are affected: Output rsun and image center should be center of undistorted image

    for (s=0; s<maxo; ++s) for (t=0; t<=(maxo-1-s); ++t) for (i=0; i<=s; ++i) for (j=0; j<=t; ++j){

      bin_si=(float)factorial(s)/factorial(i)/factorial(s-i);   
      bin_tj=(float)factorial(t)/factorial(j)/factorial(t-j);

      av_offset_x=av_offset_x+coef_x[s][t]*bin_si*bin_tj*pow(rsun/rad, s+t-(i+j))*pow(xc, i)*pow(yc, j)*intsincos(s-i,t-j)*rad; //calculate offs offset due to distortion
      av_offset_y=av_offset_y+coef_y[s][t]*bin_si*bin_tj*pow(rsun/rad, s+t-(i+j))*pow(xc, i)*pow(yc, j)*intsincos(s-i,t-j)*rad;
     
      rad_term1=rad_term1+2.0*coef_x[s][t]*bin_si*bin_tj*pow(rsun/rad,s+t-(i+j)+1)*pow(xc,i)*pow(yc,j)*intsincos(s-i+1,t-j);
      rad_term3=rad_term3+2.0*coef_y[s][t]*bin_si*bin_tj*pow(rsun/rad,s+t-(i+j)+1)*pow(xc,i)*pow(yc,j)*intsincos(s-i,t-j+1);

      for (sp=0; sp<maxo; ++sp) for (tp=0; tp<=(maxo-1-s); ++tp) for (ip=0; ip<=sp; ++ip) for (jp=0; jp<=tp; ++jp){

	bin_sip=(float)factorial(sp)/factorial(ip)/factorial(sp-ip);
	bin_tjp=(float)factorial(tp)/factorial(jp)/factorial(tp-jp);

	rad_term2=rad_term2+coef_x[s][t]*coef_x[sp][tp]*bin_si*bin_tj*bin_sip*bin_tjp*pow(rsun/rad,s+t+sp+tp-(i+j+ip+jp))*pow(xc, i+ip)*pow(yc, j+jp)*intsincos(s+sp-i-ip, t+tp-j-jp);
	rad_term4=rad_term4+coef_y[s][t]*coef_y[sp][tp]*bin_si*bin_tj*bin_sip*bin_tjp*pow(rsun/rad,s+t+sp+tp-(i+j+ip+jp))*pow(xc, i+ip)*pow(yc, j+jp)*intsincos(s+sp-i-ip, t+tp-j-jp);
      }
    }

    float dist_rad=sqrt(rsun*rsun-(rad_term1+rad_term2+rad_term3+rad_term4)*rad*rad); //calculate output rsun from different terms

   

    //calculate distortion coefficients in target image coordinates. Considers only p-angle rotation !!
     for (s=0; s<maxo; ++s) for (t=0; t<=(maxo-1-s); ++t) for (i=0; i<=s; ++i) for (j=0; j<=t; ++j){
       for (isample=0; isample<nsample; ++isample){
	coefx_rot[isample][s+t-(i+j)][i+j]=coefx_rot[isample][s+t-(i+j)][i+j]+coef_x[s][t]*bin_si*bin_tj*pow(cos(pp0[isample]-p0), s-i+j)*pow(sin(pp0[isample]-p0), i+t-j)*pow(-1.0, i); ;
	coefy_rot[isample][s+t-(i+j)][i+j]=coefy_rot[isample][s+t-(i+j)][i+j]+coef_y[s][t]*bin_si*bin_tj*pow(cos(pp0[isample]-p0), s-i+j)*pow(sin(pp0[isample]-p0), i+t-j)*pow(-1.0, i);
      }
     }
  

     //apply correction due to distortion;
     rsun=dist_rad;
     cent_x=cent_x+av_offset_x;
     cent_y=cent_y+av_offset_y;


     ///keywords for output image 
     // key_out->rsun=dist_rad; //note correction due to distortion
     //key_out->b0=b0;
     // key_out->p0=p0;                   //p-angle = averag p-angle (no "target" p-angle)
     //key_out->xx0=cent_x-av_offset_x;   //note correction due to distortion
     //key_out->yy0=cent_y-av_offset_y;   //note correction due to distortion
     // key_out->dist=dist;
     // key_out->time=tint;   //output time=target time

    ////////////////

  t1=dsecnd()-t1;



  t1=dsecnd();

  int order2_rot_coef=const_param->order2_rot_coef;   //squared order of differential rotation coefficients

  float *rot_coef;
  rot_coef=(float *)(MKL_malloc((order2_rot_coef+1)*sizeof(float), malign)); 

  for (i=0; i<=order2_rot_coef; ++i) rot_coef[i]=const_param->diffrot_coef[i];
      


  t1=dsecnd()-t1;

  
 float cosp0, sinp0;


 float *x1, *x2, *x3, *x1c;  //auxilliary valiables
 float y1,y2,y3;             //auxilliary valiables

 x1=(float *)MKL_malloc(nx*sizeof(float),malign);  //auxilliary valiables
 x2=(float *)MKL_malloc(nx*sizeof(float),malign);
 x3=(float *)MKL_malloc(nx*sizeof(float),malign);
 x1c=(float *)MKL_malloc(nx*sizeof(float),malign);
  
 float *dist_x, *dist_y;   //distortion vector field
 dist_x=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
 dist_y=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));

 for (i=0; i<nx; i++){x1[i]=((float)i-cent_x)/rad; x2[i]=x1[i]*x1[i]; x3[i]=x2[i]*x1[i];}

 derotation(rsun, cent_x, cent_y, dist, p0, b0, rot_coef, order2_rot_coef, shift, nx, ny); //calculate rotation pattern at CCD coordinate

 /////preparation for mask building
 float *ierror_out;   
 unsigned char *maskp;
 ierror_out=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));
 float thresold=2.0; // !! changed to 0.2 for test

 ///////////////////////////////////////////////////////////////
  for (isample=0; isample<nsample; ++isample){ //loop over images
 
    if (status_arr[isample] == 0){

    cosp0=cos(pp0[isample]-p0);
    sinp0=sin(pp0[isample]-p0);

    //calculate distortion vector field in target coordinates
    for (j=0; j<ny; j++){ //loop over target coordinates: calculate distortion shift
   
      y1=((float)j-cent_y)/rad; // centered at cent_y
      y2=y1*y1;
      y3=y2*y1;

     for (i=0; i<nx; i++){

       dist_x[j*nlead+i]=(coefx_rot[isample][2][0]*x2[i]+coefx_rot[isample][3][0]*x3[i]+coefx_rot[isample][1][1]*x1[i]*y1+coefx_rot[isample][2][1]*x2[i]*y1+coefx_rot[isample][0][2]*y2+coefx_rot[isample][1][2]*x1[i]*y2+coefx_rot[isample][0][3]*y3)*rad;
       dist_y[j*nlead+i]=(coefy_rot[isample][2][0]*x2[i]+coefy_rot[isample][3][0]*x3[i]+coefy_rot[isample][1][1]*x1[i]*y1+coefy_rot[isample][2][1]*x2[i]*y1+coefy_rot[isample][0][2]*y2+coefy_rot[isample][1][2]*x1[i]*y2+coefy_rot[isample][0][3]*y3)*rad;
     }
    }
    //////////////////////////////////////


   


   
    for (i=0; i<nx; i++) x1c[i]=((float)i-xx0[isample]);

    image_intout=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));

    time=(float)(tsample[isample]-tint); // difference to target time

    printf("time lag: %g \n", time);

   for (j=0; j<ny; j++){ //loop over target coordinates: calculate shift

     y1=-sinp0*((float)j-yy0[isample])+xx0[isample] + cosp0*(xx0[isample]-cent_x)-sinp0*(yy0[isample]-cent_y); //rotation around (xx0[isample], yy0[isample]) ,and shift to the center point (x term in xin[..]=.. )
     y2=cosp0*((float)j-yy0[isample])+yy0[isample] + sinp0*(xx0[isample]-cent_x)+cosp0*(yy0[isample]-cent_y);

     for (i=0; i<nx; i++){

       xin[j*nlead+i]=cosp0*x1c[i]+y1+ dist_x[j*nlead+i] + shift[j*nlead+i]*time; // distortion and rotational shift
       yin[j*nlead+i]=sinp0*x1c[i]+y2+ dist_y[j*nlead+i] + shift[ny*nlead+j*nlead+i]*time;  

      }
   }

   t1=dsecnd();

   status_arr[isample]=finterpolate(&fints,order_int,*(images+isample), xin, yin, image_intout, nx,ny,nlead,nx,ny,nlead,maxext,fillval);
   t1=dsecnd()-t1;  
   printf("spatial interpolation \t %f \n", t1);  

   
   *(int_images+isample)=image_intout;  //pointer to each output image
    

   ///////make masks out of errors
  
   finterpolate(&fints_error,1,*(ierrors+isample), xin, yin, ierror_out, nx,ny,nlead,nx,ny,nlead,maxext,2.0f);

   maskp=(unsigned char *)(MKL_malloc(nlead*ny*sizeof(unsigned char),malign));
   for (i=0; i<nx; ++i) for (j=0; j<ny; ++j) if (ierror_out[j*nlead+i] > thresold) maskp[j*nlead+i]=1; else maskp[j*nlead+i]=0;
   *(masks+isample)=maskp;
 

   
    }
  
  } // end loop over images

 
  
  


  //end interpolation///////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //count valid images
  int nvalid=0;
  printf("status: \t");
  for (isample=0; isample<nsample; ++isample){nvalid=nvalid+status_arr[isample];}
  nvalid=nsample-nvalid;

  printf("valid iamges %d\n", nvalid);

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
  printf("nvalid: %u %u \n", nvalid, count);

  t1=dsecnd();
  status=tinterpolate(nvalid,tvalid,tint,const_param->nconst,valid_images,valid_masks,image_out,nx,ny,nlead);  //call temporal interpolation function

  printf("status %u \n", status);
   t0=dsecnd()-t0;
  printf("Time to interpolate = %f seconds.\n",t0);


  printf("done\n");



	     
printf("D\n");
  // free memory
 MKL_FreeBuffers();
 
 for (isample=0;isample<nsample;isample++){if (status_arr[isample] ==0){MKL_free(*(int_images+isample)); MKL_free(*(masks+isample));}}
  
 free_finterpolate(&fints);
  
 free_finterpolate(&fints_error);
 
 MKL_free(int_images);
 
 MKL_free(masks);
 
 MKL_free(xin);
 
 MKL_free(yin);
 
 MKL_free(dist_x);
 
 MKL_free(dist_y);
 
 MKL_free(shift);
 
 MKL_free(x1c);
 
 MKL_free(x1);
 
 MKL_free(x2);
 
 MKL_free(x3);
 
 MKL_free(rot_coef);
 
 MKL_free(tsample);
 
 MKL_free(ddist);
 
 MKL_free(xx0);
 
 MKL_free(yy0);
 
 MKL_free(pp0);
 
 MKL_free(bb0);
 
 MKL_free(rrsun);
 
 MKL_free(status_arr);
 
 MKL_free(ierror_out);
 
 if (nvalid > 0){
   MKL_free(tvalid);
   MKL_free(valid_images);
   MKL_free(valid_masks);
 }
 ///////////////////////


 return status;  //status of temporal interpolation
 
 

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
   const float maxphi=atan(1.0/au);
 
   
   float sinp0=sin(p0);
   float cosp0=cos(p0);
   float sinb0=sin(b0);
   float cosb0=cos(b0);

   int nlead=nx;

   for (j=0; j<ny; ++j){

 
     for (i=0; i<ny; ++i){
   
       xy[0]=(cosp0*(i-cent_x)-sinp0*(j-cent_y))/radius;     
       xy[1]=(sinp0*(i-cent_x)+cosp0*(j-cent_y))/radius;  // normalized xy coordinate in the sky (solar radius=1)	

	  xyp[0]=cosp0*xy[0]+sinp0*xy[1];
	  xyp[1]=-sinp0*xy[0]+cosp0*xy[1];  // p-angle rotation

	  inr=sqrt(xyp[0]*xyp[0]+xyp[1]*xyp[1]); // radial coordinate in the sky
         
 
             
	  if (inr >= 1.0)
	    {
	      xyp[0]=xyp[0]/inr;
	      xyp[1]=xyp[1]/inr;  // force off-limb point to radial distance 1.0
	      inr=1.0;

	      phie=maxphi;
	      ind=(1.0-1.0/au/au)/(1.0+1.0/au/au);
            }
	  else
            {
	      phie=inr*maxphi;
	      beta=tan(phie);
              ind=(beta*au-beta*sqrt(1.0+beta*beta-au*au*beta*beta))/(1+beta*beta);
            }

	  if (inr != 0.0) ikf=ind/inr; else ikf=(au-1.0)/au; //limit for exact center


       	  xyz[0]=sqrt(1.0-ind*ind);
	  xyz[1]=xyp[0]*ikf;
	  xyz[2]=xyp[1]*ikf;

	  	      
	  slat=sinb0*xyz[0]+cosb0*xyz[2];  //sine of heliographic latitude

	   //omeg=2.0*M_PI*1.0e-9*(452.0- 49.0*pow(slat, 2) -84.0*pow(slat, 4) - 31.7)*time; // differential rotation rate [in rad/s]
           //omeg=1.0e-6*(2.838-0.301*pow(slat,2)-0.526*pow(slat, 4)-0.1991)*time; //Snodgrass 1984a
	  //omeg=1e-6*(2.913-0.405*pow(slat,2)-0.422*pow(slat,4)-0.1991); //Komm, Howard, Harvey
	
	  omeg=0.0;
	  for (k=0; k<=order2_rot_coef; ++k) omeg=omeg+1e-6*rot_coef[k]*pow(slat,2*k); 

	  sxyz[0]=-omeg*xyz[1]*cosb0;
	  sxyz[1]=omeg*(xyz[0]*cosb0-xyz[2]*sinb0); // linearized rotation
	  sxyz[2]=omeg*xyz[1]*sinb0; //

	  //singam=sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]); // projection from sphere onto observing plane
	  //cosgam=xyz[0];
	  //phi=atan(singam/(au-cosgam));
	  

	  //if (singam !=0.0) ing=1./ikf; else ing=1.0;  //limit for exact center

          sxyp[0]=sxyz[1]/ikf;
	  sxyp[1]=sxyz[2]/ikf;

	  sxy[0]=cosp0*sxyp[0]-sinp0*sxyp[1];
	  sxy[1]=sinp0*sxyp[0]+cosp0*sxyp[1];  // inverse p-angle rotation 

	  


	  shift[j*nlead+i]=sxy[0]*radius;
	  shift[ny*nlead+j*nlead+i]=sxy[1]*radius;
	   
     }
   }

}


void derotation_grid(float *xin, float *yin, float radius, float cent_x, float cent_y, float dist, float p0, float b0, float *shift, int nx, int ny){
    

  int i, j;
  float xy[2], xyp[2];
  float xyz[3], sxyz[3];
  float sxyp[2], sxy[2];
  float inr, ind, ikf, phie, ing, beta;

   
  float slat, omeg;
  float singam, cosgam, phi; 
      
    
   const float au=215.020*dist; // (1 astonomical unit in solar radii) au=149597870.691 km, R_SUN=695740 km (Kuhn 2004)
   const float maxphi=atan(1.0/au);
 
   int nlead=nx;


   for (j=0; j<ny; ++j){
     for (i=0; i<nx; ++i){
   
	 xy[0]=(xin[j*nx+i]-cent_x)/radius;
	 xy[1]=(yin[j*nx+i]-cent_y)/radius;  // normalized xy coordinate in the sky (solar radius=1)

	  xyp[0]=cos(p0)*xy[0]+sin(p0)*xy[1];
	  xyp[1]=-sin(p0)*xy[0]+cos(p0)*xy[1];  // p-angle rotation

	  inr=sqrt(xyp[0]*xyp[0]+xyp[1]*xyp[1]); // radial coordinate in the sky
         
 
             
	  if (inr >= 1.0)
	    {
	      xyp[0]=xyp[0]/inr;
	      xyp[1]=xyp[1]/inr;  // force off-limb point to radial distance 1.0
	      inr=1.0;

	      phie=maxphi;
	      ind=(1.0-1.0/au/au)/(1.0+1.0/au/au);
            }
	  else
            {
	      phie=inr*maxphi;
	      beta=tan(phie);
              ind=(beta*au-beta*sqrt(1.0+beta*beta-au*au*beta*beta))/(1+beta*beta);
            }

	  if (inr != 0.0) ikf=ind/inr; else ikf=(au-1.0)/au; //limit for exact center


       	  xyz[0]=sqrt(1.0-ind*ind);
	  xyz[1]=xyp[0]*ikf;
	  xyz[2]=xyp[1]*ikf;

	  	      
	  slat=sin(b0)*xyz[0]+cos(b0)*xyz[2];  //sine of heliographic latitude

	   //omeg=2.0*M_PI*1.0e-9*(452.0- 49.0*pow(slat, 2) -84.0*pow(slat, 4) - 31.7)*time; // differential rotation rate [in rad/s]
           //omeg=1.0e-6*(2.838-0.301*pow(slat,2)-0.526*pow(slat, 4)-0.1991)*time; //Snodgrass 1984a
	  omeg=1e-6*(2.913-0.405*pow(slat,2)-0.422*pow(slat,4)-0.1991); //Komm, Howard, Harvey

	  sxyz[0]=-omeg*xyz[1]*cos(b0);
	  sxyz[1]=omeg*(xyz[0]*cos(b0)-xyz[2]*sin(b0)); // linearized rotation
	  sxyz[2]=omeg*xyz[1]*sin(b0); //

	  //singam=sqrt(xyz[1]*xyz[1]+xyz[2]*xyz[2]); // projection from sphere onto observing plane
	  //cosgam=xyz[0];
	  //phi=atan(singam/(au-cosgam));
	  

	  //if (singam !=0.0) ing=1./ikf; else ing=1.0;  //limit for exact center

          sxyp[0]=sxyz[1]/ikf;
	  sxyp[1]=sxyz[2]/ikf;

	  sxy[0]=cos(p0)*sxyp[0]-sin(p0)*sxyp[1];
	  sxy[1]=sin(p0)*sxyp[0]+cos(p0)*sxyp[1];  // inverse p-angle rotation 

	  


	  shift[j*nx+i]=sxy[0]*radius;
	  shift[ny*nx+j*nx+i]=sxy[1]*radius;
	   
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


//icc -O3 -o target_test target_test.c tinterpolate.c finterpolate.c cholesky_down.c gapfill.c -lmkl_em64t -lguide -lpthread
