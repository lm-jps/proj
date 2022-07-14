#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include "finterpolate.h"
#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) < (y)) ? (y) : (x))

int init_finterpolate(
  struct fint_struct *pars,
  int method
)
{
  const int malign=32;
  const int nsubin0=100,nsubin=nsubin0+1,maxorder=32;
  int i,j,k;
  double *kers;
  FILE *fileptr;

  pars->method=method;

  switch (method) {
  case fint_test:

    pars->maxorder=maxorder;
    pars->fover=nsubin0;
    pars->nsub=nsubin0+1;
  
    fileptr = fopen ("/tmp20/schou/coeff3_0-32_100_double.bin", "r");
    if (fileptr == NULL) {
      printf("File not found in finterpolate!\n");
      return 2;
    }
    kers=(double *)(MKL_malloc(nsubin*maxorder*(maxorder+1)*sizeof(double),malign));
    fread ((char*)kers,sizeof(double),(maxorder+1)*nsubin*maxorder,fileptr);
    fclose(fileptr);
  
    pars->kersx=(float *)(MKL_malloc(nsubin*maxorder*(maxorder+1)*sizeof(float),malign));
    for (i=0;i<=maxorder;i++) {
      for (j=0;j<nsubin;j++) {
        for (k=0;k<maxorder;k++) {
          pars->kersx[i*nsubin*maxorder+j*maxorder+k]=kers[i*nsubin*maxorder+j*maxorder+k];
        }
      }
    }
    MKL_free(kers);
    break;

  case fint_linear:
    break;

  case fint_cubic_conv:
    break;

  default:
    printf("Unimplemented method in init_finterpolate\n");
    return 1;
    break;

  } //switch
  return 0;
}

int free_finterpolate(
  struct fint_struct *pars
)
{
  if (pars->method==fint_test) MKL_free (pars->kersx);

  return 0;
}

int winterpolate(
struct fint_struct *pars,
int order,
float *image_in,
float *xin,
float *yin,
float *image_out,
int nxin,
int nyin,
int nleadin,
int nx,
int ny,
int nlead
)

{
  int i,j,i1,j1;
  int order2,lorder;
  float *xker,*yker;
  const int ione = 1;
  int malign=32;
  int leaddiv=4;
  int ixin,iyin,ixin1,iyin1;
  float *ixins,*iyins,*ixin1s,*iyin1s;
  float fxin1,fyin1,fxin2,fyin2;
  float *fxins,*fyins,*xinp,*yinp,*helpf,*fxin1s,*fyin1s;
  float *xk1,*xk2,*yk1,*yk2;
  float *imp;
  float sum,sum1;
  float *kerso;

  ixins=(float *)(MKL_malloc(nx*sizeof(int),malign));
  iyins=(float *)(MKL_malloc(nx*sizeof(int),malign));
  fxins=(float *)(MKL_malloc(nx*sizeof(float),malign));
  fyins=(float *)(MKL_malloc(nx*sizeof(float),malign));
  ixin1s=(float *)(MKL_malloc(nx*sizeof(float),malign));
  iyin1s=(float *)(MKL_malloc(nx*sizeof(float),malign));
  fxin1s=(float *)(MKL_malloc(nx*sizeof(float),malign));
  fyin1s=(float *)(MKL_malloc(nx*sizeof(float),malign));
  helpf=(float *)(MKL_malloc(nx*sizeof(float),malign));

  xker=(float *)(MKL_malloc(order*sizeof(float),malign));
  yker=(float *)(MKL_malloc(order*sizeof(float),malign));

  order2=order/2;
  lorder=leaddiv*((order-1)/leaddiv+1);
  kerso=pars->kersx+order*(pars->maxorder)*(pars->nsub);


// This awful code basically hardcodes the common cases and allows the
// compiler to optimize the code for a significant speed improvement
  switch(order) {
  case 2:
#define OOO1 2
#define OOO2 1
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 4:
#define OOO1 4
#define OOO2 2
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 6:
#define OOO1 6
#define OOO2 3
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 8:
#define OOO1 8
#define OOO2 4
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 10:
#define OOO1 10
#define OOO2 5
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 12:
#define OOO1 12
#define OOO2 6
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 14:
#define OOO1 14
#define OOO2 7
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 16:
#define OOO1 16
#define OOO2 8
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 18:
#define OOO1 18
#define OOO2 9
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 20:
#define OOO1 20
#define OOO2 10
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 22:
#define OOO1 22
#define OOO2 11
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 24:
#define OOO1 24
#define OOO2 12
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 26:
#define OOO1 26
#define OOO2 13
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 28:
#define OOO1 28
#define OOO2 14
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 30:
#define OOO1 30
#define OOO2 15
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  case 32:
#define OOO1 32
#define OOO2 16
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
// If not hardcoded do general calculation
  default:
#define OOO1 order
#define OOO2 order2
#include "finterpolate.include"
#undef OOO1
#undef OOO2
    break;
  }

  MKL_free(ixins);
  MKL_free(iyins);
  MKL_free(fxins);
  MKL_free(fyins);
  MKL_free(ixin1s);
  MKL_free(iyin1s);
  MKL_free(fxin1s);
  MKL_free(fyin1s);
  MKL_free(helpf);
  MKL_free(xker);
  MKL_free(yker);

  return 0;
}

int linterpolate(
float *image_in,
float *xin,
float *yin,
float *image_out,
int nxin,
int nyin,
int nleadin,
int nx,
int ny,
int nlead,
float maxext,
float fillval
)

{
  int i,j;
  int malign=32;
  int ixin,iyin;
  float *ixins,*iyins;
  float fxin1,fyin1,fxin2,fyin2;
  float *xinp,*yinp;
  float *imp;
  float xmin,xmax,ymin,ymax;

  ixins=(float *)(MKL_malloc(nx*sizeof(int),malign));
  iyins=(float *)(MKL_malloc(nx*sizeof(int),malign));

  xmin=-maxext;
  xmax=nxin-1.0f+maxext;
  ymin=-maxext;
  ymax=nyin-1.0f+maxext;

  for (j=0; j<ny; j++) {
    xinp=xin+j*nlead;
    vsFloor(nx,xinp,ixins);
    yinp=yin+j*nlead;
    vsFloor(nx,yinp,iyins);
    for (i=0; i<nx; i++) {
      ixin=(int)ixins[i]; // Integer pixel to interpolate to
      iyin=(int)iyins[i]; // Integer pixel to interpolate to
      if ((ixin>=0) && (ixin <(nxin-1)) && (iyin>=0) && (iyin <(nyin-1))) {
// Normal case

        fxin1=xinp[i]-ixin;
        fyin1=yinp[i]-iyin;
        fxin2=1.0f-fxin1;
        fyin2=1.0f-fyin1;
  
        imp=image_in+ixin+nleadin*iyin;

/* Brute force addition */
        image_out[i+nlead*j]=fyin2*(fxin2*imp[0]+fxin1*imp[1])+
                             fyin1*(fxin2*imp[nleadin]+fxin1*imp[nleadin+1]);
      } // if
      else {
        if ((xinp[i]>=xmin) && (xinp[i]<=xmax) && (yinp[i]>=ymin) && (yinp[i]<=ymax)) {
// Extrapolating or point exactly at edge
// Move ixin and iyin inside of array
          ixin=maxval(0,minval(ixin,nxin-2));
          iyin=maxval(0,minval(iyin,nyin-2));

          fxin1=xinp[i]-ixin;
          fyin1=yinp[i]-iyin;
          fxin2=1.0f-fxin1;
          fyin2=1.0f-fyin1;
    
          imp=image_in+ixin+nleadin*iyin;
  
  /* Brute force addition */
          image_out[i+nlead*j]=fyin2*(fxin2*imp[0]+fxin1*imp[1])+
                               fyin1*(fxin2*imp[nleadin]+fxin1*imp[nleadin+1]);
        }
        else {
// Outside of array
          image_out[i+nlead*j]=fillval;
        }
      }
    } // i=
  } //j=

  MKL_free(ixins);
  MKL_free(iyins);

  return 0;
}

int ccinterpolate( // Cubic convolution interpolation
float *image_in,
float *xin,
float *yin,
float *image_out,
int nxin,
int nyin,
int nleadin,
int nx,
int ny,
int nlead,
float maxext,
float fillval
)

{
  int i,j,i1,j1;
  float xker[4],yker[4];
  int malign=32;
  int ixin,iyin;
  float *ixins,*iyins;
  float *xinp,*yinp;
  float *imp;
  float sum,sum1;
  float s,s2;
  float xmin,xmax,ymin,ymax;

  ixins=(float *)(MKL_malloc(nx*sizeof(int),malign));
  iyins=(float *)(MKL_malloc(nx*sizeof(int),malign));

  xmin=1-maxext;
  xmax=nxin-2.0f+maxext;
  ymin=1-maxext;
  ymax=nyin-2.0f+maxext;

  for (j=0; j<ny; j++) {
    xinp=xin+j*nlead;
    vsFloor(nx,xinp,ixins);
    yinp=yin+j*nlead;
    vsFloor(nx,yinp,iyins);
    for (i=0; i<nx; i++) {
      ixin=(int)ixins[i]; // Integer pixel to interpolate to
      iyin=(int)iyins[i]; // Integer pixel to interpolate to
      if ((ixin>=1) && (ixin <(nxin-2)) && (iyin>=1) && (iyin <(nyin-2))) {

// RML kernel calculation based on fc3kernel in winterpolate_old.c
// Kernels are unnormalized and need to be divided by 2 

        s=xinp[i]-ixin;
        s2 = s*s;
        xker[0] = -s*(1.0f - s*(2.0f - s));
        xker[1] = 2.0f - s2*(5.0f - 3.0f*s);
        xker[2] = s*(1.0f + s*(4.0f - 3.0f*s));
        xker[3] = -s2*(1.0f-s);

        s=yinp[i]-iyin;
        s2 = s*s;
        yker[0] = -s*(1.0f - s*(2.0f - s));
        yker[1] = 2.0f - s2*(5.0f - 3.0f*s);
        yker[2] = s*(1.0f + s*(4.0f - 3.0f*s));
        yker[3] = -s2*(1.0f-s);
  
        imp=image_in+ixin-1+nleadin*(iyin-1);

/* Brute force addition */
        sum=0.0f;
        for (i1=0; i1<4; i1++) {
          sum=sum+yker[i1]*(xker[0]*imp[0]+xker[1]*imp[1]+xker[2]*imp[2]+xker[3]*imp[3]);
          imp=imp+nleadin;
        }
        image_out[i+nlead*j]=sum/4;
      } // if
      else {
        if ((xinp[i]>=xmin) && (xinp[i]<=xmax) && (yinp[i]>=ymin) && (yinp[i]<=ymax)) {
// Extrapolating or point exactly at edge
// Move ixin and iyin inside of array
          ixin=maxval(1,minval(ixin,nxin-3));
          iyin=maxval(1,minval(iyin,nyin-3));
// RML kernel calculation based on fc3kernel in winterpolate_old.c
// Kernels are unnormalized and need to be divided by 2 

          s=xinp[i]-ixin;
          s2 = s*s;
          xker[0] = -s*(1.0f - s*(2.0f - s));
          xker[1] = 2.0f - s2*(5.0f - 3.0f*s);
          xker[2] = s*(1.0f + s*(4.0f - 3.0f*s));
          xker[3] = -s2*(1.0f-s);
  
          s=yinp[i]-iyin;
          s2 = s*s;
          yker[0] = -s*(1.0f - s*(2.0f - s));
          yker[1] = 2.0f - s2*(5.0f - 3.0f*s);
          yker[2] = s*(1.0f + s*(4.0f - 3.0f*s));
          yker[3] = -s2*(1.0f-s);
    
          imp=image_in+ixin-1+nleadin*(iyin-1);
  
  /* Brute force addition */
          sum=0.0f;
          for (i1=0; i1<4; i1++) {
            sum=sum+yker[i1]*(xker[0]*imp[0]+xker[1]*imp[1]+xker[2]*imp[2]+xker[3]*imp[3]);
            imp=imp+nleadin;
          }
          image_out[i+nlead*j]=sum/4;
        }
        else {
// Outside of array
          image_out[i+nlead*j]=fillval;
        }
      }
// Ought to do something here.
    } // i=
  } //j=

  MKL_free(ixins);
  MKL_free(iyins);

  return 0;
}

int finterpolate(
struct fint_struct *pars,
int order,
float *image_in,
float *xin,
float *yin,
float *image_out,
int nxin,
int nyin,
int nleadin,
int nx,
int ny,
int nlead,
float maxext,
float fillval
)
{
  int status;

  switch (pars->method) {
  case fint_test:
    status=winterpolate(pars,order,image_in,xin,yin,image_out,nxin,nyin,nleadin,nx,ny,nlead);
    break;

  case fint_linear:
    status=linterpolate(image_in,xin,yin,image_out,nxin,nyin,nleadin,nx,ny,nlead,maxext,fillval);
    break;

  case fint_cubic_conv:
    status=ccinterpolate(image_in,xin,yin,image_out,nxin,nyin,nleadin,nx,ny,nlead,maxext,fillval);
    break;

  default:
    status=1;
    break;
  }

return status;

}

