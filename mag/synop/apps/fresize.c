#include <stdio.h>

#include <stdlib.h>

#include <time.h>

#include <string.h>

#include <math.h>

#include <mkl_blas.h>

#include <mkl_service.h>

#include <mkl_lapack.h>

#include <mkl_vml_functions.h>

#include <omp.h>

#include "fresize.h"

#define minval(x,y) (((x) < (y)) ? (x) : (y))

#define maxval(x,y) (((x) < (y)) ? (y) : (x))

#define fresize_sample 1

#define fresize_bin 2

#define fresize_1d 3

#define fresize_2d 4



double sinc(double x)

{

if (fabs(x) < (1.e-10))

  return 1.;

else

  return sin(M_PI*x)/(M_PI*x);

}



int init_fresize_sample(

  struct fresize_struct *pars,

  int nsub

)

{



  pars->method=fresize_sample;

  pars->nsub=nsub;



  return 0;

}



int init_fresize_bin(

  struct fresize_struct *pars,

  int nsub

)

{



  pars->method=fresize_bin;

  pars->nsub=nsub;



  return 0;

}



int init_fresize_boxcar(

  struct fresize_struct *pars,

  int hwidth, // Half width of boxcar. Full is 2*hwidth+1.

  int nsub // Distance between sampled points

)

{

  const int malign=32;

  int fwidth;

  int i,j;



  pars->method=fresize_1d;

  pars->nsub=nsub;

  pars->hwidth=hwidth;

  fwidth=2*hwidth+1;

  pars->kerx=(float *)(MKL_malloc(fwidth*sizeof(float),malign));

  pars->kery=(float *)(MKL_malloc(fwidth*sizeof(float),malign));

  for (i=0;i<fwidth;i++) {

    pars->kerx[i]=1.0f/fwidth;

    pars->kery[i]=1.0f/fwidth;

  }



  return 0;

}



int init_fresize_gaussian(

  struct fresize_struct *pars,

  float sigma, // Shape is exp(-(d/sigma)^2/2)

  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.

  int nsub // Distance between sampled points

)

{

  const int malign=32;

  int fwidth;

  int i,j;

  float sum,x,y;



  pars->method=fresize_1d;

  pars->nsub=nsub;

  pars->hwidth=hwidth;

  fwidth=2*hwidth+1;

  pars->kerx=(float *)(MKL_malloc(fwidth*sizeof(float),malign));

  pars->kery=(float *)(MKL_malloc(fwidth*sizeof(float),malign));

  sum=0.0f;

  for (i=0;i<fwidth;i++) {

    x=(i-hwidth)/sigma;

    y=exp(-x*x/2);

    sum=sum+y;

  }

  for (i=0;i<fwidth;i++) {

    x=(i-hwidth)/sigma;

    y=exp(-x*x/2);

    pars->kerx[i]=y/sum;

    pars->kery[i]=y/sum;

  }



  return 0;

}



int free_fresize(

  struct fresize_struct *pars

)

{

  if (pars->method==fresize_1d) {

    MKL_free (pars->kerx);

    MKL_free (pars->kery);

  }



  return 0;

}



int fsample(

  float *image_in,

  float *image_out,

  int nxin,

  int nyin,

  int nleadin,

  int nxout,

  int nyout,

  int nleadout,

  int nsub,

  int xoff,

  int yoff,

  float fillval

)



{

  int i,j;

  int imin,imax,jmin,jmax;

  float *imp;



// Find first and last indices for which resulting input point is within image

  if (xoff >= 0) imin=0; else imin=((-xoff+nsub-1)/nsub);

  if (xoff <= 0) imax=nxout-1; else imax=minval(((nxin-xoff-1)/nsub),nxout-1);

  if (yoff >= 0) jmin=0; else jmin=((-yoff+nsub-1)/nsub);

  if (yoff <= 0) jmax=nyout-1; else jmax=minval(((nyin-yoff-1)/nsub),nyout-1);



  imp=image_in+yoff*nleadin+xoff;



#pragma omp parallel default(none) private(i,j) shared(image_in,image_out,imp,fillval) shared(imin,imax,jmin,jmax) shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)

  { // Needed to define parallel region



// Fill below valid region

#pragma omp for

    for (j=0; j<jmin; j++) {

      for (i=0; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



// Valid heights

#pragma omp for

    for (j=jmin; j<=jmax; j++) {

// Fill to the left

      for (i=0; i<imin; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

// Valid region

      for (i=imin; i<=imax; i++) {

//      image_out[j*nleadout+i]=image_in[j*nsub*nleadin+i*nsub];

        image_out[j*nleadout+i]=imp[j*nsub*nleadin+i*nsub];

      } // i=

// Fill to the right

      for (i=imax+1; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



// Fill above valid region

#pragma omp for

    for (j=jmax+1; j<nyout; j++) {

      for (i=0; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



  } // End of parallel region



  return 0;

}



int fbin(

  float *image_in,

  float *image_out,

  int nxin,

  int nyin,

  int nleadin,

  int nxout,

  int nyout,

  int nleadout,

  int nsub,

  int xoff,

  int yoff,

  float fillval

)



{

  int i,j,i1,j1;

  int imin,imax,jmin,jmax;

  float *imp,*impi;

  float sum;



// Find first and last indices for which resulting input point is within image

  if (xoff >= 0) imin=0; else imin=((-xoff+nsub-1)/nsub);

  if (xoff <= 0) imax=nxout-1; else imax=minval(((nxin-xoff-nsub)/nsub),nxout-1);

  if (yoff >= 0) jmin=0; else jmin=((-yoff+nsub-1)/nsub);

  if (yoff <= 0) jmax=nyout-1; else jmax=minval(((nyin-yoff-nsub)/nsub),nyout-1);



  imp=image_in+yoff*nleadin+xoff;



#pragma omp parallel default(none) private(i,j,i1,j1,impi,sum) shared(image_in,image_out,imp,fillval) shared(imin,imax,jmin,jmax) shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
  { // Needed to define parallel region



// Fill below valid region

#pragma omp for

    for (j=0; j<jmin; j++) {

      for (i=0; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



// Valid heights

#pragma omp for

    for (j=jmin; j<=jmax; j++) {

// Fill to the left

      for (i=0; i<imin; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

// Valid region

      for (i=imin; i<=imax; i++) {

        impi=imp+j*nleadin*nsub+i*nsub;

        sum=0.0f;

        for (j1=0; j1<nsub; j1++) {

          for (i1=0; i1<nsub; i1++) {

            sum=sum+impi[i1];

          }

          impi=impi+nleadin;

        }

        image_out[j*nleadout+i]=sum/(nsub*nsub);

      } // i=

// Fill to the right

      for (i=imax+1; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



// Fill above valid region

#pragma omp for

    for (j=jmax+1; j<nyout; j++) {

      for (i=0; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



  } // End of parallel region



  return 0;

}



int f1d(

  struct fresize_struct *pars,

  float *image_in,

  float *image_out,

  int nxin,

  int nyin,

  int nleadin,

  int nxout,

  int nyout,

  int nleadout,

  int xoff,

  int yoff,

  float fillval

)



{

  const int malign=32;

  char transpose[] = "t";

  char normal[] = "n";

  const int ione = 1;

  const float fone = 1.0f;

  const float fzero = 0.0f;

  int i,j,i1,j1;

  int imin,imax,jmin,jmax;

  float *inp,*inpi,*inpj;

  float sum;

  int nsub,hwidth;

  float *kerx,*kery,*work;

  int xoffl,xoffh,yoffl,yoffh;

  int nwork;

  double t1,t2,t3;

  int n1,n2,nchunk;



  nsub=pars->nsub;

  hwidth=pars->hwidth;

  kerx=pars->kerx;

  kery=pars->kery;



// Find first and last indices for which resulting input point is within image

  xoffl=xoff-hwidth;

  xoffh=xoff+hwidth;

  if (xoffl >= 0) imin=0; else imin=((-xoffl+nsub-1)/nsub);

  if (xoffh <= 0) imax=nxout-1; else imax=minval(((nxin-xoffh-1)/nsub),nxout-1);



  yoffl=yoff-hwidth;

  yoffh=yoff+hwidth;

  if (yoffl >= 0) jmin=0; else jmin=((-yoffl+nsub-1)/nsub);

  if (yoffh <= 0) jmax=nyout-1; else jmax=minval(((nyin-yoffh-1)/nsub),nyout-1);



// Get work array big enough to hold convolution in x

  nwork=nxout*nyin;

  work=(float *)(MKL_malloc(nwork*sizeof(float),malign));



nchunk=64;

n1=(imax-imin+1); // Size of matrix for sgemv

n2=2*hwidth+1; // Size of matrix for sgemv

t1=dsecnd();

 

#pragma omp parallel default(none) private(i,j,i1,j1,inpi,inpj,sum) shared(image_in,image_out,fillval) shared(imin,imax,jmin,jmax,hwidth,kerx,kery,work,xoffl,yoffl) shared(normal,transpose,n1,n2,nchunk) shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
  { // Needed to define parallel region



// First convolve in x

// Brute force loop takes longer than it ought to, but have not

// found an obvious solution.

#pragma omp for

    for (j=0;j<nyin;j++) {

      inpj=image_in+j*nleadin+xoffl; // Note offset to match kernels.

      for (i=imin; i<=imax; i++) {

        sum=0.0f;

        inpi=inpj+i*nsub;

        for (i1=0; i1<=2*hwidth; i1++) {

          sum=sum+kerx[i1]*inpi[i1];

        }

        work[j*nleadout+i]=sum;

// This works but is slower

//      work[j*nleadout+i]=sdot(&n2,kerx,&ione,inpj+i*nsub,&ione);

      }

// This clever trick gives error message

//sgemv(transpose,&n2,&n1,&fone,inpj+imin*nsub,&nsub,kerx,&ione,&fzero,work+j*nleadout+imin,&ione);

    }

/*

// This also works but is slower

#pragma omp for

for (j=0;j<nyin;j=j+nchunk) {

  for (i=imin; i<=imax; i++) {

    sgemv(transpose,&n2,&nchunk,&fone,image_in+j*nleadin+xoffl+i*nsub,&nleadin,kerx,&ione,&fzero,work+j*nleadout+i,&nleadout);

}

*/



// Then convolve in y



#pragma omp for

    for (j=jmin; j<=jmax; j++) {

      inpj=work+(yoffl+j*nsub)*nleadout;

/* Old brute force code

      for (i=imin; i<=imax; i++) {

        sum=0.0f;

        for (j1=0; j1<=2*hwidth; j1++) {

          sum=sum+inpj[j1*nleadout+i]*kery[j1];

        }

        image_out[j*nleadout+i]=sum;

      }

*/

      sgemv(normal,&n1,&n2,&fone,inpj+imin,&nleadout,kery,&ione,&fzero,image_out+j*nleadout+imin,&ione);

    }



// Fill below valid region

#pragma omp for

    for (j=0; j<jmin; j++) {

      for (i=0; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



// Valid heights

#pragma omp for

    for (j=jmin; j<=jmax; j++) {

// Fill to the left

      for (i=0; i<imin; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

// Fill to the right

      for (i=imax+1; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



// Fill above valid region

#pragma omp for

    for (j=jmax+1; j<nyout; j++) {

      for (i=0; i<nxout; i++) {

        image_out[j*nleadout+i]=fillval;

      } // i=

    } //j=



  } // End of parallel region

  MKL_free(work);

  return 0;

}





int fresize(

  struct fresize_struct *pars,

  float *image_in,

  float *image_out,

  int nxin,

  int nyin,

  int nleadin,

  int nx,

  int ny,

  int nlead,

  int xoff,

  int yoff,

  float fillval

)

{

  int status;



  switch (pars->method) {



  case fresize_sample:

    status=fsample(image_in,image_out,nxin,nyin,nleadin,nx,ny,nlead,pars->nsub,xoff,yoff,fillval);

    break;



  case fresize_bin:

    status=fbin(image_in,image_out,nxin,nyin,nleadin,nx,ny,nlead,pars->nsub,xoff,yoff,fillval);

    break;



  case fresize_1d:

    status=f1d(pars,image_in,image_out,nxin,nyin,nleadin,nx,ny,nlead,xoff,yoff,fillval);

    break;



  default:

    status=1;

    break;

  }



return status;



}




