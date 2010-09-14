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
#include <complex.h>
#include <fftw3.h>
#include "fresize.h"
#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) < (y)) ? (y) : (x))
#define fresize_sample 1
#define fresize_bin 2
#define fresize_1d 3
#define fresize_2d 4
#define fresize_1d_fft 5
#define fresize_2d_fft 6

static double sinc(double x)
{
if (fabs(x) < (1.e-10))
  return 1.;
else
  return sin(M_PI*x)/(M_PI*x);
}

int make_fft1d(
  struct fresize_struct *pars,
  int nxin,
  int nyin
)
{
  fftwf_complex *helpc,*fkernel,*fkernely;
  float *helpin,*helpout;
  fftwf_plan plan1,plan2,plan1y,plan2y;
  int nxinp,nyinp; // Complex array size
  int nin,ninp;
  int hwidth /*,fwidth */;
  int i,/*j,*/i1/*,j1*/;
  float c;

  hwidth=pars->hwidth;
  /* fwidth=2*hwidth+1; */
  nxinp=(nxin/2+1);
  nyinp=(nyin/2+1);
  nin=maxval(nxin,nyin); // Max value to use same arrays
  ninp=maxval(nxinp,nyinp); // Max value to use same arrays

  helpin=(float*) fftwf_malloc(sizeof(float) * nin);
  fkernel=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxinp);
  fkernely=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nyinp);
  helpc=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ninp);
  helpout=(float*) fftwf_malloc(sizeof(float) * nin);

  plan1=fftwf_plan_dft_r2c_1d(nxin,helpin,helpc,FFTW_ESTIMATE);
  plan2=fftwf_plan_dft_c2r_1d(nxin,helpc,helpout,FFTW_ESTIMATE);
  plan1y=fftwf_plan_dft_r2c_1d(nxin,helpin,helpc,FFTW_ESTIMATE);
  plan2y=fftwf_plan_dft_c2r_1d(nxin,helpc,helpout,FFTW_ESTIMATE);

  pars->helpin=helpin;
  pars->helpout=helpout;
  pars->fkernel=fkernel;
  pars->fkernely=fkernely;
  pars->helpc=helpc;
  pars->plan1=plan1;
  pars->plan2=plan2;
  pars->plan1y=plan1y;
  pars->plan2y=plan2y;
  pars->nxin=nxin;
  pars->nyin=nyin;
  pars->nxinp=nxinp;
  pars->nyinp=nyinp;
  pars->method=fresize_1d_fft;

/* First transform kernel in x*/

/* Zero entire array */
  for (i=0;i<nxin;i++) {
    helpin[i]=0;
  }

// Copy kernel to new array
  for (i=-hwidth;i<=hwidth;i++) {
    i1=(i+nxin) % nxin;
    helpin[i1]=pars->kerx[i+hwidth];
  }

/* Transform kernel */
  fftwf_execute(plan1);

/* Copy transform and scale */
  c=1.0f/nxin;

  for (i=0;i<nxinp;i++) {
    fkernel[i]=c*helpc[i];
  }

/* Then transform kernel in y*/

/* Zero entire array */
  for (i=0;i<nxin;i++) {
    helpin[i]=0;
  }

/* Copy kernel to new array */
  for (i=-hwidth;i<=hwidth;i++) {
    i1=(i+nyin) % nyin;
    helpin[i1]=pars->kery[i+hwidth];
  }

/* Transform kernel */
  fftwf_execute(plan1y);

/* Copy transform and scale */
  c=1.0f/nyin;

  for (i=0;i<nyinp;i++) {
    fkernely[i]=c*helpc[i];
  }

  return 0;
}

int make_fft2d(
  struct fresize_struct *pars,
  int nxin,
  int nyin
)
{
  fftwf_complex *helpc,*fkernel;
  float *helpin,*helpout;
  fftwf_plan plan1,plan2;
  int nxinp; // Complex array size
  int hwidth,fwidth;
  int i,j,i1,j1;
  float c;

  hwidth=pars->hwidth;
  fwidth=2*hwidth+1;
  nxinp=(nxin/2+1);

  helpin=(float*) fftwf_malloc(sizeof(float) * nxin*nyin);
  fkernel=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxinp*nyin);
  helpc=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxinp*nyin);
  helpout=(float*) fftwf_malloc(sizeof(float) * nxin*nyin);

  plan1=fftwf_plan_dft_r2c_2d(nyin,nxin,helpin,helpc,FFTW_ESTIMATE);
  plan2=fftwf_plan_dft_c2r_2d(nyin,nxin,helpc,helpout,FFTW_ESTIMATE);

  pars->helpin=helpin;
  pars->helpout=helpout;
  pars->fkernel=fkernel;
  pars->helpc=helpc;
  pars->plan1=plan1;
  pars->plan2=plan2;
  pars->nxin=nxin;
  pars->nyin=nyin;
  pars->nxinp=nxinp;
  pars->method=fresize_2d_fft;

/* First transform kernel */

/* Zero entire array */
  for (j=0;j<nyin;j++) {
    for (i=0;i<nxin;i++) {
      helpin[j*nxin+i]=0;
    }
  }

/*  helpin[0]=1.0f; */

/* Copy kernel to new array */
  for (j=-hwidth;j<=hwidth;j++) {
    j1=(j+nyin) % nyin;
    for (i=-hwidth;i<=hwidth;i++) {
      i1=(i+nxin) % nxin;
      helpin[j1*nxin+i1]=pars->ker[(j+hwidth)*fwidth+(i+hwidth)];
    }
  }

/* Transform kernel */
  fftwf_execute(plan1);

/* Copy transform and scale */
  c=1.0f/nxin/nyin;

  for (j=0;j<nyin;j++) { 
    for (i=0;i<nxinp;i++) {
      fkernel[j*nxinp+i]=c*helpc[j*nxinp+i];
    }
  }

  return 0;
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
  int i /*,j*/;

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

int init_fresize_boxcar_fft(
  struct fresize_struct *pars,
  int hwidth, // Half width of boxcar. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  int nxin, // Array size
  int nyin // Array size
)
{
  int status;

  status=init_fresize_boxcar(pars,hwidth,nsub);
  if (status != 0) return status;

  status=make_fft1d(pars,nxin,nyin);
  return status;
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
  int i /*,j*/;
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
    y=(float)exp(-x*x/2); /* There is expf (a float version), but not sure what impact this has. */
    sum=sum+y;
  }
  for (i=0;i<fwidth;i++) {
    x=(i-hwidth)/sigma;
    y=(float)exp(-x*x/2);
    pars->kerx[i]=y/sum;
    pars->kery[i]=y/sum;
  }

  return 0;
}

int init_fresize_gaussian_fft( // Simple square truncated Gaussian. FFT version.
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2)
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  int nxin, // Array size
  int nyin // Array size
)
{
  int status;

  status=init_fresize_gaussian(pars,sigma,hwidth,nsub);
  if (status != 0) return status;

  status=make_fft1d(pars,nxin,nyin);
  return status;
}

int init_fresize_sinc( // Sinc filter
  struct fresize_struct *pars,
  float wsinc, /* Shape is sinc(d/wsinc)*ap(d)
                  wsinc is the amount by which the Nyquist is reduced. 
                  May want wsinc=nsub. */
  int hwidth, // Half width of kernel. Full is 2*hwidth+1.
  int iap, /* Apodization method. Always ap=0 for d>nap*wsinc.
              iap=0 means no apodization ap=1 
              iap=1 uses parabola ap=1-(d/(nap*wsinc))^2
              iap=2 uses sinc ap=sinc(d/(nap*wsinc)) */
  int nap, /* Sinc apodization width in units of wsinc. 
              Normally hwidth=nap*wsinc,
              but hwidth=nap*wsinc-1 works for integer */
  int nsub // Distance between sampled points
)
{
  const int malign=32;
  int fwidth;
  int i /*,j*/;
  float sum,x,y;

  pars->method=fresize_1d;
  pars->nsub=nsub;
  pars->hwidth=hwidth;
  fwidth=2*hwidth+1;
  pars->kerx=(float *)(MKL_malloc(fwidth*sizeof(float),malign));
  pars->kery=(float *)(MKL_malloc(fwidth*sizeof(float),malign));
  sum=0.0f;
  for (i=0;i<fwidth;i++) {
    x=(float)(i-hwidth);
    if (abs(x) > (nap*wsinc)) {
      y=0.0f;
    }
    else {
       y=(float)sinc(x/wsinc);
    }
    if (iap == 1) y=y*(1-(x/(nap*wsinc))*(x/(nap*wsinc)));
    if (iap == 2) y=(float)(y*sinc(x/(nap*wsinc)));
    pars->kerx[i]=y;
    sum=sum+y;
  }
  for (i=0;i<fwidth;i++) {
    pars->kerx[i]=pars->kerx[i]/sum;
    pars->kery[i]=pars->kerx[i];
  }

  return 0;
}

int init_fresize_sinc_fft( // Sinc filter. FFT version.
  struct fresize_struct *pars,
  float wsinc, /* Shape is sinc(d/wsinc)*ap(d)
                  wsinc is the amount by which the Nyquist is reduced.
                  May want wsinc=nsub. */
  int hwidth, // Half width of kernel. Full is 2*hwidth+1.
  int iap, /* Apodization method. Always ap=0 for d>nap*wsinc.
              iap=0 means no apodization ap=1
              iap=1 uses parabola ap=1-(d/(nap*wsinc))^2
              iap=2 uses sinc ap=sinc(d/(nap*wsinc))
              all other cases give ap=1 (not guaranteed) */
  int nap, /* Sinc apodization width in units of wsinc.
              Normally hwidth=nap*wsinc,
              but hwidth=nap*wsinc-1 works for integer */
  int nsub, // Distance between sampled points
  int nxin, // Array size
  int nyin // Array size
)
{
  int status;

  status=init_fresize_sinc(pars,wsinc,hwidth,iap,nap,nsub);
  if (status != 0) return status;

  status=make_fft1d(pars,nxin,nyin);
  return status;
}

int init_fresize_gaussian2( // Circularly truncated Gaussian
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2)
  float rmax, // Truncation radius. Probably rmax<=hwidth.
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub // Distance between sampled points
)
{
  const int malign=32;
  int fwidth;
  int i,j;
  float sum,xi,xj,y,r2,rmax2;

  pars->method=fresize_2d;
  pars->nsub=nsub;
  pars->hwidth=hwidth;
  fwidth=2*hwidth+1;
  pars->ker=(float *)(MKL_malloc(fwidth*fwidth*sizeof(float),malign));
  rmax2=rmax*rmax/sigma/sigma;

  sum=0.0f;
  for (j=0;j<fwidth;j++) {
    xj=(j-hwidth)/sigma;
    for (i=0;i<fwidth;i++) {
      xi=(i-hwidth)/sigma;
      r2=xi*xi+xj*xj;
      if (r2 <= rmax2 ) {
         y=(float)exp(-r2/2);
      }
      else {
        y=0.0f;
      }
      pars->ker[fwidth*j+i]=y;
      sum=sum+y;
    }
  }
// Normalize kernel
  for (j=0;j<fwidth;j++) {
    for (i=0;i<fwidth;i++) {
      pars->ker[fwidth*j+i]=pars->ker[fwidth*j+i]/sum;
    }
  }

  return 0;
}

int init_fresize_gaussian2_fft( // Circularly truncated Gaussian. FFT version.
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2) 
  float rmax, // Truncation radius. Probably rmax<=hwidth.
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  int nxin, // Input size
  int nyin // Input size
)
{
  int status;

  status=init_fresize_gaussian2(pars,sigma,rmax,hwidth,nsub);
  if (status != 0) return status;

  status=make_fft2d(pars,nxin,nyin);
  return status;
}

int init_fresize_airy( // 2D Airy filter
  struct fresize_struct *pars,
  float cdown, /* cdown is the amount by which the Nyquist is reduced. */
  int hwidth, /* Half width of kernel. Full is 2*hwidth+1.
                 Set to <0 to make routine set appropriate value */
  int iap, /* Apodization method. Always ap=0 for d>Z_nap, where Z_nap os
              the position of the nap'th zero.
              iap=0 means no apodization ap=1
              iap=1 uses parabola ap=1-(d/Z_nap)^2
              iap=2 uses sinc ap=sinc(d/(Z_nap))
              iap=3 uses Airy with first zero at Z_nap
              all other cases give ap=1 (not guaranteed) */
  int nap, /* Apodizes to nap'th zero */
  int nsub // Distance between sampled points
)
{
  const int malign=32;
  const float beszeros[20]={0.0000000,3.8317060,7.0155867,10.173468,13.323692,16.470630,19.615859,22.760084,25.903672,29.046829,32.189680,35.332308,38.474766,41.617094,44.759319,47.901461,51.043535,54.185554,57.327525,60.469458}; // first 20 zeros in Bessel function
  float cb;
  int fwidth;
  int i,j;
  float sum,xi,xj;
  float r2,r,rc,airy,ra;
  float nzero;

  cb=(float)(cdown/M_PI); // Amount to scale Bessel function to get zeros right.
  nzero=beszeros[nap]*cb; // Position of nap'th zero
  if (hwidth < 0) hwidth=(float)floor(nzero);

  pars->method=fresize_2d;
  pars->nsub=nsub;
  pars->hwidth=hwidth;
  fwidth=2*hwidth+1;
  pars->ker=(float *)(MKL_malloc(fwidth*fwidth*sizeof(float),malign));

  sum=0.0f;
  for (j=0;j<fwidth;j++) {
     xj=(float)(j-hwidth);
    for (i=0;i<fwidth;i++) {
       xi=(float)(i-hwidth);
      r2=xi*xi+xj*xj;
      r=(float)maxval(sqrt(r2),1e-20);
      rc=r/cb;
      if (r <= nzero ) {
        airy=j1f(rc)/rc;
      }
      else {
        airy=0.0f;
      }
      ra=r/nzero;
      if (iap == 1) airy=airy*(1-ra*ra);
      if (iap == 2) airy=(float)(airy*sinc(ra));
      if (iap == 3) {
        ra=rc*beszeros[1]/beszeros[nap];
        airy=airy*j1f(ra)/ra;
      }
      pars->ker[fwidth*j+i]=airy;
      sum=sum+airy;
    }
  }
// Normalize kernel
  for (j=0;j<fwidth;j++) {
    for (i=0;i<fwidth;i++) {
      pars->ker[fwidth*j+i]=pars->ker[fwidth*j+i]/sum;
    }
  }


  return 0;
}

int init_fresize_airy_fft( // 2D Airy filter. FFT version.
  struct fresize_struct *pars,
  float cdown, /* cdown is the amount by which the Nyquist is reduced. */
  int hwidth, /* Half width of kernel. Full is 2*hwidth+1.
                 Set to <0 to make routine set appropriate value */
  int iap, /* Apodization method. Always ap=0 for d>Z_nap, where Z_nap os
              the position of the nap'th zero.
              iap=0 means no apodization ap=1
              iap=1 uses parabola ap=1-(d/Z_nap)^2
              iap=2 uses sinc ap=sinc(d/(Z_nap))
              iap=3 uses Airy with first zero at Z_nap
              all other cases give ap=1 (not guaranteed) */
  int nap, /* Apodizes to nap'th zero */
  int nsub, // Distance between sampled points
  int nxin, // Input size
  int nyin // Input size
)
{
  int status;

  status=init_fresize_airy(pars,cdown,hwidth,iap,nap,nsub);
  if (status != 0) return status;

  status=make_fft2d(pars,nxin,nyin);
  return status;
}

int free_fresize(
  struct fresize_struct *pars
)
{

  if (pars->method==fresize_1d) {
    MKL_free (pars->kerx);
    MKL_free (pars->kery);
  }

  if (pars->method==fresize_2d) {
    MKL_free (pars->ker);
  }

  if (pars->method==fresize_1d_fft) {
    MKL_free (pars->kerx);
    MKL_free (pars->kery);
    fftwf_free(pars->helpin);
    fftwf_free(pars->fkernel);
    fftwf_free(pars->fkernely);
    fftwf_free(pars->helpc);
    fftwf_free(pars->helpout);
    fftwf_destroy_plan(pars->plan1);
    fftwf_destroy_plan(pars->plan2);
    fftwf_destroy_plan(pars->plan1y);
    fftwf_destroy_plan(pars->plan2y);
  }

  if (pars->method==fresize_2d_fft) {
    MKL_free (pars->ker);
    fftwf_free(pars->helpin);
    fftwf_free(pars->fkernel);
    fftwf_free(pars->helpc);
    fftwf_free(pars->helpout);
    fftwf_destroy_plan(pars->plan1);
    fftwf_destroy_plan(pars->plan2);
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

#pragma omp parallel default(none) \
private(i,j) \
shared(image_in,image_out,imp,fillval) \
shared(imin,imax,jmin,jmax) \
shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
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

#pragma omp parallel default(none) \
private(i,j,i1,j1,impi,sum) \
shared(image_in,image_out,imp,fillval) \
shared(imin,imax,jmin,jmax) \
shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
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
  /*char transpose[] = "t"; */
  char normal[] = "n";
  const int ione = 1;
  const float fone = 1.0f;
  const float fzero = 0.0f;
  int i,j,i1 /*,j1*/;
  int imin,imax,jmin,jmax;
  float /**inp*,*/ *inpi,*inpj;
  float sum;
  int nsub,hwidth;
  float *kerx,*kery,*work;
  int xoffl,xoffh,yoffl,yoffh;
  int nwork;
  /*double t1,t2,t3;*/
  int n1,n2 /*,nchunk*/;

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

  /*nchunk=64;*/
  n1=(imax-imin+1); // Size of matrix for sgemv
  n2=2*hwidth+1; // Size of matrix for sgemv
  /*t1=dsecnd();*/
 
#pragma omp parallel default(none) \
private(i,j,i1,/*j1,*/inpi,inpj,sum)          \
shared(image_in,image_out,fillval) \
shared(imin,imax,jmin,jmax,hwidth,kerx,kery,work,xoffl,yoffl) \
shared(normal,/*transpose,*/n1,n2/*,nchunk*/)               \
shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
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

/*
// Brute force code
#pragma omp for
    for (j=jmin; j<=jmax; j++) {
      inpj=work+(yoffl+j*nsub)*nleadout;
      for (i=imin; i<=imax; i++) {
        sum=0.0f;
        for (j1=0; j1<=2*hwidth; j1++) {
          sum=sum+inpj[j1*nleadout+i]*kery[j1];
        }
        image_out[j*nleadout+i]=sum;
      }
    }
*/

// Code using sgemv
#pragma omp for
    for (j=jmin; j<=jmax; j++) {
      inpj=work+(yoffl+j*nsub)*nleadout;
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

int f1d_fft(
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
  int i,j,i1 /*,j1*/;
  int imin,imax,jmin,jmax;
  float /**inp,*inpi,*/ *inpj;
  /*float sum;*/
  int nsub,hwidth;
  float /**kerx,*kery,*/ *work;
  int xoffl,xoffh,yoffl,yoffh;
  int nwork;
  double t1, /*t2,t3,*/ t4;
  int /*n1,n2,*/ nchunk;
  fftwf_complex *helpc,*fkernel,*fkernely;
  float *helpin,*helpout;
  int nxinp,nyinp;
  float *iwork,*owork;

  if (pars->nxin != nxin) return 1;
  if (pars->nyin != nyin) return 1;

  nsub=pars->nsub;
  hwidth=pars->hwidth;
  /*kerx=pars->kerx;*/
  /*kery=pars->kery;*/
  helpc=pars->helpc;
  fkernel=pars->fkernel;
  fkernely=pars->fkernely;
  helpin=pars->helpin;
  helpout=pars->helpout;
  nxinp=pars->nxinp;
  nyinp=pars->nyinp;

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
  iwork=(float *)(MKL_malloc(nchunk*nyin*sizeof(float),malign));
  owork=(float *)(MKL_malloc(nchunk*nyout*sizeof(float),malign));

// No OpenMP yet, but have left statements in for good measure.
//#pragma omp parallel default(none) \
//private(i,j,i1,j1,inpi,inpj,sum) \
//shared(image_in,image_out,fillval) \
//shared(imin,imax,jmin,jmax,hwidth,kerx,kery,work,xoffl,yoffl) \
//shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
  { // Needed to define parallel region

  t1=dsecnd();
// First convolve in x
//#pragma omp for
    for (j=0;j<nyin;j++) {
      inpj=image_in+j*nleadin;
      for (i=0;i<nxin;i++) helpin[i]=inpj[i];
      fftwf_execute(pars->plan1);
      for (i=0; i<nxinp; i++) helpc[i]=fkernel[i]*helpc[i];
      fftwf_execute(pars->plan2);
      for (i=imin; i<=imax; i++) work[j*nleadout+i]=helpout[i*nsub+xoff];
    }
  t1=dsecnd()-t1;

// Then convolve in y

/*
// Brute force algorithm
  t2=dsecnd();
//#pragma omp for
    for (i=imin; i<=imax; i++) {
      for (j=0;j<nyin;j++) helpin[j]=work[j*nleadout+i];
      fftwf_execute(pars->plan1y);
      for (j=0; j<nyinp; j++) helpc[j]=fkernely[j]*helpc[j];
      fftwf_execute(pars->plan2y);
      for (j=jmin; j<=jmax; j++) image_out[j*nleadout+i]=helpout[j*nsub+yoff];
    }
  t2=dsecnd()-t2;

// Block the output array.
  t3=dsecnd();
//#pragma omp for
  for (i1=imin;i1<=imax;i1=i1+nchunk) {
    for (i=i1; i<=minval(i1+nchunk-1,imax); i++) {
      for (j=0;j<nyin;j++) helpin[j]=work[j*nleadout+i];
      fftwf_execute(pars->plan1y);
      for (j=0; j<nyinp; j++) helpc[j]=fkernely[j]*helpc[j];
      fftwf_execute(pars->plan2y);
      for (j=jmin; j<=jmax; j++) owork[j+(i-i1)*nyout]=helpout[j*nsub+yoff];
    }
    for (j=jmin; j<=jmax; j++) {
      for (i=i1; i<=minval(i1+nchunk-1,imax); i++) {
        image_out[j*nleadout+i]=owork[j+(i-i1)*nyout];
      }
    }
  }
  t3=dsecnd()-t3;
*/

// Block both the input output arrays.
  t4=dsecnd();
//#pragma omp for
  for (i1=imin;i1<=imax;i1=i1+nchunk) {
    for (j=0;j<nyin;j++) {
      for (i=i1; i<=minval(i1+nchunk-1,imax); i++) {
        iwork[(i-i1)*nyin+j]=work[j*nleadout+i];
      }
    }
    for (i=i1; i<=minval(i1+nchunk-1,imax); i++) {
      for (j=0;j<nyin;j++) helpin[j]=iwork[(i-i1)*nyin+j];
      fftwf_execute(pars->plan1y);
      for (j=0; j<nyinp; j++) helpc[j]=fkernely[j]*helpc[j];
      fftwf_execute(pars->plan2y);
      for (j=jmin; j<=jmax; j++) owork[j+(i-i1)*nyout]=helpout[j*nsub+yoff];
    }
    for (j=jmin; j<=jmax; j++) {
      for (i=i1; i<=minval(i1+nchunk-1,imax); i++) {
        image_out[j*nleadout+i]=owork[j+(i-i1)*nyout];
      }
    }
  }
  t4=dsecnd()-t4;

// Fill below valid region
//#pragma omp for
    for (j=0; j<jmin; j++) {
      for (i=0; i<nxout; i++) {
        image_out[j*nleadout+i]=fillval;
      } // i=
    } //j=

// Valid heights
//#pragma omp for
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
//#pragma omp for
    for (j=jmax+1; j<nyout; j++) {
      for (i=0; i<nxout; i++) {
        image_out[j*nleadout+i]=fillval;
      } // i=
    } //j=

  } // End of parallel region

  MKL_free(iwork);
  MKL_free(owork);

  return 0;
}

int f2d(
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
  /*const int malign=32;*/
  /*char transpose[] = "t";*/
  /*char normal[] = "n";*/
  /*const int ione = 1;*/
  /*const float fone = 1.0f;*/
  /*const float fzero = 0.0f;*/
  int i,j,i1,j1;
  int imin,imax,jmin,jmax;
  float *inpi,*kerp;
  float sum;
  int nsub,hwidth,fwidth;
  float *ker;
  int xoffl,xoffh,yoffl,yoffh;
  /*double t1,t2,t3;*/

  nsub=pars->nsub;
  hwidth=pars->hwidth;
  ker=pars->ker;
  fwidth=2*hwidth+1;

// Find first and last indices for which resulting input point is within image
  xoffl=xoff-hwidth;
  xoffh=xoff+hwidth;
  if (xoffl >= 0) imin=0; else imin=((-xoffl+nsub-1)/nsub);
  if (xoffh <= 0) imax=nxout-1; else imax=minval(((nxin-xoffh-1)/nsub),nxout-1);

  yoffl=yoff-hwidth;
  yoffh=yoff+hwidth;
  if (yoffl >= 0) jmin=0; else jmin=((-yoffl+nsub-1)/nsub);
  if (yoffh <= 0) jmax=nyout-1; else jmax=minval(((nyin-yoffh-1)/nsub),nyout-1);


  /*t1=dsecnd();*/
 
#pragma omp parallel default(none) \
private(i,j,i1,j1,inpi,kerp,sum) \
shared(image_in,image_out,fillval) \
shared(imin,imax,jmin,jmax,hwidth,fwidth,ker,xoffl,yoffl) \
shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
  { // Needed to define parallel region

#pragma omp for
    for (j=jmin; j<=jmax; j++) {
      for (i=imin; i<=imax; i++) {
        inpi=image_in+(j*nsub+yoffl)*nleadin+i*nsub+xoffl; // Note offset to match kernels.
        sum=0.0f;
        for (j1=0; j1<=2*hwidth; j1++) {
          kerp=ker+j1*fwidth;
          for (i1=0; i1<=2*hwidth; i1++) {
            sum=sum+kerp[i1]*inpi[i1];
          }
          inpi=inpi+nleadin;
        }
        image_out[j*nleadout+i]=sum;
      }
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

  return 0;
}

// this version uses matrix vector multiplys. Not much faster.
int f2d_mat(
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
  /*char normal[] = "n";*/
  const int ione = 1;
  const float fone = 1.0f;
  /*const float fzero = 0.0f;*/
  int i,j,/*i1,*/j1;
  int imin,imax,jmin,jmax;
  float *inpi,*kerp;
  /*float sum;*/
  int nsub,hwidth,fwidth;
  float *ker;
  int xoffl,xoffh,yoffl,yoffh;
  /*double t1,t2,t3;*/
  int /*n1,*/n2,nchunk,jc,nc;
  float *work;

  nsub=pars->nsub;
  hwidth=pars->hwidth;
  ker=pars->ker;
  fwidth=2*hwidth+1;

// Find first and last indices for which resulting input point is within image
  xoffl=xoff-hwidth;
  xoffh=xoff+hwidth;
  if (xoffl >= 0) imin=0; else imin=((-xoffl+nsub-1)/nsub);
  if (xoffh <= 0) imax=nxout-1; else imax=minval(((nxin-xoffh-1)/nsub),nxout-1);

  yoffl=yoff-hwidth;
  yoffh=yoff+hwidth;
  if (yoffl >= 0) jmin=0; else jmin=((-yoffl+nsub-1)/nsub);
  if (yoffh <= 0) jmax=nyout-1; else jmax=minval(((nyin-yoffh-1)/nsub),nyout-1);


  /*t1=dsecnd();*/
 
/*n1=jmax-jmin+1;*/
n2=nleadin*nsub;
nchunk=64;

#pragma omp parallel default(none) \
private(i,j,/*i1,*/j1,inpi,kerp,/*sum,*/work,nc)      \
shared(/*normal,*/transpose,/*n1,*/n2,nchunk)       \
shared(image_in,image_out,fillval) \
shared(imin,imax,jmin,jmax,hwidth,fwidth,ker,xoffl,yoffl) \
shared(nxin,nyin,nleadin,nxout,nyout,nleadout,nsub)
  { // Needed to define parallel region

work=(float *)(MKL_malloc(nyout*sizeof(float),malign));
#pragma omp for
    for (j=jmin; j<=jmax; j++) {
      for (i=imin; i<=imax; i++) {
        image_out[j*nleadout+i]=0.0f;
      }
    }

#pragma omp for
for (jc=jmin;jc<=jmax;jc=jc+nchunk) {
nc=minval(jc+nchunk-1,jmax)-jc+1;
    for (j1=0; j1<=2*hwidth; j1++) {
      kerp=ker+j1*fwidth;
      for (i=imin; i<=imax; i++) {
        inpi=image_in+(jc*nsub+yoffl+j1)*nleadin+i*nsub+xoffl;
        sgemv(transpose,&fwidth,&nc,&fone,inpi,&n2,kerp,&ione,&fone,image_out+jc*nleadout+i,&nleadout);
      }
    }
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

    MKL_free(work);
  } // End of parallel region

  return 0;
}

int f2d_fft(
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
  fftwf_complex *helpc,*fkernel;
  float *helpin,*helpout;
  fftwf_plan plan1,plan2;
  int nxinp; // Complex array size
  int hwidth/*,fwidth*/;
  int i,j,/*i1,*/j1;
  /*float c;*/
  int imin,imax,jmin,jmax;
  int xoffl,xoffh,yoffl,yoffh;
  int nsub;

  if (pars->nxin != nxin) return 1;
  if (pars->nyin != nyin) return 1;

  nsub=pars->nsub;
  hwidth=pars->hwidth;
  /*fwidth=2*hwidth+1;*/
  nxinp=(nxin/2+1);

  helpin=pars->helpin;
  fkernel=pars->fkernel;
  helpc=pars->helpc;
  helpout=pars->helpout;
  plan1=pars->plan1;
  plan2=pars->plan2;

/* Now copy and transform input image */

  for (j=0;j<nyin;j++) {
    for (i=0;i<nxin;i++) {
      helpin[j*nxin+i]=image_in[j*nleadin+i];
    }
  }

  fftwf_execute(plan1);

/* Multiply by kernel transform */

  for (j=0;j<nyin;j++) {
    for (i=0;i<nxinp;i++) {
      helpc[j*nxinp+i]=fkernel[j*nxinp+i]*helpc[j*nxinp+i];
    }
  }

  fftwf_execute(plan2);

// Copy data to output array

// Find first and last indices for which resulting input point is within image
  xoffl=xoff-hwidth;
  xoffh=xoff+hwidth;
  if (xoffl >= 0) imin=0; else imin=((-xoffl+nsub-1)/nsub);
  if (xoffh <= 0) imax=nxout-1; else imax=minval(((nxin-xoffh-1)/nsub),nxout-1);

  yoffl=yoff-hwidth;
  yoffh=yoff+hwidth;
  if (yoffl >= 0) jmin=0; else jmin=((-yoffl+nsub-1)/nsub);
  if (yoffh <= 0) jmax=nyout-1; else jmax=minval(((nyin-yoffh-1)/nsub),nyout-1);

  for (j=jmin; j<=jmax; j++) {
    j1=(j*nsub+yoff)*nxin+xoff;
    for (i=imin; i<=imax; i++) {
      image_out[j*nleadout+i]=helpout[j1+i*nsub];
    }
  }

// Fill below valid region
  for (j=0; j<jmin; j++) {
    for (i=0; i<nxout; i++) {
      image_out[j*nleadout+i]=fillval;
    } // i=
  } //j=

// Valid heights
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
  for (j=jmax+1; j<nyout; j++) {
    for (i=0; i<nxout; i++) {
      image_out[j*nleadout+i]=fillval;
    } // i=
  } //j=

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

  case fresize_1d_fft:
    status=f1d_fft(pars,image_in,image_out,nxin,nyin,nleadin,nx,ny,nlead,xoff,yoff,fillval);
    break;

  case fresize_2d:
    status=f2d_mat(pars,image_in,image_out,nxin,nyin,nleadin,nx,ny,nlead,xoff,yoff,fillval);
    break;

  case fresize_2d_fft:
    status=f2d_fft(pars,image_in,image_out,nxin,nyin,nleadin,nx,ny,nlead,xoff,yoff,fillval);
    break;

  default:
    status=1;
    break;
  }

return status;

}

int init_fresize_user(
  struct fresize_struct *pars,
  int hwidth, // Half width of kernel. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  float *user_ker // User specified kernel to convolve with.
                  // Must be of size (2*hwidth+1) x (2*hwidth+1).
                  // Kernel need not be and will not be normalized.
)
{
  const int malign=32;
  int fwidth;
  int i,j;

  pars->method=fresize_2d;
  pars->nsub=nsub;
  pars->hwidth=hwidth;
  fwidth=2*hwidth+1;
  pars->ker=(float *)(MKL_malloc(fwidth*fwidth*sizeof(float),malign));
  memcpy(pars->ker,user_ker,sizeof(float)*fwidth*fwidth);

  return 0;
}

