#include <complex.h>
#include <fftw3.h>

struct fresize_struct {
  int method;
  int nsub;
  int hwidth;
  float *ker,*kerx,*kery;
  fftwf_complex *helpc,*fkernel,*fkernely;
  float *helpin,*helpout;
  fftwf_plan plan1,plan2,plan1y,plan2y;
  int nxin,nyin,nxinp,nyinp;
};

int init_fresize_sample(
  struct fresize_struct *pars,
  int nsub // Distance between sampled points
);

int init_fresize_bin(
  struct fresize_struct *pars,
  int nsub // Binsize
);

int init_fresize_boxcar(
  struct fresize_struct *pars,
  int hwidth, // Half width of boxcar. Full is 2*hwidth+1.
  int nsub // Distance between sampled points
);

int init_fresize_boxcar_fft(
  struct fresize_struct *pars,
  int hwidth, // Half width of boxcar. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  int nxin, // Array size
  int nyin // Array size
);

int init_fresize_gaussian( // Simple square truncated Gaussian
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2)
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub // Distance between sampled points
);

int init_fresize_gaussian_fft( // Simple square truncated Gaussian. FFT version.
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2)
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  int nxin, // Array size
  int nyin // Array size
);

int init_fresize_sinc( // Sinc filter
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
  int nsub // Distance between sampled points
);

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
);

int init_fresize_gaussian2( // Circularly truncated Gaussian
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2)
  float rmax, // Truncation radius. Probably rmax<=hwidth.
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub // Distance between sampled points
);

int init_fresize_gaussian2_fft( // Circularly truncated Gaussian. FFT version.
  struct fresize_struct *pars,
  float sigma, // Shape is exp(-(d/sigma)^2/2)
  float rmax, // Truncation radius. Probably rmax<=hwidth.
  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  int nxin, // Input size
  int nyin // Input size
);

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
);

int init_fresize_airy_fft( // 2D Airy filter. FFT version.
  struct fresize_struct *pars,
  float cdown, /* cdown is the amount by which the Nyquist is reduced. */
  int hwidth, /* Half width of kernel. Full is 2*hwidth+1.
                 Set to <0 to make routine set appropriate value */
  int iap, /* Apodization method. Always ap=0 for d>nap*cdown.
              iap=0 means no apodization ap=1
              iap=1 uses parabola ap=1-(d/(nap*cdown))^2
              iap=2 uses sinc ap=sinc(d/(nap*cdown))
              iap=3 uses Airy with first zero at nap'th zero of main
              all other cases give ap=1 (not guaranteed) */
  int nap, /* Apodizes to nap'th zero */
  int nsub, // Distance between sampled points
  int nxin, // Input size
  int nyin // Input size
);

int init_fresize_user(
  struct fresize_struct *pars,
  int hwidth, // Half width of kernel. Full is 2*hwidth+1.
  int nsub, // Distance between sampled points
  float *user_ker // User specified kernel to convolve with.
                  // Must be of size (2*hwidth+1) x (2*hwidth+1).
                  // Kernel need not be and will not be normalized.
);

int free_fresize(
  struct fresize_struct *pars
);

int fresize(
  struct fresize_struct *pars, // Must have been initialized by init_fresize_XXX
  float *image_in,
  float *image_out,
  int nxin, // Size of input image
  int nyin, // Size of input image
  int nleadin, // Leading dimension of input image. nleadin>=nxin
  int nxout, // Size of xin, yin and image_out
  int nyout, // Size of xin, yin and image_out
  int nleadout, // Leading dimension. nlead>=nx
  int xoff, // Offset in x direction
  int yoff, // Offset in y direction
  float fillval // Value to use if outside area
);

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
);

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
);

