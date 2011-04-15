struct fresize_struct {

  int method;

  int nsub;

  int hwidth;

  float *ker,*kerx,*kery;

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



int init_fresize_gaussian(

  struct fresize_struct *pars,

  float sigma, // Shape is exp(-(d/sigma)^2/2)

  int hwidth, // Half (truncation) width of kernel. Full is 2*hwidth+1.

  int nsub // Distance between sampled points

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


