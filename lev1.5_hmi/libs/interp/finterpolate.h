#define fint_test 1
#define fint_linear 2
#define fint_cubic_conv 3

struct fint_struct {
  int method;
  int maxorder;
  int nsub;
  float fover;
  float *kersx;
};

int init_finterpolate(
  struct fint_struct *pars,
  int method // See above
);

int free_finterpolate(
  struct fint_struct *pars
);

int finterpolate(
  struct fint_struct *pars, // Must have been initialized by init_finterpolate
  int order, // Ignored in some cases (eg. for fint_linear)
  float *image_in,
  float *xin,
  float *yin,
  float *image_out,
  int nxin, // Size of input image
  int nyin, // Size of input image
  int nleadin, // Leading dimension of input image. nleadin>=nxin
  int nx, // Size of xin, yin and image_out
  int ny, // Size of xin, yin and image_out
  int nlead, // Leading dimension. nlead>=nx
  float maxext, // Max extrapolation distance
  float fillval // Value to use if outside area
);

