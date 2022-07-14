#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>

#include "tinterpolate.h"
#include "finterpolate.h"
#include "gapfill.h"



struct initial{
  float *dist_coef;
  int order_dist_coef;
  float *diffrot_coef;
  int order2_rot_coef;
  int order_int;
  int nconst;
  struct fill_struct fills;
};

struct keyword{
  float rsun;
  float xx0;
  float yy0;
  float dist;
  float b0;
  float p0;
  double time;
};



int initialize_interpol(struct initial *const_param);

int do_interpolate(float **images, float **errors, float *image_out, struct keyword *key, struct keyword *key_out, struct initial *const_param, int nsample, int nx, int ny);

int do_gapfill(float *image, unsigned char *masks, struct initial *const_param, float *ierror, int nx, int ny);

void free_interpol(struct initial *const_param);






//#define minval(x,y) (((x) < (y)) ? (x) : (y))
//#define maxval(x,y) (((x) > (y)) ? (x) : (y))

