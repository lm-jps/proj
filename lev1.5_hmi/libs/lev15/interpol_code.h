#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>

#include "tinterpolate.h"
#include "finterpolate.h"
#include "gapfill.h"


struct init_files{
  char *dist_file_front;
  char *dist_file_side;
  char *diffrot_coef;
};


struct initial{
  float *dist_coef[2][16];
  int order_dist_coef;
  float *diffrot_coef;
  int order2_rot_coef;
  int order_int;
  int nconst;
  char *code_version;
  struct fill_struct fills;
};

struct keyword{
  int camera;
  int focus;
  float rsun;
  float xx0;
  float yy0;
  float dist;
  float b0;
  float p0;
  double time;
};



int initialize_interpol(struct initial *const_param, struct init_files *infiles, int nx, int ny, const char *dpath);

int do_interpolate(float **images, char **errors, float *image_out, struct keyword *key, struct keyword *key_out, struct initial *const_param, int nsample, int nx, int ny, float average, const char *dpath);

int do_gapfill(float *image, unsigned char *masks, struct initial *const_param, char *ierror, int nx, int ny);

void free_interpol(struct initial *const_param);

char *interpol_version();



//#define minval(x,y) (((x) < (y)) ? (x) : (y))
//#define maxval(x,y) (((x) > (y)) ? (x) : (y))

