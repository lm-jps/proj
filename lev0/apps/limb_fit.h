#ifndef __LIMBFIT_H
#define __LIMBFIT_H
struct mempointer
{
  double *xrp;
  double *yrp;
  double *imrphi;
  double *rc;
  double *phic;

  double *avgphi;

  float *im;
  double *imhp;
  double *imhp_full;
  double *imro;
  float *image;
  float *image_full;
  double *parab;
  unsigned char *mask_p;
  float *cnorm;
  float *ierror;
};

#define nanval -2147483648

const int cent_err=20;

  const double high=1.02;
  const double low=0.98;

  const double limit_var=20.0;
  const double limit_cc=20.0;

  const int lim=8;
  const int parsize=2*lim+1;

//gapfill constants

const int gapfill_order=11;
const int gapfill_method=11;
const float gapfill_regular=0.0025;
const int gapfill_order2=gapfill_order/2;

const int min_imcnf=80;
const int max_imcnf=124;

char *X0_MP_key="X0_MP";
char *Y0_MP_key="Y0_MP";
char *RSUN_OBS_key="RSUN_OBS";
char *IMSCL_MP_key="IMSCL_MP";

char *HCAMID_key="HCAMID";
char *HCFTID_key="HCFTID";

int light_val1=2;
int light_val2=3;

const double percent_good=0.75;

int limb_fit(DRMS_Record_t *record, int *image_in, double *rsun_lf, double *x0_lf, double *y0_lf, int nx, int ny, int method);

#endif
