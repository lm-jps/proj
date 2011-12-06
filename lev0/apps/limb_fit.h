#ifndef __LIMB_FIT_H
#define __LIMB_FIT_H
struct mempointer
{
  double *xrp;
  double *yrp;
  double *imrphi;
  double *rc;
  double *phic;

  double *avgphi;

  double *imhp;
  double *imro;
  float *image;
  double *parab;
  unsigned char *mask_p;
  float *cnorm;
  float *ierror;
  float *imcp;
};

#define nanval -2147483648

const double rad_corr_fac=1.00130;
const double foc_corr=-1.08490;
  const double high=1.03;
  const double low=0.97;

  const double limit_var=1.8;
  const double limit_cc=15.0;

// work-around for compiler problem below
#define kLimbFit_lim 8
//const int lim=8;
const int lim=kLimbFit_lim;

// this is not legal - can't initialize global vars with variables.
// for now, replace lim with its value to fix the build.
//const int parsize=2*lim+1;
const int parsize=2*kLimbFit_lim+1;
//gapfill constants

// work-around for compiler problem below
#define kLimbFit_gapfill_order 11
//const int gapfill_order=11;
const int gapfill_order=kLimbFit_gapfill_order;
const int gapfill_method=11;
const float gapfill_regular=0.0025;
// this is not legal - can't initialize global vars with variables.
// for now, replace lim with its value to fix the build.
//const int gapfill_order2=gapfill_order/2;
const int gapfill_order2=kLimbFit_gapfill_order/2;

const int min_imcnf=80;
const int max_imcnf=1024;

char *X0_MP_key="X0_MP";
char *Y0_MP_key="Y0_MP";
char *RSUN_OBS_key="RSUN_OBS";
char *IMSCL_MP_key="IMSCL_MP";

char *HCAMID_key="HCAMID";
char *HCFTID_key="HCFTID";
char *MISSVAL_key="MISSVALS";

int light_val1=2;
int light_val2=3;

const double percent_good=0.75;



int limb_fit(DRMS_Record_t *record, float *image_in, double *rsun_lf, double *x0_lf, double *y0_lf, int nx, int ny, int method, const char *dpath);

#endif
