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

/* not used anywhere #define nanval -2147483648 */

/* Moved all these definition to limb_fit_functions.c - should not define
 * any global variables in a header file. */
extern const double rad_corr_fac;
extern const double foc_corr;
extern const double high;
extern const double low;

extern const double limit_var;
extern const double limit_cc;

extern const int lim;

extern const int parsize;
//gapfill constants

extern const int gapfill_order;
extern const int gapfill_method;
extern const float gapfill_regular;
extern const int gapfill_order2;

extern const int min_imcnf;
extern const int max_imcnf;

extern char *X0_MP_key;
extern char *Y0_MP_key;
extern char *RSUN_OBS_key;
extern char *IMSCL_MP_key;

extern char *HCAMID_key;
extern char *HCFTID_key;
extern char *MISSVAL_key;

extern int light_val1;
extern int light_val2;

extern const double percent_good;



//int limb_fit(DRMS_Record_t *record, float *image_in, double *rsun_lf, double *x0_lf, double *y0_lf, int nx, int ny, int method, const char *dpath);
int limb_fit(DRMS_Record_t *record, float *image_in, double *rsun_lf, double *x0_lf, double *y0_lf, int nx, int ny, int method);

#endif
