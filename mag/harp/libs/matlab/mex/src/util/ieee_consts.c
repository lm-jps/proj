/*
 * ieee_consts: Generate IEEE special values
 *
 * mxt_getnanf: float NaN
 * mxt_getnand: double NaN
 * mxt_getinfd: double +Infinity
 */

/*LINTLIBRARY*/

#include "ieee_consts.h"  /* ensure consistency */

/* from Harbison&Steele by way of GNU configure ...
   returns 1 for bigendian, 0 for littleendian */

static int is_bigendian(void)
{
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1L;
  return(u.c[sizeof (long) - 1] == 1);
}

/* below, set the inf/nan buffers up only once */

float mxt_getnanf(void)
{
  static int ieee_fnan_inited = 0;
  /* the union ensures the underlying float is aligned properly */
  static union {
    float f;
    unsigned char c[4];
  } nan;
  
  if (!ieee_fnan_inited) {
    int i;
    for (i = 0; i < 4; i++)
      nan.c[i] = 1;
    if (is_bigendian()) {
      nan.c[0] = 0x7F;
      nan.c[1] = 0x80;
    } else {
      nan.c[3] = 0x7F;
      nan.c[2] = 0x80;
    }
    ieee_fnan_inited = 1; /* flag: buffer is set up now */
  }
  return(nan.f);
}

double mxt_getnand(void)
{
  static int ieee_nan_inited = 0;
  /* the union ensures the underlying double is aligned properly */
  static union {
    double d;
    unsigned char c[8];
  } nan;

  if (!ieee_nan_inited) {
    int i;
    for (i = 0; i < 8; i++)
      nan.c[i] = 1;  /* have seen 0 here also */
    if (is_bigendian()) {
      nan.c[0] = 0x7F;
      nan.c[1] = 0xF0; /* with 0 above goes 0xF8 here */
    } else {
      nan.c[7] = 0x7F;
      nan.c[6] = 0xF0; /* ditto */
    }
    ieee_nan_inited = 1; /* flag: buffer is set up now */
  }
  return(nan.d);
}


double mxt_getinfd(void)
{
  static int ieee_inf_inited = 0;
  /* the union ensures the underlying double is aligned properly */
  static union {
    double d;
    unsigned char c[8];
  } inf;
  
  if (!ieee_inf_inited) {
    int i;
    for (i = 0; i < 8; i++)
      inf.c[i] = 0;
    if (is_bigendian()) {
      inf.c[0] = 0x7F;
      inf.c[1] = 0xF0;
    } else {
      inf.c[7] = 0x7F;
      inf.c[6] = 0xF0;
    }
    ieee_inf_inited = 1; /* buffer set up */
  }
  return(inf.d);
}
