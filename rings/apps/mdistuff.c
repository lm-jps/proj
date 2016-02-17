/*
 *  Local function for correction of MDI image distortion; these are copied
 *    directly from the function in ~CM/src/libMDI.d/mdi_corrections.c
 *    Only the names of the functions and defineds have been changed to
 *    avoid possible conflicts.
 *
 */
/*
 *  Assumed error in SOHO position angle w.r.t solar north (deg.)
 */
#include "math.h"
 
#define MDI_IMG_SOHO_PA	(-0.2)

void mtrack_MDI_correct_pa (double *pa) {
  *pa += MDI_IMG_SOHO_PA * M_PI / 180.0;
}

/*
 *  Correct for error in measured center coordinates due to plate tipping;
 *    the constants are as defined above
 */

#define MDI_IMG_ACPA	(1.01e-3)
#define MDI_IMG_ASPA	(-1.49e-3)

void mtrack_MDI_correct_imgctr (double *xc, double *yc, double rsun) {
  double rs2;

  rs2 = rsun * rsun / 512.0;
  *xc -= MDI_IMG_ASPA * rs2;
  *yc -= MDI_IMG_ACPA * rs2;
}

/*
 *  mtrack_MDI_image_stretch
 *
 *  Modify "plate" coordinates to account for known optical distortions
 *    in the MDI instrument
 *  It is assumed that the coordinates *x and *y are in terms of half the
 *    full plate width, i.e. that they are in the range [-1.0, 1.0] and
 *    relative to its center; this of course requires an external correction
 *    to be applied in the cases of extracted rasters and binned data.
 *    For MDI the half-plate-width is 512 * 21 um.  The 2nd-order radial
 *    correction constant is given in Kuhn et al. (Ap.J. 613, 1241; 2004)
 *    as 1.58e-3 where the radial unit is in cm. The constant used here is thus
 *    1.58e-3  * (.0021 * 512)^2
 *
 *  By ignoring the 4th-order term the function can be used for inverse
 *    as well as direct stretching.
 *
 *  The value of the integer direct specifies the direction of the
 *    transformation: +1 for a correction from "perfect" to plate coordinates, 
 *    -1 for transformation from plate to perfect
 *
 *  Bugs:
 *    There is no check that |direct| = 1; it can actually be changed as a
 *	scale factor (useful for testing)
 *
 */

#define MDI_IMG_STRETCH	(1.83e-3)

void mtrack_MDI_image_stretch (double *x, double *y, int n, int direct) {
  double f, r2, s;

  s = direct * MDI_IMG_STRETCH;
  while (n--) {
    r2 = *x * *x + *y * *y;
    f = 1.0 + s * r2;
    *x++ *= f;
    *y++ *= f;
  }
}

/*
 *  mtrack_MDI_image_tip
 *
 *  Correct for ellipticity of image due to plate tipping
 *  The constants are:
 *    TIP = 2.61 deg = 0.04555
 *    EFL = 25.3 (effective focal length in units of plate half-width)
 *    PA = -56 deg
 *    SPA = sin (PA),  CPA = cos (PA)
 *    AEP = TIP / EFL
 *    BEP = TIP^2 / 4 = 5.187e-4
 *    BCPA = BEP * CPA,  BSPA = BEP * SPA
 *
 *  Bugs:
 *    There is no check that |direct| = 1; it can actually be changed as a
 *	scale factor (useful for testing)
 *    
 */

#define MDI_IMG_SPA	(-0.8290)
#define MDI_IMG_CPA	(0.5592)
#define MDI_IMG_AEP	(1.80e-3)
#define MDI_IMG_BCPA	(2.90e-4)
#define MDI_IMG_BSPA	(-4.30e-4)

void mtrack_MDI_image_tip (double *x, double *y, int n, int direct) {
  double x0, y0, s, t;

  while (n--) {
    x0 = *x;
    y0 = *y;
    t = direct * (MDI_IMG_SPA * x0 + MDI_IMG_CPA * y0);
    s = direct * (MDI_IMG_CPA * x0 - MDI_IMG_SPA * y0);
    *x += MDI_IMG_BSPA * t;
    *y += MDI_IMG_BCPA * t;
    *x -= MDI_IMG_BCPA * s;
    *y += MDI_IMG_BSPA * s;
    t *= MDI_IMG_AEP;
    *x++ += t * x0;
    *y++ += t * y0;
  }
}

