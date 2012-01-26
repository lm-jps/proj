
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `smoothsphere' a.k.a. `ssp' as a C library.
 *
 * Made by intermediate binary `smoothsphere.out' on:
 * 	Thu Jul  7 23:02:43 2011
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * smoothsphere: smooth a projected sphere with a kernel
 * 
 *  [y,s,p]=smoothsphere(x,geom,k,kparam,kwt,bws)
 *  * Given a solar image of location given by `geom', smooth it
 *  using the radially-symmetric kernel `k'.
 *  * This version is threaded, using 8 threads by default.
 *  Set MXT_NUM_THREADS in the environment to lower this number.
 *  * Off-disk values are indicated by NaN in the output.
 *  * The kernel used is dependent only on the weighted distance
 *  between two positions, say P1 and P2, in three-dimensional
 *  coordinates, normalized to live on the unit sphere:
 *    d = P1 - P2
 *    dist = d' W d  (>= 0)
 *  for a diagonal weight matrix W.
 *  The kernel values in `k' are specified for values of `dist'
 *  in a linear range from 0 to kparam(1), typically 0.015.
 *  Given a certain value of `dist' found between an image pixel
 *  and the kernel center, the associated weight is interpolated
 *  linearly using the table.  If dot > kparam(1), a zero weight
 *  is used.
 *  * We have typically let k be a gaussian kernel of width
 *  (standard deviation) = 0.0325, i.e.,
 *    k(dist) = C * exp[ -0.5 * dist / sqr(0.0325) ],
 *  since dist is a squared quantity already, and C is a constant
 *  chosen to make k have unit norm as a spherical kernel.
 *  * As a shortcut, if k is a scalar, the kernel is taken to
 *  be the above function, but with width 0.0325 multiplied by
 *  k(1).  In this case, the LUT consists of 256 points evenly
 *  spaced over [0,kparam(1)].
 *  * The diagonal portion of the distance weight matrix W may be
 *  supplied as a triple kwt.  This weights distances between
 *  P1 and P2 in the x, y, and z directions respectively.   The
 *  P-angle (rotation in the x-y plane) is taken into account in
 *  this weighting, so that y is cross-track, and x and z are
 *  always along-track.
 *  * Typically the kernel falls to zero rather quickly.  As a
 *  computational shortcut, it is assumed that the kernel extends
 *  only kparam(2) pixels on each side of its center.  For W = I,
 *  this implies that an "on" pixel in x extends to influence at
 *  most a "swath" 2*kparam(2)+1 pixels on a side.  Both the P-angle
 *  and the W matrix are taken into account in finding the swath.
 *  For instance, a W_y > 1 will cause the y portion of the swath
 *  to shrink, because distance decreases faster.  The swath is
 *  always parallel to the (i,j) image axes, so P-angles not a
 *  multiple of 90 degrees will cause the swath to be enlarged
 *  to contain a rotated rectangle.
 *  * To decrease computational load, nearby pixels may be grouped
 *  into meta-pixels called blocks.  This is especially trouble-free
 *  away from the limb, where local geometry is nearly planar, so
 *  the blocking is given as a table that depends on z (in [0,1]).
 *  A row in the table of (z,bw) means: above the value z, use
 *  a blocking of bw (>1).  For z=0 to bws(1,1), no blocking is used.
 *  * For example, the small-image default means to use single pixels
 *  below 0.4, 2x2 blocks above 0.4, and 4x4 above 0.6.  For typical
 *  MDI images, this is 16%, 20%, and 64% of on-disk pixels, respectively.
 *  Using block sizes that do not pack together is legal, but causes
 *  inefficiency due to poor fits between abutting blocks; diagnose
 *  bad packing with the p output.
 *  * Supply bws=[] to treat all pixels singly (turn off blocking).
 *  * The optional output s is the z-coordinate (in [0,R]).
 *  The optional p output is the block (patch) number of each
 *  pixel.
 * 
 *  Inputs:
 *    real x[m,n];
 *    real geom[5];  -- [x0 y0 rsun b0 p0]
 *    real k[p] or k[1];
 *    real kparam[2];  -- [top window]
 *    opt real kwt = [1 1 1];
 *    opt real bws[bwnum,2] = [0.4 2;0.6 4];         -- m <  2048
 *                          = [0.3 2;0.5 4; 0.7 8];  -- m >= 2048
 * 
 *  Outputs:
 *    real y[m,n];
 *    opt real s[m,n];
 *    opt int  p[m,n];
 * 
 *  See also:
 * 
 *  implemented as a mex file
 * 
 * 
 */

#ifndef _mexfn_smoothsphere_h_
#define _mexfn_smoothsphere_h_

// function entry point
mexfn_lib_t main_smoothsphere;

// argument counts
#define MXT_ssp_NARGIN_MIN 	4
#define MXT_ssp_NARGIN_MAX 	6
#define MXT_ssp_NARGOUT_MIN	0
#define MXT_ssp_NARGOUT_MAX	3

// input argument numbers
#define MXT_ssp_ARG_x	0
#define MXT_ssp_ARG_geom	1
#define MXT_ssp_ARG_k	2
#define MXT_ssp_ARG_kparam	3
#define MXT_ssp_ARG_kwt	4
#define MXT_ssp_ARG_bws	5

// output argument numbers
#define MXT_ssp_ARG_y	0
#define MXT_ssp_ARG_s	1
#define MXT_ssp_ARG_p	2


#endif // _mexfn_smoothsphere_h_

// (file ends)
