
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `smoothsphere' a.k.a. `ssp' as a C library.
 *
 * Made by intermediate binary smoothsphere.out on Mon May  3 17:24:41 2010
 *
 * Source file ../Gen-include.c last modified on Mon May  3 17:04:13 2010
 */

// Original documentation block
/*
 * smoothsphere: smooth a projected sphere with a kernel
 * 
 *  [y,s]=smoothsphere(x,center,p0,k,kparam,kwt,bws)
 *  * Given a solar image of location given by center, and
 *  p-angle p0, smooth it using the radially-symmetric kernel
 *  listed in k.
 * * * Off-disk values are indicated by NaN in the output.
 *  * The kernel used is dependent only on the weighted distance
 *  between two positions, say P1 and P2, in three-dimensional
 *  coordinates, normalized to live on the unit sphere:
 *    d = P1 - P2
 *    dist = d' W d  (>= 0)
 *  for a diagonal weight matrix W.
 *  The kernel values are specified for values of "dist" taken
 *  from a linear range from 0 to kparam(1).  Given a certain
 *  value of "dist" found between an image pixel and the kernel
 *  center, the associated weight is interpolated linearly using
 *  the table.  If dot > kparam(1), a zero weight is used.
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
 *  The default means to use single pixels below 0.4, 2x2 blocks
 *  above 0.4, and 4x4 above 0.6.  For typical MDI images, this is
 *  16%, 20%, and 64% of on-disk pixels, respectively.  Supply
 *  bws=[] to treat all pixels singly (turn off blocking).
 *  * The optional output s is the z-coordinate (in [0,R]).
 *  The optional p output is the block (patch) number of each
 *  pixel.
 * 
 *  Inputs:
 *    real x[m,n];
 *    real center[3];  -- [center_x center_y r_sun]
 *    real p0;         -- p-angle (degrees)
 *    real k[p];
 *    real kparam[2];  -- [top window]
 *    opt real kwt = [1 1 1];
 *    opt real bws[bwnum,2] = [0.4 2;0.6 4];
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

// function entry point
mexfn_lib_t main_smoothsphere;

// argument counts
#define MXT_ssp_NARGIN_MIN 	5
#define MXT_ssp_NARGIN_MAX 	7
#define MXT_ssp_NARGOUT_MIN	0
#define MXT_ssp_NARGOUT_MAX	3

// input argument numbers
#define MXT_ssp_ARG_x	0
#define MXT_ssp_ARG_center	1
#define MXT_ssp_ARG_p0	2
#define MXT_ssp_ARG_k	3
#define MXT_ssp_ARG_kparam	4
#define MXT_ssp_ARG_kwt	5
#define MXT_ssp_ARG_bws	6

// output argument numbers
#define MXT_ssp_ARG_y	0
#define MXT_ssp_ARG_s	1
#define MXT_ssp_ARG_p	2


// (file ends)
