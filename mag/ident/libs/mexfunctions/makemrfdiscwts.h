
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `makemrfdiscwts' a.k.a. `mdw' as a C library.
 *
 * Made by intermediate binary `makemrfdiscwts.out' on:
 * 	Thu Jun 16 23:36:44 2011
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * makemrfdiscwts	make distance metric for mrf segmentation
 * 
 *  dist=makemrfdiscwts(n,del,ctr,rho,mode,ell)
 *  * Given image sizes n, and neighborhood size del, make a distance
 *  metric array dist, such that dist(:,m1,n1) is the distance-to-
 *  neighbors of the site m1,n1, where 1 <= m1 <= n(1) and 1 <= n1 <= n(2).
 *  * Size of del = 3 corresponds to a 3x3 neighborhood around each site.
 *  * The distance is computed for a disk of center ctr(1:2) and radius
 *  ctr(3); ctr(1) corresponds to x or n1 and ctr(2) to y or m1.
 *  ctr(1:2) = [1 1] corresponds to the first pixel.
 *  * Thus, for instance, if Cx, Cy, r is the MDI center from mdidisk(),
 *  we would specify ctr=[Cx Cy r].  Note that ctr(1) and n(2) correspond
 *  to x, while ctr(2) and n(1) correspond to y.
 *  * Note that p-angle and b-angle are not needed by this routine.
 *  For simplicity, a 5-tuple containing these as entries 4 and 5
 *  (a `geom' vector) can be given for ctr, and the tail will be ignored.
 *  * The distance scale factor is rho.  Zero corresponds to uniform
 *  distances (all one) and large positive values correspond to
 *  separation-sensitive distances.  100 (20..200) is a typical value
 *  for images of radius 500.  To gain more insight, the raw distances
 *  (theta below) are returned when rho=NaN is given.  A histogram
 *  of theta values relative to exp(-theta*rho) will show whether most
 *  theta's land in the linear part of the exponential.
 *  * If rho < 0 is given, the scale factor is taken to be abs(rho), but
 *  the distances are rescaled so that the overall max equals 1.  Otherwise,
 *  the distances diminish across the entire disk as rho increases.
 *  * More explicitly, consider two points (s, s') within the 2d disk.
 *  They map to points on the unit sphere, say
 *    s = (x y z), s' = (x' y' z'),
 *  The angle theta (>0) between them is found, and the distance
 *  is returned as:
 *    dist = exp(-theta*rho)
 *  Because theta is very close to zero, the code in fact uses the chord
 *  distance rather than the arc distance between the two points.
 *  * Because the distance is symmetric, dist only includes distances
 *  between sites s,s' where s' is less than s.  This is the case if
 *  the pixel corresponding to s' comes before s in the memory footprint
 *  of an image.
 *  * Elliptical discs are allowed for.  In this case, ctr(3) is the
 *  major axis semidiameter, and ell(1) is that of the minor axis.
 *  ell(2) is the counterclockwise rotation, in degrees, of the major
 *  axis off of the "m" or "y" axis.  If not given, ell(1)=ctr(3) and
 *  ell(2)=0.  The implied condition ell(1) <= ctr(3) is not required
 *  by the code.
 *  * If either s or s' is off-disk, zero is put in the corresponding
 *  distance value.  This is required by mrf_segment_wts.
 *  * The required mode string switches between sesw (mode = 'sesw')
 *  or transposed (mode = 'sene') pixel ordering, which works as follows.
 *  * The normal HMI (and normal MDI) pixel ordering starts in the
 *  southeast corner, and the first scan line of pixels runs toward the
 *  southwest corner.  This is `sesw' ordering.  The transposed ordering
 *  is `sene'; this ordering is what we used for the JPL MDI processing.
 *  This is implemented via internal stride parameters.
 *  * This is implemented as a MEX file.  It agrees with makemrfdiscwts2.m
 * 
 *  Inputs:
 *    int n(1) or (2); -- if n is scalar, assume a square image
 *    int del;
 *    real ctr(3) or ctr(4) or ctr(5);
 *    real rho;
 *    string mode;
 *    opt real ell(2) = [ctr(3) 0];
 * 
 *  Outputs:
 *    real dist ((del*del-1)/2,n(1),n(2));
 * 
 *  See Also:  makemrfdiscwts2, mrf_segment_wts
 * 
 *  turmon oct 2006
 * 
 * 
 */

#ifndef _mexfn_makemrfdiscwts_h_
#define _mexfn_makemrfdiscwts_h_

// function entry point
mexfn_lib_t main_makemrfdiscwts;

// argument counts
#define MXT_mdw_NARGIN_MIN 	5
#define MXT_mdw_NARGIN_MAX 	6
#define MXT_mdw_NARGOUT_MIN	0
#define MXT_mdw_NARGOUT_MAX	1

// input argument numbers
#define MXT_mdw_ARG_n	0
#define MXT_mdw_ARG_del	1
#define MXT_mdw_ARG_ctr	2
#define MXT_mdw_ARG_rho	3
#define MXT_mdw_ARG_mode	4
#define MXT_mdw_ARG_ell	5

// output argument numbers
#define MXT_mdw_ARG_dist	0


#endif // _mexfn_makemrfdiscwts_h_

// (file ends)
