
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `clean_edge_label' a.k.a. `cel' as a C library.
 *
 * Made by intermediate binary `clean_edge_label.out' on:
 * 	Thu Jul  7 22:00:02 2011
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 *  clean_edge_label: clean up labels at extreme edge of a disk
 * 
 *  [res,nclean]=clean_edge_label(img,center,delta,fill,mode);
 *  * Remove possibly-tainted labels at edge of disk by cleaning
 *  pixels in an annulus of width delta pixels.  Also, sets all
 *  off-disk values to the value of the (1,1) pixel of img.
 *  * There are two cleaning modes.  Easiest is, set values in
 *  the edge annulus to `fill' (NaN is OK).  But, if `mode'
 *  contains the string "adaptive", we instead propagate the value
 *  just inside the annulus radially outward to the limb, and the
 *  `fill' value is unused.
 *  * The number of on-disk pixels altered is in nclean.
 *  * Primary disk parameters (x-center, y-center, radius) are in
 *  `center'.  The sun is the disk defined by these three numbers.
 *  As a convenience, 'center' can be a "geom" vector, a 5-tuple
 *  with beta and p-angle at the end (these are not needed by
 *  this routine).
 *  * If delta < 0, img is propagated to res without further ado.
 *  (The off-disk clearing is not done either.)
 *  In this case, the center parameter need not be correct.
 *  * The "mode" string also switches between sesw (mode = 'sesw')
 *  or transposed (mode = 'sene') pixel ordering.  You must specify
 *  one or the other.
 *  * The normal HMI (and normal MDI) pixel ordering starts in the
 *  southeast corner, and the first scan line of pixels runs toward the
 *  southwest corner.  This is `sesw' ordering.  The transposed ordering
 *  is `sene'; this ordering is what we used for the JPL MDI processing.
 *  This is implemented via internal stride parameters.
 * 
 * 
 *  Inputs:
 *    real img(m,n);
 *    real center(3) or (5);
 *    real delta;
 *    real fill;
 *    string mode;
 * 
 *  Outputs:
 *    real res(m,n);
 *    int nclean;
 * 
 * 
 *  implemented as a mex file
 * 
 * 
 */

#ifndef _mexfn_clean_edge_label_h_
#define _mexfn_clean_edge_label_h_

// function entry point
mexfn_lib_t main_clean_edge_label;

// argument counts
#define MXT_cel_NARGIN_MIN 	5
#define MXT_cel_NARGIN_MAX 	5
#define MXT_cel_NARGOUT_MIN	1
#define MXT_cel_NARGOUT_MAX	2

// input argument numbers
#define MXT_cel_ARG_image	0
#define MXT_cel_ARG_center	1
#define MXT_cel_ARG_delta	2
#define MXT_cel_ARG_fill	3
#define MXT_cel_ARG_mode	4

// output argument numbers
#define MXT_cel_ARG_res	0
#define MXT_cel_ARG_nclean	1


#endif // _mexfn_clean_edge_label_h_

// (file ends)
