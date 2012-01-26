
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `concomponent' a.k.a. `ccp' as a C library.
 *
 * Made by intermediate binary concomponent.out on Mon May  3 17:24:30 2010
 *
 * Source file ../Gen-include.c last modified on Mon May  3 17:04:13 2010
 */

// Original documentation block
/*
 * concomponent: identify connected components
 * 
 *  [y,nr]=concomponent(x,nbr);
 *  * identifies connected components of a binary labeling x, storing
 *  the result as a like-sized image y.  The number of regions is
 *  optionally returned in nr.
 *  * On-region components of y are in the range 1..nr.  Inputs x
 *  of 0 or NaN are treated as off-region, but passed intact through
 *  to output region map as 0 or NaN.
 *  * Either the 4-pixel neighborhood (N/S/E/W) template can be used,
 *  or the 8-pixel neighborhood, which includes diagonal pixels.
 *  Default is 8.
 *  * This is implemented as a MEX file.
 * 
 *  Inputs:
 *    real x(m,n); -- 0, 1, or NaN
 *    opt int nbr = 8; -- 4 or 8
 * 
 *  Outputs:
 *    int y(m,n); -- 0, NaN, or 1..nr
 *    opt int nr;
 * 
 *  See Also: region_bb
 * 
 *  turmon may 2001
 * 
 * 
 */

// function entry point
mexfn_lib_t main_concomponent;

// argument counts
#define MXT_ccp_NARGIN_MIN 	1
#define MXT_ccp_NARGIN_MAX 	2
#define MXT_ccp_NARGOUT_MIN	0
#define MXT_ccp_NARGOUT_MAX	2

// input argument numbers
#define MXT_ccp_ARG_x	0
#define MXT_ccp_ARG_nbr	1

// output argument numbers
#define MXT_ccp_ARG_y	0
#define MXT_ccp_ARG_nr	1


// (file ends)
