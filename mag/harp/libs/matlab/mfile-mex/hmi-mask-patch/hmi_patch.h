
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `hmi_patch' a.k.a. `Hpat' as a C library.
 *
 * Made by intermediate binary `hmi_patch.out' on:
 * 	Fri Jul  8 01:30:07 2011
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * hmi_patch	driver for HMI patch finding
 * 
 *  [bb,s,yrgn,crit]=hmi_patch(y,mag,geom,active,ker,kwt,tau)
 *  * Find active-region patches in a mask image y, returning them as
 *  a list of bounding boxes bb, and a re-encoded mask image yrgn.
 *  Optionally returns the smoothed mask via crit.
 *  * If mag is given as [], the statistics s are not computed and
 *  the return value is empty.
 *  * The parameter active tells what parts of y are considered
 *  to be within-active-region: (y == active) identifies the active
 *  stuff to be combined into patches.
 *  * The parameters ker, kwt, and tau control grouping, see
 *  smoothsphere for more.
 *  * Depending on later needs, some other morphological parameters
 *  might be added so that very tiny ARs are removed.  Currently
 *  this is not needed.
 * 
 *  Inputs:
 *    real y(m,n)
 *    real mag(m,n) or (0,0)
 *    real geom(5)
 *    int active
 *    real ker(Nk)
 *    real kwt(3)
 *    real tau
 * 
 *  Outputs:
 *    int bb(nr,4)
 *    real stats(nr,28) or (0,0)
 *    real yrgn(m,n)
 *    opt real crit(m,n)
 * 
 *  See Also:  smoothsphere region_bb concomponent roi_stats_mag
 * 
 *  turmon oct 2009, june 2010, sep 2010
 * 
 * 
 */

#ifndef _mexfn_hmi_patch_h_
#define _mexfn_hmi_patch_h_

// function entry point
mexfn_lib_t main_hmi_patch;

// argument counts
#define MXT_Hpat_NARGIN_MIN 	7
#define MXT_Hpat_NARGIN_MAX 	7
#define MXT_Hpat_NARGOUT_MIN	2
#define MXT_Hpat_NARGOUT_MAX	4

// input argument numbers
#define MXT_Hpat_ARG_y	0
#define MXT_Hpat_ARG_mag	1
#define MXT_Hpat_ARG_geom	2
#define MXT_Hpat_ARG_active	3
#define MXT_Hpat_ARG_ker	4
#define MXT_Hpat_ARG_kwt	5
#define MXT_Hpat_ARG_tau	6

// output argument numbers
#define MXT_Hpat_ARG_bb	0
#define MXT_Hpat_ARG_stats	1
#define MXT_Hpat_ARG_yrgn	2
#define MXT_Hpat_ARG_crit	3


#endif // _mexfn_hmi_patch_h_

// (file ends)
