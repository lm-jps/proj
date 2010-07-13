
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `hmi_patch' a.k.a. `Hpat' as a C library.
 *
 * Made by intermediate binary `hmi_patch.out' on:
 * 	Wed Jun 23 11:25:43 2010
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * hmi_patch	driver for HMI patch finding
 * 
 *  [bb,s,yrgn]=hmi_patch(y,mag,ctr,p0,beta,active,ker,kwt,tau)
 *  * Find active-region patches in a mask image y, returning them as
 *  a list of bounding boxes bb, and a re-encoded mask image yrgn.
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
 *    real mag(m,n)
 *    real ctr(3)
 *    real p0
 *    real beta
 *    int active
 *    real ker(Nk)
 *    real kwt(3)
 *    real tau
 * 
 *  Outputs:
 *    int bb(nr,4)
 *    real stats(nr,22)
 *    real yrgn(m,n)
 * 
 *  See Also:  smoothsphere region_bb concomponent roi_stats_mag
 * 
 *  turmon oct 2009, june 2010
 * 
 * 
 */

#ifndef _mexfn_hmi_patch_h_
#define _mexfn_hmi_patch_h_

// function entry point
mexfn_lib_t main_hmi_patch;

// argument counts
#define MXT_Hpat_NARGIN_MIN 	9
#define MXT_Hpat_NARGIN_MAX 	9
#define MXT_Hpat_NARGOUT_MIN	2
#define MXT_Hpat_NARGOUT_MAX	3

// input argument numbers
#define MXT_Hpat_ARG_y	0
#define MXT_Hpat_ARG_mag	1
#define MXT_Hpat_ARG_ctr	2
#define MXT_Hpat_ARG_p0	3
#define MXT_Hpat_ARG_beta	4
#define MXT_Hpat_ARG_active	5
#define MXT_Hpat_ARG_ker	6
#define MXT_Hpat_ARG_kwt	7
#define MXT_Hpat_ARG_tau	8

// output argument numbers
#define MXT_Hpat_ARG_bb	0
#define MXT_Hpat_ARG_stats	1
#define MXT_Hpat_ARG_yrgn	2


#endif // _mexfn_hmi_patch_h_

// (file ends)
