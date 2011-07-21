
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `hmi_segment' a.k.a. `Hseg' as a C library.
 *
 * Made by intermediate binary `hmi_segment.out' on:
 * 	Fri Jul  8 01:30:07 2011
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * hmi_segment	driver for HMI segmentation
 * 
 *  [y,s,post,nclean]=hmi_segment(xm,xp,edge,iter,T,beta,alpha,geom,rho,m1,...)
 *  * Integrated routine for deriving HMI segmentations.  Uses
 *  models (m1,m2,...), plus a magnetogram-photogram pair (xm,xp),
 *  to deduce an integer labeling.  Besides images and models,
 *  it also requires some disk parameters and labeling smoothness
 *  parameters.
 *  * The defaults for iter and T are inherited from mrf_segment_wts,
 *  and should be looked up there.
 *  * The extreme limb can be a problem for labeling; we set the last
 *  `edge' pixels to quiet (1), and return the number of such pixels
 *  modified, optionally, in nclean.
 * 
 *  Inputs:
 *    real xm(m,n)
 *    real xp(m,n)
 *    real edge(2)
 *    int iter(1) or (2)
 *    real T[0] or [1] or [2] or [3] or [4]
 *    real beta[1] or [K,K]
 *    real alpha[K] or [0] = []
 *    real geom(5)
 *    real rho
 *    real m1(l,k1)
 *    ...
 *    real mR(l,kR)
 * 
 *  Outputs:
 *    int y(m,n)
 *    opt real s(R,nS)
 *    opt real post
 *    opt int nclean
 * 
 *  See Also:  makemrfdiscwts, mrf_segment_wts, mixNprobNd,
 *             clean_edge_label, roi_stats_mag
 * 
 *  turmon oct 2009, june 2010, june 2011
 * 
 * 
 */

#ifndef _mexfn_hmi_segment_h_
#define _mexfn_hmi_segment_h_

// function entry point
mexfn_lib_t main_hmi_segment;

// argument counts
#define MXT_Hseg_NARGIN_MIN 	11
#define MXT_Hseg_NARGIN_MAX 	13
#define MXT_Hseg_NARGOUT_MIN	1
#define MXT_Hseg_NARGOUT_MAX	4

// input argument numbers
#define MXT_Hseg_ARG_xm	0
#define MXT_Hseg_ARG_xp	1
#define MXT_Hseg_ARG_edge	2
#define MXT_Hseg_ARG_iter	3
#define MXT_Hseg_ARG_T	4
#define MXT_Hseg_ARG_beta	5
#define MXT_Hseg_ARG_alpha	6
#define MXT_Hseg_ARG_geom	7
#define MXT_Hseg_ARG_rho	8
#define MXT_Hseg_ARG_m1	9
#define MXT_Hseg_ARG_m2	10
#define MXT_Hseg_ARG_m3	11
#define MXT_Hseg_ARG_m4	12

// output argument numbers
#define MXT_Hseg_ARG_y	0
#define MXT_Hseg_ARG_s	1
#define MXT_Hseg_ARG_post	2
#define MXT_Hseg_ARG_nclean	3


#endif // _mexfn_hmi_segment_h_

// (file ends)
