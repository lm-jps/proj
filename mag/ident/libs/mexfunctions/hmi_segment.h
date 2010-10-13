
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `hmi_segment' a.k.a. `Hseg' as a C library.
 *
 * Made by intermediate binary `hmi_segment.out' on:
 * 	Mon Oct 11 16:21:54 2010
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * hmi_segment	driver for HMI segmentation
 * 
 *  y=hmi_segment(xm,xp,iter,T,beta,alpha,ctr,rho,m1,...)
 *  * Integrated routine for deriving HMI segmentations.  Uses
 *  models (m1,m2,...), plus a magnetogram-photogram pair (xm,xp),
 *  to deduce an integer labeling.  Besides images and models,
 *  it also requires some disk parameters and labeling smoothness
 *  parameters.
 *  * The defaults for iter and T are inherited from mrf_segment_wts,
 *  and should be looked up there.
 * 
 *  Inputs:
 *    real xm(m,n)
 *    real xp(m,n)
 *    int iter[1] or [2];
 *    real T[0] or [1] or [2] or [3] or [4];
 *    real beta[1] or [K,K];
 *    real alpha[K] or [0] = [];
 *    real ctr(3) or ctr(4) or ctr(5);
 *    real rho;
 *    real m1(l,k1)
 *    ...
 *    real mR(l,kR)
 * 
 *  Outputs:
 *    int y(m,n)
 *    opt real post
 * 
 *  See Also:  makemrfdiscwts, mrf_segment_wts, mixNprobNd, clean_edge_label
 * 
 *  turmon oct 2009, june 2010
 * 
 * 
 */

#ifndef _mexfn_hmi_segment_h_
#define _mexfn_hmi_segment_h_

// function entry point
mexfn_lib_t main_hmi_segment;

// argument counts
#define MXT_Hseg_NARGIN_MIN 	10
#define MXT_Hseg_NARGIN_MAX 	12
#define MXT_Hseg_NARGOUT_MIN	1
#define MXT_Hseg_NARGOUT_MAX	2

// input argument numbers
#define MXT_Hseg_ARG_xm	0
#define MXT_Hseg_ARG_xp	1
#define MXT_Hseg_ARG_iter	2
#define MXT_Hseg_ARG_T	3
#define MXT_Hseg_ARG_beta	4
#define MXT_Hseg_ARG_alpha	5
#define MXT_Hseg_ARG_ctr	6
#define MXT_Hseg_ARG_rho	7
#define MXT_Hseg_ARG_m1	8
#define MXT_Hseg_ARG_m2	9
#define MXT_Hseg_ARG_m3	10
#define MXT_Hseg_ARG_m4	11

// output argument numbers
#define MXT_Hseg_ARG_y	0
#define MXT_Hseg_ARG_post	1


#endif // _mexfn_hmi_segment_h_

// (file ends)
