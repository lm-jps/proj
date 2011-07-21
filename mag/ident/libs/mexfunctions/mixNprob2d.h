
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `mixNprob2d' a.k.a. `m2d' as a C library.
 *
 * Made by intermediate binary `mixNprob2d.out' on:
 * 	Thu Jul  7 23:02:43 2011
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 *  mixNprob2d: probabilities under normal mixture (2 dimensions)
 * 
 *  p=mixNprob2d(i1,i2,model,mode);
 *  * This routine has been superseded by mixNprobNd.  But we retain it
 *  because it can be much faster for 2d models.
 * 
 *  Inputs:
 *    real i1[m,n];
 *    real i2[m,n];
 *    real model[6,k];
 *    opt int logmode = 1;
 * 
 *  Outputs:
 *    real p[m,n];
 * 
 *  See Also: mixN2mixture2d, mixNprobNd
 * 
 *  implemented as a mex file
 * 
 *  MJT 20 july 1998: tested exhaustively against mixNprob,
 *  the m-file for finding gaussian mixture probabilities,
 *  and found agreement to within floating-point accuracy
 * 
 * 
 */

#ifndef _mexfn_mixNprob2d_h_
#define _mexfn_mixNprob2d_h_

// function entry point
mexfn_lib_t main_mixNprob2d;

// argument counts
#define MXT_m2d_NARGIN_MIN 	3
#define MXT_m2d_NARGIN_MAX 	4
#define MXT_m2d_NARGOUT_MIN	1
#define MXT_m2d_NARGOUT_MAX	1

// input argument numbers
#define MXT_m2d_ARG_i1	0
#define MXT_m2d_ARG_i2	1
#define MXT_m2d_ARG_model	2
#define MXT_m2d_ARG_logmode	3

// output argument numbers
#define MXT_m2d_ARG_p	0


#endif // _mexfn_mixNprob2d_h_

// (file ends)
