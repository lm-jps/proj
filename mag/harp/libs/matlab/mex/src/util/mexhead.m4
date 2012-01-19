dnl File to be expanded by m4
dnl It is OK to modify this file, but not the generated file

/*
 * mexhead.h: mex header file
 *
 * #defines and function templates useful regardless of whether
 * mex C-code is to be used standalone or from matlab.
 * 
 * Modified for Matlab 5 API, August 2000
 * 
 * Michael Turmon, 1996, 1997, 2000, 2002
 * 
 * (Note: this file was generated automatically by m4; do not edit directly)
 * 
 */

#ifndef _mexhead_h_
#define _mexhead_h_

/* 
 * Generic stuff 
 */
#define max(a,b)   (((a) > (b)) ? (a) : (b))
#define min(a,b)   (((a) < (b)) ? (a) : (b))

/* 
 * Useful function types 
 */
/* mexfn_t: a mexFunction() */
typedef void (mexfn_t)(int, mxArray **, int , const mxArray **);
/* mexfn_lib_t: mexFunction as a C library call */
typedef char * (mexfn_lib_t)(int, mxArray **, int , mxArray **);

include(mexargcheck.h)
include(mextool.h)
include(mex_api_extras.h)
include(ieee_consts.h)
include(num_strings.h)

#endif /* _mexhead_h_ */

