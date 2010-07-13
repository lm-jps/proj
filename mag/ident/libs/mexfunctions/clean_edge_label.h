
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `clean_edge_label' a.k.a. `cel' as a C library.
 *
 * Made by intermediate binary clean_edge_label.out on Mon May  3 17:24:17 2010
 *
 * Source file ../Gen-include.c last modified on Mon May  3 17:04:13 2010
 */

// Original documentation block
/*
 *  clean_edge_label: clean up labels at extreme edge of a disk
 * 
 *  [res,nclean]=clean_edge_label(img,center,delta,mode);
 *  * Remove possibly-tainted labels at edge of disk by examining an
 *  annulus of width delta.
 *  * Two modes are supported.  If mode < 0, cleaning is done
 *  by propagating the value just inside the annulus  outward.
 *  If mode >= 0, the values in the annulus are simply set to mode.
 *  * The number of pixels actually altered is in nclean.
 *  * Primary disk parameters (x-center, y-center, radius) are in center().
 *  The sun is the disk defined by these three numbers.
 *  * If delta = 0, img is propagated to res without further ado;
 *  in this case the center parameter need not be correct.
 * 
 *  Inputs:
 *    real img(m,n);
 *    real center(3);
 *    real delta;
 *    int mode;
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

// function entry point
mexfn_lib_t main_clean_edge_label;

// argument counts
#define MXT_cel_NARGIN_MIN 	4
#define MXT_cel_NARGIN_MAX 	4
#define MXT_cel_NARGOUT_MIN	1
#define MXT_cel_NARGOUT_MAX	2

// input argument numbers
#define MXT_cel_ARG_image	0
#define MXT_cel_ARG_center	1
#define MXT_cel_ARG_delta	2
#define MXT_cel_ARG_mode	3

// output argument numbers
#define MXT_cel_ARG_res	0
#define MXT_cel_ARG_nclean	1


// (file ends)
