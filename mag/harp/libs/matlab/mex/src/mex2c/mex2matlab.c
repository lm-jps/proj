 /*
  * mex2matlab
  * 
  * library for interoperability in calling mex code from matlab, 
  * and from the command line or other interpreters. 
  *
  * these routines are not in the matlab C API, but had to be added
  * to get certain useful functionality.
  *
  */

#include <stdio.h>
#include "mex.h" /* must be present, but its contents are unused */
#include "mexhead.h"



/************************************************************************
 * 
 * dispatcher shell
 * 
 * (for compatibility with pure C library dispatcher)
 *
 * In contrast to the C dispatcher, calls within a library mexfunction
 * that is called from Matlab will route thru the Mathworks version
 * of mexErrMsgTxt.  This will generate an immediate return to the
 * matlab prompt.  So, there's no need for error handling here.
 *
 ***********************************************************************/


char *
mxt_mexDispatcher(mexfn_t *mexfn, 
		  const char *progname,
		  int verbosity,
		  int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
#ifdef DEBUG
  printf("  calling library mexfunction...\n");
#endif
  (*mexfn)(nlhs, plhs, nrhs, prhs);
#ifdef DEBUG
  printf("  returned normally from library mexfunction.\n");
#endif
  return NULL; /* successful return */
}


/************************************************************************
 * 
 * RANGE SETTING 
 * this is un-used by matlab -- hence not in the Matlab API as formulated
 * by Mathworks -- but required for some aspects of the command-line API.
 * 
 ***********************************************************************/


/*
 * set the range of the Matrix structure to [datamin,datamax].
 */
void    
setrange(mxArray *pm, double datamin, double datamax)
{ }

/*
 * get the range, [datamin,datamax], of the Matrix structure 
 * here datamin > datamax so the caller knows not to bother adjusting it;
 * if called from Matlab this will have no function
 */
void    
getrange(const mxArray *pm, double *datamin, double *datamax)
{ *datamin = 1.0; *datamax = -1.0; }

