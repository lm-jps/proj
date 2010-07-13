/*
 * mex_api_extras.h: mex header file
 *
 * function templates for routines/APIs that are NOT in the standard
 * matlab API, but that are needed to make mex operable from non-matlab
 * environments, like the shell.
 * 
 *
 * These routines are distinct from the rest of mexhead.h because
 * the underlying functional definitions must change depending on the 
 * environment (ie, call from matlab or not); the other mextools
 * utilities do not.
 * Because of this, their definitions are in the mex2c directory.  These
 * are just the declarations.  
 * The declarations must appear here because they must be loaded with 
 * mexhead.h.  This is because mexhead.h is the only include file, 
 * beyond mex.h, that I can be assured an ordinary mex2c-enabled 
 * function will load.  In short:
 *  declarations in mexhead.h (it's the only non-mex.h required #include)
 *  definitions in libmex2X.c (it's the only platform-specific library)
 *
 * Michael Turmon, 2002
 * 
 */

#ifndef _mex_api_extras_h_
#define _mex_api_extras_h_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * DISPATCHER
 */

char *
mxt_mexDispatcher(mexfn_t *mexfn, 
		  const char *progname,
		  int verbosity,
		  int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs);

 /* 
  * RANGE-SETTING TEMPLATES
  */

extern
void    
setrange(mxArray *pm, double datamin, double datamax);

extern
void    
getrange(const mxArray *pm, double *datamin, double *datamax);

#ifdef __cplusplus
}	/* extern "C" */
#endif

#endif /* _mex_api_extras_h_ */

