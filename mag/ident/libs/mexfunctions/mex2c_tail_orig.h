
/***************************************************************************
 * Generic include file providing hook to mexFunction
 * 
 * This file contains one function to interface
 * a mexFunction, declared static in the current namespace, to a caller
 * This file has been written to be generic across multiple mex files.
 * It expects that PROGNAME is #defined to be a unique string.
 * It expects to be #included from the file containing mexFunction, because
 * mexFunction must be declared static to get around namespace collisions.
 ***************************************************************************/

#ifndef _generic_mex2c_tail_h_
#define _generic_mex2c_tail_h_

// time functions
#include <sys/resource.h> 

/* 
 * Make the static mexFunction() here visible to caller
 * as main_PROGNAME(...) where PROGNAME is #defined in
 * the including file.
 */

/* utility: paste main_ plus PROGNAME into one token */
#define INITNAME1(X) INITNAME(X)
#define INITNAME(X) main_ ## X

void INITNAME1(PROGNAME) (int nlhs, 
			  mxArray **plhs, 
			  int nrhs, 
			  const mxArray **prhs)  
{
  struct rusage ru0;
  struct rusage ru1;
  double dt_u, dt_s; // system and user elapsed time

  printf("\tin local mexfunction (%s)\n", progname);
  getrusage(RUSAGE_SELF, &ru0); // tic
  mexFunction(nlhs, plhs, nrhs, prhs);
  getrusage(RUSAGE_SELF, &ru1); // toc
  dt_u = ((double) ru1.ru_utime.tv_sec + ru1.ru_utime.tv_usec * 1e-6) -
         ((double) ru0.ru_utime.tv_sec + ru0.ru_utime.tv_usec * 1e-6);
  dt_s = ((double) ru1.ru_stime.tv_sec + ru1.ru_stime.tv_usec * 1e-6) -
         ((double) ru0.ru_stime.tv_sec + ru0.ru_stime.tv_usec * 1e-6);
  printf("\tfinished local mexfunction (%s) in dt = %.3f s\n", progname, dt_u + dt_s);
}

#undef INITNAME1
#undef INITNAME

#endif /* _generic_mex2c_tail_h */
