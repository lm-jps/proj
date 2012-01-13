/**************************************************************
 * 
 * mex2c_util definitions
 * 
 * generally useful utilities for mex2c, in either its command-line
 * or its xml/pleo incarnations
 *
 * input/output of numbers, strings, and matrices 
 * in text form as input to shell programs, for example
 *
 * version 3, Michael Turmon, 1999
 ****************************************************************/

#ifndef _mex2c_util_h_
#define _mex2c_util_h_

#include "mex.h"        /* mxArray data structure */

/************************************************************************
 * 
 * Function templates
 * 
 ***********************************************************************/

/* misc */
void
printerror(int status);

/* fits matrix i/o */
mxArray *
GetFITS(char *filename, int transposed);
void
PutFITS(mxArray *m, char *filename, int bitpix, int transposed, double scaling);

/* literal matrix i/o */
void
compute_range(double *data, long n, double *dmin, double *dmax);
mxArray *
GetNumericMatrix(char *value);
mxArray *
GetStringMatrix(char *value);
void
PutNumericMatrix2d(FILE *fp, 
		   int m,
		   int n,
		   double *pr, 
		   char *delim[]);
void
PutNumericArray(FILE *fp, 
		 mxArray *pm, 
		 char *head, char *tail, char *delim[]);
void
PutStringMatrix(FILE *fp,
		mxArray *pm, 
		char *head, char *tail, char *delim[]);

#endif /* _mex2c_util_h_ */
