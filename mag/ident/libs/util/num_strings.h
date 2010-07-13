/*
 * 
 * num_strings.c definitions
 * 
 * generally useful utilities for mex2c
 *
 * input/output of numbers, strings, and matrices 
 * in text form as input to shell programs, for example
 *
 * version 3, Michael Turmon, 1999
 */

#ifndef _num_strings_h_
#define _num_strings_h_

/************************************************************************
 * 
 * Constants
 * 
 ***********************************************************************/

/* codes for various FITS keyword input/output modes */
#define MXT_IO_UBRACK  0x100  /* mask, if set, means "unbracketed" */
#define MXT_IO_TBAD    0      /* bad type */
#define MXT_IO_TSHORT  1      /* short integer */
#define MXT_IO_TINT    2      /* integer */
#define MXT_IO_TFLOAT  3      /* float scalar */
#define MXT_IO_TDOUBLE 4      /* double scalar */
#define MXT_IO_TSTRING 5      /* string */
#define MXT_IO_TMATRIX 6      /* real matrix */
#define MXT_IO_TMATRIX_MAX 100  /* maximum #entries in a command-line
			     matrix (2 digit index) */
#define MXT_IO_TMATRIX_STRLEN (MXT_IO_TMATRIX_MAX*30)  /* largest string for the matrix */



/************************************************************************
 * 
 * Function templates
 * 
 ***********************************************************************/

/* misc */
int
do_transpose(double *x,
	     double *y,
	     mwSize D,
	     mwSize *By);

mwSize *
string_to_array(mwSize datamax,
		char *s, 
		mwSize *dimP, 
		double *data);
int
string_to_matrix(mwSize datamax,
		 char *s, 
		 mwSize *Mp, 
		 mwSize *Np, 
		 double *data,
		 int raw);
int
string_to_matrix_ints(int datamax,
		      char *s, 
		      int *Mp, 
		      int *Np, 
		      double *data,
		      int raw);
void
matrix_to_string(double *data,
		 mwSize M, 
		 mwSize N,
		 int bracketed, 
		 char *s);
void
matrix_to_string_ints(double *data,
		      int M, 
		      int N,
		      int bracketed, 
		      char *s);
int
mxt_io_convert_string_to_type(const char *type);

#endif /* _num_strings_h_ */

