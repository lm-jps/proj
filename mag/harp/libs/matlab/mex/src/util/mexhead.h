
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

/*
 * ARGUMENT CHECKING
 * (mexargcheck.c)
 * 
 * Michael Turmon, 2002
 * 
 */

#ifndef _mexargcheck_h_
#define _mexargcheck_h_

#ifdef __cplusplus
extern "C" {
#endif

/* simple size/type checking */
#define IsFullRealArray(pm) (mxIsDouble(pm) && \
                             !mxIsComplex(pm) && \
                             !mxIsSparse(pm))
#define IsFullRealMatrix(pm) (IsFullRealArray(pm) && \
                             (mxGetNumberOfDimensions(pm) <= 2))
#define IsRealVector(pm) (IsFullRealMatrix(pm) && \
                          (mxGetM(pm) == 1 || \
                           mxGetN(pm) == 1))
/* Vectors in broad sense */
#define IsRealLooseVector(pm) (IsFullRealMatrix(pm) && \
                          (mxGetM(pm) <= 1 || \
                           mxGetN(pm) <= 1))
#define IsLooseVector(pm) (mxGetM(pm) <= 1 || \
                           mxGetN(pm) <= 1)
#define IsRealScalar(pm) (IsFullRealMatrix(pm) && \
                           mxGetNumberOfElements(pm) == 1)
/* Scalar in broad sense */
#define IsRealLooseScalar(pm) (IsFullRealMatrix(pm) && \
                           mxGetNumberOfElements(pm) <= 1)
#define IsLooseScalar(pm) (mxGetNumberOfElements(pm) <= 1)
#define IsEmpty(pm)       (mxGetNumberOfElements(pm) == 0)

/*
 *  Function templates 
 */

int
mexargparse(int narg, 
	    const mxArray **args,
	    const char **names, const char **specs, const char **msgs, const char *gen_msg);
int start_sizechecking();
int sizeinit(const mxArray *pm);
int sizecheck(const char *msg, int num);
int sizecheck_msg(const char *msg, const char **argnames, int num);
int sizeagree(const mxArray *pm);
int sizeagreeM(const mxArray *pm);
int sizeagreeN(const mxArray *pm);
int sizeagreeMN(const mxArray *pm);
int sizesare(int *s, int d);
int sizeis(int m, int n);
int sizeis3(int m, int n, int p);
int sizeisM(int m);
int sizeisN(int n);
int sizeisMN(int m, int n);

#ifdef __cplusplus
}	/* extern "C" */
#endif

#endif /* _mexargcheck_h_ */

/*
 * ARGUMENT TRANSLATION TEMPLATES
 * (mextool.c)
 *
 * Michael Turmon, 2002
 * 
 */

#ifndef _mextool_h_
#define _mextool_h_

#ifdef __cplusplus
extern "C" {
#endif
#ifdef NOT_DEFINED 
} /* fool emacs */
#endif

/* HACK: turmon 8/2010: Stanford uses matlab 7.1, these
 * types were defined only as of matlab 7.3.
 * I'm putting them here even though they should go into
 * mex.h; this file goes into mexhead.h which is included
 * early in my mexfiles...
 */
#ifdef MATLAB_MEX_FILE
/* indexing abstractions */
#define MX_COMPAT_32  /* undef for large (>2G) arrays, cf -largeArrayDims
*/
#ifdef MX_COMPAT_32
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;
#else
typedef size_t    mwSize;
typedef size_t    mwIndex;
typedef ptrdiff_t mwSignedIndex; /* when you need a signed type */
#endif
#endif


/* If mexFunction is passed nrhs < 0, it causes special behavior, including
 * passing back information about the function signature of the mexFunction.
 * This enum lists the possibilities for returned information.
 * See also mxt_PackSignature(), which wraps the information up.
 */
typedef enum {
    mxt_SignatureNARG = 1,
    mxt_SignatureDocstring,
    mxt_SignatureInNames,
    mxt_SignatureInSpecs,
    mxt_SignatureOutNames
} mxt_Signature;

mxArray * 
mxt_PackSignature(mxt_Signature index, 
		  int in_lo, int in_hi, int out_lo, int out_hi, 
		  const char **in_names, 
		  const char **in_specs, 
		  const char **out_names, 
		  const char *docstring);

char **
mxt_ArrayToStrings(const mxArray *a);

int
mxt_put_matrix(char *name,
	       int inx,
	       double *dat,
	       mwSize m,
	       mwSize n);
double
mxt_make_scalar(
		const mxArray *pm,
		double default_val);
double *
mxt_make_vector(
		const mxArray *pm,
		mwSignedIndex M,
		double *default_val,
		mwSize default_len);
double **
mxt_make_matrix1(
		 const mxArray *pm,
		 mwSignedIndex M,
		 mwSignedIndex N,
		 double default_val);
double **
mxt_make_matrix2(
		 const mxArray *pm,
		 mwSignedIndex M,
		 mwSignedIndex N,
		 double default_val);

double **
mxt_matrix_index(double *data, mwSize m, mwSize n);
double ***
mxt_matrix3_index(double *data, mwSize m, mwSize n, mwSize p);
double **
mxt_matrix_alloc(double *data, mwSize m, mwSize n);
double ***
mxt_matrix3_alloc(double *data, mwSize m, mwSize n, mwSize p);

int
mxt_transpose_double_mxArray(mxArray *pa);

#ifdef __cplusplus
}	/* extern "C" */
#endif

#endif /* _mextool_h_ */


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


/*
 * GENERATING IEEE FLOATING-POINT CONSTANTS
 * (ieee_consts.c)
 *
 * Michael Turmon, 2002
 * 
 */

#ifndef _ieee_consts_h_
#define _ieee_consts_h_

#ifdef __cplusplus
extern "C" {
#endif

float  mxt_getnanf(void);
double mxt_getnand(void);
double mxt_getinfd(void);

#ifdef __cplusplus
}	/* extern "C" */
#endif

#endif /* _ieee_consts_h_ */

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



#endif /* _mexhead_h_ */

