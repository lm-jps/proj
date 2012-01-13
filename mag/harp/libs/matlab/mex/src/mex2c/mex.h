/*
 * mex.h
 * 
 * These are C replacements for mex data structures.  This allows
 * functions written for compilation and execution
 * under the mex environment to be used standalone from C, 
 * without the rest of the matlab environment.
 *
 * This file must be a direct replacement for the file
 * .../matlab/extern/include/mex.h
 *
 * version 3, Michael Turmon, 1999
 * version 4, Michael Turmon, 2002
 * version 4.1, Michael Turmon, 2003 (for matlab r13/v6.5)
 * version 5.0, Michael Turmon, 2008, add Python/numpy, light updates for v7
 * version 5.1, turmon, 2009, remove some matlab-specific mimicry
 * use mwSize, etc., for large arrays, Michael Turmon, 2010 (matlab R2009b)
 */

#ifndef _mex_h_
#define _mex_h_

#define USING_MEX2C /* so applications can determine environment */

/*
 * Python.h must be included before anything else.
 * mexfiles must include this mex.h *first*
 */
#ifdef USING_MEX2PY
#define PY_ARRAY_UNIQUE_SYMBOL mex2py_ARRAY_API
#include "Python.h"
#include "arrayobject.h"
#endif

#ifdef __cplusplus
extern "C" {
#ifdef NOT_DEFINED
} /* fool emacs */
#endif
#endif

/* mathworks mex.h includes these; for consistency, so do we */
#include <limits.h>   
#include <stddef.h>   /* size_t, ptrdiff_t */
#include <stdlib.h>   /* alloc, free */
#include <stdio.h>    /* for printf() below */
#include <assert.h>   /* for assert() below */
#include <float.h>    /* for DBL_EPSILON and other IEEE f.p. below */
#include <math.h>     /* ditto */

/* need functions/global variables to be static 
 * when compiling against mex2py (and anything else offering 
 * multi-mexfunction integration) because of namespace conflicts 
 * among multiple mex2* modules.  mex2c_cli does not have that issue, 
 * so here we just define StaticP to empty.
 * mnemonic: static or not-static = StaticP as in Lisp predicates
 */
#ifdef USING_MEX2PY
/* hook at end, used by mex2py to include methods table there */
#define MEX2C_TAIL_HOOK 
/* need certain things to be static when compiling against mex2py 
 * because of namespace conflicts among multiple mex2py modules */
#define StaticP static
#else /* !USING_MEX2PY */
/* allow an external definition for StaticP to override */
#ifndef StaticP
#define StaticP
#endif /* ifndef */
#endif

/*
 * constants
 */

#define mxMAXNAM 64  /* same as matlab v6 thru v7 */

/*
 * types, shadowing matlab
 */

/* boolean abstraction */
#if defined(NO_BUILT_IN_SUPPORT_FOR_BOOL)
typedef unsigned char bool;
#ifndef false
#define false (0)
#endif
#ifndef true 
#define true (1)
#endif
#else  /* defined(NO_BUILT_IN_SUPPORT_FOR_BOOL) */
#include <stdbool.h> /* typical case */
#endif /* defined(NO_BUILT_IN_SUPPORT_FOR_BOOL) */

typedef bool mxLogical;

/* indexing abstractions */
#define MX_COMPAT_32  /* undef for large (>2G) arrays, cf -largeArrayDims */
#ifdef MX_COMPAT_32
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;
#else
typedef size_t    mwSize;
typedef size_t    mwIndex;
typedef ptrdiff_t mwSignedIndex; /* when you need a signed type */
#endif

/* character abstraction */
#ifdef USING_MEX2PY
/* numpy stores its string type as char's */
typedef char mxChar;
#else /* !USING_MEX2PY */
typedef unsigned short mxChar; /* follow matlab: use a short int */
#endif

typedef enum {
  mxUNKNOWN_CLASS = 0,
  mxCELL_CLASS = 1,
  mxSTRUCT_CLASS,
  mxLOGICAL_CLASS,  /* new at r13 */
  mxCHAR_CLASS,
  mxSPARSE_CLASS,   /* obsolete by r13 */
  mxDOUBLE_CLASS,
  mxSINGLE_CLASS,
  mxINT8_CLASS,
  mxUINT8_CLASS,
  mxINT16_CLASS,
  mxUINT16_CLASS,
  mxINT32_CLASS,
  mxUINT32_CLASS,
  mxINT64_CLASS,
  mxUINT64_CLASS,
  mxFUNCTION_CLASS,
  mxOPAQUE_CLASS,
  mxOBJECT_CLASS
} mxClassID;

typedef enum {
    mxREAL,
    mxCOMPLEX
} mxComplexity;

/*
 * fundamental array data structure 
 */
#ifdef USING_MEX2PY

#define mxArray PyArrayObject /* for mex2py, it's easy */

#else /* !USING_MEX2PY */

/* this defintion is recursive */
struct _mxA {
  /* bookkeeping */
  mxClassID     tag;   /* class of contents: controls data union below */
  bool    sparse;      /* attribute: applies various content classes above */
  /* note: no fields are set up for the other sparse machinery */
  mwSize  ndim;        /* number of dimensions */
  mwSize *dims;        /* size of each dimension */
  mwSize  M;           /* dims[0] as a shortcut */
  mwSize  N;           /* prod_{1..d-1} dims[i] as a shortcut */
  mwSize  Nelem;       /* prod(dims) as a shortcut */
  /* extra info for mex2c */
  double  datamin;     /* data range -- for FITS output */
  double  datamax;
  /* data pointers, variant depending on tag above */
  /* note: no variants for function handles and opaques because
   * there are no mx* functions using them yet */
  union { 
    struct {               /* numeric array */
      void *pr;            /* real data */
      void *pi;            /* imag data */
    } numA;
    struct {               /* char array */
      mxChar *pc;          /* char data */
    } charA;
    struct {               /* logical array */
      mxLogical *pl;       /* logical data */
    } logA;
    struct {               /* struct array */
      struct _mxA **ps;    /* structure data */
      char **fname;        /* field names */
      int nfield;          /* number of fields (always int) */
    } structA;
    struct {               /* cell array */
      struct _mxA **pc;    /* cell data */
    } cellA;
  } data;
};

/* the actual type name as used */
typedef struct _mxA mxArray;

#endif /* USING_MEX2PY */


#ifdef USING_MEX2PY
/* Declaration for the C-API entry point
 * This is a convenience for use in the per-mexfunction glue file.
 */
/* this is an abbreviation for clarity here, not to be used elsewhere */
typedef void (mexfn_t_abbrev)(int, mxArray **, int , const mxArray **);

PyObject *mx_entry(mexfn_t_abbrev *, PyObject *, PyObject *,
		   int, int, int, int, 
                   const char **, const char **, const char **);
#endif

/*
 * Old version compatibility has been removed.
 * These V4/V5 functions are here #defined away from their
 * target functions.  This results in a compile error.  
 */
#define mxIsFull()         mxIsFull_is_obsolete
#define mxCreateFull()     mxCreateFull_is_obsolete
#define mxIsString()       mxIsString_is_obsolete
#define mxFreeMatrix()     mxFreeMatrix_is_obsolete
#define mxSetLogical()     mxSetLogical_is_obsolete
#define mxClearLogical()   mxClearLogical_is_obsolete
#define mxSetName()        mxSetName_is_obsolete
#define mxGetName()        mxGetName_is_obsolete

/*
 * mexPrintf and mxAssert's, as #defines 
 */
#define mexPrintf printf
#define mxAssert(expr,msg)	assert(expr)
#define mxAssertS(expr,msg)	assert(expr)

/*
 * function prototypes
 *
 * functions that were introduced in later versions
 * are defined for earlier versions anyway, to encourage migration
 */

void    mexErrMsgTxt( const char *error_msg); 
void    mexWarnMsgTxt(const char *error_msg); 
int     mexEvalString(const char *string); 
int     mexCallMATLAB(int nlhs, mxArray *plhs[], int nrhs, 
		      mxArray *prhs[], const char *name); 
void    mexLock(void);
void    mexUnlock(void);
bool    mexIsLocked(void);
int     mexAtExit(void (*exit_fcn)(void)); /* NI */

StaticP 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]); 
const char *mexFunctionName(void); 
int mexPutArray(mxArray *pa, const char *ws);
int mexPutVariable(const char *workspace, const char *name, const mxArray *parray);
/* NB: most mex* functions have been omitted for brevity */

/* memory */
void *mxMalloc(size_t n);
void *mxCalloc(size_t n, size_t size);
void  mxFree(void *ptr);
void *mxRealloc(void *ptr, size_t size);
/* get/set attribute */
mxClassID mxGetClassID(const mxArray *pa);
void    *mxGetData(const mxArray *pa);
void     mxSetData(mxArray *pa, void *newdata);
void    *mxGetImagData(const mxArray *pa);
void     mxSetImagData(mxArray *pa, void *newdata);
double  *mxGetPi(const mxArray *pm);
double  *mxGetPr(const mxArray *pm);
void     mxSetPi(mxArray *pm, double *pi);
void     mxSetPr(mxArray *pm, double *pr);
/* get contents */
double   mxGetScalar(const mxArray *pm); 
int      mxGetString(const mxArray *pa, char *buf, mwSize buflen); 
mxChar  *mxGetChars(const mxArray *pa);      /* get data, for chars */
mxLogical *mxGetLogicals(const mxArray *pa); /* get data, for logicals */
/* get/set sizes */
size_t   mxGetM(const mxArray *pm); /* matlab does use size_t here */
size_t   mxGetN(const mxArray *pm);
void     mxSetM(mxArray *pm, mwSize m);
void     mxSetN(mxArray *pm, mwSize n);
mwSize   mxGetNumberOfDimensions(const mxArray *pa);
const mwSize *mxGetDimensions(const mxArray *pa);
int      mxSetDimensions(mxArray *pa, const mwSize *size, mwSize ndims);
size_t   mxGetNumberOfElements(const mxArray *pa); /* is size_t */
size_t   mxGetElementSize(const mxArray *pa);
mwIndex  mxCalcSingleSubscript(const mxArray *pa, mwSize nsubs, const mwIndex *subs);
/* get/set sparse */
mwIndex *mxGetIr(const mxArray *pa);           /* NI */
void     mxSetIr(mxArray *pa, mwIndex *newir); /* NI */
mwIndex *mxGetJc(const mxArray *pa);           /* NI */
void     mxSetJc(mxArray *pa, mwIndex *newjc); /* NI */
mwSize   mxGetNzmax(const mxArray *pa);        /* NI */
void     mxSetNzmax(mxArray *pa, mwSize nzmax);/* NI */
/* get/set structs and cells.  field numbers are int's. */
mxArray *mxGetCell(const mxArray *pa, mwIndex i);           /* NI */
void     mxSetCell(mxArray *pa, mwIndex i, mxArray *value); /* NI */
int      mxGetNumberOfFields(const mxArray *pa);            /* NI */
int      mxGetFieldNumber(const mxArray *pa, const char *name);        /* NI */
mxArray *mxGetFieldByNumber(const mxArray *pa, mwIndex i, int f);      /* NI */
void     mxSetFieldByNumber(mxArray *pa,mwIndex i,int f, mxArray *val);/* NI */
mxArray *mxGetField(const mxArray *pa, mwIndex i, const char *f);      /* NI */
void     mxSetField(mxArray *pa,mwIndex i,const char *f,mxArray *val); /* NI */
const char *mxGetFieldNameByNumber(const mxArray *pa, int f);          /* NI */
const char *mxGetClassName(const mxArray *pa);                         /* NI */
int      mxSetClassName(mxArray *pa, const char *classname);           /* NI */
int      mxAddField(mxArray *pa, const char *fieldname);               /* NI */
void     mxRemoveField(mxArray *pa, int fnum);                         /* NI */
/* get/set objects */
mxArray *mxGetProperty(const mxArray *pa, mwIndex i, const char *name);/* NI */
void     mxSetProperty(mxArray *pa, mwIndex i, const char *name, const mxArray *val); /* NI */

bool     mxIsFromGlobalWS(const mxArray *pa);
bool     mxIsEmpty(  const mxArray *pa); 
bool     mxIsComplex(const mxArray *pa); 
bool     mxIsDouble( const mxArray *pa); 
bool     mxIsSparse( const mxArray *pa); 
bool     mxIsNumeric(const mxArray *pa); 
bool     mxIsChar(   const mxArray *pa); 
bool     mxIsCell(   const mxArray *pa);
bool     mxIsStruct( const mxArray *pa);
bool     mxIsOpaque( const mxArray *pa);
bool     mxIsFunctionHandle(const mxArray *pa);
bool     mxIsObject( const mxArray *pa);
bool     mxIsClass(const mxArray *pa, const char *name); /* NI */
bool     mxIsSingle( const mxArray *pa);
bool     mxIsLogical(const mxArray *pa);
bool     mxIsLogicalScalar(const mxArray *pa);
bool     mxIsLogicalScalarTrue(const mxArray *pa);
bool     mxIsInt8(   const mxArray *pa);
bool     mxIsUint8(  const mxArray *pa);
bool     mxIsInt16(  const mxArray *pa);
bool     mxIsUint16( const mxArray *pa);
bool     mxIsInt32(  const mxArray *pa);
bool     mxIsUint32( const mxArray *pa);
bool     mxIsInt64(  const mxArray *pa);
bool     mxIsUint64( const mxArray *pa);

/* creation, copying, destruction of mxArray's */
mxArray *mxCreateNumericArray(mwSize n, const mwSize *dims, mxClassID id, mxComplexity f);
mxArray *mxCreateNumericMatrix(mwSize m, mwSize n, mxClassID classid, mxComplexity f);
mxArray *mxCreateDoubleMatrix (mwSize m, mwSize n, mxComplexity f);
mxArray *mxCreateDoubleScalar(double value);
#define mxCreateScalarDouble(d) mxCreateDoubleScalar(d)
mxArray *mxCreateSparse(mwSize m, mwSize n, mwSize nzmax, mxComplexity f); /* NI */
mxArray *mxCreateSparseLogicalMatrix(mwSize m, mwSize n, mwSize nzmax); /* NI */
mxArray *mxCreateLogicalArray(mwSize ndim, const mwSize *dims);
mxArray *mxCreateLogicalMatrix(mwSize m, mwSize n);
mxArray *mxCreateLogicalScalar(mxLogical value);
char    *mxArrayToString(const mxArray *pa);
mxArray *mxCreateString(const char *str); 
mxArray *mxCreateCharArray(mwSize ndim, const mwSize *dims);
mxArray *mxCreateCharMatrixFromStrings(mwSize m, const char **str);
mxArray *mxCreateCellMatrix(mwSize m, mwSize n); /* NI */
mxArray *mxCreateCellArray(mwSize ndim, const mwSize *dims); /* NI */
mxArray *mxCreateStructMatrix(mwSize m, mwSize n,int nf,const char **names); /* NI */
mxArray *mxCreateStructArray(mwSize nd, const mwSize *dims, int nf, const char **names); /* NI */
mxArray *mxDuplicateArray(const mxArray *in);
void     mxDestroyArray(mxArray *pa); 

/* IEEE floating point */
#define mxGetEps() DBL_EPSILON /* in float.h */
#define mxGetInf getinfd  /* in ieee_consts.c -> libmextools.a */
#define mxGetNaN getnand  /* in ieee_consts.c -> libmextools.a */
#define mxIsFinite finite
#define mxIsInf(x) (!finite(x)) 
#define mxIsNaN isnan     /* in math.h */
int finite(double dsrc);  /* location of the declaration varies */

#ifdef __cplusplus
#ifdef NOT_DEFINED
{ /* fool emacs */
#endif
}	/* extern "C" */
#endif

#endif /* _mex_h_ */

