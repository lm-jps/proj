/*
 * mex_stubs.c
 * 
 * These are stubs for mex functions to allow functions written
 * under the mex environment to be used directly from C, 
 * without the rest of the matlab environment.
 *
 * Numeric arrays (of whatever type), string arrays, and logicals 
 * are fully supported.
 * Support for struct and cell arrays is spotty.
 * Sparse arrays are not supported.
 *
 * For function definition, please see Matlab External Interface Guide.
 * 
 */

/*LINTLIBRARY*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>

#include "mex.h"
#include "mexhead.h"  /* for mxt_getnand() */

/* for mexFunctionName(), set by main() */
extern char *mex2c_progname;

/*
 * Shared API functions
 * 
 * these are shared between mex_stubs and pymex_stubs
 * the first uses only API functions
 * the second is exotic unimplemented functions
 */

/* for runtime errors on unimplemented functions */
void 
mex2c_unimplemented(char *fn) {
  char msg[256];
  
  snprintf(msg, sizeof(msg), 
	   "%s: error: mex2c does not yet implement %s", 
	   mex2c_progname, fn);
  mexErrMsgTxt(msg);
}

#include "mex_stubs_api.c"    /* stubs using only API functions */
#include "mex_stubs_unimp.c"  /* universally unimplemented stubs */


/************************************************************************
 * 
 * mx* for memory allocation
 * 
 ***********************************************************************/

void *mxMalloc(size_t n)                { return malloc(n); }
void *mxCalloc(size_t n, size_t size)   { return calloc(n, size); }
void  mxFree(void *ptr)                 {        free(ptr); }
void *mxRealloc(void *ptr, size_t size) { return realloc(ptr, size); }

/************************************************************************
 * 
 * mx* for creation/destruction/duplication
 * 
 ***********************************************************************/

mxArray *
mxCreateNumericArray(mwSize n, const mwSize *dims, mxClassID id, mxComplexity cflag)
{
  mxArray *pm; /* result */
   
  /* make space for the basic object */
  if ((pm = (mxArray *) malloc(sizeof(mxArray))) == NULL)
    return NULL;
  /* set its type before anything else */
  pm->tag = id;
  pm->sparse = 0;
  /* set its size */
  if (mxSetDimensions(pm, dims, n)) {
    /* can fail to alloc space for the dimension vector */
    free(pm);
    return NULL;
  }
  /* set data range: indicate it is totally unknown */
  pm->datamin = mxt_getnand();
  pm->datamax = mxt_getnand();
  /* memory for real part */
  pm->data.numA.pr = (void *) 
    calloc((size_t) mxGetNumberOfElements(pm), mxGetElementSize(pm));
  if (pm->data.numA.pr == NULL) {
    free(pm->dims); /* dimension vector was allocated */
    free(pm);
    return NULL; /* calloc failed */
  }
  /* memory for imaginary part if necessary */
  if (cflag == mxCOMPLEX) {
    pm->data.numA.pi = (void *) 
      calloc((size_t) mxGetNumberOfElements(pm), mxGetElementSize(pm));
    if (pm->data.numA.pi == NULL) {
      free(pm->data.numA.pr); /* real data */
      free(pm->dims); /* dimension vector */
      free(pm);
      return NULL; /* calloc failed */
    }
  } else {
    pm->data.numA.pi = NULL;
  }
  return pm; /* return the matrix */
}

mxArray *
mxCreateNumericMatrix(mwSize m, mwSize n, mxClassID classid, mxComplexity cflag)
{
  mwSize dims[2];

  dims[0] = m; dims[1] = n;
  return mxCreateNumericArray((mwSize) 2, dims, classid, cflag);
}

mxArray *
mxCreateDoubleMatrix(mwSize m, mwSize n, mxComplexity cflag)
{
  mxArray *pm; /* result */
  mwSize dims[2]; /* space for dimensions */
   
  /* make space for the basic object */
  if ((pm = (mxArray *) malloc(sizeof(mxArray))) == NULL)
    return(NULL);
  /* set its type before anything else */
  pm->tag = mxDOUBLE_CLASS;
  pm->sparse = 0;
  /* set its size */
  dims[0] = m; dims[1] = n;
  if (mxSetDimensions(pm, dims, (mwSize) 2)) {
    /* can fail to alloc space for the dimension vector */
    free(pm);
    return NULL;
  }
  /* set data range: indicate it is totally unknown */
  pm->datamin = mxt_getnand();
  pm->datamax = mxt_getnand();
  /* memory for real part */
  pm->data.numA.pr = (void *) calloc((size_t) (m*n), sizeof(double));
  if (pm->data.numA.pr == NULL) {
    free(pm->dims); /* dimension vector was allocated */
    free(pm);
    return NULL; /* calloc failed */
  }
  /* memory for imaginary part if necessary */
  if (cflag == mxCOMPLEX) {
    pm->data.numA.pi = (void *) calloc((size_t) (m*n), sizeof(double));
    if (pm->data.numA.pi == NULL) {
      free(pm->data.numA.pr);
      free(pm->dims); /* dimension vector was allocated */
      free(pm);
      return NULL; /* calloc failed */
    }
  } else {
    pm->data.numA.pi = NULL;
  }
  return pm;
}

/*
 * Create a two-dimensional array to hold double-precision scalar.
 */
mxArray *
mxCreateDoubleScalar(double value)
{
  mxArray *pm;

  if ((pm = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 1, mxREAL))) {
    (mxGetPr(pm))[0] = value; /* set up the value */
    pm->datamin = pm->datamax = value; /* can set up the range too */
  }
  return pm; /* ok for NULL case too */
}

/* 
 * creation of logicals
 * (not adequately tested)
 */
mxArray *
mxCreateLogicalArray(mwSize n, const mwSize *dims) 
{
  mxArray *pm; /* result */
   
  /* make space for the basic object */
  if ((pm = (mxArray *) malloc(sizeof(mxArray))) == NULL)
    return NULL;
  /* set its type before anything else */
  pm->tag = mxLOGICAL_CLASS;
  pm->sparse = 0;
  /* set its size */
  if (mxSetDimensions(pm, dims, n)) {
    /* can fail to alloc space for the dimension vector */
    free(pm);
    return NULL;
  }
  /* set data range: for logical, it's known */
  pm->datamin = 0;
  pm->datamax = 1;
  /* memory for data */
  pm->data.logA.pl = (mxLogical *) 
    calloc((size_t) mxGetNumberOfElements(pm), sizeof(mxLogical));
  if (pm->data.logA.pl == NULL) {
    free(pm);
    return NULL; /* calloc failed */
  }
  return pm;
}

mxArray *
mxCreateLogicalMatrix(mwSize m, mwSize n) 
{
  mxArray *pm; /* result */
  mwSize dims[2]; /* space for dimensions */
   
  /* make space for the basic object */
  if ((pm = (mxArray *) malloc(sizeof(mxArray))) == NULL)
    return NULL;
  /* set its type before anything else */
  pm->tag = mxLOGICAL_CLASS;
  pm->sparse = 0;
  /* set its size */
  dims[0] = m; dims[1] = n;
  if (mxSetDimensions(pm, dims, (mwSize) 2)) {
    /* can fail to alloc space for the dimension vector */
    free(pm);
    return NULL;
  }
  /* set data range: for logical, it's known */
  pm->datamin = 0;
  pm->datamax = 1;
  /* memory for the data */
  pm->data.logA.pl = (mxLogical *) calloc((size_t) (m*n), sizeof(mxLogical));
  if (pm->data.logA.pl == NULL) {
    free(pm);
    return NULL; /* calloc failed */
  }
  return pm;
}

mxArray *mxCreateLogicalScalar(bool value) {
  mxArray *pm; /* result */

  if ((pm = mxCreateLogicalMatrix((mwSize) 1, (mwSize) 1))) {
    (mxGetLogicals(pm))[0] = value; /* set up the value */
  }
  return(pm);
}



/* reclaim space of *pm via free()
 * the shell of support for cells and structs is present, but
 * it needs some polishing.  as usual, no sparse support.
 */
void
mxDestroyArray(mxArray *pm) {
  mwSize i;

  if (pm == NULL)
    return; /* nothing to do */
  switch (mxGetClassID(pm)) {
  case mxDOUBLE_CLASS:
  case mxSINGLE_CLASS:
  case mxINT8_CLASS:
  case mxUINT8_CLASS:
  case mxINT16_CLASS:
  case mxUINT16_CLASS:
  case mxINT32_CLASS:
  case mxUINT32_CLASS:
  case mxINT64_CLASS:
  case mxUINT64_CLASS:
    free((void *) pm->data.numA.pr);
    free((void *) pm->data.numA.pi);
    break;
  case mxCHAR_CLASS:
    free((void *) pm->data.charA.pc);
    break;
  case mxLOGICAL_CLASS:
    free((void *) pm->data.logA.pl);
    break;
  case mxCELL_CLASS:
    for (i = 0; i < pm->Nelem; i++)
      mxDestroyArray(pm->data.cellA.pc[i]);  /* components */
    free((void *) pm->data.cellA.pc); /* array of component ptrs */
    break;
  case mxSTRUCT_CLASS:
    /* still need to decide how fields are stored in the ps[] array */
    for (i = 0; i < pm->Nelem; i++)
      mxDestroyArray(pm->data.structA.ps[i]); /* components */
    free((void *) pm->data.structA.ps); /* array of component ptrs */
    for (i = 0; i < pm->data.structA.nfield; i++)
      free((void *) pm->data.structA.fname[i]); /* field names */
    free((void *) pm->data.structA.fname); /* list of field names */
    break;
  case mxFUNCTION_CLASS:
  case mxOPAQUE_CLASS:
  case mxOBJECT_CLASS:
  default:
    break;
  }
  mxFree(pm->dims); /* dimension list */
  mxFree(pm);
}

/* copy *pm
 * the shell of support for cells and structs is present, but
 * it needs some polishing.
 */
mxArray *
mxDuplicateArray(const mxArray *pm) {
  mxArray *out;
  size_t length;
  mwSize i;

  if (pm == NULL) return(NULL); /* nothing to do */
  /* first set up the base struct */
  if ((out = malloc(sizeof(mxArray))) == NULL) return(NULL); /* alloc or die */
  memcpy((void *) out, (void *) pm, sizeof(mxArray));
  /* must malloc all pointers within the base struct */
  /* NB, mxSetDims does malloc space; returns nonzero for failure */
  if (mxSetDimensions(out, mxGetDimensions(pm), mxGetNumberOfDimensions(pm)))
    return(NULL);
  /* data part */
  length = mxGetNumberOfElements(pm) * mxGetElementSize(pm); /* bytes */
  switch (mxGetClassID(pm)) {
  case mxDOUBLE_CLASS:
  case mxSINGLE_CLASS:
  case mxINT8_CLASS:
  case mxUINT8_CLASS:
  case mxINT16_CLASS:
  case mxUINT16_CLASS:
  case mxINT32_CLASS:
  case mxUINT32_CLASS:
  case mxINT64_CLASS:
  case mxUINT64_CLASS:
    /* real */
    mxSetData(out, calloc(length, (size_t) 1));
    if (mxGetData(out) == NULL) return (NULL);
    memcpy(mxGetData(out), mxGetData(pm), length);
    /* imaginary, if necessary */
    if (mxIsComplex(pm)) {
      mxSetImagData(out, calloc(length, (size_t) 1));
      if (mxGetImagData(out) == NULL) return (NULL);
      memcpy(mxGetImagData(out), mxGetImagData(pm), length);
    }
    break;
  case mxCHAR_CLASS:
    out->data.charA.pc = (mxChar *) calloc(length, (size_t) 1);
    if (out->data.charA.pc == NULL) return (NULL);
    memcpy((mxChar *) out->data.charA.pc, 
	   (mxChar *)  pm->data.charA.pc, length);
    break;
  case mxCELL_CLASS:
    /* space for mxArray cell pointers */
    out->data.cellA.pc = (mxArray **) calloc(length, (size_t) 1);
    if (out->data.cellA.pc == NULL) return (NULL);
    for (i = 0; i < pm->Nelem; i++) 
      if ((                out->data.cellA.pc[i] = 
	   mxDuplicateArray(pm->data.cellA.pc[i])) == NULL)
	return(NULL); /* check as we go */
    break;
  case mxSTRUCT_CLASS:
    /* still need to decide how fields are stored in the ps[] array */
    /* space for mxArray struct pointers */
    out->data.structA.ps = (mxArray **) calloc(length, (size_t) 1);
    if (out->data.structA.ps == NULL) return (NULL);
    for (i = 0; i < pm->Nelem; i++) 
      if ((                out->data.structA.ps[i] = 
	   mxDuplicateArray(pm->data.structA.ps[i])) == NULL)
	return(NULL); /* check as we go */
    /* field name pointers */
    out->data.structA.fname = (char **) 
      calloc(out->data.structA.nfield, sizeof(char *));
    if (out->data.structA.fname == NULL) return (NULL);
    /* field names themselves */
    for (i = 0; i < out->data.structA.nfield; i++)
      if ((out->data.structA.fname[i] = 
	   strdup(pm->data.structA.fname[i])) == NULL)
	return(NULL);
    break;
  case mxFUNCTION_CLASS:
  case mxOPAQUE_CLASS:
  case mxOBJECT_CLASS:
  default:
    break;
  }
  return(out);
}

mxArray *
mxCreateString(const char *str)
{
  const char *sptr[1]; 

  sptr[0] = str; /* just to make sure we can take its address */
  return mxCreateCharMatrixFromStrings((mwSize) 1, sptr);
}


/*
 * NB: written to allow constituent strings str[0], str[1], etc., 
 * to have different length.  Pads the resulting rectangular
 * matrix with spaces.
 */
mxArray *
mxCreateCharMatrixFromStrings(mwSize m, const char **str)
{
  mxArray *pm; /* the returned char matrix */
  mxChar *dst; /* char data part of *pm */
  mwSize dims[2]; /* its size */
  size_t strlen_one, strlen_max; /* component lengths */
  mwSize m_ct; /* count in m (row) dimension */
  mwSize n_ct; /* count in n (col) dimension */
  int blankwrite; /* bool: filling with blanks or taking from **str */

  /* 1: Make the basic object */
  /* find its size */
  for (strlen_max = 0, m_ct = 0; m_ct < m; m_ct++) {
    strlen_one = strlen(str[m_ct]);
    if (strlen_one > strlen_max) strlen_max = strlen_one;
  }
  dims[0] = m; /* rows */
  dims[1] = strlen_max; /* columns */
  pm = mxCreateCharArray((mwSize) 2, dims);
  if (pm == NULL) return(pm); /* failure */
  /* 2: Populate it */
  dst = pm->data.charA.pc;
  for (m_ct = 0; m_ct < m; m_ct++) {
    for (blankwrite = 0, n_ct = 0; n_ct < strlen_max; n_ct++) {
      /* set to true if we see a null */
      blankwrite |= (str[m_ct][n_ct] == '\0');
      /* write mxChar or pad with blanks, note indexing */
      if (blankwrite)
	dst[m_ct + n_ct * m] = (mxChar) ' '; /* space character */
      else
	dst[m_ct + n_ct * m] = (mxChar) str[m_ct][n_ct]; /* source */
    }
  }  
  /* that's all */
  return(pm);  
}


mxArray *
mxCreateCharArray(mwSize n, const mwSize *dims)
{
  mxArray *pm; /* the returned char matrix */

  /* make space for the basic object */
  if ((pm = (mxArray *) malloc(sizeof(mxArray))) == NULL)
    return(NULL);
  /* set its type before anything else */
  pm->tag = mxCHAR_CLASS;
  /* set its size */
  if (mxSetDimensions(pm, dims, n)) {
    /* can fail to alloc space for the dimension vector */
    free(pm);
    return NULL;
  }
  /* set default data range (probably unneeded) */
  pm->datamin = (double) 0;
  pm->datamax = (double) 65535; /* 16-bit unsigned max */
  /* memory for data part */
  pm->data.charA.pc = (mxChar *) 
    calloc((size_t) mxGetNumberOfElements(pm), 
	   (size_t) mxGetElementSize(pm));
  if (pm->data.charA.pc == NULL) {
    free(pm);
    return NULL; /* calloc failed */
  }
  return pm; /* return the matrix */
}

int
mxGetString(const mxArray *pm, char *buf, mwSize buflen)
{
  mxChar *src;
  mwSize n, nxfer;

  if (!mxIsChar(pm))
    return 1; /* not a string array: failure */
  /* eliminate "can transfer nothing" case */
  if (buflen == 0)
    return 1; /* always fails because can't null-terminate */
  /* below here, buflen-1 is still >= 0 */
  nxfer = mxGetNumberOfElements(pm); /* number of MxChars to transfer */
  if (nxfer > (buflen-1))
    nxfer = buflen-1; /* can only xfer this many */
  /* transfer the characters */
  src = pm->data.charA.pc;
  for (n = 0; n < nxfer; n++)
    *buf++ = (char) *src++; /* converting from MxChar to char */
  *buf = '\0'; /* must null-terminate ourselves */
  /* get the return code right */
  if (mxGetNumberOfElements(pm) > (buflen-1))
    return 1; /* there was not enough room for everything */
  return 0; /* 0 on success (see doc) */
}

char *
mxArrayToString(const mxArray *pm)
{
  mxChar *src;
  char *dst; /* return value */
  mwSize n, nxfer;

  if (!mxIsChar(pm))
    return(NULL); /* not a string array: can do nothing */
  nxfer = mxGetNumberOfElements(pm); /* number of MxChars to transfer */
  /* make the space */
  dst = (char *) malloc((size_t) nxfer + 1); /* bytes, plus the null */
  if (dst == NULL) return NULL;
  /* transfer the characters */
  src = pm->data.charA.pc;
  for (n = 0; n < nxfer; n++)
    dst[n] = (char) src[n]; /* converting from MxChar to char */
  dst[nxfer] = '\0'; /* must null-terminate ourselves */
  return dst; 
}

/************************************************************************
 * 
 * mx* for get/set pairs
 * 
 ***********************************************************************/

/* class-id */
mxClassID mxGetClassID(const mxArray *pa) {return(pa->tag);}

/* these always return size_t's */
size_t mxGetM(const mxArray *pm) {return (size_t) pm->M;}
size_t mxGetN(const mxArray *pm) {return (size_t) pm->N;}

/* mxSetM and mxSetN must allow successive calls to correctly
 * initialize N, M, and Nelem, even if all were nonsense before.
 */
void
mxSetM(mxArray *pm, mwSize m) {
  mwSize i; 
  pm->dims[0] = pm->M = m;
  for (pm->Nelem = pm->M, pm->N = 1, i = 1; i < pm->ndim; i++) {
    pm->Nelem *= pm->dims[i];
    pm->N     *= pm->dims[i];
  }
}

void
mxSetN(mxArray *pm, mwSize n) {
  mwSize i; 
  pm->dims[1] = n;
  for (pm->Nelem = pm->M, pm->N = 1, i = 1; i < pm->ndim; i++) {
    pm->Nelem *= pm->dims[i];
    pm->N     *= pm->dims[i];
  }
}

/* this is an mwSize */
mwSize mxGetNumberOfDimensions(const mxArray *pa) {return (mwSize) pa->ndim; }

/* the spec is ambiguous, but we choose not to copy pa->dims */
const mwSize *
mxGetDimensions(const mxArray *pa) { 
  return (const mwSize *) pa->dims; 
}

/* this is a size_t */
size_t mxGetNumberOfElements(const mxArray *pa) { return (size_t) pa->Nelem;}

/* the spec is clear that a new copy should be made */
int 
mxSetDimensions(mxArray *pa, const mwSize *size, mwSize ndims) {
  mwSize *size2; /* holding area for mxCalloc */
  mwSize i;

  if (ndims < (mwSize) 2)
    return 1; /* failure */
  if ((size2 = mxCalloc((size_t) ndims, sizeof(mwSize))) == NULL)
    return 1; /* return code *is* nonzero for failure */
  memcpy((void *) size2, (void *) size, ndims*sizeof(mwSize));
  pa->dims = size2; /* could not do earlier because it might be NULL */
  pa->ndim = ndims;
  /* update convenience variables also */
  pa->M = size[0];
  for (i = 1, pa->Nelem = pa->M, pa->N = 1; i < ndims; i++) {
    pa->Nelem *= size[i];
    pa->N     *= size[i];
  }
  return 0; /* success */
}

mwIndex
mxCalcSingleSubscript(const mxArray *pa, mwSize n, const mwSize *subs) {
  mwSignedIndex k; /* allow < 0 */
  mwIndex index;

  if (n == 0) return 0; /* not actually possible in mex */
  /* this correctly handles n=1 and up */
  for (index = subs[n-1], k = n-2; k >= 0; k--)
    index = subs[k] + pa->dims[k] * index;
  return index;
}

/* real, imaginary parts */
void *
mxGetImagData(const mxArray *pm) {
  return (mxIsNumeric(pm) ? (void *) pm->data.numA.pi : NULL);
}

void *
mxGetData(const mxArray *pm) {
  return (mxIsNumeric(pm) ? (void *) pm->data.numA.pr : NULL);
}

mxChar *
mxGetChars(const mxArray *pm) {
  return (mxIsChar(pm) ? (void *) pm->data.charA.pc : NULL);
}

mxLogical *
mxGetLogicals(const mxArray *pm){
  return (mxIsLogical(pm) ? (void *) pm->data.logA.pl : NULL);
}

void mxSetImagData(mxArray *pm, void *newdata) {pm->data.numA.pi = newdata;}
void mxSetData(    mxArray *pm, void *newdata) {pm->data.numA.pr = newdata;}

double *
mxGetPi(const mxArray *pm) {
  return (mxGetClassID(pm) == mxDOUBLE_CLASS ? 
	  (double *) pm->data.numA.pi : NULL);
}

double *
mxGetPr(const mxArray *pm) {
  return (mxGetClassID(pm) == mxDOUBLE_CLASS ? 
	  (double *) pm->data.numA.pr : NULL);
}

void mxSetPi(mxArray *pm, double *pi) {pm->data.numA.pi = (void *) pi;}
void mxSetPr(mxArray *pm, double *pr) {pm->data.numA.pr = (void *) pr;}

double
mxGetScalar(const mxArray *pm) {
  if (mxGetClassID(pm) == mxDOUBLE_CLASS) {
    double *p_real = (double *) pm->data.numA.pr;
    return p_real[0];
  }
  else
    return mxt_getnand();
}

/************************************************************************
 * 
 * mxIsX functions
 * 
 * for the rest, see mex_stubs_api.c
 ***********************************************************************/

bool
mxIsComplex(const mxArray *pm)
{
  if (mxGetClassID(pm) == mxDOUBLE_CLASS)
    return (pm->data.numA.pi != NULL);
  else
    return 0;
}

bool
mxIsSparse(const mxArray *pm)
{
  /* NB: no support for sparse arrays */
  return (pm->sparse);
}


/************************************************************************
 * 
 * Stubs for mex routines
 *
 * These routines interact with Matlab, so they mostly can't be emulated
 * 
 * For function definition, please see Matlab External Interface Guide.
 * 
 ***********************************************************************/

/* this is defined by the invoker so we don't call exit() here */
void mxt_ErrHandler(const char *);

void
mexErrMsgTxt(const char *error_msg)
{
   fprintf(stderr, "%s\n", error_msg);
   mxt_ErrHandler(error_msg);  // this is back in the namespace of the invoker
}

void
mexWarnMsgTxt(const char *error_msg)
{
   fprintf(stderr, "%s\n", error_msg); /* TODO: would warn(3c) be better? */
}

void 
mexErrMsgIdAndTxt(const char *id, const char *error_msg,  ... )
{
  va_list args;

  va_start(args, error_msg);
  fprintf(stderr, "%s\n", id); /* Of the form Component:error */
  vfprintf(stderr, error_msg, args);
  fprintf(stderr, "\n");
  va_end(args);
  mxt_ErrHandler(id);  // this is back in the namespace of the invoker
}

void 
mexWarnMsgIdAndTxt(const char *id, const char *warn_msg,  ... )
{
  va_list args;

  va_start(args, warn_msg);
  fprintf(stderr, "%s\n", id); /* Of the form Component:warning */
  vfprintf(stderr, warn_msg, args);
  fprintf(stderr, "\n");
  va_end(args);
}


/* uses the external variable which must be defined in mexFunction */
const char *
mexFunctionName(void)
{
  return(strdup(mex2c_progname));
}


/************************************************************************
 * API to determine and set input and output data range
 *
 * these are needed for output to FITS files that aren't in a floating-
 * point format, i.e. their BITPIX > 0.  Floating point numbers in the
 * given range must be scaled to fit into an integer representation.
 ***********************************************************************/

void
getrange(const mxArray *pm, double *datamin, double *datamax)
{
  *datamin = pm->datamin;
  *datamax = pm->datamax;
}

void
setrange(mxArray *pm, double datamin, double datamax)
{
  pm->datamin = datamin;
  pm->datamax = datamax;
}

