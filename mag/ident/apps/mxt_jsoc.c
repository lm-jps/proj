/*
 * mxt_jsoc
 * 
 * Routines for JSOC/DRMS and mextools interaction
 * 
 * For HMI data pipeline integration.
 * 
 * Michael Turmon, 2009
 */

/********************************************************
 *
 * Error handling
 *
 *********************************************************/

#include <stdlib.h> /* for exit() */
#include <setjmp.h> /* for mexFunction exit via mexErrMsgTxt */

#ifdef NOT_DEFINED

/* hook for within-mexFunction exception handling via setjmp/longjmp */
static jmp_buf MXT_exception_env;
static int MXT_exception_valid = 0; // is the above valid


/*
 * mxt_ErrHandler: called by mexErrMsgTxt to process an exception
 * that occurs within a mexFunction.
 * 
 * This hook can be defined differently for each interface
 * mode so that errors within the mexFunction could be handled
 * as desired (not just by exit()).
 */
void
mxt_ErrHandler(void)
{
#ifdef DEBUG
  printf("  in mxt_ErrHandler, returning to %s.\n", 
	MXT_exception_valid ? "caller in main()" : "system");
#endif
  /* jump back: the matching setjmp() is at entry to mexFunction */
  if (MXT_exception_valid)
    longjmp(MXT_exception_env, 1); // allow clean up
  else
    exit(1); // nowhere else to go
}

#endif

/********************************************************
 *
 * Routines for making mxArrays from input CmdParams
 *
 *********************************************************/

/* 
 * mxt_param_to_matrix: convert FLOATS to mxArray
 * 
 * The FLOATS are vectors, and the mxArrays are matrices.
 * So, we need to say what the dimensions should be.
 *   M: -1, or preferred #rows for the array
 *   N: -1, or preferred #cols for the array
 *   If M = N = -1, we just assume M = 1.
 * If the number of rows/columns does not divide evenly
 * into the number of elements in the FLOAT list, we return NULL.
 */
static
mxArray *
mxt_param_to_matrix(CmdParams_t *params, char *name, int M, int N)
{
  mxArray *pA;
  double *x;
  int i, Ntot;
  char name1[80];

  if (M < 0 && N < 0)
    M = 1; // make it a row vector
  // this encodes number of values in the list
  snprintf(name1, sizeof(name1), "%s_nvals", name);
  Ntot = params_get_int(params, name1);
  // establish dimensions
  if (M < 0) 
    M = Ntot / N;  // figure M out
  else if (N < 0)
    N = Ntot / M;  // figure N out
  if (M*N != Ntot)
    return NULL;   // dimensions were incompatible
  pA = mxCreateDoubleMatrix(M, N, mxREAL);
  if (!pA) 
    return pA;
  // set up contents
  x = mxGetPr(pA); // data of pA
  for (i = 0; i < Ntot; i++) {
    snprintf(name1, sizeof(name1), "%s_%d_value", name, i);
    x[i] = params_get_double(params, name1);
  }
  return pA;
}

/* 
 * mxt_param_to_scalar: convert FLOAT to mxArray (1x1)
 */
static
mxArray *
mxt_param_to_scalar(CmdParams_t *params, char *name)
{
  return mxCreateDoubleScalar(params_get_double(params, name));
}

/* 
 * mxt_param_to_string: convert FLOAT to character mxArray
 *
 * The mxArray will have M = 1.
 */
static
mxArray *
mxt_param_to_string(CmdParams_t *params, char *name)
{
  return mxCreateString(params_get_str(params, name));
}

/********************************************************
 *
 * Routines for mxArray / DRMS_Array_t conversions
 *
 *********************************************************/

/*
 * convert mxArray to drms array
 *
 * (now only covers one case, but this is the right abstraction)
 */
static 
DRMS_Array_t *
mxArray2drms_array(mxArray *in_a, DRMS_Type_t out_type) {

  int status;
  DRMS_Array_t *out_a;
  int out_dims[2];
  unsigned char *out_data;
  double *in_data = mxGetPr(in_a);
  const mwSize *in_dims = mxGetDimensions(in_a);
  size_t Ntot = mxGetNumberOfElements(in_a);
  size_t i;

  // set up sizes (assume 2d case)
  out_dims[0] = (int) in_dims[0]; 
  out_dims[1] = (int) in_dims[1];
  // make empty drms_array (space for data is allocated)
  out_a = drms_array_create(out_type, 2, out_dims, NULL, &status);
  if (status)
    return NULL;
  // copy data over
  out_data = (unsigned char *) out_a->data;
  for (i = 0; i < Ntot; i++)
    if (isnan(in_data[i]))
      out_data[i] = 0;  // nan input
    else
      out_data[i] = (unsigned char) in_data[i];
  return out_a;
}
