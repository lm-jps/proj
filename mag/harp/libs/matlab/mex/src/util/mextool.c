/*
 * mextool.c: simple tools for Matlab MEX files.
 * 
 * all function names begin mxt_
 *
 * most functions here pertain to initializing, allocating or 
 * indexing of scalars, vectors, matrices, and 3d arrays.
 *
 * Michael Turmon, April 1997
 * Modified for Matlab 5 API, August 2000
 * enhanced 2002
 * 
 */

/*LINTLIBRARY*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "mextool.h"     /* function templates: current file */
#include "num_strings.h" /* for do_transpose */
#include "ieee_consts.h" /* for mxt_getnand() */
#include "mexargcheck.h" /* for IsEmpty, etc. */

/* forward declaration */
static void *
mxt_index_block_alloc(mwSize n_dbl, mwSize n_ptr, double **data);

/* mxt_PackSignature: wrap up information in the mexFunction's metadata
 * (NARG{IN,OUT}_{MIN_MAX}, in_names, in_specs, out_names, docstring)
 * into a mxArray for transmission out of the mexFunction.
 * This is used to print usage information (docstring) or to retrieve
 * input/output names (xio interface).
 * 
 * Returns NULL on failure, otherwise a mxArray * summarizing the "index".
 */

mxArray * 
mxt_PackSignature(mxt_Signature index, 
		  int in_lo, int in_hi, int out_lo, int out_hi, 
		  const char **in_names, 
		  const char **in_specs, 
		  const char **out_names, 
		  const char *docstring)
{
  mxArray *a = NULL;
 
  switch (index) {
  case mxt_SignatureNARG:
    /* return a 1x4 matrix of NARG* */
    a = mxCreateDoubleMatrix(1, 4, mxREAL);
    if (!a) break;
    mxGetPr(a)[0] = (double) in_lo;
    mxGetPr(a)[1] = (double) in_hi;
    mxGetPr(a)[2] = (double) out_lo;
    mxGetPr(a)[3] = (double) out_hi;
    break;
  case mxt_SignatureDocstring:
    a = mxCreateString(docstring);
    break;
  case mxt_SignatureInNames:
    a = mxCreateCharMatrixFromStrings(in_hi, in_names);
    break;
  case mxt_SignatureInSpecs:
    a = mxCreateCharMatrixFromStrings(in_hi, in_specs);
    break;
  case mxt_SignatureOutNames:
    a = mxCreateCharMatrixFromStrings(out_hi, out_names);
    break;
  }
  return a;
}

/* mxt_ArrayToStrings: convert mxArray to C char **
 * This is the inverse function of mxCreateCharMatrixFromStrings.
 * Note that the individual strings in "a" must be the same length
 * because they are stored as a matrix.  This routine eats whitespace
 * on the end of each string to compensage for the blank-padding done
 * when the strings are stored as a matrix.
 * As a convenience, the returned list of string pointers is one
 * entry longer than necessary, and the last character pointer is NULL.
 */
char **
mxt_ArrayToStrings(const mxArray *a)
{
  char **s;
  mxChar *sa; /* the data pointer within "a" points to mxChar's */
  mwSize M, N;   /* the size of a */
  mwSignedIndex m, n; /* allow to be < 0 */

  if (!a || !mxIsChar(a)) 
    return NULL;
  M = mxGetM(a);  /* number of strings */
  N = mxGetN(a);  /* their shared length */
  sa = (mxChar *) mxGetChars(a);  /* head of the char block in a */
  if ((s = calloc(M+1, sizeof(char *))) == NULL)
    return NULL;
  for (m = 0; m < M; m++) {
    /* make space for string #m (one extra for \0) */
    if ((s[m] = calloc(N+1, sizeof(char))) == NULL)
      return NULL;
    /* copy string m, which is not contiguous */
    for (n = 0; n < N; n++)
      s[m][n] = (char) sa[m+n*M]; 
    /* eat spaces from back of string m: N-1, ..., 0 */
    for (n = N-1; (n >= 0) && (s[m][n] == ' '); n--)
      s[m][n] = '\0';
  }
  return s;
}
  




/* mxt_put_matrix: for debugging, puts data into global workspace
 * as an mxn matrix.  takes only a pointer to the head of the
 * data, which must be contiguous, but could be of any size
 * or dimensions.
 * The data is put into the workspace as a 2d quantity named, for example,
 * x42 if name = "x" and inx = 42.  If inx < 0, inx is not used for
 * the name at all.
 */
int
mxt_put_matrix(char *name,
	      int inx,
	      double *dat,
	      mwSize m,
	      mwSize n)
{
  mxArray *arr;
  char full_name[mxMAXNAM];

  /* return 1; / * if uncommented, disables the routine */
  if ((arr = mxCreateDoubleMatrix(m, n, mxREAL)) == NULL)
    return 0; /* could not get the memory */
  memcpy(mxGetPr(arr), dat, m*n*sizeof(double)); /* array <- dat */
  if (inx >= 0)
    snprintf(full_name, sizeof(full_name), "%s%d", name, inx);
  else
    snprintf(full_name, sizeof(full_name), "%s", name);
  mexPutVariable("base", full_name, arr);
  mxDestroyArray(arr); /* eliminate the array */
  return 1;
}

 /*
  * mxt_make_scalar
  * From mxArray structure, extract a scalar.  If input is empty,
  * use default.  If not empty or a scalar, raise an error.
  */
double
mxt_make_scalar(
		const mxArray *pm,
		double default_val)

{
  if (!IsRealLooseScalar(pm)) {
    mexErrMsgTxt("! internal error: cannot mxt_make_scalar from that");
    return mxt_getnand(); /* not actually reached */
  } else if (IsEmpty(pm))
    return default_val;
  else
    return mxGetScalar(pm);
}


 /*
  * mxt_make_vector
  * From mxArray structure, extract a vector.  If input is shorter
  * than desired length M, use corresponding entry in default
  * vector.  If default runs out, fill remainder by replicating
  * last default.  Letting desired length equal -1 means use whatever
  * size the input is (except if this causes an error, below).
  * If input not a real vector (length zero OK) or has more entries
  * than desired, raise an error.
  * If output vector size does not equal requested size, a
  * block is calloc'ed , else original mxArray data is used.
  * This means that unless the caller knows for certain that
  * a size match will or will not occur, they cannot mxFree the 
  * returned vector.  It might or might not have been calloc'd
  * In particular, if M is given as -1 no calloc will occur.
  *
  * possibly useful functionality to add: if pm==NULL, make and return a 
  * copy of the default.
  */
double *
mxt_make_vector(
		const mxArray *pm,
		mwSignedIndex M,   /* can be < 0 */
		double *default_val,
		mwSize default_len)
{
  double *chunk;
  double *data_val;
  mwSize M_present;
  mwSize m;

  if (!IsRealLooseVector(pm)) {
    /* something very wrong */
    mexErrMsgTxt("! internal error: cannot make_vector from non-numerics");
    return NULL; /* not actually reached */
  } else if ((M != -1) && (mxGetNumberOfElements(pm) > M)) {
    /* more numbers supplied than needed */
    mexErrMsgTxt("! internal error: cannot make_vector, input too long");
    return NULL; /* not actually reached */
  } else if ((M < 0) || (M == mxGetNumberOfElements(pm)))
    /* just pass on what we have */
    return(mxGetPr(pm));
  else {
    /* fill in as needed (note M >= 0 now) */
    chunk = (double *) mxCalloc((size_t) M, sizeof(double));
    if (!chunk) return NULL;
    M_present = mxGetNumberOfElements(pm);
    data_val = mxGetPr(pm);
    for (m = 0; m < M; m++)
      if (m < M_present) /* normal copying */
	chunk[m] = data_val[m];
      else if (m < default_len) /* fill in from current default */
	chunk[m] = default_val[m];
      else if (default_len > 0) /* fill in from last default present */
	chunk[m] = default_val[default_len-1];
      else /* no defaults present: fill in with zero */
	chunk[m] = 0.0;
    return(chunk);
  }
}

 /*
  * mxt_make_matrix1
  * From mxArray structure, extract a matrix.  If input is scalar or 
  * null (and desired size is not scalar or null!), make a default
  * of a diagonal matrix with the supplied scalar, or the default
  * if input is empty, along the diagonal.  (Non-square target
  * outputs are OK.)  Letting desired length equal (-1,-1) means 
  * use whatever size the input is.
  * Note, this can accept a d>2 dimensional input -- in this case
  * a 2d matrix is made with all the input elements, and M = dims(1),
  * and N = prod(dims(2:end)).
  * If input not a real matrix of size equal to requested size,
  * and is neither empty nor a scalar, raise an error.
  * Returned result is a vector of pointers to the columns of
  * the input, so that m[a][b] is m(b,a) in the normal (row,col)
  * indexing scheme of mathematics.
  * In all cases, the pointer-vector is calloc'ed.  
  * If output matrix size does not equal requested size, a
  * block is calloc'ed for the matrix, else original mxArray data
  * is used.  Because of this, in either case all storage allocated
  * in this routine can be freed with one call to mxFree.
  * See also make_matrix2, which handle default_val differently.
  */
double **
mxt_make_matrix1(
		 const mxArray *pm,
		 mwSignedIndex M,
		 mwSignedIndex N,
		 double default_val)
{
  double *chunk;
  double **chunk2;
  int do_fillin, size_match;
  mwSize pm_N, pm_M, pm_MN;
  mwSize N_target, M_target;
  mwSize m, n;

  /* note, this allows making a 2d matrix from a d-dimensional input,
   * because mxGetM * mxGetN = mxGetNumberOfElements in this case */
  if (!IsFullRealArray(pm)) { /* something very wrong */
    mexErrMsgTxt("! internal error: cannot make_matrix1 from that");
    return NULL; /* not actually reached */
  }
  /* FIND USEFUL SIZES */
  /* size of input mxArray */
  pm_M = mxGetM(pm);
  pm_N = mxGetN(pm);
  pm_MN = pm_M * pm_N;
  /* size of data to be produced by this routine */
  N_target = (M == -1) ? pm_N : N; 
  M_target = (M == -1) ? pm_M : M;
  /* FILL DATA IF NEEDED */
  /* will size(in) == size(out) ? */
  size_match = (M == -1) || ((pm_M == M) && (pm_N == N));
  /* need to fill data (using scalar input, or the default_val)? */
  do_fillin = !size_match && (pm_MN <= 1);
  if (!size_match && !do_fillin) { /* bad # of numbers supplied */
    mexErrMsgTxt("error: mxt_make_matrix1: size conflict when defaulting");
    return NULL; /* not actually reached */
  }
  /* have detected all errors: now set up data */
  if (size_match) {
    chunk = mxGetPr(pm); /* easy case: use what is there */
    chunk2 = (double **) mxCalloc((size_t) N_target, 
				  sizeof(double *)); /* col ptrs */
  } else {
    /* here do_fillin == 1 => must fill in */
    /* note here N_target == N, etc. */
    /* get the block (chunk2) plus the data block (chunk) within it */
    chunk2 = (double **) mxt_index_block_alloc(M*N, N, &chunk);
    if (!chunk2) return NULL; /* failed calloc */
    /* now, fill in data */
    if (pm_MN == 1) default_val = mxGetScalar(pm); /* replace */
    for (n = 0; n < N; n++)
      for (m = 0; m < M; m++)
        chunk[n*M+m] = m == n ? default_val : 0.0; /* plug in to diags */
  }
  /* SET UP COLUMN POINTERS */
  for (n = 0; n < N_target; n++) 
    chunk2[n] = chunk + M_target*n;
  return chunk2;
}


 /*
  * mxt_make_matrix2
  * From mxArray structure, extract a matrix.  If input is scalar or 
  * null (and desired size is not scalar or null!), make a default
  * replicating the supplied scalar, or the default if input is empty, 
  * throughout the matrix.  Letting desired length equal (-1,-1) means 
  * use whatever size the input is.
  * Note, this can accept a d>2 dimensional input -- in this case
  * a 2d matrix is made with all the input elements, and M = dims(1),
  * and N = prod(dims(2:end)).
  * If input not a real matrix of size equal to requested size,
  * and is neither empty nor a scalar, raise an error (since a default
  * value is requested and a valid one cannot be found).
  * Returned result is a vector of pointers to the columns of
  * the input, so that m[a][b] is m(b,a) in the normal (row,col)
  * indexing scheme of mathematics.
  * In all cases, the pointer-vector is calloc'ed.  
  * If output matrix size does not equal requested size, a
  * block is calloc'ed for the matrix, else original mxArray data
  * is used.  Because of this, in either case all storage allocated
  * in this routine can be freed with one call to mxFree.
  * See also make_matrix1, which handle default_val differently.
  */
double **
mxt_make_matrix2(
		 const mxArray *pm,
		 mwSignedIndex M,
		 mwSignedIndex N,
		 double default_val)
{
  double *chunk;
  double **chunk2;
  int do_fillin, size_match;
  mwSize pm_N, pm_M, pm_MN;
  mwSize N_target, M_target;
  mwSize n;

  /* note, this allows making a 2d matrix from a d-dimensional input,
   * because mxGetM * mxGetN = mxGetNumberOfElements in this case */
  if (!IsFullRealArray(pm))  /* something very wrong */
    mexErrMsgTxt("! internal error: cannot mxt_make_matrix2 from that");
  /* FIND USEFUL SIZES */
  /* size of input mxArray */
  pm_M = mxGetM(pm);
  pm_N = mxGetN(pm);
  pm_MN = pm_M * pm_N;
  /* size of data to be produced by this routine */
  N_target = (M == -1) ? pm_N : N; 
  M_target = (M == -1) ? pm_M : M;
  /* FILL DATA IF NEEDED */
  /* will size(in) == size(out) ? */
  size_match = (M == -1) || ((pm_M == M) && (pm_N == N));
  /* need to fill data (using scalar input, or the default_val)? */
  do_fillin = !size_match && (pm_MN <= 1);
  if (!size_match && !do_fillin) /* bad # of numbers supplied */
    mexErrMsgTxt("error: mxt_make_matrix2: size conflict when defaulting");
  /* have detected all errors: now set up data from input or defaults */
  if (size_match) {
    chunk = mxGetPr(pm); /* easy case: use what is there */
    chunk2 = (double **) mxCalloc((size_t) N_target, 
				  sizeof(double *)); /* col ptrs */
    if (!chunk2) return NULL; /* failed calloc */
  } else  {
    /* here do_fillin == 1 => must fill in */
    /* note here N_target == N, etc. */
    /* get the block (chunk2) plus the data block (chunk) within it */
    chunk2 = (double **) mxt_index_block_alloc(M*N, N, &chunk);
    if (!chunk2) return NULL; /* failed calloc */
    /* now, fill in data */
    if (pm_MN == 1) default_val = mxGetScalar(pm); /* replace */
    for (n = 0; n < M*N; n++)
      chunk[n] = default_val; /* plug in all over the matrix */
  }
  /* SET UP COLUMN POINTERS */
  for (n = 0; n < N_target; n++) 
    chunk2[n] = chunk + M_target*n;
  return(chunk2);
}

/*
 * mxt_matrix_index: allocate column pointers for a m by n
 * matrix, returning the pointers.  They will enable indexing
 * the matrix via handle[j][i] where j < n and i < m.  (Note
 * the bracket reversal due to differences in C and Matlab indexing.)
 * Pointers are freed using a single call to mxfree(); this
 * does not free the data they point to!
 * Returns null on error (failed calloc)
 */
double **
mxt_matrix_index(double *data, mwSize m, mwSize n)
{
  double **handle;  /* whole block of memory (data + col ptrs) */
  mwSize j;

  /* room for col ptrs */
  handle = mxCalloc(n, sizeof(double *)); 
  if (!handle) return NULL; /* trouble */
  /* set up ptrs */
  for (j = 0; j < n; j++)
    handle[j] = data + j*m;
  return handle;
}


/*
 * mxt_matrix3_index: allocate slice and column pointers for a mXnXp
 * matrix, returning the pointers.  They will enable indexing
 * the matrix via handle[k][j][i] where k < p, j < n, i < m.  (Note
 * the bracket reversal due to differences in C and Matlab indexing.)
 * Both layers of pointers are freed using a single call to mxfree(); 
 * this does not free the data they point to!
 * Note that about n*p pointer slots are allocated.
 * Returns null on error (failed calloc)
 */
double ***
mxt_matrix3_index(double *data, 
		  mwSize m,
		  mwSize n, 
		  mwSize p)
{
  double ***handle;  /* whole block of memory (slice + col ptrs) */
  double  **handle2; /* start of col ptrs */
  size_t units;      /* number of bytes needed */
  mwSize j, k;          /* indexing */

  /* room for slice + col ptrs */
  units = p * sizeof(double **) + p*n * sizeof(double *);
  handle = mxCalloc(units, (size_t) 1); 
  if (!handle) return NULL;
  /* set up ptrs */
  handle2 = (double **) (handle + p); /* starts second layer of n*p pointers */
  for (k = 0; k < p; k++) {
    handle[k] = handle2 + k*n; /* leave room for n col ptrs at handle[k] */
    for (j = 0; j < n; j++)
      handle[k][j] = data + k*(n*m) + j*m;
  }
  return handle;
}


/*
 * mxt_matrix_alloc: allocate data and pointers for a m by n
 * matrix, returning the col pointers.  Both the pointers
 * and the data section are a contiguous block, so they 
 * are freed using a single call to mxfree()
 * If data is non-null, is is copied into the data segment
 * of the new matrix handle.
 * Returns null on error (failed calloc)
 */
double **
mxt_matrix_alloc(double *data, mwSize m, mwSize n)
{
  double **handle;  /* whole block of memory (data + col ptrs) */
  double *chunk;    /* points to head of numeric part */
  mwSize j;

  /* get the block (handle) plus the data block (chunk) within it */
  handle = (double **) mxt_index_block_alloc(m*n, n, &chunk);
  if (!handle) return NULL;
  /* set up ptrs */
  for (j = 0; j < n; j++)
    handle[j] = chunk + j*m;
  /* plug in data if it exists */
  if (data)
    memcpy(chunk, data, (size_t) m*n*sizeof(double));
  return handle;
}


/*
 * mxt_matrix3_alloc: allocate data and pointers for a mXnXp
 * matrix, returning the col pointers.  All the pointers
 * and the data section are a contiguous block, so they 
 * are freed using a single call to mxfree()
 * If data is non-null, is is copied into the data segment
 * of the new matrix handle.
 * Returns null on error (failed calloc)
 */
double ***
mxt_matrix3_alloc(double *data, mwSize m, mwSize n, mwSize p)
{
  double ***handle;  /* whole block of memory (data + col ptrs) */
  double **handle2;  /* points to head of column pointers (p*n ptrs) */
  double *chunk;     /* points to head of numeric part */
  mwSize j, k;

  /* room for data + col ptrs, slice ptrs */
  /* get the block (handle) plus the data block (chunk) within it */
  handle = (double ***) mxt_index_block_alloc(m*n*p, n*p+p, &chunk);
  if (!handle) return NULL; /* failed calloc */
  /* set up ptrs */
  handle2 = (double **) (handle + p);  /* go p double **'s beyond handle */
  for (k = 0; k < p; k++) {
    handle[k] = handle2 + k*n; /* leave room for n col ptrs at handle[k] */
    for (j = 0; j < n; j++)
      handle[k][j] = chunk + k*(n*m) + j*m;
  }
  /* plug in data if it exists */
  if (data)
    memcpy(chunk, data, (size_t) m*n*p*sizeof(double));
  return handle;
}

/*
 * Transposition
 */

/* 
 * transpose a matrix structure "in place"
 * 
 * The shell of the mxArray remains intact.  However, its data
 * segment is reset to point to the transposed data, and the old
 * segment is freed.
 */
int
mxt_transpose_double_mxArray(mxArray *pa)
{
  mwSize ndim;
  mwSize *old_dims;
  mwSize *new_dims;
  double *old_pr;
  double *new_pr;
  mwSize d;

  /* make flipped dimensions */
  ndim = mxGetNumberOfDimensions(pa); /* dimensionality */
  old_dims = (mwSize *) mxGetDimensions(pa);
  new_dims = (mwSize *) mxCalloc(ndim, sizeof(*new_dims)); 
  if (!new_dims) return(0); /* calloc fails */
  for (d = 0; d < ndim; d++)
    new_dims[d] = old_dims[ndim-d-1]; /* copy in reverse order */
  /* make flipped data */
  new_pr = (double *) mxCalloc(mxGetNumberOfElements(pa), sizeof(double));
  if (!new_pr) {
    mxFree(new_dims);
    return 0; /* calloc fails */
  }
  old_pr = mxGetPr(pa);
  if (!do_transpose(old_pr, new_pr, ndim, new_dims)) {
    mxFree(new_dims);
    mxFree(new_pr);
    return 0; /* trouble */
  }
  /* install flipped parameters */
  mxSetDimensions(pa, new_dims, ndim); /* new dimensions */
  mxSetPr(pa, new_pr); /* new data heap */
  /* free alloc'd memory */
  mxFree(old_pr); /* old data segment */
  mxFree(new_dims); /* mxSetDimensions alloc's its own space for dims */
  mxFree(old_dims); /* this has been overwritten */
  return 1; /* OK */
}


/*
 * mxt_index_block_alloc: utility routine, allocate ptr/data blocks
 * 
 * allocate storage for n_dbl doubles (data) that will be indexed by
 * n_ptr double* pointers.  Because the pointers may be of a different
 * length than the doubles they point to, to share the block we must be
 * careful about alignment.  After callocing the combined chunk, we
 * locate the right place within it to begin the data; this is output
 * via *data.  The head of the whole block is returned via the
 * function call (NULL for calloc failure).
 */

static void *
mxt_index_block_alloc(mwSize n_dbl, mwSize n_ptr, double **data)
{
  double **block_head; /* head of combined block */
  size_t units;        /* bytes of the combined ptr+data segment */
  size_t ptr_bytes;    /* bytes of the ptr segment alone */
  size_t bump;         /* difference to align the data part */

  ptr_bytes = n_ptr * sizeof(double *);
  /* bump ptr_bytes up to a multiple of sizeof(double) */
  bump = ptr_bytes % sizeof(double);
  if (bump != 0)
    ptr_bytes += sizeof(double) - bump;
  /* reset n_ptr to account for the bump */
  n_ptr = ptr_bytes / sizeof(double *);
  /* the requirement, in bytes */
  units = n_ptr * sizeof(double *) + n_dbl * sizeof(double);
  /* allocate the memory */
  block_head = (double **) mxCalloc(units, (size_t) 1);
  if (!block_head) return NULL;
  /* set up return values */
  *data = (double *) (block_head + n_ptr);  /* alignment IS ok */
  return (void *) block_head;    /* to a void* for return */
}
