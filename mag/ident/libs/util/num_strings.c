/*
 * num_strings.c
 * 
 * generally useful utilities involving numbers and strings
 *
 * needs mex.h for simple type definitions
 *
 * Michael Turmon, 2002, 2010
 */

/*LINTLIBRARY*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"         /* mwSize type, et al. */
#include "num_strings.h" /* function templates: current file */
#include "ieee_consts.h" /* mxt_getnand */

/**************************************************************
 * 
 * input/output of numbers, strings, and matrices 
 * in text form as input to shell programs, for example
 *
 * the subroutines do not need FITSIO declarations.
 *
 * they do need the simplest matlab "mwSize" typedef's
 *
 * string_to_array: converts a string representing an arbitrary-
 * dimensional array to the double-precision array
 * string_to_matrix: converts a string representing a matrix
 * to a double-precision matrix
 * matrix_to_string: converts a double-precision matrix to a 
 * string
 * convert_string_to_type: converts a string to type number
 * 
 ****************************************************************/

/* do_transpose: computes transpose, similar to matlab permute.
 * 
 * x: pointer to vectorized input
 * y: pointer to vectorized output (which is overwritten)
 * D: number of dimensions
 * By: bounds on each coordinate in y, as in mxGetDimensions
 *
 * Return value: y is written to.  A small malloc (size 2D+1) is
 * made, which may fail.  If so, 0 is returned, else 1.
 * 
 * This routine is nontrivial because it works for any dimensionality.
 *
 * Its behavior may be understood by letting ix and iy be the indexes
 * of corresponding elements of x and y.  If cx and cy are their
 * base-Bx and base-By digit expansions, we can write
 *
 * ix <--> cx[0],...cx[D-1]
 * iy <--> cy[0],...cy[D-1]
 *
 * where cx[0] = cy[D-1], etc., ie cy is the digit-reversed version
 * of cx.  So if you tell me one of these four quantities, I can find
 * the other three.
 * The loop unfolds by incrementing ix by 1 each time.  There's no
 * need to keep track of cx.  However, cy must be maintained to 
 * make the proper adjustments to iy.  In particular, ix++ corresponds 
 * to cx[0]++, ie to cy[D-1]++, ie to iy += delta_y where 
 * delta_y = By[0] * ... * By[D-2].  It's also useful to define
 * 
 * pi[0] = 1
 * pi[d] = By[0] * ... * By[d-1] 
 * 
 * so that pi[d] is the increment in iy due to an increment by one
 * in cy[d].  Note that pi[D] is the number of entries in the array.
 * The tricky point is carries: Some increments of ix cause a carry
 * "upward" in cx, which ripples downward in cy.  To adjust cy for
 * a carry, just set one of its digits to zero, and then increment
 * the previous one.  In iy, this corresponds to subtracting B*pi[d]
 * from iy, and adding pi[d-1] back in.
 * 
 * In July 2000 (matlab 5.3/linpack and 5.3.1/lapack) this was within 10%
 * of matlab's own ' for large, two-dimensional matrices.
 *
 * Michael Turmon, July 2000
 */
int
do_transpose(double *x,
	     double *y,
	     mwSize D,
	     mwSize *By)
{
  mwSize d; /* dimension counter */
  mwSize ix, iy; /* indexes into x and y */
  mwSize *pi = (mwSize *) calloc(D+1, sizeof(*pi));
  mwSize *cy = (mwSize *) calloc(D,   sizeof(*cy));

  if (pi == NULL || cy == NULL) return 0; /* failure */
  /* initial calculations */
  /* NB, pi[d] stores the increment in the scalar index implied by 
     incrementing digit d of the coordinate vector for y */
  for (pi[0] = 1, d = 1; d <= D; d++)
    pi[d] = pi[d-1] * By[d-1]; /* partial products */
  /* loop over all elements */
  ix = iy = 0; /* start both at origin (0,0,...0) */
  /* NB: all digits of cy are also 0 because of calloc */
  while (ix < pi[D]) {
    /* perform an assignment */
    y[iy] = x[ix];
    /* increment coordinates */
    ix++;           /* inc index into x */
    iy += pi[D-1];  /* inc index into y (increments MSD by 1) */
    cy[D-1]++;      /* inc MSD of y index, counterpart of above line */
    /* propagate carry down towards LSD if necessary */
    for (d = D-1; d > 0 && cy[d] == By[d]; d--) {
      /* there was overflow of digit-d --> digit-(d-1) */
      /* first, fix cy by carrying once towards LSD */
      cy[d] = 0; /* reset digit d */
      cy[d-1]++; /* increment digit d-1 (which may in turn overflow) */
      /* next, fix iy */
      /* step 1, subtract pi[d+1] = By[d] * pi[d].  This can be viewed
	 as resetting iy to its status before we began looping on this 
	 digit (doing By[d] iterations, each adding pi[d] to iy).  It
	 can also be viewed as the precise counterpart of the line above
	 setting cy[d] from By[d] to zero. */
      iy -= pi[d+1]; /* NB, d <= D-1, but pi[D] is valid */
      /* step 2, move up to the next element in this digit's loop by
	 adding pi[d-1] to the index.  This is the counterpart of the
	 line above, adding one to cy[d-1]. */
      iy += pi[d-1]; /* NB, d > 0, so we never reach too low into pi */
    }
  }
  return 1;
}

/* 
 * parse an array, input as a string
 *
 * string format:
 * array  := [array,array,...,array]
 *        := matrix
 * matrix := [a_11,a_12,...,a_1n;a_21,...,a_2n;...;a_m1,a_mn]
 * and in particular, the "matrix" is as indicated below in string_to_matrix.
 * Each component "array" must be the same size of all the other ones.
 * So, for example, a 1x1x1 array of the constant 42 is: [[42]].
 * Generally, the dimensionality is 1 plus the number of opening
 * left brackets (thus, always 2 or greater).
 * Suppose the dimensions of the array are d_1 x d_2 x ... d_N.
 * The number of arrays in the outermost group is d_1, and the number
 * in the innermost group is d_N-2.  Each of the latter arrays is 
 * exactly a "matrix" that is read in by string_to_matrix.  Thus,
 * in the input string, the "d_N" dimension changes fastest, while
 * in the output representation, the "d_1" dimension changes fastest.
 * The separating commas in the expansion of the array are not necessary
 * and may be omitted.
 * 
 * There are many places where errors are detected.  We use a flag and a
 * goto to handle errors: the flag is set to indicate error at the start
 * of the routine, and only set to "OK" when input is validated.  Any error
 * condition along the way branches to the error label, where all cleanup 
 * is done.  The error-free cases fall through to this cleanup stage also,
 * but with the flag set to indicate OK.  The cleanup is just freeing of
 * alloc'd variables -- they are 0 if they were never alloc'd.
 */
mwSize *
string_to_array(mwSize datamax,
		char *s, 
		mwSize *dimP, 
		double *data) 
{
  char *s1;          /* string position */
  char *s2 = 0;      /* temporary string (alloc'd) */
  double *data0 = 0; /* temporary for input values (alloc'd) */
  /* consistently use SignedIndex below because we use (-1) for flags */
  mwSignedIndex *dims = 0;  /* dimension size (to be returned) (alloc'd) */
  mwSignedIndex *dims0 = 0; /* dimension counters (alloc'd) */
  mwSignedIndex Nd;    /* number of dimensions */
  mwSignedIndex d;     /* counts dimensions */
  mwSignedIndex level; /* counter */
  mwSignedIndex m, Nm; /* number of submatrices */
  mwSignedIndex M, N;  /* submatrix dimensions */
  size_t l;            /* a string length */
  int convert_ok = 0;  /* assume error, only corrected if read OK */

  *dimP = 0; /* nonsense value */
  s2 = calloc(strlen(s)+1, sizeof(char)); /* guaranteed to be long enough */
  data0 = calloc(datamax, sizeof(double));
  if (!s2 || !data0) goto error; /* failed calloc */
  /* check for balanced brackets, and find dimensionality */
  for (Nd = level = 0, s1 = s; *s1; s1++) {
    if (*s1 == '[') level++;
    else if (*s1 == ']') level--;
    if (level < 0) goto error; /* too many ] */
    if (level > Nd) Nd = level; /* update depth */
  }
  /* level != 0 => unbalanced []. Nd == 0 => no brackets found */
  if (level != 0 || Nd == 0) goto error;
  Nd++; /* #dims = #brackets + 1.  Note: Nd >= 2 */
  /* obtain and check size */
  dims  = calloc((size_t) Nd, sizeof(*dims)); /* space for the size */
  dims0 = calloc(Nd, sizeof(*dims0)); /* counters to compare against above */
  if (!dims || !dims0) goto error;
  for (d = 0; d < Nd; d++) dims[d] = -1; /* dims = unknown */
  level = 0; /* use as a flag: was any non-bracket seen? */
  for (d = -1, s1 = s; *s1; s1++) {
    /* d: which level of bracket are we at (0 = outermost)
       dims0[d]: # bracketed groups at level d so far
       dims[d]:  # bracketed groups seen before at level d, -1 = don't know
    */
    if (*s1 == '[') {
      if (d == -1 && dims[0] != -1)
	goto error; /* only allow 1 element at the top level (disallow [][]) */
      d++;
      dims0[d] = 0; /* begin checking dimension d */
    } else if (*s1 == ']') {
      /* close the group at dimension d */
      if (dims[d] < 0)
	dims[d] = dims0[d]; /* no size yet: fill it in for later */
      if (dims[d] != dims0[d]) 
	goto error;  /* check the size */
      d--; /* move up to the enclosing group */
      if (d >= 0) /* NB: outermost group has d = -1 */
	dims0[d]++; /* we completed one more element in that group */
    } else if (*s1 == ',') {
      if (d == Nd-2)
	continue;  /* all commas in lowest grouping ok (check them later) */
      if (*(s1 - 1) == ']' && *(s1 + 1) == '[')
	continue; /* higher commas must appear so: ],[ */
      goto error; /* disallow [[1],[2],] and [,[1]] and [1], and ,[] etc. */
    } else {
      /* other characters: must be at bottom level */
      if (d != Nd - 2)
	goto error; /* other chars, not at bottom level: must be an error */
      level = 1; /* we did see other characters (=> non-empty matrix) */
    }
  }
  /* find cumulative sizes dims[0], dims[0]*dims[1], ... 
   * note that dims[Nd-1] and dims[Nd-2] are not known yet */
  for (d = 1, dims0[0] = dims[0]; d < Nd-2; d++)
    dims0[d] = dims[d] * dims0[d-1];
  /* special-case empty matrices: hard to isolate for string_to_matrix() */
  if (!level) {
    Nm = 0; /* no matrices to read */
    dims[Nd-1] = dims[Nd-2] = 0; /* set these up now */
  } else {
    Nm = Nd == 2 ? 1 : dims0[Nd-3]; /* number of matrices to read */
  }
  /* read in each component matrix */
  for (M = N = -1/* flag */, m = 0, s1 = s; m < Nm; m++) {
    /* skip s1 forward to first non-bracket */
    l = strspn(s1, "[],");
    s1 += l-1;  /* s1 = "[<matrix>]..." */
    /* extract [<matrix>] alone from s1 into s2 -- leave no trailing junk */
    l = strcspn(s1, "]");
    strncpy(s2, s1, l+1); s2[l+1] = '\0';
    /* read M*N elements into data starting at s1 */
    /* first iter: m=0 so m*M*N=0, which is good since M,N unknown then! */
    if (!string_to_matrix(datamax/Nm, s2, &M, &N, data0+m*M*N, 1))
      goto error; /* fail */
    /* check submatrix size */
    if (m == 0) {
      dims[Nd-1] = N; dims[Nd-2] = M; /* first submatrix sets size for later */
    } else {
      if (dims[Nd-1] != N || dims[Nd-2] != M) goto error; /* mismatch */
    }
    /* ready for next iter */
    l = strcspn(s1, "]");
    s1 += l; /* move up to ] */
  }
  /* check for junk at end of s1 */
  if (Nm > 0 && (strlen(s1) != strspn(s1, "]")))
    goto error; /* remainder must be only ]'s */
  /* transpose into matlab ordering */
  if (!do_transpose(data0, data, (mwSize) Nd, (mwSize *) dims)) goto error;
  *dimP = Nd; 
  convert_ok = 1; /* signal success */
 error:  /* free all memory that has been allocated, and exit */
  if (s2)    free(s2);
  if (data0) free(data0);
  if (dims0) free(dims0);
  if (dims && !convert_ok)
    free(dims); /* if error, should also free dims - it will be unused */
  /* if success: return pointer to the dimension list (it was calloc'd here) 
   *             cast it from (mwSignedIndex *) to (mwSize *)
   * if failure: return NULL 
   */
  return convert_ok ? ((mwSize *) dims) : NULL;
}

/*
 * parse a matrix, which is input as a string
 *
 * must be a <string> like:
 * <[a_11,a_12,...,a_1n;a_21,...,a_2n;...;a_m1,a_mn]>
 * where each a_ik = <whitespace><NUMBER><whitespace>
 * the whitespace can be null, and NUMBER is either
 * an ordinary floating-point constant or the
 * literal strings NaN, Inf, -Inf, which are interpreted 
 * as the corresponding IEEE special values.
 *
 * returns the computed size in Mp and Np, the 
 * numbers in data, and 0 or 1 according to 
 * success or failure.  The maximum number
 * of items allowed to convert is in datamax.
 */
int
string_to_matrix(mwSize datamax,
		 char *s, 
		 mwSize *Mp, 
		 mwSize *Np, 
		 double *data,
		 int raw) 
{
  char s1[MXT_IO_TMATRIX_STRLEN]; /* s without [ ] */
  char row1[MXT_IO_TMATRIX_STRLEN]; /* first row of s1 */
  char *delim, *delim2; /* placeholders */
  int nchar;
  int l = strlen(s);
  int row, col; /* row and column sizes */
  int ri, ci;   /* row, column indexing */
  int nscan; /* number of data read in so far */

  *Mp = *Np = 0;
  /* strip [] */
  if (l < 2)
    return 0; /* not enough chars */ 
  if (s[0] != '[' || s[l-1] != ']') 
    return 0; /* insist on no whitespace here */
  strcpy(s1, s+1); /* strips [ */
  s1[l-2] = '\0';  /* strips ] */
  /* find number of rows */
  row = 1; /* assume one row */
  delim = s1;
  while ((delim = strchr(delim, ';')) != NULL) {
    row++;
    delim++; /* advance past ; */
  }
  /* find number of columns */
  /* isolate the first row of s1 into row1 */
  for (delim2 = row1, delim = s1; *delim && *delim != ';'; )
    *delim2++ = *delim++;
  *delim2 = '\0'; /* terminate */
  /* now, as before with rows */
  col = 1; /* assume one col */
  delim = row1;
  while ((delim = strchr(delim, ',')) != NULL) {
    col++; 
    delim++; /* advance past , */
  }
  /* check for empty matrix */
  if (col == 1 && row == 1) {
    if (strspn(s1, " \t") == strlen(s1)) /* nothing but whitespace */
      col = 0; /* but let row = 1 still */
  }
  /* setup and check sizes */
  *Mp = row;
  *Np = col;
  if (col * row > datamax)
    return 0;
  /* scan entries, enforcing format */
  delim = s1; /* current pos in string */
  nscan = -1;  /* read first element */
  for (ri = 0; ri < row; ri++) 
    for (ci = 0; ci < col; ci++) {
      nscan = raw ? nscan + 1 : ci*row + ri; /* current target element */
      /* skip leading spaces */
      delim += strspn(delim, " \t");
      /* get the number */
      if (strncasecmp(delim, "nan", (size_t) 3) == 0) {
	data[nscan] = mxt_getnand();
	delim += 3; /* advance */
      } else if (strncasecmp(delim, "-inf", (size_t) 4) == 0) {
	data[nscan] = -1.0 * mxt_getinfd();
	delim += 4; /* advance */
      } else if (strncasecmp(delim, "inf", (size_t) 3) == 0) {
	data[nscan] = mxt_getinfd();
	delim += 3; /* advance */
      } else {
	if (sscanf(delim, "%lf%n", &data[nscan], &nchar) != 1)
	  return 0; /* no number present in this slot */
	delim += nchar; /* advance */
      }
      /* skip trailing spaces */
      delim += strspn(delim, " \t"); 
      /* check for proper delimiter */
      if (ci != col-1) {
	if (*delim != ',') /* need a , unless at end of row */
	  return 0;
      } else {
	if (ri != row-1)
	  if (*delim != ';') /* need a ; unless at end of matrix */
	    return 0;
      }
      delim++; /* advance past ; or , */
    }
  delim--; /* compensate for advance past nonexistent final delimiter above */
  if (*delim != '\0')
    return 0; /* junk at end */
  return 1;
}


/*
 * as above, but with int arguments, for compatibility with 
 * extract/insert_keywords
 */
int
string_to_matrix_ints(int datamax,
		      char *s, 
		      int *Mp, 
		      int *Np, 
		      double *data,
		      int raw) 
{
  mwSize Mp1, Np1;
  int result;

  result =
    string_to_matrix((mwSize) datamax, s, &Mp1, &Np1, data, raw);
  *Mp = Mp1;
  *Np = Np1;
  return result;
}


/*
 * output a string corresponding to a list of doubles (matrix)
 *
 * the string s must already be allocated by the calling program
 *
 * This routine is currently little used.
 * It has been largely replaced by PutNumericArray/PutNumericMatrix_level
 *
 * Outputs a <string> like:
 * <[a_11,a_12,...,a_1n;a_21,...,a_2n;...;a_m1,a_mn]>
 * where each a_ik = is either
 * an ordinary floating-point constant or the
 * literal strings NaN, Inf, -Inf, which are written for 
 * the corresponding IEEE special values.
 *
 * Brackets can be omitted.
 *
 */
void
matrix_to_string(double *data,
		 mwSize M, 
		 mwSize N,
		 int bracketed, /* 1 if surrounded by [] */
		 char *s) 

{
  int m, n; /* counters */
  int nchar; /* number of characters converted */

  /* opening [ */
  if (bracketed) *s++ = '[';
  /* the numbers */
  for (m = 0; m < M; m++) {
    for (n = 0; n < N; n++) {
      /* get the number */
      if (isnan(*data)) {
	strcpy(s, "NaN");
	s += 3;
      } else if (!finite(*data) && *data < 0) {
	strcpy(s, "-Inf");
	s += 4;
      } else if (!finite(*data) && *data > 0) {
	strcpy(s, "Inf");
	s += 3;
      } else {
	sprintf(s, "%.12g%n", *data, &nchar);
	s += nchar;
      }
      data++; /* move on */
      if (n != N-1)
	*s++ = ','; /* separating numbers: , */
    }
    if (m != M-1)
	*s++ = ';'; /* separating rows: ; */
  }
  /* closing ] */
  if (bracketed) *s++ = ']';
  *s++ = '\0';
}

/*
 * as above, but with ints for compatibility with extract_keywords
 */
void
matrix_to_string_ints(double *data,
		      int M, 
		      int N,
		      int bracketed,
		      char *s) 
{
  matrix_to_string(data, (mwSize) M, (mwSize) N, bracketed, s);
}

/*
 * convert a string to a type number
 */

int
mxt_io_convert_string_to_type(const char *type)
{
  int flag;

  /* account for a leading u flag */
  if (*type == 'u') {
    type++;  /* advance past it */
    flag = MXT_IO_UBRACK;   /* flag, if set, means "unbracketed" */
  } else {
    flag = 0;
  }
  /* OR the flag into the return value */
  if (strcmp(type, "TSHORT") == 0)
    return (flag | MXT_IO_TSHORT);
  else if (strcmp(type, "TINT") == 0)
    return (flag | MXT_IO_TINT);
  else if (strcmp(type, "TFLOAT") == 0)
    return (flag | MXT_IO_TFLOAT);
  else if (strcmp(type, "TDOUBLE") == 0)
    return (flag | MXT_IO_TDOUBLE);
  else if (strcmp(type, "TSTRING") == 0)
    return (flag | MXT_IO_TSTRING);
  else if (strcmp(type, "TMATRIX") == 0)
    return (flag | MXT_IO_TMATRIX); /* actually a matrixio type */
  else
    return MXT_IO_TBAD;
}

