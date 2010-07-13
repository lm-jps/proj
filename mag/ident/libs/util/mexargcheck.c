/* 
 * mexargcheck.c
 *
 * check out arguments to mex-file
 *
 * two components, arg type checker and size checker
 *
 * Michael Turmon, 1996, rev. 1997, rev. 2000 (Matlab 5)
 */

/*LINTLIBRARY*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "mex.h"
#include "mexargcheck.h"

#define ERRBUFSIZE 120

/* declarations for later static functions */
static int mexargparse1(const mxArray *arg, const char *spec);
static int mexspecmatch(const mxArray *arg, char tspec, char sspec,
	     int *sizes, int d);
static int translate_spec(char *buf, const char *spec);

/* 
 * TYPE CHECKER ENTRY POINT
 *
 * check a list of arguments, printing error message if bad arg found
 *
 * the args are matched one by one against a corresponding spec, which
 * for one arg is:
 *
 * <spec> = <one-spec>               -- only one possibility
 *        = <one-spec> `|' <spec>    -- more than one template acceptable
 * <one-spec> = <type-spec><shape-spec><size-spec>
 * <type-spec> = 'R'                 -- real numbers
 *             = 'I'                 -- integers
 *             = 'B'                 -- booleans (0/1 only)
 *             = 'S'                 -- string
 *             = 'C'                 -- cells
 *             = 'F'                 -- struct (with Fields)
 *             = 'X'                 -- any type
 * <shape-spec> = 'A'                -- array (or empty, scalar, vector, matrix)
 *              = 'M'                -- matrix (or empty, scalar, vector)
 *              = 'V'                -- vector (or empty, scalar)
 *              = 'S'                -- scalar (or empty)
 * <size-spec> = ''                  -- no specification
 *             = '(' <M> ')'         -- primarily for vector and scalar
 *             = '(' <M> ',' <N> ')' -- for matrix types
 *             = '(' <M1> ',' ... <Mp> ')' -- for array types
 *
 * To specify a scalar integer that can't be empty, you could use
 * IS(1) or IS(1,1); for a real vector of size two, use RV(2).  If
 * it must be a row vector, use RV(1,2).  A two-dimensional matrix
 * is RM.  The generic numeric array is RA,
 * which just means a numeric quantity that is not complex, in particular
 * it could be empty.
 */
int
mexargparse(
	    int narg,      /* number to check */
	    const mxArray **args, /* arg list itself */
	    const char **names,  /* names of each arg above */
	    const char **specs,  /* spec for each arg above */
	    const char **msgs,   /* special per-arg message, NULL ok */
	    const char *gen_msg) /* generic message to be printed */

{
  char errbuf[ERRBUFSIZE];
  char errmsg[ERRBUFSIZE];
  int a;
  int perr;
	
  for (a = 0; a < narg; a++) {
    if (args[a] == NULL) break;
    if ((perr = mexargparse1(args[a], specs[a])) != 0) {
      /* found an ordinary error */
      if (perr > 0) {
	translate_spec(errbuf, specs[a]);
	sprintf(errmsg, "%s: Arg %d = %s, %s %s%s%s", 
		gen_msg, a+1, names[a], 
		msgs ? (*msgs[a] == '\0' ? "" : "(") : "", 
		msgs ? msgs[a] : "", 
		msgs ? (*msgs[a] == '\0' ? "" : ")") : "", 
		errbuf);
      }
      /* found a bad specification */
      else {
	sprintf(errmsg, 
		"%s: internal error: parse problem, typespec <%s> for arg %d",
		gen_msg,
		specs[a],
		a+1);
      }
      mexErrMsgTxt(errmsg);
      return 0; /* not reached: value not important */
    }
  }
  return 1; /* value not important */
}


/*
 * parse one argument
 * 
 * returns: 0 for spec match, 
 *         +1 if argument does not match spec
 *         -1 if spec error, 
 */
static int
mexargparse1(
	     const mxArray *arg, 
	     const char *spec)
{
  int index;       /* finds single clause within spec */
  char typespec, shapespec; /* tags for type and shape */
  char *spec_size; /* copy of size part of one clause */
  char *spec1;     /* finds single size in spec_size */
  int *sizes;      /* list of identified sizes */
  int d, dim;      /* length of above list */
  int spec_matched;
  int match;
	
  index = 0;
  spec_matched = 0;
  while ((!spec_matched) && (spec[index] != '\0')) {
    /* printf("** matching <%s>\n", spec+index); */
    /* find type, shape (easy) */
    typespec = spec[index++];
    if (!spec[index]) return(-1); /* early termination of spec */
    shapespec = spec[index++];
    /* Find size if present (harder).  Remainder must be either empty,
       or the start of another spec ("|" + SPEC), or a size spec
       of form (s1,...,sD) terminated by | or NULL */
    if (spec[index] == '\0') {
      sizes = (int *) NULL; d = 0; /* no size info at all */
    } else if (spec[index] == '|') {
      index++;  /* no size, so just move to the spec beyond the | */
      sizes = (int *) NULL; d = 0; /* no size info at all */
    } else if (spec[index] == '(') {
      /* there is a size -- update index, then get it */
      /* copy size to operate comfortably */
      spec_size = strdup(spec+index+1); /* start after the opening paren */
      /* move index beyond the current spec */
      while (spec[index] != '\0' && spec[index] != '|')
	index++;
      if (spec[index] == '|') index++; /* move beyond |, sigh */
      /* isolate just the s1,s2,...,sD part */
      if ((spec1 = strchr(spec_size, (int) '|')) != 0)
	*spec1 = '\0'; /* more than one spec: kill subsequent specs */
      else
	spec1 = spec_size + strlen(spec_size); /* one spec: go to its end */
      /* either way, now spec1 points to the terminating NULL of spec_size */
      spec1--; /* now, points to the last char */
      if (*spec1 == ')')
	*spec1 = '\0'; /* kill the closing paren */
      else
	return(-1); /* bad format: no closing paren */
      /* there is a size -- scan to find how many dims it has */
      for (spec1 = spec_size, d = 1; 1 /*CONSTCOND*/; d++)
	if ((spec1 = strchr(spec1, (int) ',')) == NULL)
	  break;
	else
	  spec1++;   /* ok, move 1 beyond the , we found */
      /* make space */
      sizes = mxCalloc(d, sizeof(int));
      /* scan again to find the size of each dim */
      for (spec1 = spec_size, dim = 0; dim < d; dim++) {
	if (sscanf(spec1, "%d", &(sizes[dim])) != 1)
	  return(-1); /* could not find a number */
	spec1 = strchr(spec1, (int) ',') + 1; /* 1 beyond , */
      }
      /* no more need for spec_size */
      free(spec_size);  /* strdup uses malloc() */
      /* done finding sizes */
    } else {
      return(-1); /* neither \0, nor |, nor (...) */
    }
    /* OK, now find if there was a match */
    match = mexspecmatch(arg, typespec, shapespec, sizes, d);
    /* printf("matching %c %c %d: %d\n", typespec, shapespec, d, match); */
    if (sizes) mxFree(sizes); /* get rid of sizes, too */
    /* take care of match info */
    if (match < 0)
      return(match); /* spec error */
    else
      spec_matched |= match;
    }
  return(spec_matched ? 0 : 1);
}

/*
 * does the matrix arg match the specification given?
 * 0: no match, +1: match, -1: invalid spec
 */
static int
mexspecmatch(
	     const mxArray *arg,
	     char tspec,
	     char sspec,
	     int *sizes,  /* -1 in an entry if no match specified */
	     int d)       /* d=0 if no matches at all */
	
{
  double one_elem;
  
  /* TYPE SPEC */
  switch (toupper(tspec))  {
  case 'R': /* real */
    if (!IsFullRealArray(arg))
      return(0);
    break;
  case 'I': /* integer */
    if (!IsFullRealArray(arg))
      return(0);
    else { /* check if we can */
      if (!IsEmpty(arg)) {
	one_elem = *mxGetPr(arg);
	if (!isnan(one_elem))
	  if (one_elem != ((int) one_elem))
	    return(0);
      }
    }
    break;
  case 'B': /* boolean */
    if (!IsFullRealArray(arg))
      return(0);
    else { /* check if we can */
      if (!IsEmpty(arg)) {
	one_elem = *mxGetPr(arg);
	if (!isnan(one_elem))
	  if ((one_elem != 0) && (one_elem != 1))
	    return(0);
      }
    }
    break;
  case 'S': /* string type */
    if (!mxIsChar(arg))
      return(0);
    break;
  case 'F': /* struct (with "Fields") type */
    if (!mxIsStruct(arg))
      return(0);
    break;
  case 'C': /* cell array type */
    if (!mxIsCell(arg))
      return(0);
    break;
  case 'X': /* any type */
    break;
  default:
    return(-1); /* error in spec! */
  }
  /* SHAPE SPEC */
  switch (toupper(sspec))  {
  case 'A':
    /* anything's an array in our view */
    break;
  case 'M':
    /* a matrix has no more than two dimensions */
    if (mxGetNumberOfDimensions(arg) > 2)
      return(0);
    break;
  case 'V':
    if (!IsLooseVector(arg))
      return(0);
    break;
  case 'S':
    if (!IsLooseScalar(arg))
      return(0);
    break;
  default:
    return(-1); /* error in spec! */
  }
  /* SIZE SPEC */
  if (d <= 0) {
    /* EMPTY */
    /* omit checking and fall through */
  } else if (d == 1) {
    /* vector: 1xsize or sizex1 is OK */
    if (mxGetNumberOfDimensions(arg) > 2)
      return(0); /* need #dims at most 2 */
    /* can fail if it's not [(1xsize) or (sizex1)] */
    if (!((mxGetM(arg) == 1 && mxGetN(arg) == sizes[0]) ||
	  (mxGetN(arg) == 1 && mxGetM(arg) == sizes[0])))
      if (sizes[0] >= 0) /* ok, also need size >= 0 so it's not excluded */
	return(0); 
  } else if (d >= 2) {
    int D = mxGetNumberOfDimensions(arg);
    const mwSize *dnums = mxGetDimensions(arg);
    int dim;
    if (d != D)
      return(0); /* #dims must match */
    for (dim = 0; dim < D; dim++)
      if (sizes[dim] >= 0 && sizes[dim] != dnums[dim])
	return(0); /* each un-excluded dim must match */
  }
  return(1);
}

/* 
 * translate a character string spec, see above, to english text
 */
static int
translate_spec(
	       char *buf,
	       const char *spec)
{
  int index;
  int index2;
	
  index = 0;
  *buf = '\0';
  while (spec[index]) {
    /* type spec  */
    switch (toupper(spec[index++])) {
    case 'R':
      strcat(buf, "Real");
      break;
    case 'I':
      strcat(buf, "Integer");
      break;
    case 'B':
      strcat(buf, "Boolean");
      break;
    case 'S':
      strcat(buf, "Char");
      break;
    case 'F':
      strcat(buf, "Struct");
      break;
    case 'C':
      strcat(buf, "Cell");
      break;
    case 'X':
      strcat(buf, "Any-type");
      break;	
    default:
      strcat(buf, "!internal error, bad type spec");
      break;
    }
    /* shape spec */
    switch (toupper(spec[index++]))  {
    case 'A':
      strcat(buf, " array");
      break;
    case 'M':
      strcat(buf, " matrix");
      break;
    case 'V':
      strcat(buf, " vector");
      break;
    case 'S':
      strcat(buf, " scalar");
      break;
    default:
      strcat(buf, "!internal error, bad shape spec");
      break;
    }
    /* size spec */			
    index2 = strlen(buf);
    while (spec[index] != '\0' && spec[index] != '|')  
      buf[index2++] = spec[index++];
    if (spec[index] == '|')
      buf[index2++] = spec[index++]; /* get that one too */
    buf[index2] = '\0';
  }
  return 1;
}
	
/* 
 * SIZE CHECKER ENTRY POINT
 *
 * calling sequence for a single argument, matrix1:
 *   sizeinit(array1);          -- initialize: we are checking array1
 *   sizeagreeXX(array2);       -- any number of reps of size templates...
 *   sizeisXX(number);          -- ...of any sort
 *   ...and so on...
 *   sizecheck_msg(message, names, ARG_NO);
 *                               -- check if some check matched
 *
 * Calls are two types, sizeagreeXX for "agreement" of array sizes, and
 * sizeisXX for explicit matching of a size to a given integer.  
 * The XX phrase is replaced by '' if both sizes must agree in number and 
 * position, by 'M' if only the first size is to agree, by 'N' if the 
 * second, and by 'MN' if both must agree but either orientation is OK.
 *
 * Thus sizeis(1, 10) asks for a 1-by-10 matrix, but sizeisMN(1,10)
 * also allows a 10-by-1 matrix.  
 * 
 * Array arguments can be checked by sizeagree() and sizesare().
 *
 * If no agreement has been found among the several calls by the 
 * time sizecheck_msg is called, an eror will be printed.
 *
 * The start_sizechecking() call must be made before the first 
 * sizeinit(m) call.  The overall format if multiple arguments
 * need size checking is (ignore the [0] and [!0] for now):
 * 
 * [!0]
 * start_sizechecking()   [0]
 * sizeinit(arg1)...[!0]...sizecheck_msg [0]
 * sizeinit(arg1)...[!0]...sizecheck_msg [0]
 *   ...
 * sizeinit(arg1)...[!0]...sizecheck_msg [0]
 */
 
/* static vars to maintain state across calls */
/* p_testm is init'ed to non-NULL so the programmer is forced 
   to call start_sizechecking() at outset.  otherwise, coders would be
   liable to forget the sizeinit() opening call.  the null value acts as
   a flag to ensure that sizeinit/sizecheck calls are properly interlaced. 
   The bracketed numbers above ([0]) indicate when this variable is NULL
   and non-NULL during the sequence of calls.  as for the casts, under
   mex, the size of mxArray may not be defined, making it unable to 
   do the pointer arithmetic, so we do the arithmetic as a char.  all we
   need is a non-null value.  */
static const mxArray *p_testm = (const mxArray *)((char *)NULL+1); 
static int test_status;  /* have we found a size match yet */

/* clear the flag: 
 * needed because p_testm is preserved across mex invocations (!) */
int start_sizechecking()
{
  p_testm = NULL; /* there is no "current matrix" */
  return 1;
}

int sizeinit(const mxArray *pm)
{
  if (p_testm == NULL) {  /* calls interlaced properly: no current array */
    p_testm = pm;         /* set up current array */
    test_status = 0;      /* no size agreement yet */
  } else {
    /* we were in the middle of checking another arg! */
    mexErrMsgTxt("internal error: size checking calls in wrong order");
    return 0;
  }
  return 1;
}

/* deprecated: provides number, not name, of arg in error msg */
int sizecheck(const char *msg, int num)
{
  return sizecheck_msg(msg, NULL, num); /* NULL name handled below */
}
	
/* preferred gateway: provides name of arg in error msg */
int sizecheck_msg(const char *msg, const char **argnames, int num)
{
  char errbuf[2*ERRBUFSIZE], sbuf[ERRBUFSIZE];
   
  if (!test_status)  {
    int dim, d = mxGetNumberOfDimensions(p_testm);
    const mwSize *dnums = mxGetDimensions(p_testm);
    for (*sbuf = '\0', dim = 0; dim < d; dim++) {
      sprintf(errbuf, "%d", (int) dnums[dim]);
      strcat(sbuf, errbuf);
      /* printf("%d: %s\n", dim, sbuf); */
      if (dim < d-1) strcat(sbuf, "x"); /* as in, 4x3x2x5 */
    }
    sprintf(errbuf, "%s: arg %s has wrong size (is %s)",
	    msg, 
	    argnames ? argnames[num] : "", 
	    sbuf);
     mexErrMsgTxt(errbuf);
   }
   p_testm = NULL; /* reset to verify interlaced calls */
   return 1;
}
	
/* utility: verify that the "current matrix" is set up */
static int size_null_matrix_check()
{
  if (p_testm == NULL) {
    mexErrMsgTxt("internal error: size checking calls in wrong order");
    return 0;
  }
  return 1;
}

/* 
 * sizeagree:
 * sizes of two quantities match
 */

int sizeagree(const mxArray *pm) /* baseline size agrees with input */
{
  int d;

  size_null_matrix_check();
  d = mxGetNumberOfDimensions(pm); /* do this after above check! */
  if (d != mxGetNumberOfDimensions(p_testm)) 
    return 0; /* unequal sizes */
  test_status |= 
    (memcmp((void *) mxGetDimensions(p_testm), 
	    (void *) mxGetDimensions(pm), 
	    d*sizeof(int)) == 0); /* memcmp == 0 on match */
  return 1;
}

int sizeagreeM(const mxArray *pm) /* baseline size `M' agrees with input */
{
  size_null_matrix_check();
  /* no dimensionality check for agreement along 1 dim */
  test_status |= (mxGetM(p_testm) == mxGetM(pm));
  return 1;
}

int sizeagreeN(const mxArray *pm) /* baseline size `N' agrees with input */
{
  size_null_matrix_check();
  /* no dimensionality check for agreement along 1 dim */
  test_status |= (mxGetN(p_testm) == mxGetN(pm));
  return 1;
}
      
/* note: this only makes sense for 2-d inputs */
int sizeagreeMN(const mxArray *pm) /* baseline size agrees with input or input' */
{
  if (mxGetNumberOfDimensions(p_testm) != 2) return 0; /* insist on 2 dims */
  sizeagree(pm);
  test_status |= ((mxGetN(p_testm) == mxGetM(pm)) && 
		  (mxGetM(p_testm) == mxGetN(pm)));
  return 1;
}

/* 
 * size{is,are}:
 * size matches explicitly given dimensions 
 */

/* made for arbitrary dimensional inputs */
int sizesare(int *s, int d) /* baseline size equals s[0],...s[d-1] */
{
  size_null_matrix_check();
  if (d != mxGetNumberOfDimensions(p_testm)) return 0; /* unequal sizes */
  test_status |= 
    (memcmp((void *) mxGetDimensions(p_testm), 
	    (void *) s, 
	    d*sizeof(int)) == 0); /* memcmp == 0 on match */
  return 1;
}

/* shortcut for 2-d inputs */
int sizeis(int m, int n) /* baseline size equals [m,n] */
{
  size_null_matrix_check();
  if (mxGetNumberOfDimensions(p_testm) != 2) return 0; /* unequal sizes */
  test_status |= ((mxGetM(p_testm) == m) && 
		  (mxGetN(p_testm) == n));
  return 1;
}

/* shortcut for 3-d inputs */
int sizeis3(int m, int n, int p) /* baseline size equals [m,n,p] */
{
  int s[3]; /* hold the sizes */

  size_null_matrix_check();
  if (mxGetNumberOfDimensions(p_testm) != 3) return 0; /* unequal sizes */
  /* note, mxGetN on 3-d inputs is ambiguous */
  s[0] = m; s[1] = n; s[2] = p;
  test_status |= 
    (memcmp((void *) mxGetDimensions(p_testm), 
	    (void *) s, 
	    3*sizeof(int)) == 0); /* memcmp == 0 on match */
  return 1;
}

int sizeisM(int m) /* first baseline size equals m */
{
  size_null_matrix_check();
  test_status |= (mxGetM(p_testm) == m);
  return 1;
}

int sizeisN(int n) /* second baseline size equals n */
{
  size_null_matrix_check();
  test_status |= (mxGetN(p_testm) == n);
  return 1;
}

/* for 2d inputs only */
int sizeisMN(int m, int n) /* baseline size equals [m,n] or [n,m] */
{
  sizeis(m, n);
  sizeis(n, m);
  return 1;
}

