/*
 * mex2c_util
 * 
 * generally useful utilities for mex2c, in either its command-line
 * or its xml/pleo incarnations
 *
 * version 3, Michael Turmon, 1999, 2001
 */

/*LINTLIBRARY*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>     /* unlink() */
#include <fitsio.h>
#include "mex2c_util.h" /* function templates: this file */
#include "mex.h"        /* mxArray data structure */
#include "mexhead.h"    /* standard mex utilities */

/* global variables to hold program state */
extern char *mex2c_progname;  /* name of this program, for error messages */
extern char mex2c_phase[];     /* for error messages: where are we */

/************************************************************************
 * 
 * Miscellaneous
 * 
 ***********************************************************************/

/* 
 * Print out cfitsio error messages, if any, and exit program 
 */
void
printerror(int status)
{
  fprintf(stderr, "%s: Runtime error (%s).  Exiting.\n",
	mex2c_progname, 
	status == 0 ? "not a FITSIO error" : "FITSIO trace follows");
  if (*mex2c_phase)
    fprintf(stderr, "%s: error while: %s\n", 
	    mex2c_progname, 
	    mex2c_phase);
  fits_report_error(stderr, status); /* print whole message stack */
  exit(status != 0 ? status : 1);  /* terminate, ensure nonzero status */
}


/*
 * compute the range of supplied input
 */
void
compute_range(double *data, long n, double *dmin, double *dmax)
{
  long i;
  double min1, max1; /* local copies for faster adjustment */
  
  max1 = -1.0 * mxt_getinfd();
  min1 =        mxt_getinfd();
  for (i = 0; i < n; i++, data++) {
    if (!isnan(*data)) {
      if (*data < min1) min1 = *data;
      if (*data > max1) max1 = *data;
    }
  }
  /* if there were no valid entries, this is OK - leaving min > max
   * will turn scaling off according to our conventions.  if there
   * was just one valid entry, min = max which is also OK.  finally,
   * if either min=-Inf or max=+Inf, that's the user's problem! */
  /* set output pointers */
  *dmin = min1; *dmax = max1;
}


/************************************************************************
 * 
 * FITS I/O
 * 
 ***********************************************************************/

/*
 * get a fits file to a mxArray structure
 */
mxArray *
GetFITS(char *filename, int file_transpose)
{
  mxArray *pm;          /* returned mxArray */
  fitsfile *fptr;       /* fitsio input file */
  int status = 0;       /* fitsio status: OK to start with */
  int naxis;            /* n-dimensional image */     
  mwSize *naxesI;       /* extent along each dimension (mwSize) */
  long *naxesL;         /* extent along each dimension (long) */
  int nfound;           /* returned count */
  long nelements;       /* total number of elements */
  int i;                /* count dimensions */
  int    bitpix;
  double bscale, bzero;
  double nullval;       /* cfitsio BLANK checking */
  double datamin, datamax;

  fits_open_file(&fptr, filename, READONLY, &status);  /* open for reading */
  /* read the NAXIS keyword to get number of dimensions */
  fits_read_key(fptr, TINT, "NAXIS", &naxis, NULL, &status);
  if (status)  /* bail before the allocation */
    printerror(status);
  /* make space for the dimensions */
  naxesI = (mwSize *) calloc((size_t) naxis, sizeof(*naxesI));
  naxesL = (long *)   calloc((size_t) naxis, sizeof(*naxesL));
  if (naxesI == NULL || naxesL == NULL) { /* highly unlikely */
    fprintf(stderr, "Failed to input-convert matrix (failed calloc)\n");
    printerror(0);
  }
  /* read the NAXIS* keywords to get image size */
  fits_read_keys_lng(fptr, "NAXIS", 1, naxis, naxesL, &nfound, &status);
  for (i = 0; i < naxis; i++) 
    naxesI[i] = (mwSize) naxesL[i];
  if (status)  /* bail before the allocation */
    printerror(status);
  /* read in bitpix (mandatory KW) */
  fits_read_key(fptr, TLONG, "BITPIX", &bitpix, NULL, &status);
  if (status)
    printerror(status);
  /* set up nullval, which controls cfitsio NaN checking, based on bitpix */
  if (bitpix < 0)
    /* floating point: tell cfitsio not to check if pixel[i] is a
       special value -- just pass the values straight through.  
       If cfitsio checking is turned on, it converts all 
       IEEE special values, both NaN and +/- Inf, to nullval */
    nullval = 0.0; /* i.e., turn checking off and pass NaN/Inf through */
  else
    /* integer: let cfitsio check if pixel[i] is BLANK.
       When cfitsio checking is turned on, it converts all 
       BLANKs to nullval (NaN). */
    nullval = mxt_getnand();  /* set FITS BLANKs to this */
  /* create matrix */
  pm = mxCreateNumericArray((mwSize) naxis, naxesI, mxDOUBLE_CLASS, mxREAL);
  nelements = (long) mxGetNumberOfElements(pm);
  /* read in the image */
  fits_read_img(fptr, TDOUBLE, 1L, nelements, &nullval,
		mxGetPr(pm), &nfound, &status);
  /* read in bscale and bzero (optional KWs) */
  if (fits_read_key(fptr, TDOUBLE, "BSCALE", &bscale, NULL, &status)) {
    bscale = 1.0; status = 0; fits_clear_errmsg(); 
  }
  if (fits_read_key(fptr, TDOUBLE, "BZERO", &bzero, NULL, &status)) {
    bzero = 0.0;  status = 0; fits_clear_errmsg(); 
  }
  /* close the file */
  if (fits_close_file(fptr, &status))
    printerror(status);

  /* convert to matlab-style ordering if needed */
  if (file_transpose && !mxt_transpose_double_mxArray(pm)) {
    fprintf(stderr, "Failed to input-convert matrix (failed malloc)\n");
    printerror(0);
  }
  /* initialize range of data */
  switch (bitpix) {
  case 8:
    datamin = bscale * (  0.0) + bzero;
    datamax = bscale * (254.0) + bzero;
    break;
  case 16:
    /* signed 16-bit integers, (+/-) (2^15 - 1) */
    datamin = bscale * (-32767.0) + bzero;
    datamax = bscale * ( 32767.0) + bzero;
    break;
  case 32: 
    /* signed 32-bit integers, (+/-) (2^31 - 1) */
    datamin = bscale * (-2147483647.0) + bzero;
    datamax = bscale * ( 2147483647.0) + bzero;
    break;
  case -32:
  case -64:
    /* indicate the range is unknown -- note, if output is also
     * floating point, the range will be unneeded anyway. */
    datamin = datamax = mxt_getnand();
    break;
  default:
    datamin = datamax = mxt_getnand(); // kill compiler warning
    fprintf(stderr, "getfits: illegal bitpix = %d", bitpix);
    printerror(0);
    break;
  }
  setrange(pm, datamin, datamax);
  return(pm);  /* return the transposed matrix */
}

/*
 * put a mxArray structure to a fits file 
 */
void
PutFITS(mxArray *pm, char *filename, int bitpix, int file_transpose, double scaling)
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status;           /* cfitsio error status */
  int naxis;            /* n-dimensional image */     
  const mwSize *naxesI; /* extent along each dimension (mwSize)  */
  long *naxesL;         /* extent along each dimension (long) */
  long nelements;       /* total number of elements */
  long i;                /* count dimensions */
  double datamin, datamax;
  double min, max;
  double bscale, bzero;
  int blank;
  /* const double nan = mxt_getnand(); */

  /* transpose the mxArray if needed */
  if (file_transpose && !mxt_transpose_double_mxArray(pm)) {
    fprintf(stderr, "Failed to output-convert matrix (failed malloc)\n");
    printerror(0);
  }
  /* file creation */
  unlink(filename);   /* delete old file */
  status = 0;         /* initialize status before calling fitsio routines */
  if (fits_create_file(&fptr, filename, &status)) 
    printerror(status);
  /* initialize the header */
  naxis = (int) mxGetNumberOfDimensions(pm);
  nelements = (long) mxGetNumberOfElements(pm);
  naxesI = mxGetDimensions(pm); /* mwSize version */
  naxesL = (long *) calloc(naxis, sizeof(long));
  if (naxesL == NULL) { /* highly unlikely */
    fprintf(stderr, "Failed to output-convert matrix (failed calloc)\n");
    printerror(0);
  }
  for (i = 0; i < naxis; i++) 
    naxesL[i] = (long) naxesI[i];
  if (fits_create_img(fptr, bitpix, naxis, naxesL, &status))
    printerror(status); 
  free(naxesL);
  /* 
   * Determine scaling for image data
   */
  /* read in datamin and datamax */
  getrange(pm, &datamin, &datamax);
  if ((isnan(datamin) || isnan(datamax)) && (bitpix > 0)) {
    /* don't have range and will need it: loop again to find */
    compute_range(mxGetPr(pm), nelements, &datamin, &datamax);
  }

  /* compute blank, bscale and bzero for bitpix > 0 */
  /* Allow guidance from the caller:
   *   if scaling > 0, let bzero = 0, bscale = scaling.
   *   if scaling == 0, define scaling to be the identity.  
   * Otherwise (scaling < 0), use the stated range:
   *   if datamin < datamax, then define the scaling to be the identity;
   *   if datamin == datamax, let bzero = datamax and scale = 1.0 */
  switch (bitpix) {
  case 8:
    min = 0;
    max = 254;
    blank = 255;
    if (scaling > 0) {
      bzero = 0; bscale = scaling;
    } else if (scaling == 0) {
      bzero = 0; bscale = 1;
    } else if (datamax > datamin) { /* typical case */
      bscale = (0.5 * datamax - 0.5 * datamin) / (0.5 * max - 0.5 * min);
      bzero = datamin - min * bscale;
    } else if (datamax == datamin) { /* e.g., all numbers equal */
      bscale = 1; bzero = datamax; 
    } else { /* pass values through directly */
      bscale = 1; bzero = 0; 
    }
    break;
  case 16:
    max = 32767.0; min = -max; /* (+/-) (2^15 - 1) */
    blank = min-1;
    if (scaling > 0) {
      bzero = 0; bscale = scaling;
    } else if (scaling == 0) {
      bzero = 0; bscale = 1;
    } else if (datamax > datamin) { /* typical case */
      bscale = (0.5 * datamax - 0.5 * datamin) / (0.5 * max - 0.5 * min);
      bzero = datamin - min * bscale;
    } else if (datamax == datamin) { /* e.g., all numbers equal */
      bscale = 1; bzero = datamax; 
    } else { /* pass values through directly */
      bscale = 1; bzero = 0; 
    }
    break;
  case 32:
    max = 2147483647.0; min = -max; /* (+/-) (2^31 - 1) */
    blank = min-1;
    if (scaling > 0) {
      bzero = 0; bscale = scaling;
    } else if (scaling == 0) {
      bzero = 0; bscale = 1;
    } else if (datamax > datamin) { /* typical case */
      bscale = (0.5 * datamax - 0.5 * datamin) / (0.5 * max - 0.5 * min);
      bzero = datamin - min * bscale;
    } else if (datamax == datamin) { /* e.g., all numbers equal */
      bscale = 1; bzero = datamax; 
    } else { /* pass values through directly */
      bscale = 1; bzero = 0; 
    } 
    break;
  case -32:
  case -64:
    bscale = 1; bzero = 0; /* no need for blank */
    break;
  default:
    fprintf(stderr, "putfits: illegal bitpix = %d", bitpix);
    printerror(0);
    break;
  }
  /* insert bscale, bzero into header if necessary */
  if (bscale != 1 || bzero != 0) {
    fits_write_key(fptr, TDOUBLE, "BSCALE", &bscale, "bscale", &status);
    fits_write_key(fptr, TDOUBLE, "BZERO",  &bzero,  "bzero", &status);
  }
  /* insert blank into header for the integer-valued types */
  if (bitpix > 0) {
    fits_write_key(fptr, TINT,    "BLANK",  &blank,  "blank", &status);
  }
  /* miscellaneous info */
  fits_write_date(fptr, &status);
  fits_write_history(fptr, "Generated by mex2c version 4.0", &status);
  if (status)
    printerror(status);
  /* write the array of doubles to the FITS file */
  /* must manually change NaNs in the data to double-precision
     values which will be translated to BLANK by fits_write.
     I cannot use fits_write_imgnull, which is supposed to handle
     BLANKs.  This is because that routine inserts BLANK when
     the data == a supplied number.  If the supplied number is
     an IEEE NaN, then the == never succeeds, and BLANK is
     never written.  So we must take care of it here. */
  if (bitpix >= 0) {
    double *data = mxGetPr(pm);
    double blank_data = bscale * blank + bzero; /* transforms into BLANK */
    long i;

    for (i = 0; i < nelements; i++)
      if (isnan(data[i]))
	data[i] = blank_data;
  }
  /* write the data */
  fits_write_img(fptr, TDOUBLE, 1L, nelements, mxGetPr(pm), &status);
  /* close the file */
  fits_close_file(fptr, &status);
  /* One last check to be sure things went OK */
  if (status)
    printerror(status);    
}


/************************************************************************
 * 
 * literal matrix I/O
 * 
 ***********************************************************************/

/*
 * replaces the old version
 */
mxArray *
GetNumericMatrix(char *value)
{
  size_t N;      /* number of total data elements */
  mwSize Nd;     /* number of dimensions */
  mwSize *dimsM; /* size of each dimension (alloc'd elsewhere) */
  double *data;  /* data itself */
  mwSize n_data = 16384; /* max #data elements */
  double datamin, datamax;
  mxArray *pm;

  if (!(data = calloc((size_t) n_data, sizeof(*data))))
    return NULL;
  if (!(dimsM = string_to_array(n_data, value, &Nd, data))) {
    free(data);
    return NULL;
  }
  /* make space for data */
  if (!(pm = mxCreateNumericArray((mwSize) Nd, dimsM, mxDOUBLE_CLASS, mxREAL))) {
    free(dimsM);
    free(data);
    return NULL;
  }
  /* copy data over */
  N = (size_t) mxGetNumberOfElements(pm);
  memcpy(mxGetPr(pm), data, (size_t) N*sizeof(double));
  compute_range(data, (long) N, &datamin, &datamax); /* easy to do here */
  setrange(pm, datamin, datamax);
  /* clean up and return */
  free(dimsM);
  free(data);
  return pm;
}


/*
 * Store a string of form: ['s1';'s2';...;'sN'] into
 * an mxArray data structure.  No white space is allowed
 * outside quoted strings.  Component strings may be
 * of differing lengths; they are padded with spaces.
 * Semicolons may be present within strings.  There
 * is no method to insert a literal ' into a string.
 */
mxArray *
GetStringMatrix(char *value)
{
  mxArray *pm; /* the result */
  char **strs; /* the component strings */
  int  Nrow; /* total number of rows */
  int i;  /* counts rows */
  char *s, *s1, *s2; /* pointers within value */
  char *value_end; /* points to terminating null in value */
  int quote_level; /* status during scan of quoted string */

  value_end = value + strlen(value);
  /*
   * 1: Check format, and count rows 
   */
  /* check start and end positions */
  if (*value != '[') {
    fprintf(stderr, "invalid ([ bracket) string matrix %s\n", value);
    printerror(0);
  }
  if (*(value_end-1) != ']') {
    fprintf(stderr, "invalid (] bracket) string matrix %s\n", value);
    printerror(0);
  }
  /* check format and count number of rows; not completely trivial
     because we must allow for ; separators within quoted strings */
  quote_level = 1; /* looking for new quoted string */
  for (Nrow = 1, s = value + 1; s < value_end-1; s++) {
    if (quote_level == 2) {
      /* within quoted string: skip chars until ' */
      if (*s == '\'')
	quote_level = 0; /* out of quoted string mode */
    } else if (quote_level == 0) {
      /* out of quoted string: must see ; separator */
      if (*s != ';') {
	fprintf(stderr, "invalid (missing ;) string matrix %s\n", value);
	printerror(0);
      }
      Nrow++; /* count a new string now */
      quote_level++; /* next stage */
    } else if (quote_level == 1) {
      /* looking for new quoted string */
      if (*s != '\'') {
	fprintf(stderr, "invalid (missing ') string matrix %s\n", value);
	printerror(0);
      }
      quote_level++; /* next stage */
    }
  } /* end for */
  /* ensure correct end state */
  if (quote_level != 0) {
    fprintf(stderr, "invalid (mismatched ') string matrix %s\n", value);
    printerror(0);
  }
  /*
   * 2: Scan in component strings
   */
  /* allocate space for component strings */
  strs = (char **) calloc(Nrow, sizeof(char *));
  if (strs == NULL) { /* unlikely! */
    fprintf(stderr, "failed calloc for string matrix %s\n", value);
    printerror(0);
  }
  /* extract component strings; since format has checked out OK,
     can just look for stuff between ' pairs. */
  s2 = value; /* "end of last found string" -- first char in value */
  for (i = 0; i < Nrow; i++) {
    s1 = strchr(s2+1, '\''); /* find open quote */
    s2 = strchr(s1+1, '\''); /* find close quote */
    /* NB, one extra char below == NULL since calloc zeros its buffer */
    strs[i] = (char *) calloc(s2 - s1, sizeof(char)); 
    if (strs[i] == NULL) { /* unlikely! */
      fprintf(stderr, "failed calloc for string matrix %s\n", value);
      printerror(0);
    }
    memcpy(strs[i], s1+1, (size_t) (s2-(s1+1))); /* s1..s2 excluding ends */
  }
  /*
   * 3: Create object and clean up
   */
  /* cast below from char ** to const char ** to match function template */
  pm = mxCreateCharMatrixFromStrings((mwSize) Nrow, (const char **) strs);
  if (pm == NULL) { /* should not happen */
    fprintf(stderr, "object creation failed for string matrix %s\n", value);
    printerror(0);
  }
  /* set up its range: probably unneeded */
  setrange(pm, 0.0, 254.0);
  /* free up intermediate stuff */
  for (i = 0; i < Nrow; i++) 
    free(strs[i]);
  free(strs);
  /* that's all */
  return(pm);
}


/* print numeric matrix segment to stream: recursive helper function
 *
 * stride is the size of the overall matrix to be printed; the 
 * recursive step divides one matrix of leading dimension D and
 * stride S into D recursive calls, each of stride S/D.
 * The matrix is therefore printed in the order it appears
 * in memory.  The innermost nesting is printed as a 
 * vector of length dims[ndim-1], and the outermost nesting
 * is a group of length dims[0].
 *
 * delim arg is various delimiters, as described below.
 * 
 */
static
void
PutNumericMatrix_level(FILE *fp, 
		       int ndim,
		       const int *dims,
		       int stride,
		       double *pr, 
		       char *delim[])
{
  if (ndim > 2) {
    /* recursive case */
    int i;

    fputs(delim[0], fp);     /* matrix start */
    for (i = 0; i < dims[0]; i++) {
      if (i > 0) fputs(delim[4], fp);   /* between-matrix delimiter */
      PutNumericMatrix_level(fp, ndim-1, dims+1, stride/dims[0], 
			     pr+(stride/dims[0])*i, delim);
    }
    fputs(delim[1], fp);     /* matrix end */
  } else {
    /* base case */
    int i, j;
    double num;

    fputs(delim[0], fp);       /* matrix start */
    for (i = 0; i < dims[0]; i++) {
      if (i > 0) fputs(delim[5], fp);   /* between-row delimiter */
      fputs(delim[2], fp);     /* row start */
      for (j = 0; j < dims[1]; j++) {
	if (j > 0) fputs(delim[6], fp); /* between-element delimiter */
	/* put the number -- checking for special values */
	num = pr[i*dims[1]+j];  /* current number */
	if (isnan(num))
	  fputs("NaN", fp);     /* IEEE nan */
	else if (!finite(num) && num < 0)
	  fputs("-Inf", fp);    /* IEEE -infinity */
	else if (!finite(num) && num > 0)
	  fputs("Inf", fp);     /* IEEE +infinity */
	else
	  fprintf(fp, delim[7], num); /* a number */
      }
      fputs(delim[3], fp);      /* row end */
    }
    fputs(delim[1], fp);        /* matrix end */
  }
}

/* print (any-dimensional) numeric array to stream; need at least two dims.
 *
 * delim arg is various delimiters:
 * matrix start, matrix end,
 * row start, row end,
 * inter-matrix delimiter,
 * inter-row delimiter, inter-element delimiter,
 * element format
 * 
 * the array is printed as a nested series of 2d matrices, separated
 * by inter-matrix delimiter, and surrounded by matrix start/end.
 *
 * typical delimiter values would be:
 * {"[", "]", ",", "", "", ";", ",", "%12.9f"} for mex2c_cli style
 * {"[", "]", " ", "[", "]", "\n", " ", "%12.9f"} for mex2c_xml style
 */

void
PutNumericArray(FILE *fp, 
		mxArray *pa, 
		char *head, char *tail, char *delim[])
{
  int d, Nd;
  const mwSize *dimsM; /* dims within pa */
  int *dimsI;    /* dims as int's for the printing routine */

  /* note, dimsI refers to the original, untransposed, array */
  Nd = (int) mxGetNumberOfDimensions(pa);
  dimsI = (int *) mxMalloc((size_t) (Nd * sizeof(int)));
  dimsM = mxGetDimensions(pa);
  for (d = 0; d < Nd; d++)
    dimsI[d] = (int) dimsM[d];
  /*
  printf("<<put0: nd=%d, d=[%d,%d,..], l=%d a[0]=%x>>", Nd, dimsI[0], dimsI[1],
	 (int) mxGetNumberOfElements(pa), 
	 mxGetPr(pa)); */
  /* easiest to transpose, because elements are printed in this order */
  mxt_transpose_double_mxArray(pa); /* NB: this overwrites pa -> dims */
  fputs(head, fp);           /* starting string */
  /*
  printf("<<put1: nd=%d, d=[%d,%d,..], l=%d a[0]=%x>>", Nd, dimsI[0], dimsI[1],
	 mxGetNumberOfElements(pa), 
	 mxGetPr(pa)); */
  PutNumericMatrix_level(fp, Nd, dimsI, (int) mxGetNumberOfElements(pa), 
			 mxGetPr(pa), delim);
  fputs(tail, fp);           /* closing string */
  mxFree((void *) dimsI);
}


/* print string matrix to stream.
 *
 * delim arg is various delimiters:
 * matrix start, matrix end,
 * string start, string end,
 * inter-string delimiter
 * 
 * typical delimiter values would be:
 * {"[", "]", "\"", "\"", ";"} for mex2c_cli style
 * {"[", "]", "", "", "\n"} for mex2c_xml style
 */
void
PutStringMatrix(FILE *fp,
		mxArray *pm, 
		char *head, char *tail, char *delim[])
{
  int m = mxGetM(pm);
  int n = mxGetN(pm);
  char *pc;
  int i,j;

  pc = mxArrayToString(pm);  /* extract string data as char's */
  /* the function malloc's and so may fail */
  if (pc == NULL) {
    fprintf(stderr, "object creation failed on printing string matrix\n");
    printerror(0);
  }
  fputs(head, fp);           /* starting string */
  fputs(delim[0], fp);       /* matrix start */
  for(i = 0; i < m; i++) {
    fputs(delim[2], fp);     /* string start */
    for(j = 0; j < n; j++) {
      fprintf(fp, "%c", pc[j*m+i]);
    }
    fputs(delim[3], fp);     /* string end */
    if( i != m-1 )
      fputs(delim[4], fp);       /* matrix end */
  }
  fputs(delim[1], fp);       /* matrix end */
  fputs(tail, fp);           /* closing string */
  /* free the memory */
  free(pc);
}

