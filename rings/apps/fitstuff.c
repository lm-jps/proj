/*
 *  Simple functions for accessing plain FITS files without extensions
 *    The basic structure associated with the FITS file is kept in a
 *    fitsinfo struct defined here. The header information is just kept
 *    in a string that is not parsed beyound the essentials; it is assumed
 *    that keyword information is available elsewhere, e.g. in an RDB file
 *    or in a DRMS
 *
 *  Functions:
 *    int cpfitsinfo()
 *    fitsinfo_t *rsb_open_fits()
 *    void rsb_close_fits()
 *    int rsb_copy_fits_data()
 *    double rsb_getfitskey_float()
 *    char *rsb_getfitskey_str()
 *    int rsb_keyval_from_prime_header()
 *    int rsb_read_prime_header()
 *    int rsb_read_data()
 *    int rsb_read_data_slice()
 *    fitsinfo_t *rsb_read_fits()
 *    int rsb_write_prime_header()
 *    int rsb_write_fits_data()
 *    int rsb_write_data_slice()
 *    int rsb_write_fits ()
 *
 *  All int functions return 0 on success, non-zero on failure; fitsinfo_t
 *    functions return NULL on failure
 *
 *  Bugs:
 *    rsb_write_fits() was quickly hacked and scaling is not supported
 *    comments in header records are ignored
 *    *rsb_getfitskey_str() fails to compress paired quotes
 *
 *  Revision history is at end of file
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <drms_types.h>

typedef struct fits_info {
  FILE *fp;
  int dataunit;
  int mode;
  char *path;
  int bitpix;
  int naxis;
  int *axis;
  int scaling;
  int allvalid;
  int blank;
  double bscale;
  double bzero;
  char *history;
  void *data;
} fitsinfo_t;

int cpfitsinfo (fitsinfo_t *new, fitsinfo_t *old) {
  int n;

  if (!old || !new) {
    fprintf (stderr, "error in cpfitsinfo(): old or new struct undefined\n");
    return 1;
  }
  new->bitpix = old->bitpix;
  new->naxis = old->naxis;
  if (new->naxis) new->axis = (int *)malloc (new->naxis * sizeof (int));
  for (n = 0; n < new->naxis; n++) new->axis[n] = old->axis[n];
  new->scaling = old->scaling;
  new->allvalid = old->allvalid;
  new->blank = old->blank;
  new->bscale = old->bscale;
  new->bzero = old->bzero;
  return 0;
}

fitsinfo_t *rsb_open_fits (char *filename, char *how) {
  fitsinfo_t *finf;
  FILE *fp = fopen (filename, how);

  if (!fp) {
    fprintf (stderr, "error in rsb_open_fits(): could not open %s for %s\n",
	filename, how);
    return NULL;
  }
  finf = (fitsinfo_t *)malloc (sizeof (fitsinfo_t));
  if (!finf) {
    fprintf (stderr, "error in rsb_open_fits(): could not malloc info struct\n");
    return NULL;
  }
  finf->path = (char *)malloc (strlen (filename) + 1);
  if (!finf->path) {
    fprintf (stderr, "error in rsb_open_fits(): malloc error\n");
    free (finf);
    return NULL;
  }
  finf->fp = fp;
  strcpy (finf->path, filename);
  finf->mode = (strncmp (how, "r", 1)) ? 1 : 0;
  finf->axis = NULL;
  finf->history = NULL;
  finf->data = NULL;
  return finf;
}

void rsb_close_fits (fitsinfo_t *finf) {
  if (finf->mode) {
		/*  make sure writeable file is NULL padded to n*2880 bytes  */
    ;
  }
  if (fclose (finf->fp))
    fprintf (stderr, "Error in rsb_close_fits(): could not close file\n");
  free (finf->path);
  if (finf->axis) free (finf->axis);
  if (finf->history) free (finf->history);
  if (finf->data) free (finf->data);
  free (finf);
}

int rsb_read_prime_header (fitsinfo_t *finf) {
  FILE *fp = finf->fp;
  double dval;
  int card, n;
  int ival;
  int END_not_found = 1, block_ct = 0;
  char *line;
  char fblk[2880], keyname[10];

  if (!fp) {
    fprintf (stderr, "error in rsb_read_prime_header(): file not open\n");
    return 1;
  }
  rewind (fp);
  finf->scaling = 0;
  finf->bscale = 1.0;
  finf->bzero = 0.0;
  finf->allvalid = 1;
  finf->history = NULL;
  finf->history = calloc (2880, sizeof (char));
  while (END_not_found) {
    if (fread (fblk, sizeof (char), 2880, fp) != 2880) {
      fprintf (stderr, "error in rsb_read_prime_header(): read block failed\n");
      return 1;
    }
    card = 0;
    line = fblk;
    if (!block_ct) {
      if (strncmp (line, "SIMPLE  = ", 10)) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	fprintf (stderr, "%80.80s\n", line);
	return 1;
      }
      line += 80;
      card++;
      if (strncmp (line, "BITPIX  = ", 10)) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	fprintf (stderr, "%80s\n", line);
	return 1;
      }
      if (sscanf (line, "BITPIX  = %d", &ival) != 1) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	fprintf (stderr, "%80s\n", line);
	return 1;
      }
      if (ival != 8 && ival != 16 && ival != 32 && ival != 64 &&
          ival != -32 && ival != -64) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	fprintf (stderr, "%80s\n", line);
	return 1;
      }
      finf->bitpix = ival;  
      line += 80;
      card++;
      if (strncmp (line, "NAXIS   = ", 10)) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	fprintf (stderr, "%80s\n", line);
	return 1;
      }
      if (sscanf (line, "NAXIS   = %d", &ival) != 1) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	fprintf (stderr, "%80s\n", line);
	return 1;
      }
      if (ival < 0 || ival > 999) {
	fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	return 1;
      }
      finf->naxis = ival;
      if (finf->naxis) finf->axis = (int *)malloc (finf->naxis * sizeof (int));
      line += 80;
      card++;
      for (n = 0; n < finf->naxis; n++) {
	sprintf (keyname, "NAXIS%d", n + 1);
	if (strncmp (line, keyname, strlen (keyname))) {
	  fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	  fprintf (stderr, "%80s\n", line);
	  return 1;
	}
	if (sscanf (line, "%s = %d", keyname, &ival) != 2) {
	  fprintf (stderr, "FITS conformance error detected in rsb_read_prime_header()\n");
	  fprintf (stderr, "%80s\n", line);
	  return 1;
	}
	finf->axis[n] = ival;
        line += 80;
        card++;
	if (card == 36) {
	  if (fread (fblk, sizeof (char), 2880, fp) != 2880) {
	    fprintf (stderr,
		"error in rsb_read_prime_header(): read block failed\n");
	    return 1;
	  }
	  line = fblk;
	  card = 0;
	  block_ct++;
	}
      }
    }
    block_ct++;
    finf->history = realloc (finf->history, block_ct * 2880);
    for (; card < 36; card++, line += 80) {
      if (!strncmp (line, "BLANK   = ", 10))
        if (sscanf (line, "BLANK   = %d", &ival) == 1)
	  if (finf->bitpix > 0) {
	    finf->blank = ival;
	    finf->allvalid = 0;
	  }
      if (!strncmp (line, "BSCALE  = ", 10))
        if (sscanf (line, "BSCALE  = %lg", &dval) == 1)
	  if (finf->bitpix > 0 && dval != 0.0) {
	    finf->bscale = dval;
	    finf->scaling = finf->bitpix;
	  }
      if (!strncmp (line, "BZERO   = ", 10))
        if (sscanf (line, "BZERO   = %lg", &dval) == 1)
	  if (finf->bitpix > 0) {
	    finf->bzero = dval;
	    finf->scaling = finf->bitpix;
	  }
      if (!strncmp (line, "HISTORY ", 8)) {
	strncat (finf->history, line + 8, 72);
	strcat (finf->history, "\n");
      }
      if (!strncmp (line,
"END                                                                             ",
	  80)) {
	END_not_found = 0;
					     /*  ignore anything after this  */
	card = 36;
      }
    }
  }
  finf->dataunit = ftell (fp);
  return 0;
}

double rsb_getfitskey_float (char *keyvalstr, int *status) {
  double dval;
  char *endptr;

  if (!keyvalstr) {
    if (status) *status = -1;
    return (strtod ("NaN", NULL));
  }
  dval = strtod (keyvalstr, &endptr);
  if (status) *status = 0;
  if (endptr) {
    if (*endptr == '/') return dval;
    if (isblank (*endptr)) return dval;
    if (status) *status = 1;
  }
  return dval;
}
				   /*  bug: fails to compress paired quotes  */
char *rsb_getfitskey_str (char *keyvalstr, int *status) {
  char *c1, *c2, *line, *str;

  if (!keyvalstr) {
    if (status) *status = -1;
    return (NULL);
  }
  line = strdup (keyvalstr);
  if ((c1 = strchr (line, '\'')) && (c2 = strrchr (c1 + 1, '\'')))
    c1++;
  else {			  /*  string but no quotes, take first word  */
						    /*  skip leading blanks  */
    for (c1 = line; *c1 == ' '; c1++) ;
					      /*  get trailing blank or EOS  */
    for (c2=c1; *c2 && *c2 != ' '; c2++) ;
  }
  *c2 = '\0';
  while (c2 > c1 && *(c2-1) == ' ') *--c2 = '\0';
  str = malloc (c2 - c1 + 2);
  strcpy (str, c1);
  free (line);
  if (status) *status = 0;
  return str;
}

int rsb_keyval_from_prime_header (fitsinfo_t *finf, const char *key,
    char *valstr) {
  FILE *fp = finf->fp;
  int card;
  char *line;
  char fblk[2880], search[11];
  int END_not_found = 1, key_found = 0;

  if (!fp) {
    fprintf (stderr, "error in rsb_keyval_from_prime_header(): file not open\n");
    return 1;
  }
  rewind (fp);
  sprintf (search, "%-8.8s= ", key);
  if (fread (fblk, sizeof (char), 2880, fp) != 2880) {
    fprintf (stderr, "error in rsb_keyval_from_prime_header(): read block failed\n");
    return 1;
  }
  card = 0;
  line = fblk;
  while (END_not_found) {
    for (; card < 36; card++, line += 80) {
      if (!strncmp (line, search, 10)) {
	key_found++;
	strncpy (valstr, line+10, 70); 
      }
      if (!strncmp (line,
"END                                                                             ",
 	  80)) {
	END_not_found = 0;
					     /*  ignore anything after this  */
	card = 36;
      }
      if (card == 35) {
	if (fread (fblk, sizeof (char), 2880, fp) != 2880) {
	  fprintf (stderr,
		"error in rsb_keyval_from_prime_header(): read block failed\n");
	  return 1;
	}
	line = fblk;
	card = 0;
      }
    }
  }
  return (key_found - 1);
}

int rsb_write_prime_header (fitsinfo_t *finf) {
  FILE *fp = finf->fp;
  int n, cardct = 0;
  char card[80], ifield[21];

  rewind (fp);
  fprintf (fp, "%-80.80s", "SIMPLE  =                    T");
  cardct++;
  snprintf (card, 80, "BITPIX  = %20d", finf->bitpix);
  fprintf (fp, "%-80.80s", card);
  cardct++;
  snprintf (card, 80, "NAXIS   = %20d", finf->naxis);
  fprintf (fp, "%-80.80s", card);
  cardct++;
  for (n = 0; n < finf->naxis; n++) {
    snprintf (card, 80, "NAXIS%d", n + 1);
    if (n < 99) strcat (card, " ");
    if (n < 9) strcat (card, " ");
    strcat (card, "= ");
    sprintf (ifield, "%20d", finf->axis[n]);
    strcat (card, ifield);
    fprintf (fp, "%-80.80s", card);
    cardct++;
  }
  if (!finf->allvalid) {
    snprintf (card, 80, "BLANK   = %20d", finf->blank);
    fprintf (fp, "%-80.80s", card);
    cardct++;
  }
  if (finf->scaling) {
    snprintf (card, 80, "BSCALE  = %20.12E", finf->bscale);
    fprintf (fp, "%-80.80s", card);
    cardct++;
    snprintf (card, 80, "BZERO   = %20.12E", finf->bzero);
    fprintf (fp, "%-80.80s", card);
    cardct++;
  }
  fprintf (fp, "%-80.80s", "END");
  cardct++;
  snprintf (card, 80, "%-80.80s", " ");
  for (n = 36; n > cardct; n--) fwrite (card, 80, 1, fp);

  finf->dataunit = ftell (fp);
  return 0;
}

int rsb_copy_fits_data (fitsinfo_t *new, fitsinfo_t *old) {
  long long bytecount;
  int blocks, n, status;
  unsigned char block[2880];

  fseek (old->fp, old->dataunit, SEEK_SET);
  fseek (new->fp, new->dataunit, SEEK_SET);
  bytecount = fabs (old->bitpix);
  bytecount /= 8;
  for (n = 0; n < old->naxis; n++) bytecount *= old->axis[n];
  blocks = bytecount / 2880;
  if (bytecount % 2880) blocks++;
  for (n = 0; n < blocks; n++) {
    if ((status = fread (block, 1, 2880, old->fp)) != 2880) {
      fprintf (stderr, "error in rsb_copy_fits_data() reading data block %d\n", n);
      fprintf (stderr, "fread returned %d != 2880\n", status);
      return 1;
    }
    if ((status = fwrite (block, 1, 2880, new->fp)) != 2880) {
      fprintf (stderr, "error in rsb_copy_fits_data() writing data block %d\n", n);
      fprintf (stderr, "fwrite returned %d != 2880\n", status);
      return 1;
    }
  }

  return 0;
}

int rsb_copy_fits_data_slices (fitsinfo_t *new, fitsinfo_t *old,
    int slice_size, long long start, int slice_count) {
  long long bytecount;
  int n, status;
  unsigned char *block;

  fseek (old->fp, old->dataunit, SEEK_SET);
  fseek (new->fp, new->dataunit, SEEK_SET);
  fseek (old->fp, slice_size * start, SEEK_CUR);
  block = malloc (slice_size);
  bytecount = 0;
  for (n = 0; n < slice_count; n++) {
    if ((status = fread (block, 1, slice_size, old->fp)) != slice_size) {
      fprintf (stderr, "error in rsb_copy_fits_data_slices() reading data slice %d\n", n);
      fprintf (stderr, "fread returned %d != %d\n", status, slice_size);
      return 1;
    }
    if ((status = fwrite (block, 1, slice_size, new->fp)) != slice_size) {
      fprintf (stderr, "error in rsb_copy_fits_data_slices() writing data block %d\n", n);
      fprintf (stderr, "fwrite returned %d != %d\n", status, slice_size);
      return 1;
    }
    bytecount += slice_size;
  }
  if (bytecount % 2880)
    for (n = bytecount % 2880; n < 2880; n++) fputc (0, new->fp);
  free (block);

  return 0;
}
							/*  not implemented  */
int rsb_copy_fits_data_slices_binned (fitsinfo_t *new, fitsinfo_t *old,
    int slice_size, long long start, int slice_count, int binwidth) {
  long long bytecount;
  int n;
  int binned_slice = slice_size / binwidth / binwidth;

  fseek (old->fp, old->dataunit, SEEK_SET);
  fseek (new->fp, new->dataunit, SEEK_SET);
  fseek (old->fp, slice_size * start, SEEK_CUR);
  bytecount = 0;
  for (n = 0; n < slice_count; n++) {
    bytecount += binned_slice;
  }
  if (bytecount % 2880)
    for (n = bytecount % 2880; n < 2880; n++) fputc (0, new->fp);
  return 0;
}

static void flip_byte_order (void *data, long long vlength, int dlength) {
  unsigned long long *d;
  unsigned int *i;
  unsigned short *s;
  long long n;

  if (!data) return;
  switch (dlength) {
    case (8):
      d = (unsigned long long *)data;
      for (n = 0; n < vlength; n++, i++) *d = (*d << 32) | (*d >> 32);
      vlength *= 2;
    case (4):
      i = (unsigned int *)data;
      for (n = 0; n < vlength; n++, i++) *i = (*i << 16) | (*i >> 16);
      vlength *= 2;
    case (2):
      s = (unsigned short *)data;
      for (n = 0; n < vlength; n++, s++) *s = (unsigned short)(*s << 8) | (unsigned short)(*s >> 8);
  }
  return;
}

int rsb_read_fits_data (fitsinfo_t *fits, void **data) {
  long long bytecount, ntot, brem;
  int dsize, n, status, vlen;
  char *ndat, *buf;
  int buflen = 16777216;

  fseek (fits->fp, fits->dataunit, SEEK_SET);
  bytecount = fabs (fits->bitpix);
  bytecount /= 8;
  dsize = bytecount;
  ntot = (fits->naxis) ? 1 : 0;
  for (n = 0; n < fits->naxis; n++) ntot *= fits->axis[n];
  bytecount *= ntot;

  buf = *data = (char *)malloc (bytecount);

#if __BYTE_ORDER == __LITTLE_ENDIAN
	/*  To avoid having to malloc twice the size of the data array,
					just read 16 MB chunks to swab  */
  ndat = malloc (buflen);
  brem = bytecount;
  vlen = buflen / dsize;
  while (brem > buflen) {
    if ((status = fread (ndat, 1, buflen, fits->fp)) != buflen) {
      fprintf (stderr, "error in rsb_read_fits_data()\n");
      fprintf (stderr, "fread returned %d != %d\n", status, buflen);
      return 1;
    }
    flip_byte_order (ndat, vlen, dsize);
    memcpy (buf, ndat, buflen);
    brem -= buflen;
    buf += buflen;
  }
  if (brem) {
    if ((status = fread (ndat, 1, brem, fits->fp)) != brem) {
      fprintf (stderr, "error in rsb_read_fits_data()\n");
      fprintf (stderr, "fread returned %d != %lld\n", status, brem);
      return 1;
    }
    vlen = brem / dsize;
    flip_byte_order (ndat, vlen, dsize);
    memcpy (buf, ndat, brem);
  }
  free (ndat);
#else
  if ((status = fread (*data, 1, bytecount, fits->fp)) != bytecount) {
    fprintf (stderr, "error in rsb_read_fits_data()\n");
    fprintf (stderr, "fread returned %d != %ld\n", status, bytecount);
    return 1;
  }
#endif
  return 0;
}

fitsinfo_t *rsb_read_fits (char *filename) {
  void *data;
  fitsinfo_t *finf = rsb_open_fits (filename, "r");
  if (!finf) return finf;
  if (rsb_read_prime_header (finf)) return NULL;
  if (rsb_read_fits_data (finf, &data)) return NULL;
  finf->data = data;
  return finf;
}

int rsb_read_data_slice (fitsinfo_t *fits, void **data, int *offsets,
    int *lengths) {
  long long bytecount, ntot, offset, brem;
  int dsize, n, status, vlen;
  char *ndat, *buf;
  int buflen = 16777216;

  bytecount = fabs (fits->bitpix);
  bytecount /= 8;
  dsize = bytecount;
  ntot = (fits->naxis) ? 1 : 0;
  for (n = 0; n < fits->naxis; n++) {
    if (lengths[n] <= 0) lengths[n] = fits->axis[n];
    ntot *= lengths[n];
  }
  bytecount *= ntot;
  buf = *data = (char *)malloc (bytecount);

  offset = 0;
  ntot = (fits->naxis) ? 1 : 0;
  for (n = 0; n < fits->naxis; n++) {
    offset += ntot * offsets[n];
    ntot *= fits->axis[n];
  }
  offset *= dsize;

  fseek (fits->fp, fits->dataunit, SEEK_SET);
  fseek (fits->fp, offset, SEEK_CUR);

#if __BYTE_ORDER == __LITTLE_ENDIAN
	/*  To avoid having to malloc twice the size of the data array,
					just read 16 MB chunks to swab  */
  ndat = malloc (buflen);
  brem = bytecount;
  vlen = buflen / dsize;
  while (brem > buflen) {
    if ((status = fread (ndat, 1, buflen, fits->fp)) != buflen) {
      fprintf (stderr, "error in rsb_read_fits_data()\n");
      fprintf (stderr, "fread returned %d != %d\n", status, buflen);
      return 1;
    }
    flip_byte_order (ndat, vlen, dsize);
    memcpy (buf, ndat, buflen);
    brem -= buflen;
    buf += buflen;
  }
  if (brem) {
    if ((status = fread (ndat, 1, brem, fits->fp)) != brem) {
      fprintf (stderr, "error in rsb_read_fits_data()\n");
      fprintf (stderr, "fread returned %d != %lld\n", status, brem);
      return 1;
    }
    vlen = brem / dsize;
    flip_byte_order (ndat, vlen, dsize);
    memcpy (buf, ndat, brem);
  }
  free (ndat);
#else
  if ((status = fread (*data, 1, bytecount, fits->fp)) != bytecount) {
    fprintf (stderr, "error in rsb_read_fits_data()\n");
    fprintf (stderr, "fread returned %d != %ld\n", status, bytecount);
    return 1;
  }
#endif
  return 0;
}

int rsb_write_fits_data (fitsinfo_t *new, void *data) {
  long long bytecount, ntot, brem, status;
  int dsize, n, vlen;
  char *ndat, *buf;
  int buflen = 16777216;

  fseek (new->fp, new->dataunit, SEEK_SET);
  bytecount = fabs (new->bitpix);
  bytecount /= 8;
  dsize = bytecount;
  ntot = (new->naxis) ? 1 : 0;
  for (n = 0; n < new->naxis; n++) ntot *= new->axis[n];
  bytecount *= ntot;
#if __BYTE_ORDER == __LITTLE_ENDIAN
		/*  To avoid having to malloc an additional data array,
					just write 16 MB chunks to swab  */
  ndat = malloc (buflen);
  brem = bytecount;
  vlen = buflen / dsize;
  buf = (char *)data;
  while (brem > buflen) {
    memcpy (ndat, buf, buflen);
    flip_byte_order (ndat, vlen, dsize);
    if ((status = fwrite (ndat, 1, buflen, new->fp)) != buflen) {
      fprintf (stderr, "error in rsb_write_fits_data()\n");
      fprintf (stderr, "fwrite returned %lld != %d\n", status, buflen);
      return 1;
    }
    brem -= buflen;
    buf += buflen;
  }
  if (brem) {
    vlen = brem / dsize;
    memcpy (ndat, buf, brem);
    flip_byte_order (ndat, vlen, dsize);
    if ((status = fwrite (ndat, 1, brem, new->fp)) != brem) {
      fprintf (stderr, "error in rsb_write_fits_data()\n");
      fprintf (stderr, "fwrite returned %lld != %lld\n", status, brem);
      return 1;
    }
  }
  free (ndat);
#else
  if ((status = fwrite (data, 1, bytecount, new->fp)) != bytecount) {
    fprintf (stderr, "error in rsb_write_fits_data()\n");
    fprintf (stderr, "fwrite returned %ld != %ld\n", status, bytecount);
    return 1;
  }
#endif

  if (bytecount % 2880) {
    for (n = bytecount%2880; n < 2880; n++) fputc (0, new->fp);
  }
  return 0;
}

int rsb_write_fits (DRMS_Array_t *array, char *filename) {
  fitsinfo_t *fits;
/*
  float *odata;
  short *data;
*/
  long long n, nt;

  fits = rsb_open_fits (filename, "w");
  if (!fits) return 1;
  fits->bitpix = 8 * drms_sizeof (array->type);
  if (array->type == DRMS_TYPE_FLOAT || array->type == DRMS_TYPE_DOUBLE)
    fits->bitpix *= -1;
/*
  fits->bitpix = 16;
*/
  fits->naxis = array->naxis;
  fits->axis = (int *)malloc (fits->naxis * sizeof (int));
  nt = 1;
  for (n = 0; n < fits->naxis; n++) {
    fits->axis[n] = array->axis[n];
    nt *= fits->axis[n];
  }
  fits->allvalid = 0;
  fits->blank = -32768;
  fits->scaling = 1;
  fits->bscale = 1.0;
  fits->bzero = 0.0;
  rsb_write_prime_header (fits);
/*
  data = (short *)malloc (nt * sizeof (short));
  odata = (float *)array->data;
  for (n = 0; n < nt; n++) data[n] = odata[n];
  rsb_write_fits_data (fits, data);
*/
  rsb_write_fits_data (fits, array->data);
  rsb_close_fits (fits);
  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  09.09.17		version that went into first release under proj/rings
 *  09.12.02		fixed three icc11 compiler warnings
 */
