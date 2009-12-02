/*
 *  keystuff.c					(linked from) ~rick/src/util
 *
 *  library of miscellaneous utility functions for dealing with DRMS record keys
 *
 *  check_and_setkey_TYPE()
 *  append_keyval_to_primekeyval()
 *  create_primekey_from_keylist()
 *  drms_wcs_timestep()
 *
 *  Bugs:
 *    There is evidently a bug in append_keyval_to_primekeyval() in keystuff.c,
 *	as the ingest_track module was generating a seg fault in
 *	drms_setkey_string() with values of the PrimeKeyString like
 *  1987:000.0|22.5|-52.5|2002.03.30_00:44:30_TAI|1664|16.0|16.0|0.100000|0.000000|Postels|0.12500|0.00000|prog:mdi,level:lev1.5,series:fd_V_01h[80986-81014],sel:[53-36]
 *	I was only able to work around this by changing the initial value of
 *	buflen in append_keyval_to_primekeyval() from 128 to 256
 *    Many combinations of legitimate prefixes and unit strings for WCS time
 *	descriptions are not yet supported
 *
 *  Revision history is at end of file
 */

int check_and_set_key_int (DRMS_Record_t *new, const char *key, int val) {
  DRMS_Keyword_t *keywd;
  int vreq;
  int status;

  if (!(keywd = drms_keyword_lookup (new, key, 1))) return 0;
  if (keywd->info->recscope != 1) {
		       /*  it's not a constant, so don't worry, just set it  */
    drms_setkey_int (new, key, val);
    return 0;
  }
  vreq = drms_getkey_int (new, key, &status);
  if (status) {
    fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
    return 1;
  }
  if (vreq != val) {
    char format[256];
    sprintf (format,
	"Error:  parameter value %s for keyword %%s\n  differs from required value %s\n",
	keywd->info->format, keywd->info->format);
    fprintf (stderr, format, val, key, vreq);
    return 1;
  }
  return 0;
}

int check_and_set_key_float (DRMS_Record_t *new, const char *key, float val) {
  DRMS_Keyword_t *keywd;
  float vreq;
  int status;
  char sreq[128], sval[128];

  if (!(keywd = drms_keyword_lookup (new, key, 1))) return 0;
  if (keywd->info->recscope != 1) {
		       /*  it's not a constant, so don't worry, just set it  */
    drms_setkey_float (new, key, val);
    return 0;
  }
  vreq = drms_getkey_float (new, key, &status);
  if (status) {
    fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
    return 1;
  }
  sprintf (sreq, keywd->info->format, vreq);
  sprintf (sval, keywd->info->format, val);
/*
  if (vreq != val) {
*/
  if (strcmp (sreq, sval)) {
    char format[256];
    sprintf (format,
	"Error:  parameter value %s for keyword %%s\n  differs from required value %s\n",
	keywd->info->format, keywd->info->format);
    fprintf (stderr, format, val, key, vreq);
    return 1;
  }
  return 0;
}

int check_and_set_key_double (DRMS_Record_t *new, const char *key, double val) {
  DRMS_Keyword_t *keywd;
  double vreq;
  int status;
  char sreq[128], sval[128];

  if (!(keywd = drms_keyword_lookup (new, key, 1))) return 0;
  if (keywd->info->recscope != 1) {
		       /*  it's not a constant, so don't worry, just set it  */
    drms_setkey_double (new, key, val);
    return 0;
  }
  vreq = drms_getkey_double (new, key, &status);
  if (status) {
    fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
    return 1;
  }
  sprintf (sreq, keywd->info->format, vreq);
  sprintf (sval, keywd->info->format, val);
  if (strcmp (sreq, sval)) {
    char format[256];
    sprintf (format,
	"Error:  parameter value %s for keyword %%s\n  differs from required value %s\n",
	keywd->info->format, keywd->info->format);
    fprintf (stderr, format, val, key, vreq);
    return 1;
  }
  return 0;
}

int check_and_set_key_str (DRMS_Record_t *new, const char *key, char *val) {
  DRMS_Keyword_t *keywd;
  char *vreq;
  int status;

  if (!(keywd = drms_keyword_lookup (new, key, 1))) return 0;
  if (keywd->info->recscope != 1) {
		       /*  it's not a constant, so don't worry, just set it  */
    drms_setkey_string (new, key, val);
    return 0;
  }
  vreq = drms_getkey_string (new, key, &status);
  if (status) {
    fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
    return 1;
  }
  if (strcmp (vreq, val)) {
    char format[256];
    sprintf (format,
	"Error:  parameter value \"%s\" for keyword %%s\n  differs from required value \"%s\"\n",
	keywd->info->format, keywd->info->format);
    fprintf (stderr, format, val, key, vreq);
    return 1;
  }
  return 0;
}

int check_and_set_key_time (DRMS_Record_t *new, const char *key, TIME tval) {
  DRMS_Keyword_t *keywd;
  TIME treq;
  int status;
  char sreq[64], sval[64];

  if (!(keywd = drms_keyword_lookup (new, key, 1))) return 0;
  if (keywd->info->recscope != 1) {
		       /*  it's not a constant, so don't worry, just set it  */
    drms_setkey_time (new, key, tval);
    return 0;
  }
  treq = drms_getkey_time (new, key, &status);
  if (status) {
    fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
    return 1;
  }
  sprint_time (sreq, treq, keywd->info->unit, atoi (keywd->info->format));
  sprint_time (sval, tval, keywd->info->unit, atoi (keywd->info->format));
  if (strcmp (sval, sreq)) {
    fprintf (stderr, "Error:  parameter value %s for keyword %s\n", sval, key);
    fprintf (stderr, "  differs from required value %s\n", sreq);
    return 1;
  }
  return 0;
}

int check_and_copy_key (DRMS_Record_t *new, DRMS_Record_t *old,
    const char *key) {
  DRMS_Keyword_t *oldkey, *newkey;
  int status;

  if (!(oldkey = drms_keyword_lookup (old, key, 1))) return 0;
  if (!(newkey = drms_keyword_lookup (new, key, 1))) return 0;
		     /*  check that key types agree (this could be relaxed)  */
  if (oldkey->info->type != newkey->info->type) {
    fprintf (stderr, "Error: type for keyword %s differs in input and output\n",
	key);
    return 1;
  }
  if (newkey->info->recscope != 1) {
		      /*  it's not a constant, so don't worry, just copy it  */
    drms_copykey (new, old, key);
    return 0;
  }
  switch (newkey->info->type) {
/*
    case DRMS_TYPE_CHAR:
      return check_and_set_key_char (new, key,
	  drms_getkey_char (old, key, &status));
    case DRMS_TYPE_SHORT:
      return check_and_set_key_short (new, key,
	  drms_getkey_short (old, key, &status));
*/
    case DRMS_TYPE_INT:
      return check_and_set_key_int (new, key,
	  drms_getkey_int (old, key, &status));
/*
    case DRMS_TYPE_LONGLONG:
      return check_and_set_key_longlong (new, key,
	  drms_getkey_longlong (old, key, &status));
*/
    case DRMS_TYPE_FLOAT:
      return check_and_set_key_float (new, key,
	  drms_getkey_float (old, key, &status));
    case DRMS_TYPE_DOUBLE:
      return check_and_set_key_double (new, key,
	  drms_getkey_double (old, key, &status));
    case DRMS_TYPE_TIME:
      return check_and_set_key_time (new, key,
	drms_getkey_time (old, key, &status));
    case DRMS_TYPE_STRING:
      return check_and_set_key_str (new, key,
	  drms_getkey_string (old, key, &status));
    default:
      fprintf (stderr, "Error: Unknown type %d for keyword %s\n",
	    newkey->info->type, key);
      return 1;
  }
  return 0;
}
/*
 *  Create a prime key string value by concatenating the sprinted values of a
 *    selected list of keys in a DRMS record
 *  The values are ordered as in the key list, and separated by '|' characters;
 *    '|' characters in strings are converted to "\|"
 *  All keywords in the list must be present in the record, otherwise a NULL
 *    string is returned
 */

void string_insert_escape_char (char **str, const char esc) {
  int i, n, ct = 0, len;

  len = strlen (*str);
  for (n = 0; n < len; n++) if ((*str)[n] == esc) ct++;
  if (ct) {
    *str = realloc (*str, 2*len);
    for (n = len; n < 2*len; n++) (*str)[n] = '\0';
  }
  for (n = 0; n < len; n++) {
    if ((*str)[n] == esc) {
      for (i = len; i >  n; i--) (*str)[i] = (*str)[i-1];
      (*str)[n] = '\\';
      len++;
      n++;
    }
  }
}

int append_keyval_to_primekeyval (char **pkey, DRMS_Record_t *rec,
    const char *key) {
  DRMS_Keyword_t *keywd;
  int size;
  static int curlen;
  char **buf;
  int buflen = 256;

  if (!*pkey) {
    *pkey = calloc (buflen, sizeof (char));
    curlen = buflen;
  }
  size = strlen (*pkey);
  if (size) {
					/*  after first key, append separator  */
    size++;
    if (size >= curlen) {
      char *tmp = calloc (curlen, sizeof (char));
      strcpy (tmp, *pkey);
      curlen += buflen;
      *pkey = realloc (*pkey, curlen);
      bzero (*pkey, curlen);
      strcpy (*pkey, tmp);
      free (tmp);
    }
    strcat (*pkey, "|");
  }
  if (!(keywd = drms_keyword_lookup (rec, key, 1))) {
    fprintf (stderr, "create_primekey_from_keylist(): %s not found\n", key);
    if (*pkey) free (*pkey);
    *pkey = NULL;
    return 1;
  }
  buf = malloc (sizeof (char **));
  *buf = malloc (DRMS_DEFVAL_MAXLEN * sizeof (char *));
  drms_keyword_snprintfval (keywd, *buf, DRMS_DEFVAL_MAXLEN);
  if (keywd->info->type == DRMS_TYPE_STRING ||
      keywd->info->type == DRMS_TYPE_TIME) {
    string_insert_escape_char (buf, '|');
  }
  size += strlen (*buf);
  if (size >= curlen) {
    curlen +- buflen;
    *pkey = realloc (*pkey, curlen);
  }
  strncat (*pkey, *buf, strlen (*buf));
  free (*buf);
  free (buf);
  return 0;
}

char *create_primekey_from_keylist (DRMS_Record_t *rec, char **keylist,
    int keyct) {
  char *pkey = NULL;
  int n;
  for (n = 0; n < keyct; n++)
    if (append_keyval_to_primekeyval (&pkey, rec, keylist[n])) break;

  return pkey;
}

int drms_wcs_timestep (DRMS_Record_t *rec, int axis, double *tstep) {
  double dval;
  int status;
  char *sval;
  char delta[9], unit[9];

  *tstep = FP_NAN;
  sprintf (delta, "CDELT%d", axis);
  sprintf (unit, "CUNIT%d", axis);
  dval = drms_getkey_double (rec, delta, &status);
  if (status || isnan (dval)) return 1;
  *tstep = dval;
  sval = drms_getkey_string (rec, unit, &status);

  if (status || !sval) {
    fprintf (stderr, "Warning: no keyword %s: \"s\" assumed\n", unit);
    return 0;
  }
  if (!strcmp (sval, "s")) return 0;
  if (!strcmp (sval, "min")) {
    *tstep *= 60.0;
    return 0;
  }
  if (!strcmp (sval, "h")) {
    *tstep *= 3600.0;
    return 0;
  }
  if (!strcmp (sval, "d")) {
    *tstep *= 86400.0;
    return 0;
  }
  if (!strcmp (sval, "a") || !strcmp (sval, "yr")) {
    *tstep *= 31557600.0;
    return 0;
  }

  if (!strcmp (sval, "ns")) {
    *tstep *= 1.0e-9;
    return 0;
  }
  if (!strcmp (sval, "us")) {
    *tstep *= 1.0e-6;
    return 0;
  }
  if (!strcmp (sval, "ms")) {
    *tstep *= 0.001;
    return 0;
  }
  if (!strcmp (sval, "cs")) {
    *tstep *= 0.01;
    return 0;
  }
  if (!strcmp (sval, "ds")) {
    *tstep *= 0.1;
    return 0;
  }
  if (!strcmp (sval, "das")) {
    *tstep *= 10.0;
    return 0;
  }
  if (!strcmp (sval, "hs")) {
    *tstep *= 100.0;
    return 0;
  }
  if (!strcmp (sval, "ks")) {
    *tstep *= 1000.0;
    return 0;
  }
  if (!strcmp (sval, "Ms")) {
    *tstep *= 1.0e6;
    return 0;
  }
  if (!strcmp (sval, "gs")) {
    *tstep *= 1.0e9;
    return 0;
  }

  if (!strcmp (sval, "nmin")) {
    *tstep *= 6.0e-8;
    return 0;
  }
  if (!strcmp (sval, "umin")) {
    *tstep *= 6.0e-5;
    return 0;
  }
  if (!strcmp (sval, "mmin")) {
    *tstep *= 0.06;
    return 0;
  }
  if (!strcmp (sval, "cmin")) {
    *tstep *= 0.6;
    return 0;
  }
  if (!strcmp (sval, "dmin")) {
    *tstep *= 6.0;
    return 0;
  }
  if (!strcmp (sval, "damin")) {
    *tstep *= 600.0;
    return 0;
  }
  if (!strcmp (sval, "hmin")) {
    *tstep *= 6000.0;
    return 0;
  }
  if (!strcmp (sval, "kmin")) {
    *tstep *= 60000.0;
    return 0;
  }
  if (!strcmp (sval, "Mmin")) {
    *tstep *= 6.0e7;
    return 0;
  }
  if (!strcmp (sval, "gmin")) {
    *tstep *= 6.0e10;
    return 0;
  }

  if (!strcmp (sval, "nh")) {
    *tstep *= 3.6e-6;
    return 0;
  }
  if (!strcmp (sval, "uh")) {
    *tstep *= 3.6e-3;
    return 0;
  }
  if (!strcmp (sval, "mh")) {
    *tstep *= 0.36;
    return 0;
  }
  if (!strcmp (sval, "ch")) {
    *tstep *= 3.6;
    return 0;
  }
  if (!strcmp (sval, "dh")) {
   *tstep *= 360.0;
   return 0;
  }
  if (!strcmp (sval, "dah")) {
    *tstep *= 3.6e4;
    return 0;
  }
  if (!strcmp (sval, "hh")) {
    *tstep *= 3.6e5;
    return 0;
  }
  if (!strcmp (sval, "kh")) {
    *tstep *= 3.6e6;
    return 0;
  }
  if (!strcmp (sval, "Mh")) {
    *tstep *= 3.6e9;
    return 0;
  }
  if (!strcmp (sval, "gh")) {
    *tstep *= 3.6e12;
    return 0;
  }

  if (!strcmp (sval, "nd")) {
    *tstep *= 8.64e-5;
    return 0;
  }
  if (!strcmp (sval, "ud")) {
    *tstep *= 8.64-2;
    return 0;
  }
  if (!strcmp (sval, "md")) {
    *tstep *= 86.4;
    return 0;
  }
  if (!strcmp (sval, "cd")) {
    *tstep *= 864.0;
    return 0;
  }
  if (!strcmp (sval, "dd")) {
    *tstep *= 8640.0;
    return 0;
  }
  if (!strcmp (sval, "dad")) {
    *tstep *= 8.64e5;
    return 0;
  }
  if (!strcmp (sval, "hd")) {
    *tstep *= 8.64e6;
    return 0;
  }
  if (!strcmp (sval, "kd")) {
    *tstep *= 8.64e7;
    return 0;
  }
  if (!strcmp (sval, "Md")) {
    *tstep *= 8.64e10;
    return 0;
  }
  if (!strcmp (sval, "gd")) {
    *tstep *= 8.64e13;
    return 0;
  }

  fprintf (stderr, "Error: %s type of %s unrecognized by drms_wcs_timestep()\n",
      unit, sval);
  return 1;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  09.04.27 or earlier	created file
 *  09.05.16		added functions for creating prime key string
 *  09.11.06		added WCS parsing function
 *  09.11.09		modified check_and_set_key_* to only require matching
 *		to precision of keyword format for float, double, & time
 */
