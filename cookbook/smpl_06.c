/*
 *  smpl_06.c						$DRMS/proj/cookbook/
 *
 *  Updates a record in the data series created in smpl_05 with new data
 *
 *  Illustrates features of the DRMS_Keyword struct, and "updating" of
 *    records in DRMS via cloning and recnum; also illustrates processing
 *    of multiple records via the DRMS_RecordSet struct
 *
 *  Usage:
 *    smpl_06 [ds= cols= rows= ]
 *
 *  Revision history is at end of file.
 */


char *module_name = "CookbookRecipe:06";
char *version_id = "1.0";

#include <jsoc_main.h>

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds",   "drms.noaa_ar", "name of data series"},
  {ARG_TIME,    "date", "Not Specified", "report date"},
  {ARG_INT,     "ar", "Not Specified", "active region number"},
  {ARG_STRING,	"key",  "", "key to be modified (required)"},
  {ARG_STRING,	"value",  "", "key value to be set (required)"},
  {ARG_FLAG,	"f",	"", "force update of multiple records"}, 
  {ARG_END}
};

int setkey_general (DRMS_Record_t *rec, const char *key, char *val) {
  DRMS_Keyword_t *keywd = drms_keyword_lookup (rec, key, 0);
  double dv;
  long long lv;
  int iv;

  if (!keywd) {
    fprintf (stderr, "Error: key %s not found\n", key);
    return 1;
  }

  switch (drms_keyword_gettype (keywd)) {
    case DRMS_TYPE_STRING:
      return drms_setkey_string (rec, key, val);
    case DRMS_TYPE_TIME:
      return drms_setkey_time (rec, key, sscan_time (val));
    case DRMS_TYPE_DOUBLE:
      dv = atof (val);
      return drms_setkey_double (rec, key, dv);
    case DRMS_TYPE_FLOAT:
      dv = atof (val);
      return drms_setkey_float (rec, key, (float)dv);
    case DRMS_TYPE_INT:
      iv = atoi (val);
      return drms_setkey_int (rec, key, iv);
    case DRMS_TYPE_SHORT:
      iv = atoi (val);
      return drms_setkey_int (rec, key, (short)iv);
    case DRMS_TYPE_CHAR:
      iv = atoi (val);
      return drms_setkey_int (rec, key, (char)iv);
    case DRMS_TYPE_LONGLONG:
      lv = atoll (val);
      return drms_setkey_longlong (rec, key, lv);
    default:
      fprintf (stderr, "Error: unknown key type for key %s\n", key);
      return 1;
  }

  return 1;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ds;
  DRMS_Record_t *record;
  int rec, recct, status = 0;

  char *dspec = params_get_str (params, "ds");
  TIME date = params_get_time (params, "date");
  int arnum = params_get_int (params, "ar");
  char *keynam = params_get_str (params, "key");
  char *keyval = params_get_str (params, "value");
  int force = params_isflagset (params, "f");

  if (!(ds = drms_open_records (drms_env, dspec, &status))) {
    fprintf (stderr, "Error: unable to open input data set %s\n", dspec);
    return 1;
  }

  if (!drms_ismissing_time (date) && !drms_ismissing_int (arnum)) {
	    /*  if valid values supplied for both date and ar, form a query  */
    char newspec[DRMS_MAXSERIESNAMELEN], timestr[16];
    char *brack = strchr (dspec, '[');	       /* strip existing query spec  */
    if (brack) *brack = '\0';
    drms_close_records (ds, DRMS_FREE_RECORD);
    date += 43200.0;				    /* round to nearest day  */
    sprint_time (timestr, date, "Z", -3);
    sprintf (newspec, "%s[%s][%d]", dspec, timestr, arnum);
    if (!(ds = drms_open_records (drms_env, newspec, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", newspec);
      return 1;
    }
  }
  if ((recct = ds->n) > 1) {
    if (force) printf ("Updating %d records:\n", recct);
    else {
      fprintf (stderr, "Warning: %d records matched query.\n", recct);
      fprintf (stderr,
	  "         To update all of these records, use the -f flag\n");
      fprintf (stderr,
	  "         You may limit the query with dataset specification notation\n");
      fprintf (stderr,
	  "           or by specifying *both* the ar and date parameters\n");
      return 1;
    }
  }
  for (rec = 0; rec < recct; rec++) {
    record = drms_clone_record (ds->records[rec], DRMS_PERMANENT,
	DRMS_SHARE_SEGMENTS, &status);
    if (status = setkey_general (record, keynam, keyval)) {
      fprintf (stderr, "Error: unable to set %s = %s for record %d\n", keynam,
	  keyval, rec);
      drms_close_record (record, DRMS_FREE_RECORD);
    } else drms_close_record (record, DRMS_INSERT_RECORD);
  }
  return 0;
}
/*
 *  Revision History
 *
 *  09.08.03	file created by R Bogart
 */
