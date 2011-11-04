/*
 *  smpl_05.c						$DRMS/proj/cookbook/
 *
 *  Populates the data series created from noaa_ar.jsd with the data in
 *    noaa_ar.dat
 *  Illustrates features of the DRMS_Keyword struct, and writing to DRMS
 *    records
 *
 *  Usage:
 *    smpl_05 [ds= dat= ]
 *
 *  Revision history is at end of file.
 */


char *module_name = "CookbookRecipe:05";
char *version_id = "1.0";

#include <jsoc_main.h>

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds",   "drms.noaa_ar", "name of data series"},
  {ARG_STRING,	"data", "noaa_ar.dat",  "file containing input data"},
  {ARG_FLAG, "v", "", "verbose flag"}, 
  {ARG_END}
};


int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_Record_t *record;
  FILE *datafile;
  TIME obstime;
  int arnum, lathg, lonhg, loncm, area, extent, count;
  char timestr[64], class[16], type[16];
  size_t linelen = 256;
  char *line = (char *)malloc (linelen);
  int rct = 0, status = 0;

  char *series = strdup (params_get_str (params, "ds"));
  char *filename = strdup (params_get_str (params, "data"));
  int verbose = cmdparams_exists (&cmdparams, "v");

  datafile = fopen (filename, "r");

  while (getline (&line, &linelen, datafile) >=0) {
    sscanf (line, "%s %d %d %d %d %d %s %d %d %s", timestr, &arnum,
	&lathg, &loncm, &lonhg, &area, class, &extent, &count, type);
    obstime = sscan_time (timestr);
    record = drms_create_record (drms_env, series, DRMS_PERMANENT, &status);
    if (!record) {
      fprintf (stderr, "Error: Unable to create record %d in series %s\n",
	rct, series);
      fprintf (stderr, "       Does series exist and is it writeable by you?\n");
      return 1;
    }
    drms_setkey_time (record, "Date", obstime);
    drms_setkey_int (record, "Region", arnum);
    drms_setkey_int (record, "Lat", lathg);
    drms_setkey_int (record, "Lon", loncm);
    drms_setkey_int (record, "CarrLon", lonhg);
    drms_setkey_int (record, "Area", area);
    drms_setkey_int (record, "Extent", extent);
    drms_setkey_int (record, "Count", count);
    drms_setkey_string (record, "Class", class);
    drms_setkey_string (record, "Type", type);
    drms_close_record (record, DRMS_INSERT_RECORD);
    rct++;
  }

  fclose (datafile);
  if (verbose) printf ("wrote %d records\n", rct);
  return status;
}
/*
 *  Revision History
 *
 *  09.07.29	file created by R Bogart
 */
