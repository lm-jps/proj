/*
 *  smpl_03.c						$DRMS/proj/cookbook/
 *
 *  Prints a list of the series names known to DRMS, with record counts
 *    for each
 *  Illustrates the connection to the DRMS database and ways that SQL
 *    queries can be directly run and the results analyzed at the lowest
 *    level of the DRMS API
 *
 *  Usage:
 *    smpl_03 [nmax= ...]
 *
 *  Revision history is at end of file.
 */

char *module_name = "CookbookRecipe:03";
char *version_id = "1.0";

#include <jsoc_main.h>
#include <regex.h>

ModuleArgs_t module_args[] = {
  {ARG_INT,	"nmax",	"100", "maximum number of series to be listed"},
  {ARG_END}
};

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DB_Text_Result_t *qres, *sqres;
  int series, seriesct;
  char query[DRMS_MAXQUERYLEN];

  int nmax = params_get_int (params, "nmax");
  	/*  Query the database to get all series names from the master list  */
  sprintf (query, "select seriesname from %s()", DRMS_MASTER_SERIES_TABLE);
  if ((qres = drms_query_txt (drms_env->session, query)) == NULL) {
    fprintf (stderr, "Cant find DRMS\n");
    return 1;
  }
  seriesct = qres->num_rows;
  printf ("%d series found", seriesct);
  if (seriesct > nmax) {
    seriesct = nmax;
    printf (" (only the first %d will be listed)", seriesct);
  }
  printf ("\n");

  for (series = 0; series < seriesct; series++) {
    char *seriesname = qres->field[series][0];
    printf ("%s\t", seriesname);
    sprintf (query, "select count (recnum) from %s", seriesname);
       /*  Query the database to get the record count from the series table  */
			  /*  (every data series must have a "recnum" field) */
    if (sqres = drms_query_txt (drms_env->session, query)) {
      printf ("%s", sqres->field[0][0]);
      db_free_text_result (sqres);
    } else printf ("?");
    printf ("\n");
  }

  db_free_text_result (qres);
  return 0;
}

/*
 *  Revision History
 *
 *  09.04.20	file created by R Bogart
 */
