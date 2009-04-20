/*
 *  smpl_03.c						$DRMS/proj/cookbook/
 *
 *  Prints the number of unique records in the selected data series, and
 *    the defined number of data segments per record for the series
 *  Illustrates features of the DRMS_Record struct, and concept of uniqueness
 *    for DRMS records
 *
 *  Usage:
 *    smpl_03 ds= ...
 *
 *  Revision history is at end of file.
 */

char *module_name = "CookbookRecipe:03";
char *version_id = "1.0";

#include <jsoc_main.h>
#include <regex.h>

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds", "", "name of data series"},
  {ARG_END}
};

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DB_Text_Result_t *qres;
  DRMS_Record_t *record;
  regmatch_t pmatch[10];
  int series, seriesct;
  int status = 0;
  char query[DRMS_MAXQUERYLEN];

  char *ds = params_get_str (params, "ds");

  record = drms_template_record (drms_env, ds, &status);
  if (record && !status) {
    if (record->seriesinfo->pidx_num) {
      char qry[1024];
      sprintf (query, "select %s from %s group by %s",
	  (record->seriesinfo->pidx_keywords[0])->info->name, ds,
	  (record->seriesinfo->pidx_keywords[0])->info->name);
      for (int i=1; i<record->seriesinfo->pidx_num; i++) {
	sprintf (qry, ", %s",
	    (record->seriesinfo->pidx_keywords[i])->info->name);
	strcat (query, qry);
      }
      if ((qres = drms_query_txt (drms_env->session, query))) {
	printf ("%s contains %d unique records", ds, qres->num_rows);
	db_free_text_result (qres);
      }
    } else printf ("%s: no records found", ds);
    printf (", with %d data segment(s) per record\n", record->segments.num_total);
    drms_free_record (record);
  } else printf ("%s: no such series\n", ds);

  return status;
}

/*
 *  Revision History
 *
 *  09.04.20	file created by R Bogart
 */
