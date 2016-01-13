#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"
#include <unistd.h>
#include <strings.h>

#define NOTSPECIFIED "***NOTSPECIFIED***"
/* List of default parameter values. */
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "ds", NOTSPECIFIED, "", ""},   /* Series to store into */
  {ARG_STRING, "in", NOTSPECIFIED, "", ""},   /* File to store from */
  {ARG_STRING, "file_name", NOTSPECIFIED, "", ""}, /* Title of file being backed up */
  {ARG_STRING, "src_dir", NOTSPECIFIED, "", ""}, /* Title of file being backed up */
  {ARG_STRING, "backup_type", NOTSPECIFIED, "", ""}, /* Type of backup */
  {ARG_INT,    "backup_id", "0", "", ""}, /* backup_id, if appending */
  {ARG_INT,    "file_id", "0", "", ""}, /* file_id, index starts at 1 */
  {ARG_INT,    "chunk_id", "0", "", ""}, /* chunk_id, index starts at 1 */
  {ARG_INT,    "chunk_max", "0", "", ""}, /* chunk_id, index starts at 1 */
  {ARG_INT,    "chunk_size", "0", "", ""}, /* file size in bytes*/
  {ARG_TIME,   "backup_date", NOTSPECIFIED, "date backup was taken"},
  {ARG_TIME,   "date", NOTSPECIFIED, "date backup was put into SUMS"},
  {ARG_FLAG,   "s", "0", "", ""},   /* start a new backup session (new backup_id) */
  {ARG_FLAG,   "h", "0", "", ""},   /* print usage message and quit */
  {ARG_FLAG,   "v", "0", "", ""},   /* verbose flag, normally do not use  */
  {ARG_END}
};
//  {ARG_DOUBLE, "backup_date", "0", ""}, /* backup date */
//  {ARG_DOUBLE, "date", "0", ""}, /* archive date */
/* Module name presented to DRMS. */
char *module_name = "store_lmbackup_v2";

/* Some global variables for this module. */
int verbose = 0;

/* Check DRMS session status and calling arguments and set verbose variable */
int nice_intro () {
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  if (usage) {
    printf ("store_lmbackup ds=<series> in=<file> backup_title=<title> file_name=<name> backup_type=<type> {-h} {-v}\n"
	    "  -s: STart new backup\n"
	    "  -h: print this message and exit\n"
	    "  -v: verbose\n"
	    "ds=    - data series to store file into\n"
	    "in=    - file to store\n"
	    "sel=   - optional description for retrieval key\n");
    return (1);
  }
  if (verbose) cmdparams_printall (&cmdparams);
  return (0);
}

/* Module main function. */
int DoIt (void) {
  int status = 0;

  const char *in, *filename;
  char *series, *note, *sel, *rsp;
  DRMS_Record_t *rec, *template;
  int yes_create = cmdparams_get_int (&cmdparams, "c", NULL);


  int backup_id = 0;
  int file_id =
    cmdparams_get_int(&cmdparams, "file_id", NULL);
  int chunk_id =
    cmdparams_get_int(&cmdparams, "chunk_id", NULL);
  int chunk_max =
    cmdparams_get_int(&cmdparams, "chunk_max", NULL);
  int chunk_size =
    cmdparams_get_int(&cmdparams, "chunk_size", NULL);
  const char* file_name =
    cmdparams_get_str(&cmdparams, "file_name", NULL);
  const char* src_dir =
    cmdparams_get_str(&cmdparams, "src_dir", NULL);
  const char* backup_type =
    cmdparams_get_str(&cmdparams, "backup_type", NULL);
  const char* backup_date =
    cmdparams_get_str(&cmdparams, "backup_date", &status);
  const char* date=
    cmdparams_get_str(&cmdparams, "date", &status);


  if (nice_intro())
    return (0);

  series = strdup(cmdparams_get_str(&cmdparams, "ds", NULL));
  in = cmdparams_get_str(&cmdparams, "in", NULL);
  if (strcmp(series, NOTSPECIFIED) == 0)
    {
      fprintf(stderr, "'ds' series must be specified.  Abort\n");
      return(1);
    }
  if (strcmp(in, NOTSPECIFIED) == 0)
    {
      fprintf(stderr, "'in' file must be specified.  Abort\n");
      return(1);
    }

  /* First, check to see that the file exists */
  if (access(in, R_OK) != 0)
    {
      printf("The requested file can not be accessed.  %s\n", in);
      return(1);
    }

  /* Now, extract the filename from the path */
  filename = rindex(in, '/');
  if (filename)
    filename++;
  else
    filename = in;

  /* Check to see if series exists */
  template = drms_template_record(drms_env, series, &status);
  if (template==NULL && status == DRMS_ERROR_UNKNOWNSERIES)
    {
      printf("Series '%s' does not exist. Give up.\n", series);
      return(1);
    }
  else
    if(status)
      {
	fprintf(stderr, "DRMS problem looking up series.\n");
	return(status);
      }

  /*
   * If new backup session, automatically increment the backup_id
   * by first querying for it
   */
  if (cmdparams_get_int(&cmdparams, "s", NULL))
    {
      DRMS_RecordSet_t *recordset =
	drms_open_records(drms_env, series, &status);
      if(!recordset)
	{
	  printf("Series '%s' does not exist.\n", series);
	  return(1);
	}
      if(recordset->n == 0)
	{
	  // printf("Default backup_id to 1\n");
	  backup_id = 1;
	}
      else
	{
	  int rec_n, cur_bkid, max_bkid = 0;
	  for(rec_n = 0; rec_n < recordset->n; rec_n++)
	    {
	      rec = recordset->records[(recordset->n)-1];
	      cur_bkid = drms_getkey_int(rec, "backup_id", &status);
	      if(status)
		{
		  fprintf(stderr, "DRMS Problem getting last backup_id");
		  return(1);
		}
	      if(cur_bkid > max_bkid)
		{
		  max_bkid = cur_bkid;
		}
	    }
	  backup_id = max_bkid+1;
	  // printf("Max bkid = %i\n",max_bkid);
	  /*
	   * printf("Reading from record #%i\n", (recordset->n)-1);
	   * rec = recordset->records[(recordset->n)-1];
	  */
	  /*
	  DRMS_Segment_t *seg = drms_segment_lookupnum(rec, 0);
	  char path[DRMS_MAXPATHLEN];
	  drms_record_directory(rec,path,1);
	  drms_segment_filename(seg, path);
	  printf("Next available backup_id is %i, last seg in %s\n",backup_id, path);
	  */
	}
    }
  else
    {
      backup_id = cmdparams_get_int(&cmdparams, "backup_id", NULL);
      if (backup_id <= 0)
	{
	  fprintf(stderr, "Need to specify backup_id when continuing a backup session\n");
	  return(1);
	}
    }
  printf("%i\n",backup_id);

  if (file_id == 0)
    {
      fprintf(stderr,"Need to specify file_id\n");
      return(1);
    }
  if (chunk_id == 0)
    {
      fprintf(stderr,"Need to specify chunk_id\n");
      return(1);
    }
  if (chunk_max == 0)
    {
      fprintf(stderr,"Need to specify chunk_max\n");
      return(1);
    }
  if (chunk_size == 0)
    {
      fprintf(stderr,"Need to specify chunk_size\n");
      return(1);
    }
  if (strcmp(file_name, NOTSPECIFIED) == 0)
    {
      fprintf(stderr,"Need to specify file_name\n");
      return(1);
    }
  if (strcmp(backup_type, NOTSPECIFIED) == 0)
    {
      fprintf(stderr,"Need to specify backup_type\n");
      return(1);
    }
  if (strcmp(backup_date, NOTSPECIFIED) == 0)
    {
      fprintf(stderr,"Need to specify valid backup_date\n");
      return(1);
    }
  if (strcmp(date, NOTSPECIFIED) == 0)
    {
      fprintf(stderr,"Need to specify valid date\n");
      return(1);
    }

  /* Now ready to make a new record and set keywords */
  rec = drms_create_record(drms_env, series, DRMS_PERMANENT, &status);
  if (!rec || status)
    {
      printf("drms_create_record failed, series=%s, status=%d.  Aborting.\n", series,status);
      return(status);
    }
  if ((status = drms_setkey_int(rec, "backup_id", backup_id)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'backup_id'\n", series);
      else
	printf("ERROR: drms_setkey_int failed for 'backup_id'\n");
      return(1);
    }
  if ((status = drms_setkey_int(rec, "file_id", file_id)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'file_id'\n", series);
      else
	printf("ERROR: drms_setkey_int failed for 'file_id'\n");
      return(1);
    }
  if ((status = drms_setkey_int(rec, "chunk_id", chunk_id)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'chunk_id'\n", series);
      else
	printf("ERROR: drms_setkey_int failed for 'chunk_id'\n");
      return(1);
    }
  if ((status = drms_setkey_int(rec, "chunk_max", chunk_max)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'chunk_max'\n", series);
      else
	printf("ERROR: drms_setkey_int failed for 'chunk_max'\n");
      return(1);
    }
  if ((status = drms_setkey_int(rec, "chunk_size", chunk_size)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'chunk_size'\n", series);
      else
	printf("ERROR: drms_setkey_int failed for 'chunk_size'\n");
      return(1);
    }
  if (strcmp(src_dir, NOTSPECIFIED) == 1)
    {
      if ((status = drms_setkey_string(rec, "src_dir", src_dir)))
	{
	  if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	    printf("ERROR: series %s does not have a keyword named 'src_dir'\n", series);
	  else
	    printf("ERROR: drms_setkey_string failed for 'src_dir'\n");
	  return(1);
	}
    }
  else
    {
      if ((status = drms_setkey_string(rec, "src_dir", " ")))
	{
	  if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	    printf("ERROR: series %s does not have a keyword named 'src_dir'\n", series);
	  else
	    printf("ERROR: drms_setkey_string failed for 'src_dir'\n");
	  return(1);
	}
    }
  if ((status = drms_setkey_string(rec, "file_name", file_name)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'file_name'\n", series);
      else
	printf("ERROR: drms_setkey_string failed for 'file_name'\n");
      return(1);
    }
  if ((status = drms_setkey_string(rec, "backup_type", backup_type)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'backup_type'\n", series);
      else
	printf("ERROR: drms_setkey_string failed for 'backup_type'\n");
      return(1);
    }
  if ((status = drms_setkey_string(rec, "backup_date", backup_date)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'backup_date'\n", series);
      else
	printf("ERROR: drms_setkey_double failed for 'backup_date'\n");
      return(1);
    }
  if ((status = drms_setkey_string(rec, "date", date)))
    {
      if (status == DRMS_ERROR_UNKNOWNKEYWORD)
	printf("ERROR: series %s does not have a keyword named 'date'\n", series);
      else
	printf("ERROR: drms_setkey_double failed for 'date'\n");
      return(1);
    }
  if ((status = drms_segment_write_from_file(drms_segment_lookup(rec, "file_seg"), (char*)in)))
    {
      printf("drms_segment_write_file failed, status=%d\n",status);
      return(status);
    }

  if ((status = drms_close_record(rec, DRMS_INSERT_RECORD)))
    printf("drms_close_record failed!\n");

  return(status);
}
