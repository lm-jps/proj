#define DEBUG 1
#define DEBUG 0

/*
 *  jsoc_fetch - cgi-bin program to recieve jsoc export and export status requests
 *
 */
#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"
#include "json.h"

#define EXPORT_SERIES "jsoc.export"
#define EXPORT_SERIES_NEW "jsoc.export_new"
#define EXPORT_USER "jsoc.export_user"

static char x2c (char *what) {
  register char digit;

  digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
  digit *= 16;
  digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0'));
  return (digit);
}

static void CGI_unescape_url (char *url) {
  register int x, y;

  for (x = 0, y = 0; url[y]; ++x, ++y) {
    if ((url[x] = url[y]) == '%') {
      url[x] = x2c (&url[y+1]);
      y += 2;
    }
  }
  url[x] = '\0';
}

int drms_count_records(DRMS_Env_t *env, char *recordsetname, int *status)
  {
  int stat, filter, mixed;
  char *query=NULL, *where=NULL, *seriesname=NULL;
  int count = 0;
  DB_Text_Result_t *tres;

  stat = drms_recordset_query(env, recordsetname, &where, &seriesname, &filter, &mixed);
      if (stat)
        goto failure;

  stat = 1;
  query = drms_query_string(env, seriesname, where, filter, mixed, DRMS_QUERY_COUNT, NULL);
      if (!query)
        goto failure;

  tres = drms_query_txt(env->session,  query);

  if (tres && tres->num_rows == 1 && tres->num_cols == 1)
    count = atoi(tres->field[0][0]);
  else
    goto failure;

  free(seriesname);
  free(query);
  free(where);
  *status = DRMS_SUCCESS;
  return(count);

  failure:
  if (seriesname) free(seriesname);
  if (query) free(query);
  if (where) free(where);
  *status = stat;
  return(0);
  }


ModuleArgs_t module_args[] =
{ 
  {ARG_STRING, "op", "Not Specified", "<Operation>"},
  {ARG_STRING, "ds", "Not Specified", "<record_set query>"},
  {ARG_STRING, "requestid", "Not Specified", "JSOC export request identifier"},
  {ARG_STRING, "process", "Not Specified", "string containing program and arguments"},
  {ARG_STRING, "requestor", "Not Specified", "name of requestor"},
  {ARG_STRING, "notify", "Not Specified", "email address of requestor"},
  {ARG_STRING, "shipto", "Not Specified", "mail address of requestor"},
  {ARG_STRING, "protocol", "as-is", "exported file protocol"},
  {ARG_STRING, "format", "json", "return content type"},
  {ARG_STRING, "method", "url", "return method"},
  {ARG_FLAG, "h", "0", "help - show usage"},
  {ARG_STRING, "QUERY_STRING", "Not Specified", "AJAX query from the web"},
  {ARG_END}
};

char *module_name = "jsoc_fetch";

int nice_intro ()
  {
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\njsoc_info {-h} "
        "  details are:\n"
	"op=<command> tell which ajax function to execute\n"
	"ds=<recordset query> as <series>{[record specifier]} - required\n"
	"requestid=JSOC export request identifier\n"
	"process=string containing program and arguments\n"
	"requestor=name of requestor\n"
	"notify=email address of requestor\n"
	"shipto=mail address of requestor\n"
	"protocol=exported file protocol\n"
	"format=return content type\n"
	"method=return method\n"
	"h=help - show usage\n"
	"QUERY_STRING=AJAX query from the web"
	);
    return(1);
    }
  return (0);
  }

void drms_sprint_rec_query(char *text, DRMS_Record_t *rec)
  {
  int iprime, nprime=0;
  char **external_pkeys, *pkey;
  DRMS_Keyword_t *rec_key;
  if (!rec)
    {
    sprintf(text, "** No Record **");
    return;
    }
  strcpy(text,rec->seriesinfo->seriesname);
  external_pkeys =
        drms_series_createpkeyarray(rec->env, rec->seriesinfo->seriesname, &nprime, NULL);
  if (external_pkeys && nprime > 0)
    {
    for (iprime = 0; iprime < nprime; iprime++)
      {
      char val[1000];
      pkey = external_pkeys[iprime];
      rec_key = drms_keyword_lookup (rec, pkey, 1);
      drms_keyword_snprintfval(rec_key, val, sizeof(val));
      strcat(text, "[");
      strcat(text, val);
      strcat(text, "]");
      }
    }
  else
    sprintf(text, "[:#%lld]",rec->recnum);
  return;
  }

/* quick export of recordset - on entry it is known that all the records are online
 * so no directory file need be built.  The json structure will have 3 elements
 * added:
 *   size : size in bytes
 *   count : number of files
 *   data  : array of count objects containing
 *         record : record query with segment name suffix
 *         filename : path to segment file
 * set online=0 if known to be online, to 1 if should block until online
 * returns number of files found.
 */

int quick_export_rs( json_t *jroot, DRMS_RecordSet_t *rs, int online,  long long size)
  {
  char numval[200];
  char query[DRMS_MAXQUERYLEN];
  char record[DRMS_MAXQUERYLEN];
  char recpath[DRMS_MAXPATHLEN];
  char segpath[DRMS_MAXPATHLEN];
  int i;
  int count;
  json_t *data;
  DRMS_Record_t *rec;
  data = json_new_array();
  count = 0;
  for (i=0; i < rs->n; i++)
    {
    DRMS_Segment_t *seg;
    HIterator_t hit;
    rec = rs->records[i];
    drms_sprint_rec_query(query, rec);
    drms_record_directory(rec, recpath, online);
    hiter_new (&hit, &rec->segments);
    while ((seg = (DRMS_Segment_t *)hiter_getnext (&hit)))
      {
      json_t *recobj = json_new_object();
      char *jsonstr;
      count += 1;
      strcpy(record, query);
      strcat(record, "{");
      strcat(record, seg->info->name);
      strcat(record, "}");
      strncpy(segpath, recpath, DRMS_MAXPATHLEN);
      strncat(segpath, "/", DRMS_MAXPATHLEN);
      strncat(segpath, seg->filename, DRMS_MAXPATHLEN);
      jsonstr = string_to_json(record);
      json_insert_pair_into_object(recobj, "record", json_new_string(jsonstr));
      free(jsonstr);
      jsonstr = string_to_json(segpath);
      json_insert_pair_into_object(recobj, "filename", json_new_string(jsonstr));
      free(jsonstr);
      json_insert_child(data, recobj);
      }
    }
  sprintf(numval, "%ld", count);
  json_insert_pair_into_object(jroot, "count", json_new_number(numval));
  sprintf(numval, "%ld", size);
  json_insert_pair_into_object(jroot, "size", json_new_number(numval));
  json_insert_pair_into_object(jroot, "data", data);
  return(count);
  }

char *string_to_json(char *in)
  {
  char *new, *c = in;
  wchar_t *tmp, *work;
  tmp = work = (wchar_t *)malloc(sizeof(wchar_t)*(strlen(in)+1));
  while (c && *c) *tmp++ = *c++;
  *tmp = NULL;
  new = json_escape(work);
  free(work);
  return(new);
  }

#include <time.h>
TIME timenow()
  {
  TIME UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
  TIME now = (double)time(NULL) + UNIX_epoch;
  return(now);
  }

#define JSONDIE(msg) \
  {	\
  char *msgjson;	\
  char errval[10];	\
  char *json;	\
  json_t *jroot = json_new_object();	\
if (DEBUG) fprintf(stderr,"%s\n",msg);	\
  msgjson = string_to_json(msg);	\
  json_insert_pair_into_object(jroot, "status", json_new_number("4"));	\
  json_insert_pair_into_object(jroot, "error", json_new_string(msgjson));	\
  json_tree_to_string(jroot,&json);	\
  printf("Content-type: application/json\n\n");	\
  printf("%s\n",json);	\
  fflush(stdout);	\
  return(1);	\
  }

/* Module main function. */
int DoIt(void)
  {
						/* Get command line arguments */
  char *op;
  char *in;
  char *requestid;
  char *process;
  char *requestor;
  long requestorid;
  char *notify;
  char *format;
  char *shipto;
  char *method;
  char *protocol;
  char *errorreply;
  int size;
  TIME reqtime;
  TIME esttime;
  TIME exptime;
  TIME now;
  double waittime;
  char *web_query;
  int from_web,status;
  int dojson=1, dotxt=0, dohtml=0, doxml=0;
  DRMS_RecordSet_t *exports;
  DRMS_Record_t *export_log;
  char new_requestid[200];
  char status_query[1000];
  char *export_series; 

  if (nice_intro ()) return (0);

  web_query = strdup (cmdparams_get_str (&cmdparams, "QUERY_STRING", NULL));
  from_web = strcmp (web_query, "Not Specified") != 0;

  if (from_web)
    {
    char *getstring, *ds, *p;
    CGI_unescape_url(web_query);
    getstring = strdup (web_query);
    for (p=strtok(getstring,"&"); p; p=strtok(NULL, "&"))
      {
      char *key=p, *val=index(p,'=');
      if (!val)
	 JSONDIE("Bad QUERY_STRING");
      *val++ = '\0';
      cmdparams_set(&cmdparams, key, val);
      }
    free(getstring);
    }
  free(web_query);

  op = cmdparams_get_str (&cmdparams, "op", NULL);
  requestid = cmdparams_get_str (&cmdparams, "requestid", NULL);
  in = cmdparams_get_str (&cmdparams, "ds", NULL);
  process = cmdparams_get_str (&cmdparams, "process", NULL);
  format = cmdparams_get_str (&cmdparams, "format", NULL);
  method = cmdparams_get_str (&cmdparams, "method", NULL);
  protocol = cmdparams_get_str (&cmdparams, "protocol", NULL);
  requestor = cmdparams_get_str (&cmdparams, "requestor", NULL);
  notify = cmdparams_get_str (&cmdparams, "notify", NULL);
  shipto = cmdparams_get_str (&cmdparams, "shipto", NULL);
  requestorid = cmdparams_get_str (&cmdparams, "requestorid", NULL);

  dojson = strcmp(format, "json") == 0;
  dotxt = strcmp(format, "txt") == 0;
  dohtml = strcmp(format, "html") == 0;
  doxml = strcmp(format, "xml") == 0;

  export_series = EXPORT_SERIES;
  /*  op == exp_request  */
  if (strcmp(op,"exp_request") == 0) 
    {
    int status=0;
    int segcount = 0;
    int irec;
    int all_online = 1;
    DRMS_RecordSet_t *rs;
    export_series = EXPORT_SERIES_NEW;
    size=0;
    rs = drms_open_records(drms_env, in, &status);
    if (!rs)
	JSONDIE("Can not open RecordSet");

    // Do survey of recordset
    all_online = 1;
    for (irec=0; irec < rs->n; irec++) 
      {
      char recpath[DRMS_MAXPATHLEN];
      DRMS_Record_t *rec = rs->records[irec];
      drms_record_directory(rec,recpath,0);
      if (strncmp(recpath, "/SUM", 4) != 0)
          all_online = 0;
      size += drms_record_size(rec);
      segcount += drms_record_numsegments(rec);
      }

    // Do quick export if possible
    if (strcmp(method,"url_quick")==0 && strcmp(protocol,"as-is")==0 && segcount == 0)
      JSONDIE("There are no files in this RecordSet");
    if (strcmp(method,"url_quick")==0 && strcmp(protocol,"as-is")==0 && all_online)
      {
      if (dojson)
        {
        char *json;
        char *strval;
        int count;
        json_t *jroot = json_new_object();
        count = quick_export_rs(jroot, rs, 0, size); // add count, size, and array data of names and paths
        json_insert_pair_into_object(jroot, "requestid", json_new_string(""));
        free(strval);
        strval = string_to_json(method);
        json_insert_pair_into_object(jroot, "method", json_new_string(strval));
        free(strval);
        strval = string_to_json(protocol);
        json_insert_pair_into_object(jroot, "protocol", json_new_string(strval));
        free(strval);
        json_insert_pair_into_object(jroot, "wait", json_new_number("0"));
        json_insert_pair_into_object(jroot, "status", json_new_number("0"));
        json_tree_to_string(jroot,&json);

        printf("Content-type: application/json\n\n");
        printf("%s\n",json);
        fflush(stdout);

        free(json);
        return(0);
        }  
      else
        JSONDIE("format not implemented yet");
      }

    // Must do full export processing

    // Get RequestID
    {
    FILE *fp = popen("/home/phil/cvs/JSOC/bin/linux_ia32/GetJsocRequestID", "r");
    if (fscanf(fp, "%s", new_requestid) != 1)
      JSONDIE("Cant get new RequestID");
    pclose(fp);
    requestid = new_requestid;
    }

    now = timenow();

    // Add Requestor info to jsoc.export_user series 
    // Can not watch for new information since can not read this series.
    //   start by looking up requestor 
    if (strcmp(requestor, "Not Specified") != 0)
      {
      DRMS_Record_t *requestor_rec;
#ifdef IN_MY_DREAMS
      DRMS_RecordSet_t *requestor_rs;
      char requestorquery[2000];
      sprintf(requestorquery, "%s[? Requestor = '%s' ?]", EXPORT_USER, requestor);
      requestor_rs = drms_open_records(drms_env, requestorquery, &status);
      if (!requestor_rs)
        JSONDIE("Cant find requestor info series");
      if (requestor_rs->n == 0)
        { // First request for this user
        drms_close_records(requestor_rs, DRMS_FREE_RECORD);
#endif
        requestor_rec = drms_create_record(drms_env, EXPORT_USER, DRMS_PERMANENT, &status);
        if (!requestor_rec)
          JSONDIE("Cant create new user info record");
        requestorid = requestor_rec->recnum;
        drms_setkey_int(requestor_rec, "RequestorID", requestorid);
        drms_setkey_string(requestor_rec, "Requestor", requestor);
        drms_setkey_string(requestor_rec, "Notify", notify);
        drms_setkey_string(requestor_rec, "ShipTo", shipto);
        drms_setkey_time(requestor_rec, "FirstTime", now);
        drms_setkey_time(requestor_rec, "UpdateTime", now);
        drms_close_record(requestor_rec, DRMS_INSERT_RECORD);
#ifdef IN_MY_DREAMS
        }
      else
        { // returning user, look for updated info
        // WARNING ignore adding new info for now - XXXXXX must fix this later
        requestor_rec = requestor_rs->records[0];
        requestorid = drms_getkey_int(requestor_rec, "RequestorID", NULL);
        drms_close_records(requestor_rs, DRMS_FREE_RECORD);
        }
#endif
      }
    else
      requestorid = 0;

    // Create new record in export control series
    // This will be copied into the cluster-side series on first use.
    export_log = drms_create_record(drms_env, export_series, DRMS_PERMANENT, &status);
    if (!export_log)
      JSONDIE("Cant create new export control record");

    drms_setkey_string(export_log, "RequestID", requestid);
    drms_setkey_string(export_log, "DataSet", in);
    drms_setkey_string(export_log, "Processing", process);
    drms_setkey_string(export_log, "Protocol", protocol);
    drms_setkey_string(export_log, "Method", method);
    drms_setkey_string(export_log, "Format", format);
    drms_setkey_time(export_log, "ReqTime", now);
    drms_setkey_time(export_log, "EstTime", now+10); // Crude guess for now
    drms_setkey_longlong(export_log, "Size", size);
    drms_setkey_int(export_log, "Status", 2);
    drms_setkey_int(export_log, "Requestor", requestorid);
    drms_close_record(export_log, DRMS_INSERT_RECORD);
    }

  // Now report back to the requestor by dropping into the code for status request.

  // op = exp_status

  if (strcmp(requestid, "Not Specified") == 0)
    JSONDIE("RequestID must be provided");

  sprintf(status_query, "%s[%s]", export_series, requestid);
  exports = drms_open_records(drms_env, status_query, &status);
  if (!exports)
    JSONDIE("Cant locate export series");
  if (exports->n < 1)
    JSONDIE("Cant locate export request");
  export_log = exports->records[0];

  status     = drms_getkey_int(export_log, "Status", NULL);
  in         = drms_getkey_string(export_log, "DataSet", NULL);
  process = drms_getkey_string(export_log, "Processing", NULL);
  protocol   = drms_getkey_string(export_log, "Protocol", NULL);
  method     = drms_getkey_string(export_log, "Method", NULL);
  format     = drms_getkey_string(export_log, "Format", NULL);
  reqtime    = drms_getkey_time(export_log, "ReqTime", NULL);
  esttime    = drms_getkey_time(export_log, "EstTime", NULL); // Crude guess for now
  size       = drms_getkey_int(export_log, "Size", NULL);
  requestorid = drms_getkey_int(export_log, "Requestor", NULL);

  // Do special actions on status
  switch (status)
    {
    case 0:
            errorreply = NULL;
            waittime = 0;
            break;
    case 1:
            errorreply = NULL;
            waittime = esttime - timenow();
            break;
    case 2:
            errorreply = NULL;
            waittime = esttime - timenow();
            break;
    case 3:
            errorreply = "Request too large";
            waittime = 999999;
            break;
    case 4:
            waittime = 999999;
            errorreply = "RecordSet specified does not exist";
            break;
    case 5:
            waittime = 999999;
            errorreply = "Request was completed but is now deleted, 7 day limit exceeded";
            break;
    default:
      JSONDIE("Illegal status in export record");
    }

  // Return status information to user
  if (dojson)
    {
    char *json;
    char *strval;
    char numval[100];
    json_t *jsonval;
    json_t *jroot = json_new_object();

    if (status > 0)
      {
      sprintf(numval, "%d", status);
      json_insert_pair_into_object(jroot, "status", json_new_number(numval));
      strval = string_to_json(requestid);
      json_insert_pair_into_object(jroot, "requestid", json_new_string(strval));
      free(strval);
      strval = string_to_json(method);
      json_insert_pair_into_object(jroot, "method", json_new_string(strval));
      free(strval);
      strval = string_to_json(protocol);
      json_insert_pair_into_object(jroot, "protocol", json_new_string(strval));
      free(strval);
      sprintf(numval, "%1.0lf", waittime);
      json_insert_pair_into_object(jroot, "wait", json_new_number(numval));
      sprintf(numval, "%ld", size);
      json_insert_pair_into_object(jroot, "size", json_new_number(numval));
      if (errorreply)
        {
        strval = string_to_json(errorreply);
        json_insert_pair_into_object(jroot, "error", json_new_string(strval));
        free(strval);
        }
      if (status > 2)
        {
        strval = string_to_json("jsoc_help@jsoc.stanford.edu");
        json_insert_pair_into_object(jroot, "contact", json_new_string(strval));
        free(strval);
        }
      json_tree_to_string(jroot,&json);
      printf("Content-type: application/json\n\n");
      printf("%s\n",json);
      fflush(stdout);
      return(0);
      }
    else  // (status == 0)
      {
      // this what the user has been waiting for.  The export record segment dir should
      // contain a file containing the json to be returned to the user.  The dir will be returned as well as the file.
      char logpath[DRMS_MAXPATHLEN];
      FILE *fp;
      int c;
      drms_record_directory(export_log, logpath, 0);
      strncat(logpath, "/", DRMS_MAXPATHLEN);
      strncat(logpath, "index.json", DRMS_MAXPATHLEN);
      fp = fopen(logpath, "r");
      if (!fp)
        JSONDIE("Export should be complete but return index.json file not found");

      printf("Content-type: application/json\n\n");
      while ((c = fgetc(fp)) != EOF)
        putchar(c);
      fclose(fp);
      fflush(stdout);
      return(0);
      }
    }
  else
    JSONDIE("format not implemented yet");

  return(0);
}

