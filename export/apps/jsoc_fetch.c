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
#include "printk.h"


#include <time.h>


#define EXPORT_SERIES "jsoc.export"
#define EXPORT_SERIES_NEW "jsoc.export_new"
#define EXPORT_USER "jsoc.export_user"

int dojson, dotxt, dohtml, doxml;

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

ModuleArgs_t module_args[] =
{ 
  {ARG_STRING, "op", "Not Specified", "<Operation>"},
  {ARG_STRING, "ds", "Not Specified", "<record_set query>"},
  {ARG_STRING, "seg", "Not Specified", "<record_set segment list>"},
  {ARG_STRING, "sunum", "Not Specified", "<sunum list for SU exports>"},
  {ARG_STRING, "requestid", "Not Specified", "JSOC export request identifier"},
  {ARG_STRING, "process", "Not Specified", "string containing program and arguments"},
  {ARG_STRING, "requestor", "Not Specified", "name of requestor"},
  {ARG_STRING, "notify", "Not Specified", "email address of requestor"},
  {ARG_STRING, "shipto", "Not Specified", "mail address of requestor"},
  {ARG_STRING, "protocol", "as-is", "exported file protocol"},
  {ARG_STRING, "filenamefmt", "{seriesname}.{recnum:%d}.{segment}", "exported file filename format"},
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
	"seg=<list of segment names to append to dataset spec\n"
	"requestid=JSOC export request identifier\n"
	"process=string containing program and arguments\n"
	"requestor=name of requestor\n"
	"notify=email address of requestor\n"
	"shipto=mail address of requestor\n"
	"protocol=exported file protocol\n"
	"filenamefmt=exported file filename format\n"
	"format=return content type\n"
	"method=return method\n"
	"h=help - show usage\n"
	"QUERY_STRING=AJAX query from the web"
	);
    return(1);
    }
  return (0);
  }

char *json_text_to_string(char *in)
  {
  char *o, *new = (char *)malloc(strlen(in)+1);
  char *i;
  for (i=in, o=new; *i; )
    {
    if (*i == '\\')
      {
      i++;
      if (*i == '/' || *i == '"' || *i == '\\' )
        { *o++ = *i++; continue;}
      else if (*i == 'b')
	{ *o++ = '\b'; i++; continue;}
      else if (*i == 'f')
	{ *o++ = '\f'; i++; continue;}
      else if (*i == 'r')
	{ *o++ = '\r'; i++; continue;}
      else if (*i == 't')
	{ *o++ = '\t'; i++; continue;}
      else if (*i == 'n')
	{ *o++ = '\n'; i++; continue;}
// need to do uXXXX too!
      }
    *o++ = *i++;
    }
  *o = '\0';
  return(new);
  }

char *string_to_json(char *in)
  {
  char *new;
  new = json_escape(in);
  return(new);
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
    HIterator_t *hit = NULL;
    rec = rs->records[i];
    drms_sprint_rec_query(query, rec);
    while (seg = drms_record_nextseg(rec, &hit))
      {
      DRMS_Record_t *segrec;
      json_t *recobj = json_new_object();
      char *jsonstr;
      segrec = seg->record;
      count += 1;
      strcpy(record, query);
      strcat(record, "{");
      strcat(record, seg->info->name);
      strcat(record, "}");
      drms_record_directory(segrec, segpath, online); // just to insure SUM_get done
      drms_segment_filename(seg, segpath);
      jsonstr = string_to_json(record);
      json_insert_pair_into_object(recobj, "record", json_new_string(jsonstr));
      free(jsonstr);
      jsonstr = string_to_json(segpath);
      json_insert_pair_into_object(recobj, "filename", json_new_string(jsonstr));
      free(jsonstr);
      json_insert_child(data, recobj);
      }
    free(hit);
    }
  if (jroot) // i.e. if dojson, else will be NULL for the dotxt case.
    {
    sprintf(numval, "%ld", count);
    json_insert_pair_into_object(jroot, "count", json_new_number(numval));
    sprintf(numval, "%ld", size);
    json_insert_pair_into_object(jroot, "size", json_new_number(numval));
    json_insert_pair_into_object(jroot, "dir", json_new_string(""));
    json_insert_pair_into_object(jroot, "data", data);
    }
  else
    {
    json_t *recobj = data->child;
    printf("count=%ld\n", count);
    printf("size=%ld\n", size);
    printf("dir=/\n");
    printf("# DATA\n");
    while (recobj)
      {
      char *ascii_query, *ascii_path;
      json_t *record = recobj->child;
      json_t *filename = record->next;
      json_t *recquery = record->child;
      json_t *pathname = filename->child;
      ascii_query = json_text_to_string(recquery->text);
      ascii_path = json_text_to_string(pathname->text);
      printf("%s\t%s\n",ascii_query, ascii_path);
      free(ascii_query);
      free(ascii_path);
      recobj = recobj->next;
      }
    }
  return(count);
  }

TIME timenow()
  {
  TIME UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
  TIME now = (double)time(NULL) + UNIX_epoch;
  return(now);
  }

SUM_t *my_sum=NULL;

SUM_info_t *drms_get_suinfo(long long sunum)
  {
  int status;
  if (my_sum && my_sum->sinfo->sunum == sunum)
    return(my_sum->sinfo);
  if (!my_sum)
    {
    if ((my_sum = SUM_open(NULL, NULL, printkerr)) == NULL)
      {
      printkerr("drms_open: Failed to connect to SUMS.\n");
      return(NULL);
      }
    }
  if (status = SUM_info(my_sum, sunum, printkerr))
    {
    printkerr("Fail on SUM_info, status=%d\n", status);
    return(NULL);
    }

  return(my_sum->sinfo);
  }

#ifdef NEVER
// use at end of any section that exits Doit after using SUM_info
  if (my_sum)
    SUM_close(my_sum,printkerr);
#endif

#define JSONDIE(msg) {die(dojson,msg,"");return(1);}
#define JSONDIE2(msg,info) {die(dojson,msg,info);return(1);}

die(int dojson, char *msg, char *info)
  {
  char *msgjson;
  char errval[10];
  char *json;
  char message[10000];
  json_t *jroot = json_new_object();
if (DEBUG) fprintf(stderr,"%s%s\n",msg,info);
  strcpy(message,msg); 
  strcat(message,info); 
  if (dojson)
    {
    msgjson = string_to_json(message);
    json_insert_pair_into_object(jroot, "status", json_new_number("4"));
    json_insert_pair_into_object(jroot, "error", json_new_string(msgjson));
    json_tree_to_string(jroot,&json);
    printf("Content-type: application/json\n\n");	
    printf("%s\n",json);
    }
  else
    {
    printf("Content-type: text/plain\n\n");
    printf("status=4\nerror=%s\n", message);
    }
  fflush(stdout);
  if (my_sum)
    SUM_close(my_sum,printkerr);
  return(1);
  }

send_file(DRMS_Record_t *rec, int segno)
  {
  DRMS_Segment_t *seg = drms_segment_lookupnum(rec, 0);
  char path[DRMS_MAXPATHLEN];
  char sudir[DRMS_MAXPATHLEN];
  FILE *fp;
  int b;

  drms_record_directory(rec,sudir,0);
  strcpy(path, "/web/jsoc/htdocs");
  strncat(path, sudir, DRMS_MAXPATHLEN);
  strncat(path, "/", DRMS_MAXPATHLEN);
  strncat(path, seg->filename, DRMS_MAXPATHLEN);
fprintf(stderr,"path: %s\n",path);
  fp = fopen(path, "r");
  if (!fp)
    JSONDIE2("Can not open file for export: ",path);
  switch (seg->info->protocol)
    {
    case DRMS_FITS:
        printf("Content-Type: application/fits\n\n");	
        break;
    default:
        printf("Content-Type: application/binary\n\n");	
//  Content-Type: application/fits
//Length: 152,640 (149K) [application/fits]

        // look at file extension to guess file type here and print content type line
    }
  while ((b = fgetc(fp))!= EOF)
    fputc(b, stdout);
  fclose(fp);
  fflush(stdout);
  }

/* Module main function. */
int DoIt(void)
  {
						/* Get command line arguments */
  char *op;
  char *in;
  char *seglist;
  char *requestid;
  char *process;
  char *requestor;
  long requestorid;
  char *notify;
  char *format;
  char *shipto;
  char *method;
  char *protocol;
  char *filenamefmt;
  char *errorreply;
  char *sunumlist;
  long long size;
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

  dojson=1; dotxt=0; dohtml=0; doxml=0;

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
  sunumlist = cmdparams_get_str (&cmdparams, "sunum", NULL);
  seglist = cmdparams_get_str (&cmdparams, "seg", NULL);
  process = cmdparams_get_str (&cmdparams, "process", NULL);
  format = cmdparams_get_str (&cmdparams, "format", NULL);
  method = cmdparams_get_str (&cmdparams, "method", NULL);
  protocol = cmdparams_get_str (&cmdparams, "protocol", NULL);
  filenamefmt = cmdparams_get_str (&cmdparams, "filenamefmt", NULL);
  requestor = cmdparams_get_str (&cmdparams, "requestor", NULL);
  notify = cmdparams_get_str (&cmdparams, "notify", NULL);
  shipto = cmdparams_get_str (&cmdparams, "shipto", NULL);
  requestorid = cmdparams_get_str (&cmdparams, "requestorid", NULL);

  dojson = strcmp(format, "json") == 0;
  dotxt = strcmp(format, "txt") == 0;
  dohtml = strcmp(format, "html") == 0;
  doxml = strcmp(format, "xml") == 0;

  export_series = EXPORT_SERIES;

  /*  op == exp_su - export Storage Units */
  if (strcmp(op,"exp_su") == 0)
    {
    char *this_sunum, *sunumlist, *sunumlistptr;
    long long sunums[DRMS_MAXQUERYLEN/8];  // should be enough!
    char *paths[DRMS_MAXQUERYLEN/8];
    char *series[DRMS_MAXQUERYLEN/8];
    long long sunum;
    int count;
    int status=0;
    int all_online;

    export_series = EXPORT_SERIES_NEW;
    // Do survey of sunum list
    size=0;
    all_online = 1;
    count = 0;
    sunumlist = strdup(in);

    while (this_sunum = strtok_r(sunumlist, ",", &sunumlistptr))
      {
      SUM_info_t *sinfo;
      TIME expire;
      sunum = atoll(this_sunum);
      sunumlist = NULL;
      sinfo = drms_get_suinfo(sunum);
      if (!sinfo)
        JSONDIE("Invalid sunum, SUM_info call failed");
      size += sinfo->bytes;
      if (strcmp(sinfo->online_status,"Y")==0)
        {
        int y,m,d,hr,mn;
        char sutime[50];
        sscanf(sinfo->effective_date,"%4d%2d%2d%2d%2d", &y,&m,&d,&hr,&mn);
        sprintf(sutime, "%4d.%02d.%02d_%02d:%02d", y,m,d,hr,mn);
        expire = (sscan_time(sutime) - now)/86400.0;
        }
      if (strcmp(sinfo->online_status,"N")==0 || expire < 3)
        {  // need to stage or reset retention time
        all_online = 0;
        }
      else
        {
        sunums[count] = sunum;
        paths[count] = strdup(sinfo->online_loc);
        series[count] = strdup(sinfo->owning_series);
        }
      count += 1;
      }
    if (count==0)
      JSONDIE("There are no files in this RecordSet");

    // Do quick export if possible
    if (strcmp(method,"url_quick")==0 && strcmp(protocol,"as-is")==0  && all_online)
      {
      if (dojson)
        {
        int i;
        char *json;
        char *strval;
        char numval[50];
        json_t *jroot = json_new_object();
        json_t *data;
        data = json_new_array();
        for (i=0; i < count; i++)
          {
          json_t *suobj = json_new_object();
          char *jsonstr;
          char numval[40];
          sprintf(numval,"%lld",sunums[i]);
          jsonstr = string_to_json(numval); // send as string in case long long fails
          json_insert_pair_into_object(jroot, "sunum", json_new_string(jsonstr));
          free(jsonstr);
          jsonstr = string_to_json(series[i]);
          json_insert_pair_into_object(suobj, "series", json_new_string(jsonstr));
          free(jsonstr);
          jsonstr = string_to_json(paths[i]);
          json_insert_pair_into_object(suobj, "path", json_new_string(jsonstr));
          free(jsonstr);
          json_insert_child(data, suobj);
          }
        sprintf(numval, "%ld", count);
        json_insert_pair_into_object(jroot, "count", json_new_number(numval));
        sprintf(numval, "%lld", size);
        json_insert_pair_into_object(jroot, "size", json_new_number(numval));
        json_insert_pair_into_object(jroot, "dir", json_new_string(""));
        json_insert_pair_into_object(jroot, "data", data);
        json_insert_pair_into_object(jroot, "requestid", json_new_string(""));
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
        }  
      else
        {
        int i;
        printf("Content-type: text/plain\n\n");
        printf("# JSOC Quick Data Export of as-is files.\n");
        printf("status=0\n");
        printf("requestid=\"Not Specified\"\n");
        printf("method=%s\n", method);
        printf("protocol=%s\n", protocol);
        printf("wait=0\n");
        printf("count=%d\n", count);
        printf("size=%lld\n", size);
        printf("dir=/\n");
        printf("# DATA\n");
        for (i=0; i<count; i++)
          printf("%lld\t%s\t%s\n",sunums[i],series[i],paths[i]);
        }
      if (my_sum)
        SUM_close(my_sum,printkerr);
      return(0);
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
#ifdef SHOULD_BE_HERE
check for requestor to be valid remote DRMS site
#else // for now
      requestorid = 0;
#endif
      }
    else
      requestorid = 0;

    // FORCE process to be su_export
    process = "su_export";

    // Create new record in export control series
    // This will be copied into the cluster-side series on first use.
    export_log = drms_create_record(drms_env, export_series, DRMS_PERMANENT, &status);
    if (!export_log)
      JSONDIE("Cant create new export control record");
    drms_setkey_string(export_log, "RequestID", requestid);
    drms_setkey_string(export_log, "DataSet", in);
    drms_setkey_string(export_log, "Processing", process);
    drms_setkey_string(export_log, "Protocol", protocol);
    drms_setkey_string(export_log, "FilenameFmt", filenamefmt);
    drms_setkey_string(export_log, "Method", method);
    drms_setkey_string(export_log, "Format", format);
    drms_setkey_time(export_log, "ReqTime", now);
    drms_setkey_time(export_log, "EstTime", now+10); // Crude guess for now
    drms_setkey_longlong(export_log, "Size", size);
    drms_setkey_int(export_log, "Status", 2);
    drms_setkey_int(export_log, "Requestor", requestorid);
    drms_close_record(export_log, DRMS_INSERT_RECORD);
    } // end of exp_su
  /*  op == exp_request  */
  else if (strcmp(op,"exp_request") == 0) 
    {
    int status=0;
    int segcount = 0;
    int irec;
    int all_online = 1;
    char dsquery[DRMS_MAXQUERYLEN];
    char *p;
    DRMS_RecordSet_t *rs;
    export_series = EXPORT_SERIES_NEW;
    size=0;
    strncpy(dsquery,in,DRMS_MAXQUERYLEN);
    if (index(dsquery,'[') == NULL)
      strcat(dsquery,"[]");
    if (strcmp(seglist,"Not Specified") != 0)
      {
      if (index(dsquery,'{') != NULL)
        JSONDIE("Can not give segment list both in key and explicitly in recordset.");
      strncat(dsquery, "{", DRMS_MAXQUERYLEN);
      strncat(dsquery, seglist, DRMS_MAXQUERYLEN);
      strncat(dsquery, "}", DRMS_MAXQUERYLEN);
      }
    if ((p=index(dsquery,'{')) != NULL && strncmp(p+1, "**ALL**", 7) == 0)
      *p = '\0';
    rs = drms_open_records(drms_env, dsquery, &status);
    if (!rs)
	JSONDIE2("Can not open RecordSet: ",dsquery);

    // Do survey of recordset
    all_online = 1;
    for (irec=0; irec < rs->n; irec++) 
      {
      // Must check each segment since some may be linked and offline.
      DRMS_Record_t *rec = rs->records[irec];
      DRMS_Segment_t *seg;
      HIterator_t *segp = NULL;
      while (seg = drms_record_nextseg(rec, &segp))
        {
        DRMS_Record_t *segrec = seg->record;
        SUM_info_t *sinfo = drms_get_suinfo(segrec->sunum);
        if (!sinfo)
	  JSONDIE2("Bad sunum in a record in RecordSet: ", dsquery);
	if (strcmp(sinfo->online_status,"N") == 0)
            all_online = 0;
        else
            {
            struct stat buf;
	    char path[DRMS_MAXPATHLEN];
	    drms_record_directory(segrec, path, 0);

            drms_segment_filename(seg, path);
            if (stat(path, &buf) != 0)
	      JSONDIE2("Bad path (stat failed) in a record in RecordSet: ", dsquery);
            size += buf.st_size;
            segcount += 1;
            }
         }
      if (segp)
        free(segp); 
      }
    if (my_sum)
      SUM_close(my_sum,printkerr);

    // Do quick export if possible
    if (strcmp(method,"url_quick")==0 && (strcmp(protocol,"as-is")==0 || strcmp(protocol,"su")==0) && segcount == 0)
      JSONDIE("There are no files in this RecordSet");
    if (strcmp(method,"url_quick")==0 && (strcmp(protocol,"as-is")==0 || strcmp(protocol,"su")==0) && all_online)
      {
      if (0 && segcount == 1) // If only one file then do immediate delivery of that file.
        {
	send_file(rs->records[0], 0);
        }
      else if (dojson)
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
        }  
      else
	{
	printf("Content-type: text/plain\n\n");
	printf("# JSOC Quick Data Export of as-is files.\n");
	printf("status=0\n");
	printf("requestid=\"Not Specified\"\n");
	printf("method=%s\n", method);
	printf("protocol=%s\n", protocol);
	printf("wait=0\n");
        quick_export_rs(NULL, rs, 0, size); // add count, size, and array data of names and paths
	}
      return(0);
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
    drms_setkey_string(export_log, "DataSet", dsquery);
    drms_setkey_string(export_log, "Processing", process);
    drms_setkey_string(export_log, "Protocol", protocol);
    drms_setkey_string(export_log, "FilenameFmt", filenamefmt);
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
  // This is entry point for status request and tail of work for exp_request and exp_su
  // If data was as-is and online and url_quick the exit will have happened above.

  // op = exp_status

  if (strcmp(requestid, "Not Specified") == 0)
    JSONDIE("RequestID must be provided");

  sprintf(status_query, "%s[%s]", export_series, requestid);
  exports = drms_open_records(drms_env, status_query, &status);
  if (!exports)
    JSONDIE2("Cant locate export series: ", status_query);
  if (exports->n < 1)
    JSONDIE2("Cant locate export request: ", status_query);
  export_log = exports->records[0];

  status     = drms_getkey_int(export_log, "Status", NULL);
  in         = drms_getkey_string(export_log, "DataSet", NULL);
  process = drms_getkey_string(export_log, "Processing", NULL);
  protocol   = drms_getkey_string(export_log, "Protocol", NULL);
  filenamefmt = drms_getkey_string(export_log, "FilenameFmt", NULL);
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
  if (1)
    {
    char *json;
    char *strval;
    char numval[100];
    json_t *jsonval;
    json_t *jroot=NULL;

    if (status > 0)
      {
      if (dojson)
	{
        jroot = json_new_object();
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
        }
      else
        {
	printf("Content-type: text/plain\n\n");
        printf("# JSOC Data Export Not Ready.\n");
        printf("status=%d\n", status);
        printf("requestid=%s\n", requestid);
        printf("method=%s\n", method);
        printf("protocol=%s\n", protocol);
        printf("wait=%d\n",waittime);
	printf("size=%ld\n",size);
        if (errorreply)
	  printf("error=\"%s\"\n", errorreply);
	if (status > 2)
	  printf("contact=jsoc_help@jsoc.stanford.edu\n");
        }
      fflush(stdout);
      }
    else  // (status == 0)
      {
      // this what the user has been waiting for.  The export record segment dir should
      // contain a file containing the json to be returned to the user.  The dir will be returned as well as the file.
      char logpath[DRMS_MAXPATHLEN];
      FILE *fp;
      int c;
      char *indexfile = (dojson ? "index.json" : "index.txt");
      jroot = json_new_object();
      drms_record_directory(export_log, logpath, 0);
      strncat(logpath, "/", DRMS_MAXPATHLEN);
      strncat(logpath, indexfile, DRMS_MAXPATHLEN);
      fp = fopen(logpath, "r");
      if (!fp)
        JSONDIE2("Export should be complete but return %s file not found", indexfile);

      if (dojson)
        printf("Content-type: application/json\n\n");
      else
	printf("Content-type: text/plain\n\n");
      while ((c = fgetc(fp)) != EOF)
        putchar(c);
      fclose(fp);
      fflush(stdout);
      }
    }
  return(0);
  }
