#define DEBUG 1
#define DEBUG 0

/*
 *  jsoc_info - prints keyword information and/or file path for given recordset
 *
 *  new version of original show_keys expanded to have more than just keyword info.
 *
 *  Bugs:
 *	Fails (with a segmentation fault) if there are no records in the
 *	  requested series
 *	Fails (with a segmentation fault) if called with -p flag or with a
 *	  value for file and there are no data segments associated with the
 *	  requested series
 */

/**
\defgroup jsoc_info jsoc_info

Prints keyword, segment, and other information and/or file path for given recordset.

\ref jsoc_info can list the keyword names and values, and the segment
names and file names (full paths) for each record in a record set. It
can also list the full path to the record direcory in SUMS, which
contains the segment files. Exactly what information gets printed is
controlled by command-line flags (see below). The \a -k flag controls
the format of the output.  If it is set, then the output is in table
format, with a header row showing the keyword names.  Otherwise,
keyword name=value pairs are listed one per line.  If the \a -a flag
is set, \ref jsoc_info lists the names of all series keywords, prime
keywords, and segments, and exits.  Otherwise, it prints keyword and
segment information as specified by the other flags and arguments.  If
the \a -p flag is set and \a seglist is specified, then the full paths
for the segment files will be displayed. If the \a -p flag is set, but
\a seglist is not specified, then only the full path to the record's
storage unit will be displayed.

The number of records for which information will be printed must be
specified, either by supplying a \a record_set string that selects a
subset of records from a series, or by supplying the \a n=nrecords
argument, which indicates the number of records.

\par Synopsis:

\code
jsoc_info [-ajklpqrsDRIVER_FLAGS] ds=<record_set> [n=<nrecords>] [key=<keylist>] [seg=<seglist>]
\endcode

\b Example:
To show the storage-unit paths for a maximum of 10
records:
\code
  jsoc_info -p ds=su_arta.TestStoreFile n=10
\endcode

\b Example:
To show information, in non-table format, for all keywords,
plus the segment named file_seg, for a maximum of 10 records:
\code
  jsoc_info ds=su_arta.TestStoreFile -akr n=10 seg=file_seg
\endcode

\par Flags:
\c -a: Show all keyword names and values for each  record  specified
by \a record_set  or \a nrecords.  \a -a takes precedence over \a
keylist.  
\par
\c -j: List the names of all series keywords, prime keywords, and
segments, and links in jsd format and exit. 
\par
\c -k: List keyword name=value pairs, one per line. Otherwise print
all keyword values on a single line and print a header line containing
the keyword names (table format).
\par
\c -l: List the names of all series keywords, prime keywords,  and
segments, and exit. 
\par
\c -p: Include in the output the full storage-unit path for each record
\par
\c -q: Quiet - omit the header line listing keyword names if the -k
flag is set
\par
\c -r:  Include the record number in the output
\par
\c -s:  Include statistics of series in the output

\par Driver flags: 
\ref jsoc_main

\param record_set
A series name followed by an optional record-set specification (i.e.,
\a seriesname[RecordSet_filter]). Causes selection of a subset of
records in the series. This argument is required, and if no record-set
filter is specified, then \a n=nrecords must be present.

\param nrecords
\a nrecords specifies the maximum number of records for which
information is printed.  If \a nrecords < 0, \ref jsoc_info displays
information for the last \a nrecords records in the record set. If
\a nrecords > 0, \ref jsoc_info displays information for the first
\a nrecords records in the record set. If \a record_set contains a
record set filter, then \a nrecords can reduce the total number of
records for which information is displayed.

\param keylist
Comma-separated list of keyword names. For each keyword listed,
information will be displayed.  \a keylist is ignored in the case that
the \a -a flag is set.

\param seglist
Comma-separated list of segment names. For each segment listed, the
full path to the segment's file is displayed
(if the \a -p flag is set) or the file name of the
segment's file name is displayed (if the \a -p flag is unset).

\bug
The program will produce superflous and non-meaningful output if
called with the \a -p flag and \a seglist is provided on the command line.

\sa
retrieve_file drms_query describe_series

@{
*/
#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"
#include "json.h"
// #include "json.h"

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

/* drms_record_getlogdir */

/* Returns path of directory that contains any saved log information for the given record */
/* If log is offline, returns message, if log was not saved or otherwise not found returns NULL */
/* The returned char* should be freed after use. */

char *drms_record_getlogdir(DRMS_Record_t *rec)
  { 
  char *logpath;
  char query[DRMS_MAXQUERYLEN];
  DB_Text_Result_t *qres;

  sprintf(query,
    "select online_loc, online_status from sum_main a, %s.drms_session b where a.ds_index = b.sunum and b.sessionid=%lld",
    rec->sessionns, rec->sessionid);
  if ((qres = drms_query_txt(drms_env->session, query)) && qres->num_rows > 0)
    {
    if (qres->field[0][1][0] == 'Y')
      logpath = strdup(qres->field[0][0]);
    else
      logpath = strdup("Log Offline");
    }
  else
    logpath = NULL;
  db_free_text_result(qres);
  return logpath;
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
  {ARG_STRING, "key", "Not Specified", "<comma delimited keyword list>"},
  {ARG_STRING, "seg", "Not Specified", "<comma delimited segment list>"},
  {ARG_FLAG, "h", "0", "help - show usage"},
  {ARG_FLAG, "R", "0", "Show record query"},
  {ARG_FLAG, "z", "0", "emit JSON output"},
  {ARG_STRING, "QUERY_STRING", "Not Specified", "AJAX query from the web"},
  {ARG_END}
};

char *module_name = "jsoc_info";
/** @}*/
int nice_intro ()
  {
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\njsoc_info {-h} "
	"op=<command> ds=<recordset query> {n=0} {key=<keylist>} {seg=<segment_list>}\n"
        "  details are:\n"
	"op=<command> tell which ajax function to execute\n"
	"ds=<recordset query> as <series>{[record specifier]} - required\n"
	"key=<comma delimited keyword list>, for all use -a flag\n"
	"seg=<comma delimited segment list>\n"
	);
    return(1);
    }
  return (0);
  }

/* find first record in series that owns the given record */
DRMS_RecordSet_t *drms_find_rec_first(DRMS_Record_t *rec, int wantprime)
  {
  int iprime, nprime;
  int status;
  DRMS_RecordSet_t *rs;
  char query[DRMS_MAXQUERYLEN];
  strcpy(query, rec->seriesinfo->seriesname);
  nprime = rec->seriesinfo->pidx_num;
  if (wantprime && nprime > 0) 
    // only first prime key is used for now
     // for (iprime = 0; iprime < nprime; iprime++)
      strcat(query, "[#^]");
  else
    strcat(query, "[:#^]");
// fprintf(stderr,"test 1 query is %s\n",query);
  rs = drms_open_records(rec->env, query, &status);
// fprintf(stderr,"test 1 status is %d\n",status);
  return(rs);
  }

/* find last record in series that owns the given record */
DRMS_RecordSet_t *drms_find_rec_last(DRMS_Record_t *rec, int wantprime)
  {
  int iprime, nprime;
  int status;
  DRMS_RecordSet_t *rs;
  char query[DRMS_MAXQUERYLEN];
  strcpy(query, rec->seriesinfo->seriesname);
  nprime = rec->seriesinfo->pidx_num;
  if (wantprime && nprime > 0) 
    // only first prime key is used for now
     // for (iprime = 0; iprime < nprime; iprime++)
      strcat(query, "[#$]");
  else
    strcat(query, "[:#$]");
  rs = drms_open_records(rec->env, query, &status);
  return(rs);
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
      if (rec_key->info->type == DRMS_TYPE_TIME)
      {
	 drms_sprintfval_format(val, rec_key->info->type, &(rec_key->value), rec_key->info->unit, 0);
      }
      else
      {
	 drms_sprintfval_format(val, rec_key->info->type, &(rec_key->value), rec_key->info->format, 0);
      }
      strcat(text, "[");
      strcat(text, val);
      strcat(text, "]");
      }
    }
  else
    sprintf(text, "[:#%lld]",rec->recnum);
  return;
  }

char * string_to_json(char *in)
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


static void list_series_info(DRMS_Record_t *rec, json_t *jroot)
  {
  DRMS_Keyword_t *key;
  DRMS_Segment_t *seg;
  DRMS_Link_t *link;
  HIterator_t hit;
  int iprime;
  char *notework;

if (DEBUG) fprintf(stderr,"   starting Primekeys\n");
  /* show the prime index keywords */
  json_t *primekeys = json_new_string("primekeys");
  json_t *primearray = json_new_array();
  int npkeys = rec->seriesinfo->pidx_num;
  if (npkeys > 0)
    {
    int i;
    for (i=0; i<npkeys; i++)
        {
        DRMS_Keyword_t *skey, *pkey;
	int status;
        skey = pkey = rec->seriesinfo->pidx_keywords[i];
	if (pkey->info->recscope > 1)
            pkey = drms_keyword_slotfromindex(pkey);
        json_insert_child(primearray, json_new_string(pkey->info->name));
        }
    }
  else
    json_insert_child(primearray, json_new_null());
  json_insert_child(primekeys, primearray);
  json_insert_child(jroot,primekeys);

if (DEBUG) fprintf(stderr,"   starting DBindex\n");
  /* show DB index keywords */
  json_t *dbindex = json_new_string("dbindex");
  json_t *indexarray = json_new_array();
  if (rec->seriesinfo->dbidx_num > 0)
    {
    int i;
    for (i=0; i<rec->seriesinfo->dbidx_num; i++)
        json_insert_child(indexarray, json_new_string((rec->seriesinfo->dbidx_keywords[i])->info->name));
    }
  else
    json_insert_child(indexarray, json_new_null());
  json_insert_child(dbindex, indexarray);
  json_insert_child(jroot,dbindex);

if (DEBUG) fprintf(stderr,"   starting keywords\n");
  /* show all keywords */
  json_t *allkeys = json_new_string("keywords");
  json_t *keyarray = json_new_array();
  hiter_new (&hit, &rec->keywords);
  while ((key = (DRMS_Keyword_t *)hiter_getnext (&hit)))
    {
    json_t *keyinfo = json_new_object();
    json_t *keytype;
if (DEBUG) fprintf(stderr,"      starting %s\n",key->info->name);
    json_insert_pair_into_object(keyinfo, "name", json_new_string(key->info->name));
    if (key->info->islink)
	keytype = json_new_string("link");
    else
	keytype = json_new_string(drms_type_names[key->info->type]);
    json_insert_pair_into_object(keyinfo, "type", keytype);
if (DEBUG) fprintf(stderr,"      starting note %s\n",key->info->description);
    notework = string_to_json(key->info->description);
if (DEBUG) fprintf(stderr,"      continuing note %s\n",notework);
    json_insert_pair_into_object(keyinfo, "note", json_new_string(notework));
    free(notework);
    json_insert_child(keyarray, keyinfo);
    }
  json_insert_child(allkeys,keyarray);
  json_insert_child(jroot,allkeys);
  
if (DEBUG) fprintf(stderr,"   starting segments\n");
  /* show the segments */
  json_t *allsegs = json_new_string("segments");
  json_t *segarray = json_new_array();
  if (rec->segments.num_total)
    {
    hiter_new (&hit, &rec->segments);
    while ((seg = (DRMS_Segment_t *)hiter_getnext (&hit)))
      { /* segment name, units, protocol, dims, description */
      json_t *seginfo = json_new_object();
      json_t *keytype;
      int iaxis, naxis = seg->info->naxis;
      json_insert_pair_into_object(seginfo, "name", json_new_string(seg->info->name));
      if (seg->info->islink)
	    {
            char linkinfo[DRMS_MAXNAMELEN+10];
	    sprintf(linkinfo, "link via %s", seg->info->linkname);
            json_insert_pair_into_object(seginfo, "units", json_new_null());
            json_insert_pair_into_object(seginfo, "protocol", json_new_string(linkinfo));
            json_insert_pair_into_object(seginfo, "dims", json_new_null());
	    }
	else
	    {
            char prot[DRMS_MAXNAMELEN];
            char diminfo[160];
            int iaxis;
            strcpy(prot, drms_prot2str(seg->info->protocol));
            json_insert_pair_into_object(seginfo, "units", json_new_string(seg->info->unit));
            json_insert_pair_into_object(seginfo, "protocol", json_new_string(prot));
	    diminfo[0] = '\0';
	    for (iaxis=0; iaxis<naxis; iaxis++)
	        {
	        if (iaxis != 0)
                    strcat(diminfo,"x");
		if (seg->info->scope == DRMS_VARDIM)
                    strcat(diminfo,"VAR");
		else
		    {
	            char size[10];
		    sprintf(size,"%d",seg->axis[iaxis]);
                    strcat(diminfo,size);
		    }
		}
            json_insert_pair_into_object(seginfo, "dims", json_new_string(diminfo));
            }
      notework = string_to_json(seg->info->description);
      json_insert_pair_into_object(seginfo, "note", json_new_string(notework));
      free(notework);
      json_insert_child(segarray, seginfo);
      }
    }
//  else
//    json_insert_child(segarray, json_new_null());
  json_insert_child(allsegs,segarray);
  json_insert_child(jroot,allsegs);

if (DEBUG) fprintf(stderr,"   starting links\n");
  /* show the links */
  json_t *alllinks = json_new_string("links");
  json_t *linkarray = json_new_array();
  if (rec->links.num_total)
    {
    hiter_new (&hit, &rec->links);
    while ((link = (DRMS_Link_t *)hiter_getnext (&hit)))
      {
      json_t *linkinfo = json_new_object();
      json_insert_pair_into_object(linkinfo, "name", json_new_string(link->info->name));
      json_insert_pair_into_object(linkinfo, "target", json_new_string(link->info->target_series));
      json_insert_pair_into_object(linkinfo, "kind", json_new_string(link->info->type == STATIC_LINK ? "STATIC" : "DYNAMIC"));
      notework = string_to_json(link->info->description);
      json_insert_pair_into_object(linkinfo, "note", json_new_string(notework));
      free(notework);
      json_insert_child(linkarray,linkinfo);
      }
    }
//  else
//    json_insert_child(linkarray,json_new_null());
  json_insert_child(alllinks,linkarray);
  json_insert_child(jroot,alllinks);
  return;
  }

void get_series_stats(DRMS_Record_t *rec, json_t *jroot)
  {
  DRMS_RecordSet_t *rs;
  int iprime, nprime;
  int status;
  char query[DRMS_MAXQUERYLEN];
  json_t *interval = json_new_object();

  nprime = rec->seriesinfo->pidx_num;
  if (nprime > 0)
    sprintf(query,"%s[#^]", rec->seriesinfo->seriesname);
  else
    sprintf(query,"%s[:#^]", rec->seriesinfo->seriesname);
  rs = drms_open_records(rec->env, query, &status);
  if (!rs || rs->n < 1)
    {
    json_insert_pair_into_object(interval, "FirstRecord", json_new_string("NA"));
    json_insert_pair_into_object(interval, "FirstRecnum", json_new_string("NA"));
    json_insert_pair_into_object(interval, "LastRecord", json_new_string("NA"));
    json_insert_pair_into_object(interval, "LastRecnum", json_new_string("NA"));
    json_insert_pair_into_object(interval, "MaxRecnum", json_new_number("0"));
    if (rs) drms_free_records(rs);
    json_insert_pair_into_object(jroot, "Interval", interval);
    return(0);
    }
  else
    {
    char recquery[DRMS_MAXQUERYLEN];
    char *jsonquery;
    char val[100];
    int status, count;
    drms_sprint_rec_query(recquery,rs->records[0]);
    jsonquery = string_to_json(recquery);
    status = json_insert_pair_into_object(interval, "FirstRecord", json_new_string(jsonquery));
if (status != JSON_OK) fprintf(stderr, "json_insert_pair_into_object, status=%d, text=%s\n",status,jsonquery);
    free(jsonquery);
    sprintf(val,"%d", rs->records[0]->recnum);
    json_insert_pair_into_object(interval, "FirstRecnum", json_new_number(val));
    drms_free_records(rs);
  
    if (nprime > 0)
      sprintf(query,"%s[#$]", rec->seriesinfo->seriesname);
    else
      sprintf(query,"%s[:#$]", rec->seriesinfo->seriesname);
    rs = drms_open_records(rec->env, query, &status);
    drms_sprint_rec_query(recquery,rs->records[0]);
    jsonquery = string_to_json(recquery);
    json_insert_pair_into_object(interval, "LastRecord", json_new_string(jsonquery));
    free(jsonquery);
    sprintf(val,"%d", rs->records[0]->recnum);
    json_insert_pair_into_object(interval, "LastRecnum", json_new_number(val));
    drms_free_records(rs);
 
    sprintf(query,"%s[:#$]", rec->seriesinfo->seriesname);
    rs = drms_open_records(rec->env, query, &status);
    sprintf(val,"%d", rs->records[0]->recnum);
    json_insert_pair_into_object(interval, "MaxRecnum", json_new_number(val));
    drms_free_records(rs);
    }
  json_insert_pair_into_object(jroot, "Interval", interval);
  return(0);
  }

#define JSONDIE(msg) \
  {	\
  char *msgjson;	\
  char errval[10];	\
  char *json;	\
  json_t *jroot = json_new_object();	\
if (DEBUG) fprintf(stderr,"%s\n",msg);	\
  msgjson = string_to_json(msg);	\
  json_insert_pair_into_object(jroot, "status", json_new_number("1"));	\
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
  char *keylist;
  char *seglist;
  char *web_query;
  int from_web, keys_listed, segs_listed;

/* JSON building places */

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
    // Force JSON for now
    cmdparams_set (&cmdparams,"z", "1");
    free(getstring);
    }

  op = cmdparams_get_str (&cmdparams, "op", NULL);
  in = cmdparams_get_str (&cmdparams, "ds", NULL);
  keylist = strdup (cmdparams_get_str (&cmdparams, "key", NULL));
  seglist = strdup (cmdparams_get_str (&cmdparams, "seg", NULL));
  keys_listed = strcmp (keylist, "Not Specified");
  segs_listed = strcmp (seglist, "Not Specified");

  /*  op == series_struct  */
  if (strcmp(op,"series_struct") == 0) 
    {
    char *p, *seriesname;
    json_t *jroot;
    char *json;
    DRMS_Record_t *rec;
    int status=0;
    /* Only want keyword info so get only the template record for drms series or first record for other data */
    seriesname = strdup (in);
    if ((p = index(seriesname,'['))) *p = '\0';
    rec = drms_template_record (drms_env, seriesname, &status);
    if (status)
      JSONDIE("series not found");
    jroot = json_new_object();
    list_series_info(rec, jroot);
    get_series_stats(rec, jroot);
    json_insert_pair_into_object(jroot, "status", json_new_number("0"));
    json_tree_to_string(jroot,&json);
    /* send the output json back to client */
    printf("Content-type: application/json\n\n");
    printf("%s\n",json);
    fflush(stdout);
    free(seriesname);
    return(0);
    }

  /*  op == rs_summary  */
  if (strcmp(op,"rs_summary") == 0) 
    {
    json_t *jroot;
    char *json;
    int count=0, status=0;
    char val[100];
    jroot = json_new_object();
    /* get series count */
    count = drms_count_records(drms_env, in, &status);
    if (status)
      JSONDIE("series not found");
    /* send the output json back to client */
    sprintf(val, "%d", count);
    json_insert_pair_into_object(jroot, "count", json_new_number(val));
    json_insert_pair_into_object(jroot, "status", json_new_number("0"));
    json_tree_to_string(jroot,&json);
    printf("Content-type: application/json\n\n");
    printf("%s\n",json);
    fflush(stdout);
    return(0);
    }

  /*  op == rs_list  */
  if (strcmp(op,"rs_list") == 0) 
    {
    int wantRecInfo = cmdparams_get_int(&cmdparams, "R", NULL);
    DRMS_RecordSet_t *recordset;
    DRMS_Record_t *rec;
    char *keys[1000];
    char *segs[1000];
    int ikey, nkeys = 0;
    int iseg, nsegs = 0;
    char count[100];
    json_t *jroot, **keyvals = NULL, **segvals = NULL, *recinfo;
    json_t *json_keywords = json_new_array();
    json_t *json_segments = json_new_array();
    char *json;
    int status=0;
    int irec, nrecs;
    jroot = json_new_object();
    recinfo = json_new_array();

    /* Open record_set */
    recordset = drms_open_records (drms_env, in, &status);
    if (!recordset) 
      JSONDIE(" jsoc_info: series not found.");
  
    nrecs = recordset->n;
    if (nrecs == 0)
      {
      json_insert_pair_into_object(jroot, "count", json_new_number("0"));
      json_insert_pair_into_object(jroot, "status", json_new_number("0"));
      json_tree_to_string(jroot,&json);
      printf("Content-type: application/json\n\n");
      printf("%s\n",json);
      fflush(stdout);
      return(0);
      }
  
    /* get list of keywords to print for each record */
    nkeys = 0;
    if (keys_listed) 
      { /* get specified list */
      char *thiskey;
      CGI_unescape_url(keylist);
      for (thiskey=strtok(keylist, ","); thiskey; thiskey=strtok(NULL,","))
	{
	if (strcmp(thiskey,"**NONE**")==0)
	  {
	  nkeys = 0;
	  break;
	  }
	if (strcmp(thiskey, "**ALL**")==0)
          {
          DRMS_Keyword_t *key;
          HIterator_t hit;
          hiter_new (&hit, &recordset->records[0]->keywords);
          while ((key = (DRMS_Keyword_t *)hiter_getnext (&hit)))
            keys[nkeys++] = strdup (key->info->name);
	  }
  	else
	  keys[nkeys++] = strdup(thiskey);
	}
      }
    free (keylist);
    /* place to put an array of keyvals per keyword */
    if (nkeys)
      keyvals = (json_t **)malloc(nkeys * sizeof(json_t *));
    for (ikey=0; ikey<nkeys; ikey++)
      {
      json_t *val = json_new_array();
      keyvals[ikey] = val;
      }
  
    /* get list of segments to show for each record */
    nsegs = 0;
    if (segs_listed) 
      { /* get specified segment list */
      char *thisseg;
      CGI_unescape_url(seglist);
      for (thisseg=strtok(seglist, ","); thisseg; thisseg=strtok(NULL,","))
	{
	if (strcmp(thisseg,"**NONE**")==0)
	  {
	  nsegs = 0;
	  break;
	  }
	if (strcmp(thisseg, "**ALL**")==0)
	  {
          DRMS_Segment_t *seg;
          HIterator_t hit;
          hiter_new (&hit, &recordset->records[0]->segments);
          while ((seg = (DRMS_Segment_t *)hiter_getnext (&hit)))
            segs[nsegs++] = strdup (seg->info->name);
	  }
  	else
	  segs[nsegs++] = strdup(thisseg);
	}
      }
    free (seglist);
    /* place to put an array of segvals per segment */
    if (nsegs)
      segvals = (json_t **)malloc(nsegs * sizeof(json_t *));
    for (iseg=0; iseg<nsegs; iseg++)
      {
      json_t *val = json_new_array();
      segvals[iseg] = val;
      }

    /* loop over set of selected records */
    for (irec = 0; irec < nrecs; irec++) 
      {
      char recquery[DRMS_MAXQUERYLEN];
      char *jsonquery;
      json_t *recobj = json_new_object();
      rec = recordset->records[irec];  /* pointer to current record */
      if (wantRecInfo)
	{
        drms_sprint_rec_query(recquery,rec);
        jsonquery = string_to_json(recquery);
        json_insert_pair_into_object(recobj, "name", json_new_string(jsonquery));
        free(jsonquery);
	}
      /* now get keyword information */
      for (ikey=0; ikey<nkeys; ikey++) 
        {
        DRMS_Keyword_t *rec_key_ikey; 
        json_t *thiskeyval = keyvals[ikey]; 
        json_t *val;
        char rawval[20000];
        char *jsonval;

        if (strcmp(keys[ikey],"*recnum*") == 0)
	  {
	  sprintf(rawval,"%ld",rec->recnum);
	  val = json_new_number(rawval);
	  }
        else if (strcmp(keys[ikey], "*logdir*") == 0)
          {
	  char *logdir = drms_record_getlogdir(rec);
	  if (logdir)
	    {
	    jsonval = string_to_json(logdir);
	    free(logdir);
	    }
	  else
	    jsonval = string_to_json("NO LOG");
	  val = json_new_string(jsonval);
  	  free(jsonval);
          }
        else
	  {
          rec_key_ikey = drms_keyword_lookup (rec, keys[ikey], 1); 
	  if (rec_key_ikey->info->type == DRMS_TYPE_TIME)
	  {
	     drms_sprintfval_format(rawval, rec_key_ikey->info->type, &(rec_key_ikey->value),
				    rec_key_ikey->info->unit, 0);
	  }
	  else
	  {
	     drms_sprintfval_format(rawval, rec_key_ikey->info->type, &(rec_key_ikey->value),
				 rec_key_ikey->info->format, 0);
	  }
          switch (rec_key_ikey->info->type)
            {
            case DRMS_TYPE_STRING:
            case DRMS_TYPE_TIME:
		jsonval = string_to_json(rawval);
		val = json_new_string(jsonval);
		free(jsonval);
		break;
            default:
		if (strcmp(rawval,"nan")==0)
		    val = json_new_number("99999999.99999");
		else
		    val = json_new_number(rawval);
            }
	  }
        json_insert_child(thiskeyval, val);
        }
  
      /* now show desired segments */
      int online = 0;
      for (iseg=0; iseg<nsegs; iseg++) 
        {
        DRMS_Segment_t *rec_seg_iseg = drms_segment_lookup (rec, segs[iseg]); 
        char *jsonpath;
        char path[DRMS_MAXPATHLEN];
        json_t *thissegval = segvals[iseg]; 
        json_t *obj = json_new_object();
        json_t *val;

        drms_record_directory (rec, path, 0);
  	strncat(path, "/", DRMS_MAXPATHLEN);
  	strncat(path, rec_seg_iseg->filename, DRMS_MAXPATHLEN);
        jsonpath = string_to_json(path);
        json_insert_child(thissegval, json_new_string(jsonpath));
        free(jsonpath);
        online = strncmp(path, "/SUM",4) == 0;
        }

      /* finish record info for this record */
      if (wantRecInfo)
	{
        json_insert_pair_into_object(recobj, "online", json_new_number(online ? "1" : "0"));
        json_insert_child(recinfo, recobj);
	}
      }
  
  /* Finished.  Clean up and exit. */
      if (wantRecInfo)
        json_insert_pair_into_object(jroot, "recinfo", recinfo);

      for (ikey=0; ikey<nkeys; ikey++) 
        {
        json_t *keyname = json_new_string(keys[ikey]); 
        json_t *keyobj = json_new_object();
        json_insert_pair_into_object(keyobj, "name", keyname);
        json_insert_pair_into_object(keyobj, "values", keyvals[ikey]);
        json_insert_child(json_keywords, keyobj);
        }
      json_insert_pair_into_object(jroot, "keywords", json_keywords);

      for (iseg=0; iseg<nsegs; iseg++) 
        {
        json_t *segname = json_new_string(segs[iseg]); 
        json_t *segobj = json_new_object();
        json_insert_pair_into_object(segobj, "name", segname);
        json_insert_pair_into_object(segobj, "values", segvals[iseg]);
        json_insert_child(json_segments, segobj);
        }
      json_insert_pair_into_object(jroot, "segments", json_segments);
      sprintf(count, "%d", nrecs);
      json_insert_pair_into_object(jroot, "count", json_new_number(count));
      json_insert_pair_into_object(jroot, "status", json_new_number("0"));
    
    drms_close_records(recordset, DRMS_FREE_RECORD);
    json_tree_to_string(jroot,&json);
    printf("Content-type: application/json\n\n");
    printf("%s\n",json);
    fflush(stdout);
    return(0);
    }
}

