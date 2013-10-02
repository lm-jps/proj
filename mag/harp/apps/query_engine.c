/*
 * query_engine: send information about a recordset to a socket
 *
 * This is a socket-compatible version of show_info and jsoc_info.  It listens
 * for metadata requests on a socket, and outputs query results in the same
 * way that jsoc_info would, as strings with JSON objects.
 *
 * Each query is a single newline-terminated line consisting of space-separated 
 * name=value pairs, where names and values are below:
 *
 * op   -- operation, one of: series_struct, rs_summary, rs_list, or exit
 * ds   -- record_set query list
 * key  -- "comma-delimited keyword list. Special values: 
 *         **ALL**, **NONE**, *recnum*, *sunum*, *size*, *online*, 
 *         *retain*, *archive*, *logdir*, *dir_mtime*
 * link -- comma delimited link list links or special values: **ALL**, **NONE**
 * seg  -- comma delimited segment list segment names or special values: 
 *         **ALL**, **NONE**
 * n    -- ceiling on the number of records to return
 *
 * "op" is mandatory, and "ds" is mandatory unless op=exit.  key, link, and seg
 * are meaningful only for op=rs_list, and are optional even for rs_list.  "n" 
 * is optional.
 *
 * This program was written to improve efficiency of metadata queries by 
 * external programs.  Obtaining metadata by capturing stdout from a shell program,
 * or by using the jsoc_info interface over HTTP, was taking too long due to 
 * program startup time.  Using this interface reduces typical single-record 
 * response times from 400ms (shell) or 200ms (HTTP) to 20ms.  For query-intensive 
 * programs, this is a big improvement.
 *
 * Usage:
 *   query_engine [-mRop] logfile=... portfile=... port=...  n=... verb=...
 *
 * where:
 *   -m:  include MIME header "Content-type: application/json" (default is not)
 *   -R:  show record query
 *   -o:  add owner info to series_struct
 *   -p:  preserve logfile (below) after exit (default is not)
 *  logfile (string): log file (default is no log)
 *        One line per request is recorded, plus notes of errors and total volume.
 *  portfile (string): port filename (default is no file)
 *        The name of the port chosen (a free, high-numbered port) is written
 *        to this file for use by clients.
 *  port (int): server port (default is to auto-assign)
 *        This is the port number.  It is auto-assigned by default.
 *  verb (int): verbosity (0, 1, 2)
 *        This writes some info to stdout.  Default is 0 (totally quiet).
 *  n (int): recordSet limit
 *        default value: 0 (no limit).  This is a session-wide limit.  It is 
 *        over-ridden by the "n=" part of the query line, for that query alone.
 *
 * Example:
 *   (in one window)
 *   % query_engine logfile=/tmp/query.log portfile=/tmp/pf
 *   (in another window)
 *   % echo 'op=rs_list ds=hmi.M_720s[$]' key=T_REC,CROTA2 | netcat localhost `cat /tmp/pf`
 *   % echo 'op=exit' | netcat localhost `cat /tmp/pf`
 *
 * Remarks:
 * It is OK to terminate query_engine with ctrl-c, although the recommended way 
 * is by sending "op=exit\n".  Erroneous queries are just skipped with the 
 * appropriate JSON reply (nonzero "status", text "message" explanation).
 * The server does not fork, so it is meant to serve only one client.  Multiple
 * server processes can run at the same time, because they are not tied to one
 * particular port.  Unless verb is nonzero (zero default), there is no output
 * to stdout.  Erroneous requests will produce log messages (if logfile is given),
 * as well as output to stderr.
 * 
 * Michael Turmon, JPL, September 2013.
 */


#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"
#include "json.h"
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>


// query length (from the socket)
#define QUERY_MAX 1024
// time string length (YYYY.MM.DD_HH:MM:SS.mmm)
#define TIMESTR_MAX 32
// prefix optionally added as header to returned value
#define JSON_MIME "Content-type: application/json\n\n"
// undefined string
#define ARG_UNDEF "@"


/*
 * singletons to encapsulate the overall server and query state
 */

// query
typedef struct {
  // times
  double t0;
  double t1;
  char t0_str[TIMESTR_MAX];
  char t1_str[TIMESTR_MAX];
  // everything else
  char line[QUERY_MAX];
  char op  [QUERY_MAX];
  char ds  [QUERY_MAX];
  char emsg[QUERY_MAX];  // error message, if nonzero status
  char remote_host[MAXHOSTNAMELEN];
  int max_recs;
  int status;  // nonzero for error
} query_bag_t;

// server
typedef struct {
  char hostname[MAXHOSTNAMELEN];
  int port;         // port (0 if not open)
  char *port_fn;    // filename to put port number in (or NULL)
  FILE *log;        // server log file
  char *log_fn;     // log filename
  int preserve_log; // preserve log file after exit?
  unsigned int 
    export_num;     // number of export requests
  unsigned long 
    export_len;     // length of response
  double 
    timeout_time;   // seconds until server times out
  int timeout;      // has the server timed out?
  int running;      // 0 if we are shutting it down
  int sock;         // socket
  int connfd;       // connection fd
  int head;         // attach mime head?
} server_bag_t;


// set according to verbosity input parameter
static int verbflag;
// log file pointers (initialized to NULL)
static FILE *LOGout = NULL;
static FILE *LOGerr = NULL;

/************************************************************* 
 *
 * ERROR HANDLING/LOGGING MACROS
 *
 *************************************************************
 */

// V_printf: facilitate verbflag output
//   if flag is > 0, output is to stdout+LOGout, if < 0, to stderr+LOGerr
//   if flag is 0, no output is made at all
// If LOGout/LOGerr are NULL, they are skipped.
//   It is legal to have LOGerr = stdout, so that error messages are 
//   duplicated to stdout.
// The message is printed in the form:
//   <message>                               (if preamble is NULL)
//   <module_name>: <preamble><message>      (if preamble is non-NULL)
//   <module_name>: <preamble+1><message>\n  (if non-NULL && *preamble==\n)
// Usage:
//   V_printf(verbflag > 0, "\t", "Mask(%d) = %d\n", 2048, mask[2048]);
void
V_printf(int flag, char *preamble, char *format, ...) {
  va_list args;
  extern char *module_name;
  extern FILE *LOGout;  // NULL OK to suppress
  extern FILE *LOGerr;  // NULL OK to suppress
  FILE *fp = (flag > 0) ? stdout : stderr;
  FILE *lp = (flag > 0) ? LOGout : LOGerr;
  int append_newline;

  if (flag != 0) {
    // print to stream from stdio
    va_start(args, format);
    append_newline = 0;
    if (preamble) {
      if (preamble[0] == '\n') {
	append_newline = 1;
	preamble++;
      }
      fprintf(fp, "%s: %s", module_name, preamble);
    }
    vfprintf(fp, format, args); // va_list version of printf
    if (append_newline)
      fprintf(fp, "\n");
    fflush(fp);
    va_end(args);
    // optionally print to log file
    if (lp) { 
      va_start(args, format);
      if (preamble) 
	fprintf(lp, "%s: %s", module_name, preamble);
      vfprintf(lp, format, args);
      if (append_newline)
	fprintf(fp, "\n");  // append_newline was set above
      fflush(lp);
      va_end(args);
    }
  } else {
    va_start(args, format);
    // (do nothing)
    va_end(args);
  }
}
// exit with error
//  (standard trick to swallow the semicolon)
//  (use abbreviation where putting \n in V_printf preamble writes \n 
//   at tail of V_printf line; this trick allows us to use V_printf with
//   custom printf-style DIE messages like:
//   DIE("Could not read input file `%s'.", input_file);
// NOTE: return() comes back from DoIt, so only call within DoIt()
#define DIE(...) do { \
	fflush(stdout); \
        V_printf(-1, "\nFATAL.  ", __VA_ARGS__); \
	if (LOGout) fclose(LOGout); \
	if (LOGerr) fclose(LOGerr); \
        return 1; \
        } while (0)
// report non-fatal error
//  (standard trick to swallow the semicolon)
#define WARN(...) do { \
	fflush(stdout); \
        V_printf(-1, "\nWARNING (continuing).  ", __VA_ARGS__); \
        } while (0)

/*************************************************************
 *
 * DRMS UTILITIES
 *
 *************************************************************/

/* drms_record_getlogdir */

/* Returns path of directory that contains any saved log information for the given record */
/* If log is offline, returns message, if log was not saved or otherwise not found returns NULL */
/* The returned char* should be freed after use. */

char *drms_record_getlogdir(DRMS_Record_t *rec)
  { 
  char *logpath;
  char query[DRMS_MAXQUERYLEN];
  DB_Text_Result_t *qres;

  sprintf(query, "select sunum from %s.drms_session where sessionid=%lld", rec->sessionns, rec->sessionid);
  if ((qres = drms_query_txt(drms_env->session, query)) && qres->num_rows>0)
    {
    if (qres->field[0][0][0] == '\0')
      logpath = strdup("No log avaliable");
    else
      {
      DRMS_StorageUnit_t *su;
      int status, save_retention = drms_env->retention;
      int retrieve = 0;
      su = malloc(sizeof(DRMS_StorageUnit_t));
      su->sunum = atoll(qres->field[0][0]);
      drms_env->retention = DRMS_LOG_RETENTION;
      // FIXME: 
      // In server mode (i.e., "make query_engine_sock"), the build
      // fails.  arta 30 sep 2013 suggests:
      //   "For drms_su_getsudir() use drms_getunit() (which checks 
      //    the DRMS_CLIENT macro internally)."
      // But this replacement has a different signature.  Punt for now,
      // it is not used anyway.
#ifdef DRMS_CLIENT
      status = 1; // Copy what was done for the server.
      char *or_else = "Log offline";
#else
      status = 1; // act as if it failed
      char *or_else = "Logdir unimplemented";
#endif
      if (!status)
        logpath = strdup(su->sudir);
      else
        logpath = strdup(or_else);
      free(su);
      drms_env->retention = save_retention;
      }
    }
  else
    logpath = strdup("Log query failed");
  db_free_text_result(qres);
  return logpath;
  }

int drms_ismissing_keyval(DRMS_Keyword_t *key)
  {
  XASSERT(key);
  switch(key->info->type)
    {
    case DRMS_TYPE_CHAR:
      return(drms_ismissing_char(key->value.char_val));
    case DRMS_TYPE_SHORT:
      return(drms_ismissing_short(key->value.short_val));
    case DRMS_TYPE_INT:
      return(drms_ismissing_int(key->value.int_val));
    case DRMS_TYPE_LONGLONG:
      return(drms_ismissing_longlong(key->value.longlong_val));
    case DRMS_TYPE_FLOAT:
      return(drms_ismissing_float(key->value.float_val));
    case DRMS_TYPE_DOUBLE:
      return(drms_ismissing_double(key->value.double_val));
    case DRMS_TYPE_TIME:
      return(drms_ismissing_time(key->value.time_val));
    case DRMS_TYPE_STRING:
      return(drms_ismissing_string(key->value.string_val));
    default:
      V_printf(-1, "", "ERROR: Unhandled DRMS type %d\n",(int)key->info->type);
      XASSERT(0);
    }
  return 0;
  }


/* find first record in series that owns the given record */
DRMS_RecordSet_t *drms_find_rec_first(DRMS_Record_t *rec, int wantprime)
  {
  int nprime;
  int status;
  DRMS_RecordSet_t *rs;
  char query[DRMS_MAXQUERYLEN];
  strcpy(query, rec->seriesinfo->seriesname);
  nprime = rec->seriesinfo->pidx_num;
  if (wantprime && nprime > 0) 
    // only first prime key is used for now
     // for (iprime = 0; iprime < nprime; iprime++)
      strcat(query, "[^]");
  else
    strcat(query, "[:#^]");
  V_printf(verbflag > 1, "", "test 1 query is %s\n",query);
  rs = drms_open_nrecords(rec->env, query, 1, &status);
  V_printf(verbflag > 1, "", "test 1 status is %d\n",status);
  return(rs);
  }

/* find last record in series that owns the given record */
DRMS_RecordSet_t *drms_find_rec_last(DRMS_Record_t *rec, int wantprime)
  {
  int nprime;
  int status;
  DRMS_RecordSet_t *rs;
  char query[DRMS_MAXQUERYLEN];
  strcpy(query, rec->seriesinfo->seriesname);
  nprime = rec->seriesinfo->pidx_num;
  if (wantprime && nprime > 0) 
    // only first prime key is used for now
     // for (iprime = 0; iprime < nprime; iprime++)
      strcat(query, "[$]");
  else
    strcat(query, "[:#$]");
  rs = drms_open_nrecords(rec->env, query, -1, &status);
  return(rs);
  }

/* returns series owner as static string */
char *drms_getseriesowner(DRMS_Env_t *drms_env, char *series, int *status)
   {
   char *nspace = NULL;
   char *relname = NULL;
   int istat = DRMS_SUCCESS;
   DB_Text_Result_t *qres = NULL;
   static char owner[256];
   owner[0] = '\0';

   if (!get_namespace(series, &nspace, &relname))
      {
      char query[1024];
      strtolower(nspace);
      strtolower(relname);

      snprintf(query, sizeof(query), "SELECT pg_catalog.pg_get_userbyid(T1.relowner) AS owner FROM pg_catalog.pg_class AS T1, (SELECT oid FROM pg_catalog.pg_namespace WHERE nspname = '%s') AS T2 WHERE T1.relnamespace = T2.oid AND T1.relname = '%s'", nspace, relname);

      if ((qres = drms_query_txt(drms_env->session, query)) != NULL)
         {
         if (qres->num_cols == 1 && qres->num_rows == 1)
            strcpy(owner, qres->field[0][0]);
         db_free_text_result(qres);
         }
      }
   *status = owner[0] != 0;
   return(owner);
   }

static char * 
string_to_json(char *in)
{ 
  // for json vers 0.9 no longer uses wide chars
  char *new;
  new = json_escape(in);
  return new;
}


static void list_series_info(DRMS_Env_t *drms_env, DRMS_Record_t *rec, json_t *jroot, int wantOwner)
{
  DRMS_Keyword_t *key;
  DRMS_Segment_t *seg;
  DRMS_Link_t *link;
  HIterator_t *last = NULL;
  char intstring[100];
  char *notework;
  char *owner;
  json_t *indexarray, *primearray, *keyarray, *segarray, *linkarray;
  json_t *primeinfoarray;
  int npkeys;
  int status;
  char prevKeyName[DRMS_MAXNAMELEN] = "";
  char baseKeyName[DRMS_MAXNAMELEN];

  /* add description from seriesinfo */
  notework = string_to_json(rec->seriesinfo->description);
  json_insert_pair_into_object(jroot, "note", json_new_string(notework));
  free(notework);
  /* add retention, unitsize, archive, and tapegroup integers */
  sprintf(intstring, "%d", rec->seriesinfo->retention);
  json_insert_pair_into_object(jroot, "retention", json_new_number(intstring));
  sprintf(intstring, "%d", rec->seriesinfo->unitsize);
  json_insert_pair_into_object(jroot, "unitsize", json_new_number(intstring));
  sprintf(intstring, "%d", rec->seriesinfo->archive);
  json_insert_pair_into_object(jroot, "archive", json_new_number(intstring));
  sprintf(intstring, "%d", rec->seriesinfo->tapegroup);
  json_insert_pair_into_object(jroot, "tapegroup", json_new_number(intstring));
  /* add owner for series */
  if (wantOwner) {
    owner = string_to_json(drms_getseriesowner(drms_env, rec->seriesinfo->seriesname, &status));
    json_insert_pair_into_object(jroot, "owner", json_new_string(owner));
    free(owner);
  }
  
  /* show the prime index keywords */
  // both the original simple list and new array of objects are generated -- XXXXX REMOVE SOMEDAY                             
  // for compatibility.  The older "primekeys" array may be removed in the future.  23Nov09 -- XXXXX REMOVE SOMEDAY           
  // old lines marked with trailing XXXXX REMOVE SOMEDAY
  primearray = json_new_array(); // XXXXX REMOVE SOMEDAY
  primeinfoarray = json_new_array();
  npkeys = rec->seriesinfo->pidx_num;
  if (npkeys > 0)
    {
    int i;
    json_t *primeinfo, *val;
    char *jsonstr;
    for (i=0; i<npkeys; i++)
        {
        json_t *primeinfo = json_new_object();
        DRMS_Keyword_t *pkey, *skey;
        pkey = rec->seriesinfo->pidx_keywords[i];
	if (pkey->info->recscope > 1)
        {
           char rawval[100];
           skey = drms_keyword_slotfromindex(pkey);
           json_insert_pair_into_object(primeinfo, "name", json_new_string(skey->info->name));
           json_insert_pair_into_object(primeinfo, "slotted", json_new_number("1"));
           drms_keyword_snprintfval(drms_keyword_stepfromvalkey(skey), rawval, sizeof(rawval));
           jsonstr = string_to_json(rawval);
           val = json_new_string(jsonstr);
           free(jsonstr);
           json_insert_pair_into_object(primeinfo, "step", val);
           json_insert_child(primearray, json_new_string(skey->info->name)); // XXXXX REMOVE SOMEDAY                         

        }
        else
        {
           json_insert_pair_into_object(primeinfo, "name", json_new_string(pkey->info->name));
           json_insert_pair_into_object(primeinfo, "slotted", json_new_number("0"));
           json_insert_child(primearray, json_new_string(pkey->info->name)); // XXXXX REMOVE SOMEDAY                         
        }
        json_insert_child(primeinfoarray, primeinfo);            
        }
    }
  else
  {
     json_insert_child(primearray, json_new_null()); // XXXXX REMOVE SOMEDAY                                                   
     json_insert_child(primeinfoarray, json_new_null());
  }
  json_insert_pair_into_object(jroot, "primekeys", primearray); // XXXXX REMOVE SOMEDAY                                       
  json_insert_pair_into_object(jroot, "primekeysinfo", primeinfoarray);
 
  /* show DB index keywords */
  indexarray = json_new_array();
  if (rec->seriesinfo->dbidx_num > 0)
    {
    int i;
    for (i=0; i<rec->seriesinfo->dbidx_num; i++)
        json_insert_child(indexarray, json_new_string((rec->seriesinfo->dbidx_keywords[i])->info->name));
    }
  else
    json_insert_child(indexarray, json_new_null());
  json_insert_pair_into_object(jroot, "dbindex", indexarray);

  V_printf(verbflag > 1, "", "   starting all keywords\n");
  /* show all keywords */
  keyarray = json_new_array();

  /* We don't want to follow the link just yet - we need a combination of source and target 
   * keyword information below. For now, just get the source keyword info. */
  while ((key = drms_record_nextkey(rec, &last, 0)))
    {
    json_t *keyinfo;
    json_t *keytype;
    json_t *defval, *recscope;
    char rawval[100], *jsonstr;
    int persegment = key->info->kwflags & kKeywordFlag_PerSegment;
    V_printf(verbflag > 1, "", "   starting keyword %s\n",key->info->name);
    /* If a keyword is a linked keyword, then it cannot be a per-segment keyword also. */
    if (persegment)
      {
      char *underscore;
      strcpy(baseKeyName, key->info->name);
      underscore = rindex(baseKeyName, '_');
      if (underscore) *underscore = '\0';
      if (strcmp(prevKeyName, baseKeyName) == 0)
        continue;  // only report the first instance of persegment keywords.
      strcpy(prevKeyName, baseKeyName);
      }
    keyinfo = json_new_object();
    json_insert_pair_into_object(keyinfo, "name", json_new_string(persegment ? baseKeyName : key->info->name));
    if (key->info->islink)
    {
       /* provide link name and target keyword name */
       char linknames[100], *tmpstr;
       json_t *linkinfo;
       sprintf(linknames,"%s->%s", key->info->linkname, key->info->target_key);
       tmpstr = string_to_json(linknames);
       linkinfo = json_new_string(tmpstr);
       json_insert_pair_into_object(keyinfo, "linkinfo", linkinfo);
       free(tmpstr);

       /* Display the target keyword data type. Must follow link now. */
       int lnkstat = DRMS_SUCCESS;
       DRMS_Keyword_t *linkedkw = drms_template_keyword_followlink(key, &lnkstat);

       if (lnkstat == DRMS_SUCCESS && linkedkw)
       {
          keytype = json_new_string(drms_type_names[linkedkw->info->type]);
       }
       else
       {
          keytype = json_new_string("link");
       }
    }
    else
       keytype = json_new_string(drms_type_names[key->info->type]);
    json_insert_pair_into_object(keyinfo, "type", keytype);
    // scope                                                                                                                      
    // redundant - persegment = key->info->kwflags & kKeywordFlag_PerSegment;
    if (persegment)
      {
      recscope = json_new_string("segment");
      }
    else
      {
      jsonstr = string_to_json((char *)drms_keyword_getrecscopestr(key, NULL));
      recscope = json_new_string(jsonstr);
      free(jsonstr);
      }
    json_insert_pair_into_object(keyinfo, "recscope", recscope);
    // default                                                                                                                    
    drms_keyword_snprintfval(key, rawval, sizeof(rawval));
    jsonstr = string_to_json(rawval);
    defval = json_new_string(jsonstr);
    free(jsonstr);
    json_insert_pair_into_object(keyinfo, "defval", defval);

    notework = string_to_json(key->info->unit);
    json_insert_pair_into_object(keyinfo, "units", json_new_string(notework));
    free(notework);
    notework = string_to_json(key->info->description);
    json_insert_pair_into_object(keyinfo, "note", json_new_string(notework));
    free(notework);
    json_insert_child(keyarray, keyinfo);
    }
  json_insert_pair_into_object(jroot, "keywords", keyarray);

  V_printf(verbflag > 1, "", " done with keywords, start segments\n");
  
  /* show the segments */
  segarray = json_new_array();
  if (rec->segments.num_total)
    {
       if (last)
       {
          hiter_destroy(&last);
       }

      while ((seg = drms_record_nextseg(rec, &last, 0)))
      { /* segment name, units, protocol, dims, description */
      json_t *seginfo = json_new_object();
      int naxis = seg->info->naxis;
      json_insert_pair_into_object(seginfo, "name", json_new_string(seg->info->name));
      if (seg->info->islink)
	    {
            char linkinfo[DRMS_MAXNAMELEN+10];
	    sprintf(linkinfo, "link via %s", seg->info->linkname);
            json_insert_pair_into_object(seginfo, "type", json_new_string(drms_type_names[seg->info->type]));
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
            json_insert_pair_into_object(seginfo, "type", json_new_string(drms_type_names[seg->info->type]));
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
  json_insert_pair_into_object(jroot, "segments", segarray);

  V_printf(verbflag > 1, "", " done with segments, start links\n");
  /* show the links */
  linkarray = json_new_array();
  if (rec->links.num_total)
    {
      if (last)
	{
          hiter_destroy(&last);
	}

      while ((link = drms_record_nextlink(rec, &last)))
	{
	  V_printf(verbflag > 1, "", " link: %s\n",link->info->name);
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
  json_insert_pair_into_object(jroot, "links", linkarray);

  if (last)
    {
      hiter_destroy(&last);
    }
  return;
}


static int get_series_stats(DRMS_Record_t *rec, json_t *jroot)
  {
  DRMS_RecordSet_t *rs;
  int nprime;
  int status;
  char query[DRMS_MAXQUERYLEN];
  json_t *interval = json_new_object();

  nprime = rec->seriesinfo->pidx_num;
  if (nprime > 0)
    sprintf(query,"%s[^]", rec->seriesinfo->seriesname);
  else
    sprintf(query,"%s[:#^]", rec->seriesinfo->seriesname);
  rs = drms_open_nrecords(rec->env, query, 1, &status);

  if (status == DRMS_ERROR_QUERYFAILED)
  {
     if (rs) 
     {
        drms_free_records(rs);
     }

     return status;
  }

  if (!rs || rs->n < 1)
    {
    json_insert_pair_into_object(interval, "FirstRecord", json_new_string("NA"));
    json_insert_pair_into_object(interval, "FirstRecnum", json_new_string("NA"));
    json_insert_pair_into_object(interval, "LastRecord", json_new_string("NA"));
    json_insert_pair_into_object(interval, "LastRecnum", json_new_string("NA"));
    json_insert_pair_into_object(interval, "MaxRecnum", json_new_number("0"));
    if (rs) drms_free_records(rs);
    json_insert_pair_into_object(jroot, "Interval", interval);
    return DRMS_SUCCESS;
    }
  else
    {
    char recquery[DRMS_MAXQUERYLEN];
    char *jsonquery;
    char val[100];
    int status;
    drms_sprint_rec_query(recquery,rs->records[0]);
    jsonquery = string_to_json(recquery);
    status = json_insert_pair_into_object(interval, "FirstRecord", json_new_string(jsonquery));
    if (status != JSON_OK) 
      V_printf(-1, "", "json_insert_pair_into_object, status=%d, text=%s\n",status,jsonquery);
    free(jsonquery);
    sprintf(val,"%lld", rs->records[0]->recnum);
    json_insert_pair_into_object(interval, "FirstRecnum", json_new_number(val));
    drms_free_records(rs);
  
    if (nprime > 0)
      sprintf(query,"%s[$]", rec->seriesinfo->seriesname);
    else
      sprintf(query,"%s[:#$]", rec->seriesinfo->seriesname);
    rs = drms_open_nrecords(rec->env, query, -1, &status);

    if (status == DRMS_ERROR_QUERYFAILED)
    {
       if (rs)
       {
          drms_free_records(rs);
       }

       return status;
    }

    drms_sprint_rec_query(recquery,rs->records[0]);
    jsonquery = string_to_json(recquery);
    json_insert_pair_into_object(interval, "LastRecord", json_new_string(jsonquery));
    free(jsonquery);
    sprintf(val,"%lld", rs->records[0]->recnum);
    json_insert_pair_into_object(interval, "LastRecnum", json_new_number(val));
    drms_free_records(rs);
 
    sprintf(query,"%s[:#$]", rec->seriesinfo->seriesname);
    rs = drms_open_records(rec->env, query, &status);

    if (status == DRMS_ERROR_QUERYFAILED)
    {
       if (rs)
       {
          drms_free_records(rs);
       }

       return status;
    }

    sprintf(val,"%lld", rs->records[0]->recnum);
    json_insert_pair_into_object(interval, "MaxRecnum", json_new_number(val));
    drms_free_records(rs);
    }
  json_insert_pair_into_object(jroot, "Interval", interval);
  return 0;
}


/*************************************************************
 *
 * SERVER UTILITIES
 *
 *************************************************************/

static
void
current_time(double *t_num, char *t_str)
{
  struct timeval tv;
  time_t nowtime;
  struct tm *nowtm;
  char tmbuf[TIMESTR_MAX];

  gettimeofday(&tv, NULL);
  if (t_num)
    *t_num = tv.tv_sec + tv.tv_usec/1000000.0;
  if (t_str) {
    nowtime = tv.tv_sec;
    nowtm = localtime(&nowtime);
    strftime(tmbuf, TIMESTR_MAX, "%Y.%m.%d_%H:%M:%S", nowtm);
    // round to millisec
    snprintf(t_str, TIMESTR_MAX, "%s.%03d", tmbuf, (int) (tv.tv_usec/1000));
  }
}

static
void
insert_log_header(server_bag_t *server_bag)
{
  extern char *module_name;
  FILE *fp = server_bag->log;
  char tnow_str[TIMESTR_MAX];

  if (!fp) return;
  current_time(NULL, tnow_str);
  fprintf(fp, "# %s: transaction log\n", module_name);
  fprintf(fp, "# running as PID %d on %s:%d\n", 
	  getpid(), server_bag->hostname, server_bag->port);
  fprintf(fp, "# start time = %s\n", tnow_str);
  fflush(fp);
}

/*
 * puts in the footer and closes the file
 */
static
void
insert_log_footer(server_bag_t *server_bag)
{
  extern char *module_name;
  char tnow_str[TIMESTR_MAX];
  FILE *fp = server_bag->log;

  if (!fp) return;
  current_time(NULL, tnow_str);
  fprintf(fp, "# %s: closing transaction log\n", module_name);
  fprintf(fp, "# end time = %s\n", tnow_str);
  fprintf(fp, "# requests = %d\n", server_bag->export_num);
  fprintf(fp, "# bytes = %ld\n",   server_bag->export_len);
  fflush(fp);
  fclose(fp);
  server_bag->log = NULL;
}

// insert_log_summary - record this call of the server.
static
void
insert_log_summary(server_bag_t *Sbag,
		   query_bag_t *Qbag)
{
  FILE *log = Sbag->log;

  if (!log) return;
  if (Qbag->status == 0) {
    // success
    fprintf(log, "host='%s'\t", Sbag->hostname);
    fprintf(log, "lag=%0.3f\t", (Qbag->t1) - (Qbag->t0));
    fprintf(log, "IP='%s'\t",   Qbag->remote_host);
    fprintf(log, "op='%s'\t",   Qbag->op);
    fprintf(log, "ds='%s'\t",   Qbag->ds);
    fprintf(log, "n=%d\t",      Qbag->max_recs);
    fprintf(log, "status=%d\n", Qbag->status);
  } else {
    // always have t1_str
    fprintf(log, "# error on following line at %s:\n", Qbag->t1_str);
    fprintf(log, "# |%s|\n", Qbag->line);
    fprintf(log, "# Error was: %s\n", *(Qbag->emsg) ? Qbag->emsg : "Unknown");
  }
  fflush(log);
}

// insert_log_message -- put a message into the log
static
void
insert_log_message(server_bag_t *Sbag, char *msg)
{
  FILE *log = Sbag->log;
  if (!log) return;
  fprintf(log, "# %s\n", msg);
  fflush(log);
}

/*
 * set up server socket and log file
 */

static
char *
server_setup(int port_given, 
	     const char *port_fn, 
	     const char *log_fn, 
	     server_bag_t *server_bag)
{
  int sock = 0;
  struct sockaddr_in serv_addr; 
  socklen_t serv_addr_len = sizeof(serv_addr);
  int rc;
  int port;
  FILE *fp;

  // fill in local hostname (zero it because gethostname may not write last \0)
  memset(server_bag->hostname, 0, sizeof(server_bag->hostname));
  gethostname(server_bag->hostname, sizeof(server_bag->hostname)-1);

  sock = socket(AF_INET, SOCK_STREAM, 0);
  if (sock < 0)
    return "socket() failed";

  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  if (port_given > 0)
    serv_addr.sin_port = htons((short) port_given);
  else
    serv_addr.sin_port = htons(0);
  rc = bind(sock, (struct sockaddr*) &serv_addr, sizeof(serv_addr)); 
  if (rc < 0)
    return "bind() failed";

  rc = getsockname(sock, (struct sockaddr*) &serv_addr, &serv_addr_len); 
  if (rc < 0)
    return "getsockname() failed";
  port = ntohs(serv_addr.sin_port);

  server_bag->port = port;
  server_bag->sock = sock;

  // set up the log file, if desired
  if (strcmp(log_fn, ARG_UNDEF) == 0) {
    server_bag->log = NULL;
    server_bag->log_fn = NULL;
  } else {
    if ((fp = fopen(log_fn, "w")) != NULL) {
      server_bag->log = fp;
    } else {
      server_bag->log = NULL;
      return "Could not create log file";
    }
    server_bag->log_fn = strdup(log_fn);
  }
  insert_log_header(server_bag); // NULL log is OK

  // write the port chosen to a file, if desired
  if (strcmp(port_fn, ARG_UNDEF) == 0) {
    server_bag->port_fn = NULL;
  } else {
    if ((fp = fopen(port_fn, "w")) == NULL) {
      V_printf(-1, "", "Could not create port file `%s'\n", port_fn);
      return "Could not create port file";
    }
    fprintf(fp, "%d\n", (int) port);
    fclose(fp);
    server_bag->port_fn = strdup(port_fn);
  }
  return NULL; // OK
}

/*
 * tear down server log file and socket
 *   note: the OnSIGINT handler also traps thru this routine
 */
static
char *
server_teardown(server_bag_t *server_bag)
{
  int rc;
  char *msg = NULL;

  // remove port filename, if it was created
  if (server_bag->port_fn != NULL)
    unlink(server_bag->port_fn);  // fail silently
  // clean termination in log (closes log)
  insert_log_footer(server_bag); // NULL log is OK
  // close socket
  if (server_bag->port > 0) {
    server_bag->port = 0;
    rc = close(server_bag->sock);
    if (rc < 0)
      msg = "Error on socket close";
  }
  // remove log, unless desired
  if (server_bag->log_fn) {
    if (!server_bag->preserve_log)
      unlink(server_bag->log_fn);  // fail silently
    free(server_bag->log_fn);
  }
  return msg;
}


/*
 * blocking read on fd, but time-out after t0 seconds, returns:
 *   -1: if error
 *    0: if timed out
 *    1: if fd is ready to read
 */
static
int 
isready(int fd, double t0)
{
  int rc;
  fd_set fds;
  struct timeval tv;
  
  // socket wait time
  tv.tv_sec = (int) floor(t0);
  tv.tv_usec = (int) floor((t0 - tv.tv_sec) * 1e-6);
  
  // fd set
  FD_ZERO(&fds);
  FD_SET(fd, &fds);

  rc = select(fd+1, &fds, NULL, NULL, &tv);
  if (rc < 0)
    return -1;
  return FD_ISSET(fd,&fds) ? 1 : 0;
}

static
char *
input_query(server_bag_t *server_bag, query_bag_t *query_bag)
{
  int connfd;
  int rc;
  int n;
  struct sockaddr remote_addr;
  socklen_t addr_size = sizeof remote_addr;

  // we reset this below, but we set it up here because some code paths
  // lead out of this routine first; this time is better than nothing.
  current_time(&(query_bag->t0), query_bag->t0_str);

  server_bag->timeout = 0;
  // use of select() here allows us to close server due to inactivity
  rc = isready(server_bag->sock, server_bag->timeout_time);
  if (rc == 0) {
    server_bag->timeout = 1;
    return "server timed out";
  } else if (rc < 0) {
    return "select() failed";
  }
  
  // (since we select()'ed above, this should read immediately
  //   (when Matlab signaling fails, this can return 0, which is not an error,
  //    but will eventually cause a failure)
  connfd = accept(server_bag->sock, &remote_addr, &addr_size); 
  if (connfd < 0)
    return "accept() failed";
  else if (connfd == 0) 
    V_printf(-1, "", "input_query: got connfd == 0 from sock = %d\n",
	     server_bag->sock);

  // the first read gets it all, and further reads block
  n = sizeof(query_bag->line) - 1; // room for \0
  rc = read(connfd, query_bag->line, n);
  if (rc < 0)
    return "read failed";
  // make it a string
  query_bag->line[rc] = '\0';
  // chop at first newline
  n = strcspn(query_bag->line, "\n\r");
  if ((query_bag->line)[n] == '\0')
    return "read failed: no newline seen (line too long?)";
  (query_bag->line)[n] = '\0'; 
  
  // get remote hostname
  rc = getnameinfo(&remote_addr, addr_size, 
		   query_bag->remote_host, 
		   sizeof(query_bag->remote_host), NULL, 0, 0);
  if (rc != 0)
    strcpy(query_bag->remote_host, "unresolved_host");

  V_printf(verbflag > 1, "", "remote host = %s\n", query_bag->remote_host);
  V_printf(verbflag > 1, "", "got <%s>\n", query_bag->line);

  current_time(&(query_bag->t0), query_bag->t0_str);
  // this also means a request is in progress
  server_bag->connfd = connfd;

  return NULL;
}


/*
 * decompose a line into its pieces
 * 
 * the line looks like:
 *  op=rs_list ds=hmi.M_720s key=T_REC,CROTA2 seg=magnetogram
 * although we allow white space around the equal sign, and at line's beginning and end.
 */
static
char *
parse_line(query_bag_t *qBag, 
	   int srv_max_recs,
	   char *keylist, 
	   char *seglist, 
	   char *linklist)
{
  const char *delims = " \t\n\r";
  char *token;
  char *argname, *argvalu;
  char *line;

  // defaults
  strcpy(qBag->op, "rs_list");   // assume rs_list
  strcpy(qBag->ds, "");
  strcpy(keylist,  "");
  strcpy(seglist,  "");
  strcpy(linklist, "");
  qBag->max_recs = srv_max_recs; // start with session default
  // this routine uses "line" as a workspace
  line = strdup(qBag->line);
  // parse the string
  while ((token = strsep(&line, delims)) != NULL) {
    if (*token == '\0') continue; // allow extra whitespace
    argname = token;
    argvalu = strchr(token, (int) '=');
    if (!argvalu) {
      free(line);
      return "Missing = in arg=value pair in input line";
    }
    *argvalu++ = '\0'; // null the = sign, and move to next char
    // skip leading whitespace in argvalu
    argvalu += strspn(argvalu, delims); 
    // null out trailing whitespace in argvalu (or re-null the final null)
    *(argvalu + strcspn(argvalu, delims)) = '\0';
    // null out trailing whitespace in argname (or re-null the final null)
    *(argname + strcspn(argname, delims)) = '\0';
    // case on argname
    if (strcmp(argname, "op") == 0) {
      strcpy(qBag->op, argvalu);
    } else if (strcmp(argname, "ds") == 0) {
      strcpy(qBag->ds, argvalu);
    } else if (strcmp(argname, "key") == 0) {
      strcpy(keylist, argvalu);
    } else if (strcmp(argname, "seg") == 0) {
      strcpy(seglist, argvalu);
    } else if (strcmp(argname, "link") == 0) {
      strcpy(linklist, argvalu);
    } else if (strcmp(argname, "n") == 0) {
      if (sscanf(argvalu, "%d", &(qBag->max_recs)) != 1) {
	qBag->max_recs = srv_max_recs; // reset to something sane
	free(line);
	return "Bad number in n=max_recs parameter in input line";
      }
    } else {
      free(line);
      return "Unrecognized arg name in arg=value pair in input line";
    }
  }
  free(line);
  // ds is required, unless exiting
  if (strcmp(qBag->op, "exit") != 0 && strlen(qBag->ds) == 0)
    return "Missing ds=value pair in input line";
  // OK
  return NULL;
}


/*************************************************************
 *
 * JSON EXPORT
 *
 *************************************************************/

/*
 * includes the necessary loop to send all of buf
 */
int 
write_buf(int fd, char *buf, int len)
{
  int rc;

  while ((rc = write(fd, buf, len)) > 0) {
    buf += rc;
    len -= rc;
  }
  if (rc < 0)
    return -1; // error
  else if (rc == 0 && len == 0)
    return 0; // complete write
  else
    return 1; // incomplete write
}

/*
 * export json_t object to file
 *
 * status of 0 indicates OK
 */
static
void
export_json(server_bag_t *server_bag, query_bag_t *query_bag, json_t *jroot)
{
  char *initial_json;
  char *final_json;
  char runtime[32];  // a float
  char status_string[12]; // an integer
  int rc, ok1, ok2;
  int fd = server_bag->connfd;
  int head = server_bag->head;
  int len1, len2;

  // set up end-time
  current_time(&(query_bag->t1), query_bag->t1_str);
  // put runtime in
  snprintf(runtime, sizeof(runtime), "%0.3f", query_bag->t1 - query_bag->t0);
  json_insert_pair_into_object(jroot, "runtime", json_new_number(runtime));
  // status
  snprintf(status_string, sizeof(status_string), "%8d", query_bag->status);
  json_insert_pair_into_object(jroot, "status", json_new_number(status_string));
  // convert to text 
  json_tree_to_string(jroot, &initial_json);
  final_json = json_format_string(initial_json);
  free(initial_json);
  // write it out
  ok1 = ok2 = 0; // assume trouble
  if (head) {
    len1 = strlen(JSON_MIME);
    rc = write_buf(fd, JSON_MIME, strlen(JSON_MIME));
    if (rc < 0)
      V_printf(-1, "", "Error writing json mime\n");
    else if (rc > 0)
      V_printf(-1, "", "Incomplete write on json mime\n");
    else
      ok1 = 1;
  } else {
    ok1  = 1;
    len1 = 0;
  }
  len2 = strlen(final_json);
  rc = write_buf(fd, final_json, len2);
  if (rc < 0)
    V_printf(-1, "", "Error writing json buf\n");
  else if (rc > 0)
    V_printf(-1, "", "Incomplete write on json buf\n");
  else
    ok2 = 1;

  rc = close(fd);
  if (rc < 0)
    V_printf(-1, "", "Error closing json buf\n");

  server_bag->connfd = 0; // indicate no request being processed
  free(final_json);
  // log a response error, if any
  if (!ok1 || !ok2 || rc < 0)
    insert_log_message(server_bag, 
		       "communication error responding to request below.");
  // log the request last
  insert_log_summary(server_bag, query_bag);
  server_bag->export_num += 1;
  server_bag->export_len += len1 + len2;
}


/*
 * export an exceptional status
 * 
 *   if query_bag->status is nonzero, an error is assumed
 *   if status is zero, it's a meta-operation like exiting the server
 *
 */
static
void
export_json_exception(server_bag_t *server_bag, query_bag_t *query_bag, char *msg)
{
  char *emsg;
  char *emsg_json;
  json_t *jroot = json_new_object();

  emsg = msg ? msg : "Unspecified error";
  // set up for logging later
  snprintf(query_bag->emsg, sizeof(query_bag->emsg), "%s", emsg);
  if (query_bag->status)
    V_printf(verbflag, "", "Failed request (%s):\n%s\n",
	     emsg,
	     query_bag->line);
  emsg_json = string_to_json(emsg);
  if (query_bag->status) 
    json_insert_pair_into_object(jroot, "error",   json_new_string(emsg_json));
  else
    json_insert_pair_into_object(jroot, "message", json_new_string(emsg_json));
  export_json(server_bag, query_bag, jroot);
  json_free_value(&jroot);
}


/*************************************************************
 *
 * CONTROL WITHIN MAIN ROUTINE
 *
 *************************************************************/

char *module_name = "query_engine";

ModuleArgs_t module_args[] = { 
  {ARG_STRING, "logfile", ARG_UNDEF, "log file (default is no log)"},
  {ARG_STRING, "portfile",ARG_UNDEF, "port filename (default is no file)"},
  {ARG_INT,    "port",    "0",       "server port (default is to auto-assign)"},
  {ARG_INT,    "verb",    "0",       "verbosity (0, 1, 2)"},
  {ARG_INT,    "n",       "0",       "recordSet limit"},
  {ARG_FLOAT,  "timeout", "30.0",    "server-exit timeout (minutes, >0)"},
  {ARG_FLAG,   "m",       "0",       "include MIME header"},
  {ARG_FLAG,   "p",       "0",       "preserve logfile after exit"},
  {ARG_FLAG,   "R",       "0",       "show record query"},
  {ARG_FLAG,   "o",       "0",       "add owner info to series_struct"},
  {ARG_END}
};


// skip finding response to the current request
//   (just a convenient abbreviation)
#define SKIP_REQUEST(msg)						\
  do {									\
    export_json_exception(&server_bag, &query_bag, (msg));		\
    goto loop_end; } while (0)

/* 
 * interrupt handler for SIGINT 
 */
// static variable eliminates need to allocate memory during interrupt handler
static char *sigint_message = NULL;

/*
 * set up the json reply in the static variable
 */
static
void 
OnSIGINT_init(int setup, int head)
{
  if (setup > 0) {
    char *initial_json, *final_json;
    char *header;
    json_t *jroot = json_new_object();
    
    header = head ? JSON_MIME : "";
    // put in expected name=value pairs
    json_insert_pair_into_object(jroot, 
				 "error", 
				 json_new_string("server killed during request"));
    json_insert_pair_into_object(jroot, 
				 "status", 
				 json_new_string("-1"));
    json_insert_pair_into_object(jroot, 
				 "runtime", 
				 json_new_number("0.0"));
    // convert to text
    json_tree_to_string(jroot, &initial_json);
    final_json = json_format_string(initial_json);
    free(initial_json);
    // install in the static variable
    sigint_message = malloc(strlen(final_json) + strlen(header) + 1);
    sprintf(sigint_message, "%s%s", header, final_json);
  } else {
    // free it
    free(sigint_message);
  }

}

/*
 * sigint handler
 */
static
int 
OnSIGINT(void *data)
{
  int rc;
  server_bag_t *server_bag = (server_bag_t *) data;

  V_printf(-1, "", "Recieved sigint, shutting down.\n");
  insert_log_message(server_bag, "Received sigint, shutting down.");
  // if we were serving a request, cancel it
  //   (the standard path would be export_json_exception(), but that
  //   routine has a lot of allocation, so we export a canned message instead.)
  if (server_bag->connfd != 0) {
    rc = write(server_bag->connfd, sigint_message, strlen(sigint_message));
    if (rc < 0)
      V_printf(-1, "", "Error writing final json buf\n");
    rc = close(server_bag->connfd);
    if (rc < 0)
      V_printf(-1, "", "Error closing json output channel\n");
  }
  server_teardown(server_bag);
  return 0;
}

/* 
 * main function 
 *
 * Three kinds of errors: fatal to the app, fatal to the request,
 * and via sigint.
 * - fatal to the app: explanations to stderr, and DoIt() itself
 * returns to the caller.  The exit is through DIE().  There is
 * presently no case where this happens during a request, but if 
 * it did, the request should be cancelled.
 * - fatal to the request: explanation as a JSON result to the 
 * requester, and return to get more requests.
 * - sigint: should have an explanation to stderr, and, if in the
 * middle of a request, JSON to the requester.
 */
int 
DoIt(void)
{
  // arguments
  const char *log_fn  = cmdparams_get_str(&cmdparams, "logfile", NULL);
  const char *port_fn = cmdparams_get_str(&cmdparams, "portfile",NULL);
  int port_given      = cmdparams_get_int(&cmdparams, "port",    NULL);
  int verbflag        = cmdparams_get_int(&cmdparams, "verb",    NULL);
  float timeout_time  = cmdparams_get_float(&cmdparams, "timeout", NULL);
  int max_recs_server = cmdparams_get_int(&cmdparams, "n",       NULL);
  int wantMimeHead    = cmdparams_get_int(&cmdparams, "m",       NULL);
  int wantRecInfo     = cmdparams_get_int(&cmdparams, "R",       NULL);
  int wantOwner       = cmdparams_get_int(&cmdparams, "o",       NULL);
  int preserveLog     = cmdparams_get_int(&cmdparams, "p",       NULL);
  // pieces of input line
  char ds[QUERY_MAX];
  char keylist[QUERY_MAX];
  char seglist[QUERY_MAX];
  char linklist[QUERY_MAX];
  // parameter singletons
  server_bag_t server_bag;
  query_bag_t query_bag;
  int max_recs;
  json_t *jroot = NULL;
  char *msg;
  int rc;

  // initialize the server-param container
  memset(&server_bag, 0, sizeof(server_bag));
  server_bag.connfd = 0; // explicitly indicate no request being processed
  server_bag.preserve_log = preserveLog;
  server_bag.head = wantMimeHead;
  server_bag.timeout_time = timeout_time * 60.0;  // convert to seconds
  if (!(server_bag.timeout_time > 0 && server_bag.timeout < 1e8))
    DIE("Given timeout (%g) is out of range", timeout_time);

  // Register function that will be called if this DRMS module catches SIGINT.
  //   (will print an error message to stderr if it fails.)
  OnSIGINT_init(1, wantMimeHead);
  CleanerData_t cleaner;
  cleaner.cb = (pFn_Cleaner_t) &OnSIGINT;
  cleaner.data = (void *) &server_bag;
  /* FIXME: From arta 30 sep 2013:
     For the cleaner code, you use drms_client_registercleaner() if the macro DRMS_CLIENT is defined, and use drms_server_registercleaner otherwise (not all functions have a nice wrapper that chooses the correct version of API function). */
#ifdef DRMS_CLIENT
  drms_client_registercleaner(drms_env, &cleaner);
#else
  drms_server_registercleaner(drms_env, &cleaner);
#endif

  // binds the port
  msg = server_setup(port_given, port_fn, log_fn, &server_bag);
  if (msg) 
    DIE("Could not start up server connections: %s", msg);

  // It is possible to daemonize this server, but it turned out to be
  // iffy because exit of the parent process will close the DRMS
  // server down.  It's disabled.
#ifdef DAEMONIZE
  // daemonize
  switch (fork()) {
  case -1:
    // error
    DIE("fork failed");
    break;
  default:
    // original process -- if we return normally, the DRMS server
    // associated with this job will be shut down when DoIt returns
    close(server_bag.sock);
    _exit(1); // seems to work
    // kill(getpid(), 9); // also seems to work
    return 0;
    break;
  case 0:
    // daemon process
    break;
  }
#endif

  /*
   * begin serving queries
   */
  V_printf(verbflag, "", "Serving on port %d\n", server_bag.port);
  rc = listen(server_bag.sock, 5); 
  if (rc < 0)
    DIE("socket listen() failed");
  server_bag.running = 1;  // op=exit will reset this
  while (server_bag.running) {
    // ensure query summary is valid
    memset(&query_bag, 0, sizeof(query_bag));
    query_bag.status = 1; // assume failure so SKIP_REQUEST works

    // ensure the json message is valid
    jroot = json_new_object();

    // None of the json_* calls (and others) are error-checked.  
    // As a band-aid, alloc a large block here to see if we have
    // headroom, if not, fail.
    if ((msg = malloc(4*1024*1024)) == NULL) {
      V_printf(-1, "", "Quitting: out of memory\n");
      // log this, because the server has started
      insert_log_message(&server_bag, "Quit: out of memory");
      // close the server
      server_teardown(&server_bag);
      DIE("Out of memory");
    } else {
      free(msg);
    }
    
    // get a request
    //  (sets up server_bag to indicate request is in process)
    msg = input_query(&server_bag, &query_bag);
    if (server_bag.timeout) {
      // set up to skip request and exit
      query_bag.status = 0;   // request OK
      server_bag.running = 0; // exit now
      strcpy(query_bag.op, "timeout");  // for logging
      SKIP_REQUEST("Server exiting due to inactivity");
    } else if (msg != NULL) {
      SKIP_REQUEST(msg);
    }
    V_printf(verbflag>1, "", "Input query will go to fd = %d\n", server_bag.connfd);

    // break request (query_bag.line) into parts
    msg = parse_line(&query_bag, max_recs_server, keylist, seglist, linklist);
    if (msg)
      SKIP_REQUEST(msg);
    strcpy(ds, query_bag.ds); // we modify ds in-place below (not query_bag.ds)
    max_recs = query_bag.max_recs; // abbreviation

    V_printf(verbflag, "\t", "Running: %s ds=`%s' (K=`%s', S=`%s', L=`%s')\n", 
	     query_bag.op, query_bag.ds, keylist, seglist, linklist);

    /*******************************************************************
     *** exit
     *******************************************************************/
  if (strcmp(query_bag.op, "exit") == 0) {
    server_bag.running = 0; // set to exit the server
    query_bag.status = 0;   // mark request as OK
    SKIP_REQUEST("Exiting server normally.");
  }
    /*******************************************************************
     *** series_struct  
     *******************************************************************/
  else if (strcmp(query_bag.op, "series_struct") == 0) {
    char *p, *emsg;
    DRMS_Record_t *rec;
    int status = 0;

    /* Only want keyword info so get only the template record 
       for drms series or first record for other data */
    if (p = index(ds, '[')) *p = '\0';
    if (p = index(ds, '{')) *p = '\0';
    // (don't need to free rec)
    rec = drms_template_record(drms_env, ds, &status);
    if (status == DRMS_ERROR_QUERYFAILED) {
      emsg = (char *) DB_GetErrmsg(drms_env->session->db_handle);
      if (!emsg) emsg = "problem with database query";
      SKIP_REQUEST(emsg);
    } else if (status != 0) {
      SKIP_REQUEST("series not found");
    }

    list_series_info(drms_env, rec, jroot, wantOwner);
    if (get_series_stats(rec, jroot) == DRMS_ERROR_QUERYFAILED) {
      emsg = (char *) DB_GetErrmsg(drms_env->session->db_handle);
      if (!emsg) emsg = "problem with database query";
      SKIP_REQUEST(emsg);
    }
    // send answer
    query_bag.status = 0; // OK
    export_json(&server_bag, &query_bag, jroot);
  }

    /*******************************************************************
     *** op == rs_summary
     *******************************************************************/
  else if (strcmp(query_bag.op, "rs_summary") == 0) {
    char *emsg;
    int count=0, status=0;
    int countlimit = abs(max_recs);
    char val[32]; // an integer
    /* get series count */
    char *bracket = index(ds, '{');
    if (bracket)
      *bracket = '\0';
    if (countlimit) {
      DRMS_RecordSet_t *recordset = drms_open_nrecords(drms_env, ds, max_recs, &status);
      if (status == DRMS_ERROR_QUERYFAILED) {
	emsg = (char *) DB_GetErrmsg(drms_env->session->db_handle);
	if (!emsg) emsg = "problem with database query";
	SKIP_REQUEST(emsg);
      } else if (status != 0 || (!recordset)) {
	SKIP_REQUEST("unable to open records");
      }
      count = recordset->n;
      drms_close_records(recordset, DRMS_FREE_RECORD);
    } else {
      count = drms_count_records(drms_env, ds, &status);
    }
    if (bracket)
      *bracket = '{';
    if (status == DRMS_ERROR_QUERYFAILED) {
      emsg = (char *) DB_GetErrmsg(drms_env->session->db_handle);
      if (!emsg) emsg = "problem with database query";
      SKIP_REQUEST(emsg);
    } else if (status != 0) {
      SKIP_REQUEST("series not found");
    }
    
    /* send the output json back to client */
    snprintf(val, sizeof(val), "%d", count);
    json_insert_pair_into_object(jroot, "count", json_new_number(val));
    // send to output
    query_bag.status = 0; // OK
    export_json(&server_bag, &query_bag, jroot);
    } 

    /*******************************************************************
     *** op == rs_list
     *******************************************************************/
  else if (strcmp(query_bag.op, "rs_list") == 0) {
    DRMS_RecordSet_t *recordset;
    DRMS_Record_t *rec, *template;
    DRMS_RecChunking_t cstat = kRecChunking_None;
    char seriesname[DRMS_MAXQUERYLEN];
    char *keys[1000];
    char *segs[1000];
    char *links[1000];
    int ikey, nkeys = 0;
    int iseg, nsegs = 0;
    int ilink, nlinks = 0;
    int irec, nrecs;
    char count[32]; // integer
    json_t **keyvals = NULL, **segvals = NULL, **segdims = NULL, **linkvals = NULL;
    json_t *recinfo = NULL;
    json_t **segcparms = NULL;
    json_t **segbzeros = NULL;
    json_t **segbscales = NULL;
    json_t *json_keywords = NULL;
    json_t *json_segments = NULL;
    json_t *json_links = NULL;
    char *emsg;
    int status = 0;
    int record_set_staged = 0;
    char *lbracket;
    int requireSUMinfo;

    /* Get template record */
    strcpy(seriesname, ds);
    lbracket = index(seriesname, '[');
    if (lbracket) *lbracket = '\0';
    template = drms_template_record(drms_env, seriesname, &status);
    if (status == DRMS_ERROR_QUERYFAILED) {
      emsg = (char *) DB_GetErrmsg(drms_env->session->db_handle);
      if (!emsg) emsg = "problem with database query";
      SKIP_REQUEST(emsg);
    } else if (status != 0) {
      SKIP_REQUEST("series not found");
    }

    /* Open record_set(s) */
    if (max_recs == 0) {
      //      recordset = drms_open_recordset (drms_env, in, &status);
      // temporarily reverting to drms_open_records until I can fix the problem with
      // not passing a segment-list ot drms_open_recordset().
      recordset = drms_open_records(drms_env, ds, &status);
    } else {
      // max_recs specified via "n=" parameter. 
      recordset = drms_open_nrecords(drms_env, ds, max_recs, &status);
    }
    if (status == DRMS_ERROR_QUERYFAILED) {
      emsg = (char *) DB_GetErrmsg(drms_env->session->db_handle);
      if (!emsg) emsg = "problem with database query";
      SKIP_REQUEST(emsg);
    } else if (status == DRMS_QUERY_TRUNCATED) {
      SKIP_REQUEST("Query truncated; too many records?");
    }
    if (status != 0 || (!recordset))
      SKIP_REQUEST("series not found.");
    nrecs = recordset->n;
    if (nrecs == 0) {
      // bail early to avoid edge case
      json_insert_pair_into_object(jroot, "count", json_new_number("0"));
      query_bag.status = 0; // OK
      export_json(&server_bag, &query_bag, jroot);
      // jroot is freed at loop end
      goto loop_end;
    }
  
    /* get list of keywords to print for each record */
    /* Depending on the set of keywords to print, we will know whether or not 
     * we need to call SUM_infoEx(). Here's the list of keys that necessitate 
     * a SUM_infoEx() call:
     *  
     *   *size*
     *   *online*
     *   *retain*
     *   *archive* 
     */
    requireSUMinfo = 0;
    nkeys = 0;

    if (strlen(keylist) > 0) { 
      /* get specified list */
      char *thiskey;
      for (thiskey = strtok(keylist, ","); thiskey; thiskey = strtok(NULL, ",")) {
	if (strcmp(thiskey,"**NONE**") == 0) {
	  nkeys = 0;
	  break;
	} else if (strcmp(thiskey, "**ALL**") == 0) {
          DRMS_Keyword_t *key;
          HIterator_t *last = NULL;
          while ((key = drms_record_nextkey(template, &last, 0)))
	    if (!drms_keyword_getimplicit(key))
	      keys[nkeys++] = strdup(key->info->name);
          if (last) hiter_destroy(&last);
	} else {
	  keys[nkeys++] = strdup(thiskey);
	}
        if (strcmp(thiskey, "*size*"   ) == 0 || 
            strcmp(thiskey, "*online*" ) == 0 || 
            strcmp(thiskey, "*retain*" ) == 0 || 
            strcmp(thiskey, "*archive*") == 0)
	  requireSUMinfo = 1;
      }
    }
    /* place to put an array of keyvals per keyword */
    if (nkeys) {
      keyvals = (json_t **) malloc(nkeys * sizeof(json_t *));
      for (ikey = 0; ikey < nkeys; ikey++)
	keyvals[ikey] = json_new_array();
    } else {
      keyvals = NULL;
    }
  
    /* get list of segments to show for each record */
    nsegs = 0;
    if (strlen(seglist) > 0) { 
      char *thisseg;
      for (thisseg = strtok(seglist, ","); thisseg; thisseg = strtok(NULL, ",")) {
	if (strcmp(thisseg, "**NONE**") == 0) {
	  nsegs = 0;
	  break;
	} else if (strcmp(thisseg, "**ALL**") == 0) {
          DRMS_Segment_t *seg;
          HIterator_t *last = NULL;
          while ((seg = drms_record_nextseg(template, &last, 0)))
            segs[nsegs++] = strdup(seg->info->name);
          if (last) hiter_destroy(&last);
	} else {
	  segs[nsegs++] = strdup(thisseg);
	}
      }
    }
    /* place to put an arrays of segvals and segdims per segment */
    if (nsegs) {
      segvals    = (json_t **) malloc(nsegs * sizeof(json_t *));
      segdims    = (json_t **) malloc(nsegs * sizeof(json_t *));
      segcparms  = (json_t **) malloc(nsegs * sizeof(json_t *));
      segbzeros  = (json_t **) malloc(nsegs * sizeof(json_t *));
      segbscales = (json_t **) malloc(nsegs * sizeof(json_t *));
      for (iseg=0; iseg < nsegs; iseg++) {
	segvals[iseg]    = json_new_array();
	segdims[iseg]    = json_new_array();
	segcparms[iseg]  = json_new_array();
	segbzeros[iseg]  = json_new_array();
	segbscales[iseg] = json_new_array();
      }
    } else {
      segvals = segdims = segcparms  = segbzeros  = segbscales = NULL;
    }

    /* get list of links to print for each record */
    nlinks = 0;
    if (strlen(linklist) > 0) { 
      char *thislink;
      for (thislink = strtok(linklist, ","); thislink; thislink = strtok(NULL,",")) {
	if (strcmp(thislink, "**NONE**") == 0) {
	  nlinks = 0;
	  break;
	} else if (strcmp(thislink, "**ALL**") == 0) {
          DRMS_Link_t *link;
          HIterator_t *last = NULL;
          while ((link = drms_record_nextlink(template, &last)))
            links[nlinks++] = strdup(link->info->name);
          if (last)
	    hiter_destroy(&last);
	} else {
	  links[nlinks++] = strdup(thislink);
	}
      }
    }
    /* place to put an array of linkvals per keyword */
    if (nlinks) {
      linkvals = (json_t **) malloc(nlinks * sizeof(json_t *));
      for (ilink=0; ilink<nlinks; ilink++) 
	linkvals[ilink] = json_new_array();
    } else {
      linkvals = NULL;
    }

    /* ART - find out if we will be needing SUM info. If so, call drms_record_getinfo(). */
    if (requireSUMinfo) {
       drms_record_getinfo(recordset);
    }

    /* loop over set of selected records */
    for (irec = 0; irec < nrecs; irec++) {
      char recquery[DRMS_MAXQUERYLEN];
      char *jsonquery;
      json_t *recobj;

      if (max_recs == 0) {
        rec = drms_recordset_fetchnext(drms_env, recordset, &status, &cstat, NULL);
      } else {
        rec = recordset->records[irec];  /* pointer to current record */
        status = DRMS_SUCCESS;
      }

      if (wantRecInfo) {
        drms_sprint_rec_query(recquery,rec);
        jsonquery = string_to_json(recquery);
	recobj = json_new_object();
        json_insert_pair_into_object(recobj, "name", json_new_string(jsonquery));
        free(jsonquery);
      }
      /* now get keyword information */
      for (ikey=0; ikey<nkeys; ikey++) {
        DRMS_Keyword_t *rec_key_ikey; 
        json_t *thiskeyval = keyvals[ikey]; 
        json_t *val;
        char rawval[20000];
        char *jsonval;

        if (strcmp(keys[ikey],"*recnum*") == 0)
	  {
	  sprintf(rawval,"%lld",rec->recnum);
	  val = json_new_number(rawval);
	  }
        else if (strcmp(keys[ikey],"*sunum*") == 0)
	  {
	  sprintf(rawval,"%lld",rec->sunum);
	  val = json_new_number(rawval);
	  }
        else if (strcmp(keys[ikey],"*size*") == 0)
	  {
          char size[40];
	  SUM_info_t *sinfo = rec->suinfo;
          if (!sinfo)
	    val = json_new_string("NA");
	  else
            {
            sprintf(size,"%.0f", sinfo->bytes);
	    val = json_new_string(size);
            }
	  }
        else if (strcmp(keys[ikey],"*online*") == 0)
	  {
	  SUM_info_t *sinfo = rec->suinfo;
          if (!sinfo)
	    val = json_new_string("NA");
	  else
	    val = json_new_string(sinfo->online_status);
	  }
        else if (strcmp(keys[ikey],"*retain*") == 0)
	  {
	  SUM_info_t *sinfo = rec->suinfo;
          if (!sinfo)
	    val = json_new_string("NA");
	  else
	    {
            int y,m,d;
	    char retain[20];
            if (strcmp("N", sinfo->online_status) == 0)
              val = json_new_string("N/A");
            else
              {
              sscanf(sinfo->effective_date, "%4d%2d%2d", &y,&m,&d);
              sprintf(retain, "%4d.%02d.%02d",y,m,d);
	      val = json_new_string(retain);
              }
	    }
	  }
        else if (strcmp(keys[ikey],"*archive*") == 0)
	  {
	  SUM_info_t *sinfo = rec->suinfo;
          if (!sinfo)
	    val = json_new_string("NA");
	  else
            {
	    if(sinfo->pa_status == DAAP && sinfo->pa_substatus == DAADP)
              val = json_new_string("Pending");
            else
	      val = json_new_string(sinfo->archive_status);
            }
	  }
        else if (strcmp(keys[ikey], "*recdir*") == 0)
          { // get record directory
          char path[DRMS_MAXPATHLEN];
          if (!record_set_staged)
	    {
            drms_stage_records(recordset, 0, 0);
            record_set_staged = 1;
	    }
          drms_record_directory (rec, path, 0);
          jsonval = string_to_json(path);
          val = json_new_string(jsonval);
          free(jsonval);
          }
        else if (strcmp(keys[ikey], "*dirmtime*") == 0)
          { // get record dir last change date
	  struct stat buf;
          char path[DRMS_MAXPATHLEN];
          char timebuf[100];
          if (!record_set_staged)
	    {
            drms_stage_records(recordset, 0, 0);
            record_set_staged = 1;
	    }
          drms_record_directory (rec, path, 0);
          stat(path, &buf);
          sprint_ut(timebuf, buf.st_mtime + UNIX_EPOCH);
          jsonval = string_to_json(timebuf);
          val = json_new_string(jsonval);
          free(jsonval);
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
        else {
          rec_key_ikey = drms_keyword_lookup (rec, keys[ikey], 1); 
          if (!rec_key_ikey) {
	    V_printf(-1, "", "error, keyword not in series: %s\n",keys[ikey]);
	    // SKIP_REQUEST("Keyword not in series");
	    jsonval = string_to_json("Invalid KeyLink");
	  } else if (drms_ismissing_keyval(rec_key_ikey) && 
		     strcmp(keys[ikey],"QUALITY") != 0) {
	    jsonval = string_to_json("MISSING");
	  } else {
	    drms_keyword_snprintfval(rec_key_ikey, rawval, sizeof(rawval));
	    /* always report keyword values as strings */
	    jsonval = string_to_json(rawval);
	  }
	  val = json_new_string(jsonval);
	  free(jsonval);
	}
        json_insert_child(thiskeyval, val);
      }
  
      /* now show desired segments */
      int online = 0;
      for (iseg=0; iseg<nsegs; iseg++) 
        {
	DRMS_Segment_t *rec_seg_iseg = drms_segment_lookup(rec, segs[iseg]); 
        char *jsonpath;
        char *jsondims;
        char path[DRMS_MAXPATHLEN];
        json_t *thissegval = segvals[iseg]; 
        json_t *thissegdim = segdims[iseg]; 
        json_t *thissegcparms = segcparms[iseg];
        json_t *thissegbzero = segbzeros[iseg];
        json_t *thissegbscale = segbscales[iseg];
        if (rec_seg_iseg) {
	  int iaxis;
	  int naxis = rec_seg_iseg->info->naxis;
          char dims[100], dimval[20];

          // Get paths into segvals
          if (!record_set_staged)
	    {
            drms_stage_records(recordset, 0, 0);
            record_set_staged = 1;
	    }
          drms_record_directory (rec_seg_iseg->record, path, 0);
	  // if (!*path)
	  //   SKIP_REQUEST("Can not retrieve record path, SUMS may be offline");
          if (!*path) {
	    strcpy(path, "NoDataDirectory");
	  } else {
	    strncat(path, "/", DRMS_MAXPATHLEN);
	    strncat(path, rec_seg_iseg->filename, DRMS_MAXPATHLEN);
	  }
          jsonpath = string_to_json(path);
          json_insert_child(thissegval, json_new_string(jsonpath));
          free(jsonpath);
          online = strncmp(path, "/SUM", 4) == 0;

          // Get seg dimension info into segdims
          dims[0] = '\0';
          for (iaxis=0; iaxis<naxis; iaxis++) {
            if (iaxis)
              strcat(dims, "x");
            sprintf(dimval,"%d",rec_seg_iseg->axis[iaxis]);
            strcat(dims, dimval);
	  }
          jsondims = string_to_json(dims);
          json_insert_child(thissegdim, json_new_string(jsondims));
          free(jsondims);

          /* Print bzero and bscale (use format of the implicit keywords, %g) IFF 
           * the segment protocol implies the values (fits, fitsz, tas, etc.) */
          char keybuf[DRMS_MAXKEYNAMELEN];
          DRMS_Keyword_t *anckey = NULL;
          char *jsonkeyval = NULL;
          /* cparms */
          if (strlen(rec_seg_iseg->cparms)) {
	    jsonkeyval = string_to_json(rec_seg_iseg->cparms);
	    json_insert_child(thissegcparms, json_new_string(jsonkeyval));
	    free(jsonkeyval);
	  }
          /* bzero */
          snprintf(keybuf, sizeof(keybuf), "%s_bzero", segs[iseg]);
          anckey = drms_keyword_lookup(rec, keybuf, 1);
          if (anckey) {
	    drms_keyword_snprintfval(anckey, keybuf, sizeof(keybuf));
	    /* always report keyword values as strings */
	    jsonkeyval = string_to_json(keybuf);
	    json_insert_child(thissegbzero, json_new_string(jsonkeyval));
	    free(jsonkeyval);
	  }
          /* bscale */
          anckey = NULL;
          snprintf(keybuf, sizeof(keybuf), "%s_bscale", segs[iseg]);
          anckey = drms_keyword_lookup(rec, keybuf, 1);
          if (anckey) {
	    drms_keyword_snprintfval(anckey, keybuf, sizeof(keybuf));
	    /* always report keyword values as strings */
	    jsonkeyval = string_to_json(keybuf);
	    json_insert_child(thissegbscale, json_new_string(jsonkeyval));
	    free(jsonkeyval);
	  }
	} else {
          char *nosegmsg = "InvalidSegName";
          DRMS_Segment_t *segment = hcon_lookup_lower(&rec->segments, segs[iseg]);
          if (segment && segment->info->islink)
            nosegmsg = "BadSegLink";
          jsonpath = string_to_json(nosegmsg);
          json_insert_child(thissegval, json_new_string(jsonpath));
          free(jsonpath);
          jsondims = string_to_json("NA");
          json_insert_child(thissegdim, json_new_string(jsondims));
          free(jsondims);
	}
        }

      /* now show desired links */
      for (ilink=0; ilink<nlinks; ilink++) {
        DRMS_Link_t *rec_link = hcon_lookup_lower (&rec->links, links[ilink]); 
        DRMS_Record_t *linked_rec = drms_link_follow(rec, links[ilink], &status);
        char linkquery[DRMS_MAXQUERYLEN];
        if (linked_rec) {
          if (rec_link->info->type == DYNAMIC_LINK)
            drms_sprint_rec_query(linkquery, linked_rec);
          else
            sprintf(linkquery, "%s[:#%lld]", 
		    linked_rec->seriesinfo->seriesname, 
		    linked_rec->recnum);
          drms_close_record(linked_rec, DRMS_FREE_RECORD);
  
          json_t *thislinkval = linkvals[ilink]; 
          json_insert_child(thislinkval, json_new_string(linkquery));
	} else {
          json_t *thislinkval = linkvals[ilink]; 
          json_insert_child(thislinkval, json_new_string("Invalid_Link"));
	}
      }

      /* finish record info for this record */
      if (wantRecInfo) {
	recinfo = json_new_array();
	json_insert_pair_into_object(recobj, "online", json_new_number(online ? "1" : "0"));
        json_insert_child(recinfo, recobj);
      }
    }
  
    /* Finished.  Clean up and exit. */
    if (wantRecInfo)
      json_insert_pair_into_object(jroot, "recinfo", recinfo);

      json_keywords = json_new_array();
      for (ikey=0; ikey<nkeys; ikey++) 
        {
        json_t *keyname = json_new_string(keys[ikey]); 
        json_t *keyobj = json_new_object();
        json_insert_pair_into_object(keyobj, "name", keyname);
        json_insert_pair_into_object(keyobj, "values", keyvals[ikey]);
        json_insert_child(json_keywords, keyobj);
        }
      json_insert_pair_into_object(jroot, "keywords", json_keywords);

      json_segments = json_new_array();
      for (iseg=0; iseg<nsegs; iseg++) 
        {
        json_t *segname = json_new_string(segs[iseg]); 
        json_t *segobj = json_new_object();
        json_insert_pair_into_object(segobj, "name", segname);
        json_insert_pair_into_object(segobj, "values", segvals[iseg]);
        json_insert_pair_into_object(segobj, "dims", segdims[iseg]);
        json_insert_pair_into_object(segobj, "cparms", segcparms[iseg]);
        json_insert_pair_into_object(segobj, "bzeros", segbzeros[iseg]);
        json_insert_pair_into_object(segobj, "bscales", segbscales[iseg]);
        json_insert_child(json_segments, segobj);
        }
      json_insert_pair_into_object(jroot, "segments", json_segments);


      json_links = json_new_array();
      for (ilink=0; ilink<nlinks; ilink++) 
        {
        json_t *linkname = json_new_string(links[ilink]); 
        json_t *linkobj = json_new_object();
        json_insert_pair_into_object(linkobj, "name", linkname);
        json_insert_pair_into_object(linkobj, "values", linkvals[ilink]);
        json_insert_child(json_links, linkobj);
        }
      json_insert_pair_into_object(jroot, "links", json_links);

    drms_close_records(recordset, DRMS_FREE_RECORD);

    // insert metadata
    snprintf(count, sizeof(count), "%d", nrecs);
    json_insert_pair_into_object(jroot, "count", json_new_number(count));
    // send to output
    query_bag.status = 0; // OK
    export_json(&server_bag, &query_bag, jroot);

    // free
    if (nlinks > 0)
      free(linkvals);
    if (nkeys > 0)
      free(keyvals);
    if (nsegs > 0) {
      free(segvals);
      free(segdims);
      free(segcparms);
      free(segbzeros);
      free(segbscales);
    }
    for (iseg = 0; iseg < nsegs; iseg++)
      free(segs[iseg]);
    for (ikey = 0; ikey < nkeys; ikey++)
      free(keys[ikey]);
    for (ilink = 0; ilink < nlinks; ilink++)
      free(links[ilink]);
    } /* rs_list */

    /*******************************************************************
     *** op == unknown
     *******************************************************************/
  else {
    // unknown operator
    SKIP_REQUEST("unknown operator, skipping.");
  }
  // if we got to here, the request was handled normally
  query_bag.status = 0; 
  // point of re-entry for errors within the above steps
  loop_end:
  // these frees are applicable to all code paths
  if (jroot)
    json_free_value(&jroot);
  } // end outer loop
  // close socket, log file, etc.
  msg = server_teardown(&server_bag);
  if (msg)
    V_printf(-1, "", "Error closing server: %s\n", msg);
  // free pre-allocated error message
  OnSIGINT_init(0, wantMimeHead);
  // return OK
  return 0;
}


