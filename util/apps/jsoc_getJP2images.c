/**
   @defgroup jsoc_getimages jsoc_getimages get paths of existing jp2000 images
   @ingroup su_util

   @brief process a DRMS query for AIA data and return paths wrt JSOC /web/jsoc/htdocs/data/aia/images

   @par Synopsis:
   @code
   jsoc_getimages {-j} in=aia lev1 query
   The output is on stdout in web acceptable format.  json if -j flag, else utf8 (ascii).
   @endcode

   This is a special purpose module that takes an aia lev1 query and return the paths
   to existing jpeg2000 images that match the query.  It uses the export module method
   to enable the "@cadence" style query for aia.lev1 data.

   The json output is 3 objects, "images", "count", and "status".
   the "images" object is an array of individual image objects eahc of which
   contains 3 or 4 strings: "time", "wave", "url", and optionally "HMItag".
   "time" is of the form yyyymmdd_hhmmss
   and wave is the AIA wavelength.  "url" is a complete url to fetch the image.
   If the flag "t=1" from an html call, or the flag -t on command line then
   the HMI time tag for the nearest 15m slot is shown.

   @par Flags:
   @c
   -j  create output in json with proper html header, else ascii with html header.
   -t  add timetag for nearest HMI 15m slot.

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in  The input data series.

   @par Exit_Status:
   If json output is requested "count" and "status" keys will be provided.
   Status is non zero on failures.

   @par Example:
   get AIA data on 15m ticks.
   @code
   jsoc_getimages in='aia.lev1[2013.10.20/1d@15m]' -j 
   @endcode

   @bug
   None known so far.

*/

#include "jsoc.h"
#include "jsoc_main.h"
#include <unistd.h>

static char x2c (char *what)
  {
  char digit;
  digit = (char)(what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
  digit *= 16;
  digit = (char)(digit + (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0')));
  return (digit);
  }

static void CGI_unescape_url (char *url)
  {
  int x, y;
  for (x = 0, y = 0; url[y]; ++x, ++y)
    {
    if ((url[x] = url[y]) == '%')
      {
      url[x] = x2c (&url[y+1]);
      y += 2;
      }
    }
  url[x] = '\0';
  }

char *module_name = "jsoc_getimages";

#define DIE(msg) {printf("%s\n",(json ? "{\"status\":-1}" : ""));\
                  fflush(stdout);\
                  fprintf(stderr,"%s, status=%d\n",msg,status);\
                  return(-1);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
     {ARG_STRING, "QUERY_STRING", "NOT SPECIFIED",  "Params from web GET"},
     {ARG_FLAG, "j", "0", "Generate json formatted output, plain ascii otherwise"},
     {ARG_FLAG, "t", "0", "Generate HMI timetag nearest to each aia image found"},
     {ARG_END}
};

const char *get_input_recset(DRMS_Env_t *drms_env, const char *in);

int DoIt(void)
  {
  int status = DRMS_SUCCESS;
  DRMS_RecordSet_t *inRS;
  DRMS_Record_t *inRec;
  char date__obs[MAXSTR];
  char wave[MAXSTR];
  char HMItag[MAXSTR];
  char fname[MAXSTR];
  char timetag[MAXSTR];
  char path[MAXSTR];
  char testpath[MAXSTR];
  char showpath[MAXSTR];
  int year, month, day, hour, minute, sec, fsec;
  int irec, nrecs, ngood;
  int from_web;
  char *web_query;
  char *inStr;
  char *in_filename = NULL;
  int json, wantHMItime;
  const char *inQuery;

  web_query = strdup (cmdparams_get_str (&cmdparams, "QUERY_STRING", NULL));
  from_web = strcmp (web_query, "Not Specified") != 0;

  if (from_web)
    {
    
    char *getstring, *p;
    CGI_unescape_url(web_query);
    getstring = strdup (web_query);
    for (p=strtok(getstring,"&"); p; p=strtok(NULL, "&"))
      {
      char *key=p, *val=index(p,'=');
      if (!val)
         {
         json = 1;  // assume in case of error
         DIE("Bad QUERY_STRING");
         }
      *val++ = '\0';
      cmdparams_set(&cmdparams, key, val);
      }
    free(getstring);
    }

  json = params_get_int(&cmdparams, "j");
  wantHMItime = params_get_int(&cmdparams, "t");
  inQuery = params_get_str(&cmdparams, "in");

  inStr = strdup(get_input_recset(drms_env, (char *)inQuery));
  if (!inStr || *inStr=='\0')
    DIE("Cant make special cadence recordset list file");
  inRS = drms_open_records(drms_env, inStr, &status);
  if (strcmp(inStr, inQuery) && *inStr == '@')
    unlink(inStr+1);
  if (status || inRS->n == 0)
    DIE("No input data found");

  nrecs = inRS->n;
  if (nrecs == 0)
    DIE("No records found");

  if (json)
    {
    printf("Content-type: application/json\n\n");
    printf("{\"images\":[\n");
    }
  else
    {
    printf("Content-type: text/plain\n\n");
    }

  for (ngood = irec = 0; irec < nrecs; irec++)
    {
    inRec = inRS->records[irec];
    if (wantHMItime)
      { // get nearest HMI 15m time slot.
      TIME t_rec = drms_getkey_time(inRec, "T_REC", NULL);
      t_rec = (t_rec + 450.0) / 900.0;
      int slot = floor(t_rec);
      t_rec = 900.0 * slot;
      char HMIt_rec[MAXSTR];
      sprint_time(HMIt_rec, t_rec, "TAI", 0);
      int y,m,d,h,M,s;
      sscanf(HMIt_rec, "%4d.%02d.%02d_%02d:%02d:%02d",&y,&m,&d,&h,&M,&s);
      sprintf(HMItag, "%4d%02d%02d_%02d%02d%02d",y,m,d,h,M,s);
      }
    strncpy(date__obs, drms_getkey_string(inRec, "DATE__OBS", NULL), MAXSTR);
    strncpy(wave, drms_getkey_string(inRec, "WAVELNTH", NULL), MAXSTR);
    sscanf(date__obs, "%4d-%2d-%2dT%2d:%2d:%2d.%2d", &year,&month,&day,&hour,&minute,&sec,&fsec);
    sprintf(fname, "%4d_%02d_%02d__%02d_%02d_%02d_%02d__SDO_AIA_AIA_%s.jp2",
	year,month,day,hour,minute,sec,fsec,wave);
    sprintf(timetag, "%4d%02d%02d_%02d%02d%02d", year,month,day,hour,minute,sec);
    sprintf(path, "data/aia/images/%4d/%02d/%02d/%s/%s", year, month, day, wave, fname);
    sprintf(testpath, "/web/jsoc/htdocs/%s", path);
    sprintf(showpath, "http://jsoc.stanford.edu/%s", path);
    if (access(testpath, R_OK) == 0)
      {
      if (json)
        {
        if (wantHMItime)
            printf("%s{\"time\":\"%s\",\"wave\":\"%s\",\"url\":\"%s\",\"HMItag\":\"%s\"}\n",
		(ngood ? "," : " "), timetag, wave, showpath, HMItag);
        else
            printf("%s{\"time\":\"%s\",\"wave\":\"%s\",\"url\":\"%s\"}\n",
		(ngood ? "," : " "), timetag, wave, showpath);
        }
      else
        {
        printf("%s\n",showpath);
        }
      ngood++;
      }
    }

  if (json)
    {
    printf("], \"count\":%d, \"status\":0}\n", ngood);
    }
  else
    {
    }
  drms_close_records(inRS, DRMS_FREE_RECORD);
  return (DRMS_SUCCESS);
  } // end of DoIt

// ----------------------------------------------------------------------

const char *get_input_recset(DRMS_Env_t *drms_env, const char *inQuery)
  {
  static char newInQuery[102];
  TIME epoch = (cmdparams_exists(&cmdparams, "epoch")) ? params_get_time(&cmdparams, "epoch") : 0;
  DRMS_Array_t *data;
  TIME t_start, t_stop, t_now, t_want, t_diff, this_t_diff;
  int status = 1;
  int nrecs, irec;
  int nslots, islot;
  long long *recnums;
  TIME *t_this, half;
  TIME cadence;
  double *drecnum, *dquality;
  int quality;
  long long recnum;
  char keylist[DRMS_MAXQUERYLEN];
  static char filename[100];
  char *tmpdir;
  FILE *tmpfile;
  char newIn[DRMS_MAXQUERYLEN];
  char seriesname[DRMS_MAXQUERYLEN];
  char *lbracket;
  char *at = index(inQuery, '@');
      
  if (at && *at && (strncmp(inQuery,"aia.lev1[", 9)==0 ||
                    strncmp(inQuery,"hmi.lev1[", 9)==0 ||
                    strncmp(inQuery,"aia.lev1_nrt2[",14)==0 ||
                    strncmp(inQuery,"hmi.lev1_nrt[", 13)==0 ))
    {
    char *ip=(char *)inQuery, *op=newIn, *p;
    long n, mul;
        
    char *segSpec = NULL;
    const char *psl = NULL;
    char *pbracket = NULL;
        
    psl = strchr(inQuery, '{');
    if (psl)
    {
        segSpec = strdup(psl + 1);
        
        if (!segSpec)
        {
            fprintf(stderr, "No memory.\n");
            return NULL;
        }
        
        pbracket = strchr(segSpec, '}');
        
        if (!pbracket)
        {
            fprintf(stderr, "Invalid segment specification.\n");
            return NULL;
        }
        
        *pbracket = '\0';
        
        if (!*segSpec)
        {
            fprintf(stderr, "Invalid segment specification.\n");
            return NULL;
        }
    }
        
    while ( *ip && ip<at )
      *op++ = *ip++;
    ip++; // skip the '@'
    n = strtol(ip, &p, 10); // get digits only
    if (*p == 's') mul = 1;
    else if (*p == 'm') mul = 60;
    else if (*p == 'h') mul = 3600;
    else if (*p == 'd') mul = 86400;
    else 
      {
      fprintf(stderr,"cant make sense of @xx cadence spec for aia or hmi lev1 data");
      return(NULL);
      }
    cadence = n * mul;
    ip = ++p;  // skip cadence multiplier
    while ( *ip )
      *op++ = *ip++;
    *op = '\0';
    half = cadence/2.0;
    sprintf(keylist, "T_OBS,QUALITY,recnum");
    data = drms_record_getvector(drms_env, newIn, keylist, DRMS_TYPE_DOUBLE, 0, &status);
    if (!data || status)
      {
      fprintf(stderr,"getkey_vector failed status=%d\n", status);
      return(NULL);
      }
    nrecs = data->axis[1];
    irec = 0;
    t_this = (TIME *)data->data;
    dquality = (double *)data->data + 1*nrecs;
    drecnum = (double *)data->data + 2*nrecs;
    if (epoch > 0.0)
      {
      int s0 = (t_this[0] - epoch)/cadence;
      TIME t0 = s0*cadence + epoch;
      t_start = (t0 < t_this[0] ? t0 + cadence : t0);
      }
    else
      t_start = t_this[0];
    t_stop = t_this[nrecs-1];
    nslots = (t_stop - t_start + cadence/2)/cadence;
    recnums = (long long *)malloc(nslots*sizeof(long long));
    for (islot=0; islot<nslots; islot++)
      recnums[islot] = 0;
    islot = 0;
    t_want = t_start;
    t_diff = 1.0e9;
    for (irec = 0; irec<nrecs; irec++)
        {
        t_now = t_this[irec];
        quality = (int)dquality[irec] & 0xFFFFFFFF;
        recnum = (long long)drecnum[irec];
        this_t_diff = fabs(t_now - t_want);
        if (quality < 0)
          continue;
        if (t_now <= (t_want-half))
          continue;
        while (t_now > (t_want+half))
          {
          islot++;
          if (islot >= nslots)
             break;
          t_want = t_start + cadence * islot;
          this_t_diff = fabs(t_now - t_want);
          t_diff = 1.0e8;
          }
        if (this_t_diff <= t_diff)
          recnums[islot] = recnum;
        t_diff = fabs(t_now - t_want);
        }
    if (islot+1 < nslots)
      nslots = islot+1;  // take what we got.
    strcpy(seriesname, inQuery);
    lbracket = index(seriesname,'[');
    if (lbracket) *lbracket = '\0';
    
    tmpdir = getenv("TMPDIR");
    if (!tmpdir) tmpdir = "/tmp";
    sprintf(filename, "%s/hg_patchXXXXXX", tmpdir);
    mkstemp(filename);
    tmpfile = fopen(filename,"w");
    for (islot=0; islot<nslots; islot++)
      if (recnums[islot])
      {
          if (!segSpec)
          {
              fprintf(tmpfile, "%s[:#%lld]\n", seriesname, recnums[islot]);
          }
          else
          {
              fprintf(tmpfile, "%s[:#%lld]{%s}\n", seriesname, recnums[islot], segSpec);
          }
      }
        
    if (segSpec)
        {
        free(segSpec);
        }
        
    fclose(tmpfile);
    free(recnums);
    drms_free_array(data);
    sprintf(newInQuery,"@%s", filename);
    return(newInQuery);
    }
  else
    return(inQuery);
  }
