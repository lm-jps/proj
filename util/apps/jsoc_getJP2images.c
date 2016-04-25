/**
   @defgroup jsoc_getJP2images jsoc_getJP2images get paths of existing jp2000 images
   @ingroup su_util

   @brief process a DRMS query for AIA data and return paths wrt JSOC /web/jsoc/htdocs/data/aia/images

   @par Synopsis:
   @code
   jsoc_getJP2images {t=0} {j=0} in=aia lev1 query
   @endcode

   This is a special purpose module that takes an aia lev1 query and return the paths
   to existing jpeg2000 images that match the query.  It uses the export module method
   to enable the "@cadence" style query for aia.lev1 data.

   The json output (j=1) is 3 objects, "images", "count", and "status".
   the "images" object is an array of individual image objects each of which
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

   @param in  The AIA lev1 record-set spec.

   @par Exit_Status:
   If json output is requested "count" and "status" keys will be provided.
   Status is non zero on failures.

   @par Example:
   get AIA data on 15m ticks.
   @code
   jsoc_getJP2images in='aia.lev1[2013.10.20/1d@15m]' -j 
   @endcode

   @par Example:
   Get timetag, wavelength, URL, nearest HMI timetag for some EUV data near 0 UT on one day, via cgi-bin.
   @code
   http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_getJP2images?in=aia.lev1_euv_12s[2015.08.27/40s][171,335,211]&j=1&t=1
   @endcode
   Response is:
   @code
   {"images":[
 {"time":"20150826_235958","wave":"211","url":"http://jsoc.stanford.edu/data/aia/images/2015/08/26/211/2015_08_26__23_59_58_62__SDO_AIA_AIA_211.jp2","HMItag":"20150827_000000"}
,{"time":"20150827_000013","wave":"335","url":"http://jsoc.stanford.edu/data/aia/images/2015/08/27/335/2015_08_27__00_00_13_63__SDO_AIA_AIA_335.jp2","HMItag":"20150827_000000"}
,{"time":"20150827_000034","wave":"171","url":"http://jsoc.stanford.edu/data/aia/images/2015/08/27/171/2015_08_27__00_00_34_34__SDO_AIA_AIA_171.jp2","HMItag":"20150827_000000"}
], "count":3, "status":0}
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

char *module_name = "jsoc_getJP2images";

#define DIE(msg) {printf("%s\n",(json ? "{\"status\":-1}" : ""));\
                  fflush(stdout);\
                  fprintf(stderr,"%s, status=%d\n",msg,status);\
                  return(-1);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "times", "NOT SPECIFIED",  "AIA times desired."},
     {ARG_INT, "wave", "0",  "AIA wavelength in Angstroms"},
     {ARG_STRING, "cadence", "0",  "AIA cadence spec, e.g. '15m' or '60s' style"},
     {ARG_STRING, "QUERY_STRING", "NOT SPECIFIED",  "Params from web GET"},
     {ARG_FLAG, "j", "0", "Generate json formatted output, plain ascii otherwise"},
     {ARG_FLAG, "t", "0", "Generate HMI timetag nearest to each aia image found"},
     {ARG_END}
};

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
  char *cadence, *p;
  char *in_filename = NULL;
  int json, wantHMItime;
  const char *inQuery;
  char Query[1024];
  int iwave;
  int step, jp2step;
  int mul;
  int rec_step;

  web_query = strdup (cmdparams_get_str (&cmdparams, "QUERY_STRING", NULL));
  from_web = strcmp (web_query, "NOT SPECIFIED") != 0;

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

  iwave = params_get_int(&cmdparams, "wave");
  json = params_get_int(&cmdparams, "j");
  wantHMItime = params_get_int(&cmdparams, "t");
  inQuery = (char *)params_get_str(&cmdparams, "times");
  cadence = (char *)params_get_str(&cmdparams, "cadence");
  int n = strtol(cadence, &p, 10); // get digits only
  if (*p == 's') mul = 1;
  else if (*p == 'm') mul = 60;
  else if (*p == 'h') mul = 3600;
  else if (*p == 'd') mul = 86400;
  n *= mul;  // number of seconds per cadence

  // cadence of UV is 24s, of VIS is 1 hour, of EUV is 12s
  // but jp2 images made for every second UV, every 3rd EUV, and every VIS image.
  // step is the seconds per possible image for desired wavelength
  step = (iwave == 1600 || iwave == 1700 ? 24: (iwave == 4500 ? 3600 : 12));
  // jp2step is number of possible images skipped for jp2000 image cadence
  jp2step = (iwave == 1600 || iwave == 1700 ? 2: (iwave == 4500 ? 1 : 3));
  rec_step = n/(step); // rec_step is number of records to skip for the desired cadence


  if (json)
    {
    printf("Content-type: application/json\n\n");
    }
  else
    {
    printf("Content-type: text/plain\n\n");
    }

  char *at = index(inQuery, '@');
  if (at)
    DIE("Use cadence param instead of @ syntax");
  sprintf(Query,"aia.lev1_nrt2[%s][?WAVELNTH=%d?]", inQuery, iwave);
fprintf(stderr,"query=%s\n",Query);
  inRS = drms_open_records(drms_env, Query, &status);
  if (status || inRS->n == 0)
    DIE("No input data found");

  nrecs = inRS->n;

  if (wave == 0)
    DIE("The wave parameter must be provided");

  if (json)
    {
    printf("{\"images\":[\n");
    }

  ngood = 0;
  irec = 0;

  int trysteps[] = {0, -1, 1, -2, 2, -3, 3};
  int try;
  int found = 0;
  while (irec < nrecs)
    {
    for (try=0; try<7; try++)
      {
      int tryrec;
      tryrec = irec + trysteps[try];
      if (tryrec < 0) tryrec = 0;
      if (tryrec >= nrecs) tryrec = nrecs - 1;
      inRec = inRS->records[tryrec];
      strncpy(date__obs, drms_getkey_string(inRec, "DATE__OBS", NULL), MAXSTR);
      strncpy(wave, drms_getkey_string(inRec, "WAVELNTH", NULL), MAXSTR);
      sscanf(date__obs, "%4d-%2d-%2dT%2d:%2d:%2d.%2d", &year,&month,&day,&hour,&minute,&sec,&fsec);
      sprintf(fname, "%4d_%02d_%02d__%02d_%02d_%02d_%02d__SDO_AIA_AIA_%s.jp2",
	  year,month,day,hour,minute,sec,fsec,wave);
      sprintf(path, "data/aia/images/%4d/%02d/%02d/%s/%s", year, month, day, wave, fname);
      sprintf(testpath, "/web/jsoc/htdocs/%s", path);
      if (access(testpath, R_OK) == 0)
        {
        found = 1;
        break;
        }
      }
    if (found)
      {
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
      sprintf(timetag, "%4d%02d%02d_%02d%02d%02d", year,month,day,hour,minute,sec);
      sprintf(showpath, "http://jsoc.stanford.edu/%s", path);
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
        printf("%s\n",testpath);
        }
      ngood++;
      }
    irec += rec_step; 
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

