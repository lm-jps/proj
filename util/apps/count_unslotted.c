/*

 count_unslotted - in=recset, counts records matching query.  Query may include the '@cadence' construct
 even in case of aia.lev1 and hmi.lev1 and other mostly empty slotted or unslotted series.

 Should be modified to simply count the lines in the input_recset file, or popen show_info -c in the general case.

*/
#include "jsoc.h"
#include "jsoc_main.h"
#include "atoinc.h"

char *module_name = "count_unslotted";

#define	DIE(msg) {fprintf(stderr,"%s  Status=%d\n",msg, status); return(status?status:1);}
#define DIE_get_recset(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return NULL;}

ModuleArgs_t module_args[] =
  {
    {ARG_STRING, "in", "NOTSPECIFIED", "input series. e.g 'mdi.fd_M_96m_lev18'"},
    {ARG_TIME, "epoch", "1993.01.01_TAI", "Reference epoch for '@' slots" },
    {ARG_END}
  };

char *get_input_recset(DRMS_Env_t *drms_env, char *inQuery);

int DoIt(void)
{
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *inRS = NULL;
  int status = DRMS_SUCCESS, nrecs; 
  char inQuery[DRMS_MAXQUERYLEN];
  char in[DRMS_MAXQUERYLEN];
  const char *ingiven = params_get_str(params, "in");
  char *new_in, *in_filename;
  strcpy(in, ingiven);
  new_in = get_input_recset(drms_env, in);
  if (new_in != in)
    {
    strcpy(in, new_in);
    in_filename = new_in+1;
    }

  inRS = drms_open_records(drms_env, in, &status);
  if (status)
     {
     fprintf(stderr,"Query is: %s\n",in);
     DIE("No input data found");
     }
  nrecs = inRS->n;
  drms_close_records(inRS, DRMS_FREE_RECORD); 
  if (in_filename)
    {
    unlink(in_filename);
    }

  printf("%d\n", nrecs);
  return DRMS_SUCCESS;
}


// In cases known to not have compact slotted series and cadence is specified
// generate explicit recordset list of closest good record to desired grid
// First get vector of times and quality
// Then if vector is not OK, quit.
// then: make temp file to hold recordset list
//       start with first time to define desired grid,
//       make array of desired times.
//       make empty array of recnums
//       search vector for good images nearest desired times
//       for each found time, write record query


char *get_input_recset(DRMS_Env_t *drms_env, char *inQuery)
  {
  static char newInQuery[DRMS_MAXSERIESNAMELEN+2];
  int epoch_given = cmdparams_exists(&cmdparams, "epoch");
  TIME epoch, t_epoch;
  DRMS_Array_t *data;
  DRMS_Record_t *inTemplate;
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
  char filename[DRMS_MAXSERIESNAMELEN];
  char *tmpdir;
  FILE *tmpfile;
  char newIn[DRMS_MAXQUERYLEN];
  char seriesname[DRMS_MAXQUERYLEN];
  char *lbracket;
  char *at = index(inQuery, '@');
  int npkeys;
  char *timekeyname;
  double t_step;

  strcpy(seriesname, inQuery);
  lbracket = index(seriesname,'[');
  if (lbracket) *lbracket = '\0';
  inTemplate = drms_template_record(drms_env, seriesname, &status);
  if (status || !inTemplate) DIE_get_recset("Input series can not be found");

  // Now find the prime time keyword name
  npkeys = inTemplate->seriesinfo->pidx_num;
  timekeyname = NULL;
  if (npkeys > 0)
    {
    int i;
    for (i=0; i<npkeys; i++)
        {
        DRMS_Keyword_t *pkey, *skey;
        pkey = inTemplate->seriesinfo->pidx_keywords[i];
        if (pkey->info->recscope > 1)
           pkey = drms_keyword_slotfromindex(pkey);
        if (pkey->info->type != DRMS_TYPE_TIME)
          continue;
        if(i > 0) DIE_get_recset("Input series must have TIME keyword first, for now...");
        timekeyname = pkey->info->name;
        t_step = drms_keyword_getdouble(drms_keyword_stepfromslot(pkey), &status);
        if (status) DIE_get_recset("problem getting t_step");
        t_epoch = drms_keyword_getdouble(drms_keyword_epochfromslot(pkey), &status);
        if (status) DIE_get_recset("problem getting t_epoch");
        }
    }
  else
    DIE_get_recset("Must have time prime key");
  epoch = epoch_given ? params_get_time(&cmdparams, "epoch") : t_epoch;

  if (at && *at && ((strncmp(inQuery,"aia.lev1[", 9)==0 ||
                    strncmp(inQuery,"hmi.lev1[", 9)==0 ||
                    strncmp(inQuery,"aia.lev1_nrt2[",14)==0 ||
                    strncmp(inQuery,"hmi.lev1_nrt[", 13)==0 ) ||
                   epoch_given))
    {
    char *ip=(char *)inQuery, *op=newIn, *p;
    long n, mul;
    while ( *ip && ip<at )
      *op++ = *ip++;
    ip++; // skip the '@'
    n = strtol(ip, &p, 10); // get digits only
    if (*p == 's') mul = 1;
    else if (*p == 'm') mul = 60;
    else if (*p == 'h') mul = 3600;
    else if (*p == 'd') mul = 86400;
    else 
      DIE_get_recset("cant make sense of @xx cadence spec");
    cadence = n * mul;
    ip = ++p;  // skip cadence multiplier
    while ( *ip )
      *op++ = *ip++;
    *op = '\0';
    half = cadence/2.0;
    sprintf(keylist, "%s,QUALITY,recnum", timekeyname);
    data = drms_record_getvector(drms_env, newIn, keylist, DRMS_TYPE_DOUBLE, 0, &status);
    if (!data || status)
      {
      fprintf(stderr, "status=%d\n", status);
      DIE_get_recset("getkey_vector failed\n");
      }
    nrecs = data->axis[1];
    irec = 0;
    t_this = (TIME *)data->data;
    dquality = (double *)data->data + 1*nrecs;
    drecnum = (double *)data->data + 2*nrecs;
    if (epoch_given)
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
        if (islot < nslots && this_t_diff <= t_diff)
          recnums[islot] = recnum;
        t_diff = fabs(t_now - t_want);
        }
    if (islot+1 < nslots)
      nslots = islot+1;  // take what we got.
    tmpdir = getenv("TMPDIR");
    if (!tmpdir) tmpdir = "/tmp";
    sprintf(filename, "%s/%sXXXXXX", tmpdir, module_name);
    mkstemp(filename);
    tmpfile = fopen(filename,"w");
    for (islot=0; islot<nslots; islot++)
      if (recnums[islot])
        fprintf(tmpfile, "%s[:#%lld]\n", seriesname, recnums[islot]);
    fclose(tmpfile);
    free(recnums);
    drms_free_array(data);
    sprintf(newInQuery,"@%s", filename);
    return(newInQuery);
    }
  else
	  return(inQuery);
  }
