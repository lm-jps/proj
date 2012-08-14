#include "jsoc.h"
#include "jsoc_main.h"

char *module_name = "hmi_fixCROTA2";

/*
 * This program updates CROTA2 based on results of the 2012 transit of Venus.
 * CROTA2 is corrected by adding either CAM1_DELTA or CAM2_DELTA
 * CALVERS0 is set to 1
 * DATE is updated.
 * A HISTORY line is added.
 *
 * call with hmi_fixCROTA2 ds=<seriesname> first=<first_time> last=<last_time>
 *
 * The program will operate on blocks of about 100d, 10d, 1d, or 6h as appropriate for
 * less than 25000 records in a chunk.
 *
 * The program requires (for now) that the first primekey be a slotted TIME.
 *
 * No query arguments other than the time range for each chunk are used.
 *
 */

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

#define NOTSPECIFIED "Not Specified"

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "ds", NOTSPECIFIED,  "data series to process."},
     {ARG_TIME, "first", NOTSPECIFIED, "first time to process."},
     {ARG_TIME, "last", NOTSPECIFIED, "last time to process."},
     {ARG_END}
};

#define     CAM1_DELTA	(0.0135 - 0.082596) // i.e. change 180.082596 to 180.0135 by adding CAM1_DELTA to INST_ROT and  CROTAT2
#define     CAM2_DELTA	(-0.0702)           // i.e. change 180.0 to 179.9298 by adding CAM1_DELTA to INST_ROT and  CROTAT2
#define     BLOCKSIZE   (25000)             // approx number of records to process in one call, will be rounded down to nice time.

int DoIt(void)
  {
  int status = DRMS_SUCCESS;
  DRMS_Record_t *inTemplate;
  DRMS_RecordSet_t *inRS, *outRS;
  DRMS_Keyword_t *pkey;
  int irec, nrecs, npkeys;
  const char *dsSeries = params_get_str(&cmdparams, "ds");
  char dsQuery[DRMS_MAXQUERYLEN];
  char *timekeyname;
  TIME t_first = params_get_time(&cmdparams, "first");
  TIME t_last = params_get_time(&cmdparams, "first");
  TIME t_step, t_epoch, t_block, t_start, t_stop;

  char history[4096];

  // Now find the prime time keyword name
  inTemplate = drms_template_record(drms_env, dsSeries, &status);
  if (status || !inTemplate) DIE("series not found");
  npkeys = inTemplate->seriesinfo->pidx_num;
  timekeyname = NULL;
  if (npkeys > 0)
    {
    int i;
    for (i=0; i<npkeys; i++)
        {
        pkey = inTemplate->seriesinfo->pidx_keywords[i];
        if (pkey->info->recscope > 1)
           pkey = drms_keyword_slotfromindex(pkey);
        if (pkey->info->type != DRMS_TYPE_TIME)
          continue;
        if(i > 0) DIE("Input series must have TIME keyword first, for now...");
        timekeyname = pkey->info->name;
        t_step = drms_keyword_getdouble(drms_keyword_stepfromslot(pkey), &status);
        if (status) DIE("problem getting t_step");
        t_epoch = drms_keyword_getdouble(drms_keyword_epochfromslot(pkey), &status);
        if (status) DIE("problem getting t_epoch");
        break;
        }
    }
  else
    DIE("Must have time prime key");

  t_block = t_step * BLOCKSIZE;
  if (t_block > 100*86400)
    t_block = 100*86500;
  else if (t_block > 10*86400)
    t_block = 10*86400;
  else if (t_block > 86400)
    t_block = 86400;
  else if (t_block > 6*3600)
    t_block = 6*3600;

  for (t_start = t_first; t_start <= t_last; t_start += t_block)
    {
    char first[100], last[100];
    t_stop = t_start + t_block - t_step;
    if (t_stop > t_last)
      t_stop = t_last;
    sprint_time(first, t_start, pkey->info->unit, 0);
    sprint_time(last, t_stop, pkey->info->unit, 0);

    sprintf(dsQuery, "%s[%s-%s]", dsSeries, first, last);
    
    inRS = drms_open_records(drms_env, dsQuery, &status);
    nrecs = inRS->n;
    if (status || nrecs == 0)
      DIE("No records found");
  
    outRS = drms_clone_records(inRS, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &status);
    nrecs = outRS->n;
    if (status || nrecs == 0)
      DIE("No records cloned");
    drms_close_records(inRS, DRMS_FREE_RECORD);
  
    for (irec=0; irec<nrecs; irec++)
      {
      int crotstat, camstat;
      DRMS_Record_t *rec = outRS->records[irec]; 
      double crota2 = drms_getkey_double(rec, "CROTA2", &crotstat);
      int camera = drms_getkey_int(rec, "CAMERA", &camstat);
      if (crotstat || camstat)
        {
        fprintf(stderr, "ERROR getkey CROTA2 or CAMERA bad, irec=%d, t_start=%s, crotstat=%d, camstat=%d\n",
           irec, first, crotstat, camstat);
        continue;
        }
  
      // The action is all between here and the end of the irec loop
      if (camera == 1)
        {
        crota2 += CAM1_DELTA;
        drms_setkey_double(rec, "INST_ROT", CAM1_DELTA);
        sprintf(history, "CROTA2 corrected by adding %6.4f degrees", CAM1_DELTA);
        }
      else
        {
        crota2 += CAM2_DELTA;
        drms_setkey_double(rec, "INST_ROT", CAM2_DELTA);
        sprintf(history, "CROTA2 corrected by adding %6.4f degrees", CAM2_DELTA);
        }
      
      drms_setkey_double(rec, "CROTA2", crota2);
      drms_setkey_int(rec, "CALVERS0", 1);
      drms_appendhistory(rec, history, 1);
      drms_setkey_time(rec, "DATE", CURRENT_SYSTEM_TIME);
      } //end of "irec" loop
  
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    fprintf(stdout, "%s CROTA2 fixed for %s through %s\n", dsSeries, first, last);
    } // end of time chunk loop

  return (DRMS_SUCCESS);
  } // end of DoIt
