#include <string.h>
#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"

#define NOT_SPECIFIED "***Not Specified***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "dsinp", NOT_SPECIFIED, "Input series query"},
  {ARG_STRING, "dsout", NOT_SPECIFIED, "Output series"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

char *module_name = "aia_slot";
int verbose;

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("aia_slot {-h} {-v} dsinp=series_record_spec dsout=output_series\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        "dsinp=<recordset query> as <series>{[record specifier]} - required\n"
        "dsout=<series> - required\n");
    return(1);
  }
  return(0);
}

void sprint_time_ISO (char *tstring, TIME t)
{
  sprint_at(tstring,t);
  tstring[4] = tstring[7] = '-';
  tstring[10] = 'T';
  tstring[19] = '\0';
}

int DoIt ()
{
  int irec, nrecs, status, wl, fsn, qual, first = 1, cmdexp, explim;
  char *dsinp, *dsout, now_str[100];
  double tr_step;
  long long tr_index;
  TIME t_rec, t_obs, tr_epoch;
  DRMS_Record_t *inprec, *outrec;
  DRMS_RecordSet_t *inprs;
  DRMS_Keyword_t *inpkey = NULL, *outkey = NULL;
  DRMS_Segment_t *inpseg, *outseg;

  if (nice_intro(0)) return(0);

  dsinp = strdup(cmdparams_get_str(&cmdparams, "dsinp", NULL));
  dsout = strdup(cmdparams_get_str(&cmdparams, "dsout", NULL));
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("dsinp argument is required");
  if (strcmp(dsout, NOT_SPECIFIED)==0) DIE("dsout argument is required");

  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  for (irec=0; irec<nrecs; irec++) {
    inprec = inprs->records[irec];
    outrec = drms_create_record(drms_env, dsout, DRMS_PERMANENT, &status);
    if (status) DIE("cant create recordset");
    if (first) {
      tr_step = drms_getkey_double(outrec, "T_REC_step", &status);
      if (status) DIE("T_REC_step not found!");
      tr_epoch = drms_getkey_time(inprec, "T_REC_epoch", &status);
      if (status) DIE("T_REC_epoch not found!");
      first = 0;
    }
    sprint_time_ISO(now_str, CURRENT_SYSTEM_TIME);
    drms_setkey_string(outrec, "DATE", now_str);
    t_obs = drms_getkey_time(inprec, "T_OBS", &status);
    if (status) DIE("T_OBS not found!");
    cmdexp = drms_getkey_int(inprec, "AIMGSHCE", &status);
    if (status) DIE("cant get commanded exposure AIMGSHCE");
    wl = drms_getkey_int(inprec, "WAVELNTH", &status);
    if (!status) drms_setkey_int(outrec, "WAVELNTH", wl);
    switch (wl) {
      case 94:
      case 131:
      case 211:
      case 304:
      case 335:
      case 1600: explim = 2800; break;
      case 171:
      case 193: explim = 1931; break;
      case 1700: explim = 966; break;
      case 4500: explim = 483; break;
      default: DIE("Bad wavelength");
    }
    if (cmdexp < explim) { t_obs -= (1.0357*explim - cmdexp)*0.0005; }
    tr_index = (t_obs - tr_epoch)/tr_step;
    t_rec = tr_epoch + tr_index*tr_step;
    drms_setkey_time(outrec, "T_REC", t_rec);
    fsn = drms_getkey_int(inprec, "FSN", &status);
    if (!status) drms_setkey_int(outrec, "FSN", fsn);
    qual = drms_getkey_int(inprec, "QUALITY", &status);
    if (!status) drms_setkey_int(outrec, "QUALITY", qual);
    status = drms_link_set("lev1", outrec, inprec);
    drms_close_record(outrec, DRMS_INSERT_RECORD);
  }
  return 0;
}
