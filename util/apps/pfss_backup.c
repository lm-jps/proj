#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"

#define NOTSPECIFIED "***NOTSPECIFIED***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

ModuleArgs_t module_args[] = {
  {ARG_STRING, "ds", "lm_jps.pfss", "output data series"},
  {ARG_STRING, "ser_beg", NOTSPECIFIED, "serial number to begin"},
  {ARG_STRING, "ser_end", NOTSPECIFIED, "serial number to end (not saved)"},
  {ARG_STRING, "spath", "/archive/pfss/kitrun48/surffield-serial", ""},
  {ARG_STRING, "bpath", "/archive/pfss/kitrun48/Bfield-serial", ""},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

char *module_name = "pfss_backup";

int verbose = 0;

void sprint_time_ISO (char *tstring, TIME t)
{
  sprint_at(tstring,t);
  tstring[4] = tstring[7] = '-';
  tstring[10] = 'T';
  tstring[19] = '\0';
}

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("pfss_backup -h -v ds=output_series ser_beg=sn1 ser_end=sn_stop\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        "ds=<output series name> default is lm_jps.pfss\n"
        "ser_beg=<first serial number to backup>\n"
        "ser_end=<serial number to stop> (not backup)\n");
    return(1);
  }
  return(0);
}

int DoIt ()
{
  char *ds, *spath, *bpath, sname[512], bname[512], now_str[64];
  int ser_beg, ser_end, ser_num, status=0;
  int s_ok, b_ok; /* yr, mo, da, hr, mn; */
  time_t tnow, trec;
  DRMS_Record_t *outrec;
  DRMS_Segment_t *outseg;
  struct tm tm_rec;
  struct stat sbuf;

  if (nice_intro(0)) return(0);
  tnow = time(NULL);
  ser_beg = cmdparams_get_int(&cmdparams, "ser_beg", NULL);
  ser_end = cmdparams_get_int(&cmdparams, "ser_end", NULL);
  ds = strdup(cmdparams_get_str(&cmdparams, "ds", NULL));
  spath = strdup(cmdparams_get_str(&cmdparams, "spath", NULL));
  bpath = strdup(cmdparams_get_str(&cmdparams, "bpath", NULL));

  for (ser_num = ser_beg; ser_num < ser_end; ser_num++) {
    trec = 21600*ser_num + 836179440;
    gmtime_r(&trec, &tm_rec);
    snprintf(sname, 512, "%s/kitrun048_%5.5d.sav", spath, ser_num);
    snprintf(bname, 512, "%s/Bfield_%5.5d.sav", bpath, ser_num);
    outrec = drms_create_record(drms_env, ds, DRMS_PERMANENT, &status);
    sprint_time_ISO(now_str, trec + UNIX_EPOCH);
    drms_setkey_string(outrec, "model_date", now_str);
    drms_setkey_int(outrec, "serial", ser_num);
    if (stat(sname, &sbuf)) {
      s_ok = 0;
      fprintf(stderr, "Can not stat %s\n", sname);
    } else {
      s_ok = 1;
      sprint_time_ISO(now_str, sbuf.st_mtime + UNIX_EPOCH);
      drms_setkey_string(outrec, "s_calc_date", now_str);
    }
    if (stat(bname, &sbuf)) {
      b_ok = 0;
      fprintf(stderr, "Can not stat %s\n", bname);
    } else {
      b_ok = 1;
      sprint_time_ISO(now_str, sbuf.st_mtime + UNIX_EPOCH);
      drms_setkey_string(outrec, "b_calc_date", now_str);
    }
    sprint_time_ISO(now_str, CURRENT_SYSTEM_TIME);
    drms_setkey_string(outrec, "backup_date", now_str);
    if (s_ok) {
      outseg = drms_segment_lookup(outrec, "surffield");
      if(drms_segment_write_from_file(outseg, sname)) {
        fprintf(stderr, "Can not write surface field segment for %s.\n", sname);
      }
    }
    if (b_ok) {
      outseg = drms_segment_lookup(outrec, "Bfield");
      if(drms_segment_write_from_file(outseg, bname)) {
        fprintf(stderr, "Can not write B field segment for %s.\n", bname);
      }
    }
    status = drms_close_record(outrec, DRMS_INSERT_RECORD);
    if(status) fprintf(stderr, "Error closing rec for ser #: %d.\n", ser_num);
  }
  return status;
}
