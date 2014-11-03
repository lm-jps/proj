#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"

#define NOTSPECIFIED "***NOTSPECIFIED***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

ModuleArgs_t module_args[] = {
  {ARG_STRING, "ds", "lm_jps.pfss", "output data series"},
  {ARG_STRING, "ser_beg", NOTSPECIFIED, "serial number to begin"},
  {ARG_STRING, "ser_end", NOTSPECIFIED, "serial number to end (not saved)"},
  {ARG_STRING, "date_beg", NOTSPECIFIED, "date to begin"},
  {ARG_STRING, "date_end", NOTSPECIFIED, "date to end (not saved)"},
  {ARG_STRING, "spath", "/archive/pfss/kitrun48/surffield-serial", ""},
  {ARG_STRING, "bpath", "/archive/pfss/kitrun48/Bfield-serial", ""},
  {ARG_STRING, "sbase", "%s/kitrun048_%5.5d.sav", "surface base file name"},
  {ARG_STRING, "bbase", "%s/Bfield_%5.5d.sav", "B field base file name"},
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
        "ser_end=<serial number to stop> (not backup)\n"
        "date_beg=<first date to backup>\n"
        "date_end=<date to stop> (not backup)\n"
        "spath=<base path to surface field files>\n"
        "bpath=<base path to B field files>\n"
        "sbase=<surface field base file name>\n"
        "bbase=<B field field base file name>\n");
    return(1);
  }
  return(0);
}

int DoIt ()
{
  char *date_beg, *date_end, *ds, *spath, *bpath, *sbase, *bbase;
  char sname[512], bname[512], now_str[64];
  int ser_beg, ser_end, ser_num, status=0;
  int s_ok, b_ok, sp_ok, bp_ok; /* yr, mo, da, hr, mn; */
  time_t trec, tnow;
  DRMS_Record_t *outrec;
  DRMS_Segment_t *outseg;
  struct tm tm_rec;
  struct stat sbuf;

  if (nice_intro(0)) return(0);
  ser_beg = cmdparams_get_int(&cmdparams, "ser_beg", NULL);
  ser_end = cmdparams_get_int(&cmdparams, "ser_end", NULL);
  ds = strdup(cmdparams_get_str(&cmdparams, "ds", NULL));
  spath = strdup(cmdparams_get_str(&cmdparams, "spath", NULL));
  bpath = strdup(cmdparams_get_str(&cmdparams, "bpath", NULL));
  sbase = strdup(cmdparams_get_str(&cmdparams, "sbase", NULL));
  bbase = strdup(cmdparams_get_str(&cmdparams, "bbase", NULL));
if (verbose) fprintf(stdout, "'%s'\n'%s'\n", spath, bpath);
if (verbose) fprintf(stdout, "'%s'\n'%s'\n", sbase, bbase);
  if (ser_beg < 0) {
    time_t tr_b, tr_e;
    date_beg = strdup(cmdparams_get_str(&cmdparams, "date_beg", NULL));
    tr_b = sscan_time(date_beg) - UNIX_EPOCH + 218;
    ser_beg = (tr_b - 836179440)/21600;
    date_end = strdup(cmdparams_get_str(&cmdparams, "date_end", NULL));
    tr_e = sscan_time(date_end) - UNIX_EPOCH + 218;
    ser_end = (tr_e - 836179440)/21600;
    if (verbose) {
      fprintf(stderr, "trec: %ld %ld,ser_beg, end: %d, %d\n",tr_b,tr_e,ser_beg,ser_end);
    }
  }
  if (stat(spath, &sbuf)) {
    sp_ok = 0;
    fprintf(stderr, "Can not stat %s\n", spath);
//    return(1);
  } else sp_ok = 1;
  if (stat(bpath, &sbuf)) {
    bp_ok = 0;
    fprintf(stderr, "Can not stat %s\n", bpath);
//    return(1);
  } else bp_ok = 1;

//sleep(60);
  for (ser_num = ser_beg; ser_num < ser_end; ser_num++) {
    trec = 21600*ser_num + 836179440 - 8;
    if (verbose) {
      fprintf(stderr, "trec: %ld,ser_beg, end: %d, %d\n", trec, ser_beg, ser_end);
    }
    gmtime_r(&trec, &tm_rec);
    snprintf(sname, 512, sbase, spath, ser_num);
    snprintf(bname, 512, bbase, bpath, ser_num);
    outrec = drms_create_record(drms_env, ds, DRMS_PERMANENT, &status);
    if (status) DIE("Can't create output record");
//  drms_setkey_time(outrec, "model_date", (double) (trec + UNIX_EPOCH));
    sprint_time_ISO(now_str, trec + UNIX_EPOCH);
    drms_setkey_string(outrec, "model_date", now_str);
    if (verbose) fprintf(stderr, "model_date: %s\n", now_str);
    drms_setkey_int(outrec, "serial", ser_num);
    if (stat(sname, &sbuf)) {
      s_ok = 0;
      if (s_ok) fprintf(stderr, "Can not stat %s\n", sname);
    } else {
      s_ok = 1;
//    drms_setkey_time(outrec,"s_calc_date",(double)(sbuf.st_mtime+UNIX_EPOCH));
      sprint_time_ISO(now_str, sbuf.st_mtime + UNIX_EPOCH);
      drms_setkey_string(outrec, "s_calc_date", now_str);
    }
    if (stat(bname, &sbuf)) {
      b_ok = 0;
      if (b_ok) fprintf(stderr, "Can not stat %s\n", bname);
    } else {
      b_ok = 1;
//    drms_setkey_time(outrec,"b_calc_date",(double)(sbuf.st_mtime+UNIX_EPOCH));
      sprint_time_ISO(now_str, sbuf.st_mtime + UNIX_EPOCH);
      drms_setkey_string(outrec, "b_calc_date", now_str);
    }
    tnow = CURRENT_SYSTEM_TIME;
    drms_setkey_time(outrec, "backup_date", tnow);
    sprint_time_ISO(now_str, tnow);
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
    if (s_ok || b_ok) status = drms_close_record(outrec, DRMS_INSERT_RECORD);
    else status = drms_close_record(outrec, DRMS_FREE_RECORD);
    if(status) fprintf(stderr, "Error closing rec for ser #: %d.\n", ser_num);
  }
  return status;
}
