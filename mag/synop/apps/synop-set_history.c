#include <unistd.h>
#ifndef CVSTAG
#define CVSTAG "undefined"
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "jsoc_main.h"

// Defined in synop-saveparm.c
extern char *savestr;


char *cvsinfo_set_history = "cvsinfo: $Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/synop-set_history.c,v 1.1 2016/05/05 19:49:28 arta Exp $";

long long set_history(DRMS_Link_t *histlink)
{
  int status=0;
  int len;
  long long hold;
  TIME tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
  char path[DRMS_MAXPATHLEN+1];
  DRMS_Record_t *histrec = drms_create_record(drms_env, 
                                              histlink->info->target_series,
                                              DRMS_PERMANENT, 
                                              &status);

  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr,"ERROR: could not open a record in history dataseries %s, status = %d\n", histlink->info->target_series, status);
    return -1;
  }

  if ((len=readlink("/proc/self/exe", path, DRMS_MAXPATHLEN)) == -1)
  {
    fprintf(stderr, "ERROR: cannot locate binary\n");
    drms_close_record(histrec, DRMS_FREE_RECORD);
    return -1;
  }
  else
    path[len]='\0';

  status = drms_setkey_string(histrec, "MODPATH", path);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr,"ERROR: problem writing keyword MODPATH in history dataseries, status = %d\n", status);
    drms_close_record(histrec, DRMS_FREE_RECORD);
    return -1;
  }
  status = drms_setkey_string(histrec, "MODNAME", cmdparams.argv[0]);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr,"ERROR: problem writing keyword MODNAME in history dataseries, status = %d\n", status);
    drms_close_record(histrec, DRMS_FREE_RECORD);
    return -1;
  }
  status = drms_setkey_string(histrec, "ARGSUSED", savestr);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr,"ERROR: problem writing keyword ARGSUSED in history dataseries, status = %d\n", status);
    drms_close_record(histrec, DRMS_FREE_RECORD);
    return -1;
  }
  status = drms_setkey_string(histrec, "CVSTAG", CVSTAG);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr,"ERROR: problem writing keyword CVSTAG in history dataseries, status = %d\n", status);
    drms_close_record(histrec, DRMS_FREE_RECORD);
    return -1;
  }
  tnow = (double)time(NULL);
  tnow += UNIX_epoch;
  status = drms_setkey_time(histrec, "DATE", tnow);
  if (status != DRMS_SUCCESS)
  {
    fprintf(stderr,"ERROR: problem writing keyword DATE in history dataseries, status = %d\n", status);
    drms_close_record(histrec, DRMS_FREE_RECORD);
    return -1;
  }

  hold=histrec->recnum;
  drms_close_record(histrec, DRMS_INSERT_RECORD);
  return hold;
}
