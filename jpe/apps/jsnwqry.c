/*-----------------------------------------------------------------------------
 * /home/production/cvs/JSOC/proj/jpe/apps/jsnwqry.c
 *-----------------------------------------------------------------------------
 *
 *  jsnwqry.c to qry the database given prog, level, series and start series_num
 *  and end series_num.
 *
 * Just like jsnqry but also shows create_date, archive_tape and arch_tape_fn.
 *
 * Usage: jsnwqry -pprog -llevel -sseries -bstart_sn -eend_sn
 * where: -pprog = prog is prog name eg. mdi_eof
 *       -llevel = level is level name eg. lev0
 *       -sseries = series is series name eg. 00022000_01h
 *        For series alone wild card can be used for search eg.00022000%
 *        Note: % is the wild card character
 *       -bstart_sn -eend_sn = start_sn and end_sn are series number ranges
*/

#include <jsoc_main.h>
#include <sum_rpc.h>
#include <cmdparams.h>
#include <drms.h>
#include <drms_names.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <strings.h>
#include <soi_key.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>     //for alarm(2) among other things...
#include <printk.h>
#include <setjmp.h>
#include "pe.h"

//#define out_namespace "su_production"
#define out_namespace "dsds"

void printkey();
int soi_errno = NO_ERROR;

int start_sn, end_sn, next_sn, series_num;
int abort_active;		/* set while doing an abort */
int nocontrolc;                 /* don't honor ^C if set */
int debugflg;			/* run all pvm servers in the debug mode */
				/* also do keyiterate. Don't use w/tae */ 
int ccactive = 0;		/* ^C already in progress */
int dbxflg = 0;

char a_prog[32], a_level[32], a_series[64];
char line[1024];
char *username;
char *in_ds;
KEY *dslist;
KEY *alist;
SUM_t *sumhandle = NULL;

static struct timeval TIMEOUT = { 20, 0 };

// List of default parameter values. NOTE: use get_args().
ModuleArgs_t module_args[] = {
  {ARG_END}
};

CmdParams_t cmdparams;
char **argv = NULL;
char *pargv[12];		//copy argv[] strings due to free at end
int argc = 0;
char *module_name = "jpeq";      // Module name presented to DRMS
jmp_buf env;


void get_args(int argc, char *argv[])
{
	int c;
	if(argc != 6) {
		printf("Usage: jsnwqry -pprog -llevel -sseries -bstart_sn -eend_sn\n");
		printf("where: -pprog = prog is prog name eg. mdi_eof\n");
		printf("       -llevel = level is level name eg. lev0\n");
		printf("       -sseries = series is series name eg. 00022000_01h\n");
		printf("	For series alone wild card can be used for search eg.00022000%%\n");
 		printf("	Note: %% is the wild card character\n");
		printf("       -bstart_sn -eend_sn = start_sn and end_sn are series number ranges\n");
		exit(1);
	}
	while(--argc > 0 && (*++argv)[0] == '-') {
		while( c = *++argv[0]) 
		   switch (c) {
		   case 'p':
		     if(*++argv[0] != NULL) {
			strcpy((char *)a_prog, argv[0]);
		     }
			while(*++argv[0] != NULL);
			--argv[0];
			break;
		   case 'l':
		     if(*++argv[0] != NULL) {
			strcpy((char *)a_level, argv[0]);
		     }
			while(*++argv[0] != NULL);
			--argv[0];
			break;
		   case 's':
		     if(*++argv[0] != NULL) {
			strcpy((char *)a_series, argv[0]);
		     }
			while(*++argv[0] != NULL);
			--argv[0];
			break;
		   case 'b':
		     if(*++argv[0] != NULL) {
			start_sn = atoi(argv[0]);
		     }
		   	break;
		   case 'e':
		     if(*++argv[0] != NULL) {
			end_sn = atoi(argv[0]);
		     }
			break;	
		}
	}
}

/* Return ptr to "mmm dd hh:mm:ss". */
char *ljdatestring(void)
{
  time_t t;
  char *str;

  t = time(NULL);
  str = ctime(&t);
  str[19] = 0;
  return str+4;          /* isolate the mmm dd hh:mm:ss */
}


int DoIt()
{
  int wildflg = 0;
  int found = 0;
  int c, rstatus;
  char drmsname[128], cmd[128], buf[128], path[128], qcmd[128], rstring[128];
  char logfile[128];
  char *cptr, *cptr1;
  FILE *fin;
  DRMS_RecordSet_t *rset;
  DRMS_Record_t *rec;
  DRMS_SeriesInfo_t *seriesinfo;
  SUM_t *sumhandle = NULL;
  SUM_info_t *sinfo;
    uint64_t sunum;

  cmdparams_get_argv(&cmdparams, &argv, &argc);
  for(c=0; c < argc; c++) {
    pargv[c] = strdup(argv[c]);
  }
  printk_set(printf, printf);
  if(!(username = (char *)getenv("USER"))) username = "nouser";
  get_args(argc, pargv);
  if(strchr(a_series, '%')) { wildflg = 1; }
  //convert DSDS name to DRMS name
  sprintf(drmsname, "%s__%s__%s", a_prog, a_level, a_series);
      sprintf(cmd, "echo \"%s\" | sed 's/\\\./_/g' | sed 's/-/__/g'", drmsname);
      //printf("cmd = %s\n", cmd);
      fin = popen(cmd, "r");
      fgets(buf, sizeof buf, fin);
      fclose(fin);
      cptr = rindex(buf, '\n');
      if(cptr) *cptr = NULL;

  if(wildflg) {
    sprintf(qcmd, "select seriesname from %s.drms_series where seriesname like '%s.%s'", out_namespace, out_namespace, buf);
    sprintf(cmd, "echo \"%s\" | psql -h hmidb jsoc", qcmd); 
    fin = popen(cmd, "r");
    while(fgets(buf, sizeof buf, fin)) {
      if(strstr(buf, out_namespace)) {
        cptr = rindex(buf, '\n');
        if(cptr) *cptr = NULL;
        found = 1;
        break;  //we have a series name
      }
    }
    if(!found) return(0);
    cptr = strstr(buf, "__");
    if(!cptr) return(1);
    cptr1 = strstr(cptr+2, "__"); //start of real series name
    if(!cptr1) return(1);
    strcpy(a_series, cptr1+2);
  }

  if((sumhandle = SUM_open(NULL, NULL, printk)) == 0) {
    printk("Failed on SUM_open()\n");
    return(1);
  }

  while(1) {
    for(next_sn = start_sn; next_sn <= end_sn; next_sn++) {
      if(wildflg) sprintf(drmsname, "%s[%d]", buf, next_sn);
      else sprintf(drmsname, "%s.%s[%d]", out_namespace, buf, next_sn);
      rset = drms_open_records(drms_env, drmsname, &rstatus);
      if(rstatus) {
        printk("Can't do drms_open_records(%s)\n", drmsname);
        continue;
      }
      if(!rset || rset->n == 0) {
        //printk("No records for %s\n", drmsname); //don't do printout
        continue;
      }
      rec = rset->records[0];
      sunum = (uint64_t)(rec->sunum);
      
       /* SUM_info() is not supported by the new MT SUMS. Convert to SUM_infoArray(), 
        * providing an array of size one. A better solution would be to re-write the 
        * enclosing loop to batch the records. This way, many records could be
        * opened at once, and then a chunk of SUNUMs could be provided to 
        * SUM_infoArray(). However, this program isn't likely to be used again, 
        * so we'll leave the loop inefficient.
        */
        if(SUM_infoArray(sumhandle, &sunum, 1, printk)) {
          printk("Fail on SUM_info() for sunum=%u\n", rec->sunum);
          return(1);  
        }
        else {
            /* sumhandle->sinfo is the first node in a linked list of size one. */
          sinfo = sumhandle->sinfo;
        } 

      printk("prog:%s,level:%s,series:%s[%d] %s %s %d\n", 
	a_prog, a_level, a_series, next_sn, sinfo->creat_date, sinfo->arch_tape, sinfo->arch_tape_fn);
      drms_close_record(rec, DRMS_FREE_RECORD);
    }
    if(wildflg) {
      if(!fgets(buf, sizeof buf, fin)) {
        fclose(fin);
        break;
      }
      if(!strstr(buf, out_namespace)) {
        fclose(fin);
        break;
      }
      cptr = rindex(buf, '\n');
      if(cptr) *cptr = NULL;
      cptr = strstr(buf, "__");
      if(!cptr) break;
      cptr1 = strstr(cptr+2, "__"); //start of real series name
      if(!cptr1) break;
      strcpy(a_series, cptr1+2);
    }
    else break;
  }
  if(sumhandle)
    SUM_close(sumhandle, printk);
}
