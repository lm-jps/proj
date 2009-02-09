/* call_drms_in.c
 * Gets the wd of the ARG_DATA_IN datasets
*/
#include <stdio.h>
#include <drms.h>
#include <drms_names.h>
#include <soi_key.h>
#include <printk.h>

extern DRMS_Env_t *drms_env;
extern void printkey (KEY *key);
extern void abortit();
extern int ampexflg;		// 1 = retrieve from tape

//#define out_namespace "dsds"
#define out_namespace "su_production"

//Take the su dir and convert to the corresponding ds_index (i.e. sunum)
int path2index(char *path, ulong *dsindex) {
  char pathx[DRMS_MAXPATHLEN];
  char *cptr, *cptr2;

  strcpy(pathx, path);
  cptr = strstr(pathx, "/D");	//find first /D
  if(!cptr) {
    printk("***Fatal Error: no /D in path %s\n", pathx);
    abortit(1);
  }
  //if there's a second /D then that is the orig one that we want
  cptr2 = strstr(cptr+1, "/D");
  if(cptr2) cptr = cptr2;
  cptr2 = index(cptr+1, '/');	//elim an other '/' after sunum
  if(cptr2) *cptr = NULL;
  cptr2 = cptr+2;		//skip the /D
  *dsindex = strtoul(cptr2, NULL, 0);
  return(0);
}


/*
* Here an example of a keylist received here:
* (the args are given by arg_data_in_0, arg_data_in_1, etc) 
*
*lago_tid:	KEYTYP_INT	0
*pds_host:	KEYTYP_STRING	d00
*dsds_uid:	KEYTYP_UINT	0
*alsoin_2_basename:	KEYTYP_STRING	fd_V_01h.029795
*alsoin_2_wd:	KEYTYP_STRING	/tmp40/jim/tmp/lev1.5/fd_V_01h/29795
*alsoin_2_rule:	KEYTYP_STRING	wd:{dbase}/{level}/{series}/{#%04d#series};bn:{series}.{#%06d#series}
*alsoin_1_basename:	KEYTYP_STRING	fd_V_01h.029794
*alsoin_1_wd:	KEYTYP_STRING	/tmp40/jim/tmp/lev1.5/fd_V_01h/29794
*alsoin_1_rule:	KEYTYP_STRING	wd:{dbase}/{level}/{series}/{#%04d#series};bn:{series}.{#%06d#series}
*alsoin_0_basename:	KEYTYP_STRING	fd_V_01h.029793
*alsoin_0_wd:	KEYTYP_STRING	/tmp40/jim/tmp/lev1.5/fd_V_01h/29793
*alsoin_0_rule:	KEYTYP_STRING	wd:{dbase}/{level}/{series}/{#%04d#series};bn:{series}.{#%06d#series}
*alsoin_2_dbase:	KEYTYP_STRING	/tmp40/jim/tmp
*alsoin_1_dbase:	KEYTYP_STRING	/tmp40/jim/tmp
*alsoin_0_dbase:	KEYTYP_STRING	/tmp40/jim/tmp
*alsoin_0_series_sn:	KEYTYP_INT	29793
*alsoin_0_fmt:	KEYTYP_STRING	%04d
*alsoin_0_incr:	KEYTYP_INT	1
*alsoin_0_lsn:	KEYTYP_INT	-1
*alsoin_0_fsn:	KEYTYP_INT	0
*alsoin_0_data:	KEYTYP_STRING	prog:mdi,level:lev1.5,series:fd_V_01h[29793-29795]
*alsoin_0_prog:	KEYTYP_STRING	mdi
*alsoin_0_level:	KEYTYP_STRING	lev1.5
*alsoin_0_series:	KEYTYP_STRING	fd_V_01h
*alsoin_0_series_range:	KEYTYP_STRING	29793,29794,29795
*alsoin_1_series_sn:	KEYTYP_INT	29794
*alsoin_1_fmt:	KEYTYP_STRING	%04d
*alsoin_1_incr:	KEYTYP_INT	1
*alsoin_1_lsn:	KEYTYP_INT	-1
*alsoin_1_fsn:	KEYTYP_INT	0
*alsoin_1_data:	KEYTYP_STRING	prog:mdi,level:lev1.5,series:fd_V_01h[29793-29795]
*alsoin_1_prog:	KEYTYP_STRING	mdi
*alsoin_1_level:	KEYTYP_STRING	lev1.5
*alsoin_1_series:	KEYTYP_STRING	fd_V_01h
*alsoin_1_series_range:	KEYTYP_STRING	29793,29794,29795
*alsoin_2_series_sn:	KEYTYP_INT	29795
*alsoin_2_fmt:	KEYTYP_STRING	%04d
*alsoin_2_incr:	KEYTYP_INT	1
*alsoin_2_lsn:	KEYTYP_INT	-1
*alsoin_2_fsn:	KEYTYP_INT	0
*alsoin_2_data:	KEYTYP_STRING	prog:mdi,level:lev1.5,series:fd_V_01h[29793-29795]
*alsoin_2_prog:	KEYTYP_STRING	mdi
*alsoin_2_level:	KEYTYP_STRING	lev1.5
*alsoin_2_series:	KEYTYP_STRING	fd_V_01h
*alsoin_2_series_range:	KEYTYP_STRING	29793,29794,29795
*alsoin_nsets:	KEYTYP_INT	3
*alsoin_dbase:	KEYTYP_STRING	/tmp40/jim/tmp
*arg_data_in_1:	KEYTYP_STRING	alsoin
*in_0_basename:	KEYTYP_STRING	
*in_0_wd:	KEYTYP_STRING	.
*in_basename:	KEYTYP_STRING	
*in_wd:	KEYTYP_STRING	.
*in_0_dbase:	KEYTYP_STRING	/tmp40/jim/tmp
*in_0_series_sn:	KEYTYP_INT	-1
*in_0_fmt:	KEYTYP_STRING	%04d
*in_0_incr:	KEYTYP_INT	1
*in_0_lsn:	KEYTYP_INT	-1
*in_0_fsn:	KEYTYP_INT	0
*in_0_data:	KEYTYP_STRING	prog:mdi_klist,level:copy,series:vw_V[-1]
*in_0_prog:	KEYTYP_STRING	mdi_klist
*in_0_level:	KEYTYP_STRING	copy
*in_0_series:	KEYTYP_STRING	vw_V
*in_0_series_range:	KEYTYP_STRING	-1
*in_fmt:	KEYTYP_STRING	%04d
*in_incr:	KEYTYP_INT	1
*in_lsn:	KEYTYP_INT	-1
*in_fsn:	KEYTYP_INT	0
*in_data:	KEYTYP_STRING	prog:mdi_klist,level:copy,series:vw_V[-1]
*in_prog:	KEYTYP_STRING	mdi_klist
*in_level:	KEYTYP_STRING	copy
*in_series:	KEYTYP_STRING	vw_V
*in_series_sn:	KEYTYP_INT	-1
*in_series_range:	KEYTYP_STRING	-1
*in_nsets:	KEYTYP_INT	1
*in_dbase:	KEYTYP_STRING	/tmp40/jim/tmp
*arg_data_in_0:	KEYTYP_STRING	in
*dlocalflg:	KEYTYP_INT	0
*db:	KEYTYP_STRING	mdi_2
*anotherout:	KEYTYP_STRING	prog:mdi,level:tmp,series:test[666]
*out:	KEYTYP_STRING	prog:mdi_raw,level:tlm,series:nxt2ser[3006,3007]
*alsoin:	KEYTYP_STRING	prog:mdi,level:lev1.5,series:fd_V_01h[29793-29795]
*in:	KEYTYP_STRING	prog:mdi_klist,level:copy,series:vw_V[-1]
*T_FIRST:	KEYTYP_STRING	'1999.11.06_00:01:00_TAI'
*/

KEY *call_drms_in(KEY *list, int dbflg) 
{
  static KEY *alist;
  char drmsname[MAX_STR];
  char cmd[128], buf[128];
  char path[DRMS_MAXPATHLEN] = {0};
  char *stmp, *cptr;
  DRMS_RecordSet_t *rset;
  DRMS_Record_t *rec;
  FILE *fin;
  char argname[MAX_STR], inname[MAX_STR], ext[MAX_STR];
  double dbytes;
  ulong dsindex;
  int i, loop, innsets, rstatus, touch;
  int restoreenv = 0;
  int num_ds = 0;                       /* total # of ds queried */
  int ntmp;

  alist = newkeylist();
  for(loop = 0; ; loop++) {
    sprintf(ext, "arg_data_in_%d", loop);
    if(!findkey(list, ext)) break;      /* all done, exit for(loop)*/
    strcpy(argname, GETKEY_str(list, ext)); /* e.g. "in" */
    strcpy(inname, argname); strcat(inname, "_nsets");
    innsets = getkey_int(list, inname);
    for(i=0; i < innsets; i++) {        /* do for each dataset */
      sprintf(ext, "_%d", i);
      strcpy(inname, argname); strcat(inname, ext); /* e.g. "in_0" */
      setkey_str(&alist, "inname", inname);
      sprintf(ext, "inname_%d", num_ds);
      setkey_str(&alist, ext, inname);		/* also make unique inname*/
      sprintf(ext, "%s_prog", inname);		//e.g. in_0_prog
      stmp = GETKEY_str(list, ext);
      if(!strcmp(stmp, "")) continue;			//not a ds name
      sprintf(drmsname, "%s__", stmp); 

//      setkey_str(&alist, ext, GETKEY_str(list, ext));
//      sprintf(ext, "%s_prog_sn", inname);
//      if(findkey(list, ext))
//        setkey_int(&alist, ext, getkey_int(list, ext));
      sprintf(ext, "%s_level", inname);
      stmp = GETKEY_str(list, ext);
      sprintf(drmsname, "%s%s__", drmsname, stmp);
      sprintf(ext, "%s_series", inname);
      stmp = GETKEY_str(list, ext);
      sprintf(drmsname, "%s%s", drmsname, stmp);
      sprintf(cmd, "echo \"%s\" | sed 's/\\\./_/g' | sed 's/-/__/g'", drmsname);
      //printf("cmd = %s\n", cmd);
      fin = popen(cmd, "r");
      fgets(buf, sizeof buf, fin);
      cptr = rindex(buf, '\n');
      if(cptr) *cptr = NULL;
      sprintf(ext, "%s_series_sn", inname);
      if(findkey(list, ext)) {
        ntmp = getkey_int(list, ext);
      }
      else ntmp = -1;
      sprintf(buf, "%s[%d]", buf, ntmp);
      sprintf(drmsname, "%s.%s", out_namespace, buf); 

      if(findkey(list, ext))
        setkey_int(&alist, ext, getkey_int(list, ext));
      if(findkey(list, "touch")) {
        touch = getkey_int(list, "touch");
        setkey_int(&alist, "touch", touch);
        restoreenv = drms_env->retention;
        drms_env->retention = touch;
      }
      if(findkey(list, "ampex_tid")) 
        setkey_int(&alist, "ampex_tid", getkey_int(list, "ampex_tid"));

    //printf("\nThe keylist for call_drms_in() is:\n"); //!!TEMP
    //keyiterate(printkey, list);

    rset = drms_open_records(drms_env, drmsname, &rstatus);
    if(rstatus) {
      printk("Can't do drms_open_records(%s)\n", drmsname);
      return(NULL);         
    }
    sprintf(ext, "%s_wd", inname);
    if(!rset || rset->n == 0) {
      printk("No prev ds\n");  
      setkey_str(&alist, ext, "");	//send back empty wd
    }
    else {
      drms_stage_records(rset, 1, 1);	//retrieve from tape if needed
      rec = rset->records[0];  /* !!TBD fix 0 when more than one */
      if(!ampexflg) {
        drms_record_directory (rec, path, 0);
      } else {
        printk("Possible wait for drms retrieve from tape...\n");
        drms_record_directory (rec, path, 1);
      }
      setkey_str(&alist, ext, path);
      dbytes = du_dir(path);
      sprintf(ext, "%s_bytes", inname);
      setkey_double(&alist, ext, dbytes);
      if(rstatus = path2index(path, &dsindex)) {
        //doesn't return here. called abortit()
      }
      sprintf(ext, "%s_ds_index", inname);
      setkey_ulong(&alist, ext, dsindex);
      //printk("In call_drms_in() ext=%s dsindex=%u\n", ext, dsindex);
      //keyiterate(printkey, alist);
    }
      num_ds++;
    }
    if(restoreenv != 0) drms_env->retention = restoreenv;
  }
  return(alist);
}
