static char rcsid[] = "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/spawnds.c,v 1.1 2009/02/03 22:18:47 production Exp $";
/* Routines to spawn a local dsds_svc or a pe_rpc to get to the global
 * dsds_svc.
*/
#include <stdio.h>
#include <string.h>
#include <pvm3.h>
#include "pe.h"

int find_pe_rpc(int mytid, int (*msg)(char *fmt, ...));

/* If a dsds_svc is not running, spawn it on the local machine. If already
 * running check that it's database connection is the one that we asked for.
 * Do likewise for dsds_rm.
 *   db = name of database for dsds_svc to connect to
 *   mytid = the pvm task id of the caller
 *   dsdstid = address of current/new task id of the spawned dsds_svc
 *   debugflg = spawn in debug mode if true
 *   msgf = function to print messages in printf fashion
 * Return non-0 on error.
*/
int spawn_dsds_local(char *db, int mytid, int *dsdstid, int debugflg,
                      int (*msgf)(char *fmt, ...))
{
  KEY *alist, *blist;
  struct taskinfo *taskp;
  char thishost[MAX_STR];
  char *args[2];
  char *dsdsdb;
  int i, ntask, dsds_rm_tid;

  if(pvm_tasks(0, &ntask, &taskp)) { /* get all tasks on the virt machine */
    (*msgf)("Error on pvm_taks() call for dsds_svc info\n");
    return(1);
  }
  for(i=0; i<ntask; i++) {
    if(!strcmp(taskp[i].ti_a_out, "dsds_svc")) {
      *dsdstid=taskp[i].ti_tid;
      break;
    }
  }
  if(i == ntask) {		/* dsds_svc not there, start it */
    gethostname(thishost,MAX_STR);
    args[0] = db;
    args[1] = NULL;
    if(debugflg)
      pvm_spawn("dsds_svc",args,PvmTaskDebug,thishost,1,dsdstid);
    else
      pvm_spawn("dsds_svc",args,PvmTaskHost,thishost,1,dsdstid);
    if(*dsdstid < 0) {
      (*msgf)("Error spawing dsds_svc on %s\n", thishost);
      return(1);
    }
    (*msgf)("dsds_svc tid=%x spawned on %s for db=%s\n", *dsdstid, thishost, db);
  }
  else {			/* dsds_svc is there. ck its database */
    alist=newkeylist();
    if((blist = (KEY *)call_dsds(&alist, REQDBNAME, *dsdstid, mytid, msgf, debugflg)) == NULL) {
      (*msgf)("Error requesting database name from dsds_svc\n");
      return(1);
    }
    dsdsdb = GETKEY_str(blist, "dsds_dbname");
    (*msgf)("dsds_svc tid=%x currently running for db=%s\n",*dsdstid,dsdsdb);
    freekeylist(&alist); freekeylist(&blist);
    if(strcmp(dsdsdb, db)) {
      (*msgf)("The dsds_svc must be halted if you really want to connect to db=%s\n", db);
      return(1);
    }
  }
  /* now do the same for dsds_rm */
  for(i=0; i<ntask; i++) {
    if(!strcmp(taskp[i].ti_a_out, "dsds_rm")) {
      dsds_rm_tid=taskp[i].ti_tid;
      break;
    }
  }
  if(i == ntask) {		/* dsds_rm not there, start it */
    gethostname(thishost,MAX_STR);
    args[0] = db;
    args[1] = NULL;
    if(debugflg)
      pvm_spawn("dsds_rm",args,PvmTaskDebug,thishost,1,&dsds_rm_tid);
    else
      pvm_spawn("dsds_rm",args,PvmTaskHost,thishost,1,&dsds_rm_tid);
    if(dsds_rm_tid < 0) {
      (*msgf)("Error spawing dsds_rm on %s\n", thishost);
      return(1);
    }
    (*msgf)("dsds_rm tid=%x spawned on %s for db=%s\n", dsds_rm_tid, thishost, db);
  }
  else {			/* dsds_rm is there. ck its database */
    /* No! dsds_rm may be deleting and you'll have to wait very long*/
    /* Its db should be right as it's always started up along with dsds_svc */
/*    alist=newkeylist();
/*    if((blist = (KEY *)call_dsds(&alist, REQDBNAME, dsds_rm_tid, mytid, msgf, debugflg)) == NULL) {
/*      (*msgf)("Error requesting database name from dsds_rm\n");
/*      return(1);
/*    }
/*    dsdsdb = GETKEY_str(blist, "dsds_dbname");
/*
/*    (*msgf)("dsds_rm tid=%x currently runningfor db=%s\n",dsds_rm_tid,dsdsdb);
/*    freekeylist(&alist); freekeylist(&blist);
/*    if(strcmp(dsdsdb, db)) {
/*      (*msgf)("The dsds_rm must be halted if you really want to connect to db=%s\n", db);
/*      return(1);
/*    }
*/
  }
  return(0);
}

/* Start a pe_rpc with pvm_spawn(), and wait for it's init message.
 * Return non-0 on error.
*/
int spawn_rpc_now(int *perpctid, int debugflg, int (*msgf)(char *fmt, ...))
{
  char thishost[MAX_STR];
  char *args[2];
  int initstat, dummy;

  gethostname(thishost,MAX_STR);
  /*args[0] = db;*/
  args[0] = NULL;
  args[1] = NULL;
  if(debugflg)
    pvm_spawn("pe_rpc",args,PvmTaskDebug,thishost,1,perpctid);
  else
    pvm_spawn("pe_rpc",args,PvmTaskHost,thishost,1,perpctid);
  if(*perpctid < 0) {
    (*msgf)("Error spawing pe_rpc on %s\n", thishost);
    return(1);
  }
  (*msgf)("pe_rpc tid=%x spawned on %s\n", *perpctid, thishost);
  /* get the pe_rpc init msg that it's ready to go */
  if(pvm_recv(*perpctid, MSGDSDS) < 0) {
    (*msgf)("Error receiving initialization response from pe_rpc\n");
    return(1);
  }
  pvm_upkint(&dummy, 1, 1);           /* don't care about request# */
  pvm_upkint(&initstat, 1, 1);        /* the pe_rpc init status */
  if(initstat) {
    (*msgf)("pe_rpc init err %d. Most likely pe_rpc_svc not running.\nContact CM@shoom.stanford.edu\n", initstat);
    return(1);
  }
  return(0);
}

/* See if a pe_rpc is running in our virtual machine, if not start one
 * on the local machine. If one is already running, verify that the db
 * that it is connected to is the one that we asked for.
 * Store the pvm tid of the pe_rpc in the given vrb perpctid.
 *   db = name of database for pe_rpc to connect to
 *   mytid = the pvm task id of the caller
 *   perpctid = address of current/new task id of the spawned pe_rpc
 *   debugflg = spawn in debug mode if true
 *   msgf = function to print messages in printf fashion
 * Return non-0 on error.
*/
int spawn_pe_rpc(char *db, int mytid, int *perpctid, int debugflg,
                  int (*msgf)(char *fmt, ...))
{
  KEY *alist, *blist;
  struct taskinfo *taskp;
  char *dsdsdb;
  int i, ntask;

  if(pvm_tasks(0, &ntask, &taskp)) { /* get all tasks on the virt machine */
    (*msgf)("Error on pvm_tasks() call for pe_rpc info\n");
    return(1);
  }
  for(i=0; i<ntask; i++) {
    if(!strcmp(taskp[i].ti_a_out, "pe_rpc")) {
      *perpctid=taskp[i].ti_tid;	/* current pe_rpc tid */
      break;
    }
  }
  if(i == ntask) {		/* pe_rpc not there, start it */
    if(spawn_rpc_now(perpctid, debugflg, msgf)) 
      return(1);
  }
  /* pe_rpc is there. check the dsds_svc database */
  alist=newkeylist();
  if((blist = (KEY *)call_dsds(&alist, REQDBNAME, *perpctid, mytid, msgf, debugflg)) == NULL) {
    (*msgf)("Error requesting database name from dsds_svc:\n");
    (*msgf)("I'm restarting pe_rpc as you may have a stale connection...\n");
    pvm_kill(*perpctid);	/* kill the current one */
/*    if(find_pe_rpc(mytid, msgf)) {	/* rogue pe_rpc present */
/*      /* !!TEMP */
/*      (*msgf)("If you know how the rogue pe_rpc got there please inform the CM.\n");
/*      (*msgf)("You don't have to do anything. This is informational only.\n");
/*      /*return(1);*/
/*    }
*/
    if(spawn_rpc_now(perpctid, debugflg, msgf)) return(1);
    if((blist = (KEY *)call_dsds(&alist, REQDBNAME, *perpctid, mytid, msgf, debugflg)) == NULL) {
      (*msgf)("***Error RE-requesting database name from dsds_svc\n");
      freekeylist(&alist);
      return(1);
    }
  }
  dsdsdb = GETKEY_str(blist, "dsds_dbname");
  (*msgf)("dsds_svc currently running for db=%s\n", dsdsdb);
  freekeylist(&alist); freekeylist(&blist);
  if(strcmp(dsdsdb, db)) {
    (*msgf)("The dsds_svc must be halted if you really want to connect to db=%s\n", db);
    return(1);
  }
  return(0);
}


/* This will determine if any pe_rpc is running on the current machine
 * for the current user.
 * It does this via  a ps command and a search for pe_rpc for the user. 
 * This is intended to be called by spawn_pe_rpc() after it has killed
 * any stale pe_rpc in the current pvm machine. This will then determine
 * if a rogue pe_rpc is running. There should be only one pe_rpc 
 * running for a user in a machine.
 * For now will output warning message and return an error status.
 * (After we get some experience this can maybe kill off any rogue
 * pe_rpc's for a user.)
 *
 * Called with the pvm task id of the caller and a ptr to a function 
 * to print a formatted msg (e.g. printf).
 * Returns non-0 on error.
*/

int find_pe_rpc(int mytid, int (*msg)(char *fmt, ...))
{
  FILE *fplog;
  char *thisuser, *cptr;
  char lfile[64], acmd[96], line[128], uline[64];
  int rpid;

  sprintf(lfile, "/tmp/find_pe_rpc.%x.log", mytid);
  sprintf(acmd, "ps -ef | grep pe_rpc  1> %s 2>&1", lfile);
  if(system(acmd)) {
    (*msg)("**find_pe_rpc():Can't execute %s.\n", acmd);
    return(-1);
  }
  if((fplog=fopen(lfile, "r")) == NULL) {
    (*msg)("**find_pe_rpc():Can't open %s to find any rogue pe_rpc\n", lfile);
    return(-1);
  }
  if(!(thisuser = (char *)getenv("USER"))) thisuser = "nouser";
  while(fgets(line, 128, fplog)) {       /* get ps lines */
    if(!(strstr(line, "pe_rpc \n"))) continue;
    if(strstr(line, "grep pe_rpc")) continue;
    sscanf(line, "%s %d", uline, &rpid); /* get user name & process id */
    if(!strcmp(thisuser, uline)) {	/* it's a rogue pe_rpc */
      (*msg)("***You have a rogue pe_rpc pid=%d running outside your pvm machine\n", rpid);
      fclose(fplog);
      return(-1);
    }
  }
  fclose(fplog);
  return (0);
}

/*
$Id: spawnds.c,v 1.1 2009/02/03 22:18:47 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/spawnds.c,v $
$Author: production $
*/
/* $Log: spawnds.c,v $
 * Revision 1.1  2009/02/03 22:18:47  production
 * initial
 *
 * Revision 1.15  2002/07/02  20:47:48  jim
 * call_dsds() now uses KEY **
 *
 * Revision 1.14  1997/09/04  17:08:56  jim
 * elim rouge pe_rpc check i.e dont call find_pe_rpc()
 *
 * Revision 1.13  1997/07/14  18:17:21  jim
 * add  CM@shoom.stanford.edu
 *
 * Revision 1.12  1997/04/30  23:30:51  jim
 * take out the last change
 *
 * Revision 1.11  1997/04/30  23:03:17  jim
 * don't kill the current pe_rpc that is being replaced
 *
 * Revision 1.10  1997/04/28  18:21:15  jim
 * continue now even if see rogue pe_rpc. just give informational msgs.
 *
 * Revision 1.9  1997/03/28  23:50:58  jim
 * add find_pe_rpc() to detect rogue pe_rpc
 *
 * Revision 1.8  1996/06/18  17:55:53  jim
 * add (char *fmt, ...)
 *
 * Revision 1.7  1996/05/17  17:26:41  jim
 * kill of any old pe_rpc
 *
 * Revision 1.6  1996/01/24  01:04:18  jim
 * dont query dsds_rm for it's database name
 *
 * Revision 1.5  1995/12/15  18:28:04  jim
 * fix spelling error
 *
 * Revision 1.4  1995/12/13  19:41:42  jim
 * *** empty log message ***
 * */
