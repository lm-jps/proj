static char rcsid[] = "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/start_pvmd.c,v 1.1 2009/01/23 21:30:46 production Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <pvm3.h>
#include "pe.h"

#define MAXSTR 128
#define PVMNOT "/usr/local/logs/dsds/PVM_NOT"

int ck_pe_rpc_dsds(int (*msg)());
int ck_cluster(char *hostname, char *user, int (*msg)());
int start_pvm_now(char *hostname, int (*msg)());

/* Starts a pvm daemon on the local host. 
 * Allows for a future pvm_addhosts() call to add aditya or fault, etc.. 
 * The file to start pvmd3 is built in /tmp/hostfile.$USER.
 * Uses $STAGING for ep= to define the path for the executables.
 * Registers the calling task with pvm and returns it tid. Returns 0 on error.
 *
 * msg = pointer to a function to print a vararg message.
 *
 * NEW:4Jun2002: First check if on a cluster machine (currently n00-n07) by
 * seeing if /etc/SGE/pvm_start exists.
 * If so then see if we can reach dsds_svc through a pvm/pe_rpc. If so then
 * we already have a good pvmd running, else we start one by calling 
 * /etc/SGE/pvm_start which will clean out any remnants of an old pe_rpc/pvmd 
 * and /tmp/pvm* files and start a new one outside of the SGE environment, 
 * so that it will not be killed off when the current job completes.
*/
int start_pvm(int (*msg)())
{
  char hostname[MAXSTR];
  char *user;
  int status;

  gethostname(hostname, MAXSTR);
  if(!(user = (char *)getenv("USER"))) user = "pvm";
  status = ck_cluster(hostname, user, msg);
  if(status == 0) return(0);
  return(start_pvm_now(hostname, msg));
}

/* Returns 0 on error, else the tid of the started pvmd.
*/
int start_pvm_now(char *hostname, int (*msg)())
{
  FILE *fphost;
  struct stat stbuf;
  int mytid, rint, napcnt, status;
  char filename[MAXSTR], cmd[MAXSTR];
  char *staging, *machine, *user;

  if((mytid = pvm_mytid()) < 0) {	/* pvmd not running */
    /* don't run if PVM_NOT file present  */
    if(stat(PVMNOT, &stbuf) != -1) { /* dont try to start a pvmd */
      (*msg)("*** pvm daemon not allowed. The SSSC is not available.\n");
      return(0);
    }
    /* new 15Sep98: delay a random time and check again in case
     * someone else has already started it (usually multiple
     * scripts still running after an sssc_stop).
    */
#ifdef __sgi
    srand(getpid());			/* seed for rand() */
    rint = rand();			/* 0 to 65535 */
    napcnt = rint/4369;			/* a sginap count 0 to 15 */
    sginap(napcnt);			/* delay 0 to 5 seconds */
#else
    sleep(1);
#endif
    if((mytid = pvm_mytid()) < 0) {	/* pvmd still not running */
      (*msg)("\nStarting a pvm daemon on %s\n", hostname);
      if(!(user = (char *)getenv("USER"))) user = "pvm";
      sprintf(filename, "/tmp/hostfile.%s", user);
      if((fphost=fopen(filename, "w")) == NULL) {
        (*msg)("Can't open %s to configure pvmd\n", filename);
        return(0);
      }
      if(!(staging = (char *)getenv("STAGING"))) {
        (*msg)("You must have an env $STAGING to find the executables for pvm\n");
        return(0);
      }
      if(!(machine = (char *)getenv("MACHINE"))) {
        (*msg)("You must have an env $MACHINE to find the executables for pvm\n");
        return(0);
      }
      fprintf(fphost, "#%s built by start_pvm()\n", filename);
      fprintf(fphost, "%s ep=%s/bin/_%s\n", hostname, staging, machine);
      /* put in aditya and fault, etc. for possible pvm_addhosts() later */
      /* NEW: this is no longer used, except for production possibly */
      /* adding tarax or sonar to its virtual machine. */
      /* NOTE: on 15May03 tarax became an alias for sonar */
      fprintf(fphost, "&tarax.Stanford.EDU ep=%s/bin/_sgi4\n", staging);
      fprintf(fphost, "&sonar.Stanford.EDU ep=%s/bin/_sgi4\n", staging);
/*******
      fprintf(fphost, "&fault.Stanford.EDU ep=%s/bin/_sgi4\n", staging);
      fprintf(fphost, "&shock.Stanford.EDU ep=%s/bin/_linux4\n", staging);
      fprintf(fphost, "&rick.Stanford.EDU ep=%s/bin/_linux4\n", staging);
      fprintf(fphost, "&phil.Stanford.EDU ep=%s/bin/_sgi4\n", staging);
      fprintf(fphost, "&yeti.Stanford.EDU ep=%s/bin/_sgi\n", staging);
      fprintf(fphost, "&shoom.Stanford.EDU ep=%s/bin/_sgi\n", staging);
      fprintf(fphost, "&aditya.Stanford.EDU ep=%s/bin/_sun4\n", staging);
      fprintf(fphost, "&shiver.Stanford.EDU ep=%s/bin/_sun4\n", staging);
      fprintf(fphost, "&soidb.Stanford.EDU ep=%s/bin/_sun4\n", staging);
      fprintf(fphost, "&quake.Stanford.EDU ep=%s/bin/_mips\n", staging);
********/
      fclose(fphost);
      sprintf(cmd, "pvmd3 %s&\n", filename);
      if(system(cmd)) {
        (*msg)("Error starting pvm daemon\n");
        return(0);
      }
      while(1) {
        sleep(2);			/* let it start */
        if((mytid = pvm_mytid()) > 0) break;
      }
    }
  } 
  return(mytid);
}

/* If on a cluster machine then see if a valid pvmd is already running, and
 * if not start one with an rsh to the machine for /etc/SGE/pvm_start.
 * Make sure that 2 different machines don't get in each others way by 
 * using an exclusive lock on a /tmp file.
 * The rsh must be done from some other machine (use sonar) so that the SGE
 * doesn't end up knowing about the pvmd and killing it when the calling pe
 * exits.(NEW: Don't need to do the rsh from sonar if redirect output ok.)
 * Return -1 if not on a cluster machine.
 * Return 0 if failure.
 * Return 1 if on cluster and pvmd is now running.
*/
int ck_cluster(char *hostname, char *user, int (*msg)())
{
  char cmd[MAXSTR];
  int status, fdes;


  if(!access("/etc/SGE/pvm_start", F_OK)) {	/* on a cluster machine */
    sprintf(cmd, "/tmp/pvm_lock.%s", user);
    if((fdes=open(cmd, O_WRONLY | O_CREAT, 0644)) == -1) {
      (*msg)("Can't open %s. Proceed anyway.\n", cmd);
    }
    if(fdes != -1) {
       lockf(fdes, F_LOCK, 0);		/* lock or block */
    }
    if(status = ck_pe_rpc_dsds(msg)) {	/* can't reach dsds_svc */
      //printf("Cannot reach dsds_svc, restarting\n");
      printf("Restarting...\n");
      if(status != 1) {			/* if other then no pvmd, exit it*/
	pvm_exit();
        system("echo y | /CM/bin/_linux4/pvm_halt");
      }
      (*msg)("Starting a pvmd on cluster %s\n", hostname);
/****old way to get around SGE killing the pvmd******
      sprintf(cmd, "rsh sonar \"rsh %s /etc/SGE/pvm_start&\"&", hostname);
*****/
      sprintf(cmd, "rsh %s /etc/SGE/pvm_start \>\& /tmp/pvm_start.%s.log&",
		hostname, user);
      (*msg)("cmd: %s\n", cmd); /* !!TEMP */
      if(system(cmd)) {
        (*msg)("Error starting pvm daemon with /etc/SGE/pvm_start\n");
        return(0);
      }
      while(1) {
        sleep(2);                       /* let it start */
        if(pvm_mytid() > 0) break;
      }
      sleep(2);				/* let pe_rpc start */
    }
    if(fdes != -1) {
      lockf(fdes, F_ULOCK, 0);		/* unlock so others can run */
      close(fdes);
    }
    return(1);
  }
  return(-1); 
}

/* See if a pvmd/pe_rpc is running in our virtual machine. 
 * Give error return if not.
 * If one is already running, verify that the dsds_svc can be contacted,
 * and print the db name that the dsds_svc is connected to.
 * Return 0 on success, else:
 *  1 = pvmd not running
 *  2 = pvmd internal error
 *  3 = pe_rpc not running
 *  4 = error contacting dsds_svc
*/
int ck_pe_rpc_dsds(int (*msg)())
{
  KEY *alist, *blist;
  struct taskinfo *taskp;
  char *dsdsdb;
  int i, ntask, perpctid;

  int mytid;

  if((mytid = pvm_mytid()) < 0) {       /* pvmd not running */
    (*msg)("pvmd not running\n");
    return(1);
  }
  if(pvm_tasks(0, &ntask, &taskp)) { /* get all tasks on the virt machine */
    (*msg)("Error on pvm_tasks() call for pe_rpc info\n");
    return(2);
  }
  for(i=0; i<ntask; i++) {
    if(!strcmp(taskp[i].ti_a_out, "pe_rpc")) {
      perpctid=taskp[i].ti_tid;        /* current pe_rpc tid */
      break;
    }
  }
  if(i == ntask) {              /* pe_rpc not there */
    (*msg)("pe_rpc not there\n");
    return(3);
  }
  /* pe_rpc is there. check the dsds_svc database */
  //For jpe, don't try to contact dsds_svc
//  alist=newkeylist();
//  if((blist = (KEY *)call_dsds(&alist, REQDBNAME, perpctid, mytid, msg, 0)) == NULL) {
//    (*msg)("Error requesting database name from dsds_svc:\n");
//    return(4);
//  }
//  freekeylist(&alist); freekeylist(&blist);
  return(0);
}

/*main(int argc, char *argv[])
/*{
/*  int mstat;
/*
/*  mstat=start_pvm(printf);
/*  printf("The pvm tid = %d\n", mstat);
/*}
*/

/*
$Id: start_pvmd.c,v 1.1 2009/01/23 21:30:46 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/start_pvmd.c,v $
$Author: production $
*/
/* $Log: start_pvmd.c,v $
 * Revision 1.1  2009/01/23 21:30:46  production
 * initial
 *
 * Revision 1.27  2006/03/27  22:49:42  jim
 * add mode 0644 to open()
 *
 * Revision 1.26  2004/11/09  23:25:25  rmunk
 * Updated file handling and removed extraneous printfs.
 *
 * Revision 1.25  2004/11/06 20:02:38  rmunk
 * removed C++ style comment that the SGI compiler chokes on.
 *
 * Revision 1.24  2004/11/05  18:43:26  rmunk
 * Changed file handling to speed up file locking.
 *
 * Revision 1.23  2003/05/15 17:38:02  jim
 * add comment on tarax alias to sonar
 *
 * Revision 1.22  2003/02/27 17:55:06  jim
 * remove adding host to the /tmp/hostfile.user
 *
 * Revision 1.21  2003/02/26  19:12:29  jim
 * change rick from _sgi to _linux4
 *
 * Revision 1.20  2002/08/29 23:20:28  jim
 * output to /tmp/pvm_start.%s.log and sleep(2)
 *
 * Revision 1.19  2002/08/21  23:14:00  jim
 * sprintf(cmd, "/tmp/pvm_lock.%s", user), user instead of hostname
 *
 * Revision 1.18  2002/08/14  17:08:20  jim
 * change shock to _linux4
 *
 * Revision 1.17  2002/07/02  20:50:03  jim
 * add exclusive lock on /tmp file
 *
 * Revision 1.16  2002/06/07 21:30:56  jim
 * don't need double rsh to "fool" SGE
 *
 * Revision 1.15  2002/06/07 00:04:42  jim
 * add code for cluster check and pvm_start call
 *
 * Revision 1.14  2001/08/24 20:56:35  jim
 * remove restriction of PVM_NOT only on tarax
 *
 * Revision 1.13  2001/05/02 21:04:28  jim
 * change phil to bin/_sgi4
 *
 * Revision 1.12  2000/03/20 22:14:42  jim
 * update host list. add sonar, remove flytrap, rumble
 *
 * Revision 1.11  1998/11/05  23:37:25  CM
 * add include for mips
 *
 * Revision 1.10  1998/10/14  18:48:55  jim
 * dont start pvmd on tarax if PVM_NOT file found
 *
 * Revision 1.9  1998/09/15  18:16:34  jim
 * add sginap delay
 *
 * Revision 1.8  1996/12/09  17:20:42  jim
 * change fault to _sgi4
 *
 * Revision 1.7  1995/07/24  20:28:39  jim
 * add yeti
 *
 * Revision 1.6  1995/06/30  17:37:30  jim
 * make tarax path _sgi4. add rick, phil, shoom
 *
 * Revision 1.5  1995/05/01  21:44:01  jim
 * add host shock and soidb
 *
 * Revision 1.4  1994/11/07  23:59:32  jim
 * added tarax.Stanford.EDU to hostfile
 *
 * Revision 1.3  1994/10/19  23:43:02  jim
 * put .user with hostfile name
 *
 * Revision 1.2  1994/09/16  18:20:43  jim
 * initial
 * */
