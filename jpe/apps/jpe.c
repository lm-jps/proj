/*-----------------------------------------------------------------------------
 * /home/production/cvs/JSOC/proj/jpe/apps/jpe.c
 *-----------------------------------------------------------------------------
 *
 * This is the old pe.c from SOI MDI that is converted to run like the
 * original pe but with DRMS/SUMS interface instead of DSDS.
 * The original is found in /CM/src/pipe/pe.c
 *
 * HERE are the new jpe.c notes:
 * The following original pe options are deprecated in jpe:
 *  -d, -t, -L, db= , rcp= , pds=
 *
 ***************************************************************************
 * HERE are the original pe notes:
 * INPUTS:
 *Command line args of "[-v -d -t -L -A touch=# db=db_name
 *                      rcp=/dir pds=host_name]
 *                      map=map_file_name"
 *          -v is verbose mode, -t is tae mode with peui, -d is the debug
 *          mode where all servers are run under dbx, -L is local
 *          dsds_svc mode, or -A runs the ampex_svc in case it is needed for
 *          data retrieval from tape, db_name is the
 *          database name to connect to and overrides any DBNAME env vrbl,
 *          touch gives the #of days to set all the input datasets before they
 *          are eligible for deletion. This includes input datasets that are
 *          already online and those that are retrieved from tape.
 *          A value of "keep" will retain the input ds until a findkeep is
 *          run to release the kept datasets.
 *          If not given and no DBNAME env then uses mdi_2.
 *          The rcp= gives the dir to rcp the input ds to before running
 *          a module. This is for machines for which the NFS of the /PDS
 *          partitions is too slow.
 *          The pds= gives the name of the machine to assign local /PDS
 *          storage from. This is when there are multiple /PDS sets.
 *          If not specified, will use the env vrbl PDS_SET_HOST, else
 *          will use the local host that pe is executing on.
 *          See sample map file /home/soi/CM/src/pipe/map.examp
 * OUTPUTS: Messages to servers to process the dataset given in the map file.
 *          Messages to dsd_svc for storage and dataset name resolution.
 *          Messages to peui if -t flag.
 * RETURNS: exit(0) normal, exit(1) on fatal error
 * EXTERNALLY READ: Map file
 * EXTERNALLY MODIFIED: None
 *
 * DESCRIPTION:
 * This is the Pipeline Execution (PE) prototype. It will handle any number of
 * servers up to MAX_SERV different server types (i.e. different strategy
 * level module executables).
 * It takes as an input a map file name. The map file contains a list of
 * servers and various control directives. See src/pipe/map.txt for a spec.
 * The prog, level and series names for a dataset must all be present.
 * The servers and their order in the pipeline is determined by their order
 * in the map file.
 * PE requires pvm daemons to be running on all the hosts specified in the
 * map file. Deamons will be started automatically as required.
 * (PE once distributed the files for a dataset for each server amoungst the
 * number of host machined specified. This was taken out when a datacollection
 * was conceived and the production target machine was determined to be a
 * single SGI Power Challenge.)
 *
 */ 

#include <jsoc_main.h>
#include <sum_rpc.h>
#include <cmdparams.h>
#include <drms.h>
#include <drms_names.h>
#include <pvm3.h>
#include <soi_key.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>	//for umask(2)
#include <dirent.h>
#include <unistd.h>	//for alarm(2) among other things...
#include <printk.h>
#include <setjmp.h>
#include "soi_args.h"
#include "pe.h"
#include "dsds.h"

#define MONE (long)-1           /* needs to be long for IRIX64 compile */
#define PEMAILFILE "/tmp/pe.%s.%d.mail"
#define PKTSZ 1788		//size of VCDU pkt
#define MAXFILES 8192		//max # of file can handle in tlmdir
#define NUMTIMERS 8		//number of seperate timers avail
//#define IMAGE_NUM_COMMIT 12	//number of complete images until commit
#define IMAGE_NUM_COMMIT 2	//!!TEMP number of complete images until commit
#define TESTAPPID 0x199		//appid of test pattern packet
#define TESTVALUE 0xc0b		//first value in test pattern packet
#define MAXERRMSGCNT 10		//max # of err msg before skip the tlm file
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define mailable_users_num 6    /* #of users below that will generate mail */
static char *mailable_users[] = {"jprod","daemon","tprod","production",
                                "jeneen", "thailand"};
#define mailees_num 1   /* #below that are sent the daemon & production mail */
static char *mailees[] = {"prod2@solar2"};

// List of default parameter values. 
// NOTE: only used for cmdparams -H call. Normally use get_cmd()
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "map", NOTSPECIFIED, "map file with the pe directives"},
  {ARG_INT, "touch","-1","#of days to retain input datasets before deletable"}, 
  {ARG_FLAG, "A", "0", "retrieve from tape flag"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_END}
};

CmdParams_t cmdparams;
char **argv = NULL;
int argc = 0;
char *module_name = "jpe";	// Module name presented to DRMS
jmp_buf env;

int resp_dsds(int dsdstid);
double du_dir();

FILE *pemailfp;		// fp for pe ouput mail in /tmp/pe.pid.mail */
FILE *pelogfp;          /* fp for pe ouput log for this run */
SERVER stab[MAX_SERV+1];// one for ea server type plus a term
KEY *wd_dup_list;               /* all output wd are checked for dups */
KEY *wd_mk_list;                /* all output wd made are checked for dups */
KEY *JPElist;			/* info for archive groups */

static char datestr[32];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];
static DRMS_Record_t *rs;
extern DRMS_Env_t *drms_env;
extern KEY *call_drms_in(KEY *list, int dbflg);

unsigned int uid;               /* for jpe, assigned to pid */
unsigned int vcdu_24_cnt, vcdu_24_cnt_next;
int wdkey_int;                  /* make a uni key name in form_arg_data_out()*/
int verbose;
int whk_status;
int ampexflg;                   /* run with the ampex_svc */
int nocontrolc;                 /* don't honor ^C if set */
int archactive;                 /* 1 during call_archive(), -1 if ^C */
int ccactive;                   /* ^C already in progress */
int mailable;                   // set if we should send mail about aborts
int mailtrue;                   // set if ever call pemail()
int touchflg;			// days for touch
double reqbytes;                /* # of bytes requested by DSDSOUT= */
int rcpflg;                     /* set if rcp= given as input arg */
int dbxflg = 0;                 /* user defined to print keylists for debug */
int debugflg;                   /* !!TBD this is deprecated */
int reqmegs;                    /* set by get_cmd() if DSDSOUT=# in map file*/
int nowarn;                     /* set by get_cmd() if NOWARN=# in map file*/
int dsdsin;                     /* set by get_cmd() if any d=1 in map file */
int msg_id;                     /* assigned sequentially to server types */
int num_hosts;                  /* set by get_cmd() */
int effective_hosts;            /* set by get_cmd() changed by form_split()*/
int total_tlm_vcdu;
int total_missing_vcdu;
int dsds_tid;                   /* originally the tid of dsds_svc */
                                /* now the tid of pe_rpc that gets us to dsds*/
int pe_tid;                     /* tid of this pe */
int ampex_tid;                  /* tid of ampex_svc */
int lago_tid;                   /* tid of lago_svc (obsolete) */
int tae_tid;                    /* tid if we were spawned by tae ui */
int errmsgcnt, fileimgcnt;
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when ingest_lev0 is called for a restart
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
int firstfound = 0;		// set if see any file after startup
int firstalloc = 1;		// only assign out wd once in a archive group
int JPE_out_nsets;		// set 0 when start a new archive group
int ALRMSEC = 60;               // must get 2 in a row for no image timeout
char argvc[32], argindir[96], arglogfile[96], argoutdir[96];
char timetag[32];
char pchan[8];			// primary channel to listen to e.g. VC02 
char rchan[8];			// redundant channel to listen to e.g. VC10 
char stopfile[80];		// e.g. /usr/local/logs/lev0/VC04_stop
char tlmseriesname[128];	// e.g. hmi.tlm
char lev0seriesname[128];	// e.g. hmi.lev0
char tlmnamekey[128];		// shortened tlm file name for TLMDSNAM keyword
char tlmnamekeyfirst[128];	// for TLMDSNAM keyword for 1st file of image
char oldtlmdsnam[128];		// TLMDSNAM keyword from prev rec in db
char mailname[256];		// mail log file name
char pdshost[MAX_STR];          /* from pds= on cmd line */
char pe_map[256];		// name of the map file given as input
char dsdswd[MAX_STR];           /* working directory ret from Data Allocate */
char *username;			// from getenv("USER") 
char *rcpdir;                   /* from rcp= on cmd line */
char *t_first_arg;              /* set if modules has a T_FIRST arg */
char *t_last_arg;               /* set if modules has a T_LAST arg */


struct namesort {		// sorted file names in tlmdir 
  char *name;
};
typedef struct namesort NAMESORT;

// Setup global datestr[] like: 2008.07.14_08:29:31
char *do_datestr() {
  time_t tval;
  struct tm *t_ptr;

  tval = time(NULL);
  t_ptr = localtime(&tval);
  sprintf(datestr, "%d.%02d.%02d_%02d:%02d:%02d", 
	  (t_ptr->tm_year+1900), (t_ptr->tm_mon+1),
	  t_ptr->tm_mday, t_ptr->tm_hour, t_ptr->tm_min, t_ptr->tm_sec);
  return(datestr);
}

// Returns a time tag like  yyyy.mm.dd.hhmmss 
char *gettimetag()
{
  struct timeval tvalr;
  struct tm *t_ptr;

  gettimeofday(&tvalr, NULL);
  t_ptr = localtime((const time_t *)&tvalr);
  sprintf(timetag, "%04d.%02d.%02d.%02d%02d%02d",
        (t_ptr->tm_year+1900), (t_ptr->tm_mon+1), t_ptr->tm_mday, t_ptr->tm_hour, t_ptr->tm_min, t_ptr->tm_sec);
  return(timetag);
}


void BeginTimer(int n)
{
  gettimeofday (&first[n], NULL);
}

float EndTimer(int n)
{
  gettimeofday (&second[n], NULL);
  if (first[n].tv_usec > second[n].tv_usec) {
    second[n].tv_usec += 1000000;
    second[n].tv_sec--;
  }
  return (float) (second[n].tv_sec-first[n].tv_sec) +
    (float) (second[n].tv_usec-first[n].tv_usec)/1000000.0;
}


/* Outputs the variable format message (re: printf) to the pe log file.
*/
int pelog(char *fmt, ...)
{
  va_list args;
  char string[32768];

  va_start(args, fmt);
  vsprintf(string, fmt, args);
  if(pelogfp) {
    fprintf(pelogfp, string);
    fflush(pelogfp);
  }
  va_end(args);
  return(0);
}

/* Open the pe log file for this pe run
*/
void open_pelog(char *filename)
{
  if((pelogfp=fopen(filename, "w")) == NULL)
    fprintf(stderr, "Can't open the log file %s\n", filename);
}

void kill_pvm()
{
  int i;
  HDATA *hnext;

  for(i=0; i<MAX_SERV; i++) {           /* kill all servers on all hosts */
    if(stab[i].name==NULL)
      break;
    for(hnext=(HDATA *)gethnext(stab[i].hosts); hnext != NULL;
         hnext=(HDATA *)gethnext((HDATA *)MONE)) {
      if(hnext->tid > 0) {
        pvm_kill(hnext->tid);
        /*peklog("In kill_pvm(): tid = %x %s\n", hnext->tid, stab[i].name);*/
        pelog("%x\tn/a\tkilled\t0\t%s\t<NONE>\t<NONE>\n", hnext->tid,stab[i].name);
/* use pvm_kill() instead of sending a MSGKILL message */
/*        pvm_initsend(PvmDataDefault);
/*        if(pvm_send(hnext->tid, MSGKILL)) {
/*          printk("Error sending %s on %s its kill msg\n",
/*              stab[i].name, hnext->host_name);
/*        }
*/
      }
    }
  }
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit jpe w/ status = %d\n", stat);
  //if (pelogfp) fclose(pelogfp);
  //return(stat);
  longjmp(env, 2);	//get out of here
}

void sighandler(sig)
    int sig;
{
  KEY *list = newkeylist();
  int respcnt;

  if (signal(SIGINT, SIG_IGN) != SIG_IGN)
      signal(SIGINT, sighandler);
  if (signal(SIGTERM, SIG_IGN) != SIG_IGN)
      signal(SIGTERM, sighandler);
  if(archactive) {		/* don't abort during archive */
    printf("Ctrl-C during an archive update. Will abort when done...\n");
    archactive = -1;		/* tell call_archive() that it happened */
    return;
  }
  if(nocontrolc) {
    /* tape retrieval underway. send cancel to dsds_svc. */
    if(ccactive) {
      printf("A ^C is already active. I'm awaiting a dsds_svc reply...\n");
    }
    else {
      ccactive = 1;
      setkey_uint(&list, "dsds_uid", uid); /* give our OPEN uid */
      setkey_int(&list, "ampex_tid", ampex_tid);
      setkey_str(&list, "USER", username);
      printf("^C signal received by pe\n");
      printf("Sending TAPECANCEL request to dsds_svc...\n");
      if((list = call_dsds(&list,REQTAPECANCEL,dsds_tid,pe_tid,printf,debugflg)) == NULL) {
        printf("Can't REQTAPECANCEL with dsds_svc\n");  /* don't abort here */
      }
      else {
        if(getkeytype(list,"MSG_PEND")==KEYTYP_INT) { /* wait for completion */
          respcnt = 0;
          while(1) {
            if(resp_dsds(dsds_tid)) {     /* rets every timeout interval */
              printf("..."); respcnt++;
              if(respcnt == 20) {
                respcnt = 0;
                printf("\n");
              }
              continue;
            }
            if(respcnt != 0) printf("\n");
            break;
          }
          printf("Ampex tape read request has been canceled.\n");
          abortit(1);
        }
        else {
          printf("Failed to get TAPECANCEL msg to ampex_svc\n");
          /*abortit(1);*/
        }
      }
    }
  }
  else {
    printf("\n***PE received a termination signal\n");
    abortit(1);
  }
}

/* This is called when an expected response has timed out. Check to see if
 * any of the busy servers have exited (i.e. most likely crashed).
 * Abort if anyone who is suppose to be there is gone.
*/
void who_died()
{
  SERVER *sptr;
  HDATA *hnext;
  struct taskinfo *taskp;
  int i, j, ntasks;

  for(i=0; i<MAX_SERV; i++) {
    sptr = &stab[i];
    if(sptr->name == NULL)
      break;                            /* everyone is still there */
    if(sptr->busyall) {
      for(hnext=(HDATA *)gethnext(sptr->hosts); hnext != NULL;
         hnext=(HDATA *)gethnext((HDATA *)MONE))
      {
        if(hnext->busy) {
          if(pvm_tasks(pvm_tidtohost(hnext->tid), &ntasks, &taskp)) {
            pemail("***Error on pvm_tasks() call\n");
            abortit(1);
          }
          for(j=0; j<ntasks; j++) {
            if(taskp[j].ti_tid == hnext->tid) break;
          }
          if(j == ntasks) {
            pelog("%x\tn/a\tcrash\t0\t%s\t<NONE>\t<NONE>\n", hnext->tid,sptr->name);
            pemail("***Unexpected exit of %s tid=%x on %s **\n",
                        sptr->name, hnext->tid, hnext->host_name);
            /* I don't think the below is correct. pe doens't know how to
             * proceed if unexpected exit. Just wrap up pe for now.
            */
            if (sptr->noabort) {     /*  map file had ABORT_ACTION=continue  */
/*
              pemail ("    Not aborting as map file had ABORT_ACTION=continue override\n");
*/
              hnext->tid = 0;           /* don't try to kill it if abort*/
              hnext->busy = 0;          /* indicate this host is free */
              sptr->busyall--;          /* decrement master busy count */
              return;
            }
            hnext->tid = 0;             /* don't try to kill it when abort*/
            abortit(1);
          }
        }
      }
    }
  }
}


/*
 *  arglist.c						  R S Bogart  93.04.01
 *
 *  This function loops through the parameter keylist passed to it,
 *    checking for the first parameter keyname to match the target
 *    (arguments) name.  If no match is found, a 0 is returned.  The
 *    calling (main) program is responsible for obtaining a matching
 *    parameter.
 *  If a match is found, the keytype is checked against the desired
 *    (arguments) type.  Note that there is not a one-to-one mapping
 *    between these types.  If a match is found, the program returns
 *    1.  Otherwise, it returns 0.  The calling program is responsible
 *    for obtaining additional parameters or fixing the type.
 *  Mapping of argument type to keytype:
 *      argument type	matches
 *	ARG_DATASET	KEYTYP_STRING (special)
 *	ARG_FLAG	KEYTYP_BYTE, KEYTYP_SHORT, KEYTYP_INT, KEYTYP_LONG
 *	ARG_TIME	KEYTYP_TIME
 *	ARG_INT		KEYTYP_BYTE, KEYTYP_SHORT, KEYTYP_INT, KEYTYP_LONG
 *			KEYTYP_UBYTE, KEYTYP_USHORT, KEYTYP_UINT, KEYTYP_ULONG
 *	ARG_FLOAT	KEYTYP_FLOAT, KEYTYP_DOUBLE
 *	ARG_STRING	KEYTYP_STRING
 *	ARG_FILEPTR	KEYTYP_FILEP, KEYTYP_STRING (special)
 *
 *  The special case ARG_FILEPTR matches KEYTYP_STRING when the value
 *    is either "<" or ">".  The special case ARG_DATASET requires a
 *    string that is a Dataset Name resolved to a file path name plus
 *    template for selection and slicing.
 *
 *  In general, the parameter list will only contain types of either
 *    BYTE (for flags) or STRING.  Thus, these are really the only
 *    types that need to be checked.
 */
//#include <module.h>

int arg_available (KEY *param, char *name, int type)
{
int paramtype = getkeytype (param, name);
			       /*  For strings, key->type is string length  */
if (paramtype == KEYTYP_BYTE)
{
  if (type == ARG_FLAG)
    return (1);
  else
    return (0);
}
else if (paramtype >= KEYTYP_STRING)
  return (1);
else
  return (0);
}

/* Print the usage message and abort
*/
void usage()
{
  printk("Usage:\njpe [-v] [-A] [touch=#] map=map_file_name\n");
  printk("where: -v = verbose\n");
  printk("       -A = access the tape_svc if need to retrieve data\n");
  printk("       touch = #of days to retain input datasets before deletable\n");
  abortit(1);
}

/* A password is required for local mode. Get it or abort.
*/
void getpasswd()
{
  int getpasswd;
  char *envpasswd;
  char pline[MAX_STR];

  getpasswd = 1;
  if(envpasswd = (char *)getenv("LMODEPASSWD")) { /* try env vrbl first */
    if(!strcmp(envpasswd, "mdi4soi")) getpasswd = 0;
  }
  if(getpasswd) {
    printf("Need password to run pe in local mode: passwd = ");
    if(system("stty -echo")) {
      printf("\n  You're not the owner of this tty (you did an su?)\n");
      printf("  WARNING: Your passwd will be echoed\n");
    }
    if(gets(pline) == NULL) {
      system("stty echo");
      printf("\n***Invalid passwd\n");
      abortit(1);
    }
    system("stty echo");
    if(strcmp(pline, "mdi4soi")) {
      printf("\n***Invalid passwd\n");
      abortit(1);
    }
    printf("\n");
  }
}

/* Open the pe mail file if this is a production user. Insert the initial
 * id string.
*/
void open_pemail(char *filename, char *idstr)
{
  int i;
  char string[128];

  mailtrue = 0;         /* set if ever call pemail() */
  mailable = 0;         /* set if this user generates mail */
  for(i = 0; i < mailable_users_num; i++) {
    if(!strcmp(username, mailable_users[i])) mailable = 1;
  }
  if(mailable) {
    if((pemailfp=fopen(filename, "w")) == NULL) {
      fprintf(stderr, "Can't open the mail file %s\n", filename);
      mailable = 0;
      return;
    }
    fprintf(pemailfp, idstr);
  }
}

/* Outputs the variable format message (re: printf) to the pe mail file
 * if enabled and also does a msg() call for stdout (or tae). Sets a flag
 * so we know if ever called and therefore mailit() should mail.
*/
int pemail(char *fmt, ...)
{
  va_list args;
  char string[32768];

  va_start(args, fmt);
  vsprintf(string, fmt, args);
  printk(string);
  if(mailable) {
    fprintf(pemailfp, string);
    fflush(pemailfp);
  }
  va_end(args);
  mailtrue = 1;
  return(0);
}

/* If the pemail is enabled then send the mail file to the mailees.
*/
void mailit()
{
  int i;
  char cmd[128];

  if(mailable && mailtrue) {
    /* all mailees get daemon & production mail */
    if(!(strcmp("daemon", username)) || !(strcmp("production", username))) {
      for(i = 0; i < mailees_num; i++) {
        sprintf(cmd, "Mail -s \"pe 3* mail\" %s < %s", mailees[i], mailname);
        system(cmd);
      }
    }
    else {                      /* each user gets their own */
      sprintf(cmd, "Mail -s \"pe 3* mail\" %s < %s", username, mailname);
      system(cmd);
    }
  }
}

/* Called when get a shell escape (!) in the map file, with the cmd line after
 * the ! character. Will intercept a setenv and perform it, else will do a
 * system call with the given cmd line.
*/
void do_cmd(char *line)
{
  char cmd[256];
  char *token, *name, *value;
  char *putenvcmd;

  strcpy(cmd, line);
  if(token=(char *)strtok(line, " \t\n")) {
    if(!strcmp(token, "setenv")) {
      if(name=(char *)strtok(NULL, " \t\n")) {
        if(value=(char *)strtok(NULL, " \t\n")) {
          /*setenv(name, value, 1); */
          putenvcmd = (char *)malloc(128);
          sprintf(putenvcmd, "%s=%s", name, value);
          putenv(putenvcmd);
        }
      }
    }
    else {
      system(cmd);
    }
  }
}

/* Add the server name and message id number to the next 
 * free entry in the server definition table.
 * The server name has "_svc" appended, the message id
 * is the next sequential message id from the global variable msg_id.
*/
void setstab(char *sname)
{
  int i;
  char *s;

  for(i=0; i<MAX_SERV; i++) {
    if(stab[i].name == NULL) {
      s = (char *)malloc((strlen(sname))+5);
      strcpy(s, sname);
      strcat(s, "_svc");
      stab[i].name=s;
      msg_id++;
      stab[i].msgid=msg_id;
      /* keep a seperate ver# for ea server in case we implement this later */
      strcpy(stab[i].version, soi_version);	/* from soi_version.h */
      break;
    }
  }
  if(i == MAX_SERV) {
    pemail("***Map file exceeds pe limit of %d servers\n",MAX_SERV);
    abortit(1);
  }
}

/* Determine if a pvm daemon is running on the given host. Start one if not.
 * Returns non-0 on error.
*/
int addpvmd(char *host)
{
  struct hostinfo *hostp;
  char *hostptr[2];
  int nhost, narch, i, hinfos;

  /* get info on the present virtual machine configuration */
  if(pvm_config(&nhost, &narch, &hostp)) {
    pemail("***Can't get pvm configuration. Pvm daemon(s) running?\n");
    return(1);
  }
  for(i=0; i<nhost; i++) {
    if(!strcmp(hostp[i].hi_name, host)) {
      break;
    }
  }
  if(i == nhost) {
    printk("A pvm daemon is been added for host %s...\n", host);
    /* you must let any previous spawns settle down before you do the */
    /* pvm_addhosts(). If you don't pvm 3.0 gets into a terrible state! */
    sleep(2);
    hostptr[0] = host;
    if((pvm_addhosts(hostptr, 1, &hinfos)) < 0) {
      pemail("***Can't add a pvm daemon for host %s\n", host);
      return(1);
    }
    if(hinfos < 0) {
      pemail("***Can't add %s pvm daemon. Error %d.\n", host, hinfos);
      pemail("An independently started one may be running? Halt it.\n");
      return(1);
    }
    sleep(2);
  }
  return(0);
}


/* Sets up all the server definition table to have all the given hosts for each
 * server. The server names must already be set in stab[]. The input hname is a
 * list of comma delimited host names (e.g. rumble,quake,dragon). 
 * Also copies the servers' map file keylist to each host's keylist.
 * Sets the global variable num_hosts and effective_hosts. 
 * Will start a pvm daemon on each host if one is not already there.
 * Returns non-0 on error.
*/
int sethosts(char *hname)
{
  HDATA *hnext;
  char *host;
  int j;

  num_hosts=0;
  host=(char *)strtok(hname,",");
  while(host) {
    for(j=0; j<MAX_SERV; j++) {         /* set hosts for each server */
      if(stab[j].name == NULL)
        break;                          /* exit for(j) loop */
      sethname(&stab[j].hosts, host);   /* add host name to hdata list */
    }
    if(addpvmd(host)) return(1);	/* add a pvm daemon if needed */
    num_hosts++;
    host=(char *)strtok(NULL,",");
  }
  effective_hosts = num_hosts;
  if(effective_hosts != 1) {
    pemail("**Multiple hosts no longer allowed. Will use first host only\n");
    effective_hosts = 1;
  }

  for(j=0; j<MAX_SERV; j++) {		/* copy all map_lists to all hosts */
    if(stab[j].name==NULL)
      break;
    for(hnext=(HDATA *)gethnext(stab[j].hosts); hnext != NULL;
	 hnext=(HDATA *)gethnext((HDATA *)MONE)) {
      hnext->param_list = newkeylist();
      add_keys(stab[j].map_list, &hnext->param_list);
    }
  }
  return(0);
}


/*
 *  cmdparam.c					~soi/(version)/src/libast.d
 *
 *  This is a quick and dirty piece of the parameter evaluation
 *    stuff from the original param.c.  It is put in to use only what
 *    is needed and avoid inclusion and linking with other libraries
 *  Functions:
 *  CollectParams	creates key-list of all command line parameters
 *
 *  Responsible:  Rick Bogart			RBogart@solar.Stanford.EDU
 *
 *  Bugs:
 *    No man pages.
 *    CollectParams() ignores -Pxxx printer flag.
 */
int CollectParams (KEY **newlist, int argc, char **argv)
{
int pct, i;
char *arg, *valu, *keyname;
char isval, flagval;
static char flagname[] = "x"; /* a string with length 1 */

pct = argc;
*newlist = newkeylist ();

isval = 0;
for (i=1; i<argc; i++)
{
  arg = argv[i];
  if (isval)
		 /*  If previously flagged as value, fill remainder of key  */
  {
    setkey_str (newlist, keyname, arg);
    free (keyname);
    isval = 0;
    pct--;
    continue;
  }
  if (arg[0] == '-' || arg[0] == '+')
	   /*  Set as byte type for flags, starting with - or + characters  */
  {
    valu = arg;
    valu++;
    isval = 0;
    if (arg[0] == '-' && arg[1] == 'P')
    {
						/*  Ignore -P printer flag  */
      pct--;
      continue;
    }
    flagval = (arg[0] == '-') ? -1 : 1;
    while (*valu)
    {
      flagname[0] = valu[0];
      setkey_byte (newlist, flagname, flagval);
      valu++;
    }
    continue;
  }
  if (valu = index (arg, '='))
  {
    *valu = 0;
    keyname = key_strdup (arg);
    *valu = '=';
    valu++;
    if (strlen (valu))
    {
      setkey_str (newlist, keyname, valu);
      free (keyname);
      isval = 0;
    }
    else
    {
      isval = 1;
    }
    continue;
  }
						       /*  Default: ignore  */
  pct--;
  isval = 0;
}
return (pct - 1);
}

void get_cmd(int argc, char *argv[])
{
  int xargc, first, pfound, groupflg, groupid, gid, archgrp;
  int qon, i;
  char *toks[128];
  char *tokens1 = "*\t\n";
  char *tokens2 = " \t\n";
  char *tokens3 = "*\n";
  char *tokens4 = " \n";
  char hostname[MAX_STR], logkill[MAX_STR];
  char *inname, *sp, *hname, *cptr;
  char *line, *line2, *token, *tokenchars, *tokencharsend;
  char carch;
  FILE *mapfp;
  SERVER *sptr;
  KEY *work_list;                       /* keylist from ea line in map file */

  CollectParams(&work_list, argc, argv);/* get calling sequence */
  if(!(inname=getkey_str(work_list, "map"))) {
    usage();                            /* abort */
  }
  gethostname(hostname, MAX_STR);
  if(!findkey(work_list, "touch"))     /* see how long retain any input ds*/
    touchflg = -1;                      /* don't change current setting */
  else {
    cptr = getkey_str(work_list, "touch");
    if(!strcmp(cptr, "keep")) {   /* it's a touch "keep" */
      touchflg = 9999999;         /* pass on the default keep value */
    }
    else {
      if(!isdigit((int)*cptr)) {
        pemail("***touch= must say \"keep\" or give integer #of days\n");
        abortit(1);
      }
      touchflg = atoi(cptr);
    }
  }
  if(findkey(work_list, "pds")) {
    strcpy(pdshost, getkey_str(work_list, "pds")); /* host for /PDS storage*/
  } else if((char *)getenv("PDS_SET_HOST")) {
    strcpy(pdshost, (char *)getenv("PDS_SET_HOST"));
  } else {
    strcpy(pdshost, hostname);          /* use this host */
  }
  if(cptr = index(pdshost, '.')) {
    *cptr = NULL;                       /* remove any .Stanford.EDU */
  }
  if(findkey(work_list, "rcp")) {
    rcpdir=getkey_str(work_list, "rcp");
    if(strncmp(rcpdir, "/", 1)) {
      pemail("***The rcp= dir must start with \"/\"\n");
      abortit(1);
    }
    rcpflg = 1;
  }
  /*taeflg=(getkeytype(work_list, "t") == KEYTYP_BYTE);*/
  verbose=(getkeytype(work_list, "v") == KEYTYP_BYTE);
  //debugflg=(getkeytype(work_list, "d") == KEYTYP_BYTE);
  ampexflg=(getkeytype(work_list, "A") == KEYTYP_BYTE);
  if(getkeytype(work_list, "O") == KEYTYP_BYTE) {
    printf("\npe does not have a \"-O\" switch.\n");
    abortit(1);
  }
  if(!(mapfp=fopen(inname, "r"))) {
    pemail("***Can't open the map file %s\n", inname);
    abortit(1);
  }
  strcpy(pe_map, inname);			/* save as global */
  hname=hostname;				/* default HOST= in map file */
  sptr = &stab[0];				/* first server table */
  first=1; pfound=0; groupflg=0; gid=0; groupid=0, archgrp=0;
  line = (char *)malloc(16384);
  line2 = (char *)malloc(16384);
  while(fgets(line, 16384, mapfp)) {		/* get map file lines */
    if(line[0] == '#' || line[0] == '\n') continue;
    if(line[0] == '!') {
      do_cmd(line+1); 
      continue;
    }
    if(strchr(line, '"')) {	/* parse for quoted string */
      tokenchars = tokens1;
      tokencharsend = tokens3;
      /* substitue "*" for space except in quoted strings */
      /* Note: can't use "*" in quoted strings */
      qon = 0; line2[0] = NULL;
      for(i=0; ; i++) {
        if(line[i] == '"') {
          if(qon) qon = 0;
          else qon = 1;
          i++;			/* dont include the quote */
        }
        if((line[i] == '*') && qon) {
          pemail("Illegal \"*\" char in quoted string\n");
          abortit(1);
        }
        if((line[i] == ' ') && !qon) {
          line[i] = '*';
        }
        strncat(line2, &line[i], 1);
        if(line[i] == '\n') break;
      }
      if(qon) {
        pemail("Mismatched quotes in: %s\n", line);
        abortit(1);
      }
    } else {			/* no " in line */
      tokenchars = tokens2;
      tokencharsend = tokens4;
      strcpy(line2, line);
    }
    toks[0]=argv[0];
    xargc=1;
    token=(char *)strtok(line2, tokenchars);
    while(token) {
      if(token[0] == '\\') {
          line = (char *)malloc(16384);
          line2 = (char *)malloc(16384);
          while(fgets(line, 16384, mapfp)) {
            if(line[0] == '#' || line[0] == '\n') continue;
            break;
          }
          if(strchr(line, '"')) {	/* parse for quoted string */
            tokenchars = tokens1;
            tokencharsend = tokens3;
            /* substitue "*" for space except in quoted strings */
            /* Note: can't use "*" in quoted strings */
            qon = 0; line2[0] = NULL;
            for(i=0; ; i++) {
              if(line[i] == '"') {
                if(qon) qon = 0;
                else qon = 1;
                i++;			/* dont include the quote */
              }
              if((line[i] == '*') && qon) {
                pemail("Illegal \"*\" char in quoted string\n");
                abortit(1);
              }
              if((line[i] == ' ') && !qon) {
                line[i] = '*';
              }
              strncat(line2, &line[i], 1);
              if(line[i] == '\n') break;
            }
            if(qon) {
              pemail("Mismatched quotes in: %s\n", line);
              abortit(1);
            }
          } else {			/* no " in line */
            tokenchars = tokens2;
            tokencharsend = tokens4;
            strcpy(line2, line);
          }
          if(!(token=(char *)strtok(line2, tokenchars)))
            break;
      }
      toks[xargc]=token;
      xargc++;
      token=(char *)strtok(NULL, tokencharsend);
    }
    CollectParams(&work_list, xargc, toks);	/* get into a keylist */
    if(sp=getkey_str(work_list, "DSDSOUT"))
      reqmegs=atoi(sp);
    else if(sp=getkey_str(work_list, "NOWARN"))
      nowarn=atoi(sp);
    else if(sp=getkey_str(work_list, "HOST"))
      hname=sp;
    else if(sp=getkey_str(work_list, "START_ARCHIVE")) {
      if(atoi(sp)) {				/* START_ARCHIVE is true */
        if(archgrp) {
          pemail("***Non-matching START/END_ARCHIVE pair\n");
          abortit(1);
        }
        archgrp = 1;				/* must be 1 for below */
      }
    }
    else if(sp=getkey_str(work_list, "END_ARCHIVE")) {
      if(atoi(sp)) {
        if(!archgrp) {
          pemail("***Non-matching START/END_ARCHIVE pair\n");
          abortit(1);
        }
        archgrp = 0; 
        /* now change the archive_group of the last server from 1 to -1 */
        if(pfound) {				/* make sure there was a p= */
          --sptr;				/* go back to prev table */
          sptr->archive_group = -1;		/* indicated end of group */
          ++sptr;
        }
      }
    }
    else if(sp=getkey_str(work_list, "START_GROUP")) {
      if(atoi(sp)) {				/* START_GROUP is true */
        if(groupflg) {
          pemail("***Non-matching START/END_GROUP pair\n");
          abortit(1);
        }
        groupflg = 1;
        groupid = ++gid;
      }
    }
    else if(sp=getkey_str(work_list, "END_GROUP")) {
      if(atoi(sp)) {
        if(!groupflg) {
          pemail("***Non-matching START/END_GROUP pair\n");
          abortit(1);
        }
        groupflg = 0; 
        groupid = 0;
      }
    }
    else if(sp=getkey_str(work_list, "p")) { 	/* line w/process spec */
      pfound=1;
      deletekey(&work_list, "p");		/* remove if from keylist */
      setstab(sp);				/* put name in server table */
      sptr->cphist = 1;				/* default cp hist log */
      sptr->groupid = groupid;			/* real group id or 0 */
      sptr->archive_group = archgrp;		/* 0 or 1 if in archive group */
      if(first) {
        first=0;
        sptr->firstserver=1;
      }
      if(sp=getkey_str(work_list, "d")) {	/* ck for dsds in flag */
        if(sptr->dsin=atoi(sp))
          dsdsin=1;				/* flg if any d=1 */
        deletekey(&work_list, "d");
      }
      //NEW for jpe we must always query DRMS for input so porce flags to 1
      sptr->dsin = 1;
      dsdsin=1;
					   /*  removed 99.02.12 (see below)  */
/*
      if(sp=getkey_str(work_list, "l")) {
        sptr->cphist=atoi(sp);
        deletekey(&work_list, "l");
      }
*/
		    /*  check for copy suppression (new key as of 97.12.16)  */
		       /*  (intent is to remove "l" from special args list)  */
      if (sp = getkey_str (work_list, "COPY_HISTORY")) {
        sptr->cphist = strcasecmp (sp, "no");
        deletekey (&work_list, "COPY_HISTORY");
      }
					   /*  removed 99.02.12 (see below)  */
/*
      if(sp=getkey_str(work_list, "x")) {
        sptr->noabort=atoi(sp);
        deletekey(&work_list, "x");
      }
*/
		      /*  check for abort override (new key as of 97.11.18)  */
		       /*  (intent is to remove "x" from special args list)  */
      if (sp = getkey_str (work_list, "ABORT_ACTION")) {
        sptr->noabort = strcasecmp (sp, "continue") ? 0 : 1;
        deletekey (&work_list, "ABORT_ACTION");
      }
      if(sp=getkey_str(work_list, "a")) {	/* ck for archive flag */
        carch = *sp++;
        if(carch != 'a' && carch != 't' && carch != 'p' && carch != 'n' && carch != '0') {
          pemail("***Illegal a= spec for p=%s. Must be a, t, p or n\n",sptr->name);
          abortit(1);
        }
        if(carch == '0') {
          //sptr->archive = NULL;
          sptr->archive = 't';		// force a=0 to t for jpe
        }
        else {
          sptr->archive = carch;		/* archive char a, t,p or n */
          sptr->archive_day=atoi(sp);		/* will be 0 if fails */
        }
        deletekey(&work_list, "a");		/* remove it from keylist */
      }
      else {
        //sptr->archive = NULL;
        sptr->archive = 't';			// force a=0 to t for jpe
      }
      if(sp=getkey_str(work_list, "s")) {	/* ck for split fsn-lsn flag */
        sptr->split=atoi(sp);
        sptr->split=0;				/* !! implement split later */
        deletekey(&work_list, "s");		/* remove it from keylist */
      }
      add_keys(work_list, &sptr->map_list);	/* add to server keylist */
      sptr++;					/* next server table */
    }						/* end if("p") */
    else {
      pemail("***A map file line must start w/a control stmt, or p=\n");
      pemail("  %s\n", line);
      abortit(1);
    }
  }						/* end while(fgets()) */
  fclose(mapfp);
  if(!hname) {
    pemail("***The map file has no HOST specification\n");
    abortit(1);
  }
  if(!pfound) {
    pemail("***The map file has no p= line to define the server\n");
    abortit(1);
  }
  if(groupflg) {
    pemail("***Non-matching START/END_GROUP pair\n");
    abortit(1);
  }
  if(archgrp) {
    pemail("***Non-matching START/END_ARCHIVE pair\n");
    abortit(1);
  }
  if(sethosts(hname)) abortit(1);	/* hosts & map_lists to serv table */
  /* elim pekill log */
  /*sprintf(logkill, "/usr/local/logs/pekill/pekill.%s.log", username);
  /*open_pekill_log(logkill);
  */
}

/* Determine the pvm machine configuration and which tasks are running.
 * Spawn the pe_rpc that will bridge us to the dsds_svc via the pe_rpc_svc.
 * Spawn the servers as required at each host.
*/
void spawn_pvm()
{
  struct hostinfo *hostp;
  SERVER *sptr;
  HDATA *hnext;
  int i, nhost, narch;

  /* get info on the present virtual machine configuration */
  if(pvm_config(&nhost, &narch, &hostp)) {
    pemail("***Can't get pvm configuration. Pvm daemon(s) running?\n");
    abortit(1);
  }
  /* assume that all servers run on all hosts given in the cmd line */
  for(hnext=(HDATA *)gethnext(stab[0].hosts); hnext != NULL;
         hnext=(HDATA *)gethnext((HDATA *)MONE))
  {
    for(i=0; i<nhost; i++) {
      if(!strcmp(hostp[i].hi_name, hnext->host_name))
      /*if(!strstr(hnext->host_name, hostp[i].hi_name))*/
        break;
    }
    if(i == nhost) {
      pemail("***No pvm daemon running on %s\n", hnext->host_name);
      abortit(1);
    }
  }

  /* spawn the servers on each host (host names already set in stab[].hosts) */
  /* As of Apr 11, 1997 this is now done in kick_server(). */
/*  for(i=0; i<MAX_SERV; i++) {
/*    sptr = &stab[i];
/*    if(sptr->name == NULL)
/*      break;					/* no more servers */
/*    hnext = (HDATA *)gethnext(sptr->hosts);	/* get first hdata entry */
/*    while(hnext) {
/*      if(debugflg)
/*        pvm_spawn(sptr->name, (char **)0, PvmTaskDebug, hnext->host_name, 1,
/*		&hnext->tid);
/*      else
/*        pvm_spawn(sptr->name, (char **)0, PvmTaskHost, hnext->host_name, 1,
/*		&hnext->tid);
/*      if(hnext->tid < 0) {
/*        printf("***Can't spawn %s on %s\n",sptr->name,hnext->host_name);
/*       abortit(1);
/*      }
/*      printf("%s tid=%x spawned on %s\n",
/*		sptr->name, hnext->tid, hnext->host_name);
/*      hnext = (HDATA *)gethnext((HDATA *)MONE);	/* get next entry */
/*    }
/*  }
*/
}

/* Send the given server for the given host the message id that that server
 * type will use for completion messages to PE. This must be the first thing
 * sent to servers. Also sends the verbose flag and any history and errlog
 * files specified in the pe map file.
*/
void msgid_send(SERVER *sptr, HDATA *hnext)
{
  KEY *alist, *keybad;
  char *log;

  pvm_initsend(PvmDataDefault);
  pvm_pkint(&sptr->msgid, 1, 1);        /* pack the msg id */
  /*pvm_pkint(&verbose, 1, 1);          /* and the verbose flag */
  /* make keylist w/verbose flg and any history and errlog from map file*/
  alist=newkeylist();
  setkey_int(&alist, "verbose", verbose);
  if(log=getkey_str(sptr->map_list, "history"))
    setkey_str(&alist, "history", log);
  if(log=getkey_str(sptr->map_list, "errlog"))
    setkey_str(&alist, "errlog", log);
  if(keybad=(KEY *)pack_keylist(alist)) {
    pemail("***Err packing pvm msg, type=%d name=%s\n",keybad->type,keybad->name
);
    abortit(1);
  }
  if(pvm_send(hnext->tid, MSGMSG)) {
    pemail("***Error sending %s on %s its message id\n",
        sptr->name,hnext->host_name);
    abortit(1);
  }
}

/* Ask for and get the argument list from the given server.
*/
void arg_recv(SERVER *sptr, HDATA *hnext)
{
  argument *args;
  struct timeval tvalr;
  uint64_t tsr;
  int bufid, i, j;

    pvm_initsend(PvmDataDefault);
    if(pvm_send(hnext->tid, MSGARGS)) {
      pemail("***Error calling %s on %s for arg list\n",sptr->name,hnext->host_name);
      abortit(1);
    }
    gettimeofday(&tvalr, NULL); tsr = tvalr.tv_sec;
    while(1) {
      if(!(bufid=pvm_nrecv(hnext->tid, MSGARGS))) { /* haven't got msg yet */
        gettimeofday(&tvalr, NULL);
        if(tvalr.tv_sec - tsr > RESPWAIT) {
          pemail("***Timeout awaiting an argument list from %s\n",sptr->name);
          abortit(1);
        }
      }
      else
        break;
    }
    if(bufid < 0) {
      pemail("***Error receiving arg list from %s\n",sptr->name);
      abortit(1);
    }
    for(j=0; j<MAX_ARGS; j++) {
      args = &sptr->arguments[j];
      pvm_upkint(&args->kind,1,1);
      args->key=(char *)malloc(MAX_ARG_STR);
      pvm_upkstr(args->key);
                 /*  Check if module is using an arg that is reserved by pe  */
      if (!strcmp (args->key, "p") || !strcmp (args->key, "a")) {
        printf("**WARNING: Module has arg %s which is reserved by pe\n",
			args->key);
      }
      args->default_value=(char *)malloc(MAX_ARG_STR);
      pvm_upkstr(args->default_value);
      args->range=(char *)malloc(MAX_ARG_STR);
      pvm_upkstr(args->range);
      args->description=(char *)malloc(MAX_ARG_STR*4);
      pvm_upkstr(args->description);
      if(args->kind == ARG_END)
        break;                  /* got all args, go onto next server */
    }                           /* end for(MAX_ARGS) */
    if(j == MAX_ARGS) {
      pemail("***Never got an ARG_END from %s after %d arguments\n",
      sptr->name, MAX_ARGS);
      abortit(1);
    }
}

/* Check that the given output working directory is not a duplicate of
 * a previous output wd. Called with the wd to check and the keylist to
 * check it in (normally wd_dup_list or wd_mk_list which has all the wd's
 * so far). Returns non-0 if this is a dup wd.
*/
int wd_dup_ck(KEY *list, char *wd)
{
  KEY *walker = list;

  while(walker) {
    if(!strcmp(wd, (char *)walker->val)) {
      /* elim old abort on dups
      printf("***Error: pe detects two datasets that have the same output wd:\n");
      printf(" wd = %s\n", wd);
      printf(" You do not have your prog name template set to give unique wd's\n");
      abortit(1);
      */
      return(1);
    }
    walker = walker->next;
  }
  return(0);
}


/* Satisfy and expand the arg list for the given server and host entry.
 * Stores the resulting key list in the server's param_list.
 * Makes any output directories if output is in dsds storage.
 * Will abort if an arg list cannot be satisfied. 
*/
void ck_arglist(SERVER *sptr, HDATA *hx)
{
  DIR *dfd;
  double dblval;
  int found, outnsets, i, intval, parse_status;
  char dirstr[MAX_STR], inname[MAX_STR], ext[MAX_STR], wdkey[MAX_STR];
  char *rootdbase, *wd, *valstr, *failstr, *strck, *cptr;
  argument *arg;

  arg=sptr->arguments;
  while(arg->kind != ARG_END) {
    found=arg_available(hx->param_list, arg->key, arg->kind);
    if (!found) {	/*  Then use default from module argument table  */
      if (found = strlen (arg->default_value));
        setkey_str (&hx->param_list, arg->key, arg->default_value);
    }
    if(!found) {
      pemail("***Can't satisfy arg \"%s\" for server %s\n", arg->key,sptr->name);
      abortit(1);
    }
    if(arg->kind == 1) {	/* catch any obsolete ARG_DATASET use */
      pemail("***An obsolete ARG_DATASET found in %s\n", sptr->name);
      abortit(1);
    }
    if (arg->kind == ARG_DATA_OUT) { 
//!!TBD ck about a=a for jpe
      if(sptr->archive == 'a') { /* do mkdir case for append archive */
        strcpy(inname, arg->key); strcat(inname, "_nsets");
        outnsets = getkey_int(hx->param_list, inname);
        for(i=0; i < outnsets; i++) {   /* do for each dataset */
          sprintf(ext, "%s_%d_wd", arg->key, i);
          wd = getkey_str(hx->param_list, ext);
          if((dfd=opendir(wd)) == NULL) { /* dir doesn't exist, make it */
            strcpy(dirstr, "mkdir -p "); strcat(dirstr, wd);
            if(system(dirstr)) {
              pemail("***Cannot %s for %s\n",dirstr,sptr->name);
              abortit(1);
            }
          }
        }
      }
      else {			/* non append archive case */
        if(reqmegs) {		/* ouput is in dsds storage */
          strcpy(inname, arg->key); strcat(inname, "_nsets");
          outnsets = getkey_int(hx->param_list, inname);
          /* put any {dbase} dir in the duplicate wd list */
          sprintf(ext, "%s_0_dbase", arg->key);
          if(findkey(hx->param_list, ext)) {
            rootdbase = getkey_str(hx->param_list, ext);
            sprintf(wdkey, "wd_%d", wdkey_int++); /* form uniq name of key */
            setkey_str(&wd_mk_list, wdkey, rootdbase);
          }
          for(i=0; i < outnsets; i++) {	/* do for each dataset */
            sprintf(ext, "%s_%d_wd", arg->key, i);
            cptr = (char *)sptr->archive;
            if(cptr == (char *)'t' || cptr == (char *)'p' || cptr == (char *)'n') {
              wd = getkey_str(hx->param_list, ext);
              cptr = strstr(wd, "/SUM");
              if(cptr != wd) {		/* not a /SUM wd illegal if archive */
                pemail("Cannot archive non-/SUM storage: %s = %s\n", ext,wd);
                abortit(1);
              }
            }
            /* don't make same dir twice to avoid error msg */
            if(!wd_dup_ck(wd_mk_list, getkey_str(hx->param_list,ext))) {
//!!!TBD ck mkdir here
              strcpy(dirstr, "mkdir -p ");
              strcat(dirstr, getkey_str(hx->param_list, ext));
              if(system(dirstr)) {
                pemail("***Cannot %s for %s\n",dirstr,sptr->name);
                abortit(1);
              }
              sprintf(wdkey, "wd_%d", wdkey_int++); /* form uniq name of key */
              setkey_str(&wd_mk_list, wdkey, getkey_str(hx->param_list, ext));
            }
          }
        }			/* end ouput is in dsds storage */
      }				/* end non append archive case */
    }
    else if (arg->kind == ARG_FILEPTR) {
      pemail("***ARG_FILEPTR used by %s, is not supported by pe\n",sptr->name);
      abortit(1);
    }
    else if (arg->kind == ARG_FLOAT) {
       valstr = getkey_str (hx->param_list, arg->key);
       dblval = strtod(valstr, &failstr);
       if (valstr == failstr) dblval = D_MISSING;
       setkey_double (&hx->param_list, arg->key, dblval);
    }
    else if (arg->kind == ARG_INT) {
       valstr = getkey_str (hx->param_list, arg->key);
       intval = (int)strtol(valstr, &failstr,0);
       if (valstr == failstr) intval = I_MISSING;
       setkey_int (&hx->param_list, arg->key, intval);
    }
    else if (arg->kind == ARG_TIME) 
       setkey_time (&hx->param_list, arg->key,
        sscan_time (getkey_str (hx->param_list, arg->key)));
    else if (arg->kind == ARG_FLAG) {
       if(!(strck = getkey_str (hx->param_list, arg->key))) {
         /*pemail("Illegal ARG_FLAG %s\n", arg->key);*/
         /*abortit(1);*/
         setkey_byte(&hx->param_list, arg->key, 1); /* force to 1 */
       }
       else {
         setkey_byte(&hx->param_list, arg->key, atoi(strck));
       }
    }
    else if (arg->kind == ARG_FLOATS) {
      /*keyiterate(printkey, hx->param_list); /* !!TEMP */
      if(parse_status=parse_array (&hx->param_list, arg->key, KEYTYP_DOUBLE)) {
        pemail("***ARG_FLOATS error %d for arg %s\n", parse_status, arg->key);
        abortit(1);
      }
    } 
    else if (arg->kind == ARG_INTS) {
      if(parse_status=parse_array (&hx->param_list, arg->key, KEYTYP_INT)) {
        pemail("***ARG_INTS error %d for arg %s\n", parse_status, arg->key);
        abortit(1);
      }
    } else if (arg->kind == ARG_NUME) {
      char **names;
      int nval, nvals = parse_numerated (arg->range, &names);
      valstr = getkey_str (hx->param_list, arg->key);
      for (nval = 0; nval < nvals; nval++)
        if (!strcmp (names[nval], valstr)) break;
      if (nval >= nvals) nval = 0;
      setkey_int (&hx->param_list, arg->key, nval);
    }
    arg++;
  }						/* end while(!=ARG_END) */
}

/* Called by form_keylist() for the first input dataset specified in a map
 * file for a server. Modifies the fsn & lsn for the number of hosts if the 
 * split flag was given in the map for this server. Will not do any fsn-lsn
 * split if the dataset is a datacollection or no fsn-lsn in the keylist.
 * Puts the appropriately updated keylist into the param_list for each host.
 * Sets the global vrb effective_hosts for the number of hosts to be actually
 * used.
 * !!NOTICE: The split flag will always be 0 for now. Full implementation
 * will be done later when it is clear how we want to use it. The split
 * related code that you'll see is only a partial implementation, mostly wrt
 * the input dataset. Must figure out what to do with it in process_resp().
*/
void form_split(SERVER *sptr, KEY *xlist, char *key)
{
  HDATA *hnext;
  char keystr[MAX_STR];
  int i, fsn, lsn, delta, drem, xfsn, xlsn, nosplit;

/* !! TEMP */
/*  inname=getkey_str(xlist, key);		/* get the ds name */
/*  printf("***** At form_split() the keylist for %s is:\n", inname);
/*  keyiterate(printkey, xlist);
/* !! END TEMP */

  nosplit = 0;
  if(!sptr->split) {
    nosplit = 1;				/* no split requested */
  }
  else {
    sprintf(keystr, "%s_nsets", key);
    if(findkey(xlist, keystr)) {
      if((getkey_int(xlist, keystr)) != 1) {
        printf("**Ignoring the s=1 split flag given for %s:\n", sptr->name);
        printf("  The %s= specification is not a single dataset\n", key);
        nosplit = 1;				/* this is a datacollection */
      }
      else {
        sprintf(keystr, "%s_0_fsn", key);
        if(!findkey(xlist, keystr)) {	/* no "in"_0_fsn */
          printf("**Ignoring the s=1 split flag given for %s:\n", sptr->name);
          printf("  No sel:[fsn-lsn] for first input ds given in map file or found from dsds\n");
          nosplit = 1;
        }
      }
    }
    else
      nosplit = 1;				/* this is not a dataset */
  }
  if(nosplit) {
    effective_hosts = 1;
    hnext = (HDATA *)gethnext(sptr->hosts);	/* get first hdata entry */
    add_keys(xlist, &hnext->param_list);	/* copy over all keys */
    return;
  }

  sprintf(keystr, "%s_0_fsn", key);
  fsn=getkey_int(xlist, keystr);
  sprintf(keystr, "%s_0_lsn", key);
  lsn=getkey_int(xlist, keystr);

  effective_hosts = num_hosts;			/* num_hosts set by get_cmd*/
  delta=((lsn-fsn)+1)/effective_hosts;
  if(delta == 0)
    delta = 1;					/* fewer files than hosts */
  drem=((lsn-fsn)+1) % effective_hosts;
  hnext = (HDATA *)gethnext(sptr->hosts);	/* get first hdata entry */

  for(i=0; i<effective_hosts; i++) {
    if(((lsn-fsn)+1) == i) {			/* fewer files than hosts */
      effective_hosts = i;
      break;
    }
    add_keys(xlist, &hnext->param_list);	/* copy over all keys */
    xfsn=fsn+(delta*i);
    xlsn=xfsn+delta-1;
    if(i+1 == effective_hosts)
      xlsn=xlsn+drem;
    sprintf(keystr, "%s_0_fsn", key);
    setkey_int(&hnext->param_list, keystr, xfsn);
    sprintf(keystr, "%s_0_lsn", key);
    setkey_int(&hnext->param_list, keystr, xlsn);
    /* !!TEMP */
    /*printf("***** In form_split the keylist for %s is:\n", hnext->host_name);
    /*keyiterate(printkey, hnext->param_list);
    /* !! END TEMP */
    /* for now this hnext will be non-null for effective_hosts */
    hnext = (HDATA *)gethnext((HDATA *)MONE);	/* get next hdata entry */
  }
}

/* Called by form_keylist() for each input datacollection key in the argument 
 * table for the given server (i.e. for each ARG_DATA_IN). 
 * Returns the parsed keylist, added to the input keylist. 
 *   xlist = input keylist to add this datacollection to
 *   sptr = SERVER table of this current server (i.e. module)
 *   arg = argument list for this server
 *   seq = 0 for the first input datacollection, then 1, 2, etc.
*/
KEY *form_arg_data_in(KEY *xlist, SERVER *sptr, argument *arg, int seq)
{
  //KEY *blist;
  char dbasekey[MAX_STR], inname[MAX_STR], ext[MAX_STR];
  char *wd;
  int innsets, parse, respcnt, i;

  if(!seq) add_keys(sptr->map_list, &xlist);
  if(!arg_available(xlist, arg->key, arg->kind))  {
    if(!(strlen(arg->default_value))) {	/* no default value */
      pemail("***Can't satisfy arg \"%s\" for server %s\n",arg->key,sptr->name);
      abortit(1);
    }
    else 
      setkey_str(&xlist, arg->key, arg->default_value);
  }
  if(touchflg != -1)
    setkey_int(&xlist, "touch", touchflg);
  sprintf(ext, "arg_data_in_%d", seq);
  setkey_str(&xlist, ext, arg->key);	/* add name of arg */
  /* define {dbase} in case used in the env variable name template */
  strcpy(dbasekey, arg->key); strcat(dbasekey, "_dbase");
  /* setup if input is from previous modules output */
  if(reqmegs && !sptr->firstserver && !sptr->dsin) {
    setkey_str(&xlist, dbasekey, dsdswd);
    sprintf(ext, "%s_rule", arg->key); 	/* e.g. in_rule */
    setkey_str(&xlist, ext, "wd:{dbase}"); /* do in case parse list again */
    sprintf(ext, "%s_%d_rule", arg->key, seq);
    setkey_str(&xlist, ext, "wd:{dbase}"); /* override rule for SUMS compat */
  }
  else if(getenv("dbase"))
    setkey_str(&xlist, dbasekey, getenv("dbase"));
  else
    setkey_str(&xlist, dbasekey, "/tmp"); /* there must be some {dbase} */

  if(parse=parse_list(&xlist, arg->key)) {/* expand ds name */
      if(parse != CANNOT_FILL_TEMPLATE) { /* ignore this error for input ds */
        pemail("***Error %d in parse_list for %s\n", parse, sptr->name);
        if(dbxflg) {
          printf("\n***** The xlist after the parse is:\n");
          keyiterate(printkey, xlist);
        }
        abortit(1);
      }
  }
  if(dbxflg) {
    printf("\n***** The xlist after the parse for %s is:\n", sptr->name);
    keyiterate(printkey, xlist);
  }

  /* if this is an in ds not from dsds then a wd must */
  /* be present in the keylist for each dataset */
  if(!sptr->dsin) {
    /* don't test in the default in is "NOT SPECIFIED" */
    if(strcmp(getkey_str(xlist, arg->key), "NOT SPECIFIED")) {
      strcpy(inname, arg->key); strcat(inname, "_nsets");
      innsets = getkey_int(xlist, inname);
      for(i=0; i < innsets; i++) {          /* do for each dataset */
        sprintf(ext, "%s_%d_wd", arg->key, i);
        if(!(wd=getkey_str(xlist, ext)) || !strcmp(wd, "")) {
          pemail("***No %s for %s. Missing rule or template in map file.\n",
  	 ext, sptr->name);
          abortit(1);
        }
        if(*wd == '.') {
          pemail("***%s of \".\" not valid for %s.\n", ext, sptr->name);
          abortit(1);
          /*printf("**!!WARN: %s of \".\" not valid for %s.\n",ext, sptr->name);*/
        }
      }
    }
  }
  return(xlist);
}

/* Get all the input datasets from dsds_svc. Called with the expanded
 * keylist from form_keylist(). The args are given by arg_data_in_0,
 * arg_data_in_1, etc. keys. Returns a new keylist with the values supplied by
 * dsds_svc. Only returns when all dataset are brought online if -A (ampex)
 * mode. If not -A mode then the wd's are made null here for all ds with a 
 * non-zero status. This is what the module sees for a wd.
*/
KEY *dsds_arg_data_in(KEY *xlist) 
{
  static KEY *blist;
  char inname[MAX_STR], ext[MAX_STR];
  char *wd, *argname, *svcname, *svcversion;
  int respcnt, status, i;

  blist = newkeylist();
  /* now query dsds_svc for each dataset in the input keylist */
  /* the args are given by arg_data_in_0, arg_data_in_1, etc */
  setkey_uint(&xlist, "dsds_uid", uid);  /* give our OPEN uid */
  setkey_str(&xlist, "pds_host", pdshost); /* host for the /PDS set */
  /* send either "ampex_tid" or "lago_tid" keyword to dsds_svc */
  if(ampexflg)
    setkey_int(&xlist, "ampex_tid", ampex_tid);
  else
    setkey_int(&xlist, "lago_tid", lago_tid);
  /* only send svc_version & svc_name if it was in the map file */
  if(svcversion=getkey_str(xlist,"svc_version")) {
    /*setkey_str(&xlist, "svc_version", svcversion);*/
    printf("**Obsolete use of svc_version in map file. Ignoring...\n");
  }
  if(svcname=getkey_str(xlist,"svc_name")) {
    /*setkey_str(&xlist, "svc_name", svcname);*/
    printf("**Obsolete use of svc_name in map file. Ignoring...\n");
  }
  printf("Querying for the input datasets...\n");
  nocontrolc = 1;			/* no ^C during retrieve */

  blist = (KEY *)call_drms_in(xlist, dbxflg);
  if(blist == NULL) {
    pemail("Can't resolve/retrieve the input datasets with drms.\n");
    pemail("(Usually no -A switch specified.)\n");
    abortit(1);
  }
//printf("\n*** The blist after the call_drms_in() in dsds_arg_data_in() is:\n");
//keyiterate(printkey, blist); //!!!TEMP

  /* Need to wait for final reply in case some data was not online */
//  if(findkey(blist, "DS_ARCHIVE")) {	/* data not on-line */
//    printf("Data are being retrieved from tape for one or more input datasets.\n");
//    printf("Waiting (may be long) for data to be retrieved ...\n");
//    respcnt = 0;
//    while(1) {
//      if(resp_dsds(dsds_tid)) {		/* rets every timeout interval */
//        printf("..."); respcnt++;
//        if(respcnt == 20) {
//          respcnt = 0;
//          printf("\n");
//        }
//        continue;
//      }
//      if(respcnt != 0) printf("\n");
//      break;
//    }
//    freekeylist(&blist);
//    if((blist = (KEY *)unpack_dsds(printf)) == NULL) {
//      pemail("***Error on data retrieval\n");
//      abortit(1);
//    }
//  }
  nocontrolc = 0;			/* allow ^C now */
  /* ck the status of each returned ds and print warning if any have
   * a non-0 status. All ds are found from inname_%d.
  */
//  for(i = 0; ; i++) {
//    sprintf(ext, "inname_%d", i);
//    if(!findkey(blist, ext)) break;	/* done, exit for() loop */
//    argname = getkey_str(blist, ext);	/* gets e.g. in_0 */
//    sprintf(ext, "status_%d", i);
//    status = getkey_int(blist, ext);	/* the corresponding status */
//    if(status) {
//      printf("Dsds_svc returned error %d for %s query.\n",status,argname);
//      printf("  **Assume error means a null %s_wd\n", argname);
//      sprintf(ext, "%s_wd", argname);
//      setkey_str(&blist, ext, "");
//    }
//    /* indicate in the keylist that dsds was called for the ds */
//    sprintf(ext, "%s_dsdsin", argname); /* e.g. in_0_dsdsin */
//    setkey_int(&blist, ext, 1);
//  }
  return((KEY *)blist);
}

/* Called by form_keylist() for each output datacollection key in the argument 
 * table for the given server. Outputs the parsed keylist, in a new keylist, 
 * for the dataset. Sets the global vrbl wdkey_int.
*/
KEY *form_arg_data_out(SERVER *sptr, argument *arg)
{
  KEY *xlist, *blist;
  SERVER *sptr2;
  FILE *fin;
  DRMS_Segment_t *segment;
  char path[DRMS_MAXPATHLEN] = {0};
  char *wd, *pname, *lname, *sname, *stmp, *cptr;
  char dbasekey[MAX_STR], ext[MAX_STR], wdkey[MAX_STR], inname[MAX_STR];
  char drmsname[MAX_STR], cmd[MAX_STR], buf[128];
  int innsets, parse, levnum, i, ntmp, dstatus;

  char *prog, *level, *series, *jdata;
  char jpedata[192];
  int seriesnum, jcnt, k;
  int jix = 0;

  //printf("\n*****The sptr->map_list at the start of form_arg_data_out is:\n");
  //keyiterate(printkey, sptr->map_list);

  xlist=newkeylist();
  add_keys(sptr->map_list, &xlist);
  if(!arg_available(xlist, arg->key, arg->kind))  {
    pemail("***The ARG_DATA_OUT \"%s=\" for %s ", arg->key, sptr->name);
    pemail("is missing or\n  a non-allowed default value given\n");
    abortit(1);
  }
  strcpy(dbasekey, arg->key); strcat(dbasekey, "_dbase");
  if(strcmp(dsdswd, "")) 		/* set {dbase} to dsds storage */
    setkey_str(&xlist, dbasekey, dsdswd);
  else if(getenv("dbase"))
    setkey_str(&xlist, dbasekey, getenv("dbase"));
  else
    setkey_str(&xlist, dbasekey, "/tmp"); /* there must be some {dbase} */

  if(parse=parse_list(&xlist, arg->key)) {/* expand ds name */
    pemail("***Error %d in parse_list for %s\n", parse, sptr->name);
    if(dbxflg) {
      printf("\n***** The xlist after the parse is:\n");
      keyiterate(printkey, xlist);
    }
    abortit(1);
  }
  if(dbxflg) {
    printf("\n***** The xlist after the parse is:\n");
    keyiterate(printkey, xlist);
  }
  strcpy(inname, arg->key); strcat(inname, "_nsets");
  innsets = getkey_int(xlist, inname);

  /* If this is an appendable output (a=a) then query wd from ds_naming tbl */
  if(sptr->archive == 'a') {  //!!TBD ck use in jpe
    //NOTE: appendable output ds not supported in DRMS/SUMS.
    //See the original pe.c for how the ds_naming db table was used for this.
    //pemail("***Error: jpe does not allow appendable output dataseries\n");
    //printk("***Error: jpe does not allow appendable output dataseries\n");
    //abortit(1);
    //OK, going to try to support appendable output to /SUM0/PAS !!TEMP
    for(i=0; i < innsets; i++) {          /* do for each dataset */
      /* set this wd as the _dbase term in the keylist and remove the current */
      /* evaluation and parse again to get the final output wd. */
      setkey_str(&xlist, dbasekey, "/SUM0/PAS");
          //!!!TEMP
          //printf("\n***** The xlist before the deletekey() is:\n");
          //keyiterate(printkey, xlist);
      sprintf(ext, "%s_%d_wd", arg->key, i);
      deletekey(&xlist, ext);             /* rem the current _0_wd */
      if(parse=parse_list(&xlist, arg->key)) {/* expand ds name */
        pemail("***Error %d in parse_list for %s\n", parse, sptr->name);
        if(debugflg) {
          printf("\n***** The xlist after the parse is:\n");
          keyiterate(printkey, xlist);
        }
        abortit(1);
      }
    }
          //!!!TEMP
          //printf("\n***** The xlist after the parse is:\n");
          //keyiterate(printkey, xlist);
    return(xlist);
  }

  /* This is not an appendable ds so proceed with rule. */
  /* a rule must be present in the keylist for each dataset */
  for(i=0; i < innsets; i++) {          /* do for each dataset */
    sprintf(ext, "%s_%d_level_sn", arg->key, i);
    if(findkey(xlist, ext)) {
      levnum = getkey_int(xlist, ext);
      if(levnum != 0) {
        pemail("***Map Error: Explicit output level # other than 0 is not allowed\n");
        abortit(1);
      }
      else {
        if(sptr->archive) {
          pemail("***Map Error: Archive of level # 0 output is not allowed\n");
          abortit(1);
        }
      }
    }
    sprintf(ext, "%s_%d_prog", arg->key, i);	//e.g. out_0_prog
    stmp = GETKEY_str(xlist, ext);
    sprintf(drmsname, "%s__", stmp);
    sprintf(ext, "%s_%d_level", arg->key, i);	//e.g. out_0_level
    stmp = GETKEY_str(xlist, ext);
    sprintf(drmsname, "%s%s__", drmsname, stmp);
    sprintf(ext, "%s_%d_series", arg->key, i);	//e.g. out_0_series
    stmp = GETKEY_str(xlist, ext);
    sprintf(drmsname, "%s%s", drmsname, stmp);
    sprintf(cmd, "echo \"%s\" | sed 's/\\\./_/g' | sed 's/-/__/g'", drmsname);
    //printf("cmd = %s\n", cmd);
    fin = popen(cmd, "r");
    fgets(buf, sizeof buf, fin);
    cptr = rindex(buf, '\n');
    if(cptr) *cptr = NULL;
    //sprintf(ext, "%s_%d_series_sn", arg->key, i);
    //ntmp = getkey_int(xlist, ext);
    //sprintf(buf, "%s[%d]", buf, ntmp);
    sprintf(drmsname, "dsds.%s", buf);
    if(findkey(JPElist, "JPE_out_nsets")) {
      jix = getkey_int(JPElist,  "JPE_out_nsets");
    }

  //printf("\nout drmsname = %s\n", drmsname); //!!TEMP

  if(reqmegs) {
    if(firstalloc) {
      // drms create record to get a wd for the output
      rs = drms_create_record(drms_env, drmsname, DRMS_PERMANENT, &dstatus);
      if(dstatus) {
        printk("**ERROR %d: Can't create record for %s\n", dstatus, drmsname);
        abortit(1);
      }
  
      drms_record_directory(rs, path, 0);		//get wd
      //printk("path = %s\n", path);		//!!TEMP
      strcpy(dsdswd, path);			//if next in is prev out
      sprintf(ext, "%s_%d_rs", arg->key, i);
      setkey_fileptr(&sptr->map_list, ext, (FILE *)rs);
      //sptr->rsarray[i] = rs;
      sprintf(ext, "%s_%d_wd", arg->key, i);	//e.g. out_0_wd
      setkey_str(&xlist, ext, path);
      wd = path;
      sprintf(ext, "%s_%d_prog", arg->key, i);      //e.g. out_0_prog
      prog = getkey_str(xlist, ext);                //e.g mdi_eof
      sprintf(ext, "%s_%d_level", arg->key, i);     //e.g. out_0_level
      level = getkey_str(xlist, ext);               //e.g lev1.5
      sprintf(ext, "%s_%d_series", arg->key, i);    //e.g. out_0_series
      series = getkey_str(xlist, ext);              //e.g. vw_V_06h
      sprintf(ext, "%s_%d_series_sn", arg->key, i); //e.g. out_0_series_sn
      seriesnum = getkey_int(xlist, ext);           //e.g. 23220
      sprintf(jpedata, "prog:%s,level:%s,series:%s[%d]",
                prog, level, series, seriesnum);
      //save this info in jpelist for others in this archive group to check
      sprintf(ext, "JPE_%d_data", jix);
      setkey_str(&JPElist, ext, jpedata);
      sprintf(ext, "JPE_%d_wd", jix);
      setkey_str(&JPElist, ext, path);
      sprintf(ext, "JPE_%d_rs", jix);
      setkey_fileptr(&JPElist, ext, (FILE *)rs);
      jix++;
    }
    else {   //see if this out arg is already in JPElist for the arch grp
      sprintf(ext, "%s_%d_prog", arg->key, i);      //e.g. out_0_prog
      prog = getkey_str(xlist, ext);                //e.g mdi_eof
      sprintf(ext, "%s_%d_level", arg->key, i);     //e.g. out_0_level
      level = getkey_str(xlist, ext);               //e.g lev1.5
      sprintf(ext, "%s_%d_series", arg->key, i);    //e.g. out_0_series
      series = getkey_str(xlist, ext);              //e.g. vw_V_06h
      sprintf(ext, "%s_%d_series_sn", arg->key, i); //e.g. out_0_series_sn
      seriesnum = getkey_int(xlist, ext);           //e.g. 23220
      sprintf(jpedata, "prog:%s,level:%s,series:%s[%d]",
   	             prog, level, series, seriesnum);
      jcnt = getkey_int(JPElist, "JPE_out_nsets");
      for(k=0; k < jcnt; k++) {  //try to find this output ds in JPElist
        sprintf(ext, "JPE_%d_data", k);
        jdata = getkey_str(JPElist, ext);
        if(!strcmp(jpedata, jdata)) {	//use this ones wd
          sprintf(ext, "JPE_%d_wd", k);
          wd = getkey_str(JPElist, ext);
          break;
        }
      }
      if(k == jcnt) {			//existing wd not found
        // drms create record to get a wd for the output
        rs = drms_create_record(drms_env, drmsname, DRMS_PERMANENT, &dstatus);
        if(dstatus) {
          printk("**ERROR %d: Can't create record for %s\n", dstatus, drmsname);
          abortit(1);
        }
        drms_record_directory(rs, path, 0);		//get wd
        sprintf(ext, "%s_%d_rs", arg->key, i);
        setkey_fileptr(&sptr->map_list, ext, (FILE *)rs);
        sprintf(ext, "%s_%d_wd", arg->key, i);	//e.g. out_0_wd
        setkey_str(&xlist, ext, path);
        printk("%s path = %s\n", ext, path);		//!!TEMP
        wd = path;
  sprintf(ext, "%s_%d_prog", arg->key, i);      //e.g. out_0_prog
  prog = getkey_str(xlist, ext);                //e.g mdi_eof
  sprintf(ext, "%s_%d_level", arg->key, i);     //e.g. out_0_level
  level = getkey_str(xlist, ext);               //e.g lev1.5
  sprintf(ext, "%s_%d_series", arg->key, i);    //e.g. out_0_series
  series = getkey_str(xlist, ext);              //e.g. vw_V_06h
  sprintf(ext, "%s_%d_series_sn", arg->key, i); //e.g. out_0_series_sn
  seriesnum = getkey_int(xlist, ext);           //e.g. 23220
  sprintf(jpedata, "prog:%s,level:%s,series:%s[%d]",
                prog, level, series, seriesnum);
  //save this info in jpelist for others in this archive group to check
  sprintf(ext, "JPE_%d_data", jix);
  setkey_str(&JPElist, ext, jpedata);
  sprintf(ext, "JPE_%d_wd", jix);
  setkey_str(&JPElist, ext, path);
  jix++;
      }
      else {			//use the wd from JPElist
        sprintf(ext, "JPE_%d_wd", k);
        wd = getkey_str(JPElist, ext);
        sprintf(ext, "%s_%d_wd", arg->key, i);    //e.g. out_0_wd
        setkey_str(&xlist, ext, wd);
        sprintf(ext, "JPE_%d_rs", k);
        rs = getkey_fileptr(JPElist, ext);
        sprintf(ext, "%s_%d_rs", arg->key, i);
        setkey_fileptr(&sptr->map_list, ext, (FILE *)rs);
      }
    }
  }

    if(!wd || !strcmp(wd, "")) {
      pemail("***No %s for %s. Missing rule or template in map file.\n",
	 ext, sptr->name);
      abortit(1);
    }
    if(*wd == '.') {
      pemail("***%s of \".\" not valid for %s.\n", ext, sptr->name);
      abortit(1);
    }
    /* check if this is a duplicate dir that we already have */
    if(!wd_dup_ck(wd_dup_list, wd)) {	/* no dup wd so add it to list */
      sprintf(wdkey, "wd_%d", wdkey_int++);/* form the unique name of new key */
      setkey_str(&wd_dup_list, wdkey, wd);
    }
    else {
      if(!nowarn) {
      printf("**WARNING: pe detects two datasets that have the same output wd:\n");
      printf("  Make sure that you have run with START_ARCHIVE/END_ARCHIVE\n");
      printf("  wd = %s\n", wd);
      }
    }
  }
  setkey_int(&JPElist, "JPE_out_nsets", jix);
  firstalloc = 0;
  if(dbxflg) {
    printf("\n***** The xlist at the end of form_arg_data_out is:\n");
    keyiterate(printkey, xlist);
    printf("\n##### The JPElist at the end of form_arg_data_out is:\n");
    keyiterate(printkey, JPElist);
  }
  return(xlist);
}

/* Called by form_keylist() after the input ds have been resolved with
 * dsds_svc and if the user has requested an rcp of the datasets by 
 * specifying "rcp=/rdir" on the command line.
 * Called with the keylist of the resolved input ds where arg_data_in_0,
 * arg_data_in_1, etc. give the prefix name for each input ds. For example,
 * arg_data_in_0 = "in" and then "in_nsets" gives the number of ds for
 * this arg and "in_0_wd", "in_1_wd", etc. will have the wds.
 * These wds will be copied to the /rdir given by the caller and the
 * keylist will be modified to reflect these dirs as containing the 
 * input ds.
*/
void in_ds_rcp(KEY *list)
{
  int i, loop, innsets;
  char *argname, *wd, *cptr;
  char ext[MAX_STR], inname[MAX_STR], cmd[MAX_STR];
  char targetdir[MAX_STR];
  struct stat stbuf;

  for(loop = 0; ; loop++) {
    sprintf(ext, "arg_data_in_%d", loop);
    if(!findkey(list, ext)) break;		/* all done, exit for(loop)*/
    argname = getkey_str(list, ext);		/* e.g. "in" */
    strcpy(inname, argname); strcat(inname, "_nsets");
    innsets = getkey_int(list, inname);
    for(i=0; i < innsets; i++) {        /* do for each dataset */
      sprintf(ext, "_%d_wd", i);
      strcpy(inname, argname); strcat(inname, ext); /* e.g. "in_0_wd" */
      wd = getkey_str(list, inname);
      if(strcmp(wd, "")) {
        cptr = index(wd+1, '/');	/* pos to e.g. /PDS12/ */
        sprintf(targetdir, "%s%s", rcpdir, cptr);
        if(stat(targetdir, &stbuf)) {	/* dir doesn't exist */
          sprintf(cmd, "mkdir -p %s", targetdir);
          if(system(cmd)) {
            pemail("***Cannot %s\n", cmd);
            abortit(1);
          }
        }
        sprintf(cmd, "rcp -p %s:%s/* %s", prod_host_prime(), wd, targetdir);
        if(system(cmd)) {
          pemail("***Cannot %s\n", cmd);
          abortit(1);
        }
        printf("%s\n", cmd);
        setkey_str(&list, inname, targetdir);
      }
    }
  }
}

/* Called whenever the next server in the map file is about to be run for
 * the first time. 
 * Forms the keylist for the given server in the param_list for each host
 * from the map file keylist and the argument list from the server for all 
 * datasets. 
 * The map file key list was set in sptr->map_list by get_cmd().
 * The argument list was set in sptr->arguments by arg_recv().
 * If d=1 in the map file for this server then will query dsds_svc for the 
 * wd and fsn/lsn for each in dataset. However, if an exlicit [fsn-lsn] was 
 * given in the map file, then the values from the database will not override
 * the given values (fsn/lsn now obsolete in the db). 
 * Calls form_split() (now obsolete) for the first input dataset only in the 
 * argument list which will potentially distribute the dataset execution for 
 * the fsn-lsn over a number of hosts (now obsolete). This sets the global 
 * variable effective_hosts (now always set to 1).
 * Finally calls ck_arglist() for each host that is to run the server.
 *   sptr = SERVER table of new server to run
 *
 * NEW 11/14/96: All the ARG_DATA_IN args are built into a single keylist
 * which is sent to dsds_svc to resolve the wd and bring online if needed.
 * Formerly dsds_svc was called for each single dataset for each ARG_DATA_IN.
 * Now dsds_svc will get one msg and bring all the input datasets that this
 * server needs online before it returns the answer here with the wd's 
 * resolved. (This is to allow dsds_svc to optimize the retrieve from 
 * tape of all the datasets that a module needs. This will greatly speed
 * things up.)
 *
 * NEW Apr 11, 1997: Add the *hnext as an input arg and don't loop on it
 * anymore as the effective_hosts will always be 1.
*/
void form_keylist(SERVER *sptr, HDATA *hnext)
{
  KEY *xlist, *inlist, *dslist;
  argument *arg, *argsort[MAX_ARGS];
  char ext[MAX_STR];
  char *pname;
  int inseq, i, ax;
  int inapp = 0;

  /* make sure we see any ARG_DATA_IN first */
  arg=sptr->arguments;
  i = 0;
  while(arg->kind != ARG_END) {
    if(arg->kind == ARG_DATA_IN) argsort[i++] = arg;
    arg++;
  }
  arg=sptr->arguments;
  while(arg->kind != ARG_END) {
    if(arg->kind != ARG_DATA_IN) argsort[i++] = arg;
    arg++;
  }
  argsort[i] = arg;			/* include ARG_END */

  /* verify and manipulate the ds names given in the arg list */
  inseq = 0; inlist = newkeylist(); dslist = newkeylist();
  ax = 0; arg = argsort[ax];
  while(arg->kind != ARG_END) {		/* look at all args */
    switch(arg->kind) {
    case ARG_DATA_IN:
      inlist = (KEY *)form_arg_data_in(inlist, sptr, arg, inseq);
      sprintf(ext, "%s_prog", arg->key);    /* e.g. in_prog */
      pname = getkey_str(inlist, ext);  
      //check if this is a  dsds appendable input series
      if(!strcmp(pname, "mdi_eof_log") || !strcmp(pname, "mdi_eof_rec")
         || !strcmp(pname, "mdi_log") || !strcmp(pname, "mdi_rec") ||
         !strcmp(pname, "mdi_sim_rec") || !strcmp(pname, "soi_eof_rec")) {
        inapp = 1;
      }
      if(!inseq++) {				/* first input dc */
        /* NOTE: form_split() is vestigial. Just set effective_hosts */
        /*form_split(sptr, inlist, arg->key);	/* sets effective_hosts */
        /*effective_hosts = 1;*/
      }
      add_keys(inlist, &hnext->param_list);/* copy over all keys */
      break;
    case ARG_DATA_OUT:
      xlist = (KEY *)form_arg_data_out(sptr, arg);
      add_keys(xlist, &hnext->param_list);/* copy over all keys */
      freekeylist(&xlist);
      break;
    default:
      break;
    }
    arg = argsort[++ax];
  }					/* end arg->kind != ARG_END */
  if((inseq) && (sptr->dsin) && (!inapp)) { //ARG_DATA_IN present & call dsds
    if(dbxflg) {			/* !!TEMP for debug */
      printf("\nThe list before the dsds_arg_data_in(inlist) call:\n");
      keyiterate(printkey, inlist);
    }
    dslist = (KEY *)dsds_arg_data_in(inlist);	/* call dsds for final list */
    if(dbxflg) {			/* !!TEMP for debug */
      printf("\nThe list after the dsds_arg_data_in(inlist) call:\n");
      keyiterate(printkey, dslist);
    }
    if(rcpflg) {			/* user requested rcp of input */
      /* rcp all the ds input wd to the dir given in the rcp= arg */
      /* and then modify the wd in the list. */
      in_ds_rcp(dslist);
    }
    add_keys(dslist, &hnext->param_list);/* copy over all keys */
    //printf("\n*** The param_list after add_keys() is:\n");
    //keyiterate(printkey, hnext->param_list); //!!!TEMP
    freekeylist(&inlist); freekeylist(&dslist);
  }
  else if(inseq) {
    add_keys(inlist, &hnext->param_list);/* copy over all keys */
    freekeylist(&inlist);
  }
  /* Now satisfy and expand the arg list for all hosts for this server */
  ck_arglist(sptr, hnext);		/* expand and ck args */
  if(findkey(hnext->param_list, "T_FIRST")) {
    t_first_arg = getkey_str(hnext->param_list, "T_FIRST");
  } 
  else {
    t_first_arg = NULL;
  }
  if(findkey(hnext->param_list, "T_LAST")) {
    t_last_arg = getkey_str(hnext->param_list, "T_LAST");
  } 
  else {
    t_last_arg = NULL;
  }
}

/* Pack and send a message to a server. The message will consists of the
given doflg and the keylist for this server contained in the stab[].
Called with the stab[] address of the server, the hosts entry in stab[],
and the value of doflg.
Will pack the message, set the host busy and send the msg to the server.
*/
void send_to_serv(SERVER *sptr, HDATA *hx, int doflg)
{
  KEY *keybad;
  argument *arg;
  char cfsn[MAX_STR], clsn[MAX_STR];

  hx->busy=1;                   /* set this host busy */
  sptr->busyall++;              /* inc master busy count */
  arg=sptr->arguments;
  while(arg->kind != ARG_END) {
    if(arg->kind == ARG_DATA_IN) {
      sprintf(cfsn, "%s_0_fsn", arg->key);
      if(!sptr->split || !findkey(hx->param_list, cfsn)) {
        printf("%s %s sent:\n  %s=%s\n",hx->host_name,
                sptr->name,arg->key,getkey_str(hx->param_list,arg->key));
        break;
      }
      else {            /* put split info in the msg */
        sprintf(clsn, "%s_0_lsn", arg->key);
        printf("%s %s sent:\n  %s=%s split[%d-%d]\n",hx->host_name,
                sptr->name,arg->key,getkey_str(hx->param_list,arg->key),
                getkey_int(hx->param_list, cfsn),
                getkey_int(hx->param_list, clsn));
        break;
      }
    }
    arg++;
  }
  if(!sptr->cphist)
    setkey_int(&hx->param_list, "nocphist", 1);  /* don't cp hist log */
  setkey_int(&hx->param_list, "dsds_tid", dsds_tid); /* always add dsds_tid */
  setkey_uint(&hx->param_list, "dsds_uid", uid); /* always add dsds_uid */
  setkey_str(&hx->param_list, "pe_mapfile", pe_map); /* add pe map name */
  pvm_initsend(PvmDataDefault); /* init the send buffer */
  pvm_pkint(&doflg, 1, 1);      /* pack the do flag */
  if(keybad=(KEY *)pack_keylist(hx->param_list)) {
    pemail("***Err packing a pvm msg, type=%d name=%s\n",keybad->type,keybad->name);
    abortit(1);
  }
  printf("Sending to server:\n");  /* !!!TEMP */
  keyiterate(printkey, hx->param_list);
  if(pvm_send(hx->tid, sptr->msgid)) {/* send the msg to server */
    pemail("***Error calling %s on %s\n",sptr->name,hx->host_name);
    abortit(1);
  }
}


/* This routine is normally called by call_archive() after all the output 
 * datasets of a server have been archived. It will output into each dataset's 
 * wd, a text file that has the full 6 part name of all the actual input and 
 * output datasets that were used by the server. 
 * Also changes the files of each wd to be owned by production w/mode 644.
 *   list = keylist sent to the server and update by pe with the level numbers
 *   arg = argument table from the server
 * Returns non-0 on error.
*/
int ds_names_file(KEY *list, argument *argu)
{
  AT *atp;
  argument *arg;
  char ext[MAX_STR], rdbline[1024], sysstr[1024];
  char *prog, *level, *series, *wd;
  double bytes;
  unsigned long dsindex;
  int nsets, levnum, prognum, seriesnum;
  int i = 0;
  char DSINFO[] ="arg\tprog\tprog#\tlevel\tlev#\tseries\tser#\tdsindex\tbytes\n----\t-----\t------\t------\t-----\t-------\t-----\t-------\t------\n";
  if(!(atp = at_create(DSINFO))) {
    printf("**Can't at_create() the table for ds.rdb file\n");
    return(1);
  }
  arg = argu;
  while(arg->kind != ARG_END) {
    if(arg->kind == ARG_DATA_IN || arg->kind == ARG_DATA_OUT) {
      sprintf(ext, "%s_nsets", arg->key);
      nsets = getkey_int(list, ext);
      for(i=0; i < nsets; i++) {
        sprintf(ext, "%s_%d_prog", arg->key, i);
        if(findkey(list, ext)) prog = getkey_str(list, ext);
        else prog = "(null)";
        sprintf(ext, "%s_%d_prog_sn", arg->key, i);
        if(findkey(list, ext)) 
          prognum = getkey_int(list, ext);
        else
          prognum = -1;
        sprintf(ext, "%s_%d_level", arg->key, i);
        if(findkey(list, ext)) level = getkey_str(list, ext);
        else level = "(null)";
        sprintf(ext, "%s_%d_level_sn", arg->key, i);
        if(findkey(list, ext))
          levnum = getkey_int(list, ext);
        else
          levnum = -1;
        sprintf(ext, "%s_%d_series", arg->key, i);
        if(findkey(list, ext)) series = getkey_str(list, ext);
        else series = "(null)";
        sprintf(ext, "%s_%d_series_sn", arg->key, i);
        if(findkey(list, ext))
          seriesnum = getkey_int(list, ext);
        else
          seriesnum = -1;
        sprintf(ext, "%s_%d_ds_index", arg->key, i);
        if(findkey(list, ext))
          dsindex = getkey_ulong(list, ext);
        else
          dsindex = (unsigned long)-1;
        sprintf(ext, "%s_%d_bytes", arg->key, i);
        if(findkey(list, ext))
          bytes = getkey_double(list, ext);
        else
          bytes = 0;
        sprintf(ext, "%s_%d", arg->key, i);
        sprintf(rdbline, "%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%g",
		ext,prog,prognum,level,levnum,series,seriesnum,dsindex,bytes);
        if(at_put_info(atp, rdbline) != AT_OK) {
          printf("**Error on at_put_info() to ds.rdb file\n");
          return(1);
        }
      }
    }
    arg++;
  }
  /* now write out the ds.rdb file to each output wd */
  arg = argu;
  while(arg->kind != ARG_END) {
    if(arg->kind == ARG_DATA_OUT) {
      sprintf(ext, "%s_nsets", arg->key);
      nsets = getkey_int(list, ext);
      for(i=0; i < nsets; i++) {
        sprintf(ext, "%s_%d_wd", arg->key, i);
        wd = getkey_str(list, ext);
        if(strcmp(wd, "")) {			/* a wd is present */
          sprintf(ext, "%s/ds.rdb", wd);
          if(at_write(atp, ext) != AT_OK) {
            printf("**Error on at_write() to %s\n", ext);
            return(1);
          }
          /* make production owner 755 mode for mdi_2 only */
            /*3Nov00 Don't do it this way, leaves "." group writable*/
            /*sprintf(sysstr, "chmod -R go-w %s/*; chown -Rhf production %s",*/
            sprintf(sysstr, "chmod -R go-ws %s; chown -Rhf production %s", 
                    wd, wd);
            if(system(sysstr)) {
              printf("**Warning: Error on: %s\n", sysstr);
              /*return(1);*/
            }
        }
      }
    }
    arg++;
  }
  at_free(atp);
  return(0);
}

/* Look for an overview.fits file in the given working directory. If found
 * extract the t_first and t_last strings from it and sets them as "t_first"
 * and "t_last" keys in the given keylist. Returns 0 on error.
 * NEW: 4Dec01 If no overview.fits is present then use the values of t_first 
 * and t_last if they were sent to the modules as args as indicated by 
 * t_first_arg and t_last_arg.
*/
int get_overview(char *wd, KEY **list)
{
  SDS *sds;
  DIR *dfd;
  struct dirent *dp;
  char *value;
  char name[512];
  int found = 0;

  if((dfd=opendir(wd)) == NULL) {
    printf("**Can't opendir(%s) to find overview.fits\n", wd);
    return(0);
  }
  while((dp=readdir(dfd)) != NULL) {
    if(!strstr(dp->d_name, "overview.fits"))
      continue;
    found = 1;
    sprintf(name, "%s/%s", wd, dp->d_name);
    closedir(dfd);
/*    if((sds=sds_read_FITS_header(name)) == NULL) {  */
    if ((sds = sds_get_fits_head (name)) == NULL) {
      printf("**Error on get_fits_head(%s)\n", name);
      return(0);
    }
    if((value=sds_search_attrvalue_str(sds, "T_FIRST")) == NULL) {
      /*printf("**No t_first in %s\n", name); */
      return(0);
    }
    if(!strcmp(value, ""))
      setkey_str(list, "t_first", "-1");
    else
      setkey_str(list, "t_first", value);
    if((value=sds_search_attrvalue_str(sds, "T_LAST")) == NULL) {
      printf("**No t_last in %s\n", name);
      return(0);
    }
    if(!strcmp(value, ""))
      setkey_str(list, "t_last", "-1");
    else
      setkey_str(list, "t_last", value);
    break;
  }
  if(found)
    return(1);
  else {
    closedir(dfd);
    if(t_first_arg) {
      setkey_str(list, "t_first", t_first_arg);
      if(t_last_arg) {
        setkey_str(list, "t_last", t_last_arg);
      }
      else {
        setkey_str(list, "t_last", "-1");
      }
      return(1);
    }
    if(t_last_arg) {
      setkey_str(list, "t_first", "-1");
      setkey_str(list, "t_last", t_last_arg);
      return(1);
    }
    return(0);
  }
}

/* Set up a keylist for a call to REQUPDATE to dsds_svc.
 * Creates a new keylist and populates it and returns it.
 *
 *  basename = the base key name e.g. out_0
 *  stab =  pointer to the completing server table
 *  hdata = host data structure pointer that contains param_list which is
 *          the keylist originally passed to the server
 *  rlist = keylist that was returned from the server that has the fsn & lsn
 *          for each output dataset
 *  status= the status code returned by the completing server
*/
KEY *set_key_archive(char *basename, SERVER *stab, HDATA *hdata,  KEY *rlist,
			int status)
{
  KEY *alist;
  char *wd, *effective_date, *warnmsg;
  char ext[MAX_STR];
  double dsize;
  /*int fsn, lsn, nofsnlsn;*/

  alist=newkeylist();
  sprintf(ext, "%s_wd", basename);
  wd = getkey_str(hdata->param_list, ext);
  /* try to find t_first/t_last in overview.fits file in wd */
  get_overview(wd, &alist);
  setkey_int(&alist, "status", status);
  if((warnmsg=getkey_str(rlist, "WARNING")))
    setkey_str(&alist, "warning", warnmsg);

/*  sprintf(ext, "%s_fsn", basename);
/*  fsn = 0; lsn = 0; nofsnlsn = 0;
/*  if(!findkey(rlist, ext))
/*    nofsnlsn = 1;
/*  else
/*    fsn = getkey_int(rlist, ext);
/*  setkey_int(&alist, "fsn", fsn);
/*  sprintf(ext, "%s_lsn", basename);
/*  if(!findkey(rlist, ext))
/*    nofsnlsn = 1;
/*  else
/*    lsn = getkey_int(rlist, ext);
/* !!TBD: eventually elim setting fsn & lsn in keylist
/*  if(nofsnlsn)
/*    printf("No %s fsn and/or lsn from %s. Assume 0\n", basename, stab->name);
/*  setkey_int(&alist, "lsn", lsn);
*/
  setkey_str(&alist, "wd", wd);
  sprintf(ext, "%s_prog", basename);
  setkey_str(&alist, "prog", getkey_str(hdata->param_list, ext));
  sprintf(ext, "%s_series", basename);
  setkey_str(&alist,"series", getkey_str(hdata->param_list, ext));
  sprintf(ext, "%s_prog_sn", basename);
  if(!findkey(hdata->param_list, ext))
    setkey_int(&alist, "prog_sn", -1);
  else
    setkey_int(&alist,"prog_sn", getkey_int(hdata->param_list, ext));
  sprintf(ext, "%s_series_sn", basename);
  if(!findkey(hdata->param_list, ext))
    setkey_int(&alist, "series_sn", -1);
  else
    setkey_int(&alist,"series_sn", getkey_int(hdata->param_list, ext));
  sprintf(ext, "%s_level", basename);
  setkey_str(&alist,"level", getkey_str(hdata->param_list, ext));
  sprintf(ext, "%s_level_sn", basename);
  if(findkey(hdata->param_list, ext))
    setkey_int(&alist,"level_sn", getkey_int(hdata->param_list, ext));
  dsize = du_dir(wd);                     /* get #bytes of storage */
  /* add 10% to storage to keep in sync with what df reports */
  /*dsize = dsize + (dsize * 0.1);*/
  setkey_int(&alist, "lago_tid", lago_tid);
  setkey_double(&alist, "dsds_bytes", dsize);
  setkey_uint(&alist, "dsds_uid", uid);   /* give our OPEN uid */
  setkey_str(&alist, "svc_name", stab->name);
  setkey_str(&alist, "svc_version", stab->version);
  setkey_str(&alist, "username", username);
  return(alist);
}

/* Called by process_resp() when  data products just created by a completing
 * server are to be archived to DSDS.
 * Calls dsds_svc with a catalog update message for each output ds.
 *  stab =  pointer to the completing server table
 *  hdata = host data structure pointer that contains param_list which is
         *          the keylist originally passed to the server
         *  rlist = keylist that was returned from the server that has the fsn & lsn
 *          for each output dataset
 *  status= the status code returned by the completing server
*/
void call_archive(SERVER *stab, HDATA *hdata, KEY *rlist, int status)
{
  KEY *alist, *blist;
  DRMS_Segment_t *segment;
  DRMS_Record_t *rsx;
  argument *arg;
  char *wd, *cptr;
  char oname[MAX_STR], ext[MAX_STR];
  int outnsets, i, snum;
  unsigned long sunum;
  double dsize;

  archactive = 1;                       /* set flg for archive in progress */
  arg = stab->arguments;
  while(arg->kind != ARG_END) {         /* find all output datasets */
    if(arg->kind == ARG_DATA_OUT) {
      sprintf(oname, "%s_nsets", arg->key);
      outnsets = getkey_int(hdata->param_list, oname);
      for(i=0; i < outnsets; i++) {     /* do for each dataset */
        sprintf(oname, "%s_%d", arg->key, i);
        /* form the keylist for the REQUPDATE call */
        alist = set_key_archive(oname, stab, hdata, rlist, status);
        wd = getkey_str(alist, "wd");
        dsize = getkey_double(alist, "dsds_bytes");
        //printf("\ncall_archive() hdata->param_list=\n");
       // keyiterate(printkey, hdata->param_list); //!!TEMP

      sprintf(oname, "%s_%d_rs", arg->key, i);
      rsx = (DRMS_Record_t *)getkey_fileptr(stab->map_list, oname);
      //rsx = stab->rsarray[i];
      snum = getkey_int(alist, "series_sn");
      //printf("%%%% snum in call_archive() = %d\n", snum); //!!TEMP
      status = drms_setkey_int(rsx, "snum", snum);
      if(status) {
        printk("ERROR: can't drms_setkey_int() for snum=%u\n", snum);
        abortit(1);
      }

  //status= rsx->seriesinfo->archive;  //!!TEMP for test
  //printf("The archive flag is %d\n", status);
  sunum = rsx->sunum;
  rsx->seriesinfo->retention = stab->archive_day;
  if(!stab->archive || stab->archive == 't') 	//don't archive if 0 or temp
    rsx->seriesinfo->archive = 0;
  else {
    rsx->seriesinfo->archive = 1;
    stab->archive_complete = 1;		/* don't del wd at end */
  }

  //printk("Call: drms_close_record()\n"); //!!TEMP
  if((status = drms_close_record(rsx, DRMS_INSERT_RECORD))) {
    printk("**ERROR: drms_close_record failed status=%d\n", status);
    abortit(1);
  }
  sprintf(oname, "%s_%d_bytes", arg->key, i);
  setkey_double(&hdata->param_list, oname, dsize); /* also add # bytes */
  sprintf(oname, "%s_%d_ds_index", arg->key, i);
  setkey_ulong(&hdata->param_list, oname, sunum); // and ds_index

/* !!TEMP print the servers keylist */
/*printf("In call_archive() print the servers keylist for %s is:\n",stab->name);
/*keyiterate(printkey, hdata->param_list);
/* end TEMP */

        freekeylist(&alist); 
        printf("Archive pending: wd=%s bytes=%g\n", wd, dsize);
      }					/* end for each dataset */
    }					/* end ARG_DATA_OUT */
    arg++;
  }					/* end ARG_END */
  stab->archive_complete = 1;           /* don't del wd at end */
  ds_names_file(hdata->param_list, stab->arguments);//write out all ds names
  if(archactive == -1) {		/* ^C during the archive */
    abortit(1);				/* ok to abort now */
  }
  archactive = 0;			/* arch no longer active */
}

/* Called by process_resp() when a completing server is marked as 
 * NO archive to DSDS but DSDS storage was obtained. There might be some
 * level number 0 datasets that we want to catalog. If there are no level #0
 * datasets then this call is a noop.
 * This will check for any level number 0 datasets and make an entry in the
 * dsds_main table for them, but with no archive pending status. They will
 * have delete pending status.
 * Calls dsds_svc with a catalog update 0 message for each output ds level# 0.
 *  stab =  pointer to the completing server table
 *  hdata = host data structure pointer that contains param_list which is 
 *          the keylist originally passed to the server
 *  rlist = keylist that was returned from the server that has the fsn & lsn
 *          for each output dataset
 *  status= the status code returned by the completing server
*/
void call_archive_0(SERVER *stab, HDATA *hdata, KEY *rlist, int status)
{
  KEY *alist, *blist;
  argument *arg;
  char *wd;
  char oname[MAX_STR], ext[MAX_STR];
  int outnsets, found, i;
  double dsize;

  arg = stab->arguments; found = 0;
  while(arg->kind != ARG_END) {		/* find all output datasets */
    if(arg->kind == ARG_DATA_OUT) {
      sprintf(oname, "%s_nsets", arg->key);
      outnsets = getkey_int(hdata->param_list, oname);
      for(i=0; i < outnsets; i++) {	/* do for each dataset */
        sprintf(oname, "%s_%d", arg->key, i);
        sprintf(ext, "%s_level_sn", oname);
        if(findkey(hdata->param_list, ext)) {
          if(getkey_int(hdata->param_list,ext) != 0) /* not a lev#0 ds */
            continue;			/* continue the for() for next ds */
        }
        else continue;			/* no lev# so not lev#0 */
        found = 1;
        stab->archive_day = 2;		/* make del pend in 2 day */
        /* form the keylist for the REQUPDATE0 call */
        alist = set_key_archive(oname, stab, hdata, rlist, status);
        wd = getkey_str(alist, "wd");
        dsize = getkey_double(alist, "dsds_bytes");
        if((blist = (KEY *)call_dsds(&alist, REQUPDATE0, dsds_tid, pe_tid, printf, debugflg)) == NULL) {
          pemail("***Error making a dsds_main entry for  wd=%s\n", wd);
          abortit(1);
        }
        printf("level#0 entry: wd=%s bytes=%g\n", wd, dsize);

        /* update the servers keylist with the ds_index used (level# is 0)*/
        if(findkey(blist, "ds_index")) {
          sprintf(ext, "%s_ds_index", oname);
          setkey_ulong(&hdata->param_list,ext,getkey_ulong(blist, "ds_index"));
        }
        sprintf(ext, "%s_bytes", oname);
        setkey_double(&hdata->param_list, ext, dsize); /* also add # bytes */

/* !!TEMP print the servers keylist */
/*printf("In call_archive_0() print the servers keylist for %s is:\n",stab->name);
/*keyiterate(printkey, hdata->param_list);
/* end TEMP */

        freekeylist(&alist); freekeylist(&blist);
      }					/* end for each dataset */
    }					/* end ARG_DATA_OUT */
    arg++;
  }					/* end ARG_END */
  if(found)
    ds_names_file(hdata->param_list, stab->arguments);/* write all ds names */
}


/* Start the given server for the number of hosts. Effective_hosts is
 * set by form_keylist() (now obsolete. effective_hosts always 1).
 * New on Apr 11, 1997 is to spawn the servers here and send them the
 * message id that they will use and then get their arg list.
*/
void kick_server(SERVER *sptr)
{
  HDATA *hnext;
  int i;

  if(!sptr->archive_group) {
    firstalloc = 1;		 //!!!TBD check this
    JPE_out_nsets = 0;
    if(JPElist) freekeylist(&JPElist);
    JPElist = newkeylist();
  }
  hnext = (HDATA *)gethnext(sptr->hosts);/* get first hdata entry */
  for(i=0; i<effective_hosts; i++) {    /* spawn server on ea host */
    //if(debugflg)
    //  pvm_spawn(sptr->name, (char **)0, PvmTaskDebug, hnext->host_name, 1,
    //          &hnext->tid);
    //else
      pvm_spawn(sptr->name, (char **)0, PvmTaskHost, hnext->host_name, 1,
              &hnext->tid);
    if(hnext->tid < 0) {
      pemail("***Can't spawn %s on %s\n",sptr->name,hnext->host_name);
     abortit(1);
    }
    printk("%s tid=%x spawned on %s\n",
              sptr->name, hnext->tid, hnext->host_name);
    msgid_send(sptr, hnext);            /* send the msd id to use */
    arg_recv(sptr, hnext);              /* get the servers arg list */
    //printk("\nThe server sends back this arg list\n"); //!!!TEMP
    //keyiterate(printkey, hnext->param_list);
    form_keylist(sptr, hnext);          /* make keylist for the args in sptr */
    send_to_serv(sptr,hnext,0);         /* send the processing msg */
    hnext = (HDATA *)gethnext((HDATA *)MONE);/* get next hdata entry */
    /* the next hdata is vestigial as effective_hosts now always 1 */
  }
}

/* Wait for any server to give a completion message of its file processing.
 * Returns the stab index of the completing server, else rets -1 if timeout.
*/
int resp_any()
{
  struct timeval tvalr;
  uint64_t tsr;
  int retcode, bufid, bytes, msgtag, tid, i;

  /* Timeout in case never get a completion response */
  gettimeofday(&tvalr, NULL);
  tsr = tvalr.tv_sec;
  while(1) {
    gettimeofday(&tvalr, NULL);
    if(tvalr.tv_sec - tsr > RESPWAIT) {
      retcode = -1;                     /* give failure ret code */
      break;;                           /* exit while loop */
    }
    if(bufid=pvm_nrecv(-1, -1)) {       /* get any completion from any server */
      if(pvm_bufinfo(bufid, &bytes, &msgtag, &tid)) { /* get msg id */
        pemail("***Can't get bufinfo on bufid=%d\n", bufid);
        abortit(1);
      }
      for(i=0; i<MAX_SERV; i++) {
        if(stab[i].name == NULL) {      /* no more servers */
          pemail("***Got a bad msg type response %d from tid=%x\n", msgtag, tid)
;
          abortit(1);
        }
        if(stab[i].msgid == msgtag)
          break;                        /* exit for loop */
      }
      retcode=i;
      break;                            /* exit while loop */
    }                                   /* end if get any completion */
#ifdef __sgi
    sginap(1);                          /* give cpu 1/3 sec break */
#else
    sleep(1);                           /* give the cpu a break */
#endif
  }                                     /* end while(1) */
  return(retcode);
}

/* Process the response message for the given stab index.
 * The response is tid int, status int, host name, server output keylist.
 * Archive the output datasets if indicated in the map file.
 * Send the keylist of the next server type to the next server as input.
 * NOTE: Must be re-entrant
*/
void process_resp(int ix)
{
  SERVER *stabp, *stabn;
  HDATA *hnext;
  KEY *list;
  argument *arg;
  char remhost[MAX_STR], sysstr[MAX_STR], ext[MAX_STR];
  char *errtxt, *remoutfile, *outwd, *cptr, *warnmsg, *auxinfo;
  int j, remerrno, remtid, wait, groupid, sindex, outs, abortnum;

  stabp = &stab[ix];
  pvm_upkint(&remtid, 1,1 );
  pvm_upkint(&remerrno, 1,1 );
  pvm_upkstr(remhost);
  list=newkeylist();
  if(!(list=(KEY *)unpack_keylist(list))) {     /* unpack msg into keylist */
    printk("**Warning: no keylist returned from %s tid=%x\n",stabp->name,remtid);
  }
  /*printk("\n");                          /* seperate ea completion */
  if(!(remoutfile=getkey_str(list, "out")))
    remoutfile = "<NONE>";
  if(remerrno) {                        /* the server indicated an error */
    switch(remerrno) {                  /* will be replaced by err str rte */
    case MISSING_FILES:
      errtxt = "MISSING_FILES";
      break;
    default:
      errtxt = " ";
      break;
    }
    printk("    %s %x status %d on %s:%s\n", stabp->name, remtid, remerrno, remhost, errtxt);
  }

/*  printf("The keylist returned by %s is:\n", stabp->name);
/*  keyiterate(printkey, list);
*/

  hnext = (HDATA *)gethnext(stabp->hosts);/* get first hdata entry */
  while(hnext) {                        /* find the task that responded */
    if(hnext->tid == remtid) {
      hnext->busy = 0;                  /* indicate this host is free */
      stabp->busyall--;                 /* decrement master busy count */
      /* copy the map file to all the output datasets */
      arg = stabp->arguments;
    if(stabp->cphist) {                 /* ok to cp map to output wds */
      while(arg->kind != ARG_END) {
        if(arg->kind == ARG_DATA_OUT) {
          sprintf(ext, "%s_nsets", arg->key);
          if(findkey(hnext->param_list, ext)) {
            outs = getkey_int(hnext->param_list, ext); /* #of output ds */
            for(j=0; j < outs; j++) {
              sprintf(sysstr, "%s_%d_wd", arg->key, j);
              if(outwd = getkey_str(hnext->param_list, sysstr)) {
                if(cptr = rindex(pe_map, '/'))
                  cptr=cptr+1;          /* point to file name only */
                else
                  cptr=pe_map;
                sprintf(sysstr, "cp %s %s;chmod u+w,g+w %s/%s",
                        pe_map, outwd, outwd, cptr);
                if(system(sysstr)) {    /* copy map file to the output wd */
                  pemail("***Error on: %s\n", sysstr);
                  abortit(1);
                }
              }
            }
          }
        }
        arg++;
      }                                 /* end while(arg->kind != ARG_END) */
    }
      if(!(auxinfo=getkey_str(list, "AUXINFO")))
        auxinfo = "<NONE>";
      if(findkey(list, "abortflg")) {
        if((abortnum=getkey_int(list, "abortflg"))) {
          if(groupid = stabp->groupid) {/* module is part of a group */
            pemail("    ***Module %s abort. See %s\n",
                stabp->name, getkey_str(list, "log_file"));
            /*if((warnmsg=getkey_str(list, "WARNING"))) {
            /*  printk("    with warning msg:\n    %s\n", warnmsg);
            /*}*/
            pemail("    Not aborting as module is part of a group\n");
          }
          else {
            pemail("    ***Module %s abort. See %s\n",
                stabp->name, getkey_str(list, "log_file"));
            if (stabp->noabort) {    /*  map file had ABORT_ACTION=continue  */
/*              pemail ("    Not aborting as map file had ABORT_ACTION=CONTINUE
override\n"); */
              if(abortnum == 1) {       /* don't archive at end */
                if(stabp->archive != 'a') { /* unless appendable */
                  stabp->archive = 0;
                }
              }
            }
            else {
              if((warnmsg=getkey_str(list, "WARNING"))) {
                pemail("    with warning msg:\n    %s\n", warnmsg);
              }
              pelog("%x\t%d\tabort\t0\t%s\t%s\t%s\n",
                remtid,remerrno,stabp->name,remoutfile,auxinfo);
              /* catalog any level #0 ds if storage was obtained */
              if(!stabp->archive && reqmegs)
                call_archive_0(stabp, hnext, list, remerrno);
              /* archive anyway if abortnum was 2 */
              if(abortnum == 2) {
                if(stabp->archive && (stabp->archive_group != 1)) {
                  if(stabp->archive != 'a') { /* not arch for append */
                    if(reqmegs) {       /* DSDS must have assigned storage */
                      call_archive(stabp, hnext, list, remerrno);
                    }
                  }
                  else {                        /* arch for append */
                    stabp->archive_complete = 1;/* don't del wd at end */
                  }
                }
              }
              /* New 03Dec98. Rm from any archive_group so the wd
               * will get deleted at the end (by deregdsds())
              */
              if(abortnum == 1) {       /* rm from any archive_group */
                stabp->archive_group = 0;
              }
              abortit(-1);              /* ret w/module fatal error status */
            }
          }
        }
      }
      printk("    %s %s %x complete:\n    %s\n",
                remhost,stabp->name,remtid,remoutfile);
      if((warnmsg=getkey_str(list, "WARNING"))) {
        printk("    **with warning msg:\n    %s\n", warnmsg);
      }

      pvm_kill(hnext->tid);             /* kill this completed task */
      /*peklog("In process_resp(): tid = %x %s\n", hnext->tid, stabp->name);*/
      hnext->tid = 0;                   /* don't kill again in kill_pvm() */

      if(!stabp->busyall) {
        /* Put the following msg back in when implement split
        /*printk("    All hosts for %s have completed\n", stabp->name);
        */
        /* archive output if flaged and not in the middle of an arch group */
        if(stabp->archive && (stabp->archive_group != 1)) {
          if(stabp->archive != 'a') {   /* not arch for append */
            if(reqmegs) {               /* DSDS must have assigned storage */
              call_archive(stabp, hnext, list, remerrno);
            }
            else {
              printk("**Ignoring archive of data for %s:\n",stabp->name);
              printk("  No DSDS storage alloc request in map file.\n");
            }
          }
          else {                        /* arch for append */
            stabp->archive_complete = 1;/* don't del wd at end */
          }
        }
/* !!TEMP spot for new code?? */
        /* catalog any level #0 ds if storage was obtained */
        if(!stabp->archive && reqmegs)
          call_archive_0(stabp, hnext, list, remerrno);
      }
      if(abortnum) {
        pelog("%x\t%d\tabort\t%d\t%s\t%s\t%s\n",
        remtid,remerrno,stabp->archive_complete,stabp->name,remoutfile,auxinfo);
      }
      else {
        pelog("%x\t%d\tnormal\t%d\t%s\t%s\t%s\n",
        remtid,remerrno,stabp->archive_complete,stabp->name,remoutfile,auxinfo);
      }
      break;                            /* exit while(hnext) */
    }
    hnext = (HDATA *)gethnext((HDATA *)MONE);/* get next hdata entry */
  }                                     /* end while(hnext) */
  if(!hnext) {
    pemail("***No such host %s for server %s\n",remhost,stabp->name);
    abortit(1);
  }
  /* wait for all the hosts of a server or all the servers of a group */
  wait = 0;
  if(groupid = stabp->groupid) {        /* this is a group */
    for(j=0; j<MAX_SERV; j++) {
      if(stab[j].groupid == groupid) {/* same group as completing server */
        if(stab[j].busyall) {
           wait = 1;                    /* someone in group still busy */
          break;
        }
      }
    }
  }
  else {                                /* not a group but host may be busy*/
    if(stabp->busyall) wait = 1;
  }
  /* eventually some call to process_resp() below will fall through this wait */
  while(wait) {
    if((sindex=resp_any())== -1)        /* get anyones response */
      who_died();                       /* timeout. see if anyone died */
    else {
      process_resp(sindex);             /* re-enter to process the response */
      freekeylist(&list);
      return;
    }
  }
  if(groupid)
    printk("End Server Group #%d\n", groupid);

  /* kick off the next server type that's not in the completing group */
  stabn = &stab[ix+1];                  /* next server table */
  while(1) {
  if(stabn->name != NULL && groupid != 0 && stabn->groupid == groupid) stabn++;
  else break;
  }
  if(stabn->name != NULL) {             /* kick off the next server type */
    if(groupid = stabn->groupid) {      /* this is a group */
      printk("Start Server Group #%d\n", groupid);
      while(stabn->groupid == groupid) {
        kick_server(stabn);     /* create the keylist & start the server */
        stabn++;
      }
    }
    else
      kick_server(stabn);               /* kick the non-group server */
  }
  freekeylist(&list);
}

//!!TDB at end of all module execution

/* Deregister with dsds_svc as appropriate. If any output storage was gotten
 * from dsds_svc then any unused amount is returned. Any data sets not archived
 * with dsds_svc are deleted (except for level number 0). 
 * Reports on the storage used by each output ds.
 * Called from main() or abortit().
*/
void dereg()
{
  KEY *alist, *blist;
  SERVER *sptr;
  HDATA *hdata;
  DRMS_Record_t *rsx;
  argument *arg;
  char *wd, *cptr;
  char oname[MAX_STR], mcmd[256];
  double free_bytes, dsize;
  double mapstore = 0.0;		/* total storage for pe map run */
  int i, j, outnsets;

  if(reqmegs) {				/* output storage was obtained */
    /* first remove any directories not archived */
    for(i=0; i < MAX_SERV; i++) {
      sptr = &stab[i];
      if(sptr->name == NULL)
        break;
      /* see if output was really archived */
      arg = sptr->arguments;
      hdata = (HDATA *)gethnext(sptr->hosts);/* get first hdata entry */
      while(arg->kind != ARG_END) {	/* find all output datasets */
        if(arg->kind == ARG_DATA_OUT) {
          sprintf(oname, "%s_nsets", arg->key);
          outnsets = getkey_int(hdata->param_list, oname);
          for(j=0; j < outnsets; j++) {/* do for each dataset */
            sprintf(oname, "%s_%d_wd", arg->key, j);
            if(wd=getkey_str(hdata->param_list, oname)) {
              /* only count if a /SUM assigned dir */
              cptr = strstr(wd, "/SUM");
              if(cptr == wd) {		/* a /SUM wd */
                /* only count the last one in an archive group */
                /* and never count any appendable ds */
                if((sptr->archive_group != 1) && (sptr->archive != 'a')) {
                  dsize = du_dir(wd); /* get # of bytes */
                  mapstore += dsize;
                }
              }
            }
            /* don't delete any level # 0 */
            sprintf(oname, "%s_%d_level_sn", arg->key, j);
            if(findkey(hdata->param_list, oname)) {
              if(getkey_int(hdata->param_list, oname) == 0)
                continue;		/* do the next for() */
            }
            if(sptr->archive != 'a') { /* dont rm any appendable */
              if(!sptr->archive_complete && (sptr->archive_group != 1)) {
                //tell drms not to archive
                sprintf(oname, "%s_%d_rs", arg->key, j);
                rsx = (DRMS_Record_t *)getkey_fileptr(sptr->map_list, oname);
                rsx->seriesinfo->archive = 0;
                if(wd) {
                  if(strstr(wd, "/SUM")) {	/* dont rm a non SUM dir */
                    printk("Removing non-archived wd=%s bytes=%g\n", wd, dsize);
                    if(rmdirs(wd, dsdswd)) { /* rm dir and any other blanks */
                      printk("**Can't rm wd=%s\n", wd);
                      printk("(NOTE: May be due to an NFS delay. Usually ignore.)\n");
                    }
                  }
                }
              }
            }
          }
        }
        arg++;
      }
    }
    printk("High water storage bytes=%g\n", mapstore);
//    alist=newkeylist();			/* create an empty list */
//    /*free_bytes = reqbytes - (bytes_used * 1.07); /* add 7% for overhead */
//    //free_bytes = reqbytes - bytes_used;
//    setkey_double(&alist, "free_bytes", free_bytes);/* # of bytes not used */
//    setkey_double(&alist, "req_bytes", reqbytes);/*# of bytes orig requested*/
//    setkey_long(&alist, "dsds_uid", uid);	/* give our OPEN uid */
//!!!TBD clean this up
//    if((blist = (KEY *)call_dsds(&alist, REQDEALLOC, dsds_tid, pe_tid, msg, debugflg)) == NULL) {
//     pemail("***Can't deallocate %g bytes from dsds_svc\n",free_bytes);
//    abortit(1);
// }
//    freekeylist(&alist); freekeylist(&blist);
//    printk("Deallocated %g bytes with dsds_svc\n", free_bytes);
//    if(mapstore > reqbytes) {
//      printk("**High water storage > storage requested. Map file should request more\n");
      /* send out special mail for this */
//      sprintf(mcmd, "echo \"Storage exceeded for user=%s map=%s\nRequested=%g High Water=%g\" | Mail -s \"pe excess storage\" jim,jeneen,%s", username, pe_map, reqbytes, mapstore, username);
//      system(mcmd);
//    }
  }
//  if(reqmegs || dsdsin) {
//    alist=newkeylist();                  /* create an empty list */
//    setkey_long(&alist, "dsds_uid", uid);/* give our OPEN uid */
//    if((blist = (KEY *)call_dsds(&alist, REQCLOSE, dsds_tid, pe_tid, msg, debugflg)) == NULL) {
//      msg("***Can't CLOSE with dsds_svc\n");
//      abortit(1);
//    }
//    freekeylist(&alist); freekeylist(&blist);
//    printk("Deregistered with dsds_svc as uid=%ld\n", uid);
//    uid = 0;				/* indicate that we're closed */
//  }
}

/* This is the do_pipe to run the strategy level modules.
 * It will send the keylist of the first server type to the all of the first
 * servers and await any response. This will start off the chain where each
 * response causes the next free server type to be called.
 * This only returns when the entire pipe is done.
*/
void do_pipe()
{
  SERVER *sptr;
  int i, groupid;
  int sindex;

  /* Start the first server or servers if this is a group. All members of */
  /* a group have the same non-0 groupid as set by get_cmd(). */
  sptr = &stab[0];
  if(groupid = sptr->groupid) {
    printk("Start Server Group #%d\n", groupid);
    while(sptr->groupid == groupid) {
      kick_server(sptr);        /* create the keylist & start the server */
      sptr++;
    }
  }
  else
    kick_server(sptr);          /* kick the non-group server */
  /* Now wait for anyone still busy */
  for(i=0; i<MAX_SERV; i++) {
    if(stab[i].name == NULL)
      break;                            /* no more servers */
    while(stab[i].busyall) {
      if((sindex=resp_any())== -1)      /* get anyones response */
        who_died();                     /* timeout. see if anyone died */
      else
        process_resp(sindex);           /* clears the busy flags */
    }
  }
}


// Initial setup stuff called when main is first entered.
void setup(int argc, char *argv[])
{
  struct timeval tvalr;
  struct tm *t_ptr;
  char logname[128], string[128], cwdbuf[128], idstr[256];
  int tvalr_int, i;

  if((pe_tid=start_pvm(printf)) == 0) {
    fprintf(stderr, "Can't start a pvm daemon!!\n");
    exit(1);
  }
  if(!prod_host_set()) {       /* set primary & second host names */
    fprintf(stderr, "Error accessing file %s\n", PHNAME);
    exit (1);
  }
  if (signal(SIGINT, SIG_IGN) != SIG_IGN)
      signal(SIGINT, sighandler);
  if (signal(SIGTERM, SIG_IGN) != SIG_IGN)
      signal(SIGTERM, sighandler);

  gettimeofday(&tvalr, NULL);
  #ifndef __sparc
  tvalr_int = (int)tvalr.tv_sec; /* need int vrbl for this to work on sgi4*/
  t_ptr = localtime(&tvalr_int);
  #else
  t_ptr = localtime(&tvalr.tv_sec);
  #endif
  sprintf(datestr, "%s", asctime(t_ptr));
  datestr[24] = NULL;
  wd_dup_list = newkeylist();	/* init the list of potential dup wd's */
  wd_mk_list = newkeylist();	/* init the list of potential mkdir wd's */
  JPElist = newkeylist();	/* init the list for any archive group */
  tae_tid=pvm_parent();		/* get tid if tae spawned us */
  uid = (uint)getpid();		/* use our pid for the old dsds_uid */
  pvm_serror(1);		/* enable auto error reporting */
  if(!(username = (char *)getenv("USER"))) username = "nouser";
  sprintf(logname, PELOGFILE, username, getpid());
  open_pelog(logname);
  printk_set(printf, printf);
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
  for(i=0; i < argc; i++) { 	/* echo cmd line */
    sprintf(string, "%s%s", argv[i], (i < argc-1) ? " " : "");
    strcat(idstr, string);
  }
  strcat(idstr, "\n");
  sprintf(string, "pe started as tid=%x pid=%d user=%s\n", 
			pe_tid, getpid(), username);
  strcat(idstr, string);
  sprintf(mailname, PEMAILFILE, username, getpid());
  open_pemail(mailname, idstr);
  printk("%s", idstr);
  pelog("#%s for %s on %s\n", logname, username, datestr);
  pelog("#listed in order of completion\n");
  pelog("#tid\tstatus\tcomplet\tarchive\tmodule\t\t\"out\"returned\tauxinfo\n");
  pelog("TID\tSTATUS\tCOMPLET\tARCHIVE\tMODULE\tOUT\tAUXINFO\n");
  pelog("---\t------\t-------\t-------\t------\t---\t-------\n");
  /* note: the columns in this file are defined in dohrsc.h */
  umask(002);			/* allow group write */
}

// Module main function. 
int DoIt()
{
  pid_t pid;
  char *args[6];
  char callcmd[128];

  cmdparams_get_argv(&cmdparams, &argv, &argc);
  setup(argc, argv);
  get_cmd(argc, argv);
  spawn_pvm();
  if(setjmp(env) != 0) {	//longjmp() has been called. get out
    dereg();			//end session !!TBD need?
    kill_pvm();                   /* kill all servers on all hosts */
    pvm_exit();
    if (pemailfp) fclose(pemailfp);
    mailit();
    if (pelogfp) fclose(pelogfp);
    printk("PE Abnormal Completion\n");
    return(1);
  }
  do_pipe();
  dereg();			//end session !!TBD need?
  kill_pvm();                   /* kill all servers on all hosts */
  pvm_exit();
  if (pemailfp) fclose(pemailfp);
  mailit();
  if (pelogfp) fclose(pelogfp);
  printk("PE Normal Completion\n");
  return(0);
}
