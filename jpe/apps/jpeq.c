/*-----------------------------------------------------------------------------
 * /home/production/cvs/JSOC/proj/jpe/apps/jpeq.c
 *-----------------------------------------------------------------------------
 *
 * This is the old peq.c from SOI MDI that is converted to run like the
 * original peq but with DRMS/SUMS interface instead of DSDS.
 * The original is found in /CM/src/pipe/peq.c
 *
 * The following original peq options are deprecated in jpeq:
 *  -d, -t, -L, -O , -S
 ***************************************************************************/

/* Here are the original peq.c notes: */
/* peq.c
 * Queries the given database for the given dataset name and prints 
 * the info to stdout.
 * Note: Does not have the capability to ask for a specific svc_version 
 * or svc_name as in a pe map file.
 * Sample call:
 *   peq [-v] [-w|W] [-d] [-A] [-L] [-Oname] [-t#] [-f] [-S] dataset_name
 * where: -v = verbose
 *        -w = print out the working dirs only
 *        -W = print out the working dirs and creation dates only
 *        -d = run dsds_svc in debug mode
 *        -Ahost = access the ampex_svc if need to retrieve data
 *             Host is an optional name to assign local /PDS storage from.
 *             If not specified, will use the env vrbl PDS_SET_HOST, else
 *             will use the local host.
 *        -L = run local dsds_svc/ampex_svc. don't use pe_rpc
 *             This is for development use only and requires a password
 *        -Oname = Oracle db to query instead of the default char *dbname
 *        -t# = (touch) #of days to retain datasets before deletable or
 *		"keep" to retain until findkeep is run. 
 *        -f = the dataset_name is a file name containing dataset_names
 *	       like "prog:name,level:name[#],series:name[#]"
 *             This file is processed 50 lines at a time.
 *        -S = get the dataset from the SUMS
 *
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
//#include <signal.h>
#include <soi_key.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <rpc/rpc.h>
//#include <rpc/pmap_clnt.h>
//#include <sys/socket.h>
#include <unistd.h>     //for alarm(2) among other things...
#include <printk.h>
#include <pvm3.h>
#include <setjmp.h>
//#include <soi_error.h>
//#include <sum_pe.h>
#include "soi_args.h"
#include "pe.h"
#include "dsds.h"

#define MAXLINE 96		/* max line size for passwd */
#define MAXDSREQ 200		/* max ds in any one datacollection */
#define MAXFLINES 100		/* max lines in a -f file */
#define NOTSPECIFIED "***NOTSPECIFIED***"

extern void pepeqprog_1(struct svc_req *rqstp, SVCXPRT *transp);
extern int pepeq_wait();
extern int pepeq_poll();
extern KEY *call_drms_in(KEY *list, int dbflg);
void printkey();
void deregdsds();
int resp_dsds(int dsdstid);
int soi_errno = NO_ERROR;

long uid;			/* assigned by dsds_svc REQOPEN call */
int verbose;			/* set by get_cmd() */
int wdonly;			/* set by get_cmd() */
int wddate;			/* set by get_cmd() */
int taeflg;			/* set by get_cmd() */
int dsds_tid;			/* originally the tid of dsds_svc */
                                /* now the tid of pe_rpc that gets us to dsds*/
int pe_tid;			/* tid of this pe */
int ampex_tid;			/* tid of ampex_svc */
int tae_tid;			/* tid if we were spawned by tae ui */
int abort_active;		/* set while doing an abort */
int retrieveflg;		/* retrieve data from tape if not on-line */
int ampexflg;			/* use ampex_svc instead of lago_svc */
int nocontrolc;                 /* don't honor ^C if set */
int fileflg;			/* take ds names from the given file */
int debugflg;			/* run all pvm servers in the debug mode */
				/* also do keyiterate. Don't use w/tae */ 
int touchflg = -1;		/* #of days to retain ds before deletable */
int ccactive = 0;		/* ^C already in progress */
int dbxflg = 0;
int linenum = 0;
int ds_lines = 0;
uint32_t ouruid;

char database[MAX_STR];
char pds_set_host[MAX_STR];
char thishost[MAX_STR];
char line[1024];
char *dbname = "mdi_2";
char *pdshost;
char *username;
char *in_ds;
FILE *dsfp;
KEY *dslist;
KEY *alist;
CLIENT *clntsumpe;
SUM_t *sumhandle = NULL;

static struct timeval TIMEOUT = { 20, 0 };

// List of default parameter values.
// NOTE: only used for cmdparams -H call. Normally use get_cmd()
ModuleArgs_t module_args[] = {
  {ARG_STRING, "dsname", NOTSPECIFIED, "dataset name"},
  {ARG_FLAG, "t","0","#of days to retain input datasets before deletable"},
  {ARG_FLAG, "A", "0", "retrieve from tape flag"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "w", "0", "wd only flag"},
  {ARG_FLAG, "f", "0", "dsname is a file of dataset names"},
  {ARG_END}
};

CmdParams_t cmdparams;
char **argv = NULL;
int argc = 0;
char *module_name = "jpeq";      // Module name presented to DRMS
jmp_buf env;


/* Return ptr to "mmm dd hh:mm:ss". */
char *datestring(void)
{
  time_t t;
  char *str;

  t = time(NULL);
  str = ctime(&t);
  str[19] = 0;
  return str+4;          /* isolate the mmm dd hh:mm:ss */
}

/* Got a fatal error sometime after registering with pvm. 
 * Degregister and close with dsds_svc as approriate. 
*/
void abortit()
{
  if(abort_active) {		/* we're already aborting */
    printk("Abort an abort! dsds_svc may not have cleaned up\n");
    longjmp(env, 2);      //get out of here
  }
  abort_active=1;
  printk("Abort in progress ...\n");
  longjmp(env, 2);      //get out of here
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
  if(nocontrolc) {
    /* tape retrieval underway. send cancel to dsds_svc. */
    if(ccactive) {
      printk("A ^C is already active. I'm awaiting a dsds_svc reply...\n");
    }
    else {
      ccactive = 1;
      setkey_long(&list, "dsds_uid", uid); /* give our OPEN uid */
      setkey_int(&list, "ampex_tid", ampex_tid);
      setkey_str(&list, "USER", username);
      printk("^C signal received by peq\n");
/*************************!!!TBD****************************************
      printk("Sending TAPECANCEL request to dsds_svc...\n");
      if((list = call_dsds(&list,REQTAPECANCEL,dsds_tid,pe_tid,msg,debugflg)) == NULL) {
        printk("Can't REQTAPECANCEL with dsds_svc\n");  // don't abort here
      }
      else {
        if(getkeytype(list,"MSG_PEND")==KEYTYP_INT) { // wait for completion
          respcnt = 0;
          while(1) {
            if(resp_dsds(dsds_tid)) {	// rets every timeout interval
              printk("..."); respcnt++;
              if(respcnt == 20) {
                respcnt = 0;
                printk("\n");
              }
              continue;
            }
            if(respcnt != 0) printk("\n");
            break;
          }
          printk("Ampex tape read request has been canceled.\n");
          freekeylist(&list);
          abortit();
        }
        else {
          printk("Failed to get TAPECANCEL msg to ampex_svc\n");// don't abort
        }
      }
      freekeylist(&list);
*****************************************************************/
    }
  }
  else {
    printk("\npeq received a termination signal\n");
    abortit();
  }
}

void usage()
{
  printf("Usage:\n");
  printf(" peq [-v] [-w|W] [-A] [-t#] [-f] dataset_name\n");
  printf(" where -v = verbose mode\n");
  printf("       -w = print out the working dirs only\n");
  //printf("       -W = print out the working dirs and creation dates only\n");
  printf("       -A = access the ampex_svc if need to retrieve data\n");
  printf("       -t# = (touch) #of days to retain datasets before deletable\n");
  printf("       -f = the dataset_name is a file name containing dataset names\n");
  printf("        dataset_name like \"prog:name,level:name,series:name[#]\"\n");
  printf("        This file is processed 50 lines at a time.\n");
  printf("        The series# can be a range. A level# may be given.\n");
  printf(" Note: If you use \"[]\" in the ds name you must put the ds name in quotes\n");
  exit(1);
}

/* Get the next 50 lines from the input file containing the dataset names.
 * Will open the file on the first call.
 * Returns the number of lines actually read, else 0 and closes the file.
 * Returns the datacollection of the lines read in malloc'd memory in_ds.
 * Uses the global vrbl dsfp & linenum for the position in the input file.
 * We can only do 50 lines at a time due to error in rpc lib when
 * we try to send too long a keylist through pe_rpc/pe_rpc_svc. (OLD from peq)
*/
int get_ds_lines()
{
  int i;
  int first = 1;

  if(!fileflg) { return(0); }		/* not in -f file mode */
  if((in_ds = (char *)malloc(65536)) == NULL) {
    printk("Can't malloc storage for dataset names\n");
    longjmp(env, 2);      //get out of here
  }
  for(i = 0; i < MAXFLINES; i++) {
    if(fgets(line, 1024, dsfp)) {    /* get a ds name */
      if(line[0] == '#' || line[0] == '\n') continue;
      line[strlen(line)-1] = NULL;      /* elim cr at end */
      if(first) { first = 0; strcpy(in_ds, line); }
      else { strcat(in_ds, ";"); strcat(in_ds, line); }
    }
    else break;
  }
  if(i == 0) {
    fclose(dsfp);
    return(0);
  }
  return(i);
}

/* Gets the command line and parses the dataset name into  the global
 * dslist keylist. Called before we register with dsds.
*/
void get_cmd(int argc, char *argv[])
{
  int c, getpasswd;
  char *envpasswd, *cptr;
  char hostname[MAX_STR];

  gethostname(hostname, MAX_STR);
  while(--argc > 0 && (*++argv)[0] == '-') {
    while(c = *++argv[0])
      switch(c) {
      case 't':
        if(*++argv[0] != NULL) {	/* get # of days */
          cptr = argv[0];
          if(!strcmp(cptr, "keep")) {	/* it's a touch "keep" */
            touchflg = 9999999;		/* pass on the default keep value */
          }
          else {
            if(!isdigit((int)*cptr)) {
              printk("-t switch must say \"keep\" or give integer #of days\n");
              exit(1);
            }
            touchflg = atoi(argv[0]);
          }
        }
        while(*++argv[0] != NULL);
        --argv[0];
        break;
      case 'v':
        verbose=1;
        break;
      case 'w':
        wdonly=1;
        break;
      case 'W':
        wdonly=1; wddate = 1;
        break;
      case 'A':
        ampexflg=1;
        retrieveflg=1;
        break;
      case 'f':
        fileflg=1;
        break;
      default:
        usage();
        break;
      }
  }
  if(argc != 1) usage();
  dslist = newkeylist();
  if(!fileflg) {			/* ds are on cmd line */
    setkey_str(&dslist, "in", argv[0]);
    if (parse_list (&dslist, "in")) {
      printk("Error in parse_list() of dataset name\n");
      longjmp(env, 2);      //get out of here
    }
    c = getkey_int(dslist, "in_nsets");
    if(c > MAXDSREQ) {
      printk("Error: max number of ds in any one request is %d\n", MAXDSREQ);
      longjmp(env, 2);      //get out of here
    }
  }
  else {				/* get ds from a file */
    /* open and get the first 100 lines */
    /* we can only do 100 lines at a time due to error in rpc lib when */
    /* we try to send too long a keylist through pe_rpc/pe_rpc_svc. (OLD peq) */
    if(!(dsfp=fopen(argv[0], "r"))) {
      printk("Can't open the input file %s\n", argv[0]);
      longjmp(env, 2);      //get out of here
    }
  }
}

/* Determine the pvm machine configuration and which tasks are running. 
 * Spawn the dsds_svc on the local machine or spawn the pe_rpc that will
 * bridge us to the master dsds_svc via the pe_rpc_svc.
*/
void spawn_pvm()
{
  struct hostinfo *hostp;
  int info, nhost, narch, i;

  /* get info on the present virtual machine configuration */
  if((info=pvm_config(&nhost, &narch, &hostp))) {
    printk("Can't get pvm configuration (%d). Pvm daemon(s) running?\n", info);
    abortit();
  }
  for(i=0; i<nhost; i++) {
    if(!strcmp(hostp[i].hi_name, thishost)) {
      break;
    }
  }
  if(i == nhost) {
    printk("No pvm daemon running on %s\n", thishost);
    abortit();
  }

  //if(spawn_pe_rpc(dbname, pe_tid, &dsds_tid, debugflg, printk))
  //  abortit();
}

/*!!!!!!!!!!OBSOLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/* Take all keys on the given keylist that have a number in their name and
 * change the name to add the bias number. For example, for bias = 50
 * a key named in_1_wd will become in_51_wd in the keylist.
*/
KEY *renumber_keylist(KEY *list, int bias)
{
  KEY *walker = list;
  char *cptr, *cptr1, *wptr;
  char tmp[256], newkey[256], newnum[32];
  int i, m, n;

  if (!list) {
    return NULL;
  }
  while (walker) {
    if(!(cptr = (char *)index(walker->name, 'in_'))) {
      if(!(cptr = (char *)rindex(walker->name, '_'))) {
        walker = walker->next;
        continue;
      }
      strcpy(tmp, cptr+1);
    }
    else {
      wptr = (char *)walker->name;
      m = (int)(cptr - wptr)+1;
      if(cptr1 = (char *)index(cptr+1, '_')) {
        n = (int)(cptr1-(cptr+1));
        strncpy(tmp, cptr+1, n);
        tmp[n] = NULL;
      }
      else {
        strcpy(tmp, cptr+1);
      }
    }
    if(isdigit((unsigned char)tmp[0])) {
      /*printf("name = %s tmp = %s\n", wptr, tmp); /* !!TEMP */
      i = atoi(tmp);
      i += bias;
      strncpy((char *)newkey, wptr, m);
      newkey[m] = NULL;
      sprintf(newnum, "%d", i);
      strcat((char *)newkey+m, newnum);
      if(cptr1) {
        strcat((char *)newkey+m+strlen(newnum), cptr1);
      }
      cptr = (char *)malloc(strlen(newkey)+1);
      strcpy(cptr, newkey);
      free(walker->name);
      walker->name = cptr;
      /*printf("newkey = %s\n", newkey);*/
      if(strstr(cptr, "inname")) {      /* must change value of these keys */
        sprintf(tmp, "in_%d", i);
        setkey_str(&walker, cptr, tmp);
      }
    }
    walker = walker->next;
  }
  return walker;
}

/* Query the dataset name given on the command line (now in dslist).
 * Abort on any error.
 * 09Jan97 - Modified for REQREQALL but make sure keep backward compatible
 * as peq is used in many scripts.
*/
void queryds()
{
  char keyname[MAX_STR], ext[MAX_STR], server_name[MAX_STR], dummy[MAX_STR];
  char *wd, *cd, *instr, *argname, *cptr, *cptr1, *server_name_p, *call_err;
  int respcnt, status, i;
  int offline = 0;
  uint32_t sumpeback;
  register SVCXPRT *transp;
  KEY *blist = newkeylist();

  alist = newkeylist();
  /* give name of the arg for future REQREQALL call */
  setkey_str(&dslist, "arg_data_in_0", "in");
  setkey_long(&dslist, "dsds_uid", uid);/* give our OPEN uid */
  setkey_str(&dslist, "pds_host", pdshost); /* host for the /PDS set */
  /* first determine if off-line and need to retrieve */
  if(ampexflg)
    setkey_int(&dslist, "ampex_tid", ampex_tid);
  if(touchflg != -1)
    setkey_int(&dslist, "touch", touchflg);

  if(dbxflg) {
    printk("\nThe dslist at the REQREQALLNR:\n");
    keyiterate(printkey, dslist);
  }
  printk("Querying for the input datasets...\n");
  nocontrolc = 1;                       /* no ^C during retrieve */
  alist = (KEY *)call_drms_in(dslist, dbxflg);
  if(alist == NULL) {
    printk("Can't resolve/retrieve the input datasets with drms.\n");
    printk("(Usually no -A switch specified.)\n");
    abortit(1);
  }

  //renumber_keylist(alist, linenum);   //!!!TBD what to do with this


//!!TEMP
    //printk("\nThe alist from call_drms_in():\n");
    //keyiterate(printkey, alist);
    printk("\n");
//return; //!!TEMP

    for(i = linenum; ; i++) {
      sprintf(ext, "inname_%d", i);
      if(!findkey(alist, ext)) break;     /* done, exit for() loop */
      argname = GETKEY_str(alist, ext);   /* gets e.g. in_0 */
      if(wdonly) {
        sprintf(ext, "%s_wd", argname);   /* e.g. in_0_wd */
        wd = GETKEY_str(alist, ext);
        printk("%s\tKEYTYP_STRING\t%s\n", ext, wd);
        if(wddate) {			/* and creat_date too */
          sprintf(ext, "%s_creat_date", argname);
          cd = GETKEY_str(alist, ext);
          printk("%s\tKEYTYP_STRING\t%s\n", ext, cd);
        }
      }
    }
    if(!wdonly) keyiterate(printkey, alist);
    freekeylist(&alist);
}

/* Deregister with dsds_svc. Called from main() or abortit().
*/
void deregdsds()
{
  KEY *list, *retlist;

  list=newkeylist();                  /* create an empty list */
  setkey_long(&list, "dsds_uid", uid);/* give our OPEN uid */
  if((retlist = call_dsds(&list,REQCLOSE,dsds_tid,pe_tid,printk,debugflg)) == NULL) {
    printk("Can't CLOSE with dsds_svc\n");
    freekeylist(&list);
    abortit();
  }
  printk("Deregistered with dsds_svc as uid=%ld\n", uid);
  uid = 0;				/* indicate that we're closed */
  freekeylist(&list); freekeylist(&retlist);
  freekeylist(&dslist);
}

/* Initial setup stuff called when main if first entered.
*/
void setup()
{
  if((pe_tid=start_pvm(printk)) == 0) {
    fprintf(stderr, "Can't start a pvm daemon!!\n");
    exit(1);
  }
  printk_set(printf, printf);
  if(!(username = (char *)getenv("USER"))) username = "nouser";
  if (signal(SIGINT, SIG_IGN) != SIG_IGN)
      signal(SIGINT, sighandler);
  if (signal(SIGTERM, SIG_IGN) != SIG_IGN)
      signal(SIGTERM, sighandler);
  /*tae_tid=pvm_parent();		/* get tid if tae spawned us */
  pvm_serror(1);			/* enable auto error reporting */
  gethostname(thishost,MAX_STR);
}

int DoIt()
{
  int c;

 if(setjmp(env) != 0) {        //longjmp() has been called. get out
    pvm_exit();
    //if (pemailfp) fclose(pemailfp);
    //mailit();
    //if (pelogfp) fclose(pelogfp);
    if(sumhandle) 
      SUM_close(sumhandle, printk);
    printk("jpeq Abnormal Completion\n");
    return(1);
  }


  cmdparams_get_argv(&cmdparams, &argv, &argc);

  setup();				/* start pvm and init things */
  get_cmd(argc, argv);			/* check the calling sequence */
  spawn_pvm();				/* start the servers as needed */
  if(!wdonly) {				//open sums so can call sum_info()
    if((sumhandle = SUM_open(NULL, NULL, printk)) == 0) {
      printk("Failed on SUM_open()\n");
      longjmp(env, 2);      //get out of here
    }
  }
  while(1) {
    if(!fileflg) {
      queryds();
      break;
    }
    if(!(ds_lines = get_ds_lines())) {
      break;                              /* break while(1) */
    }
    dslist = newkeylist();
    setkey_str(&dslist, "in", in_ds);
    if (parse_list (&dslist, "in")) {
      printk("Error in parse_list() of dataset name\n");
      longjmp(env, 2);      //get out of here
    }
    c = getkey_int(dslist, "in_nsets");
    if(c > MAXDSREQ) {
      printk("Error: the next list of ds exceed the max of %d\n", MAXDSREQ);
      longjmp(env, 2);      //get out of here
    }
    queryds();				/* query for the dataset name */
    if(ds_lines >= MAXFLINES) {
      printk("Error: max of %d lines in -f file mode\n", MAXFLINES);
      longjmp(env, 2);      //get out of here
    }
    linenum += ds_lines;
  }
  //deregdsds();			/* end the session with dsds_svc */
  pvm_exit();
  //freekeylist(&dslist);
  if(sumhandle) 
    SUM_close(sumhandle, printk);
  printk("jpeq Normal Completion\n");
  return(0);
}
