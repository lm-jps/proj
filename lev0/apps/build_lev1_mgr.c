/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev1/apps/build_lev1_mgr.c
 *-----------------------------------------------------------------------------
 *
 * This is a module that runs with DRMS and processes lev0
 * filtergrams to lev1 by running the build_lev1 module. 
 *
 *build_lev1_mgr
 *	instru= hmi or aia
 *	dsin= default hmi.lev0e or aia.lev0e
 *	dsout= hmi.lev1e or aia.lev1e
 *	brec= first lev0 record# to process
 *	erec= last lev0 record# to process
 *
 * Runs build_lev1 processes to create lev1 datasets
 * Has two modes:
 * Stream Mode (one instance)
 *	brec=0, erec=0 
 *	- start at lev0 rec of highest lev1 FSN. Run forever.
 *	- fork up to 8 (MAXCPULEV1) build_lev1 for every 
 *	  12 (MAXRECLEV1) lev0 records. 
 *	- when an build_lev1 completes, fork another for next 12 rec
 *	- if no more lev0 records available, sleep and try again
 *	- if 8 CPU not enough to keep up with lev0, go to 16, etc.
 *
 * Reprocessing Mode (any number of instances)
 *	brec=1000, erec=2000 
 *	- qsub up to 128 (MAXQSUBLEV1) build_lev1 for 12 records ea
 *	- when a job completes qsub next 12 records until erec is reached
 *	- when all jobs are done, build_lev1_mgr will exit
 *
 */ 

#include <jsoc.h>
#include <cmdparams.h>
#include <drms.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>	//for umask(2)
#include <unistd.h>	//for alarm(2) among other things...
#include <printk.h>
#include <errno.h>
#include <sys/wait.h>


//default in and out data series
#define LEV0SERIESNAMEHMI "hmi.lev0e"
#define LEV0SERIESNAMEAIA "aia.lev0e"
#define LEV1SERIESNAMEHMI "su_production.hmi_lev1e"	//temp test case
#define LEV1SERIESNAMEAIA "su_production.aia_lev1e"	//temp test case

#define LEV1LOG_BASEDIR "/usr/local/logs/lev1"
#define QSUBDIR "/scr21/production/qsub"
#define DONUMRECS 120		//!!TEMP #of lev0 records to do and then exit
#define H1LOGFILE "/usr/local/logs/lev1/build_lev1_mgr.%s.log"
#define NUMTIMERS 8		//number of seperate timers avail
#define MAXRECLEV1 12	//number of lev0 to lev1 images to do at a time
#define MAXCPULEV1 8	//max number of forks to do at a time for stream mode
#define NOTSPECIFIED "***NOTSPECIFIED***"

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "instru", NOTSPECIFIED, "instrument. either hmi or aia"},
  {ARG_STRING, "dsin", NOTSPECIFIED, "dataset of lev0 filtergrams"},
  {ARG_STRING, "dsout", NOTSPECIFIED, "dataset of lev1 output"},
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_INTS, "brec", "0", "first lev0 rec# to process"},
  {ARG_INTS, "erec", "0", "last lev0 rec# to process"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_END}
};

ModuleArgs_t *gModArgs = module_args;
CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "build_lev1_mgr";

FILE *h1logfp;		// fp for h1 ouput log for this run 
FILE *qsubfp;		// fp for qsub script
static TIME sdo_epoch;
static char datestr[32];
static char open_dsname[256];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];

int verbose;
long long brec, erec;		//begin and end lev0 rec# to do
int stream_mode = 0;		//0=qsub build_lev1, 1=fork it locally
int hmiaiaflg = 0;		//0=hmi, 1=aia
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when build_lev0 is called for a restart
int abortflg = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
char logname[128];
char argdsin[128], argdsout[128], arglogfile[128], arginstru[80];
char argbrec[80], argerec[80];
char timetag[32];
char tlmseriesname[128];	// e.g. hmi.tlm
char lev0seriesname[128];	// e.g. hmi.lev0
char *username;			// from getenv("USER") 
char *logfile;			// optional log name passed in 
char *instru;			// instument. hmi or aia 
char *dsin;			// lev0 input dataset
char *dsout;			// lev1 output dataset


int nice_intro ()
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage: build_lev1_mgr [-vhr]\n"
	"instru=<hmi|aia> dsin=<lev0> dsout=<lev1> brec=<rec#> erec=<rec#>"
	"\n[logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
	"  -r: restart. only used when we restart our selves periodically\n"
	"instru= instrument. must be 'hmi' or 'aia'\n"
	"dsin= data set name of lev0 input\n"
	"      default hmi=hmi.lev0e   aia=aia.lev0e\n"
	"dsout= data set name of lev1 output\n"
	"      default hmi=su_production.hmi_lev1e   aia=su_production.aia_lev1e\n"
	"brec= first lev0 rec# to process. 0=Stream Mode if erec also 0\n"
	"erec= last lev0 rec# to process. 0=Stream Mode if brec also 0\n"
	"logfile= optional log file name. If not given uses:\n"
        "         /usr/local/logs/lev1/build_lev1_mgr.<time_stamp>.log\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  restartflg = cmdparams_get_int (&cmdparams, "r", NULL);
  return (0);
}

/* Return pointer to "Wed Jun 30 21:49:08 1993\n" */
char *get_datetime()
{
  struct timeval tvalr;
  struct tm *t_ptr;
  static char datestr[32];

  gettimeofday(&tvalr, NULL);
  t_ptr = localtime((const time_t *)&tvalr.tv_sec);
  sprintf(datestr, "%s", asctime(t_ptr));
  return(datestr);
}

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

// Outputs the variable format message (re: printf) to the log file.
int h1log(const char *fmt, ...)
{
  va_list args;
  char string[32768];

  va_start(args, fmt);
  vsprintf(string, fmt, args);
  if(h1logfp) {
    fprintf(h1logfp, string);
    fflush(h1logfp);
  } 
  else {			// couldn't open log 
    printf(string);		// also print to stdout
    fflush(stdout);
  }
  va_end(args);
  return(0);
}

int send_mail(char *fmt, ...)
{
  va_list args;
  char string[1024], cmd[1024];

  va_start(args, fmt);
  vsprintf(string, fmt, args);
  sprintf(cmd, "echo \"%s\" | Mail -s \"build_lev1 mail\" lev0_user", string);
  system(cmd);
  va_end(args);
  return(0);
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit build_lev1_mgr w/ status = %d\n", stat);
  if (h1logfp) fclose(h1logfp);
  exit(stat);
}

int forkstream(long long recn0, long long maxrecn0)
{
  pid_t pid, wpid;
  long long numofrecs, frec, lrec;
  int stat_loc, i, j, k, l;
  char *args[8], pcmd[128];
  char args1[128], args2[128], args3[128], args4[128], args5[128], args6[128];

  numofrecs = (maxrecn0 - recn0) + 1;
  //numofrecs = 192;		//!!TEMP force this for now

  i = numofrecs/MAXRECLEV1;	//this many to fork
  if(i == 0) return(0);
  lrec = recn0-1;
  for(j=i; j > 0; j = j - MAXCPULEV1) {
    if(j < MAXCPULEV1) l = j;
    else l = MAXCPULEV1;
    for(k=0; k < l; k++) { //fork this many at a time
      frec = lrec+1; lrec = (frec + MAXRECLEV1)-1;
      //printf("!!TBD for j=%d brec=%d erec=%d\n", j, frec, lrec);
      if((pid = fork()) < 0) {
        printk("***Can't fork(). errno=%d\n", errno);
        return(1);
      }
      else if(pid == 0) {                   /* this is the beloved child */
        args[0] = "build_lev1";
        sprintf(args1, "dsin=%s", dsin);
        args[1] = args1;
        sprintf(args2, "dsout=%s", dsout);
        args[2] = args2;
        sprintf(args3, "brec=%lu", frec);
        args[3] = args3;
        sprintf(args4, "erec=%lu", lrec);
        args[4] = args4;
        sprintf(args5, "instru=%s", instru);
        args[5] = args5;
        //sprintf(args6, "logfile=%s/l1_%s_%d.log", QSUBDIR, gettimetag(), k);
        sprintf(args6, "logfile=%s/l1s_%lu_%lu.log", QSUBDIR, frec, lrec);
        args[6] = args6;
        args[7] = NULL;
        printk("execvp: %s %s %s %s %s %s %s\n", 
		args[0],args[1],args[2],args[3],args[4],args[5],args[6]);
        if(execvp(args[0], args) < 0) {
          printk("***Can't execvp() build_lev1. errno=%d\n", errno);
          exit(1);
        }
      }
      printf("forked pid = %d\n", pid);
    }
    wpid = -1;			//wait for any child
    //wpid = pid;		//wait for this child
    //wpid = 0;			//wait for child in the same process group
    while(1) {
      pid = waitpid(wpid, &stat_loc, 0);
      if(pid == -1) {
        if(errno == ECHILD) printf("!!No More Children\n");errno=0;
        break;
      }
      else {
        printf("returned pid = %d stat_loc = %d\n", pid, stat_loc);
      }
    }
    //now update lev1_highest_lev0_recnum table with lrec
    sprintf(pcmd, "echo \"update lev1_highest_lev0_recnum set lev0recnum=%lu, date='%s' where lev0series='%s'\" | psql -h hmidb jsoc", lrec, get_datetime(), dsin);
    system(pcmd);
  }
  return(0); //!!TEMP
}

int do_ingest()
{
  FILE *fin;
  int rstatus, dstatus, ncnt, qcnt, i;
  long long recnum0, maxrecnum0, recnum1;
  uint64_t jid;
  char recrange[128], qlogname[128], qsubcmd[128], string[128];
  char astr[32], bstr[32], pcmd[128];

  ncnt = 0; qcnt = 0;
  if(stream_mode) {		//start past last lev0 rec# processed 
    sprintf(pcmd, "echo \"select lev0recnum from lev1_highest_lev0_recnum where lev0series='%s'\" | psql -h hmidb jsoc", dsin);
    fin = popen(pcmd, "r");
    while(fgets(string, sizeof string, fin)) {  //get psql return line
      if(strstr(string, "lev0recnum")) continue;
      if(strstr(string, "-----")) continue;
      sscanf(string, "%lu", &recnum0);		//get lev0 rec# 
      recnum0++;				//start at next rec#
      break;
    }
    fclose(fin);
    sprintf(pcmd, "echo \"select max(recnum) from %s\" | psql -h hmidb jsoc", dsin);
    fin = popen(pcmd, "r");
    while(fgets(string, sizeof string, fin)) {  //get psql return line
      if(strstr(string, "max")) continue;
      if(strstr(string, "-----")) continue;
      sscanf(string, "%lu", &maxrecnum0);	//get max lev0 rec# 
      break;
    }
    fclose(fin);
    printf("Stream Mode starting at lev0 recnum = %lu maxrecnum = %lu\n", 
		recnum0, maxrecnum0);
    if(recnum0 > maxrecnum0) return(0);		//nothing to do. go wait
    forkstream(recnum0, maxrecnum0);
    return(0);		//!!TEMP
    //!!TBD fix. put in a phoney brec for now
    //recnum0 = 1000;
  }
  else {			//reprocessing mode. use brec/erec
    recnum0 = brec;		//!!TEMP for now
  }
//!!TBD below... Ignore stream_mode for now!!
  brec = recnum0;
  erec = recnum0+(MAXRECLEV1-1);
  sprintf(recrange, ":#%lld-#%lld", brec, erec);
  sprintf(open_dsname, "%s[%s]", dsin, recrange);
  printk("open_dsname = %s\n", open_dsname);
  printk("#levnum recnum fsn\n");
  while(1) {
    sprintf(qlogname, "%s/qsub_%s_%d.csh", QSUBDIR, username, qcnt++);
    if((qsubfp=fopen(qlogname, "w")) == NULL) {
      fprintf(stderr, "**Can't open the qsub log file %s\n", qlogname);
      return(1);		//!!TBD
    }
    fprintf(qsubfp, "#!/bin/csh\n");
    fprintf(qsubfp, "echo \"TMPDIR = $TMPDIR\"\n");
    fprintf(qsubfp, "build_lev1 dsin=%s dsout=%s brec=%d erec=%d instru=%s logfile=%s/l1_%s_$JOB_ID.log\n", 
		dsin, dsout, brec, erec, instru, QSUBDIR, gettimetag()); 
    fclose(qsubfp);
    sprintf(qsubcmd, "qsub -o %s -e %s -q j.q %s", 
		QSUBDIR, QSUBDIR, qlogname);
    //sprintf(qsubcmd, "qsub -q j.q %s", qlogname);
    printf("%s\n", qsubcmd);
    printk("%s\n", qsubcmd);
    sprintf(qsubcmd, "%s | grep \"Your job\"", qsubcmd);
    fin = popen(qsubcmd, "r");
    while(fgets(string, sizeof string, fin)) {  //get qsub return line
      sscanf(string, "%s %s %d", astr, bstr, &jid); /* get job_id */
    }
    fclose(fin);
    printf("\$JOB_ID = %u\n", jid);
    ncnt = ncnt + MAXRECLEV1;
    if(ncnt >= DONUMRECS) { break; }	//!!TEMP only do this many records
    recnum0++;
    brec = brec + MAXRECLEV1;
    erec = brec + (MAXRECLEV1-1);
    sprintf(recrange, ":#%lld-#%lld", brec, erec);
    sprintf(open_dsname, "%s[%s]", dsin, recrange);
    printk("open_dsname = %s\n", open_dsname); //!!TEMP

  }
  return(0);
}

// Initial setup stuff called when main is first entered.
void setup()
{
  FILE *fin;
  char string[128], cwdbuf[128], idstr[256], lfile[128];
  int tpid;

  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  do_datestr();
  printk_set(h1log, h1log);	// set for printk calls 
  printk("%s\n", datestr);
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
  sprintf(string, "build_lev1_mgr started as pid=%d user=%s\n", getpid(), username);
  strcat(idstr, string);
  printk("%s", idstr);
  printf("%s", idstr);
  if(restartflg) printk("-r ");
  sprintf(arginstru, "instru=%s", instru);
  sprintf(argdsin, "dsin=%s", dsin);
  sprintf(argdsout, "dsout=%s", dsout);
  sprintf(argbrec, "brec=%lu", brec);
  sprintf(argerec, "erec=%lu", erec);
  sprintf(arglogfile, "logfile=%s", logname);
  printk("%s %s %s %s %s\n", arginstru, argdsin, argdsout, argbrec, argerec);
  printf("%s %s %s %s %s\n", arginstru, argdsin, argdsout, argbrec, argerec);
  if(!restartflg) {
    //printk("tlmseriesname=%s\nlev0seriesname=%s\n", 
    //		tlmseriesname, lev0seriesname);
  }
  sprintf(idstr,  "ps -ef | grep %s", LEV1VIEWERNAME);
  fin = popen(idstr, "r");
  while(fgets(string, sizeof string, fin)) {  //get ps line
    if(!(strstr(string, "perl"))) continue;
    sscanf(string, "%s %d", idstr, &tpid); /* get user name & process id */
    sprintf(lfile, "%s/build_lev1_mgr_restart_%d.touch", LEV1LOG_BASEDIR, tpid);
    sprintf(idstr, "/bin/touch %s", lfile);
    printk("%s\n", idstr);
    system(idstr);
  }
  umask(002);			// allow group write 
}

// Module main function. 
int main(int argc, char **argv)
{
  int wflg = 1;
  char line[80];

  if (cmdparams_parse(&cmdparams, argc, argv) == -1)
  {
     fprintf(stderr,"Error: Command line parsing failed. Aborting.\n");
     return 1;
  }
  if (nice_intro())
    return (0);
  if(!(username = (char *)getenv("USER"))) username = "nouser"; 
  instru = cmdparams_get_str(&cmdparams, "instru", NULL);
  if(strcmp(instru, "hmi") && strcmp(instru, "aia")) {
    printf("Error: instru= must be given as 'hmi' or 'aia'\n");
    return(0);
  }
  if(!strcmp(instru, "aia")) hmiaiaflg = 1;
  dsin = cmdparams_get_str(&cmdparams, "dsin", NULL);
  dsout = cmdparams_get_str(&cmdparams, "dsout", NULL);
  brec = cmdparams_get_int(&cmdparams, "brec", NULL);
  erec = cmdparams_get_int(&cmdparams, "erec", NULL);
  if(brec == 0 && erec == 0) stream_mode = 1;
  logfile = cmdparams_get_str(&cmdparams, "logfile", NULL);
  if (strcmp(dsin, NOTSPECIFIED) == 0) {
    if(hmiaiaflg == 0) dsin = LEV0SERIESNAMEHMI;
    else dsin = LEV0SERIESNAMEAIA;
  }
  if (strcmp(dsout, NOTSPECIFIED) == 0) {
    if(hmiaiaflg == 0) dsout = LEV1SERIESNAMEHMI;
    else dsout = LEV1SERIESNAMEAIA;
  }
  if(hmiaiaflg) {
    if(strstr(dsin, "hmi") || strstr(dsout, "hmi")) {
      printf("Warning: You said instru=aia but have 'hmi' in ds name?\n");
      printf("Do you want to abort this [y/n]? ");
      if(gets(line) == NULL) { return(0); }
      if(strcmp(line, "n")) { return(0); }
    }
  }
  else {
    if(strstr(dsin, "aia") || strstr(dsout, "aia")) {
      printf("Warning: You said instru=hmi but have 'aia' in ds name?\n");
      printf("Do you want to abort this [y/n]? ");
      if(gets(line) == NULL) { return(0); }
      if(strcmp(line, "n")) { return(0); }
    }
  }
  if (strcmp(logfile, NOTSPECIFIED) == 0) {
    sprintf(logname, H1LOGFILE, gettimetag());
  }
  else {
    sprintf(logname, "%s", logfile);
  }
  if(restartflg) {
    //sleep(30);		//make sure old ingest_lev0 is done
    if((h1logfp=fopen(logname, "a")) == NULL)
      fprintf(stderr, "**Can't open for append the log file %s\n", logname);
  }
  else {
    if((h1logfp=fopen(logname, "w")) == NULL)
      fprintf(stderr, "**Can't open the log file %s\n", logname);
  }
  setup();
  while(wflg) {
    if(do_ingest()) {        // loop to get files from the lev0
      printk("**ERROR: Some drms error  for %s\n", open_dsname);
    }
    sleep(5);		//wait for more lev0 to appear
    //wflg = 0;		//!!TBD
  }
  return(0);
}

