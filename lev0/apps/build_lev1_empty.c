/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev1/apps/build_lev1_empty.c
 *-----------------------------------------------------------------------------
 *
 * This is called by lev1_def_gui calling lev1_def_gui_called which figures
 * out the lev1 FSN that are missing for the Ord date just processed.
 * These missing FSNs are sent here and this module makes empty hmi.lev1
 * records for them. A missing record just contains the date, fsn, and
 * QUALITY keyword with the high bit set.
 */ 

#include <jsoc_main.h>
#include <cmdparams.h>
#include <drms.h>
#include <drms_names.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>	//for umask(2)
#include <unistd.h>	//for alarm(2) among other things...
#include <printk.h>
#include <astro.h>
#include <fresize.h>
#include <gapfill.h>
#include "fftw3.h" 
#include "imgdecode.h"
#include "lev0lev1.h"
#include "quallev1.h"

//default in and out data series
#define LEV1SERIESNAMEHMI "hmi.lev1"
#define LEV1SERIESNAMEAIA "su_production.aia_lev1e"	//temp test case
#define DSFFNAMEHMI "hmi.flatfield"
//#define DSFFNAMEAIA "su_production.aia_flatfield"	//temp test case
//#define DSFFNAMEAIA "aia_test.flatfield"	//temp test case
#define DSFFNAMEAIA "aia.flatfield"

#define LEV1LOG_BASEDIR "/usr/local/logs/lev1"
#define H1LOGFILE "/usr/local/logs/lev1/build_lev1_empty.%s.log"
#define NUMTIMERS 8		//number of seperate timers avail
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define LOGTEST 0
#define CAL_HCFTID 17		//image is cal mode 

int compare_rptr(const void *a, const void *b);
static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "instru", NOTSPECIFIED, "instrument. either hmi or aia"},
  {ARG_STRING, "dsout", NOTSPECIFIED, "dataset of lev1 output"},
  {ARG_STRING, "argfile", NOTSPECIFIED, "file name with FSN args from lev1_def_gui_called"},
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_INTS, "bfsn", "0", "first lev1 fsn# to process. 0=error"},
  {ARG_INTS, "efsn", "0", "last lev1 fsn# to process. 0=error"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_END}
};

CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "build_lev1_empty";

FILE *h1logfp;		// fp for h1 ouput log for this run 
static LEV0LEV1 lev0lev1;
static LEV0LEV1 *l0l1 = &lev0lev1;
static DRMS_Record_t *rs;
static DRMS_Record_t *rs0, *rs1, *rsff, *rsbad_pix, *rec_bad_aia;
static DRMS_Record_t *rptr;
static DRMS_Segment_t *segment;
static DRMS_Segment_t *segmentff;
static DRMS_Segment_t *darkseg;
static DRMS_Segment_t *badseg;
static DRMS_Segment_t *badoutpixseg;
static DRMS_Segment_t *spikeseg;
static DRMS_Array_t *segArray;
static DRMS_RecordSet_t *rset0, *rset1, *rsetff, *rsbad_aia;
static DRMS_Array_t *Array0;
static DRMS_Array_t *Arrayff;
static DRMS_Array_t *ArrayDark;
static DRMS_Array_t *ArrayBad;
static DRMS_Array_t *ArraySpike;
static TIME sdo_epoch;
//static PTINFO *ptinfo = NULL;
//static PTINFO ptdata;
static char bld_vers[16];
static char datestr[32], datestrZ[32];
static char open_dsname[256];
static char dsffname[128];
static char path[DRMS_MAXPATHLEN], bad_pix_path[DRMS_MAXPATHLEN];
static char bad_aia_path[DRMS_MAXPATHLEN];
static char rs1_path[DRMS_MAXPATHLEN];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];
static char *orbseries = "sdo.fds_orbit_vectors";
//static char *orbseries = "sdo_ground.fds_orbit_vectors";

static int nspikes, respike, fid, aiftsid, *oldvalues, *spikelocs, *newvalues;
static int hcftid, aiagp6;
static short aifcps;
double aiascale = 1.0;

IORBIT_Info_t *IOinfo = NULL;
IORBIT_Info_t IOdata;
LIBASTRO_Error_t IOstatus;
unsigned int fsnarray[NUMRECLEV1];
unsigned int fsnx = 0;
//short data1[MAXPIXELS];
//int data1[MAXPIXELS];
float data1[MAXPIXELS];		//floats for HMI
float ftmp;
int data1A[MAXPIXELS];		//ints for AIA
int array_cosmic[16777216];	//4096x4096
double tgttimes[NUMRECLEV1];

unsigned int bfsn, efsn;	//begin and end fsn.
unsigned int bnumx, enumx;	//has the bfsn/efsn pair
int verbose;
int hmiaiaflg = 0;		//0=hmi, 1=aia
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when build_lev1 is called for a restart
int abortflg = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
//global quality flags
int flatmiss[NUMRECLEV1];
int orbmiss[NUMRECLEV1];
int asdmiss[NUMRECLEV1];
int mpdmiss[NUMRECLEV1];
int limbmiss[NUMRECLEV1];
int noimage[NUMRECLEV1];
int missflg[NUMRECLEV1];

char logname[128];
char argdsout[128], arglogfile[128], arginstru[80], argargfile[128];
char argbx[80], argex[80];
char timetag[32];
char tlmseriesname[128];	// e.g. hmi.tlm
char lev0seriesname[128];	// e.g. hmi.lev0
char *username;			// from getenv("USER") 
char *logfile;			// optional log name passed in 
char *argfile;			// arg file name passed in 
char *instru;			// instument. hmi or aia 
char *dsout;			// lev1 output dataset
char *dsff;                     // flat field dataset


int get_nspikes() { return nspikes; }
int get_respike(void) { return respike; }
int *get_spikelocs() { return spikelocs; }
int *get_oldvalues() { return oldvalues; }
int *get_newvalues() { return newvalues; }
void set_nspikes(int new_nspikes) { nspikes = new_nspikes; }
void set_spikelocs(int *new_spikelocs) { spikelocs = new_spikelocs; }
void set_oldvalues(int *new_oldvalues) { oldvalues = new_oldvalues; }
void set_newvalues(int *new_newvalues) { newvalues = new_newvalues; }

int nice_intro ()
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\nbuild_lev1_empty [-vh] "
	"instru=<hmi|aia> dsout=<lev1>\n"
	"bfsn=<fsn#> efsn=<fsn#>\n"
        "argfile=<file>\n"
	"[logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
	"instru= instrument. must be 'hmi' or 'aia'\n"
	"dsout= data set name of lev1 output\n"
	"      default hmi=su_production.hmi_lev1e   aia=su_production.aia_lev1e\n"
	"bfsn= first fsn# to process. 0=error must be given by build_lev1_mgr\n"
	"efsn= last fsn# to process. 0=error must be given by build_lev1_mgr\n"
	"argfile= file with FSN args from lev1_def_gui_called\n"
	"logfile= optional log file name. If not given uses:\n"
        "         /usr/local/logs/lev1/build_lev1_empty.<time_stamp>.log\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  //restartflg = cmdparams_get_int (&cmdparams, "r", NULL);
  return (0);
}

TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss)
{
static int firstcall = 1;
if (firstcall)
  {
  firstcall = 0;
  }
/* XXX fix build 3/18/2008, arta */
return(sdo_epoch + (TIME)sdo_s + (TIME)(sdo_ss & 0xFFFF)/65536.0);
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

// Setup global datestrZ[] like: 2010-09-28T20:19:38Z
char *do_datestrZ() {
  time_t tval;
  struct tm *t_ptr;

  tval = time(NULL);
  t_ptr = gmtime(&tval);
  sprintf(datestrZ, "%d-%02d-%02dT%02d:%02d:%02dZ", 
	  (t_ptr->tm_year+1900), (t_ptr->tm_mon+1),
	  t_ptr->tm_mday, t_ptr->tm_hour, t_ptr->tm_min, t_ptr->tm_sec);
  return(datestrZ);
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
  sprintf(cmd, "echo \"%s\" | Mail -s \"build_lev1_empty mail\" lev0_user", string);
  system(cmd);
  va_end(args);
  return(0);
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit build_lev1_empty w/ status = %d\n", stat);
  if (h1logfp) fclose(h1logfp);
  exit(stat);
}

//#include "aia_despike.c"
//#include "do_flat.c"
//#include "get_image_location.c"
//#include "limb_fit_function.c"
//#include "cosmic_ray.c"
//#include "heightformation.c"

//Called with the range to do. The args are fsn#s.
//int do_ingest(long long bbrec, long long eerec)
int do_ingest(unsigned int bbrec, unsigned int eerec)
{
  int rstatus, dstatus, lstatus, ncnt, fcnt, i, j, k, qualint, nobs;

  ncnt = (eerec - bbrec) + 1;
    rset1 = drms_create_records(drms_env, ncnt, dsout, DRMS_PERMANENT,&dstatus);
    if(dstatus) {
      printk("**ERROR: Can't create records for %s\n", dsout);
      for(j=0; j < ncnt; j++) {	 //set qual bits
        noimage[j] = 1;
      }
      return(1);
    }
  for(i=0; i < ncnt; i++) {
    rs = rset1->records[i]; 
    dstatus = drms_setkey_int(rs, "FSN", bbrec++);
    dstatus = drms_setkey_string(rs, "LEV0SERIES", "N/A");
    dstatus = drms_setkey_int(rs, "QUALITY", Q_MISSALL);
    dstatus = drms_setkey_string(rs, "BLD_VERS", bld_vers);
    dstatus = drms_setkey_time(rs, "DATE", time(0) + UNIX_EPOCH);
  }
  drms_close_records(rset1, DRMS_INSERT_RECORD); //close lev1 records
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
  gethostname(idstr, 256);
  printf("Host: %s\n", idstr);
  printk("Host: %s\n", idstr);
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
  sprintf(string, "build_lev1_empty started as pid=%d ppid=%d user=%s\n", 
		getpid(), getppid(), username);
  strcat(idstr, string);
  printk("%s", idstr);
  printf("%s", idstr);
  //if(restartflg) printk("-r ");
  sprintf(arginstru, "instru=%s", instru);
  sprintf(argdsout, "dsout=%s", dsout);
    sprintf(argbx, "bfsn=%u", bfsn);
    sprintf(argex, "efsn=%u", efsn);
  sprintf(arglogfile, "logfile=%s", logname);
  sprintf(argargfile, "argfile=%s", argfile);
  printk("%s %s %s %s %s %s\n", 
	arginstru, argdsout, argbx, argex, argargfile, arglogfile);
  printf("%s %s %s %s %s %s\n", 
	arginstru, argdsout, argbx, argex, argargfile, arglogfile);
  if(!restartflg) {
    //printk("tlmseriesname=%s\nlev0seriesname=%s\n", 
    //		tlmseriesname, lev0seriesname);
  }
  sprintf(bld_vers, "%s", jsoc_version);
  /********
  sprintf(idstr,  "ps -ef | grep %s", LEV1VIEWERNAME);
  fin = popen(idstr, "r");
  while(fgets(string, sizeof string, fin)) {  //get ps line
    if(!(strstr(string, "perl"))) continue;
    sscanf(string, "%s %d", idstr, &tpid); // get user name & process id 
    sprintf(lfile, "%s/build_lev1_empty_restart_%d.touch", LEV1LOG_BASEDIR, tpid);
    sprintf(idstr, "/bin/touch %s", lfile);
    printk("%s\n", idstr);
    system(idstr);
  }
  ********/
  umask(002);			// allow group write 
}

// Module main function. 
int DoIt(void)
{
  FILE *fplog;
  long long numofrecs, frec, lrec;
  int numrec, numofchunks, i;
  unsigned int lowfsn, highfsn;
  char line[128], scr[80];

  if (nice_intro())
    return (0);
  if(!(username = (char *)getenv("USER"))) username = "nouser"; 
  instru = cmdparams_get_str(&cmdparams, "instru", NULL);
  if(strcmp(instru, "hmi") && strcmp(instru, "aia")) {
    printf("Error: instru=%s must be given as 'hmi' or 'aia'\n", instru);
    printk("Error: instru=%s must be given as 'hmi' or 'aia'\n", instru);
    return(0);
  }
  if(!strcmp(instru, "aia")) hmiaiaflg = 1;
  dsout = cmdparams_get_str(&cmdparams, "dsout", NULL);
  bfsn = cmdparams_get_int(&cmdparams, "bfsn", NULL);
  efsn = cmdparams_get_int(&cmdparams, "efsn", NULL);
/**************************************
    if(bfsn == 0 || efsn == 0) {
      fprintf(stderr, "bfsn and efsn must be given. 0 not allowed\n");
      return(0);
    }
    if(bfsn > efsn) {
      fprintf(stderr, "bfsn must be <= efsn\n");
      return(0);
    }
    bnumx = bfsn;
    enumx = efsn;
**************************************/
  logfile = cmdparams_get_str(&cmdparams, "logfile", NULL);
  argfile = cmdparams_get_str(&cmdparams, "argfile", NULL);
  if (strcmp(dsout, NOTSPECIFIED) == 0) {
    fprintf(stderr, "You must specify dsout\n");
    return(0);
  }
  if (strcmp(argfile, NOTSPECIFIED) == 0) {
    fprintf(stderr, "You must specify argfile\n");
    return(0);
  }
  if (strcmp(logfile, NOTSPECIFIED) == 0) {
    sprintf(logname, H1LOGFILE, gettimetag());
  }
  else {
    sprintf(logname, "%s", logfile);
  }
  if((h1logfp=fopen(logname, "w")) == NULL)
    fprintf(stderr, "**Can't open the log file %s\n", logname);
  setup();

  if((fplog=fopen(argfile, "r")) == NULL) {
    fprintf(stderr, "Can't open %s\n", argfile);
    return(0);
  }
  //argfile lines look like:
  //Missing FSN: went from 12207204 to 12207207
  //Missing FSN: went from 12207603 to 12207610
  while(fgets(line, 128, fplog)) {       /* get ps lines */
     if (strstr(line, "Missing FSN:")) {
       sscanf(line, "%s %s %s %s %u %s %u", scr, scr, scr, scr, &lowfsn, scr, &highfsn);
       printf("lowfsn=%u highfsn=%u\n", lowfsn, highfsn); //!!TEMP
       if(do_ingest(lowfsn+1, highfsn-1)) {
         printf("build_lev1_empty abort\nSee log: %s\n", logname); 
         send_mail("build_lev1_empty abort\nSee log: %s\n", logname); 
         return(0);
       }
     }
  }
  fclose(fplog);
  printf("build_lev1_empty done\n"); 
  return(0);
}
