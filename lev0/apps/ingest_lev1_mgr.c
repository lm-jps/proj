/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev1/apps/ingest_lev1_mgr.c
 *-----------------------------------------------------------------------------
 *
 * This is a module that runs with DRMS and continuously processes lev0
 * filtergrams to lev1 by doing a qsub to run an ingest_lev1 for each
 * of 12 (IMAGE_NUM_COMMIT) lev0 records to be processed.
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
#include "imgdecode.h"


//default in and out data series
#define LEV0SERIESNAMEHMI "hmi.lev0e"
#define LEV0SERIESNAMEAIA "aia.lev0e"
#define LEV1SERIESNAMEHMI "su_production.hmi_lev1e"	//temp test case
#define LEV1SERIESNAMEAIA "su_production.aia_lev1e"	//temp test case

#define LEV1LOG_BASEDIR "/usr/local/logs/lev1"
#define DONUMRECS 120		//!!TEMP #of lev0 records to do and then exit
#define H1LOGFILE "/usr/local/logs/lev1/ingest_lev1_mgr.%s.log"
#define NUMTIMERS 8		//number of seperate timers avail
#define IMAGE_NUM_COMMIT 12	//!!TEMP number of lev0 images to do at a time
#define NOTSPECIFIED "***NOTSPECIFIED***"

int compare_rptr(const void *a, const void *b);

/******************************************************
extern int decode_next_hk_vcdu(unsigned short *tbuf, CCSDS_Packet_t **hk, unsigned int *Fsn);
extern int write_hk_to_drms();
extern void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg);
extern int set_HMI_mech_values(DRMS_Record_t *rec);
******************************************************/

TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "instru", NOTSPECIFIED, "instrument. either hmi or aia"},
  {ARG_STRING, "dsin", NOTSPECIFIED, "dataset of lev0 filtergrams"},
  {ARG_STRING, "dsout", NOTSPECIFIED, "dataset of lev1 output"},
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_INTS, "firstfsn", "-1", "first lev0 fsn to process. -1=start from last lev1 fsn"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_END}
};

CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "ingest_lev1_mgr";

FILE *h1logfp;		// fp for h1 ouput log for this run 
FILE *qsubfp;		// fp for qsub script
static IMG Image0, Image1;
static IMG *Img0 = &Image0;
static IMG *Img1 = &Image1;
//static CCSDS_Packet_t *Hk;
static DRMS_Record_t *rs;
static DRMS_Record_t *rs0, *rs1;
static DRMS_Record_t *rptr;
static DRMS_Segment_t *segment;
static DRMS_Array_t *segArray;
static DRMS_RecordSet_t *rset0, *rset1;
static DRMS_Array_t *Array0;
static TIME sdo_epoch;
static char datestr[32];
static char open_dsname[256];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];

unsigned int fsnarray[IMAGE_NUM_COMMIT];
unsigned int fsnx = 0;

short *rdat;
int verbose;
int firstfsn;
int hmiaiaflg = 0;		//0=hmi, 1=aia
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when ingest_lev0 is called for a restart
int abortflg = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
char logname[128];
char argdsin[128], argdsout[128], arglogfile[128], argfsn[80], arginstru[80];
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
    printf ("Usage:\ningest_lev1_mgr [-vhr] "
	"instru=<hmi|aia> dsin=<lev0> dsout=<lev1> firstfsn=<fsn>\n"
	"                   [logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
	"  -r: restart. only used when we restart our selves periodically\n"
	"instru= instrument. must be 'hmi' or 'aia'\n"
	"dsin= data set name of lev0 input\n"
	"      default hmi=hmi.lev0e   aia=aia.lev0e\n"
	"dsout= data set name of lev1 output\n"
	"      default hmi=su_production.hmi_lev1e   aia=su_production.aia_lev1e\n"
	"firstfsn= first lev0 fsn to process. -1=start from last lev1 fsn\n"
	"logfile= optional log file name. If not given uses:\n"
        "         /usr/local/logs/lev1/ingest_lev1_mgr.<time_stamp>.log\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  restartflg = cmdparams_get_int (&cmdparams, "r", NULL);
  return (0);
}

static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss)
{
static int firstcall = 1;
if (firstcall)
  {
  firstcall = 0;
  }
/* XXX fix build 3/18/2008, arta */
return(sdo_epoch + (TIME)sdo_s + (TIME)(sdo_ss)/65536.0);
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
  sprintf(cmd, "echo \"%s\" | Mail -s \"ingest_lev0 mail\" lev0_user", string);
  system(cmd);
  va_end(args);
  return(0);
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit ingest_lev1_mgr w/ status = %d\n", stat);
  if (h1logfp) fclose(h1logfp);
  exit(stat);
}

int do_ingest()
{
  DRMS_Record_t *irpt;
  FILE *fin;
  int rstatus, dstatus, ncnt, i;
  long long recnum0, recnum1;
  uint64_t jid;
  char recrange[128], qlogname[128], qsubcmd[128], string[128];
  char astr[32], bstr[32];

  ncnt = 0;		//!!TEMP use of ncnt
  fsnx = firstfsn;
  if(fsnx == -1) {	//default is to start past last lev1 record
    sprintf(open_dsname, "%s[:#$]", dsout);
    rset1 = drms_open_records(drms_env, open_dsname, &rstatus); //open last lev1
    if(!rset1 || (rset1->n == 0) || rstatus) {
      printk("Can't do drms_open_records(%s)\n", open_dsname);
      return(1);		//!!TBD
      //abortit(1);          // !!!TBD
    }
    rs1 = rset1->records[0];
    rs0 = drms_link_follow(rs1, "LEV0REC", &rstatus);
    if(rstatus) {
      printk("Can't do drms_link_follow for %s\n", open_dsname);
      return(1);		//!!TBD
    }
    drms_close_records(rset1, DRMS_FREE_RECORD); //close the last lev1 rec
    recnum0 = rs0->recnum;
    recnum0++;			//start at next recnum
    printk("Starting off with %s recnum=%lld\n", dsin, recnum0);
    printf("Starting off with %s recnum=%lld\n", dsin, recnum0);
  }
  else {			//an FSN was given as an input arg
    printk("!!TEMP starting off with %s FSN=%u\n", dsin, fsnx);
    printf("!!TEMP starting off with %s FSN=%u\n", dsin, fsnx);
    sprintf(open_dsname, "%s[%u]", dsin, fsnx);	//open lev0 for this fsn
    rset0 = drms_open_records(drms_env, open_dsname, &rstatus); //open lev0
    if(!rset0 || (rset0->n == 0) || rstatus) {
      printk("Can't do drms_open_records(%s)\n", open_dsname);
      return(1);		//!!TBD
    }
    rs0 = rset0->records[0];
    recnum0 = rs0->recnum;	//start at this recnum
    drms_close_records(rset0, DRMS_FREE_RECORD); //close the fsn lev0 rec
  }
  sprintf(recrange, ":#%lld-#%lld", recnum0, recnum0+(IMAGE_NUM_COMMIT-1));
  sprintf(open_dsname, "%s[%s]", dsin, recrange);
  printk("open_dsname = %s\n", open_dsname);
  printk("#levnum recnum fsn\n");
  while(1) {
    sprintf(qlogname, "/usr/local/logs/lev1/qsub_%s_%d.csh", username, ncnt++);
    if((qsubfp=fopen(qlogname, "w")) == NULL) {
      fprintf(stderr, "**Can't open the qsub log file %s\n", qlogname);
      return(1);		//!!TBD
    }
    fprintf(qsubfp, "#!/bin/csh\n");
    fprintf(qsubfp, "echo \"TMPDIR = $TMPDIR\"\n");
    //fprintf(qsubfp, "ingest_lev1 firstfsn=%u instru=%s logfile=$TMPDIR/jim.log\n", firstfsn, instru);
    fprintf(qsubfp, "ingest_lev1 firstfsn=%u instru=%s logfile=/home/production/qsub/lev1/l1_%s_$JOB_ID.log\n", firstfsn, instru, gettimetag()); 
    fclose(qsubfp);
    sprintf(qsubcmd, "qsub -o /home/production/qsub -e /home/production/qsub -q j.q %s | grep \"Your job\"", qlogname);
    //sprintf(qsubcmd, "qsub -q j.q %s", qlogname);
    printf("%s\n", qsubcmd);
    printk("%s\n", qsubcmd);
    //system(qsubcmd);
    fin = popen(qsubcmd, "r");
    while(fgets(string, sizeof string, fin)) {  //get qsub return line
      sscanf(string, "%s %s %d", astr, bstr, &jid); /* get job_id */
    }
    fclose(fin);
    printf("\$JOB_ID = %u\n", jid);
    //sprintf(qsubcmd, "echo \"JOB_ID = $JOB_ID\"\n");
    //system(qsubcmd);




/***************ORIG lev1 stuff below ****************************************
    rset0 = drms_open_records(drms_env, open_dsname, &rstatus); //open lev0
    if(!rset0 || (rset0->n == 0) || rstatus) {
      printk("Can't do drms_open_records(%s)\n", open_dsname);
      return(1);		//!!TBD
      //abortit(1);          // !!!TBD
    }
      rptr = (DRMS_Record_t *)malloc(IMAGE_NUM_COMMIT * sizeof(DRMS_Record_t));
      if(rptr == NULL) {
        printk("Can't malloc() for DRMS_Record_t sort\n");
        return(1);		//!!TBD
      }
      //make opened records sequencial in mem for sort
      for(i=0; i < IMAGE_NUM_COMMIT; i++) {
        memcpy(&rptr[i], rset0->records[i], sizeof(DRMS_Record_t));
      }
      //Must sort to get in ascending record # order
      qsort(rptr, IMAGE_NUM_COMMIT, sizeof(DRMS_Record_t), &compare_rptr);

      //this loop is for the benefit of the lev1view dispaly to show
      //all 12 lev0 records opened
      for(i=0; i < IMAGE_NUM_COMMIT; i++) {
        rs0 = &rptr[i];
        recnum0 = rs0->recnum;
        //must get fsn in case got record from last lev1 record
        fsnx = drms_getkey_int(rs0, "FSN", &rstatus); 
        fsnarray[i] = fsnx;
        printk("*0 %u %u\n", recnum0, fsnx);
      }

    for(i=0; i < IMAGE_NUM_COMMIT; i++) { 
      rs0 = &rptr[i];
      recnum0 = rs0->recnum;
      fsnx = fsnarray[i]; 
      //printf("rec# for %d = %lld fsn=%u\n", i, recnum0, fsnx); //!!!TEMP
      segment = drms_segment_lookupnum(rs0, 0);
      Array0 = drms_segment_read(segment, DRMS_TYPE_SHORT, &rstatus);
      if(rstatus) {
        printk("Can't do drms_segment_read() %s\n", open_dsname);
        return(1);              // !!!!TBD ck
      }
      short *adata = (short *)Array0->data;
      memcpy(Img0->dat, adata, 2*MAXPIXELS);
      drms_free_array(Array0); //must free from drms_segment_read()

      sprintf(open_dsname, "%s[%u]", dsout, fsnx);
      rs = drms_create_record(drms_env, dsout, DRMS_PERMANENT, &dstatus);
      if(dstatus) {
        printk("**ERROR: Can't create record for %s\n", open_dsname);
        return(1);
      }
      dstatus = drms_setkey_int(rs, "FSN", fsnx);
      if(dstatus = drms_setlink_static(rs, "LEV0REC", recnum0)) {
        printk("**ERROR on drms_setlink_static() for %s\n", open_dsname);
        return(1);
      }
      if(!(segment = drms_segment_lookup(rs, "image"))) {
        printk("No drms_segment_lookup(rs, image) for %s\n", open_dsname);
        return(1);
      }
      rdat = Img1->dat;
      segArray = drms_array_create(DRMS_TYPE_SHORT,
                                       segment->info->naxis,
                                       segment->axis,
                                       rdat,
                                       &dstatus);
      //call_keh_cheng();	//!!TBD

      dstatus = drms_segment_write(segment, segArray, 0);
      if (dstatus) {
        printk("ERROR: drms_segment_write error=%d for fsn=%u\n", dstatus,fsnx);
      }
      recnum1 = rs->recnum;
      if((dstatus = drms_close_record(rs, DRMS_INSERT_RECORD))) {
        printk("**ERROR: drms_close_record failed for %s\n", open_dsname);
        return(1);
      }
      printk("*1 %u %u\n", recnum1, fsnx);
    }

    drms_close_records(rset0, DRMS_FREE_RECORD);
********************************END orig stuff *************************/
    ncnt = ncnt + IMAGE_NUM_COMMIT;
    if(ncnt >= DONUMRECS) { break; }	//!!TEMP only do this many records
    recnum0++;
    sprintf(recrange, ":#%lld-#%lld", recnum0, recnum0+(IMAGE_NUM_COMMIT-1));
    sprintf(open_dsname, "%s[%s]", dsin, recrange);
    printk("open_dsname = %s\n", open_dsname); //!!TEMP

  }
  return(0);
}

int compare_rptr(const void *a, const void *b)
{
  DRMS_Record_t *x=(DRMS_Record_t *)a, *y=(DRMS_Record_t *)b;

  if(x->recnum < y->recnum) return(-1);
  if(x->recnum > y->recnum) return(1);
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
  sprintf(string, "ingest_lev1_mgr started as pid=%d user=%s\n", getpid(), username);
  strcat(idstr, string);
  printk("%s", idstr);
  printf("%s", idstr);
  if(restartflg) printk("-r ");
  sprintf(arginstru, "instru=%s", instru);
  sprintf(argdsin, "dsin=%s", dsin);
  sprintf(argdsout, "dsout=%s", dsout);
  sprintf(argfsn, "firstfsn=%d", firstfsn);
  sprintf(arglogfile, "logfile=%s", logname);
  printk("%s %s %s\n", arginstru, argdsin, argdsout);
  printf("%s %s %s\n", arginstru, argdsin, argdsout);
  if(!restartflg) {
    //printk("tlmseriesname=%s\nlev0seriesname=%s\n", 
    //		tlmseriesname, lev0seriesname);
  }
  sprintf(idstr,  "ps -ef | grep %s", LEV1VIEWERNAME);
  fin = popen(idstr, "r");
  while(fgets(string, sizeof string, fin)) {  //get ps line
    if(!(strstr(string, "perl"))) continue;
    sscanf(string, "%s %d", idstr, &tpid); /* get user name & process id */
    sprintf(lfile, "%s/ingest_lev1_mgr_restart_%d.touch", LEV1LOG_BASEDIR, tpid);
    sprintf(idstr, "/bin/touch %s", lfile);
    printk("%s\n", idstr);
    system(idstr);
  }
  umask(002);			// allow group write 
  //Image.initialized = 0;	// init the two image structures 
  //ImageOld.initialized = 0;
  //Img = &Image;
  //ImgO = &ImageOld;
  //Img->initialized = 0;
  //ImgO->initialized = 0;
}

// Module main function. 
int DoIt(void)
{
  int wflg = 1;
  char line[80];

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
  firstfsn = cmdparams_get_int(&cmdparams, "firstfsn", NULL);
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
    //sleep(1);
    wflg = 0;		//!!TBD
  }
  return(0);
}

