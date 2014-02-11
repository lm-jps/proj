// Test program only. Has vestigial code from ingest_lev0.c
// Will put file in ds su_production.tlm_test via drms calls
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
#include <dirent.h>
#include <unistd.h>	//for alarm(2) among other things...
#include <printk.h>
#include "packets.h"
#include "imgdecode.h"
#include "decode_hk_vcdu.h"
#include "decode_hk.h"
#include "load_hk_config_files.h"

#define LEV0SERIESNAMEHMI "su_production.lev0_test"
#define TLMSERIESNAMEHMI "su_production.tlm_test"
#define LEV0SERIESNAMEAIA "su_production.lev0_test_aia"
#define TLMSERIESNAMEAIA "su_production.tlm_test_aia"
//#define LEV0SERIESNAMEHMI "hmi.lev0_60d"
//#define TLMSERIESNAMEHMI "hmi.tlm_60d"
//#define LEV0SERIESNAMEAIA "aia.lev0_60d"
//#define TLMSERIESNAMEAIA "aia.tlm_60d"
#define H0LOGFILE "/usr/local/logs/lev0/ingest_lev0.%s.%s.%s.log"
#define PKTSZ 1788		//size of VCDU pkt
#define MAXFILES 512		//max # of file can handle in tlmdir
#define NUMTIMERS 8		//number of seperate timers avail
#define IMAGE_NUM_COMMIT 12	//number of complete images until commit
#define TESTAPPID 0x199		//appid of test pattern packet
#define TESTVALUE 0xc0b		//first value in test pattern packet
#define MAXERRMSGCNT 10		//max # of err msg before skip the tlm file
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define ENVFILE "/home/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE"

extern int decode_next_hk_vcdu(unsigned short *tbuf, CCSDS_Packet_t **hk, unsigned int *Fsn);
extern int write_hk_to_drms();
extern void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg);
extern int set_HMI_mech_values(DRMS_Record_t *rec);
extern TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

char *module_name = "test0";

FILE *h0logfp;		// fp for h0 ouput log for this run 
IMG Image, ImageOld;
IMG *Img;
static CCSDS_Packet_t *Hk;
static DRMS_Record_t *rs;
static DRMS_Segment_t *segment;
static DRMS_Array_t *segArray;
static DRMS_RecordSet_t *rset;
static DRMS_Record_t *rs_old, *rsc;
static DRMS_Segment_t *segmentc;
static DRMS_Array_t *cArray, *oldArray;
static char datestr[32];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];

unsigned int fsn = 0;
unsigned int fsnx = 0;
unsigned int fsn_prev = 0;
unsigned int fsn_pre_rexmit = 0;
unsigned int fid = 0;

short *rdat;
long long vcdu_seq_num;
long long vcdu_seq_num_next;
long long total_missing_im_pdu;
unsigned int vcdu_24_cnt, vcdu_24_cnt_next;
int verbose;
int appid;
int hmiaiaflg;			// 0=hmi, 1=aia
int whk_status;
int total_tlm_vcdu;
int total_missing_vcdu;
int errmsgcnt, fileimgcnt;
int imagecnt = 0;		// num of images since last commit 
int tmpcnt = 0;			//!!TEMP
int callcnt = 0;
int ALRMSEC = 70;               // seconds for alarm signal 
char timetag[32];
char pchan[8];			// primary channel to listen to e.g. VC02 
char rchan[8];			// redundant channel to listen to e.g. VC10 
char tlmseriesname[96];		// e.g. hmi.tlm
char lev0seriesname[96];	// e.g. hmi.lev0
char tlmnamekey[96];		// shortened tlm file name for TLMDSNAM keyword
char tlmnamekeyfirst[96];	// for TLMDSNAM keyword for 1st file of image
char *username;			// from getenv("USER") 
char *tlmdir;			// tlm dir name passed in 
char *outdir;			// output dir for .tlm file (can be /dev/null)
char *logfile;			// optional log name passed in 
char *vc;			// virtual channel to process, e.g. VC02 
struct p_r_chans {
  char *pchan;
  char *rchan;
  int instru;		//0=hmi, 1=aia
};
typedef struct p_r_chans P_R_CHANS;

P_R_CHANS p_r_chan_pairs[] = {
{"VC01", "VC09", 1},		// AIA 
{"VC04", "VC12", 1},		// AIA 
{"VC02", "VC10", 0},		// HMI 
{"VC05", "VC13", 0},		// HMI 
{"n/a", "n/a"}
};

// List of default parameter values.
ModuleArgs_t module_args[] = {
  {ARG_STRING, "vc", NOTSPECIFIED, "Primary virt channel to listen to"},
  {ARG_STRING, "indir", NOTSPECIFIED, "directory containing the files to ingest"},
  {ARG_STRING, "outdir", NOTSPECIFIED, "directory to move the files to after the ingest"},
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_END}
};

CmdParams_t cmdparams;

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

// Got a fatal error. 
void abortit(int stat)
{
  printf("***Abort in progress ...\n");
  printf("**Exit test0 w/ status = %d\n", stat);
  exit(stat);
}

// Called 60 secs after the last .tlm file was seen.
// Will close any opened image. NOTE: a reprocessed image cannot
// be active. It is always closed at the end of get_tlm().
void alrm_sig(int sig)
{
  signal(SIGALRM, alrm_sig);
  printf("alrm_sig: drms_server_end_transaction()\n");
  drms_server_end_transaction(drms_env, 0 , 0); //commit
  printf("alrm_sig: drms_server_begin_transaction()\n");
  drms_server_begin_transaction(drms_env); //start another cycle
}

// This is the main loop that gets the .qac and .tlm files and 
// puts them into ds tlmseriesname and extracts the lev0 and puts it in 
// lev0seriesname in DRMS.
void do_ingest()
{
  DRMS_Record_t *rs_tlm;
  int i, j, status;
  pid_t pid;
  char *args[5];
  char path[DRMS_MAXPATHLEN];
  char name[128], line[128], tlmfile[128], tlmname[96];
  char cmd[128], xxname[128], vername[16];
  char *token, *filename;

/*********************************************************************
  for(j=0; j < 1024; j++) {
      printf("drms_server_end_transaction()\n");
      drms_server_end_transaction(drms_env, 0 , 0);
      printf("drms_server_begin_transaction()\n");
      drms_server_begin_transaction(drms_env); //start another cycle
      printf("Cycle end for %d calls j=%d\n", ++callcnt, j);
  }
  exit(0);
*********************************************************************/

  sprintf(name, "%s", "/dds/stage/VC05_2008_148_15_51_05_00078703ab8_1c298_00.qac");
  sprintf(tlmfile, "%s", "/dds/stage/VC05_2008_148_15_51_05_00078703ab8_1c298_00.tlm");
  sprintf(tlmname, "%s", "VC05_2008_148_15_51_05_00078703ab8_1c298_00");
  sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMI);
  for(j=0; j < 10; j++) {
    rs_tlm = drms_create_record(drms_env, tlmseriesname, 
				DRMS_PERMANENT, &status);
    if(status) {
      printf("***Can't create record for %s\n\n", tlmseriesname);
      //continue;
      abortit(3);
    }
    printf("create_record ok for %d calls\n", ++callcnt);
    strcpy((char *)tlmnamekey, (char *)tlmname);
    if((status = drms_setkey_string(rs_tlm, "filename", tlmnamekey))) {
      printf("**ERROR: drms_setkey_string failed for 'filename'\n");
    }
    drms_record_directory(rs_tlm, path, 0);
    if(!*path) {
      printf("***ERROR: No path to segment for %s\n", tlmseriesname);
      abortit(3);
    }
    sprintf(cmd, "/bin/cp %s %s", name, path);
    printf("*cp qac to %s\n", path);
    //printf("%s\n", cmd);
    if(status = system(cmd)) {
      printf("**ERROR: %d on: %s\n", status, cmd);
    }
/********************************************************
    sprintf(cmd, "/bin/cp %s %s", tlmfile, path);
    printf("*cp tlm to %s\n", path);
    printf("%s\n", cmd);
    if(system(cmd)) {
      printf("**ERROR: on: %s\n", cmd);
    }
********************************************************/
    if((status = drms_close_record(rs_tlm, DRMS_INSERT_RECORD))) {
      printf("**ERROR: drms_close_record failed for %s\n", tlmseriesname);
    }

    if(tmpcnt++ >= 2) { 
      printf("drms_server_end_transaction()\n");
      drms_server_end_transaction(drms_env, 0 , 0);
      printf("drms_server_begin_transaction()\n");
      drms_server_begin_transaction(drms_env); //start another cycle
      tmpcnt = 0;
    }
  }
  alarm(3);			//restart alarm if no more files come
  printf("sleep(6)\n");
  sleep(6); 	//wait for alarm timeout
}

// Module main function. 
int DoIt(void)
{
  signal(SIGALRM, alrm_sig);
  while(1) {
    do_ingest();
/***********************************
    if(callcnt > 64) {
      system("test0");
      break;
    }
*********************************/
  }
}

