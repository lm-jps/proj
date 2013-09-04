/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev0/apps/ingest_lev0.c
 * NOTE: This originally came from hmi_lev0.c on hmi0
 *-----------------------------------------------------------------------------
 *
 * This is a module that runs with DRMS and continuously extracts images and HK
 * data from .tlm files that appear in the given input dir. It puts the .tlm
 * and .qac files in the DRMS dataset TLMSERIESNAME and also in the outdir if
 * given.  It extracts images from the .tlm files an puts them in
 * the DRMS dataset LEV0SERIESNAME and extracts hk data to appropriate hk datasets.
 *
 * Call for testing (set TLMSERIESNAME and LEV0SERIESNAME to test datasets): 
 * ingest_lev0 vc=VC05 indir=/scr20/jim/tlm outdir=/scr20/jim/tlm/out [logfile=name]
 *
 * OPERATION:
 * This is normally run on the dedicated lev0 host. Four instances are run.
 * No outdir= is given, i.e. the .tlm and .qac files only end up in SUMS.
 *
 * ingest_lev0 vc=VC02 indir=/dds/socdc/hmi
 * ingest_lev0 vc=VC05 indir=/dds/socdc/hmi
 * ingest_lev0 vc=VC10 indir=/dds/socdc/aia
 * ingest_lev0 vc=VC13 indir=/dds/socdc/aia
 */ 

#include <jsoc_main.h>
#include <cmdparams.h>
#include <drms.h>
#include <drms_names.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <strings.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>	//for umask(2)
#include <dirent.h>
#include <unistd.h>	//for alarm(2) among other things...
#include <printk.h>
#include <math.h>
#include "packets.h"
#include "imgdecode.h"
#include "decode_hk_vcdu.h"
#include "decode_hk.h"
#include "load_hk_config_files.h"
#include "add_small_image.c"
#include "mypng.h"
#include "tdsignals.h"
#include "quallev0.h"

#define RESTART_CNT 2	//#of tlm files to process before restart

//#define LEV0SERIESNAMEHMI "su_production.lev0d_test"
//#define LEV0SERIESNAMEHMI "su_production.lev0f_test"
//#define LEV0SERIESNAMEHMI "su_production.lev0f_hmi"
////#define LEV0SERIESNAMEHMI "su_production.lev0f_hmi_test"
//For test with DDS on Jan 19, 2010
//#define LEV0SERIESNAMEHMI "su_production.lev0f_hmi_JAN2010"
//#define LEV0SERIESNAMEHMI "su_production.lev0f_hmi_junk"
//#define LEV0SERIESNAMEHMI "hmi.lev0f"
////#define TLMSERIESNAMEHMI "su_production.tlm_test"
//For test with DDS on Jan 19, 2010
//#define TLMSERIESNAMEHMI "su_production.tlm_hmi_JAN2010"
//#define TLMSERIESNAMEHMI "su_production.tlm_hmi_junk"

#define LEV0SERIESNAMEHMIGND "hmi_ground.lev0_dds"
#define TLMSERIESNAMEHMIGND "hmi_ground.tlm_dds"
//#define TLMSERIESNAMEHMI "hmi.tlm_reingest"

//#define LEV0SERIESNAMEAIA "su_production.lev0d_test_aia"
//#define LEV0SERIESNAMEAIA "su_production.lev0f_aia"
////#define LEV0SERIESNAMEAIA "aia.lev0f"
//For test with DDS on Jan 19, 2010
//#define LEV0SERIESNAMEAIA "su_production.lev0f_aia_JAN2010"
//#define LEV0SERIESNAMEAIA "su_production.lev0f_aia_junk"
////#define TLMSERIESNAMEAIA "su_production.tlm_test_aia"
//For test with DDS on Jan 19, 2010
//#define TLMSERIESNAMEAIA "su_production.tlm_aia_JAN2010"
//#define TLMSERIESNAMEAIA "su_production.tlm_aia_junk"

#define LEV0SERIESNAMEAIAGND "aia_ground.lev0_dds"
#define TLMSERIESNAMEAIAGND "aia_ground.tlm_dds"
//#define TLMSERIESNAMEAIA "aia.tlm_reingest"

//#define LEV0SERIESNAMEHMI "hmi.lev0_60d"
//#define TLMSERIESNAMEHMI "hmi.tlm_60d"
//#define LEV0SERIESNAMEAIA "aia.lev0_60d"
//#define TLMSERIESNAMEAIA "aia.tlm_60d"

//When change to these data series below to save real data.
#define TLMSERIESNAMEHMI "hmi.tlm"
#define LEV0SERIESNAMEHMI "hmi.lev0a"
#define LEV0SERIESNAMEAIA "aia.lev0"
#define TLMSERIESNAMEAIA "aia.tlm"

//#define LEV0SERIESNAMEAIA "aia.lev0d"
//#define TLMSERIESNAMEAIA "aia.tlmd"
//#define LEV0SERIESNAMEHMI "hmi.lev0d"
//#define TLMSERIESNAMEHMI "hmi.tlmd"

//#define LEV0SERIESNAMEAIA "aia.lev0e"
//#define TLMSERIESNAMEAIA "aia.tlme"
//#define LEV0SERIESNAMEHMI "hmi.lev0e"
//#define TLMSERIESNAMEHMI "hmi.tlme"

//Also, change setting in $JSOCROOT/proj/lev0/apps/SOURCE_ENV_HK_DECODE file to:
//setenv HK_LEV0_BY_APID_DATA_ID_NAME      lev0
//setenv HK_DF_HSB_DIRECTORY               /tmp21/production/lev0/hk_hsb_dayfile
//setenv HK_JSVNMAP_DIRECTORY              /home/production/cvs/TBL_JSOC/lev0/hk_jsn_map_file/prod


#define H0LOGFILE "/usr/local/logs/lev0/ingest_lev0.%s.%s.%s.log"
#define PKTSZ 1788		//size of VCDU pkt
#define MAXFILES 65535		//max # of file can handle in tlmdir
#define NUMTIMERS 8		//number of seperate timers avail
#define IMAGE_NUM_COMMIT 12	//number of complete images until commit
//#define IMAGE_NUM_COMMIT 2	//!!TEMP number of complete images until commit
#define TESTAPPID 0x199		//appid of test pattern packet
#define TESTVALUE 0xc0b		//first value in test pattern packet
#define MAXERRMSGCNT 10		//max # of err msg before skip the tlm file
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define ENVFILE "/home2/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE"
#define ENVFILE_GND "/home2/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE_GROUND"

extern int decode_next_hk_vcdu(unsigned short *tbuf, CCSDS_Packet_t **hk, unsigned int *Fsn);
extern int write_hk_to_drms();
extern void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg);
extern int set_HMI_mech_values(DRMS_Record_t *rec);
static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "vc", NOTSPECIFIED, "Primary virt channel to listen to"},
  {ARG_STRING, "indir", NOTSPECIFIED, "directory containing the files to ingest"},
  {ARG_STRING, "outdir", NOTSPECIFIED, "directory to move the files to after the ingest"}, 
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_FLAG, "g", "0", "ground data flag"},
  {ARG_END}
};

CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "ingest_lev0";

typedef enum
{
   kThreadSigErr_Success
} ThreadSigErr_t;
static int gLoop = 1;
static pthread_mutex_t mutex;
static struct timeval tv0;
static td_alarm_t talarm = 0;

FILE *h0logfp;		// fp for h0 ouput log for this run 
static IMG Image, ImageOld;
static IMG *Img, *ImgO, *ImgC;
static CCSDS_Packet_t *Hk;
static DRMS_Record_t *rs;
static DRMS_Segment_t *segment;
static DRMS_Array_t *segArray;
static DRMS_RecordSet_t *rset;
static DRMS_Record_t *rs_old, *rsc;
static DRMS_Segment_t *segmentc;
static DRMS_Array_t *cArray, *oldArray;
static TIME sdo_epoch;
static char datestr[32];
static char bld_vers[16];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];

int INVALtime;
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
int grounddata;
int appid;
int testid1, testid2;
int hmiaiaflg;			//0=hmi, 1=aia
int whk_status;
int total_tlm_vcdu;
int total_missing_vcdu;
int errmsgcnt, fileimgcnt;
int cntsleeps = 0;
int paused = 0;
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when ingest_lev0 is called for a restart
int abortflg = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
int firstfound = 0;		// set if see any file after startup
int rexmitmode = 0;		// set if this process is doing /rexmit dir
int ALRMSEC = 60;               // must get 2 in a row for no image timeout
int sleep_interval = 2;		// #of sec to sleep after do_ingest() calls
char logname[128];
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
char *username;			// from getenv("USER") 
char *tlmdir;			// tlm dir name passed in 
char *outdir;			// output dir for .tlm file (can be /dev/null)
char *logfile;			// optional log name passed in 
char *vc;			// virtual channel to process, e.g. VC02 
struct stat stbuf;
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

struct namesort {		// sorted file names in tlmdir 
  char *name;
};
typedef struct namesort NAMESORT;


//Function callback that gets called when alarm 'rings'.  This runs in the signal
//thread that receives, the alarm signal, not the main thread. */
void shandler(int sig, pthread_mutex_t *mtx)
{
   float elapsed;
   struct timeval tv1;
   gettimeofday(&tv1, NULL);

   elapsed = (float)((tv1.tv_sec * 1000000.0 + tv1.tv_usec -
                      (tv0.tv_sec * 1000000.0 + tv0.tv_usec)) / 1000000.0);

   printk("Thread '%lld' received alarm signal '%d'.\n", 
		(long long )pthread_self(), sig);
   printk("Elapsed time is %f seconds.\n", elapsed);

   /* This isn't in the main thread, so put globals in a critical region. */
   pthread_mutex_lock(mtx);
   sigalrmflg = 1;
   pthread_mutex_unlock(mtx);
}

int nice_intro ()
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\ningest_lev0 [-vh] "
	"vc=<virt chan> indir=</dir> [outdir=</dir>] [logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
        "  -g: output to hmi_ground/aia_ground\n"
	"  -r: restart. only used when we restart our selves periodically\n"
	"vc= primary virt channel to listen to e.g. VC02\n"
	"indir= directory containing the files to ingest\n"
	"outdir= optional dir to copy the files to after the ingest\n"
	"logfile= optional log file name. Will create one if not given\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  grounddata = cmdparams_get_int (&cmdparams, "g", NULL);
  restartflg = cmdparams_get_int (&cmdparams, "r", NULL);
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
int h0log(const char *fmt, ...)
{
  va_list args;
  char string[32768];

  va_start(args, fmt);
  vsprintf(string, fmt, args);
  if(h0logfp) {
    fprintf(h0logfp, string);
    fflush(h0logfp);
  } else {			// couldn't open log 
    printf(string);
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
  sprintf(cmd, "echo \"%s\" | Mail -s \"ingest_lev0 mail\" jsoc_ops", string);
  system(cmd);
  va_end(args);
  return(0);
}

//Set the QUALLEV0 keyword
void do_quallev0(DRMS_Record_t *rs, IMG *img, int fsn) 
{
  char *hwltnset, *cdark;
  char *hseqerr, *aistate;
  char wave_str[16];
  int status, hsqfgsn, asqfsn, hgp1rgst;
  float percentd;
  uint32_t missvals, datav;
  uint32_t hcf1encd, hcf2encd, hps1encd, hps2encd, hps3encd; 
  uint32_t hwt1encd, hwt2encd, hwt3encd, hwt4encd;
  uint32_t hcf1pos, hcf2pos, hpl1pos, hpl2pos, hpl3pos, hwl1pos, hwl2pos;
  uint32_t hwl3pos, hwl4pos;
  uint32_t aiawvlen, aifwen, aiasen;
  short aifiltyp;
  uint32_t quallev0 = 0;
  uint32_t qmiss = 0;

  if(img->overflow) quallev0 = quallev0 | Q_OVFL;
  if(img->headerr) quallev0 = quallev0 | Q_HDRERR;
  if(img->nerrors) quallev0 = quallev0 | Q_CMPERR;
  if(img->last_pix_err) quallev0 = quallev0 | Q_LPXERR;
  if(img->reopened) quallev0 = quallev0 | Q_REOPENED;
  missvals = img->totalvals - img->datavals;
  if(missvals > 0) quallev0 = quallev0 | Q_MISS0;
  datav = img->totalvals;
  if(missvals > (uint32_t)(datav * 0.01)) quallev0 = quallev0 | Q_MISS1;
  if(missvals > (uint32_t)(datav * 0.05)) quallev0 = quallev0 | Q_MISS2;
  if(missvals > (uint32_t)(datav * 0.25)) quallev0 = quallev0 | Q_MISS3;
  if(missvals == datav) quallev0 = quallev0 | Q_MISSI;
  if(img->datavals == 0) quallev0 = quallev0 | Q_MISSALL; //high bit, no data

//  if(missvals != 0) {
//    qmiss = (uint32_t)(0.84 * log(missvals/img->totalvals)) + 15;
//    qmiss = qmiss & 0x0f;		//need this?
//    qmiss = qmiss * 256;		//bits 8-11
//    quallev0 = quallev0 | qmiss;
//  }

if(fsn == 469769216) quallev0 = quallev0 | Q_CORRUPT; //corrupt image
if(INVALtime) quallev0 = quallev0 | Q_INVALTIME; 
if(cdark = drms_getkey_string(rs, "IMG_TYPE", &status)) {
  if(!strcmp(cdark, "DARK")) quallev0 = quallev0 | Q_DARK;
  free(cdark);
}
if(!hmiaiaflg) {		//HMI specific qual bits
  hsqfgsn = drms_getkey_int(rs, "HSQFGSN", &status);
  if(status || (fsn != hsqfgsn)) quallev0 = quallev0 | Q_NOISP;
  //Removed per Rock's email Re: lev0 quality updates 9Jun2010 10:22
  //if(hseqerr = drms_getkey_string(rs, "HSEQERR", &status)) {
  //  if(strcmp(hseqerr, "SUCCESS")) quallev0 = quallev0 | Q_SEQERR;
  //  free(hseqerr);
  //}
  if(hwltnset = drms_getkey_string(rs, "HWLTNSET", &status)) {
    if(!strcmp(hwltnset, "OPEN")) quallev0 = quallev0 | Q_ISSOPEN;
    free(hwltnset);
  }
  hcf1encd = drms_getkey_int(rs, "HCF1ENCD", &status);
  hcf1pos = drms_getkey_int(rs, "HCF1POS", &status);
  //if(!((hcf1encd == hcf1pos) || (hcf1encd == hcf1pos+1) || (hcf1encd == hcf1pos-1)))
  if(!((hcf1encd == hcf1pos) || (hcf1encd == hcf1pos-1)))
    quallev0 = quallev0 | Q_HCF1ENCD;
  hcf2encd = drms_getkey_int(rs, "HCF2ENCD", &status);
  hcf2pos = drms_getkey_int(rs, "HCF2POS", &status);
  //if(!((hcf2encd == hcf2pos) || (hcf2encd == hcf2pos+1) || (hcf2encd == hcf2pos-1)))
  if(!((hcf2encd == hcf2pos) || (hcf2encd == hcf2pos-1)))
    quallev0 = quallev0 | Q_HCF2ENCD;
  hps1encd = drms_getkey_int(rs, "HPS1ENCD", &status);
  hpl1pos = drms_getkey_int(rs, "HPL1POS", &status);
  if(!((hpl1pos == hps1encd) || (hpl1pos == (hps1encd+1) % 240)))
    quallev0 = quallev0 | Q_HPS1ENCD;
  hps2encd = drms_getkey_int(rs, "HPS2ENCD", &status);
  hpl2pos = drms_getkey_int(rs, "HPL2POS", &status);
  if(!((hpl2pos == hps2encd) || (hpl2pos == (hps2encd+1) % 240)))
    quallev0 = quallev0 | Q_HPS2ENCD;
  hps3encd = drms_getkey_int(rs, "HPS3ENCD", &status);
  hpl3pos = drms_getkey_int(rs, "HPL3POS", &status);
  if(!((hpl3pos == hps3encd) || (hpl3pos == (hps3encd+1) % 240)))
    quallev0 = quallev0 | Q_HPS3ENCD;
  hwt1encd = drms_getkey_int(rs, "HWT1ENCD", &status);
  hwl1pos = drms_getkey_int(rs, "HWL1POS", &status);
  if(!((hwl1pos == hwt1encd) || (hwl1pos == (hwt1encd+1) % 240)))
    quallev0 = quallev0 | Q_HWT1ENCD;
  hwt2encd = drms_getkey_int(rs, "HWT2ENCD", &status);
  hwl2pos = drms_getkey_int(rs, "HWL2POS", &status);
  if(!((hwl2pos == hwt2encd) || (hwl2pos == (hwt2encd+1) % 240)))
    quallev0 = quallev0 | Q_HWT2ENCD;
  hwt3encd = drms_getkey_int(rs, "HWT3ENCD", &status);
  hwl3pos = drms_getkey_int(rs, "HWL3POS", &status);
  if(!((hwl3pos == hwt3encd) || (hwl3pos == (hwt3encd+1) % 240)))
    quallev0 = quallev0 | Q_HWT3ENCD;
  hwt4encd = drms_getkey_int(rs, "HWT4ENCD", &status);
  hwl4pos = drms_getkey_int(rs, "HWL4POS", &status);
  if(!((hwl4pos == hwt4encd) || (hwl4pos == (hwt4encd+1) % 240)))
    quallev0 = quallev0 | Q_HWT4ENCD;
  //New 14Aug2013. Set Q_GPREGBIT0/1 
  hgp1rgst = drms_getkey_int(rs, "HGP1RGST", &status);
  if(hgp1rgst != DRMS_MISSING_INT) {
    hgp1rgst = (hgp1rgst << 28) & 0x30000000;
    quallev0 = quallev0 | hgp1rgst;
  }
}
else {				//AIA specific qual bits
  asqfsn = drms_getkey_int(rs, "ASQFSN", &status);
  if(fsn != asqfsn) quallev0 = quallev0 | Q_NOISP;
  if(aistate = drms_getkey_string(rs, "AISTATE", &status)) {
    if(!strcmp(aistate, "OPEN")) quallev0 = quallev0 | AQ_ISSOPEN;
    free(aistate);
  }
  strcpy(wave_str, "UNKNOWN");
  aiawvlen = drms_getkey_int(rs, "AIAWVLEN", &status);
  aifiltyp = drms_getkey_short(rs, "AIFILTYP", &status);
  aifwen = drms_getkey_int(rs, "AIFWEN", &status);
  aiasen = drms_getkey_int(rs, "AIASEN", &status);
  switch(aiawvlen) {
  case 9:			//9.4
    if(aifiltyp == 0) {
      if((aifwen != 269) && (aifwen != 270)) {
        quallev0 = quallev0 | A94Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if((aifwen != 11) && (aifwen != 12)) {
        quallev0 = quallev0 | A94Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A94Mech_Err;
      }
    }
    //All the wave_str code can be compacted, but I'm keeping it 
    //symetric with Rock's notes
    if((aifwen == 269) || (aifwen == 270)) {
      strcpy(wave_str, "94_THIN");
    }
    else if((aifwen == 11) || (aifwen == 12)) {
      strcpy(wave_str, "94_THICK");
    }
    else if((aifwen == 74) || (aifwen == 75)) {
      strcpy(wave_str, "94_OPEN");
    }
    break;
  case 1:			//13.1
    if(aifiltyp == 0) {
      if((aifwen != 269) && (aifwen != 270)) {
        quallev0 = quallev0 | A131Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if((aifwen != 11) && (aifwen != 12)) {
        quallev0 = quallev0 | A131Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A131Mech_Err;
      }
    }
    if((aifwen == 269) || (aifwen == 270)) {
      strcpy(wave_str, "131_THIN");
    }
    else if((aifwen == 11) || (aifwen == 12)) {
      strcpy(wave_str, "131_THICK");
    }
    else if((aifwen == 74) || (aifwen == 75)) {
      strcpy(wave_str, "131_OPEN");
    }
    break;
  case 7:			//17.1
    if(aifiltyp == 0) {
      if((aifwen != 203) && (aifwen != 204)) {
        quallev0 = quallev0 | A171Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if((aifwen != 11) && (aifwen != 12)) {
        quallev0 = quallev0 | A171Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      //NOTE: no 171_OPEN
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A171Mech_Err;
      }
    }
    if((aifwen == 203) || (aifwen == 204)) {
      strcpy(wave_str, "171_THIN");
    }
    else if((aifwen == 11) || (aifwen == 12)) {
      strcpy(wave_str, "171_THICK");
    }
    break;
/* !!!TBD fix below like did above after resolve MIX_THIN stuff w/Rock */
  case 3:			//19.3
    if(aifiltyp == 0) {
      if(aiasen != 6) {
        quallev0 = quallev0 | A193Mech_Err;
      }
      if((aifwen != 269) && (aifwen != 270)) {
        quallev0 = quallev0 | A193Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if(aiasen != 6) {
        quallev0 = quallev0 | A193Mech_Err;
      }
      if((aifwen != 11) && (aifwen != 12)) {
        quallev0 = quallev0 | A193Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      if(aiasen != 6) {
        quallev0 = quallev0 | A193Mech_Err;
      }
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A193Mech_Err;
      }
    }
    if(aiasen == 6) {
      if((aifwen == 269) || (aifwen == 270)) {
        strcpy(wave_str, "193_THIN");
      } else if((aifwen == 11) || (aifwen == 12)) {
        strcpy(wave_str, "193_THICK");
      } else if((aifwen == 74) || (aifwen == 75)) {
        strcpy(wave_str, "193_OPEN");
      }
    }
    else {
      if((aifwen == 269) || (aifwen == 270)) {
        strcpy(wave_str, "MIX_THIN");
      } else if((aifwen == 11) || (aifwen == 12)) {
        strcpy(wave_str, "MIX_THICK");
      } else if((aifwen == 74) || (aifwen == 75)) {
        strcpy(wave_str, "MIX_OPEN");
      }
    }
    break;
  case 2:			//21.1
    if(aifiltyp == 0) {
      if(aiasen != 24) {
        quallev0 = quallev0 | A211Mech_Err;
      }
      if((aifwen != 203) && (aifwen != 204)) {
        quallev0 = quallev0 | A211Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if(aiasen != 24) {
        quallev0 = quallev0 | A211Mech_Err;
      }
      if((aifwen != 137) && (aifwen != 138)) {
        quallev0 = quallev0 | A211Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      if(aiasen != 24) {
        quallev0 = quallev0 | A211Mech_Err;
      }
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A211Mech_Err;
      }
    }
    if(aiasen == 24) {
      if((aifwen == 203) || (aifwen == 204)) {
        strcpy(wave_str, "211_THIN");
      } else if((aifwen == 137) || (aifwen == 138)) {
        strcpy(wave_str, "211_THICK");
      } else if((aifwen == 74) || (aifwen == 75)) {
        strcpy(wave_str, "211_OPEN");
      }
    }
    else {
      if((aifwen == 203) || (aifwen == 204)) {
        strcpy(wave_str, "MIX_THIN");
      } else if((aifwen == 137) || (aifwen == 138)) {
        strcpy(wave_str, "MIX_THICK");
      } else if((aifwen == 74) || (aifwen == 75)) {
        strcpy(wave_str, "MIX_OPEN");
      }
    }
    break;
  case 8:			//30.4
    if(aifiltyp == 0) {
      if((aifwen != 203) && (aifwen != 204)) {
        quallev0 = quallev0 | A304Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if((aifwen != 137) && (aifwen != 138)) {
        quallev0 = quallev0 | A304Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A304Mech_Err;
      }
    }
    if((aifwen == 203) || (aifwen == 204)) {
      strcpy(wave_str, "304_THIN");
    } else if((aifwen == 137) || (aifwen == 138)) {
      strcpy(wave_str, "304_THICK");
    } else if((aifwen == 74) || (aifwen == 75)) {
      strcpy(wave_str, "304_OPEN");
    }
    break;
  case 0:			//33.5
    if(aifiltyp == 0) {
      if((aifwen != 203) && (aifwen != 204)) {
        quallev0 = quallev0 | A335Mech_Err;
      }
    }
    else if(aifiltyp == 1) {
      if((aifwen != 137) && (aifwen != 138)) {
        quallev0 = quallev0 | A335Mech_Err;
      }
    }
    else if(aifiltyp == 2) {
      if((aifwen != 74) && (aifwen != 75)) {
        quallev0 = quallev0 | A335Mech_Err;
      }
    }
    if((aifwen == 203) || (aifwen == 204)) {
      strcpy(wave_str, "335_THIN");
    } else if((aifwen == 137) || (aifwen == 138)) {
      strcpy(wave_str, "335_THICK");
    } else if((aifwen == 74) || (aifwen == 75)) {
      strcpy(wave_str, "335_OPEN");
    }
    break;
  case 4:			//160.0
    if((aifwen != 269) && (aifwen != 270)) {
      quallev0 = quallev0 | A160Mech_Err;
    }
    else strcpy(wave_str, "1600");
    break;
  case 5:			//170.0
    if((aifwen != 137) && (aifwen != 138)) {
      quallev0 = quallev0 | A170Mech_Err;
    }
    else strcpy(wave_str, "1700");
    break;
  case 6:			//450.0
    if((aifwen != 74) && (aifwen != 75)) {
      quallev0 = quallev0 | A450Mech_Err;
    }
    else strcpy(wave_str, "4500");
    break;
  }
  if(!strcmp(wave_str, "UNKNOWN")) quallev0 = quallev0 | AQ_INVAL_WL;
}
  //drms_setkey_int(rs, "QUALLEV0", quallev0);	//don't use this anymore
  drms_setkey_int(rs, "QUALITY", quallev0);
  drms_setkey_string(rs, "BLD_VERS", bld_vers); //build vers to every record
  drms_setkey_string(rs, "WAVE_STR", wave_str);
  percentd = (float)((100.0 * (float)img->datavals)/(float)img->totalvals);
  drms_setkey_float(rs, "PERCENTD", percentd);
}

// Close out an image.
void close_image(DRMS_Record_t *rs, DRMS_Segment_t *seg, DRMS_Array_t *array,
		IMG *img, int fsn)
{
  STAT stat;
  int status, n, k;
  uint32_t missvals;
  long long cmdx;
  char tlmdsname[128];

  printk("*Closing image for fsn = %u\n", fsn);
  if(imgstat(img, &stat)) {
    printk("**Error on imgstat() for fsn = %u\n", fsn);
  }
  else {
    drms_setkey_short(rs, "DATAMIN", stat.min);
    drms_setkey_short(rs, "DATAMAX", stat.max);
    drms_setkey_short(rs, "DATAMEDN", stat.median);
    drms_setkey_float(rs, "DATAMEAN", stat.mean);
    drms_setkey_float(rs, "DATARMS", stat.rms);
    drms_setkey_float(rs, "DATASKEW", stat.skew);
    drms_setkey_float(rs, "DATAKURT", stat.kurt);
  }
  // CAMERA set by HMI_compute_exposure_times()
  if(hmiaiaflg) {		// except for AIA use telnum 
    drms_setkey_int(rs, "CAMERA", (img->telnum)+1);
    cmdx = drms_getkey_longlong(rs, "AIMGSHCE", &status);
    if (0 == cmdx) drms_setkey_int(rs, "AIMGTYP", 1);
  }
  drms_setkey_int(rs, "IMGAPID", img->apid);
  drms_setkey_int(rs, "CROPID", img->cropid);
  if(drms_setkey_int(rs, "LUTID", img->luid)) {
    printk("ERROR on setkey_int for LUTID\n");
  }
  drms_setkey_int(rs, "TAPCODE", img->tap);
  n = (img->N) & 0x1F;
  n = n << 3;
  k = (img->K) & 0x07;
  if(img->N == 16) { 
    drms_setkey_int(rs, "COMPID", 0);
    drms_setkey_int(rs, "BITSELID", 0);
  }
  else { 
    drms_setkey_int(rs, "COMPID", n+k);
    drms_setkey_int(rs, "BITSELID", img->R);
  }
  drms_setkey_int(rs, "FID", img->fid);
  drms_setkey_int(rs, "TOTVALS", img->totalvals);
  drms_setkey_int(rs, "DATAVALS", img->datavals);
  drms_setkey_int(rs, "NPACKETS", img->npackets);
  drms_setkey_int(rs, "NERRORS", img->nerrors);
  drms_setkey_short(rs, "EOIERROR", img->last_pix_err);
  drms_setkey_short(rs, "HEADRERR", img->headerr);
  drms_setkey_short(rs, "OVERFLOW", img->overflow);
  missvals = img->totalvals - img->datavals;
  drms_setkey_int(rs, "MISSVALS", missvals);
  snprintf(tlmdsname, 128, "%s[%s]", tlmseriesname, tlmnamekeyfirst);
  drms_setkey_string(rs, "TLMDSNAM", tlmdsname);
  unsigned int pts = img->first_packet_time >> 16;
  int ptss = img->first_packet_time & 0xffff;
  TIME fpt = SDO_to_DRMS_time(pts, ptss);
  drms_setkey_double(rs, "IMGFPT", fpt);
  drms_setkey_double(rs, "DATE", CURRENT_SYSTEM_TIME);
  do_quallev0(rs, img, fsn);		//set the QUALLEV0 keyword

  status = drms_segment_write(seg, array, 0);
  if (status) {
    printk("ERROR: drms_segment_write error=%d for fsn=%u\n", status, fsn);
  }
  add_small_array(rs, array, 8, 16); //add Phil's png and small fits
  array->data = NULL;        // must do before free 
  drms_free_array(array);
  if((status = drms_close_record(rs, DRMS_INSERT_RECORD))) {
    printk("**ERROR: drms_close_record failed for %s fsn=%u\n",
                        lev0seriesname, fsn);
  }
  img->initialized = 0;		//indicate image is ready for use again
  //imgdecode_init_hack(img);
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit ingest_lev0 w/ status = %d\n", stat);
  if (h0logfp) fclose(h0logfp);
  exit(stat);
}


int compare_names(const void *a, const void *b)
{
  NAMESORT *x=(NAMESORT *)a, *y=(NAMESORT *)b;
  return(strcmp(x->name+4, y->name+4)); // skip VC02/VC05 in compare 
}

unsigned short MDI_getshort (unsigned char *c)    //  machine independent  
{
  unsigned short s = 0;

  s = (unsigned short) *c++ << 8;
  s |= (unsigned short) *c;
  return s;
}

//fsn_normal_new_image()
int fsn_normal_new_image() 
{
  int dstatus, i;

  // start a new image 
  Img->initialized = 0;
  Img->reopened = 0;
  errmsgcnt = 0;
  fsn_prev = fsnx;
  sprintf(tlmnamekeyfirst, "%s", tlmnamekey);	//save for TLMDSNAM
  rs = drms_create_record(drms_env, lev0seriesname, DRMS_PERMANENT, &dstatus);
  if(dstatus) {
    if(dstatus == DRMS_ERROR_SUMOPEN) {
      printk("**ERROR: DRMS can't open w/SUMS. Aborting...\n");
      abortit(4);
    }
    printk("**ERROR: Can't create record for %s fsn=%u\n", lev0seriesname, fsnx);
    return(1);
  }
  dstatus = drms_setkey_int(rs, "FSN", fsnx);
  if(!(segment = drms_segment_lookup(rs, "image"))) {
    printk("No drms_segment_lookup(rs, image)\n");
    return(1);
  }
  //must initialize Img in case no data pkt comes in for this image
#if 0
    Img->datavals = 0;
    Img->npackets = 0;
    Img->nerrors = 0;
    Img->last_pix_err = 0;
    Img->first_packet_time = UINT64_MAX;
    for (i = 0; i < MAXPIXELS; ++i)
        Img->dat[i] = BLANK;
    for (i = 0; i < MAXHIST; ++i)
        Img->hist[i] = 0;
#endif
    imgdecode_init_hack(Img);

  rdat = Img->dat;
  segArray = drms_array_create(DRMS_TYPE_SHORT,
                                       segment->info->naxis,
                                       segment->axis,
                                       rdat,
                                       &dstatus);
  ImgC = Img;		//set current image pointer
  return(0);
} 

// An fsn has changed in the context of a normal stream of
// tlm files. Returns 0 on success.
int fsn_change_normal()
{
  int dstatus, rstatus, compid, n, k;
  char reopen_dsname[256];
  char *cptr;

  printk("*FSN has changed from %u to %u\n", fsn_prev, fsnx);
  if(fsn_prev == 0) {	//startup mode. restore any prev image
    errmsgcnt = 0;
    fsn_prev = fsnx;
    if(hmiaiaflg)
      sprintf(reopen_dsname, "%s[%u]", LEV0SERIESNAMEAIA, fsnx);
    else
      sprintf(reopen_dsname, "%s[%u]", LEV0SERIESNAMEHMI, fsnx);
    printk("Open normal prev ds: %s\n", reopen_dsname);
    rset = drms_open_records(drms_env, reopen_dsname, &rstatus);
    if(rstatus) {
      printk("Can't do drms_open_records(%s)\n", reopen_dsname);
      return(1);          // !!!TBD
    }
    if(!rset || (rset->n == 0) || rstatus) {
      printk("No prev ds\n");     // start a new image
      fsn_normal_new_image();
    }
    else {
      Img->initialized = 1;
      Img->reopened = 1;
      Img->fsn = fsnx;
      Img->apid = appid;
      rs_old = rset->records[0];
      rs = drms_clone_record(rs_old, DRMS_PERMANENT, DRMS_COPY_SEGMENTS, &rstatus);
      if(rstatus) {
        printk("Can't do drms_clone_record()\n");
        return(1);		// !!!TBD ck 
      }
      drms_close_records(rset, DRMS_FREE_RECORD);
      rstatus = drms_setkey_int(rs, "FSN", fsnx);
      Img->telnum = drms_getkey_int(rs, "CAMERA", &rstatus);
      Img->telnum--;
      Img->cropid = drms_getkey_int(rs, "CROPID", &rstatus);
      Img->fid = drms_getkey_int(rs, "FID", &rstatus); //!!NEW 22Mar2010
      Img->luid = drms_getkey_int(rs, "LUTID", &rstatus);
      if(rstatus) {
        Img->luid = 0;	//!!!TEMP try this
        printk("ERROR on getkey_int for LUTID\n");
      }
      Img->tap = drms_getkey_int(rs, "TAPCODE", &rstatus);
      compid = drms_getkey_int(rs, "COMPID", &rstatus);
      if(compid == 0) {
        n = 16; k = 0;
      }
      else {
        k = compid & 0x07;		//low 3 bits
        n = compid >> 3;			//next 5 bits
        n = n & 0x1F;
      }
      Img->N = n;
      Img->K = k;
      Img->R = drms_getkey_int(rs, "BITSELID", &rstatus);
      Img->overflow = drms_getkey_int(rs, "OVERFLOW", &rstatus);
      Img->headerr = drms_getkey_int(rs, "HEADRERR", &rstatus);
      Img->totalvals = drms_getkey_int(rs, "TOTVALS", &rstatus);
      Img->datavals = drms_getkey_int(rs, "DATAVALS", &rstatus);
      Img->npackets = drms_getkey_int(rs, "NPACKETS", &rstatus);
      Img->nerrors = drms_getkey_int(rs, "NERRORS", &rstatus);
      Img->last_pix_err = drms_getkey_int(rs, "EOIERROR", &rstatus);

      TIME fpt = drms_getkey_double(rs, "IMGFPT", &rstatus);
      Img->first_packet_time = round(65536.0*(fpt - sdo_epoch));
      if(Img->first_packet_time == 0) 
        Img->first_packet_time = UINT64_MAX; //fix for stuck 1958 value
      snprintf(oldtlmdsnam, 128, "%s", drms_getkey_string(rs, "TLMDSNAM",&rstatus));
      cptr = strstr(oldtlmdsnam, "VC");
      if(cptr) 
        snprintf(tlmnamekeyfirst, 23, "%s", cptr);
      else
        sprintf(tlmnamekeyfirst, "%s", "UNK");
      segment = drms_segment_lookupnum(rs, 0);
      cArray = drms_segment_read(segment, DRMS_TYPE_SHORT, &rstatus);
      if(rstatus) {
        printk("Can't do drms_segment_read()\n");
        return(1);		// !!!!TBD ck 
      }
      short *adata = (short *)cArray->data;
      memcpy(Img->dat, adata, 2*MAXPIXELS);
      rdat = Img->dat;
      segArray = drms_array_create(DRMS_TYPE_SHORT, 
                                         segment->info->naxis,
                                         segment->axis,
                                         rdat,
                                         &dstatus);
      ImgC = Img;           //set current image pointer
      drms_free_array(cArray); //must free from drms_segment_read()
    }
  }
  else {
    if(rs) {		// make sure have a created record
      close_image(rs, segment, segArray, Img, fsn_prev);
      //force a commit to DB on an image boundry & not during a rexmit
      if(imagecnt++ >= IMAGE_NUM_COMMIT) {
        printk("drms_server_end_transaction()\n");
        drms_server_end_transaction(drms_env, 0 , 0);
        printk("drms_server_begin_transaction()\n");
        drms_server_begin_transaction(drms_env); //start another cycle
        imagecnt = 0;
      }
      fileimgcnt++;
    }
    else {
      printk("**ERROR: Null record ptr for an opened image fsn=%u\n",fsn_prev);
    }
    if(fsn_normal_new_image()) {
      printk("**ERROR: Can't create new image\n");
      return(1);
    }
  }
  return(0);
}

// An fsn has changed in the context of a retransmitted or higher version
// tlm file. Returns 0 on success.
int fsn_change_rexmit() 
{
  char rexmit_dsname[256];
  int rstatus, dstatus, compid, n, k;
  char *cptr;

  if(fsn_prev != 0) {   // close image of prev fsn if not 0 
    if(rsc) {		//make sure have created record
      close_image(rsc, segmentc, oldArray, ImgO, fsn_prev);
      imagecnt++;
      fileimgcnt++;
    }
    else {
      printk("**ERROR: Null record ptr for an rexmit image fsn=%u\n",fsn_prev);
    }
  }
  errmsgcnt = 0;
  fsn_prev = fsnx;
  if(hmiaiaflg) 
    sprintf(rexmit_dsname, "%s[%u]", LEV0SERIESNAMEAIA, fsnx);
  else 
    sprintf(rexmit_dsname, "%s[%u]", LEV0SERIESNAMEHMI, fsnx);
  printk("Open prev ds: %s\n", rexmit_dsname);
  rset = drms_open_records(drms_env, rexmit_dsname, &rstatus); 
  if(rstatus) {
    printk("Can't do drms_open_records(%s)\n", rexmit_dsname);
    return(1);		// !!!TBD 
  }
  if(!rset || (rset->n == 0) || rstatus) {
startnew:
    printk("No prev ds\n");	// start a new image 
    ImgO->initialized = 0;
    ImgO->reopened = 0;
    imgdecode_init_hack(ImgO);
    sprintf(tlmnamekeyfirst, "%s", tlmnamekey);	//save for TLMDSNAM
    rsc = drms_create_record(drms_env, lev0seriesname, DRMS_PERMANENT, &rstatus);
    if(rstatus) {
      if(rstatus == DRMS_ERROR_SUMOPEN) {
        printk("**ERROR: DRMS can't open w/SUMS. Aborting...\n");
        abortit(4);
      }
      printk("Can't create record for %s fsn=%u\n", lev0seriesname, fsnx);
      return(1);                     // !!!TBD ck this 
    }
    rstatus = drms_setkey_int(rsc, "FSN", fsnx);
    if(!(segmentc = drms_segment_lookup(rsc, "image"))) {
      printk("No drms_segment_lookup(rsc, image)\n");
      return(1);
    }
    rdat = ImgO->dat;
    oldArray = drms_array_create(DRMS_TYPE_SHORT,
                                       segmentc->info->naxis,
                                       segmentc->axis,
                                       rdat,
                                       &dstatus);
  }
  else {
    ImgO->initialized = 1;
    ImgO->reopened = 1;
    ImgO->fsn = fsnx;
    ImgO->apid = appid;
    rs_old = rset->records[0];
    rsc = drms_clone_record(rs_old, DRMS_PERMANENT, DRMS_COPY_SEGMENTS, &rstatus);
    if(rstatus) {
      printk("Can't do drms_clone_record()\n");
      printk("Assume this was a temp ds and the segments are now gone...\n");
      goto startnew;
      //return(1);		// !!!TBD ck 
    }
    drms_close_records(rset, DRMS_FREE_RECORD);
    rstatus = drms_setkey_int(rsc, "FSN", fsnx);
    ImgO->telnum = drms_getkey_int(rsc, "CAMERA", &rstatus);
    ImgO->telnum--;
    ImgO->cropid = drms_getkey_int(rsc, "CROPID", &rstatus);
    ImgO->fid = drms_getkey_int(rsc, "FID", &rstatus); //!!NEW 22Mar2010
    ImgO->luid = drms_getkey_int(rsc, "LUTID", &rstatus);
    ImgO->tap = drms_getkey_int(rsc, "TAPCODE", &rstatus);
    compid = drms_getkey_int(rsc, "COMPID", &rstatus);
    if(compid == 0) {
      n = 16; k = 0;
    }
    else {
      k = compid & 0x07;              //low 3 bits
      n = compid >> 3;                        //next 5 bits
      n = n & 0x1F;
    }
    ImgO->N = n;
    ImgO->K = k;
    ImgO->R = drms_getkey_int(rsc, "BITSELID", &rstatus);
    ImgO->overflow = drms_getkey_int(rsc, "OVERFLOW", &rstatus);
    ImgO->headerr = drms_getkey_int(rsc, "HEADRERR", &rstatus);
    ImgO->totalvals = drms_getkey_int(rsc, "TOTVALS", &rstatus);
    ImgO->datavals = drms_getkey_int(rsc, "DATAVALS", &rstatus);
    ImgO->npackets = drms_getkey_int(rsc, "NPACKETS", &rstatus);
    ImgO->nerrors = drms_getkey_int(rsc, "NERRORS", &rstatus);
    ImgO->last_pix_err = drms_getkey_int(rsc, "EOIERROR", &rstatus);

    TIME fpt = drms_getkey_double(rsc, "IMGFPT", &rstatus);
    ImgO->first_packet_time = round(65536.0*(fpt - sdo_epoch));
    if(ImgO->first_packet_time == 0) 
      ImgO->first_packet_time = UINT64_MAX; //fix for stuck 1958 value
    snprintf(oldtlmdsnam, 128, "%s", drms_getkey_string(rsc, "TLMDSNAM",&rstatus));
    cptr = strstr(oldtlmdsnam, "VC");
    if(cptr) 
      snprintf(tlmnamekeyfirst, 23, "%s", cptr);
    else
      sprintf(tlmnamekeyfirst, "%s", "UNK");
    segmentc = drms_segment_lookupnum(rsc, 0);
    cArray = drms_segment_read(segmentc, DRMS_TYPE_SHORT, &rstatus);
    if(rstatus) {
      printk("Can't do drms_segment_read()\n");
      return(1);		// !!!!TBD ck 
    }
    short *adata = (short *)cArray->data;
    memcpy(ImgO->dat, adata, 2*MAXPIXELS);
    rdat = ImgO->dat;
    oldArray = drms_array_create(DRMS_TYPE_SHORT, 
                                       segmentc->info->naxis,
                                       segmentc->axis,
                                       rdat,
                                       &dstatus);
    drms_free_array(cArray); //must free from drms_segment_read()
  }
  ImgC = ImgO;		//set current image pointer
  return(0);
}

// Process the tlm file to validate and to extract the lev0 image.
int get_tlm(char *file, int rexmit, int higherver)
{
  FILE *fpin;
  unsigned char cbuf[PKTSZ];
  long long gap_42_cnt;
  int status, rstatus, fpkt_cnt, i, j, sync_bad_cnt;
  int datval, eflg, firstflg;
  unsigned int cnt1, cnt2, cnt3, gap_24_cnt;
  int zero_pn;
  unsigned short pksync1, pksync2;
  float ftmp;
  int decode_status=0;
  unsigned int Fsn;

  if(!(fpin = fopen(file, "r"))) {	// open the tlm input 
    printk("*Can't open tlm file %s\n", file);
    return(1);
  }
  BeginTimer(1);			//time tlm file processing
  printk("*Processing tlm file %s\n", file);
  fpkt_cnt = sync_bad_cnt = eflg = errmsgcnt = fileimgcnt = 0;
  zero_pn = gap_24_cnt = gap_42_cnt = 0;
  firstflg = 1; 
  if(rexmit || higherver) {
    if(Img->initialized) {		// close normal mode image
      if(rs) {
        close_image(rs, segment, segArray, Img, fsn_prev);
        printk("drms_server_end_transaction()\n");
        drms_server_end_transaction(drms_env, 0 , 0);
        printk("drms_server_begin_transaction()\n");
        drms_server_begin_transaction(drms_env); //start another cycle
        imagecnt = 0;
      }
      else
        printk("**ERROR: Null record ptr for an opened image fsn=%u\n",fsn_prev);
    }
    fsn_pre_rexmit = fsn_prev;		// restore this at end of rexmit tlm file
    fsn_prev = 0;			// cause a new image 
  }
  // read a VCDU packet 
  while((status = fread(cbuf,sizeof(char),PKTSZ,fpin) ) == PKTSZ) {
    pksync1 = MDI_getshort(cbuf);
    pksync2 = MDI_getshort(cbuf+2);
    if((pksync1 == 0) && (pksync2 == 0)) { // skip 0 pn code 
      if(!zero_pn) {			// give msg for 1st one only 
        printk("*0 PN code at pkt# %d\n", fpkt_cnt);
        printk("*Subsequent ones will be ignored until non-0 again\n");
        zero_pn = 1;
      }
      fpkt_cnt++;			// count # of pkts found 
      continue;
    }
    if((pksync1 != 0x1acf) || (pksync2 != 0xfc1d)) {
      printk("*Lost sync at VCDU pkt# %d. pksync1=%x pksync2=%x\n", 
		fpkt_cnt, pksync1, pksync2);
      fpkt_cnt++;			// count # of pkts found 
      eflg++;
      if(sync_bad_cnt++ > 4) {
        printk("**Too many out of sync packets.\n");
        //return(1);
        return(0);		//changed on 01Oct2010
      }
      printk("  Will attempt to press on...\n");
      zero_pn = 0;
      continue;
    }
    if(firstflg) {		// print first good sync found 
      printk("*VCDU pkt# %d sync = %x %x\n", fpkt_cnt, pksync1, pksync2);
    }
    fpkt_cnt++;			// count # of pkts found 
    // get 24 bit VCDU counter 
    cnt1 = MDI_getshort(cbuf+6);
    cnt2 = MDI_getshort(cbuf+8);
    cnt2 = (cnt2 >> 8)& 0xFF;
    cnt2 = ((cnt1 << 8)& 0xFF00) + cnt2;
    cnt1 = (cnt1 >> 8)& 0xFF;
    vcdu_24_cnt = (cnt1*65536) + cnt2;
    if(vcdu_24_cnt_next != vcdu_24_cnt) {
      printk("*VCDU 24bit seq num out of sequence. exp: %u  rec: %u\n", 
	    vcdu_24_cnt_next, vcdu_24_cnt);
      if(vcdu_24_cnt_next > vcdu_24_cnt) {
        printk("*NOTE: VCDU 24 bit counter retarded\n"); //cntr does go thru 0
        printk("*NOTE: gap report will be inaccurate (tbd)\n");
      }
      if(!firstflg) {		// don't count gap across .tlm files 
        gap_24_cnt += vcdu_24_cnt - vcdu_24_cnt_next;
      }
    }
    vcdu_24_cnt_next = vcdu_24_cnt + 1;
    // now get the 42bit IM_PDU counter 
    cnt1 = MDI_getshort(cbuf+10);
    cnt2 = MDI_getshort(cbuf+12);
    cnt3 = MDI_getshort(cbuf+14);
    cnt1 = cnt1 & 0x03ff;
    vcdu_seq_num = (cnt1*4294967296) + (cnt2*65536) + cnt3;
    // printk("vcdu_seq_num = %lld\n", vcdu_seq_num); 
    if(vcdu_seq_num_next != vcdu_seq_num) {
      printk("*IM_PDU seq num out of sequence. exp: %lld  rec: %lld\n", 
	    vcdu_seq_num_next, vcdu_seq_num);
      if(vcdu_seq_num_next > vcdu_seq_num) {
        printk("*NOTE: IM_PDU 42 bit counter retarded\n");
        printk("*NOTE: gap report will be inaccurate\n");
      }
      if(!firstflg) {		// don't count gap across .tlm files 
        gap_42_cnt += vcdu_seq_num - vcdu_seq_num_next;
      }
      eflg++;
    }
    firstflg = 0;
    vcdu_seq_num_next = vcdu_seq_num + 1;
    // get the App ID. Low 11 bit of short at buf+18 
    appid = MDI_getshort(cbuf+18);
    appid = appid & 0x07ff;
    if(hmiaiaflg) {		//aia
      //testid1 = 509; testid2 = 519;
      testid1 = 507; testid2 = 517;
    }
    else {
      //testid1 = 409; testid2 = 419;
      testid1 = 407; testid2 = 417;
    }
    //if(appid == TESTAPPID) {	// appid of test pattern 
    if((appid == testid1) || (appid == testid2)) {
      if(errmsgcnt++ < MAXERRMSGCNT) {
        printk("*Test ApID of %0x found for IM_PDU Cntr = %lld\n", 
			appid, vcdu_seq_num);
        for(i=0, j=TESTVALUE; i < 877; i=i+2, j++) {
          datval = MDI_getshort(cbuf+32+i);	// next data value 
          if(datval != j) {
            printk("*Test data value=%0x, expected=%0x for IM_PDU Cntr=%lld\n", 
		datval, j, vcdu_seq_num);
            printk("*File = %s\n", file);
            eflg++;
            break;		// skip the rest of this packet 
          }
        }
      }
      continue; 		// go on to next packet 
    }

    //printk("$$$$$ appid found =  0x%x %d\n", appid, appid); //!!TEMP
    // Parse tlm packet headers. 
    if(appid == APID_HMI_SCIENCE_1 || appid == APID_HMI_SCIENCE_2 || 
	appid == APID_AIA_SCIENCE_1 || appid == APID_AIA_SCIENCE_2)
    {
      cnt1 = MDI_getshort(cbuf+32);
      cnt2 = MDI_getshort(cbuf+34);
      fsnx = (unsigned int)(cnt1<<16)+(unsigned int)(cnt2);
      fsnx = fsnx & 0x3fffffff;		//low 30bits for fsn */
      if(fsnx == 0) continue;		//a 0 fsn is not acceptable
      if(rexmit || higherver) {
        if(fsnx != fsn_prev) {          // the fsn has changed
          if(fsn_change_rexmit()) {	//handle old & new images
            printk("***FATAL ERROR in fsn_change_rexmit()\n");
            return(1);
          }
        }
        else {				//prev could be hk w/diff apid
          ImgC->apid = appid;
        }
      }
      else {			// continuing normal stream
        if(fsnx != fsn_prev) {          // the fsn has changed 
          if(fsn_change_normal()) {	//handle old & new images
            printk("***FATAL ERROR in fsn_change_normal()\n");
            return(1);
          }
        }
      }
      // send the sci data to Keh-Cheng. call with pointer to M_PDU_Header 
      rstatus = imgdecode((unsigned short *)(cbuf+10), ImgC);
      switch(rstatus) {
      case 0:
        // A science data VCDU was successfully decoded 
        break;
      case IMGDECODE_DECOMPRESS_ERROR:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_DECOMPRESS_ERROR\n");
        break;
      case IMGDECODE_TOO_MANY_PIXELS:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_TOO_MANY_PIXELS\n");
        break;
      case IMGDECODE_BAD_N:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_N\n");
        break;
      case IMGDECODE_BAD_APID:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_APID\n");
        break;
      case IMGDECODE_NO_LOOKUP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_NO_LOOKUP_TABLE\n");
        break;
      case IMGDECODE_LOOKUP_ID_MISMATCH:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_LOOKUP_ID_MISMATCH\n");
        break;
      case IMGDECODE_BAD_LOOKUP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_LOOKUP_TABLE\n");
        break;
      case IMGDECODE_NO_CROP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_NO_CROP_TABLE\n");
        break;
      case IMGDECODE_CROP_ID_MISMATCH:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_CROP_ID_MISMATCH\n");
        break;
      case IMGDECODE_BAD_CROP_GEOMETRY:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_CROP_GEOMETRY\n");
        break;
      case IMGDECODE_BAD_CROP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_CROP_TABLE\n");
        break;
      case IMGDECODE_BAD_CROP_SKIP_TAKE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_CROP_SKIP_TAKE\n");
        break;
      case IMGDECODE_BAD_OFFSET:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_BAD_OFFSET\n");
        break;
      case IMGDECODE_OUT_OF_MEMORY:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: IMGDECODE_OUT_OF_MEMORY\n");
        break;
      default:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode() ret: unknown err status = %d:\n", rstatus);
        break;
      }
    }
    else {			// send the HK data to Carl 
      //printk("$$$ appid assumed hk =  0x%x %d\n", appid, appid); //!!TEMP
      decode_status = decode_next_hk_vcdu((unsigned short *)(cbuf+10), &Hk, &Fsn);
      switch (decode_status) {
        case SUCCESS_HK_NEED_TO_WTD_CTD:
          printk("*ISP found for fsn = %u\n", Fsn);
          fsnx = Fsn;
          if(rexmit || higherver) {
            if(fsnx != fsn_prev) {          // the fsn has changed
              if(fsn_change_rexmit()) {	//handle old & new images
                printk("***FATAL ERROR in fsn_change_rexmit()\n");
                return(1);
              }
            }
            //calculate and setkey some values from the keywords returned
            HMI_compute_exposure_times(rsc, Hk->keywords, hmiaiaflg);
            whk_status = write_hk_to_drms(rsc, &Hk); //ISP keywords to drms
            if(!hmiaiaflg) {
              if(whk_status = set_HMI_mech_values(rsc)) {
                printk("***ERROR: mechanism position keywords are wrong!\n");
              }
            }
          }
          else {					//normal stream
            if(fsnx != fsn_prev) {          // the fsn has changed 
              if(fsn_change_normal()) {	//handle old & new images
                printk("***FATAL ERROR in fsn_change_normal()\n");
                return(1);
              }
            }
            //calculate and setkey some values from the keywords returned
            HMI_compute_exposure_times(rs, Hk->keywords, hmiaiaflg);
            whk_status = write_hk_to_drms(rs, &Hk); //ISP keywords to drms
            if(!hmiaiaflg) {
              if(whk_status = set_HMI_mech_values(rs)) {
                printk("***ERROR: mechanism position keywords are wrong!\n");
              }
            }
          }
          hk_ccsds_free(&Hk);
          rstatus=0;
          break;
        case SUCCESS_HK_NEED_TO_CTD:
          /*NO need to write  keywords's to drms for level0 data series 
            BUT need to commit to drms for level0 by APID data series*/
          break;
        case SUCCESS_SKIP_IMAGE:
          printk("decode_next_hk_vcdu() ret:Warning:SUCCESS_SKIP_IMAGE\n");
          break;
        case SUCCESS_SKIP_PROCESSING_APID:
          printk("decode_next_hk_vcdu() ret:Warning:SUCCESS_SKIP_PROCESSING_APID:\n");
          break;
        case ERROR_HK_FAILED_GETTING_FSN:
          printk("decode_next_hk_vcdu() ret:Warning:ERROR_HK_FAILED_GETTING_FSN:\n");
          break;
        case ERROR_NODATA:
          printk("decode_next_hk_vcdu() ret:Warning:ERROR_NODATA\n");
          break;
        case ERROR_HK_ENVIRONMENT_VARS_NOT_SET:
          printk("decode_next_hk_vcdu() ret:ERROR:ERROR_HK_ENVIRONMENT_VARS_NOT_SET\n");
           break;
        case ERROR_HK_FAILED_OPEN_DAYFILE:
          printk("decode_next_hk_vcdu() ret:ERROR:ERROR_HK_FAILED_OPEN_DAYFILE\n");
           break;
        default:
          printk("decode_next_hk_vcdu() ret:Warning:Unexpected return code for decode_status:<%d>:\n", decode_status);
          break;
      }

    }
  }				// end of rd of vcdu pkts 
  fclose(fpin);
  if(!eflg) {
    printk("*No errors in tlm file\n");
  }
  if(rexmit || higherver) {	// close the opened record 
    if(rsc) {			// make sure have created record
      close_image(rsc, segmentc, oldArray, ImgO, fsnx);
      printk("drms_server_end_transaction()\n");  //commit at end of file
      drms_server_end_transaction(drms_env, 0 , 0);
      printk("drms_server_begin_transaction()\n\n");
      drms_server_begin_transaction(drms_env); //start another cycle
      imagecnt = 0;
      fileimgcnt++;
      //fsn_prev = fsn_pre_rexmit; // restore orig for next normal .tlm file
      fsn_prev = 0;		//need for restart of normal mode
      rsc = 0;
    }
  }
  ftmp = EndTimer(1);
  printk("**Processed %s\n**with %d images and %d VCDUs in %f sec\n\n",
	file, fileimgcnt, fpkt_cnt, ftmp);

  //now commit to DB in case another file doesn't come in
/***************************TEMP***********************************
NOTE: Can't do this unless close current image and then set fsn=0
      so that the next tlm file will start by opening the image again.
      I don't think that we want this. Instead find a way to get the ^C
      instead of DRMS, and then do a close and commit.
  printk("drms_server_end_transaction()\n");
  drms_server_end_transaction(drms_env, 0 , 0);
  printk("drms_server_begin_transaction()\n\n");
  drms_server_begin_transaction(drms_env); //start another cycle
  imagecnt = 0;
***************************TEMP***********************************/

  if(fpkt_cnt != total_tlm_vcdu) {
    printk("**WARNING: Found #vcdu=%d; expected=%d\n", fpkt_cnt, total_tlm_vcdu);
  }
  //if(gap_24_cnt != total_missing_vcdu) {
  //  printk("**WARNING: VCDU 24bit cntr gaps=%d; expected=%d\n",
  //      gap_24_cnt, total_missing_vcdu);
  //}
  if(gap_42_cnt != total_missing_im_pdu) {
    printk("**WARNING: IM_PDU 42bit cntr gaps=%lld; expected=%lld\n",
	 gap_42_cnt, total_missing_im_pdu);
  }
  ignoresigalrmflg = 1;		//got a file. noop next alarm signal timeout
  return(0);
}

// This is the main loop that gets the .qac and .tlm files and 
// puts them into ds tlmseriesname and extracts the lev0 and puts it in 
// lev0seriesname in DRMS.
void do_ingest()
{
  FILE *fp;
  DRMS_Record_t *rs_tlm;
  DIR *dfd;
  NAMESORT *nameptr;
  struct dirent *dp;
  int i, j, status, mvstat;
  int rexmit, higherversion;
  char path[DRMS_MAXPATHLEN];
  char name[128], line[128], tlmfile[128], tlmname[96];
  char cmd[256], xxname[128], vername[16];
  char *token, *filename;

  if((dfd=opendir(tlmdir)) == NULL) {
    printk("**Can't opendir(%s) to find files\n", tlmdir);
    abortit(3);
  }
  i = 0;
  if((nameptr = (NAMESORT *)malloc(MAXFILES * sizeof(NAMESORT))) == NULL) {
    printk("***Can't alloc memory for file name sort\n");
    abortit(3);
  }

  while((dp=readdir(dfd)) != NULL) {
    // printk("%s\n", dp->d_name) ; continue; // !!TEMP 
    // Only accept our files. 
    if(strstr(dp->d_name, pchan) || strstr(dp->d_name, rchan) || strstr(dp->d_name, ".dsf")) {
      nameptr[i++].name = strdup(dp->d_name);
      if(i >= MAXFILES) {
        printk("***Fatal error. Too many (%d) files in %s\n", MAXFILES, tlmdir);
        printf("***Fatal error. Too many (%d) files in %s\n", MAXFILES, tlmdir);
        abortit(3);
      }
      if(!strstr(dp->d_name, ".dsf")) cntsleeps = 0;	//we saw a file
    }
  }
  closedir(dfd);
  qsort(nameptr, i, sizeof(NAMESORT), &compare_names);

  for(j=0; j < i; j++) {
    //printk("####QSORT FILES: %s\n", nameptr[j].name); // !!TEMP 
    // NOTE: the dsf files is moved to indir/dsf e.g. (/dds/soc2pipe/aia/dsf)
    // Currently the cron job pipefe_rm does this:
    // `/bin/mv $dsfname /dds/socdc/hmi/dsf`
    //Only do for V01/VC02. VC04/VC05 won't mv the .dsf
    if(!strcmp(pchan, "VC01") || !strcmp(pchan, "VC02")) { 
      if(strstr(nameptr[j].name, ".dsf")) {
        //sprintf(cmd, "/bin/mv %s/%s %s", tlmdir, nameptr[j].name, outdir);
        sprintf(cmd, "/bin/mv -f %s/%s %s/dsf/",tlmdir,nameptr[j].name,tlmdir);
        printk("*mv dsf file to %s/dsf/\n", tlmdir);
        printk("%s\n", cmd);
        if(system(cmd)) {
          printk("***Error on: %s\n", cmd);
        }
      }
    }
    if(!strstr(nameptr[j].name, ".qac")) {	// can be .qac or .qacx 
      free(nameptr[j].name);
      continue;
    }
    rexmit = higherversion = 0;
    if(strstr(nameptr[j].name, ".qacx")) {  	// this is a rexmit file 
      rexmit = 1;
    }
    sprintf(name, "%s/%s", tlmdir, nameptr[j].name);
    printk("%s\n*Found qac file:\n* %s\n", do_datestr(), name);
    if(!(fp=fopen(name, "r"))) {
      printk("***Can't open %s\n", name);
      free(nameptr[j].name);
      continue;
    }
    // NOTE: the qac file is already verified by the caller of ingest_lev0 
    while(fgets(line, 256, fp)) {	// get qac file lines 
      if(line[0] == '#' || line[0] == '\n') continue;
      if(strstr(line, "TLM_FILE_NAME=")) {
        token = (char *)strtok(line, "=");
        token = (char *)strtok(NULL, "\n");
        printk("tlm file is %s\n", token);
        sprintf(tlmfile, "%s/%s", tlmdir, token);
        sprintf(tlmname, "%s", token);
      }
      else if(strstr(line, "TLM_FILE_SIZE=")) {
        token = (char *)strtok(line, "=");
        token = (char *)strtok(NULL, "=");
        printk("*tlm file size is %s", token);
      }
      else if((strstr(line, "TOTAL_TLM_IM_PDU=")) || (strstr(line, "TOTAL_VCDU="))) {
        token = (char *)strtok(line, "=");
        token = (char *)strtok(NULL, "\n");
        total_tlm_vcdu = atoi(token);
      }
      else if(strstr(line, "TOTAL_MISSING_VCDU=")) {
        token = (char *)strtok(line, "=");
        token = (char *)strtok(NULL, "\n");
        total_missing_vcdu = atoi(token);
      }
      else if(strstr(line, "TOTAL_MISSING_IM_PDU=")) {
        token = (char *)strtok(line, "=");
        token = (char *)strtok(NULL, "\n");
        total_missing_im_pdu = atol(token);
        break;
      }
    }
    fclose(fp);
    //alarm(ALRMSEC);		//restart alarm if no more files come
    strcpy((char *)tlmnamekey, (char *)tlmname);
    filename = (char *)rindex(tlmname, '.');
    if(!filename) {
      printk("**ERROR: Bad formed tlm name in qac file\n");
      printk("  tlm name = %s\n", tlmname);
      continue;
    }
    *filename = 0;			// elim .tlm 
    *((char *)tlmnamekey+22) = 0;	// elim after hh_mm_ss
    rs_tlm = drms_create_record(drms_env, tlmseriesname, 
				DRMS_PERMANENT, &status);
    if(status) {
      if(status == DRMS_ERROR_SUMOPEN) {
        printk("**ERROR: DRMS can't open w/SUMS. Aborting...\n");
        abortit(4);
      }
      printk("***Can't create record for %s. Status=%d\n", tlmseriesname, status);
      continue;
    }
    if((status = drms_setkey_string(rs_tlm, "filename", tlmnamekey))) {
      printk("**ERROR: drms_setkey_string failed for 'filename'\n");
    }
    drms_record_directory(rs_tlm, path, 0);
    if(!*path) {
      printk("***ERROR: No path to segment for %s\n", tlmseriesname);
      continue;
    }
    if(outdir) {
      sprintf(cmd, "cp -p %s %s", name, outdir);
      printk("*cp qac file to %s\n", outdir);
      printk("%s\n", cmd);
      if(system(cmd)) {
        printk("***Error on: %s\n", cmd);
      }
      sprintf(cmd, "cp -p %s %s", tlmfile, outdir);
      printk("*cp tlm file to %s\n", outdir);
      printk("%s\n", cmd);
      if(system(cmd)) {
        printk("***Error on: %s\n", cmd);
      }
    }
    sprintf(cmd, "/bin/mv %s %s", name, path);
    printk("*mv qac to %s\n", path);
    printk("%s\n", cmd);
    if(status = system(cmd)) {
      printk("**ERROR: %d errno=%d on: %s\n", status, errno, cmd);
      if(WIFEXITED(status)) {
        if(mvstat = WEXITSTATUS(status)) {  //status ret by mv
          printk("**ERROR: mv exit status = %d\n", mvstat);
        }
      }
    }
    sprintf(cmd, "/bin/mv %s %s", tlmfile, path);
    printk("*mv tlm to %s\n", path);
    printk("%s\n", cmd);
    if(status = system(cmd)) {
      printk("**ERROR: %d errno=%d on: %s\n", status, errno, cmd);
      if(WIFEXITED(status)) {
        if(mvstat = WEXITSTATUS(status)) {  //status ret by mv
          printk("**ERROR: mv exit status = %d\n", mvstat);
        }
      }
      printk("**Continue after ERROR on mv of tlm file\n");
      continue;
    }
    if((status = drms_close_record(rs_tlm, DRMS_INSERT_RECORD))) {
      printk("**ERROR: drms_close_record failed for %s\n", tlmseriesname);
    }

    sprintf(xxname, "%s/%s.tlm", path, tlmname);
    filename = (char *)rindex(tlmname, '_');
    sprintf(vername, "%s", filename);	// e.g. _00 or _01 etc. 
    if(strcmp(vername, "_00")) {	// this is a higher vers # file 
      higherversion = 1;
      printk("Higher version tlm found: %s.tlm\n", tlmname);
    }
    firstfound = 1;			//a file has been seen
    if(get_tlm(xxname, rexmit, higherversion)) { // lev0 extraction of image 
      printk("***Error in lev0 extraction for %s\n", xxname);
      //printk("***Going to abort\n");
      //abortflg = 1;
    }
    if((stat(stopfile, &stbuf) == 0) || abortflg) { break; } //signal to stop
  }
  free(nameptr);
}

// Initial setup stuff called when main is first entered.
void setup()
{
  FILE *fp;
  int i;
  char string[128], cwdbuf[128], idstr[256];
  char envfile[100], s1[256],s2[256],s3[256], line[256];
  ThreadSigErr_t error = kThreadSigErr_Success;

  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  do_datestr();
  printk_set(h0log, h0log);	// set for printk calls 
  printk("%s\n", datestr);
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
  if(grounddata)
    sprintf(string, "ingest_lev0 -g started as pid=%d user=%s\n", getpid(), username);
  else
    sprintf(string, "ingest_lev0 started as pid=%d user=%s\n", getpid(), username);
  strcat(idstr, string);
  printk("*%s", idstr);
  if(restartflg) printk("-r ");
  sprintf(argvc, "vc=%s", vc);
  sprintf(argindir, "indir=%s", tlmdir);
  sprintf(arglogfile, "logfile=%s", logname);
  if(outdir) {
    sprintf(argoutdir, "outdir=%s", outdir);
    printk("%s %s %s %s\n", argvc, argindir, argoutdir, arglogfile);
  } else {
    printk("%s %s %s\n", argvc, argindir, arglogfile);
  }
  strcpy(pchan, vc);		// virtual channel primary 
  if(rexmitmode)		//new 'stopX' 11/2/2012
    sprintf(stopfile, "/usr/local/logs/lev0/%s_stopX", pchan);
  else
    sprintf(stopfile, "/usr/local/logs/lev0/%s_stop", pchan);
  //Dont rm stopfile any more (1/6/2012). 
  //With the new rexmit dir there is a second
  //ingest_lev0 running on a VC and it needs the stop file too.
  //The stop file is removed when ingest_lev0 is started by
  //doingestlev0_[HMI,AIA].pl
  //sprintf(string, "/bin/rm -f %s", stopfile);	//remove any stop file
  //system(string);
  for(i=0; ; i++) {		// ck for valid and get redundant chan 
    if(!strcmp(p_r_chan_pairs[i].pchan, pchan)) {
      strcpy(rchan, p_r_chan_pairs[i].rchan);
      hmiaiaflg = p_r_chan_pairs[i].instru;
      break;
    }
    if(!strcmp(p_r_chan_pairs[i].pchan, "n/a")) {
      printk("!!ERROR: Invalid VCid (%s) specified\n", pchan);
      fprintf(stderr, "!!ERROR: Invalid VCid (%s) specified. Abort\n", pchan);
      abortit(1);
    }
  }
  if(hmiaiaflg) {
    if(grounddata) {
      sprintf(tlmseriesname, "%s", TLMSERIESNAMEAIAGND);
      sprintf(lev0seriesname, "%s", LEV0SERIESNAMEAIAGND);
    }
    else {
      sprintf(tlmseriesname, "%s", TLMSERIESNAMEAIA);
      sprintf(lev0seriesname, "%s", LEV0SERIESNAMEAIA);
    }
  }
  else {
    if(grounddata) {
      sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMIGND);
      sprintf(lev0seriesname, "%s", LEV0SERIESNAMEHMIGND);
    }
    else {
      sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMI);
      sprintf(lev0seriesname, "%s", LEV0SERIESNAMEHMI);
    }
  }
  if(!restartflg) {
    printk("tlmseriesname=%s\nlev0seriesname=%s\n", 
		tlmseriesname, lev0seriesname);
  }
  sprintf(bld_vers, "%s", jsoc_version);
  umask(002);			// allow group write 
  Image.initialized = 0;	// init the two image structures 
  ImageOld.initialized = 0;
  Img = &Image;
  ImgO = &ImageOld;
  imgdecode_init_hack(Img);
  imgdecode_init_hack(ImgO);

  //set environment variables for hk code
  //create filename and path
  if(grounddata) strcpy(envfile, ENVFILE_GND );
  else strcpy(envfile, ENVFILE );
  //fopen file
  if(!(fp=fopen(envfile, "r"))) {
      printk("***Can't open %s. Check setting is correct\n", envfile);
      exit(0);
  }
  //read in lines
  while( fgets(line, MAXLINE_IN_FILE, fp) != NULL )
  {
    if (!strncmp(line, "#", 1)) {
      continue; /*skip comments */
    }
    else  {
      sscanf(line, "%s %s %s ", s1,s2,s3);
      // set each line to setenv();
      setenv(s2, s3, 1);
    }
  }
  fclose(fp);

/***************************!!TBD ck this out later*************************
  //setup a sigalrm for threads
  //pthread_mutex_init(&mutex, NULL);

  // Set alarm
  gettimeofday(&tv0, NULL);
  if (td_createalarm(ALRMSEC, shandler, &mutex, &talarm)) { //to shandler in ALRMSEC
    pthread_mutex_destroy(&mutex); 
    printk("Can't set an alarm signal for %dsec of no data\n", ALRMSEC);
  }
*****************************************************************************/
}

// Module main function. 
int DoIt(void)
{
  pid_t pid;
  int wflg = 1;
  char *args[6];
  char callcmd[256];

  if (nice_intro())
    return (0);
  if(!(username = (char *)getenv("USER"))) username = "nouser"; 
  vc = cmdparams_get_str(&cmdparams, "vc", NULL);
  tlmdir = cmdparams_get_str(&cmdparams, "indir", NULL);
  outdir = cmdparams_get_str(&cmdparams, "outdir", NULL);
  logfile = cmdparams_get_str(&cmdparams, "logfile", NULL);
  if(strstr(tlmdir, "rexmit")) { rexmitmode = 1; }
  if (strcmp(vc, NOTSPECIFIED) == 0) {
    fprintf(stderr, "'vc' virt channel must be specified.  Abort\n");
    return(1);
  }
  if (strcmp(tlmdir, NOTSPECIFIED) == 0) {
    fprintf(stderr, "'indir' must be specified.  Abort\n");
    return(1);
  }
  if (strcmp(outdir, NOTSPECIFIED) == 0) {
    outdir = NULL;
  }
  if (strcmp(logfile, NOTSPECIFIED) == 0) {
    sprintf(logname, H0LOGFILE, username, vc, gettimetag());
  }
  else {
    sprintf(logname, "%s", logfile);
  }
  if(restartflg) {
    //sleep(30);		//make sure old ingest_lev0 is done
    if((h0logfp=fopen(logname, "a")) == NULL)
      fprintf(stderr, "**Can't open for append the log file %s\n", logname);
  }
  else {
    if((h0logfp=fopen(logname, "w")) == NULL)
      fprintf(stderr, "**Can't open the log file %s\n", logname);
  }
  setup();
  while(wflg) {
    do_ingest();                // loop to get files from the input dir 
/************************NOOP this out for now********************************
    if(sigalrmflg) {		// process an alarm timout for no data in
      sigalrmflg = 0;
      td_destroyalarm(&talarm, &mutex);
      gettimeofday(&tv0, NULL);
      if (td_createalarm(ALRMSEC, shandler, &mutex, &talarm)) { 
        pthread_mutex_destroy(&mutex);
        printk("Can't set an alarm signal for %dsec of no data\n", ALRMSEC);
      }
      if(!ignoresigalrmflg && firstfound) {
        if(Image.initialized) {
          if(rs && fsn_prev) {	//make sure have a created record
            close_image(rs, segment, segArray, &Image, fsn_prev);
            printk("*Closed image on timeout FSN=%u\n", fsn_prev);
          }
        }
        else {
          printk("*No image to close on timeout\n");
        }
        printk("alrm_sig: drms_server_end_transaction()\n");
        drms_server_end_transaction(drms_env, 0 , 0); //commit
        printk("alrm_sig: drms_server_begin_transaction()\n");
        drms_server_begin_transaction(drms_env); //start another cycle
        fsn_prev = 0;	//make sure don't try to flush this image
        imagecnt = 0;	//and start a new transaction interval
      }
      ignoresigalrmflg = 0;
    }
*****************************************************************************/

    if((stat(stopfile, &stbuf) == 0) || abortflg) {
      printk("Abort or Found file: %s. Terminate.\n", stopfile);
      if(abortflg) send_mail("Abort for ingest_lev0 for %s\n", pchan);
      //now close any open image
      if(Image.initialized) {
        if(rs) {		//make sure have a created record
          close_image(rs, segment, segArray, &Image, fsn_prev);
        }
      }
      else {
        printk("*No image to close on exit\n");
      }
      //No need to commit here. the exit will do it
      //printk("restart: drms_server_end_transaction()\n"); 
      //drms_server_end_transaction(drms_env, 0 , 0); //commit
      //NEW 1Jun2010: leave the stopfile. 
      //Must use doingestlev0_HMI(AIA).pl to start ingest_lev0
      //sprintf(callcmd, "/bin/rm -f %s", stopfile);
      //system(callcmd);
      if(rexmitmode)		//new 11/02/21012
        sprintf(callcmd, "touch /usr/local/logs/lev0/%s_exitX", pchan);
      else
        sprintf(callcmd, "touch /usr/local/logs/lev0/%s_exit", pchan);
      system(callcmd);		//let the world know we're gone
      wflg = 0; //leave DoIt()
      continue;
    }
    sleep(sleep_interval);	//normally 2 sec
    if(cntsleeps == 0) {	//file was seen
      if(paused) {		//send resume data flow msg
        paused = 0;
        if(!rexmitmode) {
          send_mail("tlm files seen again for ingest_lev0 for %s\n", pchan);
        }
      }
    }
    cntsleeps++;		//#of 2sec sleeps w/o any files in do_ingest()
    if(cntsleeps > 300) {	// >600sec w/o any files
      if(!paused) {
        if(!rexmitmode) {
          send_mail("No files seen for ingest_lev0 for %s for 600sec\n", pchan);
        }
        paused = 1;
      }
    }
  }
  /*************!!TBD noop this out for now
  td_destroyalarm(&talarm);
  //pthread_mutex_destroy(&mutex);  //can't destroy here. used in shandler()
  pthread_mutex_destroy(&mutex);  //if don't detroy here can end up w/2 ingest_lev0
				  //on a restart by doingestlev0.pl ??
  ***************************************************/
  return(0);
}

