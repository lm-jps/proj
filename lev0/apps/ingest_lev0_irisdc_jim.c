/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev0/apps/ingest_lev0_irisdc.c
 *-----------------------------------------------------------------------------
 *
 * This is a module that runs with DRMS and continuously extracts images and HK
 * data from .tlm files that appear in the given input dir. It puts the .tlm
 * and .qac files in the DRMS dataset TLMSERIESNAME and also in the outdir if
 * given.  It extracts images from the .tlm files an puts them in
 * the DRMS dataset LEV0SERIESNAME and extracts hk data to appropriate hk datasets.
 * (NOTE: Initially uses TLMSERIESNAMEHMI for convenience. But it's for IRIS)
 * The ISP dataset name is given in the file SOURCE_ENV_FOR_HK_DECODE_IRIS 
 * defined by ENVFILE.
 *
 * Call on datacapture front end by socdciris:
 * ingest_lev0_irisdc -l vc=VC03 indir=/sds/soc2soc/iris outdir=/sds/soc2pipe/iris
 * pipedir=/sds/pipe2soc/iris logfile=
 * /usr/local/logs/soc/soc_iris_VC03_prodtest_2013.07.15_16:49:17.log 
 * JSOC_DBNAME=irisdb JSOC_DBHOST=irisdc
 *
 * NOTE: If outdir is given, then we do not remove files from it.
 * NOTE: outdir can not be used for running in the backend (cl1n001) so that 
 * we can rm files.
 *
 * This module is also run on the JSOC backend (typically cl1n001) to get the tlm
 * and lev0 data into the JSOC DB. Here it is called by doingestlev0_IRIS.pl:
 * ingest_lev0_irisdc vc=VC03 indir=/sds/soc2pipe/iris [logfile=name]
 * ingest_lev0_irisdc vc=VC03 indir=/sds/soc2pipe/iris/rexmit [logfile=name]
 *
 * The /sds/soc2pipe/iris on the irisdc machine must be NFS'd to the backend 
 * machine.
 * NOTE: !!This module is only valid for IRIS. IRIS uses hmiaiaflg=0
 *
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
#include <signal.h>
#include "packets.h"
#include "imgdecode_iris.h"
#include "decode_hk_vcdu.h"
#include "decode_hk.h"
#include "load_hk_config_files.h"
#include "add_small_image.c"
#include "mypng.h"
#include "tdsignals.h"
#include "quallev0_iris.h"

#define RESTART_CNT 2	//#of tlm files to process before restart
#define MVDSFDIR "/sds/soc2pipe/iris/dsf/"

//use the hmi names for iris. And hmiaiaflg=0
#define TLMSERIESNAMEHMIAMES  "iris_ground.tlm_ames" 
#define TLMSERIESNAMEHMI  "iris.tlm_jim" 
//#define LEV0SERIESNAMEHMI "iris_ground.lev0_ROT" 
//#define LEV0SERIESNAMEHMI "iris_ground.lev0_dc1" 
#define LEV0SERIESNAMEHMI "iris.lev0_jim" 
#define LEV0SERIESNAMEHMIAMES "iris_ground.lev0_ames" 

#define WDLOGDIR "/dds/logs"
#define WDLOGFILE "wdlog"

#define H0LOGFILE "/usr/local/logs/lev0/ingest_lev0_iris.%s.%s.%s.log"
#define PKTSZ 1788		//size of VCDU pkt
//#define DEFAULTDB "irisdb_sums" //the default db to connect to
#define DEFAULTDB "jsoc_sums" //the default db to connect to
#define MAXFILES 65535		//max # of file can handle in tlmdir
#define NUMTIMERS 8		//number of seperate timers avail
//#define IMAGE_NUM_COMMIT 12	//number of complete images until commit
#define IMAGE_NUM_COMMIT 1	//!!TEMP number of complete images until commit
#define TESTAPPID 0x199		//appid of test pattern packet
#define TESTVALUE 0xc0b		//first value in test pattern packet
#define MAXERRMSGCNT 10		//max # of err msg before skip the tlm file
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define ENVFILE "/home/prodtest/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE_IRIS_jim"
#define ENVFILEA "/home/prodtest/cvs/IRIS/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE_AMES"
#define ENVFILE_GND "/home/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE_GROUND"

extern int decode_next_hk_vcdu(unsigned short *tbuf, CCSDS_Packet_t **hk, unsigned int *Fsn);
extern int write_hk_to_drms();
extern void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg);
//extern int set_HMI_mech_values(DRMS_Record_t *rec);
extern void sprint_time_ISOX (char *tstring, TIME t);
static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

int nxlook[5][9] = {
  {1096, 0, 0, 0, 0, 0, 0, 0, 0},          //0,0 
  {0, 1096, 1096, 0, 1096, 0, 0, 0, 1096}, //1,1 1,2 1,4 1,8
  {0, 548, 548, 0, 548, 0, 0, 0, 548},     //2,1 2,2 2,4 2,8
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 274, 274, 0, 274, 0, 0, 0, 274}      //4,1 4,2 4,4 4,8
};

int nylook[5][9] = {
  {4144, 0, 0, 0, 0, 0, 0, 0, 0},         //0,0 
  {0, 4144, 2072, 0, 1036, 0, 0, 0, 518}, //1,1 1,2 1,4 1,8
  {0, 4144, 2072, 0, 1036, 0, 0, 0, 518}, //2,1 2,2 2,4 2,8
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 4144, 2072, 0, 1036, 0, 0, 0, 518}   //4,1 4,2 4,4 4,8
};

int nxlookTap[5][9] = {
  {1096, 0, 0, 0, 0, 0, 0, 0, 0},          //0,0 
  {0, 1096, 1096, 0, 1096, 0, 0, 0, 1096}, //1,1 1,2 1,4 1,8
  {0, 548, 548, 0, 548, 0, 0, 0, 548},     //2,1 2,2 2,4 2,8
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 274, 274, 0, 274, 0, 0, 0, 274}      //4,1 4,2 4,4 4,8
};

int nylookTap[5][9] = {
  {2072, 0, 0, 0, 0, 0, 0, 0, 0},         //0,0 
  {0, 2072, 1036, 0, 518, 0, 0, 0, 259}, //1,1 1,2 1,4 1,8
  {0, 2072, 1036, 0, 518, 0, 0, 0, 259}, //2,1 2,2 2,4 2,8
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 2072, 1036, 0, 518, 0, 0, 0, 259}   //4,1 4,2 4,4 4,8
};

//Use defaults if lookup gives 0. (Also output warning message)
int nxlookdefault = 1096;
int nylookdefault = 4144;
int nxlookTapdefault = 1096;
int nylookTapdefault = 2072;

struct nextimage {
  DRMS_Record_t *rsN;
  DRMS_Segment_t *segN;
  IMG *imgN;
  int fsnN;
};

struct nextimage NextImage;
static int OpenNextImg = 0;

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "vc", NOTSPECIFIED, "Primary virt channel to listen to"},
  {ARG_STRING, "indir", NOTSPECIFIED, "directory containing the files to ingest"},
  {ARG_STRING, "outdir", NOTSPECIFIED, "directory to move the files to after the ingest"}, 
  {ARG_STRING, "pipedir", NOTSPECIFIED, "directory to get parc files from JSOC backend"}, 
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_FLAG, "A", "0", "use Ames datasets"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_FLAG, "l", "0", "wdlog flag"},
  {ARG_FLAG, "c", "0", "continue processing the current fsn after a imgdecode() error"},
  {ARG_END}
  //No -g flg for ingest iris {ARG_FLAG, "g", "0", "ground data flag"},
};

CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "ingest_lev0_iris";

typedef enum
{
   kThreadSigErr_Success
} ThreadSigErr_t;
static short sumsptrl, sumspat;
static int printflg = 0;
static int errskip=0;
static int seqerror = 0;
static int retardskip=0;
static int nx, ny, npix;
static int gLoop = 1;
static pthread_mutex_t mutex;
static struct timeval tv0;
static td_alarm_t talarm = 0;
FILE *wdlogfp = NULL;	//CARL-ingest_lev0_iris
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
static DRMS_SegmentDimInfo_t DimInfo;
static TIME sdo_epoch;
static char datestr[32];
static char bld_vers[16];
//char dtime[80];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];
static char Indir[300];  //CARL-ingest_lev0_iris
char wdlogname[128],oldwdlogname[128],Incmd[128];  //CARL-ingest_lev0_iris

int INVALtime;
unsigned int fsn = 0;
unsigned int fsnx = 0;
unsigned int fsnISP = 0;
unsigned int fsnISPTOCLOSE = 0;
unsigned int fsnISPX = 0;
unsigned int fsnISP_noop = 0;
unsigned int fsn_prev = 0;
unsigned int ispfound = 0;
unsigned int fsn_pre_rexmit = 0;
unsigned int fid = 0;
unsigned int tapcode = 0;
unsigned int cropid = 0;
unsigned int isysn = 0;

short *rdat;
short tmpdat[MAXPIXELS];
long long vcdu_seq_num;
long long vcdu_seq_num_next;
long long total_missing_im_pdu;
unsigned int vcdu_24_cnt, vcdu_24_cnt_next;
int verbose;
int grounddata = 0;
int appid;
int testid1, testid2;
int hmiaiaflg;			//0=hmi, 1=aia
int whk_status;
int total_tlm_vcdu;
int total_missing_vcdu;
int errmsgcnt, fileimgcnt;
int cntsleeps = 0;
int nofiletimeout = 0;
int timeoutclose = 0;           //set if 10 sec t.o. and a close_image was forced
int paused = 0;
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when ingest_lev0 is called for a restart
int amesflg = 0;		// set when use Ames ds names
int logflg = 0;			// set if keep wdlogs in /dds
int abortflg = 0;
int conterr = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
int firstfound = 0;		// set if see any file after startup
int ALRMSEC = 60;               // must get 2 in a row for no image timeout
int sleep_interval = 2;		// #of sec to sleep after do_ingest() calls
char ispquery[256];		// query to open isp ds
char logname[128];
char argvc[32], argindir[96], arglogfile[96], argoutdir[96], argpipedir[96];
char timetag[32];
char pchan[8];			// primary channel to listen to e.g. VC02 
char rchan[8];			// redundant channel to listen to e.g. VC10 
char stopfile[80];		// e.g. /usr/local/logs/lev0/VC04_stop
char tlmseriesname[128];	// e.g. hmi.tlm
char lev0seriesname[128];	// e.g. hmi.lev0
char tlmnamekey[128];		// shortened tlm file name for TLMDSNAM keyword
char tlmnamekeyfirst[128];	// for TLMDSNAM keyword for 1st file of image
char oldtlmdsnam[128];		// TLMDSNAM keyword from prev rec in db
char *dbname = DEFAULTDB;
char *username;			// from getenv("USER") 
char *tlmdir;			// tlm dir name passed in 
char *outdir;			// output dir for .tlm file or null
char *pipedir;			// get .parc files from JSOC backend
char *logfile;			// optional log name passed in 
char *vc;			// virtual channel to process, e.g. VC02 
DRMS_RecordSet_t *RSISP, *RSISPTO;
struct stat stbuf;
struct p_r_chans {
  char *pchan;
  char *rchan;
  int instru;		//0=hmi, 1=aia
};
typedef struct p_r_chans P_R_CHANS;

P_R_CHANS p_r_chan_pairs[] = {
//Only valid for IRIS
//{"VC01", "VC09", 1},		// AIA 
//{"VC04", "VC12", 1},		// AIA 
//{"VC02", "VC10", 0},		// HMI 
//{"VC05", "VC13", 0},		// HMI 
//NOTE: for IRIS there is no redundant channel. But use VC05 so we can test w/it
{"VC03", "VC05", 0},            // IRIS. Use hmiaiaflg=0
{"n/a", "n/a"}
};

struct namesort {		// sorted file names in tlmdir 
  char *name;
};
typedef struct namesort NAMESORT;


int nice_intro ()
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\ningest_lev0_irisdc_jim [-vh] "
	"vc=<virt chan> indir=</dir> [outdir=</dir>] [pipedir=</dir>] [logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
//        "  -g: output to hmi_ground/aia_ground\n"
        "  -l: keep logs under /dds (usually do for front end only)\n"
        "  -c: continue processing the current fsn after a imgdecode() error\n"
	"  -r: restart. only used when we restart ourselves periodically\n"
        "  -A: Ames. Use Ames ds e.g. iris_ground.tlm_ames\n"
	"vc= primary virt channel to listen to e.g. VC02\n"
	"indir= directory containing the files to ingest\n"
	"outdir= optional dir to copy the files to after the ingest\n"
	"        Do not specify this when run in the backend (cl1n001) so rm can be done\n"
	"pipedir= optional dir to get .parc files from JSOC backend\n"
	"logfile= optional log file name. Will create one if not given\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  //grounddata = cmdparams_get_int (&cmdparams, "g", NULL);
  conterr = cmdparams_get_int (&cmdparams, "c", NULL);
  restartflg = cmdparams_get_int (&cmdparams, "r", NULL);
  amesflg = cmdparams_get_int (&cmdparams, "A", NULL);
  logflg = cmdparams_get_int (&cmdparams, "l", NULL);
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
  sprintf(cmd, "echo \"%s\" | Mail -s \"ingest_lev0_iris mail\" lev0_user", string);
  system(cmd);
  va_end(args);
  return(0);
}

//Set the QUALLEV0 keyword
void do_quallev0(DRMS_Record_t *rs, IMG *img, int fsn) 
{
  char *iissloop, *cdark, *cled;
  char *hseqerr, *aistate;
  char wave_str[16];
  int status, isqfsn, asqfsn;
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
if(cled = drms_getkey_string(rs, "IMG_TYPE", &status)) {
  if(!strcmp(cled, "LED")) quallev0 = quallev0 | Q_LED;
  free(cled);
}
  isqfsn = drms_getkey_int(rs, "ISQFSN", &status);
  if(status || (fsn != isqfsn)) quallev0 = quallev0 | Q_NOISP;
  if(iissloop = drms_getkey_string(rs, "IISSLOOP", &status)) {
    if(!strcmp(iissloop, "OPEN")) quallev0 = quallev0 | Q_ISSOPEN;
    free(iissloop);
  }
  drms_setkey_int(rs, "QUALITY", quallev0);
  drms_setkey_string(rs, "BLD_VERS", bld_vers); //build vers to every record
  percentd = (float)((100.0 * (float)img->datavals)/(float)img->totalvals);
  drms_setkey_float(rs, "PERCENTD", percentd);
}

// Close out an image.
void close_image(DRMS_Record_t *rs, DRMS_Segment_t *seg, DRMS_Array_t *array,
		IMG *img, int fsn)
{
  DRMS_RecordSet_t *rsisp;
  DRMS_Record_t *rsispr;
  STAT stat;
  int status, n, k, ival, i, j, nnx, nny;
  int drms_status = 0;
  short sval;
  uint32_t missvals;
  long long cmdx, lval;
  TIME cmdxd;
  char *cptr, *iname;
  char tlmdsname[128];
  char queryext[256];
  char recdir[DRMS_MAXPATHLEN];//CARL ADDED for ingest_lev0_iris

  if(!img->initialized)  {
    printk("**No image to close in close_image() for fsn = %u\n", fsn);
    return;
  }
  printk("*Closing image for fsn = %u\n", fsn);
  if(seqerror) {
    printk("**Sequence error during fsn %u. No closed image\n", fsn);
    seqerror = 0;
    return;
  }
  if(ispfound != fsn) {
    printk("**No ISP found for fsn %u. No closed image\n", fsn);
    return;
  }
  if(imgstat_iris(img, &stat)) {
    printk("**Error on imgstat_iris() for fsn = %u\n", fsn);
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
  else {
    //For IRIS!!! set the camera # from the telescope # 
    //No, this is now done in hmi_time_setting.c
    //drms_setkey_int(rs, "CAMERA", (img->telnum)+1);
    //printf("telnum is %d\n", (img->telnum)+1);
  }
  drms_setkey_int(rs, "IMGAPID", img->apid);
  drms_setkey_int(rs, "CROPID", img->cropid);
  drms_setkey_int(rs, "ISYSN", img->isysn);
  if(drms_setkey_int(rs, "LUTID", img->luid)) {
    printk("ERROR on setkey_int for LUTID\n");
  }
  status = drms_setkey_short(rs, "SUMSPTRL", sumsptrl);
  status = drms_setkey_short(rs, "SUMSPAT", sumspat);
  drms_setkey_int(rs, "TAPCODE", img->tap);
  //drms_setkey_int(rs, "TAPCODE", tapcode);  //don't use from the img struct
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
  if(!timeoutclose) {
    rsisp = RSISP;      //saved rec set pointer for ISP
  }
  else {                //must reopen
    timeoutclose = 0;
    /* open ISP for drms series. ispquery saved from write_hk_to_drms_iris.c */
    sprintf(queryext, "%s[? I_SQ_FRAME_SN=%d ?]", ispquery, fsn);
    printk("queryext = %s\n", queryext);  //!!TEMP
    rsisp = drms_open_records(drms_env, queryext, &drms_status);
  }

    if(!rsisp || !rsisp->n || !rsisp->records || (fsn != fsnISP && !fsnISP_noop)) {
      printk("ERROR: Can't open isp record to put keywords in lev0.\n");
      printk("       The ISP was not received prior to the image for fsn %lu\n", fsn);
      printk("       Look for 'seq num out of sequence' error in log\n");
      printk("drms_status=%d, rsisp=%lu, fsn=%lu, fsnISP=%lu, fsnISP_noop=%lu\n",
                drms_status, rsisp, fsn, fsnISP, fsnISP_noop); //!!TEMP
      printk("       **ERROR: Can't close_image() for fsn = %u\n", fsn);
      img->initialized = 0;     //image is ready for use again. New 16Sep2013
      return;           //New. abort. 11Sep2013
    }
    else if(rsisp->n > 2 || rsisp->n < 0) {
      printk("ERROR: Got bad drms_open_records() pointer for isp. fsn=%lu\n", fsn);
      printk("       Look for 'seq num out of sequence' error in log\n");
      printk("       **ERROR: Can't close_image() fsn=%u\n", fsn);
      img->initialized = 0;     //image is ready for use again. New 16Sep2013
      return;           //New. abort. 11Sep2013
    }
    else {
      fsnISP_noop = 0;
      rsispr = rsisp->records[0];
      //NOTE: Can't use drms_copykeys() because source is long keywords
      //and destination is short keywords.
      //drms_status = drms_copykeys(rs, rsispr, 0, kDRMS_KeyClass_Explicit);
      //if(drms_status != DRMS_SUCCESS) {
      //  printk("ERROR: in drms_copykeys() to copy isp keys to lev0. Proceed\n");
      //}
      iname = getenv("HK_ISP_IRIS_DSNAME");
      status = drms_setkey_string(rs, "ISPSNAME", iname);
      iname = drms_getkey_string(rsispr, "PACKET_VERSION_NUMBER",&status);
      status = drms_setkey_string(rs, "ISPPKTVN", iname);
      iname = drms_getkey_string(rsispr, "HK_SOURCE",&status);
      status = drms_setkey_string(rs, "HKSOURCE", iname);
      cmdxd = drms_getkey_double(rsispr, "PACKET_TIME", &status);
      status = drms_setkey_double(rs, "ISPPKTIM", cmdxd);
      cmdx = drms_getkey_longlong(rsispr, "APID56_TIMECODE_SEC", &status);
      status = drms_setkey_longlong(rs, "ITCS56", cmdx);
      cmdx = drms_getkey_longlong(rsispr, "APID56_TIMECODE_SSEC", &status);
      status = drms_setkey_longlong(rs, "ITCSS56", cmdx);
      ival = drms_getkey_int(rsispr, "APID56_APID_VERSION", &status);
      status = drms_setkey_int(rs, "IVER56", ival);
      sval = drms_getkey_short(rsispr, "I_SQ_ISYS_NUM", &status);
      status = drms_setkey_short(rs, "ISQISYSN", sval);
      ival = drms_getkey_int(rsispr, "I_SQ_FRAME_SN", &status);
      status = drms_setkey_int(rs, "ISQFSN", ival);
//printf("In close_image() I_SQ_FRAME_SN = %d\n", ival); //!!TEMP
      //start NEW 03Apr2012
      sval = drms_getkey_short(rsispr, "I_IMG_CRS_ID", &status);
      status = drms_setkey_short(rs, "IICRSID", sval);
      //sval = drms_getkey_short(rsispr, "I_IMG_OBSLIST_ID", &status);
      //status = drms_setkey_short(rs, "IIOBSLID", sval);
      ival = drms_getkey_int(rsispr, "I_IMG_OBSLIST_ID", &status);
      status = drms_setkey_int(rs, "IIOBSLID", ival);
      //sval = drms_getkey_short(rsispr, "I_IMG_FRMLIST_ID", &status);
      //status = drms_setkey_short(rs, "IIFRMLID", sval);
      ival = drms_getkey_int(rsispr, "I_IMG_FRMLIST_ID", &status);
      status = drms_setkey_int(rs, "IIFRMLID", ival);
      //sval = drms_getkey_short(rsispr, "I_IMG_SG_FUV_NUV_ID_1", &status);
      //status = drms_setkey_short(rs, "IFDBID1", sval);
      //ival = drms_getkey_int(rsispr, "I_IMG_SJI_ID_2", &status);
      //status = drms_setkey_int(rs, "IFDBID2", ival);
      //ival = drms_getkey_int(rsispr, "I_IMG_STATUS", &status);
      //status = drms_setkey_int(rs, "IIHIS5", ival);
      //sval = drms_getkey_short(rsispr, "I_IMG_FILTER_TYPE", &status);
      //status = drms_setkey_short(rs, "IIMGFTYP", sval);
      //iname = drms_getkey_string(rsispr, "I_IMG_ISS_LOOP",&status);
      iname = drms_getkey_string(rsispr, "I_ISS_LOOP",&status);
      status = drms_setkey_string(rs, "IISSLOOP", iname);
      sval = drms_getkey_short(rsispr, "I_IMG_CFG_DELAY_1", &status);
      status = drms_setkey_short(rs, "IIMGCFD1", sval);
      sval = drms_getkey_short(rsispr, "I_IMG_CFG_DELAY_2", &status);
      status = drms_setkey_short(rs, "IIMGCFD2", sval);
      sval = drms_getkey_short(rsispr, "I_IMG_CFG_DELAY_3", &status);
      status = drms_setkey_short(rs, "IIMGCFD3", sval);
      sval = drms_getkey_short(rsispr, "I_IMG_CFG_DELAY_4", &status);
      status = drms_setkey_short(rs, "IIMGCFD4", sval);
      sval = drms_getkey_short(rsispr, "I_GTP_SUNVECTOR_X", &status);
      status = drms_setkey_short(rs, "IGTPSVX", sval);
      sval = drms_getkey_short(rsispr, "I_GTP_SUNVECTOR_Y", &status);
      status = drms_setkey_short(rs, "IGTPSVY", sval);
      sval = drms_getkey_short(rsispr, "I_ISS_PZTOFFA", &status);
      status = drms_setkey_short(rs, "IISSPZTA", sval);
      sval = drms_getkey_short(rsispr, "I_ISS_PZTOFFB", &status);
      status = drms_setkey_short(rs, "IISSPZTB", sval);
      sval = drms_getkey_short(rsispr, "I_ISS_PZTOFFC", &status);
      status = drms_setkey_short(rs, "IISSPZTC", sval);
      //end NEW 03Apr2012
      //sval = drms_getkey_short(rsispr, "I_IMG_FC_POSITION", &status);
      //status = drms_setkey_short(rs, "IIFCPOS", sval);
      //sval = drms_getkey_short(rsispr, "I_IMG_FW_ENCODER", &status);
      //status = drms_setkey_short(rs, "IIFWPOS", sval);
      ival = drms_getkey_int(rsispr, "I_IMG_SH_CMDED_EXPOSURE", &status);
      status = drms_setkey_int(rs, "IIMGSHCE", ival);
      ival = drms_getkey_int(rsispr, "I_OBC_TIME_SH1_FUV_SEC", &status);
      status = drms_setkey_int(rs, "IIMGOTS1", ival);
      ival = drms_getkey_int(rsispr, "I_OBC_TIME_SH1_FUV_SS", &status);
      status = drms_setkey_int(rs, "IMGOTSS1", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_A_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "IFUVACT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_A_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "IFUVAOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_B_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "IFUVBCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_B_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "IFUVBOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_C_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "IFUVCCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_C_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "IFUVCOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_D_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "IFUVDCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_D_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "IFUVDOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_E_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "IFUVECT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_E_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "IFUVEOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_F_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "IFUVFCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_F_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "IFUVFOT", ival);
      ival = drms_getkey_int(rsispr, "I_OBC_TIME_SH2_NUV_SEC", &status);
      status = drms_setkey_int(rs, "IIMGOTS2", ival);
      ival = drms_getkey_int(rsispr, "I_OBC_TIME_SH2_NUV_SS", &status);
      status = drms_setkey_int(rs, "IMGOTSS2", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_A_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "INUVACT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_A_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "INUVAOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_B_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "INUVBCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_B_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "INUVBOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_C_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "INUVCCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_C_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "INUVCOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_D_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "INUVDCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_D_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "INUVDOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_E_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "INUVECT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_E_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "INUVEOT", ival);
      ival = drms_getkey_int(rsispr, "I_OBC_TIME_SH3_SJI_SEC", &status);
      status = drms_setkey_int(rs, "IIMGOTS3", ival);
      ival = drms_getkey_int(rsispr, "I_OBC_TIME_SH3_SJI_SS", &status);
      status = drms_setkey_int(rs, "IMGOTSS3", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_A_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "ISJIACT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_A_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "ISJIAOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_B_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "ISJIBCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_B_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "ISJIBOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_C_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "ISJICCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_C_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "ISJICOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_D_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "ISJIDCT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_D_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "ISJIDOT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_E_CLOSE_TIME", &status);
      status = drms_setkey_int(rs, "ISJIECT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_E_OPEN_TIME", &status);
      status = drms_setkey_int(rs, "ISJIEOT", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_AEC_TABLE_ID", &status);
      status = drms_setkey_int(rs, "IAECTID", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_HIST_BIN_1", &status);
      status = drms_setkey_int(rs, "IIHIS1", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_HIST_BIN_2", &status);
      status = drms_setkey_int(rs, "IIHIS2", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_HIST_BIN_3", &status);
      status = drms_setkey_int(rs, "IIHIS3", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_HIST_BIN_4", &status);
      status = drms_setkey_int(rs, "IIHIS4", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_HIST_BIN_5", &status);
      status = drms_setkey_int(rs, "IIHIS5", ival);
      ival = drms_getkey_int(rsispr, "I_IMG_HIST_FSN", &status);
      status = drms_setkey_int(rs, "IIHISFSN", ival);
      //New 27Apr2012
      ival = drms_getkey_int(rsispr, "I_MC_FM_CMDED_POS", &status);
      status = drms_setkey_int(rs, "IFMCPOS", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FM_POS", &status);
      status = drms_setkey_int(rs, "IFMPOS", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FUV_ENCODER", &status);
      status = drms_setkey_int(rs, "IFUVENC", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FW_CMDED_TGT", &status);
      status = drms_setkey_int(rs, "IFWCTGT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FW_ENCODER", &status);
      status = drms_setkey_int(rs, "IFWENC", ival);
      ival = drms_getkey_int(rsispr, "I_MC_FW_POS", &status);
      status = drms_setkey_int(rs, "IFWPOS", ival);
      sval = drms_getkey_short(rsispr, "I_GTP_BIAS_OFFSET_X", &status);
      status = drms_setkey_short(rs, "IGTPOFFX", sval);
      sval = drms_getkey_short(rsispr, "I_GTP_BIAS_OFFSET_Y", &status);
      status = drms_setkey_short(rs, "IGTPOFFY", sval);
      sval = drms_getkey_short(rsispr, "I_IMG_FRAME_TYPE", &status);
      status = drms_setkey_short(rs, "IIFRMTYP", sval);
      iname = drms_getkey_string(rsispr, "I_LED1_HW",&status);
      status = drms_setkey_string(rs, "ILED1HWS", iname);
      iname = drms_getkey_string(rsispr, "I_LED1_SW",&status);
      status = drms_setkey_string(rs, "ILED1SWS", iname);
      iname = drms_getkey_string(rsispr, "I_LED2_HW",&status);
      status = drms_setkey_string(rs, "ILED2HWS", iname);
      iname = drms_getkey_string(rsispr, "I_LED2_SW",&status);
      status = drms_setkey_string(rs, "ILED2SWS", iname);
      iname = drms_getkey_string(rsispr, "I_LED3_HW",&status);
      status = drms_setkey_string(rs, "ILED3HWS", iname);
      iname = drms_getkey_string(rsispr, "I_LED3_SW",&status);
      status = drms_setkey_string(rs, "ILED3SWS", iname);
      ival = drms_getkey_int(rsispr, "I_MC_NUV_ENCODER", &status);
      status = drms_setkey_int(rs, "INUVENC", ival);
      sval = drms_getkey_short(rsispr, "I_RT_PZTOFFA", &status);
      status = drms_setkey_short(rs, "IRTPZTA", sval);
      sval = drms_getkey_short(rsispr, "I_RT_PZTOFFB", &status);
      status = drms_setkey_short(rs, "IRTPZTB", sval);
      sval = drms_getkey_short(rsispr, "I_RT_PZTOFFC", &status);
      status = drms_setkey_short(rs, "IRTPZTC", sval);
      ival = drms_getkey_int(rsispr, "I_MC_SJI_ENCODER", &status);
      status = drms_setkey_int(rs, "ISJIENC", ival);
      sval = drms_getkey_short(rsispr, "I_HR_BIT_ID", &status);
      status = drms_setkey_short(rs, "ISQBITID", sval);
      sval = drms_getkey_short(rsispr, "I_HR_COMP_ID", &status);
      status = drms_setkey_short(rs, "ISQCMPID", sval);
      sval = drms_getkey_short(rsispr, "I_HR_TAP_CODE", &status);
      status = drms_setkey_short(rs, "IISQTAP", sval);
      sval = drms_getkey_short(rsispr, "I_WB_PZTOFFA", &status);
      status = drms_setkey_short(rs, "IWBPZTA", sval);
      sval = drms_getkey_short(rsispr, "I_WB_PZTOFFB", &status);
      status = drms_setkey_short(rs, "IWBPZTB", sval);
      sval = drms_getkey_short(rsispr, "I_WB_PZTOFFC", &status);
      status = drms_setkey_short(rs, "IWBPZTC", sval);
      ival = drms_getkey_int(rsispr, "I_MC_WM1_CMDED_DELAY", &status);
      status = drms_setkey_int(rs, "IWM1CDEL", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM1_POS", &status);
      status = drms_setkey_int(rs, "IWM1CPOS", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM1_CMDED_TGT", &status);
      status = drms_setkey_int(rs, "IWM1CTGT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM1_ENCODER", &status);
      status = drms_setkey_int(rs, "IWM1ENC", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM2_CMDED_DELAY", &status);
      status = drms_setkey_int(rs, "IWM2CDEL", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM2_POS", &status);
      status = drms_setkey_int(rs, "IWM2CPOS", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM2_CMDED_TGT", &status);
      status = drms_setkey_int(rs, "IWM2CTGT", ival);
      ival = drms_getkey_int(rsispr, "I_MC_WM2_ENCODER", &status);
      status = drms_setkey_int(rs, "IWM2ENC", ival);
      ival = drms_getkey_int(rsispr, "I_TAT_ACS_STATUS", &status);
      status = drms_setkey_int(rs, "ITACSSTS", ival);
      //added 15Aug2012
      ival = drms_getkey_int(rsispr, "I_ISS_CONTROL", &status);
      status = drms_setkey_int(rs, "IISSCTRL", ival);
      sval = drms_getkey_short(rsispr, "I_ISS_DIODES", &status);
      status = drms_setkey_short(rs, "IISSDIOD", sval);
      sval = drms_getkey_short(rsispr, "I_ISS_OFFSET_RANGE", &status);
      status = drms_setkey_short(rs, "IISSOFFR", sval);
      sval = drms_getkey_short(rsispr, "I_ISS_SOURCE", &status);
      status = drms_setkey_short(rs, "IISS_SRC", sval);
      ival = drms_getkey_int(rsispr, "I_ISS_X_ERROR", &status);
      status = drms_setkey_int(rs, "ISXER", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_X_ERROR_AVG", &status);
      status = drms_setkey_int(rs, "ISXERAVG", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_X_ERROR_MAX", &status);
      status = drms_setkey_int(rs, "ISXERMAX", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_X_ERROR_MIN", &status);
      status = drms_setkey_int(rs, "ISXERMIN", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_Y_ERROR", &status);
      status = drms_setkey_int(rs, "ISYER,", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_Y_ERROR_AVG", &status);
      status = drms_setkey_int(rs, "ISYERAVG", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_Y_ERROR_MAX", &status);
      status = drms_setkey_int(rs, "ISYERMAX", ival);
      ival = drms_getkey_int(rsispr, "I_ISS_Y_ERROR_MIN", &status);
      status = drms_setkey_int(rs, "ISYERMIN", ival);
      sval = drms_getkey_short(rsispr, "I_MC_FUV_CMDED_DELAY", &status);
      status = drms_setkey_short(rs, "IFUVCDEL", sval);
      sval = drms_getkey_short(rsispr, "I_MC_FUV_CMDED_EXP", &status);
      status = drms_setkey_short(rs, "IFUVCEXP", sval);
      sval = drms_getkey_short(rsispr, "I_MC_NUV_CMDED_DELAY", &status);
      status = drms_setkey_short(rs, "INUVCDEL", sval);
      sval = drms_getkey_short(rsispr, "I_MC_NUV_CMDED_EXP", &status);
      status = drms_setkey_short(rs, "INUVCEXP", sval);
      sval = drms_getkey_short(rsispr, "I_MC_SJI_CMDED_DELAY", &status);
      status = drms_setkey_short(rs, "ISJICDEL", sval);
      sval = drms_getkey_short(rsispr, "I_MC_SJI_CMDED_EXP", &status);
      status = drms_setkey_short(rs, "ISJICEXP", sval);
      ival = drms_getkey_int(rsispr, "I_SQ_FLT_N_REPEAT", &status);
      status = drms_setkey_int(rs, "IIFLNRPT", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_FLT_REPEAT", &status);
      status = drms_setkey_int(rs, "IIFLRPT", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_OLT_N_REPEAT", &status);
      status = drms_setkey_int(rs, "IIOLNRPT", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_OLT_REPEAT", &status);
      status = drms_setkey_int(rs, "IIOLRPT", ival);
      //ival = drms_getkey_int(rsispr, "I_SQ_FUV_FDBT_ID", &status);
      //status = drms_setkey_int(rs, "IIFUVFDB", ival);
      lval = drms_getkey_longlong(rsispr, "I_SQ_FUV_FDBT_ID", &status);
      status = drms_setkey_longlong(rs, "IIFUVFDB", lval);
      //ival = drms_getkey_int(rsispr, "I_SQ_NUV_FDBT_ID", &status);
      //status = drms_setkey_int(rs, "IINUVFDB", ival);
      lval = drms_getkey_longlong(rsispr, "I_SQ_NUV_FDBT_ID", &status);
      status = drms_setkey_longlong(rs, "IINUVFDB", lval);
      //ival = drms_getkey_int(rsispr, "I_SQ_SJI_FDBT_ID", &status);
      //status = drms_setkey_int(rs, "IISJIFDB", ival);
      lval = drms_getkey_longlong(rsispr, "I_SQ_SJI_FDBT_ID", &status);
      status = drms_setkey_longlong(rs, "IISJIFDB", lval);
      ival = drms_getkey_int(rsispr, "I_SQ_FDBT_IDX", &status);
      status = drms_setkey_int(rs, "ISQFDBTI", ival);
      sval = drms_getkey_short(rsispr, "I_SQ_OWT_ID", &status);
      status = drms_setkey_short(rs, "ISQOWTID", sval);
      sval = drms_getkey_short(rsispr, "I_SQ_SRT_ID", &status);
      status = drms_setkey_short(rs, "ISQSRTID", sval);
      //status = drms_setkey_short(rs, "SUMSPTRL", sumsptrl);
      //status = drms_setkey_short(rs, "SUMSPAT", sumspat);
      //Oct 31, 2012
      iname = drms_getkey_string(rsispr, "I_AEC_FLAG",&status);
      status = drms_setkey_string(rs, "IAECFLAG", iname);
      iname = drms_getkey_string(rsispr, "I_AEC_EV_FLAG",&status);
      status = drms_setkey_string(rs, "IAECEVFL", iname);
      iname = drms_getkey_string(rsispr, "I_AEC_FLARE_FLAG",&status);
      status = drms_setkey_string(rs, "IAECFLFL", iname);
      ival = drms_getkey_int(rsispr, "I_AEC_ENA_MSK", &status);
      status = drms_setkey_int(rs, "IAECENAM", ival);
      ival = drms_getkey_int(rsispr, "I_AEC_LAPSE_TIME", &status);
      status = drms_setkey_int(rs, "IAECLTIM", ival);
      ival = drms_getkey_int(rsispr, "I_AEC_EV_LAPSE_TIME", &status);
      status = drms_setkey_int(rs, "IAECEVLT", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_FLT_IDX", &status);
      status = drms_setkey_int(rs, "ISQFLTDX", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_FLT_N_IDX", &status);
      status = drms_setkey_int(rs, "ISQFLTNX", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_OLT_IDX", &status);
      status = drms_setkey_int(rs, "ISQOLTDX", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_OLT_N_IDX", &status);
      status = drms_setkey_int(rs, "ISQOLTNX", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_PZT_A", &status);
      status = drms_setkey_int(rs, "ISQPZTA", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_PZT_B", &status);
      status = drms_setkey_int(rs, "ISQPZTB", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_PZT_C", &status);
      status = drms_setkey_int(rs, "ISQPZTC", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_ORB_PERIOD", &status);
      status = drms_setkey_int(rs, "ISQORBPD", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_ORB_EQ_SEC", &status);
      status = drms_setkey_int(rs, "ISQORBSE", ival);
      ival = drms_getkey_int(rsispr, "I_SQ_ORB_EQ_SUB", &status);
      status = drms_setkey_int(rs, "ISQORBSU", ival);
      //NEW 20May2013
      lval = drms_getkey_longlong(rsispr, "I_SQ_OLT_ID", &status);
      status = drms_setkey_longlong(rs, "ISQOLTID,", lval);
      lval = drms_getkey_longlong(rsispr, "I_SQ_FLT_ID", &status);
      status = drms_setkey_longlong(rs, "ISQFLTID,", lval);
      //NEW 21May2013
      int isqisysn      = drms_getkey_int(rsispr, "I_SQ_ISYS_NUM", &status);
      long long iifuvfdb = drms_getkey_longlong(rsispr, "I_SQ_FUV_FDBT_ID",&status);
      long long iinuvfdb = drms_getkey_longlong(rsispr, "I_SQ_NUV_FDBT_ID",&status);
      long long iisjifdb = drms_getkey_longlong(rsispr, "I_SQ_SJI_FDBT_ID",&status);
      switch (isqisysn) {
         case 0:
            drms_setkey_longlong(rs, "IIFDBID", iifuvfdb);
            break;
         case 1:
            drms_setkey_longlong(rs, "IIFDBID", iinuvfdb);
            break;
         case 2:
            drms_setkey_longlong(rs, "IIFDBID", iisjifdb);
            break;
       }

      //drms_close_record(rsispr, DRMS_FREE_RECORD);
      drms_close_records(rsisp, DRMS_INSERT_RECORD);
      //printk("!!TEMP set rsisp = 0 in close_image()\n");
      //RSISP = 0;  
    }
  //Now rotate image counterclokwise 90 degrees
  nnx = img->nx;	//this is always nx > ny 
  nny = img->ny;
  npix = nnx * nny;
  for(i=0; i < nny; ++i) {
    for(j=0; j < nnx; ++j)
      tmpdat[i*nnx+j] = img->dat[j*nny+i];
  }
  memcpy(img->dat, tmpdat, 2*npix);
  status = drms_segment_writewithkeys(seg, array, 0);
  //status = drms_segment_write(seg, array, 0);
  if (status) {
    printk("ERROR: drms_segment_write error=%d for fsn=%u\n", status, fsn);
    printk("ABORT: close_image(). No drms_close_record() done\n");
    array->data = NULL;        // must do before free 
    if(array) drms_free_array(array);
    img->initialized = 0;		//indicate image is ready for use again
    return;				//Abort. New 11Sep2013
  }
  add_small_array(rs, array, 8, 16); //add Phil's png and small fits

  drms_record_directory(rs, recdir, 1);//CARL ADDED for ingest_lev0_iris
  //ingest_lev0_iris -CARL
  printf("recdir from ingest_lev0_iris:<%s>\n",recdir);//CARL ADDED PRINT
  if(wdlogfp)
  {
    fprintf(wdlogfp,"%s\n",recdir);
    fflush(wdlogfp);
    sprintf(Incmd,"/bin/rm -f %s/%09u.fits",Indir,fsn);/*in case there*/
    system(Incmd);
    //sprintf(Incmd,"ln -s %s  %s/",recdir,Indir);
    sprintf(Incmd,"ln -s %s/image.fits  %s/%09u.fits",recdir,Indir,fsn);
    h0log("%s\n",Incmd);/*!!TEMP*/
    system(Incmd);
  }
  

  array->data = NULL;        // must do before free 
  if(array) drms_free_array(array);
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
  printk("**Exit ingest_lev0_iris w/ status = %d\n", stat);
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
    NextImage.rsN = 0;	//when next data pkt come in, there is no record
    return(1);
  }
  dstatus = drms_setkey_int(rs, "FSN", fsnx);
  if(!(segment = drms_segment_lookup(rs, "image"))) {
    printk("No drms_segment_lookup(rs, image)\n");
    NextImage.rsN = 0;	//when next data pkt come in, there is no record
    return(1);
  }
  /* Changes to ingest_lev0 for IRIS include the vardim segment, which
   * cannot be opened until the next data packet comes in and the nx/ny code
   * can be gotten from the fid (i.e. low byte of CamHdr3). So info is saved
   * here for the new record to be setdim'd when the next packet comes in.
  */
  //if(drms_segment_setdims(segment, &DimInfo)) {
  //  printk("No drms_segment_setdims(segment, &DimInfo)\n");
  //  return(1);
  //}

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
    //imgdecode_init_hack(Img);
    imgdecode_iris_init_hack(Img);

  rdat = Img->dat;
  //segArray = drms_array_create(DRMS_TYPE_SHORT,
  //                                     segment->info->naxis,
  //                                     segment->axis,
  //                                     rdat,
  //                                     &dstatus);
  ImgC = Img;		//set current image pointer
  //set up info for when the next pkt comes in w/the new nx/ny bits.
  NextImage.rsN = rs;
  NextImage.segN = segment;
  NextImage.imgN = Img;
  NextImage.fsnN = fsnx;
  OpenNextImg = 1;
  return(0);
} 

// An fsn has changed in the context of a normal stream of
// tlm files. Returns 0 on success.
int fsn_change_normal()
{
  int dstatus, rstatus, compid, n, k;
  int i, j, nnx, nny;
  char reopen_dsname[256];
  char *cptr;

  printk("*FSN has changed from %u to %u\n", fsn_prev, fsnx);
  if(fsnx < fsn_prev) {
    retardskip = 1;
    printk("**FSN retardation. Skip %u until normal increase\n", fsnx);
  }
  else retardskip = 0;
  printflg = 0; //!!TEMP
  if(fsn_prev == 0) {	//startup mode. restore any prev image
    errmsgcnt = 0;
    errskip = 0;
    fsn_prev = fsnx;
    //if(hmiaiaflg)
      //sprintf(reopen_dsname, "%s[%u]", LEV0SERIESNAMEAIA, fsnx);
    //else
      sprintf(reopen_dsname, "%s[%u]", LEV0SERIESNAMEHMI, fsnx);
    printk("Open normal prev ds: %s\n", reopen_dsname);
    rset = drms_open_records(drms_env, reopen_dsname, &rstatus);
    if(rstatus) {
      printk("Can't do drms_open_records(%s)\n", reopen_dsname);
      return(1);          // !!!TBD
    }
    if(!rset || (rset->n == 0) || rstatus) {
startnew:
      printk("No prev ds\n");     // start a new image
      fsn_normal_new_image();
    }
    else {
      ispfound = fsnx;		//assume isp already in ds found
      Img->initialized = 1;
      Img->reopened = 1;
      Img->fsn = fsnx;
      Img->apid = appid;
      rs_old = rset->records[0];
      rs = drms_clone_record(rs_old, DRMS_PERMANENT, DRMS_COPY_SEGMENTS, &rstatus);
      if(rstatus || !rs) {
        printk("Can't do drms_clone_record()\n");
        goto startnew;		//new 30Sep2013
        //return(1);		// !!!TBD ck 
      }
      drms_close_records(rset, DRMS_FREE_RECORD);
      rstatus = drms_setkey_int(rs, "FSN", fsnx);
      Img->telnum = drms_getkey_int(rs, "CAMERA", &rstatus);
      Img->telnum--;
      Img->cropid = drms_getkey_int(rs, "CROPID", &rstatus);
      Img->isysn = drms_getkey_int(rs, "ISYSN", &rstatus);
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
      if(!segment) {
        printk("Can't do drms_segment_lookupnum()\n");
        return(1);
      }
      cArray = drms_segment_read(segment, DRMS_TYPE_SHORT, &rstatus);
      if(rstatus || !cArray) {
        printk("Can't do drms_segment_read()\n");
        return(1);		// !!!!TBD ck 
      }
      short *adata = (short *)cArray->data;
      //memcpy(Img->dat, adata, 2*MAXPIXELS);
      memcpy(Img->dat, adata, 2*segment->axis[0]*segment->axis[1]);
      rdat = Img->dat;
      //rotate image 90 deg clockwise
      nnx = segment->axis[0];
      nny = segment->axis[1];
      npix = nnx * nny;
      for(i=0; i < nnx; ++i) {
        for(j=0; j < nny; ++j)
          tmpdat[i*nny+j] = Img->dat[j*nnx+i];
      }
      memcpy(Img->dat, tmpdat, 2*npix);

      segArray = drms_array_create(DRMS_TYPE_SHORT, 
                                         segment->info->naxis,
                                         segment->axis,
                                         rdat,
                                         &dstatus);
      if(dstatus || !segArray) {
        printk("Can't do drms_array_create()\n");
        return(1);		// !!!!TBD ck 
      }
      ImgC = Img;           //set current image pointer
      drms_free_array(cArray); //must free from drms_segment_read()
    }
  }
  else {
    if(errskip) {     //prev fsn had error. don't try to close
        errskip = 0;
        printk("*Skip Closing image for fsn = %u\n", fsn_prev);
    }
    else {
      if(rs && segArray) {		// make sure have a created record
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
    }
    if(fsn_normal_new_image()) {
      printk("**ERROR: Can't create new image\n");
      return(1);
    }
  }
  //printk("*FSN %u tapcode=%d cropid=%d\n", fsnx, tapcode, cropid);
  return(0);
}

// An fsn has changed in the context of a retransmitted or higher version
// tlm file. Returns 0 on success.
int fsn_change_rexmit() 
{
  char rexmit_dsname[256];
  int rstatus, dstatus, compid, n, k;
  int i, j, nnx, nny;
  char *cptr;

  if(fsn_prev != 0) {   // close image of prev fsn if not 0 
    if(rsc) {		//make sure have created record
      if(errskip) {     //prev fsn had error. don't try to close
          errskip = 0;
          printk("*Skip Closing image for fsn = %u\n", fsn_prev);
      }
      else {
        close_image(rsc, segmentc, oldArray, ImgO, fsn_prev);
        imagecnt++;
        fileimgcnt++;
      }
    }
    else {
      printk("**ERROR: Null record ptr for an rexmit image fsn=%u\n",fsn_prev);
    }
  }
  errmsgcnt = 0;
  fsn_prev = fsnx;
  //if(hmiaiaflg) 
    //sprintf(rexmit_dsname, "%s[%u]", LEV0SERIESNAMEAIA, fsnx);
  //else 
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
    //imgdecode_init_hack(ImgO);
    imgdecode_iris_init_hack(ImgO);
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
    rstatus = drms_segment_setdims(segmentc, &DimInfo);
    rdat = ImgO->dat;
    oldArray = drms_array_create(DRMS_TYPE_SHORT,
                                       segmentc->info->naxis,
                                       segmentc->axis,
                                       rdat,
                                       &dstatus);
    if(dstatus || !oldArray) {
      printk("Can't do drms_array_create(%s)\n", rexmit_dsname);
      return(1);		// !!!TBD 
    }
  }
  else {
    ispfound = fsnx;		//assume isp already in ds found
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
    ImgO->isysn = drms_getkey_int(rsc, "ISYSN", &rstatus);
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
    if(!segmentc) {
      printk("Can't do drms_segment_lookupnum()\n");
      return(1);		// !!!!TBD ck 
    }
    cArray = drms_segment_read(segmentc, DRMS_TYPE_SHORT, &rstatus);
    if(rstatus) {
      printk("Can't do drms_segment_read()\n");
      return(1);		// !!!!TBD ck 
    }
    short *adata = (short *)cArray->data;
    //memcpy(ImgO->dat, adata, 2*MAXPIXELS);  
    memcpy(ImgO->dat, adata, 2*segmentc->axis[0]*segmentc->axis[1]); 
    rdat = ImgO->dat;
    //rotate image 90 deg clockwise
    nnx = segmentc->axis[0];
    nny = segmentc->axis[1];
    npix = nnx * nny;
    for(i=0; i < nnx; ++i) {
      for(j=0; j < nny; ++j)
        tmpdat[i*nny+j] = Img->dat[j*nnx+i];
    }
    memcpy(ImgO->dat, tmpdat, 2*npix);

    oldArray = drms_array_create(DRMS_TYPE_SHORT, 
                                       segmentc->info->naxis,
                                       segmentc->axis,
                                       rdat,
                                       &dstatus);
    if(dstatus || !oldArray) {
      printk("Can't do drms_array_create()\n");
      return(1);		// !!!!TBD ck 
    }
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
  int datval, eflg, firstflg, dstatus, ntmp;
  unsigned int cnt1, cnt2, cnt3, gap_24_cnt, nxbits, nybits;
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
      if(vcdu_24_cnt_next != 0) seqerror = 1;
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
      if(vcdu_seq_num_next != 0) seqerror = 1;
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
    //if(appid == TESTAPPID) "{" // appid of test pattern 
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
    //if(appid == APID_HMI_SCIENCE_1 || appid == APID_HMI_SCIENCE_2 || 
    //	appid == APID_AIA_SCIENCE_1 || appid == APID_AIA_SCIENCE_2)
    if(appid == APID_IRIS_SCIENCE)
    {
      cnt1 = MDI_getshort(cbuf+32);	//fsn high 2 are camera system number
      cnt2 = MDI_getshort(cbuf+34);
      isysn = (unsigned int)(cnt1<<16)+(unsigned int)(cnt2);
      isysn = isysn & 0xc0000000;	//high 2 bits
      isysn = (unsigned int)(isysn >> 30);

      cnt1 = MDI_getshort(cbuf+40);	//cropid 12 bits
      cropid = (unsigned int)(cnt1 >> 4);
      cnt1 = MDI_getshort(cbuf+42);	//tap code in high 4 bits
      tapcode = (unsigned int)(cnt1 >> 12);
      //nx,ny bits are in CamHdr3 low 8 bits
      cnt1 = MDI_getshort(cbuf+36);
      nxbits = (unsigned int)(cnt1 >> 4) & 0x0f;
      nybits = (unsigned int)cnt1 & 0x0f;
      sumsptrl = (short)nybits;	//for keyword SUMSPTRL (reversed)
      sumspat = (short)nxbits;		//for keyword SUMSPAT (reversed)
      //printk("cnt1 = %0x  nxbits/nybits = %d/%d\n", cnt1, nxbits, nybits);
      //nxbits = 0;		//force them to 0 for now !!TEMP
      //nybits = 0;		//force them to 0 for now !!TEMP
      if((tapcode == 10) || (tapcode == 9) || (tapcode == 0 && cropid>0 && isysn !=
0)) {
        nx = nxlookTap[nxbits][nybits];
        ny = nylookTap[nxbits][nybits];
        if(nx == 0) {
          nx = nxlookTapdefault;
          ny = nylookTapdefault;
          printk("WARNING: nxbits/nybits=%d/%d gave invalid nx of 0 for fsn=%u. Use default\n", 
		nxbits, nybits, fsnx);
        }
        if(ny == 0) {
          nx = nxlookTapdefault;
          ny = nylookTapdefault;
          printk("WARNING: nxbits/nybits=%d/%d gave invalid ny of 0 for fsn=%u. Use default\n", 
		nxbits, nybits, fsnx);
        }
      }
      else {
        nx = nxlook[nxbits][nybits];
        ny = nylook[nxbits][nybits];
        if(nx == 0) {
          nx = nxlookdefault;
          printk("WARNING: nxbits/nybits=%d/%d gave invalid nx of 0 for fsn=%u. Use default\n",
		 nxbits, nybits, fsnx);
        }
        if(ny == 0) {
          ny = nylookdefault;
          printk("WARNING: nxbits/nybits=%d/%d gave invalid ny of 0 for fsn=%u. Use default\n",
		 nxbits, nybits, fsnx);
        }
      }
      //This is the ROT version, so swap the nx and ny.
      ntmp = nx;
      nx = ny;
      ny = ntmp;
if(!printflg) {		//!!TEMP
  printflg = 1;
  printk("nx/ny = %d/%d for fsn = %u\n", nx, ny, fsnx);
  printk("*FSN %u tapcode=%d cropid=%d\n", fsnx, tapcode, cropid);
}
      DimInfo.naxis = 2;
      DimInfo.axis[0] = nx;
      DimInfo.axis[1] = ny;
      cnt1 = MDI_getshort(cbuf+32);
      cnt2 = MDI_getshort(cbuf+34);
      fsnx = (unsigned int)(cnt1<<16)+(unsigned int)(cnt2);
      fsnx = fsnx & 0x3fffffff;		//low 30bits for fsn */
      if(fsnx == 0) continue;		//a 0 fsn is not acceptable
      if(tapcode == 10) {
        //fsnx = fsnx + 100000000;
        fsnx++;		//treat as next fsn
        fsnISP = fsnx;  //and use last isp data
        fsnISP_noop = 1;
      }
      if(OpenNextImg) {
        rs = NextImage.rsN;
        segment = NextImage.segN;
        ImgC = NextImage.imgN;
        if(drms_segment_setdims(segment, &DimInfo)) {
          printk("No drms_segment_setdims(segment, &DimInfo)\n");
          return(1);
        }
        segArray = drms_array_create(DRMS_TYPE_SHORT,
                                     segment->info->naxis,
                                     segment->axis,
                                     rdat,
                                     &dstatus);
        if(dstatus || !segArray) {
          printk("Can't do drms_array_create()\n");
          return(1);		// !!!!TBD ck 
        }
        OpenNextImg = 0;
      }
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
          if(fsnx < fsn_prev) {		//the prev isp is not for this image
            printk("**WARNING: ISP out of order with image data\n");
            RSISP = 0;			//indicate to close_image() no isp
          }
          if(fsn_change_normal()) {	//handle old & new images
            printk("***FATAL ERROR in fsn_change_normal()\n");
            return(1);
          }
        }
      }
      // send the sci data to Keh-Cheng. call with pointer to M_PDU_Header 
      if(errskip || retardskip) { continue; }	//skip until new fsn
      ImgC->nx = nx;
      ImgC->ny = ny;
      rstatus = imgdecode_iris((unsigned short *)(cbuf+10), ImgC);
      switch(rstatus) {
      case 0:
        // A science data VCDU was successfully decoded 
        break;
      case IMGDECODE_DECOMPRESS_ERROR:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_DECOMPRESS_ERROR\n");
        break;
      case IMGDECODE_TOO_MANY_PIXELS:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_TOO_MANY_PIXELS\n");
        break;
      case IMGDECODE_BAD_N:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_N\n");
        break;
      case IMGDECODE_BAD_APID:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_APID\n");
        break;
      case IMGDECODE_NO_LOOKUP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_NO_LOOKUP_TABLE\n");
        break;
      case IMGDECODE_LOOKUP_ID_MISMATCH:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_LOOKUP_ID_MISMATCH\n");
        break;
      case IMGDECODE_BAD_LOOKUP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_LOOKUP_TABLE\n");
        break;
      case IMGDECODE_NO_CROP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_NO_CROP_TABLE\n");
        break;
      case IMGDECODE_CROP_ID_MISMATCH:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_CROP_ID_MISMATCH\n");
        break;
      case IMGDECODE_BAD_CROP_GEOMETRY:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_CROP_GEOMETRY\n");
        break;
      case IMGDECODE_BAD_CROP_TABLE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_CROP_TABLE\n");
        break;
      case IMGDECODE_BAD_CROP_SKIP_TAKE:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_CROP_SKIP_TAKE\n");
        break;
      case IMGDECODE_BAD_OFFSET:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_BAD_OFFSET\n");
        break;
      case IMGDECODE_OUT_OF_MEMORY:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: IMGDECODE_OUT_OF_MEMORY\n");
        break;
      default:
        if(errmsgcnt++ < MAXERRMSGCNT)
          printk("*imgdecode_iris() ret: unknown err status = %d:\n", rstatus);
        break;
      }
      if(rstatus) {
         if(!conterr) {
           errskip = 1;	//skip vcdu until new fsn
           printk("*FSN imgdecode() error. Skip fsn=%u\n", fsnx);
         }
         else printk("*FSN imgdecode() error. Continue processing fsn=%u\n", fsnx);
      }
    }
    else {			// send the HK data to Carl 
      //printk("$$$ appid assumed hk =  0x%x %d\n", appid, appid); //!!TEMP
      decode_status = decode_next_hk_vcdu((unsigned short *)(cbuf+10), &Hk, &Fsn);
      //printk("decode_status = %d\n", decode_status); //!!TEMP
      switch (decode_status) {
        case SUCCESS_HK_NEED_TO_WTD_CTD:
          printk("*ISP found for fsn = %u\n", Fsn);
          ispfound = fsnISPX;
          fsnISP = fsnISPX; 
          fsnISPX = Fsn;
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
printk("Call write_hk_to_drms() after SUCCESS_HK_NEED_TO_WTD_CTD for rexmit\n"); //!!TEMP
            whk_status = write_hk_to_drms(rsc, &Hk); //ISP keywords to drms
            //if(!hmiaiaflg) {
            //  if(whk_status = set_HMI_mech_values(rsc)) {
            //    printk("***ERROR: mechanism position keywords are wrong!\n");
            //  }
            //}
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
printk("Calling write_hk_to_drms() after SUCCESS_HK_NEED_TO_WTD_CTD for normal fsn change\n"); //!!TEMP
            whk_status = write_hk_to_drms(rs, &Hk); //ISP keywords to drms
            //if(!hmiaiaflg) {
            //  if(whk_status = set_HMI_mech_values(rs)) {
            //    printk("***ERROR: mechanism position keywords are wrong!\n");
            //  }
            //}
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
  printk("%s\n**Processed %s\n**with %d images and %d VCDUs in %f sec\n\n",
	do_datestr(),file, fileimgcnt, fpkt_cnt, ftmp);

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
    //if(strstr(dp->d_name, pchan) || strstr(dp->d_name, rchan)) {
    if(strstr(dp->d_name, pchan) || strstr(dp->d_name, rchan) || strstr(dp->d_name, ".dsf")) {
      nameptr[i++].name = strdup(dp->d_name);
      if(i >= MAXFILES) {
        printk("***Fatal error. Too many (%d) files in %s\n", MAXFILES, tlmdir);
        printf("***Fatal error. Too many (%d) files in %s\n", MAXFILES, tlmdir);
        abortit(3);
      }
      cntsleeps = 0;		//we saw a file
      if(timeoutclose) {
        //now reset to 0 in close_image()
        //timeoutclose = 0;
        fsnISP = fsnISPTOCLOSE;
      }
    }
  }
  closedir(dfd);
  if(i == 0) {
    free(nameptr);
    return;
  }
  qsort(nameptr, i, sizeof(NAMESORT), &compare_names);
  nofiletimeout = atoi(getenv("FLUSH_2SEC_COUNT"));

  for(j=0; j < i; j++) {
    //printk("####QSORT FILES: %s\n", nameptr[j].name); // !!TEMP 
    // NOTE: the dsf files is moved to outdir/dsf e.g. (/dds/soc2pipe/iris/dsf)
//!!NOTE: for IRIS always mv the dsf file to /dds/soc2pipe/iris/dsf
//    if(outdir) {
        if(strstr(nameptr[j].name, ".dsf")) {
          //sprintf(cmd, "/bin/mv -f %s/%s %s/dsf/", tlmdir, nameptr[j].name, outdir);
          //printk("*mv dsf file to %s/dsf/\n", outdir);
          sprintf(cmd, "/bin/mv -f %s/%s %s", tlmdir, nameptr[j].name, MVDSFDIR);
          printk("*mv dsf file to %s\n", MVDSFDIR);
          printk("%s\n", cmd);
          if(system(cmd)) {
            printk("***Error on: %s\n", cmd);
          }
        }
//    }
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
    //copy qac and tlm into /SUMs
    sprintf(cmd, "/bin/cp -p %s %s", name, path);
    printk("*cp qac to %s\n", path);
    printk("%s\n", cmd);
    if(status = system(cmd)) {
      printk("**ERROR: %d errno=%d on: %s\n", status, errno, cmd);
      if(WIFEXITED(status)) {
        if(mvstat = WEXITSTATUS(status)) {  //status ret by mv
          printk("**ERROR: cp exit status = %d\n", mvstat);
        }
      }
    }
    sprintf(cmd, "/bin/cp -p %s %s", tlmfile, path);
    printk("*cp tlm to %s\n", path);
    printk("%s\n", cmd);
    if(status = system(cmd)) {
      printk("**ERROR: %d errno=%d on: %s\n", status, errno, cmd);
      if(WIFEXITED(status)) {
        if(mvstat = WEXITSTATUS(status)) {  //status ret by mv
          printk("**ERROR: cp exit status = %d\n", mvstat);
        }
      }
      sprintf(cmd, "/bin/mv -f %s %s", name, "/sds/reject/iris");
      printk("%s\n", cmd);
      if(status = system(cmd)) {
        printk("**ERROR: on mv\n");
      }
      printk("**Continue after ERROR on cp of tlm file\n");
      continue;
    }
    if(outdir) {		//mv tlm and qac to output dir
      sprintf(cmd, "/bin/mv -f %s %s", name, outdir);
      printk("*mv qac file to %s\n", outdir);
      printk("%s\n", cmd);
      if(system(cmd)) {
        printk("***Error on: %s\n", cmd);
      }
      sprintf(cmd, "/bin/mv -f %s %s", tlmfile, outdir);
      printk("*mv tlm file to %s\n", outdir);
      printk("%s\n", cmd);
      if(system(cmd)) {
        printk("***Error on: %s\n", cmd);
      }
    }
    else {			//if no outdir then this is backend so rm files
      sprintf(cmd, "/bin/rm -f %s", name);
      printk("%s\n", cmd);
      if(status = system(cmd)) {
        printk("**ERROR: %d errno=%d on: %s\n", status, errno, cmd);
        send_mail("ERROR: %s. We must be able to rm these files\n", cmd);
        if(WIFEXITED(status)) {
          if(mvstat = WEXITSTATUS(status)) {  //status ret by rm
            printk("**ERROR: rm exit status = %d\n", mvstat);
          }
        }
      }
      sprintf(cmd, "/bin/rm -f %s", tlmfile);
      printk("%s\n", cmd);
      if(status = system(cmd)) {
        printk("**ERROR: %d errno=%d on: %s\n", status, errno, cmd);
        send_mail("ERROR: %s. We must be able to rm these files\n", cmd);
        if(WIFEXITED(status)) {
          if(mvstat = WEXITSTATUS(status)) {  //status ret by rm
            printk("**ERROR: rm exit status = %d\n", mvstat);
          }
        }
        printk("**Continue after ERROR on rm of tlm file\n");
        continue;
      }
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

/* This is called from the main loop to check if any .parc files are in
 * the pipeline to soc dir.
 * The .parc file is sent to /dds/pipe2soc/iris 
 * every time the pipeline back end system does a tapearc.
 * Actually by this cron job on d02:
 *   0 0 * * * /home/production/cvs/JSOC/base/sums/scripts/build_parc_file.pl
 * The .parc file
 * has info on storage units that were archived successfully by the backend,
 * and so can be marked as such in the data capture sum_main table.
 * The sum_main table is updated for its safe_tape info from the .parc.
 * A .parc file looks like:
 * dcs0.jsoc:/dds/pipe2soc/aia> t AIA_2007_131_11_56.parc
 * #owning_series_name                         tape_id  fn  date
 * VC01_2007_131_11_51_39_0123456789A_FFFFF_00 000000S1 666 2007-04-12 17:15:45
 * VC04_2007_131_11_52_09_0123456789A_FFFFF_00 000000S1 666 2007-04-12 17:15:45
 *
 * storage_unit_name pipeline_tape_id_archived_on tape_fn date
*/
void do_pipe2soc() {
  DIR *dfd;
  struct dirent *dp;
  FILE *fp;
  int ptape_fn, complete;
  char line[128], fname[128], cmd[128];
  char *su_name, *ptape_id, *ptape_date;

  if((dfd=opendir(pipedir)) == NULL) {
    printk("**Can't opendir(%s) to find files\n", pipedir);
    return;
  }
  while((dp=readdir(dfd)) != NULL) {
    if(strstr(dp->d_name, ".parc")) {
      sprintf(fname, "%s/%s", pipedir, dp->d_name);
      if(!(fp=fopen(fname, "r"))) {
        printk("***Can't open %s\n", fname);
        continue;
      }
      printk("Found parc file: %s\n", fname);
      complete = 1;
      while(fgets(line, 128, fp)) {       // get .parc file lines 
        if(line[0] == '#' || line[0] == '\n') continue;
        printk("%s", line);
        su_name = (char *)strtok(line, " ");
        ptape_id = (char *)strtok(NULL, " ");
        //ptape_fn = atoi((char *)strtok(NULL, " "));
        errno = 0;
        ptape_fn = (int)strtol((char *)strtok(NULL, " "), (char **) NULL, 10);
        if(errno) continue;	//!!TBD ck what to do here
        ptape_date = (char *)strtok(NULL, "\n");
        if(SUMLIB_SafeTapeUpdate(su_name,ptape_id,ptape_fn,ptape_date)) {
          printk("**ERROR in SUMLIB_SafeTapeUpdate(%s...)\n", su_name);
          complete = 0;
        }
      }
      fclose(fp);
      if(complete) {
        sprintf(cmd, "/bin/rm -f %s", fname);
      }
      else {
        sprintf(cmd, "/bin/mv -f %s %s/err/", fname, pipedir);
      }
      printk("%s\n", cmd);
      system(cmd);
    }
  }
  closedir(dfd);
}

// Initial setup stuff called when main is first entered.
void setup()
{
  FILE *fp;
  int i;
  char string[128], cwdbuf[128], idstr[256];
  char envfile[100], s1[256],s2[256],s3[256], line[256];
  ThreadSigErr_t error = kThreadSigErr_Success;

  //sprint_time_ISOX(dtime,CURRENT_SYSTEM_TIME);
  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  do_datestr();
  printk_set(h0log, h0log);	// set for printk calls 
  printk("%s\n", datestr);
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
  if(grounddata)
    sprintf(string, "ingest_lev0_irisdc_jim -g started as pid=%d user=%s\n", getpid(), username);
  else
    sprintf(string, "ingest_lev0_irisdc_jim [-A] started as pid=%d user=%s\n", getpid(), username);
  strcat(idstr, string);
  printk("*%s", idstr);
  if(restartflg) printk("-r ");
  sprintf(argvc, "vc=%s", vc);
  sprintf(argindir, "indir=%s", tlmdir);
  sprintf(arglogfile, "logfile=%s", logname);
  if(outdir) {
    sprintf(argoutdir, "outdir=%s", outdir);
    printk("%s %s %s %s\n", argvc, argindir, argoutdir, arglogfile);
  }
  else {
    printk("%s %s %s\n", argvc, argindir, arglogfile);
  }
  if(pipedir) {
    sprintf(argpipedir, "pipedir=%s", pipedir);
    printk("Also check for .parc files in %s\n", argpipedir);
  } 
  strcpy(pchan, vc);		// virtual channel primary 
  sprintf(stopfile, "/usr/local/logs/lev0/%s_stop", pchan);
  sprintf(string, "/bin/rm -f %s", stopfile);	//remove any stop file
  system(string);
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
    //if(grounddata) {
    //  sprintf(tlmseriesname, "%s", TLMSERIESNAMEAIAGND);
    //  sprintf(lev0seriesname, "%s", LEV0SERIESNAMEAIAGND);
    //}
    //else {
    //  sprintf(tlmseriesname, "%s", TLMSERIESNAMEAIA);
    //  sprintf(lev0seriesname, "%s", LEV0SERIESNAMEAIA);
    //}
  }
  else {
    if(grounddata) {
      //sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMIGND);
      //sprintf(lev0seriesname, "%s", LEV0SERIESNAMEHMIGND);
    }
    else {
      if(amesflg) {
        sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMIAMES);
        sprintf(lev0seriesname, "%s", LEV0SERIESNAMEHMIAMES);
      }
      else {
        sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMI);
        sprintf(lev0seriesname, "%s", LEV0SERIESNAMEHMI);
      }
    }
  }
  sprintf(bld_vers, "%s", jsoc_version);
  umask(002);			// allow group write 
  Image.initialized = 0;	// init the two image structures 
  ImageOld.initialized = 0;
  Img = &Image;
  ImgO = &ImageOld;
  //imgdecode_init_hack(Img);
  //imgdecode_init_hack(ImgO);
  imgdecode_iris_init_hack(Img);
  imgdecode_iris_init_hack(ImgO);

  //set environment variables for hk code
  //create filename and path
  if(grounddata) strcpy(envfile, ENVFILE_GND );
  else {
    if(amesflg) strcpy(envfile, ENVFILEA);
    else strcpy(envfile, ENVFILE);
  }
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
  if(!restartflg) {
    printk("tlmseriesname=%s\nlev0seriesname=%s\nispseriesname=%s\n", 
		tlmseriesname, lev0seriesname, getenv("HK_ISP_IRIS_DSNAME"));
  }

  //ingest_lev0_test -CARL
  if(logflg) {	//keep wdlog and links to images
    sprintf(wdlogname,"%s_%s",WDLOGDIR,pchan);
    printf("making directory: %s\n", wdlogname);
    mkdir(wdlogname,0777);
    sprintf(wdlogname,"%s_%s/%s",WDLOGDIR,pchan,WDLOGFILE);
    sprintf(oldwdlogname,"%s_%s/%s.old",WDLOGDIR,pchan,WDLOGFILE);
    rename(wdlogname,oldwdlogname);
    if((wdlogfp=fopen(wdlogname,"w")) == NULL) 
    {
      h0log("Can't open the log file %s\n", wdlogname);
    }
    else
    {
      fprintf(wdlogfp, "#image file of lev0 ds by hmi_lev0 %s\n", datestr);
      sprintf(Indir,"%s_%s/%s",WDLOGDIR,pchan,datestr);
      if(mkdir(Indir,0775))
      {
        h0log("Can't mkdir %s. Proceed without it\n",Indir);
      }
    } 
    //fclose(wdlogfp);
  }
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
  vc = (char *)cmdparams_get_str(&cmdparams, "vc", NULL);
  tlmdir = (char *)cmdparams_get_str(&cmdparams, "indir", NULL);
  outdir = (char *)cmdparams_get_str(&cmdparams, "outdir", NULL);
  pipedir = (char *)cmdparams_get_str(&cmdparams, "pipedir", NULL);
  logfile = (char *)cmdparams_get_str(&cmdparams, "logfile", NULL);
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
  if (strcmp(pipedir, NOTSPECIFIED) == 0) {
    pipedir = NULL;
  }
  else {
    if(DS_ConnectDB_Q(dbname)) {	//needed by do_pipe2soc()
      printk("**Can't connect to DB %s\n", dbname);
      printk("**ERROR: Will not be able to process .parc files\n");
      pipedir = NULL;
    }
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
    if(pipedir) do_pipe2soc();  // get .parc file from backend & update db
    if((stat(stopfile, &stbuf) == 0) || abortflg) {
      printk("Abort or Found file: %s. Terminate.\n", stopfile);
      if(abortflg) send_mail("Abort for ingest_lev0_iris for %s\n", pchan);
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
      drms_server_end_transaction(drms_env, 0 , 0); //commit
      //NEW 1Jun2010: leave the stopfile. 
      //Must use doingestlev0_HMI(AIA).pl to start ingest_lev0
      //sprintf(callcmd, "/bin/rm -f %s", stopfile);
      //system(callcmd);
      sprintf(callcmd, "touch /usr/local/logs/lev0/%s_exit", pchan);
      system(callcmd);		//let the world know we're gone
      wflg = 0; //leave DoIt()
      continue;
    }
    sleep(sleep_interval);	//normally 2 sec
    if(cntsleeps == 0) {	//file was seen by do_ingest()
      if(paused) {		//send resume data flow msg
        paused = 0;
        send_mail("tlm files seen again for ingest_lev0_iris for %s\n", pchan);
      }
    }
    else {			//cntsleeps = #of 2sec intervals w/o new file
      if(cntsleeps >= nofiletimeout) {  //flush out anything that we have
        //now close any open image
        if(Image.initialized) {
          if(rs) {		//make sure have a created record
            printf("Data flow stopped for %d sec. Closing current image.\n", 
			nofiletimeout*2);
            printk("Data flow stopped for %d sec. Closing current image.\n", 
			nofiletimeout*2);
            RSISP = RSISPTO;    //use the timeout *rs
            fsnISPTOCLOSE = fsnISP;
            fsnISP = fsnISPX;	//inc last isp
            close_image(rs, segment, segArray, &Image, fsn_prev);
            timeoutclose = 1;
            drms_server_end_transaction(drms_env, 0 , 0); //commit
            drms_server_begin_transaction(drms_env); //start another cycle
            fsn_prev = 0; 	//start over for next file
          }
        }
        else {
          //printk("*No image to close on data timeout\n");
        }
      }
    }
    cntsleeps++;		//#of 2sec sleeps
    if(cntsleeps > 300) {	// >600sec w/o any files
      if(!paused) {
        //send_mail("No files seen for ingest_lev0_iris for %s for 600sec\n", pchan);
        paused = 1;
      }
    }
  }
  if(pipedir) DS_DisConnectDB_Q();
  return(0);
}

