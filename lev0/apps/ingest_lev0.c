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
#include "add_small_image.c"

#define RESTART_CNT 2	//#of tlm files to process before restart

#define LEV0SERIESNAMEHMI "su_production.lev0_test"
#define TLMSERIESNAMEHMI "su_production.tlm_test"
#define LEV0SERIESNAMEAIA "su_production.lev0_test_aia"
#define TLMSERIESNAMEAIA "su_production.tlm_test_aia"

//#define LEV0SERIESNAMEHMI "hmi.lev0_60d"
//#define TLMSERIESNAMEHMI "hmi.tlm_60d"
//#define LEV0SERIESNAMEAIA "aia.lev0_60d"
//#define TLMSERIESNAMEAIA "aia.tlm_60d"

//When change to these data series below to save real data.
//#define TLMSERIESNAMEHMI "hmi.tlm"
//#define LEV0SERIESNAMEAIA "aia.lev0"
//#define LEV0SERIESNAMEAIA "aia.lev0a"
//#define TLMSERIESNAMEAIA "aia.tlm"
//#define LEV0SERIESNAMEHMI "hmi.lev0"
//Also, change setting in $JSOCROOT/proj/lev0/apps/SOURCE_ENV_HK_DECODE file to:
//setenv HK_LEV0_BY_APID_DATA_ID_NAME      lev0
//setenv HK_DF_HSB_DIRECTORY               /tmp21/production/lev0/hk_hsb_dayfile
//setenv HK_JSVNMAP_DIRECTORY              /home/production/cvs/TBL_JSOC/lev0/hk_jsn_map_file/prod


#define H0LOGFILE "/usr/local/logs/lev0/ingest_lev0.%s.%s.%s.log"
#define PKTSZ 1788		//size of VCDU pkt
#define MAXFILES 8192		//max # of file can handle in tlmdir
#define NUMTIMERS 8		//number of seperate timers avail
#define IMAGE_NUM_COMMIT 12	//number of complete images until commit
#define TESTAPPID 0x199		//appid of test pattern packet
#define TESTVALUE 0xc0b		//first value in test pattern packet
#define MAXERRMSGCNT 10		//max # of err msg before skip the tlm file
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define ENVFILE "/home2/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DECODE"

extern int decode_next_hk_vcdu(unsigned short *tbuf, CCSDS_Packet_t **hk, unsigned int *Fsn);
extern int write_hk_to_drms();
extern void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg);
extern int set_HMI_mech_values(DRMS_Record_t *rec);
TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "vc", NOTSPECIFIED, "Primary virt channel to listen to"},
  {ARG_STRING, "indir", NOTSPECIFIED, "directory containing the files to ingest"},
  {ARG_STRING, "outdir", NOTSPECIFIED, "directory to move the files to after the ingest"}, 
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_END}
};

CmdParams_t cmdparams;

// Module name presented to DRMS. 
char *module_name = "ingest_lev0";

FILE *h0logfp;		// fp for h0 ouput log for this run 
IMG Image, ImageOld;
IMG *Img, *ImgO, *ImgC;
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
int hmiaiaflg;			//0=hmi, 1=aia
int whk_status;
int total_tlm_vcdu;
int total_missing_vcdu;
int errmsgcnt, fileimgcnt;
int file_processed_count = 0;	// ingest_lev0 restarts when = RESTART_CNT
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when ingest_lev0 is called for a restart
int sigalrmflg = 0;             // set on signal so prog will know 
int sigtermflg = 0;             // set on signal so prog will know 
int ALRMSEC = 70;               // seconds for alarm signal 
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


int nice_intro ()
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\ningest_lev0 [-vh] "
	"vc=<virt chan> indir=</dir> [outdir=</dir>] [logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
	"  -r: restart. only used when we restart our selves periodically\n"
	"vc= primary virt channel to listen to e.g. VC02\n"
	"indir= directory containing the files to ingest\n"
	"outdir= optional dir to copy the files to after the ingest\n"
	"logfile= optional log file name. Will create one if not given\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
  restartflg = cmdparams_get_int (&cmdparams, "r", NULL);
  return (0);
}

static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss)
{
static int firstcall = 1;
static TIME sdo_epoch;
if (firstcall)
  {
  firstcall = 0;
  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  }
/* XXX fix build 3/18/2008, arta */
return(sdo_epoch + (TIME)sdo_s + (TIME)(sdo_ss)/65536.0);
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

// Close out an image.
void close_image(DRMS_Record_t *rs, DRMS_Segment_t *seg, DRMS_Array_t *array,
		IMG *img, int fsn)
{
  STAT stat;
  int status, n, k;
  char tlmdsname[128];

  printk("*Closing image for fsn = %u\n", fsn);
  if(imgstat(img, &stat)) {
    printk("**Error on imgstat() for fsn = %u\n", fsn);
  }
  else {
    drms_setkey_short(rs, "DATAMIN", stat.min);
    drms_setkey_short(rs, "DATAMAX", stat.max);
    drms_setkey_short(rs, "DATAMEDN", stat.median);
    drms_setkey_double(rs, "DATAMEAN", stat.mean);
    drms_setkey_double(rs, "DATARMS", stat.rms);
    drms_setkey_double(rs, "DATASKEW", stat.skew);
    drms_setkey_double(rs, "DATAKURT", stat.kurt);
  }
  // CAMERA set by HMI_compute_exposure_times()
  if(hmiaiaflg) {		// except for AIA use telnum 
    drms_setkey_int(rs, "CAMERA", (img->telnum)+1);
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
  drms_setkey_int(rs, "COMPID", n+k);
  drms_setkey_int(rs, "BITSELID", img->R);
  drms_setkey_int(rs, "TOTVALS", img->totalvals);
  drms_setkey_int(rs, "DATAVALS", img->datavals);
  drms_setkey_int(rs, "NPACKETS", img->npackets);
  drms_setkey_int(rs, "NERRORS", img->nerrors);
  drms_setkey_int(rs, "EOIERROR", img->last_pix_err);
  drms_setkey_int(rs, "HEADRERR", img->headerr);
  drms_setkey_int(rs, "OVERFLOW", img->overflow);
  drms_setkey_int(rs, "MISSVALS", img->totalvals - img->datavals);
  snprintf(tlmdsname, 128, "%s[%s]", tlmseriesname, tlmnamekeyfirst);
  drms_setkey_string(rs, "TLMDSNAM", tlmdsname);
  int pts = img->first_packet_time >> 16;
  int ptss = img->first_packet_time & 0xffff;
  TIME fpt = SDO_to_DRMS_time(pts, ptss);
  drms_setkey_double(rs, "IMGFPT", fpt);
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
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit ingest_lev0 w/ status = %d\n", stat);
  if (h0logfp) fclose(h0logfp);
  exit(stat);
}

// Called 60 secs after the last .tlm file was seen.
// Will close any opened image. NOTE: a reprocessed image cannot
// be active. It is always closed at the end of get_tlm().
// Can't do the code here because of re-entrancy problem.
void alrm_sig(int sig)
{
  signal(SIGALRM, alrm_sig);
  sigalrmflg = 1;		//tell main loop of signal
}

// User signal 1 is sent by ilev0_alrm1.pl to tell us to exit.
void usr1_sig(int sig)
{
  printk("%s usr1_sig received\n", datestr);
  sleep(1);
  sigtermflg = 1;              // tell main loop to exit
}

void usr2_sig(int sig)
{
  printk("%s usr2_sig received\n", datestr);
  sleep(1);
  sigtermflg = 1;              // tell main loop to exit
}

void sighandler(int sig)
{
  printk("%s signal received\n", datestr);
  sigtermflg = 1;		//tell main loop
  return;
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
  int dstatus;

  // start a new image 
  Img->initialized = 0;
  Img->reopened = 0;
  errmsgcnt = 0;
  fsn_prev = fsnx;
  sprintf(tlmnamekeyfirst, "%s", tlmnamekey);	//save for TLMDSNAM
  rs = drms_create_record(drms_env, lev0seriesname, DRMS_PERMANENT, &dstatus);
  if(dstatus) {
    printk("**ERROR: Can't create record for %s fsn=%u\n", lev0seriesname, fsnx);
    return(1);
  }
  dstatus = drms_setkey_int(rs, "FSN", fsnx);
  if(!(segment = drms_segment_lookup(rs, "image"))) {
    printk("No drms_segment_lookup(rs, image)\n");
    return(1);
  }
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
      //return(1);          // !!!TBD
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
      Img->luid = drms_getkey_int(rs, "LUTID", &rstatus);
      if(rstatus) {
        Img->luid = 0;	//!!!TEMP try this
        printk("ERROR on getkey_int for LUTID\n");
      }
      Img->tap = drms_getkey_int(rs, "TAPCODE", &rstatus);
      compid = drms_getkey_int(rs, "COMPID", &rstatus);
      k = compid & 0x07;		//low 3 bits
      n = compid >> 3;			//next 5 bits
      n = n & 0x1F;
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
  int rstatus, dstatus;
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
    //return(1);		// !!!TBD 
  }
  if(!rset || (rset->n == 0) || rstatus) {
    printk("No prev ds\n");	// start a new image 
    ImgO->initialized = 0;
    ImgO->reopened = 0;
    sprintf(tlmnamekeyfirst, "%s", tlmnamekey);	//save for TLMDSNAM
    rsc = drms_create_record(drms_env, lev0seriesname, DRMS_PERMANENT, &rstatus);
    if(rstatus) {
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
      return(1);		// !!!TBD ck 
    }
    drms_close_records(rset, DRMS_FREE_RECORD);
    rstatus = drms_setkey_int(rsc, "FSN", fsnx);
    ImgO->telnum = drms_getkey_int(rsc, "CAMERA", &rstatus);
    ImgO->telnum--;
    ImgO->cropid = drms_getkey_int(rsc, "CROPID", &rstatus);
    ImgO->luid = drms_getkey_int(rsc, "LUTID", &rstatus);
    ImgO->tap = drms_getkey_int(rsc, "TAPCODE", &rstatus);
    ImgO->N = drms_getkey_int(rsc, "N", &rstatus);
    ImgO->K = drms_getkey_int(rsc, "K", &rstatus);
    ImgO->R = drms_getkey_int(rsc, "R", &rstatus);
    ImgO->overflow = drms_getkey_int(rsc, "OVERFLOW", &rstatus);
    ImgO->headerr = drms_getkey_int(rsc, "HEADRERR", &rstatus);
    ImgO->totalvals = drms_getkey_int(rsc, "TOTVALS", &rstatus);
    ImgO->datavals = drms_getkey_int(rsc, "DATAVALS", &rstatus);
    ImgO->npackets = drms_getkey_int(rsc, "NPACKETS", &rstatus);
    ImgO->nerrors = drms_getkey_int(rsc, "NERRORS", &rstatus);
    ImgO->last_pix_err = drms_getkey_int(rsc, "EOIERROR", &rstatus);
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
      if(rs) 
        close_image(rs, segment, segArray, Img, fsn_prev);
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
        return(1);
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
    if(appid == TESTAPPID) {	// appid of test pattern 
      if(errmsgcnt++ < MAXERRMSGCNT) {
        printk("*Test ApID of %0x found for IM_PDU Cntr = %lld\n", 
			TESTAPPID, vcdu_seq_num);
        for(i=0, j=TESTVALUE; i < 877; i=i+2, j++) {
          datval = MDI_getshort(cbuf+32+i);	// next data value 
          if(datval != j) {
            printk("*Test data value=%0x, expected=%0x for IM_PDU Cntr=%lld\n", 
		datval, j, vcdu_seq_num);
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
            continue;			//get next vcdu
          }
        }
      }
      else {			// continuing normal stream
        if(fsnx != fsn_prev) {          // the fsn has changed 
          if(fsn_change_normal()) {	//handle old & new images
            continue;			//get next vcdu
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
      //printk("$$$$$ appid for Carl =  0x%x %d\n", appid, appid); //!!TEMP
      decode_status = decode_next_hk_vcdu((unsigned short *)(cbuf+10), &Hk, &Fsn);
      switch (decode_status) {
        case SUCCESS_HK_NEED_TO_WTD_CTD:
          printk("*ISP found for fsn = %u\n", Fsn);
          fsnx = Fsn;
          if(rexmit || higherver) {
            if(fsnx != fsn_prev) {          // the fsn has changed
              if(fsn_change_rexmit()) {	//handle old & new images
                continue;			//get next vcdu
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
                continue;			//get next vcdu
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
      imagecnt++;
      fileimgcnt++;
      //fsn_prev = fsn_pre_rexmit; // restore orig for next normal .tlm file
      fsn_prev = 0;		//need for restart of normal mode
      rsc = 0;
    }
  }
  ftmp = EndTimer(1);
  printk("**Processed %s\n**with %d images and %d VCDUs in %f sec\n\n",
	file, fileimgcnt, fpkt_cnt, ftmp);
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
  int i, j, status;
  int rexmit, higherversion;
  char path[DRMS_MAXPATHLEN];
  char name[128], line[128], tlmfile[128], tlmname[96];
  char cmd[128], xxname[128], vername[16];
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
    }
  }
  closedir(dfd);
  qsort(nameptr, i, sizeof(NAMESORT), &compare_names);

  for(j=0; j < i; j++) {
    //printk("####QSORT FILES: %s\n", nameptr[j].name); // !!TEMP 
    // NOTE: the dsf files stay in the indir for now 
    // Currently the cron job pipefe_rm does this:
    // `/bin/mv $dsfname /dds/socdc/hmi/dsf`
    /*******************
    if(strstr(nameptr[j].name, ".dsf")) {
      sprintf(cmd, "/bin/mv %s/%s %s", tlmdir, nameptr[j].name, outdir);
      printk("*mv dsf file to %s\n", outdir);
      printk("%s\n", cmd);
      if(system(cmd)) {
        printk("***Error on: %s\n", cmd);
      }
    }
    ********************/
    if(!strstr(nameptr[j].name, ".qac")) {	// can be .qac or .qacx 
      free(nameptr[j].name);
      continue;
    }
    rexmit = higherversion = 0;
    if(strstr(nameptr[j].name, ".qacx")) {  	// this is a rexmit file 
      rexmit = 1;
    }
    sprintf(name, "%s/%s", tlmdir, nameptr[j].name);
    printk("\n*Found qac file:\n* %s\n", name);
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
    alarm(ALRMSEC);		//restart alarm if no more files come
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
      printk("***Can't create record for %s\n", tlmseriesname);
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
      printk("**ERROR: %d on: %s\n", status, cmd);
    }
    sprintf(cmd, "/bin/mv %s %s", tlmfile, path);
    printk("*mv tlm to %s\n", path);
    printk("%s\n", cmd);
    if(status = system(cmd)) {
      printk("**ERROR: %d on: %s\n", status, cmd);
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
    if(get_tlm(xxname, rexmit, higherversion)) { // lev0 extraction of image 
      printk("***Error in lev0 extraction for %s\n", xxname);
    }
    if(stat(stopfile, &stbuf) == 0) { break; } //signal to stop
    //if(++file_processed_count >= RESTART_CNT) break; 
    //if(sigtermflg) break;
  }
  free(nameptr);
}

// Initial setup stuff called when main is first entered.
void setup()
{
  FILE *fp;
  int i;
  time_t tval;
  struct tm *t_ptr;
  char string[128], cwdbuf[128], idstr[256];
  char envfile[100], s1[256],s2[256],s3[256], line[256];

  signal(SIGALRM, alrm_sig);
  //signal(SIGUSR1, &usr1_sig);	//handle signal 16 sent by ilev0_alrm1.pl
  signal(SIGUSR2, &usr2_sig);	//handle signal 16 sent by ilev0_alrm1.pl
  if (signal(SIGINT, SIG_IGN) != SIG_IGN)
    signal(SIGINT, sighandler);
  if (signal(SIGTERM, SIG_IGN) != SIG_IGN)
    signal(SIGTERM, sighandler);

  tval = time(NULL);
  if (signal(SIGINT, SIG_IGN) != SIG_IGN)
    signal(SIGINT, sighandler);
  if (signal(SIGTERM, SIG_IGN) != SIG_IGN)
    signal(SIGTERM, sighandler);

  tval = time(NULL);
  t_ptr = localtime(&tval);
  sprintf(datestr, "%d.%02d.%02d_%02d:%02d:%02d", 
	  (t_ptr->tm_year+1900), (t_ptr->tm_mon+1),
	  t_ptr->tm_mday, t_ptr->tm_hour, t_ptr->tm_min, t_ptr->tm_sec);
  printk_set(h0log, h0log);	// set for printk calls 
  printk("%s\n", datestr);
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
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
  sprintf(stopfile, "/usr/local/logs/lev0/%s_stop", pchan);
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
    sprintf(tlmseriesname, "%s", TLMSERIESNAMEAIA);
    sprintf(lev0seriesname, "%s", LEV0SERIESNAMEAIA);
  }
  else {
    sprintf(tlmseriesname, "%s", TLMSERIESNAMEHMI);
    sprintf(lev0seriesname, "%s", LEV0SERIESNAMEHMI);
  }
  umask(002);			// allow group write 
  Image.initialized = 0;	// init the two image structures 
  ImageOld.initialized = 0;
  Img = &Image;
  ImgO = &ImageOld;
  Img->initialized = 0;
  ImgO->initialized = 0;

  //set environment variables for hk code
  //create filename and path
  strcpy(envfile, ENVFILE );
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
  //close file
  fclose(fp);
}

// Module main function. 
int DoIt(void)
{
  pid_t pid;
  int wflg = 1;
  char *args[6];
  char callcmd[128];

  if (nice_intro())
    return (0);
  if(!(username = (char *)getenv("USER"))) username = "nouser"; 
  vc = cmdparams_get_str(&cmdparams, "vc", NULL);
  tlmdir = cmdparams_get_str(&cmdparams, "indir", NULL);
  outdir = cmdparams_get_str(&cmdparams, "outdir", NULL);
  logfile = cmdparams_get_str(&cmdparams, "logfile", NULL);
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
    sleep(12);		//make sure old ingest_lev0 is done
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
/****************************************************************************
    //note alrm sig doesn't work for multi-threaded process
    if(sigalrmflg) {		// process an alarm timout for no data in
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
      sigalrmflg = 0;
    }
**************************************************************************/

/************************************************************************
    if(file_processed_count >= RESTART_CNT) {	//restart this ingest
      if((pid = fork()) < 0) {
        printk("***Can't fork(). errno=%d\n", errno);
      }
      else if(pid == 0) {                   // this is the beloved child 
        printk("\nexecvp of ingest_lev0\n");
        args[0] = "ingest_lev0";
          args[1] = "-r";
          args[2] = argvc;
          args[3] = argindir;
          args[4] = arglogfile;
          if(outdir) {
            args[5] = argoutdir;
            args[6] = NULL;
          } else {
            args[5] = NULL;
          }
        if(execvp(args[0], args) < 0) {
          printk("***Can't execvp() ingest_lev0. errno=%d\n", errno);
          exit(1);
        }
      }
      //now parent closed any open image
      if(Image.initialized) {
        if(rs) {		//make sure have a created record
          close_image(rs, segment, segArray, &Image, fsn_prev);
          printk("*Closed image on restart FSN=%u\n", fsn_prev);
        }
      }
      else {
        printk("*No image to close on restart\n");
      }
      printk("restart: drms_server_end_transaction()\n");
      drms_server_end_transaction(drms_env, 0 , 0); //commit
      wflg = 0;		//exit DoIt
      continue;
    }
**************************************************************************/

//    if(file_processed_count >= RESTART_CNT) {	//restart this ingest
//      //now parent closed any open image
//      if(Image.initialized) {
//        if(rs) {		//make sure have a created record
//          close_image(rs, segment, segArray, &Image, fsn_prev);
//          printk("*Closed image on restart FSN=%u\n", fsn_prev);
//        }
//      }
//      else {
//        printk("*No image to close on restart\n");
//      }
//      /***************************************************************
//      printk("restart: drms_server_end_transaction()\n");
//      drms_server_end_transaction(drms_env, 0 , 0); //commit left overs
//      **********************************************************/
//      sprintf(callcmd, "ingest_lev0 -r %s %s %s &", argvc,argindir,arglogfile);
//      if(system(callcmd)) {
//        printk("**ERROR on: %s\n", callcmd);
//      }
//      wflg = 0;		//leave DoIt()
//      continue;
//    }
//!!!TEMP test
    if(stat(stopfile, &stbuf) == 0) {
      printk("Found file: %s. Terminate.\n", stopfile);
      //now close any open image
      if(Image.initialized) {
        if(rs) {		//make sure have a created record
          close_image(rs, segment, segArray, &Image, fsn_prev);
        }
      }
      else {
        printk("*No image to close on exit\n");
      }
      sprintf(callcmd, "/bin/rm -f %s", stopfile);
      system(callcmd);
      sprintf(callcmd, "touch /usr/local/logs/lev0/%s_exit", pchan);
      system(callcmd);		//let the world know we're gone
      wflg = 0; //leave DoIt()
      continue;
    }
    //detect if exit. NOTE: do_ingest() can take long time if lots of files
    //!!NOTE: it may be the DRMS server that gets the ^C
    //man signal say undefined for multi-threaded process
    if(sigtermflg) {
      wflg = 0; //leave DoIt()
      continue;
    }
    sleep(2);
  }
  return(0);
}

