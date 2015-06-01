/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev1/apps/build_lev1.c
 *-----------------------------------------------------------------------------
 *
 * This is a module that runs with DRMS and processes lev0
 * filtergrams to lev1.
 * It is scheduled by build_lev1_mgr either by qsub to the cluster
 * or by fork to run on the local machine. It is called with 
 * 12 (NUMRECLEV1 defined in lev0lev1.h) lev0 images at a time.
 * See the -h nice_intro() for details of the calling options.
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
#include "list.h"
#include </home/jsoc/include/fftw3.h> 
//#include "fftw3.h" 
#include "imgdecode.h"
#include "lev0lev1.h"
#include "quallev1.h"
#include "limb_fit.h"
#include "cosmic_ray.h"
#include "get_pointing_info.c"

//default in and out data series
#define LEV0SERIESNAMEIRIS "iris_ground.lev0_dc1"
#define LEV1SERIESNAMEIRIS "iris_ground.lev1_dc1"
//#define DSCRSNAME "iris_ground.crs_table"
#define DSCRSNAME "iris.crs_table"
//#define DSCRSNAME "iris_ground.window_table"
#define DSFFNAME "iris.flatfield_prelim"
//#define DSFFNAME "iris_ground.flatfield_test"

//#define DSFFNAME "su_richard.flatfield"		//temp test case
//#define DSFFNAMEHMI "su_production.hmi_flatfield"	//temp test case
//#define DSFFNAMEHMI "hmi.flatfield"
//#define DSFFNAMEAIA "su_production.aia_flatfield"	//temp test case
//#define DSFFNAMEAIA "aia_test.flatfield"	//temp test case
//#define DSFFNAMEAIA "aia.flatfield"

#define LEV1LOG_BASEDIR "/usr/local/logs/lev1"
#define H1LOGFILE "/usr/local/logs/lev1/build_lev1_iris.%s.log"
#define NUMTIMERS 8		//number of seperate timers avail
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define LOGTEST 0
#define CAL_HCFTID 17		//image is cal mode 
#define STOP_FILE "/usr/local/logs/lev1/build_mgr_stop_iris"

// SAA-HLZ
#define SAA_HLZ_SERIES "iris.saa_hlz"

int compare_rptr(const void *a, const void *b);
static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "mode", NOTSPECIFIED, "either recnum or fsn"},
  {ARG_STRING, "dsin", NOTSPECIFIED, "dataset of lev0 filtergrams"},
  {ARG_STRING, "dsout", NOTSPECIFIED, "dataset of lev1 output"},
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_INTS, "brec", "0", "first lev0 rec# to process. 0=error must be given by build_lev1_mgr"},
  {ARG_INTS, "erec", "0", "last lev0 rec# to process. 0=error must be given by build_lev1_mgr"},
  {ARG_INTS, "bfsn", "0", "first lev0 fsn# to process. 0=error must be given by build_lev1_mgr"},
  {ARG_INTS, "efsn", "0", "last lev0 fsn# to process. 0=error must be given by build_lev1_mgr"},
  {ARG_INTS, "quicklook", "1", "1=quick look, 0 = definitive mode"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_END}
};

CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "build_lev1_iris_NEW";

FILE *h1logfp;		// fp for h1 ouput log for this run 
static struct stat stbuf;
static LEV0LEV1 lev0lev1;
static LEV0LEV1 *l0l1 = &lev0lev1;
//static CCSDS_Packet_t *Hk;
static DRMS_Record_t *rs;
static DRMS_Record_t *crsrec;
static DRMS_Record_t *rs0, *rs1, *rsff, *rsbad_pix, *rec_bad_aia, *rt, *rresp;
static DRMS_Record_t *rptr;
static DRMS_Segment_t *segment, *segment0;
static DRMS_Segment_t *segmentff;
static DRMS_Segment_t *darkseg;
static DRMS_Segment_t *badseg;
//static DRMS_Segment_t *badoutpixseg;
//static DRMS_Segment_t *spikeseg;
static DRMS_Array_t *segArray;
static DRMS_RecordSet_t *crsset;
static DRMS_RecordSet_t *rset0, *rset1, *rsetff, *rsbad_aia, *rs_t=NULL,
                        *rs_resp=NULL;
static DRMS_Array_t *Array0;
static DRMS_Array_t *Arrayff;
static DRMS_Array_t *ArrayDark;
static DRMS_Array_t *ArrayBad;
static DRMS_Array_t *ArraySpike;
static TIME sdo_epoch;
static PTINFO *ptinfo = NULL;
static PTINFO ptdata;
static char bld_vers[16];
static char datestr[32];
static char open_dsname[256];
static char dsffname[128];
static char path[DRMS_MAXPATHLEN], bad_pix_path[DRMS_MAXPATHLEN];
static char bad_aia_path[DRMS_MAXPATHLEN];
static char rs1_path[DRMS_MAXPATHLEN];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];
//static char *orbseries = "sdo.fds_orbit_vectors";
//static char *orbseries = "su_arta.orbitvectors";
static char *orbseries = "iris.orbit_vectors";

static int nspikes, respike, fid, aiftsid, *oldvalues, *spikelocs, *newvalues;
static int hcftid, aiagp6;
static short aifcps;
double aiascale = 1.0;

IORBIT_Info_t *IOinfo = NULL;

IORBIT_Info_t IOdata;
LIBASTRO_Error_t IOstatus;
unsigned int fsnarray[NUMRECLEV1];
unsigned int fsnx = 0;
//short data1[MAXPIXELS];		//iris lev1 is shorts
//int data1[MAXPIXELS];
float data1[MAXPIXELS];		//floats for HMI
float ftmp;
int data1A[MAXPIXELS];		//ints for AIA
short data1S[MAXPIXELS];	//shorts for iris
int array_cosmic[16777216];	//4096x4096
double tgttimes[NUMRECLEV1];

long long brec, erec, bfsn, efsn; //begin and end lev0 rec/fsn. must be same data type
long long bnumx, enumx;		  //has either the brec/erec of bfsn/efsn pair accord to mode
int verbose;
int modeflg = 0;		//0=fsn, 1=recnum
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when build_lev1 is called for a restart
int abortflg = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
int quicklook;
//global quality flags
int flatmiss[NUMRECLEV1];
int orbmiss[NUMRECLEV1];
int asdmiss[NUMRECLEV1];
int mpdmiss[NUMRECLEV1];
int limbmiss[NUMRECLEV1];
int noimage[NUMRECLEV1];
int missflg[NUMRECLEV1];

char logname[128];
char argdsin[128], argdsout[128], arglogfile[128];
char argbx[80], argex[80], argquick[80], argmode[80];
char timetag[32];
char tlmseriesname[128];	// e.g. hmi.tlm
char lev0seriesname[128];	// e.g. hmi.lev0
char *username;			// from getenv("USER") 
char *logfile;			// optional log name passed in 
char *mode;			// given mode. recnum or fsn
char *dsin;			// lev0 input dataset
char *dsout;			// lev1 output dataset
char *dsff;                     // flat field dataset

//!!TEMP
typedef struct {
  float rsun_lf;
  float x0_lf;
  float y0_lf;
} LIMB_SOMETHING;



//int get_nspikes() { return nspikes; }
//int get_respike(void) { return respike; }
//int *get_spikelocs() { return spikelocs; }
//int *get_oldvalues() { return oldvalues; }
//int *get_newvalues() { return newvalues; }
//void set_nspikes(int new_nspikes) { nspikes = new_nspikes; }
//void set_spikelocs(int *new_spikelocs) { spikelocs = new_spikelocs; }
//void set_oldvalues(int *new_oldvalues) { oldvalues = new_oldvalues; }
//void set_newvalues(int *new_newvalues) { newvalues = new_newvalues; }

//Set the QUALITY keyword for lev1
void do_quallev1(DRMS_Record_t *rs0, DRMS_Record_t *rs1, int inx, unsigned int fsn)
{
  int quallev1 = 0;
  int rstatus;
  char *pchar;

  quallev1 = missflg[inx];
  if(flatmiss[inx]) quallev1 = quallev1 | Q_NOFLAT;
  if(orbmiss[inx]) quallev1 = quallev1 | Q_NOORB;
  if(limbmiss[inx]) quallev1 = quallev1 | Q_NOLIMB;
  if(asdmiss[inx]) {
    quallev1 = quallev1 | Q_NOASD;
    drms_setkey_string(rs1, "ASD_REC", DRMS_MISSING_STRING);
  }
  if(mpdmiss[inx]) quallev1 = quallev1 | Q_NOMPD;
  if(noimage[inx]) quallev1 = quallev1 | Q_MISSALL;
  if(ptinfo) {
    ptdata = ptinfo[inx];
    if(strcmp(ptdata.acs_mode, "SCIENCE")) { //not in science pointing mode
      quallev1 = quallev1 | Q_NOACS_SCI;
    }
    if(strcmp(ptdata.acs_eclp, "NO")) {      //in eclipse
      quallev1 = quallev1 | Q_ACS_ECLP;
    }
    if(strcmp(ptdata.acs_sunp, "YES")) {     //not in sun presence flag
      quallev1 = quallev1 | Q_ACS_SUNP;
    }
    if(strcmp(ptdata.acs_safe, "NO")) {     //safemode set
      quallev1 = quallev1 | Q_ACS_SAFE;
    }
  }
  pchar = drms_getkey_string(rs0, "IMG_TYPE", &rstatus);
  if(rstatus) {
    printk("ERROR: in drms_getkey_string(IMG_TYPE) fsn=%u\n", fsn);
  }
  else {
    if(strcmp(pchar, "LIGHT")) {    //dark image
      quallev1 = quallev1 | Q_IMG_TYPE;
    }
  }
  if(quicklook) {
    quallev1 = quallev1 | Q_NRT;
  }
  drms_setkey_int(rs1, "QUALITY", quallev1);
  drms_setkey_string(rs1, "BLD_VERS", bld_vers); //build vers to every record
}

int nice_intro ()
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\nbuild_lev1_iris [-vhr] "
	"mode=<recnum|fsn> dsin=<lev0> dsout=<lev1>\n"
	"brec=<rec#>|bfsn=<fsn#> erec=<rec#>|efsn=<fsn#>\n"
	"quicklook=<0|1> [logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
	"  -r: restart. only used when we restart our selves periodically\n"
        "mode= recnum: brec and erec have the record # range to process \n"
        "      fsn: bfsn and efsn have the fsn # range to process\n"
        "      For safety, the mode and arg name used must be consistent\n"
	"dsin= data set name of lev0 input\n"
	"      default hmi=hmi.lev0e   aia=aia.lev0e\n"
	"dsout= data set name of lev1 output\n"
	"      default hmi=su_production.hmi_lev1e   aia=su_production.aia_lev1e\n"
	"brec= first lev0 rec# to process. 0=error must be given by build_lev1_mgr\n"
	"erec= last lev0 rec# to process. 0=error must be given by build_lev1_mgr\n"
	"bfsn= first fsn# to process. 0=error must be given by build_lev1_mgr\n"
	"efsn= last fsn# to process. 0=error must be given by build_lev1_mgr\n"
	"quicklook= 1 = quicklook mode, 0 = definitive mode\n"
	"logfile= optional log file name. If not given uses:\n"
        "         /usr/local/logs/lev1/build_lev1_aia.<time_stamp>.log\n");
    return(1);
    }
  verbose = cmdparams_get_int (&cmdparams, "v", NULL);
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
  sprintf(cmd, "echo \"%s\" | Mail -s \"build_lev1_iris mail\" lev0_user", string);
  system(cmd);
  va_end(args);
  return(0);
}

// Got a fatal error. 
void abortit(int stat)
{
  printk("***Abort in progress ...\n");
  printk("**Exit build_lev1_iris w/ status = %d\n", stat);
  if (h1logfp) fclose(h1logfp);
  exit(stat);
}

//!!TBD Keh-Cheng
int rdout_mode_correct()
{
  return(0);
}

//!!TBD Art
int orbit_calc()
{
  return(0);
}

static int GetSaaHlz(DRMS_Env_t *env, DRMS_RecordSet_t *rs, double *tgttimes, int ntimes, int (*myprintkerr)(const char *fmt, ...))
{
    IORBIT_SaaHlzInfo_t *saahlzinfo = NULL;
    int itime;
    LinkedList_t *list = NULL;
    LinkedList_t **pList = NULL;
    ListNode_t *node = NULL;
    char *eventType = NULL;
    int nhlz;
    int shlz;
    int saa;
    HContainer_t *colToList = NULL;
    HContainer_t **pColToList = NULL;
    char timeStr[IORBIT_SAAHLZINFO_TIME_KEY_LEN];
    LIBASTRO_Error_t rv;
    int ret = 0;
    
    rv = iorbit_getSaaHlzInfo(env, SAA_HLZ_SERIES, tgttimes, ntimes, &saahlzinfo);
    
    if (rv == kLIBASTRO_Success)
    {
        for (itime = 0; itime < ntimes; itime++)
        {
            nhlz = 0;
            shlz = 0;
            saa = 0;
            
            snprintf(timeStr, sizeof(timeStr), "%lf", tgttimes[itime]);
            
            if ((pColToList = hcon_lookup(saahlzinfo, timeStr)) != NULL)
            {
                colToList = *pColToList;
                
                if ((pList = hcon_lookup(colToList, IORBIT_SAAHLZINFO_KW_EVENT_TYPE)) != NULL)
                {
                    list = *pList;
                    list_llreset(list);
                    
                    while((node = list_llnext(list)) != NULL)
                    {
                        eventType = (char *)node->data;
                        
                        if (strcasecmp(eventType, "NHLZ") == 0)
                        {
                            nhlz = 1;
                        }
                        else if (strcasecmp(eventType, "SHLZ") == 0)
                        {
                            shlz = 1;
                        }
                        else if (strcasecmp(eventType, "SAA") == 0)
                        {
                            saa = 1;
                        }
                    }
                    
                    /* Set SAA-HLZ keywords */
                    if (nhlz && shlz)
                    {
                        /* ERROR - set HLZ to 0. */
                        printkerr("ERROR - For tobs %lf, both event types NHLZ and SHLZ exist.\n", tgttimes[itime]);
                        drms_setkey_int(rs->records[itime], "HLZ", 0);
                    }
                    else if (nhlz)
                    {
                        printf("  Setting HLZ to 1.\n");
                        drms_setkey_int(rs->records[itime], "HLZ", 1);
                    }
                    else if (shlz)
                    {
                        printf("  Setting HLZ to 2.\n");
                        drms_setkey_int(rs->records[itime], "HLZ", 2);
                    }
                    else 
                    {
                        drms_setkey_int(rs->records[itime], "HLZ", 0);                    
                    }
                    
                    if (saa)
                    {
                        drms_setkey_int(rs->records[itime], "SAA", 1);
                    }
                    else
                    {
                        drms_setkey_int(rs->records[itime], "SAA", 0);
                    }
                }
                else
                {
                    myprintkerr("ERROR - SAA-HLZ info for keyword %s unexecpectedly missing.\n", IORBIT_SAAHLZINFO_KW_EVENT_TYPE);
                    ret = 1;
                }
            }
            else
            {
                myprintkerr("ERROR - SAA-HLZ info for tobs %lf unexecpectedly missing.\n", tgttimes[itime]);
                ret = 1;
            }
        }
        
        iorbit_cleanSaaHlzInfo(&saahlzinfo);
    }
    else
    {
        myprintkerr("ERROR - Couldn't query db properly.\n");
        ret = 1;
    }
    
    return ret;
}

//#include "aia_despike.c"
//#include "do_flat.c"
#include "get_image_location.c"
#include "limb_fit_function.c"
#include "cosmic_ray.c"
#include "heightformation.c"

//Called with the range to do. The args are either rec#s or fsn#s accord to modeflg
int do_ingest(long long bbrec, long long eerec)
{
  //FILE *fwt;
  Image_Location *p_imageloc;
  Image_Location imageloc[NUMRECLEV1];
  TIME t_obs0;
  TIME tobs[NUMRECLEV1];
  float percentd;
  float cdelt1, rsun, crpix1, crpix2, crota2;
  double rsun_lf, x0_lf, y0_lf;
  int rstatus, dstatus, lstatus, ncnt, fcnt, i, j, k, qualint, nobs;
  int hshiexp, hcamid, nbad, n_cosmic;
  int *spikedata, status, axes[2], nbytes;
  uint32_t missvals, totvals;
  long long recnum0, recnum1, recnumff;
  char recrange[128], lev0name[128], flatrec[128], temprec[128], pointrec[128];
  //char tmpname[80];
  double scroll;

  if(modeflg) sprintf(recrange, ":#%lld-#%lld", bbrec, eerec);
  else sprintf(recrange, "%lld-%lld", bbrec, eerec);
  sprintf(open_dsname, "%s[%s]", dsin, recrange);
  printk("open_dsname = %s\n", open_dsname);
  printk("#levnum recnum fsn\n");

    t_obs0 = 0;
    rset0 = drms_open_records(drms_env, open_dsname, &rstatus); //open lev0
    if(!rset0 || (rset0->n == 0) || rstatus) {
      printk("Can't do drms_open_records(%s)\n", open_dsname);
      printf("Can't do drms_open_records(%s)\n", open_dsname);
      //return(1);
      return(0);
    }
    drms_stage_records(rset0, 1, 0);
    ncnt = rset0->n;
    rptr = (DRMS_Record_t *)malloc(ncnt * sizeof(DRMS_Record_t));
    if(rptr == NULL) {
      printk("Can't malloc() for DRMS_Record_t sort\n");
      return(1);
    }
    //make opened records sequential in mem for sort
    for(i=0; i < ncnt; i++) {
      memcpy(&rptr[i], rset0->records[i], sizeof(DRMS_Record_t));
    }
    //Must sort to get in ascending t_obs order
    qsort(rptr, ncnt, sizeof(DRMS_Record_t), &compare_rptr);

    //this loop is for the benefit of the lev1view display to show
    //info on all lev0 records opened
    for(i=0; i < ncnt; i++) {
      flatmiss[i] = 0;		//init quality flags
      orbmiss[i] = 0;
      limbmiss[i] = 0;
      asdmiss[i] = 0;
      mpdmiss[i] = 0;
      noimage[i] = 0;
      missflg[i] = 0;
      rs0 = &rptr[i];
      recnum0 = rs0->recnum;
      //must get fsn in case got record from last lev1 record
      fsnx = drms_getkey_int(rs0, "FSN", &rstatus); 
      fsnarray[i] = fsnx;
      printk("*0 %u %u\n", recnum0, fsnx);
      //also set up call for get_pointing_info() and iorbit_getinfo()
      tobs[i] = drms_getkey_time(rs0, "t_obs", &rstatus);
      if(rstatus) {
        printk("Error on drms_getkey_time() fsn=%u. Use DRMS_MISSING_TIME\n", 
		fsnx);
       tobs[i] = DRMS_MISSING_TIME;
      }
    }
    //New from Art's stuff 07Aug2013
    HContainer_t *keymap = NULL;

    keymap = hcon_create(DRMS_MAXKEYNAMELEN, DRMS_MAXKEYNAMELEN, NULL, NULL, NULL, NULL, 0);

    if (!keymap)
    {
        IOstatus = kLIBASTRO_OutOfMemory;
    }
    else {
        hcon_insert(keymap, "kXGCI", "geixobs");
        hcon_insert(keymap, "kYGCI", "geiyobs");
        hcon_insert(keymap, "kZGCI", "geizobs");
        hcon_insert(keymap, "kXHCI", "heixobs");
        hcon_insert(keymap, "kYHCI", "heiyobs");
        hcon_insert(keymap, "kZHCI", "heizobs");
        hcon_insert(keymap, "kRSUNOBS", "rsunobs");
        hcon_insert(keymap, "kOBSVR", "obsvr");
        hcon_insert(keymap, "kDSUNOBS", "dsunobs");
        hcon_insert(keymap, "kOBSDATE", "obsdate");

        IOstatus = iorbit_getinfo_ext(drms_env,
                                orbseries,
                                NULL,
                                IORBIT_Alg_Quadratic,
                                tobs,
                                ncnt,
                                kIORBIT_CacheAction_DontCache,
                                &IOinfo,
                                keymap);
        hcon_destroy(&keymap);
    }
    if(IOstatus != kLIBASTRO_Success) {
      if(IOstatus == kLIBASTRO_InsufficientData) {
        printk("***ERROR in iorbit_getinfo_ext: kLIBASTRO_InsufficientData\n");
      }
      else {
       printk("***ERROR in iorbit_getinfo_ext() status=%d\n", IOstatus);
      }
      for(j=0; j < ncnt; j++) {  //set qual bits
        orbmiss[j] = 1;
      }
      return(1);                //abort. new 2/22/2011
    }
    
    /* IOStatus == kLIBASTRO_Success */
    
    rset1 = drms_create_records(drms_env, ncnt, dsout, DRMS_PERMANENT,&dstatus);
    if(dstatus) {
      printk("**ERROR: Can't create records for %s\n", dsout);
      for(j=0; j < ncnt; j++) {	 //set qual bits
        noimage[j] = 1;
      }
      return(1);	//new 2/22/2011
    }
    
    /*****************************************************************/
    /* Obtain SAA-HLZ information, and set the SAA and HLZ keywords. */
    
    if (GetSaaHlz(drms_env, rset1, tobs, ncnt, printkerr))
    {
        printk("***Error - Unable to fetch SAA-HLZ information.\n");
        return 1;
    }
    
    /* End SAA-HLZ                                                   */
    /*****************************************************************/
    

    for(i=0; i < ncnt; i++) { 	//do for all the sorted lev0 records
      //StartTimer(2);	//!!TEMP
      rs0 = &rptr[i];
      recnum0 = rs0->recnum;
      fsnx = fsnarray[i]; 
      sprintf(lev0name, "%s[%u]", dsin, fsnx);
      if(drms_getkey_int(rs0, "QUALITY", 0) < 0) {
        printk("Bad QUALITY for %s, no lev1 made\n", lev0name);
        noimage[i] = 1;
        //continue;	//make an image anyway so the lev1 fsn will exist
      }
      segment0 = drms_segment_lookupnum(rs0, 0);
      int seg0sz = segment0->axis[0]*segment0->axis[1];
      Array0 = drms_segment_read(segment0, DRMS_TYPE_SHORT, &rstatus);
      if(!Array0) {
        printk("Can't do drms_segment_read() %s status=%d\n", 
			lev0name, rstatus);
        noimage[i] = 1;
        return(1);	//return until we learn
        continue;
      }
      l0l1->adata0 = (short *)Array0->data; //free at end
      l0l1->dat1.adata1A = &data1A; 	    //int out
      l0l1->rs0 = rs0;
      l0l1->recnum0 = recnum0;
      l0l1->fsn = fsnx;
      l0l1->nx = segment0->axis[0];
      l0l1->ny = segment0->axis[1];
      l0l1->datavals = drms_getkey_int(rs0, "DATAVALS", &rstatus);
      l0l1->missvals = drms_getkey_int(rs0, "MISSVALS", &rstatus);

      sprintf(open_dsname, "%s[%u]", dsout, fsnx);
      rs = rset1->records[i]; 

      // 
      // find closest iris.pointing_data record and set scroll
      //
      {
        if(!quicklook) {
	  DRMS_RecordSet_t *rset;
	  DRMS_Record_t *rt;
	  int st;
	  TIME t_obs, time_qbi;
	  float aeulrbrx, aeulrbry, aeulrbrz, acg_roll, ophase;
	  t_obs = drms_getkey_time(rs0, "t_obs", &st);
	  if (st) {
	      t_obs = DRMS_MISSING_TIME;
          }
          else {
	  sprintf(open_dsname, "iris.pointing_data[? date_obs > %f and date_obs <= %f ?]", t_obs-5, t_obs+5); 
          printf("%s\n", open_dsname);
	  rset = drms_open_records(drms_env, open_dsname, &st);
	  if (!rset || !rset->n || st) {
	      printk("Error in drms_open_records(%s); setting scroll to zero\n", open_dsname);
	      scroll = 0.0;
	  } else {
	      // There should be only one record returned
	      rt = rset->records[0];
              sprintf(pointrec, "iris.pointing_data[:#%lld]", rt->recnum);
              if(dstatus = drms_setkey_string(rs, "POINTREC", pointrec)) {
                printk("**ERROR on drms_setkey_string() for %s\n", pointrec);
              }
	      time_qbi = drms_getkey_time(rt, "TIME_QBI", &st);
	      aeulrbrx = drms_getkey_float(rt, "A_EULERBR_X", &st);
	      aeulrbry = drms_getkey_float(rt, "A_EULERBR_Y",  &st);
	      aeulrbrz = drms_getkey_float(rt, "A_EULERBR_Z", &st);
	      acg_roll = drms_getkey_float(rt, "A_CG_ROLL_ANGLE", &st);
	      ophase   = drms_getkey_float(rt, "OPHASE", &st);

	      drms_setkey_time(rs, "TIME_QBI", time_qbi);
	      drms_setkey_float(rs, "AEULRBRX", aeulrbrx);
	      drms_setkey_float(rs, "AEULRBRY", aeulrbry);
	      drms_setkey_float(rs, "AEULRBRZ", aeulrbrz);
	      drms_setkey_float(rs, "SAT_ROT", aeulrbrz);
	      drms_setkey_float(rs, "ACG_ROLL", acg_roll);
	      drms_setkey_float(rs, "OPHASE",   ophase);

	      scroll = aeulrbrz;
	  }

	  if (rset)
	      drms_close_records(rset, DRMS_FREE_RECORD);
          }
        }
      }
      
      drms_record_directory(rs, rs1_path, 0);
      if(!*rs1_path) {
        printk("***ERROR: No path to segment for %s\n", open_dsname);
        noimage[i] = 1;
        continue;
      }
      printf("\npath to lev1 = %s\n", rs1_path);	//!!TEMP
      if(rstatus = iris_isp2wcs(rs0, rs, scroll)) {
        printk("**ERROR: iris_isp2wcs() status = %d\n", rstatus);
        printk("Press on after error at fsn=%u...\n", fsnx);
        printf("**ERROR: at fsn %u\n", fsnx);
        continue;
        //printf("**ERROR: Abort on fsn = %u\n", fsnx);
        //return(1);
      }
      dstatus = drms_setkey_int(rs, "FSN", fsnx);
      //dstatus = drms_setkey_string(rs, "LEV0SERIES", lev0name); //no such keyword
      if(!(segment = drms_segment_lookup(rs, "image_lev1"))) {
        printk("No drms_segment_lookup(rs, image_lev1) for %s\n", open_dsname);
        noimage[i] = 1;
        continue;
      }
        segArray = drms_array_create(DRMS_TYPE_INT,
                                       segment->info->naxis,
                                       segment0->axis,
                                       &data1A,
                                       &dstatus);
//short *adata = (short *)Array0->data;
//short *bdata = &data1S;
//memcpy(bdata, adata, 2*seg0sz); //!!TEMP mv the lev0 in for now

      //transer the lev0 keywords
      rstatus = drms_copykeys(rs, rs0, 0, kDRMS_KeyClass_Explicit);
      if(rstatus != DRMS_SUCCESS) {
        printk("Error %d in drms_copykeys() for fsn %u\n", fsnx);
        return(1);		//new 2/22/2011
        continue;
      }
      qualint = drms_getkey_int(rs0, "QUALITY", &rstatus);
      drms_setkey_int(rs, "QUALLEV0", qualint);
      fid = drms_getkey_int(rs0, "FID", &rstatus);

      short isqisysn = drms_getkey_short(rs0, "ISQISYSN", &rstatus);  
      int iimgots1 = drms_getkey_int(rs0, "IIMGOTS1", &rstatus);
      int iimgots2 = drms_getkey_int(rs0, "IIMGOTS2", &rstatus);
      int iimgots3 = drms_getkey_int(rs0, "IIMGOTS3", &rstatus); 
      switch (isqisysn) {
          case 0:
            drms_setkey_int(rs, "IIMGOTS", iimgots1);
            break;
          case 1:
            drms_setkey_int(rs, "IIMGOTS", iimgots2);
            break;
          case 2:
            drms_setkey_int(rs, "IIMGOTS", iimgots3);
            break;
      } 

      drms_setkey_time(rs, "T_OBS", tobs[i]);
      printk("t_obs for lev0 = %10.5f fsn=%u\n", tobs[i], fsnarray[i]);  //!!TEMP
      drms_setkey_double(rs, "DATE", CURRENT_SYSTEM_TIME);
      if(ptinfo) {
        ptdata = ptinfo[i];
          drms_setkey_float(rs, "SAT_Y0", ptdata.sat_y0);
          drms_setkey_float(rs, "SAT_Z0", ptdata.sat_z0);
          drms_setkey_float(rs, "SAT_ROT", ptdata.sat_rot);
          drms_setkey_string(rs, "ACS_MODE", ptdata.acs_mode);
          drms_setkey_string(rs, "ACS_ECLP", ptdata.acs_eclp);
          drms_setkey_string(rs, "ACS_SUNP", ptdata.acs_sunp);
          drms_setkey_string(rs, "ACS_SAFE", ptdata.acs_safe);
          drms_setkey_string(rs, "ASD_REC", ptdata.asd_rec);
          drms_setkey_string(rs, "ACS_CGT", ptdata.acs_cgt);
      }
      if(IOinfo) {
        IOdata = IOinfo[i];
           drms_setkey_double(rs, "RSUN_OBS", IOdata.rsun_obs);
           drms_setkey_double(rs, "DSUN_OBS", IOdata.dsun_obs);
           drms_setkey_double(rs, "OBS_VR", IOdata.obs_vr);
           drms_setkey_double(rs, "GEIX_OBS", IOdata.gciX);
           drms_setkey_double(rs, "GEIY_OBS", IOdata.gciY);
           drms_setkey_double(rs, "GEIZ_OBS", IOdata.gciZ);
           drms_setkey_double(rs, "HEIX_OBS", IOdata.hciX);
           drms_setkey_double(rs, "HEIY_OBS", IOdata.hciY);
           drms_setkey_double(rs, "HEIZ_OBS", IOdata.hciZ);
           drms_setkey_string(rs, "ORB_REC", IOdata.orb_rec); 
      }
        
        

        
        
                
        
      drms_setkey_float(rs, "X0_MP", imageloc[i].x);
      drms_setkey_float(rs, "Y0_MP", imageloc[i].y);
      drms_setkey_float(rs, "INST_ROT", imageloc[i].instrot);
      //drms_setkey_float(rs, "INST_ROT", 180.0); //force this for now
      drms_setkey_float(rs, "IMSCL_MP", imageloc[i].imscale);
      drms_setkey_string(rs, "MPO_REC", imageloc[i].mpo_rec);

      int camera = drms_getkey_int(rs0, "CAMERA", &rstatus);
      if(rstatus) {
        printk("Can't do drms_getkey_int() for fsn %u\n", fsnx);
        noimage[i] = 1;
        goto TEMPSKIP;
        //return(1);
      }
      if ( 0 == drms_setkey_time(rs, "T_REC", tobs[i])) {
        int status, allstat = 0;
        double tr_step;
        long long tr_index;
        TIME t_rec, tr_epoch;
        tr_index = drms_getkey_longlong(rs, "T_REC_index", &status);
        allstat += status;
        tr_step =  drms_getkey_double(rs, "T_REC_step", &status);
        allstat += status;
        tr_epoch = drms_getkey_time(rs, "T_REC_epoch", &status);
        allstat += status;
        if (0 == allstat) {
          t_rec = tr_epoch + tr_index*tr_step;
          drms_setkey_time(rs, "T_REC", t_rec);
        }
      }
      //now get all the crs_table values according to lev0 IICRSID
      int iicrsid = drms_getkey_int(rs0, "IICRSID", &rstatus);
      if(rstatus) {
        printk("Can't do drms_getkey_int() for IICRSID for fsn %u\n", fsnx);
        printk("No CRS_TABLE keywords for lev1 fsn %u\n", fsnx);
      }
      else {
        sprintf(open_dsname, "%s[%d]", DSCRSNAME, iicrsid);
        printk("Open: %s\n", open_dsname);
        crsset = drms_open_records(drms_env, open_dsname, &rstatus);
        if(!crsset || (crsset->n == 0) || rstatus) {
          printk("Can't do drms_open_records(%s)\n", open_dsname);
        }
        else {
          crsrec = crsset->records[0];
          char *crsstr = drms_getkey_string(crsrec, "CRS_DESC", &status);
          drms_setkey_string(rs, "CRS_DESC", crsstr);
          crsstr = drms_getkey_string(crsrec, "CRS_TYPE", &status);
          drms_setkey_string(rs, "CRS_TYPE", crsstr);
          int crsint = drms_getkey_int(crsrec, "CRS_NREG", &status);
          drms_setkey_int(rs, "CRS_NREG", crsint);
          crsint = drms_getkey_int(crsrec, "TSR1", &status);
          drms_setkey_int(rs, "TSR1", crsint);
          crsint = drms_getkey_int(crsrec, "TER1", &status);
          drms_setkey_int(rs, "TER1", crsint);
          crsint = drms_getkey_int(crsrec, "TSC1", &status);
          drms_setkey_int(rs, "TSC1", crsint);
          crsint = drms_getkey_int(crsrec, "TEC1", &status);
          drms_setkey_int(rs, "TEC1", crsint);
          crsint = drms_getkey_int(crsrec, "TSR2", &status);
          drms_setkey_int(rs, "TSR2", crsint);
          crsint = drms_getkey_int(crsrec, "TER2", &status);
          drms_setkey_int(rs, "TER2", crsint);
          crsint = drms_getkey_int(crsrec, "TSC2", &status);
          drms_setkey_int(rs, "TSC2", crsint);
          crsint = drms_getkey_int(crsrec, "TEC2", &status);
          drms_setkey_int(rs, "TEC2", crsint);
          crsint = drms_getkey_int(crsrec, "TSR3", &status);
          drms_setkey_int(rs, "TSR3", crsint);
          crsint = drms_getkey_int(crsrec, "TER3", &status);
          drms_setkey_int(rs, "TER3", crsint);
          crsint = drms_getkey_int(crsrec, "TSC3", &status);
          drms_setkey_int(rs, "TSC3", crsint);
          crsint = drms_getkey_int(crsrec, "TEC3", &status);
          drms_setkey_int(rs, "TEC3", crsint);
          crsint = drms_getkey_int(crsrec, "TSR4", &status);
          drms_setkey_int(rs, "TSR4", crsint);
          crsint = drms_getkey_int(crsrec, "TER4", &status);
          drms_setkey_int(rs, "TER4", crsint);
          crsint = drms_getkey_int(crsrec, "TSC4", &status);
          drms_setkey_int(rs, "TSC4", crsint);
          crsint = drms_getkey_int(crsrec, "TEC4", &status);
          drms_setkey_int(rs, "TEC4", crsint);
          crsint = drms_getkey_int(crsrec, "TSR5", &status);
          drms_setkey_int(rs, "TSR5", crsint);
          crsint = drms_getkey_int(crsrec, "TER5", &status);
          drms_setkey_int(rs, "TER5", crsint);
          crsint = drms_getkey_int(crsrec, "TSC5", &status);
          drms_setkey_int(rs, "TSC5", crsint);
          crsint = drms_getkey_int(crsrec, "TEC5", &status);
          drms_setkey_int(rs, "TEC5", crsint);
          crsint = drms_getkey_int(crsrec, "TSR6", &status);
          drms_setkey_int(rs, "TSR6", crsint);
          crsint = drms_getkey_int(crsrec, "TER6", &status);
          drms_setkey_int(rs, "TER6", crsint);
          crsint = drms_getkey_int(crsrec, "TSC6", &status);
          drms_setkey_int(rs, "TSC6", crsint);
          crsint = drms_getkey_int(crsrec, "TEC6", &status);
          drms_setkey_int(rs, "TEC6", crsint);
          crsint = drms_getkey_int(crsrec, "TSR7", &status);
          drms_setkey_int(rs, "TSR7", crsint);
          crsint = drms_getkey_int(crsrec, "TER7", &status);
          drms_setkey_int(rs, "TER7", crsint);
          crsint = drms_getkey_int(crsrec, "TSC7", &status);
          drms_setkey_int(rs, "TSC7", crsint);
          crsint = drms_getkey_int(crsrec, "TEC7", &status);
          drms_setkey_int(rs, "TEC7", crsint);
          crsint = drms_getkey_int(crsrec, "TSR8", &status);
          drms_setkey_int(rs, "TSR8", crsint);
          crsint = drms_getkey_int(crsrec, "TER8", &status);
          drms_setkey_int(rs, "TER8", crsint);
          crsint = drms_getkey_int(crsrec, "TSC8", &status);
          drms_setkey_int(rs, "TSC8", crsint);
          crsint = drms_getkey_int(crsrec, "TEC8", &status);
          drms_setkey_int(rs, "TEC8", crsint);

          crsint = drms_getkey_int(crsrec, "WIN_FLIP", &status);
          l0l1->winflip = crsint;
          drms_setkey_int(rs, "WIN_FLIP", crsint); 
          printk("Close: %s\n", open_dsname);
          drms_close_records(crsset, DRMS_FREE_RECORD);
        }
      }

/****NEW from Rock 10Sep2013************************************************/
      if(!quicklook) {
/*  ITF1CCD1  CCD1_FUV1_OPERATING                                    */
/*  ITF2CCD2  CCD2_FUV2_OPERATING                                    */
/*  ITNUCCD3  CCD3_NUV_OPERATING                                     */
/*  ITSJCCD4  CCD4_SJI_OPERATING                                     */
/*  BT06CBPX  CEB_ON_THE_POSX_AXIS                                   */
/*  BT07CBNX  CEB_ON_THE_NEGX_AXIS                                   */
/*  BT15IEB   ELECTRONICS_BOX                                        */
/*  IT08GTWM  GUIDE_TELESCOPE_WEDGE_HC_MOTOR_TS08                    */
/*  IT14SPPX  SPECTROGRAPH_OPTICS_PACKAGE_HOZ5_CONTROL_TS14 (pos_x)  */
/*  IT16SPNX  SPECTROGRAPH_OPTICS_PACKAGE_HOZ7_CONTROL_TS16 (neg_x)  */

        float itf1ccd1, itf2ccd2, itnuccd3, itsjccd4;
        float bt06cbpx, bt07cbnx, bt15ieb, it08gtwm;
        float it14sppx, it16spnx;

        char *dstemp = "iris.temperatures_60s";
        if(fabs(tobs[i] - t_obs0) > 300.0) {
          char *selstr = "select max(date_obs) from ";
          char *whrstr = "where date_obs <= ";
          int nr;
          if(rs_t) {
            drms_close_records(rs_t, DRMS_FREE_RECORD);
            rs_t = NULL;
          }
          sprintf(open_dsname, "%s[? date_obs=(%s %s %s %f) ?]",
                  dstemp, selstr, dstemp, whrstr, tobs[i]);
          printf(" %s\n", open_dsname);
          rt = NULL;
          rs_t = drms_open_records(drms_env, open_dsname, &rstatus);
          if(rstatus) printk("Can not open temperature series.\n");
          else {
            nr = rs_t->n;
            if(nr != 1) printk("%d records != 1.\n", nr);
            rt = rs_t->records[0];
          }
          t_obs0 = tobs[i];
        }
        if (rt) {
          int st;
          sprintf(temprec, "%s[:#%lld]", dstemp, rt->recnum);
          if(dstatus = drms_setkey_string(rs, "TEMP_REC", temprec)) {
            printk("**ERROR on drms_setkey_string() for %s\n", temprec);
          }
          itf1ccd1 = drms_getkey_float(rt, "ITF1CCD1", &st);
          itf2ccd2 = drms_getkey_float(rt, "ITF2CCD2", &st);
          itnuccd3 = drms_getkey_float(rt, "ITNUCCD3", &st);
          itsjccd4 = drms_getkey_float(rt, "ITSJCCD4", &st);
          bt06cbpx = drms_getkey_float(rt, "BT06CBPX", &st);
          bt07cbnx = drms_getkey_float(rt, "BT07CBNX", &st);
          bt15ieb  = drms_getkey_float(rt, "BT15IEB",  &st);
          it08gtwm = drms_getkey_float(rt, "IT08GTWM", &st);
          it14sppx = drms_getkey_float(rt, "IT14SPPX", &st);
          it16spnx = drms_getkey_float(rt, "IT16SPNX", &st);

          drms_setkey_float(rs, "ITF1CCD1", itf1ccd1);
          drms_setkey_float(rs, "ITF2CCD2", itf2ccd2);
          drms_setkey_float(rs, "ITNUCCD3", itnuccd3);
          drms_setkey_float(rs, "ITSJCCD4", itsjccd4);
          drms_setkey_float(rs, "BT06CBPX", bt06cbpx);
          drms_setkey_float(rs, "BT07CBNX", bt07cbnx);
          drms_setkey_float(rs, "BT15IEB",  bt15ieb);
          drms_setkey_float(rs, "IT08GTWM", it08gtwm);
          drms_setkey_float(rs, "IT14SPPX", it14sppx);
          drms_setkey_float(rs, "IT16SPNX", it16spnx);
        }
      }
/****END NEW from Rock 10Sep2013************************************************/

      //For quicklook, find closest iris.timeline record & set scroll 26Nov2013
      if(quicklook) {
          DRMS_RecordSet_t *rset;
          DRMS_Record_t *rt;
          int st;
          TIME t_obs, roll_start;
          float roll_deg;
          t_obs = drms_getkey_time(rs0, "t_obs", &st);
        if (st) {
           t_obs = DRMS_MISSING_TIME;
        }
        else {
          sprintf(open_dsname, "iris.timeline_roll[? roll_start > %f and roll_start <= %f ?]", t_obs-86401, t_obs+1);
          printf("%s\n", open_dsname);
          rset = drms_open_records(drms_env, open_dsname, &st);

          if (!rset || !rset->n || st) {
              printk("Error in drms_open_records(%s); setting scroll to zero\n", open_dsname);
              scroll = 0.0;

          } else {
              // Pick last record
              //rt = rset->records[0];
              rt = rset->records[(rset->n)-1];
              sprintf(pointrec, "iris.timeline_roll[:#%lld]", rt->recnum);
              if(dstatus = drms_setkey_string(rs, "POINTREC", pointrec)) {
                printk("**ERROR on drms_setkey_string() for %s\n", pointrec);
              }
              roll_start = drms_getkey_time(rt, "ROLL_START", &st);
              roll_deg = drms_getkey_float(rt, "DEGREES", &st);

              drms_setkey_float(rs, "SAT_ROT", roll_deg);

              scroll = roll_deg;
          }

          if (rset)
              drms_close_records(rset, DRMS_FREE_RECORD);
        }
      }

      //char *wavstr = drms_getkey_string(rs0, "WAVE_STR", &rstatus);
      char *imgpath = drms_getkey_string(rs0, "IMG_PATH", &rstatus);
      if(rstatus) {
        printk("Can't do drms_getkey_string() for IMG_PATH\n");
        return(1);
      }
      if(quicklook) {
        sprintf(open_dsname, "%s[? t_start=(select max(t_start) from %s where t_start <= %10.5f and t_stop > %10.5f and IMG_PATH='%s') and IMG_PATH='%s' ?]",
      		dsffname, dsffname, tobs[i], tobs[i], imgpath, imgpath);
      }
      else {
        sprintf(open_dsname, "%s[? t_start <= %10.5f and t_stop > %10.5f and IMG_PATH='%s' ?]",
      	dsffname, tobs[i], tobs[i], imgpath);
      }
      //printf("!!TEMP Flat field query: %s\n", open_dsname); //!!TEMP
      //printk("!!TEMP Flat field query: %s\n", open_dsname); //!!TEMP
      rsetff = drms_open_records(drms_env, open_dsname, &rstatus); //open FF 
      if(!rsetff || (rsetff->n == 0) || rstatus) {
        printk("Can't do drms_open_records(%s)\n", open_dsname);
        flatmiss[i] = 1; noimage[i] = 1;
        goto TEMPSKIP;
        return(1);		//new 2/22/2011
      }
      fcnt = rsetff->n;
      if(fcnt > 1) {
        printk("More than one FF found for %s?\n", open_dsname);
        printk("Use last one of %d found\n", fcnt); //!!TEMP
        //return(1);		//!!TBD
      }
      //rsff = rsetff->records[0];
      rsff = rsetff->records[fcnt-1];
      recnumff = rsff->recnum;
      sprintf(flatrec, "%s[:#%lld]", dsffname, recnumff);
      if(dstatus = drms_setkey_string(rs, "FLAT_REC", flatrec )) {
        printk("**ERROR on drms_setkey_string() for %s\n", flatrec);
      }
      drms_record_directory(rsff, path, 1);
      if(!*path) {
        printk("***ERROR: No path to segment for %s\n", open_dsname);
        //goto TEMPSKIP;	//!!!TEMP until have good flatfield
        return(1);
      }
      //printf("\npath to FF = %s\n", path);	//!!TEMP
      segmentff = drms_segment_lookup(rsff, "flatfield");
      Arrayff = drms_segment_read(segmentff, DRMS_TYPE_FLOAT, &rstatus);
      if(!Arrayff) {
        printk("Can't do drms_segment_read() for Flat Field status=%d\n", 
			rstatus);
        return(1);
      }
      l0l1->adataff = (float *)Arrayff->data; //!!TBD free at end

/**************************No DARK for IRIS?**********************************/
      //darkseg = drms_segment_lookup(rsff, "DARK");
      darkseg = drms_segment_lookup(rsff, "dark");
      ArrayDark = drms_segment_read(darkseg, DRMS_TYPE_FLOAT, &rstatus);
      if(!ArrayDark) {
        printk("Can't do drms_segment_read() for dark. status=%d\n", rstatus);
        return(1);
      }
      l0l1->adatadark = (float *)ArrayDark->data; //free at end

      badseg = drms_segment_lookup(rsff, "bad_pixel");
      ArrayBad = drms_segment_read(badseg, DRMS_TYPE_INT, &rstatus);
      nbad = drms_array_size(ArrayBad)/sizeof(int);
      if(!ArrayBad) {
        printk("Can't do drms_segment_read() for bad_pixel. status=%d\n", 
			rstatus);
        return(1);
      }
      l0l1->adatabad = (int *)ArrayBad->data; //free at end
      //if(!(badoutpixseg = drms_segment_lookup(rs, "bad_pixel"))) {
      //    printk("No drms_segment_lookup(rs, bad_pixel) for lev1\n");
      //    return(1);
      //}
      l0l1->rs1 = rs;
      l0l1->rsff = rsff;
      l0l1->recnum1 = rs->recnum;  
      l0l1->darkflag = 0;
      l0l1->sumx = drms_getkey_short(rs0, "SUMSPTRL", &rstatus);
      l0l1->sumy = drms_getkey_short(rs0, "SUMSPAT", &rstatus);
      hshiexp = drms_getkey_int(rs, "HSHIEXP", &rstatus);
      hcamid = drms_getkey_int(rs, "HCAMID", &rstatus);
        float sumdc=0.0;
        int idc, numdc=0;
        int aimgshce = drms_getkey_int(rs, "AIMGSHCE", &rstatus);
        //if(aimgshce == 0) l0l1->darkflag = 1;
        short iifrmtyp = drms_getkey_short(rs0, "IIFRMTYP", &rstatus);
        if((iifrmtyp == 2) || (iifrmtyp == 5)) l0l1->darkflag = 1;
        if(rstatus = do_flat_iris(l0l1)) {
          printk("***ERROR in do_flat_iris() status=%d\n", rstatus);
          printf("***ERROR in do_flat_iris() status=%d\n", rstatus);
          flatmiss[i] = 1; noimage[i] = 1;
          //return(1);		//!!TBD what to do?
          goto FLATERR;
        }

        //sprintf(tmpname, "/tmp/data_lev1.%u", fsnx);
        //fwt = fopen(tmpname, "w");
        //int wsize = 4096 * 4096;
        //fwrite(l0l1->adata1, sizeof(float), wsize, fwt);
        //fclose(fwt);
        drms_setkey_float(rs, "OSCNMEAN", l0l1->oscnmean);
        drms_setkey_float(rs, "OSCNRMS", l0l1->oscnrms);
        drms_setkey_int(rs, "DATAMIN", l0l1->datamin);
        drms_setkey_int(rs, "DATAMAX", l0l1->datamax);
        drms_setkey_int(rs, "DATAMEDN", l0l1->datamedn);
        drms_setkey_float(rs, "DATAMEAN", l0l1->datamean);
        drms_setkey_float(rs, "DATARMS", l0l1->data_rms);
        drms_setkey_float(rs, "DATASKEW", l0l1->dataskew);
        drms_setkey_float(rs, "DATAKURT", l0l1->datakurt);
        drms_setkey_int(rs, "DATAVALS", l0l1->datavals);
        drms_setkey_int(rs, "MISSVALS", l0l1->missvals);
        missvals = (uint32_t)l0l1->missvals;
        totvals = (uint32_t)l0l1->datavals + missvals;
        drms_setkey_int(rs, "TOTVALS", (int)totvals);
        percentd = (float)((100.0 * (float)l0l1->datavals)/(float)totvals);
        drms_setkey_float(rs, "PERCENTD", percentd);
        if(missvals > 0) missflg[i] = missflg[i] | Q_1_MISS0;
        if(missvals > (uint32_t)(totvals * 0.01))
           missflg[i] = missflg[i] | Q_1_MISS1;
        if(missvals > (uint32_t)(totvals * 0.05)) 
           missflg[i] = missflg[i] | Q_1_MISS2;
        if(missvals > (uint32_t)(totvals * 0.25)) 
           missflg[i] = missflg[i] | Q_1_MISS3;
        if(l0l1->datavals == 0) 
           missflg[i] = missflg[i] | Q_MISSALL; //high bit, no data
        //if(nspikes) {
        //  if (spikeseg = drms_segment_lookup(rs,"spikes") ) {
        //    nbytes = nspikes*sizeof(int);
        //    axes[0] = nspikes;
        //    axes[1] = 3;
        //    spikedata = (int *) malloc (3*nspikes*sizeof(int));
        //    ArraySpike = drms_array_create(DRMS_TYPE_INT, 2, axes,
        //               (void *) spikedata, &status);
        //    memcpy((void *)spikedata, (void *)spikelocs, nbytes);
        //    memcpy((void *)(spikedata+nspikes), (void *)oldvalues, nbytes);
        //    memcpy((void *)(spikedata+2*nspikes), (void *)newvalues, nbytes);
        //    status = drms_segment_write(spikeseg, ArraySpike, 0);
        //    drms_free_array(ArraySpike);
        //  } // else { printf("spikes segment not found\n"); }
        //}

/**************************No DARK for IRIS?**********************************/

FLATERR:
      drms_close_records(rsetff, DRMS_FREE_RECORD);
      free(ArrayDark->data);
      free(Arrayff->data);
      free(ArrayBad->data);

TEMPSKIP:
  x0_lf = DRMS_MISSING_DOUBLE;
  y0_lf = DRMS_MISSING_DOUBLE;
  rsun_lf = DRMS_MISSING_DOUBLE;
  //Don't send calmode or darks to limb_fit()
  //calmode: HCFTID=17
  //dark: HSHIEXP=0 and HCAMID=0 or 1
  int skiplimb = 0;
  hcftid = drms_getkey_int(rs, "HCFTID", &rstatus);
  if(hcftid == CAL_HCFTID) {	//don't call limb_fit()
    printk("Cal mode image fsn=%u\n", fsnx);
    skiplimb = 1;
  }
  else {
    if(hshiexp == 0) {
      if(hcamid == 0 || hcamid == 1) {
        printk("Dark image fsn=%u\n", fsnx);
        skiplimb = 1;
      }
    }
  }
  lstatus = 1;				//default bad limb_fit
  //goto WCSEND;		//skip all WCS stuff for IRIS

WCSEND:

  //if(hmiaiaflg) {                       //aia
  //  int wl = drms_getkey_int(rs, "WAVELNTH", &rstatus);
  //}

  do_quallev1(rs0, rs, i, fsnx);
  //ftmp = StopTimer(2); //!!!TEMP
  //printf( "\nTime sec in inside loop for fsn=%u : %f\n\n", fsnx, ftmp );

IRISSKIP:
      dstatus = drms_segment_writewithkeys(segment, segArray, 0);
      //dstatus = drms_segment_write(segment, segArray, 0);
      if (dstatus) {
        printk("ERROR: drms_segment_write error=%d for fsn=%u\n",
dstatus,fsnx);
        noimage[i] = 1;
      }
      recnum1 = rs->recnum;
      printk("*1 %u %u\n", recnum1, fsnx);
      free(Array0->data);
//      if (rs_t) drms_close_records(rs_t, DRMS_FREE_RECORD);
      if (rs_resp) drms_close_records(rs_resp, DRMS_FREE_RECORD);
      rs_resp = NULL;
  }				//END do for all the sorted lev0 records

  drms_close_records(rset0, DRMS_FREE_RECORD);   //close lev0 records
  drms_close_records(rset1, DRMS_INSERT_RECORD); //close lev1 records
  return(0);
}

int compare_rptr(const void *a, const void *b)
{
  TIME t1, t2;
  int rstatus;
  DRMS_Record_t *x=(DRMS_Record_t *)a, *y=(DRMS_Record_t *)b;

  t1 = drms_getkey_time(x, "t_obs", &rstatus); 
  if(rstatus) t1 = DRMS_MISSING_TIME;	//treat error as missing t_obs
  t2 = drms_getkey_time(y, "t_obs", &rstatus); 
  if(rstatus) t2 = DRMS_MISSING_TIME;
  if(t1 < t2) return(-1);
  if(t1 > t2) return(1);
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
  sprintf(string, "build_lev1_iris started as pid=%d ppid=%d user=%s\n", 
		getpid(), getppid(), username);
  strcat(idstr, string);
  printk("%s", idstr);
  printf("%s", idstr);
  if(restartflg) printk("-r ");
  sprintf(argmode, "mode=%s", mode);
  sprintf(argdsin, "dsin=%s", dsin);
  sprintf(argdsout, "dsout=%s", dsout);
  if(modeflg) {
    sprintf(argbx, "brec=%lld", brec);
    sprintf(argex, "erec=%lld", erec);
  }
  else {
    sprintf(argbx, "bfsn=%lld", bfsn);
    sprintf(argex, "efsn=%lld", efsn);
  }
  sprintf(argquick, "quicklook=%d", quicklook);
  sprintf(arglogfile, "logfile=%s", logname);
  printk("%s %s %s %s %s %s %s\n", 
	argmode, argdsin, argdsout, argbx, argex, argquick, arglogfile);
  printf("%s %s %s %s %s %s %s\n", 
	argmode, argdsin, argdsout, argbx, argex, argquick, arglogfile);
  if(!restartflg) {
    //printk("tlmseriesname=%s\nlev0seriesname=%s\n", 
    //		tlmseriesname, lev0seriesname);
  }
  sprintf(bld_vers, "%s", jsoc_version);
  sprintf(idstr,  "ps -ef | grep %s", LEV1VIEWERNAME);
  fin = popen(idstr, "r");
  while(fgets(string, sizeof string, fin)) {  //get ps line
    if(!(strstr(string, "perl"))) continue;
    sscanf(string, "%s %d", idstr, &tpid); /* get user name & process id */
    sprintf(lfile, "%s/build_lev1_iris_restart_%d.touch", LEV1LOG_BASEDIR, tpid);
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
  long long numofrecs, frec, lrec;
  int numrec, numofchunks, i;
  char line[80];

  if (nice_intro())
    return (0);
  if(!(username = (char *)getenv("USER"))) username = "nouser"; 
  mode = cmdparams_get_str(&cmdparams, "mode", NULL);
  if(strcmp(mode, "recnum") && strcmp(mode, "fsn")) {
    printf("Error: mode= must be given as 'recnum' or 'fsn'\n");
    return(0);
  }
  if(!strcmp(mode, "recnum")) modeflg = 1;
  dsin = cmdparams_get_str(&cmdparams, "dsin", NULL);
  dsout = cmdparams_get_str(&cmdparams, "dsout", NULL);
  brec = cmdparams_get_int(&cmdparams, "brec", NULL);
  erec = cmdparams_get_int(&cmdparams, "erec", NULL);
  bfsn = cmdparams_get_int(&cmdparams, "bfsn", NULL);
  efsn = cmdparams_get_int(&cmdparams, "efsn", NULL);
  quicklook = cmdparams_get_int(&cmdparams, "quicklook", NULL);
  //quicklook = 1; //!!TEMP for test
  if(modeflg) {		//recnum mode
    if(brec == 0 || erec == 0) {
      fprintf(stderr, "brec and erec must be given for recnum mode. 0 not allowed\n");
      return(0);
    }
    if(brec > erec) {
      fprintf(stderr, "brec must be <= erec\n");
      return(0);
    }
    bnumx = brec;
    enumx = erec;
  }
  else {		//fsn mode
    if(bfsn == 0 || efsn == 0) {
      fprintf(stderr, "bfsn and efsn must be given for fsn mode. 0 not allowed\n");
      return(0);
    }
    if(bfsn > efsn) {
      fprintf(stderr, "bfsn must be <= efsn\n");
      return(0);
    }
    bnumx = bfsn;
    enumx = efsn;
  }
  logfile = cmdparams_get_str(&cmdparams, "logfile", NULL);
  if (strcmp(dsin, NOTSPECIFIED) == 0) {
    dsin = LEV0SERIESNAMEIRIS;
  }
  if (strcmp(dsout, NOTSPECIFIED) == 0) {
    dsout = LEV1SERIESNAMEIRIS;
  }
  sprintf(dsffname, "%s", DSFFNAME);

  if (strcmp(logfile, NOTSPECIFIED) == 0) {
    sprintf(logname, H1LOGFILE, gettimetag());
  }
  else {
    sprintf(logname, "%s", logfile);
  }
  if(restartflg || LOGTEST) {
    if((h1logfp=fopen(logname, "a")) == NULL)
      fprintf(stderr, "**Can't open for append the log file %s\n", logname);
  }
  else {
    if((h1logfp=fopen(logname, "w")) == NULL)
      fprintf(stderr, "**Can't open the log file %s\n", logname);
  }
  setup();
  numofrecs = (enumx - bnumx) + 1;
  numrec = NUMRECLEV1;	   //# of records to do at a time
  numofchunks = numofrecs/numrec;
  if((numofrecs % numrec) != 0) numofchunks++; //extra loop for partial chunk
  lrec = bnumx-1;
  for(i = 0; i < numofchunks; i++) {
    frec = lrec+1; lrec = (frec + numrec)-1;
    if(lrec > enumx) lrec=enumx;
    if(do_ingest(frec, lrec)) {  //do a chunk to get files from the lev0
      printf("build_lev1_iris abort\nSee log: %s\n", logname); 
      send_mail("build_lev1_iris abort\nSee log: %s\n", logname); 
      return(0);
    }
    if(modeflg) {               //only do for recnum mode
      if(stat(STOP_FILE, &stbuf) == 0) {
        printf("Stop file %s seen. Exit build_lev1_iris.\n", STOP_FILE);
        break;
      }
    }
  }
  printf("build_lev1_iris done last fsn=%u\n", fsnx); 
  return(0);
}
