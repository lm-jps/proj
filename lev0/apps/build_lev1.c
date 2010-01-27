/*-----------------------------------------------------------------------------
 * cvs/JSOC/proj/lev1/apps/build_lev1.c
 *-----------------------------------------------------------------------------
 *
 * This is a module that runs with DRMS and processes lev0
 * filtergrams to lev1.
 * It is scheduled by build_lev1_mgr either by qsub to the cluster
 * or by fork to run on the local machine. It is called with 
 * 17 (NUMRECLEV1 defined in lev0lev1.h) lev0 images at a time.
 */ 
//!!!TBD - in development.

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
#include "imgdecode.h"
#include "lev0lev1.h"


//default in and out data series
//#define LEV0SERIESNAMEHMI "hmi.lev0e"
//#define LEV0SERIESNAMEAIA "aia.lev0e"
#define LEV0SERIESNAMEHMI "su_production.lev0f_hmi"
#define LEV0SERIESNAMEAIA "su_production.lev0f_aia"
//#define LEV1SERIESNAMEHMI "su_production.hmi_lev1e"	//temp test case
#define LEV1SERIESNAMEHMI "hmi.lev1"
#define LEV1SERIESNAMEAIA "su_production.aia_lev1e"	//temp test case
//#define DSFFNAME "su_richard.flatfield"			//temp test case
#define DSFFNAME "su_production.hmi_flatfield"		//temp test case

#define LEV1LOG_BASEDIR "/usr/local/logs/lev1"
#define H1LOGFILE "/usr/local/logs/lev1/build_lev1.%s.log"
#define NUMTIMERS 8		//number of seperate timers avail
#define NOTSPECIFIED "***NOTSPECIFIED***"
#define LOGTEST 0

int compare_rptr(const void *a, const void *b);

/******************************************************
extern int decode_next_hk_vcdu(unsigned short *tbuf, CCSDS_Packet_t **hk, unsigned int *Fsn);
extern int write_hk_to_drms();
extern void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg);
extern int set_HMI_mech_values(DRMS_Record_t *rec);
******************************************************/

static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);

// List of default parameter values. 
ModuleArgs_t module_args[] = { 
  {ARG_STRING, "instru", NOTSPECIFIED, "instrument. either hmi or aia"},
  {ARG_STRING, "dsin", NOTSPECIFIED, "dataset of lev0 filtergrams"},
  {ARG_STRING, "dsout", NOTSPECIFIED, "dataset of lev1 output"},
  {ARG_STRING, "logfile", NOTSPECIFIED, "optional log file name. Will create one if not given"},
  {ARG_INTS, "brec", "-1", "first lev0 rec# to process. -1=error must be given by build_lev1_mgr"},
  {ARG_INTS, "erec", "-1", "last lev0 rec# to process. -1=error must be given by build_lev1_mgr"},
  {ARG_INTS, "quicklook", "1", "1=quick look, 0 = definitive mode"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "h", "0", "help flag"},
  {ARG_FLAG, "r", "0", "restart flag"},
  {ARG_END}
};

CmdParams_t cmdparams;
// Module name presented to DRMS. 
char *module_name = "build_lev1";

FILE *h1logfp;		// fp for h1 ouput log for this run 
static IMG Image0, Image1;
static IMG *Img0 = &Image0;
static IMG *Img1 = &Image1;
static LEV0LEV1 lev0lev1;
static LEV0LEV1 *l0l1 = &lev0lev1;
//static CCSDS_Packet_t *Hk;
static DRMS_Record_t *rs;
static DRMS_Record_t *rs0, *rs1, *rsff, *rsbad_pix;
static DRMS_Record_t *rptr;
static DRMS_Segment_t *segment;
static DRMS_Segment_t *segmentff;
static DRMS_Segment_t *darkseg;
static DRMS_Segment_t *badseg;
static DRMS_Array_t *segArray;
static DRMS_RecordSet_t *rset0, *rset1, *rsetff;
static DRMS_Array_t *Array0;
static DRMS_Array_t *Arrayff;
static DRMS_Array_t *ArrayDark;
static DRMS_Array_t *ArrayBad;
static TIME sdo_epoch;
static char datestr[32];
static char open_dsname[256];
static char path[DRMS_MAXPATHLEN], bad_pix_path[DRMS_MAXPATHLEN];
static char rs1_path[DRMS_MAXPATHLEN];
static struct timeval first[NUMTIMERS], second[NUMTIMERS];
static char *orbseries = "sdo.fds_orbit_vectors";

IORBIT_Info_t *IOinfo = NULL;
IORBIT_Info_t IOdata;
LIBASTRO_Error_t IOstatus;
unsigned int fsnarray[NUMRECLEV1];
unsigned int fsnx = 0;
//short data1[MAXPIXELS];
int data1[MAXPIXELS];
double tgttimes[NUMRECLEV1];

long long brec, erec;           //begin and end lev0 rec# to do
int verbose;
int hmiaiaflg = 0;		//0=hmi, 1=aia
int imagecnt = 0;		// num of images since last commit 
int restartflg = 0;		// set when ingest_lev0 is called for a restart
int abortflg = 0;
int sigalrmflg = 0;             // set on signal so prog will know 
int ignoresigalrmflg = 0;       // set after a close_image()
int quicklook;
char logname[128];
char argdsin[128], argdsout[128], arglogfile[128], arginstru[80];
char argbrec[80], argerec[80], argquick[80];
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
    printf ("Usage:\nbuild_lev1 [-vhr] "
	"instru=<hmi|aia> dsin=<lev0> dsout=<lev1> brec=<rec#>\n"
	"                 erec=<rec#> quicklook=<0|1> [logfile=<file>]\n"
	"  -h: help - show this message then exit\n"
	"  -v: verbose\n"
	"  -r: restart. only used when we restart our selves periodically\n"
	"instru= instrument. must be 'hmi' or 'aia'\n"
	"dsin= data set name of lev0 input\n"
	"      default hmi=hmi.lev0e   aia=aia.lev0e\n"
	"dsout= data set name of lev1 output\n"
	"      default hmi=su_production.hmi_lev1e   aia=su_production.aia_lev1e\n"
	"brec= first lev0 rec# to process. -1=error must be given by build_lev1_mgr\n"
	"erec= last lev0 rec# to process. -1=error must be given by build_lev1_mgr\n"
	"quicklook= 1 = quicklook mode, 0 = definitive mode\n"
	"logfile= optional log file name. If not given uses:\n"
        "         /usr/local/logs/lev1/build_lev1.<time_stamp>.log\n");
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
  printk("**Exit build_lev1 w/ status = %d\n", stat);
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

//!!TBD Keh-Cheng
int sc_pointing()
{
  return(0);
}

#include "do_flat.c"
#include "get_pointing_info.c"
#include "get_image_location.c"

int do_ingest(long long bbrec, long long eerec)
{
  PTINFO *ptinfo = NULL;
  PTINFO ptdata;
  Image_Location *p_imageloc;
  Image_Location imageloc[NUMRECLEV1];
  TIME t_obs0;
  TIME tobs[NUMRECLEV1];
  float percentd;
  int rstatus, dstatus, ncnt, fcnt, i, j, qualint, nobs, totvals;
  long long recnum0, recnum1, recnumff;
  char recrange[128], lev0name[128], flatrec[128];

  sprintf(recrange, ":#%lld-#%lld", bbrec, eerec);
  sprintf(open_dsname, "%s[%s]", dsin, recrange);
  printk("open_dsname = %s\n", open_dsname);
  printk("#levnum recnum fsn\n");

    rset0 = drms_open_records(drms_env, open_dsname, &rstatus); //open lev0
    if(!rset0 || (rset0->n == 0) || rstatus) {
      printk("Can't do drms_open_records(%s)\n", open_dsname);
      return(1);
    }
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
    if(rstatus = get_pointing_info(drms_env, tobs, ncnt, &ptinfo)) {
      printk("**ERROR: get_pointing_info() status = %d  fsn tobs ASD:\n", 
		rstatus);
      for(j=0; j < ncnt; j++) {	 //!!TEMP debuf stuff
        printk("%u %10.5f ", fsnarray[j], tobs[j]);
        ptdata = ptinfo[j];
        printk("%s\n", ptdata.asd_rec);
      }
    }
    if ((IOstatus = iorbit_getinfo(drms_env,
                       orbseries,
                       NULL,
                       IORBIT_Alg_Quadratic,
                       tobs,
                       ncnt,
                       kIORBIT_CacheAction_DontCache,
                       &IOinfo)) != kLIBASTRO_Success)
    {
      if(IOstatus == kLIBASTRO_InsufficientData) {
        printk("***ERROR in iorbit_getinfo: kLIBASTRO_InsufficientData\n");
      }
      else { 
       printk("***ERROR in iorbit_getinfo() status=%d\n", IOstatus);
      }
    }
    rset1 = drms_create_records(drms_env, ncnt, dsout, DRMS_PERMANENT,&dstatus);
    if(dstatus) {
      printk("**ERROR: Can't create records for %s\n", dsout);
      return(1);
    }
    //Now fill in info for call to Carl's get_image_location()
    for(i=0; i < ncnt; i++) {
      rs0 = &rptr[i];
      imageloc[i].tobs = tobs[i];
      imageloc[i].camera = drms_getkey_int(rs0, "CAMERA", &rstatus);
      if(rstatus) {
        printk("ERROR: in drms_getkey_int(CAMERA) fsn=%u\n", fsnarray[i]);
      }
      imageloc[i].wavelength = drms_getkey_int(rs0, "WAVELNTH", &rstatus);
      if(rstatus) {
        printk("ERROR: in drms_getkey_int(WAVELNTH) fsn=%u\n", fsnarray[i]);
      }
      snprintf(imageloc[i].telescope, 10, "%s", 
		drms_getkey_string(rs0, "TELESCOP", &rstatus));
      if(rstatus) {
        printk("ERROR: in drms_getkey_string(TELESCOP) fsn=%u\n", fsnarray[i]);
      }
    }
    p_imageloc = imageloc;
    rstatus = get_image_location(drms_env, ncnt, &p_imageloc);


    for(i=0; i < ncnt; i++) { 	//do for all the sorted lev0 records
      rs0 = &rptr[i];
      recnum0 = rs0->recnum;
      fsnx = fsnarray[i]; 
      sprintf(lev0name, "%s[%u]", dsin, fsnx);
      //printf("rec# for %d = %lld fsn=%u\n", i, recnum0, fsnx); //!!!TEMP
      segment = drms_segment_lookupnum(rs0, 0);
      Array0 = drms_segment_read(segment, DRMS_TYPE_SHORT, &rstatus);
      if(!Array0) {
        printk("Can't do drms_segment_read() %s status=%d\n", 
			open_dsname, rstatus);
        return(1);
      }
      //short *adata = (short *)Array0->data;
      //memcpy(Img0->dat, adata, 2*MAXPIXELS);
      //drms_free_array(Array0); //must free from drms_segment_read()
      l0l1->adata0 = (short *)Array0->data; //!!TBD free at end
      l0l1->adata1 = (int *)&data1;
      l0l1->rs0 = rs0;
      l0l1->recnum0 = recnum0;
      l0l1->fsn = fsnx;
      l0l1->himgcfid = drms_getkey_int(rs0, "HIMGCFID", &rstatus);
      //!!TEMP force a good himgcfid for our 'junk' data
      //if(l0l1->himgcfid < 0 || l0l1->himgcfid >= MAXHIMGCFGS) {
      //  l0l1->himgcfid = 104;
      //}
      if(rstatus) {
        printk("Can't do drms_getkey_int(HIMGCFID) for fsn %u\n", fsnx);
        return(1); 
      }

      sprintf(open_dsname, "%s[%u]", dsout, fsnx);
      rs = rset1->records[i]; 
      drms_record_directory(rs, rs1_path, 0);
      if(!*rs1_path) {
        printk("***ERROR: No path to segment for %s\n", open_dsname);
        //goto TEMPSKIP;	//!!!TEMP until have good flatfield
        return(1);
      }
      //printf("\npath to lev1 = %s\n", rs1_path);	//!!TEMP
      dstatus = drms_setkey_int(rs, "FSN", fsnx);
      dstatus = drms_setkey_string(rs, "LEV0SERIES", lev0name);
      if(!(segment = drms_segment_lookup(rs, "image_lev1"))) {
        printk("No drms_segment_lookup(rs, image_lev1) for %s\n", open_dsname);
        return(1);
      }
      segArray = drms_array_create(DRMS_TYPE_INT,
                                       segment->info->naxis,
                                       segment->axis,
                                       &data1,
                                       &dstatus);
      //transer the lev0 keywords
      rstatus = drms_copykeys(rs, rs0, 0, kDRMS_KeyClass_Explicit);
      if(rstatus != DRMS_SUCCESS) {
        printk("Error %d in drms_copykeys() for fsn %u\n", fsnx);
        return(1);
      }
      qualint = drms_getkey_int(rs0, "QUALITY", &rstatus);
      drms_setkey_int(rs, "QUALLEV0", qualint);

      drms_setkey_time(rs, "T_OBS", tobs[i]);
      printk("t_obs for lev0 = %10.5f\n", tobs[i]);  //!!TEMP
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
          drms_setkey_string(rs, "ACS_SAFE", ptdata.acs_safe);
          drms_setkey_string(rs, "ASD_REC", ptdata.asd_rec);
          drms_setkey_string(rs, "ACS_CGT", ptdata.acs_cgt);
          //note: acs_cgt not yet in jsd
      }
      if(IOinfo) {
        IOdata = IOinfo[i];
           drms_setkey_double(rs, "HCIEC_X", IOdata.hciX);
           drms_setkey_double(rs, "HCIEC_Y", IOdata.hciY);
           drms_setkey_double(rs, "HCIEC_Z", IOdata.hciZ);
           drms_setkey_double(rs, "GCIEC_X", IOdata.gciX);
           drms_setkey_double(rs, "GCIEC_Y", IOdata.gciY);
           drms_setkey_double(rs, "GCIEC_Z", IOdata.gciZ);
           drms_setkey_float(rs, "DSUN_OBS", (float)IOdata.dsun_obs);
           drms_setkey_double(rs, "OBS_VR", IOdata.obs_vr);
           drms_setkey_double(rs, "OBS_VW", IOdata.obs_vw);
           drms_setkey_double(rs, "OBS_VN", IOdata.obs_vn);
           drms_setkey_float(rs, "CRLN_OBS", (float)IOdata.crln_obs);
           drms_setkey_float(rs, "CRLT_OBS", (float)IOdata.crlt_obs);
           drms_setkey_int(rs, "CAR_ROT", (int)IOdata.car_rot);
      }
      drms_setkey_float(rs, "X0_MP", imageloc[i].x);
      drms_setkey_float(rs, "Y0_MP", imageloc[i].y);
      drms_setkey_float(rs, "INST_ROT", imageloc[i].instrot);
      drms_setkey_float(rs, "IMSCL_MP", imageloc[i].imscale);
      drms_setkey_string(rs, "MPO_REC", imageloc[i].mpo_rec);

      int camera = drms_getkey_int(rs0, "CAMERA", &rstatus);
      if(rstatus) {
        printk("Can't do drms_getkey_int() for fsn %u\n", fsnx);
        return(1);
      }
      //Now figure out what flat field to use
      if(drms_ismissing_time(tobs[i])) {
        printk("DRMS_MISSING_TIME for fsn=%u. Continue...\n", fsnx);
        goto TEMPSKIP;
      }
      //sprintf(open_dsname, "%s[? t_start <= %10.5f and t_stop > %10.5f and CAMERA=%d and FLATFIELD_FINAL=1 ?]", 
      //sprintf(open_dsname, "%s[? t_start <= %10.5f and t_stop > %10.5f and CAMERA=%d ?]", 
      //	DSFFNAME, tobs[i], tobs[i], camera);
      if(quicklook) {
        sprintf(open_dsname, "%s[? t_start=(select max(t_start) from %s where t_start <= %10.5f and t_stop > %10.5f and CAMERA=%d) and CAMERA=%d ?]",
      		DSFFNAME, DSFFNAME, tobs[i], tobs[i], camera, camera);
      }
      else {
        sprintf(open_dsname, "%s[? t_start <= %10.5f and t_stop > %10.5f and CAMERA=%d and flatfield_version=(select max(flatfield_version) from %s where t_start <= %10.5f and t_stop > %10.5f and CAMERA=%d) ?]",
      	DSFFNAME, tobs[i], tobs[i], camera, DSFFNAME, tobs[i], tobs[i], camera);
      }
      //printf("!!TEMP Flat field query: %s\n", open_dsname); //!!TEMP
      rsetff = drms_open_records(drms_env, open_dsname, &rstatus); //open FF 
      if(!rsetff || (rsetff->n == 0) || rstatus) {
        printk("Can't do drms_open_records(%s)\n", open_dsname);
        return(1);
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
      sprintf(flatrec, "%s[:#%lld]", DSFFNAME, recnumff);
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

      darkseg = drms_segment_lookup(rsff, "DARK");
      ArrayDark = drms_segment_read(darkseg, DRMS_TYPE_FLOAT, &rstatus);
      if(!ArrayDark) {
        printk("Can't do drms_segment_read() for DARK. status=%d\n", rstatus);
        return(1);
      }
      l0l1->adatadark = (float *)ArrayDark->data; //!!TBD free at end

      badseg = drms_segment_lookup(rsff, "BAD_PIXEL");
      ArrayBad = drms_segment_read(badseg, DRMS_TYPE_INT, &rstatus);
      if(!ArrayBad) {
        printk("Can't do drms_segment_read() for BAD_PIXEL. status=%d\n", 
			rstatus);
        return(1);
      }
      l0l1->adatabad = (int *)ArrayBad->data; //!!TBD free at end

      rsbad_pix = drms_link_follow(rsff, "perm_bad_pixel", &rstatus);
      if(rstatus) {
        printk("Can't do drms_link_follow(rsff, \"perm_bad_pixel\"). status=%d\n", rstatus);
        return(1);
      }
      drms_record_directory(rsbad_pix, bad_pix_path, 1);
      if(!*bad_pix_path) {
        printk("***ERROR: No path to segment for BAD_PIXEL\n");
        return(1);
      }
      //printf("\npath to BAD_PIXEL = %s\n", bad_pix_path); //!!TEMP
      //cp the bad pixel list to the lev1 output
      char cmd[128];
      sprintf(cmd, "cp %s/bad_pixel_list.fits %s", bad_pix_path, rs1_path);
      printk("%s\n", cmd);
      system(cmd);

      l0l1->rs1 = rs;
      l0l1->rsff = rsff;
      l0l1->recnum1 = rs->recnum;  

      if(rstatus = do_flat(l0l1)) {
        printk("***ERROR in do_flat() status=%d\n", rstatus);
        printf("***ERROR in do_flat() status=%d\n", rstatus);
        //return(1);		//!!TBD what to do?
      }
      else {
        drms_setkey_short(rs, "DATAMIN", (short)l0l1->datamin);
        drms_setkey_short(rs, "DATAMAX", (short)l0l1->datamax);
        drms_setkey_short(rs, "DATAMEDN", (short)l0l1->datamedn);
        drms_setkey_float(rs, "DATAMEAN", l0l1->datamean);
        drms_setkey_float(rs, "DATARMS", l0l1->data_rms);
        drms_setkey_float(rs, "DATASKEW", l0l1->dataskew);
        drms_setkey_float(rs, "DATAKURT", l0l1->datakurt);
        drms_setkey_int(rs, "DATAVALS", l0l1->datavals);
        drms_setkey_int(rs, "MISSVALS", l0l1->missvals);
        totvals = l0l1->datavals + l0l1->missvals;
        drms_setkey_int(rs, "TOTVALS", totvals);
        percentd = (float)((100.0 * (float)l0l1->datavals)/(float)totvals);
        drms_setkey_float(rs, "PERCENTD", percentd);
      }
      drms_close_records(rsetff, DRMS_FREE_RECORD);
      free(ArrayDark->data);
      free(Arrayff->data);
      free(ArrayBad->data);

TEMPSKIP:

      dstatus = drms_segment_write(segment, segArray, 0);
      if (dstatus) {
        printk("ERROR: drms_segment_write error=%d for fsn=%u\n", dstatus,fsnx);
      }
      recnum1 = rs->recnum;
      printk("*1 %u %u\n", recnum1, fsnx);
      free(Array0->data);
  
     //!!TBD add here a call for get_limb_fit() from Richard. The call
    // is TBD. It works on the lev1 image. It returns values for
    // float rsun_lf, float x0_lf, float y0_lf

    }

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
  getcwd(cwdbuf, 126);
  sprintf(idstr, "Cwd: %s\nCall: ", cwdbuf);
  sprintf(string, "build_lev1 started as pid=%d ppid=%d user=%s\n", 
		getpid(), getppid(), username);
  strcat(idstr, string);
  printk("%s", idstr);
  printf("%s", idstr);
  if(restartflg) printk("-r ");
  sprintf(arginstru, "instru=%s", instru);
  sprintf(argdsin, "dsin=%s", dsin);
  sprintf(argdsout, "dsout=%s", dsout);
  sprintf(argbrec, "brec=%lld", brec);
  sprintf(argerec, "erec=%lld", erec);
  sprintf(argquick, "quicklook=%d", quicklook);
  sprintf(arglogfile, "logfile=%s", logname);
  printk("%s %s %s %s %s %s %s\n", 
	arginstru, argdsin, argdsout, argbrec, argerec, argquick, arglogfile);
  printf("%s %s %s %s %s %s %s\n", 
	arginstru, argdsin, argdsout, argbrec, argerec, argquick, arglogfile);
  if(!restartflg) {
    //printk("tlmseriesname=%s\nlev0seriesname=%s\n", 
    //		tlmseriesname, lev0seriesname);
  }
  sprintf(idstr,  "ps -ef | grep %s", LEV1VIEWERNAME);
  fin = popen(idstr, "r");
  while(fgets(string, sizeof string, fin)) {  //get ps line
    if(!(strstr(string, "perl"))) continue;
    sscanf(string, "%s %d", idstr, &tpid); /* get user name & process id */
    sprintf(lfile, "%s/build_lev1_restart_%d.touch", LEV1LOG_BASEDIR, tpid);
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
  quicklook = cmdparams_get_int(&cmdparams, "quicklook", NULL);
  if(brec == -1 || erec == -1) {
    fprintf(stderr, "brec and/or erec must be given. -1 not allowed\n");
    return(0);
  }
  if(brec > erec) {
    fprintf(stderr, "brec must be <= erec\n");
    return(0);
  }
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
  if(restartflg || LOGTEST) {
    //sleep(30);		//make sure old ingest_lev0 is done
    if((h1logfp=fopen(logname, "a")) == NULL)
      fprintf(stderr, "**Can't open for append the log file %s\n", logname);
  }
  else {
    if((h1logfp=fopen(logname, "w")) == NULL)
      fprintf(stderr, "**Can't open the log file %s\n", logname);
  }
  setup();
  numofrecs = (erec - brec) + 1;
  numrec = NUMRECLEV1;	   //# of records to do at a time
  numofchunks = numofrecs/numrec;
  if((numofrecs % numrec) != 0) numofchunks++; //extra loop for partial chunk
  lrec = brec-1;
  for(i = 0; i < numofchunks; i++) {
    frec = lrec+1; lrec = (frec + numrec)-1;
    if(lrec > erec) lrec=erec;
    if(do_ingest(frec, lrec)) {  //do a chunk to get files from the lev0
      //printk("**ERROR: Some error after open of %s\n", open_dsname);
      //printf("**ERROR: Some error after open of %s\n", open_dsname);
      printf("build_lev1 abort\n"); //!!TEMP
      printf("See log: %s\n", logname);
      return(0);
    }
  }
  printf("build_lev1 done last fsn=%u\n", fsnx); //!!TEMP
  return(0);
}

//!!TEMP position of keywords set up for limb_figure code when it
//goes in
void tmp() {

typedef struct {
  float rsun_lf;
//  double rsun_obs;
//  float imscl_mp;
  float x0_lf;
  float y0_lf;
//  float x0_mp;
//  float y0_mp;
//  float sat_y0;
//  float sat_z0;
//  float inst_rot;
//  float sat_rot;
} LIMB_SOMETHING;

LIMB_SOMETHING limb;

  if(limb->rsun_lf == DRMS_MISSING_FLOAT) {
    drms_setkey_float(rs, "CDELT1", imageloc[i].imscale);
    drms_setkey_float(rs, "CDELT2", imageloc[i].imscale);
    drms_setkey_float(rs, "R_SUN", (float)RSUN_OBS/imageloc[i].imscale); //TBD get RSUN_OBS from iorbit
  }

  
}
