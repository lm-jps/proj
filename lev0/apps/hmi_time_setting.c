/* block of code to do the HMI shutter time calculations
   and set some keywords.  This block sets:
    T_OBS
    EXPTIME
    EXPSDEV
    DATE__OBS

 Call: HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg)

 rec is lev0 record just before ISP keywords are added.
 isp is lev0 ISP linked-list from the ISP HK telem decode.

*/

#include <drms_keyword.h>
#include "packets.h"
#include "printk.h"

static int HK_getkey_int(HK_Keyword_t *isp, char *key)
{
while (isp)
  {
  if (strcmp(key,isp->fitsname)==0)
     return((int)isp->raw_value);
  isp = isp->next;
  }
return(DRMS_MISSING_INT);
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

static void sprint_time_ISO (char *tstring, TIME t)
{
sprint_time(tstring,t,"UTC",0);
tstring[4] = tstring[7] = '-';
tstring[10] = 'T';
tstring[19] = '\0';
} 

/* Two ways to calculate counter rollovers. Take your pick.
static int nrollct(int cmdt, int clt)
{
  if (cmdt <  50332) return 0;
  if (cmdt <  83886) return ((clt > 0x800000) ? 0 : 1);
  if (cmdt < 117441) return 1;
  if (cmdt < 150995) return ((clt > 0x800000) ? 1 : 2);
  if (cmdt < 184550) return 2;
  if (cmdt < 218104) return ((clt > 0x800000) ? 2 : 3);
  if (cmdt < 251659) return 3;
  return ((clt > 0x800000) ? 3 : 4);
}
*/
static int nrollct(int cmdt, int clt)
{
  if (cmdt < 50332) return 0;
  int nr = (cmdt - 16777)/67109;
  int rem = cmdt%67109;
  if ((rem < 16777) || (rem > 50332)) return nr + (clt > 0x800000 ? 0 : 1);
  return nr;
}

void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg)
{
static float waveltab[10] = { 33.5, 13.1,
                              21.1, 19.3,
                              160.0, 170.0, 450.0, 17.1,
                               30.4, 9.4 };
static char date_obs[100]; 
//static char datestr[100];
TIME t_obs,MJD_epoch = -3727641600.000; /* 1858.11.17_00:00:00_UT  */
double exptime, expsdev;

if (flg) {	/* AIA */
  char *aia_instru[] = { "AIA_1", "AIA_2", "AIA_3", "AIA_4" };
  int aimgshce = HK_getkey_int(isp, "AIMGSHCE"); /* AIA_IMG_SH_CMDED_EXPOSURE */
  int aimgots  = HK_getkey_int(isp, "AIMGOTS");  /* AIA_IMG_OBT_TIME_SH_SEC  */
  int aimgotss = HK_getkey_int(isp, "AIMGOTSS"); /* AIA_IMG_OBT_TIME_SH_SS   */
  int aimshcbc = HK_getkey_int(isp, "AIMSHCBC"); /* AIA_IMG_SH_CLOSE_BOT_CENTR*/
  aimshcbc += 0x1000000*nrollct(aimgshce, aimshcbc);
  int aimshcbe = HK_getkey_int(isp, "AIMSHCBE"); /* AIA_IMG_SH_CLOSE_BOT_EDGE */
  aimshcbe += 0x1000000*nrollct(aimgshce, aimshcbe);
  int aimshctc = HK_getkey_int(isp, "AIMSHCTC"); /* AIA_IMG_SH_CLOSE_TOP_CENTR*/
  aimshctc += 0x1000000*nrollct(aimgshce, aimshctc);
  int aimshcte = HK_getkey_int(isp, "AIMSHCTE"); /* AIA_IMG_SH_CLOSE_TOP_EDGE */
  aimshcte += 0x1000000*nrollct(aimgshce, aimshcte);
  int aimshobc = HK_getkey_int(isp, "AIMSHOBC"); /* AIA_IMG_SH_OPEN_BOT_CENTR*/
  int aimshobe = HK_getkey_int(isp, "AIMSHOBE"); /* AIA_IMG_SH_OPEN_BOT_EDGE */
  int aimshotc = HK_getkey_int(isp, "AIMSHOTC"); /* AIA_IMG_SH_OPEN_TOP_CENTR*/
  int aimshote = HK_getkey_int(isp, "AIMSHOTE"); /* AIA_IMG_SH_OPEN_TOP_EDGE */
  double shebc = (aimshcbc - aimshobc)*4.0e-6;
  double shebe = (aimshcbe - aimshobe)*4.0e-6;
  double shetc = (aimshctc - aimshotc)*4.0e-6;
  double shete = (aimshcte - aimshote)*4.0e-6;
  double offset = ( aimshobc + aimshcbc + aimshobe + aimshcbe +
                    aimshotc + aimshctc + aimshote + aimshcte )/8.0e6;
  if(aimgshce == 0) {			//use pkt time for t_obs
    int axsec = HK_getkey_int(isp, "ATCSISP");
    int axssec = HK_getkey_int(isp, "ATCSSISP");
    if((axsec == DRMS_MISSING_INT) || (axssec == DRMS_MISSING_INT)) {
      axsec = HK_getkey_int(isp, "ATCS211");
      axssec = HK_getkey_int(isp, "ATCSS211");
    }
    if((axsec == DRMS_MISSING_INT) || (axssec == DRMS_MISSING_INT)) {
      axsec = HK_getkey_int(isp, "ATCS239");
      axssec = HK_getkey_int(isp, "ATCSS239");
    }
    if((axsec != DRMS_MISSING_INT) && (axssec != DRMS_MISSING_INT)) {
      axssec = axssec >> 16;
      t_obs = SDO_to_DRMS_time(axsec, axssec);
    }
  }
  else {
    t_obs = SDO_to_DRMS_time(aimgots, aimgotss) + offset;
  }
  exptime = (shebc + shebe + shetc + shete)/4.0;
  expsdev = sqrt((shebc*shebc + shebe*shebe + shetc*shetc +
                 shete*shete)/4.0 - exptime*exptime);
  int asqhdr = HK_getkey_int(isp, "ASQHDR");
  if (DRMS_MISSING_INT == asqhdr) asqhdr = HK_getkey_int(isp, "A827A");
  int ahtelid = (asqhdr >> 30) & 0x3;
  int ahtlfsn = asqhdr & 0x3FFFFFF;
  int asqtnum = HK_getkey_int(isp, "ASQTNUM");
  drms_setkey_string(rec, "INSTRUME", aia_instru[ahtelid]);
//  CAMERA set elsewhere, as is FSN
//  drms_setkey_int(rec, "CAMERA", ahtelid + 1);
  int wavel = HK_getkey_int(isp, "AIAWVLEN");    /* AIA_IMG_WAVELENGTH */
  if (wavel == DRMS_MISSING_INT) wavel = HK_getkey_int(isp, "A856C");
  if (wavel != DRMS_MISSING_INT) {
    if (wavel < 10) drms_setkey_float(rec, "WAVELNTH", waveltab[wavel]);
    else drms_setkey_float(rec, "WAVELNTH", (float) wavel);
    /* delete this line and above else line after testing */
    /* if wavel > 9, WAVELNTH = DRMS_MISSING_FLOAT by default */
  }
  int aecmode = HK_getkey_int(isp, "AECMODE");   /* AIA_IMG_AEC_MODE         */
  if (aecmode == DRMS_MISSING_INT) aecmode = HK_getkey_int(isp, "A8328");
  if (aecmode == DRMS_MISSING_INT) aecmode = HK_getkey_int(isp, "A8215-00");
  if (aecmode != DRMS_MISSING_INT) {
    if (aecmode) drms_setkey_string(rec, "A8215_00", "On");
    else drms_setkey_string(rec, "A8215_00", "Off");
  }
} else {	/* HMI */
  int hshiexp = HK_getkey_int(isp, "HSHIEXP"); /* HMI_FSW_IMG_CMDED_EXPOSURE */
  int hobitsec = HK_getkey_int(isp, "HOBITSEC"); /* HMI_OBT_IMG_TIME_SHM_SEC */
  int hobitss  = HK_getkey_int(isp, "HOBITSS");  /* HMI_OBT_IMG_TIME_SHM_SS  */
  int hshmiclb = HK_getkey_int(isp, "HSHMICLB"); /* HMI_SHM_IMG_CLOSE_BOTTOM */
  int hshmiclm = HK_getkey_int(isp, "HSHMICLM"); /* HMI_SHM_IMG_CLOSE_MIDDLE */
  int hshmiclt = HK_getkey_int(isp, "HSHMICLT"); /* HMI_SHM_IMG_CLOSE_TOP    */
  int hshmiopb = HK_getkey_int(isp, "HSHMIOPB"); /* HMI_SHM_IMG_OPEN_BOTTOM  */
  int hshmiopm = HK_getkey_int(isp, "HSHMIOPM"); /* HMI_SHM_IMG_OPEN_MIDDLE  */
  int hshmiopt = HK_getkey_int(isp, "HSHMIOPT"); /* HMI_SHM_IMG_OPEN_TOP     */
  double sheb = (hshmiclb - hshmiopb)*1.0e-6;
  double shem = (hshmiclm - hshmiopm)*1.0e-6;
  double shet = (hshmiclt - hshmiopt)*1.0e-6;
  double offset = ( hshmiclb + hshmiopb + hshmiclm + hshmiopm + hshmiclt + hshmiopt )/6.0e6;
  
  if(hshiexp == 0) {			//use pkt time for t_obs
    int hxsec = HK_getkey_int(isp, "HTCSISP");
    int hxssec = HK_getkey_int(isp, "HTCSSISP");
    if((hxsec == DRMS_MISSING_INT) || (hxssec == DRMS_MISSING_INT)) {
      hxsec = HK_getkey_int(isp, "HTCS1BD");
      hxssec = HK_getkey_int(isp, "HTCSS1BD");
    }
    if((hxsec == DRMS_MISSING_INT) || (hxssec == DRMS_MISSING_INT)) {
      hxsec = HK_getkey_int(isp, "HTCS1DB");
      hxssec = HK_getkey_int(isp, "HTCSS1DB");
    }
    if((hxsec != DRMS_MISSING_INT) && (hxssec != DRMS_MISSING_INT)) {
      hxssec = hxssec >> 16;
      t_obs = SDO_to_DRMS_time(hxsec, hxssec);
    }
  }
  else {
    t_obs = SDO_to_DRMS_time(hobitsec,hobitss) + offset;
  }
  exptime = (sheb + shem + shet)/3.0;
  expsdev = sqrt((sheb*sheb + shem*shem + shet*shet)/3.0 - exptime*exptime);
  
  int camid = HK_getkey_int(isp, "HCAMID");
  camid = camid & 0x01;
  camid++;
  drms_setkey_int(rec, "CAMERA", camid);
  if(camid == 1) drms_setkey_string(rec, "INSTRUME", "HMI_SIDE1");
  else drms_setkey_string(rec, "INSTRUME", "HMI_FRONT2");
}
  
//sprint_time_ISO(date_obs, t_obs - exptime/2.0);
TIME date__obs = t_obs - exptime/2.0;
TIME mjd = date__obs - MJD_epoch;
double  mjd_day = floor(mjd / 86400.0);
double  mjd_time = mjd - 86400.0 * mjd_day;
//sprint_time(datestr, CURRENT_SYSTEM_TIME, "ISO", 0);
//sprint_time(date_obs, date__obs, "ISO", 0);
  
drms_setkey_double(rec, "T_OBS", t_obs);
drms_setkey_double(rec, "EXPTIME", exptime);
drms_setkey_float(rec, "EXPSDEV", expsdev);
//drms_setkey_string(rec, "DATE__OBS", date_obs);
drms_setkey_double(rec, "DATE__OBS", date__obs);
drms_setkey_double(rec, "MJD", mjd_day);
drms_setkey_double(rec, "TIME", mjd_time);
//NOTE: DATE now set in close_image() in ingest_lev0.c
//drms_setkey_string(rec, "DATE", datestr);
//drms_setkey_double(rec, "DATE", CURRENT_SYSTEM_TIME);
}

