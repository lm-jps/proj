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
return(sdo_epoch + (TIME)sdo_s + (TIME)(sdo_ss)/65536.0;
}

static void sprint_time_ISO (char *tstring, TIME t)
{
sprint_time(tstring,t,"UTC",0);
tstring[4] = tstring[7] = '-';
tstring[10] = 'T';
tstring[19] = '\0';
} 

void HMI_compute_exposure_times(DRMS_Record_t *rec, HK_Keyword_t *isp, int flg)
{
int hobitsec = HK_getkey_int(isp, "HOBITSEC");   /* HMI_OBT_IMG_TIME_SHM_SEC */
int hobitss  = HK_getkey_int(isp, "HOBITSS");    /* HMI_OBT_IMG_TIME_SHM_SS  */
int hshmiclb = HK_getkey_int(isp, "HSHMICLB");   /* HMI_SHM_IMG_CLOSE_BOTTOM */
int hshmiclm = HK_getkey_int(isp, "HSHMICLM");   /* HMI_SHM_IMG_CLOSE_MIDDLE */
int hshmiclt = HK_getkey_int(isp, "HSHMICLT");   /* HMI_SHM_IMG_CLOSE_TOP    */
int hshmiopb = HK_getkey_int(isp, "HSHMIOPB");   /* HMI_SHM_IMG_OPEN_BOTTOM  */
int hshmiopm = HK_getkey_int(isp, "HSHMIOPM");   /* HMI_SHM_IMG_OPEN_MIDDLE  */
int hshmiopt = HK_getkey_int(isp, "HSHMIOPT");   /* HMI_SHM_IMG_OPEN_TOP     */
double sheb = (hshmiclb - hshmiopb)*1.0e-6;
double shem = (hshmiclm - hshmiopm)*1.0e-6;
double shet = (hshmiclt - hshmiopt)*1.0e-6;
double offset = ( hshmiclb + hshmiopb + hshmiclm + hshmiopm + hshmiclt + hshmiopt )/6.0;
TIME t_obs = SDO_to_DRMS_time(hobitsec,hobitss) + offset;
double exptime = (sheb + shem + shet)/3.0;
double expsdev = sqrt((sheb*sheb + shem*shem + shet*shet)/3.0 - exptime*exptime);

char date_obs[100];
sprint_time_ISO(date_obs, t_obs - exptime/2.0);
TIME MJD_epoch = -3727641600.000; /* 1858.11.17_00:00:00_UT  */
TIME date__obs = t_obs - exptime/2.0;
TIME mjd = date__obs - MJD_epoch;
double  mjd_day = floor(mjd / 86400.0);
double  mjd_time = mjd - 86400.0 * mjd_day;


//drms_setkey_TIME(rec, "T_OBS", t_obs);
drms_setkey_double(rec, "T_OBS", t_obs);
drms_setkey_double(rec, "EXPTIME", exptime);
drms_setkey_double(rec, "EXPSDEV", expsdev);
drms_setkey_string(rec, "DATE__OBS", date_obs);
drms_setkey_double(rec, "MJD", mjd_day);
drms_setkey_double(rec, "TIME", mjd_time);

if(flg == 0) {			// HMI
  int camid = HK_getkey_int(isp, "HCAMID");
  camid = camid & 0x01;
  camid++;
  drms_setkey_int(rec, "CAMERA", camid);
  if(camid == 1) drms_setkey_string(rec, "INSTRUME", "HMI_SIDE1");
  else drms_setkey_string(rec, "INSTRUME", "HMI_FRONT2");
}
}

