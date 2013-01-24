#include <drms_keyword.h>
#include "packets.h"
#include "printk.h"

extern int INVALtime;

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
return(sdo_epoch + (TIME)sdo_s + (TIME)(sdo_ss & 0xFFFF)/65536.0);
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
  static char date_obs[100];
  TIME t_obs,MJD_epoch = -3727641600.000; /* 1858.11.17_00:00:00_UT  */
  double exptime=0.0, expdev=0.0, expoff=0.0;
  double she[6];
  float int_time;
  int frmtyp;
  int i, iimgots, iimgotss, iris_camera[3] = {1, 2, 2}, status;
  int ifuvot[6], inuvot[5], isjiot[5], ifuvct[6], inuvct[5], isjict[5];
  char *iris_instru[] = { "FUV", "NUV", "SJI" };
  char *fuvotkw[] = { "IFUVAOT", "IFUVBOT", "IFUVCOT", "IFUVDOT", "IFUVEOT",
                      "IFUVFOT" };
  char *nuvotkw[] = { "INUVAOT", "INUVBOT", "INUVCOT", "INUVDOT", "INUVEOT" };
  char *sjiotkw[] = { "ISJIAOT", "ISJIBOT", "ISJICOT", "ISJIDOT", "ISJIEOT" };
  char *fuvctkw[] = { "IFUVACT", "IFUVBCT", "IFUVCCT", "IFUVDCT", "IFUVECT",
                      "IFUVFCT" };
  char *nuvctkw[] = { "INUVACT", "INUVBCT", "INUVCCT", "INUVDCT", "INUVECT" };
  char *sjictkw[] = { "ISJIACT", "ISJIBCT", "ISJICCT", "ISJIDCT", "ISJIECT" };

  int isqisysn = HK_getkey_int(isp, "ISQISYSN");
  int iimgcfd1 = HK_getkey_int(isp, "IIMGCFD1");
  int iimgcfd2 = HK_getkey_int(isp, "IIMGCFD2");
  int iimgcfd3 = HK_getkey_int(isp, "IIMGCFD3");
  int iimgcfd4 = HK_getkey_int(isp, "IIMGCFD4");
  int iimgshce = HK_getkey_int(isp, "IIMGSHCE");

  for (i=0; i<6; i++) {
    ifuvot[i] = HK_getkey_int(isp, fuvotkw[i]);
    ifuvct[i] = HK_getkey_int(isp, fuvctkw[i]);
    if(i<5) {
      inuvot[i] = HK_getkey_int(isp, nuvotkw[i]);
      inuvct[i] = HK_getkey_int(isp, nuvctkw[i]);
      isjiot[i] = HK_getkey_int(isp, sjiotkw[i]);
      isjict[i] = HK_getkey_int(isp, sjictkw[i]);
    }
  }

  switch (isqisysn) {
    case 0:
      if (ifuvot[0] < ifuvot[1]) {
        while (ifuvct[0] < ifuvct[2]) ifuvct[0] += 16777216;/*remove rollover*/
        ifuvct[0] -= 16777216;             /* A close time is < C close time */
        while (ifuvct[1] < ifuvct[0]) ifuvct[1] += 16777216;
        while (ifuvct[3] < ifuvct[2]) ifuvct[3] += 16777216;
        while (ifuvct[4] < ifuvct[2]) ifuvct[4] += 16777216;
        while (ifuvct[5] < ifuvct[2]) ifuvct[5] += 16777216;
      } else {
        ifuvct[5] += 16777216;
        while (ifuvct[5] < ifuvct[2]) ifuvct[5] += 16777216;/*remove rollover*/
        ifuvct[5] -= 16777216;             /* A close time is < C close time */
        while (ifuvct[4] < ifuvct[5]) ifuvct[4] += 16777216;
        while (ifuvct[3] < ifuvct[4]) ifuvct[3] += 16777216;
        while (ifuvct[1] < ifuvct[2]) ifuvct[1] += 16777216;
        while (ifuvct[0] < ifuvct[2]) ifuvct[0] += 16777216;
      }
      for (i=0; i<6; i++) {
        expoff += (ifuvct[i] + ifuvot[i])*4.0e-6/6.0;
        she[i] = (ifuvct[i] - ifuvot[i])*4.0e-6;
        exptime += she[i];
        expdev += she[i]*she[i];
      }
      exptime /= 6.0;
      expdev = sqrt(expdev/6.0 - exptime*exptime);
      iimgots = HK_getkey_int(isp, "IIMGOTS1");
      iimgotss = HK_getkey_int(isp, "IMGOTSS1");
      break;
    case 1:
      if (inuvot[0] < inuvot[1]) {
        while (inuvct[0] < inuvct[2]) inuvct[0] += 65536; /* remove rollover */
        inuvct[0] -= 65536;                /* A close time is < C close time */
        while (inuvct[1] < inuvct[0]) inuvct[1] += 65536;
        while (inuvct[3] < inuvct[2]) inuvct[3] += 65536;
        while (inuvct[4] < inuvct[2]) inuvct[4] += 65536;
      } else {
        inuvct[4] += 65536;
        while (inuvct[4] < inuvct[2]) inuvct[4] += 65536; /* remove rollover */
        inuvct[4] -= 65536;                /* A close time is < C close time */
        while (inuvct[3] < inuvct[4]) inuvct[3] += 65536;
        while (inuvct[1] < inuvct[2]) inuvct[1] += 65536;
        while (inuvct[0] < inuvct[2]) inuvct[0] += 65536;
      }
      for (i=0; i<5; i++) {
        expoff += (inuvct[i] + inuvot[i])*4.0e-6/5.0;
        she[i] = (inuvct[i] - inuvot[i])*4.0e-6;
        exptime += she[i];
        expdev += she[i]*she[i];
      }
      exptime /= 5.0;
      expdev = sqrt(expdev/5.0 - exptime*exptime);
      iimgots = HK_getkey_int(isp, "IIMGOTS2");
      iimgotss = HK_getkey_int(isp, "IMGOTSS2");
      break;
    case 2:
      if (isjiot[0] < isjiot[1]) {
        while (isjict[0] < isjict[2]) isjict[0] += 65536; /* remove rollover */
        isjict[0] -= 65536;                /* A close time is < C close time */
        while (isjict[1] < isjict[0]) isjict[1] += 65536;
        while (isjict[3] < isjict[2]) isjict[3] += 65536;
        while (isjict[4] < isjict[2]) isjict[4] += 65536;
      } else {
        isjict[4] += 65536;
        while (isjict[4] < isjict[2]) isjict[4] += 65536; /* remove rollover */
        isjict[4] -= 65536;                /* A close time is < C close time */
        while (isjict[3] < isjict[4]) isjict[3] += 65536;
        while (isjict[1] < isjict[2]) isjict[1] += 65536;
        while (isjict[0] < isjict[2]) isjict[0] += 65536;
      }
      for (i=0; i<5; i++) {
        expoff += (isjict[i] + isjiot[i])*4.0e-6/5.0;
        she[i] = (isjict[i] - isjiot[i])*4.0e-6;
        exptime += she[i];
        expdev += she[i]*she[i];
      }
      exptime /= 5.0;
      expdev = sqrt(expdev/5.0 - exptime*exptime);
      iimgots = HK_getkey_int(isp, "IIMGOTS3");
      iimgotss = HK_getkey_int(isp, "IMGOTSS3");
      break;
  }
  if((iimgots == 0) && (iimgshce != 0)) INVALtime = 1;
  else INVALtime = 0;
  
  //if(iimgshce == 0) {
  //  drms_setkey_string(rec, "IMG_TYPE", "DARK");
  //} else {
  //  drms_setkey_string(rec, "IMG_TYPE", "LIGHT");
  //}
  frmtyp = HK_getkey_int(isp, "IIFRMTYP");
  switch (frmtyp) {
    case 0:
      drms_setkey_string(rec, "IMG_TYPE", "CURRENT"); 
      break;
    case 1:
      drms_setkey_string(rec, "IMG_TYPE", "LIGHT"); 
      break;
    case 2:
      drms_setkey_string(rec, "IMG_TYPE", "DARK"); 
      break;
    case 3:
      drms_setkey_string(rec, "IMG_TYPE", "LED"); 
      break;
    case 4:
      drms_setkey_string(rec, "IMG_TYPE", "LTC"); 
      break;
    case 5:
      drms_setkey_string(rec, "IMG_TYPE", "SPAT"); 
      break;
    case 6:
      drms_setkey_string(rec, "IMG_TYPE", "VPAT"); 
      break;
    default:
      drms_setkey_string(rec, "IMG_TYPE", "UNKNOWN"); 
      break;
  }

  drms_setkey_double(rec, "EXPTIME", exptime);
  drms_setkey_double(rec, "EXPSDEV", expdev);
  if ((iimgcfd3 == DRMS_MISSING_SHORT) || (iimgcfd4 == DRMS_MISSING_SHORT)) {
    int_time = DRMS_MISSING_FLOAT; /* no info to calculate integration time */
  } else {
    int_time = (iimgcfd4 - iimgcfd3)*0.0078125;
  }
  drms_setkey_float(rec, "INT_TIME", int_time);
  drms_setkey_string(rec, "INSTRUME", iris_instru[isqisysn]);
  t_obs = SDO_to_DRMS_time(iimgots, iimgotss) + expoff/2.0;
  drms_setkey_double(rec, "T_OBS", t_obs);
  
  TIME date__obs = t_obs - exptime/2.0;
  drms_setkey_double(rec, "DATE__OBS", date__obs);
  drms_setkey_int(rec, "CAMERA", iris_camera[isqisysn]);
  int aecmode = HK_getkey_int(isp, "AECMODE");

}
