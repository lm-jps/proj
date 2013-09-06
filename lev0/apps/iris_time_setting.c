#include <drms_keyword.h>
#include "packets.h"
#include "printk.h"

extern int INVALtime;
//extern unsigned int cropid;
extern unsigned int fsnx;

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
  //double isppktim = HK_getkey_double(isp, "PACKET_TIME");
  //double isppktim = drms_getkey_time(rec, "PACKET_TIME", &status);

  int use_pktim;

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
  use_pktim = 0;
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
      use_pktim = 1;
      break;
    case 3:
      drms_setkey_string(rec, "IMG_TYPE", "LED"); 
      break;
    case 4:
      drms_setkey_string(rec, "IMG_TYPE", "LTC"); 
      break;
    case 5:
      drms_setkey_string(rec, "IMG_TYPE", "SPAT"); 
      use_pktim = 1;
      break;
    case 6:
      drms_setkey_string(rec, "IMG_TYPE", "VPAT"); 
      use_pktim = 1;
      break;
    default:
      drms_setkey_string(rec, "IMG_TYPE", "UNKNOWN"); 
      break;
  }

  drms_setkey_double(rec, "EXPTIME", exptime);
  drms_setkey_double(rec, "EXPSDEV", expdev);
  if ((iimgcfd1 == DRMS_MISSING_SHORT) || (iimgcfd4 == DRMS_MISSING_SHORT)) {
    int_time = DRMS_MISSING_FLOAT; /* no info to calculate integration time */
  } else {
    int_time = (iimgcfd4 - iimgcfd1)*0.0078125;
  }
  drms_setkey_float(rec, "INT_TIME", int_time);
  drms_setkey_string(rec, "INSTRUME", iris_instru[isqisysn]);
  //if (use_pktim) t_obs = isppktim;
  //else t_obs = SDO_to_DRMS_time(iimgots, iimgotss) + expoff/2.0;
  if(use_pktim) {
    iimgots = HK_getkey_int(isp, "ITCS56");
    iimgotss = HK_getkey_int(isp, "ITCSS56"); 
    iimgotss = iimgotss >> 16; 
  }
  //t_obs = SDO_to_DRMS_time(iimgots, iimgotss) + expoff/2.0;
  t_obs = SDO_to_DRMS_time(iimgots, iimgotss) - exptime/2.0;
  drms_setkey_double(rec, "T_OBS", t_obs);
  
  TIME date__obs = t_obs - exptime/2.0;
  drms_setkey_double(rec, "DATE__OBS", date__obs);
  drms_setkey_int(rec, "CAMERA", iris_camera[isqisysn]);
  int aecmode = HK_getkey_int(isp, "AECMODE");

  //int cropid = drms_getkey_int(rec, "CROPID", &status); //Not in there yet
  int iicrsid = HK_getkey_int(isp, "IICRSID"); //use instead of cropid
  int ifwpos = HK_getkey_int(isp, "IFWPOS");
  //int iifuvfdb = HK_getkey_int(isp, "IIFUVFDB");
  //int iinuvfdb = HK_getkey_int(isp, "IINUVFDB");
  //int iisjifdb = HK_getkey_int(isp, "IISJIFDB");
printk("fsnx = %d\n", fsnx);
printk("iicrsid = %d\n", iicrsid);
printk("ifwpos = %d\n", ifwpos);
//printk("iifuvfdb = %d\n", iifuvfdb);
//printk("iinuvfdb = %d\n", iinuvfdb);
//printk("iisjifdb = %d\n", iisjifdb);
//For debug. Print out entire key list
//HK_Keyword_t *isptmp;
//isptmp = isp;
//do {
//  printk("%s raw_value=%u\n", isptmp->name, isptmp->raw_value);
//  isptmp = isptmp->next;
//} while(isptmp->next);

  switch (isqisysn) {
    case 0: 
	    //drms_setkey_short(rec, "IIFDBID", iifuvfdb);
            drms_setkey_string(rec, "IMG_PATH", "FUV");
            break;
    case 1: 
	    //drms_setkey_short(rec, "IIFDBID", iinuvfdb);
            if (iicrsid < 4) drms_setkey_string(rec, "IMG_PATH", "NUV-SJI");
            else drms_setkey_string(rec, "IMG_PATH", "NUV");
            break;
    case 2: 
            //drms_setkey_short(rec, "IIFDBID", iisjifdb);
            if (iicrsid < 4) { drms_setkey_string(rec, "IMG_PATH", "NUV-SJI"); }
            else {
	      switch (ifwpos) {
              case 1: 
              case 2: 
	            drms_setkey_string(rec, "IMG_PATH", "SJI_5000W"); 
	            break;
              case 31: 
              case 32: 
	            drms_setkey_string(rec, "IMG_PATH", "SJI_1330"); 
	            break;
              case 61: 
              case 62: 
	            drms_setkey_string(rec, "IMG_PATH", "SJI_2796"); 
	            break;
              case 91: 
              case 92: 
	            drms_setkey_string(rec, "IMG_PATH", "SJI_1400"); 
	            break;
              case 121: 
              case 122: 
	            drms_setkey_string(rec, "IMG_PATH", "SJI_2832"); 
	            break;
              case 151: 
              case 152: 
	            drms_setkey_string(rec, "IMG_PATH", "SJI_1600W"); 
	            break;
              default:
	            drms_setkey_string(rec, "IMG_PATH", "SJI_UNKNOWN"); 
                    break;
              }
            }
            break;
  }
}
