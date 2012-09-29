#include "jsoc.h"
#include "jsoc_main.h"

char *module_name = "aia_fixCROTA2";

/*
 * This program updates CROTA2 based on results of the 2012 transit of Venus.
 * CROTA2 is corrected by adding CCDx_DELTA 
 * If INST_ROT is present it will be updated.
 * DATE is updated.
 * A HISTORY line is added.
 *
 * call with aia_fixCROTA2 ds=<recordset_query>
 *
 */

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

#define NOTSPECIFIED "Not Specified"
#define MP_SERIES "sdo.master_pointing"
#define GMP_MAX_MPO_REC_SIZE 200

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "ds", NOTSPECIFIED,  "data series to process."},
     {ARG_END}
};

static char mp_query[512], respq[512];

int DoIt(void)
{
  int status = DRMS_SUCCESS;
  DRMS_Record_t *rmp, *rresp;
  DRMS_RecordSet_t *inRS, *outRS, *rsmp=NULL, *rs_resp=NULL;;
  int irec, nrecs, nresp;
  const char *dsSeries = params_get_str(&cmdparams, "ds");
  TIME t_obs, mpt_end=0.0;

  inRS = drms_open_records(drms_env, dsSeries, &status);
  nrecs = inRS->n;
  if (status || nrecs == 0)
  {
    fprintf(stdout, " status=%d, no records found, skip this block\n",status);
    fflush(stdout);
    return (DRMS_SUCCESS);
  }

  outRS = drms_clone_records_nosums(inRS, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &status);
  nrecs = outRS->n;
  if (status || nrecs == 0)
    DIE("No records cloned");
  drms_close_records(inRS, DRMS_FREE_RECORD);

  for (irec=0; irec<nrecs; irec++)
  {
    DRMS_Record_t *rec = outRS->records[irec];
    float arot[10], ascl[10];
    float crota2, sat_rot;
    int i, quality, instrot_status, wl;
    char mpo_rec[GMP_MAX_MPO_REC_SIZE];
    char *rotkeys[] = { "A_094_INSTROT", "A_131_INSTROT", "A_171_INSTROT",
       "A_193_INSTROT", "A_211_INSTROT", "A_304_INSTROT", "A_335_INSTROT",
       "A_1600_INSTROT", "A_1700_INSTROT", "A_4500_INSTROT"};
    char *sclkeys[] = { "A_094_IMSCALE", "A_131_IMSCALE", "A_171_IMSCALE",
       "A_193_IMSCALE", "A_211_IMSCALE", "A_304_IMSCALE", "A_335_IMSCALE",
       "A_1600_IMSCALE", "A_1700_IMSCALE", "A_4500_IMSCALE"};
    char *dsresp = "aia.response";
    char *wavstr = drms_getkey_string(rec, "WAVE_STR", &status);

    quality = drms_getkey_int(rec, "QUALITY", &status);
    if (!status && quality < 0)
      continue;

    t_obs = drms_getkey_time(rec, "T_OBS", &status);
    sat_rot = drms_getkey_float(rec, "SAT_ROT", NULL);
    if (drms_ismissing_float(sat_rot)) sat_rot = 0.0;
    wl = drms_getkey_int(rec, "WAVELNTH", NULL);

    // The action is all between here and the end of the irec loop
    if (t_obs > mpt_end) {
      int iw, nrmp;
      sprintf(mp_query, "%s[? T_START <= %f and %f < T_STOP ?]",
              MP_SERIES, t_obs, t_obs);
      rsmp = drms_open_records(drms_env, mp_query, &status);
      nrmp = rsmp->n;
      if (status || nrmp != 1) DIE("No Master Pointing Record");
      rmp = rsmp->records[0];
      sprintf(mpo_rec, "%s[:#%lld]", MP_SERIES, rmp->recnum);
      mpt_end = drms_getkey_time(rmp, "T_STOP", &status);
      if (status) DIE("No T_STOP for Master Pointing Record");
      for (iw=0; iw<10; iw++) {
        arot[iw] = drms_getkey_float(rmp, rotkeys[iw], &instrot_status);
        if (instrot_status) DIE("Value for inst rotation not found");
        ascl[iw] = drms_getkey_float(rmp, sclkeys[iw], &status);
        if (status) DIE("Value for inst scale not found");
      }
      drms_close_records(rsmp, DRMS_FREE_RECORD);
    }
    switch (wl) {
      case   94: i = 0; break;
      case  131: i = 1; break;
      case  171: i = 2; break;
      case  193: i = 4; break;
      case  211: i = 4; break;
      case  304: i = 5; break;
      case  335: i = 6; break;
      case 1600: i = 7; break;
      case 1700: i = 8; break;
      case 4500: i = 9; break;
      default: DIE("Bad wavelength in AIA level 1 record");
    }
    drms_setkey_float(rec, "INST_ROT", arot[i]);
    drms_setkey_float(rec, "CDELT1", ascl[i]);
    drms_setkey_float(rec, "CDELT2", ascl[i]);
    drms_setkey_float(rec, "IMSCL_MP", ascl[i]);
    crota2 = sat_rot + arot[i];
    drms_setkey_float(rec, "CROTA2", crota2);
    drms_setkey_time(rec, "DATE", CURRENT_SYSTEM_TIME);
    drms_setkey_string(rec, "MPO_REC", mpo_rec);
    sprintf(respq, "%s[][%s][? t_start <= %10.5f and t_stop > %10.5f ?]",
            dsresp, wavstr, t_obs, t_obs);
    rresp = NULL;
    rs_resp = drms_open_records(drms_env, respq, &status);
    if (status) { DIE("Can not open aia.response series.\n");
    } else  {
      nresp = rs_resp->n;
      if (nresp>0) rresp = rs_resp->records[0];
    }
    if(rresp) {
      float eperdn = drms_getkey_float(rresp, "EPERDN", &status);
      drms_setkey_float(rec, "DN_GAIN", eperdn);
      float dnperpht = drms_getkey_float(rresp, "DNPERPHT", &status);
      drms_setkey_float(rec, "DNPERPHT", dnperpht);
      float eff_wl = drms_getkey_float(rresp, "EFF_WVLN", &status);
      drms_setkey_float(rec, "EFF_WVLN", eff_wl);
      float p1 = drms_getkey_float(rresp, "EFFA_P1", &status);
      float p2 = drms_getkey_float(rresp, "EFFA_P2", &status);
      float p3 = drms_getkey_float(rresp, "EFFA_P3", &status);
      TIME t_start = drms_getkey_float(rresp, "T_START", &status);
      float dt = (float) (t_obs - t_start)/86400.0;
      float factor = ((p3*dt + p2)*dt + p1)*dt + 1.0;
      float eff_area = drms_getkey_float(rresp, "EFF_AREA", &status);
      eff_area = eff_area*factor;
      drms_setkey_float(rec, "EFF_AREA", eff_area);
      int ver_num = drms_getkey_int(rresp, "VER_NUM", &status);
      drms_setkey_int(rec, "DN_GN_V", ver_num);
      drms_setkey_int(rec, "EFF_AR_V", ver_num);
    }
  } //end of "irec" loop
  if (drms_close_records(outRS, DRMS_INSERT_RECORD))
    fprintf(stderr, "drms_close_records failure for %s\n", dsSeries);

  fprintf(stdout, "%s CROTA2 fixed\n", dsSeries);
  fflush(stdout);

  return (DRMS_SUCCESS);
} // end of DoIt
