/* DRMS module to calculate disk center by limb fitting */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "jsoc_main.h"
#include "drms.h"
#include "fitsexport.h"
#include "drms_names.h"

#define NOT_SPECIFIED "***Not Specified***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

int limbcompute(float *x, int nx, int ny, float xcguess, float ycguess, float rguess, float rrange,
  int limbmode, int useprevflag, float fwhm);

float rsun_offset(int wavelength);

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "dsinp", NOT_SPECIFIED, "Series name w/optional record spec"},
  {ARG_STRING, "dsout", NOT_SPECIFIED, "Output series"},
  {ARG_STRING, "segment", "0", "segment number"},
  {ARG_STRING, "sum", "1", "number of images to sum to reduce noise"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

extern float sdisk_xc, sdisk_yc, sdisk_r;
char *module_name = "limbfit_aia";
int verbose = 0;

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    fprintf(stderr, "limbfit_aia {-h} {-v} dsinp=series_record_spec dsout=out_series\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        "dsinp=<recordset query> as <series>{[record specifier]} - required\n"
        "segment=<segment number> default is 0\n"
        "sum=<number of images to sum> to reduce noise\n"
        );
    return(1);
  }
  return(0);
}

int DoIt ()
{
  char *dsinp, *outpathbase, *date_obs, outfilename[512];
  int limbmode=0, useprevflag=0;
  int fsn, good_rec, irec, isum=0, nfits=0, nrecs, segment, status=0, sum, wl;
  int nx, ny;
  int yr, mo, da, hr, mn, sc, ss;
  float rguess, xcguess, ycguess, rrange=100.0, fwhm=3.0, *sumarr=NULL;
  DRMS_Record_t *inprec;
  DRMS_RecordSet_t *inprs;
  DRMS_Segment_t *inpseg;
  DRMS_Array_t *inparr=NULL;
  struct stat sbuf;
  TIME t_obs;
  int i, ut;

  if (nice_intro(0)) return(0);
  dsinp = strdup(cmdparams_get_str(&cmdparams, "dsinp", NULL));
  segment =  cmdparams_get_int(&cmdparams, "segment", NULL);
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("in argument is required");
  sum = cmdparams_get_int(&cmdparams, "sum", NULL);
  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  fprintf(stderr, "%d records, sum: %d\n", nrecs, sum);
  for (irec=0; irec<nrecs; irec++) {
    char *filename, *inpsegname, outpath[512];
    inprec = inprs->records[irec];
    inpseg = drms_segment_lookupnum(inprec, segment);
    inpsegname = inpseg->info->name;
    filename = inpseg->filename;
    if(strlen(filename)) {
      if (verbose) printf("rec %d, protocol %d, '%s', '%s'\n",
               irec, inpseg->info->protocol, inpsegname, filename);
      wl = drms_getkey_int(inprec, "WAVELNTH", &status);
      fsn = drms_getkey_int(inprec, "FSN", &status);
      rguess = drms_getkey_float(inprec, "R_SUN", &status);
      rguess += rsun_offset(wl);
      xcguess = drms_getkey_float(inprec, "CRPIX1", &status);
      ycguess = drms_getkey_float(inprec, "CRPIX2", &status);
      t_obs = drms_getkey_time(inprec, "T_OBS", &status);
      date_obs = drms_getkey_string(inprec, "T_OBS", &status);
      if ((t_obs<0.0) || (fsn == 0x1c001c00)) continue;
      limbmode = 0; 
      if (wl == 304) limbmode = 3;
      else if (wl > 2000) limbmode = 2;
      else if (wl > 1000) limbmode = 1;
      inparr = drms_segment_read(inpseg, DRMS_TYPE_FLOAT, &status);
      if (status) DIE("drms_segment_read failed!");
      nx = inparr->axis[0]; ny = inparr->axis[1];
      if (isum) { float *values =  inparr->data;
        for (i=0; i<nx*ny; i++) sumarr[i] += values[i];
      } else {
        if (sum == 1) sumarr = inparr->data;
        else sumarr = (float*) calloc(nx*ny, sizeof(float));
        memcpy(sumarr, inparr->data, nx*ny*sizeof(float));
      }
      isum++;
      if (isum == sum) {
        limbcompute(sumarr, nx, ny, xcguess, ycguess, rguess, rrange,
                    limbmode, useprevflag, fwhm);
        if (sum>1) { free(sumarr); sumarr = NULL; }
        nfits++; isum = 0;
        ut = t_obs + 220924763;
        printf("%.2f\t%.2f\t%.2f\t%d\t%s\t%.2f\n", 
               sdisk_xc, sdisk_yc, sdisk_r, ut, date_obs, rguess);
      }
      if (inparr) drms_free_array(inparr);
    }
  }

  return status;
}
