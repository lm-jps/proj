#include <string.h>
#include "jsoc_main.h"
#include "drms.h"
#include "fitsexport.h"
#include "drms_names.h"

#define NOT_SPECIFIED "***Not Specified***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

int image_magrotate(void *, int, int, int, float, float, float, float,
                    void **, int *, int *, int, int);

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "dsinp", NOT_SPECIFIED, "Input series query"},
  {ARG_STRING, "dsout", NOT_SPECIFIED, "Output series"},
  {ARG_STRING, "pathout", "/scr21/jsoc/aia/synoptic", "out base path"},
  {ARG_STRING, "drms", "0", "write to drms series"},
  {ARG_STRING, "file", "1", "write to hdir FITS"},
  {ARG_STRING, "rescale", "1", "rescale to fixed plate scale"},
  {ARG_STRING, "regrid", "1", "regrid type 0: nearest neighbor, 1: bicubic"},
  {ARG_STRING, "scale_to", "0.6", "rescale to fixed plate scale"},
  {ARG_STRING, "do_stretchmarks", "0", "fill in empty pixels created"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

char *module_name = "aia_synoptic", *imtype;
int rescale = 1, regridtype = 1, verbose = 0, do_stretchmarks = 0;
float scale_to = 0.6, scale_by;

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("aia_lev1p5 {-h} {-v} dsinp=series_record_spec dsout=output_series\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        "pathout=<output directory path>, default=/scr21/jsoc/aia/synoptic\n"
        "drms=1, write output to DRMS series, default=0\n"
        "file=1, write output to hdir FITS file, default=1\n"
        "rescale=0, do not rescale to fixed plate scale, default=1\n"
        "regrid=0, nearest neighbor, default=1, bicubic\n"
        "scale_to=0.6, rescale to fixed plate scale, default=0.6\n"
        "do_stretchmarks=0, fill in empty pixels created, default=0\n"
        "dsinp=<recordset query> as <series>{[record specifier]} - required\n"
        "dsout=<series> - required if output to drms series\n");
    return(1);
  }
  return(0);
}
    
void sprint_time_ISO (char *tstring, TIME t)
{ 
  sprint_at(tstring,t);
  tstring[4] = tstring[7] = '-';
  tstring[10] = 'T';
  tstring[19] = '\0';
}

int check_outpath(int yr, int mo, int da, int hr, char *outpathbase)
{
  struct stat sbuf;
  char outfilename[512];

  if (stat(outpathbase, &sbuf)) {
    fprintf(stderr, "Cant stat output base path\n"); return 1;
  }
  sprintf(outfilename, "%s/%4.4d", outpathbase, yr);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make year directory '%s'\n", outfilename);
      return 1;
    }
  }
  sprintf(outfilename, "%s/%4.4d/%2.2d", outpathbase, yr, mo);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make month directory '%s'\n", outfilename);
      return 1;
    }
  }
  sprintf(outfilename, "%s/%4.4d/%2.2d/%2.2d", outpathbase, yr, mo, da);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make day directory '%s'\n", outfilename);
      return 1;
    }
  }
  sprintf(outfilename, "%s/%4.4d/%2.2d/%2.2d/H%2.2d00",
                       outpathbase, yr, mo, da, hr);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make day directory '%s'\n", outfilename);
      return 1;
    }
  }
  return 0;
}

int DoIt ()
{
  int irec, iseg, nrecs, nsegs, status, n_used=1, wl, explim, fsn, quality;
  int file, drms, cmdexp, nogood, oy, om, od, oh, on;
  char c1, c2, c3, *date_str, *dsinp, *dsout, now_str[100], *outpathbase=NULL;
  char querystr[512], outpath[512], outfilename[512], *outcparms="compress";
  float crpix1, crpix2, cdelt1, cdelt2, crota2, x0, y0;
  float mag = 1.0, dx, dy;
  double gaex_obs, gaey_obs, gaez_obs, haex_obs, haey_obs, haez_obs, tr_step;
  long long tr_index;
  FILE *fimg_time;
  TIME t_obs, t_rec, tr_epoch, t_obs_alt;
  DRMS_Record_t *inprec, *outrec;
  DRMS_RecordSet_t *inprs, *aecrs;
  DRMS_Keyword_t *inpkey = NULL, *outkey = NULL;
  DRMS_Array_t *inparr=NULL, *outarr=NULL;
  DRMS_Segment_t *inpseg, *outseg;
  HIterator_t *hit = NULL;
  if (nice_intro(0)) return(0);

  dsinp = strdup(cmdparams_get_str(&cmdparams, "dsinp", NULL));
  dsout = strdup(cmdparams_get_str(&cmdparams, "dsout", NULL));
  outpathbase = strdup(cmdparams_get_str(&cmdparams, "pathout", NULL));
  drms = cmdparams_get_int(&cmdparams, "drms", NULL);
  file = cmdparams_get_int(&cmdparams, "file", NULL);
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("dsinp argument is required");
  if (strcmp(dsout,NOT_SPECIFIED)==0) DIE("dsout argument is required");
  rescale = cmdparams_get_int(&cmdparams, "rescale", NULL);
  scale_to = cmdparams_get_float(&cmdparams, "scale_to", NULL);
  regridtype = cmdparams_get_int(&cmdparams, "regrid", NULL);
  do_stretchmarks = cmdparams_get_int(&cmdparams, "do_stretchmarks", NULL);

  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  printf("%d records\n", nrecs);
  for (irec=0; irec<nrecs; irec++) {
    int yr, mo, da, hr, mn, sc;
    inprec = inprs->records[irec];
    wl = drms_getkey_int(inprec, "WAVELNTH", &status);
    if (status) DIE("WAVELNTH not found!");
    switch (wl) {
      case 94:
      case 131:
      case 211:
      case 304:
      case 335:
      case 1600: explim = 2800; break;
      case 171:
      case 193: explim = 1931; break;
      case 1700: explim = 966; break;
      case 4500: explim = 483; break;
      default: DIE("Bad wavelength");
    }
    imtype = drms_getkey_string(inprec, "IMG_TYPE", &status);
    if (status) DIE("cant get IMG_TYPE");
    fsn = drms_getkey_int(inprec, "FSN", &status);
    if (status) DIE("cant get FSN");
    cmdexp = drms_getkey_int(inprec, "AIMGSHCE", &status);
    if (status) DIE("cant get commanded exposure AIMGSHCE");
    quality = drms_getkey_int(inprec, "QUALITY", &status);
    if (status) DIE("cant get QUALITY");
    t_obs = drms_getkey_time(inprec, "T_OBS", &status);
    if (status) DIE("T_OBS not found!");
    if (cmdexp < explim) {
      char *pos;
      strncpy(querystr, dsinp, 512);
      t_obs_alt = t_obs - 30.0;
      pos = index(querystr, '[');
      snprintf(pos, 512, "[? T_REC between %10f and %10f ?][? WAVELNTH = %d ?]",
               t_obs - 30.0, t_obs + 30.0, wl);
      nogood = 1;
printf("%s\n", querystr);
    } else nogood = 0;
    while (nogood) {
      nogood = 0;
    }
    if (nogood) continue;
    date_str = drms_getkey_string(inprec, "T_OBS", &status);
    sscanf(date_str, "%d%c%d%c%d%c%d:%d:%d",
           &yr, &c1, &mo, &c2, &da, &c3, &hr, &mn, &sc);
    if (file) {
      if (check_outpath(yr, mo, da, hr, outpathbase)) DIE("Cant write to out");
    }
    outrec = drms_create_record(drms_env, dsout, DRMS_PERMANENT, &status);
    if (status) DIE("cant create recordset");
    status = drms_copykeys(outrec, inprec, 0, kDRMS_KeyClass_Explicit);
    if (status) DIE("Error in drms_copykeys()");

    sprint_time_ISO(now_str, CURRENT_SYSTEM_TIME);
    drms_setkey_string(outrec, "DATE", now_str);
    t_rec = t_obs + 30.0;
    if ( 0 == drms_setkey_time(outrec, "T_REC", t_rec)) {
      date_str = drms_getkey_string(outrec, "T_REC", &status);
      t_rec = (long long) (t_obs+0.5);
      drms_setkey_time(outrec, "T_REC", t_rec);
    } else {
      tr_step = 120.0; tr_epoch = 0;
      tr_index = (t_obs - tr_epoch)/tr_step;
      t_rec = tr_epoch + tr_index*tr_step;
    }
    sscanf(date_str, "%d%c%d%c%d%c%d:%d",
           &yr, &c1, &mo, &c2, &da, &c3, &hr, &mn);
    if (file) {
      if (check_outpath(yr, mo, da, hr, outpathbase)) DIE("Cant write to out");
    }
    crpix1 = drms_getkey_float(inprec, "CRPIX1", &status);
    if (status) DIE("CRPIX1 not found!");
    crpix2 = drms_getkey_float(inprec, "CRPIX2", &status);
    if (status) DIE("CRPIX2 not found!");
    cdelt1 = drms_getkey_float(inprec, "CDELT1", &status);
    if (status) DIE("CDELT1 not found!");
    cdelt2 = drms_getkey_float(inprec, "CDELT2", &status);
    if (status) DIE("CDELT2 not found!");
    if (rescale) {
      mag = fabs(cdelt1 / scale_to);
      cdelt1 = 4.0*cdelt1/mag;
      cdelt2 = 4.0*cdelt2/mag;
      drms_setkey_float(outrec, "CDELT1", cdelt1);
      drms_setkey_float(outrec, "CDELT2", cdelt2);
    }
    crota2 = drms_getkey_float(inprec, "CROTA2", &status);
    if (status) DIE("CROTA2 not found!");
    x0 = drms_getkey_float(inprec, "X0", &status);
//       if (status) DIE("X0 not found!");
    y0 = drms_getkey_float(inprec, "Y0", &status);
//       if (status) DIE("Y0 not found!");
    if (hit) hiter_destroy(&hit);
    nsegs = hcon_size(&inprec->segments);
    for (iseg=0; iseg<1; iseg++) {
      void *output_array = NULL;
      int i, j, n, m, dtyp, nx, ny, outaxis[2];
      char *filename = NULL, *inpsegname = NULL;
      inpseg = drms_segment_lookupnum(inprec, iseg);
      inpsegname = inpseg->info->name;
      filename = inpseg->filename;
      if (0 == iseg) {
        outseg = drms_segment_lookupnum(outrec, iseg);
        if (!outseg) DIE("Cant get output segment");
        inparr = drms_segment_read(inpseg, DRMS_TYPE_FLOAT, &status);
        if (status) DIE("drms_segment_read failed!");
        dtyp = 3;
        n = inparr->axis[0]; m = inparr->axis[1];
//      dx = crpix1 - (n + 1.0)*0.5; dy = crpix2 - (m + 1.0)*0.5;
        dx = (n + 1.0)*0.5 - crpix1; dy = (m + 1.0)*0.5 - crpix2;

        status = image_magrotate( inparr->data, n, m, dtyp, crota2,
                 mag, dx, dy, &output_array, &nx, &ny,
                 regridtype, do_stretchmarks);
        if (status) DIE("image_magrotate failed!");
        outaxis[0] = 1024; outaxis[1] = 1024;
        outarr = drms_array_create(DRMS_TYPE_FLOAT, 2, outaxis,
                                   NULL, &status);
        if (status) DIE("drms_array_create failed!");
//      outarr->type = DRMS_TYPE_FLOAT;
        outarr->bscale = 0.0625;
        outarr->bzero = 0.0;
        for (i=0; i<1024; i++) {
          int k0 = 4*i;
          for (j=0; j<1024; j++) {
            int k, l, nn=0, l0=4*j;
            float sum=0.0;
            for (k=k0; k<k0+4; k++) for (l=l0; l<l0+4; l++) {
              int ndx = 4096*k + l;
              float pv = *((float *)(output_array)+ndx);
              if (!isnan(pv)) {sum += pv; nn++; }
            }
            *((float *)(outarr->data)+(1024*i+j)) = sum/nn;
          }
        }
        drms_setkey_float(outrec, "CROTA2", 0.0);
        drms_setkey_float(outrec, "CRPIX1", (1024 + 1.0)*0.5);
        drms_setkey_float(outrec, "CRPIX2", (1024 + 1.0)*0.5);
        drms_setkey_float(outrec, "X0", 0.0);
        drms_setkey_float(outrec, "Y0", 0.0);
        drms_segment_write(outseg, outarr, 0);
        if (file) {
          struct stat sbuf;
          sprintf(outpath, "%s/%4.4d/%2.2d/%2.2d/H%2.2d00",
                  outpathbase, yr, mo, da, hr);
          sprintf(outfilename,
             "%s/AIA%4.4d%2.2d%2.2d_%2.2d%2.2d_%4.4d.fits",
             outpath, yr, mo, da, hr, mn, wl);
          oy = yr; om = mo; od = da; oh = hr; on = mn;
          if (stat(outfilename, &sbuf)) {
            outcparms = "compress";
            if (DRMS_SUCCESS != fitsexport_export_tofile(outseg, outcparms,
              outfilename, NULL, NULL)) DIE("cant write FITS file");
          }
        }
        if (inparr) drms_free_array(inparr);
        if (output_array) free(output_array);
        if (outarr) drms_free_array(outarr);
      }
    }
    if (drms) drms_close_record(outrec, DRMS_INSERT_RECORD);
    else drms_close_record(outrec, DRMS_FREE_RECORD);
  }
  if (wl > 1500) return 0;
  sprintf(outpath, "%s/image_times", outpathbase);
  if (fimg_time = fopen(outpath, "w")) {
    fprintf(fimg_time, "Time     %d%2.2d%2.2d_%2.2d%2.2d00\n",
                                 oy, om,   od,  oh, on);
    fprintf(fimg_time, "synoptic  http://jsoc.stanford.edu/data/aia/synoptic/");
    fprintf(fimg_time, "%d/%2.2d/%2.2d/H%2.2d00/AIA%d%2.2d%2.2d_%2.2d%2.2d00\n"
           , oy, om, od, oh, oy, om, od, oh, on);
    fclose(fimg_time);
  } else fprintf(stderr, "Can't open image_times file.\n");
  sprintf(outpath, "%s/image_times.json", outpathbase);
  if (fimg_time = fopen(outpath, "w")) {
    fprintf(fimg_time, "{\"first\":\"20100513_000000\",\"last\":\"");
    fprintf(fimg_time, "%d%2.2d%2.2d_%2.2d%2.2d00\"}\n", oy, om, od, oh, on);
    fclose(fimg_time);
  } else fprintf(stderr, "Can't open image_times.json file.\n");

  return 0;
}
