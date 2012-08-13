#include <string.h>
#include "jsoc_main.h"
#include "drms.h"
#include "fitsexport.h"
#include "drms_names.h"

#define NOT_SPECIFIED "***Not Specified***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

/*
asdf #include "imageRotate.c"
*/
int image_magrotate(void *, int, int, int, float, float, float, float,
                    void **, int *, int *, int, int);
ModuleArgs_t module_args[] =
{
  {ARG_STRING, "dsinp", NOT_SPECIFIED, "Input series query"},
  {ARG_STRING, "dsout", NOT_SPECIFIED, "Output series"},
  {ARG_STRING, "tbeg", NOT_SPECIFIED, "begin time"},
  {ARG_STRING, "pathout", "/scr21/jsoc/aia/synoptic", "out base path"},
  {ARG_STRING, "cadence", "180", "output cadence"},
  {ARG_STRING, "drms", "0", "write to drms series"},
  {ARG_STRING, "file", "1", "write to hdir FITS"},
  {ARG_STRING, "despike", "1", "remove cosmic ray hits"},
  {ARG_STRING, "rescale", "1", "rescale to fixed plate scale"},
  {ARG_STRING, "regrid", "1", "regrid type 0: nearest neighbor, 1: bicubic"},
  {ARG_STRING, "scale_to", "0.6", "rescale to fixed plate scale"},
  {ARG_STRING, "do_stretchmarks", "0", "fill in empty pixels created"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

TIME tbeg;
char *module_name = "aia_synoptic_nrt", *imtype;
int despike = 1, rescale = 1, regridtype = 1, verbose = 0, do_stretchmarks = 0,
    cadence = 180;
float scale_to = 0.6, scale_by;

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("aia_synoptic_nrt {-h} {-v} dsinp=series_record_spec dsout=output_series\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        "pathout=<output directory path>, default=/scr21/jsoc/aia/synoptic\n"
        "drms=1, write output to DRMS series, default=0\n"
        "file=1, write output to hdir FITS file, default=1\n"
        "despike=0, do not despike, default=1, despike\n"
        "rescale=0, do not rescale to fixed plate scale, default=1\n"
        "regrid=0, nearest neighbor, default=1, bicubic\n"
        "scale_to=0.6, rescale to fixed plate scale, default=0.6\n"
        "cadence=180, output cadence, default=90\n"
        "tbeg=2010.04.01_00:00:00, begin time\n"
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
  int got_all_wl=0,iwl, irec, iseg, nrecs, nsegs, status, n_used=1, wl, wls[10];
  int file, drms;
  char c1, c2, c3, *date_str, *dsinp, *dsout, now_str[100], *outpathbase=NULL;
  char linkname[512], outpath[512], outfilename[512], *outcparms="compress";
  float crpix1, crpix2, cdelt1, cdelt2, crota2, x0, y0;
  float mag = 1.0, dx, dy;
  double gaex_obs, gaey_obs, gaez_obs, haex_obs, haey_obs, haez_obs, tr_step;
  long long tr_index;
  TIME t_obs, t_rec, tr_epoch;
  DRMS_Record_t *inprec, *outrec;
  DRMS_RecordSet_t *inprs;
  DRMS_Keyword_t *inpkey = NULL, *outkey = NULL;
  DRMS_Array_t *inparr=NULL, *outarr=NULL;
  DRMS_Segment_t *inpseg, *outseg;
  HIterator_t *hit = NULL;
  if (nice_intro(0)) return(0);

  dsinp = strdup(cmdparams_get_str(&cmdparams, "dsinp", NULL));
  dsout = strdup(cmdparams_get_str(&cmdparams, "dsout", NULL));
  outpathbase = strdup(cmdparams_get_str(&cmdparams, "pathout", NULL));
  tbeg = cmdparams_get_time(&cmdparams, "tbeg", NULL);
  drms = cmdparams_get_int(&cmdparams, "drms", NULL);
  file = cmdparams_get_int(&cmdparams, "file", NULL);
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("dsinp argument is required");
  if (strcmp(dsout,NOT_SPECIFIED)==0) DIE("dsout argument is required");
/*
  if (strcmp(tbeg, NOT_SPECIFIED)==0) DIE("tbeg argument is required");
  if (strcmp(tend, NOT_SPECIFIED)==0) DIE("tend argument is required");
*/
  despike = cmdparams_get_int(&cmdparams, "despike", NULL);
  rescale = cmdparams_get_int(&cmdparams, "rescale", NULL);
  scale_to = cmdparams_get_float(&cmdparams, "scale_to", NULL);
  regridtype = cmdparams_get_int(&cmdparams, "regrid", NULL);
  do_stretchmarks = cmdparams_get_int(&cmdparams, "do_stretchmarks", NULL);

  for (iwl=0; iwl<10; iwl++) wls[iwl] = 0;
  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  printf("%d records\n", nrecs);
  for (irec=0; irec<nrecs; irec++) {
    int skip, yr, mo, da, hr, mn, sc;
    inprec = inprs->records[irec];
    imtype = drms_getkey_string(inprec, "IMG_TYPE", &status);
    if (!strncmp(imtype, "DARK", 4)) continue;
    t_obs = drms_getkey_time(inprec, "T_OBS", &status);
    if (status) DIE("T_OBS not found!");
    date_str = drms_getkey_string(inprec, "T_OBS", &status);
    sscanf(date_str, "%d%c%d%c%d%c%d:%d:%d",
           &yr, &c1, &mo, &c2, &da, &c3, &hr, &mn, &sc);
    if (file) {
      if (check_outpath(yr, mo, da, hr, outpathbase)) DIE("Cant write to out");
    }
    if (mn % 3) continue;
    if ((t_obs - tbeg) > (cadence*n_used)) {
      for (iwl=0; iwl<10; iwl++) wls[iwl] = 0;
      got_all_wl = 0;
      n_used++;
    }
    if (got_all_wl) continue;
    wl = drms_getkey_int(inprec, "WAVELNTH", &status);
    if (status) DIE("WAVELNTH not found!");
    skip =0;
    switch (wl) {
      case 94:   if (wls[0]) { skip = 1; } else { wls[0] = 1; }
                 break;
      case 131:  if (wls[1]) { skip = 1; } else { wls[1] = 1; }
                 break;
      case 171:  if (wls[2]) { skip = 1; } else { wls[2] = 1; }
                 break;
      case 193:  if (wls[3]) { skip = 1; } else { wls[3] = 1; }
                 break;
      case 211:  if (wls[4]) { skip = 1; } else { wls[4] = 1; }
                 break;
      case 304:  if (wls[5]) { skip = 1; } else { wls[5] = 1; }
                 break;
      case 335:  if (wls[6]) { skip = 1; } else { wls[6] = 1; }
                 break;
      case 1600: if (wls[7]) { skip = 1; } else { wls[7] = 1; }
                 break;
      case 1700: if (wls[8]) { skip = 1; } else { wls[8] = 1; }
                 break;
      case 4500: if (wls[9]) { skip = 1; } else { wls[9] = 1; }
                 break;
    }
    if (skip) continue;
    outrec = drms_create_record(drms_env, dsout, DRMS_PERMANENT, &status);
    if (status) DIE("cant create recordset");
    status = drms_copykeys(outrec, inprec, 0, kDRMS_KeyClass_Explicit);
    if (status) DIE("Error in drms_copykeys()");

    sprint_time_ISO(now_str, CURRENT_SYSTEM_TIME);
    drms_setkey_string(outrec, "DATE", now_str);
    if ( 0 == drms_setkey_time(outrec, "T_REC", t_obs)) {
      tr_index = drms_getkey_longlong(outrec, "T_REC_index", &status);
      if (status) DIE("T_REC_index not found!");
      tr_step =  drms_getkey_double(outrec, "T_REC_step", &status);
      if (status) DIE("T_REC_step not found!");
      tr_epoch = drms_getkey_time(outrec, "T_REC_epoch", &status);
      if (status) DIE("T_REC_epoch not found!");
      t_rec = tr_epoch + tr_index*tr_step;
      drms_setkey_time(outrec, "T_REC", t_rec);
      date_str = drms_getkey_string(outrec, "T_REC", &status);
    } else {
      tr_step = 180.0; tr_epoch = 0;
      tr_index = (t_obs - tr_epoch)/tr_step;
      t_rec = tr_epoch + tr_index*tr_step;
    }
    sscanf(date_str, "%d%c%d%c%d%c%d:%d",
           &yr, &c1, &mo, &c2, &da, &c3, &hr, &mn);
    sc = 0;
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
             "%s/AIA%4.4d%2.2d%2.2d_%2.2d%2.2d%2.2d_%4.4d.fits",
             outpath, yr, mo, da, hr, mn, sc, wl);
          if (stat(outfilename, &sbuf)) {
            outcparms = ""; // uncompressed for now
            if (DRMS_SUCCESS != fitsexport_export_tofile(outseg, outcparms,
              outfilename, NULL, NULL)) DIE("cant write FITS file");
/*
            sprintf(linkname, "%s/mostrecent/AIAsynoptic%4.4d.fits",
                    outpathbase, wl);
            unlink(linkname);
            if (symlink(outfilename, linkname)) {
              fprintf(stderr, "Can not link '%s' to '%s'.\n",
                      outfilename, linkname);
            }
*/
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

  return 0;
}
