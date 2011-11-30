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

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "dsinp", NOT_SPECIFIED, "Series name w/optional record spec"},
  {ARG_STRING, "pathout", "/cache/sdo/test_jps", "out base path"},
  {ARG_STRING, "prefix", "AIA", "output filename prefix"},
  {ARG_STRING, "segment", "0", "segment number"},
  {ARG_STRING, "suffix", NOT_SPECIFIED, "output filename suffix"},
  {ARG_STRING, "skip", "0", "skip writing existing files"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "c", "0", "Compress output image"},
  {ARG_FLAG, "s", "0", "sleep flag"},
  {ARG_FLAG, "u", "0", "Uncompress output image"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_FLAG, "o", "0", "first segment flag"},
  {ARG_END}
};

char *module_name = "aia_lev0_fits";

int de_compress=0, do_compress=0, do_sleep = 0, verbose = 0, do_1segment = 0;
int skip = 0, segment;

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("aia_lev0_fits {-h} {-v} dsinp=<recordset query> pathout=name\n"
        "  -h: print this message\n"
        "  -c: compress image\n"
        "  -s: sleep between writes\n"
        "  -u: uncompress image\n"
        "  -v: verbose\n"
        "  -o: export only first segment\n"
        "skip=<[0|1]>\n"
        "dsinp=<recordset query> as <series>{[record specifier]} - required\n"
        "prefix=<output filename prefix> default is AIA\n"
        "segment=<segment number> default is 0\n"
        "suffix=<output filename suffix> default is wavelength as a number\n"
        "pathout=<output directory path>\n");
    return(1);
  }
  do_sleep = cmdparams_get_int(&cmdparams, "s", NULL) != 0;
  do_compress = cmdparams_get_int(&cmdparams, "c", NULL) != 0;
  de_compress = cmdparams_get_int(&cmdparams, "u", NULL) != 0;
  do_1segment = cmdparams_get_int(&cmdparams, "o", NULL) != 0;
  skip = cmdparams_get_int(&cmdparams, "skip", NULL);
  return(0);
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
      fprintf(stderr, "Cant make year diectory '%s'\n", outfilename);
      return 1;
    }
  }
  sprintf(outfilename, "%s/%4.4d/%2.2d", outpathbase, yr, mo);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make month diectory '%s'\n", outfilename);
      return 1;
    }
  }
  sprintf(outfilename, "%s/%4.4d/%2.2d/%2.2d", outpathbase, yr, mo, da);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make day diectory '%s'\n", outfilename);
      return 1;
    }
  }
  sprintf(outfilename, "%s/%4.4d/%2.2d/%2.2d/H%2.2d00",
                       outpathbase, yr, mo, da, hr);
  if (stat(outfilename, &sbuf)) {
    if (mkdir(outfilename, 0777)) {
      fprintf(stderr, "Cant make day diectory '%s'\n", outfilename);
      return 1;
    }
  }
  return 0;
}

int DoIt ()
{
  int irec, iseg, length[2], nrecs, nsegs, status=0, iwvl, ni=0, fsn, good_rec;
  int yr, mo, da, hr, mn, sc;
  int wvls[10] = { 171, 193, 335, 304, 1700, 211, 131, 94, 1600, 4500 };
  char *dsinp, *outpathbase, *date_obs, outfilename[512], *outcparms=NULL, *wu;
  char c1, c2, c3, *prefix, *suffix, *imtype, ds[3];
  time_t t0, tnow;
  DRMS_Record_t *inprec;
  DRMS_RecordSet_t *inprs;
  DRMS_Segment_t *inpseg;
  struct stat sbuf;
  TIME t_obs;

  if (nice_intro(0)) return(0);
  dsinp = strdup(cmdparams_get_str(&cmdparams, "dsinp", NULL));
  outpathbase = strdup(cmdparams_get_str(&cmdparams, "pathout", NULL));
  prefix = strdup(cmdparams_get_str(&cmdparams, "prefix", NULL));
  segment =  cmdparams_get_int(&cmdparams, "segment", NULL);
  suffix = strdup(cmdparams_get_str(&cmdparams, "suffix", NULL));
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("in argument is required");

  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  printf("%d records\n", nrecs);
  t0 = time(NULL);
  for (irec=0; irec<nrecs; irec++) {
    inprec = inprs->records[irec];
    nsegs = hcon_size(&inprec->segments);
    if (do_1segment) nsegs = 1;
    for (iseg=segment; iseg<segment + nsegs; iseg++) {
      char *filename, *inpsegname, outpath[512];
      float wavl;
      inpseg = drms_segment_lookupnum(inprec, iseg);
      inpsegname = inpseg->info->name;
      filename = inpseg->filename;
      if(strlen(filename)) {
        if (verbose) printf("rec %d, seg %d, protocol %d, '%s', '%s'\n",
               irec, iseg, inpseg->info->protocol, inpsegname, filename);
        wu = drms_getkey_string(inprec, "WAVEUNIT", &status);
        if (status && strcmp(suffix, NOT_SPECIFIED)==0) suffix = ""; 
        if (!strncasecmp(wu, "nm", 2)) {
          iwvl = (int)(drms_getkey_float(inprec, "WAVELNTH", &status)*10+.5);
        } else if (!strncasecmp(wu, "angstrom", 8)) {
          iwvl = drms_getkey_int(inprec, "WAVELNTH", &status);
        }
        fsn = drms_getkey_int(inprec, "FSN", &status);
        t_obs = drms_getkey_time(inprec, "T_OBS", &status);
        date_obs = drms_getkey_string(inprec, "T_OBS", &status);
        good_rec = 1;
        if ((t_obs<0.0) || (fsn == 0x1c001c00)) good_rec = 0;
        sscanf(date_obs, "%d%c%d%c%d%c%d:%d:%d",
                        &yr, &c1, &mo, &c2, &da, &c3, &hr, &mn, &sc);
        if (check_outpath(yr, mo, da, hr, outpathbase))DIE("Cant write to out");
        if (do_1segment && !iseg) inpsegname = "";
        imtype = drms_getkey_string(inprec, "IMG_TYPE", &status);
        if (strncmp(imtype, "DARK", 4)) strcpy(ds, "");
        else strcpy(ds, "d");
        sprintf(outpath, "%s/%4.4d/%2.2d/%2.2d/H%2.2d00",
                outpathbase, yr, mo, da, hr);
        if (strcmp(suffix, NOT_SPECIFIED)==0) {
          sprintf(outfilename,
            "%s/%s%4.4d%2.2d%2.2d_%2.2d%2.2d%2.2d_%4.4d%s%s.fits",
            outpath, prefix, yr, mo, da, hr, mn, sc, iwvl, inpsegname, ds);
        } else {
          sprintf(outfilename,
            "%s/%s%4.4d%2.2d%2.2d_%2.2d%2.2d%2.2d_%s%s%s.fits",
            outpath, prefix, yr, mo, da, hr, mn, sc, suffix, inpsegname, ds);
        }
        if ((stat(outfilename, &sbuf) || !skip) && good_rec) {
          if (do_compress) outcparms = "compress";
          if (de_compress) {
            outcparms = "";
          }
          if (DRMS_SUCCESS != fitsexport_export_tofile(inpseg, outcparms, outfilename)) {
            DIE("cant write FITS file");
          }
          if (do_sleep) {
            ni++;
            tnow = time(NULL);
            if ((tnow - t0) < (1.25*ni)) sleep(1);
          }
        }
      }
    }
  }
  return status;
}
