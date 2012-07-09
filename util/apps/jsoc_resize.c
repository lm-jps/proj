/* adapted from proj/lev1.5_aia/aia_lev1p5.c */
/* restrictions:  should not resize by "too much" since does simple bicubic mapping
   and large changes will not preserve all information.  Also, the rotation is done
   about CRPIX1 and CRPIX2 which must be contained in the image, or at least
   the final location must be contained in the image after rotation.  A patch from
   HMI will result in all missing data.
 */

#include <string.h>
#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"

#define NOT_SPECIFIED "***Not Specified***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

#define NBINS 1048576
static int hist[NBINS];

int image_magrotate(void *, int nin, int min, int data_type_input, float theta, float mag,
    float dx, float dy, void **outarray, int *nx, int *ny, int regridtype_input, int stretch_lines);

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "in", NOT_SPECIFIED, "Input series query"},
  {ARG_STRING, "out", NOT_SPECIFIED, "Output series"},
  // {ARG_STRING, "despike", "0", "remove cosmic ray hits"},
  {ARG_INT, "rescale", "0", "rescale to fixed plate scale"},
  {ARG_INT, "regrid", "1", "regrid type 0: nearest neighbor, 1: bicubic"},
  {ARG_INT, "center_to", "0", "center to 0: Sun center, 1: First image"},
  {ARG_DOUBLE, "scale_to", "0.6", "rescale to new CDELT scale"},
  {ARG_STRING, "do_stretchmarks", "0", "replicate pixels created else use missing value"},
  {ARG_STRING, "requestid", NOT_SPECIFIED, "RequestID if used for export"},
  {ARG_FLAG, "c", "0", "Center Sun in frame, else use first image as reference"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

char *module_name = "jsoc_resize";

int verbose; 

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("jsoc_resize {-h} {-v} dsinp=series_record_spec dsout=output_series\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        // "despike=0, do not despike, default=0, despike\n"
        "rescale=0, do not rescale to fixed plate scale, default=0\n"
        "regrid=0, nearest neighbor, default=1, bicubic\n"
        "center_to=0, centering target, 0=Sun, 1=First_image; default=0\n"
        "scale_to=0.6, rescale to fixed plate scale, default=0.6 arcsec/pixel\n"
        "do_stretchmarks=0, fill in empty pixels created, default=0\n"
        "in=<recordset query> as <series>{[record specifier]} - required\n"
        "out=<series> - required\n");
    return(1);
  }
  return(0);
}
    
int set_statistics(DRMS_Segment_t *seg, DRMS_Array_t *arr, int mode);

int DoIt ()
{
  if (nice_intro(0)) return(0);

  int irec, iseg, nrecs, nsegs, status, is_aia=0;
  char now_str[100];
  float crpix1, crpix2, cdelt1, cdelt2, crota2, x0, y0;
  float mag = 1.0, dx, dy;
  DRMS_Record_t *inprec, *outrec;
  DRMS_RecordSet_t *inprs;
  DRMS_Keyword_t *inpkey = NULL, *outkey = NULL;
  DRMS_Array_t *inparr=NULL, *outarr=NULL;
  DRMS_Segment_t *inpseg, *outseg;

  const char *requestid = cmdparams_get_str(&cmdparams, "requestid", NULL);
  char *dsinp = strdup(cmdparams_get_str(&cmdparams, "in", NULL));
  const char *dsout = cmdparams_get_str(&cmdparams, "out", NULL);
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("\"in\" argument is required");
  if (strcmp(dsout, NOT_SPECIFIED)==0)
    {
    char newout[DRMS_MAXNAMELEN];
    char *inp, inseries[DRMS_MAXNAMELEN];
    strncpy(inseries, dsinp, DRMS_MAXNAMELEN);
    inp = index(inseries, '[');
    if (inp)
      *inp = '\0';
    strncpy(newout, inseries, DRMS_MAXNAMELEN);
    strncat(newout, "_resized", DRMS_MAXNAMELEN);
    dsout = strdup(newout);
    }
  else
    dsout = strdup(dsout);
  // int despike = cmdparams_get_int(&cmdparams, "despike", NULL);
  int rescale = cmdparams_get_int(&cmdparams, "rescale", NULL);
  double scale_to = cmdparams_get_double(&cmdparams, "scale_to", NULL);
  int regridtype = cmdparams_get_int(&cmdparams, "regrid", NULL);
  int center_to = cmdparams_get_int(&cmdparams, "center_to", NULL);
  int do_stretchmarks = cmdparams_get_int(&cmdparams, "do_stretchmarks", NULL);
  double scale_by;

  if (strstr(dsinp, "aia")) is_aia = 1;
  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  printf("%d records\n", nrecs);
float usedx = 0;
float usedy = 0;
float usemag = 1.0;
  for (irec=0; irec<nrecs; irec++)
    {
    int quality;
    TIME t_rec, t_obs;

    inprec = inprs->records[irec];
    t_obs = drms_getkey_time(inprec, "T_OBS", &status); if (status) DIE("T_OBS not found!");
    if (t_obs < 0)
      continue;
    if (is_aia)
      {
      int qualmask = 0x800100f0;
      int quallev0;
      quallev0 = drms_getkey_int(inprec, "QUALLEV0", &status); if (status) DIE("QUALLEV0 not found!");
      if (quallev0 & qualmask)
        continue;
      }
    quality = drms_getkey_int(inprec, "QUALITY", &status);
    if (quality < 0)
      continue;

    outrec = drms_create_record(drms_env, (char*)dsout, DRMS_PERMANENT, &status);
      if (status) DIE("cant create recordset");
    status = drms_copykeys(outrec, inprec, 0, kDRMS_KeyClass_Explicit);
      if (status) DIE("Error in drms_copykeys()");
    if (strcmp(requestid, NOT_SPECIFIED) != 0)
      drms_setkey_string(outrec, "requestid", requestid);

    drms_setkey_time(outrec, "DATE", CURRENT_SYSTEM_TIME);

    if (is_aia)
      {
      if (!drms_keyword_lookup(inprec, "T_REC", 1) && drms_keyword_lookup(outrec, "T_REC", 1))
        {
        double tr_step;
        long long tr_index;
        TIME tr_epoch;
        if ( 0 == drms_setkey_time(outrec, "T_REC", t_obs))
          {
          tr_index = drms_getkey_longlong(outrec, "T_REC_index", &status); if (status) DIE("T_REC_index not found!");
          tr_step =  drms_getkey_double(outrec, "T_REC_step", &status); if (status) DIE("T_REC_step not found!");
          tr_epoch = drms_getkey_time(outrec, "T_REC_epoch", &status); if (status) DIE("T_REC_epoch not found!");
          t_rec = tr_epoch + tr_index*tr_step;
          drms_setkey_time(outrec, "T_REC", t_rec);
          }
        }
      if (!drms_keyword_lookup(inprec, "GAEX_OBS", 1))
        drms_setkey_double(outrec, "GAEX_OBS", drms_getkey_double(inprec, "GCIEC_X", NULL));
      if (!drms_keyword_lookup(inprec, "GAEY_OBS", 1))
        drms_setkey_double(outrec, "GAEY_OBS", drms_getkey_double(inprec, "GCIEC_Y", NULL));
      if (!drms_keyword_lookup(inprec, "GAEZ_OBS", 1))
        drms_setkey_double(outrec, "GAEZ_OBS", drms_getkey_double(inprec, "GCIEC_Z", NULL));
      if (!drms_keyword_lookup(inprec, "HAEX_OBS", 1))
        drms_setkey_double(outrec, "HAEX_OBS", drms_getkey_double(inprec, "HCIEC_X", NULL));
      if (!drms_keyword_lookup(inprec, "HAEY_OBS", 1))
        drms_setkey_double(outrec, "HAEY_OBS", drms_getkey_double(inprec, "HCIEC_Y", NULL));
      if (!drms_keyword_lookup(inprec, "HAEZ_OBS", 1))
        drms_setkey_double(outrec, "HAEZ_OBS", drms_getkey_double(inprec, "HCIEC_Z", NULL));
      }

    crpix1 = drms_getkey_double(inprec, "CRPIX1", &status); if (status) DIE("CRPIX1 not found!");
    crpix2 = drms_getkey_double(inprec, "CRPIX2", &status); if (status) DIE("CRPIX2 not found!");
    cdelt1 = drms_getkey_double(inprec, "CDELT1", &status); if (status) DIE("CDELT1 not found!");
    if (rescale)
      {
      mag = fabs(cdelt1 / scale_to);
      cdelt1 /= mag;
      }
    crota2 = drms_getkey_double(inprec, "CROTA2", &status); if (status) DIE("CROTA2 not found!");
   
    nsegs = hcon_size(&inprec->segments);
    for (iseg=0; iseg<1; iseg++)
      {
      void *output_array = NULL;
      int i, ix, iy, n, m, dtyp, nx, ny, npix;
      inpseg = drms_segment_lookupnum(inprec, iseg);
      if (0 == iseg) // XXXXXX WARNING - for now only does first segment
        {
        inparr = drms_segment_read(inpseg, DRMS_TYPE_FLOAT, &status); if (status) DIE("drms_segment_read failed!");
        outseg = drms_segment_lookupnum(outrec, iseg); if (!outseg) DIE("Cant get output segment");
        int outdims[2];
        dtyp = 3;
        n = inparr->axis[0];
        m = inparr->axis[1];
 
        if (is_aia) // set outer rows and cols to 0
          {
          float *fval;
          fval = inparr->data;
          for (ix = 0; ix<n; ix++)
            {
            fval[ix] = 0.0;
            fval[(m - 1)*n + ix] = 0.0;
            }
          for (iy = 0; iy<m; iy++)
            {
            fval[iy*n] = 0.0;
            fval[iy*n + n - 1] = 0.0;
            }
          }

        if (center_to == 0)
          {
          usedx = (n + 1.0)*0.5;
          usedy = (m + 1.0)*0.5;
          }
        else
          {
          if (irec == 0)
            {
            usedx = crpix1;
            usedy = crpix2;
            usemag = cdelt1;
            }
          mag = cdelt1/usemag;
          cdelt1 = usemag;
          }
        dx = usedx - crpix1;
        dy = usedy - crpix2;
        status = image_magrotate( (float *)inparr->data, n, m, dtyp, crota2,
                 mag, dx, dy, &output_array, &nx, &ny,
		 regridtype, do_stretchmarks);
        if (status) DIE("image_magrotate failed!");
        outdims[0] = nx;
        outdims[1] = ny;
        outarr = drms_array_create(inparr->type, 2, outdims, output_array, &status);
        if (status) DIE("drms_array_create failed!");

        drms_setkey_float(outrec, "CROTA2", 0.0);
        drms_setkey_float(outrec, "CRPIX1", (nx + 1.0)*0.5);
        drms_setkey_float(outrec, "CRPIX2", (ny + 1.0)*0.5);
        drms_setkey_float(outrec, "CDELT1", cdelt1);
        drms_setkey_float(outrec, "CDELT2", cdelt1);
        set_statistics(outseg, outarr, 1);
        // get info for array from segment
        outarr->bzero = outseg->bzero;
        outarr->bscale = outseg->bscale;
        if (requestid == NOT_SPECIFIED)
          drms_segment_write(outseg, outarr, 0);
        else
          drms_segment_writewithkeys(outseg, outarr, 0);
        if (inparr) drms_free_array(inparr);
        if (outarr) drms_free_array(outarr);
        }
      drms_close_record(outrec, DRMS_INSERT_RECORD);
      }
    }
  return 0;
}
