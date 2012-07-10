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

int drms_appendcomment(DRMS_Record_t *rec, const char *comment, int nl);

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
  {ARG_FLAG, "c", "0", "crop image to RSUN_OBS-1/4 pixel"},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

char *module_name = "jsoc_resize";


// Definitions at top of file
// ---------------------------------------------------------------------

#define     Deg2Rad    (M_PI/180.0)
#define     Rad2arcsec (3600.0/Deg2Rad)
#define     arcsec2Rad (Deg2Rad/3600.0)
#define     Rad2Deg    (180.0/M_PI)

struct ObsInfo_struct
  {
  // from observation info
  TIME  t_obs;
  double rsun_obs, obs_vr, obs_vw, obs_vn;
  double crpix1, crpix2, cdelt1, cdelt2, crota2;
  double crval1, crval2;
  double cosa, sina;
  double obs_b0;
  // observed point
  int i,j;
  // parameters for observed point
  double x,y,r;
  double rho;
  double lon;
    double lat;
  double sinlat, coslat;
  double sig;
  double mu;
  double chi;
  double obs_v;
  };

typedef struct ObsInfo_struct ObsInfo_t;

void rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in);
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
int crop_image(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
char *get_input_recset(DRMS_Env_t *drms_env, char *in, TIME cadence);

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus);

// ---------------------------------------------------------------------

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
  int do_crop = cmdparams_get_int(&cmdparams, "c", NULL) != 0;

  double scale_by;
  char comment[4096];

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
    ObsInfo_t *ObsLoc;

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
    if (verbose)
      fprintf(stderr,"rec %d of %d, quality=%X\n",irec,nrecs,quality);
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
        else
          {
          ObsLoc = GetObsInfo(inpseg, NULL, &status);
          upNcenter(inparr, ObsLoc);
          if (do_crop)
            crop_image(inparr, ObsLoc);
          crota2 = ObsLoc->crota2;
          crpix1 = ObsLoc->crpix1;
          crpix2 = ObsLoc->crpix2;
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
        if (verbose)
          {
          fprintf(stderr,"image_magrotate called with: "
            "data=%p, n=%d, m=%d, dtyp=%d, crota2=%f, mag=%f, dx=%f, dy=%f, regridtype=%d, do_stretch=%d\n",
            (float *)inparr->data, n, m, dtyp, crota2, mag, dx, dy,regridtype, do_stretchmarks);
          fprintf(stderr,"image_magrotate returned: "
            "into out=%p, nx=%d, ny=%d, status=%d\n",  output_array, nx, ny, status);
          }
        sprintf(comment,"resize: scale_to=%f, mag=%f, dx=%f, dy=%f, regridtype=%d, do_stretch=%d",
          scale_to, mag, dx, dy,regridtype, do_stretchmarks);
        drms_appendcomment(outrec, comment, 0);
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
      else
        {
        if (verbose)
          fprintf(stderr,"WARNING - only first segment will be processed.\n");
        }
      drms_close_record(outrec, DRMS_INSERT_RECORD);
      }
    }
  return 0;
}



// definitions at bottom of file

// ----------------------------------------------------------------------

/* center whith whole pixel shifts and rotate by 180 if needed */
/* Only apply center if it will not result in an image crop.  I.e. not ever
   for AIA, and not for HMI or MDI or other if a shift of more than 20 arcsec
   is implied  */
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc)
  {
  int nx, ny, ix, iy, i, j, xoff, yoff, max_off;
  double rot, x0, y0, mid;
  float *data;
  if (!arr || !ObsLoc)
    return(1);
  data = arr->data;
  nx = arr->axis[0];
  ny = arr->axis[1];
  x0 = ObsLoc->crpix1 - 1;
  y0 = ObsLoc->crpix2 - 1;
  mid = (nx-1.0)/2.0;
  if ((rot = fabs(ObsLoc->crota2)) > 179 && rot < 181)
    {
    // rotate image by 180 degrees by a flip flip
    float val;
    int half = nx / 2;
    int odd = nx & 1;
    if (odd) half++;
    for (iy=0; iy<half; iy++)
      {
      for (ix=0; ix<nx; ix++)
        {
        i = iy*nx + ix;
        j = (ny - 1 - iy)*nx + (nx - 1 - ix);
        val = data[i];
        data[i] = data[j];
        data[j] = val;
        }
      }
    x0 = nx - 1 - x0;
    y0 = ny - 1 - y0;
    rot = ObsLoc->crota2 - 180.0;
    if (rot < -90.0) rot += 360.0;
    ObsLoc->crota2 = rot;
    }
  // Center to nearest pixel - if OK to do so
  xoff = round(x0 - mid);
  yoff = round(y0 - mid);
  max_off = 20.0 / ObsLoc->cdelt1;
  if (arr->parent_segment &&
      arr->parent_segment->record &&
      arr->parent_segment->record->seriesinfo && 
      arr->parent_segment->record->seriesinfo->seriesname && 
      strncasecmp(arr->parent_segment->record->seriesinfo->seriesname, "aia", 3) &&
      abs(xoff) < max_off && abs(yoff) < max_off) 
    {
    if (abs(xoff) >= 1)
      {
      for (iy=0; iy<ny; iy++)
        {
        float valarr[nx];
        for (ix=0; ix<nx; ix++)
          {
          int jx = ix + xoff;
          if (jx < nx && jx >= 0)
            valarr[ix] = data[iy*nx + jx];
          else
            valarr[ix] = DRMS_MISSING_FLOAT;
          }
        for (ix=0; ix<nx; ix++)
          data[iy*nx + ix] = valarr[ix];
        }
      x0 -= xoff;
      }
    if (abs(yoff) >= 1)
      {
      for (ix=0; ix<nx; ix++)
        {
        float valarr[ny];
        for (iy=0; iy<ny; iy++)
          {
          int jy = iy + yoff;
          if (jy < ny && jy >= 0)
            valarr[iy] = data[jy*nx + ix];
          else
            valarr[iy] = DRMS_MISSING_FLOAT;
          }
        for (iy=0; iy<ny; iy++)
          data[iy*nx + ix] = valarr[iy];
        }
      y0 -= yoff;
      }
    }
  // update center location
  ObsLoc->crpix1 = x0 + 1;
  ObsLoc->crpix2 = y0 + 1;
  return(0);
  }

// ----------------------------------------------------------------------

int crop_image(DRMS_Array_t *arr, ObsInfo_t *ObsLoc)
  {
  int nx, ny, ix, iy, i, j, xoff, yoff;
  double x0, y0;
  double rsun = ObsLoc->rsun_obs/ObsLoc->cdelt1;
  double scale, crop_limit2;
  float *data;
  if (!arr || !ObsLoc)
    return(1);
  data = arr->data;
  nx = arr->axis[0];
  ny = arr->axis[1];
  x0 = ObsLoc->crpix1 - 1;
  y0 = ObsLoc->crpix2 - 1;
  scale = 1.0/rsun;
  // crop_limit = 0.99975; // 1 - 1/4000, 1/2 HMI pixel.
  crop_limit2 = 0.99950; // square of 1 - 1/4000, 1/2 HMI pixel.
  for (iy=0; iy<ny; iy++)
    for (ix=0; ix<nx; ix++)
      {
      double x, y, R2;
      float *Ip = data + iy*nx + ix;
      if (drms_ismissing_float(*Ip))
        continue;
      x = ((double)ix - x0) * scale; /* x,y in pixel coords */
      y = ((double)iy - y0) * scale;
      R2 = x*x + y*y;
      if (R2 > crop_limit2)
        *Ip = DRMS_MISSING_FLOAT;
      }
  return(0);
  }

// ----------------------------------------------------------------------

#define CHECK(keyname) {if (status) {fprintf(stderr,"Keyword failure to find: %s, status=%d\n",keyname,status); *rstatus=status; return(NULL);}}

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus)
  {
  TIME t_prev;
  DRMS_Record_t *rec;
  TIME t_obs;
  double dv;
  ObsInfo_t *ObsLoc;
  int status;

  if (!seg || !(rec = seg->record))
    { *rstatus = 1; return(NULL); }

  ObsLoc = (pObsLoc ? pObsLoc : (ObsInfo_t *)malloc(sizeof(ObsInfo_t)));
  if (!pObsLoc)
    memset(ObsLoc, 0, sizeof(ObsInfo_t));

  t_prev = ObsLoc->t_obs;
  t_obs = drms_getkey_time(rec, "T_OBS", &status); CHECK("T_OBS");

  if (t_obs <= 0.0)
    { *rstatus = 2; return(NULL); }

  if (t_obs != t_prev)
    {
    ObsLoc->crpix1 = drms_getkey_double(rec, "CRPIX1", &status); CHECK("CRPIX1");
    ObsLoc->crpix2 = drms_getkey_double(rec, "CRPIX2", &status); CHECK("CRPIX2");
    ObsLoc->crval1 = drms_getkey_double(rec, "CRVAL1", &status); CHECK("CRVAL1");
    ObsLoc->crval2 = drms_getkey_double(rec, "CRVAL2", &status); CHECK("CRVAL2");
    ObsLoc->cdelt1 = drms_getkey_double(rec, "CDELT1", &status); CHECK("CDELT1");
    ObsLoc->cdelt2 = drms_getkey_double(rec, "CDELT2", &status); CHECK("CDELT1");
    ObsLoc->crota2 = drms_getkey_double(rec, "CROTA2", &status); CHECK("CROTA2");
    ObsLoc->sina = sin(ObsLoc->crota2*Deg2Rad);
    ObsLoc->cosa = sqrt (1.0 - ObsLoc->sina*ObsLoc->sina);
    ObsLoc->rsun_obs = drms_getkey_double(rec, "RSUN_OBS", &status);
    if (status)
      {
      double dsun_obs = drms_getkey_double(rec, "DSUN_OBS", &status); CHECK("DSUN_OBS");
      ObsLoc->rsun_obs = asin(696000000.0/dsun_obs)/arcsec2Rad;
      }
    ObsLoc->obs_vr = drms_getkey_double(rec, "OBS_VR", &status); CHECK("OBS_VR");
    ObsLoc->obs_vw = drms_getkey_double(rec, "OBS_VW", &status); CHECK("OBS_VW");
    ObsLoc->obs_vn = drms_getkey_double(rec, "OBS_VN", &status); CHECK("OBS_VN");
    ObsLoc->obs_b0 = drms_getkey_double(rec, "CRLT_OBS", &status); CHECK("CRLT_OBS");
    ObsLoc->t_obs = t_obs;
    }
  *rstatus = 0;
  return(ObsLoc);
  }

// ----------------------------------------------------------------------

// Generate explicit recordset list of closest good record to desired grid
// First get vector of times and quality
// Then if vector is not OK, quit.
// then: make temp file to hold recordset list
//       start with first time to define desired grid,
//       make array of desired times.
//       make empty array of recnums
//       search vector for good images nearest desired times
//       for each found time, write record query
char *get_input_recset(DRMS_Env_t *drms_env, char *in, TIME cadence)
  {
  TIME epoch = (cmdparams_exists(&cmdparams, "epoch")) ? params_get_time(&cmdparams, "epoch") : 0;
  DRMS_Array_t *data;
  TIME t_start, t_stop, t_now, t_want, t_diff, this_t_diff;
  int status;
  int nrecs, irec;
  int nslots, islot;
  long long *recnums;
  TIME *t_this, half = cadence/2.0;
  double *drecnum, *dquality;
  int quality;
  long long recnum;
  char keylist[DRMS_MAXQUERYLEN];
  static char filename[100];
  char *tmpdir;
  FILE *tmpfile;
  char *lbracket;
  char seriesname[DRMS_MAXQUERYLEN];

  sprintf(keylist, "T_OBS,QUALITY,recnum");
  data = drms_record_getvector(drms_env, in, keylist, DRMS_TYPE_DOUBLE, 0, &status);
  if (!data || status)
        {
        fprintf(stderr, "getkey_vector failed status=%d\n", status);
        return(NULL);
        }
  nrecs = data->axis[1];
  irec = 0;
  t_this = (TIME *)data->data;
  dquality = (double *)data->data + 1*nrecs;
  drecnum = (double *)data->data + 2*nrecs;
  if (epoch > 0.0)
    {
    int s0 = (t_this[0] - epoch)/cadence;
    TIME t0 = s0*cadence + epoch;
    t_start = (t0 < t_this[0] ? t0 + cadence : t0);
    }
  else
    t_start = t_this[0];
  t_stop = t_this[nrecs-1];
  nslots = (t_stop - t_start + cadence/2)/cadence;
  recnums = (long long *)malloc(nslots*sizeof(long long));
  for (islot=0; islot<nslots; islot++)
    recnums[islot] = 0;
  islot = 0;
  t_want = t_start;
  t_diff = 1.0e9;
  for (irec = 0; irec<nrecs; irec++)
      {
      t_now = t_this[irec];
      quality = (int)dquality[irec] & 0xFFFFFFFF;
      recnum = (long long)drecnum[irec];
      this_t_diff = fabs(t_now - t_want);
      if (quality < 0)
        continue;
      if (t_now <= (t_want-half))
        continue;
      while (t_now > (t_want+half))
        {
        islot++;
        if (islot >= nslots)
           break;
        t_want = t_start + cadence * islot;
        this_t_diff = fabs(t_now - t_want);
        t_diff = 1.0e8;
        }
      if (this_t_diff <= t_diff)
        recnums[islot] = recnum;
      t_diff = fabs(t_now - t_want);
      }
  if (islot+1 < nslots)
    nslots = islot+1;  // take what we got.
  strcpy(seriesname, in);
  lbracket = index(seriesname,'[');
  if (lbracket) *lbracket = '\0';
  tmpdir = getenv("TMPDIR");
  if (!tmpdir) tmpdir = "/tmp";
  sprintf(filename, "%s/hg_patchXXXXXX", tmpdir);
  mkstemp(filename);
  tmpfile = fopen(filename,"w");
  for (islot=0; islot<nslots; islot++)
    if (recnums[islot])
      fprintf(tmpfile, "%s[:#%lld]\n", seriesname, recnums[islot]);
  fclose(tmpfile);
  free(recnums);
  drms_free_array(data);
  return(filename);
  }
