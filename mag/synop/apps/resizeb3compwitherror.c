/**
   @defgroup jsoc_rebin jsoc_rebin reduce/increase image size by integer multiples
   @ingroup su_util

   @brief Reduce (increase) the dimension of the input data by an integer factor.

   @par Synopsis:
   @code
   jsoc_rebin  in=input_data out=output_data  scale=<scale> method={boxcar,gaussian}
   where scale is <1 for size reduction and >1 for size increase.  The output image scale
   will be the nearest integer to "scale" or its reciprocal for scale < 1.0. 
   @endcode

   This is a general purpose module that takes a series of input 
   data and modifies its spatial size by a factor 
   specified by "scale".
   The method for avaraging (interpolation) can be specified 
   through the input "method". The current version handles a simple 
   boxcar average and Gaussian filtered sampling. If 'scale' < 1 then the input is reduced in size.

   The image is not registered to solar center by this module.  The image will be rotated by
   a flip-flip procedure is the CROTA2 parameter is near +-180.0 unless the -u flag (unchanged) is
   present.  If the -c flag is present the image will be cropped at the solar limb before scaling.
   If the -h flag or requestid parameter is present the output segments will have full headers.

   If gaussian smoothing is specified (via method) then the FWHM and nvector parameters should also
   be provided.  These define the Gaussian smoothing vector.  nvector is the full width of the
   smoothing function.  nvector will be adjusted such that nvector is odd if 1/scale rounds to an odd integer.

   @par Flags:
   @c
   -c  Crop before scaling.  Use rsun_obs/cdelt1 for limb radius.
   -h  Write full FITS headers.
   -u  Leave unchanged, do NOT rotate by 180 degrees if CROTA2 is near 180.  default is to do flip-flip method so
       image is norths up and no pixel values are changed.

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in  The input data series.
   @param out The output series.

   @par Exit_Status:

   @par Example:
   Takes a series of 1024 X 1024 MDI full disk magnetogram and
   produces images with the resolution of HMI, 4096 X 4096.

   @code
   jsoc_rebin in='mdi.fd_M_96m_lev18[2003.10.20/1d]' out='su_phil.mdi_M_4k' scale=4 method='boxcar'
   @endcode

   @par Example:
   Reduces the resolution of HMI images 4096 X 4096 to that of MDI,
   1024 X 1024. Here the input is the HMI images and the output
   is the lower resolution HMI images.  Crop and rotate before rescaling.
   @code
   jsoc_rebin  -c in=hmi.M_45s[2011.10.20/1d]' out='su_phil.hmi_M_1k_45s' scale=0.25 method='boxcar'
   @endcode

   @bug
   None known so far.

*/

#include "jsoc.h"
#include "jsoc_main.h"
#include "fstats.h"

char *module_name = "resizeb3compwitherror";

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
     {ARG_STRING, "out", "NOT SPECIFIED",  "Output data series."},
     {ARG_FLAG, "c", "0", "Crop at rsun_obs."},
     {ARG_FLAG, "h", "0", "Include full FITS header in output segment."},
     {ARG_FLAG, "u", "0", "do not rotate by 180 if needed."},
     {ARG_FLOAT, "scale", "1.0", "Scale factor."},
     {ARG_FLOAT, "FWHM", "-1.0", "Smoothing Gaussian FWHM for method=gaussian."},
     {ARG_INT, "nvector", "-1.0", "Smoothing Gaussian vector length for method=gaussian."},
     {ARG_INT, "inseg", "0", "Input segment number"},
     {ARG_INT, "outseg", "0", "Output segment number"},
     {ARG_STRING, "method", "boxcar", "conversion type, one of: boxcar, gaussian."},
     {ARG_STRING, "requestid", "NA", "RequestID if called as an export processing step."},
     {ARG_END}
};

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
//  int mapmmax, sinbdivs; //Added by YL
  };

typedef struct ObsInfo_struct ObsInfo_t;

void rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in);
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
int crop_image(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
const char *get_input_recset(DRMS_Env_t *drms_env, const char *in);

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus);
ObsInfo_t *GetMinObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus);

int DoIt(void)
  {
  int status = DRMS_SUCCESS;
  DRMS_RecordSet_t *inRS, *outRS;
  int irec, nrecs;
  const char *inQuery = params_get_str(&cmdparams, "in");
  char *inStr;
  const char *outSeries = params_get_str(&cmdparams, "out");
  const char *method = params_get_str(&cmdparams, "method");
  const char *requestid = params_get_str(&cmdparams, "requestid");
  int nvector = params_get_int(&cmdparams, "nvector");
  float fscale = params_get_float(&cmdparams, "scale");
  float fwhm = params_get_float(&cmdparams, "FWHM");
  int crop = params_get_int(&cmdparams, "c");
  int as_is = params_get_int(&cmdparams, "u");
  int inseg = params_get_int(&cmdparams, "inseg");
  int outseg = params_get_int(&cmdparams, "outseg");
  int full_header = params_get_int(&cmdparams, "h") || strcmp(requestid, "NA");
  char *in_filename = NULL;

  int iscale, ivec, vec_half;
  double *vector;
  char history[4096];
  int mapmmax, sinbdivs; //Added by YL

  if (fscale < 1.0) // shrinking
    iscale = 1.0/fscale + 0.5;
  else  // enlarging
    iscale = fscale + 0.5;
  if (nvector < 0)
    nvector = iscale;
  // Both 1/scale and nvector must be odd or both even so add 1 to nvector if needed       
  if (((iscale & 1) && !(nvector & 1)) || ((!(iscale & 1) && (nvector & 1) )))
    nvector += 1;
  vector = (double *)malloc(nvector * sizeof(double));
  vec_half = nvector/2; // counts on truncate to int if nvector is odd.

  if (strcasecmp(method, "boxcar")==0 && fscale < 1)
    {
    for (ivec = 0; ivec < nvector; ivec++)
      vector[ivec] = 1.0;
    sprintf(history, "Boxcar bin by %d%s%s",
      iscale,
      (crop ? ", Cropped at rsun_obs" : ""),
      (!as_is ? ", North is up" : "") );
    }
  else if (strcasecmp(method, "boxcar")==0 && fscale >= 1)
    {
    if (nvector != iscale)
      DIE("For fscale>=1 nvector must be fscale");
    for (ivec = 0; ivec < nvector; ivec++)
      vector[ivec] = 1.0;
    sprintf(history, "Replicate to expand by %d%s%s",
      iscale,
      (crop ? ", Cropped at rsun_obs" : ""),
      (!as_is ? ", North is up" : "") );
    }
  else if (strcasecmp(method, "gaussian")==0) // do 2-D vector weights calculated as Gaussian
    {
    if (fwhm < 0)
      DIE("Need FWHM parameter");
    for (ivec = 0; ivec < nvector; ivec++)
      {
      double arg = (ivec - (nvector-1)/2.0) * (ivec - (nvector-1)/2.0);
      vector[ivec] = exp(-arg/(fwhm*fwhm*0.52034));
      }
    sprintf(history, "Scale by %f with Gasussian smoothing FWHM=%f, nvector=%d%s%s",
      fscale, fwhm, nvector,
      (crop ? ", Cropped at rsun_obs" : ""),
      (!as_is ? ", North is up" : "") );
    }
  else
    DIE("invalid conversion method");

  inStr = strdup(get_input_recset(drms_env, (char *)inQuery));
  if (!inStr || *inStr=='\0') DIE("Cant make special cadence recordset list file");
  inRS = drms_open_records(drms_env, inStr, &status);
  if (strcmp(inStr, inQuery) && *inStr == '@')
    unlink(inStr+1);
  if (status || inRS->n == 0)
    DIE("No input data found");

  nrecs = inRS->n;
  if (nrecs == 0)
    DIE("No records found");
  drms_stage_records(inRS, 1, 1);

  outRS = drms_create_records(drms_env, nrecs, (char *)outSeries, DRMS_PERMANENT, &status);
  if (status)
    DIE("Output recordset not created");

  for (irec=0; irec<nrecs; irec++)
    {
    ObsInfo_t *ObsLoc;
    DRMS_Record_t *outRec, *inRec;
    DRMS_Segment_t *outSeg, *inSeg;
    DRMS_Array_t *inArray, *outArray;
    float *inData, *outData;
    float val;
    int quality=0;
    int vcomp; // YL
// YL moves section to here 
      int inx, iny, outx, outy, i, j;
      int in_nx, in_ny, out_nx, out_ny, outDims[2];
// YL end of moving

    inRec = inRS->records[irec];
    outRec = outRS->records[irec]; // YL 
    quality = drms_getkey_int(inRec, "QUALITY", &status);
    if (status || (!status && quality >= 0))
      {
      for (vcomp=0; vcomp<4; vcomp++) {  // YL -- start vcomp loop for vector B
      double power = 1.0; // YL 
      inSeg = drms_segment_lookupnum(inRec, vcomp);
      inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
      if (vcomp == 3) power = 2.0; // YL

      if (status)
        {
        printf(" No data file found but QUALITY not bad, status=%d\n", status);
        drms_free_array(inArray);
        continue;
        }
      if (crop || !as_is)
        {
        ObsLoc = GetObsInfo(inSeg, NULL, &status);
        if (!as_is) upNcenter(inArray, ObsLoc);
        if (crop) crop_image(inArray, ObsLoc);
        }
      else
        {
        ObsLoc = GetMinObsInfo(inSeg, NULL, &status);
        }

      in_nx = inArray->axis[0]; //YL
      in_ny = inArray->axis[1]; //YL
      out_nx = in_nx * fscale + 0.5; //LY
      out_ny = in_ny * fscale + 0.5; //YL
      outDims[0] = out_nx; outDims[1] = out_ny; //YL
      outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
      inData = (float *)inArray->data; //YL
      outData = (float *)outArray->data; //YL

      if (fscale > 1.0)
        {
        int out_go = (iscale-1)/2.0 + 0.5;
        for (iny = 0; iny < in_ny; iny += 1)
          for (inx = 0; inx < in_nx; inx += 1)
            {
            val = inData[in_nx*iny + inx];
            for (j = 0; j < nvector; j += 1)
              {
              outy = iny*iscale + out_go + j - vec_half;
              for (i = 0; i < nvector; i += 1)
                {
                outx = inx*iscale + out_go + i - vec_half;
                if (outx >= 0 && outx < out_nx && outy >= 0 && outy < out_ny)
                  outData[out_nx*outy + outx] = val;
                }
              }
            }
        }
      else
        {
        int in_go = (iscale-1)/2.0 + 0.5;
        for (outy = 0; outy < out_ny; outy += 1)
          for (outx = 0; outx < out_nx; outx += 1)
            {
            double total = 0.0;
            double weight = 0.0;
            int nn = 0;
            for (j = 0; j < nvector; j += 1)
              {
              iny = outy*iscale + in_go + j - vec_half;
              for (i = 0; i < nvector; i += 1)
                {
                inx = outx*iscale + in_go + i - vec_half;
                if (inx >= 0 && inx < in_nx && iny >=0 && iny < in_ny)
                  {
                  val = inData[in_nx*(iny) + inx];
                  if (!drms_ismissing_float(val))
                    {
                    double w = vector[i]*vector[j];
                    total += pow(w*val, power); // YL
                    weight += w;
                    nn++;
                    }
                  }
                }
              }
            if (vcomp == 3) total = sqrt(total); // YL
            outData[out_nx*outy + outx] = (nn > 0 ? total/weight : DRMS_MISSING_FLOAT); 
            }
        }

      // Use the input array as the best guess for scale and zero
      outArray->bzero = inArray->bzero;
      outArray->bscale = inArray->bscale;
  
      drms_free_array(inArray);
  
      // write data file
        outSeg = drms_segment_lookupnum(outRec, vcomp); // YL
        outArray->parent_segment = outSeg;
        status = drms_segment_write(outSeg, outArray, 0);
      if (status)
        DIE("problem writing file");
      } // YL -- end vcomp loop

      // copy all keywords
      drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);

      // Now fixup coordinate keywords
      // Only CRPIX1,2 and CDELT1,2 and CROTA2 should need repair.
      drms_setkey_double(outRec, "CDELT1", ObsLoc->cdelt1/fscale);
      drms_setkey_double(outRec, "CDELT2", ObsLoc->cdelt2/fscale);
      drms_setkey_double(outRec, "CRPIX1", (ObsLoc->crpix1-0.5) * fscale + 0.5);
      drms_setkey_double(outRec, "CRPIX2", (ObsLoc->crpix2-0.5) * fscale + 0.5);
      drms_setkey_double(outRec, "CROTA2", ObsLoc->crota2);
      mapmmax = drms_getkey_int(inRec, "MAPMMAX", &status);
      drms_setkey_int(outRec, "MAPMMAX", (mapmmax-0.5) * fscale + 0.5);
      sinbdivs = drms_getkey_int(inRec, "SINBDIVS", &status);
      drms_setkey_int(outRec, "SINBDIVS", (sinbdivs-0.5) * fscale + 0.5);

      if (strcasecmp(method, "gaussian")==0) 
        drms_setkey_double(outRec, "FWHM", fwhm);
      drms_appendhistory(outRec, history, 1);
      drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
      if (strcmp(requestid, "NA") != 0)
        drms_setkey_string(outRec, "RequestID", requestid);
      drms_setkey_int(outRec, "TOTVALS_1", out_nx*out_ny);
      set_statistics(outSeg, outArray, 1);
      drms_free_array(outArray);

      }
    } //end of "irec" loop

  drms_close_records(inRS, DRMS_FREE_RECORD);
  drms_close_records(outRS, DRMS_INSERT_RECORD);
  return (DRMS_SUCCESS);
  } // end of DoIt

// ----------------------------------------------------------------------

/* center whith whole pixel shifts and rotate by 180 if needed */
/* Only apply center if it will not result in an image crop.  I.e. not ever
   for AIA, and not for HMI or MDI or other if a shift of more than 20 arcsec
   is implied  */
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc)
  {
  int nx, ny, ix, iy, i, j, xoff, yoff, max_off;
  double rot, x0, y0, midx, midy;
  float *data;
  float *data2;
  if (!arr || !ObsLoc)
    return(1);
  data = arr->data;
  nx = arr->axis[0];
  ny = arr->axis[1];
  x0 = ObsLoc->crpix1 - 1;
  y0 = ObsLoc->crpix2 - 1;
  midx = (nx-1.0)/2.0;
  midy = (ny-1.0)/2.0;
  if ((rot = fabs(ObsLoc->crota2)) > 179 && rot < 181)
    {
    // rotate image by 180 degrees by a flip flip
    float val;
    int half = ny / 2;
    int odd = ny & 1;
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
  xoff = round(x0 - midx);
  yoff = round(y0 - midy);
  max_off = 20.0 / ObsLoc->cdelt1;
  if (arr->parent_segment &&
      arr->parent_segment->record &&
      arr->parent_segment->record->seriesinfo && 
      arr->parent_segment->record->seriesinfo->seriesname && 
      strncasecmp(arr->parent_segment->record->seriesinfo->seriesname, "aia", 3) &&
      abs(xoff) < max_off && abs(yoff) < max_off) 
    {
    if (abs(xoff) >= 1 || abs(yoff) >= 1)
      {
      data2 = malloc(4*nx*ny);
      for (iy=0; iy<ny; iy++)
        {
        int jy = iy + yoff;
        for (ix=0; ix<nx; ix++)
          {
          int jx = ix + xoff;
          int idx = jy*nx + jx;
          int idx2 = iy*nx + ix;
          if (jx<0 || jx>=nx || jy<0 || jy>=ny)
            data2[idx2] = DRMS_MISSING_FLOAT;
          else
            data2[idx2] = data[idx];
          }
        }
      x0 -= xoff;
      y0 -= yoff;
      free(data);
      arr->data = data2;
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
    ObsLoc->crota2 = drms_getkey_double(rec, "CROTA2", &status); if (status) ObsLoc->crota2 = 0.0; // WCS default
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

/* GetMinObsInfo - gets minimum standard WCS keywords for e.g. heliographic mapped data */
ObsInfo_t *GetMinObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus)
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
    ObsLoc->crota2 = drms_getkey_double(rec, "CROTA2", &status); if (status) ObsLoc->crota2 = 0.0; // WCS default
    ObsLoc->sina = sin(ObsLoc->crota2*Deg2Rad);
    ObsLoc->cosa = sqrt (1.0 - ObsLoc->sina*ObsLoc->sina);
    ObsLoc->rsun_obs = drms_getkey_double(rec, "RSUN_OBS", &status);
    ObsLoc->obs_vr = drms_getkey_double(rec, "OBS_VR", &status);
    ObsLoc->obs_vw = drms_getkey_double(rec, "OBS_VW", &status);
    ObsLoc->obs_vn = drms_getkey_double(rec, "OBS_VN", &status);
    ObsLoc->obs_b0 = drms_getkey_double(rec, "CRLT_OBS", &status);
    ObsLoc->t_obs = t_obs;
    }
  *rstatus = 0;
  return(ObsLoc);
  }

// ----------------------------------------------------------------------

// In cases known to not have compact slotted series and cadence is specified
// generate explicit recordset list of closest good record to desired grid
// First get vector of times and quality
// Then if vector is not OK, quit.
// then: make temp file to hold recordset list
//       start with first time to define desired grid,
//       make array of desired times.
//       make empty array of recnums
//       search vector for good images nearest desired times
//       for each found time, write record query


const char *get_input_recset(DRMS_Env_t *drms_env, const char *inQuery)
  {
  static char newInQuery[102];
  TIME epoch = (cmdparams_exists(&cmdparams, "epoch")) ? params_get_time(&cmdparams, "epoch") : 0;
  DRMS_Array_t *data;
  TIME t_start, t_stop, t_now, t_want, t_diff, this_t_diff;
  int status = 1;
  int nrecs, irec;
  int nslots, islot;
  long long *recnums;
  TIME *t_this, half;
  TIME cadence;
  double *drecnum, *dquality;
  int quality;
  long long recnum;
  char keylist[DRMS_MAXQUERYLEN];
  static char filename[100];
  char *tmpdir;
  FILE *tmpfile;
  char newIn[DRMS_MAXQUERYLEN];
  char seriesname[DRMS_MAXQUERYLEN];
  char *lbracket;
  char *at = index(inQuery, '@');
  if (at && *at && (strncmp(inQuery,"aia.lev1[", 9)==0 ||
                    strncmp(inQuery,"hmi.lev1[", 9)==0 ||
                    strncmp(inQuery,"aia.lev1_nrt2[",14)==0 ||
                    strncmp(inQuery,"hmi.lev1_nrt[", 13)==0 ))
    {
    char *ip=(char *)inQuery, *op=newIn, *p;
    long n, mul;
    while ( *ip && ip<at )
      *op++ = *ip++;
    ip++; // skip the '@'
    n = strtol(ip, &p, 10); // get digits only
    if (*p == 's') mul = 1;
    else if (*p == 'm') mul = 60;
    else if (*p == 'h') mul = 3600;
    else if (*p == 'd') mul = 86400;
    else 
      {
      fprintf(stderr,"cant make sense of @xx cadence spec for aia or hmi lev1 data");
      return(NULL);
      }
    cadence = n * mul;
    ip = ++p;  // skip cadence multiplier
    while ( *ip )
      *op++ = *ip++;
    *op = '\0';
    half = cadence/2.0;
    sprintf(keylist, "T_OBS,QUALITY,recnum");
    data = drms_record_getvector(drms_env, newIn, keylist, DRMS_TYPE_DOUBLE, 0, &status);
    if (!data || status)
      {
      fprintf(stderr,"getkey_vector failed status=%d\n", status);
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
    strcpy(seriesname, inQuery);
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
    sprintf(newInQuery,"@%s", filename);
    return(newInQuery);
    }
  else
    return(inQuery);
  }
