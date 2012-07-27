#define CVSVERSION "$Id: hmi_limbdark.c,v 1.9 2012/07/27 23:07:57 phil Exp $"
/**
   @defgroup analysis
   @ingroup su_util

   @brief Compute solar limb darkening and remove.

   @par
   @code
   hmi_limbdark - emoves limb darkening.
  
   Input is expected to be hmi (or drms mdi) brightness image data with.
   Output is the residual image after removing limb-darkening with
   a polynomial fit.  The fit order defaults to 2.  The limb-darkening
   profile comes from Tech Note SOI_TN_095 but the values in use here are from
   an average of several runs with the -f flag set.
   Modified from MDI DSDS program fitlimbdark.
   Default is to rotate to solar north up, and center to nearest pixel.

   Fit formula comes from Pierce, A.K. and C. Slaughter, "Solar Limb Darkening", Solar Physics 51, 25-41, 1977.
  
   Parameters are:
     in         input recordset, expected to be an Ic product
     out        output seriesname
     f          FLAG: Compute limb-darkening fit parameters.
     l          FLAG: non-zero supresses limb darkening removal
     c          FLAG: Do not center or rotate from input data orientation
     r          FLAG: non-zero performs REVERSE limb darkening removal
     n          FLAG: normalize output by dividing by image mean
                      i.e. puts limb darkening back in
     croplimit   float crop limit for removing limb darkening, default 1.0
     requestid  optional string
  
   @endcode

   @par Synopsis:
   @code
   hmi_fitlimbdark  in=input data out=output data
   @endcode

*/

#include "jsoc.h"
#include "jsoc_main.h"
#include "fstats.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

char *module_name = "hmi_limbdark";

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

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

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *ObsLoc, int *status);

int  GetObsLocInfo(DRMS_Segment_t *seg, int i, int j, ObsInfo_t *ObsLoc);

int rm_limbdark(DRMS_Array_t *arr, DRMS_Array_t *outarr, ObsInfo_t *ObsLoc, double *coefs, int *ncrop, int do_reverse, int do_norm, double crop_limit2);

int fit_limbdark(DRMS_Array_t *arr, ObsInfo_t *ObsLoc, double* coefs);

int upNcenter(DRMS_Array_t *data, ObsInfo_t *ObsLoc);

// for lsqfit, also needs math.h
int lsqfitd(double y[],double x[],double a[],int n,int np);
double imprvd(double a[], double lu[], double x[], double b[], int n);
int ludcmpd(double a[], int n);
void bkslvd(double lu[], double b[], int n);
int matinvd(double a[], int n);
int matsold(double a[], double x[], double b[], int n);

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
     {ARG_STRING, "out", "NOT SPECIFIED",  "Output data series."},
     {ARG_FLAG, "l", "0", "disable limb darkening removal, implied by -r"},
     {ARG_FLAG, "r", "0", "Restore limb darkening"},
     {ARG_FLAG, "f", "0", "Fit limb darkening before applying"},
     {ARG_FLAG, "n", "0", "Normalize the final image by dividing by the mean"},
     {ARG_FLAG, "c", "0", "Supress center and flip 180 degrees if CROTA~180."},
     {ARG_FLAG, "x", "0", "Exclude pixels that deviate from the default LD, <0.875 or >1.25."},
     {ARG_FLAG, "h", "0", "Include full headers, set when requestid is present"},
     {ARG_DOUBLES, "coefs", "0.0", "Limb darkening coeficients, 5 needed"},
     {ARG_FLOAT, "croplimit", "1.0", "crop limit for removing limbdarkening"},
     {ARG_STRING, "requestid", "NA",  "JSOC export identifier"},
     {ARG_END}
};

int DoIt(void)
  {
  int noLD = cmdparams_isflagset(&cmdparams, "l");
  int noFlip = cmdparams_isflagset(&cmdparams, "c");
  int restoreLD = cmdparams_isflagset(&cmdparams, "r");
  int do_fit = cmdparams_isflagset(&cmdparams, "f");
  int do_normalize = cmdparams_isflagset(&cmdparams, "n");
  int do_exclude = cmdparams_isflagset(&cmdparams, "x");
  const char *inQuery = params_get_str(&cmdparams, "in");
  const char *outSeries = params_get_str(&cmdparams, "out");
  const char *requestid = params_get_str(&cmdparams, "requestid");
  int full_headers = cmdparams_isflagset(&cmdparams, "h") || strcmp(requestid, "NA");
  float crop_limit = params_get_float(&cmdparams, "croplimit");
  double crop_limit2 = crop_limit*crop_limit;
  char *p;
  int status = 0;
  ObsInfo_t *ObsLoc;
  // Coef version 1
  // static double defaultcoefs[] = {1.0, 0.443000, 0.139000, 0.041000, 0.012500, 0.001900};
  // char *CoefVersion = "1";
  // Coef version 2
  static double defaultcoefs[] = {1.0, 0.459224, 0.132395, 0.019601, 0.000802, -4.31934E-05 };
  char *CoefVersion = "2";

  double use_coefs[6];
  double n_user_coefs = cmdparams_get_int(&cmdparams, "coefs_nvals", &status);

  int nrecs, irec;
  DRMS_RecordSet_t *inRS, *outRS;

  if (strcmp(inQuery, "NOT SPECIFIED") == 0)
    DIE("Must have input series");;
  if (strcmp(outSeries, "NOT SPECIFIED") == 0)
    DIE("Must have output series");;

  printf("FitLimbDark\n");
  if (n_user_coefs == 6)
    {
    double *cmdcoefs;
    int i;
    CoefVersion = "user given";
    cmdparams_get_dblarr(&cmdparams, "coefs", &cmdcoefs, &status);
    for (i=0; i<6; i++)
      {
      use_coefs[i] = cmdcoefs[i];
      printf(" Coef%d = %0.6f\n", i, use_coefs[i]);
      }
    }
  else
    {
    int i;
    for (i=0; i<6; i++)
      use_coefs[i] = defaultcoefs[i];
    printf(" Use default coefs\n");
    }
  
  if (restoreLD) noLD = 1;
  if (noLD)
	printf("   Supress limb darkening removal\n");

  inRS = drms_open_records(drms_env, inQuery, &status);
  if (status || !inRS || (nrecs = inRS->n) == 0)
    DIE("No input records found");

  outRS = drms_create_records(drms_env, nrecs, (char *)outSeries, DRMS_PERMANENT, &status);
  if (status || !outRS)
    DIE("Failed to create output recordset");

  for (irec=0; irec<nrecs; irec++)
    {
    DRMS_Record_t *inRec = inRS->records[irec];
    DRMS_Record_t *outRec;
    int quality = drms_getkey_int(inRec, "QUALITY", NULL);
    ObsInfo_t *ObsLoc;

    outRec = outRS->records[irec];
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);
    drms_copykey(outRec, inRec, "T_REC");
    drms_copykey(outRec, inRec, "T_OBS");
    drms_copykey(outRec, inRec, "QUALITY");
    drms_setkey_time(outRec, "DATE", time(0) + UNIX_EPOCH);
    drms_setkey_string(outRec, "RequestID", requestid);
    drms_setkey_string(outRec, "CODEVER4", CVSVERSION);

    if (quality >= 0)
      {
      double mean=1.0;
      double coefs[6];
      DRMS_Segment_t *inSeg, *outSeg;
      DRMS_Array_t *inArray, *outArray;
      int ncropped = 0;
      inSeg = drms_segment_lookupnum(inRec, 0);
      inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
      if (status)
        {
        printf(" No data file found but QUALITY not bad, status=%d\n", status);
        drms_free_array(inArray);
        continue;
        }

      ObsLoc = GetObsInfo(inSeg, NULL, &status);
      if (status)
        DIE("Failed to get observatory location.");

      // if (do_exclude && do_fit)
      if (do_fit)
        {
        if (do_exclude )
          {
          DRMS_Array_t *xArray = drms_array_create(DRMS_TYPE_FLOAT, inArray->naxis, inArray->axis, NULL, &status);
          DRMS_Array_t *inxArray = drms_array_create(DRMS_TYPE_FLOAT, inArray->naxis, inArray->axis, NULL, &status);
          float *iv = (float*)inArray->data;
          float *v = (float*)xArray->data;
          float *nv = (float*)inxArray->data;
          int i, n = inArray->axis[0] * inArray->axis[1];
          int nskip = 0;
          float missval = DRMS_MISSING_FLOAT;
          rm_limbdark(inArray, xArray, ObsLoc, use_coefs, &ncropped, 0, 1, crop_limit2);
          for (i=0; i<n; i++, v++, iv++)
            if (!isnan(*v) && *v > 0.850 && *v < 1.150) 
              *nv++ = *iv;
            else
              {
              *nv++ = missval;
              if (!isnan(*v))
                nskip++;
              }
          printf("exclude %d pixels\n", nskip);
          fit_limbdark(inxArray, ObsLoc, coefs);
          drms_free_array(xArray);
          drms_free_array(inxArray);
          }
        else
          fit_limbdark(inArray, ObsLoc, coefs);
        }

      outSeg = drms_segment_lookupnum(outRec, 0);
      outArray = drms_array_create(DRMS_TYPE_FLOAT, inArray->naxis, inArray->axis, NULL, &status);

      if (rm_limbdark(inArray, outArray, ObsLoc, (do_fit ? coefs : use_coefs), &ncropped, restoreLD, do_normalize, crop_limit2) == DRMS_SUCCESS)
        {
        int totvals, datavals;

        if (drms_keyword_lookup(outRec, "COEF_VER", 0))
          {
          drms_setkey_string(outRec, "COEF_VER", CoefVersion);
          drms_setkey_float(outRec, "LDCoef0", (float)(do_fit ? coefs[0] : use_coefs[0]));
          drms_setkey_float(outRec, "LDCoef1", (float)(do_fit ? coefs[1] : use_coefs[1]));
          drms_setkey_float(outRec, "LDCoef2", (float)(do_fit ? coefs[2] : use_coefs[2]));
          drms_setkey_float(outRec, "LDCoef3", (float)(do_fit ? coefs[3] : use_coefs[3]));
          drms_setkey_float(outRec, "LDCoef4", (float)(do_fit ? coefs[4] : use_coefs[4]));
          drms_setkey_float(outRec, "LDCoef5", (float)(do_fit ? coefs[5] : use_coefs[5]));
          }
        else
          {
          char LDhistory[1000];
          sprintf(LDhistory, "hmi_limbdark: COEF_VER=%s", CoefVersion);
          drms_appendhistory(outRec, LDhistory, 0);
          sprintf(LDhistory, "LDCoefs=%0.6f,%0.6f,%0.6f,%0.6f,%0.6f,%0.6f", 
            (float)(do_fit ? coefs[0] : use_coefs[0]),
            (float)(do_fit ? coefs[1] : use_coefs[1]),
            (float)(do_fit ? coefs[2] : use_coefs[2]),
            (float)(do_fit ? coefs[3] : use_coefs[3]),
            (float)(do_fit ? coefs[4] : use_coefs[4]),
            (float)(do_fit ? coefs[5] : use_coefs[5]));
          drms_appendhistory(outRec, LDhistory, 1);
          }
        if (!noFlip)
          upNcenter(outArray, ObsLoc);
        drms_setkey_double(outRec, "CRPIX1", ObsLoc->crpix1);
        drms_setkey_double(outRec, "CRPIX2", ObsLoc->crpix2);
        drms_setkey_double(outRec, "CROTA2", ObsLoc->crota2);

        if (do_normalize)
          {
          outArray->bzero = 1.0;
          outArray->bscale= 1.0/30000.0;
          }
       else
          {
          outArray->bzero = 32768.0;
          outArray->bscale=1.5;
          }

        set_statistics(outSeg, outArray, 1);
        totvals = drms_getkey_int(outRec, "TOTVALS", &status);
        if (!status)
          {
          datavals = drms_getkey_int(outRec, "DATAVALS", &status);
          totvals -= ncropped;
          drms_setkey_int(outRec, "TOTVALS", totvals);
          drms_setkey_int(outRec, "MISSVALS", totvals - datavals);
          }

        if (full_headers)
          drms_segment_writewithkeys(outSeg, outArray, 0);
        else
          drms_segment_write(outSeg, outArray, 0);
        drms_free_array(inArray);
        drms_free_array(outArray);
        }
      // printf("Done with %s\n", drms_getkey_string(inRec, "T_OBS", NULL));
      }
    else
      printf("Skip Rec %d, Quality=%#08x\n", irec, quality);
    }
  drms_close_records(outRS, DRMS_INSERT_RECORD);
  drms_close_records(inRS, DRMS_FREE_RECORD);
  return (DRMS_SUCCESS);
  } // end of DoIt

int fit_limbdark(DRMS_Array_t *arr, ObsInfo_t *ObsLoc, double* coefs)
  {
  int status = 0;
  double x0 = ObsLoc->crpix1 - 1;
  double y0 = ObsLoc->crpix2 - 1;
  double rsun = ObsLoc->rsun_obs/ObsLoc->cdelt1;
  
  int ix, iy;
  int nx = arr->axis[0];
  int ny = arr->axis[1];
  double scale;
  double crop_limit2;
  float *data = (float*)arr->data;
  float missval = DRMS_MISSING_FLOAT;
  int n;
  int ord;
  double *f = (double *)malloc(nx*ny*sizeof(double));
  double *c = (double *)malloc(6*nx*ny*sizeof(double));
  double fitcoefs[6];
  
  if (!f || !c) DIE("malloc problem");

  scale = 1.0/rsun;
  crop_limit2 = 0.99975;

  n = 0;
  for (iy=0; iy<ny; iy++)
    for (ix=0; ix<nx; ix++)
        {
	double costheta2;
        double xi, mu, z, ld;
	double x, y, R2;
        float *Ip = data + iy*nx + ix;

	if (drms_ismissing_float(*Ip))
	  continue;

        /* get coordinates of point */
        x = ((double)ix - x0) * scale; /* x,y in pixel coords */
        y = ((double)iy - y0) * scale;

        /* only fit points within limit radius */
        R2 = x*x + y*y;

        if (R2 >= crop_limit2)
	  continue;

        costheta2 = 1.0 - R2;
        mu = sqrt(costheta2);
	xi = log(mu);
	z = 1.0;
        f[n] = *Ip;
        c[6*n + 0] = 1.0;
	for (ord=1; ord<6; ord++)
	  {
	  z *= xi;
          c[6*n + ord] = z;
	  }
        n++;
	}
    if (!lsqfitd(f, c, fitcoefs, 5, n))
      {
      fprintf(stderr,"lsqfit failure\n");
      fitcoefs[0] = fitcoefs[1] = fitcoefs[2] = fitcoefs[3] = fitcoefs[4] = fitcoefs[5] =  DRMS_MISSING_DOUBLE;
      }
    printf("Fit %d points, Coefs = %0.6f", n, fitcoefs[0]);
    for (ord=0; ord<6; ord++)
      {
      coefs[ord] = fitcoefs[ord]/fitcoefs[0];
      printf(", %8.6f", coefs[ord]);
      }
    printf("\n");
    coefs[0] = fitcoefs[0];
  free(f);
  free(c);
  return(0);
  }

int rm_limbdark(DRMS_Array_t *arr, DRMS_Array_t *outarr, ObsInfo_t *ObsLoc, double *coefs, int *ncrop, int do_reverse, int do_norm, double crop_limit2)
  {
  double x0 = ObsLoc->crpix1 - 1;
  double y0 = ObsLoc->crpix2 - 1;
  double rsun = ObsLoc->rsun_obs/ObsLoc->cdelt1;
  
  int ix, iy;
  int nx = arr->axis[0];
  int ny = arr->axis[1];
  double scale;
  float *data = (float *)arr->data;
  float *odata = (float *)outarr->data;
  float missval = DRMS_MISSING_FLOAT;
  int ncropped = 0;
  
  scale = 1.0/rsun;

  for (iy=0; iy<ny; iy++)
    for (ix=0; ix<nx; ix++)
        {
	double costheta2;
        double xi, mu, z, ld;
	double x, y, R2;
        int ord;
        float *Ip = data + iy*nx + ix;
        float *Op = odata + iy*nx + ix;

	if (drms_ismissing_float(*Ip))
          {
          *Op = missval;
	  continue;
          }

        /* get coordinates of point */
        x = ((double)ix - x0) * scale; /* x,y in pixel coords */
        y = ((double)iy - y0) * scale;

        /* only fix points within limit radius */
        R2 = x*x + y*y;

        if (R2 >= crop_limit2)
	  {
	  *Op = missval;
	  ncropped++;
	  continue;
	  }

        if (R2 < 0.99975)
          costheta2 = 1.0 - R2;
        else
          costheta2 = 1.0 - 0.99975;  /* 1/4 pixel from limb */

        mu = sqrt(costheta2);
	xi = log(mu);
	z = 1.0;
	ld = 1.0;
	for (ord=1; ord<6; ord++)
	  {
	  z *= xi;
	  ld += coefs[ord] * z;
	  }
        if (ld <= 0.0)
          {
          *Op = missval;
          ncropped++;
          }
	else if (do_reverse)
	  *Op = *Ip * ld;
	else
	  *Op = *Ip / ld;
	}
  *ncrop = ncropped;
  if (do_norm)
    {
    double mean;
    float *v = odata;
    int i, n = outarr->axis[0] * outarr->axis[1];
    if (coefs[0] == 1.0)
      {
      int nok=0;
      v = odata;
      for (i=0; i<n; i++, v++)
        if (!isnan(*v)) 
          { nok++; mean += *v; }
      if (nok)
        mean /= nok;
      else
        mean = 1.0;
      }
    else
      mean = coefs[0];
    v = odata;
    for (i=0; i<n; i++, v++)
      if (!isnan(*v)) 
        *v /= mean;
    }
  return(0);
  }

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

/* center whith whole pixel shifts and rotate by 180 if needed */
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc)
  {
  int nx, ny, ix, iy, i, j, xoff, yoff;
  double rot, x0, y0, mid;
  float *data;
  if (!arr || !ObsLoc)
    return(1);
  data = (float *)arr->data;
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
    x0 = nx - x0;
    y0 = ny - y0;
    rot = ObsLoc->crota2 - 180.0;
    if (rot < -90.0) rot += 360.0;
    ObsLoc->crota2 = rot;
    }
  xoff = round(x0 - mid);
  yoff = round(y0 - mid);
  if (abs(xoff) > 1.0)
    {
    for (iy=0; iy<ny; iy++)
      {
      float valarr[nx];
      for (ix=0; ix<nx; ix++)
        {
        int jx = ix - xoff;
        if (jx >= nx) jx -= nx;
        if (jx < 0) jx += nx;
        valarr[jx] = data[iy*nx + ix];
        }
      for (ix=0; ix<nx; ix++)
        data[iy*nx + ix] = valarr[ix];
      }
    x0 -= xoff;
    }
  if (abs(yoff) > 1.0)
    {
    for (ix=0; ix<nx; ix++)
      {
      float valarr[ny];
      for (iy=0; iy<ny; iy++)
        {
        int jy = iy - yoff;
        if (jy >= ny) jy -= ny;
        if (jy < 0) jy += ny;
        valarr[jy] = data[iy*nx + ix];
        }
      for (iy=0; iy<ny; iy++)
        data[iy*nx + ix] = valarr[iy];
      }
    y0 -= yoff;
    }
  ObsLoc->crpix1 = x0 + 1;
  ObsLoc->crpix2 = y0 + 1;
  return(0);
  }



#define ITMAX 7		/* max imprvd iterations */
#define EPS 2.8e-17	/* relative precision of arithmetic */

/*
 * Solve the linear system a x = b for x, where b and x are
 * n dimensional vectors and a is a n x n matrix of doubles,
 * uses ludcmpd (lower/upper triangular decomposition), bkslvd
 * (back solve lu system), imprvd (refine solution from ludcmpd
 * and bkslvd - calls bkslvd itself)
 */
int matsold(double a[], double x[], double b[], int n)
{
	double *lu, err, wrk, wrk1;
	extern double imprvd();
	int k;

	/*
	 * allocate storage for lu decomposition of a
	 */

	lu = (double *)malloc( (n*n*8) );

	/*
	 * copy a to lu
	 */

	{
	double *src = &a[0];
	double *dest = lu;
	int i = (n*n);

	do
		{
		*(dest++) = *(src++);
		} while(--i);
	}

	/*
	 * get lu decomposition in lu
	 */

	if(ludcmpd(lu, n) == 0)
		{	/* singular matrix - give up */
		free(lu);
		return(0);
		}

	/*
	 * calculate the solution to the linear system
	 * using a maximum of ITMAX iterations of imprvd
	 */

	{
	int i;
	double *dest = x;
	double *src = b;

	/*
	 * copy b to x
	 */

	i = n;
	do
		{
		*(dest++) = *(src++);
		}while(--i);

	/*
	 * backsolve and improve
	 */

	bkslvd(lu, x, n);

	for(k=0;k<ITMAX;k++)
		{
		err = imprvd(a,lu,x,b,n);
		wrk = 0.0;
		src = x;
		i = n;
		do
			{
			wrk1 = (*src < 0.0) ? -(*src++) : *(src++);
			wrk = (wrk < wrk1) ? wrk1 : wrk;
			}while(--i);
		wrk1 = err/wrk;
		if(wrk1 < EPS)break;
		}

	}

	/*
	 * free temporary storage
	 */

	free(lu);
	return(1);
}

/*
 * invert the n x n matrix of double precision numbers A,
 * uses ludcmpd (lower/upper triangular decomposition), bkslvd
 * (back solve lu system), imprvd (refine solution from ludcmpd
 * and bkslvd - calls bkslvd itself)
 */
int matinvd(double a[], int n)
{
	double *ainv, *lu, *x, *b, err, wrk, wrk1;
	extern double imprvd();
	int j, k;

	/*
	 * allocate storage for intermediate arrays and
	 * vectors
	 */

	ainv = (double *)malloc( (n*n*8) );
	lu = (double *)malloc( (n*n*8) );
	x = (double *)malloc( (n*8) );
	b = (double *)malloc( (n*8) );

	/*
	 * copy a to lu
	 */

	{
	double *src = &a[0];
	double *dest = lu;
	int i = (n*n);

	do
		{
		*(dest++) = *(src++);
		} while(--i);
	}

	/*
	 * get lu decomposition in lu
	 */

	if(ludcmpd(lu, n) == 0)
		{	/* singular matrix - give up */
		free(ainv);
		free(lu);
		free(x);
		free(b);
		return(0);
		}

	/*
	 * calculate the columns of the inverse one at a time
	 * using a maximum of ITMAX iterations of imprvd
	 */

	for(j = 0;j<n;j++)
		{
		int i;
		double *dest = b;
		double *src = x;

		/*
		 * construct correct right hand sides for calculating
		 * the current column of inverse
		 */

		i = n;
		do
			{
			*(dest++) = 0.0;
			*(src++) = 0.0;
			}while(--i);
		b[j] = x[j] = 1.0;

		/*
		 * find column of inverse
		 */

		bkslvd(lu, x, n);

		for(k=0;k<ITMAX;k++)
			{
			err = imprvd(a,lu,x,b,n);
			wrk = 0.0;
			src = x;
			i = n;
			do
				{
				wrk1 = (*src < 0.0) ? -(*src++) : *(src++);
				wrk = (wrk < wrk1) ? wrk1 : wrk;
				}while(--i);
			wrk1 = err/wrk;
			if(wrk1 < EPS)break;
			}

		/*
		 * save column of inverse in ainv
		 */

		src = x;
		dest = &ainv[j];
		i = n;
		do
			{
			*dest = *(src++);
			dest += n;
			}while(--i);
		}
	{
	int i = n*n;
	double *src = ainv;
	double *dest = &a[0];

	/*
	 * copy inverse onto original a
	 */

	do
		{
		*(dest++) = *(src++);
		}while(--i);
	}

	/*
	 * free temporary storage
	 */

	free(ainv);
	free(lu);
	free(x);
	free(b);
	return(1);
}

/*
 * back solve for x in Ax = b, A an n x n matrix double prescison 
 * numbers and x and b are n dimensional vectors - given the LU
 * decomposition of A.  b is overwritten with x.
 */
void bkslvd(double lu[], double b[], int n)
{
	int i, k;

	/*
	 * multiply by lower triangular
	 */

	for(i=0;i<(n-1);i++)
		{
		int j = n-1-i;
		double *lower = &lu[i*(n+1)+n];
		double *current = &b[i+1];
		double first = b[i];

		do
			{

			*(current++) -= first*(*lower);
			lower += n;
			} while(--j);
		}

	/*
	 * backsolve upper triangular
	 */

	for(i = n-1;i >= 0;i--)
		{
		double *factor = &lu[n*i + i + 1];
		double *vector = &b[i+1];
		double *current = &b[i];

		for(k=i;k<n-1;k++)
			{
			*current -= *(factor++)*(*(vector++));
			}
		*current /= lu[n*i+i];
		}
	return;
}

/*
 * LU decomposition of an n x n matrix of double precision
 * numbers by Gaussian elimination in place - original
 * matrix is overwritten - no pivoting
 *
 * clearer but slower version is in comments at the end - present
 * version is optimized for PDP 11/45
 */
int ludcmpd(double a[], int n)
{
	double wrk;	/* to save divides */
	int i, j;

	for(i=0;i<(n-1);i++)
		{
		if(a[i+n*i] != 0.0)wrk = 1.0/a[i+n*i];
			else return(0);	/* zero divide - give up */
		for(j=i+1;j<n;j++)
			{
			double *first = &a[i+1+n*i];
			double *row = &a[i+1+n*j];
			int k = n-(i+1);
			double elim = a[i+n*j]*wrk;

			/* save lower triangular */
			a[i+n*j] = elim;
			/*
			 * generate upper triangular - optimized
			 * for fast innermost loop
			 *
			 * note that initialization time
			 * may dominate for small n - but
			 * then it doesn't take long anyway
			 */
			do
				{
				*(row++) -= *(first++)*elim;
				}while(--k);
			}
		}
	return(1);
}

/*
 * improve the solution, x, to the linear equation Ax = b, where A is a
 * n x n matrix of double precision numbers, x is a provisional
 * solution to the equation (most likely from bkslv), b is the
 * original right hand side and lu is the decomposition of A into
 * a product of a lower triangular and an upper triangular matrix
 * (most likely from ludcmp). Returns the maximum absolute error
 * of the error vector (Ax - b).  The improved x overwrites the old x.
 */
double imprvd(double a[], double lu[], double x[], double b[], int n)
{
	double *berr, err;
	int i;

	berr = (double *)malloc( (n*8) );
	err = 0.0;
	/*
	 * calculate error vector and max absolute error
	 */

	for(i=0;i<n;i++)
		{
		int j = n;
		double *vector = &x[0];
		double *row = &a[n*i];
		double accum = 0.0;

		do
			{
			accum += (*(vector++))*(*(row++));
			}while(--j);
		berr[i] = b[i] - accum;

		accum = (berr[i] < 0) ? -berr[i] : berr[i];
		err = (err < accum) ? accum : err;
		}

	bkslvd(lu, berr, n);
	for(i=0;i<n;i++)x[i] += berr[i];
	free(berr);
	return(err);
}

/*
 * linear least squares fit to a vector y(i) of np vectors x(i,j), where
 * i runs from 1 to np, and j runs from 0 to n, the coeficients are returned
 * in the vector a(j). The weight of each y(i) is contained in  x(i,0) -
 * patterned on the algol routine written by Leif Svalgaard in algol.
 * The usual "normal equations" for a least square fit are solved, with
 * the following modifications:
 *
 *	1)  The sum of the weights is normalized to np, the number of data
 *	points, forming adjusted weights, w(i).
 *
 *	2)  The y(i) are normalized by subtracting the weighted mean,
 *	ym = (sum of y(i)*w(i))/np, and dividing by an estimate of the
 *	square root of the variance sqrt( (sum (y(i)-ym)**2)/(np-1)).
 *
 *	3)  The fitting vectors x(i,j) are also normalized by subttracting
 *	the weighted mean and dividing by the square root of the estimated
 *	variance.
 *
 *	4)  Both left and right sides of the linear equations are divided
 *	by (np-1) so that the expected values of the elements of the matrix
 *	and the right hand side are 1.
 *
 *	5)  The resulting solutions to the linear equations must be transformed
 *	to account for the above normalizations before being returned.
 *
 * Returns 1 on succesful completion, 0 on some error condition. From
 * C, &y[0] should point to np doubles, &x[0] should point to np*(n+1)
 * doubles and &a[0] should point to (n+1) doubles.
 */
int lsqfitd(double y[],double x[],double a[],int n,int np)
{
	double *ar, *lu, *sx, *xm, *r;
	int i, k;
	double f1, ym, sum, s, wmi;
	extern double sqrt(), imprvd();
	/*
	 * first allocate and zero working arrays and zero sums.
	 *
	 * ym is the weighted mean of the y(i)'s
	 *
	 * s is the square root of the estimated variance of the y(i)'s
	 *
	 * wmi is the one over the sum of the weights
	 *
	 * xm(j) is the vector of the weighted means of the x(i,j)'s
	 *
	 * sx(j) is the vector of the square roots of the estimated
	 * variances of the x(i,j)'s
	 *
	 */

	if (np < 3)
		return(0);
	ar = (double *)malloc( n*n*8 );
	lu = (double *)malloc( n*n*8 );
	sx = (double *)malloc( n*8 );
	xm = (double *)malloc( n*8 );
	r = (double *)malloc( n*8 );

	sum = ym = s = 0.0;
	/* zero r and ar */
	for(i=0;i<n;i++)
		{
		double *dp = &ar[i*n];
		int j = n;

		r[i] = 0.0;
		do
			{
			*(dp)++ = 0.0;
			}while(--j);
		}
	{
	/* zero sx and xm */
	double *dp1 = &sx[0];
	double *dp2 = &xm[0];
	int j = n;

	do
		{
		*(dp1++) = 0.0;
		*(dp2++) = 0.0;
		}while(--j);
	}

	/* accumulate sum of weights, weighted sum of y's and
	 * fitting vectors
	 */
	for(i=0;i<np;i++)
		{
		double *dp1 = &xm[0];
		double *dp2 = &x[(n+1)*i];
		int j = n;
		double wi;

		wi = *(dp2++);
		sum += wi;
		ym += wi*y[i];
		do
			{
			*(dp1++) += wi*(*(dp2++));
			}while(--j);
		}
	/* if the sum of the weights is 0, return 0 */
	if(sum == 0.0)
		{	
		free(ar);
		free(lu);
		free(sx);
		free(xm);
		free(r);
		return(0);
		}

	/* divide to get weighted averages from weighted sums */
	ym /= sum;
	wmi = ((double)np)/sum;
	for(i=0;i<n;i++)xm[i] /= sum;

	/* accumulate weighted square of differences from weighted means, 
	 * matrix elements and right hand side of linear equation
	 */
	for(i=0;i<np;i++)
		{
		double fy, wi;

		wi = x[i*(n+1)]*wmi;
		fy = y[i] - ym;
		s += wi*fy*fy;
		for(k=0;k<n;k++)
			{
			double *dp1 = &ar[k];
			double *dp2 = &ar[n*k];
			double *dp3 = &xm[0];
			double *dp4 = &x[i*(n+1)+1];
			int j;
			double fx;

			fx = *(dp4+k) - *(dp3+k);
			sx[k]+= wi*fx*fx;
			r[k] += wi*fx*fy;

			for(j=0;j<=k;j++)
			{
				*dp2 += wi*fx*(*(dp4++) - *(dp3++));
				*dp1 = *(dp2++);
				dp1 += n;
				}
			}
		}

	/* normalize equations so all matrix elements have
	 * expectation value 1
	 */

	f1 = np - 1.0;
	s = sqrt(s/f1);
	for(i=0;i<n;i++)
		{
		double *dp1 = &ar[i];
		double *dp2 = &ar[i*n];
		double *dp3 = &sx[i];
		int j;
		double sxj;

		*dp3 = sqrt(*dp3/f1);
		sxj = *dp3;
		r[i] /= f1*(*dp3)*s;
		dp3 = &sx[0];

		for(j=0;j<=i;j++)
			{
			if (*dp3 == 0.0)
				return(0);
			*dp2 /= f1*(*(dp3++))*sxj;
			*dp1 = *(dp2++);
			dp1 += n;
			}
		}

	/*
	 * solve linear system - using ludcmpd, bkslvd, and imprvd
	 * first make a copy of ar in lu 
	 */

	{
	double *dp1 = &ar[0];
	double *dp2 = &lu[0];
	int j = n*n;

	do
		{
		*(dp2++) = *(dp1++);
		}while(--j);
	}

	if(ludcmpd(lu, n) == 0)
		{	/* singular matrix - give up */	
		free(ar);
		free(lu);
		free(sx);
		free(xm);
		free(r);
		return(0);
		}

	/* copy r to a offset for the constant term */

	{
	double *dp1 = &a[1];
	double *dp2 = &r[0];
	int j = n;

	do
		{
		*(dp1++) = *(dp2++);
		}while(--j);
	}

	/* back solve equations */

	bkslvd(lu, &a[1], n);

	/* do a maximum of ITMAX imprvd refinements */

	for(i=0;i<ITMAX;i++)
		if(imprvd(ar, lu, &a[1], r, n) < EPS)break;

	/* transform calculated fits back to original variables */

	{
	double *dp1 = &a[0];
	double *dp2 = &sx[0];
	double *dp3 = &xm[0];
	int j;

	*(dp1++) = ym;
	for(j=0;j<n;j++)
		{
		*dp1 *= s;
		*dp1 /= *(dp2++);
		a[0] -= (*(dp1++))*(*(dp3++));
		}
	}

	/* free allocated storage */
	free(ar);
	free(lu);
	free(sx);
	free(xm);
	free(r);

	return(1);
}

/* Algol routine after which this routine is patterned 

.
external
procedure lsqfit(y,x,a,n,np); value n,np;
real array y,x,a; integer n,np;
begin
	comment:	Computes the least-squares fit to
				y(i) = a(0) + a(1)*x(i,1) + ... + a(n)*x(i,n)
				the weight for y(i) is stored in x(i,0).
				The number of points is 'np'    ;
integer yl,yh,yn,al,ah;
yl:= 1; yh:= np; yn:= yh-yl+1; ah:= n; al:= 1;

begin real array ar(1:ah,1:ah), xm,r,sx(1:ah);
real s,sum,ym,f1,wm,fx,fy,wi,sxj,aj;
integer i,j,k;

sum:= ym:= s:= 0.0;
for j:= 1 step 1 until ah do
begin xm(j):= sx(j):= a(j):= r(j):= 0;
 for k:= 1 step 1 until ah do ar(j,k):= 0;
end; 

for i:= yl step 1 until yh do
begin wi:= x(i,0);
 sum:= sum + wi; ym:= ym + wi*y(i);
 for j:= 1 step 1 until ah do xm(j):= xm(j) + wi*x(i,j);
end;

if sum=0 then a(0):= missing
else
begin
ym:= ym/sum; wm:= sum/yn;
for j:= 1 step 1 until ah do xm(j):= xm(j)/sum;

for i:= yl step 1 until yh do
begin wi:= x(i,0)/wm; fy:= y(i) - ym;
 s:= s + wi*fy**2;
 for j:= 1 step 1 until ah do
 begin fx:= x(i,j) - xm(j);
  sx(j):= sx(j) + wi*fx**2; r(j):= r(j) + wi*fx*fy;
  for k:= 1 step 1 until j do
  ar(j,k):= ar(j,k) + wi*fx*(x(i,k) - xm(k));
 end j;
end i;

f1:= yn - 1; s:= sqrt(s/f1);
for j:= 1 step 1 until ah do
begin sxj:= sx(j):= sqrt(sx(j)/f1); r(j):= r(j)/(f1*sxj*s);
 for k:= 1 step 1 until j do
 begin ar(j,k):= ar(j,k)/(f1*sxj*sx(k));
  ar(k,j):= ar(j,k);
 end;
end;

if syminv(ar) = 0 then a(0):= missing else
begin comment: normal case;   a(0):= ym;
 for j:= 1 step 1 until ah do
 begin for k:= 1 step 1 until ah do a(j):= a(j)+r(k)*ar(j,k);
  aj:= a(j):= a(j)*s/sx(j); a(0):= a(0) - aj*xm(j);
 end;
end;

end;
end;

end;
end
*/
