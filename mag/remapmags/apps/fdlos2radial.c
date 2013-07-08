/*
* This jsoc module is to convert los magnetic field to radial field
* by assuming that the field is purely radial.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "jsoc_main.h"
#include "astro.h"
#include "drms_dsdsapi.h"

#define PI              (M_PI)
#define RADSINDEG       (PI/180)

//#define MIN(a,b) (((a)<(b)) ? (a) : (b))
//#define MAX(a,b) (((a)>(b)) ? (a) : (b))

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}

/* cmd-line parameters */
#define kRecSetIn       "in"
#define kSeriesOut      "out"
#define kSegIn          "segin"
#define kSegOut         "segout"
#define kNOTSPECIFIED   "not specified"

#define RADSINDEG       (PI/180)
#define RAD2ARCSEC      (648000. / M_PI)

char *module_name = "fdlos2radial";

ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kRecSetIn, "", "Input data records."},
   {ARG_STRING,  kSeriesOut, "", "Output data series."},
   {ARG_STRING,  kSegIn, kNOTSPECIFIED, ""},
   {ARG_STRING,  kSegOut, kNOTSPECIFIED, ""},
   {ARG_INT, "XSIZE", "4096", "", ""},
   {ARG_INT, "YSIZE", "4096", "", ""},
   {ARG_END}
};

int DoIt(void)
{

    int status = DRMS_SUCCESS;
    int ds, nds;
    int xsize, ysize;
    char *inRecQuery, *outRecQuery;
    char tstr[64];
    int y,m,d;
    DRMS_RecordSet_t *inRD, *outRD;

    inRecQuery = (char *)params_get_str(&cmdparams, kRecSetIn);
    outRecQuery = (char *)params_get_str(&cmdparams, kSeriesOut);
    xsize = cmdparams_get_int(&cmdparams, "XSIZE", &status);
    ysize = cmdparams_get_int(&cmdparams, "YSIZE", &status);
    inRD = drms_open_records(drms_env, inRecQuery, &status);
    if (status || inRD->n == 0)
        DIE("No input data found");
    nds = inRD->n;

    outRD = drms_create_records(drms_env, nds, outRecQuery, DRMS_PERMANENT, &status);
    if (status)
        DIE("Output recordset not created");

    for (ds = 0; ds < nds; ds++)
      {
        int xDim = xsize, yDim = ysize;
        double xDist = 0.0;
        double yDist = 0.0;
        double xDist2 = 0.0;
        double yDist2 = 0.0;
        double xDist2PlusyDist2 = 0.0;
        double dToDiskCenter = 0.0;
        int yOff = 0;
        int iData = 0;
        float radialCorrected = 0;

        double sinrho = 0;
        double cosrho = 0;
        float dSun, rSun_ref, cdelt, asd;
        float rSun, xCenter, yCenter;
        float *bRadial;
        float *bLos;
        DRMS_Segment_t *inSeg = NULL;
        DRMS_Array_t *inArray;
        DRMS_Record_t *inRec;

        inRec = inRD->records[ds];

        dSun = drms_getkey_double(inRec, "DSUN_OBS", &status);
        rSun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
          if (status) rSun_ref = 6.96e8;
        cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
        asd = asin(rSun_ref/dSun);
        rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
        xCenter = drms_getkey_float(inRec, "CRPIX1", &status);
        yCenter = drms_getkey_float(inRec, "CRPIX2", &status);
        bRadial = (float *) malloc(xDim * yDim * sizeof(float));

        inSeg = drms_segment_lookupnum(inRec, 0);
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status)
         {
          printf("Bad data, skip! \n");
          continue;
         }
        bLos = (float *)inArray->data;
      
        int jy = 0;
        for (jy = 0; jy < yDim; jy++)
          { 
              int ix = 0;
              yDist = (double)jy - yCenter;
              yDist2 = yDist * yDist;
              yOff = jy * xDim;

              for (ix = 0; ix < xDim; ix++)
                {
                   iData = yOff + ix;
                   xDist = (double)ix - xCenter;
                   xDist2 = xDist * xDist;
                   xDist2PlusyDist2 = xDist2 + yDist2;
                   dToDiskCenter = sqrt(xDist2PlusyDist2);
                   if (dToDiskCenter >= rSun)
                     {
                         bRadial[iData] = DRMS_MISSING_FLOAT;
                         continue;
                     }

                   if (isnan(bLos[iData]))
                     {
                        bRadial[iData] = DRMS_MISSING_FLOAT;
                        continue;
                     }

                   sinrho = dToDiskCenter / rSun;
                   cosrho = sqrt(1 - sinrho * sinrho);
                   if (cosrho > 0)
                     {
                       bRadial[iData] =  (float)bLos[iData] / cosrho;
//                         bRadial[iData] =  (float)bLos[iData];
                     }
                   else
                     {
                       bRadial[iData] = DRMS_MISSING_FLOAT;
                     }
                }
          }

        DRMS_Record_t *outRec;
        DRMS_Segment_t *outSeg;
        DRMS_Array_t *outArray;
        int outDims[2] = {xDim, yDim};
        outRec = outRD->records[ds];

        drms_copykey(outRec, inRec, "T_REC");
        sprint_ut(tstr, CURRENT_SYSTEM_TIME);
        sscanf(tstr, "%d.%d.%d", &y, &m, &d);
        sprintf(tstr, "%04d-%02d-%02d", y, m, d);
        drms_setkey_string(outRec, "DATE", tstr);

//        drms_copykey(outRec, inRec, "DATE");
        drms_copykey(outRec, inRec, "DATE__OBS");
        drms_copykey(outRec, inRec, "INSTRUME");
        drms_copykey(outRec, inRec, "CAMERA");
        drms_copykey(outRec, inRec, "HISTORY");
        drms_copykey(outRec, inRec, "COMMENT");
        drms_copykey(outRec, inRec, "BLD_VERS");

        drms_copykey(outRec, inRec, "CRPIX1");
        drms_copykey(outRec, inRec, "CRPIX2");
        drms_copykey(outRec, inRec, "CRVAL1");
        drms_copykey(outRec, inRec, "CRVAL2");
        drms_copykey(outRec, inRec, "CDELT1");
        drms_copykey(outRec, inRec, "CDELT2");
        drms_copykey(outRec, inRec, "CADENCE");
//        drms_copykey(outRec, inRec, "T_REC_step");
        drms_copykey(outRec, inRec, "CROTA2");
        drms_copykey(outRec, inRec, "CRDER1");
        drms_copykey(outRec, inRec, "CRDER2");
        drms_copykey(outRec, inRec, "CSYSER1");
        drms_copykey(outRec, inRec, "CSYSER2");

        drms_copykey(outRec, inRec, "RSUN_OBS");
        drms_copykey(outRec, inRec, "DSUN_OBS");
        drms_copykey(outRec, inRec, "CRLN_OBS");
        drms_copykey(outRec, inRec, "CRLT_OBS");
        drms_copykey(outRec, inRec, "HGLN_OBS");
        drms_copykey(outRec, inRec, "CAR_ROT");
        drms_copykey(outRec, inRec, "OBS_VR");
        drms_copykey(outRec, inRec, "OBS_VW");
        drms_copykey(outRec, inRec, "OBS_VN");

        drms_copykey(outRec, inRec, "T_OBS");
        drms_copykey(outRec, inRec, "QUALITY");

        drms_copykey(outRec, inRec, "TOTVALS");
        drms_copykey(outRec, inRec, "DATAVALS");
        drms_copykey(outRec, inRec, "MISSVALS");
        drms_copykey(outRec, inRec, "CALVER64");

// write Bradial as the first segment
        outSeg = drms_segment_lookupnum(outRec, 0);
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, bRadial, &status);
        status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          {
            printf("problem writing file, skip! \n");
            continue;
          }
        drms_free_array(outArray);
        drms_free_array(inArray);
      } 

  drms_close_records(inRD, DRMS_FREE_RECORD);
  drms_close_records(outRD, DRMS_INSERT_RECORD);
  return 0;
}

/*
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/remapmags/apps/fdlos2radial.c,v $
$Author: yliu $
*/

/* $Log: fdlos2radial.c,v $
 * Revision 1.3  2013/07/08 18:09:01  yliu
 * corrected a bug
 *
Purpose: Convert los mags to radial mags by assuming that the field is purely radial.
% vectortransform in='su_yang.hmi_vector' out='su_yang.hmi_maghelio'

*/

