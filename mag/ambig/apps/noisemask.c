/*
 */

#include <math.h>
#include <sys/time.h>
//#include "drms_dsdsapi.h"
//#include "soi_ephem.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jsoc_main.h"
#include "imagefromchebyshev.c"

/* cmd-line parameters */
#define kRecSetIn       "in"
#define kSeriesOut      "out"
#define kSegIn          "segin"
#define kSegOut         "segout"
#define kNOTSPECIFIED   "not specified"
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}

int noisemask(xDim, yDim, xcen, ycen, rsun, vrcenter, image) 
    double *image;
    int xDim, yDim, xcen, ycen, rsun;
    float vrcenter;
 {

    CmdParams_t *params = &cmdparams;
    int valid, status = 0;
    DRMS_Segment_t *inSeg, *inSegfinal, *outSeg;
    DRMS_RecordSet_t *inRS, *inRSfinal, *outRS;
    DRMS_Record_t *inRec, *inRecfinal, *outRec;
    DRMS_Array_t *inArray, *outArray;
    char *inQuery="su_yang.cheby_coef", *outQuery;
    char *inRecQuery, *outRecQuery;
    char *inQueryfinal, *vr_start_str, *vr_stop_str;
//    int xDim = 4096, yDim = 4096, order = 15;
    int vr_start, vr_stop, order;
    double *coef, *mask;
    int i, j, t, s, nrecs, ii, jj, k;
    int sunSize = 2 * rsun;  

//    inQuery = (char *)params_get_str(&cmdparams, "su_yang.cheby_coef");
    vr_start = (int)(vrcenter - 50.0);
    vr_stop = (int)(vrcenter + 50.0);
    inQueryfinal = (char *)malloc(100 * sizeof(char));
    vr_start_str = (char *)malloc(100 * sizeof(char));
    vr_stop_str = (char *)malloc(100 * sizeof(char));
    mask = (double *)malloc(sunSize * sunSize * sizeof(double));

//sprintf(inQuery, "su_yang.cheby_coef");
//printf("%s\n", inQuery);

//printf(" I'm here \n");

    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0)
       DIE("No input data found_a");
    inRec = inRS->records[0];

    sprintf(vr_start_str, "%d",  vr_start);
    sprintf(vr_stop_str, "%d",  vr_stop);
    sprintf(inQueryfinal, "%s[%s-%s]", inRec->seriesinfo->seriesname, vr_start_str, vr_stop_str);
    printf("%s\n", inQueryfinal);
    drms_close_records(inRS, DRMS_FREE_RECORD);

    inRSfinal = drms_open_records(drms_env, inQueryfinal, &status);
    if (status || inRSfinal->n == 0) DIE("No input data found_b");
    nrecs = inRSfinal->n;
    inRecfinal = inRSfinal->records[0];
    inSeg = drms_segment_lookupnum(inRecfinal, 0);
    inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
    order = inArray->axis[0];
    coef = (double *)malloc(order * order * sizeof(double));

        for (i = 0; i < order; i++)
        for (j = 0; j < order; j++) {
            coef[i * order + j] = 0;
        }

//printf(" I'm here \n");

    if (nrecs == 1)
      {
        inRecfinal = inRSfinal->records[0];
        inSeg = drms_segment_lookupnum(inRecfinal, 0);
        inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
        if (status)
           {
              printf(" No data file found, status=%d\n", status);
              drms_free_array(inArray);
            }
         double *inData = (double *)inArray->data;
//                order = inArray->axis[0];

                for (jj = 0; jj < order; jj++)
                    {
                    for (ii = 0; ii < order; ii++)
                        {
                          coef[jj * order + ii] = inData[jj * order + ii];
                        }
                    }
      }

     if (nrecs >= 2) 
      {
         for (i = 0; i < 2; i++)
         {
            inRecfinal = inRSfinal->records[i];
            inSeg = drms_segment_lookupnum(inRecfinal, 0);
            inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
            if (status)
               {
                  printf(" No data file found, status=%d\n", status);
                  drms_free_array(inArray);
                  continue;
                }
             double *inData = (double *)inArray->data;
                 for (jj = 0; jj < order; jj++)
                    {
                    for (ii = 0; ii < order; ii++)
                        {
                          coef[jj * order + ii] += 0.5 * inData[jj * order + ii];
                        }
                    }
           }
         }

        imagefromchebyshev(mask, sunSize, sunSize, order, coef);

// paste the mask in the 4096x4096 array

        double yDist, xDist, yDist2, xDist2, xDist2PlusyDist2, dToDiskCenter;
        int jy = 0, yOff, iData, yOffmask, iMask, xdelta, ydelta;
        for (jy = 0; jy < yDim; jy++)
          {
              int ix = 0;
              yDist = (double)jy - ycen;
              yDist2 = yDist * yDist;
              yOff = jy * xDim;
              ydelta = 0;
              if (jy >= ycen-rsun & jy <= ycen+rsun) ydelta = jy - (ycen-rsun);
              yOffmask = ydelta * sunSize;
              for (ix = 0; ix < xDim; ix++)
                {
                   iData = yOff + ix;
                   iMask = 0; 
                   if (ix >= xcen-rsun & ix <= xcen+rsun) iMask = yOffmask + ix - (xcen-rsun); 
                   xDist = (double)ix - xcen;
                   xDist2 = xDist * xDist;
                   xDist2PlusyDist2 = xDist2 + yDist2;
                   dToDiskCenter = sqrt(xDist2PlusyDist2);
                   if (dToDiskCenter >= rsun || iMask >= sunSize * sunSize)
                     {
                         image[iData] = DRMS_MISSING_DOUBLE;
                         continue;
                     }
                   image[iData] = mask[iMask];    
                }
          }
	 free(inQueryfinal); free(vr_start_str); free(vr_stop_str);
//printf(" I'm here \n");
    return 0;
}
/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  v 0.6 08.06.13  R Bogart    first working version
 *
 */
