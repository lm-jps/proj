/*
 */

#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jsoc_main.h"
#include "imagefromchebyshev.c"

// Dec 16 Xudong
#define kMaskInquery	"hmi.lookup_ChebyCoef_TWindow"

int noisemaskimag4twindow(xDim, yDim, xcen, ycen, rsun, vrcenter, maskid, image, maskQuery) 
    double *image;
    int xDim, yDim, maskid;
    float xcen, ycen, rsun;
    float vrcenter;
		char *maskQuery;        // Sep 25
 {
    CmdParams_t *params = &cmdparams;
    int valid, status = 0;
    DRMS_Segment_t *inSeg, *inSegfinal, *outSeg;
    DRMS_RecordSet_t *inRS, *inRSfinal, *outRS;
    DRMS_Record_t *inRec, *inRecfinal, *outRec;
    DRMS_Array_t *inArray, *outArray;
    char *inQuery=kMaskInquery, *outQuery;
    char *inRecQuery, *outRecQuery;
    char *inQueryfinal, *vr_start_str, *vr_stop_str, *maskid_str;
    char ttemp[64];
    int vr_start, vr_stop, order, vr_coef;
    float weight;
    double *coef, *mask, xc_shift, yc_shift;
    int i, j, t, s, nrecs, ii, jj, k;
    int sunSize = (int)(2 * (rsun+1));  

    vr_start = (int)(vrcenter - 50.0);
    vr_stop = (int)(vrcenter + 50.0);
    inQueryfinal = (char *)malloc(100 * sizeof(char));
    vr_start_str = (char *)malloc(100 * sizeof(char));
    vr_stop_str = (char *)malloc(100 * sizeof(char));
    maskid_str = (char *)malloc(100 * sizeof(char));
    mask = (double *)malloc(sunSize * sunSize * sizeof(double));

    sprintf(vr_start_str, "%d",  vr_start);
    sprintf(vr_stop_str, "%d",  vr_stop);
    sprintf(maskid_str, "%d",  maskid);
    sprintf(inQueryfinal, "%s[%s-%s][%s]", inQuery, vr_start_str, vr_stop_str, maskid_str);			// Dec 16 Xudong
    sprintf(maskQuery, "%s", inQueryfinal);		// Dec 16 Xudong
    printf("%s\n", inQueryfinal);

    inRSfinal = drms_open_records(drms_env, inQueryfinal, &status);
    if (status || inRSfinal->n == 0) {printf("No input data found\n"); return (status ? status : -1);};		// Dec 16
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

    if (nrecs == 1)
      {
        inRecfinal = inRSfinal->records[0];
        inSeg = drms_segment_lookupnum(inRecfinal, 0);
        inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
        if (status)
           {
              printf(" No data file found, status=%d\n", status);
              drms_free_array(inArray);  return (status ? status : -1);
            }
         double *inData = (double *)inArray->data;

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
             vr_coef = drms_getkey_int(inRecfinal, "VRCENT", &status);
             weight = 1.0 - fabs((float)(vr_coef) - vrcenter)/50.0;
             double *inData = (double *)inArray->data;
                 for (jj = 0; jj < order; jj++)
                    {
                    for (ii = 0; ii < order; ii++)
                        {
                          coef[jj * order + ii] += weight * inData[jj * order + ii];
                        }
                    }
           }
         }

        xc_shift = (double)xcen - (int)xcen;
        yc_shift = (double)ycen - (int)ycen;

        imagefromchebyshev(mask, sunSize, sunSize, order, coef, xc_shift, yc_shift);
printf("xc-shift = %f, yc-shift = %f\n", xc_shift, yc_shift);

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
    free(coef); free(mask);
    drms_free_array(inArray);
    drms_close_records(inRSfinal, DRMS_FREE_RECORD);
    return 0;
}
/*
 */
