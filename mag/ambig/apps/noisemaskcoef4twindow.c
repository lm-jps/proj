/*
* This module derives the noise mask of full-disk vector magnetogram from the FD10 field strength.
* Change of phase map is taken into account.
*/

#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jsoc_main.h"
#include "fitimage_float.c"

#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include <omp.h>
#include "fresize.c"

#define PI              (M_PI)
#define RADSINDEG       (PI/180)

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}

/* cmd-line parameters */
#define kRecSetIn       "in"
#define kSeriesOut      "out"
#define kSegIn          "segin"
#define kSegOut         "segout"
#define kNOTSPECIFIED   "not specified"

#define RADSINDEG       (PI/180)
#define RAD2ARCSEC      (648000. / M_PI)

/* ################## Wrapper for Jesper's code ################## */

void frebin(float *image_in, float *image_out, int nx_in, int ny_in, int nx_out, int ny_out, int nbin)
{
  struct fresize_struct fresizes;
//  int nxout, nyout;
  int nlead_in = nx_in, nlead_out = nx_out;

//  nxout = nx / nbin; nyout = ny / nbin;
  init_fresize_gaussian2(&fresizes, nbin, (3 * nbin)/2, (3 * nbin)/2, nbin);

//  fresize(&fresizes, image_in, image_out, nx, ny, nlead, nxout, nyout, nxout, (nbin / 2) * 2, (nbin / 2) * 2, DRMS_MISSING_FLOAT);

  fresize(&fresizes, image_in, image_out, nx_in, ny_in, nlead_in, nx_out, ny_out, nlead_out, 1, 1, DRMS_MISSING_FLOAT);
// from Jesper: I think that, in general, you want offsets of 1 if nbin=3 and nbin
//              divides into the original image size. That way the output image is effectively
//              centered on the input image.
  

  free_fresize(&fresizes);

}

char *module_name = "noisemaskcoef4twindow";
float median(int n, float *x);

ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kRecSetIn, "", "Input data records."},
   {ARG_STRING,  kSeriesOut, "", "Output data series."},
   {ARG_STRING,  kSegIn, kNOTSPECIFIED, ""},
   {ARG_STRING,  kSegOut, kNOTSPECIFIED, ""},
   {ARG_INT, "ORDER", "15", "", ""},
   {ARG_INT, "MASKID", "1", "", ""},
   {ARG_FLOAT, "VR_START", "-3600.0", "", ""},
   {ARG_FLOAT, "VR_STOP", "3600.0", "", ""},
   {ARG_END}
};

int DoIt(void)
{

    int status = DRMS_SUCCESS, newchunk, ds, nds, nnds;
    float vr, vrstart, vrstop, vrstep = 50.0, vr_center;
    char *inRecQuery, *outRecQuery;
    char *phasemapid_str, phaseid[64];
    int phasemapid;
    char tstr[64];
    int y,m,d,maskid;
    int imrec[36000]; 
    int ngood, ndMax = 50;
    float vr_left, vr_right, dSun, rSun_ref, cdelt, asd, rSun, xCenter, yCenter, bCutoff = 300.0;
    int xDim = 4096, yDim = 4096;
    long long imgData;
    double xDist = 0.0, yDist = 0.0, xDist2 = 0.0, yDist2 = 0.0, xDist2PlusyDist2 = 0.0, dToDiskCenter = 0.0;
    int mm, nn, order, yOff = 0, iData = 0;
    DRMS_RecordSet_t *inRD, *outRD;
    DRMS_Record_t *inRec;
    DRMS_Record_t *outRec;
    DRMS_Segment_t *inSeg = NULL, *outSeg;
    DRMS_Array_t *inArray, *outArray;
    DRMS_RecChunking_t cstat;
    int outDims[2] = {xDim, yDim};
    TIME t_rec;

    phasemapid_str = (char *)malloc(100 * sizeof(char));
    inRecQuery = (char *)params_get_str(&cmdparams, kRecSetIn);
    outRecQuery = (char *)params_get_str(&cmdparams, kSeriesOut);
    maskid = cmdparams_get_int(&cmdparams, "MASKID", &status);
    vrstart = cmdparams_get_float(&cmdparams, "VR_START", &status);
    vrstop = cmdparams_get_float(&cmdparams, "VR_STOP", &status);
    order = cmdparams_get_int(&cmdparams, "ORDER", &status);

    nds = (int)((vrstop - vrstart)/vrstep) + 1;
//    outRD = drms_create_records(drms_env, nds, outRecQuery, DRMS_PERMANENT, &status);
//    if (status) DIE("Output recordset not created");

  for (ds = 0; ds < nds; ds++)
  {
    float *bPixel, *datacube, *bField, *bMask;
    double *coef;
    int idx = 0;
    float crpix1, crpix2, crval1, crval2, cdelt1, cdelt2, crota2;
//    outRec = outRD->records[ds];
    vr_center = vrstart + ds * vrstep; 
    vr_left = vr_center - vrstep;
    vr_right = vr_center + vrstep;

    datacube = (float *) malloc(xDim * yDim * ndMax * sizeof(float));
    bMask = (float *) malloc(xDim * yDim * sizeof(float));

    nnds = drms_count_records(drms_env, inRecQuery, &status);
    inRD = drms_open_recordset(drms_env, inRecQuery, &status);
printf("n=%d\n", nnds);
    while ((inRec = drms_recordset_fetchnext(drms_env, inRD, &status, &cstat, &newchunk)) != NULL)
      {
        vr = drms_getkey_double(inRec, "OBS_VR", &status);
        if (vr > vr_left && vr < vr_right)
          { 
            ++idx;
            if (idx == 1)
              {
                dSun = drms_getkey_double(inRec, "DSUN_OBS", &status);
                rSun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
                if (status) rSun_ref = 6.96e8;
                cdelt = drms_getkey_float(inRec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
                asd = asin(rSun_ref/dSun);
                rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
                xCenter = drms_getkey_float(inRec, "CRPIX1", &status) - 1.0;
                yCenter = drms_getkey_float(inRec, "CRPIX2", &status) - 1.0;
                phasemapid = drms_getkey_int(inRec, "INVPHMAP", &status);

                sprintf(phasemapid_str, "%d",  phasemapid);
                printf("phase ID = %s\n", phasemapid_str);

                t_rec = drms_getkey_time(inRec, "T_REC", &status);

                crpix1 = drms_getkey_float(inRec, "CRPIX1", &status);
                crpix2 = drms_getkey_float(inRec, "CRPIX2", &status);
                crval1 = drms_getkey_float(inRec, "CRVAL1", &status);
                crval2 = drms_getkey_float(inRec, "CRVAL2", &status);
                cdelt1 = drms_getkey_float(inRec, "CDELT1", &status);
                cdelt2 = drms_getkey_float(inRec, "CDELT2", &status);
                crota2 = drms_getkey_float(inRec, "CROTA2", &status);
              }

            t_rec = drms_getkey_time(inRec, "T_REC", &status);
            vr = drms_getkey_double(inRec, "OBS_VR", &status);
            sprint_at(tstr, t_rec);
            inSeg = drms_segment_lookup(inRec, "field");
            inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
            if (status)
             {
              printf("Bad data, skip! \n");
              --idx; continue;
             }
            bField = (float *)inArray->data;
            imgData = (idx - 1) * xDim * yDim;
            int jy = 0;
            for (jy = 0; jy < yDim; jy++)
              {
                int ix = 0;
                yOff = jy * xDim;
                for (ix = 0; ix < xDim; ix++)
                  {
                   iData = yOff + ix;
                   datacube[imgData + iData] = bField[iData];
                  }
              }
              printf("idx=%d, Vr=%f, T_REC = %s\n", idx, vr, tstr);
              drms_free_array(inArray);

          }
          if (idx > ndMax - 1) break; //maximum data number is set to be 50.
        }

      ngood = idx;
printf("ngood=%d\n", ngood);
      if (ngood < 10) {free(datacube); continue;}
      else
      {
        bPixel = (float *) malloc(ngood * sizeof(float));

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
                   if (dToDiskCenter >= rSun) bMask[iData] = DRMS_MISSING_FLOAT;
                   else
                     {
                       int countPixel = 0, iimg = 0;
                       imgData = 0;
                       for (iimg = 0; iimg < ngood; iimg++)
                        {
                          imgData = iimg * xDim * yDim + iData;
                          if (datacube[imgData] < bCutoff)
                             {
                               bPixel[countPixel] = datacube[imgData];
                               ++countPixel;
                             }
                         }
                        if (countPixel == 0) bMask[iData] = DRMS_MISSING_FLOAT;
                        else if (countPixel == 1) bMask[iData] = bPixel[0];
                        else
                          {
                           bMask[iData] = median(countPixel, bPixel);
                          }
                      }
                  }
             }

// smooth by convolving a gaussian function -- call fresize code (Schou)
       int nbin = 8, xout, yout, xleft, yleft, jData;
       float *outData, *outResize, xcResize, ycResize, rResize;
       double *mask;

       xout = xDim/nbin; yout = yDim/nbin;
       outData = (float *) malloc(xout * yout * sizeof(float)); 
       frebin(bMask, outData, xDim, yDim, xout, yout, nbin);

       xcResize = xCenter/nbin;
       ycResize = yCenter/nbin;
       rResize = rSun/nbin;
       xleft = (int)(xcResize - rResize) - 1;
       yleft = (int)(ycResize - rResize) - 1;
       outDims[0] = (int)(2*rResize) + 3; outDims[1] = (int)(2*rResize) + 3;
       outResize = (float *) malloc(outDims[0] * outDims[1] * sizeof(float));

       for (jy = 0; jy < outDims[1]; jy++)
         {
           int ix = 0;
           jData = (jy + yleft) * xout;
           for (ix = 0; ix < outDims[0]; ix++)
             {
               imgData = jData + ix + xleft;
               outResize[jy * outDims[0] + ix] = outData[imgData];
             }
          }

    mm = outDims[0]; nn = outDims[1];
    coef = (double *)malloc(order * order * sizeof(double));
    fitimage_float(outResize, mm, nn, order, 1, coef);

// write the output

        outRD = drms_create_records(drms_env, 1, outRecQuery, DRMS_PERMANENT, &status);
        if (status) DIE("Output recordset not created");
        outRec = outRD->records[0];
        drms_setkey_time(outRec, "T_REC_O", t_rec);

/*
                drms_copykey(outRec, inRec, "INVPHMAP");
                drms_copykey(outRec, inRec, "INSTRUME");
                drms_copykey(outRec, inRec, "CAMERA");
                drms_copykey(outRec, inRec, "HISTORY");
                drms_copykey(outRec, inRec, "COMMENT");
*/

                drms_setkey_float(outRec, "CRPIX1_O", crpix1);
                drms_setkey_float(outRec, "CRPIX2_O", crpix2);
                drms_setkey_float(outRec, "CRVAL1_O", crval1);
                drms_setkey_float(outRec, "CRVAL2_O", crval2);
                drms_setkey_float(outRec, "CDELT1_O", cdelt1);
                drms_setkey_float(outRec, "CDELT2_O", cdelt2);
                drms_setkey_float(outRec, "CROTA2_O", crota2);
                drms_setkey_double(outRec, "DSUN_ORI", dSun);
                drms_setkey_double(outRec, "RSUN_ORI", rSun);
                drms_setkey_int(outRec, "VRCENT", (int)(vr_center));

        outDims[0] = order; outDims[1] = order;      
        outSeg = drms_segment_lookupnum(outRec, 0);
        outArray = drms_array_create(DRMS_TYPE_DOUBLE, 2, outDims, coef, &status);
        status = drms_segment_write(outSeg, outArray, 0);
        if (status) printf("problem writing file, skip! \n");
        drms_setkey_int(outRec, "NGOOD", ngood);
        drms_setkey_int(outRec, "POLYDEGR", order);
        drms_setkey_int(outRec, "MASKID", maskid);
        drms_setkey_float(outRec, "BCUTOFF", bCutoff);
        drms_setkey_float(outRec, "DELTAVR", vrstep);
        sprint_ut(tstr, CURRENT_SYSTEM_TIME);
        sscanf(tstr, "%d.%d.%d", &y, &m, &d);
        sprintf(tstr, "%04d-%02d-%02d", y, m, d);
        drms_setkey_string(outRec, "DATE", tstr);
        drms_free_array(outArray);
        free(outResize); free(outData); free(bPixel); free(datacube);
        drms_close_records(outRD, DRMS_INSERT_RECORD);
    }
  drms_close_records(inRD, DRMS_FREE_RECORD);
  }
//  drms_close_records(outRD, DRMS_INSERT_RECORD);  
  return 0;
}

float median(int n, float x[]) {
    float temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }
 
    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

/*
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/ambig/apps/noisemaskcoef4twindow.c,v $
$Author: yliu $
*/

