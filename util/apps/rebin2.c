/**
   @defgroup rebin2 rebin2 reduce/increase resolution to that of MDI/HMI
   @ingroup su_util

   @brief Reduce (increase) the resolution of the input data to that of MDI (HMI) by factor of 2.

   @par Synopsis:
   @code
   rebin2  in=input data out=output data  factor=2*n mode=simple
   where n is a multiple/fraction of 2.
   @endcode

   This is a general purpose module that takes a series of input 
   data and modifies its spatial resolution by a factor (multiples 
   or fractions of 2) as required and gives out a set of output 
   data. The method for avaraging (interpolation) can be specified 
   through the input "mode". The current version handles a simple 
   boxcar average. If 'scale' < 0 then the input is reduced in size to 1/|scale|.

   Make sure you created the appropriate output series before running 
   the program. For example, su_bala.rebin2up.jsd

   @par Flags:
   @c none
   Currently it doesn't have any flags.

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in  The input data series.
   @param out The output series.

   @par Exit_Status:
   Brief description of abnormal, non-zero, exit values.

   @par Example:
   Takes a series of 1024 X 1024 MDI full disk magnetogram and
   produces images with the resolution of HMI, 4096 X 4096.

   @code
   rebin2 in='mdi.fd_M_96m_lev18[2003.10.20/1d]' out='su_bala.rebin2up' scale=4 mode='simple'
   @endcode

   @par Example:
   Reduces the resolution of HMI images 4096 X 4096 to that of MDI,
   1024 X 1024. Here the input is the HMI images and the output
   is the lower resolution HMI images.
   @code
   rebin2 in='su_bala.rebin2up[2003.10.20/1d]' out='su_bala.rebin2down' scale=0.25 mode='simple'
   @endcode

   @bug
   None known so far.

   @par Code:
   The doxygen code that makes this page is here:
   @verbinclude  ./doxygen_moduletemplate.txt
*/

#include "jsoc.h"
#include "jsoc_main.h"
#include "fstats.h"

char *module_name = "rebin2";

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
     {ARG_STRING, "out", "NOT SPECIFIED",  "Output data series."},
     {ARG_FLOAT, "scale", "NOTSPECIFIED", "Scale factor. negative for shrink factor"},
     {ARG_STRING, "mode", "simple", "conversion type."},
     {ARG_END}
};

int DoIt(void)
  {
  int status = DRMS_SUCCESS;
  DRMS_RecordSet_t *inRS, *outRS;
  int irec, nrecs;
  const char *inQuery = params_get_str(&cmdparams, "in");
  const char *outSeries = params_get_str(&cmdparams, "out");
  const char *mode = params_get_str(&cmdparams, "mode");
  int scale = params_get_int(&cmdparams, "scale");
  int iscale, factor;
  float fscale;
  if (scale < 0.0)
    { // shrinking
    iscale = -scale;
    fscale = 1.0/iscale;
    }
  else
    { // enlarging
    iscale = (int) scale;
    fscale = scale;
    }

  inRS = drms_open_records(drms_env, inQuery, &status);
  if (status || inRS->n == 0)
    DIE("No input data found");
  nrecs = inRS->n;
 if (nrecs == 0)
    DIE("No records found");
  outRS = drms_create_records(drms_env, nrecs, (char *)outSeries, DRMS_PERMANENT, &status);
  if (status)
    DIE("Output recordset not created");

  for (irec=0; irec<nrecs; irec++)
    {
    DRMS_Record_t *outRec, *inRec;
    DRMS_Segment_t *outSeg, *inSeg;
    DRMS_Array_t *inArray, *outArray;
    float *inData, *outData;
    int inx, iny, outx, outy, i, j;
    float val;
    int quality=0;
 
    inRec = inRS->records[irec];
    quality = drms_getkey_int(inRec, "QUALITY", &status);
    if (status || (!status && quality >= 0))
      {
      inSeg = drms_segment_lookupnum(inRec, 0);
      inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
      if (status)
        {
        printf(" No data file found but QUALITY not bad, status=%d\n", status);
        drms_free_array(inArray);
        continue;
        }
  
      int naxis = inArray->naxis;
      int in_nx = inArray->axis[0];
      int in_ny = inArray->axis[1];
      int in_inc = (scale < 0 ? 1 : iscale);
      int out_inc = (scale < 0 ? iscale : 1);
      int out_nx = in_nx * fscale;
      int out_ny = in_ny * fscale;
  
      int outDims[2] = {out_nx, out_ny};
      inData = (float *)inArray->data;
      outArray = drms_array_create(DRMS_TYPE_FLOAT, naxis, outDims, NULL, &status);
      outData = (float *)outArray->data;
  
      if (strcmp(mode, "simple")==0) // setup for better HMI proxy
        {
        if (scale < 0) // need to average input chunk
          {
          for (outy = 0; outy < out_ny; outy += 1)
            for (outx = 0; outx < out_nx; outx += 1)
              {
              double total = 0.0;
              int nn = 0;
              for (j = 0; j < in_inc; j += 1)
                for (i = 0; i < in_inc; i += 1)
                  {
                  val = inData[in_nx*(outy*iscale + j) + outx*iscale + i];
                  if (!drms_ismissing_float(val))
                      {
                       total = total + val; 
                       nn++;
                      }
                  }
              outData[out_nx*outy + outx] = (nn > 0 ? total/nn : DRMS_MISSING_FLOAT); 
              }
          }
        else  // need to replicate input point
          {
          for (iny = 0; iny < in_ny; iny += 1)
            for (inx = 0; inx < in_nx; inx += 1)
              {
              val = inData[in_nx*iny + inx];
              for (j = 0; j < out_inc; j += 1)
                for (i = 0; i < out_inc; i += 1)
                  outData[out_nx*(iny*iscale + j) + inx*iscale + i] = val;
              }
          }
        }
    else
      DIE("invalid conversion mode");
  
    drms_free_array(inArray);
  
    // write data file
    outRec = outRS->records[irec];
    outSeg = drms_segment_lookupnum(outRec, 0);

    // copy all keywords
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);
    set_statistics(outSeg, outArray, 1);

    // Now fixup coordinate keywords
    // Only CRPIX1,2 and CDELT1,2 should need repair.
    drms_setkey_double(outRec, "CDELT1", drms_getkey_double(inRec, "CDELT1", NULL)/fscale);
    drms_setkey_double(outRec, "CDELT2", drms_getkey_double(inRec, "CDELT2", NULL)/fscale);
    drms_setkey_double(outRec, "CRPIX1", (drms_getkey_double(inRec, "CRPIX1", NULL)-0.5) * fscale + 0.5);
    drms_setkey_double(outRec, "CRPIX2", (drms_getkey_double(inRec, "CRPIX2", NULL)-0.5) * fscale + 0.5);
    drms_setkey_int(outRec, "SCALE", scale);

  
    // get info for array from segment
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    outArray->parent_segment = outSeg;
  
    status = drms_segment_write(outSeg, outArray, 0);
    if (status)
      DIE("problem writing file");
    drms_free_array(outArray);
    }
   } //end of "irec" loop

  drms_close_records(inRS, DRMS_FREE_RECORD);
  drms_close_records(outRS, DRMS_INSERT_RECORD);
  return (DRMS_SUCCESS);
  } // end of DoIt
