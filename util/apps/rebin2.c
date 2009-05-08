/**
   @defgroup rebin2 rebin2 reduce/increase resolution to that of MDI/HMI
   @ingroup a_programs

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
   boxcar average.

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
   rebin2 in='mdi.fd_M_96m_lev18[2003.10.20/1d]' out='su_bala.rebin2up' factor=4 mode='simple'
   @endcode

   @par Example:
   Reduces the resolution of HMI images 4096 X 4096 to that of MDI,
   1024 X 1024. Here the input is the HMI images and the output
   is the lower resolution HMI images.
   @code
   rebin2 in='su_bala.rebin2up[2003.10.20/1d]' out='su_bala.rebin2down' factor=0.25 mode='simple'
   @endcode

   @bug
   None known so far.

   @par Code:
   The doxygen code that makes this page is here:
   @verbinclude  ./doxygen_moduletemplate.txt
*/

#include "jsoc_main.h"

char *module_name = "rebin2";

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
     {ARG_STRING, "out", "NOT SPECIFIED",  "Output data series."},
     {ARG_FLOAT, "nfactor", "NOTSPECIFIED", "Reduction factor."},
     {ARG_STRING, "mode", "simple", "conversion type."},
     {ARG_END}
};

int DoIt(void)
  {
  int resoln, factor;
  int status = DRMS_SUCCESS;
  DRMS_RecordSet_t *inRS, *outRS;
  int irec, nrecs;
  char *inQuery = params_get_str(&cmdparams, "in");
  char *outQuery = params_get_str(&cmdparams, "out");
  char *mode = params_get_str(&cmdparams, "mode");
  float nfactor = params_get_float(&cmdparams, "nfactor");
  if (nfactor < 1.0)
    {
     resoln = 1;
     factor = (int) (1/nfactor);
    }
  else
    {
     resoln = 0;
     factor = (int) nfactor;
    }

  inRS = drms_open_records(drms_env, inQuery, &status);
  if (status || inRS->n == 0)
    DIE("No input data found");
  nrecs = inRS->n;
  outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
  if (status)
    DIE("Output recordset not created");

 switch (resoln){
   case 1:
/* reducing the resolution (of HMI data) to that of MDI
   resoln = 1
*/ 
    for (irec=0; irec<nrecs; irec++)
     {
      DRMS_Record_t *outRec, *inRec;
      DRMS_Segment_t *outSeg, *inSeg;
      DRMS_Array_t *inArray, *outArray;
      float *inData, *outData;
      int inx, iny, outx, outy, i, j, nn;
      float val;
      float total = 0.;
 
      inRec = inRS->records[irec];
      inSeg = drms_segment_lookupnum(inRec, 0);
      inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
      if (status)
       {
        printf(" No data file found, status=%d\n", status);
        drms_free_array(inArray);
        continue;
       }
      int naxis = inArray->naxis;
      int inDims[2] = {inArray->axis[0], inArray->axis[1]};
      int outDims[2] = {inDims[0]/factor, inDims[1]/factor};
      inData = (float *)inArray->data;
      outArray = drms_array_create(DRMS_TYPE_FLOAT, naxis, outDims, NULL, &status);
      outData = (float *)outArray->data;

      if (strcmp(mode, "simple")==0) // setup for better HMI proxy
       {
        for (iny = 0; iny < inDims[1]; iny++)
         {
          outy = iny/factor;
          for (inx = 0; inx < inDims[0]; inx++)
           {
            outx = inx/factor;
            for (j = 0; j < factor; j++)
              for (i = 0; i < factor; i++)
               {
                val = inData[inDims[0]*(iny+j) + inx + i];
                if (!drms_ismissing_float(val))
                 {
                  total = total + val; 
                  nn++;
                 }
                else 
                total = drms_ismissing_float(val); 
              }       //i loop
              if (!drms_ismissing_float(total)) 
                outData[outDims[0]*outy + outx] = total/nn; 
              else 
                outData[outDims[0]*outy + outx] = drms_ismissing_float(val); 
              total = 0.;
              nn = 0;
           }         // end of "inx"
         }           // end of "iny"
      }             // end of "mode" loop
    else
      DIE("invalid conversion mode");

    drms_free_array(inArray);

    // write data file
    outRec = outRS->records[irec];
    outSeg = drms_segment_lookupnum(outRec, 0);

    // get info for array from segment
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    outArray->parent_segment = outSeg;

    status = drms_segment_write(outSeg, outArray, 0);
    if (status)
      DIE("problem writing file");
    drms_free_array(outArray);

    // copy all keywords
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);
 
    // Now fixup coordinate keywords
    // Only CRPIX1,2 and CDELT1,2 should need repair.
    drms_setkey_double(outRec, "CDELT1", drms_getkey_double(inRec, "CDELT1", NULL)*4);
    drms_setkey_double(outRec, "CDELT2", drms_getkey_double(inRec, "CDELT2", NULL)*4);
    drms_setkey_double(outRec, "CRPIX1", (drms_getkey_double(inRec, "CRPIX1", NULL)-1) /4.0 + 1.0);
    drms_setkey_double(outRec, "CRPIX2", (drms_getkey_double(inRec, "CRPIX2", NULL)-1) /4.0 + 1.0);
    } //end of "irec" loop
/*  
  -----------------  end of resoln 1 ------------------- */
     printf("completing %d \n", resoln);
      break;
    case 0:
/* --------- increasing the resolation to that of HMI ---------
   resol = 0
*/

  for (irec=0; irec<nrecs; irec++)
    {
    DRMS_Record_t *outRec, *inRec;
    DRMS_Segment_t *outSeg, *inSeg;
    DRMS_Array_t *inArray, *outArray;
    float *inData, *outData;
    int inx, iny, outx, outy, i, j;

    inRec = inRS->records[irec];
    inSeg = drms_segment_lookupnum(inRec, 0);
    inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
    if (status)
      {
      printf(" No data file found, status=%d\n", status);
      drms_free_array(inArray);
      continue;
      }
      int naxis = inArray->naxis;
      int inDims[2] = {inArray->axis[0], inArray->axis[1]};
      int outDims[2] = {inDims[0]*factor, inDims[1]*factor};
    inData = (float *)inArray->data;
    outArray = drms_array_create(DRMS_TYPE_FLOAT, naxis, outDims, NULL, &status);
    outData = (float *)outArray->data;


    if (strcmp(mode, "simple")==0) // setup for better HMI proxy
      {
       for (iny = 0; iny < inDims[1]; iny++)
        {
         outy = factor*iny;
         for (inx = 0; inx < inDims[0]; inx++)
          {
           float val = inData[inDims[0]*iny + inx];
if (iny==512 && inx == 512) fprintf(stderr,"512,512 = %f\n",val);
           outx = factor*inx;
           for (j = 0; j < factor; j++)
             for (i = 0; i < factor; i++)
               outData[outDims[0]*(outy+j) + outx + i] = val;
          }
        }
      }
    else
      DIE("invalid conversion mode");

    drms_free_array(inArray);

    // write data file
    outRec = outRS->records[irec];
    outSeg = drms_segment_lookupnum(outRec, 0);

    // get info for array from segment
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    outArray->parent_segment = outSeg;

    status = drms_segment_write(outSeg, outArray, 0);
    if (status)
      DIE("problem writing file");
    drms_free_array(outArray);

    // copy all keywords
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);

    // Now fixup coordinate keywords
    // Only CRPIX1,2 and CDELT1,2 should need repair.
    drms_setkey_double(outRec, "CDELT1", drms_getkey_double(inRec, "CDELT1", NULL)/4.0);
    drms_setkey_double(outRec, "CDELT2", drms_getkey_double(inRec, "CDELT2", NULL)/4.0);
    drms_setkey_double(outRec, "CRPIX1", (drms_getkey_double(inRec, "CRPIX1", NULL)-1) *4.0 + 1.0);
    drms_setkey_double(outRec, "CRPIX2", (drms_getkey_double(inRec, "CRPIX2", NULL)-1) *4.0 + 1.0);
    }
/*
 ------------------------- end resol = 0 -----------------*/
     printf("ciompleting %d \n", resoln);
     break;
  } //end switch case

  drms_close_records(inRS, DRMS_FREE_RECORD);
  drms_close_records(outRS, DRMS_INSERT_RECORD);
  return (DRMS_SUCCESS);
} // end of DoIt
