// module to ingest fits images deconvolved by Aimee with a PSF model into the su_couvidat.lev1 series
//EXAMPLE: _linux_x86_64/proj/lev1.5_hmi/apps/ingest_Aimee fits="/tmp20/norton/scat/deconvolved/image_lev1_01.fits" out="su_couvidat.lev1" inRec="hmi.lev1[2012.06.06_02:27:00.86_UTC]"

#include <jsoc_main.h>
#include <cmdparams.h>
#include <drms.h>
#include <drms_names.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <HMIparam.h>           //contains definitions for some HMI filter parameters
#include "drms_defs.h"
#include "drms_fitsrw.h"
#include "fstats.h"             //header for the statistics function of Keh-Cheng
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>	//for umask(2)
#include <unistd.h>	//for alarm(2) among other things...
#include <printk.h>
#include "/home/jsoc/cvs/Development/JSOC/proj/libs/astro/astro.h"
#include <fresize.h>
#include <gapfill.h>
#include "/home/jsoc/cvs/Development/JSOC/proj/lev0/apps/imgdecode.h"
#include "/home/jsoc/cvs/Development/JSOC/proj/lev0/apps/lev0lev1.h"
#include "/home/jsoc/cvs/Development/JSOC/proj/lev0/apps/limb_fit.h"

char *module_name    = "ingest_dcon_gen";      //name of the module

//arguments of the module
ModuleArgs_t module_args[] = 
{
     {ARG_STRING,"fits", "", "FITS file to be ingested into the series"},
     {ARG_STRING,"out","", "Output series"},
     {ARG_STRING,"inRec", "", "input query"},
     {ARG_END}
};

#include "/home/jsoc/cvs/Development/JSOC/proj/lev0/apps/limb_fit_function.c"

//CORRECTION OF HEIGHT FORMATION
//returns 0 if corrections were successful, 1 otherwise
int heightformation(int FID, double OBSVR, float *CDELT1, float *RSUN, float *CRPIX1, float *CRPIX2, float CROTA2)
{
  int wl=0;
  int status=0;
  float correction=0.0,correction2=0.0;
  
  wl = (FID/10)%20;  //temp is now the filter index
  
  if( (wl >= 0) && (wl < 20) )
    {
      correction  = 0.445*exp(-(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.25)*(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.25)/7.1);
      correction2 = 0.39*(-2.0*(wl-10.- (float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)/6.15)*exp(-(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)*(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)/6.15);
      printf("A CDELT1= %f CROTA2= %f OBS_VR= %f RSUN= %f correction= %f \n",*CDELT1,CROTA2,OBSVR,*RSUN,correction);
      *CDELT1 = *CDELT1*(*RSUN)/((*RSUN)-correction);
      printf("B CDELT1= %f CROTA2= %f OBS_VR= %f \n",*CDELT1,CROTA2,OBSVR);
      *RSUN   = *RSUN-correction;
      *CRPIX1 = *CRPIX1-cos(M_PI-CROTA2*M_PI/180.)*correction2;
      *CRPIX2 = *CRPIX2-sin(M_PI-CROTA2*M_PI/180.)*correction2;
	}
  else status=1;
  
  return status;
  
}

int DoIt(void) 
{
    int status = DRMS_SUCCESS;
    int nRecs;
    DRMS_Type_t type = DRMS_TYPE_FLOAT; 
    float RSUN_LF=0.0, X0_LF=0.0, Y0_LF=0.0;
    double tempRSUN=0.0, tempX0=0.0, tempY0=0.0;
    int FID;
    double OBSVR;
    float CROTA2,CDELT1;

    char *allvers = NULL; /* If 'y', then don't do a 'group by' on the primekey value.
                           * The rationale for this is to allow users to get all versions
                           * of the requested DRMS records */
    char **sets = NULL;
    DRMS_RecordSetType_t *settypes = NULL; /* a maximum doesn't make sense */
    char **snames = NULL;
    char **filts = NULL;
    int nsets = 0;
    DRMS_RecQueryInfo_t rsinfo; /* Filled in by parser as it encounters elements. */
    
    char *inputSeries = NULL;
    
    int nSegs;
    DRMS_Array_t *image = NULL;
    DRMS_Record_t *inputRec = NULL;
    DRMS_Record_t *outputRec = NULL;
    HIterator_t *lastseg = NULL;
    DRMS_Segment_t *inputSeg = NULL;
    DRMS_Segment_t *outputSeg = NULL;
    DRMS_Segment_t *tSeg = NULL;
    
    char filename[DRMS_MAXSEGFILENAME];

    const char *inDir = cmdparams_get_str(&cmdparams, "fits", NULL);
    const char *dsout = cmdparams_get_str(&cmdparams, "out", NULL);
    const char *inRecQuery = cmdparams_get_str(&cmdparams, "inRec", NULL);
  
    /* Parse input series name to extract the series name. */
    if (drms_record_parserecsetspec(inRecQuery, &allvers, &sets, &settypes, &snames, &filts, &nsets, &rsinfo) != DRMS_SUCCESS)
    {
        fprintf(stderr, "Invalid input record-set specification.\n");
        return EXIT_FAILURE;
    }

    if (nsets != 1)
    {
        fprintf(stderr, "Invalid input record-set specification. This module does not support record-set subsets.\n");
        return EXIT_FAILURE;
    }
    
    if (!snames[0] || strlen(snames[0]) == 0)
    {
        fprintf(stderr, "Invalid input record-set specification (invalid series name).\n");
        return EXIT_FAILURE;
    }
    
    inputSeries = strdup(snames[0]);

    /* Free all memory allocated with the drms_record_parserecsetspec() call. */
    drms_record_freerecsetspecarr(&allvers, &sets, &settypes, &snames, &filts, nsets);
    
    /* Check for the existence of the output series. */
    if (!drms_series_exists(drms_env, dsout, &status))
    {
        fprintf(stderr, "Output series %s doesn't exist\n",dsout);  
        return EXIT_FAILURE;
    } 
    
    printf("Output series %s exists.\n", dsout);

    /* Open input records. */
    DRMS_RecordSet_t *data = drms_open_records(drms_env, inRecQuery, &status);   //open the records from the input series
    if (status == DRMS_SUCCESS && data != NULL && data->n > 0)
    {
        nRecs = data->n;
        
        if (nRecs != 1)
        {
            fprintf(stderr, "This module operates on a single input/output record at a time. %d records were specified for input.\n", nRecs);
            exit(EXIT_FAILURE);
        }
        
        inputRec = data->records[0];
        
        printf("Number of input records satisfying the request = %d.\n", nRecs);
    }
    else
    {
        printf("Failure opening at least one input records.\n");
        return EXIT_FAILURE;
    }
    
    /* Open output record. */
    outputRec = drms_create_record(drms_env, dsout, DRMS_PERMANENT, &status);
    if (status != DRMS_SUCCESS || !outputRec)
    {
        fprintf(stderr, "Could not create output record in series %s.\n", dsout);
        return EXIT_FAILURE;
    }

    while ((inputSeg = drms_record_nextseg(inputRec, &lastseg, 0)))
    {
        if (inputSeg->info->islink)
        {
            tSeg = drms_segment_lookup(inputRec, inputSeg->info->name);

            if (!tSeg)
            {
               /* No link set for this record or missing segment struct - skip and continue. */
               continue;
            }
        }
        else
        {
            tSeg = inputSeg;
        }
        
        /* Determine the name of the input FITS file from the segment name. Use the source segment name (if the segment is linked). 
         * We do not need to read the input segment, so do not use drms_segment_filename() to obtain the segment file name. */
        if (tSeg->info->protocol == DRMS_TAS)
        {
            fprintf(stderr, "This module supports FITS segments only (the input series has a TAS segment).\n");
            return EXIT_FAILURE;
        }
        else
        {
            char fitsPath[PATH_MAX];

            snprintf(fitsPath, sizeof(fitsPath), "%s", inDir);
            
            if (rindex(inDir, '/'))
            {
                fitsPath[sizeof(fitsPath) - 1] = '\0';
            }
            
            snprintf(filename, sizeof(filename), "%s/%s%s", inDir, inputSeg->info->name, drms_prot2ext(inputSeg->info->protocol));
        }

        /* Read the input FITS file for this segment. */
        image = drms_fitsrw_read(drms_env, filename, 0, NULL, &status);
        
        if (status != DRMS_SUCCESS || !image)
        {
            fprintf(stderr, "Unable to read input FITS file %s.\n", filename);
            return EXIT_FAILURE;
        }
        
        /* Open output segment. */
        outputSeg = drms_segment_lookup(outputRec, inputSeg->info->name);
        
        if (!outputSeg)
        {
            fprintf(stderr, "Unable to find output segment %s.\n", inputSeg->info->name);
            return EXIT_FAILURE;
        }
        
        image->bzero = outputSeg->bzero;
        image->bscale = outputSeg->bscale;
        image->israw = 0;

        printf("Writing a segment to DRMS, %s.\n", outputSeg->info->name);
        if (drms_segment_write(outputSeg, image, 0) != DRMS_SUCCESS)
        {
            fprintf(stderr, "Unable to write output segment %s.\n", outputSeg->info->name);
            return EXIT_FAILURE;
        }
        
        drms_free_array(image);
    } /* loop over segments */
    
    if (lastseg)
    {
        hiter_destroy(&lastseg);
    }
    
    /* WARNING! DANGER, DANGER! This module does not currently re-calculate the output record's 
    * CRPIX1, CRPIX2, R_SUN, or CDELT1 keyword values. The deconvolution process alters
    * these values, but the values that show in the output series are the values that 
    * existed prior to deconvolution. 
    *
    * TODO - figure out how to re-calculate these keywords.
    */
    if (drms_copykeys(outputRec, inputRec, 0, kDRMS_KeyClass_Explicit) != DRMS_SUCCESS)
    {
        fprintf(stderr, "Could not copy keywords to output series record.\n");
        return EXIT_FAILURE;
    }

    drms_close_record(outputRec, DRMS_INSERT_RECORD); /* insert the record into DRMS */
 
    return 0;
}
