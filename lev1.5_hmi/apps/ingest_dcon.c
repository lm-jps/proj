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

char *module_name    = "ingest_Aimee";      //name of the module

//arguments of the module
ModuleArgs_t module_args[] = 
{
     {ARG_STRING,"fits", "","FITS file to be ingested into the series"},
     {ARG_STRING,"out","","Output series"},
     {ARG_STRING,"inRec", "","input query"},
     {ARG_STRING, "dpath","/home/jsoc/cvs/Development/JSOC","path of development code"},
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

int DoIt(void) {

  int status        = 0;
  int nRecs;
  DRMS_Type_t type = DRMS_TYPE_FLOAT; 
  float RSUN_LF=0.0, X0_LF=0.0, Y0_LF=0.0;
  double tempRSUN=0.0, tempX0=0.0, tempY0=0.0;
  int FID;
  double OBSVR;
  float CROTA2,CDELT1;

  const char *infile=cmdparams_get_str(&cmdparams, "fits", NULL);
  const char *dsout =cmdparams_get_str(&cmdparams, "out", NULL);
  const char *inRecQuery=cmdparams_get_str(&cmdparams, "inRec", NULL);
  const char *dpath = cmdparams_get_str(&cmdparams,"dpath",NULL);

  //CHECK WHETHER OUTPUT SERIES EXISTS
  drms_series_exists(drms_env, dsout, &status);                      //check whether or not output series exists
  if(status != DRMS_SUCCESS)                                        //drms_env is defined in jsoc_main.h  if (status == DRMS_ERROR_UNKNOWNSERIES)
    {
      printf("Output series %s doesn't exist\n",dsout);  //if the output series does not exit
      exit(EXIT_FAILURE);                                            //we exit the program
    } 
  if(status == DRMS_SUCCESS)
    {
      printf("Output series %s exists.\n",dsout);
    }

  //READ THE INPUT FITS FILE (lev1 image)
  DRMS_Array_t *image;
  image=drms_fitsrw_read(drms_env,infile,0,NULL,&status);
    

  //OPEN hmi.lev1 RECORDS
  DRMS_RecordSet_t *data = drms_open_records(drms_env,inRecQuery,&status);   //open the records from the input series
  if (status == DRMS_SUCCESS && data != NULL && data->n > 0)
    {
      nRecs = data->n;                                           //number of records in the input series 
      printf("Number of level 1 records satisfying the request= %d \n",nRecs);    //check if the number of records is appropriate
    }
  else
    {
      printf("Input level 1 series %s doesn't exist\n",inRecQuery);//if the input series does not exit
      exit(EXIT_FAILURE);                                            //we exit the program
    }

 //CREATE A RECORD IN OUTPUT SERIES
  DRMS_RecordSet_t *dataout = NULL;
  DRMS_Record_t    *recout  = NULL;
  DRMS_Segment_t   *segout  = NULL;
  dataout = drms_create_records(drms_env,1,dsout,DRMS_PERMANENT,&status);
  if (status != DRMS_SUCCESS)
    {
      printf("Could not create a record for the lookup tables\n");
      exit(EXIT_FAILURE);
    }
  recout = dataout->records[0];

  //READ THE INPUT FITS FILE (bad pixel list)
  DRMS_Segment_t *segin1 = NULL;
  DRMS_Array_t *arrin1   = NULL;
  segin1    = drms_segment_lookupnum(data->records[0],1);
  arrin1    = drms_segment_read(segin1,type,&status);

  //RUN THE LIMBFINDER
  FID     = drms_getkey_int(data->records[0],"FID",&status);
  CROTA2  = drms_getkey_float(data->records[0],"CROTA2",&status);
  CDELT1  = drms_getkey_float(data->records[0],"CDELT1",&status);
  OBSVR   = drms_getkey_double(data->records[0],"OBS_VR",&status);
  printf("CDELT1= %f CROTA2= %f OBS_VR= %f \n",CDELT1,CROTA2,OBSVR);
  status  = limb_fit(data->records[0],image->data,&tempRSUN,&tempX0,&tempY0,4096,4096,0);
  RSUN_LF=(double)tempRSUN;
  X0_LF=(double)tempX0;
  Y0_LF=(double)tempY0;
  printf("X0_LF= %f Y0_LF= %f RSUN_LF= %f FID= %d \n",X0_LF,Y0_LF,RSUN_LF,FID);
  //CORRECT FOR FORMATION HEIGHT
  status  = heightformation(FID,OBSVR,&CDELT1,&RSUN_LF,&X0_LF,&Y0_LF,-CROTA2);
  printf("CDELT1= %f \n",CDELT1);
  //COPY KEYWORDS
  drms_copykeys(dataout->records[0],data->records[0],0,kDRMS_KeyClass_Explicit);
  status = drms_setkey_float(dataout->records[0],"CRPIX1",X0_LF+1);
  status = drms_setkey_float(dataout->records[0],"CRPIX2",Y0_LF+1);
  status = drms_setkey_float(dataout->records[0],"R_SUN",RSUN_LF);
  status = drms_setkey_float(dataout->records[0],"CDELT1",CDELT1);

  //WRITE IMAGE SEGMENT
  printf("Writing a record on the DRMS\n");
  segout = drms_segment_lookupnum(recout, 0);
  image->bzero = segout->bzero; image->bscale = segout->bscale; image->israw = 0;
  drms_segment_write(segout,image,0);
  //WRITE BAD PIXEL LIST SEGMENT
  printf("Writing a record on the DRMS\n");
  segout = drms_segment_lookupnum(recout, 1);
  image->bzero = segout->bzero; image->bscale = segout->bscale; image->israw = 0;
  drms_segment_write(segout,arrin1,0);
  
  //CLOSE RECORDS
  drms_close_records(dataout, DRMS_INSERT_RECORD); //insert the record in DRMS
 
  return(0);
  
}
