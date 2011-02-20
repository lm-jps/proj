 /*
 * cosmic_ray_post - Post processing of cosmic ray series for false positives
 *
 */

/**
\defgroup module_flatfield module_flatfield - derive line-of-sight observables

\par Synopsis
\code
cosmic_ray_post input_series= camera= datum= hour=
\endcode

\details

cosmic_ray_post is post-processing the cosmic ray series (hmi.cosmic_rays) for false positives. The module processes 2 hours at a time, starting 
at the time given by the datum and hour argument

par Options

\par Mandatory arguments:
\li \c input_series="string" where string is the series name of the intermediate flatfield (su_production.flatfield_fid)
\li \c camera=cam,  side camera: cam=1, front camera: cam=2
\li \c datum="date" date="yyyy.mm.dd" TAI-day for which the intermediate flatfields have been calculated. End of TAI-day datum is T_START of updated flatfield
\li \c hour=hour

\par Examples

\b Example 1:

\code
cosmic_ray_post input_series="su_production.flatfield_fid" camera=2 datum="2010.10.09" hour=2
\endcode

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>
#include </home/jsoc/include/fftw3.h>
#include <omp.h>
#include <fresize.h>
#include <module_flatfield.h>

char *module_name  = "module_flatfield_combine";    //name of the module

#define kRecSetIn      "input_series" //name of the series containing the input filtergrams
#define datumn  "datum"
#define cameran "camera"
#define fsnname      "fsn"
#define fsnf_name "fsn_first"
#define fsnl_name "fsn_last"

#define minval(x,y) (((x) < (y)) ? (x) : (y))                  

ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn, "",  "Input data series."},
     {ARG_STRING, datumn, "yyyy.mm.dd", "Datum string"},
     {ARG_INT, "hour", "00", "Hour"},
     {ARG_INT, cameran, 0, "Camera"},
     {ARG_INT, fsnf_name, "0"},
     {ARG_INT, fsnl_name, "2147483647"},
     {ARG_END}
};



void printtime();
int is_element(int value, int *array, int n);

int DoIt(void)
{

#include "module_flatfield_const.h"

  //*********************
  //read input parameters

  const char *inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL); //cmdparams is defined in jsoc_main.h
  const char *datum =  cmdparams_get_str(&cmdparams, datumn, NULL);
 
 int cameraint = cmdparams_get_int(&cmdparams, cameran, NULL);
 int hour= cmdparams_get_int(&cmdparams, "hour", NULL);
 int fsn_first=cmdparams_get_int(&cmdparams, fsnf_name, NULL);
 int fsn_last=cmdparams_get_int(&cmdparams, fsnl_name, NULL);


 //**********************
 //define variables


  int  status= DRMS_SUCCESS, stat=DRMS_SUCCESS;

  DRMS_Segment_t *seg, *seg_val, *seg_sig;
  DRMS_Segment_t *segout, *segout_val, *segout_sig;
  
  DRMS_Type_t type_time = DRMS_TYPE_TIME;
  DRMS_Type_t type_int = DRMS_TYPE_INT;
  DRMS_Type_t type_float = DRMS_TYPE_FLOAT;
  DRMS_Type_t type_double = DRMS_TYPE_DOUBLE;

  DRMS_RecordSet_t *data, *dataout;
  DRMS_Record_t  *recout, *rec0;
  DRMS_Array_t *arrout, *arrout_val, *arrout_sig;

  int *cosmic_ray_data;
  float *val_data, *sig_data;

  DRMS_Type_Value_t  keyvalue_time;

  int ct,ct_prev,ct_next,ctr;
  float *sig, *lev;
  int *hits;
  int axisbad[1];

  int i, j, k, c; //loop variables


  int axisout[2]={nx,ny};

  
  //**********************
  //*CHECK WHETHER THE FLATFIELD OUTPUT SERIES EXISTS                                                                    */

    
       drms_series_exists(drms_env, filename_cosmic2_out, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",filename_cosmic2_out);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
  

      printf("START!\n");
      printtime();

      //*****************************
      //build query string

      char fnumb[2]={""};
      char ffnumb[2]={""};
      char timefirst[256]="";

      strcat(timefirst, datum);
      strcat(timefirst, "_");
      sprintf(fnumb, "%2.2d", hour);
      strcat(timefirst, fnumb);
      strcat(timefirst, ":00:00.00_TAI");
 
 

      char query0[256]="";
      TIME tfirst, tlast;

      printf("fsnfs %d %d\n", fsn_first, fsn_last);


 if (fsn_first == 0 || fsn_last == 2147483647)
   {

 char timelast[256]="";
 strcat(timelast, datum);
 strcat(timelast, "_");
 sprintf(ffnumb, "%2.2d", hour+1);
 strcat(timelast, ffnumb);
 strcat(timelast, ":59:59.99_TAI");


     tfirst=sscan_time(timefirst);
     tlast=sscan_time(timelast);

     char datefirst[256]="";
     sprint_ut(datefirst, tfirst-5.0);

 char datelast[256]="";
 sprint_ut(datelast, tlast+5.0);


 strcat(query0, inRecQuery);
 strcat(query0, "[");
 strcat(query0, datefirst);
 strcat(query0, "-");
 strcat(query0, datelast);
 strcat(query0, "][?CAMERA=");
 sprintf(fnumb, "%1.1d", cameraint);
 strcat(query0, fnumb);
 strcat(query0, "?]");

  printf("query string time: %s\n", query0);
    
  }


  if (fsn_first != 0 && fsn_last != 2147483647)
    {

      tfirst=0;
      tlast=3881520000.0;

      strcat(query0, inRecQuery);
      strcat(query0, "[][");
       char fsnf[10]={""};
       sprintf(fsnf, "%d", fsn_first-2);
       strcat(query0, fsnf);
      strcat(query0, "-");
      char fsnl[10]={""};
      sprintf(fsnl, "%d", fsn_last+2);
       strcat(query0, fsnl);
       strcat(query0, "][?CAMERA=");
      char fnumb[2]={""};
     sprintf(fnumb, "%1.1d", cameraint);
     strcat(query0, fnumb);
      strcat(query0, "?]");
      printf("query string fsn: %s\n", query0);  
    }

  
  //****************************
  //open records

 
  data     = drms_open_records(drms_env,query0,&stat);

  if (data == NULL){printf("can not open records\n"); exit(EXIT_FAILURE);}

  drms_stage_records(data, 1, 0);

  int nRecs=data->n;
  printf("number of records %d\n", nRecs);
  printtime();


  //******************************
  //read out keywords


if (stat == DRMS_SUCCESS && data != NULL && nRecs > 0)
  {
  int keyvalue_fid[nRecs];
  int keyvalue_cam[nRecs];
  DRMS_Record_t *rec0[nRecs];
  int  *cosmic_rays[nRecs];
  int keyvalue_recnum[nRecs];
  DRMS_Array_t *arr_cosmic[nRecs], *arr_level[nRecs], *arr_sig[nRecs];
  short wave[nRecs];
  short pol[nRecs];
  int count[nRecs];
  double time_fl[nRecs];
  int keyvalue_fsn[nRecs];
  int keyvalue_flag[nRecs];
  float keyvalue_factor[nRecs];
  float *significance[nRecs], *level[nRecs];
  int new_count[nRecs];
  int npairs[nRecs];


  printtime();
  printf("create_records ...");
  dataout = drms_create_records(drms_env,nRecs,filename_cosmic2_out,DRMS_PERMANENT,&stat);
  printf("done\n");
  printtime();

  printf("begin keyword  reading loop\n");
for (k=0; k<nRecs; ++k)
  {
    rec0[k]=data->records[k];

   
   
    count[k]= drms_getkey_int(rec0[k],keycount,&status);
    keyvalue_fid[k] = drms_getkey_int(rec0[k],fidkey,&status);
    keyvalue_cam[k] = drms_getkey_int(rec0[k],keycamera,&status);
    keyvalue_fsn[k]= drms_getkey_int(rec0[k],keyfsn,&status);
    keyvalue_flag[k] = drms_getkey_int(rec0[k],keyexmax,&status);
    keyvalue_factor[k] = drms_getkey_float(rec0[k],keylimit,&status);

    keyvalue_time=drms_getkey(rec0[k], keytobs, &type_time, &status);   ;   //read "T_OBS" of filtergram
    time_fl[k] = keyvalue_time.time_val;

    wave[k]=((keyvalue_fid[k]-10000)/10-5)/2;
    pol[k]=(keyvalue_fid[k] % 10);    
    
    keyvalue_recnum[k]=drms_getkey_int(rec0[k], recnumkey,&status);

      }
 printtime();

 

 //***********************************
 //data reading loop


printf("begin data  reading loop\n");

 for (k=0; k<nRecs; ++k)
   {

     seg = drms_segment_lookup(rec0[k], segmentname_cosmic);
     seg_val=drms_segment_lookup(rec0[k], segmentname_val);
     seg_sig=drms_segment_lookup(rec0[k], segmentname_sig);
	    
     arr_cosmic[k]= drms_segment_read(seg, type_int, &status);
     arr_level[k]=drms_segment_read(seg_val, type_float, &status);
     arr_sig[k]=drms_segment_read(seg_sig, type_float, &status);


     cosmic_rays[k] = arr_cosmic[k]->data;
     significance[k]=arr_sig[k]->data;
     level[k]=arr_level[k]->data;

   }

 printtime();

 //******************************
 //test loop

  int elem_prev, elem_next;

  int *cosmic_new;
  float *val_new;
  float *sig_new;

  int counter=0;
  int falsecounter=0, goodcounter=0;

  cosmic_new=(int *)(malloc(limit_cosmic*sizeof(int)));
  val_new=(float *)(malloc(limit_cosmic*sizeof(float)));
  sig_new=(float *)(malloc(limit_cosmic*sizeof(float)));




printf("begin test loop\n");

 for (k=0; k<nRecs; ++k)
   {
     

     if (time_fl[k] >= tfirst && time_fl[k] <= tlast  && keyvalue_fsn[k] >= fsn_first && keyvalue_fsn[k] <= fsn_last)
       {

	 hits=cosmic_rays[k];
	 sig=significance[k];
	 lev=level[k];

	 ct=minval(count[k], limit_cosmic);

	 if (count[k] > 0 && k >=1 && k < (nRecs-1)) //at least one count, not first or last record -> do weeding
	   {
	     
	     ct_prev=minval(count[k-1], limit_cosmic);
	     ct_next=minval(count[k+1], limit_cosmic);
  

	     ctr=0;

	     for (c=0; c<ct; ++c)
	       {
	    
	     	 ++counter;
		 elem_prev=is_element(hits[c], cosmic_rays[k-1], ct_prev);
		 elem_next=is_element(hits[c], cosmic_rays[k+1], ct_next);

		 if (elem_prev == 0 && elem_next == 0)
		   {
		     
		     cosmic_new[ctr]=hits[c];
		     val_new[ctr]=lev[c];
		     sig_new[ctr]=sig[c];
		     ++ctr;
		     ++goodcounter;
		   }
		 else
		   {
		     ++falsecounter;
		   }
	       }
	    	       
	     new_count[k]=ctr;
	   }
	 else  //else, just copy
       {
	 ctr=count[k];

	 for (c=0; c<ct; ++c)
	   {
	    		
	     cosmic_new[c]=hits[c];
	     val_new[c]=lev[c];
	     sig_new[c]=sig[c];                                  
	    
	   }
	       
	 new_count[k]=ctr;
       }

	 //*************************
	 //set keywords
     

     recout = dataout->records[k];

     status=drms_setkey_int(recout, keyfsn, keyvalue_fsn[k]);
     status=drms_setkey_time(recout, keytobs, time_fl[k]);
     
     status=drms_setkey_int(recout, fidkey, keyvalue_fid[k]);
     status=drms_setkey_int(recout, keycamera, keyvalue_cam[k]);
     status=drms_setkey_int(recout, keyexmax, keyvalue_flag[k]);
     status=drms_setkey_float(recout, keylimit, keyvalue_factor[k]); 
	
     drms_keyword_setdate(recout);

     status=drms_setkey_int(recout, keycount, new_count[k]);
    
     if (keyvalue_cam[k] == 2) status=drms_setkey_string(recout, keyinstrument, camera_str_front);
     if (keyvalue_cam[k] == 1) status=drms_setkey_string(recout, keyinstrument, camera_str_side);
	    
   	    
	 
     //************************************
     //write out data segments


     if (new_count[k] >= 0)
       {

	 segout = drms_segment_lookup(recout, segmentname_cosmic);
	 segout_val=drms_segment_lookup(recout, segmentname_val);
	 segout_sig=drms_segment_lookup(recout, segmentname_sig);

     axisbad[0]=new_count[k];
     arrout=drms_array_create(type_int,1,axisbad,NULL,&status);
     cosmic_ray_data=arrout->data;
	  	 
     arrout_val=drms_array_create(type_float,1,axisbad,NULL,&status);
     val_data=arrout_val->data;

     arrout_sig=drms_array_create(type_float,1,axisbad,NULL,&status);
     sig_data=arrout_sig->data;


     for (i=0; i<new_count[k]; ++i){cosmic_ray_data[i]=cosmic_new[i]; val_data[i]=val_new[i]; sig_data[i]=sig_new[i];} //copy new data array

       status=drms_segment_write(segout, arrout, 0);
       status=drms_segment_write(segout_val, arrout_val, 0);
       status=drms_segment_write(segout_sig, arrout_sig, 0);
	   
       drms_free_array(arrout);
       drms_free_array(arrout_val);
       drms_free_array(arrout_sig);
       }
       }
    
   }
 
 //****************
 //free arrays

 for (k=0; k<nRecs; ++k)
   {
     drms_free_array(arr_level[k]);
     drms_free_array(arr_cosmic[k]);
     drms_free_array(arr_sig[k]);
   }

   
 drms_close_records(dataout, DRMS_INSERT_RECORD);


 drms_close_records(data,DRMS_FREE_RECORD);

 free(cosmic_new);
 free(sig_new);
 free(val_new);


 printf("correctly identified cosmic rays %d out of %d %f\n", goodcounter, (goodcounter+falsecounter), (float)goodcounter/(float)(falsecounter+goodcounter));
 
  }
 else 
    {
      printf("No data records found\n");
    }     

 
   



 printtime();

 printf("COMLETED!\n");
 return 0;

}





int is_element(int value, int *array, int n)
{

  int isel, ctr;

  ctr=0;
  isel=0;
  while (ctr < n && isel == 0)
    {
      isel=(value == array[ctr]);
      ++ctr;
    }

  return isel;
}


void printtime()	
// print time
{
  time_t timer, timerd;
  char *timestring;
  int i;

	timerd=time(&timer);
	timestring=ctime(&timer);
        for (i=0; i<24; ++i) printf("%c", *(timestring+i));
	printf("\n");
}



