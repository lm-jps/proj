#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>
#include "module_flatfield.h"

int get_latest_dark(TIME t_0, int camera, const char *fname, const char *segname, float *offpoint, long long *recnum, int *focus) //function to retrive  latest dark records (the latest before t0)
{
#include "module_flatfield_const.h" 

  int status, nRec;
  int i,j;
  DRMS_RecordSet_t *ff_series;
  DRMS_Record_t *record_off;
  DRMS_Segment_t *segin  = NULL;
  DRMS_Array_t *arr_off;
  float *arr_data;

   char t_orig_str[256]="";
   sprint_ut(t_orig_str, 0.0); //one month before t_0
  
   char t0_str[256]="";
   sprint_ut(t0_str, t_0);

   char *camstr_side="[1]";
   char *camstr_front="[2]";

   char *camstr;
   if (camera == 1){camstr=camstr_side;}
   if (camera == 2){camstr=camstr_front;}


  char query_series[256]="";
  strcat(query_series, fname);
  strcat(query_series, camstr);
  strcat(query_series, "[");
  strcat(query_series, t_orig_str);
  strcat(query_series, "-");
  strcat(query_series, t0_str);
  strcat(query_series, "]");

  printf("query series %s\n", query_series);

  ff_series=drms_open_records(drms_env,query_series,&status);

  if (status == DRMS_SUCCESS && ff_series != NULL && ff_series->n != 0)
    {
      nRec=ff_series->n;

      record_off=ff_series->records[nRec-1];
      segin    = drms_segment_lookup(record_off, segname);
      arr_off = drms_segment_read(segin, segin->info->type, &status);
      arr_data=arr_off->data;
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) offpoint[j*nx+i]=arr_data[j*nx+i]; 

      *recnum=record_off->recnum;
      *focus = drms_getkey_int(record_off,keyfocusflat,&status);
      printf("record number %ld\n", *recnum);

    }
  else 
    {
      printf("can not find hmi  record\n");
      return 1;
    }

  return 0;
}


int get_latest_flat(TIME t_0, int camera, const char *fname, const char *segname, float *offpoint, long long *recnum, int *focus) //function to retrieve latest flatfield record
{
#include "module_flatfield_const.h" 

  int status, nRec;
  int i,j;
  DRMS_RecordSet_t *ff_series;
  DRMS_Record_t *record_off;
  DRMS_Segment_t *segin  = NULL;
  DRMS_Array_t *arr_off;
  float *arr_data;

   char t_orig_str[256]="";
   sprint_ut(t_orig_str, 0.0); //one month before t_0
  
   char t0_str[256]="";
   sprint_ut(t0_str, t_0);

   char *camstr_side="[1]";
   char *camstr_front="[2]";

   char *camstr;
   if (camera == 1){camstr=camstr_side;}
   if (camera == 2){camstr=camstr_front;}


  char query_series[256]="";
  strcat(query_series, fname);
  strcat(query_series, camstr);
  strcat(query_series, "[");
  strcat(query_series, t_orig_str);
  strcat(query_series, "-");
  strcat(query_series, t0_str);
  strcat(query_series, "][?PZTFLAG=1?]");

  printf("query series %s\n", query_series);

  ff_series=drms_open_records(drms_env,query_series,&status);

  if (status == DRMS_SUCCESS && ff_series != NULL && ff_series->n != 0)
    {
      nRec=ff_series->n;

      record_off=ff_series->records[nRec-1];
      segin    = drms_segment_lookup(record_off, segname);
      arr_off = drms_segment_read(segin, segin->info->type, &status);
      arr_data=arr_off->data;
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) offpoint[j*nx+i]=arr_data[j*nx+i]; 

      *recnum=record_off->recnum;
      *focus = drms_getkey_int(record_off,keyfocusflat,&status);
      printf("record number %ld\n", *recnum);

    }
  else 
    {
      printf("can not find hmi  record\n");
      return 1;
    }

  return 0;
}




int get_latest_bad(TIME t_0, int camera, const char *fname, const char *segname, short *bad, long long *recnum) //function to obtain latest permanent bad pixel list (latest before t0)
{
#include "module_flatfield_const.h" 

  int status, nRec;
  int i,j;
  DRMS_RecordSet_t *ff_series;
  DRMS_Record_t *record_off;
  DRMS_Segment_t *segin  = NULL;
  DRMS_Array_t *arr_off;
  int *bad_data;

   char t_orig_str[256]="";
   sprint_ut(t_orig_str, 0.0); //one month before t_0
  
   char t0_str[256]="";
   sprint_ut(t0_str, t_0);

   char *camstr_side="[1]";
   char *camstr_front="[2]";

   char *camstr;
   if (camera == 1){camstr=camstr_side;}
   if (camera == 2){camstr=camstr_front;}


  char query_series[256]="";
  strcat(query_series, fname);
  strcat(query_series, camstr);
  strcat(query_series, "[");
  strcat(query_series, t_orig_str);
  strcat(query_series, "-");
  strcat(query_series, t0_str);
  strcat(query_series, "]");

  printf("query series %s\n", query_series);

  ff_series=drms_open_records(drms_env,query_series,&status);

  if (status == DRMS_SUCCESS && ff_series != NULL && ff_series->n != 0)
    {
      nRec=ff_series->n;

      record_off=ff_series->records[nRec-1];
      segin    = drms_segment_lookup(record_off, segname);
      arr_off = drms_segment_read(segin, segin->info->type, &status);
      
      for (i=0; i<nx; ++i) for (j=0; j<ny; ++j) bad[j*nx+i]=1;  //initialize bad pixel array
      int nelem=arr_off->axis[0];
	
         
	  if (nelem > 0)
	    {
	      bad_data=arr_off->data;
	      for (i=0; i<nelem; ++i) bad[bad_data[i]]=0;
	    }

      *recnum=record_off->recnum;
      printf("bad pix record number %ld\n", *recnum);

    }
  else 
    {
      printf("can not find hmi record\n");
      return 1;
    }

  return 0;
}


int retrieve_offpoint(char *query, int camera, float *offpoint, long long *recnum, int *focus) //fucntion to retrieve lateset offpoint record
{
#include "module_flatfield_const.h" 
  DRMS_RecordSet_t *ff;
  int status, i, j, cam_id;
  DRMS_Record_t *record_off;
  DRMS_Segment_t *segin  = NULL;
  DRMS_Array_t *arr_off;
  float *arr_data;

  ff = drms_open_records(drms_env,query,&status);
  if (status == DRMS_SUCCESS && ff != NULL && ff->n != 0)
    {
      record_off=ff->records[0];
      cam_id = drms_getkey_int(record_off,keycamera,&status);
      *focus=drms_getkey_int(record_off,keyfocusflat,&status);

      if (cam_id != camera){printf("flatfield record inconsistent with camera\n"); return 1;}

      segin    = drms_segment_lookup(record_off, segmentname_offpoint);
      arr_off = drms_segment_read(segin, segin->info->type, &status);
      arr_data=arr_off->data;
      *recnum=record_off->recnum;
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) offpoint[j*nx+i]=arr_data[j*nx+i];
    }
else
  {
    printf("can not find reference offpoint flatfield\n");
    return 1;
  }
  return 0;

}



int get_flat(char *query, int camera, float *flatfield, float *offpoint,  
	     short *bad) //get flatfield by query link (FLAT_REC)
{

#include "module_flatfield_const.h" 

  DRMS_RecordSet_t *ff, *badrecs;
  DRMS_Record_t *record_flat;
  DRMS_Record_t *reclink_off, *reclink_dark, *reclink_bad;
  int status, status_dat;
  DRMS_Segment_t *segin  = NULL;
  DRMS_Array_t *arr_flat;
  DRMS_Array_t *arr_offpoint, *arr_bad;
  float *arr_data;
  int *bad_data;
  int i, j;
  int cam_id;

  for (i=0; i<nx; ++i) for (j=0; j<ny; ++j) bad[j*nx+i]=1;  //initialize bad pixel array

  ff = drms_open_records(drms_env,query,&status);
  status_dat=1;

if (status == DRMS_SUCCESS && ff != NULL && ff->n != 0)
	{

	  record_flat=ff->records[0];
	  segin    = drms_segment_lookup(record_flat, segmentname);
	  arr_flat = drms_segment_read(segin, segin->info->type, &status_dat);
	  arr_data=arr_flat->data;
	}

 
 if (status_dat ==0)
   {

	  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flatfield[j*nx+i]=arr_data[j*nx+i];


	  cam_id = drms_getkey_int(record_flat,keycamera,&status);

	  if (cam_id != camera){printf("flatfield record inconsistent with camera\n"); return 1;}

	  reclink_off=drms_link_follow(record_flat, linkoff, &status);;
	  if (reclink_off == NULL){printf("can not follow link to offpoint flatfield\n"); return 1;}
	  //recnum_off=reclink_off->recnum;

	  reclink_dark=drms_link_follow(record_flat, linkdark, &status);
	  if (reclink_dark == NULL){printf("can not follow link to dark image\n"); return 1;}
	  //recnum_dark=reclink_dark->recnum;


	  reclink_bad=drms_link_follow(record_flat, linkbad, &status);
	  if (reclink_bad == NULL){printf("can not follow link to bad pixel list\n"); return 1;}
	  //recnum_bad=reclink_bad->recnum;




	  //follow offpoint link
	  segin    = drms_segment_lookup(reclink_off, segmentname_offpoint);
	  arr_offpoint= drms_segment_read(segin, segin->info->type, &status);

	  arr_data=arr_offpoint->data;
	  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) offpoint[j*nx+i]=arr_data[j*nx+i];

	  //follow bad pixel link


	  segin    = drms_segment_lookup(reclink_bad, segmentname_badpix);
	  arr_bad= drms_segment_read(segin, segin->info->type, &status);
	  int nelem=arr_bad->axis[0];
	  bad_data=arr_bad->data;
         
	  if (nelem > 0)
	    {
	      for (i=0; i<nelem; ++i) bad[bad_data[i]]=0;
	    }

	  return 0;
   }
else
  {
	      printf("Flatfield %s not found\n", query);
	      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flatfield[j*nx+i]=1.0;
	      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) offpoint[j*nx+i]=1.0;

  //get latest bad data 
	      int *bad_data;
	      //built query string for bad pixel record
	      char query_bad[256]={""};
	      strcat(query_bad,filename_badpix);
	      strcat(query_bad,"[");
	      char ffnumb[1]={""};
	      sprintf(ffnumb, "%1.1d", camera);
	      strcat(query_bad, ffnumb);
	      strcat(query_bad, "]");

	      badrecs=drms_open_records(drms_env,query_bad,&status);
	      int nRec=badrecs->n;
	      if (nRec > 0)
		{
		  reclink_bad=badrecs->records[nRec-1];

		  segin    = drms_segment_lookup(reclink_bad, segmentname_badpix);
		  arr_bad= drms_segment_read(segin, segin->info->type, &status);
		  int nelem=arr_bad->axis[0];
		  bad_data=arr_bad->data;

		  if (nelem > 0)
		    {
		      for (i=0; i<nelem; ++i) bad[bad_data[i]]=0;
		    }

		 
		}
	      else 
		{
		  printf("can not find any bad pixel array for camera\n");
		  return 1;
		}
	      return 1;
  }


}



int read_flatfield_series(TIME t_0, int camera, float *flatfield, int *focus, TIME tobs_link[2], long long recnum[6],
			  struct rotpar *rot_cur) // read latest flatfield series and provide record number for links
{
#include "module_flatfield_const.h" 
  DRMS_RecordSet_t *ff_series;
 int i, j, k, nRec, status;
 DRMS_Record_t *record_flat_series, *record_clone;
 DRMS_Record_t *reclink_off, *reclink_dark, *reclink_bad;
 double tm1[4096];
 double tm2[4096];
 int version[4096];
 TIME t_start, t_stop_new;
 DRMS_Segment_t *segin  = NULL;
 DRMS_Array_t *arr;
 float *arr_data;

//built query string

 char *camstr_side="[1]";
 char *camstr_front="[2]";

 char *camstr;
 if (camera == 1){camstr=camstr_side;}
 if (camera == 2){camstr=camstr_front;}


  char t_orig_str[256]="";
  sprint_ut(t_orig_str, t_0-30.0*24.0*60.0*60.0);
  
  char t0_str[256]="";
  sprint_ut(t0_str, t_0);

  char query_series[256]="";
  strcat(query_series, filename_flatfield_series);
  strcat(query_series, camstr);
  strcat(query_series, "[");
  strcat(query_series, t_orig_str);
  strcat(query_series, "-");
  strcat(query_series, t0_str);
  strcat(query_series, "]");

  printf("query series %s\n", query_series);

  ff_series=drms_open_records(drms_env,query_series,&status);

  if (status == DRMS_SUCCESS && ff_series != NULL && ff_series->n != 0 && ff_series->n < 4096)
    {
    nRec=ff_series->n;

    for (k=0; k<nRec; ++k)
    {
      record_flat_series=ff_series->records[k];
      tm1[k]=drms_getkey_time(record_flat_series, keytstart,&status);
      tm2[k]=drms_getkey_time(record_flat_series, keytstop,&status);
      version[k]=drms_getkey_int(record_flat_series,keyversion,&status);
    }
    }
     else
      {
	printf("no flatfield series found\n");
	return 1;
      }


  for (k=0; k< nRec; ++k)
    {
      if (tm1[k] <= t_0)
	{
	  if (tm2[k] > t_0)
	    {
		    
	      if (version[k] > 0)
			{
			  t_stop_new=tm2[k];
			  t_start=tm1[k];
			  record_clone=ff_series->records[k];
			  rot_cur->flatfield_version=version[k];
			}
		      else 
			{
			  t_stop_new=0.0;
			  t_start=tm1[k];
			  record_clone=ff_series->records[k];
			  rot_cur->flatfield_version=0;

			}
		      }

		    }
		  
    }	    

  //get keyword and link information for record clone
  segin    = drms_segment_lookup(record_clone, segmentname);
  arr = drms_segment_read(segin, segin->info->type, &status);
  arr_data=arr->data;
  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flatfield[j*nx+i]=arr_data[j*nx+i];


  rot_cur->rotbad=drms_getkey_int(record_clone,keynewpix,&status);
  rot_cur->rotpairs=drms_getkey_int(record_clone,keynpairs,&status);
  rot_cur->rotcadence=drms_getkey_float(record_clone,keycadence,&status);

  *focus=drms_getkey_int(record_clone,keyfocusflat,&status);

  reclink_off=drms_link_follow(record_clone, linkoff, &status);;
  if (reclink_off == NULL){printf("can not follow link to offpoint flatfield\n"); return 1;}
  recnum[0]=reclink_off->recnum;

  reclink_dark=drms_link_follow(record_clone, linkdark, &status);
  if (reclink_dark == NULL){printf("can not follow link to dark image\n"); return 1;}
  recnum[1]=reclink_dark->recnum;


  reclink_bad=drms_link_follow(record_clone, linkbad, &status);
  if (reclink_bad == NULL){printf("can not follow link to bad pixel list\n"); return 1;}
  recnum[2]=reclink_bad->recnum;



  tobs_link[0]=t_start;
  tobs_link[1]= t_stop_new;

 return 0;
}


int write_rot_flatfield(char *filename_flatfield, DRMS_Array_t *arrout_new, int camera,
		     TIME t_0, int focus, struct rotpar rot_new)
{    
    
#include "module_flatfield_const.h"
 
  DRMS_Segment_t *segout = NULL;
  DRMS_RecordSet_t *dataout = NULL, *dataout_new=NULL;
  DRMS_Record_t *recout  = NULL;

 
  int cam_id;
  char *camera_str;

  TIME t_stop=t_0+24.0*60.0*60.0;

   if (camera == 1){ cam_id=1; camera_str=camera_str_side;}
   if (camera == 2){ cam_id=2; camera_str=camera_str_front;}


 

   int status, stat;

   dataout_new = drms_create_records(drms_env,1,filename_rot_flatfield,DRMS_PERMANENT,&stat);


 if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield);
	      return 1;
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Write rotational flatfield\n");
	     
	      recout = dataout_new->records[0];

	      status = drms_setkey_time(recout, keytstart, t_0);
	      status = drms_setkey_int(recout, keycamera, cam_id);

	      status=drms_setkey_int(recout, keyfocusflat, focus);
	      status=drms_setkey_time(recout, keytstop, t_stop); 
	      drms_setkey_string(recout, keyinstrument, camera_str);

	      status=drms_setkey_int(recout, "COPYFLAG", 0);
	      status=drms_setkey_int(recout, keynewpix, rot_new.rotbad);
	      status=drms_setkey_int(recout, keynpairs, rot_new.rotpairs);
	      status=drms_setkey_float(recout, keycadence, rot_new.rotcadence);

	      drms_keyword_setdate(recout);


	      segout = drms_segment_lookup(recout, segmentname);
	      if (segout == NULL){printf("could not find segment\n"); return 1;}
	
              //arr_flat->bzero=segout->bzero;  //for file compression
	      //arr_flat->bscale=segout->bscale;
	      //arr_flat->israw=0;

	      status=drms_segment_write(segout, arrout_new, 0);        //write the file containing the data
	      if (status != DRMS_SUCCESS){printf("could not write segment\n"); return 1;}
	      
	    }

	  if (dataout_new != NULL) status=drms_close_records(dataout_new, DRMS_INSERT_RECORD); //insert the record in DRMS
	  
	  if (status == DRMS_SUCCESS)
	    {
	      printf("done\n");
	      return 0;
	    }
	  else 
	    {
	      return 1;
	    }

}



int write_flatfields(char *filename_flatfield, DRMS_Array_t *arr_flat, DRMS_Array_t *arrout_new, int camera,
                     long long recnum[6], TIME tobs_link[2], TIME t_0, int focus, int focusclone,
                     struct rotpar rot_new, 
                     struct rotpar rot_cur) //write out rotational flatfields
{    
    
#include "module_flatfield_const.h" 
  TIME t_start=tobs_link[0];
  TIME t_stop_new=tobs_link[1];

  if (t_stop_new == 0.0) t_stop_new=t_0+365.0*24.0*60.0*60.0;


  long long recnum_off=recnum[0];
  long long recnum_dark=recnum[1];
  long long recnum_bad=recnum[2];
  long long recnum2_off=recnum[3];
  long long recnum2_dark=recnum[4];
  long long recnum2_bad=recnum[5];


   DRMS_Segment_t *segout = NULL;
   DRMS_RecordSet_t *dataout = NULL, *dataout_new=NULL;
   DRMS_Record_t *recout  = NULL;

 
   int cam_id;
   char *camera_str;

   if (camera == 1){ cam_id=1; camera_str=camera_str_side;}
   if (camera == 2){ cam_id=2; camera_str=camera_str_front;}


 

   int status, stat;

 

   int flatvers, flatvers_back, flatvers_forward;

 /////////////////////////////////////////////////////////////////////////////////
 ///write flatfield with new T_STOP /////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////

   // dataout=drms_clone_records(ff_last, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &stat);
   

   flatvers=rot_cur.flatfield_version;
   if (flatvers == 0){flatvers_back=1; flatvers_forward=0;}
   if (flatvers > 0){flatvers_back=flatvers; flatvers_forward=flatvers+1;}

   dataout = drms_create_records(drms_env,1,filename_flatfield,DRMS_PERMANENT,&stat);
   
  if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield);
	      return 1;
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Overwrite current flatfield with definite T_STOP\n");
	      printf("Writing a record on the DRMS for the series %s\n",filename_flatfield);

	      //front camera
	     
	      recout = dataout->records[0];
	 
	      status = drms_setkey_time(recout, keytstart, t_start);
	      status=drms_setkey_int(recout, keycamera, cam_id);

	      status=drms_setkey_int(recout, keyfocusflat, focusclone);
	      status=drms_setkey_time(recout, keytstop, t_0); 
	      status=drms_setkey_int(recout, keyversion, flatvers_back);

   //write flatfield version
	      status = drms_setkey_string(recout, keyinstrument, camera_str);

	      status=drms_setkey_int(recout, keynewpix, rot_cur.rotbad);
	      status=drms_setkey_int(recout, keynpairs, rot_cur.rotpairs);
	      status=drms_setkey_float(recout, keycadence, rot_cur.rotcadence);
	      printf("set links\n");
	
	      drms_keyword_setdate(recout);
              /******************************************/
	      /*set link                                */	 
              /******************************************/
	        {
		  //	 const char *linkname = "offpoint_flatfield";
	    
	         status=drms_setlink_static(recout, linkoff, recnum_off);
		 printf("offpoint: link to recum%ld\n", recnum_off);
	       	 
	          }
	          {
		    //	   const char *linkname = "bias_dark";
	       	  
		   status=drms_setlink_static(recout, linkdark, recnum_dark);
		   printf("dark: link to recum%ld\n", recnum_dark);
	      	 	       	 
	          }
	          {
		    //   const char *linkname = "perm_bad_pixel";
	        
	      	   status=drms_setlink_static(recout, linkbad, recnum_bad);
		   printf("badpix: link to recum%ld\n", recnum_bad);
	      	  
	       }

	
     
		  segout = drms_segment_lookup(recout, segmentname);

		  //arr_flat->bzero=segout->bzero;  //for file compression
		  //arr_flat->bscale=segout->bscale;
		  //arr_flat->israw=0;
	   
	         status=drms_segment_write(segout, arr_flat, 0);        //write the file containing the data
	      

	    }


	  drms_close_record(recout, DRMS_INSERT_RECORD); //insert the record in DRMS


    
    
    

    

	  ///////////////////////////////////////////////////////////////////////
	  //////////write new updated flatfield ////////////////////////////////
	  /////////////////////////////////////////////////////////////////////
	 
	 
    
 dataout_new = drms_create_records(drms_env,1,filename_flatfield,DRMS_PERMANENT,&stat);


 if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield);
	      return 1;
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Write new flatfield\n");
	     
	      recout = dataout_new->records[0];

	      status = drms_setkey_time(recout, keytstart, t_0);
	      status = drms_setkey_int(recout, keycamera, cam_id);

	      status=drms_setkey_int(recout, keyfocusflat, focus);
	      status=drms_setkey_time(recout, keytstop, t_stop_new); 
	      status=drms_setkey_int(recout, keyversion, flatvers_forward);    //update FLAT_FINAL
	      drms_setkey_string(recout, keyinstrument, camera_str);

	      status=drms_setkey_int(recout, keynewpix, rot_new.rotbad);
	      status=drms_setkey_int(recout, keynpairs, rot_new.rotpairs);
	      status=drms_setkey_float(recout, keycadence, rot_new.rotcadence);

	      drms_keyword_setdate(recout);
              /******************************************/
	      /*set link                                */	 
              /******************************************/
	      {
		    	    
		   status=drms_setlink_static(recout, linkoff, recnum2_off);
	    	   
		   if (status != DRMS_SUCCESS){printf("could not set like to offpoint flatfield\n"); return 1;}
		   else {printf("offpoint: link to recum%ld\n", recnum2_off);}
	      }
	      {
	
		   status=drms_setlink_static(recout, linkdark, recnum2_dark);
	      	   
		   if (status != DRMS_SUCCESS){printf("could not set like to dark\n"); return 1;}
		   else {printf("dark: link to recum%ld\n", recnum2_dark);}
	      }

	      {
		       
		   status=drms_setlink_static(recout, linkbad, recnum2_bad);
		   if (status != DRMS_SUCCESS){printf("could not set like to bad pixel list\n"); return 1;}
		   else {printf("badpix: link to recum%ld\n", recnum2_bad);}
	      }


	      segout = drms_segment_lookup(recout, segmentname);
	      if (segout == NULL){printf("could not find segment\n"); return 1;}
	      //arr_flat->bzero=segout->bzero;  //for file compression
	      //arr_flat->bscale=segout->bscale;
	      //arr_flat->israw=0;

	      status=drms_segment_write(segout, arrout_new, 0);        //write the file containing the data
	      if (status != DRMS_SUCCESS){printf("could not write segment\n"); return 1;}
	      
	    }
    



	  status=drms_close_records(dataout_new, DRMS_INSERT_RECORD); //insert the record in DRMS
	  
	  if (status == DRMS_SUCCESS)
	    {
	      printf("done\n");
	      return 0;
	    }
	  else 
	    {
	      return 1;
	    }

}



    
int write_flatfield_fid(DRMS_Array_t *arrout_new, int camera, long long recnum_offpoint,
			TIME t_0, int focus,int fid, 
			struct rotpar rot_new) //write out per-FID flatfields
{    
    
#include "module_flatfield_const.h"  
 


   DRMS_Segment_t *segout = NULL;
   DRMS_RecordSet_t *dataout = NULL, *dataout_new=NULL;
   DRMS_Record_t *recout  = NULL;

 
   int cam_id;
   char *camera_str;

   if (camera == 1){ cam_id=1; camera_str=camera_str_side;}
   if (camera == 2){ cam_id=2; camera_str=camera_str_front;}

 
   int status, stat;

 
   const int FLAG_PRELIMINARY=0;
   const int FLAG_FINAL=1;

   char query_flat[256]="";
   strcat(query_flat, filename_offpoint);
   strcat(query_flat, "[:#");
   char ffnumb[10]={""};
   sprintf(ffnumb, "%1.1ld", recnum_offpoint);
   strcat(query_flat, ffnumb);
   strcat(query_flat, "]");

   printf("query offpoint flatfield: %s", query_flat);

	  ///////////////////////////////////////////////////////////////////////
	  //////////write new updated flatfield ////////////////////////////////
	  /////////////////////////////////////////////////////////////////////
	 
	 
    
 dataout_new = drms_create_records(drms_env,1,filename_flatfield_fid,DRMS_PERMANENT,&stat);


 if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield_fid);
	      return 1;
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Write new flatfield with fid %d\n", fid);
	     
	      recout = dataout_new->records[0];

	      status=drms_setkey_time(recout, keyday, t_0);
	      //	      status = drms_setkey_time(recout, keytstart, t_start);
	      status=drms_setkey_int(recout, keycamera, cam_id);
	      status=drms_setkey_int(recout, fidkey, fid);
	      //status=drms_setkey_int(recout, fsnskey, fsns);

	      status=drms_setkey_int(recout, keyfocusflat, focus);
	      //	      status=drms_setkey_time(recout, keytstop, t_stop); 
	      //status=drms_setkey_int(recout, keyversion, FLAG_PRELIMINARY);    //update FLAT_FINAL
	      drms_setkey_string(recout, keyinstrument, camera_str);

	      //status=drms_setkey_int(recout, keynewpix, rot_new.rotbad);
	      status=drms_setkey_int(recout, keynpairs, rot_new.rotpairs);
	      status=drms_setkey_float(recout, keycadence, rot_new.rotcadence);
	      //status=drms_setkey_int(recout, recnumkey, (int)recnum_off);
	      status=drms_setkey_string(recout, querykey, query_flat);
	      

	      drms_keyword_setdate(recout);
            


	      segout = drms_segment_lookup(recout, segmentname);
	      if (segout == NULL){printf("could not find segment\n"); return 1;}

              //arrout_new->bzero=segout->bzero;  //for file compression
	      //arrout_new->bscale=segout->bscale;
	      //arrout_new->israw=0;
	   
	      status=drms_segment_write(segout, arrout_new, 0);        //write the file containing the data
	      if (status != DRMS_SUCCESS){printf("could not write segment\n"); return 1;}
	    }
    



	  status=drms_close_records(dataout_new, DRMS_INSERT_RECORD); //insert the record in DRMS

	  if (status == DRMS_SUCCESS)
	    {
	      printf("done\n");
	      return 0;
	    }
	  else 
	    {
	      return 1;
	    }
}




