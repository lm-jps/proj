#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>
#include "module_flatfield.h"




int get_flatfields(int camera, TIME t_0, int focus, 
                   float *flatfield, float *offpoint, 
                   short *bad, long long recnum[6], TIME tobs_link[2],
		   struct rotpar *rot_cur)

{
  printf("get_flatfield\n");
#include "module_flatfield_const.h" // !!changed to test       

  DRMS_Segment_t *segin  = NULL;
    
  DRMS_RecordSet_t *ff;
   DRMS_Record_t *record_flat;
   DRMS_Record_t *reclink_off, *reclink_dark, *reclink_bad;

   DRMS_Array_t *arr_flat;
   float *arr_data;

   int i,j;


   int version,fvers,maxversion;
   int cam_id;

   TIME t_start, t_stop; 
   //TIME t_obs_offpoint, t_obs_dark, t_obs_bad;
 
   long long recnum_off,recnum_dark,recnum_bad;
   long long recnum2_off,recnum2_dark;

   DRMS_Type_t type_time      = DRMS_TYPE_TIME;

   int status;

   DRMS_Type_t type=DRMS_TYPE_FLOAT;

 
   //const char *key1 = "T_STOP";
   //const char *key2 = "HMI_SEQ_ID_FOCUS";
   //const char *key3 = "FLATFIELD_VERSION";

   //const char *key4 = "ROTF_FLATFIELD";
   //const char *key5 = "ROTF_N_PAIRS";
   //const char *key6 = "ROTF_CADENCE";

   //const char *keylink1 = "T_OBS_OFFPOINT";
   //const char *keylink2 = "T_OBS_DARK";
   //const char *keylink3 = "T_OBS_BADPIX"; 

   //const char *linkoff = "offpoint_flatfield";
   //const char *linkdark = "bias_dark";
   //const char *linkbad = "perm_bad_pixel";

   
     char *camstr_side="[1][";
     char *camstr_front="[2][";

     char *camstr;
   if (camera == 1){camstr=camstr_side;}
   if (camera == 2){camstr=camstr_front;}


  for (i=0; i<nx; ++i) for (j=0; j<ny; ++j) bad[j*nx+i]=1;  //initialize bad pixel array
 
  char t_orig_str[256];
  sprint_ut(t_orig_str, 0.0);
  
  

  
    char t0_str[256];
    sprint_ut(t0_str, t_0-5.0);

      {
   char query[256]={""};
   strcat(query,filename_flatfield);
   strcat(query,camstr);
   // strcat(query,"[2]["); //look for most recent flatfield
   strcat(query,t_orig_str);
   strcat(query,"-");
   strcat(query,t0_str);
   strcat(query,"]");
   //char ffnumb[2]={""};
   //sprintf(ffnumb, "%2.2d", focus);
   //strcat(query, ffnumb);
   //strcat(query, "]");
   printf("query string %s\n", query);
      

   ff = drms_open_records(drms_env,query,&status);
      }

 
      //*******************************************************************
      //query for most recent flatfield
      //*******************************************************************     
     


      char date_string[256];
      if (status == DRMS_SUCCESS && ff != NULL)
	{
	  int nRec=ff->n;

	  if (nRec > 0){
	  record_flat = ff->records[nRec-1];

	  //get max flatfield version
	  maxversion=0;
	  for (i=0; i<nRec; ++i){fvers=drms_getkey_int(ff->records[i], keyversion,&status); if (fvers > maxversion) maxversion=fvers;}
	  printf("maxversion %d %d\n", maxversion, nRec);

	  segin    = drms_segment_lookup(record_flat, segmentname);
	  printf("%ld\n", segin); //debug
	  arr_flat = drms_segment_read(segin, segin->info->type, &status);

	  printf("%ld\n", arr_flat); //debug

	  type=segin->info->type; // read type of flatfield

	  arr_data=arr_flat->data;
	
	  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flatfield[j*nx+i]=arr_data[j*nx+i];   //read data
	  printf("flatfield read\n"); 
	  cam_id = drms_getkey_int(record_flat,keycamera,&status);  //read keywords
	  t_start = drms_getkey_time(record_flat,keytstart,&status);
	  t_stop = drms_getkey_time(record_flat,keytstop,&status);
	  rot_cur->rotbad=drms_getkey_int(record_flat,keynewpix,&status);
	  rot_cur->rotpairs=drms_getkey_int(record_flat,keynpairs,&status);
	  rot_cur->rotcadence=drms_getkey_float(record_flat,keycadence,&status);
      
       //  t_obs_offpoint= drms_getkey_time(record_flat,keylink1,&status);
       //  t_obs_dark=drms_getkey_time(record_flat,keylink2,&status);
       //  t_obs_bad=drms_getkey_time(record_flat,keylink3,&status);
  
       reclink_off=drms_link_follow(record_flat, linkoff, &status);
       if (reclink_off == NULL){printf("can not follow link to offpoint flatfield\n"); exit(EXIT_FAILURE);}
       recnum_off=reclink_off->recnum;

       reclink_dark=drms_link_follow(record_flat, linkdark, &status);
       if (reclink_dark == NULL){printf("can not follow link to dark image\n"); exit(EXIT_FAILURE);}
       recnum_dark=reclink_dark->recnum;


       reclink_bad=drms_link_follow(record_flat, linkbad, &status);
       if (reclink_bad == NULL){printf("can not follow link to bad pixel list\n"); exit(EXIT_FAILURE);}
       recnum_bad=reclink_bad->recnum;

       //check most recent final flatfield for version number 
       version=drms_getkey_int(record_flat,keyversion,&status);
       printf("raw version %d\n", version);

       if (version > 0)
	 {
	   rot_cur->flatfield_version=version;
	 }
       else 
	 {
	   rot_cur->flatfield_version=maxversion;
	 }
	  	       
       printf("version after %d \n", rot_cur->flatfield_version);

       sprint_ut(date_string, t_start);
       printf("Most recent flatfield read. T_START: %s\n", date_string);

	  }
	  else
	    {
	  printf("Most recent flatfield not found\n");
	  exit(EXIT_FAILURE);
	    }
	}
      else 
	{
	  
	  printf("Most recent flatfield not found\n");
	  exit(EXIT_FAILURE);
	}
  





   DRMS_Record_t *record_offpoint;
   DRMS_Array_t *arr_offpoint;
   TIME t_obs2_offpoint;
  
   
   //////////////////////////////////////////////////
   // query for latest offpoint flatfield////////////
   //////////////////////////////////////////////////


   //   {
   //  char query_offpoint[256]={""}; 
   
   // strcat(query_offpoint,filename_offpoint);
   // strcat(query_offpoint,camstr);
   //strcat(query_offpoint,"[2][");
   // strcat(query_offpoint,t_orig_str);
   // strcat(query_offpoint, "-");
   //strcat(query_offpoint,t0_str);
   //strcat(query_offpoint,"]");
   //char ffnumb[2]={""};
   //sprintf(ffnumb, "%2.2d", focus);
   //strcat(query_offpoint, ffnumb);
   //strcat(query_offpoint, "]");

   //  printf("query offp %s\n", query_offpoint);
   

   //ff_offpoint = drms_open_records(drms_env,query_offpoint,&status);
   // }





   record_offpoint=reclink_off;
   recnum2_off=record_offpoint->recnum;

   segin    = drms_segment_lookup(record_offpoint, "offpoint_flatfield");
   arr_offpoint= drms_segment_read(segin, segin->info->type, &status);

   arr_data=arr_offpoint->data;
   for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) offpoint[j*nx+i]=arr_data[j*nx+i];

   t_obs2_offpoint=drms_getkey_time(record_offpoint,keytobs,&status);

   sprint_ut(date_string, t_obs2_offpoint);
   printf("Most recent offpoint flatfield read. T_OBS: %s %f\n", date_string, offpoint[2048*4096+2048]);



 //query for latest dark



 DRMS_Record_t *record_dark;
 TIME  t_obs2_dark;

 // {
 //  char query_dark[256]={""}; 
 //   strcat(query_dark,filename_dark);
 //   strcat(query_dark,camstr);
 //   strcat(query_dark,t_orig_str);
 //   strcat(query_dark, "-");
 //   strcat(query_dark,t0_str);
 //   strcat(query_dark,"]");
   
 //   ff_dark = drms_open_records(drms_env,query_dark,&status);
 
 //}



    record_dark=reclink_dark;
  recnum2_dark=record_dark->recnum;

  t_obs2_dark=drms_getkey_time(record_dark,keytobs,&status);

  sprint_ut(date_string, t_obs2_dark);
  printf("Most recent dark read. T_OBS: %s \n", date_string);

     


//read for latest bad pixel list


 DRMS_Record_t *record_bad;
 DRMS_Array_t *arr_bad;
 TIME t_obs2_bad;
 long long recnum2_bad;

 
 
 // char query_bad[256]={""}; 
 // strcat(query_bad,filename_badpix); 
 //   strcat(query_bad,camstr);
 //   strcat(query_bad,t_orig_str);
 //   strcat(query_bad, "-");
 //   strcat(query_bad,t0_str);
 //   strcat(query_bad,"]");
 //   printf("query bad %s\n", query_bad);

   // ff_bad = drms_open_records(drms_env,query_bad,&status);
 



       
       int *bad_data;

       //record_bad = ff_bad->records[nRec-1];
       record_bad=reclink_bad;
       recnum2_bad=record_bad->recnum;

       segin    = drms_segment_lookup(record_bad, segmentname_badpix);
       arr_bad= drms_segment_read(segin, segin->info->type, &status);
       int nelem=arr_bad->axis[0];
       bad_data=arr_bad->data;
         
     	
       t_obs2_bad=drms_getkey_time(record_bad,keytobs,&status);

       sprint_ut(date_string, t_obs2_bad);
       printf("Most recent bad pixel list read. T_OBS: %s\n", date_string);





 


 recnum[0]=recnum_off;
 recnum[1]=recnum_dark;
 recnum[2]=recnum_bad;
 recnum[3]=recnum2_off;
 recnum[4]=recnum2_dark;
 recnum[5]=recnum2_bad;
 







 tobs_link[0]=t_start;
 tobs_link[1]= t_stop;

 //tobs_link[4]= t_obs_offpoint_front;
 //tobs_link[5]= t_obs_offpoint_side; 
 //tobs_link[6]=t_obs_dark_front;
 //tobs_link[7]= t_obs_dark_side;
 //tobs_link[8]= t_obs_bad_front;
 //tobs_link[9]= t_obs_bad_side;
 //tobs_link[10]= t_obs2_offpoint_front;
 //tobs_link[11]= t_obs2_offpoint_side;
 //tobs_link[12]= t_obs2_dark_front;
 //tobs_link[13]= t_obs2_dark_side;
 //tobs_link[14]= t_obs2_bad_front;
 //tobs_link[15]= t_obs2_bad_side;

 
 return 0;
     

      }

int write_flatfields(DRMS_Array_t *arr_flat, DRMS_Array_t *arrout_new, int camera,
                     long long recnum[6], TIME tobs_link[2], TIME t_0, int focus,
                     struct rotpar rot_new, 
                     struct rotpar rot_cur)
{    
    
#include "module_flatfield_const.h" // !!changed to test  
  TIME t_start=tobs_link[0];

  //TIME t_obs_offpoint_front=tobs_link[4];
  // TIME t_obs_offpoint_side=tobs_link[5];
  // TIME  t_obs_dark_front=tobs_link[6];
  // TIME  t_obs_dark_side=tobs_link[7];
  //TIME t_obs_bad_front=tobs_link[8];
  //TIME t_obs_bad_side=tobs_link[9];
  //TIME t_obs2_offpoint_front=tobs_link[10];
  //TIME t_obs2_offpoint_side=tobs_link[11];
  //TIME  t_obs2_dark_front=tobs_link[12];
  //TIME  t_obs2_dark_side=tobs_link[13];
  //TIME t_obs2_bad_front=tobs_link[14];
  //TIME t_obs2_bad_side=tobs_link[15];

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

   //const char *key1 = "T_STOP";
   //const char *key2 = "HMI_SEQ_ID_FOCUS";
   //const char *key3 = "FLATFIELD_VERSION";

   //const char *key4 = "ROTF_FLATFIELD";
   //const char *key5 = "ROTF_N_PAIRS";
   //const char *key6 = "ROTF_CADENCE";

 

   int status, stat;

 
   const int FLAG_PRELIMINARY=0;
   const int FLAG_FINAL=1;


 /////////////////////////////////////////////////////////////////////////////////
 ///write flatfield with new T_STOP /////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////

   // dataout=drms_clone_records(ff_last, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &stat);
   
 
   dataout = drms_create_records(drms_env,1,filename_flatfield_out,DRMS_PERMANENT,&stat);
   
  if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield_out);
	      exit(EXIT_FAILURE);
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Overwrite current flatfield with definite T_STOP\n");
	      printf("Writing a record on the DRMS for the series %s\n",filename_flatfield_out);

	      //front camera
	     
	      recout = dataout->records[0];

	      status = drms_setkey_time(recout, keytstart, t_start);
	      status=drms_setkey_int(recout, keycamera, cam_id);

	      status=drms_setkey_int(recout, keyfocusflat, focus);
	      status=drms_setkey_time(recout, keytstop, t_0); 
	      status=drms_setkey_int(recout, keyversion, rot_cur.flatfield_version);

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
	
	       	 
	          }
	          {
		    //	   const char *linkname = "bias_dark";
	       	  
		   status=drms_setlink_static(recout, linkdark, recnum_dark);
	
	      	 	       	 
	          }
	          {
		    //   const char *linkname = "perm_bad_pixel";
	        
	      	   status=drms_setlink_static(recout, linkbad, recnum_bad);
	
	      	  
	       }

	         segout = drms_segment_lookup(recout, segmentname);
	   
	         status=drms_segment_write(segout, arr_flat, 0);        //write the file containing the data
	      

	    }


	  drms_close_record(recout, DRMS_INSERT_RECORD); //insert the record in DRMS


    
    
    

    

	  ///////////////////////////////////////////////////////////////////////
	  //////////write new updated flatfield ////////////////////////////////
	  /////////////////////////////////////////////////////////////////////
	 
	 
    
 dataout_new = drms_create_records(drms_env,1,filename_flatfield_out,DRMS_PERMANENT,&stat);


 if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield_out);
	      exit(EXIT_FAILURE);
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Write new flatfield with definite T_STOP\n");
	     
	      recout = dataout_new->records[0];

	      status = drms_setkey_time(recout, keytstart, t_0);
	      status=drms_setkey_int(recout, keycamera, cam_id);

	      status=drms_setkey_int(recout, keyfocusflat, focus);
	      status=drms_setkey_time(recout, keytstop, t_0+365.0*24.0*60.0*60.0); //update T_STOP=T_START+365 days
	      status=drms_setkey_int(recout, keyversion, FLAG_PRELIMINARY);    //update FLAT_FINAL
	      drms_setkey_string(recout, keyinstrument, camera_str);

	      status=drms_setkey_int(recout, keynewpix, rot_new.rotbad);
	      status=drms_setkey_int(recout, keynpairs, rot_new.rotpairs);
	      status=drms_setkey_float(recout, keycadence, rot_new.rotcadence);

	      drms_keyword_setdate(recout);
              /******************************************/
	      /*set link                                */	 
              /******************************************/
	      {
		//  const char *linkname = "offpoint_flatfield";
	    	    
		   status=drms_setlink_static(recout, linkoff, recnum2_off);
	    	       	 printf("link status %d\n", status);
	      }
	      {
	      	//   const char *linkname = "bias_dark";

		   status=drms_setlink_static(recout, linkdark, recnum2_dark);
	      	   printf("link status %d\n", status);
	      }

	      {
	       	//   const char *linkname = "perm_bad_pixel";
	       
		   status=drms_setlink_static(recout, linkbad, recnum2_bad);
		   printf("link status %d\n", status);
	      }


	      segout = drms_segment_lookup(recout, segmentname);
	      drms_segment_write(segout, arrout_new, 0);        //write the file containing the data
	      printf("done\n");
	    }
    



	  drms_close_records(dataout_new, DRMS_INSERT_RECORD); //insert the record in DRMS
	  //	  drms_free_records(dataout_new);

	  return 0;
}



    
int write_flatfield_fid(DRMS_Array_t *arrout_new, int camera,
			TIME tobs_link[2], TIME t_0, int focus,int fid, int fsns,
			struct rotpar rot_new)
{    
    
#include "module_flatfield_const.h"  
  TIME t_start=tobs_link[0];
 


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



	  ///////////////////////////////////////////////////////////////////////
	  //////////write new updated flatfield ////////////////////////////////
	  /////////////////////////////////////////////////////////////////////
	 
	 
    
 dataout_new = drms_create_records(drms_env,1,filename_flatfield_fid,DRMS_PERMANENT,&stat);


 if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",filename_flatfield_fid);
	      exit(EXIT_FAILURE);
	    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Write new flatfield with fid %d\n", fid);
	     
	      recout = dataout_new->records[0];

	      status = drms_setkey_time(recout, keytstart, t_0);
	      status=drms_setkey_int(recout, keycamera, cam_id);
	      status=drms_setkey_int(recout, fidkey, fid);
	      status=drms_setkey_int(recout, fsnskey, fsns);

	      status=drms_setkey_int(recout, keyfocusflat, focus);
	      status=drms_setkey_time(recout, keytstop, t_0+365.0*24.0*60.0*60.0); //update T_STOP=T_START+365 days
	      status=drms_setkey_int(recout, keyversion, FLAG_PRELIMINARY);    //update FLAT_FINAL
	      drms_setkey_string(recout, keyinstrument, camera_str);

	      status=drms_setkey_int(recout, keynewpix, rot_new.rotbad);
	      status=drms_setkey_int(recout, keynpairs, rot_new.rotpairs);
	      status=drms_setkey_float(recout, keycadence, rot_new.rotcadence);

	      drms_keyword_setdate(recout);
            


	      segout = drms_segment_lookup(recout, segmentname);
	      drms_segment_write(segout, arrout_new, 0);        //write the file containing the data
	      printf("done\n");
	    }
    



	  drms_close_records(dataout_new, DRMS_INSERT_RECORD); //insert the record in DRMS


	  return 0;
}



void apod_circ(float rad, float nb, float offx, float offy, float *vd)               
{
  float *rarr;
  rarr=(float *)(malloc(nx*ny*sizeof(float)));
  int i, j;

  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) rarr[j*nx+i]=sqrt(((float)i-((float)nx/2+offx))*((float)i-((float)ny/2+offx))+((float)j-((float)nx/2+offy))*((float)j-((float)ny/2+offy)));
	 
 
  for (j=0; j<ny; ++j){
    for (i=0; i<nx; ++i){
      if (rarr[j*nx+i] < rad) 
	vd[j*nx+i]=1.0;

      if (rarr[j*nx+i] >= rad && rarr[j*nx+i] < (rad+nb))
	vd[j*nx+i]=0.5*cos(M_PI/nb*(rarr[j*nx+i]-rad))+0.5;

      if (rarr[j*nx+i] >= (rad+nb))
	vd[j*nx+i]=0.0;

      

    }
  }
	  
}

