
// write out flatfields


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>



char *module_name    = "write_flatfield";    //name of the module

                                      //arguments of the module
ModuleArgs_t module_args[] =        
{
  {ARG_STRING, "instrument", ""}, //instrument AIA or HMI
  {ARG_STRING, "series_offpoint", ""}, //series name of offpoint flat field
  {ARG_STRING, "series_badpix", ""},   //series name of bad pixel series
  {ARG_STRING, "series_dark", ""},     //series name of dark series
  {ARG_STRING, "series_flatfield", ""}, //series name of flatfield series
  {ARG_STRING,"file_flatfield"},    //file name for flatfield (binary file of floats) (not required for AIA)
  {ARG_INT, "camera", "0", "1-2"},      //camera for HMI (1(side) or 2(front))
  {ARG_STRING, "wave_str", "dd"},         //wave string for AIA
  {ARG_INT, "focus",  "0"},             //focus position (1-16) for HMI
  {ARG_TIME, "t_start"},                //t_start
  {ARG_TIME, "t_stop", "DRMS_MISSING_VALUE"},          //t_stop, if omitted t_stop=t_start+1 year
  {ARG_TIME, "t_obs_dark"},             //for linking T_OBS of corresponding dark
  {ARG_TIME, "t_obs_badpix"},           
  {ARG_TIME, "t_obs_offpoint"},         //for linking T_OBS of corresponding bad pixel series
  {ARG_INT,  "flatfield_version", "0"}, //flatfield version (0: preliminary, >0: final)
  {ARG_INT, "nx", "4096"},              //x dim of image
  {ARG_INT, "ny", "4096"},               //y dim of image
  {ARG_END}
};



/////////////////////////////////////////////////////////////////////////////////


/*-------------------------------------------------------------*/
/*                                                             */
/*   DoIt is the entry point of the module                     */
/*   the name MUST be DoIt for a DRMS module                   */
/*                                                             */
/*-------------------------------------------------------------*/

int DoIt(void)
{

 
  int i,j,k;   
  int status=0;

  const char *instrument;
  instrument=cmdparams_get_str(&cmdparams, "instrument", &status);

  if (strcmp(instrument,"HMI") != 0 && strcmp(instrument,"AIA") != 0){printf("Instrument invalid\n"); exit(EXIT_FAILURE);}

  int inst_hmi=!strcmp(instrument,"HMI");
  int inst_aia=!strcmp(instrument,"AIA");



  const char *input_flatfield;
  input_flatfield = cmdparams_get_str(&cmdparams, "file_flatfield", &status);
  

  int foc_id; 
  if (inst_hmi){
    foc_id=cmdparams_get_int(&cmdparams, "focus", &status);
   }

 
  const char *wavelength;
  int vvd=0, vvv=-1;
  short aia_cam[23]={1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4};

  if (inst_aia)
    {
      wavelength=cmdparams_get_str(&cmdparams, "wave_str", &status);
char **wavestrs=(char **)(malloc(23*sizeof(char *)));
 
  
  wavestrs[0]="131_THIN";
  wavestrs[1]="131_THICK";
  wavestrs[2]="131_OPEN";
  wavestrs[3]="335_THIN";
  {char *wv="335_THICK"; wavestrs[4]=wv;}
  {char *wv="335_OPEN"; wavestrs[5]=wv;}
  {char *wv="193_THIN"; wavestrs[6]=wv;}
  {char *wv="193_THICK"; wavestrs[7]=wv;}
  {char *wv="193_OPEN"; wavestrs[8]=wv;}
  {char *wv="211_THIN"; wavestrs[9]=wv;}
  {char *wv="211_THICK"; wavestrs[10]=wv;}
{char *wv="211_OPEN"; wavestrs[11]=wv;}
{char *wv="171_THIN"; wavestrs[12]=wv;}
{char *wv="171_THICK"; wavestrs[13]=wv;}
{char *wv="1600"; wavestrs[14]=wv;}
{char *wv="1700"; wavestrs[15]=wv;}
{char *wv="4500"; wavestrs[16]=wv;}
{char *wv="94_THIN"; wavestrs[17]=wv;}
{char *wv="94_THICK"; wavestrs[18]=wv;}
{char *wv="94_OPEN"; wavestrs[19]=wv;}
{char *wv="304_THIN"; wavestrs[20]=wv;}
{char *wv="304_THICK"; wavestrs[21]=wv;}
{char *wv="304_OPEN"; wavestrs[22]=wv;}

 



for (i=0; i<23; ++i){vvd += !strcmp(wavestrs[i], wavelength); if (!strcmp(wavestrs[i], wavelength)) vvv=i;}
if (vvd == 0){printf("nonexisting wavelength id\n"); exit(EXIT_FAILURE);} else {printf("number %d\n", vvv);}
    }





  int cam_id;
  if (inst_hmi) cam_id=cmdparams_get_int(&cmdparams, "camera", &status);

 char *camera_string;

 if (inst_hmi)
  {
    if (cam_id == 1) camera_string="HMI_SIDE1";
      if (cam_id == 2) camera_string="HMI_FRONT2";
      if (cam_id < 1 || cam_id > 2){printf("wrong camera id for HMI\n"); exit(EXIT_FAILURE);}
  }

 if (inst_aia){
    if (aia_cam[vvv] == 1) camera_string="AIA_ATA1";
    if (aia_cam[vvv] == 2) camera_string="AIA_ATA2";
    if (aia_cam[vvv] == 3) camera_string="AIA_ATA3";
    if (aia_cam[vvv] == 4) camera_string="AIA_ATA4";
   }


TIME t_start=cmdparams_get_time(&cmdparams, "t_start", NULL);
TIME t_stop_preliminary=cmdparams_get_time(&cmdparams, "t_stop", &status);

     if (status != 0)  {printf("in here\n"); t_stop_preliminary=t_start+365.0*24.0*60.0*60.0;}

     TIME t_dark=cmdparams_get_time(&cmdparams, "t_obs_dark", NULL);
     TIME t_bad=cmdparams_get_time(&cmdparams, "t_obs_badpix", NULL);
     TIME t_offpoint=cmdparams_get_time(&cmdparams, "t_obs_offpoint", NULL);

       
      const char *series_name_offpoint;
       series_name_offpoint=cmdparams_get_str(&cmdparams, "series_offpoint", NULL);

       drms_series_exists(drms_env, series_name_offpoint, &status);
       if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",series_name_offpoint);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
     
       const char *series_name_dark;
       series_name_dark=cmdparams_get_str(&cmdparams, "series_dark", NULL);
       if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",series_name_dark);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
     


       const char *series_name_bad;
       series_name_bad=cmdparams_get_str(&cmdparams, "series_badpix", NULL);
       if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",series_name_bad);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 

    
       const char *series_name_flatfield;
       series_name_flatfield=cmdparams_get_str(&cmdparams,"series_flatfield", NULL);
       if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",series_name_flatfield);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 

       int version;
      version=cmdparams_get_int(&cmdparams, "flatfield_version", NULL);
    

   int stat = DRMS_SUCCESS; 
  int  error       = 0;


DRMS_RecordSet_t *dataout, *dataout_off, *dataout_dark, *dataout_bad;
DRMS_RecordSet_t *ff_offpoint, *ff_bad, *ff_dark, *ff_flat;
DRMS_Record_t *recout  = NULL;
DRMS_Record_t *record_off, *record_dark, *record_bad, *record_flat;
DRMS_Segment_t *segout = NULL, *segin = NULL;
DRMS_Array_t *arrout, *arr_bad;    
long long recnum_off, recnum_dark, recnum_bad;
DRMS_Type_t type;
DRMS_Type_t type_time      = DRMS_TYPE_TIME;
size_t bytes_read;
FILE *fgram;
float *gout1;
int *gout_bad;
int nbad=0;
int *badpix;

       const int nx=cmdparams_get_int(&cmdparams, "nx", NULL);
       const int ny=cmdparams_get_int(&cmdparams, "ny", NULL);

       int axisin[2];                                                   //size of input arrays
       int axisout[2]={nx,ny};

     
 //KEYWORDS  ***********************************************************************/

      const char *primekey1 = "CAMERA";
      const char *primekey2 = "T_OBS";                            //1st prime key of output data
      const char *primekey2a="T_START";      
      const char *primekey3  = "HMI_SEQ_ID_FOCUS";                              //2nd prime key of output data
      

      const char *key1 = "FSN_FIRST";     // output keyword //copy
      const char *key2 = "FSN_LAST";  // output keyword //copy 
      const char *key3 = "FSN_INPUT";  // output keyword //copy 

      const char *key4= "T_STOP";
      const char *key5= "FLATFIELD_VERSION";
      const char *key6="INSTRUME";

      const char *key1aia="WAVE_STR";
  
      char *camstr_side="[1]";
      char *camstr_front="[2]";

      char *camstr;

      char t0_off[256];
      sprint_ut(t0_off, t_offpoint);

      char t0_flat[256];
      sprint_ut(t0_flat, t_start);

      char query[256]={""}, query_flat[256]={""}, query_dark[256]={""}, query_bad[256]={""};
      char q1_hmi[256]={""};
      char q1_aia[256]={""};
      char t_orig_str[256];
      sprint_ut(t_orig_str, 0.0);


      strcat(query,series_name_offpoint);
      strcat(query_flat, series_name_flatfield);

      if (inst_hmi)
	{
        if (cam_id == 1){camstr=camstr_side;}
	if (cam_id == 2){camstr=camstr_front;}
	strcat(q1_hmi, camstr);
	strcat(query, q1_hmi);
	strcat(query_flat, q1_hmi);
	}

      if (inst_aia)
	{
	  
	  strcat(q1_aia,"[");
	  strcat(q1_aia, wavelength);
	  strcat(q1_aia, "]");

	  strcat(query, q1_aia);
	}

strcat(query, "[");
strcat(query,t0_off);
strcat(query,"]");

if (inst_hmi)
  {
    strcat(query_flat, "[");
    strcat(query_flat,t_orig_str);
    strcat(query_flat,"-");
    strcat(query_flat,t0_flat);
    strcat(query_flat,"]");
  }

   if (inst_hmi)
     {
       strcat(query,"[");

       char ffnumb[2]={""};
       sprintf(ffnumb, "%2.2d", foc_id);

       strcat(query, ffnumb);
       strcat(query, "]");

     }

  


printf("query offpoint %s\n", query);
ff_offpoint = drms_open_records(drms_env,query,&status);

   if (status == DRMS_SUCCESS && ff_offpoint!= NULL && ff_offpoint->n > 0)
	{
	  record_off=ff_offpoint->records[0];
	  recnum_off=record_off->recnum;
	}
	else
	  {
 
	    printf("Offpoint flatfield not found\n");
	    exit(EXIT_FAILURE);
	  }

//search for dark
  
    char t0_dark[256];
    sprint_ut(t0_dark, t_dark);

if (inst_hmi)
  {
   strcat(query_dark,series_name_dark);
   strcat(query_dark,q1_hmi);
   strcat(query_dark,"[");
   strcat(query_dark,t0_dark);
   strcat(query_dark,"]");
   printf("query string %s\n", query_dark);
  }

if (inst_aia)
  {strcat(query_dark,series_name_dark);
   strcat(query_dark,q1_aia);
   strcat(query_dark,"[");
   strcat(query_dark,t0_dark);
   strcat(query_dark,"]");
   printf("query string %s\n", query_dark);
  }


   ff_dark = drms_open_records(drms_env,query_dark,&status);
      

  

   if (status == DRMS_SUCCESS && ff_dark != NULL && ff_dark->n > 0)
	{
	  record_dark = ff_dark->records[0];
	  recnum_dark=record_dark->recnum;
	}
	else
	  {
 
	    printf("Dark record not found\n");
	    exit(EXIT_FAILURE);
	  }


//search for badpix
   char t0_bad[256];
   sprint_ut(t0_bad, t_bad);

if (inst_hmi)
  {
   strcat(query_bad,series_name_bad);
   strcat(query_bad,q1_hmi);
   strcat(query_bad,"["); 
   strcat(query_bad,t0_bad);
   strcat(query_bad,"]");
   printf("query string %s\n", query_bad);
  }

if (inst_aia)
  {
    strcat(query_bad,series_name_bad);
    strcat(query_bad,q1_aia);
    strcat(query_bad,"[");
    strcat(query_bad,t0_bad);
    strcat(query_bad,"]");
    printf("query string %s\n", query_bad);
  }


   ff_bad = drms_open_records(drms_env,query_bad,&status);
      

  

   if (status == DRMS_SUCCESS && ff_bad != NULL && ff_bad->n > 0)
	{
	  record_bad = ff_bad->records[0];
	  recnum_bad=record_bad->recnum;

	  segin    = drms_segment_lookup(record_bad, "bad_pixel_list");
	  arr_bad = drms_segment_read(segin, segin->info->type, &status);
	  nbad=arr_bad->axis[0];
	  badpix=arr_bad->data;

	}
	else
	  {
 
	    printf("Bad pixel record not found\n");
	    exit(EXIT_FAILURE);
	  }

   type    = DRMS_TYPE_FLOAT;                                     //type of the output data
   float *f_flatfield;
   
     
      fgram=fopen(input_flatfield, "rb");

      if (fgram==NULL){fputs("Input file not found", stderr); printf("Input file not found\n"); exit(EXIT_FAILURE);}
       f_flatfield=(float *)(malloc(nx*ny*sizeof(float)));
       bytes_read=fread(f_flatfield,sizeof(float),nx*ny,fgram);
   
       fclose(fgram);

      
       arrout  = drms_array_create(type,2,axisout,NULL,&status);
       gout1=arrout->data;

      
       for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)  gout1[j*nx+i]=f_flatfield[j*nx+i];

      
     

  /***********************************************************************************************************/
  /*CHECK WHETHER THE FLATFIELD SERIES EXIST                                                                    */
  /***********************************************************************************************************/
  
    
   drms_series_exists(drms_env, (char *)series_name_flatfield, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Flatfield series %s doesn't exist\n",series_name_flatfield);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 


               
      /***********************************************************************************************************/
      /*WRITE FLATFIELD DATA IN THE DRMS                                                                         */
      /***********************************************************************************************************/
   
              
      dataout = drms_create_records(drms_env,1,(char *)series_name_flatfield,DRMS_PERMANENT,&stat);
       
    

	  if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n",series_name_flatfield);
	      exit(EXIT_FAILURE);
    }
	  if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Writing a record on the DRMS for the series %s\n",series_name_flatfield);
	      recout = dataout->records[0];

	      status=0;	      
	      status += drms_setkey_time(recout, primekey2a, t_start);
	      status += drms_setkey_time(recout,key4,t_stop_preliminary);        
	      status += drms_setkey_int(recout,key5,version);

	      if (inst_hmi) status += drms_setkey_int(recout, primekey1,cam_id);
	      if (inst_hmi) status += drms_setkey_int(recout,primekey3,foc_id);

	      	      status += drms_setkey_string(recout, key6, camera_string);
	   
	      if (inst_hmi)
		{
	       status += drms_setkey_int(recout, "ROTF_FLATFIELD", 0);
	       status += drms_setkey_int(recout, "ROTF_N_PAIRS", 0);
	       status += drms_setkey_float(recout, "ROTF_CADENCE", 0.0);
		}

	      if(inst_aia)
		{
		  status += drms_setkey_string(recout, key1aia, wavelength);
		
		}

	      if (status != 0){printf("error setting keywords"); exit(EXIT_FAILURE);}
	       drms_keyword_setdate(recout);
	      
              /******************************************/
	      /*set link                                */	 
              /******************************************/
	       status=0;
	      	 {
	        	   const char *linkname = "OFFPOINT_FLAT";
			   status+=drms_setlink_static(recout, linkname, recnum_off);
			   printf("link status offpoint %d\n", status);
	       	 }
	       	 {
	       	   const char *linkname = "DARK";
		   status+=drms_setlink_static(recout, linkname, recnum_dark);
		   printf("link status dark %d\n", status);

		 }
		 {
	       	   const char *linkname = "BAD_PIXEL";
		   status+=drms_setlink_static(recout, linkname, recnum_bad);
	      	   printf("link status badpix %d\n", status);
		   if (status != 0){printf("error setting links\n"); exit(EXIT_FAILURE);}
	       	 }

		 //**************************************************/
		 //write out flatfield                               /
		 //**************************************************/	
	     
		 printf("write out segment\n");
		 segout = drms_segment_lookup(recout, "flatfield");
		 if (segout == NULL){printf("could not find segment\n"); exit(EXIT_FAILURE);}

		 status=drms_segment_write(segout, arrout, 0);        //write the file containing the data
		 if (status != 0){printf("could not write segment\n"); exit(EXIT_FAILURE);}
		   }

		 printf("done\n");
	    

	  drms_close_records(dataout, DRMS_INSERT_RECORD); //insert the record in DRMS

	
	  drms_free_array(arrout);




	  printf("COMPLETED!\n");

 
    
   
      
     return 0;

  

}
















 






