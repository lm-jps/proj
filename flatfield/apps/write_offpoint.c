
// write out flatfields


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>



char *module_name    = "write_offpoint";    //name of the module

                                      //arguments of the module
ModuleArgs_t module_args[] =        
{
  {ARG_STRING, "instrument"},      //HMI or AIA
  {ARG_STRING, "series_offpoint"},   //series name of offpoint flat field
  {ARG_STRING,"file_offpoint"},          //file name of offpoint flatfield (binary file of floats)
  {ARG_INT, "camera", "0"},      //camera for HMI (1 or 2)
  {ARG_STRING, "wave_str", "dd"},  //wave string for AIA
  {ARG_INT, "focus",  "0"},              //focus position (1-16) for HMI
  {ARG_TIME, "t_obs"},                   //T_OBS
  {ARG_INTS, "fsn_list_offpoint", "-1,-1"},       //comma separated list of FSNs
  {ARG_INTS, "fsn_list_pzt", "-1,-1"}, 
  {ARG_INT, "nx", "4096"},                  //x dim of image 
  {ARG_INT, "ny", "4096"},                   //y dim of image
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

  int status=DRMS_SUCCESS;
  int i,j,k;   
  const char *input_offpoint;
  input_offpoint = cmdparams_get_str(&cmdparams, "file_offpoint", &status); 

  const char *instrument;
  instrument=cmdparams_get_str(&cmdparams, "instrument", &status);

  if (strcmp(instrument,"HMI") != 0 && strcmp(instrument,"AIA") != 0){printf("Instrument invalid\n"); exit(EXIT_FAILURE);}

  int inst_hmi=!strcmp(instrument,"HMI");
  int inst_aia=!strcmp(instrument,"AIA");

  
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

  int foc_id;
  if (inst_hmi)foc_id=cmdparams_get_int(&cmdparams, "focus", &status);

  int cam_id;
  cam_id=cmdparams_get_int(&cmdparams, "camera", &status);

  char *camera_string;

if (inst_hmi)
  {
    if (cam_id == 1) camera_string="HMI_SIDE1";
      if (cam_id == 2) camera_string="HMI_FRONT2";
      if (cam_id < 1 || cam_id > 2){printf("wrong camera id for HMI\n"); exit(EXIT_FAILURE);}
  }

 if (inst_aia)
   {
    if (aia_cam[vvv] == 1) camera_string="AIA_ATA1";
    if (aia_cam[vvv] == 2) camera_string="AIA_ATA2";
    if (aia_cam[vvv] == 3) camera_string="AIA_ATA3";
    if (aia_cam[vvv] == 4) camera_string="AIA_ATA4";
     }


     
     TIME t_offpoint=cmdparams_get_time(&cmdparams, "t_obs", &status);
    

     int *fsn_list_offset;
     fsn_list_offset=(int *)(malloc(sizeof(int)*1024));
     int noffset=cmdparams_get_intarr(&cmdparams, "fsn_list_offpoint", &fsn_list_offset, &status);

     int *fsn_list_pzt;
     fsn_list_pzt=(int *)(malloc(sizeof(int)*1024));
     int npzt=cmdparams_get_intarr(&cmdparams, "fsn_list_pzt", &fsn_list_pzt, &status);

   
     int fsn_first=fsn_list_offset[0];
     int fsn_last=fsn_list_offset[noffset-1];
   
    
     char  fsn_string[13*1024]={""};
     if (noffset > 1024){ printf("too many offpoint frames\n"); exit(EXIT_FAILURE);}

      
     for (k=0; k<noffset; ++k)
	{ 
	  char ffnumb[12]={""};
	  sprintf(ffnumb, "%d", fsn_list_offset[k]);
	  strcat(fsn_string, ffnumb);
	  if (k<(noffset-1)) strcat(fsn_string, ",");
 	}

     if (fsn_list_pzt[0] != -1)
       {
	 strcat(fsn_string, ";PZT_FSN:");
	 for (k=0; k<npzt; ++k)
	   {
	  char ffnumb[12]={""};
	  sprintf(ffnumb, "%d", fsn_list_pzt[k]);
	  strcat(fsn_string, ffnumb);
	  if (k<(npzt-1)) strcat(fsn_string, ",");
	   }
       }
 


     const char *series_name_offpoint;
     series_name_offpoint=cmdparams_get_str(&cmdparams, "series_offpoint", &status);

    

     int stat = DRMS_SUCCESS;
     int  error       = 0;


     DRMS_RecordSet_t *dataout, *dataout_off, *dataout_dark, *dataout_bad;
   
     DRMS_Record_t *recout  = NULL;
     DRMS_Segment_t *segout = NULL;
     DRMS_Array_t *arrout, *arrout_off, *arrout_dark, *arrout_bad;     ;
     long long recnum_off, recnum_dark, recnum_bad;
     DRMS_Type_t type = DRMS_TYPE_FLOAT;
     DRMS_Type_t type_time      = DRMS_TYPE_TIME;
     size_t bytes_read;
     FILE *fgram;
     float *gout1;
     int *gout_bad;

     int nx=cmdparams_get_int(&cmdparams, "nx", &status);
     int ny=cmdparams_get_int(&cmdparams, "ny", &status);

     int axisin[2];                                                   //size of input arrays
     int axisout[2]={nx,ny};

     
 //KEYWORDS  ***********************************************************************/

      const char *primekey1 = "CAMERA";
      const char *primekey2 = "T_OBS";                            //1st prime key of output data
      const char *primekey3  = "HMI_SEQ_ID_FOCUS";                              //2nd prime key of output dat
      

      const char *key1 = "FSN_FIRST";     // output keyword //copy
      const char *key2 = "FSN_LAST";  // output keyword //copy 
      const char *key3 = "FSN_INPUT";  // output keyword //copy 

      const char *key6="INSTRUME";

      const char *key1aia="WAVE_STR";
    

  printf("START!\n");



  /***********************************************************************************************************/
  /*CHECK WHETHER THE OFFPOINT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
  
    
      drms_series_exists(drms_env, series_name_offpoint, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",series_name_offpoint);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Output series %s exists.\n",series_name_offpoint);
	}


     
      //READ AND WRITE OFFPOINT FLATFIELD SERIES

      printf("Reading OFFPOINT flatfield\n"); 

      float *f_offpoint;

       fgram=fopen(input_offpoint, "rb");

       if (fgram==NULL){fputs("File error", stderr); printf("Could not read file\n"); exit(EXIT_FAILURE);}
       f_offpoint=(float *)(malloc(nx*ny*sizeof(float)));
       bytes_read=fread(f_offpoint,sizeof(float),nx*ny,fgram);
   
       fclose(fgram);
       printf("done\n");
  
  
      
       arrout_off  = drms_array_create(type,2,axisout,NULL,&status);
       if (status !=0 || arrout_off == NULL){printf("could not create array\n"); exit(EXIT_FAILURE);}

       gout1=arrout_off->data;

      
       for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)  gout1[j*nx+i]=f_offpoint[j*nx+i];
       dataout_off = drms_create_records(drms_env,1,(char *)series_name_offpoint,DRMS_PERMANENT,&stat);
      

       if (stat != DRMS_SUCCESS)
      	    {
      	      printf("Could not create a record for the series %s\n",series_name_offpoint);
      	      exit(EXIT_FAILURE);
      	    }

      	  if (stat == DRMS_SUCCESS)
      	    {	  
      	      printf("Writing a record on the DRMS for the series %s\n",series_name_offpoint);
	      recout= dataout_off->records[0];
	      recnum_off=recout->recnum;

	      //write keywords
	      status=0;
	      status += drms_setkey_time(recout, primekey2, t_offpoint);
	      if (inst_hmi) status += drms_setkey_int(recout, primekey1, cam_id);
	      if (inst_hmi) status += drms_setkey_int(recout, primekey3, foc_id);

	      if (inst_aia) status += drms_setkey_string(recout, key1aia, wavelength);
	     
	      status += drms_setkey_int(recout, key1, fsn_first);
	      status += drms_setkey_int(recout, key2, fsn_last);
	      status += drms_setkey_string(recout, key3, fsn_string);
	      status += drms_setkey_string(recout, key6, camera_string);

	      if (status != 0){printf("error setting keywords"); exit(EXIT_FAILURE);}
	      drms_keyword_setdate(recout);


	      segout = drms_segment_lookup(recout, "offpoint_flatfield");
	      if (segout == NULL){printf("could not find segment\n"); exit(EXIT_FAILURE);}

	      status=drms_segment_write(segout, arrout_off, 0);
	      if (status != 0){printf("could not write segment\n"); exit(EXIT_FAILURE);}

	    }

	  drms_close_records(dataout_off, DRMS_INSERT_RECORD);
	  printf("done\n");





	  //drms_free_array(arrout);
	  drms_free_array(arrout_off);




	  printf("COMPLETED!\n");

 
    
   
      
     return 0;

  

}


