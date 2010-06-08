
// write out flatfields


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>



char *module_name    = "write_badpix";    //name of the module

                                      //arguments of the module
ModuleArgs_t module_args[] =        
{
  {ARG_STRING, "instrument"},         //HMI or AIA
  {ARG_STRING, "series_badpix", ""},  //bad pixel series name
  {ARG_STRING, "file_badpix"},        //bad pixels file (binary file of ints)
  {ARG_INT, "camera", "0"},           //camera for HMI (1(side) or 2(front))
  {ARG_STRING, "wave_str", "dd"},     //wave string for AIA
  {ARG_TIME, "t_obs"},                //T_OBS
  {ARG_INT,  "nbad", "0", "0-16777216"},   //number of bad pixels in file
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
 

  const char *input_badpix;
  input_badpix = cmdparams_get_str(&cmdparams, "file_badpix", NULL);

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

      int cam_id;
      cam_id=cmdparams_get_int(&cmdparams, "camera", NULL);

      const char *camera_string;
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


     TIME t_start=cmdparams_get_time(&cmdparams, "t_obs", NULL);
     TIME t_offpoint=t_start;

     int nbad=cmdparams_get_int(&cmdparams, "nbad", NULL);
 
       const char *series_name_bad;
       series_name_bad=cmdparams_get_str(&cmdparams, "series_badpix", NULL);


     
      
       int stat = DRMS_SUCCESS;


       DRMS_RecordSet_t *dataout, *dataout_bad;
   
       DRMS_Record_t *recout  = NULL;
       DRMS_Segment_t *segout = NULL;
       DRMS_Array_t *arrout, *arrout_bad;      
       long long recnum_bad;
       DRMS_Type_t type;
       DRMS_Type_t type_time      = DRMS_TYPE_TIME;
       size_t bytes_read;
       FILE *fgram;
       float *gout1;
       int *gout_bad;

    
         
 //KEYWORDS  ***********************************************************************/

      const char *primekey1 = "CAMERA";
      const char *primekey2 = "T_OBS";                            //1st prime key of output data

      const char *key6="INSTRUME";

      const char *key1aia="WAVE_STR";
 

  printf("START!\n");




     
       //********************************************
      //write bad pixel list
      //********************************************
      drms_series_exists(drms_env, series_name_bad, &status);
	  if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Series %s doesn't exist\n",series_name_bad);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Series %s exists.\n",series_name_bad);
	}

      int *f_bad;
      fgram=fopen(input_badpix, "rb");
      if (fgram==NULL){fputs("File error", stderr); printf("could not read file\n"); exit(EXIT_FAILURE);}
      f_bad=(int *)(malloc(nbad*sizeof(int)));
      bytes_read=fread(f_bad,sizeof(int),nbad,fgram);
      fclose(fgram);

      int axisout_bad[1];
      axisout_bad[0]=nbad;                      
      type    = DRMS_TYPE_INT;            
                          
      arrout_bad  = drms_array_create(type,1,axisout_bad,NULL,&status);
      if (status !=0 || arrout_bad == NULL){printf("could not create array\n"); exit(EXIT_FAILURE);}

      gout_bad  = arrout_bad->data;

      for (j=0; j<nbad; ++j) gout_bad[j]=f_bad[j];
      dataout_bad = drms_create_records(drms_env,1,(char *)series_name_bad,DRMS_PERMANENT,&stat);
    
      if (stat != DRMS_SUCCESS)
	    {
	      printf("Could not create a record for the series %s\n", series_name_bad);
	      exit(EXIT_FAILURE);
	    }
      if (stat == DRMS_SUCCESS)
	    {	  
	      printf("Writing a record on the DRMS for the series %s\n",series_name_bad);
	      recout = dataout_bad->records[0];
	      recnum_bad=recout->recnum;

	      status=0;
	      status += drms_setkey_time(recout, primekey2, t_start);

	      if (inst_hmi) status += drms_setkey_int(recout, primekey1,cam_id);

	      if (inst_aia) status += drms_setkey_string(recout, key1aia, wavelength);
	      
	      status += drms_setkey_string(recout, key6, camera_string);


	      if (status != 0){printf("error setting keywords"); exit(EXIT_FAILURE);}
	      drms_keyword_setdate(recout);

	      segout = drms_segment_lookup(recout, "bad_pixel_list");
	      if (segout == NULL){printf("could not find segment\n"); exit(EXIT_FAILURE);}

	      status=drms_segment_write(segout, arrout_bad, 0);
	      if (status != 0){printf("could not write segment\n"); exit(EXIT_FAILURE);}

	    }

      drms_close_records(dataout_bad, DRMS_INSERT_RECORD);;
      drms_free_array(arrout_bad);
      printf("done\n");

      printf("COMPLETED!\n");

 
    
   
      
     return 0;

  

}

///////////////////////////////




