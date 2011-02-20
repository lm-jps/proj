/*
 * module_flatfield - Calculates rotational flatfields (first step) and produces cosmic ray records (first step)
 *
 */

/**
\defgroup module_flatfield module_flatfield - 

\par Synopsis
\code
module_flatfield input_series= camera= datum= camera= cadence= cosmic_rays= flatfield= fsn_first= fsn_last= cosmic_ray_series
\endcode

\details

module_flatfield calculates the rotational flatfield by FID, and produces the cosmic ray records. Both are written to an intermediate series. The flatfield update is produced by module_flatfield_combine, and the cosmic ray series by cosmic_ray_post. If cosmic_ray_series is set to hmi.cosmic_rays, module_flatfield produces the final cosmic ray records directly. 



\par Options

\par Mandatory arguments:

\li \c input_series="string" where string is the series name of the input level 1 data (hmi.lev1 or hmi.lev1_nrt)
\li \c camera=cam,  side camera: cam=1, front camera: cam=2
\li \c fid=fid: observable FID (10054-10159)
\li \c datum="date" date="yyyy.mm.dd" TAI-day for which the flatfields and cosmic_rays are calculated. End of TAI-day datum is T_START of updated flatfield
\li \c cadence=cadence: integer number in seconds for the cadence for the framelist that is run (Example: Framelist Mod C, front camera: cadence=45, side camera: cadence=135)
\li \c cosmic_rays=flag: flag=1 if cosmic ray records should be produced, otherwise flag=0
\li \c flatfield=flag: flag=1 if rotational flatfield should be produced, otherwise flag=0

\par Optional arguments:
\li \c cosmic_ray_series="string" where string is the series name of the output cosmic ray series (default is "su_production.cosmic_rays")
\li \c fsn_first: first FSN for which cosmic ray record is desired, which is included in the flatfield calculation, overrides argument datum: datum still needed for flatfield identification
\li \c fsn_last: last FSN or which cosmic ray record is desired, which is included in the flatfield calculation, overrides argument datum


\par Examples

\b Example 1:

\code
module_flatfield input_series="hmi.lev1" camera=2 datum="2010.10.09" fid=10059 camera=2 cosmic_rays=1 flatfield=1 cadence=45 
module_flatfield input_series="hmi.lev1" camera=2 datum="2010.10.09" fid=10059 camera=2 cosmic_rays=1 flatfield=1 cadence=45 fsn_first=12050442 fsn_last=12050459 cosmic_ray_series="hmi.cosmic_rays"
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
#include <module_flatfield.h>

char *module_name  = "module_flatfield";    //name of the module

#define kRecSetIn      "input_series" //name of the series containing the input filtergrams
#define kRecSetOut      "cosmic_ray_series" //name of the series containing the output series for cosmic rays // optional
#define kDSCosmic      "cosmic_rays" //name of the cosmic ray series
#define kDSFlatfield      "flatfield" 
#define cadence_name "cadence"  //cadence in sec
#define fid_name "fid" //name of fid argument
#define cam_name "camera" // 0 side, 1 front
#define fsnf_name "fsn_first"
#define fsnl_name "fsn_last"
#define datumn  "datum"


#define minval(x,y) (((x) < (y)) ? (x) : (y))                                        
#define maxval(x,y) (((x) < (y)) ? (y) : (x))

ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn, "",  "Input data series."},
     {ARG_STRING, datumn, "", "datum string"},
     {ARG_INT, kDSCosmic, "", "Cosmic rays flag"},
     {ARG_INT, kDSFlatfield, "", "Flatfield flag"},
     {ARG_INT, cadence_name, "", "Cadence in sec"},
     {ARG_INT, fid_name, "", "FID"},
     {ARG_INT, cam_name, "", "Camera"},
     {ARG_INT, fsnf_name, "0"},
     {ARG_INT, fsnl_name, "2147483647"},
     {ARG_STRING, kRecSetOut, "default",  "Cosmic ray output series."},
     
     {ARG_END}
};



/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






long sign(double);
double mat_rh(long[], double[], double[], int);
void nine_diag(long[], double[], int, int, double[], long[], double);
void blockiter(double[], double[], double[], double*, double[], int, int, double);
void tridag(double veca[], double vecb[], double vecc[], double vecr[], double vecu[], int n);
void mat_mult(double[], double[], double[], double[], double[], int);
void printtime();
void highpass(int, int, double, double[nx][ny]);
void derotation(double, double, double, double, double, double, double [2][nx][ny]);

void limb_darkening(double radius, double cent_x, double cent_y, double *b, int order, double *limb_dark);
int flatfield(double *rhsp, double *rhsm,short *badpix, int pairs, double *gout, double *param, struct code_param, double deltat);
void apod_circ(float rad, float nb, float offx, float offy, float *vd);


int get_flat(char *query, int camera, float *flatfield, float *offpoint, 
	     short *bad);

int get_latest_bad(TIME t_0, int camera, const char *fname, const char *segname, short *bad, long long *recnum);

int get_latest_flat(TIME t_0, int camera, const char *fname, const char *segname, float *offpoint, long long *recnum, int *focus);



int write_flatfield_fid(DRMS_Array_t *arrout_new, int camera, long long recnum_offpoint, 
		        TIME t_0, int focus,int fid,
			struct rotpar rot_new);
/*-------------------------------------------------------------*/
/*                                                             */
/*   DoIt is the entry point of the module                     */
/*   the name MUST be DoIt for a DRMS module                   */
/*                                                             */
/*-------------------------------------------------------------*/

int DoIt(void)
{

#include "module_flatfield_const.h"

  ////////////////////////////////////////
  //read out input parameters
  /////////////////////////////////////////


  const char *inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL); //cmdparams is defined in jsoc_main.h
  const char *cosmic_ray_series = cmdparams_get_str(&cmdparams, kRecSetOut, NULL); //cmdparams is defined in jsoc_main.h
  if (strcmp(cosmic_ray_series, "default") == 0){cosmic_ray_series=filename_cosmic; }

  int cosmic_flag=cmdparams_get_int(&cmdparams, kDSCosmic, NULL);
  int flatfield_flag=cmdparams_get_int(&cmdparams, kDSFlatfield, NULL);
  int cadence=cmdparams_get_int(&cmdparams, cadence_name, NULL);

  int fid=cmdparams_get_int(&cmdparams, fid_name, NULL);
  int cam=cmdparams_get_int(&cmdparams, cam_name, NULL)-1;
  int fsn_first=cmdparams_get_int(&cmdparams, fsnf_name, NULL);
  int fsn_last=cmdparams_get_int(&cmdparams, fsnl_name, NULL);

  const char *datum=cmdparams_get_str(&cmdparams, datumn, NULL);
 

  if (fid < minfid || fid > maxfid){printf("Not an observable FID\n"); return 1;}

  int  status, status1, status2; 
  int status0=DRMS_SUCCESS;
  int stat = DRMS_SUCCESS; 
  int status_flatfield=DRMS_SUCCESS;
  int status_flatfield_relative=DRMS_SUCCESS;
 

  DRMS_RecordSet_t *data, *data_flat, *data_last, *dataout;
  DRMS_Record_t  *ff_last=NULL;
  DRMS_Array_t *arrin0;
  double *arrinL0; //pointer to data
  DRMS_Segment_t *segin = NULL, *segin_flat=NULL, *segout = NULL, *segin_offpoint=NULL, *segout_val=NULL, *segout_sig=NULL;
  DRMS_Record_t *rec_flat, *recout, *reclink_off;

  DRMS_Type_t type_time = DRMS_TYPE_TIME;
  DRMS_Type_t type_int = DRMS_TYPE_INT;
  DRMS_Type_t type_float = DRMS_TYPE_FLOAT;
  DRMS_Type_t type_double = DRMS_TYPE_DOUBLE;

  // int axisin[2]={nx,ny};                   
                                 //size of input arrays

  double deltat=(double)cadence;
  printf("cadence %lf seconds\n", deltat);

  TIME t_0, t_stamp;
  int focus,camid,fsns;
  
  TIME   interntime=sscan_time("2010.03.18_22:12:17.77_UTC");
  int axisbad[1];


  int i,j,k,km1,c,ki,iii,jjj,kkk;                                                   //loop variables
  int signid=0;
  int nfr;
 
 //********************************************************************************************************************************/
  //define and initialize quantities
  //********************************************************************************************************************************/



  //sigma for each FID (cosmic rays)

 
  //  double *sigma=(double *)(malloc(nx*ny*sizeof(double)));
  float *limw=(float *)(malloc(nx*ny*sizeof(float)));


  //DRMS_arrays for output arrays
  DRMS_Array_t *arr_flat, *arr_offpoint;
  DRMS_Array_t *arrout_new;
  DRMS_Array_t *arrout, *arrout_val, *arrout_sig;

  int *cosmic_ray_data;
  float *val_data, *sig_data;

  //Array pointers for storage
 
  short *bad;
  short *badpix;
  double *flati;
 
  float *apod;

  float *flat_off, *flat_relative;
  float *app_flat;

  app_flat=(float *)(malloc(nx*ny*sizeof(float)));
  flat_off=(float *)(malloc(nx*ny*sizeof(float)));
  flat_relative=(float *)(malloc(nx*ny*sizeof(float)));
  badpix=(short *)(malloc(nx*ny*sizeof(short)));
  flati=(double *)(malloc(nx*ny*sizeof(double)));
 
  apod=(float *)(malloc(nx*ny*sizeof(float)));

  float *flatfield_new;
  float *fflat;
  
 

 

  double flat[nx][ny];
  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){flat[i][j]=0.0;}
 


  double *rhsp, *rhsm;
  double *limb_dark;
  int *count_p, *count_m;

  int count, ccount,ccount1,ccount2, count_filtergram; // counters 
  int count_flatfields;

  int last, current, next, holder;
  int index_last, index_current, index_next;
  float wlr;
  float *rr, *a1,*a2, *w1,*w2,*cmc;
  int flatfound, flatrec_index;

  rr=(float *)(malloc(nx*ny*sizeof(float)));
  a1=(float *)(malloc(nx*ny*sizeof(float)));
  a2=(float *)(malloc(nx*ny*sizeof(float)));
  w1=(float *)(malloc(nx*ny*sizeof(float)));
  w2=(float *)(malloc(nx*ny*sizeof(float)));
  cmc=(float *)(malloc(nx*ny*sizeof(float)));
  

  rhsp=(double *)(malloc(nx*ny*sizeof(double)));
  rhsm=(double *)(malloc(nx*ny*sizeof(double)));
         
  count_p=(int *)(malloc(nx*ny*sizeof(int)));
  count_m=(int *)(malloc(nx*ny*sizeof(int)));

  limb_dark=(double *)(malloc(nx*ny*sizeof(double)));

  //cosmic ray
  
  
  float arrimg[3][ny][nx];
  int *cosmicarr;
  cosmicarr=(int *)(malloc(nx*ny*sizeof(int)));
  float *cosmicval;
  cosmicval=(float *)(malloc(nx*ny*sizeof(float)));
  float *cosmicsig;
  cosmicsig=(float *)(malloc(nx*ny*sizeof(float)));
  float diff_forward, diff_backward;

  int fsnarr[3]={0 , 0, 0};

  unsigned char *cmarr=(unsigned char *)(malloc(nx*ny*sizeof(unsigned char)));

  char *query_flat, *query_flat_ref, *query_flat_relative;
  char flat_default[256]={"none"};
 
  double *param;
  param=(double *)(malloc(6*sizeof(double)));

  
  //parameters to identify records
  TIME tobs_link[2];
  long long recnum_offpoint, recnum_badpix;
  struct rotpar rot_cur;
  struct rotpar rot_new;

  int rotpairs;


  printf("START!\n");

  /***********************************************************************************************************/
  /*CHECK WHETHER THE FLATFIELD OUTPUT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
    
       drms_series_exists(drms_env, filename_flatfield_fid, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",filename_flatfield_fid);       //if the output series does not exit
	  return 1;                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Output series %s exists.\n",filename_flatfield_fid);
	}

  //***********************************************************************************************************/
  /*CHECK WHETHER THE COSMIC RAY OUTPUT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
    
      drms_series_exists(drms_env, cosmic_ray_series, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Cosmic ray series %s doesn't exist\n",cosmic_ray_series);       //if the output series does not exit
	  return 1;                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Cosmic ray series %s exists.\n",cosmic_ray_series);
	}
  //**************************************************************************************************************/


      //**********************************
      //parallelization
      //**********************************


      int nthreads; 
      nthreads=omp_get_max_threads();
      //nthreads=omp_get_num_procs();                                      //number of threads supported by the machine where the code is running
      //omp_set_num_threads(nthreads);                                     //set the number of threads to the maximum value
      printf("Number of threads run in parallel = %d \n",nthreads);

  
  //**************************************************************************************************************/
  //built query strings
  //**************************************************************************************************************/
            
 char timefirst[256]="";
 strcat(timefirst, datum);
 strcat(timefirst, "_00:00:00.00_TAI");

 char timelast[256]="";
 strcat(timelast, datum);
 strcat(timelast, "_23:59:59.99_TAI");

 int pad;
 pad=(int)(time_limit/filtergram_cadence+1.0);


 TIME tfirst, tlast;

 char query0[256]="";

if (fsn_first == 0 || fsn_last == 2147483647)
  {

    tfirst=sscan_time(timefirst);
    tlast=sscan_time(timelast);

    char datefirst[256]="";
    sprint_ut(datefirst, tfirst-time_limit);

    char datelast[256]="";
    sprint_ut(datelast, tlast+time_limit);


 strcat(query0, inRecQuery);
 strcat(query0, "[");
 strcat(query0, datefirst);
 strcat(query0, "-");
 strcat(query0, datelast);
 strcat(query0, "][?FID=");
 char ffnumb[2]={""};
 sprintf(ffnumb, "%5.5d", fid);
 strcat(query0, ffnumb);
 strcat(query0, "?][?CAMERA=");
 char fnumb[2]={""};
 sprintf(fnumb, "%1.1d", cam+1);
 strcat(query0, fnumb);
 strcat(query0, "?]");

  printf("query string: %s\n", query0);

  }



  if (fsn_first != 0 && fsn_last != 2147483647)
    {

      tfirst=0;
      tlast=3881520000.0;

      strcat(query0, inRecQuery);
      strcat(query0, "[][");
       char fsnf[10]={""};
       sprintf(fsnf, "%d", fsn_first-pad);
       strcat(query0, fsnf);
      strcat(query0, "-");
      char fsnl[10]={""};
      sprintf(fsnl, "%d", fsn_last+pad);
       strcat(query0, fsnl);
      strcat(query0, "][?FID=");
      char ffnumb[2]={""};
     sprintf(ffnumb, "%5.5d", fid);
      strcat(query0, ffnumb);
      strcat(query0, "?][?CAMERA=");
      char fnumb[2]={""};
     sprintf(fnumb, "%1.1d", cam+1);
     strcat(query0, fnumb);
      strcat(query0, "?]");
      printf("query string: %s\n", query0);  
    }

  //************************************************************************


  
  //*************************************************************************
  //open records and stage data
  //*************************************************************************


  
      data     = drms_open_records(drms_env,query0,&stat);
      if (data == NULL || stat != 0 || data->n == 0){printf("can not open records\n"); return 1;}

      drms_stage_records(data, 1, 0); // stage needed records from tape and wait until ready.


      //************************************************************************
      //define keyword arrays
      //************************************************************************


  int nRecs=data->n;
  printf("number of records %d\n", nRecs);




  DRMS_Record_t *rec0[nRecs];
 
  //array for keyvalues
  int keyvalue_cam[nRecs];
  int keyvalue_wl[nRecs];
  int keyvalue_pl[nRecs];
  DRMS_Type_Value_t  keyvalue_time;
  long tmind[nRecs+1];
  int keyvalue_fid[nRecs];
  char *keyvalue_iss[nRecs];
  char *keyvalue_flatnumb[nRecs];
  int keyvalue_flatr[nRecs];
  int keyvalue_fsn[nRecs];
  double time_fl[nRecs];
  //int cosmic_ray_check[nRecs];

  //arrays for Lev 1 solar parameters   
  float keyvalue_rsun[nRecs];
  double keyvalue_dsun_obs[nRecs];
  float keyvalue_p0[nRecs];
  float keyvalue_b0[nRecs];
  float keyvalue_X0[nRecs];
  float keyvalue_Y0[nRecs];
  float keyvalue_vrad[nRecs];

  

  int present[nRecs];        //indices
  int present_forward[nRecs];
  int present_backward[nRecs];
  int index[nRecs];
  int statarr[nRecs];

  //********************************************************************************


  //****************************************************
  //get time reference
  //****************************************************

   
 
  t_0 = drms_getkey_time(data->records[nRecs/2],keytobs,&status);
  fsns=drms_getkey_int(data->records[nRecs/2],keyfsn,&status);
  focus = drms_getkey_int(data->records[nRecs/2],keyfocus,&status);
  
  char tmstr[256]="";
  strcat(tmstr, datum);
  strcat(tmstr, "_00:00:00.00_TAI");

  t_stamp=sscan_time(tmstr)+24.0*60.0*60.0; 

 
  //************************************************************


  int update_stat=2;   // 2: not enough data // 1: no update required // 0: flatfield updated 
    

  //*************************************************************
  //do camera test
  //************************************************************

  if (stat == DRMS_SUCCESS && data != NULL && nRecs > 0)
    {

      nfr=0;
      for (k=0; k<nRecs; ++k)
	{
	  statarr[k]=0;
	  rec0[k]=data->records[k];
	  keyvalue_cam[k]=drms_getkey_int(rec0[k],keycam,&status);
	  if (status != 0 || (keyvalue_cam[k] != cam_id_front && keyvalue_cam[k] != cam_id_side)) statarr[k]=statarr[k]+1024;
	  ++nfr;
	}

     
     
   
     
      
   

  
      

      //******************************************************** 
      //read records and test keyword range
      //********************************************************
     

      for (k=0; k<nRecs; ++k){
        keyvalue_time=drms_getkey(rec0[k], keytobs, &type_time, &status);       //read "T_OBS" of filtergram
        time_fl[k] = keyvalue_time.time_val;
        tmind[k]=(long)time_fl[k]-(long)interntime;
	
      
                     //tmind: time since first filtergram (in s) 

	keyvalue_fsn[k]=drms_getkey_int(rec0[k], keyfsn, &status);

	keyvalue_wl[k]=drms_getkey_int(rec0[k],keywl,&status);               //wavelength ID
	keyvalue_pl[k]=drms_getkey_int(rec0[k],keypl,&status);       // Polarization ID
      

	//test level_1 keywords


	keyvalue_rsun[k]=drms_getkey_float(rec0[k], lev1_r_sun,&status); // /drms_getkey_float(rec0[k], lev1_imsc,&status0);
	if (keyvalue_rsun[k] < rsun_min || keyvalue_rsun[k] > rsun_max ||  isnan(keyvalue_rsun[k]))statarr[k]=statarr[k]+32; 

	keyvalue_dsun_obs[k]=drms_getkey_double(rec0[k], lev1_dist ,&status);
	if (status !=0 || keyvalue_dsun_obs[k] < dsun_obs_min || keyvalue_dsun_obs[k] > dsun_obs_max || isnan(keyvalue_dsun_obs[k])) statarr[k]=statarr[k]+16;

	keyvalue_p0[k]=drms_getkey_float(rec0[k], lev1_p0, &status); 
	if (status !=0 || keyvalue_p0[k] < p0_min || keyvalue_p0[k] > p0_max || isnan(keyvalue_p0[k])) statarr[k]=statarr[k]+4;

	keyvalue_b0[k]=drms_getkey_float(rec0[k], lev1_b0, &status);
	if (status !=0 || keyvalue_b0[k] < b0_min || keyvalue_b0[k] > b0_max ||  isnan(keyvalue_b0[k])) statarr[k]=statarr[k]+8;

	keyvalue_X0[k]=drms_getkey_float(rec0[k], lev1_x0, &status1);
	keyvalue_Y0[k]=drms_getkey_float(rec0[k], lev1_y0,&status2);
	if (status1 !=0  || status2 != 0 || keyvalue_X0[k] < X0_min || keyvalue_X0[k] > X0_max || isnan(keyvalue_X0[k]) || keyvalue_Y0[k] < Y0_min || keyvalue_Y0[k] > Y0_max || isnan(keyvalue_Y0[k])) statarr[k]=statarr[k]+64;


	keyvalue_fid[k]=drms_getkey_int(rec0[k],fidkey,&status);
	if (status !=0 || keyvalue_fid[k] < minfid || keyvalue_fid[k] > maxfid || isnan(keyvalue_fid[k])) statarr[k]=statarr[k]+256;
	keyvalue_iss[k]=drms_getkey_string(rec0[k],isskey,&status);
	if (status !=0 || strcmp(keyvalue_iss[k], "CLOSED") != 0) statarr[k]=statarr[k]+1;

	keyvalue_flatnumb[k]=drms_getkey_string(rec0[k],flatnkey,&status); 

	
	keyvalue_vrad[k]=drms_getkey_float(rec0[k], lev1_vr, &status);
	if (keyvalue_vrad[k] < vrad_min || keyvalue_vrad[k] > vrad_max) statarr[k]=statarr[k]+512;

   

      }
  

      //******************************/ 
      //calculate average for solar parameters /
      //******************************/


      float R_SUN=0.0;
      float XX0=0.0;
      float YY0=0.0;
      float P_ANG=0.0;
      float B_ANG=0.0;
      float dist=0.0;

      float dcount=0.0;
      for (k=0; k<nRecs; ++k){
	//printf("rsun %f %f %f %f %f %lf %d %d\n", keyvalue_rsun[k], keyvalue_X0[k], keyvalue_Y0[k], keyvalue_p0[k], keyvalue_b0[k], keyvalue_dsun_obs[k], statarr[k], keyvalue_fid[k]);
	if (statarr[k] == 0){

	R_SUN=R_SUN+keyvalue_rsun[k];
	XX0=XX0+keyvalue_X0[k];
	YY0=YY0+keyvalue_Y0[k];
	P_ANG=P_ANG-keyvalue_p0[k]/180.0f*M_PI;
	B_ANG=B_ANG+keyvalue_b0[k]/180.0f*M_PI;
	dist=dist+(float)keyvalue_dsun_obs[k]/oneau;
	dcount=dcount+1.0f;
	}
      }
    


      //
  
      if (dcount > 0.0){
	R_SUN=R_SUN/dcount;
	XX0=XX0/dcount;
	YY0=YY0/dcount;
	P_ANG=P_ANG/dcount;
	B_ANG=B_ANG/dcount;
	dist=dist/dcount;
     
       *(param+0)=(double)R_SUN;  //R_SUN in pixel
       *(param+1)=(double)XX0;    //center x in pixel
       *(param+2)=(double)YY0;    //center y in pixel
       *(param+3)=(double)P_ANG;  //P-angle in rad // 
       *(param+4)=(double)B_ANG;  //B-angle in rad
       *(param+5)=(double)dist;   //distance in AU
      } 
      else
	{
	  printf("Invalid lev1 keywords\n");
	  return 1;
	}
	
     
      printf("disk center and radius: %f %f %f\n", XX0, YY0,R_SUN);	   	   	;
  

      //**************************************************************
      //Start flatfield calculations
      //**************************************************************

      //************************
      //get valid frames for flatfield calculation
      //************************

      if (cam ==0)camid=cam_id_side;
      if (cam ==1)camid=cam_id_front;



      count_flatfields=0;

      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){rhsm[j*nx+i]=0.0; rhsp[j*nx+i]=0.0; count_p[j*nx+i]=0; count_m[j*nx+i]=0; flati[j*nx+i]=0.0; cmarr[j*nx+i]=1;}

	
	for (i=0; i<nRecs; ++i){present_backward[i]=0; present_forward[i]=0; present[i]=0;}
	//for (i=0; i<nRecs; ++i){present[i]=0; index[i]=-1;}

	  //count filtergrams, and number of valid pairs
	  //count_filtergram=0; 
   
	  for (k=0; k<nRecs; ++k){present[k]=1; index[k]=k;}  //loop over all filtergrams: search for valid filtergrams
	  printf("number of filtergrams %d\n", nRecs);
 	  dataout = drms_create_records(drms_env,nRecs,(char *)cosmic_ray_series,DRMS_PERMANENT,&stat); //create ALL cosmic_ray records

	//check each forward and backward pair of images for identity
	//check for wavelength, polarization,  P-angle, center, and solar radius identity, and for time difference being nominal cadence

	
	  count_filtergram=nRecs;

	  //if (count_filtergram >=3){  //if at least 3 filtergrams of this type

	  
	  printf("filtergram with FID  %d: %d\n", fid, count_filtergram);
	for (c=1; c<count_filtergram; ++c){

//backward pairs
	  present_backward[c]=present[c];
	  if (statarr[c] != 0 || statarr[c-1] != 0 || keyvalue_wl[c] != keyvalue_wl[c-1] || keyvalue_pl[c] != keyvalue_pl[c-1] || (tmind[c]-tmind[c-1]) != (long)deltat || fabs(keyvalue_p0[c]- keyvalue_p0[c-1]) > 0.2 || sqrt(pow(keyvalue_X0[c]-keyvalue_X0[c-1],2)+pow(keyvalue_Y0[c]-keyvalue_Y0[c-1],2)) > limit_centerdiff[cam] || fabs(keyvalue_rsun[c]-keyvalue_rsun[c-1]) > limit_rsundiff){present_backward[c]=0;}

	  if (sqrt(pow(keyvalue_X0[c]-keyvalue_X0[c-1],2)+pow(keyvalue_Y0[c]-keyvalue_Y0[c-1],2)) > limit_centerdiff_cosmic) statarr[c]=statarr[c]+128;
	  //  printf("%d %d %ld %f %f\n", keyvalue_wl[index[c]]-keyvalue_wl[index[c-1]],keyvalue_pl[index[c]]- keyvalue_pl[index[c-1]], (tmind[index[c]]-tmind[index[c-1]]), fabs(keyvalue_p0[index[c]]- keyvalue_p0[index[c-1]]), fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c-1]]));
	}


	//forward pairs
	ccount=count_filtergram-1; //number of pairs (forward or backward)
	for (c=0; c<(count_filtergram-1); ++c){
	  present_forward[c]=present[c];
	  if (statarr[c] != 0 || statarr[c+1] != 0 || keyvalue_wl[c] != keyvalue_wl[c+1] || keyvalue_pl[c] != keyvalue_pl[c+1] || (tmind[c+1]-tmind[c]) != (long)deltat || fabs(keyvalue_p0[c+1]- keyvalue_p0[c]) > 0.2 || sqrt(pow(keyvalue_X0[c]-keyvalue_X0[c+1],2)+pow(keyvalue_Y0[c]-keyvalue_Y0[c+1],2)) > limit_centerdiff[cam] || fabs(keyvalue_rsun[c]-keyvalue_rsun[c+1]) > limit_rsundiff){present_forward[c]=0; --ccount;}

	  if (sqrt(pow(keyvalue_X0[c]-keyvalue_X0[c+1],2)+pow(keyvalue_Y0[c]-keyvalue_Y0[c+1],2)) > limit_centerdiff_cosmic || fabs(keyvalue_X0[c]-XX0) > limit_offpoint || fabs(keyvalue_Y0[c]-YY0) > limit_offpoint) statarr[c]=statarr[c]+128;
	}

	printf("number of filtergrams, valid pairs of frames %d \t %d \n", count_filtergram, ccount);



	//	for (k=0; k < count_filtergram; ++k) printf("%d \t %d \t %s \t %d\t %d \t %d \t %ld   \t %d  \t %f \t %f \t  %f %f %d %d %d %d\n", k, keyvalue_fsn[k], keyvalue_iss[k], keyvalue_wl[k], keyvalue_pl[k], keyvalue_cam[k], tmind[k], keyvalue_fid[k], keyvalue_rsun[k], keyvalue_X0[k], keyvalue_Y0[k], keyvalue_vrad[k], statarr[k], present[k], present_forward[k], present_backward[k]);



	//*******************************************
	///read relative flatfield
	//********************************************

	status_flatfield_relative=get_latest_flat(t_stamp, cam+1, filename_offpoint, segmentname_offpoint, flat_relative, &recnum_offpoint, &focus);




	 if (status_flatfield_relative == 0)
	   {
	     printf("Relative flatfield read: record number:%ld\n", recnum_offpoint);
	   }
	 else 
	   {
	     printf("WARNING: Could not read relative flatfield");
	   }

	 ///////////////////


	 
      


	 //******************************************
	 //calculate limit parameters for cosmic ray detection
	 //******************************************


#pragma omp parallel for private(jjj,iii)	   
	      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
		{
				  
		  rr[jjj*nx+iii]=sqrt(((float)iii-XX0)*((float)iii-XX0)+((float)jjj-YY0)*((float)jjj-YY0))/R_SUN;
		  a1[jjj*nx+iii]=coef0[cam][0]+coef0[cam][1]*rr[jjj*nx+iii]+coef0[cam][2]*rr[jjj*nx+iii]*rr[jjj*nx+iii];
		  a2[jjj*nx+iii]=coef1[cam][0]+coef1[cam][1]*rr[jjj*nx+iii]+coef1[cam][2]*rr[jjj*nx+iii]*rr[jjj*nx+iii];
		  w1[jjj*nx+iii]=coef2[cam][0]+coef2[cam][1]*rr[jjj*nx+iii]+coef2[cam][2]*rr[jjj*nx+iii]*rr[jjj*nx+iii];
		  w2[jjj*nx+iii]=w1[jjj*nx+iii];
		  cmc[jjj*nx+iii]=coef4[cam][0]+coef4[cam][1]*rr[jjj*nx+iii]+coef4[cam][2]*rr[jjj*nx+iii]*rr[jjj*nx+iii];

		}







	  //**********************************************************************************
	  //loop over filtergrams: obtain cosmic ray records and get input arrays for rotational flatfield calculations
	  //**********************************************************************************

	  printf("cosmic ray detection\n");

	  ccount1=0;
	  ccount2=0;
	  last=2; current=0; next=1;
	  double avg=0.0;
	  query_flat_ref=flat_default;
	  printtime();

	  printf("count filtergram %d\n", count_filtergram);

	  for (kkk=0; kkk<(count_filtergram+1); ++kkk){ // loop starts at 0
	  
	    if (cosmic_flag || ccount >= cthreshold[cam])
	      {

		km1=kkk-1;

		//**************************
		//read in filtergram

		if (kkk < count_filtergram)
		  {
		    printf("fsn %d %d\n", kkk, keyvalue_fsn[kkk]);
	 
		    segin    = drms_segment_lookupnum(rec0[kkk], 0);
		    arrin0= drms_segment_read(segin, type_double, &status);
	    

	  
		    if (status != DRMS_SUCCESS || arrin0 == NULL)
		      {
	    
			printf("Error: there is a problem with the filtergram number %d \n", kkk); //if there is a problem with the filtergram
			--ccount;

			present_backward[kkk]=0;
			if (kkk < (count_filtergram-1)) present_backward[kkk+1]=0;

			present_forward[kkk]=0;
			if (kkk > 0 ) present_forward[kkk-1]=0;

			present[kkk]=0;
			for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)   arrimg[next][jjj][iii]=NAN;
			flatfound=1;
		      }
		    else

		      {
	     
			printf("filtergram with FSN %d read\n",keyvalue_fsn[kkk]);
			arrinL0 = arrin0->data;


			//**************************************
			//flatfield query

			query_flat =keyvalue_flatnumb[kkk];
			if (strcmp(query_flat, query_flat_ref) != 0)
			  {
			    printf("%d %d read new flatfield\n", k, keyvalue_fsn[kkk]);

		 
			    status_flatfield=get_flat(query_flat, cam+1, app_flat, flat_off,  badpix);
			    flatfound=status_flatfield;

			    if (status_flatfield ==0)
			      {
				printf("Flatfield %s read\n", query_flat);

			      }
			    else

			      {
				printf("Failure reading flatfield / Use unity flatfield\n");

				status=get_latest_bad(t_stamp, cam+1, filename_badpix, segmentname_badpix, badpix, &recnum_badpix); //retrieve bad for links a
				if (status != 0){printf("could not find bad pixel list"); return 1;}
				else {printf("recnum badpix %ld\n", recnum_badpix);}
			      }
		 


			    query_flat_ref=query_flat;
			  }

			if (status_flatfield != 0)
			  {
			    --ccount;
			    statarr[kkk]=statarr[kkk]+2;

			    present_forward[kkk]=0;
			    if (kkk > 0 ) present_forward[kkk-1]=0;

			    present_backward[kkk]=0;
			    if (kkk < (count_filtergram-1)) present_backward[kkk+1]=0;
			    printf("no flatfield for filtergram %d\n", kkk);
			  }
	
	     
			if (flatfound == 0) for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)   arrimg[next][jjj][iii]=(float)arrinL0[jjj*nx+iii]*app_flat[jjj*nx+iii]/flat_relative[jjj*nx+iii];      	    
			if (flatfound != 0) for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)   arrimg[next][jjj][iii]=(float)arrinL0[jjj*nx+iii];

			fsnarr[next]=keyvalue_fsn[kkk];

			keyvalue_flatr[kkk]=flatfound;
		      }
		    drms_free_array(arrin0);
		  }


	    
		count=-1;

		//**********************************************
		//cosmic ray detection

		if (kkk >= 2 && kkk < count_filtergram)
		  {
		    index_last=kkk-2;
		    index_current=kkk-1;
		    index_next=kkk;

		    if (statarr[index_current] < 4)
		      {		

			flatrec_index=keyvalue_flatr[index_last]+keyvalue_flatr[index_current]+keyvalue_flatr[index_next]; //==0 if all three flatfields could be retrieved, ==3 if none has been retrieved

	
	
			if (statarr[index_last] < 4 && statarr[index_next] < 4 && (tmind[index_current]-tmind[index_last]) < time_limit && (tmind[index_next]-tmind[index_current]) < time_limit && (flatrec_index == 0 || flatrec_index == 3))
			  {
			    count=0;
#pragma omp parallel 
			    {
#pragma omp for private(jjj,iii,wlr)	   
			      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
				{
	
				  if (rr[jjj*nx+iii] < rad_cosmic_ray)
				    {
				      wlr=((float)(keyvalue_wl[km1] % 20)-5.0)/2.0*lambda_sep-keyvalue_vrad[km1]/v_c*lambda0-(cos(keyvalue_p0[km1]/180.*M_PI)*((float)iii-keyvalue_X0[km1])-sin(keyvalue_p0[km1]/180.*M_PI)*((float)jjj-keyvalue_Y0[km1]))*cos(keyvalue_b0[km1]/180.0*M_PI)/keyvalue_rsun[km1]*cpa.rotcoef0*radsun_mm/v_c*lambda0;
				      limw[jjj*nx+iii]=a1[jjj*nx+iii]*exp(-(wlr-coef3[cam][0])*(wlr-coef3[cam][0])/2.0/w1[jjj*nx+iii]/w1[jjj*nx+iii])+a2[jjj*nx+iii]*exp(-(wlr-coef3[cam][1])*(wlr-coef3[cam][1])/2.0/w2[jjj*nx+iii]/w2[jjj*nx+iii])+cmc[jjj*nx+iii];
				    }
				}


	
		    
	  //head= NULL;
#pragma omp for private(iii,jjj,diff_forward,diff_backward)
			      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
				{
				  cmarr[jjj*nx+iii]=1;
	   
				  if (rr[jjj*nx+iii] < rad_cosmic_ray)
				    {
				      diff_backward=arrimg[current][jjj][iii]-arrimg[last][jjj][iii];
				      diff_forward=arrimg[current][jjj][iii]-arrimg[next][jjj][iii];

				      if (diff_forward < 0.0 && diff_backward < 0.0) signid=1; else signid=0;

				      if (fabs(diff_forward) > (factor[cam][signid]*limw[jjj*nx+iii]) && fabs(diff_backward) > (factor[cam][signid]*limw[jjj*nx+iii])){
	     		    
					if (badpix[jjj*nx+iii]){
					  if (!isnan(arrimg[last][jjj][iii]) && !isnan(arrimg[current][jjj][iii]) && !isnan(arrimg[next][jjj][iii])){
#pragma omp critical(detect)
					    {
					      cosmicarr[count]=jjj*nx+iii;
					      if (diff_forward > 0.0 && diff_backward > 0.0)
						{
						  cosmicval[count]=minval(diff_forward, diff_backward);
						  cosmicsig[count]=fabs(cosmicval[count])/limw[jjj*nx+iii];
						  ++count; //count number of cosmic rays (using first differences)
						  cmarr[jjj*nx+iii]=0;
						}
					      if (diff_forward < 0.0 && diff_backward < 0.0)
						{
						  cosmicval[count]=maxval(diff_forward, diff_backward);
						  cosmicsig[count]=fabs(cosmicval[count])/limw[jjj*nx+iii];
						  ++count; //count number of cosmic rays (using first differences)
						  cmarr[jjj*nx+iii]=0;
						}

					      

					      

					      
					    }
					  }
					}
		    
				      }
				    }
				}
		   
			    }
			  }
		      }
	
		  }
	    
	      
	      //cosmic_ray_check[index_current]=count;
	    
	




	  /////////////////////////////////////////////////////////
	  //write out cosmic ray series
	  /////////////////////////////////////////////////////////

		//	if (cosmic_flag && kkk == 0){
		//	  if (present[kkk] && keyvalue_fsn[kkk] >= fsn_first && keyvalue_fsn[kkk] <= fsn_last && time_fl[kkk] >= tfirst && time_fl[kkk] <= tlast)
		//	    {
		//	      printf("anormal writeout of first record\n");
		//	      recout = dataout->records[kkk];
		//
		//		status=drms_setkey_int(recout, keyfsn, keyvalue_fsn[kkk]);
		//		status=drms_setkey_time(recout, keytobs, time_fl[kkk]);
		//		status=drms_setkey_int(recout, keycount, -1);
	       //		status=drms_setkey_int(recout, fidkey, fid);
		//		status=drms_setkey_int(recout, keycamera, keyvalue_cam[kkk]);
		//		status=drms_setkey_int(recout, keyexmax, 0);
		//	status=drms_setkey_float(recout, keylimit, factor[cam][0]);
	
		//	drms_keyword_setdate(recout);
	    
		//		if (keyvalue_cam[kkk] == cam_id_front) status=drms_setkey_string(recout, keyinstrument, camera_str_front);
		//	if (keyvalue_cam[kkk] == cam_id_side) status=drms_setkey_string(recout, keyinstrument, camera_str_side);
		//
		//	segout = drms_segment_lookup(recout, segmentname_cosmic);
			//	segout_val=drms_segment_lookup(recout, segmentname_val);
		//	segout_sig=drms_segment_lookup(recout, segmentname_sig);
	    
	 

		//	axisbad[0]=0;
		//	arrout=drms_array_create(type_int,1,axisbad,NULL,&status);
		//	cosmic_ray_data=arrout->data;
	  	 
		//	arrout_val=drms_array_create(type_float,1,axisbad,NULL,&status);
		//	val_data=arrout_val->data;

		//	arrout_sig=drms_array_create(type_float,1,axisbad,NULL,&status);
		//	sig_data=arrout_sig->data;

		//	status=drms_segment_write(segout, arrout, 0);
		//	status=drms_segment_write(segout_val, arrout_val, 0);
		//	status=drms_segment_write(segout_sig, arrout_sig, 0);
	   
		//	drms_free_array(arrout);
		//	drms_free_array(arrout_val);
		//	drms_free_array(arrout_sig);
		//   }
		//	}



		if (cosmic_flag && kkk >= 1){
		  
		 
		  if (present[km1] && keyvalue_fsn[km1] >= fsn_first && keyvalue_fsn[km1] <= fsn_last && time_fl[km1] >= tfirst && time_fl[km1] <= tlast)
		    {
		      printf("fsn count %d %d\n", keyvalue_fsn[km1], count);
		      printf("fsn order %d %d %d\n", fsnarr[last], fsnarr[current],fsnarr[next]);

		      int cosfind, limit_flag;
		      if (count == -1) cosfind=0; else if (count < limit_cosmic){cosfind=count; limit_flag=0;} else {cosfind=limit_cosmic; limit_flag=1;}


			recout = dataout->records[km1];

			status=drms_setkey_int(recout, keyfsn, keyvalue_fsn[km1]);
			status=drms_setkey_time(recout, keytobs, time_fl[km1]);
			status=drms_setkey_int(recout, keycount, count);
			status=drms_setkey_int(recout, fidkey, fid);
			status=drms_setkey_int(recout, keycamera, cam+1);
			status=drms_setkey_int(recout, keyexmax, limit_flag);
			status=drms_setkey_float(recout, keylimit, factor[cam][0]);
	
			drms_keyword_setdate(recout);
	    
			if (keyvalue_cam[km1] == cam_id_front) status=drms_setkey_string(recout, keyinstrument, camera_str_front);
			if (keyvalue_cam[km1] == cam_id_side) status=drms_setkey_string(recout, keyinstrument, camera_str_side);
	    
			segout = drms_segment_lookup(recout, segmentname_cosmic);
			segout_val=drms_segment_lookup(recout, segmentname_val);
			segout_sig=drms_segment_lookup(recout, segmentname_sig);
	    
	 

			axisbad[0]=cosfind;
			arrout=drms_array_create(type_int,1,axisbad,NULL,&status);
			cosmic_ray_data=arrout->data;
	  	 
			arrout_val=drms_array_create(type_float,1,axisbad,NULL,&status);
			val_data=arrout_val->data;

			arrout_sig=drms_array_create(type_float,1,axisbad,NULL,&status);
			sig_data=arrout_sig->data;


			for (i=0; i<cosfind; ++i){cosmic_ray_data[i]=cosmicarr[i]; val_data[i]=cosmicval[i]; sig_data[i]=cosmicsig[i];}

			status=drms_segment_write(segout, arrout, 0);
			status=drms_segment_write(segout_val, arrout_val, 0);
			status=drms_segment_write(segout_sig, arrout_sig, 0);
	   
			drms_free_array(arrout);
			drms_free_array(arrout_val);
			drms_free_array(arrout_sig);
		    }
		}

	    
	    

	  



	   
	  //**************************************************************************************************
	  //permanent bad pixel calculation
		

		if (ccount >= cthreshold[cam] && kkk >= 1){
		  update_stat=1;

		  if (present_backward[km1] != 0) ++ccount1;
		  if (present_forward[km1] != 0) ++ccount2;

#pragma omp parallel
		  {
	 
		    if (present_backward[km1] != 0){
		    
#pragma omp for private(jjj,iii)
		      for (jjj=0; jjj<ny; ++jjj){
			for (iii=0; iii<nx; ++iii){
			  //if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0 && cmarr[jjj*nx+iii]){rhsp[jjj*nx+iii]=rhsp[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]);  // //add up frames for backward
			  if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0){rhsp[jjj*nx+iii]=rhsp[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]); //!!disable cosmic ray control
			    count_p[jjj*nx+iii]=count_p[jjj*nx+iii]+1;}
			}
		      }
		    }
	
	 
		    if (present_forward[km1] != 0){
		    
		#pragma omp for private(jjj,iii)
		      for (jjj=0; jjj<ny; ++jjj){
			for (iii=0; iii<nx; ++iii){
			  //if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0 && cmarr[jjj*nx+iii]){rhsm[jjj*nx+iii]=rhsm[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]);   // //add up frames for forward
			  if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0){rhsm[jjj*nx+iii]=rhsm[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]); //!!disable cosmic ray control
			    count_m[jjj*nx+iii]=count_m[jjj*nx+iii]+1;}
			}
		      }
		    }
		  }

	      
      	      
		}


	    
	    
	      }

	    //**************************************************	 
	    //exchange indices

	    holder=next;
	    next=last;
	    last=current;
	    current=holder;
	 


	  }//end loop over filtergrams of same type
	  
	  //************************************************************************************




	  printtime();	  
	  printf("loop done\n");
	  drms_close_records(dataout, DRMS_INSERT_RECORD);

	  printf("reference flatfield and limb darkening\n");
   

	  //*********************************
	  //calculate limb darkening function and correct for it

	  limb_darkening((double)R_SUN, (double)XX0, (double)YY0, b_coef, b_order, limb_dark);

	 
 #pragma omp parallel for private(jjj,iii)
	  for (jjj=0; jjj<ny; ++jjj){
	    for (iii=0; iii<nx; ++iii){
	       
	      if (count_p[jjj*nx+iii] > 0) rhsp[jjj*nx+iii]=(rhsp[jjj*nx+iii]/(double)count_p[jjj*nx+iii]-log(limb_dark[jjj*nx+iii])); 
	      if (count_m[jjj*nx+iii] > 0) rhsm[jjj*nx+iii]=(rhsm[jjj*nx+iii]/(double)count_m[jjj*nx+iii]-log(limb_dark[jjj*nx+iii])); 
	    }
	  }

 
	  //	  if (debug)  //
	  //{
	  //printf("write debug\n");
	  //printf("count_filtergram %d\n", count_filtergram);
          // FILE *fileptr;                                                                                                                          
          // int *aaa;                                                                                                                            
          // aaa=(int *)(malloc(count_filtergram*sizeof(int)));                                                                                              
          // fileptr = fopen ("/tmp20/richard/interpol/dd.bin", "w");                                                                                    
	   // for (i=0;i<count_filtergram;i++){aaa[i]=farr[i];} fwrite ((char*)(aaa),sizeof(int),count_filtergram,fileptr);
	   //  for (i=0;i<count_filtergram;i++){aaa[i]=barr[i];} fwrite ((char*)(aaa),sizeof(int),count_filtergram,fileptr);
          // fclose(fileptr);                                                                                                                            
	  //	free(aaa);
	  //}


      //*******************************************************************/
      //do flatfield calculations                                          /
      //*******************************************************************/

	  printf("forward / backwards pairs %d %d \n", ccount1, ccount2);
 
	  if (ccount1>= cthreshold[cam]){

  
      for (j=0; j<ny; ++j)  for (i=0; i<nx; ++i) flati[j*nx+i]=0.0;
 
      rotpairs=ccount1;

      cpa.norm=normconst[cam];
      cpa.omega=omegaconst[cam];
      cpa.convergence=convconst[cam];
      cpa.maxiter=maxiterconst[cam];

      if (flatfield_flag && status_flatfield_relative == 0)
	{
	  status=flatfield(rhsp, rhsm, badpix, ccount, flati, param, cpa, deltat);
	}
      else
	{
	  status=1;
	}
  

      if (status ==0)
	{
      //*********************************************************************//
      
	  for (j=0; j<ny; ++j)  for (i=0; i<nx; ++i)flat[i][j]=exp(flati[j*nx+i]);
	  update_stat=0;
   
      //*********************************************************************//
	}
	  }
	 


      if (status == 0) printf("flatfield for camera %d calculated \n", cam); else printf("could not calcuate flatfield for camera %d \n", cam);

    
     

      //*********************************************************************
      //write out  flatfield
      //********************************************************************

 


      int axisout[2]={nx,ny};                             //size of the output array
      arrout_new  = drms_array_create(type_float,2,axisout,NULL,&status);
      printf("array size %ld\n", drms_array_size(arrout_new)/sizeof(float)); 
  
      flatfield_new=arrout_new->data;
  
 

      //*************************
      //calculate apodization
             
       apod_circ(R_SUN*cpa.croprad*0.95, R_SUN*cpa.croprad*0.04, XX0-(float)nx/2.0, YY0-(float)ny/2.0, apod);

 
  
       //************************************
       //copy and apodize flatfield

for (j=0; j<ny; ++j){
      for (i=0; i<nx; ++i){
	flatfield_new[j*nx+i]=(flat[i][j]-1.0)*apod[j*nx+i]+1.0;
      }
 }



 if (update_stat == 0)rot_new.rotpairs=rotpairs; else rot_new.rotpairs=0;
 rot_new.rotcadence=(float)deltat;
    //rot_new.flatfield_version=rot_cur.flatfield_version;


    
    printf("update_stat %d\n", update_stat);

    printf("query flat %s\n", query_flat);

    if (flatfield_flag)
      {
	status=write_flatfield_fid(arrout_new, cam+1, recnum_offpoint, t_0, focus, fid, rot_new);;
	if (status != 0){printf("could not write rotational flatfield\n"); return 1;}
      }


	  drms_free_array(arrout_new);

	
    }
  else 
    {
      printf("No data records found\n");
      return 1;
    }     

    
  
  if (data != NULL)
    {
      printf("close DRMS session \n");
      drms_close_records(data,DRMS_FREE_RECORD);       
    }


  //************************************************************
  //free arrays
  //****************************************************


  free(param);
  free(rr);
  free(a1);
  free(a2);
  free(w1);
  free(w2);
  free(cmc);
  free(cosmicarr);
  free(cosmicsig);
  free(cosmicval);
  free(cmarr);
  free(apod);
	  free(flat_off);
	  free(flat_relative);
	  free(badpix);
	  free(app_flat);
	  free(flati);

	  free(rhsp);
	  free(rhsm);
	  free(count_p);
	  free(count_m);
	  free(limb_dark);

	  free(limw);


  printtime();  
  printf("COMPLETED!\n");

  
  
  return 0;

}








//valgrind --leak-check=full --track-origins=yes --show-reachable=yes module_flatfield_256 input_series="su_production.lev1_hmi_test[986331-986430]" cosmic_rays=1
