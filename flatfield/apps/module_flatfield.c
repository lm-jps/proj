


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
     {ARG_INT, kDSCosmic, "0", "Cosmic rays flag"},
     {ARG_INT, kDSFlatfield, "0", "Flatfield flag"},
     {ARG_INT, cadence_name, "0", "Cadence in sec"},
     {ARG_INT, fid_name, "0", "FID"},
     {ARG_INT, cam_name, "0", "Camera"},
     {ARG_INT, fsnf_name, "0"},
     {ARG_INT, fsnl_name, "0"},
     {ARG_STRING, datumn},
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

int get_flatfields(char *query_flat, TIME t_0, int camera, int focus,
                   float *flatfield, float *offpoint, 
                   short *bad, long long recnum[6], TIME tobs_link[2],
		   struct rotpar *rot_cur);


int write_flatfields(char *filename_flatfield_out, DRMS_Array_t *arr_flat, DRMS_Array_t *arrout_new, int camera,
                     long long recnum[6], TIME tobs_link[2], TIME t_0, int focus,
                     struct rotpar rot_new, 
                     struct rotpar rot_cur);

int write_flatfield_fid(DRMS_Array_t *arrout_new, int camera, long long recnum_off, char *query_flat, 
			TIME tobs_link[2], TIME t_0, int focus,int fid,
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

  const char *inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL); //cmdparams is defined in jsoc_main.h
  int cosmic_flag=cmdparams_get_int(&cmdparams, kDSCosmic, NULL);
  int flatfield_flag=cmdparams_get_int(&cmdparams, kDSFlatfield, NULL);
  int cadence=cmdparams_get_int(&cmdparams, cadence_name, NULL);

  int fid=cmdparams_get_int(&cmdparams, fid_name, NULL);
  int cam=cmdparams_get_int(&cmdparams, cam_name, NULL)-1;
  int fsn_first=cmdparams_get_int(&cmdparams, fsnf_name, NULL);
  int fsn_last=cmdparams_get_int(&cmdparams, fsnl_name, NULL);

  const char *datum=cmdparams_get_str(&cmdparams, datumn, NULL);
 

  if (fid < minfid || fid > maxfid){printf("Not an observable FID\n"); exit(EXIT_FAILURE);}

  int  status= DRMS_SUCCESS; 
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

  TIME t_0;
  int focus,camid,fsns;
  
  TIME   interntime;
  int axisbad[1];


  int i,j,k,km1,c,ki,iii,jjj,kkk;                                                   //loop variables



 
 //********************************************************************************************************************************/
  //Constants to be set: will maybe be input parameters 
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
  double *flatc;
  float *apod;

  float *flat_off, *flat_relative;
  float *app_flat;

  app_flat=(float *)(malloc(nx*ny*sizeof(float)));
  flat_off=(float *)(malloc(nx*ny*sizeof(float)));
  flat_relative=(float *)(malloc(nx*ny*sizeof(float)));
  badpix=(short *)(malloc(nx*ny*sizeof(short)));
  flati=(double *)(malloc(nx*ny*sizeof(double)));
  flatc=(double *)(malloc(nx*ny*sizeof(double)));
  apod=(float *)(malloc(nx*ny*sizeof(float)));

  float *flatfield_new;
  float *fflat;
  
 

 

  double flat[nx][ny];

 


  double *rhsp, *rhsm;
  double *limb_dark;
  int *count_p, *count_m;

  int count, ccount, count_filtergram, bad_pix; // counters 
  int count_flatfields;

  int last, current, next, holder;
  int index_last, index_current, index_next;
  float wlr;
  float *rr, *a1,*a2, *w1,*w2,*cmc;

  rr=(float *)(malloc(nx*ny*sizeof(float)));
  a1=(float *)(malloc(nx*ny*sizeof(float)));
  a2=(float *)(malloc(nx*ny*sizeof(float)));
  w1=(float *)(malloc(nx*ny*sizeof(float)));
  w2=(float *)(malloc(nx*ny*sizeof(float)));
  cmc=(float *)(malloc(nx*ny*sizeof(float)));
  //struct list *curr, *head;
  //struct list **list_pointer;

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



  unsigned char *cmarr=(unsigned char *)(malloc(nx*ny*sizeof(unsigned char)));

  
  
  //parameters to identify records
  TIME tobs_link[2];
  long long recnum[6];
  long long recnum_ref[6];
  struct rotpar rot_cur;
  struct rotpar rot_new;




  printf("START!\n");

  /***********************************************************************************************************/
  /*CHECK WHETHER THE FLATFIELD OUTPUT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
    
       drms_series_exists(drms_env, filename_flatfield_fid, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",filename_flatfield_fid);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Output series %s exists.\n",filename_flatfield_fid);
	}

  //***********************************************************************************************************/
  /*CHECK WHETHER THE COSMIC RAY OUTPUT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
    
       drms_series_exists(drms_env, filename_cosmic, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Cosmic ray series %s doesn't exist\n",filename_cosmic);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Cosmic ray series %s exists.\n",filename_cosmic);
	}
  //**************************************************************************************************************/


      //set parallelization parameters
      int nthreads=1; 
      nthreads=omp_get_num_procs();                                      //number of threads supported by the machine where the code is running
      omp_set_num_threads(nthreads);                                     //set the number of threads to the maximum value
      printf("Number of threads run in parallel = %d \n",nthreads);

  //********************************************************************************************************************************/

 
 
  //**************************************************************************************************************/
  //read records
  //**************************************************************************************************************/
            
 char timefirst[256]="";
 strcat(timefirst, datum);
 strcat(timefirst, "_00:00:00.00_TAI");

 char timelast[256]="";
 strcat(timelast, datum);
 strcat(timelast, "_23:59:59.99_TAI");

 TIME tfirst=sscan_time(timefirst)-deltat;
 TIME tlast=sscan_time(timelast)+2*deltat;

 char datefirst[256]="";
 sprint_ut(datefirst, tfirst);

 char datelast[256]="";
 sprint_ut(datelast, tlast);


 char query0[256]="";

if (fsn_first == 0 && fsn_last == 0)
  {

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



  if (fsn_first != 0 && fsn_last != 0)
    {

      strcat(query0, inRecQuery);
      strcat(query0, "[][");
       char fsnf[10]={""};
       sprintf(fsnf, "%d", fsn_first);
       strcat(query0, fsnf);
      strcat(query0, "-");
      char fsnl[10]={""};
      sprintf(fsnl, "%d", fsn_last);
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

  

  
      data     = drms_open_records(drms_env,query0,&stat);
      if (data == NULL || stat != 0 || data->n == 0){printf("can not open records\n"); exit(EXIT_FAILURE);}

      drms_stage_records(data, 0, 1); //added follwoing the recommendation of Phil

  int nRecs=data->n;
  printf("number of records %d\n", nRecs);

  int nfr;


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
  int keyvalue_fsn[nRecs];
  double time_fl[nRecs];
  int cosmic_ray_check[nRecs];

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

  char *query_flat, *query_flat_ref, *query_flat_relative;

  //double *sigmacoef;

  //sigmacoef=(double *)(malloc(5*sizeof(double)));
  //sigmacoef[0]=constsigma[(fid-10000)/10]*8.0; sigmacoef[1]=0.0; sigmacoef[2]=0.0; sigmacoef[3]=0.0; sigmacoef[4]=0.0;


  
 
 t_0 = drms_getkey_time(data->records[nRecs/2],keytobs,&status);
  fsns=drms_getkey_int(data->records[nRecs/2],keyfsn,&status);
  focus = drms_getkey_int(data->records[nRecs/2],keyfocus,&status);
  

 
       


  int update_stat=2;   // 2: not enough data // 1: no update required // 0: flatfield updated 
    

  if (stat == DRMS_SUCCESS && data != NULL && nRecs > 0)
    {

      nfr=0;
      for (k=0; k<nRecs; ++k)
	{
	  statarr[k]=0;
	  cosmic_ray_check[k]=-1;
	  rec0[k]=data->records[k];
	  keyvalue_cam[k]=drms_getkey_int(rec0[k],keycam,&status);
	  if (status != 0 || (keyvalue_cam[k] != cam_id_front && keyvalue_cam[k] != cam_id_side)) statarr[k]=statarr[k]+1024;
	  ++nfr;
	}

     
      int fvers;
          
     

   
       
 
   
     
      
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){flat[i][j]=0.0;}

  
      

      //******************************************************** 
      //read records
      //********************************************************
     
      printtime();

      for (k=0; k<nRecs; ++k){
        keyvalue_time=drms_getkey(rec0[k], keytobs, &type_time, &status);       //read "T_OBS" of filtergram
        time_fl[k] = keyvalue_time.time_val;
        tmind[k]=(long)time_fl[k]-(long)interntime;                           //tmind: time since first filtergram (in s) 

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

	keyvalue_X0[k]=drms_getkey_float(rec0[k], lev1_x0, &status);
	if (status !=0  || keyvalue_X0[k] < X0_min || keyvalue_X0[k] > X0_max || isnan(keyvalue_X0[k])) statarr[k]=statarr[k]+64;

	keyvalue_Y0[k]=drms_getkey_float(rec0[k], lev1_y0,&status);
	if (status !=0 || keyvalue_Y0[k] < Y0_min || keyvalue_Y0[k] > Y0_max || isnan(keyvalue_Y0[k])) statarr[k]=statarr[k]+128;
	keyvalue_fid[k]=drms_getkey_int(rec0[k],fidkey,&status);
	if (status !=0 || keyvalue_fid[k] < minfid || keyvalue_fid[k] > maxfid || isnan(keyvalue_fid[k])) statarr[k]=statarr[k]+256;
	keyvalue_iss[k]=drms_getkey_string(rec0[k],isskey,&status);
	if (status !=0 || strcmp(keyvalue_iss[k], "CLOSED") != 0) statarr[k]=statarr[k]+1;

	keyvalue_flatnumb[k]=drms_getkey_string(rec0[k],flatnkey,&status); 

	
	keyvalue_vrad[k]=drms_getkey_float(rec0[k], lev1_vr, &status);
	if (keyvalue_vrad[k] < vrad_min || keyvalue_vrad[k] > vrad_max) statarr[k]=statarr[k]+512;

	if (keyvalue_fsn[k] == 11986751) printf("%d \t %d \t %s \t %d\t %d \t %d \t %ld   \t %d  \t %f \t %f \t  %f %f %d\n", k, keyvalue_fsn[k], keyvalue_iss[k], keyvalue_wl[k], keyvalue_pl[k], keyvalue_cam[k], tmind[k], keyvalue_fid[k], keyvalue_rsun[k], keyvalue_X0[k], keyvalue_Y0[k], keyvalue_vrad[k], statarr[k]);

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
	//	printf("rsun %f %f %f %f %f %lf %d %d\n", keyvalue_rsun[k], keyvalue_X0[k], keyvalue_Y0[k], keyvalue_p0[k], keyvalue_b0[k], keyvalue_dsun_obs[k], statarr[k], keyvalue_fid[k]);
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
      double *param;
      param=(double *)(malloc(6*sizeof(double)));      

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
	  exit(EXIT_FAILURE);
	}
	
     
      //********************************

      int rotpairs;
          
   

  
 
  

      //**************************************************************
      //Start flatfield calculations
      //**************************************************************

     
 

	
	




		if (cam ==0)camid=cam_id_side;
		if (cam ==1)camid=cam_id_front;



	count_flatfields=0;

	
	//**************************************************************************************
        //loop over different types of filtergrams                                            // 
	//**************************************************************************************



  


	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){rhsm[j*nx+i]=0.0; rhsp[j*nx+i]=0.0; count_p[j*nx+i]=0; count_m[j*nx+i]=0; flati[j*nx+i]=0.0;}

	
	for (i=0; i<nRecs; ++i){present_backward[i]=0; present_forward[i]=0;}
	for (i=0; i<nRecs; ++i){present[i]=0; index[i]=-1;}

	  //count filtergrams, and number of valid pairs
	  count_filtergram=0; 
   
	  for (k=0; k<nRecs; ++k){present[k]=1; index[count_filtergram]=k; ++count_filtergram;}  //loop over all filtergrams: search for valid filtergrams
	  printf("number of filtergrams %d\n", count_filtergram);
 	  dataout = drms_create_records(drms_env,count_filtergram-2,filename_cosmic,DRMS_PERMANENT,&stat); //create ALL cosmic_ray records

	//check each forward and backward pair of images for identity
	//check for wavelength, polarization,  P-angle, center, and solar radius identity, and for time difference being nominal cadence

	

	  if (count_filtergram >=3){  //if at least 3 filtergrams of this type

	  
	  printf("filtergram with FID  %d: %d\n", fid, count_filtergram);
	for (c=1; c<count_filtergram; ++c){

//backward pairs
	  present_backward[index[c]]=present[index[c]];
	  if (statarr[index[c]] != 0 || statarr[index[c-1]] != 0 || keyvalue_wl[index[c]] != keyvalue_wl[index[c-1]] || keyvalue_pl[index[c]] != keyvalue_pl[index[c-1]] || (tmind[index[c]]-tmind[index[c-1]]) != (long)deltat || fabs(keyvalue_p0[index[c]]- keyvalue_p0[index[c-1]]) > 0.2 || sqrt(pow(keyvalue_X0[index[c]]-keyvalue_X0[index[c-1]],2)+pow(keyvalue_Y0[index[c]]-keyvalue_Y0[index[c-1]],2)) > limit_centerdiff[cam] || fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c-1]]) > limit_rsundiff){present_backward[index[c]]=0;}

	  if (sqrt(pow(keyvalue_X0[index[c]]-keyvalue_X0[index[c-1]],2)+pow(keyvalue_Y0[index[c]]-keyvalue_Y0[index[c-1]],2)) > limit_centerdiff_cosmic) statarr[index[c]]=statarr[index[c]]+64;
	  //  printf("%d %d %ld %f %f\n", keyvalue_wl[index[c]]-keyvalue_wl[index[c-1]],keyvalue_pl[index[c]]- keyvalue_pl[index[c-1]], (tmind[index[c]]-tmind[index[c-1]]), fabs(keyvalue_p0[index[c]]- keyvalue_p0[index[c-1]]), fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c-1]]));
	}


	//forward pairs
	ccount=count_filtergram-1; //number of pairs (forward or backward)
	for (c=0; c<(count_filtergram-1); ++c){
	  present_forward[index[c]]=present[index[c]];
	  if (statarr[index[c]] != 0 || statarr[index[c+1]] != 0 || keyvalue_wl[index[c]] != keyvalue_wl[index[c+1]] || keyvalue_pl[index[c]] != keyvalue_pl[index[c+1]] || (tmind[index[c+1]]-tmind[index[c]]) != (long)deltat || fabs(keyvalue_p0[index[c+1]]- keyvalue_p0[index[c]]) > 0.2 || sqrt(pow(keyvalue_X0[index[c]]-keyvalue_X0[index[c+1]],2)+pow(keyvalue_Y0[index[c]]-keyvalue_Y0[index[c+1]],2)) > limit_centerdiff[cam] || fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c+1]]) > limit_rsundiff){present_forward[index[c]]=0; --ccount;}

	  if (sqrt(pow(keyvalue_X0[index[c]]-keyvalue_X0[index[c+1]],2)+pow(keyvalue_Y0[index[c]]-keyvalue_Y0[index[c+1]],2)) > limit_centerdiff_cosmic) statarr[index[c]]=statarr[index[c]]+64;
	}

	printf("number of filtergrams, valid pairs of frames %d \t %d \n", count_filtergram, ccount);


	//	for (kkk=0; kkk<count_filtergram-2; ++kkk){k=index[kkk]; printf("%d %d %d %d %ld %d %d %d %f %ld\n",kkk,k,keyvalue_fid[k],keyvalue_cam[k],tmind[k]-tmind[0],present[k],present_backward[k],present_forward[k],sqrt(pow(keyvalue_X0[index[kkk]]-keyvalue_X0[index[kkk+1]],2)+pow(keyvalue_Y0[index[kkk]]-keyvalue_Y0[index[kkk+1]],2)), (tmind[index[kkk+1]]-tmind[index[kkk]])); if ((tmind[index[kkk+1]]-tmind[index[kkk]]) != 135.0) ++counter_aaa;}

      //combine flatfield correction and limb-darkening


	//get reference flatfield
	 query_flat_relative =keyvalue_flatnumb[nRecs-1];

	 status_flatfield_relative=get_flatfields(query_flat_relative, t_0, cam+1, focus,  app_flat, flat_relative,  badpix, recnum_ref, tobs_link, &rot_cur);

    
	 
	  //read in first filtergram
	printf("read first filtergram\n");

      

	  last=0; current=1; next=2;
	  index_last=index[0]; index_current=index[1];
	 
	  
	  segin    = drms_segment_lookupnum(rec0[index_last], 0);
	  arrin0= drms_segment_read(segin, type_double, &status);

	  if (status != DRMS_SUCCESS)
	    {
	      printf("Error: there is a problem with the filtergram number %d \n", k); //if there is a problem with the filtergram
	      present_backward[index_last]=0;
	      present_forward[index_last]=0;
	      present[index_last]=0;
	      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) arrimg[last][jjj][iii]=NAN; //
	    }
	  else
	    {
	     
	      printf("filtergram with FSN %d read\n",keyvalue_fsn[index_last]);
	      arrinL0 = arrin0->data;
	      
	      printf("%d %d read flatfield\n", k, keyvalue_fsn[k]);
	      query_flat =keyvalue_flatnumb[index_last];

	      status_flatfield=get_flatfields(query_flat, t_0, cam+1, focus,  app_flat, flat_off,  badpix, recnum, tobs_link, &rot_cur);

	      if (status_flatfield ==0)
		    printf("Flatfield read \n recnum offpoint %ld\n", recnum[0]);
		  else
		    printf("Failure reading flatfield / Use unity flatfield\n");

	      if (status_flatfield != 0)
		{
		  statarr[0]=statarr[0]+2;
		  present_forward[0]=0;
		  present_backward[0]=0;
		}

	      query_flat_ref=query_flat;
	      
		for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) arrimg[last][jjj][iii]=(float)arrinL0[jjj*nx+iii]*app_flat[jjj*nx+iii]/flat_relative[jjj*nx+iii]; 

	      k=index[0];

	    
	      
	
	       	      if (present_forward[k] != 0){ 
		
	         #pragma omp parallel for private(jjj,iii)
	        	for (jjj=0; jjj<ny; ++jjj){
	          for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsm[jjj*nx+iii]=rhsm[jjj*nx+iii]+log((double)arrimg[last][jjj][iii]); //add up first frame for forward
	            count_m[jjj*nx+iii]=count_m[jjj*nx+iii]+1;}
	      	  }
	      	}
		      }

	      drms_free_array(arrin0);
	    }

	  
	
//read in second filtergram
printf("read second filtergram\n");

	  segin    = drms_segment_lookupnum(rec0[index_current], 0);
	  arrin0= drms_segment_read(segin, type_double, &status);

	  
	  if (status != DRMS_SUCCESS)
	    {
	      printf("Error: there is a problem with the filtergram number %d \n", k); //if there is a problem with the filtergram
	      present_backward[index_current]=0;
	      present_forward[index_current]=0;
	      present[index_current]=0;
	    }
	  else
	    {
	      
	     
	      printf("filtergram with FSN %d read\n",keyvalue_fsn[index_current]);
	      arrinL0 = arrin0->data;

	      query_flat =keyvalue_flatnumb[index_current];
	      if (strcmp(query_flat, query_flat_ref) != 0)
		{
		 
		  printf("%d %d read new flatfield !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", k, keyvalue_fsn[k]);

		  status_flatfield=get_flatfields(query_flat, t_0, cam+1, focus,  app_flat, flat_off,  badpix, recnum, tobs_link, &rot_cur);
		  
		  if (status_flatfield ==0)
		    printf("Flatfield read \n recnum offpoint %ld\n", recnum[0]);
		  else
		    printf("Failure reading flatfield / Use unity flatfield\n");

		  query_flat_ref=query_flat;
		 
		}


		
	    
	      
	      if (status_flatfield != 0)
		    { 
		      statarr[1]=statarr[1]+2;
		      present_forward[1]=0;
		      present_backward[1]=0;
		    }


	      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) arrimg[current][jjj][iii]=(float)arrinL0[jjj*nx+iii]*app_flat[jjj*nx+iii]/flat_relative[jjj*nx+iii];

	      k=index[1];
	      
	      //printf("fsn %d %d\n", k, keyvalue_fsn[k]);

	      //	       if (present_backward[k] != 0){

	      //#pragma omp parallel for private(jjj,iii)
	      //		 for (jjj=0; jjj<ny; ++jjj){
	      //	      	  for (iii=0; iii<nx; ++iii){
	      //	      	    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsp[jjj*nx+iii]=rhsp[jjj*nx+iii]+log(arrinL0[jjj*nx+iii]); //  //add up frames (flatfielded with offpoint flatfield and limb darkening removed
	      //	            count_p[jjj*nx+iii]=count_p[jjj*nx+iii]+1;}
	      //		    }
	      //		  }
	      //		 }

	      //	      if (present_forward[k] != 0){

	      //	      	#pragma omp parallel for private(jjj,iii)
	      //			for (jjj=0; jjj<ny; ++jjj){
	      //			  for (iii=0; iii<nx; ++iii){
	      //	          if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsm[jjj*nx+iii]=rhsm[jjj*nx+iii]+log(arrinL0[jjj*nx+iii]); //add up frames (flatfielded with offpoint flatfield and limb darkening removed
		    //	            count_m[jjj*nx+iii]=count_m[jjj*nx+iii]+1;}
	      //   }
	      // 	}

	      //  }
      	      drms_free_array(arrin0);
	    }


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
	  //loop over filtergrams of the same type
	  //**********************************************************************************

	  printf("cosmic ray detection\n");


	  double avg=0.0;
	  printtime();

	  for (kkk=2; kkk<count_filtergram; ++kkk){
	  
	    if (cosmic_flag || ccount >= cthreshold)
	      {
	  
	   
	    k=index[kkk];
	    km1=index[kkk-1];

	    printf("fsn %d %d\n", k, keyvalue_fsn[k]);
	 
	    segin    = drms_segment_lookupnum(rec0[k], 0);
	    arrin0= drms_segment_read(segin, type_double, &status);
	    

	  
	    if (status != DRMS_SUCCESS || arrin0 == NULL)
	    {
	      printf("Error: there is a problem with the filtergram number %d \n", k); //if there is a problem with the filtergram
	      present_backward[k]=0;
	      present_forward[k]=0;
	      present[k]=0;
	    }
	  else
	    {
	     
	      printf("filtergram with FSN %d read\n",keyvalue_fsn[k]);
	      arrinL0 = arrin0->data;

	      query_flat =keyvalue_flatnumb[k];
	      if (strcmp(query_flat, query_flat_ref) != 0)
		{
		  printf("%d %d read new flatfield !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", k, keyvalue_fsn[k]);

		  status_flatfield=get_flatfields(query_flat, t_0, cam+1, focus,  app_flat, flat_off,  badpix, recnum, tobs_link, &rot_cur);

		  if (status_flatfield ==0)
		    {
		    printf("Flatfield read \n recnum offpoint %ld\n", recnum[0]);
		    }
		  else
		    {
		      printf("Failure reading flatfield / Use unity flatfield\n");
		    }
		 


		  query_flat_ref=query_flat;
		}

	      if (status_flatfield != 0)
		{
	      statarr[k]=statarr[k]+2;
	      present_forward[k]=0;
	      present_backward[k]=0;
		}
	      	    
	     


	      index_last=index[kkk-2];
	      index_current=index[kkk-1];
	      index_next=index[kkk];
	      count=-1;
	 
	      if (statarr[km1] < 4)
		count=0;
		{
	 #pragma omp parallel 
	      {    
#pragma omp for private(jjj,iii,wlr)	   
		for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
		{
		  arrimg[next][jjj][iii]=(float)arrinL0[jjj*nx+iii]*app_flat[jjj*nx+iii]/flat_relative[jjj*nx+iii];
		 if (rr[jjj*nx+iii] < rad_cosmic_ray)		 
		   {
		     wlr=((float)keyvalue_wl[km1]-505.)/2.0*lambda_sep-keyvalue_vrad[km1]/v_c*lambda0-(cos(keyvalue_p0[km1]/180.*M_PI)*((float)iii-keyvalue_X0[km1])-sin(keyvalue_p0[km1]/180.*M_PI)*((float)jjj-keyvalue_Y0[km1]))*cos(keyvalue_b0[km1]/180.0*M_PI)/keyvalue_rsun[km1]*cpa.rotcoef0*radsun_mm/v_c*lambda0;
		  limw[jjj*nx+iii]=a1[jjj*nx+iii]*exp(-(wlr-coef3[cam][0])*(wlr-coef3[cam][0])/2.0/w1[jjj*nx+iii]/w1[jjj*nx+iii])+a2[jjj*nx+iii]*exp(-(wlr-coef3[cam][1])*(wlr-coef3[cam][1])/2.0/w2[jjj*nx+iii]/w2[jjj*nx+iii])+cmc[jjj*nx+iii];
		   }
		}


	
		    
	  //head= NULL;
#pragma omp for private(iii,jjj)
	  for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
	    {
	      cmarr[jjj*nx+iii]=1;
	   
	      if (rr[jjj*nx+iii] < rad_cosmic_ray)
		{
		  if (fabs(arrimg[current][jjj][iii]-arrimg[last][jjj][iii]) > (factor[cam]*limw[jjj*nx+iii]) && fabs(arrimg[current][jjj][iii]-arrimg[next][jjj][iii]) > (factor[cam]*limw[jjj*nx+iii])){
	     
		    if ((tmind[index_current]-tmind[index_last]) < time_limit && (tmind[index_next]-tmind[index_current]) < time_limit){
		    if (badpix[jjj*nx+iii]){
		    if (!isnan(arrimg[last][jjj][iii]) && !isnan(arrimg[current][jjj][iii]) && !isnan(arrimg[next][jjj][iii])){
        #pragma omp critical(detect) 
		      {
		     cosmicarr[count]=jjj*nx+iii;
		     cosmicval[count]=minval(fabs(arrimg[current][jjj][iii]-arrimg[last][jjj][iii]), fabs(arrimg[current][jjj][iii]-arrimg[next][jjj][iii]));
		     cosmicsig[count]=cosmicval[count]/limw[jjj*nx+iii];

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

		printf("fsn count %d %d %d %d\n", keyvalue_fsn[km1], count, statarr[km1], status_flatfield);
	  cosmic_ray_check[index_current]=count;
	    
	




	  /////////////////////////////////////////////////////////
	  //write out cosmic ray series
	  /////////////////////////////////////////////////////////

	  if (cosmic_flag && present[km1])
	    {
	  int cosfind, limit_flag;
	  if (count == -1) cosfind=0; else if (count < limit_cosmic){cosfind=count; limit_flag=0;} else {cosfind=limit_cosmic; limit_flag=1;}


	    recout = dataout->records[kkk-2];
	    status=drms_setkey_int(recout, keyfsn, keyvalue_fsn[km1]);
	    status=drms_setkey_time(recout, keytobs, time_fl[km1]);
	    status=drms_setkey_int(recout, keycount, count);
	    status=drms_setkey_int(recout, fidkey, fid);
	    status=drms_setkey_int(recout, keycamera, cam+1);
	    status=drms_setkey_int(recout, keyexmax, limit_flag);
	    status=drms_setkey_float(recout, keylimit, factor[cam]);

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
	    
	    
	    
	  



	   
	  //**************************************************************************************************
	  //permanent bad pixel calculation
	  //*************************************************************************************************

	  
	      if (ccount >= cthreshold){
		update_stat=1;

		
#pragma omp parallel
	{
	  //if (present_backward[km1] != 0 && kkk > 2){
	  if (present_backward[km1] != 0){
#pragma omp for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0 && cmarr[jjj*nx+iii]){rhsp[jjj*nx+iii]=rhsp[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]);  // //add up frames for backward
		      count_p[jjj*nx+iii]=count_p[jjj*nx+iii]+1;}
		  }
		}
	  }
	
	  //if (present_forward[km1] != 0 && kkk < (count_filtergram-1)){
	  if (present_forward[km1] != 0){
		#pragma omp for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0 && cmarr[jjj*nx+iii]){rhsm[jjj*nx+iii]=rhsm[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]);   // //add up frames for forward
		      count_m[jjj*nx+iii]=count_m[jjj*nx+iii]+1;}
		  }
		}
	      }
	}

	      
      	      
	      }

	      drms_free_array(arrin0);
	    }
	      }


	    holder=next;
	    next=last;
	    last=current;
	    current=holder;


	  }//end loop over filtergrams of same type


	  km1=index[kkk-1];  //last filtergram 
	  #pragma omp parallel
	  	{
	  	  if (present_backward[km1] != 0){

	  #pragma omp for private(jjj,iii)
	  		for (jjj=0; jjj<ny; ++jjj){
	  		  for (iii=0; iii<nx; ++iii){
	  		    if (!isnan(arrimg[current][jjj][iii]) && arrimg[current][jjj][iii] > 0.0){rhsp[jjj*nx+iii]=rhsp[jjj*nx+iii]+log((double)arrimg[current][jjj][iii]);  // //add up last frame to backward (no cosmic ray check!)
	  	      count_p[jjj*nx+iii]=count_p[jjj*nx+iii]+1;}
	  	  }
	  	}
	    }
	  	}




	  printtime();	  
	  printf("loop done\n");
	  drms_close_records(dataout, DRMS_INSERT_RECORD);

	  printf("reference flatfield and limb darkening\n");
   
	  limb_darkening((double)R_SUN, (double)XX0, (double)YY0, b_coef, b_order, limb_dark);

	  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flatc[j*nx+i]=log(limb_dark[j*nx+i]);


 #pragma omp parallel for private(jjj,iii)
	  for (jjj=0; jjj<ny; ++jjj){
	    for (iii=0; iii<nx; ++iii){
	       
	      if (count_p[jjj*nx+iii] > 0) rhsp[jjj*nx+iii]=(rhsp[jjj*nx+iii]/(double)count_p[jjj*nx+iii]-flatc[jjj*nx+iii]); 
	      if (count_m[jjj*nx+iii] > 0) rhsm[jjj*nx+iii]=(rhsm[jjj*nx+iii]/(double)count_m[jjj*nx+iii]-flatc[jjj*nx+iii]); 
	    }
	  }


	  if (debug)  //
	{
         printf("write debug\n");
           FILE *fileptr;                                                                                                                          
           float *aaa;                                                                                                                            
           aaa=(float *)(malloc(nx*sizeof(float)));                                                                                              
           fileptr = fopen ("/tmp20/richard/interpol/dd.bin", "w");                                                                                    
	   for (j=0; j<ny; ++j){ for (i=0;i<nx;i++){aaa[i]=rhsm[j*nx+i];} fwrite ((char*)(aaa),sizeof(float),nx,fileptr);}
           fclose(fileptr);                                                                                                                            
      	free(aaa);
	}


      //*******************************************************************/
      //do flatfield calculations                                          /
      //*******************************************************************/

	  if (ccount >= cthreshold){

  
      for (j=0; j<ny; ++j)  for (i=0; i<nx; ++i) flati[j*nx+i]=0.0;
 
      rotpairs=ccount;

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
	 
	  }

	


      if (status == 0) printf("flatfield for camera %d calculated \n", cam); else printf("could not calcuate flatfield for camera %d \n", cam);

    
     
    //********************************
   

  // highpass flatfields
    
      //     highpass(nx, ny, fwhm, flat); //   
 


      //*********************************************************************
      //update flatfield
      //********************************************************************

 


   int axisout[2]={nx,ny};                             //size of the output arrays
      arrout_new  = drms_array_create(type_float,2,axisout,NULL,&status);
      printf("array size %ld\n", drms_array_size(arrout_new)/sizeof(float)); 
  
 

      flatfield_new=arrout_new->data;
  
 

      //update existing flatfield
       bad_pix=0;
 

         
       apod_circ(R_SUN*cpa.croprad*0.95, R_SUN*cpa.croprad*0.04, XX0-(float)nx/2.0, YY0-(float)ny/2.0, apod);

 
  

      //copy and apodize flatfield
for (j=0; j<ny; ++j){
      for (i=0; i<nx; ++i){
	flatfield_new[j*nx+i]=(flat[i][j]-1.0)*apod[j*nx+i]+1.0;
      }
 }



    rot_new.rotbad=bad_pix;
    rot_new.rotpairs=rotpairs;
    rot_new.rotcadence=(float)deltat;
    rot_new.flatfield_version=rot_cur.flatfield_version;


    //**************************************************************************************
    //write out flatfields (original with new T_STOP + updated
    //**************************************************************************************
  

    //    if (update_stat == 0 || update_stat == 1) write_flatfields(arr_flat, arrout_new, 2, recnum,
    //								     tobs_link, focus, rot_new, rot_cur);

    //    if (update_stat[0] == 0 || update_stat[0] == 1) write_flatfields(arr_flat_side, arrout_new_side, 1, recnum_side,
    //								   tobs_link_side, focus, rot_new_side, rot_cur_side);



    printf("update_stat %d\n", update_stat);

    printf("query flat %s\n", query_flat);
    printf("recnum %ld\n", recnum[0]);
    if (update_stat == 0) write_flatfield_fid(arrout_new, cam+1, recnum_ref[0], query_flat_relative, tobs_link, t_0, focus, fid, rot_new); 



	  drms_free_array(arrout_new);

	  free(param);
    }
  else 
    {
      printf("No data records found\n");
    }     

    
  
  if (data != NULL)
    {
      printf("close DRMS session \n");
      drms_close_records(data,DRMS_FREE_RECORD);       
    }

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
  //free(sigmacoef);
  free(apod);
	  free(flat_off);
	  free(flat_relative);
	  free(badpix);
	  free(app_flat);
	  free(flati);
	  free(flatc);
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







//module_flatfield_256 input_series="su_production.lev1_hmi_test[986331-987330]" cosmic_rays=1

//valgrind --leak-check=full --track-origins=yes --show-reachable=yes module_flatfield_256 input_series="su_production.lev1_hmi_test[986331-986430]" cosmic_rays=1
