


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
#define cadence_front_name "cadence_front"  //name of cadence for front camera
#define cadence_side_name "cadence_side"  //name of cadence for side camera
#define fid_name "fid" //name of fid argument


                                        //arguments of the module
ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn, "",  "Input data series."},
     {ARG_INT, kDSCosmic, 0, "Check for cosmic rays flag"},
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


int get_flatfields(int camera, TIME t_0, int focus,
                   float *flatfield, float *offpoint, 
                   short *bad, long long recnum[6], TIME tobs_link[2],
		   struct rotpar *rot_cur);


int write_flatfields(DRMS_Array_t *arr_flat, DRMS_Array_t *arrout_new, int camera,
                     long long recnum[6], TIME tobs_link[2], int focus,
                     struct rotpar rot_new, 
                     struct rotpar rot_cur);

int write_flatfield_fid(DRMS_Array_t *arrout_new, int camera,
			TIME tobs_link[2], int focus,int fid, int fsns,
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
  int cadence_front=cmdparams_get_int(&cmdparams, cadence_front_name, NULL);
  int cadence_side=cmdparams_get_int(&cmdparams, cadence_side_name, NULL);
  int fid=cmdparams_get_int(&cmdparams, fid_name, NULL);

  if (fid < minfid || fid > maxfid){printf("Not an observable FID\n"); exit(EXIT_FAILURE);}

  int  status= DRMS_SUCCESS; 
  int status0=DRMS_SUCCESS;
  int stat = DRMS_SUCCESS; 
 

  DRMS_RecordSet_t *data, *dataout;
  DRMS_Record_t  *ff_last_front=NULL, *ff_last_side=NULL;
  DRMS_Array_t *arrin0;
  double *arrinL0; //pointer to data
  DRMS_Segment_t *segin = NULL, *segout = NULL;
  DRMS_Record_t *recout;

  DRMS_Type_t type_time = DRMS_TYPE_TIME;
  DRMS_Type_t type_int = DRMS_TYPE_INT;
  DRMS_Type_t type_float = DRMS_TYPE_FLOAT;
  DRMS_Type_t type_double = DRMS_TYPE_DOUBLE;

  // int axisin[2]={nx,ny};                                                    //size of input arrays
  int axisout[2]={nx,ny};                             //size of the output arrays
  double deltat[2]={(double)cadence_side, (double)cadence_front};

  TIME t_0;
  int focus,camid,fsns;
  TIME   interntime;
  int axisbad[1];


  int i,j,k,km1,c,ki,iii,jjj,kkk,cam;                                                   //loop variables



 
 //********************************************************************************************************************************/
  //Constants to be set: will maybe be input parameters 
  //********************************************************************************************************************************/



  //sigma for each FID (cosmic rays)

 
  double *sigma=(double *)(malloc(nx*ny*sizeof(double)));
  float *limw=(float *)(malloc(nx*ny*sizeof(float)));


  //DRMS_arrays for output arrays
  DRMS_Array_t *arr_flat_front, *arr_flat_side;
  DRMS_Array_t *arrout_new_front, *arrout_new_side;
  DRMS_Array_t *arrout;

  int *cosmic_ray_data;

  //Array pointers for storage
  float *flatfield_front, *offpoint_front;
  float *flatfield_side, *offpoint_side;
  short *bad_front, *bad_side;
  short *badpix, *badpix_t;
  double *flati;
  double *flatc;

  offpoint_front=(float *)(malloc(nx*ny*sizeof(float)));
  offpoint_side=(float *)(malloc(nx*ny*sizeof(float)));
  flatfield_front=(float *)(malloc(nx*ny*sizeof(float)));
  flatfield_side=(float *)(malloc(nx*ny*sizeof(float)));
  bad_front=(short *)(malloc(nx*ny*sizeof(short)));
  bad_side=(short *)(malloc(nx*ny*sizeof(short)));
  badpix_t=(short *)(malloc(nx*ny*sizeof(short)));
  flati=(double *)(malloc(nx*ny*sizeof(double)));
  flatc=(double *)(malloc(nx*ny*sizeof(double)));

  float *flatfield_front_new, *flatfield_side_new;
  float *fflat_front, *fflat_side;
  
  float *flat_off, *flat_yest;

 

  double flat_front[nx][ny];
  double flat_side[nx][ny];
 


  double *rhsp, *rhsm;
  double *limb_dark;
  int *count_p, *count_m;

  int count, ccount, count_filtergram, bad_pix_front, bad_pix_side; // counters 
  int count_flatfields;

  int last, current, next, holder;
  int index_last, index_current, index_next;
  float rad, wlr, rr;
  float *a1,*a2, *w1,*w2,*cmc;

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
  //double *subarray;
  
  float arrimg[3][ny][nx];
  int *cosmicarr;
  cosmicarr=(int *)(malloc(nx*ny*sizeof(int)));

  //subarray=(double *)(malloc(nx*ny*sizeof(double)));
  //for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) subarray[jjj*nx+iii]=0.0;
  //

  
  //parameters to identify records
  TIME tobs_link_front[2], tobs_link_side[2];
  long long recnum_front[6], recnum_side[6];
  struct rotpar rot_cur_front, rot_cur_side;
  struct rotpar rot_new_front, rot_new_side;




  printf("START!\n");

  /***********************************************************************************************************/
  /*CHECK WHETHER THE FLATFIELD OUTPUT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
    
       drms_series_exists(drms_env, filename_flatfield, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",filename_flatfield);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Output series %s exists.\n",filename_flatfield);
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

 
      data     = drms_open_records(drms_env,(char *)inRecQuery,&stat);
      if (data == NULL || stat != 0 || data->n == 0){printf("can not open records\n"); exit(EXIT_FAILURE);}

  int nRecs=data->n;
  printf("number of records %d\n", nRecs);

  int nfr_front, nfr_side;


  DRMS_Record_t *rec0[nRecs];
 
  //array for keyvalues
  int keyvalue_cam[nRecs];
  int keyvalue_wl[nRecs];
  int keyvalue_pl[nRecs];
  DRMS_Type_Value_t  keyvalue_time;
  long tmind[nRecs+1];
  int keyvalue_fid[nRecs];
  char *keyvalue_iss[nRecs];
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



 double *sigmacoef;

 sigmacoef=(double *)(malloc(5*sizeof(double)));
 sigmacoef[0]=constsigma[(fid-10000)/10]*8.0; sigmacoef[1]=0.0; sigmacoef[2]=0.0; sigmacoef[3]=0.0; sigmacoef[4]=0.0;


  
  //list_pointer=(struct list **)(malloc(nRecs*sizeof(struct list *)));




  int update_stat[2]={2,2};   // 2: not enough data // 1: no update required // 0: flatfield updated 
    

  if (stat == DRMS_SUCCESS && data != NULL && nRecs > 0)
    {

      nfr_front=0;
      nfr_side=0;
      for (k=0; k<nRecs; ++k)
	{
	  cosmic_ray_check[k]=-1;
	  rec0[k]=data->records[k];
	  keyvalue_cam[k]=drms_getkey_int(rec0[k],keycam,&status);
	  if (keyvalue_cam[k] == cam_id_side) ++nfr_side;
	  if (keyvalue_cam[k] == cam_id_front) ++nfr_front;
	 
	}

      t_0 = drms_getkey_time(rec0[nRecs-1],keytobs,&status);
      fsns=drms_getkey_int(rec0[nRecs-1],keyfsn,&status);
      focus = drms_getkey_int(rec0[nRecs-1],keyfocus,&status);

         
     
      printf("get flatfields from database\n");

      get_flatfields(2, t_0, focus,  flatfield_front, 
		     offpoint_front,  bad_front, recnum_front, tobs_link_front, &rot_cur_front);  //get front camera flatfield


    //********************************

      printf("version ff %d %f\n", rot_cur_front.flatfield_version,rot_cur_front.rotcadence);
   
      
      get_flatfields(1, t_0, focus,  flatfield_side,
		     offpoint_side,  bad_side, recnum_side, tobs_link_side, &rot_cur_side);  //get side camera flatfield

      printf("version ff %d\n", rot_cur_side.flatfield_version);
     
      
      //********************************************************************************************************************************/
    
 
   
     
      
      for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){flat_front[i][j]=0.0; flat_side[i][j]=0.0;}

  
      

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
	statarr[k]=0;

	keyvalue_rsun[k]=drms_getkey_float(rec0[k], lev1_r_sun,&status); // /drms_getkey_float(rec0[k], lev1_imsc,&status0);
	if (keyvalue_rsun[k] < rsun_min || keyvalue_rsun[k] > rsun_max ||  isnan(keyvalue_rsun[k]))statarr[k]=statarr[k]+1; //!! check for status, status0

	keyvalue_dsun_obs[k]=drms_getkey_double(rec0[k], lev1_dist ,&status);
	if (status !=0 || keyvalue_dsun_obs[k] < dsun_obs_min || keyvalue_dsun_obs[k] > dsun_obs_max || isnan(keyvalue_dsun_obs[k])) statarr[k]=statarr[k]+2;

	keyvalue_p0[k]=drms_getkey_float(rec0[k], lev1_p0, &status); 
	if (status !=0 || keyvalue_p0[k] < p0_min || keyvalue_p0[k] > p0_max || isnan(keyvalue_p0[k])) statarr[k]=statarr[k]+4;

	keyvalue_b0[k]=drms_getkey_float(rec0[k], lev1_b0, &status);
	if (status !=0 || keyvalue_b0[k] < b0_min || keyvalue_b0[k] > b0_max ||  isnan(keyvalue_b0[k])) statarr[k]=statarr[k]+8;

	keyvalue_X0[k]=drms_getkey_float(rec0[k], lev1_x0, &status);
	if (status !=0  || keyvalue_X0[k] < X0_min || keyvalue_X0[k] > X0_max || isnan(keyvalue_X0[k])) statarr[k]=statarr[k]+16;

	keyvalue_Y0[k]=drms_getkey_float(rec0[k], lev1_y0,&status);
	if (status !=0 || keyvalue_Y0[k] < Y0_min || keyvalue_Y0[k] > Y0_max || isnan(keyvalue_Y0[k])) statarr[k]=statarr[k]+32;

	keyvalue_fid[k]=drms_getkey_int(rec0[k],fidkey,&status);
	if (status !=0 || keyvalue_fid[k] < minfid || keyvalue_fid[k] > maxfid || isnan(keyvalue_fid[k])) statarr[k]=statarr[k]+64;

	keyvalue_iss[k]=drms_getkey_string(rec0[k],isskey,&status);
	if (status !=0 || strcmp(keyvalue_iss[k], "CLOSED") != 0) statarr[k]=statarr[k]+128;

	keyvalue_vrad[k]=drms_getkey_float(rec0[k], lev1_vr, &status);
	if (keyvalue_vrad[k] < vrad_min || keyvalue_vrad[k] > vrad_max) statarr[k]=statarr[k]+256;


	//	printf("%d \t %d \t %s \t %d\t %d \t %d \t %ld   \t %d \t %f \t %f \t %f \t %f \t %f \t %f %f %d\n", k, keyvalue_fsn[k], keyvalue_iss[k], keyvalue_wl[k], keyvalue_pl[k], keyvalue_cam[k], tmind[k], keyvalue_fid[k], keyvalue_rsun[k], keyvalue_dsun_obs[k], keyvalue_p0[k], keyvalue_b0[k], keyvalue_X0[k], keyvalue_Y0[k], keyvalue_vrad[k], statarr[k]);

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
	B_ANG=B_ANG-keyvalue_b0[k]/180.0f*M_PI;
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
	// !! all parameters modified to account for image rotation
       *(param+0)=(double)R_SUN;  //R_SUN in pixel
       *(param+1)=(double)(4095.0-XX0);    //center x in pixel
       *(param+2)=(double)(4095.0-YY0);    //center y in pixel
       *(param+3)=(double)P_ANG+M_PI;  //P-angle in rad // !! rot img
       *(param+4)=(double)(-B_ANG);  //B-angle in rad
       *(param+5)=(double)dist;   //distance in AU
      } 
      else
	{
	  printf("Invalid lev1 keywords\n");
	  exit(EXIT_FAILURE);
	}
	
            printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", R_SUN, XX0, YY0, P_ANG, B_ANG, dist);
      //********************************

      int rotpairs[2]={0,0};
          
   

      //*********limb darkening calculation ****************************
      limb_darkening((double)R_SUN, (double)XX0, (double)YY0, b_coef, b_order, limb_dark);

 
  

      //**************************************************************
      //Start flatfield calculations
      //**************************************************************

     
      for (cam=0; cam<=1; ++cam){   //loop over cameras // 

	
	printf("camera \t %d\n", cam);


	if (cam == 0){flat_off=offpoint_side; flat_yest=flatfield_side;}
	if (cam == 1){flat_off=offpoint_front; flat_yest=flatfield_front;}

	if (cam ==0){badpix=bad_side;}
	if (cam ==1){badpix=bad_front;}

	if (cam ==0)camid=cam_id_side;
	if (cam ==1)camid=cam_id_front;

	count_flatfields=0;

	for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)badpix_t[jjj*nx+iii]=badpix[(4095-jjj)*nx+(4095-iii)];

	//**************************************************************************************
        //loop over different types of filtergrams                                            // 
	//**************************************************************************************



  


	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){rhsm[j*nx+i]=0.0; rhsp[j*nx+i]=0.0; count_p[j*nx+i]=0; count_m[j*nx+i]=0; flati[j*nx+i]=0.0;}

	
	for (i=0; i<nRecs; ++i){present_backward[i]=0; present_forward[i]=0;}
	for (i=0; i<nRecs; ++i){present[i]=0; index[i]=-1; statarr[k]=0;}

	  //count filtergrams, and number of valid pairs
	  count_filtergram=0; 
   
	  for (k=0; k<nRecs; ++k) if (statarr[k] == 0 && keyvalue_fid[k] == fid && keyvalue_cam[k] == camid){present[k]=1; index[count_filtergram]=k; ++count_filtergram;}  //loop over all filtergrams: search for valid filtergrams with fid, camera cam
	 
 

	//check each forward and backward pair of images for identity
	//check for wavelength, polarization,  P-angle, center, and solar radius identity, and for time difference being nominal cadence
	//backward pairs

	  if (count_filtergram >=3){  //if at least 3 filtergrams of this type

	  
	  printf("filtergram with FID  %d: %d\n", fid, count_filtergram);
	for (c=1; c<count_filtergram; ++c){


	  present_backward[index[c]]=present[index[c]];
	  if (keyvalue_wl[index[c]] != keyvalue_wl[index[c-1]] || keyvalue_pl[index[c]] != keyvalue_pl[index[c-1]] || (tmind[index[c]]-tmind[index[c-1]]) != (long)deltat[cam] || fabs(keyvalue_p0[index[c]]- keyvalue_p0[index[c-1]]) > 0.2 || sqrt(pow(keyvalue_X0[index[c]]-keyvalue_X0[index[c-1]],2)+pow(keyvalue_Y0[index[c]]-keyvalue_Y0[index[c-1]],2)) > limit_centerdiff || fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c-1]]) > limit_rsundiff){present_backward[index[c]]=0;}
	  //  printf("%d %d %ld %f %f\n", keyvalue_wl[index[c]]-keyvalue_wl[index[c-1]],keyvalue_pl[index[c]]- keyvalue_pl[index[c-1]], (tmind[index[c]]-tmind[index[c-1]]), fabs(keyvalue_p0[index[c]]- keyvalue_p0[index[c-1]]), fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c-1]]));
	}


	//forward pairs
	ccount=count_filtergram-1; //number of pairs (forward or backward)
	for (c=0; c<count_filtergram-1; ++c){
	  present_forward[index[c]]=present[index[c]];
	  if (keyvalue_wl[index[c]] != keyvalue_wl[index[c+1]] || keyvalue_pl[index[c]] != keyvalue_pl[index[c+1]] || (tmind[index[c+1]]-tmind[index[c]]) != (long)deltat[cam] || fabs(keyvalue_p0[index[c+1]]- keyvalue_p0[index[c]]) > 0.2 || sqrt(pow(keyvalue_X0[index[c]]-keyvalue_X0[index[c+1]],2)+pow(keyvalue_Y0[index[c]]-keyvalue_Y0[index[c+1]],2)) > limit_centerdiff || fabs(keyvalue_rsun[index[c]]-keyvalue_rsun[index[c+1]]) > limit_rsundiff){present_forward[index[c]]=0; --ccount;}
	}

	printf("number of filtergrams, valid pairs of frames %d \t %d \n", count_filtergram, ccount);

	
	for (kkk=0; kkk<count_filtergram; ++kkk){k=index[kkk]; printf("%d %d %d %d %ld %d %d %d\n",kkk,k,keyvalue_fid[k],keyvalue_cam[k],tmind[k]-tmind[0],present[k],present_backward[k],present_forward[k]);}
	
      //combine flatfield correction and limb-darkening

	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flatc[j*nx+i]=log(limb_dark[j*nx+i])+log(flat_off[j*nx+i])-
	  log((double)flat_yest[j*nx+i]);


	  //do cosmic ray detection for each filtergram


	limb_darkening((double)R_SUN, (double)XX0, (double)YY0, sigmacoef, 5, sigma);




    
	 
	  //read in first filtergram
	  

      

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
	      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) arrimg[last][jjj][iii]=NAN; // !! find proper solution
	    }
	  else
	    {
	     
	      //   printf("filtergram with FSN %d read\n",keyvalue_fsn[index_last]);
	      arrinL0 = arrin0->data;
	      //       for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) if (isnan(arrinL0[jjj*nx+iii])) arrinL0[jjj*nx+iii]=0.0;

	      if (cosmic_flag)
		for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) arrimg[last][jjj][iii]=(float)arrinL0[jjj*nx+iii]; 

	      k=index[0];

	    
	      
	
	      if (present_forward[k] != 0){

	      #pragma omp parallel for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsm[(4095-jjj)*nx+(4095-iii)]=rhsm[(4095-jjj)*nx+(4095-iii)]+log(arrinL0[jjj*nx+iii]); //!! rotate _image   //add up frames (flatfielded with offpoint flatfield and limb darkening removed
		      count_m[(4095-jjj)*nx+(4095-iii)]=count_m[(4095-jjj)*nx+(4095-iii)]+1;}
		  }
		}
	
	      }
	      drms_free_array(arrin0);
	    }

	  
	
//read in second filtergram

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
	      
	     
	      //printf("filtergram with FSN %d read\n",keyvalue_fsn[index_current]);
	      arrinL0 = arrin0->data;
	      //          for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) if (isnan(arrinL0[jjj*nx+iii])) arrinL0[jjj*nx+iii]=0.0;

	      if (cosmic_flag)
		for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) arrimg[current][jjj][iii]=(float)arrinL0[jjj*nx+iii];

	      k=index[1];
	      printf("fsn %d %d\n", k, keyvalue_fsn[k]);

	      if (present_backward[k] != 0){

		#pragma omp parallel for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsp[(4095-jjj)*nx+(4095-iii)]=rhsp[(4095-jjj)*nx+(4095-iii)]+log(arrinL0[jjj*nx+iii]); // !! rot img //add up frames (flatfielded with offpoint flatfield and limb darkening removed
		      count_p[(4095-jjj)*nx+(4095-iii)]=count_p[(4095-jjj)*nx+(4095-iii)]+1;}
		  }
		}
	      }

	      if (present_forward[k] != 0){
	printf("forward %d\n", k);
		#pragma omp parallel for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsm[(4095-jjj)*nx+(4095-iii)]=rhsm[(4095-jjj)*nx+(4095-iii)]+log(arrinL0[jjj*nx+iii]); // !! rot img   //add up frames (flatfielded with offpoint flatfield and limb darkening removed
		      count_m[(4095-jjj)*nx+(4095-iii)]=count_m[(4095-jjj)*nx+(4095-iii)]+1;}
		  }
		}

	      }
      	      drms_free_array(arrin0);
	    }


#pragma omp parallel for private(jjj,iii,rr)	   
	      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
		{
				  
		  rr=sqrt(pow((float)iii-XX0,2)+pow((float)jjj-YY0,2))/R_SUN;
		  a1[jjj*nx+iii]=coef0[0]+coef0[1]*rr+coef0[2]*rr*rr;
		  a2[jjj*nx+iii]=coef1[0]+coef1[1]*rr+coef1[2]*rr*rr;
		  w1[jjj*nx+iii]=coef2[0]+coef2[1]*rr;
		  w2[jjj*nx+iii]=coef3[0]+coef3[1]*rr;
		  cmc[jjj*nx+iii]=coef4[0]+coef4[1]*rr+coef4[2]*rr*rr;

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
	    //printf("fid value %d\n", keyvalue_fid[k]);
	    //printf("order %d %d %d\n", last, current, next);

	    segin    = drms_segment_lookupnum(rec0[k], 0);
	    arrin0= drms_segment_read(segin, type_double, &status);
	    

	  
	  if (status != DRMS_SUCCESS)
	    {
	      printf("Error: there is a problem with the filtergram number %d \n", k); //if there is a problem with the filtergram
	      present_backward[k]=0;
	      present_forward[k]=0;
	      present[k]=0;
	    }
	  else
	    {
	     
	      //   printf("filtergram with FSN %d read\n",keyvalue_fsn[k]);
	      arrinL0 = arrin0->data;

	      //#pragma omp parallel for private(jjj,iii)
		//            for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) if (isnan(arrinL0[jjj*nx+iii])) arrinL0[jjj*nx+iii]=0.0;
	               
	    

	  if (cosmic_flag)
	    {

	      index_last=index[kkk-2];
	      index_current=index[kkk-1];
	      index_next=index[kkk];
	      count=0;


	 #pragma omp parallel 
	      {    
#pragma omp for private(jjj,iii,wlr)	   
	      for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
		{
		  arrimg[next][jjj][iii]=(float)arrinL0[jjj*nx+iii];
				 
		  wlr=((float)keyvalue_wl[km1]-505.)/2.0*69e-3-keyvalue_vrad[km1]/3e8*6173.4-(cos(keyvalue_p0[km1]/180.*M_PI)*((float)iii-keyvalue_X0[km1])-sin(keyvalue_p0[km1]/180.*M_PI)*((float)jjj-keyvalue_Y0[km1]))/keyvalue_rsun[km1]*cpa.rotcoef0*698.0/3e8*6173.4;
		  limw[jjj*nx+iii]=a1[jjj*nx+iii]*exp(-pow(wlr-0.162,2)/2.0/w1[jjj*nx+iii]/w1[jjj*nx+iii])+a2[jjj*nx+iii]*exp(-pow(wlr-0.275,2)/2.0/w2[jjj*nx+iii]/w2[jjj*nx+iii])+coef4[0]+cmc[jjj*nx+iii];
		}

	   	     
	    
	  //head= NULL;
#pragma omp for private(iii,jjj,rad)
	  for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii)
	    {
	      rad=sqrt(pow((float)iii-XX0,2)+pow((float)jjj-YY0,2))/R_SUN;
	      if (rad < rad_cosmic_ray)
		{
		  if (fabs(arrimg[current][jjj][iii]-arrimg[last][jjj][iii]) > (factor*limw[jjj*nx+iii]) && fabs(arrimg[current][jjj][iii]-arrimg[next][jjj][iii]) > (factor*limw[jjj*nx+iii])){
	     
		if ((tmind[index_current]-tmind[index_last]) < time_limit && (tmind[index_next]-tmind[index_current]) < time_limit){
		  if (badpix[jjj*nx+iii]){
		    if (!isnan(arrimg[last][jjj][iii]) && !isnan(arrimg[current][jjj][iii]) && !isnan(arrimg[next][jjj][iii])){
                     #pragma omp critical(detect) 
		      {
		    cosmicarr[count]=jjj*nx+iii;
		    ++count; //count number of cosmic rays (using first differences)
		      }
		    }
		  }
		}
	      }
		}
	    }
	      }
	  printf("fsn count %d %d\n", keyvalue_fsn[km1], count);
	  cosmic_ray_check[index_current]=count;
	    
	      
	  

	  /////////////////////////////////////////////////////////
	  //write out cosmic ray series
	  /////////////////////////////////////////////////////////

	  dataout = drms_create_records(drms_env,1,filename_cosmic,DRMS_PERMANENT,&stat);
	    recout = dataout->records[0];
	    status=drms_setkey_int(recout, keyfsn, keyvalue_fsn[km1]);
	    status=drms_setkey_time(recout, keytobs, time_fl[km1]);
	    status=drms_setkey_int(recout, keycount, count);
	    status=drms_setkey_int(recout, fidkey, fid);
	    status=drms_setkey_int(recout, keycamera, cam+1);


	    drms_keyword_setdate(recout);
	    
	    if (keyvalue_cam[km1] == cam_id_side) status=drms_setkey_string(recout, keyinstrument, camera_str_front);
            if (keyvalue_cam[km1] == cam_id_front) status=drms_setkey_string(recout, keyinstrument, camera_str_side);
	    
	    segout = drms_segment_lookup(recout, segmentname_cosmic);
	    
	    axisbad[0]=count;
            arrout=drms_array_create(type_int,1,axisbad,NULL,&status);
	    cosmic_ray_data=arrout->data;
	  	 
	    //curr=list_pointer[k];
	    //while(curr){
	    //  cosmic_ray_data[ncount]=curr->val;
	    //  curr=curr->next;
	    //  ++ncount;
	    //}

	    for (i=0; i<count; ++i) cosmic_ray_data[i]=cosmicarr[i];
	    status=drms_segment_write(segout, arrout, 0);
	    drms_close_records(dataout, DRMS_INSERT_RECORD);
	    drms_free_array(arrout);
	    

	    printf("cosmic ray detection done\n");
	    }
	    
	  



	  holder=next;
	  next=last;
	  last=current;
	  current=holder;

	   
	  //**************************************************************************************************
	  //permanent bad pixel calculation
	  //*************************************************************************************************

	  
	      if (ccount >= cthreshold){
		update_stat[cam]=1;

		//	printf("begin loop\n");
	
#pragma omp parallel
	{
	  if (present_backward[k] != 0){
	    printf("backward %d\n", k);

#pragma omp for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsp[(4095-jjj)*nx+(4095-iii)]=rhsp[(4095-jjj)*nx+(4095-iii)]+log(arrinL0[jjj*nx+iii]);  // !! rot image //add up frames (flatfielded with offpoint flatfield and limb darkening removed
		      count_p[(4095-jjj)*nx+(4095-iii)]=count_p[(4095-jjj)*nx+(4095-iii)]+1;}
		  }
		}
	  }
	
	      if (present_forward[k] != 0){
	printf("forward %d\n", k);
		#pragma omp for private(jjj,iii)
		for (jjj=0; jjj<ny; ++jjj){
		  for (iii=0; iii<nx; ++iii){
		    if (!isnan(arrinL0[jjj*nx+iii]) && arrinL0[jjj*nx+iii] > 0.0){rhsm[(4095-jjj)*nx+(4095-iii)]=rhsm[(4095-jjj)*nx+(4095-iii)]+log(arrinL0[jjj*nx+iii]);   //!! rot img  //add up frames (flatfielded with offpoint flatfield and limb darkening removed
		      count_m[(4095-jjj)*nx+(4095-iii)]=count_m[(4095-jjj)*nx+(4095-iii)]+1;}
		  }
		}
	      }
	}

	      
      	      
	      }

	      drms_free_array(arrin0);
	    }
	      }
	  }//end loop over filtergrams of same type

	  printtime();	  
	  printf("loop done\n");

 #pragma omp parallel for private(jjj,iii)
	  for (jjj=0; jjj<ny; ++jjj){
	    for (iii=0; iii<nx; ++iii){
	       
	      if (count_p[jjj*nx+iii] > 0) rhsp[jjj*nx+iii]=(rhsp[jjj*nx+iii]/(double)count_p[jjj*nx+iii]-flatc[(4095-jjj)*nx+(4095-iii)]); 
	      if (count_m[jjj*nx+iii] > 0) rhsm[jjj*nx+iii]=(rhsm[jjj*nx+iii]/(double)count_m[jjj*nx+iii]-flatc[(4095-jjj)*nx+(4095-iii)]); 
	    }
	  }

	 
      //*******************************************************************/
      //do flatfield calculations                                          /
      //*******************************************************************/

	  if (ccount >= cthreshold){

  
      for (j=0; j<ny; ++j)  for (i=0; i<nx; ++i) flati[j*nx+i]=0.0;
 
      rotpairs[cam]=rotpairs[cam]+ccount;
      cpa.norm=normconst;

      status=flatfield(rhsp, rhsm, badpix_t, ccount, flati, param, cpa, deltat[cam]);
  

      if (status ==0)
	{
      //*********************************************************************//
      //linear combination of  different flatfields       (here, average!!)  //
	  if (cam == 0){for (j=0; j<ny; ++j)  for (i=0; i<nx; ++i)flat_side[i][j]=flat_side[i][j]+exp(flati[j*nx+i]);} //average flatfields 
	  if (cam == 1){for (j=0; j<ny; ++j)  for (i=0; i<nx; ++i)flat_front[i][j]=exp(flati[j*nx+i])+flat_front[i][j];   printf("just added front flatfield\n");}
   
      ++count_flatfields;
      //*********************************************************************//
	}
	  }
	 
	  }

	
	printf("count flatfields %d \n", count_flatfields);
  

	if (cam == 0){
	  if (count_flatfields != 0){ for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)  flat_side[i][j]=flat_side[i][j]/(double)count_flatfields;} 
	  else {for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flat_side[i][j]=1.0;} // just average all available flatfields -> improve to weighted average
	}

	if (cam == 1){
	  if (count_flatfields != 0){ for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)  flat_front[i][j]=flat_front[i][j]/(double)count_flatfields;} 
	  else {for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flat_front[i][j]=1.0;}
	}

      if (count_flatfields != 0) printf("flatfield for camera %d calculated \n", cam); else printf("not enough data for flatfield for camera %d \n", cam);




   


      
      } //end camera loop
    
     
      if (debug)  //
	{
         printf("write debug\n");
           FILE *fileptr;                                                                                                                          
           double *aaa;                                                                                                                            
           aaa=(double *)(malloc(nx*sizeof(double)));                                                                                              
           fileptr = fopen ("/tmp20/richard/interpol/dd.bin", "w");                                                                                    
           for (j=0; j<ny; ++j){ for (i=0;i<nx;i++){aaa[i]=flat_side[i][j];} fwrite ((char*)(aaa),sizeof(double),nx,fileptr);}   
           for (j=0; j<ny; ++j){ for (i=0;i<nx;i++){aaa[i]=flat_front[i][j];} fwrite ((char*)(aaa),sizeof(double),nx,fileptr);}   
           fclose(fileptr);                                                                                                                            
      	free(aaa);
	}
    //********************************
   

  // highpass flatfields
    
      //     highpass(nx, ny, fwhm, flat_front); //   !!remove high pass for test
      //     highpass(nx, ny, fwhm, flat_side);  // 
    
  


      //*********************************************************************
      //update flatfield
      //********************************************************************

 
      arr_flat_front=drms_array_create(type_float,2,axisout,NULL,&status);
      arr_flat_side=drms_array_create(type_float,2,axisout,NULL,&status);

 
      arrout_new_front  = drms_array_create(type_float,2,axisout,NULL,&status);
      arrout_new_side  = drms_array_create(type_float,2,axisout,NULL,&status);

    
      fflat_front=arr_flat_front->data;
      fflat_side=arr_flat_side->data;

      flatfield_front_new=arrout_new_front->data;
      flatfield_side_new=arrout_new_side->data;


    for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) fflat_front[j*nx+i]=flatfield_front[j*nx+i];
    for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)   fflat_side[j*nx+i]=flatfield_side[j*nx+i];

    bad_pix_front=0;
    bad_pix_side=0;

    for (j=0; j<ny; ++j){
      for (i=0; i<nx; ++i){

      if (flat_front[i][j] < threshold_lower|| flat_front[i][j] > threshold_upper)  
       {++bad_pix_front; 
	 flatfield_front_new[j*nx+i]=flatfield_front[j*nx+i]*(1.0f+(float)flat_front[i][j]); update_stat[1]=0;
       }
     else
       {
	 flatfield_front_new[j*nx+i]=flatfield_front[j*nx+i];
       }

     if (flat_side[i][j] < threshold_lower || flat_side[i][j] > threshold_upper)   
       {
	 ++bad_pix_side; 
	 flatfield_side_new[j*nx+i]=flatfield_side[j*nx+i]*(1.0f+(float)flat_side[i][j]); update_stat[0]=0;
       }
     else 
       {
	 flatfield_side_new[j*nx+i]=flatfield_side[j*nx+i];
       }

      }
    }



    rot_new_front.rotbad=bad_pix_front;
    rot_new_front.rotpairs=rotpairs[1];
    rot_new_front.rotcadence=(float)deltat[1];
    rot_new_front.flatfield_version=rot_cur_front.flatfield_version;

    rot_new_side.rotbad=bad_pix_side;
    rot_new_side.rotpairs=rotpairs[0];
    rot_new_side.rotcadence=(float)deltat[0];
    rot_new_side.flatfield_version=rot_cur_side.flatfield_version;

   

    //**************************************************************************************
    //write out flatfields (original with new T_STOP + updated
    //**************************************************************************************
    printf("update_stat %d %d\n", update_stat[0], update_stat[1]);

    //    if (update_stat[1] == 0 || update_stat[1] == 1) write_flatfields(arr_flat_front, arrout_new_front, 2, recnum_front,
    //								     tobs_link_front, focus, rot_new_front, rot_cur_front);

  // if (update_stat[0] == 0 || update_stat[0] == 1) write_flatfields(arr_flat_side, arrout_new_side, 1, recnum_side,
    //								   tobs_link_side, focus, rot_new_side, rot_cur_side);


    //write out full flatfield for test //!!
for (j=0; j<ny; ++j){
      for (i=0; i<nx; ++i){
	flatfield_front_new[j*nx+i]=flat_front[4095-i][4095-j];
	flatfield_side_new[j*nx+i]=flat_side[4095-i][4095-j];
      }
 }




 if (update_stat[1] == 0 || update_stat[1] == 1) write_flatfield_fid(arrout_new_front, 2, tobs_link_front, focus, fid, fsns, rot_new_front); //!! writeout unconditional

 if (update_stat[0] == 0 || update_stat[0] == 1) write_flatfield_fid(arrout_new_side, 1, tobs_link_side, focus, fid, fsns, rot_new_side);   // !!writeout unconditional

    //*************************************************************************************
    //write out records of cosmic ray hits
    //*************************************************************************************

  //  if (cosmic_flag)
  //  {


  //printf("Write cosmic ray hit series\n");

  //  for (k=0; k<nRecs; ++k){
	  	    
  //    if (cosmic_ray_check[k] != -1){
  //     count=cosmic_ray_check[k];
	   
	    
	    //write out cosmic ray list
  //    dataout = drms_create_records(drms_env,1,filename_cosmicRMS_PERMANENT,&stat);
  //    recout = dataout->records[0];
  //    status=drms_setkey_int(recout, keyfsn, keyvalue_fsn[k]);
  //    status=drms_setkey_time(recout, keytobs, time_fl[k]);
  //    status=drms_setkey_int(recout, keycount, count);

  //    drms_keyword_setdate(recout);
	    
  //    if (keyvalue_cam[k] == cam_id_side) status=drms_setkey_string(recout, keyinstrument, camera_str_front);
  //        if (keyvalue_cam[k] == cam_id_front) status=drms_setkey_string(recout, keyinstrument, camera_str_side);
	    
  //    segout = drms_segment_lookup(recout, segmentname_cosmic);
	    
  //        int axisbad[1];
  //axisbad[0]=count;
  //        arrout=drms_array_create(type_int,1,axisbad,NULL,&status);
  //    cosmic_ray_data=arrout->data;
  //    int ncount=0;
	 
	    //curr=list_pointer[k];
	    //while(curr){
	    //  cosmic_ray_data[ncount]=curr->val;
	    //  curr=curr->next;
	    //  ++ncount;
	    //}
  //    for (jjj=0; jjj<ny; ++jjj) for (iii=0; iii<nx; ++iii) if (cosmicarr[iii][jjj] == 1) {cosmic_ray_data[ncount]=jjj*nx+iii; ++ncount;}

  //    status=drms_segment_write(segout, arrout, 0);
  //    drms_close_records(dataout, DRMS_INSERT_RECORD);
  //    drms_free_array(arrout);
  //    drms_free_records(dataout);
  //    }
  //	  }
  //
  //    }
 

    //******************************************************************

	  drms_free_array(arr_flat_front);
	  drms_free_array(arr_flat_side);
	  drms_free_array(arrout_new_front);
	  drms_free_array(arrout_new_side);
	  


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

  free(a1);
  free(a2);
  free(w1);
  free(w2);
  free(cmc);
  free(cosmicarr);
  free(sigmacoef);
	  free(offpoint_front);
	  free(offpoint_side);
	  free(bad_front);
	  free(bad_side);
	  free(flati);
	  free(flatc);
	  free(flatfield_front);
	  free(flatfield_side);
	  free(rhsp);
	  free(rhsm);
	  free(count_p);
	  free(count_m);
	  free(limb_dark);
	  free(sigma);
	  free(limw);


  printtime();  
  printf("COMPLETED!\n");

  
  
  return 0;

}







//module_flatfield_256 input_series="su_production.lev1_hmi_test[986331-987330]" cosmic_rays=1

//valgrind --leak-check=full --track-origins=yes --show-reachable=yes module_flatfield_256 input_series="su_production.lev1_hmi_test[986331-986430]" cosmic_rays=1
