
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


ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn, "",  "Input data series."},
     {ARG_STRING, datumn, "", "Datum string"},
     {ARG_END}
};

//void highpass(int, int, double, double[nx][ny]);
int get_flatfields(char *query_flat, TIME t_0, int camera, int focus,
                   float *flatfield, float *offpoint, 
                   short *bad, long long recnum[6], TIME tobs_link[2],
		   struct rotpar *rot_cur);


int write_flatfields(char *filename_flatfield_out, DRMS_Array_t *arr_flat, DRMS_Array_t *arrout_new, int camera,
                     long long recnum[6], TIME tobs_link[2], TIME t_0, int focus,
                     struct rotpar rot_new, 
                     struct rotpar rot_cur);

int DoIt(void)
{

#include "module_flatfield_const.h"

  const char *inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL); //cmdparams is defined in jsoc_main.h
  const char *datum =  cmdparams_get_str(&cmdparams, datumn, NULL);
  int cameraint = cmdparams_get_int(&cmdparams, cameran, NULL);

  int  status= DRMS_SUCCESS, stat=DRMS_SUCCESS;

  DRMS_Segment_t *segin = NULL;


  DRMS_Type_t type_time = DRMS_TYPE_TIME;
  DRMS_Type_t type_int = DRMS_TYPE_INT;
  DRMS_Type_t type_float = DRMS_TYPE_FLOAT;
  DRMS_Type_t type_double = DRMS_TYPE_DOUBLE;

  DRMS_RecordSet_t *data, *dataout;
  DRMS_Record_t  *ff_last_front=NULL, *ff_last_side=NULL;
  DRMS_Array_t *arr_flat, *arr_flat_new, *arr_flat_rel;

  int axisout[2]={nx,ny};

  TIME t_0;
  int focus, camera;
  float cadence;
  int npairstot=0;

  int i,j,k;
  float h1,h2;


  float *flatfield, *offpoint, *flat, *ffhp;
  short *bad;


  float *flatfield_new;
  float *ff, *ff_new, *ff_rel;

  long long recnum[6];
  TIME tobs_link[2];

  struct rotpar rot_cur, rot_new, rot_rel;
  struct fresize_struct fresizes;

 

  flatfield=(float *)(malloc(nx*ny*sizeof(float)));
  flatfield_new=(float *)(malloc(nx*ny*sizeof(float)));
  offpoint=(float *)(malloc(nx*ny*sizeof(float)));
  bad=(short *)(malloc(nx*ny*sizeof(short)));
  
  flat=(float *)(calloc(nx*ny, sizeof(float)));
  ffhp=(float *)(malloc(nx*ny*sizeof(float)));


  //set parallelization parameters
      int nthreads=1; 
      nthreads=omp_get_num_procs();                                      //number of threads supported by the machine where the code is running
      omp_set_num_threads(nthreads);                                     //set the number of threads to the maximum value
      printf("Number of threads run in parallel = %d \n",nthreads);



  /***********************************************************************************************************/
  /*CHECK WHETHER THE FLATFIELD OUTPUT SERIES EXIST                                                                    */
  /***********************************************************************************************************/
    
       drms_series_exists(drms_env, filename_flatfield_out, &status);
      if (status == DRMS_ERROR_UNKNOWNSERIES)
	{
	  printf("Output series %s doesn't exist\n",filename_flatfield_out);       //if the output series does not exit
	  exit(EXIT_FAILURE);                                        //we exit the program
	} 
      if (status == DRMS_SUCCESS)
	{
	  printf("Output series %s exists.\n",filename_flatfield_out);
	}


      char query[256]="";
      strcat(query, inRecQuery);
      strcat(query, "[");
      char cami[1]={""};
      sprintf(cami, "%d", cameraint);
      strcat(query, cami);
      strcat(query, "][");
      strcat(query, datum);
      strcat(query, "_00:00:00_TAI-");
      strcat(query, datum);
      strcat(query, "_23:59:59_TAI");
      strcat(query, "]");

      printf("%s\n", query);
      data     = drms_open_records(drms_env,query,&stat);

 if (data == NULL){printf("can not open records\n"); exit(EXIT_FAILURE);}

  int nRecs=data->n;
  printf("number of records %d\n", nRecs);



if (stat == DRMS_SUCCESS && data != NULL && nRecs > 0)
  {
  int keyvalue_fid[nRecs];
  DRMS_Record_t *rec0[nRecs];
  float *flatfield_fid[nRecs];
  int keyvalue_recnum[nRecs];
  char *keyvalue_query[nRecs];
  DRMS_Array_t *arrin0[nRecs];
  short wave[nRecs];
  short pol[nRecs];

  int npairs[nRecs];

  float nwav[6];
  for (i=0; i<6; ++i) nwav[i]=0.0;


for (k=0; k<nRecs; ++k)
  {
    rec0[k]=data->records[k];
    keyvalue_fid[k] = drms_getkey_int(rec0[k],fidkey,&status);
    keyvalue_recnum[k]=drms_getkey_int(rec0[k], recnumkey,&status);
    keyvalue_query[k]=drms_getkey_string(rec0[k], querykey,&status);
    segin    = drms_segment_lookupnum(rec0[k], 0);
    arrin0[k]= drms_segment_read(segin, type_float, &status);
    flatfield_fid[k] = arrin0[k]->data;
    wave[k]=((keyvalue_fid[k]-10000)/10-5)/2;
    pol[k]=(keyvalue_fid[k] % 10);
 
    npairs[k]= drms_getkey_int(rec0[k],keynpairs,&status);
    npairstot=npairstot+npairs[k];
    printf("%d %d %d %d %d\n", k,  keyvalue_fid[k], pol[k], wave[k], keyvalue_recnum[k]);
    if (wave[k] < 0 || wave[k] > 5){printf("not a tabulated wavelength\n"); exit(EXIT_FAILURE);}

    nwav[wave[k]]=nwav[wave[k]]+1.0;
  }

 

 char timefirst[256]="";
 strcat(timefirst, datum);
 strcat(timefirst, "_00:00:00.00_TAI");

 TIME t_0=sscan_time(timefirst)+24.0*60.0*60.0; 

 //t_0=drms_getkey_time(rec0[nRecs-1],keytstart,&status);
 focus = drms_getkey_int(rec0[nRecs-1],keyfocusflat,&status);
 camera= drms_getkey_int(rec0[nRecs-1],keycamera,&status);
 cadence=drms_getkey_float(rec0[nRecs-1],keycadence,&status);


 printf("%lf %d %d %f\n", t_0, focus, camera, cadence);

 int fvers; 
 float sum=0.0;


char *qstr_refflat=keyvalue_query[0];
 if (cameraint ==1) qstr_refflat="su_richard.flatfield[:#168]"; // !! debug
 if (cameraint ==2) qstr_refflat="su_richard.flatfield[:#169]"; // !! debug

for (k=0; k<nRecs; ++k)
  {
    ff=flatfield_fid[k];
    
   

    if (nwav[wave[k]] > 0.0 && keyvalue_recnum[k] == keyvalue_recnum[0])
      {
	sum=sum+cof[cameraint-1][wave[k]]/nwav[wave[k]];
	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flat[j*nx+i]=flat[j*nx+i]+cof[cameraint-1][wave[k]]*(float)ff[j*nx+i]/nwav[wave[k]];
      }
	  
  }

 printf("query flat %s\n", qstr_refflat);
 get_flatfields(qstr_refflat, t_0, cameraint, focus,  flatfield, 
		offpoint,  bad, recnum, tobs_link, &rot_cur);        

   printf("recnum %ld\n", recnum[0]);



   //for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) if (bad[j*nx+i]) ffield[j][i]=(double)flat[j*nx+i]; else ffield[j][i]=1.0;
   for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) if (bad[j*nx+i] == 0) flat[j*nx+i]=1.0;

 //highpass

const float kernel_size=4.0;
 float stdd=fwhm[cameraint-1]/2.0/sqrt(2.0*log(2.0)); 
 int nwd=(int)(fwhm[cameraint-1]*kernel_size);  //kernel size
 init_fresize_gaussian(&fresizes,stdd,nwd,1);

fresize(&fresizes,flat, ffhp, nx,ny,nx,nx,ny,nx,0,0,1.0);   

 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) ffhp[j*nx+i]=flat[j*nx+i]-ffhp[j*nx+i];

 // highpass(nx, ny, fwhm[cameraint-1], ffield);


 printf("sum %f\n", sum);

 if (sum != 1.0)
 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) ffhp[j*nx+i]=ffhp[j*nx+i]/sum;

 drms_close_records(data,DRMS_FREE_RECORD);


 

 int count=0;

 if (update_flag == 2)
   {
     for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
       {
	 if (bad[j*nx+i])
	   {
	     flatfield_new[j*nx+i]=(1.0+ffhp[j*nx+i])*offpoint[j*nx+i];
	     if (ffhp[j*nx+i] != 0.0) ++count;
	   }
	 else 
	   {
	     flatfield_new[j*nx+i]=offpoint[j*nx+i];
	   }
       }
   }

 if (update_flag == 1)
   {
     for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
       if (bad[j*nx+i] && (flat[j*nx+i] < threshold_lower|| flat[j*nx+i] > threshold_upper))
	 {
	   ++count;
	   flatfield_new[j*nx+i]=(1.0f+ffhp[j*nx+i])*offpoint[j*nx+i];
	 }
       else
	 {
	   flatfield_new[j*nx+i]=offpoint[j*nx+i];
	 }
   }

 if (update_flag == 0)
   printf("update flag 0\n");
   {
     for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
       flatfield_new[j*nx+i]=offpoint[j*nx+i];
   }




 arr_flat=drms_array_create(type_float,2,axisout,NULL,&status);
 arr_flat_new=drms_array_create(type_float,2,axisout,NULL,&status);
 arr_flat_rel=drms_array_create(type_float,2,axisout,NULL,&status);

 ff_new=arr_flat_new->data;
 ff=arr_flat->data;
 ff_rel=arr_flat_rel->data;

 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)ff[j*nx+i]=flatfield[j*nx+i];
 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)ff_new[j*nx+i]=flatfield_new[j*nx+i];
 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) if (bad[j*nx+i]) ff_rel[j*nx+i]=(1.0+ffhp[j*nx+i])*offpoint[j*nx+i]; else ff_rel[j*nx+i]=offpoint[j*nx+i];

rot_new.rotbad=count;
rot_new.rotpairs=npairstot;
rot_new.rotcadence=cadence;

 rot_rel.rotbad=0;
 rot_rel.rotpairs=npairstot;
 rot_rel.rotcadence=cadence;

 write_flatfields(filename_flatfield_out, arr_flat, arr_flat_new, camera, recnum, tobs_link, t_0, focus, rot_new, rot_cur);
 write_flatfields(filename_flatfield_rel, arr_flat, arr_flat_rel, camera, recnum, tobs_link, t_0, focus, rot_rel, rot_cur);
  }
 else
 
   {
      printf("No data records found\n");
   }     



 free(flatfield);
 free(flatfield_new);
 free(offpoint);
 free(bad);
 free(flat);

 free(ffhp);
 free_fresize(&fresizes);


 printf("COMLETED!\n");
 return 0;

}











