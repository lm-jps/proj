
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>
#include </home/jsoc/include/fftw3.h>
#include <omp.h>
#include <module_flatfield.h>

char *module_name  = "module_flatfield_combine";    //name of the module

#define kRecSetIn      "input_series" //name of the series containing the input filtergrams
#define fsnname      "fsn"


ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn, "",  "Input data series."},
     {ARG_END}
};


int get_flatfields(int camera, TIME t_0, int focus,
                   float *flatfield, float *offpoint, 
                   short *bad, long long recnum[6], TIME tobs_link[2],
		   struct rotpar *rot_cur);


int write_flatfields(DRMS_Array_t *arr_flat, DRMS_Array_t *arrout_new, int camera,
                     long long recnum[6], TIME tobs_link[2], int focus,
                     struct rotpar rot_new, 
                     struct rotpar rot_cur);

int DoIt(void)
{

#include "module_flatfield_const.h"

  const char *inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL); //cmdparams is defined in jsoc_main.h
 
  int  status= DRMS_SUCCESS, stat=DRMS_SUCCESS;

  DRMS_Segment_t *segin = NULL;


  DRMS_Type_t type_time = DRMS_TYPE_TIME;
  DRMS_Type_t type_int = DRMS_TYPE_INT;
  DRMS_Type_t type_float = DRMS_TYPE_FLOAT;
  DRMS_Type_t type_double = DRMS_TYPE_DOUBLE;

  DRMS_RecordSet_t *data, *dataout;
  DRMS_Record_t  *ff_last_front=NULL, *ff_last_side=NULL;
  DRMS_Array_t *arr_flat, *arr_flat_new;

  int axisout[2]={nx,ny};

  TIME t_0;
  int focus, camera;
  float cadence;
  int npairstot=0;

  int i,j,k;
  float h1,h2;


  float *flatfield, *offpoint, *flat;
  short *bad;


  float *flatfield_new;
  float *flatx, *flaty, *ff, *ff_new, *fx, *fy;

  long long recnum[6];
  TIME tobs_link[2];

  struct rotpar rot_cur, rot_new;


  flatfield=(float *)(malloc(nx*ny*sizeof(float)));
  flatfield_new=(float *)(malloc(nx*ny*sizeof(float)));
  offpoint=(float *)(malloc(nx*ny*sizeof(float)));
  bad=(short *)(malloc(nx*ny*sizeof(short)));
  flatx=(float *)(calloc(nx*ny,sizeof(float)));
  flaty=(float *)(calloc(nx*ny,sizeof(float)));
  fx=(float *)(malloc(nx*ny*sizeof(float)));
  fy=(float *)(malloc(nx*ny*sizeof(float)));
  flat=(float *)(malloc(nx*ny*sizeof(float)));

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


      data     = drms_open_records(drms_env,(char *)inRecQuery,&stat);

 if (data == NULL){printf("can not open records\n"); exit(EXIT_FAILURE);}

  int nRecs=data->n;
  printf("number of records %d\n", nRecs);



if (stat == DRMS_SUCCESS && data != NULL && nRecs > 0)
  {
  int keyvalue_fid[nRecs];
  DRMS_Record_t *rec0[nRecs];
  float *flatfield_fid[nRecs];
  DRMS_Array_t *arrin0[nRecs];
  short wave[nRecs];
  short pol[nRecs];
  short index[nRecs];
  int npairs[nRecs];


for (k=0; k<nRecs; ++k)
  {
    rec0[k]=data->records[k];
    keyvalue_fid[k] = drms_getkey_int(rec0[k],fidkey,&status);

    segin    = drms_segment_lookupnum(rec0[k], 0);
    arrin0[k]= drms_segment_read(segin, type_float, &status);
    flatfield_fid[k] = arrin0[k]->data;
    wave[k]=((keyvalue_fid[k]-10000)/10-5)/2;
    pol[k]=(keyvalue_fid[k] % 10)-8;
    index[k]=wave[k]*2+pol[k];

    npairs[k]= drms_getkey_int(rec0[k],keynpairs,&status);
    npairstot=npairstot+npairs[k];
    printf("%d %d %d %d %d\n", k, index[k], keyvalue_fid[k], pol[k], wave[k]);
  }

 t_0=drms_getkey_time(rec0[nRecs-1],keytstart,&status);
 focus = drms_getkey_int(rec0[nRecs-1],keyfocusflat,&status);
 camera= drms_getkey_int(rec0[nRecs-1],keycamera,&status);
 cadence=drms_getkey_float(rec0[nRecs-1],keycadence,&status);

 printf("%lf %d %d %f\n", t_0, focus, camera, cadence);

 get_flatfields(camera, t_0, focus,  flatfield, 
		     offpoint,  bad, recnum, tobs_link, &rot_cur);  //get front camera flatfield


for (k=0; k<nRecs; ++k)
  {
    ff=flatfield_fid[k];

    for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
      {
    flatx[j*nx+i]=flatx[j*nx+i]+cof_combx[index[k]]*ff[j*nx+i];
    flaty[j*nx+i]=flaty[j*nx+i]+cof_comby[index[k]]*ff[j*nx+i];
    flat[j*nx+i]=flatfield[j*nx+i]/offpoint[j*nx+i];
      }
  }

 drms_close_records(data,DRMS_FREE_RECORD);

  FILE *fileptr;                                                                                                                          
         
   for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
  {
    fx[j*nx+i]=1.0;
    if(i > 0)
      {
	if (flatx[j*nx+i] != 1.0 && flatx[j*nx+i-1] != 1.0)
	  {
	    h2=flatx[j*nx+i]-flatx[j*nx+i-1];
	  }
	else
	  {
	    h2=0.0;
	  }
	    fx[j*nx+i-1]=flat[j*nx+i-1]+h2;
      }
  }
  
 

for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
  {
    fy[j*nx+i]=1.0;
    if (j>0){
    if (flaty[j*nx+i] != 1.0 && flaty[(j-1)*nx+i] != 1.0)
      {
	h1=flaty[j*nx+i]-flaty[(j-1)*nx+i];
      }
    else
      {
	h1=0.0;
      }
	fy[(j-1)*nx+i]=flat[(j-1)*nx+i]+h1;
      
    }
  }

 int count=0;
for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)
  {
    flatfield_new[j*nx+i]=0.5*(fx[j*nx+i]+fy[j*nx+i])*offpoint[j*nx+i];
   if (flatfield_new[j*nx+i] != flatfield[j*nx+i]) ++count;
  }



   fileptr = fopen ("/tmp20/richard/interpol/dd.bin", "w");                                             
	   fwrite((char*)(fy),sizeof(float),nx*ny,fileptr);
	   fclose(fileptr);


 arr_flat=drms_array_create(type_float,2,axisout,NULL,&status);
 arr_flat_new=drms_array_create(type_float,2,axisout,NULL,&status);

 ff_new=arr_flat_new->data;
 ff=arr_flat->data;

 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)ff[j*nx+i]=flatfield[j*nx+i];
 for (j=0; j<ny; ++j) for (i=0; i<nx; ++i)ff_new[j*nx+i]=flatfield_new[j*nx+i];

rot_new.rotbad=count;
rot_new.rotpairs=npairstot;
rot_new.rotcadence=cadence;

 write_flatfields(arr_flat, arr_flat_new, camera, recnum, tobs_link, focus, rot_new, rot_cur);

}
 else
 
   {
      printf("No data records found\n");
   }     



 free(flatfield);
 free(flatfield_new);
 free(offpoint);
 free(bad);
 free(flatx);
 free(flaty);
 free(fx);
 free(fy);

 printf("COMLETED!\n");
 return 0;

}











