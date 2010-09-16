/*-------------------------------------------------------------------------------------------------------*/
/* ANSI 99 C CODE TO CALCULATE THE POLYNOMIAL COEFFICIENTS TO CORRECT THE DOPPLER VELOCITIES RETURNED    */
/* BY THE MDI-LIKE ALGORITHM                                                                             */
/*-------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <HMIparam.h>           //contains definitions for some HMI filter parameters
#include <mkl.h>

char *module_name    = "correction_velocities";   //name of the module
#define kRecSetIn      "begin"        //beginning time for which an output is wanted. MANDATORY PARAMETER.
#define kRecSetIn2     "end"          //end time for which an output is wanted. MANDATORY PARAMETER.
#define kTypeSetIn     "levin"        //series name of the input data
#define kTypeSetOut    "levout"       //series name of the output data


//convention for light and dark frames for keyword HCAMID
#define LIGHT_SIDE  2   //SIDE CAMERA
#define LIGHT_FRONT 3   //FRONT CAMERA
#define DARK_SIDE   0   //SIDE CAMERA
#define DARK_FRONT  1   //FRONT CAMERA

//arguments of the module
ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn,  "",  "beginning time for which an output is wanted"},
     {ARG_STRING, kRecSetIn2, "",  "end time for which an output is wanted"},
     {ARG_STRING, kTypeSetIn, "",  "series name of input data"},
     {ARG_STRING, kTypeSetOut,"",  "series name of output data"},
     {ARG_END}
};


/*-----------------------------------------------------------------------------------------------------*/
/* Function to perform linear interpolation                                                            */
/* found on the internet                                                                               */
/* returns the values yinterp at points x of 1D function yv (at the points xv)                         */
/*-----------------------------------------------------------------------------------------------------*/

void lininterp1f(double *yinterp, double *xv, double *yv, double *x, double ydefault, int m, int minterp)
{
    int i, j; 
    int nrowsinterp, nrowsdata;
    nrowsinterp = minterp;
    nrowsdata = m;
    for (i=0; i<nrowsinterp; i++)
      {
	if((x[i] < xv[0]) || (x[i] > xv[nrowsdata-1])) yinterp[i] = ydefault;
	else
	  {   
	    for(j=1; j<nrowsdata; j++)
	      {      
		if(x[i]<=xv[j])
		  {		   
		    yinterp[i] = (x[i]-xv[j-1]) / (xv[j]-xv[j-1]) * (yv[j]-yv[j-1]) + yv[j-1];
		    break;
		  }
	      }
	  }
      }
} 


/*------------------------------------------------------------------------------------------------------*/
/*                                                                                                      */
/*  MAIN PROGRAM                                                                                        */
/*                                                                                                      */
/*------------------------------------------------------------------------------------------------------*/



int DoIt(void) {

  int errbufstat    =setvbuf(stderr, NULL, _IONBF, BUFSIZ);                     //for debugging purpose when running on the cluster
  int outbufstat    =setvbuf(stdout, NULL, _IONBF, BUFSIZ);

  char *inRecQuery  = cmdparams_get_str(&cmdparams, kRecSetIn,      NULL);      //beginning time
  char *inRecQuery2 = cmdparams_get_str(&cmdparams, kRecSetIn2,     NULL);      //end time
  char *inLev       = cmdparams_get_str(&cmdparams, kTypeSetIn,     NULL);      //input series
  char *outLev      = cmdparams_get_str(&cmdparams, kTypeSetOut,    NULL);      //output series

  int error=0;


  int malign=32;
  a=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  coeffd=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  char uplo[] = "U";
  int ngood,info;
  int ione = 1;

  dpotrf(uplo,&ngood,a,&ngood,&info); // Cholesky decomposition
  dpotrs(uplo,&ngood,&ione,a,&ngood,coeffd,&ngood,&info);

  MKL_free(a);
  MKL_free(coeffd);

  return error;
  
}
