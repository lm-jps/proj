/*-----------------------------------------------------------------------------------------*/
/*                                                                                         */
/* Program to compute Dopplergrams using a MDI-like algorithm                              */
/*(compute the phases of the 1st and 2nd Fourier coefficients)                             */
/* Author: S. Couvidat (based on a code by J. Schou)                                       */
/* Version 1.8 April 20, 2009                                                              */
/*                                                                                         */
/* uses a MDI-like algorithm with 5 or 6 tuning positions                                  */
/* averages the velocities returned by 1st and 2nd Fourier                                 */
/* coefficient, and by LCP and RCP                                                         */
/* the code also estimates the Fe I linewidth, linedepth, and                              */
/* the continuum intensity                                                                 */
/*                                                                                         */
/*-----------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <jsoc_main.h>
#include <omp.h>                      //OpenMP header

#undef I                              //I is the complex number (0,1). We un-define it to avoid confusion

struct parameterDoppler {             //structure to provide some parameters defined in HMIparam.h to Dopplergram()
  float FSRNB;
  float FSRWB;
  float FSRE1;
  float FSRE2;
  float FSRE3;
  float FSRE4;
  float FSRE5;
  float dlamdv;
  int maxVtest;
  int maxNx;
  int ntest;
  float dvtest;
  float MISSINGDATA;
  float MISSINGRESULT;
};


int Dopplergram(DRMS_Array_t **arrLev1p,DRMS_Array_t **arrLev15,int framelistSize,DRMS_Array_t *Lookuptable,float Rsun,float X0,float Y0,struct parameterDoppler DopplerParameters)
{

  float FSR[7];
  float dlamdv;
  int maxVtest;
  int maxNx;
  int ntest;
  float dvtest;
  float MISSINGDATA;
  float MISSINGRESULT;
  int  i,j,iii;                                                       //loop variables


  FSR[0]=DopplerParameters.FSRNB;
  FSR[1]=DopplerParameters.FSRWB;
  FSR[2]=DopplerParameters.FSRE1;
  FSR[3]=DopplerParameters.FSRE2;
  FSR[4]=DopplerParameters.FSRE3;
  FSR[5]=DopplerParameters.FSRE4;
  FSR[6]=DopplerParameters.FSRE5;
  dlamdv=DopplerParameters.dlamdv;
  maxVtest=DopplerParameters.maxVtest;
  maxNx=DopplerParameters.maxNx;
  ntest=DopplerParameters.ntest;
  dvtest=DopplerParameters.dvtest;
  MISSINGDATA=DopplerParameters.MISSINGDATA;
  MISSINGRESULT=DopplerParameters.MISSINGRESULT;
  float vtest[ntest];
  float poly[ntest],poly2[ntest];                                     //for the interpolation of the look-up tables

  int  status      = 0; 
  int  error       = 0;                                               //error code returned by the routine: 0=no error, !=0 means error
  int  nthreads;                                                      //for OpenMP
  int  nRows,nColumns;
  DRMS_Type_t type;
  float distance;

  //check whether or not framelistSize is an even number
  if(framelistSize != 12 && framelistSize != 10)
    {
      printf("Error in subroutine computing Dopplergrams: the framelist size is not 12 or 10\n");
      error=3;
      goto endroutine;
    }

  nthreads=omp_get_num_procs();                                      //number of threads supported by the machine where the code is running
  omp_set_num_threads(nthreads);                                     //set the number of threads to the maximum value
  printf("number of threads= %d\n",nthreads);

  int N=framelistSize/2;                                             //number of wavelengths
  float *cosi;
  float *sini;
  float *cos2i;
  float *sin2i;
  cosi = (float *)malloc(N*sizeof(float));  
  sini = (float *)malloc(N*sizeof(float));  
  cos2i= (float *)malloc(N*sizeof(float));  
  sin2i= (float *)malloc(N*sizeof(float));  
  if(cosi == NULL || sini == NULL || cos2i == NULL || sin2i == NULL)
    {
      printf("Error: memory could not be allocated in Dopplergram() function\n");
      exit(EXIT_FAILURE);
    }
  
  //check that all level 1p data have the same data type
  for(i=0;i<framelistSize;++i)
    {
      type    = arrLev1p[i]->type;                                  //float if produced by Jesper's routine
      if (type != DRMS_TYPE_FLOAT)                                  //if the type is not FLOAT
	{
	  printf("Error in subroutine computing Dopplergrams: data type of level 1p data is not FLOAT\n");
	  error = 1;
	  free(cosi);
	  free(sini);
	  free(cos2i);
	  free(sin2i);
	  goto endroutine;
	}
    }
  

  nRows    = arrLev1p[0]->axis[1];
  nColumns = arrLev1p[0]->axis[0];

  //check that all filtergrams have the same dimensions
  for(i=1;i<framelistSize;++i)
    {
      if (arrLev1p[i]->axis[1] != nRows || arrLev1p[i]->axis[0] != nColumns )//if filtergrams are not nRows*nColumns
	{
	  printf("Error in subroutine computing Dopplergrams: dimensions of level 1p data are not %d x %d \n",nRows,nColumns);
	  error = 2;
	  free(cosi);
	  free(sini);
	  free(cos2i);
	  free(sin2i);
	  goto endroutine;
	}
    }


  float *lam0g  = arrLev15[0]->data ;                                 //Dopplergram
  float *B0g    = arrLev15[1]->data;                                  //magnetogram
  float *Idg    = arrLev15[2]->data;                                  //linedepth
  float *widthg = arrLev15[3]->data;                                  //linewidth
  float *I0g    = arrLev15[4]->data;                                  //continuum

  memset(lam0g,  0.0, drms_array_size(arrLev15[0]));                    //fill the observable arrays with 0
  memset(B0g  ,  0.0, drms_array_size(arrLev15[1]));
  memset(Idg  ,  0.0, drms_array_size(arrLev15[2]));
  memset(widthg, 0.0, drms_array_size(arrLev15[3]));
  memset(I0g   , 0.0, drms_array_size(arrLev15[4]));
  
  
  //variables for the MDI-like algorithm
  float FSRNB   = FSR[0];                       //FSR Narrow-Band Michelson, in Angstroms (WILL CHANGE ONCE THE VALUE IS ACCURATELY MEASURED)
  float dtune   = FSRNB/2.5;                    //wavelength separation between each tuning position, about 69 mA  (SHOULD BE THE SAME FOR 5 OR 6 WAVELENGTHS)
  float dv      = 1.0/dlamdv;                   //conversion factor from wavelength to velocity
  float dvtune  = dtune*dv;
  float *tune   = NULL, *angle=NULL;
  tune=(float *)malloc(N*sizeof(float));
  if(tune == NULL)
    {
      printf("Error: unable to allocate memory to tune\n");
      exit(EXIT_FAILURE);
    }
  angle=(float *)malloc(N*sizeof(float));
  if(angle == NULL)
    {
      printf("Error: unable to allocate memory to angle\n");
      exit(EXIT_FAILURE);
    }


  if(N == 6)
    {
      tune[0]=-2.5;//tuning positions
      tune[1]=-1.5;
      tune[2]=-0.5;
      tune[3]= 0.5;
      tune[4]= 1.5;
      tune[5]= 2.5;
      for(i=0;i<N;++i) angle[i]=tune[i];
    }
  if(N == 5)
    {
      tune[0]=-2.0;
      tune[1]=-1.0;
      tune[2]= 0.0;
      tune[3]= 1.0;
      tune[4]= 2.0;
      for(i=0;i<N;++i) angle[i]=tune[i];
    }

  float period = (float)(N-1)*dtune;
  float *L=NULL,*R=NULL;
  L=(float *)malloc(N*sizeof(float));
  R=(float *)malloc(N*sizeof(float));
  if(L == NULL || R == NULL)
    {
      printf("Error: unable to allocate memory for L and/or N\n");
      exit(EXIT_FAILURE);
    }
  float pv1,pv2;
  float magnetic;
  float f1LCPc,f1RCPc,f1LCPs,f1RCPs,f2LCPc,f2RCPc,f2LCPs,f2RCPs;
  float vLCP,vRCP,v2LCP,v2RCP;
  float temp,temp2,temp3,tempbis,temp2bis,temp3bis;
  float meanL=0.0,meanR=0.0;

  pv1 = dvtune*(float)(N-1);
  pv2 = pv1/2.;
  magnetic = 1.0/(2.0*4.67e-5*0.000061733433*2.5*299792458.0);      //Lande factor=2.5 for Fe I line
  
  for(i=0;i<N;++i) 
    {
      tune[i] = tune[i]*dtune;
      angle[i]= angle[i]*2.0*M_PI/(float)N;
      cosi[i] = cos(angle[i]);
      sini[i] = sin(angle[i]);
      cos2i[i]= cos(2.0*angle[i]);
      sin2i[i]= sin(2.0*angle[i]);
    }
  
  

  /*  cosi[0]=-0.84580903;
  cosi[1]=0.026937725;
  cosi[2]= 0.86905351;
  cosi[3]= 0.86905351;
  cosi[4]=0.026937725;
  cosi[5]=-0.84580903;


  sini[0]=-0.53348578;
  sini[1]=-0.99963711;
  sini[2]=-0.49471809;
  sini[3]=0.49471809;
  sini[4]=0.99963711;
  sini[5]=0.53348578;

  cos2i[0]=0.43078584;
  cos2i[1]=-0.99854872;
  cos2i[2]=0.51050802;
  cos2i[3]=0.51050802;
  cos2i[4]=-0.99854872;
  cos2i[5]=0.43078584;

  sin2i[0]=0.90245419;
  sin2i[1]=-0.053855899;
  sin2i[2]=-0.85987299;
  sin2i[3]=0.85987299;
  sin2i[4]=0.053855899;
  sin2i[5]=-0.90245419;*/



  //array containing the look-up table (type FLOAT)
  float *lookupt = Lookuptable->data ;
  int axist[3];    //dimensions of the look-up tables
  axist[0]=Lookuptable->axis[0];
  axist[1]=Lookuptable->axis[1];     
  axist[2]=Lookuptable->axis[2];
  if(axist[0] > maxVtest || axist[1] > maxNx || axist[2] > maxNx)
    {
      printf("Error in subroutine computing Dopplergrams: dimensions of the look-up tables must not exceed 1002x256x256\n"); //if there is a problem
      error = 4;
      free(cosi);
      free(sini);
      free(cos2i);
      free(sin2i);
      free(tune);
      free(angle);
      free(L);
      free(R);
      goto endroutine;        
    }
  printf("Dimensions of the look-up tables: %d %d %d\n",axist[0],axist[1],axist[2]);

  //we reconstruct the test (input) velocities used to produce the look-up tables
  //WARNING: MUST BE THE SAME AS IN lookup.c
  for(i=0;i<ntest;++i) vtest[i] = dvtest*((float)i-((float)ntest-1.0)/2.0); 

  int   index_lo=0;
  int   index_hi=ntest-1;
  float RR1,RR2;     //variables for the bilinear interpolation of the look-up tables
  int   x0,y0,x1,y1;
  float x,y,xa,xb,ya,yb;
  int   ratio,indexL,indexR,indexL2,indexR2,row,column,step=10;
  long  loc1,loc2,loc3,loc4; //coded 8 bytes
  float Kfourier;
  float *tempvec=NULL;
  Kfourier = dtune/period*2.0;
  ratio = nRows/axist[1];

  /***********************************************************************************************************/
  /*LOOP OVER ALL THE PIXELS OF THE FILTERGRAMS                                                              */
  /***********************************************************************************************************/

  FILE *fp;
  fp=fopen("temp.txt","w");

#pragma omp parallel for default(none) shared(step,arrLev1p,cosi,sini,cos2i,sin2i,pv1,pv2,index_lo,index_hi,vtest,period,dtune,dv,I0g,B0g,Idg,lam0g,widthg,tunem,magnetic,axist,ratio,lookupt,nRows,nColumns,MISSINGDATA,MISSINGRESULT,Kfourier,Rsun,X0,Y0,ntest) private(tempvec,iii,L,R,f1LCPc,f1RCPc,f1LCPs,f1RCPs,vLCP,vRCP,f2LCPc,f2RCPc,f2LCPs,f2RCPs,temp,tempbis,temp2,temp2bis,temp3,temp3bis,meanL,meanR,v2LCP,v2RCP,x0,y0,x,y,x1,y1,RR1,RR2,i,loc1,loc2,loc3,loc4,xa,xb,ya,yb,indexL,indexR,indexL2,indexR2,poly,poly2,row,column,distance)
      for(iii=0;iii<nRows*nColumns;++iii)
	{
	  ///iii=column+row*nColumns
	  row   =iii / nColumns;
	  column=iii % nColumns;

	  distance = sqrt(((float)row-Y0)*((float)row-Y0)+((float)column-X0)*((float)column-X0));
	  //printf("%f %f\n",distance,Rsun);

	  tempvec=arrLev1p[0]->data;
	      if (tempvec[iii] != MISSINGDATA && tempvec[iii] > 0.0 && distance <= Rsun) //only check the 1st filtergram for missing data? 
		{

		  /*-------------------------------------------------------------*/
		  /* MDI-like algorithm                                          */
		  /*-------------------------------------------------------------*/

		  for(i=0;i<N;++i) 
		    {
		      tempvec=arrLev1p[i*2]->data;   //LCP
		      L[i]= tempvec[iii] ;
		      tempvec=arrLev1p[i*2+1]->data; //RCP
		      R[i]= tempvec[iii] ;
		    }

		  //First and Second Fourier coefficients
		  f1LCPc=0.0;
		  f1RCPc=0.0;
		  f1LCPs=0.0;
		  f1RCPs=0.0;
		  f2LCPc=0.0;
		  f2RCPc=0.0;
		  f2LCPs=0.0;
		  f2RCPs=0.0;
		  for(i=0;i<N;++i) 
		    {
		      f1LCPc += cosi[i]*L[i];
		      f1RCPc += cosi[i]*R[i];
		      f1LCPs += sini[i]*L[i];
		      f1RCPs += sini[i]*R[i];
		      f2LCPc += cos2i[i]*L[i];
		      f2RCPc += cos2i[i]*R[i];
		      f2LCPs += sin2i[i]*L[i];
		      f2RCPs += sin2i[i]*R[i];
		    }
 		  
		  //printf("%f %f %f %f\n",f1LCPc,f1RCPc,f1LCPs,f1RCPs);
		  vLCP    = atan2(-f1LCPs,-f1LCPc)*pv1/2.0/M_PI; //-f1LCPs and -f1LCPc because the first Fourier coefficient is -2/T*...
		  vRCP    = atan2(-f1RCPs,-f1RCPc)*pv1/2.0/M_PI;
      		  v2LCP   = atan2(-f2LCPs,-f2LCPc)*pv2/2.0/M_PI;
		  v2RCP   = atan2(-f2RCPs,-f2RCPc)*pv2/2.0/M_PI;
		  v2LCP   = fmod((v2LCP-vLCP+10.5*pv2),pv2)-pv2/2.0+vLCP; //we use the uncorrected velocity, i.e. phase, of the 1st Fourier coefficient to correct for the estimate of v2LCP and v2RCP, because the range of velocities obtained with the second Fourier coefficient is half the range of the first Fourier coefficient
		  v2RCP   = fmod((v2RCP-vRCP+10.5*pv2),pv2)-pv2/2.0+vRCP; 
			
	  
		  /*-------------------------------------------------------------*/
		  /* bilinear interpolation of the look-up tables at pixel (x,y) */
		  /*-------------------------------------------------------------*/
		  		  
		  //NB: it depends on how the filtergrams rebinning from 4096*4096 to axist[1]*axist[2] was done in phasemaps.c
		  //find the 4 neighbors (x0,y0), (x0,y1), (x1,y0), and (x1,y1) of (column,row) on the grid of the look-up tables, and deal with boundary problems

		  x0 = (column/ratio); //a column number
		  y0 = (row/ratio);    //a row number
		  x1 = x0+1;
		  y1 = y0+1;
		  y  = (float)(row % ratio)   /(float)ratio+(float)y0;
		  x  = (float)(column % ratio)/(float)ratio+(float)x0;

		  if(x1 >= axist[1])
		    {
		      x0 = x0-1;
		      x1 = x1-1;
		    }
		  if(y1 >= axist[2])
		    {
		      y0 = y0-1;
		      y1 = y1-1;
		    }
		  
		  //perform the bilinear interpolation
		  //NB: x0*axist[0]+y0*axist[0]*axist[1] cannot be larger than 256x256x1002 and therefore can be coded on 4 bytes 
		  loc1=x0*axist[0]+y0*axist[0]*axist[1];
		  loc2=x1*axist[0]+y0*axist[0]*axist[1];
		  loc3=x0*axist[0]+y1*axist[0]*axist[1];
		  loc4=x1*axist[0]+y1*axist[0]*axist[1];
		  
		  xa=((float)x1-x);
		  xb=(x-(float)x0);
		  ya=((float)y1-y);
		  yb=(y-(float)y0);
		  
		  indexL =ntest+1;
		  indexR =ntest+1;
		  indexL2=ntest+1;
		  indexR2=ntest+1;

		  //initializes arrays (NOT NEEDED ACTUALLY)
		  memset(poly,0.0,sizeof(poly));
		  memset(poly2,0.0,sizeof(poly2));

		  poly[0]  =ya*(lookupt[loc1]*xa+lookupt[loc2]*xb)+yb*(lookupt[loc3]*xa+lookupt[loc4]*xb);                          //for 1st Fourier coefficient
		  poly2[0] =ya*(lookupt[loc1+ntest]*xa+lookupt[loc2+ntest]*xb)+yb*(lookupt[loc3+ntest]*xa+lookupt[loc4+ntest]*xb);  //for 2nd Fourier coefficient
		  for(i=step;i<ntest;i=i+step) //make sure (ntest-1)/step is an integer
		    {
		      poly[i]  =ya*(lookupt[loc1+i]*xa+lookupt[loc2+i]*xb)+yb*(lookupt[loc3+i]*xa+lookupt[loc4+i]*xb);                          //for 1st Fourier coefficient
		      poly2[i] =ya*(lookupt[loc1+i+ntest]*xa+lookupt[loc2+i+ntest]*xb)+yb*(lookupt[loc3+i+ntest]*xa+lookupt[loc4+i+ntest]*xb);  //for 2nd Fourier coefficient
		      if(poly[i] > vLCP && poly[i-step]   <= vLCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {
			      poly[j]  =ya*(lookupt[loc1+j]*xa+lookupt[loc2+j]*xb)+yb*(lookupt[loc3+j]*xa+lookupt[loc4+j]*xb);                          //for 1st Fourier coefficient
			      if(poly[j] > vLCP && poly[j-1]   <= vLCP) indexL = j-1;
			    }
			}
		      if(poly[i] > vRCP && poly[i-step]   <= vRCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {
			      poly[j]  =ya*(lookupt[loc1+j]*xa+lookupt[loc2+j]*xb)+yb*(lookupt[loc3+j]*xa+lookupt[loc4+j]*xb);                          //for 1st Fourier coefficient
			      if(poly[j] > vRCP && poly[j-1]   <= vRCP) indexR = j-1;
			    }
			}
		      if(poly2[i]> v2LCP && poly2[i-step] <= v2LCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {
			      poly2[j] =ya*(lookupt[loc1+j+ntest]*xa+lookupt[loc2+j+ntest]*xb)+yb*(lookupt[loc3+j+ntest]*xa+lookupt[loc4+j+ntest]*xb);  //for 2nd Fourier coefficient
			      if(poly2[j] > v2LCP && poly2[j-1]   <= v2LCP) indexL2 = j-1;
			    }
			}
		      if(poly2[i]> v2RCP && poly2[i-step] <= v2RCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {
			      poly2[j] =ya*(lookupt[loc1+j+ntest]*xa+lookupt[loc2+j+ntest]*xb)+yb*(lookupt[loc3+j+ntest]*xa+lookupt[loc4+j+ntest]*xb);  //for 2nd Fourier coefficient
			      if(poly2[j] > v2RCP && poly2[j-1]   <= v2RCP) indexR2 = j-1;
			    } 
			}
		    }

		  if(indexL == ntest+1 || indexR == ntest+1 || indexL2 == ntest+1 || indexR2 == ntest+1)
		    {
		      lam0g[iii]  = MISSINGRESULT;
		      B0g[iii]    = MISSINGRESULT;
		      widthg[iii] = MISSINGRESULT;
		      Idg[iii]    = MISSINGRESULT;
		      I0g[iii]    = MISSINGRESULT;
		    }
		  else
		    {

		      if(column%128 == 0 || row%128 == 0)
			{
			  fprintf(fp,"%d %d %f %f %f %f",column,row,vLCP,vRCP,v2LCP,v2RCP);
			}


		      //We linearly interpolate in the look-up table for the 1st Fourier coefficient to retrieve the actual velocities
		      vLCP   = vtest[indexL] +(vLCP-poly[indexL])   *(vtest[indexL+1] -vtest[indexL]) /(poly[indexL+1]  -poly[indexL]  );
		      vRCP   = vtest[indexR] +(vRCP-poly[indexR])   *(vtest[indexR+1] -vtest[indexR]) /(poly[indexR+1]  -poly[indexR]  );
		      //We linearly interpolate in the look-up table for the 2nd Fourier coefficient to retrieve the actual velocities
		      v2LCP  = vtest[indexL2]+(v2LCP-poly2[indexL2])*(vtest[indexL2+1]-vtest[indexL2])/(poly2[indexL2+1]-poly2[indexL2]);
		      v2RCP  = vtest[indexR2]+(v2RCP-poly2[indexR2])*(vtest[indexR2+1]-vtest[indexR2])/(poly2[indexR2+1]-poly2[indexR2]);
		    

		      if(column%128 == 0 || row%128 == 0)
			{
			  fprintf(fp," %f %f %f %f\n",vLCP,vRCP,v2LCP,v2RCP);
			}

		      /*-------------------------------------------------------------*/
		      /* Calculation of observables                                  */
		      /*-------------------------------------------------------------*/
		      
		      f1LCPc = f1LCPc*Kfourier;
		      f1RCPc = f1RCPc*Kfourier;
		      f1LCPs = f1LCPs*Kfourier;
		      f1RCPs = f1RCPs*Kfourier;
		      f2LCPc = f2LCPc*Kfourier;
		      f2RCPc = f2RCPc*Kfourier;
		      f2LCPs = f2LCPs*Kfourier;
		      f2RCPs = f2RCPs*Kfourier;
		      
		      //We compute the Doppler velocity
		      lam0g[iii]  =-(vLCP+vRCP+v2LCP+v2RCP)/4.;//simple average. Need weights? SIGN CONVENTION: HERE, v<0 FOR MOTION TOWARD THE OBSERVER (BLUESHIFT)
		      //printf("lam0g= %f\n",lam0g[iii]);
		      
		      //We compute the l.o.s. magnetic field
		      B0g[iii]    = (vLCP-vRCP+v2LCP-v2RCP)/2.0*magnetic;
		      
		      //We compute the linewidth (in Angstroms)
		      temp        = period/M_PI*sqrt(1.0/6.0*log((f1LCPc*f1LCPc+f1LCPs*f1LCPs)/(f2LCPc*f2LCPc+f2LCPs*f2LCPs)));
		      tempbis     = period/M_PI*sqrt(1.0/6.0*log((f1RCPc*f1RCPc+f1RCPs*f1RCPs)/(f2RCPc*f2RCPc+f2RCPs*f2RCPs)));
		      widthg[iii] = (temp+tempbis)*sqrt(log(2.0)); //we want the FWHM not the sigma of the Gaussian	      
		      //We compute the linedepth
		      temp2       = period/2.0*sqrt(f1LCPc*f1LCPc+f1LCPs*f1LCPs)/sqrt(M_PI)/temp*exp(M_PI*M_PI*temp*temp/period/period);
		      temp2bis    = period/2.0*sqrt(f1RCPc*f1RCPc+f1RCPs*f1RCPs)/sqrt(M_PI)/tempbis*exp(M_PI*M_PI*tempbis*tempbis/period/period);
		      Idg[iii]    = (temp2+temp2bis)/2.0;
		      
		      //We compute the continuum intensity
		      temp3       = (vLCP+v2LCP)/2.0/dv;
		      temp3bis    = (vRCP+v2RCP)/2.0/dv;
		      meanL=0;
		      meanR=0;
		      for(i=0;i<N;++i)
			{
			  meanL      += L[i];
			  meanR      += R[i];
			}
		      meanL=meanL/(float)N;
		      meanR=meanR/(float)N;
		      for(i=0;i<N;++i)
			{
			  meanL      += temp2   /(float)N*exp(-(tune[i]-temp3)   *(tune[i]-temp3)   /temp   /temp);
			  meanR      += temp2bis/(float)N*exp(-(tune[i]-temp3bis)*(tune[i]-temp3bis)/tempbis/tempbis);
			}
		      
		      I0g[iii]    = (meanL+meanR)/2.0;
		    }		  
		  
		}//if (tempvec[iii] != MISSINGDATA && tempvec[iii] != 0.0 && distance <= Rsun)
	      else
		{
		  lam0g[iii]  = MISSINGRESULT;
		  B0g[iii]    = MISSINGRESULT;
		  widthg[iii] = MISSINGRESULT;
		  Idg[iii]    = MISSINGRESULT;
		  I0g[iii]    = MISSINGRESULT;
		}
	      
	}//for iii
  
      free(cosi);
      free(sini);
      free(cos2i);
      free(sin2i);
      free(tune);
      free(angle);
      free(L);
      free(R);

 endroutine:

      fclose(fp);

      return error;
      
}//end routine
