/*-----------------------------------------------------------------------------------------*/
/*                                                                                         */
/* Program to compute Dopplergrams using a MDI-like algorithm                              */
/*(compute the phases of the 1st and 2nd Fourier coefficients)                             */
/* Author: S. Couvidat (based on a code by J. Schou)                                       */
/* Version 1.8 April 20, 2009                                                              */
/*                                                                                         */
/* USES BICUBIC INTERPOLATION UNLIKE Dopplergram.c                                         */
/* uses a MDI-like algorithm with 5 or 6 tuning positions                                  */
/* averages the velocities returned by 1st and 2nd Fourier                                 */
/* coefficient, and by LCP and RCP                                                         */
/* the code also estimates the Fe I linewidth, linedepth, and                              */
/* the continuum intensity                                                                 */
/*                                                                                         */
/* FOR THE MEANING OF ROWS AND COLUMNS WE FOLLOW THE CONVENTION:                           */
/* A 2D ARRAY A IN C IS DEFINED AS A[row][column]                                          */
/*                                                                                         */
/* CONVENTION FOR THE VELOCITIES:                                                          */
/* POSITIVE VELOCITIES CORRESPOND TO REDSHIFT (MOVEMENTS AWAY FROM THE OBSERVER)           */
/*-----------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <jsoc_main.h>
#include <omp.h>                      //OpenMP header

#undef I                              //I is the complex number (0,1). We un-define it to avoid confusion

struct parameterDoppler {             //structure to provide some parameters defined in HMIparam.h to Dopplergram()
  double FSRNB;
  double FSRWB;
  double FSRE1;
  double FSRE2;
  double FSRE3;
  double FSRE4;
  double FSRE5;
  double dlamdv;
  int maxVtest;
  int maxNx;
  int ntest;
  double dvtest;
  float MISSINGDATA;
  float MISSINGRESULT;
};


/*-------------------------------------------------------------------*/
/* bicubic interpolation routine derived from Numerical Recipes in C */
/*-------------------------------------------------------------------*/

void bcucof(float y[4],float y1[4],float y2[4],float y12[4],float d1,float d2,float c[4][4])
 {
   static int wt[16][16]= {
     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
   
   int l,k,j,i;
   float xx,d1d2,cl[16],x[16];

   d1d2=d1*d2;
   for(i=0;i<=3;i++){
     x[i]=y[i];
     x[i+4]=y1[i]*d1;
     x[i+8]=y2[i]*d2;
     x[i+12]=y12[i]*d1d2;
   }
   for(i=0;i<=15;i++){
     xx=0.0;
     for(k=0;k<=15;k++) xx += wt[i][k]*x[k];
     cl[i]=xx;
   }
   l=0;
   for(i=0;i<=3;i++)
     for(j=0;j<=3;j++) c[i][j]=cl[l++];

 }


/*-----------------------------------------------------------------------------------------*/
/* MDI-like algorithm                                                                      */
/*-----------------------------------------------------------------------------------------*/


int Dopplergram2(DRMS_Array_t **arrLev1p,DRMS_Array_t **arrLev15,int framelistSize,DRMS_Array_t *Lookuptable,float Rsun,float X0,float Y0,struct parameterDoppler DopplerParameters,int *MISSVALS,int *SATVALS,float cdelt1)
{

  double FSR[7];
  double dlamdv;
  int maxVtest;
  int maxNx;
  int ntest;
  double dvtest;
  float MISSINGDATA;
  float MISSINGRESULT;
  int  i,j,iii,k;                                                      //loop variables
  float ExtraCrop = 50.;//we crop the data at Rsun+ExtraCrop in pixels

  double magnetic = 1.0/(2.0*4.67E-5*0.000061733433*2.5*299792458.0); //Lande factor=2.5 for Fe I line

  //parameters of the Fe I line profiles observed by Roger Ulrich (OBTAINED BY CODE centertolimb.pro; MAKE SURE IT'S THE SAME AS IN lookup.c)
  double cost[3]={1.0,0.70710678,0.5};                    //cos(theta) where theta is the angular distance from disk center
  double minimumCoeffs[2]={0.41922611,0.24190794};        //minimum intensity, assuming the continuum is at 1 (result from a Gaussian fit in the range [-0.055,0.055] Angstroms)
  double FWHMCoeffs[2]={151.34559,-58.521771};            //FWHM of the solar line in Angstrom (result from a Gaussian fit in the range [-0.055,0.055] Angstroms)
  double FWHM,minimum,angulardistance,minimumR,minimumL;

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
  double vtest[ntest];                                                //internally calculations are done in double precision
  float poly[ntest],poly2[ntest];                                     //for the spatial interpolation of the look-up tables (THE TABLES ARE ASSUMED TO BE OF TYPE FLOAT)

  int  status      = 0; 
  int  error       = 0;                                               //error code returned by the routine: 0=no error, !=0 means error
  int  nRows,nColumns;
  DRMS_Type_t type;
  float distance;

  //check whether or not framelistSize is an even number
  if(framelistSize != 12 && framelistSize != 10)
    {
      printf("Error in subroutine computing the Dopplergrams: the framelist size is not 12 or 10\n");
      error=3;
      return error;
    }

  int N=framelistSize/2;                                             //number of different wavelengths
  double *cosi;
  double *sini;
  double *cos2i;
  double *sin2i;
  cosi = (double *)malloc(N*sizeof(double));  
  sini = (double *)malloc(N*sizeof(double));  
  cos2i= (double *)malloc(N*sizeof(double));  
  sin2i= (double *)malloc(N*sizeof(double));  
  if(cosi == NULL || sini == NULL || cos2i == NULL || sin2i == NULL)
    {
      printf("Error: memory could not be allocated in the Dopplergram() function\n");
      exit(EXIT_FAILURE);
    }
  float L[N],R[N]; //NO MALLOCS BECAUSE THEY ARE WRITTEN INSIDE THE OMP LOOP, THE LEVEL 1P DATA ARE ASSUMED TO BE OF TYPE FLOAT
 

  //check that all level 1p data have the same data type (IS THAT NECESSARY?)
  for(i=0;i<framelistSize;++i)
    {
      type    = arrLev1p[i]->type;                                  //float if produced by Jesper's routine
      if (type != DRMS_TYPE_FLOAT)                                  //if the type is not FLOAT
	{
	  printf("Error in subroutine computing the Dopplergrams: data type of level 1p data is not FLOAT\n");
	  error = 1;
	  free(cosi);
	  free(sini);
	  free(cos2i);
	  free(sin2i);
	  return error;
	}
    }
  
  //with C convention Array[row][column], axis[0] is the number of columns 
  nRows    = arrLev1p[0]->axis[1];
  nColumns = arrLev1p[0]->axis[0];

  //check that all filtergrams have the same dimensions (IS THAT NECESSARY?)
  for(i=1;i<framelistSize;++i)
    {
      if (arrLev1p[i]->axis[1] != nRows || arrLev1p[i]->axis[0] != nColumns )//if filtergrams are not nRows*nColumns
	{
	  printf("Error in subroutine computing the Dopplergrams: dimensions of level 1p data are not %d x %d \n",nRows,nColumns);
	  error = 2;
	  free(cosi);
	  free(sini);
	  free(cos2i);
	  free(sin2i);
	  return error;
	}
    }


  //THE LEVEL 1.5 DATA ARE ASSUMED TO BE OF TYPE FLOAT
  float *lam0g  = arrLev15[0]->data ;                                 //Dopplergram
  float *B0g    = arrLev15[1]->data;                                  //magnetogram
  float *Idg    = arrLev15[2]->data;                                  //linedepth
  float *widthg = arrLev15[3]->data;                                  //linewidth
  float *I0g    = arrLev15[4]->data;                                  //continuum

  memset(lam0g,  0.0, drms_array_size(arrLev15[0]));                  //fill the observable arrays with 0
  memset(B0g  ,  0.0, drms_array_size(arrLev15[1]));
  memset(Idg  ,  0.0, drms_array_size(arrLev15[2]));
  memset(widthg, 0.0, drms_array_size(arrLev15[3]));
  memset(I0g   , 0.0, drms_array_size(arrLev15[4]));
  
  
  //variables for the MDI-like algorithm
  double FSRNB   = FSR[0];                       //FSR Narrow-Band Michelson, in Angstroms (WILL CHANGE ONCE THE VALUE IS ACCURATELY MEASURED)
  double dtune   = FSRNB/2.5;                    //wavelength separation between each tuning position, nominally 68.8 mA  (SHOULD BE THE SAME FOR 5 OR 6 WAVELENGTHS)
  double dv      = 1.0/dlamdv;                   //conversion factor from wavelength to velocity
  double dvtune  = dtune*dv;
  double *tune   = NULL, *angle=NULL;
  tune=(double *)malloc(N*sizeof(double));
  if(tune == NULL)
    {
      printf("Error: unable to allocate memory to tune\n");
      exit(EXIT_FAILURE);
    }
  angle=(double *)malloc(N*sizeof(double));
  if(angle == NULL)
    {
      printf("Error: unable to allocate memory to angle\n");
      exit(EXIT_FAILURE);
    }


  if(N == 6)
    {
      printf("USING 6 WAVELENGTHS\n");
      tune[0]=+2.5;//I0 (location of the filter)
      tune[1]=+1.5;//I1
      tune[2]=+0.5;//I2
      tune[3]=-0.5;//I3
      tune[4]=-1.5;//I4
      tune[5]=-2.5;//I5
      for(i=0;i<N;++i) angle[i]=tune[i];
    }
  if(N == 5)
    {
      printf("USING 5 WAVELENGTHS\n");
      tune[0]=+2.0;
      tune[1]=+1.0;
      tune[2]= 0.0;
      tune[3]=-1.0;
      tune[4]=-2.0;
      for(i=0;i<N;++i) angle[i]=tune[i];
    }

  double period = (double)(N-1)*dtune;
  double pv1,pv2;
  double f1LCPc,f1RCPc,f1LCPs,f1RCPs,f2LCPc,f2RCPc,f2LCPs,f2RCPs;
  double vLCP,vRCP,v2LCP,v2RCP;
  double temp,temp2,temp3,tempbis,temp2bis,temp3bis;
  double meanL=0.0,meanR=0.0;

  pv1 = dvtune*(double)(N-1);
  pv2 = pv1/2.;
  
  for(i=0;i<N;++i) 
    {
      tune[i] = tune[i]*dtune;
      angle[i]= angle[i]*2.0*M_PI/(double)N;
      cosi[i] = cos(angle[i]);
      sini[i] = sin(angle[i]);
      cos2i[i]= cos(2.0*angle[i]);
      sin2i[i]= sin(2.0*angle[i]);
    }
  
  
  //array containing the look-up table (type FLOAT)
  float *lookupt = Lookuptable->data ;
  int axist[3];                  //dimensions of the look-up tables
  axist[0]=Lookuptable->axis[0]; //number of input velocities (ntest*2)
  axist[1]=Lookuptable->axis[1]; //number of columns (256)    
  axist[2]=Lookuptable->axis[2]; //number of rows (256)
  printf("Dimensions of the look-up tables: %d %d %d\n",axist[0],axist[1],axist[2]);

  if(axist[0] > maxVtest || axist[1] > maxNx || axist[2] > maxNx)
    {
      printf("Error in subroutine computing the Dopplergrams: dimensions of the look-up tables exceed what is allowed: %d %d %d\n",maxVtest,maxNx,maxNx); //if there is a problem
      error = 4;
      free(cosi);
      free(sini);
      free(cos2i);
      free(sin2i);
      free(tune);
      free(angle);
      return error;        
    }

  //to accommodate the bilinear interpolation code we format the 1D lookup tables into 2d arrays
  float lookupt2D[Lookuptable->axis[1]][Lookuptable->axis[2]];


  //we reconstruct the test (input) velocities used to produce the look-up tables
  //WARNING: MUST BE THE SAME AS IN lookup.c
  for(i=0;i<ntest;++i) vtest[i] = dvtest*((double)i-((double)ntest-1.0)/2.0); 

  int   index_lo=0;
  int   index_hi=ntest-1;
  float RR1,RR2;     //variables for the bilinear interpolation of the look-up tables
  int   x0,y0,x1,y1;
  float xa,xb,ya,yb;
  int   ratio,indexL,indexR,indexL2,indexR2,row,column,step=10;
  long  loc0,loc1,loc2,loc3,loc4,loc5,loc6,loc7,loc8,loc9,loc10,loc11,loc12,loc13,loc14,loc15; //coded 8 bytes
  double Kfourier;
  float *tempvec=NULL;
  Kfourier = dtune/period*2.0;
  ratio = nRows/axist[2];
  float minlookupt1=0.0,maxlookupt1=0.0,minlookupt2=0.0,maxlookupt2=0.0;
  
  int MISSVALS2=0; //the number of missing values (NaN)
  int SATVALS2 =0; //number of staurated values
  float coeffs[4][4];
  float u,t,yy[4],yy1[4],yy2[4],yy12[4];
  float offset=((float)ratio-1.0)/2.0; //because the phase maps were rebinned using rebin() which means that pixel 0 is actually (ratio-1)/2 on the initial grid

  /***********************************************************************************************************/
  /*LOOP OVER ALL THE PIXELS OF THE FILTERGRAMS                                                              */
  /***********************************************************************************************************/

#pragma omp parallel default(none) reduction(+:MISSVALS2,SATVALS2) shared(step,arrLev1p,cosi,sini,cos2i,sin2i,pv1,pv2,index_lo,index_hi,vtest,period,dtune,dv,I0g,B0g,Idg,lam0g,widthg,magnetic,axist,ratio,lookupt,nRows,nColumns,MISSINGDATA,MISSINGRESULT,Kfourier,Rsun,X0,Y0,ntest,tune,N,cost,minimumCoeffs,FWHMCoeffs,offset,ExtraCrop,cdelt1) private(tempvec,iii,L,R,f1LCPc,f1RCPc,f1LCPs,f1RCPs,vLCP,vRCP,f2LCPc,f2RCPc,f2LCPs,f2RCPs,temp,tempbis,temp2,temp2bis,temp3,temp3bis,meanL,meanR,v2LCP,v2RCP,x0,y0,x1,y1,RR1,RR2,i,loc0,loc1,loc2,loc3,loc4,xa,xb,ya,yb,indexL,indexR,indexL2,indexR2,poly,poly2,row,column,distance,j,minlookupt1,maxlookupt1,minlookupt2,maxlookupt2,FWHM,minimum,angulardistance,minimumR,minimumL,coeffs,u,t,k,yy,yy1,yy2,yy12,loc5,loc6,loc7,loc8,loc9,loc10,loc11,loc12,loc13,loc14,loc15)
 {

#pragma omp for
   for(iii=0;iii<nRows*nColumns;++iii)
	{
	  //with the convention adopted for a 2D array, the index iii is defined as iii=column+row*nColumns
	  row   =iii / nColumns;
	  column=iii % nColumns;

	  distance = sqrt(((float)row-Y0)*((float)row-Y0)+((float)column-X0)*((float)column-X0)); //distance in pixels
	  if(distance <= (Rsun+ExtraCrop))
		{

		  //CALCULATE THE SOLAR LINE PARAMETER AT THE DISTANCE FROM DISK CENTER (USING A LAW DERIVED FROM THE 3 Fe I PROFILES PROVIDED BY ROGER ULRICH)
		  //if(distance <= Rsun) angulardistance=cos(asin(distance/Rsun)); //convert projected angular distance from disk center to cos(angle) between the solar surface normal and the l.o.s. to the observer (INFINITE DISTANCE APPROXIMATION)
		  //else angulardistance=0.0;
		  //FWHM=FWHMCoeffs[0]+FWHMCoeffs[1]*angulardistance;
		  //FWHM=FWHM/2.0/sqrt(log(2.0))/1000.; //we convert from FWHM in mA to sigma in A
		  //minimum=minimumCoeffs[0]+minimumCoeffs[1]*angulardistance; //minimum intensity Id

		  //CALCULATE THE SOLAR LINE PARAMETER AT THE DISTANCE FROM DISK CENTER (USING A LAW DERIVED FROM A LINEWIDTH MAP OBTAINED WITH THE MDI-LIKE ALGORITHM)
		  distance=distance*cdelt1; //convert from pixels to arcsecs
		  FWHM=100.67102+0.015037016*distance-0.00010128197*distance*distance+3.1548385E-7*distance*distance*distance-3.7298102E-10*distance*distance*distance*distance+1.7275788E-13*distance*distance*distance*distance*distance;
		  FWHM=FWHM/2.0/sqrt(log(2.0))/1000.; //we convert from FWHM in mA to sigma in A

		  /*-------------------------------------------------------------*/
		  /* MDI-like algorithm                                          */
		  /*-------------------------------------------------------------*/

		  for(i=0;i<N;++i) 
		    {
		      tempvec=(float *)arrLev1p[i*2]->data;   //LCP (ASSUMES THE LEVEL 1p DATA ARE STORED LCP/RCP BACK-TO-BACK) IN THE ORDER I0, I1, I2, I3, I4, and I5
		      L[i]= tempvec[iii] ;
		      tempvec=(float *)arrLev1p[i*2+1]->data; //RCP
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
		      f1LCPc += cosi[i] *(double)L[i];
		      f1RCPc += cosi[i] *(double)R[i];
		      f1LCPs += sini[i] *(double)L[i];
		      f1RCPs += sini[i] *(double)R[i];
		      f2LCPc += cos2i[i]*(double)L[i];
		      f2RCPc += cos2i[i]*(double)R[i];
		      f2LCPs += sin2i[i]*(double)L[i];
		      f2RCPs += sin2i[i]*(double)R[i];
		    }
 		  
		  vLCP    = atan2(-f1LCPs,-f1LCPc)*pv1/2.0/M_PI; //-f1LCPs and -f1LCPc so that the jump is at 180 degrees and not 0 degrees
		  vRCP    = atan2(-f1RCPs,-f1RCPc)*pv1/2.0/M_PI;
      		  v2LCP   = atan2(-f2LCPs,-f2LCPc)*pv2/2.0/M_PI;
		  v2RCP   = atan2(-f2RCPs,-f2RCPc)*pv2/2.0/M_PI;
		  v2LCP   = fmod((v2LCP-vLCP+10.5*pv2),pv2)-pv2/2.0+vLCP; //we use the uncorrected velocity, i.e. phase, of the 1st Fourier coefficient to correct for the estimate of v2LCP and v2RCP, because the range of velocities obtained with the second Fourier coefficient is half the range of the first Fourier coefficient
		  v2RCP   = fmod((v2RCP-vRCP+10.5*pv2),pv2)-pv2/2.0+vRCP; 

		  if(isnan(f1LCPc) || isnan(f1RCPc) || isnan(f1LCPs) || isnan(f1RCPs))
		    {
		      MISSVALS2+=1;
		      lam0g[iii]  = MISSINGRESULT;
		      B0g[iii]    = MISSINGRESULT;
		      widthg[iii] = MISSINGRESULT;
		      Idg[iii]    = MISSINGRESULT;
		      I0g[iii]    = MISSINGRESULT;
		      continue;
		    }

		  /*-------------------------------------------------------------*/
		  /* bicubic interpolation of the look-up tables at pixel (x,y)  */
		  /*-------------------------------------------------------------*/
		  		  
		  //NB: it depends on how the filtergrams rebinning from 4096*4096 to axist[1]*axist[2] was done in phasemaps.c
		  //find the 4 neighbors (x0,y0), (x0,y1), (x1,y0), and (x1,y1) of (column,row) on the grid of the look-up tables, and deal with boundary problems

		  x0 = floor( ((float)column-offset)/(float)ratio); //a column number
		  y0 = floor( ((float)row-offset)/(float)ratio);    //a row number
		  x1 = x0+1;
		  y1 = y0+1;

		  //boundary issues
		  if(x1 >= axist[1]-1)
		    {
		      x0 = x0-2;
		      x1 = x1-2;
		    }
		  if(y1 >= axist[2]-1)
		    {
		      y0 = y0-2;
		      y1 = y1-2;
		    }
		  if(x1 < 2)
		    {
		      x0 = 1;
		      x1 = 2;
		    }
		  if(y1 < 2)
		    {
		      y0 = 1;
		      y1 = 2;
		    }

		  xb = (((float)column-offset)-(float)x0*(float)ratio)/(float)ratio;
		  yb = (((float)row-offset)-(float)y0*(float)ratio)/(float)ratio;

		  //for bicubic interpolation
		  loc0= x0   *axist[0]+ y0   *axist[0]*axist[1]; //x0
		  loc1= x1   *axist[0]+ y0   *axist[0]*axist[1]; //x1
		  loc2= x1   *axist[0]+ y1   *axist[0]*axist[1]; //x2
		  loc3= x0   *axist[0]+ y1   *axist[0]*axist[1]; //x3
		  loc4=(x0-1)*axist[0]+ y0   *axist[0]*axist[1]; //x4
		  loc5= x0   *axist[0]+(y0-1)*axist[0]*axist[1]; //x5
		  loc6= x1   *axist[0]+(y0-1)*axist[0]*axist[1]; //x6
		  loc7=(x1+1)*axist[0]+ y0   *axist[0]*axist[1]; //x7
		  loc8=(x1+1)*axist[0]+ y1   *axist[0]*axist[1]; //x8
		  loc9= x1   *axist[0]+(y1+1)*axist[0]*axist[1]; //x9
		  loc10= x0   *axist[0]+(y1+1)*axist[0]*axist[1]; //x10
		  loc11=(x0-1)*axist[0]+ y1   *axist[0]*axist[1]; //x11
		  loc12=(x0-1)*axist[0]+(y0-1)*axist[0]*axist[1]; //x12
		  loc13=(x1+1)*axist[0]+(y0-1)*axist[0]*axist[1]; //x13
		  loc14=(x1+1)*axist[0]+(y1+1)*axist[0]*axist[1]; //x14
		  loc15=(x0-1)*axist[0]+(y1+1)*axist[0]*axist[1]; //x15

		  indexL =ntest+1;
		  indexR =ntest+1;
		  indexL2=ntest+1;
		  indexR2=ntest+1;

		  //initializes arrays (NOT NEEDED ACTUALLY)
		  memset(poly,0.0,sizeof(poly));
		  memset(poly2,0.0,sizeof(poly2));


		  //minimum velocity
		  i=0;
		  yy[0]=lookupt[loc0+i];
		  yy[1]=lookupt[loc1+i];
		  yy[2]=lookupt[loc2+i];
		  yy[3]=lookupt[loc3+i];
		  yy1[0]=(lookupt[loc1+i]-lookupt[loc4+i])/2.0; //derivative in x direction
		  yy1[1]=(lookupt[loc7+i]-lookupt[loc0+i])/2.0;
		  yy1[2]=(lookupt[loc8+i]-lookupt[loc3+i])/2.0;
		  yy1[3]=(lookupt[loc2+i]-lookupt[loc11+i])/2.0;
		  yy2[0]=(lookupt[loc3+i] -lookupt[loc5+i])/2.0; //derivative in y direction
		  yy2[1]=(lookupt[loc2+i] -lookupt[loc6+i])/2.0;
		  yy2[2]=(lookupt[loc9+i] -lookupt[loc1+i])/2.0;
		  yy2[3]=(lookupt[loc10+i]-lookupt[loc0+i])/2.0;
		  yy12[0]=(lookupt[loc2+i] -lookupt[loc6+i] -lookupt[loc11+i]+lookupt[loc12+i])/4.0; //cross-derivatives
		  yy12[1]=(lookupt[loc8+i] -lookupt[loc13+i]-lookupt[loc3+i] +lookupt[loc5+i] )/4.0;
		  yy12[2]=(lookupt[loc14+i]-lookupt[loc7+i] -lookupt[loc10+i]+lookupt[loc0+i] )/4.0;
		  yy12[3]=(lookupt[loc9+i] -lookupt[loc1+i] -lookupt[loc15+i]+lookupt[loc4+i] )/4.0;
		  bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
		  u=yb;
		  t=xb;
		  poly[0]=0.0;
		  for (i = 3; i >= 0; i--) poly[0] = t*poly[0] + ((coeffs[i][3]*u+coeffs[i][2])*u + coeffs[i][1])*u + coeffs[i][0];

		  i=ntest;
		  yy[0]=lookupt[loc0+i];
		  yy[1]=lookupt[loc1+i];
		  yy[2]=lookupt[loc2+i];
		  yy[3]=lookupt[loc3+i];
		  yy1[0]=(lookupt[loc1+i]-lookupt[loc4+i])/2.0; //derivative in x direction
		  yy1[1]=(lookupt[loc7+i]-lookupt[loc0+i])/2.0;
		  yy1[2]=(lookupt[loc8+i]-lookupt[loc3+i])/2.0;
		  yy1[3]=(lookupt[loc2+i]-lookupt[loc11+i])/2.0;
		  yy2[0]=(lookupt[loc3+i] -lookupt[loc5+i])/2.0; //derivative in y direction
		  yy2[1]=(lookupt[loc2+i] -lookupt[loc6+i])/2.0;
		  yy2[2]=(lookupt[loc9+i] -lookupt[loc1+i])/2.0;
		  yy2[3]=(lookupt[loc10+i]-lookupt[loc0+i])/2.0;
		  yy12[0]=(lookupt[loc2+i] -lookupt[loc6+i] -lookupt[loc11+i]+lookupt[loc12+i])/4.0; //cross-derivatives
		  yy12[1]=(lookupt[loc8+i] -lookupt[loc13+i]-lookupt[loc3+i] +lookupt[loc5+i] )/4.0;
		  yy12[2]=(lookupt[loc14+i]-lookupt[loc7+i] -lookupt[loc10+i]+lookupt[loc0+i] )/4.0;
		  yy12[3]=(lookupt[loc9+i] -lookupt[loc1+i] -lookupt[loc15+i]+lookupt[loc4+i] )/4.0;
		  bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
		  u=yb;
		  t=xb;
		  poly2[0]=0.0;
		  for (i = 3; i >= 0; i--) poly2[0] = t*poly2[0] + ((coeffs[i][3]*u+coeffs[i][2])*u + coeffs[i][1])*u + coeffs[i][0];

		  minlookupt1 = poly[0];
		  minlookupt2 = poly2[0];


		  //maximum velocity
		  i=ntest-1;
		  yy[0]=lookupt[loc0+i];
		  yy[1]=lookupt[loc1+i];
		  yy[2]=lookupt[loc2+i];
		  yy[3]=lookupt[loc3+i];
		  yy1[0]=(lookupt[loc1+i]-lookupt[loc4+i])/2.0; //derivative in x direction
		  yy1[1]=(lookupt[loc7+i]-lookupt[loc0+i])/2.0;
		  yy1[2]=(lookupt[loc8+i]-lookupt[loc3+i])/2.0;
		  yy1[3]=(lookupt[loc2+i]-lookupt[loc11+i])/2.0;
		  yy2[0]=(lookupt[loc3+i] -lookupt[loc5+i])/2.0; //derivative in y direction
		  yy2[1]=(lookupt[loc2+i] -lookupt[loc6+i])/2.0;
		  yy2[2]=(lookupt[loc9+i] -lookupt[loc1+i])/2.0;
		  yy2[3]=(lookupt[loc10+i]-lookupt[loc0+i])/2.0;
		  yy12[0]=(lookupt[loc2+i] -lookupt[loc6+i] -lookupt[loc11+i]+lookupt[loc12+i])/4.0; //cross-derivatives
		  yy12[1]=(lookupt[loc8+i] -lookupt[loc13+i]-lookupt[loc3+i] +lookupt[loc5+i] )/4.0;
		  yy12[2]=(lookupt[loc14+i]-lookupt[loc7+i] -lookupt[loc10+i]+lookupt[loc0+i] )/4.0;
		  yy12[3]=(lookupt[loc9+i] -lookupt[loc1+i] -lookupt[loc15+i]+lookupt[loc4+i] )/4.0;
		  bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
		  u=yb;
		  t=xb;
		  poly[ntest-1]=0.0;
		  for (i = 3; i >= 0; i--) poly[ntest-1] = t*poly[ntest-1] + ((coeffs[i][3]*u+coeffs[i][2])*u + coeffs[i][1])*u + coeffs[i][0];

		  i=2*ntest-1;
		  yy[0]=lookupt[loc0+i];
		  yy[1]=lookupt[loc1+i];
		  yy[2]=lookupt[loc2+i];
		  yy[3]=lookupt[loc3+i];
		  yy1[0]=(lookupt[loc1+i]-lookupt[loc4+i])/2.0; //derivative in x direction
		  yy1[1]=(lookupt[loc7+i]-lookupt[loc0+i])/2.0;
		  yy1[2]=(lookupt[loc8+i]-lookupt[loc3+i])/2.0;
		  yy1[3]=(lookupt[loc2+i]-lookupt[loc11+i])/2.0;
		  yy2[0]=(lookupt[loc3+i] -lookupt[loc5+i])/2.0; //derivative in y direction
		  yy2[1]=(lookupt[loc2+i] -lookupt[loc6+i])/2.0;
		  yy2[2]=(lookupt[loc9+i] -lookupt[loc1+i])/2.0;
		  yy2[3]=(lookupt[loc10+i]-lookupt[loc0+i])/2.0;
		  yy12[0]=(lookupt[loc2+i] -lookupt[loc6+i] -lookupt[loc11+i]+lookupt[loc12+i])/4.0; //cross-derivatives
		  yy12[1]=(lookupt[loc8+i] -lookupt[loc13+i]-lookupt[loc3+i] +lookupt[loc5+i] )/4.0;
		  yy12[2]=(lookupt[loc14+i]-lookupt[loc7+i] -lookupt[loc10+i]+lookupt[loc0+i] )/4.0;
		  yy12[3]=(lookupt[loc9+i] -lookupt[loc1+i] -lookupt[loc15+i]+lookupt[loc4+i] )/4.0;
		  bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
		  u=yb;
		  t=xb;
		  poly2[ntest-1]=0.0;
		  for (i = 3; i >= 0; i--) poly2[ntest-1] = t*poly2[ntest-1] + ((coeffs[i][3]*u+coeffs[i][2])*u + coeffs[i][1])*u + coeffs[i][0];

		  maxlookupt1 = poly[ntest-1];
		  maxlookupt2 = poly2[ntest-1];

		  for(i=step;i<ntest;i=i+step) //make sure (ntest-1)/step is an integer
		    {
		      //poly[i]  =ya*(lookupt[loc1+i]*xa+lookupt[loc2+i]*xb)+yb*(lookupt[loc3+i]*xa+lookupt[loc4+i]*xb);                          //for 1st Fourier coefficient
		      //poly2[i] =ya*(lookupt[loc1+i+ntest]*xa+lookupt[loc2+i+ntest]*xb)+yb*(lookupt[loc3+i+ntest]*xa+lookupt[loc4+i+ntest]*xb);  //for 2nd Fourier coefficient

		      yy[0]=lookupt[loc0+i];
		      yy[1]=lookupt[loc1+i];
		      yy[2]=lookupt[loc2+i];
		      yy[3]=lookupt[loc3+i];
		      yy1[0]=(lookupt[loc1+i]-lookupt[loc4+i])/2.0; //derivative in x direction
		      yy1[1]=(lookupt[loc7+i]-lookupt[loc0+i])/2.0;
		      yy1[2]=(lookupt[loc8+i]-lookupt[loc3+i])/2.0;
		      yy1[3]=(lookupt[loc2+i]-lookupt[loc11+i])/2.0;
		      yy2[0]=(lookupt[loc3+i] -lookupt[loc5+i])/2.0; //derivative in y direction
		      yy2[1]=(lookupt[loc2+i] -lookupt[loc6+i])/2.0;
		      yy2[2]=(lookupt[loc9+i] -lookupt[loc1+i])/2.0;
		      yy2[3]=(lookupt[loc10+i]-lookupt[loc0+i])/2.0;
		      yy12[0]=(lookupt[loc2+i] -lookupt[loc6+i] -lookupt[loc11+i]+lookupt[loc12+i])/4.0; //cross-derivatives
		      yy12[1]=(lookupt[loc8+i] -lookupt[loc13+i]-lookupt[loc3+i] +lookupt[loc5+i] )/4.0;
		      yy12[2]=(lookupt[loc14+i]-lookupt[loc7+i] -lookupt[loc10+i]+lookupt[loc0+i] )/4.0;
		      yy12[3]=(lookupt[loc9+i] -lookupt[loc1+i] -lookupt[loc15+i]+lookupt[loc4+i] )/4.0;
		      bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
		      u=yb;
		      t=xb;
		      poly[i]=0.0;
		      for (j = 3; j >= 0; j--) poly[i] = t*poly[i] + ((coeffs[j][3]*u+coeffs[j][2])*u + coeffs[j][1])*u + coeffs[j][0];

		      yy[0]=lookupt[loc0+i+ntest];
		      yy[1]=lookupt[loc1+i+ntest];
		      yy[2]=lookupt[loc2+i+ntest];
		      yy[3]=lookupt[loc3+i+ntest];
		      yy1[0]=(lookupt[loc1+i+ntest]-lookupt[loc4+i+ntest])/2.0; //derivative in x direction
		      yy1[1]=(lookupt[loc7+i+ntest]-lookupt[loc0+i+ntest])/2.0;
		      yy1[2]=(lookupt[loc8+i+ntest]-lookupt[loc3+i+ntest])/2.0;
		      yy1[3]=(lookupt[loc2+i+ntest]-lookupt[loc11+i+ntest])/2.0;
		      yy2[0]=(lookupt[loc3+i+ntest] -lookupt[loc5+i+ntest])/2.0; //derivative in y direction
		      yy2[1]=(lookupt[loc2+i+ntest] -lookupt[loc6+i+ntest])/2.0;
		      yy2[2]=(lookupt[loc9+i+ntest] -lookupt[loc1+i+ntest])/2.0;
		      yy2[3]=(lookupt[loc10+i+ntest]-lookupt[loc0+i+ntest])/2.0;
		      yy12[0]=(lookupt[loc2+i+ntest] -lookupt[loc6+i+ntest] -lookupt[loc11+i+ntest]+lookupt[loc12+i+ntest])/4.0; //cross-derivatives
		      yy12[1]=(lookupt[loc8+i+ntest] -lookupt[loc13+i+ntest]-lookupt[loc3+i+ntest] +lookupt[loc5+i+ntest] )/4.0;
		      yy12[2]=(lookupt[loc14+i+ntest]-lookupt[loc7+i+ntest] -lookupt[loc10+i+ntest]+lookupt[loc0+i+ntest] )/4.0;
		      yy12[3]=(lookupt[loc9+i+ntest] -lookupt[loc1+i+ntest] -lookupt[loc15+i+ntest]+lookupt[loc4+i+ntest] )/4.0;
		      bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
		      u=yb;
		      t=xb;
		      poly2[i]=0.0;
		      for (j = 3; j >= 0; j--) poly2[i] = t*poly2[i] + ((coeffs[j][3]*u+coeffs[j][2])*u + coeffs[j][1])*u + coeffs[j][0];

		      if(poly[i] > vLCP && poly[i-step]   <= vLCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {

			      yy[0]=lookupt[loc0+j];
			      yy[1]=lookupt[loc1+j];
			      yy[2]=lookupt[loc2+j];
			      yy[3]=lookupt[loc3+j];
			      yy1[0]=(lookupt[loc1+j]-lookupt[loc4+j])/2.0; //derivative in x direction
			      yy1[1]=(lookupt[loc7+j]-lookupt[loc0+j])/2.0;
			      yy1[2]=(lookupt[loc8+j]-lookupt[loc3+j])/2.0;
			      yy1[3]=(lookupt[loc2+j]-lookupt[loc11+j])/2.0;
			      yy2[0]=(lookupt[loc3+j] -lookupt[loc5+j])/2.0; //derivative in y direction
			      yy2[1]=(lookupt[loc2+j] -lookupt[loc6+j])/2.0;
			      yy2[2]=(lookupt[loc9+j] -lookupt[loc1+j])/2.0;
			      yy2[3]=(lookupt[loc10+j]-lookupt[loc0+j])/2.0;
			      yy12[0]=(lookupt[loc2+j] -lookupt[loc6+j] -lookupt[loc11+j]+lookupt[loc12+j])/4.0; //cross-derivatives
			      yy12[1]=(lookupt[loc8+j] -lookupt[loc13+j]-lookupt[loc3+j] +lookupt[loc5+j] )/4.0;
			      yy12[2]=(lookupt[loc14+j]-lookupt[loc7+j] -lookupt[loc10+j]+lookupt[loc0+j] )/4.0;
			      yy12[3]=(lookupt[loc9+j] -lookupt[loc1+j] -lookupt[loc15+j]+lookupt[loc4+j] )/4.0;
			      bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
			      u=yb;
			      t=xb;
			      poly[j]=0.0;
			      for (k = 3; k >= 0; k--) poly[j] = t*poly[j] + ((coeffs[k][3]*u+coeffs[k][2])*u + coeffs[k][1])*u + coeffs[k][0];
			      
			      //for 1st Fourier coefficient
			      if(poly[j] > vLCP && poly[j-1]   <= vLCP) indexL = j-1;
			    }
			}
		      if(poly[i] > vRCP && poly[i-step]   <= vRCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {

			      yy[0]=lookupt[loc0+j];
			      yy[1]=lookupt[loc1+j];
			      yy[2]=lookupt[loc2+j];
			      yy[3]=lookupt[loc3+j];
			      yy1[0]=(lookupt[loc1+j]-lookupt[loc4+j])/2.0; //derivative in x direction
			      yy1[1]=(lookupt[loc7+j]-lookupt[loc0+j])/2.0;
			      yy1[2]=(lookupt[loc8+j]-lookupt[loc3+j])/2.0;
			      yy1[3]=(lookupt[loc2+j]-lookupt[loc11+j])/2.0;
			      yy2[0]=(lookupt[loc3+j] -lookupt[loc5+j])/2.0; //derivative in y direction
			      yy2[1]=(lookupt[loc2+j] -lookupt[loc6+j])/2.0;
			      yy2[2]=(lookupt[loc9+j] -lookupt[loc1+j])/2.0;
			      yy2[3]=(lookupt[loc10+j]-lookupt[loc0+j])/2.0;
			      yy12[0]=(lookupt[loc2+j] -lookupt[loc6+j] -lookupt[loc11+j]+lookupt[loc12+j])/4.0; //cross-derivatives
			      yy12[1]=(lookupt[loc8+j] -lookupt[loc13+j]-lookupt[loc3+j] +lookupt[loc5+j] )/4.0;
			      yy12[2]=(lookupt[loc14+j]-lookupt[loc7+j] -lookupt[loc10+j]+lookupt[loc0+j] )/4.0;
			      yy12[3]=(lookupt[loc9+j] -lookupt[loc1+j] -lookupt[loc15+j]+lookupt[loc4+j] )/4.0;
			      bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
			      u=yb;
			      t=xb;
			      poly[j]=0.0;
			      for (k = 3; k >= 0; k--) poly[j] = t*poly[j] + ((coeffs[k][3]*u+coeffs[k][2])*u + coeffs[k][1])*u + coeffs[k][0];

			      //for 1st Fourier coefficient
			      if(poly[j] > vRCP && poly[j-1]   <= vRCP) indexR = j-1;
			    }
			}
		      if(poly2[i]> v2LCP && poly2[i-step] <= v2LCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {

			      yy[0]=lookupt[loc0+j+ntest];
			      yy[1]=lookupt[loc1+j+ntest];
			      yy[2]=lookupt[loc2+j+ntest];
			      yy[3]=lookupt[loc3+j+ntest];
			      yy1[0]=(lookupt[loc1+j+ntest]-lookupt[loc4+j+ntest])/2.0; //derivative in x direction
			      yy1[1]=(lookupt[loc7+j+ntest]-lookupt[loc0+j+ntest])/2.0;
			      yy1[2]=(lookupt[loc8+j+ntest]-lookupt[loc3+j+ntest])/2.0;
			      yy1[3]=(lookupt[loc2+j+ntest]-lookupt[loc11+j+ntest])/2.0;
			      yy2[0]=(lookupt[loc3+j+ntest] -lookupt[loc5+j+ntest])/2.0; //derivative in y direction
			      yy2[1]=(lookupt[loc2+j+ntest] -lookupt[loc6+j+ntest])/2.0;
			      yy2[2]=(lookupt[loc9+j+ntest] -lookupt[loc1+j+ntest])/2.0;
			      yy2[3]=(lookupt[loc10+j+ntest]-lookupt[loc0+j+ntest])/2.0;
			      yy12[0]=(lookupt[loc2+j+ntest] -lookupt[loc6+j+ntest] -lookupt[loc11+j+ntest]+lookupt[loc12+j+ntest])/4.0; //cross-derivatives
			      yy12[1]=(lookupt[loc8+j+ntest] -lookupt[loc13+j+ntest]-lookupt[loc3+j+ntest] +lookupt[loc5+j+ntest] )/4.0;
			      yy12[2]=(lookupt[loc14+j+ntest]-lookupt[loc7+j+ntest] -lookupt[loc10+j+ntest]+lookupt[loc0+j+ntest] )/4.0;
			      yy12[3]=(lookupt[loc9+j+ntest] -lookupt[loc1+j+ntest] -lookupt[loc15+j+ntest]+lookupt[loc4+j+ntest] )/4.0;
			      bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
			      u=yb;
			      t=xb;
			      poly2[j]=0.0;
			      for (k = 3; k >= 0; k--) poly2[j] = t*poly2[j] + ((coeffs[k][3]*u+coeffs[k][2])*u + coeffs[k][1])*u + coeffs[k][0];

			      //for 2nd Fourier coefficient
			      if(poly2[j] > v2LCP && poly2[j-1]   <= v2LCP) indexL2 = j-1;
			    }
			}
		      if(poly2[i]> v2RCP && poly2[i-step] <= v2RCP)
			{
			  for(j=i-step+1;j<=i;j++)
			    {

			      yy[0]=lookupt[loc0+j+ntest];
			      yy[1]=lookupt[loc1+j+ntest];
			      yy[2]=lookupt[loc2+j+ntest];
			      yy[3]=lookupt[loc3+j+ntest];
			      yy1[0]=(lookupt[loc1+j+ntest]-lookupt[loc4+j+ntest])/2.0; //derivative in x direction
			      yy1[1]=(lookupt[loc7+j+ntest]-lookupt[loc0+j+ntest])/2.0;
			      yy1[2]=(lookupt[loc8+j+ntest]-lookupt[loc3+j+ntest])/2.0;
			      yy1[3]=(lookupt[loc2+j+ntest]-lookupt[loc11+j+ntest])/2.0;
			      yy2[0]=(lookupt[loc3+j+ntest] -lookupt[loc5+j+ntest])/2.0; //derivative in y direction
			      yy2[1]=(lookupt[loc2+j+ntest] -lookupt[loc6+j+ntest])/2.0;
			      yy2[2]=(lookupt[loc9+j+ntest] -lookupt[loc1+j+ntest])/2.0;
			      yy2[3]=(lookupt[loc10+j+ntest]-lookupt[loc0+j+ntest])/2.0;
			      yy12[0]=(lookupt[loc2+j+ntest] -lookupt[loc6+j+ntest] -lookupt[loc11+j+ntest]+lookupt[loc12+j+ntest])/4.0; //cross-derivatives
			      yy12[1]=(lookupt[loc8+j+ntest] -lookupt[loc13+j+ntest]-lookupt[loc3+j+ntest] +lookupt[loc5+j+ntest] )/4.0;
			      yy12[2]=(lookupt[loc14+j+ntest]-lookupt[loc7+j+ntest] -lookupt[loc10+j+ntest]+lookupt[loc0+j+ntest] )/4.0;
			      yy12[3]=(lookupt[loc9+j+ntest] -lookupt[loc1+j+ntest] -lookupt[loc15+j+ntest]+lookupt[loc4+j+ntest] )/4.0;
			      bcucof(yy,yy1,yy2,yy12,1.0,1.0,coeffs);
			      u=yb;
			      t=xb;
			      poly2[j]=0.0;
			      for (k = 3; k >= 0; k--) poly2[j] = t*poly2[j] + ((coeffs[k][3]*u+coeffs[k][2])*u + coeffs[k][1])*u + coeffs[k][0];

			      //for 2nd Fourier coefficient
			      if(poly2[j] > v2RCP && poly2[j-1]   <= v2RCP) indexR2 = j-1;
			    } 
			}
		    }


		  //TO DEAL WITH SATURATION
		  if(vLCP < minlookupt1)
		    { 
		      vLCP=vtest[0];
		      indexL=ntest+1;
		      //printf("take1 %f %f %f %f\n",vLCP,minlookupt1,vRCP,distance);
		    }
		  if(vLCP > maxlookupt1)
		    {
		      vLCP=vtest[ntest-1];
		      indexL=ntest+1;
		      //printf("take2 %f %f %f %f\n",vLCP,maxlookupt1,vRCP,distance);
		    }
		  if(vRCP < minlookupt1)
		    { 
		      vRCP=vtest[0];
		      indexR=ntest+1;
		      //printf("take3 %f %f %f %f\n",vRCP,minlookupt1,vLCP,distance);
		    }
		  if(vRCP > maxlookupt1)
		    {
		      vRCP=vtest[ntest-1];
		      indexR=ntest+1;
		      //printf("take4 %f %f %f %f\n",vRCP,maxlookupt1,vLCP,distance);
		    }
		  if(v2LCP < minlookupt2)
		    { 
		      v2LCP=vtest[0];
		      indexL2=ntest+1;
		    }
		  if(v2LCP > maxlookupt2)
		    {
		      v2LCP=vtest[ntest-1];
		      indexL2=ntest+1;
		    }
		  if(v2RCP < minlookupt2)
		    { 
		      v2RCP=vtest[0];
		      indexR2=ntest+1;
		    }
		  if(v2RCP > maxlookupt2)
		    {
		      v2RCP=vtest[ntest-1];
		      indexR2=ntest+1;
		    }

		  if(indexL == ntest+1 || indexR == ntest+1 || indexL2 == ntest+1 || indexR2 == ntest+1)
		    {
		      SATVALS2   += 1;
		    }
		  //else
		  //{

		  //We linearly interpolate in the look-up table for the 1st Fourier coefficient to retrieve the actual velocities
		  if(indexL != ntest+1)  vLCP   = vtest[indexL] +(vLCP-(double)poly[indexL])   *(vtest[indexL+1] -vtest[indexL]) /((double)poly[indexL+1]  -(double)poly[indexL]  );
		  if(indexR != ntest+1)  vRCP   = vtest[indexR] +(vRCP-(double)poly[indexR])   *(vtest[indexR+1] -vtest[indexR]) /((double)poly[indexR+1]  -(double)poly[indexR]  );
		  //We linearly interpolate in the look-up table for the 2nd Fourier coefficient to retrieve the actual velocities
		  if(indexL2 != ntest+1) v2LCP  = vtest[indexL2]+(v2LCP-(double)poly2[indexL2])*(vtest[indexL2+1]-vtest[indexL2])/((double)poly2[indexL2+1]-(double)poly2[indexL2]);
		  if(indexR2 != ntest+1) v2RCP  = vtest[indexR2]+(v2RCP-(double)poly2[indexR2])*(vtest[indexR2+1]-vtest[indexR2])/((double)poly2[indexR2+1]-(double)poly2[indexR2]);
		  
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
		  //lam0g[iii]  = (float)((vLCP+vRCP+v2LCP+v2RCP)/4.);//simple average. Need weights? REMINDER: SIGN CONVENTION: v<0 FOR MOTION TOWARD THE OBSERVER (BLUESHIFT)
		  lam0g[iii]  = (float)((vLCP+vRCP)/2.);
		  
		  //We compute the l.o.s. magnetic field
		  //B0g[iii]    = (float)((vLCP-vRCP+v2LCP-v2RCP)/2.0*magnetic);
		  B0g[iii]    = (float)((vLCP-vRCP)*magnetic);
		  
		  //We compute the linewidth (in Angstroms)
		  temp        = period/M_PI*sqrt(1.0/6.0*log((f1LCPc*f1LCPc+f1LCPs*f1LCPs)/(f2LCPc*f2LCPc+f2LCPs*f2LCPs)));
		  tempbis     = period/M_PI*sqrt(1.0/6.0*log((f1RCPc*f1RCPc+f1RCPs*f1RCPs)/(f2RCPc*f2RCPc+f2RCPs*f2RCPs)));
		  widthg[iii] = (float)((temp+tempbis)*sqrt(log(2.0)))*1000.; //we want the FWHM not the sigma of the Gaussian, in milliAngstoms	      
		  
		  //We compute the linedepth
		  //temp2       = period/2.0*sqrt(f1LCPc*f1LCPc+f1LCPs*f1LCPs)/sqrt(M_PI)/temp*exp(M_PI*M_PI*temp*temp/period/period);
		  //temp2bis    = period/2.0*sqrt(f1RCPc*f1RCPc+f1RCPs*f1RCPs)/sqrt(M_PI)/tempbis*exp(M_PI*M_PI*tempbis*tempbis/period/period);
		  //Idg[iii]    = (float)((temp2+temp2bis)/2.0);
		  
		  //We compute the linedepth (alternative formula)
		  //temp2       = period/2.0/sqrt(M_PI)/temp   *pow((f1LCPc*f1LCPc+f1LCPs*f1LCPs),2./3.)/pow((f2LCPc*f2LCPc+f2LCPs*f2LCPs),1./6.);
		  //temp2bis    = period/2.0/sqrt(M_PI)/tempbis*pow((f1RCPc*f1RCPc+f1RCPs*f1RCPs),2./3.)/pow((f2RCPc*f2RCPc+f2RCPs*f2RCPs),1./6.);
		  temp2       = period/2.0*sqrt(f1LCPc*f1LCPc+f1LCPs*f1LCPs)/sqrt(M_PI)/FWHM*exp(M_PI*M_PI*FWHM*FWHM/period/period);  //to not use the second Fourier coefficient
		  temp2bis    = period/2.0*sqrt(f1RCPc*f1RCPc+f1RCPs*f1RCPs)/sqrt(M_PI)/FWHM*exp(M_PI*M_PI*FWHM*FWHM/period/period);
		  Idg[iii]    = (float)((temp2+temp2bis)/2.0);
		  
		  //We compute the continuum intensity
		  //temp3       = (vLCP+v2LCP)/2.0/dv;
		  //temp3bis    = (vRCP+v2RCP)/2.0/dv;
		  //meanL=0.0;
		  //meanR=0.0;
		  //for(i=0;i<N;++i)
		  //{
		  //  meanL      += L[i];
		  //	  meanR      += R[i];
		  //}
		  //meanL=meanL/(double)N;
		  //meanR=meanR/(double)N;
		  //for(i=0;i<N;++i)
		  //{
		  //  meanL      += (temp2   /(double)N*exp(-(tune[i]-temp3)   *(tune[i]-temp3)   /temp   /temp));
		  //  meanR      += (temp2bis/(double)N*exp(-(tune[i]-temp3bis)*(tune[i]-temp3bis)/tempbis/tempbis));
		  //}
		  //
		  //I0g[iii]    = (float)((meanL+meanR)/2.0);
		  
		  temp3       = vLCP/dv;
		  temp3bis    = vRCP/dv;
		  meanL=0.0;
		  meanR=0.0;
		  for(i=0;i<N;++i)
		    {
		      meanL      += L[i];
		      meanR      += R[i];
		    }
		  meanL=meanL/(double)N;
		  meanR=meanR/(double)N;
		  //minimumL=(L[0]+L[N-1])/2.*minimum/(double)N; //(L[0]+L[N-1])/2. estimate of the continuum
		  //minimumR=(R[0]+R[N-1])/2.*minimum/(double)N;
		  for(i=0;i<N;++i)
		    {
		      //meanL      += (minimumL*exp(-(tune[i]-temp3)   *(tune[i]-temp3)   /FWHM/FWHM));
		      //meanR      += (minimumR*exp(-(tune[i]-temp3bis)*(tune[i]-temp3bis)/FWHM/FWHM));
		      meanL      += (temp2   /(double)N*exp(-(tune[i]-temp3)   *(tune[i]-temp3)   /FWHM/FWHM));
		      meanR      += (temp2bis/(double)N*exp(-(tune[i]-temp3bis)*(tune[i]-temp3bis)/FWHM/FWHM));
		    }
		  
		  I0g[iii]    = (float)((meanL+meanR)/2.0);
		  
		  // }		  
		  
		}//if (distance <= Rsun)
	  else
	    {
	      lam0g[iii]  = MISSINGRESULT;
	      B0g[iii]    = MISSINGRESULT;
	      widthg[iii] = MISSINGRESULT;
	      Idg[iii]    = MISSINGRESULT;
	      I0g[iii]    = MISSINGRESULT;
	    }
	  
	}//for iii
 }//end #pragma omp parallel
 
 
 free(cosi);
 free(sini);
 free(cos2i);
 free(sin2i);
 free(tune);
 free(angle);
 
 //fclose(fp);
 
 *MISSVALS=MISSVALS2;
 *SATVALS=SATVALS2;
 return error;
 
}//end routine
