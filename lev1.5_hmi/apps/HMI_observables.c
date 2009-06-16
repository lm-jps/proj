/*----------------------------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                                        */
/* OBSERVABLES MODULE FOR THE HMI PIPELINE OF STANFORD UNIVERSITY                                                                         */
/*                                                                                                                                        */
/* Authors:                                                                                                                               */
/* COUVIDAT, SCHOU, WACHTER                                                                                                               */
/* Version 1.0 April 10, 2009                                                                                                             */
/*                                                                                                                                        */
/* Assumption: in a cotune sequence, the filtergrams corresponding to the same wavelength                                                 */
/* are grouped together (example: I2LCP, I2RCP, I1LCP, I1RCP OK.; but I2LCP I1RCP I2RCP I1LCP Not OK!                                     */
/*                                                                                                                                        */
/* with a 6-position cotune sequence (assuming a separation of 68 mA):                                                                    */
/* I0 is centered at +170 mA                                                                                                              */
/* I1 is centered at +102. mA                                                                                                             */
/* I2 is centered at +34. mA                                                                                                              */
/* I3 is centered at -34. mA                                                                                                              */
/* I4 is centered at -102. mA                                                                                                             */
/* I5 is centered at -170 mA                                                                                                              */
/*                                                                                                                                        */
/* with a 5-position cotune sequence:                                                                                                     */
/* I0 is centered at +136 mA                                                                                                              */
/* I1 is centered at +68 mA                                                                                                               */
/* I2 is centered at  0 mA                                                                                                                */
/* I3 is centered at -68 mA                                                                                                               */
/* I4 is centered at -136 mA                                                                                                              */
/*                                                                                                                                        */
/* the level 1p data have 2 series: one with 24 segments named I0,Q0,U0,V0,..., one with 12 segments: LCP0,RCP0...                        */
/* the level 1.5 data have 5 series, each with 1 segment named dopplergram, magnetogram, linedepth, linewidth, continuum                  */
/*                                                                                                                                        */
/*                                                                                                                                        */
/* NAMING CONVENTIONS:                                                                                                                    */
/* LEVEL 1 FILTERGRAMS = flat-fielded, dark-subtracted filtergrams                                                                        */
/* LEVEL 1d FILTERGRAMS = Gapfilled, derotated, undistorted and temporally interpolated level 1                                           */
/*                        filtergrams (LEVEL 1d SERIES HAVE TWO PRIME KEYS: T_REC AND FID)                                                */
/* LEVEL 1p FILTERGRAMS = Polarization calibrated data IQUV and/or LCP+RCP                                                                */
/* LEVEL 1.5 DATA       = observables                                                                                                     */
/*                                                                                                                                        */
/*----------------------------------------------------------------------------------------------------------------------------------------*/

//Currently works on a limited amount of test data:
//for testing purposes: HMI_observables times="[2007.10.14_23:25:00-2007.10.14_23:25:10]" levin="lev1" levout="lev15" wavelength=2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>                      //OpenMP header
#include "interpol_code.h"            //from Richard
#include "polcal.h"                   //from Jesper
#include "HMIparam.h"

#undef I                              //I is the complex number (0,1). We un-define it to avoid confusion

char *module_name    = "HMI_observables"; //name of the module

#define kRecSetIn      "times"        //time range for which an output is wanted. MANDATORY PARAMETER.
                                      //the output will be UNIFORM IN EARTH TIME, NOT SDO TIME
#define kTypeSetIn     "levin"        //data level of the input filtergrams (lev1,lev1d,lev1p) LEV1 BY DEFAULT
#define kTypeSetOut    "levout"       //data level of the output series (lev1d,lev1p, or lev1.5) LEV1.5 BY DEFAULT
#define WaveLengthIn   "wavelength"   //filtergram Ii starting the framelist (i ranges from 0 to 5). MANDATORY PARAMETER.
#define PolarizationIn "polarization" //type of data returned by the polarization function (IQUV or LCP+RCP). LCP+RCP by default

                                      //arguments of the module
ModuleArgs_t module_args[] =        
{
     {ARG_STRING, kRecSetIn, ""       ,  "Time range (UT times) for which an output is wanted, in the format [2008.12.25_00:00:00-2008.12.25_01:00:00]"},
     {ARG_STRING, kTypeSetIn, "lev1.0",  "Level of input series: 1.0,1d,1p"},
     {ARG_STRING, kTypeSetOut,"lev1.5",  "Level of output series, combination of: 1d,1p,1.5. For example: 1p,1.5"},
     {ARG_INT   , WaveLengthIn,""     ,  "Index of the wavelength starting the framelist. FROM 0 TO 5"},
     {ARG_INT   , PolarizationIn, "3" ,  "Type of data returned by the polarization function: 1=IQUV; 2=LCP+RCP from 4 polarizations; 3=LCP+RCP from 2 polarizations. Must be 2 or 3 if lev1.5 is the desired output"},
     {ARG_END}
};

/*------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                  */
/*from Jesper: function that gives the sign of v                                                                    */
/*      Gives -1 for negative, 0 for 0 and 1 for positive                                                           */
/*                                                                                                                  */
/*------------------------------------------------------------------------------------------------------------------*/

int signj(int v)
{
return v > 0 ? 1 : (v < 0 ? -1 : 0);
}


/*------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                  */
/*function that uses the FID of a filtergram to tell what wavelength it is                                          */
/* returns -1 if there is an error                                                                                  */
/*                                                                                                                  */
/*------------------------------------------------------------------------------------------------------------------*/

int WhichWavelength(int FID)
{
  int result;

  switch (FID)
    {
    case 0: result=0; //LCP
      break;
    case 1: result=1;
      break;
    case 2: result=2;
      break;
    case 3: result=3;
      break;
    case 4: result=4;
      break;
    case 5: result=5;
      break;
    case 6: result=0;  //RCP
      break;
    case 7: result=1;
      break;
    case 8: result=2;
      break;
    case 9: result=3;
      break;
    case 10: result=4;
      break;
    case 11: result=5;
      break;
    default: result=-1;
    }

  return result;

}

/*------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                  */
/*function from R. Wachter                                                                                          */
/*                                                                                                                  */
/*------------------------------------------------------------------------------------------------------------------*/

void printtime()    // print time: use to debug code by identifying the part that take more time to run
{
  time_t timer, timerd;
  char *timestring;
  short i;

    timerd=time(&timer);
    timestring=ctime(&timer);
        for (i=0; i<24; ++i) printf("%c", *(timestring+i));
    printf("\n");
}


/*------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                  */
/*function to obtain data regarding the framelist that was used for a cotune sequence                               */
/*framelistInfo returns the number of filtergrams in the framelist (or 0 if a problem occured)                      */
/*PHWPLPOS is a pointer to an array of the HWLPOS and HPLPOS values of all the filtergrams of the framelist         */
/*WavelengthIndex is a pointer to an array of size framelistSize containing the order of the different wavelengths. */
/*For instance:                                                                                                     */
/*WavelengthIndex={0,0,2,2,4,4,1,1,5,5,3,3} means that the framelist is I0,I2,I4,I1,I5,I3 (with 2 polarizations)    */
/*format of PHWPLPOS:                                                                                               */
/*PHWPLPOS[i*7  ]=HWL1POS[filtergram index i of WavelengthIndex]                                                    */
/*PHWPLPOS[i*7+1]=HWL2POS[filtergram index i of WavelengthIndex]                                                    */
/*PHWPLPOS[i*7+2]=HWL3POS[filtergram index i of WavelengthIndex]                                                    */
/*PHWPLPOS[i*7+3]=HWL4POS[filtergram index i of WavelengthIndex]                                                    */
/*PHWPLPOS[i*7+4]=HPL1POS[filtergram index i of WavelengthIndex]                                                    */
/*PHWPLPOS[i*7+5]=HPL2POS[filtergram index i of WavelengthIndex]                                                    */
/*PHWPLPOS[i*7+6]=HPL3POS[filtergram index i of WavelengthIndex]                                                    */ 
/*                                                                                                                  */
/*------------------------------------------------------------------------------------------------------------------*/

int framelistInfo(int HFLID,int *PHWPLPOS,int *Wave)
{
  int framelistSize;

  if(HFLID == 1)
    {
      framelistSize=12;
      //PHWPLPOS=(int *)malloc(framelistSize*7*sizeof(int));
      //Wave=(int *)malloc(framelistSize*sizeof(int));

      //first filtergram (LCP)
      PHWPLPOS[0]=22;PHWPLPOS[1]=29;PHWPLPOS[2]=0;PHWPLPOS[3]=82;PHWPLPOS[4]=62;PHWPLPOS[5]=91;PHWPLPOS[6]=24;
      //second filtergram (RCP)
      PHWPLPOS[7]=22;PHWPLPOS[8]=29;PHWPLPOS[9]=0;PHWPLPOS[10]=82;PHWPLPOS[11]=62;PHWPLPOS[12]=91;PHWPLPOS[13]=54;
      //third filtergram (LCP)
      PHWPLPOS[14]=28;PHWPLPOS[15]=77;PHWPLPOS[16]=0;PHWPLPOS[17]=58;PHWPLPOS[18]=62;PHWPLPOS[19]=91;PHWPLPOS[20]=24;
      //4th filtergram (RCP)
      PHWPLPOS[21]=28;PHWPLPOS[22]=77;PHWPLPOS[23]=0;PHWPLPOS[24]=58;PHWPLPOS[25]=62;PHWPLPOS[26]=91;PHWPLPOS[27]=54;
      //5th filtergram
      PHWPLPOS[28]=34;PHWPLPOS[29]=65;PHWPLPOS[30]=0;PHWPLPOS[31]=94;PHWPLPOS[32]=62;PHWPLPOS[33]=91;PHWPLPOS[34]=24;
      //6th filtergram
      PHWPLPOS[35]=34;PHWPLPOS[36]=65;PHWPLPOS[37]=0;PHWPLPOS[38]=94;PHWPLPOS[39]=62;PHWPLPOS[40]=91;PHWPLPOS[41]=54;
      //7th filtergram
      PHWPLPOS[42]=40;PHWPLPOS[43]=53;PHWPLPOS[44]=0;PHWPLPOS[45]=70;PHWPLPOS[46]=62;PHWPLPOS[47]=91;PHWPLPOS[48]=24;
      //8th filtergram
      PHWPLPOS[49]=40;PHWPLPOS[50]=53;PHWPLPOS[51]=0;PHWPLPOS[52]=70;PHWPLPOS[53]=62;PHWPLPOS[54]=91;PHWPLPOS[55]=54;
      //9th filtergram
      PHWPLPOS[56]=46;PHWPLPOS[57]=41;PHWPLPOS[58]=0;PHWPLPOS[59]=106;PHWPLPOS[60]=62;PHWPLPOS[61]=91;PHWPLPOS[62]=24;
      //10th filtergram
      PHWPLPOS[63]=46;PHWPLPOS[64]=41;PHWPLPOS[65]=0;PHWPLPOS[66]=106;PHWPLPOS[67]=62;PHWPLPOS[68]=91;PHWPLPOS[69]=54;
      //11th filtergram
      PHWPLPOS[70]=52;PHWPLPOS[71]=29;PHWPLPOS[72]=0;PHWPLPOS[73]=82;PHWPLPOS[74]=62;PHWPLPOS[75]=91;PHWPLPOS[76]=24;
      //12th filtergram
      PHWPLPOS[77]=52;PHWPLPOS[78]=29;PHWPLPOS[79]=0;PHWPLPOS[80]=82;PHWPLPOS[81]=62;PHWPLPOS[82]=91;PHWPLPOS[83]=54;

      Wave[0]=0;Wave[1]=0;Wave[2]=1;Wave[3]=1;Wave[4]=2;Wave[5]=2;Wave[6]=3;Wave[7]=3;Wave[8]=4;Wave[9]=4;Wave[10]=5;Wave[11]=5;
      //0=I0,1=I1, etc...
    }
  else
    {
      printf("Error: HFLID should be 1\n");
      exit(EXIT_FAILURE);
    }

 return framelistSize;
 //framelistSize=0 if a problem occured
}



/*------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                  */
/*function that locates all the characters of a string s2 in a string s1                                            */
/*used to read the levout parameter                                                                                 */
/*returns 0 if all the characters of s2 are present in s1, 1 otherwise                                              */
/*                                                                                                                  */
/*------------------------------------------------------------------------------------------------------------------*/

int StringLocate(const char *s1,const char *s2)
{
  int i=0,result=0;
  char *ptr;
  char temp[2]="a\0";
  
  while(s2[i] != '\0')
    {
      temp[0]=s2[i];
      ptr = strpbrk(s1,temp);
      if(ptr == NULL) result=1;
      ++i;
    }

  return result;
}


/*------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                  */
/* Function that produces a mask for the gapfilling and temporal interpolation functions                            */
/*                                                                                                                  */
/*------------------------------------------------------------------------------------------------------------------*/


//for do_gapfill
int MaskCreation(unsigned char *Mask, int nColumn, int nRow)
{
  int status=0,i;
  int nElem;
  float row,column,distance;
  nElem=nRow*nColumn;
  float X0=2048.,Y0=2048.;

  for(i=0;i<nElem;++i)
    {
      
      row = (float)(i / nColumn);
      column = (float)(i % nColumn);
      distance=sqrt(  (row-Y0)*(row-Y0) + (column-X0)*(column-X0) ); //distance to Sun's center
      if(distance <= 1800.)
	{
	  if(row == 2047. || row == 2048. || column == 2047. || column == 2048.)
	    {
	      Mask[i] = 1;
	    }
	  else  Mask[i]=0;      //0 means pixel not missing
	}
      else Mask[i] = 2; //1 means pixel is missing and needs to be gap-filled
                        //2 means pixel is missing but does not need to be filled
      
    }


  return status;
}




/*---------------------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                                 */
/*                                                                                                                                 */
/*                                                                                                                                 */
/*                                                                                                                                 */
/*   DoIt is the entry point of the module                                                                                         */
/*                                                                                                                                 */
/*                                                                                                                                 */
/*                                                                                                                                 */
/*                                                                                                                                 */
/*---------------------------------------------------------------------------------------------------------------------------------*/

int DoIt(void)
{


  //Reading the command line parameters
  //*****************************************************************************************************************

  char *inRecQuery         = cmdparams_get_str(&cmdparams, kRecSetIn, NULL);         //time range
  char *inLev              = cmdparams_get_str(&cmdparams, kTypeSetIn, NULL);        //level of inout series
  char *outLev             = cmdparams_get_str(&cmdparams, kTypeSetOut, NULL);       //level of output series
  int   WavelengthID       = cmdparams_get_int(&cmdparams,WaveLengthIn , NULL);      //wavelength of the target filtergram
  int   PolarizationType   = cmdparams_get_int(&cmdparams,PolarizationIn, NULL);     //type of data returned by polarization function


  printf("COMMAND LINE= %s %s %s %d %d\n",inRecQuery,inLev,outLev,WavelengthID,PolarizationType);

  // Parameters                                                                                                    
  //*****************************************************************************************************************

#define  MaxNString 256                   //maximum length of strings in character number
  int   NumWavelengths=6;                 //number of possible values for the input WaveLengthID parameter
  int   MaxNumFiltergrams=36;             //maximum number of filtergrams in a framelist: 6 wavelengths*6 polarizations
  int   TempIntNum  = 6;                  //number of points requested for temporal interpolation (MUST BE EVEN !!!!!)
  TIME  TimeCaution = 50.;                //extra time in seconds padded to the beginning and ending input times
  const int   nRecmax     = 23040;        //maximum number of records that can be opened at once by the program (roughly 1 day of filtergrams: 23040=86400/3.75)
  char  HMISeriesLev1[MaxNString]  = "su_couvidat.HMISeriesLev1";    //name of the level 1 data series 
  char  HMISeriesLev1d[MaxNString] = "su_couvidat.HMISeriesLev1d";   //name of the level 1d data series
  char  HMISeriesLev1pa[MaxNString]= "su_couvidat.HMISeriesLev1pa";  //name of the level 1p data series FOR I+Q+U+V
  char  HMISeriesLev1pb[MaxNString]= "su_couvidat.HMISeriesLev1pb";  //name of the level 1p data series FOR LCP+RCP
  char  HMISeriesLev15a[MaxNString]= "su_couvidat.TestData";         //name of the level 1.5 data series FOR DOPPLERGRAM
  char  HMISeriesLev15b[MaxNString]= "su_couvidat.TestData2";        //name of the level 1.5 data series FOR MAGNETOGRAM
  char  HMISeriesLev15c[MaxNString]= "su_couvidat.TestData3";        //name of the level 1.5 data series FOR LINEDEPTH
  char  HMISeriesLev15d[MaxNString]= "su_couvidat.TestData4";        //name of the level 1.5 data series FOR LINEWIDTH
  char  HMISeriesLev15e[MaxNString]= "su_couvidat.TestData5";        //name of the level 1.5 data series FOR CONTINUUM
  char  HMISeriesLookup[MaxNString]= "su_couvidat.TestData7";        //name of the series containing the look-up tables for the MDI-like algorithm
  int   MaxSearchDistance0=4;             //maximum distance in framelist size at which we look for the filtergrams needed by the temporal interpolation code
  TIME  DopplergramCadence = 48.;         //Cadence in seconds (MUST MATCH THE T_REC_step OF LEV 1d, LEV 1p, and LEV 1.5 JSD FILES!!!)

  //Miscellaneous variables
  /******************************************************************************************************************/

  TIME  MaxSearchDistance;
  TIME  TREC_EPOCH = 0.;                   //Base epoch for T_REC keyword. Center of slot 0 for level 1d, 1p, and 1.5 data series. MUST BE THE SAME AS JSD FILES
  TIME  TREC_STEP = 0.;
  char  HMISeries[MaxNString];
  char  HMILookup[MaxNString];
  char  AcceptedLev[16][4] = {"0","d","p","5","dp","pd","d5","5d","dp5","d5p","pd5","p5d","5dp","5pd","p5","5p"};
  int   TestLevIn[3], TestLevOut[15];
  int   i,temp,TotalIn,TotalOut,observables,k,ii,iii;
  int   status = DRMS_SUCCESS, statusA[16], TotalStatus;
  DRMS_RecordSet_t *recLev1  = NULL;  //records for the level 1 data
  DRMS_RecordSet_t *recLev1d = NULL;  //records for the level 1d data
  DRMS_RecordSet_t *recLev1p = NULL;  //record for the level 1p data
  DRMS_RecordSet_t *recLev15a= NULL;  //record for the level 1.5 DOPPLERGRAM
  DRMS_RecordSet_t *recLev15b= NULL;  //record for the level 1.5 MAGNETOGRAM
  DRMS_RecordSet_t *recLev15c= NULL;  //record for the level 1.5 LINEDEPTH
  DRMS_RecordSet_t *recLev15d= NULL;  //record for the level 1.5 LINEWIDTH
  DRMS_RecordSet_t *recLev15e= NULL;  //record for the level 1.5 CONTINUUM
  DRMS_RecordSet_t *lookup    =NULL;  //record for the look-up tables for the MDI-like algorithm
  int   nRecs1,nRecs1d,nRecs1p,nRecs15;//number of records for the different types of data
  TIME  temptime=0.0;
  char  timeBegin[MaxNString] ="2000.12.25_00:00:00";
  char  timeEnd[MaxNString]   ="2000.12.25_00:00:00";
  char  timeBegin2[MaxNString]="2000.12.25_00:00:00";
  char  timeEnd2[MaxNString]  ="2000.12.25_00:00:00";
  TIME  TimeBegin,TimeEnd,TimeBegin2,TimeEnd2,TargetTime;
  int   ActualTempIntNum;                        //actual number of filtergrams used for temporal interpolation (if some are missing, ActualTempIntNum < TempIntNum)
  int   framelistSize=0;                         //size of the framelist for level 1 data
  TIME *internTOBS ;                             //array containing the time T_OBS of each filtergram opened, in seconds elapsed since a given date
  TIME  trec;
  TIME  tobs;
  int  *HWL1POS;                                 //Commanded wavelengths of each level 1 filtergram
  int  *HWL2POS;
  int  *HWL3POS;
  int  *HWL4POS;
  int  *HPL1POS;                                 //Commanded polarization of each level 1 filtergram
  int  *HPL2POS;
  int  *HPL3POS;
  int  *FID;                                     //filtergram ID
  int  *HFLID;                                   //framelist ID
  int  *FSN;
  float *RSUN;                                   //Radius of the Sun's image in pixels
  float *CROTA1;                                 //negative of solar P angle
  float *CRLTOBS;                                //solar B angle
  float *DSUNOBS;                                //Distance from Sun's center to SDO in meters
  float *X0;                                     //x-axis location of solar disk center in pixels 
  float *Y0;                                     //y-axis location of solar disk center in pixels 
  int  *SegmentRead;                             //Array that provides the status of the segments of the level 1 filtergrams
                                                 //SegmentRead[i]=0 if the segment of the filtergram i is not in memory
                                                 //SegmentRead[i]=1 if the segment of the filtergram i is in memory and will be used again
                                                 //SegmentRead[i]=-1 if the segment of the filtergram i is missing or corrupt                                       
  int  *Badkeyword;
  float RSUNAVG,X0AVG,Y0AVG,DSUNOBSAVG,CRLTOBSAVG,CROTA1AVG;
  char *TOBS = "T_OBS";                   //SHOULD I USE ANOTHER KEYWORD THAN T_OBS ??????????
  char *TuningPos1 = "HWL1POS";                  //tuning position of the 1st HCM
  char *TuningPos2 = "HWL2POS";                  //tuning position of the 2nd HCM
  char *TuningPos3 = "HWL3POS";                  //tuning position of the 3rd HCM
  char *TuningPos4 = "HWL4POS";                  //tuning position of the 4th HCM
  char *PolarizationPos1 = "HPL1POS";            //commanded polarization
  char *PolarizationPos2 = "HPL2POS";            //commanded polarization 
  char *PolarizationPos3 = "HPL3POS";            //commanded polarization 
  char *FiltergramID= "FID";                     //filtergram ID
  char *FramelistID= "HFLID";                    //framelist ID
  char *Rsun="R_SUN";
  char *Crota1="CROTA1";
  char *Crltobs="CRLT_OBS";
  char *Dsunobs="DSUN_OBS";
  char *x0="X0";
  char *y0="Y0";
  char *TREC="T_REC";                            //nominal time for the level 1d, 1p, and 1.5 data
  char *TRECEPOCH="T_REC_epoch";
  char *TRECSTEP="T_REC_step";
  DRMS_Array_t **Segments;                       //pointer to pointers to structures that will contain the segments of the level 1 filtergrams
  DRMS_Array_t **Ierror;                         //for gapfilling code
  int   TargetWavelength=0;                      //index of the filtergram level 1 with the wavelength WavelengthID and that is closest to TargetTime
  int  *IndexFiltergram;                         //an array that will contain the indeces of level 1 filtergrams with wavelength=WavelengthID 
  int  nIndexFiltergram;                         //size of array IndexFiltergram
  int   TargetHFLID=0;
  int   TargetHWLPOS[4];
  int   TargetHPLPOS[3];
  DRMS_Segment_t *segin  = NULL;
  int   axisin[2] ;                              //size of input filtergrams (axisin[0]=Ncolumns; axisin[1]=Nrows)
  DRMS_Array_t *arrin[TempIntNum];               //arrays that will contain pointers to the segments of the filtergrams needed for temporal interpolation
  DRMS_Array_t *arrerrors[TempIntNum];           //arrays that will contain pointers to the error maps returned by the gapfilling code
  int   axisout[2];                              //size of output filtergrams
  int   Nelem=0;                                 //total number of elements in a level 1 filtergram
  int   PHWPLPOS[MaxNumFiltergrams*7];
  int   WavelengthIndex[MaxNumFiltergrams], *OrganizeFramelist,  *OrganizeFramelist2, *FramelistArray, *SegmentStatus;
  int   FiltergramLocation, Needed;
  struct  initial const_param;                   //structure containing the parameters for Richard's functions
  unsigned char *Mask;                           //pointer to a 4096x4096 mask signaling which pixels are missing and which need to be filled
  float *image;                                  //for gapfilling code
  float *ierror;                                 //for gapfilling code
  float **images;                                //for temporal interpolation function
  float **ierrors;                               //for temporal interpolation function
  float **imagesout;                             //for polarization calibration function
  int   MaxDistance=0;
  DRMS_Array_t **arrLev1d= NULL;                 //pointer to pointer to an array that will contain a lev1d data produced by Richard's function
  DRMS_Array_t **arrLev1p= NULL;                 //pointer to pointer to an array that will contain a lev1p data produced by Jesper's function
  DRMS_Array_t **arrLev15= NULL;                 //pointer to pointer to an array that will contain a lev1.5 data produced by Seb's function
  DRMS_Array_t *TempSegment=NULL;
  DRMS_Array_t *arrintable= NULL;
  DRMS_Type_t type1d = DRMS_TYPE_FLOAT;          //type of the level 1d data produced by Richard's function
  DRMS_Type_t type1p = DRMS_TYPE_FLOAT;          //type of the level 1p data produced by Jesper's function
  DRMS_Type_t type15 = DRMS_TYPE_FLOAT;          //type of the level 1.5 data produced by Seb's function
  DRMS_Segment_t *segout = NULL;
  DRMS_Type_t typet = DRMS_TYPE_TIME;            //type of the keywords of type TIME!!!
  int Lev1pWanted=0;                             //do we need to produce and save level 1p data?
  int Lev1dWanted=0;                             //do we need to produce and save level 1d data?
  int Lev15Wanted=0;                             //do we need to produce and save level 1.5 data?
  int Segments1d=0;
  int Segments1p=0;
  struct keyword *KeyInterp=NULL;                //pointer to a list of structures containing some keywords needed by the temporal interpolation code
  struct keyword KeyInterpOut;
  struct parameterDoppler DopplerParameters;     //structure to provide some parameters defined in HMIparam.h to Dopplergram()
  int *ps1=NULL,*ps2=NULL,*ps3=NULL,*fid=NULL,*Wavelengths=NULL;
  struct polcal_struct pars;                     //for initialization of Jesper's routine
  int method=1;                                  //1 is the currently only implemented method WARNING: TO MODIFY
  float tsel;                                    //polarization selector temperature (in degree Celsius)
  float tfront;                                  //front window temperature (in degree Celsius)
  int npol;
  int npolout;
  int nSegs1p;
  int fsn;
  float SUNRADIUS, SUNCENTERX, SUNCENTERY;

  char Lev1pSegName[36][5]={"I0","Q0","U0","V0","I1","Q1","U1","V1","I2","Q2","U2","V2","I3","Q3","U3","V3","I4","Q4","U4","V4","I5","Q5","U5","V5","LCP0","RCP0","LCP1","RCP1","LCP2","RCP2","LCP3","RCP3","LCP4","RCP4","LCP5","RCP5"};            //names of the segments of the level 1 p records
                                                 //[0:23] are the segments for IQUV data, [24:35] are the segments for LCP+RCP data
  int Lev1pOffset;                               //this offset is 0 for IQUV segments, 24 for LCP+RCP  segments


  //Checking the values of the parameters inLev and outLev
  /******************************************************************************************************************/

  TotalIn=0;
  for(i=0;i<=2;++i)  
    {
      //temp = strcmp(inLev,AcceptedLev[i]);
      temp = StringLocate(inLev,AcceptedLev[i]); //the use of StringLocate allows a loose syntax for inLev and outLev (e.g. spaces can be present, commas can be forgotten...) 

      if(temp == 0) TestLevIn[i]=1; else  TestLevIn[i]=0;
      TotalIn+=TestLevIn[i];      
    }

  TotalOut=0;
  for(i=0;i<=14;++i) 
    {
      //temp = strcmp(outLev,AcceptedLev[i+1]); //i+1 because outLev cannot be "lev1"
      temp = StringLocate(outLev,AcceptedLev[i+1]);
      if(temp == 0) TestLevOut[i]=1; else TestLevOut[i]=0;
      TotalOut+=TestLevOut[i];
    }

  //Checking that the command line parameters are valide
  /******************************************************************************************************************/

  //checking that parameter for input filtergrams are valid
  if(TotalIn != 1)
    {
      if(StringLocate(inLev,"1") == 0) //special case: the user typed lev1 instead of 1.0
	{
	  TestLevIn[0]=1;
	  TotalIn=1;
	}
      else
	{
	  printf("The parameter levin must be one of the following strings (select only one): lev1, lev1d, or lev1p\n");
	  exit(EXIT_FAILURE);
	}
    }

  if(TotalOut == 0)
    {
      printf("The parameter levout must be a combination of the following strings: lev1d, lev1p, or lev1.5\n");
      exit(EXIT_FAILURE);
    }

  //Check for consistency: lowest output level>input level  
  if(TestLevIn[1]==1 && TestLevOut[13]==0 && TestLevOut[14]==0 && TestLevOut[1]==0 && TestLevOut[2]==0) //case levin=lev1d, lowest levout !=lev1p
    {
      printf("The parameter levin is the string lev1d, therefore levout can only be a combination of the strings lev1p and lev1.5\n");
      exit(EXIT_FAILURE);
    }
  if(TestLevIn[2]==1 && TestLevOut[2]==0)                      //case levin=lev1p, levout !=lev1.5
    {
      printf("The parameter levin is the string lev1p, therefore levout can only be the string lev1.5\n");
      exit(EXIT_FAILURE);
    }

  //check for valid values of PolarizationType
  if(PolarizationType !=1 && PolarizationType !=2 && PolarizationType !=3)
    {
      printf("The parameter polarization must be equal to 1, 2, or 3\n");
      exit(EXIT_FAILURE);
    }


  //Check for consistency: if outputs=level 1.5 AND level 1p data are not the input (i.e. they must be calculated) then PolarizationType must be 2 or 3 (LCP+RCP)
  if( (TestLevOut[0]!=1 && TestLevOut[1] !=1 && TestLevOut[3]!=1 && TestLevOut[4] !=1) && (TestLevIn[2] == 0) && (PolarizationType != 2 && PolarizationType != 3) )                      //case levin==lev1p, levout !=lev1.5
    {
      printf("The parameter levout includes the string lev1.5, but the level 1p data have be computed because the parameter levin is not lev1p. Therefore PolarizationType must be set to 1 ---LCP+RCP---.\n");
      exit(EXIT_FAILURE);
    }

  //check that the requested time string is in the correct format ([2008.12.25_00:00:00-2008.12.25_01:00:00])
  if(inRecQuery[0] != '[' || inRecQuery[40] != ']' || inRecQuery[20] != '-' || inRecQuery[5] != '.' || inRecQuery[8] != '.' || inRecQuery[25] != '.' || inRecQuery[28] != '.' || inRecQuery[11] != '_' || inRecQuery[31] != '_' || inRecQuery[14] != ':' || inRecQuery[17] != ':' || inRecQuery[34] != ':' || inRecQuery[37] != ':')
    {
      printf("The input parameter times does not have the correct format [2008.12.25_00:00:00-2008.12.25_01:00:00]\n");
      exit(EXIT_FAILURE);
    }

  //check that the command line parameter for the target filtergram is valid (must be in the range [0,5] for I0, I1, I2, I3, I4, or I5)
  if(WavelengthID > NumWavelengths-1 || WavelengthID < 0)
    {
      printf("The parameter filtergram is not in the range 0-5\n");
      exit(EXIT_FAILURE);  
    }


  //displaying on the screen the input and output data levels selected
  /******************************************************************************************************************/


  if(TestLevIn[0]==1) printf("Input data are level 1.0\n");
  if(TestLevIn[1]==1) printf("Input data are level 1d\n");
  if(TestLevIn[2]==1) printf("Input data are level 1p\n");
  if(TestLevOut[0]==1 || TestLevOut[3]==1 || TestLevOut[4]==1 || TestLevOut[5]==1 || TestLevOut[6]==1 || TestLevOut[7]==1 || TestLevOut[8]==1 || TestLevOut[9]==1 || TestLevOut[10]==1 || TestLevOut[11]==1 || TestLevOut[12]==1)
    {
      printf("Output data are level 1d\n");
      Lev1dWanted=1;
    }
  if(TestLevOut[1]==1 || TestLevOut[3]==1 || TestLevOut[4]==1 || TestLevOut[7]==1 || TestLevOut[8]==1 || TestLevOut[9]==1 || TestLevOut[10]==1 || TestLevOut[11]==1 || TestLevOut[12]==1 || TestLevOut[13]==1 || TestLevOut[14]==1)
    {  
      printf("Output data are level 1p\n");
      Lev1pWanted=1;
    }
  if(TestLevOut[2]==1 || TestLevOut[5]==1 || TestLevOut[6]==1 || TestLevOut[7]==1 || TestLevOut[8]==1 || TestLevOut[9]==1 || TestLevOut[10]==1 || TestLevOut[11]==1 || TestLevOut[12]==1 || TestLevOut[13]==1 || TestLevOut[14]==1)
    {
      printf("Output data are level 1.5\n");
      Lev15Wanted=1;
    }
  


  //converting the string containing the requested times into the DRMS TIME type data
  /******************************************************************************************************************/

  for(i=0;i<19;++i)
    {
      timeBegin[i]=inRecQuery[i+1];
      timeEnd[i]  =inRecQuery[i+21];
    }

  TimeBegin=sscan_time(timeBegin);                             //number of seconds elapsed since a given date
  TimeEnd  =sscan_time(timeEnd);   


  if(TimeBegin > TimeEnd)
    {
      printf("Error: the ending time must be after the beginning time!\n");
      exit(EXIT_FAILURE);
    }


  //initialization of Richard's and Jesper's codes
  //*************************************************************************************

  if(TestLevIn[0]==1)
    {
      status = initialize_interpol(&const_param);
      if(status != 0)
	{
	  printf("Error: could not initialize the gapfilling, derotation, and temporal interpolation routines\n");
	  exit(EXIT_FAILURE);
	}
      //printf("INITIALIZATION INTERPOL: %d %d %d %d %f %f %f\n",const_param.order_dist_coef,const_param.order2_rot_coef,const_param.order_int,const_param.nconst,const_param.diffrot_coef[0],const_param.diffrot_coef[1],const_param.diffrot_coef[2],); //CHECK THE FILLS STRUCTURE
    }
  if(Lev1pWanted || (Lev15Wanted && TestLevIn[2]==0))
    {
      status = init_polcal(&pars,method);
      if(status != 0)
	{
	  printf("Error: could not initialize the polarization calibration routine\n");
	  exit(EXIT_FAILURE);
	}
    }




  /***********************************************************************************************************************/
  /*                                                                                                                     */
  /*                                                                                                                     */
  /* IF INPUT IS LEVEL 1 FILTERGRAMS                                                                                     */
  /*                                                                                                                     */  
  /*                                                                                                                     */
  /***********************************************************************************************************************/


  if ( TestLevIn[0]==1 )                                       //input data are level 1 filtergrams (just flat-fielded)
    {

      //the requested time range [timeBegin,timeEnd] must be increased to take into account
      //the fact that the temporal interpolation scheme requires data points before and after
      //the first and last times wanted, and also we must add an extra time because of the
      //small difference there will be between SDO time and Earth time and because of the framelist length
      //*************************************************************************************

      TimeBegin2=TimeBegin-(TIME)TempIntNum*DopplergramCadence/2.-TimeCaution;
      TimeEnd2  =TimeEnd  +(TIME)TempIntNum*DopplergramCadence/2.+TimeCaution;
      sprint_ut(timeBegin2,TimeBegin2);                        //convert the time from TIME format to a string with UTC type
      sprint_ut(timeEnd2,TimeEnd2);
      strcat(HMISeriesLev1,"[]");                              //NB: I ASSUME THIS SERIES HAS AT LEAST 2 PRIMEKEYS, AND THAT THE TIME IS THE SECOND ONE, HENCE THIS FIRST [] !
      strcat(HMISeriesLev1,"[");                               //append the times at the end of the dataseries name for level 1 filtergrams
      strcat(HMISeriesLev1,timeBegin2);
      strcat(HMISeriesLev1,"-");
      strcat(HMISeriesLev1,timeEnd2);
      strcat(HMISeriesLev1,"]");                               //HMISeriesLev1 is in the format: seriesname[2000.12.25_00:00:00_UT-2000.12.25_00:00:00_UT]

      //opening the records in the range [TimeBegin2,TimeEnd2] and reading their keywords
      //*************************************************************************************

      recLev1 = drms_open_records(drms_env,HMISeriesLev1,&status); 
      
      if (status == DRMS_SUCCESS && recLev1 != NULL && recLev1->n > 0)  //successful opening of the input records (all these conditions are needed because the DRMS may claim it managed to open some records but the number of records might actually be 0. BUG?)
	{
	  nRecs1 = recLev1->n;                                     //number of records opened 

	  if(nRecs1 >= nRecmax)                                    //make sure this number of records does not exceed the maximum value allowed
	    {
	      printf("Too many records requested\n");
	      exit(EXIT_FAILURE);
	    }
	  
	  printf("Number of level 1 records opened: %d\n",nRecs1);

	  
	  FSN    = (int *)malloc(nRecs1*sizeof(int));          //Commanded wavelengths of each filtergram
	  if(FSN == NULL)
	    {
	      printf("Error: memory could not be allocated to FSN\n");
	      exit(EXIT_FAILURE);
	    }
	  internTOBS = (TIME *)malloc(nRecs1*sizeof(TIME));        //array containing the time T_OBS of each filtergram opened, in seconds elapsed since a given date
	  if(internTOBS == NULL)
	    {
	      printf("Error: memory could not be allocated to internTOBS\n");
	      exit(EXIT_FAILURE);
	    }
	  HWL1POS    = (int *)malloc(nRecs1*sizeof(int));          //Commanded wavelengths of each filtergram
	  if(HWL1POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HWL1POS\n");
	      exit(EXIT_FAILURE);
	    }
	  HWL2POS    = (int *)malloc(nRecs1*sizeof(int));
	  if(HWL2POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HWL2POS\n");
	      exit(EXIT_FAILURE);
	    }
	  HWL3POS    = (int *)malloc(nRecs1*sizeof(int));
	  if(HWL3POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HWL3POS\n");
	      exit(EXIT_FAILURE);
	    }
	  HWL4POS    = (int *)malloc(nRecs1*sizeof(int));
	  if(HWL4POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HWL4POS\n");
	      exit(EXIT_FAILURE);
	    }
	  HPL1POS    = (int *)malloc(nRecs1*sizeof(int));          //Commanded polarization of each filtergram
	  if(HPL1POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HPL1POS\n");
	      exit(EXIT_FAILURE);
	    }
	  HPL2POS    = (int *)malloc(nRecs1*sizeof(int));
	  if(HPL2POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HPL2POS\n");
	      exit(EXIT_FAILURE);
	    }
	  HPL3POS    = (int *)malloc(nRecs1*sizeof(int));
	  if(HPL3POS == NULL)
	    {
	      printf("Error: memory could not be allocated to HPL3POS\n");
	      exit(EXIT_FAILURE);
	    }
	  FID        = (int *)malloc(nRecs1*sizeof(int));          //filtergram ID
	  if(FID == NULL)
	    {
	      printf("Error: memory could not be allocated to FID\n");
	      exit(EXIT_FAILURE);
	    }
	  HFLID      = (int *)malloc(nRecs1*sizeof(int));          //framelist ID
	  if(HFLID == NULL)
	    {
	      printf("Error: memory could not be allocated to HFLID\n");
	      exit(EXIT_FAILURE);
	    }
	  RSUN       = (float *)malloc(nRecs1*sizeof(float));      //Radius of the Sun's image in pixels
	  if(RSUN == NULL)
	    {
	      printf("Error: memory could not be allocated to RSUN\n");
	      exit(EXIT_FAILURE);
	    }
	  CROTA1     = (float *)malloc(nRecs1*sizeof(float));      //negative of solar P angle
	  if(CROTA1 == NULL)
	    {
	      printf("Error: memory could not be allocated to CROTA1\n");
	      exit(EXIT_FAILURE);
	    }
	  CRLTOBS    = (float *)malloc(nRecs1*sizeof(float));      //solar B angle
	  if(CRLTOBS == NULL)
	    {
	      printf("Error: memory could not be allocated to CRLTOBS\n");
	      exit(EXIT_FAILURE);
	    }
	  DSUNOBS    = (float *)malloc(nRecs1*sizeof(float));      //Distance from Sun's center to SDO in meters
	  if(DSUNOBS == NULL)
	    {
	      printf("Error: memory could not be allocated to DSUNOBS\n");
	      exit(EXIT_FAILURE);
	    }
	  X0         = (float *)malloc(nRecs1*sizeof(float));      //x-axis location of solar disk center in pixels 
	  if(X0 == NULL)
	    {
	      printf("Error: memory could not be allocated to X0\n");
	      exit(EXIT_FAILURE);
	    }
	  Y0         = (float *)malloc(nRecs1*sizeof(float));      //y-axis location of solar disk center in pixels
 	  if(Y0 == NULL)
	    {
	      printf("Error: memory could not be allocated to Y0\n");
	      exit(EXIT_FAILURE);
	    }
	  SegmentRead= (int *)malloc(nRecs1*sizeof(int));
	  if(SegmentRead == NULL)
	    {
	      printf("Error: memory could not be allocated to SegmentRead\n");
	      exit(EXIT_FAILURE);
	    }
	  Badkeyword= (int *)malloc(nRecs1*sizeof(int));
	  if(Badkeyword == NULL)
	    {
	      printf("Error: memory could not be allocated to Badkeyword\n");
	      exit(EXIT_FAILURE);
	    }
	  Segments   = (DRMS_Array_t **)malloc(nRecs1*sizeof(DRMS_Array_t *));
	  if(Segments == NULL)
	    {
	      printf("Error: memory could not be allocated to Segments\n");
	      exit(EXIT_FAILURE);
	    }
	  Ierror   = (DRMS_Array_t **)malloc(nRecs1*sizeof(DRMS_Array_t *));
	  if(Ierror == NULL)
	    {
	      printf("Error: memory could not be allocated to Ierror\n");
	      exit(EXIT_FAILURE);
	    }
	  IndexFiltergram = (int *)malloc(nRecs1*sizeof(int));     //array that will contain the record index of the filtergrams with the same wavelength as WavelengthID
	  if(IndexFiltergram == NULL)
	    {
	      printf("Error: memory could not be allocated to IndexFiltergram\n");
	      exit(EXIT_FAILURE);
	    }

	  //reading some keyword values for all the open records (PUT -1 IF THE KEYWORD IS MISSING; MAKE SURE THE KEYWORDS ARE UNSIGNED INT !!!!) and
	  //create an array IndexFiltergram with the record index of all the filtergrams with the wavelength WavelengthID
	  //***********************************************************************************************************************

	  k=0;
	  for(i=0;i<nRecs1;++i) //loop over all the opened records
	    {	  
	      Badkeyword[i] = 0; //initialization: the keywords needed by do_interpolate are fine
	      //"Important" keywords: if one is missing, the filtergram is rejected
	      FSN[i] = drms_getkey_int(recLev1->records[i],"FSN",&statusA[0]); //not actually needed, just for debugging purpose
	      internTOBS[i] = drms_getkey_time(recLev1->records[i],TOBS,&statusA[0]);
	      HWL1POS[i]    = drms_getkey_int(recLev1->records[i],TuningPos1,&statusA[1]);
	      HWL2POS[i]    = drms_getkey_int(recLev1->records[i],TuningPos2,&statusA[2]); 
	      HWL3POS[i]    = drms_getkey_int(recLev1->records[i],TuningPos3,&statusA[3]); 
 	      HWL4POS[i]    = drms_getkey_int(recLev1->records[i],TuningPos4,&statusA[4]);
	      HPL1POS[i]    = drms_getkey_int(recLev1->records[i],PolarizationPos1,&statusA[5]);
	      HPL2POS[i]    = drms_getkey_int(recLev1->records[i],PolarizationPos2,&statusA[6]);
	      HPL3POS[i]    = drms_getkey_int(recLev1->records[i],PolarizationPos3,&statusA[7]);
	      FID[i]        = drms_getkey_int(recLev1->records[i],FiltergramID,&statusA[8]);
	      HFLID[i]      = drms_getkey_int(recLev1->records[i],FramelistID,&statusA[9]);
	      //The following are keywords needed by do_interpolate: if one is missing, the filtergram is not rejected but is tagged as having a problem
	      RSUN[i]       = drms_getkey_float(recLev1->records[i],Rsun,&statusA[10]);
	      if(statusA[10] != DRMS_SUCCESS)
		{
		  RSUN[i]=MISSINGKEYWORD;
		  Badkeyword[i]=1;
		}
	      CROTA1[i]     = drms_getkey_float(recLev1->records[i],Crota1,&statusA[11]);
	      if(statusA[11] != DRMS_SUCCESS)
		{
		  CROTA1[i]=MISSINGKEYWORD;
		  Badkeyword[i]=1;
		}
	      CRLTOBS[i]    = drms_getkey_float(recLev1->records[i],Crltobs,&statusA[12]);
	      if(statusA[12] != DRMS_SUCCESS)
		{
		  CRLTOBS[i]=MISSINGKEYWORD;
		  Badkeyword[i]=1;
		} 
	      DSUNOBS[i]    = drms_getkey_float(recLev1->records[i],Dsunobs,&statusA[13]);
	      if(statusA[13] != DRMS_SUCCESS)
		{
		  DSUNOBS[i]=MISSINGKEYWORD;
		  Badkeyword[i]=1;
		}
	      X0[i]         = drms_getkey_float(recLev1->records[i],x0,&statusA[14]);
	      if(statusA[14] != DRMS_SUCCESS)
		{
		  X0[i]=MISSINGKEYWORD;
		  Badkeyword[i]=1;
		}
	      Y0[i]         = drms_getkey_float(recLev1->records[i],y0,&statusA[15]);
	      if(statusA[15] != DRMS_SUCCESS)
		{
		  Y0[i]=MISSINGKEYWORD;
		  Badkeyword[i]=1;
		}
	      SegmentRead[i]= 0; //initialization: segment not in memory

	      //Now we test whether any important keyword is missing and we act accordingly
	      TotalStatus=0;
	      if(statusA[0] !=0)  //the time is not known
		{
		  internTOBS[i] = MISSINGKEYWORD;  //SHOULD I SET IT TO 0 INSTEAD? NANs ARE SLOW TO DELA WITH...
		  TotalStatus =1;
		}
	      for(ii=1;ii<=9;++ii) TotalStatus+=statusA[ii];
	      if(TotalStatus != 0)  //one of the "important" keyword is missing
		{
		  printf("Error: level 1 filtergram is missing at least one keyword\n");
		  //we set some keywords to unrealistic values so that this record cannot be considered a valid record later in the program and will be rejected
		  HWL1POS[i]    = MISSINGKEYWORDINT;
		  HWL2POS[i]    = MISSINGKEYWORDINT;
		  HWL3POS[i]    = MISSINGKEYWORDINT;
		  HWL4POS[i]    = MISSINGKEYWORDINT;
		  HPL1POS[i]    = MISSINGKEYWORDINT;
		  HPL2POS[i]    = MISSINGKEYWORDINT;
		  HPL3POS[i]    = MISSINGKEYWORDINT;
		}
	      else
		{
		  if(WhichWavelength(FID[i]) == WavelengthID) //function requestedWavelength returns 1 if FID of the filtergram i is one of those corresponding to WavelengthID
		    {
		      IndexFiltergram[k]=i;
		      ++k;
		    }
		}
	    }
	  nIndexFiltergram=k;
	  if(nIndexFiltergram == 0) //no filtergram was found with the target wavelength in the opened records
	    {
	      printf("Error: no filtergram was found with the wavelength %d in the requested level 1 records %s\n",WavelengthID,HMISeriesLev1);
	      exit(EXIT_FAILURE);
	    }	  
	  /*void *TempPoint;
	  TempPoint= realloc(IndexFiltergram,nIndexFiltergram*sizeof(int)); //to save some memory
	  if(TempPoint == NULL)
	    {
	      printf("Error: cannot reallocate memory for Indexfiltergram\n");
	      exit(EXIT_FAILURE);
	    }
	    else IndexFiltergram = (int *)TempPoint;*/


	} 
      else
	{
	  printf("Error: unable to open the requested level 1 records %s\n",HMISeriesLev1);
	  exit(EXIT_FAILURE);
	}    

    }



  /******************************************************************************************************************************************/
  /*                                                                                                                                        */
  /*                                                                                                                                        */
  /*                                                                                                                                        */
  /* LOOP OVER OBSERVABLE TIMES                                                                                                             */
  /*                                                                                                                                        */
  /*                                                                                                                                        */
  /*                                                                                                                                        */
  /******************************************************************************************************************************************/

  TargetTime = TimeBegin; //NEED TO ADD A FUNCTION TO SWITCH FROM SDO TIME TO EARTH TIME
  while(TargetTime <= TimeEnd)
    {
      

      sprint_ut(timeBegin2,TargetTime);                             //convert the time TargetTime from TIME format to a string with UTC type
      printf("\n Target time %s\n",timeBegin2);
      printf("-----------------------------------------------------------------------------------\n");


      /****************************************************************************************************************************/
      /*                                                                                                                          */
      /*                                                                                                                          */
      /* IF INPUT IS LEVEL 1d FILTERGRAMS                                                                                         */
      /*                                                                                                                          */
      /*                                                                                                                          */
      /****************************************************************************************************************************/
      
      if (TestLevIn[1]==1)                                          //input data are level 1d filtergrams (flat-fielded+derotated+un-distorted+temporally interpolated+gapfilled)
	{
	  strcpy(HMISeries,HMISeriesLev1d);
	//strcat(HMISeries,"[]");                                   //NB: I ASSUME THIS SERIES HAS AT LEAST 2 PRIMEKEYS, AND THAT THE TIME IS THE SECOND ONE, HENCE THIS FIRST [] !
	  strcat(HMISeries,"[");                                    //append the time at the end of the dataseries name for level 1d filtergrams
	  strcat(HMISeries,timeBegin2);
	  strcat(HMISeries,"]");                                    //HMISeriesLev1d is in the format: seriesname[2000.12.25_00:00:00_UT]
	  
	  recLev1d = drms_open_records(drms_env,HMISeries,&status); //open ALL the lev1d records whose T_OBS = target time (T_OBS is ASSUMED TO BE A PRIME KEY OF lev1d DATA)
	        
	  if (status == DRMS_SUCCESS && recLev1d != NULL && recLev1d->n > 0)           //successful opening of the input record
	    {
	      nRecs1d = recLev1d->n;                                //number of records in the level 1d series that have the same T_OBS value
	      
	      if(nRecs1d >= MaxNumFiltergrams)                      //If too many records were opened
		{
		  printf("Number of open record is larger than %d\n",MaxNumFiltergrams);
		  exit(EXIT_FAILURE);
		}
	      printf("Number of level 1d records opened= %d\n",nRecs1d);

	      trec = drms_getkey_time(recLev1d->records[0],TREC,&status); //T_REC, the slot time
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: unable to read the %s keyword of level 1d data\n",TREC);
		  exit(EXIT_FAILURE);
		}
	      tobs = drms_getkey_time(recLev1d->records[0],TOBS,&status); //T_REC, the nominal time
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: unable to read the %s keyword of level 1d data\n",TOBS);
		  exit(EXIT_FAILURE);
		}
	      TREC_STEP= drms_getkey_time(recLev1d->records[0],TRECSTEP,&status);
	      if(TREC_STEP != DopplergramCadence)
		{
		  printf("Error: the cadence is not equal to the T_REC_step keyword of the level 1d data\n");
		  exit(EXIT_FAILURE);
		}

	      DRMS_SegmentDimInfo_t di;
	      segin  = drms_segment_lookupnum(recLev1d->records[0],0);    //locating the first segment of the level 1p filtergram (either I0 or LCP0)
	      status = drms_segment_getdims(segin,&di);
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: unable to read the dimensions of the data segment of level 1d data\n");
		  exit(EXIT_FAILURE);
		}
	      axisin[0]= di.axis[0];                            //dimensions of the level 1p input data
	      axisin[1]= di.axis[1];
	      axisout[0]=axisin[0];                             //dimensions of the level 1.5 data
	      axisout[1]=axisin[1];

	      Segments1d=0;
	      Segments1p=0;
	    } 
	  else
	    {
	      printf("Unable to open the series %s for target time %s\n",HMISeries,timeBegin2);
	      goto NextTargetTime;
	    }
	} 





      /****************************************************************************************************************************/
      /*                                                                                                                          */
      /*                                                                                                                          */
      /* IF INPUT IS LEVEL 1 FILTERGRAMS                                                                                          */
      /*                                                                                                                          */
      /*                                                                                                                          */
      /****************************************************************************************************************************/

      if ( TestLevIn[0]==1 )                                        //input data are level 1 filtergrams (just flat-fielded)
	{

	  //We look for the filtergram with the wavelength WavelengthID that is closest to the target time
	  //*************************************************************************************

	  temptime = 100000000.0;                                   //in seconds; roughly 3.2 years
	  i=TargetWavelength;                                       //starting from the index of the filtergram with the target wavelength and closest to the previous target time
	  while(fabs(internTOBS[IndexFiltergram[i]]-TargetTime) <= temptime) //while the time difference decreases
	    {                                                                                        //when it increases, we know we reached the minimum
	      temptime=fabs(internTOBS[IndexFiltergram[i]]-TargetTime);
	      if(i <= nIndexFiltergram-2) ++i;
	      else break;
	    }
	  if(temptime >= DopplergramCadence )//if time difference between target time and time of the closest filtergram with the wavelength WavelengthID is larger than the cadence
	    {
	      goto NextTargetTime;           //should I really use a GOTO? It's kind of dangerous...
	    }
	  TargetWavelength= i-1;             //index of the filtergram with WavelengthID and that is closest to the target time: now called the TARGET FILTERGRAM
	  temp            = IndexFiltergram[TargetWavelength]; //index of the target filtergram
	  printf("TARGET FILTERGRAM = %d\n",FSN[temp]);
	  TargetHFLID     =   HFLID[temp];
	  TargetHWLPOS[0] = HWL1POS[temp];
	  TargetHWLPOS[1] = HWL2POS[temp];
	  TargetHWLPOS[2] = HWL3POS[temp];
	  TargetHWLPOS[3] = HWL4POS[temp];
	  TargetHPLPOS[0] = HPL1POS[temp];
	  TargetHPLPOS[1] = HPL2POS[temp];
	  TargetHPLPOS[2] = HPL3POS[temp];

	  if(SegmentRead[temp] != 1) //segment not already in memory
	    {
	      segin           = drms_segment_lookupnum(recLev1->records[temp],0);     //locating the first segment of the level 1 filtergram (SHOULD HAVE ONLY 1 SEGMENT)
	      Segments[temp]  = drms_segment_read(segin,type1d, &status); //reading the segment into memory (and converting it into type1d data)
	      if (status != DRMS_SUCCESS || Segments[temp] == NULL)
		{
		  Segments[temp]=NULL;
		  Ierror[temp]=NULL;
		  SegmentRead[temp]=-1; 
		  goto NextTargetTime;
		} 
	    }
	  axisin[0]= Segments[temp]->axis[0] ;                             //dimensions of the level 1 target filtergram
	  axisin[1]= Segments[temp]->axis[1] ;
	  axisout[0]=axisin[0];                                            //dimensions of the level 1d filtergram
	  axisout[1]=axisin[1];
	  Nelem=axisin[0]*axisin[1];
	  if(axisin[0] <= 0 || axisin[1] <= 0)
	    {
	      printf("Error: dimensions of segment of level 1 data record %d are 0x0\n",temp);
	      drms_free_array(Segments[temp]);
	      Segments[temp]=NULL;
	      SegmentRead[temp]=-1;                                        //-1 indicates a problem with the segment
	      goto NextTargetTime;
	    } 

	  //gapfilling of the filtergram just read
	  if(SegmentRead[temp] != 1)
	    {
	      Mask = (unsigned char *)malloc(Nelem*sizeof(unsigned char));
	      if(Mask == NULL)
		{
		  printf("Error: cannot allocate memory for Mask\n");
		  exit(EXIT_FAILURE);
		}
	      status = MaskCreation(Mask,axisin[0],axisin[1]) ;//first find the mask of missing pixels
	      if(status != 0)
		{
		  printf("Error: unable to create a mask for the gap filling function\n");
		  exit(EXIT_FAILURE);
		}
	      image  = Segments[temp]->data;
	      Ierror[temp] = drms_array_create(type1d,2,axisout,NULL,&status);
	      if(status != DRMS_SUCCESS || Ierror[temp] == NULL)
		{
		  drms_free_array(Segments[temp]);
		  Segments[temp]=NULL;
		  Ierror[temp]=NULL;
		  SegmentRead[temp]=-1; 
		  goto NextTargetTime;
		}     
	      ierror = Ierror[temp]->data;
	      status = do_gapfill(image,Mask,&const_param,ierror,axisin[0],axisin[1]); //then call the gapfilling function
	      if(status != 0)                               //gapfilling failed
		{
		  printf("Error: gapfilling code did not work on a level 1 filtergram at target time %s\n",timeBegin2);
		  SegmentRead[temp]=-1; 
		  drms_free_array(Segments[temp]);
		  drms_free_array(Ierror[temp]);
		  Segments[temp]=NULL;
		  Ierror[temp]=NULL;
		}
	      SegmentRead[temp]=1;
	      free(Mask);
	    }


	  //We call the function framelistInfo to obtain data regarding the framelist that was used for the cotune sequence
	  //to which the target filtergram belongs
	  //****************************************************************************************
	  framelistSize   = framelistInfo(TargetHFLID,PHWPLPOS,WavelengthIndex); //framelistSize is the number of filtergrams in the framelist
	                                                                //PHWPLPOS is a pointer to an array of the HWLPOS and HPLPOS values of all the filtergrams of the framelist
	                                                                //WavelengthIndex is a pointer to an array of integers of size framelistSize containing the order of the different wavelengths. 
	                                                                //For instance:
	                                                                //WavelengthIndex={0,0,2,2,4,4,1,1,5,5,3,3} means that the framelist is I0,I2,I4,I1,I5,I3 (with 2 different polarizations)
                                                            	        //format of PHWPLPOS: 
                                                            	        //PHWPLPOS[i*7  ]=HWL1POS[filtergram index i of WavelengthIndex]
                                                            	        //PHWPLPOS[i*7+1]=HWL2POS[filtergram index i of WavelengthIndex]
                                                            	        //PHWPLPOS[i*7+2]=HWL3POS[filtergram index i of WavelengthIndex]
                                                            	        //PHWPLPOS[i*7+3]=HWL4POS[filtergram index i of WavelengthIndex]
                                                            	        //PHWPLPOS[i*7+4]=HPL1POS[filtergram index i of WavelengthIndex]
                                                            	        //PHWPLPOS[i*7+5]=HPL2POS[filtergram index i of WavelengthIndex]
                                                            	        //PHWPLPOS[i*7+6]=HPL3POS[filtergram index i of WavelengthIndex] 
	                                                                //framelistInfo returns 0 if there is an error

	  if(framelistSize == 0)
	    {
	      printf("Error: cannot obtain information regarding the framelist for the level 1 data at target time %s\n",timeBegin2);
	      goto NextTargetTime;
	    }


	  //Creating DRMS records for the level 1d data
	  //**************************************************************************************************************

	  if (Lev1dWanted) recLev1d = drms_create_records(drms_env,framelistSize,HMISeriesLev1d,DRMS_PERMANENT,&status);  //record will be saved
	  else recLev1d = drms_create_records(drms_env,framelistSize,HMISeriesLev1d,DRMS_TRANSIENT,&status); //record will be discarded at the end of session
	  if(status != DRMS_SUCCESS || recLev1d == NULL || recLev1d->n < framelistSize)
	    {
	      printf("Error: cannot create records for the level 1d data %s at target time %s\n",HMISeriesLev1d,timeBegin2);
	      goto NextTargetTime;
	    }
	  nRecs1d= framelistSize; //number of level 1d data records to produce
	  Segments1d=0;           //segments for the level 1d data do not exist
	  Segments1p=0;           //segments for the level 1p data do not exist
	  TREC_EPOCH=drms_getkey_time(recLev1d->records[0],TRECEPOCH,&status); //CHANGE: I DON'T NEED TO DO THAT AT EACH ITERATION...
	  if(status != DRMS_SUCCESS)
	    {
	      printf("Error: unable to read the %s keyword for level 1d data\n",TRECEPOCH);
	      exit(EXIT_FAILURE);
	    }
	  TREC_STEP= drms_getkey_time(recLev1d->records[0],TRECSTEP,&status);
	  if(status != DRMS_SUCCESS)
	    {
	      printf("Error: unable to read the %s keyword for level 1d data\n",TRECSTEP);
	      exit(EXIT_FAILURE);
	    }
	  if(TREC_STEP != DopplergramCadence)
	    {
	      printf("Error: the cadence is not equal to the T_REC_step keyword of the level 1d data\n");
	      exit(EXIT_FAILURE);
	    }

	  tobs = TargetTime; //nominal time
	  trec = floor((tobs-TREC_EPOCH+TREC_STEP/2.0)/TREC_STEP)*TREC_STEP+TREC_EPOCH;  //slot time

	  //We decide how to group the filtergrams together to select the other target filtergrams (those close to TargetTime)
	  //**************************************************************************************************************
 
	  if( (framelistSize % 2) == 0) MaxDistance=framelistSize/2; //even number
	  else MaxDistance=(framelistSize-1)/2;                      //odd number
	  for(i=0;i<framelistSize;++i) if(WavelengthIndex[i] == WavelengthID && PHWPLPOS[i*7+4] == TargetHPLPOS[0] && PHWPLPOS[i*7+5] == TargetHPLPOS[1] && PHWPLPOS[i*7+6] == TargetHPLPOS[2] ) temp=i; //temp is the index of WavelengthIndex corresponding to the target filtergram
	  OrganizeFramelist =(int *)malloc(framelistSize*sizeof(int));//contains the index (relative to the target filtergram) of the filtergrams of the framelist
	  if(OrganizeFramelist == NULL)
	    {
	      printf("Error: memory could not be allocated to OrganizeFramelist\n");
	      exit(EXIT_FAILURE);
	    }
	  OrganizeFramelist2=(int *)malloc(framelistSize*sizeof(int));
	  if(OrganizeFramelist2 == NULL)
	    {
	      printf("Error: memory could not be allocated to OrganizeFramelist2\n");
	      exit(EXIT_FAILURE);
	    }

	  OrganizeFramelist[0]=-temp;
	  OrganizeFramelist2[0]=-1;
	  for(i=1;i<framelistSize;++i)
	    {
	      OrganizeFramelist[i]=i-temp;
	      if(WavelengthIndex[i] == WavelengthIndex[i-1]) OrganizeFramelist2[i]= OrganizeFramelist2[i-1]; //we want to make sure that same wavelengths are taken together
	      else OrganizeFramelist2[i]=signj(i-temp); //sign of i-temp
	    }

	  //creation of the array that will contain the location of the TempIntNum filtergrams of each type needed for the temporal interpolation algorithm
	  //format of FramelistArray:
	  //FramelistArray[i+framelistSize*k] = record index of the filtergram number k (out of the TempNumInt needed for the temporal interpolation)
	  //and of type WavelengthIndex[i] 
	  //***************************************************************************************************************************** 

	  MaxSearchDistance=(TIME)MaxSearchDistance0*DopplergramCadence;
	  FramelistArray = (int *)malloc(framelistSize*TempIntNum*sizeof(int));
	  if(FramelistArray == NULL)
	    {
	      printf("Error: memory could not be allocated to FramelistArray\n");
	      exit(EXIT_FAILURE);
	    }

	  for(i=0;i<framelistSize;++i)                                                                  //we loop over all filtergram types
	    {


	      //We fill the first framelistSize memory blocks of FramelistArray with IndexFiltergram[TargetWavelength]+OrganizeFramelist[i]
	      //*************************************************************************************

	      FiltergramLocation=IndexFiltergram[TargetWavelength]+OrganizeFramelist[i];                //location of the sought after filtergram in relation to the target filtergram
	      k=FiltergramLocation;      

	      FramelistArray[i]=k;
	      if(HWL1POS[k] != PHWPLPOS[i*7] || HWL2POS[k] != PHWPLPOS[i*7+1] || HWL3POS[k] !=  PHWPLPOS[i*7+2] || HWL4POS[k] != PHWPLPOS[i*7+3] || HPL1POS[k] != PHWPLPOS[i*7+4] || HPL2POS[k] != PHWPLPOS[i*7+5] || HPL3POS[k] !=  PHWPLPOS[i*7+6]) //there is an error somewhere: the filtergram we found is not what it is supposed to be!
		{ 
		  FramelistArray[i]=-1;                                                                 //-1 to indicate that the filtergram is missing
		}
	      

	      //We look for the TempIntNum-1 other filtergrams we need for the temporal interpolation
	      //*************************************************************************************

	      if(OrganizeFramelist2[i] < 0)   //the filtergram of index i in WavelengthIndex is located at or before the target filtergram
		{

		  //I need TempIntNum/2-1 filtergrams earlier than  FramelistArray[i], and TempIntNum/2 later
		  k=FiltergramLocation-1;
		  for(ii=1;ii<=TempIntNum/2-1;++ii)
		    {
		      while(HWL1POS[k] != PHWPLPOS[i*7] || HWL2POS[k] != PHWPLPOS[i*7+1] || HWL3POS[k] !=  PHWPLPOS[i*7+2] || HWL4POS[k] != PHWPLPOS[i*7+3] || HPL1POS[k] != PHWPLPOS[i*7+4] || HPL2POS[k] != PHWPLPOS[i*7+5] || HPL3POS[k] !=  PHWPLPOS[i*7+6])
			{
			  if ((internTOBS[k]-internTOBS[FiltergramLocation]) >= -MaxSearchDistance) --k; 
			  else break;
			}
		      if((internTOBS[k]-internTOBS[FiltergramLocation]) >= -MaxSearchDistance)
			{
			  FramelistArray[i+framelistSize*ii]=k;
			}
		      else  FramelistArray[i+framelistSize*ii]=-1; //missing filtergram
		      --k;
		    }


		  k=FiltergramLocation+1;
		  for(ii=TempIntNum/2;ii<TempIntNum;++ii)
		    {
		      while(HWL1POS[k] != PHWPLPOS[i*7] || HWL2POS[k] != PHWPLPOS[i*7+1] || HWL3POS[k] !=  PHWPLPOS[i*7+2] || HWL4POS[k] != PHWPLPOS[i*7+3] || HPL1POS[k] != PHWPLPOS[i*7+4] || HPL2POS[k] != PHWPLPOS[i*7+5] || HPL3POS[k] !=  PHWPLPOS[i*7+6])
			{
			  if ((internTOBS[k]-internTOBS[FiltergramLocation]) <= MaxSearchDistance)  ++k;
			  else break;
			}
		      if((internTOBS[k]-internTOBS[FiltergramLocation]) <= MaxSearchDistance) FramelistArray[i+framelistSize*ii]=k;
		      else FramelistArray[i+framelistSize*ii]=-1; //missing filtergram
		      ++k;
		    }
		  
		}
	      else   //in WavelengthIndex, the filtergram i is located after the target filtergram
		{


		  //I need TempIntNum/2 filtergrams earlier than  FramelistArray[i], and TempIntNum/2-1 later
		  k=FiltergramLocation-1;
		  for(ii=1;ii<=TempIntNum/2;++ii)
		    {
		      while(HWL1POS[k] != PHWPLPOS[i*7] || HWL2POS[k] != PHWPLPOS[i*7+1] || HWL3POS[k] !=  PHWPLPOS[i*7+2] || HWL4POS[k] != PHWPLPOS[i*7+3] || HPL1POS[k] != PHWPLPOS[i*7+4] || HPL2POS[k] != PHWPLPOS[i*7+5] || HPL3POS[k] !=  PHWPLPOS[i*7+6])
			{
			  if ((internTOBS[k]-internTOBS[FiltergramLocation]) >= -MaxSearchDistance) --k; 
			  else break;
			}
		      if((internTOBS[k]-internTOBS[FiltergramLocation]) >= -MaxSearchDistance) FramelistArray[i+framelistSize*ii]=k;
		      else FramelistArray[i+framelistSize*ii]=-1; //missing filtergram
		      --k;
		    }


		  k=FiltergramLocation+1;
		  for(ii=TempIntNum/2+1;ii<TempIntNum;++ii)
		    {
		      while(HWL1POS[k] != PHWPLPOS[i*7] || HWL2POS[k] != PHWPLPOS[i*7+1] || HWL3POS[k] !=  PHWPLPOS[i*7+2] || HWL4POS[k] != PHWPLPOS[i*7+3] || HPL1POS[k] != PHWPLPOS[i*7+4] || HPL2POS[k] != PHWPLPOS[i*7+5] || HPL3POS[k] !=  PHWPLPOS[i*7+6])
			{
			  if ((internTOBS[k]-internTOBS[FiltergramLocation]) <= MaxSearchDistance)  ++k;
			  else break;
			}
		      if((internTOBS[k]-internTOBS[FiltergramLocation]) <= MaxSearchDistance) FramelistArray[i+framelistSize*ii]=k;
		      else FramelistArray[i+framelistSize*ii]=-1; //missing filtergram
		      ++k;
		    }
		}
	      
	    }//end of for(i=0;i<framelistSize;++i)

	  free(OrganizeFramelist);
	  free(OrganizeFramelist2);
	  OrganizeFramelist=NULL;
	  OrganizeFramelist2=NULL;

	  //looking for filtergrams whose segment is already in memory but that are not needed anymore, and delete them
	  //**************************************************************************************************************


	  for(ii=0;ii<nRecs1;++ii) //loop over all the opened records (TRY TO FIND A FASTER WAY TO DO THAT...) 
	    {
	      if(SegmentRead[ii] == 1) //segment has already been read and is currently in memory
		{
		  Needed=0; //segment not needed a priori
		  for (k=0;k<framelistSize;++k)
		    {
		      for(i=0;i<TempIntNum;++i) if (FramelistArray[k+framelistSize*i] == ii) Needed=1; //Ah, my bad!!! The segment is actually needed
		    }
		  if(Needed == 0)//we delete the segment
		    {
		      drms_free_array(Segments[ii]);
		      drms_free_array(Ierror[ii]);
		      Segments[ii]=NULL;
		      Ierror[ii]=NULL;
		      SegmentRead[ii] = 0;
		    }

		}
	      
	    }
	  
	  /******************************************************************************************************************************/ 
	  /*for each type of filtergram                                                                                                 */
	  /* WE PRODUCE THE LEVEL 1D FILTERGRAMS                                                                                        */
	  /******************************************************************************************************************************/
	  
	  //creating the arrays that will contain the level 1d data
	  arrLev1d = (DRMS_Array_t **)malloc(nRecs1d*sizeof(DRMS_Array_t *));
	  if(arrLev1d == NULL)
	    {
	      printf("Error: memory could not be allocated to arrLev1d\n");
	      exit(EXIT_FAILURE);
	    }

	  for(k=0;k<nRecs1d;++k)
	    {
	      arrLev1d[k]= drms_array_create(type1d,2,axisout,NULL,&status);         
	      if(status != DRMS_SUCCESS || arrLev1d[k] == NULL)
		{
		  printf("Error: cannot create a DRMS array for a level 1d filtergram with index %d at target time %s\n",k,timeBegin2);
		  free(FramelistArray);
		  FramelistArray=NULL;
		  goto NextTargetTime;
		}
	    }
	  Segments1d=1; //segments for the level 1d data exist



	  //Compute the average values of R_SUN, X0, Y0, CRLT_OBS, CROTA1, and DSUNOBS, for the do_interpolate function
	  //**********************************************************************************************************
	  
	  RSUNAVG=0.0;
	  X0AVG=0.0;
	  Y0AVG=0.0;
	  DSUNOBSAVG=0.0;
	  CRLTOBSAVG=0.0;
	  CROTA1AVG=0.0;
		 
	  for(k=0;k<framelistSize;++k)
	    {
		  for(i=0;i<TempIntNum;++i) 
		    {
		      temp=FramelistArray[k+framelistSize*i];
		      if(temp != -1)  //if the filtergram is not missing
			{ 
			  if(RSUN[temp] != MISSINGKEYWORD && RSUN[temp] >= MINRSUN && RSUN[temp] <= MAXRSUN) RSUNAVG=RSUN[temp];
			  else Badkeyword[temp] = 1; //tag the filtergram as having a bad keyword
			  if(X0[temp] != MISSINGKEYWORD && X0[temp] >= XYMIN && X0[temp] <= XYMAX) X0AVG=X0[temp];
			  else Badkeyword[temp] = 1;
			  if(Y0[temp] != MISSINGKEYWORD && Y0[temp] >= XYMIN && Y0[temp] <= XYMAX) Y0AVG=Y0[temp];
			  else Badkeyword[temp] = 1;
			  if(DSUNOBS[temp] != MISSINGKEYWORD && DSUNOBS[temp] >= DSUNOBSMIN && DSUNOBS[temp] <= DSUNOBSMAX) DSUNOBSAVG=DSUNOBS[temp]/AstroUnit;
			  else Badkeyword[temp] = 1;
			  if(CRLTOBS[temp] != MISSINGKEYWORD && CRLTOBS[temp] >= CRLTOBSMIN && CRLTOBS[temp] <= CRLTOBSMAX) CRLTOBSAVG=CRLTOBS[temp];
			  else Badkeyword[temp] = 1;
			  if(CROTA1[temp] != MISSINGKEYWORD && CROTA1[temp] >= CROTA1MIN && CROTA1[temp] <= CROTA1MAX) CROTA1AVG=CROTA1[temp];
			  else Badkeyword[temp] = 1;
			}
		    }

	    }


	  for(k=0;k<framelistSize;++k)
	    {
	      
	      ActualTempIntNum=TempIntNum;

	      //Read the segments of the level 1 filtergrams needed and do the gapfilling
	      //***************************************************************************************************

	      for(i=0;i<TempIntNum;++i)
		{
		  temp=FramelistArray[k+framelistSize*i]; //index of the record
		  if(temp != -1)                          //if the filtergram is not missing 
		    {
		      if(Badkeyword[temp] != 1)           //if a keyword needed by do_interpolate is not missing
			{

			  if(SegmentRead[temp] == 0)      //if its segment is NOT already in memory and needs to be read
			    {
			      segin   = drms_segment_lookupnum(recLev1->records[temp], 0);
			      Segments[temp] = drms_segment_read(segin,type1d, &status); //pointer toward the segment (convert the data into type1d)
			      if (status != DRMS_SUCCESS || Segments[temp] == NULL)
				{
				  printf("Error: could not read the segment of level 1 record index %d at target time %s\n",temp,timeBegin2); //if there is a problem  
				  ActualTempIntNum-=1; //we will use one less filtergram for the temporal interpolation
				  arrin[i] = NULL;
				  arrerrors[i] = NULL;
				  Segments[temp] = NULL;
				  Ierror[temp] = NULL;
				  SegmentRead[temp]=-1;
				}  
			      else
				{
				  Ierror[temp] = drms_array_create(type1d,2,axisout,NULL,&status);
				  if(status != DRMS_SUCCESS || Ierror[temp] == NULL)
				    {
				      printf("Error: could not create an array for Ierror at target time %s\n",timeBegin2); //if there is a problem
				      drms_free_array(Segments[temp]);
				      Segments[temp]=NULL;
				      Ierror[temp]=NULL;
				      SegmentRead[temp]=-1; 
				      ActualTempIntNum-=1; //we will use one less filtergram for the temporal interpolation
				      arrin[i] = NULL;
				      arrerrors[i] = NULL;
				    }    
				  else
				    {
				      arrin[i] = Segments[temp];
				      arrerrors[i] = Ierror[temp];
				      if( arrin[i]->axis[0] != axisin[0]  || arrin[i]->axis[1] != axisin[1]) //segment does not have the same size as the segment of the target filtergram (PROBLEM HERE: I CURRENTLY DON'T CHECK IF A SEGMENT ALREADY IN MEMORY HAS THE SAME SIZE AS THE TARGET FILTERGRAM)
					{
					  printf("Error: level 1 record index %d at target time %s has a segment with dimensions %d x %d instead of %d x %d\n",temp,timeBegin2,arrin[i]->axis[0],arrin[i]->axis[1],axisin[0],axisin[1]);
					  drms_free_array(Segments[temp]);
					  drms_free_array(Ierror[temp]);
					  ActualTempIntNum-=1; //we will use one less filtergram for the temporal interpolation
					  arrin[i] = NULL;
					  arrerrors[i] = NULL;
					  Segments[temp]=NULL;
					  Ierror[temp]=NULL;
					  SegmentRead[temp]=-1;
					}
				      else
					{
					  SegmentRead[temp]=1; //now the segment for this record is in memory
					  //call gapfilling code of Richard and Jesper
					  //*******************************************************************
					  Mask = (unsigned char *)malloc(Nelem*sizeof(unsigned char));
					  if(Mask == NULL)
					    {
					      printf("Error: cannot allocate memory for Mask\n");
					      exit(EXIT_FAILURE);
					    }
					  status = MaskCreation(Mask,axisin[0],axisin[1]) ;//first find the mask of missing pixels
					  image  = arrin[i]->data;
					  ierror = arrerrors[i]->data;
					  status = do_gapfill(image,Mask,&const_param,ierror,axisin[0],axisin[1]); //then call the gapfilling function
					  free(Mask);
					  if(status != 0)                               //gapfilling failed
					    {
					      printf("Error: gapfilling code did not work on a level 1 filtergram at target time %s\n",timeBegin2);
					      ActualTempIntNum-=1; //we will use one less filtergram for the temporal interpolation
					      arrin[i] = NULL;     //IS THAT NECESSARY?
					      arrerrors[i] = NULL;
					      SegmentRead[temp]=-1; 
					      drms_free_array(Segments[temp]);
					      drms_free_array(Ierror[temp]);
					      Segments[temp]=NULL;
					      Ierror[temp]=NULL;
					    }
					}
				    }
				}
			    }//if(SegmentRead[temp] == 0)
			  
			  if(SegmentRead[temp] == 1)//if its segment is already in memory and is needed
			    {
			      arrin[i] = Segments[temp];
			      arrerrors[i] = Ierror[temp];
			    }
			  
			  if(SegmentRead[temp] == -1)//if its segment cannot be read (needed because of the TargetFiltergram that is read and gapfilled earlier
			    {
			      ActualTempIntNum-=1;
			      arrin[i] = NULL;
			      arrerrors[i] = NULL;
			    }
			}//if(Badkeyword != 1)
		      else //at least one of the keyword needed by do_interpolate is missing, so we just discard this filtergram so that it's not used by do_interpolate
			{
			  printf("Error: at least one of the keyword needed by the temporal interpolation function is missing, at target time %s\n",timeBegin2);
			  ActualTempIntNum-=1; //we will use one less filtergram for the temporal interpolation
			  arrin[i] = NULL;
			  arrerrors[i] = NULL;
			}
		    }//if(temp != -1)
		  else
		    {
		      ActualTempIntNum-=1; //the filtergram is missing: we will use one less filtergram for the temporal interpolation
		      arrin[i] = NULL;
		      arrerrors[i] = NULL;
		    }
		}//end for(i=0;i<=TempIntNum-1;++i)

	      //Do the temporal interpolation (Now the TempIntNum data segments needed for the temporal interpolation are in arrin)
	      //*******************************************************************************************************************


	      if(ActualTempIntNum >= 2) //if we have enough level 1 filtergrams to perform the temporal interpolation (IS 2 OK?)
		{

		  images = (float **)malloc(ActualTempIntNum*sizeof(float *));
		  if(images == NULL)
		    {
		      printf("Error: memory could not be allocated to images\n");
		      exit(EXIT_FAILURE);
		    }

		  ierrors = (float **)malloc(ActualTempIntNum*sizeof(float *));
		  if(ierrors == NULL)
		    {
		      printf("Error: memory could not be allocated to ierrors\n");
		      exit(EXIT_FAILURE);
		    }


		  KeyInterp = (struct keyword *)malloc(ActualTempIntNum*sizeof(struct keyword));
		  if(KeyInterp == NULL)
		    {
		      printf("Error: memory could not be allocated to KeyInterp\n");
		      exit(EXIT_FAILURE);
		    }
		  

		  //look for the available arrays
		  ii=0;
		  for(i=0;i<TempIntNum;++i) if (arrin[i] != NULL && arrerrors[i] != NULL) 
		    {
		      temp=FramelistArray[k+framelistSize*i];
		      printf("FSN filtergram used: %d\n",FSN[temp]);
		      images[ii]=arrin[i]->data;
		      ierrors[ii]=arrerrors[i]->data;
		      KeyInterp[ii].rsun=RSUN[temp];
		      KeyInterp[ii].xx0=X0[temp];
		      KeyInterp[ii].yy0=Y0[temp];
		      KeyInterp[ii].dist=DSUNOBS[temp]/AstroUnit; //Richard's code expect distance in AU (CHECK VALUE OF 1 AU)
		      KeyInterp[ii].b0=CRLTOBS[temp];
		      KeyInterp[ii].p0=CROTA1[temp];
		      KeyInterp[ii].time=internTOBS[temp];
		      ii+=1;
		    } //ii should be equal to ActualTempIntNum
		  

		  //provide the target values of some keywords for the temporal interpolation
		  KeyInterpOut.rsun=RSUNAVG;
		  KeyInterpOut.xx0=X0AVG;
		  KeyInterpOut.yy0=Y0AVG;
		  KeyInterpOut.dist=DSUNOBSAVG/AstroUnit;
		  KeyInterpOut.b0=CRLTOBSAVG;
		  KeyInterpOut.p0=CROTA1AVG;
		  KeyInterpOut.time=TargetTime;
		  
		  //calling Richard's code: temporal interpolation, de-rotation, un-distortion
		  status=do_interpolate(images,ierrors,arrLev1d[k]->data,KeyInterp,&KeyInterpOut,&const_param,ActualTempIntNum,axisin[0],axisin[1]);
		  if (status != 0)
		    {
		      printf("Error: temporal interpolation failed at target time %s\n",timeBegin2);
		      drms_free_array(arrLev1d[k]);
		      arrLev1d[k] = NULL;		      
		    }
		  else
		    {
		      for(i=0;i<TempIntNum;++i)
			{
			  temp=FramelistArray[k+framelistSize*i];
			  if(temp != -1) break;
			}

		      drms_copykeys(recLev1d->records[k],recLev1->records[temp],1, kDRMS_KeyClass_Explicit); //we copy all the keywords from the level 1 data to this new level 1d record			    		      
		      statusA[0] = drms_setkey_float(recLev1d->records[k],"R_SUN",KeyInterpOut.rsun); 
		      statusA[1] = drms_setkey_float(recLev1d->records[k],"X0",KeyInterpOut.xx0);
		      statusA[2] = drms_setkey_float(recLev1d->records[k],"Y0",KeyInterpOut.yy0);
		      statusA[3] = drms_setkey_float(recLev1d->records[k],"DSUN_OBS",KeyInterpOut.dist*AstroUnit);
		      statusA[4] = drms_setkey_float(recLev1d->records[k],"CRLT_OBS",KeyInterpOut.b0);
		      statusA[5] = drms_setkey_float(recLev1d->records[k],"CROTA1",KeyInterpOut.p0);
		      statusA[6] = drms_setkey_time(recLev1d->records[k],TOBS,KeyInterpOut.time); //KeyInterpOut.time should be equal to TargetTime and is the nominal time
		      statusA[7] = drms_setkey_time(recLev1d->records[k],TREC,trec); //TREC is the slot time

		      TotalStatus=0;
		      for(i=0;i<8;++i) TotalStatus+=statusA[i];
		      if(TotalStatus != 0)
			{
			  printf("Error: could not set some keyword for the level 1d data\n");
			  exit(EXIT_FAILURE);
			}

		      //if needed, write the segment
		      if (Lev1dWanted)
			{
			  segout = drms_segment_lookupnum(recLev1d->records[k], 0);
			  drms_segment_write(segout,arrLev1d[k],0); //write the file containing the data
			}
		      
		    }

		  free(images);
		  free(ierrors);
		  images=NULL;
		  ierrors=NULL;
		  free(KeyInterp);
		  KeyInterp=NULL;
		}//if(ActualTempIntNum >= 2)
	      else
		{
		  printf("Error: not enough valid level 1 filtergrams to produce a level 1d filtergram at target time %s\n",timeBegin2);
		  drms_free_array(arrLev1d[k]);
		  arrLev1d[k]=NULL;
		  //WHAT ELSE TO DO????
		}
	      
	      
	    }//end of for(k=0;k<framelistSize;++k)
      
      
	  free(FramelistArray);
	  FramelistArray=NULL;
	}//end of if(TestLevIn[0] == 1)
      
  




      /****************************************************************************************************************************/
      /*                                                                                                                          */
      /*                                                                                                                          */
      /* IF INPUT IS LEVEL 1P FILTERGRAMS                                                                                         */
      /*                                                                                                                          */
      /*                                                                                                                          */
      /****************************************************************************************************************************/
      

      if (TestLevIn[2]==1) 
	{
	  if(PolarizationType ==3 || PolarizationType ==2) strcpy(HMISeries,HMISeriesLev1pb); //LCP+RCP
	  else strcpy(HMISeries,HMISeriesLev1pa);                   //I+Q+U+V
	//strcat(HMISeries,"[]");                                   //NB: I ASSUME THIS SERIES HAS AT LEAST 2 PRIMEKEYS, AND THAT THE TIME IS THE SECOND ONE, HENCE THIS FIRST [] !
	  strcat(HMISeries,"[");                                    //append the time at the end of the dataseries name for level 1d filtergrams
	  strcat(HMISeries,timeBegin2);
	  strcat(HMISeries,"]");                                    //HMISeriesLev1d is in the format: seriesname[2000.12.25_00:00:00_UT]
	  
	  recLev1p = drms_open_records(drms_env,HMISeries,&status);
	  
	  if (status == DRMS_SUCCESS && recLev1p != NULL && recLev1p->n > 0)  //successful opening of the input records
	    {	      
	      if(recLev1p->n > 1)                            //make sure this number of records does not exceed the maximum value allowed
		{
		  printf("Too many records: %d\n",recLev1p->n);
		  exit(EXIT_FAILURE);
		}
	      printf("Number of level 1p records opened= %d\n",recLev1p->n);

	      nSegs1p = drms_record_numsegments(recLev1p->records[0]); //SHOULD HAVE A KEYWORD IN THE LEV1P SERIES INSTEAD !!!
	      printf("NUMBER OF SEGMENTS %d\n",nSegs1p);


	      if (PolarizationType ==1)
		{
		  nRecs1p = nSegs1p/4;               //5 or 6 wavelengths, 4 polarizations (I+Q+U+V)

		}
	      if (PolarizationType ==2 || PolarizationType ==3)
		{
		  nRecs1p = nSegs1p/2;  //5 or 6 wavelengths, 2 polarizations (LCP+RCP)
		}

	      trec = drms_getkey_time(recLev1p->records[0],TREC,&status);
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: cannot read the keyword %s\n",TREC);
		  exit(EXIT_FAILURE);
		}
	      tobs = drms_getkey_time(recLev1p->records[0],TOBS,&status);
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: cannot read the keyword %s\n",TOBS);
		  exit(EXIT_FAILURE);
		}
	      TREC_STEP= drms_getkey_time(recLev1p->records[0],TRECSTEP,&status);
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: cannot read the keyword %s\n",TRECSTEP);
		  exit(EXIT_FAILURE);
		}

	      if(TREC_STEP != DopplergramCadence)
		{
		  printf("Error: the cadence is not equal to the T_REC_step keyword of the level 1p data\n");
		  exit(EXIT_FAILURE);
		}

	      DRMS_SegmentDimInfo_t di;
	      segin  = drms_segment_lookupnum(recLev1p->records[0],0);     //locating the first segment of the level 1p filtergram (either I0 or LCP0)
	      status = drms_segment_getdims(segin,&di);
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: cannot read the dimensions of the data segment of level 1p data\n");
		  exit(EXIT_FAILURE);
		}
	      axisin[0]= di.axis[0];                            //dimensions of the level 1p input data
	      axisin[1]= di.axis[1];
	      axisout[0]=axisin[0];                             //dimensions of the level 1.5 data
	      axisout[1]=axisin[1];

	      Segments1p=0;//segments for level 1p data not read

	    } 
	  else
	    {
	      printf("Unable to open the series %s for time %s\n",HMISeries,timeBegin2);
	      goto NextTargetTime;
	    }
	}
      





      /****************************************************************************************************************************/
      /*                                                                                                                          */
      /*                                                                                                                          */
      /* IF REQUESTED OUTPUT INCLUDES LEVEL 1P AND/OR LEVEL 1.5 DATA                                                              */
      /* WE PRODUCE LEVEL 1P DATA                                                                                                 */
      /*                                                                                                                          */
      /*                                                                                                                          */
      /****************************************************************************************************************************/


      if (Lev1pWanted == 1 || (Lev15Wanted == 1 && TestLevIn[2]!=1)) 
	{
	  // Given consistency tests and the preceeding code we know that lev1d exists
	  //the level 1d records are the recLev1d[k], with the number of k = nRecs1d
	  
	  int ActualnSegs1d=0;       //number of level 1d segments that are not a NULL pointer
	  
	  //reading the level 1d segments 
	  //******************************************************************************
	  if(Segments1d == 0)        //segments are not in memory, then we populate the arrLev1d DRMS arrays
	    {
	      ActualnSegs1d=nRecs1d; 
	      
	      arrLev1d = (DRMS_Array_t **)malloc(nRecs1d*sizeof(DRMS_Array_t *));
	      if(arrLev1d == NULL)
		{
		  printf("Error: memory could not be allocated to arrLev1d\n");
		  exit(EXIT_FAILURE);
		}

	      for(i=0;i<nRecs1d;++i)
		{
		  arrLev1d[i]= drms_array_create(type1d,2,axisout,NULL,&status);         
		  if(status != DRMS_SUCCESS || arrLev1d[i] == NULL)
		    {
		      printf("Error: cannot create a DRMS array for a level 1d filtergram with index %d at target time %s\n",k,timeBegin2);
		      goto NextTargetTime;
		    }	      
		  if(recLev1d->records[i] != NULL)
		    {
		      segin   = drms_segment_lookupnum(recLev1d->records[i], 0);
		      arrLev1d[i] = drms_segment_read(segin, segin->info->type, &status); //pointer toward the segment
		      if(status != DRMS_SUCCESS || arrLev1d[i] == NULL)
			{
			  printf("Error: could not read the segment for level 1d data index %d at target time %s \n",i,timeBegin2);
			  arrLev1d[i] = NULL;
			  ActualnSegs1d-=1;
			}
		    }
		  else
		    {
		      arrLev1d[i] = NULL;
		      ActualnSegs1d-=1;
		    }
		}
	  
	      Segments1d = 1;        //now segments are in memory
	      status=1;
	      i=0;
	      while(status != DRMS_SUCCESS)
		{
		  TargetHFLID     = drms_getkey_int(recLev1d->records[i],FramelistID,&status);
		  if(i < nRecs1d-1) ++i;
		  else break;
		}
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: HFLID keyword cannot be read on level 1d data at target time %s\n",timeBegin2);
		  Segments1d = 0; 
		  goto NextTargetTime;
		}
	      else framelistSize = framelistInfo(TargetHFLID,PHWPLPOS,WavelengthIndex);
	    }
	  else                       //segments are already in memory
	    {
	      ActualnSegs1d=nRecs1d;
	      for(i=0;i<nRecs1d;++i)
		{
		  if(arrLev1d[i] == NULL) ActualnSegs1d-=1;
		}	      
	    }//end if(Segments1d == 0)
	  
      
      
	  //checking if some records are missing, and then if some segments are missing
	  //******************************************************************************
	  
	  if(ActualnSegs1d != framelistSize)
	    {
	      printf("Error: some records for the level 1d filtergrams are missing to produce level 1p data at target time %s\n",timeBegin2);
	      Segments1d=0;
	      goto NextTargetTime; //MAYBE DO SOMETHING ELSE THAN JUST DROPPING ALL THE LEV1P CALCULATIONS ?
	    }
      
      	  
	  //creating the record for the level 1p data (NB: 1 RECORD HOLDS ALL THE WAVELENGTHS)
	  //**********************************************************************************
	 
	  if(PolarizationType ==3) //producing LCP+RCP from 2 polarizations
	    {
	      nRecs1p=ActualnSegs1d/2; //nrecs1p is actually a number of groups of data segments (i.e. number of different wavelengths), not a number of records (NAME POORLY CHOSEN!)
	      npol=2;
	      npolout=2;
	      nSegs1p=npolout*nRecs1p;
	      Lev1pOffset=24;
	      if (Lev1pWanted) recLev1p = drms_create_records(drms_env,1,HMISeriesLev1pb,DRMS_PERMANENT,&status); //we create just one record that will contain several segments  
	      else recLev1p = drms_create_records(drms_env,1,HMISeriesLev1pb,DRMS_TRANSIENT,&status);
	    }
	  if(PolarizationType ==2)
	    {
	      nRecs1p=ActualnSegs1d/4; //producing LCP+RCP from 4 polarizations
	      npol=4;
	      npolout=2;
	      nSegs1p=npolout*nRecs1p;
	      Lev1pOffset=24;
	      if (Lev1pWanted) recLev1p = drms_create_records(drms_env,1,HMISeriesLev1pb,DRMS_PERMANENT,&status); //we create just one record that will contain several segments  
	      else recLev1p = drms_create_records(drms_env,1,HMISeriesLev1pb,DRMS_TRANSIENT,&status);
	    }
	  if(PolarizationType ==1) //producing I+Q+U+V
	    {
	      nRecs1p=ActualnSegs1d/4;
	      npol=4; 
	      npolout=4;
	      nSegs1p=npolout*nRecs1p;
	      Lev1pOffset=0;
	      if (Lev1pWanted) recLev1p = drms_create_records(drms_env,1,HMISeriesLev1pa,DRMS_PERMANENT,&status); //we create just one record that will contain several segments 
	      else recLev1p = drms_create_records(drms_env,1,HMISeriesLev1pa,DRMS_TRANSIENT,&status);
	    }
	  
	  if (status != DRMS_SUCCESS || recLev1p == NULL)
	    {
	      printf("Could not create a record for the level 1p series at target time %s\n",timeBegin2);
	      Segments1d=0;
	      goto NextTargetTime;
	    }	  
	  
	  //Creating arrays of input and output images for Jesper's code
	  //**************************************************************
	  
	  arrLev1p = (DRMS_Array_t **)malloc(nSegs1p*sizeof(DRMS_Array_t *));
	  if(arrLev1p == NULL)
	    {
	      printf("Error: memory could not be allocated to arrLev1p\n");
	      exit(EXIT_FAILURE);
	    }
	  
	  for(i=0;i<nSegs1p;++i)
	    {
	      arrLev1p[i] = drms_array_create(type1p,2,axisout,NULL,&status);
	      if(status != DRMS_SUCCESS || arrLev1p[i] == NULL)
		{
		  printf("Error: cannot create an array for a level 1p data at target time %s\n",timeBegin2);
		  for(ii=0;ii<nRecs1d;++ii)
		    {
		      if(arrLev1d[ii] != NULL)
			{
			  drms_free_array(arrLev1d[ii]);
			  arrLev1d[ii]=NULL;
			  Segments1d=0;
			}
		    }
		  Segments1p=0;
		  goto NextTargetTime;
		}
	    }
	  
	  images = (float **)malloc(npol*sizeof(float *));
	  if(images == NULL)
	    {
	      printf("Error: memory could not be allocated to images\n");
	      exit(EXIT_FAILURE);
	    }

	  imagesout = (float **)malloc(npolout*sizeof(float *));
	  if(imagesout == NULL)
	    {
	      printf("Error: memory could not be allocated to imagesout\n");
	      exit(EXIT_FAILURE);
	    }
	  
	  ps1 = (int *)malloc(npol*sizeof(int *));
	  if(ps1 == NULL)
	    {
	      printf("Error: memory could not be allocated to ps1\n");
	      exit(EXIT_FAILURE);
	    }
	  ps2 = (int *)malloc(npol*sizeof(int *));
	  if(ps2 == NULL)
	    {
	      printf("Error: memory could not be allocated to ps2\n");
	      exit(EXIT_FAILURE);
	    }
	  ps3 = (int *)malloc(npol*sizeof(int *));
	  if(ps3 == NULL)
	    {
	      printf("Error: memory could not be allocated to ps3\n");
	      exit(EXIT_FAILURE);
	    }
	  
	  fid = (int *)malloc(nRecs1d*sizeof(int *));
	  if(fid == NULL)
	    {
	      printf("Error: memory could not be allocated to fid\n");
	      exit(EXIT_FAILURE);
	    }
	  
	  Wavelengths = (int *)malloc(nRecs1p*sizeof(int *));
	  if(Wavelengths == NULL)
	    {
	      printf("Error: memory could not be allocated to Wavelengths\n");
	      exit(EXIT_FAILURE);
	    }


          //look for the available arrays and what are the wavelengths of the level 1d filtergrams
	  ii=0;
	  for(i=0;i<nRecs1d;++i)
	    {
	      if (arrLev1d[i] != NULL) 
		{
		  fid[i]=drms_getkey_int(recLev1d->records[i],FiltergramID,&status);     //to know which wavelength the filtergrams represent
		  if( status != 0)
		    {
		      printf("Error: unable to read the keyword fid in a level 1d record\n");
		      goto NextTargetTime;
		    }
		  if(i == 0) Wavelengths[0]=WhichWavelength(fid[0]);
		  temp=0;
		  for(k=0;k<=ii;++k) if(WhichWavelength(fid[i]) == Wavelengths[k]) temp=1;
		  if(temp == 0)
		    {
		      ii+=1;
		      Wavelengths[ii]=WhichWavelength(fid[i]);
		    }
		}
	      else
		{
		  goto NextTargetTime; 
		}
	    }

	  if(ii+1 != nRecs1p) //the number of wavelengths is not what it should be
	    {
	      printf("Error: the number of wavelengths in level 1d data: %d; is not what it should be: %d\n",ii+1,nRecs1p);
	      goto NextTargetTime; 
	    }

	  
	  //propagate keywords from level 1d data to level 1p	  
	  drms_copykeys(recLev1p->records[0],recLev1d->records[0],1, kDRMS_KeyClass_Explicit); //we copy all the keywords from the level 1d data to this level 1p record
	  
	  for(k=0;k<nRecs1p;++k) //loop over the groups of level 1p data segments (LCP+RCP or I+Q+U+V). Should be the number of different wavelengths
	    {

	      i=0;
	      for(ii=0;ii<nRecs1d;++ii) if (WhichWavelength(fid[ii]) == Wavelengths[k])       //find out which wavelength
		{
		  printf("wavelength=%d, polarization %d\n",Wavelengths[k],i);
		  images[i]=arrLev1d[ii]->data;
		  ps1[i]=drms_getkey_int(recLev1d->records[ii],PolarizationPos1,&statusA[0]); //WARNING: MODIFIY TO ACCOUNT FOR POTENTIAL ERRORS
		  ps2[i]=drms_getkey_int(recLev1d->records[ii],PolarizationPos2,&statusA[1]);
		  ps3[i]=drms_getkey_int(recLev1d->records[ii],PolarizationPos3,&statusA[2]);
		  if( (statusA[0]+statusA[1]+statusA[2]) != 0)
		    {
		      printf("Error: unable to read one or several keyword(s) in level 1d data\n");
		      goto NextTargetTime;
		    }
		  i+=1;
		}		      	      
	      
	      for(i=0;i<npolout;++i) imagesout[i]=arrLev1p[k*npolout+i]->data;


	      //Calling Jesper's code
	      //**************************************************************
	      tsel=28.;  //polarization selector temperature//WARNING: NEED TO MODIFY
	      tfront=28.;//front window temperature //WARNING: NEED TO MODIFY
	      
	      printf("Producing level 1p data\n");
	      polcal(&pars,npol,PolarizationType,images,imagesout,ps1,ps2,ps3,tsel,tfront,axisout[0],axisout[1],axisout[1]);	        
	      
	      //Putting output images in the proper records
	      //**************************************************************
	      	
	      if(Lev1pWanted) //if required, write the segment on file
		{
		  for(i=0;i<npolout;++i)
		    {		      
		      segout = drms_segment_lookup(recLev1p->records[0],Lev1pSegName[i+Lev1pOffset+Wavelengths[k]*npolout]);
		      drms_segment_write(segout,arrLev1p[k*npolout+i], 0);        //write the file containing the data (WE ASSUME THAT imagesout ARE IN THE ORDER I,Q,U,V OR LCP,RCP)		  
		    }	      
		  
		}
	      
	    }//end of for(k=0;k<nRecs1p;++k)
	 	  
	  Segments1p=1; //data segments for level 1p data are in memory  
	  
	}//end of producing level 1p data
  


      /****************************************************************************************************************************/
      /*                                                                                                                          */
      /*                                                                                                                          */
      /* IF LEVEL 1.5 OUTPUT DESIRED                                                                                              */
      /*                                                                                                                          */
      /*                                                                                                                          */
      /****************************************************************************************************************************/
   

      if (Lev15Wanted)
	{
	  // Given consistency tests and the preceeding code we know that lev1p exists and that PolarizationType = 1 (if the level 1p data were computed rather than read. If read, we don't care about the value of PolarizationType)

	  //reading the level 1p segments
	  //******************************************************************************

	  int ActualnSegs1p=nSegs1p;

	  if(Segments1p == 0)
	    {

	      if(PolarizationType ==2 || PolarizationType ==3) nSegs1p=nRecs1p*2; //number of level 1p data segments that are not a NULL pointer
	      if(PolarizationType ==1) nSegs1p=nRecs1p*4;

	      arrLev1p = (DRMS_Array_t **)malloc(nSegs1p*sizeof(DRMS_Array_t *));
	      if(arrLev1p == NULL)
		{
		  printf("Error: memory could not be allocated to arrLev1p\n");
		  exit(EXIT_FAILURE);
		}
	      for (i=0;i<nSegs1p;++i)
		{
		  arrLev1p[i] = drms_array_create(type1p,2,axisout,NULL,&status);
		  if(status != DRMS_SUCCESS || arrLev1p[i] == NULL)
		    {
		      printf("Error: cannot create an array for a level 1p data at target time %s\n",timeBegin2);
		      goto NextTargetTime;
		    }
		}

	      for(i=0;i<nSegs1p;++i)
		{
		  segin   = drms_segment_lookupnum(recLev1p->records[0], i);
		  arrLev1p[i] = drms_segment_read(segin, segin->info->type, &status); //pointer toward the segment
		  if(status != DRMS_SUCCESS || arrLev1p[i] == NULL)
		    {
		      printf("Error: could not read the segment for level 1p data index %d at target time %s \n",i,timeBegin2);
		      arrLev1p[i] = NULL;
		      ActualnSegs1p-=1;
		    }
		}
	      //ActualnSegs1p is now the actual number of segments

	      TargetHFLID = drms_getkey_int(recLev1p->records[0],FramelistID,&status);
	      if(status != DRMS_SUCCESS)
		{
		  printf("Error: HFLID keyword cannot be read on level 1p data at target time %s\n",timeBegin2);
		  goto NextTargetTime;
		}
 
	      Segments1p=1;  //now the data segments for level 1p data are in memory
	      framelistSize = framelistInfo(TargetHFLID,PHWPLPOS,WavelengthIndex);
	      
	    }//end of if(Segments1p == 0)

	  //checking if some level 1p segments are missing
	  //******************************************************************************
	  
	  if(ActualnSegs1p != nSegs1p) //some segments are missing
	    {
	      printf("Error: some level 1p data are missing to produce level 1.5 data: %d\n",ActualnSegs1p);
	      goto NextTargetTime; //MAYBE DO SOMETHING ELSE THAN JUST DROPPING ALL THE LEV1P CALCULATIONS ?
	    }
	  
	  //reading the appropriate look-up table for the MDI-like algorithm
	  //******************************************************************************

	  i=0;
	  while(WavelengthIndex[i] != 2) i++; //i contains the index of I2
	  int HCMNBT,HCMWBT,HCMPOLT,HCME1T;
	  if(framelistSize == 12) // 6 wavelengths
	    {
	      HCMNBT = PHWPLPOS[7*i+3]-12;
	      HCMWBT = PHWPLPOS[7*i+1]-6;
	      HCMPOLT= PHWPLPOS[7*i+2];
	      HCME1T = PHWPLPOS[7*i+0]+3;
	    }
	  if(framelistSize == 10) // 5 wavelengths
	    {
	      HCMNBT = PHWPLPOS[7*i+3];
	      HCMWBT = PHWPLPOS[7*i+1];
	      HCMPOLT= PHWPLPOS[7*i+2];
	      HCME1T = PHWPLPOS[7*i+0];
	    }
	  printf("%d %d %d %d\n",HCMNBT,HCMWBT,HCME1T,HCMPOLT);

	  char *filter=NULL;
	  char *where       = "FSN_REC BETWEEN 706315 AND 800000"; //WARNING, CHANGE THAT!!!
	  char *HCMNBs      = "HWL4POS";          //keyword for the position of the HCM of NB Michelson
	  char *HCMPOLs     = "HWL3POS";          //keyword for the position of the HCM of the tuning polarizer 
	  char *HCMWBs      = "HWL2POS";          //keyword for the position of the HCM of WB Michelson
	  char *HCME1s      = "HWL1POS";          //keyword for the position of the HCM of Lyot E1
	  char *Ns          = "NWL";              //keyword for the number of wavelength

	  int mixed=1;
	  int allvers = 0;
	  DB_Text_Result_t *tres,*tres1,*tres2,*tres3,*tres4,*tres5;
	  char *query,*query1,*query2,*query3,*query4,*query5;
	  double *count=NULL;
	  int row,col,row1,col1,row2,col2,row3,col3,row4,col4,row5,col5;
	  int NBC,WBC,E1C,POLC,NC;
	  query =drms_query_string(drms_env,HMISeriesLookup, where, filter, mixed,DRMS_QUERY_FL, NULL,TREC,    allvers); //get the data records as a function of TREC
	  query1=drms_query_string(drms_env,HMISeriesLookup, where, filter, mixed,DRMS_QUERY_FL, NULL,HCMNBs,  allvers); //get the data records as a function of HCMNB
	  query2=drms_query_string(drms_env,HMISeriesLookup, where, filter, mixed,DRMS_QUERY_FL, NULL,HCMWBs,  allvers); //get the data records as a function of HCMWB
	  query3=drms_query_string(drms_env,HMISeriesLookup, where, filter, mixed,DRMS_QUERY_FL, NULL,HCMPOLs, allvers); //get the data records as a function of HCMPOL
	  query4=drms_query_string(drms_env,HMISeriesLookup, where, filter, mixed,DRMS_QUERY_FL, NULL,HCME1s,  allvers); //get the data records as a function of HCME1
	  query5=drms_query_string(drms_env,HMISeriesLookup, where, filter, mixed,DRMS_QUERY_FL, NULL,Ns,  allvers);     //get the data records as a function of N (number of wavelengths)

	  tres  = drms_query_txt(drms_env->session,query);
	  tres1 = drms_query_txt(drms_env->session,query1);
	  tres2 = drms_query_txt(drms_env->session,query2);
	  tres3 = drms_query_txt(drms_env->session,query3);
	  tres4 = drms_query_txt(drms_env->session,query4);
	  tres5 = drms_query_txt(drms_env->session,query5);

	  row =tres->num_rows;
	  col =tres->num_cols;
	  row1=tres1->num_rows;
	  col1=tres1->num_cols;
	  row2=tres2->num_rows;
	  col2=tres2->num_cols;
	  row3=tres3->num_rows;
	  col3=tres3->num_cols;
	  row4=tres4->num_rows;
	  col4=tres4->num_cols;
	  row5=tres5->num_rows;
	  col5=tres5->num_cols;

	  if(col != col1 || col != col2 || col != col3 || col != col4 || col != col5 || row != row1 || row != row2 || row != row3 || row != row4 || row != row5)
	    {
	      printf("Error: the records of the look-up table series %s do not have the same number of keywords needed by this program\n",HMISeriesLookup);
	      exit(EXIT_FAILURE);
	    }

	  printf("row= %d, col=%d\n",row,col);
	  if(col != 1)
	    {
	      printf("Error: problem with the query to retrieve the look-up tables\n");
	      exit(EXIT_FAILURE);
	    }

	  count = (double *)malloc(row*sizeof(double));
	  if(count == NULL)
	    {
	      printf("Error: memory could not be allocated to count\n");
	      exit(EXIT_FAILURE);
	    }

	  temptime = 1000000000.0;
	  temp = 0;
	  for(i=0;i<row;++i)
	    {
	      count[i] = fabs(atof(tres->field[i][0])-TargetTime); //absolute value of the difference between the T_REC of the look-up tables and the TargetTime
	      NBC = atoi(tres1->field[i][0])-HCMNBT; 
	      WBC = atoi(tres2->field[i][0])-HCMWBT;
	      POLC= atoi(tres3->field[i][0])-HCMPOLT;
	      E1C = atoi(tres4->field[i][0])-HCME1T;
	      NC  = atoi(tres5->field[i][0])-framelistSize/2;
	      if(count[i] < temptime && NBC+WBC+POLC+E1C+NC == 0)
		{
		  temptime=count[i];
		  temp=i; //temp will contain the index at which the T_REC value of the look-up table is closest to the TargetTime
		}
	    }
	  sprint_ut(query,atof(tres->field[temp][0])); 
	  strcpy(HMILookup,HMISeriesLookup); 
	  strcat(HMILookup,"[]["); //assumes that the look-up table series has two prime keywords, and T_REC is the second one!!!
	  strcat(HMILookup,query);
	  strcat(HMILookup,"]");

	  printf("Look-up table query= %s\n",HMILookup);

	  free(count);

	  lookup  = drms_open_records(drms_env,HMILookup,&status); 
	  if (status == DRMS_SUCCESS && lookup != NULL)
	    {
	      if (lookup->n > 1) 
		{
		  printf("Error: more than 1 lookup table record was downloaded.\n");
		  exit(EXIT_FAILURE);
		}
	      if (lookup->n <= 0) 
		{
		  printf("Error:no record for the look-up tables were downloaded.\n");
		  exit(EXIT_FAILURE);
		}
	    }
	  else
	    {
	      printf("Error: can't open the look-up table series.\n");
	      exit(EXIT_FAILURE);
	    }
	  
	  segin     = drms_segment_lookupnum(lookup->records[0], 0);
	  arrintable= drms_segment_read(segin, segin->info->type, &status);
	  if (status != DRMS_SUCCESS)
	    {
	      printf("Error: unable to read the data segment of the look-up table record\n"); //if there is a problem
	      exit(EXIT_FAILURE);               
	    } 
	  else printf("look-up table record read\n");
	  
	  //producing the level 1.5 data
	  //******************************************************************************
	  
	  nRecs15   = 5; //Dopplergram, l.o.s. magnetic field, linewidth, linedepth, continuum
	  
	  recLev15a = drms_create_records(drms_env,1,HMISeriesLev15a,DRMS_PERMANENT,&statusA[0]); //RECORD FOR DOPPLERGRAM
	  recLev15b = drms_create_records(drms_env,1,HMISeriesLev15b,DRMS_PERMANENT,&statusA[1]); //RECORD FOR MAGNETOGRAM
	  recLev15c = drms_create_records(drms_env,1,HMISeriesLev15c,DRMS_PERMANENT,&statusA[2]); //RECORD FOR LINEDEPTH
	  recLev15d = drms_create_records(drms_env,1,HMISeriesLev15d,DRMS_PERMANENT,&statusA[3]); //RECORD FOR LINEWIDTH
	  recLev15e = drms_create_records(drms_env,1,HMISeriesLev15e,DRMS_PERMANENT,&statusA[4]); //RECORD FOR CONTINUUM

	  if ( (statusA[0]+statusA[1]+statusA[2]+statusA[3]+statusA[4]) != DRMS_SUCCESS || recLev15a == NULL || recLev15b == NULL || recLev15c == NULL|| recLev15d == NULL || recLev15e == NULL)
	    {
	      printf("Could not create a record for one or several level 1.5 data series, at target time %s\n",timeBegin2); 
	      goto NextTargetTime;
	    }


	  arrLev15 = (DRMS_Array_t **)malloc(nRecs15*sizeof(DRMS_Array_t *));
	  if(arrLev15 == NULL)
	    {
	      printf("Error: memory could not be allocated to arrLev15\n");
	      exit(EXIT_FAILURE);
		}
	  for (i=0;i<nRecs15;++i)
	    {
	      arrLev15[i] = drms_array_create(type15,2,axisout,NULL,&status);
	      if(status != DRMS_SUCCESS || arrLev15[i] == NULL)
		{
		  printf("Error: cannot create an array for a level 1.5 data at target time %s\n",timeBegin2);
		  goto NextTargetTime;
		}
		  
	    }
	      
	  //read the R_SUN, X0, and Y0 keywords
	  SUNRADIUS = drms_getkey_float(recLev1p->records[0],Rsun,&status); //all the Lev1p records for a same observable should have the same values of R_SUN, X0, and Y0
	  if(status != DRMS_SUCCESS)
	    {
	      printf("Error: %s keyword cannot be read on level 1p data at target time %s\n",Rsun,timeBegin2);
	      goto NextTargetTime;
	    }

	  SUNCENTERX = drms_getkey_float(recLev1p->records[0],x0,&status);
	  if(status != DRMS_SUCCESS)
	    {
	      printf("Error: %s keyword cannot be read on level 1p data at target time %s\n",x0,timeBegin2);
	      goto NextTargetTime;
	    }

	  SUNCENTERY = drms_getkey_float(recLev1p->records[0],y0,&status);
	  if(status != DRMS_SUCCESS)
	    {
	      printf("Error: %s keyword cannot be read on level 1p data at target time %s\n",y0,timeBegin2);
	      goto NextTargetTime;
	    }


	  //read the parameters defined in HMIparam.h
	  DopplerParameters.FSRNB=FSR[0];
	  DopplerParameters.FSRWB=FSR[1];
	  DopplerParameters.FSRE1=FSR[2];
	  DopplerParameters.FSRE2=FSR[3];
	  DopplerParameters.FSRE3=FSR[4];
	  DopplerParameters.FSRE4=FSR[5];
	  DopplerParameters.FSRE5=FSR[6];
	  DopplerParameters.dlamdv=dlamdv;
	  DopplerParameters.maxVtest=maxVtest;
	  DopplerParameters.maxNx=maxNx;
	  DopplerParameters.ntest=ntest;
	  DopplerParameters.dvtest=dvtest;
	  DopplerParameters.MISSINGDATA=MISSINGDATA;
	  DopplerParameters.MISSINGRESULT=MISSINGRESULT;

	  Dopplergram(arrLev1p,arrLev15,nSegs1p,arrintable,SUNRADIUS,SUNCENTERX,SUNCENTERY,DopplerParameters);

	  //setting keywords
	  fsn = drms_getkey_int(recLev1p->records[0],"FSN",&status);
	  if(status != DRMS_SUCCESS)
	    {
	      printf("Error: FSN keyword cannot be read on level 1p data at target time %s\n",timeBegin2);
	      goto NextTargetTime;
	    }
	  statusA[0]  = drms_setkey_int(recLev15a->records[0],"FSN_START",fsn); //Dopplergram
	  statusA[1]  = drms_setkey_int(recLev15b->records[0],"FSN_START",fsn); //Magnetogram
	  statusA[2]  = drms_setkey_int(recLev15c->records[0],"FSN_START",fsn); //Linedepth
	  statusA[3]  = drms_setkey_int(recLev15d->records[0],"FSN_START",fsn); //Linewidth
	  statusA[4]  = drms_setkey_int(recLev15e->records[0],"FSN_START",fsn); //Continuum
	  statusA[5]  = drms_setkey_time(recLev15a->records[0],TREC,trec);
	  statusA[6]  = drms_setkey_time(recLev15b->records[0],TREC,trec);
	  statusA[7]  = drms_setkey_time(recLev15c->records[0],TREC,trec);
	  statusA[8]  = drms_setkey_time(recLev15d->records[0],TREC,trec);
	  statusA[9]  = drms_setkey_time(recLev15e->records[0],TREC,trec);
	  statusA[10] = drms_setkey_time(recLev15a->records[0],TOBS,tobs);
	  statusA[11] = drms_setkey_time(recLev15b->records[0],TOBS,tobs);
	  statusA[12] = drms_setkey_time(recLev15c->records[0],TOBS,tobs);
	  statusA[13] = drms_setkey_time(recLev15d->records[0],TOBS,tobs);
	  statusA[14] = drms_setkey_time(recLev15e->records[0],TOBS,tobs);
	  TotalStatus=0;
	  for(i=0;i<15;++i) TotalStatus+=statusA[i];
	  if(TotalStatus != 0)
	    {
	      printf("Error: unable to set some keywords for the level 1.5 data\n");
	      goto NextTargetTime;
	    }


	  //writing lev1.5 data
	  segout = drms_segment_lookupnum(recLev15a->records[0], 0);
	  drms_segment_write(segout,arrLev15[0], 0);
	  //printf("%ld %s\n",recLev15a->records[0]->sunum,recLev15a->records[0]->su->sudir);
	  segout = drms_segment_lookupnum(recLev15b->records[0], 0);
	  drms_segment_write(segout,arrLev15[1], 0);
	  segout = drms_segment_lookupnum(recLev15c->records[0], 0);
	  drms_segment_write(segout,arrLev15[2], 0);
	  segout = drms_segment_lookupnum(recLev15d->records[0], 0);
	  drms_segment_write(segout,arrLev15[3], 0);
	  segout = drms_segment_lookupnum(recLev15e->records[0], 0);
	  drms_segment_write(segout,arrLev15[4], 0);  
	  
	}//end of producing the level 1.5 data
      
      
      
    NextTargetTime:
      
      
      /****************************************************************************************************************************/
      /*                                                                                                                          */
      /*                                                                                                                          */
      /* FREEING THE RECORDS                                                                                                      */
      /*                                                                                                                          */
      /*                                                                                                                          */
      /****************************************************************************************************************************/

      printf("FREEING RECORD\n");

      if(recLev1d != NULL) //if input data level is either level 1 or level 1d
	{
	  if(recLev1d->n > 0)
	    {
	      //if lev1d data required as output, then we save them in the DRMS
	      if (Lev1dWanted) status=drms_close_records(recLev1d,DRMS_INSERT_RECORD);     //insert the record in DRMS if the record has been marked as PERMANENT
	      else status=drms_close_records(recLev1d,DRMS_FREE_RECORD);
	      recLev1d=NULL;
	      for(i=0;i<nRecs1d;++i) if(arrLev1d[i] != NULL)
		{
		  drms_free_array(arrLev1d[i]); //also frees images
		}
	      if(arrLev1d != NULL) free(arrLev1d);
	      arrLev1d=NULL;
	      Segments1d=0;
	    }
	}


      if(recLev1p != NULL)
	{
	  if(recLev1p->n > 0)
	    {
	      if (Lev1pWanted) status=drms_close_records(recLev1p,DRMS_INSERT_RECORD);	//if lev1p data are a requested output, then they need to be recorded
	      else  status=drms_close_records(recLev1p,DRMS_FREE_RECORD);
	      recLev1p=NULL;
	      for(i=0;i<nSegs1p;++i) if(arrLev1p[i] != NULL) //also frees imagesout
		{
		  drms_free_array(arrLev1p[i]);
		}
	      if(arrLev1p != NULL) free(arrLev1p);
	      arrLev1p=NULL;
	      if(ps1 != NULL) free(ps1);
	      if(ps2 != NULL) free(ps2);
	      if(ps3 != NULL) free(ps3);
	      if(fid != NULL) free(fid);
	      if(Wavelengths != NULL) free(Wavelengths);
	      Wavelengths=NULL;
	      images=NULL;
	      imagesout=NULL;
	      ps1=NULL;
	      ps2=NULL;
	      ps3=NULL;
	      fid=NULL;
	      Segments1p=0;
	    }
	}
      

      if(recLev15a != NULL)
	{
	  if(recLev15a->n > 0)
	    {
	      status=drms_close_records(recLev15a,DRMS_INSERT_RECORD);
	      status=drms_close_records(recLev15b,DRMS_INSERT_RECORD);
	      status=drms_close_records(recLev15c,DRMS_INSERT_RECORD);
	      status=drms_close_records(recLev15d,DRMS_INSERT_RECORD);
	      status=drms_close_records(recLev15e,DRMS_INSERT_RECORD);
	      recLev15a=NULL;
	      recLev15b=NULL;
	      recLev15c=NULL;
	      recLev15d=NULL;
	      recLev15e=NULL;
	      for (i=0;i<nRecs15;++i) if(arrLev15[i] != NULL)
		{
		  drms_free_array(arrLev15[i]);
		}
	      if(arrLev15 != NULL) free(arrLev15);
	      arrLev15=NULL;
	      if(arrintable != NULL) drms_free_array(arrintable);
	      arrintable=NULL;
	      status=drms_close_records(lookup,DRMS_FREE_RECORD); 
	    }
	}



      TargetTime+=DopplergramCadence; 
 
   }//end while(TargetTime <= TimeEnd)


  if (TestLevIn[0]==1) //input data are level 1 filtergrams
    {
      if(recLev1->n > 0)
	{
	  status=drms_close_records(recLev1,DRMS_FREE_RECORD);  
	  recLev1=NULL;
	  free(internTOBS);
	  free(HWL1POS); 
	  free(HWL2POS);
	  free(HWL3POS);
	  free(HWL4POS);
	  free(HPL1POS); 
	  free(HPL2POS);
	  free(HPL3POS);
	  free(FID);
	  free(HFLID);
	  free(RSUN);
	  free(CROTA1);
	  free(CRLTOBS);
	  free(DSUNOBS);
	  free(X0);
	  free(Y0);
	  free(SegmentRead);
	  free(Segments);
	  free(Badkeyword);
	  free(Ierror);  
	  free(IndexFiltergram);
	}
    }

  if(TestLevIn[0]==1)
    {
      free_interpol(&const_param);
    }

  if(Lev1pWanted || (Lev15Wanted && TestLevIn[2]==0))
    {
      status = free_polcal(&pars);
    }

  status=0;
  return status;

}


