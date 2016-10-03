// USING A VOIGT PROFILE FOR THE Fe I LINE

/*-------------------------------------------------------------------------------------------------------*/

/*  Modified from JSOC/proj/lev1.5_hmi/apps/phasemaps_voight.c

/*-------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------
 * Changes needed for OBS_MODE  -- Phil Jan 2016
 * OBS_MODE detunes use the same sequence, 27 exposures with a dark at each end.
 * In OBS_MODE the WCS keywords are populated.
 * Due to CRPIX2 jitter of 0.75 pixel and about half that for CRPIX1 the image must
 * be smoothed perhaps as well as binning to eliminate high gradients in LOS velocity
 * Changes:
 *   0.  Use HCFTID to determine CAL_MODE or OBS_MODE
 *   1.  Must compute LOS velocity for each location
 *   2.  Must set solarradiusmax to R_SUN in right units...
 *-------------------------------------------------------------------------------------------------------
 */
/* XXXXXX Changes made in Jan 2016 by Phil
 * 1. changed comment for HCAMID
 * 2. Changes all instances in DoIt of 'exit(' to 'return('
 * 3. Add code from sunrotation.c for orbit removal and average sun doppler removal
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort_double.h>//GSL = Gnu Scientifc Library
#include <gsl/gsl_blas.h>       //GNU implementation of the BLAS routines
#include <gsl/gsl_linalg.h>     //SVD algorithm from the GNU Scientific Library
#include <gsl/gsl_statistics.h> //mean and sigma

#include "mdi.h"           //contains definitions for some HMI filter functions
#include "mdi_params.h"      //contains definitions for some HMI filter parameters
#include "mdi.c"           //contains some HMI filter functions

double Ni_Line   =  LAMBDA_Ni; //Ni I line central wavelength in vacuum at rest, 6767.78 A 

#include <omp.h>                //Open MP header
#include <fresize.h>            //from Jesper: to rebin the 1024x1024 images
#include "interpol_code.h"      //from Richard, for de-rotation and gap-filling

char *module_name    = "MDI_phasemaps";   //name of the module

#define kRecSetIn      "input_series" //names of the arguments of the module
#define kDSOut         "phasemap_series"
#define kreduced       "reduced"      //maps in 64x64 or 256x256?             
#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) < (y)) ? (y) : (x))

#define Deg2Rad  (M_PI/180.0)

//arguments of the module
ModuleArgs_t module_args[] =        
{
     {ARG_STRING, "in", "" ,  "Input data series."},
     {ARG_STRING, "out",    "" ,  "Phase Maps series."},
     {ARG_INT   , kreduced,  "0",  "64x64 (1), 32x32 (2), 128x128 (3), or standard (256x256) resolution (0)?"},
     {ARG_DOUBLE , "FSRNB", "-9999.0" , "FSR NB"},
     {ARG_DOUBLE , "FSRWB", "-9999.0" , "FSR WB"},
     {ARG_FLOAT  , "center", "" , "center of blocker filter"},
     {ARG_FLOAT  , "shift", "" , "wavelength shift of the non-tunable profile"},
     {ARG_DOUBLE , "thresh", "", "threshold"},
     {ARG_END}
};

// function to get LOS velocity of observed point
double Solar_velocity(double x,double y,double rsun2,double OBS_B0,double CROTA2, double CDELT, double obs_vr, double obs_vw, double obs_vn);
int img2sphere (double x, double y, double ang_r, double latc, double lonc,
    double pa, double *p_rho, double *p_lat, double *p_lon, double *p_sinlat,
    double *p_coslat, double *p_sig, double *p_mu, double *p_chi);
// rebin function that deals with missing data.
int rebinArraySF(float *outData, float *inData, int outNx, int inNx, float rsun, float x0, float y0, float fillvalue);

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


/*-----------------------------------------------------------------------------------------*/
/*									                   */
/* FUNCTION TO PRODUCE A LOOK-UP TABLE CORRESPONDING TO THE PHASES	                   */
/*									                   */
/* the convention for the filter names is (with 4 filters for Doppler):                    */
/* F1 is at -113.4 mA   M1 @ -27= -0.6 FSR1, M2 @ -13=-0.2889 FSR2                         */
/* F2 is at  -37.8 mA   M1 @ -9 = -0.2 FSR1, M2 @ -4 =-0.0889 FSR2                         */
/* F3 is at   37.8 mA   M1 @  9 =  0.2 FSR1, M2 @  5 = 0.1111 FSR2                         */
/* F4 is at  113.4 mA   M1 @  27=  0.6 FSR1, M2 @ 14 = 0.3111 FSR2                         */
/* Fc is at -184.3 mA   M1 @ -45= -1.0 FSR1, M2 @-22 =-0.4889 FSR2                         */
/*                                                                                         */
/* the relative phases of NB and WB increase when the steps of the HCMs                    */
/* increase                                                                                */
/*-----------------------------------------------------------------------------------------*/

void cotunetable(double NBphase,double WBphase,int table[2][20])
{

  int HCMNB, HCMWB, i;

  HCMNB    = (int)round(-NBphase/6.0)+60;            //for convention Intensity = 0.5*[1+cos(lam+phase)]
  HCMWB    = (int)round(-WBphase/6.0)+60;


  table[0][0]=-30; table[1][0]=0;
  table[0][1]=-27; table[1][1]=-6;
  table[0][2]=-24; table[1][2]=-12;
  table[0][3]=-21; table[1][3]=-18;
  table[0][4]=-18; table[1][4]=-24;
  table[0][5]=-15; table[1][5]=-30;    //I5 for N=6 tuning positions
  table[0][6]=-12; table[1][6]=24;   //I4 for N=5
  table[0][7]=-9;  table[1][7]=18;   //I4 for N=6
  table[0][8]=-6;  table[1][8]=12;   //I3 for N=5
  table[0][9]=-3;  table[1][9]=6;    //I3 for N=6
  table[0][10]=0;  table[1][10]=0;   //I2 for N=5
  table[0][11]=3;  table[1][11]=-6; //I2 for N=6
  table[0][12]=6;  table[1][12]=-12;//I1 for N=5
  table[0][13]=9;  table[1][13]=-18; //I1 for N=6
  table[0][14]=12; table[1][14]=-24;//I0 for N=5
  table[0][15]=15; table[1][15]=-30; //I0 for N=6
  table[0][16]=18; table[1][16]=24; 
  table[0][17]=21; table[1][17]=18;
  table[0][18]=24; table[1][18]=12;
  table[0][19]=27; table[1][19]=6;


  for(i=0;i<20;++i)
    {
      table[0][i] = table[0][i] + HCME1;
      table[1][i] = table[1][i] + HCMWB;
    }

}
  

/*------------------------------------------------------------------------------------------------------*/
/*                                                                                                      */
/*  MAIN PROGRAM                                                                                        */
/*                                                                                                      */
/*------------------------------------------------------------------------------------------------------*/


#define UNK_MODE 0
#define CAL_MODE 1
#define OBS_MODE 2

#define DIE(msg) {fprintf(stderr,msg); return(EXIT_FAILURE); }
#define DIE1(msg,val) {fprintf(stderr,msg,val); return(EXIT_FAILURE); }

int DoIt(void)
  {
  const char *inRecQuery     = cmdparams_get_str(&cmdparams, "in", NULL); //cmdparams is defined in jsoc_main.h
  const char *dsout          = cmdparams_get_str(&cmdparams, "out" ,NULL);    //series name of the output phase maps series
  int reduced          = cmdparams_get_int(&cmdparams  , kreduced ,NULL);  //if reduced=1 then we only compute 64x64 maps (nx2=64), =2 it's 32x32, =3 128x128, instead of the usual 256x256 (reduced =0)

  // Parameters defined in HMIparam.h with values for MDI
  float centerblocker2 =  cmdparams_get_float(&cmdparams,"center", NULL);
  float shiftw         = cmdparams_get_float(&cmdparams, "shift" , NULL);
  double FSR[2] = {FSR1,FSR2};
  int nx2 = 128;		// phasemap dimensions
  int ny2 = 128;		// phasemap dimensions
  int nelement = 2;		// number tunable elements
  int ntest = 801;		// number of test velocities
  double dvtest = 10.0		// m/s, so -4000:4000 m/s 
  double vtestmax = 4000.0
  int nfront = 201;		// number of wavelengths in front window transmission profile
  int nblocker = 401;		// number of wavelengths in blocker and Lyot profile, +- 1A so 5 mA resolution.
  double centerblocker = LCEN	// Lyot center wavelength WRT Ni line 0 velocity.
  double dlamdv = Ni_Line/C;	// 2.2574e-5 Angstroms per m/s
  int numberoffiltergrams_detune = 20 ; // 9 for detune, 1 FC, 10 for Dop sequence.
  double solar_radius = 696000000.0;
  double Astrounit = 149597870691.0;
  float phaseguess[3] = {272.0, 24.0};	// initial guess for cotune phase, assuming motor positions 79, 48
  double depth0 = 0.64;		// Ni line depth (I0 = 0.4)
  double widtho = 0.116;	// Ni Line FWHM Angstroms, (Norton et al)
  double thresh = 60000;	// guess for continuum
  float MINRSUN = 0.0;
  float MAXRSUN = 512;
  float XYMIN = 0;
  float XYMAX = 1024;
  double NOMINALSCALE = 2.0;	// scale for calmode, arcsec/pixel

  // Defaults for the following 8 params are in HMIparam.h
  double try;
  try = cmdparams_get_double(&cmdparams,"FSRNB" , NULL);
    if (try != -9999.0) FSR[0] = try;
  try = cmdparams_get_double(&cmdparams,"FSRWB" , NULL);
    if (try != -9999.0) FSR[1] = try;
  try = cmdparams_get_double(&cmdparams,"thresh", NULL);
    if (try != -9999.0) thresh = try;
  printf("FSR[0]=%f, FSR[1] = %f\n",FSR[0],FSR[1]);

  char COMMENT[256];
  
  char *COMMENTS= "COMMENT";

  printf("FSR PARAMETERS = %f %f %f %f %f\n",FSR[0],centerblocker2,shiftw,thresh);
  printf("COMMENT: %s\n",COMMENT);

  if(reduced == 1)
    {
      nx2=64;
      ny2=64;
    }
  if(reduced == 2)
    {
      nx2=32;
      ny2=32;
    }
  if(reduced == 3)
    {
      nx2=128;
      ny2=128;
    }

  int  status       = DRMS_SUCCESS; 
  int  error        = 0;
  int  nx           = 1024;                                           //number of columns of input lev 1 filtergrams
  int  ny           = 1024;                                           //number of rows of input lev 1 filtergrams
  int  Nelem        = nx*ny;
  int  ratio;                                                         //ratio nx/nx2 (MUST BE AN EVEN NUMBER)
  double iterstop   = 1.e-7;                                //criterion to stop the iterations: we stop when the maximum relative change in the parameters is less than iterstop

#define nseq  9                                                      //number of wavelength positions in the detune sequences (EXLUDING THE DARK FRAMES)
#define maxsteps 155                                                //maximum steps allowed for the iterative Least-Squares algorithm
#define nparam 6                                                      //number of parameters we fit for

  int nimages       = nseq+1;                                         //number of images in the detune sequence
  int indeximage[nimages];
  double factor     = 1000.;                                          //to make sure all the fitted parameters are, roughly, of the same order of magnitude
  double dpi        = 2.0*M_PI;
  int i,j,k,iii,jjj;                                                  //variables for the loops
  float distance2[nx*ny];
  float distance[ny2*nx2];                                            //format: [row][column] (in idl notation: distance[nx2,ny2])
                                                                      //nx2 is a number of column, ny2 is a number of rows
  double Inten[ny2][nx2][nseq];
  float solarradiusmax;                                               //maximum radius at which the phases will be computed

  //FOR DO_INTERPOLATE()
#define TABLEPATH "/home/jsoc/cvs/Development/JSOC/proj/lev1.5_hmi/libs/lev15/"

  struct init_files initfiles;
  char MDISeriesTemp[256];

  //WAVELENGTH GRID PARAMETERS

#define  nlam 20002                                                    //MUST BE AN EVEN NUMBER
double dlam       = 0.0005;

  double lam[nlam];
  for(i=0;i<nlam;i++) lam[i] = ((double)i-((double)nlam-1.0)/2.0)*dlam;//symmetrical around 0 but NO 0.
 
  double blockerint[nlam];
  double templine[nlam];
  double line[nlam];
  double linea[nlam];
  double lineb[nlam];
  double linec[nlam];
  double lined[nlam];
  double gaussian[nlam];
  double l;
  double a=0.03; //damping for Voigt profile (value from Foy, 1972. Should be 0.0315 actually)

  double profilef[nlam];
  double profileg[nlam];
  double dlinedI0[nlam];
  double dlinedIc[nlam];
  double dlinedw[nlam];
  double err[maxsteps];
  double history[maxsteps][nparam];
  double ydefault    = 0.;
  int    converg;
  int    compteur;
  double dIntendI0[nseq];
  double dIntendw[nseq];
  double dIntendIc[nseq];
  double dIntendPhi0[nseq];
  double dIntendPhi1[nseq];
  double Icg[nx2][nx2];
  double fwidthg[nx2][nx2];
  double fdepthg[nx2][nx2];
  double residual;
  double jose;
  thresh=thresh/factor;   //we will divide the intensities by factor, thresh is defined in HMIparam.h
  size_t minimum=0;
  FILE *fp=NULL;
  int row, column;
  int indexref;

  int nthreads=omp_get_num_procs();         //number of threads supported by the machine where the code is running
  if (nthreads > 8) nthreads = 8;
  omp_set_num_threads(nthreads);            //set the number of threads to the maximum value
  printf("number of threads for OpenMP = %d\n",nthreads);

  //TEST THE REBINNING RATIO
  if( (nx % nx2) == 0)  ratio =nx/nx2;                                 //typically ratio=16, because nx2=256 and nx=1024
  else
    {
      printf("Error: nx = %d must be a multiple of nx2 = %d\n",nx,nx2);
      DIE("");
    }

  //for(i=0;i<nlam;++i) gaussian[i]=0.015*exp(-(lam[i]+0.225)*(lam[i]+0.225)/0.2/0.2)-0.004*exp(-(lam[i]-0.150)*(lam[i]-0.150)/0.22/0.22); //GAUSSIAN LINE ADDED TO VOIGT PROFILE TO SIMULATE ASYMMETRY
  //for(i=0;i<nlam;++i) gaussian[i]=0.015*exp(-(lam[i]+0.225)*(lam[i]+0.225)/0.2/0.2);
  //for(i=0;i<nlam;++i) gaussian[i]=0.0;

  /***********************************************************************************************************/
  /*DRMS ACCESS                                                                                              */
  /***********************************************************************************************************/

  printf("QUERY = %s\n",inRecQuery);
  DRMS_RecordSet_t *data = drms_open_records(drms_env,inRecQuery,&status);//open the records from the input series

  if (status == DRMS_SUCCESS && data != NULL && data->n > 0)
    {
      int nRecs = data->n;                                               //number of records in the input series 
      printf("Number of filtergrams to be read= %d \n",nRecs);           //check if the number of records is appropriate
      if (nRecs != numberoffiltergrams_detune)  
	  DIE1("Problem: program requires %d filtergrams to produce the phase maps\n",numberoffiltergrams_detune);

      //the format of the detune9 sequence is assumed to be front camera, then side camera, then front camera, and so on...
      for(i=0;i<nimages;++i) indeximage[i]=i

      DRMS_Record_t *rec[nimages];                                       //we only work with 1 camera (9 detune positions)
      DRMS_Record_t *rec2[2]     ;                                       //records for the first and last filtergrams of the detune sequence
      for(i=0;i<nimages;++i) rec[i] = data->records[indeximage[i]];      //records that will be used to compute the phase maps
      rec2[0] = data->records[0];                                        //records that will be used to set the T_REC and FSN_REC keywords for the phasemaps
      rec2[1] = data->records[nimages-1];

      DRMS_Segment_t *segin  = NULL;
      int  Frecnum          = 0;
      int  Lrecnum         = 0;
      int  Mrecnum         = 0;  // average recnum
      DRMS_Type_t type       = DRMS_TYPE_FLOAT;
      DRMS_Type_t typeEr     = DRMS_TYPE_CHAR;

      const char *HCMNB      = "MTM1";                            //keyword for the position of the HCM of NB Michelson
      const char *HCMWB      = "MTM2";                            //keyword for the position of the HCM of WB Michelson
      const char *CRPIX1     = "CRPIX1";                             //x-axis location of solar disk center (IN OBS_MODE) in pixels, starting at 1
      const char *CRPIX2     = "CRPIX2";                             //y-axis location of solar disk center in pixels
      const char *SCALE      = "CDELT1";                             //image scale in the x direction
      const char *VR         = "OBS_VR";                             //SDO orbital velocity away from Sun m/s 
      const char *VW         = "OBS_VW";                             //SDO orbital velocity in direction of increasing HGLN_OBS 
      const char *VN         = "OBS_VN";                             //SDO orbital velocity in direction of increasing HGLT_OBS
      const char *RADIUS     = "R_SUN";                              //Image radius in pixels
      const char *EXPTIME    = "EXPTIME";                            //mean shutter open time
      const char *HCMNB0     = "HCMNB";
      const char *HCMWB0     = "HCMWB";
      const char *NXs        = "NX";
      const char *FSRNBs     = "FSRNB";
      const char *FSRWBs     = "FSRWB";
      const char *DATEs      = "DATE";

      int Mode	             = UNK_MODE;
      TIME   interntime;
      TIME   interntime2;
      TIME   interntime3;
      float  X0[nseq],Y0[nseq];                                      //location of Sun's center in pixels (starts 0.0)
      float  CDELT1[nseq];                                           //image scale in the x direction
      float  RSUN[nseq];
      int    FSN,NBADPERM,imcfg;
      double VELOCITY[nseq];                                         //Sun-SDO radial velocity
      double EXPOSURE[nseq];
      double OBS_VR, OBS_VW, OBS_VN, lam0, dlam_LOS[nx2*ny2];

      //initialization of Richard's code for de-rotation and gapfilling, and Jesper's code for rebinning
      //************************************************************************************************
      
      struct  initial const_param;                                   //structure containing the parameters for Richard's functions
 
      char dpath[]="/home/jsoc/cvs/Development/JSOC/";
      status = initialize_interpol(&const_param,&initfiles,nx,ny,dpath);
      if(status != 0)
	  DIE("Error: could not initialize the gapfilling routine\n");
      
      Mask = (unsigned char *)malloc(Nelem*sizeof(unsigned char));
      if(Mask == NULL)
	  DIE("Error: cannot allocate memory for Mask\n");

      struct fresize_struct fresizes;
      init_fresize_bin(&fresizes,ratio);

      
      /***********************************************************************************************************/
      /*READING, GAPFILLING, AND REBINNING THE FILTERGRAMS                                                       */
      /***********************************************************************************************************/

      int tuningint[nseq][2]; //HCM positions of the tuning waveplates
      int axisout22[2] = {nx,ny}; //column, row
 
      float  tempframe[ny2*nx2];
      double tuning[nseq][];
      float  *arrinL = NULL;       
      char  *ierror  = NULL;

      DRMS_Array_t *arrin     = NULL;
      DRMS_Array_t *Ierror    = NULL;
 
      Ierror = drms_array_create(typeEr,2,axisout22,NULL,&status);
      if(status != DRMS_SUCCESS || Ierror == NULL)
	  DIE("Error: could not create Ierror\n");


      //LOOP OVER ALL THE FILTERGRAMS OF THE DETUNE SEQUENCE
      printf("READING, GAPFILLING, AND REBINNING THE FILTERGRAMS in %dx%d\n",nx2,ny2);
      for(i=0;i<nseq;i++)                                            //read and rebin the filtergrams
	{
	  //READING IMAGE
	  segin     = drms_segment_lookupnum(rec[i],0);            //assume first record is first detune filtergram
	  arrin     = drms_segment_read(segin,type,&status);
	  if(status != DRMS_SUCCESS)
	      DIE("Error: unable to read a data segment\n");
	  arrinL = arrin->data;

	  //READING SOME KEYWORDS
	  FSN = drms_getkey_int(rec[i],"recnum",&status);
	  printf("FSN IMAGE = %d\n",FSN);
	  NBADPERM =  drms_getkey_int(rec[i],"NBADPERM",&status);
	  if(status != DRMS_SUCCESS) NBADPERM=-1;

	  tuningint[i][0] = drms_getkey_int(rec[i],"RLC_MTM1" ,&status); //tuning positions of the HCM of NB
	  if(status != DRMS_SUCCESS)
	      DIE1("Error: unable to read the keyword %s\n","MDI_MTM1");
	  else
              printf(" %d ",tuningint[i][0]);
	  tuningint[i][1] = drms_getkey_int(rec[i],"RLC_MTM2",&status); //tuning positions of the HCM of WB
	  if(status != DRMS_SUCCESS)
	      DIE1("Error: unable to read the keyword %s\n","MDI_MTM2");
	  else
              printf(" %d ",tuningint[i][1]);

          DPC = drms_getkey_int(rec[i], "DPC", &status);
          if (i==0)
            {
	    Mode = (DPC == 0x404e4400 ? CAL_MODE : OBS_MODE);
            printf(" %0X ",DPC);
            }
          if (Mode == OBS_MODE)
            {
	    X0[i]         = drms_getkey_float(rec[i],"MDI_X0",&status);
	    if(status != DRMS_SUCCESS || isnan(X0[i]) || X0[i] == -1)
	       fprintf(stderr,"Error: unable to read the keyword %s\n","MDI_X0");
            else
	       printf(" %f ",X0[i]);
	    Y0[i]         = drms_getkey_float(rec[i],CRPIX2,&status);
	    if(status != DRMS_SUCCESS || isnan(Y0[i]) || Y0[i] == -1)
	        fprintf(stderr,"Error: unable to read the keyword %s\n",CRPIX2);
            else
	       printf(" %f ",Y0[i]);
            }
          else
            {
	    X0[i]= 511.5;
	    Y0[i]= 511.5;
            }

	  CDELT1[i]  = drms_getkey_float(rec[i],SCALE,&status); //I don't think I really need CDELT...
	  if(status != DRMS_SUCCESS || isnan(CDELT1[i]))
	    {
	    fprintf(stderr,"Warning CDELT1 not found, use 0.5\n");
	    CDELT1[i]=0.5;
	    }
	  printf(" %f ",CDELT1[i]);

          if (Mode == OBS_MODE)
            {
	    RSUN[i]  = drms_getkey_float(rec[i],RADIUS,&status);
	    if(status != DRMS_SUCCESS || isnan(RSUN[i]))
	          DIE1("Error: unable to read the keyword %s\n",RADIUS);
            }
          else
	     RSUN[i]=511.5;
	  printf(" %f ",RSUN[i]);

	  VELOCITY[i]= drms_getkey_double(rec[i],VR,&status);
	  if(status != DRMS_SUCCESS || isnan(VELOCITY[i]))
	      DIE1("Error: unable to read the keyword %s\n",VR);
          printf(" %f ",VELOCITY[i]);

	  EXPOSURE[i]= drms_getkey_double(rec[i],EXPTIME,&status);
	  if(status != DRMS_SUCCESS || isnan(EXPOSURE[i]))
	      DIE1("Error: unable to read the keyword %s\n",EXPTIME);
	   printf(" %f ",EXPOSURE[i]);

        // Rebin the input array
        if (Mode == OBS_MODE)
          {
          // Need to rebin allowing missing data
          rebinArraySF(tempframe, arrinL, nx2, nx, RSUN[i], X0[i], Y0[i], 0.0);
          }
        else
          {
          //REBINNING, USING THE FUNCTION OF JESPER
          fresize(&fresizes,arrinL,tempframe,nx,ny,nx,nx2,ny2,nx2,0,0,0.0f);
          }

        //put the image contained in tempframe into the Inten array [IN DOUBLE PRECISION]
        for(jjj=0;jjj<ny2;++jjj)
           for(iii=0;iii<nx2;iii++)
              Inten[jjj][iii][i] = (double)tempframe[iii+jjj*nx2]; ///jjj=row,iii=column
        drms_free_array(arrin);
        }

      free_interpol(&const_param);

      //CREATE OUTPUT ARRAYS
      DRMS_Array_t *arrout   = NULL;                                 //array that will contain the phase maps AND Fe I LINE CHARACTERISTICS (LINEWIDTH AND LINEDEPTH)
      int axisout[3]         = {4,nx2,ny2};                          //size of the output arrays (nx2=number of columns, ny2=number of rows)
      type    = DRMS_TYPE_FLOAT;                                     //type of the output data: FLOAT
      arrout  = drms_array_create(type,3,axisout,NULL,&status);      //create the array that will contain the phase maps
      if(status != DRMS_SUCCESS || arrout == NULL)
	  DIE("Error: cannot create array arrout\n");

      DRMS_Array_t *arrout2  = NULL;                                 //array that will contain the reconstructed intensities 
      int axisout2[3]        = {nseq,nx2,ny2};                       //size of the output arrays (nx2=number of columns, ny2=number of rows)
      type    = DRMS_TYPE_DOUBLE; 
      arrout2 = drms_array_create(type,3,axisout2,NULL,&status); 
      if(status != DRMS_SUCCESS || arrout2 == NULL)
	  DIE("Error: cannot create array arrout\n");

      DRMS_Array_t *arrout3  = NULL;                                 //array that will contain the quality flag (did the code converge or not?)
      int axisout3[2]        = {nx2,ny2};                            //size of the output arrays (nx2=number of columns, ny2=number of rows)
      type    = DRMS_TYPE_SHORT; 
      arrout3 = drms_array_create(type,2,axisout3,NULL,&status); 
      if(status != DRMS_SUCCESS || arrout3 == NULL)
	  DIE("Error: cannot create array arrout\n");

      float  *Phig   = arrout->data;                                  //relative phases of the tunable elements
      double *IntenR = arrout2->data;                                 //reconstructed intensities
      short  *quality= arrout3->data;                                 //quality value
      memset(arrout->data,0.0,drms_array_size(arrout));               //initializes the phase values
      memset(arrout2->data,0.0,drms_array_size(arrout2));             //initializes the reconstructed intensities values
      memset(arrout3->data,0,drms_array_size(arrout3));               //initializes the quality flag values


/* MDI Keywords are:

 T_REC		1996.04.29_21:58:16_TAI – 2011.03.25_21:55:44_TAI
 ITYPE		"FD", "DARK", ""
 T_OBS		1996.04.29_21:58:16_TAI – 2011.03.25_21:55:44_TAI
 RECORD		int 0-19
 OBSMODE	"STD", "CAL", ""
 MDI_CAL1	0, 1, 2, 3, -1
 MDI_CAL2	0, 1, 3, -1
 MDI_MTM1	range 15 : 114, -1
 MDI_MTM2	range 41 : 169, -1
 RLC_MTM1	Relative to LC steps or -1
 RLC_MTM2	relative to LC steps or -1
 TRG_MTM1	e.g. "F1","F2","F3","F4","FC","LC-15",LC",LC+15","FC+1","F1+1","F2+1","F3+1","F4+1"
 TRG_MTM2	e.g. "F1","F2","F3","F4","FC","LC-15",LC",LC+15","FC+1","F1+1","F2+1","F3+1","F4+1"
Target and relative correlations:
MTM	 Fc	F1	F2	F3	F4	LC-15	LC	LC+15
MTM1	-45	-27	-9	9	27	-15	0	15
MTM2	-22	-13	-4	5	14	-15	0	15
NOTE MDI_MTMx = LC(x) + TRG_MTMx
 MDI_EXP	1 or 5 or -1
 OBSNAME	"cal_detune2b","cal_detune2a","obs_detune2b","obs_detune2a", "obs_detune1b","obs_detune1a"
		most are 2a or 2b or missing
 MDI_X0
 MDI_Y0
 CROTA2
 CAR_ROT
 CRLN_OBS
 CRLT_OBS
 DSUN_OBS
 OBS_VR
 OBS_VW
 OBS_VN
 RSUN_OBS
 EARTH_DT
 DATAMIN
 DATAMAX
 DATAMEDN
 DATAMEAN
 DATARMS
 DATASKEW
 DATAKURT
 DPC		ox404e4400 -> OBSMODE, 0x404d4400 -> CALMODE
 T_REC_epoch = 1977.01.01_00:33:12_TAI   ???
 T_REC_step = 1
*/


      //READ SOME KEYWORDS OF THE FIRST AND LAST FILTERGRAMS OF THE FULL DETUNE SEQUENCE
      Frecnum   = drms_getkey_int(rec2[0],"recnum",&status);          //value of the FSN of the 1st filtergram
      if(status != DRMS_SUCCESS)
	  DIE1("Error: unable to read keyword %s of the first filtergram\n","recnum");

      interntime = drms_getkey_time(rec2[0],"T_OBS",&status);        //value of the T_OBS of the 1st filtergram
      if(status != DRMS_SUCCESS)
	  DIE1("Error: unable to read keyword %s of the first filtergram\n","T_OBS");

      //we take the average value of the radial velocity (in m/s) of observer from 1st and last filtergrams
      OBS_VR = drms_getkey_double(rec2[0],VR,&status);            //value of disk center radial velocity for the 1st filtergram
      if(status != DRMS_SUCCESS)
	{
	  printf("Error: unable to read the keyword %s of the first filtergram\n",VR);
	  return(EXIT_FAILURE);
	}
      OBS_VR += drms_getkey_double(rec2[1],VR,&status);            //value of radial velocity for the last filtergram
      if(status != DRMS_SUCCESS)
	{
	  printf("Error: unable to read the keyword %s of the last filtergram\n",VR);
	  return(EXIT_FAILURE);
	}
      OBS_VR /= 2.0;

      Lrecnum  = drms_getkey_int(rec2[1],"recnum",&status);          //value of the FSN of the last filtergram
      if(status != DRMS_SUCCESS)
	  DIE1("Error: unable to read keyword %s of the last filtergram\n","recnum");
      interntime2= drms_getkey_time(rec2[1],"T_OBS",&status);        //value of the T_OBS of the last filtergram
      if(status != DRMS_SUCCESS)
	  DIE1("Error: unable to read keyword %s of the last filtergram\n","T_OBS");

      printf("recnum value of first filtergram %d\n",Frecnum);          //print value of the FSN of the 1st filtergram
      printf("recnum value of last filtergram %d\n",Lrecnum);          //print value of the FSN of the last filtergram


      //CONVERT THE POSITIONS OF TUNING HCM INTO PHASES (AND CHECK THAT ANGLES ARE 0, 120, OR 240 DEGREES)
      for(i=0;i<nseq;++i) 
	{
	  //FOR NB
	  tuning[i][0]  =  (double)((tuningint[i][0] *   8) % 360);
	  if(tuning[i][0] != 0.0 && tuning[i][0] != 120.0 && tuning[i][0] != 240.0 && tuning[i][0] != -120.0 && tuning[i][0] != -240.0)
	    {
	      printf("Error: detune sequence should have phases of 0, 120, or 240 degrees only: %f\n",tuning[i][0]);
	      return(EXIT_FAILURE);
	    }
	  tuning[i][0]  *=M_PI/180.0;

	  //FOR WB
	  tuning[i][1]  =  (double)((tuningint[i][1] *   8) % 360);
	  if(tuning[i][1] != 0.0 && tuning[i][1] != 120.0 && tuning[i][1] != 240.0 && tuning[i][1] != -120.0 && tuning[i][1] != -240.0)
	      {
		printf("Error: detune sequence should have phases of 0, 120, or 240 degrees only: %f\n",tuning[i][1]);
		return(EXIT_FAILURE);
	      }
	  tuning[i][1]  *=M_PI/180.0;

	}


      //FIND MINIMUM OF RSUN AMONG ALL THE IMAGES OF THE DETUNE SEQUENCE (NB: IN CAL_MODE, RSUN IS NOT THE SOLAR RADIUS)
      solarradiusmax = 10000.; //pixels
      indexref=0;
      for(k=0;k<nseq;k++) if(RSUN[k] < solarradiusmax)
	{
	  solarradiusmax = RSUN[k];
	  indexref=k;
	}
      // WARNING PHS - should this be solarradiusmax -= ratio/2 ?? XXXXXXX if problems change back to -= ratio or something...
      //solarradiusmax -= ratio/2;  //to avoid problems at the edge, because we average by groups of ratio pixels
      solarradiusmax -= ratio;
				// (IF WE DON'T DO THAT, THERE ARE A LOT OF NON-CONVERGENCE AT
                                // THE EDGE OF THE IMAGE DISK, AND THE CODE IS NOTICEABLY SLOWER)
      printf("Maximum radius (in pixels) at which phases are computed %f\n",solarradiusmax);

      //MAP OF RADIAL DISTANCES IN PIXELS
      printf("Pixel position of the image center used to compute the phase maps: %f %f\n",X0[indexref],Y0[indexref]);
      // WARNING PHS - This is where we need to decide if co-registration is required. XXXXXX for obsmode

      for(iii=0;iii<ny*nx;iii++)
	{
	  row    = iii / nx ;
	  column = iii % nx ;
	  distance2[iii]=sqrt( ((float)row-Y0[indexref])*((float)row-Y0[indexref])+((float)column-X0[indexref])*((float)column-X0[indexref]) ); 
	}

      if (Mode == OBS_MODE)
        {
        rebinArraySF(distance, distance2, nx2, nx, RSUN[indexref], X0[indexref], Y0[indexref], 512.0);
        }
      else
        {
        fresize(&fresizes,distance2,distance,nx,ny,nx,nx2,ny2,nx2,0,0,0.0f);
        free_fresize(&fresizes); 
        }



      /***********************************************************************************************************/
      /*BLOCKER FILTER + FRONT WINDOW                                                                            */
      /***********************************************************************************************************/
 
      mdi_filters(nlam, lam, blockerint);

      for(i=0;i<2;++i) FSR[i] = dpi/FSR[i];


      /***********************************************************************************************************/
      /*WAVELENGTH SHIFT DUE TO THE l.o.s. VELOCITY (DEPENDS ON THE LOCATION AT THE SOLAR SURFACE IN OBS_MODE)    */
      /* convention: + = away from Sun (redshift)                                                                */
      /***********************************************************************************************************/


      if (Mode == CAL_MODE)
        {
	for(jjj=0; jjj<ny2; ++jjj)
          for(iii=0; iii<nx2; ++iii)
	    dlam_LOS[jjj*nx2+iii] = dlamdv*OBS_VR;
        printf("VELOCITY = %f\n",OBS_VR);
        }
      else
        {
        int stat;
// int dims[2] = {nx2,nx2};
// DRMS_Array_t *arr = drms_array_create(DRMS_TYPE_FLOAT,2,dims,NULL,&stat);
// float *arrdata = arr->data;
        double DSUN_OBS = drms_getkey_double(rec2[0], "DSUN_OBS", &stat);
        double OBS_B0 = drms_getkey_double(rec2[0], "CRLT_OBS", &status); stat += status;
        double CROTA2 = drms_getkey_double(rec2[0], "CROTA2", &status); stat += status;
        OBS_VW = drms_getkey_double(rec2[0], VW, &status); stat += status;
        OBS_VN = (drms_getkey_double(rec2[0], VN, &status)); stat += status;
        if (stat || isnan(DSUN_OBS) || isnan(OBS_B0) || isnan(CROTA2) || isnan(OBS_VW) || isnan(OBS_VN))
          {
	  printf("one of CROTA2, CRLT_OBS, DSUN_OBS, OBS_VW, or OBS_VW keywords not available.");
	  return(EXIT_FAILURE);
          }
// printf("OBS_VR=%f, OBS_VW=%f, OBS_VN=%f, CROTA2=%f\n",OBS_VR,OBS_VW,OBS_VN,CROTA2);
        double X02 = X0[indexref]/ratio;
        double Y02 = Y0[indexref]/ratio;
        double rsun2 = RSUN[indexref]/ratio;
        double CDELT12 = CDELT1[0]*ratio;
        //l.o.s. velocity at pixel location [iii,jjj] 
	for(jjj=0; jjj<ny2; ++jjj)
          {
          for(iii=0; iii<nx2; ++iii)
            if (jjj<10 && iii==0)
{
	      dlam_LOS[jjj*nx2+iii] = dlamdv * OBS_VR; // for bullseye disk averages
// arrdata[jjj*nx2+iii] = 0.0;
}
            else
	      {
              double x = iii - X02; // binned pixels
              double y = jjj - Y02;
	      double V_sun = Solar_velocity(x,y,rsun2,OBS_B0,CROTA2,CDELT12,OBS_VR,OBS_VW,OBS_VN);
	      dlam_LOS[jjj*nx2+iii] = dlamdv * V_sun;
// arrdata[jjj*nx2+iii] = dlamdv * V_sun;
	      }
          }
// DRMS_RecordSet_t *ors = drms_create_records(drms_env, 1, "su_phil.dlam_LOS", DRMS_PERMANENT, &stat);
// DRMS_Record_t *orec = ors->records[0];
// DRMS_Segment_t *oseg = drms_segment_lookupnum(orec,0);
// drms_segment_write(oseg,arr,0);
// drms_close_records(ors,DRMS_INSERT_RECORD);
// return(0);
        }

// SPECIAL BULLS'EYE AVERAGING TEST
// Fill lower left part of rebinned array with averages from 10% to 100% in 10 steps or radius
{
for (jjj=0; jjj<10; jjj++)
  dlam_LOS[jjj*nx2+0] = dlamdv*OBS_VR;

for(i=0;i<nseq;++i) 
  {
  int ib;
  double bullseye[10];
  int nbullseye[10];
  for (ib=0; ib<10; ib++)
    {
    bullseye[ib] = 0.0;
    nbullseye[ib] = 0;
    }
  for(jjj=0;jjj<ny2;++jjj)
    for(iii=0;iii<nx2;iii++)
      {
      double val = Inten[jjj][iii][i];
      float dist = distance[jjj*nx2 + iii];
      if (!isnan(val))
        {
        for (ib=0; ib<10; ib++)
          {
          if (dist < round(solarradiusmax*(ib+1.0)/10.0))
	    {
            bullseye[ib] += val;
            nbullseye[ib] += 1;
            }
          }
        }
      }
  for (ib=0; ib<10; ib++)
    {
    Inten[ib][0][i] = (nbullseye[ib] ? bullseye[ib]/nbullseye[ib] : 0.0);
    }
  }
}
// END BULLSEYE GENERATION

      
      //OPENING OF BINARY FILE CONTAINING CONTRASTS OF TUNABLE ELEMENTS

XXXXXXX use CON1 and CON2 from mdi_param.c

      
      gsl_vector *Residual = NULL;
      gsl_matrix *Jac      = NULL;
      gsl_matrix *Weights  = NULL;
      gsl_matrix *VV       = NULL;
      gsl_vector *SS       = NULL;
      gsl_vector *work     = NULL;
      gsl_vector *tempvec  = NULL;
      gsl_vector *tempvec2 = NULL;
      
      int    location,location1[8];

      /***********************************************************************************************************/
      /*LOOP OVER ALL THE PIXELS OF THE CCD                                                                      */
      /***********************************************************************************************************/      

      printf("START LOOP\n");
      printf("NB: IF THERE ARE MANY NON CONVERGENCE, REMEMBER TO CHECK THE VALUE OF thresh IN HMIparam.h\n");

#pragma omp parallel default(none) shared(phaseguess,solarradiusmax,factor,Inten,axisout,axisout2,depth0,thresh,width0,nx2,ny2,FSR,CON1,CON2,dpi,lam,tuning,dlam,Phig,distance,blockerint,IntenR,quality,iterstop,shiftw,Icg,fdepthg,fwidthg,a, lam0, dlam_LOS, dlamdv)   private(location,profilef,residual,Residual,i,j,iii,jjj,k,converg,compteur,line,templine,dlinedI0,dlinedw,tempval,profileg,dIntendI0,dIntendw,dIntendIc,dIntendPhi0,dIntendPhi1,dIntendPhi2,tempres,history,err,minimum,Jac,Weights,VV,SS,work,tempvec,tempvec2,jose,phaseE2,phaseE3,phaseE4,phaseE5,contrasts,l,linea,lineb,linec,lined,gaussian,dlinedIc)
      {

	Residual = gsl_vector_alloc(nseq);
	Jac      = gsl_matrix_alloc(nseq,nparam);
	Weights  = gsl_matrix_alloc(nparam,nparam);
	VV       = gsl_matrix_alloc(nparam,nparam);
	SS       = gsl_vector_alloc(nparam);
	work     = gsl_vector_alloc(nparam);
	tempvec  = gsl_vector_alloc(nparam);
	tempvec2 = gsl_vector_alloc(nparam);
	
#pragma omp for
	for(jjj=0;jjj<ny2;jjj++) ///jjj=row
	  {
	    printf("%d/%d\n",jjj,ny2);
	    for(iii=0;iii<nx2;iii++) ///iii=column
	      {	   
 
		// if(distance[jjj*nx2+iii] <= solarradiusmax) // Changed for BULLSEYE test
		if(distance[jjj*nx2+iii] <= solarradiusmax || (jjj<10 && iii==0))
		  {
                  lam0 = dlam_LOS[jjj*nx2+iii];
		    // for(j=0;j<nseq;j++) printf(" %f ",Inten[jjj][iii][j]);
		    converg = 0;
		    compteur= 0;
		    
		    Icg[jjj][iii] = thresh;       //estimate of the solar continuum
		    fdepthg[jjj][iii] = depth0*thresh;//depth of the solar line according to Stenflo & Lindegren (1977)
		    fwidthg[jjj][iii] = width0;       //width of the solar line according to Stenflo & Lindegren (1977)
		    
		    Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]] = (float)((double)phaseguess[0]*M_PI/180./FSR[0]); //guess values for the phases; NB
		    Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]] = (float)((double)phaseguess[1]*M_PI/180./FSR[1]); // WB
		    
		    //NON-TUNABLE HMI TRANSMISSION PROFILE
		    location    =iii+jjj*nx2;

		    //iii = column, jjj = row


		    while(converg == 0)
		      {
			for (k=0;k<nlam;k++)  //BUILD THE SOLAR LINE PROFILE, ASSUMING TWO GAUSSIAN PROFILES
			  {

			    gaussian[k]= -0.0074*exp(-(lam[k]-lam0+0.200)*(lam[k]-lam0+0.200)/0.13/0.13)-0.021*exp(-(lam[k]-lam0-0.05)*(lam[k]-lam0-0.05)/0.18/0.18); //calibration 13
XXXXXXXXXXX line profile here, shifted for OBS_VR+V

			    dlinedIc[k]= 1.0-gaussian[k];
			    l=(lam[k]-lam0)/fwidthg[jjj][iii];
			    //if(fabs(lam[k]-lam0) > 1.5) //no need to calculate far from the line center, because it's NaN 
			    if(fabs(l) > 26.5) //sinh((x>26.5)^2) = NAN
			      {
				line[k]=Icg[jjj][iii]-gaussian[k]*Icg[jjj][iii];
				dlinedw[k] =0.0;
				dlinedI0[k] =0.0;
			      }
			    else
			      {
				line[k]  = Icg[jjj][iii] - fdepthg[jjj][iii] * exp(-l*l) * (1.0-a*2.0/sqrt(M_PI)*(1.0/2.0/l/l)*((4.0*l*l+3.0)*(l*l+1.0)*exp(-l*l)-1.0/l/l*(2.0*l*l+3.0)*sinh(l*l)))-gaussian[k]*Icg[jjj][iii];
				
				//for the derivative relative to fwidthg
				l=(lam[k]-lam0)/(fwidthg[jjj][iii]+0.0001);
				linea[k]= Icg[jjj][iii] - fdepthg[jjj][iii]*exp(-l*l)*(1.0-a*2.0/sqrt(M_PI)*(1.0/2.0/l/l)*((4.0*l*l+3.0)*(l*l+1.0)*exp(-l*l)-1.0/l/l*(2.0*l*l+3.0)*sinh(l*l)))-gaussian[k]*Icg[jjj][iii];
				l=(lam[k]-lam0)/(fwidthg[jjj][iii]-0.0001);
				lineb[k]= Icg[jjj][iii]-fdepthg[jjj][iii]*exp(-l*l)*(1.0-a*2.0/sqrt(M_PI)*(1.0/2.0/l/l)*((4.0*l*l+3.0)*(l*l+1.0)*exp(-l*l)-1.0/l/l*(2.0*l*l+3.0)*sinh(l*l)))-gaussian[k]*Icg[jjj][iii]; 
				dlinedw[k] = (linea[k]-lineb[k])/(2.0*0.0001);//derivative relative to fwidthg
				
				//for the derivative relative to fdepthg
				l=(lam[k]-lam0)/fwidthg[jjj][iii];
				linec[k]=Icg[jjj][iii]-(fdepthg[jjj][iii]+0.001)*exp(-l*l)*(1.0-a*2.0/sqrt(M_PI)*(1.0/2.0/l/l)*((4.0*l*l+3.0)*(l*l+1.0)*exp(-l*l)-1.0/l/l*(2.0*l*l+3.0)*sinh(l*l)))-gaussian[k]*Icg[jjj][iii];
				lined[k]=Icg[jjj][iii]-(fdepthg[jjj][iii]-0.001)*exp(-l*l)*(1.0-a*2.0/sqrt(M_PI)*(1.0/2.0/l/l)*((4.0*l*l+3.0)*(l*l+1.0)*exp(-l*l)-1.0/l/l*(2.0*l*l+3.0)*sinh(l*l)))-gaussian[k]*Icg[jjj][iii];
				dlinedI0[k] = (linec[k]-lined[k])/(2.0*0.001);//derivative relative to fdepthg
			      }				
			  }
for (k=0;k<nlam;k++)  //CHECK THE SOLAR LINE PROFILE, etc.
{
int nansfound=0;
if (isnan(blockerint[k])) {nansfound++;printf("nan found in blockerint, k=%d\n",k);}
if (isnan(lam[k])) {nansfound++;printf("nan found in lam, k=%d\n",k);}
if (isnan(linea[k])) {nansfound++;printf("nan found in linea, k=%d\n",k);}
if (isnan(lineb[k])) {nansfound++;printf("nan found in lineb, k=%d\n",k);}
if (isnan(linec[k])) {nansfound++;printf("nan found in linec, k=%d\n",k);}
if (isnan(lined[k])) {nansfound++;printf("nan found in lined, k=%d\n",k);}
if (isnan(profilef[k])) {nansfound++;printf("nan found in profilef, k=%d\n",k);}
if (isnan(gaussian[k])) {nansfound++;printf("nan found in gaussian\n");}
if (nansfound) {printf("jjj=%d, iii=%d, distance=%f\n",jjj,iii,distance[jjj*nx2+iii]);exit(1);}
}

			for(j=0;j<nseq;j++)
			  {
			    residual      = 0.0;
			    dIntendI0[j]  = 0.0;
			    dIntendw [j]  = 0.0;
			    dIntendIc[j]  = 0.0;
			    dIntendPhi0[j]= 0.0;
			    dIntendPhi1[j]= 0.0;
			    
			    for (k=0;k<nlam;k++) 
			      {
				profileg[k]    = 0.125 * profilef[k] *
                                  (1.0+CON1*cos(FSR[0]*(lam[k]+(double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][0])) *
				  (1.0+CON2*cos(FSR[1]*(lam[k]+(double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][1])); 
if (isnan(profileg[k]))
{
printf("nan found in profileg, k=%d, jjj=%d, iii=%d, seq=%d\n",k,jjj,iii,j);
printf("contrasts={%f,%f}, tuning[%d]={%f, %f}\n",contrasts[0],contrasts[1],j,tuning[j][0],tuning[j][1]);
printf("Phig[this jjj,iii]=%f, lam[k]=%f, FSR[0-1]={%f, %f}, phaseguess[0-1]={%f, %f}\n",Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]],lam[k],FSR[0],FSR[1],phaseguess[0],phaseguess[1]);
exit(1);
}

				residual      += line[k]*profileg[k]*dlam;
				dIntendI0[j]  += (dlinedI0[k]*profileg[k])*dlam;
				dIntendw [j]  += (dlinedw[k] *profileg[k])*dlam;
				dIntendIc[j]  += (dlinedIc[k]*profileg[k])*dlam;
				dIntendPhi0[j]+= (line[k]*(-CON1*FSR[0]*sin(FSR[0]*(lam[k]+(double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][0])) * 
					(1.0+CON2*cos(FSR[1]*(lam[k]+(double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][1])) / 8.0*profilef[k])*dlam;
				dIntendPhi1[j]+= (line[k]*(-CON2*FSR[1]*sin(FSR[1]*(lam[k]+(double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][1])) *
					(1.0+contrasts[0]*cos(FSR[0]*(lam[k]+(double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][0])) / 8.0*profilef[k])*dlam;
			      }
			    jose=Inten[jjj][iii][j]/factor-residual;
			    //printf("%d %d %f %f %f %f %f\n",iii,jjj,dIntendI0[j],dIntendIc[j],dIntendw[j],dIntendPhi0[j],dIntendPhi1[j]);
			    gsl_vector_set(Residual,j,jose);                               //residual is computed	  
			    gsl_matrix_set(Jac,(size_t)j,0,dIntendI0[j]);                  //Jacobian matrix
			    gsl_matrix_set(Jac,(size_t)j,1,dIntendIc[j]);
			    gsl_matrix_set(Jac,(size_t)j,2,dIntendw[j]);
			    gsl_matrix_set(Jac,(size_t)j,3,dIntendPhi0[j]);
			    gsl_matrix_set(Jac,(size_t)j,4,dIntendPhi1[j]);
			    gsl_matrix_set(Jac,(size_t)j,5,dIntendPhi2[j]);
			    
			  }//for(j=0;j<nseq;++j)

// printf("D\n");
			//SVD algorithm from the GNU Scientific Library
			//A is MxN. A = U S V^T for M >= N. On output the matrix A is replaced by U.
			gsl_linalg_SV_decomp(Jac,VV,SS,work);
			gsl_matrix_set_zero(Weights);
			for(i=0;i<nparam;i++) gsl_matrix_set(Weights,(size_t)i,(size_t)i,1.0/gsl_vector_get(SS,(size_t)i)); //creates the diagonal matrix
			gsl_blas_dgemv(CblasTrans  ,1.0,Jac,Residual   ,0.0,tempvec );
			gsl_blas_dgemv(CblasNoTrans,1.0,Weights,tempvec,0.0,tempvec2);
			gsl_blas_dgemv(CblasNoTrans,1.0,VV,tempvec2    ,0.0,tempvec );
			
			//compute relative changes in the parameters
			tempres[0]   = fabs(gsl_vector_get(tempvec,0)/fdepthg[jjj][iii]);
			tempres[1]   = fabs(gsl_vector_get(tempvec,1)/Icg[jjj][iii]);
			tempres[2]   = fabs(gsl_vector_get(tempvec,2)/fwidthg[jjj][iii]);
			tempres[3]   = fabs(gsl_vector_get(tempvec,3)/(double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]);
			tempres[4]   = fabs(gsl_vector_get(tempvec,4)/(double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]);
			

// printf("D0: tempres[x], x0: %g, X1: %g, X2: %g, X3: %g, X4: %g\n",tempres[0],tempres[1],tempres[2],tempres[3],tempres[4]);
			err[compteur]= fmax(tempres[0],fmax(tempres[1],fmax(tempres[2],fmax(tempres[3],fmax(tempres[4])))));
					      
// printf("D1: compteur=%d, err=%g, iterstop=%g\n",compteur,err[compteur],iterstop);
			if( (compteur == maxsteps-1) && (err[compteur] > iterstop) ) converg = 2; //no convergence 
			if(err[compteur] <= iterstop)                                converg = 1; //convergence (we stop the iteration when the maximum relative change is less than iterstop)
// printf("E, converg=%d\n",converg);
			
			if(converg == 2)
			  {
			    gsl_sort_smallest_index(&minimum,1,err,1,maxsteps);                //finds the minimum of err
			    fdepthg[jjj][iii] = history[minimum][0];
			    Icg[jjj][iii]     = history[minimum][1];
			    fwidthg[jjj][iii] = history[minimum][2];
			    Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]] = (float)history[minimum][3];
			    Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]] = (float)history[minimum][4];
			    quality[location]=1;                                               //the quality flag is raised to mention that the code did not converge
			  }
			
			if(converg != 2)
			  {
			    fdepthg[jjj][iii] += gsl_vector_get(tempvec,0);
			    if( (fdepthg[jjj][iii] < (0.2*thresh*depth0)) || (fdepthg[jjj][iii] > (3.0*depth0*thresh)) || isnan(fdepthg[jjj][iii])) fdepthg[jjj][iii] = depth0*thresh;
			    Icg[jjj][iii]     += gsl_vector_get(tempvec,1);
			    if( (Icg[jjj][iii] < (0.2*thresh)) || (Icg[jjj][iii] > (3.0*thresh)) || isnan(Icg[jjj][iii])) Icg[jjj][iii] = thresh;
			    fwidthg[jjj][iii] += gsl_vector_get(tempvec,2);
			    if( (fwidthg[jjj][iii] > (5.0*width0)) || (fwidthg[jjj][iii] < (0.3*width0)) || isnan(fwidthg[jjj][iii])) fwidthg[jjj][iii] = width0;
			    Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]] += (float)gsl_vector_get(tempvec,3);
			    if( fabsf(Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]) > (1.2*M_PI/FSR[0]) ) Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]] = 0.;
			    Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]] += (float)gsl_vector_get(tempvec,4);
			    if( fabsf(Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]) > (1.2*M_PI/FSR[1]) ) Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]] = 0.;
			  }
// printf("G, converg=%d\n",converg);
			
			history[compteur][0]=fdepthg[jjj][iii];
			history[compteur][1]=Icg[jjj][iii];
			history[compteur][2]=fwidthg[jjj][iii];
			history[compteur][3]=(double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]];
			history[compteur][4]=(double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]];
					   
			compteur += 1;
			
			if(converg == 2)
                          {
			  // for (j=maxsteps-4; j<maxsteps;j++)
			  j = maxsteps-1;
                          printf("no convergence %d,%d: %d: fd=%f, Ic=%f, fw=%f, NB=%f, WB=%f, E1=%f\n",jjj,iii,j,history[j][0],history[j][1],history[j][2],history[j][3]*FSR[0]*Deg2Rad,history[j][4]*FSR[1]*Deg2Rad,history[j][5]*FSR[2]*Deg2Rad);
                          }
		      }//while(converg==0)
		


		    /***********************************************************************************************************/
		    /*RECONSTRUCTING THE DETUNE SEQUENCE                                                                       */
		    /***********************************************************************************************************/       
  
		    for(k=0;k<nlam;++k)
		      {
			l=(lam[k]-lam0)/fwidthg[jjj][iii];
			//if(fabs(lam[k]-lam0) <= 1.5) //no need to calculate far from the line center
			if(fabs(l) <= 26.5)
			  {
			    line[k]  = Icg[jjj][iii]-fdepthg[jjj][iii]*exp(-l*l)*(1.0-a*2.0/sqrt(M_PI)*(1.0/2.0/l/l)*((4.0*l*l+3.0)*(l*l+1.0)*exp(-l*l)-1.0/l/l*(2.0*l*l+3.0)*sinh(l*l)))-gaussian[k]*Icg[jjj][iii];			    
			  }
			else line[k]=Icg[jjj][iii]-gaussian[k]*Icg[jjj][iii];
		      }
		    for(j=0;j<nseq;j++)
		      {
			residual   = 0.0;
			for(k=0;k<nlam;++k)
			  {
			    profileg[k]=profilef[k]*0.125 *
                              (1.0+contrasts[0]*cos(FSR[0]*(lam[k]+(double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][0])) *
                              (1.0+contrasts[1]*cos(FSR[1]*(lam[k]+(double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]])+tuning[j][1]));
			    residual      += line[k]*profileg[k]*dlam;
			  }
			IntenR[j+iii*axisout2[0]+jjj*axisout2[0]*axisout2[1]] = residual;
		      }
		    
		    Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=(float)((double)Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]*(FSR[0]*180.0/M_PI));
		    Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=(float)((double)Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]*(FSR[1]*180.0/M_PI));
		    Phig[2+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=fwidthg[jjj][iii];//we save sigma in Angstroms
		    Phig[3+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=fdepthg[jjj][iii]/Icg[jjj][iii];//we save the linedepth, for a continuum of 1
		  }//if(distance <= solarradiusmax)
		else
		  {
		    Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0; //NB Michelson
		    Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0; //WB Michelson
		    Icg[jjj][iii]=0.0;
		    fwidthg[jjj][iii]=0.0;
		    fdepthg[jjj][iii]=0.0;
		    Phig[2+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0;
		    Phig[3+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0;
		  }
	      }//for jjj
	  }//for iii
	
	
	gsl_vector_free(Residual);
	gsl_matrix_free(Jac);
	gsl_matrix_free(VV);
	gsl_vector_free(SS);
	gsl_vector_free(work);
	gsl_vector_free(tempvec);
	gsl_vector_free(tempvec2);
	gsl_matrix_free(Weights);	  
      }//pragma omp parallel    
      

      /***********************************************************************************************************/
      /*COMPUTE THE SPATIAL AVERAGES OF THE DETUNE SEQUENCES (ORIGINAL AND RECONSTRUCTED)                        */
      /***********************************************************************************************************/       
      
      
      double meanInten[nseq], meanIntenR[nseq];
      double NBphase=0.0, WBphase=0.0, E1phase=0.0;
      double meanIcg=0.0,meanwidthg=0.0,meandepthg=0.0;
      compteur = 0;

      for(j=0;j<nseq;++j)
	{
	  meanInten[j] =0.0;
	  meanIntenR[j]=0.0;
	}

      for(jjj=0;jjj<ny2;jjj++) ///jjj=row
	{
	  for(iii=0;iii<nx2;iii++) ///iii=column
	    {	    
	      //if(distance[jjj*nx2+iii] <= (solarradiustable/CDELT1[indexref]) && quality[iii+jjj*nx2] == 0)
	      if(distance[jjj*nx2+iii] <= (900.0/CDELT1[indexref]) && quality[iii+jjj*nx2] == 0) //900 arcsec instead of the full disk 976"
		{
		  NBphase += Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]];
		  WBphase += Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]];
		  meanIcg += Icg[jjj][iii];
		  meanwidthg += fwidthg[jjj][iii]; 
		  meandepthg += fdepthg[jjj][iii]; 
		  compteur+= 1;
		  for(j=0;j<nseq;++j)
		    {
		      // meanInten[j] +=Inten[jjj][iii][j]/factor;                             //original intensities
		      meanInten[j] +=Inten[jjj][iii][j]/nseq;                                  //original intensities
		      meanIntenR[j]+=factor/nseq*IntenR[j+iii*axisout2[0]+jjj*axisout2[0]*axisout2[1]];    //reconstructed intensities
		    }
		}
	    }
	}
      
      NBphase = NBphase/(double)compteur;
      WBphase = WBphase/(double)compteur;
      meanIcg = meanIcg/(double)compteur;
      meanwidthg = meanwidthg/(double)compteur;
      meandepthg = meandepthg/(double)compteur;


      drms_free_array(arrout2);



      /***********************************************************************************************************/
      /*ATTEMPT TO CORRECT THE NON-CONVERGENCE                                                                   */
      /***********************************************************************************************************/      

      printf("CORRECTION OF PIXELS WHERE THERE WAS NO CONVERGENCE, USING NEIGHBOR AVERAGING\n");
      for(jjj=0;jjj<ny2;jjj++) ///jjj=row
	{
	  for(iii=0;iii<nx2;iii++) ///iii=column
	    {
	      location=iii+jjj*nx2;
	      if(quality[location] ==1) //there was no convergence
		{
		  //location of 8 neighboring pixels
		  location1[0]=(jjj+1)*nx2+iii;     //(jjj+1,iii) (row,column)
		  location1[1]=(jjj+1)*nx2+(iii+1); //(jjj+1,iii+1)
		  location1[2]=jjj*nx2+(iii+1);     //(jjj,iii+1)
		  location1[3]=(jjj-1)*nx2+(iii+1); //(jjj-1,iii+1)
		  location1[4]=(jjj-1)*nx2+iii;     //(jjj-1,iii)
		  location1[5]=(jjj-1)*nx2+(iii-1); //(jjj-1,iii-1)
		  location1[6]=jjj*nx2+(iii-1);     //(jjj,iii-1)
		  location1[7]=(jjj+1)*nx2+(iii-1); //(jjj+1,iii-1)
		  compteur=0;
		  Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0;
		  Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0;
		  Phig[2+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0;
		  Phig[3+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=0.0;
		  for(j=0;j<8;++j)
		    {
		      if(quality[location1[j]] == 0) //there was convergence
			{
			  Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]+=Phig[0+axisout[0]*(location1[j])];
			  Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]+=Phig[1+axisout[0]*(location1[j])];
			  Phig[2+iii*axisout[0]+jjj*axisout[0]*axisout[1]]+=Phig[2+axisout[0]*(location1[j])];
			  Phig[3+iii*axisout[0]+jjj*axisout[0]*axisout[1]]+=Phig[3+axisout[0]*(location1[j])];
			  compteur+=1;
			}		  
		    }
		  if(compteur != 0)
		    {
		      Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=Phig[0+iii*axisout[0]+jjj*axisout[0]*axisout[1]]/(float)compteur;
		      Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=Phig[1+iii*axisout[0]+jjj*axisout[0]*axisout[1]]/(float)compteur;
		      Phig[2+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=Phig[2+iii*axisout[0]+jjj*axisout[0]*axisout[1]]/(float)compteur;
		      Phig[3+iii*axisout[0]+jjj*axisout[0]*axisout[1]]=Phig[3+iii*axisout[0]+jjj*axisout[0]*axisout[1]]/(float)compteur;
		    }
		}
	    }
	}

      /***********************************************************************************************************/
      /*DISPLAY RESULTS                                                                                          */
      /***********************************************************************************************************/      
 
      printf("RESULTS\n");
      printf("-----------------------------------------------------------------------------\n");
      printf("\n SPATIALLY AVERAGED INTENSITIES \n");
      for(j=0;j<nseq;++j) printf("%f %f\n",meanInten[j],meanIntenR[j]);
      printf("\n SPATIALLY AVERAGED PHASES (IN DEGREES)\n");
      printf("NB Michelson: %f\n",NBphase);
      printf("WB Michelson: %f\n",WBphase);
      printf("\n SPATIALLY AVERAGED WIDTH, DEPTH, AND CONTINUUM\n");
      printf("WIDTH: %f \n",meanwidthg);
      printf("DEPTH: %f \n",meandepthg);
      printf("CONTINUUM: %f \n",meanIcg);
      printf("\n COTUNE TABLE \n");
      printf("WB    NB\n");
      int    table[2][20];
      
      cotunetable(NBphase,WBphase,table);
      
      for(j = 0; j < 20; ++j)
      {
        printf("%d   %d\n",table[0][j],table[1][j]);      
      } 

      /***********************************************************************************************************/
      /*CREATE A RECORD FOR THE PHASE MAPS                                                                       */
      /***********************************************************************************************************/      
      
      DRMS_RecordSet_t *dataout = NULL;
      DRMS_Record_t  *recout    = NULL;
      DRMS_Segment_t *segout    = NULL;

      dataout = drms_create_records(drms_env,1,dsout,DRMS_PERMANENT,&status);
      
      if (status != DRMS_SUCCESS)
	{
	  printf("Could not create a record for the phase maps\n");
	  return(EXIT_FAILURE);
	}
      if (status == DRMS_SUCCESS)
	{
	printf("Writing a record on the DRMS for the phase maps\n");
	recout = dataout->records[0];

	  //WRITE KEYWORDS
	  
XXXXXXXXXXXXX fix keywords
	status += drms_setkey_int(recout,"START_rec",Frecnum);              //FSN_START
	status += drms_setkey_time(recout,"T_START",interntime);          //T_START
	status += drms_setkey_int(recout,"STOP_rec",Lrecnum); 
	status += drms_setkey_time(recout,"T_STOP",interntime2);          //T_STOP
	Mrecnum= (Frecnum+Lrecnum)/2;                           //FSN_REC (PRIME KEY)
	interntime3 = (interntime+interntime2)/2.0;                  //T_REC   (PRIME KEY)
	status += drms_setkey_int(recout,"RECNUM",keyvalue3);
	status += drms_setkey_time(recout,"T_REC",interntime3);
	status += drms_setkey_string(recout,COMMENTS,COMMENT);        //COMMENT ON THE DATA 
	status += drms_setkey_float(recout,CRPIX1,X0[indexref]+1.0);  //CRPIX1
	status += drms_setkey_float(recout,CRPIX2,Y0[indexref]+1.0);  //CRPIX2
	status += drms_setkey_float(recout,SCALE,CDELT1[indexref]);   //CDELT1
	status += drms_setkey_int(recout,FOCUS,HCFTID[indexref]);     //HCFTID
	status += drms_setkey_float(recout,RADIUS,solarradiusmax);    //R_SUN (NOT NECESSARILY SOLAR RADIUS, BUT SIZE OF THE IMAGE. solarradiusmax differs from RSUN[indexref] by -ratio pixels)
	status += drms_setkey_int(recout,HCMNB0,(int)round(-NBphase/6.0)+60); // NB motor position for average phase
	status += drms_setkey_int(recout,HCMWB0,(int)round(-WBphase/6.0)+60); // WB motor position for average phase
	status += drms_setkey_int(recout,NXs,nx2);
	char  DATEOBS[256];
	sprint_time(DATEOBS,CURRENT_SYSTEM_TIME,"UTC",1);
	status += drms_setkey_string(recout,DATEs,DATEOBS);           //DATE AT WHICH THE FILE WAS CREATED
	status += drms_setkey_float(recout,FSRNBs,2.0*M_PI/FSR[0]);
	status += drms_setkey_float(recout,FSRWBs,2.0*M_PI/FSR[1]);
	status += drms_setkey_float(recout,CBLOCKERs,centerblocker2);
        status += drms_setkey_float(recout, "PHASENBM", NBphase);
        status += drms_setkey_float(recout, "PHASEWBM", WBphase);
        status += drms_setkey_float(recout, "LW_MEAN", meanwidthg);
        status += drms_setkey_float(recout, "LD_MEAN", meandepthg);
        status += drms_setkey_float(recout, "CON_MEAN", meanIcg);
        status += drms_setkey_string(recout,"MODE",(Mode==OBS_MODE ? "OBS_MODE" : "CAL_MODE"));
	if (status)
          {
          printf("Failed to write one or more keywords.");
          return(EXIT_FAILURE);
          }

	  //WRITE DATA SEGMENTS
	  segout = drms_segment_lookupnum(recout, 0);
	  drms_segment_write(segout, arrout, 0);                       //write the file containing the phase maps
	  segout = drms_segment_lookupnum(recout, 1);
	  drms_segment_write(segout, arrout3,0);                       //write the file containing the quality flag
          printf("Phases and quality segments written\n");
          // Write phasemaps as separate files
          char *segnames[] = {"NB", "WB", "LW", "LD"};
          for (k=0; k<4; k++)
            {
// printf("Segment %s to be written.\n",segnames[k]);
            int status;
            int axis[2]         = {ny2,nx2};           //size of the output arrays (nx2=number of columns, ny2=number of rows)
            DRMS_Array_t *arr = drms_array_create(DRMS_TYPE_FLOAT,2,axis,NULL,&status);  //create the array that will contain the phase maps
            float *data = arr->data;
            for (jjj=0; jjj<ny2; jjj++)
              for (iii=0; iii<nx2; iii++)
		data[jjj*nx2 + iii] = Phig[k+iii*axisout[0]+jjj*axisout[0]*axisout[1]];
	    segout = drms_segment_lookup(recout, segnames[k]);
if (!segout){printf("lookup for %s segment failed\n", segnames[k]);continue;}
	    status = drms_segment_write(segout, arr,0);                       //write the file containing a phasemap or value
	    printf("Segment for %s %s written.\n",segnames[k],(status? "FAILED TO BE" : ""));
            drms_free_array(arr);
            }

	  //CLOSE THE RECORDS
	  drms_close_records(dataout, DRMS_INSERT_RECORD);             //insert the record in DRMS
	}      
      drms_free_array(arrout);   
      drms_free_array(arrout3); 
    }
    
  return error;
  
}


// Functions needed for OBS_MODE LOS velocity

double Solar_velocity(double x,double y,double rsun2,double OBS_B0,double CROTA2, double CDELT, double obs_vr, double obs_vw, double obs_vn)
  {
  double Msun = 200.0;   // Mean velocity
  double A = 1995.0;
  double B = -315.0;
  double LS0 = -455.0;
  double LS1 = 140.0;
  double LS2 = 720.0;
  OBS_B0 *= Deg2Rad;
  CROTA2 *= Deg2Rad;
  double arcsec2Rad = Deg2Rad/3600;
  double P = -CROTA2;
  double cosa = cos(CROTA2);
  double sina = sin(CROTA2);
  double X = x*cosa - y*sina; // pixels
  double Y = y*cosa + x*sina;
  double r = sqrt(x*x + y*y)/rsun2;
  double lon, sinlat, coslat, sig;
  double sinB2, rfac, r2, z, val;
  if (r < 1.0)
      img2sphere(X/rsun2, Y/rsun2, rsun2*CDELT*arcsec2Rad, OBS_B0, 0.0, P, NULL, NULL, &lon, &sinlat, &coslat, &sig, NULL, NULL);
  else
    return(0.0);
  // away from Sun is positive, away from East is positive, away from South is positive, red shift is positive.
  double obs_v = obs_vr*cos(sig) - obs_vw*sin(X*CDELT*arcsec2Rad) - obs_vn*sin(Y*CDELT*arcsec2Rad);
  sinB2 = sinlat * sinlat;
  rfac = coslat * cos(OBS_B0) * sin(lon);
  r2 = r * r;
  z = 1 - sqrt(1.0 - r2); // 1 - costheta
//printf("x=%f, y=%f, lon=%f, lat=%f\n",X/rsun2,Y/rsun2,lon/Deg2Rad,asin(sinlat)/Deg2Rad);
  // val = Msun + rfac * (A + B* sinB2 * (1.0 + sinB2)) + z * (LS0 + z * (LS1 + z * LS2)) ;
  // val = obs_vr;
  // val = obs_v;
  // val = obs_vr + Msun + rfac * A;
  val = obs_v + Msun + rfac * (A + B* sinB2 * (1.0 + sinB2)) + z * (LS0 + z * (LS1 + z * LS2)) ;
  return(val);
  }

/*
 *                     FUNCTIONS
 *
 *          img2sphere()
 *
 *  Written: Rick Bogart
 *
 */

int img2sphere (double x, double y, double ang_r, double latc, double lonc,
    double pa, double *p_rho, double *p_lat, double *p_lon, double *p_sinlat,
    double *p_coslat, double *p_sig, double *p_mu, double *p_chi)
  {
  /*
   *  Taken from cartography.c, written by Rick Bogart
   *  Map projected coordinates (x, y) to (lon, lat) and (rho | sig, chi)
   *
   *  Arguments:
   *    x            Plate locations, in units of the image radius, relative
   *    y                to the image center
   *                      These x,y are 0:1 not arcsecs and not rotated by crota2.
   *    ang_r         Apparent semi-diameter of sun (angular radius of sun at
   *                      the observer's tangent line)
   *    latc          Latitude of disc center, uncorrected for light travel time
   *    lonc          Longitude of disc center
   *    pa            Position angle of solar north on image, measured eastward
   *                      from north (sky coordinates)
   *  Return values:
   *    rho           Angle point:sun_center:observer
   *    lon           Heliographic longitude
   *    lat           Heliographic latitude
   *    sinlat        sine of heliographic latitude
   *    coslat        cosine of heliographic latitude
   *    sig           Angle point:observer:sun_center
   *    mu            cosine of angle between the point:observer line and the
   *                      local normal
   *    chi           Position angle on image measured westward from solar
   *                      north
   *
   *  All angles are in radians.
   *  Return value is 1 if point is outside solar radius (in which case the
   *    heliographic coordinates and mu are meaningless), 0 otherwise.
   *  It is assumed that the image is direct; the x or y coordinates require a
   *    sign change if the image is inverted.
   *
   * modified to ignore setting return values for NULL pointers.
   *
   */

  double w_rho, w_lat, w_lon, w_sinlat, w_coslat, w_sig, w_mu, w_chi;
  double *rho= &w_rho, *lat= &w_lat, *lon= &w_lon, *sinlat= &w_sinlat, *coslat= &w_coslat, *sig= &w_sig, *mu= &w_mu, *chi= &w_chi;
  static double ang_r0 = 0.0, sinang_r = 0.0, tanang_r = 0.0;
  static double latc0 = 0.0, coslatc = 1.0, sinlatc = 0.0;
  double cosr, sinr, sinlon, sinsig;

  if (p_rho) rho=p_rho;
  if (p_lat) lat=p_lat;
  if (p_lon) lon=p_lon;
  if (p_sinlat) sinlat=p_sinlat;
  if (p_coslat) coslat=p_coslat;
  if (p_sig) sig=p_sig;
  if (p_mu) mu=p_mu;
  if (p_chi) chi=p_chi;

  if (ang_r != ang_r0) {
      sinang_r = sin(ang_r);
      tanang_r = tan(ang_r);
      ang_r0 = ang_r;
    }

  if (latc != latc0) {
      sinlatc = sin(latc);
      coslatc = cos(latc);
      latc0 = latc;
    }

  *chi = atan2 (x, y) + pa;
  while (*chi > 2 * M_PI) *chi -= 2 * M_PI;
  while (*chi < 0.0) *chi += 2 * M_PI;
  /*  Camera curvature correction, no small angle approximations  */
  *sig = atan (hypot (x, y) * tanang_r);
  sinsig = sin (*sig);
  *rho = asin (sinsig / sinang_r) - *sig;
  if (*sig > ang_r) return (-1);
  *mu = cos (*rho + *sig);
  sinr = sin (*rho);
  cosr = cos (*rho);
  *sinlat = sinlatc * cosr + coslatc * sinr * cos (*chi);
  *coslat = sqrt (1.0 - *sinlat * *sinlat);
  *lat = asin (*sinlat);
  sinlon = sinr * sin (*chi) / *coslat;
  *lon = asin (sinlon);
  if (cosr < (*sinlat * sinlatc)) *lon = M_PI - *lon;
  *lon += lonc;
  while (*lon < 0.0) *lon += 2 * M_PI;
  while (*lon >= 2 * M_PI) *lon -= 2 * M_PI;
  return (0);
  }

int rebinArraySF(float *outData, float *inData, int outNx, int inNx, float rsun, float x0, float y0, float fillvalue)
  { //  Rebin only using data inside the Sun, else use fillvalue
  float rsun2 = 0.999*rsun*rsun; // require whole pixels for 4K images
  int nmiss = 0;
  float *inp, *outp;
  int inNy = inNx;
  int outNy = outNx;
  int inI, inJ, outI, outJ, i;
  int binsize = inNx/outNx;
  double inRow[inNx], outRow[outNx];
  int inRowN[inNx], outRowN[outNx];

  for (outJ=0; outJ < outNy; outJ++)
    {
    int inRow0 = outJ * binsize;
    for (outI=0; outI<outNx; outI++)
      {
      outRowN[outI] = 0;
      outRow[outI] = 0.0;
      }
    for (inJ = inRow0; inJ < inRow0+binsize; inJ++)
      {
      inp = inData + inJ*inNx;
      for (outI=0; outI<outNx; outI++)
        for (i=0; i<binsize; i++, inp++)
          {
          inI = outI*binsize + i;
          if (*inp != DRMS_MISSING_FLOAT && (((inI-x0)*(inI-x0) + (inJ-y0)*(inJ-y0)) < rsun2))
            {
            outRow[outI] += *inp;
            outRowN[outI]++;
            }
          }
       }
    for (outI=0; outI < outNx; outI++)
      if (outRowN[outI])
        outData[outJ*outNx + outI] = outRow[outI]/outRowN[outI];
      else
        {
        outData[outJ*outNx + outI] = fillvalue;
        nmiss++;
        }
    }  // end outJ
  }
