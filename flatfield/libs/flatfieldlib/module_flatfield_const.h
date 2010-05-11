

const int cam_id_front=3, cam_id_side=2; // !! check if HCAMID of CAMERA
const float fwhm=5.0;                   // fwhm for highpass filter
const int cthreshold=500; // minimum number of pairs of frames (minimum is 3)
const float threshold_lower=-10000.0;  //threshold for lower limit; //!!correction turned off 
const float threshold_upper=10000.0; //threshold for upper limit; // !! correction turned off

const int debug=1; //write out debug information

struct code_param cpa;
cpa.convergence=(double)1e-3;                //convergence threshold
cpa.maxiter=100;                       //maximum iterations
cpa.omega=1.4;                         //overrelaxation parameter

cpa.croprad=0.90;   //crop radius 
//cpa.rotcoef0= 2.913-0.1991;            //differential rotation coefficients (Komm Howard Harvey)
//cpa.rotcoef1=-0.405;
//cpa.rotcoef2=-0.422;
cpa.rotcoef0= 2.0*M_PI*1e-3*(452.- 31.7);         //differential rotation coefficients
cpa.rotcoef1= -2.0*M_PI*1e-3*49.0;
cpa.rotcoef2= -2.0*M_PI*1e-3*84.0;

const double normconst=0.01/600.;                   //regularization parameter
//other possible differential rotation models
 
 //omeg=2.0*M_PI*1.0e-9*(452.0- 49.0*pow(slat, 2) -84.0*pow(slat, 4) - 31.7)*time; // differential rotation rate [in rad/s]
 //omeg=1.0e-6*(2.838-0.301*pow(slat,2)-0.526*pow(slat, 4)-0.1991)*time; //Snodgrass 1984a
 //omeg=1e-6*(2.913-0.405*pow(slat,2)-0.422*pow(slat,4)-0.1991); //Komm, Howard, Harvey


//Limb darkening constants

const int b_order=5;
double b_coef[6]={0.34,1.37,-2.04,2.70,-1.94,0.559};   // Neckel and Labs (MDI wavelength), not good for extreme limb;
//b_coef[0]=0.34; 
//b_coef[1]=1.37;
//b_coef[2]=-2.04;
//b_coef[3]=2.70;
//b_coef[4]=-1.94;
//b_coef[5]=0.559;

//define fids

int kiconst;

const int minfid=10000;
const int maxfid=10199;
const int nfid=maxfid-minfid+1;      // number of FIDs 

int fid_list[nfid];
for (kiconst=minfid; kiconst<=maxfid; ++kiconst) fid_list[kiconst-minfid]=kiconst; //define FIDs

//cosmic ray detection
//set sigma for different FIDs, here constant sigma !! (cosmic rays)

double **sigma_fid;
sigma_fid=(double **)(malloc(nfid*sizeof(double*)));

// limits for lev1 keywords 
const float rsun_min=1850.0;
const float rsun_max=1940.0;

      const double dsun_obs_min=149597870691.0*0.9;
      const double dsun_obs_max=149597870691.0*1.1;

      const float p0_min=0;
      const float p0_max=360.0;

      const float b0_min=-7.4;
      const float b0_max=7.4;

      const float X0_min=2047.5-118.0;
      const float X0_max=2047.0+118.0;

      const float Y0_min=2047.5-118.0;
      const float Y0_max=2047.5+118.0;

      const float limit_centerdiff=0.2;
      const float limit_rsundiff=0.5;
/////


  //************************************************************************
  //constants for cosmic ray detection
  //************************************************************************

//define sigma (here:constant for each fid)
double constsigma[20]={281, 281, 281, 281, 281, 281.923,281.923,272.457,272.457, 257.217,257.217,270.393,270.393,381.843, 381.843,338.651, 338.651, 339, 339, 339};
float rad_cosmic_ray=0.9;

double *sigmacoef;
for (kiconst=0; kiconst<nfid; ++kiconst)
  {
 sigmacoef=(double *)(malloc(5*sizeof(double)));
 sigmacoef[0]=constsigma[(fid_list[kiconst]-10000)/10]*8.0; sigmacoef[1]=0.0; sigmacoef[2]=0.0; sigmacoef[3]=0.0; sigmacoef[4]=0.0;
 sigma_fid[kiconst]=sigmacoef;
  }



  const float factor=6.0;
  const long time_limit=120;

 
  //keyword names
   const char *keyfsn    = "FSN";                                 //1st prime key of input data
   const char *keytobs   = "T_OBS";                                //2nd prime key of input data
   const char *keycamera = "CAMERA";             //1st prime key of output data
   const char *keycam = "HCAMID"; //"HMI_SEQ_ID_EXP_PATH";
   const char *keyfocus  = "HCFTID"; 
   const char *keyfocusflat="HMI_SEQ_ID_FOCUS"; 
   const char *keyinstrument="INSTRUME";

   const char *keywl="HWLTID"; //"HMI_SEQ_ID_WLT";
   const char *keypl="HPLTID"; //"HMI_SEQ_ID_PST";

   const char *keytstart = "T_START";            //2nd prime key of output data               
   const char *keytstop = "T_STOP";
   const char *fidkey="FID";
   const char *isskey="HWLTNSET";
   const char *keyversion = "FLATFIELD_VERSION";

   const char *keynewpix = "ROTF_FLATFIELD";
   const char *keynpairs = "ROTF_N_PAIRS";
   const char *keycadence = "ROTF_CADENCE";

   const char *keylink1 = "T_OBS_OFFPOINT";
   const char *keylink2 = "T_OBS_DARK";
   const char *keylink3 = "T_OBS_BADPIX"; 

   const char *keycount = "COUNT";

   const char *linkoff = "OFFPOINT_FLAT";
   const char *linkdark = "DARK";
   const char *linkbad = "BAD_PIXEL";

   const char *segmentname="flatfield";
   const char *segmentname_badpix="bad_pixel_list";
   const char *segmentname_cosmic="cosmic_ray_hits";
   
   
   char *camera_str_front="HMI_FRONT2";
   char *camera_str_side="HMI_SIDE1";

//lev 1 keywords

const char *lev1_r_sun="RSUN_LF";
const char *lev1_imsc="CDELT1";
const char *lev1_dist="DSUN_OBS";
const char *lev1_p0="CROTA2";
const char *lev1_b0="CRLT_OBS";
const char *lev1_x0="X0_LF";
const char *lev1_y0="Y0_LF";

  //series names
char *filename_offpoint="su_richard.offpoint_flatfield"; // !! all su_richard for test
char *filename_dark="su_richard.dark";
char *filename_badpix="su_richard.bad_pixel_list";
char *filename_flatfield="su_richard.flatfield";
char *filename_flatfield_out="su_richard.flatfield";
char *filename_cosmic="su_richard.cosmic_rays";
  



