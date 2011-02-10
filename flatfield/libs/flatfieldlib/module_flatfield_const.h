

const int cam_id_front=3, cam_id_side=2; // 
const double fwhm[2]={1.0, 1.0};                   // fwhm for highpass filter  
const int cthreshold=300; // minimum number of pairs of frames (minimum is 3) 
const float threshold_lower=-10000.0;  //threshold for lower limit; //!!correction turned off 
const float threshold_upper=10000.0; //threshold for upper limit; // !! correction turned off

const int update_flag=0; //0: no update, 1: single pixel update, 2:full update // 

const int debug=0; //write out debug information  

struct code_param cpa;

cpa.croprad=0.99;   //crop radius !changed from 0.98 09.30.2010 
//cpa.rotcoef0= 2.913-0.1991;            //differential rotation coefficients (Komm Howard Harvey)
//cpa.rotcoef1=-0.405;
//cpa.rotcoef2=-0.422;
cpa.rotcoef0= 2.0*M_PI*1e-3*(452.- 31.7);         //differential rotation coefficients
cpa.rotcoef1= -2.0*M_PI*1e-3*49.0;
cpa.rotcoef2= -2.0*M_PI*1e-3*84.0;

const double omegaconst[2]={1.7, 1.2};                         //overrelaxation parameter // 
const double normconst[2]={0.03/600., 0.01/600.};           //regularization parameter
const int maxiterconst[2]={300,200}; 
const double convconst[2]={1e-4, 1e-4};
//other possible differential rotation models
 
 //omeg=2.0*M_PI*1.0e-9*(452.0- 49.0*pow(slat, 2) -84.0*pow(slat, 4) - 31.7)*time; // differential rotation rate [in rad/s]
 //omeg=1.0e-6*(2.838-0.301*pow(slat,2)-0.526*pow(slat, 4)-0.1991)*time; //Snodgrass 1984a
 //omeg=1e-6*(2.913-0.405*pow(slat,2)-0.422*pow(slat,4)-0.1991); //Komm, Howard, Harvey


//Limb darkening constants

const int b_order=5;
double b_coef[6]={0.34,1.37,-2.04,2.70,-1.94,0.559};   // Neckel and Labs (MDI wavelength), not good for extreme limb;

//static double LDCoef[] = {0.0, 0.4430, 0.1390, 0.041, 0.0125, 0.0019};

//b_coef[0]=0.34; 
//b_coef[1]=1.37;
//b_coef[2]=-2.04;
//b_coef[3]=2.70;
//b_coef[4]=-1.94;
//b_coef[5]=0.559;

//define fids

//int kiconst;
const int minfid=10000;
const int maxfid=10199;
const int nfid=(maxfid-minfid+1);      // number of FIDs 

//int fid_list[nfid];
//for (kiconst=minfid; kiconst<=maxfid; ++kiconst) fid_list[kiconst-minfid]=kiconst; //define FIDs

//cosmic ray detection



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

      const float vrad_min=-5000.0;
      const float vrad_max=5000.0;

      const float limit_centerdiff[2]={0.1, 0.05};
const float limit_centerdiff_cosmic=1.0;
      const float limit_rsundiff=0.1;
/////


  //************************************************************************
  //constants for cosmic ray detection
  //************************************************************************

//define sigma (here:constant for each fid)
//double constsigma[20]={281., 281., 281., 281., 281., 281.923,281.923,272.457,272.457, 257.217,257.217,270.393,270.393,381.843, 381.843,338.651, 338.651, 339., 339., 339.};
float rad_cosmic_ray=0.98;
long time_limit=200; //time distance in sec for filtergrams with same FID to be used in cosmic ray detection
float limit_offpoint=20.0; //limit for offpoint of individual filtergram from average offpoint
float filtergram_cadence=1.8; //minimum filtergram cadence in sec

//constants for cosmic ray detection as a function of effective wavelength

//set derive from standard deviation of time series
//float coef0[2][3]={{559.763 , 466.619, -840.870}, {516.289, 406.209, -730.977}};
//float coef1[2][3]={{1354.98, -70.6893, -1014.83}, {1281.43,230.564,-1225.94}};
//float coef2[2][2]={{0.0156957, 0.0167262},{0.0124592, 0.0212838}};
//float coef3[2][2]={{0.0200183, 0.0277104},{0.0211192,0.0272843}};
//float coef4[2][3]={{1585.98+200.0, -21.0609, -1044.63},{1577.13, 163.140, -1170.68}};

//set derived from standard deviation of derivative of time series
//float coef0[2][3]={{1395.02, -491.606, -606.292}, {1012.51, -374.209, -442.236}};
//float coef1[2][3]={{2519.27, -298.148, -2029.38}, {1370.68, -249.615, -982.241}};
//float coef2[2][3]={{0.0216513, 0.00287712, 0.0198303}, {0.0282894, -0.00629716, 0.0302188}};
//float coef3[2][2]={{0.17, 0.27}, {0.17, 0.27}};
//float coef4[2][3]={{2332.02, -0.489186,  -1549.49}, {1028.30,   75.0354, -662.894}};

float coef0[2][3]={{1151.30, -130.552, -704.380},{916.720,     -291.264  ,   -436.917}};
float coef1[2][3]={{2517.30, -339.439, -1943.22},{1369.02 ,    -424.949,     -808.191}};
float coef2[2][3]={{0.0262323, -0.00716024, 0.0240938},{0.0302306, -0.0134902, 0.0346718}};
float coef3[2][2]={{0.17, 0.27}, {0.17, 0.27}};
float coef4[2][3]={{2255.87, 117.511, -1541.76},{1034.88,  95.1960, -659.753}};


const float lambda0=6173.3433;
const float lambda_sep=68.8e-3;
const float v_c=2.99792e8;
const float radsun_mm=695.5;

float cof_combx[12]={0.145544,0.00812206, 0.233276, 0.121743, 0.190385,0.105851, 0.135802,0.0402982, 0.0571681,-0.0916504,0.0853409,-0.0318835};
float cof_comby[12]={0.283009, 0.0545081, 0.235930, 0.0957727, 0.245011, 0.0815215, 0.124488, 0.0380576, 0.0556949, -0.0611975, -0.0492513, -0.103545};


//float cof[2][6]={{1.22201 ,    0.362697 ,    0.359817 ,   -0.142182 ,   -0.433347 ,   -0.368998}, {0.300682 ,    0.370525,     0.316930 ,    0.164742,   -0.0441823 ,   -0.108697}}; 
float cof[2][6]={{0.187514, 0.542523, 0.284918, 0.143162,-0.136262,-0.0218554}, {0.311833, 0.400153, 0.205032,0.0928222,-0.00851838, -0.00132237}}; //for filter fwhm = 1.0

const float factor[2][2]={{6.2, 9.3}, {6.2, 9.3}}; 
const int limit_cosmic=100000;
  

 
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
   const char *flatnkey="FLAT_REC";
   const char *keyversion = "FLATFIELD_VERSION";
   const char *fsnskey="FSN_START";
   const char *keyday="DAY";
   const char *recnumkey="RECNUMOFF";
   const char *querykey="FLATQ";

   const char *keynewpix = "ROTF_FLATFIELD";
   const char *keynpairs = "ROTF_N_PAIRS";
   const char *keycadence = "ROTF_CADENCE";

   const char *keylink1 = "T_OBS_OFFPOINT";
   const char *keylink2 = "T_OBS_DARK";
   const char *keylink3 = "T_OBS_BADPIX"; 

   const char *keycount = "COUNT";
   const char *keyexmax = "EXMAX";
   const char *keylimit = "DETLIM";
   const char *linkoff = "OFFPOINT_FLAT";
   const char *linkdark = "DARK";
   const char *linkbad = "BAD_PIXEL";

   const char *segmentname="flatfield";
   const char *segmentname_offpoint="offpoint_flatfield";
   const char *segmentname_badpix="bad_pixel_list";
   const char *segmentname_dark="dark";
   const char *segmentname_cosmic="cosmic_ray_hits";
   const char *segmentname_val="level";
   const char *segmentname_sig="significance";
   
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
const char *lev1_vr="OBS_VR";


// !! checked in version

  //series names

//char *filename_flatfield_out="su_richard.flatfield";  //output of module_flatfield_combine
char *filename_flatfield_out="hmi.flatfield"; //for checked in version

char *filename_flatfield_series="hmi.flatfield";

char *filename_badpix="hmi.bad_pixel_list";
char *filename_offpoint="hmi.offpoint_flatfield";
char *filename_dark="hmi.dark";

//char *filename_flatfield_rel="su_richard.flatfield_rel";  //output of module_flatfield_combine
char *filename_flatfield_rel="su_production.flatfield_rel"; //for checked in version

//char *filename_flatfield_fid="su_richard.flatfield_fid_b"; //output of module_flatfield // input to module_flatfield_combine
char *filename_flatfield_fid="su_production.flatfield_fid"; //for checked in version

//char *filename_cosmic="su_richard.cosmic_rays"; //output of module_flatfield
char *filename_cosmic="su_production.cosmic_rays"; //for checked in version  


//char *filename_cosmic2_out="su_richard.cosmic_rays_b";
char *filename_cosmic2_out="hmi.cosmic_rays";


