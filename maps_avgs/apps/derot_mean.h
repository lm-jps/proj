/*    derot_mean.h  */

/* **************************************************************
 *  
 * Derot_mean header file contains listings of 
 *    keywords and defined parameters 
 *    structures
 *    error codes
 *  
 *  
 *  **************************************************************
 */
#ifndef DEROT_MEAN_INCL
#define DEROT_MEAN_INCL
/* ****** Some Structures Used in derot_mean ****** */
/* Working Array and Geometry Value Structure */
typedef struct 
{
  DRMS_Array_t *DatArray;
  DRMS_Array_t *WgtArray;
  float *dat;    //
  float *wgt;
  float *lon; 
  float *lat;
  short  *osun;
  TIME   timctr;
  double LN;             // need to add the rest of the geometry stuff
  double LT;             // need to add the rest of the geometry stuff
  float Dmax;
  float Dmin;
  float Wmax;
  float Wmin;
  DRMS_Record_t *rec;
  double  irecdt;
  int   irecno;
  double  crln;    /* Carrington Lon for Observer - from ephem */
  int     crot;    /* Carrington Rotation */
} OutImgs;

/* Parameter List Structure */
typedef struct 
{
  char    oser[DRMS_MAXSERIESNAMELEN]; /* Output series name */
  char    iser[DRMS_MAXSERIESNAMELEN]; /* Input series name */
  TIME    Tstart;     /* Time of first output image */
  TIME    Tend;       /* Defines the time of last image <= Tend */
  int     Tstep;      /* Cadence of output images - in seconds */
  int     cols;       /* Output image number of columns */
  int     rows;       /* Output image number of rows */
  int     size;       /* Output image rows * cols columns */
  double  A0;         /* Solid body roation rate in (micro-rad/sec) */
  double  A2;         /* Rotation coef for Sin^2(B) term (micro-rad/sec) */
  double  A4;         /* Rotation coef for Sin^4(B) term (micro-rad/sec) */
  double  Meri_V;     /* Meridional Flow term in ??? units */
  char    WtFunc[10]; /* String specifying the type of weighting function */
  double  WtParm;     /* Parameter to set weighting function */
  int     WtLen;      /* Length of weighting function in seconds */
  double  *Wt;        /* Pointer to Weighting Function array */
  int     projcode;   /* Code for projection used for Rick's cartography */
  char    projname[100];  /* Type of projection for output mapping */
  double  rsun;       /* semi-diameter of Sun in pixels (required for AERIAL)*/
  double  dist;       /* Distance to Sun in A.U. (required for AERIAL) */
  double  ang_rad;    /* semi-diameter of Sun in radians (used for AERIAL)*/
  double  pangle;     /* Sets regular p-angle (in radians)=0 */
  double  CROTA2;     /* Sets to minus "P-angle" for output, usually=0 */
  double  CRLN;       /* Sets the L0 angle for output, usually=NaN */
  double  CRLT;       /* Sets the B0 angle for output, usually=0 */
  char    CTYPE1[100];     /* Type for Coordinate axis 1 */
  char    CTYPE2[100];     /* Type for Coordinate axis 2 */
  double  CRPIX1;     /* Fiducial Pixel for axis 1 */
  double  CRPIX2;     /* Fiducial Pixel for axis 2 */
  double  CRVAL1;     /* Cordinate for Fiducial pixel 1 */
  double  CRVAL2;     /* Cordinate for Fiducial pixel 2 */
  double  CDELT1;     /* Pixel size, axis 1 - defined at Fiducial Pixel */
  double  CDELT2;     /* Pixel size, axis 2 - defined at Fiducial Pixel */
  short   verbose;    /* Flag set if verbose */
  int     Nmaps;      /* Number of output maps (size of OutImg structure) */
  //  This is a poor man's error handling system, so errors can be pass up
  int     status;          /* to pass back error status codes from functions */
  char    Msg[2000];       /* to pass back error messages from functions */
  int     signal;          /* to pass back signals from functions */
  char    SigMsg[2000];    /* to pass back signal messages from functions */
} ParmList;

/* ******  Status Values & Error Codes ****** */
#define kMyMod_Success   (0)        // Operation was a success - no err code
#define kMyMod_InitErr   (101)      // Could not initialize array or variable 
#define kMyMod_ValErr    (102)      // A variable has invalid value        
#define kMyMod_MallocErr (103)      // Could not allocate memory              
#define kMyMod_Missing   (104)      // A keyword is missing        
#define kMyMod_WrongType (105)      // A variable is not the right type        
#define kMyMod_EndRec    (110)      // Signals that the last record
#define kMyMod_DONE      (999)      // Signals that the module is done


#endif
