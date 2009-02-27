/* procFDSData.h */

#ifndef _LIBASTRO_H
#define _LIBASTRO_H

#include "drms.h"
#include "cmdparams.h"

#define kALLDATAMISSING         0

typedef enum 
{
   kLIBASTRO_InterBilinear = 0,
   kLIBASTRO_InterCubic = 1,
} LIBASTRO_Interpolation_t;

typedef enum
{
   kLIBASTRO_Success = 0,
   kLIBASTRO_BadDimensionality,
   kLIBASTRO_BadData,
   kLIBASTRO_CouldntCreateData,
   kLIBASTRO_DimensionMismatch,
   kLIBASTRO_CantDoOddNumLats,
   kLIBASTRO_UnsupportedInterp,
   kLIBASTRO_UnsupportedMCOR,
   kLIBASTRO_UnsupportedVCOR,
   kLIBASTRO_InconsistentConstraints,
   kLIBASTRO_InvalidArgs,
   kLIBASTRO_Interpolation,
   kLIBASTRO_InsufficientData
} LIBASTRO_Error_t;

typedef enum
{
   kLIBASTRO_MCORLevel0,
   kLIBASTRO_MCORLevel1
} LIBASTRO_MCOR_t;

typedef enum
{
   kLIBASTRO_VCORSignOnly,
   kLIBASTRO_VCORLevel1,
   kLIBASTRO_VCORLevel2
} LIBASTRO_VCOR_t;


typedef struct LIBASTRO_Dist_Struct 
{
  int disttype;
  double scale; /* Image scale relative to FD. 5 for vw. */
  double xc, yc; /* Nominal image center in pixels. */
  double cdist; /* Cubic distortion constant for FD images */
  double feff, alpha, cosbeta, sinbeta; /* Tilt constants for FD. */
} LIBASTRO_Dist_t;

typedef struct LIBASTRO_RotRate_struct
{
  float lat;
  float r;
} LIBASTRO_RotRate_t;


int ConvAndInterpFDS(DRMS_Env_t *drmsEnv, char *seriesName, char *dateRange);
int SetDistort(CmdParams_t *params, 
	       char *distortP,
	       char *cubicP,
	       char *tiltaP,
	       char *tiltbP,
	       char *tiltfeffP,
	       LIBASTRO_Dist_t *dOut);
int Obs2helio(
	       float	*V,		/* input velocity array				*/
	       /* assumed declared V[ypixels][xpixels]		*/
	       float	*U,		/* output projected array			*/
	       int	xpixels,	/* x width of input array			*/
	       int	ypixels,	/* y width of input array			*/
	       double   x0, 	   	/* x pixel address of disk center		*/
	       double   y0, 	   	/* y pixel address of disk center		*/ 	  
	       double   BZero,          /* heliographic latitude of disk center         */
               double   P,		/* angle between CCD y-axis and solar vertical	*/       
	       /* positive to the east (CCW)			*/     
	       double   S,		/* subtended radius of the sun	radians    	*/
	       double   rsun,		/* pixels */
	       double   Rmax,		/* maximum disk radius to use (e.g. 0.95)	*/
	       
	       int	interpolation,	/* option */
	       int	cols,		/* width of output array			*/
	       int	rows,		/* height of output array			*/
	       double   Lmin,		/* start of longitude range 			*/
	       double   Ldelta,		/* increment of longitude per col	        */
	       double 	Ladjust,        
	       double   sinBdelta,	/* increment of sin latitude per row          */
	       
	       double	smajor,
	       double	sminor,
	       double	sangle,
	       double   xscale,
	       double   yscale,
	       const char *orientation,
	       
	       int	mag_correction,		/* option */
	       int	velocity_correction,	/* option */
	       double   obs_vr,
	       double   obs_vw,
	       double   obs_vn,
	       double   vsign,
	       int	NaN_beyond_rmax,
	       int      carrStretch,
	       const LIBASTRO_Dist_t *distpars,
	       float   diffrotA,
	       float   diffrotB,
	       float   diffrotC,
	       LIBASTRO_RotRate_t *rRates, /* rotation rates passed by reference*/
	       int        size);
int apodize(float *data,		/* input/output data array */
	    double b0,		/* heliographic latitude of disk center */
	    int cols, int rows,	/* width and height of input array */
	    double Lmin,		/* start of longitude range */
	    double Ldelta,		/* increment in longitude pre col */
	    double sinBdelta,	/* increment in sin latitude per row */
	    int apodlevel,		/* type of apodization */
	    /* 0 to do no apodization */
	    /* 1 to apodize in true solar coordinates */
	    /* 2 to apodize in ideal solar coordinates */
	    double apinner,		/* fractional radius to start apodization at */
	    double apwidth,		/* width of apodization */
	    int apel,		/* do elliptical apodization as described by 
				 * apx and apy */
	    double apx,		/* divide the x position by this before applying 
				   apodization */
	    double apy);	/* divide the y position by this before applying 
				   apodization */

float Ccint2(float *f, int nx, int ny, double x, double y);
double Ccint2d(double *f, int nx, int ny, double x, double y);
float Linint2(float *f, int nx, int ny, double x, double y);
double Linint2d(double *f, int nx, int ny, double x, double y);
int Regrid(DRMS_Array_t **dataArr, int *new_length, LIBASTRO_Interpolation_t scheme);
float Imaginterp(DRMS_Segment_t *img, double x, double y);

/* iorbit */
enum IORBIT_Alg_enum
{
   IORBIT_Alg_Linear = 0,
   IORBIT_Alg_Quadratic
};

typedef enum IORBIT_Alg_enum IORBIT_Alg_t;

/* Structure to return information caller of iorbit_getinfo() */
struct IORBIT_Info_struct
{
  double obstime;
  double hciX;
  double hciY;
  double hciZ;
  double hciVX;
  double hciVY;
  double hciVZ;
  double dsun_obs;
  double obs_vr;
  double obs_vw;
  double obs_vn;
  double crln_obs;
  double crlt_obs;
  double car_rot;
};

typedef struct IORBIT_Info_struct IORBIT_Info_t;

LIBASTRO_Error_t iorbit_test(DRMS_Env_t *env, const char *orbseries);
LIBASTRO_Error_t iorbit_getinfo(DRMS_Env_t *env, 
                                const char *srcseries, 
                                IORBIT_Alg_t alg,
                                const double *tgttimes, 
                                int nitems, 
                                const char *optfilter,
                                int flush,
                                IORBIT_Info_t **info);
void iorbit_cleanup();


#endif // _DRMS_LIBASTRO_H

