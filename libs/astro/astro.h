/* procFDSData.h */

#ifndef _LIBASTRO_H
#define _LIBASTRO_H

#include "drms.h"
#include "cmdparams.h"

#define kALLDATAMISSING         0
#define kMaxRecQuery            128

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


int ConvAndInterpFDS(DRMS_Env_t *drmsEnv, char *seriesName, char *dateRange);

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

enum IORBIT_CacheAction_enum
{
   kIORBIT_CacheAction_None = 0,
   kIORBIT_CacheAction_Flush,
   kIORBIT_CacheAction_DontCache
};

typedef enum IORBIT_CacheAction_enum IORBIT_CacheAction_t;

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
  double gciX;
  double gciY;
  double gciZ;
  double gciVX;
  double gciVY;
  double gciVZ;
  double dsun_obs;
  double rsun_obs;
  double obs_vr;
  double obs_vw;
  double obs_vn;
  double crln_obs;
  double crlt_obs;
  int car_rot;
  char orb_rec[kMaxRecQuery];
};

typedef struct IORBIT_Info_struct IORBIT_Info_t;

LIBASTRO_Error_t iorbit_test(DRMS_Env_t *env, const char *orbseries);
void iorbit_carrcoords(TIME t, double obsdist, double b, double hci_long, int *crot, double *L, double *B);
LIBASTRO_Error_t iorbit_getinfo(DRMS_Env_t *env, 
                                const char *srcseries, 
                                const char *optfilter, 
                                IORBIT_Alg_t alg,
                                const double *tgttimes, 
                                int nitems, 
                                IORBIT_CacheAction_t ctype,
                                IORBIT_Info_t **info);
void iorbit_cleanup();

void HeliographicLocation(TIME t, int *crot, double *L, double *B);
TIME HeliographicTime(int crot, double L);


#endif // _DRMS_LIBASTRO_H

