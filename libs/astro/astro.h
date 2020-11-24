/* procFDSData.h */

#ifndef _LIBASTRO_H
#define _LIBASTRO_H

#include "drms.h"
#include "cmdparams.h"

#define kALLDATAMISSING         0
#define kMaxRecQuery            128

/* These defines are suitable for sdo.fds_orbit_vectors, but not for other series. The
 * caller should pass in to the API functions a map (these defines --> series keywords)
 * of orbit-series keywords if they differ from these.
 */
#define IORBIT_STRINGIFY(x) #x

#define IORBIT_KEYWORD_GAE_X "GCIEC_X"
#define IORBIT_KEYWORD_GAE_Y "GCIEC_Y"
#define IORBIT_KEYWORD_GAE_Z "GCIEC_Z"
#define IORBIT_KEYWORD_GAE_VX "GCIEC_VX"
#define IORBIT_KEYWORD_GAE_VY "GCIEC_VY"
#define IORBIT_KEYWORD_GAE_VZ "GCIEC_VZ"
#define IORBIT_KEYWORD_HAE_X "HCIEC_X"
#define IORBIT_KEYWORD_HAE_Y "HCIEC_Y"
#define IORBIT_KEYWORD_HAE_Z "HCIEC_Z"
#define IORBIT_KEYWORD_HAE_VX "HCIEC_VX"
#define IORBIT_KEYWORD_HAE_VY "HCIEC_VY"
#define IORBIT_KEYWORD_HAE_VZ "HCIEC_VZ"
#define IORBIT_KEYWORD_OBS_DATE "OBS_DATE"

/* Some other series have these keywords. The union of this set and the one immediately above
 * form the set of all known keywords from any orbit series. */
#define IORBIT_KEYWORD_RSUN_OBS "RSUN_OBS"
#define IORBIT_KEYWORD_OBS_VR "OBS_VR"
#define IORBIT_KEYWORD_DSUN_OBS "DSUN_OBS"

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
    kLIBASTRO_InsufficientData,
    kLIBASTRO_Internal,
    kLIBASTRO_OutOfMemory,
    kLIBASTRO_DbStatement
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
  double hae_x;
  double hae_y;
  double hae_z;
  double hae_vx;
  double hae_vy;
  double hae_vz;
  double gae_x;
  double gae_y;
  double gae_z;
  double gae_vx;
  double gae_vy;
  double gae_vz;
  double dsun_obs;
  double rsun_obs;
  double obs_vr;
  double obs_vw;
  double obs_vn;
  double crln_obs;
  double crlt_obs;
  double hgln_obs;
  int car_rot;
  char orb_rec[kMaxRecQuery];
};

typedef struct IORBIT_Info_struct IORBIT_Info_t;

/* SAA-HLZ */
#define IORBIT_SAAHLZINFO_TIME_KEY_LEN             64
#define IORBIT_SAAHLZINFO_KW_EVENT_TYPE            "eventType"
#define IORBIT_SAAHLZINFO_KW_LEN                   32
#define IORBIT_SAAHLZINFO_KW_EVENT_TYPE_VALUE_LEN  8

typedef HContainer_t IORBIT_SaaHlzInfo_t;


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
LIBASTRO_Error_t iorbit_getinfo_ext(DRMS_Env_t *env,
                                    const char *srcseries,
                                    const char *optfilter,
                                    IORBIT_Alg_t alg,
                                    const double *tgttimes,
                                    int nitems,
                                    IORBIT_CacheAction_t ctype,
                                    IORBIT_Info_t **info,
                                    HContainer_t *keymap);
LIBASTRO_Error_t iorbit_getSaaHlzInfo(DRMS_Env_t *env,
                                      const char *series,
                                      const double *tgttimes,
                                      int nitems,
                                      IORBIT_SaaHlzInfo_t **info);

void iorbit_cleanup();
void iorbit_cleanSaaHlzInfo(IORBIT_SaaHlzInfo_t **info);

void HeliographicLocation(TIME t, int *crot, double *L, double *B);
TIME HeliographicTime(int crot, double L);


#endif // _DRMS_LIBASTRO_H
