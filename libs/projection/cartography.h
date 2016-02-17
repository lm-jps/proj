#ifndef __CARTOGRAPHY_H__
#define __CARTOGRAPHY_H__

#ifdef __cplusplus
extern "C" {
#endif

	    /*  The following code is adapted from SOI functions/keywords.c  */
#define PLATFORM_UNKNOWN	(0)
#define PLATFORM_UNRECOGNIZED	(1)
#define PLATFORM_SOHO		(2)
#define PLATFORM_GONGPLUS	(3)
#define PLATFORM_MWO60		(4)
#define PLATFORM_BBSO		(5)
#define PLATFORM_TRACE		(6)
#define PLATFORM_SPOLE_JSO	(10)
#define PLATFORM_GONG		(30)
#define PLATFORM_OBSPM		(40)

#define INSTRUMENT_UNKNOWN	(0)
#define INSTRUMENT_UNRECOGNIZED	(1)
#define INSTRUMENT_SOHO_MDI	(10)
#define INSTRUMENT_SOHO_EIT	(11)
#define INSTRUMENT_GONG_TD	(20)
#define INSTRUMENT_GONG_CT	(21)
#define INSTRUMENT_GONG_TC	(22)
#define INSTRUMENT_GONG_BB	(23)
#define INSTRUMENT_GONG_ML	(24)
#define INSTRUMENT_GONG_LE	(25)
#define INSTRUMENT_GONG_UD	(26)
#define INSTRUMENT_GONG_MERGE	(29)
#define INSTRUMENT_MWO60_MOF	(30)
#define INSTRUMENT_BBSO_SINGER	(40)
#define INSTRUMENT_TRACE	(50)
#define INSTRUMENT_MOTH		(60)
#define INSTRUMENT_OBSPM_SPHG	(70)

#define NO_DATA_DICT	(0x0001)
#define NO_SEMIDIAMETER	(0x0002)
#define NO_XSCALE	(0x0004)
#define NO_YSCALE	(0x0008)
#define NO_XCENTERLOC	(0x0010)
#define NO_YCENTERLOC	(0x0020)
#define NO_HELIO_LATC	(0x0040)
#define NO_HELIO_LONC	(0x0080)
#define NO_HELIO_PA	(0x0100)
#define NO_XUNITS	(0x0200)
#define NO_YUNITS	(0x0400)
#define NO_OBSERVER_LAT	(0x0002)
#define NO_OBSERVER_LON	(0x0004)

#define KEYSCOPE_VARIABLE	(0x80000000)
#define LOCALHS_IMGINFO_VERSION	("1.0")

/*
 *  process image info (attitude, plate scale, distortion)
 *  stub function based on SOI solar_image_params()
 */
typedef struct paramdef {
  double scale;
  double offset;
  double defval;
  unsigned int statusbit;
  char name[32];
} ParamDef;


int img2sphere(double x, double y, double ang_r, double latc, double lonc, double pa, double *rho, double *lat, double *lon, double *sinlat, double *coslat, double *sig, double *mu, double *chi);
int plane2sphere (double x, double y, double latc, double lonc, double *lat, double *lon, int projection);
int sphere2img (double lat, double lon, double latc, double lonc, double *x, double *y, double xcenter, double ycenter, double rsun, double peff, double ecc, double chi,int xinvrt, int yinvrt);
int sphere2plane (double lat, double lon, double latc, double lonc, double *x, double *y, int projection);

void mtrack_MDI_correct_pa(double *pa);
void mtrack_MDI_correct_imgctr(double *xc, double *yc, double rsun);
void mtrack_MDI_image_stretch(double *x, double *y, int n, int direct);
void mtrack_MDI_image_tip(double *x, double *y, int n, int direct);

int solar_image_info(DRMS_Record_t *img, double *xscl, double *yscl, double *ctrx, double *ctry, double *apsd, const char *rsun_key, const char *apsd_key, double *pang, double *ellipse_e, double *ellipse_pa, int *x_invrt, int *y_invrt, int *need_ephem, int AIPS_convention);

int check_and_set_key_short(DRMS_Record_t *new, const char *key, short val);
int check_and_set_key_int(DRMS_Record_t *new, const char *key, int val);
int check_and_set_key_longlong (DRMS_Record_t *new, const char *key, long long val);
int check_and_set_key_float(DRMS_Record_t *new, const char *key, float val);
int check_and_set_key_double(DRMS_Record_t *new, const char *key, double val);
int check_and_set_key_str(DRMS_Record_t *new, const char *key, char *val);
int check_and_set_key_time(DRMS_Record_t *new, const char *key, TIME tval);
int check_and_copy_key(DRMS_Record_t *new, DRMS_Record_t *old, const char *key);
int drms_appendstr_tokey(DRMS_Record_t *rec, const char *key, const char *str, int addline);
void append_args_tokey(DRMS_Record_t *rec, const char *key);
int construct_stringlist(const char *request, char token, char ***stringlist);
int copy_prime_keys(DRMS_Record_t *new, DRMS_Record_t *old);
void string_insert_escape_char(char **str, const char esc);
int append_keyval_to_primekeyval(char **pkey, DRMS_Record_t *rec, const char *key);
char *create_primekey_from_keylist (DRMS_Record_t *rec, char **keylist, int keyct);
int propagate_keys(DRMS_Record_t *to, DRMS_Record_t *from, char **keylist, int keyct);
char *iau_units_parse_unit(char *unit, double *scale);
int drms_iau_units_scale(char *unit, double *scale);
int drms_wcs_timestep(DRMS_Record_t *rec, int axis, double *tstep);

#ifdef __cplusplus
}
#endif

#endif /* __CARTOGRAPHY_H__ */
