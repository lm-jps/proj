/******************************************************************************/
/*
 *  maicalc.c						~rick/hmi/rings/src/
 *
 *  Module description
 *    Integrates functions of values in a selected set of regions from sets
 *	of line-of-sight and/or vector magnetograms after mapping in order
 *	to form time-region indices.
 *    The 
 *
 *  Responsible: Rick Bogart				rick@sun.stanford.edu
 *
 *  Usage:
 *    maicalc [-Mflnrsv] los= length= rec_step= scale= extent= [arg= ...]
 *
 *  Arguments: (type	default	description)
 *	los	str	-	input LOS magnetogram data series (or
 *				data set)
 *	los_seg	str	NotSpec	segment name(s) for field in LOS series
 *	vec	str	NotSpec	input vector magnetogram data series
 *				(or data set)
 *	mai	str	NotSpec	output data series
 *	length	int	-	length of sampling interval, in units of data cadence
 *	rec_step int	-	sampling step size [rec]
 *	cr	int	current	Carrington rotation
 *	cl	float	current or 180 Carrington longitude
 *	scale	float	-	map scale in heliographic deg / pixel
 *	extent	float	-	map extent in degrees
 *	lat	float*	[0.0]	target heliographic latitude(s)
 *	lon	float*	[0.0]	target heliographic longitude(s)
 *	mask_in	float	0.9765625 inner edge for 1-r^2 apodization
 *	mask_ex	float	1.000	outer edge for 1-r^2 apodization
 *	apodize	float	0.96875	temporal apodization edge
 *	floor	float	50.0	"noise" floor for inclusion in integration
 *				[gauss]
 *	reject	string	NotSpec	If specified, name of a file with a
 *				list of input images to be rejected unconditionally
 *	qmask	int	0x80000000 Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	trec_key string	T_REC		Keyname of time type prime keyword for
 *				input data series
 *	tobs_key string	T_OBS		Keyname of time type keyword describing
 *				observation time (midpoint) of each input image
 *	tstp_key string	Cadence		Keyname of float type keyword describing
 *				observing cadence of input data
 *	qual_key string	Quality	Keyname of int type keyword describing
 *				record quality as bitmask
 *	clat_key string	CRLT_OBS
 *	clon_key string	CRLN_OBS
 *	crot_key string	CAR_ROT
 *	rsun_key string	R_SUN
 *	apsd_key string	RSUN_OBS
 *	dsun_key string	DSUN_OBS Keyname of double type keyword for
 *				distance from sun in m for of each image
 *	cvkey	 string	CalVer64	Key name of 64-bit mask type keyword
 *				describing calibration version used
 *	max_reach float	0.5	maximum distance from target sample to search
 *				when target is unaaceptable, in units of
 *				sample spacing
 *
 *  Flags
 *	-f	force filling or reporting of records for which there are
 *		no valid data
 *	-l	list MAI's on stdout
 *	-n	turn off tracking
 *	-r	turn off apodization
 *	-s	turn off despiking (elimination of outliers)
 *	-v	run verbose
 *	-M	correct for MDI distortion using new model
 *	-O	correct for MDI distortion using old model
 *	-S	suppress correction of SDO Carrington Longitude
 *
 *  Bugs:
 *    Variables los_series and vec_series are unused
 *    The separate options for "old" and "new" distortion corrections of MDI
 *	data should eventually be removed in favor of a single acceptable
 *	one. The "new" correction is currently included as a local function;
 *	with parameters based on analysis of Dopplergrams that is not really
 *	applicable to magnetograms. It is not clear whether the position
 *	angle correction needs to be applied when the "new" distortion
 *	correction is used; currently it is.
 *    Data from time samples outside the temporal apodization window, if used,
 *	are ignored.
 *    There is no option for differential tracking; regions are referred to
 *	either their Carrington coordinates at the time of analysis or are
 *	"derotated" to those at the midtime of the interval if the no-track
 *	option is selected.
 *    Not all arg parameter values are propagated to output records: in
 *	particular, length, rec_step, reject, qmask, and max_reach, though
 *	their cumulative effect can be inferred from the list of contributing
 *	magnetograms in the HISTORY
 *    Processing of vector-field data is not implemented
 *    Need to look at effect of eliminate_outliers
 *    Will need separate quality and rejection keys for LOS and Vector if
 *	processed at same time
 *    Only a few of the essential output keywords are checked; the series
 *	structure check really only checks existence
 *    Calculation of maximum distance from center of image may not be correct
 *	in case of no tracking
 *    Having a default for target time is a bad idea; should be required
 *
 *  Revision history is at the end of the file
 */
/******************************************************************************/

/******************** defines, includes, and global declarations **************/
#define MODULE_VERSION_NUMBER	("2.4")
#define KEYSTUFF_VERSION "keystuff_v11.c"

#include <jsoc_main.h>

#include KEYSTUFF_VERSION

char *module_name = "maicalc";
char *module_desc = "integration of mapped and tracked magnetogram data";
char *version_id = MODULE_VERSION_NUMBER;

#define OBS_CRLN_FIXED	(0x10000000)
#define OBS_CROT_FIXED	(0x00000010)

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"los", "",
      "input LOS magnetogram data series or dataset"}, 
  {ARG_STRING,	"los_seg", "Not Specified",
      "LOS magnetogram series segment (default: first rank 2)"}, 
  {ARG_STRING,	"vec", "Not Specified",
      "input vector magnetogram data series or dataset"}, 
  {ARG_DATASERIES, "mai", "Not Specified", "output data series for MAI\'s"},
  {ARG_STRING,	"tmid", "Not Specified", "midpoint of tracking interval"}, 
  {ARG_INT,	"cr", "Not Specified", "Carrington Rotation (default: current)"},
  {ARG_FLOAT,	"cl", "Not Specified",
      "Carrington longitude of central meridian (default: current or 180"},
  {ARG_INT,   "length", "",
      "length of sampling interval, in units of input cadence", "(1,)"},
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_INT,     "rec_step", "", "sampling step, record steps", "[1,)"},
  {ARG_FLOAT,	"max_reach", "0.5", "maximum reach to acceptable record",
      "[0.0,1.0]"},
  {ARG_FLOAT,	"scale", "",
      "integrating map scale for interpolation [deg/pixel]", "(0.0,)"},
  {ARG_FLOAT,	"extent", "", "extent of integrating maps [degrees]"},
  {ARG_FLOATS,	"lat", "",
      "heliographic latitude(s) of tracking center(s) [deg]"},
  {ARG_FLOATS,	"lon", "",
      "heliographic longitude(s) of tracking center(s) [deg]"},
  {ARG_FLOAT,   "mask_in", "0.9765625", "inner edge of 1 - r^2 apodization",
      "[0.0,)"},
  {ARG_FLOAT,   "mask_ex", "1.0", "outer edge of 1 - r^2 apodization", "[0.0,)"},
  {ARG_FLOAT,	"apodize", "0.96875", "temporal apodization edge",
	"[0.0,1.0]"},
  {ARG_FLOAT,   "floor", "50.0", "noise floor (gauss)", "[0.0,)"},
  {ARG_STRING,	"tobs_key", "T_OBS", "keyname for observation time"}, 
  {ARG_STRING,	"trec_key", "T_REC", "keyname of (slotted) prime key"}, 
  {ARG_STRING,	"tstp_key", "CADENCE",  "keyname for image observation time"}, 
  {ARG_STRING,	"qual_key", "Quality",  "keyname for 32-bit image quality field"}, 
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"crot_key", "CAR_ROT", "keyname for image Carrington rotation"}, 
  {ARG_STRING,	"rsun_key", "R_SUN", "keyname for image semi-diameter (pixel)"}, 
  {ARG_STRING,	"apsd_key", "RSUN_OBS",
      "keyname for apparent solar semi-diameter (arcsec)"}, 
  {ARG_STRING,	"dsun_key", "DSUN_OBS", "keyname for observer distance"}, 
  {ARG_STRING,	"cvkey", "CalVer64", "keyname for Calibration Version key"}, 
  {ARG_FLAG,	"f",	"", "force filling/reporting of records with no data"}, 
  {ARG_FLAG,	"l",	"", "list MAI values only on stdout"}, 
  {ARG_FLAG,	"n",	"", "turns off tracking; target cl only determines time"}, 
  {ARG_FLAG,	"r",	"", "turns off apodization"}, 
  {ARG_FLAG,	"s",	"", "turns off despiking (elimination of outliers)"}, 
  {ARG_FLAG,	"v",	"", "verbose mode (overrides -l flag)"},
  {ARG_FLAG,	"M",	"", "correct for MDI distortion using new model"}, 
  {ARG_FLAG,	"O",	"", "correct for MDI distortion using old model"}, 
  {ARG_FLAG,	"S",	"", "suppress correction of SDO Carrington longitude"}, 
  {}
};

#include "selstuff.c"
#include "soho_ephem.c"
#include "cartography.c"
#include "imginfo.c"
#include "mdistuff.c"

typedef struct {
  TIME trec;
  TIME tobs;
  double clon;
  double clat;
  double dist;
  double apsd;
  int crot;
} OBS_EPHEMERIS;

/***************************** internal functions *****************************/
/*
 *  Interpolate Carrington longitude of disc center, and if necessary
 *    Carrington Rotation number, from tobs to trec
 */
int interp_crlfix (double *tobs, double trec, int nt, double **crcl,
    double *clfix, int *crfix) {
  double clt, dt, dl, fac;
  int n = 0;
			      /*  if target out of bounds, can't interpolate  */
  if (trec < tobs[0] || trec > tobs[nt-1]) return 1;
  while (trec >= tobs[n] && n < nt) n++;
  if (n == nt) return 1;
					  /*  likewise if obstime is missing  */
  if (time_is_invalid ((TIME)(tobs[n])) ||
      time_is_invalid ((TIME)(tobs[n-1]))) return 1;
  dt = tobs[n] - tobs[n-1];
  if (isnan (dt)) return 1;
  fac = (trec - tobs[n-1]) / dt;
  dl = crcl[1][n] - crcl[1][n-1];
  *crfix = crcl[0][n-1] + 0.001;
  if (dl > 180.0) {
    dl -= 360.0;
    *crfix += 1;
  }
  if (dl < -180.0) {
    dl += 360.0;
    *crfix -= 1;
  }

  clt = crcl[1][n-1] + fac * dl;
  if (clt < 0.0) {
    clt += 360.0;
    *crfix += 1;
  }
  *clfix = clt;
  return 0;
}
/*
 *  Interpolate standard, non-periodic variables from tobs to trec
 */
int interp_hmi_fixes (double *x, double x0, int nx, double **yt, double *y,
    int yct) {
  double dx, dy, fac;
  int i;
  int n = 0;
			      /*  if target out of bounds, can't interpolate  */
  if (x0 < x[0] || x0 > x[nx-1]) return 1;
  while (x0 >= x[n] && n < nx) n++;
  if (n == nx) return 1;
					  /*  likewise if obstime is missing  */
  if (time_is_invalid ((TIME)(x[n])) ||
      time_is_invalid ((TIME)(x[n-1]))) return 1;
  dx = x[n] - x[n-1];
  fac = (x0 - x[n-1]) / dx;
  for (i = 0; i < yct; i++) {
    dy = yt[i][n] - yt[i][n-1];
    y[i] = yt[i][n-1] + fac * dy;
  }
  return 0;
}
/*
 *  Adjust ephemeris values applicable to T_OBS to T_REC by interpolation
 */
void adjust_obs2rec (OBS_EPHEMERIS *rec, DRMS_RecordSet_t *ids, int ct,
    char *trec_key, char *tobs_key, char *crot_key, char *clon_key,
    char *clat_key, char *dsun_key, char *apsd_key, TIME tstrt, TIME tstop,
    double data_cadence, char *source_series) {
  DRMS_Record_t *irec;
  DRMS_Array_t *sdoephvec, *sdocrlvec;
  double **efixp, **crlfix;
  double *tfix, *efix;
  double loncfix;
  int l, m, n, tblct, valct;
  int nr, status;
  int fixct, efixct, crfix;
  char rec_query[DRMS_MAXQUERYLEN];
  char keylist[256], tbuf0[64], tbuf1[64];
				   /*  fill ephemeris vector from key values  */
  for (nr = 0; nr < ct; nr++) {
    irec = ids->records[nr];
    rec[nr].trec = drms_getkey_time (irec, trec_key, &status);
    rec[nr].tobs = drms_getkey_time (irec, tobs_key, &status);
    rec[nr].crot = drms_getkey_int (irec, crot_key, &status);
    rec[nr].clon = drms_getkey_double (irec, clon_key, &status);
    rec[nr].clat = drms_getkey_double (irec, clat_key, &status);
    rec[nr].dist = drms_getkey_double (irec, dsun_key, &status);
    rec[nr].apsd = drms_getkey_double (irec, apsd_key, &status);
  }
  sprint_time (tbuf0, tstrt - data_cadence, "TAI", 0);
  sprint_time (tbuf1, tstop + data_cadence, "TAI", 0);
  sprintf (rec_query, "%s[%s-%s]\n", source_series, tbuf0, tbuf1);
  sprintf (keylist, "%s,%s,%s,%s", tobs_key, clat_key,
        dsun_key, apsd_key);
  sdoephvec = drms_record_getvector (drms_env, rec_query, keylist,
	DRMS_TYPE_DOUBLE, 0, &status);
  tblct = sdoephvec->axis[1];
  valct = sdoephvec->axis[0];
  tfix = (double *)malloc (tblct * sizeof (double));
  for (n = 0; n < tblct; n++)
    tfix[n] = ((double *)sdoephvec->data)[n];
  efixp = (double **)malloc (valct * sizeof (double *));
  efix = (double *)malloc (valct * sizeof (double));
  for (l = 0; l < valct; l++) {
    efixp[l] = (double *)malloc (tblct * sizeof (double));
    for (m = 0; m < tblct; m++, n++)
      efixp[l][m] = ((double *)sdoephvec->data)[n];
  }
  drms_free_array (sdoephvec);
  fixct = tblct;
  efixct = valct;
  sprint_time (tbuf0, tstrt - data_cadence, "TAI", 0);
  sprint_time (tbuf1, tstop + data_cadence, "TAI", 0);
  sprintf (rec_query, "%s[%s-%s]\n", source_series, tbuf0, tbuf1);
  sprintf (keylist, "%s,%s,%s", tobs_key, crot_key, clon_key);
  sdocrlvec = drms_record_getvector (drms_env, rec_query, keylist,
      DRMS_TYPE_DOUBLE, 0, &status);
  valct = sdocrlvec->axis[0];
  crlfix = (double **)malloc (valct * sizeof (double *));
  n = tblct;
  for (l = 0; l < valct; l++) {
    crlfix[l] = (double *)malloc (tblct * sizeof (double));
    for (m = 0; m < tblct; m++, n++)
      crlfix[l][m] = ((double *)sdocrlvec->data)[n];
  }
  drms_free_array (sdocrlvec);
  for (nr = 0; nr < ct; nr++) {
    status = interp_crlfix (tfix, (double)(rec[nr].trec), fixct,
	crlfix, &loncfix, &crfix);
    if (status) continue;
    rec[nr].crot = crfix;
    rec[nr].clon = loncfix;
    status = interp_hmi_fixes (tfix, (double)(rec[nr].trec), fixct,
	efixp, efix, efixct);
    if (status) continue;
    rec[nr].clat = efix[0];
    rec[nr].dist = efix[1];
    rec[nr].apsd = efix[2];
  }
}
		/*  new local function for MDI distortion correction,
					eventually goes into mdistuff.c */
/*
 * mtrack_MDI_image_distort
 *
 * Correct MDI image for distortion as determined by comparison with HMI
 *   Dopplergrams.
 *
 */

void mtrack_MDI_image_distort (double *x, double *y, int n, int direct) {
  double x0, y0, dx, dy;
  double x2, y2, x3, y3, xy, xxy, yyx;
/*
  double coefx[10] = {5.80955670e-09, -2.89719337e-10, -1.22504236e-05,
      -5.14509007e-06,  8.20720819e-03, 4.72616540e-03, -1.89092055e+00,
      -1.35952890e-09,  5.75679518e-09, -4.40762222e-06};
  double coefy[10] = {2.29602983e-10, 5.10638195e-09, -2.79825322e-06,
      -6.08401833e-06,  3.02299760e-03,  3.32798869e-03, -5.60494963e-01,
       6.26347786e-09, -6.08012807e-10, -7.44956991e-06};
*/
  double coefx[10] = {1.52294042e-03, -7.59482821e-05, -2.06390160e-03,
      -1.35427259e-03, -1.22480162e-03, -3.62530274e-04, 6.78226392e-04,
      -3.56392399e-04, 1.50910931e-03, 4.64799760e-05};
  double coefy[10] = {6.01888856e-05, 1.33860741e-03, 3.88014840e-04,
      7.37651827e-04, -3.51475902e-04, -1.37785819e-03, -2.22602054e-04,
      1.64193314e-03, -1.59386922e-04, -8.51982926e-04};
  while (n--) {
    x0 = *x;
    y0 = *y;
    x2 = x0*x0;
    y2 = y0*y0;
    x3 = x2*x0;
    y3 = y2*y0;
    xy = x0*y0;
    xxy = x2*y0;
    yyx = y2*x0;
    dx = coefx[0]*x3 + coefx[1]*y3 + coefx[2]*x2 + coefx[3]*y2 + coefx[4]*x0 +
	coefx[5]*y0 + coefx[6] + coefx[7]*xxy + coefx[8]*yyx + coefx[9]*xy;
    dy = coefy[0]*x3 + coefy[1]*y3 + coefy[2]*x2 + coefy[3]*y2 + coefy[4]*x0 +
	coefy[5]*y0 + coefy[6] + coefy[7]*xxy + coefy[8]*yyx + coefy[9]*xy;
    *x += dx*direct;
    *y += dy*direct;
    x++;
    y++;
  }
}

typedef enum {LOC_UNKNOWN, LOC_GEOC, LOC_MWO, LOC_GONG_MR, LOC_GONG_LE,
    LOC_GONG_UD, LOC_GONG_TD, LOC_GONG_CT, LOC_GONG_TC, LOC_GONG_BB,
    LOC_GONG_ML, LOC_SOHO, LOC_SDO} platloc;

#define RSUNM		(6.96e8)

#define OUTLIER_RATIO	(6.0)
#define OUTLIER_BASE	(400.0)
		 /*  global declaration of missing to be initialized as NaN  */
float missing_val;
				/*  adapted from SOI libM/interpolate.c  */
/*
 *  Cubic convolution interpolation, based on GONG software, received Apr-88
 *  Interpolation by cubic convolution in 2 dimensions with 32-bit float data
 *
 *  Reference:
 *    R.G. Keys in IEEE Trans. Acoustics, Speech, and Signal Processing, 
 *    Vol. 29, 1981, pp. 1153-1160.
 *
 *  Usage:
 *    float ccint2 (float *f, int nx, int ny, double x, double y)
 *    double ccint2d (double *f, int nx, int ny, double x, double y)
 *
 *  Bugs:
 *    Currently interpolation within one pixel of edge is not
 *      implemented.  If x or y is in this range or off the picture, the
 *      function value returned is NaN.
 */
				/*  Calculate the interpolation kernel.  */
static void ccker (double *u, double s) {
  double s2, s3;

  s2= s * s;
  s3= s2 * s;
  u[0] = s2 - 0.5 * (s3 + s);
  u[1] = 1.5*s3 - 2.5*s2 + 1.0;
  u[2] = -1.5*s3 + 2.0*s2 + 0.5*s;
  u[3] = 0.5 * (s3 - s2);
}
					/*  Cubic convolution interpolation  */
float ccint2 (float *f, int nx, int ny, double x, double y) {
  double  ux[4], uy[4];
  double sum;
  int ix, iy, ix1, iy1, i, j;

  if (x < 1.0 || x >= (float)(nx-2) || y < 1.0 || y >= (float)(ny-2))
    return missing_val;

  ix = (int)x;
  ccker (ux,  x - (double)ix);
  iy = (int)y;
  ccker (uy,  y - (double)iy);

  ix1 = ix - 1;
  iy1 = iy - 1;
  sum = 0.;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      sum = sum + f[(iy1+i) * nx + ix1 + j] * uy[i] * ux[j];
  return (float)sum;
}
					  /*  nearest value "interpolation"  */
float nearest_val (float *f, int cols, int rows, double x, double y) {
  int col, row;
  if (x < -0.5 || x >= (cols - 0.5) || y < -0.5 || y >= (rows - 0.5))
    return missing_val;
  col = x + 0.5;
  row = y + 0.5;
  return f[col + row * cols];
}

#define INTERP_NEAREST_NEIGHBOR	(1)
				/*  adapted from SOI functions/sds_interp.c  */
float array_imaginterp (DRMS_Array_t *img, double x, double y, int schema) {
/*
 *  Interpolate to an arbitrary grid location {x, y} from a DRMS Array
 *    containing a projected solar image.  The aim of this function is
 *    is to provide an ideal interpolation weighted by foreshortening,
 *    limb darkening, and vector projection, but for now this is simply
 *    a stub function that extracts information from the attributes of
 *    the dataset and calls a simple two-dimensional cubic convolutional
 *    interpolation function.
 *
 *  x and y are in the range [-1,-1] at the "lower left" of the first pixel
 *    to [1,1] at the "upper right" of the last pixel in the image.
 *  (These are converted to the ccint2 conventions, with x and y in
 *    the range [0,0] at the "center" of the first pixel to
 *    [cols-1, rows-1] at the "center" of the last pixel.)
 *  Interpolation near the edges is not allowed.
 *
 *  Bugs:
 *    Interpolation within one pixel of edge is not implemented.  If x or y
 *	is in this range or off the image, the function returns zero.
 *    The function assumes a fixed scale in both directions, so that if one
 *	dimension is larger than another the scale is applied to the larger.
 *    Only floating point data are supported by the function, and there is
 *	not even any testing for validity.
 */
  double xs, ys;
  int cols, rows, mdim;

  cols = img->axis[0];
  rows = img->axis[1];
  mdim = (cols > rows) ? cols : rows;
  xs = 0.5 * (x + 1.0) * mdim - 0.5;
  ys = 0.5 * (y + 1.0) * mdim - 0.5;
  if (schema == INTERP_NEAREST_NEIGHBOR)
    return nearest_val (img->data, cols, rows, xs, ys);
  else return ccint2 (img->data, cols, rows, xs, ys);
}
/*
 *  Functions to support reading from a specified rejection list of the
 *    appropriate format
 */
int fgetline (FILE *in, char *line, int max) {
  if (fgets (line, max, in) == NULL) return 0;
  else return (strlen (line));
}

int read_reject_list (FILE *file, int **list) {
  int ds, sn, rec, last_rec;
  int allocd = 1024, ct = 0, gap = 0;
  char line[1024], t_str[64], estr[16];

  *list = (int *)malloc (allocd * sizeof (int));
  while (fgetline (file, line, 1024)) {
    if (strlen (line) == 1) continue;
    if (line[0] == '#') continue;
    if (sscanf (line, "%d %d %d %s", &ds, &sn, &rec, t_str) != 4) {
      if (sscanf (line, "%d %s", &rec, t_str) != 2) {
        sscanf (line, "%s", estr);
        if (strcmp (estr, "...")) continue;
        gap = 1;
        last_rec = rec;
        continue;
      }
    }
    if (gap) {
      while (rec > last_rec) {
	last_rec++;
	(*list)[ct++] = last_rec;
	if (ct >= allocd) {
	  allocd += 1024;
	  *list = (int *)realloc (*list, allocd * sizeof (int));
	}
      }
      gap = 0;
      continue;
    }
    (*list)[ct++] = rec;
    if (ct >= allocd) {
      allocd += 1024;
      *list = (int *)realloc (*list, allocd * sizeof (int));
    }
  }
  return ct;
}

double apodization (double x, double y, double inner, double outer) {
  double f, r = hypot (x, y), r2;

  if (inner >= outer) r = (r < outer) ? 0.0 : 1.0;
  else r = (r - inner) / (outer - inner);
  r2 = 1.0 - r * r;
  f = (r >= 1.0) ? 0.0 : (r <= 0.0) ? 1.0 : r2 * r2;
  return f;
}

int eliminate_outliers (DRMS_Array_t *img, double accept, double baseval) {
/*
 *  Eliminate values differing by a ratio of "accept" from the average of
 *    their neighbors.  Values of 0 are always accepted, likewise values
 *    for which the average of the neighbor values is within the range of
 *    noise "baseval"
 *  Bugs:
 *    edge pixels are ignored
 *    no account is taken of sign of value
 */
  double s, w;
  double amn, amx;
  float *z;
  float v, blank = 0.0 / 0.0;
  long long n;
  int cl[8];
  int col, cols, row, rows, i;
  int count = 0;
  if (accept <= 0.0) return -1;
  if (accept < 1.0) accept = 1.0 / accept;
  amx = accept;
  amn = 1.0 / accept;
  cols = drms_array_nth_axis (img, 0);
  rows = drms_array_nth_axis (img, 1);
  cl[0] = -cols - 1;
  cl[1] = -cols;
  cl[2] = -cols + 1;
  cl[3] = -1;
  cl[4] = 1;
  cl[5] = cols - 1;
  cl[6] = cols;
  cl[7] = cols + 1;

  z = (float *)img->data;
		/*  Go in reverse read-out order, since last pixel before gap
							   is often corrupt  */
  n = drms_array_count (img) - 1;
  n -= cols;
  for (row = rows - 2; row; row--) {
    for (col = cols - 2; col; col--, n--) {
      v = z[n];
      if (isnan (v)) continue;
      if (v == 0.0) continue;
      s = w = 0.0;
      for (i = 0; i < 8; i++) {
	v = z[n + cl[i]];
	if (isnan (v)) continue;
	w += 1.0;
	s += v;
      }
      if (w != 0.0) {
        v = fabs (z[n]);
	s /= w;
	s = fabs (s);
	if ((v > amx * s) && (v > baseval)) {
	  z[n] = blank;
	  count++;
	} else if ((s > baseval) && (v < amn * s)) {
	  z[n] = blank;
	  count++;
	}
      }
    }
    n -= 2;
  }
  return count;
}

int truncate_noise (DRMS_Array_t *img, double accept) {
  long long n = drms_array_count (img);
  float *z;
  int count = 0;

  z = (float *)img->data;
  while (n--) {
    if (fabs (*z) < accept) {
      *z = 0.0;
      count++;
    }
    z++;
  }
  return count;
}

int *ndimsegments (DRMS_Record_t *rec, int ndim, int *ct) {
/*  
 *  Returns a list of the segnums of the actually available segments (in case
 *    segments have been named in the dataset specification) of rank ndim
 *    in the record rec
 */
  DRMS_Record_t *temprec;
  DRMS_Segment_t *seg;
  static int *seglist = NULL;
  int found, n, segct, status;

  temprec = drms_template_record (drms_env, rec->seriesinfo->seriesname, &status);
  segct = drms_record_numsegments (temprec);
  if (!seglist) seglist = (int *)realloc (seglist, segct * sizeof (int));
  found = 0;
  for (n = 0; n < segct; n++) {
    seg = drms_segment_lookupnum (rec, n);
    if (!(seg = drms_segment_lookupnum (rec, n))) continue;
    if (seg->info->naxis == ndim) seglist[found++] = n;
  }
  *ct = found;
  return seglist;
}

int good_record (DRMS_RecordSet_t *ds, int nrt, int rstp, float max_reach,
    char *qual_key, unsigned int qmask, int *rejects, int *reject_list,
    int *shift) {
  DRMS_Record_t *rec;
  TIME t, tobs, trec;
  unsigned int quality;
  int idrec, match;
  int n, nr, nrr, nrmn, nrmx, offset, status;
  static int qcheck = 1;
  int found = 1;
  static unsigned char *ok = NULL;

  rec = ds->records[nrt];
  if (qcheck) {
	/*  check for data missing or otherwise unacceptable bits in quality  */
    quality = drms_getkey_int (rec, qual_key, &status);
    if (status) qcheck = 0;
    else if (quality & qmask) found = 0;
  }
  if (*rejects && found) {
    idrec = drms_getkey_int (rec, "T_REC_index", &status);
    match = 0;
    if (status) {
      fprintf (stderr, "Warning: \"T_REC_index\" keyword not found\n");
      fprintf (stderr, "         up to %d bad images could be processed\n",
	  *rejects);
      *rejects = 0;
    }
    for (n = 0; n < *rejects; n++) {
      if (idrec == reject_list[n]) {
	match = 1;
	break;
      }
    }
    if (match) found = 0;
  }
  if (found) {
    *shift = 0;
    return nrt;
  }
				   /*  target record unacceptable, try range  */
  nrr = -1;
  *shift = offset = max_reach * rstp;
  if (!offset) return nrr;
  nrmx = nrt + offset;
  nrmn = nrt - offset;
  ok = (unsigned char *)realloc (ok, ds->n * sizeof (char));
  if (nrmx >= ds->n) nrmx = ds->n - 1;
  if (nrmn < 0) nrmn = 0;
  for (nr = nrmn; nr <=  nrmx; nr++) {
    rec = ds->records[nr];
    if (qcheck) {
	/*  check for data missing or otherwise unacceptable bits in quality  */
      quality = drms_getkey_int (rec, qual_key, &status);
      if (status) qcheck = 0;
      else if (quality & qmask) {
	ok[nr] = 0;
	continue;
      }
    }
    if (*rejects) {
				   /*  check for inclusion in rejection list  */
      idrec = drms_getkey_int (rec, "T_REC_index", &status);
      match = 0;
      for (n = 0; n < *rejects; n++) {
	if (idrec == reject_list[n]) {
	  match = 1;
	  break;
	}
      }
      if (match) {
	ok[nr] = 0;
	continue;
      }
    }
    ok[nr] = 1;
  }
  for (nr = nrmn; nr <=  nrmx; nr++) {
    if (ok[nr] && (abs (nr - nrt) < *shift)) {
      nrr = nr;
      *shift = abs (nr - nrt);
    }
  }
  return nrr;
}

int drms_key_is_slotted (DRMS_Env_t *drms_env, const char *keyname,
    const char *dsname) {
  DRMS_Record_t *rec;
  DRMS_Keyword_t *key;
  int status;

  rec = drms_template_record (drms_env, dsname, &status);
  if (status) return 0;
  key = drms_keyword_lookup (rec, keyname, 1);
  if (!key) return 0;
  status = (key->info->recscope > 1);
  drms_free_keyword_struct (key);
  drms_free_record (rec);
  return status;
}

int get_cadence (DRMS_Record_t *rec, const char *source, const char *tstp_key,
    const char *trec_key, double *cadence) {
			      /*  attempt to determine uniform input cadence  */
/*
 *  This is a modified version of the similar function in mtrack; it neither
 *    uses nor set an output step size
 *  If the data cadence keyword tstp_key ("Cadence" by default) is present
 *    in the input series and is constant, use it
 *  If it is variable, or missing, and the time keyword trec_key ("T_REC" by
 *    default) exists and is slotted, use its step length as the cadence
 */
  DRMS_Keyword_t *keywd;
  int status;

  if ((keywd = drms_keyword_lookup (rec, tstp_key, 1))) {
					      /*  cadence should be constant  */
    *cadence = drms_getkey_double (rec, tstp_key, &status);
    if (keywd->info->recscope != 1) {
      fprintf (stderr, "Warning: %s is variable in input series %s\n",
	  tstp_key, source);
      drms_free_keyword_struct (keywd);
      if ((keywd = drms_keyword_lookup (rec, trec_key, 1))) {
	if (drms_key_is_slotted (drms_env, trec_key, source)) {
	  char *stepkey = malloc (strlen (trec_key) + 8);
	  sprintf (stepkey, "%s_step", trec_key);
	  *cadence = drms_getkey_double (rec, stepkey, &status);
	  free (stepkey);
	} else {
	  fprintf (stderr, "         and %s is not slotted\n", trec_key);
	  drms_free_keyword_struct (keywd);
	  return 1;
	}
	drms_free_keyword_struct (keywd);
      } else {
	fprintf (stderr, "         and %s is not present\n", trec_key);
	drms_free_keyword_struct (keywd);
	return 1;
      }
    }
  } else {
			     /* cadence key missing, use trec_key if slotted  */
    drms_free_keyword_struct (keywd);
    if ((keywd = drms_keyword_lookup (rec, trec_key, 1))) {
      if (drms_key_is_slotted (drms_env, trec_key, source)) {
	char *stepkey = malloc (strlen (trec_key) + 8);
	sprintf (stepkey, "%s_step", trec_key);
	*cadence = drms_getkey_double (rec, stepkey, &status);
	free (stepkey);
      } else {
	fprintf (stderr,
            "Error: data cadence keyword %s not in input series %s\n",
	    tstp_key, source);
	fprintf (stderr, "       and %s is not slotted\n", trec_key);
	return 1;
      }
      drms_free_keyword_struct (keywd);
    } else {
								 /*  give up  */
      fprintf (stderr,
          "Error: data cadence keyword %s not in input series %s\n",
	  tstp_key, source);
      fprintf (stderr, "       and %s is not present\n", trec_key);
      return 1;
    }
  }
  return 0;
}

TIME time_from_crcl (int cr, double cl, int from_soho) {
  if (from_soho)
    return SOHO_meridian_crossing (cl, cr);
  else
    return earth_meridian_crossing (cl, cr);
}

/******************************************************************************/
			/*  module body begins here  */
/******************************************************************************/
int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ds = NULL, *maids;
  DRMS_Record_t *rec;
  DRMS_Segment_t *seg;
  DRMS_Array_t *mgram;
  DRMS_Keyword_t *keywd;
  struct maploc {
    double lat;
    double lon;
    double wt;
  } *mloc;
  OBS_EPHEMERIS *obsval;
  TIME *tmin, *tmax;
  TIME tcm, tbase, tmid, tstrt, tstop, avgtobs;
  double *sum, *suma, *sum_wt, *area;
  double *sb, *sbt, *sf, *sft, *s1, *st, *st2, *det;
  double *f0, *f1, *fs0, *fs1;
  double *bmin, *bmax, *lonmin, *lonmax, *latmin, *latmax, *rmax;
  double *muavg, *xavg, *yavg;
  double lat, lon;
  double twt, ttwt, tfac, wtfac, img_weight;
  double sumt, sumtwt;
  double img_radius, ellipse_e, ellipse_pa;
  double x, y, x0, y0, xstp, ystp, mscale;
  double xc, yc, xscale, yscale, imgscale, latc, lonc, peff;
  double coslatc, sinlatc, cos_lat, cos_lon, cos_lat_lon, cos_cang, r;
  double data_cadence, t_eps, phase;
  float *clat, *clon;
  float b;
  int *seglist, *misspixct, *muct, *skipct, *skiprec;
  int *reject_list = NULL;
  int expct, nrec, recct, segct, valct;
  int mcol, mcols, mrow, mrows, map, mapn, maps, maprec, mapskip;
  int n, ntot, s, size;
  int plate_cols, plate_rows, plate_width;
  int x_invrt, y_invrt;
  int badqual, rejects, shifted, shiftct, tnegct, tposct, use;
  int usefit0, usefit1;
  int need_ephem, geo_times;
  int rank, off, count;
  int kstat, status;
  int MDI_correct, MDI_correct_distort, old_MDI_correct, old_MDI_correct_distort;
  int SDO_fixclon, SDO_fixcrot;
  short *adjust;
  char *short_time, *source_series, *los_series, *vec_series, *eoser;
  char rec_query[DRMS_MAXQUERYLEN];
  char reclist[8192], histbuf[128];
  char module_ident[64], time_str[64], key[64];
						    /*  constants & counters  */
  double raddeg = M_PI / 180.0;
  double degrad = 1.0 / raddeg;
  int projection = LAMBERT;
  int tmid_unknown = 1;
  platloc platform = LOC_UNKNOWN;
  char *mapname[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
      "LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
      "stereographic", "orthographic", "LambertAzimuthal"};
  missing_val = 0.0 / 0.0;

  char *losset = strdup (params_get_str (params, "los"));
  char *los_segname = strdup (params_get_str (params, "los_seg"));
  char *vecset = strdup (params_get_str (params, "vec"));
  char *maiser = strdup (params_get_str (params, "mai"));
  char *tmid_str = strdup (params_get_str (params, "tmid"));
  int tmid_cr = params_get_int (params, "cr");
  double tmid_cl = params_get_double (params, "cl");
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  int latct = params_get_int (params, "lat_nvals");
  int lonct = params_get_int (params, "lon_nvals");
  float map_scale = params_get_float (params, "scale");
  float map_size = params_get_float (params, "extent");
  double apode_inner = params_get_double (params, "mask_in");
  double apode_outer = params_get_double (params, "mask_ex");
  double apode_edge = params_get_double (params, "apodize");
  double noise_floor = params_get_double (params, "floor");
  int interval = params_get_int (params, "length");
  int rec_step = params_get_int (params, "rec_step");
  float max_reach = params_get_float (params, "max_reach");
  char *rejectfile = strdup (params_get_str (params, "reject"));
  char *tobs_key = strdup (params_get_str (params, "tobs_key"));
  char *trec_key = strdup (params_get_str (params, "trec_key"));
  char *tstp_key = strdup (params_get_str (params, "tstp_key"));
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *crot_key = strdup (params_get_str (params, "crot_key"));
  char *rsun_key = strdup (params_get_str (params, "rsun_key"));
  char *apsd_key = strdup (params_get_str (params, "apsd_key"));
  char *dsun_key = strdup (params_get_str (params, "dsun_key"));
  char *calverkey = strdup (params_get_str (params, "cvkey"));

  int fillbadrec = params_isflagset (params, "f");
  int listmais = params_isflagset (params, "l");
  int no_track = params_isflagset (params, "n");
  int no_apode = params_isflagset (params, "r");
  int despike = params_isflagset (params, "s") ? 0 : 1;
  int verbose = params_isflagset (params, "v");
  int MDI_proc = params_isflagset (params, "M");
  int old_MDI_proc = params_isflagset (params, "O");
  int SDO_fixkeys = (params_isflagset (params, "S")) ? 0 : 1;

  int uselos = strcmp (losset, "Not Specified");
  int usevec = strcmp (vecset, "Not Specified");
  int no_out = strcmp (maiser, "Not Specified") ? 0 : 1;

  int obs2recfix = strcasecmp (trec_key, tobs_key);

  if (MDI_proc && old_MDI_proc) {
    fprintf (stderr, "Error: -M or -O flag may be specified, but not both\n");
    return 1;
  }

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s: JSOC version %s\n", module_ident, jsoc_version);

  if (no_apode) apode_edge = 1.0;
  
  if (strcmp (tmid_str, "Not Specified")) {
    if (sscanf (tmid_str, "%d:%lf", &tmid_cr, &tmid_cl) != 2) {
      tmid = sscan_time (tmid_str);
      tmid_unknown = 0;
    }
  } else {
    if (tmid_cr < 0) {
		       /*  just approximate to previous carrington longitude  */
      double rsun, lat, lon, vr, vn, vw;
      tcm = CURRENT_SYSTEM_TIME;
      earth_ephemeris (tcm, &rsun, &lat, &lon, &vr, &vn, &vw);
      tmid_cr = carrington_rots (tcm, 0);
      tmid_cr--;
      if (isnan (tmid_cl)) tmid_cl = lon;
    } else if (isnan (tmid_cl)) tmid_cl = 0.0;
  }

  if (!uselos && !usevec) {
    fprintf (stderr, "Error: at least one of los and vec must be specified\n");
    return 1;
  }
  los_series = strdup (losset);
  if (uselos && key_params_from_dspec (losset)) {
    fprintf (stderr,
	"Warning: explicit data set specification is not supported\n");
    fprintf (stderr,
	"         record range will be ignored\n");
    eoser = strchr (los_series, '[');
    if (eoser) *eoser = '\0';
    eoser = strchr (los_series, '{');
    if (eoser) *eoser = '\0';
  }
  vec_series = strdup (vecset);
  if (usevec && key_params_from_dspec (vecset)) {
    fprintf (stderr,
	"Warning: explicit data set specification is not supported\n");
    fprintf (stderr,
	"         record range will be ignored\n");
    eoser = strchr (vec_series, '[');
    if (eoser) *eoser = '\0';
    eoser = strchr (vec_series, '{');
    if (eoser) *eoser = '\0';
  }
  if (usevec) {
    fprintf (stderr, "Warning: use of vector data is not supported!\n");
    usevec = 0;
    if (!uselos) return 1;
  }
  if (no_out) {
    if (verbose) {
      fprintf (stderr,
	"No output series specified, results will only be reported\n");
/*
      printf (" Lon   Lat       |Bz|   d|Bz|/dt    Bz     dBz/dt   FillFac [Gauss(/sec)]:\n");
*/
    }
    else listmais = 1;
  } else {
				      /*  check output data series structure  */
    rec = drms_template_record (drms_env, maiser, &status);
    if (status) {
      fprintf (stderr, "Error: unable to locate output series %s\n", maiser);
      return 1;
    }
    drms_close_record (rec, DRMS_FREE_RECORD);
  }
				       /*  determine data segment(s) to used  */
  if (uselos) {
    if (strcasecmp (los_segname, "Not Specified")) {
fprintf (stderr, "Sorry, segment name parsing is not implemented\n");
return 1;
		   /*  use template record to determine segment if necessary  */
    } else {
      rec = drms_template_record (drms_env, losset, &status);
      if (status || !rec) {
	fprintf (stderr, "Error: unable to open template record in series %s\n",
	    losset);
	return 1;
      }
      seglist = ndimsegments (rec, 2, &segct);
      if (segct > 1) {
	fprintf (stderr, "Warning: series %s contains %d segments of rank 2\n",
	    rec->seriesinfo->seriesname, segct);
	seg = drms_segment_lookupnum (rec, seglist[0]);
	fprintf (stderr, "         segment %d (%s) will be used\n", seglist[0],
	    seg->info->name);
      }
      segct = 1;
      drms_close_record (rec, DRMS_FREE_RECORD);
    }
    source_series = losset;
  }
		    /*  get required series info from first record in series  */
						/*  platform, cadence, phase  */
  snprintf (rec_query, 256, "%s[:#^]", source_series);
  if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
    fprintf (stderr,
	"Error: unable to open input data set \"%s\"\n", rec_query);
    fprintf (stderr, "       status = %d\n", status);
    return 1;
  }
  rec = ds->records[0];
  if (get_cadence (rec, source_series, tstp_key, trec_key, &data_cadence)) {
    drms_close_records (ds, DRMS_FREE_RECORD);
    return 1;
  }
  t_eps = 0.5 * data_cadence;
  geo_times = 1;
  if ((keywd = drms_keyword_lookup (rec, "TELESCOP", 1))) {
						     /*  should be constant  */
    if (keywd->info->recscope != 1)
      fprintf (stderr, "Warning: TELESCOP is variable in input series %s\n",
	  source_series);
    if (!strcmp (drms_getkey_string (rec, "TELESCOP", &status), "SDO/HMI"))
      platform = LOC_SDO;
    else if (!strcmp (drms_getkey_string (rec, "TELESCOP", &status), "SOHO")) {
      platform = LOC_SOHO;
      geo_times = 0;
    }
  }
  if (platform == LOC_UNKNOWN) fprintf (stderr,
	"Warning: observing location unknown, assumed geocenter\n");
  if (platform != LOC_SDO) SDO_fixkeys = 0;
  SDO_fixclon = SDO_fixcrot = SDO_fixkeys;
  if (tmid_unknown) {
			   /*  tmid specified as CR:CL : need ephemeris info  */
    tmid = (geo_times) ?
	time_from_crcl (tmid_cr, tmid_cl, 0) :
	time_from_crcl (tmid_cr, tmid_cl, platform == LOC_SOHO);
  };
			   /*  adjust to phase of input, within data cadence  */
  tbase = drms_getkey_time (rec, trec_key, &status);
  phase = fmod ((tmid - tbase), data_cadence);
  tcm = tmid;
  tmid -= phase;
  tstrt = tmid - 0.5 * interval * data_cadence;
  tstop = tstrt + (interval - 1) * data_cadence;
  drms_close_records (ds, DRMS_FREE_RECORD);
		 /*  determine appropriate range for and open input data set  */
  if (drms_key_is_slotted (drms_env, trec_key, source_series)) {
    DRMS_Record_t *rec = drms_template_record (drms_env, source_series,
	&status);
    TIME pkeybase;
    double pkeystep;
    int tstrt_ind, tstop_ind;
    char *pkindx = malloc (strlen (trec_key) + 8);
    char *pkbase = malloc (strlen (trec_key) + 8);
    char *pkstep = malloc (strlen (trec_key) + 8);

    sprintf (pkindx, "%s_index", trec_key);
    sprintf (pkbase, "%s_epoch", trec_key);
    sprintf (pkstep, "%s_step", trec_key);
    pkeybase = drms_getkey_time (rec, pkbase, &status);
    pkeystep = drms_getkey_double (rec, pkstep, &status);
    tstrt_ind = (tstrt - pkeybase) / pkeystep;
    tstop_ind = (tstop - pkeybase) / pkeystep;
    snprintf (rec_query, 256, "%s[?%s between %d and %d?]", source_series,
	pkindx, tstrt_ind, tstop_ind);

    free (pkindx);
    free (pkbase);
    free (pkstep);
    drms_free_record (rec);
  } else snprintf (rec_query, 256, "%s[?%s between %23.16e and %23.16e?]",
      source_series, trec_key, tstrt - t_eps, tstop + t_eps);
  if (!(ds = drms_open_records (drms_env, rec_query, &status))) {
    fprintf (stderr, "Error: unable to open input data set \"%s\\n", rec_query);
    fprintf (stderr, "       status = %d\n", status);
    drms_close_records (ds, DRMS_FREE_RECORD);
    return 1;
  }
  if ((recct = ds->n) < 2) {
    fprintf (stderr, "Error: <2 records in selected input set:\n");
    fprintf (stderr, "       %s with selection\n", source_series);
    fprintf (stderr, "       %s\n", rec_query);
    drms_close_records (ds, DRMS_FREE_RECORD);
    return 1;
  }
		  /*  support special hack of reading of rejection list file  */
  rejects = 0;
  if (strcmp (rejectfile, "Not Specified")) {
    FILE *rejectfp = fopen (rejectfile, "r");
    if (rejectfp) {
      rejects = read_reject_list (rejectfp, &reject_list);
      fclose (rejectfp);
    }
    else fprintf (stderr,
	"Warning: could not open rejection list %s; ignored\n", rejectfile);
  }
					   /*  allocate tile tracking arrays  */
  maps = (latct > lonct) ? latct : lonct;
  clat = (float *)malloc (maps * sizeof (float));
  clon = (float *)malloc (maps * sizeof (float));
  for (map = 0; map < latct; map++) {
    snprintf (key, 64, "lat_%d_value", map);
    clat[map] = params_get_float (params, key) * raddeg;
  }
  for (map = 0; map < lonct; map++) {
    snprintf (key, 64, "lon_%d_value", map);
    clon[map] = params_get_float (params, key) * raddeg;
  }
					      /*  extend arrays as necessary  */
  for (map = latct; map < maps; map++) clat[map] = clat[map-1];
  for (map = lonct; map < maps; map++) clon[map] = clon[map-1];
  reclist[0] = '\0';

  tfac = 2.0 / (data_cadence * interval);
					     /*  Initialize target locations  */
  mcols = mrows = map_size / map_scale + 0.5;
  size = mcols * mrows;
  xstp = ystp = map_scale * raddeg;
  x0 = 0.5 * (1.0 - mcols) * xstp;
  y0 = 0.5 * (1.0 - mrows) * ystp;
  mscale = (map_size > map_scale) ?
      2.0 / (map_size - map_scale) / raddeg : 0.0;
			      /*  pixel area in cm^2: gives flux in Maxwells  */
/*
  pxl_area = xstp * ystp * RSUNCM * RSUNCM;
*/
  old_MDI_correct = old_MDI_correct_distort = old_MDI_proc;
  MDI_correct = MDI_correct_distort = MDI_proc;
  ntot = mrows * mcols;
  ntot *= maps;
  mloc = (struct maploc *)malloc (ntot * sizeof (struct maploc));
  f0 = (double *)malloc (maps * sizeof (double));
  f1 = (double *)malloc (maps * sizeof (double));
  fs0 = (double *)malloc (maps * sizeof (double));
  fs1 = (double *)malloc (maps * sizeof (double));
  sum = (double *)malloc (maps * sizeof (double));
  suma = (double *)malloc (maps * sizeof (double));
  sum_wt = (double *)calloc (maps, sizeof (double));
  area = (double *)calloc (maps, sizeof (double));
  det = (double *)malloc (maps * sizeof (double));
  s1 = (double *)calloc (maps, sizeof (double));
  st = (double *)calloc (maps, sizeof (double));
  st2 = (double *)calloc (maps, sizeof (double));
  sb = (double *)calloc (maps, sizeof (double));
  sbt = (double *)calloc (maps, sizeof (double));
  sf = (double *)calloc (maps, sizeof (double));
  sft = (double *)calloc (maps, sizeof (double));
  bmin = (double *)malloc (maps * sizeof (double));
  bmax = (double *)malloc (maps * sizeof (double));
  lonmin = (double *)malloc (maps * sizeof (double));
  lonmax = (double *)malloc (maps * sizeof (double));
  latmin = (double *)malloc (maps * sizeof (double));
  latmax = (double *)malloc (maps * sizeof (double));
  rmax = (double *)calloc (maps, sizeof (double));
  misspixct = (int *)calloc (maps, sizeof (int));
  skipct = (int *)calloc (maps, sizeof (int));
  skiprec = (int *)calloc (maps, sizeof (int));
  tmin = (TIME *)malloc (maps * sizeof (TIME));
  tmax = (TIME *)malloc (maps * sizeof (TIME));
  muavg = (double *)calloc (maps, sizeof (double));
  muct = (int *)calloc (maps, sizeof (int));
  xavg = (double *)calloc (maps, sizeof (double));
  yavg = (double *)calloc (maps, sizeof (double));
  adjust = (short *)calloc (maps, sizeof (short));
	/*  calculate spatial apodization weights for each pixel in each map  */
  for (n = 0, map = 0; map < maps; map++) {
    for (mrow = 0, y = y0; mrow < mrows; mrow++, y += ystp) {
      for (mcol = 0, x = x0; mcol < mcols; mcol++, x += xstp, n++) {
	off = plane2sphere (x, y, clat[map], clon[map], &lat, &lon, projection);
	mloc[n].lat = lat;
	mloc[n].lon = lon;
	mloc[n].wt = (off) ? 0.0 : (no_apode) ? 1.0 :
	    apodization (x * mscale, y * mscale, apode_inner, apode_outer);
	sum_wt[map] += mloc[n].wt;
      }
    }
  }
							   /*  and normalize  */
  for (map = 0, n = 0; map < maps; map++) {
    if (sum_wt[map] != 0.0) {
      wtfac = 1.0 / sum_wt[map];
      for (s = 0; s < size; s++, n++)
	mloc[n].wt *= wtfac;
    } else n += size;
    bmax[map]= -HUGE_VAL;
    bmin[map] = HUGE_VAL;
  }
  valct = 0;
  sumt = sumtwt = 0.0;
  tnegct = tposct = 0;
  usefit0 = usefit1 = 1;
		   /*  adjust ephemeris values for mapping from tobs to trec  */
  if (obs2recfix) {
    obsval = (OBS_EPHEMERIS *)malloc (recct * sizeof (OBS_EPHEMERIS));
    adjust_obs2rec (obsval, ds, recct, trec_key, tobs_key, crot_key, clon_key,
        clat_key, dsun_key, apsd_key, tstrt, tstop, data_cadence, source_series);
  }
					      /*  loop through input records  */
  for (nrec = rec_step / 2; nrec < recct; nrec += rec_step) {
    TIME tobs, trec;

    use = good_record (ds, nrec, rec_step, max_reach, qual_key, qmask,
	&rejects, reject_list, &shifted);
    if (use < 0) continue;
    rec = ds->records[use];
    trec = drms_getkey_time (rec, trec_key, &status);
    tobs = drms_getkey_time (rec, tobs_key, &status);
    if (tobs < tcm) tnegct++;
    else tposct++;
    seg = drms_segment_lookupnum (rec, seglist[0]);
    mgram = drms_segment_read (seg, DRMS_TYPE_FLOAT, &status);
    if (!mgram || status) continue;
    if ((rank = drms_array_naxis (mgram)) != 2) {
      fprintf (stderr, "improper format for record %d: rank = %d\n", use,
	  rank);
      drms_free_array (mgram);
      continue;
    }
    plate_cols = drms_array_nth_axis (mgram, 0);
    plate_rows = drms_array_nth_axis (mgram, 1);
			 /*  time-weight value for time derivs: from -1 to 1  */
    twt = tfac * (tobs - tcm);
    ttwt = fabs (twt);
						    /*  temporal apodization  */
    if (apode_edge < 1.0) {
      if (ttwt > apode_edge) {
	double t = (ttwt - apode_edge) / (1.0 - apode_edge);
	if (t > 1.0) t = 1.0;
	t *= t;
	ttwt = (1.0 - t);
	ttwt *= ttwt;
      } else ttwt = 1.0;
    } else ttwt = 1.0;
    if (ttwt <= 0.0) continue;
    count = drms_getkey_int (rec, "DATAVALS", &status);
    if (status) count = plate_cols * plate_rows;
			    /*  get needed info from record keys for mapping  */
    status = solar_image_info (rec, &xscale, &yscale, &xc, &yc,
	&img_radius, rsun_key, apsd_key, &peff, &ellipse_e, &ellipse_pa,
	&x_invrt, &y_invrt, &need_ephem, 0);
    if (status & NO_SEMIDIAMETER) {
      int keystat = 0;
      double dsun_obs = drms_getkey_double (rec, dsun_key, &keystat);
      if (keystat) {
	fprintf (stderr, "Error: one or more essential keywords or values missing; skipped\n");
	fprintf (stderr, "solar_image_info() returned %08x\n", status);
	continue;
      }
				/*  set image radius from scale and distance  */
      img_radius = asin (RSUNM / dsun_obs);
      img_radius *= 3600.0 * degrad;
      img_radius /= (xscale <= yscale) ? xscale : yscale;
      status &= ~NO_SEMIDIAMETER;
    }
    if (status) {
      fprintf (stderr, "Error: one or more essential keywords or values missing; skipped\n");
      fprintf (stderr, "solar_image_info() returned %08x\n", status);
      continue;
    }
    if (SDO_fixcrot) {
      keywd = drms_keyword_lookup (rec, calverkey, 1);
      if (keywd) {
        long long calver = drms_getkey_longlong (rec, calverkey, &status);
	if (calver & OBS_CROT_FIXED) ;
        else peff -= 0.00702 * raddeg;
      } else {
        fprintf (stderr, "Warning: required keyword %s not found\n", calverkey);
	fprintf (stderr, "         no calibration version filtering applied\n");
	fprintf (stderr, "         %s and %s values assumed valid\n", clon_key,
	    crot_key);
	SDO_fixclon = SDO_fixcrot = SDO_fixkeys = 0;
      }
    }
    if (obs2recfix) img_radius = obsval[use].apsd / xscale;
    if (old_MDI_correct) {
      mtrack_MDI_correct_imgctr (&xc, &yc, img_radius);
      mtrack_MDI_correct_pa (&peff);
    } else if (MDI_correct) {
			     /*  not sure whether this is appropriate or not  */
      mtrack_MDI_correct_pa (&peff);
    }
			       /*  append to list of records used in history  */
    sprint_time (time_str, tobs, "TAI", 0);
    time_str[19] = '\0';
    if (strlen (reclist) < 8170) {
      if (strlen (reclist)) {
	strcat (reclist, ", ");
	short_time = time_str + 8;
	strcat (reclist, short_time);
      } else strcat (reclist, time_str);
    }
    if (obs2recfix) {
      lonc = raddeg * obsval[use].clon;
      latc = raddeg * obsval[use].clat;;
    } else {
      lonc = drms_getkey_double (rec, clon_key, &status) * raddeg;
      latc = drms_getkey_double (rec, clat_key, &status) * raddeg;
    }
    if (SDO_fixclon) {
      keywd = drms_keyword_lookup (rec, calverkey, 1);
      if (keywd) {
        long long calver = drms_getkey_longlong (rec, calverkey, &status);
	if (calver & OBS_CRLN_FIXED) ;
        else lonc -= 0.081894 * raddeg;
      } else {
        fprintf (stderr, "Warning: required keyword %s not found\n", calverkey);
	fprintf (stderr, "         no calibration version filtering applied\n");
	fprintf (stderr, "         %s and %s values assumed valid\n", clon_key,
	    crot_key);
	SDO_fixclon = SDO_fixcrot = SDO_fixkeys = 0;
      }
    }
    coslatc = cos (latc);
    sinlatc = sin (latc);
    if (no_track) lonc = tmid_cl * raddeg;
    imgscale = (xscale > yscale) ? xscale : yscale;
						 /*  Truncate "noise" values  */
    status = truncate_noise (mgram, noise_floor);
    count -= status;
    if (despike) eliminate_outliers (mgram, OUTLIER_RATIO, OUTLIER_BASE);
				      /*  form average for effective obstime  */
    sumt += (tobs - tcm) * ttwt;
    sumtwt += ttwt;

    plate_width = (plate_cols > plate_rows) ? plate_cols : plate_rows;
    xc -= 0.5 * (plate_cols - 1);
    yc -= 0.5 * (plate_rows - 1);
    xc *= 2.0 / plate_width;
    yc *= 2.0 / plate_width;
    img_radius *= 2.0 / plate_width;
				     /*  map field values for each location  */
    for (map = 0, n = 0; map < maps; map++) {
      double bigwt, totwt;
      sum[map] = suma[map] = 0.0;
      ntot = 0;
      totwt = bigwt = 0.0;
      img_weight = ttwt;
      for (s = 0; s < size; s++, n++) {
	if (mloc[n].wt == 0.0) {
	  skipct[map]++;
	  continue;
	}
	if (sphere2img (mloc[n].lat, mloc[n].lon, latc, lonc, &x, &y, xc, yc,
	      img_radius, peff, ellipse_e, ellipse_pa, x_invrt, y_invrt)) {
	  img_weight -= mloc[n].wt * ttwt;
	  misspixct[map]++;
	  continue;
	}
	if (old_MDI_correct_distort) {
	  mtrack_MDI_image_tip (&x, &y, 1, 1);
	  mtrack_MDI_image_stretch (&x, &y, 1, 1);
	} else if (MDI_correct_distort) {
	  mtrack_MDI_image_distort (&x, &y, 1, 1);
	}
	b = array_imaginterp (mgram, x, y, 0);
	if (isnan (b)) {
	  img_weight -= mloc[n].wt * ttwt;
	  misspixct[map]++;
	  continue;
	}
	if (b > bmax[map]) {
	  bmax[map] = b;
	  lonmax[map] = mloc[n].lon;
	  latmax[map] = mloc[n].lat;
	  tmax[map] = trec;
	}
	if (b < bmin[map]) {
	  bmin[map] = b;
	  lonmin[map] = mloc[n].lon;
	  latmin[map] = mloc[n].lat;
	  tmin[map] = trec;
	}
						 /*  populate averaging sums  */
	sum[map] += b * mloc[n].wt;
	suma[map] += fabs (b) * mloc[n].wt;
	if (b != 0.0) bigwt += mloc[n].wt;
	totwt +=  mloc[n].wt;
	ntot++;
				/*  calculate map pixel distance from center  */
	cos_lat = cos (mloc[n].lat);
	cos_lon = cos (mloc[n].lon - lonc);
	cos_lat_lon = cos_lat * cos_lon;
	cos_cang  = sin (mloc[n].lat) * sinlatc + coslatc * cos_lat_lon;
	r = (cos_cang < 0.0) ? 1.0 : sqrt (1.0 - cos_cang * cos_cang);
	if (r > rmax[map]) rmax[map] = r;
        if (cos_cang < 0.0) continue;
	muct[map]++;
	muavg[map] += cos_cang;
	xavg[map] += x - xc;
	yavg[map] += y - yc;
      }
      if (totwt > 0.0) area[map] += img_weight * bigwt / totwt;
      if (img_weight < 1.0) {
	sum[map] *= img_weight;
	suma[map] *= img_weight;
      }
      s1[map] += img_weight;
      st[map] += img_weight * twt;
      st2[map] += img_weight * twt * twt;
      sb[map] += sum[map];
      sbt[map] += twt * sum[map];
      sf[map] += suma[map];
      sft[map] += twt * suma[map];
    }
/*
sprint_time (time_str, tobs, "TAI", -1);
printf ("%s %.4f %.4f %.4f %.4f %.4f\n", time_str, twt, suma[0], sum[0], sf[0], sb[0]);
*/
    drms_free_array (mgram);
    valct++;
  }
  drms_close_records (ds, DRMS_FREE_RECORD);
  if (valct < 2) {
    fprintf (stderr, "Error: <2 usable magnetograms in interval\n");
    return 1;
  }
	/*  trap cases with no valid data either before or after target time  */
  if (!tnegct || !tposct) usefit0 = 0;
		     /*  trap cases where linear fit should not be attempted  */
  if (valct < 3) {
    usefit0 = usefit1 = 0;
  }
  expct = valct * size;
  sumt /= sumtwt;
  avgtobs = tcm + sumt;

  if (verbose)
    printf (" Lon   Lat        MAI    dMAIdt    <Bz>   d<Bz>dt    Area\n");
  for (map = 0; map < maps; map++) {
    det[map] = 1.0 / (s1[map] * st2[map] - st[map] * st[map]);
    area[map] /= s1[map];
  }
  mapskip = 0;
  for (map = 0; map < maps; map++) {
    f0[map] = det[map] * (st2[map] * sf[map] - st[map] * sft[map]);
    f1[map] = det[map] * tfac * (s1[map] * sft[map] - st[map] * sf[map]);
    fs0[map] = det[map] * (st2[map] * sb[map] - st[map] * sbt[map]);
    fs1[map] = det[map] * tfac * (s1[map] * sbt[map] - st[map] * sb[map]);
    if (!usefit0) {
      f0[map] = sf[map] / valct;
      fs0[map] = sb[map] / valct;
    }
    if (!usefit1) {
      f1[map] = fs1[map] = 0.0 / 0.0;
    }
    if (misspixct[map] + skipct[map] >= expct) {
      f0[map] = f1[map] = fs0[map] = fs1[map] = 0.0 / 0.0;
      lonmin[map] = lonmax[map] = latmin[map] = latmax[map] = 0.0 / 0.0;
      bmin[map] = bmax[map] = tmax[map] = rmax[map] = area[map] = 0.0 / 0.0;
      if (!fillbadrec) {
					      /*  delete this map for output  */
	mapskip++;
	skiprec[map] = 1;
      }
    }
/*
printf ("sum_wt = %f, sum(B) = %f, sum(|B|) = %f\n", s1[0], sb[0], sf[0]);
printf ("det = %f, st2*sf = %f, st*sft = %f\n", det[0], st2[0]*sf[0], st[0]*sft[0]);
printf ("sf = %f, sb = %f\n", sf[0], sb[0]);
*/
	   /*  trap remaining cases with negative intercepts for linear fits  */
    if (f0[map] < 0.0) {
      fprintf (stderr,
	  "%d:%03.0f %05.1f%+05.1f MAI adjusted from %.3f to %.3f\n",
	  tmid_cr, tmid_cl, clon[map] * degrad, clat[map]  * degrad, f0[map],
	  sf[map] / valct);
      if (verbose) {
        printf ("%05.1f %+05.1f : %8.3f %8.3f %8.3f %8.3f %9.4f\n",
	  clon[map] * degrad, clat[map]  * degrad, sf[map] / valct,
	  1.0e6 * f1[map], fs0[map], 1.0e6 * fs1[map], area[map]);
      } else if (listmais) printf ("%8.3f", sf[map] / valct);
      adjust[map] = 1;
    } else {
      if (verbose) {
	printf ("%05.1f %+05.1f : %8.3f %8.3f %8.3f %8.3f %9.4f\n",
	    clon[map] * degrad, clat[map]  * degrad, f0[map], 1.0e6 * f1[map],
	    fs0[map], 1.0e6 * fs1[map], area[map]);
      } else if (listmais) printf ("%8.3f", f0[map]);
    }
  }
  if (listmais && !verbose) printf ("\n");
  if (mapskip == maps) {
    fprintf (stderr, "Warning: no valid records to populate\n");
    no_out = 1;
  }
  if (no_out) return 0;
				/*  create and populate output series records  */
  if (!(maids = drms_create_records (drms_env, maps - mapskip, maiser,
      DRMS_PERMANENT, &status))) {
    fprintf (stderr, "Error: unable to create %d records in series %s\n",
	maps - mapskip, maiser);
    fprintf (stderr, "       drms_create_records() returned status %d\n",
	status);
    return 1;
  }
						/*  set output record values  */
  kstat = 0;
  for (maprec = map = 0; map < maps; map++) {
    float loncm = clon[map] * degrad - tmid_cl;
    float lonminhg = lonmin[map] * degrad;
    float lonmaxhg = lonmax[map] * degrad;
    if (skiprec[map]) continue;
    while (loncm > 180.0) loncm -= 360.0;
    while (loncm < -180.0) loncm += 360.0;
    while (lonminhg < 0.0) lonminhg += 360.0;
    while (lonminhg > 360.0) lonminhg -= 360.0;
    while (lonmaxhg < 0.0) lonmaxhg += 360.0;
    while (lonmaxhg > 360.0) lonmaxhg -= 360.0;
    rec = maids->records[maprec];
    kstat += check_and_set_key_int (rec, "CarrRot", tmid_cr);
    kstat += check_and_set_key_float (rec, "CMLon", tmid_cl);
    kstat += check_and_set_key_float (rec, "LonHG", clon[map] * degrad);
    kstat += check_and_set_key_float (rec, "LatHG", clat[map] * degrad);
    kstat += check_and_set_key_float (rec, "LonCM", loncm);
    if (adjust[map]) {
      char comment[80];
      sprintf (comment, "MAI adjusted from %.3f to %.3f\n", f0[map],
          sf[map] / valct);
      kstat += drms_appendstr_tokey (rec, "COMMENT", comment, 1);
      f0[map] = sf[map] / valct;
    }
    kstat += check_and_set_key_float (rec, "MAI", f0[map]);
    if (usefit1)
      kstat += check_and_set_key_float (rec, "dMAIdt", 1.0e6 * f1[map]);
    kstat += check_and_set_key_float (rec, "MeanBigBz", fs0[map]);
    if (usefit1)
      kstat += check_and_set_key_float (rec, "dMnBigBdt", 1.0e6 * fs1[map]);
    kstat += check_and_set_key_float (rec, "BigBArea", area[map]);
    kstat += check_and_set_key_float (rec, "BMin", bmin[map]);
    kstat += check_and_set_key_float (rec, "LonBMin", lonminhg);
    kstat += check_and_set_key_float (rec, "LatBMin", latmin[map] * degrad);
    kstat += check_and_set_key_time (rec, "t_BMin", tmin[map]);
    kstat += check_and_set_key_float (rec, "BMax", bmax[map]);
    kstat += check_and_set_key_float (rec, "LonBMax", lonmaxhg);
    kstat += check_and_set_key_float (rec, "LatBMax", latmax[map] * degrad);
    kstat += check_and_set_key_time (rec, "T_OBS", avgtobs);
    kstat += check_and_set_key_time (rec, "t_BMax", tmax[map]);
    kstat += check_and_set_key_float (rec, "rMax", rmax[map]);
    kstat += check_and_set_key_int (rec, "Samples", valct);
    sprintf (histbuf, "%s magnetograms used:", source_series);
    kstat += drms_appendstr_tokey (rec, "HISTORY", histbuf, 1);
    kstat += drms_appendstr_tokey (rec, "HISTORY", reclist, 1);
    kstat += check_and_set_key_float(rec, "BadPix", (double)misspixct[map] / expct);
    kstat += check_and_set_key_str (rec, "MapProj", mapname[projection]);
    kstat += check_and_set_key_float (rec, "MapScale", map_scale);
    kstat += check_and_set_key_float (rec, "Size", map_size);
    kstat += check_and_set_key_float (rec, "BFloor", noise_floor);
    kstat += check_and_set_key_float (rec, "Apode_t", apode_edge);
    kstat += check_and_set_key_float (rec, "Apod_Min", apode_inner);
    kstat += check_and_set_key_float (rec, "Apod_Max", apode_outer);
    if (muct[map]) {
      kstat +=
	  check_and_set_key_float (rec, "MeanMu", muavg[map] / muct[map]);
      xavg[map] /= muct[map];
      yavg[map] /= muct[map];
      kstat += check_and_set_key_float (rec, "MeanPA",
	  degrad * atan2 (xavg[map], yavg[map]));
    }
    kstat += check_and_set_key_str (rec, "Source", source_series);
    kstat += check_and_set_key_time (rec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str (rec, "Module", module_ident);
    kstat += check_and_set_key_str (rec, "BLD_VERS", jsoc_version);
    maprec++;
  }
  if (kstat) {
    fprintf (stderr, "Error: unable to set all key values in %s\n", maiser);
    drms_close_records (maids, DRMS_FREE_RECORD);
    return 1;
  }
  drms_close_records (maids, DRMS_INSERT_RECORD);								/*  clean up  */
  return 0;
}
/******************************************************************************/
/*
 *  Revision history
 *  all mods by Rick Bogart unless otherwise noted
 *
 *  2016.10.26	Rick Bogart	created this file, based on v 1.1 of maicalc.c
 *  	For previous revision history, see
 *	~rick/hmi/rings/src/versions/maicalc/v11.c
 *  2017.08.03	Added reporting of maximum radial distance of input data,
 *	fraction of magnetogram pixels skipped due to mapping or despiking;
 *	Moved reporting of magnetograms used after validity checks;
 *	Changed default behavior to include despiking (as originally)
 *  2017.09.14	Added option of specifying tmid as time (date_time or
 *`	Carrington) instead of CR and CL separately;
 *	Modified algorithm for determination of sampling interval to agree
 *	with that for mtrack, requiring change of argument from "interval"
 *	to "length";
 *	Implemented temporal apodization weighting to agree with the form
 *	for pspec3, requiring additional optional argument "apodize", with
 *	same default as for pspec3; note that the default for mask_in is not
 *	the same as that for pspec3, but rather the value used in the HMI
 *	rings pipeline
 *	Added setting of effective (sampling-weighted) T_OBS key in output
 *	series, also keys for spatial and temporal apodization limits
 *  2017.09.26	Added verbose header
 *	Added optional -l flag for just listing. for consistency with output
 *	of earlier versions
 *  v 2.0 frozen 2017.09.25
 *      Note: Binaries built prior to 2017.10.11 had version ID 2.08x
 *  2017.11.02	Fixed minor bug affecting area determination when partial
 *  	magnetograms contribute
 *  2017.11.13	Added -f option to force filling of records with no valid data
 *	with NaN values, otherwise skip insertion of these records
 *		Added -O flag for including (old) MDI distortion model
 *		Added -M flag for including new MDI distortion model
 *  2018.01.25	Close reject list when done
 *		Added argument tobs_key; consistent use of tobs for
 *	time-weighting and trend fitting
 *		Fixed bug causing non-zero (and potentially large) weights
 *	assigned to time samples beyond the (assumed) outer temporal apodization
 *	edge
 *		Modified so that when all valid time samples are on one side
 *	of the interval from the target time, use the unweighted means of Bz
 *	and |Bz| for the MAI rather than the constant term in the linear fit
 *	(an extrapolation to the target time)
 *		Modified so that linear fit in time is not attempted at all
 *	when there are fewer than 3 valid time samples (i.e. 2, since no fit
 *	results are returned when there is only a single sample
 *		Modified so that any remaining negative MAI values are converted
 *	to the unweighted means of |Bz|, but with no other changes to reported
 *	parameters in those cases
 *		Added setting of keywords Samples, MeanMu, MeanPA, COMMENT
 *  v 2.1 frozen 2018.02.09
 *  2018.07.24	Fixed "bug" in setting of Apodization limits keywords -
 *  	keyword names corrected to those in use
 *  2018.09.20	Added check for validity of input data series  
 *  v 2.2 frozen 2020.03.17
 *  2020.09.07	Added calibration version-sensitive correction of SDO CMLon
 *	and position angle, flag for suppression of corrections
 *  v 2.3 frozen 2020.09.22
 *  2020.12.18	Fixed possible bug introduced in previous version for
 *	correction of SDO CMLon (but probably never reached)
 *  2021.01.08	Fixed additional possible bug introduced in previous version for
 *	correction of SDO CMLon (but probably never reached)
 *  v 2.4 frozen 2021.02.15
 *		
 */
/******************************************************************************/
