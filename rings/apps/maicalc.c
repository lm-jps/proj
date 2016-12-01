/******************************************************************************/
/*
 *  maicalc.c						~rick/hmi/rings/src/
 *
 *  Integrate functions of values in a selected set of regions from a set
 *    of magnetograms after mapping in order to form time-region indices
 *
 *  Responsible:
 *      Rick Bogart					rick@sun.stanford.edu
 *
 *  Parameters: (type   default         description)
 *	ds	str	-		input data series name
 *	cr	int	current		Carrington rotation
 *	cl	float	0 or current	Carrington longitude
 *	interval float	-		length of sampling interval [deg]
 *	rec_step int	64		sampling step size [rec]
 *	scale	float	-		map scale in heliographic deg / pixel
 *	extent	float	15.36		map extent in degrees
 *	lat	float*	[0.0]		Target heliographic latitude(s)
 *	lon	float*	[0.0]		Target heliographic longitude(s)
 *	mask_in	float	0.9765625	inner edge for 1-r^2 apodization
 *	mask_ex	float	1.000		outer edge for 1-r^2 apodization
 *	floor	float	50.0		"noise" floor for inclusion in integration [gauss]
 *	reject	string	-		If specified, name of a file with a
 *				list of input images to be rejected regardless
 *	qmask	int	0x80000000	Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	qual_key string	Quality		Keyname of int type keyword describing
 *				record quality as bitmask
 *	clat_key string	CRLT_OBS
 *	clon_key string	CRLN_OBS
 *	crot_key string	CAR_ROT
 *	rsun_key string	R_SUN
 *	apsd_key string	RSUN_OBS
 *	dsun_key string	DSUN_OBS	Keyname of double type keyword for
 *				distance from sun in m for of each image
 *	max_reach float	0.5		maximum distance from target sample
 *				to search when target is unaaceptable, in units
 *				of sample spacing
 *
 *  Flags
 *	-n	turn off tracking
 *	-r	turn off apodization
 *	-t	include time range check (for backward compatability)
 *	-v	run verbose
 *
 *  Bugs:
 *    No check for multiple segment, just reads segment 0 unconditionally
 *    There is no correction for foreshortening
 *    The query used for opening records by the function
 *`	select_dataset_from_time_interval() is inefficient
 *
 *  Possible Improvements:
 *
 *  Revision history is at end of file.
 */
/******************************************************************************/

#include <jsoc_main.h>

char *module_name = "maicalc";
char *module_desc = "integration of mapped and tracked magnetogram data";
char *version_id = "1.1";

#define RSUNM		(6.96e8)

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds", "", "input data series or dataset"}, 
  {ARG_INT,	"cr", "Not Specified", "Carrington Rotation (default: current)"},
  {ARG_FLOAT,	"cl", "Not Specified",
      "Carrington longitude of central meridian (default: current or 180"},
  {ARG_FLOAT,   "interval", "", "length of sampling interval, Carr. deg",
      "(0.0,)"},
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_INT,     "rec_step", "", "sampling step, record steps", "[1,)"},
  {ARG_FLOAT,	"max_reach", "0.5", "maximum reach to acceptable record",
      "[0.0,1.0]"},
  {ARG_FLOAT,	"scale", "",
      "integrating map scale for interpolation [deg/pixel]", "(0.0,)"},
  {ARG_FLOAT,	"extent", "", "extent of integrating maps [degrees]"},
  {ARG_FLOATS,	"lat", "[0.0]",
      "heliographic latitude(s) of tracking center(s) [deg]"},
  {ARG_FLOATS,	"lon", "[0.0]",
      "heliographic longitude(s) of tracking center(s) [deg]"},
  {ARG_FLOAT,   "mask_in", "0.9765625", "inner edge of 1 - r^2 apodization",
      "[0.0,)"},
  {ARG_FLOAT,   "mask_ex", "1.0", "outer edge of 1 - r^2 apodization", "[0.0,)"},
  {ARG_FLOAT,   "floor", "50.0", "noise floor (gauss)", "[0.0,)"},
  {ARG_STRING,  "file", "Not Specified", "output file name for diagnostics"},
  {ARG_STRING,	"qual_key", "Quality",  "keyname for 32-bit image quality field"}, 
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"crot_key", "CAR_ROT", "keyname for image Carrington rotation"}, 
  {ARG_STRING,	"rsun_key", "R_SUN", "keyname for image semi-diameter (pixel)"}, 
  {ARG_STRING,	"apsd_key", "RSUN_OBS",
      "keyname for apparent solar semi-diameter (arcsec)"}, 
  {ARG_STRING,	"dsun_key", "DSUN_OBS", "keyname for observer distance"}, 
  {ARG_FLAG,	"n",	"", "turns off tracking; target cl only determines time"}, 
  {ARG_FLAG,	"r",	"", "turns off apodization"}, 
  {ARG_FLAG,	"t",	"", "turns on time range checks"}, 
  {ARG_FLAG,	"v",	"", "runs in verbose mode"}, 
  {}
};

#include "selstuff.c"
#include "solephem.c"
#include "cartography.c"
#include "imginfo.c"

#define OUTLIER_RATIO	(6.0)
#define OUTLIER_BASE	(400.0)

		 /*  global declaration of missing to be initialized as NaN  */
float missing_val;

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
  float v, blank = 0.0/ 0.0;
  long long n;
  int cl[8];
  int col, cols, row, rows, i;
  int count = 0;
  if (accept <= 0.0) return -1;
  if (accept < 1.0) accept = 1.0 / accept;
  amx = accept;
  amn = 1.0 / accept;
/*
  MDI_image_reconstruct (img);
*/
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
/*
if (count < 10)
printf ("[%d, %d] : z[%d] = %f out of range %f (%.1f)\n",
col, row, n, z[n], s, w);
*/
	  z[n] = blank;
	  count++;
	} else if ((s > baseval) && (v < amn * s)) {
/*
if (count < 10)
printf ("[%d, %d] : z[%d] = %f out of range %f (%.1f)\n",
col, row, n, z[n], s, w);
*/
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
			   /*  local estimator of Carrington rotation number
			       (only used for default when cr not specified)  */
#define CARRINGTON_EPOCH    ("1853.11.09_12:00")
#define CARR_ROT_SYNODIC    (27.275311 * 86400.0)
			   /*  Synodic period from sidereal period of 25.38d
					      and tropical year of 365.2422d  */
double carrington_rots (TIME obs_time) {
  TIME upd;
  double ephem[EPHEM_SIZE];
  double car_rot, clong, clest;
  double r, lat, lon, vr, vn, vw;
  int carr_ct;

  car_rot = 1.0 + (obs_time - sscan_time (CARRINGTON_EPOCH)) / CARR_ROT_SYNODIC;
  carr_ct = car_rot;
  clest = car_rot - carr_ct;
  calc_sun_ephemeris (obs_time, ephem, 0.0, 0.0);
  lon = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
  clong = 1.0 - lon / 360.0;
  if ((clong - clest) > 0.5) carr_ct--;
  if ((clest - clong) > 0.5) carr_ct++;
  return (carr_ct + clong);
}
/*
 *  Function to determine the closest acceptable record, if any,
 *    to the target nrt, within 0.5 of rstp
 *  Acceptability is based on the presence of of a T_OBS value within
 *   the range (t0, t1), a quality keyword value matching bits in qmask
 *   and presence of the T_REC_index value in a rejection list
 */
int good_record (DRMS_RecordSet_t *ds, int nrt, int rstp, float max_reach,
    TIME t0, TIME t1, int time_check, char *qual_key, unsigned int qmask,
    int *rejects, int *reject_list, int *shift) {
  DRMS_Record_t *rec;
  TIME t, tobs;
  unsigned int quality;
  int idrec, match;
  int n, nr, nrr, nrmn, nrmx, offset, status;

  static int qcheck = 1;
  int found = 1;
  static unsigned char *ok = NULL;

  rec = ds->records[nrt];
  if (time_check) {
    tobs = drms_getkey_time (rec, "T_OBS", &status);
    t = tobs - MISSION_EPOCH;
    if (t < t0 || t > t1 || status) found = 0;
  }
  if (qcheck) {
	/*  check for data missing or otherwise unacceptable bits in quality  */
    quality = drms_getkey_int (rec, qual_key, &status);
    if (status) qcheck = 0;
    else if (quality & qmask) found = 0;
  }
  if (*rejects) {
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
    tobs = drms_getkey_time (rec, "T_OBS", &status);
    t = tobs - MISSION_EPOCH;
    if (t < t0 || t > t1 || status) {
      ok[nr] = 0;
      continue;
    }
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
/******************************************************************************/
			/*  module body begins here  */
/******************************************************************************/
int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ds;
  DRMS_Record_t *rec;
  DRMS_Segment_t *seg;
  DRMS_Array_t *mgram;
  struct maploc {
    double lat;
    double lon;
    double wt;
  } *mloc;
  FILE *out;
  TIME t, t0, t1, tcm;
  double *sum, *suma, *sb, *sbt, *sf, *sft, *s1, *st, *st2, *det;
  double img_radius, ellipse_e, ellipse_pa;
  double x, y, x0, y0, xstp, ystp, mscale;
  double xc, yc, xscale, yscale, imgscale, latc, lonc, peff;
  double lat, lon;
  double twt, tfac, wtfac, img_weight;
  double f0, f1;
  float *clat, *clon;
  float b;
  int *reject_list = NULL;
  int nrec, recct, use, valct;
  int mcol, mcols, mrow, mrows, map, maps, n, ntot, s, size;
  int plate_cols, plate_rows, plate_width;
  int x_invrt, y_invrt;
  int need_ephem;
  int rank, off, count;
  int badqual, rejects, shifted, shiftct;
  char *dpc_str;
  char module_ident[64], time_str[64], key[64];
						    /*  constants & counters  */
  double raddeg = M_PI / 180.0;
  double degrad = 1.0 / raddeg;
  int projection = LAMBERT;
  int status = 0;
  missing_val = 0.0 / 0.0;
							      /*  parameters  */
  char *inset = strdup (params_get_str (params, "ds"));
  char *filename = strdup (params_get_str (params, "file"));
  int cr = params_get_int (params, "cr");
  float cl = params_get_float (params, "cl");
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  int latct = params_get_int (params, "lat_nvals");
  int lonct = params_get_int (params, "lon_nvals");
  float map_scale = params_get_float (params, "scale");
  float map_size = params_get_float (params, "extent");
  double apode_inner = params_get_double (params, "mask_in");
  double apode_outer = params_get_double (params, "mask_ex");
  double noise_floor = params_get_double (params, "floor");
  double interval = params_get_double (params, "interval");
  int rec_step = params_get_int (params, "rec_step");
  float max_reach = params_get_float (params, "max_reach");
  char *rejectfile = strdup (params_get_str (params, "reject"));
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *crot_key = strdup (params_get_str (params, "crot_key"));
  char *rsun_key = strdup (params_get_str (params, "rsun_key"));
  char *apsd_key = strdup (params_get_str (params, "apsd_key"));
  char *dsun_key = strdup (params_get_str (params, "dsun_key"));

  int no_track = params_isflagset (params, "n");
  int no_apode = params_isflagset (params, "r");
  int time_check = params_isflagset (params, "t");
  int verbose = params_isflagset (params, "v");
  int filereq = strcasecmp (filename, "Not Specified");

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) fprintf (stderr, "%s: JSOC version %s\n", module_ident, jsoc_version);

  if (cr < 0) {
		       /*  just approximate to previous carrington longitude  */
    double rsun, lat, lon, vr, vn, vw;
    tcm = CURRENT_SYSTEM_TIME;
    earth_ephemeris (tcm, &rsun, &lat, &lon, &vr, &vn, &vw);
    cr = carrington_rots (tcm);
    cr--;
    if (isnan (cl)) cl = lon;
  } else if (isnan (cl)) cl = 0.0;

  if (key_params_from_dspec (inset)) {
fprintf (stderr, "Error: data set specification not supported\n");
fprintf (stderr, "       use data series and params cr, cl, and length\n");
return 0;
  } else {
				/*  only the input data series is named,
				    get record specifications from arguments  */
    char ttarget[64];
    sprintf (ttarget, "%d:%06.2f", cr, cl);
/*  hack to approximate longitude length from length in minutes  */
/*  synodic rotation rate varies from ~ 13.17 - 13.23 deg/day  */
    if (!(ds = select_dataset_from_time_interval (inset, ttarget, interval)))
/*
	interval * 13.2 / 1440.0)))
*/
      return 1;
    if ((recct = ds->n) < 2) {
      printf ("<2 records in selected input set\n");
      drms_close_records (ds, DRMS_FREE_RECORD);
      return 1;
    }
  }
		 /*  support special hack of reading of rejection list file  */
  rejects = 0;
  if (strcmp (rejectfile, "Not Specified")) {
    FILE *rejectfp = fopen (rejectfile, "r");
    if (rejectfp) rejects = read_reject_list (rejectfp, &reject_list);
    else fprintf (stderr,
	"Warning: could not open rejection list %s; ignored\n", rejectfile);
  }

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
  for (map = latct; map < maps; map++) clat[map] = clat[map-1];
  for (map = lonct; map < maps; map++) clon[map] = clon[map-1];

  if (filereq) {
    if (!(out = fopen (filename, "w"))) {
      fprintf (stderr, "Error - unable to open file %s for table output\n", filename);
      return 1;
    }
    if (verbose) fprintf (out, "CR %d:%03.0f\n", cr, cl);
  } else verbose = 0;
/*
 *  determine mission time corresponding to meridian crossing of given
 *    longitude for given Carrington rotation and limits of tracking interval
 */
  tcm = earth_meridian_crossing (cl, cr);
  tcm -= MISSION_EPOCH;
/*  hack to approximate interval length in min from longitude interval  */
/*  synodic rotation rate varies from ~ 13.17 - 13.23 deg/day  */
  interval *= 1440 / 13.2;
  t0 = tcm - interval * 30.0;
  t1 = tcm + interval * 30.0;
  tfac = 2.0 / (t1 - t0);
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
  ntot = 0;
  for (mrow=0, y=y0; mrow<mrows; mrow++, y+=ystp) {
    for (mcol=0, x=x0; mcol<mcols; mcol++, x+=xstp, ntot++) ;
  }
  ntot *= maps;
  mloc = (struct maploc *)malloc (ntot * sizeof (struct maploc));
  sum = (double *)malloc (maps * sizeof (double));
  suma = (double *)malloc (maps * sizeof (double));
  det = (double *)malloc (maps * sizeof (double));
  s1 = (double *)calloc (maps, sizeof (double));
  st = (double *)calloc (maps, sizeof (double));
  st2 = (double *)calloc (maps, sizeof (double));
  sb = (double *)calloc (maps, sizeof (double));
  sbt = (double *)calloc (maps, sizeof (double));
  sf = (double *)calloc (maps, sizeof (double));
  sft = (double *)calloc (maps, sizeof (double));

  for (n = 0, map = 0; map < maps; map++) {
    sum[map] = 0.0;
    for (mrow = 0, y = y0; mrow < mrows; mrow++, y += ystp) {
      for (mcol = 0, x = x0; mcol < mcols; mcol++, x += xstp, n++) {
	off = plane2sphere (x, y, clat[map], clon[map], &lat, &lon, projection);
	mloc[n].lat = lat;
	mloc[n].lon = lon;
	mloc[n].wt = (off) ? 0.0 : (no_apode) ? 1.0 :
	    apodization (x * mscale, y * mscale, apode_inner, apode_outer);
	sum[map] += mloc[n].wt;
      }
    }
  }
  for (map = 0, n = 0; map < maps; map++) {
    if (sum[map] != 0.0) {
      wtfac = 1.0 / sum[map];
      for (s = 0; s < size; s++, n++)
	mloc[n].wt *= wtfac;
    } else n += size;
  }
  valct = 0;
					      /*  loop through input records  */
  badqual = shiftct = shifted = 0;
  for (nrec = 0; nrec < recct; nrec += rec_step) {
    float bmax = -HUGE_VAL, bmin = HUGE_VAL;
    double lonmax, lonmin, latmax, latmin;
		/*  call function to determine the closest acceptable record,
			  if any, to the target nrec within +/- 0.5 rec_step  */
    
    use = good_record (ds, nrec, rec_step, max_reach, t0, t1, time_check,
	qual_key, qmask, &rejects, reject_list, &shifted);
    if (use < 0) {
/*
fprintf (stderr, "record %d skipped after trying all within %d\n", nrec,
shifted);
*/
      if (shifted) badqual++;
      continue;
    }
    if (shifted) shiftct++;
    rec = ds->records[use];
    t = drms_getkey_time (rec, "T_OBS", &status) - MISSION_EPOCH;
    seg = drms_segment_lookupnum (rec, 0);
    mgram = drms_segment_read (seg, DRMS_TYPE_FLOAT, &status);
    if (!mgram || status) {
/*
      sprint_time (time_str, tobs, "TAI", 0);
      fprintf (stderr, "Error: segment read failed for record %s\n", time_str);
*/
      if (mgram) drms_free_array (mgram);
      continue;
    }
    if ((rank = drms_array_naxis (mgram)) != 2) {
      fprintf (stderr, "improper format for record %d: rank = %d\n", nrec,
	  rank);
      drms_free_array (mgram);
      continue;
    }
    plate_cols = drms_array_nth_axis (mgram, 0);
    plate_rows = drms_array_nth_axis (mgram, 1);
    twt = tfac * (t - tcm);
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
	fprintf (stderr,
	    "Error: one or more essential keywords or values missing; skipped\n");
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
      fprintf (stderr,
          "Error: one or more essential keywords or values missing; skipped\n");
      fprintf (stderr, "solar_image_info() returned %08x\n", status);
      continue;
    }
    lonc = drms_getkey_double (rec, clon_key, &status) * raddeg;
    latc = drms_getkey_double (rec, clat_key, &status) * raddeg;
    if (no_track) lonc = cl * raddeg;

    imgscale = (xscale > yscale) ? xscale : yscale;
    if (verbose && filereq) {
      sprint_time (time_str, t + MISSION_EPOCH, "UT", -1);
      dpc_str = drms_getkey_string (rec, "DPC", &status);
      if (!status)
	fprintf (out, "%.16s : %08lx\n", time_str, strtol (dpc_str, NULL, 16));
      else fprintf (out, "%.16s\n", time_str);
    }
						 /*  Truncate "noise" values  */
    status = truncate_noise (mgram, noise_floor);

    if (verbose && filereq) fprintf (out,
	"  zeroed %d (of %d) values within limits +/- %.1f\n",
	status, count, noise_floor);

    count -= status;
						      /*  Eliminate outliers  */
    status = eliminate_outliers (mgram, OUTLIER_RATIO, OUTLIER_BASE);
    if (status && verbose && filereq)
      fprintf (out, "  removed %d (of %d) outlier(s)\n", status, count);
      /*  Correct values for geometric effects in observed coordinate system  */
/*
    status = correct_foreshort (mgram);
*/
    status = 0;
    plate_width = (plate_cols > plate_rows) ? plate_cols : plate_rows;
    xc -= 0.5 * (plate_cols - 1);
    yc -= 0.5 * (plate_rows - 1);
    xc *= 2.0 / plate_width;
    yc *= 2.0 / plate_width;
    img_radius *= 2.0 / plate_width;
				      /*  map field values for each location  */
    for (map = 0, n = 0; map < maps; map++) {
      sum[map] = suma[map] = 0.0;
      ntot = 0;
      img_weight = 1.0;
      for (s = 0; s < size; s++, n++) {
	if (mloc[n].wt == 0.0) continue;
	if (sphere2img (mloc[n].lat, mloc[n].lon, latc, lonc, &x, &y, xc, yc,
	      img_radius, peff, 0.0, 0.0, 0, 0)) {
	  img_weight -= mloc[n].wt;
	  continue;
	}
	b = array_imaginterp (mgram, x, y, 0);
	if (isnan (b)) {
	  img_weight -= mloc[n].wt;
	  continue;
	}
	if (b > bmax) {
	  bmax = b;
	  lonmax = mloc[n].lon;
	  latmax = mloc[n].lat;
	}
	if (b < bmin) {
	  bmin = b;
	  lonmin = mloc[n].lon;
	  latmin = mloc[n].lat;
	}
	sum[map] += b * mloc[n].wt;
	suma[map] += fabs (b) * mloc[n].wt;
	ntot++;
      }
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

    if (verbose && filereq) {
      fprintf (out, "  min = %7.1f at [%6.2f, %6.2f]  ", bmin,
	  lonmin/raddeg, latmin/raddeg);
      fprintf (out, "  max = %7.1f at [%6.2f, %6.2f]\n", bmax,
	  lonmax/raddeg, latmax/raddeg);
    }
    drms_free_array (mgram);
    valct++;
  }
  drms_close_records (ds, DRMS_FREE_RECORD);
  if (valct < 2) {
    fprintf (stderr, "Error: <2 usable magnetograms in interval\n");
    if (filereq) fclose (out);
    return 1;
  }
  if (verbose && filereq) fprintf (out, "\n");

  for (map = 0; map < maps; map++)
    det[map] = 1.0 / (s1[map] * st2[map] - st[map] * st[map]);
  if (filereq) fprintf (out, " Lon   Lat       |Bz|   d|Bz|/dt    Bz     dBz/dt [Gauss(/sec)]:\n");

  for (map = 0; map < maps; map++) {
    f0 = det[map] * (st2[map] * sf[map] - st[map] * sft[map]);
    printf ("%8.3f", f0);
    if (filereq) {
      double fs0 = det[map] * (st2[map] * sb[map] - st[map] * sbt[map]);
      double fs1 = det[map] * tfac * (s1[map] * sbt[map] - st[map] * sb[map]);
      f1 = det[map] * tfac * (s1[map] * sft[map] - st[map] * sf[map]);
      fprintf (out, "%05.1f %+05.1f : %8.3f %8.3f %8.3f %8.3f\n",
	  clon[map] / raddeg, clat[map] / raddeg, f0, 1.0e6 * f1, fs0,
	  1.0e6 * fs1);
    }
  }
  printf ("\n");
  if (filereq) fclose (out);
  if (verbose && badqual) fprintf (stderr,
      "    %d of %d records rejected for quality matching %08x\n",
      badqual, valct, qmask);
  if (verbose && shiftct) fprintf (stderr,
      "%d records shifted from  targets for quality issues\n", shiftct);

  return status;
}

/*
 *  Revision history (all mods by Rick Bogart unless otherwise noted)
 *
 *  04.04.27	created this file, based on SOI module maginteg
 *  v 0.0	compiles and builds as standalone program
 *  v 0.1	builds as DRMS module
 *  v 0.2	runs with minimal functionality
 *  v 0.3	full argument list and parsing, fixed maginteg bug (cr unspec)!
 *  v 0.4	opens records from DRMS
 *  v 0.5	does everything but the actual mapping (interpolation)
 *  v 0.6	full functionality of original maginteg, though output values
 *	differ drastically; clearly buggy
 *  v 0.7	fixed various bugs relating to keyword interpretation and
 *	uninitialized variables; results agree approximately, though not
 *	exactly, with those from maginteg run on same data with comparable
 *	parameters
 *  v 0.7 frozen 2010.04.29
 *  v 0.8	Added stdout of MAI (|B|) values for selected regions as a
 *	single line, appropriate for parameter passing to mtrack
 *		Modified file output to only include per record info when
 *	running verbose, and not to send any extra info to stdout; made
 *	output to file optional
 *		Take target list of latitudes and longitudes a la mtrack
 *	rather than a fixed tabular range
 *		Fixed various bugs relating to selection of the integration
 *	interval in terms of longitude rotation; optionally print out both
 *	signed and unsigned means and time derivs (and verbose output) to named
 *	file
 *		Removed default in series; removed check on image size
 *	(seg faults)
 *		Added dsun_key argument for observer distance keyword, and
 *	mechanism for extracting semi-diameter from observer distance if missing;
 *	other minimal fixes to get it to work with HMI test LOS magnetograms
 *	(without verifying that it still works with MDI data!)
 *		Modified defaults to be appropriate for HMI 15 deg tiles
 *		Added rec_step argument with default corresponding to 48 min
 *	for HMI; provided local version of Carrington rots estimator
 *  v 0.8 frozen 2010.08.20
 *  v 0.9	Modified for new call to solar_image_info
 *  v 0.9 frozen 2012.04.24
 *  v 1.0	Make sure 0 status is returned on successful completion;
 *		Clean up unused commented code; more complete error checking,
 *	including status return when there are no accepted magnetograms
 *  v 1.0 frozen 2016.03.10
 *  v 1.1	Removed default values for several arguments: interval,
 *	rec_step, scale, extent
 *		Added argument max_reach, with default value of 0.5 to control
 *	how far from target times to search for acceptable images
 *		Added ability to reject input magnetograms on basis of a
 *	a quality mask, and to search for acceptable "nearby" records (within
 *	max_reach * rec_step (which no longer defaults to 64)
 *		Made time range check for acceptability optional, set by
 *	default, but with intent that it be option requiring flag, preserved
 *	only for backward compatability
 *  v 1.1 frozen 2016.11.09
 */
