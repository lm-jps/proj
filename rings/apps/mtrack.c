/*
 *  mtrack.c						~rick/src/mtrack
 *
 *  Responsible:  Rick Bogart				RBogart@spd.aas.org
 *
 *  Generate multiple mapped tracked data cubes at different locations from
 *    a common time sequence of solar images
 *
 *  (Main module body begins around line 560)
 *
 *  Parameters: (type   default         description)
 *	in	DataSet none            Input dataset
 *				A series of images of all or part of the
 *				solar disc in a "plate" coordinate system
 *				(helioprojective geometry). It may be
 *				specified as either a data set query or
 *				just a data series name, in which case
 *				additional selection parameters are required
 *				The dataset is assumed to include in its
 *				ancillary data the parameters required for
 *				heliographic mapping, namely the observation
 *				time and heliographic location.
 *      segment string	-		Name of the input data segment (only
 *				required if the input series has multiple
 *				segments)
 *      out     DataSer none            Output data series name
 *				The output series must have prime keys that
 *				can distinguish individual output records by
 *				LatHG, LonHG, and either MidTime or CarrTime
 *				(= CarrRot:CMLon)
 *      bckgn	string	-		If specified, pointer to image to be
 *				be pre-subtracted from each input image. Can
 *				be either a data record segment specifier or
 *				a file name
 *	reject	string	-		If specified, name of a file with a
 *				list of input images to be rejected regardless
 *				of quality mask values
 *	qmask	int	0x00000000	Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	tmid	string	-		Midpoint of target tracking interval;
 *				ignored if input specified as data set; defined
 *				as string because it can be specified in either
 *				regular time format or as CR:CL
 *	length	int	-		Length of tracking interval, in units
 *				of data cadence; ignored if input specified as
 *				data set
 *	tstart	time	-		Start of target tracking interval;
 *				ignored if input specified as data set, or if
 *				tmid and length are specified
 *	tstop	time	-		End of target tracking interval;
 *				ignored if input specified as data set, or if
 *				tmid and length are specified
 *	lat	float*	[0.0]		Target heliographic latitude(s)
 *	lon	float*	[0.0]		Target heliographic longitude(s)
 *      map     enum  "Postel"               Mapping option:
 *				recognized values are "carree", "Cassini",
 *				"Mercator", "cyleqa", "sineqa", "gnomonic",
 *				"Postel", "stereographic", "orthographic",
 *				and "Lambert" (and possibly others).
 *      interp  enum  "cubiconv"           Interpolation option:
 *				recognized values are "cubiconv" (and possibly
 *				others)
 *      scale   float   0.0             Scale of map (heliographic degrees /
 *				pixel) at location appropriate for mapping
 *				option; a 0 value implies autoscaling to best
 *				scale of image.  Typical values would be
 *				about 0.057 * resolution [arcsec / pixel]
 *				for solar radii of 1000 arcsec, or about
 *				0.12 for MDI or GONG or MtWilson
 *				0.03 for HMI or AIA
 *	cols	int	0		Columns in output maps; 0 -> rows
 *	rows	int	0		Rows in output maps; 0 -> cols
 *      map_pa  float   0.0             The angle between heliographic north
 *				and "up" on the output map (in the direction
 *				of increasing rows) [deg[, in the sense that a
 *				positive position angle represents a clockwise
 *				displacement of the north axis.
 *      a0      float   -0.02893	Coefficients in sin^2 (latitude)
 *      a2      float   -0.3441		expansion of rotation rate minus
 *      a4      float   -0.5037		Carrington rotation (urad/sec)
 *				ignored if -c or -n flag
 *         Some representative values to use:
 *                 a0       a2      a4
 *		-0.02893 -0.3441 -0.5037	Snodgrass 82/84 (default)       
 *		-0.02713 -0.3014 -0.5263	Snodgrass 67/84
 *		-0.02533 -0.3079 -0.5278	used by Haber & Hill
 *		-0.02833 -0.4100 -0.4189	Ulrich et al. SP 117, 291
 *		-0.0853	 -0.351  -0.443		Howard & Harvey
 *		 0.03749 -0.5596  0		Newton & Nunn
 *		 0.0137  -0.339  -0.485		Snodgrass & Ulrich ApJ 351, 309
 *		 0.1067  -0.484  -0.361		Snodgrass & Ulrich " (66-87)
 *      merid_v float   0.0		Meridional rate for tracking
 *				urad/sec, +ve northward
 *	bscale	float	0.0		Value scaling parameter for output
 *	bzero	Float	0.0		value offset parameter for output
 *	trec_key string	T_REC		Keyname of time type prime keyword for
 *				input data series
 *	tobs_key string	T_OBS		Keyname of time type keyword describing
 *				observation time (midpoint) of each input image
 *	tstp_key string	Cadence		Keyname of float type keyword describing
 *				observing cadence of input data
 *	qual_key string	Quality		Keyname of int type keyword describing
 *				record quality as bitmask
 *	clon_key string	CRLN_OBS	Keyname of float type record describing
 *				centre Carrington longitude of each input image
 *	clat_key string	CRLT_OBS	Keyname of float type record describing
 *				centre Carrington latitude of each input image
 *	rsun_key string	R_SUN		Keyname of float type record describing
 *				apparent solar semidiameter of each image
 *      mai	 float*  NaN	MAI (Magnetic Activity Index) value to
 *				set as keyword
 *
 *  Flags
 *	-c	track at Carrington rate (a0 = a2 = a4 = merid_v = 0)
 *	-n	turn off tracking, just correct for anomalous observer motion
 *		(a0 = -2.8653291, a2 = a4 = merid_v = 0
 *	-o	remove line-of-sight component of observer velocity
 *	-r	remove line-of-sight component of solar rotation
 *	-v	run verbose
 *	-x	experimental run (do not insert output records)
 *	-G	use GONG keywords for input
 *	-M	use MDI keywords for input and correct for MDI distortion
 *
 *  Notes:
 *    This module is a DRMS-based version of fastrack, which it is intended
 *	to replace.
 *
 *  Bugs:
 *    There was a bug in NetDRMS 2.0b drms_segment_writeslice that can cause
 *	the code to compress the wrong tile size, or even fail with a floating-
 *	point exception, when the output segment has a vardim protocol. It
 *	is not known if this bug is still present in 2.0
 *    The code does not check that the output data segments match those of
 *	the output series when the segment protocol is variable
 *    Checks for validity of tstart and tstop parameters are unreliable, due
 *	to bugs in DRMS treatment of scans of invalid strings for times
 *	(ticket #177)
 *    The Log segment, if written, is only written to the first record, but
 *	is overwritten for each output record!
 *    Should free log, map
 *    Will not accept an input dataset specification of @*
 *    The image foreshortening corrections are appropriate for 1 AU, independent
 *	of DSUN_OBS
 *    The input of the ellipse position angle has not been verified, but then
 *	neither has the correction for ellipticity of the image altogether
 *    If there are multiple data segments of rank 3 in the output series, a
 *	warning is generated and the first one is used
 *    The protocol of the Log segment, if present in the output series, is not
 *	verified to be generic
 *    There is no verification that numerous essential key data are actually
 *	in the input dataset, in particular trec_key, tobs_key, and
 *	crot_key. (clon_key and clat_key are checked)
 *    Uses private local FITS code rather than cfitsio for reading in
 *	background image
 *    When a target time and length are given, the actual number of records
 *	may differ by one from the expected number; this only occurs if
 *	the target time differs by a small but non-zero amount (< 0.1 sec)
 *	from either the data record time or the midway point between data
 *	record times. This happens for example, if the length is odd and
 *	the target time differs from an actual observation time by more than
 *	0.1 usec but less than 0.05 sec, with a 1-minute cadence
 *    The input data are unconditionally read in as floats (single-precision),
 *	and the output data written as scaled shorts.
 *    There is evidently no WCS conventional name for the Cassini-Soldner
 *	(transverse plate carree) projection; CAS is arbitrarily used; the
 *	alternative would be to interchange HGLT and HGLN, but that would
 *	necessitate a change in the position angle
 *
 *  Future Updates
 *    0.9 incorporate writing of output by slices to avoid memory limits;
 *	reorganize code, adding functions and simplifying main module;
 *	fix up keywords; add stubs for initial input oversampling
 *	(including array limits)
 *    Add target cadence arg, modify length to be in units of target cadence
 *    Add option for tracking from interpolated set of target locations
 *    Add ability to track from remapped images, not just helioprojective
 *    Interpolate in time to target cadence
 *    Create new data series as needed
 *    Check for MDI source input data and consistency of flag
 *
 *  Revision history is at end of file
 */

#include <jsoc_main.h>

						      /*  module identifier  */
char *module_name = "mtrack";
char *module_desc = "track multiple regions from solar image sequences";
char *version_id = "0.9";

#define CARR_RATE       (2.86532908457)
#define RSUNM		(6.96e8)
#define INTERP_NEAREST_NEIGHBOR	(1)

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"in", "", "input data series or dataset"}, 
  {ARG_STRING,	"segment", "Not Specified",
      "input data series segment; ignored if series only has one segment"}, 
  {ARG_STRING,	"out",	"", "output data series"}, 
  {ARG_STRING,	"bckgn", "Not Specified",
      "background image to be subtracted from input data"}, 
  {ARG_STRING,	"reject", "Not Specified", "file containing rejection list"}, 
  {ARG_INT,	"qmask", "0x00000000", "quality bit mask for image rejection"},
  {ARG_STRING,	"tmid", "Not Specified", "midpoint of tracking interval"}, 
  {ARG_INT,	"length", "0",
      "target length of tracking interval [input cadence]"}, 
				      /*  necessitated by bug (ticket #177)  */
  {ARG_TIME,	"tstart", "JD_0", "start of tracking interval"}, 
  {ARG_TIME,	"tstop", "JD_0", "end of tracking interval"}, 
  {ARG_FLOAT,	"tstep", "Not specified", "temporal cadence for output"},
  {ARG_FLOATS,	"lat", "[0.0]",
      "heliographic latitude(s) of tracking center(s) [deg]"},
  {ARG_FLOATS,	"lon", "[0.0]",
      "heliographic longitude(s) of tracking center(s) [deg]"},
  {ARG_NUME, "map", "Postel", "map projection",
      "carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
  {ARG_NUME, "interp", "cubiconv", "interpolation option",
      "cubiconv, nearest"},
  {ARG_FLOAT,	"scale", "Not specified", "map scale at center [deg/pxl]"},
  {ARG_INT,     "cols", "0", "columns in output map(s)"},
  {ARG_INT,     "rows", "0", "rows in output map(s)"},
  {ARG_FLOAT,	"map_pa", "0.0", "position angle of north on output map"},
  {ARG_FLOAT,	"a0", "-0.02893",
      "equatorial rotation - Carrington rate [urad/sec]"},
  {ARG_FLOAT,	"a2", "-0.3441", "solar rotation parameter A2 [urad/sec]"},
  {ARG_FLOAT,	"a4", "-0.5037", "solar rotation parameter A4 [urad/sec]"},
  {ARG_FLOAT,	"merid_v", "0.0",
      "solar meridional velocity [urad/sec]; 0.0014368 * rate in m/s"},
  {ARG_FLOAT,	"bscale", "0.0", "output scaling factor"},
  {ARG_FLOAT,	"bzero", "0.0", "output offset"},
  {ARG_STRING,	"trec_key", "T_REC", "keyname of (slotted) prime key"}, 
  {ARG_STRING,	"tobs_key", "T_OBS", "keyname for image observation time"}, 
  {ARG_STRING,	"tstp_key", "CADENCE",  "keyname for image observation time"}, 
  {ARG_STRING,	"qual_key", "Quality",  "keyname for 32-bit image quality field"}, 
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"crot_key", "CAR_ROT", "keyname for image Carrington rotation"}, 
  {ARG_STRING,	"rsun_key", "R_SUN", "keyname for apparent solar semi-diameter"}, 
  {ARG_FLOATS,	"mai", "NaN", "(MDI) Magnetic Activity Indices"},
  {ARG_FLAG,	"c",	"", "track at Carrington rate"}, 
  {ARG_FLAG,	"n",	"", "turn off tracking"}, 
  {ARG_FLAG,	"o",	"",
      "remove line-of-sight component of observer velocity"}, 
  {ARG_FLAG,	"r",	"",
      "remove line-of-sight component of solar rotation"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {ARG_FLAG,	"x",	"", "experimental mode (do not save output)"}, 
  {ARG_FLAG,	"G",	"", "use GONG keywords for input"}, 
  {ARG_FLAG,	"M",	"",
      "use MDI keywords for input and correct for MDI distortion"}, 
  {}
};

#include "keystuff.c"
#include "earth_ephem.c"
#include "soho_ephem.c"
#include "fitstuff.c"
#include "cartography.c"
#include "imginfo.c"
#include "mdistuff.c"

char *prime_keylist[] = {"CarrRot", "CMLon", "LonHG", "LatHG", "Duration",
    "Width", "Height", "Source", "ZonalTrk", "MeridTrk", "MapProj", "MapScale"};

static char *create_prime_key (DRMS_Record_t *rec) {
  int pct = sizeof (prime_keylist) / sizeof (char *);

  return create_primekey_from_keylist (rec, prime_keylist, pct);
}

/*
 *  Functions to support reading from a specified rejection list of the
 *    appropriate format
 */

static int fgetline (FILE *in, char *line, int max) {
  if (fgets (line, max, in) == NULL) return 0;
  else return (strlen (line));
}

static int read_reject_list (FILE *file, int **list) {
  int ds, sn, rec, last_rec;
  int allocd = 1024, ct = 0, gap = 0;
  char line[1024], t_str[64], estr[16];

  *list = (int *)malloc (allocd * sizeof (int));
  while (fgetline (file, line, 1024)) {
    if (strlen (line) == 1) continue;
    if (line[0] == '#') continue;
    if (sscanf (line, "%d %d %d %s", &ds, &sn, &rec, t_str) != 4) {
     sscanf (line, "%s", estr);
      if (strcmp (estr, "...")) continue;
      gap = 1;
      last_rec = rec;
      continue;
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

static int key_params_from_dspec (const char *dspec) {
/*
 *  Establish whether target times are determined from dataset specifier
 *  assume that if a bracket is found in the dataset specifier it is a
 *  set of well-specified records containing a properly ordered input dataset;
 *  otherwise, a dataseries to be queried
 */
  int n, nt = strlen (dspec);

  for (n = 0; n < nt; n++) if (dspec[n] == '[') return 1;
  return 0;
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
    return FP_NAN;

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
    return FP_NAN;
  col = x + 0.5;
  row = y + 0.5;
  return f[col + row * cols];
}

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

static int fsphere2img (double sin_lat, double cos_lat, double lon,
	double latc, double lonc, double *x, double *y,
	double xcenter, double ycenter, double rsun,
	double cospa, double sinpa, double ecc, double chi, int ellipse,
	int xinvrt, int yinvrt) {
/*
 *  Optimized version of sphere2img() for case when lat and pa do not
 *    change frequently
 */
  static double sin_asd = 0.004660, cos_asd = 0.99998914;
                                                   /*  appropriate to 1 AU  */
  static double last_latc = 0.0, cos_latc = 1.0, sin_latc = 0.0;
  double r, cos_cang, xr, yr;
  double cos_lat_lon;
  double squash, cchi, schi, c2chi, s2chi, xp, yp;
  int hemisphere;

  if (latc != last_latc) {
      sin_latc = sin (latc);
      cos_latc = cos (latc);
      last_latc = latc;
  }
  cos_lat_lon = cos_lat * cos (lon - lonc);

  cos_cang   = sin_lat * sin_latc + cos_latc * cos_lat_lon;
  hemisphere = (cos_cang < 0.0) ? 1 : 0;
  r = rsun * cos_asd / (1.0 - cos_cang * sin_asd);
  xr = r * cos_lat * sin (lon - lonc);
  yr = r * (sin_lat * cos_latc - sin_latc * cos_lat_lon);
                                        /*  Change sign for inverted images  */
  if (xinvrt) xr *= -1.0;
  if (yinvrt) yr *= -1.0;
                          /*  Correction for ellipsoidal squashing of image  */
  if (ellipse) {
    squash = sqrt (1.0 - ecc * ecc);
    cchi = cos (chi);
    schi = sin (chi);
    s2chi = schi * schi;
    c2chi = 1.0 - s2chi;
    xp = xr * (s2chi + squash * c2chi) - yr * (1.0 - squash) * schi * cchi;
    yp = yr * (c2chi + squash * s2chi) - xr * (1.0 - squash) * schi * cchi;
    xr = xp;
    yr = yp;
  }
  *x = xr * cospa - yr * sinpa;
  *y = xr * sinpa + yr * cospa;

  *y += ycenter;
  *x += xcenter;
  return (hemisphere);
}

float missing_val = FP_NAN;

static void perform_mappings (DRMS_Array_t *img, float *maps, double *delta_rot,
    int mapct, double *maplat, double *maplon,
    double *map_coslat, double *map_sinlat, int pixct,
    double delta_time, unsigned char *offsun,
    double latc, double lonc, double xc, double yc,
    double a0, double a2, double a4, double merid_v, double radius, double pa,
    double ellipse_e, double ellipse_pa, int x_invrt, int y_invrt,
    int interpolator, int MDI_correct_distort) {
/*
 *  Perform the mappings from the target heliographic coordinate sets
 *    appropriate to each output cube into the image coordinates (as
 *    corrected) for spatial interpolation of the data values
 */
		/*  STUB function: checked so far: delta_time, delta_lat,
	latc, lonc, xc, yc, x_invrt, y_invrt, ellipse_e, delta_lon, radius, pa  */
  static double sin_asd = 0.004660, cos_asd = 0.99998914;
                                                   /*  appropriate to 1 AU  */
  double lat, lon, cos_lat, sin_lat;
  double xx, yy;
  float interpval;
  int m, n, mset;

  double delta_lat = merid_v * delta_time;
  int plate_cols = img->axis[0];
  int plate_rows = img->axis[1];
  double plate_width = (plate_cols > plate_rows) ? plate_cols : plate_rows;
  int no_merid_v = (merid_v == 0.0);
  double cos_pa = cos (pa);
  double sin_pa = sin (pa);
  int fit_ellipse = (ellipse_e > 0.0 && ellipse_e < 1.0);

  mset = 0;
  xc *= 2.0 / plate_width;
  yc *= 2.0 / plate_width;
  radius *= 2.0 / plate_width;

  if (no_merid_v) {
    if (x_invrt || y_invrt || fit_ellipse) {
      for (m = 0; m < mapct; m++) {
	double delta_lon = delta_rot[m] * delta_time;
	for (n = 0; n < pixct; n++) {
	  if (offsun[n + mset]) {
	    maps[n + mset] = missing_val;
            continue;
	  }
       /*  Calculate heliographic coordinates corresponding to map location  */
	  sin_lat = map_sinlat[n + mset];
	  cos_lat = map_coslat[n + mset];
	  lon = maplon[n + mset] + delta_lon;
     /*  Calculate plate coordinates corresponding to heliocentric location  */
	  if (fsphere2img (sin_lat, cos_lat, lon, latc, lonc, &xx, &yy, xc, yc,
	      radius, cos_pa, sin_pa, ellipse_e, ellipse_pa, fit_ellipse,
	      x_invrt, y_invrt)) {
	    maps[n + mset] = missing_val;
	    continue;
	  }
	  if (plate_cols > plate_rows) yy -= 1.0 - plate_rows / plate_width;
	  if (plate_rows > plate_cols) xx -= 1.0 - plate_cols / plate_width;
				  /*  Correction for MDI distortion and tip  */
       /*  should be replaced by call to MDI_correct_plateloc when verified  */
	  if (MDI_correct_distort) {
	    mtrack_MDI_image_tip (&xx, &yy, 1, 1);
	    mtrack_MDI_image_stretch (&xx, &yy, 1, 1);
	  }
	  interpval = array_imaginterp (img, xx, yy, interpolator);
	  maps[n + mset] = (isnan (interpval)) ? missing_val : interpval;
	}
	mset += pixct;
      }
    } else {
      double r, cos_cang, xr, yr, cos_lat_lon;
      double cos_latc = cos (latc);
      double sin_latc = sin (latc);
      for (m = 0; m < mapct; m++) {
	double delta_lon = delta_rot[m] * delta_time;
	for (n = 0; n < pixct; n++) {
	  if (offsun[n + mset]) {
	    maps[n + mset] = missing_val;
	    continue;
	  }
       /*  Calculate heliographic coordinates corresponding to map location  */
	  sin_lat = map_sinlat[n + mset];
	  cos_lat = map_coslat[n + mset];
	  lon = maplon[n + mset] + delta_lon;
	  cos_lat_lon = cos_lat * cos (lon - lonc);
	  cos_cang  = sin_lat * sin_latc + cos_latc * cos_lat_lon;
	  if (cos_cang < 0.0) {
	    maps[n + mset] = missing_val;
	    continue;
	  }
     /*  Calculate plate coordinates corresponding to heliocentric location  */
	  r = radius * cos_asd / (1.0 - cos_cang * sin_asd);
	  xr = r * cos_lat * sin (lon - lonc);
	  yr = r * (sin_lat * cos_latc - sin_latc * cos_lat_lon);
	  xx = xr * cos_pa - yr * sin_pa;
	  yy = xr * sin_pa + yr * cos_pa;
	  yy += yc;
	  xx += xc;
		  /*  should take tests outside loop, just modify xc and yc  */
	  if (plate_cols > plate_rows) yy -= 1.0 - plate_rows / plate_width;
	  if (plate_rows > plate_cols) xx -= 1.0 - plate_cols / plate_width;
				  /*  Correction for MDI distortion and tip  */
       /*  should be replaced by call to MDI_correct_plateloc when verified  */
	  if (MDI_correct_distort) {
	    mtrack_MDI_image_tip (&xx, &yy, 1, 1);
	    mtrack_MDI_image_stretch (&xx, &yy, 1, 1);
	  }
	  interpval = array_imaginterp (img, xx, yy, interpolator);
	  maps[n + mset] = (isnan (interpval)) ? missing_val : interpval;
        }
	mset += pixct;
      }
    }
    return;
  }
/*
 *  only executed if there is a meridional component to the tracking rate
 */
  for (m = 0; m < mapct; m++) {
    double delta_lon = delta_rot[m] * delta_time;
    for (n = 0; n < pixct; n++) {
      if (offsun[n + mset]) {
        maps[n + mset] = missing_val;
        continue;
      }
       /*  Calculate heliographic coordinates corresponding to map location  */
      lat = maplat[n + mset] + delta_lat;
      lon = maplon[n + mset] + delta_lon;
     /*  Calculate plate coordinates corresponding to heliocentric location  */
      if (sphere2img (lat, lon, latc, lonc, &xx, &yy, xc, yc, radius, pa,
          ellipse_e, ellipse_pa, x_invrt, y_invrt)) {
	maps[n + mset] = missing_val;
	continue;
      }
      if (plate_cols > plate_rows) yy -= 1.0 - plate_rows / plate_width;
      if (plate_rows > plate_cols) xx -= 1.0 - plate_cols / plate_width;
				  /*  Correction for MDI distortion and tip  */
      if (plate_cols > plate_rows) yy -= 1.0 - plate_rows / plate_width;
      if (plate_rows > plate_cols) xx -= 1.0 - plate_cols / plate_width;
				  /*  Correction for MDI distortion and tip  */
      if (MDI_correct_distort) {
	mtrack_MDI_image_tip (&xx, &yy, 1, 1);
	mtrack_MDI_image_stretch (&xx, &yy, 1, 1);
      }
      interpval = array_imaginterp (img, xx, yy, interpolator);
      maps[n + mset] = (isnan (interpval)) ? missing_val : interpval;
    }
    mset += pixct;
  }
}

static int verify_keys (DRMS_Record_t *rec, const char *clon,
    const char *clat, double *keyscale) {
  DRMS_Keyword_t *keywd;
  
  keywd = drms_keyword_lookup (rec, clon, 1);
  if (!keywd) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for Carrington longitude of observer not found\n",
	clon);
    fprintf (stderr, "       Must supply an appropriate value for clon_key\n");
    fprintf (stderr, "       (Carrington longitude of disc center in deg)\n");
    return -1;
  }
  if (keywd->info->type == DRMS_TYPE_TIME ||
      keywd->info->type == DRMS_TYPE_STRING ||
      keywd->info->type == DRMS_TYPE_RAW) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for observer Carrington longitude is of wrong type\n",
	clon);
    fprintf (stderr, "       Must supply an appropriate value for clon_key\n");
    fprintf (stderr, "       (Carrington longitude of disc centre in deg)\n");
    return -1;
  }
  if (strncasecmp (keywd->info->unit, "deg", 3)) {
    fprintf (stderr,
	"Warning: Keyword \"%s\" for observer Carrington longitude has unit of \"%s\"\n",
	clon, keywd->info->unit);
    fprintf (stderr, "         ignored, \"deg\" assumed\n");
  }

  keywd = drms_keyword_lookup (rec, clat, 1);
  if (!keywd) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for observer heliographic latitude not found\n",
	clat);
    fprintf (stderr, "       Must supply an appropriate value for clat_key\n");
    fprintf (stderr, "       (heliographic latitude of disc centre in deg)\n");
    return -1;
  }
  if (keywd->info->type == DRMS_TYPE_TIME ||
      keywd->info->type == DRMS_TYPE_STRING ||
      keywd->info->type == DRMS_TYPE_RAW) {
    fprintf (stderr,
	"Error: Keyword \"%s\" for observer heliographic latitude is of wrong type\n",
	clat);
    fprintf (stderr, "       Must supply an appropriate value for clat_key\n");
    fprintf (stderr, "       (heliographic latitude of disc centre in deg)\n");
    return -1;
  }
  if (strncasecmp (keywd->info->unit, "deg", 3)) {
    fprintf (stderr,
	"Warning: Keyword \"%s\" for observer heliographic latitude has unit of \"%s\"\n",
	clat, keywd->info->unit);
    fprintf (stderr, "         ignored, \"deg\" assumed\n");
  }

  return 0;
}

static int cleanup (int error, DRMS_RecordSet_t *ids, DRMS_RecordSet_t *ods,
    DRMS_Array_t *data, DRMS_Array_t *map) {
  if (data) drms_free_array (data);
  if (map) drms_free_array (map);
  if (ids) drms_close_records (ids, DRMS_FREE_RECORD);
  if (ods) {
    if (error) drms_close_records (ods, DRMS_FREE_RECORD);
    else drms_close_records (ods, DRMS_INSERT_RECORD);
  }
  return error;
}

static void free_all (float *clat, float *clon, double *delta_rot, float *maps,
    float *last, double *maplat, double *maplon, double *map_sinlat,
    unsigned char *offsun) {
  free (clat);
  free (clon);
  free (delta_rot);
  free (maps);
  free (last);
  free (maplat);
  free (maplon);
  free (map_sinlat);
  free (offsun);
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids = NULL, *ods = NULL;
  DRMS_Record_t **orecs;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *record_segment, *oseg, *logseg;
  DRMS_Array_t *data_array = NULL, *map_array = NULL;
  DRMS_Keyword_t *keywd;
  fitsinfo_t *bckfits = NULL;
  FILE **log;
  TIME trec, tobs, tmid, tbase, tfirst, tlast, ttrgt;
  double *maplat, *maplon, *map_coslat, *map_sinlat = NULL;
  double *delta_rot;
  double carr_lon, cm_lon_start, cm_lon_stop, lon_span;
  double cm_lon_first, cm_lon_last;
  double cadence, data_cadence, coverage, t_eps, phase;
  double lat, lon, img_lat, img_lon;
  double tmid_cl, lon_cm;
  double x, y, x0, y0, xstp, ystp;
  double cos_phi, sin_phi;
  double img_xc, img_yc, img_xscl, img_yscl, img_radius, img_pa;
  double ellipse_e, ellipse_pa;
  double kscale;
  float *maps = NULL, *last = NULL;
  float *data, *clat, *clon, *mai, *bck = NULL, *zero;
  enum platloc {LOC_UNKNOWN, LOC_MWO, LOC_GONG_MR, LOC_GONG_LE, LOC_GONG_UD,
      LOC_GONG_TD, LOC_GONG_CT, LOC_GONG_TC, LOC_GONG_BB, LOC_GONG_ML, LOC_SOHO,
      LOC_SDO} platform = LOC_UNKNOWN;
  long long nn;
  unsigned int quality;
  int *reject_list;
  int axes[3], slice_start[3], slice_end[3];
  int recct, rgn, rgnct, segct, valid;
  int tmid_cr, cr_start, cr_stop;
  int cr_first, cr_last;
  int col, row, pixct, dxy, i, found, n, nr, or;
  int blankvals, no_merid_v, rejects, status, verbose_logs, setmais;
  int x_invrt, y_invrt;
  int need_cadence, need_ephem, MDI_correct, MDI_correct_distort;
  unsigned char *offsun;
  char *input, *source, *eoser, *osegname, *pkeyval;
  char logfilename[DRMS_MAXPATHLEN], ctimefmt[DRMS_MAXFORMATLEN];
  char rec_query[256];
  char module_ident[64], key[64], tbuf[64], ptbuf[64], ctime_str[16];

  double raddeg = M_PI / 180.0;
  double degrad = 1.0 / raddeg;
  int need_ephem_from_time = 0;
  int need_crcl = 1;
  int check_platform = 0;
  int scaling_override = 0;
  int segnum = 0;
  int extrapolate = 1;
  int found_first = 0, found_last = 0;
  char *mapname[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
      "LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
      "stereographic", "orthographic", "LambertAzimuthal"};
  char *interpname[] = {"Cubic Convolution", "Nearest Neighbor"};
  char *wcscode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
      "SIN", "ZEA"};
						 /*  process command params  */
  char *inset = params_get_str (params, "in");
  char *outser = params_get_str (params, "out");
  char *bckgnfile = params_get_str (params, "bckgn");
  char *rejectfile = params_get_str (params, "reject");
  char *seg_name = params_get_str (params, "segment");
/*
  unsigned int qmask = params_get_int (params, "qmask");
*/
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  char *tmid_str = params_get_str (params, "tmid");
  TIME tstrt = params_get_time (params, "tstart");
  TIME tstop = params_get_time (params, "tstop");
  int length = params_get_int (params, "length");
  double tstep = params_get_double (params, "tstep");
  int latct = params_get_int (params, "lat_nvals");
  int lonct = params_get_int (params, "lon_nvals");
  int maict = params_get_int (params, "mai_nvals");
  int proj = params_get_int (params, "map");
  int intrpopt = params_get_int (params, "interp");
  double map_scale = params_get_double (params, "scale");
  int map_cols = params_get_int (params, "cols");
  int map_rows = params_get_int (params, "rows");
  double map_pa = params_get_double (params, "map_pa") * raddeg;
  double a0 = params_get_double (params, "a0");
  double a2 = params_get_double (params, "a2");
  double a4 = params_get_double (params, "a4");
  double merid_v = params_get_double (params, "merid_v");
  double bscale = params_get_double (params, "bscale");
  double bzero = params_get_double (params, "bzero");

  char *trec_key = params_get_str (params, "trec_key");
  char *tobs_key = params_get_str (params, "tobs_key");
  char *tstp_key = params_get_str (params, "tstp_key");
  char *qual_key = params_get_str (params, "qual_key");
  char *clon_key = params_get_str (params, "clon_key");
  char *clat_key = params_get_str (params, "clat_key");
  char *crot_key = params_get_str (params, "crot_key");
  char *rsun_key = params_get_str (params, "rsun_key");

  int carr_track = params_isflagset (params, "c");
  int no_track = params_isflagset (params, "n");
  int remove_obsvel = params_isflagset (params, "o");
  int remove_rotation = params_isflagset (params, "r");
  int verbose = params_isflagset (params, "v");
  int dispose = (params_isflagset (params, "x")) ? DRMS_FREE_RECORD :
      DRMS_INSERT_RECORD;
  int MDI_proc = params_isflagset (params, "M");
/*
printf ("time_is_invalid (%s -> %23.16e) = %d\n",
  params_get_str (params, "tstart"), tstrt, time_is_invalid (tstrt));
printf ("drms_ismissing_time (%s) = %d\n", params_get_str (params, "tstart"),
  drms_ismissing_time (tstrt));
return 0;
*/
  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);
  verbose_logs = (dispose == DRMS_INSERT_RECORD) ? verbose : 0;
	  /*  get lists of latitudes and longitudes defining region centers  */
  rgnct = (latct > lonct) ? latct : lonct;
  setmais = (maict == rgnct);
  need_cadence = isnan (tstep);
  if (isnan (map_scale) || map_scale == 0.0) {
    fprintf (stderr,
	"Error: auto scaling from image resolution not implemented;\n");
    fprintf (stderr, "       scale parameter must be set.\n");
    return 1;
  }

  if (no_track) {
    a0 = -(CARR_RATE);
    a2 = a4 = merid_v = 0.0;
    if (carr_track) {
      fprintf (stderr, "Error: inconsistent flags -c and -n\n");
      return 1;
    }
  }
  if (carr_track) a0 = a2 = a4 = merid_v = 0.0;
                                     /*  input expected in microradian/sec   */
  a0 *= 1.0e-6;
  a2 *= 1.0e-6;
  a4 *= 1.0e-6;
  merid_v *= 1.0e-6;
  no_merid_v = (merid_v == 0.0);
  if (map_cols < 1) map_cols = map_rows;
  if (map_rows < 1) map_rows = map_cols;
  if (map_rows < 1) {
    fprintf (stderr, "Error: at least one of \"cols\" or \"rows\" must be set\n");
    return 1;
  }
  MDI_correct = MDI_correct_distort = MDI_proc;
  pixct = map_rows * map_cols;
  cos_phi = cos (map_pa);
  sin_phi = sin (map_pa);
  xstp = ystp = map_scale * raddeg;
  x0 = 0.5 * (1.0 - map_cols) * xstp;
  y0 = 0.5 * (1.0 - map_rows) * ystp;

  if (!(ods = drms_create_records (drms_env, rgnct, outser, DRMS_PERMANENT,
      &status))) {
    fprintf (stderr, "Error: unable to create %d records in series %s\n",
	rgnct, outser);
    fprintf (stderr, "       drms_create_records() returned status %d\n", status); 
    return 1;
  }
  if (verbose) printf ("creating %d record(s) in series %s\n", rgnct, outser);
  if (verbose && dispose == DRMS_FREE_RECORD)
      printf ("experimental run, output records will not be saved\n");
					/*  check output data series struct  */
  orec = drms_recordset_getrec (ods, 0);
  if (!orec) {
    fprintf (stderr, "Error accessing record %d in series %s\n", 0, outser);
    drms_close_records (ods, DRMS_FREE_RECORD);
    return 1;
  }
  segct = drms_record_numsegments (orec);
  found = 0;
  for (n = 0; n < segct; n++) {
    record_segment = drms_segment_lookupnum (orec, n);
    if (record_segment->info->naxis != 3) continue;
    if (record_segment->info->scope == DRMS_CONSTANT) continue;
    if (record_segment->info->scope == DRMS_VARIABLE) {
      if (record_segment->axis[0] != map_cols ||
	  record_segment->axis[1] != map_rows ||
	  record_segment->axis[2] != length) continue;
    }
    if (!found) osegname = strdup (record_segment->info->name);
    found++;
  }
  if (found < 1) {
    fprintf (stderr,
	"Error: no data segment of dimension 3 and appropriate size in output series %s\n", outser);
    drms_close_records (ods, DRMS_FREE_RECORD);
    return 1;
  }
  record_segment = drms_segment_lookup (orec, osegname);
  if (found > 1) {
    fprintf (stderr,
	"Warning: multiple data segments of dimension 3 and appropriate size in output series %s\n", outser);
    fprintf (stderr, "       using \"%s\"\n", osegname);
  }
	    /*  use output series default segment scaling if not overridden  */
  scaling_override = 0;
  if (bscale == 0.0) {
    bscale = record_segment->bscale;
    scaling_override = 1;
  }
  if (isnan (bzero)) {
    bzero = record_segment->bzero;
    scaling_override = 1;
  }
					  /*  check for segment named "Log"  */
  logseg = drms_segment_lookup (orec, "Log");
  if (logseg) drms_segment_filename (logseg, logfilename);
  else if (verbose) {
    fprintf (stderr,
	"Warning: segment \"Log\" not present in output series %s\n", outser);
    fprintf (stderr, "         verbose logging turned off\n");
    verbose_logs = 0;
  }
  if ((keywd = drms_keyword_lookup (orec, "CMLon", 1)))
    sprintf (ctimefmt, "%%d:%s", keywd->info->format);
  else sprintf (ctimefmt, "%%d:%%07.3f");

  if (key_params_from_dspec (inset)) {
				   /*  input specified as specific data set  */
    if (!(ids = drms_open_records (drms_env, inset, &status))) {
      fprintf (stderr, "Error: (%s) unable to open input data set %s\n",
        module_ident, inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    if ((recct = ids->n) < 2) {
      fprintf (stderr, "Error: (%s) <2 records in selected input set\n",
	  module_ident);
      fprintf (stderr, "       %s\n", inset);
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    input = strdup (inset);
    source = strdup (inset);
    eoser = strchr (source, '[');
    if (eoser) *eoser = '\0';
    irec = ids->records[0];
    status = verify_keys (irec, clon_key, clat_key, &kscale);
    if (status) {
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    if (need_cadence) {
			 /*  output cadence not specified, use data cadence  */
      if ((keywd = drms_keyword_lookup (irec, tstp_key, 1))) {
					     /*  cadence should be constant  */
	if (keywd->info->recscope != 1) {
	  fprintf (stderr, "Warning: cadence is variable in input series %s\n",
	      source);
	  fprintf (stderr, "         output cadence will be %f\n",
	       drms_getkey_double (irec, tstp_key, &status));
	}
      } else {
			 /*  could infer from slotting info as well, but...  */
	fprintf (stderr,
          "Error: data cadence keyword %s not in input series %s\n",
	  tstp_key, source);
	fprintf (stderr, "       Specify desired output cadence as tstep\n");
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
      data_cadence = drms_getkey_double (irec, tstp_key, &status);
      tstep = data_cadence;
    }
    tstrt = drms_getkey_time (irec, trec_key, &status);
    tstop = drms_getkey_time (ids->records[recct - 1], trec_key, &status);
    length = (tstop - tstrt + 1.01 * tstep) / tstep;
    tmid = 0.5 * (tstrt + tstop);
    segct = drms_record_numsegments (irec);
  } else {
				/*  only the input data series is named,
				   get record specifications from arguments  */
    source = strdup (inset);
    if (!strcmp (tmid_str, "Not Specified") &&
	(time_is_invalid (tstrt) || time_is_invalid (tstop))) {
      fprintf (stderr,
	  "Error: either a specific data record set must be selected as input\n");
      fprintf (stderr, "       or (tmid and length) or (tstart and tstop) must be\n");
      fprintf (stderr, "       specified\n");
      return 1;
    }
		   /*  get required series info from first record in series  */
					       /*  platform, cadence, phase  */
    snprintf (rec_query, 256, "%s[#^]", inset);
    if (!(ids = drms_open_records (drms_env, rec_query, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", inset);
      fprintf (stderr, "       status = %d\n", status);
      return 1;
    }
    irec = ids->records[0];
    status = verify_keys (irec, clon_key, clat_key, &kscale);
    if (status) {
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    if ((keywd = drms_keyword_lookup (irec, tstp_key, 1))) {
					     /*  cadence should be constant  */
      if (keywd->info->recscope != 1) {
	fprintf (stderr, "Warning: cadence is variable in input series %s\n",
	    source);
	if (need_cadence)
	  fprintf (stderr, "         output cadence will be %f\n",
	      drms_getkey_double (irec, tstp_key, &status));
      }
    } else {
			 /*  could infer from slotting info as well, but...  */
      fprintf (stderr,
          "Error: data cadence keyword %s not in input series %s\n",
	  tstp_key, source);
      fprintf (stderr, "       Specify desired output cadence as tstep\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    data_cadence = drms_getkey_double (irec, tstp_key, &status);
    t_eps = 0.5 * data_cadence;
    if (need_cadence)
			 /*  output cadence not specified, use data cadence  */
      tstep = data_cadence;

    if ((keywd = drms_keyword_lookup (irec, "TELESCOP", 1))) {
						     /*  should be constant  */
      if (keywd->info->recscope != 1)
	fprintf (stderr, "Warning: TELESCOP is variable in input series %s\n",
	    source);
      if (!strcmp (drms_getkey_string (irec, "TELESCOP", &status), "SOHO"))
	platform = LOC_SOHO;
      else if (!strcmp (drms_getkey_string (irec, "TELESCOP", &status), "NSO-GONG")) {
	if ((keywd = drms_keyword_lookup (irec, "SITE", 1))) {
	  if (!strcmp (drms_getkey_string (irec, "SITE", &status), "MR"))
	    platform = LOC_GONG_MR;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "LE"))
	    platform = LOC_GONG_LE;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "UD"))
	    platform = LOC_GONG_UD;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "TD"))
	    platform = LOC_GONG_TD;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "CT"))
	    platform = LOC_GONG_CT;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "BB"))
	    platform = LOC_GONG_BB;
	  else if (!strcmp (drms_getkey_string (irec, "SITE", &status), "ML"))
	    platform = LOC_GONG_ML;
	  else {
	    platform = LOC_GONG_MR;
	  }
	} else {
          fprintf (stderr, "Warning: unspecified GONG site: MR assumed\n");
	  platform = LOC_GONG_MR;
	}
      }
    }
    if (platform == LOC_UNKNOWN) {
      fprintf (stderr, "Warning: observing location unknown, assumed geocenter\n");
    }
    segct = drms_record_numsegments (irec);

    if (strcmp (tmid_str, "Not Specified")) {
/*  determine start and stop times from length (in units of tstep) and midtime
				   (which can be CR:CL as well as date_time) */
       if (sscanf (tmid_str, "%d:%lf", &tmid_cr, &tmid_cl) == 2) {
			  /*  tmid specified as CR:CL : need ephemeris info  */
	need_crcl = 0;
	if (platform == LOC_SOHO ||
	    (platform >= LOC_GONG_MR && platform <= LOC_GONG_ML) ||
	    platform == LOC_UNKNOWN) {
	  if (platform == LOC_SOHO)
	    tmid = SOHO_meridian_crossing (tmid_cl, tmid_cr);
	  else
	    tmid = earth_meridian_crossing (tmid_cl, tmid_cr);
	  sprint_time (ptbuf, tmid, "", 0);
			  /*  adjust to phase of input, within data cadence  */
	  tbase = drms_getkey_time (irec, trec_key, &status);
	  phase = fmod ((tmid - tbase), data_cadence);
	  tmid -= phase;
	  if (phase > 0.5 * data_cadence) tmid += data_cadence;
	  if (verbose) {
	    sprint_time (tbuf, tmid, "", 0);
	    printf ("Target time %d:%05.1f = %s adjusted to %s\n",
	        tmid_cr, tmid_cl, ptbuf, tbuf);
	  }
	} else {
	  fprintf (stderr,
	      "Time specification in CR:CL for observing platform not supported\n");
	  drms_close_records (ids, DRMS_FREE_RECORD);
	  return 1;
	}
      } else {
			      /*  tmid specified as normal date-time string  */
        tmid = sscan_time (tmid_str);
      }
      tbase = drms_getkey_time (irec, trec_key, &status);
      phase = fmod ((tmid - tbase), data_cadence);
      tstrt = tmid - 0.5 * length * tstep;
      tstop = tstrt + (length - 1) * tstep;
			  /*  adjust stop time to reflect sampling symmetry  */
      if ((fabs (phase) < 0.001 * t_eps) && length % 2)
	tstop += tstep;
      if ((fabs (phase - t_eps) < 0.001 * t_eps) && (length % 2 == 0))
	tstop += tstep;
    } else {
	       /*  tstart and tstop specified, determine midtime and length  */
      tmid = 0.5 * (tstrt + tstop);
      length = (tstop - tstrt + 1.01 * tstep) / tstep;
    }
    drms_close_records (ids, DRMS_FREE_RECORD);

    snprintf (rec_query, 256, "%s[?%s > %23.16e and %s < %23.16e?]", inset,
	trec_key, tstrt - t_eps, trec_key, tstop + t_eps);
    if (!(ids = drms_open_records (drms_env, rec_query, &status))) {
      fprintf (stderr, "Error: unable to open input data set %s\n", rec_query);
      fprintf (stderr, "       status = %d\n", status);
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
    if ((recct = ids->n) < 2) {
      fprintf (stderr, "Error: (%s) <2 records in selected input set\n",
	  module_ident);
      fprintf (stderr, "       %s\n", inset);
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }

    input = strdup (rec_query);
    source = strdup (inset);
			    /*  end determination of record set from params  */
  }

  if (verbose) {
    sprint_time (ptbuf, tstrt, "", -1);
    sprint_time (tbuf, tstop, "", -1);
    printf ("tracking data from %s - %s at cadence of %.1f s\n", ptbuf, tbuf,
	tstep);
  }

  if (segct > 1) {
    if (strcmp (seg_name, "Not Specified")) {
      int n;
      segnum = -1;
      for (n = 0; n < segct; n++) {
	record_segment = drms_segment_lookupnum (irec, n);
	if (!strcmp (record_segment->info->name, seg_name)) {
	  segnum = n;
	  break;
	}
      }
      if (segnum < 0) {
	fprintf (stderr,
	    "Error: requested segment %s not found in input series\n",
	    seg_name);
	drms_close_records (ids, DRMS_FREE_RECORD);
	return 1;
      }
    } else {
      fprintf (stderr, "Error: input data series contains multiple segments\n");
      fprintf (stderr, "       segment must be specified\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 1;
    }
  }
  cadence = data_cadence;
  tobs = drms_getkey_time (ids->records[0], tobs_key, &status);
  if (fabs (tobs - tstrt) < cadence) {
    irec = ids->records[0];
    cm_lon_start = drms_getkey_double (irec, clon_key, &status);
    cr_start = drms_getkey_int (irec, crot_key, &status);
    found_first = 1;
  } else need_ephem_from_time = 1;

  tobs = drms_getkey_time (ids->records[recct - 1], tobs_key, &status);
  if (fabs (tobs - tstop) < cadence) {
    irec = ids->records[recct - 1];
    cm_lon_stop = drms_getkey_double (irec, clon_key, &status);
    cr_stop = drms_getkey_int (irec, crot_key, &status);
    found_last = 1;
  } else need_ephem_from_time = 1;

  found = 0;
  for (nr = 0; nr < recct; nr++) {
    tobs = drms_getkey_time (ids->records[nr], tobs_key, &status);
    if (fabs (tobs - tmid) < cadence) {
      irec = ids->records[nr];
      if (need_crcl) {
        tmid_cl = drms_getkey_double (irec, clon_key, &status);
        tmid_cr = drms_getkey_int (irec, crot_key, &status);
      }
      found = 1;
      found_last = 1;
    }
    if (!found_first && !time_is_invalid (tobs)) {
      if (tobs > tmid) found_first = 1;
      else {
	tfirst = tobs;
	irec = ids->records[nr];
	cm_lon_first = drms_getkey_double (irec, clon_key, &status);
	cr_first = drms_getkey_int (irec, crot_key, &status);
      }
    }
    if (!found_last && !time_is_invalid (tobs)) {
      if (tobs > tmid) {
	tlast = tobs;
	irec = ids->records[nr];
	cm_lon_last = drms_getkey_double (irec, clon_key, &status);
	cr_last = drms_getkey_int (irec, crot_key, &status);
	found_last = 1;
      }
    }
  }
  if (!found) need_ephem_from_time = 1;
/*
printf ("recct = %d, length = %d\n", ids->n, length);
*/
  axes[0] = map_cols;
  axes[1] = map_rows;
  axes[2] = 1;
  slice_start[0] = slice_start[1] = 0;
  slice_end[0] = axes[0] - 1;
  slice_end[1] = axes[1] - 1;

  if (need_ephem_from_time) {
    double rsun, vr, vn, vw;
    TIME table_mod_time;
printf ("need ephem from time\n");
    if (platform == LOC_SOHO) {
      soho_ephemeris (tstrt, &rsun, &img_lat, &cm_lon_start, &vr, &vn, &vw,
	  &table_mod_time);
      cr_start = carrington_rots (tstrt, 1);
      soho_ephemeris (tmid, &rsun, &img_lat, &tmid_cl, &vr, &vn, &vw,
	  &table_mod_time);
      tmid_cr = carrington_rots (tmid, 1);
      soho_ephemeris (tstop, &rsun, &img_lat, &cm_lon_stop, &vr, &vn, &vw,
	  &table_mod_time);
      cr_stop = carrington_rots (tstop, 1);
    } else {
      if (!found_first || !found_last) {
        fprintf (stderr, "Error: Carrington ephemeris from time not supported\n");
        fprintf (stderr, "         and no valid times in data for estimation!\n");
	return 1;
      }
		/*  estimate midpoint ephemeris by linear interpolation
					 of closest observations to midtime  */
 /*  This code has not been well tested in crossover of Carrington rotation  */
      fprintf (stderr, "Warning: Carrington ephemeris from time not supported\n");
      tmid_cr = cr_first;
    	  /*  assume at most one rotation between first and last estimators  */
      if (cr_last != cr_first) {
        cm_lon_last -= 360.0;
      }
      tmid_cl = cm_lon_first +
          (cm_lon_last - cm_lon_first) * (tmid - tfirst) / (tlast - tfirst);
    	      /*  assume at most one rotation between first and last estimators  */
      if (tmid_cl < 0.0) {
	tmid_cr++;
	tmid_cl += 360.0;
      }
      fprintf (stderr, "         estimating midpoint as %d:%08.4f\n",
	  tmid_cr, tmid_cl);
	/*  long extrapolations, and lazy correction for change of rotation
					number,but only needed for lon_span  */
      cr_start = cr_first;
      cm_lon_start = cm_lon_first +
          (cm_lon_last - cm_lon_first) * (tstrt - tfirst) / (tlast - tfirst);
      cr_stop = cr_last;
      cm_lon_stop = cm_lon_first +
          (cm_lon_last - cm_lon_first) * (tstop - tfirst) / (tlast - tfirst);
      if (cm_lon_stop > cm_lon_start) cr_stop++;
    }
  }
  lon_span = cm_lon_start - cm_lon_stop;
  while (cr_start < cr_stop) {
    cr_start++;
    lon_span += 360.0;
  }
					  /*  allocate map parameter arrays  */
  clat = (float *)malloc (rgnct * sizeof (float));
  clon = (float *)malloc (rgnct * sizeof (float));
  mai = (float *)malloc (rgnct * sizeof (float));
  delta_rot = (double *)malloc (rgnct * sizeof (double));
  for (i = 0; i < latct; i++) {
    snprintf (key, 64, "lat_%d_value", i);
    clat[i] = params_get_float (params, key) * raddeg;
  }
  for (i = 0; i < lonct; i++) {
    snprintf (key, 64, "lon_%d_value", i);
    clon[i] = params_get_float (params, key) * raddeg;
  }
  for (i = latct; i < rgnct; i++) clat[i] = clat[i-1];
  for (i = lonct; i < rgnct; i++) clon[i] = clon[i-1];
  if (setmais) {
    for (i = 0; i < latct; i++) {
      snprintf (key, 64, "mai_%d_value", i);
      mai[i] = params_get_float (params, key);
    }
  }
  for (rgn = 0; rgn < rgnct; rgn++) {
    double sin_lat = sin (clat[rgn]);
    sin_lat *= sin_lat;
    delta_rot[rgn] =  a0 + sin_lat * (a2 + a4 * sin_lat);
  }
  log = (FILE **)malloc (rgnct * sizeof (FILE *));
	/*  create output memory structures: data arrays and file pointers,
					     one per region (output record)  */
  orecs = (DRMS_Record_t **)malloc (rgnct * sizeof (DRMS_Record_t *));
  maps = (float *)malloc (pixct * rgnct * sizeof (float));
  last = (float *)malloc (pixct * rgnct * sizeof (float));
		      /*  allocate mapping info arrays for all output cubes  */
  maplat = (double *)malloc (pixct * rgnct * sizeof (double));
  maplon = (double *)malloc (pixct * rgnct * sizeof (double));
  if (no_merid_v) {
    map_coslat = maplat;
    map_sinlat = (double *)malloc (pixct * rgnct * sizeof (double));
  }
  offsun = (unsigned char *)malloc (pixct * rgnct * sizeof (char));
    /*  Calculate heliographic coordinates corresponding to map location(s)  */
  n = 0;
  for (rgn = 0; rgn < rgnct; rgn++) {
    double xrot, yrot;
    for (row=0, y =y0; row < map_rows; row++, y +=ystp) {
      for (col=0, x =x0; col < map_cols; col++, x +=xstp, n++) {
	xrot = x * cos_phi - y * sin_phi;
	yrot = y * cos_phi + x * sin_phi;
	offsun[n] = plane2sphere (xrot, yrot, clat[rgn], clon[rgn], &lat, &lon,
	    proj);
	maplat[n] = lat;
	maplon[n] = lon;
	if (no_merid_v) {
	  map_coslat[n] = cos (lat);
	  map_sinlat[n] = sin (lat);
	}
      }
    }
  }
				    /*  this should not really be necessary  */
  zero = (float *)calloc (pixct, sizeof (float));
  map_array = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, (void *)zero, &status);
  map_array->bscale = bscale;
  map_array->bzero = bzero;

  dxy = 0;
  if (strcmp (bckgnfile, "Not Specified")) {
                                       /*  read in average background image  */
    bckfits = rsb_read_fits (bckgnfile);
    if (bckfits) {
      if (bckfits->naxis) dxy = 1;
      for (n = 0; n < bckfits->naxis; n++)
        dxy *= bckfits->axis[n];
      if (bckfits->bitpix == -32) bck = (float *)bckfits->data;
      else {
        bck = (float *)malloc (dxy * sizeof (float));
	if (bckfits->bitpix == 8) {
	  char *buf = (char *)bckfits->data;
	  for (nn = 0; nn < dxy; nn++) bck[nn] = (buf[nn] == bckfits->blank) ?
	      FP_NAN : bckfits->bscale * buf[nn] + bckfits->bzero;
	} else if (bckfits->bitpix == 16) {
	  short *buf = (short *)bckfits->data;
	  for (nn = 0; nn < dxy; nn++) bck[nn] = (buf[nn] == bckfits->blank) ?
	      FP_NAN : bckfits->bscale * buf[nn] + bckfits->bzero;
	} else if (bckfits->bitpix == 32) {
	  int *buf = (int *)bckfits->data;
	  for (nn = 0; nn < dxy; nn++) bck[nn] = (buf[nn] == bckfits->blank) ?
	      FP_NAN : bckfits->bscale * buf[nn] + bckfits->bzero;
	} else if (bckfits->bitpix == -64) {
	  double *buf = (double *)bckfits->data;
	  for (nn = 0; nn < dxy; nn++)
	    bck[nn] = buf[nn];
	}
      }
    } else fprintf (stderr,
	  "Warning: unable to open background file %s; ignored\n", bckgnfile);
  }
		 /*  support special hack of reading of rejection list file  */
  rejects = 0;
  if (strcmp (rejectfile, "Not Specified")) {
    FILE *rejectfp = fopen (rejectfile, "r");
    if (rejectfp) rejects = read_reject_list (rejectfp, &reject_list);
    else fprintf (stderr,
	"Warning: could not open rejection list %s; ignored\n", rejectfile);
  }

  found = 0;
  for (rgn = 0; rgn < rgnct; rgn++) {
/*
    map[rgn] = (float *)malloc (pixct * sizeof (float));
*/
    orec = orecs[rgn] = drms_recordset_getrec (ods, rgn);
    if (!orec) {
      fprintf (stderr, "Error accessing record %d in series %s\n", rgn, outser);
      drms_close_records (ids, DRMS_FREE_RECORD);
      drms_close_records (ods, DRMS_FREE_RECORD);
      free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
      if (bckfits) {
	rsb_close_fits (bckfits);
	if (bck) free (bck);
      }
      return 1;
    }
    if (verbose_logs) {
      log[rgn] = fopen (logfilename, "w");
      if (!log[rgn]) {
        fprintf (stderr, "Error: unable to open log file %s\n", logfilename);
	drms_close_records (ids, DRMS_FREE_RECORD);
	drms_close_records (ods, DRMS_FREE_RECORD);
	free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
	if (bckfits) {
	  rsb_close_fits (bckfits);
	  if (bck) free (bck);
	}
	return 1;
      }
    }
		      /*  set segment axis lengths and scaling if necessary  */
    oseg = drms_segment_lookup (orec, osegname);
    if (oseg->info->scope == DRMS_VARDIM) {
      oseg->axis[0] = map_cols;
      oseg->axis[1] = map_rows;
      oseg->axis[2] = length;
    }
    if (scaling_override) {
      oseg->bscale = bscale;
      oseg->bzero = bzero;
    }
  }
					     /*  loop through input records  */
  valid = 0;
  ttrgt = tstrt;
  or = 0;
  for (nr = 0; nr < recct; nr++) {
    irec = ids->records[nr];
    tobs = drms_getkey_time (irec, tobs_key, &status);
			      /*  replace with call to solar_ephemeris_info  */
    img_lon = drms_getkey_double (irec, clon_key, &status);
    img_lat = drms_getkey_double (irec, clat_key, &status);
    if (time_is_invalid (tobs)) {
							  /*  no data, skip  */
      continue;
    }
    quality = drms_getkey_int (irec, qual_key, &status);
    if ((quality & qmask) && !status) {
						    /*  partial image, skip  */
      continue;
    }
    blankvals = drms_getkey_int (irec, "MISSVALS", &status);
    if (blankvals && !status) {
						    /*  partial image, skip  */
      continue;
    }
			 /*  check for record quality, reject as applicable  */
    if (rejects) {
      int idrec = drms_getkey_int (irec, "T_REC_index", &status);
      int match = 0;
      if (status) {
	fprintf (stderr, "Warning: \"T_REC_index\" keyword not found\n");
	fprintf (stderr, "         up to %d bad images could be processed\n",
	  rejects);
	rejects = 0;
      }
      for (n = 0; n < rejects; n++) {
	if (idrec == reject_list[n]) {
	  match = 1;
	  break;
	}
      }
      if (match) {
/*
	if (verbose) 
	  printf ("  rejecting img %d -- in list %s\n", idrec, rejectfile);
*/
        continue;
      }
    }
			   /*  get needed info from record keys for mapping  */
    status = solar_image_info (irec, &img_xscl, &img_yscl, &img_xc, &img_yc,
	&img_radius, rsun_key, &img_pa, &ellipse_e, &ellipse_pa, &x_invrt,
	&y_invrt, &need_ephem);
    if (status) {
      fprintf (stderr, "Error: one or more essential keywords or values missing; skipped\n");
      fprintf (stderr, "solar_image_info() returned %08x\n", status);
      continue;
    }
    if (MDI_correct) {
      mtrack_MDI_correct_imgctr (&img_xc, &img_yc, img_radius);
      mtrack_MDI_correct_pa (&img_pa);
    }
						  /*  read input data image  */
    record_segment = drms_segment_lookupnum (irec, segnum);
    data_array = drms_segment_read (record_segment, DRMS_TYPE_FLOAT, &status);
    data = (float *)data_array->data;
				   /*  remove average background from image  */
    if (dxy) {
      for (nn = 0; nn < dxy; nn++) {
        if (isnan (data[nn])) continue;
	data[nn] -= bck[nn];
      }
    }

    sprint_time (tbuf, tobs, "TAI", 0);
    tbuf[strlen (tbuf) - 4] = '\0';
    img_xc -= 0.5 * (data_array->axis[0] - 1);
    img_yc -= 0.5 * (data_array->axis[1] - 1);
		 	/*  should be taken care of in solar_ephemeris_info  */
img_lon *= raddeg;
img_lat *= raddeg;
    if (need_ephem) {
      ; 
    }

    if (!extrapolate)
      memcpy (last, maps, rgnct * pixct * sizeof (float));
		      /*  loop through output records, appending time slice  */
    perform_mappings (data_array, maps, delta_rot, rgnct, maplat,
	maplon, map_coslat, map_sinlat, pixct, tobs - tmid, offsun,
	img_lat, img_lon, img_xc, img_yc, a0, a2, a4, merid_v, img_radius,
	img_pa, ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt,
	MDI_correct_distort);
/*
printf ("%4d: %s %.4f %.5f %7.2f %7.2f %7.2f %7.2f\n", nr, tbuf, img_lon, img_lat,
maps[0], maps[pixct - 1], maps[(rgnct-1)*pixct], maps[rgnct*pixct - 1]);
    if (add_before) {
printf ("adding %d records before first of %d\n", add_before, recct);
      while (add_before--) {
	for (rgn = 0; rgn < rgnct; rgn++) {
	  memcpy (&map[rgn][(length - nr - add_before - 1)*pixct], &maps[rgn*pixct],
	      pixct * sizeof (float));
	  if (verbose) {
	    fprintf (log[rgn], "record %d extrapolated from image %s\n",
		length - nr - add_before - fill_ct - 1, tbuf);
	  }
	}
      }
      add_before = 0;
    }
*/
	    /*  extrapolate first image backward to start time if necessary  */
    if (extrapolate) {
      while (ttrgt < tobs) {
	if (verbose) printf ("step %d extrapolated from image %s\n", or, tbuf);
	for (rgn = 0; rgn < rgnct; rgn++) {
	  orec = orecs[rgn];
	  oseg = drms_segment_lookup (orec, osegname);
	  memcpy (map_array->data, &maps[rgn*pixct], pixct * sizeof (float));
	  slice_start[2] = slice_end[2] = or;
	  status = drms_segment_writeslice (oseg, map_array, slice_start,
	      slice_end, 0);
	  if (status) {
	    fprintf (stderr, "Error writing data to record %d in series %s\n",
		rgn, outser);
	    fprintf (stderr, "      series may not have appropriate structure\n");
	    drms_free_array (map_array);
	    drms_close_records (ods, DRMS_FREE_RECORD);
	    free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
	    return 1;
	  }
	  if (verbose_logs) {
	    fprintf (log[rgn], "step %d extrapolated from image %s\n", or,
		tbuf);
	  }
	}
	or++;
	if (or >= length) {
	  fprintf (stderr, "Error: reached output length limit\n");
	  drms_close_records (ids, DRMS_FREE_RECORD);
	  drms_close_records (ods, DRMS_FREE_RECORD);
	  free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
	  if (bckfits) {
	    rsb_close_fits (bckfits);
	    if (bck) free (bck);
	  }
	  return 1;
	}
	ttrgt += tstep;
      }
      extrapolate = 0;
    }
    if (ttrgt < tobs) {
	/*  linearly interpolate individual pixels between last valid map
								and current  */
      double f, g;
      int ct, ntot = rgnct * pixct;
      float *val = (float *)malloc (ntot * sizeof (float));
      char *skip = (char *)malloc (ntot * sizeof (char));
      
      for (n = 0; n < ntot; n++) {
	skip[n] = (isnan (last[n]) || isnan (maps[n]));
	val[n] = missing_val;
      }
      while (ttrgt < tobs) {
        f = (ttrgt - tlast) / (tobs - tlast);
	g = 1.0 - f;
  	for (n = 0; n < ntot; n++)
	  if (!skip[n]) val[n] = g * last[n] + f * maps[n];
	if (verbose) printf ("step %d interpolated from images %s and %s\n",
	      or, ptbuf, tbuf);
	for (rgn = 0; rgn < rgnct; rgn++) {
	  orec = orecs[rgn];
	  oseg = drms_segment_lookup (orec, osegname);
	  memcpy (map_array->data, &val[rgn*pixct], pixct * sizeof (float));
	  slice_start[2] = slice_end[2] = or;
	  status = drms_segment_writeslice (oseg, map_array, slice_start,
	      slice_end, 0);
	  if (status) {
	    fprintf (stderr, "Error writing data to record %d in series %s\n",
		rgn, outser);
	    fprintf (stderr, "      series may not have appropriate structure\n");
	    drms_free_array (map_array);
	    drms_close_records (ods, DRMS_FREE_RECORD);
	    free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
	    return 1;
	  }
	  if (verbose_logs) fprintf (log[rgn],
	      "step %d interpolated from images %s and %s\n", or, ptbuf, tbuf);
	}
	or++;
	if (or >= length) {
	  if (nr < recct - 1)
	    fprintf (stderr,
	      "Warning: reached output limit before last input record processed\n");
	  ttrgt = tstop;
	}
	ttrgt += tstep;
      }
      free (val);
      free (skip);
    }

    if (ttrgt == tobs) {
      if (verbose) printf ("step %d mapped from image %s\n", or, tbuf);
      for (rgn = 0; rgn < rgnct; rgn++) {
/*
	memcpy (&map[rgn][nr*pixct], &maps[rgn*pixct], pixct * sizeof (float));
*/
	orec = orecs[rgn];
	oseg = drms_segment_lookup (orec, osegname);
	memcpy (map_array->data, &maps[rgn*pixct], pixct * sizeof (float));
	slice_start[2] = slice_end[2] = or;
	status = drms_segment_writeslice (oseg, map_array, slice_start,
	    slice_end, 0);
	if (status) {
	  fprintf (stderr, "Error writing data to record %d in series %s\n",
	      rgn, outser);
	  fprintf (stderr, "      series may not have appropriate structure\n");
	  drms_free_array (map_array);
	  drms_close_records (ods, DRMS_FREE_RECORD);
	  free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
	  return 1;
	}
	if (verbose_logs)
          fprintf (log[rgn], "step %d mapped from image %s\n", or, tbuf);
      }
      or++;
      ttrgt += tstep;
    }
    tlast = tobs;
    strcpy (ptbuf, tbuf);
    drms_free_array (data_array);
    data_array = NULL;
    valid++;
  }
  drms_close_records (ids, DRMS_FREE_RECORD);
  if (bckfits) {
    rsb_close_fits (bckfits);
    if (bck) free (bck);
  }
					 /*  extend last image if necessary  */
  while (ttrgt <= tstop && or < length) {
    if (verbose) printf ("step %d extrapolated from image %s\n", or, tbuf);
    for (rgn = 0; rgn < rgnct; rgn++) {
/*
      memcpy (&map[rgn][or*pixct], &maps[rgn*pixct], pixct * sizeof (float));
*/
      orec = orecs[rgn];
      oseg = drms_segment_lookup (orec, osegname);
      memcpy (map_array->data, &maps[rgn*pixct], pixct * sizeof (float));
      slice_start[2] = slice_end[2] = or;
      status = drms_segment_writeslice (oseg, map_array, slice_start,
	  slice_end, 0);
      if (status) {
	fprintf (stderr, "Error writing data to record %d in series %s\n",
	    rgn, outser);
	fprintf (stderr, "      series may not have appropriate structure\n");
	drms_free_array (map_array);
	drms_close_records (ods, DRMS_FREE_RECORD);
	free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
	return 1;
      }
      if (verbose_logs)
	  fprintf (log[rgn], "step %d extrapolated from image %s\n", or, tbuf);
    }
    ttrgt += tstep;
    or++;
  }

  coverage = (double)valid / (double)length;
						      /*  write out records  */
  for (rgn = 0; rgn < rgnct; rgn++) {
							 /*  set key values  */
    int kstat = 0;
    orec = orecs[rgn];
    kstat += check_and_set_key_str   (orec, "WCSNAME", "Carrington Heliographic/Time");
    kstat += check_and_set_key_int   (orec, "WCSAXES", 3);
    snprintf (key, 64, "CRLN-%s", wcscode[proj]);
    kstat += check_and_set_key_str   (orec, "CTYPE1", key);
    snprintf (key, 64, "CRLT-%s", wcscode[proj]);
    kstat += check_and_set_key_str   (orec, "CTYPE2", key);
    kstat += check_and_set_key_str   (orec, "CTYPE3", "TIME");
    kstat += check_and_set_key_str   (orec, "CUNIT1", "deg");
    kstat += check_and_set_key_str   (orec, "CUNIT2", "deg");
    kstat += check_and_set_key_str   (orec, "CUNIT3", "s");
    kstat += check_and_set_key_float (orec, "CRPIX1", 0.5 * map_cols + 0.5);
    kstat += check_and_set_key_float (orec, "CRPIX2", 0.5 * map_rows + 0.5);
    kstat += check_and_set_key_float (orec, "CRPIX3", 0.5 * length + 0.5);
    kstat += check_and_set_key_float (orec, "CRVAL1", clon[rgn] * degrad);
    kstat += check_and_set_key_float (orec, "CRVAL2", clat[rgn] * degrad);
    kstat += check_and_set_key_double(orec, "CRVAL3", tmid);
    kstat += check_and_set_key_float (orec, "CDELT1", map_scale);
    kstat += check_and_set_key_float (orec, "CDELT2", map_scale);
    kstat += check_and_set_key_float (orec, "CDELT3", tstep);
    kstat += check_and_set_key_float (orec, "LonHG", clon[rgn] * degrad);
    kstat += check_and_set_key_float (orec, "LatHG", clat[rgn] * degrad);
    kstat += check_and_set_key_str   (orec, "MapProj", mapname[proj]);
    kstat += check_and_set_key_str   (orec, "Interp", interpname[intrpopt]);
    kstat += check_and_set_key_time  (orec, "MidTime", tmid);
    snprintf (ctime_str, 16, ctimefmt, tmid_cr, tmid_cl);
    kstat += check_and_set_key_str   (orec, "CarrTime", ctime_str);
    kstat += check_and_set_key_int   (orec, "CarrRot", tmid_cr);
    kstat += check_and_set_key_float (orec, "CMLon", tmid_cl);
    carr_lon = (360.0 * tmid_cr) + 360.0 - tmid_cl;
    kstat += check_and_set_key_double(orec, "CarrLon", carr_lon);
    lon_cm = clon[rgn] * degrad - tmid_cl;
    while (lon_cm < -180.0) lon_cm += 360.0;
    while (lon_cm > 180.0) lon_cm -= 360.0;
    kstat += check_and_set_key_float (orec, "LonCM", lon_cm);
    kstat += check_and_set_key_time  (orec, "T_START", tstrt);
    kstat += check_and_set_key_time  (orec, "T_STOP", tstop);
    kstat += check_and_set_key_float (orec, "LonSpan", lon_span);
    kstat += check_and_set_key_float (orec, "Coverage", coverage);
    kstat += check_and_set_key_float (orec, "ZonalTrk", delta_rot[rgn] * 1.0e6);
    kstat += check_and_set_key_float (orec, "ZonalVel",
	(delta_rot[rgn] + CARR_RATE * 1.0e-6) * RSUNM * cos (clat[rgn]));
    kstat += check_and_set_key_float (orec, "LonDrift",
	delta_rot[rgn] * degrad * length * tstep);
    kstat += check_and_set_key_float (orec, "MeridTrk", merid_v * 1.0e6);
    kstat += check_and_set_key_float (orec, "MeridVel", merid_v * RSUNM);
    kstat += check_and_set_key_float (orec, "LatDrift",
	merid_v * degrad * length * tstep);
    kstat += check_and_set_key_str   (orec, "Module", module_ident);
    kstat += check_and_set_key_str   (orec, "Source", source);
    kstat += check_and_set_key_str   (orec, "Input", input);
    if (strcmp (bckgnfile, "Not Specified"))
      kstat += check_and_set_key_str (orec, "Backgrnd", bckgnfile);
    if (strcmp (rejectfile, "Not Specified"))
      kstat += check_and_set_key_str (orec, "RejectList", rejectfile);
    kstat += check_and_set_key_float (orec, "MapScale", map_scale);
    kstat += check_and_set_key_float (orec, "Cadence", tstep);
    kstat += check_and_set_key_float (orec, "Duration", length * tstep);
    kstat += check_and_set_key_float (orec, "Width", map_cols * map_scale);
    kstat += check_and_set_key_float (orec, "Height", map_rows * map_scale);
    kstat += check_and_set_key_float (orec, "Size", sqrt (map_rows * map_cols) * map_scale);
    kstat += check_and_set_key_float (orec, "PosAng", map_pa / raddeg);
    kstat += check_and_set_key_float (orec, "RSunRef", 1.0e-6 * RSUNM);
    if (setmais && isfinite (mai[rgn]))
      kstat += check_and_set_key_float (orec, "MAI", mai[rgn]);
    if (scaling_override) {
      sprintf (key, "%s_bscale", osegname);
      kstat += check_and_set_key_double(orec, key, bscale);
      sprintf (key, "%s_bzero", osegname);
      kstat += check_and_set_key_double(orec, key, bzero);
    }
    if (MDI_correct)
      kstat += check_and_set_key_float (orec, "MDI_PA_Corr", MDI_IMG_SOHO_PA);
    if (pkeyval = create_prime_key (orec))
      kstat += check_and_set_key_str (orec, "PrimeKeyString", pkeyval);
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record %d in series %s\n",
	  rgn, outser);
      fprintf (stderr, "      series may not have appropriate structure\n");
      drms_free_array (map_array);
      drms_close_records (ods, DRMS_FREE_RECORD);
      free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
      return 1;
    }
    if (verbose_logs) fclose (log[rgn]);
  }
  drms_close_records (ods, dispose);
  drms_free_array (map_array);
  free_all (clat, clon, delta_rot, maps, last, maplat, maplon, map_sinlat, offsun);
  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  08.04.22 	created this file, based on fastrack
 *  v 0.1	add all args, modify as appropriate for DRMS;
 *	add all output keywords, fill as possible from args (except dataset)
 *    0.2	add input dataset processing, fill output keys as necessary
 *    0.3	process input record list (dataset) info, reflect in keys
 *	add platform info and ephemeris for SOHO
 *    0.4	open and close input segments, add background image processing
 *    0.5	create and write to output segments
 *    0.6	create interpolation targets, fix some key value bugs
 *    0.7	interpolate input data to output maps using crude keyword
 *	search; add interpolation of skipped images; close most memory leaks;
 *	fixed possible fastrack bug in calculation of f in interpolate_images()
 *    0.8	added keywords for output cadence, start/stop specification,
 *	interpolation option;
 *	extensive modifications to allow for time interpolation to different
 *		cadence from input, missing records and non-uniformly spaced
 *		records in input;
 *	removed restriction on output segment name;
 *	added correction for MDI image distortion and PA and image center
 *		(with -M flag), but not yet correction for plate scale if
 *		necessary;
 *	added processing of personal MDI image quality rejection lists
 *	added setting of optional PrimeKeyString
 *	adjust midtime to phase of input when specified as CR:CL
 *    0.9	write output by slices, avoiding massive memory allocation
 *	added code for nearest-value "interpolation" (untested)
 *	added processing of optional OBS_ASD key in place of R_SUN if missing
 *	changed defaults for output scaling to use output series defaults
 *	change Duration key type from int to float
 *	reject images with MISSVALS != 0 (if present)
 *	fixed setting of ZonalTrk and MeridTrk keys, added corresponding Vel
 *		keys and RSunRef; added MDI_PA_Corr key; fixed setting of
 *		Duration and CDELT3 keys
 *	adjusted Carrington rate (by 5 parts in 10^9)
 *	fixed setting of scaling when overrides in effect; support vardim
 *		output by setting segment axis
 *	fixed a few icc11 compiler warnings
 *	added argument rsun_key; verification of some essential keys
 *	added recognition of GONG site data origins, support for JPL table-based
 *		ephemerides, processing of quality keyword against qmask
 *	added optional processing of MAI value(s) into keyword
 *	changed constituent of primekeystring from CarrTime to CarrRot and
 *		CMLon
 *	added setting of Size keyword (geometric mean of Height and Width) and
 *		Cadence (= CDELT3)
 *  v 0.9 frozen 2010.02.08
 *
 */
