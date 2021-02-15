/******************************************************************************/
/*
 *  maproj.c						~rick/src/maproj
 *
 *  Responsible:  Rick Bogart				rick@sun.stanford.edu
 *
 *  Map a set of input solar images into a set of output maps
 *
 *  Parameters:	(type	default		description)
 *	in	DataSet -		Input dataset
 *				A set of images of all or part of the
 *				solar disc in a "plate" coordinate system
 *				(helioprojective geometry).
 *      out     DataSer -		Output data series name
 *	clat	double	NaN		Map central heliographic latitude
 *				(defaults to lat of [first] image center)
 *	clon	double	NaN		Map central heliographic longitude
 *				(defaults to lon of [first] image center)
 *      scale   double	0.0		Scale of map (heliographic degrees /
 *				pixel) at location appropriate for mapping
 *				option; a 0 value implies autoscaling to best
 *	cols	int	0		Columns in output maps; 0 -> rows
 *	rows	int	0		Rows in output maps; 0 -> cols
 *      map     enum  "Postel"		Mapping option:
 *				recognized values are "carree", "Cassini",
 *				"Mercator", "cyleqa", "sineqa", "gnomonic",
 *				"Postel", "stereographic", "orthographic",
 *				and "Lambert" (and possibly others).
 *      interp  enum  "cubiconv"	Interpolation option:
 *				recognized values are "cubiconv", "nearest",
 *				and "bilinear" (and possibly others)
 *	grid	float  Unspec		If supplied, the spacing in degrees of
 *				a latitude/longitude grid to be overlain on the
 *				output map(s). The overlay value is -Inf where
 *				there would be valid data, +Inf where there is
 *				not. Points are considered on a grid line if
 *				they are within grid_ratio * the grid spacing from it
 *	grid_ratio float 0.01		Ratio of grid line width to spacing
 *      map_pa  float   0.0             The angle between heliographic north
 *				and "up" on the output map (in the direction
 *				of increasing rows) [deg[, in the sense that a
 *				positive position angle represents a clockwise
 *				displacement of the north axis.
 *	bscale	float	0.0		Value scaling parameter for output
 *				defaults to output series default
 *	bzero	float	NaN		Value offset parameter for output
 *				defaults to output series default
 *	ldcoef	float*	(1.0}		coefficients for polynomial
 *				limb-darkening function
 *	qmask	int	0x80000000	Quality mask for data acceptability;
 *				records rejected if (qmask & qkey:value) != 0
 *	qual_key str	Quality		Key name of uint type keyword describing
 *			data quality
 *	clon_key string	CRLN_OBS	Keyname of float type keyword describing
 *				centre Carrington longitude of each input image
 *	clat_key string	CRLT_OBS	Keyname of float type keyword describing
 *				centre Carrington latitude of each input image
 *	rsun_key string	R_SUN		Keyname of float type keyword describing
 *				apparent solar semidiameter of image [pixel]
 *	apsd_key string	RSUN_OBS	Keyname of float type keyword describing
 *				apparent solar semidiameter of image [arcsec]
 *	dsun_key string	DSUN_OBS	Keyname of double type keyword describing
 *				r distance from sun for of each image
 *
 *  Flags
 *	-c	center map(s) at image center(s)
 *	-n	no output: diagnostics only
 *	-s	interpret clon as Stonyhurst rather than Carrington longitude
 *	-v	run verbose
 *	-M	correct for MDI distortion
 *	-R	convert scalar values assuming that they are the line-of-sight
 *		components of radial vectors (i.e. multiply by secant of
 *		center-to-limb angle)
 *
 *  Bugs:
 *    No provision for propagation of default or selected keywords; only tries
 *	to copy all common keys before explicit setting
 *    Uses considerable replicated code from mtrack, esp. array_imaginterp()
 *	function and its dependencies; should be consolidated
 *    No provision for anisotropic scaling (CDELT1 != CDELT2); PCi_j not
 *	adjusted either
 *    There is evidently no WCS conventional name for the Cassini-Soldner
 *	(transverse plate carree) projection; CAS is arbitrarily used; the
 *	alternative would be to interchange HGLT and HGLN, but that would
 *	necessitate a change in the position angle
 *
 *  Future Updates
 *    Transformation from one mapping to another
 *
 *  Revision history is at end of file
 */
/******************************************************************************/

#define MODULE_VERSION_NUMBER	("1.4")
#define CARTOGRAPHY_VERSION "cartography_v11.c"
#define KEYSTUFF_VERSION "keystuff_v11.c"

#include <jsoc_main.h>
#include CARTOGRAPHY_VERSION
#include KEYSTUFF_VERSION
#include "imginfo.c"
#include "mdistuff.c"

#define RSUNM		(6.96e8)
#define INTERP_NEAREST_NEIGHBOR	(1)
#define INTERP_BILINEAR	(2)
						       /*  module identifier  */
char *module_name = "maproj";
char *module_desc = "mapping from solar images";
char *version_id = MODULE_VERSION_NUMBER;

ModuleArgs_t module_args[] = {
  {ARG_DATASET,	"in", "", "input data set"}, 
  {ARG_DATASERIES, "out", "", "output data series"}, 
  {ARG_DOUBLE,	"clat", "Not Specified", "heliographic latitude of map center [deg]"},
  {ARG_DOUBLE,	"clon", "Not Specified", "Carrington longitude of map center [deg]"},
  {ARG_DOUBLE,	"scale", "Not specified", "map scale at center [deg/pxl]"},
  {ARG_INT,     "cols", "0", "columns in output map"},
  {ARG_INT,     "rows", "0", "rows in output map"},
  {ARG_NUME,	"map", "orthographic", "map projection",
      "carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
  {ARG_NUME,	"interp", "cubiconv", "interpolation option",
      "cubiconv, nearest, bilinear"},
  {ARG_FLOAT,	"grid", "Not Specified",
      "if specified, spacing of grid overlay [deg]"},
  {ARG_FLOAT,	"grid_ratio", "0.01", "ratio of grid line width to spacing"},
  {ARG_FLOAT,	"map_pa", "0.0", "position angle of north on output map [deg]"},
  {ARG_FLOAT,	"bscale", "0.0", "output scaling factor"},
  {ARG_FLOAT,	"bzero", "Default", "output offset"},
  {ARG_FLOATS,	"ldcoef", "{1.0}",
      "coefficients in limb darkening polynomial in mu"},
  {ARG_INT,	"qmask", "0x80000000", "quality bit mask for image rejection"},
  {ARG_STRING,	"qual_key", "Quality", "keyname for 32-bit image quality field"},
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"rsun_key", "R_SUN", "keyname for image semi-diameter (pixel)"}, 
  {ARG_STRING,	"apsd_key", "RSUN_OBS", "keyname for apparent solar semi-diameter (arcsec)"}, 
  {ARG_STRING,	"dsun_key", "DSUN_OBS", "keyname for observer distance"}, 
  {ARG_STRING,	"RequestID", "none", "RequestID for jsoc export management"},
  {ARG_FLAG,	"c",	"", "center map at center of image"}, 
  {ARG_FLAG,	"n",	"", "no output records produced; diagnostics only"}, 
  {ARG_FLAG,	"s",	"", "clon is Stonyhurst longitude"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {ARG_FLAG,	"M",	"", "correct for MDI distortion"},
  {ARG_FLAG,	"R",	"",
      "convert scalars assuming line-of-sight components of radial vectors"}, 
  {}
};

		 /*  global declaration of missing to be initialized as NaN  */
float missing_val;

double limbdark (double mu, double *coef, int order) {
/*
 *  Approximation to a limb-darkening function in powers of mu
 *    [cos (central angle)]; order >= 0 with array coef[order+1]
 *  For HMI continuum, the first order function 0.5 + 0.5 * mu does a better
 *    job of at least visually flattening than the canonical second order
 *    function 0.3 + 0.93 * mu - 0.23 * mu^2. A piecewsie linear fit:
 *	0.5 + 0.5 * mu (mu >= 0.5),
 *	0.375 + 0.75 * mu (0 <= mu < 0.5)
 *    does an even better job.
 */
  double fc, mup;
  int n;

  fc = 0.0;
  mup = 1.0;
  for (n = 0; n <= order; n++) {
    fc += coef[n] * mup;
    mup *= mu;
  }
  return fc;
}
  /*  adjust intensity values (presumably) for limb darkening before mapping  */
int adjust_for_limb_dark (DRMS_Array_t *img, double img_xc, double img_yc,
    double rsun, double *ldcoef, int ldorder) {
  double rsun2, xp, xp2, xctr, yp, yp2, yctr;
  double mu, mu2;
  int n, col, cols, row, rows;

  float *v = (float *)img->data;
  int status = 0;
			      /*  no point in just multiplying by a constant  */
  if (ldorder < 1) return status;

  cols = img->axis[0];
  rows = img->axis[1];
  xctr = img_xc + 0.5 * cols - 0.5;
  yctr = img_yc + 0.5 * rows - 0.5;
  rsun2 = rsun * rsun;
  for (row = 0, n = 0; row < rows; row++) {
    yp = row - yctr;
    yp2 = yp * yp;
    for (col = 0; col < cols; col++, n++) {
      xp = col - xctr;
      xp2 = xp * xp;
      mu2 = 1.0 - (xp2 + yp2) / rsun2;
      if (mu2 < 0) continue;
      mu = sqrt (mu2);
      v[n] /= limbdark (mu, ldcoef, ldorder);
    }
  }
  return status;
}
				/*  Calculate the interpolation kernel.  */
void ccker (double *u, double s) {
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
  if (isnan (x) || isnan (y)) return missing_val;

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
						 /*  Bilinear interpolation  */
float linint2 (float *f, int cols, int rows, double x, double y) {
  double p, q, val;
  int col = (int)x, row = (int)y;
  int onerow = cols * row;
  int colp1 = col + 1, onerowp1 = onerow + cols;

  if (x < 0.0 || x > cols  || y < 0.0 || y >= rows)
    return missing_val;
  if (isnan (x) || isnan (y)) return missing_val;
  p = x - col;
  q = y - row;
  val = (1 - p) * (1 - q) * f[col + onerow]
      + p * (1 - q) * f[colp1 + onerow]
      + (1 - p) * q * f[col + onerowp1]
      + p * q * f[colp1 + onerowp1];
  return val;
}
					  /*  nearest value "interpolation"  */
float nearest_val (float *f, int cols, int rows, double x, double y) {
  int col, row;
  if (x < -0.5 || x >= (cols - 0.5) || y < -0.5 || y >= (rows - 0.5))
    return missing_val;
  if (isnan (x) || isnan (y)) return missing_val;
  col = x + 0.5;
  row = y + 0.5;
  return f[col + row * cols];
}

float array_imaginterp (DRMS_Array_t *img, double x, double y,
    int schema) {
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
  else if (schema == INTERP_BILINEAR)
    return linint2 (img->data, cols, rows, xs, ys);
  else return ccint2 (img->data, cols, rows, xs, ys);
}

void perform_mapping (DRMS_Array_t *img, float *map,
    double *maplat, double *maplon, double *map_coslat, double *map_sinlat,
    int pixct, unsigned char *offsun, double latc, double lonc,
    double xc, double yc, double radius, double pa, double ellipse_e,
    double ellipse_pa, int x_invrt, int y_invrt, int interpolator,
    int cvlostor, int MDI_correct_distort) {
/*
 *  Perform the mappings from the target heliographic coordinate sets
 *    appropriate to each output cube into the image coordinates (as
 *    corrected) for spatial interpolation of the data values
 */
  static double sin_asd = 0.004660, cos_asd = 0.99998914;
						   /*  appropriate to 1 AU  */
  double r, cos_cang, xr, yr;
  double cos_lat, sin_lat, lon, cos_lat_lon;
  double xx, yy;
  double cosrho, sinrho;
  float interpval;
  int n;

  double cos_pa = cos (pa);
  double sin_pa = sin (pa);
  double cos_latc = cos (latc);
  double sin_latc = sin (latc);
  int plate_cols = img->axis[0];
  int plate_rows = img->axis[1];
  double plate_width = (plate_cols > plate_rows) ? plate_cols : plate_rows;

  xc *= 2.0 / plate_width;
  yc *= 2.0 / plate_width;
  radius *= 2.0 / plate_width;

  for (n = 0; n < pixct; n++) {
       /*  Calculate heliographic coordinates corresponding to map location  */
    if (offsun[n]) {
      map[n] = missing_val;
      continue;
    }
    sin_lat = map_sinlat[n];
    cos_lat = map_coslat[n];
    lon = maplon[n];
    cos_lat_lon = cos_lat * cos (lon - lonc);
    cos_cang  = sin_lat * sin_latc + cos_latc * cos_lat_lon;
    if (cos_cang < 0.0) {
      map[n] = missing_val;
      continue;
    }
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
    interpval = array_imaginterp (img, xx, yy, interpolator);
    if (cvlostor) {
      sinrho = sqrt (xr * xr + yr * yr);
      cosrho = sqrt (1.0 - sinrho * sinrho);
      interpval /= cosrho;
    }
				  /*  Correction for MDI distortion and tip  */
       /*  should be replaced by call to MDI_correct_plateloc when verified  */
    if (MDI_correct_distort) {
      mtrack_MDI_image_tip (&xx, &yy, 1, 1);
      mtrack_MDI_image_stretch (&xx, &yy, 1, 1);
    }
    map[n] = (isnan (interpval)) ? missing_val : interpval;
  }
}

int near_grid_line (double lon, double lat, double grid, double near) {
/*
 *  Return 1 if a target point (lon, lat) is within (near) deg. of a grid line
 *    with spacing (grid) deg.
 */
  static double degrad = 180.0 / M_PI;
  double g2 = 0.5 * grid;

  lon *= degrad;
  lat *= degrad;

  while (lon < 0.0) lon += grid;
  while (lon > g2) lon -= grid;
  if (fabs (lon) < near) return 1;
  while (lat < 0.0) lat += grid;
  while (lat > g2) lat -= grid;
  if (fabs (lat) < near) return 1;
  return 0;
}

int parseseglist (char *list, const char *delim, char ***parsed_list) {
  int n, count = 0;
  int maxct = strlen (list);
  char **tmplist = (char **)malloc (maxct * sizeof (char *));
  char *entry;
  char *orig = strdup (list);
  char c;

  while (*orig != '{' && *orig != '\0') orig++;
  if (orig[0] == '{') orig++;
  else {
    free (tmplist);
    return count;
  }
  for (n = 0; n < strlen (orig); n++) if (orig[n] == '}') orig[n] = '\0';
  entry = strtok (orig, delim);
  while (entry) {
    tmplist[count++] = strdup (entry);
    entry = strtok (NULL, delim);
  }
  *parsed_list = (char **)malloc (count * sizeof (char **));
  for (n = 0; n < count; n++) {
    char *subs = (*parsed_list)[n] = strdup (tmplist[n]);
    free (tmplist[n]);
					     /*  remove leading white space  */
    while (isspace (c = *subs)) subs++;
    (*parsed_list)[n] = subs;
					    /*  remove trailing white space  */
    if (strlen (subs)) {
      subs += strlen (subs) - 1;
      while (isspace (c = *subs)) {
	*subs = '\0';
	subs--;
      }
    }
  }
  free (tmplist);
  return count;
}
/*
 *  Check for permission to write to output data series
 */
static int check_output_series_permission (char *series, int warnonly) {
  DB_Text_Result_t *qres;
  char *dbuser;
  char cmd[DRMS_MAXQUERYLEN];

  int status = 0;
						     /*  get current db role  */
#ifndef DRMS_CLIENT
  dbuser = drms_env->session->db_handle->dbuser;
#else
  drms_send_commandcode (drms_env->session->sockfd, DRMS_GETDBUSER);
  dbuser = receive_string(drms_env->session->sockfd);
#endif
					   /*  check for existence of series  */
  if (!drms_series_exists (drms_env, series, &status)) {
    if (warnonly) fprintf (stderr, "Warning: ");
    else {
      fprintf (stderr, "Error: ");
      status = 1;
    }
    fprintf (stderr, "output data series %s does not exist\n", series);
    return status;
  }
						    /*  check for permission  */
  sprintf (cmd, "SELECT has_table_privilege(\'%s\', \'insert\')", series);
  if ((qres = drms_query_txt (drms_env->session, cmd)) == NULL) {
    fprintf (stderr, "Error: can\'t connect to DRMS database\n");
    return 1;
  }
  if (strcmp (qres->field[0][0], "t")) {
    if (warnonly) fprintf (stderr, "Warning: ");
    else {
      fprintf (stderr, "Error: ");
      status = 1;
    }
    fprintf (stderr, "You do not have insert permission on output series %s\n",
	series);
  }
  return status;
}
/*
 *  Check output data series structure for the existence of segments of
 *    appropriate protocol for array segment writes
 */
static int check_output_record_structure (DRMS_Record_t *rec, int rank,
    int *axis, int *segnum, int warnonly) {
  DRMS_Segment_t *seg;
  int i, n, segct;

  int found = 0;

  segct = drms_record_numsegments (rec);
  for (n = 0; n < segct; n++) {
    seg = drms_segment_lookupnum (rec, n);
    if (seg->info->protocol == DRMS_PROTOCOL_INVALID ||
	  seg->info->protocol == DRMS_GENERIC) continue;
    if (seg->info->naxis != rank) continue;
    if (seg->info->scope != DRMS_VARDIM)
      for (i = 0; i < rank; i++) if (seg->axis[i] != axis[i]) continue;
    if (!found) *segnum = n;
    found++;
  }
  if (found < 1) {
    if (warnonly) fprintf (stderr, "Warning: ");
    else fprintf (stderr, "Error:   ");
    fprintf (stderr,
	"no segment of rank %d and appropriate size in output series\n",
	rank);
    fprintf (stderr, "         %s\n", rec->seriesinfo->seriesname);
  }
  if (found > 1) {
    fprintf (stderr,
	"Warning: multiple segments of rank %d and appropriate size in output series\n",
	rank);
    fprintf (stderr, "         %s; ", rec->seriesinfo->seriesname);
    seg = drms_segment_lookupnum (rec, *segnum);
    fprintf (stderr, "using \"%s\"\n", seg->info->name);
  }
  return found;
}

/******************************************************************************/
			/*  module body begins here  */
/******************************************************************************/
int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids;
  DRMS_Record_t *irec, *orec;
  DRMS_Segment_t *iseg, *oseg;
  DRMS_Array_t *image = NULL, *map = NULL;
  double *maplat, *maplon, *map_coslat, *map_sinlat;
  double *ldcoef;
  double x, y, x0, y0, xstp, ystp, xrot, yrot;
  double lat, lon, cos_phi, sin_phi;
  double img_lat, img_lon;
  double img_xscl, img_yscl, img_xc, img_yc, img_radius, img_pa;
  double grid_width;
  double ellipse_e, ellipse_pa;
  float *data;
  unsigned int quality;
  int axes[2];
  int img, imgct, pixct, segct;
  int osegnum;
  int found, kstat, status;
  int need_ephem, need_lat, need_lon;
  int x_invrt, y_invrt;
  int n, col, row;
  int MDI_correct, MDI_correct_distort;
  char **segnames;
  unsigned char *offsun, *ongrid;
  char *input, *isegname;
  char source[DRMS_MAXQUERYLEN], recid[DRMS_MAXQUERYLEN];
  char module_ident[64], key[64], tbuf[64];

  double raddeg = M_PI / 180.0;
  double degrad = 1.0 / raddeg;
  char *mapname[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
      "LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
      "stereographic", "orthographic", "LambertAzimuthal"};
  char *interpname[] = {"Cubic Convolution", "Nearest Neighbor", "Bilinear"};
  char *wcscode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
      "SIN", "ZEA"};
  missing_val = 0.0 / 0.0;
  float bblank = -1.0 / 0.0;
  float wblank = 1.0 / 0.0;
  int adjust_values = 0;
  int fix_limb_dark = 0;
						  /*  process command params  */
  char *inset = strdup (params_get_str (params, "in"));
  char *outser = strdup (params_get_str (params, "out"));
  double clat = params_get_double (params, "clat") * raddeg;
  double clon = params_get_double (params, "clon") * raddeg;
  double map_scale = params_get_double (params, "scale");
  double map_pa = params_get_double (params, "map_pa") * raddeg;
  float bscale = params_get_float (params, "bscale");
  float bzero = params_get_float (params, "bzero");
  float grid_spacing = params_get_float (params, "grid");
  float grid_ratio = params_get_float (params, "grid_ratio");
  unsigned int qmask = cmdparams_get_int64 (params, "qmask", &status);
  int map_cols = params_get_int (params, "cols");
  int map_rows = params_get_int (params, "rows");
  int proj = params_get_int (params, "map");
  int intrpopt = params_get_int (params, "interp");
  int ldcvals = params_get_int (params, "ldcoef_nvals");
  char *qual_key = strdup (params_get_str (params, "qual_key"));
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *rsun_key = strdup (params_get_str (params, "rsun_key"));
  char *apsd_key = strdup (params_get_str (params, "apsd_key"));
  char *dsun_key = strdup (params_get_str (params, "dsun_key"));
  char *RequestID = strdup (params_get_str (params, "RequestID"));
  int center = params_isflagset (params, "c");
  int no_save = params_isflagset (params, "n");
  int stonyhurst = params_isflagset (params, "s");
  int verbose = params_isflagset (params, "v");
  int overlay = (isfinite (grid_spacing));
  int MDI_proc = params_isflagset (params, "M");
  int cvlostor = params_isflagset (params, "R");

  int dispose = (no_save) ? DRMS_FREE_RECORD : DRMS_INSERT_RECORD;
  int ldorder = ldcvals - 1;

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s: JSOC version %s\n", module_ident, jsoc_version);
						/*  check calling parameters  */
  if (map_cols < 1) map_cols = map_rows;
  if (map_rows < 1) map_rows = map_cols;
  if (map_rows < 1) {
    fprintf (stderr, "Error: at least one of \"cols\" or \"rows\" must be set\n");
    return 1;
  }
  if (isnan (map_scale) || map_scale == 0.0) {
    fprintf (stderr,
	"Error: auto scaling from image resolution not implemented;\n");
    fprintf (stderr, "       scale parameter must be set.\n");
    return 1;
  }
		     /*  optional parameters for preprocessing of image data  */
  if (ldorder) {
    adjust_values = fix_limb_dark = 1;
    ldcoef = (double *)malloc (ldcvals * sizeof (double));
    for (n = 0; n < ldcvals; n++) {
      sprintf (key, "ldcoef_%d_value", n);
      ldcoef[n] = params_get_double (params, key);
    }
  }

  need_lat = center || isnan (clat);
  need_lon = center || isnan (clon);
  MDI_correct = MDI_correct_distort = MDI_proc;
  cos_phi = cos (map_pa);
  sin_phi = sin (map_pa);
  xstp = ystp = map_scale * raddeg;
  x0 = 0.5 * (1.0 - map_cols) * xstp;
  y0 = 0.5 * (1.0 - map_rows) * ystp;
  grid_width = grid_ratio * grid_spacing;
							     /*  check input  */
  if (!(ids = drms_open_records (drms_env, inset, &status))) {
    fprintf (stderr, "Error: (%s) unable to open input data set \"%s\"\n",
        module_ident, inset);
    fprintf (stderr, "       status = %d\n", status);
    return 1;
  }
  if ((imgct = ids->n) < 1) {
    fprintf (stderr, "Error: (%s) no records in selected input set\n",
	module_ident);
    fprintf (stderr, "       %s\n", inset);
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
  input = strdup (inset);
	   /*  determine appropriate input record segment (if more than one)  */
  irec = ids->records[0];
  segct = drms_record_numsegments (irec);
  if (segct) {
    int namect = parseseglist (inset, ",", &segnames);
    found = 0;
    if (namect) {
      for (n = 0; n < namect; n++) {
	iseg = drms_segment_lookup (irec, segnames[n]);
	if (iseg->info->naxis != 2) continue;
        if (!found) isegname = segnames[n];
	found++;
      }
    } else {
      for (n = 0; n < segct; n++) {
        iseg = drms_segment_lookupnum (irec, n);
        if (iseg->info->naxis != 2) continue;
        if (!found) isegname = strdup (iseg->info->name);
        found++;
      }
    }
    if (found > 1) {
      fprintf (stderr,
	  "Warning: multiple data segments of dimension 2 in input series %s\n",
	  input);
      fprintf (stderr, "       using \"%s\"\n", isegname);
    }
  } else {
    fprintf (stderr, "Error: no data segments in input series %s\n", input);
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
  if (found < 1 || iseg->info->naxis != 2) {
    fprintf (stderr,
	"Error: no data segment of dimension 2 in input series %s\n", input);
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
						     /*  check output series  */
  if (check_output_series_permission (outser, no_save)) return 1;
  orec = drms_template_record (drms_env, outser, &status);
  if (status || !orec) {
    if (no_save) fprintf (stderr, "Warning: ");
    else {
      drms_close_records (ids, DRMS_FREE_RECORD);
      fprintf (stderr, "Error: ");
    }
    fprintf (stderr,
	" unable to verify appropriate structure of output series %s\n",
	outser);
    if (!no_save) return 1;
  }
  axes[0] = map_cols;
  axes[1] = map_rows;
  if (check_output_record_structure (orec, 2, axes, &osegnum, no_save) < 1)
    return 1;
						 /*  create output map array  */
  map = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, NULL, &status);
  if (status) {
    fprintf (stderr, "Error: couldn't create output array\n");
    return 1;
  }
  pixct = map_cols * map_rows;
  maplat = (double *)malloc (pixct * sizeof (double));
  maplon = (double *)malloc (pixct * sizeof (double));
  map_coslat = (double *)malloc (pixct * sizeof (double));
  map_sinlat = (double *)malloc (pixct * sizeof (double));
  offsun = (unsigned char *)malloc (pixct * sizeof (char));
  if (overlay) ongrid = (unsigned char *)malloc (pixct * sizeof (char));
	     /*  use output series default segment scaling if not overridden  */
  oseg = drms_segment_lookupnum (orec, osegnum);
  if (bscale == 0.0) {
    bscale = oseg->bscale;
    if (verbose) printf ("bscale set to output default: %g\n", bscale);
  }
  if (isnan (bzero)) {
    bzero = oseg->bzero;
    if (verbose) printf ("bzero set to output default: %g\n", bzero);
  }
  map->bscale = bscale;
  map->bzero = bzero;
  if (need_lat || need_lon) {
    if (need_lon) {
      img_lon = drms_getkey_double (irec, clon_key, &status);
      if (status || isnan (img_lon)) fprintf (stderr,
	  "Error: no valid value for key %s in first input record\n", clon_key);
      else if (!center) need_lon = 0;
      clon = img_lon * raddeg;
    }
    if (need_lat) {
      img_lat = drms_getkey_double (irec, clat_key, &status);
      if (status || isnan (img_lat)) fprintf (stderr,
	  "Error: no valid value for key %s in first input record\n", clat_key);
      else if (!center) need_lat = 0;
      clat = img_lat * raddeg;
    }
    if (verbose) printf ("map(s) centred at latitude %+.2f, longitude %.2f\n",
	clat / raddeg, clon/raddeg);
  }
     /*  Calculate heliographic coordinates corresponding to map location(s)  */
  for (n=0, row=0, y=y0; row < map_rows; row++, y += ystp) {
    for (col=0, x=x0; col < map_cols; col++, x += xstp, n++) {
      xrot = x * cos_phi - y * sin_phi;
      yrot = y * cos_phi + x * sin_phi;
      offsun[n] = plane2sphere (xrot, yrot, clat, clon, &lat, &lon, proj);
      maplat[n] = lat;
      maplon[n] = lon;
      map_coslat[n] = cos (lat);
      map_sinlat[n] = sin (lat);
      if (overlay) ongrid[n] =
	  near_grid_line (maplon[n], maplat[n], grid_spacing, grid_width);
    }
  }

  data = (float *)map->data;
					 /*  process individual input images  */
  for (img = 0; img < imgct; img++) {
    irec = ids->records[img];
			  /*  check for record quality, reject as applicable  */
    quality = drms_getkey_int (irec, qual_key, &status);
    if ((quality & qmask) && !status) {
      if (verbose) printf ("record %d skipped: quality = %08x\n", img, quality);
      continue;
    }
 							 /*  get input image  */
    drms_sprint_rec_query (source, irec);
    iseg = drms_segment_lookup (irec, isegname);
    image = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
						/*  trap failed segment read  */
    if (status || !image) {
      fprintf (stderr, "Error reading segment %s from record\n      %s\n",
	  iseg->info->name, source);
      if (image) drms_free_array (image);
      continue;
    }
			    /*  get needed info from record keys for mapping  */
			      /*  replace with call to solar_ephemeris_info?  */
    img_lon = drms_getkey_double (irec, clon_key, &status);
    img_lat = drms_getkey_double (irec, clat_key, &status);
    if (img && (need_lon || need_lat)) {
      if (need_lon) clon = img_lon * raddeg;
      if (need_lat) clat = img_lat * raddeg;
      if (verbose) printf ("map centred at latitude %+.2f, longitude %.2f\n",
	clat / raddeg, clon/raddeg);
   /*  Recalculate heliographic coordinates corresponding to map location(s)  */
      for (n=0, row=0, y=y0; row < map_rows; row++, y += ystp) {
	for (col=0, x=x0; col < map_cols; col++, x += xstp, n++) {
	  xrot = x * cos_phi - y * sin_phi;
	  yrot = y * cos_phi + x * sin_phi;
	  offsun[n] = plane2sphere (xrot, yrot, clat, clon, &lat, &lon, proj);
	  maplat[n] = lat;
	  maplon[n] = lon;
	  map_coslat[n] = cos (lat);
	  map_sinlat[n] = sin (lat);
	  if (overlay) ongrid[n] =
	      near_grid_line (maplon[n], maplat[n], grid_spacing, grid_width);
	}
      }
    }

    status = solar_image_info (irec, &img_xscl, &img_yscl, &img_xc, &img_yc,
	&img_radius, rsun_key, apsd_key, &img_pa, &ellipse_e, &ellipse_pa,
	&x_invrt, &y_invrt, &need_ephem, 0);
    if (status & NO_SEMIDIAMETER) {
      int keystat = 0;
      double dsun_obs = drms_getkey_double (irec, dsun_key, &keystat);
      if (keystat) {
	fprintf (stderr, "Error: one or more essential keywords or values missing; skipped\n");
	fprintf (stderr, "solar_image_info() returned %08x :\n", status);
	if (status & NO_DATA_DICT)
	  fprintf (stderr, "  internal data dictionary for platform\n");
	if (status & NO_SEMIDIAMETER) fprintf (stderr, "  semidiameter\n");
	if (status & NO_XSCALE) fprintf (stderr, "  x-scale\n");
	if (status & NO_YSCALE) fprintf (stderr, "  y-scale\n");
	if (status & NO_XCENTERLOC) fprintf (stderr, "  x-center\n");
	if (status & NO_YCENTERLOC) fprintf (stderr, "  y-center\n");
	if (status & NO_HELIO_LATC) fprintf (stderr, "  disc-center latitude\n");
	if (status & NO_HELIO_LONC) fprintf (stderr, "  disc-center longitude\n");
	if (status & NO_XUNITS) fprintf (stderr, "  x units\n");
	if (status & NO_YUNITS) fprintf (stderr, "  y units\n");
	continue;
      }
			       /*  set image radius from scale and distance  */
      img_radius = asin (RSUNM / dsun_obs);
      img_radius *= 3600.0 * degrad;
      img_radius /= (img_xscl <= img_yscl) ? img_xscl : img_yscl;
      status &= ~NO_SEMIDIAMETER;
    }
    if (status == KEYSCOPE_VARIABLE) {
      fprintf (stderr, "Warning: one or more keywords expected constant are variable\n");
    } else if (status) {
      fprintf (stderr, "Error: one or more essential keywords or values missing; skipped\n");
      fprintf (stderr, "solar_image_info() returned %08x :\n", status);
      if (status & NO_DATA_DICT)
	fprintf (stderr, "  internal data dictionary for platform\n");
      if (status & NO_SEMIDIAMETER) fprintf (stderr, "  semidiameter\n");
      if (status & NO_XSCALE) fprintf (stderr, "  x-scale\n");
      if (status & NO_YSCALE) fprintf (stderr, "  y-scale\n");
      if (status & NO_XCENTERLOC) fprintf (stderr, "  x-center\n");
      if (status & NO_YCENTERLOC) fprintf (stderr, "  y-center\n");
      if (status & NO_HELIO_LATC) fprintf (stderr, "  disc-center latitude\n");
      if (status & NO_HELIO_LONC) fprintf (stderr, "  disc-center longitude\n");
      if (status & NO_XUNITS) fprintf (stderr, "  x units\n");
      if (status & NO_YUNITS) fprintf (stderr, "  y units\n");
     continue;
    }
    if (MDI_correct) {
      mtrack_MDI_correct_imgctr (&img_xc, &img_yc, img_radius);
      mtrack_MDI_correct_pa (&img_pa);
    }
    img_xc -= 0.5 * (image->axis[0] - 1);
    img_yc -= 0.5 * (image->axis[1] - 1);
		 	/*  should be taken care of in solar_ephemeris_info  */
    img_lon *= raddeg;
    img_lat *= raddeg;
    if (adjust_values) {
      if (fix_limb_dark) {
        adjust_for_limb_dark (image, img_xc, img_yc, img_radius, ldcoef,
	    ldorder);
      }
    }
						    /*  set up output record  */
    orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    oseg = drms_segment_lookupnum (orec, osegnum);
 						     /*  perform the mapping  */
    perform_mapping (image, data, maplat, maplon, map_coslat, map_sinlat,
        pixct, offsun, img_lat, img_lon, img_xc, img_yc, img_radius, img_pa,
	ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt, cvlostor,
	MDI_correct_distort);
    if (overlay) {
      for (n = 0; n < pixct; n++) 
        if (ongrid[n]) data[n] = (isfinite (data[n])) ? bblank : wblank;
    }
				/*  write map array to output record segment  */
    status = drms_segment_write (oseg, map, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      return 1;
    }
    drms_copykeys (orec, irec, 0, kDRMS_KeyClass_Explicit);
					    /*  set output record key values  */
    kstat = 0;
    kstat += check_and_set_key_str   (orec, "WCSNAME", "Carrington Heliographic");
    kstat += check_and_set_key_int   (orec, "WCSAXES", 2);
    snprintf (key, 64, "CRLN-%s", wcscode[proj]);
    kstat += check_and_set_key_str   (orec, "CTYPE1", key);
    snprintf (key, 64, "CRLT-%s", wcscode[proj]);
    kstat += check_and_set_key_str   (orec, "CTYPE2", key);
    kstat += check_and_set_key_str   (orec, "CUNIT1", "deg");
    kstat += check_and_set_key_str   (orec, "CUNIT2", "deg");
    kstat += check_and_set_key_float (orec, "CRPIX1", 0.5 * map_cols + 0.5);
    kstat += check_and_set_key_float (orec, "CRPIX2", 0.5 * map_rows + 0.5);
    kstat += check_and_set_key_float (orec, "CRVAL1", clon * degrad);
    kstat += check_and_set_key_float (orec, "CRVAL2", clat * degrad);
    kstat += check_and_set_key_float (orec, "CDELT1", map_scale);
    kstat += check_and_set_key_float (orec, "CDELT2", map_scale);
    if (map_pa != 0.0) {
      kstat += check_and_set_key_float (orec, "PC1_1", cos (map_pa));
/*  PC1_2 should be multiplied by CDELT2 / CDELT1  */
      kstat += check_and_set_key_float (orec, "PC1_2", sin (map_pa));
/*  PC2_1 should be multiplied by CDELT1 / CDELT2  */
      kstat += check_and_set_key_float (orec, "PC2_1", sin (map_pa));
      kstat += check_and_set_key_float (orec, "PC2_2", cos (map_pa));
    }

    kstat += check_and_set_key_float (orec, "LonHG", clon * degrad);
    kstat += check_and_set_key_float (orec, "LatHG", clat * degrad);
    kstat += check_and_set_key_str   (orec, "MapProj", mapname[proj]);
    kstat += check_and_set_key_float (orec, "MapScale", map_scale);
    kstat += check_and_set_key_float (orec, "Width", map_cols * map_scale);
    kstat += check_and_set_key_float (orec, "Height", map_rows * map_scale);
    kstat += check_and_set_key_float (orec, "Size", sqrt (map_rows * map_cols) * map_scale);
    kstat += check_and_set_key_float (orec, "Map_PA", map_pa / raddeg);
    kstat += check_and_set_key_float (orec, "RSunRef", 1.0e-6 * RSUNM);
    kstat += check_and_set_key_str   (orec, "Interp", interpname[intrpopt]);
    kstat += check_and_set_key_short (orec, "MDICorr", MDI_proc);
    kstat += check_and_set_key_short (orec, "LoS2Rad", cvlostor);

    kstat += check_and_set_key_str   (orec, "Module", module_ident);
    kstat += check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
    kstat += check_and_set_key_time  (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_time  (orec, "DATE", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str   (orec, "Source", source);
    kstat += check_and_set_key_str   (orec, "Input", input);
    kstat += check_and_set_key_str   (orec, "RequestID", RequestID);

    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record %d in series %s\n",
	  img, outser);
      fprintf (stderr, "      series may not have appropriate structure\n");
      drms_close_record (orec, DRMS_FREE_RECORD);
      return 1;
    }
    drms_close_record (orec, dispose);
    drms_free_array (image);
  }
  drms_free_array (map);
  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  2008.04.22 	created this file, based on fastrack
 *  10.10.04	updated output keyword set and fleshed out arguments, commented
 *		added setting of WCS PC matrix elements as necessary; added
 *		input data processing, mapping with interpolation options
 *  v 0.8	frozen 2010.10.05
 *  11.05.02	added option for grid overlay; removed inappropriate defaults
 *  12.01.17	added option for removal of MDI distortion
 *  12.03.07	added promiscuous copy to target keys before careful setting
 *  v 0.9 frozen 2012.08.03
 *  13.02.01	removed unused variable
 *  13.04.23	added detail for failures when geometry keywords missing
 *  14.03.27	added option for copying of selected keywords from input to
 *		output
 *  14.08.01	made clat and clon optional arguments, with default being
 *		centering at image center coordinates (or just the unspecified
 *		coordinate) of first image; implemented option for centering
 *		maps at center of each image; added setting of DATE (same as
 *		Created); added option for converting scalar values assuming
 *		that they are the line-of-sight components of radial vectors;
 *		Added optional string argument RequestID
 *  14.11.17	changed copykeys to avoid copying implicit indexing keys
 *		added setting of optional keywords LoS2Rad and MDICorr depending
 *		on flag values
 *  v 1.0 frozen 2015.01.13
 *  15.07.27	free each input image array
 *  v 1.1 frozen 2015.09.21
 *  16.02.09	trap failed segment reads
 *  16.02.26	added nomap option -n for testing
 *  v 1.2 frozen 2016.06.20
 *  16.09.07	added function for parsing optional segment list in input
 *		specification; changed to use segment names rather than internal
 *		numbers which get mangled by links
 *  18.07.30	added optional argument grid_ratio to specify ratio of overlay
 *		grid line width to spacing, needed for small spacings
 *  18.09.14	added option and related arguments for checking quality key
 *		value against bit mask; improved checking of acceptability of
 *		output series; some version control of include files
 *  v 1.3 frozen 2019.02.11
 *  19.08.09	added option for removal of limb darkening before mapping
 *  20.07.13	fixed bug occurring when no appropriate input segments found
 *  20.12.08	added trap for NaN values in map target locations which can be
 *	  	returned in rare circumstances by version 10 of plane2sphere()
 *	  	without setting the offsun flag to the interpolation routines
 *		use v 11 of cartographic routines
 *  v 1.4 frozen 2021.02.15
 *
 */
