/*
 *  maproj.c						~rick/src/maproj
 *
 *  Responsible:  Rick Bogart				RBogart@spd.aas.org
 *
 *  Map a set of input solar images into a set of output maps
 *
 *  Parameters:	(type	default		description)
 *	in	DataSet TBD		Input dataset
 *				A set of images of all or part of the
 *				solar disc in a "plate" coordinate system
 *				(helioprojective geometry).
 *      out     DataSer TBD		Output data series name
 *	clat	double	0.0		Map central heliographic latitude
 *	clon	double	0.0		Map central heliographic longitude
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
 *				they are within 0.01 * the grid spacing from it
 *      map_pa  float   0.0             The angle between heliographic north
 *				and "up" on the output map (in the direction
 *				of increasing rows) [deg[, in the sense that a
 *				positive position angle represents a clockwise
 *				displacement of the north axis.
 *	bscale	float	0.0		Value scaling parameter for output
 *	bzero	float	NaN		Value offset parameter for output
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
 *	-c	center map9s) at image center(s)
 *	-s	interpret clon as Stonyhurst rather than Carrington longitude
 *	-v	run verbose
 *	-M	correct for MDI distortion
 *
 *  Bugs:
 *    Basic functionality is present, but with several fixed and inappropriate
 *	defaults and some missing arguments
 *    No provision for propagation of default or selected keywords
 *    Uses considerable replicated code from mtrack, esp. array_imaginterp()
 *	function and its dependencies; should be consolidated
 *    Values for Input and Source are inappropriate, refer to whole input
 *	data set
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

#include <jsoc_main.h>
#include "fstats.h"  // XXXX PHS
#include "/home/jsoc/cvs/Development/JSOC/proj/rings/apps/cartography.c" // XXX PHS path change
#include "/home/jsoc/cvs/Development/JSOC/proj/rings/apps/imginfo.c" // XXX PHS path change
#include "/home/jsoc/cvs/Development/JSOC/proj/rings/apps/mdistuff.c" // XXX PHS path change

#define RECTANGULAR	(0)
#define CASSINI		(1)
#define MERCATOR	(2)
#define CYLEQA		(3)
#define SINEQA		(4)
#define GNOMONIC	(5)
#define POSTEL		(6)
#define STEREOGRAPHIC	(7)
#define ORTHOGRAPHIC	(8)
#define LAMBERT		(9)

#define RSUNM		(6.96e8)
#define INTERP_NEAREST_NEIGHBOR	(1)
#define INTERP_BILINEAR	(2)
						       /*  module identifier  */
char *module_name = "maproj";
char *module_desc = "mapping from solar images";
char *version_id = "0.9";

ModuleArgs_t module_args[] = {
  {ARG_DATASET,	"in", "", "input data set"}, 
  {ARG_DATASERIES, "out", "", "output data series"}, 
  {ARG_DOUBLE,	"clat", "0.0", "heliographic latitude of map center [deg]"},
  {ARG_DOUBLE,	"clon", "0.0", "Carrington longitude of map center [deg]"},
  {ARG_DOUBLE,	"scale", "Not specified", "map scale at center [deg/pxl]"},
  {ARG_NUME, "map", "orthographic", "map projection",
      "carree, Cassini, Mercator, cyleqa, sineqa, gnomonic, Postel, stereographic, orthographic, Lambert"},
  {ARG_NUME, "interp", "cubiconv", "interpolation option",
      "cubiconv, nearest, bilinear"},
  {ARG_FLOAT, "grid", "Not Specified",
      "if specified, spacing of grid overlay [deg]"},
  {ARG_INT,     "cols", "0", "columns in output map"},
  {ARG_INT,     "rows", "0", "rows in output map"},
  {ARG_FLOAT,	"map_pa", "0.0", "position angle of north on output map [deg]"},
  {ARG_FLOAT,	"bscale", "0.0", "output scaling factor"},
  {ARG_FLOAT,	"bzero", "Default", "output offset"},
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"rsun_key", "R_SUN", "keyname for image semi-diameter (pixel)"}, 
  {ARG_STRING,	"apsd_key", "RSUN_OBS", "keyname for apparent solar semi-diameter (arcsec)"}, 
  {ARG_STRING,	"dsun_key", "DSUN_OBS", "keyname for observer distance"}, 
  {ARG_STRING,	"requestid", "none", "RequestID for jsoc export management"},  // XXX PHS
  {ARG_FLAG, "A", NULL, "Generate output for all input segments"},
  {ARG_FLAG,	"c",	"", "center map at center of image"}, 
  {ARG_FLAG,	"s",	"", "clon is Stonyhurst longitude"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {ARG_FLAG,	"x",	"", "write full headers, for export purposes."}, // XXX PHS
  {}
};

		 /*  global declaration of missing to be initialized as NaN  */
float missing_val;

	/*  Calculate the interpolation kernel.  */
	/* Robert G. Keys, "Cubic Convolution Interpolation for digital Image Processing", IEEE
         * Transactions on Acoustics, Speech, and Signal Processing, Vol ASSP-29, No. 6, December 1981. 
         */
void ccker (double *u, double s) {
  double s2, s3;

  s2= s * s;
  s3= s2 * s;
  u[0] = s2 - 0.5 * (s3 + s);
  u[1] = 1.5*s3 - 2.5*s2 + 1.0;
  u[2] = -1.5*s3 + 2.0*s2 + 0.5*s;
  u[3] = 0.5 * (s3 - s2);
}

        /* note that the end points do not need to be discarded using:
         * equiv of for f[0] = 3*data[0] - 3*data[1] + data[2]; 
         * and      for f[n-1] = 3*data[n-1] - 3*data[n-2] + data[n-3]; 
         * for target points from x=0 to x=1, and x=n-2 to x=n-2.
         * see equation after eqn 25 in Keys' paper 
         */

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
						 /*  Bilinear interpolation  */
float linint2 (float *f, int cols, int rows, double x, double y) {
  double p, q, val;
  int col = (int)x, row = (int)y;
  int onerow = cols * row;
  int colp1 = col + 1, onerowp1 = onerow + cols;

  if (x < 0.0 || x > cols  || y < 0.0 || y >= rows)
    return missing_val;
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
    int MDI_correct_distort) {
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

/* These were copied from keystuff.c. I made a number of changes - the functions 
 * used to return success if the key did not exist (not sure why). They also 
 * used to follow links and attempt to set linked keywords, which was not correct. 
 * And they did not check the return value of the lower-level DRMS set-key functions. */
static int check_and_set_key_str(DRMS_Record_t *outRec, const char *key, const char *val) 
{
    DRMS_Keyword_t *keywd = NULL;
    char *vreq;
    int status;
    int errOut;

    errOut = 0;

    if (!(keywd = drms_keyword_lookup(outRec, key, 1))) return 1;

    if (!drms_keyword_isconstant(keywd)) 
    {
                   /*  it's not a constant, so don't worry, just set it  */
        return (drms_setkey_string(outRec, key, val) != DRMS_SUCCESS);
    }

    vreq = drms_getkey_string(outRec, key, &status);
    if (status || !vreq) 
    {
        fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
        errOut = 1;
    }
    
    if (strcmp(vreq, val)) 
    {
        char format[256];
        
        snprintf (format, sizeof(format), "Error:  new value \"%s\" for constant keyword %%s\n differs from the existing value \"%s\"\n", keywd->info->format, keywd->info->format);
        fprintf(stderr, format, val, key, vreq);
        errOut = 1;
    }
    
    if (vreq)
    {
        /* drms_getkey_string() allocates mem. */
        free(vreq);
        vreq = NULL;
    }
    
    if (errOut)
    {
        return 1;
    }
    
    return 0;
}

static int check_and_set_key_time(DRMS_Record_t *outRec, const char *key, TIME tval) 
{
    DRMS_Keyword_t *keywd = NULL;
    TIME treq;
    int status;
    char sreq[64], sval[64];

    if (!(keywd = drms_keyword_lookup(outRec, key, 1))) return 1;

    if (!drms_keyword_isconstant(keywd)) 
    {
                   /*  it's not a constant, so don't worry, just set it  */
        return (drms_setkey_time(outRec, key, tval) != DRMS_SUCCESS);
    }
    
    treq = drms_getkey_time(outRec, key, &status);
    if (status) 
    {
        fprintf(stderr, "Error retrieving value for constant keyword %s\n", key);
        return 1;
    }
    
    sprint_time(sreq, treq, keywd->info->unit, atoi(keywd->info->format));
    sprint_time(sval, tval, keywd->info->unit, atoi(keywd->info->format));
    
    if (strcmp(sval, sreq)) 
    {
        fprintf (stderr, "Error:  value %s for constant keyword %s\n", sval, key);
        fprintf (stderr, "  differs from existing value %s\n", sreq);
        return 1;
    }

    return 0;
}

static int check_and_set_key_int(DRMS_Record_t *outRec, const char *key, int val) 
{
    DRMS_Keyword_t *keywd = NULL;
    int vreq;
    int status;

    if (!(keywd = drms_keyword_lookup(outRec, key, 1))) return 1;

    if (!drms_keyword_isconstant(keywd))
    {
                   /*  it's not a constant, so don't worry, just set it  */
        return (drms_setkey_int(outRec, key, val) != DRMS_SUCCESS);
    }

    vreq = drms_getkey_int(outRec, key, &status);
    if (status) 
    {
        fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
        return 1;
    }
    
    if (vreq != val) 
    {
        char format[256];
        
        snprintf (format, sizeof(format), "Error:  new value \"%s\" for constant keyword %%s\n differs from the existing value \"%s\"\n", keywd->info->format, keywd->info->format);
        fprintf(stderr, format, val, key, vreq);
        return 1;
    }
    
    return 0;
}

static int check_and_set_key_float(DRMS_Record_t *outRec, const char *key, float val) 
{
    DRMS_Keyword_t *keywd = NULL;
    float vreq;
    int status;
    char sreq[128], sval[128];

    if (!(keywd = drms_keyword_lookup(outRec, key, 1))) return 1;

    if (!drms_keyword_isconstant(keywd)) 
    {
               /*  it's not a constant, so don't worry, just set it  */
        return (drms_setkey_float(outRec, key, val) != DRMS_SUCCESS);
    }
    
    vreq = drms_getkey_float(outRec, key, &status);
    
    if (status) 
    {
        fprintf (stderr, "Error retrieving value for constant keyword %s\n", key);
        return 1;
    }
    
    sprintf(sreq, keywd->info->format, vreq);
    sprintf(sval, keywd->info->format, val);
    
    if (strcmp(sreq, sval)) 
    {
        char format[256];
        
        snprintf(format, sizeof(format), "Error:  parameter value %s for keyword %%s\n  differs from required value %s\n", keywd->info->format, keywd->info->format);
        fprintf(stderr, format, val, key, vreq);
        return 1;
    }

    return 0;
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids, *ods;
  DRMS_Record_t *inrec, *orec;
  DRMS_Segment_t *inseg = NULL;
  DRMS_Segment_t *oseg = NULL;
  DRMS_Array_t *image = NULL, *map = NULL;
  int irec;
  double *maplat, *maplon, *map_coslat, *map_sinlat;
  double x, y, x0, y0, xstp, ystp, xrot, yrot;
  double lat, lon, cos_phi, sin_phi;
  double img_lat, img_lon;
  double img_xscl, img_yscl, img_xc, img_yc, img_radius, img_pa;
  double grid_width;
  double ellipse_e, ellipse_pa;
  float *data;
  int axes[2];
  int img, pixct, segct;
  int recCount;
  int found, kstat, status;
  int need_ephem;
  int x_invrt, y_invrt;
  int n, col, row;
  int MDI_correct, MDI_correct_distort;
  unsigned char *offsun, *ongrid;
  char *input;
  char source[DRMS_MAXQUERYLEN], recid[DRMS_MAXQUERYLEN];
  char module_ident[64], key[64], tbuf[64];

  double raddeg = M_PI / 180.0;
  double degrad = 1.0 / raddeg;
  int scaling_override = 0;
  char *mapname[] = {"PlateCarree", "Cassini-Soldner", "Mercator",
      "LambertCylindrical", "Sanson-Flamsteed", "gnomonic", "Postel",
      "stereographic", "orthographic", "LambertAzimuthal"};
  char *interpname[] = {"Cubic Convolution", "Nearest Neighbor", "Bilinear"};
  char *wcscode[] = {"CAR", "CAS", "MER", "CEA", "GLS", "TAN", "ARC", "STG",
      "SIN", "ZEA"};
  missing_val = 0.0 / 0.0;
  float bblank = -1.0 / 0.0;
  float wblank = 1.0 / 0.0;
						  /*  process command params  */
  const char *inset = params_get_str (params, "in");
  char *outser = strdup (params_get_str (params, "out"));
  double clat = params_get_double (params, "clat") * raddeg;
  double clon = params_get_double (params, "clon") * raddeg;
  double map_scale = params_get_double (params, "scale");
  double map_pa = params_get_double (params, "map_pa") * raddeg;
  float bscaleIn = params_get_float (params, "bscale");
  float bzeroIn = params_get_float (params, "bzero");
  float grid_spacing = params_get_float (params, "grid");
  int map_cols = params_get_int (params, "cols");
  int map_rows = params_get_int (params, "rows");
  int proj = params_get_int (params, "map");
  int intrpopt = params_get_int (params, "interp");
  char *RequestID = strdup (params_get_str (params, "requestid")); // XXX PHS
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *rsun_key = strdup (params_get_str (params, "rsun_key"));
  char *apsd_key = strdup (params_get_str (params, "apsd_key"));
  char *dsun_key = strdup (params_get_str (params, "dsun_key"));
  int center = params_isflagset (params, "c");
  int stonyhurst = params_isflagset (params, "s");
  int verbose = params_isflagset (params, "v");
// XXX PHS added -x
  int want_headers = params_isflagset (params, "x");
  int overlay = (isfinite (grid_spacing));
  int MDI_proc = params_isflagset (params, "M");
  
    int hasSegList = 0;
    int doingAllSegs = params_isflagset(params, "A");
    float bscale;
    float bzero;

  snprintf (module_ident, sizeof(module_ident), "%s v %s", module_name, version_id);
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
  MDI_correct = MDI_correct_distort = MDI_proc;
  cos_phi = cos (map_pa);
  sin_phi = sin (map_pa);
  xstp = ystp = map_scale * raddeg;
  x0 = 0.5 * (1.0 - map_cols) * xstp;
  y0 = 0.5 * (1.0 - map_rows) * ystp;
  grid_width = 0.01 * grid_spacing;
  
    char *testSegList = NULL;
    char *pCh = NULL;

    testSegList = rindex(inset, '}');
    if (testSegList)
    {
        for (pCh = testSegList + 1, hasSegList = 1; *pCh; pCh++)
        {
            if (!isspace(*pCh))
            {
                /* There is a non-ws char after '}' - this is not a valid seglist. drms_open_records() will fail below. */
                hasSegList = 0;
                break;
            }
        }
    }

							     /*  check input  */
    if (!(ids = drms_open_records(drms_env, inset, &status))) 
    {
        fprintf (stderr, "Error: (%s) unable to open input data set \"%s\"\n", module_ident, inset);
        fprintf (stderr, "       status = %d\n", status);
        return 1;
    }
  
    recCount = ids->n;

    if (recCount < 1) 
    {
        fprintf(stderr, "Error: (%s) no records in selected input set\n", module_ident);
        fprintf(stderr, "       %s\n", inset);
        drms_close_records(ids, DRMS_FREE_RECORD);
        return 1;
    }

						  /*  open output record set  */
  if (verbose) printf ("creating %d records in series %s\n", recCount, outser);
  ods = drms_create_records(drms_env, recCount, outser, DRMS_PERMANENT, &status);
  if (!ods || status != DRMS_SUCCESS) 
  {
    fprintf (stderr, "Error: unable to create %d records in series %s\n", recCount, outser);
    fprintf (stderr, "       drms_create_records() returned status %d\n", status); 
    return 1;
  }

          /* determine if output series is <input series>_mod as used for processing on export
             from exportdata.html */
  char *runderscore = rindex(outser, '_');
  int running_as_export = (runderscore && strcmp(runderscore, "_mod")==0);

						 /*  create output map array  */
  axes[0] = map_cols;
  axes[1] = map_rows;
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
	     

     /*  Calculate heliographic coordinates corresponding to map location(s)  */
    for (n=0, row=0, y=y0; row < map_rows; row++, y += ystp) 
    {
        for (col=0, x=x0; col < map_cols; col++, x += xstp, n++) 
        {
            xrot = x * cos_phi - y * sin_phi;
            yrot = y * cos_phi + x * sin_phi;
            offsun[n] = plane2sphere (xrot, yrot, clat, clon, &lat, &lon, proj);
            maplat[n] = lat;
            maplon[n] = lon;
            map_coslat[n] = cos (lat);
            map_sinlat[n] = sin (lat);
            if (overlay) ongrid[n] = near_grid_line (maplon[n], maplat[n], grid_spacing, grid_width);
        }
    }

    data = (float *)map->data;




					  /*  process individual input mages  */
					  
    /* This is really a loop over records (wherein it is assumed that there is a single image 
     * per record). */
  for (irec = 0; irec < recCount; irec++) 
  {
							 /*  get input image  */
    inrec = ids->records[irec];
    drms_sprint_rec_query(source, inrec);
    
    if (drms_record_numsegments(inrec) < 1)
    {
        fprintf(stderr, "Error: no data segments in input record-set %s\n", inset);
        drms_close_records(ids, DRMS_FREE_RECORD);
        return 1;
    }
    
    /* Segment loop */
    DRMS_Segment_t *firstSeg = NULL;
    HIterator_t *segIter = NULL;
    DRMS_Segment_t *orig = NULL;
    int validSegsCount;

    validSegsCount = 0;
    while ((inseg = drms_record_nextseg2(inrec, &segIter, 1, &orig)) != NULL)
    {
        if (!firstSeg)
        {
            firstSeg = inseg;
        }
        else
        {
            if (!hasSegList && !doingAllSegs)
            {
                /* We are processing only the first segment of each record. */
                break;
            }
        }

        /* Exclude segments that do not meet the requirements of this program. */
        if (drms_segment_getnaxis(inseg) == 2)
        {
            ++validSegsCount;
        }
        else
        {
            continue;
        }    
    
        if (validSegsCount == 0)
        {
            fprintf (stderr, "Error: no data segment of dimension 2 in input series %s\n", inset);
            drms_close_records(ids, DRMS_FREE_RECORD);
            drms_close_records(ods, DRMS_FREE_RECORD);
            return 1;    
        }

        fprintf(stderr, "Maproj reading rec %d\n", irec);
        image = drms_segment_read(inseg, DRMS_TYPE_FLOAT, &status);
                   /*  get needed info from record keys for mapping  */
                     /*  replace with call to solar_ephemeris_info?  */
                 
        if (!image || status != DRMS_SUCCESS)
        {
            fprintf(stderr, "Unable to read input segment (%s).\n", inseg->info->name);
            drms_close_records(ids, DRMS_FREE_RECORD);
            drms_close_records(ods, DRMS_FREE_RECORD);
            return 1; 
        }
	
        img_lon = drms_getkey_double(inrec, clon_key, &status);
        img_lat = drms_getkey_double(inrec, clat_key, &status);

        status = solar_image_info(inrec, &img_xscl, &img_yscl, &img_xc, &img_yc, &img_radius, rsun_key, apsd_key, &img_pa, &ellipse_e, &ellipse_pa, &x_invrt, &y_invrt, &need_ephem, 0);
        if (status & NO_SEMIDIAMETER)
        {
            int keystat = 0;
            double dsun_obs = drms_getkey_double(inrec, dsun_key, &keystat);

            if (keystat) 
            {
                fprintf (stderr, "Error: one or more essential keywords or values missing, needed %s; skipped\n",dsun_key);
                fprintf (stderr, "solar_image_info() returned %08x\n", status);
                continue;
            }
                       /*  set image radius from scale and distance  */
            img_radius = asin (RSUNM / dsun_obs);
            img_radius *= 3600.0 * degrad;
            img_radius /= (img_xscl <= img_yscl) ? img_xscl : img_yscl;
            status &= ~NO_SEMIDIAMETER;
        }
        
        if (status == KEYSCOPE_VARIABLE) 
        {
            fprintf(stderr, "Warning: one or more keywords expected constant are variable\n");
        } 
        else if (status) 
        {
            fprintf(stderr, "Warning: one or more essential keywords or values missing, 2nd message; skipped\n");
            fprintf(stderr, "solar_image_info() returned %08x\n", status);
            continue;
        }
    
        if (MDI_correct) 
        {
            mtrack_MDI_correct_imgctr(&img_xc, &img_yc, img_radius);
            mtrack_MDI_correct_pa(&img_pa);
        }
        
        img_xc -= 0.5 * (image->axis[0] - 1);
        img_yc -= 0.5 * (image->axis[1] - 1);
                /*  should be taken care of in solar_ephemeris_info  */
        img_lon *= raddeg;
        img_lat *= raddeg;
    
						    /*  set up output record  */
        orec = ods->records[irec];
        oseg = drms_segment_lookupnum(orec, orig->info->segnum);
        
        
        /* All of this code is segment-dependent. It used to reside outside any loop
         * (before the record loop). It was moved here, inside the segment loop. The
         * map array was allocated outside any loop, and freed outside any loop. */
        scaling_override = 0;
        bscale = bscaleIn;
        if (bscaleIn == 0.0) 
        {
            bscale = oseg->bscale;
            if (verbose) printf ("bscale set to output default: %g\n", bscale);
        }
        else
        {
            scaling_override = 1;
        }
    
        bzero = bzeroIn;
        if (isnan (bzero)) 
        {
            bzero = oseg->bzero;
            if (verbose) printf ("bzero set to output default: %g\n", bzero);
        }
        else
        {
            scaling_override = 1;
        }
    
        map->bscale = bscale;
        map->bzero = bzero;    
        /* End moved scaling code. */
        
                                 /*  perform the mapping  */
        perform_mapping(image, data, maplat, maplon, map_coslat, map_sinlat, pixct, offsun, img_lat, img_lon, img_xc, img_yc, img_radius, img_pa, ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt, MDI_correct_distort);
    
        /* We are guaranteed that image is not NULL at this point. */
        drms_free_array(image); // PHS

        if (overlay) 
        {
            for (n = 0; n < pixct; n++) 
            {
                if (ongrid[n]) data[n] = (isfinite(data[n])) ? bblank : wblank;
            }
        }

        // XXX PHS moved segment write to end to allow writing with headers.

                            /*  set output record key values  */
        // XXX PHS different copykeys flag, do statistics after copykeys.
        drms_copykeys(orec, inrec, 0, kDRMS_KeyClass_Explicit); // XXX PHS
        drms_setkey_string(orec, "requestid", RequestID); // XXX PHS
        set_statistics(oseg, map, 1); // XXX PHS

        kstat = 0;
        kstat += check_and_set_key_str   (orec, "WCSNAME", "Carrington Heliographic");
        kstat += check_and_set_key_int   (orec, "WCSAXES", 2);
        snprintf (key, sizeof(key), "CRLN-%s", wcscode[proj]);
        kstat += check_and_set_key_str   (orec, "CTYPE1", key);
        snprintf (key, sizeof(key), "CRLT-%s", wcscode[proj]);
        kstat += check_and_set_key_str   (orec, "CTYPE2", key);
        kstat += check_and_set_key_str   (orec, "CUNIT1", "deg");
        kstat += check_and_set_key_str   (orec, "CUNIT2", "deg");
        kstat += check_and_set_key_float (orec, "CRPIX1", 0.5 * map_cols + 0.5);
        kstat += check_and_set_key_float (orec, "CRPIX2", 0.5 * map_rows + 0.5);
        kstat += check_and_set_key_float (orec, "CRVAL1", clon * degrad);
        kstat += check_and_set_key_float (orec, "CRVAL2", clat * degrad);
        kstat += check_and_set_key_float (orec, "CDELT1", map_scale);
        kstat += check_and_set_key_float (orec, "CDELT2", map_scale);
        // XXX PHS, CROTA2 will be different after projection.
        kstat += check_and_set_key_float (orec, "CROTA2", map_pa * degrad);
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

        kstat += check_and_set_key_str   (orec, "Module", module_ident);
        kstat += check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
        kstat += check_and_set_key_time  (orec, "Created", CURRENT_SYSTEM_TIME);
        kstat += check_and_set_key_str   (orec, "Source", source);
        kstat += check_and_set_key_str   (orec, "Input", inset);

        // XXX PHS if running_as_export ignore setkey failures and add HISTORY for mapping params since
        // input_series"_mod" will not have the maproj specific keywords.
        if (running_as_export) 
        { 
            char buf[1024];
        
            drms_appendhistory(orec, "MapProj=", 1); drms_appendhistory(orec, mapname[proj], 0);
            sprintf(buf, "LonHG=%f, LatHG=%f", clon * degrad, clat * degrad);
            drms_appendhistory(orec, buf, 1); 
            drms_appendhistory(orec, "Interp=", 1); drms_appendhistory(orec, interpname[intrpopt], 0);
            if (overlay) 
            {
                sprintf(buf, "%f", grid_spacing);
                drms_appendhistory(orec, "GridSpacing=", 1); drms_appendhistory(orec, buf, 0);
            }
            drms_setkey_time(orec, "DATE", CURRENT_SYSTEM_TIME);
            drms_appendhistory(orec, "Source=", 1); drms_appendhistory(orec, source, 0);
        } 
        else if (kstat)
        {
            fprintf (stderr, "Error writing key value(s) to record %d in series %s\n", irec, outser);
            fprintf (stderr, "      series may not have appropriate structure\n");
            drms_close_records(ids, DRMS_FREE_RECORD);
            drms_close_records(ods, DRMS_FREE_RECORD);
            return 1;
        }

        // XXX PHS if scaling not given use input segment for bscale and bzero.
        if (!scaling_override) 
        {
            map->bscale = inseg->bscale;
            if (verbose) printf("bscale set to input record value: %g\n", bscale);
            map->bzero = inseg->bzero;
            if (verbose) printf("bzero set to input record value: %g\n", bzero);
        }

        // XXX PHS moved segment write to here to allow writing with headers.
                    /*  write map array to output record segment  */
        if (want_headers)
        {
            status = drms_segment_writewithkeys(oseg, map, 0);
        }
        else
        {
            status = drms_segment_write(oseg, map, 0);
        }

        if (status) 
        {
            drms_sprint_rec_query (recid, orec);
            fprintf (stderr, "Error writing data to record %s\n", recid);
            fprintf (stderr, "      series may not have appropriate structure\n");
            drms_close_records(ids, DRMS_FREE_RECORD);
            drms_close_records(ods, DRMS_FREE_RECORD);
            return 1;
        }
    } /* End segment loop */
    
    if (segIter)
    {
        hiter_destroy(&segIter);
    }
    
    if (!firstSeg)
    {
        /* Never found a segment for this record.*/
        fprintf(stderr, "Error: no data segment of dimension 2 in input record-set %s\n", inset);
        drms_close_records(ids, DRMS_FREE_RECORD);
        drms_close_records(ods, DRMS_FREE_RECORD);
        return 1;
    }
  } /* End record loop. */
  
  drms_close_records(ods, DRMS_INSERT_RECORD);
  drms_free_array(map);
  return 0;
}
