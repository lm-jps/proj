/*
 *  Modified from maproj.c						~rick/src/maproj
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
 *    T
 *  Revision history is at end of file
 */

#include <jsoc_main.h>
#include "fstats.h"
#include "cartography.h"
#include "magutils.h"

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

#define PI              (M_PI)
#define RADSINDEG       (PI/180)
#define ARCSECSINRAD    (3600*180/PI)
#define RAD2ARCSEC      (648000. / M_PI)
#define RSUNM		(6.96e8)
#define INTERP_NEAREST_NEIGHBOR	(1)
#define INTERP_BILINEAR	(2)

// Defined in synop-imginfo.c:
extern int synop_solar_image_info(DRMS_Record_t *img, double *xscl, double *yscl, double *ctrx, double *ctry, double *apsd, const char *rsun_key, const char *apsd_key, double *pang, double *ellipse_e, double *ellipse_pa, int *x_invrt, int *y_invrt, int *need_ephem, int AIPS_convention);

// Defined in synop-cartography.c:
extern int synop_plane2sphere (double x, double y, double latc, double lonc, double *lat, double *lon, int projection);
extern int synop_img2sphere(double x, double y, double ang_r, double latc, double lonc, double pa, double *rho, double *lat, double *lon, double *sinlat, double *coslat, double *sig, double *mu, double *chi);
extern int synop_sphere2img(double lat, double lon, double latc, double lonc, double *x, double *y, double xcenter, double ycenter, double rsun, double peff, double ecc, double chi, int xinvrt, int yinvrt);
extern int synop_sphere2plane(double lat, double lon, double latc, double lonc, double *x, double *y, int projection);

void do_boxcar(float *image_in, float *image_out, int in_nx, int in_ny, float fscale, int power);

						       /*  module identifier  */
char *module_name = "maproj3comperrorlonat02deg";
char *module_desc = "mapping from solar images rounded at 0.2 degree at longitude";
char *version_id = "1.0";

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
  {ARG_FLOAT,   "RESCALE",     "0.1",  "Scale factor."}, // YLiu
  {ARG_STRING,	"clon_key", "CRLN_OBS", "keyname for image central longitude"}, 
  {ARG_STRING,	"clat_key", "CRLT_OBS", "keyname for image central latitude"}, 
  {ARG_STRING,	"rsun_key", "R_SUN", "keyname for image semi-diameter (pixel)"}, 
  {ARG_STRING,	"apsd_key", "RSUN_OBS", "keyname for apparent solar semi-diameter (arcsec)"}, 
  {ARG_STRING,	"dsun_key", "DSUN_OBS", "keyname for observer distance"}, 
  {ARG_FLAG,	"c",	"", "center map at center of image"}, 
  {ARG_FLAG,	"s",	"", "clon is Stonyhurst longitude"}, 
  {ARG_FLAG,	"v",	"", "verbose mode"}, 
  {}
};

		 /*  global declaration of missing to be initialized as NaN  */
float missing_val;

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

                                        /*  Error propogation for cubic convolution interpolation  */
float ccint2_error (float *f, int nx, int ny, double x, double y) {
  double  ux[4], uy[4];
  double sum, squareValue;
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
  for (i = 0; i < 4; i++){
    for (j = 0; j < 4; j++){
      squareValue = f[(iy1+i) * nx + ix1 + j] * uy[i] * ux[j];
      sum = sum + squareValue * squareValue;
                           }
                         }
  return (float)sqrt(sum);
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

                                                 /*  Error propagation for bilinear interpolation  */
float linint2_error (float *f, int cols, int rows, double x, double y) {
  double p, q, val;
  int col = (int)x, row = (int)y;
  int onerow = cols * row;
  int colp1 = col + 1, onerowp1 = onerow + cols;

  if (x < 0.0 || x > cols  || y < 0.0 || y >= rows)
    return missing_val;
  p = x - col;
  q = y - row;
  val = (1 - p) * (1 - q) * f[col + onerow] * (1 - p) * (1 - q) * f[col + onerow]
      + p * (1 - q) * f[colp1 + onerow] * p * (1 - q) * f[colp1 + onerow]
      + (1 - p) * q * f[col + onerowp1] * (1 - p) * q * f[col + onerowp1]
      + p * q * f[colp1 + onerowp1] * p * q * f[colp1 + onerowp1];
  return sqrt(val);
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

float array_imaginterp_error (DRMS_Array_t *img, double x, double y,
    int schema) { 
/*
 *  Error propogation for this image interpolation.
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
 *      is in this range or off the image, the function returns zero.
 *    The function assumes a fixed scale in both directions, so that if one
 *      dimension is larger than another the scale is applied to the larger.
 *    Only floating point data are supported by the function, and there is
 *      not even any testing for validity.
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
    return linint2_error (img->data, cols, rows, xs, ys);
  else return ccint2_error (img->data, cols, rows, xs, ys);
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

void perform_mapping_error (DRMS_Array_t *img, float *map,
    double *maplat, double *maplon, double *map_coslat, double *map_sinlat,
    int pixct, unsigned char *offsun, double latc, double lonc,
    double xc, double yc, double radius, double pa, double ellipse_e,
    double ellipse_pa, int x_invrt, int y_invrt, int interpolator,
    int MDI_correct_distort) {
/*
 *  Perform the mappings from the target heliographic coordinate sets
 *    appropriate to each output cube into the image coordinates (as
 *    corrected) for error propagation for spatial interpolation of the error values
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
    interpval = array_imaginterp_error (img, xx, yy, interpolator);
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

// -- YLiu

void do_boxcar(float *image_in, float *image_out, int in_nx, int in_ny, float fscale, int power)
{
  int iscale, nvector, vec_half;
  int inx, iny, outx, outy, i, j;
  float val;

      iscale = 1.0/fscale + 0.5;
      nvector = iscale;
      vec_half = nvector/2;

      int in_go = (iscale-1)/2.0 + 0.5;
      int out_nx = in_nx * fscale + 0.5;
      int out_ny = in_ny * fscale + 0.5;

        for (outy = 0; outy < out_ny; outy += 1)
          for (outx = 0; outx < out_nx; outx += 1)
            {
            double total = 0.0;
            double weight = 0.0;
            int nn = 0;
            for (j = 0; j < nvector; j += 1)
              {
              iny = outy*iscale + in_go + j - vec_half;
              for (i = 0; i < nvector; i += 1)
                {
                inx = outx*iscale + in_go + i - vec_half;
                if (inx >= 0 && inx < in_nx && iny >=0 && iny < in_ny)
                  {
                  val = image_in[in_nx*(iny) + inx];
                  if (!drms_ismissing_float(val))
                    {
                    double w = 1.0;
                    total += pow(w*val, power); // -- YLiu
                    weight += w;
                    nn++;
                    }
                  }
                }
              }
            if (power == 2) total = sqrt(total); //-- YL
            image_out[out_nx*outy + outx] = (nn > 0 ? total/weight : DRMS_MISSING_FLOAT);
            }

}


// -- YLiu

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids, *ods;
  double *maplat, *maplon, *map_coslat, *map_sinlat;
  double x, y, x0, y0, xstp, ystp, xrot, yrot;
  double lat, lon, cos_phi, sin_phi;
  double img_lat, img_lon;
  double img_xscl, img_yscl, img_xc, img_yc, img_radius, img_pa;
  double grid_width;
  double ellipse_e, ellipse_pa;
  int axes[2], outaxes[2];
  int img, imgct, pixct, segct;
  int isegnum, osegnum;
  int found, kstat, status;
  int need_ephem;
  int x_invrt, y_invrt;
  int n, col, row;
  int MDI_correct, MDI_correct_distort;
  unsigned char *offsun, *ongrid;
  char *input, *isegname, *osegname;
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
  char *inset = strdup (params_get_str (params, "in"));
  char *outser = strdup (params_get_str (params, "out"));
  double clat = params_get_double (params, "clat") * raddeg;
  double clon = params_get_double (params, "clon") * raddeg;
  double map_scale = params_get_double (params, "scale");
  double map_pa = params_get_double (params, "map_pa") * raddeg;
  float bscale = params_get_float (params, "bscale");
  float bzero = params_get_float (params, "bzero");
  float grid_spacing = params_get_float (params, "grid");
  float rescale = params_get_float (params, "RESCALE"); // YLiu
  int map_cols = params_get_int (params, "cols");
  int map_rows = params_get_int (params, "rows");
  int proj = params_get_int (params, "map");
  int intrpopt = params_get_int (params, "interp");
  char *clon_key = strdup (params_get_str (params, "clon_key"));
  char *clat_key = strdup (params_get_str (params, "clat_key"));
  char *rsun_key = strdup (params_get_str (params, "rsun_key"));
  char *apsd_key = strdup (params_get_str (params, "apsd_key"));
  char *dsun_key = strdup (params_get_str (params, "dsun_key"));
  int center = params_isflagset (params, "c");
  int stonyhurst = params_isflagset (params, "s");
  int verbose = params_isflagset (params, "v");
  int overlay = (isfinite (grid_spacing));
  int MDI_proc = params_isflagset (params, "M");

  TIME trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */

// additional variables
    double vrcenter;
    float xcen, ycen, rsun;
    int sunSize, vrcent;
    char *tobsstring, ttemp[64], tstr[64];
    int maskid;
    int xMask = 4096, yMask = 4096;
// end of define additional variables

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
  MDI_correct = MDI_correct_distort = MDI_proc;
  cos_phi = cos (map_pa);
  sin_phi = sin (map_pa);
  xstp = ystp = map_scale * raddeg;
  x0 = 0.5 * (1.0 - map_cols) * xstp;
  y0 = 0.5 * (1.0 - map_rows) * ystp;
  grid_width = 0.01 * grid_spacing;
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

  if (!(ods = drms_create_records (drms_env, imgct, outser, DRMS_PERMANENT,
      &status))) {
    fprintf (stderr, "Error: unable to create %d records in series %s\n",
	imgct, outser);
    fprintf (stderr, "       drms_create_records() returned status %d\n",
	status); 
    return 1;
  }
	  /*  determine appropriate output record segment (if more than one)  */

						 /*  create output map array  */
  axes[0] = map_cols;
  axes[1] = map_rows;
  outaxes[0] = map_cols * rescale + 0.5;
  outaxes[1] = map_rows * rescale + 0.5;

  pixct = map_cols * map_rows;
  maplat = (double *)malloc (pixct * sizeof (double));
  maplon = (double *)malloc (pixct * sizeof (double));
  map_coslat = (double *)malloc (pixct * sizeof (double));
  map_sinlat = (double *)malloc (pixct * sizeof (double));
  offsun = (unsigned char *)malloc (pixct * sizeof (char));
  if (overlay) ongrid = (unsigned char *)malloc (pixct * sizeof (char));
	     /*  use output series default segment scaling if not overridden  */

					  /*  process individual input mages  */
  for (img = 0; img < imgct; img++) {
							 /*  get input image  */
    DRMS_Record_t *irec, *orec;
    DRMS_Segment_t *iseg, *oseg;
    DRMS_Array_t *image = NULL, *outimage = NULL, *map = NULL, *outmap = NULL;
    float *bTotal, *bIncl, *bAzim, *disamb;
    int length[2], outDim[2];
    double xcrpix1, ycrpix1, b0, obsl0, p, S, rsunobs, imagescale, rsun, dsunobs;
    int xDim, yDim, iData, ix, jy, yOff;
    float *data, *outdata;

// additional variables
    double dSun, rSun_ref, asd, rSun;
    float  *mask;
    float cdelt;
    mask = (float *)malloc(xMask * yMask * sizeof(float));
// end of define additional variables

    irec = ids->records[img];
    drms_sprint_rec_query (source, irec);
    iseg = drms_segment_lookup(irec, "inclination");
    image = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
			   /*  get needed info from record keys for mapping  */
			     /*  replace with call to solar_ephemeris_info?  */
    img_lon = drms_getkey_double (irec, clon_key, &status);
    img_lat = drms_getkey_double (irec, clat_key, &status);

    status = synop_solar_image_info (irec, &img_xscl, &img_yscl, &img_xc, &img_yc,
	&img_radius, rsun_key, apsd_key, &img_pa, &ellipse_e, &ellipse_pa,
	&x_invrt, &y_invrt, &need_ephem, 0);
    if (status & NO_SEMIDIAMETER) {
      int keystat = 0;
      double dsun_obs = drms_getkey_double (irec, dsun_key, &keystat);
      if (keystat) {
	fprintf (stderr, "Error: one or more essential keywords or values missing; skipped\n");
	fprintf (stderr, "solar_image_info() returned %08x\n", status);
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
      fprintf (stderr, "solar_image_info() returned %08x\n", status);
      continue;
    }
    if (MDI_correct) {
      mtrack_MDI_correct_imgctr (&img_xc, &img_yc, img_radius);
      mtrack_MDI_correct_pa (&img_pa);
    }

    img_xc -= 0.5 * (image->axis[0] - 1);
    img_yc -= 0.5 * (image->axis[1] - 1);
		 	/*  should be taken care of in solar_ephemeris_info  */
    img_lon = 0.2 * round(10.0 * img_lon/2.0);
printf("lon=%f\n", img_lon);
    img_lon *= raddeg; clon = img_lon;
    img_lat *= raddeg; //clat = 0.0;

    bIncl = (float *)image->data;

    iseg = drms_segment_lookup(irec, "field");
    image = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
                           /*  get needed info from record keys for mapping  */
                             /*  replace with call to solar_ephemeris_info?  */
    bTotal = (float *)image->data;

    iseg = drms_segment_lookup(irec, "azimuth");
    image = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
                           /*  get needed info from record keys for mapping  */
                             /*  replace with call to solar_ephemeris_info?  */
    bAzim = (float *)image->data;

    iseg = drms_segment_lookup(irec, "disambig");
    image = drms_segment_read (iseg, DRMS_TYPE_FLOAT, &status);
                           /*  get needed info from record keys for mapping  */
                             /*  replace with call to solar_ephemeris_info?  */
    disamb = (float *)image->data;

// start calculate the grids

  for (n=0, row=0, y=y0; row < map_rows; row++, y += ystp) {
    for (col=0, x=x0; col < map_cols; col++, x += xstp, n++) {
      xrot = x * cos_phi - y * sin_phi;
      yrot = y * cos_phi + x * sin_phi;
      offsun[n] = synop_plane2sphere (xrot, yrot, clat, clon, &lat, &lon, proj);
      maplat[n] = lat;
      maplon[n] = lon;
      map_coslat[n] = cos (lat);
      map_sinlat[n] = sin (lat);
      if (overlay) ongrid[n] =
          near_grid_line (maplon[n], maplat[n], grid_spacing, grid_width);
    }
  }

// end the computation

// start computing noise mask

        vrcenter = drms_getkey_double(irec, "OBS_VR", &status);
        dSun = drms_getkey_double(irec, "DSUN_OBS", &status);
        rSun_ref = drms_getkey_double(irec, "RSUN_REF", &status);
        if (status) rSun_ref = 6.96e8;
        cdelt = drms_getkey_float(irec, "CDELT1", &status);  // in arcsec, assumimg dx=dy
        asd = asin(rSun_ref/dSun);
        rSun = asin(rSun_ref / dSun) * RAD2ARCSEC / cdelt;
        xcen = drms_getkey_float(irec, "CRPIX1", &status) - 1.0;
        ycen = drms_getkey_float(irec, "CRPIX2", &status) - 1.0;
        tobsstring = drms_getkey_string(irec, "INVPHMAP", &status);
        trec = drms_getkey_time(irec, "T_REC", &status);
        strcpy(ttemp, tobsstring);
        maskid = obstime2maskid(trec);

printf("vrcenter=%f, rSun=%f, xcen=%f, ycen=%f, INVPHMAP=%s\n", vrcenter, rSun, xcen, ycen, ttemp);
        noisemaskimag4twindow(xMask, yMask, xcen, ycen, rSun, vrcenter, maskid, mask);

// end of noise mask computation

// start transform of vectB to Brtp

    xcrpix1 = drms_getkey_double(irec, "CRPIX1", &status) - 1.0;
    ycrpix1 = drms_getkey_double(irec, "CRPIX2", &status) - 1.0;
    b0 = drms_getkey_double(irec, "CRLT_OBS", &status);
    obsl0 = drms_getkey_double(irec, "CRLN_OBS", &status);
    p = drms_getkey_double(irec, "CROTA2", &status);
    rsunobs = drms_getkey_double(irec, "RSUN_OBS", &status);
    imagescale = drms_getkey_double(irec, "CDELT1", &status);
    dsunobs = drms_getkey_double(irec, "DSUN_OBS", &status);
    rsun = rsunobs/imagescale;
    S = RSUNM / dsunobs;

      xDim = image->axis[0];
      yDim = image->axis[1];
      double bxHel, byHel, bzHel;
      float *bRadial, *bTheta, *bPhi;
      bRadial = (float *)malloc(xDim * yDim * sizeof(float));
      bTheta = (float *)malloc(xDim * yDim * sizeof(float));
      bPhi = (float *)malloc(xDim * yDim * sizeof(float));

        double vxx, vyy, vlon, vlat, vcoslat, vsinlat;
        double vmu, vrho, vsig, vchi;
        jy = 0; iData = 0; yOff = 0;
        for (jy = 0; jy < yDim; jy++)
          {
            ix = 0;
            vyy = (double)jy - ycrpix1;
            vyy /= img_radius;
            yOff = jy * xDim;

            for (ix = 0; ix < xDim; ix++)
              {
                iData = yOff + ix;
                vxx = (double)ix - xcrpix1;
                vxx /= img_radius;
                if (isnan(bTotal[iData]) || isnan(bAzim[iData]) || isnan(bIncl[iData]))
                  {
                    bRadial[iData] = DRMS_MISSING_FLOAT;
                    bTheta[iData] = DRMS_MISSING_FLOAT;
                    bPhi[iData] = DRMS_MISSING_FLOAT;
                    continue;
                  }

                synop_img2sphere (vxx, vyy, asin(S), b0 * RADSINDEG, obsl0 * RADSINDEG, p * RADSINDEG, &vrho, &vlat, &vlon,
                    &vsinlat, &vcoslat, &vsig, &vmu, &vchi);

                double bx = 0.0, by = 0.0, bz = 0.0;
//                if (disamb[iData] >= 4.0) bAzim[iData] += 180.0; //-- radial acute 
                if ((int)(disamb[iData]/2)%2 == 1) bAzim[iData] += 180.0; //-- random assumption
//                if ((int)(disamb[iData])%2 == 1) bAzim[iData] += 180.0; // -- potential field solution
                bx = -bTotal[iData] * sin(bIncl[iData] * RADSINDEG)
                     * sin(bAzim[iData] * RADSINDEG);
                by = bTotal[iData] * sin(bIncl[iData] * RADSINDEG)
                     * cos(bAzim[iData] * RADSINDEG);
                bz = bTotal[iData] * cos(bIncl[iData] * RADSINDEG);
                     // Azimuth angle is defined here to increase counter-clockwisely. The zero angle points 
                     // to the right. 
                     // transform the magnetic vector from image coordinates to the heliographic coordinates.

                img2helioVector (bx, by, bz, &bxHel, &byHel, &bzHel, vlon, vlat, obsl0 * RADSINDEG, b0 * RADSINDEG, p * RADSINDEG);
                bPhi[iData] = bxHel;
                bTheta[iData] = -byHel;
                bRadial[iData] = bzHel;
                }
            }

// end transform vectB

						    /*  set up output record  */
    orec = ods->records[img];
    map = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, NULL, &status);
    data = (float *)map->data;

    outmap = drms_array_create (DRMS_TYPE_FLOAT, 2, outaxes, NULL, &status);
    outdata = (float *)outmap->data;

    oseg = drms_segment_lookup(orec, "Br");
    length[0] = xDim; length[1] = yDim;
    outimage = drms_array_create(DRMS_TYPE_FLOAT, 2, length, NULL, &status);

    outimage = drms_array_create(DRMS_TYPE_FLOAT, 2, length, bRadial, &status);
 						     /*  perform the mapping  */
    perform_mapping (outimage, data, maplat, maplon, map_coslat, map_sinlat,
        pixct, offsun, img_lat, img_lon, img_xc, img_yc, img_radius, img_pa,
	ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt, MDI_correct_distort);
    if (overlay) {
      for (n = 0; n < pixct; n++) 
        if (ongrid[n]) data[n] = (isfinite (data[n])) ? bblank : wblank;
    }

    do_boxcar(data, outdata, axes[0], axes[1], rescale, 1); // YLiu
				/*  write map array to output record segment  */
    outmap->israw = 0;            // always compressed
    outmap->bzero = oseg->bzero;
    outmap->bscale = oseg->bscale;

    status = drms_segment_write (oseg, outmap, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      return 1;
    }

    double statMin, statMax, statMedn, statMean, statSig, statSkew, statKurt;
    int statNgood, ipixels;
    if (fstats(outaxes[0] * outaxes[1], outdata, &statMin, &statMax, &statMedn, &statMean, &statSig,
        &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");
// image statistics
    ipixels = outaxes[0] * outaxes[1];
    drms_setkey_int(orec, "TOTVALS_1", ipixels);
    drms_setkey_int(orec, "DATAVALS_1", statNgood);
    ipixels = outaxes[0] * outaxes[1]-statNgood;
    drms_setkey_int(orec, "MISSVALS_1", ipixels);
    drms_setkey_double(orec, "DATAMIN_1", statMin);
    drms_setkey_double(orec, "DATAMAX_1", statMax);
    drms_setkey_double(orec, "DATAMEDN_1", statMedn);
    drms_setkey_double(orec, "DATAMEAN_1", statMean);
    drms_setkey_double(orec, "DATARMS_1", statSig);
    drms_setkey_double(orec, "DATASKEW_1", statSkew);
    drms_setkey_double(orec, "DATAKURT_1", statKurt);
    drms_free_array(outimage);

    oseg = drms_segment_lookup(orec, "Bt");
    outimage = drms_array_create(DRMS_TYPE_FLOAT, 2, length, bTheta, &status);

                                                     /*  perform the mapping  */
    perform_mapping (outimage, data, maplat, maplon, map_coslat, map_sinlat,
        pixct, offsun, img_lat, img_lon, img_xc, img_yc, img_radius, img_pa,
        ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt, MDI_correct_distort);
    if (overlay) {
      for (n = 0; n < pixct; n++) 
        if (ongrid[n]) data[n] = (isfinite (data[n])) ? bblank : wblank;
    }

    do_boxcar(data, outdata, axes[0], axes[1], rescale, 1); // YLiu
                                /*  write map array to output record segment  */
    outmap->israw = 0;            // always compressed
    outmap->bzero = oseg->bzero;
    outmap->bscale = oseg->bscale;

    status = drms_segment_write (oseg, outmap, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      return 1;
    }

    if (fstats(outaxes[0] * outaxes[1], outdata, &statMin, &statMax, &statMedn, &statMean, &statSig,
        &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");
// image statistics 
    ipixels = outaxes[0] * outaxes[1];
    drms_setkey_int(orec, "TOTVALS_2", ipixels);
    drms_setkey_int(orec, "DATAVALS_2", statNgood);
    ipixels = outaxes[0] * outaxes[1]-statNgood;
    drms_setkey_int(orec, "MISSVALS_2", ipixels);
    drms_setkey_double(orec, "DATAMIN_2", statMin);
    drms_setkey_double(orec, "DATAMAX_2", statMax);
    drms_setkey_double(orec, "DATAMEDN_2", statMedn);
    drms_setkey_double(orec, "DATAMEAN_2", statMean);
    drms_setkey_double(orec, "DATARMS_2", statSig);
    drms_setkey_double(orec, "DATASKEW_2", statSkew);
    drms_setkey_double(orec, "DATAKURT_2", statKurt);
    drms_free_array(outimage);

    oseg = drms_segment_lookup(orec, "Bp");
    outimage = drms_array_create(DRMS_TYPE_FLOAT, 2, length, bPhi, &status);
                                                     /*  perform the mapping  */
    perform_mapping (outimage, data, maplat, maplon, map_coslat, map_sinlat,
        pixct, offsun, img_lat, img_lon, img_xc, img_yc, img_radius, img_pa,
        ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt, MDI_correct_distort);
    if (overlay) {
      for (n = 0; n < pixct; n++) 
        if (ongrid[n]) data[n] = (isfinite (data[n])) ? bblank : wblank;
    }

    do_boxcar(data, outdata, axes[0], axes[1], rescale, 1); // YLiu
                                /*  write map array to output record segment  */
    outmap->israw = 0;            // always compressed
    outmap->bzero = oseg->bzero;
    outmap->bscale = oseg->bscale;

    status = drms_segment_write (oseg, outmap, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      return 1;
    }

    if (fstats(outaxes[0] * outaxes[1], outdata, &statMin, &statMax, &statMedn, &statMean, &statSig,
        &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");
// image statistics 
    ipixels = outaxes[0] * outaxes[1];
    drms_setkey_int(orec, "TOTVALS_3", ipixels);
    drms_setkey_int(orec, "DATAVALS_3", statNgood);
    ipixels = outaxes[0] * outaxes[1]-statNgood;
    drms_setkey_int(orec, "MISSVALS_3", ipixels);
    drms_setkey_double(orec, "DATAMIN_3", statMin);
    drms_setkey_double(orec, "DATAMAX_3", statMax);
    drms_setkey_double(orec, "DATAMEDN_3", statMedn);
    drms_setkey_double(orec, "DATAMEAN_3", statMean);
    drms_setkey_double(orec, "DATARMS_3", statSig);
    drms_setkey_double(orec, "DATASKEW_3", statSkew);
    drms_setkey_double(orec, "DATAKURT_3", statKurt);
    drms_free_array(outimage);

//  Noise mask segment
    oseg = drms_segment_lookup(orec, "B_error");
    outimage = drms_array_create(DRMS_TYPE_FLOAT, 2, length, mask, &status);

                                                     /*  perform the mapping  */
    perform_mapping (outimage, data, maplat, maplon, map_coslat, map_sinlat,
        pixct, offsun, img_lat, img_lon, img_xc, img_yc, img_radius, img_pa,
        ellipse_e, ellipse_pa, x_invrt, y_invrt, intrpopt, MDI_correct_distort);
    if (overlay) {
      for (n = 0; n < pixct; n++)
        if (ongrid[n]) data[n] = (isfinite (data[n])) ? bblank : wblank;
    }

    do_boxcar(data, outdata, axes[0], axes[1], rescale, 2); // YLiu
                                /*  write map array to output record segment  */
    outmap->israw = 0;            // always compressed
    outmap->bzero = oseg->bzero;
    outmap->bscale = oseg->bscale;

    status = drms_segment_write (oseg, outmap, 0);
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      return 1;
    }

    if (fstats(outaxes[0] * outaxes[1], outdata, &statMin, &statMax, &statMedn, &statMean, &statSig,
        &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");
// image statistics 
    ipixels = outaxes[0] * outaxes[1];
    drms_setkey_int(orec, "TOTVALS_4", ipixels);
    drms_setkey_int(orec, "DATAVALS_4", statNgood);
    ipixels = outaxes[0] * outaxes[1]-statNgood;
    drms_setkey_int(orec, "MISSVALS_4", ipixels);
    drms_setkey_double(orec, "DATAMIN_4", statMin);
    drms_setkey_double(orec, "DATAMAX_4", statMax);
    drms_setkey_double(orec, "DATAMEDN_4", statMedn);
    drms_setkey_double(orec, "DATAMEAN_4", statMean);
    drms_setkey_double(orec, "DATARMS_4", statSig);
    drms_setkey_double(orec, "DATASKEW_4", statSkew);
    drms_setkey_double(orec, "DATAKURT_4", statKurt);
    drms_free_array(outimage);

//    drms_copykeys (orec, irec, 0, kDRMS_KeyClass_All);

     drms_copykey(orec, irec, "T_REC");
     drms_copykey(orec, irec, "T_OBS");
     drms_copykey(orec, irec, "CADENCE");
     drms_copykey(orec, irec, "DATASIGN");
     drms_copykey(orec, irec, "OBS_VR");
     drms_copykey(orec, irec, "OBS_VW");
     drms_copykey(orec, irec, "OBS_VN");
     drms_copykey(orec, irec, "DATE__OBS");
     drms_copykey(orec, irec, "QUALITY");
     drms_copykey(orec, irec, "CAMERA");    
     drms_copykey(orec, irec, "CRLN_OBS");  
     drms_copykey(orec, irec, "CRLT_OBS");
     drms_copykey(orec, irec, "CAR_ROT");
     drms_copykey(orec, irec, "DSUN_OBS");
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
    kstat += check_and_set_key_float (orec, "CRPIX1", 0.5 * map_cols * rescale + 0.5);
    kstat += check_and_set_key_float (orec, "CRPIX2", 0.5 * map_rows * rescale + 0.5);
    kstat += check_and_set_key_float (orec, "CRVAL1", clon * degrad);
    kstat += check_and_set_key_float (orec, "CRVAL2", clat * degrad);
    kstat += check_and_set_key_float (orec, "CDELT1", map_scale/rescale);
    kstat += check_and_set_key_float (orec, "CDELT2", map_scale/rescale);
    kstat += check_and_set_key_float (orec, "CROTA2", 0.0);
    if (map_pa != 0.0) {
      kstat += check_and_set_key_float (orec, "PC1_1", cos (map_pa));
/*  PC1_2 should be multiplied by CDELT2 / CDELT1  */
      kstat += check_and_set_key_float (orec, "PC1_2", sin (map_pa));
/*  PC2_1 should be multiplied by CDELT1 / CDELT2  */
      kstat += check_and_set_key_float (orec, "PC2_1", sin (map_pa));
      kstat += check_and_set_key_float (orec, "PC2_2", cos (map_pa));
    }

//    kstat += check_and_set_key_float (orec, "LonHG", img_lon * degrad);
//    kstat += check_and_set_key_float (orec, "LatHG", img_lat * degrad);
    kstat += check_and_set_key_str   (orec, "MapProj", mapname[proj]);
    kstat += check_and_set_key_float (orec, "MapScale", map_scale);
    kstat += check_and_set_key_float (orec, "Width", map_cols * rescale * map_scale);
    kstat += check_and_set_key_float (orec, "Height", map_rows * rescale * map_scale);
    kstat += check_and_set_key_float (orec, "Size", sqrt (map_rows * map_cols) * rescale * map_scale);
    kstat += check_and_set_key_float (orec, "Map_PA", map_pa / raddeg);
    kstat += check_and_set_key_float (orec, "RSunRef", 1.0e-6 * RSUNM);
    kstat += check_and_set_key_str   (orec, "Interp", interpname[intrpopt]);
    kstat += check_and_set_key_str   (orec, "Module", module_ident);
    kstat += check_and_set_key_str   (orec, "BLD_VERS", jsoc_version);
//    kstat += check_and_set_key_time  (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str   (orec, "Source", source);
    kstat += check_and_set_key_str   (orec, "Input", input);

    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record %d in series %s\n",
	  img, outser);
      fprintf (stderr, "      series may not have appropriate structure\n");
      drms_close_records (ods, DRMS_FREE_RECORD);
      return 1;
    }
        tnow = (double)time(NULL);
        tnow += UNIX_epoch;
        drms_setkey_time(orec, "DATE", tnow);

       free(bIncl); free(bTotal); free(bAzim);
       drms_free_array(image);
       drms_free_array (map);
       drms_free_array (outmap);
  }
  drms_close_records(ids, DRMS_FREE_RECORD);
  drms_close_records (ods, DRMS_INSERT_RECORD);
  return 0;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  08.04.22 	created this file, based on fastrack
 *  10.10.04	updated output keyword set and fleshed out arguments, commented
 *		added setting of WCS PC matrix elements as necessary; added
 *		input data processing, mapping with interpolation options
 *  v 0.8	frozen 2010.10.05
 *  11.05.02	added option for grid overlay; removed inappropriate defaults
 *  12.01.17	added option for removal of MDI distortion
 *  12.03.07	added promiscuous copy to target keys before careful setting
 *  v 0.9 frozen 12.08.03
 *
 *  15.09.02    modified maproj (from Rick) so that the module can remap 
 *              vector B (3 components)
 *                               -- YL 
 *  15.09.08    add error estimate based on noise mask that is related with orbital velocity.
 *                               -- YL
 *
 *
 */
