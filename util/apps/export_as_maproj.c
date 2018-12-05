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

/* From cartography.c. */
static double arc_distance (double lat, double lon, double latc, double lonc) 
{
  double cosa = sin (lat) * sin (latc) + cos (lat) * cos (latc) * cos (lon - lonc);
  return acos (cosa);
}

static int plane2sphere (double x, double y, double latc, double lonc, double *lat, double *lon, int projection) 
{
/*
 *  Perform the inverse mapping from rectangular coordinates x, y on a map
 *    in a particular projection to heliographic (or geographic) coordinates
 *    latitude and longitude (in radians).
 *  The map coordinates are first transformed into arc and azimuth coordinates
 *    relative to the center of the map according to the appropriate inverse
 *    transformation for the projection, and thence to latitude and longitude
 *    from the known heliographic coordinates of the map center (in radians).
 *  The scale of the map coordinates is assumed to be in units of radians at
 *    the map center (or other appropriate location of minimum distortion).
 *
 *  Arguments:
 *      x }         Map coordinates, in units of radians at the scale
 *      y }           appropriate to the map center
 *      latc        Latitude of the map center (in radians)
 *      lonc        Longitude of the map center (in radians)
 *      *lat        Returned latitude (in radians)
 *      *lon        Returned longitude (in radians)
 *      projection  A code specifying the map projection to be used: see below
 *
 *  The following projections are supported:
 *      RECTANGULAR     A "rectangular" mapping of x and y directly to
 *                      longitude and latitude, respectively; it is the
 *			normal cylindrical equidistant projection (plate
 *			carrÈe) tangent at the equator and equidistant
 *			along meridians. Central latitudes off the equator
 *			merely result in displacement of the map in y
 *			Also known as CYLEQD
 *	CASSINI		The transverse cylindrical equidistant projection
 *			(Cassini-Soldner) equidistant along great circles
 *			perpendicular to the central meridian
 *      MERCATOR        Mercator's conformal projection, in which paths of
 *                      constant bearing are straight lines
 *      CYLEQA          Lambert's normal equal cylindrical (equal-area)
 *                      projection, in which evenly-spaced meridians are
 *                      evenly spaced in x and evenly-spaced parallels are
 *                      separated by the cosine of the latitude
 *      SINEQA          The Sanson-Flamsteed sinusoidal equal-area projection,
 *                      in which evenly-spaced parallels are evenly spaced in
 *                      y and meridians are sinusoidal curves
 *      GNOMONIC        The gnomonic, or central, projection, in which all
 *                      straight lines correspond to geodesics; projection
 *                      from the center of the sphere onto a tangent plane
 *      POSTEL          Postel's azimuthal equidistant projection, in which
 *                      straight lines through the center of the map are
 *                      geodesics with a uniform scale
 *      STEREOGRAPHIC   The stereographic projection, mapping from the
 *                      antipode of the map center onto a tangent plane
 *      ORTHOGRAPHIC    The orthographic projection, mapping from infinity
 *                      onto a tangent plane
 *      LAMBERT         Lambert's azimuthal equivalent projection
 *
 *  The function returns -1 if the requested point on the map does not project
 *    back to the sphere or is not a principal value, 1 if it projects to a
 *    point on a hidden hemisphere (if that makes sense), 0 otherwise
 */
  static double latc0 = 0.0, sinlatc = 0.0, coslatc = 1.0 ;
  double r, rm, test;
  double cosr, sinr, cosp, sinp, coslat, sinlat, sinlon;
  double sinphi, cosphi, phicom;
  int status = 0;

  if (latc != latc0) {
    coslatc = cos (latc);
    sinlatc = sin (latc);
  }
  latc0 = latc;

  switch (projection) {
    case (RECTANGULAR):
      *lon = lonc + x;
      *lat = latc + y;
      if (arc_distance (*lat, *lon, latc, lonc) > M_PI_2) status = 1;
      if (fabs (x) > M_PI || fabs (y) > M_PI_2) status = -1;
      return status;
    case (CASSINI): {
      double sinx = sin (x);
      double cosy = cos (y + latc);
      double siny = sin (y + latc);
      *lat = acos (sqrt (cosy * cosy + siny * siny * sinx * sinx));
      if (y < -latc) *lat *= -1;
      *lon = (fabs (*lat) < M_PI_2) ? lonc + asin (sinx / cos (*lat)) : lonc;
      if (y > (M_PI_2 - latc) || y < (-M_PI_2 - latc)) 
        *lon = 2 * lonc + M_PI - *lon;
      if (*lon < -M_PI) *lon += 2* M_PI;
      if (*lon > M_PI) *lon -= 2 * M_PI;
      if (arc_distance (*lat, *lon, latc, lonc) > M_PI_2) status = 1;
      if (fabs (x) > M_PI || fabs (y) > M_PI_2) status = -1;
      return status;
      }
    case (CYLEQA):
      if (fabs (y) > 1.0) {
        y = copysign (1.0, y);
	status = -1;
      }
      cosphi = sqrt (1.0 - y*y);
      *lat = asin ((y * coslatc) + (cosphi * cos (x) * sinlatc));
      test = (cos (*lat) == 0.0) ? 0.0 : cosphi * sin (x) / cos (*lat);
      *lon = asin (test) + lonc;
      if (fabs (x) > M_PI_2) {
        status = 1;
        while (x > M_PI_2) {
          *lon = M_PI - *lon;
	  x -= M_PI;
        }
        while (x < -M_PI_2) {
          *lon = -M_PI - *lon;
	  x += M_PI;
        }
      }
      if (arc_distance (*lat, *lon, latc, lonc) > M_PI_2) status = 1;
      return status;
    case (SINEQA):
      cosphi = cos (y);
      if (cosphi <= 0.0) {
        *lat = y;
        *lon = lonc;
	if (cosphi < 0.0) status = -1;
        return status;
      }
      *lat = asin ((sin (y) * coslatc) + (cosphi * cos (x/cosphi) * sinlatc));
      coslat = cos (*lat);
      if (coslat <= 0.0) {
        *lon = lonc;
	if (coslat < 0.0) status = 1;
        return status;
      }
      test = cosphi * sin (x/cosphi) / coslat;
      *lon = asin (test) + lonc;
      if (fabs (x) > M_PI * cosphi) return (-1);
/*
      if (fabs (x) > M_PI_2) {
        status = 1;
        while (x > M_PI_2) {
          *lon = M_PI - *lon;
	  x -= M_PI;
        }
        while (x < -M_PI_2) {
          *lon = -M_PI - *lon;
	  x += M_PI;
        }
*/
      if (fabs (x) > M_PI_2 * cosphi) {
        status = 1;
        while (x > M_PI_2 * cosphi) {
          *lon = M_PI - *lon;
	  x -= M_PI * cosphi;
        }
        while (x < -M_PI_2 * cosphi) {
          *lon = -M_PI - *lon;
	  x += M_PI * cosphi;
        }
      }
/*
      if (arc_distance (*lat, *lon, latc, lonc) > M_PI_2) status = 1;
*/
      return status;
    case (MERCATOR):
      phicom = 2.0 * atan (exp (y));
      sinphi = -cos (phicom);
      cosphi = sin (phicom);
      *lat = asin ((sinphi  * coslatc) + (cosphi * cos (x) * sinlatc));
      *lon = asin (cosphi * sin (x) / cos (*lat)) + lonc;
      if (arc_distance (*lat, *lon, latc, lonc) > M_PI_2) status = 1;
      if (fabs (x) > M_PI_2) status = -1;
      return status;
  }
                                          /*  Convert to polar coordinates  */
  r = hypot (x, y);
  cosp = (r == 0.0) ? 1.0 : x / r;
  sinp = (r == 0.0) ? 0.0 : y / r;
                                                        /*  Convert to arc  */
  switch (projection) {
    case (POSTEL):
      rm = r;
      if (rm > M_PI_2) status = 1;
      break;
    case (GNOMONIC):
      rm = atan (r);
      break;
    case (STEREOGRAPHIC):
      rm = 2 * atan (0.5 * r);
      if (rm > M_PI_2) status = 1;
      break;
    case (ORTHOGRAPHIC):
      if ( r > 1.0 ) {
	r = 1.0;
        status = -1;
      }
      rm = asin (r);
      break;
    case (LAMBERT):
      if ( r > 2.0 ) {
        r = 2.0;
        status = -1;
      }
      rm = 2 * asin (0.5 * r);
      if (rm > M_PI_2 && status == 0) status = 1;
      break;
  }
  cosr = cos (rm);
  sinr = sin (rm);
                                          /*  Convert to latitude-longitude  */
  sinlat = sinlatc * cosr + coslatc * sinr * sinp;
  *lat = asin (sinlat);
  coslat = cos (*lat);
  sinlon = (coslat == 0.0) ? 0.0 : sinr * cosp / coslat;
/*  This should never happen except for roundoff errors, but just in case:   */
/*
  if (sinlon + 1.0 <= 0.0) *lon = -M_PI_2;
  else if (sinlon - 1.0 >= 0.0) *lon = M_PI_2;
  else
*/
  *lon = asin (sinlon);
                                     /*  Correction suggested by Dick Shine  */
  if (cosr < (sinlat * sinlatc)) *lon = M_PI - *lon;
  *lon += lonc;
  return status;
}

/* From imginfo.c. */
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

typedef struct paramdef {
  double scale;
  double offset;
  double defval;
  unsigned int statusbit;
  char name[32];
} ParamDef;


static double lookup (DRMS_Record_t *rec, ParamDef key, int *status) {
  double value = key.defval;
  int lookupstat = 0;

  value = drms_getkey_double (rec, key.name, &lookupstat);
  value = value * key.scale + key.offset;
  if (lookupstat) *status |= key.statusbit;
  if (isnan (value)) *status |= key.statusbit;
  return value;
}

static char *lookup_str (DRMS_Record_t *rec, ParamDef key, int *status) {
  DRMS_Keyword_t *keywd;
  int lstat;
  char *value;

  value = drms_getkey_string (rec, key.name, &lstat);
  if (lstat) *status |= key.statusbit;
					     /*  cadence should be constant  */
  if ((keywd = drms_keyword_lookup (rec, key.name, 1))) {
    if (keywd->info->recscope != 1) *status |= KEYSCOPE_VARIABLE;
  } else *status |= key.statusbit;
  return value;
}

static int solar_image_info (DRMS_Record_t *img, double *xscl, double *yscl,
    double *ctrx, double *ctry, double *apsd, const char *rsun_key,
    const char *apsd_key, double *pang, double *ellipse_e, double *ellipse_pa,
    int *x_invrt, int *y_invrt, int *need_ephem, int AIPS_convention) {
/*
 *  Provides the following values from the DRMS record:
 *    xscl  scale in the image column direction (arc-sec/pixel)
 *    yscl  scale in the image row direction (arc-sec/pixel)
 *    ctrx  (virtual) fractional pixel column of the center of the
 *	solar disc
 *    ctry  (virtual) fractional pixel row of the center of the
 *	solar disc
 *    apsd  apparent semi-diameter (semimajor-axis) of the solar disc, in
 *	pixel units
 *    pang  position angle of solar north relative to image vertical (y-axis,
 *	[0,0] -> [0,1]), measured westward (clockwise), in radians
 *    eecc  eccentricity of best-fit ellipse describing limb
 *    eang  position angle of best-fit ellipse describing limb, relative
 *	to direction of solar north, measured westward (clockwise), in radians
 *    xinv  0 if image is direct, 1 if flipped by columns
 *    yinv  0 if image is direct, 1 if flipped by rows
 *  If AIPS_convention is true (!=0), it is assumed that the input keywords
 *    representing position and ellipse angles are measured westward
 *    (clockwise) relative to their nominal axes; otherwise they are measured
 *    eastward (counter-clockwise) relative to their nominal axes.
 *  NO!
 *  The following data types are supported: SOHO-MDI, GONG+, Mt. Wilson MOF,
 *    SOHO-EIT, TRACE, BBSO Ha
 *  NO!
 *  If the data are not recognizably of one of these types, the function
 *    returns NO_DATA_DICT; if one or more required keywords are missing
 *    the function returns a status mask indicating which values could not
 *    be filled.  No flag is set if the image ellipse parameters are not
 *    available (unless they are required for other parameters), but the
 *    eccentricity and pericenter angle are quietly set to 0.
 */
  enum param {
    XSCL, YSCL, XUNI, YUNI, LATC, LONC, CTRX, CTRY, PANG, APSD, RSUN,
    ESMA, ESMI, EECC, EANG, XASP, YASP, PCT
  };
  enum known_plat {
    UNKNOWN
  };
  static ParamDef param[PCT];
  static double raddeg =  M_PI / 180.0;
  static double degrad = 180.0 / M_PI;
  double ella, ellb;
  int n, status = 0;
  static int scale_avail, xinv_type, yinv_type;
  static int hdrtype = UNKNOWN, lasthdr = UNKNOWN - 1;
  char *strval;
/*
 *  Set up the appropriate dictionary for interpretation of keywords
 */
  if (lasthdr != hdrtype) {
    if (lasthdr >= UNKNOWN)
      fprintf (stderr,
	  "Warning from solar_image_info(): record keywords may change\n");
    for (n = 0; n < PCT; n++) {
      sprintf (param[n].name, "No Keyword");
      param[n].scale = 1.0;
      param[n].offset = 0.0;
      param[n].defval = NAN;
      param[n].statusbit = 0;
    }
    param[RSUN].statusbit = NO_SEMIDIAMETER;
    param[APSD].statusbit = NO_SEMIDIAMETER;
    param[XSCL].statusbit = NO_XSCALE;
    param[YSCL].statusbit = NO_YSCALE;
    param[XUNI].statusbit = NO_XUNITS;
    param[YUNI].statusbit = NO_YUNITS;
    param[CTRX].statusbit = NO_XCENTERLOC;
    param[CTRY].statusbit = NO_YCENTERLOC;
    param[CTRX].statusbit = NO_XCENTERLOC;
    param[CTRY].statusbit = NO_YCENTERLOC;
    param[LATC].statusbit = NO_HELIO_LATC;
    param[LONC].statusbit = NO_HELIO_LONC;
    param[PANG].statusbit = NO_HELIO_PA;
    param[EECC].defval = 0.0;
    param[EANG].defval = 0.0;
    param[XASP].defval = 1.0;
    param[YASP].defval = 1.0;

    switch (hdrtype) {
      default:					      /*  Assume WCS HPLN/T  */
					 /*  WITH CERTAIN MDI SPECIFIC ONES  */
	scale_avail = 1;
	xinv_type = yinv_type = 0;
	sprintf (param[XUNI].name, "CUNIT1");
	sprintf (param[YUNI].name, "CUNIT2");
	sprintf (param[XSCL].name, "CDELT1");
	sprintf (param[YSCL].name, "CDELT2");
	sprintf (param[CTRX].name, "CRPIX1");
	param[CTRX].offset = -1.0;
	sprintf (param[CTRY].name, "CRPIX2");
	param[CTRY].offset = -1.0;
	*need_ephem = 0;
	strval = lookup_str (img, param[XUNI], &status);
	if (!(status & NO_XUNITS)) {
	  if (!strcmp (strval, "arcsec")) param[XSCL].scale = 1.0;
	  else if (!strcmp (strval, "arcmin")) param[XSCL].scale = 1.0 / 60.0;
	  else if (!strcmp (strval, "deg")) param[XSCL].scale = 1.0 / 3600.0;
	  else if (!strcmp (strval, "mas")) param[XSCL].scale = 1000.0;
	  else if (!strcmp (strval, "rad")) param[XSCL].scale = degrad * 3600.0;
/*
	  need_units = status & KEYSCOPE_VARIABLE;
*/
	}
	if (strval) free (strval);
	strval = lookup_str (img, param[YUNI], &status);
	if (!(status & NO_YUNITS)) {
	  if (!strcmp (strval, "arcsec")) param[YSCL].scale = 1.0;
	  else if (!strcmp (strval, "arcmin")) param[YSCL].scale = 1.0 / 60.0;
	  else if (!strcmp (strval, "deg")) param[YSCL].scale = 1.0 / 3600.0;
	  else if (!strcmp (strval, "mas")) param[YSCL].scale = 1000.0;
	  else if (!strcmp (strval, "rad")) param[YSCL].scale = degrad * 3600.0;
/*
	  need_units = status & KEYSCOPE_VARIABLE;
*/
	}
	if (strval) free (strval);
   /*  the following are appropriate for MDI, but not strictly based on WCS  */
/*
	sprintf (param[RSUN].name, "R_SUN");
	sprintf (param[APSD].name, "OBS_ASD");
*/
	strncpy (param[RSUN].name, rsun_key, 31);
	strncpy (param[APSD].name, apsd_key, 31);
	sprintf (param[PANG].name, "CROTA2");
	param[PANG].scale = -raddeg;
	if (AIPS_convention) param[PANG].scale *= -1;
	sprintf (param[ESMA].name, "S_MAJOR");
	sprintf (param[ESMI].name, "S_MINOR");
	sprintf (param[EANG].name, "S_ANGLE");
	param[EANG].scale = -raddeg;
	if (AIPS_convention) param[EANG].scale *= -1;
    }
  }
  lasthdr = hdrtype;
				    /*  Plate info: image scale, distortion  */
  *apsd = lookup (img, param[RSUN], &status);
  if (scale_avail) {
    *xscl = lookup (img, param[XSCL], &status);
    *yscl = lookup (img, param[YSCL], &status);
    if (status & NO_SEMIDIAMETER) {
      status &= ~NO_SEMIDIAMETER;
      *apsd = lookup (img, param[APSD], &status);
      if (status & NO_SEMIDIAMETER) {
        *need_ephem = 1;
      } else {
	if (!(status & (NO_XSCALE | NO_YSCALE))) {
	  *apsd /= (*xscl <= *yscl) ? *xscl : *yscl;
	  status &= ~NO_SEMIDIAMETER;
	}
      }
    }
  }
  ella = lookup (img, param[ESMA], &status);
  ellb = lookup (img, param[ESMI], &status);
  *ellipse_e = sqrt ((ella - ellb) * (ella + ellb)) / ella;
  *ellipse_pa = lookup (img, param[EANG], &status);
				    /*  Pointing (attitude: image location)  */
  *ctrx = lookup (img, param[CTRX], &status);
  *ctry = lookup (img, param[CTRY], &status);
					  /*  Position angle of solar north  */
  *pang = lookup (img, param[PANG], &status);
                                                      /*  Image orientation  */
  *x_invrt = xinv_type;
  *y_invrt = yinv_type;

  return status;
}

/* From mdistuff.c. */
#define MDI_IMG_ACPA	(1.01e-3)
#define MDI_IMG_ASPA	(-1.49e-3)

static void mtrack_MDI_correct_imgctr (double *xc, double *yc, double rsun) {
  double rs2;

  rs2 = rsun * rsun / 512.0;
  *xc -= MDI_IMG_ASPA * rs2;
  *yc -= MDI_IMG_ACPA * rs2;
}

/*
 *  mtrack_MDI_image_tip
 *
 *  Correct for ellipticity of image due to plate tipping
 *  The constants are:
 *    TIP = 2.61 deg = 0.04555
 *    EFL = 25.3 (effective focal length in units of plate half-width)
 *    PA = -56 deg
 *    SPA = sin (PA),  CPA = cos (PA)
 *    AEP = TIP / EFL
 *    BEP = TIP^2 / 4 = 5.187e-4
 *    BCPA = BEP * CPA,  BSPA = BEP * SPA
 *
 *  Bugs:
 *    There is no check that |direct| = 1; it can actually be changed as a
 *	scale factor (useful for testing)
 *    
 */

#define MDI_IMG_SPA	(-0.8290)
#define MDI_IMG_CPA	(0.5592)
#define MDI_IMG_AEP	(1.80e-3)
#define MDI_IMG_BCPA	(2.90e-4)
#define MDI_IMG_BSPA	(-4.30e-4)

static void mtrack_MDI_image_tip (double *x, double *y, int n, int direct) {
  double x0, y0, s, t;

  while (n--) {
    x0 = *x;
    y0 = *y;
    t = direct * (MDI_IMG_SPA * x0 + MDI_IMG_CPA * y0);
    s = direct * (MDI_IMG_CPA * x0 - MDI_IMG_SPA * y0);
    *x += MDI_IMG_BSPA * t;
    *y += MDI_IMG_BCPA * t;
    *x -= MDI_IMG_BCPA * s;
    *y += MDI_IMG_BSPA * s;
    t *= MDI_IMG_AEP;
    *x++ += t * x0;
    *y++ += t * y0;
  }
}

/*
 *  mtrack_MDI_image_stretch
 *
 *  Modify "plate" coordinates to account for known optical distortions
 *    in the MDI instrument
 *  It is assumed that the coordinates *x and *y are in terms of half the
 *    full plate width, i.e. that they are in the range [-1.0, 1.0] and
 *    relative to its center; this of course requires an external correction
 *    to be applied in the cases of extracted rasters and binned data.
 *    For MDI the half-plate-width is 512 * 21 um.  The 2nd-order radial
 *    correction constant is given in Kuhn et al. (Ap.J. 613, 1241; 2004)
 *    as 1.58e-3 where the radial unit is in cm. The constant used here is thus
 *    1.58e-3  * (.0021 * 512)^2
 *
 *  By ignoring the 4th-order term the function can be used for inverse
 *    as well as direct stretching.
 *
 *  The value of the integer direct specifies the direction of the
 *    transformation: +1 for a correction from "perfect" to plate coordinates, 
 *    -1 for transformation from plate to perfect
 *
 *  Bugs:
 *    There is no check that |direct| = 1; it can actually be changed as a
 *	scale factor (useful for testing)
 *
 */

#define MDI_IMG_STRETCH	(1.83e-3)

static void mtrack_MDI_image_stretch (double *x, double *y, int n, int direct) {
  double f, r2, s;

  s = direct * MDI_IMG_STRETCH;
  while (n--) {
    r2 = *x * *x + *y * *y;
    f = 1.0 + s * r2;
    *x++ *= f;
    *y++ *= f;
  }
}

#define MDI_IMG_SOHO_PA	(-0.2)

static void mtrack_MDI_correct_pa (double *pa) {
  *pa += MDI_IMG_SOHO_PA * M_PI / 180.0;
}


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

static void perform_mapping (DRMS_Array_t *img, float *map,
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

/* These check_and_set_key_* functions were copied from keystuff.c. I made a number of changes - 
 * the functions used to return success if the key did not exist. They also 
 * used to follow links and attempt to set linked keywords, but linked keywords
 * are generally in read-only records, so I think this was incorrect.
 * And they did not check the return value of the lower-level DRMS set-key functions. */
static int check_and_set_key_str(DRMS_Record_t *outRec, const char *key, const char *val) 
{
    DRMS_Keyword_t *keywd = NULL;
    char *vreq;
    int status;
    int errOut;

    errOut = 0;

    if (!(keywd = drms_keyword_lookup(outRec, key, 0))) return 1;

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

    if (!(keywd = drms_keyword_lookup(outRec, key, 0))) return 1;

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

    if (!(keywd = drms_keyword_lookup(outRec, key, 0))) return 1;

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

    if (!(keywd = drms_keyword_lookup(outRec, key, 0))) return 1;

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
  int pixct;
  int recCount;
  int kstat, status;
  int need_ephem;
  int x_invrt, y_invrt;
  int n, col, row;
  int MDI_correct, MDI_correct_distort;
  unsigned char *offsun, *ongrid;
  char source[DRMS_MAXQUERYLEN], recid[DRMS_MAXQUERYLEN];
  char module_ident[64], key[64];

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
    
    fprintf(stderr, "Maproj reading rec %d\n", irec);
    
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
        bscale = bscaleIn;
        if (bscaleIn == 0.0) 
        {
            bscale = oseg->bscale;
            if (verbose) printf ("bscale set to output default: %g\n", bscale);
        }
    
        bzero = bzeroIn;
        if (isnan (bzero)) 
        {
            bzero = oseg->bzero;
            if (verbose) printf ("bzero set to output default: %g\n", bzero);
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
