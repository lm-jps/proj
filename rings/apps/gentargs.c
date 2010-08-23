/*
 *  gentargs.c						~rick/hmi/rings/src
 *
 *  Calculate appropriate limits for tracking intervals based on.time,
 *    tile sizes, target latitudes, and B0 at time of meridian crossing
 *    and generate an appropriate set of targets
 *
 *  The program behaves in two different modes:
 *
 *    If one of the non-default grid values is selected, it prints two*n lines
 *	of output in pairs with one corresponding to the target list of
 *	latitudes for the selected time and the second to the corresponding
 *	longitudes for the selected grid. If the list contains more than about
 *	200 target points, it will be broken up into multiple pairs.
 *
 *    Otherwise, the program prints out a set of lines corresponding to the
 *	target latitudes for the selected time and the corresponding
 *	longitude ranges for tracking across a disc passage in the HMI
 *	ring-diagrams grid. In this case, the target longitudes will be those
 *	nearest a central meridian longitude in steps of 2.5 deg.
 *    The normal mode output consists of a set of lines, one for each tile
 *	size and duration. Each line contains the following numbers:
 *	  tile size (5, 15, or 30)
 *	  starting CM longitude
 *	  ending CM longitude
 *	  number of latitudes to be tracked (N) (>= 1)
 *	  N target latitudes
 *      If the starting CM longitude is smaller than the ending CM longitude,
 *	  it is assumed to refer to the previous Carrington rotation.
 *
 *    The "dense-pack" is a set of 195 regions arranged in a roughly circular
 *	grid spaced at 7.5 deg centres, centred at Lat 0 and CM Lon, and
 *	extending 105 deg in longitude and 120 deg in latitude; the targets
 *	are expected to be tracked for about 15 deg of rotation (27 hr).
 *	The targets will always be printed for the central meridian longitude
 *	at 15 deg nearest that for the target time.
 *    The "structure-pack" is a set of 39 regions arranged in three longitude
 *	columns spaced at 7.5 deg centres, centred at Lat 0 and CM Lon, and
 *	extending 30 deg in longitude and 90 deg in latitude; the targets
 *	are expected to be tracked for about 55 deg of rotation (100 hr).
 *	The targets will be printed for the central meridian longitude at
 *	15 deg nearest that for the target time.
 *    The HMI ring-diagrams grid is described at
 *	http://hmi.stanford.edu/teams/rings/pipelines.html
 *	These targets will be printed for the central meridian longitudes
 *	at 5 deg, 15 deg, or 30 deg nearest that for the target time, depending
 *	on the grid size.
 *
 *  Responsible:  Rick Bogart				RBogart@spd.aas.org
 *
 *  Usage:
 *    gentargs [-vM] [time= ] [grid= ]
 *
 *  Parameters: (type   default 	description)
 *	time	Time	(curtime-CR/3)	midtime for target calculation (midtime
 *				of tracking interval for example)
 *	grid	enum	discross	choice of output
 *				recognized values are "discross", "rdsyn05",
 *				"rdsyn15", "rdsyn30", "timed", "mdidp", "mdisp",
 *				"rd+30" (and possibly others).
 *		discross:	default mode, print values for trackings of
 *				appropriate latitude targets for disc passages
 *		rdsyn05		HMI ring-diagrams synoptic set of 5-deg tiles
 *				(2748 or 2727, varying with B0)
 *		rdsyn15		HMI ring-diagrams synoptic set of 15-deg tiles
 *				(284 or 281, varying with B0)
 *		rdsyn30		HMI ring-diagrams synoptic set of 15-deg tiles
 *				(69, varying with B0)
 *		timed20		HMI time-distance synoptic set of 25 tiles at
 *				20 deg spacings
 *		timed24		HMI time-distance synoptic set of 25 tiles at
 *				24 deg spacings
 *		mdidp		MDI "extended dense-pack" set of 195 tiles
 *		mdisp		MDI "structure-pack" set of 39 tiles
 *		rdcm05		central meridian strip from rdsyn05
 *				(64 or 65, varying with B0)
 *		rdcm15		central meridian strip from rdsyn15
 *				(20 or 19, varying with B0)
 *		rdcm30		central meridian strip from rdsyn30 (9)
 *		rdeq05		equatorial strip from rdsyn05 (63)
 *		rdeq15		equatorial strip from rdsyn15 (19)
 *		rdeq30		equatorial strip from rdsyn30 (9)
 *		rd+05		central meridian plus equator from rdsyn05
 *				(127 or 126, varying with B0)
 *		rd+15		central meridian plus equator from rdsyn15
 *				(38 or 37, varying with B0)
 *		rd+30		central meridian plus equator from rdsyn30 (17)
 *		rdx05		- not implemented -
 *		rdx15		- not implemented -
 *		rdx30		"X" pattern from rdsyn30 (lat = +/- lonCM) (17)
 *
 *  Flags
 *	M		use SOHO ephemeris for B0 calculation or conversion of
 *		time in calendar-clock format to Carrington time
 *	c		report longitudes relative to CM rather than HG
 *	v		run in verbose mode (with commentary to stderr)
 *	
 *
 *  Bugs:
 *    The 60-deg latitude list for 30-deg high B0 discross appears incomplete
 *    The "Stonyhurst" option of printing longitudes relative to CM is only
 *	implemented for the rd+30 grid
 *    A number of the options have not been implemented at all
 *    There may be problems parsing regular calendar-clock time strings
 *
 *  Revision history is at end of file.
 */

#include <drms.h>
#include <math.h>
		/*  include static versions of functions for soho_ephemeris  */
#include "soho_ephem.c"
						      /*  module identifier  */
char *module_name = "track_target_list";
char *version_id = "1.0";

ModuleArgs_t module_args[] = {
  {ARG_TIME,   "time",  "now-120deg", "midpoint of desired tracking interval"},
  {ARG_NUME, "grid", "discross", "target grid",
      "discross, rdsyn05, rdsyn15, rdsyn30, timed20, timed24, mdidp, mdisp, rd+05, rd+15, rd+30, rdx05, rdx15, rdx30, rdcm05, rdcm15, rdcm30, rdeq05, rdeq15, rdeq30"},
  {ARG_INTS,	"ar", "{}",
      "target active regions for ring diagrams structure grid"},
  {ARG_FLAG,   "M", "", "midtime Carrington longitude appropriate for MDI"}, 
  {ARG_FLAG,   "c", "", "report longitudes relative to CM rather than HG"}, 
  {ARG_FLAG,   "v", "", "run verbose (with commentary to stderr)"}, 
  {}
};

ModuleArgs_t *gModArgs = module_args;
CmdParams_t cmdparams;

enum {
  L150_60S, L150_45S, L150_30S, L150_15S,
  L150_000,
  L150_15N, L150_30N, L150_45N, L150_60N
};

enum {
  L75_75S,  L75_675S, L75_60S,  L75_525S, L75_45S,
  L75_375S, L75_30S,  L75_225S, L75_15S,  L75_075S,
  L75_000,
  L75_075N, L75_15N,  L75_225N, L75_30N,  L75_375N,
  L75_45N,  L75_525N, L75_60N,  L75_675N, L75_75N
};

enum {
  L25_825S, L25_80S,  L25_775S, L25_75S,  L25_725S,
  L25_70S,  L25_675S, L25_65S,  L25_625S, L25_60S,  L25_575S, L25_55S,
  L25_525S, L25_50S,  L25_475S, L25_45S,  L25_425S, L25_40S,  L25_375S,
  L25_35S,  L25_325S, L25_30S,  L25_275S, L25_25S,  L25_225S, L25_20S,
  L25_175S, L25_15S,  L25_125S, L25_10S,  L25_075S, L25_05S,  L25_025S,
  L25_000,
  L25_025N, L25_05N,  L25_075N, L25_10N,  L25_125N, L25_15N,  L25_175N,
  L25_20N,  L25_225N, L25_25N,  L25_275N, L25_30N,  L25_325N, L25_35N,
  L25_375N, L25_40N,  L25_425N, L25_45N,  L25_475N, L25_50N,  L25_525N,
  L25_55N,  L25_575N, L25_60N,  L25_625N, L25_65N,  L25_675N, L25_70N,
  L25_725N, L25_75N,  L25_775N, L25_80N,  L25_825N
};

enum {
  X25_85S,  X25_825S, X25_80S,  X25_775S, X25_75S,  X25_725S,
  X25_70S,  X25_675S, X25_65S,  X25_625S, X25_60S,  X25_575S, X25_55S,
  X25_525S, X25_50S,  X25_475S, X25_45S,  X25_425S, X25_40S,  X25_375S,
  X25_35S,  X25_325S, X25_30S,  X25_275S, X25_25S,  X25_225S, X25_20S,
  X25_175S, X25_15S,  X25_125S, X25_10S,  X25_075S, X25_05S,  X25_025S,
  X25_000,
  X25_025N, X25_05N,  X25_075N, X25_10N,  X25_125N, X25_15N,  X25_175N,
  X25_20N,  X25_225N, X25_25N,  X25_275N, X25_30N,  X25_325N, X25_35N,
  X25_375N, X25_40N,  X25_425N, X25_45N,  X25_475N, X25_50N,  X25_525N,
  X25_55N,  X25_575N, X25_60N,  X25_625N, X25_65N,  X25_675N, X25_70N,
  X25_725N, X25_75N,  X25_775N, X25_80N,  X25_825N, X25_85N
};

static TIME earth_meridian_crossing (double lonn, int crr) {
/*
 *  Provide the time at which a selected heliographic longitude (in degrees)
 *    is on the geocentric meridian in a given Carrington Rotation
 */
  double ephem[EPHEM_SIZE];
  double lon;
  double phi, secapp;

  secapp = ((crr - CRR0) * 360.0 - (lonn - LONN0)) * K;
  calc_sun_ephemeris (secapp, ephem, 0.0, 0.0);
		  /*  soho_ephemeris always returns a longitude in [0, 360)  */
  while (lonn >= 360.0) lonn -= 360.0;
  while (lonn < 0.0) lonn += 360.0;
  lon = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
  phi = lon - lonn;
  while (phi > 180.0) phi -= 360.0;
  while (phi < -180.0) phi += 360.0;
  while (fabs (phi) > ERR) {
    secapp += phi * K;       
    calc_sun_ephemeris (secapp, ephem, 0.0, 0.0);
    lon = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
    phi = lon - lonn;
    while (phi > 180.0) phi -= 360.0;
    while (phi < -180.0) phi += 360.0;
  }
  return secapp;
}

int struc_pack_list (double clon) {
  float lat, lon, loncm;
					      /*  generate latitude targets  */
  for (loncm = -15.0; loncm <= 15.0; loncm += 15.0) {
    for (lat = -45.0; lat <= 45.0; lat += 7.5)
      printf (" %+05.1f", lat);
  }
  printf ("\n");
					     /*  generate longitude targets  */
  for (loncm = -15.0; loncm <= 15.0; loncm += 15.0) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    for (lat = -45.0; lat <= 45.0; lat += 7.5)
      printf (" %05.1f", lon);
  }
  printf ("\n");

  return 0;
}

int dense_pack_list (double clon) {
  float lat, lon, loncm;

  for (lat = -22.5; lat <= 22.5; lat += 7.5)
    printf (" %+05.1f", lat);
  for (loncm = -45.0; loncm <= -37.5; loncm += 7.5) {
    for (lat = -37.5; lat <= 37.5; lat += 7.5)
      printf (" %+05.1f", lat);
  }
  for (lat = -45.0; lat <= 45.0; lat += 7.5)
    printf (" %+05.1f", lat);
  for (loncm = -22.5; loncm <= 22.5; loncm += 7.5) {
    for (lat = -52.5; lat <= 52.5; lat += 7.5)
      printf (" %+05.1f", lat);
  }
  for (lat = -45.0; lat <= 45.0; lat += 7.5)
    printf (" %+05.1f", lat);
  for (loncm = 37.5; loncm <= 45.0; loncm += 7.5) {
    for (lat = -37.5; lat <= 37.5; lat += 7.5)
      printf (" %+05.1f", lat);
  }
  for (lat = -22.5; lat <= 22.5; lat += 7.5)
    printf (" %+05.1f", lat);
  for (lat = -60.0; lat <= 60.0; lat += 120.0) {
    for (loncm = -15.0; loncm <= 15.0; loncm += 15.0)
      printf (" %+05.1f", lat);
  }
  printf ("\n");

  loncm = -52.5;
  lon = clon + loncm;
  while (lon < 0.0) lon += 360.0;
  while (lon >= 360.0) lon -= 360.0;
  for (lat = -22.5; lat <= 22.5; lat += 7.5)
    printf (" %05.1f", lon);
  for (loncm = -45.0; loncm <= -37.5; loncm += 7.5) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    for (lat = -37.5; lat <= 37.5; lat += 7.5)
      printf (" %05.1f", lon);
  }
  loncm = -30.0;
  lon = clon + loncm;
  while (lon < 0.0) lon += 360.0;
  while (lon >= 360.0) lon -= 360.0;
  for (lat = -45.0; lat <= 45.0; lat += 7.5)
    printf (" %05.1f", lon);
  for (loncm = -22.5; loncm <= 22.5; loncm += 7.5) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    for (lat = -52.5; lat <= 52.5; lat += 7.5)
      printf (" %05.1f", lon);
  }
  loncm = 30.0;
  lon = clon + loncm;
  while (lon < 0.0) lon += 360.0;
  while (lon >= 360.0) lon -= 360.0;
  for (lat = -45.0; lat <= 45.0; lat += 7.5)
    printf (" %05.1f", lon);
  for (loncm = 37.5; loncm <= 45.0; loncm += 7.5) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    for (lat = -37.5; lat <= 37.5; lat += 7.5)
      printf (" %05.1f", lon);
  }
  loncm = 52.5;
  lon = clon + loncm;
  while (lon < 0.0) lon += 360.0;
  while (lon >= 360.0) lon -= 360.0;
  for (lat = -22.5; lat <= 22.5; lat += 7.5)
    printf (" %05.1f", lon);

  for (lat = -60.0; lat <= 60.0; lat += 120.0) {
    for (loncm = -15.0; loncm <= 15.0; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  return 0;
}

int timed_grid_list (double clon, float step) {
  float lat, lon, loncm;
  float edge = 2.0 * step;

  for (lat = -edge; lat <= edge; lat += step) {
    for (loncm = -edge; loncm <= edge; loncm += step)
      printf (" %+05.1f", lat);
  }
  printf ("\n");
  for (lat = -edge; lat <= edge; lat += step) {
    for (loncm = -edge; loncm <= edge; loncm += step) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");
  return 0;
}

int hmi30_pack_list (double clon, double b0) {
  float lat, lon, loncm, lim;

  lat = -60.0;
  lim = (b0 < -3.625) ? 60.0 : 30.0;
  for (loncm = -lim; loncm <= lim; loncm += 30.0)
    printf (" %+05.1f", lat);
  lat = -45.0;
  lim = (b0 <= 3.625) ? 60.0 : 45.0;
  for (loncm = -lim; loncm <= lim; loncm += 15.0)
    printf (" %+05.1f", lat);
  for (lat = -30.0; lat <= 30.0; lat += 15.0) {
    for (loncm = -60.0; loncm <= 60.0; loncm += 15.0)
      printf (" %+05.1f", lat);
  }
  lat = 45.0;
  lim = (b0 >= -3.625) ? 60.0 : 45.0;
  for (loncm = -lim; loncm <= lim; loncm += 15.0)
    printf (" %+05.1f", lat);
  lat = 60.0;
  lim = (b0 > 3.625) ? 60.0 : 30.0;
  for (loncm = -lim; loncm <= lim; loncm += 30.0)
    printf (" %+05.1f", lat);
  printf ("\n");

  lat = -60.0;
  lim = (b0 < -3.625) ? 60.0 : 30.0;
  for (loncm = -lim; loncm <= lim; loncm += 30.0) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    printf (" %05.1f", lon);
  }
  lat = -45.0;
  lim = (b0 <= 3.625) ? 60.0 : 45.0;
  for (loncm = -lim; loncm <= lim; loncm += 15.0) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    printf (" %05.1f", lon);
  }
  for (lat = -30.0; lat <= 30.0; lat += 15.0) {
    for (loncm = -60.0; loncm <= 60.0; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  lat = 45.0;
  lim = (b0 >= -3.625) ? 60.0 : 45.0;
  for (loncm = -lim; loncm <= lim; loncm += 15.0) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    printf (" %05.1f", lon);
  }
  lat = 60.0;
  lim = (b0 > 3.625) ? 60.0 : 30.0;
  for (loncm = -lim; loncm <= lim; loncm += 30.0) {
    lon = clon + loncm;
    while (lon < 0.0) lon += 360.0;
    while (lon >= 360.0) lon -= 360.0;
    printf (" %05.1f", lon);
  }
  printf ("\n");

  return 0;
}

int hmi30_cm_list () {
  float lat;

  for (lat = -60.0; lat <= 60.0; lat += 15.0) printf (" %+05.1f", lat);
  printf ("\n");
  return 0;
}

int hmi30_eq_list (double clon, int stony) {
  float lon, loncm;

  if (stony) {
    for (loncm = -60.0; loncm <= 60.0; loncm += 15.0) printf (" %+05.1f", loncm);
  } else {
    for (loncm = -60.0; loncm <= 60.0; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");
  return 0;
}

int hmi30_stgrg_list (double clon, int stony) {
  float lat, lon, loncm;

  for (lat = -60.0; lat < 0.0; lat += 15.0) printf (" %+05.1f", lat);
  lat = 0.0;
  for (loncm = -60.0; loncm <= 60.0; loncm += 15.0) printf (" %+05.1f", lat);
  for (lat = 15.0; lat <= 60.0; lat += 15.0) printf (" %+05.1f", lat);
  printf ("\n");

  if (stony) {
    for (lat = -60.0; lat < 0.0; lat += 15.0) printf (" %+05.1f", 0.0);
    for (loncm = -60.0; loncm <= 60.0; loncm += 15.0) printf (" %+05.1f", loncm);
    for (lat = 15.0; lat <= 60.0; lat += 15.0) printf (" %+05.1f", 0.0);
  } else {
    for (lat = -60.0; lat < 0.0; lat += 15.0) printf (" %05.1f", clon);
    for (loncm = -60.0; loncm <= 60.0; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    for (lat = 15.0; lat <= 60.0; lat += 15.0) printf (" %05.1f", clon);
  }
  printf ("\n");
  return 0;
}

int hmi30_stand_list (double clon, int stony) {
  float lat, lon, loncm;

  for (lat = -60.0; lat < 0.0; lat += 15.0) printf (" %+05.1f %+05.1f", lat, lat);
  lat = 0.0;
  printf (" %+05.1f", lat);
  for (lat = 15.0; lat <= 60.0; lat += 15.0) printf (" %+05.1f %+05.1f", lat, lat);
  printf ("\n");

  if (stony) {
    for (loncm = -60.0; loncm < 0.0; loncm += 15.0) printf (" %+05.1f %+05.1f",
	loncm, -loncm);
    printf (" %+05.1f", 0.0);
    for (loncm = 15; loncm <= 60.0; loncm += 15.0) printf (" %+05.1f %+05.1f",
	loncm, -loncm);
  } else {
    for (loncm = -60.0; loncm < 0.0; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
      lon = clon - loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    printf (" %05.1f", clon);
    for (loncm = 15; loncm <= 60.0; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
      lon = clon - loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");
  return 0;
}

int hmi15_pack_list (double clon, double b0) {
  float lat, lon, loncm, lim;

  if (b0 < -3.625) {
    lat = -75.0;
    printf (" %+05.1f", lat);
    lat = -67.5;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0)
      printf (" %+05.1f", lat);
    lat = -60.0;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0)
      printf (" %+05.1f", lat);
    lat = -52.5;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5)
      printf (" %+05.1f", lat);
    lim = 70.0;
    for (lat = -45.0; lat <= -37.5; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0)
        printf (" %+05.1f", lat);
    }
    lim = 75.0;
    for (lat = -30.0; lat <= -15.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5)
        printf (" %+05.1f", lat);
    }
    lat = -7.5;
    lim = 67.5;
    for (loncm = -lim; loncm <= lim; loncm += 7.5)
      printf (" %+05.1f", lat);
    printf ("\n");

    lat = -75.0;
    printf (" %05.1f", clon);
    lat = -67.5;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = -60.0;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = -52.5;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lim = 70.0;
    for (lat = -45.0; lat <= -37.5; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lim = 75.0;
    for (lat = -30.0; lat <= -15.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lat = -7.5;
    lim = 67.5;
    for (loncm = -lim; loncm <= lim; loncm += 7.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    printf ("\n");

    lim = 67.5;
    for (lat = 0.0; lat <= 30.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5)
        printf (" %+05.1f", lat);
    }
    lim = 60.0;
    for (lat = 37.5; lat <= 45.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0)
        printf (" %+05.1f", lat);
    }
    lat = 52.5;
    lim = 50.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5)
      printf (" %+05.1f", lat);
    lat = 60.0;
    lim = 45.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0)
      printf (" %+05.1f", lat);
    lat = 67.5;
    lim = 20.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0)
      printf (" %+05.1f", lat);
    printf ("\n");

    lim = 67.5;
    for (lat = 0.0; lat <= 30.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lim = 60.0;
    for (lat = 37.5; lat <= 45.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lat = 52.5;
    lim = 50.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 60.0;
    lim = 45.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 67.5;
    lim = 20.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    printf ("\n");
  } else if (b0 > 3.625) {
    lat = -67.5;
    lim = 20.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0)
      printf (" %+05.1f", lat);
    lat = -60.0;
    lim = 45.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0)
      printf (" %+05.1f", lat);
    lat = -52.5;
    lim = 50.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5)
      printf (" %+05.1f", lat);
    lim = 60.0;
    for (lat = -45.0; lat <= -37.5; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0)
        printf (" %+05.1f", lat);
    }
    lim = 67.5;
    for (lat = -30.0; lat <= 0.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5)
        printf (" %+05.1f", lat);
    }
    printf ("\n");

    lat = -67.5;
    lim = 20.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = -60.0;
    lim = 45.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = -52.5;
    lim = 50.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lim = 60.0;
    for (lat = -45.0; lat <= -37.5; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lim = 67.5;
    for (lat = -30.0; lat <= 0.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    printf ("\n");

    lat = 7.5;
    lim = 67.5;
    for (loncm = -lim; loncm <= lim; loncm += 7.5)
      printf (" %+05.1f", lat);
    lim = 75.0;
    for (lat = 15.0; lat <= 30.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5)
        printf (" %+05.1f", lat);
    }
    lim = 70.0;
    for (lat = 37.5; lat <= 45.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0)
        printf (" %+05.1f", lat);
    }
    lat = 52.5;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5)
      printf (" %+05.1f", lat);
    lat = 60.0;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0)
      printf (" %+05.1f", lat);
    lat = 67.5;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0)
      printf (" %+05.1f", lat);
    lat = 75.0;
    printf (" %+05.1f", lat);
    printf ("\n");

    lat = 7.5;
    lim = 67.5;
    for (loncm = -lim; loncm <= lim; loncm += 7.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lim = 75.0;
    for (lat = 15.0; lat <= 30.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lim = 70.0;
    for (lat = 37.5; lat <= 45.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lat = 52.5;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 12.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 60.0;
    lim = 75.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 67.5;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 75.0;
    printf (" %05.1f", clon);
    printf ("\n");
  } else {
    lat = -67.5;
    lim = 40.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0)
      printf (" %+05.1f", lat);
    lat = -60.0;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0)
      printf (" %+05.1f", lat);
    lat = -52.5;
    lim = 62.5;
    for (loncm = -lim; loncm <= lim; loncm += 12.5)
      printf (" %+05.1f", lat);
    lim = 70.0;
    for (lat = -45.0; lat <= -37.5; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0)
        printf (" %+05.1f", lat);
    }
    lim = 67.5;
    for (lat = -30.0; lat <= 0.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5)
        printf (" %+05.1f", lat);
    }
    printf ("\n");

    lat = -67.5;
    lim = 40.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = -60.0;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = -52.5;
    lim = 62.5;
    for (loncm = -lim; loncm <= lim; loncm += 12.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lim = 70.0;
    for (lat = -45.0; lat <= -37.5; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lim = 67.5;
    for (lat = -30.0; lat <= 0.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    printf ("\n");

    lim = 67.5;
    for (lat = 7.5; lat <= 30.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5)
        printf (" %+05.1f", lat);
    }
    lim = 70.0;
    for (lat = 37.5; lat <= 45.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0)
        printf (" %+05.1f", lat);
    }
    lat = 52.5;
    lim = 62.5;
    for (loncm = -lim; loncm <= lim; loncm += 12.5)
      printf (" %+05.1f", lat);
    lat = 60.0;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0)
      printf (" %+05.1f", lat);
    lat = 67.5;
    lim = 40.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0)
      printf (" %+05.1f", lat);
    printf ("\n");

    lim = 67.5;
    for (lat = 7.5; lat <= 30.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 7.5) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lim = 70.0;
    for (lat = 37.5; lat <= 45.0; lat += 7.5) {
      for (loncm = -lim; loncm <= lim; loncm += 10.0) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    lat = 52.5;
    lim = 62.5;
    for (loncm = -lim; loncm <= lim; loncm += 12.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 60.0;
    lim = 60.0;
    for (loncm = -lim; loncm <= lim; loncm += 15.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    lat = 67.5;
    lim = 40.0;
    for (loncm = -lim; loncm <= lim; loncm += 20.0) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    printf ("\n");
  }

  return 0;
}

int hmi15_cm_list (double b0) {
  float lat, latmin, latmax;

  latmin = -67.5;
  latmax = 67.5;
  if (b0 < -3.625) latmin = -75.0;
  if (b0 > 3.625) latmax = 75.0;
  for (lat = latmin; lat <= latmax; lat += 7.5) printf (" %+05.1f", lat);
  printf ("\n");
  return 0;
}

int hmi15_eq_list (double clon, int stony) {
  float lon, loncm;

  if (stony) {
    for (loncm = -67.5; loncm <= 67.5; loncm += 7.5) printf (" %+05.1f", loncm);
  } else {
    for (loncm = -67.5; loncm <= 67.5; loncm += 7.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");
  return 0;
}

int hmi15_stgrg_list (double clon, double b0, int stony) {
  float lat, latmin, latmax, lon, loncm;

  latmin = -67.5;
  latmax = 67.5;
  if (b0 < -3.625) latmin = -75.0;
  if (b0 > 3.625) latmax = 75.0;

  for (lat = latmin; lat < 0.0; lat += 7.5) printf (" %+05.1f", lat);
  lat = 0.0;
  for (loncm = -67.5; loncm <= 67.5; loncm += 7.5) printf (" %+05.1f", lat);
  for (lat = 7.5; lat <= latmax; lat += 7.5) printf (" %+05.1f", lat);
  printf ("\n");

  if (stony) {
    for (lat = latmin; lat < 0.0; lat += 7.5) printf (" %+05.1f", 0.0);
    for (loncm = -67.5; loncm <= 67.5; loncm += 7.5) printf (" %+05.1f", loncm);
    for (lat = 7.5; lat <= latmax; lat += 7.5) printf (" %+05.1f", 0.0);
  } else {
    for (lat = latmin; lat < 0.0; lat += 7.5) printf (" %05.1f", clon);
    for (loncm = -67.5; loncm <= 67.5; loncm += 7.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    for (lat = 7.5; lat <= latmax; lat += 7.5) printf (" %05.1f", clon);
  }
  printf ("\n");
  return 0;
}

int hmi5_pack_list (double clon, double b0) {
  float limit[X25_85N + 1], latv[X25_85N + 1], step[X25_85N + 1];
  float lat, lon, loncm, lim;
  int nlat;

  for (nlat = X25_85S; nlat <= X25_85N; nlat++)
    latv[nlat] = -85.0 + 2.5 * nlat;
  for (nlat = X25_85S; nlat <= X25_75S; nlat++)
    step[nlat] = 10.0;
  for (nlat = X25_725S; nlat <= X25_675S; nlat++)
    step[nlat] = 7.5;
  for (nlat = X25_65S; nlat <= X25_425S; nlat++)
    step[nlat] = 5.0;
  for (nlat = X25_40S; nlat <= X25_40N; nlat++)
    step[nlat] = 2.5;
  for (nlat = X25_425N; nlat <= X25_65N; nlat++)
    step[nlat] = 5.0;
  for (nlat = X25_675N; nlat <= X25_725N; nlat++)
    step[nlat] = 7.5;
  for (nlat = X25_75N; nlat <= X25_85N; nlat++)
    step[nlat] = 10.0;

  if (b0 < -3.625) {
    limit[X25_85S] = 50.0;
    for (nlat = X25_825S; nlat <= X25_75S; nlat++)
      limit[nlat] = 70.0;
    for (nlat = X25_725S; nlat <= X25_425S; nlat++)
      limit[nlat] = 75.0;
    for (nlat = X25_40S; nlat <= X25_125S; nlat++)
      limit[nlat] = 80.0;
    for (nlat = X25_10S; nlat <= X25_075N; nlat++)
      limit[nlat] = 77.5;
    for (nlat = X25_10N; nlat <= X25_20N; nlat++)
      limit[nlat] = 75.0;
    for (nlat = X25_225N; nlat <= X25_30N; nlat++)
      limit[nlat] = 72.5;
    for (nlat = X25_325N; nlat <= X25_375N; nlat++)
      limit[nlat] = 70.0;
    limit[X25_40N] = 67.5;
    limit[X25_425N] = 65.0;
    for (nlat = X25_45N; nlat <= X25_55N; nlat++)
      limit[nlat] = 60.0;
    limit[X25_575N] = 55.0;
    for (nlat = X25_60N; nlat <= X25_625N; nlat++)
      limit[nlat] = 50.0;
    limit[X25_65N] = 40.0;
    limit[X25_675N] = 30.0;
    limit[X25_70N] = 22.5;
    limit[X25_725N] = 7.5;
    for (nlat = X25_75N; nlat <= X25_85N; nlat++)
      limit[nlat] = -1;
  } else if (b0 > 3.625) {
    for (nlat = X25_85S; nlat <= X25_75S; nlat++)
      limit[nlat] = -1;
    limit[X25_725S] = 7.5;
    limit[X25_70S] = 22.5;
    limit[X25_675S] = 30.0;
    limit[X25_65S] = 40.0;
    for (nlat = X25_625S; nlat <= X25_60S; nlat++)
      limit[nlat] = 50.0;
    limit[X25_575S] = 55.0;
    for (nlat = X25_55S; nlat <= X25_45S; nlat++)
      limit[nlat] = 60.0;
    limit[X25_425S] = 65.0;
    limit[X25_40S] = 67.5;
    for (nlat = X25_375S; nlat <= X25_325S; nlat++)
      limit[nlat] = 70.0;
    for (nlat = X25_30S; nlat <= X25_225S; nlat++)
      limit[nlat] = 72.5;
    for (nlat = X25_20S; nlat <= X25_10S; nlat++)
      limit[nlat] = 75.0;
    for (nlat = X25_075S; nlat <= X25_10N; nlat++)
      limit[nlat] = 77.5;
    for (nlat = X25_125N; nlat <= X25_40N; nlat++)
      limit[nlat] = 80.0;
    for (nlat = X25_425N; nlat <= X25_725N; nlat++)
      limit[nlat] = 75.0;
    for (nlat = X25_75N; nlat <= X25_825N; nlat++)
      limit[nlat] = 70.0;
    limit[X25_85N] = 50.0;
  } else {
    for (nlat = X25_85S; nlat <= X25_825S; nlat++)
      limit[nlat] = -1;
    limit[X25_80S] = limit[X25_80N] = 0.0;
    limit[X25_775S] = limit[X25_775N] = 30.0;
    limit[X25_75S] = limit[X25_75N] = 40.0;
    limit[X25_725S] = limit[X25_725N] = 45.0;
    for (nlat = X25_70S; nlat <= X25_675S; nlat++)
      limit[nlat] = 52.5;
    limit[X25_65S] = limit[X25_65N] = 60.0;
    for (nlat = X25_625S; nlat <= X25_525S; nlat++)
      limit[nlat] = 65.0;
    for (nlat = X25_50S; nlat <= X25_425S; nlat++)
      limit[nlat] = 70.0;
    for (nlat = X25_40S; nlat <= X25_275S; nlat++)
      limit[nlat] = 75.0;
    for (nlat = X25_25S; nlat <= X25_25N; nlat++)
      limit[nlat] = 77.5;
    for (nlat = X25_275N; nlat <= X25_40N; nlat++)
      limit[nlat] = 75.0;
    for (nlat = X25_425N; nlat <= X25_50N; nlat++)
      limit[nlat] = 70.0;
    for (nlat = X25_525N; nlat <= X25_625N; nlat++)
      limit[nlat] = 65.0;
    for (nlat = X25_675N; nlat <= X25_70N; nlat++)
      limit[nlat] = 52.5;
    for (nlat = X25_825N; nlat <= X25_85N; nlat++)
      limit[nlat] = -1;
  }

  if (b0 > 3.625) {
    for (nlat = X25_85S; nlat <= X25_55S; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
	printf (" %+05.1f", latv[nlat]);
    }
    printf ("\n");
    for (nlat = X25_85S; nlat <= X25_55S; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    printf ("\n");
  } else {
    for (nlat = X25_85S; nlat <= X25_675S; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
	printf (" %+05.1f", latv[nlat]);
    }
    printf ("\n");
    for (nlat = X25_85S; nlat <= X25_675S; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    printf ("\n");

    for (nlat = X25_65S; nlat <= X25_55S; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
	printf (" %+05.1f", latv[nlat]);
    }
    printf ("\n");
    for (nlat = X25_65S; nlat <= X25_55S; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    printf ("\n");
  }

  for (nlat = X25_525S; nlat <= X25_425S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_525S; nlat <= X25_425S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_40S; nlat <= X25_35S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_40S; nlat <= X25_35S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_325S; nlat <= X25_275S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_325S; nlat <= X25_275S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_25S; nlat <= X25_20S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_25S; nlat <= X25_20S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_175S; nlat <= X25_125S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_175S; nlat <= X25_125S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_10S; nlat <= X25_05S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_10S; nlat <= X25_05S; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_025S; nlat <= X25_025N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_025S; nlat <= X25_025N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_05N; nlat <= X25_10N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_05N; nlat <= X25_10N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_125N; nlat <= X25_175N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_125N; nlat <= X25_175N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_20N; nlat <= X25_25N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_20N; nlat <= X25_25N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_275N; nlat <= X25_325N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_275N; nlat <= X25_325N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_35N; nlat <= X25_40N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_35N; nlat <= X25_40N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  for (nlat = X25_425N; nlat <= X25_525N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
      printf (" %+05.1f", latv[nlat]);
  }
  printf ("\n");
  for (nlat = X25_425N; nlat <= X25_525N; nlat++) {
    for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");

  if (b0 < -3.625) {
    for (nlat = X25_55N; nlat <= X25_85N; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
	printf (" %+05.1f", latv[nlat]);
    }
    printf ("\n");
    for (nlat = X25_55N; nlat <= X25_85N; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
    		}
    }
    printf ("\n");
  } else {
    for (nlat = X25_55N; nlat <= X25_65N; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
	printf (" %+05.1f", latv[nlat]);
    }
    printf ("\n");
    for (nlat = X25_55N; nlat <= X25_65N; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
    		}
    }
    printf ("\n");

    for (nlat = X25_675N; nlat <= X25_85N; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat])
	printf (" %+05.1f", latv[nlat]);
    }
    printf ("\n");
    for (nlat = X25_675N; nlat <= X25_85N; nlat++) {
      for (loncm = -limit[nlat]; loncm <= limit[nlat]; loncm += step[nlat]) {
	lon = clon + loncm;
	while (lon < 0.0) lon += 360.0;
	while (lon >= 360.0) lon -= 360.0;
	printf (" %05.1f", lon);
      }
    }
    printf ("\n");
  }

  return 0;
}

int hmi5_cm_list (double b0) {
  float lat, latmin, latmax;

  latmin = -67.5;
  latmax = 67.5;
  if (b0 < -3.625) {
    latmin = -85.0;
    latmax = 72.5;
  } else if (b0 > 3.625) {
    latmin = -72.5;
    latmax = 85.0;
  } else {
    latmin = -80.0;
    latmax = 80.0;
  }
  for (lat = latmin; lat <= latmax; lat += 2.5) printf (" %+05.1f", lat);
  printf ("\n");
  return 0;
}

int hmi5_eq_list (double clon, int stony) {
  float lon, loncm;

  if (stony) {
    for (loncm = -77.5; loncm <= 77.5; loncm += 2.5) printf (" %+05.1f", loncm);
  } else {
    for (loncm = -77.5; loncm <= 77.5; loncm += 2.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
  }
  printf ("\n");
  return 0;
}

int hmi5_stgrg_list (double clon, double b0, int stony) {
  float lat, latmin, latmax, lon, loncm;

  latmin = -67.5;
  latmax = 67.5;
  if (b0 < -3.625) {
    latmin = -85.0;
    latmax = 72.5;
  } else if (b0 > 3.625) {
    latmin = -72.5;
    latmax = 85.0;
  } else {
    latmin = -80.0;
    latmax = 80.0;
  }

  for (lat = latmin; lat < 0.0; lat += 2.5) printf (" %+05.1f", lat);
  lat = 0.0;
  for (loncm = -77.5; loncm <= 77.5; loncm += 2.5) printf (" %+05.1f", lat);
  for (lat = 2.5; lat <= latmax; lat += 2.5) printf (" %+05.1f", lat);
  printf ("\n");

  if (stony) {
    for (lat = latmin; lat < 0.0; lat += 2.5) printf (" %+05.1f", 0.0);
    for (loncm = -77.5; loncm <= 77.5; loncm += 2.5) printf (" %+05.1f", loncm);
    for (lat = 2.5; lat <= latmax; lat += 2.5) printf (" %+05.1f", 0.0);
  } else {
    for (lat = latmin; lat < 0.0; lat += 2.5) printf (" %05.1f", clon);
    for (loncm = -77.5; loncm <= 77.5; loncm += 2.5) {
      lon = clon + loncm;
      while (lon < 0.0) lon += 360.0;
      while (lon >= 360.0) lon -= 360.0;
      printf (" %05.1f", lon);
    }
    for (lat = 2.5; lat <= latmax; lat += 2.5) printf (" %05.1f", clon);
  }
  printf ("\n");
  return 0;
}

void hmi_transit_list (int size, float *lat, int *trlat, int latct, float clon,
    int lonspan, int ct) {
  float lonmax, lonmin;
  int n;

  lonmax = clon + 0.5 * lonspan;
  lonmin = clon - 0.5 * lonspan;
  while (lonmax > 360.0) lonmax -= 360.0;
  while (lonmin < 0.0) lonmin += 360.0;
  printf ("%d %.2f %.2f %d", size, lonmax, lonmin, ct);
  for (n = 0; n < latct; n++) {
    if (trlat[n] == lonspan) printf (" %+05.1f", lat[n]);
  }
  printf ("\n");
}

int main (int argc, char **argv) {
  ModuleArgs_t *arg = module_args;
  CmdParams_t *params = &cmdparams;
  TIME modt;
  double ephem[EPHEM_SIZE];
  double b0, cl;
  double au, vr, vn, vw;
  double degrad = 180.0 / M_PI;
  float lat25[L25_825N + 1], lat75[L75_75N + 1], lat150[L150_60N + 1];
  int trlat25[L25_825N + 1], trlat75[L75_75N + 1], trlat150[L150_60N + 1];
  int ct[181];
  int cr, latct, lonspc, lonstp, n, nmin, nmax;
  int dosp25 = 0, dosp75 = 0, dosp150  = 0;
  char tbuf[64];
  enum grid_choice {DISCROSS, RD_SYN05, RD_SYN15, RD_SYN30, TD_SYN20, TD_SYN24,
      MDI_DP, MDI_SP, RD_PLUS05, RD_PLUS15, RD_PLUS30, RD_X05, RD_X15, RD_X30,
      RD_CM05, RD_CM15, RD_CM30, RD_EQ05, RD_EQ15, RD_EQ30,
      TARGETS, LAST_CHOICE};
  int trackct = 0, track5ct = 0, track16ct = 0, track32ct = 0;

  int status = cmdparams_parse (params, argc, argv);
  if (status == CMDPARAMS_QUERYMODE) {
    cmdparams_usage (argv[0]);
    return 0;
  } else if (status == CMDPARAMS_NODEFAULT) {
    fprintf (stderr, "For usage, type %s [-H|--help]\n", argv[0]);
    return 0;
  } else if (status < 0) {
    fprintf (stderr, "Error: Command line parsing failed. Aborting.\n");
    fprintf (stderr, "For usage, type %s [-H|--help]\n", argv[0]);
    return 1;
  }

  TIME tmid = params_get_time (params, "time");
  int target_grid = params_get_int (params, "grid");
  int formdi = cmdparams_isflagset (params, "M");
  int stony = cmdparams_isflagset (params, "c");
  int verbose = cmdparams_isflagset (params, "v");
  char *tstr = strdup (params_get_str (params, "time"));

  int do_dense_pack = target_grid == MDI_DP;
  int do_struc_pack = target_grid == MDI_SP;
  int do_hmi30_pack = target_grid == RD_SYN30;
  int do_hmi15_pack = target_grid == RD_SYN15;
  int do_hmi5_pack = target_grid == RD_SYN05;
  int do_timed_grid20 = target_grid == TD_SYN20;
  int do_timed_grid24 = target_grid == TD_SYN24;
  int do_cm30 = target_grid == RD_CM30;
  int do_cm15 = target_grid == RD_CM15;
  int do_cm05 = target_grid == RD_CM05;
  int do_eq30 = target_grid == RD_EQ30;
  int do_eq15 = target_grid == RD_EQ15;
  int do_eq05 = target_grid == RD_EQ05;
  int do_stgrg30 = target_grid == RD_PLUS30;
  int do_stgrg15 = target_grid == RD_PLUS15;
  int do_stgrg05 = target_grid == RD_PLUS05;
  int do_stand30 = target_grid == RD_X30;
  int do_stand15 = target_grid == RD_X15;
  int do_stand05 = target_grid == RD_X05;

  lonstp = (do_dense_pack || do_struc_pack || do_hmi15_pack|| do_cm15 || do_eq15) ? 150 :
      (do_hmi30_pack|| do_eq30 || do_cm30) ? 300 : \
      (do_timed_grid20 || do_timed_grid24) ? 1 : 25;
					   /*  parse the target time string  */
  if (sscanf (tstr, "%d:%lf", &cr, &cl) != 2) {
				      /* requested time not in format CR:CL  */
    if (sscanf (tstr, "CT%lf", &cl) != 1) {
				        /* requested time not in format CT*  */
      if (time_is_invalid (tmid)) {
	if (strcmp (tstr, "now-120deg")) {
          printf ("unrecognized time \"%s\"\n", tstr);
	  return 1;
	}
					 /*  default time: current - 1/3 CR  */
	tmid = CURRENT_SYSTEM_TIME;
	if (verbose)
	  fprintf (stderr, "finding time corresponding to 120 deg in past: ");
	sprint_time (tbuf, tmid, "Z", -1);
	if (formdi) {
	    sprint_time (tbuf, tmid, "Z", -1);
	    fprintf (stderr, "\nWarning: requested time %s outside range of SOHO ephemris\n",
		tbuf);
	} else {
	  calc_sun_ephemeris (tmid, ephem, 0.0, 0.0);
	  cl = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
	}
	cr = carrington_rots (tmid, formdi);
	cl += 120.0;
	if (cl >= 360.0) {
          cl -= 360.0;
	  cr--;
	}
	if (verbose) fprintf (stderr, "%d:%05.1f\n", cr, cl);
      }
    } else {
						    /*  convert CT to CR:CL  */
      cr = cl;
      cl = 360.0 * (1.0 + cr - cl);
		      /*  time corresponding to requested meridian crossing  */
      tmid = (formdi) ? SOHO_meridian_crossing (cl, cr) :
	    earth_meridian_crossing (cl, cr);
      if (verbose) {
	sprint_time (tbuf, tmid, "Z", -1);
	fprintf (stderr, "requested crossing occurs at %s\n", tbuf);
      }
    }
  } else {
		      /*  time corresponding to requested meridian crossing  */
    tmid = (formdi) ? SOHO_meridian_crossing (cl, cr) :
	earth_meridian_crossing (cl, cr);
    if (verbose) {
      sprint_time (tbuf, tmid, "Z", -1);
      fprintf (stderr, "requested crossing occurs at %s\n", tbuf);
    }
  }
  if (formdi) soho_ephemeris (tmid, &au, &b0, &cl, &vr, &vn, &vw, &modt);
  else {
    calc_sun_ephemeris (tmid, ephem, 0.0, 0.0);
    cl = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
    b0 = degrad * ephem[EPHEM_B0];
  }
  cr = carrington_rots (tmid, formdi);
    /*  adjust time to next longitude slot appropriate for requested blocks  */
  lonspc = 10.0 * (cl + 0.05 * lonstp);
  if (lonspc % lonstp > 0) lonspc -= lonspc % lonstp;
  cl = 0.1 * lonspc;
  if (cl >= 360.0) {
    cl -= 360.0;
    cr--;
  }

  if (do_dense_pack) return dense_pack_list (cl);
  if (do_struc_pack) return struc_pack_list (cl);
  if (do_hmi30_pack) return hmi30_pack_list (cl, b0);
  if (do_hmi15_pack) return hmi15_pack_list (cl, b0);
  if (do_hmi5_pack) return hmi5_pack_list (cl, b0);
  if (do_timed_grid20) return timed_grid_list (cl, 20.0);
  if (do_timed_grid24) return timed_grid_list (cl, 24.0);
  if (do_cm30) return hmi30_cm_list ();
  if (do_cm15) return hmi15_cm_list (b0);
  if (do_cm05) return hmi5_cm_list (b0);
  if (do_eq30) return hmi30_eq_list (cl, stony);
  if (do_eq15) return hmi15_eq_list (cl, stony);
  if (do_eq05) return hmi5_eq_list (cl, stony);
  if (do_stgrg30) return hmi30_stgrg_list (cl, stony);
  if (do_stgrg15) return hmi15_stgrg_list (cl, b0, stony);
  if (do_stgrg05) return hmi5_stgrg_list (cl, b0, stony);
  if (do_stand30) return hmi30_stand_list (cl, stony);
  if (do_stand15 || do_stand05) {
    fprintf (stderr, "Error: selected grid option is unimplemented\n");
    return 1;
  }

  if (verbose) {
    fprintf (stderr, "target time adjusted to %d:%05.1f ", cr, cl);
    sprint_time (tbuf, tmid, "Z", -1);
    fprintf (stderr, ": %s, B0 = %.3f\n", tbuf, b0);
  }

  for (n = L25_825S; n <= L25_825N; n++) {
    trlat25[n] = 0;
    lat25[n] = -82.5 + n * 2.5;
  }
  for (n = L75_75S; n <= L75_75N; n++) {
    trlat75[n] = 0;
    lat75[n] = -75.0 + n * 7.5;
  }
  for (n = L150_60S; n <= L150_60N; n++) {
    trlat150[n] = 0;
    lat150[n] = -60.0 + n * 15;
  }

  if (lonspc % 450 == 0) {
					/*  45 deg spacing, cadence 8  */
    dosp75 = 1;
    if (b0 < -3.625) trlat75[L75_75S] = 90;
    else if (b0 > 3.625) trlat75[L75_75N] = 90;
  }
  if (lonspc % 300 == 0) {
					/*  30 deg spacing, cadence 12  */
    dosp150 = 1;
    if (b0 < -3.625) trlat150[L150_60S] = 100;
    else if (b0 <= 3.625) trlat150[L150_60S] = trlat150[L150_60N] = 65;
    else trlat150[L150_60N] = 100;
  }
  if (lonspc % 200 == 0) {
					/*  20 deg spacing, cadence 18  */
    dosp75 = 1;
    if (b0 < -3.625) trlat75[L75_675S] = 120;
    else if (b0 <= 3.625) trlat75[L75_675S] = trlat75[L75_675N] = 75;
    else trlat75[L75_675N] = 120;
  }
  if (lonspc % 150 == 0) {
					/*  15 deg spacing, cadence 24  */
    dosp150 = dosp75 = 1;
    trlat150[L150_000] = 130;
    if (b0 < -3.625) {
      trlat150[L150_15S] = 130;
      trlat150[L150_15N] = 120;
      trlat150[L150_30S] = 130;
      trlat150[L150_30N] = 110;
      trlat150[L150_45S] = 125;
      trlat150[L150_45N] = 85;
      trlat75[L75_60S] = 135;
      trlat75[L75_60N] = 70;
    } else if (b0 <= 3.625) {
      trlat150[L150_15S] = trlat150[L150_15N] = 130;
      trlat150[L150_30S] = trlat150[L150_30N] = 120;
      trlat150[L150_45S] = trlat150[L150_45N] = 105;
      trlat150[L150_60S] = trlat150[L150_60N] = 65;
      trlat75[L75_60S] = trlat75[L75_60N] = 105;
    } else {
      trlat150[L150_15N] = 130;
      trlat150[L150_15S] = 120;
      trlat150[L150_30N] = 130;
      trlat150[L150_30S] = 110;
      trlat150[L150_45N] = 125;
      trlat150[L150_45S] = 85;
      trlat75[L75_60N] = 135;
      trlat75[L75_60S] = 70;
    }
  }
  if (lonspc % 120 == 0) {
					/*  12 deg spacing, cadence 30  */
    dosp75 = 1;
    if (b0 < -3.625) {
      trlat75[L75_525S] = 140;
      trlat75[L75_525N] = 95;
    } else if (b0 <= 3.625) {
      trlat75[L75_525S] = trlat75[L75_525N] = 120;
    } else {
      trlat75[L75_525S] = 95;
      trlat75[L75_525N] = 140;
    }
  }
  if (lonspc % 100 == 0) {
					/*  10 deg spacing, cadence 36  */
    dosp75 = dosp25 = 1;
    if (b0 < -3.625) {
      trlat75[L75_45S] = trlat75[L75_375S] = 145;
      trlat75[L75_375N] = 120;
      trlat75[L75_45N] = 110;
    } else if (b0 <= 3.625) {
      trlat75[L75_45S] = trlat75[L75_45N] = 130;
      trlat75[L75_375S] = trlat75[L75_375N] = 135;
    } else {
      trlat75[L75_45S] = 110;
      trlat75[L75_375S] = 120;
      trlat75[L75_375N] = trlat75[L75_45N] = 145;
    }

    if (b0 < -3.625) {
      trlat25[L25_825S] = 90;
      trlat25[L25_80S] = 115;
      trlat25[L25_775S] = 130;
      trlat25[L25_75S] = 135;
    } else if (b0 <= 3.625) {
      trlat25[L25_75S] = trlat25[L25_75N] = 65;
    } else {
      trlat25[L25_825N] = 90;
      trlat25[L25_80N] = 115;
      trlat25[L25_775N] = 130;
      trlat25[L25_75N] = 135;
    }
  }
  if (lonspc % 75 == 0) {
					/*  7.5 deg spacing, cadence 48  */
    dosp75 = dosp25 = 1;
    for (n = L75_075S; n <= L75_075N; n++) trlat75[n] = 150;
    if (b0 < -3.625) {
      trlat75[L75_30S] = 145;
      trlat75[L75_225S] = 150;
      trlat75[L75_15S] = trlat75[L75_075S] = 145;
      trlat75[L75_000] = 145;
      trlat75[L75_075N] = trlat75[L75_15N] = 140;
      trlat75[L75_225N] = 135;
      trlat75[L75_30N] = 130;
    } else if (b0 <= 3.625) {
      for (n = L75_15S; n <= L75_15N; n++) trlat75[n] = 145;
      trlat75[L75_30S] = trlat75[L75_225S] = 140;
      trlat75[L75_225N] = trlat75[L75_30N] = 140;
    } else {
      trlat75[L75_30N] = 145;
      trlat75[L75_225N] = 150;
      trlat75[L75_15N] = trlat75[L75_075N] = 145;
      trlat75[L75_000] = 145;
      trlat75[L75_075S] = trlat75[L75_15S] = 140;
      trlat75[L75_225S] = 135;
      trlat75[L75_30S] = 130;
    }

    if (b0 < -3.625) {
      trlat25[L25_725S] = 140;
      trlat25[L25_70S] = 145;
      trlat25[L25_675S] = 150;
      trlat25[L25_675N] = 55;
      trlat25[L25_70N] = 20;
    } else if (b0 <= 3.625) {
      trlat25[L25_725S] = trlat25[L25_725N] = 85;
      trlat25[L25_70S] = trlat25[L25_70N] = 100;
      trlat25[L25_675S] = trlat25[L25_675N] = 110;
    } else {
      trlat25[L25_70S] = 20;
      trlat25[L25_675S] = 55;
      trlat25[L25_675N] = 150;
      trlat25[L25_70N] = 145;
      trlat25[L25_725N] = 140;
    }
  }
  if (lonspc % 50 == 0) {
					      /*  5 deg spacing, cadence 72  */
    dosp25 = 1;
    if (b0 < -3.625) {
      trlat25[L25_65S] = 150;
      for (n = L25_625S; n <= L25_575S; n++) trlat25[n] = 155;
      for (n = L25_55S; n <= L25_425S; n++) trlat25[n] = 160;
      trlat25[L25_425N] = trlat25[L25_45N] = 130;
      trlat25[L25_475N] = 125;
      trlat25[L25_50N] = 120;
      trlat25[L25_525N] = 115;
      trlat25[L25_55N] = 110;
      trlat25[L25_575N] = 105;
      trlat25[L25_60N] = 95;
      trlat25[L25_625N] = 90;
      trlat25[L25_65N] = 75;
    } else if (b0 <= 3.625) {
      trlat25[L25_65S] = trlat25[L25_65N] = 120;
      trlat25[L25_625S] = trlat25[L25_625N] = 125;
      trlat25[L25_60S] = trlat25[L25_60N] = 130;
      trlat25[L25_575S] = trlat25[L25_575N] = 130;
      trlat25[L25_55S] = trlat25[L25_55N] = 135;
      for (n = L25_525S; n <= L25_475S; n++) trlat25[n] = 140;
      for (n = L25_45S; n <= L25_425S; n++) trlat25[n] = 145;
      for (n = L25_425N; n <= L25_425N; n++) trlat25[n] = 145;
      for (n = L25_475N; n <= L25_525N; n++) trlat25[n] = 140;
    } else {
      trlat25[L25_65S] = 75;
      trlat25[L25_625S] = 90;
      trlat25[L25_60S] = 95;
      trlat25[L25_575S] = 105;
      trlat25[L25_55S] = 110;
      trlat25[L25_525S] = 115;
      trlat25[L25_50S] = 120;
      trlat25[L25_475S] = 125;
      trlat25[L25_425S] = trlat25[L25_45S] = 130;
      for (n = L25_425N; n <= L25_55N; n++) trlat25[n] = 160;
      for (n = L25_575N; n <= L25_625N; n++) trlat25[n] = 155;
      trlat25[L25_65N] = 150;
    }
/*
for (n = L25_825S; n <= L25_825N; n++) printf ("%5.1f: %d\n", lat25[n], trlat25[n]);
*/
  }
  if (lonspc % 25 == 0) {
					   /*  2.5 deg spacing, cadence 144  */
    dosp25 = 1;
    for (n = L25_05S; n <= L25_05N; n++) trlat25[n] = 155;
    if (b0 < -3.625) {
      for (n = L25_40S; n <= L25_175S; n++) trlat25[n] = 160;
      for (n = L25_15S; n <= L25_05N; n++) trlat25[n] = 155;
      for (n = L25_075N; n <= L25_175N; n++) trlat25[n] = 150;
      for (n = L25_20N; n <= L25_275N; n++) trlat25[n] = 145;
      for (n = L25_30N; n <= L25_35N; n++) trlat25[n] = 140;
      trlat25[L25_375N] = trlat25[L25_40N] = 135;
    } else if (b0 <= 3.625) {
      trlat25[L25_40S] = trlat25[L25_40N] = 145;
      for (n = L25_375S; n <= L25_225S; n++) trlat25[n] = 150;
      for (n = L25_20S; n <= L25_20N; n++) trlat25[n] = 155;
      for (n = L25_225N; n <= L25_375N; n++) trlat25[n] = 150;
    } else {
      for (n = L25_40N; n >= L25_175N; n--) trlat25[n] = 160;
      for (n = L25_15N; n >= L25_05S; n--) trlat25[n] = 155;
      for (n = L25_075S; n >= L25_175S; n--) trlat25[n] = 150;
      for (n = L25_20S; n >= L25_275S; n--) trlat25[n] = 145;
      for (n = L25_30S; n >= L25_35S; n--) trlat25[n] = 140;
      trlat25[L25_375S] = trlat25[L25_40S] = 135;
    }
  }
	    /*  this should never happen with the possible values of lonstp  */
  if (!dosp25 && !dosp75 && !dosp150) {
    fprintf (stderr, "requested longitude %.1f not on any standard cadence\n",
	cl);
    return 1;
  }
				       /*  32 deg regions at 15 deg spacing  */
  latct = L150_60N + 1;
  nmin = 180;
  nmax = 0;
  for (n = nmax; n <= nmin; n++) ct[n] = 0;
  for (n = 0; n < latct; n++) {
    if (trlat150[n]) {
      if (trlat150[n] < nmin) nmin = trlat150[n];
      if (trlat150[n] > nmax) nmax = trlat150[n];
      ct[trlat150[n]]++;
    }
  }
  for (n = nmin; n <= nmax; n++) {
    if (ct[n]) {
      hmi_transit_list (30, lat150, trlat150, latct, cl, n, ct[n]);
      trackct += ct[n];
      track32ct += ct[n];
    }
  }
				      /*  16 deg regions at 7.5 deg spacing  */
  latct = L75_75N + 1;
  nmin = 180;
  nmax = 0;
  for (n = nmax; n <= nmin; n++) ct[n] = 0;
  for (n = 0; n < latct; n++) {
    if (trlat75[n]) {
      if (trlat75[n] < nmin) nmin = trlat75[n];
      if (trlat75[n] > nmax) nmax = trlat75[n];
      ct[trlat75[n]]++;
    }
  }
  for (n = nmin; n <= nmax; n++) {
    if (ct[n]) {
      hmi_transit_list (15, lat75, trlat75, latct, cl, n, ct[n]);
      trackct += ct[n];
      track16ct += ct[n];
    }
  }
				    /*  5.12 deg regions at 2.5 deg spacing  */
  latct = L25_825N + 1;
  nmin = 180;
  nmax = 0;
  for (n = nmax; n <= nmin; n++) ct[n] = 0;
  for (n = 0; n < latct; n++) {
    if (trlat25[n]) {
      if (trlat25[n] < nmin) nmin = trlat25[n];
      if (trlat25[n] > nmax) nmax = trlat25[n];
      ct[trlat25[n]]++;
    }
  }  

  for (n = nmin; n <= nmax; n++) {
    if (ct[n]) {
      hmi_transit_list (5, lat25, trlat25, latct, cl, n, ct[n]);
      trackct += ct[n];
      track5ct += ct[n];
    }
  }

  if (trackct && verbose)
    fprintf (stderr, "%d:%05.1f\t%3d %3d %3d\n", cr, cl, track32ct, track16ct,
	track5ct);
  return 0;
}

/*
 *  Revision History: all mods by R. Bogart unless otherwise noted
 *
 *  	09.11.28	Created program, based on genmtscr and
 *		~rick/mdi/rings/gendptargs
 *  v 0.7
 *  	10.02.16	Added option for MDI "structure_pack" (a few longitudes
 *		tracked for longer durations)
 *		Fixed separation of stdout and stderr for verbose output
 *		Fixed bug in determination of target CM longitude
 *  v 0.7 frozen 2010.02.17
 *  v 0.8
 *	10.02.18	Defined new format for output of hmi transit info and
 *		implemented
 *  v 0.8 frozen 2010.03.10
 *  v 0.9
 *	10.04.13	Added options for "daily" HMI-packs of 5-, 15- and
 *		30-deg tiles
 *	10.04.21	Added option for time-distance grid
 *		Changed option selection to enumerated "grid" argument rather
 *		than flags; added several options for future expansion (crosses,
 *		strips)
 *  v 0.9 frozen 2010.04.23
 *  v 1.0
 *	10.05.12	Implemented low B0 targets for daily HMI 5-deg and
 *		15-deg tilesets, CM strips, EQ strip for 30-deg tileset
 *	10.05.17	Fixed bug in 15-deg tileset for |lat| 52.5 strips
 *		implemented other EQ strips, + cross for 30
 *	10.05.22	Fixed CM strip for 5-deg tiles
 *	10.05.29	Fixed error in high-lat lon lists for 15-deg grid at
 *		"equinox"
 *	10.07.01	Changed timed option to timed20, added timed24 option
 *	10.07.20	Added -c option for reporting of longitudes relative
 *		to central meridian (for rd+30 only so far)
 *	10.08.10	Implemented + cross for 5 and 15-deg tiles and x cross
 *		for 30-deg tiles (all with -c option); made -c option available
 *		for equatorial strips
 *  v 1.0 frozen 2010.08.23
 */
