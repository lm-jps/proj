#define EARTH_EPHEM_VERSION_NUM     (1.2)

#define JPL_EPHEM_VERSION "jpleph_v11.c"

#include JPL_EPHEM_VERSION

#define AU_METERS	(1.49597870691e11)

void hc_earth_ephem (double jd, double pos[6], char *table) {
/*
 *  Calls suitably modified JPL-provided FORTRAN code to get J2000.0 positions
 *    and velocities from table.  Table units are AU and AU/day; returned
 *    units are in m and m/s
 */
  long i;
  int targ = 3, ctr = 11;
  if (strlen (table)) pleph_ext (jd, targ, ctr, 0, pos, table);
  else pleph (jd, targ, ctr, pos);
  for (i=0; i<3; i++) pos[i] *= AU_METERS;
  for (i=3; i<6; i++) pos[i] *= AU_METERS / 86400.0;
}

int earth_ephemeris (TIME obs_time, double *r, double *b, double *l,
    double *vr, double *vn, double *vw) {
/*
 *  Interpolate heliocentric position-velocity vector in summary SOHO orbit
 *    table
 *
 *  Bugs:
 *    Location of tables is hardcoded at compile time, should be passed as
 *      parameter
 *    For times outside the bounds of tables, or if table is unreadable,
 *	returns 1
 *    Function pleph_ext() leaves file open for possible next use.
 */
  FILE *tstbl;
  double gcitdt, jd, posvel[6];
  double ex, ey, ez, evx, evy, evz;
  double evr, evw, evn;
  double lapp, w, lx, rx, bx, cx, ce, vx, vy, vz;
  int icar;
  char tabledir[256], tablename[256], tstr[64];

  static double txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz;
  static int first_call = 1;
  static char lastable[256];

  double deg2rad = M_PI / 180;
  double rad2deg = 1.0 / deg2rad;
        /*  Longitude of ascending node of solar equator, 2000.0 = 16.13 deg  */
  double alpha = 16.13 * deg2rad;
     /*  Inclination of solar equator on mean equator of 2000.0 = 26.16 deg  */
   /*  Is this a typo? 26.13 deg has been used in code and seems to be right  */
  double delta = 26.13 * deg2rad;
  double car0 = 1650.0;
  double carrate = 4.2434255e-7;			   /*  rotations/sec  */

#ifdef JPL_EPHEM_TABLEDIR
  strcpy (tabledir, JPL_EPHEM_TABLEDIR);
#else
  fprintf (stderr, "Error (earth_ephem): unknown location for JPL tables\n");
  return 1;
#endif
				 /*  Check that JD is within ephemeris range  */
  gcitdt = obs_time + 32.184;
  jd = gcitdt / 86400.0 + 2443144.5;
  if (jd < 2378448.5) {
    sprint_time (tstr, gcitdt, "Z", -1);
    fprintf (stderr, "Error: time %s is earlier than any tabulated values\n",
	tstr);
    return 1;
  } else if (jd < 2433264.5)
    sprintf (tablename, "%s/unxp1800.406", tabledir);
  else if (jd < 2451536.5)
    sprintf (tablename, "%s/unxp1950.405", tabledir);
  else if (jd < 2469808.5)
    sprintf (tablename, "%s/unxp2000.405", tabledir);
  else {
    sprint_time (tstr, gcitdt, "Z", -1);
    fprintf (stderr, "Error: time %s is later than any tabulated values\n",
	tstr);
    return 1;
  }
      /*  Set up matrix for transformation from J2000.0 to solar coordinates  */
  if (first_call) {
    txx = cos (alpha);
    txy = sin (alpha);
    txz = 0.0;
    tyx = -sin (alpha) * cos (delta);
    tyy = cos (alpha) * cos (delta);
    tyz = sin (delta);
    tzx = sin (alpha) * sin (delta);
    tzy = -cos (alpha) * sin (delta);
    tzz = cos (delta);
					    /*  also check table readability  */
    if (!(tstbl = fopen (tablename, "r"))) {
      fprintf (stderr, "Error: can\'t read ephemeris table %s\n", tablename);
      return 1;
    }
    fclose (tstbl);
    strcpy (lastable, tablename);
    first_call = 0;
  } else {
						 /*  has time range changed?  */
    if (strcmp (tablename, lastable)) {
						      /*  if so, check again  */
      if (!(tstbl = fopen (tablename, "r"))) {
	fprintf (stderr, "Error: can\'t read ephemeris table %s\n", tablename);
	return 1;
      }
      fclose (tstbl);
      strcpy (lastable, tablename);
    }
  }
	     /*  Get the desired earth center quantities for the target time  */
  hc_earth_ephem (jd, posvel, tablename);
  ex = posvel[0];
  ey = posvel[1];
  ez = posvel[2];
  evx = posvel[3];
  evy = posvel[4];
  evz = posvel[5];

  rx = sqrt (ex*ex + ey*ey + ez*ez);
  *r = rx / AU_METERS;

  bx = asin ((tzx*ex + tzy*ey + tzz*ez) / rx);
  *b = rad2deg * bx;
	/*  lx is the Carrington rotation number + (1-fractional longitude)
				this is a continously increasing number  */
  lapp = car0 + carrate * gcitdt;
  w = (84.10 + 14.1844*(jd - 2451545.0))*deg2rad;
  lx = atan2 ((tyx*ex + tyy*ey + tyz*ez),(txx*ex + txy*ey + txz*ez));
  cx = 1.0 - 0.5 * fmod (lx-w, 2 * M_PI) / M_PI;
  icar = lapp - cx + 0.5;			/*  round off lapp - lx  */
  ce = cx + icar;
  icar = ce + 1;
  *l = 360.0 * (icar - ce);

  vx = txx*evx + txy*evy + txz*evz;
  vy = tyx*evx + tyy*evy + tyz*evz;
  vz = tzx*evx + tzy*evy + tzz*evz;
  *vr = cos(lx)*cos(bx)*vx + sin(lx)*cos(bx)*vy + sin(bx)*vz;
  *vw = -sin(lx)*vx + cos(lx)*vy;
  *vn = -cos(lx)*sin(bx)*vx - sin(lx)*sin(bx)*vy + cos(bx)*vz;

  return 0;
}

#define E_CRR0    (1650)
#define E_ERR     (0.0000770)
#define E_LONN0   (357.894)
			/*  At 0 of internal time (1977.01.01_00:00:00_TAI),
						       CR:CL = 1650:357.894  */
#define E_K       (6546.07452)
		/*  K is the average time in seconds for the sun to rotate
			    1 degree in longitude: the synodic period / 360
				    based on a tropical year of 365.24220 d  */
		
TIME earth_meridian_crossing (double lonn, int crr) {
/*
 *  Provide the time at which a selected heliographic longitude (in degrees)
 *    is on the geocentric meridian in a given Carrington Rotation
 */
  double au, lon, lat, vr, vn, vw;
  double phi, secapp;

#ifndef JPL_EPHEM_TABLEDIR
  return 0.0 / 0.0;
#endif
  secapp = ((crr - E_CRR0) * 360.0 - (lonn - E_LONN0)) * E_K;
  earth_ephemeris (secapp, &au, &lat, &lon, &vr, &vn, &vw);
		 /*  earth_ephemeris always returns a longitude in [0, 360)  */
  while (lonn >= 360.0) lonn -= 360.0;
  while (lonn < 0.0) lonn += 360.0;
  phi = lon - lonn;
  while (phi > 180.0) phi -= 360.0;
  while (phi < -180.0) phi += 360.0;
  while (fabs (phi) > E_ERR) {
    secapp += phi * E_K;       
    earth_ephemeris (secapp, &au, &lat, &lon, &vr, &vn, &vw);
    phi = lon - lonn;
    while (phi > 180.0) phi -= 360.0;
    while (phi < -180.0) phi += 360.0;
  }
  return secapp;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  09.12.08 or earlier	created file
 *  v 1.0 frozen 2019.12.08
 *  20.05.30	Added version number definition; replaced define of table
 *		directory with DRMS serverdefs define; trap undefined case
 */
