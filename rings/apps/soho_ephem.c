/*
 *  static versions of functions for soho_ephemeris for direct inclusion
 */

#ifndef LOCAL_SOLEPHEM_INCL
#include "solephem.c"
#endif

static int sds_flip (void *data, long length, int size) {
  long i,j;
  char *c, t;
  unsigned short *s;
  unsigned int *l;
  
  if (data == NULL) return 1;

  if (size == 2) {
    s = (unsigned short *)data;
    for (i=0; i<length; i++, s++) *s = (*s<<8) | (*s>>8);
  } else if (size == 4) {
    l = (unsigned int *)data;
    for (i=0; i<length; i++, l++)
      *l = (*l<<24) | ((*l&0xff00)<<8) | ((*l&0xff0000)>>8) | (*l>>24);
  } else {
    c = (char *)data;
    for (i=0; i<length; i++, c+=size)
      for (j=0; j<(size/2); j++) {
	t = *(c+j);
	*(c+j) = *(c+size-j-1);
	*(c+size-j-1) = t;
      }
  }

  return 0;
}

#define SUMMARY_FILE	("/home/soi/CM/tables/ephemeris/summary")
#define TABLE_DT    (600.0)
#define SUMMARY_REC (6)

int soho_ephemeris (TIME obs_time, double *r, double *b, double *l,
    double *vr, double *vn, double *vw, TIME *table_mod_time) {
/*
 *  Interpolate heliocentric position-velocity vector in summary SOHO orbit
 *    table
 *
 *  Bugs:
 *    Doesn't deal with gaps.
 *    For times outside the bounds of table, or if table is unreadable,
 *      earth (L0L0) ephemeris is returned,
 *    Leaves file open for possible next use.
 */
  static FILE *sumfile = NULL;
  static int firstcall = 1;
  static TIME tab_start, tab_stop;
  struct stat stat_buf;
  double obs_vel0[SUMMARY_REC], obs_vel1[SUMMARY_REC];
  double dt, f0, f1;
  int gap, offset, reclen = SUMMARY_REC * sizeof (double);

  if (firstcall) {
    stat (SUMMARY_FILE, &stat_buf);
    *table_mod_time = stat_buf.st_mtime + UNIX_EPOCH;
    sumfile = fopen (SUMMARY_FILE, "r");
    if (sumfile) {
      fread (&tab_start, sizeof (TIME), 1, sumfile);
      fread (&tab_stop, sizeof (TIME), 1, sumfile);
#if __BYTE_ORDER == __LITTLE_ENDIAN
      sds_flip ((void *)&tab_start, 1, sizeof (TIME));
      sds_flip ((void *)&tab_stop, 1, sizeof (TIME));
#endif
    }
    firstcall = 0;
  }
    
  if (!sumfile || obs_time < tab_start || obs_time > (tab_stop + TABLE_DT)) {
    double p1, p2;
    double ephem[EPHEM_SIZE];
    calc_sun_ephemeris (obs_time, ephem, 0.0, 0.0);
    *r = 0.99 * ephem[EPHEM_DIST];
    *l = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
    *b = ephem[EPHEM_B0] * 180.0 / M_PI;
    *vw = ephem[EPHEM_VEARTH];
    *vn = ephem[EPHEM_B0V];
    calc_sun_ephemeris (obs_time - 1.0e5, ephem, 0.0, 0.0);
    p1 = ephem[EPHEM_DIST];
    calc_sun_ephemeris (obs_time + 1.0e5, ephem, 0.0, 0.0);
    p2 = ephem[EPHEM_DIST];
    *vr = (1.496e11)*(p2 - p1)/2.0e5;
    return (sumfile ? -1 : -2);
  }
  if (obs_time >= tab_stop) {
    fseek (sumfile, -reclen, SEEK_END);
    fread (obs_vel0, sizeof (double), SUMMARY_REC, sumfile);
#if __BYTE_ORDER == __LITTLE_ENDIAN
    sds_flip ((void *)obs_vel0, SUMMARY_REC, sizeof (double));
#endif
    *vw = obs_vel0[4];
    *vn = obs_vel0[5];
    *r = obs_vel0[0];
    *vr = obs_vel0[3];
    *b = obs_vel0[2];
    *l = obs_vel0[1];
    return (1);
  }
  offset = reclen * floor ((obs_time - tab_start) / TABLE_DT) + 
    2 * sizeof (TIME);
  fseek (sumfile, offset, SEEK_SET);
  fread (obs_vel0, sizeof (double), SUMMARY_REC, sumfile);
  fread (obs_vel1, sizeof (double), SUMMARY_REC, sumfile);
#if __BYTE_ORDER == __LITTLE_ENDIAN
  sds_flip ((void *)obs_vel0, SUMMARY_REC, sizeof (double));
  sds_flip ((void *)obs_vel1, SUMMARY_REC, sizeof (double));
#endif
  dt = fmod ((obs_time - tab_start), TABLE_DT) / TABLE_DT;
  f1 = dt; f0 = 1.0 - f1;
  *r = f0 * obs_vel0[0] + f1 * obs_vel1[0];
  while (obs_vel1[1] > obs_vel0[1]) obs_vel1[1] -= 360.0;
  *l = f0 * obs_vel0[1] + f1 * obs_vel1[1];
  while (*l < 0.0) *l += 360.0;
  *b = f0 * obs_vel0[2] + f1 * obs_vel1[2];
  *vr = f0 * obs_vel0[3] + f1 * obs_vel1[3];
  *vw = f0 * obs_vel0[4] + f1 * obs_vel1[4];
  *vn = f0 * obs_vel0[5] + f1 * obs_vel1[5];
  return (0);
}

#define CARRINGTON_EPOCH    ("1853.11.09_12:00")
#define CARR_ROT_SYNODIC    (27.275311 * 86400.0)
			  /*  Synodic period from sidereal period of 25.38d
					     and tropical year of 365.2422d  */
double carrington_rots (TIME obs_time, int forsoho) {
  TIME upd;
  double ephem[EPHEM_SIZE];
  double car_rot, clong, clest;
  double r, lat, lon, vr, vn, vw;
  int carr_ct;

  car_rot = 1.0 + (obs_time - sscan_time (CARRINGTON_EPOCH)) / CARR_ROT_SYNODIC;
  carr_ct = car_rot;
  clest = car_rot - carr_ct;
  if (forsoho)
    soho_ephemeris (obs_time, &r, &lat, &lon, &vr, &vn, &vw, &upd);
  else {
    calc_sun_ephemeris (obs_time, ephem, 0.0, 0.0);
    lon = 360.0 - fmod (ephem[EPHEM_L0], 360.0);
  }
  clong = 1.0 - lon / 360.0;
  if ((clong - clest) > 0.5) carr_ct--;
  if ((clest - clong) > 0.5) carr_ct++;
  return (carr_ct + clong);
}

#define CRR0    (1650)
#define ERR     (0.0000770)
#define LONN0   (357.895)
#define K       (6546.07464)
		/*  K is the average time in seconds for the sun to rotate
			    1 degree in longitude: the synodic period / 360  */

TIME SOHO_meridian_crossing (double lonn, int crr) {
/*
 *  Provide the time at which a selected heliographic longitude (in degrees)
 *    is on the SOHO meridian in a given Carrington Rotation
 */
  TIME modt;
  double au, lon, lat, vr, vn, vw;
  double phi, secapp;

  secapp = ((crr - CRR0) * 360.0 - (lonn - LONN0)) * K;
  soho_ephemeris (secapp, &au, &lat, &lon, &vr, &vn, &vw, &modt);
		  /*  soho_ephemeris always returns a longitude in [0, 360)  */
  while (lonn >= 360.0) lonn -= 360.0;
  while (lonn < 0.0) lonn += 360.0;
  phi = lon - lonn;
  while (phi > 180.0) phi -= 360.0;
  while (phi < -180.0) phi += 360.0;
  while (fabs (phi) > ERR) {
    secapp += phi * K;       
    soho_ephemeris (secapp, &au, &lat, &lon, &vr, &vn, &vw, &modt);
    phi = lon - lonn;
    while (phi > 180.0) phi -= 360.0;
    while (phi < -180.0) phi += 360.0;
  }
  return secapp;
}
