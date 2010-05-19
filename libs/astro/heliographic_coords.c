/*
 * void HeliographicLocation(TIME t, int *crot, double *L, double *B);
 * TIME HeliographicTime(int crot, double L);
 *
 * heritage from solephem (CM), WSO, and iorbit (JSOC)
 */

#include <math.h>
#include <timeio.h>

#define EPS 3.e-17
#define OMEGA_0 (75.76)     // degrees, solar equator crosses ecliptic, pre 2007 value with no aberration correction
#define OMEGA_1 (0.01397)   // degrees per Julian year (365.2500d) precession
#define SUN_INCLINATION (7.25) // Carrington solar inclination to ecliptic
#define SININCL  (0.126198969135829761)
#define COSINCL  (0.992004949679714999)

#define C  (299792.458)                  // c in km/s
#define AU  (149597870.69)   // AU in km
#define T2000           (725760032.0)    // sscan_time("2000.01.01_00") J2000 base time
#define TCARR           (-3881476800.0)  // sscan_time("1854.01.01_12:00_TAI")   Carr ref epcoh est. 
#define CARR_DEGDAY   (14.1844000)       // Adopted degrees per day, includes precession
#define CARR_ROT_SYNODIC    (27.275311 * 86400.0) // estimate of synodic carrington rotation period

#define PI      (3.141592653589793)
#define TWO_PI  (2*PI)
#define DEGRAD  (180.0/PI)
#define SID     (86400.0)

// HeliographicLocation(TIME t, int *crot, double *L, double *B);
// time t in TAI.
// It returns Carrington rotation, latitude, and longitude, as int for rotation and in degrees for angles.
// It corrects for light travel time and has coefs that correct for aberration.
// crot, L, B are Carrington rotation, longitude, and latitude of sub-observer point
// The method has been tested for 2000 to 2010 with errors never more than 0.005 degrees compared to JPL tables.
// BUT there may be long term drifts of up to 0.015 degrees per year depending on the input coordinate system
// assumptions compared to the test series.

/*
 * HeliographicLocation - returns Carrington rotation, longitude, and latitude of disk center at time t
 */
void HeliographicLocation(TIME t, int *crot, double *L, double *B)
  {
  static int firstcall = 1;
  static TIME t1900 = 0.0;
  static TIME t1950 = 0;
  static TIME t2000 = 0;
  double ecanom(double e, double manom);
  double tc19, tc1950, tc20, e, eps, perigee, lo, g, ge, v, se_long;
  double r, Omega, se_long_Omega, xs, xc, sinB, sunB0, hci_long;
  double solrots;       // solar rotations by time t
  double carr_rots;     // Rotations since Carrington Epoch
  double l;             // longitude working var (rotations)
  double tlight;        // Light travel time from Sun center to observer
  double clest;         // estimate of carrington longitude
  int rot;              // Carrington rotation number of sub-observer point
                                                                                      
  /* time in Julian centuries since 1900. */
  if (firstcall)
    {
    t1900 = sscan_time("1900.01.01_UTC");
    t1950 = sscan_time("JD_2433282.423459"); // B1950.0=1949.12.31_22:09:15_UTC
    t2000 = sscan_time("2000.01.01_00:00_UTC"); // J2000 ref time
    }
  if (isnan(t) || t < 0) return;
  tc19 = (t - t1900)/(SID*36525.0);
  tc20 = (t - t2000)/(SID*36525.0);
  tc1950 = (t - t1950)/(SID*36525.0);
  /* eccentricity of earth's orbit (B1950.0) */
  e = 0.01673011 - tc1950*(4.193e-5 + tc1950*1.26e-7);
  /* mean obliquity of the ecliptic (B1950.0) */
  eps = 0.40920619245 - tc1950*(2.27132786e-4 + tc1950*(1.5466e-8 - tc1950*8.77516e-9));
  /* mean longitude of perigee, mean equinox of date (B1950.0) */
  perigee=4.9232342127 + tc1950*( 3.0013215e-2 + tc1950*( 7.999426e-6 + tc1950*5.8178e-8));
  /* geometric mean longitude, mean equinox of date (AENA p. 539) */
  lo = 4.8815279341118791 + tc19*(628.33195094248177 + tc19*5.2796209873e-6);
  /* mean anomaly - (c.f. AENA p. 540) */
  g = lo - perigee;
  /* eccentric anomaly, mean equinox of date */
  ge = ecanom(e, g);
  /* true anomaly, mean equinox of date (SA p. 118) */
  v = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*ge));
  /* geometric longitude, mean equinox of date  (Green p 253) */
  se_long = v + perigee;
  if(se_long > TWO_PI)se_long -= TWO_PI;
  else if (se_long < 0)se_long += TWO_PI;

  /* distance in AU */
  r=(1.0-e*e)/(1.0+e*cos(se_long-perigee));
  /* longitude of ascending node of solar equator on ecliptic (AENA p. 553) */
  Omega  = (OMEGA_0 + OMEGA_1*tc20*100.0)/DEGRAD;  // Solar equator eclitic longitude (radians), no precession used
  // se_long_Omega = se_long - Omega;                 // heliocentric inertial coordinate longitude of observer (radians) 
  se_long_Omega = se_long - Omega + PI;                 // heliocentric inertial coordinate longitude of observer (radians) 
  if (se_long_Omega > TWO_PI) se_long_Omega -= TWO_PI;
  else if (se_long_Omega < 0.0) se_long_Omega += TWO_PI;
  xs = sin(se_long_Omega);
  xc = cos(se_long_Omega);
  sinB = -xs * SININCL;
  sunB0 = asin(sinB);
  // sunB0 -= hsec_z/obsdist;  // correct for Z distance, max value about 0.01 degrees, ignore inclination.
  hci_long = atan2(xs*COSINCL, xc);
  if (hci_long < 0) hci_long += TWO_PI;

  tlight = r*AU/C;
  solrots = (CARR_DEGDAY * (t - tlight - TCARR)/86400.0 -0.125)/DEGRAD;
  l = hci_long - solrots;
  l = modf(l/TWO_PI, &carr_rots);  // carrington longitude in rotations with bias
  if(l < 0) l += 1;

  // Get rotation number, first from "guess" then adjust based on accurate longitude.
  // guess is good to a few degrees.  2 rotations happened before the Carrington ref time
  carr_rots = 3.0 + (t - tlight - TCARR)/CARR_ROT_SYNODIC;
  rot = carr_rots;
  clest = (1 + rot - carr_rots);
  if ((l - clest) > 0.5)
    rot++;
  if ((clest - l) > 0.5)
    rot--;

  if (crot) *crot = rot;
  if (L) *L = 360*l;
  if (B) *B = sunB0 * DEGRAD;
  return;
  }

/* 
 * HeliographicTime returns the time at which the given location is at disk center
 */

#define C1      0.0757646111 /* days/degree at 27.27526 */
#define T1853   "1853.11.09_22:27:24_UTC"


TIME HeliographicTime(int crot, double L)
  {
  static TIME t1853=0.0;
  double i, b, l;
  int rot;
  double CThave, CTgot, clong;
  double err, CTp50m, CTm50m;
  TIME t;
  if (t1853 == 0.0) t1853 = sscan_time(T1853);
  CThave = 360.0 * crot - L;
  t = (t1853 + CThave*C1*SID);
  HeliographicLocation(t, &rot, &l, &b);
  CTgot = 360.0 * rot - l;
  err = CThave - CTgot;
  HeliographicLocation(t+3000.0, &rot, &l, &b);
  CTp50m = 360.0 * rot - l;
  HeliographicLocation(t-3000.0, &rot, &l, &b);
  CTm50m = 360.0 * rot - l;
  t += 6000 * (err/(CTp50m - CTm50m));
  return(t);
  }

/*
 * ecanom calculates the eccentric anomaly from the mean anomaly
 * and the eccentricity of date.  The method used is sufficient
 * applications of Newton's method to the equation (c.f. SA p. 117)
 *
 *	E - M - eps * sin( E ) = 0
 *
 * where E is the eccentric anomaly and M is the mean anomaly. The
 * diagnostic counter and message can be removed when this has be
 * thoroughly tested.
 */

static double ecanom(double e, double manom)
  {
  double ecm_new, ecm_old;
  ecm_old = manom + e*sin(manom);    /* first approximation */
  while (1)
    {
    ecm_new = ecm_old - (ecm_old - manom - e * sin( ecm_old))/ (1 - e * cos( ecm_old));
    if(fabs((ecm_new - ecm_old)/(ecm_new + ecm_old)) < EPS || fabs((ecm_new - ecm_old)) < EPS )
      break;
    ecm_old = ecm_new;
    }
  return(ecm_new);
  }


