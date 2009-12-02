/*
 *  static versions of functions for solephem for direct inclusion
 */

#ifndef LOCAL_SOLEPHEM_INCL
#define LOCAL_SOLEPHEM_INCL
#endif

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

#define EPS 3.e-17
static double ecanom (double e, double manom) {
  int i;
  double ecm_new, ecm_old;

  ecm_old = manom + e*sin(manom);    /* first approximation */
  i = 0;	/* diagnostic counter */
  for (;;) {
	++i;
	ecm_new = ecm_old - (ecm_old - manom - e * sin( ecm_old)) /
			(1 - e * cos( ecm_old));
	if (fabs((ecm_new - ecm_old)/(ecm_new + ecm_old)) < EPS ||
	  fabs((ecm_new - ecm_old)) < EPS ) break;
	ecm_old = ecm_new;
    if (i > 16) {
      fprintf (stderr, "ecanom: excessive iterations: %d\nold,new", i);
      fprintf (stderr, " %g %g %g\n", ecm_old, ecm_new, ecm_new-ecm_old);
      break;
    }
  }
  return (ecm_new);
}

#define EPHEM_OFFSET	(-11865398400.0) /*  (sscan_time ("1601.01.01_00:00:00"))  */

#define EPHEM_SIZE	35	/* size of array to declare for solephem etc. */

#define EPHEM_T		1	/* time used for calculations		*/
#define EPHEM_TCENT	2	/* time in centuries since 1900.0	*/
#define EPHEM_JD	3	/* time as Julian date			*/
#define EPHEM_ST	4	/* sidereal time			*/
#define EPHEM_RA	5	/* right ascention (geometric)		*/
#define EPHEM_DEC	6	/* declination (geometric)		*/
#define EPHEM_CLONG	7	/* geometric longitude			*/
#define EPHEM_L0	8	/* carrington degrees of central point of disk */
#define EPHEM_B0	9	/* latitude of central point of disk	*/
#define EPHEM_P		10	/* position angle of northern extremity of rotation */
				/* axis. + to east from north (in sky)	*/
#define EPHEM_DIST	11	/* true distance to sun (AU)		*/
#define EPHEM_RSUN	12	/* true semi-diameter of sun (arc-sec)	*/
#define EPHEM_VEARTH	13	/* angular velocity of earth (m/s)	*/
#define EPHEM_EOT	14	/* equation of time			*/
#define EPHEM_B0V	15	/* b0 wobble (m per s)			*/
#define EPHEM_PHI	16	/* phi = aux angle for deltav		*/
#define EPHEM_SINP	17	/* sin(p)				*/
#define EPHEM_COSP	18	/* cos(p)				*/
#define EPHEM_PAROT	19	/* pa of rotation axis wrt ecliptic     */
				  /*  the following defined in delta_v  */
#define EPHEM_RHO	20	/* rho					*/
#define EPHEM_L		21	/* Longitude on disk			*/
#define EPHEM_B		22	/* heliographic latitude		*/
#define EPHEM_R		23	/* r/Rsun				*/
#define EPHEM_SECZ	24	/* 1/cos(zentih angle)			*/
#define EPHEM_SINL	25	/* sin(L)				*/
#define EPHEM_SINB	26	/* sin(B)				*/
				   /*  the following defined in deltaI  */
#define EPHEM_IA	27	/* I/I_zenith				*/
#define EPHEM_LD	28	/* I/I_center - limb darkening 5250	*/
#define EPHEM_AIRM	29	/* airmass				*/
#define EPHEM_DV	30	/* deltav				*/
			    /*  the following defined in sun_ephemeris  */
#define EPHEM_OBS_LON	31	/* Observer longitude			*/
#define EPHEM_OBS_LAT	32	/* Observer latitude			*/
				     /*  following used in mdi_point()  */

static void solephem (double time, double *ephem) {
/*
  This procedure finds the ephemerides for the sun at the time "time".
  several quantities are computed and returned in the array ephem
  These quantities are (except as noted) the same as found in the
  American Ephemeris and Nautical Almanac (AENA).  See that document and
  the Explanatory Supplement to the Astronomical Ephemeris and the American
  Ephemeris and Nautical Almanac (ESAEAENA - 1973 edition) as well as         
  Smart "Textbook on Spherical Astronomy" (SA) for the methods used.
  all angles (including sideral time) in radians.  Note that ephem is
  assumed to point to 19 doubles of storage.
  ephem[1]:	 time used for computation 
	2:	 t = time in centuries of 36525 days since 1900.0
	3:	 jd = Julian date
	4:	 st = siderial time (nearest 0 UT)
	5:	 ra = right ascention (geometric)
	6:	 dec = declination (geometric)
	7:	 clong = geometric longitude
	8:	 l0 = carrington degrees of central point of disk
	9:	 b0 = latitude of central point of disk
       10:	 p = position angle of northern extremity of rotation
			axis. + to east from north (in sky)
       11:	 r = true distance to sun (AU)[
       12:	 rsun = true semi-diameter of sun (arc-sec)
		 uncorrected for irradiance (c.f. AENA p. 541) 
       13:	 vearth = angular velocity of earth * 1 AU. (m/s)
       14:	 eot = equation of time
       15:	 b0 wobble (m per s)
       16:	 phi = aux angle for deltav
       17:	 sin(p)
       18:	 cos(p)
       19:       position angle of Sun wrt the ecliptic

	In some AENA formulae, seemingly superfluous digits are the result
	of the use of the AENA numbers given in degrees, minutes and seconds
	and converted to radians with a full double precision representation
	of pi.

	Note that because of the AENA formulae used, some of the numbers do
	not make sense if a time before the middle of the 20th century is given.

	References to SA are to 6th edition (1977) and to AENA the 1979 issue.
	Adapted and modified from version on rc4000 - josh

  **    References to Green are to "Spherical Astronomy", l985, by RM
	Green.     11/86  rob
			*/

  double dwork1, dwork2, dsd;
  long year,modsd,ly,sd,doy;
  double tenths,tenths_sid_mean,tt;
  double days_1900,t,mean_sid_time,hams,rams,ra,lo,eot,pi,two_pi,
	clong,omega,xx,xs,xc,cosclo,b0,l0,l00,perigee,e,rday,
	clong_omega,p,r,rsun,sinincl,cosincl,tanincl,
	sinb0,dec,g,ge,v,eps,y;


/* note that internally some calculations are performed in tenths of seconds.
 * the origin of this is an older definition of time used on a previous
 * machine - when you know all about it, you fix it
 */

	/* constants, rday is radians per tenth of second */

  pi=3.14159265358979324;
  two_pi=6.28318530717995865;
  rday=two_pi/864000;
/* sine, cosine and tangent of the inclination of the sun's rotation axis to 
 * the ecliptic - all those digits are flourish - 7.25 degrees, AENA, p. 553
 */
	sinincl = 0.126198969135829761;
	cosincl = 0.992004949679714999;
	tanincl = 0.127216068001046929;
/* tenths of seconds in UT day and UT day number from 1601 */
/*  Assume times refer to old baseline; the explicit correction is not made  */
/*
  time -= sscan_time ("1601.01.01_00:00:00");
*/
	tenths = modf((double)time/86400.0, &dsd);
	sd = dsd;
	tenths *= 864000.0;
/* year and day of year - for carrington degrees */
	modf((double)sd/36524.0, &dwork1);
	modf((double)sd/146096.0, &dwork2);
	modsd = (long)(sd + dwork1 - dwork2);
	ly=modsd/1461;
	dwork1 = (modsd+1-1461*ly)/365.3;
	if(dwork1 < 0)
		{
		modf(fabs(dwork1), &dwork1);
		dwork1 = -1.0 - dwork1;
		}
	else
		modf(dwork1, &dwork1);
	year=4*ly+dwork1;
	modf((double)(year*365.25), &dwork1);
	doy= modsd + 1 - dwork1;
	year=1601+year;
/* days since 1900, and the unit of time used for AENA formulae, t,
 * the time in Julian centuries of 36525 mean solar days from Greenwich
 * mean noon (12h UT) on 1900 January 0 (AENA p. 528).
 * t is an integer # days + 0.5 (for calculation of quantities at 0h UT
 * on the current date) and tt is the current time.
 */
	days_1900=dsd-109206.0-0.5;
/* 109206.0 = atodate("1899:12:31:0:0:0")/86400.0 */
	t=days_1900/36525.0;
	tt = t + tenths/(864000.0*36525.0);
/* eccentricity of earth's orbit (AENA p. 540) */
	e = 0.01675104 - tt*(4.18e-5 + 1.26e-7*tt);
/* mean obliquity of the ecliptic (AENA p. 540) */
	eps = 0.40931975520272993 - tt*(2.271109689158e-4 +
		tt*(2.86040072e-8 - tt*8.7751276e-9));
/* mean longitude of perigee, mean equinox of date (AENA p. 539) */
	perigee=4.90822946686888692 + tt*( 3.000526416797e-2 +
		tt*( 7.9024630021e-6 + tt*5.81776417e-8));
/* geometric mean longitude, mean equinox of date (AENA p. 539) */
	lo = 4.8815279341118791 + tt*(628.33195094248177 + 5.2796209873e-6*tt);
/* where did all these ^ sig. figures come from?   rob */
/* y, used in calculation of right ascension (ra) */
	y = tan( (double)(eps/2.0) );
	y = y*y;
/* siderial time in tenths of seconds (AENA p. 528) */
	tenths_sid_mean=239258.36+86401845.42*t+
				0.929*t*t+1.0027379093*tenths;
/* mean sideral time within day in radians */
	mean_sid_time = two_pi*modf(tenths_sid_mean/864000.0, &dwork2);
/* hams is the hour angle of the mean sun */
	hams=(tenths-432000.0)*rday; if(  hams < 0  ) hams=hams+two_pi;
/* rams is the right ascension of the mean sun */
	rams=mean_sid_time-hams; if(  rams < 0  ) rams=rams+two_pi;
/* mean anomaly - (c.f. AENA p. 540) */
	g = lo - perigee;
/* eccentric anomaly, mean equinox of date */
	ge = ecanom(e, g);
/* true anomaly, mean equinox of date (SA p. 118) */
	v = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*ge));
/* geometric longitude, mean equinox of date  (Green p 253) */
	clong = v + perigee;
	if(clong > two_pi)clong -= two_pi;
	  else if (clong < 0)clong += two_pi;
/* right ascension, mean equinox of date (SA p. 148) */
	ra = atan2( (1.0-y)*sin(clong),(1+y)*cos(clong));
	if(ra < 0)ra += two_pi;
/* equation of time */
	eot = rams - ra;
	if(eot > pi)eot -= two_pi;
	  else if(eot < -pi)eot += two_pi;
/* declination (SA p. 40) */
	dec = asin(sin(eps)*sin(clong));
/* longitude of ascending node of solar equator on ecliptic (AENA p. 553) */
      omega=1.2857258823024895+2.436188747575e-4*(time-7857648000.0)/3.15576e7;
/* 7,857,648,000 = atodate("1850:1:1:0:0:0") */
	clong_omega=clong-omega;
	if( clong_omega < 0 )
		clong_omega= clong_omega+two_pi;
	  else if( clong_omega >= two_pi)
		clong_omega= clong_omega-two_pi;
	xs=sin(clong_omega); xc= cos(clong_omega);
	cosclo = xc= cos(clong_omega);
/* heliocentric latitude of center of disk, ESAEAENA p. 308   */
	sinb0= xs*sinincl;
	b0 = asin(sinb0);
/* heliocentric longitude of the center of the disk ESAEAENA p. 308
 * currently used to correct older version below.
 */
	l00 = atan2(-xs*cosincl,-xc)+5.1097298952004202 + 0.247564432906997103
	  * (2430000.5 - 2415020.0 - days_1900 - (tenths/864000.0));
	l00 = modf(l00/two_pi, &dwork1);
	l00 *= two_pi;
	if(l00 < 0)l00 += two_pi;
	xx = xc * sinincl;
/*
to within 0.25 degrees, or 4 arc-sec, (i.e. ignore inclination
  of the solar equator) carrington degrees is degrees sun has
  rotated since 1854.0 + degrees before 1854.0 (==720) - degrees
  earth has moved in its orbit since 1854.0 + omega
  calculated in rotations, then converted to degrees
*/
	l0= (time - 7983921600.0)/(25.38*86400.0)
/* 7,983,921,600 = atodate("1854:1:1:12:0:0") */
		 + 2.5 - clong_omega/two_pi - (year-1854);
	if(  doy>90  && clong>omega  ) l0=l0-1;
/* refine to account for the inclination of the sun's rotation axis */
	dwork1 = l0;
	modf(l0, &l0);
	l0 = l0 + (two_pi - l00)/two_pi;
/* correct for possible error due to truncation of l0 at Carr Rot boundary */
	if (dwork1 - l0 > 0.9) l0 = l0 + 1.0;
	else if (l0 - dwork1 > 0.9) l0 = l0 - 1.0;
/* position angle of the northern extremity of the axis of rotation
 * measured eastward from the north point of the disk, ESAEAENA p. 309 
 */
	p=atan(-tan(eps)*cos(clong))+atan(-tanincl*xc);
	xc=cos(clong-perigee);
	r=(1.0-e*e)/(1.0+e*xc);
	rsun=959.63/r;

  ephem[EPHEM_T] = time;
  ephem[EPHEM_TCENT] = t;
  ephem[EPHEM_JD] = 2415020.0 + days_1900;
  ephem[EPHEM_ST] = mean_sid_time;
  ephem[EPHEM_RA] = ra;
  ephem[EPHEM_DEC] = dec;
  ephem[EPHEM_CLONG] = clong;
  ephem[EPHEM_L0] = l0*360;
  ephem[EPHEM_B0] = b0;
  ephem[EPHEM_P] = p;
  ephem[EPHEM_DIST] = r;
  ephem[EPHEM_RSUN] = rsun;
  ephem[EPHEM_VEARTH] = 29785.9*(1+e*xc);
		       /*  1.496e11 * two_pi / 3.155692597474e7 = 29786.31454  */
			  /*  This is consistent with 1 AU = 1.4959792e13 cm  */
  ephem[EPHEM_EOT] = eot;
/*
					original formula prior to 2004.01.29
  ephem[EPHEM_B0V] = xx * ephem[EPHEM_RSUN] / 1.496e11;
		formula in /home/soi/src/libastro.d/solephem.c as of 2004.01.29
  ephem[EPHEM_B0V] = -xx * sin (ephem[EPHEM_RSUN]*pi/(3600*180.0)) * 1.496e11 *
      two_pi / (365.2425*86400);
*/
  ephem[EPHEM_B0V] = xx * ephem[EPHEM_VEARTH];
  ephem[EPHEM_PHI] = atan (0.4336*cos(clong)+p);
  ephem[EPHEM_SINP] = sin (p);
  ephem[EPHEM_COSP] = cos (p);
  ephem[EPHEM_PAROT] = atan (-tanincl*cosclo);
  return;
}

static void calc_sun_ephemeris (TIME time, double ephem[], double obs_lon, double obs_lat) {
  solephem (time - EPHEM_OFFSET, ephem);
  ephem[EPHEM_OBS_LON] = obs_lon * M_PI / 180.0;
  ephem[EPHEM_OBS_LAT] = obs_lat * M_PI / 180.0;
}

