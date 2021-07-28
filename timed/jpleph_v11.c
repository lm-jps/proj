/*****************************************************************************
*        *****    jpl planetary and lunar ephemerides    *****       ver.1.1 *
*                 originally modified: 4 March 2004 by RSB                   *
*****************************************************************************/
/*
 *  jpleph.c - function pleph() for reading from JPL binary ephemeris table
 *    and returning relative position-velocity coordinates for selected pair
 *    of solar system objects
 *
 *  Author: Richard S. Bogart (rbogart@spd.aas.org), Stanford University
 *
 *  This code is based on the JPL routines included in the file testeph.f
 *    available at ftp://ssd.jpl.nasa.gov/pub/eph/export/fortran/
 *    (last modified 23 July 2001) and translated into C by
 *    Piotr A. Dybczynski (dybol@phys.amu.edu.pl), version 1.2 (23 July 1997)
 *    available at ftp://ftp.astro.amu.edu.pl/pub/jpleph (see also
 *    http://www.projectpluto.com/jpl_eph.htm
 *
 *  The following substantive changes have been made to the original code:
 *
 *  Version 0.9 (frozen 2007.10.31):
 *    An internal byte-reordering function has been added in order to
 *	support reading of the big-endian JPL-distributed tables on
 *	little-endian platforms that support the IEEE-754 numeric data
 *	representations (notably Intel processors).  This allows the
 *	sharing of the same binary tables by machines with multiple
 *	architectures and obviates the necessity of doing private ASCII
 *	to binary table conversions.  The type of platform is automatically
 *	detected by comparing redundant binary and ASCII information in the
 *	JPL tables.
 *    The type of the JPL table is automatically detected and the read buffer
 *	sizes and record lengths allocated appropriately.
 *    An additional function pleph_ext() has been provided with added arguments
 *	for supporting options that were previously fixed on compilation:
 *	  - the option for returning the results in km and km/sec rather than
 *	    AU and AU/day
 *	  - specification of the filename of the binary ephemeris table
 *	The basic function pleph() uses the (as-distributed) pre-compiled
 *	values for these options.
 *    The compiler option for non-barycentric values was removed, as it was
 *	in fact unconditionally overridden.
 *    Removed option for positions only from interp()
 *    Moved checks for nutations and librations in pleph_ext() to follow
 *	state() call, where they must logically be in order to guarantee
 *	valid values
 *  Version 1.0 (frozen 2009.12.06):
 *    Modified et to split integer and fractional parts for allegedly greater
 *	interpolation precision, though it does not seem to make much difference
 *    Added static variables to (a) ensure that the earth-moon mass ratio
 *	is properly set when there are repeated calls to pleph_ext within
 *	a program; and (b) allow for reading from different tables in calls
 *	from the same calling program
 *    Replaced several initialization loops with memset and memcpy calls
 *  Version 1.1 (created 2021.07.27):
 *    Replaced explicit path for default file to define created by (and
 *	dependent on) DRMS site localization; in the absence of that define,
 *	an error message from state() is guaranteed (unless the default
 *	table happens to be in the running directory!)
 *
 */
#define JPL_EPHEMERIS_VERSION_NUM     (1.1)

#define DEFAULT_EPHEM_FILE  "unxp2000.405"
/*
#define DEFAULT_EPHEM_FILE  "/home/rick/src/ephem/tables/unxp2000.405"
#define DEFAULT_EPHEM_FILE  "/home/soi/CM/tables/jpleph.403.1950-2050"
*/
#include <stdio.h>
#include <math.h>

/*
 *  Byte reordering function added to support reading of big-endian binary
 *    tables on little-endian architectures that still use the IEEE-754
 *    floating point representation
 */
static void bytereorder (void *data, int size, int num) {
  int n, s;
  unsigned char *c = (unsigned char *)data;
  unsigned char t;

  for (n=0; n<num; n++, c+=size)
    for (s=0; s<(size/2); s++) {
    t = *(c + s);
    *(c + s) = *(c + size - s - 1);
    *(c + size - s - 1) = t;
  }
}

/****************************************************************************
****                       split(tt,fr)                                  ****
*****************************************************************************
****  this subroutine breaks a d.p. number into a d.p. integer           ****
****  and a d.p. fractional part.                                        ****
****                                                                     ****
****  calling sequence parameters:                                       ****
****                                                                     ****
****    tt = d.p. input number                                           ****
****                                                                     ****
****    fr = d.p. 2-word output array.                                   ****
****         fr(1) contains integer part                                 ****
****         fr(2) contains fractional part                              ****
****                                                                     ****
****         for negative input numbers, fr(1) contains the next         ****
****         more negative integer; fr(2) contains a positive fraction.  ****
****************************************************************************/
static void split (double tt, double fr[2]) {
  fr[1] = modf (tt, &fr[0]);
  if (tt >= 0.0 || fr[1] == 0.0) return;
			    /*  make adjustments for negative input number  */
  fr[0] = fr[0] - 1.0;
  fr[1] = fr[1] + 1.0;
  return;
}

/*****************************************************************************
**                     interp(buf,t,ncf,ncm,na,pv)                          **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**        buf   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**         tt   fractional time in interval covered by coefficients at      **
**              which interpolation is wanted  (0 <= tt <= 1)               **
**                                                                          **
**	intrvl	length of whole interval in input time units                **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**      output:                                                             **
**                                                                          **
**        pv   interpolated quantities requested.  dimension                **
**              expected is pv(ncm,ifl), dp.                                **
**                                                                          **
*****************************************************************************/

static void interp (double *coef, double tt, double intrvl, int ncf, int ncm,
    int na, double posvel[6]) {
  static double pc[18], vc[18];
  static double twot = 0.0;
  static int first_call = 1, np = 2, nv = 3;
  double dna, dt1, tc, temp, temp1, vfac;
  int i, j, l, lcfcm, offset;

  if (first_call ) {
					     /*  initialize static vectors  */
    pc[0] = 1.0;
    pc[1] = 0.0;
    vc[1] = 1.0;
    first_call = 0;
  }
/*
    entry point. get correct sub-interval number for this set
    of coefficients and then get normalized chebyshev time
    within that subinterval.
*/
  dna = (double)na;
  modf (tt, &dt1);
  temp = dna * tt;
  l = (int)(temp - dt1);
  lcfcm = l * ncf * ncm;
		   /*  tc is the normalized chebyshev time (-1 <= tc <= 1)  */
  tc = 2.0 * (modf (temp, &temp1) + dt1) - 1.0;
/*
    check to see whether chebyshev time has changed,
    and compute new polynomial values if it has.
    (the element pc[1] is the value of t1[tc] and hence
    contains the value of tc on the previous call.)
*/
  if (tc != pc[1]) {
    np = 2;
    nv = 3;
    pc[1] = tc;
    twot = tc + tc;
  }
/*
    be sure that at least 'ncf' polynomials have been evaluated
    and are stored in the array 'pc'.
*/
  if (np < ncf) {
    for (i=np; i<ncf; i++)  pc[i] = twot*pc[i-1] - pc[i-2];
    np = ncf;
  }
			/*  interpolate to get position for each component  */
  for (i=0; i<ncm; i++){		  /* ncm is a number of coordinates */
    posvel[i] = 0.0;
    offset = i*ncf + lcfcm;
    for (j=ncf-1; j>=0; j--)
      posvel[i] += pc[j] * coef[j + offset];
  }
  /*  be sure enough derivative polynomials have been generated and stored  */
					    /*  for velocity interpolation  */
  vfac = (dna + dna) / intrvl;
  vc[2] = twot + twot;
  if (nv < ncf) {
    for (i=nv; i<ncf; i++) vc[i] = twot*vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
    nv = ncf;
  }
			/*  interpolate to get velocity for each component  */
  for (i=0; i<ncm; i++) {
    posvel[i+ncm] = 0.0;
    offset = i*ncf + lcfcm;
    for (j=ncf-1; j>0; j--)
      posvel[i+ncm] += vc[j] * coef[j + offset];
    posvel[i+ncm] *= vfac;
  }
  return;
}

/*****************************************************************************
**  static int state (double *et2, int *list, char *table, int in_kms,
    double pv[][6], double *nut, double *emrat, double pvsun[6],
    int *nut_flag, int *lib_flag)
******************************************************************************
** This subroutine reads and interpolates the jpl planetary ephemeris file  **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**    Input:                                                                **
**                                                                          **
**        et2[] double, 2-element JED epoch at which interpolation          **
**              is wanted.  Any combination of et2[0]+et2[1] which falls    **
**              within the time span on the file is a permissible epoch.    **
**                                                                          **
**               a. for ease in programming, the user may put the           **
**                  entire epoch in et2[0] and set et2[1]=0.0               **
**                                                                          **
**               b. for maximum interpolation accuracy, set et2[0] =        **
**                  the most recent midnight at or before interpolation     **
**                  epoch and set et2[1] = fractional part of a day         **
**                  elapsed between et2[0] and epoch.                       **
**                                                                          **
**               c. as an alternative, it may prove convenient to set       **
**                  et2[0] = some fixed epoch, such as start of integration,**
**                  and et2[1] = elapsed interval between then and epoch.   **
**                                                                          **
**       list   12-element integer array specifying what interpolation      **
**              is wanted for each of the "bodies" on the file.             **
**                                                                          **
**                        list[i]=0, no interpolation for body i            **
**                               !=0, position and velocity                 **
**                                                                          **
**              The designation of the astronomical bodies by i is:         **
**                                                                          **
**                        i = 0: Mercury                                    **
**                          = 1: Venus                                      **
**                          = 2: Earth-Moon barycenter                      **
**                          = 3: Mars                                       **
**                          = 4: Jupiter                                    **
**                          = 5: Saturn                                     **
**                          = 6: Uranus                                     **
**                          = 7: Neptune                                    **
**                          = 8: Pluto                                      **
**                          = 9: geocentric Moon                            **
**                          =10: nutations in longitude and obliquity       **
**                          =11: lunar librations (if on file)              **
**                                                                          **
 *	table	character string containing name of file containing binary
 *		ephemeris data
 *
 *	in_kms	integer flag: if non-zero, values returned in km and km/s
 *		otherwise AU and AU/day
 *
**    output:                                                               **
**                                                                          **
**	pv[][6]	double array that will contain requested interpolated       **
**              quantities.  The body specified by list[i] will have its    **
**              state in the array starting at pv[i][0]  (on any given      **
**              call, only those words in 'pv' which are affected by the    **
**              first 10 'list' entries (and by list[11] if librations are  **
**              on the file) are set.  The rest of the 'pv' array           **
**              is untouched.)  The order of components in pv[][] is:       **
**              pv[][0]=x,pv[][1]=y,....pv[][5]=dz.                                   **
**                                                                          **
**              All output vectors are referenced to the Earth mean         **
**              equator and equinox of epoch. The Moon state is always      **
**              geocentric; the other nine states are either heliocentric   **
**              or solar-system barycentric, depending on the setting of    **
**              global variables (see below).                               **
**                                                                          **
**              Lunar librations, if on file, are put into pv[10][k] if     **
**              list[11] is 1 or 2.                                         **
**                                                                          **
**	nut	dp 4-element array that will contain nutations and rates,   **
**              depending on the setting of list[10].  the order of         **
**              quantities in nut is:                                       **
**                                                                          **
**                       d psi  (nutation in longitude)                     **
**                       d epsilon (nutation in obliquity)                  **
**                       d psi dot                                          **
**                       d epsilon dot                                      **
 *
 *	emrat	pointer to variable containing earth-moon mass ratio
 *
 *    pvsun[6]	double 6-element array containing the barycentric position
 *	and velocity of the Sun
 *
 *  nut_flag	flag value telling whether table contains nutations
 *
 *  lib_flag	flag value telling whether table contains librations
 *
 ****************************************************************************/
static int state (double *et2, int *list, char *table, int in_kms,
    double pv[][6], double *nut, double *emrat, double pvsun[6],
    int *nut_flag, int *lib_flag) {
  static FILE *binfile;
  static double *buf;
  static double val[400];
  static double au, tstart, tstop, tstep;
  static double earth_moon_mass_ratio;
  double pefau[6], pjd[4], ss[3];
  double aufac, jda0, jda1, s, t, tstp;
  static int ipt[13][3];
  static int first_call = 1, reorder = 0, recl = 0;
  static int coeffs, recsize;
  int i, j, n, ncon, numde, rec, series_num;
  static char *last_table = NULL;
  char title[3][84], name[400][6];

  if (last_table) {
    if (strcmp (table, last_table)) {
      free (last_table);
      fclose (binfile);
      first_call = 1;
    }
  }
  if (first_call) {
    if (sizeof (double) != 8) {
      fprintf (stderr, "** Error: Local architecture does not use 8-byte DP floats\n");
      return (1);
    }
    if (!(binfile = fopen (table, "r"))) {
      fprintf (stderr, "** Error: Unable to open file \"%s\" for read\n", table);
      return (1);
    }
    for (n = 0; n < 3; n++) fread (title[n], 84, 1, binfile);
    if (sscanf (title[0], "JPL Planetary Ephemeris DE%d/LE", &series_num) != 1) {
      fprintf (stderr, "** Error: Unexpected title format in file\n  %s\n",
          table);
      return (1);
    }
    if (sscanf (title[1], "Start Epoch: JED= %lf", &jda0) != 1) {
      fprintf (stderr, "** Error: Unexpected title format in file\n  %s\n",
          table);
      return (1);
    }
    if (sscanf (title[2], "Final Epoch: JED= %lf", &jda1) != 1) {
      fprintf (stderr, "** Error: Unexpected title format in file\n  %s\n",
          table);
      return (1);
    }
    if (series_num == 200 || series_num == 202) coeffs = 826;
    else if (series_num == 403 || series_num == 405)  coeffs = 1018;
    else if (series_num == 404 || series_num == 406)  coeffs = 728;
    else {
      fprintf (stderr, "** Error: Unknown ephemeris series %d in file\n  %s\\n",
          series_num, table);
      return (1);
    }
    recsize = coeffs * sizeof (double);
    buf = (double *)malloc (recsize);
						/*  May not need these  */
    for (n = 0; n < 400; n++) fread (name[n], 6, 1, binfile);

    fread (ss, sizeof (ss), 1, binfile);
    if ((fabs (ss[0] - jda0) > 0.5) || (fabs (ss[1] - jda1) > 0.5)) {
      bytereorder (ss, sizeof (double), 3);
      if ((fabs (ss[0] - jda0) > 0.5) || (fabs (ss[1] - jda1) > 0.5)) {
        fprintf (stderr, "Error: Local architecture does not support IEEE-754\n");
        fprintf (stderr, "       or binary epoch limits inconsistent with header\n");
        return (1);
      }
      reorder = 1;
    }
    tstart = ss[0];
    tstop = ss[1];
    tstep = ss[2];
				     /*  ncon is 144 in DE403, 156 in DE405  */
    fread (&ncon, sizeof (int), 1, binfile);
    fread (&au, sizeof (double), 1, binfile);
    if (reorder) bytereorder (&au, sizeof (double), 1);
    fread (emrat, sizeof (double), 1, binfile);
    if (reorder) bytereorder (emrat, sizeof (double), 1);
    earth_moon_mass_ratio = *emrat;
    for (n = 0; n < 12; n++) fread (ipt[n], sizeof (int), 3, binfile);
						/*  Ignore?  */
    fread (&numde, sizeof (int), 1, binfile);
    fread (ipt[12], sizeof (int), 3, binfile);
    if (reorder) bytereorder (ipt, sizeof (int), 39);
    *nut_flag = ipt[11][1];
    *lib_flag = ipt[12][1];
    fseek (binfile, recsize, SEEK_SET);
						/*  May not need these  */
    fread (val, sizeof (val), 1, binfile);
    if (reorder) bytereorder (val, sizeof (double), 400);
    last_table = malloc (strlen (table) + 1);
    strcpy (last_table, table);
    first_call = 0;
  }

/****************************** main entry point *****************************/

  if (et2[0] == 0.0) return 0;
  s = et2[0] - 0.5;
  split (s, &pjd[0]);
  split (et2[1], &pjd[2]);
  pjd[0] = pjd[0] + pjd[2] + 0.5;
  pjd[1] = pjd[1] + pjd[3];
  split (pjd[1], &pjd[2]);
  pjd[0] = pjd[0] + pjd[2];

/* here pjd[0] contains last midnight before epoch desired (in JED: *.5)
   and pjd[3] contains the remaining, fractional part of the epoch         */
				   /*  error return for epoch out of range  */
  if ((pjd[0] + pjd[3]) < tstart || (pjd[0] + pjd[3]) > tstop) {
    fprintf (stderr, "Requested JED not within ephemeris limits:\n");
    fprintf (stderr, "%f < %f || > %f\n", pjd[0] + pjd[3], tstart, tstop);
    return 1;
  }
  *emrat = earth_moon_mass_ratio;
		      /*  Calculate record # and relative time in interval  */
      /*  add 2 to adjust for the first two records containing header data  */
  rec = (long)((pjd[0] - tstart) / tstep) + 2;
  if (pjd[0] == tstop) rec--;
  t = (pjd[0] - ((rec - 2) * tstep + tstart) + pjd[3]) / tstep;
	      /*  read correct record if not in core (static vector buf[])  */
  if (rec != recl) {
    if (fseek (binfile, rec*recsize, SEEK_SET)) {
      fprintf (stderr, "Error seeking to new record within file \"%s\"\n",
	  table);
      return 1;
    }
    fread (buf, sizeof (*buf), coeffs, binfile);
    if (reorder) bytereorder (buf, sizeof (double), coeffs);
    recl = rec;
  }
  tstp = tstep;
  aufac = 1.0;
  if (in_kms) tstp *= 86400.0;
  else aufac /= au;
			      /*  interpolate Solar System barycentric Sun  */
  interp (&buf[ipt[10][0]-1], t, tstp, ipt[10][1], 3, ipt[10][2], pefau);
  for (n=0; n<6; n++) pvsun[n] = pefau[n] * aufac;
		  /*  check and interpolate whichever bodies are requested  */
  for (i=0; i<10; i++) {
    if (list[i] == 0) continue;
    interp (&buf[ipt[i][0]-1], t, tstp, ipt[i][1], 3, ipt[i][2], pefau);
    for (j=0; j<6; j++) pv[i][j] = pefau[j]*aufac;
  }
			    /*  do nutations if requested (and if on file)  */
  if (list[10] && ipt[11][1] > 0)
    interp (&buf[ipt[11][0]-1], t, tstp, ipt[11][1], 2, ipt[11][2], nut);
			  /*  get librations if requested (and if on file)  */
  if( list[11] && ipt[12][1] > 0) {
    interp (&buf[ipt[12][0]-1], t, tstp, ipt[12][1], 3, ipt[12][2], pefau);
    for (j=0; j<6; j++) pv[10][j] = pefau[j];
  }
  return 0;
}

/*
 *  pleph_ext (double jed, int target, int center, int in_kms, double *rrd
 *    char *table)
 *
 *  Added function to provide extended option capability via additional
 *    arguments; the traditional function pleph merely calls this function
 *    with the extra arguments set to their (as-distributed) pre-compiled
 *    values.
 *
 *  Parameters:
 *
 *    jed = julian ephemeris date at which interpolation is wanted
 *
 *    target = code number of 'target' point
 *
 *    center = code number of center point
 *
**    The numbering convention for 'target' and 'center' is:                **
**                                                                          **
**            1 = mercury           8 = neptune                             **
**            2 = venus             9 = pluto                               **
**            3 = earth            10 = moon                                **
**            4 = mars             11 = sun                                 **
**            5 = jupiter          12 = solar-system barycenter             **
**            6 = saturn           13 = earth-moon barycenter               **
**            7 = uranus           14 = nutations (longitude and obliq)     **
**                                 15 = librations, if on eph. file         **
**                                                                          **
**            (If nutations are wanted, set target = 14.                    **
**             For librations, set target = 15. set center= 0)              **
**
 *    in_kms	units for position velocity array rrd
 *		0 -> AU, AU/day; else km, km/sec
 *
 *	rrd	output 6-element, double array of position and velocity
 *		of point 'target' relative to 'center', in units requested.
 *		For librations the units are radians and radians per day.
 *		In the case of nutations the first four words of
 *		rrd will be set to nutations and rates, having units of
 *		radians and radians/day.
 *
 */
void pleph_ext (double jed, int target, int center, int in_kms, double *rrd,
    char *table) {
  double et2[2], pv[13][6];		/* pv is the position/velocity array
				NUMBERED FROM ZERO: 0:Mercury,1:Venus,...
				8=Pluto,9=Moon,10=Sun,11=SSBary,12=EMBary
			First 10 elements (0-9) are affected by state(),
			all are adjusted here.				     */
  double pvsun[6];

  double emrat;
  int i, k, nflag, lflag;
  int list[12];          /* list is a vector denoting, for which "body"
                            ephemeris values should be calculated by state():
                            0=Mercury,1=Venus,2=EMBary,...,8=Pluto,
                            9=geocentric Moon, 10=nutations in long. & obliq.
                            11= lunar librations  */

		 /*  initialize et2 for 'state' and set up component count   */
  et2[0] = trunc (jed);
  et2[1] = fmod (jed, 1.0);
  memset (rrd, 0, 6 * sizeof (double));
/*
  for (i=0; i<6; i++) rrd[i] = 0.0;
*/

  if (target == center) return;

  memset (list, 0, 12 * sizeof (int));
/*
  for (i=0; i<12; ++i) list[i] = 0;
*/
		   /*  set up proper entries in 'list' array for state call  */
  for (i=0; i<2; i++) {
    k = (i) ? center - 1 : target - 1;
    if (k <= 9) list[k] = 2;				  /*  Major planets  */
    if (k == 9) list[2] = 2;	/*  for moon state earth state is necessary  */
    if (k == 2) list[9] = 2;	/*  for earth state moon state is necessary  */
    if (k == 12) list[2] = 2;		        /*  EMBary state additional  */
  }
						     /*  make call to state  */
  state (et2, list, table, in_kms, pv, rrd, &emrat, pvsun, &nflag, &lflag);
						    /*  check for nutations  */
  if (target == 14) {
    if (nflag > 0) {
      list[10] = 2;
      state (et2, list, table, in_kms, pv, rrd, &emrat, pvsun, &nflag, &lflag);
    }
    else fprintf (stderr, "***** no nutations on the ephemeris file  ******\n");
    return;
  }
						   /*  check for librations  */
  if (target == 15) {
    if (lflag > 0) {
				 /*  there are librations on ephemeris file  */
      list[11] = 2;
      state (et2, list, table, in_kms, pv, rrd, &emrat, pvsun, &nflag, &lflag);
      memcpy (rrd, pv[10], 6 * sizeof (double));
/*
      for (i=0; i<6; i++) rrd[i] = pv[10][i];
*/
    }
    else fprintf (stderr, "***** no librations on the ephemeris file  ******\n");
    return;
  }
		    /*  Solar System barycentric Sun state goes to pv[10][]  */
  if (target == 11 || center == 11)
    memcpy (pv[10], pvsun, 6 * sizeof (double));
/*
    for (i=0; i<6; i++) pv[10][i] = pvsun[i];
*/
	 /*  Solar System Barycenter coordinates & velocities equal to zero  */
  if (target == 12 || center == 12)
    memset (pv[11], 0, 6 * sizeof (double));
/*
    for (i=0; i<6; i++) pv[11][i] = 0.0;
*/
				 /*  Solar System barycentric EMBary state:  */
  if (target == 13 || center == 13)
    memcpy (pv[12], pv[2], 6 * sizeof (double));
/*
    for (i=0; i<6; i++) pv[12][i] = pv[2][i];
*/
			    /*  if Moon from Earth or Earth from Moon .....  */
  if ((target * center) == 30 && (target + center) == 13)
    memset (pv[2], 0, 6 * sizeof (double));
/*
    for(i=0; i<6; i++) pv[2][i] = 0.0;
*/
  else {
    if (list[2] == 2)		      /*  calculate earth state from EMBary  */
      for (i=0; i<6; i++) pv[2][i] -= pv[9][i] / (1.0 + emrat);
    if (list[9] == 2)	  /*  calculate Solar System barycentric Moon state  */
      for (i=0; i<6; i++) pv[9][i] += pv[2][i];
  }

  for (i=0; i<6; i++) rrd[i] = pv[target-1][i] - pv[center-1][i];
  return;
}

/*****************************************************************************
**    pleph (jed, target, center, rrd)                                      **
**                                                                          **
**    This subroutine reads the jpl planetary ephemeris                     **
**    and gives the position and velocity of the point 'target'             **
**    with respect to 'center'.                                             **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**    jed = (double) julian ephemeris date at which interpolation           **
**           is wanted.                                                     **
**                                                                          **
**    target = integer number of 'target' point.                            **
**                                                                          **
**    center = integer number of center point.                              **
**                                                                          **
**    The numbering convention for 'target' and 'center' is:                **
**                                                                          **
**            1 = mercury           8 = neptune                             **
**            2 = venus             9 = pluto                               **
**            3 = earth            10 = moon                                **
**            4 = mars             11 = sun                                 **
**            5 = jupiter          12 = solar-system barycenter             **
**            6 = saturn           13 = earth-moon barycenter               **
**            7 = uranus           14 = nutations (longitude and obliq)     **
**                                 15 = librations, if on eph. file         **
**                                                                          **
**            (If nutations are wanted, set target = 14.                    **
**             For librations, set target = 15. set center= 0)              **
**                                                                          **
**     rrd = output 6-element, double array of position and velocity        **
**           of point 'target' relative to 'center'. The units are au and   **
**           au/day. For librations the units are radians and radians       **
**           per day. In the case of nutations the first four words of      **
**           rrd will be set to nutations and rates, having units of        **
**           radians and radians/day.                                       **
**                                                                          **
*****************************************************************************/
void pleph (double jed, int target, int center, double *rrd) {
  char deftbloc[DRMS_MAXPATHLEN];
#ifdef JPL_EPHEM_TABLEDIR
  sprintf (deftbloc, "%s/%s", JPL_EPHEM_TABLEDIR, DEFAULT_EPHEM_FILE);
#else
  sprintf (deftbloc, "%s", DEFAULT_EPHEM_FILE);
#endif
  pleph_ext (jed, target, center, 0, rrd, deftbloc);
}
