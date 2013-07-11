/*
 * soi_ephem.h 
 *
 *  Defines and declarations for ephemeris functions in libastro.
 *
 *  Responsible:  Rick Bogart                   RBogart@solar.Stanford.EDU
 *
 *  Revision history is at the end of the file.
 */

#ifndef SOI_EPHEM_INCL

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

#include <timeio.h>

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define SOI_EPHEM_VERSION_NUM	(2.9)
#define SOI_EPHEM_INCL	1

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
#define MDIPT_TIME_BEFORE_START	(-1)
#define MDIPT_TIME_BEYOND_END	(1)
#define MDIPT_CANT_INTERPOLATE	(2)
#define MDIPT_CANT_READ_TABLE	(-2)
#define MDIPT_FROM_QUICKLOOK	(3)

#define MDIPT_USE_BEST		(0)
#define MDIPT_USE_QUICKLOOK	(1) /* ftp data from EOF */
#define MDIPT_USE_DISTRIB	(2) /* cdrom data */ 
#define MDIPT_USE_DEFINITIVE	(2) /* GSFC name for cdrom data */
#define MDIPT_USE_SMOOTH	(3) /* smoothed data */
#define MDIPT_USE_DEFAULT	(4) /* default value 511.5 */

#define MDIROLL_TIME_BEFORE_START	(-1)
#define MDIROLL_TIME_BEYOND_END		(1)
#define MDIROLL_CANT_INTERPOLATE	(2)
#define MDIROLL_CANT_READ_TABLE		(-2)
#define MDIROLL_BAD_TABLE		(-3)

#define MDIROLL_USE_DEFAULT		(0) /* rollangle = 0.0 */
#define MDIROLL_USE_CONSTANT		(1) /* rollangle = constant */
#define MDIROLL_USE_SMOOTH		(2) /* "smoothed" table */
#define MDIROLL_USE_RAW			(3) /* raw roll table from CDHF */
#define MDIROLL_USE_TABLE		(4) /* user-specified roll table */

/*
 *  List of observatories and platforms with known location information.  
 *  The values here are longitude west and latitude in degrees.
 */

#define	OBS_LOC_WSO	122.167,36.967	/* Wilcox Solar Obs, Stanford	*/
	    /*  values picked so numbers in solephem will be unchanged  */
#define OBS_LOC_MWO	118.060,34.217	/* Mt Wilson Observatory 	*/
#define OBS_LOC_BBSO	000.000,00.000	/* Big Bear Solar Observatory 	*/
#define OBS_LOC_NSO	110.948,32.233	/* NSO office in Tucson 	*/
#define OBS_LOC_KPNO	111.595,31.958	/* NSO at Kitt Peak		*/
#define OBS_LOC_SPO	 75.357,39.905	/* NSO at Sac Peak		*/
#define OBS_LOC_SPOLE	000.000,-90.00	/* South Pole			*/
#define OBS_LOC_CRAO	000.000,00.000	/* Crimean Astrophys Obs	*/

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/
  

/****************************************************************************/
/***************************  FUNCTION PROTOTYPES  **************************/
/****************************************************************************/

/*  libastro  */
							/*  mdi_roll.c  */
extern int mdi_roll (char *rolltbl, TIME tobs, double *roll, TIME *tmod);
							/*  mdi_pointing.c  */
extern int mdi_point (TIME obs_time, double *xc, double *yc, TIME *mod_time,
    int source);
extern int mdi_pointing (TIME obs_time, double *xc, double *yc, TIME *mod_time);
						/*  sun_ephemeris.c  */
extern int soho_ephemeris (TIME obs_time, double *r, double *lat, double *lon,
    double *vr, double *vn, double *vw, TIME *mod_time);
extern double soho_obs_dist (TIME t);
extern double soho_latitude (TIME t);
extern double soho_longitude (TIME t);
extern void sun_ephemeris (TIME time, double *ephem, double longi, double lat);
extern double sun_delta_v (TIME time, double *ephem, double x, double y);
extern double carrington_rots (TIME time);
extern TIME SOHO_meridian_crossing (double lon, int car_rot);
extern TIME earth_meridian_crossing (double lon, int car_rot);
extern TIME meridian_crossing (double lon, int car_rot);

/* old functions, should be reserved for WSO programs */

extern void solephem (double time, double *ephem);
extern double delta_v (double time, double *ephem, double x, double y);


#endif
/*
 *  Revision History
 *  96.08.08	Rick Bogart	added prototypes for mdi_pointing()
	and soho_ephemeris(); repaired declarations for sun_ephemeris()
	and sun_delta_v();    
 *  96.09.11	R Bogart	added table mod time argument to
	mdi_pointing() and soho_ephemeris()
 *  97.01.03	R Bogart	added declaration for carrington_rots()
 *  97.04.21	R Bogart	add declaration for meridian_crossing()
 *  98.02.24	R Bogart	new prototype of mdi_pointing and defs
 *	to support multiple data sources
 *  98.03.13	Kay Leibrand	added defines for MDIPT_USE_DISTRIB and
	MDIPT_USE_SMOOTH, making MDIPT_USE_DISTRIB equivalent to
	MDIPT_USE_DEFINITIVE
 *  98.11.23	Kecheng Chu	added MDIROLL_*
 *  01.02.28	R Bogart	added helper functions for aspects of
 	soho_ephemeris
 *  01.05.07	K Chu		added MDIPT_USE_DEFAULT
 *  04.01.29	R Bogart	added one more define for solephem
 *  06.08.10	R Bogart	added declarations for additional meridian
 *	crossing functions; changed soi_time.h to timeio.h
 */

/*
$Id: soi_ephem.h,v 1.1 2013/07/11 15:13:59 yliu Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/soi_ephem.h,v $
$Author: yliu $
 * $Log: soi_ephem.h,v $
 * Revision 1.1  2013/07/11 15:13:59  yliu
 * new module for producing synoptic charts for MDI line-of-sight mags
 *
 * Revision 1.16  2007/05/04  21:18:23  rick
 * see above
 *
 * Revision 1.15  2004/01/30  01:15:45  rick
 * see above
 *
 * Revision 1.14  2001/05/07 22:07:07  kehcheng
 * added MDIPT_USE_DEFAULT
 *
 * Revision 1.13  2001/02/28 22:35:34  rick
 * added helper functions for aspects of soho_ephemeris
 *
 * Revision 1.12  1998/11/23  18:41:20  kehcheng
 * added MDIROLL_*
 *
 * Revision 1.11  1998/03/13 23:05:21  kay
 * added defines for MDIPT_USE_DISTRIB and MDIPT_USE_SMOOTH
 * making MDIPT_USE_DISTRIB equivalent to MDIPT_USE_DEFINITIVE
 *
 * Revision 1.10  1998/02/25  19:31:49  rick
 * see above
 *
 * Revision 1.9  1997/05/01  22:43:50  rick
 * see above
 *
 * Revision 1.8  1997/04/16  21:51:36  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.7  1997/01/13  18:43:06  rick
 * see above
 *
 * Revision 1.6  1996/09/12  09:48:40  rick
 * see above
 *
 * Revision 1.5  1996/08/15  15:41:14  CM
 * change stdtime.h to soi_time.h
 *
 * Revision 1.4  1996/08/15  00:18:43  CM
 * added include <stdtime.h>
 *
 * Revision 1.3  1996/08/14  20:55:55  rick
 * see above
 *
 * Revision 1.2  1995/07/13  22:31:04  phil
 * *** empty log message ***
 * */
