/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/stdtime.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: stdtime.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.6  90/10/23  13:15:57  phil
 * ascii_time and ascii_now added
 * 
 * Revision 1.5  90/10/23  12:16:13  phil
 * timezone added
 * 
 * Revision 1.4  90/10/23  11:20:03  phil
 * add _Mon declarations
 * 
 * Revision 1.3  90/10/22  14:45:12  phil
 * include timeparam()
 * 
 * Revision 1.2  90/09/20  12:03:24  phil
 * declare sprintt
 * 
 * Revision 1.1  90/09/17  16:40:53  phil
 * Initial revision
 * 
*/

#ifndef STDTIME
#define STDTIME 1
#define SOI_TIME_INCL   1

#ifndef TIME_DEFINED
typedef double TIME;
#define TIME_DEFINED
#endif

extern TIME timeparam();
extern int fprintinc(), printinc();
extern TIME date(), stdtime(), ymdtime(), bartime(), atodate();
extern TIME timenow(), atoinc();
extern TIME carrtime(), Ctime();
extern char *ascii_time(), *ascii_now();
extern char *printt(), *fprintt(), *sprintt(); 
extern char *sprintinc(), *getinc();
extern char *signif2(), *signif(), *signiffrac();
extern char *sindex(), *stindex(), *strlow(), *strup(), *mprefix();
extern char *upper(), *lower();
extern double divint();
extern double atofreq();
extern double floor2(), ceil2();
extern double atocarr();
extern double angle_norm();
extern double ask(), askt(), askct(), askinc(), askfr(), question();
extern char * asks();
extern double round();

#ifndef __linux
extern char *strdup();
#endif

typedef struct
	{TIME	_time;		/* standard time		*/
	 short	_year;		/* year as 1980			*/
	 short	_doy;		/* day of year: Jan 1 = 1	*/
	 short	_month;		/* month: Jan = 1		*/
	 short	_day;		/* day in month			*/
	 short	_hour;		/* hour in day: 0 to 23		*/
	 short	_minute;	/* minute in hour: 0 to 59	*/
	 short	_second;	/* second: 0 to 59		*/
	 short	_dow;		/* day of week: Monday = 1	*/
	 short	_dim;		/* days in current month	*/
	 short	_leap;		/* 1 if a leap year		*/
	 short	_spare1;	/* for later implimentation	*/
	 short	_spare2;	/* for later implimentation	*/
	 short	_brot;		/* Bartels rotation		*/
	 short	_bday;		/* day within rotation		*/
	 short	_Crot;		/* Carrington rotation		*/
	 short	_Clong;		/* Carrington longitude		*/
	 TIME	_Ctime;		/* Carrington time		*/
	 TIME	_SPARETIME;	/* What we all have so much of	*/
	 char	_tzone[4];	/* TimeZone from timenow()	*/
	} DATE;

extern DATE curtime;

static char * _Months[13] =
	{
	"BAD",
	"January", "February", "March", "April", "May", "June",
	"July", "August", "September", "October", "November", "December"
	} ;

static char * _MON[13] =
	{
	"BAD",
	"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
	"JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
	} ;

static char * _Mon[13] =
	{
	"BAD",
	"Jan", "Feb", "Mar", "Apr", "May", "Jun",
	"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
	} ;
#endif
