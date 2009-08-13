/*
 * time_convert - change TAI time since 1977 into ascii with a given timezone
 * and from an ascii time to seconds.  Takes one argument of time or seconds and an optional
 * flag argument of time zone.  Default is UTC for conversion to ascii.
 */

/**
   @defgroup time_convert time_convert
   @ingroup su_util

   @brief Convert among the internal DRMS time representation and other time representations.

   @par Synopsis:
   @code
   time_convert -h
   time_convert s=<secondsJSOC> | sdo=<secondsSDO> | egse=<secondsEGSE> | time=<calenderTime> |
                ord=<ordinalDate> [o=jsoc | o=sdo | o=egse | o=ord | o=cal] [zone=<zone>] [p=<precision>]

   @endcode

   Converts internal DRMS time (seconds since 15 seconds 
   before January 1, 1977 UTC) to an external ascii string 
   representation. The exact external representation depends on the command-line
   arguments provided to this utility. Possiblities include: <tt>SDO time</tt>, 
   <tt>EGSE time</tt>, <tt>Ordinal date</tt>, and <tt>Calendar time</tt>.  
   <tt>SDO time</tt> is the number of seconds that have elapsed since 
   January 1, 1958 TAI, i.e. SDO onboard time. The full format is the string
   representation of a double data type. <tt>EGSE time</tt> is 
   the number of seconds that have elapsed since @e APPROXIMATELY January 1, 2004 UTC
   (actual epoch is 2003.12.30_23:59:36.000_UTC). The full format is also the string
   representation of a double data type. <tt>Ordinal date</tt> is 
   the day number of the year (starting at day 1 on January 1). The full format
   is YYYY.DDD[_ZZZ], where YYYY is the year, DDD is the day number, and ZZZ is the
   zone (e.g., UTC, TDT, PDT, etc.). <tt>Calendar time</tt> gives the year, month, date, hour, 
   minutes, and seconds for a given system of time (like @e UTC). The full format is
   specified in <tt>JSOC TN 07-001</tt> (http://jsoc.stanford.edu/doc/timerep.html), 
   but in short it looks like:

   @par
   1.  year.month.fracday[_type]
   @par
   2.  year.month.day_hour:minute[:second][_type]
   @par
   3.  {MJD|JD}_julday[_type]

   where @a type refers to the time system or time zone (e.g., @e UTC, or @e PST), @e MJD and
   @e JD refer to Modified Julian Date and Julian Date, and @e julday refers to a Julian 
   day number.

   The precision of the seconds field is specified with the <tt>p</tt> parameter.  The default is 0.  Setting p=3 will produce output identical to the original version of time_convert.

   Alternatively, and with the appropriate command-line parameters, @ref time_convert 
   converts from any supported time representation to any other representation.

   If multiple input format strings are specified as arguments 
   (e.g., time_convert s=234235235.35 ord=1982.035), only one will be used. A descriptive 
   string will be printed describing which input string was used. Given a definitive
   input format, @ref time_convert chooses a default output format. If the input is 
   an internal time, an SDO time, an EGSE time, or an ordinal date, the default output 
   is a calendar time. If the definitive input format is a calendar time, the default
   output is an internal time. The default output format can be overwritten
   by providing the @a o=format argument. 

   When the output format is either an ordinal time or a calendar time, the 
   time can be expressed in any of several supported time systems (e.g., @e UTC,
   @e TDT, @e TAI, or even a time zone, like @e PDT). This is accomplished by
   supplying the appropriate @a zone=system argument. Refer to 
   <tt>JSOC TN 07-001</tt> (http://jsoc.stanford.edu/doc/timerep.html) for a 
   complete list of the supported "zones".

   The default time format is with no fractions of seconds.  If o=jsoc is specified the seconds are provided to the nearest ms.

   @par Flags:
   @c -h: Show usage message.

   @param s An input time formatted as an internal time. 
   <secondsJSOC> is seconds since 15 seconds before January 1, 1977 UTC.
   @param sdo  An input time formatted as an SDO time.
   <secondsSDO> is seconds since January 1, 1958 TAI.
   @param egse  An input time formatted as an EGSE time.
   <secondsEGSE> is seconds since 2003.12.30_23:59:36.000_UTC
   @param time  An input time formatted as a calendar time.
   <calenderTime> is as specified in 
   <tt>JSOC TN 07-001</tt> (http://jsoc.stanford.edu/doc/timerep.html)
   @param ord An input time formatted as an ordinal date.
   <ordinalDate> is yyyy.ddd[_zone], where @a zone is any supported 
   time system as specified in  <tt>JSOC TN 07-001</tt> 
   (http://jsoc.stanford.edu/doc/timerep.html).
   @param o The output format to be used. "jsoc" refers to internal time; 
   "sdo" refers to SDO time; "egse" refers to EGSE time; 
   "ord" refers to ordinal date; and "cal" refers to calendar time.

   @par Example to convert an internal time to the default output format (UTC calendar time):
   @code
   time_convert s=234253535.23
   @endcode

   @par Example to convert an internal time to the default output format (calendar time), but in TAI time:
   @code
   time_convert s=234253535.23 zone=TAI
   @endcode

   @par Example to convert an internal time to an EGSE time:
   @code
   time_convert s=234253535.23 o=egse
   @endcode

   @par Example to convert an EGSE time to a calendar time in the PDT time zone:
   @code
   time_convert egse=232533636.362 o=cal zone=PDT
   @endcode

   @par Example to convert an ordinal time to a calendar time in the TDT system:
   @code
   time_convert ord=2007.352 o=cal zone=TDT
   @endcode

   @par Example to convert a calendar time to an internal time:
   @code
   time_convert time=1998.02.04_06:00:17.230_UTC
   @endcode

   @par Example to convert a calendar time to an internal time (explicitly):
   @code
   time_convert time=1998.02.04_06:00:17.230_UTC o=jsoc
   @endcode

   @par Example to convert a calendar time to an SDO time
   @code
   time_convert time=1998.02.04_06:00:17.230_UTC o=sdo
   @endcode
*/
/* @{ */
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
#include "timeio.h"
#include "jsoc.h"
#include "cmdparams.h"

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "s", "NOT SPECIFIED", "<DSDS/JSOC time in seconds>"},
  {ARG_STRING, "sdo", "NOT SPECIFIED", "<SDO time in seconds>"},
  {ARG_STRING, "egse", "NOT SPECIFIED", "<EGSE time in seconds>"},
  {ARG_STRING, "time", "NOT SPECIFIED", "<Time as yyyy.mm...>"},
  {ARG_STRING, "ord", "NOT SPECIFIED", "<Time as yyyy.ddd...>"},
  {ARG_STRING, "zone", "UTC", "<Time zone>"},
  {ARG_STRING, "o", "NOT SPECIFIED", "format of time output"},
  {ARG_INT, "p", "0", "precision of seconds for time output"},
  {ARG_FLAG, "h", "0", "help message"},
  {ARG_END}
};

ModuleArgs_t *gModArgs = module_args;
/* @} */

CmdParams_t cmdparams;

#define SEC1970TO2004   1072828800
#define SECSPERDAY         86400.0

typedef enum
{
   kTimeFormat_None = 0,
   kTimeFormat_Internal,
   kTimeFormat_SDO,
   kTimeFormat_EGSE,
   kTimeFormat_Ordinal,
   kTimeFormat_Calendar
} TimeFormat_t;

char *TimeFormatStr[] = 
{
   "None",
   "Internal/JSOC",
   "SDO",
   "EGSE",
   "Ordinal",
   "Calendar"
};

static void PrintTime(TIME t, TimeFormat_t f, int precision);
static TIME ZoneAdjustment(TIME t, char *zone, int precision);

int nice_intro ()
  {
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\ntime_convert -h\n"
	"time_convert s=<secondsJSOC>|sdo=<secondsSDO>|egse=<secondsEGSE>|time=<calender time>|ord=<ordinal date> [o=jsoc | o=sdo | o=egse | o=ord | o=cal] [zone=<zone>]\n"
        "<secondsJSOC> = JSOC standard internal time, i.e. secs since 1997.01.01_TAI \n"
        "<secondsSDO> = seconds since January 1, 1958 TAI, i.e. SDO onboard time\n"
	"<secondsEGSE> = seconds since APPROXIMATELY January 1, 2004 \n"
        "<calender time> = yyyy.mm.dd_hh:mm:ss<zone>\n"
	"<ordinal date> = yyyy.ddd<zone>\n"
        "<zone> = time zone as UT, TAI, PST, etc. - default is UTC\n"
	"h - show usage message\n"
        "o - output time in format specified\n"
        "p - precision of seconds printed\n");
    return(1);
    }
  return (0);
  }

int main(int argc, char **argv)
{
int status;
TIME sscan_time(char *s);

char *s, *sdo, *time;
char *ord; /* Users can now supply the time argument as an ordinal date ('day-of-year' format - ISO 8601): YYYY.DDD */
char *egse = NULL;
char *oArg = NULL;
int precision = 0;

TIME t;

/* Parse command line parameters */
status = cmdparams_parse (&cmdparams, argc, argv);
if (status == CMDPARAMS_QUERYMODE)
  {
  cmdparams_usage (argv[0]);
  return 0;
  }

if (nice_intro ())
  return (0);

/* The are the different types of input. */
s = cmdparams_get_str(&cmdparams, "s", NULL);
sdo = cmdparams_get_str(&cmdparams, "sdo", NULL);
time = cmdparams_get_str(&cmdparams, "time", NULL);
ord = cmdparams_get_str(&cmdparams, "ord", NULL);
egse = cmdparams_get_str(&cmdparams, "egse", NULL);

/* Determine seconds precision */
precision = cmdparams_get_int(&cmdparams, "p", NULL);

int err = 0;

TIME timeToPrint;
TimeFormat_t inputFormat = kTimeFormat_None;
TimeFormat_t outputFormat = kTimeFormat_None;
TimeFormat_t outputOverride = kTimeFormat_None;
int mFormats = 0;

/* These are the different types of output. */
oArg = cmdparams_get_str(&cmdparams, "o", NULL);
if (strcmp(oArg, "NOT SPECIFIED") != 0)
{
   if (strcmp(oArg, "jsoc") == 0)
   {
      outputOverride = kTimeFormat_Internal;
   }
   else if (strcmp(oArg, "sdo") == 0)
   {
      outputOverride = kTimeFormat_SDO;
   }
   else if (strcmp(oArg, "egse") == 0)
   {
      outputOverride = kTimeFormat_EGSE;
   }
   else if (strcmp(oArg, "ord") == 0)
   {
      outputOverride = kTimeFormat_Ordinal;
   }
   else if (strcmp(oArg, "cal") == 0)
   {
      outputOverride = kTimeFormat_Calendar;
   }
   else
   {
      fprintf(stderr, "Output format %s is not recognized.  Using default.\n", oArg);
   }
}

if (strcmp(s, "NOT SPECIFIED")!= 0)
{ 
   /* internal seconds was specified */
   if (inputFormat == kTimeFormat_None)
   {
      sscanf(s, "%lf", &t);
      timeToPrint = t;
      outputFormat = kTimeFormat_Calendar;
      inputFormat =  kTimeFormat_Internal;
   }
   else
   {
      mFormats = 1;
   }
}
if (strcmp(sdo, "NOT SPECIFIED")!= 0)
{ 
   /* SDO seconds was specified */
   if (inputFormat == kTimeFormat_None)
   {
      sscanf(sdo, "%lf", &t);
      timeToPrint = t + sscan_time("1958.01.01_00:00:00_TAI");
      outputFormat = kTimeFormat_Calendar;
      inputFormat = kTimeFormat_SDO;
   }
   else
   {
      mFormats = 1;
   }
}
if (strcmp(egse, "NOT SPECIFIED") != 0)
{
   /* Time since ~ 1/1/2004 was provided as input */
   if (inputFormat == kTimeFormat_None)
   {
      sscanf(egse, "%lf", &t);
      timeToPrint = t + SEC1970TO2004 + sscan_time("1970.01.01_00:00_UTC");
      outputFormat = kTimeFormat_Calendar;
      inputFormat = kTimeFormat_EGSE;
   }
   else
   {
      mFormats = 1;
   }
}
if (strcmp(time, "NOT SPECIFIED")!= 0)
{ 
   /* calendar time string was specified */
   if (inputFormat == kTimeFormat_None)
   {
      t = sscan_time(time);
      timeToPrint = t;
      outputFormat = kTimeFormat_Internal;
      inputFormat = kTimeFormat_Calendar;
   }
   else
   {
      mFormats = 1;
   }
}
if (ord != NULL && *ord != NULL && strcmp(ord, "NOT SPECIFIED") != 0)
{
   /* day-of-year was specified - convert to internal representation */
   if (inputFormat == kTimeFormat_None)
   {
      err = 1;
      
      long year = 0;
      int day = 0;
      char *tFormat = NULL;
      
      char *ordDate = strdup(ord);
      int timeStrLen = strlen(ordDate);
      
      if (ordDate != NULL)
      {
	 char *loc = strchr(ordDate, '.');
	 char *loc2 = strchr(ordDate, '_');
	 
	 if (loc != NULL)
	 {
	    *loc = '\0';
	    
	    if (sscanf(ordDate, "%ld", &year) != 0)
	    {
	       if (loc2 != NULL && loc2 > loc)
	       {
		  *loc2 = '\0';
		  sscanf(loc + 1, "%d", &day);
		  
		  if (loc2 + 1 <= &ordDate[timeStrLen - 1])
		  {
		     tFormat = (char *)malloc(sizeof(char) * strlen(loc2 + 1) + 1);
		     if (tFormat)
		     {
			sscanf(loc2 + 1, "%s", tFormat);
		     }
		  }
	       }
	       else if (loc + 1 <= &ordDate[timeStrLen - 1])
	       {
		  sscanf((loc + 1), "%d", &day);
	       }
	       
	       err = (day < 1 || day > 366);
	       
	       if (!err)
	       {
		  /* first convert to calendar date */
		  char timeBuf[64];
		  snprintf(timeBuf, sizeof(timeBuf), "%ld.01.%d_00:00_%s", year, day, tFormat ? tFormat : "UT");
		  /* then use sscan_time to convert to internal time */
		  // t = sscan_time(timeBuf) - sscan_time("1977.01.01_00:00_TAI");
		  // sscan_time("1977.01.01_00:00_TAI") == 0
		  t = sscan_time(timeBuf);
		  timeToPrint = t;
		  outputFormat = kTimeFormat_Internal;
		  inputFormat = kTimeFormat_Ordinal;
	       }
	       else
	       {
		  fprintf(stderr, "invalid day of year\n");
	       }
	    }
	    else
	    {
	       fprintf(stderr, "invalid year\n");
	    }
	 }
	 else
	 {
	    fprintf(stderr, "missing day of year\n");
	 }
	 
	 free(ordDate);
      }
      
      if (tFormat)
      {
	 free(tFormat);
      }
   }
   else
   {
      mFormats = 1;
   }
}

if (inputFormat == kTimeFormat_None)
{
   fprintf(stderr, "No input format specified.\n");
   err = 1;
}
else if (mFormats)
{
   fprintf(stderr, 
	   "Multiple input formats specified, using %s time\n", 
	   TimeFormatStr[inputFormat]);
}

if (outputOverride != kTimeFormat_None)
{
   outputFormat = outputOverride;
}

if (!err)
{
   PrintTime(timeToPrint, outputFormat, precision);
}

if (err)
{
   fprintf(stderr,"time_convert call error\n");
   return(1);
}
else
{
   return(0);
}
}

/* t is time in TAI seconds since 1/1/1977, f is the format in which to print the time string */
void PrintTime(TIME t, TimeFormat_t f, int precision)
{
   if (f == kTimeFormat_Calendar)
   {
      char at[128];
      sprint_time(at, t, cmdparams_get_str(&cmdparams, "zone", NULL), precision);
      // broken for precision == 3, causes seg fault, fixed maybe
      // sprint_time(at, t, cmdparams_get_str(&cmdparams, "zone", NULL), 0);
      printf("%s\n", at);
   }
   else if (f == kTimeFormat_SDO)
   {
      /* no zone associated with sdo time, so zone parameter is ignored */
      printf("%12.*f\n", precision, t - sscan_time("1958.01.01_00:00:00_TAI"));
   }
   else if (f == kTimeFormat_EGSE)
   {
      /* no zone associated with egse time, so zone parameter is ignored */
      printf("%12.*f\n", precision, t - SEC1970TO2004 - sscan_time("1970.01.01_00:00_UTC"));
   }
   else if (f == kTimeFormat_Ordinal)
   {
      /* This does NOT provide a fractional day - any fractional part is discarded. */
      int ordTimeDays = 0;
      char zoneStr[32];
      char at[128];
      TIME jan1 = 0;
      TIME secsSinceJan1 = 0;

      char *zone =  cmdparams_get_str(&cmdparams, "zone", NULL);
      
      /* First, get internal time for January 1 of the year we are going to output. */
      sprint_time(at, t, zone, precision);
      char *tz = strrchr(at, '_');

      if (tz)
      {
	 snprintf(zoneStr, sizeof(zoneStr), "_%s", tz + 1);
      }

      char *dot = strchr(at, '.');
      if (dot != NULL)
      {
	 *dot = '\0';
	 char timeStr[64];
	 snprintf(timeStr, sizeof(timeStr), "%s.01.01_00:00:00%s", at, zoneStr);
	 jan1 = sscan_time(timeStr);
      }
      
      secsSinceJan1 = t - jan1;

      /* Now, get zone adjustment for jan1 and t.  Subtract Adj(t) - Adj(jan1) from 
       * secsSinceJan1. */
      TIME adjJan1 = ZoneAdjustment(jan1, zone, precision);
      TIME adjT = ZoneAdjustment(t, zone, precision);
      TIME ordTimeSecs = secsSinceJan1 - (adjT - adjJan1);

      /* Divide by TAI secs per day (86400). */
      ordTimeDays = 1 + ordTimeSecs / SECSPERDAY;      

      /* Output */
      printf("%s.%03d%s\n", at, ordTimeDays, zoneStr);
   }
   else
   {
      if (f != kTimeFormat_Internal)
      {
         // format error
         fprintf(stderr, "Invalid output format, defaulting to JSOC time.\n");
      }
         printf("%0.*f\n", precision, t);
   }
}

/* <t> is the TAI-seconds equivalent of a time in the <zone> time zone. */
TIME ZoneAdjustment(TIME t, char *zone, int precision)
{
   char timestr[128];
   sprint_time(timestr, t, zone, precision);

   char *tz = strrchr(timestr, '_');
   tz[1] = 'T';
   tz[2] = 'A';
   tz[3] = 'I';
   tz[4] = '\0';
   
   TIME taiDateSecs = sscan_time(timestr);
   return t - taiDateSecs;
}
