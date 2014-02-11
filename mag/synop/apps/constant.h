/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/constant.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: constant.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.2  90/12/18  15:51:56  phil
 * local.h included
 * 
 * Revision 1.1  90/09/17  16:40:20  phil
 * Initial revision
 * 
*/

#include <local.h>

#define PI	(3.14159265358979323844)
#define TWOPI	(6.28318530717958647688)
#define HALFPI  (PI*0.5)
#define DEGRAD	(57.29577951308232087721)
#define SECRAD  (3600.0*DEGRAD)
#define E	(2.71828182845904523536)
#define MISSING	(-8388608.0e10)
#define IMISSING (0x80000000) /* new value 9/6/90 phs */
#define RCMISS	(-8388608)
#define SMISS	(-32768)
#define BMISS	(-128)
#define TINY	(1.0e-30)
#define MAXLONG	(0x7FFFFFFF)
#define MAXSHORT (0x7FFF)
#define AU	(1.496e11)  /* 1 AU in meters  */
#define C	(2.9979e8)  /* c in m/s  */
#define DAYSECONDS (86400.0)
#define SID 	(86400.0)
#define DAYS_PER_CROT	(27.27527)
#define CDEG_PER_DAY	(360.0/27.27527)
#define RSUN	(6.95e8)	/* radius sun in meters */
