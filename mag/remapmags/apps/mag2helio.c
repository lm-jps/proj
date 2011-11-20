/* this JSOC module is a port of the SOI module described in the following comment. 
 * ported by Art Amezcua, revised by Tim Larson
 */

/*
 *  v2helio.c                            ~soi/(version)/src/modules/v2helio.c
 *
 *  This module interpolates CCD velocity data to estimate values that
 *    would be obtained at specified equal increments of heliographic 
 *    longitude and sine of latitude.  Apodization and corrections for
 *    solar rotation and limbshift are included. 
 *
 *  Responsible:  Kay Leibrand                   KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *    This module is under development.  Look for ??? in comments.
 *
 *  Planned updates:
 *    Decide when and where and if default values should be used.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "jsoc_main.h"
#include "astro.h"
#include "drms_dsdsapi.h"

#include "saveparm.c"
#include "obs2helio.c"
#include "obs2heliodb.c"
#include "fstats.h"

#define PI		(M_PI)
#define RADSINDEG 	(PI/180)
#define ARCSECSINRAD 	(3600*180/PI)
#define DAYSINYEAR	(365.2425)
#define SECSINDAY	(86400)
#define TAU_A           (499.004783806) // light travel time in seconds, = 1 AU/c
//#define TAU_A		(499.004782)	// this value used in old v2helio
#define TCARR		(25.38)		// days
#define RTRUE		(6.96000000e8)	// meters
#define AU		(149597870691)	// meters/au
//#define AU		(1.49597870e11)	// this value used in old v2helio
#define MAXLEN		(256)
#define kMAXROWS        65536

#define ARRLENGTH(ARR)  (sizeof(ARR)/sizeof(ARR[0]))
#define QUAL_NODATA	(0x80000000)
#define kNOTSPECIFIED   "not specified"

typedef enum
{
   V2HStatus_Success,
   V2HStatus_MissingParameter,
   V2HStatus_IllegalTimeFormat,
   V2HStatus_TimeConvFailed,
   V2HStatus_Unimplemented,
   V2HStatus_IllegalOrientation
} V2HStatus_t;

char *module_name = "mag2helio";

ModuleArgs_t module_args[] = 
{
   {ARG_STRING,  "in", "", "input data records"},
   {ARG_STRING,  "out", "", "output data series"},
   {ARG_STRING,  "segin", kNOTSPECIFIED, "name of input segment if not using segment 0"},
   {ARG_STRING,  "segout", kNOTSPECIFIED, "name of output segment if not using segment 0"},
//   {ARG_STRING,  "histlink", "HISTORY", "name of link to ancillary dataseries for processing metadata"},
//   {ARG_STRING,  "srclink",  "FDMr12m", "name of link to source data"},
   {ARG_INT,     "PERM", "1", "set to 0 for transient records, nonzero for permanent records"},
   {ARG_INT,     "VERB", "1", "option for level of verbosity: 0=only error and warning messages; >0=print messages outside of loop; >1=print messages inside loop; >2=debugging output", ""},  /* debug messages */
   {ARG_FLOAT,   "MAPRMAX",  "0.994", "maximum image radius", ""},
   {ARG_INT,     "RMAXFLAG", "0", "set to nonzero to set output to DRMS_MISSING outside MAPRMAX, otherwise uses 0.0 outside MAPRMAX"},  /* Non 0 sets data outside RMAX MISSING */   
   {ARG_INT,     "MAPMMAX",  "1800", "determines mapcols", ""},	/* determines mapcols, default value is 3*512 */
   {ARG_FLOAT,   "MAPLGMAX", "90.0", "longitude maximum, degrees", ""},	/* degrees */    
   {ARG_FLOAT,   "MAPLGMIN", "-90.0", "longitude minimum, degrees", ""},
   {ARG_FLOAT,   "MAPBMAX",  "90.0", "latitude maximum, degrees, also used for minimum", ""},  
   {ARG_INT,     "SINBDIVS", "720", "number of increments in sin latitude from 0 to 1", ""},	/* # of = increments in sinB from sin(0) to sin(PI/2) */
   {ARG_INT,     "LGSHIFT",  "0", "option for longitude shift: 0=none; 1=fixed rate; 2=nearest degree; 3=nearest tenth of a degree", ""}, /* 0=none; 1=fixed rate; 2=nearest Degree */
   {ARG_TIME,    "REF_T0",   "1987.01.03_17:31:12_TAI", "reference time for computing fixed rate longitude shift", ""},
   {ARG_FLOAT,   "REF_L0",   "0.0", "reference longitude for computing fixed rate longitude shift ", ""},
   {ARG_FLOAT,   "SOLAR_P",  "999.0", "P-angle; if unset, taken from keywords", ""},	/* can't use D_MISSING here */
   {ARG_FLOAT,   "PSIGN",    "1.0", "sign of SOLAR_P", ""},	/* Sign of P. For MWO data. */
   {ARG_FLOAT,   "PERR",     "0.0", "fixed P-angle error, likely -0.22", ""},	/* Fixed P-angle error. Maybe -0.22. */
   {ARG_FLOAT,   "IERR",     "0.0", "error in Carrington inclination, likely -0.10", ""},	/* Error in Carrington inclination. Maybe -0.10. */
   {ARG_TIME,    "REF_TB0",  "2001.06.06_06:57:22_TAI", "reference time for computing correction to P and B angles, roughly when B0=0", ""},
   {ARG_INT,     "INTERPO",  "1", "option for interpolation: 0=bilinear; 1=cubic convolution", ""},	/* 2 methods - see soi_fun.h */
   {ARG_INT,     "APODIZE",  "0", "option for apodization: 0=none; 1=use true solar coordinates; 2=use ideal solar coordinates (b0=0)", ""},	/* see soi_fun.h or apodize.c */
   {ARG_FLOAT,   "APINNER",  "0.90", "start of apodization in fractional image radius", ""},	/* start of apodization */
   {ARG_FLOAT,   "APWIDTH",  "0.05", "width of apodization in fractional image radius", ""},	/* width of apodization */
   {ARG_INT,     "APEL",     "0", "set to nonzero for elliptical apodization described by APX and APY", ""},	/* do elliptical apodization described by apx and apy */
   {ARG_FLOAT,   "APX",      "1.00", "divide the x position by this before applying apodization", ""},	/* divide the x position by this before applying apodization */
   {ARG_FLOAT,   "APY",      "1.00", "divide the y position by this before applying apodization", ""},	/* divide the y position by this before applying apodization */
   {ARG_INT,     "VCORLEV",  "0", "option for velocity correction: 0=none; 1=subtract a model of differential rotation; 2=also divide by line of sight projection factor for purely radial velocities", ""}, 	/* 3 levels - see soi_fun.h*/
   {ARG_INT,     "MCORLEV",  "0", "option for magnetic correction: 0=none; 1=line of sight", ""}, 	/* 2 levels - see soi_fun.h*/
   {ARG_INT,     "MOFFSETFLAG",  "0", "set to nonzero to get BFITZERO from input record and subtract from data before interpolating", ""}, 	/* 1=apply BFITZERO correction*/
   {ARG_FLOAT,   "OUTSCALE", "1.0", "bscale to use for output", ""},    /* scale for output */
   {ARG_FLOAT,   "OUTBIAS",  "0.0", "bzero to use for output", ""},    /* bias for scaled output */
   {ARG_INT,     "DISTORT",  "0", "option for distortion correction: 0=none; 1=full disk(fd) data; 2=vector-weighted(vw) data", ""}, /* 0 for none, 1 for FD, 2 for vw */
   {ARG_FLOAT,   "CUBIC",    "7.06E-9", "cubic distortion in fd units", ""}, /* Cubic distortion in FD units */
   {ARG_FLOAT,   "TILTALPHA","2.59", "tilt of CCD, degrees", ""}, /* TILT of CCD in degrees */
   {ARG_FLOAT,   "TILTBETA", "56.0", "direction of CCD tilt, degrees", ""}, /* Direction of TILT in degrees */
   {ARG_FLOAT,   "TILTFEFF", "12972.629", "effective focal length", ""}, /* Effective focal length */
   {ARG_INT,     "OFLAG",  "0", "set to nonzero to force reading orientation from keyword, otherwise \"SESW\" is assumed)", ""},  /* Non 0 skips checko (SESW assumed) */
   {ARG_INT,     "DATASIGN", "0", "value to multiply data; set to 0 to take DATASIGN from keyword, or 1.0 if not found", ""}, 	/* Non 0 forces datasign to value*/
   {ARG_INT,     "MAXMISSVALS", "0", "if >0, this becomes threshold on MISSVALS from keyword", ""},  /* max. allowed MISSING pixels */
   {ARG_INT,     "CARRSTRETCH", "0", "set to nonzero to correct for differential rotation according to DIFROT_[ABC]", ""},  /* 0 - don't correct for diff rot, 1 - correct */
   {ARG_FLOAT,   "DIFROT_A",   "13.562", "A coefficient in differential rotation adjustment (offset)", ""}, /* A coefficient in diff rot adj (offset) */
   {ARG_FLOAT,   "DIFROT_B",   "-2.04", "B coefficient (to sin(lat) ^ 2)", ""},  /* B coefficient (to sin(lat) ^ 2) */
   {ARG_FLOAT,   "DIFROT_C",   "-1.4875", "C coefficient (to sin(lat) ^ 4)", ""},  /* C coefficient (to sin(lat) ^ 4) */
   {ARG_END}
};


double getwalltime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000.0 + tv.tv_usec/1000.0;
}

double getcputime(double *utime, double *stime)
{

  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  *utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec/1000.0;
  *stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec/1000.0;
  return *utime + *stime;
}


static void CheckO(const char *orientation, V2HStatus_t *stat);

static inline void ParameterDef(int status, 
				char *pname, 
				double defaultp, 
				double *p, 
				int iRec,
                                int verbflag)
{
   /* logs warning and sets parameter to default value */
   if (status != 0)
   {
      *p = defaultp;
      if (verbflag)
         fprintf(stderr, "WARNING: default value %g used for %s, iRec = %d, status = %d\n", defaultp, pname, iRec, status);
   }
}


#define PARAMETER_ERROR(PNAME)                                     \
  if (status != DRMS_SUCCESS)                                      \
  {                                                                \
    CreateBlankRecord(inrec, outrec, quality);                     \
    fprintf(stderr,                                                \
	    "SKIP: error getting keyword %s: iRec = %d, status = %d, T_REC = %s, recnum = %lld\n",  \
	    PNAME,                                                 \
            iRec,                                                  \
  	    status,                                                \
            trecstr,                                               \
            inrec->recnum);                                        \
    if (inarr)  drms_free_array(inarr);                            \
    if (outarr) drms_free_array(outarr);                           \
    if (orientation) free(orientation);                            \
    continue;                                                      \
  }


/* Segment  will be empty, but there will be a record! */
static void CreateBlankRecord(DRMS_Record_t *inrec, DRMS_Record_t *outrec, int quality)
{
  /* create 'blank' data */
// might insert 'quality = quality | MASK' here.
   quality = quality | QUAL_NODATA;
   drms_copykey(outrec, inrec, "T_REC");
   drms_setkey_int(outrec, "QUALITY", quality);
   drms_close_record(outrec, DRMS_INSERT_RECORD);
}


int DoIt(void)
{

   int newstat = 0;
   int status = DRMS_SUCCESS;
   V2HStatus_t vstat = V2HStatus_Success;

   const char *orientationdef = "SESW    ";
   char *orientation = NULL;
   int paramsign;
   int longitude_shift, velocity_correction, interpolation, apodization;
   int mag_correction;
   int mag_offset;
   int row;
   int mapped_lmax, sinb_divisions, mapcols, maprows, nmax, nmin;
   int length[2];
   int carrStretch = 0;
   float diffrotA = 0.0;
   float diffrotB = 0.0;
   float diffrotC = 0.0;
   double tobs, tmearth, tref, trefb0;
   double smajor, sminor, sangle;
   double xscale, yscale, imagescale;
   int xpixels, ypixels, pixels;
   double obs_vr, obs_vw, obs_vn;
   double b0, bmax, bmin, sinBdelta;
   double longmax, longmin, longmax_adjusted, longmin_adjusted, longinterval;
   double p0, p, rmax;
   double ierr, perr, psign;
   double x0, y0;
   double obsdist, longshift, obsl0, refl0, mapl0, longrate, rtrue, rsun, S;
   double rsunDef;
   int obsCR;
   int apel;
   double apinner, apwidth, apx, apy;
   double scale, bias;
   double colsperdeg;

   int quality;
   double satrot, instrot;
   double dsignout, vsign;
   int distsave;
   double cubsave, tiltasave, tiltbsave, tiltfsave;
   TIME trec, tnow, UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
   char trecstr[100];

  double statMin, statMax, statMedn, statMean, statSig, statSkew, statKurt;
  int statNgood;

   LIBASTRO_RotRate_t rRates[kMAXROWS];
   LIBASTRO_Dist_t distP;
   DRMS_Segment_t *segin = NULL;
   DRMS_Segment_t *segout = NULL;
   DRMS_Array_t *inarr = NULL;
   DRMS_Array_t *outarr = NULL;
   DRMS_Record_t *inrec = NULL;
   DRMS_Record_t *outrec = NULL;
   DRMS_Type_t usetype = DRMS_TYPE_FLOAT;
   DRMS_RecLifetime_t lifetime;

   long long histrecnum=-1;

   int errbufstat=setvbuf(stderr, NULL, _IONBF, BUFSIZ);
   int outbufstat=setvbuf(stdout, NULL, _IONBF, BUFSIZ);

   double wt0, wt1, wt;
   double ut0, ut1, ut;
   double st0, st1, st;
   double ct0, ct1, ct;

   wt0=getwalltime();
   ct0=getcputime(&ut0, &st0);

   char *inRecQuery = cmdparams_save_str(&cmdparams, "in", &newstat);
   char *outSeries = cmdparams_save_str(&cmdparams, "out", &newstat);

   int checko = cmdparams_save_int(&cmdparams, "OFLAG", &newstat);
   int verbflag = cmdparams_save_int(&cmdparams, "VERB", &newstat);
   int NaN_beyond_rmax = cmdparams_save_int(&cmdparams, "RMAXFLAG", &newstat);
   int maxmissvals = cmdparams_save_int(&cmdparams, "MAXMISSVALS", &newstat);
   int permflag = cmdparams_save_int(&cmdparams, "PERM", &newstat);
   if (permflag)
      lifetime = DRMS_PERMANENT;
   else
      lifetime = DRMS_TRANSIENT;

   // read Carrington coordinate correction for differential sun rotation 
   carrStretch = cmdparams_save_int(&cmdparams, "CARRSTRETCH", &newstat);
   diffrotA = cmdparams_save_float(&cmdparams, "DIFROT_A", &newstat);
   diffrotB = cmdparams_save_float(&cmdparams, "DIFROT_B", &newstat);
   diffrotC = cmdparams_save_float(&cmdparams, "DIFROT_C", &newstat);

   longrate = 360.0 / TCARR - 360.0 / DAYSINYEAR; // degrees per day 
   longrate /= SECSINDAY; // degrees per sec 
   rtrue = RTRUE/AU;  // au 

   apodization = cmdparams_save_int(&cmdparams, "APODIZE", &newstat);
   apinner = cmdparams_save_double(&cmdparams, "APINNER", &newstat);  
   apwidth = cmdparams_save_double(&cmdparams, "APWIDTH", &newstat);
   apel = cmdparams_save_int(&cmdparams, "APEL", &newstat);
   apx = cmdparams_save_double(&cmdparams, "APX", &newstat);
   apy = cmdparams_save_double(&cmdparams, "APY", &newstat);
   longitude_shift = cmdparams_save_int(&cmdparams, "LGSHIFT", &newstat);
   mag_correction = cmdparams_save_int(&cmdparams, "MCORLEV", &newstat);
   mag_offset = cmdparams_save_int(&cmdparams, "MOFFSETFLAG", &newstat);
   velocity_correction = cmdparams_save_int(&cmdparams, "VCORLEV", &newstat);
   interpolation = cmdparams_save_int(&cmdparams, "INTERPO", &newstat);
   paramsign = cmdparams_save_int(&cmdparams, "DATASIGN", &newstat);
   rmax = cmdparams_save_double(&cmdparams, "MAPRMAX", &newstat);
   refl0 = cmdparams_save_double(&cmdparams, "REF_L0", &newstat);

   distsave = cmdparams_save_int(&cmdparams, "DISTORT", &newstat);
   cubsave = cmdparams_save_double(&cmdparams, "CUBIC", &newstat);
   tiltasave = cmdparams_save_double(&cmdparams, "TILTALPHA", &newstat);
   tiltbsave = cmdparams_save_double(&cmdparams, "TILTBETA", &newstat);
   tiltfsave = cmdparams_save_double(&cmdparams, "TILTFEFF", &newstat);

   scale = cmdparams_save_double(&cmdparams, "OUTSCALE", &newstat);
   bias = cmdparams_save_double(&cmdparams, "OUTBIAS", &newstat);
   p0 = cmdparams_save_double(&cmdparams, "SOLAR_P", &newstat);
   psign = cmdparams_save_double(&cmdparams, "PSIGN", &newstat);
   perr = cmdparams_save_double(&cmdparams, "PERR", &newstat);
   ierr = cmdparams_save_double(&cmdparams, "IERR", &newstat);
   trefb0 = cmdparams_save_time(&cmdparams, "REF_TB0", &newstat);

   SetDistort(distsave, cubsave, tiltasave, tiltbsave, tiltfsave, &distP);

   tref = cmdparams_save_time(&cmdparams, "REF_T0", &newstat);

   // determine mapcols and adjust longmin and longmax */
   mapped_lmax = cmdparams_save_int(&cmdparams, "MAPMMAX", &newstat); 
   longmax = cmdparams_save_double(&cmdparams, "MAPLGMAX", &newstat); /* degrees */
   longmin = cmdparams_save_double(&cmdparams, "MAPLGMIN", &newstat); /* degrees */
   longinterval = (180.0) / mapped_lmax;	                     /* degrees */
 
   // This does not always handle the case where 1/longinterval is an integer correctly.

   // the next two statement do nothing, right? 
   // why do nmin and max keep getting set with different RHSs? 
   nmin = (int)(longmin / longinterval); // round towards 0 
   nmax = (int)(longmax / longinterval); // round towards 0 
   colsperdeg = mapped_lmax / 180.0;
   nmin = (int)(longmin * colsperdeg); // round towards 0 
   nmax = (int)(longmax * colsperdeg); // round towards 0 
   mapcols = nmax - nmin + 1;
   longmin_adjusted = nmin * longinterval;
   longmax_adjusted = nmax * longinterval;

   // determine maprows, bmax, bmin, and sinBdelta
   sinb_divisions = cmdparams_save_int(&cmdparams, "SINBDIVS", &newstat);
   sinBdelta = 1.0/sinb_divisions;
   bmax = cmdparams_save_double(&cmdparams, "MAPBMAX", &newstat);     // degrees
   bmin = -bmax; 
   nmax = (int)(sin(RADSINDEG*bmax)*sinb_divisions); // round towards 0 
   maprows = 2*nmax;

   length[0] = mapcols;
   length[1] = maprows;

   char *segnamein = cmdparams_save_str(&cmdparams, "segin", &newstat);
   char *segnameout = cmdparams_save_str(&cmdparams, "segout", &newstat);
   int seginflag = strcmp(kNOTSPECIFIED, segnamein);
   int segoutflag = strcmp(kNOTSPECIFIED, segnameout);

//   char *histlinkname = cmdparams_save_str(&cmdparams, "histlink", &newstat);
//   char *srclinkname = cmdparams_save_str(&cmdparams, "srclink", &newstat);

   // set up ancillary dataseries for processing metadata
   DRMS_Record_t *tempoutrec = drms_create_record(drms_env, 
					       outSeries,
					       DRMS_TRANSIENT, 
					       &status);

   if (status != DRMS_SUCCESS) 
   {
     fprintf(stderr,"ERROR: couldn't open a record in output dataseries %s, status = %d\n", outSeries, status);
     return 1;
   }

//   DRMS_Link_t *histlink = hcon_lookup_lower(&tempoutrec->links, histlinkname);
//   DRMS_Link_t *srclink = hcon_lookup_lower(&tempoutrec->links, srclinkname);

/*
   if (histlink != NULL) 
   {
      DRMS_Record_t *histrec = drms_create_record(drms_env, 
						  histlink->info->target_series,
						  DRMS_PERMANENT, 
						  &status);
      if (status != DRMS_SUCCESS)
      {
         fprintf(stderr,"ERROR: could not open a record in history dataseries %s, status = %d\n", histlink->info->target_series, status);
         return 1;
      }

      histrecnum = histrec->recnum;
      tnow = (double)time(NULL);
      tnow += UNIX_epoch;
      status = drms_setkey_time(histrec, "DATE", tnow);
      if (status != DRMS_SUCCESS)
      {
         fprintf(stderr,"ERROR: problem writing keyword DATE in history dataseries, status = %d\n", status);
         drms_close_record(histrec, DRMS_FREE_RECORD);
         return 1;
      }
      status = drms_setkey_string(histrec, "MODNAME", cmdparams.argv[0]);
      if (status != DRMS_SUCCESS)
      {
         fprintf(stderr,"ERROR: problem writing keyword MODNAME in history dataseries, status = %d\n", status);
         drms_close_record(histrec, DRMS_FREE_RECORD);
         return 1;
      }
      status = drms_setkey_string(histrec, "ARGSUSED", savestr);
      if (status != DRMS_SUCCESS)
      {
         fprintf(stderr,"ERROR: problem writing keyword DATE in history dataseries, status = %d\n", status);
         drms_close_record(histrec, DRMS_FREE_RECORD);
         return 1;
      }
      drms_close_record(histrec, DRMS_INSERT_RECORD);

   }
   else
   {
      fprintf(stderr,"WARNING: could not find history link in output dataseries\n");
   }

*/

   if (newstat) 
   {
     fprintf(stderr, "ERROR: problem with input arguments, status = %d, diagnosis follows\n", newstat);
     cpsave_decode_error(newstat);
     return 0;
   }  
   else if (savestrlen != strlen(savestr)) 
   {
     fprintf(stderr, "ERROR: problem with savestr, savestrlen = %d, strlen(savestr) = %d\n", savestrlen, (int)strlen(savestr));
     return 0;
   }


// these must be present in the output dataseries and variable, not links or constants   
   char *outchecklist[] = {"T_REC", "QUALITY", "MAPMMAX", "SINBDIVS", "DATASIGN",
                           "CRPIX1", "CRVAL1", "CDELT1",
                           "CRPIX2", "CRVAL2", "CROTA2", "CDELT2" };

   int itest;
   for (itest=0; itest < ARRLENGTH(outchecklist); itest++)
   {
      DRMS_Keyword_t *outkeytest = hcon_lookup_lower(&tempoutrec->keywords, outchecklist[itest]);
      if (!outkeytest || outkeytest->info->islink || outkeytest->info->recscope == 1)
      {
         fprintf(stderr, "ERROR: output keyword %s is either missing, constant, or a link\n", outchecklist[itest]);
         return 0;
      }
   }
   drms_close_record(tempoutrec, DRMS_FREE_RECORD);


   DRMS_RecordSet_t *inRecSet = drms_open_records(drms_env, inRecQuery, &status);
   int nRecs = inRecSet->n;
   if (status != DRMS_SUCCESS || inRecSet == NULL)
   {
      fprintf(stderr, "ERROR: problem opening input recordset: status = %d\n", status);
      return 0;
   }

   if (nRecs == 0)
   {
      printf("WARNING: input recordset contains no records\nmodule %s successful completion\n", cmdparams.argv[0]);
      return 0;
   }

   if (verbflag) 
      printf("input recordset opened, nRecs = %d\n",nRecs);


// go ahead and check for the presence of these input keywords as well
   char *inchecklist[] = {"T_REC", "QUALITY", "T_OBS", "CRLT_OBS", "CRLN_OBS", 
                          //"SAT_ROT", "INST_ROT", "IM_SCALE",
                          "CDELT1", "CDELT2"};

   DRMS_Record_t *tempinrec = inRecSet->records[0];
   for (itest=0; itest < ARRLENGTH(inchecklist); itest++)
   {
      DRMS_Keyword_t *inkeytest = hcon_lookup_lower(&tempinrec->keywords, inchecklist[itest]);
      if (!inkeytest)
      {
         fprintf(stderr, "ERROR: required input keyword %s is missing\n", inchecklist[itest]);
         return 0;
      }
   }


   int iRec;
   int error=0; // only set error before a break
   int nsuccess=0;
   for (iRec = 0; iRec < nRecs; iRec++)
   {
      if (verbflag > 1) 
      {
         wt1=getwalltime();
         ct1=getcputime(&ut1, &st1);
         printf("processing record %d\n", iRec);
      }
      inrec = inRecSet->records[iRec];

      // create an output record
      outrec = drms_create_record(drms_env, 
                                  outSeries,
                                  lifetime, 
                                  &status);

      if (status != DRMS_SUCCESS || outrec==NULL)
      {
         fprintf(stderr,"ERROR: unable to open record in output dataseries, status = %d\n", status);
         error=2;
         break;
      }

//      if (histlink)
//         drms_setlink_static(outrec, histlinkname,  histrecnum);

//      if (srclink)
//         drms_setlink_static(outrec, srclinkname,  inrec->recnum);


      if (segoutflag)
         segout = drms_segment_lookup(outrec, segnameout);
      else
         segout = drms_segment_lookupnum(outrec, 0);

      // create an output array
      outarr = drms_array_create(DRMS_TYPE_FLOAT, 2, length, NULL, &status);

      if (!segout || !outarr || status != DRMS_SUCCESS)
      {
         fprintf(stderr, "ERROR: problem with output segment or array: iRec = %d, status = %d\n", iRec, status);
         if (outarr)
            drms_free_array(outarr);
         error=1;
         break;
      }

      trec = drms_getkey_time(inrec, "T_REC", &status);
      if (status != DRMS_SUCCESS)
      {
         fprintf(stderr,
	    "SKIP: error getting prime keyword %s: iRec = %d, status = %d, recnum = %lld\n",
	    "T_REC",
            iRec,
  	    status,
            inrec->recnum);
         drms_free_array(outarr);
         continue;
      }
      sprint_time(trecstr, trec, "TAI", 0); 

      quality = drms_getkey_int(inrec, "QUALITY", &status);
      if (status != DRMS_SUCCESS)
      {
        CreateBlankRecord(inrec, outrec, quality);
        fprintf(stderr,
	    "SKIP: error getting keyword %s: iRec = %d, status = %d, T_REC = %s, recnum = %lld\n",
	    "QUALITY",
            iRec,
  	    status,
            trecstr,
            inrec->recnum);
        drms_free_array(outarr);
        continue;
      }


      // insert tests on quality here.
      // if we encounter a keyword error causing a continue within the loop, we should set a bit in quality for this.
      // in  other words, CreateBlankRecord should always be preceded by setting quality, or could move this inside CreateBlankRecord.
      if (quality & QUAL_NODATA)
      {
         CreateBlankRecord(inrec, outrec, quality);
         fprintf(stderr,"SKIP: record rejected based on quality = %d: iRec = %d, T_REC = %s, recnum = %lld\n", quality, iRec, trecstr, inrec->recnum);
         drms_free_array(outarr);
         continue;
      }

      if (seginflag)
         segin = drms_segment_lookup(inrec, segnamein);
      else
         segin = drms_segment_lookupnum(inrec, 0);

      if (segin)
         inarr = drms_segment_read(segin, DRMS_TYPE_FLOAT, &status);

      if (!segin || !inarr || status != DRMS_SUCCESS)
      {
//         CreateBlankRecord(inrec, outrec, quality);
         fprintf(stderr, "ERROR: problem with input segment or array: iRec = %d, status = %d, T_REC = %s, recnum = %lld \n", iRec, status, trecstr, inrec->recnum);
         drms_free_array(outarr);
         if (inarr)
            drms_free_array(inarr);
//         error=1;
         continue;
//         break;
      }

      if (checko)
      {
         orientation = drms_getkey_string(inrec, "ORIENT", &status);
         if (status == DRMS_ERROR_UNKNOWNKEYWORD)
         {
            fprintf(stderr, "ERROR: required keyword %s missing, iRec = %d, T_REC = %s, recnum = %lld\n", "ORIENT", iRec, trecstr, inrec->recnum);
            drms_free_array(inarr);
            drms_free_array(outarr);
            error=1;
            break;
         }
         PARAMETER_ERROR("ORIENT")
         CheckO(orientation, &vstat);
         if (vstat != V2HStatus_Success)
         { 
            CreateBlankRecord(inrec, outrec, quality);
            fprintf(stderr,"SKIP: illegal %s: iRec = %d, T_REC = %s, recnum = %lld\n", "ORIENT", iRec, trecstr, inrec->recnum);
            drms_free_array(inarr);
            drms_free_array(outarr);
            free(orientation);
            continue;
         }
      }
      else
      {
         orientation = strdup(orientationdef);
      }

      if (maxmissvals > 0) 
      {
         int missvals = drms_getkey_int(inrec, "MISSVALS", &status);
         if (status == DRMS_ERROR_UNKNOWNKEYWORD)
         {
            fprintf(stderr, "ERROR: required keyword %s missing, iRec = %d, T_REC = %s, recnum = %lld\n", "MISSVALS", iRec, trecstr, inrec->recnum);
            drms_free_array(inarr);
            drms_free_array(outarr);
            free(orientation);
            error=1;
            break;
         }
         PARAMETER_ERROR("MISSVALS")
         if (missvals > maxmissvals) 
         {
            CreateBlankRecord(inrec, outrec, quality);
            fprintf(stderr, "SKIP: %d pixels MISSING, max allowed is %d: iRec = %d, T_REC = %s, recnum = %lld\n", missvals, maxmissvals, iRec, trecstr, inrec->recnum);
            drms_free_array(inarr);
            drms_free_array(outarr);
            free(orientation);
            continue;
         }
      }

      // assemble arguments

      tobs = drms_getkey_time(inrec, "T_OBS", &status);
      PARAMETER_ERROR("T_OBS")

      // MDI keyword was OBS_B0
      b0 = drms_getkey_double(inrec, "CRLT_OBS", &status);
      PARAMETER_ERROR("CRLT_OBS")

      // MDI keyword was OBS_L0
      obsl0 = drms_getkey_double(inrec, "CRLN_OBS", &status);
      PARAMETER_ERROR("CRLN_OBS")

      // MDI keyword was OBS_CR
      obsCR = drms_getkey_int(inrec, "CAR_ROT", &status);
   
      if (status != DRMS_SUCCESS)
      {
         // May 2, 2005. Commented this out to allow for no OBS_CR in MWO data. Still prints warning. OBS_CR is not used by this program or by helio2mlat.
         fprintf(stderr, "WARNING: missing %s, iRec = %d, status = %d, T_REC = %s, recnum = %lld\n", "CAR_ROT", iRec, status, trecstr, inrec->recnum);
      }

      if (p0 == 999.0) 
      {
         // MDI keyword was SOLAR_P = -(SAT_ROT + INST_ROT)
/*
         satrot = drms_getkey_double(inrec, "SAT_ROT", &status);
         PARAMETER_ERROR("SAT_ROT")
         instrot = drms_getkey_double(inrec, "INST_ROT", &status);
         PARAMETER_ERROR("INST_ROT")
         p=-(satrot+instrot);
*/
         double crota = drms_getkey_double(inrec, "CROTA2", &status);
         PARAMETER_ERROR("CROTA2")
         p=-crota;
      } 
      else 
      {
         p = p0;
      }

      // fix for 1988 MWO
      // p = 180.0 - p ;

      p = psign * p ;

      // 991839442. corresponds to hour 73878 minute 57 second 22
      // or 73878.956 or day 3078.2898, roughly when B0 is 0
      // b0=b0 * 0.986207; One way of correcting. 
      // The following is pretty good 
      // b0=b0-0.1*sin((tobs-991839442.)/31557600.*2*PI);
      // b0 = b0 + ierr * sin((tobs - 991839442.) / 31557600. * 2 * PI);
      b0 = b0 + ierr * sin((tobs - trefb0) / 31557600. * 2 * PI);

      // p=p-0.2; 
      // p=p-0.2+0.1*cos((tobs-991839442.)/31557600.*2*PI); 
      // p = p + perr - ierr * cos((tobs - 991839442.) / 31557600. * 2 * PI);
      p = p + perr - ierr * cos((tobs - trefb0) / 31557600. * 2 * PI);

      // S_MAJOR, S_MINOR, S_ANGLE, X_SCALE, Y_SCALE were MDI keywords, not used for HMI, but still necessary for GONG data
      smajor = drms_getkey_double(inrec, "S_MAJOR", &status);
      ParameterDef(status, "S_MAJOR", 1.0, &smajor, iRec, 0);

      sminor = drms_getkey_double(inrec, "S_MINOR", &status);
      ParameterDef(status, "S_MINOR", 1.0, &sminor, iRec, 0);

      sangle = drms_getkey_double(inrec, "S_ANGLE", &status);
      ParameterDef(status, "S_ANGLE", 0.0, &sangle, iRec, 0);

      xscale = drms_getkey_double(inrec, "CDELT1", &status);
      PARAMETER_ERROR("CDELT1")
      yscale = drms_getkey_double(inrec, "CDELT2", &status);
      PARAMETER_ERROR("CDELT2")

      // use xscale and yscale for the following check, then set to 1.0 for the call to obs2helio
      if (xscale != yscale)
      {
         CreateBlankRecord(inrec, outrec, quality);
         fprintf(stderr, "SKIP: %s != %s not supported, iRec = %d, T_REC = %s, recnum = %lld \n", "CDELT1", "CDELT2", iRec, trecstr, inrec->recnum);
         drms_free_array(inarr);
         drms_free_array(outarr);
         free(orientation);
         continue;
      }
imagescale=xscale;
      xscale=1.0;
      yscale=1.0;

/*
      imagescale = drms_getkey_double(inrec, "IM_SCALE", &status);
      PARAMETER_ERROR("IM_SCALE")
*/

      if (paramsign != 0)
      {
         vsign = paramsign;
      }
      else
      {
         vsign = drms_getkey_double(inrec, "DATASIGN", &status);
         ParameterDef(status, "DATASIGN", 1.0, &vsign, iRec, 1);
      }

      if (velocity_correction) 
      {
         obs_vr = drms_getkey_double(inrec, "OBS_VR", &status);
         ParameterDef(status, "OBS_VR", 0.0, &obs_vr, iRec, 1);

         obs_vw = drms_getkey_double(inrec, "OBS_VW", &status);
         ParameterDef(status, "OBS_VW", 0.0, &obs_vw, iRec, 1);

         obs_vn = drms_getkey_double(inrec, "OBS_VN", &status);
         ParameterDef(status, "OBS_VN", 0.0, &obs_vn, iRec, 1);
      }

      // MDI keyword was OBS_DIST, in AU
      obsdist = drms_getkey_double(inrec, "DSUN_OBS", &status);
      obsdist /= AU;  // note that an incorrect value of 1.49597892e11 has sometimes been used to convert between OBS_DIST and DSUN_OBS when porting data from DSDS to DRMS 
      ParameterDef(status, "DSUN_OB", 1.0, &obsdist, iRec, 1);
      S = rtrue / obsdist; // radians - approx. arcsin(rtrue/obsdist)

/*
      rsun = drms_getkey_double(inrec, "R_SUN", &status);
      rsunDef = ARCSECSINRAD * S / sqrt(1.0 - S * S) / imagescale;	   
      ParameterDef(status, "R_SUN", rsunDef, &rsun, iRec, 1);
*/
      rsun = drms_getkey_double(inrec, "RSUN_OBS", &status);
      rsun /= imagescale;


      if (longitude_shift == 1) 
      {
         tmearth = tobs+TAU_A*(1.0-obsdist); 
         longshift = (obsl0-refl0)+longrate*(tmearth-tref); // degrees 
         while (longshift > 180.0) longshift-=360.0; 
         while (longshift < -180.0) longshift+=360.0;
      }
      else if (longitude_shift == 2) // Shift center to nearest Carrington Degree 
      {
         longshift =  obsl0 - (int)(obsl0);
         if (longshift > 0.5) longshift -= 1.0;
      }
      else if (longitude_shift == 3)  // Shift center to nearest tenth of a degree 
      {
         longshift = (obsl0 * 10 - (int)(obsl0 * 10)) / 10;
         if (longshift > 0.5) longshift -= 1.0;
      }
      else
      {
         longshift = 0.0;
      }

      mapl0 = obsl0 - longshift;

      xpixels = inarr->axis[0];
      ypixels = inarr->axis[1];
      pixels  = xpixels * ypixels;

      // MDI keyword was X0
      x0 = drms_getkey_double(inrec, "CRPIX1", &status);
      ParameterDef(status, "CRPIX1", xpixels / 2, &x0, iRec, 1);
      x0 -= 1.0;

      // MDI keyword was Y0
      y0 = drms_getkey_double(inrec, "CRPIX2", &status);
      ParameterDef(status, "CRPIX2", ypixels / 2, &y0, iRec, 1);
      y0 -= 1.0;

      if (mag_offset) 
      {
         float *dat = (float *)inarr->data;
         double bfitzero = drms_getkey_double(inrec, "BFITZERO", &status);
         PARAMETER_ERROR("BFITZERO")
         drms_setkey_double(outrec, "MOFFSET", bfitzero);
         int i;

         if (status == DRMS_SUCCESS || isnan(bfitzero)) 
         {

         } 
         else 
         {
            for (i = 0; i < pixels; ++i) 
            {
               dat[i] -= (float)bfitzero;
            }

         }
      }

      if (status =       Obs2heliodb((float *)inarr->data, 
	  			   (float *)outarr->data,
				   xpixels, 
				   ypixels, 
				   x0, 
				   y0, 
				   b0 * RADSINDEG, 
				   p * RADSINDEG, 
				   S, 
				   rsun, 
				   rmax,
				   interpolation, 
				   mapcols, 
				   maprows, 
				   longmin_adjusted * RADSINDEG, 
				   longinterval * RADSINDEG,
				   longshift * RADSINDEG, 
				   sinBdelta, 
				   smajor, 
				   sminor, 
				   sangle * RADSINDEG, 
				   xscale, 
				   yscale, 
				   orientation, 
				   mag_correction,
				   velocity_correction, 
				   obs_vr, 
				   obs_vw, 
				   obs_vn, 
				   vsign, 
				   NaN_beyond_rmax,
				   carrStretch,
				   &distP,
				   diffrotA,
				   diffrotB,
				   diffrotC,
				   rRates,
				   sizeof(rRates))) 
      {
         CreateBlankRecord(inrec, outrec, quality);
         fprintf(stderr, "here failure in obs2helio: iRec = %d, status = %d, T_REC = %s, recnum = %lld\n", iRec, status, trecstr, inrec->recnum);
         drms_free_array(inarr);
         drms_free_array(outarr);
         free(orientation);
         continue; // go to next image 
      }

      int xDims = outarr->axis[0], yDims = outarr->axis[1];
      float *odata = outarr->data;

      if (fstats(xDims*yDims, odata, &statMin, &statMax, &statMedn, &statMean, &statSig,
      &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");


      // Output to log rotation rate (deg/day) vs. latitude 
      if (verbflag > 2)
      {
         fprintf(stdout, "latitude\tsidereal rotation rate (deg/day)\n");
         for (row = 0; row < maprows && row < sizeof(rRates); row++)
         {
            fprintf(stdout, "%10f\t%32f\n", rRates[row].lat / RADSINDEG, rRates[row].r);
         }
      }

      if (status =      apodize((float *)outarr->data,
				 b0 * RADSINDEG, 
				 mapcols, 
				 maprows,
				 longmin_adjusted * RADSINDEG,
				 longinterval * RADSINDEG,
				 sinBdelta,
				 apodization, 
				 apinner, 
				 apwidth, 
				 apel, 
				 apx, 
				 apy)) 
      { 
         CreateBlankRecord(inrec, outrec, quality);
         fprintf(stderr, "failure in apodize: iRec = %d, status = %d, T_REC = %s, recnum = %lld\n", iRec, status, trecstr, inrec->recnum);
         drms_free_array(inarr);
         drms_free_array(outarr);
         free(orientation);
         continue; // go to next image
      }

      int i,j,nanoffset, nancount=0;
      float *nanptr=outarr->data;

      for (i=0;i<maprows;i++)
      {
        nanoffset=i*mapcols;
        for (j=0;j<mapcols;j++)
        {
          if (isnan(*(nanptr+nanoffset+j)))
            nancount++;
        }
      }

      // it would be nice to calculate image statistics at this point and write keywords
      //  DATAMIN, DATAMAX, DATARMS, DATAMEAN, DATAMEDN, DATASKEW, DATAKURT
/*
      // Write to the output record
      outarr->bzero=bias;
      outarr->bscale=scale;
      segout->axis[0]=outarr->axis[0]; // required for vardim output
      segout->axis[1]=outarr->axis[1];
*/

      drms_segment_write(segout, outarr, 0);
      drms_copykey(outrec, inrec, "T_REC");
      drms_setkey_int(outrec, "QUALITY", quality);

      drms_setkey_int(outrec, "MAPMMAX", mapped_lmax);
      drms_setkey_int(outrec, "SINBDIVS", sinb_divisions);
      drms_setkey_int(outrec, "MISSVALS", nancount);

      drms_setkey_float(outrec, "CRPIX1", mapcols/2.0 + 0.5);
      drms_setkey_float(outrec, "CRVAL1", mapl0);
//      drms_setkey_float(outrec, "CROTA1", 0.0);
      drms_setkey_float(outrec, "CDELT1", longinterval);
      drms_setkey_float(outrec, "CRPIX2", maprows/2.0 + 0.5);
      drms_setkey_float(outrec, "CRVAL2", 0.0);
      drms_setkey_float(outrec, "CROTA2", 0.0);
      drms_setkey_float(outrec, "CDELT2", sinBdelta);

      drms_setkey_float(outrec, "MAPLGMAX", longmax);
      drms_setkey_float(outrec, "MAPLGMIN", longmin);
      drms_setkey_float(outrec, "MAPBMAX", bmax);
      drms_setkey_float(outrec, "MAPRMAX", rmax);
      drms_setkey_int(outrec, "LGSHIFT", longitude_shift);
      drms_setkey_int(outrec, "INTERPO", interpolation);
      drms_setkey_int(outrec, "MCORLEV", mag_correction);
      drms_setkey_int(outrec, "MOFFSET", mag_offset);
      drms_setkey_int(outrec, "CARSTRCH", carrStretch);
      drms_setkey_float(outrec, "DIFROT_A", diffrotA);
      drms_setkey_float(outrec, "DIFROT_B", diffrotB);
      drms_setkey_float(outrec, "DIFROT_C", diffrotC);


      // these could be links or constant, but set them just in case
      drms_setkey_string(outrec, "CTYPE1", "CRLN_CEA");
      drms_setkey_string(outrec, "CTYPE2", "CRLT_CEA");
      drms_setkey_string(outrec, "CUNIT1", "deg");
      drms_setkey_string(outrec, "CUNIT2", "sinlat");

      drms_copykey(outrec, inrec, "DATE");
      drms_copykey(outrec, inrec, "DATE__OBS");
      drms_copykey(outrec, inrec, "TELESCOP");
      drms_copykey(outrec, inrec, "INSTRUME");
      drms_copykey(outrec, inrec, "ORIGIN");
      drms_copykey(outrec, inrec, "WAVELNTH");
      drms_copykey(outrec, inrec, "CAMERA");

      drms_copykey(outrec, inrec, "T_OBS");
      drms_copykey(outrec, inrec, "CADENCE");
//      drms_copykey(outrec, inrec, "T_REC_step");
      drms_copykey(outrec, inrec, "EXPTIME");
      drms_copykey(outrec, inrec, "CRLT_OBS");
      drms_copykey(outrec, inrec, "CRLN_OBS");
      drms_copykey(outrec, inrec, "HGLN_OBS");
      drms_copykey(outrec, inrec, "CAR_ROT");
      drms_copykey(outrec, inrec, "DATE_OBS");
      drms_copykey(outrec, inrec, "DSUN_OBS");
      drms_copykey(outrec, inrec, "R_SUN");
      drms_copykey(outrec, inrec, "R_SUN_REF");
      drms_copykey(outrec, inrec, "OBS_VR");
      drms_copykey(outrec, inrec, "OBS_VW");
      drms_copykey(outrec, inrec, "OBS_VN");

// write the statistical results
        drms_setkey_int(outrec, "TOTVALS", xDims*yDims);
        drms_setkey_int(outrec, "DATAVALS", statNgood);
        drms_setkey_int(outrec, "MISSVALS", xDims*yDims-statNgood);
        drms_setkey_double(outrec, "DATAMIN", statMin);
        drms_setkey_double(outrec, "DATAMAX", statMax);
        drms_setkey_double(outrec, "DATAMEDN", statMedn);
        drms_setkey_double(outrec, "DATAMEAN", statMean);
        drms_setkey_double(outrec, "DATARMS", statSig);
        drms_setkey_double(outrec, "DATASKEW", statSkew);
        drms_setkey_double(outrec, "DATAKURT", statKurt);

      dsignout=vsign*drms_getkey_double(inrec, "DATASIGN", &status);
      if (status != DRMS_SUCCESS) dsignout=vsign;
      dsignout/=fabs(dsignout);
      drms_setkey_int(outrec, "DATASIGN", (int)dsignout);

      tnow = (double)time(NULL);
      tnow += UNIX_epoch;
      drms_setkey_time(outrec, "DATE", tnow);

      drms_close_record(outrec, DRMS_INSERT_RECORD);

      if (verbflag > 1) 
      {
        wt=getwalltime();
        ct=getcputime(&ut, &st);
        fprintf(stdout, 
                "record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                iRec,
                wt-wt1,
                ct-ct1);

      }

      drms_free_array(inarr);
      drms_free_array(outarr);
      free(orientation);
      nsuccess++;

   } // end loop

   if (inRecSet)
      drms_close_records(inRecSet, DRMS_FREE_RECORD);

   if (error == 1)
      drms_close_record(outrec, DRMS_FREE_RECORD);

   wt=getwalltime();
   ct=getcputime(&ut, &st);

   if (verbflag) 
   {
      printf("number of records processed = %d\n", nsuccess);
      fprintf(stdout, "total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
              wt-wt0,
              ct-ct0);
      if (!error)
         printf("module %s successful completion\n", cmdparams.argv[0]);
   }

   return 0;

} // end DoIt


static void CheckO(const char *orientation, V2HStatus_t *stat)
{
   /* check for legal orientation string */
   static char o[MAXLEN];
   char *c = o;

   if(!orientation)
   {
      *stat = V2HStatus_MissingParameter;
   } 
   else if (4 != sscanf(orientation, "%[NS]%[EW]%[NS]%[EW]", c++, c++, c++, c++))
   {
      *stat = V2HStatus_IllegalOrientation;
   }
   else if ((o[0] == o[2]) && (o[1] == o[3]))
   { 
      *stat = V2HStatus_IllegalOrientation;
   }
   else if ((o[0] !=o [2]) && (o[1] != o[3])) 
   {
      *stat = V2HStatus_IllegalOrientation;
   }
   else 
   {
      *stat = V2HStatus_Success;   
   }
}
