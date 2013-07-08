 /* static char rcsid[] = "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/mdisynop.c,v 1.1 2013/07/08 17:36:04 yliu Exp $"; */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "jsoc_main.h"
#include "astro.h" 
#include "drms_dsdsapi.h"
#include "soi_ephem.h"
#include "/home0/yliu/cvs/JSOC/proj/myproj/apps/src/fstats.c"
//#include "/home0/yliu/cvs/JSOC/proj/libs/astro/heliographic_coords.c"
#include "/home/wso/src/libastro.d/solephem.c"

#define C1      0.0757644 /* days/degree at 27.27527 */
#define C2      92353.9357 /* day of 1853:11:09_22h:27m:24s */
#define RADSINDEG       (PI/180)

/* cmd-line parameters */
#define kRecSetIn       "in"
#define kSeriesOut      "out"
#define kSegIn          "segin"
#define kSegOut         "segout"
#define kNOTSPECIFIED   "not specified"

#define EQPOINTS               "EQPOINTS"
#define HALFWIN                "HALFWIN"
#define LOSFLAG                "l"
#define FORCEFLAG              "f" /* force - accept data without checking version */
#define DLOGFLAG               "d" /* debug log - verbose stdout */
#define NOISESIGMA             "NOISESIGMA"  /* std dev of noise */
#define MAXNOISEADJ            "MAXNOISEADJ" /* maximum adjustment due to increases
					      * with latitude of the noise level
					      * (radial data only) */
#define MINOUTLIERPTS          "MINOUTLIERPTS" /* minimum number of points that must
						  exist before outliers can be 
					          discarded when calulating summary
					          statistics */
#define DPC_OBSR               "DPC_OBSR"
#define FD_Magnetogram_Sum     "FD_Magnetogram_Sum"
#define FD_Magnetogram         "FD_Magnetogram"
#define HWNWIDTH               "HWNWIDTH"

/* keywords in the output data */
#define kCARSTRCH              "CARSTRCH"
#define kDIFROT_A              "DIFROT_A"
#define kDIFROT_B              "DIFROT_B"
#define kDIFROT_C              "DIFROT_C"
#define CR                     "CR"

const float kOutScale = 0.1;
const int kOutBPP = 16;
const float kNOISE_EQ = 20.0;
const float sensAdjVal = 1.8;

enum MagAggregate_t
{
     kOneMinuteValue = 0,
     kFiveMinuteValue = 1
};

float gMagAgg[] =
{
   1.0,
   2.75
};

char *module_name = "mdisynop";

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}

/* Describes a column of magnetogram data.  Fields other than slot refer to 
 * the magnetogram in which this data originated.  For example, dist is the 
 * distance (in pixels) between the CM - center of the magnetogram and the 
 * magnetogram's column of data.
 */

typedef struct MagCol_struct
{
  float dist;     /* Distance, in magnetogram pixels, from CM to column. */
  float *datacol; /* Ptr to column of data elonging to one magnetogram */
  float equivPts; /* 1-minute noise-equivalent points of this column */
  int ds;         /* ds number, for id purposes */
  int sn;         /* sn number, for id purposes */
  int col;        /* original magnetogram column */
  char valid;     /* 1 if in use */
} MagCol_t;

/*
double local_earth_meridian_crossing(double L, int cr);
double earth_B(TIME t);
*/
int CalcSynCols(int start,
		int end,
		int incr,
		MagCol_t **smc,
		float *synop,
		int *wt,
		char *ww,
		float *epts,
		float *losSynop,
		int *len,
		float center,
		float nEquivPtsReq, 
		int los,
		int radialFound,
		float sensAdj,
		float noiseS, 
		float nsig,
		float maxNoiseAdj, 
		int minOutPts,
		int dlog);

int magStats(float **val, int npts, double sum, int outThreshold, double *avg, float *med);
void SortMagCols(MagCol_t *mc, int nelem);
void FreeMagColsData(int start, int end, int incr, int *n, MagCol_t **mc);
TIME CarringtonTime(int crot, double L);

ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kRecSetIn, "", "Input data records."},
   {ARG_STRING,  kSeriesOut, "", "Output data series."},
   {ARG_STRING,  kSegIn, kNOTSPECIFIED, ""},
   {ARG_STRING,  kSegOut, kNOTSPECIFIED, ""},
   {ARG_INT, CR, "", "", ""},
   {ARG_FLOAT, "nsig", "3.0", "", ""},
   {ARG_INT, "MAPMMAX", "1800", "", ""},
   {ARG_INT, "SINBDIVS", "720", "", ""},
   {ARG_FLOAT, "MAPLGMAX", "90.0", "", ""},
   {ARG_FLOAT, "MAPLGMIN", "-90.0", "", ""},
   {ARG_FLOAT, EQPOINTS, "20.0", "", ""},
   {ARG_FLOAT, HALFWIN, "30.0", "", ""},
   {ARG_FLOAT, "center", "0.0", "", ""},
   {ARG_INT, "checkqual", "1", "", ""},
   {ARG_STRING, "qualmask", "0x0", ""},
   {ARG_FLAG, LOSFLAG, "0", "-1", "1"},
   {ARG_FLAG, FORCEFLAG, "0", "-1", "1"},
   {ARG_FLAG, DLOGFLAG, "0", "-1", "1"},
   {ARG_FLOAT, NOISESIGMA, "2.0", "", ""},
   {ARG_FLOAT, MAXNOISEADJ, "3.0", "", ""},
   {ARG_INT, MINOUTLIERPTS, "4", "", ""},
   {ARG_END}
};

int DoIt(void)
{

  int newstat = 0;
  int status = DRMS_SUCCESS;
  float *synop;
  int *wt;
  char *ww;
  float *epts = NULL;
  float *losSynop = NULL;
  float *ps, *pm;
  char *odir;
  struct stat stbuf;
  float nsig;
  int cr, mapmmax, sinbdivs;
  int checkqual;
  unsigned int quality, qualvld, qualmask, qualmask1, qualmask2, qualmask3;
  float lgmax, lgmin;
  float center;
  float halfWindow;
  int len[2]; /* 0 - x pixels, 1 - y pixels */
  int ngood, idx;
  int started, ended;
  int firstidx, lastidx;
  int mapmidcol, synmidcol;
  int cols2right, cols2left;
  int mrc, mlc, src, slc;
  int mapcols, col, row, syncol;
  double maprightct, mapleftct;
  double synleftct, synrightct;
  double synstart; /* CT, in degrees, of the left-most pixel of the synoptic chart to output */
  double synend; /* CT, in degrees, of the right-most pixel of the synoptic chart to output */
  double synstep; /* degrees per pixel */
  char tstr[64];
  int mag_num = 0;
  TIME tstart, tstop, trot, tearth;
  TIME tdiff, magtime;
  double r,b,l,vr,vn,vw;
  int y,m,d;
  TIME tmod, delta_T = 0.0;
  double bearth;
  double carrtime;
  char dsname[1024], key[1024], logfn[1024], synfn[1024], wtfn[1024], losFn[1024];
  char eptsfn[1024];
  int ds, nds, noutRec;
  int fsn, lsn, sn;
  int i,j;
  double statMin, statMax, statMedn, statMean, statSig, statSkew, statKurt;
  int statNgood;
  double eph[30];

  struct {
      int recno;
      double mapct; /* Carrington time of center pixel of remapped image (in degrees) MINUS center. */
      double mapCM; /* Carrington time of CM of remapped image (in degrees). */
      float mapdev; /* deviation between center pixel and either end of the map (in degrees) */
      float tmin;
      float tmax;
      TIME tobs;
          } imrec[1024];

  int nStackMags;
  float nEquivPtsReq;
  MagCol_t **sortedMagCol = NULL; 

  long long modVers = soi_vers_num;
  long long imgVers = 0;
  int radialFound = 0;
  int losFound = 0;
  int carrStretch = -1;
  float diffrotA = DRMS_MISSING_FLOAT;
  float diffrotB = DRMS_MISSING_FLOAT;
  float diffrotC = DRMS_MISSING_FLOAT;

  const char *synopFileName = NULL;
  const char *wtFileName = NULL;
  const char *eptsFileName = NULL;

  int los;
  int force;
  int dlog;
  float noiseS = cmdparams_get_float(&cmdparams, NOISESIGMA, &status);      /* noise sigma */
  float maxNoiseAdj = cmdparams_get_float(&cmdparams, MAXNOISEADJ, &status); /* maximum adj due to increases
							* with latitude of the noise level
							* (radial data only) */
  int minOutPts = cmdparams_get_int(&cmdparams, MINOUTLIERPTS, &status);     /* minimum number of points that 
							* must exist before outliers can be
							* discarded when calulating summary
							* statistics */
  int sensAdjDone = -1;
  float sensAdj;

  int nxtSyncol = -1; /* the first synop column that hasn't been calculated */
  int calcsynret = 0;
  delta_T = sscan_time("1977.01.01_00:00:00_TAI") - sscan_time("1601.01.01_00:00:00_UT");
           // Difference between DRMS_EPOCH_S (1977.01.01_00:00:00_TAI) and WSO_EPOCH_S (1601.01.01_00:00:00_UT)
           // CarringtonTime converts time and Carrington coordinates based on Phil's code ctimes.c
           // (/home/wso/src/misc/ctimes.c) which is based on WSO obs. Thus this time difference needs to
           // be corrected to consist with DRMS time.

  cr = cmdparams_get_int(&cmdparams, CR, &status);
  nsig = cmdparams_get_float(&cmdparams, "nsig", &status);
  mapmmax = cmdparams_get_int(&cmdparams, "MAPMMAX", &status);
  sinbdivs = cmdparams_get_int(&cmdparams, "SINBDIVS", &status);
  checkqual = cmdparams_get_int(&cmdparams, "checkqual", &status);
  qualmask = strtoul(cmdparams_get_str(&cmdparams, "qualmask", &status), (char **) 0, 0);
  center = cmdparams_get_float(&cmdparams, "center", &status);
  halfWindow = cmdparams_get_float(&cmdparams, HALFWIN, &status);
  lgmax = cmdparams_get_float(&cmdparams, "MAPLGMAX", &status);
  lgmin = cmdparams_get_float(&cmdparams, "MAPLGMIN", &status);

  qualmask3 = strtoul("0x00000200", (char **) 0, 0);
  qualmask2 = strtoul("0x00000201", (char **) 0, 0);
  qualmask1 = strtoul("0x00000000", (char **) 0, 0);

  if (lgmax <= lgmin) 
    {
      printf("MAPLGMAX must be greater than MAPLGMIN\n");
      DIE("MAKES_NO_SENSE");
    }

  nEquivPtsReq = cmdparams_get_float(&cmdparams, EQPOINTS, &status);

/* synoptic char dimensions*/
  len[0] = rint(360.0 / (lgmax - lgmin) * mapmmax);
  len[1] = 2 * sinbdivs;

  synstep = 360.0 / len[0];
  synstart = (cr - 1) * 360.0;
  synend = cr * 360.0 - synstep;

  nStackMags = rint(2 * halfWindow * 1.5); /* consecutive magnetograms are 
  					  * shifted about 1 degree apart */

//  tstart = meridian_crossing(360.0, cr);
//  tstop = meridian_crossing(0.0, cr);
//  trot = meridian_crossing(180.0, cr);

//  tstart = HeliographicTime(cr, 360.0);
//  tstop = HeliographicTime(cr, 0.0);
//  trot = HeliographicTime(cr, 180.0);
//  tearth = DRMS_MISSING_TIME;
//  bearth = DRMS_MISSING_FLOAT;

  tstart = CarringtonTime(cr, 360.0);
  tstop = CarringtonTime(cr, 0.0);
  trot = CarringtonTime(cr, 180.0);

//  sprint_at(tstr, tstart);
//printf("TSTART=%s\n", tstr);
//  sprint_at(tstr, tstop);
//printf("TSTOP=%s\n", tstr);

  los = (cmdparams_isflagset(&cmdparams, LOSFLAG) != 0);
  force = (cmdparams_isflagset(&cmdparams, FORCEFLAG) != 0);
  dlog = (cmdparams_isflagset(&cmdparams, DLOGFLAG) != 0);

printf("\n######### FIRST PASS here #################\n");

/* Determine how many of the input series are valid.  Store this in idx. */
/* Also, determine the minimum and maximum CT covered by each image. */

  if (modVers < 0)
    {
       modVers = -modVers;
    }
    char *inRecQuery, *outRecQuery;
    inRecQuery = (char *)params_get_str(&cmdparams, kRecSetIn);
    outRecQuery = (char *)params_get_str(&cmdparams, kSeriesOut);
//    char *inRecQuery = cmdparams_get_str(&cmdparams, "in", &status);
//    char *outRecQuery = cmdparams_get_str(&cmdparams, "out", &status);
    DRMS_RecordSet_t *inRD, *outRD;
//    DRMS_RecordSet_t *outRec;
    DRMS_Record_t *outRec;    

    inRD = drms_open_records(drms_env, inRecQuery, &status);
    if (status || inRD->n == 0)
        DIE("No input data found");
    nds = inRD->n;
    idx = 0;
    for (ds=0; ds < nds; ds++) 
      {
        double clong;
        double cmLong;
        double mapct = 0.0;
        double mapctnew;
        double mapCM = 0.0;
        double mapCMnew;
        double remapLgmax = 0.0;
        double remapLgmin = 0.0;
        int orot;
        char t_obs[64];
        char *imgVersStr = NULL;
        int rkey = 0;
        int csKey = 0;
        float csCoeffKey;
        int sensAdjKey;
        DRMS_Record_t *inRec;
        inRec = inRD->records[ds];

        if (checkqual) 
        {
        quality = drms_getkey_int(inRec, "QUALITY", &status);
        if (drms_ismissing_int(quality))
        { 
          printf("  Can't find QUALITY in log file; quality check disabled\n"); 
          continue;
        }
          else if ((quality != qualmask1) & (quality != qualmask2) & (quality != qualmask3)) 
          {
              printf("  Bad QUALITY = 0x%08x; rejected\n", quality);
//printf("check=%d, status=%d\n", checkqual, status); 
              continue;
           }
        }

/*
        if (checkqual) 
        {
        quality = drms_getkey_int(inRec, "QUALITY", &status);
        qualvld = drms_getkey_int(inRec, "QUAL_VLD", &status);

	if (drms_ismissing_int(quality) || drms_ismissing_int(qualvld)) 
	  printf("  Can't find QUALITY in log file; quality check disabled\n"); 
	  else if (quality & qualvld & qualmask) 
          {
	      printf("  Bad QUALITY = 0x%08x; rejected\n", quality); 
	      continue;
	   }
        }
*/
      /* Don't create synoptic charts from obsolete data. */
      if (!force)
        {
	 imgVers = 0;
	 imgVersStr = drms_getkey_string(inRec, "BLDVER18", &status);

	  if (imgVersStr)
	    {
	      char *endptr = NULL;
	      imgVers = strtoll(imgVersStr, &endptr, 10);
	      if (*endptr != '\0' || endptr == imgVersStr)
	        {
	           /* Invalid build version */
	           printf("  Invalid build version information.\n"); 
	           printf("    Rejecting ds=%d.\n", ds); 
	           continue;
	        }
	    }
	    else
	      {
	        /* Missing build number, bail on this image. */
	        printf("    Rejecting ds=%d \n", ds); 
	        continue;
	      }

	    if (imgVers < 0)
	      {
	        imgVers = -imgVers;
	      }

	    if (imgVers < modVers)
	      { 
	        printf("  Attempt to use obsolete image (img lev18 vers %lld, mod lev18 vers %lld).\n", 
		    imgVers, modVers);
	        printf("    Rejecting ds=%d.\n", ds); 
	        continue;
	       }
          }

            /* Ensure that all images used are either radial values or 
            * line-of-sight values, not a mixture of the two. */
          rkey = drms_getkey_int(inRec, "FDRADIAL", &status);
          if (rkey > 0)
            {
	       printf("  Found radial keyword %d, ds=%d.\n", rkey, ds); 
	       if (losFound)
	         {
	           printf("  Attempt to use a mixture of radial and line-of-sight images.\n"); 
	           printf("    Rejecting ds=%d.\n", ds); 
	           continue;
	         }
	       else
	         {
	           radialFound = 1;
	         }
             }
           else
             {
	       if (radialFound)
	         {
	           printf("  Attempt to use a mixture of radial and line-of-sight images.\n"); 
	           printf("    Rejecting ds=%d.\n", ds); 
	           continue;
	         }
	       else
	         {
	           losFound = 1;
	         }
             }

             csKey = drms_getkey_int(inRec, kCARSTRCH, &status);
             csKey = (drms_ismissing_int(csKey)) ? 0 : csKey;

             if (idx == 0)
               {
	         carrStretch = csKey;
               }
             else
               {
	         if (csKey != carrStretch)
	           {
	             printf("  Attempt to mix carr-streched and non-carr-stretched images.\n");
	             printf("    Rejecting ds=%d.\n", ds); 
	             continue;
	           }
               }

              if (carrStretch > 0)
                {
	          csCoeffKey = (float)drms_getkey_double(inRec, kDIFROT_A, &status);
	          if (idx == 0)
	            {
	              diffrotA = (float)csCoeffKey;
	            }
	          else
	            {
	              if (fabsf(csCoeffKey - diffrotA) > 0.001)
	                {
		          printf("  Attempt to use inconsistent carr stretch parameters.\n");
		          printf("    Rejecting ds=%d.\n", ds); 
		          continue;
	                }
	            }

	           csCoeffKey = (float)drms_getkey_double(inRec, kDIFROT_B, &status);
	           if (idx == 0)
	             {
	               diffrotB = (float)csCoeffKey;
	             }
	           else
	             {
	               if (fabsf(csCoeffKey - diffrotB) > 0.001)
	                 {
		           printf("  Attempt to use inconsistent carr stretch parameters.\n");
		           printf("    Rejecting ds=%d.\n", ds); 
		           continue;
	                 }
	             }

	             csCoeffKey = (float)drms_getkey_double(inRec, kDIFROT_C, &status);
	             if (idx == 0)
	               {
	                 diffrotC = (float)csCoeffKey;
	               }
	             else
	               {
	                 if (fabsf(csCoeffKey - diffrotC) > 0.001)
	                   {
		             printf("  Attempt to use inconsistent carr stretch parameters.\n");
		             printf("    Rejecting ds=%d.\n", ds); 
		             continue;
	                    }
	               }
                   }

                   sensAdjKey = drms_getkey_int(inRec, "SENSADJ", &status);
                   if (sensAdjKey < 0)
                     {
	               sensAdjKey = 0;
                     }

                   if (sensAdjDone == -1)
                     {
	               sensAdjDone = sensAdjKey;
	               sensAdj = sensAdjDone ? sensAdjVal : 1.0;
                     }
                   else
                     {
	               /* ensure all images are the same */
	               if (sensAdjKey != sensAdjDone)
	                 { 
	                   printf(" Attempt to use a mixture of sensitivity adjusted and non-adj imgs.\n");
	                   printf("    Rejecting ds=%d.\n", ds); 
	                   continue;
	                 }
                      }

                      clong = drms_getkey_double(inRec, "CRVAL1", &status);
                      cmLong = drms_getkey_double(inRec, "CRLN_OBS", &status);
                      orot = drms_getkey_int(inRec, "CAR_ROT", &status);
                      if (drms_ismissing_float(clong) || drms_ismissing_int(orot)) 
                        {
	                  printf("  Bad MAP_L0 or OBS_CR in header; skipped"); 
	                  continue;
                        }

                      mapctnew = (double)(orot) * 360.0 - clong - center;
                      mapCMnew = (double)(orot) * 360.0 - cmLong - center;
                      if (mapCMnew < mapCM)
                        {
	                  printf("  Source image (ds=%d) has a smaller CT than previous img; skipped.\n", ds);
	                  continue;
                        }
                      else
                        {
	                  mapct = mapctnew;
	                  mapCM = mapCMnew;
                        }

                   remapLgmax = drms_getkey_double(inRec, "MAPLGMAX", &status);
                   remapLgmin = drms_getkey_double(inRec, "MAPLGMIN", &status);
                   imrec[idx].mapdev = (remapLgmax - remapLgmin) / 2.0;
                   imrec[idx].recno = drms_getkey_int(inRec, "I_DREC", &status);
                   imrec[idx].mapct = mapct;
                   imrec[idx].mapCM = mapCM;
//                   imrec[idx].ds = ds;
                   strcpy(t_obs, drms_getkey_string(inRec, "T_OBS", &status));
                   imrec[idx].tobs = sscan_time(t_obs);
//                   imrec[idx].tmin = MAX(mapCM - halfWindow, mapct + center - imrec[idx].mapdev);
//                   imrec[idx].tmax = MIN(mapCM + halfWindow, mapct + center + imrec[idx].mapdev);
                   ++idx;

               } /* for each dataset */

/* Every image pixel has a Carrington time.  Calculate the CT range for each 
 * magnetogram.
 */
    ngood = idx;

    for (idx = 0; idx < ngood; ++idx) 
      {
         imrec[idx].tmin = MAX(imrec[idx].mapCM - halfWindow, imrec[idx].mapct + center - imrec[idx].mapdev);
         imrec[idx].tmax = MIN(imrec[idx].mapCM + halfWindow, imrec[idx].mapct + center + imrec[idx].mapdev);
    }

printf("\n######### SECOND PASS #################\n"); 

      if (!(wt = (int *) calloc(len[0],4)))
        DIE("MALLOC_FAILED");
      if (!(ww = (char *) malloc(len[0]*len[1]*2)))
        DIE("MALLOC_FAILED");
      if (!(epts = (float *) malloc(len[0] * len[1] * sizeof(float))))
        DIE("MALLOC_FAILED");
      if (!(synop = (float *) malloc(len[0]*len[1]*4)))
        DIE("MALLOC_FAILED");
      if (los && radialFound && !(losSynop = (float *)malloc(len[0] * len[1] * sizeof(float))))
        {
          DIE("MALLOC_FAILED");
        }
      sortedMagCol = (MagCol_t **)malloc(sizeof(MagCol_t *) * len[0]);
      if (!sortedMagCol)
        {
          DIE("MALLOC_FAILED");
        }
      memset(sortedMagCol, 0, sizeof(MagCol_t *) * len[0]);
    idx = 0;
    started = 0;
    ended = 0;
    nxtSyncol = len[0] - 1;

/* Step through successive magnetograms that overlap the synoptic chart's CT range.  
 * Stop as soon as the focal magnetogram ceases to overlap this range.
 */

    for (ds=0; ds < nds; ds++) 
      {
        int recno;
        int col;
        float magtype = 0.0;
        char *mType = NULL;
        float equivPts = 0.0;
        char t_mag[64];
        TIME tmag;
        DRMS_Record_t *inRec;
        DRMS_Segment_t *inSeg = NULL;
        DRMS_Array_t *inArray;
        inRec = inRD->records[ds];

        strcpy(t_mag, drms_getkey_string(inRec, "T_OBS", &status));
        tmag = sscan_time(t_mag);
        if (tmag != imrec[idx].tobs) continue;

        /* Skip magnetograms whose CT range lies below the CT range of the synoptic chart. */

        if (imrec[idx].tmax < synstart) 
          {
	    ++idx;
	    continue;
         }
        /* Skip magnetograms whose CT range lies above the CT range of the synoptic chart. */
        if (imrec[idx].tmin > synend) 
          {
	    if (!ended) 
             {
	        ended = 1;
	        lastidx = idx;
	      }
	    ++idx;
	    continue;
         }

        if (!started) 
          {
	    started = 1;
	    firstidx = idx;
          }


        /* truncate the magnetogram's CT range at synstart or synend if necessary */
        maprightct = MAX(synstart, imrec[idx].tmin);
        mapleftct = MIN(synend, imrec[idx].tmax);
        mapcols = drms_getkey_int(inRec, "MAPMMAX", &status) + 1;
        mapmidcol = rint((mapcols - 1.0) / 2 + center / synstep);
        cols2right = rint((imrec[idx].mapct - maprightct) / synstep);
        cols2left = rint((imrec[idx].mapct - mapleftct) / synstep);
        mrc = mapmidcol + cols2right;
        mlc = mapmidcol + cols2left;
        synmidcol = rint((synend - imrec[idx].mapct) / synstep); /* NOT the middle column of the
								* synoptic chart, but essentially
								* the column of the synoptic chart
								* that corresponds to the middle
								* column of the remapped img. */

        /* src - right column of synoptic chart that matches with mrc */
        /* slc - left column of synoptic chart that matches with mlc */
        src = synmidcol + cols2right;
        slc = synmidcol + cols2left;

        /* If src < nxtSyncol, then it is okay to calculate all synoptic chart pixels
         * from src + 1 to nxtSynCol. */

        if (src < nxtSyncol)
          {
	    calcsynret = CalcSynCols(nxtSyncol,
	      			    src + 1,
				    -1,
				    sortedMagCol, 
				    synop,
				    wt, 
				    ww,
				    epts,
				    losSynop,
				    len,
				    center,
				    nEquivPtsReq, 
				    los,
				    radialFound,
				    sensAdj,
				    noiseS,
				    nsig,
				    maxNoiseAdj, 
				    minOutPts,
				    dlog);
	    FreeMagColsData(nxtSyncol, src + 1, -1, wt, sortedMagCol);

	    nxtSyncol = src;
	    if (calcsynret) DIE("Fail to sort data");
          }

      /* For each column of the synoptic chart, push the corresponding input magnetogram's 
       * column of data to the stack at that column.  The stack is an array of 2D 
       * input images that will be averaged to make the synoptic chart.  The dimensions 
       * of these images match the final dimensions of the synoptic chart.  Upon 
       * completion of all input series, wt is the depth of the stack at each synoptic column. 
       */
 
        magtype = drms_getkey_float(inRec, "EXPTIME", &status);
        equivPts = 0.0;

        if (magtype == 300.0)
          {
            equivPts = gMagAgg[kFiveMinuteValue];
          }
        else if (magtype == 30.0)
          {
            equivPts = gMagAgg[kOneMinuteValue];
          }
        else
          {
             equivPts = gMagAgg[kOneMinuteValue];
          }

        inSeg = drms_segment_lookupnum(inRec, 0);
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);

        if (status)
            {
            printf(" No data file found, status=%d\n", status);
            drms_free_array(inArray);
            continue;
            }

        for (col = mlc; col <= mrc; ++col) 
          {
	    MagCol_t mMagCol;

	    mMagCol.dist = ((col - mapmidcol) * synstep + imrec[idx].mapct) - imrec[idx].mapCM;
	    mMagCol.datacol = (float *)malloc(sizeof(float) * len[1]);
	    mMagCol.equivPts = equivPts;
	    mMagCol.ds = ds;
	    mMagCol.col = col;
	    mMagCol.valid = 1;
	   
	    pm = (float *)inArray->data + col;
	    ps = mMagCol.datacol;

	    for (row = 0; row < len[1]; ++row) 
	      {
	        *ps = *pm;
	        pm += mapcols;
	        ps++;
	      }
	   
	    syncol = col + slc - mlc;

	    if (!sortedMagCol[syncol])
	      {
		sortedMagCol[syncol] = (MagCol_t *)malloc(sizeof(MagCol_t) * nStackMags);
		if (!sortedMagCol[syncol]) DIE("MALLOC_FAILED");
		memset(sortedMagCol[syncol], 0, sizeof(MagCol_t) * nStackMags);
	      }

	    if (wt[syncol] > 0 && (wt[syncol] % nStackMags == 0))
	      {
	        int frames = 1 + (wt[syncol] / nStackMags);
	        MagCol_t *mc = (MagCol_t *)malloc(sizeof(MagCol_t) * frames * nStackMags);
	        memcpy(mc, 
		       sortedMagCol[syncol], 
		       sizeof(MagCol_t) * (frames - 1) * nStackMags);
	        free(sortedMagCol[syncol]);
	        sortedMagCol[syncol] = mc;
	      }

	    sortedMagCol[syncol][wt[syncol]] = mMagCol;	   
	    wt[syncol]++;

       } /* for each series */
    drms_free_array(inArray);
    ++idx;
  }

/* finish off synoptic chart - the above algorithm */
  if (nxtSyncol >= 0)
    {
       calcsynret = CalcSynCols(nxtSyncol,
  			        0,
			        -1,
			        sortedMagCol, 
			        synop,
			        wt, 
			        ww,
			        epts,
			        losSynop,
			        len,
			        center,
			        nEquivPtsReq, 
			        los,
			        radialFound,
			        sensAdj,
			        noiseS,
			        nsig,
			        maxNoiseAdj, 
			        minOutPts,
			        dlog);

       FreeMagColsData(nxtSyncol, 0, -1, wt, sortedMagCol);
       nxtSyncol = -1;
       if (calcsynret) DIE("Fail to sort data");
    }

    if (!ended) lastidx = ngood-1;
    magtime = imrec[(int)(ngood/2)].tobs;

/*
    tdiff = 1e99;
    for (idx=firstidx; idx<=lastidx; ++idx) 
        {
          if (fabs(imrec[idx].tobs - trot) < tdiff) 
            {
	      tdiff = fabs(imrec[idx].tobs - trot);
	      magtime = imrec[idx].tobs;
            }
        }
*/
  
  if (fstats(len[0]*len[1], synop, &statMin, &statMax, &statMedn, &statMean, &statSig,
      &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed...\n");

    outRD = drms_create_records(drms_env, 1, outRecQuery, DRMS_PERMANENT, &status);
    if (status)
        DIE("Output record not created");
        DRMS_Segment_t *outSeg;
        DRMS_Array_t *outArray;
        int outDims[2] = {len[0], len[1]};
        outRec = outRD->records[0];

	// ww (first segment)
	outSeg = drms_segment_lookupnum(outRec, 0);
	outArray = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, ww, &status);
	status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file");
        drms_free_array(outArray);


	// synop (second segment)
	outSeg = drms_segment_lookupnum(outRec, 1);
        outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, synop, &status);
	status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file");
        drms_free_array(outArray);

	// epts (third segment)
	outSeg = drms_segment_lookupnum(outRec, 2);
	outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, epts, &status);
	status = drms_segment_write(outSeg, outArray, 0);
        if (status)
          DIE("problem writing file");
        drms_free_array(outArray);

/* write the keywords */

//for reserved keywords

        sprint_ut(tstr, CURRENT_SYSTEM_TIME);
        sscanf(tstr, "%d.%d.%d", &y, &m, &d);
        sprintf(tstr, "%04d-%02d-%02d", y, m, d);
        drms_setkey_string(outRec, "DATE", tstr);
        drms_setkey_string(outRec, "BUNIT", "Gauss");

//for image coordinate mapping keywords

        drms_setkey_string(outRec, "CTYPE1", "CRLN-CEA");
        drms_setkey_string(outRec, "CTYPE2", "CRLT-CEA");
        i=len[0]/2+0.0;
        drms_setkey_double(outRec, "CRPIX1", i);
           // origin is at the left corner of the first pixel
        l=len[1]/2+0.0;
        drms_setkey_double(outRec, "CRPIX2", l);
           // origin is at the left corner of the first pixel
        carrtime = (cr-1)*360.0 + 180.0 - 0.5 * 360.0/len[0];;
        drms_setkey_double(outRec, "CRVAL1", carrtime);
        i=0.0;
        drms_setkey_double(outRec, "CRVAL2", i);
        l=-360.0/len[0];
        drms_setkey_double(outRec, "CDELT1", l);
        drms_setkey_double(outRec, "CDELT2", 1.0/sinbdivs);
        drms_setkey_string(outRec, "CUNIT1", "degree");
        drms_setkey_string(outRec, "CUNIT2", "Sine Latitude");

//HMI observables keywords
        sprint_at(tstr, trot - delta_T);
//        drms_setkey_string(outRec, "T_REC", tstr);
// image statistics

        i = len[0]*len[1];
        drms_setkey_int(outRec, "TOTVALS", i);
        drms_setkey_int(outRec, "DATAVALS", statNgood);
        i = len[0]*len[1]-statNgood;
        drms_setkey_int(outRec, "MISSVALS", i);
        drms_setkey_double(outRec, "DATAMIN", statMin);
        drms_setkey_double(outRec, "DATAMAX", statMax);
        drms_setkey_double(outRec, "DATAMEDN", statMedn);
        drms_setkey_double(outRec, "DATAMEAN", statMean);
        drms_setkey_double(outRec, "DATARMS", statSig);
        drms_setkey_double(outRec, "DATASKEW", statSkew);
        drms_setkey_double(outRec, "DATAKURT", statKurt);

//Synoptic map keywords
        drms_setkey_string(outRec, "T_OBS", tstr);
        drms_setkey_string(outRec, "T_ROT", tstr);
        sprint_at(tstr, tstart - delta_T);
        drms_setkey_string(outRec, "T_START", tstr);
        sprint_at(tstr, tstop - delta_T);
        drms_setkey_string(outRec, "T_STOP", tstr);
//        sprint_at(tstr, tearth);
//        drms_setkey_string(outRec, "T_EARTH", tstr);
        drms_setkey_int(outRec, "CAR_ROT", cr);
        drms_setkey_double(outRec, "CARRTIME", carrtime);
        solephem(trot, eph);
//        HeliographicLocation(trot, &i, &l, &b);
        drms_setkey_double(outRec, "B0_ROT", eph[9]/RADSINDEG);
        solephem(tstart, eph);
//        HeliographicLocation(tstart, &i, &l, &b);
        drms_setkey_double(outRec, "B0_FRST", eph[9]/RADSINDEG);
        solephem(tstop, eph);
//        HeliographicLocation(tstop, &i, &l, &b);
        drms_setkey_double(outRec, "B0_LAST", eph[9]/RADSINDEG);
//        drms_setkey_double(outRec, "EARTH_B0", bearth);
        l=(cr-1)*360.0;
        drms_setkey_double(outRec, "LON_FRST", l);
        l=cr*360.0-360.0/len[0];
        drms_setkey_double(outRec, "LON_LAST", l);
        l=-360.0/len[0];
        drms_setkey_double(outRec, "LON_STEP", l);
        l=center;
        drms_setkey_double(outRec, "W_OFFSET", l);
        drms_setkey_string(outRec, "W_WEIGHT", "Even");
        i=lastidx-firstidx+1;
        drms_setkey_int(outRec, "IMG_NUM", i);
        sprint_at(tstr, imrec[firstidx].tobs);
        drms_setkey_string(outRec, "IMG_FRST", tstr);
        sprint_at(tstr, imrec[lastidx].tobs);
        drms_setkey_string(outRec, "IMG_LAST", tstr);
        sprint_at(tstr, magtime);
        drms_setkey_string(outRec, "IMG_ROT", tstr);
        drms_setkey_float(outRec, HWNWIDTH, halfWindow);
        drms_setkey_float(outRec, EQPOINTS, nEquivPtsReq);
        drms_setkey_float(outRec, "NSIGMA", nsig);

// save the differential rotation correction parameters 
        if (carrStretch > 0)
          {
            drms_setkey_int(outRec, kCARSTRCH, carrStretch);
            drms_setkey_float(outRec, kDIFROT_A, diffrotA);
            drms_setkey_float(outRec, kDIFROT_B, diffrotB);
            drms_setkey_float(outRec, kDIFROT_C, diffrotC);
          }

/*
        i=len[0]/2;
        drms_setkey_int(outRec, "CRPIX1", i);
        l=len[1]/2+0.5;
        drms_setkey_double(outRec, "CRPIX2", l);
        carrtime = cr*360.0 - 180.0;
        drms_setkey_double(outRec, "CRVAL1", carrtime);
        i=0;
        drms_setkey_int(outRec, "CRVAL2", i);
        l=-360.0/len[0];
        drms_setkey_double(outRec, "CDELT1", l);
        drms_setkey_double(outRec, "CDELT2", 1.0/sinbdivs);
        drms_setkey_string(outRec, "CTYPE1", "Carrington Time");
        drms_setkey_string(outRec, "CTYPE2", "Sine Latitude");
        drms_setkey_string(outRec, "BUNIT", "GAUSS");
        drms_setkey_string(outRec, "DTYPE", "B");
        drms_setkey_string(outRec, "TELESCOP", "SOHO");
        drms_setkey_string(outRec, "INSTRUME", "MDI");
        drms_setkey_string(outRec, "ORIGIN", "SOHO-SOI");
        drms_setkey_string(outRec, "SITE", "md");
        drms_setkey_string(outRec, "DNAME0", "synoptic chart");
        drms_setkey_int(outRec, "CAR_ROT", cr);
        drms_setkey_double(outRec, "CARRTIME", carrtime);

        drms_setkey_time(outRec, "T_OBS", trot);
        drms_setkey_time(outRec, "T_START", tstart);
        drms_setkey_time(outRec, "T_STOP", tstop);
        drms_setkey_time(outRec, "T_REC", trot);
        drms_setkey_time(outRec, "T_ROT", trot);
        drms_setkey_time(outRec, "T_EARTH", tearth);

//        sprint_at(tstr, trot);
//        drms_setkey_string(outRec, "T_OBS", tstr);
//        drms_setkey_time(outRec, "T_OBS", trot);
//        sprint_at(tstr, tstart);
//        drms_setkey_string(outRec, "T_START", tstr);
//        sprint_at(tstr, tstop);
//        drms_setkey_string(outRec, "T_STOP", tstr);
//        drms_setkey_string(outRec, "T_REC", tstr);
//        drms_setkey_string(outRec, "T_ROT", tstr);
//        sprint_at(tstr, tearth);
//        drms_setkey_string(outRec, "T_EARTH", tstr);

//        soho_ephemeris(trot, &r, &b, &l, &vr, &vn, &vw, &tmod); 
        drms_setkey_double(outRec, "B0_ROT", b);
//        soho_ephemeris(tstart, &r, &b, &l, &vr, &vn, &vw, &tmod); 
        drms_setkey_double(outRec, "B0_FRST", b);
//        soho_ephemeris(tstop, &r, &b, &l, &vr, &vn, &vw, &tmod); 
        drms_setkey_double(outRec, "B0_LAST", b);
        drms_setkey_double(outRec, "EARTH_B0", bearth);
        l=(cr-1)*360.0;
        drms_setkey_double(outRec, "LON_FRST", l);
        l=cr*360.0-360.0/len[0];
        drms_setkey_double(outRec, "LON_LAST", l);
        l=-360.0/len[0];
        drms_setkey_double(outRec, "LON_STEP", l);
        l=center;
        drms_setkey_double(outRec, "W_OFFSET", l);
        drms_setkey_string(outRec, "W_WEIGHT", "Even");
        i=lastidx-firstidx+1;
        drms_setkey_int(outRec, "MAG_NUM", i);
        sprint_at(tstr, imrec[firstidx].tobs);
        drms_setkey_string(outRec, "MAG_FRST", tstr);
        sprint_at(tstr, imrec[lastidx].tobs);
        drms_setkey_string(outRec, "MAG_LAST", tstr);
        sprint_at(tstr, magtime);
        drms_setkey_string(outRec, "MAG_ROT", tstr);

        drms_setkey_string(outRec, "SYN_VER", soi_version);
        sprint_ut(tstr, CURRENT_SYSTEM_TIME);
        drms_setkey_string(outRec, "RUNTIME", tstr);
        sscanf(tstr, "%d.%d.%d", &y, &m, &d);
        sprintf(tstr, "%04d-%02d-%02d", y, m, d);
        drms_setkey_string(outRec, "DATE", tstr);
        drms_setkey_float(outRec, HWNWIDTH, halfWindow);
        drms_setkey_float(outRec, EQPOINTS, nEquivPtsReq);

// save the differential rotation correction parameters 
        if (carrStretch > 0)
          {
            drms_setkey_int(outRec, kCARSTRCH, carrStretch);
            drms_setkey_float(outRec, kDIFROT_A, diffrotA);
            drms_setkey_float(outRec, kDIFROT_B, diffrotB);
            drms_setkey_float(outRec, kDIFROT_C, diffrotC);
          }
*/

    if (sortedMagCol)
      {
        for (i = 0; i < len[0]; i++)
          {
	    if (sortedMagCol[i])
	      {
	        free(sortedMagCol[i]);
	      }
          }

         free(sortedMagCol);
       }
  drms_close_records(inRD, DRMS_FREE_RECORD);
  drms_close_records(outRD, DRMS_INSERT_RECORD);
  return 0;
}

int CalcSynCols(int start,
                int end,
                int incr,
                MagCol_t **smc,
                float *synop,
                int *wt,
                char *ww,
                float *epts,
                float *losSynop,
                int *len,
                float center,
                float nEquivPtsReq,
                int los,
                int radialFound,
                float sensAdj,
                float noiseS,
                float nsig,
                float maxNoiseAdj,
                int minOutPts,
                int dlog)
{
   float midSynRow;
   float cosrho; /* equal cos(magnetogram latitude) * cos(magnetogram dlatitude) 
                  *   where dlatitude is positive delta between CM and point.
                  * approximate dlatitude with 'center' (assume for all 
                  * magnetograms, data comes from CM - center, even though
                  *   they will come from points surrounding this column.
                  */
   float cosCenter;
   char *rejDescFormat;
   int col;
   int row;
   int i;
   int j;
   static int nrej = 0;
   int marked;
   float nEquivPts;

   midSynRow = (float)(len[1] - 1) / 2.0;
   cosCenter = cos(center * M_PI / 180);
   rejDescFormat = "      magnetogram: ds=%d, sn=%d, magcol=%d, magrow=%d\n";
   if (minOutPts < 2) minOutPts = 2;

   /* loop over synoptic chart column */
   for (col = start; incr > 0 ? col <= end : col >= end; col += incr)
   {
      float *val = (float *)malloc(sizeof(float) * wt[col]);
      float *ept = (float *)malloc(sizeof(float) * wt[col]);
      MagCol_t *pcols = NULL;
      marked = 0;

      if (dlog)
      {
         pcols = (MagCol_t *)malloc(sizeof(MagCol_t) * wt[col]);
         memset(pcols, 0, sizeof(MagCol_t) * wt[col]);
      }

      if (!val || !ept) printf("MALLOC FAILED");

      /* sort MagCol_ts */
      SortMagCols(smc[col], wt[col]);

      /* print out all mags for this column */

/*
      nEquivPts = 0.0;
      for (i = 0; i < wt[col]; i++) 
      {
         nEquivPts += smc[col][i].equivPts;

         if (i == 0) printf("it's begin\n");
         if (i >= 0 && i < 20) printf("col=%d, i=%d, ds=%d, dist=%f\n", col, i, smc[col][i].ds, smc[col][i].dist);
         if (i == wt[col]-1) printf("it's end\n");

         if (nEquivPts >= nEquivPtsReq && !marked) marked = 1;
      }
*/

      /* loop over synoptic chart row */
      for (row = 0; row < len[1]; ++row)
      {
         float sinLat = (row - midSynRow) / midSynRow;
         float cosLat = sqrt(1 - sinLat * sinLat);
         float sum;
         int npts;
         float sumfinal = 0.0;
         int nptsfinal = 0;
         float eptfinal = 0.0;
         int Maxnepts = nEquivPtsReq + 10;

         if (radialFound) cosrho = cosLat * cosCenter;
         sum = 0.0;
         npts = 0;
         j = 0;
         nEquivPts = 0.0;

         /* Loop over synoptic chart stack above col, row, averaging the values. */
         /* Take the magnetograms with the best values first (smallest dist to CM) 
          * and take only the fewest possible until nEquivPtsReq one-minute-equivalent 
          * points are attained. */
         for (i = 0; i < wt[col] && nEquivPts < Maxnepts; i++)
         {
            float magVal = *(smc[col][i].datacol + row);
            if (!drms_ismissing_float(magVal))
            {
               nEquivPts += smc[col][i].equivPts;
               sum += magVal;
               val[j] = magVal;
               ept[j] = smc[col][i].equivPts;

               if (dlog && pcols) pcols[j] = smc[col][i];
               j++;
               npts++;
            }
         }

         /* Calculate summary statistics */
         if (npts > 1 && nsig > 0.0) {
            double avg = 0.0;
            double ssqr = 0.0;
            double sig;
            float med;
            int nStatPts = 0;
            float *statVals;
            float noiseLevel;
            float *outlier = NULL;
            int statsDone = 0;
            int sIdx = 0;
            float dev = 0.0;
            int nPtsB4Rej = npts;

            /* calculate threshold - below this, point values are within noise levels -
               we don't want to reject them under any circumstances */
            noiseLevel = kNOISE_EQ * noiseS * sensAdj; // noiseLevel not be used so far.
            if (radialFound) noiseLevel = noiseLevel * MIN(1 / cosrho, maxNoiseAdj);
            /* copy the vals - will be destructively analyzing them */
            statVals = (float *)malloc(sizeof(float) * npts);
            if (!statVals) printf("MALLOC_FAILED");

            memcpy(statVals, val, sizeof(float) * npts);
            nStatPts = magStats(&statVals, npts, sum, minOutPts, &avg, &med);
            /* Reject outliers whose values exceeds the noise threshold */
            for (j = 0; j < nPtsB4Rej; ++j)
            {
                  if (!statsDone)
                  {
                     for (sIdx = 0; sIdx < nStatPts; sIdx++)
                     {
                        ssqr += statVals[sIdx] * statVals[sIdx];
                     }

                     sig = sqrt((ssqr - nStatPts * avg * avg) / (nStatPts - 1));
                     statsDone = 1;
                  }

                  dev = fabs(val[j] - med);
                  if (nptsfinal >= nEquivPtsReq) break;
                  if (dev < nsig * sig)
                  {
                     sumfinal += val[j];
                     nptsfinal += 1;
                     eptfinal += ept[j];
                     --npts;
                     ++nrej;
                     if (dlog)
                     {
                        char rejDescBuf[128];

                        snprintf(rejDescBuf,
                                 sizeof(rejDescBuf),
                                 rejDescFormat,
                                 pcols[j].ds,
                                 pcols[j].col,
                                 row);

                     }

                  }
            }

            free(statVals);
         } /* npts > 1 && nsig > 0.0 */

         /* Calcuate the average value for each x,y in the stack */
//       if (npts)
         if (nptsfinal)
         {

            float synVal = sumfinal/nptsfinal;
//            float minVal = (SHRT_MIN + 1);  // * kOutScale;
//            float maxVal = (SHRT_MAX - 1);  // * kOutScale;

               synop[row * len[0] + col] = synVal;
         }
         else
         {
            synop[row * len[0] + col] = DRMS_MISSING_FLOAT;
         }

//         ww[row * len[0] + col] = npts;
//       epts[row * len[0] + col] = nEquivPts;
         ww[row * len[0] + col] = nptsfinal;
         epts[row * len[0] + col] = eptfinal;
         if (los && radialFound)
         {
            int offset = row * len[0] + col;
            if (drms_ismissing_float(synop[offset]))
            {
               losSynop[offset] = DRMS_MISSING_FLOAT;
            }
            else
            {
               losSynop[offset] = synop[offset] * cosLat;
            }
         }
      } /* row loop */

      free(val);
      val = NULL;
      free(ept);
      ept = NULL;

      if (pcols)
      {
         free(pcols);
         pcols = NULL;
      }
   } /* col loop */
   return 0;
}

int cmp(const void *a, const void *b)
{
    float aa = *(float *)a;
    float bb = *(float *)b;
    return (aa<bb) ? -1 : (aa==bb) ? 0 : 1;
}

/* To simplify, assumes one outlier only */
int magStats(float **val, int npts, double sum, int outThreshold, double *avg, float *med)
{
   int actPts = npts;
   float *out = NULL;
   float *in = *val;
   int idx = 0;
   int iOut = 0;
   float medianInt = 0.0;
   float first;
   float last;

   qsort((void *)in, npts, sizeof(float), cmp);

   if (npts > outThreshold && npts > 1)
   {
      actPts = npts - 1;
      out = malloc(sizeof(float) * actPts);

      medianInt = in[npts/2];

      first = in[0];
      last = in[npts - 1];

      if (fabs(first - medianInt) <= fabs(last - medianInt))
      {
	 idx = 0;
	 sum -= last;
      }
      else
      {
	 idx = 1;
	 sum -= first;
      }

      for (iOut = 0; idx < npts && iOut < actPts; idx++, iOut++)
      {
	 out[iOut] = in[idx];
      }

      free(in);
      *val = out;
      in = *val;
   }

   *avg = sum / actPts;
   *med = in[actPts / 2];
   
   return actPts;
}

int CmpMagCol(const void *a, const void *b)
{
   MagCol_t *aa = (MagCol_t *)a;
   MagCol_t *bb = (MagCol_t *)b;

   if (aa->valid && bb->valid)
   {
      return (fabsf(aa->dist) < fabsf(bb->dist)) ? -1 : 
	((fabsf(aa->dist) == fabsf(bb->dist)) ? 0 : 1);
   }
   else if (aa->valid)
   {
      return -1;
   }
   else if (bb->valid)
   {
      return 1;
   }
   else
   {
      return 0;
   }
}

void SortMagCols(MagCol_t *mc, int nelem)
{
   /* qsort((void *)mc, nelem, sizeof(MagCol_t), CmpMagCol); */
   
   /* Try merge sort, since mc will be essentially two sorted lists. */
   int topNeg = -1;
   int topPos = -1;
   MagCol_t *res = NULL;
   int idx;

   for (idx = 0; idx < nelem; idx++)
   {
      if (mc[idx].dist >= 0)
      {
	 topPos = idx;
	 break;
      }
   }

   if (idx == nelem)
   {
      /* only negative dists */
      topNeg = nelem - 1;
   }
   else if (topPos > 0)
   {
      /* both negative and positive dists */
      topNeg = topPos - 1;
   }

   if (topNeg >= 0 || topPos >= 0)
   {
      res = (MagCol_t *)malloc(sizeof(MagCol_t) * nelem);
      memset(res, 0, sizeof(MagCol_t) * nelem);
   }

   if (res)
   {
      idx = 0;
      while (topNeg >=0 && topPos >= 0 && topPos < nelem)
      {
	 if (fabsf(mc[topNeg].dist) < fabsf(mc[topPos].dist))
	 {
	    res[idx] = mc[topNeg];
	    topNeg--;
	 }
	 else
	 {
	    res[idx] = mc[topPos];
	    topPos++;
	 }

	 idx++;
      }

      while (topNeg >= 0)
      {
	 res[idx] = mc[topNeg];
	 topNeg--;
	 idx++;
      }

      while (topPos >= 0 && topPos < nelem)
      {
	 res[idx] = mc[topPos];
	 topPos++;
	 idx++;
      }

      memcpy(mc, res, sizeof(MagCol_t) * nelem);
   }
}

void FreeMagColsData(int start, int end, int incr, int *n, MagCol_t **mc)
{
   if (mc)
   {
      float *data = NULL;
      int cIdx;
      int idx;

      for (cIdx = start; incr > 0 ? cIdx <= end : cIdx >= end; cIdx += incr)
      {
	 for (idx = 0; idx < n[cIdx]; idx++)
	 {
	    data = mc[cIdx][idx].datacol;
	    if (data)
	    {
	       free(data);
	       mc[cIdx][idx].datacol = NULL;
	       mc[cIdx][idx].valid = 0;
	    }
	 }
      }
   }
}

TIME CarringtonTime(int crot, double L)
  {
  double eph[30];
  double CT, clong=0.0;
  double err, CTp50m, CTm50m;
  TIME t;
  char *stime;
  char tstr[64];
  CT = 0.0;
  static TIME T1853 = 0.0, delta_T = 0.0;
  T1853 = sscan_time("1853.11.09_22:27:24_UTC");

  delta_T = sscan_time("1977.01.01_00:00:00_TAI") - sscan_time("1601.01.01_00:00:00_UT"); 

  CT = 360.0 * crot - L;
  t = (C2 + CT * C1) * SID;

//  t = T1853 + CT * C1 * SID;
  solephem(t, eph);
  err = CT - eph[8];
  solephem(t+6*3600.0, eph);
  CTp50m = eph[8];
  solephem(t-6*3600.0, eph);
  CTm50m = eph[8];
        // interpolate to correct t
  t += 12*3600 * (err/(CTp50m - CTm50m));
  solephem(t, eph);
  err = CT - eph[8];
  solephem(t+300.0, eph);
  CTp50m = eph[8];
  solephem(t-300.0, eph);
  CTm50m = eph[8];
        // interpolate to correct t
  t += 600 * (err/(CTp50m - CTm50m));

/*
  sprint_at(tstr, t-delta_T);
printf("TSTART=%s\n", tstr);

solephem(t, eph); 
CT = eph[8];
crot = floor(CT/360.0);
clong = 360 - (CT - 360*crot);
crot += 1;
printf("CT%04d:%06.2f\n", crot, clong);
        printf("CT=%10.3f, B0=%f\n", CT, eph[9]/RADSINDEG);
*/
  return(t);
  }

/*
double local_earth_meridian_crossing(double L, int cr)
{
    double ta, tb, t, la, lb, l;
    double ephem[EPHEM_SIZE];
    double trot;
    
    trot = meridian_crossing(L, cr);
    ta = trot - 1800.0;
    tb = trot + 1800.0;
    sun_ephemeris(ta, ephem, 0.0, 0.0);
    la = fmod(ephem[EPHEM_L0],360.0)-L;
    sun_ephemeris(tb, ephem, 0.0, 0.0);
    lb = fmod(ephem[EPHEM_L0],360.0)-L;
    while (fabs(tb-ta) > 0.01) {
	t = ta - (tb - ta) * la / (lb - la);
	sun_ephemeris(t, ephem, 0.0, 0.0);
	l = fmod(ephem[EPHEM_L0],360.0)-L;
	ta = tb;
	la = lb;
	tb = t;
	lb = l;
    }

    return t;
}

double earth_B(TIME t)
{
    double ephem[EPHEM_SIZE];
    sun_ephemeris(t, ephem, 0.0, 0.0);
    return ephem[EPHEM_B0]*180.0/M_PI;
}

*/
/*
$Id: mdisynop.c,v 1.1 2013/07/08 17:36:04 yliu Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/mdisynop.c,v $
$Author: yliu $
*/
/* $Log: mdisynop.c,v $
 * Revision 1.1  2013/07/08 17:36:04  yliu
 * new module for producing synoptic charts for MDI line-of-sight mags
 *
 * Revision 1.24  2007/10/26 17:51:39  arta
 * Fix bug where for loop limit was changed within loop.
 *
 * Revision 1.23  2007/10/18 16:45:58  arta
 * Fix bug where there was an attempt to read pixels just off the left or right edge of a remapped image when center was either 60 or -60.
 *
 * Revision 1.22  2007/09/27 17:34:25  arta
 * Change name of degrees per halfwindow parameter to HWNWIDTH.
 *
 * Revision 1.21  2007/09/27 17:12:28  arta
 * Change the W_HWIDTH attribute type to a float, to match the variable used to populate it.
 *
 * Revision 1.20  2007/09/26 16:06:41  arta
 * Make log more accurately show the stack of mags for each synoptic chart column.
 *
 * Revision 1.19  2007/09/26 00:28:38  arta
 * Use dist from CM, not dist from mid col.
 *
 * Revision 1.18  2007/08/20 16:59:41  arta
 * Add output to help identify data set name.  Fix sensAdj detection.
 *
 * Revision 1.17  2007/08/10 23:58:00  arta
 * Write out carr stretch keywords.
 *
 * Revision 1.16  2007/08/09 15:10:49  arta
 * Do proper check for out-of-bounds array index.
 *
 * Revision 1.15  2007/08/06 18:49:03  arta
 * Fix for crash in sgi4.
 *
 * Revision 1.14  2007/07/18  14:53:27  arta
 * Improve performance.  Use merge-sort, not quicksort, as data were essentially read in into two sorted lists.  Decrease memory footprint - do synop columen the data are ready, don't wait until the end.  Clear memory right after each column is done.
 *
 * Revision 1.13  2007/07/16 21:01:49  arta
 * new outlier detection.  Calculate noise at various cos(rho).
 *
 * Revision 1.12  2007/07/11 16:30:56  arta
 * magsynop - add line-of-sight flag (to output line-of-sight images in addition to radial ones), add force flag (to allow processing of obsolete data), fix bug in radial image detection
 * fixplatescale_v0 - write out soi_vers_num in BLDVER18 attribute
 * fdmagcal and fdradial - set pixels outside of rsun to missing (do not leave them unmodified), do not write BLDVER18 in fdmagcal.
 *
 * Revision 1.11  2007/05/14 18:44:49  arta
 * Ensure v2helio retains the BLDVER18 keyword for magsynop use.
 *
 * Revision 1.10  2007/05/14 16:57:39  arta
 * Use BLDVER18, not BLDVER15 when assessing how current data are.
 *
 * Revision 1.9  2007/05/14 05:10:19  arta
 * Fix sgi build.
 *
 * Revision 1.8  2007/05/08 23:39:18  arta
 * Check for consistent MDI CM build in magsynop.  Requires passing along keywords in v2helio.  Also, ensure that all images are either radial or line-of-sight, but not a mixture.
 *
 * Revision 1.7  2007/05/08 17:23:40  arta
 * Store epts file in a more compact form.
 *
 * Revision 1.6  2007/05/05 03:31:05  arta
 * Cap values larger than can be contained by max value.
 *
 * Revision 1.5  2007/03/26 19:09:49  arta
 * Track magnetogram row, col, dataseries, recordnum info to diagnose missing data.
 *
 * Revision 1.4  2007/03/16 16:35:10  arta
 * Modify algorithm to use the BEST magnetograms for each synoptic column.  BEST is closest to CM.  Use as many magnetograms as required to reach 20 1-minute magnetogram equivalence points.  Crop magnetograms at +/- 20 of CM.
 *
 * Revision 1.3  2007/01/26 18:25:28  rick
 * changed name of earth_merid... to local_earth_merid to avoid lib conflict
 *
 * Revision 1.2  2000/05/30 16:44:14  kehcheng
 * initial revision
 * */
