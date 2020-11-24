#include "astro.h"
#include "jsoc_main.h"
#include "list.h"

#define DEBUG 0


/* NOTE that all "GCI" or "HCI" are really "GAE" or "HAE" coordiantes.
 * The "AE" part of thes names is "Aries Ecliptic" where the X-axis points to Aires,
 * the equinox, and the Z axis is normal to the ecliptic.  The "CI" coordiantes are
 * Heliocentric Inertial where the x-axis is pointed to the ascending node of the solar
 * equator on the ecliptic.  We convert from HAE to HCI by first rotating about the Z axis
 * by Omega, the distance along the ecliptic from the equinox to the ascending node, about
 * 76 degrees.  Then by rotating by 7.25 degrees about the Y axis to a system that has
 * the Z axis along solar rotation axis.  PHS, Oct 2020.
 *
 * Code changes in OCt 2020 correct CRLN_OBS and calculate HGLN_OBS.
 */


#define kMaxRecIdSize 128

#define kCACHEKEYSIZE 64
#define kCACHESIZE 128

#define kPI M_PI
#define kDEG2RAD (kPI / 180)

/* Ecliptic values */
#define kALPHA (75.76 * kDEG2RAD)
#define kDELTA (7.25 * kDEG2RAD)

#define kRSUNREF 696000000.00		// Sun radius in m

/* For SDOCarringtonCoords() */
#define TWO_PI 		(2*kPI)
#define C  		(299792458)      // c in m/s
#define AU		(149597870700)   // AU in m
#define CARR_DEGDAY     (14.1844000)      // Adopted degrees per day, includes precession
// #define J2000           (725803167.816)   // sscan_time("2000.01.01_00") J2000 base time
#define TCARR           (-3881476832.184) // sscan_time("1854.01.01_12:00_TT")   Carr ref epcoh est.
#define CARR_SYNODIC    (CARR_DEGDAY - 360.0/365.25) // estimate of synodic carrington rotation rate

enum IORBIT_Slotpos_enum
{
   kMISSING = 0,
   kLT,
   kGT,
   kEQ,
   kLTE,
   kGTE
};

typedef enum IORBIT_Slotpos_enum IORBIT_Slotpos_t;

struct IORBIT_Vector_struct
{
    long long slot; /* slot number - missing if vector is not for a grid time */
    double obstime;
    double gae_x;
    double gae_y;
    double gae_z;
    double gae_vx;
    double gae_vy;
    double gae_vz;
    double hae_x;
    double hae_y;
    double hae_z;
    double hae_vx;
    double hae_vy;
    double hae_vz;
    double rsunobs;
    double obsvr;
    double dsunobs;
};

typedef struct IORBIT_Vector_struct IORBIT_Vector_t;

/* If not chunking, then gGridCache never gets updated during module run (it has all the values it
 * will ever have right from the beginning).
 */
static int gCacheChunking = 1;
static IORBIT_Vector_t *gGridCache = NULL;
static DRMS_RecordSet_t *gCacheRS = NULL;
static int gGridNItems = 0;
static int gGridCurrPos = -1;
static HContainer_t *gKeyMap = NULL;

/* for !gCacheChuunking case */
const int gMinAbove = 16; /* min no. of vectors above max tgttime */
const int gMinBelow = 16; /* min no. of vectors below min tgttime */


// SDOCarringtonCoords takes time, distance to Sun, B0 in radians, and the heliocentric ecliptic longitude in radians.
// t in TAI and obsdist in km.
// It returns Carrington rotation, latitude, and longitude, as int for rotation and in degrees for angles.
// It corrects for light travel time difference between actual distance and 1 AU, corrects for distance to Sun surface. 
// crot, L, B are Carrington rotation, longitude, and latitude of sub-observer point
// The method has been tested for 2010 to 2020 with negligible differences to JPL or SunPy. 

static void SDOCarringtonCoords(TIME t, double obsdist, double b, double hci_long, int *crot, double *L, double *B)
{
   double solrots;       // solar rotations by time t
   double carr_rots;     // Rotations since Carrington Epoch
   double l;             // longitude working var (rotations)
   double tlight;        // Light travel time from Sun center to observer
   double d_since_TCARR; // days at Sun since TCARR
   double clest;         // estimate of carrington longitude
   double correction = -0.130693;    // Correction mostly from change in CARR_DEGDAY definition in c. 1976
   int rot;              // Carrington rotation number of sub-observer point

   tlight = (obsdist*1000 -  kRSUNREF - AU)/C;
   d_since_TCARR = (t - tlight - TCARR)/86400;
   solrots = (CARR_DEGDAY * d_since_TCARR + correction)*kDEG2RAD;
   l = hci_long - solrots;
   l = modf(l/TWO_PI, &carr_rots);  // carrington longitude in rotations with bias
   if(l < 0) l += 1;

   // Get rotation number, first from "guess" then adjust based on accurate longitude.
   // guess is good to a few degrees.  2 rotations happened before the Carrington ref time
   carr_rots = 3.0 + CARR_SYNODIC * (d_since_TCARR + correction)/360.0;
   rot = (int)carr_rots;
   clest = (1 + rot - carr_rots);
   if ((l - clest) > 0.5)
     rot++;
   if ((clest - l) > 0.5)
     rot--;

   *crot = rot;
   *L = 360*l;
   *B = b / kDEG2RAD;
}

/* function to compare times - used by SortTimes() */
static int CmpTimes(const void *a, const void *b)
{
   double *aa = (double *)a;
   double *bb = (double *)b;

   if (aa && bb)
   {
      if (drms_ismissing_time(*aa) && drms_ismissing_time(*bb))
      {
         return 0;
      }
      else if (drms_ismissing_time(*aa))
      {
         return -1;
      }
      else if (drms_ismissing_time(*bb))
      {
         return 1;
      }
      else
      {
         return (*aa < *bb ? -1 : (*aa == *bb ? 0 : 1));
      }
   }
   else if (aa)
   {
      return 1;
   }
   else if (bb)
   {
      return -1;
   }
   else
   {
      return 0;
   }
}

/* assumes requested times are all larger than the missing time value (or they are
 * missing time values).
*/
static int SortTimes(const double *unsorted, int nitems, double **sorted)
{
    int nret = 0;
    double *sortedint = malloc(nitems * sizeof(double));
    int itime;

    if (sortedint)
    {
        memcpy(sortedint, unsorted, nitems * sizeof(double));
        qsort(sortedint, nitems, sizeof(double), CmpTimes);

        if (sorted)
        {
            itime = 0;
            while (drms_ismissing_time(sortedint[itime]) && itime < nitems)
            {
                itime++;
            }

            /* itime is index of first valid time */
            if (itime < nitems)
            {
                /* at least one valid time */
                nret = nitems - itime;
                *sorted = malloc(nret * sizeof(double));
                memcpy(*sorted, &(sortedint[itime]), nret * sizeof(double));
            }
        }

        free(sortedint);
        sortedint = NULL;
    }
    else
    {
        fprintf(stderr, "SortTime() - out of memory.\n");
    }

    return nret;
}

/* Like SortTimes, but don't strip missing values. */
static int SortTimesB(const double *unsorted, int nitems, double **sorted)
{
    int nret = 0;
    double *sortedint = malloc(nitems * sizeof(double));

    if (sortedint)
    {
        memcpy(sortedint, unsorted, nitems * sizeof(double));
        qsort(sortedint, nitems, sizeof(double), CmpTimes);

        if (sorted)
        {
            nret = nitems;
            *sorted = malloc(nret * sizeof(double));
            memcpy(*sorted, sortedint, nret * sizeof(double));
        }

        free(sortedint);
        sortedint = NULL;
    }
    else
    {
        fprintf(stderr, "SortTime() - out of memory.\n");
    }

    return nret;
}

/* return the relationship of slottime to tgttime */
static IORBIT_Slotpos_t GetSlotPos(double slottime, double tgttime)
{
   IORBIT_Slotpos_t slotpos;

   if (!drms_ismissing_time(slottime))
   {
      if (slottime < tgttime)
      {
         slotpos = kLT;
      }
      else if (slottime > tgttime)
      {
         slotpos = kGT;
      }
      else
      {
         slotpos = kEQ;
      }
   }
   else
   {
      slotpos = kMISSING;
   }

   return slotpos;
}

/* Returns 1 if the posA ==> posB */
static int IsPos(IORBIT_Slotpos_t posA, IORBIT_Slotpos_t posB)
{
   return ((posA == kLT && (posB == kLTE || posB == kLT)) ||
           (posA == kGT && (posB == kGTE || posB == kGT)) ||
           (posA == kEQ && posB == kLTE) ||
           (posA == kEQ && posB == kGTE) ||
           (posA == posB));
}

/* Cache functions
 *
 * record-chunking will actually make this obsolete, but for now, just
 * cache the grid times.
 */

static inline int CreateLHashKey(char *hashkey, int size, long long slot)
{
   return snprintf(hashkey, size, "%lld", slot);
}

static inline int CreateDHashKey(char *hashkey, int size, double val)
{
   return snprintf(hashkey, size, "%f", val);
}

static IORBIT_Vector_t *CreateCache()
{
   if (!gGridCache)
   {
      gGridCache = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * kCACHESIZE);
   }

   return gGridCache;
}

static void DestroyCache()
{
   if (gCacheRS)
   {
      drms_close_records(gCacheRS, DRMS_FREE_RECORD);
      gCacheRS = NULL;
   }

   if (gGridCache)
   {
      free(gGridCache);
      gGridCache = NULL;
   }
}

/* bsearch for a vec grid slot that satisfies slotpos = GetSlotPos(vectime, tgttime); IsPos(slotpos, pos) */
static IORBIT_Vector_t *LookupInCache(double tgttime, IORBIT_Slotpos_t pos)
{
   IORBIT_Vector_t *val = NULL;
   IORBIT_Vector_t *vec = NULL;
   int icache;
   int top;
   int bottom;
   int direction;
   IORBIT_Slotpos_t slotpos;

   if (gGridCache)
   {
      top = gGridNItems - 1;
      bottom = 0;

      if (IsPos(pos, kLTE))
      {
         direction = -1;
      }
      else
      {
         direction = 1;
      }

      if (direction == -1)
      {
         while (top != bottom)
         {
            icache = (int)ceil((double)(top + bottom) / 2);
            vec = &(gGridCache[icache]);
            slotpos = GetSlotPos(vec->obstime, tgttime);
            if (IsPos(slotpos, pos))
            {
               /* Got a vector that satisfies pos, but perhaps not the final one. Don't
                * exclude this vector from the possible final result, but exclude every
                * vector with a smaller obstime than this vector.
                */
               bottom = icache;
            }
            else
            {
               /* This vector doesn't satisfy pos.  It and all vectors with a larger
                * obstime can be excluded.
                */
               top = icache - 1;
            }
         }
      }
      else
      {
         while (top != bottom)
         {
            icache = (int)floor((double)(top + bottom) / 2);
            vec = &(gGridCache[icache]);
            slotpos = GetSlotPos(vec->obstime, tgttime);
            if (IsPos(slotpos, pos))
            {
               /* Got a vector that satisfies pos, but perhaps not the final one. Don't
                * exclude this vector from the possible final result, but exclude every
                * vector with a larger obstime than this vector.
                */
               top = icache;
            }
            else
            {
               /* This vector doesn't satisfy pos.  It and all vectors with a smaller
                * obstime can be excluded.
                */
               bottom = icache + 1;
            }
         }
      }

      /* top == bottom */
      vec = &(gGridCache[top]);

      /* It is possible that gGridCache[top] isn't the correct result - it never gets tested if
       * the bottom consistently failed until bottom == top. */
      slotpos = GetSlotPos(vec->obstime, tgttime);
      if (IsPos(slotpos, pos))
      {
         val = vec;
         gGridCurrPos = top;
      }
      else
      {
         gGridCurrPos = -1;
      }
   }

   return val;
}

static double GetKey(DRMS_Record_t *rec, const char *keyID)
{
    const char *keyname = NULL;

    keyname = (const char *)hcon_lookup(gKeyMap, keyID);
    if (keyname)
    {
        return drms_getkey_double(rec, keyname, NULL);
    }
    else
    {
        return DRMS_MISSING_DOUBLE;
    }
}

static void PopulateVector(DRMS_Record_t *rec, IORBIT_Vector_t *vec)
{
    const char *keyname = NULL;
    char indxkey[DRMS_MAXKEYNAMELEN] = {0};

    if (gKeyMap)
    {
        /* Orbit-series keyword names are custom. Map to the custom names from the keyword IDs. It is possible that
         * some of these keywords do not exist in the current orbit series. If so, then the value assigned will
         * be missing. */
        vec->obstime = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_DATE));
        vec->gae_x = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_X));
        vec->gae_y = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_Y));
        vec->gae_z = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_Z));
        vec->gae_vx = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_VX));
        vec->gae_vy = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_VY));
        vec->gae_vz = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_VZ));
        vec->hae_x = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_X));
        vec->hae_y = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_Y));
        vec->hae_z = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_Z));
        vec->hae_vx = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_VX));
        vec->hae_vy = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_VY));
        vec->hae_vz = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_VZ));
        vec->rsunobs = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_RSUN_OBS));
        vec->obsvr = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_VR));
        vec->dsunobs = GetKey(rec, IORBIT_STRINGIFY(IORBIT_KEYWORD_DSUN_OBS));

        keyname = (const char *)hcon_lookup(gKeyMap, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_DATE));
        if (keyname)
        {
            snprintf(indxkey, sizeof(indxkey), "%s_index", keyname);
        }
    }
    else
    {
        /* Orbit-series keyword names are standard. Use the defines to locate keyword names. */
        vec->obstime = drms_getkey_double(rec, IORBIT_KEYWORD_OBS_DATE, NULL);
        vec->gae_x = drms_getkey_double(rec, IORBIT_KEYWORD_GAE_X, NULL);
        vec->gae_y = drms_getkey_double(rec, IORBIT_KEYWORD_GAE_Y, NULL);
        vec->gae_z = drms_getkey_double(rec, IORBIT_KEYWORD_GAE_Z, NULL);
        vec->gae_vx = drms_getkey_double(rec, IORBIT_KEYWORD_GAE_VX, NULL);
        vec->gae_vy = drms_getkey_double(rec, IORBIT_KEYWORD_GAE_VY, NULL);
        vec->gae_vz = drms_getkey_double(rec, IORBIT_KEYWORD_GAE_VZ, NULL);
        vec->hae_x = drms_getkey_double(rec, IORBIT_KEYWORD_HAE_X, NULL);
        vec->hae_y = drms_getkey_double(rec, IORBIT_KEYWORD_HAE_Y, NULL);
        vec->hae_z = drms_getkey_double(rec, IORBIT_KEYWORD_HAE_Z, NULL);
        vec->hae_vx = drms_getkey_double(rec, IORBIT_KEYWORD_HAE_VX, NULL);
        vec->hae_vy = drms_getkey_double(rec, IORBIT_KEYWORD_HAE_VY, NULL);
        vec->hae_vz = drms_getkey_double(rec, IORBIT_KEYWORD_HAE_VZ, NULL);

        /* These keywords do not exist in the standard orbit series. */
        vec->rsunobs = DRMS_MISSING_DOUBLE;
        vec->obsvr = DRMS_MISSING_DOUBLE;
        vec->dsunobs = DRMS_MISSING_DOUBLE;

        snprintf(indxkey, sizeof(indxkey), "%s_index", IORBIT_KEYWORD_OBS_DATE);
    }

    vec->slot = drms_getkey_longlong(rec, indxkey, NULL);
}

/*
 * srcseries - datseries containing orbit-data information (like sdo.fds_orbit_vectors)
 * optfilter - a record-set query that causes a subset of all orbit data to be used; if NULL
 *    all orbit data used.
 * nkeep - if the cache doesn't contain a suitable grid vector, the cache must be
 *   rehydrated. nkeep is the number of existing cache items to retain when rehydrating.
 *   So the old cache isn't completed blown away. This is done so that FetchNext()
 *   and FetchPrevious() can find the next and previous items without causing a cache miss.
 *
 * returns number of items successfully cached */
static int RehydrateCache(DRMS_Env_t *env,
                          const char *srcseries,
                          const char *optfilter,
                          int nkeep)
{
   int nitems = 0;
   int ret = 0;

   IORBIT_Vector_t *vec;
   DRMS_Record_t *rec = NULL; /* record of orbit data from FDS */
   int iitem = 0;
   char query[256];
   int drmsstatus;
   DRMS_RecChunking_t cstat = kRecChunking_None;
   int newchunk;
   int firstitem = 0;
   int lastitem = 0;

   snprintf(query,
            sizeof(query),
            "%s%s",
            srcseries, optfilter ? optfilter : "");

   if (!gCacheChunking)
   {
      if (gCacheRS)
      {
         /* Cannot truly RE-hydrate if caching is turned off */
         fprintf(stderr, "Cannot rehydrate cache if caching is turned off.\n");
      }
      else
      {
         /* OK to populate cache one time when caching is turned off */
         gCacheRS = drms_open_records(env, query, &drmsstatus);

         if (gCacheRS && gCacheRS->n > 0)
         {
            gGridCache = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * gCacheRS->n);
            firstitem = 0;
            lastitem = gCacheRS->n;

            for (iitem = firstitem; iitem < lastitem; iitem++)
            {
               rec = gCacheRS->records[iitem];
               vec = &(gGridCache[iitem]);

               /* These are in J2000.0, ecliptic coordinates
                *   -gcipos and gcivel are values relative to earth
                *   -hcipos and hcivel are values relative to the sun
                */
                PopulateVector(rec, vec);
            }
         }
      }

      nitems = lastitem - firstitem;
      gGridNItems = nitems;
      ret = nitems;
   }
   else
   {
      /* gCacheChunking == 1 */
      int newcache = 0;

      if (!gGridCache)
      {
         gGridCache = CreateCache();
         newcache = 1;
      }

      if (!gCacheRS)
      {
         /* Use cursor to retrieve just a chunk */
         drms_recordset_setchunksize(kCACHESIZE);
         gCacheRS = drms_open_recordset(env, query, &drmsstatus);
      }

      /* may be fewer items than this, but this will handled in the loop below. */
      if (newcache)
      {
         nkeep = 0;
      }

      lastitem = kCACHESIZE - 1;
      firstitem = nkeep;

      /* must keep the last nkeep items since higher-level calls might
       * want to move backwards in cache up to nkeep items */
      for (iitem = 0; iitem < nkeep; iitem++)
      {
         gGridCache[iitem] = gGridCache[gGridNItems - nkeep + iitem];
      }

      /* flush the rest of the cache */
      memset(gGridCache + nkeep, 0, sizeof(IORBIT_Vector_t) * (kCACHESIZE - nkeep));

      for (iitem = firstitem; iitem <= lastitem; iitem++)
      {
         rec = drms_recordset_fetchnext(env, gCacheRS, &drmsstatus, &cstat, &newchunk);
         if (!rec)
         {
            /* last time, read last record in entire set. */
            break;
         }

         vec = &(gGridCache[iitem]);

         /* These are in J2000.0, ecliptic coordinates
          *   -gcipos and gcivel are values relative to earth
          *   -hcipos and hcivel are values relative to the sun
          */
          PopulateVector(rec, vec);
      }

      nitems = iitem - firstitem + nkeep;
      gGridNItems = nitems;
      ret = nitems - nkeep;
   }

   return ret;
}

/*
 * srcseries - datseries containing orbit-data information (like sdo.fds_orbit_vectors)
 * optfilter - a record-set query that causes a subset of all orbit data to be used; if NULL
 *    all orbit data used.
 * tgttime - the interpolated time for which we'd like find a nearby grid vector
 * nkeep - if the cache doesn't contain a suitable grid vector, the cache must be
 *   rehydrated. nkeep is the number of existing cache items to retain when rehydrating.
 *   So the old cache isn't completed blown away. This is done so that FetchNext()
 *   and FetchPrevious() can find the next and previous items without causing a cache miss.
 * pos - this parameter indicates the relationship between the tgttime and the grid vector
 *   obstime. If pos == kLT, then the Fetch() request looks for a grid vector with an
 *   obstime that is LT the tgttime.
 *
 * Fetch() assumes that there will be a vec-grid very close to tgttime, which should be
 * true, unless a user is requesting a tgttime that lies outside the existing
 * vector grid.
 *
 * Fetch() always searches the orbit series in increasing OBS_DATE.
 */
static IORBIT_Vector_t *Fetch(DRMS_Env_t *env,
                              const char *srcseries,
                              const char *optfilter,
                              double tgttime,
                              int nkeep,
                              IORBIT_Slotpos_t pos)
{
   IORBIT_Vector_t *vec = NULL;
   int again;

   /* THIS IS TRICKY! It is possible to find a vector that is kLTE than tgttime, but
    * is not the immediately-kLTE one. LookupInCache() could have selected the vector
    * that had the largest obstime in the CACHED vectors, but in fact other recordset
    * chunks had larger obstimes that were smaller than tgttime.
    *
    * The situation is asymmetrical.  If pos == kGT, we won't have vectors in the cache
    * that are GT the tgttime while there exist recordset chunks that have vectors
    * whose obstimes are less than the ones in the cache, but greater than the tgttime.
    * This won't be the case because we make Fetch() calls in increases tgttime order.
    */
   again = 1;

   while (again)
   {
      vec = LookupInCache(tgttime, pos);

      if (vec)
      {
         if (!(gGridCurrPos == gGridNItems - 1 && (pos == kLT || pos == kLTE)))
         {
            /* We have a grid vector, from the cache, that is LT the target time.
             * And grid vector does not have the largest obstime in the cache.
             * If the latter were not true, and vec had the largest obstime of
             * any vector in the cache, then this leaves open the possibility
             * that there are recordset chunks that contain grid vectors whose
             * obstimes are larger than vec's, but still LT the target time.
             */
            again = 0;
         }
         else if (!gCacheChunking)
         {
            /* The vec selected is the one with the largest obs time in the cache.
             * So it may be the case that the optfilter didn't correctly select
             * a sufficient number of records with obstimes above the tgttime.
             */
            fprintf(stderr,
                    "WARNING: There might be an insufficient number of records in '%s%s' to span the range needed for interpolation.\n",
                    srcseries,
                    optfilter);
            again = 0;
         }
      }
      else
      {
         again = 1;
      }

      if (again)
      {
         /* nkeep should be sufficiently large to accommodate the code calling
          * Fetch(), which may call FetchNext() or  FetchPrevious(). */
         TIMER_t *tmr = NULL;

         if (env->verbose)
         {
            /* time fetching records from disk */
            tmr = CreateTimer();
         }

         if (RehydrateCache(env, srcseries, optfilter, nkeep) == 0)
         {
            /* No more data - either didn't find a grid vector with an obstime
             * suitable for pos and tgttime, or it was
             * the very last item in the entire series. If the latter is true
             * then vec will now contain that last item, which is the correct
             * result. */
            again = 0;
         }

         if (env->verbose)
         {
            fprintf(stdout, "RehydrateCache() seconds elapsed: %f\n", GetElapsedTime(tmr));
            DestroyTimer(&tmr);
         }
      }
   }

   return vec;
}

static IORBIT_Vector_t *FetchPrevious()
{
   IORBIT_Vector_t *vec = NULL;

   /* Can't got backward in gridvec array */
   if (gGridCurrPos > 0)
   {
      vec = &(gGridCache[--gGridCurrPos]);
   }

   return vec;
}

static IORBIT_Vector_t *FetchNext(DRMS_Env_t *env,
                                  const char *srcseries,
                                  const char *optfilter,
                                  int nkeep)
{
   IORBIT_Vector_t *vec = NULL;

   if (gGridCurrPos < gGridNItems - 1)
   {
      vec = &(gGridCache[++gGridCurrPos]);
   }
   else
   {
      /* even if optfilter != NULL, it is okay to fetch forward in the gridvec array */
      vec = Fetch(env, srcseries, optfilter, gGridCache[gGridNItems - 1].obstime, nkeep, kGT);
   }

   return vec;
}

/* If not caching grid vectors between calls to iorbit_getinfo(), then all needed vectors
 * must be in gGridCache. The "needed" vectors include a buffer of gMinAbove and gMinBelow
 * above and below the min and max obstimes. In the top-level function, 1 hour above and
 * below is added to the min and max obstimes, so if there are at least 16 observations
 * within each hour, this function will succeed. */
static int CheckCache(double mintgt, double maxtgt)
{
   int ok = 1;

   if (!gCacheChunking)
   {
      /* If not caching, ensure that there are a sufficient number of vectors above and
       * below tgttime */
      int itime;
      itime = gGridNItems - 1;

      while (gGridNItems - 1 - itime < gMinAbove)
      {
         if (gGridCache[itime].obstime <= maxtgt)
         {
            ok = 0;
            break;
         }

         itime--;
      }

      itime = 0;
      while (itime < gMinBelow)
      {
         if (gGridCache[itime].obstime >= mintgt)
         {
            ok = 0;
            break;
         }

         itime++;
      }
   }

   return ok;
}

/*
 * indices - For each tgttime, returns the index into gvectors of vector whose
 * obstime is the largest time that is less than or equal to the tgttime.
 *
 * returns the grid vectors needed in order to interpolate values for tgttimes.
 * posbelow indicates how many grid vectors below each tgttime are needed
 * in order to interpolate to the tgttime. posabove indicates how many grid vectors
 * above each tgttime are needed for the interpolation.
 */
static int GetGridVectors(DRMS_Env_t *env,
                          const char *srcseries,
                          double *tgttimes,
                          int ntgts,
                          int nbelow,
                          IORBIT_Slotpos_t posbelow,
                          int nabove,
                          IORBIT_Slotpos_t posabove,
                          LinkedList_t **gvectors,
                          int *ngvecs,
                          HContainer_t **indices,
                          const char *optfilter)
{
   int err = 0;
   int actpoints;
   int totpoints = 0;
   IORBIT_Vector_t *vec = NULL;
   IORBIT_Slotpos_t slotpos;

   /* For each tgttime's interpolation, vecsbelow holds the grid vectors with obstimes
    * LT the tgttime and vecsabove holds the grid vectors with obstimes GT the tgttime
    *  necessary to do the interpolation*/
   IORBIT_Vector_t *vecsbelow = NULL;
   IORBIT_Vector_t *vecsabove = NULL;
   int itgt;
   int ibelow;
   int iabove;
   char lhashkey[128];
   char dhashkey[128];

   /* key is the grid vector slot (one grid vector per series slot)
    * and value is the index into gvectors. */
   HContainer_t *slottovecindex = NULL;

   if (gvectors)
   {
      /* The grid of FDS orbit vectors is NOT regular (because of leap seconds). Usually,
       * the times are 60 seconds apart, but sometimes they are 61 seconds apart. So, you can't
       * assume that observation time of the data in a slotted key's slot is the time
       * at the beginning of that slot.
       *
       * So, make the grid-time vector cache keyed by slot number. Then
       * map the tgttime to a slot, and check the cache for that slot. If it exists, then
       * check the obsdate in the slot (it could fall anywhere within the slot because of
       * the leap-seconds complication). If the obsdate is less than the tgttime, then
       * you have the closest value below the tgttime. If the obsdate is greater than the
       * tgttime, then you have the closest value above the tgttime. In this manner,
       * you can get as many grid values as needed above and below the tgttime.
       */
      *gvectors = list_llcreate(sizeof(IORBIT_Vector_t), NULL);
      vecsbelow = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * nbelow);
      vecsabove = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * nabove);
      slottovecindex = hcon_create(sizeof(int), 128, NULL, NULL, NULL, NULL, 0);

      if (indices)
      {
         *indices = hcon_create(sizeof(int), 128, NULL, NULL, NULL, NULL, 0);
      }

      for (itgt = 0; itgt < ntgts; itgt++)
      {
         /* Cache. If the cache is already full, and tgttimes[itgt] is a miss
          * (there is no grid vector in the cache whose obs time is less than tgttimes[itgt])
          * then  keep the last (nbelow + 4) cached vectors currently in the cache. Blow away
          * the rest of the cache, and rehydrate with the next chunk of vectors.
          * Keep these old vectors because the next while loop is going to walk backward in
          * the cache trying to find nbelow vectors that are LTE the tgttimes[itgt] time.
          * 4 is a buffer. Also, Fetch() and FetchNext() below could have caused
          * rehydration of a newer recordset chunk, causing this Fetch() to fail if
          * we didn't keep some of the previous chunk cached. */
         vec = Fetch(env, srcseries, optfilter, tgttimes[itgt], nbelow + 4, posbelow);

         if (!vec)
         {
            fprintf(stderr,
                    "There is no grid vector smaller than %f in '%s'.\n",
                    tgttimes[itgt],
                    srcseries);
            err = 1;
            break;
         }

         if (itgt == 0 && !gCacheChunking)
         {
            if (!CheckCache(tgttimes[0], tgttimes[ntgts - 1]))
            {
               fprintf(stderr, "Not enough data records in '%s%s' to perform interpolation.\n", srcseries, optfilter);
               err = 1;
               break;
            }
         }

         /* get the nbelow points less than tgttimes[itgt] (in reverse order of course) */
         actpoints = 0;
         while (1)
         {
            /* By definition, Fetch() returns the grid vector with an obstime that returns
             * true for IsPos(slotpos, posbelow).  It returns the vector with the largest
             * obstime that satisfies this. */
            slotpos = GetSlotPos(vec->obstime, tgttimes[itgt]);

            if (IsPos(slotpos, posbelow))
            {
               vecsbelow[nbelow - actpoints - 1] = *vec;
               actpoints++;
            }

            if (actpoints == nbelow)
            {
               break;
            }

            /* need to fetch previous grid vector */
            vec = FetchPrevious();
            if (!vec)
            {
               fprintf(stderr,
                       "Series '%s' missing a grid vector smaller than %f.\n",
                       srcseries,
                       tgttimes[itgt]);

               err = 1;
               break;
            }
         }

         if (err)
         {
            break;
         }

         /* vec now points to the vector that is nbelow points below the tgttime. Insert
          * that vector into list of vectors to be returned if it is not already there.
          * Use the fact that the list of vectors must be monotonically according to
          * obstime.
          */

         /* okay to insert the nbelow tgttimes */
         for (ibelow = 0; ibelow < nbelow; ibelow++)
         {
            CreateLHashKey(lhashkey, sizeof(lhashkey), vecsbelow[ibelow].slot);

            if (hcon_lookup(slottovecindex, lhashkey))
            {
               continue;
            }

            list_llinserttail(*gvectors, &(vecsbelow[ibelow]));
            hcon_insert(slottovecindex, lhashkey, &totpoints); /* map vec->slot to index in gvectors */
            totpoints++;
         }

         /* vecsbelow[nbelow - 1], the last vector in this array, is the vector that is
          * immediately
          * smaller than tgttimes[itgt].  So, hcon_lookup(indices, tgttimes[itgt]) -->
          * index into gvectors of vector that is immediately smaller than tgttimes[itgt]
          */
         if (indices)
         {
            CreateLHashKey(lhashkey, sizeof(lhashkey), vecsbelow[nbelow - 1].slot);
            int *pvindex = (int *)hcon_lookup(slottovecindex, lhashkey);
            CreateDHashKey(dhashkey, sizeof(dhashkey), tgttimes[itgt]);

            if (pvindex)
            {
               hcon_insert(*indices, dhashkey, pvindex); /* store slot of vec that is smaller than
                                                          * tgttimes[itgt] */
            }
         }

         /* get the nabove points greater than tgttimes[itgt] */
         actpoints = 0;

         /* Must call Fetch(), not FetchNext() since the global cache item pointer may not
          * be at the grid vector immediately below tgttimes[itgt] */
         vec = Fetch(env, srcseries, optfilter, tgttimes[itgt], nbelow + 4, posabove);

         if (!vec)
         {
            fprintf(stderr,
                    "There is no grid vector greater than %f in '%s'.\n",
                    tgttimes[itgt],
                    srcseries);
            err = 1;
            break;
         }

         while (1)
         {

            slotpos = GetSlotPos(vec->obstime, tgttimes[itgt]);

            if (IsPos(slotpos, posabove))
            {
               vecsabove[actpoints] = *vec;
               actpoints++;
            }

            if (actpoints == nabove)
            {
               break;
            }

            /* need to fetch next grid vector */
            //vec = Fetch(env, srcseries, slotkey, optfilter, vec->obstime, 0, kGT);
            vec = FetchNext(env, srcseries, optfilter, nbelow + 4);
            if (!vec)
            {
               fprintf(stderr,
                       "Series '%s' missing a grid vector smaller than %f.\n",
                       srcseries,
                       tgttimes[itgt]);
               err = 1;
               break;
            }
         }

         if (err)
         {
            break;
         }

         /* okay to insert the nabove tgttimes */
         for (iabove = 0; iabove < nabove; iabove++)
         {
            CreateLHashKey(lhashkey, sizeof(lhashkey), vecsabove[iabove].slot);

            if (hcon_lookup(slottovecindex, lhashkey))
            {
               continue;
            }

            list_llinserttail(*gvectors, &(vecsabove[iabove]));
            hcon_insert(slottovecindex, lhashkey, &totpoints);
            totpoints++;
         }
      } /* for itgt */

      if (vecsbelow)
      {
         free(vecsbelow);
      }

      if (vecsabove)
      {
         free(vecsabove);
      }

      if (slottovecindex)
      {
         hcon_destroy(&slottovecindex);
      }
   }

   if (ngvecs)
   {
      *ngvecs = totpoints;
   }

   return err;
}

/* indices - for each tgttime, index into gridVecs of vector closest (LTE) to that tgttime */
static int iorbit_interpolate(IORBIT_Alg_t alg,
                              IORBIT_Vector_t *gridVecs, /* abcissae + ordinates */
                              const double *tgttimes,
                              HContainer_t *indexmap,
                              int ntgtpts,
                              IORBIT_Vector_t **interp)
{
   int err = 0;
   int ipt;
   int vindex = 0;
   char dhashkey[128];

   if (interp)
   {
      *interp = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * ntgtpts);

      switch (alg)
      {
         case IORBIT_Alg_Linear:
           {
              fprintf(stderr, "Unsupported interpolation algorithm '%d'.\n", (int)alg);
           }
           break;
         case IORBIT_Alg_Quadratic:
           {
              double ptlow;
              double ptmid;
              double pthii;
              IORBIT_Vector_t *veclow = NULL;
              IORBIT_Vector_t *vecmid = NULL;
              IORBIT_Vector_t *vechii = NULL;
              double factorlow;
              double factormid;
              double factorhii;

              for (ipt = 0; ipt < ntgtpts; ipt++)
              {
                  if (drms_ismissing_time(tgttimes[ipt]))
                  {
                      /* skip missing times */
                      ((*interp)[ipt]).obstime = DRMS_MISSING_TIME;
                      continue;
                  }

                  CreateDHashKey(dhashkey, sizeof(dhashkey), tgttimes[ipt]);
                  vindex = *((int *)hcon_lookup(indexmap, dhashkey));

                  ptlow = gridVecs[vindex - 1].obstime;
                  ptmid = gridVecs[vindex].obstime;
                  pthii = gridVecs[vindex + 1].obstime;
                  veclow = &(gridVecs[vindex - 1]);
                  vecmid = &(gridVecs[vindex]);
                  vechii = &(gridVecs[vindex + 1]);
                  factorlow = (tgttimes[ipt] - ptmid) * (tgttimes[ipt] - pthii) /
                  ((ptlow - ptmid) * (ptlow - pthii));
                  factormid = (tgttimes[ipt] - ptlow) * (tgttimes[ipt] - pthii) /
                  ((ptmid - ptlow) * (ptmid - pthii));
                  factorhii = (tgttimes[ipt] - ptlow) * (tgttimes[ipt] - ptmid) /
                  ((pthii - ptlow) * (pthii - ptmid));

                  ((*interp)[ipt]).obstime = tgttimes[ipt];
                  ((*interp)[ipt]).slot = DRMS_MISSING_LONGLONG;

                  /* Assume all orbit series have gae_x, gae_y, gae_z. */
                  ((*interp)[ipt]).gae_x =
                  veclow->gae_x * factorlow + vecmid->gae_x * factormid + vechii->gae_x * factorhii;
                  ((*interp)[ipt]).gae_y =
                  veclow->gae_y * factorlow + vecmid->gae_y * factormid + vechii->gae_y * factorhii;
                  ((*interp)[ipt]).gae_z =
                  veclow->gae_z * factorlow + vecmid->gae_z * factormid + vechii->gae_z * factorhii;

                  /* Not all orbit series have gae_vx. */
                  if (!drms_ismissing_double(veclow->gae_vx) && !drms_ismissing_double(vecmid->gae_vx) && !drms_ismissing_double(vechii->gae_vx))
                  {
                      ((*interp)[ipt]).gae_vx =
                      veclow->gae_vx * factorlow + vecmid->gae_vx * factormid + vechii->gae_vx * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).gae_vx = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all orbit series have gae_vy. */
                  if (!drms_ismissing_double(veclow->gae_vy) && !drms_ismissing_double(vecmid->gae_vy) && !drms_ismissing_double(vechii->gae_vy))
                  {
                      ((*interp)[ipt]).gae_vy =
                      veclow->gae_vy * factorlow + vecmid->gae_vy * factormid + vechii->gae_vy * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).gae_vy = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all orbit series have gae_vz. */
                  if (!drms_ismissing_double(veclow->gae_vz) && !drms_ismissing_double(vecmid->gae_vz) && !drms_ismissing_double(vechii->gae_vz))
                  {
                      ((*interp)[ipt]).gae_vz =
                      veclow->gae_vz * factorlow + vecmid->gae_vz * factormid + vechii->gae_vz * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).gae_vz = DRMS_MISSING_DOUBLE;
                  }

                  /* Assume all orbit series have hae_x, hae_y, hae_z. */
                  ((*interp)[ipt]).hae_x =
                  veclow->hae_x * factorlow + vecmid->hae_x * factormid + vechii->hae_x * factorhii;
                  ((*interp)[ipt]).hae_y =
                  veclow->hae_y * factorlow + vecmid->hae_y * factormid + vechii->hae_y * factorhii;
                  ((*interp)[ipt]).hae_z =
                  veclow->hae_z * factorlow + vecmid->hae_z * factormid + vechii->hae_z * factorhii;

                  /* Not all orbit series have hae_vx. */
                  if (!drms_ismissing_double(veclow->hae_vx) && !drms_ismissing_double(vecmid->hae_vx) && !drms_ismissing_double(vechii->hae_vx))
                  {
                      ((*interp)[ipt]).hae_vx =
                      veclow->hae_vx * factorlow + vecmid->hae_vx * factormid + vechii->hae_vx * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).hae_vx = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all orbit series have hae_vy. */
                  if (!drms_ismissing_double(veclow->hae_vy) && !drms_ismissing_double(vecmid->hae_vy) && !drms_ismissing_double(vechii->hae_vy))
                  {
                      ((*interp)[ipt]).hae_vy =
                      veclow->hae_vy * factorlow + vecmid->hae_vy * factormid + vechii->hae_vy * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).hae_vy = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all orbit series have hae_vz. */
                  if (!drms_ismissing_double(veclow->hae_vz) && !drms_ismissing_double(vecmid->hae_vz) && !drms_ismissing_double(vechii->hae_vz))
                  {
                      ((*interp)[ipt]).hae_vz =
                      veclow->hae_vz * factorlow + vecmid->hae_vz * factormid + vechii->hae_vz * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).hae_vz = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all orbit series have rsunobs. */
                  if (!drms_ismissing_double(veclow->rsunobs) && !drms_ismissing_double(vecmid->rsunobs) && !drms_ismissing_double(vechii->rsunobs))
                  {
                      ((*interp)[ipt]).rsunobs =
                      veclow->rsunobs * factorlow + vecmid->rsunobs * factormid + vechii->rsunobs * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).rsunobs = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all orbit series have obsvr. */
                  if (!drms_ismissing_double(veclow->obsvr) && !drms_ismissing_double(vecmid->obsvr) && !drms_ismissing_double(vechii->obsvr))
                  {
                      ((*interp)[ipt]).obsvr =
                      veclow->obsvr * factorlow + vecmid->obsvr * factormid + vechii->obsvr * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).obsvr = DRMS_MISSING_DOUBLE;
                  }

                  /* Not all series have dsunobs. */
                  if (!drms_ismissing_double(veclow->dsunobs) && !drms_ismissing_double(vecmid->dsunobs) && !drms_ismissing_double(vechii->dsunobs))
                  {
                      ((*interp)[ipt]).dsunobs =
                      veclow->dsunobs * factorlow + vecmid->dsunobs * factormid + vechii->dsunobs * factorhii;
                  }
                  else
                  {
                      ((*interp)[ipt]).dsunobs = DRMS_MISSING_DOUBLE;
                  }
              }
           }
           break;
         default:
           err = 1;
           fprintf(stderr, "Unsupported interpolation algorithm '%d'.\n", (int)alg);
      }
   }
   else
   {
      err = 1;
   }

   return err;
}

/* vec - array of interpolated vectors */
/* hb, hl in radians */
static int CalcSolarVelocities(IORBIT_Vector_t *vec,
                               int nvecs,
                               double **hr,      /* DSUN_OBS */
                               double **hvr,     /* OBS_VR */
                               double **hvw,     /* OBS_VW */
                               double **hvn,     /* OBS_VN */
                               double **hgln,   /* HGLN_OBS */
                               double **hb,      /* CRLT_OBS */
                               double **hl)      /* L0 */
{
   int err = 0;

   double hci_vx;
   double hci_vy;
   double hci_vz;

   double gci_x;
   double gci_y;

   const double txx = cos(kALPHA);
   const double txy = sin(kALPHA);
   const double txz = 0.0;
   const double tyx = -sin(kALPHA) * cos(kDELTA);
   const double tyy = cos(kALPHA) * cos(kDELTA);
   const double tyz = sin(kDELTA);
   const double tzx = sin(kALPHA) * sin(kDELTA);
   const double tzy = -cos(kALPHA) * sin(kDELTA);
   const double tzz = cos(kDELTA);

   *hr = (double *)malloc(sizeof(double) * nvecs);

   if (hvr)
   {
      *hvr = (double *)malloc(sizeof(double) * nvecs);
   }

   if (hvw)
   {
      *hvw = (double *)malloc(sizeof(double) * nvecs);
   }

   if (hvn)
   {
      *hvn = (double *)malloc(sizeof(double) * nvecs);
   }

   if (hgln)
   {
      *hgln = (double *)malloc(sizeof(double) * nvecs);
   }

   *hb = (double *)malloc(sizeof(double) * nvecs);
   *hl = (double *)malloc(sizeof(double) * nvecs);

   IORBIT_Vector_t *cvec = NULL;
   int ivec;
   for (ivec = 0; ivec < nvecs; ivec++)
   {
      cvec = &(vec[ivec]);

      /* distance to sun */
      (*hr)[ivec] = sqrt(cvec->hae_x * cvec->hae_x +
                         cvec->hae_y * cvec->hae_y +
                         cvec->hae_z * cvec->hae_z);


      /* beta and l0 angles */
      (*hb)[ivec] = asin((tzx * cvec->hae_x + tzy * cvec->hae_y + tzz * cvec->hae_z) / (*hr)[ivec]);
      (*hl)[ivec] = atan2((tyx * cvec->hae_x + tyy * cvec->hae_y + tyz * cvec->hae_z),
                        (txx * cvec->hae_x + txy * cvec->hae_y + txz * cvec->hae_z));

      /* velocities transformed to solar-rotation-axis coordinates */
      hci_vx = txx * cvec->hae_vx + txy * cvec->hae_vy + txz * cvec->hae_vz;
      hci_vy = tyx * cvec->hae_vx + tyy * cvec->hae_vy + tyz * cvec->hae_vz;
      hci_vz = tzx * cvec->hae_vx + tzy * cvec->hae_vy + tzz * cvec->hae_vz;

      /* position transformed to solar-rotation-axis coordinates */
      gci_x = txx * cvec->gae_x + txy * cvec->gae_y + txz * cvec->gae_z;
      gci_y = tyx * cvec->gae_x + tyy * cvec->gae_y + tyz * cvec->gae_z;

      if (hvr)
      {
         (*hvr)[ivec] = cos((*hb)[ivec]) * (cos((*hl)[ivec]) * hci_vx + sin((*hl)[ivec]) * hci_vy) + sin((*hb)[ivec]) * hci_vz;
      }

      if (hvw)
      {
         (*hvw)[ivec] = -sin((*hl)[ivec]) * hci_vx + cos((*hl)[ivec]) * hci_vy;
      }

      if (hvn)
      {
         (*hvn)[ivec] = -sin((*hb)[ivec]) * (cos((*hl)[ivec]) * hci_vx + sin((*hl)[ivec]) * hci_vy) + cos((*hb)[ivec]) * hci_vz;
      }

      if (hgln)
      {
          (*hgln)[ivec] = atan2(-sin((*hl)[ivec]) * gci_x + cos((*hl)[ivec]) * gci_y, (*hr)[ivec]);
      }
   }

   return err;
}

void iorbit_carrcoords(TIME t, double obsdist, double b, double hci_long, int *crot, double *L, double *B)
{
   SDOCarringtonCoords(t, obsdist, b, hci_long, crot, L, B);
}

static int SetUpKeymap(DRMS_Env_t *env, const char *series, HContainer_t *keymap)
{
    int err = 0;
    const char *keyid = NULL;
    const char *keyname = NULL;
    HIterator_t *hit = NULL;
    int dstat = DRMS_SUCCESS;
    DRMS_Record_t *template = NULL;

    if (keymap)
    {
        /* Iterate through elements, verifying that each key name is one of the expected names, and each value
         * is the name of a keyword in the orbit series. */
        hit = hiter_create(keymap);

        if (hit)
        {
            template = drms_template_record(env, series, &dstat);

            if (template && dstat == DRMS_SUCCESS)
            {
                while ((keyname = hiter_extgetnext(hit, &keyid)) != NULL)
                {
                    if (strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_X)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_Y)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_Z)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_VX)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_VY)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_VZ)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_X)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_Y)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_Z)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_VX)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_VY)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_VZ)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_DATE)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_RSUN_OBS)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_VR)) != 0 &&
                        strcasecmp(keyid, IORBIT_STRINGIFY(IORBIT_KEYWORD_DSUN_OBS)) != 0)
                    {
                        fprintf(stderr, "Invalid key ID %s.\n", keyid);
                        err = 1;
                        break;
                    }

                    if (!drms_keyword_lookup(template, keyname, 1))
                    {
                        fprintf(stderr, "Invalid keyword name %s.\n", keyname);
                        err = 1;
                        break;
                    }

                    if (!gKeyMap)
                    {
                        gKeyMap = hcon_create(DRMS_MAXKEYNAMELEN, DRMS_MAXKEYNAMELEN, NULL, NULL, NULL, NULL, 0);
                    }

                    if (!gKeyMap)
                    {
                        fprintf(stderr, "Out of memory.\n");
                        err = 1;
                    }
                    else
                    {
                        hcon_insert(gKeyMap, keyid, keyname);
                    }
                }

                if (err)
                {
                    hcon_destroy(&gKeyMap);
                    gKeyMap = NULL;
                }
            }
            else
            {
                fprintf(stderr, "Unable to get template record for series %s.\n", series);
                err = 1;
            }

            hiter_destroy(&hit);
        }
        else
        {
            fprintf(stderr, "Out of memory.\n");
            err = 1;
        }
    }

    return err;
}

LIBASTRO_Error_t iorbit_getinfo_internal(DRMS_Env_t *env,
                                         const char *srcseries,
                                         const char *optfilter,
                                         IORBIT_Alg_t alg,
                                         const double *tgttimes,
                                         int nitems,
                                         IORBIT_CacheAction_t ctype,
                                         IORBIT_Info_t **info,
                                         int calcSolarV)
{
    LIBASTRO_Error_t err = kLIBASTRO_Success;
    int ntimes = 0;
    double *tgtsorted = NULL;
    IORBIT_Info_t *retvec = NULL;

    if (info && nitems > 0)
    {
        *info = calloc(nitems, sizeof(IORBIT_Info_t));

        /* Ensure tgttimes are sorted in increasing value - missing values are removed */
        ntimes = SortTimes(tgttimes, nitems, &tgtsorted);
    }
    else
    {
        err = kLIBASTRO_InvalidArgs;
    }

    if (err == kLIBASTRO_Success && ntimes <= 0)
    {
        /* All target times are missing, and there are ntimes target times with missing values. */
        int ivec;

        for (ivec = 0; ivec < nitems; ivec++)
        {
            retvec = &((*info)[ivec]);

            /* Initialize with missing values */
            retvec->obstime = DRMS_MISSING_TIME;
            retvec->hae_x = DRMS_MISSING_DOUBLE;
            retvec->hae_y = DRMS_MISSING_DOUBLE;
            retvec->hae_z = DRMS_MISSING_DOUBLE;
            retvec->hae_vx = DRMS_MISSING_DOUBLE;
            retvec->hae_vy = DRMS_MISSING_DOUBLE;
            retvec->hae_vz = DRMS_MISSING_DOUBLE;
            retvec->gae_x = DRMS_MISSING_DOUBLE;
            retvec->gae_y = DRMS_MISSING_DOUBLE;
            retvec->gae_z = DRMS_MISSING_DOUBLE;
            retvec->gae_vx = DRMS_MISSING_DOUBLE;
            retvec->gae_vy = DRMS_MISSING_DOUBLE;
            retvec->gae_vz = DRMS_MISSING_DOUBLE;
            retvec->dsun_obs = DRMS_MISSING_DOUBLE;
            retvec->rsun_obs = DRMS_MISSING_DOUBLE;
            retvec->obs_vr = DRMS_MISSING_DOUBLE;
            retvec->obs_vw = DRMS_MISSING_DOUBLE;
            retvec->obs_vn = DRMS_MISSING_DOUBLE;
            retvec->crln_obs = DRMS_MISSING_DOUBLE;
            retvec->crlt_obs = DRMS_MISSING_DOUBLE;
            retvec->car_rot = DRMS_MISSING_INT;
            snprintf(retvec->orb_rec, sizeof(retvec->orb_rec), "%s", DRMS_MISSING_STRING);
        }
    }
    else if (err == kLIBASTRO_Success)
    {
        LinkedList_t *listvecs = NULL;
        HContainer_t *indexmap = NULL;
        IORBIT_Vector_t *vecs = NULL;
        IORBIT_Vector_t *interp = NULL;

        double *dsun_obs = NULL;
        double *hvr = NULL;
        double *hvw = NULL;
        double *hvn = NULL;
        double *hgln = NULL;

        int nbelow = 0;
        int nabove = 0;
        char filter[DRMS_MAXQUERYLEN];

        TIMER_t *tmr = NULL;

        gCacheChunking = 1;

        if (ctype == kIORBIT_CacheAction_Flush)
        {
            DestroyCache();
        }
        else if (ctype == kIORBIT_CacheAction_DontCache)
        {
            gCacheChunking = 0;
            DestroyCache(); /* The cache (which actually doesn't contain a chunk, but instead contains
                             * all data within the time range specified by the filter) should have been
                             * flushed at the end of each call to iorbit_getinfo(). But, just in case
                             * flush again. */
            if (optfilter != NULL)
            {
                fprintf(stderr, "WARNING: record-query filter '%s' is not applicable when caching is turned off and will not be used.\n", optfilter);
            }
            else
            {
                /* if !gCacheChunking, then must modify optfilter to restrict search to records that will likely
                 * contain all target points, plus extra points on either end. */
                double min = tgtsorted[0];
                double max = tgtsorted[ntimes - 1];
                char mintime[128];
                char maxtime[128];
                sprint_time(mintime, min - 3600, "UTC", 0);
                sprint_time(maxtime, max + 3600, "UTC", 0);

                snprintf(filter, sizeof(filter), "[%s-%s]", mintime, maxtime);
                optfilter = filter;
            }
        }

        switch (alg)
        {
            case IORBIT_Alg_Linear:
            {
                fprintf(stderr, "Unsupported interpolation algorithm '%d'.\n", (int)alg);
            }
                break;
            case IORBIT_Alg_Quadratic:
            {
                nbelow = 2;
                nabove = 1;
            }
                break;
            default:
                fprintf(stderr, "Unsupported interpolation algorithm '%d'.\n", (int)alg);
        }

        /* Use array of time strings to check for presence of gci/hci data in
         * grid cache - make some assumptions here
         *   1. The caller will ask for times in order.
         *   2. The user won't request the same target times under normal
         *      circumstances.
         *
         * Given these assumptions, the best way to refresh the cache is to
         * go in chronological order, using cached values until a miss happens.
         * Then at that point, the cache is no longer useful - so wipe it out, get
         * another chunk of records from DRMS, and use those values to populate the cache.
         */

        /* Get all grid vectors (data from the source series - each vector of information
         * is associated with an actual observation time.  The information is not
         * interpolated), in ascending order, such that the first grid vector's abscissa
         * (observation time) is the largest abscissa smaller than than the
         * smallest result abscissa (interpolated time), and the last grid vector's abscissa is
         * the smallest abscissa greater than the largest result abscissa.
         */

        if (env->verbose)
        {
            /* time obtaining the grid vectors */
            tmr = CreateTimer();
        }

        if (GetGridVectors(env,
                           srcseries, /* series containing observed data */
                           tgtsorted, /* times we want to interpolate to */
                           ntimes,    /* number of items in tgtsorted */
                           nbelow,    /* number of grid abscissae below first result abscissa */
                           kLTE,
                           nabove,    /* number of grid abscissae above last result abscissa */
                           kGT,
                           &listvecs,
                           NULL,
                           &indexmap, /* indices into listvecs; each target time is mapped to
                                       * an index into listvecs. The resulting vector is the
                                       * grid vector that has an obstime that is the largest
                                       * time smaller than the target time. */
                           optfilter))
        {
            err = kLIBASTRO_InsufficientData;
        }

        /* tgtsorted no longer needed */
        free(tgtsorted);

        if (env->verbose)
        {
            fprintf(stdout, "GetGridVectors() seconds elapsed: %f\n", GetElapsedTime(tmr));
            DestroyTimer(&tmr);
        }

        if (!err && listvecs && indexmap)
        {
            int ivec;
            ListNode_t *ln = NULL;

            /* Make an array out of the list of grid vectors */
            vecs = malloc(sizeof(IORBIT_Vector_t) * listvecs->nitems);
            list_llreset(listvecs);
            ivec = 0;

            while ((ln = list_llnext(listvecs)) != NULL)
            {
                vecs[ivec++] = *((IORBIT_Vector_t *)(ln->data));
            }

            if (env->verbose)
            {
                /* time interpolation */
                tmr = CreateTimer();
            }

            /* This operates on the original set of target times, unlike GetGridVectors. */
            if (!iorbit_interpolate(alg,
                                    vecs, /* grid vectors (contains both abscissa and ordinate values) */
                                    tgttimes, /* result abscissae */
                                    indexmap,
                                    nitems,
                                    &interp))
            {
                double *hb = NULL;
                double *hl = NULL;
                int crot;
                double l0;
                double b0;
                char dhashkey[128];
                int *pvindex = NULL;
                int vindex;
                char tbuf[48];

                if (env->verbose)
                {
                    fprintf(stdout, "iorbit_interpolate() seconds elapsed: %f\n", GetElapsedTime(tmr));
                    DestroyTimer(&tmr);

                    /* time solar velocity calculation */
                    tmr = CreateTimer();
                }

                if (calcSolarV)
                {
                    if (!CalcSolarVelocities(interp, nitems, &dsun_obs, &hvr, &hvw, &hvn, &hgln, &hb, &hl))
                    {
                        if (env->verbose)
                        {
                            fprintf(stdout, "CalcSolarVelocities() seconds elapsed: %f\n", GetElapsedTime(tmr));
                            DestroyTimer(&tmr);
                        }

                        /* loop through resulting vectors to create output */
                        for (ivec = 0; ivec < nitems; ivec++)
                        {
                            retvec = &((*info)[ivec]);

                            /* Carrington coordinates */
                            SDOCarringtonCoords(interp[ivec].obstime, dsun_obs[ivec], hb[ivec], hl[ivec], &crot, &l0, &b0);

                            if (!drms_ismissing_time(interp[ivec].obstime))
                            {
                                retvec->obstime = interp[ivec].obstime;

                                /* Return all positions and velocity scalars in m and m/s (input is km and km/s). */
                                retvec->hae_x = interp[ivec].hae_x * 1000;
                                retvec->hae_y = interp[ivec].hae_y * 1000;
                                retvec->hae_z = interp[ivec].hae_z * 1000;
                                retvec->hae_vx = interp[ivec].hae_vx * 1000;
                                retvec->hae_vy = interp[ivec].hae_vy * 1000;
                                retvec->hae_vz = interp[ivec].hae_vz * 1000;
                                retvec->gae_x = interp[ivec].gae_x * 1000;
                                retvec->gae_y = interp[ivec].gae_y * 1000;
                                retvec->gae_z = interp[ivec].gae_z * 1000;
                                retvec->gae_vx = interp[ivec].gae_vx * 1000;
                                retvec->gae_vy = interp[ivec].gae_vy * 1000;
                                retvec->gae_vz = interp[ivec].gae_vz * 1000;
                                retvec->dsun_obs = dsun_obs[ivec] * 1000;
                                retvec->obs_vr = hvr[ivec] * 1000;
                                retvec->obs_vw = hvw[ivec] * 1000;
                                retvec->obs_vn = hvn[ivec] * 1000;
                                retvec->hgln_obs = hgln[ivec] / kDEG2RAD;
                                retvec->rsun_obs = (retvec->dsun_obs < 0.001 && retvec->dsun_obs > -0.001 ? DRMS_MISSING_DOUBLE :
                                                    648000 * asin(kRSUNREF / retvec->dsun_obs) / kPI);
                                retvec->crln_obs = l0;
                                retvec->crlt_obs = b0;
                                retvec->car_rot = crot;

                                /* indexmap contains the indices into vecs array such that hcon_lookup(indexmap, dhashkey),
                                 * where dhashkey is the string equivalent of tgttimes[itgt],
                                 * returns the index into vecs that yields the vector whose tobs is immediately smaller. */

                                /* ok to use ivec - there is one tgttime per vector */
                                CreateDHashKey(dhashkey, sizeof(dhashkey), tgttimes[ivec]);

                                if ((pvindex = (int *)hcon_lookup(indexmap, dhashkey)) != NULL)
                                {
                                    vindex = *pvindex;
                                    sprint_time(tbuf, vecs[vindex].obstime, "UTC", 0);
                                    snprintf(retvec->orb_rec, sizeof(retvec->orb_rec), "%s[%s]", srcseries, tbuf);
                                }
                                else
                                {
                                    fprintf(stderr, "Unable to locate expected target time '%s' in indexmap hash.\n", dhashkey);
                                    err = kLIBASTRO_Internal;
                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    /* No velocities in the orbit series. */
                    /* loop through resulting vectors to create output */
                    for (ivec = 0; ivec < nitems; ivec++)
                    {
                        retvec = &((*info)[ivec]);

                        retvec->obstime = interp[ivec].obstime;

                        /* Set all velocity fields to missing - no velocity information in this branch. */
                        retvec->hae_vx = DRMS_MISSING_DOUBLE;
                        retvec->hae_vy = DRMS_MISSING_DOUBLE;
                        retvec->hae_vz = DRMS_MISSING_DOUBLE;
                        retvec->gae_vx = DRMS_MISSING_DOUBLE;
                        retvec->gae_vy = DRMS_MISSING_DOUBLE;
                        retvec->gae_vz = DRMS_MISSING_DOUBLE;
                        retvec->obs_vw = DRMS_MISSING_DOUBLE;
                        retvec->obs_vn = DRMS_MISSING_DOUBLE;
                        retvec->crln_obs = DRMS_MISSING_DOUBLE;
                        retvec->crlt_obs = DRMS_MISSING_DOUBLE;
                        retvec->car_rot = DRMS_MISSING_INT;

                        if (drms_ismissing_time(retvec->obstime))
                        {
                            retvec->hae_x = DRMS_MISSING_DOUBLE;
                            retvec->hae_y = DRMS_MISSING_DOUBLE;
                            retvec->hae_z = DRMS_MISSING_DOUBLE;
                            retvec->gae_x = DRMS_MISSING_DOUBLE;
                            retvec->gae_y = DRMS_MISSING_DOUBLE;
                            retvec->gae_z = DRMS_MISSING_DOUBLE;
                            retvec->rsun_obs = DRMS_MISSING_DOUBLE;
                            retvec->obs_vr = DRMS_MISSING_DOUBLE;
                            retvec->dsun_obs = DRMS_MISSING_DOUBLE;
                        }
                        else
                        {
                            /* Convert km to m. */
                            retvec->hae_x = interp[ivec].hae_x * 1000;
                            retvec->hae_y = interp[ivec].hae_y * 1000;
                            retvec->hae_z = interp[ivec].hae_z * 1000;
                            retvec->gae_x = interp[ivec].gae_x * 1000;
                            retvec->gae_y = interp[ivec].gae_y * 1000;
                            retvec->gae_z = interp[ivec].gae_z * 1000;
                            retvec->rsun_obs = interp[ivec].rsunobs;
                            retvec->obs_vr = interp[ivec].obsvr * 1000;
                            retvec->dsun_obs = interp[ivec].dsunobs * 1000;

                            /* indexmap contains the indices into vecs array such that hcon_lookup(indexmap, dhashkey),
                             * where dhashkey is the string equivalent of tgttimes[itgt],
                             * returns the index into vecs that yields the vector whose tobs is immediately smaller. */

                            /* ok to use ivec - there is one tgttime per vector */
                            CreateDHashKey(dhashkey, sizeof(dhashkey), tgttimes[ivec]);

                            if ((pvindex = (int *)hcon_lookup(indexmap, dhashkey)) != NULL)
                            {
                                vindex = *pvindex;
                                sprint_time(tbuf, vecs[vindex].obstime, "UTC", 0);
                                snprintf(retvec->orb_rec, sizeof(retvec->orb_rec), "%s[%s]", srcseries, tbuf);
                            }
                            else
                            {
                                fprintf(stderr, "Unable to locate expected target time '%s' in indexmap hash.\n", dhashkey);
                                err = kLIBASTRO_Internal;
                                break;
                            }
                        }
                    }
                }

                if (hb)
                {
                    free(hb);
                }

                if (hl)
                {
                    free(hl);
                }
            }
            else
            {
                /* error performing interpolation */
                fprintf(stderr, "Error performing interpolation.\n");
                err = kLIBASTRO_Interpolation;
            }
        }

        if (dsun_obs)
        {
            free(dsun_obs);
        }

        if (hvr)
        {
            free(hvr);
        }

        if (hvw)
        {
            free(hvw);
        }

        if (hvn)
        {
            free(hvn);
        }

        if (vecs)
        {
            free(vecs);
        }

        if (interp)
        {
            free(interp);
        }

        if (listvecs)
        {
            list_llfree(&listvecs);
        }

        if (indexmap)
        {
            hcon_destroy(&indexmap);
        }
    }

    if (ctype == kIORBIT_CacheAction_DontCache)
    {
        DestroyCache();
    }

    return err;
}

/*
 * srcseries - Datseries containing orbit-data information (like sdo.fds_orbit_vectors)
 * alg - interpolation algorithm
 * tgttimes - array of times (doubles) for which information is desired.
 * nitems - number of times in tgttimes.
 * optfilter - a record-set query that causes a subset of all orbit data to be used; if NULL
 *    all orbit data used.
 * flush - force a flush of the grid vector cache (XXX - if the caller wants a grid vector
 *    whose obstime is smaller than any item's in the cache, then the cache needs to be
 *    flushed)
 * info - the resulting list of IORBIT_Info_t structs (each contains position and velocity data
 *    for a given observation time.
 *
 */

/*
 * This version of the API function is used for SDO.
 */
LIBASTRO_Error_t iorbit_getinfo(DRMS_Env_t *env,
                                const char *srcseries,
                                const char *optfilter,
                                IORBIT_Alg_t alg,
                                const double *tgttimes,
                                int nitems,
                                IORBIT_CacheAction_t ctype,
                                IORBIT_Info_t **info)
{
    return iorbit_getinfo_internal(env, srcseries, optfilter, alg, tgttimes, nitems, ctype, info, 1);
}

/* Added for IRIS. IRIS has no velocity keywords. */

/*
 * For IRIS, the velocities (hae_vx, hae_vy, hae_vz, gae_vx, gae_vy, gae_vz) are not used. The values
 * for these quantities will set to missing.
 *
 * To use this API function, the caller must set-up a HContainer_t that maps the following strings
 *
 * IORBIT_KEYWORD_GAE_X
 * IORBIT_KEYWORD_GAE_Y
 * IORBIT_KEYWORD_GAE_Z
 * IORBIT_KEYWORD_GAE_VX
 * IORBIT_KEYWORD_GAE_VY
 * IORBIT_KEYWORD_GAE_VZ
 * IORBIT_KEYWORD_HAE_X
 * IORBIT_KEYWORD_HAE_Y
 * IORBIT_KEYWORD_HAE_Z
 * IORBIT_KEYWORD_HAE_VX
 * IORBIT_KEYWORD_HAE_VY
 * IORBIT_KEYWORD_HAE_VZ
 * IORBIT_KEYWORD_RSUN_OBS
 * IORBIT_KEYWORD_OBS_VR
 * IORBIT_KEYWORD_DSUN_OBS
 * IORBIT_KEYWORD_OBS_DATE
 *
 * to the corresponding series keywords names. Some series may not have one or more corresponding keywords.
 * For those keywords, either set the keyword name to the empty string, or omit the map entry altogether. If
 * keymap is NULL or the container it points to is empty, then the resulting behavior is identical to
 * calling iorbit_getinfo().
 */
LIBASTRO_Error_t iorbit_getinfo_ext(DRMS_Env_t *env,
                                    const char *srcseries,
                                    const char *optfilter,
                                    IORBIT_Alg_t alg,
                                    const double *tgttimes,
                                    int nitems,
                                    IORBIT_CacheAction_t ctype,
                                    IORBIT_Info_t **info,
                                    HContainer_t *keymap)
{
    if (SetUpKeymap(env, srcseries, keymap))
    {
        /* The only real error would be if the caller provided a map that contained an entry with
         * an unknown keyname, or a value that is not a keyword name.*/
        return kLIBASTRO_InvalidArgs;
    }

    return iorbit_getinfo_internal(env, srcseries, optfilter, alg, tgttimes, nitems, ctype, info, 0);
}

/* */
LIBASTRO_Error_t iorbit_getSaaHlzInfo(DRMS_Env_t *env, const char *series, const double *tgttimes, int nitems, IORBIT_SaaHlzInfo_t **info)
{
    LIBASTRO_Error_t err = kLIBASTRO_Success;
    int ntimes = 0;
    double *tgtsorted = NULL;
    int ipt;

    if (info && nitems > 0)
    {
        /* Ensure tgttimes are sorted in increasing value - missing values are removed */
        ntimes = SortTimesB(tgttimes, nitems, &tgtsorted);
    }
    else
    {
        err = kLIBASTRO_InvalidArgs;
    }

    *info = NULL;

    if (err == kLIBASTRO_Success)
    {
        /* drms_dmsv will create a db prepared statement, then it will */
        char stmnt[8192];
        DB_Binary_Result_t **res = NULL;
        void **values = NULL;
        unsigned int nargs;
        DB_Type_t dbtypes[2];

        snprintf(stmnt, sizeof(stmnt), "SELECT eventtype FROM %s WHERE recnum in (SELECT max(recnum) FROM %s WHERE ? >= starttime AND ? <= stoptime GROUP BY starttime,stoptime)", series, series);

        nargs = 2;
        values = calloc(1, nargs * sizeof(void *));

        if (!values)
        {
            err = kLIBASTRO_OutOfMemory;
        }

        if (!err)
        {
            values[0] = calloc(1, ntimes * sizeof(double));
            if (!values[0])
            {
                err = kLIBASTRO_OutOfMemory;
            }
        }

        if (!err)
        {
            values[1] = values[0]; /* Since we're passing the same column of target times twice to iorbit, cheat and make
                                    * the second column point to the first. */

            dbtypes[0] = DB_DOUBLE;
            dbtypes[1] = DB_DOUBLE;

            for (ipt = 0; ipt < ntimes; ipt++)
            {
                ((double *)(values[0]))[ipt] = tgtsorted[ipt];
            }
        }

        if (!err)
        {
            if ((res = drms_query_bin_ntuple(env->session, stmnt, ntimes, nargs, dbtypes, values)) == NULL)
            {
                err = kLIBASTRO_DbStatement;
            }
        }

        if (!err)
        {
            /* Fill-up return container. */
            *info = hcon_create(sizeof(HContainer_t *), IORBIT_SAAHLZINFO_TIME_KEY_LEN, (void (*)(const void *))hcon_destroy, NULL, NULL, NULL, 0);
            if (!*info)
            {
                err = kLIBASTRO_OutOfMemory;
            }
        }

        if (!err)
        {
            int nrows;
            int irow;
            char eventTypeStr[IORBIT_SAAHLZINFO_KW_EVENT_TYPE_VALUE_LEN];
            HContainer_t *colToList = NULL; /* Maps column (e.g., eventtype) to list. */
            LinkedList_t *listEventType = NULL;
            char timeStr[IORBIT_SAAHLZINFO_TIME_KEY_LEN];

            /* Loop over target times. */
            for (ipt = 0; ipt < ntimes; ipt++)
            {
                snprintf(timeStr, sizeof(timeStr), "%lf", tgtsorted[ipt]);
                nrows = res[ipt]->num_rows;

                /* Create a list for each SAA-HLZ keyword for each target time. */

                /* eventtype */
                listEventType = list_llcreate(IORBIT_SAAHLZINFO_KW_EVENT_TYPE_VALUE_LEN, NULL);

                /* others */

                if (!listEventType)
                {
                    /* Out of memory. */
                    err = kLIBASTRO_OutOfMemory;
                    break;
                }

                /* Create colToList map */
                colToList = hcon_create(sizeof(LinkedList_t *), IORBIT_SAAHLZINFO_TIME_KEY_LEN, (void (*)(const void *))list_llfree, NULL, NULL, NULL, 0);
                if (!colToList)
                {
                    /* Out of memory. */
                    err = kLIBASTRO_OutOfMemory;
                    break;
                }

                /* Loop over matching records in the SAA-HLZ table. */
                for (irow = 0; irow < nrows; irow++)
                {
                    /* eventtype (column 0) */
                    db_binary_field_getstr(res[ipt], irow, 0, sizeof(eventTypeStr), eventTypeStr);
                    list_llinserttail(listEventType, eventTypeStr);

                    /* other keywords (columns > 0) */
                }

                /* Insert all lists into colToList map (one element per column/list). */
                hcon_insert(colToList, IORBIT_SAAHLZINFO_KW_EVENT_TYPE, &listEventType);

                /* Insert lists into return-value container. */
                hcon_insert(*info, timeStr, &colToList);
            }
        }

        free(values[0]);
        free(values);
        db_free_binary_result_tuple(&res, ntimes);
    }

    /* tgtsorted no longer needed */
    free(tgtsorted);

    return err;
}


void iorbit_cleanup()
{
   DestroyCache();
}

void iorbit_cleanSaaHlzInfo(IORBIT_SaaHlzInfo_t **info)
{
    if (info && *info)
    {
        hcon_destroy(info);
    }
}
