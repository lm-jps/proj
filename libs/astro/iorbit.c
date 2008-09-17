#include "astro.h"
#include "jsoc_main.h"

#define DEBUG 0

#define kXGCI "X_GEO"
#define kYGCI "Y_GEO"
#define kZGCI "Z_GEO"
#define kVXGCI "VX_GEO"
#define kVYGCI "VY_GEO"
#define kVZGCI "VZ_GEO"
#define kXHCI "X_HELIO"
#define kYHCI "Y_HELIO"
#define kZHCI "Z_HELIO"
#define kVXHCI "VX_HELIO"
#define kVYHCI "VY_HELIO"
#define kVZHCI "VZ_HELIO"
#define kOBSDATE "OBS_DATE"
#define kOBSDATE_INDEX "OBS_DATE_index"

#define kCACHEKEYSIZE 64
#define kCACHESIZE 128

HContainer_t *gGridCache = NULL;

enum IORBIT_Slotpos_enum
{
   kMISSING = 0, 
   kLTE,
   kGT
};

typedef enum IORBIT_Slotpos_enum IORBIT_Slotpos_t;

struct IORBIT_Vector_struct
{
  long long slot; /* slot number */
  double obstime; 
  double gciX;
  double gciY;
  double gciZ;
  double gciVX;
  double gciVY;
  double gciVZ;
  double hciX;
  double hciY;
  double hciZ;
  double hciVX;
  double hciVY;
  double hciVZ;
};

typedef struct IORBIT_Vector_struct IORBIT_Vector_t;

/* Gets J2000.0 positions and velocities from cdf file.
   Units are km and km/s */
static void get_earth_ephem(double jd, double pos[6])
{
   long i;
   int targ=3;
   int ctr=11;
   double au=0.1495978706910000e09;

   pleph_(&jd,&targ,&ctr,pos);
   for (i=0;i<3;i++) pos[i]=pos[i]*au;
   for (i=3;i<6;i++) pos[i]=pos[i]*au/86400;
}

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
   double *sortedint = malloc(nitems * sizeof(double));
   memcpy(sortedint, unsorted, nitems * sizeof(double));
   qsort(sortedint, nitems, sizeof(double), CmpTimes);
   int itime = 0;
   int nret = 0;

   if (sorted)
   {
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

   return nret;
}

/* Cache functions 
 *
 * record-chunking will actually make this obsolete, but for now, just
 * cache the grid times.
 */

static inline int CreateHashKey(char *hashkey, int size, long long slot)
{
   return snprintf(hashkey, size, "%lld", slot);
}

static HContainer_t *CreateCache()
{
   if (!gGridCache)
   {
      gGridCache = hcon_create(sizeof(IORBIT_Vector_t), kCACHEKEYSIZE, NULL, NULL, NULL, NULL, 0);
   }

   return gGridCache;
}

static void DestroyCache()
{
   if (gGridCache)
   {
      hcon_destroy(&gGridCache);
   }
}

static IORBIT_Vector_t *LookupInCache(const char *key)
{
   IORBIT_Vector_t *val = NULL;

   if (gGridCache)
   {
      val = hcon_lookup(gGridCache, key);
   }

   return val;
}

static void FlushCache()
{
   DestroyCache();
   CreateCache();
}

static int RehydrateCache(DRMS_Env_t *env,
                          const char *srcseries, 
                          long long stslot,
                          int npoints)
{
   int nitems = 0; 

   if (!gGridCache)
   {
      gGridCache = CreateCache();
   }

   if (gGridCache)
   {
      IORBIT_Vector_t vec;
      long long obsslot;
      double obstime;
      DRMS_RecordSet_t *rs = NULL;
      DRMS_Record_t *rec = NULL; /* record of orbit data from FDS */
      int iitem = 0;
      //char sttimestr[64];
      //double duration;
      char query[256];
      int drmsstatus;
      char hashkey[128];
  
      memset(&vec, sizeof(IORBIT_Vector_t), 0);

      nitems = (npoints > kCACHESIZE) ? npoints : kCACHESIZE;
      //duration = nitems * step;
      //sprint_time(sttimestr, sttime, "TAI", 0);
      //snprintf(query, sizeof(query), "%s[%s/%fs]", srcseries, sttimestr, duration);
      snprintf(query, sizeof(query), "%s[%s=%lld/%d]", srcseries, kOBSDATE_INDEX, stslot, nitems);
      rs = drms_open_records(env, query, &drmsstatus);

      nitems = rs->n;

      XASSERT(nitems >= npoints);

      if (nitems >= npoints)
      {
         for (iitem = 0; iitem < nitems; iitem++)
         {
            rec = rs->records[iitem];

            obsslot = drms_getkey_longlong(rec, kOBSDATE_INDEX, NULL);
            vec.slot = obsslot;

            obstime = drms_getkey_double(rec, kOBSDATE, NULL);
            vec.obstime = obstime;

            /* These are in J2000.0, ecliptic coordinates 
             *   -gcipos and gcivel are values relative to earth
             *   -hcipos and hcivel are values relative to the sun
             */
            vec.gciX = drms_getkey_double(rec, kXGCI, NULL);
            vec.gciY = drms_getkey_double(rec, kYGCI, NULL);
            vec.gciZ = drms_getkey_double(rec, kZGCI, NULL);
            vec.gciVX = drms_getkey_double(rec, kVXGCI, NULL);
            vec.gciVY = drms_getkey_double(rec, kVYGCI, NULL);
            vec.gciVZ = drms_getkey_double(rec, kVZGCI, NULL);

            vec.hciX = drms_getkey_double(rec, kXHCI, NULL);
            vec.hciY = drms_getkey_double(rec, kYHCI, NULL);
            vec.hciZ = drms_getkey_double(rec, kZHCI, NULL);
            vec.hciVX = drms_getkey_double(rec, kVXHCI, NULL);
            vec.hciVY = drms_getkey_double(rec, kVYHCI, NULL);
            vec.hciVZ = drms_getkey_double(rec, kVZHCI, NULL);

            CreateHashKey(hashkey, sizeof(hashkey), obsslot);
            hcon_insert(gGridCache, hashkey, &vec);
         }
      }
      else
      {
         fprintf(stderr, "Insufficient data in '%s' to process slot '%lld'.\n", srcseries, stslot);
      }
   }

   return nitems;
}

static IORBIT_Slotpos_t GetSlotPos(double slottime, double tgttime)
{ 
   IORBIT_Slotpos_t slotpos;

   if (!drms_ismissing_time(slottime))
   {
      if (slottime <= tgttime)
      {
         slotpos = kLTE;
      }
      else
      {
         slotpos = kGT;
      }
   }
   else
   {
      slotpos = kMISSING;
   }

   return slotpos;
}

/* Returns the grid vectors needed in order to interpolate values to tgttime */
static int GetGridVectors(DRMS_Env_t *env,
                          const char *srcseries,
                          double tgttime, 
                          IORBIT_Alg_t alg, 
                          double epoch,
                          double step,
                          IORBIT_Vector_t **gvectors)
{
   int err = 0;
   int npoints = 0;
   int ipoint = 0;
   //double *gpointsint = NULL;
   char hashkey[128];
   IORBIT_Vector_t *vec = NULL;
   long long inslot;
   long long actualslot;
   //double closepoint;
   //double firstpoint;
   //int popcache = 0;
   IORBIT_Slotpos_t slotpos;
   int rehydrated = 0;

   if (gvectors)
   {
      switch (alg)
      {
         case IORBIT_Alg_Linear:
           npoints = 2;
           //gpointsint = malloc(sizeof(double) * npoints);
           break;
         default:
           fprintf(stderr, "Unsupported interpolation algorithm '%d'.\n", (int)alg);
      }

      /* The grid of FDS orbit vectors is NOT regular (because of leap seconds). Usually, 
       * the times are 60 seconds apart, but sometimes they are 61 seconds apart. So, you can't 
       * assume that observation time of the data in a slotted key's slot is the time at the beginning
       * of that slot.
       *
       * So, scrap what I was doing. Instead, make the cache keyed by slot number. Then 
       * map the tgttime to a slot, and check the cache for that slot. If it exists, then 
       * check the obsdate in the slot (it could fall anywhere within the slot because of 
       * the leap-seconds complication). If the obsdate is less than the tgttime, then
       * you have the closest value below the tgttime. If the obsdate is greater than the
       * tgttime, then you have the closest value above the tgttime. In this manner, 
       * you can get as many grid values as needed above and below the tgttime.


       */

      /* The grid is a TSEQ slotted series, where the epoch is in the middle of the 
       * zeroth slot, and the observation time is at the beginning of the slot. So, 
       * if you find the slot that the tgt time lies within, that slot contains the
       * data that is the closest, but smaller, time. */
      inslot = floor((tgttime - epoch + (step / 2.0)) / step);
      //      closepoint = inslot * step - (step / 2.0) + epoch; /* epoch lies in middle of zeroth slot */
   
      /* find first point */
      //firstpoint = closepoint;
      //ipoint = 1;

      //while (ipoint < npoints / 2)
      //{
      //   firstpoint -= step;
      //}

      //gpointsint[0] = firstpoint;
      //ipoint = 1;

      //while (ipoint < npoints)
      //{
      //   gpointsint[ipoint] = gpointsint[ipoint - 1] + step;
      //   ipoint++;
      //}

      /* with doubles in hand, convert to time strings so we can look them up in cache 
       * and get the corresponding gci and hci values.
       */     
      *gvectors = malloc(sizeof(IORBIT_Vector_t) * npoints);

      /* get the npoints / 2 points below tgttime (in reverse order of course) */
      for (ipoint = npoints / 2 - 1; ipoint >= 0 && !err; ipoint--)
      {
         slotpos = kMISSING;
         actualslot = inslot;

         while (slotpos != kLTE)
         {
            CreateHashKey(hashkey, sizeof(hashkey), actualslot);
            vec = LookupInCache(hashkey);

            if (!vec)
            {
               if (rehydrated)
               {
                  /* Just updated the cache, but got a miss - bail */
                  fprintf(stderr, "Required grid point '%s' (slot) not in '%s'.\n", hashkey, srcseries);
                  err = 1;
                  break;
               }

               FlushCache();
               /* Cache, starting with a slot (2 * npoints / 2) slots before the actualslot.
                * Do this because we're walking backward in time in this for loop. */
               if (RehydrateCache(env, srcseries, actualslot - 2 * npoints / 2, npoints / 2) < npoints / 2)
               {
                  fprintf(stderr, "Not enough data in '%s'.\n", srcseries);
                  err = 1;
                  break;
               }

               rehydrated = 1;
               
               /* Try again */
               continue;
            }
            else if (rehydrated)
            {
               rehydrated = 0;
            }

            slotpos = GetSlotPos(vec->obstime, tgttime);
            actualslot--;
         }

         if (!err)
         {
            /* We have a slot that has data, and the obsdate is LT tgttime */
            (*gvectors)[ipoint] = *vec;
         }
      }

      /* get the npoints / 2 points above tgttime (in chronological order) */
      for (ipoint = npoints / 2 ; ipoint < npoints && !err; ipoint++)
      {
         slotpos = kMISSING;
         actualslot = inslot;

         while (slotpos != kGT)
         {
            CreateHashKey(hashkey, sizeof(hashkey), actualslot);
            vec = LookupInCache(hashkey);

            if (!vec)
            {
               if (rehydrated)
               {
                  /* Just updated the cache, but got a miss - bail */
                  fprintf(stderr, "Required grid point '%s' (slot) not in '%s'.\n", hashkey, srcseries);
                  err = 1;
                  break;
               }

               FlushCache();
               /* Cache, starting with actualslot.
                * Do this because we're walking forward in time in this for loop. */
               if (RehydrateCache(env, srcseries, actualslot, npoints / 2) < npoints / 2)
               {
                  fprintf(stderr, "Not enough data in '%s'.\n", srcseries);
                  err = 1;
                  break;
               }
               
               rehydrated = 1;

               /* Try again */
               continue;
            }
            else if (rehydrated)
            {
               rehydrated = 0;
            }

            slotpos = GetSlotPos(vec->obstime, tgttime);
            actualslot++;
         }

         if (!err)
         {
            /* We have a slot that has data, and the obsdate is LT tgttime */
            (*gvectors)[ipoint] = *vec;
         }
      }
   }

   return npoints;
}

static int iorbit_interpolate(double time, 
                              IORBIT_Alg_t alg, 
                              IORBIT_Vector_t *gridVecs, 
                              IORBIT_Vector_t *interpvec)
{
   int err = 0;


   if (interpvec)
   {
      switch (alg)
      {
         case IORBIT_Alg_Linear:
           {
              int igrid = 0;
              double frac1;
              double frac2;

              frac2 = (time - gridVecs[igrid].obstime) /  
                (gridVecs[igrid + 1].obstime - gridVecs[igrid].obstime);
              frac1 = 1 - frac2;

              interpvec->gciX = frac1 * gridVecs[igrid].gciX + frac2 * gridVecs[igrid + 1].gciX;
              interpvec->gciY = frac1 * gridVecs[igrid].gciY + frac2 * gridVecs[igrid + 1].gciY;
              interpvec->gciZ = frac1 * gridVecs[igrid].gciZ + frac2 * gridVecs[igrid + 1].gciZ;
              interpvec->gciVX = frac1 * gridVecs[igrid].gciVX + frac2 * gridVecs[igrid + 1].gciVX;
              interpvec->gciVY = frac1 * gridVecs[igrid].gciVY + frac2 * gridVecs[igrid + 1].gciVY;
              interpvec->gciVZ = frac1 * gridVecs[igrid].gciVZ + frac2 * gridVecs[igrid + 1].gciVZ;
              interpvec->hciX = frac1 * gridVecs[igrid].hciX + frac2 * gridVecs[igrid + 1].hciX;
              interpvec->hciY = frac1 * gridVecs[igrid].hciY + frac2 * gridVecs[igrid + 1].hciY;
              interpvec->hciZ = frac1 * gridVecs[igrid].hciZ + frac2 * gridVecs[igrid + 1].hciZ;
              interpvec->hciVX = frac1 * gridVecs[igrid].hciVX + frac2 * gridVecs[igrid + 1].hciVX;
              interpvec->hciVY = frac1 * gridVecs[igrid].hciVY + frac2 * gridVecs[igrid + 1].hciVY;
              interpvec->hciVZ = frac1 * gridVecs[igrid].hciVZ + frac2 * gridVecs[igrid + 1].hciVZ;
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

/*
 * tgttimes - array of times (doubles) for which information is desired.
 * nitems - number of times in tgttimes.
 * 
 */
LIBASTRO_Error_t iorbit_getinfo(DRMS_Env_t *env, 
                                const char *srcseries, 
                                IORBIT_Alg_t alg,
                                const double *tgttimes, 
                                int nitems, 
                                LinkedList_t **info)
{
   LIBASTRO_Error_t err = kLIBASTRO_Success;
   int drmsstat = DRMS_SUCCESS;

   if (info)
   {
      double *tgtsorted = NULL;
      int ntimes = 0;
     
      char **pkarr = NULL;
      int npkeys = 0;
      DRMS_Record_t *template = NULL;
      DRMS_Keyword_t *slotkw = NULL;
      TIME epoch;
      double step;

      IORBIT_Vector_t *vecs = NULL;
      IORBIT_Vector_t interpvec;
      DRMS_SlotKeyUnit_t unit;
      int itime;
      int ngrid;

      /* Ensure tgttimes are sorted in increasing value - missing values are removed */
      ntimes = SortTimes(tgttimes, nitems, &tgtsorted);

      pkarr = drms_series_createpkeyarray(env, srcseries, &npkeys, &drmsstat);

      /* There should be only one prime key in the source series */
      XASSERT(pkarr && npkeys == 1);

      template = drms_template_record(env, srcseries, &drmsstat);
      slotkw = drms_keyword_lookup(template, pkarr[0], 0);
      epoch = drms_keyword_getslotepoch(slotkw, &drmsstat);
      step = drms_keyword_getslotstep(slotkw, &unit, &drmsstat);

      drms_series_destroypkeyarray(&pkarr, npkeys);

      /* Use array of time strings to check for presence of gci/hci data in 
       * grid cache - make some assumptions here 
       *   1. The caller will ask for times in order.
       *   2. The user won't request the same target times under normal
       *      circumstances.
       *
       * Given these assumptions, the best way to refresh the cache is to 
       * go in chronological order, using cached values until a miss happens.
       * Then at that point, the cache is no longer useful - so wipe it out, get all
       * subsequent values from DRMS, and use those values to populate the cache.
       */
      for (itime = 0; itime < ntimes; itime++)
      {
         ngrid = GetGridVectors(env, srcseries, tgtsorted[itime], alg, epoch, step, &vecs);

         if (vecs)
         {
            iorbit_interpolate(tgtsorted[itime], alg, vecs, &interpvec); 
         }
      }
      
   }
   else
   {
      err = kLIBASTRO_InvalidArgs;
   }
   
   return err;
}

LIBASTRO_Error_t iorbit_test(DRMS_Env_t *env, const char *rsquery)
{
   LIBASTRO_Error_t err = kLIBASTRO_Success;
   int drmsstat = DRMS_SUCCESS;
   char *query = strdup(rsquery);
   DRMS_RecordSet_t *rset = drms_open_records(env, query, &drmsstat);

   if (rset && rset->n > 0)
   {
      /* Get GCI values */
      DRMS_Record_t *rec = NULL;
      int irec;
      TIME obsdate;
      double **gcipos = malloc(sizeof(double *) * rset->n);
      double **gcivel = malloc(sizeof(double *) * rset->n);
      double **hcipos = malloc(sizeof(double *) * rset->n);
      double **hcivel = malloc(sizeof(double *) * rset->n);

      double pi=3.14159265358979e0;
      double deg2rad=pi/180;
      double gcitdt,jd,help[6];
      double ex,ey,ez,evx,evy,evz; /* earth relative to sun */
      double sx,sy,sz,svx,svy,svz; /* s/c relative to sun */

      double coordrot = 0; /* angular diff between equatorial and ecliptic */
      double ey2;
      double ez2;
      double evy2;
      double evz2;

      for (irec = 0; irec < rset->n; irec++)
      {
         rec = rset->records[irec];

         gcipos[irec] = malloc(sizeof(double) * 3);
         gcivel[irec] = malloc(sizeof(double) * 3);
         hcipos[irec] = malloc(sizeof(double) * 3);
         hcivel[irec] = malloc(sizeof(double) * 3);

         obsdate = drms_getkey_time(rec, kOBSDATE, NULL);

         gcipos[irec][0] = drms_getkey_double(rec, kXGCI, NULL);
         gcipos[irec][1] = drms_getkey_double(rec, kYGCI, NULL);
         gcipos[irec][2] = drms_getkey_double(rec, kZGCI, NULL);
         gcivel[irec][0] = drms_getkey_double(rec, kVXGCI, NULL);
         gcivel[irec][1] = drms_getkey_double(rec, kVYGCI, NULL);
         gcivel[irec][2] = drms_getkey_double(rec, kVZGCI, NULL);

         hcipos[irec][0] = drms_getkey_double(rec, kXHCI, NULL);
         hcipos[irec][1] = drms_getkey_double(rec, kYHCI, NULL);
         hcipos[irec][2] = drms_getkey_double(rec, kZHCI, NULL);
         hcivel[irec][0] = drms_getkey_double(rec, kVXHCI, NULL);
         hcivel[irec][1] = drms_getkey_double(rec, kVYHCI, NULL);
         hcivel[irec][2] = drms_getkey_double(rec, kVZHCI, NULL);

         gcitdt = obsdate + 32.184;
         jd=gcitdt/86400+2443144.5;

         get_earth_ephem(jd,help);
         /* results for the earth are for the same time as for SOHO, not
            after the appropriate time delay */
         ex=help[0];
         ey=help[1];
         ez=help[2];
         evx=help[3];
         evy=help[4];
         evz=help[5];

         /* eXX and eXXX (not x vals though) are all equatorial coordinates -
            convert to ecliptic coords. */
         coordrot = atan(evz/evy) - 23.45 * deg2rad;
         ey2 = ey * ey;
         ez2 = ez * ez;
         evy2 = evy * evy;
         evz2 = evz * evz;

         ey = sqrt(ey2 + ez2) * cos(coordrot);
         ez = sqrt(ey2 + ez2) * sin(coordrot);
         evy = sqrt(evy2 + evz2) * cos(coordrot);
         evz = sqrt(evy2 + evz2) * sin(coordrot);

         /* ecliptic coords */
         sx=ex+gcipos[irec][0];
         sy=ey+gcipos[irec][1];
         sz=ez+gcipos[irec][2];
         svx=evx+gcivel[irec][0];
         svy=evy+gcivel[irec][1];
         svz=evz+gcivel[irec][2];
       

#if DEBUG
         fprintf(stdout, "%-15s%-15s%-15s%-15s%-15s%-15s\n", 
                 "calc_evx", "eph_evx", "calc_evy", "eph_evy", "calc_evz", "eph_evz");
         fprintf(stdout, "%-15.8f%-15.8f%-15.8f%-15.8f%-15.8f%-15.8f\n", 
                 hcivel[irec][0] - gcivel[irec][0], 
                 evx, 
                 hcivel[irec][1] - gcivel[irec][1],
                 evy, 
                 hcivel[irec][2] - gcivel[irec][2],
                 evz);
         fprintf(stdout, "\n");
         fprintf(stdout, "%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                 "calc_ex", "eph_ex", "calc_ey", "eph_ey", "calc_ez", "eph_ez");
         fprintf(stdout, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                 hcipos[irec][0] - gcipos[irec][0], 
                 ex, 
                 hcipos[irec][1] - gcipos[irec][1],
                 ey, 
                 hcipos[irec][2] - gcipos[irec][2],
                 ez);
#endif

         /* print out hci values (calculated and from MOC) */
         fprintf(stdout, "Position of s/c relative to sun (J2000.0) == calc_XX - derived from gci, fds_XX - from FDS data products\n");
         fprintf(stdout, "%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                 "calc_x", "fds_x", "calc_y", "fds_y", "calc_z", "fds_z");
         fprintf(stdout, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                 sx, hcipos[irec][0], hcipos[irec][1], sy, sz, hcipos[irec][2]);
         
         fprintf(stdout, "\n");

         fprintf(stdout, "Velocity of s/c relative to sun (J2000.0) == calc_XX - derived from gci, fds_XX - from FDS data products\n");
         fprintf(stdout, "%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                 "calc_vx", "fds_vx", "calc_vy", "fds_vy", "calc_vz", "fds_vz");
         fprintf(stdout, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                 svx, hcivel[irec][0], svy, hcivel[irec][1], svz, hcivel[irec][2]);
      }

      if (gcipos)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (gcipos[irec])
            {
               free(gcipos[irec]);
            }
         }

         free(gcipos);
      }
      
      if (gcivel)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (gcivel[irec])
            {
               free(gcivel[irec]);
            }
         }

         free(gcivel);
      }

      if (hcipos)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (hcipos[irec])
            {
               free(hcipos[irec]);
            }
         }

         free(hcipos);
      }

      if (hcivel)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (hcivel[irec])
            {
               free(hcivel[irec]);
            }
         }

         free(hcivel);
      }
   }

   if (query)
   {
      free(query);
   }

   return err;
}

void iorbit_cleanup()
{
   DestroyCache();
}
