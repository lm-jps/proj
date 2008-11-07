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

#define kPI M_PI
#define kDEG2RAD (kPI / 180)

/* Ecliptic values */
#define kALPHA (75.76 * kDEG2RAD)
#define kDELTA (7.25 * kDEG2RAD)

HContainer_t *gGridCache = NULL;

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

static inline int CreateLHashKey(char *hashkey, int size, long long slot)
{
   return snprintf(hashkey, size, "%lld", slot);
}

static inline int CreateDHashKey(char *hashkey, int size, double val)
{
   return snprintf(hashkey, size, "%f", val);
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
                          long long stslot)
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
      char query[256];
      int drmsstatus;
      char lhashkey[128];
  
      memset(&vec, sizeof(IORBIT_Vector_t), 0);

      snprintf(query, 
               sizeof(query), 
               "%s[%s=%lld/%d]", 
               srcseries, 
               kOBSDATE_INDEX, 
               stslot, 
               kCACHESIZE);

      rs = drms_open_records(env, query, &drmsstatus);

      nitems = rs->n;
  
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

         CreateLHashKey(lhashkey, sizeof(lhashkey), obsslot);
         hcon_insert(gGridCache, lhashkey, &vec);
      }
   }

   return nitems;
}

/* Determining the slot without having to go to psql, fetch records, and get it from a record 
 * involves calling the DRMS function that knows how to calculate slot numbers.  There is
 * some complexity involved when calculating the slot, so use drms_keyword_slotval2indexval().
 */
static long long GetSlot(DRMS_Keyword_t *slotkey, double tgttime)
{
#if 0
   /* This worked previously - keep until sure about drms_keyword_slotval2indexval*/
   long long slot1 = floor((tgttime - epoch + (step / 2.0)) / step);
#endif

   DRMS_Value_t valin;
   DRMS_Value_t valout;

   valin.type = DRMS_TYPE_TIME;
   valin.value.time_val = tgttime;

   if (drms_keyword_slotval2indexval(slotkey, &valin, &valout, NULL) != DRMS_SUCCESS)
   {
      fprintf(stderr, "Problem calculating the slot number from time '%f'.\n", tgttime);
      valout.value.longlong_val = DRMS_MISSING_TIME;
   }

   return valout.value.longlong_val;
}

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
           (posA == kEQ && posB == kGTE));
}


static IORBIT_Vector_t *Fetch(DRMS_Env_t *env, const char *srcseries, long long slot, int offset)
{
   IORBIT_Vector_t *vec = NULL;
   char lhashkey[128];

   CreateLHashKey(lhashkey, sizeof(lhashkey), slot);
   vec = LookupInCache(lhashkey);

   if (!vec)
   {
      FlushCache();
      
      if (RehydrateCache(env, srcseries, slot + offset) == 0)
      {
         fprintf(stderr, "Not enough data in '%s'.\n", srcseries);
      }
      else
      {
         vec = LookupInCache(lhashkey);
      }
   }

   return vec;
}

/* Returns the grid vectors needed in order to interpolate values for tgttimes.
 * indices - For each tgttime, returns the index into gvectors of vector whose 
 * grid point is the largest point that is less than or equal to the tgttime.
 * 
 */

static int GetGridVectors(DRMS_Env_t *env,
                          const char *srcseries,
                          DRMS_Keyword_t *slotkey,
                          double *tgttimes,
                          int ntgts,
                          int nbelow,
                          IORBIT_Slotpos_t posbelow,
                          int nabove,
                          IORBIT_Slotpos_t posabove,
                          LinkedList_t **gvectors,
                          int *ngvecs,
                          HContainer_t **indices)
{
   int err = 0;
   int actpoints;
   int totpoints = 0;
   IORBIT_Vector_t *vec = NULL;
   long long inslot;
   long long actualslot;
   IORBIT_Slotpos_t slotpos;
   IORBIT_Vector_t *vecsbelow = NULL;
   IORBIT_Vector_t *vecsabove = NULL;
   int itgt;
   int ibelow;
   int iabove;
   char lhashkey[128];
   char dhashkey[128];
   HContainer_t *slottovecindex = NULL;

   if (gvectors)
   {
      /* The grid of FDS orbit vectors is NOT regular (because of leap seconds). Usually, 
       * the times are 60 seconds apart, but sometimes they are 61 seconds apart. So, you can't 
       * assume that observation time of the data in a slotted key's slot is the time at the beginning
       * of that slot.
       *
       * So, make the grid-time vector cache keyed by slot number. Then 
       * map the tgttime to a slot, and check the cache for that slot. If it exists, then 
       * check the obsdate in the slot (it could fall anywhere within the slot because of 
       * the leap-seconds complication). If the obsdate is less than the tgttime, then
       * you have the closest value below the tgttime. If the obsdate is greater than the
       * tgttime, then you have the closest value above the tgttime. In this manner, 
       * you can get as many grid values as needed above and below the tgttime.
       */
 
      *gvectors = list_llcreate(sizeof(IORBIT_Vector_t));
      vecsbelow = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * nbelow);
      vecsabove = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * nabove);
      slottovecindex = hcon_create(sizeof(int), 128, NULL, NULL, NULL, NULL, 0);

      if (indices)
      {
         *indices = hcon_create(sizeof(int), 128, NULL, NULL, NULL, NULL, 0);
      }

      for (itgt = 0; itgt < ntgts; itgt++)
      {
         inslot = GetSlot(slotkey, tgttimes[itgt]);
         actualslot = inslot;

         /* get the nbelow points less than tgttimes[itgt] (in reverse order of course) */
         actpoints = 0;
         while (actpoints != nbelow)
         {
            /* Cache, starting with a slot nbelow + 8 slots before the actualslot.
             * Do this because we're walking backward in time in this for loop,
             * and we should walk backward about nbelow slots (8 is a buffer). */
            vec = Fetch(env, srcseries, actualslot, - nbelow - 8);
            
            if (!vec)
            {
               /* Just updated the cache, but got a miss - bail */
               fprintf(stderr, "Required grid point, slot number %lld, not in '%s'.\n", actualslot, srcseries);
               err = 1;
               break;
            }

            slotpos = GetSlotPos(vec->obstime, tgttimes[itgt]);

            if (IsPos(slotpos, posbelow))
            {
               vecsbelow[nbelow - actpoints - 1] = *vec;
               actpoints++;
            }

            actualslot--;
         }

         if (err)
         {
            break;
         }

         /* vec now points to the vector that is nbelow points below the tgttime. Insert
          * that vector into list of vectors to be returned if it is not already there.
          * Use the fact that the list of vectors must be monotonically increasing by
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

         /* vecsbelow[nbelow - 1], the last vector in this array, is the vector that is immediately 
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
         actualslot = inslot;
         actpoints = 0;
         while (actpoints != nabove)
         {
            vec = Fetch(env, srcseries, actualslot, 0);
            
            if (!vec)
            {
               /* Just updated the cache, but got a miss - bail */
               fprintf(stderr, "Required grid point, slot number %lld, not in '%s'.\n", actualslot, srcseries);
               err = 1;
               break;
            }

            slotpos = GetSlotPos(vec->obstime, tgttimes[itgt]);

            if (IsPos(slotpos, posabove))
            {
               vecsabove[actpoints] = *vec;
               actpoints++;
            }

            actualslot++;
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
                              int *indices, 
                              int ntgtpts,
                              IORBIT_Vector_t **interp)
{
   int err = 0;
   int ipt;

   if (interp)
   {
      *interp = (IORBIT_Vector_t *)malloc(sizeof(IORBIT_Vector_t) * ntgtpts);

      switch (alg)
      {
         case IORBIT_Alg_Linear:
           {
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
                 ptlow = gridVecs[indices[ipt] - 1].obstime;
                 ptmid = gridVecs[indices[ipt]].obstime;
                 pthii = gridVecs[indices[ipt] + 1].obstime;
                 veclow = &(gridVecs[indices[ipt] - 1]); 
                 vecmid = &(gridVecs[indices[ipt]]); 
                 vechii = &(gridVecs[indices[ipt] + 1]); 
                 factorlow = (tgttimes[ipt] - ptmid) * (tgttimes[ipt] - pthii) / 
                   ((ptlow - ptmid) * (ptlow - pthii));
                 factormid = (tgttimes[ipt] - ptlow) * (tgttimes[ipt] - pthii) / 
                   ((ptmid - ptlow) * (ptmid - pthii));
                 factorhii = (tgttimes[ipt] - ptlow) * (tgttimes[ipt] - ptmid) / 
                   ((pthii - ptlow) * (pthii - ptmid));

                 ((*interp)[ipt]).obstime = tgttimes[ipt];
                 ((*interp)[ipt]).slot = DRMS_MISSING_LONGLONG;

                 ((*interp)[ipt]).gciX = 
                   veclow->gciX * factorlow + vecmid->gciX * factormid + vechii->gciX * factorhii;
                 ((*interp)[ipt]).gciY = 
                   veclow->gciY * factorlow + vecmid->gciY * factormid + vechii->gciY * factorhii;
                 ((*interp)[ipt]).gciZ = 
                   veclow->gciZ * factorlow + vecmid->gciZ * factormid + vechii->gciZ * factorhii;

                 ((*interp)[ipt]).gciVX = 
                   veclow->gciVX * factorlow + vecmid->gciVX * factormid + vechii->gciVX * factorhii;
                 ((*interp)[ipt]).gciVY = 
                   veclow->gciVY * factorlow + vecmid->gciVY * factormid + vechii->gciVY * factorhii;
                 ((*interp)[ipt]).gciVZ = 
                   veclow->gciVZ * factorlow + vecmid->gciVZ * factormid + vechii->gciVZ * factorhii;

                 ((*interp)[ipt]).hciX = 
                   veclow->hciX * factorlow + vecmid->hciX * factormid + vechii->hciX * factorhii;
                 ((*interp)[ipt]).hciY = 
                   veclow->hciY * factorlow + vecmid->hciY * factormid + vechii->hciY * factorhii;
                 ((*interp)[ipt]).hciZ = 
                   veclow->hciZ * factorlow + vecmid->hciZ * factormid + vechii->hciZ * factorhii;

                 ((*interp)[ipt]).hciVX = 
                   veclow->hciVX * factorlow + vecmid->hciVX * factormid + vechii->hciVX * factorhii;
                 ((*interp)[ipt]).hciVY = 
                   veclow->hciVY * factorlow + vecmid->hciVY * factormid + vechii->hciVY * factorhii;
                 ((*interp)[ipt]).hciVZ = 
                   veclow->hciVZ * factorlow + vecmid->hciVZ * factormid + vechii->hciVZ * factorhii;
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
int CalcSolarVelocities(IORBIT_Vector_t *vec, int nvecs, double **hr, double **hvr, double **hvw, double **hvn)
{
   int err = 0;
   
   double *hb = NULL; /* radians */
   double *hl = NULL; /* radians */
   double shvx;
   double shvy;
   double shvz;

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

   hb = (double *)malloc(sizeof(double) * nvecs);
   hl = (double *)malloc(sizeof(double) * nvecs);

   IORBIT_Vector_t *cvec = NULL;
   int ivec;
   for (ivec = 0; ivec < nvecs; ivec++)
   {
      cvec = &(vec[ivec]);

      /* solar radius */
      (*hr)[ivec] = sqrt(cvec->hciX * cvec->hciX +
                         cvec->hciY * cvec->hciY +
                         cvec->hciZ * cvec->hciZ);


      /* beta and l0 angles */
      hb[ivec] = asin((tzx * cvec->hciX + tzy * cvec->hciY + tzz * cvec->hciZ) / (*hr)[ivec]);
      hl[ivec] = atan2((tyx * cvec->hciX + tyy * cvec->hciY + tyz * cvec->hciZ),  
                        (txx * cvec->hciX + txy * cvec->hciY + txz * cvec->hciZ));

      /* velocities transformed to solar-rotation-axis coordinates */
      shvx = txx * cvec->hciVX + txy * cvec->hciVY + txz * cvec->hciVZ;
      shvy = tyx * cvec->hciVX + tyy * cvec->hciVY + tyz * cvec->hciVZ;
      shvz = tzx * cvec->hciVX + tzy * cvec->hciVY + tzz * cvec->hciVZ;

      if (hvr)
      {
         (*hvr)[ivec] = cos(hb[ivec]) * (cos(hl[ivec]) * shvx + sin(hl[ivec]) * shvy) + sin(hb[ivec]) * shvz;
      }

      if (hvw)
      {
         (*hvw)[ivec] = -sin(hl[ivec]) * shvx + cos(hl[ivec]) * shvy;
      }

      if (hvn)
      {
         (*hvn)[ivec] = -sin(hb[ivec]) * (cos(hl[ivec]) * shvx + sin(hl[ivec]) * shvy) + cos(hb[ivec]) * shvz;
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

      LinkedList_t *listvecs = NULL;
      int *indices = NULL;
      HContainer_t *indexmap = NULL;
      IORBIT_Vector_t *vecs = NULL;
      IORBIT_Vector_t *interp = NULL;

      double *dsun_obs = NULL;
      double *hvr = NULL;
      double *hvw = NULL;
      double *hvn = NULL;

      IORBIT_Info_t retvec;

      /* Ensure tgttimes are sorted in increasing value - missing values are removed */
      ntimes = SortTimes(tgttimes, nitems, &tgtsorted);

      pkarr = drms_series_createpkeyarray(env, srcseries, &npkeys, &drmsstat);

      /* There should be only one prime key in the source series */
      XASSERT(pkarr && npkeys == 1);

      if (pkarr)
      {
         template = drms_template_record(env, srcseries, &drmsstat);
         slotkw = drms_keyword_lookup(template, pkarr[0], 0);
      }

      drms_series_destroypkeyarray(&pkarr, npkeys);

      if (slotkw)
      {
         *info = list_llcreate(sizeof(IORBIT_Info_t));

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


         /* Get all grid vectors, in ascending order, such that the first grid vector's abscissa
          * is the largest abscissa smaller than than the smallest result abscissa, 
          * and the last grid vector's abscissa is the smallest abscisaa greater than 
          * the largest result abscissa. 
          *
          *
          */

         if (GetGridVectors(env, 
                            srcseries, 
                            slotkw, 
                            tgtsorted,
                            ntimes,
                            2, /* number of grid abscissae below first result abscissa */
                            kLTE,
                            1, /* number of grid abscissae above last result abscissa */
                            kGT,
                            &listvecs,
                            NULL,
                            &indexmap))
         {
            err = kLIBASTRO_InsufficientData;
         }

         if (!err && listvecs && indexmap)
         {
            int ivec;
            ListNode_t *gv = NULL;
            int *vindex = NULL;
            char dhashkey[128];
            int itgt;

            /* Make an array out of the list of grid vectors */
            vecs = malloc(sizeof(IORBIT_Vector_t) * listvecs->nitems);
            list_llreset(listvecs);
            ivec = 0;

            while ((gv = list_llnext(listvecs)) != NULL)
            {
               vecs[ivec] = *((IORBIT_Vector_t *)(gv->data));
               ivec++;
            }

            /* Make an array of indices to gridvecs for tgttimes (the unsorted tgttimes) */
            indices = (int *)malloc(sizeof(int) * ntimes);

            for (itgt = 0; itgt < ntimes; itgt++)
            {
               CreateDHashKey(dhashkey, sizeof(dhashkey), tgttimes[itgt]);
               vindex = ((int *)hcon_lookup(indexmap, dhashkey));
               if (vindex)
               {
                  indices[itgt] = *vindex;
               }
            }

            if (!iorbit_interpolate(alg,
                                    vecs, /* grid vectors (contains both abscissa and ordinate values) */
                                    tgttimes, /* result abscissae */
                                    indices, 
                                    ntimes,
                                    &interp))
            {
               if (!CalcSolarVelocities(interp, ntimes, &dsun_obs, &hvr, &hvw, &hvn))
               {
                  /* loop through resulting vectors to create output */
                  for (ivec = 0; ivec < ntimes; ivec++)
                  {
                     retvec.obstime = interp[ivec].obstime;
                     retvec.hciX = interp[ivec].hciX;
                     retvec.hciY = interp[ivec].hciY;
                     retvec.hciZ = interp[ivec].hciZ;
                     retvec.hciVX = interp[ivec].hciVX;
                     retvec.hciVY = interp[ivec].hciVY;
                     retvec.hciVZ = interp[ivec].hciVZ;
                     retvec.dsun_obs = dsun_obs[ivec];
                     retvec.obs_vr = hvr[ivec];
                     retvec.obs_vw = hvw[ivec];
                     retvec.obs_vn = hvn[ivec];

                     list_llinserttail(*info, &retvec);
                  }
               }
            }
            else
            {
               /* error performing interpolation */
               fprintf(stderr, "Error performing interpolation.\n");
               err = kLIBASTRO_Interpolation;
            }
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
