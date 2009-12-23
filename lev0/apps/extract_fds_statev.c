/* extractFdsStateV.c */

/* Extracts heliocentric state data, for the input date and from the input data series.
 *   Inputs:
 *     Data series in - the fully qualified name of the data series that contains the FDS data
 *     Time range - the range of dates to which the extraction should be restricted
 *     Data series out - the fully qualified name of the data series in which to store 
 *       the extracted data
 */

#include <stdio.h>
#include <stdlib.h>

#include "jsoc_main.h"
#include "cmdparams.h"
#include "drms_env.h"

#define DEBUG 0

#define kDefaultNamespace "sdo"
#define kDefaultSeriesIn "fds"
#define kDefaultSeriesOut "fds_orbit_vectors"
#define kSeriesHistory "fds_orbit_ingesthist"
#define kObsDateKey "OBS_DATE"
#define kFdsDataProductKey "FDS_DATA_PRODUCT"
#define kFdsProductCompKey "FDS_PRODUCT_COMP"
#define kFdsDataFormatKey "DATA_FORMAT"
#define kFileVersionKey "FILE_VERSION"
#define kFdsProductCompHELIO "ORBIT_HELIO"
#define kFdsProductCompGEO "ORBIT_GEO"

#define kSVFileSegName "FILENAME"

#define kSVLineMax 1024
#define kMaxHashKey 1024

/* State-vector data keys */
#define kSVKeyPrimary "OBS_DATE"
#define kSVKeyXHELIO "HCIEC_X"
#define kSVKeyYHELIO "HCIEC_Y"
#define kSVKeyZHELIO "HCIEC_Z"
#define kSVKeyVxHELIO "HCIEC_VX"
#define kSVKeyVyHELIO "HCIEC_VY"
#define kSVKeyVzHELIO "HCIEC_VZ"
#define kSVKeyXGEO "GCIEC_X"
#define kSVKeyYGEO "GCIEC_Y"
#define kSVKeyZGEO "GCIEC_Z"
#define kSVKeyVxGEO "GCIEC_VX"
#define kSVKeyVyGEO "GCIEC_VY"
#define kSVKeyVzGEO "GCIEC_VZ"

#define kSVKey_idHELIO "id_HELIO"
#define kSVKey_idGEO "id_GEO"

#define gChunkSize 16384
#define gCacheKeySize 128

#define kTOBSMaxLen 64
#define kIDMaxLen 256

struct VectorNode_struct 
{
  char tobs[kTOBSMaxLen];
  double hcix;
  double hciy;
  double hciz;
  double hcivx;
  double hcivy;
  double hcivz;
  char hciID[kIDMaxLen];
  double gcix;
  double gciy;
  double gciz;
  double gcivx;
  double gcivy;
  double gcivz;
  char gciID[kIDMaxLen];
};

typedef struct VectorNode_struct VectorNode_t;

typedef HContainer_t VectorCache_t;

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "ns", kDefaultNamespace, "working namespace (sdo, sdo_ground, sdo_dev)"},
  {ARG_STRING, "seriesIn", kDefaultSeriesIn, "name of series containing FDS data"},
  {ARG_STRING, "timeRange", "NOT SPECIFIED", "DRMS time interval encompassing data file dates"},
  {ARG_STRING, "seriesOut", kDefaultSeriesOut, "name of series in which to save extracted data"},
  {ARG_FLAG, "c", "0", "create the series only"},
  {ARG_END}
};

char *module_name = "extractFdsStateV";

int nice_intro ()
{
     int usage = cmdparams_get_int(&cmdparams, "h", NULL);
     if (usage)
     {
	  printf ("Usage:\n\textractFdsStateV [-h] "
		  "[seriesIn=<seriesname>] [timeRange=<timerange>] [seriesOut=<seriesname>]\n"
		  "  details are:\n"
		  "  -h: help - show this message then exit\n"
		  "  <seriesname> - fully qualified series name.\n"
		  "  <timerange> - time value range set.\n"
		  "  seriesIn defaults to sdo.fds.\n"
		  "  timeRange defaults to all records in seriesIn.\n"
		  "  seriesOut defaults to sdo.fdsStateVectors.\n"
		  "  example - extractFdsStateV seriesIn=su_arta.TestFDSData timeRange=2006.11.20_22:38:00-2006.11.20_22:45:00,2006.11.20_22:52:00. seriesOut=su_arta.TestFDSHelio\n");
	  return(1);
     }
     return (0);
}

static VectorCache_t *CreateVectorCache()
{
   VectorCache_t *cache = (VectorCache_t *)malloc(sizeof(VectorCache_t));
   hcon_init_ext((HContainer_t *)cache, 49999, sizeof(VectorNode_t), kTOBSMaxLen, NULL, NULL);
   return cache;
}

static void DestroyVectorCache(VectorCache_t **cache)
{
   if (cache)
   {
      if (*cache)
      {
         hcon_destroy((HContainer_t **)cache);
      }
   }
}

static void VCacheCache(VectorCache_t *cache,
                        char *tobs,
                        double *hcix,
                        double *hciy,
                        double *hciz,
                        double *hcivx,
                        double *hcivy,
                        double *hcivz,
                        char *hciID,
                        double *gcix,
                        double *gciy,
                        double *gciz,
                        double *gcivx,
                        double *gcivy,
                        double *gcivz,
                        char *gciID)
{
   if (cache && tobs && *tobs)
   {
#if DEBUG
      /* the time to alloc WAS growing */
      TIMER_t *timer = CreateTimer();
      VectorNode_t *node = (VectorNode_t *)hcon_allocslot_lower(cache, tobs);
      fprintf(stdout, "alloc slot: %f seconds.\n", GetElapsedTime(timer));
      DestroyTimer(&timer);
#else
      VectorNode_t *node = (VectorNode_t *)hcon_allocslot_lower(cache, tobs);
#endif
      snprintf(node->tobs, kTOBSMaxLen, "%s", tobs);

      if (hcix)
      {
         node->hcix = *hcix;
      }
      else
      {
         node->hcix = DRMS_MISSING_DOUBLE;
      }
      if (hciy)
      {
         node->hciy = *hciy;
      }
      else
      {
         node->hciy = DRMS_MISSING_DOUBLE;
      }
      if (hciz)
      {
         node->hciz = *hciz;
      }
      else
      {
         node->hciz = DRMS_MISSING_DOUBLE;
      }
      if (hcivx)
      {
         node->hcivx = *hcivx;
      }
      else
      {
         node->hcivx = DRMS_MISSING_DOUBLE;
      }
      if (hcivy)
      {
         node->hcivy = *hcivy;
      }
      else
      {
         node->hcivy = DRMS_MISSING_DOUBLE;
      }
      if (hcivz)
      {
         node->hcivz = *hcivz;
      }
      else
      {
         node->hcivz = DRMS_MISSING_DOUBLE;
      }

      if (hciID && *hciID)
      {
         snprintf(node->hciID, kIDMaxLen, "%s", hciID);
      }
      else
      {
         snprintf(node->hciID, kIDMaxLen, "%s", "unknown");
      }

      if (gcix)
      {
         node->gcix = *gcix;
      }
      else
      {
         node->gcix = DRMS_MISSING_DOUBLE;
      }
      if (gciy)
      {
         node->gciy = *gciy;
      }
      else
      {
         node->gciy = DRMS_MISSING_DOUBLE;
      }
      if (gciz)
      {
         node->gciz = *gciz;
      }
      else
      {
         node->gciz = DRMS_MISSING_DOUBLE;
      }
      if (gcivx)
      {
         node->gcivx = *gcivx;
      }
      else
      {
         node->gcivx = DRMS_MISSING_DOUBLE;
      }
      if (gcivy)
      {
         node->gcivy = *gcivy;
      }
      else
      {
         node->gcivy = DRMS_MISSING_DOUBLE;
      }
      if (gcivz)
      {
         node->gcivz = *gcivz;
      }
      else
      {
         node->gcivz = DRMS_MISSING_DOUBLE;
      }

      if (gciID && *gciID)
      {
         snprintf(node->gciID, kIDMaxLen, "%s", gciID);
      }
      else
      {
         snprintf(node->gciID, kIDMaxLen, "%s", "unknown");
      }
   }
}

static VectorNode_t *VCacheLookup(VectorCache_t *cache, const char *tstr)
{
   VectorNode_t *node = NULL;

   if (tstr && *tstr)
   {
      node = (VectorNode_t *)hcon_lookup_lower(cache, tstr);
   }

   return node;
}

static DRMS_Keyword_t *AddKey(DRMS_Record_t *prototype, 
			      DRMS_Type_t type, 
			      const char *name,
			      const char *format,
			      const char *unit,
			      DRMS_RecScopeType_t scope,
			      const char *desc,
			      int isprime)
{
   DRMS_Keyword_t *tKey = NULL;

   XASSERT(tKey = 
	   hcon_allocslot_lower(&(prototype->keywords), name));
   memset(tKey, 0, sizeof(DRMS_Keyword_t));
   XASSERT(tKey->info = malloc(sizeof(DRMS_KeywordInfo_t)));
   memset(tKey->info, 0, sizeof(DRMS_KeywordInfo_t));
	    
   if (tKey && tKey->info)
   {
      /* record */
      tKey->record = prototype;

      /* keyword info */
      snprintf(tKey->info->name, 
	       DRMS_MAXKEYNAMELEN,
	       "%s",
	       name);
      tKey->info->type = type;
      snprintf(tKey->info->format, DRMS_MAXFORMATLEN, "%s", format);
      snprintf(tKey->info->unit, DRMS_MAXUNITLEN, "%s", unit);
      tKey->info->recscope = scope;

      if (isprime)
      {
         drms_keyword_setintprime(tKey);
      }
      else
      {
         drms_keyword_unsetintprime(tKey);
      }

      snprintf(tKey->info->description, DRMS_MAXCOMMENTLEN, "%s", desc);

      /* default value - missing */
      drms_missing(type, &(tKey->value));
   }

   return tKey;
}

static void CreateOutSeries(DRMS_Env_t *drmsEnv, char *outSeries, int *status)
{
   int stat = DRMS_SUCCESS;

   if (!drms_series_exists(drmsEnv, outSeries, &stat))
   {
      DRMS_Record_t *prototype = (DRMS_Record_t *)calloc(1, sizeof(DRMS_Record_t));
	    
      if (prototype)
      {	       
	 prototype->seriesinfo = calloc(1, sizeof(DRMS_SeriesInfo_t));

	 if (prototype->seriesinfo)
	 {
	    DRMS_Keyword_t *pkey = NULL;
	    DRMS_Keyword_t *akey = NULL;
	    char keyname[DRMS_MAXKEYNAMELEN];

	    prototype->env = drmsEnv;
	    prototype->recnum = 0;
	    prototype->sunum = -1;		 
	    prototype->init = 1;
	    prototype->sessionid = 0;
	    prototype->sessionns = NULL;
	    prototype->su = NULL;
      
	    /* Initialize container structure. */
	    hcon_init(&prototype->segments, sizeof(DRMS_Segment_t), DRMS_MAXHASHKEYLEN, 
		      (void (*)(const void *)) drms_free_segment_struct, 
		      (void (*)(const void *, const void *)) drms_copy_segment_struct);
	    /* Initialize container structures for links. */
	    hcon_init(&prototype->links, sizeof(DRMS_Link_t), DRMS_MAXHASHKEYLEN, 
		      (void (*)(const void *)) drms_free_link_struct, 
		      (void (*)(const void *, const void *)) drms_copy_link_struct);
	    /* Initialize container structure. */
	    hcon_init(&prototype->keywords, sizeof(DRMS_Keyword_t), DRMS_MAXHASHKEYLEN, 
		      (void (*)(const void *)) drms_free_keyword_struct, 
		      (void (*)(const void *, const void *)) drms_copy_keyword_struct);

	    /* series info */
	    char *user = getenv("USER");
	    snprintf(prototype->seriesinfo->seriesname,
		     DRMS_MAXSERIESNAMELEN,
		     "%s",
		     outSeries);

	    strcpy(prototype->seriesinfo->author, "unknown");
	    strcpy(prototype->seriesinfo->owner, "unknown");

	    if (user)
	    {
	       if (strlen(user) < DRMS_MAXCOMMENTLEN)
	       {
		  strcpy(prototype->seriesinfo->author, user);
	       }
		    
	       if (strlen(user) < DRMS_MAXOWNERLEN)
	       {
		  strcpy(prototype->seriesinfo->owner, user);
	       }
	    }

	    /* discard "Owner", fill it with the dbuser */
	    if (drmsEnv->session->db_direct) 
	    {
	       strcpy(prototype->seriesinfo->owner, drmsEnv->session->db_handle->dbuser);
	    }

	    prototype->seriesinfo->unitsize = 0;
            prototype->seriesinfo->tapegroup = 1;

	    snprintf(prototype->seriesinfo->description,
		     DRMS_MAXCOMMENTLEN,
		     "%s",
		     "Helio- and Geo-centric orbit position and velocity vectors.  Prime key is observation date.");

	    /* slotted keyword isdrmsprime, but not really prime */
	    akey = AddKey(prototype, DRMS_TYPE_TIME, kSVKeyPrimary, "0", "UTC", kRecScopeType_TS_EQ, "slotted observation date", 1);

	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "index");
	    pkey = AddKey(prototype, kIndexKWType, keyname, kIndexKWFormat, "none", kRecScopeType_Index, "index for OBS_DATE", 1);
            drms_keyword_setimplicit(pkey);

	    /* epoch */
	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "epoch");
	    akey = AddKey(prototype, DRMS_TYPE_TIME, keyname, "0", "UTC", kRecScopeType_Constant, "epoch for OBS_DATE", 0);
	    akey->value.time_val = sscan_time("1993.01.01_00:00:30_UT");

	    /* step */
	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "step");
	    akey = AddKey(prototype, DRMS_TYPE_STRING, keyname, "%s", "none", kRecScopeType_Constant, "time step for OBS_DATE", 0);
	    copy_string(&(akey->value.string_val), "1m");

	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyXHELIO, "%f", "km", kRecScopeType_Variable, "X position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyYHELIO, "%f", "km", kRecScopeType_Variable, "Y position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyZHELIO, "%f", "km", kRecScopeType_Variable, "Z position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVxHELIO, "%f", "km/sec", kRecScopeType_Variable, "velocity along X", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVyHELIO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Y", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVzHELIO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Z", 0);

	    akey = AddKey(prototype, DRMS_TYPE_STRING, kSVKey_idHELIO, "%s", "unique id", kRecScopeType_Variable, "sdo.fds prime key values to locate HELIO source", 0);
	    copy_string(&(akey->value.string_val), "unknown");

	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyXGEO, "%f", "km", kRecScopeType_Variable, "X position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyYGEO, "%f", "km", kRecScopeType_Variable, "Y position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyZGEO, "%f", "km", kRecScopeType_Variable, "Z position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVxGEO, "%f", "km/sec", kRecScopeType_Variable, "velocity along X", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVyGEO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Y", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVzGEO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Z", 0);

	    akey = AddKey(prototype, DRMS_TYPE_STRING, kSVKey_idGEO, "%s", "unique id", kRecScopeType_Variable, "sdo.fds prime key values to locate GEO source", 0);
	    copy_string(&(akey->value.string_val), "unknown");

	    prototype->seriesinfo->pidx_keywords[0] = pkey;
	    prototype->seriesinfo->pidx_num = 1;
	 }

	 stat = drms_create_series_fromprototype(&prototype, outSeries, 0);
      }
   }

   if (status)
   {
      *status = stat;
   }
}

static void CreateHistSeries(DRMS_Env_t *drmsEnv, char *histSeries, int *status)
{
   int stat = DRMS_SUCCESS;

   if (!drms_series_exists(drmsEnv, histSeries, &stat))
   {
      DRMS_Record_t *prototype = (DRMS_Record_t *)calloc(1, sizeof(DRMS_Record_t));
	    
      if (prototype)
      {	       
	 prototype->seriesinfo = calloc(1, sizeof(DRMS_SeriesInfo_t));

	 if (prototype->seriesinfo)
	 {
	    DRMS_Keyword_t *pkey1 = NULL;
	    DRMS_Keyword_t *pkey2 = NULL;
	    DRMS_Keyword_t *pkey3 = NULL;
	    DRMS_Keyword_t *pkey4 = NULL;
	    DRMS_Keyword_t *akey = NULL;
	    char keyname[DRMS_MAXKEYNAMELEN];

	    prototype->env = drmsEnv;
	    prototype->recnum = 0;
	    prototype->sunum = -1;		 
	    prototype->init = 1;
	    prototype->sessionid = 0;
	    prototype->sessionns = NULL;
	    prototype->su = NULL;
      
	    /* Initialize container structure. */
	    hcon_init(&prototype->segments, sizeof(DRMS_Segment_t), DRMS_MAXHASHKEYLEN, 
		      (void (*)(const void *)) drms_free_segment_struct, 
		      (void (*)(const void *, const void *)) drms_copy_segment_struct);
	    /* Initialize container structures for links. */
	    hcon_init(&prototype->links, sizeof(DRMS_Link_t), DRMS_MAXHASHKEYLEN, 
		      (void (*)(const void *)) drms_free_link_struct, 
		      (void (*)(const void *, const void *)) drms_copy_link_struct);
	    /* Initialize container structure. */
	    hcon_init(&prototype->keywords, sizeof(DRMS_Keyword_t), DRMS_MAXHASHKEYLEN, 
		      (void (*)(const void *)) drms_free_keyword_struct, 
		      (void (*)(const void *, const void *)) drms_copy_keyword_struct);

	    /* series info */
	    char *user = getenv("USER");
	    snprintf(prototype->seriesinfo->seriesname,
		     DRMS_MAXSERIESNAMELEN,
		     "%s",
		     histSeries);

	    strcpy(prototype->seriesinfo->author, "unknown");
	    strcpy(prototype->seriesinfo->owner, "unknown");

	    if (user)
	    {
	       if (strlen(user) < DRMS_MAXCOMMENTLEN)
	       {
		  strcpy(prototype->seriesinfo->author, user);
	       }
		    
	       if (strlen(user) < DRMS_MAXOWNERLEN)
	       {
		  strcpy(prototype->seriesinfo->owner, user);
	       }
	    }

	    /* discard "Owner", fill it with the dbuser */
	    if (drmsEnv->session->db_direct) 
	    {
	       strcpy(prototype->seriesinfo->owner, drmsEnv->session->db_handle->dbuser);
	    }

	    prototype->seriesinfo->unitsize = 0;
	    prototype->seriesinfo->tapegroup = 1;

	    snprintf(prototype->seriesinfo->description,
		     DRMS_MAXCOMMENTLEN,
		     "%s",
		     "Helio- and Geo-centric orbit position and velocity vectors ingestion history.");

	    /* slotted keyword isdrmsprime, but not really prime */
	    akey = AddKey(prototype, DRMS_TYPE_TIME, kObsDateKey, "0", "UTC", kRecScopeType_TS_EQ, "slotted observation date", 1);

	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "index");
	    pkey1 = AddKey(prototype, kIndexKWType, keyname, kIndexKWFormat, "none", kRecScopeType_Index, "index for OBS_DATE", 1);
            drms_keyword_setimplicit(pkey1);

	    /* epoch */
	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "epoch");
	    akey = AddKey(prototype, DRMS_TYPE_TIME, keyname, "0", "UTC", kRecScopeType_Constant, "epoch for OBS_DATE", 0);
	    akey->value.time_val = sscan_time("1993.01.01_12:00:00_UT");

	    /* step */
	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "step");
	    akey = AddKey(prototype, DRMS_TYPE_STRING, keyname, "%s", "none", kRecScopeType_Constant, "time step for OBS_DATE", 0);
	    copy_string(&(akey->value.string_val), "1440m");

            pkey2 = AddKey(prototype, DRMS_TYPE_STRING, kFdsProductCompKey, "%s", "none", kRecScopeType_Variable, "FDS data product component", 1);
            pkey3 = AddKey(prototype, DRMS_TYPE_STRING, kFdsDataFormatKey, "%s", "none", kRecScopeType_Variable, "Format of data file contained in segment", 1);
            pkey4 = AddKey(prototype, DRMS_TYPE_INT, kFileVersionKey, "%d", "none", kRecScopeType_Variable, "FDS data product", 1);

	    prototype->seriesinfo->pidx_keywords[0] = pkey1;
	    prototype->seriesinfo->pidx_keywords[1] = pkey2;
	    prototype->seriesinfo->pidx_keywords[2] = pkey3;
	    prototype->seriesinfo->pidx_keywords[3] = pkey4;
	    prototype->seriesinfo->pidx_num = 4;
	 }

	 stat = drms_create_series_fromprototype(&prototype, histSeries, 0);
      }
   }

   if (status)
   {
      *status = stat;
   }
}

static void CloseCachedRecords(HContainer_t **cache)
{
   if (cache && *cache)
   {
      DRMS_RecordSet_t *rs = *((DRMS_RecordSet_t **)hcon_lookup(*cache, "bobmould"));
      drms_close_records(rs, DRMS_FREE_RECORD);
      rs = NULL;

      /* destory cache */
      hcon_destroy(cache);
   }
}

/* This uses record-chunking to find the correct record. Assumes records are in increasing time order,
 * and calls to this function have increaing tbuf */
/* Caller must clean up records by calling CloseCachedRecords() */
static int FetchCachedRecord(DRMS_Env_t *drmsEnv, 
                             char *series, 
                             char *tbuf, 
                             HContainer_t **cache,
                             DRMS_Record_t **recout)
{
   int error = 0;
   DRMS_RecordSet_t *rs = NULL;
   int stop;
   int rehydrated;
   DRMS_Record_t **prec = NULL;
   DRMS_Record_t *rec = NULL;
   char *timestr = NULL;
   int status = DRMS_SUCCESS;
   DRMS_RecChunking_t cstat = kRecChunking_None;

   if (!cache || !recout)
   {
      error = 1;
   }
   else
   {
      *recout = NULL;

      if (*cache == NULL)
      {
         drms_recordset_setchunksize(gChunkSize);
         rs = drms_open_recordset(drmsEnv, series, &status);

         if (status != DRMS_SUCCESS || rs == NULL) 
         {
            fprintf(stderr, "drms_open_recordset() failed, query=%s, error=%d.  Aborting.\n", series, error);
         }
      }
      /* attempt to find the record that matches tbuf */
      else
      { 
         /* retrieve saved rs */
         rs = *((DRMS_RecordSet_t **)hcon_lookup(*cache, "bobmould"));

         if (rs->n > 0)
         {
            if ((prec = (DRMS_Record_t **)hcon_lookup(*cache, tbuf)) != NULL)
            {
               *recout = *prec;
            }
         }
      }

      if (rs->n > 0)
      {
         /* Have to loop on ALL chunks if can't find time in cache */
         int morerecs = 1;
         while (!prec && !error && morerecs)
         {
            /* either no cache, or cache, but miss - get more records, then try to find tbuf */
            stop = 0;
            rehydrated = 0;

            /* must not call fetchnext() if previous call retrieved last rec in chunk, otherwise
             * fetchnext() will blow away chunk, but the *cache will still point to recs 
             * in the cache. */
            while (!stop && (rec = drms_recordset_fetchnext(drmsEnv, rs, &status, &cstat)) != NULL)
            {
               if (cstat == kRecChunking_NewChunk)
               {
                  if (*cache)
                  {
                     hcon_destroy(cache);
                  }

                  *cache = hcon_create(sizeof(DRMS_Record_t *), gCacheKeySize, NULL, NULL, NULL, NULL, 0);

                  /* insert recordset */
                  hcon_insert(*cache, "bobmould", &rs);
               }
               else if (cstat == kRecChunking_LastInChunk || cstat == kRecChunking_LastInRS)
               {
                  /* this record was the last in chunk - stop caching */
                  stop = 1;
               }

               /* insert, using the obs_date keyword value */
               timestr = drms_getkey_string(rec, kObsDateKey, &status);
               if (timestr)
               {
                  hcon_insert(*cache, timestr, &rec);
                  rehydrated = 1;
               }
               else
               {
                  error = 1;
                  break;
               }
            }

            if (cstat == kRecChunking_NoMoreRecs)
            {
               morerecs = 0;
            }

            /* Try to find tbuf again */
            if (rehydrated)
            {
               if ((prec = (DRMS_Record_t **)hcon_lookup(*cache, tbuf)) != NULL)
               {
                  *recout = *prec;
               }
            }
         }
      }
      else if (*cache == NULL)
      {
         *cache = hcon_create(sizeof(DRMS_Record_t *), gCacheKeySize, NULL, NULL, NULL, NULL, 0);

         /* insert recordset */
         hcon_insert(*cache, "bobmould", &rs);
      }
   }

#if DEBUG
   if (!prec)
   {
      fprintf(stderr, "Record for time '%s' not found in series '%s'.\n", tbuf, series);
   }
#endif

   return error;
}

static int ParseSVRecFields(char *recBuf, char **date, double *xVal, double *yVal, double *zVal, 
			    double *vxVal, double *vyVal, double *vzVal)
{
     int error = 0;
     char *token;
     char *line = strdup(recBuf);

     if (line != NULL)
     {
	  token = strtok(line, " ");

	  if (token != NULL)
	  {
	       /* Must convert date to something that drms recognizes (calendar date) */
	       char year[8];
	       char day[8];
	       char hour[8];
	       char minute[8];

	       strncpy(year, token, 4);
	       year[4] = '\0';
	       strncpy(day, &token[4], 3);
	       day[3] = '\0';
	       strncpy(hour, &token[8], 2);
	       hour[2] = '\0';
	       strncpy(minute, &token[10], 2);
	       minute[2] = '\0';

	       char timeBuf[64];
	       snprintf(timeBuf, sizeof(timeBuf), "%s.01.%s_%s:%s_UT", year, day, hour, minute);

	       *date = strdup(timeBuf);
	  
	       if ((token = strtok(NULL, " ")) != NULL)
	       {
		    sscanf(token, "%lf", xVal);
	       }
	       else
	       {
		    error = 1;
	       }

	       if (error == 0)
	       {
		    if ((token = strtok(NULL, " ")) != NULL)
		    {
			 sscanf(token, "%lf", yVal);
		    }
		    else
		    {
			 error = 1;
		    }
	       }
	       
	       if (error == 0)
	       {
		    if((token = strtok(NULL, " ")) != NULL)
		    {
			 sscanf(token, "%lf", zVal);
		    }
		    else
		    {
			 error = 1;
		    }
	       }

	       if (error == 0)
	       {
		    if ((token = strtok(NULL, " ")) != NULL)
		    {
			 sscanf(token, "%lf", vxVal);
		    }
		    else
		    {
			 error = 1;
		    }
	       }
	       
	       if (error == 0)
	       {
		    if ((token = strtok(NULL, " ")) != NULL)
		    {
			 sscanf(token, "%lf", vyVal);
		    }
		    else
		    {
			 error = 1;
		    }
	       }
	       
	       if (error == 0)
	       {
		    if ((token = strtok(NULL, " ")) != NULL)
		    {
			 sscanf(token, "%lf", vzVal);
		    }
		    else
		    {
			 error = 1;
		    }
	       }
	  }
     }
     else
     {
	  error = 1;
     }

     return error;
}

static int IsDiff(double v1, double v2)
{
   if (!drms_ismissing_double(v1) && !drms_ismissing_double(v2))
   {
      return (fabs(v1 - v2) > 1.0e-11 * (fabs(v1) + fabs(v2)));
   }
   else if (!drms_ismissing_double(v1) || !drms_ismissing_double(v2))
   {
      return 1;
   }

   return 0;
}

static int ExtractStateVectors(DRMS_Env_t *drmsEnv, 
			       char *filePathHELIO, 
			       char *idHELIO,
			       char *filePathGEO, 
			       char *idGEO,
			       char *outSeries)
{
   int stat = DRMS_SUCCESS;
   int error = 0;
   char *obsDate = NULL;
   double xValHel;
   double yValHel;
   double zValHel;
   double vxValHel;
   double vyValHel;
   double vzValHel;
   double xValGeo;
   double yValGeo;
   double zValGeo;
   double vxValGeo;
   double vyValGeo;
   double vzValGeo;

   int addedRecsHELIO = 0;
   int addedRecsGEO = 0;

   HContainer_t *existRCache = NULL;
   DRMS_Record_t *existRec = NULL;
   VectorCache_t *helioOutVCache = NULL;
   VectorCache_t *outVCache = NULL;

   outVCache = CreateVectorCache();

#if DEBUG
   Hash_Table_t hash;
   hash_init(&hash, 49999, 0, (int (*)(const void *, const void *))strcmp, hash_universal_hash);
#endif

#if DEBUG
   int throttle = 0;
#endif

   /* Read heliocentric data from filePathHELIO one line at a time */
   FILE *datafp = NULL;

   if (filePathHELIO);
   {
      datafp = fopen(filePathHELIO, "r");
      if (datafp == NULL)
      {
         error = 1;
         fprintf(stderr, "Could not open %s for reading\n", filePathHELIO);
      }
   }

   if (datafp)
   {
      char lineBuf[kSVLineMax];
      int oneMore = -1;

      CreateOutSeries(drmsEnv, outSeries, &stat);

      if (filePathGEO)
      {
         /* If we are ingesting both geo and helio data, then we need to 
          * pair up the two into records (pair according to prime key 
          * value). Put these heliocentric data into a cache. 
          * and then when the geocentric data are available
          * fetch from the cache and match them up. */
         //snprintf(sbox, sizeof(sbox), "%s_sbox", outSeries);
         //CreateOutSeries(drmsEnv, sbox, &stat);
         helioOutVCache = CreateVectorCache();
      }

      /* Parsing HELIO data. */
#if DEBUG
      TIMER_t *timer = CreateTimer();
      int nitems = 0;
#endif

      while (!error && fgets(lineBuf, kSVLineMax, datafp) != NULL)
      {
	 if (oneMore == -1)
	 { 
	    if (strlen(lineBuf) >= 4)
	    {
	       if (strncmp(lineBuf, "Time", 4) == 0)
	       {
		  oneMore = 1;
	       }
	    }

	    continue;
	 }

	 if (oneMore > 0)
	 {
	    oneMore--;
	    continue;
	 }

	 ParseSVRecFields(lineBuf, 
                          &obsDate, 
                          &xValHel, 
                          &yValHel, 
                          &zValHel, 
                          &vxValHel, 
                          &vyValHel, 
                          &vzValHel);

         TIME od = sscan_time(obsDate);
         char tbuf[128];
         sprint_time(tbuf, od, "UTC", 0);

	 if (helioOutVCache)
	 {
#if DEBUG
            if (nitems < 10000)
            {
#endif 

#if DEBUG
               /* try just hashing, and don't worry about leaking the key */
               char *buff = malloc(64);
               snprintf(buff, 64, "test%d", nitems);
               TIMER_t *timer10 = CreateTimer();
               hash_insert(&hash, buff, (void *)(buff));
               hash_lookup(&hash, buff);
               fprintf(stdout, "hash insert + lookup %f seconds.\n", GetElapsedTime(timer10));
               DestroyTimer(&timer10);
#else
               VCacheCache(helioOutVCache, 
                           tbuf, 
                           &xValHel, 
                           &yValHel, 
                           &zValHel, 
                           &vxValHel, 
                           &vyValHel, 
                           &vzValHel,
                           idHELIO,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL);
#endif
             
#if DEBUG  
            }
#endif
	 }
	 else
	 {
            /* Before creating a new output record, check to see if data have changed. There 
             * exist HELIO input data only. */
           
            /* Use FetchCachedRecord() directly on the orbit-vector series -
             * this avoids multiple requests to psql; must request times in
             * ascending order to use FetchCachedRecord(). */

            /* The first time this call is made, a db cursor is created. There is
             * no need to fetch records older than the time string in tbuf (which
             * corresponds to the time - double - in od), provided that tbuf
             * is monotonically increasing in this loop. Provide a filter that 
             * limits the number of DRMS records cached to speed things up. */
            char filtered[DRMS_MAXQUERYLEN];
            snprintf(filtered, sizeof(filtered), "%s[? %s >= $(%s) ?]", outSeries, kObsDateKey, tbuf);

            if (FetchCachedRecord(drmsEnv, filtered, tbuf, &existRCache, &existRec))
            {
               double xValHelSav = drms_getkey_double(existRec, kSVKeyXHELIO, &stat);
               double yValHelSav = drms_getkey_double(existRec, kSVKeyYHELIO, &stat);
               double zValHelSav = drms_getkey_double(existRec, kSVKeyZHELIO, &stat);
               double vxValHelSav = drms_getkey_double(existRec, kSVKeyVxHELIO, &stat);
               double vyValHelSav = drms_getkey_double(existRec, kSVKeyVyHELIO, &stat);
               double vzValHelSav = drms_getkey_double(existRec, kSVKeyVzHELIO, &stat);

               char *idh = drms_getkey_string(existRec, kSVKey_idHELIO, &stat);
               
               double xValGeoSav = drms_getkey_double(existRec, kSVKeyXGEO, &stat);
               double yValGeoSav = drms_getkey_double(existRec, kSVKeyYGEO, &stat);
               double zValGeoSav = drms_getkey_double(existRec, kSVKeyZGEO, &stat);
               double vxValGeoSav = drms_getkey_double(existRec, kSVKeyVxGEO, &stat);
               double vyValGeoSav = drms_getkey_double(existRec, kSVKeyVyGEO, &stat);
               double vzValGeoSav = drms_getkey_double(existRec, kSVKeyVzGEO, &stat);

               char *id = drms_getkey_string(existRec, kSVKey_idGEO, &stat);

               if (IsDiff(xValHelSav, xValHel) || IsDiff(yValHelSav, yValHel) ||  
                   IsDiff(zValHelSav, zValHel) || IsDiff(vxValHelSav, vxValHel) ||
                   IsDiff(vyValHelSav, vyValHel) || IsDiff(vzValHelSav, vzValHel))
               {
                  /* Difference in HELIO values exists - make record with new helio
                   * values, but old geo values. */
                  VCacheCache(outVCache, 
                              tbuf,
                              &xValHel,
                              &yValHel,
                              &zValHel,
                              &vxValHel,
                              &vyValHel,
                              &vzValHel,
                              idHELIO,
                              &xValGeoSav,
                              &yValGeoSav,
                              &zValGeoSav,
                              &vxValGeoSav,
                              &vyValGeoSav,
                              &vzValGeoSav,
                              id);
               }

               if (idh)
               {
                  free(idh);
               }

               if (id)
               {
                  free(id);
               }
            }
            else
            {
               /* No geo input file, and no previously existing record in outSeries - 
                * add helio data. */
               VCacheCache(outVCache, 
                           tbuf,
                           &xValHel,
                           &yValHel,
                           &zValHel,
                           &vxValHel,
                           &vyValHel,
                           &vzValHel,
                           idHELIO,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL);
            }
	 }

         if (!filePathGEO)
         {
            addedRecsHELIO++;
         }

#if DEBUG
         nitems++;
#endif

#if DEBUG
         throttle++;
#endif

#if DEBUG
	 /* Just to speed things up during debugging */
	 if (throttle == 10)
	 {
	    break;
	 }
#endif

      } /* end while */

#if DEBUG
      fprintf(stdout, "Helio ingest time: %f seconds.\n", GetElapsedTime(timer));
      DestroyTimer(&timer);
#endif

      fclose(datafp);
      datafp = NULL;
   } /* end helio-file processing */

   /* Cannot use existing existRCache since it was populated from times 
    * from the helio data file, and the geo data file may have different
    * times. */
   CloseCachedRecords(&existRCache);
   existRec = NULL;

   /* Geocentric orbit files - save these data in the real series, along with the heliocentric
    * data from sbox series. */
   if (error == 0);
   {
      if (filePathGEO)
      {
         datafp = fopen(filePathGEO, "r");
         if (datafp == NULL)
         {
            error = 1;
            fprintf(stderr, "Could not open %s for reading\n", filePathGEO);
         }
      }
   }

   if (datafp)
   {
      char lineBuf[kSVLineMax];
      int oneMore = -1;
      int helioexist;

#if DEBUG
      TIMER_t *timer = CreateTimer();
#endif 

      while (!error && fgets(lineBuf, kSVLineMax, datafp) != NULL)
      {
	 if (oneMore == -1)
	 { 
	    if (strlen(lineBuf) >= 4)
	    {
	       if (strncmp(lineBuf, "Time", 4) == 0)
	       {
		  oneMore = 1;
	       }
	    }

	    continue;
	 }

	 if (oneMore > 0)
	 {
	    oneMore--;
	    continue;
	 }

	 ParseSVRecFields(lineBuf, 
                          &obsDate, 
                          &xValGeo, 
                          &yValGeo, 
                          &zValGeo, 
                          &vxValGeo, 
                          &vyValGeo, 
                          &vzValGeo);

         /* Find corresponding heliocentric data, if they exist */
         helioexist = 0;
         char *idHELIOtmp = NULL;
         VectorNode_t *hrec = NULL;

         TIME od = sscan_time(obsDate);
         char tbuf[128];
         sprint_time(tbuf, od, "UTC", 0);

         if (helioOutVCache)
         {
            hrec = VCacheLookup(helioOutVCache, tbuf);
         }

         /* Get heliocentric data, if they are also being ingested */
         if (hrec)
         {
            xValHel = hrec->hcix;
            yValHel = hrec->hciy;
            zValHel = hrec->hciz;
            vxValHel = hrec->hcivx;
            vyValHel = hrec->hcivy;
            vzValHel = hrec->hcivz;
            idHELIOtmp = hrec->hciID;
            helioexist = 1;
         }

         /* Before creating a new output record, check to see if data have changed. 
          * There is GEO input, and there is HELIO input if helioexist == 1. */
 
         // DRMS_RecordSet_t *rssav = drms_open_records(drmsEnv, query, &stat);
         /* Use new record-chunking way. */
         //DRMS_RecordSet_t *rssav = drms_open_recordset(drmsEnv, query, &stat);
         //FetchCachedRecord(drmsEnv, outSeries, tbuf, &orbitcache, &savrec);


         /* The first time this call is made, a db cursor is created. There is
          * no need to fetch records older than the time string in tbuf (which
          * corresponds to the time - double - in od), provided that tbuf
          * is monotonically increasing in this loop. Provide a filter that 
          * limits the number of DRMS records cached to speed things up. */
         char filtered[DRMS_MAXQUERYLEN];
         snprintf(filtered, sizeof(filtered), "%s[? %s >= $(%s) ?]", outSeries, kObsDateKey, tbuf);

         if (FetchCachedRecord(drmsEnv, filtered, tbuf, &existRCache, &existRec))
         {
            int heliodiff = 0;
            int geodiff = 0;

            double xValHelSav = drms_getkey_double(existRec, kSVKeyXHELIO, &stat);
            double yValHelSav = drms_getkey_double(existRec, kSVKeyYHELIO, &stat);
            double zValHelSav = drms_getkey_double(existRec, kSVKeyZHELIO, &stat);
            double vxValHelSav = drms_getkey_double(existRec, kSVKeyVxHELIO, &stat);
            double vyValHelSav = drms_getkey_double(existRec, kSVKeyVyHELIO, &stat);
            double vzValHelSav = drms_getkey_double(existRec, kSVKeyVzHELIO, &stat);

            char *id = drms_getkey_string(existRec, kSVKey_idHELIO, &stat);
            
            double xValGeoSav = drms_getkey_double(existRec, kSVKeyXGEO, &stat);
            double yValGeoSav = drms_getkey_double(existRec, kSVKeyYGEO, &stat);
            double zValGeoSav = drms_getkey_double(existRec, kSVKeyZGEO, &stat);
            double vxValGeoSav = drms_getkey_double(existRec, kSVKeyVxGEO, &stat);
            double vyValGeoSav = drms_getkey_double(existRec, kSVKeyVyGEO, &stat);
            double vzValGeoSav = drms_getkey_double(existRec, kSVKeyVzGEO, &stat);

            char *idg = drms_getkey_string(existRec, kSVKey_idGEO, &stat);

            /* Some of these HELIO and GEO data could be missing. */
            geodiff = (IsDiff(xValGeoSav, xValGeo) || IsDiff(yValGeoSav, yValGeo) ||  
                       IsDiff(zValGeoSav, zValGeo) || IsDiff(vxValGeoSav, vxValGeo) ||
                       IsDiff(vyValGeoSav, vyValGeo) || IsDiff(vzValGeoSav, vzValGeo));            

            if (helioexist)
            {
               /* compare existRec to content of filePathHELIO */
               heliodiff = (IsDiff(xValHelSav, xValHel) || IsDiff(yValHelSav, yValHel) ||  
                            IsDiff(zValHelSav, zValHel) || IsDiff(vxValHelSav, vxValHel) ||
                            IsDiff(vyValHelSav, vyValHel) || IsDiff(vzValHelSav, vzValHel));
            }

            if (geodiff)
            {
               /* Difference in values exists */
               if (!helioexist || !heliodiff)
               {
                  /* Difference in one or more GEO values, and no HELIO input file */
                  VCacheCache(outVCache, 
                              tbuf,
                              &xValHelSav,
                              &yValHelSav,
                              &zValHelSav,
                              &vxValHelSav,
                              &vyValHelSav,
                              &vzValHelSav,
                              id,
                              &xValGeo,
                              &yValGeo,
                              &zValGeo,
                              &vxValGeo,
                              &vyValGeo,
                              &vzValGeo,
                              idGEO);
               }
               else
               {
                  /* There is a helio input file, and one or more values in that file 
                   * differ from the values in outSeries. */
                  VCacheCache(outVCache, 
                              tbuf,
                              &xValHel,
                              &yValHel,
                              &zValHel,
                              &vxValHel,
                              &vyValHel,
                              &vzValHel,
                              idHELIOtmp,
                              &xValGeo,
                              &yValGeo,
                              &zValGeo,
                              &vxValGeo,
                              &vyValGeo,
                              &vzValGeo,
                              idGEO);
               }              
            }
            else
            {
               /* The geo file's value are identical to the ones in outSeries - 
                * save a record if the helio ones differ. */
               if (helioexist && heliodiff)
               {
                  VCacheCache(outVCache, 
                              tbuf,
                              &xValHel,
                              &yValHel,
                              &zValHel,
                              &vxValHel,
                              &vyValHel,
                              &vzValHel,
                              idHELIOtmp,
                              &xValGeoSav,
                              &yValGeoSav,
                              &zValGeoSav,
                              &vxValGeoSav,
                              &vyValGeoSav,
                              &vzValGeoSav,
                              idg);
               }
            }
            
            if (id)
            {
               free(id);
            }

            if (idg)
            {
               free(idg);
            }
         }
         else
         {
            /* No previously existing record in outSeries */
            if (helioexist)
            {
               VCacheCache(outVCache, 
                           tbuf,
                           &xValHel,
                           &yValHel,
                           &zValHel,
                           &vxValHel,
                           &vyValHel,
                           &vzValHel,
                           idHELIOtmp,
                           &xValGeo,
                           &yValGeo,
                           &zValGeo,
                           &vxValGeo,
                           &vyValGeo,
                           &vzValGeo,
                           idGEO);
            }
            else
            {
               VCacheCache(outVCache, 
                           tbuf,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           &xValGeo,
                           &yValGeo,
                           &zValGeo,
                           &vxValGeo,
                           &vyValGeo,
                           &vzValGeo,
                           idGEO);
            }
         }

         addedRecsGEO++;

         if (helioexist)
         {
            addedRecsHELIO++;
         }

#if DEBUG
	 throttle--;
	 if (throttle == 0)
	 {
	    break;
	 }
#endif        
      } /* while */

#if DEBUG
      fprintf(stdout, "Geo ingest time: %f seconds.\n", GetElapsedTime(timer));
      DestroyTimer(&timer);
#endif

      fclose(datafp);
      datafp = NULL;
   } /* end geo-file processing */

   CloseCachedRecords(&existRCache);

   /* Now, write all cached vectors into the output series. */
   HIterator_t *cachehit = hiter_create((HContainer_t *)outVCache);
   DRMS_RecordSet_t *rsout = drms_create_records(drmsEnv, 
                                                 outVCache->num_total, 
                                                 outSeries, 
                                                 DRMS_PERMANENT, 
                                                 &stat);
   DRMS_Record_t *recout = NULL;
   int irec = 0;
   VectorNode_t *node = NULL;

   if (rsout && rsout->n > 0)
   {
      while ((node = (VectorNode_t *)hiter_getnext(cachehit)) != NULL)
      {
         recout = rsout->records[irec];

         drms_setkey_string(recout, kSVKeyPrimary, node->tobs);
         drms_setkey_double(recout, kSVKeyXHELIO, node->hcix);
         drms_setkey_double(recout, kSVKeyYHELIO, node->hciy);
         drms_setkey_double(recout, kSVKeyZHELIO, node->hciz);
         drms_setkey_double(recout, kSVKeyVxHELIO, node->hcivx);
         drms_setkey_double(recout, kSVKeyVyHELIO, node->hcivy);
         drms_setkey_double(recout, kSVKeyVzHELIO, node->hcivz);
         drms_setkey_string(recout, kSVKey_idHELIO, node->hciID);

         drms_setkey_double(recout, kSVKeyXGEO, node->gcix);
         drms_setkey_double(recout, kSVKeyYGEO, node->gciy);
         drms_setkey_double(recout, kSVKeyZGEO, node->gciz);
         drms_setkey_double(recout, kSVKeyVxGEO, node->gcivx);
         drms_setkey_double(recout, kSVKeyVyGEO, node->gcivy);
         drms_setkey_double(recout, kSVKeyVzGEO, node->gcivz);
         drms_setkey_string(recout, kSVKey_idGEO, node->gciID);

         irec++;
      }

      error = (drms_close_records(rsout, DRMS_INSERT_RECORD) != DRMS_SUCCESS);
   }
   else
   {
      fprintf(stderr, "Failed to create records in series '%s'.\n", outSeries);
      error = 1;
   }

   DestroyVectorCache(&outVCache);
   DestroyVectorCache(&helioOutVCache);

   if (!error)
   {
      if (filePathHELIO && filePathGEO)
      {
	 fprintf(stdout, 
		 "Ingested '%s' (%d recs) and '%s' (%d recs).\n", 
		 filePathHELIO, 
		 addedRecsHELIO,
		 filePathGEO,
		 addedRecsGEO);
      }
      else if (filePathHELIO)
      {
	 fprintf(stdout, "Ingested '%s' (%d recs).\n", filePathHELIO, addedRecsHELIO);
      }
      else
      {
	 fprintf(stdout, "Ingested '%s' (%d recs).\n", filePathGEO, addedRecsGEO);
      }
   }

   return error;
}

int DoIt(void)
{
   int error = 0;
   int stat = DRMS_SUCCESS;

   if (drms_env == NULL)
   {
      error = 1;
   }
   else
   {
      char seriesin[DRMS_MAXSERIESNAMELEN];
      char seriesout[DRMS_MAXSERIESNAMELEN];
      char serieshist[DRMS_MAXSERIESNAMELEN];
      char *fdsTRange = NULL;

      const char *ns = cmdparams_get_str(&cmdparams, "ns", NULL);
      const char *si = cmdparams_get_str(&cmdparams, "seriesIn", NULL);
      const char *so = cmdparams_get_str(&cmdparams, "seriesOut", NULL);
      const char *timeRange = cmdparams_get_str(&cmdparams, "timeRange", NULL);

      snprintf(seriesin, sizeof(seriesin), "%s.%s", ns, si);
      snprintf(seriesout, sizeof(seriesout), "%s.%s", ns, so);
      snprintf(serieshist, sizeof(serieshist), "%s.%s", ns, kSeriesHistory);

      if (cmdparams_isflagset(&cmdparams, "c"))
      {
	 CreateOutSeries(drms_env, seriesout, &stat);
	 return 0;
      }

      /* Create the record set query */
      if (error == 0)
      {
	 if (strcmp(timeRange, "NOT SPECIFIED") != 0)
	 {
	    size_t len = strlen(kObsDateKey) + strlen(timeRange) + 3;
	    fdsTRange = malloc(sizeof(char) * len + 1);
	    if (fdsTRange != NULL)
	    {
	       snprintf(fdsTRange, len + 1, "[%s=%s]", kObsDateKey, timeRange);
	    }
	    else
	    {
	       error = 1;
	    }
	 }
      }
	  
      if (error == 0)
      {
	 /* Query database to locate the record(s) for a given time range */
	 DRMS_RecordSet_t *rsInput = NULL;
	 DRMS_RecordSet_t *rsHistory = NULL;
	 DRMS_RecordSet_t *rsHELIO = NULL;
	 DRMS_RecordSet_t *rsGEO = NULL;
	 char queryInput[DRMS_MAXQUERYLEN];
	 char queryHistory[DRMS_MAXQUERYLEN];
	 char *filePathHELIO = NULL;
	 char *filePathGEO = NULL;
	 char queryCORE[DRMS_MAXQUERYLEN];
	 char queryHELIO[DRMS_MAXQUERYLEN];
	 char queryGEO[DRMS_MAXQUERYLEN];
	 char hashbuf[kMaxHashKey];
	 char path[DRMS_MAXPATHLEN];

	 char *obs = NULL;
	 char *prodcomp = NULL;
	 char *format = NULL;
	 int filevers;

	 int history = 0;

	 HContainer_t *proclist = NULL;

	 /* Need to get the date of the last helio record ingested */
	 if (drms_series_exists(drms_env, serieshist, &stat))
	 {
	    /* kSeriesHistory keeps track of all orbit data already ingested into 
	     * seriesout.  So, open all records in seriesin, iterate through each
	     * record and skip any records that are already in kSeriesHistory. */

            /* Need this history series because not all files in seriesin will be 
             * ingested.  If a file is an updated version of a data file,
             * but the data don't change, then we add a record in serieshist, but
             * not in seriesout. */

	    history = 1;
	 }

	 if (fdsTRange)
	 {
	    snprintf(queryInput, 
		     sizeof(queryInput), 
		     "%s[%s=%s]%s,%s[%s=%s]%s", 
		     seriesin, 
		     kFdsProductCompKey,
		     kFdsProductCompHELIO,
		     fdsTRange,
		     seriesin, 
		     kFdsProductCompKey,
		     kFdsProductCompGEO,
		     fdsTRange);
	 }
	 else
	 {
	    snprintf(queryInput, 
		     sizeof(queryInput), 
		     "%s[%s=%s],%s[%s=%s]", 
		     seriesin, 
		     kFdsProductCompKey,
		     kFdsProductCompHELIO,
		     seriesin, 
		     kFdsProductCompKey,
		     kFdsProductCompGEO);
	 }

         /* If no files to ingest, this will remain empty */
         proclist = hcon_create(kMaxHashKey, 
                                kMaxHashKey,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                0);

	 rsInput = drms_open_records(drms_env, queryInput, &stat);
	 if (rsInput && rsInput->n > 0)
	 {
	    DRMS_Record_t *recin = NULL;
	    int irec;

	    for (irec = 0; irec < rsInput->n && !error; irec++)
	    {
	       recin = rsInput->records[irec];

	       if (recin)
	       {
		  obs = drms_getkey_string(recin, kObsDateKey, &stat);
		  prodcomp = drms_getkey_string(recin, kFdsProductCompKey, &stat);
		  format = drms_getkey_string(recin, kFdsDataFormatKey, &stat);
		  filevers = drms_getkey_int(recin, kFileVersionKey, &stat);

		  /* Skip if this record has already been ingested */
		  if (history)
		  {
		     snprintf(queryHistory, 
			      sizeof(queryHistory), 
			      "%s[%s=%s][%s=%s][%s=%s][%s=%d]",
			      serieshist,
                              kObsDateKey,
                              obs,
                              kFdsProductCompKey,
			      prodcomp,
                              kFdsDataFormatKey,
			      format,
                              kFileVersionKey,
			      filevers);
			      
		     rsHistory = drms_open_records(drms_env, queryHistory, &stat);
		     if (rsHistory && rsHistory->n > 0)
		     {
			/* skip - already ingested */
			drms_close_records(rsHistory, DRMS_FREE_RECORD);
			continue;			
		     }
                     else
                     {
                        drms_close_records(rsHistory, DRMS_FREE_RECORD);
                     }
		  }

		  /* Skip if this record has been processed in this session */
		  snprintf(hashbuf,
			   sizeof(hashbuf), 
			   "%s[%s=%s][%s=%s][%s=%d]", 
			   seriesin,
			   kObsDateKey,
			   obs,
			   kFdsDataFormatKey,
			   format,
			   kFileVersionKey,
			   filevers);

		  if (hcon_member(proclist, hashbuf))
		  {
		     continue;
		  }

                  /* record is okay to process (and possible ingest) - but first collect
                   * HELIO and GEO in pairs */
                  snprintf(queryCORE, 
                           sizeof(queryCORE), 
                           "%s[%s=%s][%s=%s][%s=%d]",
                           seriesin,
                           kObsDateKey,
                           obs,
                           kFdsDataFormatKey,
                           format,
                           kFileVersionKey,
                           filevers);

		  hcon_insert(proclist, queryCORE, queryCORE);
	       } 

	       if (obs)
	       {
		   free(obs);
	       }

	       if (prodcomp)
	       {
		   free(prodcomp);
	       }

	       if (format)
	       {
		   free(format);
	       }

	    } /* record loop */

	    drms_close_records(rsInput, DRMS_FREE_RECORD);
	 }

	 /* Now, open all HELIO and GEO files that match the core queries saved
	  * in proclist */
	 HIterator_t *hit = hiter_create(proclist);
	 char *iquery = NULL;
	   
	 while ((iquery = (char *)hiter_getnext(hit)) != NULL && !error)
	 {
#if 0
            /* If iquery starts with "###", then either the helio or geo (or both) 
             * record does not need to be ingested - the data have already been 
             * ingested from a file of an earlier version (but the data have not changed). */
            if (strncmp(iquery, "###", 3) == 0)
            {
               /* If there is another identical query, except that it does NOT start
                * with "###", then this means either the helio or geo data DID change 
                * so we have to reingest the data (and it may have been ingested in
                * this session already). */
               if (hcon_member(proclist, iquery[3]))
               {
                  /* The other matching record had data that did change.  When iquery 
                   * points to that record, data for both helio and geo will be ingested.
                   * Ignore this query. */
                  continue;
               }
               else
               {
                  /* There was a file with an updated version that had data that does not
                   * differ from data that has already been ingested. But we have to 
                   * mark in the history that these data have been processed. */
                  oktoadd = 0;
                  iquery = iquery[3];
               }
            }
#endif
	    snprintf(queryHELIO, 
		     sizeof(queryHELIO), 
		     "%s[%s=%s]",
		     iquery,
		     kFdsProductCompKey,
		     kFdsProductCompHELIO);

	    snprintf(queryGEO, 
		     sizeof(queryGEO), 
		     "%s[%s=%s]",
		     iquery,
		     kFdsProductCompKey,
		     kFdsProductCompGEO);

            DRMS_Record_t *recsv = NULL;
            TIME obsvalHELIO = 0;
            char *dprodvalHELIO = NULL;
            char *prodcompvalHELIO = NULL;
            char *formatvalHELIO = NULL;
            int fileversvalHELIO = 0;
            TIME obsvalGEO = 0;
            char *dprodvalGEO = NULL;
            char *prodcompvalGEO = NULL;
            char *formatvalGEO = NULL;
            int fileversvalGEO = 0;

	    rsHELIO = drms_open_records(drms_env, queryHELIO, &stat);

	    if (rsHELIO && stat != DRMS_ERROR_UNKNOWNSERIES)
	    {
	       if (rsHELIO->n != 1)
	       {
		  fprintf(stderr, "Expected only one record with query '%s'\n", queryHELIO);
		  error = 1;
	       }
	       else
	       {
		  /* Get the HELIO filepath */
		  DRMS_Record_t *recHELIO = rsHELIO->records[0];
			   
		  DRMS_Segment_t *helioseg = drms_segment_lookup(recHELIO, kSVFileSegName);

		  if (helioseg != NULL)
		  {
		     drms_record_directory(recHELIO, path, 1);

                     if (*path)
                     {
                        size_t len = strlen(path) + strlen(helioseg->filename) + 1;
                        filePathHELIO = malloc(sizeof(char) * len + 1);
                        if (filePathHELIO)
                        {
                           snprintf(filePathHELIO, len + 1, "%s/%s", path, helioseg->filename);
                        }
                     }
		  }

                  if (filePathHELIO)
                  {
                     recsv = rsHELIO->records[0];
                     obsvalHELIO = drms_getkey_time(recsv, kObsDateKey, &stat);
                     dprodvalHELIO = drms_getkey_string(recsv, kFdsDataProductKey, &stat);
                     prodcompvalHELIO = drms_getkey_string(recsv, kFdsProductCompKey, &stat);
                     formatvalHELIO = drms_getkey_string(recsv, kFdsDataFormatKey, &stat);
                     fileversvalHELIO = drms_getkey_int(recsv, kFileVersionKey, &stat);
                  }
	       }
	    }

            if (!error)
            {
               rsGEO = drms_open_records(drms_env, queryGEO, &stat);
            }

	    if (rsGEO && stat != DRMS_ERROR_UNKNOWNSERIES)
	    {
	       if (rsGEO->n != 1)
	       {
		  fprintf(stderr, "Expected only one record with query '%s'\n", queryGEO);
		  error = 1;
	       }
	       else
	       {
		  /* Get the GEO filepath */
		  DRMS_Record_t *recGEO = rsGEO->records[0];
			   
		  DRMS_Segment_t *geoseg = drms_segment_lookup(recGEO, kSVFileSegName);
			   
		  if (geoseg != NULL)
		  {
		     drms_record_directory(recGEO, path, 1);

                     if (*path)
                     {
                        size_t len = strlen(path) + strlen(geoseg->filename) + 1;
                        filePathGEO = malloc(sizeof(char) * len + 1);
                        if (filePathGEO)
                        {
                           snprintf(filePathGEO, len + 1, "%s/%s", path, geoseg->filename);
                        }
                     }
		  }

                  if (filePathGEO)
                  {
                     recsv = rsGEO->records[0];
                     obsvalGEO = drms_getkey_time(recsv, kObsDateKey, &stat);
                     dprodvalGEO = drms_getkey_string(recsv, kFdsDataProductKey, &stat);
                     prodcompvalGEO = drms_getkey_string(recsv, kFdsProductCompKey, &stat);
                     formatvalGEO = drms_getkey_string(recsv, kFdsDataFormatKey, &stat);
                     fileversvalGEO = drms_getkey_int(recsv, kFileVersionKey, &stat);
                  }
	       }
	    }

            if (!error)
            {
               /* Create HELIO and GEO ids */
               char idHELIO[256];
               char idGEO[256];

               if (rsHELIO && filePathHELIO)
               {
                  char *datestr = drms_getkey_string(rsHELIO->records[0], 
                                                     kObsDateKey, 
                                                     &stat);

                  /* Date must come first - the indexing code forces this */
                  snprintf(idHELIO, 
                           sizeof(idHELIO), 
                           "%s[%s][%s][%s][%d]", 
                           rsHELIO->records[0]->seriesinfo->seriesname,
                           datestr,
                           prodcompvalHELIO,
                           formatvalHELIO,
                           fileversvalHELIO);

                  if (datestr)
                  {
                     free(datestr);
                  }
               }

               if (rsGEO && filePathGEO)
               {
                  char *datestr = drms_getkey_string(rsGEO->records[0], 
                                                     kObsDateKey, 
                                                     &stat);

                  /* Date must come first - the indexing code forces this */
                  snprintf(idGEO, 
                           sizeof(idGEO), 
                           "%s[%s][%s][%s][%d]", 
                           rsGEO->records[0]->seriesinfo->seriesname,
                           datestr,
                           prodcompvalGEO,
                           formatvalGEO,
                           fileversvalGEO);

                  if (datestr)
                  {
                     free(datestr);
                  }
               }

               if (filePathHELIO || filePathGEO)
               {
                  error = ExtractStateVectors(drms_env, 
                                              filePathHELIO, 
                                              idHELIO,
                                              filePathGEO,
                                              idGEO,
                                              seriesout);
               }
            }
	    
	    if (!error)
	    {
              /* Add a record to the history series so that this record doesn't get
               * reingested - need to create 2 records, one for helio and one for
               * geo */
              if (!history)
              {
                 /* Create history series. */
                 CreateHistSeries(drms_env, serieshist, &stat);
              }
              
              int nhist = 0;
              int ihist = 0;
              
              if (filePathHELIO)
              {
                 nhist++;
              }

              if (filePathGEO)
              {
                 nhist++;
              }

              if (nhist > 0)
              {
                 DRMS_RecordSet_t *rshist = drms_create_records(drms_env, 
                                                                nhist, 
                                                                serieshist, 
                                                                DRMS_PERMANENT,
                                                                &stat);
                 DRMS_Record_t *rechist = NULL;

                 if (rshist)
                 {
                    if (filePathHELIO)
                    {
                       rechist = rshist->records[ihist++];
                       drms_setkey_time(rechist, kObsDateKey, obsvalHELIO);
                       drms_setkey_string(rechist, kFdsDataProductKey, dprodvalHELIO);
                       drms_setkey_string(rechist, kFdsProductCompKey, prodcompvalHELIO);
                       drms_setkey_string(rechist, kFdsDataFormatKey, formatvalHELIO);
                       drms_setkey_int(rechist, kFileVersionKey, fileversvalHELIO);
                    }

                    if (filePathGEO)
                    {
                       rechist = rshist->records[ihist++];
                       drms_setkey_time(rechist, kObsDateKey, obsvalGEO);
                       drms_setkey_string(rechist, kFdsDataProductKey, dprodvalGEO);
                       drms_setkey_string(rechist, kFdsProductCompKey, prodcompvalGEO);
                       drms_setkey_string(rechist, kFdsDataFormatKey, formatvalGEO);
                       drms_setkey_int(rechist, kFileVersionKey, fileversvalGEO);
                    }

                    drms_close_records(rshist, DRMS_INSERT_RECORD);
                 }
              }
	    }

	    if (rsHELIO)
	    {
	       drms_close_records(rsHELIO, DRMS_FREE_RECORD);
               rsHELIO = NULL;
	    }

	    if (rsGEO)
	    {
	       drms_close_records(rsGEO, DRMS_FREE_RECORD);
               rsGEO = NULL;
	    }

            if (dprodvalHELIO)
            {
               free(dprodvalHELIO);
            }

            if (prodcompvalHELIO)
            {
               free(prodcompvalHELIO);
            }

            if (formatvalHELIO)
            {
               free(formatvalHELIO);
            }

            if (dprodvalGEO)
            {
               free(dprodvalGEO);
            }

            if (prodcompvalGEO)
            {
               free(prodcompvalGEO);
            }

            if (formatvalGEO)
            {
               free(formatvalGEO);
            }

	 } /* iquery - core queries */

	 hiter_destroy(&hit);
	 hcon_destroy(&proclist);
      }
	  
      if (fdsTRange)
      {
	 free(fdsTRange);
      }	  
   }
     
   return error;
}
