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

#define DEBUG 1

#define kDefaultNamespace "sdo"
#define kDefaultSeriesIn "moc_fds"
#define kDefaultSeriesOut "fds_orbit_vectors"
#define kSeriesHistory "fds_orbit_ingesthist"
#define kObsDateKey "OBS_DATE"
#define kFdsDataProductKey "FDS_DATA_PRODUCT"
#define kFdsProductCompKey "FDS_PRODUCT_COMP"
#define kFdsDataFormatKey "DATA_FORMAT"
#define kFileVersionKey "FILE_VERSION"
#define kFdsDataProductHELIO "predHelioOrb"
#define kFdsDataProductGEO "predGeoOrb"
#define kFdsProductCompHELIO "ORBIT_HELIO"
#define kFdsProductCompGEO "ORBIT_GEO"

#define kSVFileSegName "FILENAME"

#define kSVLineMax 1024
#define kMaxHashKey 1024

/* State-vector data keys */
#define kSVKeyPrimary "OBS_DATE"
#define kSVKeyXHELIO "X_HELIO"
#define kSVKeyYHELIO "Y_HELIO"
#define kSVKeyZHELIO "Z_HELIO"
#define kSVKeyVxHELIO "Vx_HELIO"
#define kSVKeyVyHELIO "Vy_HELIO"
#define kSVKeyVzHELIO "Vz_HELIO"
#define kSVKeyXGEO "X_GEO"
#define kSVKeyYGEO "Y_GEO"
#define kSVKeyZGEO "Z_GEO"
#define kSVKeyVxGEO "Vx_GEO"
#define kSVKeyVyGEO "Vy_GEO"
#define kSVKeyVzGEO "Vz_GEO"

#define kSVKey_idHELIO "id_HELIO"
#define kSVKey_idGEO "id_GEO"

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
      tKey->info->isdrmsprime = isprime;
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

	    snprintf(prototype->seriesinfo->description,
		     DRMS_MAXCOMMENTLEN,
		     "%s",
		     "Helio- and Geo-centric orbit position and velocity vectors.  Prime key is observation date.");

	    /* slotted keyword isdrmsprime, but not really prime */
	    akey = AddKey(prototype, DRMS_TYPE_TIME, kSVKeyPrimary, "0", "UTC", kRecScopeType_TS_EQ, "slotted observation date", 1);

	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "index");
	    pkey = AddKey(prototype, kIndexKWType, keyname, kIndexKWFormat, "none", kRecScopeType_Index, "index for OBS_DATE", 1);

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

/* Caller must clean up records by calling drms_close_records() */
static int FetchRecords(DRMS_Env_t *drmsEnv, char *query, DRMS_RecordSet_t **recordSet, int *nRecs)
{
     int error = 0;

     *recordSet = drms_open_records(drmsEnv, query, &error);

     if (error != 0 || recordSet == NULL) 
     {
	  printf("drms_open_records failed, recSetQuery=%s, error=%d.  Aborting.\n", query, error);
     }
     else
     {
	  *nRecs = (*recordSet)->n;
     }

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

static int SetSVKey(DRMS_Record_t *rec, char *keyName, DRMS_Type_Value_t value)
{
   int error = 0;

   if (rec != NULL && keyName != NULL)
   {
      error = drms_setkey(rec, keyName, DRMS_TYPE_DOUBLE, &value);
      if (error != 0)
      {
	 printf("Failed to set key '%s' to record\n", keyName);
      }
   }
   else
   {
      error = 1;
      printf("Invalid record or keyword name.\n");
   }

   return error;
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
#if DEBUG
   int throttle = 0;
#endif

   char sbox[DRMS_MAXSERIESNAMELEN];

   /* Read heliocentric data from filePathHELIO one line at a time */
   FILE *datafp = NULL;

   if (error == 0 && filePathHELIO);
   {
      datafp = fopen(filePathHELIO, "r");
   }

   if (datafp != NULL)
   {
      char lineBuf[kSVLineMax];
      int oneMore = -1;

      CreateOutSeries(drmsEnv, outSeries, &stat);

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

	 ParseSVRecFields(lineBuf, &obsDate, &xValHel, &yValHel, &zValHel, &vxValHel, &vyValHel, &vzValHel);

	 DRMS_RecordSet_t *recSet = NULL;
	 if (filePathGEO)
	 {
	    /* Put these heliocentric data into temporary records in a sandbox series. */
	    snprintf(sbox, sizeof(sbox), "%s_sbox", outSeries);

	    CreateOutSeries(drmsEnv, sbox, &stat);

	    /* Since these are DRMS_TRANSIENT, they should never get saved between runs of 
	     * this module. */
	    recSet = drms_create_records(drms_env, 1, sbox, DRMS_TRANSIENT, &error);
	 }
	 else
	 {
	    recSet = drms_create_records(drms_env, 1, outSeries, DRMS_PERMANENT, &error);
	 }

	 if (recSet != NULL && error == 0)
	 {
	    DRMS_Record_t *rec = recSet->records[0];
	    int nPrime = rec->seriesinfo->pidx_num;
		    
	    if (nPrime != 1)
	    {
	       error = 1;
	       printf("The specified dataseries %s has the incorrect number of primary keys\n", outSeries);
	    }
	    else
	    {
	       DRMS_Type_Value_t value;
			      
	       /* Add primary key */
	       value.time_val = sscan_time(obsDate);
	       error = drms_setkey(rec, kSVKeyPrimary, DRMS_TYPE_TIME, &value);
	       if (error != 0)
	       {
		  printf("Failed to set key '%s' to record\n", kSVKeyPrimary);
	       }
			      
	       /* Add state vector keys */
	       value.double_val = xValHel;
	       SetSVKey(rec, kSVKeyXHELIO, value);
			      
	       value.double_val = yValHel;
	       SetSVKey(rec, kSVKeyYHELIO, value);
			      
	       value.double_val = zValHel;
	       SetSVKey(rec, kSVKeyZHELIO, value);
			      
	       value.double_val = vxValHel;
	       SetSVKey(rec, kSVKeyVxHELIO, value);
			      
	       value.double_val = vyValHel;
	       SetSVKey(rec, kSVKeyVyHELIO, value);
			      
	       value.double_val = vzValHel;
	       SetSVKey(rec, kSVKeyVzHELIO, value);

	       drms_setkey_string(rec, kSVKey_idHELIO, idHELIO);
	    }    
	 }
	 else
	 {
	    error = 1;
	    printf("The specified output series %s does not exist\n", outSeries);
	 }

	 if (error == 0)
	 {
	    error = drms_close_records(recSet, DRMS_INSERT_RECORD);
	    if (error == 0)
	    {

	       if (!filePathGEO)
	       {
		  addedRecsHELIO++;
	       }
#if DEBUG
	       throttle++;
#endif
	    }
	 }
	 else
	 {
	    error = drms_close_records(recSet, DRMS_FREE_RECORD);
	    break;
	 }
#if DEBUG
	 /* Just to speed things up during debugging */
	 if (throttle == 10)
	 {
	    break;
	 }
#endif       
      } /* end while */

      fclose(datafp);
   }
   else
   {
      /* Couldn't open file */
      error = 1;
      printf("Could not open %s for reading\n", filePathHELIO);
   }

   /* Geocentric orbit files - save these data in the real series, along with the heliocentric
    * data from sbox series. */
   if (error == 0 && filePathGEO);
   {
      datafp = fopen(filePathGEO, "r");
   }

   if (datafp != NULL)
   {
      char lineBuf[kSVLineMax];
      int oneMore = -1;
      DRMS_RecordSet_t *rs = NULL;
      char rsquery[DRMS_MAXQUERYLEN];
      int nrecs;
      DRMS_Type_t type;
      DRMS_Type_Value_t value;

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

	 ParseSVRecFields(lineBuf, &obsDate, &xValGeo, &yValGeo, &zValGeo, &vxValGeo, &vyValGeo, &vzValGeo);

	 DRMS_RecordSet_t *recSet = drms_create_records(drms_env, 1, outSeries, DRMS_PERMANENT, &error);
	 char *idHELIOtmp = NULL;

	 if (recSet != NULL && error == 0)
	 {
	    DRMS_Record_t *outrec = recSet->records[0];
	    DRMS_Record_t *inrec = NULL;
	    char tbuf[128];
			      
	    /* Set primary key */
	    value.time_val = sscan_time(obsDate);
	    error = drms_setkey(outrec, kSVKeyPrimary, DRMS_TYPE_TIME, &value);
	    if (error != 0)
	    {
	       printf("Failed to set key '%s' to record\n", kSVKeyPrimary);
	    }
	    
	    /* Find corresponding heliocentric data, if it exists */
	    TIME od = sscan_time(obsDate);
	    sprint_time(tbuf, od, "UTC", 0);
	    snprintf(rsquery, sizeof(rsquery), "%s[%s]", sbox, tbuf);
	    FetchRecords(drmsEnv, rsquery, &rs, &nrecs);
	    if (rs && rs->n == 1)
	    {
	       inrec = rs->records[0]; /* heliocentric */

	       /* Get heliocentric data */
	       if (inrec)
	       {
		  value = drms_getkey(inrec, kSVKeyXHELIO, &type, &stat);
		  SetSVKey(outrec, kSVKeyXHELIO, value);

		  value = drms_getkey(inrec, kSVKeyYHELIO, &type, &stat);
		  SetSVKey(outrec, kSVKeyYHELIO, value);

		  value = drms_getkey(inrec, kSVKeyZHELIO, &type, &stat);
		  SetSVKey(outrec, kSVKeyZHELIO, value);

		  value = drms_getkey(inrec, kSVKeyVxHELIO, &type, &stat);
		  SetSVKey(outrec, kSVKeyVxHELIO, value);

		  value = drms_getkey(inrec, kSVKeyVyHELIO, &type, &stat);
		  SetSVKey(outrec, kSVKeyVyHELIO, value);
	       
		  value = drms_getkey(inrec, kSVKeyVzHELIO, &type, &stat);
		  SetSVKey(outrec, kSVKeyVzHELIO, value);
		  
		  idHELIOtmp = drms_getkey_string(inrec, kSVKey_idHELIO, &stat);
		  if (idHELIOtmp)
		  {
		     drms_setkey_string(outrec, kSVKey_idHELIO, idHELIOtmp);
		     free(idHELIOtmp);
		     idHELIOtmp = NULL;
		  }
	       }
	      
	       drms_close_records(rs, DRMS_FREE_RECORD);
	       
	       addedRecsHELIO++;
	    }

	    /* Add geocentric state vector keys */
	    value.double_val = xValGeo;
	    SetSVKey(outrec, kSVKeyXGEO, value);
			      
	    value.double_val = yValGeo;
	    SetSVKey(outrec, kSVKeyYGEO, value);
			      
	    value.double_val = zValGeo;
	    SetSVKey(outrec, kSVKeyZGEO, value);
			      
	    value.double_val = vxValGeo;
	    SetSVKey(outrec, kSVKeyVxGEO, value);
			      
	    value.double_val = vyValGeo;
	    SetSVKey(outrec, kSVKeyVyGEO, value);
			      
	    value.double_val = vzValGeo;
	    SetSVKey(outrec, kSVKeyVzGEO, value);

	    drms_setkey_string(outrec, kSVKey_idGEO, idGEO);

	    error = (drms_close_records(recSet, DRMS_INSERT_RECORD) != DRMS_SUCCESS);

	    if (!error)
	    {
	       addedRecsGEO++;
	    }
	 }
	 
#if DEBUG
	 throttle--;
	 if (throttle == 0)
	 {
	    break;
	 }
#endif        
      } /* while */
   }
   else
   {
      /* Couldn't open file */
      error = 1;
      printf("Could not open %s for reading\n", filePathGEO);
   }

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

      char *ns = cmdparams_get_str(&cmdparams, "ns", NULL);
      char *si = cmdparams_get_str(&cmdparams, "seriesIn", NULL);
      char *so = cmdparams_get_str(&cmdparams, "seriesOut", NULL);
      char *timeRange = cmdparams_get_str(&cmdparams, "timeRange", NULL);

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
	 char *dataprod = NULL;
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

	    history = 1;
	 }

	 if (fdsTRange)
	 {
	    snprintf(queryInput, 
		     sizeof(queryInput), 
		     "%s[%s=%s]%s,%s[%s=%s]%s", 
		     seriesin, 
		     kFdsDataProductKey,
		     kFdsDataProductHELIO,
		     fdsTRange,
		     seriesin, 
		     kFdsDataProductKey,
		     kFdsDataProductGEO,
		     fdsTRange);
	 }
	 else
	 {
	    snprintf(queryInput, 
		     sizeof(queryInput), 
		     "%s[%s=%s],%s[%s=%s]", 
		     seriesin, 
		     kFdsDataProductKey,
		     kFdsDataProductHELIO,
		     seriesin, 
		     kFdsDataProductKey,
		     kFdsDataProductGEO);
	 }

	 rsInput = drms_open_records(drms_env, queryInput, &stat);
	 if (rsInput && rsInput->n > 0)
	 {
	    DRMS_Record_t *recin = NULL;
	    int irec;

	    proclist = hcon_create(kMaxHashKey, 
				   kMaxHashKey,
				   NULL,
				   NULL,
				   NULL,
				   NULL,
				   0);

	    for (irec = 0; irec < rsInput->n && !error; irec++)
	    {
	       recin = rsInput->records[irec];

	       if (recin)
	       {
		  obs = drms_getkey_string(recin, kObsDateKey, &stat);
		  dataprod = drms_getkey_string(recin, kFdsDataProductKey, &stat);
		  prodcomp = drms_getkey_string(recin, kFdsProductCompKey, &stat);
		  format = drms_getkey_string(recin, kFdsDataFormatKey, &stat);
		  filevers = drms_getkey_int(recin, kFileVersionKey, &stat);

		  /* Skip if this record has already been ingested */
		  if (history)
		  {
		     snprintf(queryHistory, 
			      sizeof(queryHistory), 
			      "%s[%s=%s][%s=%s][%s=%s][%s=%s][%s=%d]",
			      serieshist,
			      kObsDateKey,
			      obs,
			      kFdsDataProductKey,
			      dataprod,
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

		  /* record is okay to ingest - but first collect HELIO and GEO in pairs */
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
	    } /* record loop */

	    drms_close_records(rsInput, DRMS_FREE_RECORD);
	 }

	 /* Now, open all HELIO and GEO files that match the core queries saved
	  * in proclist */
	 HIterator_t *hit = hiter_create(proclist);
	 char *iquery = NULL;
	   
	 while ((iquery = (char *)hiter_getnext(hit)) != NULL)
	 {
	    snprintf(queryHELIO, 
		     sizeof(queryHELIO), 
		     "%s[%s=%s][%s=%s]",
		     iquery,
		     kFdsDataProductKey,
		     kFdsDataProductHELIO,
		     kFdsProductCompKey,
		     kFdsProductCompHELIO);

	    snprintf(queryGEO, 
		     sizeof(queryGEO), 
		     "%s[%s=%s][%s=%s]",
		     iquery,
		     kFdsDataProductKey,
		     kFdsDataProductGEO,
		     kFdsProductCompKey,
		     kFdsProductCompGEO);

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
		     size_t len = strlen(path) + strlen(helioseg->filename) + 1;
		     filePathHELIO = malloc(sizeof(char) * len + 1);
		     if (filePathHELIO)
		     {
			snprintf(filePathHELIO, len + 1, "%s/%s", path, helioseg->filename);
		     }
		  }
	       }

	       drms_close_records(rsHELIO, DRMS_FREE_RECORD);
	    }

	    rsGEO = drms_open_records(drms_env, queryGEO, &stat);

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
		     size_t len = strlen(path) + strlen(geoseg->filename) + 1;
		     filePathGEO = malloc(sizeof(char) * len + 1);
		     if (filePathGEO)
		     {
			snprintf(filePathGEO, len + 1, "%s/%s", path, geoseg->filename);
		     }
		  }
	       }

	       drms_close_records(rsGEO, DRMS_FREE_RECORD);
	    }
		  
	    error = ExtractStateVectors(drms_env, 
					filePathHELIO, 
					queryHELIO,
					filePathGEO,
					queryGEO,
					seriesout);
	 }

	 hiter_destroy(&hit);
	 hcon_destroy(&proclist);

	 if (obs)
	 {
	    free(obs);
	 }

	 if (dataprod)
	 {
	    free(dataprod);
	 }

	 if (prodcomp)
	 {
	    free(prodcomp);
	 }

	 if (format)
	 {
	    free(format);
	 }
      }
	  
      if (fdsTRange)
      {
	 free(fdsTRange);
      }	  
   }
     
   return error;
}
