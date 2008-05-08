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

#define kDefaultNamespace "sdo"
#define kDefaultSeriesIn "moc_fds"
#define kDefaultSeriesOut "fds_orbit_vectors"
#define kDefaultSeriesHistory "sdo.fds_orbit_ingesthist"
#define kObsTimeKey "OBS_DATE"
#define kFdsDataProductKeyValHel "[FDS_DATA_PRODUCT=predHelioOrb]"
#define kFdsDataProductKeyValGeo "[FDS_DATA_PRODUCT=predGeoOrb]"
#define kSVFileSegName "FILENAME"

#define kSVLineMax 1024

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

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "ns", kDefaultNamespace, "working namespace (sdo, sdo_ground, sdo_dev)"},
  {ARG_STRING, "seriesIn", kDefaultSeriesIn, "name of series containing FDS data"},
  {ARG_STRING, "timeRange", "NOT SPECIFIED", "time ranges separated by commas"},
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
	    akey->value.time_val = sscan_time("1993.01.01_12:00:00_UT");

	    /* step */
	    snprintf(keyname, sizeof(keyname), "%s_%s", kSVKeyPrimary, "step");
	    akey = AddKey(prototype, DRMS_TYPE_STRING, keyname, "%s", "none", kRecScopeType_Constant, "time step for OBS_DATE", 0);
	    copy_string(&(akey->value.string_val), "1d");

	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyXHELIO, "%f", "km", kRecScopeType_Variable, "X position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyYHELIO, "%f", "km", kRecScopeType_Variable, "Y position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyZHELIO, "%f", "km", kRecScopeType_Variable, "Z position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVxHELIO, "%f", "km/sec", kRecScopeType_Variable, "velocity along X", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVyHELIO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Y", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVzHELIO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Z", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyXGEO, "%f", "km", kRecScopeType_Variable, "X position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyYGEO, "%f", "km", kRecScopeType_Variable, "Y position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyZGEO, "%f", "km", kRecScopeType_Variable, "Z position", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVxGEO, "%f", "km/sec", kRecScopeType_Variable, "velocity along X", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVyGEO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Y", 0);
	    AddKey(prototype, DRMS_TYPE_DOUBLE, kSVKeyVzGEO, "%f", "km/sec", kRecScopeType_Variable, "velocity along Z", 0);

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

static void FreeRecords(DRMS_RecordSet_t *recSet)
{
     drms_close_records(recSet, DRMS_FREE_RECORD);
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

static int ParseDuration(char *durationStr, double *duration)
{
     int error = 0;

     char *end;
     char *p = durationStr;
     double dval = (int)strtod(p, &end);

     if ((dval == 0 && end==p) || 
	 ((dval == HUGE_VALF || dval == -HUGE_VALF) && errno == ERANGE))
     {
	  fprintf(stderr,"Syntax Error: Expected finite floating point value at "
		  "beginning of time duration, found '%s'.\n", p);
	  error = 1;
     }
     else
     {
	  p = end;
     }

     if (error == 0)
     {
	  switch(*p++)
	  {
	  case 's':
	    *duration = 1.0 * dval;
	    break;
	  case 'm':
	    *duration = 60.0 * dval;
	    break;
	  case 'h':
	    *duration = 3600.0 * dval;
	    break;
	  case 'd':
	    *duration = 86400.0 * dval;
	    break;
	  default:
	    fprintf(stderr,"Syntax Error: Time duration unit must be one of 's', "
		    "'m', 'h', 'd', found '%c', '%s'.\n", *(p-1), p);
	    error = 1;
	  }
     }

     return error;
}

/* If <timeQuery> is within the range specified by <rs>, *bAns is set to 1, otherwise it is set to 0 */
/* The data types in <rs> must be of type DRMS_TYPE_TIME */
/* Adapted from sql_primekey_value_set() */
static int TimeInRange(char *time, ValueRangeSet_t *rs, int *bAns)
{
     int error = 0;
     int bInRange = 0;

     /* Convert time string into drmsTime */
     TIME timeQuery = sscan_time(time) - sscan_time("1977.01.01_00:00_TAI");

     do 
     {
	  if (rs->type == SINGLE_VALUE)
	  {
	       bInRange = (rs->start.time_val == timeQuery) ? 1 : 0;
	       if (bInRange == 1)
	       {
		    break;
	       }
	  }
	  else
	  {
	       if (rs->type == START_END)
	       {
		    /* For the interval case, rs->x is the end of the interval */
		    if (timeQuery >= rs->start.time_val && timeQuery <= rs->x.time_val)
		    {
			 bInRange = 1;
			 break;
		    }
	       }
	       else if (rs->type == START_DURATION)
	       {
		    /* For the duration case, rs->x is the increment from rs->start */
		    TIME endTime = rs->start.time_val + rs->x.time_val;
		    if (timeQuery >= rs->start.time_val && timeQuery < endTime)
		    {
			 bInRange = 1;
			 break;
		    }
	       }

	       if (rs->has_skip)
	       {
		    fprintf(stderr,"(has_skip == 1) not yet supported in value queries.\n");
		    error = 1;
		    break;
	       }
	  }
	  
	  rs = rs->next;
     } while (rs);     
     
     if (error == 0)
     {
	  *bAns = bInRange;
     }
     else
     {
	  *bAns = 0;
     }
     
     return error;
}

static int CreateTimeRange(char *timeRange, ValueRangeSet_t **rs)
{
     int error = 0;

     if (timeRange == NULL)
     {
	  error = 1;
     }
     else
     {
	  int n;
	  char *p = timeRange;
	  ValueRangeSet_t *vr = NULL;
	  ValueRangeSet_t *head = NULL;
	  
	  do 
	  {
	       if (vr)
	       {
		    XASSERT(vr->next = (ValueRangeSet_t *)malloc(sizeof(ValueRangeSet_t)));
		    vr = vr->next;
		    memset(vr, 0, sizeof(ValueRangeSet_t));
	       }
	       else 
	       {
		    XASSERT(vr = (ValueRangeSet_t *)malloc(sizeof(ValueRangeSet_t)));
		    head = vr;
		    memset(vr,0,sizeof( ValueRangeSet_t));
	       }
	       
	       /* Get start */
	       if ((n = drms_sscanf(p, DRMS_TYPE_TIME, &vr->start)) == 0)    
	       {
		    fprintf(stderr,"Syntax Error: Expected start value of type %s in"
			    " value range, found '%s'.\n", drms_type2str(DRMS_TYPE_TIME), p);
		    error = 1;
		    break;
	       }
	       else
	       {
		    p += n;     
	       }
	       
	       if (*p == '-')
	       {
		    vr->type = START_END;
	       }
	       else if (*p == '/')
	       {
		    vr->type = START_DURATION;
	       }
	       else
	       {
		    vr->type = SINGLE_VALUE;
	       }
	       
	       /* Get end or duration "x" */
	       if (vr->type != SINGLE_VALUE)
	       {
		    ++p;
		    
		    if (vr->type == START_END)
		    {
			 if ((n = drms_sscanf(p, DRMS_TYPE_TIME, &vr->x)) == 0)    
			 {
			      fprintf(stderr,"Syntax Error: Expected end value of"
				      " type %s in value range, found '%s'.\n", drms_type2str(DRMS_TYPE_TIME), p);
			      error = 1;
			      break;
			 }
			 else
			 {
			      p += n;
			 }
		    }
		    else
		    {
			 if (ParseDuration(p, &vr->x.time_val))
			 {
			      fprintf(stderr,"Syntax Error: Expected time duration "
				      " in value range, found '%s'.\n", p);
			      error = 1;
			      break;
			 }	  
		    }
		    /* Get skip */
		    if (*p == '@')
		    {
			 ++p;
			 vr->has_skip = 1;
			 if (ParseDuration(p, &vr->skip.time_val))    
			 {
			      fprintf(stderr,"Syntax Error: Expected skip (time duration)"
				      " in value range, found '%s'.\n", p);
			      error = 1;
			      break;
			 }
		    }
		    else
		    {
			 vr->has_skip = 0;
		    }
	       }
	  } while (*p++ == ',');
	  
	  
	  if (error == 0)
	  {
	       *rs = head;
	  }
	  else
	  {
	       *rs = NULL;
	  }
     }

     return error;
}

static void DestroyTimeRange(ValueRangeSet_t *head)
{
     /* Clean up range structures */
     ValueRangeSet_t *vr = NULL;
     
     while (head)
     { 
	  vr = head->next; 
	  free(head); 
	  head = vr; 
     }
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
			       char *filePathGEO, 
			       char *timeRange, 
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

   int skippedRecs = 0;
   int addedRecs = 0;

   ValueRangeSet_t *tRangeSet = NULL;

   char sbox[DRMS_MAXSERIESNAMELEN];

   if (timeRange != NULL)
   {
      printf("Time range is %s.\n", timeRange);
      error = CreateTimeRange(timeRange, &tRangeSet);
   }

   /* Read heliocentric data from filePathHELIO one line at a time */
   FILE *datafp = NULL;

   if (error == 0);
   {
      datafp = fopen(filePathHELIO, "r");
   }

   if (datafp != NULL)
   {
      char lineBuf[kSVLineMax];
      int oneMore = -1;

      CreateOutSeries(drmsEnv, outSeries, &stat);

      while (fgets(lineBuf, kSVLineMax, datafp) != NULL)
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
	       
	 /* Strip out records that do not fall within the range of times specified in fdsTRange */
	 int bAns = 0;
	 if (TimeInRange(obsDate, tRangeSet, &bAns) != 0)
	 {
	    error = 1;
	 }
	 else
	 {
	    if (bAns == 0)
	    {
	       skippedRecs++;
	       continue;
	    }
	 }

	 /* Put these heliocentric data into temporary records in a sandbox series. */
	 snprintf(sbox, sizeof(sbox), "%s_sbox", outSeries);

	 CreateOutSeries(drmsEnv, sbox, &stat);

	 /* Since these are DRMS_TRANSIENT, they should never get saved between runs of 
	  * this module. */
	 DRMS_RecordSet_t *recSet = drms_create_records(drms_env, 1, sbox, DRMS_TRANSIENT, &error);
	  
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
	       // printf("Added obsDate=%s, x=%f, y=%f, z=%f, vx=%f, vy=%f, vz=%f\n", obsDate, xValHel, yValHel, zValHel, vxValHel, vyValHel, vzValHel);
	       addedRecs++;
	    }
	 }
	 else
	 {
	    error = drms_close_records(recSet, DRMS_FREE_RECORD);
	    break;
	 }
	       
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
   if (error == 0);
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
      DRMS_Record_t *rec = NULL;
      DRMS_Type_t type;
      DRMS_Type_Value_t value;

      while (fgets(lineBuf, kSVLineMax, datafp) != NULL)
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
	       }
	      
	       drms_close_records(rs, DRMS_FREE_RECORD);
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

	    drms_close_records(recSet, DRMS_INSERT_RECORD);
	 }
      }
   }
   else
   {
      /* Couldn't open file */
      error = 1;
      printf("Could not open %s for reading\n", filePathGEO);
   }

   if (tRangeSet)
   {
      DestroyTimeRange(tRangeSet);
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
      char *fdsSeries = NULL;
      char *fdsTRange = NULL;
      char *svSeries = NULL;
	  
      char *seriesIn = cmdparams_get_str(&cmdparams, "seriesIn", NULL);
      char *timeRange = cmdparams_get_str(&cmdparams, "timeRange", NULL);
      char *seriesOut = cmdparams_get_str(&cmdparams, "seriesOut", NULL);

      if (cmdparams_isflagset(&cmdparams, "c"))
      {
	 CreateOutSeries(drms_env, seriesOut, &stat);
	 return 0;
      }

      fdsSeries = strdup(seriesIn);
      svSeries = strdup(seriesOut);
	  
      XASSERT(fdsSeries && svSeries);
	  
      if (!fdsSeries || !svSeries)
      {
	 error = 1;
      }

      /* Create the record set query */
      if (error == 0)
      {
	 if (strcmp(timeRange, "NOT SPECIFIED") != 0)
	 {
	    size_t len = strlen(kObsTimeKey) + strlen(timeRange) + 3;
	    fdsTRange = malloc(sizeof(char) * len + 1);
	    if (fdsTRange != NULL)
	    {
	       snprintf(fdsTRange, len + 1, "[%s=%s]", kObsTimeKey, timeRange);
	    }
	    else
	    {
	       error = 1;
	    }
	 }
      }
	  
      if (error == 0)
      {
	 char *recSetQuery = NULL;
	 size_t len = strlen(fdsSeries) + strlen(kFdsDataProductKeyValHel);
	       
	 if (fdsTRange != NULL)
	 {
	    len += strlen(fdsTRange);
	    recSetQuery = malloc(sizeof(char) * len + 1);
	    if (recSetQuery != NULL)
	    {
	       snprintf(recSetQuery, len + 1, "%s%s%s", fdsSeries, fdsTRange, kFdsDataProductKeyValHel);
	    }
	    else
	    {
	       error = 1;
	    }
	 }
	 else
	 {
	    recSetQuery = malloc(sizeof(char) * len + 1);
	    snprintf(recSetQuery, len + 1, "%s%s", fdsSeries, kFdsDataProductKeyValHel);
	    if (recSetQuery == NULL)
	    {
	       error = 1;
	    }
	 }

	 printf("Heliocentric record set query is %s\n", recSetQuery);

	 if (recSetQuery == NULL)
	 {
	    error = 1;
	 }
	       
	 if (error == 0)
	 {
	    /* Query database to locate the record(s) for a given time range */
	    int nRecs = 0;
	    DRMS_RecordSet_t *recordSet = NULL;
	    DRMS_Record_t *rec = NULL;
	    int irec;
		
	    error = FetchRecords(drms_env, recSetQuery, &recordSet, &nRecs);
	    free(recSetQuery);    

	    for (irec = 0; error == 0 && irec < nRecs; irec++)
	    {   	     
	       printf("Number of records fetched is %d\n", nRecs);
	      
	       rec = recordSet->records[irec];

	       if (rec)
	       {
		  char path[DRMS_MAXPATHLEN];
				   
		  DRMS_Segment_t *helioseg = drms_segment_lookup(rec, kSVFileSegName);

		  if (helioseg != NULL)
		  {
		     drms_record_directory(rec, path, 1);
		     len = strlen(path) + strlen(helioseg->filename) + 1;
		     char *filePathHel = malloc(sizeof(char) * len + 1);
		     if (filePathHel)
		     {
			snprintf(filePathHel, len + 1, "%s/%s", path, helioseg->filename);
			char *tRange = (strcmp(timeRange, "NOT SPECIFIED") == 0) ? NULL : timeRange;


			/* Geocentric vectors */
			char rsqueryGEO[DRMS_MAXQUERYLEN];
			DRMS_RecordSet_t *rsGEO = NULL;
			int nrecsGEO = 0;
			char *filePathGeo = NULL;
			DRMS_Segment_t *geoseg = NULL;

			snprintf(rsqueryGEO, 
				 sizeof(rsqueryGEO),
				 "%s%s%s", 
				 fdsSeries, 
				 fdsTRange, 
				 kFdsDataProductKeyValGeo);

			error = FetchRecords(drms_env, rsqueryGEO, &rsGEO, &nrecsGEO);

			if (nrecsGEO != 1)
			{
			   fprintf(stderr, "Should only be one record returned.\n");
			}

			if (rsGEO)
			{
			   geoseg = drms_segment_lookup(rsGEO->records[0], kSVFileSegName);

			   if (geoseg != NULL)
			   {
			      drms_record_directory(rsGEO->records[0], path, 1);
			      len = strlen(path) + strlen(geoseg->filename) + 1;
			      filePathGeo = malloc(sizeof(char) * len + 1);

			      if (filePathGeo)
			      {
				 snprintf(filePathGeo, len + 1, "%s/%s", path, geoseg->filename);
				 error = ExtractStateVectors(drms_env, filePathHel, filePathGeo, tRange, svSeries);
				 free(filePathGeo);
			      }
			   }

			   FreeRecords(rsGEO);
			}

			free(filePathHel);

		     }
		     else
		     {
			error = 1;
		     }		     
		  }
	       }
	      
	    }

	    FreeRecords(recordSet);
	 }
      }
	  
      if (fdsSeries)
      {
	 free(fdsSeries);
      }
	  
      if (fdsTRange)
      {
	 free(fdsTRange);
      }
	  
      if (svSeries)
      {
	 free(svSeries);
      }
   }
     
   return error;
}
