/* demo_td08062007.c */

/*  The purpose of this module is to demonstrate how JSOC can input data cubes (like
 *  those output by fastrack) and transform it into a dataseries that can be used
 *  to generate power spectra.
 *
 *  The demonstration was held on August 6-9, 2007 at Stanford University.
 *
 *  Assumptions:
 *    1. One segment per record.
 *
 *
 *  Caveats:
 *    1. Run "limit stacksize unlimited" or else you'll likely get a crash.
 *
 *
 *  --Art Amezcua
 */

#include <complex.h>
#include <sys/time.h>
#include "jsoc_main.h"
#include "drms_types.h"
#include "fftw3.h"

#define PRINTJSD 0
#define DEBUGMSGS 0

char *module_name = "demo_td08062007";

#define kRecSetIn      "recsin"
#define kDSOut         "dsout"
#define kOutSeries     "su_arta.TestDemoTD"

ModuleArgs_t module_args[] =
{
     {ARG_STRING, kRecSetIn, "",          "Input data series."},
     {ARG_STRING, kDSOut,    kOutSeries,  "Output data series."},
     {ARG_END}
};

typedef enum
{
   kDemoError_Success,
   kDemoError_CouldntCheckSeries,
   kDemoError_CouldntCreateSeries,
   kDemoError_InputIncompatibleWithSeries,
   kDemoError_CouldntCreateRecord,
   kDemoError_SegmentNotFound,
   kDemoError_UnsupportedDimensionality,
   kDemoError_UnsupportedDataType,
   kDemoError_BadParameter,
   kDemoError_CouldntCreateArray,
   kDemoError_UnknownSeries,
   kDemoError_BadQuerySyntax,
   kDemoError_DrmsOpen
} DemoError_t;

typedef enum
{
   kOutputSeriesDisp_Unknown,
   kOutputSeriesDisp_DoesntExist,
   kOutputSeriesDisp_Exists,
} OutputSeriesDisp_t;

/* For testing - convert record structure into jsd-file format */

#if PRINTJSD
void drms_keyword_print_jsd(DRMS_Keyword_t *key) {
    printf("Keyword:%s",key->info->name);
    if (key->info->islink) {
      printf(", link, %s, %s, %s\n", key->info->linkname, 
	     key->info->target_key,
	     key->info->description);
    } else {
      printf(", %s", drms_type2str(key->info->type));
      printf(", %s", drms_keyword_recscopestr(key, NULL));
      if (key->info->per_segment) 
	printf(", segment");
      else 
	printf(", record");
      printf(", ");
      if (key->info->type == DRMS_TYPE_STRING) {
	char qf[DRMS_MAXFORMATLEN+2];
	sprintf(qf, "\"%s\"", key->info->format);
	printf(qf, key->value.string_val);
      }
      else 
	drms_keyword_printval(key);      
      if (key->info->unit[0] != ' ') {
	printf(", %s, %s, \"%s\"", key->info->format,
	       key->info->unit,
	       key->info->description);
      } else {
	printf(", %s, none, \"%s\"", key->info->format,
	       key->info->description);
      }
    }
    printf("\n");
}

void drms_segment_print_jsd(DRMS_Segment_t *seg) {
  int i;
  printf("Data: %s, ", seg->info->name);
  if (seg->info->islink) {
    printf("link, %s, %s", seg->info->linkname, seg->info->target_seg);
    if (seg->info->naxis) {
      printf(", %d", seg->info->naxis);
      printf(", %d", seg->axis[0]);
      for (i=1; i<seg->info->naxis; i++) {
	printf(", %d", seg->axis[i]);
      }
    }
  } else {
    switch(seg->info->scope)
      {
      case DRMS_CONSTANT:
	printf("constant");
	break;
      case DRMS_VARIABLE:
	printf("variable");
	break;
      case DRMS_VARDIM:
	printf("vardim");
	break;
      default:
	printf("Illegal value: %d", (int)seg->info->scope);
      }
    printf(", %s, %d", drms_type2str(seg->info->type), seg->info->naxis);
    if (seg->info->naxis) {
      printf(", %d", seg->axis[0]);
      for (i=1; i<seg->info->naxis; i++) {
	printf(", %d", seg->axis[i]);
      }
    }
    printf(", %s, ", seg->info->unit);  
    switch(seg->info->protocol)
      {
      case DRMS_GENERIC:
	printf("generic");
	break;
      case DRMS_BINARY:
	printf("binary");
	break;
      case DRMS_BINZIP:
	printf("binzip");
	break;
      case DRMS_FITZ:
	printf("fitz");
	break;
      case DRMS_FITS:
	printf("fits");
	break;
      case DRMS_MSI:
	printf("msi");
	break;
      case DRMS_TAS:
	printf("tas");
	if (seg->info->naxis) {
	  printf(", %d", seg->blocksize[0]);      
	  for (i=1; i<seg->info->naxis; i++)
	    printf(", %d", seg->blocksize[i]);
	}
	break;
      default:
	printf("Illegal value: %d", (int)seg->info->protocol);
      }
  }
  printf(", \"%s\"\n", seg->info->description);
}

void drms_link_print_jsd(DRMS_Link_t *link) {
  printf("Link: %s, %s, ", link->info->name, link->info->target_series);
  if (link->info->type == STATIC_LINK)
    printf("static");
  else
    printf("dynamic");
  printf(", \"%s\"\n", link->info->description);
}

void print_jsd(DRMS_Record_t *rec) {
  const int fwidth=17;
  int i;
  HIterator_t hit;
  DRMS_Link_t *link;
  DRMS_Keyword_t *key;
  DRMS_Segment_t *seg;

  printf("#=====General Series Information=====\n");
  printf("%-*s\t%s\n",fwidth,"Seriesname:",rec->seriesinfo->seriesname);
  printf("%-*s\t\"%s\"\n",fwidth,"Author:",rec->seriesinfo->author);
  printf("%-*s\t%s\n",fwidth,"Owner:",rec->seriesinfo->owner);
  printf("%-*s\t%d\n",fwidth,"Unitsize:",rec->seriesinfo->unitsize);
  printf("%-*s\t%d\n",fwidth,"Archive:",rec->seriesinfo->archive);
  printf("%-*s\t%d\n",fwidth,"Retention:",rec->seriesinfo->retention);
  printf("%-*s\t%d\n",fwidth,"Tapegroup:",rec->seriesinfo->tapegroup);
  if (rec->seriesinfo->pidx_num) {
    printf("%-*s\t%s",fwidth,"Index:",rec->seriesinfo->pidx_keywords[0]->info->name);
    for (i=1; i<rec->seriesinfo->pidx_num; i++)
      printf(", %s", (rec->seriesinfo->pidx_keywords[i])->info->name);
    printf("\n");
  }
  printf("%-*s\t%s\n",fwidth,"Description:",rec->seriesinfo->description);
  printf("\n#=====Links=====\n");
  hiter_new(&hit, &rec->links); 
  while( (link = (DRMS_Link_t *)hiter_getnext(&hit)) )
    drms_link_print_jsd(link);

  printf("\n#=====Keywords=====\n");
  hiter_new(&hit, &rec->keywords);
  while( (key = (DRMS_Keyword_t *)hiter_getnext(&hit)) )
    drms_keyword_print_jsd(key);

  printf("\n#=====Segments=====\n");
  hiter_new(&hit, &rec->segments);
  while( (seg = (DRMS_Segment_t *)hiter_getnext(&hit)) )
    drms_segment_print_jsd(seg);
}
#endif // PRINTJSD

/* CreateOutSeries
 *
 *   Using an input record as a prototype, create an output series.  Modifies
 *   segment dimensionality.
 *
 *   env - DRMS session object.
 *   rec - Record representing input data.
 */
static int CreateOutSeries(DRMS_Env_t *env, DRMS_Record_t *rec)
{
   int error = kDemoError_CouldntCreateSeries;
   int status = DRMS_SUCCESS;

   DRMS_Record_t *sourceRec = drms_template_record(env, rec->seriesinfo->seriesname, &status);
   DRMS_Segment_t *segproto = NULL;

   if (sourceRec)
   {
      DRMS_Record_t *prototype = drms_create_recproto(sourceRec, &status);

      if (prototype)
      {
	 DRMS_SegmentDimInfo_t di;
	 
	 if (hcon_size(&(prototype->segments)) != 1)
	 {
	    XASSERT(0);
	    fprintf(stderr, 
		    "Warning: more than one segment in input series; using first segment.\n");
	 }

	 if ((segproto = drms_segment_lookupnum(prototype, 0)) != NULL)
	 {
	    if (segproto->info->naxis != 3)
	    {
	       error = kDemoError_UnsupportedDimensionality;
	    }
	    else if (segproto->info->type != DRMS_TYPE_FLOAT && 
		     segproto->info->type != DRMS_TYPE_DOUBLE)
	    {
	       error = kDemoError_UnsupportedDataType;
	    }
	    else 
	    {
	       if (segproto->info->protocol == DRMS_DSDS ||
		   segproto->info->protocol == DRMS_LOCAL)
	       {
		  /* Input is a DSDS series. */
		  segproto->info->protocol = DRMS_FITS;
	       }
	       else
	       {
		  /* Input is a true DRMS series. */
	       }

	       drms_segment_getdims(segproto, &di);
	       di.naxis = 2;
	       di.axis[0] = (segproto->axis)[1] / 2;
	       di.axis[1] = (segproto->axis)[2] / 2;
	       drms_segment_setdims(segproto, &di);

#if PRINTJSD
	       print_jsd(prototype);
#endif 
	    
	       status = drms_create_series_fromprototype(&prototype, kOutSeries, 0);
	       if (status == DRMS_SUCCESS)
	       { 
		  error = kDemoError_Success;
	       }
	    }
	 }
	 else
	 {
	    drms_destroy_recproto(&prototype);
	 }
      }
   }

   return error;
}

/* CheckCompat
 *
 *   Indicate if a record of input data is compatible with the output series.  Returns 1 
 *   if compatible, 0 otherwise.
 *	 
 *   env   - DRMS session object.
 *   rec   - Record representing data.
 *   dsout - Output series name.
 */
static int CheckCompat(DRMS_Env_t *env, DRMS_Record_t *rec, const char *dsout)
{
   int compat = 0;
   int status = DRMS_SUCCESS;

   DRMS_Record_t *prototype =  drms_create_recproto(rec, &status);
   DRMS_Segment_t *segproto = NULL;
   HContainer_t *matchSegNames = NULL;
   
   if (prototype)
   {
      if ((segproto = drms_segment_lookupnum(prototype, 0)) != NULL)
      {
	 DRMS_SegmentDimInfo_t di;

	 if (segproto->info->protocol == DRMS_DSDS ||
	     segproto->info->protocol == DRMS_LOCAL)
	 {
	    /* Input is a DSDS series. */
	    segproto->info->protocol = DRMS_FITS;
	 }
	 else
	 {
	    /* Input is a true DRMS series. */
	 }

	 /* Data saved has a different format than data input. Modify prototype 
	  * to reflect this. */
	 drms_segment_getdims(segproto, &di);
	 di.naxis = 2;
	 di.axis[0] = (segproto->axis)[1] / 2;
	 di.axis[1] = (segproto->axis)[2] / 2;
	 drms_segment_setdims(segproto, &di);

	 XASSERT((matchSegNames = (HContainer_t *)malloc(sizeof(HContainer_t))) != NULL);
	 compat = drms_series_checkrecordcompat(env,
						dsout, 
						prototype, 
						matchSegNames, 
						&status);
	 
	 hcon_destroy(&matchSegNames);
	 drms_destroy_recproto(&prototype);

	 if (!compat)
	 {
	    fprintf(stderr, 
		    "Output series %s is not compatible with output data.\n", 
		    dsout);
	 }
      }     
   }

   return compat;
}

/* Crunch
 *
 *   Manipluate the 3D input data to produce a 2D power diagram as output.
 *  
 *   arrin  - 3D array of input data.
 *   arrout - 2D array of output data.
 */
DemoError_t Crunch(DRMS_Array_t *arrin, DRMS_Array_t **arrout)
{
   DemoError_t error = kDemoError_Success;
   int status = DRMS_SUCCESS;

   if (arrin && arrout)
   {
      long long arrsize = drms_array_count(arrin);
      DRMS_Type_t type = arrin->type;
      int laxis1 = drms_array_nth_axis(arrin, 0);
      int laxis2 = drms_array_nth_axis(arrin, 1);
      int laxis3 = drms_array_nth_axis(arrin, 2);
      int axisout[2] = {laxis2 / 2, laxis3 / 2};
      long long nElem = drms_array_count(arrin);

      *arrout = drms_array_create(type,
				  2,
				  axisout,
				  NULL,
				  &status);

      if (*arrout)
      {
	 fftwf_plan p;
	 fftwf_complex datc[laxis3][laxis1][laxis1 / 2 + 1];

	 /* Iterate through data values.  Works on double or float. */
	 long long iData;

	 if (type == DRMS_TYPE_FLOAT)
	 {
	    float *pDataIn = arrin->data;

	    for (iData = 0; iData < nElem; iData++)
	    {
	       pDataIn[iData] = pDataIn[iData] / (float)arrsize;
	    }
	 }
	 else if (type == DRMS_TYPE_DOUBLE)
	 {
	    double *pDataIn = arrin->data;

	    for (iData = 0; iData < nElem; iData++)
	    {
	       pDataIn[iData] = pDataIn[iData] / (double)arrsize;
	    }	   
	 }

	 p = fftwf_plan_dft_r2c_3d(laxis3,
				   laxis2,
				   laxis1,
				   arrin->data,
				   &datc[0][0][0],
				   FFTW_ESTIMATE);

	 fftwf_execute(p);
	 fftwf_destroy_plan(p);

	 int i;
	 int j;
	 int k;
	 int ii;

	 int iLim = laxis3 / 2;
	 int jLim = laxis2 / 2;
	 int kLim = laxis1 / 2;

	 memset((*arrout)->data, 0, drms_array_size(*arrout));

	 if (type == DRMS_TYPE_FLOAT)
	 {
	    float *power = (*arrout)->data;

	    for (j = 0; j < jLim; j++)
	    {
	       for (k = 0; k < kLim; k++)
	       {
		  ii = (int)(sqrtf(powf((float)j, 2.0) + powf((float)k, 2.0)));
			
		  if (ii < jLim)
		  {	    
		     for (i = 0; i < iLim; i++)
		     {
			power[i * axisout[0] + ii] +=
			  powf(cabsf(datc[i][j][k]), 2.0);
		     }
		  }
	       }
	    }
	 }
	 else if(type == DRMS_TYPE_DOUBLE)
	 {
	    double *power = (*arrout)->data;

	    for (j = 0; j < jLim; j++)
	    {
	       for (k = 0; k < kLim; k++)
	       {
		  ii = (int)(sqrtf(powf((float)j, 2.0) + powf((float)k, 2.0)));
			
		  if (ii < jLim)
		  {	    
		     for (i = 0; i < iLim; i++)
		     {
			power[i * axisout[0] + ii] +=
			  pow(cabs(datc[i][j][k]), 2.0);
		     }
		  }
	       }
	    }
	 }	   
      } /* arrout */
      else
      {
	 error = kDemoError_CouldntCreateArray;
      }
   }
   else
   {
      error = kDemoError_BadParameter;
   }

   return error;
}

/* LogTime
 *
 *   Print out to stdout a message showing elapsed time since the last time LogTime
 *   was called with the clear flag set.
 *
 *   clear - Flag to establish baseline.  Call with clear == 1 to reset timer.
 *     Then call with clear == 0 to obtain elapsed time since first call.
 *   msg   - Message to print in output.
 */
void LogTime(int clear, const char *msg)
{
   static double stTime = 0.0;
   struct timeval tv;
   double newTime = 0.0;
   double elapsed = 0.0;

   gettimeofday(&tv, NULL);
   newTime = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;

   if (clear)
   {
      stTime = newTime;
      elapsed = 0.0;
   }
   else
   {
      elapsed = newTime - stTime;
      fprintf(stdout, 
	      "%s: done (%.3g seconds elapsed)\n", 
	      msg,
	      elapsed);
      fflush(stdout);
   }
}

/* DoIt
 *
 *   Module entry point.
 */
int DoIt(void) 
{
   int status = DRMS_SUCCESS;
   DemoError_t error = kDemoError_Success;
   char *inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL);
   char rquery[DRMS_MAXQUERYLEN];
   char *dsout = cmdparams_get_str(&cmdparams, kDSOut, NULL);
   DRMS_RecordSet_t *inRecSet = NULL;
   DRMS_RecordSetType_t rqueryType;
   OutputSeriesDisp_t disp;

   rqueryType = drms_record_getquerytype(inRecQuery);

   if (rqueryType == kRecordSetType_DSDS)
   {
      fprintf(stdout, "Fetching records from DSDS");
      fflush(stdout);
      LogTime(1, NULL);
   }

   drms_series_exists(drms_env, dsout, &status);
   disp = (status == DRMS_ERROR_UNKNOWNSERIES) ? kOutputSeriesDisp_DoesntExist : 
     ((status != DRMS_SUCCESS) ? kOutputSeriesDisp_Unknown : kOutputSeriesDisp_Exists);

   if (disp == kOutputSeriesDisp_Unknown)
   {
      error = kDemoError_CouldntCheckSeries;
   }
   else
   {
      if (disp == kOutputSeriesDisp_Exists && 
	  !strchr(inRecQuery, '[') && 
	  rqueryType == kRecordSetType_PlainFile)
      {
	 /* series exists, so use its prime keys when creating records from local
	  * data files (unless user specified prime keys on cmd-line) */
	 int nPKeys = 0;
	 char **pkArray = NULL;
   
	 pkArray = drms_series_createpkeyarray(drms_env, dsout, &nPKeys, &status);
	 if (status == DRMS_SUCCESS)
	 {
	    char *buf = (char *)malloc(sizeof(char) * DRMS_MAXKEYNAMELEN * nPKeys + 32);
	    *buf = '\0';
	    int iKey;

	    for (iKey = 0; iKey < nPKeys; iKey++)
	    {
	       if (!iKey == 0)
	       {
		  strcat(buf, ",");
	       }

	       strcat(buf, pkArray[iKey]);
	    }

	    snprintf(rquery, sizeof(rquery), "%s[%s]", inRecQuery, buf);

	    if (buf)
	    {
	       free(buf);
	    }

	    drms_series_destroypkeyarray(&pkArray, nPKeys);
	 }
	 else
	 {
	    error = kDemoError_CouldntCheckSeries;
	 }
      }
      else
      {
	 /* drms query */
	 snprintf(rquery, sizeof(rquery), "%s", inRecQuery);
      }

      if (!error)
      {
	 inRecSet = drms_open_records(drms_env, rquery, &status);

	 if (rqueryType == kRecordSetType_DSDS)
	 {
	    LogTime(0, "");
	 }

	 if (status == DRMS_ERROR_UNKNOWNSERIES)
	 {
	    error = kDemoError_UnknownSeries;
	    fprintf(stderr, "Unknown series.\n");
	 }
	 else if (status == DRMS_ERROR_INVALIDDATA)
	 {
	    error = kDemoError_BadQuerySyntax;
	 }
	 else if (status == DRMS_INEXACT)
	 {
	    error = kDemoError_BadQuerySyntax;
	    fprintf(stderr, "Error in drms_open_records(), err num %d.\n"
		    "Probably bad series query syntax, or invalid local file/directory name.\n", 
		    status);
	 }
	 else if (status != DRMS_SUCCESS)
	 {
	    error = kDemoError_DrmsOpen;
	    fprintf(stderr, "Error in drms_open_records(), err num %d.\n"
		    "Probably bad series query syntax, or invalid local file/directory name.\n", 
		    status);
	 }
      }
   }

   if (!error)
   {
      int nRecs = inRecSet->n;
      int iRec;

      /* create output series rec prototype */
      DRMS_Record_t *recout = NULL;
      DRMS_Segment_t *segout = NULL;
      DRMS_Segment_t *segin = NULL;

      if (nRecs > 0)
      {
	 if (disp == kOutputSeriesDisp_DoesntExist)
	 {
#if DEBUGMSGS
	    fprintf(stdout, "env cache before\n\n");
	    hcon_print(&(drms_env->series_cache));
#endif

	    if (CreateOutSeries(drms_env, inRecSet->records[0]))
	    {
	       error = kDemoError_CouldntCreateSeries;
	    }

#if DEBUGMSGS
	    fprintf(stdout, "env cache after\n\n");
	    hcon_print(&(drms_env->series_cache));
#endif	    
	 }
      }

      for (iRec = 0; !error && iRec < nRecs; iRec++)
      {
	 DRMS_Record_t *rec = inRecSet->records[iRec];
	 if (rec)
	 {
	    fprintf(stdout, "Record %d\n", iRec);
	    fflush(stdout);
	    
	    /* Check for compatibility of each record with existing series. */
	    if (!CheckCompat(drms_env, rec, dsout))
	    {
	       error = kDemoError_InputIncompatibleWithSeries;
	    }
	    
	    if (!error)
	    {
	       segin = drms_segment_lookupnum(rec, 0);
	       if (segin)
	       {
		  DRMS_Array_t *arrin = NULL;
		  DRMS_Array_t *arrout = NULL;

		  fprintf(stdout, "  reading data");
		  fflush(stdout);
		  LogTime(1, NULL);

		  arrin = drms_segment_read(segin, segin->info->type, &status);

		  LogTime(0, "");

		  if (arrin)
		  {
		     int naxis = drms_array_naxis(arrin);
		     DRMS_Type_t type = arrin->type;

		     if (naxis != 3)
		     {
			error = kDemoError_UnsupportedDimensionality;
		     }
		     else if (type != DRMS_TYPE_FLOAT && type != DRMS_TYPE_DOUBLE)
		     {
			error = kDemoError_UnsupportedDataType;
		     }
		   
		     if (error == kDemoError_Success)
		     {
			DRMS_RecordSet_t *rs = NULL;

			fprintf(stdout, "  processing data");
			fflush(stdout);
			LogTime(1, NULL);

			error = Crunch(arrin, &arrout);

			LogTime(0, "");

			if (error == kDemoError_Success)
			{
			   rs = drms_create_records(drms_env, 
						    1, 
						    dsout, 
						    DRMS_PERMANENT, 
						    &status);
			   
			   if (status != DRMS_SUCCESS)
			   {
			      error = kDemoError_CouldntCreateRecord;
			   }
			}

			if (!error)
			{
			   recout = rs->records[0];
			}

			if (recout)
			{
			   DRMS_Keyword_t *key = NULL;
			   segout = drms_segment_lookupnum(recout, 0);
			   if (segout)
			   {
			      HIterator_t *hit = hiter_create(&(rec->keywords));
			      if (hit)
			      {
				 while ((key = hiter_getnext(hit)) != NULL)
				 {
				    status = drms_setkey(recout, 
							 key->info->name,
							 key->info->type,
							 &(key->value));
				 }

				 hiter_destroy(&hit);
			      }

			      drms_segment_write(segout, arrout, 0);
			   }

			   drms_close_records(rs, DRMS_INSERT_RECORD);
			}

			drms_free_array(arrin);
		     }
		  } /* arrin */
	       }
	       else
	       {
		  error = kDemoError_SegmentNotFound;
	       }
	    } /* !error */
	 } /* rec */
      } /* iRec */

      if (inRecSet)
      {
	 drms_close_records(inRecSet, DRMS_FREE_RECORD);
      }
   }

   return error;
}
