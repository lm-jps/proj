/* arithTool.c */

/* Mapping from vds/sds to drms:
 *   VDS (a set of datasets)     -> Set of record sets (not implemented in DRMS yet)
 *   Dataset (a set of records)  -> Record set
 *   Record (a set of SDSs)      -> Record
 *   SDS (a FITS file)           -> Segment
 *
 * VDS is implemented as an array (indexed by var number) of datasets
 *
 * <recSetIn> - Record set on which to operate
 * <withRecSet> - For binary operations, the second set of records on which to operate.
 *    The data series containing these two record sets must have equivalent primary keyword sets
 *    and overlapping segment sets.  Otherwise, the two sets are too dissimilar to operate on.
 *    The number of records in this set must be either equal to the number of records in
 *    <recSetIn>, or it must be equal to one record.  If it is the former, then
 *    there is a one-to-one correspondence between the records in <recSetIn> and
 *    <withRecSet>.  The keyword-tuple that selects one record in <recSetIn> should
 *    also select one record in <withRecSet>.  Each record in <recSetIn> is paired with
 *    its corresponding record in <withRecSet>.  The binary operation is then performed
 *    on each pair of records.  If the number of records in <withRecSet> is equal to one,
 *    then each
 *    record in <recSetIn> is paired with the single record in <withRecSet>.  Then
 *    the binary operation is performed on each pair.  The binary operation is
 *    only performed between equivalently named segments.  In other words, if a record
 *    from <recSetIn> contains three segments named SegA, SegB, and SegC, and the corresponding
 *    record from <withRecSet< contains two semgnets named SegA and SegB, the binary operation is
 *    performed between <recSetIn>:SegA and <withRecSet>:SegA and between <recSetIn>:SegB
 *    and <withRecSet>:SegB.
 *
 */

/* Adding a comment for testing. */

#include "jsoc_main.h"

#define kRecSetIn "recSetIn"
#define kWithRecSet "withRecSet"
#define kSeriesOut "seriesOut"
#define kSegList "segments"
#define kOp "op"
#define kBzero "bzero"
#define kBscale "bscale"
#define kDatatype "dtype"
#define kNoSegsFlag "n"
#define kNotSpecified "NOT SPECIFIED"
#define kMaxSegs 1024
#define kMaxQuery 2048
#define kOutSeriesDesc "Output series for arithTool calculation."
#define kInSliceLower "sl"
#define kInSliceUpper "su"
#define kPosOut "pout"

#define TESTER 0

typedef enum
{
     kArithOpUnknown = 0,
     kArithOpMul,
     kArithOpDiv,
     kArithOpAdd,
     kArithOpSub,
     kArithOpBinaryDummy, /* Delimits binary from unary ops */
     kArithOpAbs,
     kArithOpSqrt,
     kArithOpLog,
     kArithOpPow,
     kArithOpSqr,
     kArithOpRecip,
     kArithOpSubmean,
     kArithOpNop,
     kArithOpLastDummy
} ArithOp_t;

ModuleArgs_t module_args[] =
{
  {ARG_FLAG, "h", "0", "help - print usage info"},
  {ARG_STRING, kRecSetIn, "", "A set of records with segments upon which to operate."},
  {ARG_INTS, kInSliceLower, "-1", "Bottom-left corner of slice (-1 means no slicing)"},
  {ARG_INTS, kInSliceUpper, "-1", "Upper-right corner of slice (-1 means no slicing)"},
  {ARG_STRING, kWithRecSet, kNotSpecified, "For binary operations, the second operand."},
  {ARG_STRING, kSeriesOut, kNotSpecified, "Name of series in which to save extracted data."},
  {ARG_INTS, kPosOut, "-1", "Bottom-left corner of the output data array at which data are to be written."},
  {ARG_STRING, kSegList, kNotSpecified, "Comma-separated list of segments on which to operate."},
  {ARG_STRING, kOp, "", "The operation to perform."},
  {ARG_DOUBLE, kBzero, "0.0", "For integer output, the bzero to use."},
  {ARG_DOUBLE, kBscale, "1.0", "For integer output, the bscale to use."},
  {ARG_NUME, kDatatype, "double", "Data type of in-memory data array.", "char,short,int,longlong,float,double,raw"},
  {ARG_FLAG, kNoSegsFlag, NULL, "Don't create an output segmetn file - just copy keywords.", NULL},
  {ARG_END}
};

char *module_name = "arithTool";

typedef struct DRMSContainer_struct
{
  HContainer_t *items;
  void (*Free)(HContainer_t *);
  HIterator_t iter;
} DRMSContainer_t;



static int NiceIntro()
{
     int usage = cmdparams_get_int(&cmdparams, "h", NULL);
     if (usage)
     {
	  printf ("Usage:\n\tarithTool [-h] "
		  "op=<operation> recSetIn=<recspec> [withRecSet=<recspec>] [seriesOut=<seriesname>] [segList=<segments>]\n"
		  "  details are:\n"
		  "  -h: help - show this message then exit\n"
		  "  <seriesname> - fully qualified series name.\n"
		  "  <timerange> - time value range set.\n"
		  "  seriesIn defaults to sdo.fds.\n"
		  "  timeRange defaults to all records in seriesIn.\n"
		  "  seriesOut defaults to sdo.fdsStateVectors.\n"
		  "  example - extractFdsStateV seriesIn=su_arta.TestFDSData timeRange=2006.11.20_22:38:00-2006.11.20_22:45:00,2006.11.20_22:52:00. seriesOut=su_arta.TestFDSHelio\n");
	  return 1;
     }

     return 0;
}

/* No need for an associative array since this mapping will be performed one time. */
static ArithOp_t MapOp(char *opStr)
{
     ArithOp_t op = kArithOpUnknown;

     if (strncmp(opStr, "mul", 3) == 0)
     {
	  op = kArithOpMul;
     }
     else if (strncmp(opStr, "div", 3) == 0)
     {
	  op = kArithOpDiv;
     }
     else if (strncmp(opStr, "add", 3) == 0)
     {
	  op = kArithOpAdd;
     }
     else if (strncmp(opStr, "sub", 3) == 0)
     {
	  op = kArithOpSub;
     }
     else if (strncmp(opStr, "abs", 3) == 0)
     {
	  op = kArithOpAbs;
     }
     else if (strncmp(opStr, "sqrt", 4) == 0)
     {
	  op = kArithOpSqrt;
     }
     else if (strncmp(opStr, "log", 3) == 0)
     {
	  op = kArithOpLog;
     }
     else if (strncmp(opStr, "pow", 3) == 0)
     {
	  op = kArithOpPow;
     }
     else if (strncmp(opStr, "sqr", 3) == 0)
     {
	  op = kArithOpSqr;
     }
     else if (strncmp(opStr, "recip", 5) == 0)
     {
	  op = kArithOpRecip;
     }
     else if (strncmp(opStr, "nomean", 6) == 0)
     {
	  op = kArithOpSubmean;
     }
     else if (strncmp(opStr, "nop", 3) == 0)
     {
	op = kArithOpNop;
     }

     return op;
}

/* hconFree performs non-standard cleaning up of HContainer_t. */
static int CreateDRMSPrimeKeyContainer(DRMS_Record_t *recTemplate, DRMSContainer_t *conPrimeKeys,
				       void (*hconFree)(HContainer_t *hc))
{
     int error = 0;
     int drmsst = DRMS_SUCCESS;

     /* There is no HContainer_t of primary keys - make one. */

     if (conPrimeKeys != NULL)
     {
	  /* Create new HContainer_t. */
	  conPrimeKeys->items = (HContainer_t *)malloc(sizeof(HContainer_t));

	  if (conPrimeKeys->items != NULL)
	  {
	       hcon_init(conPrimeKeys->items, sizeof(DRMS_Keyword_t *), DRMS_MAXKEYNAMELEN, NULL, NULL);
	       int iPrimeKeys = 0;
               int nkeys = 0;

               char **keyarr =
                 drms_series_createpkeyarray(drms_env, recTemplate->seriesinfo->seriesname, &nkeys, &drmsst);

	       while (iPrimeKeys < nkeys)
	       {
                  DRMS_Keyword_t *keyword = drms_keyword_lookup(recTemplate, keyarr[iPrimeKeys], 1);

		    if (keyword != NULL)
		    {
                       DRMS_Keyword_t **newKeyword =
                         (DRMS_Keyword_t **)hcon_allocslot_lower(conPrimeKeys->items,
                                                                 keyword->info->name);

                       if (newKeyword != NULL)
                       {
                          *newKeyword = keyword;
                       }
                       else
                       {
                          error = 1;
                          break;
                       }
		    }
		    else
		    {
                       error = 1;
                       break;
		    }

		    iPrimeKeys++;
	       }

               drms_series_destroypkeyarray(&keyarr, nkeys);

	       /* Create iterator too. */
	       conPrimeKeys->Free = hconFree;
	       hiter_new(&(conPrimeKeys->iter), conPrimeKeys->items);
	  }
     }
     else
     {
	  error = 1;
     }

     return error;
}

/* hconFree performs non-standard cleaning up of HContainer_t. */
static void CreateDRMSSegmentContainer(DRMS_Record_t *recTemplate, DRMSContainer_t *conSegs,
				       void (*hconFree)(HContainer_t *hc))
{
     /* Create iterator too. */
     conSegs->items = &(recTemplate->segments);
     conSegs->Free = hconFree;
     hiter_new(&(conSegs->iter), conSegs->items);
}

static void DestroyDRMSContainer(DRMSContainer_t *pCon)
{
   int wasinit = (pCon->items != NULL);

   if (pCon != NULL)
   {
      if (pCon->Free != NULL && pCon->items != NULL)
      {
         (*(pCon->Free))(pCon->items);
      }

      pCon->items = NULL;

      if (wasinit)
      {
         hiter_free(&(pCon->iter));
      }
   }
}

static void ReleaseHContainer(HContainer_t *hcon)
{
     if (hcon != NULL)
     {
	  /* Standard cleaning up of HContainer_t. */
	  hcon_free(hcon);
	  free(hcon);
     }
}

static int IsValidOp(ArithOp_t op)
{
     return (op >= kArithOpUnknown && op < kArithOpLastDummy);
}

static int IsBinaryOp(ArithOp_t op)
{
     if (IsValidOp(op))
     {
	  return op < kArithOpBinaryDummy;
     }

     return 0;
}

/* Returns 0 on error. */
/* Returns 1 if keys2 is equal to keys1. */
static int KeysEqual(DRMSContainer_t *keys1, DRMSContainer_t *keys2)
{
     int ret = 1;

     /*  Iterate through keys1, looking up the keyword name in the keys2 collection.
      *  This should be relatively efficient since the lookup is a hash function.
      */

     if (hcon_size(keys1->items) != hcon_size(keys2->items))
     {
	  ret = 0;
     }
     else
     {
	  hiter_rewind(&(keys1->iter));
	  DRMS_Keyword_t **currKey = NULL;
	  DRMS_Keyword_t **currKey2 = NULL;
	  while ((currKey = (DRMS_Keyword_t **)hiter_getnext(&(keys1->iter))) != NULL)
	  {
	       if ((currKey2 = hcon_lookup_lower(keys2->items, (*currKey)->info->name)) == NULL)
	       {
		    ret = 0;
		    break;
	       }
	       else
	       {
		    /* Ensure matching keys are the same type. */
		    if (drms_keyword_type(*currKey) != drms_keyword_type(*currKey2))
		    {
			 ret = 0;
			 break;
		    }
	       }
	  }
     }

     return ret;
}


/* Returns the number of matching segments.  If matchSet != NULL,
 * returns the actual segs1 segments
 * matching in matchSet.  Returns -1 on error.
 */
static int CreateMatchingSegs(DRMSContainer_t *segs1, DRMSContainer_t *segs2, DRMSContainer_t *matchSet)
{
     int nMatch = 0;

     hiter_rewind(&(segs1->iter));
     DRMS_Segment_t *currSeg = NULL;

     while ((currSeg = (DRMS_Segment_t *)hiter_getnext(&(segs1->iter))) != NULL)
     {
	  if (segs1 == segs2 ||
	      hcon_lookup_lower(segs2->items, currSeg->info->name) != NULL)
	  {
	       nMatch++;

	       if (matchSet != NULL)
	       {
		    if (nMatch == 1)
		    {
			 matchSet->items = (HContainer_t *)malloc(sizeof(HContainer_t));
			 if (matchSet->items != NULL)
			 {
			      hcon_init(matchSet->items, sizeof(DRMS_Segment_t *),
					DRMS_MAXSEGNAMELEN,
					NULL,
					NULL);
			 }
			 else
			 {
			      nMatch = -1;
			      break;
			 }
		    }

		    DRMS_Segment_t **newSeg =
		      (DRMS_Segment_t **)hcon_allocslot_lower(matchSet->items, currSeg->info->name);

		    if (newSeg != NULL)
		    {
			 *newSeg = currSeg;
		    }
		    else
		    {
			 nMatch = -1;
			 break;
		    }
	       }
	  }
     }

     if (nMatch > 0 && matchSet != NULL)
     {
	  matchSet->Free = ReleaseHContainer;
	  hiter_new(&(matchSet->iter), matchSet->items);
     }

     return nMatch;
}

static int ValidateBinaryOperands(char *recSetIn,
				  char **segNameArr,
				  int nSegs,
				  DRMSContainer_t *inKeys, DRMSContainer_t *inSegs,
				  DRMSContainer_t *withKeys, DRMSContainer_t *withSegs,
				  DRMSContainer_t *segsToProc)
{
     int error = 0;
     int status = 0;

     int nWithRecs = 0;
     DRMS_RecordSet_t *recSet = drms_open_records(drms_env, recSetIn, &status);
     error = (status != DRMS_SUCCESS);

     if (error == 0)
     {
	  nWithRecs = recSet->n;
	  drms_close_records(recSet, DRMS_FREE_RECORD);
     }

     /* If the number of records in the withSeries != 1, then the
      * primary keys must match.
      */
     if (!error)
     {
	  if (nWithRecs > 1)
	  {
	       if (!KeysEqual(inKeys, withKeys))
	       {
		    error = 1;
	       }
	  }
     }

     if (!error)
     {
	  int nMatch = 0;

	  if (nSegs > 0)
	  {

	       /* If a segment list is specified, then both the 'in' and 'with' series
		* must contain all items in the segment list
		*/
	       int iSeg = 0;

	       for (; iSeg < nSegs; iSeg++)
	       {
		    char *aSegName = segNameArr[iSeg];

		    if (hcon_lookup_lower(inSegs->items, aSegName) == NULL ||
			hcon_lookup_lower(withSegs->items, aSegName) == NULL)
		    {
			 error = 1;
			 fprintf(stderr, "Segment %s not present in an input recordset.\n",
				 aSegName);
			 break;
		    }
	       }

	       /* Create segsToProc list - a copy of segs in inSegs. */
	       if (!error)
	       {
		    nMatch = CreateMatchingSegs(inSegs, inSegs, segsToProc);
	       }
	  }
	  else
	  {
	       /* If no segments are specified, at least one segment must match */
	       nMatch = CreateMatchingSegs(inSegs, withSegs, segsToProc);
	  }

	  if (nMatch == -1)
	  {
	       DestroyDRMSContainer(segsToProc);
	       error = 1;
	  }

	  if (!error && nMatch == 0)
	  {
	       error = 1;
	  }
     }

     return error;
}

/* Must free malloc'd memory if ret != NULL. */
static char *CreateStrFromDRMSValue(DRMS_Record_t *rec, char *keyName, int *error)
{
     return drms_getkey_string(rec, keyName, error);
}

static int MultChar(const char *pInData, const char *pWithData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *rop = NULL;
   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]) &&
	     (!!drms_ismissing_char(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop * *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int MultShort(const short *pInData, const short *pWithData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *rop = NULL;
   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]) &&
	     (!!drms_ismissing_short(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop * *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int MultInt(const int *pInData, const int *pWithData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *rop = NULL;
   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]) &&
	     (!!drms_ismissing_int(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop * *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int MultLongLong(const long long *pInData,
			const long long *pWithData,
			arraylen_t nElements,
			long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *rop = NULL;
   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]) &&
	     (!!drms_ismissing_longlong(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop * *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int MultFloat(const float *pInData, const float *pWithData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *rop = NULL;
   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]) &&
	     (!!drms_ismissing_float(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop * *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int MultDouble(const double *pInData, const double *pWithData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *rop = NULL;
   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]) &&
	     (!!drms_ismissing_double(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop * *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int DivChar(const char *pInData, const char *pWithData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *rop = NULL;
   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]) &&
	     (!!drms_ismissing_char(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop / *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int DivShort(const short *pInData, const short *pWithData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *rop = NULL;
   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]) &&
	     (!!drms_ismissing_short(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop / *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int DivInt(const int *pInData, const int *pWithData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *rop = NULL;
   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]) &&
	     (!!drms_ismissing_int(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop / *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int DivLongLong(const long long *pInData,
		       const long long *pWithData,
		       arraylen_t nElements,
		       long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *rop = NULL;
   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]) &&
	     (!!drms_ismissing_longlong(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop / *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int DivFloat(const float *pInData, const float *pWithData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *rop = NULL;
   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]) &&
	     (!!drms_ismissing_float(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop / *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int DivDouble(const double *pInData, const double *pWithData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *rop = NULL;
   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]) &&
	     (!!drms_ismissing_double(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop / *rop;

	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AddChar(const char *pInData, const char *pWithData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *rop = NULL;
   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]) &&
	     (!!drms_ismissing_char(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop + *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AddShort(const short *pInData, const short *pWithData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *rop = NULL;
   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]) &&
	     (!!drms_ismissing_short(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop + *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AddInt(const int *pInData, const int *pWithData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *rop = NULL;
   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]) &&
	     (!!drms_ismissing_int(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop + *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AddLongLong(const long long *pInData,
		       const long long *pWithData,
		       arraylen_t nElements,
		       long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *rop = NULL;
   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]) &&
	     (!!drms_ismissing_longlong(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop + *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AddFloat(const float *pInData, const float *pWithData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *rop = NULL;
   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]) &&
	     (!!drms_ismissing_float(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop + *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AddDouble(const double *pInData, const double *pWithData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *rop = NULL;
   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]) &&
	     (!!drms_ismissing_double(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop + *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubChar(const char *pInData, const char *pWithData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *rop = NULL;
   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]) &&
	     (!!drms_ismissing_char(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop - *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubShort(const short *pInData, const short *pWithData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *rop = NULL;
   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]) &&
	     (!!drms_ismissing_short(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop - *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubInt(const int *pInData, const int *pWithData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *rop = NULL;
   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]) &&
	     (!!drms_ismissing_int(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop - *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubLongLong(const long long *pInData,
		       const long long *pWithData,
		       arraylen_t nElements,
		       long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *rop = NULL;
   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]) &&
	     (!!drms_ismissing_longlong(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop - *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubFloat(const float *pInData, const float *pWithData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *rop = NULL;
   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]) &&
	     (!!drms_ismissing_float(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop - *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubDouble(const double *pInData, const double *pWithData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *rop = NULL;
   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (!pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]) &&
	     (!!drms_ismissing_double(pWithData[index])))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    rop = &(pWithData[index]);
	    *result = *lop - *rop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AbsChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = abs(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AbsShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = abs(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AbsInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = abs(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AbsLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = llabs(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AbsFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = fabsf(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int AbsDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = fabs(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrtChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = sqrt(abs(*lop));
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrtShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = sqrt(abs(*lop));
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrtInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = sqrt(abs(*lop));
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrtLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = sqrt(llabs(*lop));
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrtFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = sqrt(fabsf(*lop));
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrtDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = sqrt(fabs(*lop));
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int LogChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = log10(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int LogShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = log10(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int LogInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = log10(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int LogLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = log10(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int LogFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = log10f(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int LogDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = log10(*lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int PowChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = pow(10.0, *lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int PowShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = pow(10.0, *lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int PowInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = pow(10.0, *lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int PowLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = pow(10.0, *lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int PowFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = powf(10.0, *lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int PowDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = pow(10.0, *lop);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop * *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop * *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop * *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop * *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop * *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SqrDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop * *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int RecipChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = ((*lop != 0.0) ? 1.0 / *lop : DRMS_MISSING_CHAR);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int RecipShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = ((*lop != 0.0) ? 1.0 / *lop : DRMS_MISSING_SHORT);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int RecipInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = ((*lop != 0.0) ? 1.0 / *lop : DRMS_MISSING_INT);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int RecipLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = ((*lop != 0.0) ? 1.0 / *lop : DRMS_MISSING_LONGLONG);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int RecipFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = ((*lop != 0.0) ? 1.0 / *lop : DRMS_MISSING_FLOAT);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int RecipDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = ((*lop != 0.0) ? 1.0 / *lop : DRMS_MISSING_DOUBLE);
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubmeanChar(const char *pInData, arraylen_t nElements, char *pOutData, double mean)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop - mean;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubmeanShort(const short *pInData, arraylen_t nElements, short *pOutData, double mean)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop - mean;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubmeanInt(const int *pInData, arraylen_t nElements, int *pOutData, double mean)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop - mean;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubmeanLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData, double mean)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop - mean;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubmeanFloat(const float *pInData, arraylen_t nElements, float *pOutData, double mean)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop - mean;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int SubmeanDouble(const double *pInData, arraylen_t nElements, double *pOutData, double mean)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop - mean;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int NopChar(const char *pInData, arraylen_t nElements, char *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const char *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      char *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_CHAR;

	 if (!drms_ismissing_char(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int NopShort(const short *pInData, arraylen_t nElements, short *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const short *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      short *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_SHORT;

	 if (!drms_ismissing_short(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int NopInt(const int *pInData, arraylen_t nElements, int *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const int *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      int *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_INT;

	 if (!drms_ismissing_int(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int NopLongLong(const long long *pInData, arraylen_t nElements, long long *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const long long *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      long long *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_LONGLONG;

	 if (!drms_ismissing_longlong(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int NopFloat(const float *pInData, arraylen_t nElements, float *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const float *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      float *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_FLOAT;

	 if (!drms_ismissing_float(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

static int NopDouble(const double *pInData, arraylen_t nElements, double *pOutData)
{
   int error = 0;

   int first = 1;
   int percentDone = 0;
   int update = 0;

   const double *lop = NULL;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      double *result = NULL;
      arraylen_t index = 0;

      for (; index < nElements; index++)
      {
	 result = &(pOutData[index]);
	 *result = DRMS_MISSING_DOUBLE;

	 if (!drms_ismissing_double(pInData[index]))
	 {
	    percentDone = index * 100/nElements;

	    if (first)
	    {
	       fprintf(stdout, "%03d%% done", percentDone);
	       first = 0;
	    }
	    else if (percentDone == update)
	    {
	       update++;
	       fprintf(stdout,
		       "\b\b\b\b\b\b\b\b\b%03d%% done",
		       percentDone);
	       fflush(stdout);
	    }

	    lop = &(pInData[index]);
	    *result = *lop;
	 }
      }

      fprintf(stdout, "\b\b\b\b\b\b\b\b\b100%% done\n");
      fflush(stdout);
   }

   return error;
}

/* Gets as far as it can, filling in with MISSINGS if necessary.  Returns 1
* if any error happens. */
static int PerformOperation(ArithOp_t op, DRMS_Type_t dtype,
			    const void *pInData, const void *pWithData,
			    arraylen_t nElements, void *pOutData, double mean)
{
   int error = 0;

   if (!pInData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else if (IsBinaryOp(op) && !pWithData)
   {
      error = 1;
      fprintf(stderr, "Missing input data.\n");
   }
   else
   {
      if (IsBinaryOp(op))
      {
	 switch (op)
	 {
	    case kArithOpMul:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return MultChar(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return MultShort(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return MultInt(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return MultLongLong(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return MultFloat(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return MultDouble(pInData, pWithData, nElements, pOutData);
	      }
	      break;
	    case kArithOpDiv:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return DivChar(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return DivShort(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return DivInt(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return DivLongLong(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return DivFloat(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return DivDouble(pInData, pWithData, nElements, pOutData);
	      }
	      break;
	    case kArithOpAdd:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return AddChar(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return AddShort(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return AddInt(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return AddLongLong(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return AddFloat(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return AddDouble(pInData, pWithData, nElements, pOutData);
	      }
	      break;
	    case kArithOpSub:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return SubChar(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return SubShort(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return SubInt(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return SubLongLong(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return SubFloat(pInData, pWithData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return SubDouble(pInData, pWithData, nElements, pOutData);
	      }
	      break;
	    default:
	      error = 1;
	      fprintf(stderr, "Expecting a binary operator, got %d instead.\n", (int)op);
	 }
      }
      else
      {
	 switch (op)
	 {
	    case kArithOpAbs:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return AbsChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return AbsShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return AbsInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return AbsLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return AbsFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return AbsDouble(pInData, nElements, pOutData);
	      }
	      break;
	    case kArithOpSqrt:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return SqrtChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return SqrtShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return SqrtInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return SqrtLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return SqrtFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return SqrtDouble(pInData, nElements, pOutData);
	      }
	      break;
	    case kArithOpLog:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return LogChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return LogShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return LogInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return LogLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return LogFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return LogDouble(pInData, nElements, pOutData);
	      }
	      break;
	    case kArithOpPow:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return PowChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return PowShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return PowInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return PowLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return PowFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return PowDouble(pInData, nElements, pOutData);
	      }
	      break;
	    case kArithOpSqr:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return SqrChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return SqrShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return SqrInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return SqrLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return SqrFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return SqrDouble(pInData, nElements, pOutData);
	      }
	      break;
	    case kArithOpRecip:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return RecipChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return RecipShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return RecipInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return RecipLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return RecipFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return RecipDouble(pInData, nElements, pOutData);
	      }
	      break;
	    case kArithOpSubmean:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return SubmeanChar(pInData, nElements, pOutData, mean);
		 case DRMS_TYPE_SHORT:
		   return SubmeanShort(pInData, nElements, pOutData, mean);
		 case DRMS_TYPE_INT:
		   return SubmeanInt(pInData, nElements, pOutData, mean);
		 case DRMS_TYPE_LONGLONG:
		   return SubmeanLongLong(pInData, nElements, pOutData, mean);
		 case DRMS_TYPE_FLOAT:
		   return SubmeanFloat(pInData, nElements, pOutData, mean);
		 case DRMS_TYPE_DOUBLE:
		   return SubmeanDouble(pInData, nElements, pOutData, mean);
	      }
	      break;
	    case kArithOpNop:
	      switch (dtype)
	      {
		 case DRMS_TYPE_CHAR:
		   return NopChar(pInData, nElements, pOutData);
		 case DRMS_TYPE_SHORT:
		   return NopShort(pInData, nElements, pOutData);
		 case DRMS_TYPE_INT:
		   return NopInt(pInData, nElements, pOutData);
		 case DRMS_TYPE_LONGLONG:
		   return NopLongLong(pInData, nElements, pOutData);
		 case DRMS_TYPE_FLOAT:
		   return NopFloat(pInData, nElements, pOutData);
		 case DRMS_TYPE_DOUBLE:
		   return NopDouble(pInData, nElements, pOutData);
	      }
	      break;
	    default:
	      error = 1;
	      fprintf(stderr, "Expecting a unary operator, got %d instead.\n", (int)op);
	 }
      }
   }

   return error;
}

/* XXX - WARNING: Not modified to work with non-double data!!! But DoUnaryOp was */
int DoBinaryOp(DRMS_Env_t *drmsEnv, ArithOp_t op, DRMS_Type_t dtype,
	       DRMSContainer_t *inPrimeKeys, DRMSContainer_t *segsToProc,
	       char *inSeriesQuery,
	       char *withSeriesName, char *withSeriesQuery,
	       char *seriesOut, double bzero, double bscale)
{
     int error = 0;
     int status = 0;
     int nWithRecs = 0;
     int nInRecs = 0;
     double **pOutData = NULL;
     char **outSegNames = NULL;

     /* Open the inSeries. */
     DRMS_RecordSet_t *inRecSet = drms_open_records(drmsEnv, inSeriesQuery, &status);
     error = (status != DRMS_SUCCESS);
     DRMS_RecordSet_t *withRecSet = NULL;

     DRMS_Record_t *inRec = NULL;
     DRMS_Record_t *withRec = NULL;

     int nSegs = segsToProc->items->num_total;

     if (!error)
     {
	  nInRecs = inRecSet->n;
	  withRecSet = drms_open_records(drmsEnv, withSeriesQuery, &status);
	  error = (status != DRMS_SUCCESS);
	  nWithRecs = withRecSet->n;

	  XASSERT(nWithRecs >= 1);

	  if (nWithRecs == 0)
	  {
	       error = 1;
	  }
	  else if (nWithRecs == 1)
	  {
	       withRec = withRecSet->records[0];
	  }
     }

     int iRec = 0;
     for (; !error && iRec < nInRecs; iRec++)
     {
	  inRec = inRecSet->records[iRec];

	  fprintf(stdout, "Processing record %lld\n", inRec->recnum);

	  if (nWithRecs != 1)
	  {
	       /* We want to choose one record from all withSeries recs.  This
		* chosen record must match the prime keys of the current inSeries
		* record.  Since the inSeries and withSeries have the same prime
		* keys, because the prime key values select one record from the inSeries,
		* they also select one record from the withSeries.
		*
		* First close the withSeries record that matched
		* the previous inSeries record. */
	       drms_close_records(withRecSet, DRMS_FREE_RECORD);
	       withRecSet = NULL;

	       /* Query the withSeries to find the record whose prime keys match
		* the current inSeries record's prime key values.
		*/
	       char query[kMaxQuery];
	       char buf[kMaxQuery];
	       char *pQuery = query;
	       int maxQuery = kMaxQuery;

	       snprintf(buf, sizeof(buf), "%s", withSeriesName);
	       snprintf(pQuery, maxQuery, "%s", buf);
	       pQuery += strlen(buf);
	       maxQuery -= strlen(buf);

	       hiter_rewind(&(inPrimeKeys->iter));
	       DRMS_Keyword_t **primeKey = NULL;

	       while (maxQuery > 0 &&
		      (primeKey = hiter_getnext(&(inPrimeKeys->iter))) != NULL)
	       {
		    char *val = CreateStrFromDRMSValue(inRec, (*primeKey)->info->name, &error);
		    if (!error)
		    {
			 if (val != NULL)
			 {
			      snprintf(buf, sizeof(buf), "[%s=%s]",
					    (*primeKey)->info->name,
				       val);
			      snprintf(pQuery, maxQuery, "%s", buf);
			      pQuery += strlen(buf);
			      maxQuery -= strlen(buf);

			      if (maxQuery <= 0)
			      {
				   error = 1;
				   break;
			      }

			      free(val);
			 }
			 else
			      {
				   error = 1;
			      }
		    }
	       }

	       if (!error)
	       {
		    withRecSet = drms_open_records(drmsEnv, query, &status);
		    error = (status != DRMS_SUCCESS);

		    if (!error && withRecSet != NULL)
		    {
			 if (withRecSet -> n != 1)
			 {
			      error = 1;
			 }
			 else
			 {
			      withRec = withRecSet->records[0];
			 }
		    }
		    else
		    {
			 error = 1;
		    }
	       }
	  }

	  if (!error)
	  {
	       /* Have a single record from inSeries and a matching single
		* record from withSeries.  Get matching segments from each record. */
	       hiter_rewind(&(segsToProc->iter));
	       DRMS_Segment_t **seg = NULL;
	       DRMS_Segment_t *inSeg = NULL;
	       DRMS_Segment_t *withSeg = NULL;
	       DRMS_Array_t *inSegArray = NULL;
	       DRMS_Array_t *withSegArray = NULL;

	       pOutData = (double **)malloc(sizeof(double *) * nSegs);
	       outSegNames = (char **)malloc(sizeof(char *) * nSegs);

	       unsigned int iSeg = 0;
	       while(pOutData != NULL &&
		     outSegNames != NULL &&
		     !error &&
		     (seg = (DRMS_Segment_t **)hiter_getnext(&(segsToProc->iter))) != NULL)
	       {
		    fprintf(stdout, "  processing segment %s:  ", (*seg)->info->name);

		    inSeg = drms_segment_lookup(inRec, (*seg)->info->name);
		    withSeg = drms_segment_lookup(withRec, (*seg)->info->name);


		    XASSERT(inSeg != NULL && withSeg != NULL);
		    if (inSeg != NULL && withSeg != NULL)
		    {
			 inSegArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
			 error = (status != DRMS_SUCCESS);
			 if (!error)
			 {
			      withSegArray = drms_segment_read(withSeg, DRMS_TYPE_DOUBLE, &status);
			      error = (status != DRMS_SUCCESS);
			 }
		    }

		    if (!error)
		    {
			 /* Have inSeries seg data and withSeries seg data - do binary op. */
			 arraylen_t nElementsIn = drms_array_count(inSegArray);
			 arraylen_t nElementsWith = drms_array_count(withSegArray);
			 double *pInData = inSegArray->data;
			 double *pWithData = withSegArray->data;

			 if (nElementsIn == nElementsWith)
			 {
			      pOutData[iSeg] = (double *)malloc(sizeof(double) * nElementsIn);
			      outSegNames[iSeg] = strdup((*seg)->info->name);

			      if (pOutData[iSeg] != NULL && outSegNames[iSeg] != NULL)
			      {
				 /* Shouldn't return an error - all params checked */
				 PerformOperation(op,
						  dtype,
						  pInData,
						  pWithData,
						  nElementsIn,
						  pOutData[iSeg],
						  0.0);
				 iSeg++;
			      }
			      else
			      {
				   error = 1;
			      }
			 }
			 else
			 {
			      error = 1;
			 }
		    }

		    if (inSegArray)
		    {
		       drms_free_array(inSegArray);
		    }

		    if (withSegArray)
		    {
		       drms_free_array(withSegArray);
		    }

	       } /* while */
	  }

	  if (!error)
	  {
	       // CreateRecord(drmsEnv, inRec, outSegNames, pOutData, seriesOut, seriesOutExists);
	       /* If a record with the same primary key already exists, the new record will be
		* a newer version of the original one.
		*/
	       DRMS_Record_t *rec = drms_create_record(drmsEnv, seriesOut, DRMS_PERMANENT, &status);
	       error = (status != DRMS_SUCCESS);
	       if (rec != NULL && !error)
	       {
		    /* Write out outSeries keys - fill them in with values from inSeries. */
		    HIterator_t hit;
		    hiter_new(&hit, &(rec->keywords));
		    DRMS_Keyword_t *outKey = NULL;
		    DRMS_Keyword_t *inKey = NULL;

		    while ((outKey = (DRMS_Keyword_t *)hiter_getnext(&hit)) != NULL)
		    {
			 inKey = drms_keyword_lookup(inRec, outKey->info->name, 1);
			 if (inKey != NULL && !drms_keyword_isindex(inKey))
			 {
			      status = drms_setkey(rec,
						   outKey->info->name,
						   outKey->info->type,
						   &(inKey->value));

			      error = (status != DRMS_SUCCESS);
			 }
		    }

		    /* Write segments. */
		    unsigned int index = 0;
		    for (; index < nSegs; index++)
		    {
			 if (pOutData != NULL && outSegNames != NULL)
			 {

			      DRMS_Segment_t *inSeg = drms_segment_lookup(inRec,
									  outSegNames[index]);
			      DRMS_Segment_t *seg = drms_segment_lookup(rec, outSegNames[index]);
			      if (inSeg != NULL && seg != NULL)
			      {
				   /* segArray now owns pOutData[index]. */
				   DRMS_Array_t *segArray = drms_array_create(DRMS_TYPE_DOUBLE,
									      inSeg->info->naxis,
									      inSeg->axis,
									      pOutData[index],
									      &status);
				   error = (status != DRMS_SUCCESS);

				   if (!error)
				   {
				      segArray->bzero = bzero;
				      segArray->bscale = bscale;
				      segArray->israw = 0;
				      drms_segment_writewithkeys(seg, segArray, 0);
				      drms_free_array(segArray);
				   }
			      }

			 }
		    }

		    if (!error)
		    {
			 drms_close_record(rec, DRMS_INSERT_RECORD);
		    }
		    else
		    {
			 drms_close_record(rec, DRMS_FREE_RECORD);
		    }
	       }
	       else if (rec != NULL)
	       {
		    drms_close_record(rec, DRMS_FREE_RECORD);
	       }
	  }

	  /* Clean up seg names. DO NOT clean up pOutData - this got transferred to segArray
	   * just prior to writing the segments.
	   */
	  unsigned int index = 0;
	  for (; index < nSegs; index++)
	  {
	       if (outSegNames[index] != NULL)
	       {
		    char *pName = outSegNames[index];
		    free(pName);
	       }
	  }

	  if (outSegNames != NULL)
	  {
	       free(outSegNames);
	  }
     } /* foreach (record) */

     if (inRecSet != NULL)
     {
	  drms_close_records(inRecSet, DRMS_FREE_RECORD);
     }

     if (withRecSet != NULL)
     {
	  drms_close_records(withRecSet, DRMS_FREE_RECORD);
     }

     return error;
}

static int DoUnaryOp(DRMS_Env_t *drmsEnv, ArithOp_t op, DRMS_Type_t dtype,
		     DRMSContainer_t *segsToProc,
		     char *inSeriesQuery, int *slicelower, int *sliceupper,
		     char *seriesOut, int *pout, double bzero, double bscale, int nosegs)
{
     int error = 0;
     int status = 0;
     void *pOutData = NULL;
     char **outSegNames = NULL;
     int nSegs = 0;
     DRMS_Type_t actualtype;
    double *insegBZERO = NULL;
    double *insegBSCALE = NULL;

     /* Open the inSeries. */
     DRMS_RecordSet_t *inRecSet = drms_open_recordset(drmsEnv, inSeriesQuery, &status);
     DRMS_Record_t *inRec = NULL;

     DRMS_RecChunking_t cstat = kRecChunking_None;
     DRMS_Record_t *targetrec = NULL;

     error = (status != DRMS_SUCCESS);

     if (!error)
     {
         if (segsToProc && segsToProc->items)
         {
             nSegs = segsToProc->items->num_total;
         }
         else
         {
             nSegs = 0;
         }
     }

     if (nSegs > 0)
     {
         insegBZERO = (double *)malloc(sizeof(double) * nSegs);
         insegBSCALE = (double *)malloc(sizeof(double) * nSegs);
     }

     while ((inRec = drms_recordset_fetchnext(drmsEnv, inRecSet, &status, &cstat, NULL)) != NULL)
     {
	  fprintf(stdout, "Processing record %lld\n", inRec->recnum);

          /* If not processing segments, then simply copy DRMS keywords. */
          if (nosegs)
          {
             targetrec = drms_create_record(drmsEnv,
                                            seriesOut,
                                            DRMS_PERMANENT,
                                            &status);
             error = (status != DRMS_SUCCESS);

             if (!error)
             {
                drms_copykeys(targetrec, inRec, 0, kDRMS_KeyClass_Explicit);

                /* write DATE keyword */
                drms_keyword_setdate(targetrec);
                drms_close_record(targetrec, DRMS_INSERT_RECORD);
             }
             else
             {
                drms_close_record(targetrec, DRMS_FREE_RECORD);
             }

             continue;
          }

	  hiter_rewind(&(segsToProc->iter));
	  DRMS_Segment_t **seg = NULL;
	  DRMS_Segment_t *inSeg = NULL;
	  DRMS_Array_t *inSegArray = NULL;

	  actualtype = dtype;

	  pOutData = malloc(sizeof(void *) * nSegs);
	  outSegNames = (char **)malloc(sizeof(char *) * nSegs);
          memset(outSegNames, 0, sizeof(char *) * nSegs);

	  unsigned int iSeg = 0;
	  while(pOutData != NULL &&
		outSegNames != NULL &&
		!error &&
		(seg = (DRMS_Segment_t **)hiter_getnext(&(segsToProc->iter))) != NULL)
	  {
	       fprintf(stdout, "  processing segment %s:  ", (*seg)->info->name);
	       inSeg = drms_segment_lookup(inRec, (*seg)->info->name);

	       XASSERT(inSeg != NULL);
	       if (inSeg != NULL)
	       {
                  if (slicelower && sliceupper)
                  {
                     inSegArray = drms_segment_readslice(inSeg,
                                                         dtype,
                                                         slicelower,
                                                         sliceupper,
                                                         &status);
                  }
                  else
                  {
		    inSegArray = drms_segment_read(inSeg, dtype, &status);
                  }

		    if (dtype == DRMS_TYPE_RAW)
		    {
		       actualtype = inSeg->info->type;
		    }

                    if (status != DRMS_SUCCESS)
                    {
                       /* Error reading segment file - just skip and go to the next. */
                       fprintf(stderr, "Error reading segment - skipping to the next segment.\n");
                       continue;
                       iSeg++;
                    }

                    /* Have inSeries seg data - do unary op. */
                    arraylen_t nElementsIn = drms_array_count(inSegArray);
                    void *pInData = inSegArray->data;
                    void *pOut = NULL;

                    insegBZERO[iSeg] = inSegArray->bzero;
                    insegBSCALE[iSeg] = inSegArray->bscale;

                    switch (actualtype)
                    {
                       case DRMS_TYPE_CHAR:
                         ((char **)(pOutData))[iSeg] = (char *)malloc(sizeof(char) * nElementsIn);
                         pOut = (void *)(((char **)(pOutData))[iSeg]);
                         break;
                       case DRMS_TYPE_SHORT:
                         ((short **)(pOutData))[iSeg] = (short *)malloc(sizeof(short) * nElementsIn);
                         pOut = (void *)(((short **)(pOutData))[iSeg]);
                         break;
                       case DRMS_TYPE_INT:
                         ((int **)(pOutData))[iSeg] = (int *)malloc(sizeof(int) * nElementsIn);
                         pOut = (void *)(((int **)(pOutData))[iSeg]);
                         break;
                       case DRMS_TYPE_LONGLONG:
                         ((long long **)(pOutData))[iSeg] = (long long *)malloc(sizeof(long long) * nElementsIn);
                         pOut = (void *)(((long long **)(pOutData))[iSeg]);
                         break;
                       case DRMS_TYPE_FLOAT:
                         ((float **)(pOutData))[iSeg] = (float *)malloc(sizeof(float) * nElementsIn);
                         pOut = (void *)(((float **)(pOutData))[iSeg]);
                         break;
                       case DRMS_TYPE_DOUBLE:
                         ((double **)(pOutData))[iSeg] = (double *)malloc(sizeof(double) * nElementsIn);
                         pOut = (void *)(((double **)(pOutData))[iSeg]);
                         break;
                    }

                    outSegNames[iSeg] = strdup((*seg)->info->name);

                    if (pOut != NULL && outSegNames[iSeg] != NULL)
                    {
                       double mean = 0.0;

                       /* calc mean, if op == submean */
                       if (op == kArithOpSubmean)
                       {
                          int n = 0;
                          arraylen_t iData = 0;

                          switch (actualtype)
                          {
                             case DRMS_TYPE_CHAR:
                               {
                                  char *pData = pInData;
                                  for (; iData < nElementsIn; iData++, pData++)
                                  {
                                     if (!drms_ismissing_char(*pData))
                                     {
                                        mean += (double)*pData;
                                        n += 1;
                                     }
                                  }
                               }
                               break;
                             case DRMS_TYPE_SHORT:
                               {
                                  short *pData = pInData;
                                  for (; iData < nElementsIn; iData++, pData++)
                                  {
                                     if (!drms_ismissing_short(*pData))
                                     {
                                        mean += (double)*pData;
                                        n += 1;
                                     }
                                  }
                               }
                               break;
                             case DRMS_TYPE_INT:
                               {
                                  int *pData = pInData;
                                  for (; iData < nElementsIn; iData++, pData++)
                                  {
                                     if (!drms_ismissing_int(*pData))
                                     {
                                        mean += (double)*pData;
                                        n += 1;
                                     }
                                  }
                               }
                               break;
                             case DRMS_TYPE_LONGLONG:
                               {
                                  long long *pData = pInData;
                                  for (; iData < nElementsIn; iData++, pData++)
                                  {
                                     if (!drms_ismissing_longlong(*pData))
                                     {
                                        mean += (double)*pData;
                                        n += 1;
                                     }
                                  }
                               }
                               break;
                             case DRMS_TYPE_FLOAT:
                               {
                                  float *pData = pInData;
                                  for (; iData < nElementsIn; iData++, pData++)
                                  {
                                     if (!drms_ismissing_float(*pData))
                                     {
                                        mean += (double)*pData;
                                        n += 1;
                                     }
                                  }
                               }
                               break;
                             case DRMS_TYPE_DOUBLE:
                               {
                                  double *pData = pInData;
                                  for (; iData < nElementsIn; iData++, pData++)
                                  {
                                     if (!drms_ismissing_double(*pData))
                                     {
                                        mean += (double)*pData;
                                        n += 1;
                                     }
                                  }
                               }
                               break;
                          }

                          if (n > 0)
                          {
                             mean /= n;
                          }
                       }

                       PerformOperation(op,
                                        actualtype,
                                        pInData,
                                        NULL,
                                        nElementsIn,
                                        pOut,
                                        mean);
                    }
                    else
                    {
                       error = 1;
                    }


		    drms_free_array(inSegArray);
	       }

	       iSeg++;
	  } /* while */

	  if (!error)
	  {
               DRMS_Record_t *rec = NULL;

               /* See if we have the record-creation permission */
               if (drms_series_cancreaterecord(drmsEnv, seriesOut))
               {

                  /* If a record with the same primary key already exists, the new record will be
                   * a newer version of the original one.
                   */
                  rec = drms_create_record(drmsEnv,
                                           seriesOut,
                                           DRMS_PERMANENT,
                                           &status);
                  error = (status != DRMS_SUCCESS);
               }
               else
               {
                  fprintf(stderr, "Can't create new records in '%s'; permission denied.\n", seriesOut);
                  error = 1;
               }

	       if (rec != NULL && !error)
	       {
		    /* Write out outSeries keys - fill them in with values from inSeries. */
		    HIterator_t hit;
		    hiter_new(&hit, &(rec->keywords));
		    DRMS_Keyword_t *outKey = NULL;
		    DRMS_Keyword_t *inKey = NULL;

		    while ((outKey = (DRMS_Keyword_t *)hiter_getnext(&hit)) != NULL)
		    {
			 inKey = drms_keyword_lookup(inRec, outKey->info->name, 1);
			 if (inKey != NULL && !drms_keyword_isindex(inKey))
			 {
			      status = drms_setkey(rec,
						   outKey->info->name,
						   outKey->info->type,
						   &(inKey->value));

			      error = (status != DRMS_SUCCESS);
			 }
		    }

                    hiter_free(&hit);

		    /* Write segments, and then clean up arrays holding data and seg names. */

        /* write DATE keyword (BEFORE writing segment); do this before segment loop so all segment FITS files will have
         * the same DATE keyword value that matches the value is in the DRMS DB
         */
        drms_keyword_setdate(rec);

		    unsigned int index = 0;
		    for (; index < nSegs; index++)
		    {
			 if (pOutData != NULL && outSegNames != NULL)
			 {
			    DRMS_Segment_t *inseg = NULL;
			    DRMS_Segment_t *outseg = NULL;


                            if (outSegNames[index])
                            {
                               /* If outSegNames[index] == NULL this means that there was
                                * a problem reading the source segment file, so we're
                                * skipping this segment. */
                               inseg = drms_segment_lookup(inRec, outSegNames[index]);
                               outseg = drms_segment_lookup(rec, outSegNames[index]);

                               if (nSegs == 1 && (inseg == NULL || outseg == NULL))
                               {
                                  /* Perhaps the two series have just one segment, but the names
                                   * don't match - go ahead and use the one outseries segment
                                   * to hold the results. */
                                  inseg = drms_segment_lookupnum(inRec, 0);
                                  outseg = drms_segment_lookupnum(rec, 0);
                               }
                            }

			    if (inseg != NULL && outseg != NULL)
			    {
			       /* segArray now owns pOutData[index]. */
			       void *pOut = NULL;

			       switch (actualtype)
			       {
				  case DRMS_TYPE_CHAR:
				    pOut = (void *)(((char **)(pOutData))[index]);
				    break;
				  case DRMS_TYPE_SHORT:
				    pOut = (void *)(((short **)(pOutData))[index]);
				    break;
				  case DRMS_TYPE_INT:
				    pOut = (void *)(((int **)(pOutData))[index]);
				    break;
				  case DRMS_TYPE_LONGLONG:
				    pOut = (void *)(((long long **)(pOutData))[index]);
				    break;
				  case DRMS_TYPE_FLOAT:
				    pOut = (void *)(((float **)(pOutData))[index]);
				    break;
				  case DRMS_TYPE_DOUBLE:
				    pOut = (void *)(((double **)(pOutData))[index]);
				    break;
			       }

			       DRMS_Array_t *segArray = NULL;
                               int actaxis[DRMS_MAXRANK];
                               int end[DRMS_MAXRANK];
                               int iaxis;

                               if (slicelower && sliceupper)
                               {
                                  /* Must calculate new axis lengths */
                                  for (iaxis = 0; iaxis < inseg->info->naxis; iaxis++)
                                  {
                                     actaxis[iaxis] = sliceupper[iaxis] - slicelower[iaxis] + 1;
                                  }
                               }
                               else
                               {
                                  memcpy(actaxis, inseg->axis, sizeof(int) * inseg->info->naxis);
                               }

                               if (pout)
                               {
                                  for (iaxis = 0; iaxis < inseg->info->naxis; iaxis++)
                                  {
                                     end[iaxis] = actaxis[iaxis] + pout[iaxis] - 1;
                                  }
                               }

                               segArray = drms_array_create(actualtype,
                                                            inseg->info->naxis,
                                                            actaxis,
                                                            pOut,
                                                            &status);

			       error = (status != DRMS_SUCCESS);

			       if (!error)
			       {
				  if (dtype == DRMS_TYPE_RAW)
				  {
                                     /* We read the input array as RAW, which means the input data are scaled values
                                      * with bzero and bscale as the scaling parameters. The output array will have
                                      * been some manipulation of the input data, so the input's scaling parameters
                                      * are still relevant. */
				     segArray->bzero = insegBZERO[index];
				     segArray->bscale = insegBSCALE[index];
				     segArray->israw = 1;
				  }
				  else
				  {
                                     /* The input array data are in physical values (not RAW). If the output segment's
                                      * data type is an int, then the inverse bzero/bscale conversion will be applied,
                                      * using the values of segArray->bzero and segArray->bscale. If the output segment's
                                      * data type is a float, then segArray->bzero and segArray->bscale are ignored. */
				     segArray->bzero = bzero;
				     segArray->bscale = bscale;
				     segArray->israw = 0;
				  }

          if (pout == NULL)
          {
             drms_segment_writewithkeys(outseg, segArray, 0);
          }
          else
          {
             drms_segment_writeslice(outseg, segArray, pout, end, 0);
          }

				  drms_free_array(segArray);
			       }
			    }
			 }
		    }

		    if (!error)
		    {
		      drms_close_record(rec, DRMS_INSERT_RECORD);
		    }
		    else
		    {
			 drms_close_record(rec, DRMS_FREE_RECORD);
		    }
	       }
	       else if (rec != NULL)
	       {
		    drms_close_record(rec, DRMS_FREE_RECORD);
	       }
	  }

	  /* Clean up seg names. DO NOT clean up pOutData - this got transferred to segArray
	   * just prior to writing the segments.
	   */
	  unsigned int index = 0;
	  for (; index < nSegs; index++)
	  {
	       if (outSegNames[index] != NULL)
	       {
		    char *pName = outSegNames[index];
		    free(pName);
	       }
	  }

	  if (outSegNames != NULL)
	  {
	       free(outSegNames);
	  }

          if (pOutData)
          {
             free(pOutData);
             pOutData = NULL;
          }
     } /* foreach (record) */

     if (insegBZERO)
     {
	free(insegBZERO);
     }

     if (insegBSCALE)
     {
	free(insegBSCALE);
     }

     if (inRecSet != NULL)
     {
	  drms_close_records(inRecSet, DRMS_FREE_RECORD);
     }

     if (pOutData)
     {
        free(pOutData);
     }

     return error;
}

#if TESTER
     typedef struct yammer
     {
       int a;
       char *b;
     } yammer_t;

     void FreeF(const void *v)
     {
	  free(((yammer_t *)(v))->b);
     };

     void CopyF(const void *dst, const void *src)
     {
	  ((yammer_t *)dst)->b = strdup(((yammer_t *)src)->b);
     };
#endif

int DoIt(void)
{

#if TESTER
     char *nameArr[4] = {"name1", "name2", "name3", "name4"};
     char *strArr[4] = {"yammer1", "yammer2", "yammer3", "yammer4"};
     yammer_t valArr[4];

     valArr[0].a = 1;
     valArr[0].b = strArr[0];
     valArr[1].a = 2;
     valArr[1].b = strArr[1];
     valArr[2].a = 3;
     valArr[2].b = strArr[2];
     valArr[3].a = 4;
     valArr[3].b = strArr[3];

     yammer_t *valArrI[4] = {&valArr[0], &valArr[1], &valArr[2], &valArr[3]};

     HContainer_t *cont = hcon_create(sizeof(yammer_t),
				      128,
				      FreeF,
				      CopyF,
				      (void **)valArrI,
				      nameArr,
				      4);
     if (cont)
     {
	  HIterator_t *it = hiter_create(cont);
	  if (it)
	  {
	       yammer_t *yp = NULL;
	       while ((yp = (yammer_t *)hiter_getnext(it)) != NULL)
	       {
		    printf("int val %d, str val %s\n", yp->a, yp->b);
	       }

	       hiter_destroy(&it);
	  }

	  hcon_destroy(&cont);
     }

     printf("Aok\n");

     return 1;
#endif
     int error = 0;
     int status = 0;

     if (NiceIntro())
     {
	  return 0;
     }

     if (drms_env == NULL)
     {
	  error = 1;
     }
     else
     {
	  /* Check for a valid combination of parameters. */
	  int bBinaryOp = 0;
	  int bSeriesOut = 0;       /* user specified an out series, not just the default */
	  int bWithRecSet = 0;
	  int bSeriesOutExists = 0; /* the out series specified by the user exists
				     * and matches the inSeries in terms of
				     * keywords. */
	  int bSeriesOutEqSeriesIn = 0;
	  ArithOp_t op = kArithOpUnknown;
	  char inSeriesName[DRMS_MAXSERIESNAMELEN];
	  char withSeriesName[DRMS_MAXSERIESNAMELEN];

	  const char *recSetIn = cmdparams_get_str(&cmdparams, kRecSetIn, NULL);
          int *slicelower = NULL;
          int *sliceupper = NULL;
          int *posout = NULL;

          int slicedim = cmdparams_get_intarr(&cmdparams, kInSliceLower, &slicelower, NULL);
          if (cmdparams_get_intarr(&cmdparams, kInSliceUpper, &sliceupper, NULL) != slicedim)
          {
             error = 1;
             fprintf(stderr, "Slice dimension mismatch.\n");
          }

          if (cmdparams_get_intarr(&cmdparams, kPosOut, &posout, NULL) != slicedim)
          {
             error = 1;
             fprintf(stderr, "Slice dimension mismatch.\n");
          }

          if (!error && slicedim > 0)
          {
             if (slicelower[0] == -1 || sliceupper[0] == -1)
             {
                slicedim = 0;
                free(slicelower);
                slicelower = NULL;
                free(sliceupper);
                sliceupper = NULL;
             }

             if (posout[0] == -1)
             {
                free(posout);
                posout = NULL;
             }
          }

	  const char *withRecSet = cmdparams_get_str(&cmdparams, kWithRecSet, NULL);
	  const char *seriesOut = cmdparams_get_str(&cmdparams,
					      kSeriesOut,
					      NULL); /* actual param string
						      * could be 'NOT SPECIFIED' */
          int pfile = 0; /* is the input series a plain file? */

	  const char *segList = cmdparams_get_str(&cmdparams, kSegList, NULL);
	  const char *opStr = cmdparams_get_str(&cmdparams, kOp, NULL);

	  double bzerov = cmdparams_get_double(&cmdparams, kBzero, NULL);
	  double bscalev = cmdparams_get_double(&cmdparams, kBscale, NULL);
	  DRMS_Type_t dtype = (DRMS_Type_t)cmdparams_get_int(&cmdparams, kDatatype, NULL);
	  if (dtype > DRMS_TYPE_DOUBLE)
	  {
	     dtype = DRMS_TYPE_RAW;
	  }

          int nosegs = cmdparams_isflagset(&cmdparams, kNoSegsFlag);

	  char actualOutputSeries[DRMS_MAXSERIESNAMELEN];
	  /* xxx DO SCALILNG/OFFSET PARAMS LATER */

	  char *segNameArr[kMaxSegs];
	  bzero(segNameArr, sizeof(segNameArr));

	  int nSegs = 0;

	  DRMS_Record_t *inRecTemplate = NULL;
	  DRMS_Record_t *outRecTemplate = NULL;

	  /* op is required */
	  if (strcmp(opStr, kNotSpecified) == 0)
	  {
	       error = 1;
	       fprintf(stderr, "Operation parameters is required.\n");
	  }
	  else
	  {
	       op = MapOp(opStr);

	       if (op == kArithOpUnknown)
	       {
		    error = 1;
		    fprintf(stderr, "Invalid op %s.\n", opStr);
	       }

	       bBinaryOp = IsBinaryOp(op);
	  }

	  /* recSetIn is required */
	  if (!error && strcmp(recSetIn, kNotSpecified) == 0)
	  {
	       error = 1;
	       fprintf(stderr, "Input record set required.\n");
	  }

	  if (!error)
	  {
             if (drms_record_getquerytype(recSetIn) == kRecordSetType_PlainFile)
             {
                pfile = 1;
             }
             else
             {
                snprintf(inSeriesName, sizeof(inSeriesName), recSetIn);
                char *pChar = strchr(inSeriesName, '[');
                if (pChar != NULL)
                {
                   *pChar = '\0';
                }
             }
	  }

	  if (!error)
	  {
             /* Ensure inSeries exists in the database. */
             if (!pfile)
             {
                inRecTemplate = drms_template_record(drms_env, inSeriesName, &status);
                error = (status != 0);
                if (error)
                {
                   fprintf(stderr, "BAILING: inSeries does not exist\n");
                }
             }
             else
             {
                DRMS_RecordSet_t *pfileRS = drms_open_records(drms_env, recSetIn, &status);
                if (pfileRS && pfileRS->n > 0)
                {
                   const char *pfileSeries = pfileRS->records[0]->seriesinfo->seriesname;
                   inRecTemplate = drms_template_record(drms_env, pfileSeries, &status);
                }
                else
                {
                   error = 1;
                   fprintf(stderr, "BAILING: plain file query does not resolve into any fits files.\n");
                }
             }
	  }

	  if (!error && strcmp(seriesOut, kNotSpecified) != 0)
	  {
	       bSeriesOut = 1;

	       /* Does the specified outSeries actually exist? */
	       if (strcmp(inSeriesName, seriesOut) == 0)
	       {
		    bSeriesOutExists = 1; /* because we already verified that inSeries exists. */
		    bSeriesOutEqSeriesIn = 1;
	       }
	       else
	       {
		    outRecTemplate = drms_template_record(drms_env,
							  seriesOut,
							  &status);

		    if (status == DRMS_SUCCESS)
		    {
			 bSeriesOutExists = 1;
		    }
	       }
	  }

	  if (!error && strcmp(withRecSet, kNotSpecified) != 0)
	  {
	       bWithRecSet = 1;
	       snprintf(withSeriesName, sizeof(withSeriesName), withRecSet);
	       char *pChar = strchr(withSeriesName, '[');
	       if (pChar != NULL)
	       {
		    *pChar = '\0';
	       }
	  }

	  if (!error && strcmp(segList, kNotSpecified) != 0)
	  {
	       char *aSeg;

	       for (aSeg = strtok(segList, ","); aSeg && nSegs < kMaxSegs; aSeg = strtok(NULL, ","))
	       {
		    segNameArr[nSegs++] = strdup(aSeg);
	       }
	  }

	  if (!error)
	  {
	       /* main code */
	       DRMSContainer_t inSeriesPrimeKeys;
	       DRMSContainer_t inSeriesSegs;
	       DRMSContainer_t segsToProc;
	       DRMSContainer_t outSeriesPrimeKeys;
	       DRMSContainer_t outSeriesSegs;

               inSeriesPrimeKeys.items = NULL;
               inSeriesSegs.items = NULL;
               segsToProc.items = NULL;
               outSeriesPrimeKeys.items = NULL;
               outSeriesSegs.items = NULL;
               outSeriesSegs.Free = NULL;

	       if (!error)
	       {

                  error = CreateDRMSPrimeKeyContainer(inRecTemplate,
                                                      &inSeriesPrimeKeys,
                                                      ReleaseHContainer);
	       }

	       if (!error)
	       {
                  CreateDRMSSegmentContainer(inRecTemplate, &inSeriesSegs, NULL);
	       }

	       if (!error && bSeriesOut && bSeriesOutExists && !bSeriesOutEqSeriesIn)
	       {
		    error = CreateDRMSPrimeKeyContainer(outRecTemplate,
							&outSeriesPrimeKeys,
							ReleaseHContainer);

		    if (!error)
		    {
			 CreateDRMSSegmentContainer(outRecTemplate,
						    &outSeriesSegs,
						    NULL);
		    }
	       }

	       /* Branch on whether the op is binary or not */
	       if (!error && bBinaryOp)
	       {
		    /* withRecSet is required*/
		    error = (bWithRecSet != 1);

		    DRMSContainer_t withSeriesPrimeKeys;
		    DRMSContainer_t withSeriesSegs;
		    DRMS_Record_t *withRecTemplate = NULL;

		    if (!error)
		    {
			 withRecTemplate = drms_template_record(drms_env,
								withSeriesName,
								&status);
			 error = (status != DRMS_SUCCESS);
		    }

		    if (!error)
		    {
			 error = CreateDRMSPrimeKeyContainer(withRecTemplate,
							     &withSeriesPrimeKeys,
							     ReleaseHContainer);
		    }

		    if (!error)
		    {
			 CreateDRMSSegmentContainer(withRecTemplate,
                                                    &withSeriesSegs,
                                                    NULL);
		    }

		    if (!error)
		    {
			 error = ValidateBinaryOperands(recSetIn,
							segNameArr,
							nSegs,
							&inSeriesPrimeKeys,
							&inSeriesSegs,
							&withSeriesPrimeKeys,
							&withSeriesSegs,
							&segsToProc);
		    }

		    DestroyDRMSContainer(&withSeriesPrimeKeys);
		    DestroyDRMSContainer(&withSeriesSegs);
                    free(withSeriesSegs.items);
	       }
	       else if (!error)
	       {
		    /* Unary op */

		    int nSegsToProc = CreateMatchingSegs(&inSeriesSegs,
							 &inSeriesSegs,
							 &segsToProc);

		    if (nSegsToProc == -1)
		    {
			 error = 1;
		    }
	       }

	       /* Validate output series */
	       if (!error)
	       {
                  if (pfile)
                  {
                     if (!bSeriesOut)
                     {
                        /* No output series specified, and only a pfile input specified */
                        error = 1;
                     }
                     else
                     {
                        strncpy(actualOutputSeries, seriesOut, sizeof(actualOutputSeries));
                     }
                  }
                  else if (bSeriesOut && bSeriesOutExists)
                  {
                     if(!bSeriesOutEqSeriesIn)
                     {
                        /* Ensure that the existing outSeries' prime keys match
                         * those of inSeries. */
                        if (!KeysEqual(&inSeriesPrimeKeys, &outSeriesPrimeKeys))
                        {
                           error = 1;
                        }

                        /* Ensure that the existing outSeries' segments are a superset of
                         * those of segToProc.
                         */
                        if (!error)
                        {
                           int nSegsToProc = CreateMatchingSegs(&outSeriesSegs,
                                                                &segsToProc,
                                                                NULL);

                           if (nSegsToProc == -1)
                           {
                              error = 1;
                           }
                           else if (nSegsToProc > 0 && nSegsToProc != segsToProc.items->num_total)
                           {
                              error = 1;
                              fprintf(stderr, "Target series is missing one or more required segments.\n");
                           }
                        }
                     }

                     strncpy(actualOutputSeries, seriesOut, sizeof(actualOutputSeries));
                  }
                  else if (bSeriesOut && !bSeriesOutExists)
                  {
                     /*
                       error = (DRMS_SUCCESS !=
                       drms_create_seriesfromseries(drms_env,
                       inSeriesName,
                       seriesOut,
                       1, 0, 2, 1,
                       kOutSeriesDesc,
                       1));
                     */
                     error = 1;
                     fprintf(stdout, "Need to call drms_create_seriesfromtemplate()!!!\n");
                     XASSERT(0);

                     if (!error)
                     {
                        strncpy(actualOutputSeries, seriesOut, sizeof(actualOutputSeries));
                     }
                  }
                  else
                  {
                     /* outseries not specified at all */
                     strncpy(actualOutputSeries, inSeriesName, sizeof(actualOutputSeries));
                  }
	       }

	       /* Do operation */
	       if (!error)
	       {
		    if (bBinaryOp == 1)
		    {
			 error = DoBinaryOp(drms_env,
					    op,
					    dtype,
					    &inSeriesPrimeKeys,
					    &segsToProc,
					    recSetIn,
					    withSeriesName,
					    withRecSet,
					    actualOutputSeries,
					    bzerov,
					    bscalev);
		    }
		    else
		    {
                       /* XXX */
                       /* If the caller didn't provide bzero/bscale on the cmd-line, then the output segment's
                        * bzero/bscale values should be assumed. If the caller didn't provide dtype on the cmd-line,
                        * the the output segment's dtype should be assumed. */
                       error = DoUnaryOp(drms_env,
                                         op,
                                         dtype,
                                         &segsToProc,
                                         recSetIn,
                                         slicelower,
                                         sliceupper,
                                         actualOutputSeries,
                                         posout,
                                         bzerov,
                                         bscalev,
                                         nosegs);
		    }

	       } /* Do operation */

	       DestroyDRMSContainer(&inSeriesPrimeKeys);
               DestroyDRMSContainer(&outSeriesPrimeKeys);
	       DestroyDRMSContainer(&inSeriesSegs);
	       DestroyDRMSContainer(&outSeriesSegs);
               free(inSeriesSegs.items);
	       DestroyDRMSContainer(&segsToProc);
	  } /* main code */

          if (slicelower)
          {
             free(slicelower);
          }

          if (sliceupper)
          {
             free(sliceupper);
          }

          if (posout)
          {
             free(posout);
          }
     } /* drms env exists */

     return error;
}
