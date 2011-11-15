#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/libs/stats/set_statistics.c,v 1.4 2011/11/15 20:20:38 phil Exp $"

#include <jsoc.h>
#include <jsoc_main.h>
#include "fstats.h"

 // seg is target segment that must contain a pointer to the owning record.
 // data is the array for which the statistics will be computed.
 // mode is 0 for best method, 1 for quicker estimator of median.

 // set_statistics must be called after any use of drms_copykeys.

 // The return value is 0 for OK and the return value
 // from fstats or dstats or fstats2 or dstats2 whichever was called by
 // set_statistics depending on mode and the data type.

 // The set_statistics function will compute AND populate keywords for:
 //   DATAMEAN, DATAMIN, DATAMAX, DATAMEDN, DATARMS, DATASKEW, DATAKURT,
 //   DATAVALS and MISSVALS if they are present.

 // If the keyword scope is per segment, then the appropriate segment specific
 // keywords will be set.

#define SETKEY(str,val) 				\
  {							\
  if (thiskey = drms_keyword_lookup(rec, str, 0))	\
    drms_setkey_double(rec, str, val);  		\
  else							\
    {							\
    char name[DRMS_MAXNAMELEN];				\
    strcpy(name, str);	 				\
	    strcat(name, segstr); 				\
    if ((thiskey = drms_keyword_lookup(rec, name, 0)) && ((thiskey->info->kwflags & 1) == kKeywordFlag_PerSegment))      \
        drms_setkey_double(rec, name, val); 		\
    }                                                 \
  }

int set_statistics(DRMS_Segment_t *seg, DRMS_Array_t *data, int mode)
  {
  int status;
  DRMS_Record_t *rec;
  int i, n, nok, ntot, segnum, nsegs;
  double min, max, medn, mean, sig, skew, kurt;
  char segstr[10];
  DRMS_Keyword_t *thiskey;

  if (!seg || !data)
    return(DRMS_ERROR_INVALIDDATA);
  rec = seg->record;
  if (!rec)
    return(DRMS_ERROR_INVALIDDATA);
  nsegs = rec->segments.num_total;
  segnum = seg->info->segnum;

  sprintf(segstr, "_%03d", segnum);

  for (n=1, i=0; i<data->naxis; i++)
    n *= data->axis[i];

  ntot = drms_getkey_int(rec, "TOTVALS", &status);
  if (status)
    {
    char name[DRMS_MAXNAMELEN];
    strcpy(name, "TOTVALS");
    strcat(name, segstr);
    if ((thiskey = drms_keyword_lookup(rec, name, 0)) && ((thiskey->info->kwflags & 1) == kKeywordFlag_PerSegment))
        ntot = drms_getkey_int(rec, name, &status);
    if (status) 
      ntot = n;
    }

  if (data->type == DRMS_TYPE_FLOAT)
    {
    if (mode)
      status = fstats2(n, (float *)data->data, &min, &max, &medn, &mean, &sig, &skew, &kurt, &nok);
    else
      status = fstats(n, (float *)data->data, &min, &max, &medn, &mean, &sig, &skew, &kurt, &nok);
    }
  else if (data->type == DRMS_TYPE_DOUBLE)
    {
    if (mode)
      status = dstats2(n, (double *)data->data, &min, &max, &medn, &mean, &sig, &skew, &kurt, &nok);
    else
      status = dstats(n, (double *)data->data, &min, &max, &medn, &mean, &sig, &skew, &kurt, &nok);
    }
  else return(-1);

  if (status)
    return(status);

  SETKEY("DATAMEAN", mean);
  SETKEY("DATAMIN", min);
  SETKEY("DATAMAX", max);
  SETKEY("DATAMEDN", medn);
  SETKEY("DATARMS", sig);
  SETKEY("DATASKEW", skew);
  SETKEY("DATAKURT", kurt);
  SETKEY("DATAVALS", (double)nok);
  SETKEY("MISSVALS", (double)(ntot - nok));
  return(0);
  }
