#include <jsoc.h>
#include <jsoc_main.h>
#include "fstats.h"

 // seg is target segment that must contain a pointer to the owning record.
 // data is the array for which the statistics will be computed.
 // mode is 0 for best method, 1 for quicker estimator of median.

 // set_statistics must be called after any use of drms_copykeys.

int set_statistics(DRMS_Segment_t *seg, DRMS_Array_t *data, int mode)
  {
  int status;
  DRMS_Record_t *rec;
  int i, n, nok, ntot, segnum, nsegs;
  float *fdata;
  double *ddata;
  double min, max, medn, mean, sig, skew, kurt;
  char name[DRMS_MAXNAMELEN];
  char segstr[10];

  if (!seg || !data) return(DRMS_ERROR_INVALIDDATA);
  rec = seg->record;
  if (!rec) return(DRMS_ERROR_INVALIDDATA);
  nsegs = rec->segments.num_total;
  segnum = seg->info->segnum;

  if (nsegs>1)
    sprintf(segstr, "[%d]", segnum);
  else
    segstr[0] = '\0';

  for (n=1, i=0; i<data->naxis; i++)
    n *= data->axis[i];

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
      status = fstats2(n, (float *)data->data, &min, &max, &medn, &mean, &sig, &skew, &kurt, &nok);
    else
      status = fstats(n, (float *)data->data, &min, &max, &medn, &mean, &sig, &skew, &kurt, &nok);
    }
  if (status) return(status);

#define SETKEY(str,val) {				\
  DRMS_Keyword_t *thiskey = drms_keyword_lookup(rec, str, 0); 		\
  if (thiskey) 						\
    {							\
    strcpy(name, str); 					\
    if (thiskey->info->kwflags == kKeywordFlag_PerSegment) 	\
      strcat(name, segstr); 				\
    drms_setkey_double(rec, name, val); 		\
    } 							\
  }
 
  SETKEY("DATAMEAN", mean);
  SETKEY("DATAMIN", min);
  SETKEY("DATAMAX", max);
  SETKEY("DATAMEDN", medn);
  SETKEY("DATARMS", sig);
  SETKEY("DATASKEW", skew);
  SETKEY("DATAKURT", kurt);
  if (drms_keyword_lookup(rec, "DATAVALS", 1))
    ntot = drms_getkey_int(rec, "DATAVALS", 0) + (n - nok);
  else
    ntot = n;
  SETKEY("MISSVALS", (double)(ntot-nok));
  SETKEY("DATAVALS", (double)(nok));
  return(0);
  }
