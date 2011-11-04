/*
 *  smpl_08.c						$DRMS/proj/cookbook/
 *
 *  Calculates primary statistics for the selected segment(s) of selected
 *    records in a DRMS data series
 *
 *  Populates a data series consisting of variable size images with random
 *    data from a selected distribution function
 *  Illustrates the use of drms_segment_read() and its returned DRMS_Array
 *    struct to manipulate data in the data "cubes" of individual record
 *    segments
 *  The default record set is the (current) last record in the series
 *    created and populated in Recipe 07
 *
 *  Usage:
 *    smpl_08 [ds= ]
 *
 *  Revision history is at end of file.
 */

char *module_name = "CookbookRecipe:08";
char *version_id = "1.0";

/*
 *  Bugs:
 *    drms_segment_read seems to fail if storage protocol is bin (or bin.gz)
 *
 */

					/*  required include, declarations  */
#include <jsoc_main.h>
#include <math.h>

ModuleArgs_t module_args[] = {
  {ARG_STRING, "ds", "drms.images[:#$]", "data set"},
  {ARG_STRING, "seg", "Not Specified", "data segment name (default: all)"},
  {}
};

void calc_and_print_stats (DRMS_Segment_t *record_segment) {
  DRMS_Array_t *data_array;		       /*  array struct from segment  */
  double *data;
  double sum, sum2, avg, stdev;
  long long ntot, nv;
  int n, rank;
  int status;
			/*  read data from segment into array as doubles
				        (will block until segment is staged)  */
  data_array = drms_segment_read (record_segment, DRMS_TYPE_DOUBLE, &status);
  if (status) {
/*      if drms_segment_read() returns a non-0 status, it probably could
	not find the file, which may have aged out of the archive, or may
	not be accessible from the current system, or the segment may have an
	unsupported protocol. In any case, there's nothing to be done...      */
    printf ("Not found\n");
    return;
  }
  rank = data_array->naxis;
  ntot = 1;
  for (n = 0; n < rank; n++) ntot *= data_array->axis[n];
  data = (double *)data_array->data;
  nv = 0;
  sum = sum2 = 0;
  for (n = 0; n < ntot; n++) {
    if (isnan (data[n])) continue;
    nv++;
    sum += data[n];
    sum2 += data[n] * data[n];
  }
  printf ("%ld of %ld valid", nv, ntot);
  if (nv) {
    avg = sum / nv;
    sum2 /= nv;
    stdev = sqrt (sum2 - avg * avg);
    printf ("; mean = %11.4e std dev = %11.4e", avg, stdev);
  }
  printf ("\n");
  drms_free_array (data_array);
}
							/*  main module body  */
int DoIt () {
  DRMS_RecordSet_t *drs;					/*  data set  */
  DRMS_Record_t *record;   	       /*  individual record of the data set  */
  DRMS_Segment_t *record_segment;	 /*  single data segment of a record  */
  int recn, rec_ct, segn, seg_ct;
  int n, seg_selected;
  int status;
						/* Get command line arguments */
  char *dsspec = strdup (params_get_str (&cmdparams, "ds"));
  char *seg_name = strdup (params_get_str (&cmdparams, "seg"));

  seg_selected = strcmp (seg_name, "Not Specified");
  drs = drms_open_records (drms_env, dsspec, &status);
  if (!drs) {
    fprintf (stderr, "Error: unable to open record set %s\n", dsspec);
    return 0;
  }
  rec_ct = drs->n;
  if (rec_ct < 1) {
    fprintf (stderr, "No records in selected set %s\n", dsspec);
    return 0;
  }
		     /*  loop over all the records in the selected data set  */
  for (recn = 0; recn < rec_ct; recn++) {
    if (rec_ct > 1) printf ("record %d:\n", recn);
    record = drs->records[recn];
    seg_ct = drms_record_numsegments (record);
    if (seg_ct < 1) printf ("  no segments in selected record\n");
    if (seg_selected) {
      record_segment = drms_segment_lookup (record, seg_name);
      if (!record_segment) {
	fprintf (stderr, "Error, unable to locate segment %s of record %d\n",
	    seg_name, recn);
	continue;
      }
      calc_and_print_stats (record_segment);
    } else {
			      /*  loop over all segments defined for series  */
      for (segn = 0; segn < seg_ct; segn++) {
	record_segment = drms_segment_lookupnum (record, segn);
	if (!record_segment) {
	  fprintf (stderr, "Error, unable to locate segment %d of record %d\n",
	      segn, recn);
	  continue;
	}
	printf ("  %s\t", record_segment->info->name);
	if (record_segment->info->protocol == DRMS_PROTOCOL_INVALID ||
	    record_segment->info->protocol == DRMS_GENERIC) {
	  printf ("unsupported data protocol; skipped\n");
	  continue;
	}
        calc_and_print_stats (record_segment);
      }
    }
  }

  return 0;
}

/*
 *  Revision History
 *
 *  08.07.29	Rick Bogart	created this file as a heuristic exercise
 *  09.09.20	"		incorporated in cookbook with minor mods
 */
