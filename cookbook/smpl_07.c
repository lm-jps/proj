/*
 *  smpl_07.c						$DRMS/proj/cookbook/
 *
 *  Populates a data series consisting of variable size images with random
 *    data from a selected distribution function
 *  Illustrates features of the DRMS_Array and DRMS_Segment structs, and
 *    writing to DRMS record segments
 *
 *  Usage:
 *    smpl_07 [ds= cols= rows= ]
 *
 *  Revision history is at end of file.
 */

char *module_name = "CookbookRecipe:07";
char *version_id = "1.0";

#include <jsoc_main.h>

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"ds",   "drms.images", "name of data series"},
  {ARG_INT,	"cols", "0",  "number of columns in image (default (rows)"},
  {ARG_INT,	"rows", "0",  "number of rows in image (default: cols)"},
  {ARG_NUME,	"dist", "uniform", "distribution for pseudo-randoms",
      "uniform, dblexp, normal"},
  {ARG_INT,	"seed", "0",  "initial seed for pseudo-random number generator"},
  {ARG_FLAG, "v", "", "verbose flag"}, 
  {ARG_END}
};

enum dist_type {UNIFORM, DBL_EXP, NORMAL};

float next_rand (int distrib) {
  double val0, val1, rmod, gf;
  static double sav;
  static int avail = 0;

  switch (distrib) {
    case (UNIFORM):
      return 2 * drand48 () - 1.0;
    case (DBL_EXP):
      val0 = -log (drand48 ());
      val1 = (drand48 () < 0.5) ? -0.5 : +0.5;
      return val0 * val1;
    case (NORMAL):
      if (avail) {
	avail = 0;
	return sav;
      }
      rmod = 2.0;
      while (rmod >= 1.0) {
	val0 = 2 * drand48 () - 1.0;
	val1 = 2 * drand48 () - 1.0;
	rmod = val0 * val0 + val1 * val1;
      }
      gf = sqrt (-2.0 * log (rmod) / rmod);
      sav = gf * val0;
      avail = 1;
      return gf * val1;
    default:
      return 0.0;
  }
}

void rand_fill_array (float *data, int cols, int rows, int df) {
  double xr;
  int n, col, row, pr;
  int colsM1 = cols - 1;
  double kf1 = 1.0, kf2 = 0.5, kf3 = 1.0 / 3.0, kf4 = 0.25, kf5 = 0.2;

  n = row = 0;
  data[n++] = next_rand (df);
  for (col = 1; col < colsM1; col++, n++) {
    data[n] = kf1 * data[n-1] + next_rand (df);
  }
  data[n++] = kf2 * (data[n-1] + data[n+1-cols]) + next_rand (df);

  row++;
  for (; row < rows; row++) {
    pr = n - cols;
    data[n++] = kf3 * (data[n-1] + data[pr] + data[pr+1]) + next_rand (df);
    pr++;
    for (col = 1; col < colsM1; col++, n++, pr++) {
      data[n] = kf4 * (data[pr-1] + data[pr] + data[pr+1] + data[n-1]) +
	  next_rand (df);
    }
    data[n++] = kf5 * (data[pr-1] + data[pr] + data[pr+1] + data[n-1] + data[n+1-cols]) +
	next_rand (df);
  }
}

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_Record_t *record;
  DRMS_Segment_t *record_segment;
  DRMS_Array_t *data_array;
  float *data;
  int axes[2];
  int status = 0;

  char *series = params_get_str (params, "ds");
  int cols = params_get_int (params, "cols");
  int rows = params_get_int (params, "rows");
  int distrib = params_get_int (params, "dist");
  long int seed = params_get_int (params, "seed");
  int verbose = params_isflagset (params, "v");
	/*  make sure that at least one of cols and rows has been specified  */
  if (cols <= 0) cols = rows;
  if (rows <= 0) rows = cols;
  if (rows <= 0) {
    fprintf (stderr, "Error: at least one of cols and rows must be specified\n");
    return 1;
  }
					   /*  create the output data array  */
  axes[0] = cols;
  axes[1] = rows;
  data = (float *)malloc (cols * rows * sizeof (float));
  data_array = drms_array_create (DRMS_TYPE_FLOAT, 2, axes, (void *)data,
      &status);
			       /*  fill the output array (with random data)  */
  srand48 (seed);
  rand_fill_array (data, cols, rows, distrib);

  if (verbose) {
    int col, row, n;
    int mincol, maxcol, minrow, maxrow;
    float minval, maxval;

    minval = maxval = data[0];
    mincol = minrow = maxcol = maxrow = 0;
    for (n = 0, row = 0; row < rows; row++) {
      for (col = 0; col < cols; col++, n++) {
	if (data[n] < minval) {
	  minval = data[n];
	  mincol = col;
	  minrow = row;
	}
	if (data[n] > maxval) {
	  maxval = data[n];
	  maxcol = col;
	  maxrow = row;
	}
      }
    }
    printf ("min @ [%d, %d] : %f\n", mincol, minrow, minval);
    printf ("max @ [%d, %d] : %f\n", maxcol, maxrow, maxval);
  }

  record = drms_create_record (drms_env, series, DRMS_PERMANENT, &status);
  if (!record) {
    fprintf (stderr, "Error: Unable to create record in series %s\n", series);
    fprintf (stderr, "       Does series exist and is it writeable by you?\n");
    return 1;
  }
  record_segment = drms_segment_lookupnum (record, 0);
			/*  set the output scaling parameters for the array  */
  data_array->bscale = 1.0;
  data_array->bzero = 0.0;
  status = drms_segment_write (record_segment, data_array, 0);
  if (status) {
    fprintf (stderr, "Error writing data to segment 0 of series %s\n", series);
    fprintf (stderr, "      series may not have appropriate structure\n");
    drms_free_array (data_array);
    drms_close_record (record, DRMS_FREE_RECORD);
    return 1;
  }
  drms_close_record (record, DRMS_INSERT_RECORD);
  return 0;
}
/*
 *  Revision History
 *
 *  09.08.02	file created by R Bogart
 */
