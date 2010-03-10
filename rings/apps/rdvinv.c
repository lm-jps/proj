/*****************************************************************
 *  rdvinv.c						~rick/src/rdvinv
 *
 *  Responsible: Charles Baldner, Rick Bogart, Sarbani Basu
 *							RBogart@spd.aas.org
 *
 *  Velocity inversion code for ring diagram power spectra fits
 *
 *  Parameters: (type   default         description)
 *	in	string	"Not Specified"	Input data set descriptor
 *	out	string	fort.10.hmixy	Output filename
 *	seg	string	fit.out		Input data segment name
 *	car	int	Not Specified	Carrington rotation number for input
 *				record set? unused
 *	lat	double	Not Specified	Latitude for input record set spec?
 *				unused
 *	lon	double	Not Specified	Longitude for input record set spec?
 *				unused
 *	amu	double	0.005		Error trade-off parameter
 *	ob	double	1.0		Lower frequency limit (mHz)
 *	oe	double	5.2		Upper frequency limit (mHz)
 *	rb	double	0.97		Lower radius limit
 *	re	double	1.00		Upper radius limit
 *	num	int	40		Number of target inversion points
 *	kernel	string	/tmp21/baldner/kernels/ringker.kur_cm_opal_78_mod.ascii
 *				Pathname of file containing inversion kernels
 *	ave	string	Not Specified	Output file for averaging kernels
 *				(if not set, kernels are not written)
 *	coef	string	Not Specified	Output file for ?? coefficients
 *				(if not set, coefficients are not written)
 *
 *  Flags
 *	-v	run verbose
 *	
 *  Notes:
 *    Output is specified by argment 'out'.  If 'out' is a drms data series
 *	that can be opened, the inversions will be  saved to the appropriate
 *	drms records; otherwise, 'out' will be used as the base filename for
 *	the inversion outputs.
 *
 *  Bugs:
 *    If the input data is not in a drms series, an output series can't be
 *	created - rdvinv will exit with an error message.
 *    Arguments 'car', 'lat', and 'lon' are not currently used.
 *    Argument 'out', when a file name, cannot contain any path information.
 *    When writing output to a file, rather than a drms record set, only one
 *	set of files can be specified; if the input contains more than one
 *	data record, each inversion will overwrite the last.
 *
 *  Revision history is at end of file
 ******************************************************************/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jsoc_main.h>
#include "keystuff.c"
#include "rdutil.c"

/* prototypes */
extern void ola_(double *, int *, double *, double *, double *,
      double *, double *, int *, char *, char *, char *, char *, char *,
      int *, int *, int *, double *, double *,
      int *, double *, double *, double *,
      int, int, int, int, int);

char *module_name = "rdvinv";
char *version_id = "0.5";

ModuleArgs_t module_args[] = {
  //{ARG_STRING, "ds", "su_rsb.rdfits", "data set"},
  {ARG_STRING, "seg", "fit.out", "data segment name"},
  {ARG_INT, "car", "Not Specified", "Carrington rotation"},
  {ARG_DOUBLE, "lat", "Not Specified", "AR latitude"},
  {ARG_DOUBLE, "lon", "Not Specified", "AR longitude"},
  {ARG_DOUBLE, "amu", "0.005", "Error trade-off parameter"},
  {ARG_DOUBLE, "ob", "1.0", "Lower frequency limit (mHz)"},
  {ARG_DOUBLE, "oe", "5.2", "Upper frequency limit (mHz)"},
  {ARG_DOUBLE, "rb", "0.97", "Lower radius limit"},
  {ARG_DOUBLE, "re", "1.00", "Upper radius limit"},
  {ARG_INT, "num", "40", "Number of target inversion points"},
  {ARG_STRING, "in", "Not Specified", "Input filename - if set, overides drms lookup"},
  {ARG_STRING, "out", "fort.10.hmixy", "Output filename"},
  {ARG_STRING, "kernel", "/tmp21/baldner/kernels/ringker.kur_cm_opal_78_mod.ascii", ""},
  {ARG_STRING, "ave", "Not Specified", "output file for averaging kernels (if not set, kernels are not written)"},
  {ARG_STRING, "coef", "Not Specified", ""},
  {ARG_FLAG,    "v",    "", "run in verbose mode"},
  {ARG_END}
};

int process_record (const char *filename, int drms_output, char *outfilex,
    char *outfiley, DRMS_Segment_t *osegment_ux, DRMS_Segment_t *osegment_uy,
    char *kernel, char *ave, char *coef, int qave, int qcoef, int verbose,
    double ob, double oe, int num, double rb, double re, double amu) {
  FILE *infile = fopen (filename, "r");
  double *l, *f, *ef, *ux, *eux, *uy, *euy;
  int *n, *mask;
  int npts, i, j, status;
  int lenoutx, lenouty, lenkern, lenave, lencoef;

  status = read_fit_v (infile, &npts, &n, &l, &f, &ef, &ux, &eux, &uy, &euy);
  fclose (infile);
  if(status)	{
    fprintf(stderr, "File %s could not be read\n", filename);
    return 1;
  }
  mask = (int *)malloc (npts * sizeof(int));
  interp_vel (n, l, f, ef, ux, eux, uy, euy, npts);
  autoweed_vel (n, l, ux, uy, mask, npts);
  j = 0;
  for (i=0; i<npts; i++) {
    if (mask[i]) {
      n[j] = n[i];
      l[j] = l[i];
      f[j] = f[i];
      ef[j] = ef[i];
      ux[j] = ux[i];
      eux[j] = eux[i];
      uy[j] = uy[i];
      euy[j] = euy[i];
      j++;
    }
  }
  npts = j;
  if (drms_output) {
    drms_segment_filename (osegment_ux, outfilex);
    lenoutx = strlen (outfilex);
    drms_segment_filename(osegment_uy, outfiley);
    lenouty = strlen (outfiley);
  }
  lenkern = strlen (kernel);
  lenave = (qave) ? strlen (ave) : 0;
  lencoef = (qcoef) ? strlen (coef) : 0;

  ola_ (l, n, f, ux, eux, uy, euy, &npts, kernel, outfilex, outfiley, ave,
      coef, &qave, &qcoef, &verbose, &ob, &oe, &num, &rb, &re, &amu, lenkern,
      lenoutx, lenouty, lenave, lencoef);
  free (n);
  free (mask);
  free (l);
  free (f);
  free (ef);
  free (ux);
  free (eux);
  free (uy);
  free (euy);
  return 0;
}

/* main loop */
int DoIt(void)	{
  CmdParams_t *params = &cmdparams;
  int status = 0;
  int drms_input, drms_output;
  int nrec, rec_i, seg_ct;
  DRMS_RecordSet_t *recordSet = NULL;
  DRMS_Record_t *record, *orecord;
  DRMS_Segment_t *segment, *osegment_ux, *osegment_uy;
  //DRMS_Array_t *keytable_coord;
  char filename[DRMS_MAXPATHLEN+5], *recnum_querry;
  char buffer[1024], outfilex[DRMS_MAXPATHLEN], outfiley[DRMS_MAXPATHLEN], *suffix;
  char *keylist_coord = "recnum"; //"LatHG,LonHG,CarrTime,recnum";
  int *open_record;
  double *data_ptr;
  int *n, *mask, npts, i, j;
  double *l, *f, *ef, *ux, *eux, *uy, *euy;
  int lends, filename_set;

  char *ds = strdup (params_get_str (params, "in"));
  int car = params_get_int (params, "car");
  double lat = params_get_double (params, "lat");
  double lon = params_get_double (params, "lon");
  char *seg = strdup (params_get_str  (params, "seg"));
  char *kernel = strdup (params_get_str (params, "kernel"));
  char *out = strdup (params_get_str (params, "out"));
  char *ave = strdup (params_get_str (params, "ave"));
  char *coef = strdup (params_get_str (params, "coef"));
  double ob = params_get_double (params, "ob");
  double oe = params_get_double (params, "oe");
  double rb = params_get_double (params, "rb");
  double re = params_get_double (params, "re");
  double amu = params_get_double (params, "amu");
  int num = params_get_int (params, "num");
  int verbose = params_isflagset (params, "v");

  int qave = strcmp (ave, "Not Specified");
  int qcoef = strcmp (coef, "Not Specified");

  if(!strcmp(ds, "Not Specified"))	{
    fprintf(stderr, "Error - input must be specified using 'in' argument\n");
    return 1;
  }
  //Attempt to create output drms record
  orecord = drms_create_record (drms_env, out, DRMS_PERMANENT, &status);
  if(orecord)	{	// if record creation was successful
    seg_ct = drms_record_numsegments(orecord);
    if(seg_ct < 1)	{
      fprintf(stderr, "Error: no data segment in output series %s\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    osegment_ux = drms_segment_lookup(orecord, "Ux_inv.out");
    osegment_uy = drms_segment_lookup(orecord, "Uy_inv.out");
    // check that the appropriate data segments exist and are of type generic
    if(!osegment_ux || !osegment_uy)	{
      fprintf(stderr, "Could not find 2 output segments in output series %s\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    if(osegment_ux->info->protocol != DRMS_GENERIC)	{
      fprintf(stderr, "Ux segment is not of type GENERIC in output series %s,\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    if(osegment_uy->info->protocol != DRMS_GENERIC)	{
      fprintf(stderr, "Uy segment is not of type GENERIC in output series %s,\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    drms_output = 1;
    fprintf(stderr, "Using drms series for output\n");
  }
  else	{	// if record creation was not successful, assume `out' is a filename
    sprintf(outfilex, "%s.ux", out);
    sprintf(outfiley, "%s.uy", out);
    drms_output = 0;
    fprintf(stderr, "Not using drms series for output\n");
  }

  recordSet = drms_open_records (drms_env, ds, &status);
  if(status)  {
    //fprintf(stderr, "Call to drms_open_records returned status %i\n", status);
    //return 1;
    fprintf(stderr, "Not using drms series for input\n");
    drms_input = 0;
    nrec = 1;
  }
  else	{
    drms_input = 1;
    nrec = recordSet->n;
    if(nrec<1)	{
      fprintf(stderr, "Number of records opened was 0 or less\n");
      return 1;
    }
  }

  for(rec_i=0; rec_i<nrec; rec_i++)	{
    if(drms_input)	{
      record = recordSet->records[0];
      if (record->sunum != -1LL && record->su == NULL)	{
	record->su = drms_getunit (record->env, record->seriesinfo->seriesname,
	    record->sunum, 1, &status);
      }
      segment = drms_segment_lookup(record, seg);
      sprintf (filename, "%s/S%05i/%s", segment->record->su->sudir, segment->record->slotnum, seg);
    }
    else strcpy(filename, ds);

    if (drms_output) {
      if (rec_i) {
	orecord = drms_create_record (drms_env, out, DRMS_PERMANENT, &status);
	osegment_ux = drms_segment_lookup (orecord, "Ux_inv.out");
	osegment_uy = drms_segment_lookup (orecord, "Uy_inv.out");
	drms_segment_filename (osegment_ux, outfilex);
	suffix = strstr (outfilex, ".generic");
	if(suffix) *suffix = '\0';
	drms_segment_filename (osegment_uy, outfiley);
	suffix = strstr (outfiley, ".generic");
	if (suffix) *suffix = '\0';
      }
    }
    status = process_record (filename, drms_output, outfilex, outfiley, osegment_ux,
	osegment_uy, kernel, ave, coef, qave, qcoef, verbose, ob, oe, num,
	rb, re, amu);
    if(status && drms_output)	{
      drms_close_record(orecord, DRMS_FREE_RECORD);
    }
    else	{
		/* If a drms output segment was opened, write the inversions to it */
      if (drms_output)	{
	int kstat = 0;
	char *key_str;
	double key_dbl;
	int key_int;
	TIME key_time;
	key_str = drms_getkey_string (record, "PrimeKeyString", &status);
	if(!status && check_and_set_key_str (orecord, "PrimeKeyString", key_str)) kstat = 1;
	key_str = drms_getkey_string (record, "CarrTime", &status);
	if(!status && check_and_set_key_str (orecord, "CarrTime", key_str)) kstat = 1;
	key_dbl = drms_getkey_double (record, "LonHG", &status);
	if(!status && check_and_set_key_float (orecord, "LonHG", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "LatHG", &status);
	if(!status && check_and_set_key_float (orecord, "LatHG", key_dbl)) kstat = 1;
	if(check_and_set_key_float(orecord, "amu", amu)) kstat = 1;
	if(check_and_set_key_float(orecord, "freqmin", ob)) kstat = 1;
	if(check_and_set_key_float(orecord, "freqmax", oe)) kstat = 1;
	if(check_and_set_key_float(orecord, "radmin", rb)) kstat = 1;
	if(check_and_set_key_float(orecord, "radmax", re)) kstat = 1;
	key_time = drms_getkey_time (record, "MidTime", &status);
	if(!status && check_and_set_key_time (orecord, "MidTime", key_time)) kstat = 1;
	key_int = drms_getkey_int (record, "Duration", &status);
	if(!status && check_and_set_key_int (orecord, "Duration", key_int)) kstat = 1;
	key_str = drms_getkey_string (record, "MapProj", &status);
	if(!status && check_and_set_key_str (orecord, "MapProj", key_str)) kstat = 1;
	key_dbl = drms_getkey_double (record, "MapScale", &status);
	if(!status && check_and_set_key_float (orecord, "MapScale", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "Map_PA", &status);
	if(!status && check_and_set_key_float (orecord, "Map_PA", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "Width", &status);
	if(!status && check_and_set_key_float (orecord, "Width", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "Height", &status);
	if(!status && check_and_set_key_float (orecord, "Height", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "ZonalTrk", &status);
	if(!status && check_and_set_key_float (orecord, "ZonalTrk", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "MeridTrk", &status);
	if(!status && check_and_set_key_float (orecord, "MeridTrk", key_dbl)) kstat = 1;
	if(check_and_set_key_str(orecord, "Module", module_name)) kstat = 1;
	key_str = drms_getkey_string(record, "Module", &status);
	if(check_and_set_key_str(orecord, "Source", key_str)) kstat = 1;
	drms_close_record(orecord, DRMS_INSERT_RECORD);
      }
    }
//   drms_close_records(recordSet, DRMS_FREE_RECORD);
  }

  return 0;
}

/*
 *  Revision History
 *
 *  June 2009 - v 0.1 module created, no real drms functionality
 *  09.10.2009 - v 0.2 drms functionality added - input and output from drms records now 
 *  	possible, but input/output to rooted filenames also possible
 *  0.3 2009.10.22	fixed format of name of output files
 *  0.4 2010.02.11	unified processing calls for named file and drms input
 *	branches
 *  0.5 2010.02.12	Combined branches into one, added some very obvious
 *	error handling (no given input argument, input cannot be opened)
 *			removed rdutil.h inclusion
 *  v 0.5 frozen 2010.03.08
 *  
 *
 */
