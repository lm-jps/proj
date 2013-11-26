/*
 *  invert_td_hr.c						~rick/hmi/timed
 *
 *  Invert travel times from Gabor wavelet and/or Gizon-Birch fitting for
 *    sound-speed and velocity using multichannel deconvolution
 *    
 *  Responsible:  Junwei Zhao & Rick Bogart	RBogart@solar.stanford.edu
 *
 *  Usage:
 *    invert_td_hr [-v] in= kernels= out=
 *
 *  Arguments: (type		default         description)
 *	in	DataSet
 *      copy	String		"+"	Comma separated list of keys to
 *			propagate forward; if the list includes "+", all
 *			keys in the default list of propagation keys, plus
 *			all prime keys in the input records will be copied
 *
 *  Bugs:
 *    Output array sizes are hardcoded and fixed; no idea what would happen if
 *	input sizes did not match, nor checking to make sure they do
 *    Some assumed key information is hardcoded in Fortran subroutines (e.g. vv)
 *    Only a single input record is processed, though code is very fast,
 *	requiring kernel files to be read repeatedly
 *    Dual processing of travel-time segments from both fitting methods is not
 *	implemented
 *    No documentation
 *
 *  Notes:
 *    All the work is done by the FORTRAN subroutines inver_cs_sub_() &
 *	inver_v_sub_(), whose arguments are:
 *    inver_cs_sub_()
 *	float *oi_mn	input mean out - in (?) travel times
 *	kern *kern_file
 *	    name of a FORTRAN unformatted file containing the following vectors
 *	    of length n1 * n2 * num_annuli:
 *	float *cs	output sound speed values
 *    inver_v_sub_()
 *	float *oi	input out - in (?) travel times
 *	float *we	input west - east (?) travel times
 *	float *ns	input north - south (?) travel times
 *	float *lonc	input Stonyhurst longitude of center of region
 *	float *latc	input latitude of center of region
 *	kern *kern_file
 *	    name of a FORTRAN unformatted file containing the following vectors
 *	    of length n1 * n2 * num_annuli:
 *	float *vx	output V-west (?) values
 *	float *vy	output V-north (?) values
 *	float *vz	output V-outward (?) values
 *	int *kerntyp	input flag for approximation used in kernel calculation;
 *	    currently unused
 *
 *  Revision history is at end of file.
 */

#include <jsoc_main.h>
#include "keystuff.c"

#define KERNTYP_RAY	(1)
#define KERNTYP_BORN	(2)
                                                      /*   module identifier  */
char *module_name = "timed_invert";
char *version_id = "0.81";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "in", "", "input data record"},
  {ARG_STRING, "kernels", "", "input kernel"},
  {ARG_STRING, "out", "", "output series"}, 
  {ARG_NUME, "fitting", "both", "travel-time fitting method used",
      "Gabor, GizonBirch, both"},      
  {ARG_STRING, "copy",  "+",
      "comma separated list of keys to propagate forward"},
  {ARG_FLAG, "v", "0", "verbose output"},      
  {}
};

        /*  list of keywords to propagate (if possible) from input to output  */
char *propagate[] = {"CarrRot", "CMLon", "LonHG", "LatHG", "LonCM", "MidTime",
    "PxLocVer", "WCSAXES", "WCSNAME", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2",
    "CDELT1", "CDELT2", "CUNIT1", "CUNIT2", "CRVAL1", "CRVAL2",
    "MapScale", "MapProj", "Map_PA"};

extern void inver_cs_sub_ (float *oi_mn, float *kern_cs, float *cs);
extern void inver_v_sub_ (float *oi, float *we, float *ns, float *lonc,
    float *latc, float *kern_vx, float *kern_vy, float *kern_vz,
    float *vx, float *vy, float *vz, int *kerntyp);

enum fittypes {GABOR_WAVELET, GIZON_BIRCH, ALL_FITS};

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *ids, *kds;
  DRMS_Record_t *ttrec, *krec, *orec;
  DRMS_Segment_t *ttseg, *seg_cs, *seg_vx, *seg_vy, *seg_vz;
  DRMS_Segment_t *kcs_seg, *kvx_seg, *kvy_seg, *kvz_seg;
  DRMS_Array_t *kcs_arr, *kvx_arr, *kvy_arr, *kvz_arr;
  DRMS_Array_t *ttdat, *cs_res, *vx_res, *vy_res, *vz_res;
  float *kcs, *kvx, *kvy, *kvz;
  float *input, *oi_mn, *oi, *we, *ns, *cs, *vx, *vy, *vz;
  float *inpmn, *inpoi, *inpwe, *inpns;
  int *axes;
  int col, cols, row, rows, pln, pln0, pln1, pln2, pln3;
  int i, j, k, n, ntot, rgn, rgnct, outok;
  int pxlocver, kerntyp;
  int kstat, key_n, propct, propstd, status;
  char **copykeylist;
  char *type_segname[] = {"gabor", "GB"};
  char *fitting_key[] = {"Gabor", "GizonBirch", "both"};
  char *sval;
  char recid[DRMS_MAXQUERYLEN], source[DRMS_MAXQUERYLEN];
  char module_ident[64];

  int keyct = sizeof (propagate) / sizeof (char *);
  int docsinv = 1, dovinv = 1;

  char *ttimds = strdup (params_get_str (params, "in"));
  char *kernds = strdup (params_get_str (params, "kernels"));
  char *outser =  strdup (params_get_str (params, "out"));
  char *propagate_req = strdup (params_get_str (params, "copy"));
  int fitopt = params_get_int (params, "fitting");
  int verbose = params_isflagset (params, "v");

  snprintf (module_ident, 64, "%s v %s", module_name, version_id);
  if (verbose) printf ("%s:\n", module_ident);
						 /*  open kernel data record  */
  kds = drms_open_records (drms_env, kernds, &status);
  if (status) {
    fprintf (stderr, "Error: unable to open input record set %s\n", kernds);
    return 1;
  }
  if (kds->n != 1) {
    fprintf (stderr,
	"Error: kernel record set %s does not contain single record\n",
	kernds);
    drms_close_records (kds, DRMS_FREE_RECORD);
    return 1;
  }
  krec = kds->records[0];
  kcs_seg = drms_segment_lookup (krec, "cs");
  kvx_seg = drms_segment_lookup (krec, "vx");
  kvy_seg = drms_segment_lookup (krec, "vy");
  kvz_seg = drms_segment_lookup (krec, "vz");
  sval = drms_getkey_string (krec, "Approx", &status);
  if (status) {
    fprintf (stderr, "Warning: kernel approximation type not specified in\n");
    fprintf (stderr, "         %s ; Ray path assumed\n", kernds);
    kerntyp = KERNTYP_RAY;
  } else {
    if (!strcasecmp (sval, "Raypath")) kerntyp = KERNTYP_RAY;
    else if (!strcasecmp (sval, "Born")) kerntyp = KERNTYP_BORN;
    else {
      fprintf (stderr, "Warning: kernel approximation type %s not recognized in\n",
	  sval);
      fprintf (stderr, "         %s ; Ray path assumed\n", kernds);
      kerntyp = KERNTYP_RAY;
   }
  }
							/*  read the kernels  */
  if (kcs_seg) {
    kcs_arr = drms_segment_read (kcs_seg, DRMS_TYPE_FLOAT, &status);
    if (status) {
      fprintf (stderr, "Warning: unable to read \"cs\" segment in kernel record;\n");
      fprintf (stderr, "         sound-speed inversions will not be attempted\n");
      docsinv = 0;
    } else
    kcs = (float *)kcs_arr->data;
  } else {
    fprintf (stderr, "Warning: no \"cs\" segment in kernel record;\n");
    fprintf (stderr, "         sound-speed inversions will not be attempted\n");
    docsinv = 0;
  }

  if (kvx_seg && kvy_seg && kvz_seg) {
    kvx_arr = drms_segment_read (kvx_seg, DRMS_TYPE_FLOAT, &status);
    if (status) dovinv = 0;
    kvy_arr = drms_segment_read (kvy_seg, DRMS_TYPE_FLOAT, &status);
    if (status) dovinv = 0;
    kvz_arr = drms_segment_read (kvz_seg, DRMS_TYPE_FLOAT, &status);
    if (status) dovinv = 0;
    if (!dovinv) {
      fprintf (stderr,
	  "Warning: unable to read one or more v_i segments in kernel record;\n");
      fprintf (stderr, "         velocity inversions will not be attempted.\n");
      if (!docsinv) {
	fprintf (stderr, "         nothing to do!\n");
	drms_close_records (kds, DRMS_FREE_RECORD);
	return 1;
      }
    }
    kvx = (float *)kvx_arr->data;
    kvy = (float *)kvy_arr->data;
    kvz = (float *)kvz_arr->data;    
  } else {
    fprintf (stderr,
	"Warning: one or more v_i segments missing in kernel record;\n");
    fprintf (stderr, "         velocity inversions will not be attempted.\n");
    if (!docsinv) {
      fprintf (stderr, "         nothing to do!\n");
      drms_close_records (kds, DRMS_FREE_RECORD);
      return 1;
    }
    dovinv = 0;
  }
 						  /*  open input data series  */
  ids = drms_open_records (drms_env, ttimds, &status);
  if (status) {
    fprintf (stderr, "Error: unable to open input record set %s\n", ttimds);
    drms_close_records (kds, DRMS_FREE_RECORD);
    return 1;
fprintf (stderr, "input data set is open\n");
  }
  rgnct = ids->n;
  if (rgnct < 1) {
    fprintf (stderr, "Error: input dataset %s contains no records\n", ttimds);
    drms_close_records (kds, DRMS_FREE_RECORD);
    drms_close_records (ids, DRMS_FREE_RECORD);
    return 1;
  }
					/*  allocate input and output arrays  */
  axes = (int *)malloc (3 * sizeof (int));
  cols = axes[0] = rows = axes[1] = 256;
  axes[2] = 11;
  ntot = 1;
  for (n = 0; n < 3; n++) ntot *= axes[n];
  oi_mn = (float *)malloc (ntot * sizeof (float));
  oi = (float *)malloc (ntot * sizeof (float));
  we = (float *)malloc (ntot * sizeof (float));
  ns = (float *)malloc (ntot * sizeof (float));
  vx = (float *)malloc (ntot * sizeof (float));
  vy = (float *)malloc (ntot * sizeof (float));
  vz = (float *)malloc (ntot * sizeof (float));
  if (docsinv) {
    cs = (float *)malloc (ntot * sizeof (float));
    cs_res = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, cs, &status);
    if (status) {
      fprintf (stderr, "Error: unable to create required ouput array cs\n");
      return 1;
    }
  }

  vx_res = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, vx, &status);
  vy_res = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, vy, &status);
  vz_res = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, vz, &status);

  outok = 0;
  rgn = 0;
		    /*  read travel times from input record segment "gabor"  */
  for (rgn = 0; rgn < rgnct; rgn++) {
    if (verbose) printf ("processing record %d\n", rgn);
    ttrec = ids->records[rgn];
    drms_sprint_rec_query (source, ttrec);
		   /*  check PxLocVer, must be consistent with code version  */
    pxlocver = drms_getkey_int (ttrec, "PxLocVer", &status);
    if (pxlocver < 1 && verbose) {
      printf ("Warning: Pixel locator file used for travel time map in record\n");
      printf ("         %s has unknown version\n", source);
    }
    if (fitopt == ALL_FITS) {
      fprintf (stderr, "Warning: inversions from dual fits not supported\n");
      fprintf (stderr, "         only Gabor wavelet fits will be used\n");
      fitopt = GABOR_WAVELET;
    }
    ttseg = drms_segment_lookup (ttrec, type_segname[fitopt]);
    if (!ttseg) {
      fprintf (stderr, "Error: no segment \"%s\" in record %s; skipped\n",
	  type_segname[fitopt], source);
      continue;
    }
    ttdat = drms_segment_read (ttseg, DRMS_TYPE_FLOAT, &status);
    if (!ttdat) {
      fprintf (stderr, "Error: unable to read input data record-segment in %s\n",
          source);
      continue;
    }
    input = (float *)ttdat->data;
						/*  create an output record  */
    orec = drms_create_record (drms_env, outser, DRMS_PERMANENT, &status);
    if (!orec) {
      fprintf (stderr, "Error: unable to create record %d in series %s\n", rgn,
	  outser);
      fprintf (stderr, "       remainder of processing skipped\n");
      drms_close_records (ids, DRMS_FREE_RECORD);
      return 0;
    }
    seg_cs = drms_segment_lookup (orec, "cs");
    seg_vx = drms_segment_lookup (orec, "vx");
    seg_vy = drms_segment_lookup (orec, "vy");
    seg_vz = drms_segment_lookup (orec, "vz");
/*  check the output data series struct (only needs to be done once)  */
    if (!outok) {
				     /*  check for the appropriate segments  */
      if (!seg_cs) {
	docsinv = 0;
	if (!dovinv) {
	  fprintf (stderr, "Error: output data series %s\n", outser);
	  fprintf (stderr,
	      "       lacks segment \"cs\", and no flow inversion kernels available\n");
	  drms_close_records (ids, DRMS_FREE_RECORD);
	  return 1;
	}
	fprintf (stderr,
	    "Warning: cs segment missing in output record;\n");
	fprintf (stderr, "         sound-speed inversions will not be attempted.\n");
      }
      if (!seg_vx) {
        dovinv = 0;
	if (!docsinv) {
	  fprintf (stderr, "Error: output data series %s\n", outser);
	  fprintf (stderr,
	      "       lacks one or more segments \"vx\", \"vy\", \"vz\"\n");
	  fprintf (stderr,
	      "       and no sound-speed inversion kernels available\n");
	  drms_close_records (ids, DRMS_FREE_RECORD);
	  return 1;
	}
	fprintf (stderr,
	    "Warning: one or more v_i segments missing in output record;\n");
	fprintf (stderr, "         velocity inversions will not be attempted.\n");
      }
     /*  check structure of segments: Fortran code has hard-coded parameters  */
      if (dovinv) {
        if (seg_vx->info->naxis != 3 || seg_vy->info->naxis != 3 ||
	    seg_vz->info->naxis != 3) {
	  fprintf (stderr, "Error: output data series %s\n", outser);
	  fprintf (stderr,
	      "       has wrong rank for one or more segments \"vx\", \"vy\", \"vz\"\n");
	  drms_close_records (ids, DRMS_FREE_RECORD);
	  return 1;
	}
        if (seg_vx->info->scope != DRMS_VARDIM) {
	  if (seg_vx->axis[0] != 256 || seg_vx->axis[1] != 256 ||
	      seg_vx->axis[2] != 11) {
	    fprintf (stderr, "Error: in output data series %s\n", outser);
	    fprintf (stderr,
		"       dimensions for segment \"vx\" differ from hard-coded values\n");
	    drms_close_records (ids, DRMS_FREE_RECORD);
	    return 1;
	  }
	}
        if (seg_vy->info->scope != DRMS_VARDIM) {
	  if (seg_vy->axis[0] != 256 || seg_vy->axis[1] != 256 ||
	      seg_vy->axis[2] != 11) {
	    fprintf (stderr, "Error: in output data series %s\n", outser);
	    fprintf (stderr,
		"       dimensions for segment \"vy\" differ from hard-coded values\n");
	    drms_close_records (ids, DRMS_FREE_RECORD);
	    return 1;
	  }
	}
        if (seg_vz->info->scope != DRMS_VARDIM) {
	  if (seg_vz->axis[0] != 256 || seg_vz->axis[1] != 256 ||
	      seg_vz->axis[2] != 11) {
	    fprintf (stderr, "Error: in output data series %s\n", outser);
	    fprintf (stderr,
		"       dimensions for segment \"vx\" differ from hard-coded values\n");
	    drms_close_records (ids, DRMS_FREE_RECORD);
	    return 1;
	  }
	}
      }
      if (docsinv) {
        if (seg_cs->info->naxis != 3) {
	  fprintf (stderr, "Error: output data series %s\n", outser);
	  fprintf (stderr,
	      "       has wrong rank for segment \"cs\"\n");
	  drms_close_records (ids, DRMS_FREE_RECORD);
	  return 1;
	}
        if (seg_cs->info->scope != DRMS_VARDIM) {
	  if (seg_cs->axis[0] != 256 || seg_cs->axis[1] != 256 ||
	      seg_cs->axis[2] != 11) {
	    fprintf (stderr, "Error: in output data series %s\n", outser);
	    fprintf (stderr,
		"       dimensions for segment \"cs\" differ from hard-coded values\n");
	    drms_close_records (ids, DRMS_FREE_RECORD);
	    return 1;
	  }
	}
      }
      outok = 1;
    }
	  /*  check input dimensions: Fortran code has hard-coded parameters  */

			   /*  copy designated keywords from input to output  */
    kstat = 0;
    propct = construct_stringlist (propagate_req, ',', &copykeylist);
    propstd = 0;
					/*  copy specifically requested keys  */
    for (key_n = 0; key_n < propct; key_n++) {
      if (strcmp (copykeylist[key_n], "+"))
        kstat += check_and_copy_key (orec, ttrec, copykeylist[key_n]);
      else propstd = 1;
    }
    if (propstd) {
						      /*  copy standard keys  */
      kstat = propagate_keys (orec, ttrec, propagate, keyct);
					   /*  and any additional prime keys  */
      kstat += copy_prime_keys (orec, ttrec);
    }
    kstat += check_and_set_key_str (orec, "Module", module_ident);
    kstat += check_and_set_key_str (orec, "BLD_VERS", jsoc_version);
    kstat += check_and_set_key_time (orec, "Created", CURRENT_SYSTEM_TIME);
    kstat += check_and_set_key_str (orec, "Source", source);
    kstat += check_and_set_key_str (orec, "Input", ttimds);
    kstat += check_and_set_key_str (orec, "Kernels", kernds);
    if (kerntyp == KERNTYP_RAY)
      kstat += check_and_set_key_str (orec, "KrnAppr", "Raypath");
    else if (kerntyp == KERNTYP_BORN)
      kstat += check_and_set_key_str (orec, "KrnAppr", "Born");
    kstat += check_and_set_key_str (orec, "Fitting", fitting_key[fitopt]);
    if (kstat) {
      fprintf (stderr, "Error writing key value(s) to record in series %s:\n",
	  outser);
      fprintf (stderr, "      series may not have appropriate structure;\n");
      fprintf (stderr, "      record %d skipped\n", rgn);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
					     /*  use default output scalings  */
	     /*  this is a no-op as long as output data are stored as floats  */
    if (docsinv) {
      cs_res->bscale = seg_cs->bscale;
      cs_res->bzero = seg_cs->bzero;
    }
    if (dovinv) {
      vx_res->bscale = seg_vx->bscale;
      vx_res->bzero = seg_vx->bzero;
      vy_res->bscale = seg_vy->bscale;
      vy_res->bzero = seg_vy->bzero;
      vz_res->bscale = seg_vz->bscale;
      vz_res->bzero = seg_vz->bzero;
    }
		  /*  need to transpose input array to individual components  */
    for (k=0; k<11; k++) {
      pln = k * cols * rows;
      pln0 = 4 * k * cols * rows;
      pln1 = (4 * k + 1) * cols * rows;
      pln2 = (4 * k + 2) * cols * rows;
      pln3 = (4 * k + 3) * cols * rows;
      for (row = 0; row < rows; row++){
	pln += cols;
	pln0 += cols;
	pln1 += cols;
	pln2 += cols;
	pln3 += cols;
	for (col = 0; col < cols; col++){
	  oi_mn[pln + col] = input[pln0 + col]; 
	  oi[pln + col] = input[pln1 + col];
	  we[pln + col] = input[pln2 + col];
	  ns[pln + col] = input[pln3 + col];
	}
      }
    }

    if (docsinv) inver_cs_sub_ (oi_mn, kcs, cs);
    if (dovinv) {
      float latc = drms_getkey_float (ttrec, "LatHG", &status);
      float lonc = drms_getkey_float (ttrec, "LonCM", &status);
      inver_v_sub_ (oi, we, ns, &lonc, &latc, kvx, kvy, kvz, vx, vy, vz,
	  &kerntyp);
    }

    status = 0;
    if (docsinv) status += drms_segment_write (seg_cs, cs_res, 0);
    if (dovinv) {
      status += drms_segment_write (seg_vx, vx_res, 0);
      status += drms_segment_write (seg_vy, vy_res, 0);
      status += drms_segment_write (seg_vz, vz_res, 0);
    }
    if (status) {
      drms_sprint_rec_query (recid, orec);
      fprintf (stderr, "Error writing data to record %s\n", recid);
      fprintf (stderr, "      series may not have appropriate structure\n");
      fprintf (stderr, "      record %d skipped\n", rgn);
      drms_close_record (orec, DRMS_FREE_RECORD);
      continue;
    }
// fprintf (stderr, "inserting output record\n");
    drms_close_record (orec, DRMS_INSERT_RECORD);
// fprintf (stderr, "output record inserted\n");
  }
  return 0;
}
/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  v 0.2 frozen 11.03.11 as working reference
 *  11.03.12	added verbose flag, some reporting and documentation; travel-
 *		times input from DRMS record set; removed supefluous read fits
 *		function
 *  11.03.15	output to DRMS
 *  11.03.18	kernels read from DRMS by driver 
 *  v 0.5 frozen 11.03.23 as reference version
 *  11.03.23	processing of multiple records, some code simplification;
 *		setting of a few keywords; removed include of fftw.f from
 *		mcd_cs_hmi_sub.f in favor of declaration of sole needed
 *		parameter
 *  v 0.6 frozen 11.03.25
 *  v 0.7	Check version of pixel locator used in travel-time map and
 *		propagate by default
 *		Add arguments for region location to call to inver_v_sub_
 *		Add argument for kernel type to call to "
 *		Add module argument for fitting method to be used, and setting
 *		of keyword Fitting
 *		Add setting (and checking) of keyword KrnlAppr
 *		Check for existence of required segments in output record, and
 *		modify attempted inversions accordingly; check for agreement of
 *		ouput segment structure with hard-coded dimensions
 *  v 0.8 frozen 2012.08.04
 *  12.09.28	added option of adding to or overriding standard list of
 *		keywords from input to output
 *		removed defaults for input and output sets, changed default
 *		fitting method to both
 *  13.04.08	added WCS keys to propagation list
 *  v 0.81 frozen 2013.04.08
 */
