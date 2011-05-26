/**
   @defgroup hg_patch hg_patch extracts a rectangular 
   region or patch from an input series.  The patch is tracked across
   the solar disk from time  t_start to t_stop.  The patch is defined by either its
   Carrington rotation, longitude, latitude and size when it crosses central meridian or
   by its location in arcsec at a reference time.  If defined by Carrington coordinates the
   bounding box size is specified in degrees of longitude and latitude at disk center.  If the box is
   specified by its reference image time and location in arcsec from disk center,
   its size is also specified in arcsec.

   @ingroup a_programs

   @brief hg_patch extracts a rectangular 
   heliographic region or patch from an input series.  The patch is tracked across
   the solar disk from t_start to t_stop times.  The patch may be defined by Carrington coordinates or
   in arcsec.  The extracted patch is rectangular in original image pixels.

      Written by: Bala Poduval as extract_region, modified by Phil Scherrer

   @par Synopsis:

   @code
   hg_patch in=<recordset> out=<series_out>  \
   log=<logfile> \
   t_start=<start time> \
   t_stop=<end time> \
   cadence=<step time> \
   car_rot=<Carrington Rotation> \
   t_ref=<ref_time> \
   locunits=<location units> \
   x=<ew location> \
   y=<ns location> \
   boxunits= <patch spec units> \
   height=<height> \
   width=<width> \
   requestid=<RequestID> \
   [-t] \

   @endcode

   Cadence is the desired time between extracted frames.  It will be modified to be a multiple of the dataset's step size.
   The extracted region location may be specified by position at a reference time or by Carrington coordinates.
   The position may be specified in arcsec or Stonyhurst coordinates.
   If the car_rot parameter is specified, the Carrington location method is used.
   If the t_ref parameter is present, location is specified by a method given in <locunits>.
   The location is given in the "x" and "y" arguments where the units of x and y are determined by locunits.
   The choices for locunits are: "arcsec", "pixels", "stony", "carrlong".
   For "stony" (Stonyhurst) x is the CM distance in degrees, for "carrlong" x is Carrington longitude.
   The y value is in degrees for either "stony" or "carrlong".
   locunits can be "arcsec" or "stony". 

   The size of the patch, given in the width and height parameters can be expressed in "pixels", "arcsec", or "degrees"
   as determined by the "boxunits" parameter.

   If the logfile is present, a RecordSet query will be written to it for each image created  containing a query
   that will returnthat image.

   RequestID, if present, is the ID of an on-requst processing export request.

   If t_start or t_stop are not present, those values will be inferred from
   the on-disk time span for the center of the patch -90 degrees to +90 degrees from CM.
   If cadence is not specified, the full cadence of the dataset will be used.
   This module computes the beginning and ending time for tracking the 
   region identified at time = t_ref, x and y or by the Carrington coordinates of the box center.
   In the Carrington case the module extracts a rectangular region 
   with user defined size. the height and width of the box are specified in degrees
   of latitude and longitude with the pixel size of the box is computed from the
   projection of the box at CM.  The pixel size of the rectangular box remains constant 
   during tracking.  The box height defaults to the width.  The width defaults to 10 degrees.
  
   If the input recordset spec is exactly "[$]" it will be discarded and the time limits will be taken from t_start, t_stop, car_rot, and/or t_ref as appropriate.  This is to allow exports via jsoc_fetch to not require a bounding recordset to be specified.

   If t_start and/or t_stop are specified the referecne time or disk center location might be not
   included.
   
   If the input is only a seriesname, it must have a prime key of type time.

   If the t_ref parameter is specified, it must refer to a non-missing image in the input series.

   The output seriesname defaults to the input seriesname with a suffix of "_hgpatch".

   The -t, 'no tracking' flag casues the extracted region to remain fixed with respect to disc center.  I.e. the Carrington
   tracking is disabled.

   @par Flags:
   -t   Disable tracking.

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in The input data series.
   @param out The output series.

   @par Exit_Status:

   @par Example:
   Extracts a rectangular region of size 20 x 10 degrees and tracks it from 
   when the region is just appearing at the E-limb to t_stop.

   @code
   hg_patch in='mdi.fd_M_96m_lev18' out='su_bala.extractreg' locunits=carrlong boxunits=pixels car_rot=2009 x=295 y=5 width=30 height=20

   @endcode

   @par Example:

   @code

   @endcode

   @bug
   No doubt some bugs.


*/
#include "jsoc.h"
#include "jsoc_main.h"
#include "astro.h"
#include "fstats.h"
#include "atoinc.h"
#include <math.h>

// void HeliographicLocation(TIME t, int *crot, double *L, double *B);
// TIME HeliographicTime(int crot, double L);

char *module_name = "hg_patch";

#define NOTSPECIFIED "***NOTSPECIFIED***"
#define	DIE(msg) {fprintf(stderr,"%s  Status=%d\n",msg, status); return(status?status:1);}
#define	DIE2(msg,val) {fprintf(stderr,"%s %s,  Status=%d\n",msg, val,  status); return(status?status:1);}
#define TEST_PARAM(param) {if (status) DIE2("Required keyword missing: ", param);}

#define     Deg2Rad    (M_PI/180.0)
#define	    Rad2arcsec (3600.0/Deg2Rad)
#define	    Rad2Deg    (180.0/M_PI)

ModuleArgs_t module_args[] =
  {
    {ARG_STRING, "in", "NOTSPECIFIED", "input series. e.g 'mdi.fd_M_96m_lev18'"},
    {ARG_STRING, "out", "NOTSPECIFIED", "output seies. e.g. 'su_bala.extractreg'"},
    {ARG_STRING, "log", "NOTSPECIFIED", "output log file of records made"},
    {ARG_STRING, "requestid", "NOTSPECIFIED", "RequestID if hg_patch call originated in an export request."},
    {ARG_INT, "car_rot", "-1", "Carrington Rotation when the region crosses CM"},
    {ARG_FLOAT, "width", "0", "width of box in degrees of longitude"},
    {ARG_FLOAT, "height", "0", "height of box in degrees of latitude when it corsses CM"},
    {ARG_STRING, "boxunits", "pixels", "units of patch, 'pixels', 'arcsecs', or 'degrees'"},
    {ARG_TIME, "t_start", "JD_0", "Start time, defaults to time at 90E"},
    {ARG_TIME, "t_stop", "JD_0", "End time, defauolts to 90W"},
    {ARG_TIME, "cadence", "NOTSPECIFIED", "Cadence of product, defaults to input cadence."},
    {ARG_TIME, "t_ref", "JD_0", "Time for which x and y apply, implies ref image."},
    {ARG_STRING, "locunits", "pixels", "Location units in 'pixels', 'arcsec', or 'degrees'"},
    {ARG_STRING, "where", "1=1", "Additional 'where' clause if needed"},
    {ARG_FLOAT, "x", "0", "Location of extract box center"},
    {ARG_FLOAT, "y", "0", "Location of extract box center"},
    {ARG_FLAG, "t", "0", "Disable Carrington rate tracking"},
    {ARG_END}
  };

int img2sphere (double x, double y, double ang_r, double latc, double lonc,
        double pa, double *rho, double *lat, double *lon, double *sinlat,
        double *coslat, double *sig, double *mu, double *chi);
int sphere2img (double lat, double lon, double latc, double lonc,
        double *x, double *y, double xcenter, double ycenter,
        double rsun, double peff, double ecc, double chi,
        int xinvrt, int yinvrt);

// experimental make recordset list at cadence
char *get_input_recset(DRMS_Env_t *env, char *in, TIME cadence);

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
// These are in units where the first pixel is 1 not 0.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

#define BOXBAD 0
#define BOXARCSEC 1
#define BOXPIXELS 2
#define BOXDEGREE 3

#define LOCBAD 0
#define LOCARCSEC 1
#define LOCPIXELS 2
#define LOCSTONY 3
#define LOCCARR 4

#define TIMES_RECSET 0   // recordset spec contains time range
#define TIMES_GIVEN 1    // explicit t_start and t_stop provided
#define TIMES_IMPLICIT 2 // times to be deduced for disk transit of specified box
int DoIt(void)
{
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *inRS, *outRS;
  DRMS_Record_t *inRec, *outRec, *inTemplate, *outTemplate;
  DRMS_Segment_t *inSeg, *outSeg;
  DRMS_Array_t *inArray, *outArray;
  int i, ii, status = DRMS_SUCCESS, nrecs; 
  int  irec;
  double center_x, center_y, crpix1, crpix2, x0, y0; 
  double rsun_ref, dsun_obs, rsun, rsun_rad, rsunpix;
  double crln_obs, crlt_obs;
  double crln_obs_rad, crlt_obs_rad;
  double crln_rad, crlt_rad;
  double pa_rad;
  double cdelt;
  double urx, ury, llx, lly;
  int inAxis[2];
  int this_car_rot;
  char *ctype1, *ctype2;
  TIME t_rec;
  double box_x, box_y;
  double crln, crlt;
  char outseries[DRMS_MAXNAMELEN];
  char inseries[DRMS_MAXNAMELEN];
  char inQuery[DRMS_MAXQUERYLEN];
  char in[DRMS_MAXQUERYLEN];

  const char *ingiven = params_get_str(params, "in");
  char *inparam;
  char *lbracket;
  char *moreQuery = NULL;
  const char *outparam = params_get_str(params, "out");
  const char *logfile = params_get_str(params, "log");
  const char *requestid = params_get_str(params, "requestid");
  const char *where = params_get_str(params, "where");
  TIME t_start = params_get_time(params, "t_start");
  TIME t_stop = params_get_time(params, "t_stop");
  const char *cadence_str = params_get_str(params, "cadence");
  double cadence;
  TIME t_ref = params_get_time(params, "t_ref");
  double width = params_get_double(params, "width");
  double height = params_get_double(params, "height");
  int car_rot = params_get_int(params, "car_rot");
  double x = params_get_double(params, "x");
  double y = params_get_double(params, "y");
  const char *boxunits = params_get_str(params, "boxunits");
  const char *locunits = params_get_str(params, "locunits");
  int NoTrack = params_isflagset(params, "t");
  TIME tNotSpecified = sscan_time("JD_0");
  int do_reftime = 0;  // Target specified by reftime vs Carr rotation mode.
  double crvalx = 0.0;
  double crvaly = 0.0;
  double crota, sina, cosa;
  double width_arcsec, height_arcsec;
  double pa, deltlong;
  double center_x_first, center_y_first;  // used for NoTrack option
  int firstimage = 1;
  int boxtype, loctype;
  char *timekeyname;
  TIME t_step, t_epoch;
  int npkeys;
  int times_source = TIMES_RECSET;
  char timebuf[20];
  char *in_filename = NULL;

  FILE *log = NULL;

  if (strncasecmp(locunits, "arcsec", 6) == 0) loctype = LOCARCSEC;
  else if (strncasecmp(locunits, "pixels", 3) == 0) loctype = LOCPIXELS;
  else if (strncasecmp(locunits, "stony", 5) == 0) loctype = LOCSTONY;
  else if (strncasecmp(locunits, "carrlong", 4) == 0) loctype = LOCCARR;
  else loctype = LOCBAD;

  if (strncasecmp(boxunits, "arcsec", 6) == 0) boxtype = BOXARCSEC;
  else if (strncasecmp(boxunits, "pixels", 3) == 0) boxtype = BOXPIXELS;
  else if (strncasecmp(boxunits, "degrees", 3) == 0) boxtype = BOXDEGREE;
  else loctype = BOXBAD;

fprintf(stderr,"starting boxtype=%d, loctype=%d\n",boxtype,loctype);
  if (loctype == LOCCARR)
    {
    if (car_rot < 0) DIE("Carrington rotation number must be provided for locunits=carrlong");
    if (NoTrack) DIE("Can not use locunits=carrlong for no-tracking mode");
    do_reftime = 0;
    }
  else 
    {
    if (t_ref < 0) DIE("t_ref must be provided for locunits!=carrlong");
    do_reftime = 1;
    }

  if (loctype==LOCBAD)
    DIE2("Location units not detected, locunits", locunits)
  if (boxtype==BOXBAD)
    DIE2("patch dimensions not understood,", boxunits);

  if (strcmp(logfile, "NOTSPECIFIED") != 0)
    {
    log = fopen(logfile, "w");
    if (!log) DIE2("Can not create log file.",logfile);
    }

  inparam = strdup(ingiven);
  if (strcmp(inparam, "NOTSPECIFIED") == 0) DIE("Input series must be specified.");
  lbracket = index(inparam, '[');
  // first, get input series names.
  if (lbracket)
    {
    int n = lbracket - inparam;
    strncpy(inseries, inparam, n);
    inseries[n] = '\0';
    if (strncmp(lbracket, "[$]", 3) == 0) // Special case, discard the explicit last-record spec to allow exports via jsoc_fetch.
      {
      moreQuery = lbracket + 3;
      *lbracket = '\0';
      lbracket = NULL;
      if (t_start == tNotSpecified || t_stop == tNotSpecified)
          times_source = TIMES_IMPLICIT;
      else
          times_source = TIMES_GIVEN;
      }
    else
      moreQuery = lbracket;
    }
  else
    {
    strcpy(inseries, inparam);
    if (t_start == tNotSpecified || t_stop == tNotSpecified)
        times_source = TIMES_IMPLICIT;
    else
        times_source = TIMES_GIVEN;
    }
    
// XXXXXXXXXXX Get box location in Carrington Coords for all location types XXXXXXXXXXXXXX
  if (do_reftime) // Arc-sec, pixel, or Stonyhurst specification, get ref image information
    { // the image for ref_time is supposed to be present in the series.
    char t_ref_text[100];
    sprint_at(t_ref_text, t_ref-7200);
    sprintf(in, "%s[%s/4h][? QUALITY >=0 ?]", inseries, t_ref_text);
    inRS = drms_open_records(drms_env, in, &status); if (status || inRS->n == 0) DIE2("No input data found for t_ref",in);
    int irec, nrecs = inRS->n;
    TIME tdiff = 10000;
    for (irec=0; irec<nrecs; irec++) // find record close to t_ref
      {
      TIME newdiff;
      inRec = inRS->records[irec];
      if (drms_getkey_int(inRec,"QUALITY",NULL) < 0)
        continue;
      if ((newdiff=(fabs(drms_getkey_time(inRec,"T_OBS",NULL) - t_ref))) > tdiff)
        break;
      tdiff = newdiff;
      }
    inRec = inRS->records[irec];
    this_car_rot = drms_getkey_int(inRec, "CAR_ROT", &status); TEST_PARAM("CAR_ROT");
    crlt_obs = drms_getkey_double(inRec, "CRLT_OBS", &status); TEST_PARAM("CRLT_OBS");
    crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status); TEST_PARAM("CRLN_OBS");
    crpix1 = drms_getkey_double(inRec, "CRPIX1", &status); TEST_PARAM("CRPIX1");
    x0 = crpix1 - 1;
    crpix2 = drms_getkey_double(inRec, "CRPIX2", &status); TEST_PARAM("CRPIX2");
    y0 = crpix2 - 1;
    crvalx = drms_getkey_double(inRec, "CRVAL1", &status); TEST_PARAM("CRVAL1");
    crvaly = drms_getkey_double(inRec, "CRVAL2", &status); TEST_PARAM("CRVAL2");
    crota = drms_getkey_double(inRec, "CROTA2", &status); TEST_PARAM("CROTA2");
    pa = -crota;
    crlt_obs_rad = Deg2Rad * crlt_obs;
    crln_obs_rad = Deg2Rad * crln_obs;
    pa_rad = Deg2Rad * pa;
    sina = sin(crota*Deg2Rad);
    cosa = cos(crota*Deg2Rad);
    ctype1 = strdup(drms_getkey_string(inRec, "CTYPE1", &status)); TEST_PARAM("CTYPE1");
      if (strcmp(ctype1, "HPLN-TAN") != 0) DIE2("CTYPE1 not HPLN-TAN as required, is: ", ctype1);
    ctype2 = strdup(drms_getkey_string(inRec, "CTYPE2", &status)); TEST_PARAM("CTYPE2");
      if (strcmp(ctype2, "HPLT-TAN") != 0) DIE2("CTYPE2 not HPLT-TAN as required, is: ", ctype2);
    rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
      if (status) rsun_ref = 6.96e8;
    dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status); TEST_PARAM("DSUN_OBS");
    rsun_rad = asin(rsun_ref/dsun_obs); 
    rsun = rsun_rad*Rad2arcsec; 
    cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
    rsunpix = rsun/cdelt;
    if (loctype == LOCPIXELS || loctype==LOCARCSEC)
      { /* box ref location is in arcsec - find Carrington equivalent */
      if (loctype==LOCARCSEC)
        {
        center_x = PIX_X(x,y) - 1 - x0;
        center_y = PIX_Y(x,y) - 1 - y0;
        }
      else
        {
        center_x = x - 1 - x0;
        center_y = y - 1 - y0;
        }
      inSeg = drms_segment_lookupnum(inRec, 0);
      inAxis[0] = inSeg->axis[0];
      inAxis[1] = inSeg->axis[1];
      if ( img2sphere(center_x/rsunpix, center_y/rsunpix, rsun_rad, crlt_obs_rad, crln_obs_rad, pa_rad,
           NULL, &crlt_rad, &crln_rad, NULL, NULL, NULL, NULL, NULL) < 0)
          DIE("Starting location is off the solar disk.");
      crln = crln_rad * Rad2Deg;
      crlt = crlt_rad * Rad2Deg;
      }
    else /* loctype==LOCSTONY */
      {
      crlt = y;
      crln = crln_obs + x;
      }
    deltlong = crln - crln_obs;
    car_rot = this_car_rot;
    if (deltlong >= 360.0)
      { car_rot--; crln -= 360.0; }
    else if (deltlong < 0)
      { car_rot++; crln += 360.0; }
sprint_at(timebuf,drms_getkey_time(inRec,"T_OBS",NULL));
fprintf(stderr,"t_ref specified, box center is at %4d:%05.1f by %05.1f on %s\n",car_rot, crln,crlt, timebuf);
    drms_close_records(inRS, DRMS_FREE_RECORD);
    }
  else // Carrington specification
    {
    crln = x;
    crlt = y;
    if (crln < 0) DIE("Box longitude must be specified.");
    if (crlt < -990) DIE("box latitude must be specified.");
fprintf(stderr,"car_rot specified, box center is at %4d:%05.1f by %05.1f \n",car_rot, crln,crlt);
    }
  crln_rad = crln * Deg2Rad;
  crlt_rad = crlt * Deg2Rad;

  // Now we have crln, crlt, crln_rad, crlt_rad, car_rot for the target box. 

  // Get implied time limits from car_rot and crln
  t_ref = HeliographicTime(car_rot, crln);
  if (times_source == TIMES_IMPLICIT)
      {
      if (t_start == tNotSpecified)
        t_start = HeliographicTime(car_rot, crln + 90);
      if (t_stop == tNotSpecified)
        t_stop = HeliographicTime(car_rot, crln - 90);
      }
// XXXXXXXXXXXXXXXXX End of get target box location

// XXXXXXXXXXXXXXXX get input seriesname and output seriesname
// first, get input and output series names.

  inTemplate = drms_template_record(drms_env, inseries, &status);
  if (status || !inTemplate) DIE2("Input series can not be found: ", inseries);

  if (strcmp(outparam, "NOTSPECIFIED") == 0)
    {
    strncpy(outseries, inseries, DRMS_MAXNAMELEN);
    strncat(outseries, "_hgpatch", DRMS_MAXNAMELEN);
    }
  else
   strncpy(outseries, outparam, DRMS_MAXNAMELEN);

  // Now, make sure output series exists and get template record.
  outTemplate = drms_template_record(drms_env, outseries, &status);
  if (status || !outTemplate) DIE2("Output series can not be found: ", outseries);

  // Now find the prime time keyword name
  npkeys = inTemplate->seriesinfo->pidx_num;
  timekeyname = NULL;
  if (npkeys > 0)
    {
    int i;
    for (i=0; i<npkeys; i++)
        {
        DRMS_Keyword_t *pkey, *skey;
        pkey = inTemplate->seriesinfo->pidx_keywords[i];
        if (pkey->info->recscope > 1)
           pkey = drms_keyword_slotfromindex(pkey);
        if (pkey->info->type != DRMS_TYPE_TIME)
          continue;
	if(i > 0) DIE("Input series must have TIME keyword first, for now...");
        timekeyname = pkey->info->name;
        t_step = drms_keyword_getdouble(drms_keyword_stepfromslot(pkey), &status);
	if (status) DIE("problem getting t_step");
        t_epoch = drms_keyword_getdouble(drms_keyword_epochfromslot(pkey), &status);
	if (status) DIE("problem getting t_epoch");
        }
    }
  else
    DIE("Must have time prime key");

  // Get cadence for output data.  The output series should be slotted on a multiple of
  // the input series, multiple may be 1.
  
  if (strcmp(cadence_str, "NOTSPECIFIED") != 0)
    {
    int ratio;
    cadence = atoinc((char *)cadence_str);
    ratio = round(cadence/t_step);
    cadence = ratio * t_step;
    }
  else
    cadence = t_step;

  // Finally, get input recordset spec
  if (lbracket && times_source == TIMES_RECSET) // RecordSet query specified
    strncpy(in, inparam, DRMS_MAXQUERYLEN);
  else // need to generate query from limit times
    { 
    char t_start_text[100], t_stop_text[100], cadence_text[100];
    if (cadence > t_step)
      {
      sprintf(cadence_text,"@%fs", cadence);
      if (times_source == TIMES_IMPLICIT) // round start and stop to cadence slots for auto limits
        {
        int nsteps;
        nsteps = round((t_start - t_epoch)/cadence);
        t_start = t_epoch + nsteps * cadence;
        nsteps = round((t_stop - t_epoch)/cadence);
        t_stop = t_epoch + nsteps * cadence;
        }
      }
    else
      cadence_text[0] = '\0';
    // strncpy(inseries, inparam, DRMS_MAXNAMELEN);
    sprint_at(t_start_text, t_start);
    sprint_at(t_stop_text, t_stop);
    if (strcmp(inseries,"aia.lev1")==0 && t_step == 1.0 && cadence > 1.0) // special case for AIA slots
        {
        // experimental, AIA is not tseq slots so get vector of times and convert to list of records
        // rounded to nearest slots.  Must put this list in a temp file since may be big.
        sprintf(in, "%s[%s-%s]%s[? %s ?]", inseries, t_start_text, t_stop_text, (moreQuery ? moreQuery : ""), where);
        in_filename = get_input_recset(drms_env, in, cadence);
        if (!in_filename) DIE("Cant make AIA cadence recordset list file");
        sprintf(in, "@%s", in_filename);
        }
    else // normal case
        sprintf(in, "%s[%s-%s%s]%s[? %s ?]", inseries, t_start_text, t_stop_text, cadence_text, (moreQuery ? moreQuery : ""), where);
    }

// XXXXXXXXXXXXX End of get input and output series information

// XXXXXXXXXXXXX Now read data and extract boxes.

  inRS = drms_open_records(drms_env, in, &status);
  if (status || inRS->n == 0)
           DIE("No input data found");
  nrecs = inRS->n;

  // extract patches from each record
  for (irec = 0; irec < nrecs; irec ++)
    {
    inRec = inRS->records[irec];
    if (status || !inRec) DIE("Record read failed.");
    TIME trec = drms_getkey_time(inRec, timekeyname, &status); TEST_PARAM(timekeyname);
    if (trec < t_start)
      continue;
    if (trec > t_stop)
      break;
    // skip records without images.
    if (drms_keyword_lookup(inRec, "QUALITY", 1))
      {
      int quality = drms_getkey_int(inRec, "QUALITY", &status);
      if (quality < 0)
         continue;
      }
    
    // get coordinate information for this image
    this_car_rot = drms_getkey_int(inRec, "CAR_ROT", &status); TEST_PARAM("CAR_ROT");
    crlt_obs = drms_getkey_double(inRec, "CRLT_OBS", &status); TEST_PARAM("CRLT_OBS");
    crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status); TEST_PARAM("CRLN_OBS");
    crpix1 = drms_getkey_double(inRec, "CRPIX1", &status); TEST_PARAM("CRPIX1");
    x0 = crpix1 - 1;
    crpix2 = drms_getkey_double(inRec, "CRPIX2", &status); TEST_PARAM("CRPIX2");
    y0 = crpix2 - 1;
    pa = -drms_getkey_double(inRec, "CROTA2", &status); TEST_PARAM("CROTA2");
    // convert to radians
    crlt_obs_rad = Deg2Rad*crlt_obs;;
    crln_obs_rad = Deg2Rad*crln_obs;;
    pa_rad = Deg2Rad*pa;

//  XXXXXXXXXXXXXXXXX get coordinate mapping info using keywords from first record
    if (firstimage)
      {
      firstimage = 0;
      if (loctype==LOCCARR)
        {
        ctype1 = strdup(drms_getkey_string(inRec, "CTYPE1", &status)); TEST_PARAM("CTYPE1");
        if (strcmp(ctype1, "HPLN-TAN") != 0) DIE2("CTYPE1 not HPLN-TAN as required, is: ", ctype1);
        ctype2 = strdup(drms_getkey_string(inRec, "CTYPE2", &status)); TEST_PARAM("CTYPE2");
        if (strcmp(ctype2, "HPLT-TAN") != 0) DIE2("CTYPE2 not HPLT-TAN as required, is: ", ctype2);
        rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
        if (status) rsun_ref = 6.96e8;
        dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status); TEST_PARAM("DSUN_OBS");
        rsun = asin(rsun_ref/dsun_obs)*Rad2arcsec; 
        cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
        // in principle use deltlong to get time to use for correcting crlt_obs, ignore for now
        // get ur, ll coords of box at CM.
        rsunpix = rsun/cdelt;
        }
      if (boxtype == BOXDEGREE)
        {
        int nx, ny;
        // get corners of sample box at zero longitude wrt center.
        sphere2img(Deg2Rad*(crlt+height/2), Deg2Rad*(width/2), crlt_obs_rad, 0.0,
          &urx, &ury, 0.0, 0.0, rsunpix, pa_rad, 0, 0, 0, 0);
        sphere2img(Deg2Rad*(crlt-height/2), -Deg2Rad*(width/2), crlt_obs_rad, 0.0,
          &llx, &lly, 0.0, 0.0, rsunpix, pa_rad, 0, 0, 0, 0);
        nx = urx-llx;
        ny = ury-lly;
        // fix for 180 flipped coords.
        if (nx < 0) nx = -nx;
        if (ny < 0) ny = -ny;
        llx = -nx/2;  lly = -ny/2;  urx = nx/2; ury = ny/2;
        urx = round(urx); ury = round(ury); llx = round(llx); lly = round(lly);
        // make size odd, increase smaller distance from center
        nx = urx - llx + 1;  ny = ury - lly + 1;
        if (!(nx&1)) {if (urx > -llx) llx-=1; else urx+=1;}
        if (!(ny&1)) {if (ury > -lly) lly-=1; else ury+=1;}
        width_arcsec = nx * cdelt;
        height_arcsec = nx * cdelt;
        }
      else if (boxtype==BOXARCSEC)
        {
        int nx, ny;
        urx = width/(2.0*cdelt); llx = -urx; ury = height/(2.0*cdelt); lly = -ury;
        urx = round(urx); ury = round(ury); llx = round(llx); lly = round(lly);
        nx = urx - llx + 1;  ny = ury - lly + 1;
        width_arcsec = nx * cdelt;
        height_arcsec = nx * cdelt;
        }
      else /* boxtype == BOXPIXELS */
        {
        int nx, ny;
        int iurx, iury, illx, illy;
        iurx = round(width/2.0); illx = -iurx; iury = round(height/2.0); illy = -iury;
        if (((int)(round(width)) & 1) == 0 && ((iurx-illx+1) & 1) == 1) iurx--;
        if (((int)(round(height)) & 1) == 0 && ((iury-illy+1) & 1) == 1) iury--;
        urx = iurx; ury = iury; llx = illx; lly = illy;
        nx = urx - llx + 1;  ny = ury - lly + 1;
        width_arcsec = nx * cdelt;
        height_arcsec = nx * cdelt;
        }
      if (NoTrack)
          {
          sphere2img(crlt_rad, crln_rad, crlt_obs_rad, crln_obs_rad, &center_x, &center_y, x0, y0, rsunpix, pa_rad, 0, 0, 0, 0);
          center_x_first = center_x - x0;
          center_y_first = center_y - y0;
fprintf(stderr,"NoTrack, center_x_first=%f, center_y_first=%f, pa=%f\n",center_x_first,center_y_first,pa);
          }
      }

    if (NoTrack)
        {
        center_x = center_x_first + x0;
        center_y = center_y_first + y0;
        }
    else
        {
        sphere2img(crlt_rad, crln_rad, crlt_obs_rad, crln_obs_rad, &center_x, &center_y, x0, y0, rsunpix, pa_rad, 0, 0, 0, 0);
        }
    int x1 = round(center_x + llx);
    int y1 = round(center_y + lly);
    int x2 = round(center_x + urx);
    int y2 = round(center_y + ury);
    crpix1 = 1 + x0 - x1;
    crpix2 = 1 + y0 - y1;

    int start1[2] = {x1, y1};
    int end1[2] = {x2, y2};

    inSeg = drms_segment_lookupnum(inRec, 0);
    inAxis[0] = inSeg->axis[0];
    inAxis[1] = inSeg->axis[1];

    if (x1 >= inAxis[0] || y1 >= inAxis[1] || x2 < 0 || y2 < 0)
      continue; // slice outside image

    if (x1>=0 && y1>=0 && x2<inAxis[0] && y2<inAxis[1])
      {
      outArray = drms_segment_readslice(inSeg, DRMS_TYPE_FLOAT, start1, end1, &status);
      if (status) DIE("Cant read input record");
      }
    else
      {
      int dims[2] = {x2-x1+1,y2-y1+1};
      int start2[2], end2[2];
      int x,y;
      outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, NULL, &status);
      drms_array2missing(outArray);
      start2[0] = x1 < 0 ? 0 : x1;
      start2[1] = y1 < 0 ? 0 : y1;
      end2[0] = x2 >= inAxis[0] ? inAxis[0]-1 : x2;
      end2[1] = y2 >= inAxis[1] ? inAxis[1]-1 : y2;
      inArray = drms_segment_readslice(inSeg, DRMS_TYPE_FLOAT, start2, end2, &status);
      if (status) DIE("Cant read input record");
      int start3[2] = {start2[0]-start1[0],start2[1]-start1[1]}; // fetched slice offset in outarray
      int end3[2] = {end2[0]-end1[0],end2[1]-end1[1]};           // fetched slice offset in outarray
      int n3x = end2[0] - start2[0] + 1;
      int n3y = end2[1] - start2[1] + 1;
      for (y=0; y< n3y; y++)
        for (x=0; x < n3x; x++)
          *((float *)outArray->data + ((start3[1]+y)*dims[0] + x + start3[0])) =
            *((float *)inArray->data + (y*n3x + x));
      drms_free_array(inArray);
      }
    // Handle special case where |pa| == 180
    if (fabs(fabs(pa)-180.0) < 1.0)
      {
      // Tracking should be OK for center, but patch rotated.  Flip on both axes.
      float *data = (float *)outArray->data;
      int i,j;
      int nx = outArray->axis[0];
      int ny = outArray->axis[1];
      int midrow = (ny)/2;
      int midcol = (nx)/2;
      float val;
      for (j=0; j<midrow; j++)
        for (i=0; i<nx; i++)
          {
          val = data[j*nx + i];
          data[j*nx + i] = data[(ny - 1 - j)*nx + nx - 1 - i];
	  data[(ny - 1 - j)*nx + nx - 1 - i] = val;
	  }
      if ((nx & 1) != 0)
        {
        for (i=0; i<midcol; i++)
          {
          val = data[midrow*nx + i];
          data[midrow*nx + i] = data[midrow*nx + nx - 1 - i];
	  data[midrow*nx + nx - 1 - i] = val;
          }
        }
      pa -= 180.0;
      crpix1 = 1 + x2 - x0;
      crpix2 = 1 + y2 - y0;
      cosa = 1.0; sina = 0.0;
      }

/*
 *               writing the extracted region data file
*/
    outRS = drms_create_records(drms_env, 1, outseries, DRMS_PERMANENT, &status);
    if (status) DIE("Cant make outout record");
    outRec = outRS->records[0];
    drms_copykeys(outRec, inRec, 1, kDRMS_KeyClass_Explicit);
    outSeg = drms_segment_lookupnum(outRec, 0);
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    set_statistics(outSeg, outArray, 1);
    drms_setkey_float(outRec, "XCEN", WX((outArray->axis[0]+1)/2, (outArray->axis[1]+1)/2));
    drms_setkey_float(outRec, "YCEN", WY((outArray->axis[0]+1)/2, (outArray->axis[1]+1)/2));
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) DIE("problem writing file");
    drms_free_array(outArray);
    drms_sprint_rec_query(inQuery, inRec);
    drms_setkey_string(outRec, "SOURCE", inQuery);
    drms_setkey_time(outRec, "T_REC", trec);
    drms_setkey_double(outRec, "CRPIX1", crpix1);
    drms_setkey_double(outRec, "CRPIX2", crpix2);
    drms_setkey_double(outRec, "CRVAL1", 0.0);
    drms_setkey_double(outRec, "CRVAL2", 0.0);
    drms_setkey_double(outRec, "CRDELT1", cdelt);
    drms_setkey_double(outRec, "CRDELT2", cdelt);
    drms_setkey_double(outRec, "CROTA2", 0.0);
    drms_setkey_string(outRec, "CONTENT", "Tracked Extracted Patches, made by hg_patch");
    drms_setkey_string(outRec, "COMMENTS", "Patches");
    drms_setkey_string(outRec, "RequestID", requestid);
    drms_setkey_string(outRec, "HGBOXUNITS", boxunits);
    drms_setkey_string(outRec, "HGLOCUNITS", locunits);
    drms_setkey_float(outRec, "HGWIDE", width);
    drms_setkey_float(outRec, "HGHIGH", height);
    drms_setkey_float(outRec, "HGASWIDE", width_arcsec);
    drms_setkey_float(outRec, "HGASHIGH", height_arcsec);
    drms_setkey_float(outRec, "HGCRLN", crln);
    drms_setkey_float(outRec, "HGCRLT", crlt);
    drms_setkey_int(outRec, "HGCARROT", car_rot);
    drms_setkey_time(outRec, "HGTSTART", t_start);
    drms_setkey_time(outRec, "HGTSTOP", t_stop);
    drms_setkey_time(outRec, "DATE", time(0) + UNIX_EPOCH);
    drms_setkey_string(outRec, "HGQUERY", in);
    // drms_copykey(outRec, inRec, "EXPTIME");
    // drms_copykey(outRec, inRec, "QUALITY");
    // drms_copykey(outRec, inRec, "DSUN_OBS");
    // drms_copykey(outRec, inRec, "CRLN_OBS");
    // drms_copykey(outRec, inRec, "CRLT_OBS");
    // drms_copykey(outRec, inRec, "OBS_VR");
    // drms_copykey(outRec, inRec, "OBS_VW");
    // drms_copykey(outRec, inRec, "OBS_VN");
    // drms_copykey(outRec, inRec, "CAR_ROT");
    // drms_copykey(outRec, inRec, "T_OBS");
    // drms_copykey(outRec, inRec, "DATE__OBS");
    // drms_copykey(outRec, inRec, "WAVELNTH");
fprintf(stderr,"Box %04d ",irec);
drms_fprint_rec_query(stderr, outRec);
fprintf(stderr,"  crpix1=%f, crpix2=%f\n",crpix1,crpix2);

    if (log)
      {
      drms_fprint_rec_query(log, outRec);
      fprintf(log, "\n");
      }
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    } 

    drms_close_records(inRS, DRMS_FREE_RECORD); 
    if (in_filename)
      {
      unlink(in_filename);
      }

  printf("  DONE! \n \n");

  if (log)
    fclose(log);
  return DRMS_SUCCESS;
}    // 

/*
 *
 *                     FUNCTIONS
 *
 *          img2sphere()
 *
 *  Written: Rick Bogart 
 *
*/
int img2sphere (double x, double y, double ang_r, double latc, double lonc,
    double pa, double *p_rho, double *p_lat, double *p_lon, double *p_sinlat,
    double *p_coslat, double *p_sig, double *p_mu, double *p_chi)
  {
	/*
 *  Taken from cartography.c, written by Rick Bogart
 *  Map projected coordinates (x, y) to (lon, lat) and (rho | sig, chi)
 *  
 *  Arguments:
 *    x }           Plate locations, in same units as the image radius, relative
 *    y }               to the image center
 *    ang_r         Apparent semi-diameter of sun (angular radius of sun at
 *                      the observer's tangent line)
 *    latc          Latitude of disc center, uncorrected for light travel time
 *    lonc          Longitude of disc center
 *    pa            Position angle of solar north on image, measured eastward
 *                      from north (sky coordinates)
 *  Return values:
 *    rho           Angle point:sun_center:observer
 *    lon           Heliographic longitude
 *    lat           Heliographic latitude
 *    sinlat        sine of heliographic latitude
 *    coslat        cosine of heliographic latitude
 *    sig           Angle point:observer:sun_center
 *    mu            cosine of angle between the point:observer line and the
 *                      local normal
 *    chi           Position angle on image measured westward from solar
 *                      north
 *
 *  All angles are in radians.
 *  Return value is 1 if point is outside solar radius (in which case the
 *    heliographic coordinates and mu are meaningless), 0 otherwise.
 *  It is assumed that the image is direct; the x or y coordinates require a
 *    sign change if the image is inverted.
 *
 * modified to ignore setting return values for NULL pointers.
 *
 */


  double w_rho, w_lat, w_lon, w_sinlat, w_coslat, w_sig, w_mu, w_chi;
  double *rho= &w_rho, *lat= &w_lat, *lon= &w_lon, *sinlat= &w_sinlat, *coslat= &w_coslat, *sig= &w_sig, *mu= &w_mu, *chi= &w_chi;
  static double ang_r0 = 0.0, sinang_r = 0.0, tanang_r = 0.0;
  static double latc0 = 0.0, coslatc = 1.0, sinlatc = 0.0;
  double cosr, sinr, sinlon, sinsig;

  if (p_rho) rho=p_rho;
  if (p_lat) lat=p_lat;
  if (p_lon) lon=p_lon;
  if (p_sinlat) sinlat=p_sinlat;
  if (p_coslat) coslat=p_coslat;
  if (p_sig) sig=p_sig;
  if (p_mu) mu=p_mu;
  if (p_chi) chi=p_chi;

  if (ang_r != ang_r0) {
      sinang_r = sin(ang_r);
      tanang_r = tan(ang_r);
      ang_r0 = ang_r;
    }

  if (latc != latc0) {
      sinlatc = sin(latc);
      coslatc = cos(latc);
      latc0 = latc;
    }

  // position angle
  *chi = atan2 (x, y) + pa;
  while (*chi > 2 * M_PI) *chi -= 2 * M_PI;
  while (*chi < 0.0) *chi += 2 * M_PI;
       /*  Camera curvature correction, no small angle approximations  */
  *sig = atan (hypot (x, y) * tanang_r);
  sinsig = sin (*sig);
  *rho = asin (sinsig / sinang_r) - *sig;
if (*sig > ang_r) return (-1);
  *mu = cos (*rho + *sig);
  sinr = sin (*rho);
  cosr = cos (*rho);
  *sinlat = sinlatc * cosr + coslatc * sinr * cos (*chi);
  *coslat = sqrt (1.0 - *sinlat * *sinlat);
  *lat = asin (*sinlat);

  sinlon = sinr * sin (*chi) / *coslat;
  *lon = asin (sinlon);
  if (cosr < (*sinlat * sinlatc)) *lon = M_PI - *lon;
  *lon += lonc;
  while (*lon < 0.0) *lon += 2 * M_PI;
  while (*lon >= 2 * M_PI) *lon -= 2 * M_PI;
  return (0);
}
/*
 *
 *          sphere2img()
 *
 *  Written: Rick Bogart 
 *
*/

int sphere2img (double lat, double lon, double latc, double lonc,
        double *x, double *y, double xcenter, double ycenter,
        double rsun, double peff, double ecc, double chi,
        int xinvrt, int yinvrt) {
/*
 *  Taken from cartography.c, written by Rick Bogart
 *  Perform a mapping from heliographic coordinates latitude and longitude
 *    (in radians) to plate location on an image of the sun.  The plate
 *    location is in units of the image radius and is given relative to
 *    the image center.  The function returns 1 if the requested point is
 *    on the far side (>90 deg from disc center), 0 otherwise.
 *
 *  Arguments:
 *      lat         Latitude (in radians)
 *      lon         Longitude (in radians)
 *      latc        Heliographic latitude of the disc center (in radians)
 *      lonc        Heliographic longitude of the disc center (in radians)
 *      *x }        Plate locations, in same units as the image radius, NOT relative
 *      *y }               to the image center
 *      xcenter }   Plate locations of the image center, in units of the
 *      ycenter }     image radius, and measured from an arbitrary origin
 *                    (presumably the plate center or a corner)
 *      rsun        Apparent semi-diameter of the solar disc, in plate
 *                    coordinates
 *      peff        Position angle of the heliographic pole, measured
 *                    eastward from north, relative to the north direction
 *                    on the plate, in radians
 *      ecc         Eccentricity of the fit ellipse presumed due to image
 *                    distortion (no distortion in direction of major axis)
 *      chi         Position angle of the major axis of the fit ellipse, 
 *                    measure eastward from north,  relative to the north
 *                    direction on the plate, in radians (ignored if ecc = 0)
 *      xinvrt}     Flag parameters: if not equal to 0, the respective
 *      yinvrt}       coordinates on the image x and y are inverted
 *
 *  The heliographic coordinates are first mapped into the polar coordinates
 *    in an orthographic projection centered at the appropriate location and
 *    oriented with north in direction of increasing y and west in direction
 *    of increasing x.  The radial coordinate is corrected for foreshortening
 *    due to the finite distance to the Sun. If the eccentricity of the fit
 *    ellipse is non-zero the coordinate of the mapped point is proportionately
 *    reduced in the direction parallel to the minor axis.
 *
 *  Bugs:
 *    The finite distance correction uses a fixed apparent semi-diameter
 *    of 16'01'' appropriate to 1.0 AU.  In principle the plate radius could
 *    be used, but this would require the plate scale to be supplied and the
 *    correction would probably be erroneous and in any case negligible.
 *
 *    The ellipsoidal correction has not been tested very thoroughly.
 *
 *    The return value is based on a test which does not take foreshortening
 *    into account.
 */
  static double sin_asd = 0.004660, cos_asd = 0.99998914;
                                                   /*  appropriate to 1 AU  */
  static double last_latc = 0.0, cos_latc = 1.0, sin_latc = 0.0;
  double r, cos_cang, xr, yr;
  double sin_lat, cos_lat, cos_lat_lon, cospa, sinpa;
  double squash, cchi, schi, c2chi, s2chi, xp, yp;
  int hemisphere;

  if (latc != last_latc) {
      sin_latc = sin (latc);
      cos_latc = cos (latc);
      last_latc = latc;
    }
  sin_lat     = sin (lat);
  cos_lat     = cos (lat);
  cos_lat_lon = cos_lat * cos (lon - lonc);

  cos_cang   = sin_lat * sin_latc + cos_latc * cos_lat_lon;
  hemisphere = (cos_cang < 0.0) ? 1 : 0;
  r = rsun * cos_asd / (1.0 - cos_cang * sin_asd);
  xr = r * cos_lat * sin (lon - lonc);
  yr = r * (sin_lat * cos_latc - sin_latc * cos_lat_lon);
                        /*  Change sign for inverted images  */
  if (xinvrt) xr *= -1.0;
  if (yinvrt) yr *= -1.0;
                        /*  Correction for ellipsoidal squashing of image  */
  if (ecc > 0.0 && ecc < 1.0) {
    squash = sqrt (1.0 - ecc * ecc);
    cchi = cos (chi);
    schi = sin (chi);
    s2chi = schi * schi;
    c2chi = 1.0 - s2chi;
    xp = xr * (s2chi + squash * c2chi) - yr * (1.0 - squash) * schi * cchi;
    yp = yr * (c2chi + squash * s2chi) - xr * (1.0 - squash) * schi * cchi;
    xr = xp;
    yr = yp;
  }

  cospa = cos (peff);
  sinpa = sin (peff);
  *x = xr * cospa - yr * sinpa;
  *y = xr * sinpa + yr * cospa;

  *y += ycenter;
  *x += xcenter;

  return (hemisphere);
}


// Generate explicit recordset list of closest good record to desired grid
// First get vector of times and quality
// Then if vector is not OK, quit.
// then: make temp file to hold recordset list
//       start with first time to define desired grid, 
//       make array of desired times.
//       make empty array of recnums
//       search vector for good images nearest desired times
//       for each found time, write record query
char *get_input_recset(DRMS_Env_t *drms_env, char *in, TIME cadence)
  {
  DRMS_Array_t *data;
  TIME t_start, t_stop, t_now, t_want, t_diff;
  int status;
  int nrecs, irec;
  int nslots, islot;
  long long *recnums;
  TIME *slot_times;
  TIME *t_this;
  double *drecnum, *dquality;
  int quality;
  long long recnum;
  char keylist[DRMS_MAXQUERYLEN];
  static char filename[100];
  char *tmpdir;
  FILE *tmpfile;
  char *lbracket;
  char seriesname[DRMS_MAXQUERYLEN];

  sprintf(keylist, "T_OBS,QUALITY,recnum");
fprintf(stderr,"Get vector for: %s\n",in);
  data = drms_record_getvector(drms_env, in, keylist, DRMS_TYPE_DOUBLE, 0, &status);
  if (!data || status)
	{
	fprintf(stderr, "getkey_vector failed status=%d\n", status);
	return(NULL);
	}
  nrecs = data->axis[1];
fprintf(stderr,"A nrecs=%d\n",nrecs);
  irec = 0;
  t_this = (TIME *)data->data;
  dquality = (double *)data->data + 1*nrecs;
  drecnum = (double *)data->data + 2*nrecs;
  t_start = t_this[0];
  t_stop = t_this[nrecs-1];
  nslots = (t_stop - t_start + cadence/2)/cadence;
  recnums = (long long *)malloc(nslots*sizeof(long long));
  slot_times = (TIME *)malloc(nslots*sizeof(TIME));
fprintf(stderr,"nslots=%d\n",nslots);
  islot = 0;
  t_want = t_start;
  t_diff = 1.0e8; // 3+ years
  for (irec = 0; irec<nrecs; irec++)
      {
      t_now = t_this[irec];
      quality = (int)dquality[irec] & 0xFFFFFFFF;
      recnum = (long long)drecnum[irec];
      if (quality < 0)
        continue;
      if (fabs(t_now - t_want) <= t_diff)
        {
        slot_times[islot] = t_now;
        recnums[islot] = recnum;
        }
// fprintf(stderr,"irec=%d, t_want=%f, t_now=%f, t_diff=%f quality=%x, recnum=%lld, slot=%d, slottime=%f, slotrecnum=%lld\n",irec,t_want,t_now,t_diff,quality,recnum,islot,slot_times[islot],recnums[islot]);
      t_diff = fabs(t_now - t_want);
      if ( t_now > t_want) // past this slot, time to move to next slot
        {
        islot++;
        if (islot >= nslots)
           break;
        t_want = t_start + cadence * islot;
        slot_times[islot] = t_now;
        recnums[islot] = recnum;
        t_diff = fabs(t_now - t_want);
        }
      }
  if (islot+1 < nslots)
    nslots = islot+1;  // take what we got.
  strcpy(seriesname, in);
  lbracket = index(seriesname,'[');
  if (lbracket) *lbracket = '\0';
  tmpdir = getenv("TMPDIR");
  if (!tmpdir) tmpdir = "/tmp";
  sprintf(filename, "%s/hg_patchXXXXXX", tmpdir);
  mkstemp(filename);
  tmpfile = fopen(filename,"w");
  for (islot=0; islot<nslots; islot++)
    fprintf(tmpfile, "%s[:#%lld]\n", seriesname, recnums[islot]);
  fclose(tmpfile);
  free(recnums);
  free(slot_times);
  drms_free_array(data);
// char cmd[2000];sprintf(cmd,"echo %s\ncat %s",filename,filename);
// system(cmd);
  return(filename);
  }
