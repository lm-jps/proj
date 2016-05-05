/**
   @defgroup im_patch im_patch (image_patch) extracts a rectangular 
   region or patch from an input series.  The patch is tracked across
   the solar disk from time  t_start to t_stop.  The patch size in pixels is defined by either its
   Carrington rotation, longitude, latitude and size when it crosses central meridian or
   by its location in arcsec at a reference time.  If defined by Carrington coordinates the
   bounding box size is specified in degrees of longitude and latitude at disk center.  Else the box is
   specified by its reference image time and location in locunits.  No reprojection is done, the patch is
   a simple cutout of the original image.

   @ingroup a_programs

   @brief im_patch extracts a rectangular 
   image from an input series.  The patch is tracked across
   the solar disk from t_start to t_stop times.  The patch may be defined by Carrington or Syonyhurst coordinates or
   in arcsec or pixels.  The extracted patch is rectangular in original image pixels.


   @par Synopsis:

   @code
   im_patch in=<recordset> out=<series_out>  \
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
   If by position at a reference time, the position may be specified in ipixels, arcsec, or Stonyhurst coordinates.
   If the car_rot parameter is required if the Carrington location method is used.
   If the t_ref parameter is present, location is specified by a method given in <locunits>.
   The location is given in the "x" and "y" arguments where the units of x and y are determined by locunits.
   The choices for locunits are: "arcsec", "pixels", "stony", "carrlong".
   For "stony" (Stonyhurst) x is the CM distance in degrees, for "carrlong" x is Carrington longitude.
   The y value is in degrees of latitude for either "stony" or "carrlong".
   Arcsec locations are measured from disk center, pixel locations are measured from the lower left corner of
   the image, with the first pixel labeled (1,1).

   The size of the patch, given in the width and height parameters can be expressed in "pixels", "arcsec", or "degrees"
   as determined by the "boxunits" parameter.

   If the logfile parameter is present, a RecordSet query will be written to logfile for each image created  containing a query
   that will return that image.

   RequestID, if present, is the ID of an on-requst processing export request.

   If neither a limiting recordset specification nor t_start or t_stop are present, those values will be inferred from
   the on-disk time span for the center of the patch -90 degrees to +90 degrees from CM.
   If cadence is not specified explicitly on the command line or in the recordset spec, the full cadence of the dataset will be used.
   This module computes the beginning and ending time for tracking the 
   region identified at time = t_ref, x and y or by the Carrington coordinates of the box center.
   In the Carrington case the module extracts a rectangular region 
   with user defined size. the height and width of the box are specified in degrees
   of latitude and longitude with the pixel size of the box is computed from the
   projection of the box at CM.  The pixel size of the rectangular box remains constant 
   during tracking.  The box height defaults to the width.  The width defaults to 10 degrees.
  
   If the input recordset spec is exactly "[$]" that spec will be discarded and the time limits will be taken from t_start, t_stop, car_rot, and/or t_ref as appropriate.  This is to allow exports via jsoc_fetch to not require a bounding recordset to be specified.

   If t_start and/or t_stop are specified the referecne time or disk center location need not be
   included.
   
   If the input is only a seriesname with [$], it must have a prime key of type time.

   If the t_ref parameter is specified, there must be a non-missing image within +- 2 hours of t_ref.

   The output seriesname defaults to the input seriesname with a suffix of "_mod".

   The -t, 'no tracking' flag casues the extracted region to remain fixed with respect to disc center.  I.e. the Carrington
   tracking is disabled.  If locunits are "stony" or "carrlong" the patch center must be on the disk.
  
   The -c, 'crop' flag causes the off-limb data to be replaced with MISSING (NaNs for floating types).

   The -h "header keywords" flag causes the metadata keywords to be included in all of the output FITS files.

   If -r or -R are specified, the images will be sub-pixel registered to the target location.  If -r or -R are given
   along with -t (no tracking) then -r means register to the first frame and -R means register to Sun center.
 
   -f, -F flags and FDS are special rarely used parameters that allow tracking to a location given in
   a SDO flight dynamics system generated trajectory of e.g. a transitting moon or planet.  The -f or -F flags
   in this case will generate a fake object with value either 0 or datamax. For example, the Venus transit
   locations are found in sdo.fds[2012.06.05_00:00:00_UTC][SOLAR_TRANSIT][S][2]

   If the input series is AIA lev1 or HMI lev1 where the time slots are not fully populated then special
   code is used to produce the requested cadence.  In these cases the DRMS "time@cadence" recordset notation
   will fail so use the cadence parameter instead.

   @par Flags:
   -A   Generate output for all segments.
   -t   Disable tracking.
   -h   Include header keywords
   -r   register to first frame
   -R   register to Sun center
   -f   make black fake tracked object
   -F   make light fake tracked object


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
   im_patch in='hmi.M_45s' locunits=carrlong boxunits=pixels car_rot=2009 x=295 y=5 width=30 height=20

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

int get_tracking_xy(char *FDSfile, TIME want, double *xp, double *yp, double *radius);
void HeliographicLocation(TIME t, int *crot, double *L, double *B);
TIME HeliographicTime(int crot, double L);

char *module_name = "im_patch";

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
    {ARG_STRING, "requestid", "NOTSPECIFIED", "RequestID if im_patch call originated in an export request."},
    {ARG_INT, "car_rot", "-1", "Carrington Rotation when the region crosses CM"},
    {ARG_FLOAT, "width", "0", "width of box in boxunits (conversion to pixels made when at CM)"},
    {ARG_FLOAT, "height", "0", "height of box in boxunits (conversion to pixels made when at CM)"},
    {ARG_STRING, "boxunits", "NOTSPECIFIED", "units of patch, 'pixels', 'arcsecs', or 'degrees'"},
    {ARG_TIME, "t_start", "JD_0", "Start time, defaults to time at 90E"},
    {ARG_TIME, "t_stop", "JD_0", "End time, defauolts to 90W"},
    {ARG_TIME, "cadence", "NOTSPECIFIED", "Cadence of product, defaults to input cadence."},
    {ARG_TIME, "t_ref", "JD_0", "Time for which x and y apply, implies ref image."},
    {ARG_STRING, "locunits", "NOTSPECIFIED", "Location units in 'pixels', 'arcsec', 'stony', or 'carrlong'"},
    {ARG_STRING, "where", "NOTSPECIFIED", "Additional 'where' clause if needed"},
    {ARG_FLOAT, "x", "0", "Location of extract box center in pixels in array (1,1), arcsec from Sun center, degrees from Sun center, or Carrington longitude"},
    {ARG_FLOAT, "y", "0", "Location of extract box center, same units as x"},
    {ARG_FLAG, "A", NULL, "Generate output for all input segments"},
    {ARG_FLAG, "t", "0", "Disable Carrington rate tracking"},
    {ARG_FLAG, "h", "0", "Export keywords, i.e. include full metadata in FITS files"},
    {ARG_FLAG, "c", "0", "Crop at limb, useful for HMI M data."},
    {ARG_FLAG, "f", "0", "FAKE tracking target added, val=0"},
    {ARG_FLAG, "F", "0", "FAKE tracking target added, val=datamax"},
    {ARG_FLAG, "r", "0", "Register to fractional CRPIXi, all frames to first frame round(crpix)"},
    {ARG_STRING, "FDS", "NOTSPECIFIED", "FDS file for tracking coords"},
    {ARG_END}
  };

int img2sphere (double x, double y, double ang_r, double latc, double lonc,
        double pa, double *rho, double *lat, double *lon, double *sinlat,
        double *coslat, double *sig, double *mu, double *chi);
int sphere2img (double lat, double lon, double latc, double lonc,
        double *x, double *y, double xcenter, double ycenter,
        double rsun, double peff, double ecc, double chi,
        int xinvrt, int yinvrt);
int image_magrotate(void *, int nin, int min, int data_type_input, float theta, float mag,
    float dx, float dy, void **outarray, int *nx, int *ny, int regridtype_input, int stretch_lines);

// experimental make recordset list at cadence
//  OLD VERSION char *get_input_recset(DRMS_Env_t *env, char *in, TIME cadence, char *timekeyname);
char *get_input_recset(DRMS_Env_t *drms_env, char *inQuery);

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
  DRMS_RecordSet_t *inRS = NULL;
  DRMS_RecordSet_t *outRS = NULL;
  DRMS_Record_t *inRec = NULL;
  DRMS_Record_t *outRec = NULL;
  DRMS_Record_t *inTemplate = NULL;
  DRMS_Record_t *outTemplate = NULL;
  DRMS_Segment_t *inSeg = NULL;
  DRMS_Segment_t *outSeg = NULL;
  DRMS_Array_t *inArray = NULL;
  DRMS_Array_t *outArray = NULL;
  int i, ii, status = DRMS_SUCCESS, nrecs; 
  int  irec;
  int  OK_recs = 0;
  double center_x, center_y, crpix1, crpix2, x0, y0; 
  double rsun_ref, dsun_obs, rsun, rsun_rad, rsunpix;
  double r2;
  double this_x0;
  double this_y0;
  double crln_obs, crlt_obs;
  double crln_obs_rad, crlt_obs_rad;
  double crln_rad, crlt_rad;
  double pa_rad;
  double cdelt;
  double llx, lly;
  int inAxis[2];
  int this_car_rot;
  char *ctype1, *ctype2;
  TIME t_rec;
  double box_x, box_y;
  double crln, crlt;
  char outseries[DRMS_MAXSERIESNAMELEN];
  char inseries[DRMS_MAXSERIESNAMELEN];
  char extractedSeries[DRMS_MAXSERIESNAMELEN];
  int is_mod;
  char testseries[DRMS_MAXSERIESNAMELEN];
  char inQuery[DRMS_MAXQUERYLEN];
  char in[DRMS_MAXQUERYLEN];

  const char *ingiven = params_get_str(params, "in");
  char *inparam;
  char *lbracket, *rbracket;
  char *moreQuery = NULL;
  const char *outparam = params_get_str(params, "out");
  const char *logfile = params_get_str(params, "log");
  const char *requestid = params_get_str(params, "requestid");
  const char *where = params_get_str(params, "where");
  const char *FDSfile = params_get_str(params, "FDS");
  TIME t_start = params_get_time(params, "t_start");
  TIME t_stop = params_get_time(params, "t_stop");
  const char *cadence_str = params_get_str(params, "cadence");
  double cadence;
  TIME t_ref = params_get_time(params, "t_ref");
  double width = params_get_double(params, "width");
  double height = params_get_double(params, "height");
  int pixwidth, pixheight;
  double midx, midy;
  int car_rot = params_get_int(params, "car_rot");
  double x = params_get_double(params, "x");
  double y = params_get_double(params, "y");
  const char *boxunits = params_get_str(params, "boxunits");
  const char *locunits = params_get_str(params, "locunits");
  int NoTrack = params_isflagset(params, "t");
  int wantFAKEwhite = params_isflagset(params, "F");
  int wantFAKE = params_isflagset(params, "f") || wantFAKEwhite;
  int do_register = params_isflagset(params, "r");
  int export_keys = params_isflagset(params, "h");
  int do_crop = params_isflagset(params, "c");
  TIME tNotSpecified = sscan_time("JD_0");
  int do_reftime = 0;  // Target specified by reftime vs Carr rotation mode.
  double crvalx = 0.0;
  double crvaly = 0.0;
  double crota, sina, cosa;
  double pa, deltlong;
  double target_x, target_y;  // used for NoTrack option
  int firstimage = 1;
  int boxtype, loctype;
  char *timekeyname;
  TIME t_step, t_epoch;
  int npkeys;
  int times_source = TIMES_RECSET;
  char timebuf[20];
  char *in_filename = NULL;
  double trackx, tracky;
  double track_radius;
  double crpix1_0, crpix2_0;
  int register_padding = 5; // later this should be based on crota2 and box size.
  char *tmpstr = NULL;
  int hasSegList = 0;
  int doingAllSegs = params_isflagset(params, "A");

  FILE *log = NULL;

  if (strcmp(FDSfile,"NOTSPECIFIED")==0)
    FDSfile = NULL;

  if (strncasecmp(locunits, "arcsec", 6) == 0) loctype = LOCARCSEC;
  else if (strncasecmp(locunits, "pixels", 3) == 0) loctype = LOCPIXELS;
  else if (strncasecmp(locunits, "stony", 5) == 0) loctype = LOCSTONY;
  else if (strncasecmp(locunits, "carrlong", 4) == 0) loctype = LOCCARR;
  else loctype = LOCBAD;

  if (strncasecmp(boxunits, "arcsec", 6) == 0) boxtype = BOXARCSEC;
  else if (strncasecmp(boxunits, "pixels", 3) == 0) boxtype = BOXPIXELS;
  else if (strncasecmp(boxunits, "degrees", 3) == 0) boxtype = BOXDEGREE;
  else loctype = BOXBAD;

  if (loctype == LOCCARR)
    {
    if (car_rot < 0) DIE("Carrington rotation number must be provided for locunits=carrlong");
    // if (NoTrack) DIE("Can not use locunits=carrlong for no-tracking mode");
    do_reftime = 0;
    }
  else 
    {
    if (t_ref < 0)
      DIE("t_ref must be provided for locunits!=carrlong");
    do_reftime = 1;
    }

  if (loctype==LOCBAD)
    DIE2("Location units not detected, locunits", locunits)
  if (boxtype==BOXBAD)
    DIE2("patch dimensions not understood,", boxunits);

  if (strcasecmp(logfile, "NOTSPECIFIED") != 0)
    {
    log = fopen(logfile, "w");
    if (!log) DIE2("Can not create log file.",logfile);
    }

  if (strcasecmp(where, "NOTSPECIFIED") == 0 || *where == '\0')
    {
    where = "";
    }
  else
    {
    char wherework[4096];
    sprintf(wherework,"[? %s ?]", where);
    where = strdup(wherework);
    }

    /* We need to extract the series name from the input record-set specification. There */

    /* Determine if there is a segment list specifier. I believe this is true if the last non-ws char of the specification is a '}'. 
     * Of course, the specification could be invalid, but in that event the drms_open_records() call will fail and catch the error and
     * the module run will abort. */
     char *testSegList = NULL;
     char *pCh = NULL;
     
     testSegList = rindex(ingiven, '}');
     if (testSegList)
     {
        for (pCh = testSegList + 1, hasSegList = 1; *pCh; pCh++)
        {
            if (!isspace(*pCh))
            {
                /* There is a non-ws char after '}' - this is not a valid seglist. */
                hasSegList = 0;
                break;
            }
        }
     }

  inparam = strdup(ingiven);
  if (strcasecmp(inparam, "NOTSPECIFIED") == 0) DIE("Input series must be specified.");
  
  /* This program works only if the series' first prime-key constituent is a time keyword. Then, we want to
   * ignore the first prime key and use the t_start/t_stop argument to select records instead. 
   * I believe the program will not work with record-set specifications that contain more than one
   * subset either. So...
   * 1. Extract the series from the record-set specification. 
   * 2. Bail if there is more than one record-set subset.
   * 3. Obtain the first keyword from the prime-key set. 
   * 4. Bail if the first keyword is not a time keyword.
   * 5. Remove the time filter provided by the user in the record-set specification.
   */
   
   /* Use drms_record_parserecsetspec() to extract the series and filter strings. The previous 
    * code wasn't quite working 100% . */
    char *allvers = NULL;
    char **sets = NULL;
    DRMS_RecordSetType_t *settypes = NULL; /* a maximum doesn't make sense */
    char **snames = NULL;
    char **filts = NULL;
    int nsets = 0;
    DRMS_RecQueryInfo_t rsinfo; /* Filled in by parser as it encounters elements. */
    DRMS_Keyword_t *firstPKey = NULL;

    if (drms_record_parserecsetspec(inparam, &allvers, &sets, &settypes, &snames, &filts, &nsets, &rsinfo) == DRMS_SUCCESS)
    {
        if (nsets > 1)
        {
            status = DRMS_ERROR_INVALIDDATA;
            DIE("im_patch supports single-set record-set specifications only.");
        }
    }
    else
    {
        status = DRMS_ERROR_INVALIDDATA;
        DIE("Invalid record-set specification provided.");
    }

    snprintf(inseries, sizeof(inseries), "%s", snames[0]);
    lbracket = index(inparam, '[');
    
    inTemplate = drms_template_record(drms_env, inseries, &status);
    if (status || !inTemplate) DIE2("Input series can not be found: ", inseries);
    
    if (drms_keyword_isindex(inTemplate->seriesinfo->pidx_keywords[0]))
    {
        firstPKey = drms_keyword_slotfromindex(inTemplate->seriesinfo->pidx_keywords[0]);
    }
    else
    {
        firstPKey = inTemplate->seriesinfo->pidx_keywords[0];
    }

    if (firstPKey->info->type != DRMS_TYPE_TIME)
    {
        status = DRMS_ERROR_INVALIDDATA;
        DIE("The first prime-key constituent keyword must be of type TIME.");
    }
    
    if (filts && filts[0])
    {
        if (strlen(filts[0]) == 3 && filts[0][1] == '$')
        {
            /* The user specified the last (according to time) record. Remove that filter. */
            moreQuery = lbracket + 3;
            *lbracket = '\0';
            lbracket = NULL;
            
            if (t_start == tNotSpecified || t_stop == tNotSpecified)
                times_source = TIMES_IMPLICIT;
            else
                times_source = TIMES_GIVEN;
        }
        else
        {
            rbracket = index(lbracket, ']');
            moreQuery = rbracket ? rbracket + 1 : lbracket;
        }
    }
    else
    {
        /* No filter provided at all. */
        if (t_start == tNotSpecified || t_stop == tNotSpecified)
            times_source = TIMES_IMPLICIT;
        else
            times_source = TIMES_GIVEN;
    }
    
    drms_record_freerecsetspecarr(&allvers, &sets, &settypes, &snames, &filts, nsets);
    
// XXXXX Get box location in Carrington Coords for all location types, also initial center_X, center_y in NoTrack case  XXXXX
  if (do_reftime) // Arc-sec, pixel, or Stonyhurst specification, get ref image information
    { // the image for ref_time is supposed to be present in the series.
fprintf(stderr,"doing reftime\n");
    int okrec;
    char t_ref_pre[100];
    char t_ref_post[100];
    TIME search_width = 15;
    while (search_width < 10000)
      { 
      sprint_at(t_ref_pre, t_ref-search_width);
      sprint_at(t_ref_post, t_ref+search_width);
      sprintf(in, "%s[%s-%s][? QUALITY >= 0 ?]%s%s", inseries, t_ref_pre, t_ref_post, moreQuery, where);
fprintf(stderr,"searchwidth=%f, t_ref query is: %s\n",search_width,in);
      inRS = drms_open_records(drms_env, in, &status);
      if (!status && inRS->n > 0)
        {
        int irec, nrecs = inRS->n;
        TIME tdiff = 100000;
        for (irec=0, okrec=-1; irec<nrecs; irec++) // find record close to t_ref
          {
          TIME newdiff;
          inRec = inRS->records[irec];
fprintf(stderr,"irec=%d, tdiff=%lf\n",irec,tdiff);
          if ((newdiff=(fabs(drms_getkey_time(inRec,"T_OBS",NULL) - t_ref))) > tdiff)
            break;
          okrec = irec;
          tdiff = newdiff;
          }
        }
      if (okrec >= 0)
        break;
      search_width *= 2;
      drms_close_records(inRS, DRMS_FREE_RECORD);
      }
    if (search_width >= 10000)
      DIE2("No input data found within 2-hours of t_ref",in);
    inRec = inRS->records[okrec];
    this_car_rot = drms_getkey_int(inRec, "CAR_ROT", &status); TEST_PARAM("CAR_ROT");
    crlt_obs = drms_getkey_double(inRec, "CRLT_OBS", &status); TEST_PARAM("CRLT_OBS");
    crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status); TEST_PARAM("CRLN_OBS");
    crpix1 = drms_getkey_double(inRec, "CRPIX1", &status); TEST_PARAM("CRPIX1");
    x0 = crpix1 - 1;
    crpix2 = drms_getkey_double(inRec, "CRPIX2", &status); TEST_PARAM("CRPIX2");
    y0 = crpix2 - 1;
    crvalx = drms_getkey_double(inRec, "CRVAL1", &status); TEST_PARAM("CRVAL1");
    crvaly = drms_getkey_double(inRec, "CRVAL2", &status); TEST_PARAM("CRVAL2");
    crota = drms_getkey_double(inRec, "CROTA2", &status); if (status) crota = 0.0; // WCS default
    pa = -crota;
    crlt_obs_rad = Deg2Rad * crlt_obs;
    crln_obs_rad = Deg2Rad * crln_obs;
    pa_rad = Deg2Rad * pa;
    sina = sin(crota*Deg2Rad);
    cosa = cos(crota*Deg2Rad);
    tmpstr = drms_getkey_string(inRec, "CTYPE1", &status);
    ctype1 = strdup(tmpstr); TEST_PARAM("CTYPE1");
    if (tmpstr)
      {
      free(tmpstr);
      }
    if (strcmp(ctype1, "HPLN-TAN") != 0)
      DIE2("CTYPE1 not HPLN-TAN as required, is: ", ctype1);
    if (ctype1)
      {
      free(ctype1);
      ctype1 = NULL;
      }
    tmpstr = drms_getkey_string(inRec, "CTYPE2", &status);
    ctype2 = strdup(tmpstr); TEST_PARAM("CTYPE2");
    if (tmpstr)
      {
      free(tmpstr);
      }
    if (strcmp(ctype2, "HPLT-TAN") != 0)
      DIE2("CTYPE2 not HPLT-TAN as required, is: ", ctype2);
    if (ctype2)
      {
      free(ctype2);
      ctype2 = NULL;
      }
    rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
    if (status)
      rsun_ref = 6.96e8;
    dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status); TEST_PARAM("DSUN_OBS");
    rsun_rad = asin(rsun_ref/dsun_obs); 
    rsun = rsun_rad*Rad2arcsec; 
    cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
    rsunpix = rsun/cdelt;
    if (loctype == LOCPIXELS || loctype==LOCARCSEC)
      { /* box ref location is in arcsec - find Carrington equivalent */
      // Note x and y and PIX_X and PIX_Y have first pixel at (1,1)
      // center_x and center_y are center of target box where first pixel is (0,0)
      if (loctype==LOCARCSEC)
        {
        center_x = PIX_X(x,y) - 1;
        center_y = PIX_Y(x,y) - 1;
        }
      else
        {
        center_x = x - 1;
        center_y = y - 1;
        }
      // center_x and center_y are box center in pixels, first pixel at (0,0)
      inSeg = drms_segment_lookupnum(inRec, 0);
      inAxis[0] = inSeg->axis[0];
      inAxis[1] = inSeg->axis[1];
      if ( img2sphere((center_x - x0)/rsunpix, (center_y -y0)/rsunpix, rsun_rad, crlt_obs_rad, crln_obs_rad, pa_rad,
           NULL, &crlt_rad, &crln_rad, NULL, NULL, NULL, NULL, NULL) < 0)
          fprintf(stderr, "Starting location is off the solar disk.");
      crln = crln_rad * Rad2Deg;
      crlt = crlt_rad * Rad2Deg;
fprintf(stderr,"at do reftime, center_x=%f\n", center_x);
      }
    else /* loctype==LOCSTONY */
      {
      crlt = y;
      crln = crln_obs + x;
      sphere2img(crlt/Rad2Deg, crln/Rad2Deg, crlt_obs_rad, crln_obs_rad, &center_x, &center_y, x0, y0, rsunpix, pa_rad, 0, 0, 0, 0);
      }
    deltlong = crln - crln_obs;
    car_rot = this_car_rot;
    if (crln >= 360.0)
      { car_rot--; crln -= 360.0; }
    else if (crln < 0)
      { car_rot++; crln += 360.0; }
    drms_close_records(inRS, DRMS_FREE_RECORD);
    }
  else // Carrington specification
    {
    crln = x;
    crlt = y;
    if (crln < 0) DIE("Box longitude must be specified.");
    if (crlt < -990) DIE("box latitude must be specified.");
    sphere2img(crlt/Rad2Deg, crln/Rad2Deg, crlt_obs_rad, crln_obs_rad, &center_x, &center_y, x0, y0, rsunpix, pa_rad, 0, 0, 0, 0);
    }
  crln_rad = crln * Deg2Rad;
  crlt_rad = crlt * Deg2Rad;

  // Now we have crln, crlt, crln_rad, crlt_rad, car_rot for the target box. 
  // unless the box center is off the disk.
  // For NoTrack case we have center_x and center_y in pixels

  // Get implied time limits from car_rot and crln
  t_ref = HeliographicTime(car_rot, crln);
  if (times_source == TIMES_IMPLICIT)
      {
      if (t_start == tNotSpecified)
        {
        if (NoTrack) DIE("Must have t_start or explicit times for NoTrack");
        t_start = HeliographicTime(car_rot, crln + 90);
        }
      if (t_stop == tNotSpecified)
        {
        if (NoTrack) DIE("Must have t_stop or explicit times for NoTrack");
        t_stop = HeliographicTime(car_rot, crln - 90);
        }
      }
// XXXXXXXXXXXXXXXXX End of get target box location

// XXXXXXXXXXXXXXXX get input seriesname and output seriesname
// first, get input and output series names.
  is_mod = 0;
  if (strcasecmp(outparam, "NOTSPECIFIED") == 0)
    {
    strncpy(outseries, inseries, DRMS_MAXSERIESNAMELEN);
    strncat(outseries, "_mod", DRMS_MAXSERIESNAMELEN);
    is_mod = 1;
    }
  else
    {
    strncpy(outseries, outparam, DRMS_MAXSERIESNAMELEN);
    strncpy(testseries, inseries, DRMS_MAXSERIESNAMELEN);
    strcat(testseries, "_mod");
    if (strcmp(testseries, outseries) == 0)
      is_mod = 1;
    }

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
  
  if (strcasecmp(cadence_str, "NOTSPECIFIED") == 0)
    cadence = t_step;
  else
    {
    int ratio;
    double initial_cadence;
    initial_cadence = atoinc((char *)cadence_str);
    ratio = round(initial_cadence/t_step);
    cadence = ratio * t_step;
    if (cadence != initial_cadence)
      fprintf(stderr,"Cadence rounded from %f to %f, since t_step == %f\n", initial_cadence, cadence, t_step);
    }

  // Finally, get input recordset spec
  if (lbracket && times_source == TIMES_RECSET) // RecordSet query specified
    strncpy(in, inparam, DRMS_MAXQUERYLEN);
  else // need to generate query from limit times
    { 
    char t_start_text[100], t_stop_text[100], cadence_text[100];
    if (cadence > t_step)
      {
      sprintf(cadence_text,"@%0.0fs", cadence);
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
    sprint_at(t_start_text, t_start);
    sprint_at(t_stop_text, t_stop);
    sprintf(in, "%s[%s-%s%s]%s%s", inseries, t_start_text, t_stop_text, cadence_text, where, (moreQuery ? moreQuery : ""));
    }

  if (inparam)
  {
     free(inparam);
     inparam = NULL;
  }

  // Check for specified epoch or cadence specified and not compact slotted series, replace query if needed
  if (cmdparams_exists(&cmdparams, "epoch") || (cadence > 1.0 &&
     (strncmp(inseries, "aia.lev1", 8)==0 || strncmp(inseries, "hmi.lev1", 8)==0 )))
    { // special case for AIA and HMI lev1 slots
    char *new_in;
    new_in = get_input_recset(drms_env, in);
    if (new_in != in)
      {
      strcpy(in, new_in);
      in_filename = new_in+1;
      }
    }

// XXXXXXXXXXXXX End of get input and output series information

// XXXXXXXXXXXXX Now read data and extract boxes.

  inRS = drms_open_records(drms_env, in, &status);
  if (status || inRS->n == 0)
     {
     fprintf(stderr,"Query is: %s\n",in);
     DIE("No input data found");
     }
  nrecs = inRS->n;

  // extract patches from each record
  int nextRec = 0;
  for (irec = 0; irec < nrecs; irec ++)
    {
    double use_bzero, use_bscale;
    char history[4096];
    
    nextRec = 0;
    
    *history = '\0';
    inRec = inRS->records[irec];
    if (status || !inRec) DIE("Record read failed.");
    
    
    TIME trec = drms_getkey_time(inRec, timekeyname, &status); TEST_PARAM(timekeyname);
    if (t_start != tNotSpecified && trec < t_start)
      {
      continue;
      }
    if (t_stop != tNotSpecified && trec > t_stop)
      break;
    // skip records without images.
    if (drms_keyword_lookup(inRec, "QUALITY", 1))
      {
      int quality = drms_getkey_int(inRec, "QUALITY", &status);
      if (quality < 0)
         {
         continue;
         }
      }

    // get venus or other transit object location
    if (FDSfile)
      {
      TIME t_obs = drms_getkey_time(inRec, "T_OBS", NULL);
      if(get_tracking_xy((char *)FDSfile, t_obs, &trackx, &tracky, &track_radius))
        DIE("Failed to open FDS file");
      }
    
    // get coordinate information for this image
    this_car_rot = drms_getkey_int(inRec, "CAR_ROT", &status); TEST_PARAM("CAR_ROT");
    crlt_obs = drms_getkey_double(inRec, "CRLT_OBS", &status); TEST_PARAM("CRLT_OBS");
    crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status); TEST_PARAM("CRLN_OBS");
    crpix1 = drms_getkey_double(inRec, "CRPIX1", &status); TEST_PARAM("CRPIX1");
    x0 = crpix1 - 1;
    crpix2 = drms_getkey_double(inRec, "CRPIX2", &status); TEST_PARAM("CRPIX2");
    y0 = crpix2 - 1;
    crota = drms_getkey_double(inRec, "CROTA2", &status); TEST_PARAM("CROTA2");
    cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
    pa = -crota;
    
    int paIs180 = 0;
    if (fabs(fabs(pa)-180.0) < 1.0)
    {
        paIs180 = 1;
    }
    
    crlt_obs_rad = Deg2Rad * crlt_obs;
    crln_obs_rad = Deg2Rad * crln_obs;
    pa_rad = Deg2Rad * pa;
    sina = sin(crota*Deg2Rad);
    cosa = cos(crota*Deg2Rad);

fprintf(stderr,"T_OBS=%s, crpix1=%f, crpix2=%f\n", drms_getkey_string(inRec,"T_OBS", NULL), crpix1, crpix2);

//  XXXXXXXXXXXXXXXXX get coordinate mapping info using keywords from first record
    if (firstimage)
      {
      firstimage = 0;
      // firstimage values for register to first image option.
      int idx, idy;
      idx = crpix1 + 0.5;
      idy = crpix2 + 0.5;
      crpix1_0 = idx;
      crpix2_0 = idy;

      if (loctype==LOCCARR)
        {
           tmpstr = drms_getkey_string(inRec, "CTYPE1", &status);
            ctype1 = strdup(tmpstr); TEST_PARAM("CTYPE1");
        if (tmpstr)
        {
           free(tmpstr);
        }

        tmpstr = drms_getkey_string(inRec, "CTYPE2", &status);
        if (strcmp(ctype1, "HPLN-TAN") != 0) DIE2("CTYPE1 not HPLN-TAN as required, is: ", ctype1);
        if (ctype1)
        {
           free(ctype1);
           ctype1 = NULL;
        }
        ctype2 = strdup(tmpstr); TEST_PARAM("CTYPE2");
        if (tmpstr)
        {
           free(tmpstr);
        }

        if (strcmp(ctype2, "HPLT-TAN") != 0) DIE2("CTYPE2 not HPLT-TAN as required, is: ", ctype2);
        if (ctype2)
        {
           free(ctype2);
           ctype2 = NULL;
        }
        rsun_ref = drms_getkey_double(inRec, "RSUN_REF", &status);
        if (status) rsun_ref = 6.96e8;
        dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status); TEST_PARAM("DSUN_OBS");
        rsun_rad = asin(rsun_ref/dsun_obs); 
        rsun = rsun_rad*Rad2arcsec; 
        cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
        // in principle use deltlong to get time to use for correcting crlt_obs, ignore for now
        // get ur, ll coords of box at CM.
        rsunpix = rsun/cdelt;
        }

      if (boxtype == BOXDEGREE)
        {
        int nx, ny;
        double urx, ury;
        // get corners of sample box at zero longitude wrt center.
        sphere2img(Deg2Rad*(crlt+height/2), Deg2Rad*(width/2), crlt_obs_rad, 0.0,
          &urx, &ury, 0.0, 0.0, rsunpix, pa_rad, 0, 0, 0, 0);
        sphere2img(Deg2Rad*(crlt-height/2), -Deg2Rad*(width/2), crlt_obs_rad, 0.0,
          &llx, &lly, 0.0, 0.0, rsunpix, pa_rad, 0, 0, 0, 0);
        // llx, urx, lly, ury are image coordinates for Sun, want llx and lly
        // on plate so need to swap if image rotated.
        nx = urx - llx ;
        ny = ury - lly ;
        if (nx < 0) { nx = -nx; llx = urx; }
        if (ny < 0) { ny = -ny; lly = ury; }
        pixwidth = nx + 1;
        pixheight = ny + 1;
        }
      else if (boxtype==BOXARCSEC)
        {
        int nx, ny;
        llx = -width/(2.0*cdelt);
        lly = -height/(2.0*cdelt);
        pixwidth = 1 - 2 * llx;
        pixheight = 1 - 2 * lly;
        }
      else /* boxtype == BOXPIXELS */
        {
        pixwidth = round(width);
        pixheight = round(height);
        llx = -pixwidth / 2.0;
        lly = -pixheight / 2.0;
        }
      }

    if (!NoTrack)
        {
        if (FDSfile)
          {
          center_x = PIX_X(trackx,tracky) - 1;
          center_y = PIX_Y(trackx,tracky) - 1;
          sprintf(history+strlen(history), "Image center tracked as per %s\n", FDSfile);
fprintf(stderr,"center_x=%f,center_y=%f\n",center_x,center_y);
          }
        else
          sphere2img(crlt_rad, crln_rad, crlt_obs_rad, crln_obs_rad, &center_x, &center_y, x0, y0, rsunpix, pa_rad, 0, 0, 0, 0);
        }
    // center_x and center_y are target locations for the center of the desired patch, pixels from 0 of as-is image.
    // note HMI has not been flipped at this point

    target_x = center_x;
    target_y = center_y;
    int x1Orig;
    int y1Orig;
    
    int x1 = round(target_x + llx);
    int y1 = round(target_y + lly);
    int x2 = x1 + pixwidth - 1;
    int y2 = y1 + pixheight - 1;
    if (do_register)
      {
      x1 -= register_padding;
      y1 -= register_padding;
      x2 += register_padding;
      y2 += register_padding;
      }
    crpix1 = 1 + x0 - x1;
    crpix2 = 1 + y0 - y1;
fprintf(stderr,"at box define, crpix1=%f, x0=%f, x1=%d\n",crpix1,x0,x1);

    int start1[2] = {x1, y1};
    int end1[2] = {x2, y2};
    
    outRS = drms_create_records(drms_env, 1, outseries, DRMS_PERMANENT, &status);
    if (status) {fprintf(stderr,"Output series is %s, ",outseries); DIE("Cant make outout record");}
    outRec = outRS->records[0];
      
    
    if (do_crop)
    {
        /* BEFORE SEG LOOP */
        int retStat;
        dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &retStat);
        rsun_rad = asin(rsun_ref/dsun_obs); 
        rsun = rsun_rad*Rad2arcsec/cdelt;
        r2 = rsun*rsun;
        r2 *= 0.9995;
        
        /* WE HAVE TO CALCULATE this_x0 and this_y0 now, before paIs180 modifies crpix1 and crpix2. */
        this_x0 = crpix1 - 1;
        this_y0 = crpix2 - 1;
    }
      
    /* We need to adjust pa and crota if |pa| == 180. Originally, this was done in code that was mixed with 
     * per-segment code. But now that we made a per-segment loop, we do not want to execute this code more than 
     * once. So I moved it here, right before the segment loop. */

    /* These might be modified in the paIs180 block of code, but we need to use the original values 
     * later. */     
    x1Orig = x1;
    y1Orig = y1;

    if (paIs180)
    {
        pa -= 180.0; /* This actually does not get used again during the rest of the record loop, and is not used by the segment loop. */
        crota -= 180.0;

        crpix1 = 1 + x2 - x0;
        crpix2 = 1 + y2 - y0;
fprintf(stderr,"after rotate, crpix1=%f\n",crpix1);

        // adjust internal quantities for after the flip, may be needed by do_register. 
        /* But we will need the original values of x1 and y1 in the code that handles the patch-image intersection.
         * These quantities are saved in x1orig and y1orig.
         */
        x1 = inAxis[0] - 1 - x2;
        y1 = inAxis[1] - 1 - y2;
        target_x = inAxis[0] - 1 - center_x;
        target_y = inAxis[1] - 1 - center_y;
        // cosa = 1.0; sina = 0.0;
        strcat(history, "Image rotated 180 degrees.");
fprintf(stderr,"after flip x1=%d, target_x=%lf, y1=%d, target_y=%lf\n",x1,target_x,y1,target_y);
    }
    
    /* At this point, all the modifications to x1, y1, target_x, and target_y have been performed, 
     * so it is safe to calculate dx and dy. */
     
    /*
     * To simplify things, error-out if the input segment's patches dimensions are not equivalent.
     * Also, we need to calculate dx and dy, which require that we know the intersection of the
     * patch box with the image.
     */
    DRMS_Segment_t *firstSeg = NULL;
    HIterator_t *segIter = NULL;
    DRMS_Segment_t *orig = NULL;

    float regShiftX = 0;
    float regShiftY = 0;
    char patchIntersection;
    
    while ((inSeg = drms_record_nextseg2(inRec, &segIter, 1, &orig)) != NULL)
    {
        if (!firstSeg)
        {
            firstSeg = inSeg;
            
            /* Calculate the ultimate patch dimensions. */
            /* We need to use the original x1 and y1, not the ones modified by paIs180. */
            if (x1Orig >= inSeg->axis[0] || y1Orig >= inSeg->axis[1] || x2 < 0 || y2 < 0)
            {  
                // patch completely outside image
                patchIntersection = 'N';
                /* The record needs to be skipped. Downstream code will handle this. */
            }
            else if (x1Orig >= 0 && y1Orig >= 0 && x2 < inSeg->axis[0] && y2 < inSeg->axis[1])
            {  
                // patch entirely in image.
                patchIntersection = 'F';
            }
            else
            { 
                // patch partly outside image.
                patchIntersection = 'P';
            }            
        }
        else
        {
            if (!hasSegList && !doingAllSegs)
            {
                /* No need to ensure post-first segment matches the first segment since we are processing only the first segment. */
                break;
            }

            /* Ensure that both segments are of the same dimensions (both in number and size of each dimension). */
            int iSeg;
            int mismatch;
            
            mismatch = 0;
            if (firstSeg->info->naxis == inSeg->info->naxis)
            {
                for (iSeg = 0; iSeg < firstSeg->info->naxis; iSeg++)
                {
                    if (firstSeg->axis[iSeg] != inSeg->axis[iSeg])
                    {
                       mismatch = 1;
                       break;
                    }
                }
            }
            else
            {
                mismatch = 1;
            }

            if (!mismatch)
            {
                if (firstSeg->info->protocol == DRMS_TAS && inSeg->info->protocol == DRMS_TAS)
                {
                    for (iSeg = 0; iSeg < firstSeg->info->protocol; iSeg++)
                    {
                        if (firstSeg->blocksize[iSeg] != inSeg->blocksize[iSeg])
                        {
                            mismatch = 1;
                            break;
                        }
                    }
                }
            }
            
            if (mismatch)
            {
                DIE("When processing more than one segment, all segments' dimensions must match.");
            }        
        }
    } /* Not the real segment loop. */
    
    if (segIter)
    {
        hiter_destroy(&segIter);
    }
     
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);
    drms_sprint_rec_query(inQuery, inRec);
    drms_setkey_string(outRec, "SOURCE", inQuery);
    drms_appendhistory(outRec, inQuery, 1);
    drms_setkey_time(outRec, "T_REC", trec);
    drms_setkey_double(outRec, "CRVAL1", 0.0);
    drms_setkey_double(outRec, "CRVAL2", 0.0);
    drms_setkey_double(outRec, "CRDELT1", cdelt);
    drms_setkey_double(outRec, "CRDELT2", cdelt);
    
    if (do_register)
    {
        /* It used to be the case that in the per-segment code, if do_register was set, then 
         * crota got set to 0.0, and they the crota2 keyword got assigned that value. 
         * But now that we have a segment loop, we cannot set crota to 0.0 in the per-segment
         * code, otherwise we lose the value of crota needed by the second and greater 
         * segment. */
        drms_setkey_double(outRec, "CROTA2", 0.0);
    }
    else
    {
        drms_setkey_double(outRec, "CROTA2", crota);
    }
        
    drms_setkey_string(outRec, "CONTENT", "Tracked Extracted Patches, made by im_patch");
    drms_appendcomment(outRec, "Patches", 1);
    drms_setkey_string(outRec, "RequestID", requestid);
    drms_setkey_string(outRec, "HGBOXUNITS", boxunits);
    drms_setkey_string(outRec, "HGLOCUNITS", locunits);
    drms_setkey_float(outRec, "HGWIDE", width);
    drms_setkey_float(outRec, "HGHIGH", height);
    drms_setkey_float(outRec, "HGCRLN", crln);
    drms_setkey_float(outRec, "HGCRLT", crlt);
    drms_setkey_int(outRec, "HGCARROT", car_rot);
    drms_setkey_time(outRec, "HGTSTART", t_start);
    drms_setkey_time(outRec, "HGTSTOP", t_stop);
    drms_setkey_time(outRec, "DATE", time(0) + UNIX_EPOCH);
    drms_setkey_string(outRec, "HGQUERY", in);
fprintf(stderr,"DATAMIN=%f, DATAMAX=%f\n",drms_getkey_float(outRec,"DATAMIN",NULL),drms_getkey_float(outRec,"DATAMAX",NULL));

/*
 *               writing the extracted region data file
*/
    /* Make sure the output series is compatible with the input series. 
     * The input series and output series must have the same-named segments in
     * the same order, with one exception. If there is only a single segment in 
     * the input series and only a single segment in the output segment, then 
     * the names can differ. In that case, this code will assume that the output
     * patch will go in the single output series.
     */
    int numInSegs = hcon_size(&(inTemplate->segments));
    int numOutSegs = hcon_size(&(outTemplate->segments));
    
    if (numInSegs != numOutSegs)
    {
        DIE("Input and output series are incompatible (the number of segments differs).\n");
    }
    
    if (numInSegs < 1)
    {
        DIE("The input series must have at least one segment.\n");
    }
    
    /* inSeg - 
     * The third argument to drms_record_nextseg2() indicates that a link should be followed.
     * 
     */
    int iSeg = 0;
    char msgBuf[512];

    while ((inSeg = drms_record_nextseg2(inRec, &segIter, 1, &orig)) != NULL)
    {
        /* THIS IS THE REAL SEGMENT LOOP. */
        if (!hasSegList && !doingAllSegs && iSeg > 0)
        {
            /* We've already processed the first segment, and there is no seglist specifier, so we are not processing segments other 
             * than the first one. */
             break;
        }
        
        /* Pin the first segment that meets these conditions:
         * + The segment protocol is either FITS or TAS.
         * + NAXIS == 2.
         * + CTYPE1 == CTYPE2 == 'HPLN-TAN' (but the keyword-reading code
         *   already rejects records for which this isn't true).
         * + The existence of these keywords for the segment:
         *   CAR_ROT
         *   CRLT_OBS
         *   CRLN_OBS
         *   CRPIX1
         *   CRPIX2
         *   CRVAL1
         *   CRVAL2
         *   CROTA2
         *   CTYPE1
         *   CTYPE2
         *   DSUN_OBS
         *   CDELT1
         *   T_REC or T_OBS
         *   (but the keyword-reading code already rejects records for which 
         *   this isn't true).
         * 
         * Once a segment is pinned, then all following segments must have WCS
         * keyword values that match the WCS keyword values of the pinned segment.
         * However, since this module rejects series with per-segment WCS keywords, 
         * this will also already be true. 
         * 
         * So as long as the current segment has NAXIS == 2 and is a FITS or TAS segment,
         * then we can process it, even if it is a linked segment.
         */
        
        /* The input series and output series must have the same-named segments in
         * the same order, with one exception. If there is only a single segment in 
         * the input series and only a single segment in the output segment, then 
         * the names can differ. In that case, this code will assume that the output
         * patch will go in the single output series.
         *
         * A check for the two series having the same number of segments was already performed, 
         * before the segment loop.
         *
         * Must use orig, not inSeg, since inSeg could be a followed link.
         */
        outSeg = drms_segment_lookupnum(outRec, orig->info->segnum);
        
        if (outSeg == NULL || (numInSegs != 1 && strncasecmp(orig->info->name, outSeg->info->name, DRMS_MAXSEGNAMELEN - 1) != 0))
        {
            /* Incompatible output series. */    
            snprintf(msgBuf, sizeof(msgBuf), "Input and output series are incompatible (there is no segment named %s in the output series).\n", orig->info->name);
            DIE(msgBuf);
        }
        
        if (inSeg->info->naxis != 2 || (inSeg->info->protocol != DRMS_FITZ && inSeg->info->protocol != DRMS_FITS && inSeg->info->protocol != DRMS_TAS))
        {
            /* The input series' segment must be a 2D FITS-type of segment. */
            snprintf(msgBuf, sizeof(msgBuf), "Input segment %s is not a 2D FITS segment.\n", orig->info->name);
        }
        
        inAxis[0] = inSeg->axis[0];
        inAxis[1] = inSeg->axis[1];

        /* We need to use the original x1 and y1, not the ones modified in paIs180. */
        if (patchIntersection == 'N')
        {  // patch completely outside image
fprintf(stderr, "patch completely outside image, x1=%d, y1=%d, x2=%d, y2=%d\n", x1Orig, y1Orig, x2, y2); 
            /* If the patch is completely outside of one image, it is completely outside all images
             * since all images must be of the same dimensions. So break out of segment loop. */
            nextRec = 1;
            break;
        }
        else if (patchIntersection == 'F')
        {  // patch entirely in image.
            status = 0;
            outArray = drms_segment_readslice(inSeg, DRMS_TYPE_FLOAT, start1, end1, &status);
            if (status) DIE("Cant read input record");
fprintf(stderr,"$$$$$$$ outArray bzero, bscale are %f, %f\n", outArray->bzero, outArray->bscale);
        }
        else
        { // patch partly outside image.
            int dims[2] = {x2-x1Orig+1,y2-y1Orig+1};
            int start2[2], end2[2];
            
            int x,y;
            outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, dims, NULL, &status);
            drms_array2missing(outArray);
            
            start2[0] = x1Orig < 0 ? 0 : x1Orig;
            start2[1] = y1Orig < 0 ? 0 : y1Orig;
            end2[0] = x2 >= inAxis[0] ? inAxis[0]-1 : x2;
            end2[1] = y2 >= inAxis[1] ? inAxis[1]-1 : y2;
            
            status = 0;
            inArray = drms_segment_readslice(inSeg, DRMS_TYPE_FLOAT, start2, end2, &status);        
            if (status || !inArray) DIE("Cant read input record");
        
//fprintf(stderr,"$$$$$$$ inArray bzero, bscale are %f, %f\n",inArray->bzero, inArray->bscale);
            outArray->bzero = inArray->bzero;
            outArray->bscale = inArray->bscale;
            
            int start3[2] = {start2[0]-start1[0],start2[1]-start1[1]}; // fetched slice offset in outarray
            int end3[2] = {end2[0]-end1[0],end2[1]-end1[1]};           // fetched slice offset in outarray
            int n3x = end2[0] - start2[0] + 1;
            int n3y = end2[1] - start2[1] + 1;
        
            for (y=0; y< n3y; y++)
            {
                for (x=0; x < n3x; x++)
                {
                    *((float *)outArray->data + ((start3[1]+y)*dims[0] + x + start3[0])) =
                        *((float *)inArray->data + (y*n3x + x));
                }
            }
        
            drms_free_array(inArray);
        }

        use_bzero = outArray->bzero;
        use_bscale = outArray->bscale;
    
        // Handle crop if wanted
        if (do_crop)
        {
            int stat;
            // rsun_ref and cdelt defined above            
            float *data = (float *)outArray->data;
            int i,j;
            int nx = outArray->axis[0];        
            int ny = outArray->axis[1];
        
            for (j=0; j<ny; j++)
            {
                for (i=0; i<nx; i++)
                {
                    double x = i - this_x0;
                    double y = j - this_y0;
                
                    if ((x*x + y*y) >= r2)
                    {
                        data[j*nx + i] = DRMS_MISSING_FLOAT;
                    }
                }
            }
        }
        
        // Always handle special case where |pa| == 180
        if (paIs180)
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
            {
                for (i=0; i<nx; i++)
                {
                    val = data[j*nx + i];
                    data[j*nx + i] = data[(ny - 1 - j)*nx + nx - 1 - i];
                    data[(ny - 1 - j)*nx + nx - 1 - i] = val;
                }
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
        }
    
        /*
         *  Register fraction of pixel to target location, if desired
         */
        if (do_register)
        {
            /* IN SEG LOOP */
            DRMS_Array_t *tmpArray;
            int status;
            int i,j;
            float *data = (float *)outArray->data;
            float *newdata = NULL;
            int newnx, newny;
            int nx = outArray->axis[0];
            int ny = outArray->axis[1];
            float midx=(nx-1)/2.0;
            float midy = (ny-1)/2.0;
            float dx = (x1 + midx) - target_x;
            float dy = (y1 + midy) - target_y;

            /* We do not want to adjust crpix1 and crpix2 more than once per record. */
            if (iSeg == 0)
            {
                sprintf(history+strlen(history), "\nImage registered by shift of (%0.3f,%0.3f) pixels.", dx, dy);
                crpix1 += dx - register_padding;
                crpix2 += dy - register_padding;
            }
        
            regShiftX = dx;
            regShiftY = dy;
        
            int dtyp = 3;
        
            if (status=image_magrotate((void *)data, nx, ny, dtyp, crota, 1.0, regShiftX, regShiftY, &(void *)newdata, &newnx, &newny, 1, 0))
            {
                DIE("XXXX Failure in call to image_magrotate.  Either malloc error or nx or ny too large. Must quit.\n");
            }
        
            if (newnx != nx || newny != ny)
            {
                fprintf(stderr,"image_magrotate changed dimensions: nx want %d got %d, ny want %d got %d\n",nx,newnx,ny,newny);
            }
        
            drms_free_array(outArray);
            int outdims[2] = {pixwidth, pixheight};
            outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outdims, NULL, &status);
            data = (float *)outArray->data;
            for (j=0; j<pixheight; j++)
            {
                for (i=0; i< pixwidth; i++)
                {
                    data[j*pixwidth + i] = newdata[(j+register_padding)*nx + i+register_padding];
                }
            }

            /* This was being leaked - over 4MB per call. */
            if (newdata)
            {
                free(newdata);
                newdata = NULL;
            }

            outArray->bzero = use_bzero;
            outArray->bscale = use_bscale;
        }
    
        if (wantFAKE && FDSfile)
        {
            double r;
            int ix0, iy0;
            int nx = outArray->axis[0];
            int ny = outArray->axis[1];
            double x0, y0, r2;
            float *data = (float *)outArray->data;
            float datamax = drms_getkey_float(outRec, "DATAMAX", NULL);
            int ix,iy,ir,ir2;
        
            if (NoTrack)
            {
                x0 = trackx/cdelt + crpix1;
                y0 = tracky/cdelt + crpix2;
            }
            else
            {
                x0 = nx/2.0;
                y0 = ny/2.0;
            }
    
            ix0 = x0+0.5;
            iy0 = y0+0.5;
            // r = 28.9/cdelt;  // venus radius in pixels
            r = track_radius/cdelt;  // venus radius in pixels
            ir = r + 0.5;
            r2 = r*r;
        
            for (ix = -ir; ix <= ir; ix++)
            {
                for (iy = -ir; iy <= ir; iy++)
                {
                    x = x0 - ix0 + ix;
                    y = y0 - iy0 + iy;
                    if (x*x + y*y < r2 && ix0+ix >=0 && ix0+ix < nx && iy0+iy >=0 && iy0+iy < ny)
                    {
                        if (wantFAKEwhite)
                        {
                            data[(iy+iy0)*nx + (ix+ix0)] = datamax;
                        }
                        else
                        {
                            if (ix==0 || iy==0)
                            {
                                data[(iy+iy0)*nx + (ix+ix0)] = datamax;
                            }
                            else
                            {
                                data[(iy+iy0)*nx + (ix+ix0)] = 0.0;
                            }
                         }
                    }
                }
            }
        }

        if (iSeg == 0)
            {
            drms_setkey_double(outRec, "CRPIX1", crpix1);
            drms_setkey_double(outRec, "CRPIX2", crpix2);
            }
    
        set_statistics(outSeg, outArray, 1);
        
        if (export_keys)
          {
          status = drms_segment_writewithkeys(outSeg, outArray, 0);
          }
        else
          {
          status = drms_segment_write(outSeg, outArray, 0);
          }
        
        if (status)
          {
          DIE("problem writing file");
          }
        
        if (outArray)
          {
          drms_free_array(outArray);
          outArray = NULL;
          }

        iSeg++;
    } /* end segment loop */

    if (segIter)
      {
      hiter_destroy(&segIter);
      }

    if (log)
      {
      drms_fprint_rec_query(log, outRec);
      fprintf(log, "\n");
      }
      
    if (*history)
      drms_appendhistory(outRec, history, 1);
      
    if (nextRec)
      {    
      drms_close_records(outRS, DRMS_FREE_RECORD);
      continue;
      }
    else
      {      
      drms_close_records(outRS, DRMS_INSERT_RECORD);
      drms_free_array(outArray);
      }
    
    OK_recs += 1;
    } /* end record loop */

  drms_close_records(inRS, DRMS_FREE_RECORD); 
  if (in_filename)
    {
    unlink(in_filename);
    }

  if (OK_recs)
    printf("  DONE, %d records made! \n \n",OK_recs);
  else
    printf("  DONE but no useful records created! \n \n");

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

// In cases known to not have compact slotted series and cadence is specified
// generate explicit recordset list of closest good record to desired grid
// First get vector of times and quality
// Then if vector is not OK, quit.
// then: make temp file to hold recordset list
//       start with first time to define desired grid,
//       make array of desired times.
//       make empty array of recnums
//       search vector for good images nearest desired times
//       for each found time, write record query


#define DIE_get_recset(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return NULL;}
char *get_input_recset(DRMS_Env_t *drms_env, char *inQuery)
  {
  static char newInQuery[DRMS_MAXSERIESNAMELEN+2];
  int epoch_given = cmdparams_exists(&cmdparams, "epoch");
  TIME epoch, t_epoch;
  DRMS_Array_t *data;
  DRMS_Record_t *inTemplate;
  TIME t_start, t_stop, t_now, t_want, t_diff, this_t_diff;
  int status = 1;
  int nrecs, irec;
  int nslots, islot;
  long long *recnums;
  TIME *t_this, half;
  TIME cadence;
  double *drecnum, *dquality;
  int quality;
  long long recnum;
  char keylist[DRMS_MAXQUERYLEN];
  char filename[DRMS_MAXSERIESNAMELEN];
  char *tmpdir;
  FILE *tmpfile;
  char newIn[DRMS_MAXQUERYLEN];
  char seriesname[DRMS_MAXQUERYLEN];
  char *lbracket;
  char *at = index(inQuery, '@');
  int npkeys;
  char *timekeyname;
  double t_step;

  strcpy(seriesname, inQuery);
  lbracket = index(seriesname,'[');
  if (lbracket) *lbracket = '\0';
  inTemplate = drms_template_record(drms_env, seriesname, &status);
  if (status || !inTemplate) DIE_get_recset("Input series can not be found");

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
        if(i > 0) DIE_get_recset("Input series must have TIME keyword first, for now...");
        timekeyname = pkey->info->name;
        t_step = drms_keyword_getdouble(drms_keyword_stepfromslot(pkey), &status);
        if (status) DIE_get_recset("problem getting t_step");
        t_epoch = drms_keyword_getdouble(drms_keyword_epochfromslot(pkey), &status);
        if (status) DIE_get_recset("problem getting t_epoch");
        }
    }
  else
    DIE_get_recset("Must have time prime key");
  epoch = epoch_given ? params_get_time(&cmdparams, "epoch") : t_epoch;

  if (at && *at && ((strncmp(inQuery,"aia.lev1[", 9)==0 ||
                    strncmp(inQuery,"hmi.lev1[", 9)==0 ||
                    strncmp(inQuery,"aia.lev1_nrt2[",14)==0 ||
                    strncmp(inQuery,"hmi.lev1_nrt[", 13)==0 ) ||
                   epoch_given))
    {
    char *ip=(char *)inQuery, *op=newIn, *p;
    long n, mul;
    while ( *ip && ip<at )
      *op++ = *ip++;
    ip++; // skip the '@'
    n = strtol(ip, &p, 10); // get digits only
    if (*p == 's') mul = 1;
    else if (*p == 'm') mul = 60;
    else if (*p == 'h') mul = 3600;
    else if (*p == 'd') mul = 86400;
    else 
      DIE_get_recset("cant make sense of @xx cadence spec");
    cadence = n * mul;
    ip = ++p;  // skip cadence multiplier
    while ( *ip )
      *op++ = *ip++;
    *op = '\0';
    half = cadence/2.0;
    sprintf(keylist, "%s,QUALITY,recnum", timekeyname);
    data = drms_record_getvector(drms_env, newIn, keylist, DRMS_TYPE_DOUBLE, 0, &status);
    if (!data || status)
      {
      fprintf(stderr, "status=%d\n", status);
      DIE_get_recset("getkey_vector failed\n");
      }
    nrecs = data->axis[1];
    irec = 0;
    t_this = (TIME *)data->data;
    dquality = (double *)data->data + 1*nrecs;
    drecnum = (double *)data->data + 2*nrecs;
    if (epoch_given)
      {
      int s0 = (t_this[0] - epoch)/cadence;
      TIME t0 = s0*cadence + epoch;
      t_start = (t0 < t_this[0] ? t0 + cadence : t0);
      }
    else
      t_start = t_this[0];
    t_stop = t_this[nrecs-1];
    nslots = (t_stop - t_start + cadence/2)/cadence;
    recnums = (long long *)malloc(nslots*sizeof(long long));
    for (islot=0; islot<nslots; islot++)
      recnums[islot] = 0;
    islot = 0;
    t_want = t_start;
    t_diff = 1.0e9;
    for (irec = 0; irec<nrecs; irec++)
        {
        t_now = t_this[irec];
        quality = (int)dquality[irec] & 0xFFFFFFFF;
        recnum = (long long)drecnum[irec];
        this_t_diff = fabs(t_now - t_want);
        if (quality < 0)
          continue;
        if (t_now <= (t_want-half))
          continue;
        while (t_now > (t_want+half))
          {
          islot++;
          if (islot >= nslots)
             break;
          t_want = t_start + cadence * islot;
          this_t_diff = fabs(t_now - t_want);
          t_diff = 1.0e8;
          }
        if (islot < nslots && this_t_diff <= t_diff)
          recnums[islot] = recnum;
        t_diff = fabs(t_now - t_want);
        }
    if (islot+1 < nslots)
      nslots = islot+1;  // take what we got.
    tmpdir = getenv("TMPDIR");
    if (!tmpdir) tmpdir = "/tmp";
    sprintf(filename, "%s/%sXXXXXX", tmpdir, module_name);
    mkstemp(filename);
    tmpfile = fopen(filename,"w");
    for (islot=0; islot<nslots; islot++)
      if (recnums[islot])
        fprintf(tmpfile, "%s[:#%lld]\n", seriesname, recnums[islot]);
    fclose(tmpfile);
    free(recnums);
    drms_free_array(data);
    sprintf(newInQuery,"@%s", filename);
    return(newInQuery);
    }
  else
	  return(inQuery);
  }

// FDS file is in format of SDO FDS "transit" files.
// returns xp, yp, and radius of target all in arcsec
int get_tracking_xy(char *FDSfile, TIME want, double *xp, double *yp, double *radius)
  {
  int iline;
  char buf[200];
  TIME now0;
  FILE *data = fopen(FDSfile,"r");
  double lambda, phi, x, y, dist, obj_radius;
  double frac, x0, y0, obj_radius0;
  char object[100], now_txt[100];
  if (!data)
    return(1);

  // skip first 13 lines, for lunar transit data
  for (iline=0; iline<13; iline++)
    fgets(buf, 200, data);

  // read info until end, stop at first time after want
  while (fscanf(data, "%s %s %lf %lf %lf %lf %lf %lf", object, now_txt, &lambda, &phi, &x, &y, &dist, &obj_radius) == 8)
    {
    int yyyy, doy, hh, mm, ss;
    int m,d;
    int dim[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    double distance;
    TIME now;
    // fix doordinate system to match HMI view
// not for ISON
    x = -x;
    y = -y;
    obj_radius *= 3600;
    // extract time
    sscanf(now_txt,"%4d%3d.%2d%2d%2d", &yyyy, &doy, &hh, &mm, &ss);
    if (yyyy % 4 == 0) dim[1] = 29; else dim[1] = 28;
    for (m=1; m<=12; m++)
      {
      if (doy > dim[m-1])
        doy -= dim[m-1];
      else
        break;
      }
    d = doy;

    sprintf(now_txt, "%4d.%02d.%02d_%02d:%02d:%02d_UTC", yyyy, m, d, hh, mm, ss);
    now = sscan_time(now_txt);
    // report closest x,y for time want
    if (now < want)
      {
      now0 = now;
      x0 = x;
      y0 = y;
      obj_radius0 = obj_radius;
      }
    else
      {
      if (fabs(now-want) >= fabs(now0-want))
        {
        obj_radius = obj_radius0;
        if (now == want)
          {
          x = x0;
          y = y0;
          }
        else
          {
          frac = (want-now0)/(now-now0); 
          x = x0 + (x-x0)*frac;
          y = y0 + (y-y0)*frac;
          }
        }
      break;
      }
    }
  fclose(data);
  fprintf(stderr,"Found %s, x=%f, y=%f, radius=%f (all arcsecs)\n",now_txt,x,y,obj_radius);
  *xp = x;
  *yp = y;
  *radius = obj_radius;
  return(0);
  }
