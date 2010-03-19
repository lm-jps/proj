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
   in arcsec.

      Written by: Bala Poduval as extract_region, modified by Phil Scherrer

   @par Synopsis:

   @code
   hg_patch in=<recordset> out=<series_out>  \
   log=<logfile> \
   t_start=<start time> \
   t_stop=<end time> \
   car_rot=<Carrington Rotation> \
   crln=<longitude> \
   crlt=<latitude> \
   height=<height>
   width=<width>
   t_ref=<ref_time> \
   aswide=<arcsec_width> \
   ashigh=<arcsec_height> \
   asx=<arcsec_center_x> \
   asy=<arcsec_center_y> \
   requestid=<RequestID>

   @endcode

   The extracted region may be defined in one of two ways, either by its disk-center longitude and latitude
   in a reference Carrington rotation or by its disk location at a reference image time, ref_time.
   If the car_rot parameter is specified, the Carrington location method is used.
   If the t_ref parameter is present, the box location and size are specified in arcsec.

   If the logfile is present, a RecordSet query will be written to it for each image created  containing a query
   that will returnthat image.

   RequestID, if present, is the ID of an on-requst procssing export request.

   If t_start or t_stop are not present, those values will be inferred from
   the on-disk time span for the center of the patch -90 degrees to +90 degrees from CM.
   This module computes the beginning and ending time for tracking the 
   region identified at time = t_ref, asx and asy or by the Carrington coordinates of the box center.
   In the Carrington case the module extracts a rectangular region 
   with user defined size. the height and width of the box are specified in degrees
   of latitude and longitude with the pixel size of the box is computed from the
   projection of the box at CM.  The pixel size of the rectangular box remains constant 
   during tracking.  The box height defaults to the width.  The width defaults to 10 degrees.
  
   If the input recordset spec is exactly "[$]" it will be discarded and the time limits will be taken from t_start, t_stop, car_rot, and/or t_ref as appropriate.  This is to allow exports via jsoc_fetch to not require a bounding recordset to be specified.

   If t_start and/or t_stop are specified the referecne time or disk center location might be not
   included.
   
   Option: the width, height, latitude, and longitude are rounded to the nearest 5 degrees.
   This will allow more frequent re-use of tracked boxes.

   If the input is only a seriesname, it must have a first prime key of type time.

   If the t_ref parameter is specified, it must refer to a non-missing image in the input series.

   The output seriesname defaults to the input seriesname with a suffix of "_hgpatch".

   At present, it doesn't have any flags.

   @par Flags:
   @c none

   @par GEN_FLAGS:
   Ubiquitous flags present in every module.
   @ref jsoc_main

   @param in The input data series.
   @param out The output series.

   @par Exit_Status:
   Brief description of abnormal, non-zero, exit values.
      Arguemnts are to be in the prescribed units and formats. Otherwise,
      the program stops, with a segmentation fault or similar error message.

   @par Example:
   Extracts a rectangular region of size 20 x 10 degrees and tracks it from 
   when the region is just appearing at the E-limb to t_stop.

   @code
   hg_patch in='mdi.fd_M_96m_lev18' out='su_bala.extractreg' car_rot=2009 crln=295 crlt=5 width=30 height=20

   @endcode

   @par Example:

   @code

   @endcode

   @bug
   No doubt some bugs.

Fails when SOHO upside down
needs more prime keys in output series
fails for CR mode, too high by maybe factor of two.
in CR mode, box dims not centered.  need to find image at disk center - or at least calc for
that image...


*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "jsoc_main.h"

void HeliographicLocation(TIME t, int *crot, double *L, double *B);
TIME HelioographicTime(int crot, double L);

#include "astro.h"

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
    {ARG_FLOAT, "crln", "-1", "Carrington longitude of the box center when it crosses CM"},
    {ARG_FLOAT, "crlt", "-999", "Carrington latitude of the box center when it crosses CM"},
    {ARG_FLOAT, "width", "0", "width of box in degrees of longitude"},
    {ARG_FLOAT, "height", "0", "height of box in degrees of latitude when it corsses CM"},
    {ARG_TIME, "t_start", "JD_0", "Start time, defaults to time at 90E"},
    {ARG_TIME, "t_stop", "JD_0", "End time, defauolts to 90W"},
    {ARG_TIME, "t_ref", "JD_0", "Time for which asx and asy apply, implies ref image."},
    {ARG_FLOAT, "aswide", "0", "Extract box width in arcsec"},
    {ARG_FLOAT, "ashigh", "0", "Extract box height in arcsec"},
    {ARG_FLOAT, "asx", "0", "Location of extract box center at t_ref"},
    {ARG_FLOAT, "asy", "0", "Location of extract box center at t_ref"},
    {ARG_END}
  };

int img2sphere (double x, double y, double ang_r, double latc, double lonc,
        double pa, double *rho, double *lat, double *lon, double *sinlat,
        double *coslat, double *sig, double *mu, double *chi);
int sphere2img (double lat, double lon, double latc, double lonc,
        double *x, double *y, double xcenter, double ycenter,
        double rsun, double peff, double ecc, double chi,
        int xinvrt, int yinvrt);

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)

int DoIt(void)
{
  CmdParams_t *params = &cmdparams;
  DRMS_RecordSet_t *inRS, *outRS;
  DRMS_Record_t *inRec, *outRec, *inTemplate, *outTemplate;
  DRMS_Segment_t *inSeg, *outSeg;
  DRMS_Array_t *inArray, *outArray;
  int i, ii, status = DRMS_SUCCESS, nrecs; 
  int  irec;
  double center_x, center_y, crpix1, crpix2; 
  double rsun_ref, dsun_obs, rsun, rsunpix;
  double crln_obs, crlt_obs;
  double cdelt;
  double urx, ury, llx, lly;
  int inAxis[2];
  int this_car_rot;
  char *ctype1, *ctype2;
  TIME t_rec;
  double box_x, box_y;
  char outseries[DRMS_MAXNAMELEN];
  char inseries[DRMS_MAXNAMELEN];
  char in[DRMS_MAXQUERYLEN];

  char *ingiven = (char *)params_get_str(params, "in");
  char *inparam;
  char *lbracket;
  char *outparam = (char *)params_get_str(params, "out");
  char *logfile = (char *)params_get_str(params, "log");
  char *requestid = (char *)params_get_str(params, "requestid");
  TIME t_start = params_get_time(params, "t_start");
  TIME t_stop = params_get_time(params, "t_stop");
  TIME t_ref = params_get_time(params, "t_ref");
  double width = params_get_double(params, "width");
  double height = params_get_double(params, "height");
  int car_rot = params_get_int(params, "car_rot");
  double crln = params_get_double(params, "crln");
  double crlt = params_get_double(params, "crlt");
  double aswide = params_get_double(params, "aswide");
  double ashigh = params_get_double(params, "ashigh");
  double asx = params_get_double(params, "asx");
  double asy = params_get_double(params, "asy");
  TIME tNotSpecified = sscan_time("JD_0");
  int do_arcsec = 0;  // do-arc-sec mode.
  double crvalx = 0.0;
  double crvaly = 0.0;
  double crota, sina, cosa;
  double pa, deltlong;
  int firstimage = 1;

  FILE *log = NULL;

  if (strcmp(logfile, "NOTSPECIFIED") != 0)
    {
    log = fopen(logfile, "w");
    if (!log) DIE2("Can not create log file.",logfile);
    }

  inparam = strdup(ingiven);
  lbracket = index(inparam, '[');
  // first, get input series names.
  if (lbracket)
    {
    int n = lbracket - inparam;
    strncpy(inseries, inparam, n);
    inseries[n] = '\0';
    if (strcmp(lbracket, "[$]") == 0) // Special case, discard the explicit last-record spec to allow exports via jsoc_fetch.
      {
      *lbracket = '\0';
      lbracket = NULL;
      }
    }
  else
    strcpy(inseries, inparam);

  if (car_rot > 0 && t_ref == tNotSpecified) do_arcsec = 0;
  else if (car_rot <= 0 && t_ref != tNotSpecified) do_arcsec = 1;
  else DIE("Must provide either car_rot or t_ref, but not both.");

  if (do_arcsec) // Arc-sec specification, get ref image information
    { // the image for ref_time is supposed to be present in the series.
    char t_ref_text[100];
    sprint_at(t_ref_text, t_ref);
    sprintf(in, "%s[%s]", inseries, t_ref_text);
    inRS = drms_open_records(drms_env, in, &status); if (status || inRS->n == 0) DIE2("No input data found for t_ref",in);
    inRec = inRS->records[0];;
    this_car_rot = drms_getkey_int(inRec, "CAR_ROT", &status); TEST_PARAM("CAR_ROT");
    crlt_obs = drms_getkey_double(inRec, "CRLT_OBS", &status); TEST_PARAM("CRLT_OBS");
    crln_obs = drms_getkey_double(inRec, "CRLN_OBS", &status); TEST_PARAM("CRLN_OBS");
    crpix1 = drms_getkey_double(inRec, "CRPIX1", &status)-1; TEST_PARAM("CRPIX1");
    crpix2 = drms_getkey_double(inRec, "CRPIX2", &status)-1; TEST_PARAM("CRPIX2");
    crvalx = drms_getkey_double(inRec, "CRVAL1", &status); TEST_PARAM("CRVAL1");
    crvaly = drms_getkey_double(inRec, "CRVAL2", &status); TEST_PARAM("CRVAL2");
    crota = drms_getkey_double(inRec, "CROTA2", &status); TEST_PARAM("CROTA2");
    pa = -crota;
    sina = sin(crota*Deg2Rad);
    cosa = cos(crota*Deg2Rad);
    ctype1 = strdup(drms_getkey_string(inRec, "CTYPE1", &status)); TEST_PARAM("CTYPE1");
      if (strcmp(ctype1, "HPLN-TAN") != 0) DIE2("CTYPE1 not HPLN-TAN as required, is: ", ctype1);
    ctype2 = strdup(drms_getkey_string(inRec, "CTYPE2", &status)); TEST_PARAM("CTYPE2");
      if (strcmp(ctype2, "HPLT-TAN") != 0) DIE2("CTYPE2 not HPLT-TAN as required, is: ", ctype2);
    rsun_ref = drms_getkey_double(inRec, "DSUN_REF", &status);
      if (status) rsun_ref = 6.96e8;
    dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status); TEST_PARAM("DSUN_OBS");
    rsun = asin(rsun_ref/dsun_obs)*Rad2arcsec; 
    cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
    rsunpix = rsun/cdelt;
    urx = aswide/(2.0*cdelt); llx = -urx; ury = ashigh/(2.0*cdelt); lly = -ury;
    urx = round(urx); ury = round(ury); llx = round(llx); lly = round(lly);
    center_x = PIX_X(asx,asy);
    center_y = PIX_Y(asx,asy);
    inSeg = drms_segment_lookupnum(inRec, 0);
    inAxis[0] = inSeg->axis[0];
    inAxis[1] = inSeg->axis[1];
    drms_close_records(inRS, DRMS_FREE_RECORD);
    firstimage=0;

    if(img2sphere (asx/rsun, asy/rsun, rsunpix, crlt_obs*Deg2Rad, crln_obs*Deg2Rad, pa*Deg2Rad,
      NULL, &crlt, &crln, NULL, NULL, NULL, NULL, NULL) < 0) DIE("Starting location is off the solar disk.");
    crln *= Rad2Deg;
    crlt *= Rad2Deg;
    deltlong = crln - crln_obs;
    car_rot = this_car_rot;
    if (deltlong >= 360.0)
      {
      car_rot--;
      crln -= 360.0;
      }
    else if (deltlong < 0)
      {
      car_rot++;
      crln += 360.0;
      }
    // now we have urx, ury, llx, lly, crln, crlt, car_rot for the box given in arcsec.
    }
  else // Carrington specification
    {
    // Examine parameters for required information.
    if (crln < 0) DIE("Box longitude must be specified.");
    if (crlt < -990) DIE("box latitude must be specified.");
    if (width <= 0) width = 10.0;
    if (height <= 0) height = width;
    }
// fprintf(stderr,"Carr box location, car_rot=%d, crln=%f, crlt=%f\n",car_rot, crln, crlt);

    // Get implied time limits from car_rot and crln
    // Use quickie estimate of carrtimes for now.
    // XXXXXX need to fix this later.
    t_ref = HeliographicTime(car_rot, crln);
    if (t_start == tNotSpecified)
      t_start = HeliographicTime(car_rot, crln + 90);
    if (t_stop == tNotSpecified)
      t_stop = HeliographicTime(car_rot, crln - 90);

  // get input seriesname and recordset, and output seriesname
  if (strcmp(inparam, "NOTSPECIFIED") == 0) DIE("Input series must be specified.");
  // first, get input and output series names.
  if (lbracket) // RecordSet query specified
    strncpy(in, inparam, DRMS_MAXQUERYLEN);
  else // need to generate query from limit times
    { 
    char t_start_text[100], t_stop_text[100];
    strncpy(inseries, inparam, DRMS_MAXNAMELEN);
    sprint_at(t_start_text, t_start);
    sprint_at(t_stop_text, t_stop);
    sprintf(in, "%s[%s-%s]", inseries, t_start_text, t_stop_text);
    }
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

// HERE before open records, get the step if needed, check key name for time, etc.
// for now blindly go on, but need to get some base numbers to convert box size to pixels.

  inRS = drms_open_records(drms_env, in, &status);
  if (status || inRS->n == 0)
           DIE("No input data found");
  nrecs = inRS->n;
// fprintf(stderr,"opened, found %d records\n", nrecs);

  // extract patches from each record
  for (irec = 0; irec < nrecs; irec++)
    {
    inRec = inRS->records[irec];
    if (status || !inRec) DIE("Record read failed.");
    TIME trec0 = drms_getkey_time(inRec, "T_REC", &status); TEST_PARAM("T_REC");
    if (trec0 < t_start)
      continue;
    if (trec0 > t_stop)
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
    crpix1 = drms_getkey_double(inRec, "CRPIX1", &status)-1; TEST_PARAM("CRPIX1");
    crpix2 = drms_getkey_double(inRec, "CRPIX2", &status)-1; TEST_PARAM("CRPIX2");
    pa = -drms_getkey_double(inRec, "CROTA2", &status); TEST_PARAM("CROTA2");

    // get coordinate mapping info using keywords from first record
    if (firstimage) // true for first image of Carrington type specs, not for do_arcsec specs.
      {
      firstimage = 0;
      ctype1 = strdup(drms_getkey_string(inRec, "CTYPE1", &status)); TEST_PARAM("CTYPE1");
      if (strcmp(ctype1, "HPLN-TAN") != 0) DIE2("CTYPE1 not HPLN-TAN as required, is: ", ctype1);
      ctype2 = strdup(drms_getkey_string(inRec, "CTYPE2", &status)); TEST_PARAM("CTYPE2");
      if (strcmp(ctype2, "HPLT-TAN") != 0) DIE2("CTYPE2 not HPLT-TAN as required, is: ", ctype2);
      rsun_ref = drms_getkey_double(inRec, "DSUN_REF", &status);
      if (status) rsun_ref = 6.96e8;
      dsun_obs = drms_getkey_double(inRec, "DSUN_OBS", &status); TEST_PARAM("DSUN_OBS");
      rsun = asin(rsun_ref/dsun_obs)*Rad2arcsec; 
      cdelt = drms_getkey_double(inRec, "CDELT1", &status); TEST_PARAM("CDELT1");
      // in principle use deltlong to get time to use for correcting crlt_obs, ignore for now
      // get ur, ll coords of box at CM.
      rsunpix = rsun/cdelt;

      sphere2img(Deg2Rad*(crlt+height/2), Deg2Rad*(width/2), Deg2Rad*crlt_obs, 0.0, &urx, &ury, 0.0, 0.0, rsunpix, pa*Deg2Rad, 0, 0, 0, 0);
      sphere2img(Deg2Rad*(crlt-height/2), -Deg2Rad*(width/2), Deg2Rad*crlt_obs, 0.0, &llx, &lly, 0.0, 0.0, rsunpix, pa*Deg2Rad, 0, 0, 0, 0);
      urx = round(urx); ury = round(ury); llx = round(llx); lly = round(lly);
      inSeg = drms_segment_lookupnum(inRec, 0);
      inAxis[0] = inSeg->axis[0];
      inAxis[1] = inSeg->axis[1];
      aswide = round((urx - llx)*cdelt);
      ashigh = round((ury - lly)*cdelt);
// fprintf(stderr,"box limits from (%f,%f) to (%f,%f)\n",llx,lly,urx,ury);
// fprintf(stderr,"Now Carr box location, car_rot=%d, crln=%f, crlt=%f, rsun(as)=%f, rsunpix=%f\n",car_rot, crln, crlt, rsun, rsunpix);
      }

    sphere2img(Deg2Rad*crlt, Deg2Rad*crln, Deg2Rad*crlt_obs, Deg2Rad*crln_obs, &center_x, &center_y, crpix1, crpix2, rsunpix, pa*Deg2Rad, 0, 0, 0, 0);
    int x1 = center_x + llx;
    int y1 = center_y + lly;
    int x2 = center_x + urx;
    int y2 = center_y + ury;

    int start1[2] = {x1, y1};
    int end1[2] = {x2, y2};

    if (x1 >= inAxis[0] || y1 >= inAxis[1] || x2 < 0 || y2 < 0)
{
// fprintf(stderr,"For %s, slice completely outside image\n",drms_getkey_string(inRec,"T_REC",NULL));
      continue;
}

    inSeg = drms_segment_lookupnum(inRec, 0);
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
      }
    crpix1 -= center_x;
    crpix2 -= center_y;
/*
 *
 *               writing the extracted region data file
 *
*/
    outRS = drms_create_records(drms_env, 1, outseries, DRMS_PERMANENT, &status);
    if (status) DIE("Cant make outout record");
    outRec = outRS->records[0];
    outSeg = drms_segment_lookupnum(outRec, 0);
    outArray->bzero = outSeg->bzero;
    outArray->bscale = outSeg->bscale;
    drms_setkey_float(outRec, "XCEN", WX((outArray->axis[0]+1)/2, (outArray->axis[1]+1)/2));
    drms_setkey_float(outRec, "YCEN", WY((outArray->axis[0]+1)/2, (outArray->axis[1]+1)/2));
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) DIE("problem writing file");
    drms_free_array(outArray);
    drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);
    drms_setkey_int(outRec, "CRPIX1", crpix1);
    drms_setkey_int(outRec, "CRPIX2", crpix2);
    drms_setkey_string(outRec, "CONTENT", "Tracked Extracted Patches, made by hg_patch");
    drms_setkey_string(outRec, "COMMENTS", "Patches");
    drms_setkey_string(outRec, "RequestID", requestid);
    drms_setkey_float(outRec, "HGASWIDE", aswide);
    drms_setkey_float(outRec, "HGASHIGH", ashigh);
    drms_setkey_float(outRec, "HGCRLN", crln);
    drms_setkey_float(outRec, "HGCRLT", crlt);
    drms_setkey_int(outRec, "HGCARROT", car_rot);
    drms_setkey_time(outRec, "HGTSTART", t_start);
    drms_setkey_time(outRec, "HGTSTOP", t_stop);
    drms_setkey_time(outRec, "DATE", time(0) + UNIX_EPOCH);

    if (log)
      {
      drms_fprint_rec_query(log, outRec);
      fprintf(log, "\n");
      }
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    } 

    drms_close_records(inRS, DRMS_FREE_RECORD); 

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
 *    x }           Plate locations, in units of the image radius, relative
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
 *      *x }        Plate locations, in units of the image radius, relative
 *      *y }          to the image center
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
