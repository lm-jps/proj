	    /*  The following code is adapted from SOI functions/keywords.c  */
#define PLATFORM_UNKNOWN	(0)
#define PLATFORM_UNRECOGNIZED	(1)
#define PLATFORM_SOHO		(2)
#define PLATFORM_GONGPLUS	(3)
#define PLATFORM_MWO60		(4)
#define PLATFORM_BBSO		(5)
#define PLATFORM_TRACE		(6)
#define PLATFORM_SPOLE_JSO	(10)
#define PLATFORM_GONG		(30)
#define PLATFORM_OBSPM		(40)

#define INSTRUMENT_UNKNOWN	(0)
#define INSTRUMENT_UNRECOGNIZED	(1)
#define INSTRUMENT_SOHO_MDI	(10)
#define INSTRUMENT_SOHO_EIT	(11)
#define INSTRUMENT_GONG_TD	(20)
#define INSTRUMENT_GONG_CT	(21)
#define INSTRUMENT_GONG_TC	(22)
#define INSTRUMENT_GONG_BB	(23)
#define INSTRUMENT_GONG_ML	(24)
#define INSTRUMENT_GONG_LE	(25)
#define INSTRUMENT_GONG_UD	(26)
#define INSTRUMENT_GONG_MERGE	(29)
#define INSTRUMENT_MWO60_MOF	(30)
#define INSTRUMENT_BBSO_SINGER	(40)
#define INSTRUMENT_TRACE	(50)
#define INSTRUMENT_MOTH		(60)
#define INSTRUMENT_OBSPM_SPHG	(70)

#define NO_DATA_DICT	(0x0001)
#define NO_SEMIDIAMETER	(0x0002)
#define NO_XSCALE	(0x0004)
#define NO_YSCALE	(0x0008)
#define NO_XCENTERLOC	(0x0010)
#define NO_YCENTERLOC	(0x0020)
#define NO_HELIO_LATC	(0x0040)
#define NO_HELIO_LONC	(0x0080)
#define NO_HELIO_PA	(0x0100)
#define NO_XUNITS	(0x0200)
#define NO_YUNITS	(0x0400)
#define NO_OBSERVER_LAT	(0x0002)
#define NO_OBSERVER_LON	(0x0004)

#define KEYSCOPE_VARIABLE	(0x80000000)
#define LOCALHS_IMGINFO_VERSION	("1.0")

/*
 *  process image info (attitude, plate scale, distortion)
 *  stub function based on SOI solar_image_params()
 */
typedef struct paramdef {
  double scale;
  double offset;
  double defval;
  unsigned int statusbit;
  char name[32];
} ParamDef;

static double lookup (DRMS_Record_t *rec, ParamDef key, int *status) {
  double value = key.defval;
  int lookupstat = 0;

  value = drms_getkey_double (rec, key.name, &lookupstat);
  value = value * key.scale + key.offset;
  if (lookupstat) *status |= key.statusbit;
  if (isnan (value)) *status |= key.statusbit;
  return value;
}

static char *lookup_str (DRMS_Record_t *rec, ParamDef key, int *status) {
  DRMS_Keyword_t *keywd;
  int lstat;
  char *value;

  value = drms_getkey_string (rec, key.name, &lstat);
  if (lstat) *status |= key.statusbit;
					     /*  cadence should be constant  */
  if ((keywd = drms_keyword_lookup (rec, key.name, 1))) {
    if (keywd->info->recscope != 1) *status |= KEYSCOPE_VARIABLE;
  } else *status |= key.statusbit;
  return value;
}

static int solar_image_info (DRMS_Record_t *img, double *xscl, double *yscl,
    double *ctrx, double *ctry, double *apsd, const char *rsun_key,
    const char *apsd_key, double *pang, double *ellipse_e, double *ellipse_pa,
    int *x_invrt, int *y_invrt, int *need_ephem, int AIPS_convention) {
/*
 *  Provides the following values from the DRMS record:
 *    xscl  scale in the image column direction (arc-sec/pixel)
 *    yscl  scale in the image row direction (arc-sec/pixel)
 *    ctrx  (virtual) fractional pixel column of the center of the
 *	solar disc
 *    ctry  (virtual) fractional pixel row of the center of the
 *	solar disc
 *    apsd  apparent semi-diameter (semimajor-axis) of the solar disc, in
 *	pixel units
 *    pang  position angle of solar north relative to image vertical (y-axis,
 *	[0,0] -> [0,1]), measured westward (clockwise), in radians
 *    eecc  eccentricity of best-fit ellipse describing limb
 *    eang  position angle of best-fit ellipse describing limb, relative
 *	to direction of solar north, measured westward (clockwise), in radians
 *    xinv  0 if image is direct, 1 if flipped by columns
 *    yinv  0 if image is direct, 1 if flipped by rows
 *  If AIPS_convention is true (!=0), it is assumed that the input keywords
 *    representing position and ellipse angles are measured westward
 *    (clockwise) relative to their nominal axes; otherwise they are measured
 *    eastward (counter-clockwise) relative to their nominal axes.
 *  NO!
 *  The following data types are supported: SOHO-MDI, GONG+, Mt. Wilson MOF,
 *    SOHO-EIT, TRACE, BBSO Ha
 *  NO!
 *  If the data are not recognizably of one of these types, the function
 *    returns NO_DATA_DICT; if one or more required keywords are missing
 *    the function returns a status mask indicating which values could not
 *    be filled.  No flag is set if the image ellipse parameters are not
 *    available (unless they are required for other parameters), but the
 *    eccentricity and pericenter angle are quietly set to 0.
 */
  enum param {
    XSCL, YSCL, XUNI, YUNI, LATC, LONC, CTRX, CTRY, PANG, APSD, RSUN,
    ESMA, ESMI, EECC, EANG, XASP, YASP, PCT
  };
  enum known_plat {
    UNKNOWN
  };
  static ParamDef param[PCT];
  static double raddeg =  M_PI / 180.0;
  static double degrad = 180.0 / M_PI;
  double ella, ellb;
  int n, status = 0;
  static int scale_avail, xinv_type, yinv_type;
  static int hdrtype = UNKNOWN, lasthdr = UNKNOWN - 1;
  char *strval;
/*
 *  Set up the appropriate dictionary for interpretation of keywords
 */
  if (lasthdr != hdrtype) {
    if (lasthdr >= UNKNOWN)
      fprintf (stderr,
	  "Warning from solar_image_info(): record keywords may change\n");
    for (n = 0; n < PCT; n++) {
      sprintf (param[n].name, "No Keyword");
      param[n].scale = 1.0;
      param[n].offset = 0.0;
      param[n].defval = NAN;
      param[n].statusbit = 0;
    }
    param[RSUN].statusbit = NO_SEMIDIAMETER;
    param[APSD].statusbit = NO_SEMIDIAMETER;
    param[XSCL].statusbit = NO_XSCALE;
    param[YSCL].statusbit = NO_YSCALE;
    param[XUNI].statusbit = NO_XUNITS;
    param[YUNI].statusbit = NO_YUNITS;
    param[CTRX].statusbit = NO_XCENTERLOC;
    param[CTRY].statusbit = NO_YCENTERLOC;
    param[CTRX].statusbit = NO_XCENTERLOC;
    param[CTRY].statusbit = NO_YCENTERLOC;
    param[LATC].statusbit = NO_HELIO_LATC;
    param[LONC].statusbit = NO_HELIO_LONC;
    param[PANG].statusbit = NO_HELIO_PA;
    param[EECC].defval = 0.0;
    param[EANG].defval = 0.0;
    param[XASP].defval = 1.0;
    param[YASP].defval = 1.0;

    switch (hdrtype) {
      default:					      /*  Assume WCS HPLN/T  */
					 /*  WITH CERTAIN MDI SPECIFIC ONES  */
	scale_avail = 1;
	xinv_type = yinv_type = 0;
	sprintf (param[XUNI].name, "CUNIT1");
	sprintf (param[YUNI].name, "CUNIT2");
	sprintf (param[XSCL].name, "CDELT1");
	sprintf (param[YSCL].name, "CDELT2");
	sprintf (param[CTRX].name, "CRPIX1");
	param[CTRX].offset = -1.0;
	sprintf (param[CTRY].name, "CRPIX2");
	param[CTRY].offset = -1.0;
	*need_ephem = 0;
	strval = lookup_str (img, param[XUNI], &status);
	if (!(status & NO_XUNITS)) {
	  if (!strcmp (strval, "arcsec")) param[XSCL].scale = 1.0;
	  else if (!strcmp (strval, "arcmin")) param[XSCL].scale = 1.0 / 60.0;
	  else if (!strcmp (strval, "deg")) param[XSCL].scale = 1.0 / 3600.0;
	  else if (!strcmp (strval, "mas")) param[XSCL].scale = 1000.0;
	  else if (!strcmp (strval, "rad")) param[XSCL].scale = degrad * 3600.0;
/*
	  need_units = status & KEYSCOPE_VARIABLE;
*/
	}
	if (strval) free (strval);
	strval = lookup_str (img, param[YUNI], &status);
	if (!(status & NO_YUNITS)) {
	  if (!strcmp (strval, "arcsec")) param[YSCL].scale = 1.0;
	  else if (!strcmp (strval, "arcmin")) param[YSCL].scale = 1.0 / 60.0;
	  else if (!strcmp (strval, "deg")) param[YSCL].scale = 1.0 / 3600.0;
	  else if (!strcmp (strval, "mas")) param[YSCL].scale = 1000.0;
	  else if (!strcmp (strval, "rad")) param[YSCL].scale = degrad * 3600.0;
/*
	  need_units = status & KEYSCOPE_VARIABLE;
*/
	}
	if (strval) free (strval);
   /*  the following are appropriate for MDI, but not strictly based on WCS  */
/*
	sprintf (param[RSUN].name, "R_SUN");
	sprintf (param[APSD].name, "OBS_ASD");
*/
	strncpy (param[RSUN].name, rsun_key, 31);
	strncpy (param[APSD].name, apsd_key, 31);
	sprintf (param[PANG].name, "CROTA2");
	param[PANG].scale = -raddeg;
	if (AIPS_convention) param[PANG].scale *= -1;
	sprintf (param[ESMA].name, "S_MAJOR");
	sprintf (param[ESMI].name, "S_MINOR");
	sprintf (param[EANG].name, "S_ANGLE");
	param[EANG].scale = -raddeg;
	if (AIPS_convention) param[EANG].scale *= -1;
    }
  }
  lasthdr = hdrtype;
				    /*  Plate info: image scale, distortion  */
  *apsd = lookup (img, param[RSUN], &status);
  if (scale_avail) {
    *xscl = lookup (img, param[XSCL], &status);
    *yscl = lookup (img, param[YSCL], &status);
    if (status & NO_SEMIDIAMETER) {
      status &= ~NO_SEMIDIAMETER;
      *apsd = lookup (img, param[APSD], &status);
      if (status & NO_SEMIDIAMETER) {
        *need_ephem = 1;
      } else {
	if (!(status & (NO_XSCALE | NO_YSCALE))) {
	  *apsd /= (*xscl <= *yscl) ? *xscl : *yscl;
	  status &= ~NO_SEMIDIAMETER;
	}
      }
    }
  }
  ella = lookup (img, param[ESMA], &status);
  ellb = lookup (img, param[ESMI], &status);
  *ellipse_e = sqrt ((ella - ellb) * (ella + ellb)) / ella;
  *ellipse_pa = lookup (img, param[EANG], &status);
				    /*  Pointing (attitude: image location)  */
  *ctrx = lookup (img, param[CTRX], &status);
  *ctry = lookup (img, param[CTRY], &status);
					  /*  Position angle of solar north  */
  *pang = lookup (img, param[PANG], &status);
                                                      /*  Image orientation  */
  *x_invrt = xinv_type;
  *y_invrt = yinv_type;

  return status;
}

/*
 *  Revision History (all mods by Rick Bogart unless otherwise noted)
 *
 *  09.10.06		version that went into first release under proj/rings
 *  09.12.03		fixed two icc11 compiler warnings; added keyword name
 *		passing to solar_image_info for selected cases
 *  10.04.01		fixed up processing for semidiameter when apsd_key
 *		is missing
 *  10.08.19		make sure status is initialized to 0
 *  10.09.03		fixed initialization of default values for ellipse
 *		elements
 *  10.10.05		fixed define of keyscope_variable (removed ;)
 *  12.04.17		changed return value of pang to conform to Thompson
 *		(solar WCS) convention rather than the AIPS convention for
 *		CROTA2; added argument AIPS_convention to solar_image_info to
 *		be set if the CROTA2 keyword conforms to the AIPS convention
 *  14.10.04		added AIPS-convention correction for ellipse angle
 *  15.06.19		corrected erroneous overwrite of apsd_key
 *			added definition of code version number
 */

