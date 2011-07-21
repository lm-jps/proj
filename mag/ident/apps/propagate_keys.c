/*
 *	Propagate keywords from:
 *        M_720s -> Marmask 
 *        Ic_720s -> Marmask 
 *        M_720s -> Mharp
 *
 *	Upon success, return 0
 *	Otherwise return number of keywords failed to populate
 *	If record fails, return -1
 *
 *	Prime keys like T_REC, PNUM need to be set upfront in the main code
 *
 *	Usage:
 *	For Marmask, use:
 *		propagate_keys_mag2mask()
 *		propagate_keys_intensity2mask()
 *	For Mharp population, use:
 *		propagate_keys_harp()
 *
 *	Written by X. Sun, Jun 27 2011
 *      Contributions by M. Turmon, Jul 2011
 *
 */

// (We #includ the present .c in another .c that already includes the files below)
// #include "jsoc_main.h"
// #include "drms_types.h"

/*
 * These keys arrays are sequences of pairs,
 *   <source_key_name>, <dest_key_name>
 * where if dest == NULL, source is mapped to the same name.
 * The arrays must be null-terminated.
 */
static const char *mag2mask_keys[] = {
  // info
  "QUALITY",  NULL,
  // "QUAL_S",   NULL, // not in M_720s (7/2011)
  "QUALLEV1", NULL,
  "INSTRUME", NULL,
  "CAMERA",   NULL,
  "SATVALS",  NULL,
  // "CADENCE", // omit since it's fixed for the sets we use
  // time
  "T_OBS",    NULL,
  // "DATE_S",   NULL, // omit; this key is not in M_720s (7/2011)
  "DATE__OBS",NULL,
  // observer's geometry
  "DSUN_OBS", NULL,
  "RSUN_OBS", NULL,
  "CRLN_OBS", NULL,
  "CRLT_OBS", NULL,
  "CAR_ROT",  NULL,
  "OBS_VR",   NULL,
  "OBS_VW",   NULL,
  "OBS_VN",   NULL,
  "RSUN_OBS", NULL,
  // wcs
  "CRPIX1",   NULL,
  "CRPIX2",   NULL,
  "CRVAL1",   NULL,
  "CRVAL2",   NULL,
  "CDELT1",   NULL,
  "CDELT2",   NULL,
  "CROTA2",   NULL,
  "CRDER1",   NULL,
  "CRDER2",   NULL,
  "CSYSER1",  NULL,
  "CSYSER2",  NULL,
  // lev1 calib
  "HFLID",    NULL,
  "HCFTID",   NULL,
  "QLOOK",    NULL,
  "TINTNUM",  NULL,
  "SINTNUM",  NULL,
  "DISTCOEF", NULL,
  "ROTCOEF",  NULL,
  "POLCALM",  NULL,
  // "SOURCE",   NULL, // omit, it's huge and you can track back to mgram
  // versioning
  "CODEVER0", NULL,
  "CODEVER1", NULL,
  "CODEVER2", NULL,
  "CODEVER3", NULL,
  "CODEVER4", NULL, // not in M_720s, but maybe in future?
  // "CODEVER5", NULL,
  // end sentinel -- required
  NULL,       NULL,
  };

static const char *intensity2mask_keys[] = {
  // code versioning
  "CODEVER0", "ICCODEV0",
  "CODEVER1", "ICCODEV1",
  "CODEVER2", "ICCODEV2",
  "CODEVER3", "ICCODEV3",
  "CODEVER4", "ICCODEV4",
  // limb darkening function
  "NORMALIZE",NULL,
  "COEF_VER", NULL,
  "LDCoef1",  NULL,
  "LDCoef2",  NULL,
  "LDCoef3",  NULL,
  "LDCoef4",  NULL,
  "LDCoef5",  NULL,
  // end sentinel -- required
  NULL,       NULL,
  };

static const char *harp_keys[] = {
  "ARMCODEV", NULL,
  "ARMDOCU",  NULL,
  // end sentinel -- required
  NULL,       NULL,
  };



/*
 * generic key propagation
 */
static
int 
propagate_keys_generic(const char **keys, DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  int failCount;
  int status;
  int iKey;
  const int verbose = 0;
  DRMS_Value_t srcval; /* owns contained string */


  // keys is a sequence of pairs <source_key_name>, <dest_key_name>
  // where if dest == NULL, source is mapped to the same name.
  for (failCount = iKey = 0; keys[iKey] != NULL; iKey += 2) {
    status = DRMS_SUCCESS;
    srcval = drms_getkey_p(inRec, keys[iKey], &status);
    if (status == DRMS_SUCCESS) {
      status = drms_setkey_p(outRec, 
                             keys[iKey+1] ? keys[iKey+1] : keys[iKey], 
			     &srcval);
    }
    drms_value_free(&srcval);
    /*
    status = drms_copykey(outRec, inRec, keys[iKey]);
    */
    failCount += (status == DRMS_SUCCESS) ? 0 : 1;
    if (verbose) {
      if (status != DRMS_SUCCESS)
	printf("prop-keys: error for %s was %d\n", keys[iKey], status);
      else
	printf("prop-keys: OK for %s\n", keys[iKey]);
    }
  }
  return failCount;
}

/*
 * propagate keys in the magnetogram into a mask
 */
static
int 
propagate_keys_mag2mask(DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  return propagate_keys_generic(mag2mask_keys, inRec, outRec);
}

/*
 * propagate keys in the intensitygram into a mask
 */
static
int 
propagate_keys_intensity2mask(DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  return propagate_keys_generic(intensity2mask_keys, inRec, outRec);
}


/*
 * propagate keys in the magnetogram into a HARP
 */
static
int 
propagate_keys_harp(DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  int failCount;
  int status;
  int iKey;
	
  if (inRec == NULL || outRec == NULL)
    return -1;
  failCount = propagate_keys_mag2mask(inRec, outRec);
  if (failCount == -1) 
    return -1;
  for (iKey = 0; harp_keys[iKey] != NULL; iKey++) {
    status = drms_copykey(outRec, inRec, harp_keys[iKey]);
    failCount += (status == DRMS_SUCCESS) ? 0 : 1;
  }
	
  status = 0;
  // FIXME: can use prop_keys_generic here
  failCount += drms_setkey_float(outRec, "IMCRPIX1", 
				 drms_getkey_float(inRec, "CRPIX1", &status));
  failCount += drms_setkey_float(outRec, "IMCRPIX2", 
				 drms_getkey_float(inRec, "CRPIX2", &status));
  failCount += drms_setkey_float(outRec, "IMCRVAL1", 
				 drms_getkey_float(inRec, "CRVAL1", &status));
  failCount += drms_setkey_float(outRec, "IMCRVAL2", 
				 drms_getkey_float(inRec, "CRVAL2", &status));
  return failCount;
}
