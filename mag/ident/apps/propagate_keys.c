/*
 * Propagate keywords from:
 *    M_720s  -> Marmask 
 *    Ic_720s -> Marmask 
 *    M_720s  -> Mharp
 *
 * Upon success, return 0
 *   Otherwise return number of keywords failed to populate
 *   If record fails, return -1
 *
 * Prime keys, e.g., T_REC, need to be set upfront in the main code
 *
 * Usage:
 *   For Marmask, use:
 *	propagate_keys_mag2mask()
 *	propagate_keys_intensity2mask()
 *   For Mharp population, use:
 *	propagate_keys_harp()
 *
 *  Written by X. Sun, Jun 27 2011
 *  Rewritten, M. Turmon, Jul 2011, Aug 2013
 *
 */

// (We #include the present .c in another .c that already includes the files below)
// #include "jsoc_main.h"
// #include "drms_types.h"

/*
 * These keys arrays are sequences of pairs,
 *   <source_key_name>, <dest_key_name>
 * where if dest == NULL, source is mapped to the same name.
 * The arrays must be null-terminated.
 *
 * Keys that are constant in DRMS should not be present here.
 */

// shared keys:
//   hmi.M_720s -> hmi.Marmasm_720s
//   mdi.fd_M_96m_lev182 -> mdi.armask
// note: keep this really simple: just time and geometry
static const char *mag2mask_shared_keys[] = {
  // time
  "T_OBS",    NULL,
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
  // end sentinel -- required
  NULL,       NULL,
  };

// hmi.M_720s -> hmi.Marmask_720s
static const char *mag2mask_hmi_keys[] = {
  "INSTRUME", NULL,    // variable
  "CAMERA",   NULL,    // variable
  // "CADENCE",  NULL, // omit since it's fixed for the sets we use
  "QUALITY",  NULL,
  // "QUAL_S",   NULL, // not in M_720s (7/2011)
  "QUALLEV1", NULL,
  "SATVALS",  NULL,
  // wcs extras
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
  "CALVER64", NULL, // new 11/2012
  // "SOURCE",   NULL, // omit, it's huge and you can track back to mgram
  // versioning
  "CODEVER0", NULL,
  "CODEVER1", NULL,
  "CODEVER2", NULL,
  "CODEVER3", NULL,
  // "CODEVER4", NULL, // not in M_720s (7/2013)
  // "CODEVER5", NULL,
  // end sentinel -- required
  NULL,       NULL,
  };

// mdi.fd_M_96m_lev182 -> mdi.fd_armask
static const char *mag2mask_mdi_keys[] = {
  "QUALITY",  NULL,
  "R_SUN",    NULL,
  "X0",       NULL,
  "Y0",       NULL,
  "EARTH_DT", NULL,
  "RUNTIME",  "MRUNTIME",
  // end sentinel -- required
  NULL,       NULL,
  };

// hmi and mdi (none at present)
static const char *intensity2mask_shared_keys[] = {
  // end sentinel -- required
  NULL,       NULL,
  };

// hmi.Ic_noLimbDark_720s -> hmi.Marmask_720s
static const char *intensity2mask_hmi_keys[] = {
  // code versioning
  "CODEVER0", "ICCODEV0",
  "CODEVER1", "ICCODEV1",
  "CODEVER2", "ICCODEV2",
  "CODEVER3", "ICCODEV3",
  "CODEVER4", "ICCODEV4",
  // limb darkening function
  // "NORMALIZE",NULL, // not in Marmask_720s jsd
  "COEF_VER", NULL,
  "LDCoef0",  NULL,
  "LDCoef1",  NULL,
  "LDCoef2",  NULL,
  "LDCoef3",  NULL,
  "LDCoef4",  NULL,
  "LDCoef5",  NULL,
  // end sentinel -- required
  NULL,       NULL,
  };

// mdi.fd_InterpolatedIc -> mdi.fd_armask
//   note, for MDI, this is *interpolated* Ic's, not ordinary Ic's
//   see the mdi_interpolate_ic module
static const char *intensity2mask_mdi_keys[] = {
  // end sentinel -- required
  NULL,       NULL,
  };


/*
 *  keyword map: mag -> harp
 */

// shared keys:
//   hmi.M_720s -> hmi.Mharp_720s
//   mdi.fd_M_96m_lev182 -> mdi.tarp
// note: keep this really simple: just time and geometry
static const char *mag2harp_shared_keys[] = {
  // time
  "T_OBS",    NULL,
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
  //   the CRPIXn offsets, for fulldisk WCS, are not correct for HARP WCS
  //   the CRVALn values are the same (?)
  "CRVAL1",   NULL,
  "CRVAL2",   NULL,
  "CDELT1",   NULL,
  "CDELT2",   NULL,
  "CROTA2",   NULL,
  // save the full-disk WCS, but do not confuse with HARP WCS
  "CRPIX1",  "IMCRPIX1", 
  "CRPIX2",  "IMCRPIX2",
  "CRVAL1",  "IMCRVAL1",
  "CRVAL2",  "IMCRVAL2", 
  // end sentinel -- required
  NULL,       NULL,
  };

// hmi-only:  hmi.M_720s -> hmi.Mharp_720s
static const char *mag2harp_hmi_keys[] = {
  "INSTRUME", NULL,
  "CAMERA",   NULL,
  // "CADENCE",  NULL, // omit: it's fixed for the sets we use
  "QUALITY",  NULL,
  // "QUAL_S",   NULL, // not in M_720s (7/2011)
  "QUALLEV1", NULL,
  // WCS extras
  "CRDER1",   NULL,
  "CRDER2",   NULL,
  "CSYSER1",  NULL,
  "CSYSER2",  NULL,
  // hmi-specific lev1 calib -- still valid to include in HARP
  "HFLID",    NULL,
  "HCFTID",   NULL,
  "QLOOK",    NULL,
  "TINTNUM",  NULL,
  "SINTNUM",  NULL,
  "DISTCOEF", NULL,
  "ROTCOEF",  NULL,
  "POLCALM",  NULL,
  "CALVER64", NULL, // new 11/2012
  // "SOURCE",   NULL, // omit, it's huge and you can track back to mgram
  // versioning
  "CODEVER0", NULL,
  "CODEVER1", NULL,
  "CODEVER2", NULL,
  "CODEVER3", NULL,
  // "CODEVER4", NULL, // not in M_720s, but maybe in future?
  // "CODEVER5", NULL,
  // end sentinel -- required
  NULL,       NULL,
  };

// mdi-only: mdi.fd_M_96m_lev182 -> mdi.tarp
static const char *mag2harp_mdi_keys[] = {
  "EARTH_DT", NULL,
  "RUNTIME",  "MRUNTIME",
  // end sentinel -- required
  NULL,       NULL,
  };

/*
 *   mask -> harp
 *
 * shared, HMI-only, MDI-only
 */
// shared
//   hmi.Marmask_720s -> hmi.Mharp_720s
//   mdi.armask -> mdi.tarp
static const char *mask2harp_shared_keys[] = {
  // mask values
  "OFFDISK",  NULL,
  "QUIET",    NULL,
  "ACTIVE",   NULL,
  "NCLASS",   NULL,
  "ON_PATCH", NULL,
  "ON_PATCH", "MASK", // two versions of this
  // mask quality
  "ARM_QUAL", NULL,
  "ARM_NCLN", NULL,
  // code versions
  "ARMCODEV", NULL,
  "ARMDOCU",  NULL,
  // most important mask parameters
  "ARM_MODL", NULL,
  "ARM_BETA", NULL,
  "ARM_EDGE", NULL,
  // end sentinel -- required
  NULL,       NULL,
  };

// hmi-only: hmi.Marmask_720s -> hmi.Mharp_720s
static const char *mask2harp_hmi_keys[] = {
  // (empty for now)
  // end sentinel -- required
  NULL,       NULL,
  };

// mdi-only: mdi.armask -> mdi.tarp
static const char *mask2harp_mdi_keys[] = {
  // TBD
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

  // should not happen
  if (inRec == NULL || outRec == NULL)
    return -1;
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
 *   mode: "HMI" or "MDI"
 */
static
int 
propagate_keys_mag2mask(char *mode, DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  int s1, s2;

  s1 = propagate_keys_generic(mag2mask_shared_keys, inRec, outRec);
  if (strcmp(mode, "HMI") == 0)
    s2 = propagate_keys_generic(mag2mask_hmi_keys, inRec, outRec);
  else if (strcmp(mode, "MDI") == 0)
    s2 = propagate_keys_generic(mag2mask_mdi_keys, inRec, outRec);
  else
    s2 = 0;
  if (s1 < 0 || s2 < 0)
    return -1; // error
  else
    return s1 + s2; // missing key count
}

/*
 * propagate keys in the intensitygram into a mask
 *   mode: "HMI" or "MDI"
 */
static
int 
propagate_keys_intensity2mask(char *mode, DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  int s1, s2;

  s1 = propagate_keys_generic(intensity2mask_shared_keys, inRec, outRec);
  if (strcmp(mode, "HMI") == 0)
    s2 = propagate_keys_generic(intensity2mask_hmi_keys, inRec, outRec);
  else if (strcmp(mode, "MDI") == 0)
    s2 = propagate_keys_generic(intensity2mask_mdi_keys, inRec, outRec);
  else
    s2 = 0;
  if (s1 < 0 || s2 < 0)
    return -1; // error
  else
    return s1 + s2; // missing key count
}


/*
 * propagate keys in the mag to a HARP
 *   (not intended as an entry point)
 */
static
int 
propagate_keys_mag2harp(char *mode, DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  int s1, s2;

  s1 = propagate_keys_generic(mag2harp_shared_keys, inRec, outRec);
  if (strcmp(mode, "HMI") == 0)
    s2 = propagate_keys_generic(mag2harp_hmi_keys, inRec, outRec);
  else if (strcmp(mode, "MDI") == 0)
    s2 = propagate_keys_generic(mag2harp_mdi_keys, inRec, outRec);
  else
    s2 = 0;
  if (s1 < 0 || s2 < 0)
    return -1; // error
  else
    return s1 + s2; // missing key count
}


/*
 * propagate keys in the mask to a HARP
 *   (not intended as an entry point)
 */
static
int 
propagate_keys_mask2harp(char *mode, DRMS_Record_t *inRec, DRMS_Record_t *outRec)
{
  int s1, s2;

  s1 = propagate_keys_generic(mask2harp_shared_keys, inRec, outRec);
  if (strcmp(mode, "HMI") == 0)
    s2 = propagate_keys_generic(mask2harp_hmi_keys, inRec, outRec);
  else if (strcmp(mode, "MDI") == 0)
    s2 = propagate_keys_generic(mask2harp_mdi_keys, inRec, outRec);
  else
    s2 = 0;
  if (s1 < 0 || s2 < 0)
    return -1; // error
  else
    return s1 + s2; // missing key count
}


/*
 * propagate keys in the (mag, mask) into a HARP
 *
 * (intended entry point)
 */
static
int 
propagate_keys_harp(char *mode, DRMS_Record_t *magRec, DRMS_Record_t *maskRec, DRMS_Record_t *outRec)
{
  int s1, s2;
	
  // harp keys += mag keys
  s1 = propagate_keys_mag2harp( mode, magRec,  outRec);
  // harp keys += mask keys
  s2 = propagate_keys_mask2harp(mode, maskRec, outRec);
  // (can also other keys as necessary)
  if (s1 < 0 || s2 < 0)
    return -1; // error
  else
    return s1 + s2; // missing key count
}

