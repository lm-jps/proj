#include "jsoc_main.h"
#include "drms_types.h"

/*	Populate keywords from ME inversion results to 
 *	ME_patch, B, B_patch, ME_HARP, B_HARP, etc.
 *
 *	Upon success, return 0
 *	Otherwise return number of keywords failed to populate
 *	If record fails, return -1
 *
 *	Prime keys like T_REC, PNUM needs to be set upfront in the main code
 *
 *	Usage:
 *	For full disk image population, use:
 *		copy_me_keys() and copy_geo_keys()
 *	For patch population without remapping, use:
 *		copy_me_keys() and copy_geo_keys() and copy_patch_keys()
 *	For patch with remapping, use:
 *		copy_me_keys() and copy_patch_keys()
 *	For B series, additionally use:
 *		copy_ambig_keys()
 *
 *	Written by X. Sun, Mar 01 2011
 */

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

const char *meKeys[] =
{
	"QUALITY", "QUAL_S", "QUALLEV1", "INSTRUME", "CAMERA",		// info
	"T_OBS", "CADENCE",											// time
	"DSUN_OBS", "CRLN_OBS", "CRLT_OBS", "CAR_ROT",
	"OBS_VR", "OBS_VW", "OBS_VN", "RSUN_OBS",					// geometry
	"INVCODEV", "INVITERA", "INVLMBDA", "INVLMBDF",
	"INVTUNEN", "INVSVDTL", "INVCHIST", "INVPOLTH",
	"INVPJUMP", "INVLMBDM", "INVLMBD0", "INVLMBDB",
	"INVDLTLA", "INVLYOTW", "INVWNARW", "INVWSPAC",
	"INVINTTH", "INVNOISE", "INVCONTI", "INVWGHTI",
	"INVWGHTQ", "INVWGHTU", "INVWGHTV", "INVSTLGT",
	"INVFLPRF", "INVPHMAP", "INVBLAVE", "INVBBAVE",
	"INVNPRCS", "INVNCNVG",										// inversion
	"HFLID", "HCFTID", "QLOOK", "TINTNUM", "SINTNUM",
	"DISTCOEF", "ROTCOEF", "POLCALM", "SOURCE",
	"CODEVER0", "CODEVER1", "CODEVER2", "CODEVER3", "CODEVER4"	// misc
};

const char *patchKeys[] = 
{
	"PATCHNUM", "OFFDISK", "QUIET", "MASK"
};

const char *geoKeys[] =
{
	"CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
	"CDELT1", "CDELT2", "CROTA2",
	"CRDER1", "CRDER2", "CSYSER1", "CSYSER2"
};

const char *ambigKeys[] =
{
	"AMBGMTRY", "AMBWEAK", "AMBNEROD", "AMBNGROW",
	"AMBNPAD", "AMBNAP", "AMBNTX", "AMBNTY",
	"AMBBTHR0", "AMBBTHR1", "AMBSEED", "AMBNEQ",
	"AMBLMBDA", "AMBTFCT0", "AMBTVCTR"
};



int copy_me_keys (DRMS_Record_t *inRec, DRMS_Record_t *outRec)

{
	int failCount = 0;
	int iKey, nKeys = ARRLENGTH(meKeys);
	
	if (inRec == NULL || outRec == NULL) {
		return -1;
	}
	
	for (iKey = 0; iKey < nKeys; iKey++)
	{
		failCount += drms_copykey(inRec, outRec, meKeys[iKey]);
	}
	
	return failCount;
}



int copy_patch_keys (DRMS_Record_t *inRec, DRMS_Record_t *outRec)

{
	int failCount = 0;
	int iKey, nKeys = ARRLENGTH(patchKeys);
	
	if (inRec == NULL || outRec == NULL) {
		return -1;
	}
	
	for (iKey = 0; iKey < nKeys; iKey++)
	{
		failCount += drms_copykey(inRec, outRec, patchKeys[iKey]);
	}
	
	return failCount;
}



int copy_geo_keys (DRMS_Record_t *inRec, DRMS_Record_t *outRec)

{
	int failCount = 0;
	int iKey, nKeys = ARRLENGTH(geoKeys);
	
	if (inRec == NULL || outRec == NULL) {
		return -1;
	}
	
	for (iKey = 0; iKey < nKeys; iKey++)
	{
		failCount += drms_copykey(inRec, outRec, geoKeys[iKey]);
	}
	
	return failCount;
}



int copy_ambig_keys (DRMS_Record_t *inRec, DRMS_Record_t *outRec)

{
	int failCount = 0;
	int iKey, nKeys = ARRLENGTH(ambigKeys);
	
	if (inRec == NULL || outRec == NULL) {
		return -1;
	}
	
	for (iKey = 0; iKey < nKeys; iKey++)
	{
		failCount += drms_copykey(inRec, outRec, ambigKeys[iKey]);
	}
	
	return failCount;
}
