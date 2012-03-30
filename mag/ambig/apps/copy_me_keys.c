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
	"T_OBS", "CADENCE", "DATE_S", "DATE__OBS",											// time
	"DSUN_OBS", "CRLN_OBS", "CRLT_OBS", "CAR_ROT",
	"OBS_VR", "OBS_VW", "OBS_VN", "RSUN_OBS",					// geometry
	"INVCODEV", "INVDOCU", "INVITERA", "INVLMBDA", "INVLMBDF",
	"INVTUNEN", "INVSVDTL", "INVCHIST", "INVPOLTH",
	"INVPJUMP", "INVLMBDM", "INVLMBD0", "INVLMBDB",
	"INVDLTLA", "INVLMBDS", "INVLMBMS",
  "INVLYOTW", "INVWNARW", "INVWSPAC",
	"INVINTTH", "INVNOISE", "INVCONTI", "INVWGHTI",
	"INVWGHTQ", "INVWGHTU", "INVWGHTV", "INVSTLGT", "INVFREEP",
	"INVFLPRF", "INVPHMAP", "INVVLAVE", "INVBLAVE", "INVBBAVE",
	"INVNPRCS", "INVNCNVG",
	"INVKEYS1", "INVKEYS2", "INVKEYS3", 
	"INVKEYI1", "INVKEYI2", "INVKEYI3", 
	"INVKEYD1", "INVKEYD2", "INVKEYD3", 								// inversion
	"BUNIT",
	"BUNIT_000", "BUNIT_001", "BUNIT_002", "BUNIT_003", "BUNIT_004",
	"BUNIT_005", "BUNIT_006", "BUNIT_007", "BUNIT_008", "BUNIT_009",
	"BUNIT_010", "BUNIT_011", "BUNIT_012", "BUNIT_013", "BUNIT_014",
	"BUNIT_015", "BUNIT_016", "BUNIT_017", "BUNIT_018", "BUNIT_019",
	"BUNIT_020", "BUNIT_021", "BUNIT_022", "BUNIT_023", "BUNIT_024",
  "BUNIT_025", "BUNIT_026", "BUNIT_027",
	"HFLID", "HCFTID", "QLOOK", "TINTNUM", "SINTNUM",
	"DISTCOEF", "ROTCOEF", "POLCALM", "SOURCE",
	"CODEVER0", "CODEVER1", "CODEVER2", "CODEVER3", "CODEVER4", "CODEVER5", "CODEVER6"	// misc
};

const char *patchKeys[] = 
{
	"PATCHNUM", "OFFDISK", "QUIET", "MASK",
	"LON_MIN", "LAT_MIN", "LON_MAX", "LAT_MAX", "OMEGA"
};

const char *geoKeys[] =
{
	"CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
	"CDELT1", "CDELT2", "CROTA2",
	"CRDER1", "CRDER2", "CSYSER1", "CSYSER2",
	"IMCRPIX1", "IMCRPIX2"
};

const char *ambigKeys[] =
{
  "AMBCODEV", "AMBDOCU",
	"AMBGMTRY", "AMBWEAK", "AMBNEROD", "AMBNGROW",
	"AMBNPAD", "AMBNAP", "AMBNTX", "AMBNTY",
	"AMBBTHR0", "AMBBTHR1", "AMBSEED", "AMBNEQ",
	"AMBLMBDA", "AMBTFCT0", "AMBTFCTR", "AMBPATCH"
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
		failCount += drms_copykey(outRec, inRec, meKeys[iKey]);
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
		failCount += drms_copykey(outRec, inRec, patchKeys[iKey]);
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
		failCount += drms_copykey(outRec, inRec, geoKeys[iKey]);
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
		failCount += drms_copykey(outRec, inRec, ambigKeys[iKey]);
	}
	
	return failCount;
}
