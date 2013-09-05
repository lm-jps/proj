#ifndef __QUALLEV1_H
#define __QUALLEV1_H
//These are the meaning of the bits put into the QUALITY keyword for a drms
//record from the build_lev1 processing.

//Bit0 is the low bit (0x01)
#define Q_NOFLAT 0x01	//flatfield not available or error
#define Q_NOORB  0x02	//orbit data not available or error
#define Q_NOASD  0x04	//ancillary sci data not available or error
#define Q_NOMPD  0x08	//master pointing data not available or error
#define Q_NOLIMB 0x10	//limb fit error

//#bit 5
//#bit 6
//#bit 7

//missvals is from Img struct totalvals-datavals
#define Q_1_MISS0 0x100	//missvals > 0
#define Q_1_MISS1 0x200	//missvals > 0.01*totalvals
#define Q_1_MISS2 0x400	//missvals > 0.05*totalvals
#define Q_1_MISS3 0x800	//missvals > 0.25*totalvals

#define Q_NOACS_SCI 0x1000	//ACS_MODE != 'SCIENCE'
#define Q_ACS_ECLP 0x2000	//ACS_ECLP == 'YES'
#define Q_ACS_SUNP 0x4000	//ACS_SUNP == 'NO' no sun presence
#define Q_ACS_SAFE 0x8000	//ACS_SAFE == 'YES' . safemode flag set
#define Q_IMG_TYPE 0x10000	//Dark image
#define Q_LOOP_OPEN 0x20000	//HWLTNSET == "OPEN" or AISTATE == "OPEN"
#				//ISS Loop Open
//#18    Calibration Image	//based on FID range. See Rock
#define Q_CAL_IMG 0x40000	//Calibration image
#define Q_CALM_IMG 0x80000	//HMI cal mode image
#define Q_AIA_FOOR 0x100000	//AIA focus out of range
#define Q_AIA_REGF 0x200000	//AIA register flag

#define Q_THERM_RECOV 0x400000	//HMI thermal recovery
#define Q_LUNAR_TRAN 0x800000	//HMI lunar transit

#define Q_NRT	  0x40000000	//near real time mode (formerly quicklook)
#define Q_MISSALL 0x80000000    //Image not available. high bit

#endif

