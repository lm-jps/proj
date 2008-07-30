#ifndef __QUALLEV0_H
#define __QUALLEV0_H
//These are the meaning of the bits put into the QUALLEV0 keyword for a drms
//record from the ingest_lev0 processing.

//Bit0 is the low bit (0x01)
//The first 4 bits are determined from the Img sturct passed back by
//imgdecode()
#define Q_OVFL 0x01	//overflow flag set
#define Q_HDRERR 0x02	//header error flag set
#define Q_CMPERR 0x04	//compression error in image
#define Q_LPXERR 0x08	//last pixel error

//image status packet is missing if FSN != HSQFGSN
#define Q_NOISP 0x10	//no ISP

//missvals is from Img struct totalvals-datavals
#define Q_MISS0 0x100	//missvals > 0
#define Q_MISS1 0x200	//missvals > 0.01*datavals
#define Q_MISS2 0x400	//missvals > 0.05*datavals
#define Q_MISS3 0x800	//missvals > 0.25*datavals

//HMI sepecific
#define Q_SEQERR 0x10000	//sequencer error HSEQERR  != 'SUCCESS'
#define Q_ISSOPEN 0x20000	//ISS loop open HWLTNSET = 'OPEN'
#define Q_HCF1ENCD 0x40000	//Focus/Cal Motor 1 Error
				//HCF1ENCD ne HCF1POS +/- 1
#define Q_HCF2ENCD 0x80000	//Focus/Cal Motor 2 Error
				//HCF2ENCD ne HCF2POS +/- 1
#define Q_HPS1ENCD 0x100000	//Polarization MTR 1 Error
				//HPS1ENCD ne HPL1POS +/- 1
#define Q_HPS2ENCD 0x200000	//Polarization MTR 2 Error
				//HPS2ENCD ne HPL2POS +/- 1
#define Q_HPS3ENCD 0x400000	//Polarization MTR 3 Error
				//HPS3ENCD ne HPL3POS +/- 1
#define Q_HWT1ENCD 0x800000	//Wavelength Motor 1 Error
				//HWT1ENCD ne HWL1POS +/- 1
#define Q_HWT2ENCD 0x1000000	//Wavelength Motor 2 Error
				//HWT2ENCD ne HWL2POS +/- 1
#define Q_HWT3ENCD 0x2000000	//Wavelength Motor 3 Error
				//HWT3ENCD ne HWL3POS +/- 1
#define Q_HWT4ENCD 0x4000000	//Wavelength Motor 4 Error
				//HWT4ENCD ne HWL4POS +/- 1


//AIA sepecific !!!TBD!!!

#endif

