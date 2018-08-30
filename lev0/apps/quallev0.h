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
#define Q_MISSI 0x20	//missing image
#define Q_CORRUPT 0x40	//corrupt image (FSN=469769216 0x1c001c00)
#define Q_INVALTIME 0x80//HOBITSEC = 0 (T_OBS = 1958.01.01_00:00:00_UTC)

//missvals is from Img struct totalvals-datavals
#define Q_MISS0 0x100	//missvals > 0
#define Q_MISS1 0x200	//missvals > 0.01*totalvals
#define Q_MISS2 0x400	//missvals > 0.05*totalvals
#define Q_MISS3 0x800	//missvals > 0.25*totalvals

// Instrument misconfiguration, e.g. worng Aperture Selector target.
#define Q_MISCONFIG 0x1000   // likely set manually
#define Q_INSTR_ANOM 0x2000   //Instrument Anomaly - Entered by hand or ingest_lev0
#define Q_CAM_ANOM 0x8000   //Camera Anomaly - Entered by hand - Hao QUALITY=0x00008000

//HMI sepecific
#define Q_DARK  0x10000 //dark image (bit 16)
//#define Q_SEQERR 0x10000	//sequencer error HSEQERR  != 'SUCCESS'
#define Q_ISSOPEN 0x20000	//ISS loop open HWLTNSET = 'OPEN'
#define Q_HCF1ENCD 0x40000	//Focus/Cal Motor 1 Error
				//HCF1ENCD ne HCF1POS +/- 1
#define Q_HCF2ENCD 0x80000	//Focus/Cal Motor 2 Error
				//HCF2ENCD ne HCF2POS +/- 1
#define Q_HPS1ENCD 0x100000	//Polarization MTR 1 Error
				//HPS1ENCD ne HPL1POS +/- 1 %240
#define Q_HPS2ENCD 0x200000	//Polarization MTR 2 Error
				//HPS2ENCD ne HPL2POS +/- 1 %240
#define Q_HPS3ENCD 0x400000	//Polarization MTR 3 Error
				//HPS3ENCD ne HPL3POS +/- 1 %240
#define Q_HWT1ENCD 0x800000	//Wavelength Motor 1 Error
				//HWT1ENCD ne HWL1POS +/- 1 %240
#define Q_HWT2ENCD 0x1000000	//Wavelength Motor 2 Error
				//HWT2ENCD ne HWL2POS +/- 1 %240
#define Q_HWT3ENCD 0x2000000	//Wavelength Motor 3 Error
				//HWT3ENCD ne HWL3POS +/- 1 %240
#define Q_HWT4ENCD 0x4000000	//Wavelength Motor 4 Error
				//HWT4ENCD ne HWL4POS +/- 1 %240

#define Q_GPREGBIT0 0x10000000
#define Q_GPREGBIT1 0x20000000
#define Q_REOPENED 0x40000000	//image reopened during reconstruction; no impact on data quality except in rare cases where telemetry retransmission may have caused inflated NPACKETS value
#define Q_MISSALL 0x80000000	//data is completely missing. high bit

//Fits keyword and Image Status Packet (ISP) keyword translation:
//HWLTNSET = HMI_IMG_ISS_LOOP
//HSEQERR  = HMI_SEQ_ERROR
//
//HCF1ENCD = HMI_CF1_ENCODER
//HCF2ENCD = HMI_CF1_ENCODER
//HPS1ENCD = HMI_PL1_ENCODER
//HPS2ENCD = HMI_PL2_ENCODER
//HPS3ENCD = HMI_PL3_ENCODER
//HWT1ENCD = HMI_WT1_ENCODER
//HWT2ENCD = HM2_WT1_ENCODER
//HWT3ENCD = HM3_WT1_ENCODER
//HWT4ENCD = HM4_WT1_ENCODER
//
//HCF1POS  = HMI_FSW_CF1_CMDED_TARGET
//HCF2POS  = HMI_FSW_CF2_CMDED_TARGET
//HPL1POS  = HMI_FSW_PL1_CMDED_TARGET
//HPL2POS  = HMI_FSW_PL2_CMDED_TARGET
//HPL3POS  = HMI_FSW_PL3_CMDED_TARGET
//HWL1POS  = HMI_FSW_WT1_CMDED_TARGET
//HWL2POS  = HMI_FSW_WT2_CMDED_TARGET
//HWL3POS  = HMI_FSW_WT3_CMDED_TARGET
//HWL4POS  = HMI_FSW_WT4_CMDED_TARGET

//AIA sepecific 
#define AQ_ISSOPEN 0x20000	//ISS loop open AISTATE = 'OPEN'
#define A94Mech_Err 0x40000	//AIAWVLEN == 94 &&
                            //{(AIFILTYP == 0 && AIFWEN != 269 && AIFWEN != 270)
                            //|| (AIFILTYP == 1 && AIFWEN != 11 && AIFWEN != 12)
                            //|| (AIFILTYP == 2 && AIFWEN != 74 && AIFWEN != 75)
#define A131Mech_Err 0x80000	//AIAWVLEN == 131 &&
                            //{(AIFILTYP == 0 && AIFWEN != 269 && AIFWEN != 270)
                            //|| (AIFILTYP == 1 && AIFWEN != 11 && AIFWEN != 12)
                            //|| (AIFILTYP == 2 && AIFWEN != 74 && AIFWEN != 75)
#define A171Mech_Err 0x100000	//AIAWVLEN == 171 &&
                            //{(AIFILTYP == 0 && AIFWEN != 203 && AIFWEN != 204)
                            //|| (AIFILTYP == 1 && AIFWEN != 11 && AIFWEN != 12)
#define A193Mech_Err 0x200000	//AIAWVLEN == 193 && {AIASEN != 6
                          //|| (AIFILTYP == 0 && AIFWEN != 269 && AIFWEN != 270)
                          //|| (AIFILTYP == 1 && AIFWEN !=  11 && AIFWEN !=  12)
                          //|| (AIFILTYP == 2 && AIFWEN != 74 && AIFWEN != 75)
#define A211Mech_Err 0x400000	//AIAWVLEN == 211 && {AIASEN != 24
                          //|| (AIFILTYP == 0 && AIFWEN != 203 && AIFWEN != 204)
                          //|| (AIFILTYP == 1 && AIFWEN != 137 && AIFWEN != 138)
                          //|| (AIFILTYP == 2 && AIFWEN != 74 && AIFWEN != 75)
#define A304Mech_Err 0x800000	//AIAWVLEN == 304 &&
                          //  {(AIFILTYP == 0 && AIFWEN != 203 && AIFWEN != 204)
                          //|| (AIFILTYP == 1 && AIFWEN != 137 && AIFWEN != 138)
                          //|| (AIFILTYP == 2 && AIFWEN != 74 && AIFWEN != 75)
#define A335Mech_Err 0x1000000	//AIAWVLEN == 335 &&
                          //  {(AIFILTYP == 0 && AIFWEN != 203 && AIFWEN != 204)
                          //|| (AIFILTYP == 1 && AIFWEN != 137 && AIFWEN != 138)
                          //|| (AIFILTYP == 2 && AIFWEN != 74 && AIFWEN != 75)
#define A160Mech_Err 0x2000000	//AIAWVLEN == 1600 && AIFWEN != 269 && 
			  // AIFWEN != 270
#define A170Mech_Err 0x4000000	//AIAWVLEN == 1700 && AIFWEN != 137 && 
			  // AIFWEN != 138
#define A450Mech_Err 0x8000000	//AIAWVLEN == 4500 && AIFWEN !=  74 && 
			  // AIFWEN !=  75

#define AQ_INVAL_WL 0x10000000 //invalid wavelength WAVE_STR == "UNKNOWN"

//AIA Mechanism position definitions from Paul Boerner
//WAVELEN FILTER_TYPE     FW_ENCODER      AS_ENCODER"
//1600    "Don't check"   "269 or 270"    "Don't check"
//        "Don't check"   "269 or 270"    "Don't check"
//        "Don't check"   "269 or 270"    "Don't check"
//
//1700    "Don't check"   "137 or 138"    "Don't check"
//        "Don't check"   "137 or 138"    "Don't check"
//        "Don't check"   "137 or 138"    "Don't check"
//
//4500    "Don't check"   "74 or 75"      "Don't check"
//        "Don't check"   "74 or 75"      "Don't check"
//        "Don't check"   "74 or 75"      "Don't check"
//
//WAVELEN FILTER_TYPE     FW_ENCODER      AS_ENCODER"
//94      0               "269 or 270"    "Don't check"
//        1               "11 or 12"      "Don't check"
//        2               "74 or 75"      "Don't check"
//
//133     0               "269 or 270"    "Don't check"
//        1               "11 or 12"      "Don't check"
//        2               "74 or 75"      "Don't check"
//
//171     0               "203 or 204"    "Don't check"
//        1               "11 or 12"      "Don't check"
//        2               "Don't Check"   "Don't check"
//
//304     0               "203 or 204"    "Don't check"
//        1               "137 or 138"    "Don't check"
//        2               "74 or 75"      "Don't check"
//
//335     0               "203 or 204"    "Don't check"
//        1               "137 or 138"    "Don't check"
//        2               "74 or 75"      "Don't check"
//
//WAVELEN FILTER_TYPE     FW_ENCODER      AS_ENCODER"
//193     0               "269 or 270"    6
//        1               "11 or 12"      6
//        2               "74 or 75"      6
//
//211     0               "203 or 204"    24
//        1               "137 or 138"    24
//        2               "74 or 75"      24
//
// ***************************************************************
//
//Fits keyword and Image Status Packet (ISP) keyword translation:
//
//ASQFSN          AIA_SEQ_FRAME_SN        longlong
//AISTATE         AIA_IMG_ISS_LOOP        string
//AIAWVLEN        AIA_IMG_WAVELENGTH      int
//AIASEN          AIA_IMG_AS_ENCODER      int
//AIFILTYP        AIA_IMG_FILTER_TYPE     short
//AIFWEN          AIA_IMG_FW_ENCODE       int
//AIFOENFL        AIA_IMG_FOCUS_ENA_FLAG  short
//
#endif

