#ifndef __QUALLEV0_H
#define __QUALLEV0_H
//These are the meaning of the bits put into the QUALLEV0 keyword for a drms
//record from the ingest_lev0 processing.

//Bit0 is the low bit (0x01)
//The first 4 bits are determined from the Img sturct passed back by
//imgdecode()
#define Q_OVFL      0x01	//overflow flag set
#define Q_HDRERR    0x02	//header error flag set
#define Q_CMPERR    0x04	//compression error in image
#define Q_LPXERR    0x08	//last pixel error

//image status packet is missing if FSN != HSQFGSN
#define Q_NOISP     0x10	//no ISP
#define Q_MISSI     0x20	//missing image
#define Q_CORRUPT   0x40	//corrupt image (FSN=469769216 0x1c001c00)
#define Q_INVALTIME 0x80	//T_OBS = 1958.01.01_00:00:00_UTC

//missvals is from Img struct totalvals-datavals
#define Q_MISS0     0x100	//missvals > 0
#define Q_MISS1     0x200	//missvals > 0.01*totalvals
#define Q_MISS2     0x400	//missvals > 0.05*totalvals
#define Q_MISS3     0x800	//missvals > 0.25*totalvals

//IRIS specific
#define Q_DARK      0x10000 	//dark image;    IIFRMTYP = 'DARK'
#define Q_LED       0x40000 	//led image;    IIFRMTYP = 'LED'
#define Q_ISSOPEN   0x20000	//ISS loop open; IISSLOOP = 'OPEN'

#define Q_REOPENED  0x40000000	//image reopened during reconstruction; NPACKETS value may be incorrect
#define Q_MISSALL   0x80000000	//data is completely missing. high bit

//Fits keyword and Image Status Packet (ISP) keyword translation:
//
// ISQFSN   = I_SQ_FRAME_SN
// IISSLOOP = I_ISS_LOOP
// IIFRMTYP = I_IMG_FRAME_TYPE 
//
#endif

