#ifndef __DECODE_HK_VCDU_H
#define __DECODE_HK_VCDU_H

#include "packets.h"

/* Defines for known HMI apids for packets */
#define APID_HMI_SCIENCE_1              400
#define APID_HMI_SCIENCE_2              410
#define APID_HMI_TIME_1                 405
#define APID_HMI_TIME_2                 415
#define APID_HMI_TEST1768_1             407
#define APID_HMI_TEST1768_2             417
#define APID_HMI_TEST256x256_1          408
#define APID_HMI_TEST256x256_2          418
#define APID_HMI_TEST4096x4096_1        409
#define APID_HMI_TEST4096x4096_2        419
#define APID_HMI_IMSTAT_1               445
#define APID_HMI_IMSTAT_2               475
#define APID_HMI_SEQ_1                  451
#define APID_HMI_SEQ_2                  481
#define APID_HMI_OBT_1                  448
#define APID_HMI_OBT_2                  478

/* Defines for known AIA apids for packets */
#define APID_AIA_SCIENCE_1              500
#define APID_AIA_SCIENCE_2              510
#define APID_AIA_TIME_1                 505
#define APID_AIA_TIME_2                 515
#define APID_AIA_TEST1768_1             507
#define APID_AIA_TEST1768_2             517
#define APID_AIA_TEST256x256_1          508
#define APID_AIA_TEST256x256_2          518
#define APID_AIA_IMSTAT_1               529
#define APID_AIA_IMSTAT_2               569
#define APID_AIA_SEQ_1                  536
#define APID_AIA_SEQ_2                  576
#define APID_AIA_OBT_1                  540
#define APID_AIA_OBT_2                  586

/* Return status send back to from 
   decode_next_hk_vcdu to Jim's top lev0 module */
/*********NOTE: CTD: COMMIT TO DRMSi*********/ 
/*********NOTE: WTD: WRITE TO DRMS***********/ 
#define SUCCESS_HK_NEED_TO_CTD            (0)
#define SUCCESS_HK_NEED_TO_WTD_CTD        (1)
#define SUCCESS_SKIP_IMAGE                (4)
#define SUCCESS_SKIP_PROCESSING_APID      (5)
#define ERROR_NODATA                      (-13)
#define ERROR_HK_NO_CONFIG_DATA           (-14)
#define ERROR_HK_CANNOT_FIND_VER_NUM      (-15)
#define ERROR_HK_CANNOT_LOAD_HK_VALUES    (-16)
#define ERROR_HK_CANNOT_LOAD_ENGR_VALUES  (-17)
#define ERROR_HK_INVALID_BITFIELD_LENGTH  (-18)
#define ERROR_HK_UNHANDLED_TYPE           (-19)
#define ERROR_HK_NOSUCHDIR                (-20)
#define ERROR_HK_CANNOT_LOAD_CONFIG       (-21)
#define ERROR_HK_FAILED_WRITE_DAYFILE     (-22)
#define ERROR_HK_ENVIRONMENT_VARS_NOT_SET (-23)
#define ERROR_HK_FAILED_TO_FIND_TIMECODES (-24)
#define ERROR_HK_FAILED_CLOSE_DRMS_RECORD (-25)
#define ERROR_HK_FAILED_OPEN_DRMS_RECORD  (-26)
#define ERROR_HK_FAILED_GETTING_FSN       (-27)

/*  Overall decode_next_hk_vcdu status after saving 
    packet data to files and decode of hk packets for
    each vcdu passed to it from lev 0 top module*/
#define HK_SUCCESS_HK_ALL                 (0)
#define HK_SUCCESS_HK_SOME                (1)

/* return hk_status returned back to decode_next-hk_vcdu 
   from calls to decode_hk(), save_dayfile(), etc. */ 
#define HK_SUCCESS_HKTIME                 (0)
#define HK_SUCCESS_DECODING               (1)
#define HK_SUCCESS_WRITE_DAYFILE          (2)
#define HK_SUCCESS_REACHED_END_VCDU       (3)
#define HK_SUCCESS_SKIP_IMAGE             SUCCESS_SKIP_IMAGE 
#define HK_SUCCESS_SKIP_PROCESSING_APID   SUCCESS_SKIP_PROCESSING_APID

/* return status from write_hk_to_drms */
#define HK_SUCCESS_WROTE_TO_DRMS          (0)

/* When to write day files after getting number of VCDU 
   If set to 1 will write data to day files after every 
   vcdu received. If 100 will write data to day files
   after every 100 vcdu's */
#define HK_WRITE_AFTER_VCDU_COUNT           1
#define HK_INIT_WRITE_FLAG                  1

/* Functions */
void  hk_ccsds_free(CCSDS_Packet_t **p);
int   decode_next_hk_vcdu( unsigned short vcdu[PACKETWORDS], CCSDS_Packet_t **hk_packets, unsigned int **fsn);

#endif
