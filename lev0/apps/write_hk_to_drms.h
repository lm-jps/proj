#ifndef __WRITE_HK_TO_DRMS_H
#define __WRITE_HK_TO_DRMS_H

/* General size defines for arrays */
#define HK_LEV0_MAX_PVN_STR              10
#define HK_LEV0_MAX_JVN_STR              10
#define HK_LEV0_MAX_FILE_NAME           300
#define HK_LEV0_MAXLINE_IN_FILE         200
#define HK_LEV0_MAX_DSNAME_STR          100
#define HK_LEV0_MAX_PROJECT_NAME         50
#define HK_LEV0_MAX_DATATYPE_NAME        50
#define HK_LEV0_MAX_LONG_KW_NAME         50
#define HK_LEV0_MAX_PACKET_NAME          25
#define HK_LEV0_MAX_HKS_DATE             50
#define HK_LEV0_MAX_HKS                 100
#define HK_LEV0_MAX_LU_PREFIX_FILENAME   25

/* Identify ranges of apids for hmi,aia, and sdo */
#define HK_HSB_LOWEST_HMI_APID       400
#define HK_HSB_HIGHEST_HMI_APID      499
#define HK_HSB_LOWEST_AIA_APID       500
#define HK_HSB_HIGHEST_AIA_APID      599
#define HK_LR_LOWEST_SDO_APID         96
#define HK_LR_HIGHEST_SDO_APID       399
#define HK_LR_LOWEST_HMI_APID          1
#define HK_LR_HIGHEST_HMI_APID        31
#define HK_LR_LOWEST_AIA_APID         32
#define HK_LR_HIGHEST_AIA_APID        63

/* identify PACKET APID for low rate(LR) and high speed bus(HSB) */
#define HK_LR_HMI_ISP                 29
#define HK_LR_HMI_SEQ                 21
#define HK_LR_HMI_OBT                 18
#define HK_LR_AIA_ISP                 39
#define HK_LR_AIA_SEQ                 46
#define HK_LR_AIA_OBT                 50
#define HK_LR_SDO_ASD                129
#define HK_HSB_HMI_ISP_1             445
#define HK_HSB_HMI_ISP_2             475
#define HK_HSB_HMI_SEQ_1             451
#define HK_HSB_HMI_SEQ_2             481
#define HK_HSB_HMI_OBT_1             448
#define HK_HSB_HMI_OBT_2             478
#define HK_HSB_AIA_ISP_1             529
#define HK_HSB_AIA_ISP_2             569
#define HK_HSB_AIA_SEQ_1             536
#define HK_HSB_AIA_SEQ_2             576
#define HK_HSB_AIA_OBT_1             540
#define HK_HSB_AIA_OBT_2             580

/* Identify packet names */
#define HK_HMI_PKT_NAME_ISP        "isp"
#define HK_HMI_PKT_NAME_SEQ        "seq"
#define HK_HMI_PKT_NAME_OBT        "obt"
#define HK_AIA_PKT_NAME_ISP        "isp"
#define HK_AIA_PKT_NAME_SEQ        "seq"
#define HK_AIA_PKT_NAME_OBT        "obt"
#define HK_SDO_PKT_NAME_ASD        "asd"

/* other parameters */
#define HK_SECONDS_PER_DAY        86400
#define HK_INIT_MERGE_MAP_FILE_FLAG (0)
#define HK_MERGE_MAP_FILE_FLAG      (1)
#define HK_NON_MERGE_MAP_FILE_FLAG  (2)

/* check hk record exits parameter */
#define HK_HIGH_QUERY_RANGE         (1)
#define HK_LOW_QUERY_RANGE          (0)
#define HK_MAX_SIZE_RANGE_TIME     (50)
#define HK_MAX_SIZE_QUERY         (200)


/* Structures */

/***************** Map Data struct *********************************/
typedef struct Map_Data_struct 
{
  char pvn[HK_LEV0_MAX_PVN_STR];
  char jvn[HK_LEV0_MAX_JVN_STR];
  char dsn[HK_LEV0_MAX_DSNAME_STR];
  struct Map_Data_struct *next;
} Map_Data_t;
/***************** JSOC_Version_Map struct *************************/
typedef struct JSOC_Version_Map_struct 
{
  short apid;
  Map_Data_t *mdata;
  struct JSOC_Version_Map_struct *next;
} JSOC_Version_Map_t;
/***************  HK_timecode_inmemory struct **********************/
/*Link List of second and subsecond timecodes cache here from data series */
typedef struct HK_timecode_inmemory
{
  long int sec;
  long int subsec;
  struct HK_timecode_inmemory *next;
}HK_Timecode_t;
/***************  HK_dsn_range_inmemory struct *********************/
/*Link list structure for each data series name and contains low range*/ 
/*timecode in seconds and pointer to cached time code in data series  */
typedef struct HK_dsn_range_inmemory
{
  char dsname[HK_LEV0_MAX_DSNAME_STR];
  long int timecode_lrsec;
  long int timecode_hrsec;
  struct HK_dsn_range_inmemory *next;
  HK_Timecode_t *tcnode;
}HK_DSN_RANGE_t;
#endif
