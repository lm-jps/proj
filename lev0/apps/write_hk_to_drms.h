#ifndef __WRITE_HKOB_TO_DRMS_H
#define __WRITE_HK_TO_DRMS_H

#define HK_LEV0_MAX_PVN_STR           10
#define HK_LEV0_MAX_JVN_STR           10
#define HK_LEV0_MAX_FILE_NAME         300
#define HK_LEV0_MAXLINE_IN_FILE       200
#define HK_LEV0_MAX_DSNAME_STR        100
#define HK_LEV0_MAX_PROJECT_NAME      50
#define HK_LEV0_MAX_DATATYPE_NAME     50
#define HK_LEV0_MAX_LONG_KW_NAME      50
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


/* Structures */


/***************** Map Data struct *************************/
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


#endif
