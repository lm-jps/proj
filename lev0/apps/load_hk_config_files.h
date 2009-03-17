#ifndef __LOAD_HK_CONFIGS_H
#define __LOAD_HK_CONFIGS_H

/*****************************************************************************
 * Filename: load_hk_configs.h                                               *
 * Author: Carl Cimilluca                                                    *
 * Create Date: August, 5, 2005                                              *
 * Description: This file contains the structures and link list structures   *
 *              used to save the configuration formats for the Housekeeping  *
 *              Packet format from the high rate bit stream coming from      *
 *              the HMI-EGSE-SS telemetry simulator and DDS.                 *
 *****************************************************************************/
/******************************defines****************************************/
#define MAXLINE_IN_FILE         150
#define MAX_NUM_KW_LINES        3000
#define MAX_DATE_SIZE           15
#define MAX_HK_MNM              30
#define MAX_HK_KYWD             10
#define MAX_HK_VALUE_TYPE       10
#define MAX_APID_POINTERS       4095
#define MAX_DIRECTORY_NAME      200
#define MAX_FILE_NAME           300
#define ERROR_LOADING_HK_DATA   1
#define ERROR_SAVING_FILE_LINE  1
#define ERROR_SAVING_APID_LINE  1
#define SUCESSFUL               0
#define HK_SUCCESSFUL              0
#define MAX_SIZE_CHANGE_TIME    15
#define MAX_SIZE_CHANGE_DATE    15
#define MAX_CHAR_VERSION_NUMBER 10
#define MAX_NUM_ACON_LINES      1000
#define MAXLINE_ACON_IN_FILE    200
#define MAX_NUM_DCON_LINES      4000
#define MAXLINE_DCON_IN_FILE    200
#define MAX_VALUE_DSC           100
#define MAX_NUMBER_COFFICIENTS  6
#define HMI_ID_TYPE             "HMI"
#define AIA_ID_TYPE             "AIA"
#define OTHER_ID_TYPE           "SSIM"
#define MAX_PACKET_ID_TYPE      100
#define HK_MAX_TLM_NAME         200
#define HK_LR_HMI_ISP           29
#define HK_LR_HMI_SEQ           21
#define HK_LR_HMI_OBT           18
#define HK_LR_AIA_ISP           39
#define HK_LR_AIA_SEQ           46
#define HK_LR_AIA_OBT           50
#define HK_LR_SDO_ASD           129
#define HK_HSB_HMI_ISP_1        445
#define HK_HSB_HMI_ISP_2        475
#define HK_HSB_HMI_SEQ_1        451
#define HK_HSB_HMI_SEQ_2        481
#define HK_HSB_HMI_OBT_1        448
#define HK_HSB_HMI_OBT_2        478
#define HK_HSB_AIA_ISP_1        529
#define HK_HSB_AIA_ISP_2        569
#define HK_HSB_AIA_SEQ_1        536
#define HK_HSB_AIA_SEQ_2        576
#define HK_HSB_AIA_OBT_1        540
#define HK_HSB_AIA_OBT_2        580


/*************************typedef structures**********************************/
/*****************************************************************************/
/* Structure used to contain  pointer to the HK Configuration                */
/*  data based on the APID of packet.  This is a linked list .               */
/*                                                                           */
/*****************************************************************************/
typedef struct APID_Pointer_HK_Configs_struct    
{
  char apid_name[10]; /*#1:CC: added 12-15-2008*/
  int apid;                          /*Make either hex or decimal value*/
  struct HK_Config_Files_struct *ptr_hk_configs;/*Pointer to HK 
						  Configurations struct*/
  struct APID_Pointer_HK_Configs_struct *next;/* Link List's next node */
} APID_Pointer_HK_Configs;
/*****************************************************************************/
/* Structure used to contain file parameters and pointer to the              */
/* many keyword parameters each with 6 or so parameters .                    */
/* The structure is a link list of hk_config-file structures                 */
/*****************************************************************************/
typedef struct HK_Config_Files_struct  
{
  char apid_name[10];/*#2:CC: added 12-15-2008*/ 
  int apid_number;
  char packet_id_type[MAX_PACKET_ID_TYPE];
  char file_version_number[MAX_CHAR_VERSION_NUMBER];/*Value for the apid-version<version-
				                      number>.txt file */
  char parameter_version_number[MAX_CHAR_VERSION_NUMBER];/*Values in GTCIDS map file for 
				   HMI_VER_NUM_SEQ_STATUS */
  char date[MAX_DATE_SIZE];      /*Example: 2005/07/06 or 10 characters*/
  int number_bytes_used;        /* Example: 112 in apid file and stanford file*/
  struct Keyword_Parameter_struct *keywords;/*Pointer to array or vector of
					      keyword_parameter structures*/
  struct HK_Config_Files_struct  *next;
} HK_Config_Files;
                                                                                                           
/****************************************************************************/
/*   Structures used to contain values for analog and digital               */
/*   conversion values to use to set engineering values for                 */
/*   A and D type of keywords. This structure is a linked list.             */
/****************************************************************************/
typedef struct ALG_Conversion_struct
{
  int     number_of_coeffs;
  double  coeff [MAX_NUMBER_COFFICIENTS];
} ALG_Conversion;

typedef struct DSC_Conversion_struct
{
  int     low_range;
  int     high_range;
  char    dsc_value[MAX_VALUE_DSC];
  struct  DSC_Conversion_struct  *next;
} DSC_Conversion;

/****************************************************************************/
/*   Structure used to contain values of hk keyword parameters              */
/*   This structure is a linked list of keyword_parameter structures        */
/****************************************************************************/
typedef struct Keyword_Parameter_struct  {
  char telemetry_mnemonic_name[MAX_HK_MNM];/*Make 30-original name
					     from Lockheed Database and
					     STANFORD_TLM_HMI_AIA.txt files*/
  char keyword_name[MAX_HK_KYWD]; /*Make size at least 10 - 8 chararacters
				    of less parameter passed to L0P Decode
				    component */
  int start_byte;       /*L0P Decode will use this to parse packet data */
  int start_bit_number;
  int bit_length;       /*L0P Decode will use this to determine bit length
			  to parse packet data */
  char conv_type;       /* type of conversion facTtor */
  char type[MAX_HK_VALUE_TYPE];     /* Make at least 5 -Examples: UB, UL1, i
				       UL2, etc  or 3 characters*/
  struct ALG_Conversion_struct *alg_ptr;
  struct DSC_Conversion_struct *dsc_ptr;
  struct Keyword_Parameter_struct *next;
} Keyword_Parameter;
/****************************************************************************/
/*   Structure used to contain values of filename from directory            */
/****************************************************************************/
typedef struct APID_HKPFD_Files_struct  
{
  char version_number[MAX_CHAR_VERSION_NUMBER];
  int apid; 
  char apid_name[10]; 
  char directory_name[MAX_DIRECTORY_NAME];
  char filename[MAX_FILE_NAME];
  struct APID_HKPFD_Files_struct *next;
}   APID_HKPFD_Files;
/****************************************************************************/
/*   Structure used to contain values of GTCIDS file                        */
/****************************************************************************/
typedef  struct GTCIDS_Version_Number_struct   
{
  char file_version_number[MAX_CHAR_VERSION_NUMBER];
  char hmi_id_version_number[MAX_CHAR_VERSION_NUMBER];
  char aia_id_version_number[MAX_CHAR_VERSION_NUMBER];/*Values in GTCIDS map file */
  char change_date[MAX_SIZE_CHANGE_DATE];
  char change_time[MAX_SIZE_CHANGE_TIME];
  struct GTCIDS_Version_Number_struct  *next;
}   GTCIDS_Version_Number;
/****************************************************************************/
/*   Structure used to contain values of SHCIDS file                        */
/****************************************************************************/
typedef  struct SHCIDS_Version_Number_struct   
{
  char file_version_number[MAX_CHAR_VERSION_NUMBER];
  int apid;
  char change_date[MAX_SIZE_CHANGE_DATE];
  int  date; //change_date's yyyymmdd as integer 
  char change_time[MAX_SIZE_CHANGE_TIME];
  struct SHCIDS_Version_Number_struct  *next;
}   SHCIDS_Version_Number;
/*************************function prototypes*********************************/

APID_HKPFD_Files* read_all_hk_config_files(int apid, char file_version_number[],char packet_version_number[]);

void load_config_data( APID_HKPFD_Files* , APID_Pointer_HK_Configs* );
int save_hdpf_new_formats(FILE* file_ptr,APID_Pointer_HK_Configs *ptr);
int update_hdpf_new_formats(FILE* file_ptr, HK_Config_Files *ptr_config_node,
			    HK_Config_Files *top_ptr_config_node); 
int load_hdpf_keyword_lines(char keyword_lines[][], int i, HK_Config_Files *ptr_config_node);
int load_hdpf_dsc_lines(char dsc_lines[][], int j, HK_Config_Files *ptr_config_node);  
int load_hdpf_alg_lines(char alg_lines[][], int k, HK_Config_Files *ptr_config_node);
DSC_Conversion* create_hdpf_dsc_nodes(Keyword_Parameter *ptr_hk_keyword);
Keyword_Parameter* create_hdpf_keyword_nodes(HK_Config_Files *ptr_hk_config_node,int number_of_lines);  
GTCIDS_Version_Number * read_gtcids_hk_file();
SHCIDS_Version_Number * read_shcids_hk_file();
int load_all_apids_hk_configs(int apid, char version_number[], char pkt_date[]);
void deallocate_apid_ptr_hk_config_nodes(void);
void deallocate_GTCIDS_Version_Number(void);
void deallocate_SHCIDS_Version_Number(void);
HK_Config_Files* check_packet_version_number( HK_Config_Files *top_ptr_to_configs,
						char *version_number );
HK_Config_Files* check_file_version_number( HK_Config_Files *top_ptr_to_configs,
						char *version_number );
HK_Config_Files * reread_all_files(APID_Pointer_HK_Configs *apid_ptr, char version_number[], char pkt_date[]);
int check_free_configs_flag(void);
int (*printkerr)(const char *fmt, ...);


#ifndef  LOAD_HK_C
extern APID_Pointer_HK_Configs *gb_top_apid_ptr_hk_configs;
#endif

#endif
