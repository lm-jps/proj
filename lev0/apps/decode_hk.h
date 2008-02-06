#ifndef __DECODE_HK_H
#define __DECODE_HK_H
#include "load_hk_config_files.h"
#include "packets.h"
#include "printk.h"
/*****************************************************************************
 * Filename: decode_hk.h                                                     *
 * Author: Carl Cimilluca                                                    *
 * Create Date: August, 5, 2005                                              *
 * Description: This file contains structures and link list structures       *
 *              used by the modules in decode_hk.c.                          *
 * (C) Stanford University, 2005                                             *
 * Revision History:                                                         *
 * Draft/Revision  Date           Description                                *
 * --------------  -------------  ----------------------------------------   *
 * 0.0               8/5/2005     Created                                    *
 * 0.1               10/31/2005   Added version number fields as string      *
 *                                buffer rather than as a float value        *
 * 0.2               11/04/2005   Changed HK_Keyword to HK_Keywords_Format   *
 *                                structure. Added Rasmus's HK_Keyword struct*
 *                                definition to integrate with L0 code       *
 *                                Added packets.h to file                    *
 * 0.3               01/19/2006   Added HK_DECODER_ERROR_UNKNOWN_ENV_VARIABLE*
 *                                define.                                    *
 * 0.4               05/10/2006   Updates for decode of analog and digital   *
 *                                type keywords. Added DSC_Conversion_struct *
 *                                and ALG_Conversion_struct to HK_Keyword_   *
 *                                Format structure.                          *
 *****************************************************************************/
/*****************************includes****************************************/
/*****************************defines*****************************************/
/* Maximum sizes of arrays in HMI_HK_Keyword_struct */
#define KEYWORD_NAME_SIZE                          10
#define KEYWORD_TYPE_SIZE                          10
#define TELEMETRY_MNEMONIC_SIZE                    50
#define HK_DECODER_SUCCESSFUL                      SUCCESS
#define HK_DECODER_ERROR_UNKNOWN_APID   ERROR_HK_UNKNOWN_APID
#define HK_DECODER_ERROR_NO_CONFIG_DATA            ERROR_HK_UNKNOWN_APID
#define HK_DECODER_ERROR_CANNOT_FIND_VER_NUM       ERROR_HK_CANNOT_FIND_VER_NUM
#define HK_DECODER_ERROR_CANNOT_LOAD_HK_VALUES     ERROR_HK_CANNOT_LOAD_HK_VALUES
#define HK_DECODER_ERROR_CANNOT_LOAD_ENGR_VALUES   ERROR_HK_CANNOT_LOAD_ENGR_VALUES
#define HK_DECODER_ERROR_UNKNOWN_ENV_VARIABLE      ERROR_HK_CANNOT_LOAD_ENGR_VALUES
/***************************type definitions**********************************/
/* Use to store a link list of pointers to the link list of HK Keywords Format*/
typedef struct Pointers_HK_Keywords_Format_struct 
{
  int apid;
  char  version_number[10];
  int fsn;
  int  fid;
  int number_bytes_in_packet;
  struct HK_Keywords_Format_struct  *ptr_hk_kw;
  struct Pointers_HK_Keywords_Format_struct  *next;
} Pointers_HK_Keywords_Format;
/* Use to store a link list of HK Keywords */
typedef struct HK_Keywords_Format_struct 
{
  char keyword_name[KEYWORD_NAME_SIZE]; /*Size is limited to 8 
                                         *from HMI-EGSE-FS document 
                                         */
  unsigned int keyword_value;           /*Actual value of keyword from 
                                         *hk bit stream packet
                                         */
  char type[KEYWORD_TYPE_SIZE] ;        /*Used by Decoder to parse keyword
                                         *values. Size is limited to 3 from 
                                         *HMI-EGSE-SS-FS document 
                                         */
  char telemetry_mnemonic_name[TELEMETRY_MNEMONIC_SIZE];/*Used to debug  */
  int start_byte ;                      /*Used by Decoder Parsing keyword 
                                         *values. For example: 8th byte, 
                                         *14th byte, 18th byte, etc. Used 
                                         *to debug
                                         */
  int start_bit;
  int bit_length;
  char conv_type;                       /*keyword can be conv type R,A,D */
  struct ALG_Conversion_struct  *alg;   /*analog config data to decode   */
  struct DSC_Conversion_struct  *dsc;   /*digital config data to decode  */
  struct HK_Keywords_Format_struct *next;/*Currently there are about 50
                                           keywords 
                                          */
} HK_Keywords_Format;

int decode_hk_keywords(unsigned short *ptr, int apid, HK_Keyword_t **kw_head);
void deallocate_hk_keyword(HK_Keyword_t *head);
HK_Keyword_t *copy_hk_keywords(HK_Keyword_t *head);


#endif
