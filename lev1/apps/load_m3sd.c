/*############################################################################
# Name:        load_m3sd.c - load mean max min and sd in series              #
# Description: Load Minimum, Maximum, Mean and Standard Deviation keyword    #
#              values over interval of time. Get keyword value by sending    #
#              sending hk packets to functions in decode_hk.c file. Keywords #
#              are decoded using the  same functions that decodes hk keywords#
#              on high speed bus. load_m3sd.c loops through all hk packets & #
#              sends data to be decoded to decode_hk.c and then gets all     #
#              keyword values defined in instruction file for interval of    #
#              packet time and then finds  minimum, maximun, mean and        #
#              standard deviation and write values to drms series. The       #
#              keywords and interval of time, and data series are defined in #
#              the instruction file which has a predetermined format.        #
#              The "in" file is a hk dayfile.                                #
# Execution:   load_m3sd  in=<day-file> isf=<instruction-file>               #
#              (1)To load minimum,maximum,mean and standard deviation value: #
#                    load_m3sd in=<day-file> isf=<instruction file>          #
# Example  :  load_m3sd in=/tmp20/20070202.0x001d  isf=./temp_instr_file.txt #
# Limitation:  Setup required for environment variables in file              #
#              SOURCE_HK_LM3S. Currently the code read packet size           #
#              limit of 1000 bytes (HKLMS_MAX_PKT_SIZE ). Need to setup local#
#              version of SOURCE_HK_LM3S file since file checked in          #
#              are production version of file.                               #
# Author:      Carl                                                          #
# Date:        September,09, 2009                                            #
############################################################################*/

/**
   @defgroup load_m3sd load_m3sd 
   @ingroup su_lev1

   @brief Loads minimum, maximum, mean and standard deviation for keywords over a interval of time into data series.

   @par Synopsis:
   @code
   load_m3sd -h
   load_m3sd  in=<day filename> isf=<instruction file> 
   @endcode

   Load minimum, maximum, mean and standard deviation(m3sd) values by first sending packet
   to functions in decode_hk.c file to get keyword values. Keywords are decoded using
   the the same functions that decodes hk keywords on high speed bus. This load_m3sd executable
   loops through all hk packets and sends data to be decoded to decode_hk.c's functions
   which returns a structure containing keyword names, keyword values, keyword type
   and keyword variable. load_m3sd executable then gets keyword values per interval of
   time and calculates the minimum, maximum, mean and standard deviation for keywords. 
   Finally, load_m3sd executable writes to user defined keywords names with MIN, MAX, MEAN, and SD
   suffixed to each keyword name,  the values to a DRMS data series. DRMS data series is based on value
   defined in the instruction file. The interval of time and the keywords to get values for
   are contained in the instruction file. The keywords in m3sd data series are loaded as floats.
   A requirement is to have a data series created before running this executable. This can be done
   by running script located in lev1/scripts called cm3sd_jsd_file.pl using the instruction file.
   This script will create the jsd file. After creating jsd file, then use create_series to implement 
   data series in DRMS. An example on running this script is as follows: 
   cm3sd_jsd_file.pl isf=/home/carl/cvs/TBL_JSOC/lev1/instruction_file/su_carl/hmitest1200_thermal_template.txt


   The in parameter is a mandatory argument which should contain the directory and 
   filename of the input dayfile. Currently only one dayfile is allowed for the in 
   parameter.

   The isf parameter is mandatory argument. This is the full path to the instruction file. The
   instruction contains data series name, owner and author. The instruction file contains
   interval of time to use to calculate the minimum, maximum, mean and standard deviation for
   keywords. The current interval ranges tested are 60(seconds) to 7200(seconds). The
   instruction file contains the list of keywords. Each keyword is represented by
   a line containing the packet apid from which the keyword is contained. The line contains the
   Keyword:,apid value in decimal,long keyword name, and user defined name. Currently the user 
   defined name should be 8 chararcters. User's can use short keyword name for user defined names. 
   It is required to use the proper long keyword name that is contained in the
   HK Configuration file and STANFORD file. The line contains a user created keyword name.
   The user created names will have a MAX, MIN, MEAN or SD and suffixed to names(i.e.,TEMP1_MIN,
   TEMP1_MAX, TEMP1_MEAN, TEMP1_SD). View example instruction files at directory in CVS at  
   /home/production/cvs/TBL_JSOC/lev1/instruction_file/su_carl.

   @par Flags:
   @c -h: Shows usage message.
   @par

   @param in The full directoy path and file name to the input dayfile(required field).

   @param isf The full directory path and file name to instruction file(required field).

 
   @par Example of running script:
   @code
   load_m3sd  in=/home/production/dfile/20080918.0x0013 isf=/home/production/inst_file/instruction_file.txt
   @endcode

   @par Example of running script debug conditional compile flag turned on:
   @code
   load_m3sd  in=/home/production/df/20080918.0x0013 isf=/home/production/inst_file/instr_file.txt > Debug-Log
   @endcode

   @par Example of running help:
   @code
   load_m3sd -h
   @endcode

*/
/* Defined constants */
/******************** defines ***********************************************/
#define HKLMS_LONG_KEYWORD_NAME_SIZE   100
#define HKLMS_MAX_AU_OW_NAME           50
#define HKLMS_MAX_DSNAME_STR           100
#define HKLMS_MAX_FILE_NAME            100
#define HKLMS_MAX_INSTR_FILENAME       100
#define HKLMS_MAX_PACKET_DESCRIPTION   100
#define HKLMS_MAX_PKT_SIZE             1000
#define HKLMS_MAX_PVN_SIZE             50
#define HKLMS_PACKET_TIME_STR          100
#define HKLMS_READ_ARRAY_SIZE          (25000001)
#define HKLMS_SHORT_KEYWORD_NAME_SIZE  100
/*#define ENVFILE  "/home/production/cvs/JSOC/proj/lev1/apps/SOURCE_ENV_FOR_LOAD_M3SD" */
/*#define ENVFILE    "/home3/carl/cvs/JSOC/proj/lev1/apps/SOURCE_ENV_FOR_LOAD_M3SD"*/
#define ENVFILE  "/home/production/cvs/JSOC/proj/lev1/apps/SOURCE_ENV_FOR_LOAD_M3SD"

/******************** includes ******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include "drms.h"
#include "drms_names.h"
#include "jsoc_main.h"
#include "load_hk_config_files.h"
#include "write_hk_to_drms.h"
#include "packets.h"
#include "decode_hk.h"
#include "printk.h"

/*  @{ */
/************* modules definitions **************************************/
ModuleArgs_t module_args[] =
{
  {ARG_STRING, "in", "Not Specified", "full path to day file"},
  {ARG_STRING, "isf", "Not Specified", "full path to instruction file"},
  {ARG_END}
};
ModuleArgs_t   *ggModArgs=module_args;
char* module_name = "load_m3sd";

/********************* structures   *************************************/

/* structures holding instruction file */
typedef struct Instruction_File_Keywords_struct_t  Instruction_File_Keywords_t;
struct Instruction_File_Keywords_struct_t
{
  char kw_longname[HKLMS_LONG_KEYWORD_NAME_SIZE];
  char kw_shortname[HKLMS_SHORT_KEYWORD_NAME_SIZE];
  int dec_apid;
  int hex_apid;
  Instruction_File_Keywords_t *next;
};


typedef struct Instruction_File_Data_struct_t
{
  char inst_filename[HKLMS_MAX_INSTR_FILENAME];
  char series_name[HKLMS_MAX_DSNAME_STR];
  char author[HKLMS_MAX_AU_OW_NAME];
  char owner[HKLMS_MAX_AU_OW_NAME];
  char packet_description[HKLMS_MAX_PACKET_DESCRIPTION];
  int interval;
  Instruction_File_Keywords_t *kw_if;
} Instruction_File_Data_t;

/* structure holding min,max,mean and std dev(m3sd) values to write to drms */
typedef struct HK_Keyword_M3SD_Values_struct_t
{
  char kw_longname[HKLMS_LONG_KEYWORD_NAME_SIZE];
  char kw_shortname[HKLMS_SHORT_KEYWORD_NAME_SIZE];
  float mean;
  float min;
  float max;
  float stdev;
  struct  HK_Keyword_M3SD_Values_struct_t *next;
} HK_Keyword_M3SD_Values_t;

typedef struct HK_Keyword_M3SD_Data_struct_t
{
  TIME start_pkt_time;
  int number_points;
  int dec_apid;
  int hex_apid;
  struct  HK_Keyword_M3SD_Data_struct_t *next;
  struct  HK_Keyword_M3SD_Values_struct_t *m3sd_values;
} HK_Keyword_M3SD_Data_t;

/* struction hold temp data per interval of time which is used to calculate m3sd */
typedef struct HK_KW_Data_Values_struct_t
{
  TIME pkt_time;
  KW_Type_t       eng_type;  /* Engineering value type. */
  KW_Type_Value_t eng_value; /* Engineering value. */
  struct HK_KW_Data_Values_struct_t  *next;
} HK_KW_Data_Values_t;

typedef struct HK_KW_Data_struct_t
{
  TIME start_pkt_time;
  char  longname[MAX_KEYWORD_NAME_SIZE];
  char  shortname[MAX_FITS_NAME_SIZE];
  struct HK_KW_Data_struct_t  *next;
  struct HK_KW_Data_Values_struct_t  *kw_values;
} HK_KW_Data_t;

/******************* function prototypes  *******************************/
static void  set_env_variables();
static TIME  SDO_to_DRMS_time(int sdo_s, int sdo_ss);
void   get_next_pkt_time(TIME p_time, int intval, TIME *start, TIME *end);
int    get_minute_from_pkttime(double tc_sec);
int   get_seconds_from_pkttime(double tc_sec);
static int   get_packet_time_for_df(HK_Keyword_t *hk, TIME *ptime);
HK_KW_Data_t * save_kw_data(TIME start_pkt_time, Instruction_File_Data_t *ifs, HK_Keyword_t *kw );
void  save_kw_data_values(HK_KW_Data_t *top_kw_data_ptr, TIME pkt_time, HK_Keyword_t *kw );
void  free_all_kw_data(HK_KW_Data_t *top_kw_data_ptr);
void  save_m3sd_data(HK_KW_Data_t *top_kw_data_ptr, HK_Keyword_M3SD_Data_t **kw_m3sd_head);
void  save_m3sd_values(HK_KW_Data_t *top_kw_data_ptr, HK_Keyword_M3SD_Data_t **kw_data_head);
float get_min_value(HK_KW_Data_Values_t *ptr);
float get_max_value(HK_KW_Data_Values_t *ptr);
float get_mean_value(HK_KW_Data_Values_t *ptr);
float get_stdev_value(HK_KW_Data_Values_t *ptr);
void  write_m3sd_to_drms(HK_Keyword_M3SD_Data_t *top_m3sd_data_ptr, Instruction_File_Data_t *ifp);
int   get_number_points(HK_KW_Data_Values_t *ptr);
Instruction_File_Data_t * read_isf_data(char *inf);
void  get_keyword_values(unsigned char *read_in_buffer, Instruction_File_Data_t *ifp);
void check_status_drms_set(int status, char *kwn);
int get_day_from_pkttime(double p_time);
int get_month_from_pkttime(double p_time);
int get_yr_from_pkttime(double p_time);
void my_usage (void);

/********************* extern functions  *********************************/
extern int  get_hour_from_pkttime(double p_time);
extern SHCIDS_Version_Number *global_shcids_vn;
extern char * find_file_version_number(GTCIDS_Version_Number *top,char f_version_number[MAX_CHAR_VERSION_NUMBER], int apid);
extern char * find_fvn_from_shcids(SHCIDS_Version_Number *top,char pkt_date[MAX_SIZE_PKT_DATE],int apid);
extern double  get_packet_time(unsigned short *word_ptr);
extern int DoIt(void);
extern int nice_intro (void);
extern int check_for_sdo_apid(int apid);
extern void sprint_time (char *at, TIME t, char *zone, int precision);
extern char *get_data_packet_name(int apid ) ;
extern char *get_lookup_filename(int apid);
extern int check_hk_record_exists(char* ds_name, HK_Keyword_t *kw, int apid);
extern int check_hk_record_within_time_range( HK_Keyword_t *kw);

/********************* extern globals  *********************************/
extern GTCIDS_Version_Number *global_gtcids_vn;
extern GTCIDS_Version_Number *global_gtcids_vn;

/* @} */



/*************************************************************************
 * my_usage                                                              *
 * FUNCTION: myusage(void)                                               *
 * DESCRIPTION: Use function to display usage                            *
 *************************************************************************/
void my_usage (void)
{
  
  printf ("Usage:\nload_m3sd  [-h] "
    "in=<day filename> "
    "isf=<instruction file>  \n"
    "  details are:\n"
    "  -h: help - show this message then exit(optional field)\n"
    "  in=<day file name> -use full path to day file(required field)\n"
    "  isf=<instruction file name> -use full path to instruction file(required field)\n"
    "  Need data series already created for this to program to load keywords in data series.\n"
    "  Setup interval time in this file where range of intervals tested are 60 seconds to 7200 seconds.\n"
    "  Setup keywords to calculate mean, max, min and standard deviation values in this instruction file.\n"
    "  Setup seriesname in this instruction file for writing keyword min,max,mean and standard deviation values.\n"
    "  Example of running:\n"
    "  load_m3sd  in=/home/carl/cvs/myprod/JSOC/proj/lev1_hmi/apps/20080918.0x0013 isf=/home/carl/cvs/TBL_JSOC/lev1/instruction_file/hmitest1200_thermal_template.txt\n" );

}


/*************************************************************************
 * Nice Intro                                                            *
 * FUNCTION: nice_intro(void)                                            *
 * DESCRIPTION: Use function to display help                             *
 *************************************************************************/
int nice_intro (void)
{
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage == 1)
  {
    my_usage();
    return(1);
  }
  return (0);
}



/*************************************************************************
 * DoIT                                                                  *
 * FUNCTION: DoIt(void                                                   *
 * DESCRIPTION: Use function to create jsoc module                       *
 *************************************************************************/
int DoIt(void)
{
  /* variables */
  FILE *file_ptr;
  Instruction_File_Data_t *instr_file_ptr;
  Instruction_File_Data_t * read_isf_data(char *inf);
  char *hk_df_fn;
  char  hk_directory_filename[HKLMS_MAX_FILE_NAME];
  unsigned char *ptr_read_in_buffer;
  unsigned long int i;

  /* set environment variables */
  set_env_variables();

  /* parameter initialization */
  hk_df_fn= hk_directory_filename;

  /* Get command line arguments */
  char *in = cmdparams_get_str (&cmdparams, "in", NULL);
  char *isf = cmdparams_get_str (&cmdparams, "isf", NULL);

  /* check arguments used */
  if (nice_intro ()) return (0);

  /* check if entered day file name */
  if (in == NULL) 
  {
    printkerr("ERROR at %s, line %d: Need to enter day file name."
              "Exiting program.\nExample format for -in- file: "
              "$HOME/data/dayfiles/20060719.0x000f\n",
               __FILE__,__LINE__);
    my_usage();
    return (0);
  }

  /* get isf filename and load instruction file structure */
  instr_file_ptr= read_isf_data(isf);
  if (instr_file_ptr == NULL) 
  {
    printkerr("ERROR at %s, line %d: Note need to enter correct path "
              "and filename for instruction file <%s>. Example format for -isf- file: "
              "isf=/home/data/instruction_file1.txt. Exiting this program.\n ",
               __FILE__,__LINE__, isf);
    my_usage();
    return (0);
  }

 /*PRINT RESULTS FOR INSTRUCTION FILES DATA STRUCTURE */
#ifdef DEBUG_LM3S
 Instruction_File_Keywords_t *ptr_if_kws;
 printf("DoIt:Print Structure\n");
 printf("DoIt:interval:%d des:%s \n",instr_file_ptr->interval, instr_file_ptr->packet_description);
 printf("DoIt:owner:%s author:%s \n",instr_file_ptr->author,instr_file_ptr->owner);
 printf("DoIt:seriesname:%s  \n",instr_file_ptr->series_name);
 printf("DoIt:templatename:%s  \n",instr_file_ptr->inst_filename);
 for(ptr_if_kws =instr_file_ptr->kw_if;ptr_if_kws;ptr_if_kws=ptr_if_kws->next )
 {
   printf("DoIt: ptr_if_kws->kw_longname=<%d>\n",ptr_if_kws->dec_apid );
   printf("DoIt: ptr_if_kws->kw_longname=<%s>\n",ptr_if_kws->kw_longname );
   printf("DoIt: ptr_if_kws->kw_shortname=<%s>\n\n",ptr_if_kws->kw_shortname );
 }
#endif

  /* get  in filename and open hk dayfile*/
  strcpy(hk_df_fn, in);
  strcat(hk_df_fn, "\0");
  file_ptr=fopen(hk_df_fn,"r");
  if (!file_ptr)
  {
    printkerr("ERROR at %s, line %d: Please check filename and directory is correct. "
              " Could not get -in- directory and filename: "
              "<%s>. Example format for -in- file: in=/home/rock/20080909.0x013. Exiting execution.\n", 
              __FILE__,__LINE__, hk_df_fn);
    return(0);
  }

  /* malloc memory for holding data in file in memory */
   ptr_read_in_buffer = (unsigned char *) malloc(sizeof(unsigned char) * HKLMS_READ_ARRAY_SIZE);
 
  /*read lines in file into buffer representing packet data in file*/
  for(i=0; i < HKLMS_READ_ARRAY_SIZE;i++) ptr_read_in_buffer[i]=0; ;
  for(i = 0 ; fread(ptr_read_in_buffer + i,1,1,file_ptr) ; i++) 
  {  
    ;/* do nothing*/
    if( i == HKLMS_READ_ARRAY_SIZE - 1)
    {
      printkerr("ERROR at %s, line %d: Array for reading dayfile is too small. :"
                "<%lu>\n", __FILE__,__LINE__, i + 1);
      printkerr("Break up large dayfile using dd command. :   dd if<large-file> "
                " of=<small-file-1> bs=<packet_size + 7> count=<i.e., 1000> skip=<0,1,etc>\n");
      return (0);
    }
  }
  /* set last value in buffer to null */
  *(ptr_read_in_buffer + i) = '\0';

  /* close dayfile */
  fclose(file_ptr);

  /* get keyword values from dayfile for keywords outlined in instruction file,calculate m3sd, write m3sd to drms */
  (void)get_keyword_values(ptr_read_in_buffer, instr_file_ptr);
  return 0;  
}


/*************************************************/
/* Get Packet Time reused from decode_dayfile    */
/*************************************************/
static int get_packet_time_for_df(HK_Keyword_t *hk,  TIME *ptime)
{
  /* variables */
  HK_Keyword_t *t_hk;
  int sec;
  int subsec;
  int SEC_FOUND_FLAG=0;
  int SUBSEC_FOUND_FLAG=0;

  /* init variables */
  t_hk= hk;

  /*loop until get TIMECODE Seconds and Subseconds values */
  while (t_hk && ( !SEC_FOUND_FLAG || !SUBSEC_FOUND_FLAG ))
  {
    if(strstr(t_hk->name,"TIMECODE_SECONDS")) 
    {
      /*set found flag*/
      SEC_FOUND_FLAG=1;
      /* create time string based on some hard coded parameters for now */
      sec = t_hk->eng_value.uint32_val;
    }
    if(strstr(t_hk->name,"TIMECODE_SUBSECS")) 
    {
      /*set found flag*/
      SUBSEC_FOUND_FLAG=1;
      /* create time string based on some hard coded parameters for now */
      subsec = t_hk->eng_value.uint32_val;
    }
    t_hk= t_hk->next;
  }
  /* check if found TIMECODE_SECONDS and TIMECODE_SUBSECS */
  if (!SEC_FOUND_FLAG)
  {
    printkerr("ERROR at %s, line %d: Did not find TIMECODE_SECONDS value for"
              "calculating the PACKET_TIME keyword and index. Returning error"
              "status.\n",__FILE__,__LINE__);
    return 0;
  }
  else
  {
    if (!SUBSEC_FOUND_FLAG)
    {
      printkerr("ERROR at %s, line %d: Did not find TIMECODE_SUBSECS value for"
                "calculating the PACKET_TIME keyword and index. Returning error"
                "status, but setting time using seconds only.\n",__FILE__,__LINE__);
      *ptime =SDO_to_DRMS_time(sec, 0);
      return 0;
    }
    else
    { 
      int shifted_ss=(subsec >> 16) & 0xFFFF;
      *ptime =SDO_to_DRMS_time(sec, shifted_ss);
      return 1;
    }
  }
}

/****************************************************/
/*  SDO_to_DRMS_time                                */
/****************************************************/
static TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss) 
{
static int firstcall = 1;
static TIME sdo_epoch;
if (firstcall)
  { 
  firstcall = 0;
  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  } 
return(sdo_epoch + (TIME)sdo_s + (TIME)(sdo_ss)/65536.0);
}

/****************************************************/
/* Get Keyword Values                               */
/****************************************************/
void get_keyword_values(unsigned char *read_in_buffer, Instruction_File_Data_t *ifp)
{ 
  CCSDS_Packet_t ccsds;
  HK_Keyword_M3SD_Data_t *top_m3sd_data_ptr;
  HK_KW_Data_t *top_kw_data_ptr;
  HK_Keyword_t *kw_head,*kw;
  TIME pkt_time;
  TIME next_pkt_time;
  TIME start,end;
  TIME *st_ptr, *ed_ptr;
  char file_version_number[HKLMS_MAX_PVN_SIZE];
  char packet_version_number[HKLMS_MAX_PVN_SIZE];
  char pkt_date[MAX_SIZE_PKT_DATE]; //ascii time
  char *ptr_fvn;
  int apid;
  int factor;
  int interval_count; 
  int interval;
  int i,j,k,y,s;
  int packet_length;
  unsigned char hk_pkt_buffer[HKLMS_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKLMS_MAX_PKT_SIZE];
  unsigned short *word_ptr;
#ifdef DEBUG_LM3S
  char at[HKLMS_PACKET_TIME_STR];
  char atn[HKLMS_PACKET_TIME_STR];
#endif

  /* initialize variables */
  top_m3sd_data_ptr=NULL;
  top_kw_data_ptr=NULL;
  interval_count=0; 
  st_ptr= &start;
  ed_ptr= &end;

  /* get interval of time to use in seconds -where 600s = 10minutes */
  interval=ifp->interval;

  /* set next packet to for first packet handling to 0.0 */
  next_pkt_time=0.0;

  /* go thru each packet and save keywords to DRMS */
  for(k=0,factor=0,packet_length=0;  ; k++)
  {
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;

    /* Check if at end of all pkts */
    if (*(read_in_buffer+5+factor) == '\0')
    {
      /* at end of all packets but need to process last interval */

      /* At end of block of interval values, calculate m3sd values based on interval block & save in m3sd structures*/
      if( !top_m3sd_data_ptr)
      {
        /* first time thru, creates all HK_Keyword_M3SD_Data_t values*/
        (void)save_m3sd_data( top_kw_data_ptr,&top_m3sd_data_ptr);

        /* calculate and add HK_Keyword_M3SD_Values_t values */
        (void)save_m3sd_values(top_kw_data_ptr, &top_m3sd_data_ptr);
      }
      else
      {
        /* save m3sd data for second interval time thru nth time, just add to HK_Keyword_M3SD_Values_t values*/
        (void)save_m3sd_data( top_kw_data_ptr,&top_m3sd_data_ptr);
        (void)save_m3sd_values(top_kw_data_ptr, &top_m3sd_data_ptr);
      }

      /* now free HK_KW_Data_Values to start new block of interval values */
      free_all_kw_data(top_kw_data_ptr);

      /* since completed processing packets and added data to M3S structure, break and write data to series */
      break;
    }

    /* get packet lenght */
    packet_length=  *(read_in_buffer+5+factor);

    /* get apid */
    /* set 0th to 7th bits */
    apid =  (unsigned  short int)( (*(read_in_buffer+1+factor)) & 0x00FF );
    /* set 8th to 15th bits */
    apid |= (unsigned  short int)( (*(read_in_buffer+0+factor) << 8) & 0xFF00 );
    apid &= 0x07FF;

    /* get packet version number */
    /* if sdo-hk type apid,packet version number not used,else get packet version number for hmi or aia hk packets */
    if(check_for_sdo_apid(apid))
    {
       sprintf(packet_version_number,"%s","not applicable");
    }
    else
    {
      /*check if packet version is 0.0 -skip-print warning*/
      if (( *(read_in_buffer+14+factor) == 0) && (*(read_in_buffer+15+factor) == 0))
      {
        printkerr("Warning at %s, line %d: Getting 0.0 for packet version number."
                  "in packet data. Skip processing this packet, don't load to "
                  "DRMS and get next packet", __FILE__,__LINE__);
        continue;
      }
      /* check for version number set to zero and print warning messages. */
      if ( *(read_in_buffer+14+factor) == 0 && *(read_in_buffer+15+factor))
      {
        printkerr("Warning at %s, line %d: Getting 0 for whole number(i.e.,0.1) for "
                  "packet version number in packet data. Skip processing this packet,"
                  "don't load to DRMS and get next packet", __FILE__,__LINE__);
        continue;
      }
      /* if passed two tests above, then set packet version number */
      sprintf(packet_version_number,"%03d.%03d",*(read_in_buffer+14+factor), *(read_in_buffer+15+factor));
      strcat(packet_version_number,"\0");
    }

     /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKLMS_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
    /* set buffer to values in packet and set packet version number to value
    /* in packet or some default value */
    for (i =0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
        /* set values in array */
        hk_pkt_buffer[j++]= *(read_in_buffer+i);
    } /* end for loop setting values for one extracted packet in buffer */


    /* put in format for decode_hk to read by adjusting values  */
    /* in buffer from unsigned char to unsigned short           */
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      s_hk_pkt_buffer[y] = (unsigned short)(hk_pkt_buffer[i + 1] << 8  & 0xFF00 );
      s_hk_pkt_buffer[y] = (unsigned short)((hk_pkt_buffer[i] & 0x00FF) + s_hk_pkt_buffer[y]) ;
    }
    /* send hk_pkt_buffer to decoder function */
    word_ptr = s_hk_pkt_buffer;
    s = decode_hk_keywords(word_ptr, apid,  &ccsds.keywords);
    if (s) 
    {
      printkerr("ERROR at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n\n", __FILE__,__LINE__, s);
      printf("ERROR at %s, line %d: decode_hk_keyword function returned error status number <%d>\n\n", __FILE__,__LINE__, s);
      continue;
    }
    kw_head= ccsds.keywords;

    /*find file version if apid is sdo hk type */
    if(check_for_sdo_apid(apid))
    {
      /* after load of shcids file during call to decode_hk_keywords -
         get packet time to lookup sdo hk config file version to use */
      (void)sprint_time (pkt_date, get_packet_time(word_ptr), "TAI", 0);
      /* get file version number for ADP apid */
      ptr_fvn=find_fvn_from_shcids(global_shcids_vn, pkt_date, apid);
    }
    else
    {
      /* get fvn for hk apids */
      ptr_fvn=find_file_version_number(global_gtcids_vn, packet_version_number,apid); 
    }
    /* set file version number  */
    strcpy(file_version_number,ptr_fvn);

    /*************************/
    /* check for packet time */
    /*************************/
    /* set kw ptr to top of kw link list */
    kw=kw_head;

    /* get packet time of kw packet */ 
    (void)get_packet_time_for_df(kw, &pkt_time);

#ifdef DEBUG_LM3S
    (void)sprint_time (at,  pkt_time, "UTC", 0);
    strcat(at,"\0");
    (void)sprint_time (atn, next_pkt_time, "UTC", 0);
    strcat(atn,"\0");
#endif

    /* check if next pkt time is null then doing first set of packets over the interval of time */
    /* save kw data to structure, define head or top node and then save kw values to structure */
    if(!(int)next_pkt_time)
    {
#ifdef DEBUG_LM3S
      printf("get_keyword_values:interval case1:pkt_time is      %-12.12f & %-s\n",pkt_time,at);
      printf("get_keyword_values:interval case1:next_pkt_time is %-12.12f & %-s\n",next_pkt_time,atn);
#endif
      /* get start and end pkt time */
      (void)get_next_pkt_time(pkt_time, interval, st_ptr,  ed_ptr); 

      /* set next pkt time variable to end ptr time to trigger next block interval */
      next_pkt_time= *ed_ptr;

      /* use interval count to debug number of interval processing per dayfile */
      interval_count++;

#ifdef DEBUG_LM3S
      (void)sprint_time (atn, next_pkt_time, "UTC", 0);
      strcat(atn,"\0");
      printf("get_keyword_values:interval case1:new next_pkt_time is %-12.12f & %-s\n",next_pkt_time,atn);
#endif

      /* save keyword data for this packet to HK_Data_t nodes- i.e.,shortname,longname,start time,etc */
      top_kw_data_ptr=save_kw_data(*st_ptr, ifp, kw );

      /* save keyword data values for this packet to HK_Data_Values_t nodes- i.e.,pkt time,value,value type,etc */
      (void) save_kw_data_values(top_kw_data_ptr, pkt_time, kw );
    }
    else if (pkt_time >  next_pkt_time) 
    {

#ifdef DEBUG_LM3S
      printf("get_keyword_values:interval case2:k packet:%d  pkt_time is %-12.12f & %-s\n",k,pkt_time,at);
      printf("get_keyword_values:interval case2:next_pkt_time is %-12.12f & %-s\n",next_pkt_time, atn);
      printf("get_keyword_values:interval case2:at end of interval block, so calculate and save m3sd data.\n");
#endif

      /* At end of block of interval values, calculate m3sd values based on interval block & save in m3sd structures*/
      if( !top_m3sd_data_ptr)
      {
        /* first time thru, creates all HK_Keyword_M3SD_Data_t values*/
        (void)save_m3sd_data( top_kw_data_ptr,&top_m3sd_data_ptr);

        /* calculate and add HK_Keyword_M3SD_Values_t values */
        (void)save_m3sd_values(top_kw_data_ptr, &top_m3sd_data_ptr);
      }
      else
      {
        /* save m3sd data for second interval time thru nth time, just add to HK_Keyword_M3SD_Values_t values*/
        (void)save_m3sd_data( top_kw_data_ptr,&top_m3sd_data_ptr);
        (void)save_m3sd_values(top_kw_data_ptr, &top_m3sd_data_ptr);
      }

      /* now free HK_KW_Data_Values to start new block of interval values */
      free_all_kw_data(top_kw_data_ptr);

      /* now do next interval of packets, get start and end pkt time */
      (void)get_next_pkt_time(pkt_time, interval, st_ptr,  ed_ptr); 

      /* set next pkt time variable to end ptr time to trigger next block interval */
      next_pkt_time= *ed_ptr;

#ifdef DEBUG_LM3S
      (void)sprint_time (atn, next_pkt_time, "UTC", 0);
      strcat(atn,"\0");
      printf("get_keyword_values:interval case2:new next_pkt_time is %-12.12f & %-s\n",next_pkt_time,atn);
#endif

      /* for debug information count number of intervals doing per dayfile. */
      interval_count++;

      /* save keyword data to HK_Data_t nodes- i.e.,longname,start time,etc */
      top_kw_data_ptr=save_kw_data(*st_ptr, ifp, kw );

      /* save keyword data values to HK_Data_Values_t nodes- i.e.,pkt time,value,value type,etc */
      (void) save_kw_data_values(top_kw_data_ptr, pkt_time, kw );
    }
    else
    {
#ifdef DEBUG_LM3S
      printf("get_keyword_values:interval case3:k packet:%d pkt_time is %-12.12f & %-s\n",k,pkt_time,at);
      printf("get_keyword_values:interval case3:next_pkt_time is %-12.12f & %-s\n",next_pkt_time,atn);
#endif

      /* save kw values to HK_KW_Data_Values for all kw's in Instruction_File_keyword_t struct*/
      /* save keyword data values to HK_Data_Values_t nodes- i.e.,pkt time,value,value type,etc */
      (void) save_kw_data_values(top_kw_data_ptr, pkt_time, kw );
    }
  

  } /* big for loop throught next packet and do all packets in file*/

  /* ok looped thru all packets and found packet intervals to calculate the m3sd values/
  /* write to drms  what's in HK_Keyword_M3SD_Data_t link list */
#ifdef DEBUG_LM3S
  HK_Keyword_M3SD_Values_t *ptr_m3sd_values;
  HK_Keyword_M3SD_Data_t *ptr_m3sd_data;
  for ( ptr_m3sd_data=top_m3sd_data_ptr;ptr_m3sd_data; ptr_m3sd_data=ptr_m3sd_data->next)
  {
    printf("HK_Keyword_M3SD_Data_t node start pkt time<%f\n",ptr_m3sd_data->start_pkt_time);
    for(ptr_m3sd_values = ptr_m3sd_data->m3sd_values;ptr_m3sd_values;ptr_m3sd_values=ptr_m3sd_values->next)
    {
      printf("get_keyword_values:HK_Keyword_M3SD_Values_t node longname<%s>\n",ptr_m3sd_values->kw_longname);
      printf("get_keyword_values:HK_Keyword_M3SD_Values_t node shortname<%s>\n",ptr_m3sd_values->kw_shortname);
      printf("get_keyword_values:HK_Keyword_M3SD_Values_t node MIN:<%f>\n",ptr_m3sd_values->min);
      printf("get_keyword_values:HK_Keyword_M3SD_Values_t node MAX:<%f>\n",ptr_m3sd_values->max);
      printf("get_keyword_values:HK_Keyword_M3SD_Values_t node MEAN:<%f>\n",ptr_m3sd_values->mean);
      printf("get_keyword_values:HK_Keyword_M3SD_Values_t node SD:<%f>\n",ptr_m3sd_values->stdev);
    }
  }
  printf("Interval looped thru: %d\n",interval_count);
#endif

  /* write min,max,mean, and sd values in link list structure to drms */
  (void)write_m3sd_to_drms(top_m3sd_data_ptr, ifp);

}/*end of get_keyword_values() */



/*******************************************************************/
/* read_isf_data: reads data in instruction file or template file. */
/*******************************************************************/
Instruction_File_Data_t * read_isf_data(char *inf)
{
    void remove_blanks( char string[HKLMS_SHORT_KEYWORD_NAME_SIZE], char *new_str);
    Instruction_File_Data_t *ptr_isf;
    Instruction_File_Keywords_t *top_if_kws,*ptr_if_kws;
    FILE *file_ptr;
    char line[MAXLINE_IN_FILE];
    char kwlong[HKLMS_LONG_KEYWORD_NAME_SIZE];
    char kwshort[HKLMS_SHORT_KEYWORD_NAME_SIZE];
    char new_string[HKLMS_LONG_KEYWORD_NAME_SIZE];
    int exit_flag=0;

    /* get directory and file name */
    if(!inf)
    {
      return ((Instruction_File_Data_t *)(NULL));
    }

    /* malloc Instruction_File_Data_t structure */
    ptr_isf=(Instruction_File_Data_t*)malloc (sizeof (Instruction_File_Data_t));
    ptr_isf->kw_if=NULL; //set keyword link list to null

    /*open file & read  data in structure*/
    file_ptr = fopen(inf, "r");
    if(!file_ptr)
   {
     /*printkerr("ERROR at %s, line %d: Couldn't open instruction file <%s>.\n", __FILE__,__LINE__, inf);*/
     return NULL;
   }
   top_if_kws=NULL;

   while( fgets(line, MAXLINE_IN_FILE, file_ptr) != NULL )
   {
     if(line[0] == '#')
     {
       continue; /* skip comments */
     }
     else if (!strncmp(line, "templatename:", 13) || !strncmp(line, "Templatename:", 13))
     {
       (void)strtok(line, ":");
       strcpy( ptr_isf->inst_filename, strtok(NULL,"\n"));
       strcat( ptr_isf->inst_filename,"\0");
     }
     else if (!strncmp(line, "seriesname:", 11) || !strncmp(line, "Seriesname:", 11))
     {
       (void)strtok(line, ":");
       strcpy( ptr_isf->series_name, strtok(NULL,"\n"));
       strcat( ptr_isf->inst_filename,"\0");
     }
     else if (!strncmp(line, "author:", 7) || !strncmp(line, "Author:", 7))
     {
       (void)strtok(line, ":");
       strcpy( ptr_isf->author, strtok(NULL,"\n"));
     }
     else if (!strncmp(line, "owner:", 6) || !strncmp(line, "Owner:", 6))
     {
       (void)strtok(line, ":");
       strcpy( ptr_isf->owner, strtok(NULL,"\n"));
     }
     else if (!strncmp(line, "description:", 12) || !strncmp(line, "Description:", 12))
     {
       (void)strtok(line, ":");
       strcpy( ptr_isf->packet_description, strtok(NULL,"\n"));
     }
     else if (!strncmp(line, "interval:", 9) || !strncmp(line, "Interval:", 9))
     {
       (void)strtok(line, ":");
       sscanf( strtok(NULL,"\n"), "%d", &(ptr_isf->interval));
     }
     else if (!strncmp(line, "keyword:", 8) || !strncmp(line, "Keyword:", 8))
     {
       if (top_if_kws == NULL)
       {
         top_if_kws = ptr_if_kws = malloc(sizeof(Instruction_File_Keywords_t));
         ptr_if_kws->next = (Instruction_File_Keywords_t*)NULL;
         ptr_isf->kw_if = top_if_kws;/*hook link list to Instruction_File_t struct*/
       }
       else
       {
         ptr_if_kws->next = malloc(sizeof(Instruction_File_Keywords_t));//link between nodes set here
         ptr_if_kws = ptr_if_kws->next;
         ptr_if_kws->next = (Instruction_File_Keywords_t*)NULL;
       }
       ptr_if_kws->next= NULL;

       /* set apid, shortname and longname values in Instruction_File_Keywords_t structure */
       /*parse line this->keyword line got:Keyword:19, HMI_TS09_TILT_MIRROR,      TILTMIRR */
       /* remove keyword: */
       (void)strtok(line, ":");
       /* Get apid */
       sscanf( strtok(NULL,","), "%d", &(ptr_if_kws->dec_apid));

       /* Get longname string */
       strcpy( kwlong, strtok(NULL,","));
       /* remove blanks in longname string */
       remove_blanks( kwlong, new_string);
       strcpy(ptr_if_kws->kw_longname, new_string);
       strcat(ptr_if_kws->kw_longname,"\0");

       /* Get short name */
       strcpy( kwshort, strtok(NULL,"\n"));
       /* remove blanks in strings */
       remove_blanks( kwshort, new_string);
       strcpy(ptr_if_kws->kw_shortname,new_string);
       strcat(ptr_if_kws->kw_shortname,"\0");

#ifdef  DEBUG_LM3S
       printf("read_isf_data:keywordshort : <%-s>\n", ptr_if_kws->kw_shortname);
       printf("read_isf_data:keywordlong  : <%-s>\n",ptr_if_kws->kw_longname);
#endif
        
     } /*end-else if keyword */
     else
     {
       printkerr("WARNING:Miss parsing line:<%s> from instruction file. Skip line and continue processing.\n",line);
     }

   }/* end while */

  /* close file and return structure which contains instruction file data */
  fclose(file_ptr);

  /*check got basic values instruction file */
  if (!strcmp(ptr_isf->inst_filename,""))
  {
    printkerr("ERROR:bad templatename setting <%s> from instruction file.\n",ptr_isf->inst_filename);
    exit_flag=1;
  }
  if (!strcmp(ptr_isf->series_name,""))
  {
    printkerr("ERROR:bad seriesname setting <%s> from instruction file.\n",ptr_isf->series_name);
    exit_flag=1;
  }
  if (!strcmp(ptr_isf->owner,""))
  {
    printkerr("ERROR:bad owner setting <%s> from instruction file.\n",ptr_isf->owner);
    exit_flag=1;
  }
  if (!strcmp(ptr_isf->author, ""))
  {
    printkerr("ERROR:bad author setting <%s> from instruction file.\n",ptr_isf->author);
    exit_flag=1;
  }
  if (ptr_isf->interval <= 0)
  {
    printkerr("ERROR:bad interval setting <%d> from instruction file.\n",ptr_isf->interval);
    exit_flag=1;
  }
  

  if (exit_flag)
  {
    printkerr("ERROR: Got exit flag. Exiting from executable.\n",line);
    return ((Instruction_File_Data_t *)(NULL));
  }
#ifdef  DEBUG_LM3S
  printf("read_isf_data:templatename      : <%s>\n",ptr_isf->inst_filename);
  printf("read_isf_data:seriesname        : <%s>\n",ptr_isf->series_name);
  printf("read_isf_data:owner             : <%s>\n", ptr_isf->owner);
  printf("read_isf_data:author            : <%s>\n",ptr_isf->author);
  printf("read_isf_data:description       : <%s>\n",ptr_isf->packet_description);
  printf("read_isf_data:interval          : <%d>\n",ptr_isf->interval);
#endif

  /* Return instruction filename */
  return(ptr_isf);
}/*end function:read_isf_data */



/*******************************************************************/
/* remove_blanks : remove blanks from string from instruction file */
/*******************************************************************/
void remove_blanks( char string[HKLMS_LONG_KEYWORD_NAME_SIZE], char *new_str)
{
  char new_string[HKLMS_LONG_KEYWORD_NAME_SIZE];
  int i,j;
  for(j=0; j < HKLMS_LONG_KEYWORD_NAME_SIZE; new_string[j]='\0',j++);
  /* remove blanks in strings */
  for(j=0, i=0;string[i];i++)
  {
    if(string[i] == ' ')
    {
      ;//skip
    }
    else
    {
      new_string[j]=string[i];
      j++;
    }
  }
  strcat(new_string,"\0");
  strcpy(new_str,new_string);
}
/*************************************/
/*      get minutes from pkt time    */
/*************************************/
int get_minute_from_pkttime(double tc_sec)
{
  /* declarations */
  int min;
  char at[HKLMS_PACKET_TIME_STR];

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UTC", 0);
  strcat(at,"\0");

  //2008.09.18_16:09:54_UTC
  /* get minute and return value */
  sscanf(at,"%*d.%*d.%*d_%*hd:%d",&min);
  return (min);
}
/************************************/
/*      get seconds from pkt time   */
/************************************/
int get_seconds_from_pkttime(double tc_sec)
{
  /* declarations */
  int sec;
  char at[HKLMS_PACKET_TIME_STR];

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UTC", 0);
  strcat(at,"\0");

  //2008.09.18_16:09:54_UTC
  /* get seconds and return value */
  sscanf(at,"%*d.%*d.%*d_%*hd:%*d:%d",&sec);
  return (sec);
}
/**********************************************/
/*  Save Keyword Data for Interval in struct */
/**********************************************/
HK_KW_Data_t * save_kw_data(TIME start_pkt_time, Instruction_File_Data_t *ifs, HK_Keyword_t *kw)
{

  HK_KW_Data_t *ptr_kw_data,*top_kw_data;
  HK_Keyword_t *hk_kw;
  Instruction_File_Keywords_t *kwif;

  /*initialize variables */
  top_kw_data=NULL;
  kwif= ifs->kw_if;

  /* loop thru list of keyword to lookup up values for in HK_Keyword_t structure */
  while (kwif)
  {
#ifdef DEBUG_LM3S
     printf("save_kw_data: longname:%s values in HK_KW_Data_t structure.\n",kwif->kw_longname);
#endif
     /* create HK_KW_Data node using longname and start_pkt_time */
     if(top_kw_data == NULL)
     {
       top_kw_data = ptr_kw_data = malloc(sizeof(HK_KW_Data_t));
       ptr_kw_data->next = (HK_KW_Data_t*)NULL;
       ptr_kw_data->kw_values = (HK_KW_Data_Values_t*)NULL;
     }
     else
     {
       ptr_kw_data->next = malloc(sizeof(HK_KW_Data_t));//link between nodes set here
       ptr_kw_data = ptr_kw_data->next;
       ptr_kw_data->next = (HK_KW_Data_t*)NULL;
       ptr_kw_data->kw_values = (HK_KW_Data_Values_t*)NULL;
     }

     /* set values in HK_KW_Data_t node */
     ptr_kw_data->start_pkt_time = start_pkt_time;
     strcpy(ptr_kw_data->longname,kwif->kw_longname);
     strcpy(ptr_kw_data->shortname,kwif->kw_shortname);//new shortname
     
     /* load in data values in HK_KW_Data_Values_t structure */
     for(hk_kw=kw;hk_kw; hk_kw=hk_kw->next)
     {
        if (!strcmp(hk_kw->name, kwif->kw_longname))
        {
           break;
        }
     } 

     /* go to next HK_KW_Data_t node */
     kwif=kwif->next;
  }
  return (top_kw_data);
}
/****************************************************/
/*  Save Keyword Data Values for Interval in struct */
/****************************************************/
void save_kw_data_values(HK_KW_Data_t *top_kw_data, TIME pkt_time, HK_Keyword_t *kw )
{
  // ####logical steps####
  // get first kw 
  // lookup value
  // save value in HK_KW_Data_t-->HK_KW_Data_Values_t
  // get next kw and repeat steps

  /* variables */
  HK_KW_Data_Values_t *ptr_kw_data_values;
  HK_KW_Data_t *ptr_kw_data; 
  HK_Keyword_t *ptr_kw; 
   
  /* loop thru each longname in HK_Data_t list and find longname in HK_Keyword_t list */
  for( ptr_kw_data= top_kw_data; ptr_kw_data;ptr_kw_data=ptr_kw_data->next)
  {
    for(ptr_kw=kw; ptr_kw; ptr_kw=ptr_kw->next)
    {
      /* check if found match when looked up longname in HK_Keyword_t link list */
      if(!strcmp(ptr_kw->name, ptr_kw_data->longname))
      {
#ifdef DEBUG_LM3S
        printf("save_kw_data_values:looking for longname in HK_KW_Data_Values_t:  %s \n",ptr_kw_data->longname);
        printf("save_kw_data_values:found match in HK_Keyword_t: %s \n",ptr_kw->name);
#endif
         
        /* load HK_Data_Values nodes */
        ptr_kw_data_values=ptr_kw_data->kw_values;
        if(ptr_kw_data_values == NULL)
        {
          ptr_kw_data_values = malloc(sizeof(HK_KW_Data_Values_t));
          ptr_kw_data_values->next = (HK_KW_Data_Values_t*)NULL;
          ptr_kw_data->kw_values=ptr_kw_data_values; //link to HK_Data_t node
        }
        else
        {
          /* go to last node in HK_KW_Data_Values_t link list to add new node */
          while(ptr_kw_data_values->next)
          {
            ptr_kw_data_values=ptr_kw_data_values->next;
          }

          /*link between nodes set here */
          ptr_kw_data_values->next = malloc(sizeof(HK_KW_Data_Values_t));
          ptr_kw_data_values = ptr_kw_data_values->next;
          ptr_kw_data_values->next = (HK_KW_Data_Values_t*)NULL;
        }

        /* set kw data values based on variable type of keyword */
        ptr_kw_data_values->pkt_time=pkt_time;
        ptr_kw_data_values->eng_type=ptr_kw->eng_type;

        if(ptr_kw_data_values->eng_type == KW_TYPE_DOUBLE)
        {
           ptr_kw_data_values->eng_value.double_val=ptr_kw->eng_value.double_val;
        }
        else if(ptr_kw_data_values->eng_type ==  KW_TYPE_UINT8) 
        {
           ptr_kw_data_values->eng_value.double_val= (double)ptr_kw->eng_value.uint8_val; 
        }
        else if(ptr_kw_data_values->eng_type ==  KW_TYPE_INT8) 
        {
           ptr_kw_data_values->eng_value.double_val= (double)ptr_kw->eng_value.int8_val; 
        }
        else if(ptr_kw_data_values->eng_type ==  KW_TYPE_UINT16) 
        {
           ptr_kw_data_values->eng_value.double_val= (double)ptr_kw->eng_value.uint16_val; 
        }
        else if(ptr_kw_data_values->eng_type ==  KW_TYPE_INT16) 
        {
           ptr_kw_data_values->eng_value.double_val= (double)ptr_kw->eng_value.int16_val; 
        }
        else if(ptr_kw_data_values->eng_type == KW_TYPE_UINT32)
        {
           ptr_kw_data_values->eng_value.double_val=(double)ptr_kw->eng_value.uint32_val;
        }
        else if(ptr_kw_data_values->eng_type == KW_TYPE_INT32)
        {
           ptr_kw_data_values->eng_value.double_val=(double)ptr_kw->eng_value.int32_val;
        }
        else if(ptr_kw_data_values->eng_type == KW_TYPE_FLOAT)
        {
           ptr_kw_data_values->eng_value.double_val=(double)ptr_kw->eng_value.float_val;
        }
        else 
        {
          ptr_kw_data_values->eng_value.double_val=0.00;
          printkerr("WARNING at %s, line %d: Type '%d' not handled. Probably String or Time keywords used. This is not allowed.\n",
                     __FILE__, __LINE__, ptr_kw_data_values->eng_type);
          printf("WARNING at %s, line %d: Type '%d' not handled.\n", __FILE__, __LINE__, ptr_kw_data_values->eng_type);
        }

#ifdef DEBUG_LM3S
        printf("save_kw_data_values:ptr_kw_data_values->pkt_time is %f\n",ptr_kw_data_values->pkt_time);
        printf("save_kw_data_values:ptr_kw_data_values->eng.value. is %f\n",ptr_kw_data_values->eng_value.double_val);
#endif

        /* break from looking and get need keyword longname to lookup*/
        break;

      }/* if found match */
    }/* inner for loop to find longname in HK_Keyword_t link list */
  }/* outer for loop to traverse thru each longname in HK_Data_t list */
} /* save_kw_data_values*/

/****************************************************/
/*  Free Keyword Data stucture for Interval         */
/****************************************************/
void free_all_kw_data(HK_KW_Data_t *top_kw_data_ptr)
{
  HK_KW_Data_t *ptr_kw_data, *temp_ptr_kw_data;
  HK_KW_Data_Values_t *ptr_kw_data_values, *temp_ptr_kw_data_values;

  /* got to first HK_Data_t node and got thru and free each HK_Data_Value_t node, then free HK_Data_t node */ 
  for ( ptr_kw_data= top_kw_data_ptr; ptr_kw_data;)
  {
#ifdef DEBUG_LM3S
      printf("free_all_kw_data:1:at HK_Data_t node with start time  :%f\n",  ptr_kw_data->start_pkt_time);
      printf("free_all_kw_data:2:at HK_Data_t node with longname    :%s\n",  ptr_kw_data->longname);
#endif

      for(ptr_kw_data_values = ptr_kw_data->kw_values;ptr_kw_data_values;)
      {
#ifdef DEBUG_LM3S
         printf("free_all_kw_data:3:will free HK_Data_Values_t node:%p\n",ptr_kw_data_values);
         printf("free_all_kw_data:4:where HK_Data_Values_t node pkt time:%f\n", ptr_kw_data_values->pkt_time);
#endif

         temp_ptr_kw_data_values = ptr_kw_data_values->next;  
         free(ptr_kw_data_values);
         ptr_kw_data_values=temp_ptr_kw_data_values;
      }/* inner for loop- loop thru HK_Data_Values_t node  */

      /* get next HK_Data_t node to traverse thru */
      temp_ptr_kw_data=ptr_kw_data->next;

#ifdef DEBUG_LM3S
      printf("free_all_kw_data:5:will free HK_Data_t node:%p\n",ptr_kw_data);
#endif

      free(ptr_kw_data);
      ptr_kw_data=temp_ptr_kw_data;
  }/* outer for loop - loop thru HK_Data_t nodes */

}
/****************************************************/
/*  Save Mean,Max,Min and Standard deviation Data   */
/****************************************************/
void save_m3sd_data(HK_KW_Data_t *top_ptr_kw_data ,HK_Keyword_M3SD_Data_t **kw_data_head)
{
  /* local variables */
  HK_KW_Data_t *ptr_kw_data ;
  HK_Keyword_M3SD_Data_t *ptr_m3sd_data;
  HK_Keyword_M3SD_Data_t *prev_ptr_m3sd_data;

  /* initalize HK_Keyword_M3SD_Data_t pointers to null */
  ptr_m3sd_data=*kw_data_head;

  /* loop thru each keyword & add names HK_Keyword_M3SD_Data_t * 
   * nodes and values to HK_Keyword_M3SD_Values_t nodes        */
  for(ptr_kw_data=top_ptr_kw_data; ptr_kw_data; ptr_kw_data= ptr_kw_data->next)
  {
    /* check if first node of HK_Keyword_M3SD_Data_t link list */
    if(!ptr_m3sd_data)
    {
      *kw_data_head = ptr_m3sd_data = malloc(sizeof(HK_Keyword_M3SD_Data_t));
      ptr_m3sd_data->next = (HK_Keyword_M3SD_Data_t*)NULL;
      ptr_m3sd_data->m3sd_values = (HK_Keyword_M3SD_Values_t*)NULL;
     }
     else
     {
       /* go to end of link list NEW-8-26-2009:4:30PM*/
       for(prev_ptr_m3sd_data =ptr_m3sd_data = *kw_data_head; ptr_m3sd_data;prev_ptr_m3sd_data=ptr_m3sd_data,ptr_m3sd_data=ptr_m3sd_data->next);
       prev_ptr_m3sd_data->next = ptr_m3sd_data = malloc(sizeof(HK_Keyword_M3SD_Data_t));//link between nodes set here
       ptr_m3sd_data->next = (HK_Keyword_M3SD_Data_t*)NULL;
       ptr_m3sd_data->m3sd_values = (HK_Keyword_M3SD_Values_t*)NULL;
     }

     /* set values in HK_Keyword_M3SD_Data_t node */
     ptr_m3sd_data->start_pkt_time = ptr_kw_data->start_pkt_time; //NEW-8-26
     ptr_m3sd_data->number_points = (int)get_number_points(ptr_kw_data->kw_values);
#ifdef DEBUG_LM3S
     printf("save_m3sd_data:set start time %f\n", ptr_m3sd_data->start_pkt_time);
     printf("save_m3sd_data:set number of points %d\n", ptr_m3sd_data->number_points);
#endif
     break; 
    }/* end-forloop*/

}
/****************************************************/
/*  Save Mean,Max,Min and Standard deviation Values  */
/****************************************************/
void save_m3sd_values(HK_KW_Data_t *top_ptr_kw_data, HK_Keyword_M3SD_Data_t **kw_data_head)
{
  HK_KW_Data_t *ptr_kw_data;
  HK_Keyword_M3SD_Data_t *ptr_m3sd_data;
  HK_Keyword_M3SD_Values_t *ptr_m3sd_values;

  /* initalize local H_Keyword_M3SD_Data_t pointer */
  ptr_m3sd_data=*kw_data_head;

   for( ptr_kw_data= top_ptr_kw_data; ptr_kw_data;ptr_kw_data=ptr_kw_data->next)
   {
     for(ptr_m3sd_data=*kw_data_head; ptr_m3sd_data; ptr_m3sd_data=ptr_m3sd_data->next)
     {
       /* check if found match when looked up longname in HK_Keyword_t link list */
       if(((int)ptr_m3sd_data->start_pkt_time) == ((int)ptr_kw_data->start_pkt_time))
       { 
            /* get start time from HK_Data_t node  */
            /* go to first HK_Keyword_M3SD_Values_t node and set      */
            /* set ptr to top HK_Keyword_M3SD_Values_t node */
            ptr_m3sd_values= ptr_m3sd_data->m3sd_values;

            if(!ptr_m3sd_values)
            {
              /* got first node for HK_Keyword_M3SD_Values_t link list */
              ptr_m3sd_values = malloc(sizeof(HK_Keyword_M3SD_Values_t));
              ptr_m3sd_values->next = (HK_Keyword_M3SD_Values_t*)NULL;

              /* link HK_Keyword_M3SD_Data_t node to  HK_Keyword_M3SD_Values_t node */
              ptr_m3sd_data->m3sd_values= ptr_m3sd_values;//added FRIDAY-8-19-2009 at 3:30PM
            }
            else
            {
              /* go to last node in HK_M3SD_Data_Values_t link list to add new node */
              while(ptr_m3sd_values->next)
              {
                ptr_m3sd_values=ptr_m3sd_values->next;
              }
 
              ptr_m3sd_values->next = malloc(sizeof(HK_Keyword_M3SD_Values_t));//link between nodes set here
              ptr_m3sd_values = ptr_m3sd_values->next;
              ptr_m3sd_values->next = (HK_Keyword_M3SD_Values_t*)NULL;
            }

            /* set longname and shortname keywords in HK_Keyword_M3SD_Values_t node */
            strcpy(ptr_m3sd_values->kw_longname,ptr_kw_data->longname);
            strcpy(ptr_m3sd_values->kw_shortname, ptr_kw_data->shortname);
            ptr_m3sd_values->min = get_min_value(ptr_kw_data->kw_values);
            ptr_m3sd_values->max = get_max_value(ptr_kw_data->kw_values);
            ptr_m3sd_values->mean = get_mean_value(ptr_kw_data->kw_values);
            ptr_m3sd_values->stdev = get_stdev_value(ptr_kw_data->kw_values);
            break;
       }/* if found match of start pkt time between M3SD Struct and DATA struct*/
     }/* for  HK_Keyword_M3SD_Data_t nodes to find longname to set m3sd and start packet values */
   }/* for loop thru HK_Data_t nodes for each longname*/
}

/**************************************************************************
 * Set Environment Variables                                              *
 * FUNCTION: set_env_variables()                                          *
 * DESCRIPTION: Sets Data series project name and data type name.         *
 * Sets path to JSOC Version number map files to lookup JSOC Version      *
 * number.                                                                *
 **************************************************************************/
void set_env_variables()
{
  FILE * fp;
  char envfile[500], s1[256],s2[256],s3[256], line[256];
  /* set environment variables for hk code that decodes dayfiles*/
  /* create filename and path */
  strcpy(envfile, ENVFILE );
  /* fopen file */
  if(!(fp=fopen(envfile, "r"))) 
  {
      printf("ERROR:Can't open environment variable file <%s>. Check setting is correct.\n", envfile);
      exit(0);
  }
  /* read in lines */
  while( fgets(line, MAXLINE_IN_FILE, fp) != NULL )
  {
    if (!strncmp(line, "#", 1)) 
    {
      continue; /*skip comments */
    }
    else  
    {
      sscanf(line, "%s %s %s ", s1,s2,s3);
      /* set each line to setenv(); */
      setenv(s2, s3, 1);
    }
  }
  /* close file */
  fclose(fp);
}
            
/************************/
/* GET NUMBER OF POINTS */
/************************/
int get_number_points(HK_KW_Data_Values_t *ptr_data_values)
{
  int num_pts;
  for(num_pts=0; ptr_data_values; num_pts++, ptr_data_values=ptr_data_values->next);
  return(num_pts);
}

/******************/
/* GET MIN VALUE  */
/******************/
float get_min_value(HK_KW_Data_Values_t *ptr_data_values)
{
  float min_value=0.0;
  for(min_value=(float)ptr_data_values->eng_value.double_val, ptr_data_values=ptr_data_values->next;
      ptr_data_values; ptr_data_values=ptr_data_values->next)
  {
      if(min_value > ptr_data_values->eng_value.double_val)
      {
         min_value= (float)ptr_data_values->eng_value.double_val;
      }
  }
#ifdef DEBUG_LM3S
  printf("get_min_value(): MIN->eng.value. is %f\n",min_value);
#endif
  return(min_value);
}

/******************/
/* GET MAX VALUE  */
/******************/
float get_max_value(HK_KW_Data_Values_t *ptr_data_values)
{
  float max_value=0.0;
  for(max_value = (float)ptr_data_values->eng_value.double_val, ptr_data_values=ptr_data_values->next;
      ptr_data_values; ptr_data_values=ptr_data_values->next)
  {
      if( ptr_data_values->eng_value.double_val > max_value)
      {
         max_value= (float)(ptr_data_values->eng_value.double_val);
      }
  }
#ifdef DEBUG_LM3S
  printf("get_max_value(): MAX->eng.value. is %f\n",max_value);
#endif
  return(max_value);
}

/******************/
/* GET MEAN VALUE */
/******************/
float get_mean_value(HK_KW_Data_Values_t *ptr_data_values)
{
  float mean_value=0.0;
  float count=0.0;
  for(; ptr_data_values; ptr_data_values=ptr_data_values->next)
  {
      mean_value += (float)ptr_data_values->eng_value.double_val;
      count++;
#ifdef DEBUG_LM3S
      printf("get_mean_value(): HK_KW_Data_Values->eng.value. is %f\n",ptr_data_values->eng_value.double_val);
#endif
  }
#ifdef DEBUG_LM3S
  printf("get_mean_value(): MEAN->eng.value. is %f\n",mean_value/count);
#endif
  return(mean_value/count);
}

/**************************/
/* GET STANDARD DEVIATION */
/**************************/
float get_stdev_value(HK_KW_Data_Values_t *top_ptr_data_values)
{
  float sum_value=0.0;
  float sum_sq_value=0.0;
  float numpts=0.0;
  int int_numpts=0;
  float mean_value=0.0;
  float mean_sq_value=0.0;
  float cal_value=0.0;
  HK_KW_Data_Values_t *ptr_data_values;
  float stdev_value=0.0;

  /* get mean */
  for(ptr_data_values=top_ptr_data_values, numpts=0.0; ptr_data_values; ptr_data_values=ptr_data_values->next)
  {
      sum_sq_value += (float)(powf((float)ptr_data_values->eng_value.double_val, 2.0 ));
      sum_value    += (float)ptr_data_values->eng_value.double_val;
      numpts++;
      int_numpts++;
  }
  /* get mean of values */
  mean_value=sum_value/numpts;

  /* get mean to power of 2 squared */
  mean_sq_value = (float)(powf(mean_value, 2.0 ));

  /* get standard deviation value */
  cal_value= (((float)(sum_sq_value/(numpts - 1.0)) ) - ((float)((numpts/(numpts - 1.0)) * mean_sq_value)));
  if (cal_value < 0.0 || int_numpts == 1)
  {
    stdev_value=0.0;
  }
  else
  {
    stdev_value = sqrtf( (float)cal_value );
  }

#ifdef DEBUG_LM3S
  printf("get_stdev_value(): STANDARD-DEV->eng.value. is %15.15f\n",stdev_value);
#endif
  return(stdev_value);
}

/**********************/
/* WRTIE M3SD TO DRMS */
/**********************/
void write_m3sd_to_drms(HK_Keyword_M3SD_Data_t *top_m3sd_data_ptr, Instruction_File_Data_t *ptr_inst_file)
{
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  HK_Keyword_M3SD_Data_t *ptr_m3sd_data;
  HK_Keyword_M3SD_Values_t *ptr_m3sd_values;
  TIME pkt_time;
  char keyname[HKLMS_LONG_KEYWORD_NAME_SIZE];
  char query[HKLMS_MAX_DSNAME_STR];
  int  status, status_create, status_closed;
  int cntrec;
  int i;

  /* if data series name passed as argument use this to write keyword values to */
  if(!ptr_inst_file->series_name)
  {
    printkerr("ERROR at %s, line %d: Bad value for data series name:"
              "<%s>. Exiting program.\n", __FILE__,__LINE__, ptr_inst_file->series_name);
    return;
  }
  else
  {
     strcpy(query, ptr_inst_file->series_name) ; 
#ifdef DEBUG_LM3S
     printf("write_m3sd_to_drms:ptr_inst_file->series_name:<%s> query:<%s>\n", ptr_inst_file->series_name,query);
#endif
  }

  /*new: count number of records to write */
  for ( cntrec=0, ptr_m3sd_data=top_m3sd_data_ptr;ptr_m3sd_data; cntrec++,ptr_m3sd_data=ptr_m3sd_data->next);

  /* message to user */
  printf(". . . writing %d records to drms data series <%s>\n",cntrec, query);

  /* create records will need */
  /* create record in drms */
  rs = drms_create_records( drms_env, cntrec, query, DRMS_PERMANENT, &status_create);
  if (status_create) 
  {
    printkerr("ERROR at %s, line %d: drms_create_records "
              "returned a bad status. Status returned:<%d>. \n",
               __FILE__,__LINE__,status_create);
    if (DRMS_ERROR_UNKNOWNSERIES == status_create)
    {
      printkerr("ERROR at %s, line %d: Unknown DRMS Series. "
                "Please create series:<%s> and try again. \n",
                 __FILE__,__LINE__,query);
    }
    return;//new 9-28
  }

  for ( i=0, ptr_m3sd_data=top_m3sd_data_ptr;ptr_m3sd_data; ptr_m3sd_data=ptr_m3sd_data->next, i++)
  {
#ifdef DEBUG_LM3S
    printf("write_m3sd_to_drms:HK_Keyword_M3SD_Values_t node start_pkt_time<%f>\n",ptr_m3sd_data->start_pkt_time);
#endif

    /* get first and next records */
    rec = rs->records[i];

    pkt_time= ptr_m3sd_data->start_pkt_time;
    /* set T_START keyword using data from TIMECODE keywords and interval value*/
    /* set drms type and long telemetry name  */
    keytype= DRMS_TYPE_TIME;
    strcpy(keyname, "T_START");
    /* set packet time */
    key_anyval.time_val= pkt_time;
    /* set record */
    status = drms_setkey(rec, keyname, keytype, &key_anyval);
    /* check status */
    check_status_drms_set(status, keyname);
    
    /* set numpts keyword*/
    /* set drms type and long telemetry name  */
    keytype= DRMS_TYPE_INT;
    strcpy(keyname, "NUMPTS");
    /* set numpts keyword */
    key_anyval.int_val = (int32_t)(ptr_m3sd_data->number_points);
    /* set record */
    status = drms_setkey(rec, keyname, keytype, &key_anyval);
    /* check status */
    check_status_drms_set(status, keyname);

#ifdef DEBUG_LM3S
    printf("write_m3s_drms: NUMPTS is %d\n",(int32_t)(ptr_m3sd_data->number_points));
#endif

    /* load in min,max,mean and standard deviation keywords */
    for(ptr_m3sd_values = ptr_m3sd_data->m3sd_values;ptr_m3sd_values;ptr_m3sd_values=ptr_m3sd_values->next)
    {
      /* set longname,shortname,m3,sd */
      /* set drms type but promote up for unsigned values */
      keytype= DRMS_TYPE_FLOAT;

      /* set min */
      sprintf(keyname,"%s_%s", ptr_m3sd_values->kw_shortname,"MIN");
      strcat(keyname,"\0");
      key_anyval.float_val = ptr_m3sd_values->min;
      status = drms_setkey(rec, keyname, keytype, &key_anyval);
      /* check status */
      check_status_drms_set(status, keyname);

      /* set max */
      sprintf(keyname,"%s_%s", ptr_m3sd_values->kw_shortname,"MAX");
      strcat(keyname,"\0");
      key_anyval.float_val = ptr_m3sd_values->max;
      status = drms_setkey(rec, keyname, keytype, &key_anyval);
      /* check status */
      check_status_drms_set(status, keyname);

      /* set mean */
      sprintf(keyname,"%s_%s", ptr_m3sd_values->kw_shortname,"MEAN");
      strcat(keyname,"\0");
      key_anyval.float_val = ptr_m3sd_values->mean;
      status = drms_setkey(rec, keyname, keytype, &key_anyval);
      /* check status */
      check_status_drms_set(status, keyname);

      /* set standard deviation */
      sprintf(keyname,"%s_%s", ptr_m3sd_values->kw_shortname,"SD");
      strcat(keyname,"\0");
      key_anyval.float_val = ptr_m3sd_values->stdev;
      status = drms_setkey(rec, keyname, keytype, &key_anyval);
      /* check status */
      check_status_drms_set(status, keyname);

    }/* for-inner loop */

  } /* for loop thru each start time in HK_Keyword_M3SD_Data_t node */

  /* close record */
  status_closed = drms_close_records(rs, DRMS_INSERT_RECORD);
  if (status_closed) 
  {
    printkerr("ERROR at %s, line %d: drms_close_record "
              "returned a bad status. Status returned:<%d>. \n",
               __FILE__,__LINE__,status_closed);
  }

}

/**********************************************/
/*  CHECK STATUS DRMS SET                     */
/**********************************************/
void check_status_drms_set(int status, char *kwn)
{
  if (status)
  {
    printkerr("ERROR at %s, line %d: drms_setkey for keyword "
              "returned a bad status. Status returned:<%d> when setting "
              "keyword:<%s>\n", __FILE__,__LINE__,status, kwn);
  }
  return;
}
/**********************************************/
/*  Get Next Packet Time                      */
/**********************************************/
void get_next_pkt_time(TIME p_time, int intval , TIME *start, TIME *end)
{
  TIME start_day, start_pt_time, end_pt_time, find_time;
  char at[HKLMS_PACKET_TIME_STR];
  char start_str[HKLMS_PACKET_TIME_STR];
  int year, month, day;


  /*get year,month and day */
  year=get_yr_from_pkttime(p_time); 
  month=get_month_from_pkttime(p_time);
  day=get_day_from_pkttime(p_time);

  /* create start day string */
  sprintf(start_str,"%4.4d.%02.2d.%02.2d_00:00:00_UTC", year,month,day);
  strcat(start_str,"\0");

  /* create start time double(TIME) */
  start_day = sscan_time(start_str);

  /* find next packet time */
  find_time=start_day;
  while( find_time < p_time)
  {
    find_time += intval;
  }

  /*found end point of interval */
  end_pt_time=find_time;
  start_pt_time= find_time - intval;
  (void)sprint_time (at,  start_pt_time, "UTC", 0);
  strcat(at,"\0");
  (void)sprint_time (at,  end_pt_time, "UTC", 0);
  strcat(at,"\0");

  /* set return pointer values for start and end p_times to return back */
  *start= start_pt_time;
  *end=end_pt_time;
}

