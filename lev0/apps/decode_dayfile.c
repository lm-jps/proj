#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/lev0/apps/decode_dayfile.c,v 1.24 2011/03/15 21:01:16 carl Exp $"
/*############################################################################
# Name:        decode_dayfile.c - Decode Dayfiles                            #
# Description: Decode dayfile decodes hk packet by sending hk packets to     #
#              functions in decode_hk.c file. Keywords are decoded using the #
#              the same functions that decodes hk keywords on high speed bus.#
#              This decode_dayfiles loops through all hk packets and sends   #
#              data to be decoded to decode_hk.c and writes data to DRMS data#
#              series based on apid, project name, data type name, and jsoc  #
#              version number.When run without out parameter this code looks #
#              up name of data series name to put keyword names and values in#
# Execution:   decode_dayfile [-p] src=<src> in=<day-file>                   #
#                             [ out=<data series name>]                      #
#              (1)To create data out series on automatically:                #
#                     decode_dayfile src=<source>  in=<day-file>             #
#              (2)To specify data series :                                   #
#                     decode_dayfile src=<src> in=<day-file>                 #
#                                    out=<data series name>                  #
#              (3)To specify out data series and print to standard out :     #
#                  decode_dayfile -p src=<src> in=<day-file>                 #
#                                  out=<data series name>                    #
#              (4)create data out series automatically and print standard    #
#                 output to file.                                            #
#                  decode_dayfile -p  src=<src> in=<day-file>  > View-OUTPUT #
# Example 1:  decode_dayfile src=rtmon in=/tmp20/20070202.0x001d(Best Option)#
# Example 2:  decode_dayfile src=rtmon in=/tmp20/20070202.0x001d             #
#                             out=su_carl.lev0_0029_0001                     #
# Limitation:  Setup required for environment variables in file              #
#              SETENV_HK_DAYFILE_DECODE. Currently the code read packet size #
#              limit of 1000 bytes (HKDDF_MAX_PKT_SIZE ). Need to setup local#
#              version of SETENV_HK_DAYFILE_DECODE file since file checked in#
#              is production version of file. The source setting determines  #
#              parsing of file format therefore must match each src to its   #
#              standard file format to parse correctly. For example:         #
#                 src=moc expects file like 0129_2008_191_01.hkt             # 
#                 src=rtmon expects files like  2008_02_01.0x081             #
#                 src=hsb expects files like hsb_0445_2008_02_30_01.hkt      #
#                 src=egsefm expects files like  2008_02_01.0x081            #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on April 10, 2008 #
############################################################################*/

/**
   @defgroup decode_dayfile decode_dayfile
   @ingroup su_lev0

   @brief Decodes housekeeping packet data from dayfiles and writes decoded keywords and values to DRMS series.

   @par Synopsis:
   @code
   decode_dayfile -h
   decode_dayfile [-p] src=<moc|hsb|rtmon|egsefm> in=<day filename> [ out=<data series name> ] 
   @endcode

   Decode dayfile decodes hk packet by sending packet to functions in 
   decode_hk.c file. Keywords are decoded using the the same functions that 
   decodes hk keywords on high speed bus. This decode_dayfiles loops through 
   all hk packets and sends data to be decoded to decode_hk.c's functions 
   which returns a structure containing keyword names, keyword values, keyword
   and keyword variable types. Decode dayfile code then writes keywords names and 
   values to a DRMS data series. DRMS data series is based on apid(i.e.,0445,0529), 
   project name(i.e.,hmi,aia, sdo,etc.), data type name(i.e., lev0), and jsoc 
   version number(0001, 0002, etc.). Decode dayfile will not load records that 
   already exist in data series. Decode dayfile checks the packet time per packet 
   is within a time range of less than current date and time plus 12 hours. 
   Decode dayfile will not load data for packets that do not have HK configuration files.

   The src parameter is a mandatory argument which should be either moc, hsb, rtmon, or
   egsefm. This value is required to be set to correspond to the -in- parameter file's 
   format. If set src=moc the code expects file formats  
   <apid>_<year>_<day of year>_<version>.hkt (i.e., 0129_2008_191_01.hkt). If set 
   src=rtmon, the code expects file format <year>_<month>_<day of month>.0x<hex apid value>
   (i.e., 2008_02_01.0x081).  If set src=hsb, the code expects file format 
   hsb_<decimal apid value>_<year>_<month>_<day in month>_<version>.hkt 
   (i.e., hsb_0445_2008_02_30_01.hkt). If set src=egsefm, the code expects file format  
   <year>_<month>_<day in month>.0x<hex apid value> (i.e., 2008_02_30.0x081). Any other src
  setting will throw an error. The src value needs to match the -in- file format. LMSAL
  (src=rtmon and src=egsefm), Nasa(src=moc) and Stanford(src=hsb) each produce different 
  file formats for day files. The src is also used to set source value in the HK_SOURCE 
  keyword in Level 0 by apid data series.
   
   The in parameter is a mandatory argument which should contain the directory and 
   filename of the input dayfile. Currently only one dayfile is allowed for the in 
   parameter.

   The out parameter is optional. When run without out parameter, this 
   code creates the data series name to put keyword names and values in. 
   When run with the out parameter, this code will write keywords and values 
   to the data series name defined by the out parameter value. For either case
   there needs to be a data series already created. This code does not create
   a data series but determines and then creates if necessary the data series 
   name to write data to. 

   When use the -p flag, a report on each packet is send to standard output. 
   By using -p flag, the report output can be redirected to a file
   (i.e., > Report_ISP_DAYFILE_OCTOBER_29_2008).

   @par Flags:
   @c -h: Shows usage message.
   @c -p: Prints report on each packets keyword name and value to standout out
   @par

   @param src The source of dayfiles(required field).

   @param in The full directoy path and file name to the input dayfile(required field).

   @param out The data series name to write keyword names and values to(optional).
   Normally do not run with this option. If not present, the code will
   automatically create data series name using the apid and packet 
   version number(for hmi and aia packets) in packets or packet time
   for sdo packets in day files.  Also required to set environment 
   variables for project names and data identifer name: 
   HK_DDF_HKBYAPID_PROJECT_NAME_xxx & HK_DDF_HKBYAPID_DATA_ID_NAME. These variables are in the the
   SOURCE_ENV_FOR_HK_DAYFILE_DECODE file which is read in based on the setting of
   define ENVFILE . The data series needs to be already created for this to program 
   to load keywords in data series.  This code only creates the name of the data 
   series to write to, it does not create data series. 

   @par Example of running without out parameter using as in parameter an ISP dayfile from MOC Product Server:
   @code
   decode_dayfile src=moc in=/tmp21/production/lev0/hk_moc_dayfile/0029_2008_275_01.hkt 
   @endcode

   @par Example of running with -p flag, without out parameter using as in parameter an ISP dayfile from MOC Product Server:
   @code
   decode_dayfile -p src=moc in=/tmp21/production/lev0/hk_moc_dayfile/0029_2008_275_01.hkt  > Report_ISP_0029_2008_275_01.hkt 
   @endcode

   @par Example of running with out parameter using as in parameter an ISP dayfile from MOC Product Server:
   @code
   decode_dayfile src=moc in=/tmp21/production/lev0/hk_moc_dayfile/0029_2008_275_01.hkt out=hmi.lev0_0029_0008
   @endcode

   @par Example of running with out parameter using as in parameter an ISP dayfile file from MOC Product Server and sending error logsto file:
   @code
   decode_dayfile src=moc in=/tmp21/production/lev0/hk_moc_dayfile/0029_2008_275_01.hkt out=hmi.lev0_0029_0008 > ERROR_LOG
   @endcode

   @par Example of running with -p print flag with out parameter using as in parameter an ISP dayfile from MOC Product Server:
   @code
   decode_dayfile -p src=moc  in=/tmp21/production/lev0/hk_moc_dayfile/0029_2008_275_01.hkt out=hmi.lev0_0029_0008 >  Report_ISP_0029_2008_275_01.hkt
   @endcode

   @par Example of running help:
   @code
   decode_dayfile -h
   @endcode

*/
/* Defined constants */
/******************** defines ***********************************************/
#define HKDDF_CNTRL_L            ""
#define HKDDF_MAXLINE_IN_FILE    200
#define HKDDF_MAXNUM_PJNAMES     3
#define HKDDF_MAX_APID_STR       5
#define HKDDF_MAX_DSNAME_STR     100
#define HKDDF_MAX_FILE_NAME      100
#define HKDDF_MAX_JSVN_STR       5
#define HKDDF_MAX_LONG_KW_NAME   50
#define HKDDF_MAX_PKT_SIZE       1000
#define HKDDF_MAX_PJNAME         50
#define HKDDF_MAX_PVN_SIZE       50
#define HKDDF_MAX_PVN_STR        10
#define HKDDF_MAX_HK_SRC_KW      100
#define HKDDF_MAX_VERSION_LINES  1000
#define HKDDF_START_MERGED_PVNW  (1)
#define HKDDF_START_MERGED_PVND  (194)
#define HKDDF_MAX_PACKET_NAME    50
#define HKDDF_MAX_LOOKUP_NAME    100
#define HKDDF_MAX_SRC_ARG        10
#define HKDDF_MAX_IN_OUT_ARG     200
#define HKDDF_ENVIRONMENT_VARS_NOT_SET  (-1)
/*#define ENVFILE    "/home/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DAYFILE_DECODE"*/
/*#define ENVFILE    "/home/carl/cvs/myprod/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DAYFILE_DECODE"*/
#define ENVFILE      "/home/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DAYFILE_DECODE"
#define HKDDF_READ_ARRAY_SIZE  (25000001)

/******************** includes ******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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
  {ARG_STRING, "src", "Not Specified", "day file source"},
  {ARG_STRING, "in", "Not Specified", "full path to day file"},
  {ARG_STRING, "out", "Not Specified", "Series name"},
  {ARG_FLAG, "p", "0", "print values of keywords to standard out"},
  {ARG_END}
};
ModuleArgs_t   *ggModArgs=module_args;
char* module_name = "decode_dayfile";

/******************* function prototypes  *******************************/
static char* get_ds_pjname(int apid, char pjn[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME]);
static int   get_packet_time_for_df(HK_Keyword_t *hk, TIME *ptime);
static int   load_ds_names(char *pn, char* didn, int apid, char *pvn);
static char* lookup_dsn( int apid, char *pvn,char *fvn, int count);
static void  print_packet_time(HK_Keyword_t *kwh);
static void  save_packet_values1(unsigned char *read_in_buffer, char *out, char *in, char *src);
static void  save_packet_values2(unsigned char *read_in_buffer, char projname[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME],char *didn, char *in, char *src);
static void  saveprint_packet_values1(unsigned char *read_in_buffer, char *in, char *out, char *src);
static void  saveprint_packet_values2(unsigned char *read_in_buffer, char *in, char projname[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME],char *didn, char *src);
static void  set_env_variables();
static void  write_to_drms( int apid, char packet_version_number[MAX_CHAR_VERSION_NUMBER], char file_version_number[MAX_CHAR_VERSION_NUMBER],char *data_series, HK_Keyword_t *kw_head, char *in, char *src);
static TIME  SDO_to_DRMS_time(int sdo_s, int sdo_ss);
static char *get_hk_src_keyword(int apid, char *in_file, char *src);
static int get_date_in_filename(char *inf, char *src, int *year, int *month, int *day);
static void print_hk_src_keyword(char *ifn, int apid, char *src);
static int check_df_env_variable(void);
static void get_month_day(int year, int yearday, int *pmonth, int *pday);

/********************* extern functions  *********************************/
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
int INVALtime;/*added to help compile on production*/

/********************* structures   *************************************/
/*structure used to save values from JSVN-TO-PVN map files */
struct jsvn_map_data
{
  char apid[HKDDF_MAX_APID_STR];
  char jvn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_JSVN_STR];
  char pvn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_PVN_STR];
  char fvn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_PVN_STR];
  char dsn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_DSNAME_STR];
} jsvn_list[1], *jmap=jsvn_list;

/********************* enums      *************************************/
/* used for projname array to lookup project name */
enum project{pjnHMI,pjnAIA,pjnSDO};

/* @} */



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
    printf ("Usage:\ndecode_dayfile [-p | -h] "
      "src=<day filename source> "
      "in=<day filename> "
      "[ out=<data series name> ] \n"
      "  details are:\n"
      "  -h: help - show this message then exit(optional field)\n"
      "  -p: print all keyword names and values to standard out(optional field)\n"
      "  src=<moc|hsb|rtmon|egsefm> -use to determine input filename source (required field).\n"
      "  user needs to set src value to correspond to file format of the -in- dayfile name.\n"
      "  in=<day file name> -use full path to day file(required field)\n"
      "  out=<data series name> -data series to write keyword values to in DRMS(optional field)\n"
      "  if don't enter out value, data series name will be created automatically using the apid and\n"
      "  packet version number in packets in day file. Also required to set environment variables\n"
      "  for project names and data identifer name: HK_DDF_HKBYAPID_PROJECT_NAME_xx,\n"
      "  HK_DDF_HKBYAPID__DATA_ID_NAME, HK_DDF_DF_PROJECT_NAME_xxx, and HK_DDF_DF__DATA_ID_NAME,\n"
      "  Also need data series already created for this to program to load keywords in data series.\n"
      "  Example of running using automatic lookup of data series name and without print on:\n"
      "           decode_dayfile src=moc in=/tmp20/production/hmi_hk/20070707.0x0029\n" );
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
  char *hk_df_fn;
  char hk_directory_filename[HKDDF_MAX_FILE_NAME];
  char *pn1,*pn2,*pn3;
  char projname[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME];
  char *didn;
  int out_flag;
  unsigned char *ptr_read_in_buffer;
  unsigned long int i;
  int print_flag;
  char src[HKDDF_MAX_SRC_ARG];
  char in[HKDDF_MAX_IN_OUT_ARG];
  char out[HKDDF_MAX_IN_OUT_ARG];

  /* set environment variables */
  set_env_variables();

  /* check environment variables are set */
  if(!check_df_env_variable())
  {
     exit(1);
  }

  /* parameter initialization */
  hk_df_fn= hk_directory_filename;
  out_flag=0;

  /* Get command line arguments */
  const int c_print_flag =cmdparams_get_int (&cmdparams, "p", NULL) != 0;
  const char *c_src = cmdparams_get_str (&cmdparams, "src", NULL);
  const char *c_in = cmdparams_get_str (&cmdparams, "in", NULL);
  const char *c_out = cmdparams_get_str (&cmdparams, "out", NULL);
  print_flag = c_print_flag;
  strcpy(src ,c_src);
  strcpy(in,c_in);
  strcpy(out,c_out);

  /* check arguments used */
  if (nice_intro ()) return (0);

  /* check if entered day file name */
  if (in == NULL) 
  {
    printkerr("ERROR at %s, line %d: Need to enter day file name.",
              "Exiting program.\nExample format for out file: "
              "$HOME/EGSE/tables/hk_sim_file/20060719.0x000f\n",
               __FILE__,__LINE__);
    return (0);
  }

  if (!strcmp(out, "Not Specified")) 
  {
    /*set out_flag to 0 to trigger creating automatically data series name  */
    out_flag=0;

    /* get series name project and data type name for creating HK BY APID series name */
    pn1 = getenv("HK_DDF_HKBYAPID_PROJECT_NAME_HMI");
    pn2 = getenv("HK_DDF_HKBYAPID_PROJECT_NAME_AIA");
    pn3 = getenv("HK_DDF_HKBYAPID_PROJECT_NAME_SDO");
    didn =getenv("HK_DDF_HKBYAPID_DATA_ID_NAME");

    /* set each project name in array */
    strcpy(projname[pjnHMI],pn1);
    strcpy(projname[pjnAIA],pn2);
    strcpy(projname[pjnSDO],pn3);
  }
  else
  {
    /* if have out argument set by user set out_flag to 1 */
    out_flag=1; 
  }

 
  /* get  in filename and open file*/
  strcpy(hk_df_fn, in);
  strcat(hk_df_fn, "\0");
  file_ptr=fopen(hk_df_fn,"r");
  if (!file_ptr)
  {
    printkerr("ERROR at %s, line %d: Please check filename and directory is correct. "
              " Could not get -in- directory and filename: "
              "<%s>. Exiting execution.\n", __FILE__,__LINE__, hk_df_fn);
    return(0);
  }

   ptr_read_in_buffer = (unsigned char *) malloc(sizeof(unsigned char) * HKDDF_READ_ARRAY_SIZE);
  /*read lines in file into buffer representing a packet*/
  for(i=0; i < HKDDF_READ_ARRAY_SIZE;i++) ptr_read_in_buffer[i]= '\0'; ;
  for(i = 0 ; fread(ptr_read_in_buffer + i,1,1,file_ptr) ; i++) 
  {  
    ;/* do nothing*/
    if( i == HKDDF_READ_ARRAY_SIZE - 1)
    {
      printkerr("ERROR at %s, line %d: Array for reading dayfile is too small. :"
                "<%lu>\n", __FILE__,__LINE__, i + 1);
      printkerr("Break up large dayfile using dd command. :   dd if<large-file> "
                " of=<small-file-1> bs=<packet_size + 7> count=<i.e., 1000> skip=<0,1,etc>\n");
      return (0);
    }
  }
  *(ptr_read_in_buffer + i)= '\0';
  fclose(file_ptr);
  if(print_flag && out_flag)
  {
    /* start printing packet information and process keywords to drms*/
    /*  and use command line argument for data series name           */
    saveprint_packet_values1(ptr_read_in_buffer, in, out, src);
    return 0;
  }
  else if(print_flag && !out_flag)
  {
    /* start printing packet information and process keywords to drms and */
    /* automatically create JSVN and use as data series name              */
    saveprint_packet_values2(ptr_read_in_buffer, in, projname, didn, src );
    return 0;
  }
  else if (!print_flag && out_flag)
  {
    /* do not print packet information and use command line argument for */
    /* data series name. this is case where: if(out_flag && !print_flag) */
    /* write to drms only - do not print values to standard out          */
    save_packet_values1(ptr_read_in_buffer, out, in, src);
    return 0;
  }
  else if (!print_flag && !out_flag)
  {
    /* do not print packet information and use automatic lookup of  */
    /* data series name. Environment variables need to be set for   */
    /* project name and data identifier name. write to drms only    */
    /* - do not print values to standard out                        */
    save_packet_values2(ptr_read_in_buffer, projname, didn, in, src);
    return 0;
  }
  else 
  {
    printkerr("ERROR at %s, line %d: Save and Print Flags are not set correctly."
              "This case should not occur. Check logic. Exiting Program.\n" ,
               __FILE__,__LINE__ );
    return 1;
  }
}



/*************************************************************************
 * WRITE TO DRMS                                                         *
 * FUNCTION: void write_to_drms(char [],char *,HK_Keyword_t *,char *in)  *
 * DESCRIPTION: Gets Packet time based on Time Codes. Set packet version *
 *              number or file version number and keywords.              *
 *************************************************************************/
void write_to_drms(int apid, char pkt_ver_num[MAX_CHAR_VERSION_NUMBER], char file_ver_num[MAX_CHAR_VERSION_NUMBER], char *ds_name, HK_Keyword_t *kw_head, char *in, char *src)
{
  /* variable definitions */
  char at[200];
  char query[HKDDF_MAX_DSNAME_STR];
  int status;
  int ck_status;
  int ck_rwtr_status;
  /* keyword variables and structure */
  HK_Keyword_t *kw;
  TIME pkt_time;
  TIME pkt_timestamp;
  kw= kw_head;

  /* variable to set drms record */
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  char keyname[HKDDF_MAX_LONG_KW_NAME];
  /* drms record create variables */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;

  /* if data series name passed as argument use this to write keyword values to */
  if(ds_name)
  {
    strcpy(query, ds_name) ;  
  }
  else
  {
    printkerr("ERROR at %s, line %d: Bad value for data series name:"
              "<%s>. Exiting program.\n", __FILE__,__LINE__, ds_name);
    return;
  }

  /* check if record is within 12 hour range to remove processing of future dated packet data */
  kw=kw_head;


  ck_rwtr_status = check_hk_record_within_time_range(kw);
  if (ck_rwtr_status == 1)
  {
    /* now check if record is exits in drms series */
    kw=kw_head;
    ck_status = check_hk_record_exists(ds_name, kw, apid);
  }
  else
  {
    (void)get_packet_time_for_df(kw, &pkt_timestamp);
    (void)sprint_time (at,  pkt_timestamp, "UTC", 0);
    strcat(at,"\0");
    printkerr("WARNING:decode_dayfile: Skipping write of record to  <%s> data series. Packet time:<%s> Record not within time range.\n", ds_name, at);
    return;
  }

  /* if already exits skip writing record to DRMS data series */
  if (ck_status == 1)
  {
    ;/* skip setting */
#ifdef DEBUG_DECODE_DAYFILE
      printf("DEBUG:Message:decode_dayfile: Skipping write of record to  <%s> data series - record already exits\n", ds_name);
#endif
  }
  else
  { /* big else - create drms record and populate values */
#ifdef DEBUG_DECODE_DAYFILE
      printkerr("DEBUG:Message at %s, line %d: Prepare write of record to  <%s> data series.\n" ,
               __FILE__, __LINE__, ds_name);
#endif

    /* create record in drms */
    rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);
    if (status)
    {
      printkerr("ERROR at %s, line %d: Cannot create record using this data"
                " series name:<%s>.Exiting write to drms.\n",__FILE__,__LINE__, query);
      return;
    }
    else
    {
      ;//printf("Sucessfully created record for series <%s>\n", query);
    }
    rec = rs->records[0];

    /* get and set packet_time in record*/
    if (!get_packet_time_for_df(kw, &pkt_time))
    {
      /*did not set PACKET_TIME */
      printkerr("ERROR at %s, line %d: Could not set PACKET_TIME keyword.",
                "Exiting program.\n", __FILE__,__LINE__);
      return;
    }
    else 
    { 
      /* set PACKET_TIME keyword using data from TIMECODE keywords*/
      /* set drms type and long telemetry name  */
      keytype= DRMS_TYPE_TIME;
      strcpy(keyname, "PACKET_TIME");
      /* set packet time */
      key_anyval.time_val= pkt_time;
      /* set record */
      status = drms_setkey(rec, keyname, keytype, &key_anyval);
    }

    /* if sdo hk packet then set file version number, otherwise set packet version number */
    if(check_for_sdo_apid(apid))
    {
      /* for sdo hk packets */
      /* set FILE_VERSION_NUMBER keyword using passed value */
      /* set drms type and long telemetry name  */
      keytype= DRMS_TYPE_STRING;
      strcpy(keyname, "FILE_VERSION_NUMBER");
      /*allocate memory for string value */
      key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
      /* set packet version number */
      strcpy(key_anyval.string_val, file_ver_num);
    }
    else
    {
      /* set PACKET_VERSION_NUMBER keyword using data from HMI_VER.. keyword */
      /* set drms type and long telemetry name  */
      keytype= DRMS_TYPE_STRING;
      strcpy(keyname, "PACKET_VERSION_NUMBER");
      /*allocate memory for string value */
      key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
      /* set packet version number */
      strcpy(key_anyval.string_val, pkt_ver_num);
    }
    /* set record */
    status = drms_setkey(rec, keyname, keytype, &key_anyval);
    /* free memory */
    free (key_anyval.string_val);

    /* set HK_SOURCE keyword */
    /* set drms type */
    keytype= DRMS_TYPE_STRING;
    /* set  -long- telemetry name (HK_SOURCE) keyword for lev0 by APID data series */
    strcpy(keyname, "HK_SOURCE");
    /* allocate memory for string value */
    key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
    /* set string:<dayfile_name>.[<date>][<packet-apid>][<source>] */
    strcpy(key_anyval.string_val, get_hk_src_keyword(apid, in, src));
    status = drms_setkey(rec, keyname, keytype, &key_anyval);
    /* free memory for string */
    free (key_anyval.string_val);
    
    /* loop through keyword struct and load in record */
    kw= kw_head;
    while (kw)
    {
      if(kw->eng_type == KW_TYPE_UINT8)
      {
        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_SHORT;
        strcpy(keyname, kw->name);
        /* set drms value by casting up for unsigned values */
        key_anyval.short_val = (int16_t)kw->eng_value.uint8_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_UINT16)
      {
        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_INT;
        strcpy(keyname, kw->name);
        /* set drms value by casting up for unsigned values */
        key_anyval.int_val = (int32_t)kw->eng_value.uint16_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_UINT32)
      {
        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_LONGLONG;
        strcpy(keyname, kw->name);
        /* set drms value by casting up for unsigned values */
        key_anyval.longlong_val = (int64_t)kw->eng_value.uint32_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_INT8)
      {
        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_SHORT;
        strcpy(keyname, kw->name);
        key_anyval.short_val = (int16_t)kw->eng_value.int8_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_DOUBLE)
      {
        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_DOUBLE;
        strcpy(keyname, kw->name);
        key_anyval.double_val = kw->eng_value.double_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_FLOAT)
      {
        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_FLOAT;
        strcpy(keyname, kw->name);
        key_anyval.float_val = kw->eng_value.float_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_INT16)
      {
        /* set drms type  */
        keytype= DRMS_TYPE_SHORT;
        strcpy(keyname, kw->name);
        key_anyval.short_val = kw->eng_value.int16_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_INT32)
      {
        /* set drms type  */
        keytype= DRMS_TYPE_INT;
        strcpy(keyname, kw->name);
        key_anyval.int_val= kw->eng_value.int32_val;
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
      else if(kw->eng_type == KW_TYPE_STRING)
      {
        /* set drms type  */
        keytype= DRMS_TYPE_STRING;
        strcpy(keyname, kw->name);
        /*allocate memory for string */
        key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
        strcpy(key_anyval.string_val, kw->eng_value.string_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
        free(key_anyval.string_val);
      }
      else
      {
        printkerr("Warning at %s, line %d: Found Unknown KW_TYPE for this keyword:"
                  "<%d>. Check if new KW_TYPE has been created in HKPDF files"
                  "Skipping recording keyword in structure. Keyword is <%s>.",
                   __FILE__,__LINE__, (int)kw->eng_type, kw->name);
        /**TO DO:SKIP setting record but increment to next keyword node */
      }
      /* Next Keyword */
      kw=kw->next;
    }
    /* close record */
    status = drms_close_records(rs, DRMS_INSERT_RECORD);
    if (status)
    {
      printkerr("ERROR at %s, line %d: Cannot close drms record.\n", 
                 __FILE__,__LINE__);
      return;
    }
    else
    {
      ;//printf("Successfully closed record\n");
#ifdef DEBUG_DECODE_DAYFILE
      printkerr("DEBUG:Message at %s, line %d: Completed close of write of record to  <%s> data series.\n" ,
               __FILE__, __LINE__, ds_name);
#endif
    }
  }/* big else- create new drms records and populate */
  return;
}



/**********************************************************************************
 * TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);
 **********************************************************************************/
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



/*************************************************************************
 * GET PACKET TIME FOR DAY FILE                                          *
 * FUNCTION: int get_packet_time_for_df(HK_Keyword_t *,  char *, TIME )  *
 * DESCRIPTION: Gets Packet time based on Time Codes.                    *
 *************************************************************************/
int get_packet_time_for_df(HK_Keyword_t *hk,  TIME *ptime)
{
  /* variables */
  int sec;
  int subsec;
  HK_Keyword_t *t_hk;
  int SEC_FOUND_FLAG=0;
  int SUBSEC_FOUND_FLAG=0;
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


/*************************************************************************
 * Name:Lookup data series name Function                                 *
 * Description: Looks up data series JSOC version number and returns     *
 *              value. Looks up values in JSOC version number map files. *
 * Function: lookup_dsn(char* pvn, int apid, char *pvn, int cnt)         *
 *************************************************************************/
char* lookup_dsn( int apid, char *pvn, char *fvn, int cnt)
{
  /* variables */
  int i;
  int found_jvn=0;

  /* lookup jsvn associated with pvn */
  for(i=0; i < cnt ;i++)
  {

    if(check_for_sdo_apid(apid))
    {

      /* check if found pvn in packet */
      if (!strcmp (jmap->fvn[i], fvn ))
      {
        found_jvn = 1;
        break;
      }
    }
    else
    {
      /* check if found pvn in packet */
      if (!strcmp (jmap->pvn[i], pvn ))
      {
        found_jvn = 1;
        break;
      }
    }
  }
  if ( found_jvn != 1)
  {
    if(check_for_sdo_apid(apid))
    {
      printkerr("ERROR at %s, line %d: Did not find the data series name "
                "for this file version number. Maybe this is new configuration "
                "and a new JVN-TO-PVN map files need to updated using script.\n",
                __FILE__,__LINE__);
      return ((char *)"");
    }
    else
    {
      printkerr("ERROR at %s, line %d: Did not find the data series name "
                "for this packet version number. Maybe this is new configuration "
                "and a new JVN-TO-PVN map files need to updated using script.\n",
                __FILE__,__LINE__);
      return ((char *)"");
    }
  }
  return ((char *)jmap->dsn[i]);
}



/**************************************************************************
 * PRINT PACKET VALUES 1                                                  *
 * FUNCTION: saveprint_packet_values1(unsigned char [])                   *
 * DESCRIPTION: Prints translated and decoded packet values and writes    *
 * values to DRMS using the data series name passed in argument list      *
 **************************************************************************/
void saveprint_packet_values1(unsigned char *read_in_buffer, char *in, char *out, char *src)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  char file_version_number[HKDDF_MAX_PVN_SIZE];
  char data_seriesname[HKDDF_MAX_DSNAME_STR];
  char *ptr_fvn;
  int i,j,k,y,x,s;
  int apid;
  int factor;
  int packet_length;
  int pvnw,pvnd;
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;
  char pkt_date[MAX_SIZE_PKT_DATE]; //ascii time

  /* go thru each packet and print to standout and save to drms */
  for(k=0,factor=0,packet_length=0;  ; k++)
  {

    /* init dsn */
    for(x=0; x < HKDDF_MAX_DSNAME_STR; data_seriesname[x]='\0', x++);

    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;

    /* Check if at end of all pkts */
    if (*(read_in_buffer+5+factor) == '\0')
    {
       /* printf for report information */
       printf("\n\n*****At end of file. Found %d packets in file.\n", k);
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
    /* if ADP type apid, packet version number not used */
    if(check_for_sdo_apid(apid))
    {
       sprintf(packet_version_number,"%s","not applicable");
    }
    else
    {
      /*check if packet version is 0.0 -skip-print warning*/
      if (( *(read_in_buffer+14+factor) == 0) && (*(read_in_buffer+15+factor) == 0))
      {
        printkerr("Warning at %s, line %d: Getting 0.0 for packet version number "
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
      strcat(packet_version_number,"\0");/*added 10-2-2008*/
    }

    /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;

    /* set buffer to values in packet and set packet version number to value in packet or some default value*/
    for (i =0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
        /* set values in array */
        hk_pkt_buffer[j++]= *(read_in_buffer+i);
    } /* end for loop setting values for one extracted packet in buffer */


    /* format bytes in format required for decode_hk_keyword() function to process */
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      s_hk_pkt_buffer[y] = (unsigned short)(hk_pkt_buffer[i + 1] << 8  & 0xFF00 );
      s_hk_pkt_buffer[y] = (unsigned short)((hk_pkt_buffer[i] & 0x00FF) + s_hk_pkt_buffer[y]) ;
    }

    /* send hk_pkt_buffer to decoder function and check status returned */
    word_ptr = s_hk_pkt_buffer;
    s = decode_hk_keywords(word_ptr, apid,  &ccsds.keywords);
    if (s)
    {
      printkerr("ERROR at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n\n", __FILE__,__LINE__, s);
      continue;/*added-10-7-2008*/
      
    }
    kw_head= ccsds.keywords;

    /* if sdo hk packet then lookup file version number using packet date */
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
      /* get fvn for hk apids add 10-6-2008*/
      ptr_fvn=find_file_version_number(global_gtcids_vn, packet_version_number,apid); /*added 10-6-2008*/
    }

    /* get data name by automatically creating JSVN value based on */
    /* packet version number                                       */
    strcpy(file_version_number,ptr_fvn);
    strcpy(data_seriesname, out);


    /* print packet keywords */
    printf ( "------------------------------------------------------------------------------\n");
    printf("Packet Number in day file = %-3d       Packet Lenght          = %d \n", k, packet_length);
    if(check_for_sdo_apid(apid))
    {
      printf("Apid                      = %-3.3x       File Version Number  = %s\n", apid,file_version_number);
    }
    else
    {
      printf("Apid                      = %-3.3x       Packet Version Number  = %s\n", apid,packet_version_number);
    }
    printf("In Day File               = %s\n", in);
    printf("Out Series Name           = %s\n", out);
    printf ( "------------------------------------------------------------------------------\n");
    printf("Packet Values in unsigned shorts: \n"  );
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      /*print information contained in packet being sent to decode_hk function */
      printf("%2.2x",  s_hk_pkt_buffer[y]  & 0x00FF);
      printf("%2.2x ",  s_hk_pkt_buffer[y] >> 8 & 0x00FF);
      if ( i ==  14 || i == 30  || i == 46  ||
           i == 62  || i == 78  || i == 94  ||
           i == 110 || i == 126 || i == 142 ||
           i == 158 || i == 174 || i == 190 ||
           i == 206 || i == 222 || i == 238 ||
           i == 254 || i == 270 || i == 286 )
        printf("\n");
    }
    printf("\n");
 
    /* start printing report on keywords */
    printf ( "-----------------------------------------------------------------------------\n");
    printf ("LONG TELEM MNEMONIC NAME           SHORT    RAW VALUE  ENGR VALUE  ENGR VALUE\n");
    printf ("                                   KEYWORD    (HEX)     (DECIMAL)    (HEX)   \n");

    /* loop throught all keywords in single packet*/
    for ( i=0; ccsds.keywords ; i++)
    {
      printf ( "------------------------------------------------------------------------------\n");
      printf ( "%-35.35s", ccsds.keywords->name);
      printf ( "%-9.8s",  ccsds.keywords->fitsname);
      printf ( "%-8.8lx",  (long int )ccsds.keywords->raw_value);
      

      if( ccsds.keywords->eng_type == KW_TYPE_UINT8)
      {
	printf("    %-12d %-10.2x\n",
	ccsds.keywords->eng_value.uint8_val, ccsds.keywords->eng_value.uint8_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_UINT16)
      {
	printf("    %-12hu %-10.4x\n",
	ccsds.keywords->eng_value.uint16_val, ccsds.keywords->eng_value.uint16_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_UINT32)
      {
	printf("    %-12u %-10.8x\n",
	ccsds.keywords->eng_value.uint32_val, ccsds.keywords->eng_value.uint32_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_STRING)
      {
	printf("    %-12s %-10s\n", ccsds.keywords->eng_value.string_val, ccsds.keywords->eng_value.string_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_DOUBLE)
      {
	printf("    %-12lf \n", ccsds.keywords->eng_value.double_val );
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_INT8)
      {
        printf("    %-12hhd %-.2x\n", ccsds.keywords->eng_value.int8_val, ccsds.keywords->eng_value.int8_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_INT16)
      {
        printf("    %-12d %-10.4x\n", ccsds.keywords->eng_value.int16_val, ccsds.keywords->eng_value.int16_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_INT32)
      {
        printf("    %-12d %-10.8x\n", ccsds.keywords->eng_value.int32_val, ccsds.keywords->eng_value.int32_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_FLOAT)
      {
        printf("    %-20.20f \n", ccsds.keywords->eng_value.float_val );
      }
      else
      {
        printkerr("Warning at %s, line %d: Found Unknown KW_TYPE for this keyword:"
                  "<%d>. Check if new or not valid KW_TYPE has been created in HKPDF files"
                  "Skipping printing keyword. Keyword is <%s>.",
                   __FILE__,__LINE__, (int)ccsds.keywords->eng_type, ccsds.keywords->name);
      }

       /*get next packet */
       ccsds.keywords= ccsds.keywords->next;

    }/*end of loop throught all keywords in single packet*/ 

    /* write values to DRMS with data series name passed as argument*/
    write_to_drms( apid, packet_version_number, file_version_number, out, kw_head, in, src );
    /*now print packet time to report*/
    print_packet_time(kw_head);

    /* check if using merged packet version number- if is merged,do not print HK_SOURCE keyword since does not exist */
    sscanf(packet_version_number,"%3d.%3d",&pvnw, &pvnd);
    if(check_for_sdo_apid(apid))
    {
      print_hk_src_keyword(in,apid,src);
    }
    else if(pvnw  ==  HKDDF_START_MERGED_PVNW && pvnd >= HKDDF_START_MERGED_PVND)
    {
      print_hk_src_keyword(in,apid,src);
    }
    //else skip doing

  }/* loop throught next packet and do all packets in file*/
}/*saveprint_packet_values1*/


/**************************************************************************
 * Get Data Series name                                                   *
 * FUNCTION:  get_ds_pjname                                               *
 * DESCRIPTION: Gets data series project name based on apid  using        *
 *              values loaded in array. c                                 *
 **************************************************************************/
char* get_ds_pjname(int apid, char pjn[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME])
{
  char *pn;

  /* get project name */
 if(apid <= HK_HSB_HIGHEST_HMI_APID && apid >= HK_HSB_LOWEST_HMI_APID  || (apid <= HK_LR_HIGHEST_HMI_APID  && apid >= HK_LR_LOWEST_HMI_APID))
  {
    pn = pjn[pjnHMI];
  }
  else if((apid <= HK_HSB_HIGHEST_AIA_APID && apid >= HK_HSB_LOWEST_AIA_APID) || (apid <= HK_LR_HIGHEST_AIA_APID && apid >= HK_LR_LOWEST_AIA_APID))
  {
    pn = pjn[pjnAIA];
  }
  else if((apid <= HK_LR_HIGHEST_SDO_APID) && (apid >= HK_LR_LOWEST_SDO_APID))
  {
    pn = pjn[pjnSDO];
  }
  else
  {
    printkerr("Warning at %s, line %d: APID is not in range of ",
              "HMI, AIA, or SDO. APID: <%d>\n",
               __FILE__, __LINE__, apid);
    pn="\0";
  }
  return((char*)pn);
}


/**************************************************************************
 * Load Data Series Names                                                 *
 * FUNCTION:  load_ds_name                                                *
 * DESCRIPTION: Loads data series names in array with packet version in   *
 *              corresponding array.                                      *
 **************************************************************************/
int load_ds_names(char* pn, char* didn, int apid, char *pvn)
{
  /* variables */
  FILE *file_ptr;
  char *directory ;
  char directory_filename[HKDDF_MAX_FILE_NAME];
  char fn[HKDDF_MAX_FILE_NAME];
  char lookup_name[HKDDF_MAX_LOOKUP_NAME];
  char sn[HKDDF_MAX_FILE_NAME];
  char line[HKDDF_MAXLINE_IN_FILE];
  int count;
  int fwn,fdn, pwn, pdn;
  int pvnw,pvnd;

  /* initialize variables */
  char *filename=fn;
  char *suffix_filename=sn;
  int k=0;

  /* get directory and file name  with jsoc version numbers*/
  directory = (char *)getenv("HK_JSVNMAP_DIRECTORY");
  suffix_filename= (char *)getenv("HK_SUFFIX_JSVNMAP_FILENAME");
  if(suffix_filename == NULL)
  {
    printkerr("ERROR at %s, line %d: Set environment variable <HK_SUFFIX_JSVNMAP_FILENAME>\n",
                __FILE__,__LINE__);
  }
  if(directory == NULL)
  {
    printkerr("ERROR at %s, line %d: Set environment variable <HK_JSVNMAP_DIRECTORY>\n",
                __FILE__,__LINE__);
  }

  /* make filename to use to lookup JSOC Version Number */
  /* if sdo apid just use merged lookup file format */
  if(check_for_sdo_apid(apid))
  {
    /* new for merged jsd -make filename to use to lookup JSOC Version Number */
    strcpy(lookup_name, get_lookup_filename(apid));
    sprintf(filename,"%s%s",lookup_name,suffix_filename);
  }
  else
  { /* else if aia or hmi check packet version to decide lookup file format to use */
    /* check if using merged lookup filename - used merged if for example the pvn=1.197 then fvn=1.163 */
    sscanf(pvn,"%3d.%3d",&pvnw, &pvnd);
    if(pvnw  <=  HKDDF_START_MERGED_PVNW && pvnd < HKDDF_START_MERGED_PVND)
    {
      sprintf(filename,"%4.4d%s",apid,suffix_filename);
    }
    else
    {
      /* new for merged jsd -make filename to use to lookup JSOC Version Number */
      strcpy(lookup_name, get_lookup_filename(apid));
      sprintf(filename,"%s%s",lookup_name,suffix_filename);
    }
  }

  /* put together directory and filename */
  sprintf(directory_filename, "%s/%s", directory, filename);

  /* get apid */
  sprintf(jmap->apid,"%d",apid);

  /* open file and read */
  file_ptr=fopen(directory_filename,"r");
  if (!file_ptr)
  {
    printkerr("ERROR at %s, line %d: Check if directory and filename are",
              "correct. Could not open filename: <%s>\n",
               __FILE__, __LINE__, directory_filename);
    return(0);
  }
  /* get lines in JSVN-TO-PVN files */
  k=0;
  while( fgets(line, HKDDF_MAXLINE_IN_FILE, file_ptr) != NULL )
  {
    if(line[0] == '#')
    {
      continue; /* skip comments */
    }
    else
    {
      sscanf( line,
        "%s | %d.%d | %d.%d | %*s", jmap->jvn[k], &pwn,&pdn,&fwn,&fdn);
      sprintf(jmap->pvn[k],"%03d.%03d", pwn,pdn);

      /*process file version number without leading zero like is done for packet version number*/
      sprintf(jmap->fvn[k],"%d.%d", fwn,fdn);

      /* if jsvn is 0 make dsn null else set data series name in array*/
      if (!strcmp("0000",  jmap->jvn[k]))
      {
        strcpy(jmap->dsn[k], "\0");
      }
      else
      {
        /* check for sdo apid -if so then use merged name for output data series */
        if(check_for_sdo_apid(apid))
        {
          sprintf(jmap->dsn[k], "%s.%s_%s_%s", pn, didn, get_data_packet_name(apid), jmap->jvn[k]);
        }
        else
        {
          /* check if using merged data series name - used merged if pvn >= 1.197 then fvn >= 1.163 */
          if(pwn  <=  HKDDF_START_MERGED_PVNW && pdn < HKDDF_START_MERGED_PVND)
          {
            sprintf(jmap->dsn[k], "%s.%s_%04d_%s", pn, didn, apid, jmap->jvn[k]);
          }
          else
          {
            sprintf(jmap->dsn[k], "%s.%s_%s_%s", pn, didn, get_data_packet_name(apid), jmap->jvn[k]);
          }
        }
      }
    }
    k++;
  }
  fclose( file_ptr);
  count=k;
  return (count);
}



/**************************************************************************
 * SAVE PACKET VALUES 2                                                   *
 * FUNCTION: save_packet_values2()                                        *
 * DESCRIPTION: Translates and decoded packet values and writes           *
 * values to DRMS. Also automatically creates data series name  based on  *
 * the packet version number in packet                                    *
 * STATUS: in progress                                                    *
 **************************************************************************/
void  save_packet_values2(unsigned char *read_in_buffer, char projname[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME],char *didn, char *in, char *src)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  char file_version_number[HKDDF_MAX_PVN_SIZE];
  char data_seriesname[HKDDF_MAX_DSNAME_STR];
  char *dname;
  char *ptr_fvn;
  int i,j,k,y,x,s;
  int apid;
  int count;
  int factor;
  int packet_length;
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;
  char pkt_date[MAX_SIZE_PKT_DATE]; //ascii time

  /* go thru each packet and save to DRMS */
  for(k=0,factor=0, packet_length=0;  ; k++)
  {
    /* init dsn */
    for(x=0; x < HKDDF_MAX_DSNAME_STR; data_seriesname[x]='\0', x++);
    
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;

    /* Check if at end of all pkts */
    if (*(read_in_buffer+5+factor) == '\0')
    {
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
    /* if ADP type apid, packet version number not used */ 
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
      strcat(packet_version_number,"\0");/*added 10-2-2008*/

    }

    /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
    /* set buffer to values in packet and set packet version number to value 
    /* in packet or some default value */
    for (i =0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
	/* set values in array */
	hk_pkt_buffer[j++]= *(read_in_buffer+i);
    } /* end for loop setting values for one extracted packet in buffer */

    /* put in format for decode_hk to read by adjusting values in buffer  */
    /* from unsigned char to unsigned short and write bytes to report     */
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
 /* TAKE OUT TEMP CC 2-4-2009
 */
      printkerr("ERROR at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n\n", __FILE__,__LINE__, s);
      continue;/*added-10-7-2008*/
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
      /* get fvn for hk apids add 10-2-2008*/
      ptr_fvn=find_file_version_number(global_gtcids_vn, packet_version_number,apid); /*added 10-2-2008*/
    }
    strcpy(file_version_number,ptr_fvn);/*move to here -10-2-2008*/

    /* get data name by automatically creating JSVN value based on */
    /* packet version number                                       */
    count = load_ds_names(get_ds_pjname(apid,projname), didn, apid, packet_version_number);
    dname = lookup_dsn( apid, packet_version_number, file_version_number, count);
    strcpy(data_seriesname, dname);

    /* write values to DRMS with auto generated data series name*/
    write_to_drms( apid, packet_version_number, file_version_number, data_seriesname, kw_head, in,src );

  }/* loop throught next packet and do all packets in file*/
}/*save_packet_values2 */


/**************************************************************************
 * PRINT PACKET VALUES 2                                                  *
 * FUNCTION: saveprint_packet_values2()                                   *
 * DESCRIPTION: Prints translated and decoded packet values and writes    *
 * values to DRMS. Also automatically creates data series name  based on  *
 * the packet version number in packet                                    *
 **************************************************************************/
void  saveprint_packet_values2(unsigned char *read_in_buffer, char *in,
                               char projname[HKDDF_MAXNUM_PJNAMES][HKDDF_MAX_PJNAME],
                               char *didn, char *src)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  char file_version_number[HKDDF_MAX_PVN_SIZE];
  char data_seriesname[HKDDF_MAX_DSNAME_STR];
  char *dname;
  char *ptr_fvn;
  int i,j,k,y,x,s;
  int apid;
  int count;
  int factor;
  int packet_length;
  int pvnw,pvnd;
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;
  char pkt_date[MAX_SIZE_PKT_DATE]; //ascii time

  /* go thru each packet and print to standout and save to drms */
  for(k=0,factor=0, packet_length=0;  ; k++)
  {
    /* init dsn */
    for(x=0; x < HKDDF_MAX_DSNAME_STR; data_seriesname[x]='\0', x++);
    
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;

    /* Check if at end of all pkts */
    if (*(read_in_buffer+5+factor) == '\0')
    {
       printf("\n\n*****At end of file. Found %d packets in file.\n", k);
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
    /* if ADP type apid, packet version number not used */
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
                  "DRMS and get next packet\n", __FILE__,__LINE__);
        continue;
      }
      /* check for version number set to zero and print warning messages. */
      if ( *(read_in_buffer+14+factor) == 0 && *(read_in_buffer+15+factor))
      {
        printkerr("Warning at %s, line %d: Getting 0 for whole number(i.e.,0.1) for "
                  "packet version number in packet data. Skip processing this packet,"
                  "don't load to DRMS and get next packet\n", __FILE__,__LINE__);
        continue;
      }
      /* if passed two tests above, then set packet version number */
      sprintf(packet_version_number,"%03d.%03d",*(read_in_buffer+14+factor), *(read_in_buffer+15+factor));
      strcat(packet_version_number,"\0");
    }

    /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
    /* set buffer to values in packet and set packet version number to value 
    /* in packet or some default value */
    for (i =0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
	/* set values in array */
	hk_pkt_buffer[j++]= *(read_in_buffer+i);
    } /* end for loop setting values for one extracted packet in buffer */


    /* format bytes in format required for decode_hk_keyword() function to process */ 
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      s_hk_pkt_buffer[y] = (unsigned short)(hk_pkt_buffer[i + 1] << 8  & 0xFF00 );
      s_hk_pkt_buffer[y] = (unsigned short)((hk_pkt_buffer[i] & 0x00FF) + s_hk_pkt_buffer[y]) ;
    }

    /* send hk_pkt_buffer to decoder function and check status returned */
    word_ptr = s_hk_pkt_buffer;
    s = decode_hk_keywords(word_ptr, apid,  &ccsds.keywords);
    if (s) 
    {
      printkerr("ERROR at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n\n", __FILE__,__LINE__, s);
      continue;/*added-10-7-2008*/
    }
    kw_head= ccsds.keywords;

    /* find file version number using packet date for sdo hk packets */
    if(check_for_sdo_apid(apid))
    {
      /* after load of shcids file during call to decode_hk_keywords -
         get packet time to lookup sdo hk config file version to use */
      (void)sprint_time (pkt_date, get_packet_time(word_ptr), "TAI", 0);
      /* get file version number for sdo-hk apids */
      ptr_fvn=find_fvn_from_shcids(global_shcids_vn, pkt_date,apid);
    }
    else
    {
      /* get fvn for hk apids */
      ptr_fvn=find_file_version_number(global_gtcids_vn, packet_version_number,apid);
    }
    strcpy(file_version_number,ptr_fvn);

    /* get data name by automatically creating JSVN value based on */
    /* packet version number                                       */
    count = load_ds_names(get_ds_pjname(apid,projname), didn, apid,  packet_version_number);
    dname = lookup_dsn( apid, packet_version_number, file_version_number, count);
    strcpy(data_seriesname, dname);

    /* Move Print stuff here after call to decode_hk()!! */
    /* print packet keywords */
    printf ( "------------------------------------------------------------------------------\n");
    printf("Packet Number in day file = %-3d       Packet Lenght          = %d \n", k, packet_length);
    if(check_for_sdo_apid(apid))
    {
      printf("Apid                      = %-3.3x       File Version Number  = %s\n", apid,file_version_number);
    }
    else
    {
      printf("Apid                      = %-3.3x       Packet Version Number  = %s\n", apid,packet_version_number);
    }
    printf("In Day File               = %s\n", in);
    printf("Out Series Name           = %s\n", data_seriesname);
    printf ( "------------------------------------------------------------------------------\n");
    printf("Packet Values in unsigned shorts: \n"  );
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      /*print information contained in packet being sent to decode_hk function */
      printf("%2.2x",  s_hk_pkt_buffer[y]  & 0x00FF); 
      printf("%2.2x ",  s_hk_pkt_buffer[y] >> 8 & 0x00FF); 
      if ( i ==  14 || i == 30  || i == 46  || 
           i == 62  || i == 78  || i == 94  || 
           i == 110 || i == 126 || i == 142 ||
           i == 158 || i == 174 || i == 190 ||
           i == 206 || i == 222 || i == 238 ||
           i == 254 || i == 270 || i == 286 )
        printf("\n");
    }
    printf("\n");

    /* start print report on keywords */
    printf ( "-----------------------------------------------------------------------------\n");
    printf ("LONG TELEM MNEMONIC NAME           SHORT    RAW VALUE  ENGR VALUE  ENGR VALUE\n");
    printf ("                                   KEYWORD    (HEX)     (DECIMAL)    (HEX)   \n");

    /* loop throught all keywords in single packet*/
    for ( i=0; ccsds.keywords ; i++)
    {
      printf ( "------------------------------------------------------------------------------\n");
      printf ( "%-35.35s", ccsds.keywords->name);
      printf ( "%-9.8s",  ccsds.keywords->fitsname);
      printf ( "%-8.8lx", (long int) ccsds.keywords->raw_value);

      if( ccsds.keywords->eng_type == KW_TYPE_UINT8)
      {
	printf("    %-12d %-10.2x\n",
	ccsds.keywords->eng_value.uint8_val, ccsds.keywords->eng_value.uint8_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_UINT16)
      {
	printf("    %-12hu %-10.4x\n",
	ccsds.keywords->eng_value.uint16_val, ccsds.keywords->eng_value.uint16_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_UINT32)
      {
	printf("    %-12u %-10.8x\n",
	ccsds.keywords->eng_value.uint32_val, ccsds.keywords->eng_value.uint32_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_STRING)
      {
	printf("    %-12s %-10s\n", ccsds.keywords->eng_value.string_val, ccsds.keywords->eng_value.string_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_DOUBLE)
      {
	printf("    %-12lf \n", ccsds.keywords->eng_value.double_val );
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_INT8)
      {
        printf("    %-12hhd %-.2x\n", ccsds.keywords->eng_value.int8_val, ccsds.keywords->eng_value.int8_val);
        /* did update of printf statement on 12-18-2006 */
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_INT16)
      {
        printf("    %-12d %-10.4x\n", ccsds.keywords->eng_value.int16_val, ccsds.keywords->eng_value.int16_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_INT32)
      {
        printf("    %-12d %-10.8x\n", ccsds.keywords->eng_value.int32_val, ccsds.keywords->eng_value.int32_val);
      }
      else if(ccsds.keywords->eng_type == KW_TYPE_FLOAT)
      {
        printf("    %-20.20f \n", ccsds.keywords->eng_value.float_val);
      }
      else
      {
        printf("WARNING(decode_dayfile):****Found*** NO KW_TYPE\n");
        printf("Eng type is %d\n", (int)ccsds.keywords->eng_type);
      }

       /*get next packet */
       ccsds.keywords= ccsds.keywords->next;

    }/*end of loop throught all keywords in single packet*/ 


    /* write values to DRMS with auto generated data series name*/
    write_to_drms( apid, packet_version_number, file_version_number, data_seriesname, kw_head, in, src );
    /*now print packet time to report*/
    print_packet_time(kw_head);

    /* check if using merged packet version number- if is merged, do not print HK_SOURCE keyword since does not exist */
    sscanf(packet_version_number,"%3d.%3d",&pvnw, &pvnd);
    if(check_for_sdo_apid(apid))
    {
      print_hk_src_keyword(in,apid,src);
    }
    else if(pvnw  ==  HKDDF_START_MERGED_PVNW && pvnd >= HKDDF_START_MERGED_PVND)
    {
      print_hk_src_keyword(in,apid,src);
    }
    //else skip doing
    
  }/* loop throught next packet and do all packets in file*/
} /*  saveprint_packet_values2 */


/**************************************************************************
 * PRINT PACKET TIME                                                      *
 * FUNCTION: print_packet_time()                                          *
 * DESCRIPTION: Get values and prints keyword PACKET_TIME  to standard out*
 * STATUS: completed                                                      *
 **************************************************************************/
void print_packet_time(HK_Keyword_t *kwh)
{
  HK_Keyword_t *kw;
  TIME pkt_time;
  kw= kwh;
  pkt_time=(TIME)0.0;
  /* get and set packet_time in record*/
  if (!get_packet_time_for_df(kw, &pkt_time))
  {
    /*did not set PACKET_TIME */
    printf ( "------------------------------------------------------------------------------\n");
    printf("ERROR: COULD NOT SET PACKET_TIME KEYWORD\n");
    return;
  }
  else 
  { 
    printf ( "------------------------------------------------------------------------------\n");
    printf("PACKET_TIME                                     %15.5lf\n", pkt_time);
  }
  return;
}


/**************************************************************************
 * SAVE PACKET VALUES 1                                                   *
 * FUNCTION: save_packet_values1()                                        *
 * DESCRIPTION: Translates and decoded packet values and writes           *
 * values to DRMS. Creates data series name using argument passed at      *
 * command line.                                                          *
 * STATUS:completed                                                       *
 **************************************************************************/
void save_packet_values1(unsigned char *read_in_buffer, char *out, char *in, char *src)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  char file_version_number[HKDDF_MAX_PVN_SIZE];
  int i,j,k,y,s;
  int apid;
  int factor;
  int packet_length;
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;
  char pkt_date[MAX_SIZE_PKT_DATE]; //ascii time
  char *ptr_fvn;

  /* go thru each packet and save keywords to DRMS */
  for(k=0,factor=0,packet_length=0;  ; k++)
  {
    //printf("Packet=%d\n",k);
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;

    /* Check if at end of all pkts */
    if (*(read_in_buffer+5+factor) == '\0')
    {
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
    /* if sdo-hk type apid, packet version number not used, else get packet version number for hmi or aia hk packets */
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
      strcat(packet_version_number,"\0");/*added 10-2-2008*/
    }

     /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
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
      continue;/*added-10-7-2008*/
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
      /* get fvn for hk apids add 10-6-2008*/
      ptr_fvn=find_file_version_number(global_gtcids_vn, packet_version_number,apid); /*added 10-6-2008*/
    }
    /* set file version number  */
    strcpy(file_version_number,ptr_fvn);

    /* write values to DRMS with data series name from command line argument*/
    write_to_drms( apid, packet_version_number, file_version_number, out, kw_head, in,src );

  } /* loop throught next packet and do all packets in file*/
}/*end of save_packet_values1*/


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
      printf("ERROR:Can't open environment variable file <%s>. Check setting is correct\n", envfile);
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



/************************************************************/
/* NAME:Get HK_SOURCE Value to Set in Lev0 by APID Series   */
/* FUNCTION:   get_hk_src_keyword()                         */
/* DESCRIPTION:Get value to set HK_SOURCE keyword in Lev0   */ 
/*             by APID Series by using PACKET_TIME to set   */
/*             Date indx value, Apid to set apid index and  */
/*             hsb to source index value. The project name  */
/*             to use is based on choosing three environment*/
/*             HK_LEV0_DF_MAX_PROJECT_NAME_xxx where xxx is */
/*             either SDO,HMI or AIA.                       */
/************************************************************/
char *get_hk_src_keyword(int apid, char *in_file, char *src)
{
  /* declarations */
  char date[HK_LEV0_MAX_HKS_DATE];
  char df_datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char df_project_name[HK_LEV0_MAX_PROJECT_NAME];
  char *df_dtname, *df_pjname;
  char hk_source[HK_LEV0_MAX_HKS];
  char *hks;
  int  len_hk_source;
  int  pkt_yr, pkt_month, pkt_day;
  int status;

  /* initial pointers */
  df_dtname=df_datatype_name;
  df_pjname=df_project_name;
  hks = hk_source;

  /* get dayfile project name */
  if(apid <= HK_HSB_HIGHEST_HMI_APID && apid >= HK_HSB_LOWEST_HMI_APID  || (apid <= HK_LR_HIGHEST_HMI_APID  && apid >= HK_LR_LOWEST_HMI_APID))
  {
    df_pjname = getenv("HK_DDF_DF_PROJECT_NAME_HMI");
  }
  else if((apid <= HK_HSB_HIGHEST_AIA_APID && apid >= HK_HSB_LOWEST_AIA_APID) || (apid <= HK_LR_HIGHEST_AIA_APID && apid >= HK_LR_LOWEST_AIA_APID))
  {
    df_pjname = getenv("HK_DDF_DF_PROJECT_NAME_AIA");
  }
  else if((apid <= HK_LR_HIGHEST_SDO_APID) && (apid >= HK_LR_LOWEST_SDO_APID))
  {
    df_pjname = getenv("HK_DDF_DF_PROJECT_NAME_SDO");
  }
  else
  {
    printkerr("Warning at %s, line %d: APID is not in range of ",
              "HMI, AIA, or SDO. APID: <%d>\n",
               __FILE__, __LINE__, apid);
  }

  /* get dayfile data type name */
  df_dtname = getenv("HK_DDF_DF_DATA_ID_NAME");

  /* get date from current packet to know date to use for index to dayfile series */
  /* parse in filename to get year, month and day */
  status=get_date_in_filename(in_file, src, &pkt_yr, &pkt_month, &pkt_day);
  if(!status)
  {
    printkerr("ERROR at %s, line %d: HK_SOURCE keyword value could not be set. "
              "Setting hk source keyword value to null for input filename <%s>  "
              "with source value <%s>. \n", __FILE__, __LINE__, in_file, src);
  }

  /* Put dayfile time in this format like this: 2007.11.08_00:00:00.00_TAI */
  sprintf(date,"%d.%-02.2d.%-02.2d_00:00:00.00_TAI",pkt_yr, pkt_month, pkt_day);

  /* check length of HK_SOURCE keyword does not exceed max allowed value */
  len_hk_source=strlen(df_pjname) + strlen(df_dtname) + strlen(date) + strlen("hsb") + strlen("475");
  if(len_hk_source >= HK_LEV0_MAX_HKS - 1)
  {
    printkerr("ERROR at %s, line %d: HK_SOURCE keyword value exceeds maximum size of array. "
              "Length is %d. Setting source value to null so that we do not exceed hks array.\n",
               __FILE__, __LINE__, len_hk_source);
    strcpy(hks,"" );
  }
  else
  {
    /* create final string for hk source using format:hmi.hk_dayfile[2007.11.08_00:00:00.00_TAI][475][hsb] */
    sprintf(hks,"%s.%s[%s][%d][%s]\0", df_pjname, df_dtname, date, apid, src);
  }
  return(hks);
}


/**************************************************************************/
/* NAME:GET DATE IN FILENAME                                              */
/* FUNCTION:   int get_date_in_filename(char *inf,char *src, int *yr,     */
/*                                      int *mo, int *day)                */
/* DESCRIPTION: Get date in filename for date value in HK_SOURCE keyword. */
/*              Pass input filename and source(src) to determine file     */
/*              format. Pass back year, month and day.                    */
/**************************************************************************/
int get_date_in_filename(char *inf,char *src, int *yr,int *mo, int *day)
{
  /*declarations*/
  char fn_str[HKDDF_MAX_FILE_NAME];
  char final_fn_str[HKDDF_MAX_FILE_NAME];
  char *fn,*final_fn;
  int ret_status, year_day;

  /* intialize variables */
  fn=fn_str;
  final_fn=final_fn_str;
  ret_status=1;

  /* move to beginning of filename */
  fn= strrchr(inf, '/');

  /* move to beginning of filename without slash mark */
  final_fn=strtok(fn,"/");
  final_fn=strtok(fn,"/");

  if(!strcmp(src,"moc"))
  {
    /* process moc formatted files:0129_2008_191_01.hkt */
    /* get year and year day */
    sscanf(final_fn,"%*4s_%4d_%3d_%*2d.%*3s",yr,&year_day);
    
    /* get month and day based on year day */
    (void)get_month_day(*yr, year_day, mo, day);

  }
  else if(!strcmp(src,"rtmon"))
  {
    /* process rtmon formatted files:20080101.0x01d */
     sscanf(final_fn,"%4d%2d%2d",yr,mo,day);
  }
  else if(!strcmp(src,"hsb"))
  {
    /* process hsb formatted files:hsb_0445_2008_04_15_10_41_00.hkt */
    sscanf(final_fn,"%*3s_%*4s_%4d_%2d_%2d_%*8s.%*3s",yr,mo,day);

  }
  else if(!strcmp(src,"egsefm"))
  {
    /* process egsefm formatted files:20080101.0x01d */
     sscanf(final_fn,"%4d%2d%2d",yr,mo,day);
  }
  else
  {
    printkerr("ERROR at %s, line %d:  Unknown source value. Unexpected value:<%s>\n",
               __FILE__, __LINE__, src);
    ret_status=0;
  }
  return(ret_status);
}

/**************************************************************************
 * PRINT HK SOURCE KEYWORD                                                *
 * FUNCTION: print_hk_src_keyword()                                       *
 * DESCRIPTION: Get values and prints keyword HK_SOURCE  to standard out  *
 **************************************************************************/
void print_hk_src_keyword(char *ifn, int apid, char *src)
{
  char hk_source[HKDDF_MAX_HK_SRC_KW];

  /* get hk source keyword value */
  strcpy(hk_source, get_hk_src_keyword(apid, ifn, src));
  if (!hk_source)
  {
    /*did not set PACKET_TIME */
    printf ( "------------------------------------------------------------------------------\n");
    printf("ERROR: COULD NOT SET PACKET_TIME KEYWORD\n");
    printf ( "------------------------------------------------------------------------------\n");
    printf ("%s\n",HKDDF_CNTRL_L);
    return;
  }
  else 
  { 
    printf ( "------------------------------------------------------------------------------\n");
    printf("HK_SOURCE  %-100s\n", hk_source);
    printf ( "------------------------------------------------------------------------------\n");
    printf ("%s\n",HKDDF_CNTRL_L);
  }
  return;
}

/**************************************************************************/
/* Get Month and Day using Day of Year                                    */
/* FUNCTION: get_month_day(int year, int yearday, int *pmonth, int *pday) */
/* DESCRIPTION: Get values of month and day based on year and year day    */
/**************************************************************************/
void get_month_day(int year, int yearday, int *pmonth, int *pday)
{
  static char daytab[2][13] = 
  {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
  };
  int i, leap /* 1 for leap year, 0 for non-leap */;
  leap = ((year%4 == 0) && year%100 != 0) || year%400 == 0;
  for(i = 1; yearday > daytab[leap][i]; i++)
  {
    yearday -= daytab[leap][i];
  }
  *pmonth = i;
  *pday = yearday;
}
/************************************************************/
/* NAME:Check Environment Variable                          */
/* FUNCTION:   check_env_variable()                         */
/* DESCRIPTION: Checks all environment variables are set    */ 
/*              SOURCE_ENV_HK_DECODE file.                  */
/************************************************************/
int check_df_env_variable()
{
  /* declarations */
  char hkbyapid_project_name[HK_LEV0_MAX_PROJECT_NAME];
  char hkbyapid_datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char df_project_name[HK_LEV0_MAX_PROJECT_NAME];
  char df_datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char *hk_dtname, *hk_pjname, *df_dtname, *df_pjname;

  /* initial pointers */
  hk_pjname=hkbyapid_project_name;
  df_dtname=df_datatype_name;
  df_pjname=df_project_name;
  hk_dtname=hkbyapid_datatype_name;

  /* Check 4 environment variables are set for dayfile series name -Use to set HK_KEYWORD in HK by APID series*/

  /* get project name for HK by APID data series */
  df_pjname = getenv("HK_DDF_DF_PROJECT_NAME_HMI");
  if(df_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_DDF_DF_PROJECT_NAME_HMI>. Set the env variable "
              "HK_DDF_DF_PROJECT_NAME_HMI to dayfile data series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing dayfile data series name.\n",
              __FILE__,__LINE__);
    return  HKDDF_ENVIRONMENT_VARS_NOT_SET ;
  }
  df_pjname = getenv("HK_DDF_DF_PROJECT_NAME_AIA");
  if(df_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_DDF_DF_PROJECT_NAME_AIA>. Set the env variable "
              "HK_DDF_DF_PROJECT_NAME_AIA to dayfile data series project name"
              "(aia, aia_ground, su_carl,etc.) for a existing data series name.\n",
              __FILE__,__LINE__);
    return  HKDDF_ENVIRONMENT_VARS_NOT_SET ;
  }
  df_pjname = getenv("HK_DDF_DF_PROJECT_NAME_SDO");
  if(df_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_DDF_DF_PROJECT_NAME_SDO>. Set the env variable "
              "HK_DDF_DF_PROJECT_NAME_SDO to dayfile data series project name"
              "(sdo, sdo_ground, su_carl,etc.) for a existing data series name.\n",
              __FILE__,__LINE__);
    return  HKDDF_ENVIRONMENT_VARS_NOT_SET ;
  }

  /* get data type name for HK BY APID data series */
  df_dtname = getenv("HK_DDF_DF_DATA_ID_NAME");
  if(df_dtname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get data type name environment "
              "variable:<HK_DDF_DF_DATA_ID_NAME>. Set the env variable "
              "HK_DDF_DF_DATA_ID_NAME to dayfile data series data type name"
              "(hk_dayfile, etc.) for a existing dayfile data series name\n",
              __FILE__ , __LINE__ );
     return (HKDDF_ENVIRONMENT_VARS_NOT_SET) ;
  }

  /* Check 4 environment variables are set- Used for creating HK BY APID output series name to write to*/

  /* get dayfile project name for setting HK_SOURCE value in HK by APID data series */
  hk_pjname = getenv("HK_DDF_HKBYAPID_PROJECT_NAME_HMI");
  if(hk_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_DDF_HKBYAPID_PROJECT_NAME_HMI>. Set the env variable "
              "HK_DDF_HKBYAPID_PROJECT_NAME_HMI to series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing HK by APID series name.\n",
              __FILE__,__LINE__);
    return  HKDDF_ENVIRONMENT_VARS_NOT_SET ;
  }
  hk_pjname = getenv("HK_DDF_HKBYAPID_PROJECT_NAME_AIA");
  if(hk_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_DDF_HKBYAPID_PROJECT_NAME_AIA>. Set the env variable "
              "HK_DDF_HKBYAPID_PROJECT_NAME_AIA to series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing HK by APID series name.\n",
              __FILE__,__LINE__);
    return  HKDDF_ENVIRONMENT_VARS_NOT_SET ;
  }
  hk_pjname = getenv("HK_DDF_HKBYAPID_PROJECT_NAME_SDO");
  if(hk_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_DDF_HKBYAPID_PROJECT_NAME_SDO>. Set the env variable "
              "HK_DDF_HKBYAPID_PROJECT_NAME_SDO to series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing HK By APID series name.\n",
              __FILE__,__LINE__);
    return  HKDDF_ENVIRONMENT_VARS_NOT_SET ;
  }
  
  /* get data type name for HK BY APID data series */
  hk_dtname = getenv("HK_DDF_HKBYAPID_DATA_ID_NAME");
  if(hk_dtname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get data type name environment "
              "variable:<HK_DDF_HKBYAPID_DATA_ID_NAME>. Set the env variable "
              "HK_DDF_HKBYAPID_DATA_ID_NAME to data series data type name"
              "(lev0,lev0test, etc.) for a existing HK By APID data series name\n",
              __FILE__ , __LINE__ );
     return (HKDDF_ENVIRONMENT_VARS_NOT_SET) ;
  }
  /* all environment variable are at least set, so return 1 */
  return (1);
}

