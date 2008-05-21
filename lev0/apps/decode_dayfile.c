/*############################################################################
# Name:        decode_dayfile.c - Decode Dayfile                             #
# Description: Decode dayfile decodes hk packet by sending packet to         #
#              functions in decode_hk.c file. Keywords are decoded using the #
#              the same functions that decodes hk keywords on high speed bus.#
#              This decode_dayfiles loops through all hk packets and sends   #
#              data to be decoded to decode_hk.c and write data to DRMS data #
#              series based on apid, project name, data type name, and jsoc  #
#              version number.When run without out parameter this code looks #
#              up name of data series name to put keyword names and values in#
# Execution:   decode_dayfile [-p] in=<day-file> [ out=<data series name>]   #
#              (1)To create data out series on automatically:                #
#                     decode_dayfile  in=<day-file>                          #
#              (2)To specify data series :                                   #
#                     decode_dayfile  in=<day-file>  out=<data series name>  #
#              (3)To specify out data series and print to standard out :     #
#                  decode_dayfile -p in=<day-file>  out=<data series name>   #
#              (4)create data out series automatically and print    standard #
#                 output to file.                                            #
#                  decode_dayfile -p  in=<day-file>  > View-OUTPUT           #
# Example 1:   decode_dayfile in=/tmp20/20070202.0x001d (Best Option)        #
# Example 2:   decode_dayfile in=/tmp20/20070202.0x001d                      #
#                             out=su_carl.lev0_0029_0001                     #
# Limitation:  Setup required for environment variables in                   #
#              SETENV_HK_DAYFILE_DECODE. Currently the code read packet size #
#              limit of 1000 bytes (HKDDF_MAX_PKT_SIZE ). Need to setup local#
#              version of SETENV_HK_DAYFILE_DECODE file since file checked in#
#              is production version of file.                                #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on April 10, 2008 #
############################################################################*/

/******************** defines ***********************************************/
#define HKDDF_CNTRL_L            ""
#define HKDDF_MAX_FILE_NAME      100
#define HKDDF_MAXLINE_IN_FILE    200
#define HKDDF_MAX_VERSION_LINES  1000
#define HKDDF_MAX_APID_STR       5
#define HKDDF_MAX_JSVN_STR       5
#define HKDDF_MAX_PVN_STR        10
#define HKDDF_MAX_DSNAME_STR     100
#define HKDDF_MAX_PKT_SIZE       1000
#define HKDDF_MAX_PVN_SIZE       50
#define ENVFILE      "/home/production/cvs/JSOC/proj/lev0/apps/SOURCE_ENV_FOR_HK_DAYFILE_DECODE"

/******************** includes ******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include "drms.h"
#include "drms_names.h"
#include "jsoc_main.h"
#include "load_hk_config_files.h"
#include "packets.h"
#include "decode_hk.h"
#include "printk.h"


/************* modules definitions **************************************/
ModuleArgs_t module_args[] =
{
  {ARG_STRING, "in", "Not Specified", "full path to day file"},
  {ARG_STRING, "out", "Not Specified", "Series name"},
  {ARG_FLAG, "p", "0", "print values of  keywords to standard out"},
  {ARG_END}
};
ModuleArgs_t   *ggModArgs=module_args;
char* module_name = "decode_dayfile";

/******************* function prototypes  *******************************/
static int   get_packet_time_for_df(HK_Keyword_t *hk, TIME *ptime);
static int   load_ds_names(char *pn, char* didn, int apid);
static char* lookup_dsn( char *pvn, int count);
static void  print_packet_time(HK_Keyword_t *kwh);
static void  save_packet_values1(unsigned char read_in_buffer[], char *out);
static void  save_packet_values2(unsigned char read_in_buffer[], char *pn, char *didn );
static void  saveprint_packet_values1(unsigned char read_in_buffer[], char *in, char *out);
static void  saveprint_packet_values2(unsigned char read_in_buffer[], char *in, char *pn, char *didn );
static void  set_env_variables();
static void  write_to_drms( char packet_version_number[], char *data_series, HK_Keyword_t *kw_head);
static TIME  SDO_to_DRMS_time(int sdo_s, int sdo_ss);

/********************* extern functions  *********************************/
extern int DoIt(void);
extern int nice_intro (void);

/********************* structures   *************************************/
/*structure used to save values from JSVN-TO-PVN map files */
struct jsvn_map_data
{
  char apid[HKDDF_MAX_APID_STR];
  char jvn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_JSVN_STR];
  char pvn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_PVN_STR];
  char dsn[HKDDF_MAX_VERSION_LINES][HKDDF_MAX_DSNAME_STR];
} jsvn_list[1], *jmap=jsvn_list;



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
    printf ("Usage:\ndecode_dayfile [-p] "
      "in=<day filename> "
      "[ out=<data series name> ] \n"
      "  details are:\n"
      "  -h: help - show this message then exit\n"
      "  -p: print all keyword names and values to standard out(optional field)\n"
      "  in=<day file name>  -use full path to day file(required field).       \n"
      "  out=<data series name>  -data series to write keyword values to in DRMS(optional field). \n"
      "  if don't enter out value, data series name will be created automatically using the apid and \n"
      "  packet version number in packets in day file. Also required to set environment variables \n"
      "  for project name and data identifer name: HK_DDF_PROJECT_NAME & HK_DDF_DATA_ID_NAME.     \n"
      "  Also need data series already created for this to program to load keywords in data series.\n"
      "  Example of running using automatic lookup of data series name and without print on:      \n"
      "           decode_dayfile in=/tmp20/production/hmi_hk/20070707.0x0029                      \n" );
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
  char *pn;
  char *didn;
  int i;
  int out_flag;
  unsigned char read_in_buffer[6588020];
  unsigned char *ptr_read_in_buffer;

  /* set environment variables */
  set_env_variables();

  /* parameter initialization */
  hk_df_fn= hk_directory_filename;
  ptr_read_in_buffer = read_in_buffer;
  out_flag=0;

  /* Get command line arguments */
  int print_flag =cmdparams_get_int (&cmdparams, "p", NULL) != 0;
  char *in = cmdparams_get_str (&cmdparams, "in", NULL);
  char *out = cmdparams_get_str (&cmdparams, "out", NULL);

  /* check arguments used */
  if (nice_intro ()) return (0);

  /* check if entered day file name */
  if (in == NULL) 
  {
    printkerr("Error at %s, line %d: Need to enter day file name.",
              "Exiting program.\nExample format for out file: "
              "$HOME/EGSE/tables/hk_sim_file/20060719.0x000f\n",
               __FILE__,__LINE__);
    return (0);
  }
  if (!strcmp(out, "Not Specified")) 
  {
    /*set out_flag to 0 to trigger creating automatically data series name  */
    out_flag=0;

    /* check if environment variables setup for automatic creation of data  */
    /* series name project name */
    pn = getenv("HK_DDF_PROJECT_NAME");
    didn = getenv("HK_DDF_DATA_ID_NAME");
    if(pn == NULL) 
    {
      printkerr("Error at %s, line %d: Could not get project name environment\n"
                "variable:<HK_DDF_PROJECT_NAME>. Either set the env variable "
                "HK_DDF_PROJECT_NAME  or set out argument to a value.\n"
                "Exiting execution.\n-Example format for using out argument:"
                "<project-name>.<data specifier>_<4 digit apid>_<4 digit jsvn number>\n"
                "-Example name to use for out argument: < hmi_ground.lev0_0029_0001>\n"
                "-Example of setting environment variable: setenv HK_DDF_PROJECT_NAME hmi_ground\n",
                 __FILE__,__LINE__, pn);
      return (0);
    }

    /* series name data type name */
    if(didn == NULL) 
    {
      printkerr("Error at %s, line %d: Could not get data ID name environment\n"
                "variable:<HK_DDF_DATA_ID_NAME>. Either set the env variable "
                "HK_DDF_DATA_ID_NAME  or set out argument to a value.\n"
                "Exiting execution.\n-Example format for using out argument:"
                "<project-name>.<data specifier>_<4 digit apid>_<4 digit jsvn number>\n"
                "-Example name to use for out argument: < hmi_ground.lev0_0029_0001>\n"
                "-Example of setting environment variable: setenv HK_DDF_DATA_ID_NAME lev0\n",
                  __FILE__,__LINE__, pn);
      return (0);
    }
  }
  else
  {
    /* if have out argument set by user set out_flag to 1 */
    out_flag=1; 
  }

  /* get  in file name and open file*/
  strcpy(hk_df_fn, in);
  strcat(hk_df_fn, "\0");
  file_ptr=fopen(hk_df_fn,"r");
  if (!file_ptr)
  {
    printkerr("Error at %s, line %d: Please check filename is correct. Could not create filename:\n"
              "<%s>\n", __FILE__,__LINE__, hk_df_fn);
    return(0);
  }

  /*read lines in file into buffer representing a packet*/
  for(i=0; i < 6588020;i++) read_in_buffer[i]=0; ;
  for(i = 0 ; fread(ptr_read_in_buffer,1,1,file_ptr) ; i++,ptr_read_in_buffer++) 
  {  
    ;/* do nothing*/
  }
  ptr_read_in_buffer=NULL;
  fclose(file_ptr);
  
  if(print_flag && out_flag)
  {
    /* start printing packet information and process keywords to drms*/
    /*  and use command line argument for data series name           */
    saveprint_packet_values1(read_in_buffer, in, out);
    return 0;
  }
  else if(print_flag && !out_flag)
  {
    /* start printing packet information and process keywords to drms and */
    /* automatically create JSVN and use as data series name              */
    saveprint_packet_values2(read_in_buffer, in, pn, didn );
    return 0;
  }
  else if (!print_flag && out_flag)
  {
    /* do not print packet information and use command line argument for */
    /* data series name. this is case where: if(out_flag && !print_flag) */
    /* write to drms only - do not print values to standard out          */
    save_packet_values1(read_in_buffer, out);
    return 0;
  }
  else if (!print_flag && !out_flag)
  {
    /* do not print packet information and use automatic lookup of  */
    /* data series name. Environment variables need to be set for   */
    /* project name and data identifier name. write to drms only    */
    /* - do not print values to standard out                        */
    save_packet_values2(read_in_buffer, pn, didn);
    return 0;
  }
  else 
  {
    printkerr("Error at %s, line %d: Save and Print Flags are not set correctly."
              "This case should not occur. Check logic. Exiting Program.\n" ,
               __FILE__,__LINE__ );
    return 1;
  }
}



/*************************************************************************
 * WRITE TO DRMS                                                         *
 * FUNCTION: void write_to_drms(char [],char *,HK_Keyword_t *            *
 * DESCRIPTION: Gets Packet time based on Time Codes.                    *
 *************************************************************************/
void write_to_drms(char pkt_ver_num[],char *ds_name, HK_Keyword_t *kw_head)
{
  /* variable definitions */
  /* create series name variables*/
  char query[HKDDF_MAX_DSNAME_STR];
  int status;
  /* keyword variables and structure */
  HK_Keyword_t *kw;
  kw= kw_head;
  TIME pkt_time;

  /* variable to set drms record */
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  char keyname[50];
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
    printkerr("Error at %s, line %d: Bad value for data series name:"
              "<%s>. Exiting program.\n", __FILE__,__LINE__, ds_name);
    return;
  }

  /* create record in drms */
  rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);
  if (status)
  {
    printkerr("Error at %s, line %d: Cannot create record using this data"
              " series name:<%s>.Existing program.\n",__FILE__,__LINE__, query);
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
    printkerr("Error at %s, line %d: Could not set PACKET_TIME keyword.",
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

  /* set PACKET_VERSION_NUMBER keyword using data from HMI_VER.. keyword */
  /* set drms type and long telemetry name  */
  keytype= DRMS_TYPE_STRING;
  strcpy(keyname, "PACKET_VERSION_NUMBER");
  /*allocate memory for string value */
  key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
  /* set packet version number */
  strcpy(key_anyval.string_val, pkt_ver_num);
  /* set record */
  status = drms_setkey(rec, keyname, keytype, &key_anyval);
  /* free memory */
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
    printkerr("Error at %s, line %d: Cannot close drms record.\n", 
               __FILE__,__LINE__);
    return;
  }
  else
  {
    ;//printf("Successfully closed record\n");
  }
}



/**********************************************************************************
 * TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);
 * Note on time codes.
 * SDO/HMI,AIA keeps time in a 48-bit counter in units of 1/(2^16) seconds.  Thus
 * the top 32 bits is a seconds counter and the bottom 16 bits is a sub-seconds
 * counter.  The epoch is 1958.01.01_00:00:00.
 * Thus to convert HMI,AIA instrument time in two variables, e.g. SHS and SHSS to
 * a DRMS time the conversion is:  t_drms = SDO_EPOCH + SHS + SHSS/65536.0
 * where SDO_EPOCH = sscan_time("1958.01.01_00:00:00");
 * TAI and UTC are same at 1 Jan 1958.
 **********************************************************************************/
TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss)
{
/*changes done:
1.changed args from float to int 
2.added line below...int ss.. 
3.changed return statement with ss parameter!
*/
int ss=(sdo_ss >> 16) & 0xFFFF;
static int firstcall = 1;
static TIME sdo_epoch;
if (firstcall)
  { /* time_1958 - time_1977_TAI, to be added to SDO time to get DRMS time */
  firstcall = 0;
  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  }
return(sdo_epoch + (TIME)sdo_s + (TIME)ss/65536.0);
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
    printkerr("Error at %s, line %d: Did not find TIMECODE_SECONDS value for"
              "calculating the PACKET_TIME keyword and index. Returning error"
              "status.\n",__FILE__,__LINE__);
    return 0;
  }
  else
  {
    if (!SUBSEC_FOUND_FLAG)
    {
      printkerr("Error at %s, line %d: Did not find TIMECODE_SUBSECS value for"
                "calculating the PACKET_TIME keyword and index. Returning error"
                "status, but setting time using seconds only.\n",__FILE__,__LINE__);
      *ptime =SDO_to_DRMS_time(sec, 0);
      return 0;
    }
    else
    { 
      *ptime =SDO_to_DRMS_time(sec, subsec);
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
char* lookup_dsn( char *pvn, int cnt)
{
  /* variables */
  int i;
  int found_jvn=0;

  /* lookup jsvn associated with pvn */
  for(i=0; i < cnt ;i++)
  {
    /* check if found pvn in packet */
    if (!strcmp (jmap->pvn[i], pvn ))
    {
      found_jvn = 1;
      break;
    }
  }
  if ( found_jvn != 1)
  {
    printkerr("Error at %s, line %d: Did not find the data series name "
              "for this packet version number. Maybe this is new configuration"
              "and a new JVN-TO-PVN map files need to updated using script.\n",
              __FILE__,__LINE__);
  }
  return ((char *)jmap->dsn[i]);
}



/**************************************************************************
 * PRINT PACKET VALUES 1                                                  *
 * FUNCTION: saveprint_packet_values1(unsigned char [])                   *
 * DESCRIPTION: Prints translated and decoded packet values and writes    *
 * values to DRMS using the data series name passed in argument list      *
 **************************************************************************/
void saveprint_packet_values1(unsigned char read_in_buffer[], char *in, char *out)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  int i,j,k,y,s;
  int apid;
  int factor;
  int packet_length;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;

  for(k=0,factor=0,packet_length=0;  ; k++)
  {
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;
    /* Check if at end of all pkts */
    if (read_in_buffer[5 + factor] == (int)NULL)
    {
       printf("\n\n*****At end of file. Found %d packets in file.\n", k);
       break;
    }
    /* get packet lenght */
    packet_length=  read_in_buffer[5 + factor];

    /* apid */
    /* set 0th to 7th bits */
    apid =  (unsigned  short int)( (read_in_buffer[ 1 + factor]) & 0x00FF );
    /* set 8th to 15th bits */
    apid |= (unsigned  short int)( (read_in_buffer[ 0 + factor] << 8) & 0xFF00 );
    apid &= 0x07FF;

    /* version number */
    /*check if packet version is 0.0 -skip-print warning*/
    if (( read_in_buffer[14 + factor ] == 0) && (read_in_buffer[15 + factor ] == 0))
    {
      printkerr("Warning at %s, line %d: Getting 0.0 for packet version number."
                "in packet data. Skip processing this packet, don't load to "
                "DRMS and get next packet", __FILE__,__LINE__);
      continue;
    }

    /* version number */
    /*check if packet version is 0.0 -skip-print warning*/
    if (( read_in_buffer[14 + factor ] == 0) && (read_in_buffer[15 + factor ] == 0))
    {
      printkerr("Warning at %s, line %d: Getting 0.0 for packet version number."
                "in packet data. Skip processing this packet, don't load to "
                "DRMS and get next packet", __FILE__,__LINE__);
      continue;
    }
    sprintf(packet_version_number,"%03d.%03d",read_in_buffer[14 + factor ], read_in_buffer[15 + factor]);
    /* check for version number set to zero and print warning messages. */
    if ( read_in_buffer[14 + factor ] == 0 && read_in_buffer[15 + factor ])
    {
	printf("decode_dayfile: bad version number:going to set version number to 0.00 by default for now\n");
        continue;
    }
    /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
    /* set buffer to values in packet and set packet version number to value in packet or some default value*/
    for (i = 0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
      if ( (i == (14 + factor ))  &&  (read_in_buffer[14 + factor ] == 0)) 
      { 
       	/*can set default value here if need to */
	//hk_pkt_buffer[j++]= 0x01;
	hk_pkt_buffer[j++]= 0x00; /*can set default value here if need to */
      }
      else if  ( (i == (15 + factor)) && (read_in_buffer[15 + factor] == 0))
      {
	/*can set default value here if need to */
	//hk_pkt_buffer[j++]= 0x49; /*0x49 = 73(dec) */
	hk_pkt_buffer[j++]= 0x00;
      }
      else if  ( i == (15 + factor) ||  i == (14 + factor))  
      {
	/* using packet version value found in packet */
        hk_pkt_buffer[j++]= read_in_buffer[i];
      }
      else  
      { 
	/* set other values in array */
	hk_pkt_buffer[j++]= read_in_buffer[i];
      }
    } /* end for loop setting values for one extracted packet in buffer */

    /* print packet keywords */
    printf ( "--------------------------------------------------------------------------\n");
    printf("Packet Number in day file = %-3d       Packet Lenght          = %d \n", k, packet_length);
    printf("Apid                      = %-3.3x       Packet Version Number  = %s\n", apid,packet_version_number);
    printf("In Day File               = %s\n", in);
    printf("Out Series Name           = %s\n", out);

    /* put in format for decode_hk to read by adjusting values in buffer *
     * from unsigned char to unsigned short                              */
    printf ( "--------------------------------------------------------------------------\n");
    printf("Packet Values in unsigned shorts: \n"  );
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      s_hk_pkt_buffer[y] = (unsigned short)(hk_pkt_buffer[i + 1] << 8  & 0xFF00 );
      s_hk_pkt_buffer[y] = (unsigned short)((hk_pkt_buffer[i] & 0x00FF) + s_hk_pkt_buffer[y]) ;
      /*print information contained in packet being sent to decode_hk function */
      printf("%2.2x",  s_hk_pkt_buffer[y]  & 0x00FF); 
      printf("%2.2x ",  s_hk_pkt_buffer[y] >> 8 & 0x00FF); 
      if ( i ==  14 || i == 30  || i == 46  || 
           i == 62  || i == 78  || i == 94  || 
           i == 110 || i == 126 || i == 142 ||
           i == 158 || i == 174 || i == 190 ||
           i == 206 || i == 222 || i == 238 )
        printf("\n");
   }
   printf("\n");

    /* send hk_pkt_buffer to decoder function */
    word_ptr = s_hk_pkt_buffer;
    s = decode_hk_keywords(word_ptr, apid,  &ccsds.keywords); 
    if (s) 
    {
      printkerr("Error at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n", __FILE__,__LINE__, s);
    }
    kw_head= ccsds.keywords;

    /* start printing report on keywords */
    printf ( "--------------------------------------------------------------------------\n");
    printf ("LONG TELEM MNEMONIC NAME       SHORT    RAW VALUE  ENGR VALUE  ENGR VALUE\n");
    printf ("                               KEYWORD    (HEX)     (DECIMAL)    (HEX)   \n");

    /* loop throught all keywords in single packet*/
    for ( i=0; ccsds.keywords ; i++)
    {
      printf ( "---------------------------------------------------------------------------\n");
      printf ( "%-31.29s", ccsds.keywords->name);
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

    /* write values to DRMS using the data series name pass in argument list*/
    write_to_drms( packet_version_number, out, kw_head );
    print_packet_time(kw_head);

  }/* loop throught next packet and do all packets in file*/
}



/**************************************************************************
 * Load Data Series Names                                                 *
 * FUNCTION:  load_ds_name                                                *
 * DESCRIPTION: Loads data series names in array with packet version in   *
 *              corresponding array.                                      *
 * STATUS: Completed.                                                     *
 **************************************************************************/
int load_ds_names(char* pn, char* didn, int apid)
{
  /* variables */
  FILE *file_ptr;
  char directory_filename[HKDDF_MAX_FILE_NAME];
  char fn[HKDDF_MAX_FILE_NAME];
  char sn[HKDDF_MAX_FILE_NAME];
  char *directory ;
  char line[HKDDF_MAXLINE_IN_FILE];
  int wn, dn;
  int count;

  /* initialize variables */
  char *filename=fn;
  char *suffix_filename=sn;
  int k=0;

  /* get directory and file name  with jsoc version numbers*/
  directory = (char *)getenv("HK_JSVNMAP_DIRECTORY");
  suffix_filename= (char *)getenv("HK_SUFFIX_JSVNMAP_FILENAME");
  if(suffix_filename == NULL)
  {
    printf("Error: Set environxment variable <HK_SUFFIX_JSVNMAP_FILENAME>\n");
    printkerr("Error at %s, line %d: Set environment variable <HK_SUFFIX_JSVNMAP_FILENAME>\n",
                __FILE__,__LINE__);
  }
  if(directory == NULL)
  {
    printf("Error: Set environment variable <HK_JSVNMAP_DIRECTORY>\n");
    printkerr("Error at %s, line %d: Set environment variable <HK_JSVNMAP_DIRECTORY>\n",
                __FILE__,__LINE__);
  }
  /* make filename */
  sprintf(filename,"%4.4d%s",apid,suffix_filename);
  sprintf(directory_filename, "%s/%s", directory, filename);
  /* get apid */
  sprintf(jmap->apid,"%d",apid);
  /* open file and read */
  file_ptr=fopen(directory_filename,"r");
  if (!file_ptr)
  {
    printkerr("Error at %s, line %d: Check if directory and filename are",
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
        "%s | %d.%d | %*s | %*s", jmap->jvn[k], &wn,&dn);
      sprintf(jmap->pvn[k],"%03d.%03d", wn,dn);
      /* if jsvn is 0 make dsn null else set data series name in array*/
      if (!strcmp("0000",  jmap->jvn[k]))
      {
        strcpy(jmap->dsn[k], "\0");
      }
      else
      {
        sprintf(jmap->dsn[k], "%s.%s_%04d_%s", pn, didn, apid, jmap->jvn[k]);
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
void  save_packet_values2(unsigned char read_in_buffer[], char *pn,char *didn)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  int i,j,k,y,x,s;
  int apid;
  int count;
  int factor;
  int packet_length;
  char *dname;
  char data_seriesname[HKDDF_MAX_DSNAME_STR];
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;

  for(k=0,factor=0, packet_length=0;  ; k++)
  {
    /* init dsn */
    for(x=0; x < HKDDF_MAX_DSNAME_STR; data_seriesname[x]='\0', x++);
    
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;
    /* Check if at end of all pkts */
    if (read_in_buffer[5 + factor] == (int)NULL)
    {
       break;
    }
    /* get packet lenght */
    packet_length=  read_in_buffer[5 + factor];
    /* apid */
    /* set 0th to 7th bits */
    apid =  (unsigned  short int)( (read_in_buffer[ 1 + factor]) & 0x00FF );
    /* set 8th to 15th bits */
    apid |= (unsigned  short int)( (read_in_buffer[ 0 + factor] << 8) & 0xFF00 );
    apid &= 0x07FF;
    /* version number */
    /*check if packet version is 0.0 -skip-print warning*/
    if (( read_in_buffer[14 + factor ] == 0) && (read_in_buffer[15 + factor ] == 0))
    {
      printkerr("Warning at %s, line %d: Getting 0.0 for packet version number."
                "in packet data. Skip processing this packet, don't load to "
                "DRMS and get next packet", __FILE__,__LINE__);
      continue;
    }
    sprintf(packet_version_number,"%03d.%03d",read_in_buffer[14 + factor ], read_in_buffer[15 + factor]);
    /* check for version number set to zero and print warning messages. */
    if ( read_in_buffer[14 + factor ] == 0 && read_in_buffer[15 + factor ])
    {
	printf("decode_dayfile: bad version number:going to set version number to 0.00 by default for now\n");
        continue;
    }
    /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
    /* set buffer to values in packet and set packet version number to value 
    /* in packet or some default value */
    for (i =0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
      if ( (i == (14 + factor ))  &&  (read_in_buffer[14 + factor ] == 0)) 
      { 
       	/*can set default value here if need to */
	//hk_pkt_buffer[j++]= 0x01;
	hk_pkt_buffer[j++]= 0x00; /*can set default value here if need to */
      }
      else if  ( (i == (15 + factor)) && (read_in_buffer[15 + factor] == 0))
      {
	/*can set default value here if need to */
	//hk_pkt_buffer[j++]= 0x49; /*0x49 = 73(dec) */
	hk_pkt_buffer[j++]= 0x00;
      }
      else if  ( i == (15 + factor) ||  i == (14 + factor))  
      {
	/* using packet version value found in packet */
        hk_pkt_buffer[j++]= read_in_buffer[i];
      }
      else  
      { 
	/* set other values in array */
	hk_pkt_buffer[j++]= read_in_buffer[i];
      }
    } /* end for loop setting values for one extracted packet in buffer */
    /* get data name by automatically creating JSVN value based on */
    /* packet version number                                       */
    count = load_ds_names(pn, didn, apid );
    dname = lookup_dsn( packet_version_number, count);
    strcpy(data_seriesname, dname);
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
      printkerr("Error at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n", __FILE__,__LINE__, s);
    }
    kw_head= ccsds.keywords;
    /* write values to DRMS with auto generated data series name*/
    write_to_drms(packet_version_number, data_seriesname, kw_head );
  }/* loop throught next packet and do all packets in file*/
}


/**************************************************************************
 * PRINT PACKET VALUES 2                                                  *
 * FUNCTION: saveprint_packet_values2()                                       *
 * DESCRIPTION: Prints translated and decoded packet values and writes    *
 * values to DRMS. Also automatically creates data series name  based on  *
 * the packet version number in packet                                    *
 * STATUS: completed                                                      *
 **************************************************************************/
void  saveprint_packet_values2(unsigned char read_in_buffer[], char *in, 
                           char *pn,char *didn)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  int i,j,k,y,x,s;
  int apid;
  int count;
  int factor;
  int packet_length;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;
  char data_seriesname[HKDDF_MAX_DSNAME_STR];
  char *dname;

  for(k=0,factor=0, packet_length=0;  ; k++)
  {
    /* init dsn */
    for(x=0; x < HKDDF_MAX_DSNAME_STR; data_seriesname[x]='\0', x++);
    
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;
    /* Check if at end of all pkts */
    if (read_in_buffer[5 + factor] == (int)NULL)
    {
       printf("\n\n*****At end of file. Found %d packets in file.\n", k);
       break;
    }
    /* get packet lenght */
    packet_length=  read_in_buffer[5 + factor];
    /* apid */
    /* set 0th to 7th bits */
    apid =  (unsigned  short int)( (read_in_buffer[ 1 + factor]) & 0x00FF );
    /* set 8th to 15th bits */
    apid |= (unsigned  short int)( (read_in_buffer[ 0 + factor] << 8) & 0xFF00 );
    apid &= 0x07FF;
    /* version number */
    /*check if packet version is 0.0 -skip-print warning*/
    if (( read_in_buffer[14 + factor ] == 0) && (read_in_buffer[15 + factor ] == 0))
    {
      printkerr("Warning at %s, line %d: Getting 0.0 for packet version number."
                "in packet data. Skip processing this packet, don't load to "
                "DRMS and get next packet", __FILE__,__LINE__);
      continue;
    }
    sprintf(packet_version_number,"%03d.%03d",read_in_buffer[14 + factor ], read_in_buffer[15 + factor]);
    /* check for version number set to zero and print warning messages. */
    if ( read_in_buffer[14 + factor ] == 0 && read_in_buffer[15 + factor ])
    {
	printf("decode_dayfile: bad version number:going to set version number to 0.00 by default for now\n");
        continue;
    }
    /* Extract hk packets - initialize array with zeros */
    for(i=0; i < HKDDF_MAX_PKT_SIZE;i++) hk_pkt_buffer[i]=0x00;
    /* set buffer to values in packet and set packet version number to value 
    /* in packet or some default value */
    for (i =0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
      if ( (i == (14 + factor ))  &&  (read_in_buffer[14 + factor ] == 0)) 
      { 
       	/*can set default value here if need to */
	//hk_pkt_buffer[j++]= 0x01;
	hk_pkt_buffer[j++]= 0x00; /*can set default value here if need to */
      }
      else if  ( (i == (15 + factor)) && (read_in_buffer[15 + factor] == 0))
      {
	/*can set default value here if need to */
	//hk_pkt_buffer[j++]= 0x49; /*0x49 = 73(dec) */
	hk_pkt_buffer[j++]= 0x00;
      }
      else if  ( i == (15 + factor) ||  i == (14 + factor))  
      {
	/* using packet version value found in packet */
        hk_pkt_buffer[j++]= read_in_buffer[i];
      }
      else  
      { 
	/* set other values in array */
	hk_pkt_buffer[j++]= read_in_buffer[i];
      }
    } /* end for loop setting values for one extracted packet in buffer */
    /* get data name by automatically creating JSVN value based on */
    /* packet version number                                       */
    count = load_ds_names(pn, didn, apid);
    dname = lookup_dsn( packet_version_number, count);
    strcpy(data_seriesname, dname);
    /* print packet keywords */
    printf ( "--------------------------------------------------------------------------\n");
    printf("Packet Number in day file = %-3d       Packet Lenght          = %d \n", k, packet_length);
    printf("Apid                      = %-3.3x       Packet Version Number  = %s\n", apid,packet_version_number);
    printf("In Day File               = %s\n", in);
    printf("Out Series Name           = %s\n", data_seriesname);

    /* put in format for decode_hk to read by adjusting values in buffer  */
    /* from unsigned char to unsigned short and write bytes to report     */
    printf ( "--------------------------------------------------------------------------\n");
    printf("Packet Values in unsigned shorts: \n"  );
    for (i=0, y=0 ; i < packet_length + 6 + 1   ; i += 2, y++)
    {
      s_hk_pkt_buffer[y] = (unsigned short)(hk_pkt_buffer[i + 1] << 8  & 0xFF00 );
      s_hk_pkt_buffer[y] = (unsigned short)((hk_pkt_buffer[i] & 0x00FF) + s_hk_pkt_buffer[y]) ;
      /*print information contained in packet being sent to decode_hk function */
      printf("%2.2x",  s_hk_pkt_buffer[y]  & 0x00FF); 
      printf("%2.2x ",  s_hk_pkt_buffer[y] >> 8 & 0x00FF); 
      if ( i ==  14 || i == 30  || i == 46  || 
           i == 62  || i == 78  || i == 94  || 
           i == 110 || i == 126 || i == 142 ||
           i == 158 || i == 174 || i == 190 ||
           i == 206 || i == 222 || i == 238 )
        printf("\n");
    }
    printf("\n");

    /* send hk_pkt_buffer to decoder function */
    word_ptr = s_hk_pkt_buffer;
    s = decode_hk_keywords(word_ptr, apid,  &ccsds.keywords);
    if (s) 
    {
      printkerr("Error at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n", __FILE__,__LINE__, s);
    }
    kw_head= ccsds.keywords;

    /* start print report on keywords */
    printf ( "--------------------------------------------------------------------------\n");
    printf ("LONG TELEM MNEMONIC NAME       SHORT    RAW VALUE  ENGR VALUE  ENGR VALUE\n");
    printf ("                               KEYWORD    (HEX)     (DECIMAL)    (HEX)   \n");

    /* loop throught all keywords in single packet*/
    for ( i=0; ccsds.keywords ; i++)
    {
      printf ( "---------------------------------------------------------------------------\n");
      printf ( "%-31.29s", ccsds.keywords->name);
      printf ( "%-9.8s",  ccsds.keywords->fitsname);
      //printf ( "%-8.8x",  ccsds.keywords->raw_value);
      printf ( "%-8.8lx", (long int) ccsds.keywords->raw_value);

      //printf ( "%d",ccsds.keywords->raw_value);
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
      else
      {
        printf("WARNING(decode_dayfile):****Found*** NO KW_TYPE\n");
        printf("Eng type is %d\n", (int)ccsds.keywords->eng_type);
      }

       /*get next packet */
       ccsds.keywords= ccsds.keywords->next;

    }/*end of loop throught all keywords in single packet*/ 


    /* write values to DRMS with auto generated data series name*/
    write_to_drms( packet_version_number, data_seriesname, kw_head );
    /*now print packet time */
    print_packet_time(kw_head);
    
  }/* loop throught next packet and do all packets in file*/
}


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
    printf ( "--------------------------------------------------------------------------\n");
    printf("ERROR: COULD NOT SET PACKET_TIME KEYWORD\n");
    printf ( "--------------------------------------------------------------------------\n");
    printf ("%s\n",HKDDF_CNTRL_L);
    return;
  }
  else 
  { 
    printf ( "--------------------------------------------------------------------------\n");
    printf("PACKET_TIME                                     %15.5lf\n", pkt_time);
    printf ( "--------------------------------------------------------------------------\n");
    printf ("%s\n",HKDDF_CNTRL_L);
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
void save_packet_values1(unsigned char read_in_buffer[], char *out)
{
  /* variables */
  CCSDS_Packet_t ccsds;
  HK_Keyword_t *kw_head;
  int i,j,k,y,s;
  int apid;
  int factor;
  int packet_length;
  char packet_version_number[HKDDF_MAX_PVN_SIZE];
  unsigned char hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short s_hk_pkt_buffer[HKDDF_MAX_PKT_SIZE];
  unsigned short *word_ptr;

  for(k=0,factor=0,packet_length=0;  ; k++)
  {
    /* set pointer to beginning of packets in buffer using factor parameter */
    factor = k * (packet_length + 6 + 1 ) ;
    /* Check if at end of all pkts */
    if (read_in_buffer[5 + factor] == (int)NULL)
    {
      break;
    }
    /* get packet lenght */
    packet_length=  read_in_buffer[5 + factor];
    /* apid */
    /* set 0th to 7th bits */
    apid =  (unsigned  short int)( (read_in_buffer[ 1 + factor]) & 0x00FF );
    /* set 8th to 15th bits */
    apid |= (unsigned  short int)( (read_in_buffer[ 0 + factor] << 8) & 0xFF00 );
    apid &= 0x07FF;
    /* get packet version number */
    /*check if packet version is 0.0 */
    if (( (int)read_in_buffer[14 + factor ] == 0 ) && ((int)read_in_buffer[15 + factor ] == 0))
    {
      printkerr("Warning at %s, line %d: Getting 0.0 for packet version number."
                "in packet data. Skip processing this packet, don't load to "
                "DRMS and get next packet", __FILE__,__LINE__);
      continue;
    }
    sprintf(packet_version_number, "%03d.%03d",read_in_buffer[14 + factor ], read_in_buffer[15 + factor]);
    /* set buffer to values in packet and set packet version number to value in packet */
    for (i = 0 + factor, j=0 ; i <  packet_length + 6 + 1 + factor; i++)
    {
      if  ( i == (15 + factor) ||  i == (14 + factor))  
      {
        /* using packet version value found in packet */
        hk_pkt_buffer[j++]= read_in_buffer[i];
      }
      else  
      { 
        /* set other values in array */
	hk_pkt_buffer[j++]= read_in_buffer[i];
      }
    }/* end for loop setting values for one extracted packet in buffer */
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
      printkerr("Error at %s, line %d: decode_hk_keyword function returned"
                " error status number <%d>\n", __FILE__,__LINE__, s);
    }
    kw_head= ccsds.keywords;
    /* write values to DRMS*/
    write_to_drms( packet_version_number, out, kw_head );
  } /* loop throught next packet and do all packets in file*/
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
      printf("Error:Can't open environment variable file <%s>. Check setting is correct\n", envfile);
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
