/*############################################################################
# Name:        load_mp.c - Loads image location, leg, iss status data in     #
#                          master pointing drms series.                      #
# Description: Load image location, leg, and iss status  keyword values      #
#              in master pointing drms series. Gets image location keyword   #
#              values by reading the image location file defined in argument #
#              in input argument to executable using ilf parameter. The hmi  #
#              leg keyword data is read from hmi.leg_status drms series and  #
#              loaded in master pointing drms series. The hmi iss and aia iss#
#              keyword data is read from the hmi.iss_status,                 #
#              aia.iss_1_2_status, and aia.iss_3_4_status drms series.       #
# Execution:   load_mp  in=<day-file> isf=<instruction-file>                 #
#              (1)To load master pointing drms series keyword values:        #
#                    load_mp ilf=<image location file>                       #
# Example  :  load_mp ilf=/tmp20/image_location.txt                          #
# Limitation:  Setup required for environment variables in file              #
#              SOURCE_ENV_FOR_LOAD_MP. This file contains the name of the    #
#              master pointing series, names of leg,  and names of iss series#
#              to use.                                                       #
# Author:      Carl                                                          #
# Date:        January,29, 2010                                              #
############################################################################*/

/**
   @defgroup load_mp load_mp
   @ingroup su_lev1

   @brief Loads master pointing data into keyword names for master pointing DRMS series.

   @par Synopsis:
   @code
   load_mp -h
   load_mp  ilf=<full path to image location file text file> 
   @endcode

   Loads master pointing data into keyword names for master pointing DRMS series.
   The data is loaded from a image location text file and 4 DRMS series. The image location text file
   contains both HMI and AIA image location data. The HMI image location keyword data  will be x, y,
   imscale and instrot data for carmera 1 and 2. The AIA image location keyword data will be x, y, 
   imscale and instrot data for the different wavelenghts. The 4 DRMS series used to get keyword data
   and load into the master pointing DRMS series includes the following series:hmi.leg_status, hmi.iss_status,
   aia.iss_1_2_status and aia.iss_3_4_status. The record to get data for each series is the 
   HK_VALS keyword value which is found in the image location text file. The prime index keyword
   for the master pointer series is T_START keyword. This T_START values is set with the value
   in the image location text file(i.e.,image_location.txt file). The DATE keyword value is updated when updating
   a record's values or when create new record. The DATE is set to current time when update occurred. 

   The VERSION keyword is used to determine if the record in master pointing DRMS series has been updated 
   except for the latest record. For example if a record has been updated 5 times the version value will be 5. 
   The latest record's version is always set to 0 even when updated.  When a newer latest record is added,
   then the previous latest record is set to version 1 and the newer latest record's version is 0. 

   The T_STOP keyword for the first record in the master pointing DRMS series is intially set using the
   value in image location text file. As new records get added the T_STOP values is determined by where
   this new record is added in the master pointing series. The master pointing DRMS series is sorted by
   index T_START time. So when insert update to record or add a new record, the T_STOP values is determined
   by the next records T_START time or if this record is the latest record then use the T_STOP time in the
   image location text file. 

   The debug messages for C code can be turned on by compiling with -DDEBUG_LOAD_MP. Here is example of adding
   flag to ICC_CF_ICCCOMP line in make_basic.mk. This will turn on debug messages.
   ICC_CF_ICCCOMP  = -DDEBUG_LOAD_MP

   The setup of the environment variables in file SOURCE_ENV_FOR_LOAD_MP tells executable the name of the 
   master pointing series to write keyword values. The environment variables tell executable the 
   names of series to get data from to use to set aia iss data and hmi iss and leg data in the master pointing series. 
   The are a few variable used to point to housekeeping configuration information. The setup of the SOURCE file to 
   use is in executable code using define variable ENVFILE. Here are example settings of the variables in file.

   setenv  HK_CONFIG_DIRECTORY                  /home/production/cvs/TBL_JSOC/lev0/hk_config_file/
   setenv  HK_GTCIDS_FILE                       gtcids.txt
   setenv  HK_SHCIDS_DIRECTORY                  /home/production/cvs/TBL_JSOC/lev0/sdo_hk_config_file/
   setenv  HK_SHCIDS_FILE                       shcids.txt
   setenv  HK_LMP_MP_SERIESNAME                 sdo.master_pointing
   setenv  HK_LMP_HMI_ISS_SERIESNAME            hmi.iss_status
   setenv  HK_LMP_HMI_LEG_SERIESNAME            hmi.leg_status
   setenv  HK_LMP_AIA_ISS12_SERIESNAME          aia.iss_1_2_status
   setenv  HK_LMP_AIA_ISS34_SERIESNAME          aia.iss_3_4_status



   The ilf parameter is a mandatory argument which should contain the directory and filename of the image 
   location text file. The image location text file T_START, T_STOP, T_HKVALS UTC times(i.e.,
   2009.08.11_17:00:00.00_UTC). The image location text files contains the list of keywords to retrieve 
   from drms data series to the set in the master pointing drms series. The image location text files contain
   the VERSION setting.


   @par Flags:
   @c -h: Shows usage message.
   @par

   @param ilf The full directory path and file name to image location text file(required field).
 
   @par Example of running script:
   @code
   load_mp  ilf=/home/production/cvs/TBL_JSOC/lev1/image_loc_file/image_location.txt
   @endcode

   @par Example of running script debug conditional compile flag turned on:
   @code
   load_mp   ilf=/home/production/cvs/TBL_JSOC/lev1/image_loc_file/image_location.txt > Debug-Log
   @endcode

   @par Example of running help:
   @code
   load_mp -h
   @endcode

*/
/* Defined constants */
/******************** defines ***********************************************/
#define HKLMP_FAILED_STATUS            1
#define HKLMP_PASSED_STATUS            0
#define HKLMP_MAX_DSNAME_STR           100
#define HKLMP_MAX_DIR_FILE_NAME        200
#define HKLMP_MAX_KEYWORD_NAME_STR     100
#define HKLMP_MAX_QUERY_STR            200
#define HKLMP_MAX_CHARS_IN_LINE        200
#define HKLMP_PACKET_TIME_STR          100
#define HKLMP_VAR_VERSION              8
#define HKLMP_VAR_T_HKVALS             4
#define HKLMP_VAR_T_STOP               2
#define HKLMP_VAR_T_START              1
#define HKLMP_NO_IMAGE_LOC_FILE        100
/*#define ENVFILE  "/home/production/cvs/JSOC/proj/lev1/apps/SOURCE_ENV_FOR_LOAD_M3SD" */
/*#define ENVFILE    "/home3/carl/cvs/JSOC/proj/lev1/apps/SOURCE_ENV_FOR_LOAD_M3SD"*/
#define ENVFILE    "/home/production/cvs/JSOC/proj/lev1/apps/SOURCE_ENV_FOR_LOAD_MP"

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
  {ARG_STRING, "ilf", "Not Specified", "full path to image location file"},
  {ARG_END}
};
ModuleArgs_t   *ggModArgs=module_args;
char* module_name = "load_mp";


/******************* function prototypes  *******************************/
TIME set_tstop_value(char *seriesname, char tstart_str[HKLMP_PACKET_TIME_STR],TIME *ptstart);
int set_mp_version(char *seriesname, char tstart_str[HKLMP_PACKET_TIME_STR]);
int get_hmilegstatus_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR]);
int get_hmiissstatus_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR]);
int get_aiaiss12status_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR]);
int get_aiaiss34status_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR]);
static void  set_env_variables();
void update_tstop_version(char query[HKLMP_MAX_QUERY_STR], char mp_tstart_value[HKLMP_PACKET_TIME_STR]);
void check_oldrec_tstop_time(char query[HKLMP_MAX_QUERY_STR], char mp_tstart_value[HKLMP_PACKET_TIME_STR]);
void my_usage (void);
int set_image_loc_values(DRMS_Record_t *record, char *seriesname, TIME *ptr_tstart, char mp_tstart_str[HKLMP_PACKET_TIME_STR] , TIME *ptr_hkvals, char hk_time_val_str[HKLMP_PACKET_TIME_STR], char im_loc_fn[HKLMP_MAX_DIR_FILE_NAME]);

/********************* extern functions  *********************************/
extern int DoIt(void);
extern int nice_intro (void);
extern void sprint_time (char *at, TIME t, char *zone, int precision);


/* @} */



/*************************************************************************
 * my_usage                                                              *
 * Function: myusage(void)                                               *
 * Description: Use function to display usage                            *
 *************************************************************************/
void my_usage (void)
{
  
  printf ("Usage:\nload_mp  [-h] "
    "ilf=<image location file>  \n"
    "  details are:\n"
    "  -h: help - show this message then exit(optional field)\n"
    "  ilf=<full directory path to image location file name> -use full path to image location file file(required field)\n"
    "  note:need master pointing series already created for this program to load keywords in data series.\n"
    "  note:need to create image_location using set_mp.pl script or by hand editing file using correct file format.\n"
    "  Example of running:\n"
    "  load_mp  ilf=/home/carl/cvs/TBL_JSOC/lev1/image_location_file/image_location_2009_10_29.txt\n");
}



/*************************************************************************
 * Nice Intro                                                            *
 * Function: nice_intro(void)                                            *
 * Description: Use function to display help                             *
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
 * Function: DoIt(void                                                   *
 * Description: Use function to create jsoc module                       *
 *************************************************************************/
int DoIt(void)
{
  /* variables */
 /* drms record create variables */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  FILE *file_ptr;
  TIME *p_hkvals;
  TIME *p_tstart;
  TIME hk_vals;
  TIME tstart_vals;
  char query[HKLMP_MAX_QUERY_STR];
  char hk_directory_filename[HKLMP_MAX_DIR_FILE_NAME];
  char hk_time_value[HKLMP_PACKET_TIME_STR];
  char mp_tstart_value[HKLMP_PACKET_TIME_STR];
  char mp_seriesname[HKLMP_MAX_DSNAME_STR];
  char *hk_df_fn;
  char *sn;
  int set_image_status;
  unsigned int t_status=0;
  int status;

  /* initialize variables */
  p_hkvals=&hk_vals;
  p_tstart=&tstart_vals;
  sn=mp_seriesname;

  /* set environment variables */
  set_env_variables();

  /* parameter initialization */
  hk_df_fn= hk_directory_filename;

  /* Get command line arguments */
  char *ilf = cmdparams_get_str (&cmdparams, "ilf", NULL);

  /* check arguments used */
  if (nice_intro ()) return (0);

  /* check if entered day file name */
  if (ilf == NULL) 
  {
    printkerr("ERROR at %s, line %d: Need to enter image location file name! "
              "Exiting program.\n ", __FILE__,__LINE__);
    my_usage();
    return (0);
  }

  /* get  ilf(Image Location Filename) filename and open hk dayfile*/
  strcpy(hk_df_fn, ilf);
  strcat(hk_df_fn, "\0");
  file_ptr=fopen(hk_df_fn,"r");

  /* check if image location file exists */
  if (!file_ptr)
  {
    printkerr("ERROR at %s, line %d: Please check filename and directory is correct. "
              " Could not find or open -ilf- directory and filename: "
              "<%s>. Example format for -ilf- file: ilf=/home/rock/image_location.txt. Exiting execution.\n", 
              __FILE__,__LINE__, hk_df_fn);
    return(0);
  }
  /* close file after quick check if there */
  fclose(file_ptr);

  /* get master pointing series name */
  sn = (char *)getenv("HK_LMP_MP_SERIESNAME"); 
  if(sn == NULL)
  {
    printkerr("ERROR at %s, line %d: no master pointing seriesname environment variable set. "
              " Set envirionment HK_LMP_MP_SERIESNAME variable in sourced file "
              " SOURCE_ENV_FOR_LOAD_MP.\n", __FILE__,__LINE__);
    return (HKLMP_FAILED_STATUS);
  }
  strcpy(query,sn);
  strcat(query,"\0");
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(1): Master pointing series name is <%s>\n",query);
#endif
 
  /* create record in drms for master_pointing table */
  rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);
  if (status)
  {
    printkerr("ERROR at %s, line %d: Cannot create record using this data"
              " series name:<%s>.Existing program.\n",__FILE__,__LINE__, query);
    return 1;
  }
#ifdef DEBUG_LOAD_MP
  else
  {
    printf("load_mp:(2):Sucessfully created record for series <%s>\n", query);
  }
#endif

  /* set master pointing series rec to pass to following function to set master pointing series image location data */
  rec = rs->records[0];

#ifdef DEBUG_LOAD_MP
    printf("load_mp:(3):Calling set_image_loc_values to set image location values in  <%s>\n", query);
#endif
  /* parse data in image_location.txt file and save keyword name and values to master pointing drms series */
  /* set the T_START and T_STOP times based T_HKVALS time in master pointing table */
  set_image_status=set_image_loc_values(rec, sn, p_tstart, mp_tstart_value, p_hkvals, hk_time_value, ilf);
  if (set_image_status > 0 && set_image_status < 15)
  {
     printkerr("ERROR at %s, line %d: Returned error status <%d> from set_image_loc_values function. "
              "Exiting from load_mp executable. Please add in the required keywords listed.\n",
               __FILE__,__LINE__, set_image_status);
     /* check bit 4 is set, if not display message */
     if (!((unsigned int)set_image_status>>3 & 0x00000001))
         printkerr("-->>Add VERSION keyword to image location file.\n");
     /* check bit 1 is set, if not display message */
     if (!((unsigned int)set_image_status & 0x00000001))
         printkerr("-->>Add T_START keyword to image location file.\n");
     /* check bit 3 is set, if not display message */
     if (!((unsigned int)set_image_status>>2 & 0x00000001))
         printkerr("-->>Add T_HKVALS keyword to image location file.\n");
     /* check bit 2 is set, if not display message */
     if (!((unsigned int)set_image_status>>1 & 0x00000001))
         printkerr("-->>Add T_STOP keyword to image location file.\n");

     /* exit load_mp because of one or more errors above */
     return(0);
  }
  else if ( set_image_status == HKLMP_NO_IMAGE_LOC_FILE )
  {
     printkerr("ERROR at %s, line %d: Returned error status <%d> from set_image_loc_values function. "
              "Exiting from load_mp executable. Please check image location file name's directory and "
              "filename are correct.\n", __FILE__,__LINE__, set_image_status);
     /* exit load_mp because of one or more errors above */
     return(0);
  }
  /*else ok */

  /* get keyword data from hmi.leg_status and save  name and values to master pointing table  */
  if( get_hmilegstatus_keywords(rec, p_hkvals, hk_time_value))
  {
    printkerr("ERROR at %s, line %d: Returned error status from get_hmilegstatus_keywords function. "
              "Exiting from load_mp executable. Please use correct HKVALS setting in image loc file and rerun.\n",
               __FILE__,__LINE__);
    status = drms_close_records(rs, DRMS_FREE_RECORD);
    return(0);
  }


  /* get keyword data from hmi.iss_status and save  name and values to master pointing table  */
  if(get_hmiissstatus_keywords(rec, p_hkvals, hk_time_value))
  {
    printkerr("ERROR at %s, line %d: Returned error status from get_hmiissstatus_keywords function. "
              "Exiting from load_mp executable.  Please use correct HKVALS setting in image loc file and rerun.\n",
              __FILE__,__LINE__);
    status = drms_close_records(rs, DRMS_FREE_RECORD);
    return(0);
  }

  /* get keyword data from aia.iss_1_2_status series and save  name and values to master pointing table  */
  if(get_aiaiss12status_keywords(rec, p_hkvals, hk_time_value))
  {
    printkerr("ERROR at %s, line %d: Returned error status from get_iss12status_keywords function. "
              "Exiting from load_mp executable.  Please use correct HKVALS setting in image loc file and rerun.\n",
              __FILE__,__LINE__);
    status = drms_close_records(rs, DRMS_FREE_RECORD);
    return(0);
  }

  /* get keyword data from aia.iss_3_4_status series and save  name and values to master pointing table  */
  if(get_aiaiss34status_keywords(rec, p_hkvals, hk_time_value))
  {
    printkerr("ERROR at %s, line %d: Returned error status from get_iss34status_keywords function. "
              "Exiting from load_mp executable.  Please use correct HKVALS setting in image loc file and rerun.\n",
               __FILE__,__LINE__);
    status = drms_close_records(rs, DRMS_FREE_RECORD);
    return(0);
  }

  /* finally COMMIT values to master pointing - close record */ 
  status = drms_close_records(rs, DRMS_INSERT_RECORD);
  if (status)
  {
    printkerr("ERROR at %s, line %d: Cannot close drms record.\n",
               __FILE__,__LINE__);
    return 0;
  }
  else
  {
    ;//printf("Successfully closed record\n");
#ifdef DEBUG_LOAD_MP
      printf("load_mp:(42): Completed close of write of record of master pointing series\n" );
#endif
  }

  /* check if new T_START in image_location.txt file is greater than all other T_STARTS in master pointing table */
  (void)update_tstop_version(query, mp_tstart_value);
  (void)check_oldrec_tstop_time(query, mp_tstart_value);

  printf(". . . successfully added or updated record in master pointing series<%s> for T_START record <%s>\n",query, mp_tstart_value); 

  return 0;  
}



/**************************************************************************
 * Set Environment Variables                                              *
 * Function: set_env_variables()                                          *
 * Description: Sets Data series project name and data type name.         *
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



/**************************************************************************
 * SET TIMES IMAGE LOCATION                                               *
 * Function:set_image_loc_values()                                        *          
 * Description:Set the image location values stored in image location file*
 * to keywords in master pointing series.                                 *
 **************************************************************************/
int set_image_loc_values(DRMS_Record_t *record, char *seriesname, TIME *ptr_tstart, char mp_tstart_str[HKLMP_PACKET_TIME_STR] , TIME *ptr_hkvals, char hk_time_val_str[HKLMP_PACKET_TIME_STR], char im_loc_fn[HKLMP_MAX_DIR_FILE_NAME])
{            
  /* Pass record to set in master pointing series, seriesname for master pointing series, pass image location filename in im_loc_fn */
  /* Return back ptr_hkvals set using file data */
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  FILE *file_ptr;
  TIME now_time;
  TIME t_start_value;
  TIME t_stop_value;
  char *hk_df_fn;
  char hk_directory_filename[HKLMP_MAX_DIR_FILE_NAME];
  char keyname[HKLMP_MAX_KEYWORD_NAME_STR];
  char line[HKLMP_MAX_CHARS_IN_LINE];
  char time_str[HKLMP_PACKET_TIME_STR];
  char t_start_str[HKLMP_PACKET_TIME_STR];
  char t_stop_str[HKLMP_PACKET_TIME_STR];
  int filevar_status;
  int  ret_ver_value;
  int  status;
  struct tm timestorage;
  struct tm *time_ptr;
  time_t tvalue;

  /* initialize variables */
  time_ptr= &timestorage;
  hk_df_fn= hk_directory_filename;
  filevar_status=0;

  /*open file - get ilf(image location file) file*/
  strcpy(hk_df_fn,  im_loc_fn);
  strcat(hk_df_fn, "\0");
  file_ptr=fopen(hk_df_fn,"r");
  if (!file_ptr)
  {
    printkerr("ERROR at %s, line %d: Please check filename and directory is correct. "
              " Could not get image_location.txt  directory and filename: %s "
              " Exiting execution.\n", 
              __FILE__,__LINE__, hk_df_fn);
    return(HKLMP_NO_IMAGE_LOC_FILE);
  }

  /* get all lines and use Keyword values in file to set master pointing data series with */
  /* image location values, T_START, T_STOP and T_HKVALS values */
  while( fgets(line, MAXLINE_IN_FILE, file_ptr) != NULL )
  {
     if(line[0] == '#')
     {
       continue; /* skip comments */
     } 
     else if (!strncmp(line, "KWD", 3))
    {
       /* get keyword name */
       sscanf(line,"%*s%s%*s%*s",keyname);     
       strcat(keyname,"\0");

       /* check if one of these keywords for special processing based on logic */
       if (!strcmp(keyname, "T_START") || !strcmp(keyname, "T_STOP") || !strcmp(keyname, "T_HKVALS"))
       {
         if(!strcmp(keyname, "T_STOP"))
         {
           /* set filevar_status */
           filevar_status+=HKLMP_VAR_T_STOP;

           /* if T_START records exists then use current T_STOP value in record */
           /* else if new T_START record then find closest time to use for T_STOP */
           t_stop_value=set_tstop_value(seriesname, mp_tstart_str,ptr_tstart);
           if(t_stop_value > 0)
           {
#ifdef DEBUG_LOAD_MP
             printf("load_mp:(5):set_image_loc_values():setting T_STOP <%f>\n",t_stop_value);
#endif
             /* set drms type */
             keytype= DRMS_TYPE_TIME;
             /* set packet time */
             key_anyval.time_val= t_stop_value;
             /* set record */
             status = drms_setkey(record, keyname, keytype, &key_anyval);
             continue;
           }
           /* else skip setting here and set T_STOP to value with value in file below */

         }
         /* set T_START or T_HKVALS keyword using times in image_location.txt file */
         /* get time string */
         sscanf(line,"%*s%*s%s%*s",t_start_str);     
         /* get time as double */
         t_start_value = sscan_time( t_start_str);
         keytype= DRMS_TYPE_TIME;
         key_anyval.time_val= t_start_value;
         status = drms_setkey(record, keyname, keytype, &key_anyval);

         /* if T_HKVALS keyword set ptr_hkval variable to return back */
         if(!strcmp(keyname, "T_HKVALS"))
         {
           /* set filevar_status - got keyword T_HKVALS! */
           filevar_status+=HKLMP_VAR_T_HKVALS;

           /* set hkvals time as TIME variable and as TIME string variable */
           *ptr_hkvals=t_start_value;
           strcpy( hk_time_val_str,t_start_str);

#ifdef DEBUG_LOAD_MP
           printf("load_mp:(15):set_image_loc_values():Setting T_HKVALS keyword to float value to:<%f>\n", *ptr_hkvals);
#endif
         }
         if(!strcmp(keyname, "T_START"))
         {
           /* set filevar_status - got keyword T_START! */
           filevar_status+=HKLMP_VAR_T_START;

           /* set T_START time as TIME variable and as TIME string variable */
           *ptr_tstart=t_start_value;
           strcpy( mp_tstart_str,t_start_str);

#ifdef DEBUG_LOAD_MP
           printf("load_mp:(4):set_image_loc_values():Setting T_START float value using data in image loc file:<%f>\n", *ptr_tstart);
           printf("load_mp:(5):set_image_loc_values():String display of T_START setting<%s>\n", mp_tstart_str);
#endif
         }
       }
       else if(!strcmp(keyname, "VERSION"))
       {
         /* set filevar_status - got keyword VERSION! */
         filevar_status+=HKLMP_VAR_VERSION;

         /* if T_START records exists then use current T_STOP value in record */
         /* before setting version check if record already exists with a version value */
         ret_ver_value=set_mp_version(seriesname, mp_tstart_str);

#ifdef DEBUG_LOAD_MP
         printf("load_mp:(21):set_image_loc_values(): from set_mp_version: <%d>\n",ret_ver_value);
#endif
 
         /* set VERSION keyword value in record in master pointing drms series */
         keytype= DRMS_TYPE_INT;
         key_anyval.int_val= ret_ver_value;
         status = drms_setkey(record, keyname, keytype, &key_anyval);
       }
       else
       {
         /* set the other keyword values set up in the image_location file that define image location data */
         keytype= DRMS_TYPE_FLOAT;
         sscanf(line,"%*s%*s%f%*s",&key_anyval.float_val);     
         status = drms_setkey(record, keyname, keytype, &key_anyval);
       }
    } /* end of else */
  }/* end of while */

  /* close image location file */
  fclose(file_ptr);

  /* set DATE keyword - get CURRENT TIME "now" in year, month, day and hour format */
  tvalue = time(NULL);
  /* get current date-time in UTC */
  time_ptr = gmtime(&tvalue);
  /* initalize time_str to nulls */
  for(int i=0; i < HKLMP_PACKET_TIME_STR;time_str[i]='\0',i++);
  /* convert current time to time_str */
  sprintf(time_str, "%04d.%02d.%02d_%02d:%02d:%02d_UTC",(time_ptr->tm_year+1900), 
  (time_ptr->tm_mon+1), time_ptr->tm_mday, time_ptr->tm_hour, time_ptr->tm_min, time_ptr->tm_sec);
  /* get current time  double value*/
  now_time=sscan_time(time_str);
  //printf("time_str<%s>\n",time_str);
  // printf("current datetime is <%f>\n",now_time);
  /* set drms type and long telemetry name  */
  keytype= DRMS_TYPE_TIME;
  strcpy(keyname, "DATE");
  /* set packet time */
  key_anyval.time_val= now_time;
  /* set record */
  status = drms_setkey(record, keyname, keytype, &key_anyval);

  /* we have finished processing file now check file variable_status for missing key variables in file */
  if(filevar_status == HKLMP_VAR_VERSION+HKLMP_VAR_T_START+HKLMP_VAR_T_STOP+HKLMP_VAR_T_HKVALS)
  {
    return(0);
  }
  else 
  {
    return (filevar_status);
  }
}/* end of  set_image_loc_values() */ 



/**************************************************************************
 * GET HMI LEG STATUS KEYWORDS                                            *
 * Function:get_hmilegstatus_keywords()                                   *          
 * Description:Sets the hmi leg status keyword values listed in image     *
 * location file to keywords in master pointing series by lookimg up      *
 * keywords in the hmi.leg_status DRMS series using the HKVALS time in    *
 * in the image location file.                                            *
 **************************************************************************/
int get_hmilegstatus_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR] )
{
  /* Passed in the pointer to record to write to in master pointing series "record" and the            *
   * TIME value of the HK_VALS "p_hkvals" and the string value of TIME of the HK_VALS "hk_time_val"    *
   * Return back pass or error status.                                                                  */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  DRMS_Type_Value_t key_anyval;
  TIME adj_hkval;
  char adj_hkval_str[HKLMP_PACKET_TIME_STR];
  char nquery[HKLMP_MAX_QUERY_STR];
  char dsname[HKLMP_MAX_DSNAME_STR];
  char series[HKLMP_MAX_DSNAME_STR];
  char *sn;
  float leg_step;
  int status;

  /* initialize variables */
  sn=series;
  rs=NULL; //points to leg_status series
  rec=NULL;//point to leg_status series

  /*get series name for hmi iss status */
  sn = (char *)getenv("HK_LMP_HMI_LEG_SERIESNAME"); 
  if(sn == NULL)
  {
    printkerr("ERROR at %s, line %d: no hmi leg status environment variable set. "
              " Set envirionment HK_LMP_HMI_LEG_SERIESNAME variable in sourced file "
              " SOURCE_ENV_FOR_LOAD_MP.\n", __FILE__,__LINE__);
    return (HKLMP_FAILED_STATUS);
  }

  /* get leg status series name */
  strcpy(dsname,sn);
  strcat(dsname,"\0");

  /* open records for leg_status series and get step size */
  rs = drms_open_records(drms_env, dsname, &status);
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(22):get_hmilegstatus():number of records retrieved was %d\n",rs->n);
#endif
    leg_step=drms_getkey_float(rec, "T_START_step", &status);
    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  else
  {
    printkerr("ERROR at %s, line %d: Could not find records for <%s> series. "
              "Could not get step value needed to load this series values in master pointing series.\n",
              __FILE__,__LINE__,dsname);
    return(HKLMP_FAILED_STATUS);
  }

  /* hkval time passed in arguments */
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(23):get_hmilegstatus_keywords():hkval_str is <%s>\n",hk_time_val);
#endif

  /* adjust hkval by -7200 second for second part of query with T_START in leg_status series*/
  adj_hkval= *p_hkvals - leg_step;

  /* get string value of adj_hkval */
  (void)sprint_time (adj_hkval_str,  adj_hkval, "UTC", 0);
  strcat(adj_hkval_str,"\0");
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(24):get_hmilegstatus_keywords():adj_hkval_str is <%s>\n",adj_hkval_str);
#endif

  /* create query for hmi.leg_status record */
  /* where T_START is time string keyword in hmi.leg_status series */
  /* where hk_time_val is time string keyword from master_pointing series for looking up hk values in  hmi.leg_status series */
  /* where adj_hkval_str is adjust time string value shown above. use to check upper limit by subtracting step from hk_time_val */
  sprintf(nquery,"%s[? $(%s) >= T_START AND $(%s) <  T_START ?]",dsname, hk_time_val, adj_hkval_str);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(25):get_hmilegstatus():nquery is<%s>\n",nquery);
#endif

  /* open record using query*/
  rs = drms_open_records(drms_env, nquery, &status);
  if(!rs->n)
  {
    printkerr("ERROR at %s, line %d: No Records found using query <%s>. "
              "Could not set values from <%s> series in master pointing series.\n", 
              __FILE__,__LINE__,nquery,dsname);
    return(HKLMP_FAILED_STATUS);
  }
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(26):get_hmilegstatus():number of records retrieved was %d\n",rs->n);
#endif
    /*  do getkey and setkey for master pointing series */
    drms_setkey_float(record, "HMI_FSW_AL1_POSITION", drms_getkey_float(rec, "HAL1POS_MEAN", &status));
    drms_setkey_float(record, "HMI_FSW_AL2_POSITION", drms_getkey_float(rec, "HAL2POS_MEAN", &status));
    drms_setkey_float(record, "HMI_AL1_STATUS",drms_getkey_float(rec, "HAL1STAT_MEAN", &status));
    drms_setkey_float(record, "HMI_AL2_STATUS",drms_getkey_float(rec, "HAL2STAT_MEAN", &status));

    /* close leg_status series */
    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  return(HKLMP_PASSED_STATUS);
} /* end of  get_hmilegstatus_keywords()  */



/**************************************************************************
 * GET HMI ISS STATUS KEYWORDS                                            *
 * Function:get_hmiissstatus_keywords()                                   *          
 * Description:Sets the hmi iss status keyword values listed in image     *
 * location file to keywords in master pointing series by lookimg up      *
 * keywords in the hmi.iss_status DRMS series using the HKVALS time in    *
 * in the image location file.                                            *
 **************************************************************************/
int get_hmiissstatus_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR])
{
  /* Passed in the pointer to record to write to in master pointing series "record" and the            *
   * TIME value of the HK_VALS "p_hkvals" and the string value of TIME of the HK_VALS "hk_time_val"    *
   * Return back pass or error status.                                                                  */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  DRMS_Type_Value_t key_anyval;
  TIME adj_hkval;
  char adj_hkval_str[HKLMP_PACKET_TIME_STR];
  char nquery[HKLMP_MAX_QUERY_STR];
  char dsname[HKLMP_MAX_DSNAME_STR];
  char series[HKLMP_MAX_DSNAME_STR];
  float leg_step;
  int status;
  char *sn;
  sn=series;

  /* initialize variables */
  sn=series;
  rs=NULL;
  rec=NULL;

  /*get series name for hmi iss status */
  sn = (char *)getenv("HK_LMP_HMI_ISS_SERIESNAME"); 
  if(sn == NULL)
  {
    printkerr("ERROR at %s, line %d: no hmi iss status environment variable set. "
              " Set envirionment HK_LMP_HMI_ISS_SERIESNAME variable in sourced file "
              " SOURCE_ENV_FOR_LOAD_MP.\n", __FILE__,__LINE__);
    return (HKLMP_FAILED_STATUS);
  }

  /* get leg status series name */
  strcpy(dsname,sn);

  /* open records for iss_status series and get step size */
  rs = drms_open_records(drms_env, dsname, &status);
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(27):get_hmissstatus():number of records retrieved was %d\n",rs->n);
#endif
    leg_step=drms_getkey_float(rec, "T_START_step", &status);
    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  else
  {
    printkerr("ERROR at %s, line %d: Could not find records for <%s> series. "
              "Could not get step value needed to load this series values in master pointing series.\n",
              __FILE__,__LINE__,dsname);
    return(HKLMP_FAILED_STATUS);
  }

  /* hkval time passed in arguments */
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(28):get_hmiissstatus_keywords:hkval_str is <%s>\n",hk_time_val);
#endif

  /* adjust hkval by -7200 second for second part of query with T_START in leg_status series*/
  adj_hkval= *p_hkvals - leg_step;

  /* get string value of adj_hkval */
  (void)sprint_time (adj_hkval_str,  adj_hkval, "UTC", 0);
  strcat(adj_hkval_str,"\0");
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(29):get_hmiissstatus_keywords():adj_hkval_str is <%s>\n",adj_hkval_str);
#endif

  /* create query for hmi.iss_status record */
  /* where T_START is time string keyword in hmi.iss_status series */
  /* where hk_time_val is time string keyword from master_pointing series for looking up hk values in  hmi.iss_status series */
  /* where adj_hkval_str is adjust time string value shown above. use to check upper limit by subtracting step from hk_time_val */
  sprintf(nquery,"%s[? $(%s) >= T_START AND $(%s) <  T_START ?]",dsname,hk_time_val, adj_hkval_str);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(30):get_hmiissstatus():nquery is<%s>\n",nquery);
#endif

  /* open record using query*/
  rs = drms_open_records(drms_env, nquery, &status);
  if(!rs->n)
  {
    printkerr("ERROR at %s, line %d: No Records found using query <%s>. "
              "Could not set values from <%s> series in master pointing series.\n", 
              __FILE__,__LINE__,nquery,dsname);
    return(HKLMP_FAILED_STATUS);
  }
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(31):get_hmissstatus():number of records retrieved was %d\n",rs->n);
#endif
    /* set master pointing */
    drms_setkey_float(record, "HMI_ISS_ERRGAINY",  drms_getkey_float(rec, "HIERRGNY_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_ERRGAINZ",  drms_getkey_float(rec, "HIERRGNZ_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_ERROFFY", drms_getkey_float(rec, "HIERROFY_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_ERROFFZ", drms_getkey_float(rec, "HIERROFZ_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PZTOFFA", drms_getkey_float(rec, "HIPZTOFA_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PZTOFFB",drms_getkey_float(rec, "HIPZTOFB_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PZTOFFC",drms_getkey_float(rec, "HIPZTOFC_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PKT_YCOEF_A",drms_getkey_float(rec, "HIYCOEFA_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PKT_YCOEF_B", drms_getkey_float(rec, "HISYCOEF_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PKT_YCOEF_C", drms_getkey_float(rec, "HIYCOEFC_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PKT_ZCOEF_A", drms_getkey_float(rec, "HIZCOEFA_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PKT_ZCOEF_B", drms_getkey_float(rec, "HIZCOEFB_MEAN", &status));
    drms_setkey_float(record, "HMI_ISS_PKT_ZCOEF_C", drms_getkey_float(rec, "HIZCOEFC_MEAN", &status));

    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  return(HKLMP_PASSED_STATUS);

}/* end of get_hmiissstatus_keywords */



/**************************************************************************
 * GET AIA ISS STATUS 12 KEYWORDS                                         *
 * Function:get_aiaiss12status_keywords()                                 *          
 * Description:Sets the aia iss 12 status keyword values listed in image  *
 * location file to keywords in master pointing series by lookimg up      *
 * keywords in the aia.iss12_status DRMS series using the HKVALS time in  *
 * in the image location file.                                            *
 **************************************************************************/
int get_aiaiss12status_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR])
{
  /* Passed in the pointer to record to write to in master pointing series "record" and the            *
   * TIME value of the HK_VALS "p_hkvals" and the string value of TIME of the HK_VALS "hk_time_val"    *
   * Return back pass or error status.                                                                  */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  DRMS_Type_Value_t key_anyval;
  TIME adj_hkval;
  char adj_hkval_str[HKLMP_PACKET_TIME_STR];
  char nquery[HKLMP_MAX_QUERY_STR];
  char dsname[HKLMP_MAX_DSNAME_STR];
  char series[HKLMP_MAX_DSNAME_STR];
  char *sn;
  float leg_step;
  int status;

  /* initialize variables */
  sn=series;
  rs=NULL;
  rec=NULL;

  /*get series name for hmi iss status */
   sn = (char *)getenv("HK_LMP_AIA_ISS12_SERIESNAME"); 
  if(sn == NULL)
  {
    printkerr("ERROR at %s, line %d: no aia iss12 status environment variable set. "
              " Set envirionment HK_LMP_AIA_ISS12_SERIESNAME variable in sourced file "
              " SOURCE_ENV_FOR_LOAD_MP.\n", __FILE__,__LINE__);
    return (HKLMP_FAILED_STATUS);
  }

  /* get leg status series name */
  strcpy(dsname,sn);

  /* open records for iss_status series and get step size */
  rs = drms_open_records(drms_env, dsname, &status);
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(32):get_aiaiss12status():number of records retrieved was %d\n",rs->n);
#endif
    leg_step=drms_getkey_float(rec, "T_START_step", &status);
    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  else
  {
    printkerr("ERROR at %s, line %d: Could not find records for <%s> series. "
              "Could not get step value needed to load this series values in master pointing series.\n",
              __FILE__,__LINE__,dsname);
    return(HKLMP_FAILED_STATUS);
  }

  /* hkval time passed in arguments */
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(33):get_aiaiss12sstatus_keywords():hkval_str is <%s>\n",hk_time_val);
#endif

  /* adjust hkval by -7200 second for second part of query with T_START in leg_status series*/
  adj_hkval= *p_hkvals - leg_step;

  /* get string value of adj_hkval */
  (void)sprint_time (adj_hkval_str,  adj_hkval, "UTC", 0);
  strcat(adj_hkval_str,"\0");
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(34):get_aiaiss12status_keywords():adj_hkval_str is <%s>\n",adj_hkval_str);
#endif

  /* create query for aia.iss_1_2_status record */
  /* where T_START is time string keyword in aia.iss_1_2_status series */
  /* where hk_time_val is time string keyword from master_pointing series for looking up hk values in  aia.iss_1_2_status series */
  /* where adj_hkval_str is adjust time string value shown above. use to check upper limit by subtracting step from hk_time_val */
  sprintf(nquery,"%s[? $(%s) >= T_START AND $(%s) <  T_START ?]",dsname,hk_time_val, adj_hkval_str);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(35):get_aiaiss12status():nquery is<%s>\n",nquery);
#endif

  /* open record using query*/
  rs = drms_open_records(drms_env, nquery, &status);
  if(!rs->n)
  {
    printkerr("ERROR at %s, line %d: No Records found using query <%s>. "
              "Check the HKVALS value selected has value in aia.iss_1_2_status series. "
              "Could not set values from <%s> series in master pointing series.\n", 
              __FILE__,__LINE__,nquery,dsname);
    return(HKLMP_FAILED_STATUS);
  }
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(36):get_aiaiss12status():number of records retrieved was %d\n",rs->n);
#endif

    /* set master pointing by getting keywords values and  setting master pointing series */
    drms_setkey_float(record, "AIA_IS1_ERRGAINY", drms_getkey_float(rec, "A1ERRGNY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_ERRGAINZ",drms_getkey_float(rec, "A1ERRGNZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_ERROFFY", drms_getkey_float(rec, "A1ERROFY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_ERROFFZ", drms_getkey_float(rec, "A1ERROFZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_PZTGAINA", drms_getkey_float(rec, "A1PZTGNA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_PZTGAINB", drms_getkey_float(rec, "A1PZTGNB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_PZTGAINC", drms_getkey_float(rec, "A1PZTGNC_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_PZTOFFA", drms_getkey_float(rec, "A1PZTOFA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_PZTOFFB", drms_getkey_float(rec, "A1PZTOFB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS1_PZTOFFC", drms_getkey_float(rec, "A1PZTOFC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT1_PKT_YCOEF_A",drms_getkey_float(rec, "AGT1_YCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT1_PKT_YCOEF_B",drms_getkey_float(rec, "AGT1_YCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT1_PKT_YCOEF_C",drms_getkey_float(rec, "AGT1_YCC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT1_PKT_ZCOEF_A",drms_getkey_float(rec, "AGT1_ZCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT1_PKT_ZCOEF_B",drms_getkey_float(rec, "AGT1_ZCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT1_PKT_ZCOEF_C",drms_getkey_float(rec, "AGT1_ZCC_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_ERRGAINY",drms_getkey_float(rec, "A2ERRGNY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_ERRGAINZ",drms_getkey_float(rec, "A2ERRGNZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_ERROFFY", drms_getkey_float(rec, "A2ERROFY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_ERROFFZ",drms_getkey_float(rec, "A2ERROFZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_PZTGAINA", drms_getkey_float(rec, "A2PZTGNA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_PZTGAINB",drms_getkey_float(rec, "A2PZTGNB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_PZTGAINC",drms_getkey_float(rec, "A2PZTGNC_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_PZTOFFA",drms_getkey_float(rec, "A2PZTOFA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_PZTOFFB",drms_getkey_float(rec, "A2PZTOFB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS2_PZTOFFC",drms_getkey_float(rec, "A2PZTOFC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT2_PKT_YCOEF_A", drms_getkey_float(rec, "AGT2_YCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT2_PKT_YCOEF_B",drms_getkey_float(rec, "AGT2_YCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT2_PKT_YCOEF_C",drms_getkey_float(rec, "AGT2_YCC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT2_PKT_ZCOEF_A",drms_getkey_float(rec, "AGT2_ZCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT2_PKT_ZCOEF_B",drms_getkey_float(rec, "AGT2_ZCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT2_PKT_ZCOEF_C",drms_getkey_float(rec, "AGT2_ZCC_MEAN", &status));

    /* close aiaiss12status records */
    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  return(HKLMP_PASSED_STATUS);
}/* end of  get_aiaiss12status_keywords */



/**************************************************************************
 * GET AIA ISS STATUS 34 KEYWORDS                                         *
 * Function:get_aiaiss34status_keywords()                                 *          
 * Description:Sets the aia iss 34 status keyword values listed in image  *
 * location file to keywords in master pointing series by lookimg up      *
 * keywords in the aia.iss34_status DRMS series using the HKVALS time in  *
 * in the image location file.                                            *
 **************************************************************************/
int get_aiaiss34status_keywords(DRMS_Record_t *record, TIME *p_hkvals,char hk_time_val[HKLMP_PACKET_TIME_STR])
{
  /* Passed in the pointer to record to write to in master pointing series "record" and the            *
   * TIME value of the HK_VALS "p_hkvals" and the string value of TIME of the HK_VALS "hk_time_val"    *
   * Return back pass or error status                                                                  */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  DRMS_Type_Value_t key_anyval;
  TIME adj_hkval;
  char adj_hkval_str[HKLMP_PACKET_TIME_STR];
  char nquery[HKLMP_MAX_QUERY_STR];
  char dsname[HKLMP_MAX_DSNAME_STR];
  char series[HKLMP_MAX_DSNAME_STR];
  char *sn;
  float leg_step;
  int status;

  /* initialize variables */
  sn=series;
  rs=NULL;
  rec=NULL;

  /*get series name for hmi iss status */
  sn = (char *)getenv("HK_LMP_AIA_ISS34_SERIESNAME"); 
  if(sn == NULL)
  {
    printkerr("ERROR at %s, line %d: no aia iss34 status environment variable set. "
              " Set envirionment HK_LMP_AIA_ISS34_SERIESNAME variable in sourced file "
              " SOURCE_ENV_FOR_LOAD_MP.\n", __FILE__,__LINE__);
    return (HKLMP_FAILED_STATUS);
  }

  /* get leg status series name */
  strcpy(dsname,sn);

  /* open records for iss_status series and get step size */
  rs = drms_open_records(drms_env, dsname, &status);
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(37):get_aiaiss34status():number of records retrieved was %d\n",rs->n);
#endif
    leg_step=drms_getkey_float(rec, "T_START_step", &status);
    drms_close_records(rs, DRMS_FREE_RECORD);
  }
  else
  {
    printkerr("ERROR at %s, line %d: Could not find records for <%s> series. "
              "Could not get step value needed to load this series values in master pointing series.\n",
              __FILE__,__LINE__,dsname);
    return(HKLMP_FAILED_STATUS);
  }

  /* hkval time passed in arguments */
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(38):get_aiaiss34sstatus_keywords():hkval_str is <%s>\n",hk_time_val);
#endif

  /* adjust hkval by -7200 second or -2 hours(use leg_step!) for second part of query with T_START in leg_status series*/
  adj_hkval= *p_hkvals - leg_step;

  /* get string value of adj_hkval */
  (void)sprint_time (adj_hkval_str,  adj_hkval, "UTC", 0);
  strcat(adj_hkval_str,"\0");
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(39):get_aiaiss34status_keywords():adj_hkval_str is <%s>\n",adj_hkval_str);
#endif

  /* create query for aia.iss_3_4_status record */
  /* where T_START is time string keyword in aia.iss_3_4_status series */
  /* where hk_time_val is time string keyword from master_pointing series for looking up hk values in  aia.iss_3_4_status series */
  /* where adj_hkval_str is adjust time string value shown above. use to check upper limit by subtracting step from hk_time_val */
  sprintf(nquery,"%s[? $(%s) >= T_START AND $(%s) <  T_START ?]",dsname,hk_time_val, adj_hkval_str);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(40):get_aiaiss34status():nquery is<%s>\n",nquery);
#endif

  /* open record using query*/
  rs = drms_open_records(drms_env, nquery, &status);
  if(!rs->n)
  {
    printkerr("ERROR at %s, line %d: No Records found using query <%s>. "
              "Could not set values from <%s> series in master pointing series.\n", 
              __FILE__,__LINE__,nquery,dsname);
    return(HKLMP_FAILED_STATUS);
  }
  if(rs->n)
  {
    rec=rs->records[0];
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(41):get_aiaiss34status():number of records retrieved was %d\n",rs->n);
#endif

    /* get keyword values and set master pointing series */
    drms_setkey_float(record, "AIA_IS3_ERRGAINY", drms_getkey_float(rec, "A3ERRGNY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_ERRGAINZ", drms_getkey_float(rec, "A3ERRGNZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_ERROFFY", drms_getkey_float(rec, "A3ERROFY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_ERROFFZ", drms_getkey_float(rec, "A3ERROFZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_PZTGAINA", drms_getkey_float(rec, "A3PZTGNA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_PZTGAINB", drms_getkey_float(rec, "A3PZTGNB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_PZTGAINC", drms_getkey_float(rec, "A3PZTGNC_MEAN",&status));
    drms_setkey_float(record, "AIA_IS3_PZTOFFA", drms_getkey_float(rec, "A3PZTOFA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_PZTOFFB", drms_getkey_float(rec, "A3PZTOFB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS3_PZTOFFC", drms_getkey_float(rec, "A3PZTOFC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT3_PKT_YCOEF_A", drms_getkey_float(rec, "AGT3_YCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT3_PKT_YCOEF_B", drms_getkey_float(rec, "AGT3_YCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT3_PKT_YCOEF_C", drms_getkey_float(rec, "AGT3_YCC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT3_PKT_ZCOEF_A", drms_getkey_float(rec, "AGT3_ZCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT3_PKT_ZCOEF_B", drms_getkey_float(rec, "AGT3_ZCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT3_PKT_ZCOEF_C", drms_getkey_float(rec, "AGT3_ZCC_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_ERRGAINY", drms_getkey_float(rec, "A4ERRGNY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_ERRGAINZ", drms_getkey_float(rec, "A4ERRGNZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_ERROFFY", drms_getkey_float(rec, "A4ERROFY_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_ERROFFZ", drms_getkey_float(rec, "A4ERROFZ_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_PZTGAINA", drms_getkey_float(rec, "A4PZTGNA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_PZTGAINB", drms_getkey_float(rec, "A4PZTGNB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_PZTGAINC", drms_getkey_float(rec, "A4PZTGNC_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_PZTOFFA", drms_getkey_float(rec, "A4PZTOFA_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_PZTOFFB", drms_getkey_float(rec, "A4PZTOFB_MEAN", &status));
    drms_setkey_float(record, "AIA_IS4_PZTOFFC", drms_getkey_float(rec, "A4PZTOFC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT4_PKT_YCOEF_A", drms_getkey_float(rec, "AGT4_YCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT4_PKT_YCOEF_B", drms_getkey_float(rec, "AGT4_YCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT4_PKT_YCOEF_C", drms_getkey_float(rec, "AGT4_YCC_MEAN", &status));
    drms_setkey_float(record, "AIA_GT4_PKT_ZCOEF_A", drms_getkey_float(rec, "AGT4_ZCA_MEAN", &status));
    drms_setkey_float(record, "AIA_GT4_PKT_ZCOEF_B", drms_getkey_float(rec, "AGT4_ZCB_MEAN", &status));
    drms_setkey_float(record, "AIA_GT4_PKT_ZCOEF_C", drms_getkey_float(rec, "AGT4_ZCC_MEAN", &status));

    /* close records for aia.iss34_status series */
    drms_close_records(rs, DRMS_FREE_RECORD);
  }

  return(HKLMP_PASSED_STATUS);

} /*end of get_iss34status*/

/***************************************************************************
 * UPDATE T_STOP AND VERSION                                               *
 * Function:update_tstop_version()                                         *
 * Description: Updates the T_STOP and VERSION keyword of a record that use*
 * to be the lastest record in master pointing series. If this is true,    *
 * need to update the VERSION keyword from 0 to 1. If this is true, need to*
 * update the T_STOP time from "far out" value to T_START value of the     *
 * latest record in the master pointing series.                            *
 ***************************************************************************/
void update_tstop_version(char query[HKLMP_MAX_QUERY_STR], char mp_tstart_value[HKLMP_PACKET_TIME_STR])
{
  DRMS_RecordSet_t *old_rs;
  DRMS_Record_t *old_rec;
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  TIME new_tstop;
  TIME prev_tstart;
  char keyname[HKLMP_MAX_KEYWORD_NAME_STR];
  char nquery[HKLMP_MAX_QUERY_STR];
  int status;

  /* check if there is a record that use to be the latest record in series */
  sprintf(nquery,"%s[? T_START < $(%s) AND VERSION = 0 ?]",query,mp_tstart_value);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(43):update_tstop_version(): query is:<%s> mp_tstart_value:<%s>\n",query, mp_tstart_value); 
  printf("load_mp:(44):update_tstop_version():nquery is<%s>\n",nquery);
#endif

  /* open record using query*/
  old_rs = drms_open_records(drms_env, nquery, &status);

  if (status)
  {
    printkerr("ERROR at %s, line %d:update_tstop_version(): Cannot open record using this data"
              " series name:<%s>.Existing program.\n",__FILE__,__LINE__, query);
    status = drms_close_records(old_rs, DRMS_FREE_RECORD);
    return;
  }
  else
  {
    ;//printf("Sucessfully status for open record for series <%s>\n", query);
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(45):update_tstop_version():Sucessfully created record for series <%s>\n", query);
#endif
  }
  
  /* check if got 1 record that use to be the latest record in series */
  if (old_rs->n == 1)
  {
     old_rec = old_rs->records[0];
#ifdef DEBUG_LOAD_MP
     printf("load_mp:(46):update_tstop_version():in if:got %d record\n",old_rs->n);
#endif
  }
  else if  (old_rs->n == 0)
  {
     /* did not find any record that use to be latest-so no need to update old record */
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(47):update_tstop_version():if 0 records then no need to update old T_START Record. note got <%d> records\n",old_rs->n);
#endif
    status = drms_close_records(old_rs, DRMS_FREE_RECORD);
    return;
  }
  else
  {
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(48):update_tstop_version():in else:got not 1 or 0 records. This should not occur. Got <%d> records\n",old_rs->n);
#endif
    printkerr("ERROR at %s, line %d:update_tstop_version(): Got more than one record with"
              " verion equal to 0 for query:<%s>.Returning.\n",__FILE__,__LINE__, nquery);
    status = drms_close_records(old_rs, DRMS_FREE_RECORD);
    return;
  }
  
  /* got 1 record to reach here - so close query check */
  /* get old T_START value in mp series that was previously set to version 0 and then close*/
  prev_tstart= drms_getkey_time(old_rec, "T_START", &status);
  
  /* create record to write to since found 1 record. This record with overwrite the old rec with updates to T_STOP time and version */
  rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);

  /* set master pointing series rec to pass to following function to set master pointing series */
  rec = rs->records[0];

  /*copy all values in old rec to new rec- note using set class=kDRMS_KeyClass_All for drms_copykeys function */
  (void)drms_copykeys(rec, old_rec, 1, kDRMS_KeyClass_All);
 
  /* now set items that need new values like T_START, T_STOP, and VERSION */
  /* set previous T_START value set to version = 0 */
  strcpy(keyname, "T_START");
  keytype= DRMS_TYPE_TIME;
  key_anyval.time_val= prev_tstart;
  status = drms_setkey(rec, keyname, keytype, &key_anyval);

  /* set T_STOP to new T_START value */
  /* get time as double using T_START time in new just updated record */
  new_tstop = sscan_time( mp_tstart_value);
  /* better name for  mp_tstart_value is  mp_tstart_str*/
#ifdef DEBUG_LOAD_MP
  printf ("load_mp:(49):update_tstop_version:mp_tstart_value<%s>\n", mp_tstart_value);
  printf ("load_mp:(50):update_tstop_version:new_tstop_float<%f>\n", new_tstop);
#endif
  strcpy(keyname, "T_STOP");
  keytype= DRMS_TYPE_TIME;
  key_anyval.time_val= new_tstop;
  status = drms_setkey(rec, keyname, keytype, &key_anyval);

  /* set old rec VERSION to 1 */
  strcpy(keyname, "VERSION");
  keytype= DRMS_TYPE_INT;
  key_anyval.int_val= 1;
  status = drms_setkey(rec, keyname, keytype, &key_anyval);

  /* COMMIT values for "new" old record to master pointing - close record */ 
  status = drms_close_records(rs, DRMS_INSERT_RECORD);
  if (status != 0)
  {
    printkerr("ERROR at %s, line %d: Cannot close drms record. Status ret:<%d>  "
              "Failed to update rec that was the latest rec, with new T_STOP and VERSION values.\n",
               __FILE__,__LINE__,status);
  }
  else
  {
     ;/*added to use else case to show debug message */
#ifdef DEBUG_LOAD_MP
      printf("load_mp:(51):update_tstop_version(): Completed close of write of record for master pointing series\n");
#endif
  }

  /* close record used to copy values from */
  status = drms_close_records(old_rs, DRMS_FREE_RECORD);
  if (status != 0)
  {
    printkerr("ERROR at %s, line %d: Cannot close drms record. Status ret:<%d>\n",
               __FILE__,__LINE__,status);
  }

}/* end of update_tstop_version() */



/***************************************************************************
 * SET MASTER POINTING VERSION                                             *
 * Function:set_mp_version()                                               *
 * Description: Sets the VERSION keyword. Check if record exists, if exists*
 * and not last record, increment VERSION keyword by one.                  *
 ***************************************************************************/
int set_mp_version(char *sname, char tstart_str[HKLMP_PACKET_TIME_STR])
{
  /* Pass in master pointing series name sname and T_START Time string value tstart_str*/
  DRMS_RecordSet_t *exists_rs;
  DRMS_Record_t *exists_rec;
  int status;
  int curr_ver;
  char rec_exists_query[HKLMP_MAX_QUERY_STR];

#ifdef DEBUG_LOAD_MP
  printf("load_mp:(16):set_mp_version(): Check if not last record and if there increment VERSION.\n");
#endif

  /*build query to find record that has tstart value and has version equal to 1 or greater*/
  sprintf(rec_exists_query,"%s[? T_START=$(%s) AND VERSION>=1 ?]",sname,tstart_str);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(17):set_mp_version():Check rec exists with this query:<%s>\n", rec_exists_query);
#endif
    
  exists_rs = drms_open_records(drms_env, rec_exists_query, &status);
  if(exists_rs->n > 0)
  {
     /* find record therefore get increment value and increment by 1 and return value */
     exists_rec=exists_rs->records[0];
     curr_ver=drms_getkey_int(exists_rec, "VERSION", &status);
     status = drms_close_records(exists_rs, DRMS_FREE_RECORD);
     return(++curr_ver);
  }/* end outer if */
  else
  {
     status = drms_close_records(exists_rs, DRMS_FREE_RECORD);
     /* this case can be for latest record or for new record(Have not done) */
     /* check for other version 0 records and if this is latest or new record */
     sprintf(rec_exists_query,"%s[? T_START>$(%s) AND VERSION=0 ?]",sname,tstart_str);
     exists_rs = drms_open_records(drms_env, rec_exists_query, &status);
#ifdef DEBUG_LOAD_MP
     printf("load_mp:(18):set_mp_version():check rec exists using this query:<%s>\n",rec_exists_query);
#endif
     if(exists_rs->n > 0)
     {  
       /* this is a new record since there is one record in series greater than t_start_str */
#ifdef DEBUG_LOAD_MP
       printf("load_mp:(19):set_mp_version():Rec exists:<%s>\n", rec_exists_query);
#endif
       /*then got record with version=0 and is probably latest record, so this is new rec*/
       status = drms_close_records(exists_rs, DRMS_FREE_RECORD);
       return(1);
     }
     else
     {
       /*this is latest record if found no records greater than tstart_str, so set version=0*/
#ifdef DEBUG_LOAD_MP
       printf("load_mp:(20):set_mp_version():Rec does not exists:<%s>\n", rec_exists_query);
#endif
       status = drms_close_records(exists_rs, DRMS_FREE_RECORD);
       return(0);
     }
  }/* end outer else */
}/* end of set_mp_version() */



/***************************************************************************
 * SET T_STOP VALUE                                                        *
 * Function:set_tstop_value()                                              *
 * Description: Returns closest records T_START time to use to set the     *
 * the T_STOP value by the calling function.                               *
 ***************************************************************************/
TIME set_tstop_value(char *sname, char tstart_str[HKLMP_PACKET_TIME_STR],TIME *ptstart)
{
  /* Pass in the master pointing series name sname, the T_START string time value, and the T_START TIME value */
  DRMS_RecordSet_t *exists_rs;
  DRMS_RecordSet_t *find_rs;
  DRMS_Record_t *exists_rec;
  DRMS_Record_t *find_rec;
  TIME curr_tstop;
  TIME find_closest;
  TIME find_current;
  char find_rec_query[HKLMP_MAX_QUERY_STR];
  char rec_exists_query[HKLMP_MAX_QUERY_STR];
  int find_rs_num;
  int i;
  int status;

#ifdef DEBUG_LOAD_MP
  printf("load_mp:(6)set_tstop_value():Setting T_STOP value\n");
  printf("load_mp:(7)set_tstop_value():tstart value passed is <%s>\n", tstart_str);
#endif

  /*build query  to search for all records that are not last record(VERSION=0) and have T_START time equal to passed in t_start_str */
  sprintf(rec_exists_query,"%s[? T_START=$(%s) AND VERSION>=1 ?]",sname,tstart_str);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(8)set_tstop_value() CASE 1-if rec not last but record there:check if rec exists using query:<%s>\n", rec_exists_query);
#endif
    
  exists_rs = drms_open_records(drms_env, rec_exists_query, &status);
  if(exists_rs->n > 0)
  {
     /* found rec */
     exists_rec=exists_rs->records[0];

     /* get tstop value */
     curr_tstop=drms_getkey_time(exists_rec, "T_STOP", &status);
#ifdef DEBUG_LOAD_MP
     printf("load_mp:(9)set_tstop_value() Found rec for above Case 1-Retrieve T_STOP value from record and return back value:%f\n",curr_tstop);
#endif
    
     /* close exists_rs */
     status = drms_close_records(exists_rs, DRMS_FREE_RECORD);

     /* return back t_stop value */
     return(curr_tstop);
  }
  else
  {
      /* this is case where record is not latest record with version = 0 *
       * and  is new record and t_start is less than all records.        *
       * find the closest records T_START time and return                */
      /* close exists_rs */
      status = drms_close_records(exists_rs, DRMS_FREE_RECORD);

      /*build query to get all records with T_START time greater than passed in tstart_str to find closest time*/
      sprintf(find_rec_query,"%s[? T_START>$(%s) AND VERSION>=0 ?]",sname,tstart_str); /*VERSION=0 for 0805 rec!*/
      find_rs = drms_open_records(drms_env, find_rec_query, &status);

      /*get record set number retrieved */
      find_rs_num=find_rs->n;
      
#ifdef DEBUG_LOAD_MP
      printf("load_mp:(10)set_tstop_value():CASE 2- if rec not last but less than all T_START times. Retrieve T_STOP value using this query:<%s>\n", find_rec_query);
      printf("load_mp:(11)set_tstop_value():Rec(s) found for above query. if zero use file value. find_rs_num:<%d>\n", find_rs_num);
      printf("load_mp:(12)set_tstop_value():new tstart value:ptstart:<%f>\n", *ptstart);
#endif

      /* loop thru records and find record with closet tstart value */
      for(i=0, find_closest=0.0;i<find_rs_num;i++)
      {
        find_rec=find_rs->records[i];
#ifdef DEBUG_LOAD_MP
        printf("load_mp:(13)set_tstop_value:Found rec for above CASE 2: T_START TIME got <%f>\n",drms_getkey_time(find_rec, "T_START", &status));
#endif
        find_current=drms_getkey_time(find_rec, "T_START", &status);
        if(i == 0)
        {
          /*init to first value */
          find_closest=  find_current;
        }
        else
        {
           if( (fabs(find_closest) - fabs(*ptstart)) > (fabs(find_current) - fabs(*ptstart)) ) 
           {
                find_closest=  find_current;
           }
        }/*else*/
      }
      status = drms_close_records(find_rs, DRMS_FREE_RECORD);
#ifdef DEBUG_LOAD_MP
      printf("load_mp:(14):set_tstop_value:closest value:find_closest:<%f>\n", find_closest);
#endif

      /* return closest T_START time to use later to set the tstop value*/
      return(find_closest);
  }
}/*end  set_tstop_value() */



/***************************************************************************
 * CHECK OLD RECORD T_STOP TIME                                            *
 * Function:check_oldrec_tstop_time()                                      *
 * Description: Checks if the previous record needs to be updated with new *
 * new values for T_STOP and VERSION. Create new previous record and copies*
 * most of the value in the previous record to new record and updates      *
 * T_STOP and VERSION values.                                              *
 ***************************************************************************/
void check_oldrec_tstop_time(char query[HKLMP_MAX_QUERY_STR], char mp_tstart_value[HKLMP_PACKET_TIME_STR])
{
  /* Pass in master pointing series name query and master pointing T_START time of new record updated or added */
  DRMS_RecordSet_t *old_rs;
  DRMS_Record_t *old_rec;
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  TIME mp_tstop_value;
  TIME new_new_tstop;//fix this new_new
  TIME new_tstop;
  TIME prev_tstart;
  char keyname[HKLMP_MAX_KEYWORD_NAME_STR];
  char nquery[HKLMP_MAX_QUERY_STR];
  int status;
  int old_rec_version;
  char new_new_tstop_str[HKLMP_PACKET_TIME_STR];

  /* get t_stop */
  sprintf(nquery,"%s[? T_START = $(%s) ?]",query,mp_tstart_value);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(52):check_oldrec_tstop_time:query is:<%s> mp_tstart_value:<%s>\n",query, mp_tstart_value); 
  printf("load_mp:(53):check_oldrec_tstop_time:nquery is<%s>\n",nquery);
#endif

  /* open record using query*/
  rs = drms_open_records(drms_env, nquery, &status);
  rec=rs->records[0];

  /* set T_STOP time to T_START value of new record just loaded */
  new_new_tstop= drms_getkey_time(rec, "T_STOP", &status);

  /* convert to string time from double time for debugging */
#ifdef DEBUG_LOAD_MP
  (void)sprint_time (new_new_tstop_str,  new_new_tstop, "UTC", 0);
  printf("load_mp:(54):new_new_tstop_str is %s\n", new_new_tstop_str);
#endif
 

  /* close connection */
  status = drms_close_records(rs, DRMS_FREE_RECORD);

  /* create query to find  previous record to use to get most values for creating new previous record */
  sprintf(nquery,"%s[? T_START < $(%s) AND T_STOP = %.0f ?]",query,mp_tstart_value,  new_new_tstop);
#ifdef DEBUG_LOAD_MP
  printf("load_mp:(55):check_oldrec_tstop_time:query is:<%s> mp_tstart_value:<%s>\n",query, mp_tstart_value); 
  printf("load_mp:(56):check_oldrec_tstop_time:nquery is<%s>\n",nquery);
#endif

  /* open old or previous record using query to update tstop time */
  old_rs = drms_open_records(drms_env, nquery, &status);

  if (status)
  {
    printkerr("ERROR at %s, line %d:check_oldrec_tstop_time(): Cannot create record using this data"
              " series name:<%s>.Existing program.\n",__FILE__,__LINE__, query);
    return;
  }
  
  /* check if got 1 record */
  if (old_rs->n == 1)
  {
     old_rec = old_rs->records[0];
#ifdef DEBUG_LOAD_MP
     printf("load_mp:(57):check_oldrec_tstop_time():in if:got %d record\n",old_rs->n);
#endif
  }
  else if  (old_rs->n == 0)
  {
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(58):check_oldrec_tstop_time():if 0 records then no need to update old T_START Record. note got <%d> records\n",old_rs->n);
#endif
    return;
  }
  else
  {
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(59):check_oldrec_tstop_time():in else-check_oldrec_tstop_time:got not 1 or 0 records. Got <%d> records\n",old_rs->n);
#endif
    printkerr("ERROR at %s, line %d:check_oldrec_tstop_time(): Got more than one record with"
              " verion equal to 0 for query:<%s>.Returning.\n",__FILE__,__LINE__, nquery);
    return;
  }
  
  /* got 1 record to reach here - so close query check */
  /* get old T_START value in mp series that was previously set to version 0 and then close*/
  prev_tstart= drms_getkey_time(old_rec, "T_START", &status);
  
  /* create record to write to since found 1 record less then new T_START */
  rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);

  /* set master pointing series rec to pass to following function to set master pointing series */
  rec = rs->records[0];

  /*copy all values in old rec to new rec to create new previous rec with most of the old rec values- set class=kDRMS_KeyClass_All*/
  (void)drms_copykeys(rec, old_rec, 1, kDRMS_KeyClass_All);
 
  /* now set items that need new values like T_START, T_STOP, and VERSION */
  /* set T_START */
  strcpy(keyname, "T_START");
  keytype= DRMS_TYPE_TIME;
  key_anyval.time_val= prev_tstart;
  status = drms_setkey(rec, keyname, keytype, &key_anyval);

  /* set T_STOP to new T_START value */
  /* get time as double */
  new_tstop = sscan_time( mp_tstart_value);
  /* better name for  mp_tstart_value is  mp_tstart_str*/
#ifdef DEBUG_LOAD_MP
  printf ("load_mp:(60):check_oldrec_tstop_time():mp_tstart_value<%s>\n", mp_tstart_value);
  printf ("load_mp:(61):check_oldrec_tstop_time():new_tstop_float<%f>\n", new_tstop);
#endif
  strcpy(keyname, "T_STOP");
  keytype= DRMS_TYPE_TIME;
  key_anyval.time_val= new_tstop;
  status = drms_setkey(rec, keyname, keytype, &key_anyval);

  /* set VERSION:check version value of old_rec, increment value, and set new VERSION value in rs's rec */
  old_rec_version= drms_getkey_time(old_rec, "VERSION", &status);
  strcpy(keyname, "VERSION");
  keytype= DRMS_TYPE_INT;
  key_anyval.int_val= old_rec_version + 1;
  status = drms_setkey(rec, keyname, keytype, &key_anyval);

  /* COMMIT values to master pointing - close record and check status */
  status = drms_close_records(rs, DRMS_INSERT_RECORD);
  if (status != 0)
  {
    printkerr("ERROR at %s, line %d: Cannot close rs drms record. Status ret:<%d>\n",
               __FILE__,__LINE__,status);
  }
  else
  {
    ;/* print debug message */
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(62):check_oldrec_tstop_time():Successfully inserted update previous record and closed rs record. status <%d>\n",status);
#endif
  }

  /* close old_rs  and check status */
  status = drms_close_records(old_rs, DRMS_FREE_RECORD);
  if (status != 0)
  {
    printkerr("ERROR at %s, line %d: Cannot close old_rs drms record. Status ret:<%d>\n",
               __FILE__,__LINE__,status);
  }
  else
  {
    ;/* print debug message */
#ifdef DEBUG_LOAD_MP
    printf("load_mp:(63):check_oldrec_tstop_time():Successfully closed old_rs record. status <%d>\n",status);
#endif
  }
} /* end of check_oldrec_tstop_time() */
