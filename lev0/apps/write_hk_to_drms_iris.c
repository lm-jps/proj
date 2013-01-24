#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/lev0/apps/write_hk_to_drms_iris.c,v 1.1 2013/01/24 19:02:33 jim Exp $"
/*****************************************************************************
 * Filename: write_hk_to_drms.c                                              *
 * Author: Carl                                                              *
 * Create Date: February, 2, 2008                                            *
 * Description: This file contains modules to write housekeeping keywords    *
 *              to DRMS.                                                     *
 * (C) Stanford University, 2008                                             *
 ****************************************************************************/
#include "jsoc_main.h"
#include "packets.h"
#include "decode_hk_vcdu.h" 
#include "decode_hk.h" 
#include "write_hk_to_drms.h"
#include "printk.h"
#include "timeio.h"
/***********************   Protoype Functions  *******************************/
int   write_hk_to_drms(DRMS_Record_t *record, CCSDS_Packet_t **ccsds_pkt);
static TIME  SDO_to_DRMS_time(int sdo_s, int sdo_ss);
char *get_packet_version_number( HK_Keyword_t *kw);
char *lookup_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t **jmap, char *packet_version_number);
void load_map_data(int apid, JSOC_Version_Map_t  *jm, char* packet_version_number);
int find_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t  *jm, char *dsn);
int   get_pkt_time_from_timecodes(HK_Keyword_t *hk, TIME *ptime);
int   check_for_apid(int apid, JSOC_Version_Map_t  *jm);
void  free_jsvn_map( JSOC_Version_Map_t  *top_jm);
char *get_data_packet_name(int apid ) ;
char *get_hk_source(CCSDS_Packet_t *ccsds_ptr, TIME *pt);
int  check_env_variable();
char *get_lookup_filename(int apid ) ;
char * make_strupr (char *in_string);
int check_hk_record_exists(char* ds_name, HK_Keyword_t *kw, int apid);
short create_low_range(HK_Keyword_t *kw_head, long int *low_range);
short check_ds_cache_init(char *dsname); 
short check_curr_map_file(char* pv_num, short* current_map_flag);
int get_query_range(int range_type, HK_Keyword_t *kw, char *qr);
int initialize_timecodes_cache(char* ds_name, HK_Keyword_t *kw, int apid);
int check_hk_record_within_time_range( HK_Keyword_t *kw);

/***********************   Extern Functions  *******************************/
extern int  get_day_from_pkttime(double tc_sec);
extern int  get_month_from_pkttime(double tc_sec);
extern int  get_yr_from_pkttime(double tc_sec);
extern int  check_for_sdo_apid(int apid);
extern int  get_day_from_pkttime(double ptime);
extern int  get_hour_from_pkttime(double ptime);
extern int  get_month_from_pkttime(double ptime);
extern int  get_yr_from_pkttime(double ptime);

/***********************   global variable  *******************************/
//extern char *dtime;
extern char ispquery[];
extern DRMS_RecordSet_t *RSISP, *RSISPTO;
extern unsigned int fsnx;

HK_DSN_RANGE_t *dsr_head=NULL;
char IRIS_ISP_DS[80];	//dataset name. set by lookup_data_series_name()
static DRMS_RecordSet_t *rs = 0;



/*****************************************************************************
 * Name:        WRITE HK TO DRMS                                             *
 * Function :   int write_hk_to_drms(DRMS_Record_t *record,                  *
 *                                   CCSDS_Packet_t **ccsds_pkt)             *
 * Description: Receives link list of CCSDS_Packet_t structures and for each *
 *              CCSDS_Packet_t node, checks HK_Keyword_t link list for the   *
 *              packet_version_number, create data series name based on apid *
 *              and value looked up in the PVN-TO-JSVN map files based apid  *
 *              and packet_version_number and finally used data series name  *
 *              to loop thru HK_Keyword_t structure and write keyword long   *
 *              names and keyword values to Level 0 by APID data series name.*
 *              Or write keyword short names and keyword values to Level 0   *
 *              data series.
 *****************************************************************************/
int write_hk_to_drms(DRMS_Record_t *record, CCSDS_Packet_t **ccsds_pkt)
{
  /* declaration to set drms record */
  DRMS_Type_t keytype;
  DRMS_Type_Value_t key_anyval;
  DRMS_Record_t *rec;
  TIME pkt_time;
  char keyname[HK_LEV0_MAX_LONG_KW_NAME];
  char query[HK_LEV0_MAX_DSNAME_STR];
  char pkt_ver_num[HK_LEV0_MAX_PVN_STR];
  char *pvn, *pvtmp;
  char *directory ;
  char *suffix_filename ;
  int  status, env_status;
  int  rec_alreadycreated_flag;
  int ck_rec_status;
  int ck_rwtr_status;
  static short haveNotChecked=1;

  /* CCSDS and keyword variables and structure to get 
     value to use to set in drms */
  CCSDS_Packet_t *ccsds;
  HK_Keyword_t *kw;

  /* global pointer to structure containing map of packet version numbers
     and JSOC Version numbers to use to lookup values for data series name */
  static JSOC_Version_Map_t *jmap;

  /* initialize variables */
  TIME *ptime= &pkt_time;
  ccsds=*ccsds_pkt;
  pvn = pkt_ver_num;


  /* check if record is set like for case of Lev 0 data series*/
  if(record)
  {
    //printk("write_hk_to_drms() called w/valid record\n"); //!!TEMP
    rec_alreadycreated_flag = 1;

#ifdef DEBUG_WRITE_HK_TO_DRMS
      printkerr("DEBUG:Message at %s, line %d: Preparing to write to DRMS "
                " for Level 0 Series.\n", __FILE__,__LINE__);
#endif
  }
  else 
  {
    //printk("write_hk_to_drms() called w/null record\n"); //!!TEMP
#ifdef DEBUG_WRITE_HK_TO_DRMS
      printkerr("DEBUG:Message at %s, line %d: Preparing to write to DRMS "
                " for Level 0 BY APID Series.\n", __FILE__,__LINE__);
#endif
    /*setting for writing to DRMS Level 0 by APID data series*/
    rec_alreadycreated_flag=0; 
  }

  /* check environments variables have values each time this function is called */
  if(haveNotChecked )
  { 
    if(((env_status = check_env_variable()) != 1))
    {
      /* if not a good value- returns value not equal to one */
      return (env_status);
    }
    else
    {
      /* set to checked only if passed setting all env variables */
      haveNotChecked=0;
    }
  }

  /* get directory and get file name suffix for file containing jsoc version numbers*/
  directory = (char *)getenv("HK_JSVNMAP_DIRECTORY");
  suffix_filename= (char *)getenv("HK_SUFFIX_JSVNMAP_FILENAME");
  if(suffix_filename == NULL)
  {
    printkerr("ERROR at %s, line %d: Set environment variable <HK_SUFFIX_JSVNMAP_FILENAME>\n",
               __FILE__,__LINE__);
    return  (ERROR_HK_ENVIRONMENT_VARS_NOT_SET);
  }
  if(directory == NULL)
  {
    printkerr("ERROR at %s, line %d: Set environment variable <HK_JSVNMAP_DIRECTORY>\n",
                __FILE__,__LINE__);
    return  (ERROR_HK_ENVIRONMENT_VARS_NOT_SET) ;
  }

  /* while loop threw CCSDS_Packet_t link list */
  while (ccsds)
  { /* outer-while thru CCSDS_Packet_t nodes */

    if (rec_alreadycreated_flag)
    {
      /* use record passed by main lev0 module to write to drms */
      rec= record;
      //Added by JA for IRIS. May be !!TEMP
      drms_setkey_double(rec, "DATE", CURRENT_SYSTEM_TIME);
      pvtmp = get_packet_version_number(ccsds->keywords);
      if(pvtmp) strcpy(pvn, pvtmp);
      //strcpy(pvn , get_packet_version_number(ccsds->keywords));   

      /* lookup data series name and set in DRMS*/
      keytype= DRMS_TYPE_STRING;
      strcpy(keyname, "ISPSNAME");

      /*allocate memory for string value */
      key_anyval.string_val = (char *)malloc(sizeof(char) * 100);

      /* get isp series name */
      strcpy(key_anyval.string_val, lookup_data_series_name(ccsds,&jmap,pvn));
      status = drms_setkey(rec, keyname, keytype, &key_anyval);

      /* free memory */
      free (key_anyval.string_val);
    }
    else
    {
      /* get packet version number */
      pvtmp = get_packet_version_number(ccsds->keywords);
      if(pvtmp) strcpy(pvn, pvtmp);
      //strcpy(pvn , get_packet_version_number(ccsds->keywords));   

      /*  lookup data series name */
      strcpy(query , lookup_data_series_name(ccsds, &jmap, pvn)); 
#ifdef DEBUG_WRITE_HK_TO_DRMS
      printkerr("DEBUG:Message at %s, line %d: Writing to DRMS data"
                " series <%s> for apid <%d>\n", __FILE__,__LINE__, query, ccsds->apid);
#endif

      /* check if record is within 12 hour range to remove processing of future dated packet data */
      /* get keyword structure */
      kw = ccsds->keywords;

      ck_rwtr_status = check_hk_record_within_time_range(kw);
      if (ck_rwtr_status == 1)
      {
        /* check if record exists based on timecode values */
        /* need ds name(above), apid and kw structure(below)  to check if record exists */
        kw = ccsds->keywords;
        //ck_rec_status = check_hk_record_exists(query, kw, ccsds->apid);
        ck_rec_status = 0;	//always do
#ifdef DEBUG_WRITE_HK_TO_DRMS
        printkerr("DEBUG:Message at %s, line %d: After check if record exists. "
                  "Status:<%d> where 0 means write record and 1 means skip "
                  "write of record.\n", __FILE__,__LINE__, ck_rec_status);
#endif
      }
      else
      {
        /*set to not write record since not within time range.*/
        ck_rec_status= 1; 
      }  

      /* if already exits or not within time range, skip writing record to DRMS data series */
      if (ck_rwtr_status == 0)
      {
        ;/* skip setting */
#ifdef DEBUG_WRITE_HK_TO_DRMS
        printkerr("DEBUG:Message at %s, line %d: Check if record within time range. "
                " Status:<%d>, so skip writing record\n",
                 __FILE__,__LINE__, ck_rwtr_status);
#endif
        /* set return status for ingest_lev0 to successful but skip writing since rec not within time range */
        /* of current date and time + 12 hours. The idea is to not write data the has made up future dates  */
        status=SUCCESS_HK_SKIP_WTD_REC_TIME_OOR;

        /* go to next node or next vcdu packet */
        ccsds= ccsds->next;
        continue;
      }
      else if (ck_rec_status == 1)
      {
        ;/* skip setting */
#ifdef DEBUG_WRITE_HK_TO_DRMS
        printkerr("DEBUG:Message at %s, line %d: Check if record exists. "
                " Status:<%d>, so skip writing record\n",
                 __FILE__,__LINE__, ck_rec_status);
#else
#endif
        /* set return status for ingest_lev0 to successful but skip writing since rec exits */
        status=SUCCESS_HK_SKIP_WTD_REC_EXISTS;

        /* go to next node or next vcdu packet */
        ccsds= ccsds->next;
        continue;
      }
      else
      {
#ifdef DEBUG_WRITE_HK_TO_DRMS
        printkerr("DEBUG:Message at %s, line %d: Check if record exists. "
                  "Status:<%d>, Check if record within range.Status:<%d> "
                  "so writing record\n", __FILE__,__LINE__, ck_rec_status,ck_rwtr_status);
#endif
      }

      /* create record in drms */
      RSISP = rs;		//close_image() must use last rs
      rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);
      RSISPTO = rs;		//if tlm file timeout must use this rs
      if (status < 0)
      {
        //printk("Set RSISP = 0 in write_hk_to_drms() fsn=%lu\n", fsnx); //!!TEMP
        rs = 0;
        //Probably getting one of these error codes.
        //DRMS_ERROR_UNKNOWNRECORD    (-10004)
        //DRMS_ERROR_RECORDREADONLY   (-10015)
        //DRMS_ERROR_BADRECORDCOUNT   (-10025)
        printkerr("ERROR at %s, line %d: Could not create record for data "
                  "series <%s> because of error code <%d>. Check code in "
                  "drms_statuscodes.h file. Skipping write to drms of "
                  "packet data to record in drms data series.\n",
                  __FILE__ , __LINE__ , query,status);
       ccsds= ccsds->next;
       continue;
      }
      else
      {
        //printk("Set RSISP = rs in write_hk_to_drms() fsn=%lu\n", fsnx); //!!TEMP
        rec = rs->records[0];
        //Added by JA for IRIS. May be !!TEMP
        drms_setkey_double(rec, "DATE", CURRENT_SYSTEM_TIME);
      }
    }

    /* get and set packet_time in record*/
    kw = ccsds->keywords;
    if (!get_pkt_time_from_timecodes(kw, ptime))
    {
      /*did not set PACKET_TIME */
      if (rec_alreadycreated_flag)
      {
        ;/* don't control record closing -skip*/
      }
      else
      {
       //status = drms_close_records(rs, DRMS_INSERT_RECORD); //!!new. no close
        //printk("drms_server_end_transaction() for ISP\n");
        //drms_server_end_transaction(drms_env, 0 , 0);
        //printk("drms_server_begin_transaction() for ISP\n");
        //drms_server_begin_transaction(drms_env); //start another cycle

      }
      /* assume need to have packet time set for both cases so do this */
      printkerr("ERROR at %s, line %d: Could not set PACKET_TIME keyword. "
                "Because the TIME CODE values were not found in keywords. "
                "Exiting write to drms of data series <%s> and closing record.\n",
                  __FILE__ , __LINE__ , query);
      return (ERROR_HK_FAILED_TO_FIND_TIMECODES);
    }
    else
    {
      /*if (fsn < FSN_TAI_FIXED) added based on input from Phil*/
      /* check if need to get actual fsn and check if less than  1,000,000*/
      //*ptime += 33.0; removed based on input from Phil 3-1-2008
                                                                                 
      /* set PACKET_TIME or ISPPKTIM keyword using data from TIMECODE keywords*/
      /* set drms type and long telemetry name  */
      keytype= DRMS_TYPE_TIME;
      strcpy(keyname, "PACKET_TIME");
      /* set packet time */
      key_anyval.time_val= *ptime;
      /* set record */
      if (rec_alreadycreated_flag)
      {
        status = drms_setkey(rec, "ISPPKTIM", keytype, &key_anyval);
      }
      else
      {
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
      }
    }


    /* set PACKET_VERSION_NUMBER(xxx.yyy) keyword */
    /* set drms type */
    keytype= DRMS_TYPE_STRING;
    if (rec_alreadycreated_flag)
    {
      /* set  defined keyword name in lev0 data series  */
      strcpy(keyname, "ISPPKTVN");
    }
    else
    {
      /* set  -long- telemetry name (PACKET_VERSION_NUMBER) keyword for lev0 by APID data series */
      strcpy(keyname, "PACKET_VERSION_NUMBER");
    }
    /*allocate memory for string value */
    key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
    /* set packet version number */
    strcpy(key_anyval.string_val, pvn);
    status = drms_setkey(rec, keyname, keytype, &key_anyval);
    /* free memory for string */
    free (key_anyval.string_val);
     

    /* set HK_SOURCE keyword */
    /* set drms type */
    keytype= DRMS_TYPE_STRING;
    if (rec_alreadycreated_flag)
    {
      ;//skip!
    }
    else
    {
      /* set  -long- telemetry name (HK_SOURCE) keyword for lev0 by APID data series */
      strcpy(keyname, "HK_SOURCE");
      /*allocate memory for string value */
      key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
      /* set string:<dayfile_name>.[<date>][<packet-apid>][<source>] */
      strcpy(key_anyval.string_val, get_hk_source(ccsds, ptime));
      status = drms_setkey(rec, keyname, keytype, &key_anyval);
      /* free memory for string */
      free (key_anyval.string_val);
    }
                                                                                 
    /* loop through keyword struct and load in record */
    kw = ccsds->keywords;
    while (kw)
    { /* inner while */
      if(kw->eng_type == KW_TYPE_UINT8)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }

        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_SHORT;

        /* set drms value by casting up for unsigned values */
        key_anyval.short_val = (int16_t)kw->eng_value.uint8_val;
    //!!TEMP
    printf("KEYWORD: %s %d\n", keyname, key_anyval.short_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_UINT16)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }

        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_INT;

        /* set drms value by casting up for unsigned values */
        key_anyval.int_val = (int32_t)kw->eng_value.uint16_val;
    //!!TEMP
    printf("KEYWORD: %s %d\n", keyname, key_anyval.int_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_UINT32)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }

        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_LONGLONG;

        /* set drms value by casting up for unsigned values */
        key_anyval.longlong_val = (int64_t)kw->eng_value.uint32_val;
    //!!TEMP
    printf("KEYWORD: %s %ld\n", keyname, key_anyval.longlong_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_INT8)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
	}

         /* set drms type but promote up drms short since byte is a number value that can contain non-valid chars */
         keytype= DRMS_TYPE_SHORT;
         key_anyval.short_val = (int16_t)kw->eng_value.int8_val;
    //!!TEMP
    printf("KEYWORD: %s %d\n", keyname, key_anyval.short_val);
         status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_DOUBLE)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }

        /* set drms type but promote up for unsigned values */
        keytype= DRMS_TYPE_DOUBLE;
        key_anyval.double_val = kw->eng_value.double_val;
    //!!TEMP
    printf("KEYWORD: %s %g\n", keyname, key_anyval.double_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_INT16)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }

        /* set drms type  */
        keytype= DRMS_TYPE_SHORT;
        key_anyval.short_val = kw->eng_value.int16_val;
    //!!TEMP
    printf("KEYWORD: %s %d\n", keyname, key_anyval.short_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_INT32)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }
        /* set drms type  */
        keytype= DRMS_TYPE_INT;
        key_anyval.int_val= kw->eng_value.int32_val;
    //!!TEMP
    printf("KEYWORD: %s %d\n", keyname, key_anyval.int_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);

      }
      else if(kw->eng_type == KW_TYPE_STRING)
      {
        if(rec_alreadycreated_flag) 
        {
          strcpy(keyname, kw->fitsname);
        }
        else   
        {
          strcpy(keyname, kw->name);
        }

        /* set drms type  */
        keytype= DRMS_TYPE_STRING;

        /*allocate memory for string */
        key_anyval.string_val = (char *)malloc(sizeof(char) * 100);
        strcpy(key_anyval.string_val, kw->eng_value.string_val);
    //!!TEMP
    printf("KEYWORD: %s %s\n", keyname, key_anyval.string_val);
        status = drms_setkey(rec, keyname, keytype, &key_anyval);
        free(key_anyval.string_val);

      }
      else
      {
        printkerr("Warning at %s, line %d: Found Unknown KW_TYPE for this keyword: "
                  "<%d>. Check if new KW_TYPE has been created in HKPDF files "
                  "Skipping recording keyword in structure. Keyword is <%s>.",
                 __FILE__,__LINE__, (int)kw->eng_type, kw->name);
      }
//!!TEMP
//if(!strcmp(keyname, "ISQFSN")) {
//  printk("In write_hk_to_drms() ISQFSN = %d\n", key_anyval.int_val);
//}


#ifdef DEBUG_WRITE_HK_TO_DRMS_MORE
      /* check got good return status */
      if (status)
      {
        /* compile in by adding -DDEBUG_WRITE_HK_TO_DRMS in makefile */
        if (rec_alreadycreated_flag)
        {
          printkerr("DEBUG Messsage at %s, line %d: Warning-cannot setkey "
                    "in drms record when writing to <Level 0 series>. "
                    "Return status from drms_setkey was <%d>. "
                    "For keyword <%s>\n", __FILE__,__LINE__, status, keyname);
        }
        else
        {
          printkerr("DEBUG Messsage at %s, line %d: Warning-cannot setkey "
                    "in drms record when writing to <Level 0 By APID series>. "
                    "Return status from drms_setkey was <%d>. "
                    "For keyword <%s>\n", __FILE__,__LINE__, status, keyname);
        }
      }
#else
#endif
        
      /* Next Keyword Node */
      kw=kw->next;

    } /* inner while */

    /* close record */
    if (rec_alreadycreated_flag)
    {
       /* don't control record closing -skip*/
    }
    else
    {
      //status = drms_close_records(rs, DRMS_INSERT_RECORD); //!!new no close
        //printk("drms_server_end_transaction() for ISP\n");
        //drms_server_end_transaction(drms_env, 0 , 0);
        //printk("drms_server_begin_transaction() for ISP\n");
        //drms_server_begin_transaction(drms_env); //start another cycle
      //if (status)
      //{
      //  printkerr("ERROR at %s, line %d: Cannot close drms record.\n",
      //             __FILE__ , __LINE__ );
      //}
    }
    /* don't expect more than one node to rec_alreadycreated_flag case,
       so don;t need to go to next CCSDS_Packet_t node. If do get
       more than one node a DRMS_Record_Set_t need to be used as argument  */
    ccsds= ccsds->next;
#ifdef DEBUG_WRITE_HK_TO_DRMS
    printkerr("DEBUG:Message at %s, line %d: Complete CCSDS packet. "
              "Do next CCSDS packet.\n", __FILE__,__LINE__);
#else
#endif
    if (rec_alreadycreated_flag)
    {
       /* check if more than one CCSDS_Packet_t nodes are available */
       if (ccsds)
       {
          printkerr("ERROR at %s, line %d: Getting more than one node in CCSDS_Packet_t "
                    "link list for Lev0.1 case where one record is passed. "
                    "To handle not to pass DRMS_RECORD_SET_T as argument.\n",
                   __FILE__ , __LINE__ );
          break;
       }
    }
  }/* outer-while*/
  return status;
}



/******************************************************************************
 * Name:        GET PACKET VERSION NUMBER                                     *
 * Function :   char *get_packet_version_number( HK_Keyword_t *kw)            *
 * Description: Gets packet version number in HK_Keyword_t link list          *
 *              structure. Returns "string" in format 001.120, or 001.001,etc.*             
 *****************************************************************************/
char *get_packet_version_number( HK_Keyword_t *kw)
{
  /* declarations */
  char vn[HK_LEV0_MAX_PVN_STR];
  char *ptr_vn;

  /* initialized variables */
  ptr_vn = vn;

  while (kw)
  {
     //if( strstr(kw->name, "VER_NUM") || strstr(kw->name, "VER_TEMPERATURES")) //Fix 11-30-2009
     if(strstr(kw->name, "APID56_APID_VERSION"))  //This is the IRIS keyword name
     {
       sprintf(ptr_vn, "%03d.%03d", kw->eng_value.uint16_val >> 8  & 0x00FF, 
               kw->eng_value.uint16_val & 0x00FF );
       strcat(ptr_vn,"\0"); 
       return (ptr_vn);
     }
     kw= kw->next;
  }
  return NULL;
}


/*****************************************************************************
 * Name:        Get Data Packet Name                                         *
 * Function :   char *get_data_packet_name(int apid )                        *
 * Description: Get data packet name based on apid.Return a pointer to string* 
 *****************************************************************************/
char *get_data_packet_name(int apid ) 
{
  /* declarations */
  char *dpn;
  char data_packet_name[HK_LEV0_MAX_PACKET_NAME];

  /* initialized variables */
  dpn = data_packet_name;

  /* determine value of packet name to use for data series */
  /* set packet name based on the APID. */
  switch(apid)
  {  
  /*****************************************/
  /***  Case of Finding HMI ISP APIDS    ***/
  /*****************************************/
   case HK_LR_HMI_ISP:
   case HK_HSB_HMI_ISP_1:
   case HK_HSB_HMI_ISP_2:

        /* ISP apids */      
        strcpy(dpn, HK_HMI_PKT_NAME_ISP);
        break;

   case HK_HSB_IRIS_ISP:
        //IRIS ISP
        strcpy(dpn, HK_HSB_IRIS_PKT_NAME_ISP);
        break;

  /*****************************************/
  /***  Case of Finding AIA ISP APIDS    ***/
  /*****************************************/
   case HK_LR_AIA_ISP:
   case HK_HSB_AIA_ISP_1:
   case HK_HSB_AIA_ISP_2:

        strcpy(dpn, HK_AIA_PKT_NAME_ISP);
        break;

  /*****************************************/
  /***Case of Finding HMI SEQUENCER APIDS***/
  /*****************************************/
   case HK_LR_HMI_SEQ:
   case HK_HSB_HMI_SEQ_1:
   case HK_HSB_HMI_SEQ_2:

        strcpy(dpn, HK_HMI_PKT_NAME_SEQ);
        break;

  /*****************************************/
  /***Case of Finding AIA SEQUENCER APIDS***/
  /*****************************************/
   case HK_LR_AIA_SEQ:
   case HK_HSB_AIA_SEQ_1:
   case HK_HSB_AIA_SEQ_2:

        strcpy(dpn, HK_AIA_PKT_NAME_SEQ);
        break;

  /*****************************************/
  /***Case of Finding HMI OBT APIDS      ***/
  /*****************************************/
   case HK_LR_HMI_OBT:
   case HK_HSB_HMI_OBT_1:
   case HK_HSB_HMI_OBT_2:

        strcpy(dpn, HK_HMI_PKT_NAME_OBT);
        break;

  /*****************************************/
  /***Case of Finding AIA OBT APIDS      ***/
  /*****************************************/
   case HK_LR_AIA_OBT:
   case HK_HSB_AIA_OBT_1:
   case HK_HSB_AIA_OBT_2:

        strcpy(dpn, HK_AIA_PKT_NAME_OBT);
        break;

  /******************************************************/
  /***Case of Finding SDO ANCILLORY SCIENCE DATA APIDS***/
  /******************************************************/
   case HK_LR_SDO_ASD:

        strcpy(dpn, HK_SDO_PKT_NAME_ASD);
        break;

  /*************************************************/
  /*** All Other cases of undefined packet names ***/
  /*************************************************/
   default:
      /* use apid to name */
      sprintf(dpn,"%04d", apid);
      break;
  }/* end switch */

  return( dpn);
}



/*****************************************************************************
 * Name:        LOOKUP DATA SERIES NAME                                      *
 * Function :   char *lookup_data_series_name(CCSDS_Packet_t *ccsds_ptr,     *
 *                                            JSOC_Version_Map_t **jmap_ptr, *
 *                                             char packet_version_number)   *
 * Description: lookup the data series name based on apid and packet         * 
 *              packet version number. Returns string value                  *
 *****************************************************************************/
char *lookup_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t **jmap_ptr, char *pkt_ver_num) 
{
  //!!NOTE !!TMP?? For IRIS use a fixed ISP data set name defined in
  //SOURCE_ENV_FOR_HK_DECODE file
  sprintf(IRIS_ISP_DS, "%s", getenv("HK_ISP_IRIS_DSNAME"));
  return IRIS_ISP_DS;

   /* declarations */
   JSOC_Version_Map_t *jm, *tmp_jm, *last_jm, *prev_jm;
   char data_series_name[HK_LEV0_MAX_DSNAME_STR];
   char *dsn;
   int  apid;
   static short curr_mf_flag=0;
   short new_mf_flag;

   /* initialized variables */
   dsn = data_series_name;
   jm = *jmap_ptr;

   /* get apid */
   apid = ccsds_ptr->apid;

   /* check if using old map files or new merged map files */
   new_mf_flag = check_curr_map_file(pkt_ver_num, &curr_mf_flag);
#ifdef DEBUG_WRITE_HK_TO_DRMS
   printkerr("DEBUG:Message:newmfflag:%d curr:%d \n",new_mf_flag,curr_mf_flag);
#else
#endif

   /*check if structure has node with data series name */
   if (!jm || new_mf_flag)
   {
      /* if empty structure then create top node */
      jm = (JSOC_Version_Map_t *)malloc(sizeof(JSOC_Version_Map_t)); 
      if (!jm)
      {
        printkerr("ERROR at %s, line %d: Cannot malloc space!\n", __FILE__ , __LINE__ );
      }
      *jmap_ptr=jm;
      jm->apid=(short)apid;
      jm->next=NULL;

      /* load Map_Data nodes and create data series names 
         for each packet version number*/
      load_map_data(apid, jm, pkt_ver_num);
  
      /* find data series name in Map_Data nodes based on 
         packet version number and return*/
      (void)find_data_series_name(ccsds_ptr, jm, dsn);

   }
   else
   {
     /*find data series name based on packet version number and apid */
     (void)find_data_series_name(ccsds_ptr, jm, dsn);

     /* if don't find then read map file and load into 
        JSOC Version Number Map structure */
     if (!strcmp(dsn, "\0"))
     { 
       /* if did not find data series name in jm structure then
          check if apid value is there in jm, if there, then need to load
          probably  new PVN-TO-JSVN file with new Map_Data node with 
          needed packet version number line*/
       /*TEST This CASE **/
       if (check_for_apid(ccsds_ptr->apid, jm))
       {
           printkerr("Warning at %s, line %d: Could not find find "
              "packet version number and JSOC version number line for " 
              "APID: <%d> in Data Map structure. Need to reread files.\n", 
              __FILE__, __LINE__, ccsds_ptr->apid);

           /* NOTE NEED TO TEST THIS CASE !!! */
           /* found in JSOC_Version_Map node apid needed */
           /* free entire jm structure */
           (void)free_jsvn_map(jm);
           /* recreate jm structure by loading on demand PVN-TO-JSVN file data*/
           /* create top node */
           jm = (JSOC_Version_Map_t *)malloc(sizeof(JSOC_Version_Map_t));
           if (!jm)
           {
             printkerr("ERROR at %s, line %d: Cannot malloc space!\n", __FILE__ , __LINE__ );
           } 
           *jmap_ptr=jm;
           jm->apid=(short)apid;
           jm->next=NULL;
                                                                                                
           /* load Map_Data nodes and create data series names for each 
              packet version numberi and return data series name*/
           load_map_data(apid, jm, pkt_ver_num);
                                                                                                
           /* find data series name in Map_Data nodes based on packet version number*/
           (void)find_data_series_name(ccsds_ptr, jm, dsn);
       } 
       else 
       { 
         /* if cannot find node for apid then create one */ 
         /* go to end of JSOC Version Map nodes and add node */
         /* first create JSOC Version Map node */
         tmp_jm = (JSOC_Version_Map_t *)malloc(sizeof(JSOC_Version_Map_t)); 
         if (!tmp_jm)
         {
           printkerr("ERROR at %s, line %d: Cannot malloc space!\n", __FILE__ , __LINE__ );
         } 
         tmp_jm->apid= (short)apid;
         tmp_jm->next= NULL;

         /* link in list of map data to JSOC Version Map node */
         load_map_data(apid, tmp_jm, pkt_ver_num);     

         /*find and connect to JSOC Version Map's last node */
         for( last_jm= *jmap_ptr; last_jm; prev_jm=last_jm, last_jm=last_jm->next);
         prev_jm->next =tmp_jm;

         /* find data series name based on apid and packet version number */
         (void)find_data_series_name(ccsds_ptr, jm, dsn);
       } /* else if did not find JSOC_Version_Map node with needed apid */

     } /* end if did not find dsn value */
   } /*else if have data in jmap structure case */
   return dsn;
}



/*****************************************************************************
 * Name     :   LOAD MAP DATA                                                *
 * Function :   void load_map_data(int apid, JSOC_Version_Map_t  *jm,        *
 *                                 char *packet_version_number)              *
 * Description: Loads PVN-TO-JSVN-<apid> files data in structures. This      *
 *              structure will be used to lookup data series name.           *
 *****************************************************************************/
void load_map_data(int apid, JSOC_Version_Map_t  *jm, char *pvn)
{
  /*declarations */
  FILE *file_ptr;
  Map_Data_t  *tmp_dm, *prev_dm;
  char directory_filename[HK_LEV0_MAX_FILE_NAME];
  char dirname[HK_LEV0_MAX_FILE_NAME];
  char fn[HK_LEV0_MAX_FILE_NAME];
  char *filename;
  char *directory;
  char sn[HK_LEV0_MAX_FILE_NAME];
  char *didn;
  char *suffix_filename;
  char line[HK_LEV0_MAXLINE_IN_FILE];
  char *pn;
  int  wn, dn;
  int curr_pvn_wn,curr_pvn_dn;

  /* initialized variables */
  filename=fn;
  directory=dirname;

  /* initialized variables */
  suffix_filename=sn;


  /* get project name */
 if(apid <= HK_HSB_HIGHEST_HMI_APID && apid >= HK_HSB_LOWEST_HMI_APID  || (apid <= HK_LR_HIGHEST_HMI_APID  && apid >= HK_LR_LOWEST_HMI_APID))
  {
    pn = getenv("HK_LEV0_BY_APID_PROJECT_NAME_HMI");
  }
  else if((apid <= HK_HSB_HIGHEST_AIA_APID && apid >= HK_HSB_LOWEST_AIA_APID) || (apid <= HK_LR_HIGHEST_AIA_APID && apid >= HK_LR_LOWEST_AIA_APID))
  {
    pn = getenv("HK_LEV0_BY_APID_PROJECT_NAME_AIA");
  }
  else if((apid <= HK_LR_HIGHEST_SDO_APID) && (apid >= HK_LR_LOWEST_SDO_APID))
  {
    pn = getenv("HK_LEV0_BY_APID_PROJECT_NAME_SDO");
  }
  else
  {
    printkerr("Warning at %s, line %d: APID is not in range of ",
              "HMI, AIA, or SDO. APID: <%d>\n",
               __FILE__, __LINE__, apid);
  }

  /* get data type name */
  didn = getenv("HK_LEV0_BY_APID_DATA_ID_NAME");

  /* get directory and file name  with jsoc version numbers*/
  directory = (char *)getenv("HK_JSVNMAP_DIRECTORY");
  suffix_filename= (char *)getenv("HK_SUFFIX_JSVNMAP_FILENAME");

  /* make filename for looking up of jsoc version number */
  /* get number value of packet version number for comparison later */
  sscanf(pvn,"%3d.%3d",&curr_pvn_wn, &curr_pvn_dn);

  /* check if merge data series packet version number */
  if(curr_pvn_wn  >=  HK_LEV0_START_MERGED_PVNW && curr_pvn_dn >= HK_LEV0_START_MERGED_PVND)
  {
    strcpy(filename,get_lookup_filename(apid));
    strcat(filename,suffix_filename);
#ifdef DEBUG_WRITE_HK_TO_DRMS
    printkerr("DEBUG:Message at %s, line %d: Triggered merge case. "
              "Getting JSOC map file with new name: <%s>\n",
               __FILE__, __LINE__, filename);
#else
#endif
  }
  else
  {
    sprintf(filename,"%4.4d%s",apid,suffix_filename);
#ifdef DEBUG_WRITE_HK_TO_DRMS
    printkerr("DEBUG:Message at %s, line %d: Triggered non-merge case. "
              "Getting JSOC map file with old name: <%s>\n",
               __FILE__, __LINE__, filename);
#else
#endif
  }

  /* put together filename and directory path for lookup jsoc version file */
  sprintf(directory_filename, "%s/%s", directory, filename);

  /* get apid and set to JSOC Version Map node structure */
  jm->apid=(short)apid;

  /* open JSOC Version Map file and read */
  file_ptr=fopen(directory_filename,"r");
  if (!file_ptr)
  {
    printkerr("Error at %s, line %d: Check if directory and filename are",
              "correct. Could not open filename: <%s>\n",
               __FILE__, __LINE__, directory_filename);
    return;
  }

  /* get lines in JSVN-TO-PVN files */
  tmp_dm=NULL;
  while( fgets(line, HK_LEV0_MAXLINE_IN_FILE, file_ptr) != NULL )
  {
    if(line[0] == '#')
    {
      continue; /* skip comments */
    }
    else
    { /* got good data to load in structure */
      /* if first mdata node */
      /* for each line load map data in Map_Data_t structure */
      if (!tmp_dm)
      {
        tmp_dm = (Map_Data_t*)malloc(sizeof(Map_Data_t)); 
        if (!tmp_dm)
        {
          printkerr("ERROR at %s, line %d: Cannot malloc space!\n", __FILE__ , __LINE__ );
        } 
        tmp_dm->next= NULL;
        jm->mdata = tmp_dm;
        prev_dm =tmp_dm;
      }
      else
      {
        tmp_dm = (Map_Data_t*)malloc(sizeof(Map_Data_t)); 
        if (!tmp_dm)
        {
          printkerr("ERROR at %s, line %d: Cannot malloc space!\n", __FILE__ , __LINE__ );
        } 
        tmp_dm->next= NULL;
        prev_dm->next =tmp_dm;
        prev_dm= tmp_dm;
      }

      /* set jsoc version number field */
      sscanf( line,
         "%s | %d.%d | %*s | %*s", tmp_dm->jvn, &wn,&dn);

      /* set packet version number field */
      sprintf(tmp_dm->pvn,"%03d.%03d", wn,dn);

      /* if jsvn is 0 make dsn null else set data series name in array*/
      if (!strcmp("0000",  tmp_dm->jvn))
      {
        /* if jsv is 0 then set data series name to NULL */
        strcpy(tmp_dm->dsn, "\0");
      }
      else
      {

        /* check if pvn is for merge series names */
        if(curr_pvn_wn  >=  HK_LEV0_START_MERGED_PVNW && curr_pvn_dn >= HK_LEV0_START_MERGED_PVND)
        {
          /* if pvn => 1.197 then used merged name  */
          sprintf(tmp_dm->dsn, "%s.%s_%s_%s", pn, didn, get_data_packet_name(apid), tmp_dm->jvn);
        }
        else
        {
          /* else use non-merged name */
          sprintf(tmp_dm->dsn, "%s.%s_%4.4d_%s", pn, didn, apid, tmp_dm->jvn);
        }
      }
    } /*end else -got good data */
  } /* while loop thru lines in PVN-TO-JSVN files */

  fclose( file_ptr); 
  return ;
}



/*************************************************************************
 * FIND DATA SERIES NAME                                                 *
 * FUNCTION: char *find_data_series_name(CCSDS_Packet_t *ccsds_ptr,      *
 *                                       JSOC_Version_Map_t  *jm)        *
 * DESCRIPTION: find data series name in JSOC Version Map structure      *
 *************************************************************************/
int find_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t  *jm, char *dsn)
{
  /* declarations */
  JSOC_Version_Map_t  *tmp_jm;
  Map_Data_t *dm;
  char data_series_name[HK_LEV0_MAX_DSNAME_STR];
  char  pkt_ver_num[HK_LEV0_MAX_PVN_STR];
  char *pvn;

  /* initialized varibles */
  pvn = pkt_ver_num;
  tmp_jm= jm;

  /* set to NULL if cannot find in structure return NULL */
  strcpy (data_series_name,"\0");

  // get packet version number from ccsds 
  strcpy(pvn , get_packet_version_number(ccsds_ptr->keywords));   

  // lookup apid in jm structure
  // lookup packet version number in jm structure's Map_Data_t link list
  // Then return data series name if found and return null if did not find.
  if (!jm) 
  {
    return 0;
  }
  else
  {
    while (tmp_jm)
    {
      /* loop thru JSOC_Version_Data node to find if node has apid looking for */
      if(tmp_jm->apid == ccsds_ptr->apid)
      {
        /* found apid */
        /* Loop thru Map_Data_t nodes to get data series name based on 
           packet version number */
        dm= tmp_jm->mdata;
        while (dm)
        {
           /* find packet version number */
           if ( !strcmp(dm->pvn, pvn))
           {
              /* found pvn-now copy in data series name to return */
              strcpy (data_series_name, dm->dsn);
              break;
           }
            dm = dm->next;
        }
        break; 
      }
      else
      {
        /* go to next node looking for apid */
        tmp_jm= tmp_jm->next;
      }
    } /* while loop thru jm nodes */
  }
  strcpy(dsn, data_series_name);
#ifdef DEBUG_WRITE_HK_TO_DRMS
  printkerr("DEBUG:Message:FOUND Lookup series name:%s\n",data_series_name);
#else
#endif
  return 1;
}


                                                                                               
/*************************************************************************
 * GET PACKET TIME FROM TIME CODES                                       *
 * FUNCTION: int get_pkt_time_from_timecodes(HK_Keyword_t *,char *,TIME) *
 * DESCRIPTION: Gets Packet time based on Time Code values.              *
 *************************************************************************/
int get_pkt_time_from_timecodes(HK_Keyword_t *hk,  TIME *ptime)
{
  /* declarations */
  HK_Keyword_t *t_hk;
  int sec;
  int subsec;
  int SEC_FOUND_FLAG=0;
  int SUBSEC_FOUND_FLAG=0;

  /* initialize variables */
  t_hk= hk;
                                                                                                                    
  /*loop until get TIMECODE Seconds and Subseconds values */
  while (t_hk && ( !SEC_FOUND_FLAG || !SUBSEC_FOUND_FLAG ))
  {
    if(strstr(t_hk->name,"TIMECODE_SEC"))	//name for IRIS
    {
      /*set found flag*/
      SEC_FOUND_FLAG=1;
      /* create time string based on some hard coded parameters for now */
      sec = t_hk->eng_value.uint32_val;
    }
    if(strstr(t_hk->name,"TIMECODE_SSEC"))  //name for IRIS
    {
      /*set found flag*/
      SUBSEC_FOUND_FLAG=1;
      /* create time string based on some hard coded parameters for now */
      subsec = t_hk->eng_value.uint32_val;
    }
    t_hk= t_hk->next;
  }
                                                                                                
  if (!SEC_FOUND_FLAG)
  {
    printkerr("ERROR at %s, line %d: Did not find TIMECODE_SECONDS value for "
              "calculating the PACKET_TIME keyword and index. Returning error "
              "status.\n",__FILE__,__LINE__);
    return 0;
  }
  else
  {
    if (!SUBSEC_FOUND_FLAG)
    {
      printkerr("ERROR at %s, line %d: Did not find TIMECODE_SUBSECS value for "
                "calculating the PACKET_TIME keyword and index. Returning error "
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
 * NAME: Check for Apid                                                  *
 * FUNCTION: int check_for_apid(int apid, JSOC_Version_Map_t  *jm )      *
 * DESCRIPTION: Check for apid  in JSOC_Version_Map_t structure          *
 *************************************************************************/
int check_for_apid(int apid, JSOC_Version_Map_t  *top_jm)
{
  JSOC_Version_Map_t *tmp_jm;
  tmp_jm= top_jm;

  while( tmp_jm)
  {
     if(tmp_jm->apid == apid)
     {
       return 1; //FOUND
     }
     tmp_jm= tmp_jm->next;
  }
  return 0;//DID NOT FIND
}

/*************************************************************************
 * NAME: free JSOC VERSION NUMBER MAP                                    *
 * FUNCTION: int free_jsvn_map( JSOC_Version_Map_t  *jm )                *
 * DESCRIPTION: Free JSOC_Version_Map_t nodes and Map_Data_t nodes       *
 *************************************************************************/
void free_jsvn_map( JSOC_Version_Map_t  *top_jm)
{
  /* declarations */
  JSOC_Version_Map_t  *tmp_jm, *prev_jm;
  Map_Data_t *prev_md, *tmp_md;

  /* initialized variables */
  tmp_jm= top_jm;

  while( tmp_jm)
  {
     tmp_md= tmp_jm->mdata;
     while(tmp_md)
     {
       prev_md = tmp_md;
       tmp_md= tmp_md->next;
       free(prev_md);
     }
     prev_jm = tmp_jm;
     tmp_jm= tmp_jm->next;
     free(prev_jm);
  }
  return;
}

/*********************/
/* SDO_to_DRMS_time  */
/*********************/
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



/************************************************************/
/* NAME:Get HK_SOURCE Value to Set in Lev0 by APID Series   */
/* FUNCTION:   get_hk_source()                              */
/* DESCRIPTION:Get value to set HK_SOURCE keyword in Lev0   */ 
/*             by APID Series by using PACKET_TIME to set   */
/*             Date indx value, Apid to set apid index and  */
/*             hsb to source index value. The project name  */
/*             to use is based on choosing three environment*/
/*             HK_LEV0_DF_MAX_PROJECT_NAME_xxx where xxx is */
/*             either SDO,HMI or AIA.                       */
/************************************************************/
char *get_hk_source(CCSDS_Packet_t *ccsds_ptr, TIME *pt)
{
  /* declarations */
  char date[HK_LEV0_MAX_HKS_DATE];
  char df_datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char df_project_name[HK_LEV0_MAX_PROJECT_NAME];
  char *df_dtname, *df_pjname;
  char hk_source[HK_LEV0_MAX_HKS];
  char *hks;
  int  apid;
  int  len_hk_source;
  int  pkt_yr, pkt_month, pkt_day;

  /* initial pointers */
  df_dtname=df_datatype_name;
  df_pjname=df_project_name;
  hks = hk_source;

  /* get apid */
  apid = ccsds_ptr->apid;

  /* get dayfile project name */
  if(apid <= HK_HSB_HIGHEST_HMI_APID && apid >= HK_HSB_LOWEST_HMI_APID  || (apid <= HK_LR_HIGHEST_HMI_APID  && apid >= HK_LR_LOWEST_HMI_APID))
  {
    df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_HMI");
  }
  else if((apid <= HK_HSB_HIGHEST_AIA_APID && apid >= HK_HSB_LOWEST_AIA_APID) || (apid <= HK_LR_HIGHEST_AIA_APID && apid >= HK_LR_LOWEST_AIA_APID))
  {
    df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_AIA");
  }
  else if((apid <= HK_LR_HIGHEST_SDO_APID) && (apid >= HK_LR_LOWEST_SDO_APID))
  {
    df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_SDO");
  }
  else
  {
    printkerr("Warning at %s, line %d: APID is not in range of "
              "HMI, AIA, or SDO. APID: <%d>\n",
               __FILE__, __LINE__, apid);
  }
    df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_IRIS"); //override for iris

  /* get dayfile data type name */
  df_dtname = getenv("HK_LEV0_DF_DATA_ID_NAME");

  /* get date from current packet to know date to use for index to dayfile series */
  pkt_yr=get_yr_from_pkttime(*pt);
  pkt_month=get_month_from_pkttime(*pt);
  pkt_day=get_day_from_pkttime(*pt);

  /* Put dayfile time in this format: 2007.11.08_00:00:00.00_TAI */
  sprintf(date,"%d.%-02.2d.%-02.2d_00:00:00.00_TAI",pkt_yr, pkt_month, pkt_day);

  /* create string for hk source using format:hmi.hk_dayfile[2007.11.08_00:00:00.00_TAI][475][hsb] */
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
    sprintf(hks,"%s.%s[%s][%d][%s]\0", df_pjname, df_dtname, date, apid, "hsb");
  }

  return(hks);
}



/************************************************************/
/* NAME:Check Environment Variable                          */
/* FUNCTION:   check_env_variable()                         */
/* DESCRIPTION: Checks all environment variables are set    */ 
/*              SOURCE_ENV_HK_DECODE file. Return -23 if    */
/*              variable is not set.                        */
/************************************************************/
int check_env_variable()
{
  /* declarations */
  char datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char df_project_name[HK_LEV0_MAX_PROJECT_NAME];
  char df_datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char project_name[HK_LEV0_MAX_PROJECT_NAME];
  char jvn_lu_hmi_prefix[HK_LEV0_MAX_LU_PREFIX_FILENAME];
  char jvn_lu_aia_prefix[HK_LEV0_MAX_LU_PREFIX_FILENAME];
  char jvn_lu_sdo_prefix[HK_LEV0_MAX_LU_PREFIX_FILENAME];
  char *dtname, *pjname, *df_dtname, *df_pjname;
  char *hmiprefix,*sdoprefix,*aiaprefix;

  /* initialize pointers */
  pjname=project_name;
  df_dtname=df_datatype_name;
  df_pjname=df_project_name;
  dtname=datatype_name;
  hmiprefix=jvn_lu_hmi_prefix;
  aiaprefix=jvn_lu_aia_prefix;
  sdoprefix=jvn_lu_sdo_prefix;

  /* check 3 environment variables are set for HK By APID series name */
  /* get project name for HK by APID data series */
  pjname = getenv("HK_LEV0_BY_APID_PROJECT_NAME_HMI");
  if(pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_LEV0_BY_APID_PROJECT_NAME_HMI>. Set the env variable "
              "HK_LEV0_BY_APID_PROJECT_NAME_HMI to data series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing data series name.\n",
              __FILE__,__LINE__);
    return  ERROR_HK_ENVIRONMENT_VARS_NOT_SET ;
  }
  pjname = getenv("HK_LEV0_BY_APID_PROJECT_NAME_AIA");
  if(pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_LEV0_BY_APID_PROJECT_NAME_AIA>. Set the env variable "
              "HK_LEV0_BY_APID_PROJECT_NAME_HMI to data series project name"
              "(aia, aia_ground, su_carl,etc.) for a existing data series name.\n",
              __FILE__,__LINE__);
    return  ERROR_HK_ENVIRONMENT_VARS_NOT_SET ;
  }
  pjname = getenv("HK_LEV0_BY_APID_PROJECT_NAME_SDO");
  if(pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_LEV0_BY_APID_PROJECT_NAME_SDO>. Set the env variable "
              "HK_LEV0_BY_APID_PROJECT_NAME_SDO to data series project name"
              "(sdo, sdo_ground, su_carl,etc.) for a existing data series name.\n",
              __FILE__,__LINE__);
    return  ERROR_HK_ENVIRONMENT_VARS_NOT_SET ;
  }

  /* get and check data type name for HK BY APID data series */
  dtname = getenv("HK_LEV0_BY_APID_DATA_ID_NAME");
  if(dtname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get data type name environment "
              "variable:<HK_LEV0_BY_APID_DATA_ID_NAME>. Set the env variable "
              "HK_LEV0_BY_APID_DATA_ID_NAME to data series data type name"
              "(lev0,lev0test, etc.) for a existing data series name\n",
              __FILE__ , __LINE__ );
     return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET) ;
  }

  /* get and check 3 environment variables are setting dayfile reference name   */

  /* get dayfile project name for setting HK_SOURCE value in HK by APID data series */
  df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_HMI");
  if(df_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_LEV0_DF_PROJECT_NAME_HMI>. Set the env variable "
              "HK_LEV0_DF_PROJECT_NAME_HMI to dayfile series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing dayfile series name.\n",
              __FILE__,__LINE__);
    return  ERROR_HK_ENVIRONMENT_VARS_NOT_SET ;
  }
  df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_AIA");
  if(df_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_LEV0_DF_PROJECT_NAME_AIA>. Set the env variable "
              "HK_LEV0_DF_PROJECT_NAME_AIA to dayfile series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing dayfile series name.\n",
              __FILE__,__LINE__);
    return  ERROR_HK_ENVIRONMENT_VARS_NOT_SET ;
  }
  df_pjname = getenv("HK_LEV0_DF_PROJECT_NAME_SDO");
  if(df_pjname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get project name environment "
              "variable:<HK_LEV0_DF_PROJECT_NAME_SDO>. Set the env variable "
              "HK_LEV0_DF_PROJECT_NAME_SDO to dayfile series project name"
              "(hmi, hmi_ground, su_carl,etc.) for a existing dayfile series name.\n",
              __FILE__,__LINE__);
    return  ERROR_HK_ENVIRONMENT_VARS_NOT_SET ;
  }
  
  /* get data type name for HK BY APID data series */
  df_dtname = getenv("HK_LEV0_DF_DATA_ID_NAME");
  if(df_dtname == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get data type name environment "
              "variable:<HK_LEV0_DF_DATA_ID_NAME>. Set the env variable "
              "HK_LEV0_DF_DATA_ID_NAME to data series data type name"
              "(lev0,lev0test, etc.) for a existing data series name\n",
              __FILE__ , __LINE__ );
     return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET) ;
  }

  /* get and check each of the three prefix values for JSOC Version Number lookup file */
  hmiprefix=getenv("HK_PREFIX_JSVNMAP_FILENAME_HMI");
  if(hmiprefix == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get lookup file environment "
              "variable:<HK_PREFIX_JSVNMAP_FILENAME_HMI>. Set the env variable "
              "HK_PREFIX_JSVNMAP_FILENAME_HMI to correct prefix name"
              "(HMI-ISP or HMI-SEQ, etc.). Current value is null. \n", 
              __FILE__ , __LINE__ );
    return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET) ;
  }

  aiaprefix=getenv("HK_PREFIX_JSVNMAP_FILENAME_AIA");
  if(aiaprefix == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get lookup file environment "
              "variable:<HK_PREFIX_JSVNMAP_FILENAME_AIA>. Set the env variable "
              "HK_PREFIX_JSVNMAP_FILENAME_AIA to correct prefix name"
              "(AIA-ISP or AIA-SEQ, etc.). Current value is null.\n", 
              __FILE__ , __LINE__ );
    return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET) ;
  }

  sdoprefix=getenv("HK_PREFIX_JSVNMAP_FILENAME_SDO");
  if(sdoprefix == NULL) 
  {
    printkerr("ERROR at %s, line %d: Could not get lookup file environment "
              "variable:<HK_PREFIX_JSVNMAP_FILENAME_SDO>. Set the env variable "
              "HK_PREFIX_JSVNMAP_FILENAME_SDO to correct prefix name"
              "(SDO-ADP). Current value is null.\n", 
              __FILE__ , __LINE__ );
    return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET) ;
  }

  /* all environment variable are at least set, so return 1 */
  return (1);
}



/*****************************************************************************
 * Name:        Get LOOKUP FILENAME    - get JSVN lookup filename            *
 * Function :   char *get_lookup_filename(int apid )                         *
 * Description: Get lookup filename to based on apid. Return a pointer to    *
 *              string. String is <instrument-name>-<packet-name> where      * 
 *              instrument is HMI,AIA,or SDO and packet-name is ISP, ADP, or *
 *              SEQ.                                                         *
 *****************************************************************************/
char *get_lookup_filename(int apid ) 
{
   /* declarations */
   char *lufn;
   char lookup_filename[HK_LEV0_MAX_DSNAME_STR];
   char pkt_name[HK_LEV0_MAX_PACKET_NAME];

   /* initialized variables */
   lufn = lookup_filename;

   /* determine value of packet name to use for data series */
   /* set packet name based on the APID. */
   switch(apid)
   {  
   /*****************************************/
   /***  Case of Finding HMI ISP APIDS    ***/
   /*****************************************/
   case HK_LR_HMI_ISP:
   case HK_HSB_HMI_ISP_1:
   case HK_HSB_HMI_ISP_2:

        /* create prefix HMI-ISP using environment variable and define value */
        strcpy(pkt_name, HK_HMI_PKT_NAME_ISP);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_HMI"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   case HK_HSB_IRIS_ISP:
        //IRIS 
        strcpy(pkt_name, HK_HSB_IRIS_PKT_NAME_ISP);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_IRIS"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /*****************************************/
   /***  Case of Finding AIA ISP APIDS    ***/
   /*****************************************/
   case HK_LR_AIA_ISP:
   case HK_HSB_AIA_ISP_1:
   case HK_HSB_AIA_ISP_2:

        /* create prefix AIA-ISP */
        strcpy(pkt_name,HK_AIA_PKT_NAME_ISP);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_AIA"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /*****************************************/
   /***Case of Finding HMI SEQUENCER APIDS***/
   /*****************************************/
   case HK_LR_HMI_SEQ:
   case HK_HSB_HMI_SEQ_1:
   case HK_HSB_HMI_SEQ_2:

        /* create prefix HMI-SEQ */
        strcpy(pkt_name,HK_HMI_PKT_NAME_SEQ);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_HMI"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /*****************************************/
   /***Case of Finding AIA SEQUENCER APIDS***/
   /*****************************************/
   case HK_LR_AIA_SEQ:
   case HK_HSB_AIA_SEQ_1:
   case HK_HSB_AIA_SEQ_2:

        /* create prefix AIA-SEQ */
        strcpy(pkt_name,HK_AIA_PKT_NAME_SEQ);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_AIA"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /*****************************************/
   /***Case of Finding HMI OBT APIDS      ***/
   /*****************************************/
   case HK_LR_HMI_OBT:
   case HK_HSB_HMI_OBT_1:
   case HK_HSB_HMI_OBT_2:

        /* create prefix HMI-OBT */
        strcpy(pkt_name,HK_HMI_PKT_NAME_OBT);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_HMI"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /*****************************************/
   /***Case of Finding AIA OBT APIDS      ***/
   /*****************************************/
   case HK_LR_AIA_OBT:
   case HK_HSB_AIA_OBT_1:
   case HK_HSB_AIA_OBT_2:

        /* create prefix AIA-OBT */
        strcpy(pkt_name,HK_AIA_PKT_NAME_OBT);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_AIA"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /******************************************************/
   /***Case of Finding SDO ANCILLORY SCIENCE DATA APIDS***/
   /******************************************************/
   case HK_LR_SDO_ASD:

        /* create prefix SDO-ASD */
        strcpy(pkt_name,HK_SDO_PKT_NAME_ASD);
        sprintf(lufn,"%s-%s",getenv("HK_PREFIX_JSVNMAP_FILENAME_SDO"),make_strupr(pkt_name)); 
        strcat(lufn,"\0");
        break;

   /*************************************************/
   /*** All Other cases of undefined packet names ***/
   /*************************************************/
   default:
        /* create default prefix using apid in name */
        if(apid <= HK_HSB_HIGHEST_HMI_APID && apid >= HK_HSB_LOWEST_HMI_APID  || (apid <= HK_LR_HIGHEST_HMI_APID  && apid >= HK_LR_LOWEST_HMI_APID))
        {
          sprintf(lufn,"%s-%-03.3d",getenv("HK_PREFIX_JSVNMAP_FILENAME_HMI"),apid); 
          strcat(lufn,"\0");
        }
        else if((apid <= HK_HSB_HIGHEST_AIA_APID && apid >= HK_HSB_LOWEST_AIA_APID) || (apid <= HK_LR_HIGHEST_AIA_APID && apid >= HK_LR_LOWEST_AIA_APID))
        {
          sprintf(lufn,"%s-%-03.3d",getenv("HK_PREFIX_JSVNMAP_FILENAME_AIA"),apid); 
          strcat(lufn,"\0");
        }
        else if((apid <= HK_LR_HIGHEST_SDO_APID) && (apid >= HK_LR_LOWEST_SDO_APID))
        {
          sprintf(lufn,"%s-%-03.3d",getenv("HK_PREFIX_JSVNMAP_FILENAME_SDO"),apid); 
          strcat(lufn,"\0");
        }
        else
        {
          printkerr("Warning at %s, line %d: APID is not in range of "
                    "HMI, AIA, or SDO. APID: <%d>\n",
                    __FILE__, __LINE__, apid);
          strcpy(lufn,"\0");
        }
        break;

   }// switch

   return(lufn);
}



/*****************************************************************************
 * Name:        MAKE STRING UPPER CASE - make characters in string upper case*
 * Function :   char *make_strupr(char *a )                                  *
 * Description: Make string passed in function upper case. Return a pointer  *
 *              to upper cased string value.                                 * 
 *****************************************************************************/
char * make_strupr (char *a)
{
  char *ret = a;

  while (*a != '\0')
  {
    if (islower (*a))
     *a = (char)toupper (*a);
    ++a;
  }

  return ret;
}



/************************************************************/
/* NAME:CHECK HK RECORD EXISTS                              */
/* FUNCTION:   check_hk_record_exists                       */
/*             (char* ds_name, HK_Keyword_t *kw, int apid)  */
/* DESCRIPTION: Checks if record already exists in data     */ 
/*              series by check time code value keywords    */
/*              Also keeps in-memory the time codes for data*/
/*              series in link list                         */
/************************************************************/
int check_hk_record_exists(char* ds_name, HK_Keyword_t *kw, int apid)
{
  /* variable definitions */
  HK_DSN_RANGE_t *dsr;//current link list pointer
  HK_DSN_RANGE_t *pdsr;//previous link list pointer
  HK_DSN_RANGE_t *dsrn;//next link list pointer
  HK_Keyword_t *tkw;  //kw link list pointer
  HK_Timecode_t *tc;  //current link list pointer
  HK_Timecode_t *ptc; //previous link list pointer
  HK_Timecode_t *tcn; //next link list pointer
  int ck_status;
  int found_status=0;
  int j;
  long int sec_pkt, subsec_pkt;
  long int trigger_lr, trigger_hr;
  static int tc_last_index;
  int found_dsr_node=0;
  int ckinit=0;
  int istatus=0;

  /* check if should initialize cache structure for this data series */
  ckinit= check_ds_cache_init(ds_name); 

  /* open records once for complete dayfile of packets */
  /* if  should initialize then create cache structure for this data series */
  if(ckinit)
  {
    istatus=initialize_timecodes_cache(ds_name, kw, apid);
    if(istatus == 0)
    {
      return (0);
    }
  }

  /* get timecode seconds and subseconds in current packet */
  /* loop until get TIMECODE Seconds and Subseconds values */
  tkw=kw;
  while (tkw)
  {
    if(strstr(tkw->name,"TIMECODE_SECONDS"))
    {
      /*set found flag*/
      /* create time string based on some hard coded parameters for now */
      sec_pkt = tkw->eng_value.uint32_val;
    }
    if(strstr(tkw->name,"TIMECODE_SUBSECS"))
    {
      /*set found flag*/
      /* create time string based on some hard coded parameters for now */
      subsec_pkt = tkw->eng_value.uint32_val;
      break;
    }
    tkw= tkw->next;
  }



  /* check for cache link list dsr node for current packets dsname then check low range values in dsr node */
  for(dsr=dsr_head; dsr;dsr=dsr->next)
  {
     /* check for matching data series names */
    if (!strcmp(dsr->dsname, ds_name))
    {

      /* found ds name in dsr cache link list-now check low range threshold for holding cache values is okay */
      /* check if want to force free of cache, then query  and reload of timecodes with new range of timecodes */
      trigger_lr=dsr->timecode_lrsec + 3600; //Trigger reload an hour later for low range
      trigger_hr=dsr->timecode_hrsec - 3600; //Trigger reload an hour early for high range
      if ( (trigger_lr <= sec_pkt) && (trigger_hr >= sec_pkt) )
      {
        /* if got timecode seconds value less than low range -then values in cache ok */
        /* if got packet timecode seconds less than high range and greater than low range  -then values in cache ok */
        break;
      }
      else
      {
        /* if current time code is greater than low range or if current time code is greater than high range */
        /* then free link list in order to reload cached values. this frees cache every two days */
#ifdef DEBUG_WRITE_HK_TO_DRMS
        printkerr("Warning Monitor Message:check_hk_rec_exist:Watch this does "
                  "not occur too often. The low and high range of timecode cache needs "
                  "to be reset and reloaded.This warning message should occur only about four "
                  "times every 24 hours per HK by APID series.\n");
#endif
        /* reset range of actual query of timecodes to database -we currently doing 12 hour period*/
        /* HK_SECONDS_RANGE was changed to 6 hrs to fix bug on retrieving too large amount of data*/
        dsr->timecode_lrsec= sec_pkt - (HK_SECONDS_RANGE);
        dsr->timecode_hrsec= sec_pkt + (HK_SECONDS_RANGE);

        /* free link */
        for(j=0,tc=dsr->tcnode;tc;j++)
        {
          ptc=tc;
          tc=tc->next;
          free(ptc);
        }
        dsr->tcnode=NULL;
#ifdef DEBUG_WRITE_HK_TO_DRMS
       printkerr("DEBUG:MESSAGE: Freed <%d>  TIMECODE nodes in cache for <%s> before reloading cache TIMECODE values\n", j,ds_name);
#endif
        tc_last_index=0;
        /* reload link list done below with new values */
        //dsr_head=NULL;Took out since might cause memory leaks
        istatus=initialize_timecodes_cache(ds_name, kw, apid);
        if(istatus == 0)
        {
            return (0);
        }
        break;
      } /* end-else curent timecode value not less the low range in dsr cache node */
    }/* end-if found match for current data series name in cache link list */
  }/*end-for loop - look for dsname node and check low range */



  /* find match of time codes in current packet with timecodes in timecodes cache link list */
  /* if match set flag to not create new record in DRMS, else set flag to create new record */
  /*if found in link list -set to 1- don't create new record in drms*/

  /* loop thru data series cache link list */
  for(ck_status=found_status=0, dsr=dsr_head; dsr; dsr=dsr->next)
  {
    /* check if found matching data series name in current packet and cache link list */
    if(!strcmp(dsr->dsname, ds_name))
    {

      /* found dsr(data series range) cache node with ds_name */
      /* go to first timecode node in dsr node and loop thru timecodes in cache link list*/
      for(found_status=0, j=0, tc=dsr->tcnode;  tc;j++)
      { 
        /* check if found match for seconds value */
        if( sec_pkt == tc->sec)
        { 
          /* found second, check if found subseconds */
          if( subsec_pkt == tc->subsec)
          {
            /*its there set found_status to skip writing */
            ck_status=1;
            found_status=1;
            tc=tcn;
            break;
          }
          else
          {
            /*did not find yet-check next node in timecodes link list */
            /*its not there- set found_status to create new record */
            found_status=0;
            ck_status=0;
            tcn=tc->next;
            tc=tcn;
          }
        }
        else
        {  
          /* did not find-check next node in timecodes link list */
          found_status=0;
          ck_status=0;
          tcn=tc->next;
          tc=tcn;
        }
      }/*end-for-loop- to loop thru timecode link list looking for second and subsec match current packet's timecode */

    }/*if found ds name in dsr node */
  }/*loop thru all dsr nodes */



  /* if did not find timecode then will add to drms-therefore add values to timecodes link list */
  if(!found_status)
  { 
    /* cannot find record with these timecodes */
    /* add to timecodes link list since will be writing to DRMS */
    ck_status=found_status;

    /* check if need to create new dsr node*/
    for(dsr=dsr_head; dsr; dsr=dsr->next)
    {
      if(!strcmp(dsr->dsname, ds_name))
      {
        /* got dsr node matching data series name */
        found_dsr_node=1;
        break;
      }
    } /*end-for-loop thru dsr node looking for data series name match*/

    /* check if found node with ds_name  */
    if (!found_dsr_node)
    {
      /* if did not find dsr node then create new dsr node!! */
      /*create node for new dsr node */
      dsrn= (HK_DSN_RANGE_t *)malloc(sizeof(HK_DSN_RANGE_t));

      /* set next (ndsr) node value */
      strcpy(dsrn->dsname, ds_name);
      dsrn->timecode_lrsec= sec_pkt - HK_SECONDS_RANGE;
      dsrn->timecode_hrsec= sec_pkt + HK_SECONDS_RANGE;
      dsrn->next=NULL;
      dsrn->tcnode=NULL;

      /* check if other dsr node and go to last dsr node and add new dsr node*/
      for(dsr=dsr_head; dsr; pdsr=dsr, dsr=dsr->next);//pdsr is last

      /* if did not find because no dsr nodes with dsname value */
      if(dsr_head)
      {
        /* link previous dsr node with new dsr node*/
        pdsr->next= dsrn;
      }
      else
      {
        /* else did not find because no dsr node there- no head dsr node */
        dsr_head=dsrn;
      }
      tc=dsrn->tcnode;

    }/* end-if did not dsr node!! */
    else
    {
      /* if found dsr node set tc using dsr */
      tc=dsr->tcnode;
      dsrn=dsr;/* fix for ticket #278*/
    }

    /* check if head tc node there, if not there create head tc node */
    if (!tc)
    {
       /* malloc new space for head tc node and set values in timecode node */
       tc= (HK_Timecode_t *)malloc(sizeof(HK_Timecode_t));
       tc->sec = sec_pkt;
       tc->subsec = subsec_pkt;
       tc->next = NULL;
       dsrn->tcnode=tc;
       tc_last_index++;
    } /*end-if no head timecode node then create one */
    else
    {
       /* malloc space for new timecode node -  go to end of last node */
       while(tc) 
       {
           ptc=tc; /*save previous node to link nodes below */
           tcn=tc->next;
           tc=tcn;
       }

       /* malloc new space for next node */
       tc= (HK_Timecode_t *)malloc(sizeof(HK_Timecode_t)) ;
       tc->sec = sec_pkt;
       tc->subsec = subsec_pkt;
       tc->next = NULL;
       ptc->next=tc; /*link nodes in list */
       tc_last_index++;//use to count nodes created to compare with freed
    }/*end-else malloc timecode node and add to end of timecode link list */
  }
  else
  {
    /*found in link list -return ck_status=1 to skip creating new record */
    ck_status=found_status;
  }
  /* return check status */
  return (ck_status);
} /*end check_hk_records_exist */



/************************************************************/
/* NAME:CREATE CACHE LOW RANGE                              */
/* FUNCTION:   short create_low_range                       */
/*                   (HK_Keyword_t *kw_head,int *lr)        */
/* DESCRIPTION: Create low range for cached value for data  */ 
/*              series. Create low range by using current   */
/*              timecode in seconds - HK_SECONDS_RANGE.     */
/*              Pass in current hk structure to find time   */
/*              of first packet. Pass back low range value. */
/************************************************************/
short create_low_range(HK_Keyword_t *kw_head, long int *lr)
{
  HK_Keyword_t *tkw;
  tkw=kw_head;

  /* loop thru kw structure of keywords to find timecode seconds values */
  while (tkw)
  {
    if(strstr(tkw->name,"TIMECODE_SECONDS"))
    {
      /*if found timecode seconds keyword then set found flag by returning 1 */
      /* create time string based on some hard coded parameters for now */
      *lr = tkw->eng_value.uint32_val - HK_SECONDS_RANGE;
      return(1);
    }
    tkw= tkw->next;
  }
  /* if did not set return not found by returning 0 value */
  return(0);
}
/************************************************************/
/* NAME:CREATE CACHE HIGH RANGE                             */
/* FUNCTION:   short create_high_range                      */
/*                   (HK_Keyword_t *kw_head,int *hr)        */
/* DESCRIPTION: Create high range for cached value for data */ 
/*              series. Create high range by using current  */
/*              timecode in seconds + (HK_SECONDS_RANGE)*/
/*              Pass in current hk structure to find time   */
/*              of last packet. Pass back high range value. */
/************************************************************/
short create_high_range(HK_Keyword_t *kw_head, long int *hr)
{
  HK_Keyword_t *tkw;
  tkw=kw_head;

  /* loop thru kw structure of keywords to find timecode seconds values */
  while (tkw)
  {
    if(strstr(tkw->name,"TIMECODE_SECONDS"))
    {
      /*if found timecode seconds keyword then set found flag by returning 1 */
      /* create time string based on some hard coded parameters for now */
      *hr = tkw->eng_value.uint32_val + (HK_SECONDS_RANGE);
      return(1);
    }
    tkw= tkw->next;
  }
  /* if did not set return not found by returning 0 value */
  return(0);
}

/************************************************************/
/* NAME:CHECK DATA SERIES CACHE INITIALIZE                  */
/* FUNCTION:   check_ds_cache_init(char* ds_name            */
/* DESCRIPTION: Checks if cache has been initialize for data*/ 
/*              series by first checking the array flag for */
/*              for data series name. If there, then return */
/*              should_init value of 0. If not there, return*/
/*              should_init value of 1 and add data series  */
/*              name to array                               */
/************************************************************/
short check_ds_cache_init(char *dsname) 
{
  short i;
  static char init_cache[100][50];
  static short init_count=0;
  short should_init=0;

  /* if first time running- return back should_init cache */
  if (init_count == 0)
  {
    should_init=1;
  }
  else
  {
    /* check if data series name in array- check init_count array idexes */
    for (i=0,should_init=1; i < init_count;i++) 
    {

       if (!strcmp(init_cache[i], dsname))
       {
          /* if there- this should_init is no or return 0 */
          should_init=0;
          break; 
       }
    }
  }
  /*add value to init_cache array if never initialized before */
  if (should_init)
  {
    strcpy(init_cache[init_count],dsname); 
    init_count++;
  } 
  return(should_init);
}

/************************************************************/
/* NAME:CHECK CURRENT MAP FILES To USE                      */
/* FUNCTION:   check_curr_map_file(char* pv_str, short      */
/*                                 *curr_map_flag)          */
/* DESCRIPTION: Check packet version number to determine    */ 
/*              whether to use the old JSOC map files or    */
/*              new JSOC map by passing back flag  set to   */
/*              1 for reload map files or 0 to just use     */
/*              map files currently loaded. This is to      */
/*              ensure backward capatability                */
/************************************************************/
short check_curr_map_file(char* pv_str, short *curr_map_flag)
{
  float pv_num;
  short new_map_flag;
  float merge_threshold;

  /* get MERGE threshold packet version number value */
  merge_threshold= (float)HK_LEV0_START_MERGED_PVNW;
  merge_threshold += (float)(HK_LEV0_START_MERGED_PVND/1000.0);

  /* get packet version number float value from string value */
  sscanf(pv_str,"%f",&pv_num);
  //if(pv_num < 1.194)
  if(pv_num < merge_threshold)
  {
     new_map_flag= HK_NON_MERGE_MAP_FILE_FLAG;
  }
  //else if (pv_num > 1.194)
  else if (pv_num >= merge_threshold)
  {
     new_map_flag= HK_MERGE_MAP_FILE_FLAG;
  }
  else
  {
    printkerr("DEBUG:Warning Message:check_curr_map_file: Should not have reached else case:\n");
  }

  /* check if switch of map flag occurred- if did return 1 trigger to init map files */
  if (new_map_flag == *curr_map_flag)
  {

    return(0);
  }
  else
  {
   *curr_map_flag = new_map_flag; 
    return(1);
  }
}



/************************************************************/
/* NAME:GET QUERY RANGE TO USE TO CHECK RECORDS EXIST       */
/* FUNCTION:   get_query_range()                            */
/* DESCRIPTION: Get query range to use to do query of       */ 
/*              of records to load into cache memory which  */
/*              before loading to dataseries                */
/************************************************************/
int get_query_range(int range_type, HK_Keyword_t *kw, char qr[HK_MAX_SIZE_RANGE_TIME])
{
  int gpt_status;
  TIME pkt_time;
  TIME *ptime=&pkt_time;
  TIME pkt_time1;
  TIME *ptime1=&pkt_time1;
  int yr,mr,dr,hr;

  /* get current packet time in current packet */
  gpt_status = get_pkt_time_from_timecodes(kw, ptime);
  if(!gpt_status)
  {
    printkerr("WARNING at %s, line %d: Could not find timecode inorder "
              "to get packet time. gpt_status:<%d>.\n", __FILE__,__LINE__, gpt_status);
  }
  /* using current packet time create new low or high range */
  if(range_type == HK_LOW_QUERY_RANGE)
  {
    *ptime1= *ptime - (HK_SECONDS_RANGE);
  }
  else if (range_type == HK_HIGH_QUERY_RANGE)
  {
    *ptime1= *ptime + (HK_SECONDS_RANGE);
  }
  else
  {
    sprintf(qr, "0000.00.00_00:00:00.00_UTC");

#ifdef DEBUG_WRITE_HK_TO_DRMS
      printkerr("DEBUG:Message at %s, line %d: Bad Query range "
                " is <%s>.\n", __FILE__,__LINE__, qr);
#endif

    return (0);
  }
  hr=get_hour_from_pkttime(*ptime1);
  dr=get_day_from_pkttime( *ptime1);
  mr=get_month_from_pkttime(*ptime1);
  yr=get_yr_from_pkttime(*ptime1);
  sprintf(qr, "%04d.%02d.%02d_%02d:00:00.00_UTC",yr,mr,dr,hr);

#ifdef DEBUG_WRITE_HK_TO_DRMS
      printkerr("DEBUG:Message at %s, line %d: Query range "
                " is <%s>.\n", __FILE__,__LINE__, qr);
#endif

  return (1); 
}

/*###################################################################
#  initialize_timecodes_cache()                                   #
###################################################################*/
int initialize_timecodes_cache(char* ds_name, HK_Keyword_t *kw, int apid)
{
  /* variable definitions */
  HK_DSN_RANGE_t *dsr;//current link list pointer
  HK_DSN_RANGE_t *pdsr;//previous link list pointer
  HK_Timecode_t *tc;  //current link list pointer
  HK_Timecode_t *tcn; //next link list pointer
  char pkt_ver_num[HK_LEV0_MAX_PVN_STR];
  char tc_sec[HK_LEV0_MAX_LONG_KW_NAME];
  char tc_subsec[HK_LEV0_MAX_LONG_KW_NAME];
  char strapid[HK_LEV0_MAX_PACKET_NAME];
  char nquery[HK_MAX_SIZE_QUERY];
  char qr1[HK_MAX_SIZE_RANGE_TIME];
  char qr2[HK_MAX_SIZE_RANGE_TIME];
  char *pvn;
  int curr_pvn_dn;
  int curr_pvn_wn;
  int drms_status;
  int j;
  int lr_status;
  int hr_status;
  int gqr1_status;
  int gqr2_status;
  static long int tclr=0;
  static long int tchr=0;
  /* drms record create variables */
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;

  /* initialize variables */
  rs=NULL;
  rec=NULL;
  pvn= pkt_ver_num;

  /* get packet version number to use for TIMECODE keyword name */
  if(check_for_sdo_apid(apid))
  {
    /* for sdo 129 set pvn to merged case to pickup SDO-ASD-PVN-TO-JSVN file */
    /* also to go into merged case logic below */
    curr_pvn_wn= HK_LEV0_START_MERGED_PVNW;
    curr_pvn_dn= HK_LEV0_START_MERGED_PVND;
  }
  else
  {
    pvn=get_packet_version_number(kw);
    sscanf(pvn,"%3d.%3d",&curr_pvn_wn, &curr_pvn_dn);
  }

  /* drms record create variables */
  /* set first time run timecode low range */
  lr_status=create_low_range( kw, &tclr);
  hr_status=create_high_range( kw, &tchr);
/************************************************By JA for IRIS**********
  if( !lr_status || !hr_status)
  {
    printkerr("Warning at %s, line %d: Could not set low and/or high range value in cache memory. "
              "This is use to save in cache a limited amount of records in cache. "
              "This cache is used to test if record already exists in data series. \n",
              __FILE__, __LINE__);
  }
***************************************************************************/

  /* setup query get narrow range of records in data series to avoid trucated query error */
  /* use current packet time found in first packet to get low and high range for drms_open_records query */
  /* set query range low current day - 1 day and set query range high to current day + 1 day */
  gqr1_status=get_query_range(HK_LOW_QUERY_RANGE,kw,qr1);
  gqr2_status=get_query_range(HK_HIGH_QUERY_RANGE,kw,qr2);
  if(!gqr1_status || !gqr2_status)
  {
    printkerr("Warning at %s, line %d: Could not create query range successfully. "
              "qr1:<%s> qr2:<%s>\n", __FILE__, __LINE__);
  }
  /* create full query statement */
  sprintf(nquery,"%s[? PACKET_TIME < $(%s) AND PACKET_TIME >  $(%s) ?]",ds_name,qr2,qr1);
  sprintf(ispquery, "%s", nquery);	//save for opening to get kw for lev0
  /* open records for drms series */
  rs = drms_open_records(drms_env, nquery, &drms_status);
  if(drms_status) 
  {
    printkerr("ERROR at %s, line %d: Can't open records using drms_open_records for drms query <%s>  "
              "This series may not exist. Please create series and rerun\n",
                __FILE__,__LINE__, nquery);
    return(0);
  }

#ifdef DEBUG_WRITE_HK_TO_DRMS
  if(!rs)
  {
    printkerr("DEBUG:Message at %s, line %d: Record query to "
              "dataseries is <%s>. Return status of query:<%d> \n",
              __FILE__,__LINE__, nquery, drms_status);
  }
  else
  {
    printkerr("DEBUG:Message at %s, line %d: Record query to "
              "dataseries is <%s>. Return status of query:<%d> "
              "Return count:<%d>\n", __FILE__,__LINE__, nquery, drms_status, rs->n);
  }
#endif

  /* open records for drms series */
  /* if have null no record in series so return 0 or did not find any records in series so can write to series */
  if (!rs) 
  { 
    /* if no records there, break from if, create dsr node and timecode node,then return 0 but got 0 records maybe because hit query limit */
    printkerr("Warning at %s, line %d: There are no records for series. "
              "Therefore returning 0 for check if hk record exists. When return "
              "zero then the records will get written to series. \n", __FILE__, __LINE__);
    /* if get query 10001 this is drms query trucated therefore size of query too large so adjust time range */
    /* This could be because series was just created and has no records */
    printkerr("Warning at %s, line %d: DRMS data series <%s> returned status<%d>.Query "
              "could be too big and reached limit. The query was <%s>\n", 
              __FILE__,__LINE__, ds_name, drms_status, nquery);
    return(0); 
  }

  /* ds has records, so fill records in timecode link list and check new records in link list */
  if(rs->n)
  {

    /* have  already existing records in data series so set timecount link list seconds and subsecs from each record */
    /* if dsr head node exists */
    if(dsr_head)
    {
      /* go to last node via pdsr */
      for(dsr=dsr_head;dsr; pdsr=dsr,dsr=dsr->next)
      {
         /* new code to fix TRAC ticket 278*/
         /* check if HK_DSN_RANGE_t node exists */
         if(!strcmp(dsr->dsname, ds_name))
         {
            /*found HK_DSN_RANGE_t node with series name */
            tc=dsr->tcnode; /*set tc to null since dsr->node is set to null during free of node above*/
            break;
         }
      }/*end for loop threw HK_DSN_RANGE_t nodes */
    }
    else
    {
       /* if dsr_head does NOT exist, then assign dsr null value */
       dsr=dsr_head;
    }

    /* get apid packet value(isp,seq or 481) based on apid -use to create TIMECODE keywords names later*/
    if(curr_pvn_wn  >=  HK_LEV0_START_MERGED_PVNW && curr_pvn_dn >= HK_LEV0_START_MERGED_PVND)
    {
      /* doing setting of packet value name for merged TIMECODE names- ISP,SEQ,etc.*/
      if(apid == HK_HSB_HMI_ISP_1 || apid == HK_HSB_HMI_ISP_2 || apid == HK_LR_HMI_ISP)
      {
        strcpy(strapid, HK_HMI_PKT_NAME_ISP);
      }
      else if(apid == HK_HSB_IRIS_ISP) {
        strcpy(strapid, HK_HSB_IRIS_PKT_NAME_ISP);
      }
      else if(apid == HK_HSB_AIA_ISP_1 || apid == HK_HSB_AIA_ISP_2 || apid == HK_LR_AIA_ISP )
      {
        strcpy(strapid, HK_AIA_PKT_NAME_ISP);
      }
      else if(apid == HK_HSB_HMI_SEQ_1 || apid == HK_HSB_HMI_SEQ_2 || apid ==HK_LR_HMI_SEQ )
      {
        strcpy(strapid, HK_HMI_PKT_NAME_SEQ);
      }
      else if(apid == HK_HSB_AIA_SEQ_1 || apid == HK_HSB_AIA_SEQ_2 || apid == HK_LR_AIA_SEQ )
      {
        strcpy(strapid, HK_AIA_PKT_NAME_SEQ);
      }
      else if(apid == HK_HSB_HMI_OBT_1 || apid == HK_HSB_HMI_OBT_2 || apid ==HK_LR_HMI_OBT )
      {
        strcpy(strapid, HK_HMI_PKT_NAME_OBT);
      }
      else if(apid == HK_HSB_AIA_OBT_1 || apid == HK_HSB_AIA_OBT_2 || apid == HK_LR_AIA_OBT )
      {
        strcpy(strapid, HK_AIA_PKT_NAME_OBT);
      }
      else if(apid == HK_LR_SDO_ASD)
      {
        strcpy(strapid, HK_SDO_PKT_NAME_ASD);
      }
      else
      {
        /* for apid using hex number for timecode keywords like apid 129 and non-isp and non-seq apids */
        sprintf(strapid,"%03.3x", apid);
      }
    }
    else
    {
      /* doing setting of packet name based on non-merge TIMECODE hex names - 1BD,1DB,1E1,etc. */
      sprintf(strapid,"%03.3x", apid);
    }
        
    /* go thru each record in data series and get and set timecode sec and subsec values */
    for(j=0; j < rs->n;j++)
    {

      /* if no head node create first head node -which was checked andn set above*/
      if(!dsr)
      {
        /* create head node for data series name and range link list */
        dsr= (HK_DSN_RANGE_t *)malloc(sizeof(HK_DSN_RANGE_t));

        /* create head node for time code link list */
        tc= (HK_Timecode_t *)malloc(sizeof(HK_Timecode_t));
       
        /* set dsr node value to allocated space for node */
        if(!dsr_head) 
        {
          /* if dsr_head null then set dsr_head to newly created dsr node */
          dsr_head=dsr;
        }
        else 
        {
          /* else dsr_head does exist, so use the last node in dsr link list to add this newly created node */
          pdsr->next=dsr;
        }

        /* set values in dsr node, name, low range value, etc */
        strcpy(dsr->dsname, ds_name);
        dsr->next=NULL;
        dsr->tcnode=tc;
        dsr->timecode_lrsec=tclr;
        dsr->timecode_hrsec=tchr;

        /* read record and fill in value or cache value in link list in order to not go to db to get data again */
        rec = rs->records[j];

        /* get time code keyword strings to use to lookup value in rec */
        sprintf(tc_sec, "%s%s_%s","APID",make_strupr(strapid),"TIMECODE_SECONDS");
        sprintf(tc_subsec, "%s%s_%s","APID",make_strupr(strapid),"TIMECODE_SUBSECS");

        /* get and set time code second and subsecond values */
        tc->sec = drms_getkey_longlong(rec, tc_sec, &drms_status);
        tc->subsec = drms_getkey_longlong(rec,tc_subsec , &drms_status);

        /* set next timecode node */
        tc->next = NULL;

      }
      else /* else add newly created timecode node to end of dsr link list */
      {
        /* create next node for time code link list */
        tcn= (HK_Timecode_t *)malloc(sizeof(HK_Timecode_t));

        /* new code for fix of TRAC ticket 278 */
        /* check if HK_Timecode_t node is null-if is set top node else link last node to newly created node*/
        if(dsr->tcnode == NULL)
        {
           /* set top node for HK_Timecode_t nodes */
           dsr->tcnode=tcn;
        }
        else
        {
          /* assign next value to new node and reset tc pointer */
          tc->next=tcn;
        }
        /* end of new code for fix of TRAC ticket 278 */
        /* set tc to node going to add timecode values to */
        tc=tcn;

        /* get next record and set values in timecode link list */
        rec = rs->records[j];

        /* get time code keyword strings */
        sprintf(tc_sec, "%s%s_%s","APID",make_strupr(strapid),"TIMECODE_SECONDS");
        sprintf(tc_subsec, "%s%s_%s","APID",make_strupr(strapid),"TIMECODE_SUBSECS");

        /* get time code values */
        tc->sec = drms_getkey_longlong(rec, tc_sec, &drms_status);
        tc->subsec = drms_getkey_longlong(rec,tc_subsec , &drms_status);

        /* set next value */
        tc->next = NULL;
      }
    }/*for-loop -go thru each rec and get timecodes and set in cache link list*/
  } /* if-end of if record in data series are  > 0 */
  /*else - don't have any records in data series to load in dsr cache link list */

  /* close once for complete dayfile of packets */
  drms_status = drms_close_records(rs, DRMS_FREE_RECORD);

  return (1);

} 

/************************************************************/
/* NAME:CHECK PACKET TIME RECORD IS WITH RANGE TO PROCESS   */
/* FUNCTION:   check_hk_record_within_time_range()          */
/* DESCRIPTION: Check PACKET_TIME of packet is within time  */
/*              of current-time + 12 hours. If out of range */
/*              then return 0 for do not process.           */
/************************************************************/
int check_hk_record_within_time_range( HK_Keyword_t *kw)
{
  TIME pkt_time;
  TIME *ptime= &pkt_time;
  HK_Keyword_t *tkw;  //kw link list pointer
  int curr_packet_number;
  int curr_time_range;
  int curr_hr_range;
  int y,m,d,h;
  int gpt_status;
  struct tm timestorage;
  struct tm *time_ptr;
  time_t tvalue;
  time_ptr= &timestorage;

  /* get year, month and hour for CURRENT packet */
  tkw=kw;
  gpt_status = get_pkt_time_from_timecodes(tkw, ptime);
  if(!gpt_status)
  {
    printkerr("WARNING at %s, line %d: Could not find timecode inorder "
              "to get packet time. gpt_status:<%d>.\n", __FILE__,__LINE__, gpt_status);
  }
  d=get_day_from_pkttime(*ptime);
  h=get_hour_from_pkttime(*ptime);
  m=get_month_from_pkttime(*ptime);
  y=get_yr_from_pkttime(*ptime);

  /* set up packet time found as a integer number  to compare later -yyyymmdd*/
  curr_packet_number = (y * 10000) + (m * 100) + d ; 
#ifdef DEBUG_WRITE_HK_TO_DRMS
  printkerr("DEBUG:Message at %s, line %d: Packet time is::: y:<%d>m:<%d>d:<%d>h<%d>.\n",__FILE__,__LINE__,y,m,d,h);
#endif

  /* get CURRENT TIME now in year,month, day and hour  */
  tvalue = time(NULL);
  time_ptr = gmtime(&tvalue);
#ifdef DEBUG_WRITE_HK_TO_DRMS
  printkerr("DEBUG:Message at %s, line %d: Current time is:::y:<%d>m:<%d>d:<%d>h<%d>.\n", 
            __FILE__,__LINE__, time_ptr->tm_year+1900,time_ptr->tm_mon+1,time_ptr->tm_mday,time_ptr->tm_hour);
#endif

  /*adjust by 12 hours the Current time and setup as a integer number,yyyymmdd, to compare later with packet time*/
  if((time_ptr->tm_mon+1) == 12 &&  (time_ptr->tm_mday == 31) && (time_ptr->tm_hour >= 12))
  {
    /* adjust year,month,day if end of year range */
    curr_time_range = (((time_ptr->tm_year+1900+1) * 10000) + ( 1 * 100) + 1);
    curr_hr_range=time_ptr->tm_hour - 12;
  }
  else if( ( (time_ptr->tm_mon+1) == 1 || (time_ptr->tm_mon+1) == 3 ||  (time_ptr->tm_mon+1) == 5 ||  (time_ptr->tm_mon+1) == 7 ||  (time_ptr->tm_mon+1) == 8  ||  (time_ptr->tm_mon+1) == 10  ) && time_ptr->tm_mday == 31 && (time_ptr->tm_hour >= 12))
  {
    /* adjust month,day if end of month range */
    curr_time_range = (((time_ptr->tm_year+1900) * 10000) + ((time_ptr->tm_mon+2) * 100) + 1);
    curr_hr_range=time_ptr->tm_hour - 12;
  }
  else if( ( (time_ptr->tm_mon+1) == 4 || (time_ptr->tm_mon+1) == 6 ||  (time_ptr->tm_mon+1) == 5 ||  (time_ptr->tm_mon+1) == 9 ||  (time_ptr->tm_mon+1) == 11 ) && time_ptr->tm_mday == 30 && (time_ptr->tm_hour >= 12))
  {
    /* adjust month,day if end of month range */
    curr_time_range = ((time_ptr->tm_year+1900) * 10000) + ((time_ptr->tm_mon+2) * 100) + 1;
    curr_hr_range=time_ptr->tm_hour - 12;
  }
  else if( ( (time_ptr->tm_mon+1) == 2 && (time_ptr->tm_year+1900 == 2009) && time_ptr->tm_mday == 28) &&  (time_ptr->tm_hour >= 12))
  {
    /* adjust month,day if end of month range */
    curr_time_range = ((time_ptr->tm_year+1900) * 10000) + ((time_ptr->tm_mon+2) * 100) + 1;
    curr_hr_range=time_ptr->tm_hour - 12;
  }
  else if( ( (time_ptr->tm_mon+1) == 2 && (time_ptr->tm_year+1900 == 2010) && time_ptr->tm_mday == 29) &&  (time_ptr->tm_hour >= 12))
  {
    /* adjust month,day if end of month range */
    curr_time_range = ((time_ptr->tm_year+1900) * 10000) + ((time_ptr->tm_mon+2) * 100) + 1;
    curr_hr_range=time_ptr->tm_hour - 12;
  }
  else if((time_ptr->tm_hour >= 12))
  {
    /* adjust day if end of day range */
    curr_time_range= ((time_ptr->tm_year+1900) * 10000) + ((time_ptr->tm_mon+1) * 100) + ( time_ptr->tm_mday + 1);
    curr_hr_range=time_ptr->tm_hour - 12;
  }
  else if((time_ptr->tm_hour < 12))
  {
    /* adjust hour if same day range */
    curr_time_range= ((time_ptr->tm_year+1900) * 10000) + ((time_ptr->tm_mon+1) * 100) + ( time_ptr->tm_mday);
    curr_hr_range=time_ptr->tm_hour + 12;
  }
  else
  {
    printkerr("ERROR:check_hk_record_within_time_range:Bad else case in function\n");
  }

  /* compare packet time to current time range calculated above */
  if ( curr_packet_number < curr_time_range)
  {
#ifdef DEBUG_WRITE_HK_TO_DRMS
     printkerr("DEBUG:MESSAGE:check_hk_record_within_range:Received packets --within-- range(yyyy.mm.dd_hh):%d.%02d.%02d_%02d\n", y,m,d,h);
#endif
     return (1);
  }
  else if ( curr_packet_number == curr_time_range && h < curr_hr_range)
  {
     /* check if same day for curent packet and current time range then check hour is with current hour range */
#ifdef DEBUG_WRITE_HK_TO_DRMS
     printkerr("DEBUG:MESSAGE:check_hk_record_within_range:Received packets --within-- range(yyyy.mm.dd_hh):%d.%02d.%02d_%02d\n", y,m,d,h);
#endif
     return (1);
  }
  else
  {
#ifdef DEBUG_WRITE_HK_TO_DRMS
     printkerr("DEBUG:WARNING MESSAGE:check_hk_record_within_range:Received packets not within range(yyyy.mm.dd_hh):%d.%02d.%02d_%02d\n", y,m,d,h);
#endif
     return (0);
  }
 
}
