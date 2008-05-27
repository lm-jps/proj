/*****************************************************************************
 * Filename: write_hk_to_drms.c                                              *
 * Author: Carl                                                              *
 * Create Date: February, 2, 2008                                            *
 * Description: This file contains modules to write housekeeping keywords    *
 *              to DRMS.                                                     *
 * (C) Stanford University, 2008                                             *
 ****************************************************************************/
#define  DEBUG_WHTD  0
#include "jsoc_main.h"
#include "packets.h"
#include "decode_hk_vcdu.h" 
#include "write_hk_to_drms.h"
#include "printk.h"
/***********************   Protoype Functions  *******************************/
int   write_hk_to_drms(DRMS_Record_t *record, CCSDS_Packet_t **ccsds_pkt);
TIME  SDO_to_DRMS_time(int sdo_s, int sdo_ss);
char *get_packet_version_number( HK_Keyword_t *kw);
char *lookup_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t **jmap);
void load_map_data(int apid, JSOC_Version_Map_t  *jm);
int find_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t  *jm, char *dsn);
int   get_pkt_time_from_timecodes(HK_Keyword_t *hk, TIME *ptime);
int   check_for_apid(int apid, JSOC_Version_Map_t  *jm);
void  free_jsvn_map( JSOC_Version_Map_t  *top_jm);



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
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec;
  TIME pkt_time;
  char keyname[HK_LEV0_MAX_LONG_KW_NAME];
  char query[HK_LEV0_MAX_DSNAME_STR];
  char pkt_ver_num[HK_LEV0_MAX_PVN_STR];
  char project_name[HK_LEV0_MAX_PROJECT_NAME];
  char  datatype_name[HK_LEV0_MAX_DATATYPE_NAME];
  char *pvn, *dtname, *pjname;
  char *directory ;
  char *suffix_filename ;
  int  status;
  int rec_alreadycreated_flag;

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
  pjname=project_name;
  dtname=datatype_name;


  /* check if record is set like for case of Lev 0 data series*/
  if(record)
  {
    rec_alreadycreated_flag = 1;
  }
  else 
  {
    /*setting for writing to DRMS Level 0 by APIDi data series*/
    rec_alreadycreated_flag=0; 
  }

  /* check 3 environment variables are set */
  /* get project name */
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

  /* get data type name */
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
      strcpy(pvn , get_packet_version_number(ccsds->keywords));   

      /* lookup data series name and set in DRMS*/
      keytype= DRMS_TYPE_STRING;
      strcpy(keyname, "ISPSNAME");

      /*allocate memory for string value */
      key_anyval.string_val = (char *)malloc(sizeof(char) * 100);

      /* get isp series name */
      strcpy(key_anyval.string_val, lookup_data_series_name(ccsds,&jmap));
      status = drms_setkey(rec, keyname, keytype, &key_anyval);

      /* free memory */
      free (key_anyval.string_val);
    }
    else
    {
      /* get packet version number */
      strcpy(pvn , get_packet_version_number(ccsds->keywords));   

      /*  lookup data series name */
      strcpy(query , lookup_data_series_name(ccsds,&jmap)); 
      if (DEBUG_WHTD) printf("query is %s\n", query);

      /* create record in drms */
      rs = drms_create_records( drms_env, 1, query, DRMS_PERMANENT, &status);
      if (!rs)
      {
        printkerr("ERROR at %s, line %d: Could create record for data series "
                  "<%s>. Because the data series needs to be created. "
                  "Skipping write to drms of this data series.\n",
                  __FILE__ , __LINE__ , query);
       ccsds= ccsds->next;
       continue;
      }
      else
      {
        rec = rs->records[0];
      }
    }

    /* get and set packet_time in record*/
    kw = ccsds->keywords;
    if (!get_pkt_time_from_timecodes(kw, ptime))
    {
      /*did not set PACKET_TIME */
      if (rec_alreadycreated_flag)
      {
         /* don't control record closing -skip*/
      }
      else
      {
         status = drms_close_records(rs, DRMS_INSERT_RECORD);
      }
      /* assume need to have packet time set for both cases so do this */
      printkerr("Warning at %s, line %d: Could not set PACKET_TIME keyword. "
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
                                                                                 
    /* set packet version number(xxx.yyy) keyword */
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


#ifdef DEBUG_WRITE_HK_TO_DRMS
      /* check got good return status */
      if (status)
      {
        /* compile in by adding -DDEBUG_WRITE_HK_TO_DRMS in makefile */
        printkerr("Carl testing:Warning at %s, line %d: Cannot setkey in drms record."
                  " For keyword <%s>\n",
                 __FILE__,__LINE__, keyname);
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
      status = drms_close_records(rs, DRMS_INSERT_RECORD);
      if (status)
      {
        printkerr("ERROR at %s, line %d: Cannot close drms record.\n",
                   __FILE__ , __LINE__ );
      }
    }
    /* don't expect more than one node to rec_alreadycreated_flag case,
       so don;t need to go to next CCSDS_Packet_t node. If do get
       more than one node a DRMS_Record_Set_t need to be used as argument  */
    ccsds= ccsds->next;
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
     if( strstr(kw->name, "VER_NUM"))
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
 * Name:        LOOKUP DATA SERIES NAME                                      *
 * Function :   char *lookup_data_series_name(CCSDS_Packet_t *ccsds_ptr,     *
 *                                            JSOC_Version_Map_t **jmap_ptr) *
 * Description: lookup the data series name based on apid and packet         * 
 *              packet version number. Returns string value                  *
 *****************************************************************************/
char *lookup_data_series_name(CCSDS_Packet_t *ccsds_ptr, JSOC_Version_Map_t **jmap_ptr) 
{
   /* declarations */
   JSOC_Version_Map_t *jm, *tmp_jm, *last_jm, *prev_jm;
   char data_series_name[HK_LEV0_MAX_DSNAME_STR];
   char *dsn;
   int  apid;

   /* initialized variables */
   dsn = data_series_name;
   jm = *jmap_ptr;

   /* get apid */
   apid = ccsds_ptr->apid;

   /*check if structure has node with data series name */
   if (!jm)
   {
      /* if empty structure then create top node */
      assert(jm = (JSOC_Version_Map_t *)malloc(sizeof(JSOC_Version_Map_t))); 
      *jmap_ptr=jm;
      jm->apid=(short)apid;
      jm->next=NULL;

      /* load Map_Data nodes and create data series names 
         for each packet version number*/
      load_map_data(apid, jm);
  
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
           assert(jm = (JSOC_Version_Map_t *)malloc(sizeof(JSOC_Version_Map_t)));
           *jmap_ptr=jm;
           jm->apid=(short)apid;
           jm->next=NULL;
                                                                                                
           /* load Map_Data nodes and create data series names for each 
              packet version numberi and return data series name*/
           load_map_data(apid, jm);
                                                                                                
           /* find data series name in Map_Data nodes based on packet version number*/
           (void)find_data_series_name(ccsds_ptr, jm, dsn);
       } 
       else 
       { 
         /* if cannot find node for apid then create one */ 
         /* go to end of JSOC Version Map nodes and add node */
         /* first create JSOC Version Map node */
         assert(tmp_jm = (JSOC_Version_Map_t *)malloc(sizeof(JSOC_Version_Map_t))); 
         tmp_jm->apid= (short)apid;
         tmp_jm->next= NULL;

         /* link in list of map data to JSOC Version Map node */
         load_map_data(apid, tmp_jm);     

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
 * Function :   void load_map_data(int apid, JSOC_Version_Map_t  *jm)       *
 * Description: loads PVN-TO-JSVN-<apid> files data                          * 
 *              in structures. This structure will                           * 
 *              be used to lookup data series name                           * 
 *****************************************************************************/
void load_map_data(int apid, JSOC_Version_Map_t  *jm)
{
  /*declarations */
  FILE *file_ptr;
  Map_Data_t  *tmp_dm, *prev_dm;
  char directory_filename[HK_LEV0_MAX_FILE_NAME];
  char fn[HK_LEV0_MAX_FILE_NAME];
  char sn[HK_LEV0_MAX_FILE_NAME];
  char *didn;
  char *directory ;
  char *suffix_filename;
  char line[HK_LEV0_MAXLINE_IN_FILE];
  char *pn;
  char *pn1,*pn2,*pn3;
  int  wn, dn;

  /* initialized variables */
  char *filename=fn;

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

  /* make filename for looking up jsoc version number */
  sprintf(filename,"%4.4d%s",apid,suffix_filename);
  sprintf(directory_filename, "%s/%s", directory, filename);

  /* get apid  and set to JSOC Version Map node*/
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
        assert(tmp_dm = (Map_Data_t*)malloc(sizeof(Map_Data_t))); 
        tmp_dm->next= NULL;
        jm->mdata = tmp_dm;
        prev_dm =tmp_dm;
      }
      else
      {
        assert(tmp_dm = (Map_Data_t*)malloc(sizeof(Map_Data_t))); 
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
        /* else set data series name field*/
        sprintf(tmp_dm->dsn, "%s.%s_%04d_%s", pn, didn, apid, tmp_dm->jvn);
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
  return 1;
}


                                                                                               
/*************************************************************************
 * GET PACKET TIME FROM TIME CODES                                       *
 * FUNCTION: int get_pkt_time_from_timecodes(HK_Keyword_t *,char *,TIME) *
 * DESCRIPTION: Gets Packet time based on Time Codes.                    *
 *************************************************************************/
int get_pkt_time_from_timecodes(HK_Keyword_t *hk,  TIME *ptime)
{
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
                                                                                                
  if (!SEC_FOUND_FLAG)
  {
    printkerr("Error at %s, line %d: Did not find TIMECODE_SECONDS value for "
              "calculating the PACKET_TIME keyword and index. Returning error "
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
 * NAME: Check for Apid                                                  *
 * FUNCTION: int check_for_apid(int apid, JSOC_Version_Map_t  *jm )      *
 * DESCRIPTION: Check for apid  in structure                             *
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
