#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/lev0/apps/get_image_location.c,v 1.7 2011/04/08 19:51:34 carl Exp $" 
/* GET IMAGE LOCATION to be merged into Jim's Lev1 code */
/* NOTE1:Jim:Jim's code needs to free Image_Location after calling get_location_information. */
/* NOTE2:Jim:Example main file used to test at:/home3/carl/cvs/JSOC/proj/lev1/apps/gif_main.c*/
/******************** defines ***********************************************/
/* use by Jim's pre-productions test-> #define GMP_MASTER_POINTING_SERIES  "su_carl.master_pointing" */
/* use by Carl to test ->#define GMP_MASTER_POINTING_SERIES  "su_carl.test99_master_pointing"*/
/* use for production after launch ->#define GMP_MASTER_POINTING_SERIES  "sdo._master_pointing"*/
#define GMP_MASTER_POINTING_SERIES  "sdo.master_pointing"
#define GMP_DRMS_OPEN_FAILED           1
#define GMP_MAX_DSNAME_STR           100
#define GMP_MAX_KEYWORD_NAME_STR     100
#define GMP_MAX_QUERY_STR            200
#define GMP_MAX_MPO_REC_SIZE         200
#define GMP_MAX_TELESCOPE_STR         10
#define GMP_MALLOC_FAILED              2
#define GMP_PACKET_TIME_STR           50
#define GMP_PASSED_STATUS              0
#define GMP_TOBS_NOT_SORTED_FAILED     3

/******************** structures ********************************************/
typedef struct Image_Location_struct {
   //*******
   // Inputs
   //*******
   // Obs time
   TIME tobs;

   // For AIA:1,2,3,4. For HMI:1 or 2
   int camera;

   // For AIA:SDO/AIA. For HMI:SDO/HMI
   char telescope[GMP_MAX_TELESCOPE_STR]; 

   //wavelength
   int wavelength;

   //*******************************************
   // Returned values from get_image_location()
   //*******************************************
   float x;
   float y;
   float instrot;
   float imscale;
   float yinrtb;              //value in mp is SC_Y_INRT_BIAS;
   float zinrtb;              // value in mp is SC_Z_INRT_BIAS;
   // MPO_REC,  "Master Pointing series record pointer" 
   char mpo_rec[GMP_MAX_MPO_REC_SIZE];

} Image_Location;

typedef struct TIME_MP_struct {
   TIME tstart;
   TIME tstop;
} TIME_MP;

/******************** functions *********************************************/
int find_mp_record(TIME_MP *time_mp, int nrec, TIME tobs);

/****************************************************************************/
/*************  get_image_location function *********************************/
/****************************************************************************/
int get_image_location(DRMS_Env_t *drms_env, int ncnt, Image_Location **ptr_imageloc)
{

  /* drms record create variables */
  DRMS_RecordSet_t *rs;
  DRMS_RecordSet_t *rset;
  DRMS_Record_t *rec;
  /* local variables */
  Image_Location *tptr;
  TIME tstart_range,tend_range;
  TIME_MP *time_mp;
  char dsname[GMP_MAX_DSNAME_STR];
  char nquery[GMP_MAX_QUERY_STR];
  int i,j,idx;
  int nrec;
  int status;

  /* initialize variables */
  tptr= *ptr_imageloc;
  
  /*get master pointing series name */
  strcpy(dsname, GMP_MASTER_POINTING_SERIES);

  /* set tstart and  tend range */
  tstart_range= (tptr)->tobs;
  tend_range= (tptr+(ncnt-1))->tobs ;

  /* check for not sort list of tobs */
  if (tstart_range > tend_range)
  {
    printkerr("ERROR at %s, line %d: Failed because TOBS list is not sorted. "
              "Returning error. Not setting return values.\n",__FILE__,__LINE__);
    return (GMP_TOBS_NOT_SORTED_FAILED);
  }

  /* create query statement */
  sprintf(nquery,"%s[? T_START < %.0f AND T_STOP > %.0f ?]",dsname,tend_range,tstart_range);
  //printf("test:get_location_information:<%s>\n",nquery);

  /* open records with query */
  rset = drms_open_records(drms_env, nquery, &status);
  if(!rset || !rset->n || status)
  {
     /* return error PT_DRMS_OPEN_FAILED if cannot open */
     printkerr("ERROR at %s, line %d: Failed to open  master pointing series:<%s>. "
               "Check envirionment  GMP_MASTER_POINTING_SERIES variable is set "
               "correctly. Or check time range in master pointing series is available.\n",
               __FILE__,__LINE__, dsname);
     return (GMP_DRMS_OPEN_FAILED);
  }
  nrec= rset->n; 

  /*allocate space for time records found in query */
  time_mp= (TIME_MP *) malloc(nrec*sizeof(TIME_MP)); 
  if(!time_mp) return (GMP_MALLOC_FAILED);

  /*put time records in allocated time records */
  for(i=0; i<nrec;i++)
  {
    (time_mp+i)->tstart=drms_getkey_time(rset->records[i],"T_START",&status);
    (time_mp+i)->tstop=drms_getkey_time(rset->records[i],"T_STOP",&status);
  }

  /* loop thru array of structures and fill in values including tobs value */
  for(i=0; i < ncnt; i++)
  {
    /* find_mp_record where tobs values is between T_START and T_STOP values of record */
    idx=find_mp_record(time_mp,  nrec, (tptr+i)->tobs);
    if (idx == -1)
    {
      printkerr("WARNING at %s, line %d: Setting return values to DRMS_MISSING_FLOAT "
                "because could not find record. TOBS time passed "
                "is <%f>.\n", __FILE__,__LINE__,(tptr+i)->tobs);
      (tptr+i)->x= DRMS_MISSING_FLOAT;
      (tptr+i)->y= DRMS_MISSING_FLOAT;
      (tptr+i)->instrot=DRMS_MISSING_FLOAT;
      (tptr+i)->imscale=DRMS_MISSING_FLOAT;
      (tptr+i)->zinrtb=DRMS_MISSING_FLOAT;
      (tptr+i)->yinrtb=DRMS_MISSING_FLOAT;
      continue;
    }

    /* get record found where found tobs value is between T_START and T_STOP values in mp */
    rec=rset->records[idx];

    /* set spacecraft z inrt and y inrt bias value */
    (tptr+i)->zinrtb =drms_getkey_float(rec, "SC_Z_INRT_BIAS", &status);
    if(status == -10006) 
    {
      (tptr+i)->zinrtb=DRMS_MISSING_FLOAT;
      printkerr("WARNING at %s, line %d: Setting zinrtb failed because of error code DRMS_ERROR_UNKNOWNKEYWORD. "
                "Status returned from drms_getkey_float was %d\n",  __FILE__,__LINE__, status);
    }
    
    (tptr+i)->yinrtb=drms_getkey_float(rec, "SC_Y_INRT_BIAS", &status);
    if(status == -10006) 
    {
      (tptr+i)->yinrtb=DRMS_MISSING_FLOAT;
      printkerr("WARNING at %s, line %d: Setting yinrtb failed because of error code DRMS_ERROR_UNKNOWNKEYWORD. "
                "Status returned from drms_getkey_float was %d\n",__FILE__,__LINE__, status);
    }


    /* load parameters in structure using values from this mp record */
    if(!strcmp((tptr+i)->telescope, "SDO/HMI"))
    {
        if((tptr+i)->camera == 1)
        {
           (tptr+i)->x =drms_getkey_float(rec, "H_CAM1_X0", &status);
           (tptr+i)->y=drms_getkey_float(rec, "H_CAM1_Y0", &status);
           (tptr+i)->instrot=drms_getkey_float(rec, "H_CAM1_INSTROT", &status);
           (tptr+i)->imscale=drms_getkey_float(rec, "H_CAM1_IMSCALE", &status);
           sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
        }
        else if((tptr+i)->camera == 2)
        {
           (tptr+i)->x=drms_getkey_float(rec, "H_CAM2_X0", &status);
           (tptr+i)->y=drms_getkey_float(rec, "H_CAM2_Y0", &status);
           (tptr+i)->instrot=drms_getkey_float(rec, "H_CAM2_INSTROT", &status);
           (tptr+i)->imscale=drms_getkey_float(rec, "H_CAM2_IMSCALE", &status);
           sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
        }
        else
        {
           printkerr("WARNING at %s, line %d: Setting returns values set to DRMS_MISSING_FLOAT "
                     "because could not find camera value for hmi. Camera value passed "
                     "is <%d>.\n", __FILE__,__LINE__,(tptr+i)->camera);
           (tptr+i)->x= DRMS_MISSING_FLOAT;
           (tptr+i)->y=DRMS_MISSING_FLOAT;
           (tptr+i)->instrot=DRMS_MISSING_FLOAT;
           (tptr+i)->imscale=DRMS_MISSING_FLOAT;
        }
    }/*end if HMI data to retrieve and set */
    else if(!strcmp((tptr+i)->telescope, "SDO/AIA"))
    {
      if ( (tptr+i)->wavelength == 94)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_094_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_094_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_094_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_094_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 131)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_131_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_131_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_131_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_131_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 171)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_171_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_171_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_171_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_171_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      //added after ok from jps to treat 193 and 195 identical
      else if (((tptr+i)->wavelength == 193) || ((tptr+i)->wavelength == 195))
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_193_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_193_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_193_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_193_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 211)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_211_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_211_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_211_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_211_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 304)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_304_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_304_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_304_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_304_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 335)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_335_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_335_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_335_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_335_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 1600)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_1600_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_1600_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_1600_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_1600_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 1700)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_1700_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_1700_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_1700_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_1700_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else if ( (tptr+i)->wavelength == 4500)
      {
        (tptr+i)->x=drms_getkey_float(rec, "A_4500_X0", &status);
        (tptr+i)->y=drms_getkey_float(rec, "A_4500_Y0", &status);
        (tptr+i)->instrot=drms_getkey_float(rec, "A_4500_INSTROT", &status);
        (tptr+i)->imscale=drms_getkey_float(rec, "A_4500_IMSCALE", &status);
        sprintf((tptr+i)->mpo_rec,  "%s[:#%lld]", GMP_MASTER_POINTING_SERIES, rec->recnum);
      }
      else 
      {
        printkerr("WARNING at %s, line %d: Setting return values to DRMS_MISSING_FLOAT "
               "because could not find wavelength value for aia. Wavelength value passed "
               "is <%d>.\n", __FILE__,__LINE__,(tptr+i)->wavelength);
        (tptr+i)->x= DRMS_MISSING_FLOAT;
        (tptr+i)->y= DRMS_MISSING_FLOAT;
        (tptr+i)->instrot=DRMS_MISSING_FLOAT;
        (tptr+i)->imscale=DRMS_MISSING_FLOAT;
      }
    }/*end else if AIA data to retrieve and set */
  }/*for loop thru tobs values passed */

  /* close drms record */
  drms_close_records(rset, DRMS_FREE_RECORD);

  /*free MP_TIME structure */
  free(time_mp);
  return(GMP_PASSED_STATUS);
}

/****************************************************************************/
/*************  find_mp_record function *************************************/
/****************************************************************************/
int find_mp_record(TIME_MP *time_mp, int nrec, TIME tobs)
{
  int i;
  int index=-1;
  
  /* check if tobs time + or - range factor is between master pointing record T_START/T_STOP */ 
  for(i=0; i<nrec; i++)
  {
    if( (fabs((time_mp + i)->tstart) < fabs(tobs) ) &&
        (fabs((time_mp + i)->tstop)  > fabs(tobs) ) )
    {
       index=i;
       return index;
    }
  }
  /* never found match- return index= -1 */
  return index;
}

