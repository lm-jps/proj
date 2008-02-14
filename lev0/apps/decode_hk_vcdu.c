#include <stdio.h>
#include <assert.h>
#include "decode_hk.h"
#include "packets.h"
#include "decode_hk_vcdu.h"
#include "save_packet_to_dayfile.h"
#include <sys/time.h>
#include <time.h>
#include "timeio.h"
#include <dirent.h>
#include "printk.h"

/*********** static function prototypes **************/
static unsigned short *decode_im_pdu(unsigned short *w, IM_PDU_Packet_t *p);
static unsigned short *decode_ccsds(unsigned short *w, CCSDS_Packet_t *p);

/*************  function prototypes ******************/
TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss);
int lookup_fsn(CCSDS_Packet_t **ptr, unsigned int *Fsn);
int save_packet_to_dayfile(unsigned short *word_ptr, int apid, HK_Dayfile_Data_t **df_head); 
int write_packet_to_dayfile(HK_Dayfile_Data_t **df_head); 
int filter_out_non_isps(CCSDS_Packet_t **ptr);
int write_hk_to_drms(DRMS_Record_t *record, CCSDS_Packet_t **ccsds_pkt);
int get_status( int lev0_status[], int status_count);
 
/*********** extern function prototypes **************/
extern int free_dayfile_data( HK_Dayfile_Data_t **df_head); 
extern int free_dayfile_pkt_data( HK_Dayfile_Data_t **df_head); 
/********** High level hk telemetry routines **************/

/* Receive the next VCDU in the telemetry stream. The following
   action is taken depending on the contents of the VCDU.
    
   A) If the vcdu contains a single science data packet
      it skipped

   B) If the vcdu contains one or more housekeeping packets each
      one is APID value checked and based on APID the following
      three cases of processing is carried out

      a) decode hk keyword in structure and pass back structure

      b) skip decode of hk kywords and pass back null hk structure
         and save packet to dayfile (i.e., OBT and Sequencer packets)
   
      c) decode hk keyword in structure and pass back structure and
         also save packet to dayfile(i.e., ISP, etc)

      d) skip decode of hk keywords and skip save of packet to dayfile 
         for APID not requiring processing. 

    Overall Status value:

    0: HK_SUCCESS_HK_ALL:            All hk packets where decoded successfully 
                                     or saved to dayfiles successfully with no 
                                     errors.  Successfully skipped processing of 
                                     hk packets that are not on list to be 
                                     processed. Successfully skipped if found 
                                     Image Packets. When this occurs a CCSDS 
                                     link list is probably returned. So probably 
                                     need to check structure and if available 
                                     write HK Structure (within CCSDS structure) 
                                     of keyword names and values to DRMS.
     
    1: HK_SUCCESS_HK_SOME:           Some errors found when doing HK Processing
                                     but did successfully decode and return a 
                                     structure or successfully saved dayfiles.
                                     When this occurs a CCSDS link list is
                                     probably returned. So probably need to check
                                     structure and if available write HK Structure 
                                     (within CCSDS structure) of keyword names and 
                                     values to DRMS. 

    4: SUCCESS_SKIP_IMAGE:           Detected Image and successfully skipped processing.

    5: SUCCESS_SKIP_PROCESSING_APID: Detected ALL APIDS in VCDU that are not 
                                     normally processed. Successfully skipped 
                                     processing and returned this value.

    < 0: See Macro definitions:      All errors found when doing HK Processing for
                                     a single VCDU. The last error code occurrance 
                                     is returned, since cannot return all error 
                                     codes for every hk packet in VCDU. Consult 
                                     log file for full display of all errors that 
                                     occurred.  Consult decode_hk_vcdu.h for macro 
                                     definitions of negative error codes. The
                                     returned CCSDS structure has no items to
                                     write to DRMS when this occurs.

    Returned Values to Level 0 Top Module: See .h file
*/

/*************************  decode next hk vcdu  ****************************/
int decode_next_hk_vcdu(unsigned short vcdu[PACKETWORDS],  CCSDS_Packet_t **hk_packets, unsigned int *Fsn)
{
  /* declarations */
  CCSDS_Packet_t ccsds, *p, *p_hk;
  IM_PDU_Packet_t im_pdu;
  int spfd;
  int i, j, k; 
  unsigned char *p1,*p2;
  static unsigned short buffer[PACKETWORDS];
  unsigned short *w, *hkstart;
  int lz_status[100];
  int hk_status = ERROR_NODATA;
  int wdf_status;
  int wd_status;
  int foni_status;
  int decode_status;
  int lev0_status;
  int fsn_status=1;
  DRMS_Record_t *record=NULL;
  static short writeFlag= HK_INIT_WRITE_FLAG;
  static HK_Dayfile_Data_t *dfd=NULL;

  /* init variables */
  j=0;
  lev0_status=ERROR_NODATA; /* -13 */

  /* init lz_status array */
  for (k=0;k<100;k++) lz_status[k]=99;

  /****** 
  Copy the raw telem data to a buffer. The incoming telemetry stream is in 
  big endian format. Swap the byte order if we are doing the processing on 
  a little endian machine. 
  ******/
#if __BYTE_ORDER == __LITTLE_ENDIAN
  p1 = (unsigned char *)vcdu;
  p2 = (unsigned char *)buffer;
  for (i=0; i<PACKETWORDS; i++)
  {
    *(p2+1) = *p1;
    *p2 = *(p1+1);
    p1 += 2;
    p2 += 2;
  }
#else  
  memcpy(buffer, vcdu, PACKETBYTES);
#endif

  /* Decode IM_PDU headers. */
  w = decode_im_pdu(buffer, &im_pdu);

  /* Loop over all CCSDS packet in the packet zone. */
  *hk_packets = p_hk = NULL;
  do
  {
    /* Decode the CCSDS header. */
    hkstart = vcdu + (w-buffer);
    w = decode_ccsds(w, &ccsds);
    
    /* Branch depending on the APID. */
    switch(ccsds.apid)
    {  
    /*************************************/
    /***  Case of Finding End of VCDU  ***/
    /*************************************/
    case 0:
      /* We have reached the end of the VCDU */      
      if (hk_status != ERROR_NODATA) 
      {  
        hk_status=HK_SUCCESS_REACHED_END_VCDU;
      }
      lz_status[j++] = hk_status;
      break;

    /******************************************/
    /*** Case of Finding Image Packet-skip  ***/
    /******************************************/
    case APID_HMI_SCIENCE_1:
    case APID_HMI_SCIENCE_2:
    case APID_AIA_SCIENCE_1:
    case APID_AIA_SCIENCE_2:
      hk_status=HK_SUCCESS_SKIP_IMAGE;
      lz_status[j++] = hk_status;
      break;
      
    /*********************************************/
    /*** Case of Finding HK Time packet-skip  ***/
    /********************************************/
    case APID_HMI_TIME_1:
    case APID_HMI_TIME_2:
    case APID_AIA_TIME_1:
    case APID_AIA_TIME_2:
      /* This is an empty timestamp packet.  Do nothing. */
      hk_status=HK_SUCCESS_HKTIME;
      lz_status[j++] = hk_status;
      break;

    /***********************************************************************/
    /*** Case of Decoding to CCSDS_Packet_t struct and Writing Day File  ***/
    /***********************************************************************/
    case APID_HMI_IMSTAT_1:
    case APID_HMI_IMSTAT_2:
    case APID_AIA_IMSTAT_1:
    case APID_AIA_IMSTAT_2:

      /* This is an image status, or maybe sequencer or maybe obt packet. */
      decode_status = decode_hk_keywords(hkstart, ccsds.apid, &ccsds.keywords);

      /* set status based on decode status */
      if (decode_status != 0) 
      { 
        /* if did have error or warning use value in decode_status*/
        hk_status=decode_status;
      }
      else  
      {
        /* if passed then map to success values */
        hk_status=HK_SUCCESS_DECODING;
      }
      lz_status[j++] = hk_status;

      /* create ccsds note and add HK_Keyword_t structure to CCSDS node */
      if (hk_status==HK_SUCCESS_DECODING)
      {      
	/* Allocate a CCSDS packet. */
	assert(p = malloc(sizeof(CCSDS_Packet_t)));
	memcpy(p, &ccsds, sizeof(CCSDS_Packet_t));

	/* Append to output list. */
	if (*hk_packets == NULL)
        {
	  *hk_packets = p_hk = p;
        }
	else
	{
	  p_hk->next = p;
	  p_hk = p_hk->next;
	}
	p_hk->next = NULL;
      }

      /* write packet to dayfile */
      spfd = save_packet_to_dayfile(hkstart ,ccsds.apid, &dfd) ;
      hk_status=spfd;
      if (hk_status == ERROR_HK_ENVIRONMENT_VARS_NOT_SET) 
      {
         /* never got a chance to save  packets */
         lev0_status= ERROR_HK_ENVIRONMENT_VARS_NOT_SET;
         lz_status[j++] = hk_status;
         return (lev0_status);
       }

     // take out if SOME means only decoded some and 
     // does  not mean  decoded and wrote-dayfile some 
      lz_status[j++] = hk_status;
      break;
    /*************************************/
    /*** Case of Writing Day File Only ***/
    /*************************************/
    case APID_HMI_SEQ_1:
    case APID_HMI_SEQ_2:
    case APID_HMI_OBT_1:
    case APID_HMI_OBT_2:
    case APID_AIA_SEQ_1:
    case APID_AIA_SEQ_2:
    case APID_AIA_OBT_1:
    case APID_AIA_OBT_2:
      if (hk_status == HK_SUCCESS_HKTIME   || 
          hk_status == HK_SUCCESS_DECODING || 
          hk_status == HK_SUCCESS_WRITE_DAYFILE)
      {      

        /* write packet to dayfile */
        spfd = save_packet_to_dayfile(hkstart ,ccsds.apid, &dfd) ;
        hk_status=spfd;
        if (hk_status == ERROR_HK_ENVIRONMENT_VARS_NOT_SET) 
        {
          /* never got a chance to save packets */
          lev0_status= ERROR_HK_ENVIRONMENT_VARS_NOT_SET;
          lz_status[j++] = hk_status;
          return (lev0_status);
        }
      }

      // take out if SOME means only decoded some and 
      // does  not mean  decoded and wrote-dayfile some 
      lz_status[j++] = hk_status;
      break;      

    /*********************************/
    /*** All Other Non-Valid cases ***/
    /*********************************/
    default:
      /* This is not known HK PACKET to process */
      hk_status=HK_SUCCESS_SKIP_PROCESSING_APID;
      lz_status[j++] = hk_status;
      break;
    }// switch
    w += ccsds.length/2;
  } while ( ccsds.apid>0 && (int)(w-buffer) < PACKETWORDS);

  /* finished VCDU, now write saved dayfile packets to filesystem 
     and free dayfile structure. Currentily write out dayfile after each VCDU.
     Also keeps using initial packet times filename to write packets too.*/
  if(dfd &&  hk_status != HK_SUCCESS_SKIP_IMAGE )
  {
    if (writeFlag == HK_WRITE_AFTER_VCDU_COUNT)
    {
      wdf_status = write_packet_to_dayfile( &dfd); 
      if (HK_SUCCESS_WRITE_DAYFILE != wdf_status ) 
      {
        printkerr("ERROR at %s, line %d: Return status failed for "
                  "write_packet_to_dayfile. Probably can't open dayfile.\n", 
                   __FILE__, __LINE__);
        return (ERROR_HK_FAILED_OPEN_DAYFILE);
      }
      wdf_status=free_dayfile_pkt_data( &dfd); 
      if (HK_SUCCESS_WRITE_DAYFILE != wdf_status) 
      {
        printkerr("Warning at %s, line %d: Return status failed for "
                  "free_dayfile_pkt_data\n", __FILE__, __LINE__);
      }
      writeFlag= HK_INIT_WRITE_FLAG;
    }
    else writeFlag++;
  }
  if(*hk_packets &&  hk_status != HK_SUCCESS_SKIP_IMAGE) 
  {
      /* write Keywords to DRMS in Level 0 Data Series by APID */
      wd_status = write_hk_to_drms(record, hk_packets);

      if (wd_status == HK_SUCCESS_WROTE_TO_DRMS) 
      {
        /* wrote successfully to drms -now free CCSDS_Packet_t nodes */
        lev0_status= SUCCESS_HK_NEED_TO_CTD; /* set to 0 */
      }
      else if (wd_status == ERROR_HK_ENVIRONMENT_VARS_NOT_SET) 
      {
        /* never got a chance to write to DRMS */
        lev0_status= ERROR_HK_ENVIRONMENT_VARS_NOT_SET;
      }
      else if  (wd_status == ERROR_HK_FAILED_TO_FIND_TIMECODES)
      {
        /* never got a chance to write to DRMS */
        lev0_status= ERROR_HK_FAILED_TO_FIND_TIMECODES;
      }
      else if (wd_status == ERROR_HK_FAILED_CLOSE_DRMS_RECORD)
      {
        lev0_status= ERROR_HK_FAILED_CLOSE_DRMS_RECORD;
      }
      else if (wd_status == ERROR_HK_FAILED_OPEN_DRMS_RECORD)
      {
        lev0_status= ERROR_HK_FAILED_OPEN_DRMS_RECORD;
      }
      else 
      {
        printkerr("Warning at %s, line %d: Return status unexpected "
                  "for write_hk_to_drms\n", __FILE__, __LINE__);
        
      }
  }/* if  time to write hk packets to DRMS */

  if (get_status( lz_status, j) != 4) 
  { 
    /*printf("OVERALL Status on saving packets,decoding packets  <%d> and writing to drms lev0 hk series <%d>\n", 
            get_status( lz_status, j), lev0_status);
    */
  }
 
  /* determine structure to pass back to LEV0 top module -note only want ISPs */
  /* first should not free ccsds nodes if ISPs there */
  /* Filter out non-isp nodes in CCSDS_Packet_t link list */
  /* if there is a isp node return back success values */
  /* if there are no isp nodes return back status meaning no isp's found */
  if(hk_packets &&  hk_status != HK_SUCCESS_SKIP_IMAGE  && lev0_status >= ERROR_NODATA) 
  {
    /* filters out non-isps and free node that are not isp's */
    foni_status = filter_out_non_isps(hk_packets);  
    if(foni_status == 1)
    {
      /*got a ISP node in CCSDS_Packet_t link list so set to value for
        lev0 top module to know to write to level0 data series*/
      /* lookup FSN value to pass back */
      fsn_status= lookup_fsn(hk_packets, Fsn);
      if(fsn_status)
      {
         lev0_status = ERROR_HK_FAILED_GETTING_FSN; /* set to -27 */
      }
      else
      {
         /* successfully got FSN value */
         lev0_status = SUCCESS_HK_NEED_TO_WTD_CTD; /* set to 1 */
      }
    }
  }

  /* Set status returned for lev0 top module by evaluation of  */
  /* errors when writing to dayfile, decoding hk, etc.        */
  return (lev0_status);
}


/*******************   Packet header routines      ********************/

/* Extract fields from IM_PDU packet header. */
static unsigned short *decode_im_pdu(unsigned short *w, IM_PDU_Packet_t *p)
{
    /* Word 0 */
    p->im_pdu_count = ((long long) (*w & 0x3ff)) << 32;
    p->im_pdu_id = (unsigned char)(*w++ >> 10);
    /* Word 1 */
    p->im_pdu_count |= *w++ << 16;
    /* Word 2 */
    p->im_pdu_count |= *w++;
    /* Word 3 */  
    p->fhp = (unsigned short)((*w++) & 0x7ff);
    return w;
}

/************     Extract fields from CCSDS packet header    **************/
static unsigned short *decode_ccsds(unsigned short *w, CCSDS_Packet_t *p)
{
    /* Word 4 */  
    p->version = (unsigned char)((*w >> 13) & 0x7);
    p->type = (unsigned char)((*w >> 12) & 1);
    p->shf =  (unsigned char)((*w >> 11) & 1);
    p->apid = (unsigned short)((*w++) & 0x7ff);
    /* Word 5 */  
    p->sf = (unsigned char)((*w >> 14) & 0x3);
    p->ssc = (unsigned short)((*w++) & 0x3fff);
    /* Word 6 */  
    /* According to the CCSDS spec the packet length field contains 
       "total number of octets - header octets - 1". 
       Add 1 to get length of data field in bytes. */
    p->length = (unsigned short)(*w++ + 1);
    return w;
}

/********************** hk_ccsds_free function ******************/
void  hk_ccsds_free(CCSDS_Packet_t **ptr)
{
  CCSDS_Packet_t *tmp, *p;
  p= *ptr;
  while(p)
  {
    if (p->keywords)
    {
      deallocate_hk_keyword(p->keywords);
    }
    tmp = p->next;
    free(p);
    p = tmp;
  }
}


/********************** filter out non isps function ******************/
int filter_out_non_isps(CCSDS_Packet_t **ptr)
{
  /* declarations */
  CCSDS_Packet_t *tmp, *p, *prev;
  int isp_count=0;

  /* initialized variables */
  p= *ptr;
  prev=p;

  /* locate non-isps and free node*/
  while(p)
  {
    
    if (p->apid != APID_HMI_IMSTAT_1 && p->apid != APID_HMI_IMSTAT_2 &&
        p->apid != APID_AIA_IMSTAT_1 && p->apid != APID_AIA_IMSTAT_2)
    {  
      if (p->keywords)
      {
        deallocate_hk_keyword(p->keywords);
      }
      prev->next=p->next;
      tmp = p->next;
      free(p);
      p = tmp;
    }
    else
    {  /* if is isp */
       isp_count++;
       if(isp_count == 1)
       {
         /*set new head once*/
         *ptr=p;
       }
       prev=p;
       p=p->next;
    }
  }
  return (isp_count);
}


/********************** lookup FSN function ******************/
int lookup_fsn(CCSDS_Packet_t **ptr, unsigned int *Fsn)
{
  /* declarations */
  CCSDS_Packet_t *p;
  HK_Keyword_t *kw;
  unsigned int fsn=0;
  int status=1;  //not found case

  /* initialized variables */
  p= *ptr;
  kw=p->keywords;

  /* locate Fsn */
  while(kw)
  {
    
    if (strstr(kw->name, "FILTERGRAM_SN") || strstr(kw->name, "FRAME_SN"))
    {  
       fsn=kw->eng_value.uint32_val;
       status = 0;
       *Fsn=fsn;
       break;
    }
    else
    {  
       kw=kw->next;
    }
  }
  return (status);
}



/***********************   get status ************************************/
int get_status( int lz[], int jj)
{
  /* declarations */
  int k, count_errors, last_error, skip_apids;

  /* Set status returned to lev0 top module by evaluation of  */
  /* errors when writing to dayfile, decoding hk, etc.        */
  for (k=0, count_errors=0,skip_apids=0, last_error=0; k < jj; k++)
  {
     /* count errors  and apids skipped*/
     if (lz[k] < 0 )
     {
        count_errors+=1;
        last_error= lz[k];
     }
     else if (lz[k] ==  SUCCESS_SKIP_PROCESSING_APID)
     {
       skip_apids+=1;
     }
  }
 
  /* if no errors return status HK_SUCCESS_SKIP_IMAGE, 
     HK_SUCCESS_SKIP_PROCESSING_APID, or HK_SUCCESS_HK_ALL */ 
  if (!count_errors) 
  {
    if (lz[0] == HK_SUCCESS_SKIP_IMAGE)  return HK_SUCCESS_SKIP_IMAGE;
    else if (skip_apids == (jj - 2 ))
    {
      return HK_SUCCESS_SKIP_PROCESSING_APID;
    }
    else  return HK_SUCCESS_HK_ALL;
  }
  /* if count errors equal to total process then return 
     HK error-status of last item processed */
  else if ( count_errors == jj - 2)
  {
    return last_error;
  }
  /* if got some errors and some sucessful processing return HK_SUCCESS_HK_SOME */
  else if ( count_errors < jj - 2)
  {
     return HK_SUCCESS_HK_SOME;
  }
  else return lz[0];
}

/********************************************************/
/* SDO_to_DRMS_time --Need to hook up to library to get */
/********************************************************/
TIME SDO_to_DRMS_time(int sdo_s, int sdo_ss) 
{
/*changes done:
1.changed args from float to int 
2.added line below...int ss.. 
3.changed return statement with ss parameter!
*/
static TIME sdo_epoch;
int ss=(sdo_ss >> 16) & 0xFFFF;
static int firstcall = 1;
if (firstcall)
  { /* time_1958 - time_1977_TAI, to be added to SDO time to get DRMS time */
  firstcall = 0;
  sdo_epoch = sscan_time("1958.01.01_00:00:00_TAI");
  } 
return(sdo_epoch + (TIME)sdo_s + (TIME)ss/65536.0);
}  
