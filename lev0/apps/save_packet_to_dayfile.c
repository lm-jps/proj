/*****************************************************************************
 * Filename: save_packet_to_dayfile.c                                        *
 * Author: Carl                                                              *
 * Create Date: February, 2, 2008                                            *
 * Description: This file contains modules to save housekeeping packets      *
 *              from high speed bus to dayfiles.                             *
 * (C) Stanford University, 2008                                             *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "decode_hk.h"
#include "packets.h"
#include "save_packet_to_dayfile.h"
#include "decode_hk_vcdu.h"
#include <sys/time.h>
#include <time.h>
#include "timeio.h"
#include <dirent.h>
#include "printk.h"

/*************  function prototypes ******************/
int check_filename_apid(char *fn, int apid); 
int  check_filename_pkt_time(char *fn, unsigned short *word_ptr);
int  check_for_pkt_time( HK_Dayfile_Data_t **fdfd, unsigned short *word_ptr);
int  check_dfd_file( unsigned short *word_ptr, int apid);
int check_for_dfd(int apid, HK_Dayfile_Data_t *dfd,  HK_Dayfile_Data_t **fdfd);
int  free_dayfile_data( HK_Dayfile_Data_t **df_head); 
int  free_dayfile_pkt_data( HK_Dayfile_Data_t **df_head); 
double  get_packet_time(unsigned short *word_ptr);
char   *get_dayfilename(int apid, double tc_sec);
int  get_day_from_pkttime(double tc_sec);
int  get_hour_from_pkttime(double tc_sec);
int  get_month_from_pkttime(double tc_sec);
int  get_yr_from_pkttime(double tc_sec);
int  load_dfd_node(unsigned short *word_ptr,HK_Dayfile_Data_t **df_head);
int  load_packet_data( unsigned short *word_ptr, HK_Dayfile_Packet_t **pkt);
void set_time_values(HK_Dayfile_Data_t **dn, double tcsec);
int  write_packet_to_dayfile(HK_Dayfile_Data_t **df_head); 

/*********** extern function prototypes **************/
extern TIME  SDO_to_DRMS_time(int sdo_s, int sdo_ss);
extern void sprint_time (char *at, TIME t, char *zone, int precision);


/*********************** save_hkpkt_dayfile ************************/
/*Pass One Packet to function from VCDU*/
int save_packet_to_dayfile(unsigned short *word_ptr, int apid, HK_Dayfile_Data_t **df_head) 
{
  /* declarations */
  HK_Dayfile_Packet_t *pkt, *temp_pkt,*prev_pkt;
  HK_Dayfile_Data_t *dd, *dd_tmp, *found_dfd;
  int wdf_status;
  int ldfd_status;

  /* load packet in HK_Dayfile_Packet_t node */
  if(!load_packet_data( word_ptr, &pkt))
  {
    printkerr("ERROR at %s, line %d: Could not load packet data to "
              "day file data structure.\n"
               __FILE__, __LINE__);
  }

  /* Handle restart case where stopped data flow and restarted
     so that the dayfile names are not any longer in memory */
  if (*df_head == NULL)
  {
    ldfd_status = load_dfd_node(word_ptr, &dd_tmp);
    if (ldfd_status == ERROR_HK_ENVIRONMENT_VARS_NOT_SET) 
    {
           return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET);
    }
    /* set head of data dayfile link list */
    if(ldfd_status == HK_DF_DAYFILE_LOAD_SUCCESS) 
    {
      *df_head= dd_tmp;
    }
  }

  /* check if data-dayfile node exists or create new one and 
     attach packet node to data-dayfile node  */
  if (*df_head == NULL)
  {
    /* create top data-dayfile node */
    assert(dd = malloc(sizeof(HK_Dayfile_Data_t)));

    /* set apid value and time values in HK_DAYFILE_DATA struct*/
    dd->apid=(short)apid;
    (void)set_time_values(&dd,  get_packet_time(word_ptr));

    /* set dayfile name based on first packet time */
    strcpy(dd->dayfile, get_dayfilename(dd->apid, get_packet_time(word_ptr)));

    /* set null for end record */
    dd->next=NULL;

    /* link packet node to dayfile-data node */
    dd->pkt=pkt;

    /* set head reference to return back */
    *df_head = dd_tmp = dd;
  } /* if top dayfile-data node is null */
  else
  {
    /* check if dfd node for apid exists if so use node to attach new packet node 
       else create new dfd node and attach packet nodeo */
    dd_tmp=  *df_head;
    if (check_for_dfd(apid, dd_tmp, &found_dfd))
    {
      /* check if packet time for this packet is for same day- 
         if not need to create new dayfile filename  */
      if (check_for_pkt_time( &found_dfd, word_ptr))
      {
         ;//printf("save_packet_for_dayfile: dfd node is within pkt time \n");
      }
      else //No match-reset packet time to first packet time  and new dayfile name or use file
      {
        /* check if dayfile exists in dayfile directory and set to dfd node */
        if (check_dfd_file(word_ptr, apid))
        {
          /* write packets there and free dayfile node and pkt nodes */
          wdf_status= write_packet_to_dayfile( &dd_tmp);
          if (wdf_status == 0) 
          {
            return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET);
          }
          wdf_status= free_dayfile_data( &dd_tmp);

          /* load dfd nodes with existing dayfiles for this packet time */
          ldfd_status = load_dfd_node(word_ptr, &dd_tmp);
          if (ldfd_status == ERROR_HK_ENVIRONMENT_VARS_NOT_SET) 
          {
            return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET);
          }

          /* set head node */
          *df_head=dd_tmp;

          /* call to get found_dfd node for this packet*/
          if (check_for_dfd(apid, dd_tmp, &found_dfd)) 
          {
            ;
          }
        }  
        else
        {   
          /* if files not there for given time and apid then reuse dfd node */
          /* write packets there and free -only- pkt nodes */
          wdf_status = write_packet_to_dayfile( &dd_tmp);
          if (wdf_status == 0) 
          { 
            return (ERROR_HK_ENVIRONMENT_VARS_NOT_SET);
          }
          wdf_status = free_dayfile_pkt_data( &dd_tmp); 

          /* reuse existing node and reset time and dayfile name */
          strcpy(found_dfd->dayfile, get_dayfilename(found_dfd->apid, get_packet_time(word_ptr)));
          (void)set_time_values(&found_dfd,  get_packet_time(word_ptr));
        }  
      }/*else no match load dayfile there in dfd node or reuse dfd node by reset time and df name */

      /* set pkt node to pkt node where found apid in dayfile-data node */
      temp_pkt = found_dfd->pkt;// DEC 5 

      /* check where to add pkt data */
      if (temp_pkt)
      {
        /* if 2th to N packet node there-set here*/
        /* search to last pkt node to use to add new pkt node*/
        for(prev_pkt=temp_pkt; temp_pkt; prev_pkt=temp_pkt,temp_pkt = temp_pkt->next);

        /*add new packet node to end of dfd node's packet nodes */
        prev_pkt->next=pkt;        
      }
      else
      {
        /* if first packet node there-set here*/
        found_dfd->pkt=pkt;
      }
    }
    else
    { 
      /* if does not exist make new node for apid */
      assert(dd = malloc(sizeof(HK_Dayfile_Data_t)));
      dd->apid=(short)apid;
      (void)set_time_values(&dd,  get_packet_time(word_ptr));
      strcpy(dd->dayfile, get_dayfilename(dd->apid, get_packet_time(word_ptr)));
      dd->next=NULL;

      /* add to end of dayfile-data nodes link list*/
      dd_tmp=  *df_head;
      while ( dd_tmp->next) 
      { 
        dd_tmp= dd_tmp->next;
      }
      dd_tmp->next = dd;

      /* connect pkt node to data-dayfile node*/
      dd->pkt=pkt;
    }/*else if did not find dfd node for this apid*/
  } /*else if dfd_head is NOT null */
  return (HK_SUCCESS_WRITE_DAYFILE);
}



/******************* get packet time ****************************/
double get_packet_time(unsigned short *word_ptr)
{
  /*declarations */
  TIME pkt_time;
  double tc;
  unsigned int *wptr;
  unsigned int pkt_secs;
  unsigned int pkt_subsecs;

  /* get time in packet seconds */
  wptr = (unsigned int*)(word_ptr+3);
  /* set 0th to 7th bits */
  pkt_secs = (unsigned int) ( (*wptr) >> 24 & 0x000000FF);
  /* set 8th to 15th bits */
  pkt_secs |= (unsigned int) ( (*wptr) >> 8 & 0x0000FF00);
  /* set 24th to 31th bits */
  pkt_secs |= (unsigned int) ( (*wptr) << 24 & 0xFF000000);
  /* set 16th to 23th bits */
  pkt_secs |= (unsigned int) ( (*wptr) << 8 & 0x00FF0000);
                                                                                                    
  /* get time in packet subseconds */
  wptr = (unsigned int*)(word_ptr+5);
  /* set 0th to 7th bits */
  pkt_subsecs =  (unsigned int) ( (*wptr) >> 24 & 0x000000FF);
  /* set 8th to 15th bits */
  pkt_subsecs |= (unsigned int) ( (*wptr) >> 8 & 0x0000FF00);
  /* set 24th to 31th bits */
  pkt_subsecs |= (unsigned int) ( (*wptr) << 24 & 0xFF000000);
  /* set 16th to 23th bits */
  pkt_subsecs |= (unsigned int) ( (*wptr) << 8 & 0x00FF0000);

  /* call function to convert secs and subs to double seconds */ 
  pkt_time = SDO_to_DRMS_time(pkt_secs, pkt_subsecs);
  //pkt_time += 33.0;
  tc=pkt_time;
  return tc;
}



/*********************** load packet data ************************/
int load_packet_data( unsigned short *word_ptr, HK_Dayfile_Packet_t **pkt)
{
  HK_Dayfile_Packet_t *pp;
  char vn[10];
  char *ptr_vn;
  int i;
  unsigned short w;

  /* initialized variables */
  ptr_vn =vn;

  /* create packet data node */
  assert(pp = malloc(sizeof(HK_Dayfile_Packet_t)));

  /* get packet version number -Probably don't need*/
  w = word_ptr[7];
  sprintf(ptr_vn, "%d.%d",  w & 0x00FF, w >> 8  & 0x00FF );

  /* get and set packet size to write to pkt node data struct */
  w = word_ptr[2]; 
  pp->length = (unsigned short)((w  >> 8 & 0x00FF) + (w & 0x00FF));

  /* load packet data in pkt node data struct */
  for (i = 0 ; i < pp->length/2 + 4; i++)
  {
     pp->value[i]=word_ptr[i] ; 
  }
  pp->next=NULL;

  /* test print of content in pkt node's array -comment out*/
  /*****
  for (i=0 ; i < pp->length/2 + 4; i++)
  {
     printf("%-04.4x ",pp->value[i]);
     if (i == 7 || i == 15 || i == 23 ||
         i == 31 || i == 39 || i == 47 || i == 55 || i == 63 ||
         i == 71 || i == 79 || i == 87 || i == 95 || i == 103 ||
         i == 111 || i == 119)
       printf("\n");
  }
  printf("\n");
  ***/

  /*set pkt pointer to created node with data so can send back*/
  *pkt =pp;
  return (HK_DF_LOAD_PACKET_SUCCESS);
}



/***************** check for apid in data-dayfile structure *****************/
int check_for_dfd(int apid, HK_Dayfile_Data_t *dfd,  HK_Dayfile_Data_t **fdfd)
{
  /* declarations */
  HK_Dayfile_Data_t *dfd_t, *found_dfd;

  /* init variables */
  dfd_t=dfd;

  /* loop thru link list for dfd node for given apid */
  while (dfd_t)
  {
    /* if apid equal return pointer to found dfd node */
    if (dfd_t->apid == apid)
    {
      found_dfd=dfd_t;
      *fdfd=found_dfd; 
      return (HK_DF_CHECK_FOUND);
    }
    dfd_t=dfd_t->next;
  }
  found_dfd=NULL;
  return (HK_DF_CHECK_NOT_FOUND);
}



/********************* write_packet_to_dayfile ********************/
int write_packet_to_dayfile(HK_Dayfile_Data_t **df_head) 
{
  /*declarations*/
  FILE *file_ptr;
  HK_Dayfile_Packet_t *pp;
  HK_Dayfile_Data_t *dd;
  char filename[MAX_FILE_NAME];
  char dn[MAX_DIRECTORY_NAME];
  char *p_dn;

  /* get directory */
  p_dn = getenv("HK_DF_HSB_DIRECTORY");
  if(!p_dn) 
  {
    printkerr("ERROR at %s, line %d: Could not get directory environment "
              "variable:<HK_DF_HSB_DIRECTORY>. Set the env variable "
              "HK_DF_HSB_DIRECTORY to existing directory name. \n",
              __FILE__,__LINE__);
    return ( ERROR_HK_ENVIRONMENT_VARS_NOT_SET );
  }
  strcpy(dn,p_dn);
  
  /* loop thru data dayfile structure and write packets to dayfile based on apid */
  dd= *df_head;
  while (dd)  
  { 
    /* create filename */
    /* get apid in decimal */
    /* Use filename based on syntax: hsb_apid_yyyy_ddd_hh_mm_ss_ver.hkt
       and  form complete path including directory name */
    sprintf(filename,"%s/%s",dn, dd->dayfile);

    /* open filename */
    file_ptr = fopen( filename ,"a");
    if(!file_ptr) 
    {
      printkerr("ERROR at %s, line %d: Could not open file at:"
                "<%s>. Check have permission  "
                "to write to directory or check if directory exists. \n",
                __FILE__,__LINE__, filename);
      return ( ERROR_HK_FAILED_OPEN_DAYFILE);
    }

    /* loop thru each packet node in data-dayfile node and write 
       packet to file using packet length stored in packet node*/
    for( pp = dd->pkt; pp ; pp=pp->next)
    {
      fwrite(pp->value,1, pp->length + 7,file_ptr);
    }

    /* close dayfile */
    if (file_ptr) 
    {
      fclose(file_ptr);
    }

    /* go to next data dayfile nodes packet nodes */
    dd=dd->next;
  }
  return (HK_SUCCESS_WRITE_DAYFILE);
}



/*********************** free_dayfile_pkt_data ************************/
int free_dayfile_pkt_data( HK_Dayfile_Data_t **df_head) 
{
  /* declarations */
  HK_Dayfile_Data_t *dfd, *t_dfd;
  HK_Dayfile_Packet_t *pp, *t_pp;

  dfd= *df_head;
  if (dfd)
  {
    /*traverse thru structure and print */
    t_dfd= dfd;
    while (t_dfd)
    {

      /* free all pkt nodes attached to data dayfile node*/
      for ( pp=t_dfd->pkt; pp ; ) 
      {
        t_pp= pp->next;  
        free((HK_Dayfile_Packet_t*) pp); 
        pp=t_pp; 
      }

      /* go to next dayfile node */
      dfd= t_dfd->next;

      /* dayfile data node */
      t_dfd->pkt=NULL;  
      t_dfd=dfd;
    }
  }
  return (HK_SUCCESS_WRITE_DAYFILE);
}



/*********************** free_dayfile_data ************************/
int free_dayfile_data( HK_Dayfile_Data_t **df_head) 
{
  /* declarations */
  HK_Dayfile_Data_t *dfd, *t_dfd;
  HK_Dayfile_Packet_t *pp;

  dfd= *df_head;
  if (dfd)
  {
    /*traverse thru structure and print */
    t_dfd= dfd;
    while (t_dfd)
    {
      /* free all pkt nodes attached to data dayfile node*/
      for ( pp=t_dfd->pkt; pp ; pp=pp->next)
      {
        free((HK_Dayfile_Packet_t*) t_dfd->pkt);
      }

      /* go to next dayfile node */
      dfd= t_dfd->next;

      /* dayfile data node */
      free(t_dfd);
      t_dfd=dfd;
    }
  }
  return (HK_SUCCESS_WRITE_DAYFILE);
}



/******************** get year from pkt time *****************/
int get_yr_from_pkttime(double tc_sec)
{
  /* declarations */
  short year;
  char at[200];

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UT", 0);
  strcat(at,"\0");

  /* get year and return value */
  sscanf(at,"%hd",&year);
  return (year);
}
/******************** get month from pkt time *****************/
int get_month_from_pkttime(double tc_sec)
{
  /* declarations */
  short month;
  char at[200];

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UT", 0);
  strcat(at,"\0");

  /* get month and return value */
  sscanf(at,"%*d.%hd",&month);
  return (month);
}
/******************** get day from pkt time *****************/
int get_day_from_pkttime(double tc_sec)
{
  /* declarations */
  short day;
  char at[200];

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UT", 0);
  strcat(at,"\0");

  /* get month and return value */
  sscanf(at,"%*d.%*d.%hd",&day);
  return (day);
}
/******************** get hour from pkt time *****************/
int get_hour_from_pkttime(double tc_sec)
{
  /* declarations */
  short hour;
  char at[200];

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UT", 0);
  strcat(at,"\0");

  /* get hour and return value */
  sscanf(at,"%*d.%*d.%*d_%hd",&hour);
  return (hour);
}



/**************** get dayfilename ******************/
char * get_dayfilename(int apid, double tc_sec)
{
  /* declarations */
  int year,month,day, hour, minute, second, version;
  char filename[200];
  char at[200];
  char *fp;

  /* init variables */
  version=0;
  fp = filename;

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tc_sec, "UT", 0);
  strcat(at,"\0");

  /* make filename based on what have */
  sscanf(at,"%d.%d.%d_%d:%d:%d.%*s",&year,&month,&day,&hour,&minute,&second);
  sprintf(filename,"hsb_%4.4d_%4.4d_%2.2d_%2.2d_%2.2d_%2.2d_%2.2d_%2.2d.hkt",apid,year,month,day,hour,minute,second,version);
  strcat(filename,"\0");
  return fp;
}

/***************** check for packet time *********************/
int check_for_pkt_time( HK_Dayfile_Data_t **fdfd, unsigned short *word_ptr)
{
  /* declarations */
  HK_Dayfile_Data_t *t_dfd;
  int pkt_yr,pkt_month,pkt_day;

  /* initialized variables */
  t_dfd = *fdfd;

  /* get times from current packet */
  pkt_yr=get_yr_from_pkttime(get_packet_time(word_ptr));
  pkt_month=get_month_from_pkttime(get_packet_time(word_ptr));
  pkt_day=get_day_from_pkttime(get_packet_time(word_ptr));

  /* check if in range of current dayfile */
  if (( t_dfd->year == pkt_yr) && (t_dfd->month == pkt_month) && (t_dfd->day == pkt_day)) 
  {
    /* case of within range of current dayfile-data node's filename */
    return (HK_DF_CHECK_FOUND) ;
  }
  else
  {
    /* case of NOT within range of current dayfile-data node's filename */
    /* case of new day for packet time */
    return (HK_DF_CHECK_NOT_FOUND) ;
  }
}



/**************** set time values ******************/
void set_time_values(HK_Dayfile_Data_t **dn, double tcsec)
{
  /* declarations */
  HK_Dayfile_Data_t *td;
  char at[200];

  /* init variables */
  td= *dn;

  /* convert time code secs and subsecs into year,month,day,hour,minute,second */
  (void)sprint_time (at,  tcsec, "UT", 0);
  strcat(at,"\0");

  /* set time values in dayfile node*/
  sscanf(at,"%hd.%hd.%hd_%hd:%hd:%hd", &td->year, &td->month,
         &td->day, &td->hour, &td->minute, &td->second);
  return;
}


/******************* load dfd node ************************/
int check_dfd_file( unsigned short *word_ptr, int apid)
{
  /*declarations*/
  DIR *dir_p;
  char dn[MAX_DIRECTORY_NAME];
  char *p_dn;
  struct dirent *dir_entry_p;

  /* get directory name for HSB dayfiles */
  p_dn = getenv("HK_DF_HSB_DIRECTORY");
  strcpy(dn,p_dn);

  /* open directory - if some files there do below*/
  if ((dir_p = opendir(dn)) == NULL)
  {
    printkerr("Error at %s, line %d: Could not open directory <%s>. "
	      "Check if environment variable <HK_DF_HSB_DIRECTORY> is set .\n",
              __FILE__,__LINE__,dn);
    return ERROR_HK_FAILED_OPEN_DAYFILE;
  }

  /* read dot hkt files names based on current time of word_ptr filter out some filename.*/
  while( (dir_entry_p = readdir(dir_p)) != NULL ) 
  {
    if( strncmp(dir_entry_p->d_name,"hsb_",4) )
    {
      continue; /* not an dayfile filename - skip*/
    }
    else
    {
      /* For each filename create dfd node and
         set apid,filename,year,month,day, hour for each dfd node */
      /* check if dayfile loading is within range of latest incoming packet */
      if ( (!check_filename_pkt_time(dir_entry_p->d_name, word_ptr)) || (!check_filename_apid(dir_entry_p->d_name,apid)) )
      {
           continue; /* not an dayfile filename within time range- skip*/
      }
      else 
      {
            return (HK_DF_CHECK_FOUND );
      }
    }
  }
  return (HK_DF_CHECK_NOT_FOUND);
}


   
/******************* load dfd node ************************/
int load_dfd_node( unsigned short *word_ptr, HK_Dayfile_Data_t **dfd)
{
  /*declarations*/
  DIR *dir_p;
  HK_Dayfile_Data_t *dd,*t_dd;
  char dn[MAX_DIRECTORY_NAME];
  char *p_dn;
  struct dirent *dir_entry_p;

  /* initialize variables */
  dd=NULL;
  p_dn=dn;

  /* get directory name for HSB dayfiles */
  p_dn = getenv("HK_DF_HSB_DIRECTORY");
  if(!p_dn){
    printkerr("Error at %s, line %d: Could not set environment variable. "
	      "Check if environment variable <HK_DF_HSB_DIRECTORY> is set .\n",
              __FILE__,__LINE__);
    return ERROR_HK_ENVIRONMENT_VARS_NOT_SET;
  }

  /* open directory - if some files there do below*/
  if ((dir_p = opendir(p_dn)) == NULL)
  {
    printkerr("Error at %s, line %d: Could not open directory <%s>. "
	      "Check if environment variable <HK_DF_HSB_DIRECTORY> is set .\n",
              __FILE__,__LINE__,dn);
    return ERROR_HK_FAILED_OPEN_DAYFILE;
  }

  /* read dot hkt files names based on current time of word_ptr filter out some filename.*/
  while( (dir_entry_p = readdir(dir_p)) != NULL ) 
  {
    if( strncmp(dir_entry_p->d_name,"hsb_",4) )
    {
      continue; /* not an dayfile filename - skip*/
    }
    else
    {
      /* For each filename create dfd node and
         set apid,filename,year,month,day, hour for each dfd node */
      /* check if dayfile loading is within range of latest incoming packet */
      if ( !check_filename_pkt_time(dir_entry_p->d_name, word_ptr) )
      {
           continue; /* not an dayfile filename within time range- skip*/
      }

      if (!dd)
      {
        /* create top data-dayfile node */
        assert(dd = malloc(sizeof(HK_Dayfile_Data_t)));

        /* set top pointer **dfd */
        *dfd= dd;
        t_dd=dd;

        /* set apid value and time values in HK_DAYFILE_DATA struct*/
         sscanf(dir_entry_p->d_name, "hsb_%hd_%hd_%hd_%hd_%hd_%hd_%hd", 
                &dd->apid,&dd->year,&dd->month,&dd->day,&dd->hour,&dd->minute, &dd->second );

        /* set filename -12-18-2007 */
        strcpy(dd->dayfile, dir_entry_p->d_name);
        strcat(dd->dayfile, "\0");
        dd->next=NULL;
        dd->pkt=NULL;
      }
      else
      {
        /* create next data-dayfile node */
        assert(t_dd = malloc(sizeof(HK_Dayfile_Data_t)));

        /* set apid value and time values in HK_DAYFILE_DATA struct*/
        sscanf(dir_entry_p->d_name, "hsb_%hd_%hd_%hd_%hd_%hd_%hd_%hd",
               &t_dd->apid,&t_dd->year,&t_dd->month,&t_dd->day,&t_dd->hour,&t_dd->minute, &t_dd->second);

        /* set filename -12-18-2007 */
        /* create next data-dayfile node */
        strcpy(t_dd->dayfile, dir_entry_p->d_name);
        strcat(t_dd->dayfile,"\0");
        t_dd->next=NULL;
        t_dd->pkt=NULL;

        /* move to point to next node*/
        dd->next=t_dd;
        dd=t_dd;
      }
    }
  }
  // Check value returned here
  return HK_DF_DAYFILE_LOAD_SUCCESS;
}



/*********** check filename pkt time ************/
int check_filename_pkt_time(char *fn, unsigned short *word_ptr) 
{
  /* declarations */
  int pkt_yr,pkt_month,pkt_day;
  int fn_yr,fn_month,fn_day;

  /* get times from current packet */
  strcat(fn,"\0");
  pkt_yr=get_yr_from_pkttime(get_packet_time(word_ptr));
  pkt_month=get_month_from_pkttime(get_packet_time(word_ptr));
  pkt_day=get_day_from_pkttime(get_packet_time(word_ptr));
  sscanf(fn, "hsb_%*d_%d_%d_%d_%*d_%*d_%*d.%*s",&fn_yr, &fn_month,&fn_day);

  /* check filename time is within range */
  if ( (pkt_yr == fn_yr) && (pkt_month == fn_month)  && (pkt_day == fn_day) )
  {
    return (HK_DF_CHECK_FOUND);
  }
  else 
  {
    return (HK_DF_CHECK_NOT_FOUND) ;
  }

}


/*********** check filename apid ************/
int check_filename_apid(char *fn, int apid) 
{
  /* declarations */
  int fn_yr,fn_month,fn_day;
  int fn_apid;

  /* get times from current packet */
  strcat(fn,"\0");
  sscanf(fn, "hsb_%d_%d_%d_%d_%*d_%*d_%*d.%*s",&fn_apid,&fn_yr, &fn_month,&fn_day);

  /* check filename time is within range */
  if (apid == fn_apid)
  {
    return (HK_DF_CHECK_FOUND);
  }
  else 
  {
    return (HK_DF_CHECK_NOT_FOUND) ;
  }
}

