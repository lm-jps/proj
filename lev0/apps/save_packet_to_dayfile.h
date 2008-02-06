#ifndef __SAVE_PACKET_TO_DAYFILE_H
#define __SAVE_PACKET_TO_DAYFILE:w


/* Structures */

/***************** HK Dayfile Packet struct *************************/
typedef struct HK_Dayfile_Packet_struct 
{
  unsigned short length;
  unsigned short value[2000];
  struct HK_Dayfile_Packet_struct *next;

} HK_Dayfile_Packet_t;

/***************** HK Dayfile Data struct *************************/
typedef struct HK_Dayfile_Data_struct 
{
  /* current filename fields and information
     (i.e, hsb_<apid-in-dec>_<yyyy>_<MM>_<dd>_<hh>_<mm>_<ss>_<vv>.hkt) */
  int short apid;
  char dayfile[100];
  short year;
  short day_yr;
  short month;
  short day;
  short hour;
  short minute;
  short second;
  short version;
  /* Note:Lev 0 documents filename:hsb_<apid-in-dec>_<yyyy>_<ddd>_<hh>_<mm>_<ss>_<vv>.hkt)

  /* packet lenght of bytes per packet */
  unsigned short length;

  /* packet values  */
  struct HK_Dayfile_Packet_struct *pkt;

  /* Next pointer */
  struct HK_Dayfile_Data_struct *next;

} HK_Dayfile_Data_t;


#endif
