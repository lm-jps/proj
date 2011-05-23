// This file contains a bunch of DSDS code, ripped from the EGSE CVS tree and needed by ingest_tlm and soc_pipe_scp.
// Art
// April 19, 2011

#include <inttypes.h>
#include <dirent.h>
#include "jsoc.h"
#include "egsehmicomp.h"

typedef struct APID_Pointer_HK_Configs_struct    
{
  int apid;                              /*Make either hex or i
					   decimal value*/
  struct HK_Config_Files_struct *ptr_hk_configs;/*Pointer to HK 
						  Configurations struct*/
  struct APID_Pointer_HK_Configs_struct *next;/* Link List's next node */
} APID_Pointer_HK_Configs;

typedef struct HK_Config_Files_struct  
{
  int apid_number;
  char packet_id_type[MAX_PACKET_ID_TYPE];
  char file_version_number[MAX_CHAR_VERSION_NUMBER];/*Value for the apid-version<version-
				                      number>.txt file */
  char parameter_version_number[MAX_CHAR_VERSION_NUMBER];/*Values in GTCIDS map file for 
				   HMI_VER_NUM_SEQ_STATUS */
  char date[MAX_DATE_SIZE];      /*Example: 2005/07/06 or 10 characters*/
  int number_bytes_used;        /* Example: 112 in apid file and stanford file*/
  struct Keyword_Parameter_struct *keywords;/*Pointer to array or vector of
					      keyword_parameter structures*/
  struct HK_Config_Files_struct  *next;
} HK_Config_Files;

typedef struct APID_HKPFD_Files_struct  
{
  char version_number[MAX_CHAR_VERSION_NUMBER];
  int apid; 
  char directory_name[MAX_DIRECTORY_NAME];
  char filename[MAX_FILE_NAME];
  struct APID_HKPFD_Files_struct *next;
}   APID_HKPFD_Files;

typedef  struct GTCIDS_Version_Number_struct   
{
  char file_version_number[MAX_CHAR_VERSION_NUMBER];
  char hmi_id_version_number[MAX_CHAR_VERSION_NUMBER];
  char aia_id_version_number[MAX_CHAR_VERSION_NUMBER];/*Values in GTCIDS map file */
  char change_date[MAX_SIZE_CHANGE_DATE];
  char change_time[MAX_SIZE_CHANGE_TIME];
  struct GTCIDS_Version_Number_struct  *next;
}   GTCIDS_Version_Number;

typedef struct ALG_Conversion_struct
{
  int     number_of_coeffs;
  double  coeff [MAX_NUMBER_COFFICIENTS];
} ALG_Conversion;

typedef struct Keyword_Parameter_struct  {
  char telemetry_mnemonic_name[MAX_HK_MNM];/*Make 30-original name
					     from Lockheed Database and
					     STANFORD_TLM_HMI_AIA.txt files*/
  char keyword_name[MAX_HK_KYWD]; /*Make size at least 10 - 8 chararacters
				    of less parameter passed to L0P Decode
				    component */
  int start_byte;       /*L0P Decode will use this to parse packet data */
  int start_bit_number;
  int bit_length;       /*L0P Decode will use this to determine bit length
			  to parse packet data */
  char conv_type;       /* type of conversion facTtor */
  char type[MAX_HK_VALUE_TYPE];     /* Make at least 5 -Examples: UB, UL1, i
				       UL2, etc  or 3 characters*/
  struct ALG_Conversion_struct *alg_ptr;
  struct DSC_Conversion_struct *dsc_ptr;
  struct Keyword_Parameter_struct *next;
} Keyword_Parameter;

/******************* VCDU Packet struct ********************/
typedef struct IM_PDU_Packet_struct 
{
  /* IM_PDU Header fields. */
  uint8_t  im_pdu_id;
  uint64_t im_pdu_count;
  /* M_PDU header */
  uint8_t  spare;
  uint16_t fhp;

  /* Linked list of CCSDS packets. */
  struct CCSDS_Packet_struct *packets;

} IM_PDU_Packet_t;

typedef struct HK_Keywords_Format_struct 
{
  char keyword_name[KEYWORD_NAME_SIZE]; /*Size is limited to 8 
                                         *from HMI-EGSE-FS document 
                                         */
  unsigned int keyword_value;           /*Actual value of keyword from 
                                         *hk bit stream packet
                                         */
  char type[KEYWORD_TYPE_SIZE] ;        /*Used by Decoder to parse keyword
                                         *values. Size is limited to 3 from 
                                         *HMI-EGSE-SS-FS document 
                                         */
  char telemetry_mnemonic_name[TELEMETRY_MNEMONIC_SIZE];/*Used to debug  */
  int start_byte ;                      /*Used by Decoder Parsing keyword 
                                         *values. For example: 8th byte, 
                                         *14th byte, 18th byte, etc. Used 
                                         *to debug
                                         */
  int start_bit;
  int bit_length;
  char conv_type;                       /*keyword can be conv type R,A,D */
  struct ALG_Conversion_struct  *alg;   /*analog config data to decode   */
  struct DSC_Conversion_struct  *dsc;   /*digital config data to decode  */
  struct HK_Keywords_Format_struct *next;/*Currently there are about 50
                                           keywords 
                                          */
} HK_Keywords_Format;

typedef struct DSC_Conversion_struct
{
  int     low_range;
  int     high_range;
  char    dsc_value[MAX_VALUE_DSC];
  struct  DSC_Conversion_struct  *next;
} DSC_Conversion;


#define ERRMSG(__msg__) printkerr("ERROR at %s, line %d: " #__msg__"\n",__FILE__,__LINE__);

// globals
/* Decompression contexts */
int numcontexts=0;
Decompress_Context_t *Context[MAXCONTEXTS];
APID_Pointer_HK_Configs *global_apid_configs;
GTCIDS_Version_Number *global_gtcids_vn;
static int isHMI = 0;

// From fits.c
static int fits_copy_header(char *out, char *in, int headlen)
{
  int newlen=0;

  while (headlen>0)
  {
    if (strncmp(in, "SIMPLE  =", 9) &&
	strncmp(in, "BITPIX  =", 9) &&	
	strncmp(in, "NAXIS   =", 9) &&
	strncmp(in, "NAXIS1  =", 9) &&
	strncmp(in, "NAXIS2  =", 9) &&
	strncmp(in, "DATE    =", 9) &&
	strncmp(in, "END     ", 8))
    {
      memcpy(out, in, 80);
      out += 80;
      newlen += 80;
    }
    in += 80;
    headlen -= 80;
  }
  return newlen;
}

#define WRITE(fh,buf,len)  (fwrite(buf, len, 1, fh)==1)
#define CLOSE(fh) fclose(fh)

int fits_writeimage_raw(const char* file, int compress, int headlen, 
			char *header, int width, int height,short *data)
{
  time_t t;
  struct tm tm;
  FILE *fh;
  int i, rem, hl,nzhead, nzdata;
  unsigned char *zdata;
  char *head;
  char date[20], *p, tmp;
  unsigned int npix,buflen;
  int ndata_blocks, nhead_blocks, len;
  char buf[1024];

  if (headlen % 80)
  {
    printkerr("Error: Fits header length must be a multiple of 80.\n");
    return 1;
  }
    
  if (!(fh=fopen(file,"w")))
  {
    printkerr("Couldn't open file %s.\n",file);
    return 1;
  }
  

  /* Write standard FITS keywords. */
  head = malloc(headlen+2*2880);
  memset(head,' ',headlen+2*2880);
  p = head;
  snprintf(p,80,"SIMPLE  =                    T");
  p += 80;
  snprintf(p,80,"BITPIX  =                   16 / Number of bits per data pixel");
  p += 80;
  snprintf(p,80,"NAXIS   =                    2 / Number of data axes");
  p += 80;
  snprintf(p,80,"NAXIS1  =               %6d",width);  
  p += 80;
  snprintf(p,80,"NAXIS2  =               %6d",height);
  p += 80;
  snprintf(p,80,"BLANK   =               %6d",-32768);
  p += 80;
  /* Set up a string with the present UTC time and date. */
  t = time(NULL);
  gmtime_r(&t, &tm);
  strftime(date,20,"%F %T",&tm);
  
  snprintf(p,80,"DATE    = '%s'",date);
  p += 80;
  if (headlen>0 && header !=NULL) 
    p += fits_copy_header(p,header,headlen);
  snprintf(p,80,"END");
  p += 80;
    
  while ((p-head)%2880)
  {
    memset(p,' ',80);
    p += 80;
  }
  hl = p-head;
  nhead_blocks = hl/2880;
  for(i=0; i<hl; i++)
    if (head[i]=='\0')
      head[i] = ' ';

  /* Write header block. */ 
  if (!WRITE(fh,head,hl))
  {
    printkerr( "Failed to write header data to file %s\n",file);
    goto bailout;
  }

  npix = width*height;

  if (compress)
  {
    // allow 50% growth in size in case of bad compression
    buflen = 3*npix;
    if ((zdata= malloc(buflen))==NULL)
    {
      printkerr("Could't allocate zdata buffer\n");
      goto bailout;
    }
      
    /* Generate compressed data block. */ 
    nzdata = rice_encode2(data, npix, zdata, buflen, 64);
    if (nzdata < 0) 
    {
      printkerr( "Data compression failed with error code %d\n", 
		 nzhead);
      free(zdata);
      goto bailout;
    }

    /* Generate uncompressed mini header block. */ 
    fprintf(fh, "%d %d %d %d\n", 16, 64, 2*npix, nzdata);

    /* Write compressed data. */ 
    if (!WRITE(fh, zdata, nzdata))
    {
      printkerr( "Failed to write compressed data to file %s\n",file);
      free(zdata);
      goto bailout;
    }
    free(zdata);
  }
  else
  {
    /* The FITS standard stipulates storing binary data in big endian format. */
#if __BYTE_ORDER == __LITTLE_ENDIAN
    p = (char *)data;
    for (i=0; i<npix; i++)
    {
      tmp = *(p+1); *(p+1) = *p;  *p = tmp;
      p += 2;
    }
#endif
    
    /* Write data */
    if (!WRITE(fh,data,sizeof(unsigned short)*npix))
    {
      printkerr( "Failed to write data to file %s\n",file);
      goto bailout;
    }
  
    /* Pad with zeros to get an integer number of 2880 byte records. */
    memset(head,0,2880);
    rem = 2880 - ((sizeof(unsigned short)*npix) % 2880);
    if (!WRITE(fh, head, rem))
    {
      printkerr( "Failed to write zero padding to file %s\n",file);
      goto bailout;
    }      
  
    /* Restore data buffer to native format. */
#if __BYTE_ORDER == __LITTLE_ENDIAN
    p = (char *)data;
    for (i=0; i<npix; i++)
    {
      tmp = *(p+1); *(p+1) = *p; *p = tmp;
      p += 2;
    }
#endif
  }
  /* File successfully written to disk. Close it and return 0 to signal
     success. */
  free(head);
  CLOSE(fh);
  return 0;


 bailout:
  /* Something went wrong. Close the file, remove it and return 
     an error code. */
  free(head);
  CLOSE(fh);
  unlink(file);
  return 1;  
}
#undef WRITE
#undef CLOSE


// From printk.c

/* Set printk output functions. */
void printk_set(int (*std)(const char *fmt, ...),
                int (*err)(const char *fmt, ...))
{
  printk = std;
  printkerr = err;
}

/* Default error message output function. */
static int fprintf_wrap(const char *fmt, ...)
{
  int val;
  va_list args;

  va_start(args,fmt);
  val = vfprintf(stderr,fmt,args);
  va_end(args);
  return val;
}

/* Default error message output function. */
static int printf_wrap(const char *fmt, ...)
{
  int val;
  va_list args;

  va_start(args,fmt);
  val = vprintf(fmt,args);
  va_end(args);
  return val;
}

int (*printk)(const char *fmt, ...) = &printf_wrap; // default = output to stdout.
int (*printkerr)(const char *fmt, ...) = &fprintf_wrap;  //default = output to stderr. 

// From decode_hk.c

static void get_version_number(unsigned short *wd_ptr, char *ptr_vn) 
{
  /* declarations */
  unsigned short w;
  w = wd_ptr[7];
  /* parse high and low bytes of version number and convert from hex to decimal value */
  sprintf(ptr_vn, "%d.%d",  w & 0x00FF, w >> 8  & 0x00FF );
}

static int load_engr_values(HK_Keywords_Format *hk_keywords, HK_Keyword_t **kw_head)
{
  /* declarations */
  HK_Keyword_t  *kwt=NULL;
  DSC_Conversion *dsc;
  int i;
  /* Loop thru HK_Keywords_Format link list structure and find 
     first value to load */
  while(hk_keywords)
  {
    if (!kwt)
      *kw_head = kwt = (HK_Keyword_t*) malloc (sizeof (HK_Keyword_t));
    else
    {
      kwt->next = (HK_Keyword_t*) malloc (sizeof (HK_Keyword_t));
      kwt = kwt->next;
    }
    memset(kwt,0,sizeof(HK_Keyword_t));
    kwt->next = NULL;
    /* load keyword name */ 
    strcpy(kwt->fitsname, hk_keywords->keyword_name); 
    strcpy(kwt->name, hk_keywords->telemetry_mnemonic_name); 
    /* load raw value */
    kwt->raw_value = hk_keywords->keyword_value;
    /* load engr value and type  */
    if ( hk_keywords->conv_type == 'R')
    {
      if (!strcmp( hk_keywords->type , "UB"))
      {
        kwt->eng_type = DRMS_TYPE_UINT8;
        kwt->eng_value.uint8_val = hk_keywords->keyword_value;
      }
      else if (!strcmp( hk_keywords->type , "SB"))
      {
        kwt->eng_type = DRMS_TYPE_INT8;
        kwt->eng_value.int8_val = hk_keywords->keyword_value;
      }
      else if (!strcmp( hk_keywords->type , "IU1"))
      {
        kwt->eng_type = DRMS_TYPE_UINT16;
        kwt->eng_value.uint16_val = hk_keywords->keyword_value;

      }
      else if (!strcmp( hk_keywords->type , "IS1"))
      {
        kwt->eng_type = DRMS_TYPE_INT16;
        kwt->eng_value.int16_val = hk_keywords->keyword_value;
      }
      else if (!strcmp( hk_keywords->type , "IL1"))
      {
        kwt->eng_type = DRMS_TYPE_INT32;
        kwt->eng_value.int32_val = hk_keywords->keyword_value;
      }
      else if (!strcmp( hk_keywords->type , "UL1"))
      {
        kwt->eng_type = DRMS_TYPE_UINT32;
        kwt->eng_value.uint32_val = hk_keywords->keyword_value;
      }
      else 
      {
        printkerr("ERROR at %s, line %d: Type '%s' not handled.\n", __FILE__, __LINE__, hk_keywords->type);
        return ERROR_HK_UNHANDLED_TYPE;
      }
    }/* end-if R type */
    else if ( hk_keywords->conv_type == 'D')
    {
      if (!strcmp( hk_keywords->type,"UB") || !strcmp( hk_keywords->type ,"IU1") || !strcmp( hk_keywords->type,"UL1"))
      {
        kwt->eng_type = DRMS_TYPE_STRING;
        for (dsc = hk_keywords->dsc; dsc ; dsc= dsc->next)
        {
          if  (  hk_keywords->keyword_value  >= dsc->low_range &&
                 hk_keywords->keyword_value  <=  dsc->high_range )
          {
	    kwt->eng_value.string_val =(char *) malloc(sizeof(dsc->dsc_value) + 1);
            strcpy( kwt->eng_value.string_val , dsc->dsc_value);
            break;
          }
        }
        if ( !kwt->eng_value.string_val)
        {
	   kwt->eng_value.string_val =(char *) malloc(sizeof("NO_VALUE") + 1);
           strcpy( kwt->eng_value.string_val , "NO_VALUE");
           printkerr("WARNING at %s, line %d: There are no DSC data lines to set digital type keyword '%s'.\n"
                     "Setting value to NO_VALUE. Check config files and run script to update config files.\n",
                     __FILE__, __LINE__, hk_keywords->keyword_name);
        }
      }
      else 
      {
        printkerr("ERROR at %s, line %d: Type '%s' not handled.\n", __FILE__, __LINE__, hk_keywords->type);
        return ERROR_HK_UNHANDLED_TYPE;
      }
    } /* end else if D type */
    else if (hk_keywords->conv_type == 'A')
    {
      /* Calculation is as follows: 
         (1) get raw value
         (2) get number of coeffs
         (3) get each coeff value
         (4) engr value = coeff[0] * pow(raw, 0) + coeff[1] * pow(raw, 1) + .. coeff[4] * pow(raw, 4);
        */
      if (!strcmp( hk_keywords->type,"UB") || !strcmp( hk_keywords->type ,"IU1") || !strcmp( hk_keywords->type,"UL1"))
      {
        if(hk_keywords->alg)
        { 
          for (i=0; i  <  hk_keywords->alg->number_of_coeffs; i++)
          {
            kwt->eng_value.double_val += hk_keywords->alg->coeff[i] * pow(hk_keywords->keyword_value, i);
          }
        }
        else
        { /* handle case when no ACON line which contains coeff. values in config data */
          kwt->eng_value.double_val = 0;
          printkerr("WARNING at %s, line %d: Missing ACON line for keyword '%s'."
                    " Default engr.value set to zero.\n"
                    "Check config files for missing ACON lines for keyword.\n",
                     __FILE__, __LINE__, hk_keywords->keyword_name);
        }
        kwt->eng_type = DRMS_TYPE_DOUBLE;/*set to given type or double */
      }
      else
      {
        printkerr("ERROR at %s, line %d: Engr Type '%s' not handled.\n", __FILE__, __LINE__, hk_keywords->type);
        return ERROR_HK_UNHANDLED_TYPE;
      }
    }/* else if A type */
    else
    {
          printkerr("ERROR at %c, line %d: Conversion Type '%s' not handled.\n", __FILE__, __LINE__, hk_keywords->conv_type);
          return ERROR_HK_UNHANDLED_TYPE;
    }
    hk_keywords = hk_keywords->next;
  } /*while-loop*/
  return SUCCESSFUL; /* successful*/
}

HK_Keywords_Format *load_hk_configs(HK_Config_Files *config) 
{
  HK_Keywords_Format *new_kw, *new_kw_head;
  Keyword_Parameter *config_kw;
  DSC_Conversion *config_dsc, *prev_new_dsc, *new_dsc;
  ALG_Conversion *config_alg, *new_alg;
  int i;
  /* Initialized variables */
  config_kw=config->keywords;
  /* Allocate memory to HK_Keywords_Format structure */
  if (config_kw == NULL)      
  {
    ERRMSG("Null pointer input.");
    return NULL;
  }
  /* load values to HK_Keywords_Format structure while not null */
  new_kw_head = new_kw = NULL;
  while (config_kw)  
  {
    if (!new_kw)
      new_kw = new_kw_head = malloc(sizeof(HK_Keywords_Format));
    else
    {
      new_kw->next = malloc(sizeof(HK_Keywords_Format));
      new_kw = new_kw->next;
    }
    new_kw->next = NULL;
    /* Load config values */
    strcpy(new_kw->keyword_name, config_kw->keyword_name);
    /*for debugging only */
    strcpy(new_kw->telemetry_mnemonic_name,
	   config_kw->telemetry_mnemonic_name);
    /* for parsing later */
    new_kw->start_byte = config_kw->start_byte;  
    new_kw->start_bit = config_kw->start_bit_number;  
    /* for parsing later */
    strcpy( new_kw->type, config_kw->type);
    /* for parsing later */
    new_kw->bit_length = config_kw->bit_length;
    new_kw->conv_type = config_kw->conv_type;
    /* get DSC values for parsing later and setting engr values */
    if( config_kw->conv_type == 'D')
    {
      new_kw->dsc = (DSC_Conversion*)NULL;/*ADDED 6-26-2006 */
      config_dsc=config_kw->dsc_ptr;
      prev_new_dsc=NULL;
      while (config_dsc)
      {
        /* Set low and high range and string value */
        new_dsc=(DSC_Conversion*) malloc(sizeof(DSC_Conversion));
        if (!prev_new_dsc) 
        {
          new_kw->dsc= new_dsc;
        }
        new_dsc->low_range = config_dsc->low_range;
        new_dsc->high_range = config_dsc->high_range;
        strcpy(new_dsc->dsc_value,config_dsc->dsc_value);
        new_dsc->next = (DSC_Conversion *)NULL;
        if(prev_new_dsc)
        {
          prev_new_dsc->next = new_dsc;
        }
        /* Go to next node for values */
        config_dsc= config_dsc->next;
        prev_new_dsc=new_dsc; 
        new_dsc= new_kw->dsc->next;
      }
    }
    else
    {
       new_kw->dsc = (DSC_Conversion*)NULL;
    }
    /* get ALG values for setting engr values */
    if( config_kw->conv_type == 'A')
    {
      config_alg=config_kw->alg_ptr;
      if (config_alg)
      {
        /* create node & set low and high range and string value */
        new_alg = (ALG_Conversion*) malloc(sizeof(ALG_Conversion));
        new_kw->alg= new_alg;
        new_alg->number_of_coeffs = config_alg->number_of_coeffs;
        for(i=0; i < new_alg->number_of_coeffs; i++)
        {
          new_alg->coeff[i] = config_alg->coeff[i];
        }
      }
      else 
      {
        /* set to null -missing ACON line case */
        new_kw->alg= (ALG_Conversion*)NULL;
      }

    }
    else
    {
       new_kw->alg = (ALG_Conversion*)NULL;
    }

    /* Go to next keyword-parameter link list 
       node in HK_Config_Files structure */
    config_kw = config_kw->next;
  }
  return new_kw_head;
} 

static int load_hk_values(unsigned short *word_ptr, HK_Keywords_Format *ptr_hk_keywords)    
{
  /* declarations */
  HK_Keywords_Format *top_hk_keywords;
  unsigned int *w;
  unsigned char *byte_ptr;
  top_hk_keywords = ptr_hk_keywords;
  unsigned int keep_bits; /* 32 bits */
  int i;
  /* Initialize variable */
  byte_ptr= (unsigned char*)(word_ptr+3);
  /* Loop  thru HK_Keywords_Format link list structure and find  first value to load */
  while ( ptr_hk_keywords != NULL )  
  {
    if(!strcmp(ptr_hk_keywords->type, "UB") || !strcmp(ptr_hk_keywords->type, "SB") ) 
    {
      /* set keep bits to zero */
      keep_bits=0;
      /* get raw value */
      w = (unsigned int*) &(byte_ptr[ptr_hk_keywords->start_byte]);
      /* calculate keep bits mask to use*/
      keep_bits = (unsigned int) (pow( 2, ptr_hk_keywords->bit_length) - 1) ;
      /* adjust based on bit length  and bit position */
      ptr_hk_keywords->keyword_value = (unsigned  char)(( *w >> (8 - (ptr_hk_keywords->start_bit + ptr_hk_keywords->bit_length))) & keep_bits);
    } 
    else if( (!strcmp(ptr_hk_keywords->type, "IS1") || !strcmp(ptr_hk_keywords->type, "IU1")) && ptr_hk_keywords->bit_length <= 16) 
    {
      /* set keep bits to zero */
      keep_bits=0;
      /* get raw value */
      w = (unsigned int*) &(byte_ptr[ptr_hk_keywords->start_byte]);
      /* calculate keep bits mask to use*/
      keep_bits = (unsigned int) (pow( 2, ptr_hk_keywords->bit_length) - 1) ;
      /* set 0th to 7th bits */
      ptr_hk_keywords->keyword_value =  (unsigned  short int) (( (*w  ) >> 8) & 0x00FF );
      /* set 8th to 15th bits */
      ptr_hk_keywords->keyword_value |= (unsigned  short int) (( (*w ) << 8 ) & 0xFF00 );
      /* adjust based on bit length  and bit position */
      ptr_hk_keywords->keyword_value =  (unsigned int) (( ptr_hk_keywords->keyword_value >> (16 - (ptr_hk_keywords->start_bit + ptr_hk_keywords->bit_length))) &  keep_bits);
    }
    else if (!strcmp(ptr_hk_keywords->type, "UL1") || !strcmp(ptr_hk_keywords->type, "IL1")) 
    {  
      /* set keep bits to zero */
      keep_bits=0;
      /* calculate keep bits mask to use*/
      keep_bits = (unsigned int) (pow( 2, ptr_hk_keywords->bit_length ) - 1 ) ;
      /* get raw value */
      w = (unsigned int*) &(byte_ptr[ptr_hk_keywords->start_byte]);
      /* set 0th to 7th bits */
      ptr_hk_keywords->keyword_value =	(unsigned int) ( (*w) >> 24 & 0x000000FF);
      /* set 8th to 15th bits */
      ptr_hk_keywords->keyword_value |=	(unsigned int) ( (*w) >> 8 & 0x0000FF00);
      /* set 24th to 31th bits */
      ptr_hk_keywords->keyword_value |=	(unsigned int) ( (*w) << 24 & 0xFF000000);
      /* set 16th to 23th bits */
      ptr_hk_keywords->keyword_value |=	(unsigned int) ( (*w) << 8 & 0x00FF0000);
      /* adjust for bit length */
      ptr_hk_keywords->keyword_value =  (ptr_hk_keywords->keyword_value >> (32 - (ptr_hk_keywords->start_bit + ptr_hk_keywords->bit_length ))) & keep_bits ; 
    }                 
    else  
    {
      printkerr("ERROR at %s, line %d: Did not find this bit length for keyword %s\n",
                 __FILE__, __LINE__, ptr_hk_keywords->keyword_name );
      //return ERROR_HK_INVALID_BITFIELD_LENGTH;
    }
  ptr_hk_keywords = ptr_hk_keywords->next;           
  } 
  ptr_hk_keywords= top_hk_keywords;
  return (SUCCESSFUL); /*no errors */
} /*END-module :load_hk_values()*/

static void deallocate_hk_keywords_format(HK_Keywords_Format *head)   
{
  HK_Keywords_Format *tmp;
  DSC_Conversion *d_tmp, *d_head;

  while(head)
  {
    /*clear dsc nodes*/
    d_head=head->dsc;
    while (d_head)
    {
      d_tmp=d_head->next;
      free(d_head);
      d_head= d_tmp;
    }

    /*clear alg nodes*/
    if (head->alg)
    {
      free(head->alg);    
    }

    /*clear HK_Keyword_Format nodes */
    tmp = head->next;
    free(head);
    head = tmp;
  }
}

void deallocate_hk_keyword(HK_Keyword_t *head)   
{
  HK_Keyword_t *tmp;

  while(head)
  {
    tmp = head->next;
    if (head->eng_type == DRMS_TYPE_STRING && head->eng_value.string_val)
    {
      free(head->eng_value.string_val);
    }
    free(head);
    head = tmp;
  }
}

GTCIDS_Version_Number * read_gtcids_hk_file(APID_Pointer_HK_Configs *top_apid_ptr)
{
  /*declarations*/
  char *hk_gtcids_directory ;
  char *hk_gtcids_filename;
  char hk_gtcids_directory_filename[MAX_FILE_NAME];
  FILE *file_ptr;
  char line[MAXLINE_IN_FILE];
  GTCIDS_Version_Number *top_ptr_gtcids_vn;
  GTCIDS_Version_Number *ptr_gtcids_vn;

  /* get directory and file name */
  hk_gtcids_directory = getenv("HK_CONFIG_DIRECTORY");
  hk_gtcids_filename= getenv("HK_GTCIDS_FILE");
  if(hk_gtcids_filename == NULL) 
    hk_gtcids_filename = "gtcids.txt";
  if(hk_gtcids_directory == NULL) 
    hk_gtcids_directory = "../../tables/hk_config_file";
  sprintf(hk_gtcids_directory_filename, "%s/%s", hk_gtcids_directory,
	  hk_gtcids_filename); 
  /*open file & read  data into GTCIDS_Version_Number link list structure*/
  file_ptr = fopen(hk_gtcids_directory_filename, "r");
  if(!file_ptr)
  {
    printkerr("Error:Couldn't open Ground To Code IDS  file %s.\n",hk_gtcids_directory_filename);
    printkerr("Set the enviroment variables HK_CONFIG_DIRECTORY"
              " to point to config directory and HK_GTCIDS_FILE"
              " environment variable to point to the correct file name\n");
    //exit (1);
    return NULL;
  }
  top_ptr_gtcids_vn = NULL;
  while( fgets(line, MAXLINE_IN_FILE, file_ptr) != NULL )
  {
    if(line[0] == '#') 
      continue; /* skip comments */
    else  
    {
      if (top_ptr_gtcids_vn == NULL)    
      {
	top_ptr_gtcids_vn = ptr_gtcids_vn = malloc(sizeof(GTCIDS_Version_Number));
	ptr_gtcids_vn->next = (GTCIDS_Version_Number*)NULL;
      } 
      else  
      {
	ptr_gtcids_vn->next = malloc(sizeof(GTCIDS_Version_Number));
	ptr_gtcids_vn = ptr_gtcids_vn->next;
        ptr_gtcids_vn->next = (GTCIDS_Version_Number*)NULL;
      }
      ptr_gtcids_vn->next= NULL;
      /*Parse key values- HMI ID and Stanford ID */
      /*Locate column for HMI_ID value and STANFORD_TLM_HMI_AIA values*/
      /*HMI_ID Column is 4th column with | as delimiter*/
      /*STANFORD_TLM_HMI_AIA values  are in 7th column*/
      /*Assume always above, otherwise more code required here*/
      sscanf( line,
	      "%s %s |%*s | %*s | %s | %s | %*s | %s | %*s",
	      ptr_gtcids_vn->change_date, 
	      ptr_gtcids_vn->change_time, 
	      ptr_gtcids_vn->hmi_id_version_number,
	      ptr_gtcids_vn->aia_id_version_number,
	      ptr_gtcids_vn->file_version_number);
    } 
  } /* END-for fgets line */
  fclose(file_ptr);

  global_gtcids_vn = top_ptr_gtcids_vn; /* set global */
  return(top_ptr_gtcids_vn);
}/* END-Module: read_gtcids_hk_file */

char * find_file_version_number(GTCIDS_Version_Number *top, char p_version_number[])
{
  /*declarations*/
  GTCIDS_Version_Number  *tmp_ptr;

  for(tmp_ptr=top;tmp_ptr;tmp_ptr=tmp_ptr->next)    
  {
    if( !strcmp(tmp_ptr->hmi_id_version_number, p_version_number)) 
    {
      return tmp_ptr->file_version_number;
    }
  } /* End-for  */
  return ("");
}

APID_HKPFD_Files* read_all_hk_config_files(char f_version_number[])
{
  /*declarations */
  char dirname[200];
  char *dn;
  APID_HKPFD_Files *top, *p; 
  DIR *dir_p;
  struct dirent *dir_entry_p;

  /* intialize variables */
  int file_loaded_flag=0;
  memset(dirname, 0x0, sizeof(dirname));

  /* get directory name */
  dn = getenv("HK_CONFIG_DIRECTORY");
  if(dn == NULL) 
  {
    printkerr("Error at %s, line %d: Could not get directory environment\n"
              "variable:<HK_CONFIG_DIRECTORY>. Set the env variable "
	      "HK_CONFIG_DIRECTORY to point to the config file directory.\n",
              __FILE__,__LINE__,dn);
    return NULL;
  }
  /* Add file version number to directory path */
  strcpy( dirname, dn);
  strcat(dirname, f_version_number);
  strcat(dirname,"\0");

  /* open directory */
  if ((dir_p = opendir(dirname)) == NULL)
  {
    printkerr("Error at %s, line %d: Could not open directory <%s>. "
	      "The directory with < %s > version number does not exist. "
	      "Run make_hkpdf.pl script to create directory and config files.\n",
              __FILE__,__LINE__,dirname, f_version_number);
    return NULL;
  }
  /* read each entry until NULL.*/
  top = NULL;
  while( (dir_entry_p = readdir(dir_p)) != NULL ) 
  {
    if( strncmp(dir_entry_p->d_name,"apid",4) )
      continue; /* not an APID config file - skip*/
    if( top == NULL )   
      top = p = malloc(sizeof(APID_HKPFD_Files));
    else                
    {
      p->next = malloc(sizeof(APID_HKPFD_Files));
      p = p->next;
    } 
    p->next = NULL;
    /* load dir and filename in structure link list */
    strcpy(p->filename, dir_entry_p->d_name);
    strcpy(p->directory_name, dirname);
    /* parse filename to get apid and version number*/
    sscanf(p->filename, "%*4s-%x-%*7s-%d",
	   &p->apid, p->version_number); 
    file_loaded_flag=1;
  }
  closedir(dir_p);
  /* check if at least one file loaded and exists in directory*/
  if ( !file_loaded_flag) 
  {
    printkerr("Error at %s, line %d: Could not find config file(s) "
	      "in directory < %s >. Check if file(s) exist, if "
	      "don't exist, then run make_hkpdf.pl script to "
              "create config files.\n", __FILE__,__LINE__,dirname);
  }
  return top;
}/* END-Module: load_hk_filenames_from_directory  */

Keyword_Parameter* create_hdpf_keyword_nodes(HK_Config_Files *ptr_hk_config,
                                             int number_of_lines)   
{
  /* declarations */
  Keyword_Parameter *p;
  int k;

  /* initialize variables */
  p = ptr_hk_config->keywords = NULL;

  for( k=0; k<number_of_lines;  k++)
  {
    if(ptr_hk_config->keywords == NULL) 
    {
      /* Create first node for keyword */
      p = ptr_hk_config->keywords = malloc( sizeof(Keyword_Parameter) );
    } 
    else     
    {
      /*Create all other keyword_parameter nodes required */
      p->next = malloc(sizeof(Keyword_Parameter ));
      p = p->next;
    } 
    p->next = NULL;
  }
  return (ptr_hk_config->keywords);
}

DSC_Conversion* create_hdpf_dsc_nodes(Keyword_Parameter *ptr_hk_keyword)
{
  /* declarations */
  DSC_Conversion *p;
  int k;

  /* initializations */
  p = ptr_hk_keyword->dsc_ptr ;

  if ( p == NULL)
  {
    /*create first node for dsc data */
    p =  (DSC_Conversion *) malloc(sizeof(DSC_Conversion)); 
    ptr_hk_keyword->dsc_ptr= p; /* add connection 5-17-2006*/
  }
  else
  {
    while ( p->next  != (DSC_Conversion *)(NULL))
    {
      p= p->next;
    }
    /*create all other DSC Conversion nodes */
    p->next = (DSC_Conversion *) malloc (sizeof(DSC_Conversion));
    p=p->next;
  }
  p->next= NULL;
  return(p);
}/*end-create_hdpf_dsc_nodes module */

int load_hdpf_keyword_lines(char keyword_lines[MAX_NUM_KW_LINES][MAXLINE_IN_FILE],int number_of_lines,HK_Config_Files *ptr_hk_config_node) 
{
  /* declarations */
  Keyword_Parameter *keyword_node;
  int k=0;

  /*create space for keyword_parameter link list nodes */
  keyword_node= create_hdpf_keyword_nodes( ptr_hk_config_node, number_of_lines );
  /* load keyword values in Keyword_Parameter nodes */
  for(k=0;  keyword_lines[k][0] != '\0'; 
      k++,keyword_node=keyword_node->next) 
  {
    sscanf(keyword_lines[k]," %*s %s %s %d %d %d %s %c %*s",
	   keyword_node->keyword_name,
	   keyword_node->telemetry_mnemonic_name,
	   &(keyword_node->start_byte),
	   &(keyword_node->start_bit_number),
	   &(keyword_node->bit_length),
	   keyword_node->type,
	   &keyword_node->conv_type);
    keyword_node->dsc_ptr = (DSC_Conversion *)NULL;
    keyword_node->alg_ptr = (ALG_Conversion *)NULL;
  }
  return SUCCESSFUL;
} /*END-module: load_hdpf_keyword_lines */

int load_hdpf_dsc_lines(char dsc_lines[MAX_NUM_DCON_LINES][MAXLINE_DCON_IN_FILE ], int number_of_lines, HK_Config_Files *ptr_hk_config_node) 
{
  /* declarations */
  DSC_Conversion *dsc_node;
  int i=0;
  int k=0;
  char tlm_name[HK_MAX_TLM_NAME];
  Keyword_Parameter *kw;

  /* loop throught each of the telemetry names in keywords node link list 
     with conv type equal to D. Search for line in array with telemetry 
     name equal to keyword node in link list. When found, create DSC node 
     and load values from line. 
   */
  for ( kw= ptr_hk_config_node->keywords ; kw; kw=kw->next)
  {
    /* find D type KWDs */
    if ( kw->conv_type == 'D')
    {
      kw->dsc_ptr=(DSC_Conversion*)NULL;
      /* find line with tlm name the same as Keywords tlm name */
      for (i=0; i < number_of_lines; i++)
      {
        sscanf ( dsc_lines[i], "%*s  %s  %*d  %*d  %*d  %*d  %*d  %*d", tlm_name);
        strcat( tlm_name, "");
        if( !strcmp( kw->telemetry_mnemonic_name, tlm_name))
        {
          /* create dsc node memory */
          dsc_node  = create_hdpf_dsc_nodes( kw);
          /* load values in dsc node*/
          sscanf(dsc_lines[i]," %*s %*s %d %d %s %*s", 
            &(dsc_node->low_range), &(dsc_node->high_range),
            dsc_node->dsc_value); 
          dsc_node->next = NULL; /*set link list to null at end node*/
        } /* end of if tlm name equal */
      } /* end of for number of lines */
    } /* end of if D type keyword */
  } /* end of for each keyword */
  return 0;
} /*END-module: load_hdpf_dsc_lines */

int load_hdpf_alg_lines(char alg_lines[MAX_NUM_ACON_LINES][MAXLINE_ACON_IN_FILE],int number_of_lines,HK_Config_Files *ptr_hk_config_node) 
{
  /* declarations */
  ALG_Conversion *alg_node;
  int i, k;
  char tlm_name[HK_MAX_TLM_NAME];
  Keyword_Parameter *kw;

  /* Loop throught each of the telemetry names in keywords node link list  
     with conv type equal to A.Search for line in array with telemetry 
     name equal to keyword node in link list. When found, create ALG node
     and load values from liner.
    */
  for ( kw= ptr_hk_config_node->keywords ; kw; kw=kw->next)
  {
    /* find "A" type KWDs */
    if ( kw->conv_type == 'A')
    {
      /* find line with tlm name the same as Keywords tlm name */
      for (i=0; i < number_of_lines; i++)
      {
        /* get tlm name from line */
        sscanf ( alg_lines[i], "%*s  %s  %*d  %*d  %*d  %*d  %*d  %*d  %*d", tlm_name);
        strcat( tlm_name, "");
        if( !strncmp( kw->telemetry_mnemonic_name, tlm_name, sizeof(tlm_name)))
        {
          /* create alg node memory */
          kw->alg_ptr = (ALG_Conversion*)malloc (sizeof (ALG_Conversion));
          alg_node= kw->alg_ptr;

          /* load values in alg node-assume number of coeffs to be 5*/
          k=0;
          sscanf(alg_lines[i]," %*s %*s %d %lf %lf %lf %lf %lf %*lf %*s",
            &(alg_node->number_of_coeffs), &(alg_node->coeff[k++]),
            &(alg_node->coeff[k++]), &(alg_node->coeff[k++]),
            &(alg_node->coeff[k++]), &(alg_node->coeff[k]) );
          break;
        } /* end of if tlm name equal */
      } /* end of for number of lines */
    } /* end of if "A" type keyword */
  } /* end of for each keyword */
  return SUCCESSFUL;
} /*END-module: load_hdpf_alg_lines */

int save_hdpf_new_formats(FILE* file_ptr,APID_Pointer_HK_Configs *p_apid_ptr_hk_configs) 
{
  /* declarations */
  int err_status[3];
  HK_Config_Files *ptr_config_node;
  HK_Config_Files *previous_ptr_config_node;
  int status;
  int i, j, k;
  char line[MAXLINE_IN_FILE];
  char keyword_lines[MAX_NUM_KW_LINES][MAXLINE_IN_FILE];
  char acon_lines[MAX_NUM_ACON_LINES][MAXLINE_ACON_IN_FILE];
  char dcon_lines[MAX_NUM_DCON_LINES][MAXLINE_DCON_IN_FILE];

  /*Assume format contents of file using functional spec formats */
  /* Check if first node exists*/
  if (p_apid_ptr_hk_configs->ptr_hk_configs == NULL)   
  {
    /* create first HK Config node */
    p_apid_ptr_hk_configs->ptr_hk_configs =malloc( sizeof(HK_Config_Files));
    ptr_config_node = p_apid_ptr_hk_configs->ptr_hk_configs;
    ptr_config_node->next = NULL; 
    ptr_config_node->keywords = NULL;
  } /*END-if first node exists */
  else    
  {
    /* search for new null node in HK_Config Link list */
    for(ptr_config_node = p_apid_ptr_hk_configs->ptr_hk_configs,
	previous_ptr_config_node = ptr_config_node;
	ptr_config_node;
	previous_ptr_config_node = ptr_config_node,
	ptr_config_node = ptr_config_node->next) 
    {
      ;/* continue loop*/
    }         
    /* create next node */
    ptr_config_node = malloc( sizeof(HK_Config_Files));
    previous_ptr_config_node->next = ptr_config_node;
    ptr_config_node->next = NULL;
    ptr_config_node->keywords = NULL;
  }
  /* load keyword, dsc, or alg lines in hk configuration structures */
  for(i=0, j=0, k=0; fgets(line, MAXLINE_IN_FILE, file_ptr) != NULL; ) 
  {
    if (!strncmp(line, "#", 1))
    {
      continue; /*skip comments */
    }
    else if (!strncmp(line, "KWD", 3))
    {
      strcpy( keyword_lines[i++], line);
    }
    else if (!strncmp(line, "DCON", 4))
    {
      
      if(j >= MAX_NUM_DCON_LINES)
      {
        printkerr("WARNING: Maximum lines for array exceeded.\n"
                  "         Skipping saving DCON config data since array too small.\n"
                  "         Update define variable MAX_NUM_DCON_LINES. \n");
      }
      else
      {
        strncpy( dcon_lines[j++], line, MAXLINE_DCON_IN_FILE);
      }
    }
    else if (! strncmp(line, "ACON", 4))
    {
      if(k >= MAX_NUM_ACON_LINES)
      {
        printkerr("WARNING: Maximum lines for array exceeded.\n"
                  "         Skipping saving ACON config data since array too small.\n"
                  "         Update define variable MAX_NUM_ACON_LINES. \n");
      }
      else
      {
        strcpy( acon_lines[k++], line);
      }
    }
    else if (!strncmp(line, "FILE", 4))
    {
      sscanf(line," %*s %*s  %s ", ptr_config_node->file_version_number);
    }
    else if (!strncmp(line, "APID", 4))
    {
      if( strstr(line, "HMI")  )
      {
        sscanf( line,"%*s %x %d %s %s", 
                &(ptr_config_node->apid_number), 
                &(ptr_config_node->number_bytes_used), 
                ptr_config_node->packet_id_type, 
                ptr_config_node->date);
      }
      else if( strstr(line, "AIA") )
      {
        sscanf( line,"%*s %x %d %s %s", 
                &(ptr_config_node->apid_number), 
                &(ptr_config_node->number_bytes_used), 
                ptr_config_node->packet_id_type, 
                ptr_config_node->date);
      }
      else if( strstr(line, "SSIM") )
      {
        sscanf( line,"%*s %x %d %*s %s", 
                &(ptr_config_node->apid_number), 
                &(ptr_config_node->number_bytes_used), 
                ptr_config_node->date);
        /* set to default HMI */
        strcpy(ptr_config_node->packet_id_type,"HMI");
      }
      else
      {
        /* Backward compatible case for apid-#-version-# created*
         * using older(before 6-30-2006) make_hpf.pl script     */ 
        sscanf( line,"%*s %x %d %s", 
                &(ptr_config_node->apid_number), 
                &(ptr_config_node->number_bytes_used), 
                ptr_config_node->date);
        /* set to default HMI */
        strcpy(ptr_config_node->packet_id_type,"HMI");
      }
    }
    else
    {
      printkerr("WARNING: Could not find line with this keyword in apid-#-version-# file\n");
      printkerr("WARNING: Line tried to parse is < %s >\n", line);
    }
  }
  /* Set to null the end of array*/
  strcpy(keyword_lines[i],"");
  strcpy(dcon_lines[j],"");
  strcpy(acon_lines[k],"");
  /*Load keyword, alg and dsc lines values in configuration structures */
  err_status[0] = load_hdpf_keyword_lines(keyword_lines,i,ptr_config_node);
  err_status[1] = load_hdpf_dsc_lines(dcon_lines,j,ptr_config_node);  
  err_status[2] = load_hdpf_alg_lines(acon_lines,k,ptr_config_node);  
  /* return values */
  /* return top_hk_config_nodes by setting value if first node*/
  /* return error status when do error checks. Where
     status 0 = no errors and  status 1 equals errors */
  for (int i=0; i < sizeof(err_status); i++)  
  {
    if ( err_status[i] != ERROR_LOADING_HK_DATA )  
      continue;
    else 
      status= err_status[i];
  }/* END-for */
  return status;
}

void load_config_data(APID_HKPFD_Files *hkpfd_files,
		      APID_Pointer_HK_Configs *hk_configs)  
{
  /*declarations*/
  APID_HKPFD_Files* p;
  int apid;
  char filename[MAX_FILE_NAME];
  FILE* file_ptr; 

  /*intialized variables*/
  int err_status=0;

  /* FOR LOOP through all files in directory */
  apid = hk_configs->apid;
  p = hkpfd_files;
  while(p)
  {
    if ( p->apid == apid )
    {
      sprintf(filename,"%s/%s",p->directory_name, p->filename);
      file_ptr = fopen( filename ,"r");
      err_status = save_hdpf_new_formats(file_ptr, hk_configs);
      fclose(file_ptr);
    }
    p = p->next;
  }    
}

void load_gtcids_data( GTCIDS_Version_Number* top_ptr_gtcids_data,
                       APID_Pointer_HK_Configs *top_apid_ptr_hk_configs) 
{
  /* declarations */
  GTCIDS_Version_Number* tmp_ptr_gtcids_data;
  APID_Pointer_HK_Configs *tmp_apid_ptr_hk_configs;
  HK_Config_Files *tmp_ptr_hk_configs;

  /* Initialize data */
  tmp_ptr_gtcids_data= top_ptr_gtcids_data;
  tmp_apid_ptr_hk_configs= top_apid_ptr_hk_configs;

  /* load gtcids data in APID-PTR and HK-Config-Files structures */
  for(tmp_apid_ptr_hk_configs = top_apid_ptr_hk_configs; 
      tmp_apid_ptr_hk_configs ; 
      tmp_apid_ptr_hk_configs= tmp_apid_ptr_hk_configs->next) 
  {

    for( tmp_ptr_hk_configs = tmp_apid_ptr_hk_configs->ptr_hk_configs;
	 tmp_ptr_hk_configs ;
	 tmp_ptr_hk_configs= tmp_ptr_hk_configs->next) 
    {
      for (tmp_ptr_gtcids_data=top_ptr_gtcids_data;tmp_ptr_gtcids_data ;
	   tmp_ptr_gtcids_data=tmp_ptr_gtcids_data->next)    
      {
	if( !strcmp(tmp_ptr_gtcids_data->file_version_number,
		    tmp_ptr_hk_configs->file_version_number)) 
	{
          if(!strcmp(tmp_ptr_hk_configs->packet_id_type, HMI_ID_TYPE))
          {
	     strcpy(tmp_ptr_hk_configs->parameter_version_number,
		 tmp_ptr_gtcids_data->hmi_id_version_number);
          }
          else if(!strcmp(tmp_ptr_hk_configs->packet_id_type, AIA_ID_TYPE))
          {
	     strcpy(tmp_ptr_hk_configs->parameter_version_number,
		 tmp_ptr_gtcids_data->aia_id_version_number);
          }
          else
          {  /*default set to HMI type */
	     strcpy(tmp_ptr_hk_configs->parameter_version_number,
		 tmp_ptr_gtcids_data->hmi_id_version_number);
          }
          break;
	}
      } /* End-for thru gtcid data */
    } /*End-for loop thru hk config nodes */
  } /*End-for loop thru apid-ptr */
} /* END-module : load_gtcids_data*/ 

void deallocate_apid_hkpfd_files_nodes(APID_HKPFD_Files *ptr) 
{
  /*declarations*/
  APID_HKPFD_Files *tmp;

  while(ptr)
  {
    tmp = ptr->next;
    free(ptr);
    ptr = tmp;
  }
  return;
}

int load_all_apids_hk_configs(char version_number[])    
{
  /*  declarations  */
  APID_Pointer_HK_Configs *apid_ptr_hk_configs;
  APID_Pointer_HK_Configs *top_apid_ptr_hk_configs;
  APID_HKPFD_Files *top_apid_hkpfd_files; 
  APID_HKPFD_Files *hkpfd_files;
  GTCIDS_Version_Number *top_ptr_gtcids_data;
  int apid_array[MAX_APID_POINTERS];
  int number_of_apids, i;
  char file_version_number[50];
  char *ptr_fvn;
  HK_Config_Files *hk_configs;  /* used to print results*/
  Keyword_Parameter *hk_keyword;/* used to print results*/

  /* Init parameter */
  ptr_fvn=file_version_number;
  /* Load data from ground to code ids file */
  top_ptr_gtcids_data = read_gtcids_hk_file(top_apid_ptr_hk_configs); 
  /* Check for file version number in ground to code ids file */
  ptr_fvn=find_file_version_number(top_ptr_gtcids_data, version_number);
  /* load HK_Config_Files structures for each apid */
  if ((top_apid_hkpfd_files = read_all_hk_config_files(ptr_fvn)) == NULL)
    return ERROR_HK_NOSUCHDIR;
  /* get list of unique apids to read and allocate space for*/
  number_of_apids = 0;
  memset(apid_array, 0, sizeof(apid_array));
  for(hkpfd_files=top_apid_hkpfd_files; hkpfd_files; 
      hkpfd_files=hkpfd_files->next)
  {
    /* See if we already have this apid. */
    for (i=0; i<number_of_apids; i++)
    {
      if (apid_array[i] == hkpfd_files->apid)
        break;
    }
    if (i >= number_of_apids) \
    {
      /* Nope this was a new one. Insert it in the list. */
      apid_array[i]= hkpfd_files->apid;
      number_of_apids++;
    }
  }
  /* load apid pointer to HK Configs structure apid-value and allocate nodes */  
  for ( i=0; i < number_of_apids; i++)     
  {
    if ( i == 0)    
    { 
      /* create and load first node */
      apid_ptr_hk_configs = (APID_Pointer_HK_Configs *)malloc(sizeof(APID_Pointer_HK_Configs));
      top_apid_ptr_hk_configs = apid_ptr_hk_configs;
    }
    else    
    { 
      /* create and load next node */
      apid_ptr_hk_configs->next = (APID_Pointer_HK_Configs *)malloc(sizeof(APID_Pointer_HK_Configs));
      apid_ptr_hk_configs = apid_ptr_hk_configs->next;      
    } 
    apid_ptr_hk_configs->apid = apid_array[i];
    apid_ptr_hk_configs->next= NULL;
    apid_ptr_hk_configs->ptr_hk_configs= (HK_Config_Files*)NULL;
  }             
  /* set moving pointer to top of list */
  apid_ptr_hk_configs = top_apid_ptr_hk_configs;
  for(i=0; i < number_of_apids; i++)  
  {
    load_config_data( top_apid_hkpfd_files, apid_ptr_hk_configs);
    apid_ptr_hk_configs = apid_ptr_hk_configs->next;
  }  /*End for*/
  /* load file_version_number based on packet version number */
  load_gtcids_data(top_ptr_gtcids_data, top_apid_ptr_hk_configs);
  /* set global variable to use for decode_hk modules */
  global_apid_configs = top_apid_ptr_hk_configs;
  /* deallocate  top_apid_hkpfd_files link list */
  deallocate_apid_hkpfd_files_nodes(top_apid_hkpfd_files) ;
/******PRINT RESULTS **********/
/***
for (apid_ptr_hk_configs=global_apid_configs; apid_ptr_hk_configs;
     apid_ptr_hk_configs=apid_ptr_hk_configs->next)
{
   if( apid_ptr_hk_configs->apid == 445)
   {
      printf("APID_Pointer_HK_Configs apid is %d\n", apid_ptr_hk_configs->apid);
      for (hk_configs=apid_ptr_hk_configs->ptr_hk_configs; hk_configs;
            hk_configs=hk_configs->next)
      {
          printf("HK_Configs file version number  is %s\n", hk_configs->file_version_number );
          if (!strcmp(hk_configs->file_version_number, "1.38"))
          {
            printf("HK_Configs file version number  is %s\n", hk_configs->file_version_number );
            for ( hk_keyword = hk_configs->keywords;hk_keyword;  
                  hk_keyword = hk_keyword->next)
            {
              printf("Keywords  is %s\n", hk_configs->keywords->keyword_name );
            }
          }
      }
   }
}
***/
/******END PRINT RESULTS*******/
  return SUCCESS;  
}/*End Module: load_all_apids_hk_configs */

HK_Config_Files* check_packet_version_number( HK_Config_Files *ptr_to_configs,
				       char *version_number )       
{
  while(ptr_to_configs != NULL )  
  {
    if(!strcmp(ptr_to_configs->parameter_version_number, version_number))
    {
      break; 
    }
    ptr_to_configs = ptr_to_configs->next;
  }
  return ptr_to_configs;
}

int check_free_configs_flag(void)
{
  /*declarations */
  int status;
  FILE *file_ptr;
  char *directory ;
  char *filename;
  char directory_filename[MAX_FILE_NAME];

  /* get directory and file name */
  directory = getenv("HK_CONFIG_DIRECTORY");
  filename= getenv("HK_FREE_FLAG_FILE");
  if(filename == NULL) 
    filename = "HK_FREE_CONFIG_FLAG_FILE";
  if(directory == NULL) 
    directory = "../../tables/hk_config_file";
  sprintf(directory_filename, "%s/%s", directory, filename); 
  /*open file-if exists then clear configs each time*/
  file_ptr = fopen(directory_filename, "r");
  if (file_ptr == NULL)
  {
    status =0;
  }
  else 
  {
    status = 1;
    fclose(file_ptr);
  }
  return status;
}

void deallocate_apid_ptr_hk_config_nodes(void)     
{
  extern APID_Pointer_HK_Configs *global_apid_configs;
  APID_Pointer_HK_Configs *tmp_apid_ptr_hk_configs;
  APID_Pointer_HK_Configs *prev_apid_ptr_hk_configs;
  HK_Config_Files *tmp_ptr_hk_configs;
  HK_Config_Files *prev_ptr_hk_configs;
  Keyword_Parameter *tmp_keyword_node;
  Keyword_Parameter *prev_keyword_node;
  DSC_Conversion *tmp_dsc_node;
  DSC_Conversion *prev_dsc_node;
  ALG_Conversion *tmp_alg_node;
  ALG_Conversion *prev_alg_node;
  int free_flag=0;

  /* check if want to free all stored configurations */
  if (!(free_flag= check_free_configs_flag()))
  {
     /*skip deallocating configuration data kept in memory*/
     return;
  }

  /* assign top of link list of config data to tmp ptr */
  tmp_apid_ptr_hk_configs=global_apid_configs;

  /* Loop throught and free nodes */
  for(prev_apid_ptr_hk_configs=tmp_apid_ptr_hk_configs;
      tmp_apid_ptr_hk_configs;
      prev_apid_ptr_hk_configs = tmp_apid_ptr_hk_configs )  
  {
    for(prev_ptr_hk_configs=tmp_ptr_hk_configs=tmp_apid_ptr_hk_configs->ptr_hk_configs;
	tmp_ptr_hk_configs;
	prev_ptr_hk_configs =tmp_ptr_hk_configs)  
    {
      for(prev_keyword_node=tmp_keyword_node=
	    tmp_ptr_hk_configs->keywords;
	  tmp_keyword_node;
	  prev_keyword_node =tmp_keyword_node)  
      {
        /* free dsc nodes */
        for(prev_dsc_node=tmp_dsc_node=tmp_keyword_node->dsc_ptr;
            tmp_dsc_node; prev_dsc_node = tmp_dsc_node)
        {
          tmp_dsc_node= tmp_dsc_node->next;
          free( prev_dsc_node);
        }
        /* free alg node */
        tmp_alg_node=tmp_keyword_node->alg_ptr;
        if (tmp_keyword_node)
        {
          free( tmp_alg_node );
        }
	tmp_keyword_node=tmp_keyword_node->next;
	free(prev_keyword_node);
      }
      tmp_ptr_hk_configs =tmp_ptr_hk_configs->next;
      free(prev_ptr_hk_configs);
    }
    tmp_apid_ptr_hk_configs = tmp_apid_ptr_hk_configs->next ;
    free(prev_apid_ptr_hk_configs);
  }
  //ADDED 6-28 to test deallocate
  global_apid_configs=NULL;

}

HK_Config_Files*  reread_all_files(APID_Pointer_HK_Configs *apid_ptr_configs,char version_number[]) 
{
  /*declarations*/
  int apid;
  APID_Pointer_HK_Configs *p, *prev_p;
  char file_version_number[50];
  char *ptr_fvn;
  int found_flag;
  GTCIDS_Version_Number *ptr_gtcids_data;
  APID_HKPFD_Files *top_hkpfd_files;
  APID_HKPFD_Files *hkpfd_files;

  /* init values */
  ptr_fvn=file_version_number;
  ptr_gtcids_data = global_gtcids_vn;

  /*save apid value to look up later */
  apid= apid_ptr_configs->apid;
  /* find file version number directory to read in files */
  ptr_fvn=find_file_version_number(ptr_gtcids_data, version_number);
  /* load HK_Config_Files structures for each apid */
  if ((top_hkpfd_files = read_all_hk_config_files(ptr_fvn)) == NULL)
    return (HK_Config_Files*)NULL;
  /* check which apid node to add */
  for(hkpfd_files=top_hkpfd_files; hkpfd_files;
      hkpfd_files=hkpfd_files->next)
  {
    found_flag=0;
    for(p=global_apid_configs; p ; p=p->next)
    {
      if ( p->apid == hkpfd_files->apid)
      {
        found_flag=1;
        break;
      }
      prev_p=p;
    }
    if (!found_flag) 
    {
      /* add apid node to link list */
      prev_p->next = (APID_Pointer_HK_Configs *) malloc (sizeof (APID_Pointer_HK_Configs));
      prev_p->next->next= NULL;
      prev_p->next->ptr_hk_configs= (HK_Config_Files*) NULL;
      prev_p->next->apid=hkpfd_files->apid;
    }
  }
  /*set moving pointer to top of list and load config data */
  for(p=global_apid_configs; p ; p=p->next)
  {
    load_config_data(top_hkpfd_files, p);
  }
  /* load data for packet version numbers in hk config format link list*/
  p=global_apid_configs; 
  load_gtcids_data(ptr_gtcids_data, p);
  /* deallocate  top_apid_hkpfd_files link list */
  hkpfd_files=top_hkpfd_files;
  deallocate_apid_hkpfd_files_nodes(hkpfd_files) ;
  /* locate top HK_Config_Files Node for APID_Pointer node equal to lookup apid value */
  p = global_apid_configs;
  while (p)
  {
    if( p->apid == apid )  
    {
      return p->ptr_hk_configs;
    }/* if apid equal p->apid */
    p = p->next; //added 1-19-2006 
  } 
  /* APID not found, return NULL. */
  return (HK_Config_Files*)NULL;
}

int decode_hk_keywords(unsigned short *word_ptr, int apid, HK_Keyword_t **kw_head) 
{
  /* declarations */
  HK_Keywords_Format *ptr_hk_keywords;
  HK_Keywords_Format *top_ptr_hk_keywords;
  APID_Pointer_HK_Configs *apid_configs;
  HK_Config_Files *config_files;
  int status;
  char version_number[MAX_CHAR_VERSION_NUMBER]; 
  HK_Config_Files *matching_config;
  char *init_hdr_apid;
  char *init_hdr_pkt_version;

  /* init values */
  matching_config= (HK_Config_Files *)NULL;

  /* get environmental variables */
  init_hdr_apid= getenv("HK_INITIAL_HEADER_APID");
  init_hdr_pkt_version= getenv("HK_INITIAL_HEADER_PKT_VERSION");
  if( !init_hdr_apid || !init_hdr_pkt_version)
  {
    ERRMSG("Please set. Could not find environment variables:  <HK_INITIAL_HEADER_APID> <HK_INITIAL_HEADER_PKT_VERSION>");
    return HK_DECODER_ERROR_UNKNOWN_ENV_VARIABLE;
  }

  /* Get Version Number  from byte 8 and byte 9  bit stream for hk packet  */
  if (apid == 431) /*inital packet header is 8 bytes and called apid 431 */ 
  {
    strcpy(version_number,init_hdr_pkt_version);/*parameter version # from packet*/ 
  }
  else
  {
    get_version_number(word_ptr, version_number); 
  }

  /* check if global pointer to configuration link list exists*/
  if (!global_apid_configs )
  {
    load_all_apids_hk_configs(version_number);
  }

  /* get pointer to HK_Config_File structure based on apid # */
  apid_configs = global_apid_configs;
  while(apid_configs)   
  {
    if (apid_configs->apid == apid)  
      break; /*Found pointer apid pointer correct HK_Config_File structures */
    else
      apid_configs = apid_configs->next;
  }

  if ( apid_configs == NULL) 
  {
  
    printkerr("ERROR at %s, line %d: This apid <%x> does not have valid "
              "config data to decode packet. Check if config files exist. "
              "If don't exist, run make_hkpdf.pl script to create files.\n",
               __FILE__, __LINE__, apid);
    return HK_DECODER_ERROR_UNKNOWN_APID;
  }

  /*set pointer to configs to tmp pointer to get config data for this apid*/
  config_files = apid_configs->ptr_hk_configs;     

  /**************************************************************
   * Check if version number exists for hk-configurations       *
   * There are two version numbers:                             *
   * Version number for STANFORD_TLM_HMI_AIA file               *
   * Version number for parameter                               *
   **************************************************************/
  matching_config = check_packet_version_number(config_files,version_number);
  if ( matching_config == NULL )
  {
    /* Could not find matching config. Try to reload the config files. */
    /*printkerr("WARNING at %s, line %d: For apid <%x> and version number <%s> "
	      "cannot find config data to decode packet!!. Reread config files"
              " in folder %s.\n", __FILE__, __LINE__,apid, version_number,version_number);
    */
    config_files = reread_all_files(apid_configs, version_number);
    matching_config = check_packet_version_number(config_files,version_number);
    if ( matching_config != NULL ) 
    {
      /*printkerr("WARNING at %s, line %d: Found config for apid <%x> and version number <%s>\n",
                 __FILE__, __LINE__,apid, version_number);
       */;
    }
  }
  if ( matching_config == NULL )
  {
    /* Still no matching config information. Return an error code. */
    ERRMSG("Could not find  version number even after re-reading hk config files.");
    return HK_DECODER_ERROR_CANNOT_FIND_VER_NUM; 
  }
  /* Load hk configs in HK_Keywords_Format structure  */
  if (ptr_hk_keywords = load_hk_configs( matching_config ))
  {
    /* Load hk values in HK_Keyword_struct  */
    if (status = load_hk_values( word_ptr, ptr_hk_keywords))
      return status;
    /* set engr values and raw values in HK_Keyword_t structure which is used for Lev0 Processing*/
    if (!( status = load_engr_values(ptr_hk_keywords, kw_head)))
    {
      /* deallocate HK_Keyword_Format structure */
      deallocate_hk_keywords_format(ptr_hk_keywords);
      status = HK_DECODER_SUCCESSFUL; /*successful!!*/
      deallocate_apid_ptr_hk_config_nodes();/*check setting before doing*/
      return status;
    }
  }
  else  
  {
    ERRMSG("Could not find config data for this version number");
    status = HK_DECODER_ERROR_CANNOT_LOAD_HK_VALUES;
  }
  return status;
}/*Module:decode_hk_keywords*/



/***************************************************************************** 
 * COPY HK Keywords 
 *
 *****************************************************************************/
HK_Keyword_t *copy_hk_keywords(HK_Keyword_t *head)   
{
  HK_Keyword_t *newhead, *p;

  if (head)
  {
    newhead = malloc(sizeof(HK_Keyword_t));
    assert(newhead);
    p = newhead;
    memcpy(p, head, sizeof(HK_Keyword_t));
    if (head->eng_type == DRMS_TYPE_STRING)
      p->eng_value.string_val = strdup(head->eng_value.string_val);
    head = head->next;
    while(head)
    {      
      p->next = malloc(sizeof(HK_Keyword_t));
      assert(p->next);
      p = p->next;
      memcpy(p, head, sizeof(HK_Keyword_t));
      if (head->eng_type == DRMS_TYPE_STRING)
	p->eng_value.string_val = strdup(head->eng_value.string_val);
      head = head->next;      
    }   
  }
  else
    newhead = NULL;

  return newhead;
}

// From decompress.c

CropTable_t croptables[MAX_CROPID+1];
unsigned short *lutables[MAX_LUTID+1];
unsigned short *ilutables[MAX_LUTID+1];
static int table_initialized = 0;
static unsigned char numzero16[1<<16];

static unsigned short *decode_im_pdu(unsigned short *w, IM_PDU_Packet_t *p)
{
    /* Word 0 */
    p->im_pdu_count = ((long long) (*w & 0x3ff)) << 32;
    p->im_pdu_id = *w++ >> 10;
    /* Word 1 */
    p->im_pdu_count |= *w++ << 16;
    /* Word 2 */
    p->im_pdu_count |= *w++;
    /* Word 3 */  
    p->fhp = (*w++) & 0x7ff;
    return w;
}

static unsigned short *decode_ccsds(unsigned short *w, CCSDS_Packet_t *p)
{
    /* Word 4 */  
    p->version = (*w >> 13) & 0x7;
    p->type = (*w >> 12) & 1;
    p->shf =  (*w >> 11) & 1;
    p->apid =  (*w++) & 0x7ff;
    /* Word 5 */  
    p->sf = (*w >> 14) & 0x3;
    p->ssc =  (*w++) & 0x3fff;
    /* Word 6 */  
    /* According to the CCSDS spec the packet length field contains 
       "total number of octets - header octets - 1". 
       Add 1 to get length of data field in bytes. */
    p->length = *w++ + 1;
    return w;
}

/* Extract fields from HMI/AIA science data packet header. */
static unsigned short *decode_scidata(unsigned short *w, SciDataPacket_t *p)
{
  int i;
  /* Word 7 */  
    p->shs = (*w++)<<16;
    /* Word 8 */  
    p->shs |= *w++;
    /* Word 9 */  
    p->shss = (*w++)<<16;
    /* Word 10 */  
    p->shss |= *w++;
    /* Word 11-14 */  
    for (i=0; i<4; i++)
      p->ccdhead[i] = *w++;
    /* Word 15 */  
    p->cropid = *w >> 4;
    p->romode = (*w >> 2) & 0x3;
    p->headererr = (*w >> 1) & 0x1;
    p->oflow = *w++ & 0x1;  
    /* Word 16 */  
    p->tapcode = *w >> 12;
    p->bitselectid = (*w >> 8) & 0xf;
    p->compid = *w++ & 0xff;
    /* Word 17 */  
    p->lutid = *w >> 8;
    p->offset = (*w++ & 0xff) << 16;
    /* Word 18  */  
    p->offset |= *w++;
    
    p->data = w;
    return w;
}

/* 
   Create a new image decompression context corresponding to an image
   that is in the process of being decompressed and assembled from its
   constituent telemetry packets.  An image buffer (of type Image)
   is allocated and a pointer to it will be returned by
   decompress_next_packet when the image is complete. It is the
   responsibility of the caller of decompress_next_packet to free this
   memory.
*/
static int decompress_new_context(unsigned int FSN, unsigned int FID)
{
  Decompress_Context_t *C;
  
  if (numcontexts>=MAXCONTEXTS)
  {
    printkerr("Maximum number of contexts exceeded.\n" 
	      "More than %d images are being decoded simultaneously.\n",
	      MAXCONTEXTS);
    return -1;
  }
  else
  {
    C = malloc(sizeof(Decompress_Context_t));
    memset(C,0,sizeof(Decompress_Context_t));
    C->ID = UNIQUEID(FSN,FID);
    
    /* Init image struct. */
    C->image = malloc(sizeof(Image_t));
    C->image = memset(C->image,0,sizeof(Image_t));
    C->image->ID = C->ID;
    C->image->keywords = NULL;
    C->image->next = NULL;
    /* Init firstpacket struct. cropid = -1 is used to indicates that 
       no image packets have been received yet for this context. 
       When the first packet is received the cropid
       is set to the value in the packet header.
       ASSUMPTION: -1 is not a valid cropid. */
    C->image->firstpacket.cropid = -1; 

    /* Init stat struct. */
    C->image->stat.ID = C->ID;
    C->image->stat.starttime = time(NULL);

    /* Insert new structure in the global context table. */
    Context[numcontexts] = C;
    ++numcontexts;

    return numcontexts-1;
  }
}

/* 
   Return the number and status of images currently being
   decompressed. decompress_status returns the number of images being
   decompressed and on return *stat point to an array of structs of
   type Decompress_Stat (see hmi_compression.h), one for each
   image.
 */
int decompress_status_all(Decompress_Stat_t **stat)
{
  int i;
  if (numcontexts<=0)
  {
    *stat = NULL;
    return 0;
  }
  else
  {
    *stat = malloc(numcontexts*sizeof(Decompress_Stat_t));    
    for (i=0; i<numcontexts; i++)
      (*stat)[i] = Context[i]->image->stat;
    return numcontexts;
  }
}

/* 
   Print status information for a decompression context as returned by 
   decompress_status.
*/
void decompress_print_status( Decompress_Stat_t *stat)
{
  if (stat==NULL)
    return;

  printk("  ID = %llu = 0x%0llx\n (FSN,FID) = (%d, %d) = (0x%0x, 0x%0x).\n",
	 stat->ID,stat->ID,ID2FSN(stat->ID),ID2FID(stat->ID),
	 ID2FSN(stat->ID),ID2FID(stat->ID));
  printk("  %hu packets received.\n", stat->npackets);
  printk("  %u pixels (%d%%) out of %u received.\n",
	 stat->numpix, (100*stat->numpix)/stat->totalpix,
	 stat->totalpix);
  printk("  First packet received at %s",ctime(&(stat->starttime)));
  printk("  backup_occured = %hhu\n",stat->backup_occured);
  printk("  skip_occured = %hhu\n",stat->skip_occured);
}

/* Find the decompression context corresponding to
   a given unique ID. */
static int decompress_findcontext(unsigned int FSN, unsigned int FID)
{
  int ctx;
  long long ID;

  ID = UNIQUEID(FSN,FID);
  for (ctx=0; ctx<numcontexts; ctx++)
    if (Context[ctx] != NULL && Context[ctx]->ID == ID)
      break;
  return ctx;
}

/* Undo pixel value transformation based on lookup table. */
void decompress_undotransform(const unsigned int N, const unsigned int R, 
			      const unsigned int ILUTID, 
			      const unsigned int numpix, unsigned short *pixels)
{
  unsigned int i;
  unsigned short nmask;
  unsigned short *ILUT = ilutables[ILUTID];

  assert(N>=NMIN && N<=NMAX);
  assert(R<=RMAX);
  nmask = (1<<N)-1;

  if (ILUTID==0) {
      for (i=0; i<numpix; i++)
	  if (pixels[i] != 0xffff)
	      pixels[i] = (pixels[i] & nmask) << R;
	  else
	      pixels[i] = 0x8000; // MISSING = -32768
  } else {
      for (i=0; i<numpix; i++)
	if (pixels[i] != 0xffff) 
	    pixels[i] = (ILUT[pixels[i]] & nmask) << R;
	else
	    pixels[i] = 0x8000; // MISSING = -32768
  }
}

/* Insert pixel values from a 1-D pixel buffer into an image according 
   the a crop table. */
void decompress_uncrop(const unsigned int CROPID, const unsigned int numpix, 
		       unsigned short *pixels, Image_t *image)
{
  unsigned int row, skip, take, total, width, height;
  unsigned short *table;
  short *data, *t;
  int tapcode;
  int i,j;
  
  width = croptables[CROPID].width;
  height = croptables[CROPID].height;
  table = croptables[CROPID].table;
  data = image->data;
  tapcode = image->firstpacket.tapcode;

  if (!CROPID) {	// no crop
      memcpy(data, pixels, 2*numpix);
  } else {
      total = 0;
      for (row=0; row<height && total<numpix; row++)
      {
	skip = table[2*row];
	take = table[2*row+1];
	if (total+take>numpix)
	  take = numpix-total; /* Incomplete image. */
	memcpy(&data[row*width + skip], &pixels[total],
	       take*sizeof(unsigned short));  
	total += take;
      }
  }

  //  H--G
  //  |  |
  //  E--F
  width = image->width;
  height = image->height;
  switch (tapcode) {
  case 0:	// 4-port, corners E,F,G,H
      t = malloc(2*width*height);
      for (i=0; i<height/2; ++i) {
	  memcpy(t+i*width, data+i*width/2, width);
	  memcpy(t+(height-1-i)*width, data+(height+height/2+i)*width/2, width);
	  for (j=0; j<width/2; ++j) {
	      t[width-1-j+i*width] = data[j+(height/2+i)*width/2];
	      t[width-1-j+(height-1-i)*width] = data[j+(height+i)*width/2];
	  }
      }
      free(data);
      image->data = t;
      break;
  case 1:	// 2-port, corners E,F
      t = malloc(2*width*height);
      for (i=0; i<height; ++i) {
	  memcpy(t+i*width, data+i*width/2, width);
	  for (j=0; j<width/2; ++j)
	      t[width-1-j+i*width] = data[j+(height+i)*width/2];
      }
      free(data);
      image->data = t;
      break;
  case 2:	// 2-port, corners F,G
      t = malloc(2*width*height);
      for (i=0; i<height/2; ++i)
	  for (j=0; j<width; ++j) {
	      t[j+i*width] = data[width-1-j+i*width];
	      t[j+(height-1-i)*width] = data[width-1-j+(i+height/2)*width];
	  }
      free(data);
      image->data = t;
      break;
  case 3:	// 2-port, corners G,H
      t = malloc(2*width*height);
      for (i=0; i<height; ++i) {
	  memcpy(t+(height-1-i)*width, data+(i+height)*width/2, width);
	  for (j=0; j<width/2; ++j)
	      t[width-1-j+(height-1-i)*width] = data[j+i*width/2];
      }
      free(data);
      image->data = t;
      break;
  case 4:	// 2-port, corners E,H (counterintuitive!)
      t = malloc(2*width*height);
      for (i=0; i<height/2; ++i) {
	  memcpy(t+(height-1-i)*width, data+(i+height/2)*width, 2*width);
	  memcpy(t+i*width, data+i*width, 2*width);
      }
      free(data);
      image->data = t;
      break;
  case 5:	// 1-port, corner E
      // nothing to do
      break;
  case 6:	// 1-port, corner F
      t = malloc(2*width*height);
      for (i=0; i<height; ++i)
	  for (j=0; j<width; ++j)
	      t[j+i*width] = data[width-1-j+i*width];
      free(data);
      image->data = t;
      break;
  case 7:	// 1-port, corner G
      t = malloc(2*width*height);
      for (i=0; i<height; ++i)
	  for (j=0; j<width; ++j)
	      t[j+i*width] = data[width-1-j+(height-1-i)*width];
      free(data);
      image->data = t;
      break;
  case 8:	// 1-port, corner H
      t = malloc(2*width*height);
      for (i=0; i<height; ++i)
	  memcpy(t+i*width, data+(height-1-i)*width, 2*width);
      free(data);
      image->data = t;
      break;
  }
}

/* Free a context structure. If discard_image is 0 then do
   not free the image data associated with the context, but
   pass it back to the caller after freeing everything else. 
*/
Image_t *decompress_free_context(unsigned int ctx, int discard_image)
{
  Image_t *im;
  Decompress_Context_t *C;
  HK_Keyword_t *kw, *tmp;


  C = Context[ctx];
  Context[ctx] = Context[numcontexts-1];
  --numcontexts;
  Context[numcontexts] = NULL;

  free(C->pixelbuf);
  if (discard_image)
  {
    deallocate_hk_keyword(C->image->keywords);
    free(C->image->data);
    free(C->image);    
    im = NULL;
  }
  else
    im = C->image;
  free(C);

  return im;
}

/* 
   Prepare the image in the context no. ctx for output
   and free the associated context. 
   On input the image pixels are stored linearly in the
   context pixelbuffer. The preparation process involves 
   3 steps:

   1. Undo table lookup transformation of pixels values.
   2. Undo cropping, i.e. put the pixels into their proper
      locations in the 2D image array according to the crop table.
   3. Free context and the associated pixel buffer.
*/

static Image_t *decompress_prepare_image(int ctx)
{
  Image_t *im;
  unsigned int N, R;
  SciDataPacket_t *scipack;
  Decompress_Stat_t *stat;
  
  im = Context[ctx]->image; 
  scipack = &(im->firstpacket);
  stat = &im->stat;
  N = scipack->compid >> 3;
  if (N==0) N=16;		// uncompressed mode
  R = scipack->bitselectid;
  
  decompress_undotransform(N,R,scipack->lutid,stat->totalpix,Context[ctx]->pixelbuf);    
  decompress_uncrop(scipack->cropid,stat->totalpix,Context[ctx]->pixelbuf, im);
  
  // Get status information for this image.
  return decompress_free_context(ctx, 0);
}


/* 
   Force reconstruction of a partially received image. If a
   decompression context with the corresponding FSN and FID exists,
   the partially decompressed image it contains is untransformed,
   uncropped and returned in *image, and a value of SUCCESS (0) is
   returned. If the given FSN and FID does not match those of an image
   currently being reconstructed, a value of ERROR_NOSUCHIMAGE (-7) is
   returned.
 */
int decompress_flush_image(unsigned int FSN, unsigned int FID, 
			   Image_t **image)
{
  int ctx;  

  ctx = decompress_findcontext(FSN,FID);
  if (ctx>=numcontexts)
  {
    *image = NULL;
    return ERROR_NOSUCHIMAGE;
  }
  if (!Context[ctx]->pixelbuf) {
      *image = NULL;
      return ERROR_EMPTY_IMAGE_CONTEXT;
  }
  *image = decompress_prepare_image(ctx);
  return SUCCESS;
}

/* Write the image as a FITS file with the telemetry header values
   as FITS header keywords. */
int decompress_writefitsimage(const char* file, Image_t *image, int compress)
{
  int n_kw, hlen;
  char *p, *head;
  SciDataPacket_t *h;
  HK_Keyword_t *kw;
  int npkts, missvals;

  missvals = image->stat.totalpix - image->stat.numpix;
  npkts = image->stat.npackets;

  kw = image->keywords;
  n_kw = 0;
  while(kw)
  {
    n_kw++;
    kw = kw->next;
  }
  hlen = ((n_kw+15+35)/36)*36*80;
  head = malloc(hlen);
  assert(head);
  memset(head,' ',hlen);
  p = head;
  /* Write special filtergram keywords. */
  h = &image->firstpacket;
  snprintf(p,80,"SHS     = %u",h->shs);
  p += 80;
  snprintf(p,80,"SHSS    = %u",h->shss);
  p += 80;
  snprintf(p,80,"FSN     = %u",ID2FSN(image->ID));
  p += 80;
  // FID in scicence packet is not reliable, so hide it
  //snprintf(p,80,"FID     = %u",ID2FID(image->ID));
  //p += 80;
  snprintf(p,80,"CROPID  = %d",h->cropid);
  p += 80;
  snprintf(p,80,"ROMODE  = %d",h->romode);
  p += 80;
  snprintf(p,80,"HEADRERR= %d",h->headererr);
  p += 80;
  snprintf(p,80,"OVERFLOW= %d",h->oflow);
  p += 80;
  snprintf(p,80,"TAPCODE = %d",h->tapcode);
  p += 80;
  snprintf(p,80,"BITSELID= %d",h->bitselectid);
  p += 80;
  snprintf(p,80,"LUTID   = %d",h->lutid);
  p += 80;
  snprintf(p,80,"COMPID  = %d",h->compid);
  p += 80;
  snprintf(p,80,"OFFSET  = %d",h->offset);
  p += 80;
  snprintf(p,80,"MISSVALS= %d",missvals);
  p += 80;
  snprintf(p,80,"NUMPKTS = %d",npkts);
  p += 80;
  kw = image->keywords;
  while(kw)
  {
    if(kw->eng_type == DRMS_TYPE_STRING)
    {
      snprintf(p,80,"%-8s= %s",kw->fitsname, kw->eng_value.string_val);
    }            
    else  if(kw->eng_type == DRMS_TYPE_UINT8)
    {
      snprintf(p,80,"%-8s= %"PRIu8,kw->fitsname, kw->eng_value.uint8_val);
    }
    else if (kw->eng_type == DRMS_TYPE_UINT16)
    {
      snprintf(p,80,"%-8s= %"PRIu16,kw->fitsname, kw->eng_value.uint16_val);
    }
    else if (kw->eng_type == DRMS_TYPE_UINT32)
    {
      snprintf(p,80,"%-8s= %"PRIu32,kw->fitsname, kw->eng_value.uint32_val);
    }
    else if (kw->eng_type == DRMS_TYPE_DOUBLE)
    {
      snprintf(p,80,"%-8s= %lf",kw->fitsname, kw->eng_value.double_val);
    }
    else if (kw->eng_type == DRMS_TYPE_INT8)
    {
      snprintf(p,80,"%-8s= %"PRId8,kw->fitsname, kw->eng_value.int8_val);
    }
    else if (kw->eng_type == DRMS_TYPE_INT16)
    {
      snprintf(p,80,"%-8s= %"PRId16,kw->fitsname, kw->eng_value.int16_val);
    }
    else if (kw->eng_type == DRMS_TYPE_INT32)
    {
      snprintf(p,80,"%-8s= %"PRId32,kw->fitsname, kw->eng_value.int32_val);
    }
    else
    {
       printkerr("Unexpected execution of else. All keywords should have type values. This one does not < %s>\n",kw->name);
    }
    //rasmus's old code-snprintf(p,80,"%-8s= %"PRId64,kw->fitsname, kw->raw_value);
    kw = kw->next;
    p += 80;    
  }
  

  if (compress<=1)
    return fits_writeimage_raw(file, compress, (int) (p-head), head, 
			       image->width, image->height, image->data);
  else 
  {
    printkerr("Invalid value (%d) for compress.\n",compress);
    return 1;
  }
}

void decompress_free_hk(CCSDS_Packet_t *p)
{
  CCSDS_Packet_t *tmp;
  while(p)
  {
    if (p->keywords)
      deallocate_hk_keyword(p->keywords);
    tmp = p->next;
    free(p);
    p = tmp;
  }
}

void decompress_free_images(Image_t *image)
{
  Image_t *ip;
  HK_Keyword_t *kw, *ktmp;
  while (image)
  {
    deallocate_hk_keyword(image->keywords);
    free(image->data);
    ip = image->next;
    free(image);
    image = ip;
  }
}

static void decompress_initnonzerotables(void)
{
  unsigned int i,nz,b;
    
  for (i=0; i<=0xffff; i++)
  {
    b = i;
    nz = 0; 
    while ((b&1) == 0 && nz<16)
    {
      ++nz;
      b = b>>1;
    }
    numzero16[i] = nz;
  }
  table_initialized = 1;
}

static void decompress_initcroptables(void)
{
  unsigned int i;
    
  // TODO: Must initialize lookup- and crop tables properly here.
  for (i=0; i<MAX_CROPID; i++)
  {
    croptables[i].width = 0;
    croptables[i].height = 0;
    croptables[i].table = NULL;
  }
#if 0
  // Set up simple builtin tables.
  assert(MAX_CROPID>=5);
  decompress_square_croptable(&croptables[0],256);
  decompress_square_croptable(&croptables[1],1024);
  decompress_square_croptable(&croptables[2],4096);
  decompress_circular_croptable(&croptables[3],256);
  decompress_circular_croptable(&croptables[4],1024);
  decompress_circular_croptable(&croptables[5],4096);
#endif
}

static void decompress_initlookuptables(void)
{
  unsigned int i;
  unsigned short maxval;
  double x, finv;
    
  // TODO: Must initialize lookup- and crop tables properly here.
  for (i=0; i<MAX_LUTID; i++)
  {    
    lutables[i] = NULL;
    ilutables[i] = NULL;
  }
#if 0
  assert(MAX_LUTID>=2);
  maxval = (1<<NMAX)-1;
  lutables[0] = malloc((maxval+1)*sizeof(short));  
  lutables[1] = malloc((maxval+1)*sizeof(short));
  lutables[2] = malloc((maxval+1)*sizeof(short));
  ilutables[0] = malloc((maxval+1)*sizeof(short));  
  ilutables[1] = malloc((maxval+1)*sizeof(short));
  ilutables[2] = malloc((maxval+1)*sizeof(short));
  for (i=0, x=0; i<(1<<NMAX); i++, x += 1.0)
  {
    
    lutables[0][i]  = i;
    ilutables[0][i] = i;

    lutables[1][i]  = floor(sqrt(128*x)+0.5);
    finv = floor((x*x)/128+0.5);
    ilutables[1][i] = finv > (float)maxval ? maxval : (unsigned short) finv;

    lutables[2][i]  = floor(sqrt(1024*x)+0.5);
    finv = floor((x*x)/1024+0.5);
    ilutables[2][i] = finv > (float)maxval ? maxval : (unsigned short) finv;
  }  
#endif
}

void decompress_inittables(void)
{
  decompress_initnonzerotables();
  decompress_initcroptables();
  decompress_initlookuptables();
}

int decompress_read_croptable(const char *filename, const int cropid, CropTable_t *C)
{
  unsigned int row;
  FILE *fh;
  unsigned short skip, take;
  int id;

  fh = fopen(filename,"r");
  if (fh == NULL)
  {
    printkerr("Failed to open crop table file %s.\n",filename);
    return ERROR_BAD_OR_MISSING_CROP_TABLE;
  }
  C->totalpix = 0;
  fscanf(fh,"%d",&id);
  if (id != cropid) {
      printkerr("Crop table with the wrong ID %d supplied.\n", id);
      return ERROR_BAD_OR_MISSING_CROP_TABLE;
  }
  fscanf(fh,"%u",&C->width);
  fscanf(fh,"%u",&C->height);  
  C->table = malloc(2*C->height*sizeof(unsigned short));
  for (row=0; row<C->height; row++)
  {
    if (fscanf(fh,"%hu",&skip) != 1)
    {
      printk("Crop table file %s is too short.\n",filename);
      return ERROR_BAD_OR_MISSING_CROP_TABLE;
    }
    C->table[2*row] = skip;
    if (fscanf(fh,"%hu",&take) != 1)
    {
      printk("Crop table file %s is too short.\n",filename);
      return ERROR_BAD_OR_MISSING_CROP_TABLE;
    }
    C->table[2*row+1] = take;
    C->totalpix += take;
  }
  fclose(fh);
  return 0;
}

int decompress_read_lutable(const char *filename, const int lutid, unsigned short *ILUT)
{
  int i, id;
  FILE *fh;

  fh = fopen(filename,"r");
  if (fh == NULL)
  {
    printkerr("Failed to open lookup table file %s.\n",filename);
    return ERROR_BAD_OR_MISSING_LOOKUP_TABLE;
  }
  fscanf(fh,"%d", &id);
  if (id != lutid) {
      printkerr("ILU table with the wrong ID %d supplied.\n", id);
      return ERROR_BAD_OR_MISSING_LOOKUP_TABLE;
  }
  for (i=0; i<(1<<14); i++)
  {
    fscanf(fh,"%hu",ILUT++);
  }
  fclose(fh);
  return 0;
}

static int decompress_context_initimage(int ctx, SciDataPacket_t *scipack)
{
  Decompress_Context_t *C;
  int cropid, lutid;
  unsigned int buflen;      
  int pad=0;
  int errcode;

  C = Context[ctx];

  /* Save first packet. */
  memcpy(&(C->image->firstpacket), scipack, sizeof(SciDataPacket_t));

  /* Get cropid. */
  cropid = scipack->cropid;
  if (cropid < 0 || cropid > MAX_CROPID)
    return ERROR_INVALIDCROPID;

  if (cropid == 0) {		// no crop
      croptables[cropid].totalpix = 4096*4096;
      switch (scipack->tapcode) {
	  case 2:
	  case 4:
	  case 5:
	  case 6:
	  case 7:
	  case 8:
	      croptables[cropid].width = 4096;
	      croptables[cropid].height = 4096;
	      break;
	  default:
	      croptables[cropid].width = 2048;
	      croptables[cropid].height = 8192;
      }
  } else if (croptables[cropid].width == 0) {	// crop table not yet read
      char filename[64];
      if (isHMI)
	  snprintf(filename, 64, "/home/production/cvs/EGSE/tables/crop/hmi/crop%d", cropid);
      else
	  snprintf(filename, 64, "/home/production/cvs/EGSE/tables/crop/aia/crop%d", cropid);
      errcode=decompress_read_croptable(filename, cropid, &croptables[cropid]);
      if (errcode) return errcode;
  }

  lutid = scipack->lutid;
  if (lutid && ilutables[lutid] == NULL) {	// ILU table not yet read
      char filename[64];
      ilutables[lutid] = malloc(16384*2);
      if (isHMI)
	  snprintf(filename, 64, "/home/production/cvs/EGSE/tables/lu/hmi/ilu%d", lutid);
      else
	  snprintf(filename, 64, "/home/production/cvs/EGSE/tables/lu/aia/ilu%d", lutid);
      errcode=decompress_read_lutable(filename, lutid, ilutables[lutid]);
      if (errcode) return errcode;
  }

  /* Get total pixel count from crop table. */
  C->image->stat.totalpix = croptables[cropid].totalpix;
  buflen = croptables[cropid].width*croptables[cropid].height
           * sizeof(unsigned short);  

  /* Fill pixelbuf with MISSING = 0xffff. */
  C->pixelbuf = malloc(buflen);
  memset(C->pixelbuf, 0xff, buflen);

  C->image->width = croptables[cropid].width;
  C->image->height = croptables[cropid].height;
  if (scipack->tapcode ==0 || scipack->tapcode ==1 || scipack->tapcode ==3) {
      C->image->width *= 2;
      C->image->height /= 2;
  }
  /* Fill pixels with MISSING = -32768. */
  if (buflen % 2880)
      pad = 2880 - (buflen % 2880);
  C->image->data = malloc(buflen + pad);
  for (int i=0;i<buflen/2;++i) C->image->data[i] = -32768;
  /* Pad with zeros. */
  memset(((char *) C->image->data)+buflen, 0, pad);
  return SUCCESS;
}

static int decompress_packet(const unsigned int N, const unsigned int K, 
			     const unsigned int SAT, unsigned short *pixels, 
			     int *new_count, unsigned short *indata)
{
  unsigned int *in;
  unsigned int b, nz, sign, nbits, c31mk, c32mnmk;
  unsigned kmask, nmask;
  unsigned int diff, error;	
  unsigned int pixcnt, nmk, wordcnt, newcnt=0;
  unsigned short oldval, val;
 
  assert(N>=NMIN && N<=NMAX);
  assert(K>=KMIN && K<=KMAX);
  assert(SAT<=8);
  assert(table_initialized);
  
  /* Precalculate various constants that depend on n and k. */
  kmask = (1<<K)-1;
  nmask = (1<<N)-1;
  nmk = N-K;
  c31mk = 31-K;
  c32mnmk = 32-nmk;

  /* Decode first pixel uncompressed. */ 
  in = (unsigned int *)indata;
  wordcnt = 0; /* */
  b = *in++; /* b contains words [wordcnt:wordcnt+1] */
  oldval = (unsigned short) b;
  if (*pixels == 0xffff)
    newcnt++;
  *pixels++ = oldval;
  b = b >> 16;
  nbits = 16; /* This many "live" bits are left in b. */ 
  pixcnt = 1; /* This many pixels have been decoded. */
  

  /* Rice decompression loop with first differencing. */
  while ((wordcnt <= PACKETDATAWORDS-3))
  {  
    /* Decode k low order bits. */
    diff = b & kmask;
    if (unlikely(nbits<K))
    {
      wordcnt += 2;
      if (wordcnt >= PACKETDATAWORDS-1)
	  break;
      b = *in++;
      diff += (b << nbits) & kmask;
      b = b >> (K-nbits);
      nbits += c31mk;

      /* Decode sign bit */
      sign = b & 1;    
      b = b>>1;

      /* Decode fundamental sequence. */
      /* We know: k < 14 => nbits = 31-k > 17 so we can
         go ahead and decode the fundamental sequence 
         with a table lookup. */
      nz = numzero16[(unsigned short) b];

      /* ADD ERROR CHECK HERE!!! */
      b = b >> (nz+1);
      nbits -= nz+1;          
    }
    else
    {
      if (likely(nbits>K+1))
      {
	nbits -= K+1;
	b = b >> K;
	/* Decode sign bit */
	sign = b & 1;    
	b = b>>1;

	/* Decode fundamental sequence */
	if (likely(nbits<17))
	{
	  nz = numzero16[b];
	  
	  if (likely(nz<nbits))
	  {
	    b = b>>(nz+1);
	    nbits -= nz+1;
	    /* ADD ERROR CHECK HERE!!! */
	  }	
	  else
	  {
	    unsigned int tmp;
	    wordcnt += 2;
	    if (wordcnt >= PACKETDATAWORDS-1)
		break;
	    b = *in++;
	    tmp = numzero16[(unsigned short) b];
	    //	    hist[nz] += 1;

	    b = b>>(tmp+1);
	    nz = nbits + tmp;
	    nbits = 31-tmp;
	    /* ADD ERROR CHECK HERE!!! */
	  }
	}
	else
	{
	  nz = numzero16[(unsigned short) b];
	  b = b>>(nz+1);
	  nbits -= nz+1;
	    /* ADD ERROR CHECK HERE!!! */
	}	
      }
      else 
      {
	if (nbits==K)
	{
	  wordcnt += 2;
	  if (wordcnt >= PACKETDATAWORDS-1)
	      break;
	  b = *in++;
	  /* Decode sign bit */
	  sign = b & 1;    
	  b = b>>1;
	  /* Decode fundamental sequence */
	  nz = numzero16[(unsigned short) b];
	  b = b>>(nz+1);
	  nbits = 30 - nz;
	}
	else 
	{
	  /* Decode sign bit */
	  sign = (b >> K) & 1;
	  wordcnt += 2;
	  if (wordcnt >= PACKETDATAWORDS-1)
	      break;
	  b = *in++;
	  /* Decode fundamental sequence */
	  nz = numzero16[(unsigned short) b];
	  b = b>>(nz+1);
	  nbits = 31 - nz;
	}
      }
    }

    /* Decode the remaining (n-k) bits of diff. */
    if (unlikely(nz >= 16))
      break;
    else if (unlikely(nz==SAT)) 
    {
      diff = (diff | (b<<K))  & nmask;      
      
      if (unlikely(nbits<nmk)) 
      {
	wordcnt += 2;
	if (wordcnt >= PACKETDATAWORDS-1)
	    break;
	b = *in++;
	diff += (b<<(K+nbits)) & nmask;
	b = b >> (nmk-nbits);
	nbits += c32mnmk;        
      }
      else 
      {
	b = b >> nmk;
	nbits -= nmk;        
      }
    }
    else
      diff += nz << K;

    /* Reconstruct pixel value. */
    if (*pixels == 0xffff)
      newcnt++;
    *pixels = (sign==1 ? oldval-diff : oldval+diff);  
    oldval = *pixels++;
    ++pixcnt;
  }   

  error = 0;
  if (unlikely(nz==16))
  {
    /* 
       Illegal FS detected. This can either mean that 
       we have encountered the end of a partial packet
       or that the packet was corrupted. The remaining part of 
       a partial packet must be filled with zeros, except for
       the final word, which must match the latest decoded value.
    */
    error = (diff != 0) || (sign!=0);
    ++wordcnt;
    while (wordcnt < PACKETDATAWORDS-1)
    {
      error |= (indata[wordcnt++]!= 0);
    }
  }

  if (pixcnt==1)
    oldval=0;

  /* Check that the last word matches the last decoded
     value. */
  val = indata[PACKETDATAWORDS-1];
  error |= (val!=oldval);

  
  if (unlikely(error))
  {
    printk("ERROR: val != oldval: %u != %u at pixcnt = %d\n",
	   val,oldval,pixcnt);
        pixcnt = -pixcnt; 
  }

  *new_count =  newcnt;
  return pixcnt;
}

int decompress_reconstruct(SciDataPacket_t *scipack, int ctx, Image_t **image)
{
  int cnt, new_cnt;
  unsigned int N, R, K, SAT;
  int errcode;
  Decompress_Stat_t *stat;
  Image_t *im;

  if (!table_initialized)   
    decompress_inittables();

  im = Context[ctx]->image;
  stat = &(im->stat);
  /* If this is the first image packet received for this context initialize
     various image info in the context structure. */
  if (im->firstpacket.cropid == -1)
  {
    if ((errcode = decompress_context_initimage(ctx, scipack)))
      return errcode;
  }
  

  /********
  Get the compression parameters from the header and decompress the
  pixels in the packet.
  *********/
  N = scipack->compid >> 3;
  R = scipack->bitselectid;
  K = scipack->compid & 0x7;
  SAT = 8;

  /*******
  Check if the current packet continues where we got to in 
  this image. If not, a packet was lost or packets got out of order.
  ********/
  if ( scipack->offset != stat->numpix )
  {
    if (scipack->offset > stat->totalpix)
    {
      printkerr("Warning: Offset count (%d) is larger than " \
		"the total number of pixels (%d) for image in context %d " \
		"with FSN=%lu, FID=%lu. Discarding image!\n", 
		scipack->offset,stat->totalpix,ctx,ID2FSN(stat->ID),
		ID2FID(stat->ID)); 
      errcode = ERROR_BADOFFSET;
      goto failpartial;
    }
    else
    {
      if ( scipack->offset < stat->numpix )
      {
	/* We are receiving pixels we have already seen. Maybe 
	   we are reprocessing a package that was retransmitted
	   from the DDS? */
	stat->backup_occured = 1;
      }
      else if ( scipack->offset > stat->numpix )
      {
	/* The offset field pixel count skips ahead of what we have 
	   decoded so far. Maybe a packet was lost or packets
	   got out of order? */
	stat->skip_occured = 1;
      }	
    }      
  }

  /* Special uncompressed mode */
  if (scipack->compid == 0 || (N==16&&R==0&&K==0)) {
      int npix,i;
      npix = stat->totalpix - scipack->offset;
      if (npix > PACKETDATAWORDS)
	  npix = PACKETDATAWORDS;
      for (i=npix; i<PACKETDATAWORDS; ++i)
	  if (scipack->data[i] != 0) {
	      errcode = ERROR_RAW_MODE_TRAILING_GARBAGE;
	      goto failpartial;
	  }
      memcpy(&Context[ctx]->pixelbuf[scipack->offset], scipack->data, 2*npix);
      stat->numpix += npix;
      if (stat->numpix > stat->totalpix) {
	errcode = ERROR_TOOMANYPIXELS;
	goto failpartial; 
      }
      ++(stat->npackets);
      return SUCCESS;
  }
      
  /* Now do the actual decompression of the packet payload. */
  cnt = decompress_packet(N, K, SAT,  &Context[ctx]->pixelbuf[scipack->offset],
			  &new_cnt, scipack->data);  
  ++(stat->npackets);

  if (cnt < 0)
  {
    errcode = ERROR_CORRUPTDATA;
    goto failpartial;
  }
  else if ( cnt!=new_cnt && (new_cnt>0))
  {
    errcode = ERROR_PARTIALOVERWRITE;
    goto failpartial; 
  } 
  else if (stat->numpix+new_cnt > stat->totalpix)
  {    
    /* AAARGH!!! Very bad boy! We probably stepped on some memory... */
    // Get status information for this image.
    errcode = ERROR_TOOMANYPIXELS;
    goto failpartial; 
  }
  else
    stat->numpix += new_cnt;

  //  decompress_print_status(stat);
  return SUCCESS;

 failpartial:
  /* 
     Decoding failed after decoding part of an image. Report status of
     the image when the error occured and return what was decoded so
     far.
  */
  *image = decompress_prepare_image(ctx);
  return errcode;
}

int check_completeness(Image_t **image)
{
  int status;
  unsigned int ctx;
  Image_t *im;
  Decompress_Stat_t *stat;

  /********
  Check if we have completed the image in the current context.
  Completion occures either 
    a) by decoding as many unique pixels as expected by the croptable
  or 
    b) by reaching the last pixel.  
  *********/

  *image = NULL;
  status = SUCCESS; /* Signal that we successfully appended the pixels
		     or added HK keywords an image which is still 
		     incomplete. */

  for (ctx=0; ctx<numcontexts; ctx++)
  {
    im = Context[ctx]->image; 
    stat = &im->stat;
    /**** Change definition of "image complete"
    if (im->firstpacket.cropid != -1 && stat->numpix == stat->totalpix 
       	   && im->keywords != NULL  )
    ****/
    if (im->firstpacket.cropid != -1 && im->keywords != NULL &&
	    Context[ctx]->pixelbuf[stat->totalpix-1] != 0xffff)
    {
      /* Woohoo! We got all the pixels for this image. Now undo the
	 table lookup, bit select and cropping and hand the final image
	 back to the caller. */
      *image = decompress_prepare_image(ctx);
      status = SUCCESS_IMAGECOMPLETE; /* Signal that an image was completed. */
    }
  }
  return status;
}

int decompress_next_vcdu(unsigned short vcdu[PACKETWORDS], 
			 Image_t **image, CCSDS_Packet_t **hk_packets)
{
  int i, FSN, FID; 
  int ctx;
  int status = ERROR_NODATA;
  unsigned char *p1,*p2;
  static unsigned short buffer[PACKETWORDS];
  SciDataPacket_t scipack;
  unsigned short *w, *hkstart;
  IM_PDU_Packet_t im_pdu;
  CCSDS_Packet_t ccsds, *p, *p_hk;
  HK_Keyword_t *kw, *kw2;


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
  *image = NULL;
  do
  {
    /* Decode the CCSDS header. */
    hkstart = vcdu + (w-buffer);
    w = decode_ccsds(w, &ccsds);
    
    /* Branch depending on the APID. */
    switch(ccsds.apid)
    {  
    case 0:
      /* We have reached the end of the VCDU */      
      break;

    case APID_HMI_SCIENCE_1:
    case APID_HMI_SCIENCE_2:
      isHMI = 1;
    case APID_AIA_SCIENCE_1:
    case APID_AIA_SCIENCE_2:
      /* This is a high-rate science data packet. */
      decode_scidata(w, &scipack);

      /*******
	Extract unique FSN and FID from the header and find the context 
	for that image.
      ********/
      FSN = UNPACK_FSN(&scipack);
      FID = UNPACK_FID(&scipack);
      
      ctx = decompress_findcontext(FSN, FID);
      if (ctx>=numcontexts)
      {
	/* Not previously seen ID => create a new context. */
	if ((ctx = decompress_new_context(FSN, FID)) == -1)
	  return ERROR_CTXOVERFLOW;
      }
      if ((status = decompress_reconstruct(&scipack, ctx, image)) < 0)
	return status;      
      break;
      
    case APID_HMI_TIME_1:
    case APID_HMI_TIME_2:
    case APID_AIA_TIME_1:
    case APID_AIA_TIME_2:
    case 460:	/******* TEMP ***********/
      /* This is an empty timestamp packet.  Do nothing. */
      status = SUCCESS;
      break;

    case 417:
    case 407:
    case 517:
    case 507:
      // DCHRI test packet
      // verify test pattern
      // return bogus SUCCESS_HK so hmi_lev0 no-ops it
      for (i=11; i<888; ++i) {
	  if (buffer[i] != 0xc0b+(i-11)) {
	      printk("DC/HRI test pattern error at IM_PDU ID %d Counter %lu\n",
		      im_pdu.im_pdu_id, im_pdu.im_pdu_count);
	      break;
	  }
      }
      status = SUCCESS_HK;
      break;

    case APID_HMI_IMSTAT_1:
    case APID_HMI_IMSTAT_2:
      /* This is an image status packet. */
      status = decode_hk_keywords(hkstart, ccsds.apid, &ccsds.keywords);
      if (status==SUCCESS)
      {      
	/* Extract FSN and FID keywords. */
	FSN = -1;
	FID = -1;
	kw = ccsds.keywords;      
	
	/* ASSUMPTION: These packets contains keywords of the form
	   XXX_SEQ_FILTERGRAM_ID, where XXX = "HMI" or "AIA" */
	while(kw)
	{
	  if (!strncasecmp("SEQ_FILTERGRAM_ID", kw->name + 4, 
			   strlen("SEQ_FILTERGRAM_ID")))
	    FID = kw->raw_value;
	  if (!strncasecmp("SEQ_FILTERGRAM_SN", kw->name + 4, 
			   strlen("SEQ_FILTERGRAM_SN")))
	    FSN = kw->raw_value;
	  kw = kw->next;
	}

	if (FSN == -1)
	  return ERROR_MISSING_FSN;
	if ( FID == -1)
	  return ERROR_MISSING_FID;
	
	/* find image context and attach keywords to it. */
	ctx = decompress_findcontext(FSN,FID);  
	if (ctx>=numcontexts)
	{
	  /*   Not previously seen ID => create a new context. */
	  if ((ctx = decompress_new_context(FSN,FID)) == -1)
	    return ERROR_CTXOVERFLOW;
	}

	/* Copy keywords. */
	Context[ctx]->image->keywords = copy_hk_keywords(ccsds.keywords);
	
	/* Allocate a CCSDS packet. */
        p = malloc(sizeof(CCSDS_Packet_t));
        assert(p);
	memcpy(p, &ccsds, sizeof(CCSDS_Packet_t));
	/* Append to output list. */
	if (*hk_packets == NULL)
	  *hk_packets = p_hk = p;
	else
	{
	  p_hk->next = p;
	  p_hk = p_hk->next;
	}
	p_hk->next = NULL;
	status = SUCCESS_HK;
	goto DONE_WITH_THIS_VCDU; /**********TEMP***********/
      }
      /*      else
	      status =  ERROR_DECODEHKFAIL;*/
      break;      

    case APID_AIA_IMSTAT_1:
    case APID_AIA_IMSTAT_2:
      /* This is an image status packet. */
      status = decode_hk_keywords(hkstart, ccsds.apid, &ccsds.keywords);
      if (status==SUCCESS)
      {      
	/* Extract FSN and FID keywords. */
	FSN = -1;
	FID = -1;
	kw = ccsds.keywords;      
	
	while(kw)
	{
	  if (!strncasecmp("SEQ_HEADER", kw->name + 4, 
			   strlen("SEQ_HEADER")))
	    FSN = kw->raw_value;
	  kw = kw->next;
	}

	if (FSN == -1)
	  return ERROR_MISSING_FSN;
	
	/* find image context and attach keywords to it. */
	ctx = decompress_findcontext(FSN,FID);  
	if (ctx>=numcontexts)
	{
	  /*   Not previously seen ID => create a new context. */
	  if ((ctx = decompress_new_context(FSN,FID)) == -1)
	    return ERROR_CTXOVERFLOW;
	}

	/* Copy keywords. */
	Context[ctx]->image->keywords = copy_hk_keywords(ccsds.keywords);
	
	/* Allocate a CCSDS packet. */
        p = malloc(sizeof(CCSDS_Packet_t));
        assert(p);
	memcpy(p, &ccsds, sizeof(CCSDS_Packet_t));
	/* Append to output list. */
	if (*hk_packets == NULL)
	  *hk_packets = p_hk = p;
	else
	{
	  p_hk->next = p;
	  p_hk = p_hk->next;
	}
	p_hk->next = NULL;
	status = SUCCESS_HK;
	goto DONE_WITH_THIS_VCDU; /**********TEMP***********/
      }
      /*      else
	      status =  ERROR_DECODEHKFAIL;*/
      break;      
    default:
     /* This should be a housekeeping packet. */
      status = decode_hk_keywords(hkstart, ccsds.apid, &(ccsds.keywords));
      if (status == SUCCESS)
      {
	/* Allocate a CCSDS packet. */
        p = malloc(sizeof(CCSDS_Packet_t));
        assert(p);
	memcpy(p, &ccsds, sizeof(CCSDS_Packet_t));
	/* Append to output list. */
	if (*hk_packets == NULL)
	  *hk_packets = p_hk = p;
	else
	{
	  p_hk->next = p;
	  p_hk = p_hk->next;
	}
	p_hk->next = NULL;
	status = SUCCESS_HK;
      }
      break;
    }
    w += ccsds.length/2;
  } while ( /* status >= 0 && */ /********TEMP*******/ ccsds.apid>0 && (int)(w-buffer) < PACKETWORDS);

DONE_WITH_THIS_VCDU: /**********TEMP***********/

  /* Check if any pictures are complete as a result of
     what was decoded from this VCDU. */
  if (status >= 0)
    if (SUCCESS_IMAGECOMPLETE == check_completeness(image))
	return SUCCESS_IMAGECOMPLETE;

  return status;
}

// rice_encode2.c

#ifdef KMAX
#undef KMAX
#endif

#define KMAX 13		// max number of unencoded low-order bits
#define KBITS 4		// number of bits required to represent K

#define CHECK_OVERRUN \
    if (c > out+bufsz) { free(d); return RICE_ENCODE_OVERRUN; }

#define FLUSH \
    *c++=buf>>24; *c++=buf>>16; *c++=buf>>8; *c++=buf; buf = 0; CHECK_OVERRUN

#define PUT_N_BITS(n, val) \
    if (bits2go > n) { \
	bits2go -= n; \
	buf += val << bits2go; \
    } else if (bits2go == n) { \
	buf += val; \
	FLUSH \
	bits2go = 32; \
    } else { \
	buf += val >> (n - bits2go); \
	FLUSH \
	bits2go = 32 - (n - bits2go); \
	buf = val << bits2go; \
    }

#define PUT_FS(n) \
    if (bits2go > n+1) { \
	bits2go -= n+1; \
	buf += 1 << bits2go; \
    } else if (bits2go == n+1) { \
	++buf; \
	FLUSH \
	bits2go = 32; \
    } else { \
	FLUSH \
	bits2go += 32 - (n+1); \
	while (bits2go < 0) { \
	    c += 4; \
	    bits2go += 32; \
	} \
	buf = 1 << bits2go; \
    }

int rice_encode2(
    const short *in,		// input array of pixels
    int nin,			// number of input pixels
    unsigned char *out,		// encoded output byte array
    int bufsz,			// output buffer size
    int blksz)			// number of pixels in a coding block
{
    unsigned short *d;		// array of 1st differences to be encoded
    unsigned char *c;		// pointer to current byte in output
    unsigned buf = 0;		// bit buffer
    unsigned val;		// value to insert into bit buffer
    int bits2go = 32;		// vacancy in bit buffer
    int k;			// number of low order bits to split
    int i, j;
    short curr, last, delta;
    unsigned short top, bottom, kmask, tmp;
    double pixsum;

    d = (unsigned short *) malloc(blksz * sizeof(unsigned short));
    if (!d) return RICE_ENCODE_OUT_OF_MEMORY;

    memset(out, 0, bufsz);

    // output first pixel raw (most significant byte first)
    last = in[0];
    val = last;
    out[0] = val >> 8;
    out[1] = val;

    c = out + 2;

    // Main loop.  Note that index starts at 0.  (d[0] is always 0.)
    for (i = 0; i < nin; i += blksz) {
	pixsum = 0.0;
	if (nin - i < blksz) blksz = nin - i;
	for (j = 0; j < blksz; ++j) {
	    curr = in[i+j];
	    if (last == curr) {
		d[j] = 0;
		continue;
	    }
	    delta = curr - last;	// overflow OK
	    last = curr;
	    d[j] = (delta < 0) ? ~(delta << 1) : (delta << 1);
	    pixsum += d[j];
	}

	// Zero entropy case: output KBITS zero bits
	if (pixsum == 0.0) {
	    if (bits2go > KBITS)
		bits2go -= KBITS;
	    else {
		FLUSH
		bits2go += 32 - KBITS;
	    }
	    continue;
	}

	// Find k. 
	tmp = .5*pixsum/blksz;
	for (k=0; tmp; ++k) tmp >>= 1;

	// Small entropy case, k == 0
	if (k == 0) {
	    val = 1;
	    PUT_N_BITS(KBITS, val)
	    for (j = 0; j < blksz; ++j) {
		val = d[j];
		PUT_FS(val)
	    }
	// Medium entropy case, k <= KMAX
	} else if (k <= KMAX) {
	    val = k+1;
	    PUT_N_BITS(KBITS, val)
	    kmask = (1u<<k) - 1u;
	    for (j = 0; j < blksz; ++j) {
		val = d[j];
		top = val >> k;
		bottom = val & kmask;
		PUT_FS(top)
		PUT_N_BITS(k, bottom)
	    }
	// High entropy case, k > KMAX
	} else {
	    val = KMAX+2;
	    PUT_N_BITS(KBITS, val)
	    for (j = 0; j < blksz; ++j) {
		val = d[j];
		PUT_N_BITS(16, val);
	    }
	}
    }

    free(d);

    // final flush
    i = 4 - bits2go/8;
    if (i-- > 0) *c++ = buf >> 24;
    if (i-- > 0) *c++ = buf >> 16;
    if (i-- > 0) *c++ = buf >> 8;
    if (i-- > 0) *c++ = buf;

    return (c-out <= bufsz) ? c-out : RICE_ENCODE_OVERRUN;
}

#ifdef KMAX
#undef KMAX
#endif
