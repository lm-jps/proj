#ifndef _EGSEHMICOMP_H
#define _EGSEHMICOMP_H

/* Maximum number of images that can be decoded simultaneously. */
#define MAXCONTEXTS (256)

/* Limits for compression parameters. */
#define NMAX      (16)
#define NMIN       (6)
#define RMAX       (8)
#define RMIN       (0)
#define KMAX       (7)
#define KMIN       (0)

#define UNPACK_FSN(H) (((unsigned int)(H)->ccdhead[0]<<16)+((unsigned int)(H)->ccdhead[1]))
#define UNPACK_FID(H) (((unsigned int)(H)->ccdhead[2]<<16)+((unsigned int)(H)->ccdhead[3]))
#define IMAGE_FSN(I) UNPACK_FSN(&(I->firstpacket))
#define IMAGE_FID(I) UNPACK_FID(&(I->firstpacket))
#define ID2FSN(ID) ((uint32_t) (ID) & 0xffffffff)
#define ID2FID(ID) ((uint32_t) ((ID)>>32) & 0xffffffff)
#define UNIQUEID(FSN,FID)  ((int64_t)(FSN))
#define unlikely(a) __builtin_expect((a), 0)
#define likely(a) __builtin_expect((a), 1)

/* Maximum number of different crop tables. */
#define MAX_CROPID ((1<<12)-1)

/* Maximum number of different lookup tables. */
#define MAX_LUTID  ((1<<8)-1)

#define SUCCESS_HKCOMPLETE         (3)
#define SUCCESS_HK                 (2)
#define SUCCESS_IMAGECOMPLETE      (1)
#define SUCCESS                    (0) // collides with a system define
#define ERROR_BADOFFSET           (-1)
#define ERROR_CORRUPTDATA         (-2)
#define ERROR_BADHEADER           (-3)
#define ERROR_CTXOVERFLOW         (-4)
#define ERROR_INVALIDID           (-5)
#define ERROR_TOOMANYPIXELS       (-6)
#define ERROR_NOSUCHIMAGE         (-7)
#define ERROR_PARTIALOVERWRITE    (-8)
#define ERROR_WRONGPACKET         (-9)
#define ERROR_MISSING_FSN        (-10)
#define ERROR_MISSING_FID        (-11)
#define ERROR_INVALIDCROPID      (-12)
#define ERROR_NODATA             (-13)
#define ERROR_HK_UNKNOWN_APID              (-14)
#define ERROR_HK_CANNOT_FIND_VER_NUM       (-15)
#define ERROR_HK_CANNOT_LOAD_HK_VALUES     (-16)
#define ERROR_HK_CANNOT_LOAD_ENGR_VALUES   (-17)
#define ERROR_HK_INVALID_BITFIELD_LENGTH   (-18)
#define ERROR_HK_UNHANDLED_TYPE            (-19)
#define ERROR_HK_NOSUCHDIR                 (-20)
#define ERROR_HK_CANNOT_LOAD_CONFIG        (-21)
#define ERROR_RAW_MODE_TRAILING_GARBAGE    (-22)
#define ERROR_BAD_OR_MISSING_CROP_TABLE    (-23)
#define ERROR_BAD_OR_MISSING_LOOKUP_TABLE  (-24)
#define ERROR_EMPTY_IMAGE_CONTEXT          (-25)

#define PACKETHEADERWORDS   (19) 
#define PACKETHEADERBYTES (2*19) 
#define PACKETHEADERBITS (16*19) 

#define PACKETDATAWORDS   (869) 
#define PACKETDATABYTES (2*869) 
#define PACKETDATABITS (16*869) 

#define PACKETWORDS  (PACKETHEADERWORDS+PACKETDATAWORDS) 
#define PACKETBYTES  (PACKETHEADERBYTES+PACKETDATABYTES) 
#define PACKETBITS   (PACKETHEADERBITS+PACKETDATABITS)

#define APID_HMI_SCIENCE_1		400
#define APID_HMI_SCIENCE_2		410
#define APID_HMI_TIME_1			405
#define APID_HMI_TIME_2			415
#define APID_HMI_IMSTAT_1		445
#define APID_HMI_IMSTAT_2		475

#define APID_AIA_SCIENCE_1		500
#define APID_AIA_SCIENCE_2		510
#define APID_AIA_TIME_1			505
#define APID_AIA_TIME_2			515
#define APID_AIA_IMSTAT_1		529
#define APID_AIA_IMSTAT_2		569

#define KEYWORD_NAME_SIZE                          10
#define KEYWORD_TYPE_SIZE                          10
#define TELEMETRY_MNEMONIC_SIZE                    50
#define HK_DECODER_SUCCESSFUL                      SUCCESS
#define HK_DECODER_ERROR_UNKNOWN_APID              ERROR_HK_UNKNOWN_APID
#define HK_DECODER_ERROR_CANNOT_FIND_VER_NUM       ERROR_HK_CANNOT_FIND_VER_NUM
#define HK_DECODER_ERROR_CANNOT_LOAD_HK_VALUES     ERROR_HK_CANNOT_LOAD_HK_VALUES
#define HK_DECODER_ERROR_CANNOT_LOAD_ENGR_VALUES   ERROR_HK_CANNOT_LOAD_ENGR_VALUES
#define HK_DECODER_ERROR_UNKNOWN_ENV_VARIABLE      ERROR_HK_CANNOT_LOAD_ENGR_VALUES
#define MAX_CHAR_VERSION_NUMBER 10
#define MAX_PACKET_ID_TYPE      100
#define MAX_DATE_SIZE           15
#define MAX_APID_POINTERS       4095
#define MAX_FILE_NAME           300
#define MAX_DIRECTORY_NAME      200
#define MAX_SIZE_CHANGE_TIME    15
#define MAX_SIZE_CHANGE_DATE    15
#define MAX_HK_MNM              30
#define MAX_HK_KYWD             10
#define MAX_HK_VALUE_TYPE       10
#define MAXLINE_IN_FILE         150
#define MAX_VALUE_DSC           100
#define MAX_NUMBER_COFFICIENTS  6
#define SUCCESSFUL              0
#define HMI_ID_TYPE             "HMI"
#define AIA_ID_TYPE             "AIA"
#define OTHER_ID_TYPE           "SSIM"
#define MAX_NUM_KW_LINES        3000
#define MAX_NUM_ACON_LINES      1000
#define MAXLINE_ACON_IN_FILE    200
#define MAX_NUM_DCON_LINES      4000
#define MAXLINE_DCON_IN_FILE    200
#define HK_MAX_TLM_NAME         200
#define ERROR_LOADING_HK_DATA   1


/* Crop table in skip-take format. */
typedef struct CropTable_struct
{
  unsigned int totalpix;
  unsigned int width, height;
  unsigned short *table;      /* Skip-take table with 2*totalpix entries. */
} CropTable_t;

/*************** Global table definitions *******************/
extern CropTable_t croptables[MAX_CROPID+1];
extern unsigned short *lutables[MAX_LUTID+1];
extern unsigned short *ilutables[MAX_LUTID+1];

/*
  Rice compression
*/

#define RICE_ENCODE_OUT_OF_MEMORY       -1
#define RICE_ENCODE_OVERRUN             -2
#define RICE_DECODE_OVERRUN             -3
#define RICE_DECODE_TRUNCATED_INPUT     -4
#define RICE_DECODE_TRAILING_GARBAGE    -5

extern int rice_encode1(const char *in, int nin, unsigned char *out, int bufsz, int blksz);
extern int rice_encode2(const short *in, int nin, unsigned char *out, int bufsz, int blksz);
extern int rice_encode4(const int *in, int nin, unsigned char *out, int bufsz, int blksz);
extern int rice_decode1(const unsigned char *in, int nin, char *out, int nout, int blksz);
extern int rice_decode2(const unsigned char *in, int nin, short *out, int nout, int blksz);
extern int rice_decode4(const unsigned char *in, int nin, int *out, int nout, int blksz);

// WARNING - these definitions of DRMS_Type_t and DRMS_Type_Value_t collide
// with the definitions in JSOC/base/drms/libs/aps/drms_types.h.
// It looks like lev0 uses definitions in JSOC/proj/lev0/apps/packets.h and 
// JSOC/base/drms/libs/aps/drms_types.h. The definitions were renamed in
// packets.h so that they don't collide with the ones in drms_types.h, but
// they do collide with the ones below. So, I moved the ones below
// from egsehmicomp.h to egsehmicomp.c.
/* Keyword value types. */
#ifndef LEV0SLOP
typedef enum {DRMS_TYPE_INT8, DRMS_TYPE_INT16, DRMS_TYPE_INT32,
              DRMS_TYPE_INT64, DRMS_TYPE_UINT8, DRMS_TYPE_UINT16,
              DRMS_TYPE_UINT32, DRMS_TYPE_UINT64,
              DRMS_TYPE_FLOAT, DRMS_TYPE_DOUBLE,
              DRMS_TYPE_TIME, DRMS_TYPE_STRING} DRMS_Type_t;

typedef union DRMS_Type_Value
{
  int8_t   int8_val;
  int16_t  int16_val;
  int32_t  int32_val;
  int64_t  int64_val;
  uint8_t  uint8_val;
  uint16_t uint16_val;
  uint32_t uint32_val;
  uint64_t uint64_val;
  float   float_val;
  double  double_val;
  double  time_val;
  char   *string_val;
} DRMS_Type_Value_t;
#endif

/* Decompression status for an image. */
typedef struct Decompress_Stat_struct
{
  /* Unique Image ID. */
  int64_t ID;
  time_t starttime;        /* Timestamp indicating when processing of the
                              first packet for this image began. */
  unsigned short npackets; /* Number of packets received. */
  unsigned int numpix;     /* Number of unique pixels decompressed. */
  unsigned int totalpix;   /* Total number of unique pixels expected for
                              image. */
  unsigned char backup_occured; /*  At some point the offset counter decreased,
                                    possibly due to reprocessing of part of an
                                    image.*/
  unsigned char skip_occured; /*  At some point the offset counter skipped
                                  forward, indicating lost or out of order
                                  packets.*/
} Decompress_Stat_t;

/* Telemetry packet consisting of (possibly byteswapped) raw data and
   parsed header. */
typedef struct SciDataPacket_struct
{
  /* Science packet headers. */
  unsigned int   shs;
  unsigned int   shss;
  unsigned short ccdhead[4];
  short cropid;
  unsigned char  romode;
  unsigned char  headererr;
  unsigned char  oflow;
// TAP codes
//  H--G     0 - EFGH
//  |  |     1 - EF   2 - FG   3 - GH   4 - HE
//  E--F     5 - E    6 - F    7 - G    8 - H
  unsigned char  tapcode;
  unsigned char  bitselectid;
  unsigned char  compid;
  unsigned char  lutid;
  unsigned int   offset;

  unsigned short *data;
} SciDataPacket_t;

#ifndef LEV0SLOP
/************* Keyword struct ****************/
#define MAX_KEYWORD_NAME_SIZE 64
typedef struct HK_Keyword_struct {
     char           name[MAX_KEYWORD_NAME_SIZE];
     char           fitsname[9];
     int64_t        raw_value;    /* 64 bit integer should be able to
                                     hold any header field */
     DRMS_Type_t       eng_type;  /* Engineering value type. */
     DRMS_Type_Value_t eng_value; /* Engineering value. */
     struct HK_Keyword_struct  *next;
} HK_Keyword_t;

typedef struct Image_struct
{
  /* Unique Image ID. */
  int64_t ID;

  /* Image data. */
  short width, height;
  short *data;

  /* Decompressor statistics. */
  struct Decompress_Stat_struct stat;

  /* Meta data. */
  SciDataPacket_t firstpacket;
  HK_Keyword_t *keywords;

  /* Linked list pointer in case we need to return
     multiple completed images at a time. */
  struct Image_struct *next;
} Image_t;

/******************* CCSDS Packet struct ********************/
typedef struct CCSDS_Packet_struct 
{
  /* CCSDS Header fields. */
  unsigned char  version;
  unsigned char  type;
  unsigned char  shf;
  unsigned short apid;
  unsigned char  sf;
  unsigned short ssc;
  unsigned short length;

  /* Decoded Keywords. */
  HK_Keyword_t *keywords;

  /* Next pointer */
  struct CCSDS_Packet_struct *next;
} CCSDS_Packet_t;
#endif

/* Decompression context for an image. */
typedef struct Decompress_Context_struct
{
  /* Unique Image ID. */
  int64_t ID;
  unsigned short *pixelbuf; /* Buffer to hold decoded pixels before
                               uncropping. */
  struct Image_struct *image;           /* Buffer for final image. */
} Decompress_Context_t;

int (*printk)(const char *fmt, ...);
int (*printkerr)(const char *fmt, ...);
void printk_set(int (*std)(const char *fmt, ...),
                int (*err)(const char *fmt, ...));

#ifndef LEV0SLOP
int decompress_status_all(Decompress_Stat_t **stat);
void decompress_print_status( Decompress_Stat_t *stat);


/* Undo pixel value transformation based on lookup table. */
void decompress_undotransform(const unsigned int N, const unsigned int R, 
			      const unsigned int ILUTID, 
			      const unsigned int numpix, unsigned short *pixels);

/* Insert pixel values from a 1-D pixel buffer into an image according 
   the a crop table. */
void decompress_uncrop(const unsigned int CROPID, const unsigned int numpix, 
		       unsigned short *pixels, Image_t *image);
Image_t *decompress_free_context(unsigned int ctx, int discard_image);
int decompress_flush_image(unsigned int FSN, unsigned int FID, 
			   Image_t **image);
int decompress_writefitsimage(const char* file, Image_t *image, int compress);
void decompress_free_hk(CCSDS_Packet_t *p);
void decompress_free_images(Image_t *image);
void decompress_inittables(void);
int decompress_read_croptable(const char *filename, const int cropid, CropTable_t *C);
int decompress_read_lutable(const char *filename, const int lutid, unsigned short *ILUT);
int check_completeness(Image_t **image);
int decompress_next_vcdu(unsigned short vcdu[PACKETWORDS], 
			 Image_t **image, CCSDS_Packet_t **hk_packets);
#endif


#endif
