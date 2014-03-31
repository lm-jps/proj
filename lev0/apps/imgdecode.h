#ifndef IMGDECODE_INCL
#define IMGDECODE_INCL 1

#ifndef PACKETHEADERWORDS
#define PACKETHEADERWORDS               19
#endif
#ifndef PACKETDATAWORDS
#define PACKETDATAWORDS                 869
#endif
#ifndef PACKETWORDS
#define PACKETWORDS                     888
#endif

#define MAXPIXELS			16777216	// 4096*4096
#define MAXHIST				32768

#ifndef TABLE_DIR
#define TABLE_DIR 			"/home/production/cvs/TBL_JSOC/lev0"
#endif

#define APID_HMI_SCIENCE_1 		400
#define APID_HMI_SCIENCE_2 		410
#define APID_AIA_SCIENCE_1 		500
#define APID_AIA_SCIENCE_2 		510

#define isAIA(x)	((x) == APID_AIA_SCIENCE_1 || (x) == APID_AIA_SCIENCE_2)
#define isHMI(x)	((x) == APID_HMI_SCIENCE_1 || (x) == APID_HMI_SCIENCE_2)

#define IMGDECODE_DECOMPRESS_ERROR	(-1)
#define IMGDECODE_TOO_MANY_PIXELS	(-2)

#define IMGDECODE_BAD_N			(-101)
#define IMGDECODE_BAD_APID		(-102)
#define IMGDECODE_NO_LOOKUP_TABLE	(-103)
#define IMGDECODE_LOOKUP_ID_MISMATCH	(-104)
#define IMGDECODE_BAD_LOOKUP_TABLE	(-105)
#define IMGDECODE_NO_CROP_TABLE		(-106)
#define IMGDECODE_CROP_ID_MISMATCH	(-107)
#define IMGDECODE_BAD_CROP_GEOMETRY	(-108)
#define IMGDECODE_BAD_CROP_TABLE	(-109)
#define IMGDECODE_BAD_CROP_SKIP_TAKE	(-110)
#define IMGDECODE_BAD_OFFSET		(-111)

#define IMGDECODE_OUT_OF_MEMORY		(-201)

#define	BLANK				(-32768)

typedef struct {
    unsigned short width, height;
    unsigned totalpix;
    unsigned short *skip;
    unsigned short *take;
    unsigned *offset;
} CROPTABLE;

typedef struct {
    int initialized;
    int reopened;
    unsigned telnum;
    unsigned fsn;
    unsigned fid;
    unsigned apid;
    unsigned cropid;
    unsigned luid;
    unsigned tap;
    unsigned N, K, R;
    unsigned totalvals;
    unsigned datavals;
    unsigned npackets;
    unsigned nerrors;
    int last_pix_err;
    int overflow;
    int headerr;
    uint64_t first_packet_time;
    short dat[MAXPIXELS];
    unsigned hist[MAXHIST];
} IMG;

typedef struct {
    short min;
    short max;
    short median;
    double mean;
    double rms;
    double skew;
    double kurt;
} STAT;

int imgdecode(unsigned short *impdu, IMG *img);
int imgdecode_init_hack(IMG *img);
int imgstat(IMG *img, STAT *stat);

// Image corruption patches
// Patch 1: crop table corruption Dec 2011 - Jan 2012
#define NEED_PATCH1(fsn) (!(fsn%2) && (fsn<=33190636) && (fsn>=32249884))
// Patch 2: camera 1 lookup table corruption 30 March 2014 
#define NEED_PATCH2(fsn) (!(fsn%2) && (fsn<=70405129) && (fsn>=70397305))

#endif
