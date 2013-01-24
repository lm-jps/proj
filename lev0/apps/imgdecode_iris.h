#ifndef IMGDECODE_IRIS_INCL
#define IMGDECODE_IRIS_INCL 1

#ifndef PACKETHEADERWORDS
#define PACKETHEADERWORDS               19
#endif
#ifndef PACKETDATAWORDS
#define PACKETDATAWORDS                 869
#endif
#ifndef PACKETWORDS
#define PACKETWORDS                     888
#endif

#if 0
#ifndef IMGX
#define IMGX 1096
#endif

#ifndef IMGY
#define IMGY 4144
#endif

#define MAXPIXELS			(IMGX*IMGY)
#endif

#define MAXPIXELS			(1096*4144)

#define MAXHIST				32768

#ifndef TABLE_DIR
#define TABLE_DIR 			"/home/prodtest/cvs/TBL_JSOC/lev0"
#endif

#define APID_IRIS_SCIENCE		768

#define isIRIS(x)	((x) == APID_IRIS_SCIENCE)

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
    unsigned isysn;
    unsigned apid;
    unsigned cropid;
    unsigned luid;
    unsigned tap;
    unsigned N, K, R;
    unsigned totalvals;
    unsigned datavals;
    unsigned npackets;
    unsigned nerrors;
    unsigned nx;
    unsigned ny;
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

int imgdecode_iris(unsigned short *impdu, IMG *img);
int imgdecode_iris_init_hack(IMG *img);
int imgstat_iris(IMG *img, STAT *stat);

#endif
