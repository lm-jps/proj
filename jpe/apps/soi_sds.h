/*
 *  soi_sds.h							~CM/include
 *
 *  SDS is a generalized internal data structure supporting I/O among
 *    a number of standard external data representations.  It provides
 *    a rich structure modeled on the Common Data Format (CDF).  sds can
 *    be used as a common platform for extracting single variable data 
 *    from and saving data to a variety of external representations, 
 *    including CDF, FITS, TIFF (?), and WSO Datasets.
 *  C++ programs should include the corresponding file "soi_sds.hxx".
 *  Additional information is in the following man pages:
 *
 *  Responsible:
 *	Rick Bogart				    RBogart@solar.Stanford.EDU
 *
 *  Bugs:
 *    Only the FITS external format is supported, and not very completely
 *    The macros sIsImvalid and usIsImvalid (note spelling) are probably not
 *	used.  The GNU c compiler doesn't like them.  None of the validity
 *	macros are really correct.
 *
 *  Revision history is at the end of the file.
 */

#ifndef SOI_SDS_INCL 

/****************************************************************************/
/**************************  INCLUDE STATEMENTS  ****************************/
/****************************************************************************/

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "soi_error.h"
#include "soi_key.h"
#include "soi_machine.h"
#include "soi_missing.h"
#include "soi_NaN.h"
#include "soi_args.h"


/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define    SDS_QUIT(errcode, i) {soi_errno=errcode; return(i);}

#define MAX_TYPE       3
#define MAX_DATATYPE   5
#define MAX_BYTES     20
#define SDS_MAXSTRING	(255)
#define SDS_MAXPATH	(255)

#define FITSID        51
#define CDFID         52
#define MAXRANK      100
						       /*  fits file format  */
#define FITS_KWSIZE    8
#define FITS_VALSIZE  70
#define FITS_CARDSIZE 80
#define FITS_NCARDS   36
						      /*  datatypes defines  */
#define	SDS_VOID	(0)
#define	SDS_BYTE	(1)
#define SDS_UBYTE	(2)
#define SDS_SHORT	(3)
#define SDS_USHORT	(4)
#define SDS_INT		(5)
#define SDS_UINT	(6)
#define SDS_LONG	(7)
#define SDS_ULONG	(8)
#define SDS_FLOAT	(9)
#define SDS_DOUBLE	(10)
#define SDS_COMPLEX	(11)
#define SDS_STRING	(12)
#define SDS_TIME	(13)
#define SDS_LOGICAL	(14)
#define SDS_ANY		(15)
#define SDS_ASIS	(16)
						   /*  supported file types  */
#define SDS_NOTYPE	(0)
#define SDS_RAW         (1)
#define SDS_SDS         (2)
#define SDS_MDI_HRTLM   (3)
#define SDS_MDI_LRTLM   (4)
#define SDS_TIFF        (5)
#define SDS_FITS        (6)
#define SDS_CDF         (7)
#define SDS_FITS_TABLE	(8)
#define SDS_RDB		(9)
#define SDS_GIF		(10)
					     /*  SDS data pointer ownership  */
#define SDS_NEVER_FREE	(0)
#define SDS_OK_TO_FREE	(1)

/*  ALL_DATA_MISSING is defined in soi_error.h  */

#define SOI_SDS_VERSION_NUM	("4.5")
#define SOI_SDS_INCL

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/

typedef struct attribute {
  struct attribute   *next;      /*  Next attribute struct            */
  char               *attrname;  /*  Name of attribute                */
  void               *attrvalue; /*  Pointer to attribute value       */
  char               *comment;   /*  Holds a comment                  */
  int                datatype;   /*  Attribute value's type           */
} ATTRIBUTES;


typedef struct history {
  struct history    *next;	       /*  Pointer to next line of history  */
  char		    *text;	       /*  Contains one line of history     */
} HISTORY;

typedef struct record {	       /*  Stores info on what records in the CDF   */
			       /*  the data values are from                 */ 
  long            recvariance;	    /*  Record variance                     */
  long            startrecord;
  long            numrecords;
  long            recinterval;		/*  Interval between records.       */
} RECORD; 

typedef struct dimension {            /*  Stores what indices of what       */
	                              /*  dimensions the data are from      */
  int            rank;       /*  dimension     */
  int            *length;    /*  Num indices from each */
  int            *dimvary;   /*  array of Dimensional variance       */
  int            *startind;  /*  starting indices of each dimension  */
  int            *indinter;  /*  intervals between indices           */
} DIMENSION;

typedef struct sds_stats_struct {
  union	{
    double dm;
    float fm;
    long int lm;	unsigned int ulm;
    int im;		unsigned int uim;
    short sm;	unsigned short usm;
    char cm;	unsigned char ucm;
  } min;                   /*  Minimum of the data in native type  */
  union	{
    double dm;
    float fm;
    long int lm;
    int im;
    short sm;
    char cm;
  } min_abs;                /*  Minimum of absolute value           */
  union	{
    double dm;
    float fm;
    long int lm;	unsigned int ulm;
    int im;		unsigned int uim;
    short sm;	unsigned short usm;
    char cm;	unsigned char ucm;
  } max;                    /*  Maximum of the data in native type  */
  double        minact;     /*  Actual min of the data.             */
  double        maxact;     /*  Actual max of the data.             */
  double        rangemin;   /*  Legal min of the data.              */
  double        rangemax;   /*  Legal max of the data.              */
  double	median;	    /*  median value of the data.           */
  double        mean;       
  double        stdev;
  double        skewness;
  double	kurtosis;
  double	minabs;
  char          *var_name;    /*  Variable's name 		      */
  char          *format;      /*  C's printf format 		      */
  unsigned long   numvals;    /*  Number of elements of valid data    */
  unsigned long   missvals;   /*  Number of elements of missing data  */
  unsigned short  statvalid;  /*  Are the statistics valid?           */
} SDS_STATS;

typedef struct axis {
	 int		*axis_datatype;   /* datatype of independent vars    */
	 void           *start;     /* starting values for independent vars  */
	 void           *stepsize;  /* stepsizes of independent vars         */
	 char           **units;    /* strings with unit information         */
         char           **format;   /* strings with significant format info  */
	 char		**title;    /* the names of the independent vars     */
} AXIS;


typedef struct color {
  int		BitsPerPixel;
  int		*Red;
  int		*Green;
  int		*Blue;
} SDS_COLOR;


typedef struct sds {
  struct sds    *next;      /*  Pointer to next sds struct          */
  RECORD        *recinfo;   /*  Info on record numbers of the data  */
  DIMENSION     *diminfo;   /*  Info on dimension layout of data    */
  SDS_STATS     *stats;     /*  Statistics of sds's data            */
  SDS_COLOR	*color_table;
  char          *varname;   /*  Name of sds                         */
  void          *fillvalue; /*  Empty spot in data symbol           */
  void          *data;      /*  Pointer to data of the sds          */
  int		data_quality; /*  are data good? - sliding  scale   */
  int           datatype;   /*  Type of data                        */
  int           numbytes;   /*  Bytes per data element              */
  AXIS		*axisinfo;  /*  holds info of independent vars      */ 
  ATTRIBUTES    *attrinfo;  /*  holds local attribute info          */
  FILE		*fp;	    /*  used only when file is opened       */
  char          *filename;  /*  filename to the fits file           */
  int		own_data;   /*  set to disable "free" of data pointer  */
  double	bscale;	    /*  advisory scaling for binary output  */
  double	bzero;	    /*  advisory offset for binary output   */
  int		scale_on_write;	/*  advice on binary output	    */
  char		*comment;   /*  comment text			    */
  char		*history;   /*  history text			    */
} SDS;


/* TYPE: IDS
 * ---------
 * Type definition of IDS structure.
 *
 * The ids definition follows 
 * the ids is a list of sds's.  
 * and  a set of functions which operate on the list 
 * and functions that slice thru the data 
 */
			      /*  maximum SDSs in a chunk -- see chunk list  */
#define MAX_SDS 128

/* CHUNK LIST
 * ----------
 * Each chunk consist of an array of objects of generic type
 * In clist, elemSize determines the size of each object to be
 * stored in the chunks.
 * Chunks are linked in a single linked list manner, order of
 * objects are preserved.
 * The size of each chunk is limited by chunkSize, to be declared
 * by user.
 */

typedef struct _chunk {
  struct _chunk *next;
  int           nElem;
  void          *base;        /* generic type */
} chunk;


typedef struct _clist {
  int   listSize;   /* length of the list */
  int   elemSize;   /* size of each element */
  int   chunkSize;  /* max no. of elements in a chunk */
  chunk *head;      /* points to the first chunk */
} clist;

/*  note:  chunk and clist are TYPEs, not variables  */
/*  INSTANTIATE TO SDSlist  */
typedef clist SDSlist;  /*  this defines what form SDSlist takes on  */

typedef struct {
  SDSlist *sdslist;
  ATTRIBUTES *Attrinfo;
  AXIS *PrimaryAxis;	/* the SDS index axis is implicit in the IDS */ 
} IDS;



/* Type llist supports all linked lists with this structure.
 * The member next must be the first member of the structure.
 */
typedef struct _ll {
  struct _ll *next;
} _llist;

/* Type _llistFn is used as a mapping function, ie. a pointer 
 * to a function.  This, passed to sds_MapList, will be performed
 * on each element in the linked list.
 */
typedef int (*_llistFn)(_llist *ll, void *obj);

/****************************************************************************/
/***************************  FUNCTION PROTOTYPES  **************************/
/****************************************************************************/

					 /*  malloc routines, sds_malloc.c  */
extern SDS *sds_create (void );
extern int sds_create_data (SDS* sdsptr);
extern SDS_STATS *sds_create_stats (void);		/* sds_stats_inf.c */
extern void sds_free (SDS **sdsptrptr);
extern void sds_free_attr (ATTRIBUTES *attr);
extern void sds_free_color (SDS *sdsptr);
extern void sds_free_data (SDS *sdsptr);
extern void sds_free_dim (DIMENSION *dim);
extern void sds_free_filename (SDS *sdsptr);
extern void sds_free_ptrs (SDS *sdsptr);
extern void sds_free_varname (SDS *sdsptr);
extern ATTRIBUTES *sds_malloc_attribute (void);
extern SDS_COLOR *sds_malloc_color (void);
extern char *sds_malloc_data (long numobytes);
extern DIMENSION *sds_malloc_dimension (void);
extern void *sds_malloc_fillvalue (int datatype);
extern int *sds_malloc_length (int rank);
extern char *sds_malloc_string (int numochars);
extern void *sds_std_fillvalue (int sds_type);
							  /*  sds_utility.c  */
	    /*  sds utilities for making sds's more more or less completely  */
extern int sds_append_args_tohist (SDS *sds, argument *args, KEY *list);
extern int sds_append_attrs (SDS *sdsout, SDS *sdsin);
extern int sds_append_comment (SDS *sds, char *str);
extern int sds_append_history (SDS *sds, char *str);
extern SDS *sds_construct (int rank, int *length, int datatype, void *fillval);
extern int sds_copy_attrs (SDS *sdsout, SDS *sdsin);
extern int sds_copy_attr (SDS *sdsout, SDS *sdsin, char *key);
extern SDS *sds_dup_header (SDS *sds);               
extern int sds_fill_data (SDS *sdsptr, void *fillvalue);
extern SDS *sds_make_sds (int rank, int *length, int datatype, void *data); 
extern SDS *sds_replicate (SDS *sds);               
extern SDS *sds_steal_data (SDS *sds);
	      /*  sds_short_template makes sds from sdsin up to data malloc  */
extern SDS *sds_short_template (SDS *sds);
				     /*  debugging tools (in sds_utility.c)  */ 
extern char *sds_datatypename (int datatype);
extern char *xsds (SDS *sds);
extern void DEBUGsds(SDS *sds, char *msg);
							  /*  sds_convert.c  */
extern int sds_data_convert (SDS *sds, int datatype);
extern int sds_scale_data (SDS *sds, double scale, double bias);
extern int sds_set_scaling (SDS *sds, int bits, double scale, double bias);
						 /*  obsolete_sds_convert.c  */
extern int sds_convert (SDS *sds_in, int type_out);
							    /*  sds_slice.c  */
extern SDS *sds_slice (SDS *sds, int *start, int *end);
SDS *sds_slice_file (FILE *fp, int *start, int *end);
		       /*  C++ style access routines into the sds structure  */
							    /*  sds_query.c  */
extern int sds_datatype (SDS *sds);
extern int sds_numbytes (SDS *sds);
extern FILE *sds_stream (SDS *sds);
extern void *sds_data (SDS *sds);
extern void *sds_fillvalue (SDS *sds);
extern int sds_rank (SDS *sds);
extern int sds_dim_n (SDS *sds, int n);		       /*  length of axis n  */
extern int sds_dim0 (SDS *sds);				/*  first dimension  */
extern int sds_dim1 (SDS *sds);			       /*  second dimension  */
extern int sds_dim2 (SDS *sds);				/*  third dimension  */
extern long sds_data_length (SDS *sds);		   /*  length of data array  */
extern SDS_STATS *sds_stats_ptr (SDS *sds);
extern int sds_sizeof (int datatype);			    /*  sds_malloc.c */
extern int *sds_length (SDS *sds);
extern char *sds_filename (SDS *sds);
extern char *sds_varname (SDS *sds);
extern int sds_data_quality (SDS *sds);
extern double sds_bscale (SDS *sds);
extern double sds_bzero (SDS *sds);
							  /*  Axis routines  */
extern AXIS *sds_axisinfo(SDS *sdsptr);
extern AXIS *sds_nth_axis(SDS *sdsptr, int axis_number);
extern void free_axis(AXIS *axisinfo);
							    /*  sds_stats.c  */
extern int sds_calculate_stats (SDS *sdsptr, int *nmissing, int *npresent,
        double *minimum, double *maximum, double *mean,
        double *stdv, double *skew, double *kurt);
							/*  sds_stats_inf.c  */
extern int sds_copy_stats(SDS *sdsin, SDS *sdsout);
extern void sds_stats_info(SDS *sdsptr, int (*message)(const char *, ...));
							      /*  sds_set.c  */
						     /*  set the attributes  */
extern int sds_set_attrvalue (ATTRIBUTES *attr, void *value, int datatype);
extern int sds_set_attrname (ATTRIBUTES *attr, char *attrname);
extern int sds_set_comment (ATTRIBUTES *attr, char *comment);
extern int sds_set_attribute (SDS *sds, char *name, void *value, int attrtype, char *comment);
extern int sds_remove_attribute (SDS *sds, char *name);
extern void sds_set_stdhead (SDS *sds, char *prog_name);
					  /*  set the sds structure members  */
extern int sds_set_data (SDS *sds, void *data);
extern int sds_set_datatype (SDS *sds, int datatype);
extern int sds_set_dim (SDS *sdsptr, int dim, int dimvalue);
extern int sds_set_filename (SDS *sdsptr, char *filename);
extern int sds_set_fillval (SDS *sds, void *fillval);
extern int sds_set_fillvalue (SDS *sds, void *fillvalue);
extern int sds_set_length (SDS *sds, int *length);
extern int sds_set_numbytes (SDS *sds, int numbytes);
extern int sds_set_rank (SDS *sds, int rank);
extern int sds_set_startind (SDS *sds, int *startind);
extern int sds_set_stream (SDS *sdsptr, FILE *fp);
extern int sds_set_varname (SDS *sdsptr, char *varname);
						  /*  set the AXIS routines  */
extern int sds_set_axis(SDS *sdsptr, int *axis_datatype, 
			void *start, void *stepsize,
                    	char **units, char **format, char **title);

extern ATTRIBUTES *sds_attributes(SDS *sds);


/* query the attribute list */
extern ATTRIBUTES *sds_first_attr(SDS *sds);
extern ATTRIBUTES *sds_last_attr(SDS *sds);
extern ATTRIBUTES *sds_next_attr(ATTRIBUTES *attr);
extern ATTRIBUTES *sds_search_attr(SDS *sds, char *key);
		/* sds_attr.c */
extern void *sds_search_attrvalue(SDS *sds, char *key);
extern char *sds_search_attrvalue_str(SDS *sds, char *key);
extern double sds_search_attrvalue_double(SDS *sds, char *key);


extern char *sds_attrname(ATTRIBUTES *attrptr);
extern char *sds_attrcomment(ATTRIBUTES *attrptr);
extern void *sds_attrvalue(ATTRIBUTES *attrptr);
extern char *sds_attrvalue_str(ATTRIBUTES *attrptr);
extern int sds_attribute_type(SDS *sdsptr, char *key);
extern int sds_attrtype(ATTRIBUTES *attrptr);

/* Common keylist interface routines */
extern char *sds_getkey_str(SDS *sds, char *key);
extern char *SDS_getkey_str(SDS *sds, char *key);
extern int sds_getkey_int(SDS *sds, char *key);
extern double sds_getkey_double(SDS *sds, char *key);
extern TIME sds_getkey_time(SDS *sds, char *key);
extern TIME sds_getkey_time_interval(SDS *sds, char *key);
extern SDS_STATS *SDS_getkey_stats(SDS *sds); /* functions/data_stats.c */

extern int sds_setkey_str(SDS *sds, char *key, char *str);
extern int sds_setkey_int(SDS *sds, char *key, int val);
extern int sds_setkey_double(SDS *sds, char *key, double val);
extern int sds_setkey_time(SDS *sds, char *key, TIME time);
extern int sds_setkey_time_interval(SDS *sds, char *key, TIME time);
extern int sds_setkey_stats(SDS *sds, SDS_STATS *stats); /* functions/data_stats.c */
						       /*  SDS I/O routines  */
							    /*  general I/O  */
							       /*  sds_in.c  */
extern SDS *sds_in (char *infile, int type, int filetype);
extern SDS *sds_read (FILE *infp, int type, int filetype);
							      /*  sds_out.c  */
extern int sds_out (SDS *sdsptr, char *outfile, int filetype); 
extern int sds_write (SDS *sds, FILE *outfp, int filetype);

extern SDS *sds_read_sds(FILE *in);
extern int sds_write_sds(SDS *in, FILE *out);

extern SDS *sds_read_raw(FILE *in, int datatype);
extern int sds_write_raw(SDS *in, FILE *out);

							     /*  sds_fits.c  */
extern SDS *sds_get_fits (char *filename);
extern SDS *sds_read_fits (FILE *in);
extern SDS *sds_get_fits_head (char *filename);
extern SDS *sds_read_fits_head (FILE *in);
extern int sds_read_fits_data (SDS *sds);
extern int sds_put_fits (SDS *sds, char *filename);
extern int sds_write_fits (SDS *sds, FILE *out);
extern int sds_write_fits_head (SDS *sds, FILE *out);

							      /*  sds_gif.c  */
extern int sds_write_gif (SDS *sds, char *filename, int interlace,
    double scale, double offset, int transparent, int delay);
					      /*  sds_rfits.c : obsolescent  */
/*
extern SDS *sds_get_one_fits(char *filename);
extern int sds_get_fits_header(SDS *sdsptr, FILE *fp);
extern int sds_get_fits_data(SDS *sds, FILE *fp);
extern SDS *sds_read_FITS_header(char *filename);
*/
							     /*  sds_flip.c  */
extern int sds_flip_data(SDS *sdsptr);
extern int sds_flip(void *data, long length, int num_bytes);
					      /*  sds_wfits.c : obsolescent  */
/*
extern int sds_put_one_fits(SDS *sdsptr, char *outfile);
extern int sds_put_fits_header(SDS *sdsptr, FILE *fp);
extern int sds_put_fits_data(SDS *sdsptr, FILE *fp);
*/

/* query the stats structure  */

extern int sds_stats_valid (SDS *sdsptr);
extern char *sds_min (SDS *sdsptr);
extern char *sds_max (SDS *sdsptr);
extern char *sds_min_abs (SDS *sdsptr);
extern char *sds_mean (SDS *sdsptr);
extern char *sds_skewness (SDS *sdsptr);
extern char *sds_stdev (SDS *sdsptr);
extern char *sds_kurtosis (SDS *sdsptr);
extern char *sds_numvals (SDS *sdsptr);
extern void sds_update_stats (SDS *sdsptr, SDS_STATS *statptr);
							  /*  IDS Functions  */
						/*  User Interface routines  */
extern IDS *ids_init( void );
extern void ids_free(IDS **idsptrptr);
extern SDS *ids_nth_sds(IDS *ids, int n);
							   /*  ids_series.c  */
extern int ids_get_FITS_series_headers(IDS *ids,char *filename[],int numFiles);

extern IDS *ids_slice_FITS_series_file(IDS *ids, int *start, int *end);
extern SDS *ids_merge_series(IDS *ids, int start, int end);
						/*  Query the ids structure  */
extern int ids_1st_non_missing_rec(IDS *idsptr);
extern int ids_length(IDS *ids);
extern int ids_rank(IDS *idsptr);
extern int *ids_axis_lengths(IDS *idsptr);
extern int ids_dimn(IDS *idsptr, int axis_number);
						   /*  Low level operations  */
extern int ids_add_sds(IDS *ids, SDS *sds);
extern int ids_insert_sds(IDS *ids, SDS *sds, int n);
extern int ids_delete_sds(IDS *ids, int n);

    /*  these are extremely low-level routines not for general consumption   */
/* CHUNK LIST Helper Functions 
 * ---------------------------
 */
extern clist *CL_InitList(int elemSize, int chunkSize);
extern void *CL_DataAt(clist *list, int pos);
extern int CL_Length(clist *list);
extern void CL_FreeList(clist *list);
extern int CL_InsertAt(clist *list, void *obj, int pos);
extern int CL_DeleteAt(clist *list, int pos);

/*  the next routines are for the support of the linked lists  */
extern void sds_AddList(_llist **l, _llist *new);
extern void sds_FreeList(_llist **l);
extern int sds_ListSize(_llist **l);
extern _llist *sds_NextElement(_llist **l);
extern _llist *sds_LastElement(_llist **l);
extern _llist *sds_FirstElement(_llist **l);
extern _llist *sds_GetNthElement(_llist **l, int n);
extern _llist *sds_RemoveNthElement(_llist **l, int n);
extern _llist *sds_DeList(_llist **l);
extern _llist *sds_InsertNthElement(_llist **l, _llist *new, int n);
extern _llist *sds_MapList(_llist **l, _llistFn fn, void *obj);
extern char *ConvertToUpper(char *line);
extern char *ConvertToLower(char *line);
extern char *FirstNonWhite(char *line);
extern char *LastNonWhite(char *line);
extern int numType(char *num);
extern int sds_type(char *str);
extern void g_max_min(void *base, long n, int nbytes,
	       int (*cmp)(void *, void *),
	       void *maxptr, void *minptr);
#endif

/*
 *  Revision History
 *
 *  96.11.26	R Bogart	added advisory scaling information to sds
 *	struct for future use
 *  96.12.08	R Bogart	changed name of scale info members (the B in
 *	BSCALE stands for "brightness", not "binary"); added prototype for
 *	xsds()
 *  97.03.07	P Scherrer	added sds_getkey_stats and sds_setkey_stats
 *  97.04.16	K Chu		added #include <soi_version.h>
 *  98.01.15	R Bogart	removed obsolete declarations; removed change
 *	history prior to v 2.0 (96.11.26); added declarations for new scaling
 *	and conversion routines and new FITS interface
 *  98.04.16	R Bogart	changed version number to string; added
 *	declaration for sds_put_fits_head and removed declarations from
 *	sds_rfits.c and sds_wfits.c
 *  98.05.03	R Bogart	included history and comment members in sds
 *	struct, declared append functions for them
 *  99.01.02	K Leibrand	added minabs to SDS_STATS in implementing
 *	fix of SSTR 27
 *  99.07.23	R Bogart	modified prototype for sds_gif funct
 *  01.07.24	R Bogart	added prototype for sds_append_args_tohist()
 *  02.07.08	R Bogart	removed declarations commented out 98.01;
 *	removed RCS comments below from before 2000 (all info in above)
 */

/*
$Id: soi_sds.h,v 1.1 2009/04/24 21:53:17 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/soi_sds.h,v $
$Author: production $
*/
/* $Log: soi_sds.h,v $
 * Revision 1.1  2009/04/24 21:53:17  production
 * *** empty log message ***
 *
 * Revision 1.69  2002/07/08  21:36:56  rick
 * see above
 *
 * Revision 1.68  2001/07/24  23:15:47  rick
 * see above
 *
 */
