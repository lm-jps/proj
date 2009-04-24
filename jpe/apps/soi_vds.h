/*
 *  soi_vds.h				$(DIST)/include/soi_vds.h
 *
 *   Based on code by Phil Scherrer for handling multivariant data
 *   and dataset names
 *
 *  Bugs: yes
 *
 *  Planned updates:
 *
 */

#ifndef SOI_VDS_INCL 

/****************************************************************************/
/**************************  INCLUDE STATEMENTS  ****************************/
/****************************************************************************/

#ifndef SOI_VERSION_INCL
#include "soi_version.h"
#endif

#include "soi_sds.h"
#include "soi_NaN.h"
#include "soi_key.h"
#include "soi_str.h"
#include "soi_at.h"

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define VDS_CADENCE	(1)

#define VDS_READ	(1)
#define VDS_WRITE	(2)
#define VDS_APPEND	(3)

#define VDS_INFO	(1)
#define VDS_DATA	(2)

#define DEFAULT_PROTOCOL	"RDB.FITS"

#define VDS_CONFORMANCE_VERSION	"2003.11.18"

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/

typedef struct name_list {
  struct name_list *next;
  char		*name;
} NAMELIST;

/* add slicing info or keep the keylist around?  */
typedef struct vds_description {
  char		*wd;			/* dataset directory name */
  char		*bn;			/* basename for filename building */
  char		*fmt;			/* format of record part of datafiles */
  char		*ov;			/* Full pathname of overview file */
  char		*info;			/* */
  char		*prog;			/* prog name or NULL if not dataset */
  char		*filename;		/* for single file or datacube */
  int		 fsn;			/* requested first record */
  int		 lsn;			/* requested last record */
  int		 first_record;		/* actual first record */
  int		 last_record;		/* actual last record */
  int		 mode;			/* read write or ? append ? */
  int		 vds_extern_type;	/* external data format */
  int		 vds_intern_type;	/* internal data format */
  void		*vds_pointer;		/* for record info and fits records */
  SDS		*sds;			/* for data cube */
  int		state;			/* Flag for state of vds. */
  SDS 		*global_attributes;
  int		 nvars;
  NAMELIST	*var_names;
} VDS;

typedef int VDS_TYPE;


/****************************************************************************/
/*********************  GLOBAL & EXTERN DECLARATIONS  ***********************/
/****************************************************************************/


/* state flag to control known state of dataset. */

enum state_types {
VDS_EMPTY=0,	/* after vds_create */
VDS_CLEAN,	/* after vds_opne or vds_new */
VDS_CUBED,	/* after convert_to_cube but before use */
VDS_USED,	/* after data inserted or selected */
		/* new state types here */
VDS_NUM_STATES};

/* index used to process different types of conforming datasets */

enum conformance_types {
VDS_NOT_CONFORMING,	/* =0, must be first */
VDS_TS,			/* time series */
VDS_TS_EQ,		/* time series - equally spaced */
VDS_TS_BLOCKED,		/* time series - blocked */
VDS_L_NU,		/* l-nu diagram */
VDS_MDICAL,		/* calibration */
VDS_MISC,
VDS_TS_EQ_TILE,		/* ts_eq but sn units not T_BLOCK */
			/* insert additional types here,
                           edit conformance_names accordingly */
VDS_NUM_TYPES};		/* count of conformance types, must be last */


/* descriptions of conforming datasets indexed by conformance type */
/* vds attribute value for the key CONFORMS should be one of these names */

extern char *conformance_names[VDS_NUM_TYPES];	/* definition in libvds.d */


/* index which enumerates the way data and attributes are stored externally */

enum protocol_types {
VDS_NO_PROTOCOL,	/* =0, must be first */
VDS_CDF,
VDS_FITS,
VDS_FITS_MERGE,
VDS_FITS_TABLE,
VDS_RDB,
VDS_NUM_PROTOCOLS};	/* count of external protocol types, must be last */


/* info arrays corresponding to protocol types used to determine protocols
   and generate names - definitions in libvds.d */

extern char *protocol_names[VDS_NUM_PROTOCOLS];
extern char *protocol_suffixes[VDS_NUM_PROTOCOLS];


/* index which enumerates the way data and attributes are stored internally */
/* use mutually exclusive index sets for external and internal protocols */

enum internal_protocol_types {
VDS_START_INTERNAL=VDS_NUM_PROTOCOLS,
VDS_CDF_INTERNAL,
VDS_FILE,
VDS_CUBE,
VDS_FILES,
VDS_END_INTERNAL};

#define VDS_NO_INTERNAL_PROTOCOL 0
#define VDS_FIRST_INTERNAL_PROTOCOL VDS_START_INTERNAL+1
#define VDS_NUM_INTERNAL_PROTOCOLS VDS_END_INTERNAL-VDS_START_INTERNAL

/****************************************************************************/
/****************************  MACRO DEFINITIONS  ***************************/
/****************************************************************************/

#define    VDS_QUIT(errcode) {soi_errno=errcode; return(NULL);}



/****************************************************************************/
/***************************  FUNCTION PROTOTYPES  **************************/
/****************************************************************************/

/* vds constructors */
extern VDS *vds_open (KEY *keylist, char *root_key);
extern VDS *vds_open_all (KEY *keylist, char *root_key);
extern VDS *vds_new (KEY *keylist, char *root_key, int protocol);
extern VDS *vds_create_std(KEY *keylist, char *root_key, int mode);
extern VDS *vds_create(void);	/* create an empty VDS structure */
					    /*  vds destructors: vds_new.c  */
extern int vds_close (VDS **vdsptrptr);
extern int vds_flush_data (VDS *, int var);
extern int vds_flush_data_rec (VDS *, int var, int rec);
extern void vds_free(VDS **vdsptrptr);

/* data insertion */
extern int vds_insert_sds (VDS *vds, SDS *sds, int var_num);
extern int vds_insert_sds_rec (VDS *vds, SDS *sds, int var_num, int rec_num);

/* data selection */
extern SDS *vds_select (VDS *vds, int var_number, int recStart, int recCount, 
                  int recInterval, int *indices, int *counts, int *intervals);

extern SDS *vds_select_hdr (VDS *vds, int var_number, int rec_number);
extern SDS *vds_select_rec(VDS *vds, int var_number, int rec_number);
extern SDS *vds_select_drec(VDS *vds, int var_number, int rec_number);
extern SDS *vds_select_frec(VDS *vds, int var_number, int rec_number);

extern SDS *VDS_select_hdr (VDS *vds, int var_number, int rec_number);
extern SDS *VDS_select_rec(VDS *vds, int var_number, int rec_number);
extern SDS *VDS_select_drec(VDS *vds, int var_number, int rec_number);
extern SDS *VDS_select_frec(VDS *vds, int var_number, int rec_number);

extern SDS *vds_select_key_double(VDS *vds, int var, char *key);
extern SDS *vds_select_vars (VDS *vds, int n_vars, int *var_numbers);

/* data conversion */
extern int vds_convert_to_cube (VDS *vdsptr, int nrecs);

/* global attributes */
extern int vds_conform_type(VDS *vds);
extern int vds_protocol_type (VDS *vds, int which);

extern int vds_check_info(VDS *vds);
extern int vds_copy_conform_info(VDS *vdsout, VDS *vdsin);
extern int vds_set_conform_info(VDS *vds, char *conforms, char *protocol,
	char *t_epoch, char *t_block, char *ext_info, char *series_type);

extern char *vds_getkey_str(VDS *vds, char *key);
extern char *VDS_getkey_str(VDS *vds, char *key);
extern int vds_getkey_int(VDS *vds, char *key);
extern double vds_ketkey_double(VDS *vds, char *key);
extern TIME vds_getkey_time(VDS *vds, char *key);
extern TIME vds_getkey_time_interval(VDS *vds, char *key);

extern int vds_setkey_str(VDS *vds, char *key, char *val);
extern int vds_setkey_int(VDS *vds, char *key, int val);
extern int vds_setkey_double(VDS *vds, char *key, double val);
extern int vds_setkey_time(VDS *vds, char *key, TIME val);
extern int vds_setkey_time_interval(VDS *vds, char *key, TIME val);

extern char *vds_getkey_rec_str(VDS *vds, int var, int rec, char *key);
extern char *VDS_getkey_rec_str(VDS *vds, int var, int rec, char *key);
extern int vds_getkey_rec_int(VDS *vds, int var, int rec, char *key);
extern double vds_getkey_rec_double(VDS *vds, int var, int rec, char *key);
extern TIME vds_getkey_rec_time(VDS *vds, int var, int rec, char *key);
extern TIME vds_getkey_rec_time_interval(VDS *vds, int var, int rec, char *key);

extern int vds_attribute_type(VDS *vds, char *key);
extern int vds_set_attribute(VDS *vds, char *name, void *value, 
			     int attrtype, char *comment);
extern int vds_remove_attribute(VDS *vds, char *name);
extern int vds_copy_attribute(VDS *vdsout, VDS *vdsin, char *key);

extern void *vds_search_attrvalue(VDS *vds, char *key);
extern char *vds_search_attrvalue_str(VDS *vds, char *key);
extern double vds_search_attrvalue_double(VDS *vds, char *key);

extern int vds_check_record_file(char *filename, int conformance_type);
extern int vds_check_records(AT *recinfo, int conformance_type);

/* routines that set the structure  */
extern int vds_set_last_record(VDS *vdsptr, int last_record);
extern int vds_set_vds_ptr(VDS *vdsptr, void *vds_ptr);
extern int vds_set_intern_ptr(VDS *vdsptr, void *intern_ptr);
extern int vds_set_extern_type(VDS *vdsptr, int vds_extern_type);
extern int vds_set_intern_type(VDS *vdsptr, int vds_intern_type);
extern int vds_set_filename(VDS *vdsptr, char *filename);
extern int vds_set_var_namelist(VDS *vdsptr, NAMELIST *var_names);

/* routines that query the vds  */
extern int vds_nrecs(VDS *vds);
extern int vds_last_record(VDS *vds);

extern int vds_extern_type(VDS *vds);
extern int vds_intern_type(VDS *vds);
extern void *vds_intern_ptr(VDS *vds);

extern int vds_var_rank(VDS *vds, int var_number);
extern int vds_var_dimn(VDS *vds, int var_number, int axis_number);
extern int *vds_var_lengths(VDS *vds, int var_number);
extern int vds_var_datatype(VDS *vds, int var_number);
extern void *vds_var_fillvalue(VDS *vds, int var_number);
extern void *vds_default_fillvalue(int datatype);

extern char *vds_wd(VDS *vds);
extern char *vds_get_info_filename(VDS *vds, int var);
extern char *vds_get_data_filename (VDS *vds, int var);

extern char *vds_title (VDS *vds);

extern void DEBUGvds(VDS *vds, char *msg);

/* variable routines  */
extern int vds_nvars (VDS *vds);
extern char *vds_var_name (VDS *vds, int var_number);
extern int vds_add_var_name (VDS *vds, char *varname);
extern int vds_var_number (VDS *vds, char *var_name);

/* AXIS routines */
extern AXIS *vds_axis (VDS *vds);
extern int vds_promote_axis (VDS *vds);

							     /*  vds_util.c  */
extern VDS *vds_get_sn_range (KEY *params, char *name, int ds,
    int *fsn, int *lsn, int *dsn, int discontinuous);
extern int vds_record_count (KEY *params, char *vds_name, int discontinuous,
    int *ds_count);

#define SOI_VDS_INCL
#endif

/* Revision History
 *
 *  v 2.8  98.02.18	Rick Bogart	declared obsolete versions of vds_close
 *	and vds_flush_data_rec for compatability; removed history prior to v 1
 *  v 4.3  99.06.10	R Bogart	added declarations from vds_util
 *  v 4.5  99.10.13	R Bogart	modified "
 *  v 4.8  01.07.27	R Bogart	removed declarations of unused old_vds*
 *	from v 2.8
 *  v 5.4  03.11.18	R Bogart	modified conformance version date to
 *	reflect tightened rules
 */

/*
$Id: soi_vds.h,v 1.1 2009/04/24 21:53:31 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/soi_vds.h,v $
$Author: production $
 * $Log: soi_vds.h,v $
 * Revision 1.1  2009/04/24 21:53:31  production
 * *** empty log message ***
 *
 * Revision 1.50  2003/11/18  23:15:08  rick
 * see above
 *
 * Revision 1.49  2001/08/28 00:10:15  rick
 * see above
 *
 * Revision 1.48  1999/12/14  00:30:53  rick
 * see above
 *
 * Revision 1.47  1999/07/21  22:31:26  rick
 * see above
 *
 * Revision 1.46  1999/01/20 00:48:57  kay
 * added functions vds_check_records and vds_check_record_file
 * to implement fix to SSTR 18
 *
 * Revision 1.45  1998/02/19  00:32:22  rick
 * see above
 *
 * Revision 1.44  1997/12/15  18:17:17  phil
 * added more vds_getkey_rec_XXX types
 *
 * Revision 1.43  1997/05/24  01:28:14  phil
 * added TS_EQ_TILE conforming type.
 *
 * Revision 1.42  1997/04/16  21:58:17  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.41  1996/11/13  23:05:49  kay
 * recategorized function externs, added extern for DEBUGvds,
 * moved externs for GETKEY_str and getkey_time_interval to soi_key.h
 * moved externs for sds_getkey*, sds_setkey*, SDS_getkey_str, and DEBUGsds
 * to soi_sds.h
 *
 * Revision 1.40  1996/10/09  22:52:30  kay
 * mutually exclusive index sets for external and internal protocols
 *
 * Revision 1.39  1996/07/10  23:04:54  phil
 * added vds_flush_data_rec
 *
 * Revision 1.38  1996/06/03  20:46:28  phil
 * added vds_open_all
 *
 * Revision 1.37  1996/05/31  00:28:42  phil
 * *** empty log message ***
 *
 * Revision 1.36  1996/05/06  21:04:19  kay
 * added sds field to VDS structure and added extern for vds_flush_data
 *
 * Revision 1.35  1996/04/19  23:28:30  kay
 * redefined internal_protocol_types
 *
 *
 */
