/*
 * dsds.h
 *
 * Contains the definitions needed for the dsds_svc module.
 *
*/
#ifndef DSDS_INCL

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

#include <soi_key.h>
#include <sum_rpc.h>

#define DSDS_VERSION_NUM    (0.8)

#define DSDSMASTER "sonar.Stanford.EDU"	/* location of global dsds/ampex_svc */

#define MAX_PART 64	/* max # +1 of dedicated DSDS partitions */

#define NUM_GROUPS 256

#define MAX_RETRIEVE 196 /* max # of simultaneous lago/ampex_svc retrievals */


/* Note: Do NOT change the following. They are in the database */
#define DARW 1		/* data allocate assigned as read/write */
#define DADP 2		/* data allocate deletion pending when effective date */
#define DAAP 4		/* data allocate archive pending */
#define DARO 8		/* data request assigned as read only */
#define DAPERM 16	/* data allocate is permanent */
#define DASAP 32	/* data allocate safe archive pending */
/* Here are sub-statuses for DAAP */
#define DAAEDDP 32	/* don't archive and when effective date mark DADP */
#define DAAPERM 64	/* archive then mark DAPERM */
#define DAADP 128	/* archive then mark DADP when effective date */
/* Here are sub-statuses for DADP */
#define DADMVA 256	/* active ds move, don't delete wd */
#define DADMVC 512	/* ds has been moved, don't mark offline when rm wd */

/* Pe/uid assignment table. One of these is put onto the peuid_hdr pointer
 * each time a pe registers (opens) with dsds_svc, and removed when pe
 * deregisters (closes).
*/
//now in sum_rpc.h
//struct peuid {
//  struct peuid *next;
//  long uid;
//  int petid;
//};
//typedef struct peuid PEUID;

/* Partition assignment table. For working directory assignments made 
 * by dsds_svc, there will be a number of these tables linked onto one of 
 * the pahdr_xx pointers.
*/
//now in sum_rpc.h
//struct padata {
//  struct padata *next;
//  char *wd;
//  char *effective_date;
//  long uid;
//  double bytes;
//  int status;
//  int archsub;		/* archive pend substatuses */
//  int group_id;		/* for grouping in tape archives */
//  int safe_id;		/* for grouping in safe tape archives */
//  unsigned long ds_index;
//};
//typedef struct padata PADATA;

/* Partition definition table. One for each dedicated DSDS partition.
 * Initialized by dsds_svc from the partn_avail data base table. */
//struct partition {
//  char *name;		/* name of the partition */
//  double bytes_total;	/* total number of bytes of the partition */
//  double bytes_left;	/* bytes unassigned */
//  double bytes_alloc;	/* bytes allocated by DS_Allocate() */
//  int pds_set_num;	/* which /PDS set the partition belongs to */
//};
//typedef struct partition PART;

/* This is a linked list of keylist pointers and their status.
 * Saves the keylist here and the status returned by a query
 * to DS_DataRequest() for each dataset passed to DataRequestAll() for a pe.
*/
struct key_link {
  struct key_link *next;
  int status;
  KEY *list;
};
typedef struct key_link KEY_LINK;

/* Used by DataRequestAll() for tape retrieval for a pe.
 * An array of these end up sorted by increasing tapeid and filenum.
*/
struct tapesort {
  int tapeid;
  int filenum;
  int status;
  KEY *list;
};
typedef struct tapesort TAPESORT;

/* An array of these for all the pe's that have an active retrieval request
 * with DataRequestAll().
*/
struct pe_active {
  long dsds_uid;        /* unique id for each pe. 0 if free */
  int retnum;		/* # of retrievals to do */
  int retcurrent;	/* current retrieval # (0 to retnum-1) */
  KEY *origlist;	/* original keylist to DataRequestAll */
  TAPESORT *tapeptr;    /* the sorted tape list for ampex retrieval */
  KEY_LINK *keylink;    /* start of the key_link chain for this pe */
};
typedef struct pe_active PE_ACTIVE;


//typedef struct tape {
//   struct tape *next;
//   int  barcode;	/* lago barcode is barcode number on Exabyte tape
//		 		 * ampex barcode is combination of barcode on
//		 		 * ampex cartridge and tape partition number
//		 		 */
//   int  lastfn;		/* file number of the last successfully written
//		 		 * file on the tape 
//		 		 */
//   int  group_id;	/* integer between 0 and NUM_GROUPS which defines
//		 		 * which datasets can be archived on this tape
//		 		 */
//   long availblocks;/* number of available blocks remaining on the tape
//		 		 */
//} TAPE;

/* Defines each /PDS set. 
 * Initialized by src/libpe.d/pds_set_def.c.
*/
struct pds_set_def {
  char *host_name;	/* host that has the /PDS local */
  char *first_parti;	/* 1st partition in the set, e.g. "/PDS0" */
  char *last_parti;	/* last partition in the set, e.g. "/PDS19" */
  int pds_set_num;	/* unique number for the entire set */
};
typedef struct pds_set_def PDS_SET_DEF;

#ifdef __sgi
/* Pthread variables for thread syncronization in dsds_svc */
typedef struct p_struct_tag {
  pthread_mutex_t mutex;
  pthread_cond_t  cond;
  KEY *list;
  int req_num;
  int callertid; 
} P_STRUCT_T;
#endif




/* utilities for handling archive group id's */
extern int is_arch (int gpid);
extern int is_safe (int gpid);
extern int is_export (int gpid);
extern int get_gpid(char *prog, char *level, char *series);
extern int get_safe_gpid(char *prog, char *level, char *series);
extern int get_rgpid(char *prog, char *level, char *series);
extern char *get_export_site (int xgpid);

/* utilities for handling DataRequestAll() for REQREQALL call */
extern void ret_key_all(KEY_LINK *link, TAPESORT *tapeptr, KEY **retlist, int num);
extern TAPESORT *sort_tapes(KEY_LINK **link, int num);
extern int compare_tape(const void *a, const void *b);
extern int ret_key_link(KEY_LINK *link, KEY **retlist);
extern struct RETTAB *get_rettab(KEY *list, int *status);
extern int read_ds(struct RETTAB *ptab);
extern KEY *read_ds_ack(struct RETTAB *ptab, int *stat);
extern int startret(struct RETTAB *pret, KEY *list, PE_ACTIVE *peaptr);
extern int retrieve_open_lago(struct RETTAB *pret, int tapeid, int reqnum);
extern void log_key_link(KEY_LINK *link);
extern void free_all(PE_ACTIVE *peaptr);
extern void free_tapeptr(TAPESORT *tapeptr, int num);
extern void reqreqstat(int offflg, int (*log)(char *fmt, ...));

#ifdef PVM33	/* for pvm version 3.3 and above */
#define hostinfo pvmhostinfo
#define taskinfo pvmtaskinfo
#define pvm_serror(x) pvm_setopt(PvmAutoErr,x)
#endif		/* PVM33 */
 
#define DSDS_INCL
#endif
/*
 * Revision History
 * v0.3	93.08.09	ja	change from Mb to bytes of storage
 * v0.2	93.08.04	ja
 * v0.1	93.07.26	ja
 *      created original file
 *
*/

/*
$Id: dsds.h,v 1.1 2009/04/24 21:50:28 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/dsds.h,v $
$Author: production $
*/
/* $Log: dsds.h,v $
 * Revision 1.1  2009/04/24 21:50:28  production
 * *** empty log message ***
 *
 * Revision 1.34  2003/05/15  17:14:50  jim
 * change DSDSMASTER from tarax to sonar. Tarax has gone away. This may
 * not be entirely necessary as it's just used for spawning ampex in
 * local (-L) mode and tarax is now a sonar alias...
 *
 * Revision 1.33  2002/11/22  21:33:55  jim
 * add edate to DS_UpArchsub()
 *
 * Revision 1.32  2000/05/16 22:46:34  CM
 * add #ifdef __sgi
 *
 * Revision 1.31  2000/04/27 23:54:36  jim
 * add P_STRUCT_T
 *
 * Revision 1.30  2000/03/15  00:30:20  jim
 * add sub-statuses for DSDP and update protocols
 *
 * Revision 1.29  1999/09/07  18:46:15  jim
 * add struct pds_set_def and pds_set_num in struct partition
 *
 * Revision 1.28  1998/12/10  01:07:53  jim
 * add bytes_alloc to PART
 *
 * Revision 1.27  1998/08/03  22:08:14  kay
 * add get_rgpid, delete search_gpid, match
 *
 * Revision 1.26  1998/06/16  22:14:03  kehcheng
 * added PVM33 macros for pvm versions 3.3 and above
 *
 * Revision 1.25  1997/04/16 21:37:34  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.24  1997/04/07  23:40:51  jim
 * add reqreqstat()
 *
 * Revision 1.23  1997/01/19  19:11:49  jim
 * adjust for REQREQALL
 *
 * Revision 1.22  1997/01/14  19:24:03  jim
 * update to support REQREQALL
 *
 * Revision 1.21  1996/09/24  23:16:36  kay
 * changed *Stat* and CatDelete to use ds_index instead of wd and uid
 *
 * Revision 1.20  1996/09/23  17:06:31  kay
 * added ds_index to PADATA, DS_PaUpdate, NC_PaUpdate
 * added functions NC_PaAdd, NC_PaDelete, NC_PaReplace
 *
 * Revision 1.19  1996/08/29  19:04:32  kay
 * added NC_PaRequest_AP
 *
 * Revision 1.18  1996/08/02  00:49:51  kay
 * added NC_ExportRequest
 *
 * Revision 1.17  1996/07/18  21:11:22  kay
 * added define for DASAP, added status and safe_id to struct padata
 * added gpid to arguments of *PaUpdate and NC_StatArchive
 * added new functions NC_PaRequest, NC_SafeID, and group id utilities
 * deleted externs for *PaUpdate_AP
 *
 * Revision 1.16  1996/06/13  21:00:50  jim
 * change history prototype
 *
 * Revision 1.15  1996/03/29  17:33:53  jim
 * change meaning of DAAEDDP
 *
 * Revision 1.14  1996/01/04  16:14:35  jayasree
 * added extern TAPE *NC_TapeByDate
 *
 * Revision 1.13  1995/12/28  17:31:47  kay
 * added next pointer to TAPE structure
 *
 * Revision 1.12  1995/12/01  19:53:25  jim
 * changed DSDSMASTER to tarax
 *
 * Revision 1.11  1995/11/07  16:53:39  jim
 * added DSDSMASTER
 *
 * Revision 1.10  1995/08/04  17:12:31  jim
 * changed MAX_PART from 32 to 64
 *
 * Revision 1.9  1995/06/19  15:56:09  jayasree
 * added group_id to tape structure.
 *
 * Revision 1.8  1995/04/20  23:00:34  jim
 * made all prototype bytes a double.
 *
 * Revision 1.7  1994/09/02  22:56:48  jim
 * put back lago.h include
 *
 * Revision 1.6  1994/09/02  22:53:36  jim
 * removed include <lago.h>
 *
 * Revision 1.5  1994/09/02  20:41:34  jim
 * add DS_ConnectDB and DS_DisConnectDB prototypes
 *
 * Revision 1.4  1994/08/30  23:10:01  jayasree
 * add eff_date and archsub to DS_PaUpdate
 *
 * Revision 1.3  1994/08/29  17:29:34  jayasree
 * removed dbname from DS_ functions prototyps.
 *
 * Revision 1.2  1994/08/29  17:12:39  jim
 * added effective date and archsub to padata
 *
 * Revision 1.1  1994/02/16  23:21:14  CM
 * Initial revision
 * */
