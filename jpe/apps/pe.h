/*
 * pe.h
 *
 * Contains the definitions needed for the Pipeline Execution module
 *
*/
#ifndef PE_INCL

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

#include "dsds.h"

#define PE_VERSION_NUM    (0.9)

#define SUMSERVER "d02.Stanford.EDU"

#define PELOGFILE "/tmp/pe.%s.%d.log"
#define PEMAILFILE "/tmp/pe.%s.%d.mail"

#define MAX_STR 256	/* max size of a string */
#define MAX_SERV 128	/* max # of different server progs, e.g. remap_svc */
#define MAX_ARGS 64	/* max # of entries in argument list from a server */
#define MAX_ARG_STR 256	/* max size of an arg string */
#define MSGTIMEOUT 0	/* simulated msg # for a timeout */
#define MSGUNREGPERPC  241
#define MSGAMPEX  242
#define MSGAMPEX1 243
#define MSGAMPEX2 244 
#define MSGKILL 245	/* msg number to tell a server to exit */
#define MSGTAEL 246	/* msg number for msg from lago servers to lui */
#define MSGLAGOA 247	/* msg number for requests to lagoa_svc */
#define MSGLAGOB 248	/* msg number for requests to lagob_svc */
#define MSGTAEDS 249	/* msg number for msg from dsds_svc to tae peui */
#define MSGTAEPE 250	/* msg number for msg from pe to tae peui */
#define MSGLAGO 251	/* msg number for requests to lago_svc */
#define MSGDSDS 252	/* msg number for requests to dsds_svc */
#define MSGDSDR 253	/* msg number for request to dsdr_svc (obsolete) */
#define MSGMSG  254	/* msg number to specify a msg number to a server */
#define MSGARGS 255	/* msg number for send/receive of arg list */
#define REQOPEN 1	/* DSDS Open from pe to dsds_svc */
#define REQREQ 2	/* DSDS Data_Request from pe to dsds_svc */
#define REQUPDATE 3	/* DSDS Cat_Update from pe to dsds_svc */
#define REQCLOSE 4	/* DSDS Close from pe to dsds_svc */
#define REQALLOC 5	/* DSDS Allocate from pe to dsds_svc */
#define REQDEALLOC 6	/* DSDS Deallocate from pe to dsds_svc */
#define REQTAEON 7	/* DSDS start send log msgs to tae peui */
#define REQTAEOFF 8	/* DSDS end send log msgs to tae peui */
#define REQIMGUPDATE 9	/* DSDS update image info from cvtlm to dsds_svc */
#define REQDBNAME 10	/* database name request */
#define TAPEARCH 11	/* to ampex server to start an archive set */
#define TAPEWRT 12	/* to ampex server to write an archive */
#define TAPERD 13	/* to ampex server to read an archive */
#define TAPERET 14	/* to ampex server to start a retrieve */
#define TAPEDONE 15	/* to ampex server to end an archive set */
#define TAPEWRTDONE 16	/* to ampex server to indicate write completion */
#define TAPERDDONE 17	/* to ampex server to indicate read completion */
#define LAGOA 18	/* from lagoa server to lui to show info */ 
#define LAGOB 19 	/* from lagob server to lui to show info */ 
#define LAGOINFO 20	/* from lago server to lui to show info */
#define LAGOINIT 21	/* from lago server to lui to give init status */
#define LAGOREQTAPE 22	/* from lago server to lui to request tape */
#define LAGOREQBLANK 23	/* from lago server to lui to request blank tape */
#define LAGOTAEON 24	/* from lui to lago servers to start log msgs */
#define LAGOTAEOFF 25	/* from lui to lago servers to stop log msgs */
#define REQREQNR 26	/* DSDS Data_Request to dsds_svc with no retrieval */
#define REQILOGU 27	/* DSDS ingest_hrsc_log update fr ingest_hrsc_log_svc*/
#define REQMLOGU 28	/* DSDS merge_hrsc_log update fr merge_hrsc_log_svc*/
#define LAGOINVENTORY 29  
#define LAGOOPENDOOR 30  
#define LAGOCLOSEDOOR 31
#define LAGOACCEPT 32  
#define LAGOREMOVE 33 
#define LAGONOSLOT 34
#define LAGONOBARCODE 35
#define LAGONOTINDB 36 
#define REQSQL 37	/* DSDS dynamic sql command */
#define LAGOKILL 38
#define REQ5KLOGU 39	/* DSDS ingest_5k_log update fr ingest_5k_log_svc*/
#define AMPEXINIT       40
#define AMPEXINFO       41
#define AMPEXKILL       42
#define AMPEXTAEON      43
#define AMPEXTAEOFF     44
#define AMPEXOPENDOOR   45
#define AMPEXCLOSEDOOR  46
#define AMPEXACCEPT     47
#define AMPEXREMOVE     48
#define AMPEXINVENTORY  49
#define AMPEXREQBLANK   50
#define AMPEXREQTAPE    51
#define AMPEX1		52
#define AMPEX2		53
#define REQUPDATE0 	54	/* DSDS Cat_Update0 from pe to dsds_svc */
#define REQAPUP 	55 /* DSDS arch pend change from ampexarc to dsds_svc */
#define TAPEIN		56	/* to ampex server to insert cartridge */
#define TAPEINDONE	57	/* to ampex server to indicate that a cartridge
has been inserted */
#define REQRMRUN	58	/* run msg to dsds_rm */
#define REQRM		59	/* msg to dsds_svc to find storage for rm */
#define REQALLOCUP	60	/* msg to dsds_svc to update in-mem lists */
#define TAPERESERVE	61	/* to ampex server to reserve unit */
#define TAPERELEASE	62	/* to ampex server to release reserved unit */ 
#define AMPEXQON	63	/* to ampex server to enable request queue */
#define AMPEXQOFF	64	/* to ampex server to disable request queue */
#define REQDPUP 	65	/* DSDS delete pend change to dsds_svc */
#define AMPEXOPON	66	/* to ampex server to indicate operator on duty */ 
#define AMPEXOPOFF	67	/* to ampex server to indicate operator off duty */ 
#define TAPERETNOW	68	/* to ampex server to start a retrieve - NO Wait 
				   if cartridge not in unit and no operator */
#define TAPEARCHNOW	69	/* to ampex server to start an archive - NO Wait 
				   if cartridge not in unit and no operator */
#define REQREQNOW	70	/* DSDS data req w/retrieval unless on-shelf
				   and no operator */
#define AMPEXOP		71	/* to ampex server to ask if operator on duty */ 
#define TAPEHAVE	72	/* to ampex server to ask if tape is nearline */ 
#define REQREQALL	73	/* DSDS Data_Request for all ds in keylist*/
#define REQREQALLNR	74	/* as above but w/no Ampex retrieve */
#define TAPERETMORE	75	/* to ampex server to start another retrieve -
				   equivalent to TAPEDONE followed by TAPERET
				   but without interruption by another req */
#define AMPEXEQ		76	/* to ampex server to require that archives and 
				   retrieves have equal priority */
#define AMPEXARCH	77	/* to ampex server to require that archives
				   have priority over retrieves */
#define AMPEXRET	78	/* to ampex server to require that retrieves
				   have priority over archives */
#define TAPECANCEL      79      /* to ampex server to cancel all requests
                                   for the given dsds_uid */
#define REQTAPECANCEL	80	/* to dsds_svc to do ampex TAPECANCEL */
#define AMPEXCANCEL     81      /* to ampex server to cancel all queued requests
                                   with the given req_num */
#define REQDFFLG	82      /* to dsds_svc w/value to set dfflg */
#define REQREQMV	83      /* to dsds_svc w/datasets to mv to PDS set */
#define REQMVDPUP	84      /* to dsds_svc to update DP for moved ds */

/* column locations in pe log file /tmp/pe.pid.log */
#define C_PE_TID 0
#define C_PE_STATUS 1
#define C_PE_COMPLET 2
#define C_PE_ARCHIVE 3
#define C_PE_MODULE 4
#define C_PE_OUT 5
#define C_PE_AUXINFO 6

/* file with primary and secondary Ampex host names */
/*#define PHNAME "/home/jim/STAGING/tables/prod_host.cfg"*/
#define PHNAME "/home/soi/CM/tables/prod_host.cfg"


#include "module.h"
#include "pe_hdata.h"

struct serv_result {
	int serv_errno;
        int tid;
        char *host;
	char *out_file;
};
typedef struct serv_result serv_result;

/* Server definition table. One for each server type. */
struct server {
  char *name;			/* name of the remote server program */
  char version[16];		/* version # string such as V0R8B2 */
  int msgid;			/* msg # assigned to this server's messages */
  int busyall;			/* true if at least one host is busy */
  int firstserver;		/* true if first server in the pipe */
  int dsin;			/* =1 means to request input from dsds */
  int split;                    /* =1 to split the fsn-lsn among the hosts */
  int groupid;			/* non-0 and same for all servers of a group */
  int cphist;			/* cp hist log and map file to output wds */
  int noabort;                  /* don't abort map file if mod sets abortflg*/
  char archive;			/* t=temp, p=perm, n=normal, NULL=no archive */
  int archive_day;		/* days form now to archive (t) or del (n) */
  int archive_group;		/* 0=no grp, 1=in grp, -1=end of grp */
  int archive_complete;		/* 1=archive occured, no crash or abortflg */
  argument arguments[MAX_ARGS]; /* argument list sent from the server */
  KEY *map_list;		/* keylist from the map file for this serv */
  HDATA *hosts;			/* host dependent info */
};
typedef struct server PSERVER; /* originally called SERVER before Rick used it*/

/* prototypes */
//extern PADATA *getpadata(PADATA *list, char *wd, long uid);
//extern void setpadata(PADATA **list, char *wd, long uid, double bytes, 
//int stat, int substat, char *eff_date, 
//int group_id, int safe_id, unsigned long ds_index);
//extern void uidpadata(PADATA *new, PADATA **start, PADATA **end);
//extern void remuidpadata(PADATA **start, PADATA **end, char *wd, long uid);
//extern PADATA *getpawd(PADATA *list, char *wd);
//extern PADATA *getpauid(PADATA *list, long uid);
//extern PADATA *getpanext(PADATA *list);
//extern void rempadata(PADATA **list, char *wd, long uid);
//extern void updpadata(PADATA **list, char *wd, long uid, char *eff_date);
//extern void setpeuid(PEUID **list, long uid, int petid);
//extern PEUID *getpeuid(PEUID *list, long uid);
//extern PEUID *getpeuidnext(PEUID *list);
//extern void rempeuid(PEUID **list, long uid);
extern void pack_args(argument *args);
extern KEY *pack_keylist(KEY *list);
extern KEY *unpack_keylist(KEY *list);
extern HDATA *gethname (HDATA *list, char *host_name);
extern HDATA *gethtid (HDATA *list, int tid);
extern void sethname (HDATA **list, char *host_name);
extern HDATA *sethtid (HDATA *list, char *host_name, int tid);
extern HDATA *gethnext (HDATA *list);
extern HDATA *gethrel (HDATA *oldlist, HDATA *hloc, HDATA *newlist);
extern char *get_effdate(int plusdays);
extern void logkey (KEY *key);
extern void printkey (KEY *key);
extern int start_dsds(char *datab, int *dsdstid, int mytid, 
	int (*msg)(char *fmt, ...));
extern KEY *call_dsds(KEY **list, int reqtype, int dsdstid, int mytid, 
	int (*msgf)(char *fmt, ...), int dbflg);
extern KEY *unpack_dsds(int (*msgf)(char *fmt, ...));
extern KEY *call_sql(char *scmd, int dsdstid, int mytid, 
	int (*msg)(char *fmt, ...));
extern long reg_dsds(int dsdstid, int mytid, int (*msg)(char *fmt, ...));
extern int dereg_dsds(long uid, int dsdstid, int mytid, 
	int (*msg)(char *fmt, ...));
extern int spawn_dsds_local(char *db, int mytid, int *dsdstid, int debugflg,
	int (*msgf)(char *fmt, ...));
extern int spawn_pe_rpc(char *db, int mytid, int *perpctid, int debugflg,
	int (*msgf)(char *fmt, ...));
extern int start_ampex(char *db, int mytid, int perpctid, int *ampextid,
	int debugflg, int localflg, int (*msgf)(char *fmt, ...));
extern int send_ampex_request(KEY *list, int req_num, int mytid, int perpctid,
        int (*msgf)(char *fmt, ...));
extern int send_ampex_local_request(KEY *list, int req_num, int ampextid,
        int mytid, int (*msgf)(char *fmt, ...));
extern void free_key_link(KEY_LINK *link);
extern int set_key_link(KEY_LINK **link, KEY *list, int status);
extern void rm_key_link(KEY_LINK **alink, KEY_LINK *thisone);
extern int pds_wd_2_host(char *wd, char *hostname);
extern int pds_host_2_num(char *hostname);
extern int pds_wd_2_set(char *wd);
extern char *prod_host_prime();
extern char *prod_host_second();

#ifdef PVM33	/* for pvm version 3.3 and above */
#define hostinfo pvmhostinfo
#define taskinfo pvmtaskinfo
#define pvm_serror(x) pvm_setopt(PvmAutoErr,x)
#endif		/* PVM33 */

#define PE_INCL
#endif

/*
$Id: pe.h,v 1.2 2009/04/24 21:51:59 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/pe.h,v $
$Author: production $
*/
/* $Log: pe.h,v $
 * Revision 1.2  2009/04/24 21:51:59  production
 * change SERVER to PSERVER
 *
/* Revision 1.1  2009/02/23 22:44:22  production
/* initial
/*
 * Revision 1.70  2008/08/20 19:53:18  jim
 * *** empty log message ***
 *
 * Revision 1.69  2007/01/04  23:46:04  jim
 * add SUMSERVER
 *
 * Revision 1.68  2002/07/02  20:50:39  jim
 * call_dsds() now uses KEY **
 *
 * Revision 1.67  2001/08/24  17:28:49  jim
 * add PHNAME and prod_host_prime/second
 *
 * Revision 1.66  2000/04/04  19:22:41  jim
 * add pds_wd_2_set()
 *
 * Revision 1.65  2000/03/15  00:31:25  jim
 * add REQREQMV
 *
 * Revision 1.64  1999/10/25  18:41:23  jim
 *  MAX_ARG_STR 256
 *
 * Revision 1.63  1999/09/07  18:50:58  jim
 * add prototype for pds_wd_2_host & pds_host_2_num
 *
 * Revision 1.62  1999/04/28  17:18:38  jim
 * add REQDFFLG
 *
 * Revision 1.61  1998/09/09  22:41:06  kay
 * added AMPEXCANCEL
 *
 * Revision 1.60  1998/09/09  22:33:58  jim
 * add REQTAPECANCEL
 *
 * Revision 1.59  1998/09/02  17:16:30  kay
 * added TAPECANCEL
 *
 * Revision 1.58  1998/06/16  22:12:46  kehcheng
 * added PVM33 macros for pvm versions 3.3 and above
 *
 * Revision 1.57  1998/06/01 22:09:40  kehcheng
 * changed errno to serv_errno in struct serv_result
 *
 * Revision 1.56  1997/05/09 23:24:12  jim
 * add MSGUNREGPERPC
 *
 * Revision 1.55  1997/04/16  21:48:10  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.54  1997/03/05  01:30:08  kay
 * added AMPEXEQ, AMPEXARCH, AMPEXRET to control archive and retrieve priorities
 *
 * Revision 1.53  1997/02/10  17:20:17  kay
 * added TAPERETMORE, a way to allow uninterrupted sequential retrieves
 *
 * Revision 1.52  1997/01/14  19:24:44  jim
 * add prototypes
 *
 * Revision 1.51  1996/11/19  01:00:31  jim
 * add REQREQALL
 *
 * Revision 1.50  1996/09/23  17:10:02  kay
 * added ds_index to setpadata arguments
 *
 * Revision 1.49  1996/09/10  01:10:00  kay
 * added AMPEXOP and TAPEHAVE
 *
 * Revision 1.48  1996/08/30  17:51:00  jim
 * add REQREQNOW
 *
 * Revision 1.47  1996/08/20  21:27:53  kay
 * added request types AMPEXOPON, AMPEXOPOFF, TAPERETNOW, TAPEARCHNOW
 *
 * Revision 1.46  1996/08/06  21:34:22  jim
 * add void remuidpadata(PADATA **start, PADATA **end, char *wd, long uid)
 *
 * Revision 1.45  1996/07/18  21:09:51  kay
 * added stat and safe_id to arguments of setpadata
 *
 * Revision 1.44  1996/07/09  17:56:22  kay
 * added REQDPUP
 *
 * Revision 1.43  1996/06/21  21:51:59  jim
 * add PEMAILFILE
 *
 * Revision 1.42  1996/06/13  21:01:42  jim
 * change history prototype
 *
 * Revision 1.41  1996/05/08  00:50:46  jim
 * changed PELOGFILE
 *
 * Revision 1.40  1996/05/06  21:41:08  jim
 * up MAX_SERV from 64 to 128
 *
 * Revision 1.39  1996/03/29  00:39:54  kay
 * added request types AMPEXQON and AMPEXQOFF to control ampex_svc request queue
 *
 * Revision 1.38  1996/03/13  16:25:38  kay
 * added request types TAPERESERVE and TAPERELEASE
 *
 * Revision 1.37  1996/03/11  21:41:05  jim
 * up MAX_SERV
 *
 * Revision 1.36  1996/03/04  19:58:17  jim
 * add noabort to SERVER
 *
 * Revision 1.35  1996/02/29  22:57:38  jim
 * add cphist field
 *
 * Revision 1.34  1996/01/24  01:10:30  jim
 * add REQALLOCUP
 *
 * Revision 1.33  1995/12/28  17:51:32  kay
 * added TAPEIN and TAPEINDONE request types
 *
 * Revision 1.32  1995/12/08  00:01:55  jim
 * up MAX_SERV to 48 from 32
 *
 * Revision 1.31  1995/12/02  00:20:45  jim
 * add the pe log AT defs
 *
 * Revision 1.30  1995/12/01  19:52:55  jim
 * add PELOGFILE
 *
 * Revision 1.29  1995/11/15  17:07:36  jim
 * added REQAPUP and some prototypes
 *
 * Revision 1.28  1995/08/16  23:25:24  jim
 * added archive_complete to SERVER struct
 *
 * Revision 1.27  1995/06/19  15:56:50  jayasree
 * added group_id to tape structure.
 *
 * Revision 1.26  1995/06/05  18:04:42  jim
 * add REQUPDATE0
 *
 * Revision 1.25  1995/04/24  18:05:28  kay
 * added AMPEX message and request numbers
 *
 * Revision 1.24  1995/04/20  23:01:12  jim
 * changed setpadata prototype for double bytes
 *
 * Revision 1.23  1995/03/27  22:30:56  jim
 * fixed bad comment on last checkin
 *
 * Revision 1.22  1995/03/27  22:28:06  jim
 * added to SERVER  int archive_group;  0=no grp, 1=in grp, -1=end of grp
 *
 * Revision 1.21  1995/03/07  21:38:05  jim
 * added REQ5KLOGU
 *
 * Revision 1.20  1995/01/19  22:07:27  jim
 * add LAGOKILL def
 *
 * Revision 1.19  1995/01/19  18:03:47  jim
 * add unpack_dsds() prototype
 *
 * Revision 1.18  1995/01/17  18:59:36  jim
 * changed call_dsds() prototype for int (*msg)
 *
 * Revision 1.17  1995/01/10  18:30:57  jim
 * new prototypes
 *
 * Revision 1.16  1994/11/11  17:31:41  kay
 * eliminated duplicate message number 34
 *
 * Revision 1.15  1994/11/11  17:19:50  jim
 * add REQSQL
 *
 * Revision 1.14  1994/10/17  22:00:01  jim
 * no change
 *
 * Revision 1.13  1994/10/17  19:05:09  jim
 * changed MAX_SERV to 24
 *
 * Revision 1.12  1994/09/16  16:19:23  kay
 * restored previous numbers for DSDS requests
 *
 * Revision 1.11  1994/09/16  00:26:01  kay
 * added prototypes for logkey and printkey and
 * added some lago message types
 *
 * Revision 1.10  1994/09/02  17:55:44  jim
 * add prototype for uidpadata
 *
 * Revision 1.9  1994/09/02  17:48:03  jim
 * add groupid, archive, archive_day to SERVER. Elim elbow and oneserver.
 *
 * Revision 1.8  1994/07/19  21:07:37  jim
 * added REQILOGU and REQMLOGU
 *
 * Revision 1.7  1994/06/07  15:36:58  CM
 * CM changed MAX_ARGS to 64
 *
 * Revision 1.6  1994/06/01  18:53:52  jim
 * made version in server table a string. added reqreqnr.
 *
 * Revision 1.5  1994/04/15  18:50:46  jim
 * add REQREQNR
 *
 * Revision 1.4  1994/04/13  22:08:47  jim
 * added MSGKILL
 *
 * Revision 1.3  1994/03/08  21:07:28  kay
 * renumbered request codes to be unique
 *
 * Revision 1.2  1994/02/28  20:37:41  kay
 * added message and request types to support lui (Lago User Interface)
 *
 * Revision 1.1  1994/02/16  23:21:14  CM
 * Initial revision
 * */
