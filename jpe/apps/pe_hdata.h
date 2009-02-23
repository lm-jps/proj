/*
 * pe_hdata.h
 *
 * HDATA is a linked list storage structure for each pipeline server that 
 * is running on a set of hosts.
 * Note: A host can have multiple servers of the same type, so host names are
 *       not necessarily unique. However, tid's are guaranteed unique.
 *
 * The following functions are available for use with hdata lists:
 *	HDATA *gethname (HDATA *list, char *host_name)
 *		return pointer to entry containing host_name
 *	HDATA *gethtid (HDATA *list, int tid)
 *		return pointer to entry containing tid
 *	HDATA *sethtid (HDATA *list, char *host_name, int tid)
 *		adds the tid for the entry with host_name and tid=0
 *	void sethname (HDATA **list, char *host_name)
 *		adds an enry, unconditionally, for the host_name and sets
 *		tid and busy = 0 and param_list = null
 *      HDATA *gethnext (HDATA *list)
 *              get first entry on list, if list = null gets next entry
 *              Note: Do not interleave calls for two different lists
 *      HDATA *gethrel (HDATA *oldlist, HDATA *hloc, HDATA *newlist)
 *              get the same relative entry on the new list as the one given
 *              on the old list
*/
#ifndef PE_HDATA_INCL

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

#define PE_HDATA_VERSION_NUM    (0.8)

/* Host data for a server. Linked list for all hosts for a server type */
struct hdata {
  struct hdata *next;
  char *host_name;		/* host name that the server is running on */
  int tid;			/* task id on the host */
  int busy;			/* busy flag on the host */
  KEY *param_list;		/* param list sent to the server on the host */
};
typedef struct hdata HDATA;

#define PE_HDATA_INCL
#endif
/*
 * Revision History
 * v0.1		93.05.26	ja
 *	created original file
 *
*/

/*
$Id: pe_hdata.h,v 1.1 2009/02/23 22:45:40 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/pe_hdata.h,v $
$Author: production $
*/
/* $Log: pe_hdata.h,v $
 * Revision 1.1  2009/02/23 22:45:40  production
 * initial
 *
 * Revision 1.2  1997/04/16  21:48:57  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.1  1994/02/16  23:21:14  CM
 * Initial revision
 * */
