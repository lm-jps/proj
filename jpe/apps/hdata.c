/* hdata.c
 *
 * HDATA is a linked list storage structure for each pipeline server that
 * is running on a set of host.
 * Note: A host can have multiple servers of the same type, so host names are
 *       not necessarily unique. However, tid's are guaranteed non-0 & unique.
 *
 * The following functions are available for use with hdata lists:
 *      void sethname (HDATA **list, char *host_name)
 *              adds an entry, unconditionally, for the host_name and sets
 *		tid and busy = 0 and param_list = null
 *      HDATA *sethtid (HDATA *list, char *host_name, int tid)
 *              adds the tid for the entry with host_name and tid=0
 *      HDATA *gethname (HDATA *list, char *host_name)
 *              return pointer to entry containing host_name
 *      HDATA *gethtid (HDATA *list, int tid)
 *              return pointer to entry containing tid
 *	HDATA *gethnext (HDATA *list)
 *		get first entry on list, if list = -1 gets next entry
 *		Note: Do not interleave calls for two different lists
 *	HDATA *gethrel (HDATA *oldlist, HDATA *hloc, HDATA *newlist)
 *		get the same relative entry on the new list as the one given
 *		on the old list
*/

#include <stdlib.h>
#include <stdio.h>
#include <soi_key.h>
#include "pe_hdata.h"

#define MONE (long)-1	/* needs to be long for IRIX64 compile */

int kludgey = 2;	/* workaround for malloc / alignment problem */
/*int kludgey = 1;	/* workaround for malloc / alignment problem */

HDATA *hnext; 		/* used between gethnext() calls */

HDATA *gethname (HDATA *list, char *host_name)
{
  HDATA *walk = list;

  if(!host_name) return NULL;
  while(walk) {
    if(strcmp(walk->host_name, host_name))
      walk = walk->next;
    else
      return walk;
  }
  return walk;
}

HDATA *gethtid (HDATA *list, int tid)
{
  HDATA *walk = list;

  while(walk) {
    if(walk->tid == tid)
      return walk;
    else
      walk = walk->next;
  }
  return walk;
}

char *h_strdup (char *s)
{
  char *duped;

  if(!s) return NULL;
  duped = (char *)malloc(kludgey*(strlen(s)+1));
  strcpy(duped, s);
  return duped;
}

void sethname (HDATA **list, char *host_name)
{
  HDATA *newone;

  if(!host_name) return;
  newone = (HDATA *)malloc(kludgey*sizeof(HDATA));
  newone->next = *list;
  newone->host_name = h_strdup(host_name);
  newone->tid = 0;
  newone->busy = 0;
  newone->param_list = NULL;
  *list = newone;
}

HDATA *sethtid (HDATA *list, char *host_name, int tid)
{
  HDATA *walk = list;

  if(!host_name) return NULL;
  while(walk) {
    if(strcmp(walk->host_name, host_name))
      walk = walk->next;
    else {
      if(walk->tid == 0) {
        walk->tid = tid;
        return walk;
      }
      else
        walk = walk->next;
    }
  }
  return walk;
}

/* get first entry on list, if list = -1 gets next entry 
 * Note: Do not interleave calls for two different lists.
*/
HDATA *gethnext (HDATA *list)
{
  if(list != (HDATA *)MONE) {
    hnext = list;
    return list;
  }
  if(hnext)
    hnext = hnext->next;
  return hnext;
}

/* get the same relative entry on the new list as given on the old list */
HDATA *gethrel (HDATA *oldlist, HDATA *hloc, HDATA *newlist)
{
  HDATA *hstep, *hstep2;

  if(!oldlist || !hloc || !newlist)
    return NULL;
  hstep2=newlist;
  for(hstep=(HDATA *)gethnext(oldlist); hstep != NULL; 
	hstep=(HDATA *)gethnext((HDATA *)MONE)) {
    if(hstep == hloc) break;
    hstep2=hstep2->next;
  }
  if(hstep)
    return hstep2;
  else
    return NULL;
}
