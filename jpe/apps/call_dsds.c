static char rcsid[] = "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/call_dsds.c,v 1.1 2009/02/03 22:17:57 production Exp $";
/* Routines to do a call to dsds_svc and get a respose.
*/
#include <sys/time.h>
#include <pvm3.h>
#include <soi_error.h>
#include <unistd.h>
#include "pe.h"

#define RESPDSDS 3600		/* sec to wait for dsds_svc response */

void printkey (KEY *key);
void print_asctime(int (*msgf)(char *fmt, ...));
KEY *unpack_dsds(int (*msgf)(char *fmt, ...));
int resp_dsds(int dsdstid);

static int reqack;

/* Send the given keylist to dsds_svc for the given request type and await its 
 * response. Returns the keylist returned by dsds_svc as a new keylist, or 
 * NULL for error.
 * Args:
 *   list = addr of keylist to send to dsds_svc. mytid gets added to it
 *   reqtype = the dsds_svc request type, e.g. open, allocate, etc.
 *   dsdstid = the pvm task id of dsds_svc
 *   mytid = the pvm task id of the caller, so dsds_svc knows who to answer
 *   msgf = function to print messages with as we go
 *   dbflg = if non-0 then will print the keylist going in and coming out, 
 *           this is a special use flag, normally it will be 0.
*/
KEY *call_dsds(KEY **list, int reqtype, int dsdstid, int mytid, 
		int (*msgf)(char *fmt, ...), int dbflg)
{
  KEY *keybad, *retlist;
  char *username;
  int status;

  if(!dsdstid) {
    (*msgf)("Called call_dsds() w/null dsdstid\n");
    return(NULL);
  }
  pvm_initsend(PvmDataDefault); /* init the send buffer */
  pvm_pkint(&reqtype, 1, 1);	/* pack the reqtype flag */
  status=0;			/* status not applicable here */
  pvm_pkint(&status, 1, 1);
  setkey_int(list, "pe_tid", mytid);	/* always supply our tid */
  if(!(username = (char *)getenv("USER"))) username = "nouser";
  setkey_str(list, "USER", username); /* always supply our user */
  if(dbflg) {			/* don't use with tae ui */
    printf("\nThe keylist before the dsds_svc call %d is:\n", reqtype);
    keyiterate( printkey, *list);
  }
  if(keybad=(KEY *)pack_keylist(*list)) {
    (*msgf)("Err packing a pvm msg to dsds_svc, type=%d name=%s\n",
                keybad->type, keybad->name);
    return(NULL);
  }
  if(pvm_send(dsdstid, MSGDSDS)) {
    (*msgf)("Error calling dsds_svc\n");
    return(NULL);
  }
  if(resp_dsds(dsdstid)) {
    (*msgf)("Timeout receiving response from dsds_svc (req=%d): ", reqtype);
    print_asctime(msgf);
    return(NULL);
  }
  if(retlist = unpack_dsds(msgf)) {
    if(dbflg) {			/* don't use with tae ui */
      printf("\nThe keylist after the dsds_svc call %d is:\n", reqack);
      keyiterate( printkey, retlist);
    }
  }
  return(retlist);
}

/* Unpack the message just received from dsds. Returns the unpacked keylist
 * as a new keylist or NULL for error.
*/
KEY *unpack_dsds(int (*msgf)(char *fmt, ...))
{
  KEY *list;
  int status;
  char *errtxt;

  pvm_upkint(&reqack, 1, 1);			/* request type N/A here */
  pvm_upkint(&status, 1, 1);			/* get dsds_svc status */
  if(status) {
    switch(status) {			/* will be replaced by err str rte */
    case DS_NO_LAGO:
      errtxt = "Can't retrieve. No ampex_svc";
      break;
    case DS_ALLOC_ERR:
      errtxt = "Can't allocate storage";
      break;
    case DS_DATA_QRY:
      errtxt = "DB query failed. No such ds?";
      break;
    case DS_MAX_LAGO_ERR:
      errtxt = "Max# of active Ampex retrieves";
      break;
    case AMPEX_READ_ERROR:
      errtxt = "Ampex read error";
      break;
    default:
      errtxt = "see soi_error.h";
      break;
    }
    (*msgf)("Dsds_svc returned error code %d: %s\n", status, errtxt);
    return(NULL);
  }
  list=newkeylist();				/* create empty keylist */
  if(!(list=(KEY *)unpack_keylist(list))) {	/* unpack msg into keylist */
    (*msgf)("Error unpacking keylist from dsds_svc\n");
    return(NULL);
  }
  return(list);
}

/* Wait for dsds_svc to give a response. Returns 0 on completion, else
 * returns -1 if response timeout */
int resp_dsds(int dsdstid)
{
  struct timeval tvalr;
  long tsr;
  int retcode, bufid;

  /* Timeout in case never get a response */
  gettimeofday(&tvalr, NULL);
  tsr = tvalr.tv_sec;
  while(1) {
    if(bufid=pvm_nrecv(dsdstid, MSGDSDS)) {
      if(bufid < 0) 
        retcode = -1;
      else
        retcode = 0;
      break;				/* exit while loop */
    }					/* end if get any completion */
    gettimeofday(&tvalr, NULL);
    if(tvalr.tv_sec - tsr > RESPDSDS) {
      retcode = -1;			/* give failure ret code */
      break;;				/* exit while loop */
    }
    usleep(200);			/* give the cpu a break */
  }					/* end while(1) */
  return(retcode);
}

/* print out the current date/time in the asctime() format
*/
void print_asctime(int (*msgf)(char *fmt, ...))
{
  struct timeval tvalr;
  struct tm *t_ptr;
  int tvalr_int;

  gettimeofday(&tvalr, NULL);
  tvalr_int = (int)tvalr.tv_sec; /* need int vrbl for this to work on sgi4*/
  t_ptr = localtime(&tvalr_int);
  (*msgf)("%s", asctime(t_ptr));
}

/*
$Id: call_dsds.c,v 1.1 2009/02/03 22:17:57 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/call_dsds.c,v $
$Author: production $
*/
/* $Log: call_dsds.c,v $
 * Revision 1.1  2009/02/03 22:17:57  production
 * initial
 *
 * Revision 1.22  2004/11/05  19:21:30  jim
 * change usleep() to 200
 *
 * Revision 1.21  2004/06/22 16:46:28  jim
 * add usleep()
 *
 * Revision 1.20  2004/06/17 17:56:36  jim
 * elim sleep(1)
 *
 * Revision 1.19  2002/07/02  20:28:41  jim
 * arg now KEY **
 *
 * Revision 1.18  1998/11/19  17:45:14  jim
 * add AMPEX_READ_ERROR
 *
 * Revision 1.17  1998/09/15  21:39:06  jim
 * add another msg
 *
 * Revision 1.16  1997/04/10  20:19:15  jim
 * add reqtype to timeout msg
 *
 * Revision 1.15  1997/04/07  18:04:54  jim
 * add print_asctime()
 *
 * Revision 1.14  1997/01/08  01:08:55  jim
 * add more errtxt
 *
 * Revision 1.13  1996/07/11  23:41:14  jim
 * use sginap()
 *
 * Revision 1.12  1996/06/18  17:02:00  jim
 * make (char *fmt, ...)
 *
 * Revision 1.11  1996/04/15  16:23:33  jim
 * change "lago_svc" msg to "ampex_svc"
 *
 * Revision 1.10  1996/03/11  21:53:42  jim
 * increase response timeout
 *
 * Revision 1.9  1996/01/25  18:18:18  jim
 * change response t.o. to 15 minutes
 *
 * Revision 1.8  1995/11/15  17:51:13  jim
 * put USER in keylist
 *
 * Revision 1.7  1995/09/18  21:37:47  jim
 * changed dsds timeout
 *
 * Revision 1.6  1995/05/26  23:26:07  jim
 * changed response timeout to 10 minutes
 *
 * Revision 1.5  1995/01/19  18:09:42  jim
 * increase dsds timeout due to DS_SqlDo() running too long. Temp fix.
 * Also change void msg() to int msg().
 *
 * Revision 1.4  1994/09/19  18:41:14  jim
 * call unpack_dsds w/o KEY *list
 *
 * Revision 1.3  1994/06/24  23:40:59  jim
 * changed prototype
 *
 * Revision 1.2  1994/06/24  23:03:57  jim
 * original
 * */
