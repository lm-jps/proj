//#include "module.h"
#include "pe.h"
#include <soi_key.h>

/* Packs the fields of an argument data type into the current pvm buffer.
*/
void pack_args(argument *args)
{
  pvm_pkint(&args->kind,1,1);
  pvm_pkstr(args->key);
  pvm_pkstr(args->default_value);
  pvm_pkstr(args->range);
  pvm_pkstr(args->description);
}

/* Packs a linked list keylist into a pvm send buffer and adds an entry with
the name equal END_NAME at the end of the send buffer.
Called with the pointer of the keylist to pack.
Returns NULL if successful, else returns the KEY * of the failing entry.
*/
KEY *pack_keylist(KEY *list)
{
  KEY *walker = list;
  int info;

  while(walker) {
    pvm_pkstr(walker->name);
    pvm_pkint(&walker->type, 1, 1);
    switch(walker->type) {		/* pack value according to type */
    case KEYTYP_STRING:
      info=pvm_pkstr((char *)walker->val);
      break;
    case KEYTYP_BYTE:
      info=pvm_pkbyte((char *)walker->val, 1, 1);
      break;
    case KEYTYP_INT:
      info=pvm_pkint((int *)walker->val, 1, 1);
      break;
    case KEYTYP_FLOAT:
      info=pvm_pkfloat((float *)walker->val, 1, 1);
      break;
    case KEYTYP_DOUBLE:
      info=pvm_pkdouble((double *)walker->val, 1, 1);
      break;
    case KEYTYP_TIME:
      info=pvm_pkdouble((double *)walker->val, 1, 1);
      break;
    case KEYTYP_SHORT:
      info=pvm_pkshort((short *)walker->val, 1, 1);
      break;
    case KEYTYP_LONG:
      info=pvm_pklong((long *)walker->val, 1, 1);
      break;
    case KEYTYP_UBYTE:
      info=pvm_pkbyte((unsigned char *)walker->val, 1, 1);
      break;
    case KEYTYP_USHORT:
      info=pvm_pkshort((unsigned short *)walker->val, 1, 1);
      break;
    case KEYTYP_UINT:
      info=pvm_pkuint((unsigned int *)walker->val, 1, 1);
      break;
    case KEYTYP_UINT32:
      info=pvm_pkuint((unsigned int *)walker->val, 1, 1);
      break;
    case KEYTYP_ULONG:
      info=pvm_pkulong((unsigned long *)walker->val, 1, 1);
      break;
    }
    if(info)
      return(walker);
    walker=walker->next;
  }
  pvm_pkstr("END_NAME");
  return(walker);
}

/* Unpacks the pvm receive buffer into a linked list keylist.
Called with the pointer of the keylist to unpack into.
Returns the new keylist pointer if successful, else returns NULL.
*/

KEY *unpack_keylist(KEY *list)
{
  int info, type;
  char name[4096], val[131072];

  while(1) {
    pvm_upkstr(name);
    if(!strcmp(name, "END_NAME"))
      break;
    pvm_upkint(&type, 1, 1);
    switch(type) {		/* unpack value according to type */
    case KEYTYP_STRING:
      info=pvm_upkstr((char *)val);
      break;
    case KEYTYP_BYTE:
      info=pvm_upkbyte((char *)val, 1, 1);
      break;
    case KEYTYP_INT:
      info=pvm_upkint((int *)val, 1, 1);
      break;
    case KEYTYP_FLOAT:
      info=pvm_upkfloat((float *)val, 1, 1);
      break;
    case KEYTYP_DOUBLE:
      info=pvm_upkdouble((double *)val, 1, 1);
      break;
    case KEYTYP_TIME:
      info=pvm_upkdouble((double *)val, 1, 1);
      break;
    case KEYTYP_SHORT:
      info=pvm_upkshort((short *)val, 1, 1);
      break;
    case KEYTYP_LONG:
      info=pvm_upklong((long *)val, 1, 1);
      break;
    case KEYTYP_UBYTE:
      info=pvm_upkbyte((unsigned char *)val, 1, 1);
      break;
    case KEYTYP_USHORT:
      info=pvm_upkshort((unsigned short *)val, 1, 1);
      break;
    case KEYTYP_UINT:
      info=pvm_upkuint((unsigned int *)val, 1, 1);
      break;
    case KEYTYP_UINT32:
      info=pvm_upkuint((unsigned int *)val, 1, 1);
      break;
    case KEYTYP_ULONG:
      info=pvm_upkulong((unsigned long *)val, 1, 1);
      break;
    }
    if(info)
      return(NULL);
    /*setkey(&list, name, val, type);	/* add this entry to the keylist */
    addkey(&list, name, val, type);	/* add entry w/o cking for dups */
  }
  return(list);
}
