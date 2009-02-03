/*
 *  parse_arg.c							CM/src/libast.d
 *
 *  Function(s) for parsing argument values supplied to strategy modules
 *
 *  Contents:
 *    int parse_array (KEY **params, char *root, int type)
 *    int parse_numerated (char *klist, char *init_val, int **vals,
 *	  char ***names)
 *
 *  Responsible: Rick Bogart			RBogart@solar.Stanford.EDU
 *		Kay Leibrand			KLeibrand@solar.Stanford.EDU    
 *
 *  Notes:
 *    parse array() accepts strings of the form
 *	[DL] [token0 [, token1 [, ...]]] [DR]
 *	with optional whitespace and left and right delimiters (), [], or {},
 *	and fills the params keylist with entries describing an array of values
 *	of the given type based on the parsed tokens.
 *    Strings that cannot be parsed according to the given type are ignored;
 *	it is thus possible to supply an empty array, in which case root_nvals
 *	is set to 0 and no values of root_i_value are set.  It is up to the
 *	module to deal with this situation.  The argument value
 *	  [0, foo, 1] would produce a two-element array of integers, while the
 *	value
 *	  {xxx} would produce an empty array (of numerics).
 *    Integer tokens are parsed according to standard rules for strtol (, , 0),
 *	so that octal, decimal, and hexadecimal representations may be used
 *	and mixed freely
 *
 *  Bugs:
 *    parse_array() is only implemented for arrays of ints & floats (doubles);
 *	if an unknown keytype is presented it sets soi_errno and sets the
 *	number of values in the array to 0
 *    Although the left and right delimiters are optional, if a left delimiter
 *	is present there must be a matching right delimiter of the same type
 *    parse_numerated() is only a stub function so far; it is designed to
 *	support arguments of type ARG_NUME, with enumerated sets of permissible
 *	values
 *
 *  Revision history is at end of file.
 */

#include <soi_error.h>
#include "soi_names.h"
#include <stdlib.h>
#include <string.h>

static char *strip_delimiters (char *s) {
  char *match;
  switch (*s) {
    case '{': {
      if (match = strchr (s, '}')) *match = 0;
      else return NULL;
      return s+1;
    }
    case '[': {
      if (match = strchr (s, ']')) *match = 0;
      else return NULL;
      return s+1;
    }
    case '(': {
      if (match = strchr (s, ')')) *match = 0;
      else return NULL;
      return s+1;
    }
    default:
      return s;
  }
}

int parse_array (KEY **params, char *root, int type) {
/*
 *  Finds an entry in the params list with the key root and parses it as
 *    an array of values of given type.  Entries for each of the parsed
 *    values are inserted into the params list.  Keys for these entries
 *    are of the form root'_'n'_value', and an additional entry with
 *    key_name root'_nvals' is added.
 */
  double dval;
  long ival;
  int nvals = 0;
  char *endptr;
  char *name = getkey_str (*params, root);
  char *next, *nptr;
  static char key[MAX_STRLEN];

  soi_errno = NO_ERROR; 
  if (!name) return soi_errno = MISSING_DATA_NAME;
  str_compress (name);				     /*  remove white space  */
  name = strip_delimiters (name);			  /*  (), [], or {}  */
  if (!name) return soi_errno = MISSING_DATA_NAME;
/*
 *  name should contain a comma separated list of entities of the given type
 */
  next = name;
  while (next) {
    nptr = next;
    if (next = (char *)strchr (nptr, ',')) {
      *next = 0;
      next++;
    }
    if (!strlen (nptr)) continue;
					      /*  nptr now points to entity  */
    sprintf (key, "%s_%d_value", root, nvals);
    switch (type) {
      case KEYTYP_INT:
      case KEYTYP_SHORT:
      case KEYTYP_USHORT:
      case KEYTYP_BYTE:
      case KEYTYP_UBYTE:
        ival = strtol (nptr, &endptr, 0);
	if (endptr == nptr) continue;
        setkey_int (params, key, ival);
        break;
      case KEYTYP_DOUBLE:
      case KEYTYP_FLOAT:
        dval = strtod (nptr, &endptr);
	if (endptr == nptr) continue;
        setkey_double (params, key, dval);
        break;
						   /**  ... more cases ...  **/
      default:
        soi_errno = KEYTYPE_NOT_IMPLEMENTED;
    }
						      /**  check soi_errno  **/
    nvals++;
  }
  sprintf (key, "%s_nvals", root);
  setkey_int (params, key, nvals);
  return soi_errno;
}

int parse_numerated (char *klist, char ***names) {
/*
 *  Parses an entry in the args list of type ARG_NUME and returns the number
 *    of possible values, determined from the range entry interpreted as a
 *    comma-spearated set of strings, and a mallocd array of strings
 *    corresponding to the enumeration choices
 */
  int found = 0, maxlen;
  char *next, *nptr, *tmp;

  if (!klist) return found;
  maxlen = strlen (klist);
  tmp = malloc (maxlen + 1);
  strcpy (tmp, klist);
  str_compress (tmp);				     /*  remove white space  */
  next = tmp;

  *names = (char **)malloc (maxlen * sizeof (char *));
  while (next) {
    nptr = next;
    if (next = (char *)strchr (nptr, ',')) {
      *next = 0;
      next++;
    }
    if (!strlen (nptr)) continue;
					      /*  nptr now points to entity  */
    (*names)[found] = (char *)malloc (strlen (nptr) + 1);
    strcpy ((*names)[found], nptr);
    found++;
  }
  free (tmp);

  return found;
}
