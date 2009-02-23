/*
 *  names.c					~soi/(version)/src/libast.d
 *
 *  functions to parse data names
 *  according to Kay Leibrand 94.09.01 revision/extension of
 *  Phil Scherrer's 93.03.31 memo on Use of DataSet names
 *
 *  All parsing routines which return keylists expect to parse strings with
 *  no white space.  They set soi_errno and return NULL when errors 
 *  are encountered.  The following syntax is assumed.
 *
 *  Notation for Formal Syntax of Data Names

 *  This is revised from TN-112.
 *  The main changes are extended notation for data collections and
 *  different notation for array slices.

 *  The notation is as in TN-112 where ::= means consists of,
 *  <> indicates concatenation, and | means or.  Any characters
 *  not on the left hand side of ::= should be included as is.  For example,
 *  series: should literally appear as the first characters of a series_selector.
 *  The curly braces {} indicate optional parts except in the definition of
 *  placeholder where they should be included as is.

 *  datacollection_name::=datasubset_name{;datacollection_name}
 *  datasuperset_name::=project_selector{,class_selector},series_collector
 *  dataset_name::=project_selector{,class_selector},series_selector
 *  datasubset_name::=dataset_name{,sel:data_selector}|
                  datasuperset_name{,sel:data_selector}|data_selector
 *  project_selector::=proj:member_selector|prog:member_selector
 *  class_selector::=class_name:member_selector{,class_selector}
 *  series_collector::=series:member_collector
 *  series_selector::=series:member_selector
 *  member_selector::=member_name{[member_number]}
 *  member_collector::=member_name[member_range]
 *  data_selector::=directory_name{file_name}|file_name|
                {var_name}record_range{array_slice}{,fmt:format}
 *  member_range::=range{,member_range}
 *  range::=number{-number}
 *  array_slice::=[slice_range]{array_slice}
 *  record_range::=[{slice_range}]
 *  slice_range::=range{,increment}
 *  increment::=number
 *  number::=digit|number<>digit
 *  member_name::=name
 *  member_number::=number
 *  class_name::=name
 *  directory_name::=name/{directory_name}
 *  file_name::=name
 *  var_name::=name
 *  name::=alphabetic|name<>alphanumeric
 *  format::=%<>number<>d|%<>number<>T

 *  template::=subtemplate|wd:subtemplate;bn:subtemplate
 *  subtemplate::={subtemplate<>}literal|{subtemplate<>}placeholder
 *  placeholder::={name}|{#name}|{#format#name}
 *  literal::=character|literal<>character

 *  Contents: are in soi_names.h
 *
 *  Responsible:  Kay Leibrand			KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *
 *  Planned updates:
 *    Add functions to produce other kinds of directory/file names besides
 *	FITS files (e.g. CDF).
 *
 *  Revision history is at end of file.
 */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <ctype.h>
#include "soi_key.h"
#include "soi_names.h"
//#include "timeio.h"
#include "soi_error.h"

#define VERSION_NUM	(0.8)

/* global variables used in keyiterate */
static KEY *__collection_list;
static char __tailfmt[MAX_STRLEN];
static char __keyfmt[MAX_STRLEN];

static int add_env_template (KEY **tolist, char *root) {
/*
 * tries to add the template from the environment corresponding to 
 * the proj or prog name IF there is not a template already in the list 
 * i.e. does not override an existing template 
 * failure to find a template in the environment is not considered an error
 */

   char *pro;
   static char key[MAX_STRLEN];
   static char pro_key[MAX_STRLEN];

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return soi_errno;} 

   sprintf (key, "%s_rule", root);
   if (getkeytype (*tolist, key) == KEYTYP_STRING) /* there is a template */
      return NO_ERROR;

   sprintf (pro_key, "%s_proj", root);
   if (pro = GETKEY_str (*tolist, pro_key))
      setkey_str (tolist, key, getenv (pro)); /* try to get it from env */

   if (getkeytype (*tolist, key) == KEYTYP_STRING)  /* there is a template */
      return NO_ERROR;

   sprintf (pro_key, "%s_prog", root);
   if (pro = GETKEY_str (*tolist, pro_key))
      setkey_str (tolist, key, getenv (pro)); /* try to get it from env */
   return NO_ERROR;
}

static int add_basename (KEY **tolist, char *root) {
/*
 *  if there is a template in the given list or if a template can be
 *  obtained from the environment, it is applied to the list
 *  entries to form a basename which is then added to the list.
 *  if there is no template the basename is an empty string.
 *  if an error occurs, the list is unchanged and status is returned.
 */
   char *buffer, *charptr;
   static char key[MAX_STRLEN];
   char *template = NULL;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return soi_errno;} 

   add_env_template (tolist, root); /* if not in list */

   sprintf (key, "%s_rule", root);
   if (template = getkey_str (*tolist, key)) {
      if (charptr=(char*)strstr(template,"bn:")) {
         template=(char*)strchr(charptr,':'); template++;
      }
      else template=NULL;
   }
   /* template points to bn: portion of in_rule, or is NULL */

   sprintf (key, "%s_basename", root);
   if (!template) 
      setkey_str (tolist, key, "");
   else {
      if (!(buffer=fill_template(template, *tolist, root))) return soi_errno;
      str_collapse(buffer,'/');
      setkey_str (tolist, key, buffer);
   }
   /* if(template) free(template); /* This may cause a crash?? */
   return NO_ERROR;
}

static int add_directory (KEY **tolist, char *root) {
/*
 *  if there is a template in the given list or if a template can be
 *  obtained from the environment, it is applied to the list
 *  entries to form a directory name which is then added to the list.
 *  if there is no template the directory name is an empty string.
 *  if an error occurs, the list is unchanged and status is returned.
 *  will NOT OVERRIDE and existing root_wd entry in the list
 */
   char *buffer, *charptr;
   static char key[MAX_STRLEN];
   char *template = NULL;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return soi_errno;} 

   sprintf (key, "%s_wd", root);
   if (getkeytype (*tolist, key) == KEYTYP_STRING) return NO_ERROR;

   add_env_template (tolist, root); /* if not in list */

   sprintf (key, "%s_rule", root);
   if (template = getkey_str (*tolist, key)) {
      if (charptr=(char*)strstr(template,"wd:")) {
         template=(char*)strchr(template,':'); template++;
      } 
      if (charptr=(char*)strstr(template,"bn:")) {
         *(--charptr)= 0;
      }
   }
       /*  template points to entire in_rule, or to wd: portion, or is NULL  */
   sprintf (key, "%s_wd", root);
   if (!template) 
      setkey_str (tolist, key, ".");
   else {
      if (!(buffer=fill_template(template, *tolist, root))) return soi_errno;
      str_collapse(buffer,'/');
      setkey_str (tolist, key, buffer);
   }
   /* if(template) free(template); /* This may cause a crash?? */
   return NO_ERROR;
}

static int add_paths (KEY **tolist, char *root) {
/* adds directory and basename entries to list for multiple datasets */
/*
 *  takes a pointer to a keylist and a rootkey as arguments.
 *  If there is a template in the given list or if a template can be obtained
 *  from the environment, it is applied to the list entries to form
 *  a directory name which is then added to the list. If there is no template,
 *  the directory name is an empty string. If an error occurs, the list
 *  is unchanged and status is returned. It will NOT OVERRIDE and existing
 *  root_wd entry in the list. 
 */
   static char key[MAX_STRLEN];
   static char root_n[MAX_STRLEN];
   int nsets = 1;
   int n, status;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return soi_errno;} 

   sprintf (key, "%s_nsets", root);
   if (getkeytype (*tolist, key) == KEYTYP_INT) 
      nsets = getkey_int (*tolist, key);

   if (nsets == 1){ 
      if (status = add_directory(tolist, root)) return status;
      if (status = add_basename(tolist, root)) return status;
   }

   for (n=0; n<nsets; n++) { /* for each set */
      sprintf (root_n,"%s_%d",root,n);
      if (status = add_directory(tolist, root_n)) return status;
      if (status = add_basename(tolist, root_n)) return status;
   }
   return NO_ERROR;
}

static void add_member (KEY *member) {
/*
 *  takes a pointer to a keylist member as an argument and adds this member
 *  to a global keylist called __collection_list.  This function add_member
 *  is used by parse_datacollection as an argument to the keylist function
 *  keyiterate.  It should not be used in any other context. 
 */
   static char tail[MAX_STRLEN];
   static char key[MAX_STRLEN];

   soi_errno = NO_ERROR; 
   sscanf (member->name, __tailfmt, tail);
   sprintf (key, __keyfmt, tail);
   setkey_any (&__collection_list, key, member->val, member->type);
}

static char *expand_range (char *member_range) {
/*
 *  expands a member range into a comma-separated list of numbers
 *    returns a NULL pointer if there is a parsing error
 *  WARNING - the returned pointer is a pointer to a static char array 
 *    whose contents are destroyed by the next call to this function.
 *  There is (allegedly) no check for maximum string length. 
 */
   static char expanded_range[MAX_STRLEN];
   static char local[MAX_STRLEN];
   static char lsnstr[MAX_STRLEN];
   char *charptr=local;
   char *range;
   int range_limit;
   int n, sn, fsn, lsn;

   soi_errno = NO_ERROR; 
   if (!member_range) {soi_errno = MISSING_DATA_NAME; return NULL;}
   if (strlen (member_range) >= MAX_STRLEN) {
     soi_errno=NAME_TOO_LONG; return NULL;
   }
   strcpy (local, member_range);
   expanded_range[0] = 0;			/*  start with empty string  */

   while (charptr) {
      range = charptr;
      if (charptr = (char*)strchr (range,',')) {
         *charptr = 0; charptr++;
      }
      if ((n = sscanf (range, "%d-%d", &fsn, &lsn)) == 1) lsn = fsn;
      else if (n !=2) {soi_errno = RANGE_SYNTAX_ERROR; return NULL;}

      sprintf (lsnstr, "%d,", lsn); 
      range_limit = MAX_STRLEN - strlen (lsnstr);
      for (sn = fsn; sn <= lsn; sn++) {
         if (strlen (expanded_range) >= range_limit) {
	    errstk ("error in parse_list:expand_range(): strlen (range list) > %d\n  probably too many indices\n",
	        range_limit);
            soi_errno = NAME_TOO_LONG;
	    return NULL;
         }
         sprintf (expanded_range, "%s%d,", expanded_range, sn);
      }
   }
   *((char*)strrchr (expanded_range, ',')) = 0;	  /*  remove the last comma  */
   return expanded_range;
}

static KEY *parse_class (char *collector, char *root) {
/*
 *  takes two character strings, a
 *  [project|class|series|member]_[selector|collector] and a rootkey,
 *  as arguments and returns a pointer to a keylist which contains entries
 *  for the parsed components.
 */
   static char local[MAX_STRLEN];
   static char key[MAX_STRLEN];
   KEY *return_list = newkeylist();
   char *range, *rightbracket;
   char *expanded_range = NULL;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}
   if (!collector) {soi_errno = MISSING_DATA_NAME; return NULL;}
   if (strlen(collector)>=MAX_STRLEN) {soi_errno=NAME_TOO_LONG; return NULL;}
   strcpy (local, collector);

   if (range = (char*)strchr (local, '[')) {
      *range = 0; range ++;
      if (!(rightbracket = (char*)strchr (range, ']'))) { 
         soi_errno = CLASS_SYNTAX_ERROR;
	 return NULL;
      }
      *rightbracket = 0;
      if (!(expanded_range = expand_range(range))) return NULL;
      sprintf (key, "%s_range", root);
      setkey_str (&return_list, key, expanded_range);
      if (!(strchr (expanded_range,','))) {
         sprintf(key, "%s_sn", root);
         setkey_int (&return_list, key, atoi(expanded_range));
      }
   }
   setkey_str (&return_list, root, local);
   return return_list;
}

static KEY *parse_range (char *range, char *root) {
/*
 *  takes two character strings, a record_range and a rootkey, as arguments 
 *  and returns a pointer to a keylist which contains entries for the parsed
 *  components. These entries can include first number, last number,
 *  and increment. 
*/
/* parses slice_range (with optional increment) e.g.0-15,3 */

   KEY *return_list = newkeylist();
   static char key[MAX_STRLEN];
   int n, fsn, lsn, incr=1;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}
   if (!range) {soi_errno = MISSING_DATA_NAME; return NULL;}

   if ((n=sscanf (range, "%d-%d,%d", &fsn, &lsn, &incr)) == 1) lsn = fsn;
   else if (n == 3) ; /* an increment was specified */
   else if (n != 2) {soi_errno = RANGE_SYNTAX_ERROR; return NULL;}

   sprintf (key, "%s_fsn", root);
   setkey_int (&return_list, key, fsn);

   sprintf (key, "%s_lsn", root);
   setkey_int (&return_list, key, lsn);

   sprintf (key, "%s_incr", root);
   setkey_int (&return_list, key, incr);

   return (return_list);
}

static KEY *parse_data_selector (char *selector, char *root) {
/*
 *  takes two character strings, a data_selector and a rootkey, as arguments
 *  and returns a pointer to a keylist which contains entries for the parsed
 *  components.
 *
 *  returns a keylist which can contain entries for
 *  data_selector with key root_selector, format with key root_fmt, 
 *  working directory with key root_wd,
 *  var_name or file_name with key root_select,
 *  first and last record numbers with keys root_fsn and root_lsn, and
 *  array_slice with key root_slice
 */
   KEY *return_list = newkeylist();
   KEY *range_list;
   static char local[MAX_STRLEN];
   static char key[MAX_STRLEN];
   static char default_range[] = "0--1"; /* 0 to -1 */
   char *fmt, *records;
   char *charptr = local;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}

   if (!selector) {
      sprintf (key, "%s_fmt", root);
      setkey_str (&return_list, key, DEFAULT_SN_FMT);
      range_list = parse_range (default_range, root);
      add_keys (range_list, &return_list);
      freekeylist (&range_list);
      return return_list;
   }

   if (strlen(selector)>=MAX_STRLEN) {soi_errno=NAME_TOO_LONG; return NULL;}
   strcpy (local, selector);

   sprintf(key,"%s_selector", root);
   setkey_str (&return_list, key, selector);

   sprintf (key, "%s_fmt", root);
   if (fmt=(char*)strstr (local,",fmt:")) {
      *fmt=0; fmt++; fmt=(char*)strchr(fmt,':'); fmt++;
      setkey_str (&return_list, key, fmt);
   }
   else setkey_str (&return_list, key, DEFAULT_SN_FMT);

   if (records = (char*)strchr (local, '[')) {
      *records=0; records++;
      if (charptr = (char*)strchr(records,']')) {*charptr=0; charptr++;}
      else {soi_errno = DATA_SELECTOR_SYNTAX_ERROR; return NULL;}
      
      if (strlen(records)) range_list = parse_range (records, root);
      else range_list = parse_range (default_range, root);
      if (range_list) {
         add_keys (range_list, &return_list);
         freekeylist (&range_list);
      }
      else return NULL;

      if (*(charptr) == '[') { /* there is an array_slice */
         sprintf (key, "%s_slice", root);
         setkey_str (&return_list, key, charptr);
      }
   }
   else {
      range_list = parse_range (default_range, root);
      add_keys (range_list, &return_list);
      freekeylist (&range_list);
   }

   sprintf (key, "%s_select", root);
   if (charptr=(char*)strrchr(local,'/')) { /* data selector contains a slash */
      charptr++; setkey_str (&return_list, key, charptr); *charptr=0;
      sprintf (key, "%s_wd", root); 
      setkey_str (&return_list, key, local);
   }
   else setkey_str (&return_list, key, local); /* could be empty string */

   return return_list;
}

static KEY *parse_datasuperset (char *superset, char *root) {
/*
 *  takes two character strings, a data*set_name and a rootkey,
 *  as arguments and returns a pointer to a keylist which contains entries
 *  for the parsed components.
 */
   static char local[MAX_STRLEN];
   static char key[MAX_STRLEN];

   KEY *return_list = newkeylist();
   KEY *sublist;
   char *classname, *selector;
   char *nextptr = local;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}
   if (!superset) {soi_errno = MISSING_DATA_NAME; return NULL;}
   if (strlen(superset)>=MAX_STRLEN) {soi_errno=NAME_TOO_LONG; return NULL;}
   strcpy (local, superset);

   sprintf (key, "%s_data", root);
   setkey_str (&return_list, key, local);

   while (nextptr) {
      classname = nextptr;
      if (nextptr = (char*)strchr (classname, ':')) {
        *nextptr = 0;
	nextptr++;
      } else {
        soi_errno = DATASUPERSET_SYNTAX_ERROR;
	return NULL;
      }

      if (!strlen (classname)) {
        soi_errno = CLASS_SYNTAX_ERROR;
	return NULL;
      }
      sprintf (key, "%s_%s", root, classname);

      selector = nextptr;
      if (strcmp (classname, "series")) { /* not = */
         if (nextptr = (char*)strchr(selector, ',')) {*nextptr = 0; nextptr++;}
      }
      else /* assume the rest of the name is the series collector */
      /* there can be commas in a collector so can't just look for comma */
         nextptr = NULL;

      if (!(sublist = parse_class (selector, key))) return NULL;
      add_keys (sublist, /*to*/ &return_list);
      freekeylist (&sublist);
   }
   return return_list;
}

KEY *parse_dataname (char *data, char *root) {
/*
 * takes two character strings, a datasubset_name and a rootkey, as arguments.
 * It separates the data and selector portions of the datasubset_name
 * and in turn parses whichever of these (or both) is found.
 * It returns a pointer to a keylist which contains entries
 * for the parsed components.
*/
   KEY *return_list = newkeylist();
   KEY *sublist;
   static char local[MAX_STRLEN];
   char *selector;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}
   if (!data) {soi_errno = MISSING_DATA_NAME; return NULL;}
   if (strlen(data)>=MAX_STRLEN) {soi_errno=NAME_TOO_LONG; return NULL;}
   strcpy (local, data);

   if ((local == (char*)strstr(local, "proj:")) ||
       (local == (char*)strstr(local, "prog:"))) {
      if (selector = (char*)strstr(local, ",sel:")) {
         *selector=0;
	 selector++; 
         selector=(char*)strchr(selector,':');
	 selector++;
      }

      if (sublist = parse_datasuperset (local, root)) {
         add_keys (sublist, /*to*/ &return_list);
	 freekeylist (&sublist);
      } else return NULL;
   }
   else selector = local;
   
   if (sublist = parse_data_selector(selector,root)) { 
      add_keys (sublist, /*to*/ &return_list);
      freekeylist (&sublist);
   } else return NULL;

   return return_list;
}

#define MAXXLEN	(16384)
/*
 *  The length of this constant controls the maximum length of the original
 *    dataset specifier, e.g. strlen ("prog:a,series:b[1000-1891],sel:[0-3,2]")
 *    It is only significant when the specifier is constructed as a long list
 *    of comma separated individual dataset name specifiers, rather than
 *    making use of indices
 */
static KEY * parse_datacollection (char *collection, char *root) {
/*
 *  takes two character strings, a datacollection_name and a rootkey,
 *  as arguments.  It parses the datacollection_name and returns a pointer
 *  to a keylist which contains entries for a ordered set of components
 *  of the datacollection_name.  
 */
   int nsets = 0;		       /*  number of datasets in collection  */
   static char local[MAXXLEN];
   static char key[MAXXLEN];
   static char sn_key[MAXXLEN];
   static char range_key[MAXXLEN];
   char *nextptr, *nptr, *datasubset, *range;
   char *charptr = local;
   KEY *data_list;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}
   if (!collection) {soi_errno = MISSING_DATA_NAME; return NULL;}

   /* globals required by keyiterate function */
   __collection_list=newkeylist(); 
   sprintf (__tailfmt, "%s_%s", root, "%s"); /* e.g "in_%s" */

   if (strlen (collection) >= MAXXLEN) {
     soi_errno = NAME_TOO_LONG;
     return NULL;
   }
   strcpy (local, collection);
   sprintf (range_key, "%s_series_range", root);

   while (charptr) {				    /*  for each datasubset  */
      datasubset = charptr;
      if (charptr = (char*)strchr (datasubset,';')) {
         *charptr = 0;
	 charptr++;
      }

      if (!strlen (datasubset)) {
        soi_errno = DATACOLLECTION_SYNTAX_ERROR;
	return NULL;
      }

      if(!(data_list = parse_dataname (datasubset, root))) return NULL;

      if (range=getkey_str (data_list, range_key)) {
         /* range should contain a comma separated list of series numbers */
         nextptr = range; 
         while (nextptr) { /* for each series number */
            nptr = nextptr;
            if (nextptr = (char*)strchr (nptr, ',')) {
               *nextptr = 0; nextptr++;
            }
            /* nptr now points to series number */
            sprintf(__keyfmt,"%s_%d_%s", root, nsets, "%s"); /* e.g. in_0_%s */
            sprintf(sn_key, __keyfmt, "series_sn");
            setkey_int (&__collection_list, sn_key, atoi(nptr));
            keyiterate (add_member, /* from */ data_list);
            nsets++;
         }
         free (range);
      }
      else { /* there is no series number */
         sprintf(__keyfmt,"%s_%d_%s", root, nsets, "%s"); /* e.g. in_0_%s */
         keyiterate (add_member, /* from */ data_list);
         nsets++;
      }
      if (nsets==1) add_keys (data_list, &__collection_list);
      freekeylist (&data_list);
   }

   sprintf (key, "%s_nsets", root);
   setkey_int (&__collection_list, key, nsets);

   return __collection_list;
}

static void replicate_name (char *name, KEY **inlist, char *root) {
/*
 *  takes a character string (name), a pointer to a keylist (inlist),
 *  and a rootkey as arguments. If an integer value for root_nsets exists
 *  in inlist and either a string or integer value for root_name exist in
 *  inlist the value is replicated in the list nsets times with keys
 *  root_0_name, ...  It is used by parse_list to insert a dbase entry
 *  for each element of a datacollection.
 *
 *  if an integer value for root_nsets exists in inlist and 
 *  either a string or integer value for root_name exist in inlist 
 *  the value is replicated in the list nsets times with keys root_0_name, ...
 */
   static char key[MAX_STRLEN];
   int n, nsets, keytype, intval;
   char *stringval;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return;}
   if (!name) return;

   sprintf (key, "%s_nsets", root);
   if (KEYTYP_INT != getkeytype (*inlist, key)) return;
   nsets = getkey_int (*inlist, key);

   sprintf (key, "%s_%s", root, name);
   if (KEYTYP_VOID == (keytype = getkeytype (*inlist, key))) return;

   if (keytype == KEYTYP_STRING) {
      stringval = getkey_str(*inlist, key);
      for (n=0; n<nsets; n++) {
         sprintf (key, "%s_%d_%s", root, n, name);
         setkey_str (inlist, key, stringval);
      }
      free (stringval);
   }
   else if (keytype == KEYTYP_INT) {
      intval = getkey_int (*inlist, key);
      for (n=0; n<nsets; n++) { 
         sprintf (key, "%s_%d_%s", root, n, name);
         setkey_int (inlist, key, intval);
      }
   }
   else return;  
}

void str_compress (char *s) {
   char *source, *dest;
   source = dest = s;
   while (*source) {
      if (!isspace (*source)) {
         *dest = *source;
         dest++;
      }
      source++;
   }
   *dest = 0;
}

void str_collapse (char *s, char c) {
/* collapses multiple occurrences of char c in string s to single occurrences */
   char *source, *dest;
   int skipit=0;
   source = dest = s;
   while (*source) {
      if (*source!=c) {
         skipit=0;
         *dest = *source;
         dest++;
      }
      else if (!skipit) {
         skipit=1;
         *dest = *source;
         dest++;
      }
      source++;
   }
   *dest = 0;
}

void int2numstr (char *numstr, char *fmt, int num) {
/* converts integer to string according to the given format 
 *   DOES NOT malloc the string
 */
   int nchars;
   char c;

   numstr[0] = '\0'; 
   if (2==sscanf(fmt,"%*c%d%c",&nchars,&c) && (c=='T')) { /* convert to time */
      TIME t = MISSION_EPOCH;
      if (num<3624) /* days */ t += 86400*num;
      else if (num>219150) /* minutes */ t += 60*num;
      else /* hours */ t += 3600*num; 
      sprint_at(numstr, t);
      if (nchars<strlen(numstr)) numstr[nchars] = '\0';
   }
   else /* assume a sprintf type format */
      sprintf(numstr,fmt,num);
}

int parse_list (KEY **params, char *root) {
/*
 *Finds an entry in the params list with the key root and parses it as
 * a datacollection. Entries for each of the parsed parts are inserted 
 * into the params list. Keys for these entries are of the form root_partname.
 */
   char *dataname;
   KEY *data_list = NULL;
   soi_errno = NO_ERROR; 

   if (!(dataname = getkey_str (*params, root))) 
      {soi_errno = MISSING_DATA_NAME; return soi_errno;}
			  /*  dataname now contains the string to be parsed  */

   str_compress (dataname);			     /*  remove white space  */
   if (!strlen(dataname)) 
      {soi_errno = MISSING_DATA_NAME; return soi_errno;}

   if (!(data_list = parse_datacollection (dataname, root))) return soi_errno;
   free (dataname);

   add_keys (data_list, /*to*/ params); 
   freekeylist (&data_list);

   replicate_name ("dbase", params, root);
   return add_paths (params, root);
}

char *fill_template (char *template, KEY *fromlist, char *root){

   /* template is applied to the root entries from the list
      to form the return string. if the template is NULL an empty string
      is returned.  if there is an error soi_errno is set and 
      a NULL pointer is returned  */

   static char local[MAX_STRLEN];
   static char buffer[MAX_STRLEN];
   static char key[MAX_STRLEN];
   static char fmt[MAX_STRLEN];
   static char member_number[MAX_STRLEN];
   char *tmp = local;
   char *another, *member_name, *left, *right;

   soi_errno = NO_ERROR; 
   if (!root) {soi_errno = MISSING_ROOT_KEY; return NULL;}
   if (template) 
      if (strlen(template)<MAX_STRLEN) strcpy (local,template);
      else {soi_errno = NAME_TOO_LONG; return NULL;}
   else tmp=NULL; 
   buffer[0] = 0; /*  start with empty string  */
/* note: strtok doesn't work because an empty string is not a token.
 *  thus a beginning { or back to back }{ in template causes problems. */
   while (tmp) {
      left = index (tmp, '{');
      if (left) {
         *left = 0;
         if (strlen(buffer)+strlen(tmp) >= MAX_STRLEN) 
            {soi_errno = NAME_TOO_LONG; return NULL;}
         strcat (buffer, tmp);
         tmp = ++left;
         right = index (tmp, '}');
         if (right) {
            *right = 0;
            if (tmp[0] == '#') { /*  name begins with #  */
	    /*  example:  {#series} or {#%03d#series}  */
               if (another = (char*)strchr (++tmp, '#')) {
                  *another = 0;			    /*  end format  */
                  sprintf (fmt, "%s", tmp);
                  tmp = ++another;
               }
               else { /*  use default format  */
                  sprintf (fmt, "%s", "%d");
               }
	       /*  tmp now points to class names  */
               sprintf (key, "%s_%s_sn", root, tmp);
               if (getkeytype (fromlist, key) == KEYTYP_INT) {
                  int2numstr (member_number, fmt, getkey_int (fromlist, key));
                  if (strlen(buffer)+strlen(member_number) >= MAX_STRLEN) 
                     {soi_errno = NAME_TOO_LONG; return NULL;}
                  strcat (buffer,member_number);
               }
               else { /* missing number */
                  soi_errno = CANNOT_FILL_TEMPLATE;
                  return NULL;
               }
            }
            else {
               sprintf (key, "%s_%s", root, tmp);
               if (member_name = getkey_str(fromlist, key)) {
                  if (strlen(buffer)+strlen(member_name) >= MAX_STRLEN) 
                     {soi_errno = NAME_TOO_LONG; return NULL;}
                  strcat (buffer, member_name); 
                  free (member_name);
               }
               else {  /* missing value  */
                  soi_errno = CANNOT_FILL_TEMPLATE;
                  return NULL;
               }
            }
            tmp = ++right;
         } 
         else {
            soi_errno = TEMPLATE_SYNTAX_ERROR;
            return NULL;  /*  syntax error - no matching bracket  */
         }
      }
      else {
         if (strlen(buffer)+strlen(tmp) >= MAX_STRLEN) 
            {soi_errno = NAME_TOO_LONG; return NULL;}
         strcat (buffer, tmp);
         tmp = NULL;
      }
   }
   return buffer;  
}

char *logname (KEY *params, char *root) {
/*
 *  Forms a .log filename from keylist parameters for the given root
 *    and returns a pointer to the malloced name
 */
static char tmp[MAX_STRLEN];
static char key[MAX_STRLEN];
char *extension;
char *name;

sprintf (key, "%s_select", root);		 /*  variable name  */
if (getkeytype (params, key) == KEYTYP_STRING) {
   sprintf (tmp, "%s", getkey_str(params, key));
   if (((extension = (char *)strrchr (tmp, '.')) != 0) &&
      (strcmp (extension, ".log") == 0)) {
      ;			 /*  all is well  -- tmp contains the logfile name  */
   } else {
      soi_errno = CANNOT_MAKE_NAME;
      return NULL;
   }
} else {
   sprintf (key, "%s_wd", root);	     /*  working directory  */
   if (getkeytype (params, key) == KEYTYP_STRING)
      sprintf (tmp, "%s", getkey_str(params, key));
   else
     tmp[0] = 0;

   sprintf (key, "%s_basename", root);		      /*  basename  */
   if (getkeytype (params, key) == KEYTYP_STRING)
      if (strlen (getkey_str (params, key)))
         sprintf (tmp, "%s%s", tmp, getkey_str(params, key));

   if (strlen (tmp)) {
      sprintf (tmp, "%s.log", tmp);
   } else {
      soi_errno = CANNOT_MAKE_NAME;
      return NULL;
   }
}

if (name = (char *)malloc (2*(strlen(tmp) + 1)))
   strcpy (name, tmp);
return name;
}

char *fitsname (KEY *params, char *root, int sn) {
/*
 *  Forms a .fits filename from keylist parameters for the given root and
 *    sn and returns a pointer to the malloced name
 */
static char tmp[MAX_STRLEN];
static char key[MAX_STRLEN];
static char format[MAX_STRLEN];
char *name;
char *fmt;
char *extension;
char *lastslash;

sprintf (key, "%s_wd", root);		     /*  working directory  */
if (getkeytype (params, key) == KEYTYP_STRING)
   if (strlen (getkey_str (params, key)))
      sprintf (tmp, "%s/", getkey_str(params, key));
else tmp[0] = 0;

sprintf (key, "%s_basename", root);		      /*  basename  */
if (getkeytype (params, key) == KEYTYP_STRING)
   if (strlen (getkey_str (params, key)))
      sprintf (tmp, "%s%s.", tmp, getkey_str(params, key));

sprintf (key, "%s_select", root);		 /*  variable name  */
if (getkeytype (params, key) == KEYTYP_STRING)
   if (strlen (getkey_str (params, key)))
      sprintf (tmp, "%s%s", tmp, getkey_str(params, key));

sprintf (key, "%s_fmt", root);
if (getkeytype (params, key) == KEYTYP_STRING)	    /* sn format is in list */
   fmt = getkey_str (params, key);
else
   fmt = DEFAULT_SN_FMT;

extension = (char *)strrchr (tmp, '.');
lastslash = (char *)strrchr (tmp, '/');

if ((extension != NULL) && (strcmp (extension, ".fits") == 0)) {
				      /*  basepath ends in .fits  -- use it  */
   ;
} else if (
   ((extension != NULL) && (strcmp (extension, ".") == 0)) ||
						      /*  basepath ends in . */
       ((lastslash != NULL) && (strcmp (lastslash, "/") == 0))) {
				       /*  basepath ends in / : add sn.fits  */
   sprintf (format, "%s%s.fits", "%s", fmt);
   sprintf (tmp, format, tmp, sn);
} else if (strlen (tmp)) {				   /*  add .sn.fits  */
   sprintf (format, "%s.%s.fits", "%s", fmt);
   sprintf (tmp, format, tmp, sn);
} else {					    /*  the name is sn.fits  */
   sprintf (format, "%s.fits", fmt);
   sprintf (tmp, format, sn);
}

if (name = (char *)malloc (2*(strlen(tmp) + 1)));
   strcpy (name, tmp);
return name;
}

char *fitsname_nopath (KEY *params, char *root, int sn) {
/*
 *  same as fitsname, but without the directory pathname
 */
  static char tmp[MAX_STRLEN], key[MAX_STRLEN], format[MAX_STRLEN];
  char *name, *fmt, *extension, *lastslash;

  sprintf (key, "%s_basename", root);
  if (getkeytype (params, key) == KEYTYP_STRING)
    if (strlen (getkey_str (params, key)))
      sprintf (tmp, "%s.", getkey_str(params, key));
    else tmp[0] = '\0';

  sprintf (key, "%s_select", root);
  if (getkeytype (params, key) == KEYTYP_STRING)
    strcat (tmp, getkey_str (params, key));

  sprintf (key, "%s_fmt", root);
  fmt = (getkeytype (params, key) == KEYTYP_STRING) ?
      getkey_str (params, key) : DEFAULT_SN_FMT;

  extension = (char *)strrchr (tmp, '.');
  lastslash = (char *)strrchr (tmp, '/');

  if ((extension != NULL) && (strcmp (extension, ".fits") == 0)) {
    ;
  } else if (((extension != NULL) && (strcmp (extension, ".") == 0)) ||
      ((lastslash != NULL) && (strcmp (lastslash, "/") == 0))) {
    sprintf (format, "%s%s.fits", "%s", fmt);
    sprintf (tmp, format, tmp, sn);
  } else if (strlen (tmp)) {
    sprintf (format, "%s.%s.fits", "%s", fmt);
    sprintf (tmp, format, tmp, sn);
  } else {
    sprintf (format, "%s.fits", fmt);
    sprintf (tmp, format, sn);
  }

  if (name = (char *)malloc (2*(strlen (tmp) + 1)));
    strcpy (name, tmp);
  return name;
}

char *fitsname_noseries (KEY *params, char *root, int sn)
/*
 *  Forms a .fits filename from keylist parameters for the given root and
 *    sn and returns a pointer to the malloced name
 *  This function is the same as fitsname except that 
 *    the basename (which is formed using the {series}.{#series} template)
 *    is not added to the filename.
 */
{
static char tmp[MAX_STRLEN];
static char key[MAX_STRLEN];
static char format[MAX_STRLEN];
char *name;
char *fmt;
char *extension;
char *lastslash;

sprintf (key, "%s_wd", root);		     /*  working directory  */
if (getkeytype (params, key) == KEYTYP_STRING)
   if (strlen (getkey_str (params, key)))
      sprintf (tmp, "%s/", getkey_str(params, key));
else tmp[0] = 0;

sprintf (key, "%s_select", root);		 /*  variable name  */
if (getkeytype (params, key) == KEYTYP_STRING)
   if (strlen (getkey_str (params, key)))
      sprintf (tmp, "%s%s", tmp, getkey_str(params, key));

sprintf (key, "%s_fmt", root);
if (getkeytype (params, key) == KEYTYP_STRING)	    /* sn format is in list */
   fmt = getkey_str (params, key);
else
   fmt = DEFAULT_SN_FMT;

extension = (char *)strrchr (tmp, '.');
lastslash = (char *)strrchr (tmp, '/');

if ((extension != NULL) && (strcmp (extension, ".fits") == 0))
/*  basepath ends in .fits  -- use it */
{
   ;
}
else if (
   ((extension != NULL) && (strcmp (extension, ".") == 0)) ||
/*  basepath ends in . */
   ((lastslash != NULL) && (strcmp (lastslash, "/") == 0))
/*  basepath ends in /  */
)
/*  add sn.fits */
{
   sprintf (format, "%s%s.fits", "%s", fmt);
   sprintf (tmp, format, tmp, sn);
}
else if (strlen (tmp))
/* add .sn.fits */
{
   sprintf (format, "%s.%s.fits", "%s", fmt);
   sprintf (tmp, format, tmp, sn);
}
else
/* the name is sn.fits */
{
   sprintf (format, "%s.fits", fmt);
   sprintf (tmp, format, sn);
}

if (name = (char *)malloc (2*(strlen(tmp) + 1)));
   strcpy (name, tmp);
return name;
}


char *tlmname (KEY *params, char *root)
/*
 *  if prog exists in params then
 *  forms a $wd/tfr filename from keylist parameters for the given root 
 *    and returns a pointer to the malloced name
 *  otherwise it assumes filename is in the selector
 */
{
static char tmp[MAX_STRLEN];
static char key[MAX_STRLEN];
char *name;

sprintf (key, "%s_prog", root);	           /*  progam name  */
if (getkeytype (params, key) == KEYTYP_STRING)
  {if (strlen (getkey_str (params, key)))	    /* program name exists  */
     {
	sprintf (key, "%s_wd", root);	     /*  working directory  */
	if (getkeytype (params, key) == KEYTYP_STRING)
   	if (strlen (getkey_str (params, key)))
      		sprintf (tmp, "%s/tfr", getkey_str(params, key));
	else tmp[0] = 0;
     }
  }
 else						/* use selector as filename  */
     {
	sprintf (key, "%s_selector", root);           /*  selector  */
	if (getkeytype (params, key) == KEYTYP_STRING)
   	if (strlen (getkey_str (params, key)))
      		sprintf (tmp, "%s", getkey_str(params, key));
	else tmp[0] = 0;
     }
 if (name = (char *)malloc (2*(strlen(tmp) + 1)));
     strcpy (name, tmp);

 return name;
}


char *datasetname (KEY *params, char *root)
{
/*
 *  Forms a name string from keylist parameters for the given root
 *    and returns a pointer to the malloced name; the datasetname is
 *    just the normal full data set name (prog:...) without the possible
 *    fmt: field at the end.
 */
  static char key[MAX_STRLEN];
  char *name, *tmp;

  sprintf (key, "%s_fmt", root);
  if (getkeytype (params, key) == KEYTYP_STRING)  /*  sn format is in list  */
  {
    tmp = getkey_str (params, root);
    if (name = (char *)malloc (2 * (strlen (tmp) + 1)))
    {
      strcpy (name, tmp);
      tmp = (char*)strstr (name, ",fmt:");
      if (tmp)
        tmp[0] = 0;
    }
    return (name);
  }
  else
    return (getkey_str (params, root));
}
