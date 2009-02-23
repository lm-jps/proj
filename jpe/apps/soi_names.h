/*
 *  soi_names.h			~soi/(version)/include/soi_names.h
 *
 *  The programs referenced here provide naming services that act as
 *    intermediaries between keys provided by catalog database queries
 *    or other information and pointers to actual files.
 *  C++ programs should include the corresponding file "soi_names.hxx".
 *
 *  The functions provided implement the SOI dataset naming scheme,
 *    described in SOI Technical Notes TN-93-104.
 *
 *  Additional information is in the following man pages:
 *	names (3)
 *
 *  Responsible:  Kay Leibrand		KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *	soi_names.hxx does not exist
 *	fitsname_noseries should become the default; naming conventions
 *	  need to be formalized
 *
 *  Revision history is at the end of the file.
 */

#ifndef SOI_NAMES_INCL
/****************************************************************************/
/**************************  INCLUDE STATEMENTS  ****************************/
/****************************************************************************/

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

#include <soi_key.h>

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define SOI_NAMES_VERSION_NUM	(5.0)

#define MAX_STRLEN 16384
/*
 *  The length of this constant controls the maximum number of dataset
 *    numbers that can be parsed by creating a string of the form
 *    "n0,n1,n2,..,nn".  The maximum number of parseable numbers is:
 *		4096	8192   16384
 *    <4 digit:	1041
 *    4 digit:	 819	1638	3276
 *    5 digit:	 682	1365	2730
 *    6 digit:	 585	1170	2340
 */
#define DEFAULT_SN_FMT "%04d"
#define DEFAULT_RECORD_RANGE "0"

#define SOI_NAMES_INCL	1

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/
  
/****************************************************************************/
/*********************  GLOBAL & EXTERN DECLARATIONS  ***********************/
/****************************************************************************/

/****************************************************************************/
/****************************  MACRO DEFINITIONS  ***************************/
/****************************************************************************/

/****************************************************************************/
/***************************  FUNCTION PROTOTYPES  **************************/
/****************************************************************************/

				  /*  source file: ~CM/src/libast.d/names.c  */
extern int parse_list (KEY  **params, char * root);
extern char *fitsname (KEY *params, char *root, int sn); 
extern char *fitsname_nopath (KEY *params, char *root, int sn); 
extern char *fitsname_noseries (KEY *params, char *root, int sn); 
extern char *logname (KEY *params, char *root); 
extern char *tlmname (KEY *params, char *root);
extern char *datasetname (KEY *params, char *root); 
extern void str_compress (char *s);
extern void str_collapse (char *s, char c);
extern void int2numstr (char *numstr, char *fmt, int num);
extern char *fill_template (char *template, KEY *fromlist, char *root);
			      /*  source file: ~CM/src/libast.d/parse_arg.c  */
extern int parse_array (KEY **params, char *root, int type);
extern int parse_numerated (char *klist, char ***names);

#endif

/*
 *  Revision History
 *  V 0.0   93.02.25	R Bogart	created this file
 *  V 0.1   93.04.27	K Leibrand	created original function set			
 *  V 0.2   93.08	K Leibrand	deleted form_basepath,
 *		form_name, and form_fitsname; added fitsname and logname
 *  V 0.7   93.09.10	R Bogart	defined version number
 *	    93.12.10	R Bogart	added datasetname()
 *	    94.01.25	R Bogart	added declaration for
 *		fitsname_noseries, changed all declarations (except last)
 *		to extern
 *	    94.02.08	Katie Scott     added tlmname
 *  V 0.8   94.02.08	K Leibrand	changed VERSION_NUM
 *	    94.08.01	K Leibrand	added functions to parse data
 *		collections
 *	    94.08.03	K Leibrand	added replicate_name function so that
 *		non-dataname parts can be used in forming working directories
 *		for multiple datasets, e.g. dbase; added function add_paths
 *		to simplify directory & basename handling
 *	    94.08.24	K Leibrand	added function add_env_template
 *  v 0.9   94.11.14	R Bogart	changed VERSION_NUM
 *  v 4.5   99.11.16	R Bogart	added declaration from parse_arg
 *	99.12.14	R Bogart	added declaration for fitsname_nopath
 *	00.06.17	R Bogart	added declaration of parse_numerated
 * 	01.11.13	R Bogart	removed decalaration of obsolete
 *		function parse_list_force; quadrupled value of MAX_STRLEN to 
 *		16384; removed several declarations of functions internal
 *		to names.c
 *  v 5.0  04.08.05	R Bogart	modified declaration of parse_numerated()
 */

/*
$Id: soi_names.h,v 1.1 2009/02/23 22:49:03 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/soi_names.h,v $
$Author: production $
 * $Log: soi_names.h,v $
 * Revision 1.1  2009/02/23 22:49:03  production
 * initial
 *
 * Revision 1.21  2004/08/16  16:59:35  rick
 * *** empty log message ***
 *
 * Revision 1.20  2004/08/16 16:59:06  rick
 * see above
 *
 * Revision 1.19  2001/11/13 22:41:38  rick
 * see above
 *
 * Revision 1.18  2001/11/13  17:31:48  rick
 * see above
 *
 * Revision 1.17  2001/04/11  19:07:44  rick
 * see above
 *
 * Revision 1.16  1999/12/14  19:36:50  rick
 * see above
 *
 * Revision 1.15  1999/11/18  01:00:45  rick
 * see above
 *
 * Revision 1.14  1999/02/12  20:34:56  kay
 * change collapse to str_collapse and compress to str_compress
 * in response to SSTR 39
 *
 * Revision 1.13  1998/07/29  22:27:27  jim
 * noop
 *
 * Revision 1.12  1997/04/16  21:54:27  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.11  1996/07/02  22:02:57  kay
 * added prototype for collapse
 *
 * Revision 1.10  1996/04/01  22:44:37  kay
 * added externs for previously static functions in names.c
 *
 * Revision 1.9  1995/07/17  20:12:58  kay
 * set MAX_STRLEN to 4096
 *
 * Revision 1.8  1995/06/29  21:10:55  kay
 * increased MAX_STRLEN to 512
 *
 * Revision 1.7  1994/11/14  23:05:40  rick
 * V 0.9
 *
 * Revision 1.6  1994/08/24  20:36:10  kay
 * added function add_env_template
 *
 * Revision 1.5  1994/08/08  19:03:15  kay
 * added function add_paths to simplify directory & basename handling
 *
 * Revision 1.4  1994/08/03  20:30:26  kay
 * added replicate_name function so that non-dataname parts can be used in
 * forming working directories for multiple datasets, e.g. dbase
 *
 * Revision 1.3  1994/08/03  20:28:26  kay
 * added template error codes
 *
 * Revision 1.2  1994/08/01  22:32:22  kay
 * added functions to parse data collections
 *
 * Revision 1.1  1994/02/16  23:21:14  CM
 * Initial revision
 *
 */
