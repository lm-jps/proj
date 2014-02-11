/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/keywords.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: keywords.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.1  90/10/12  11:58:02  phil
 * Initial revision
 * 
 * Revision 1.3  90/10/03  22:39:57  phil
 * convert kw_insert to function from macro, change call sequence
 * 
 * Revision 1.2  90/09/19  14:13:49  phil
 * name change from gdds to gds
 * 
 * Revision 1.1  90/09/17  15:47:31  phil
 * Initial revision
 * 
*/

/* keywords.h */

#ifndef KEYWORDS_DEF
#define KEYWORDS_DEF

typedef struct keywords
	{
	struct keywords *next;
	char		*key;
	char		*value;
	} KEYWORDS;

extern KEYWORDS *kw_search();
extern int	 kw_delete();
extern char	*kw_value();
extern void	kw_insert();

#ifndef NULL
#define	NULL	((char *)0)
#endif

/*
 * KEYWORDS is a general token-value storage structure.  The set of
 * routines kwXXX provides tools to create and use these lists.
 * The contents of key and value are of no interest to the kw routines.
 * 
 * next	Pointer to the next token-value pair in the list.
 * 	"next" is (char *)NULL in the last item in a list.
 * 
 * key	Pointer to a null terminated printable string of ASCII characters.
 * 	The string may not contain white space.
 * 
 * value	Pointer to a null terminated string of ASCII characters.
 * 
 */

#endif
