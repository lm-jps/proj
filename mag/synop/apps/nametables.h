/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/nametables.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: nametables.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.1  90/09/17  14:32:22  phil
 * Initial revision
 * 
*/

#ifndef NAMETABLES_DEF

#define NAMETABLES_DEF

typedef struct
{
	char	*nm;
	int	action;
} nametable;

extern nametable yesno[],kindnames[],accnames[],statenames[],
	errnames[],typenames[],resnames[];

/* See nametables.c for table initializations. */

#endif
