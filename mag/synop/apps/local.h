/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/local.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: local.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.3  90/12/18  16:51:07  phil
 * _ to __
 * 
 * Revision 1.2  90/12/18  16:40:04  phil
 * ANSI symbol removed.
 * 
 * Revision 1.1  90/12/18  15:49:49  phil
 * Initial revision
 * 
*/

#ifndef _LOCAL_INCLUDED
#define _LOCAL_INCLUDED

/********************************************************************

  The following macro symbols are predefined by the compilers:
  On Banneker and other vax ultrix 3.x systems:

	unix           Any UNIX system
	bsd4_2         Berkeley UNIX Version 4.2
	ultrix         ULTRIX only
	vax            VAX only (as opposed to PDP-11)
	sparc          SUN sparc processor

  On Quake and other mips ultrix 4.x systems:

	unix           Any UNIX system
	bsd4_2         Berkeley UNIX Version 4.2
	ultrix         ULTRIX only
	mips           Any MIPS architecture
	MIPSEL         Little endian variant of MIPS architecture
	host_mips      Native compilation environment (as opposed to
                         cross-compiler)

  On mips big endian machines:

	MIPSEB         Big endian - e.g. sgi

  On the NeXT:

	NeXT           The NeXT machine
	unix
	BSD
	mc68000

  On the system V 386 machines:

	i386	       intel 80386 machine

  On the sgi irix machines:

        __sgi		sgi machines

  On the Sun sparc machines:

        __sun
        __sparc

  On ANSI standard C:

	__FILE__		text file name - char*
	__LINE__		line number - int
	__DATE__		compile date - char*
	__TIME__		compile time - char*
	__STDC__		ANSI standard C

********************************************************************/

#ifndef MACHINE

#ifdef MIPSEL
#define MACHINE mipsel
#define MACHINE_TXT "mipsel"
#endif
#ifdef MIPSEB
#define MACHINE mipseb
#define MACHINE_TXT "mipseb"
#endif
#ifdef vax
#define MACHINE vax
#define MACHINE_TXT "vax"
#endif
#ifdef NeXT
#define MACHINE next
#define MACHINE_TXT "next"
#endif
#ifdef i386
#define MACHINE i386
#define MACHINE_TXT "i386"
#endif
#ifdef __sparc
#define MIPSEB
#define MACHINE mipseb
#define MACHINE_TXT "sparc"
#endif
#ifndef MACHINE
#define MACHINE unknown
#define MACHINE_TXT "unknown"
#endif

#endif

#ifndef MACHINE_TXT
#define LITERAL_STRING(X) # X
#define MACHINE_TXT LITERAL_STRING(MACHINE)
#endif

#endif
