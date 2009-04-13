# @(#)Makevars_mac_osx_ia32.mk	
#  Revision history is at the end of the file.

# generic variable and macro definitions for Mac OS-X (10.5.6) with Intel CPU
#
#  Revision history is at the end of the file.

# make
MAKE =		/usr/bin/make

# the Shell: this must be a valid shell, e.g. sh or csh or tcsh; echo also
#   works as a no-op
SHELL =		/bin/sh

# commands and arguments

# ar - archive, ranlib
AR =		/usr/bin/ar
ARFLAGS =	crv
RANLIB =	/usr/bin/ranlib

# as - assembler
AS =		/usr/bin/as
ASFLAGS =

# cc - C compiler
NCC = 		/usr/bin/cc
GCC =		/usr/bin/gcc
ICC =		/usr/local/bin/icc
CC =		$(GCC)
CDEBUG =	-g
CDEFINES =	-DBETA
CFLAGS =	-O3 -std=c99 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
CCFLAGS =       -c $(CFLAGS)
CCGFLAGS =	-g -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
CINCLUDES =	-I$(JSOC)/base/include
DBCC =		$(CC) $(CDEBUG)
FW =		-fullwarn
MULTI =
HFILES =
CFILES =

# f77 - FORTRAN
FC =		/usr/bin/f77
FFLAGS =	-O -w -c
FLIBS =		-lCm -Xlinker -rpath -Xlinker /opt/intel_fc_80/lib -L/opt/intel_fc_80/lib -lifport -lifcore -lmkl_ia32 -lguide

# ld - link editor
LD =		/bin/ld
LDFLAGS =	
LDCMD =		$(CC) $(LDFLAGS)
MALCHK = 
LOADLIBES =
PARLIBS =	
OBJS =
AOUT =		noaout
AOUTS =		noaouts

# lex - lexical analyser
LEX =		/usr/bin/lex
LFLAGS =

# yacc
YACC =		/usr/bin/yacc
YFLAGS =

# oracle for platform
ORACLE_HOME=
ORAADD=

# rpc
RPCLIB = 

# miscellaneous commands
AWK =		/usr/bin/awk
CD =		cd
CHGRP =		/usr/bin/chgrp
CHMOD =		/bin/chmod
CHOWN =		/usr/sbin/chown
CP =		/bin/cp
CPP =		/usr/bin/cpp
DATE =		/bin/date
ECHO =		echo
INSTALL =	/usr/bin/install
LN =		/bin/ln
LS =		/bin/ls
MKDIR =		/bin/mkdir
MV =		/bin/mv
RM =		/bin/rm -f
STRIP =		/usr/bin/strip
TOUCH =		/usr/bin/touch
XSTR =		$(ECHO)

# Revision History:
#
# 09.04.13	created this file, based on Makevars_linux_ia32.mk
#
