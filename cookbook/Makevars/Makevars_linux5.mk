# @(#)Makevars_linux_x86_64.mk
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
ICC =		/usr/local/bin/icc -mcmodel=medium
CC =		$(ICC)
CDEBUG =	-g
CDEFINES =	-DBETA
CFLAGS =	-O3 -std=c99 -D_FILE_OFFSET_BITS=64 -xW -D_GNU_SOURCE
#CFLAGS =	-std=c99 -xW -D_GNU_SOURCE
CCFLAGS =       -c $(CFLAGS)
CCGFLAGS =	-g -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
CINCLUDES =	-I$(JSOC)/base/include -I/usr/include
DBCC =		$(CC) $(CDEBUG)
FW =		-fullwarn
MULTI =
HFILES =
CFILES =

# f77 - FORTRAN
FC =		/usr/bin/f77
#FC =		/usr/local/bin/ifort
FFLAGS =	-O -c
FLIBS =		-lCm -lftn -lF77 -lm -lU77 -lI77 -lblas -lisam -lm

# ld - link editor
LD =		/usr/bin/ld
LDFLAGS =	-xW
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
AWK =		/bin/awk
CD =		cd
CHGRP =		/bin/chgrp
CHMOD =		/bin/chmod
CHOWN =		/bin/chown
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
TOUCH =		/bin/touch
XSTR =		$(ECHO)

# Revision History:
#
# 07.04.19	created this file, based on ~rick/genmake.d/Makevars_linux4.mk
# 07.05.29	modified to support JSOC module builds
# 09.01.20	modified for linux5
#
