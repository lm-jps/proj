#
#	Generic per-directory makefile
#
# mjt 7 oct 1996, nov 2002, oct 2005, sep 2009

# needs MEXEXT defined, e.g. thru environment or make invocation
# needs JSOC_MACHINE defined

JSOC_ROOT = /home/jsoc

# clear out rules
.SUFFIXES:

# Put this target first so the default target is "all"
# (otherwise, targets within MexDocRules.mk will appear first)
default: all

install: all

# generic mex rules
include ../MexCodeRules.mk
# don't bother to remake the makefiles
Makefile Local.mk ../MexCodeRules.mk : ;

# common code
COMMON = $(OUTDIR)/$(CSDIR)/../../mex

# Flags
# CMEXPPFLAGS is passed to mex, which places the -I$(COMMON)/include
# *after* its own -I$(MATLABPATH)/extern/include, ensuring that 
# the Mathworks mex.h, not mine, is picked up.  The $$CFLAGS tells
# mex to add its own default CFLAGS at the end of mine.  See
# mex -help and mex -v for more.  To override the annoying -ansi
# put in by some mex scripts, can add -std=gnu89 at end of CMEXPPFLAGS.
#
# JSOC_ROOT/ dependences are for cfitsio libraries and includes
CMEXPPFLAGS += 'CFLAGS=-I$(COMMON)/src/util $$CFLAGS' -I$(JSOC_ROOT)/include
CMEXLDFLAGS = -L$(COMMON)/src/util -L$(COMMON)/src/mex2c -L$(JSOC_ROOT)/lib/$(JSOC_MACHINE) $(LIBSUSED)
LIBSUSED = -lmexrng -lmextoolsMW -lmex2matl 
MEX2C_LIBS = -lmexrng -lmextools -lm  # cf LIBSUSED, for mex2c only
LDFLAGS += -L$(COMMON)/src/util -L$(COMMON)/src/mex2c
CPPFLAGS += -I$(COMMON)/src/util -I$(COMMON)/src/mex2c

# Targets
# find mex source code
# need trickery to avoid trouble when there is no .c file 
#  ...in which case the glob done by sh in "egrep -l mexFunction *.c" 
#  matches nothing and "*.c" expands to the literal filename "*.c" which
#  does not exist, causing egrep to give an error...sigh
# later, modified to work with gnu make - 
#  would like to match mexFunction( but ( causes trouble
MEX_SRC := $(shell egrep -l '^mexFunction' `find . -maxdepth 1 -name '*.c'` /dev/null;  exit 0)
MEX_ABBRS = $(MEX_SRC:.c=.m)
MEX_FILES = $(MEX_ABBRS:.m=.$(MEXEXT))
INC_FILES = $(MEX_ABBRS:.m=.h)
CLEAN_FILES =

# per-directory modifications (`-' means no complaints if file is not there)
-include Local.mk

all: $(MEX_FILES)

distclean: clean

clean:
	rm -f $(MEX_FILES) $(CLEAN_FILES)

# these are not filename targets
.PHONY: all install clean distclean doc default 

