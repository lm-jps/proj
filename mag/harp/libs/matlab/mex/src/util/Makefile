# 
# Makefile for numerical utilities
# 
# requires MEXEXT set: via the environment, or by make MEXEXT=mexFOO
# this will done from top-level makes

# boilerplate
SHELL = /bin/sh
TOP = ../..
# 64-bit versus 32-bit compilation
ifneq (,$(findstring 64,$(MEXEXT)))
        MODE=64
else
        MODE=32
endif
# suppress Makefile re-make
# need -fexceptions to support C++ error handling when called from Matlab
CFLAGS = -m$(MODE) -fPIC -Wall -O -fexceptions
# CFLAGS = -m$(MODE) -fPIC -Wall -g -fexceptions
# need "mxArray" definitions for mextool.c
CPPFLAGS = -I../mex2c
# but for mex, do *not* want our own mxArray defs
CMEXPPFLAGS = 

# turmon -- I built this, and it works with Matlab r14 (old matlab)
# CC = $(HOME)/unix/gcc-3.3.1/bin/gcc
# turmon -- this might work with Matlab r2010
CC = /usr/local/gcc43/bin/gcc


# installation locations
MISCLIBDIR = $(TOP)/lib/$(MEXEXT)  # arch-dependent libs
MISCINCDIR = $(TOP)/include

# targets (generated includes not handled the same as ordinary ones)
LIBEXRNG	:= $(OUTDIR)/$(CSDIR)/libmexrng.a
LIBMEXTOOLS	:= $(OUTDIR)/$(CSDIR)/libmextoolsMW.a
METAHEADER	:= $(OUTDIR)/$(CSDIR)/mexhead.h
LIBRARY    	= $(LIBEXRNG) $(LIBMEXTOOLS)
INCLUDE_GEN= mexhead.h
INCLUDE    = $(INCLUDE_GEN) rng.h

# some includes are generated automatically
TARGETS    = $(LIBRARY) $(INCLUDE_GEN)

# default target (some includes are generated automatically)
all: $(TARGETS)

install: all
	ranlib $(LIBRARY)
	mkdir -p $(MISCINCDIR)
	mkdir -p $(MISCLIBDIR)
	cp -f $(INCLUDE) $(MISCINCDIR)
	cp -f $(LIBRARY) $(MISCLIBDIR)

clean:
	-	rm -f *.o
	rm -f $(LIBRARY)

distclean:   clean
	( cd $(MISCLIBDIR); rm -f $(LIBRARY) )
	( cd $(MISCINCDIR); rm -f $(INCLUDE) )
	rm -f $(INCLUDE_GEN)

# libraries
LIBEXRNG_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, rng.o)
$(LIBEXRNG):	$(LIBEXRNG_OBJ)
	ar crus $@ $^

MEXCOMP_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, mextoolMW.o mexargcheckMW.o num_stringsMW.o)
GCCCOMP_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, ieee_consts.o)
LIBMEXTOOLS_OBJ	:= $(GCCCOMP_OBJ) $(MEXCOMP_OBJ)

$(LIBMEXTOOLS): $(LIBMEXTOOLS_OBJ)
	ar crus $@ $^

$(LIBEXRNG_OBJ) :	$(notdir $(LIBEXRNG_OBJ:%.o=%.c))
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

$(GCCCOMP_OBJ) :	$(notdir $(GCCCOMP_OBJ:%.o=%.c))
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

# this is a little raw, but it works for the purpose of making
# a pure-mex library set.
# The mex'ing lines below must not pick up my own mex.h!
# In my mathworks mexfiles, I presently use: 
#   -lmexrng -lmextools -lmex2matl
# -lmextoolsMW should -lmextools in mex compilations
# -lmexrng is 100% selfcontained, so does not need special handling
#           for the mex environment vs. elsewhere.
# -lmex2matl must include mex.h, but does not rely on it.

#mextoolMW.o: mextool.c
#	mex $(CMEXPPFLAGS) -c $< 
#	@ mv `echo $< | sed 's/\.c/.o/'` $@

#mexargcheckMW.o: mexargcheck.c
#	mex $(CMEXPPFLAGS) -c $< 
#	@ mv `echo $< | sed 's/\.c/.o/'` $@

#num_stringsMW.o: num_strings.c
#	mex $(CMEXPPFLAGS) -c $< 
#	@ mv `echo $< | sed 's/\.c/.o/'` $@

$(MEXCOMP_OBJ): %MW.o:	%.o
	mv `echo $@ | sed 's/\MW.o/.o/'` $@

$(MEXCOMP_OBJ:%MW.o=%.o): $(OUTDIR)/$(CSDIR)/%.o:	%.c
	mex $(CMEXPPFLAGS) -outdir $(OUTDIR)/$(CSDIR) -c $<

# mexhead.m is a concatenation of several header files.
$(METAHEADER): mexhead.m4 \
		mexargcheck.h mextool.h mex_api_extras.h \
		ieee_consts.h num_strings.h
	m4 mexhead.m4 > $(METAHEADER)
