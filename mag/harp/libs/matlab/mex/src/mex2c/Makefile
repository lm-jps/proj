# 
# Makefile for mex2c libraries
# 
# requires MEXEXT set: via the environment, or by make MEXEXT=mexFOO
# this will done from top-level makes

# turmon aug 99, dec 2002, oct 2005

JSOC_ROOT = /home/jsoc

# boilerplate
SHELL = /bin/sh
INSTALL = install -cv
INSTALLDATA = $(INSTALL) -m 644
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

# turmon -- I built this, and it works with Matlab r14 (old matlab)
# CC = $(HOME)/unix/gcc-3.3.1/bin/gcc
# turmon -- this might work with Matlab r2010
CC = /usr/local/gcc43/bin/gcc

# installation locations
MISCLIBDIR = $(TOP)/lib/$(MEXEXT)  # arch-dependent libs
MISCINCDIR = $(TOP)/include

# targets
LIBMEX2C := $(OUTDIR)/$(CSDIR)/libmex2c_lib.a
LIBMEX2MATL := $(OUTDIR)/$(CSDIR)/libmex2matl.a
LIBRARY = $(LIBMEX2C) $(LIBMEX2MATL)
INCLUDE = mex.h

# includes (need other package dirs)
#CPPFLAGS += -I$(JSOC_ROOT)/include -I../util
CPPFLAGS += $(CFITSIOH) -I$(OUTDIR)/$(CSDIR)/../util

# default target
all:$(LIBRARY)


# libraries
#libmex2c_lib.a: libmex2c_lib.a(mex2c_main_lib.o mex_stubs.o mex2c_util.o)
LIBMEX2C_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, mex2c_main_lib.o mex_stubs.o mex2c_util.o)
$(LIBMEX2C):	$(LIBMEX2C_OBJ)
	ar crus $@ $^

#libmex2matl.a:  libmex2matl.a(mex2matlab.o)
LIBMEX2MATL_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, mex2matlab.o)
$(LIBMEX2MATL):	$(LIBMEX2MATL_OBJ)
	ar crus $@ $^

$(LIBMEX2C_OBJ)	: $(notdir $(LIBMEX2C_OBJ:%.o=%.c)) $(OUTDIR)/$(CSDIR)/../util/mexhead.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

$(LIBMEX2MATL_OBJ)	: $(notdir $(LIBMEX2MATL_OBJ:%.o=%.c))
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

# instructions for install
install: $(LIBRARY)
	ranlib $(LIBRARY)
	mkdir -p $(MISCINCDIR)
	mkdir -p $(MISCLIBDIR)
	$(INSTALLDATA) $(INCLUDE) $(MISCINCDIR)
	$(INSTALLDATA) $(LIBRARY) $(MISCLIBDIR)

# extra dependences for mex_stubs.c
$(addprefix $(OUTDIR)/$(CSDIR)/, mex_stubs.o): mex_stubs.c mex_stubs_api.c mex_stubs_unimp.c

$(addprefix $(OUTDIR)/$(CSDIR)/, mex2c_main_lib.o): mex2c_main_lib.c mex2c_dispatcher.c

# instructions for clean
clean:
	-	rm -f *.o
	rm -f $(LIBRARY)

distclean:   clean
	( cd $(MISCLIBDIR); rm -f $(LIBRARY) )
	( cd $(MISCINCDIR); rm -f $(INCLUDE) )

