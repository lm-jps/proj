#
#  Common Makefile requiring re-compilation
#

# MAKEROOT = $(DRMS)/proj/cookbook/Makevars
include $(MAKEROOT)/Makevars_$(JPLAT).mk

INCL = -I$(DRMS)/base/include
LIBD = -L$(DRMS)/lib/$(JPLAT) -L/home/jsoc/lib/$(JPLAT)

SLIBS = $(LIBD) -lcmdparams -ldstruct -ltimeio -lmisc -lm
DLIBS = $(LIBD) -ldrms -lpq -lz -lpthread -lcfitsio -lfftw3 -lfftw3f -lm

DEST = ../bin/$(JPLAT)

PROG =	gentargs

MODS =	mtrack pspec3 rdfitc

FPLS =	$(addprefix $(DEST)/, $(MODS)) $(addprefix $(DEST)/, $(PROG)) \
	$(addprefix $(DEST)/, rdvinv)

all:	$(MODS) $(PROG) rdvinv

.SECONDEXPANSION:

$(PROG):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $? $(SLIBS)
	$(RM) $@.o

$(MODS):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $? $(DLIBS)
	$(RM) $@.o

rdvinv:		$$@.o ola_xy.o
	$(FC) -nofor_main -openmp -o $(DEST)/$@ $? $(DLIBS) -lmkl_em64t
	$(RM) $?

clean:
	$(RM)  *.o
	$(RM) $(FPLS)

.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL)

.f.o:
	$(FC) $(FFLAGS) $(@:.o=.f)
