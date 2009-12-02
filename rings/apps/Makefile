#
#  Common Makefile requiring re-compilation
#

include $(MAKEROOT)/Makevars_$(JPLAT).mk

INCL = -I$(DRMS)/base/include
LIBD = -L$(DRMS)/lib/$(JPLAT) -L/home/jsoc/lib/$(JPLAT)

DLIBS = $(LIBD) -ldrms -lpq -lz -lpthread -lcfitsio -lfftw3 -lfftw3f -lm

DEST = ../bin/$(JPLAT)

MODS =	mtrack pspec3 rdfitc

FPLS =	$(addprefix $(DEST)/, $(MODS)) $(addprefix $(DEST)/, rdvinv)

all:	$(MODS) rdvinv

.SECONDEXPANSION:

$(MODS):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $? $(DLIBS)
	$(RM) $@.o

rdvinv:		$$@.o ola_xy.o
	$(FC) -nofor_main -openmp -o $(DEST)/$@ $? $(DLIBS) -lmkl_em64t
	$(RM) $@.o

clean:
	$(RM)  *.o
	$(RM) $(FPLS)

.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL)

.f.o:
	$(FC) $(FFLAGS) $(@:.o=.f)
