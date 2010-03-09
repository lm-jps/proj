#
#  Common Makefile requiring re-compilation
#

# MAKEROOT = $(DRMS)/proj/cookbook/Makevars
include $(MAKEROOT)/Makevars_$(JPLAT).mk

INCL = -I$(DRMS)/base/include
LIBD = -L$(DRMS)/lib/$(JPLAT) -L/home/jsoc/lib/$(JPLAT)

DLIBS = $(LIBD) -ldrms -lpq -lz -lpthread -lcfitsio -lm

DEST = ../bin/$(JPLAT)

MODS =	drms_rebin

FPLS =	$(addprefix $(DEST)/, $(MODS))

all:	$(MODS)

.SECONDEXPANSION:

$(MODS):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $? $(DLIBS)
	$(RM) $@.o

clean:
	$(RM)  *.o
	$(RM) $(FPLS)

.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL)

.f.o:
	$(FC) $(FFLAGS) $(@:.o=.f)
