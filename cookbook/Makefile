include Makevars/Makevars_$(JPLAT).mk

DEST =	.

INCL = -I$(DRMS)/base/include
LIBD = -L$(DRMS)/lib/$(JPLAT)

DLIBS = $(LIBD) -ldrms -lpq -lz -lpthread -lcfitsio -lm

.SECONDEXPANSION:

# simple modules, no extra libraries
MODS =        smpl_00

FPLS =	$(addprefix $(DEST)/, $(MODS))

all:	$(MODS)

$(MODS):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $@.o $(DLIBS)
	$(RM) $@.o

clean:
	$(RM) *.o

Clean:
	$(RM) *.o
	$(RM) $(FPLS)
.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL)
