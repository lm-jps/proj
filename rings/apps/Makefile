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

MODS =	datavg maicalc mtrack pspec3 rdcover rdfitc

FOBJ =	cart_to_polar.o fourier_filter.o ring_pass.o polyval.o \
	multfactor.o func1.o first_deriv_1.o second_deriv.o \
	dogleg3.o hmifits.o ringanalysis.o ccint2.o fourierlibn32.o 

FPLS =	$(addprefix $(DEST)/, $(MODS)) $(addprefix $(DEST)/, $(PROG)) \
	$(addprefix $(DEST)/, rdvinv) $(addprefix $(DEST)/, rdsinv) \
	$(addprefix $(DEST)/, rdfitf)

all:	$(MODS) $(PROG) rdfitf rdvinv

FC = ifort11
F90F = -c -O3

.SECONDEXPANSION:

$(PROG):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $? $(SLIBS)
	$(RM) $@.o

$(MODS):	$$@.o
	$(LDCMD) -o $(DEST)/$@ $? $(DLIBS)
	$(RM) $@.o

rdfitf:		$$@.o $(FOBJ)
	$(FC)  -nofor-main -openmp -o $(DEST)/$@ $? $(FOBJ) $(DLIBS) -lmkl_em64t
	$(RM) $?

rdsinv:		$$@.o h1h2.o
	$(FC) -nofor_main -openmp -o $(DEST)/$@ $? $(DLIBS) -lmkl_em64t
	$(RM) $?

rdvinv:		$$@.o ola_xy_v11.o ola_subs.o
	$(FC) -nofor_main -openmp -o $(DEST)/$@ $? $(DLIBS) -lmkl_em64t
	$(RM) $?

clean:
	$(RM)  *.o
	$(RM)  *.mod
	$(RM) $(FPLS)

cart_to_polar.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

dogleg3.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

first_deriv_1.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

fourier_filter.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

func1.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

hmifits.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

multfactor.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

ola_xy_v10.o:
	$(FC) $(FFLAGS) $(@:.o=.f)

polyval.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

read_guess.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

ring_pass.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

ringanalysis.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

second_deriv.o:
	$(FC) $(F90F) -o $(@) $(@:.o=.f90)

.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL)

.f.o:
	$(FC) $(F90F) $(@:.o=.f)
