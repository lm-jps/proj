JPLAT=linux_x86_64
include Makevars/Makevars_$(JPLAT).mk
ICC =           icc 
CC =            $(ICC)
#CFLAGS =        -O3 -std=c99 -Wall -g
CFLAGS =        -O3 -std=c99 -Wall
CCFLAGS =       -c $(CFLAGS)

#DRMS=/home/production/cvs/JSOC
DRMS=/home/jsoc/cvs/JSOC
#DRMS=/home/jsoc

#DEST =	~/bin/$(JPLAT)
DEST =	.

INCL = -I$(DRMS)/base/include -I$(JSOCROOT)/base/include -I/home/jsoc/include  
LIBD = -L$(DRMS)/lib/$(JPLAT) -L/home/jsoc/lib/$(JPLAT) -L/home/production/cvs/JSOC/lib_third_party/lib/$(JPLAT)
DLIBS = $(LIBD) -ldrms -lpq -lcfitsio -lz -lpthread -lm -lg2c -lgsl -lgslcblas

FC = f77
F90F = -c

OBJ =	lfwrp.o do_one_limbfit.o limbfit.o limb.o expmax.o expfit.o nrutil.o indexx.o sort.o


# DRMS module with FORTRAN routine
MODULE =	lfwrp

all:	$(MODULE)

$(MODULE):	%:	%.o $(OBJ)
	$(ICC) -mcmodel=medium -xW -openmp -o $(DEST)/$@ $@.o $(OBJ) $(DLIBS)
	$(RM) $@.o

clean:
	$(RM) *.o
	$(RM) *.mod
	rm $(MODULE)

.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL) 

.f.o:
	$(FC) $(F90F) $(@:.o=.f)

limb.o: limb.f
	$(FC) $(F90F) limb.f

expfit.o: expfit.c
	$(CC) $(CCFLAGS) expfit.c $(INCL) 

expmax.o: expmax.c
	$(CC) $(CCFLAGS) expmax.c $(INCL)  

nrutil.o: nrutil.c
	$(CC) $(CCFLAGS) nrutil.c $(INCL)

limbfit.o: limbfit.c sort.c indexx.c expmax.o expfit.o nrutil.o 
	$(CC) $(CCFLAGS) limbfit.c $(INCL)  -DHMI

do_one_limbfit.o: do_one_limbfit.c
	$(CC) $(CCFLAGS) do_one_limbfit.c $(INCL) 

lfwrp.o: lfwrp.c do_one_limbfit.o limbfit.o 
	$(CC) $(CCFLAGS) lfwrp.c $(INCL) 
