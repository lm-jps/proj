JPLAT=linux_x86_64
ICC =           icc11
CC =            $(ICC)
CFLAGS =        -O3 -std=c99 -D_GNU_SOURCE
CCFLAGS =       -c $(CFLAGS)

# added by K.H. ... for some reasons, we need Rick's file...
DRMS=/home/rick/drms


#DEST =	~/bin/$(JPLAT)
DEST =	.

LDCMD = ifort -nofor-main -openmp
MPLDCMD = /home/jsoc/mpich2/bin/mpif90 -nofor-main -openmp
INCL = -I$(DRMS)/base/include -I/home/jsoc/mpich2/include -I$(JSOCROOT)/base/include
LIBD = -L$(DRMS)/lib/$(JPLAT) -L/home/jsoc/lib/$(JPLAT)  -L/home/jsoc/mpich2/lib  -L/home/production/cvs/JSOC/lib_third_party/lib/$(JPLAT)
DLIBS = $(LIBD) -ldrms -lpq -lcfitsio -lmkl_em64t -lz -lpthread -lm

FC =	ifort
MPCC =	/home/jsoc/mpich2/bin/mpicc
F90F = -c -O3
MPFC =	/home/jsoc/mpich2/bin/mpif90

OBJ =	cons_param.o line_param.o filt_param.o inv_param.o line_param.o \
	svd_param.o filt_init.o free_init.o inv_init.o line_init.o svd_init.o \
	voigt_init.o wave_init.o voigt_data.o voigt_taylor.o free_memory.o \
	inv_utils.o voigt.o svdcmp.o svbksb.o forward.o invert.o wfa_guess.o

# DRMS module with FORTRAN routines and MPI
MODULE =	vfisv

all:	$(MODULE)

$(MODULE):	%:	%.o $(OBJ)
	$(MPLDCMD) -o $(DEST)/$@ $@.o $(OBJ) $(DLIBS)
	$(RM) $@.o

clean:
	$(RM) *.o
	$(RM) *.mod
	rm $(MODULE)

.c.o:
	$(MPCC) $(CCFLAGS) $(@:.o=.c) $(INCL)

.f.o:
	$(FC) $(F90F) $(@:.o=.f)

cons_param.o: cons_param.f90
	$(FC) $(F90F) cons_param.f90

filt_init.o: filt_init.f90 filt_param.o cons_param.o line_param.o
	$(FC) $(F90F) filt_init.f90

filt_param.o: filt_param.f90 line_param.o cons_param.o
	$(FC) $(F90F) filt_param.f90

forward.o: forward.f90 line_param.o cons_param.o filt_param.o inv_param.o\
	voigt.o voigt_taylor.o
	$(FC) $(F90F) forward.f90

free_init.o: free_init.f90 inv_param.o cons_param.o filt_param.o
	$(FC) $(F90F) free_init.f90

free_memory.o: free_memory.f90 filt_param.o svd_param.o line_param.o
	$(FC) $(F90F) free_memory.f90

inv_init.o: inv_init.f90 inv_param.f90
	$(FC) $(F90F) inv_init.f90

inv_param.o: inv_param.f90 cons_param.o
	$(FC) $(F90F) inv_param.f90

inv_utils.o: inv_utils.f90 cons_param.o filt_param.o \
	svd_param.o inv_param.o line_param.o svdcmp.o svbksb.o
	$(FC) $(F90F) inv_utils.f90

wfa_guess.o: wfa_guess.f90 cons_param.o filt_param.o inv_param.o
	$(FC) $(F90F) wfa_guess.f90

invert.o: invert.f90 forward.o line_param.o cons_param.o filt_param.o \
	inv_utils.o svd_param.o inv_param.o svdcmp.o svbksb.o wfa_guess.o
	$(FC) $(F90F) invert.f90

line_init.o: line_init.f90  line_param.o cons_param.o
	$(FC) $(F90F) line_init.f90

line_param.o: line_param.f90 cons_param.o
	$(FC) $(F90F) line_param.f90

svbksb.o: svbksb.f90 cons_param.o
	$(FC) $(F90F) svbksb.f90

svdcmp.o: svdcmp.f90 cons_param.o
	$(FC) $(F90F) svdcmp.f90

svd_init.o: svd_init.f90 svd_param.o inv_param.o cons_param.o
	$(FC) $(F90F) svd_init.f90

svd_param.o: svd_param.f90
	$(FC) $(F90F) svd_param.f90

#voigt.o: voigt.f
#	$(FC) $(F90F) voigt.f

voigt_data.o: voigt_data.f90 cons_param.o
	$(FC) $(F90F) voigt_data.f90

voigt_init.o: voigt_init.f90 voigt_data.o cons_param.o
	$(FC) $(F90F) voigt_init.f90

voigt_taylor.o: voigt_taylor.f90 voigt_data.o voigt_init.o \
	voigt.o cons_param.o line_param.o
	$(FC) $(F90F) voigt_taylor.f90

wave_init.o: wave_init.f90 line_param.o cons_param.o
	$(FC) $(F90F) wave_init.f90
