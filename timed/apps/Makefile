#
#  Common Makefile requiring re-compilation
#

# MAKEROOT = $(DRMS)/proj/cookbook/Makevars
include $(MAKEROOT)/Makevars_$(JPLAT).mk

FC = ifort11
LDCMD = $(FC) -nofor_main -no-ipo 

INCL = -I$(DRMS)/base/include
LIBD = -L$(DRMS)/lib/$(JPLAT) -L/home/jsoc/lib/$(JPLAT)

DLIBS = $(LIBD) -ldrms -lpq -lz -lcfitsio -lfftw3f -lmkl_em64t -lguide -lpthread

DEST = ../bin/$(JPLAT)

NVOBJ =	mcd_cs_hmi_sub.o mcd_v_hmi_sub.o

TTOBJ =	cal_HMI_2x2ave_sub.o cal_HMI_noave_sub.o sub_correlate_BLAS.o \
	sub_filter_HMI_ppline.o sub_lmfit.o sub_shift.o sub_GB_fitting_2002.o \
	sub_do_fitting.o

FPLS =	$(addprefix $(DEST)/, invert_td_hr) $(addprefix $(DEST)/, travel_times)

all:	invert_td_hr travel_times

.SECONDEXPANSION:

invert_td_hr:	$$@.o $(NVOBJ)
	$(LDCMD) -o $(DEST)/$@ $? $(DLIBS)
	$(RM) $?

travel_times:	$$@.o $(TTOBJ)
	$(LDCMD) -o $(DEST)/$@ $? $(DLIBS)
	$(RM) $?

clean:
	$(RM)  *.o
	$(RM) $(FPLS)

.c.o:
	$(CC) $(CCFLAGS) $(@:.o=.c) $(INCL)

.f.o:
	$(FC) $(FFLAGS) $(@:.o=.f)
