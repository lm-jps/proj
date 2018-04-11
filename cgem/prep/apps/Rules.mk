# Standard things
sp := $(sp).x
dirstack_$(sp) := $(d)
d := $(dir)

## never touch above ---------------------------------------------------------------

## edit below -----------------------------------------------------

# module name(s) to be generated and file name(s) to be compiled

# Common utilities
EXTRADEPS_$(d)		:= $(addprefix $(d)/, flctsubs.o)

## C-wrapper name (name must end with .c)
MODEXE_USEF_$(d) := $(addprefix $(d)/, cgem_prep cgem_cutout cgem_map cgem_flct cgem_doppcal)
MODEXE_USEF := $(MODEXE_USEF) $(MODEXE_USEF_$(d))

## wrapped Fortran codes
WRAPPEDF_OBJ_$(d)    := $(addprefix $(d)/, azim_mindiff_jsoc.o dilate.o get_pils_los.o get_pils_rad.o median_sub.o pix2helio.o doppcal_estimate.o )

# flags for compiling and linking
MYCMPFLG_$(d) := -O2 -nofor-main -fp-model strict -fomit-frame-pointer
MYLNKFLG_$(d) := -lmkl_em64t -openmp -lfftw3 -lm

## edit above -----------------------------------------------------

## may not edit below -----------------------------------------------------

## compile and link etc.

# list .o .o.d and executables to be generated
#

OBJ_$(d) := $(MODEXE_USEF_$(d):%=%.o) $(WRAPPEDF_OBJ_$(d))
DEP_$(d) := $(OBJ_$(d):%=%.d)
CLEAN := $(CLEAN) \
         $(OBJ_$(d)) \
         $(MODEXE_USEF_$(d)) \
         $(DEP_$(d)) \
         $(EXTRADEPS_$(d))


# name(?) of target executable binary
TGT_BIN := $(TGT_BIN) $(MODEXE_USEF_$(d))

S_$(d) := $(notdir $(MODEXE_USEF_$(d)))

# making object/executable with Rules file and extra options
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -fp-model strict -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/../../../libs/interpolate -I/home/jsoc/lib/$(JSOC_MACHINE) $(FFTWH) $(FFTW3LIBS) $(FMATHLIBSH)

$(EXTRADEPS_$(d)):	CF_TGT := $(CF_TGT) $(FFTWH) $(GSLH) -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/src/
$(MODEXE_$(d)):			$(EXTRADEPS_$(d))
$(MODEXE_USEF_$(d)):    $(EXTRADEPS_$(d))

$(WRAPPEDF_OBJ_$(d)): FF_TGT := $(FF_TGT) $(MYCMPFLG_$(d)) -free 
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(MYLNKFLG_$(d)) -free
$(MODEXE_USEF_$(d)): $(WRAPPEDF_OBJ_$(d))

ALL_$(d)	:= $(MODEXE_USEF_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBSTATS) $(LIBINTERP)# $(LIBFFTW3)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(FFTW3LIBS) $(CFITSIOLIBS) $(MKL)

## never touch below ---------------------------------------------------------------------------

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)): %: $(d)/%

# Standard things
-include $(DEP_$(d))

d := $(dirstack_$(sp))
sp := $(basename $(sp))

# end of this file ---------------------------------------------------------------------------------
