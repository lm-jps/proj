# Standard things
sp := $(sp).x
dirstack_$(sp) := $(d)
d := $(dir)

## never touch above ---------------------------------------------------------------

## edit below -----------------------------------------------------

# module name(s) to be generated and file name(s) to be compiled

## C-wrapper name (name must end with .c)
MODEXE_USEF_$(d) := $(addprefix $(d)/, cgem_pdfi)
MODEXE_USEF := $(MODEXE_USEF) $(MODEXE_USEF_$(d))

## wrapped Fortran codes
WRAPPEDF_OBJ_$(d)    := $(addprefix $(d)/, pdfi_wrapper4jsoc_test.o)

# flags for compiling and linking
MYCMPFLG_$(d) := -O2 -nofor-main -fp-model strict
MYLNKFLG_$(d) := -lmkl_em64t

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
         $(DEP_$(d))


# name(?) of target executable binary
TGT_BIN := $(TGT_BIN) $(MODEXE_USEF_$(d))

S_$(d) := $(notdir $(MODEXE_USEF_$(d)))

# making object/executable with Rules file and extra options
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/../../../libs/interpolate $(FFTWH) $(FFTW3LIBS) $(FMATHLIBSH) -I/home/jsoc/lib/$(JSOC_MACHINE)

$(WRAPPEDF_OBJ_$(d)): FF_TGT := $(FF_TGT) $(MYCMPFLG_$(d))
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(MYLNKFLG_$(d)) -lfishpack_r8
$(MODEXE_USEF_$(d)): $(WRAPPEDF_OBJ_$(d))

ALL_$(d)	:= $(MODEXE_USEF_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBSTATS) $(LIBINTERP)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(GSLLIBS)

## never touch below ---------------------------------------------------------------------------

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)): %: $(d)/%

# Standard things
-include $(DEP_$(d))

d := $(dirstack_$(sp))
sp := $(basename $(sp))

# end of this file ---------------------------------------------------------------------------------
