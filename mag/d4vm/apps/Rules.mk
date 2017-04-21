# Standard things
sp := $(sp).x
dirstack_$(sp) := $(d)
d := $(dir)

## never touch above ---------------------------------------------------------------

## edit below -----------------------------------------------------

# module name(s) to be generated and file name(s) to be compiled

## C-wrapper name (name must end with .c)
MODEXE_USEF_$(d) := $(addprefix $(d)/, test_d4vm dave4vm4velocity)

## wrapped Fortran codes
WRAPPEDF_$(d)    := $(addprefix $(d)/, d4vm_preproc.o d4vm_matrix.o d4vm_derivs.o d4vm_2dconv.o d4vm_solver.o d4vm_kernel.o)

# flags for compiling and linking
# MYCMPFLG := -O2 -free -nofor-main
MYLNKFLG := -lmkl_em64t_sequential -lmkl_core -lifport -lifcore -limf

## edit above -----------------------------------------------------

## may not edit below -----------------------------------------------------

## compile and link etc.

# list .o .o.d and executables to be generated
#
MODEXE_USEF := $(MODEXE_USEF) $(MODEXE_USEF_$(d))
#
OBJ_$(d) := $(MODEXE_USEF_$(d):%=%.o) $(WRAPPEDF_$(d))
#
DEP_$(d) := $(OBJ_$(d):%=%.d)

CLEAN := $(CLEAN) $(OBJ_$(d)) $(MODEXE_USEF_$(d)) $(DEP_$(d))

# name(?) of target executable binary
TGT_BIN := $(TGT_BIN) $(MODEXE_USEF_$(d))
# what is this ?
S_$(d) := $(notdir $(MODEXE_USEF_$(d)))

# making object/executable with Rules file and extra options
$(OBJ_$(d)): $(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)): CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../../libs/stats

$(WRAPPEDF_$(d)): FF_TGT := $(FF_TGT) # $(MYCMPFLG)
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(MYLNKFLG)
$(MODEXE_USEF_$(d)): $(WRAPPEDF_$(d))

## never touch below ---------------------------------------------------------------------------

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)): %: $(d)/%

# Standard things
-include $(DEP_$(d))

d := $(dirstack_$(sp))
sp := $(basename $(sp))

# end of this file ---------------------------------------------------------------------------------
