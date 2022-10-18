# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

ITH_$(d)        := $(addprefix $(d)/, invert_td_hr)
TT_$(d)         := $(addprefix $(d)/, travel_times)

MODEXE_USEF_$(d)     := $(ITH_$(d)) $(TT_$(d))
MODEXE_USEF          := $(MODEXE_USEF) $(MODEXE_USEF_$(d))

EXE_$(d)        := $(MODEXE_USEF_$(d))

OBJ_ITH_$(d)    := $(addprefix $(d)/, mcd_cs_hmi_sub.o mcd_v_hmi_sub.o)
OBJ_TT_$(d)     := $(addprefix $(d)/, cal_HMI_2x2ave_sub.o cal_HMI_noave_sub.o sub_correlate_BLAS.o sub_filter_HMI_ppline.o sub_lmfit.o sub_shift.o sub_GB_fitting_2002.o sub_do_fitting.o)

OBJ_$(d)        := $(EXE_$(d):%=%.o) $(OBJ_ITH_$(d)) $(OBJ_TT_$(d))
DEP_$(d)        := $(OBJ_$(d):%=%.d)
CLEAN           := $(CLEAN) $(OBJ_$(d)) $(EXE_$(d)) $(DEP_$(d))

TGT_BIN         := $(TGT_BIN) $(EXE_$(d))

S_$(d)          := $(notdir $(EXE_$(d)))

# Local rules
$(OBJ_$(d)):    $(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):    CF_TGT := $(CF_TGT) $(FFTWH)

ifeq ($(JSOC_MACHINE), linux_avx2)
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) -lmkl_rt
else
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(FFTW3LIBS) $(FFTW3FLIBS) -mkl
endif

$(ITH_$(d)):    $(OBJ_ITH_$(d))
$(TT_$(d)):     $(OBJ_TT_$(d))


# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)):      %:      $(d)/%

# Standard things
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))


