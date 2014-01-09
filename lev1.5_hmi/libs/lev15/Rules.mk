# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

LIBLEV15	:= $(d)/liblev15.a


OBJ_$(d)	:= $(addprefix $(d)/, interpol_code.o polcal.o Dopplergram.o Dopplergram2.o fstats.o fstats2.o dstats.o dstats2.o Dopplergram_test.o Dopplergram_test2.o Dopplergram_duvall.o Dopplergram_Jesper.o Dopplergram_largercrop.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBLEV15) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBLEV15)

S_$(d)		:= $(notdir $(LIBLEV15))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):    CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../../libs/interpolate -I$(FFTW_INCS)
$(OBJ_$(d)):    LL_TGT := $(LL_TGT) -L$(SRCDIR)/$(d)/../../../libs/interpolate -linterp
$(LIBLEV15):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
