# For the interpolate library, the make system will generate 64-bit binaries only.

# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

LIBINTERP	:= $(d)/libinterp.a

OBJ_$(d)	:= $(addprefix $(d)/, cholesky_down.o finterpolate.o gapfill.o tinterpolate.o fresize.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBINTERP) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBINTERP)

S_$(d)		:= $(notdir $(LIBINTERP))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -I$(FFTW_INCS) -DCDIR="\"$(SRCDIR)/$(d)\""


# For modules linking against this library, please 
# consult proj/lev1.5_hmi/apps/Rules.mk


$(LIBINTERP):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
