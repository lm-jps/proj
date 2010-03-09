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


# Add something like this when linking against this library
# ifeq ($(COMPILER), icc)
#   MKL     := -lmkl
#
#   ifeq ($(JSOC_MACHINE), linux_x86_64) 
#     MKL     := -static-intel -lmkl_em64t
#   endif
#
# MODPATH:	LL_TGT := $(LL_TGT) -L$(CFITSIO_LIBS) -l$(FFTW3_LIB)
#
# endif

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
