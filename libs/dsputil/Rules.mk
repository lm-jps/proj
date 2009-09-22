# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
LIBDSPUTIL	:= $(d)/libdsputil.a
OBJ_$(d)	:= $(addprefix $(d)/, bspline_interp.o bspline_coeff.o interpolate.o fftw_convolve.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBDSPUTIL) $(DEP_$(d))

S_$(d)		:= $(notdir $(LIBDSPUTIL))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk

ifeq ($(JSOC_MACHINE), linux_ia32)
  FFTW_$(d) = /home/jsoc/lib/linux-ia32
endif
ifeq ($(JSOC_MACHINE), linux_x86_64)
  FFTW_$(d) = /home/jsoc/lib/linux-x86_64
endif

$(OBJ_$(d)):	CF_TGT := -I/home/jsoc/include

$(LIBDSPUTIL):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
