# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.
dir	:= $(d)/astro
-include		$(SRCDIR)/$(dir)/Rules.mk
dir	:= $(d)/dr
-include		$(SRCDIR)/$(dir)/Rules.mk
dir	:= $(d)/dsputil
-include		$(SRCDIR)/$(dir)/Rules.mk
dir	:= $(d)/gapfiller
-include		$(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/stats
-include                $(SRCDIR)/$(dir)/Rules.mk

ifeq ($(COMPILER), icc)
ifneq ($(JSOC_MACHINE), linux_ia32) 
    # It looks like our 32-bit machines do not have recent versions of icc and mkl
    # so don't build on 32-bit machines. Also, don't trust gcc to build the code
    # that links to mkl.
dir	:= $(d)/interpolate
-include		$(SRCDIR)/$(dir)/Rules.mk
endif
endif

# Standard things
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))

