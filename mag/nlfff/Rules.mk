# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

ifeq ($(JSOC_MACHINE), linux_avx)
# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.
dir     := $(d)/apps
-include                $(SRCDIR)/$(dir)/Rules.mk
endif

# Standard things
d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
