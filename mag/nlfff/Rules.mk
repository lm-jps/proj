# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

ifneq ($(JSOC_MACHINE), linux_x86_64)
# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.
dir     := $(d)/apps
-include                $(SRCDIR)/$(dir)/Rules.mk
endif

# Standard things
d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
