# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.
# mag team probably isn't going to use apps dir
# dir     := $(d)/apps
# -include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/lorentz
-include                $(SRCDIR)/$(dir)/Rules.mk

# Standard things
d               := $(dirstack_$(sp))
sp              := $(basename $(sp))

