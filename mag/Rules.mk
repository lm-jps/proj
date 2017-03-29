# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.

dir     := $(d)/libs
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/pfss
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/ambig
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/harp
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/ident
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/patch
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/d4vm
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/nlfff
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/remapmags
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/synop
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/polarfield
-include                $(SRCDIR)/$(dir)/Rules.mk
dir     := $(d)/polcorr
-include                $(SRCDIR)/$(dir)/Rules.mk

# Standard things
d               := $(dirstack_$(sp))
sp              := $(basename $(sp))

