# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Subdirectories. Directory-specific rules are optional here. The
# order does NOT matter.
dir	:= $(d)/mexfunctions
-include		$(SRCDIR)/$(dir)/Rules.mk
dir	:= $(d)/util
-include		$(SRCDIR)/$(dir)/Rules.mk
dir	:= $(d)/mex2c
-include		$(SRCDIR)/$(dir)/Rules.mk
# more dirs just go like this:
#dir	:= $(d)/dr
#-include		$(SRCDIR)/$(dir)/Rules.mk

# Standard things
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))

