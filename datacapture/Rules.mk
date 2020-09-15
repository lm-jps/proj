# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

USE_RPC :=

ifneq ($(SUMS_USEMTSUMS),1)
	USE_RPC := yes
endif

ifneq ($(SUMS_USEMTSUMS_ALL),1)
	USE_RPC := yes
endif

# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.
ifdef USE_RPC
	dir	:= $(d)/apps
	-include		$(SRCDIR)/$(dir)/Rules.mk
endif

# Standard things
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
