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

# Subdirectories. Order matters - put libs first since apps will refer to libs.
ifdef USE_RPC
	dir	:= $(d)/libs
	-include		$(SRCDIR)/$(dir)/Rules.mk
	dir	:= $(d)/apps
	-include		$(SRCDIR)/$(dir)/Rules.mk
endif

# Standard things
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
