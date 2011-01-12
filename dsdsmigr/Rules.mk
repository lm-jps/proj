# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Subdirectories. Order matters - put libs first since apps will refer to libs.
dir	:= $(d)/libs
-include		$(SRCDIR)/$(dir)/Rules.mk
dir	:= $(d)/apps
-include		$(SRCDIR)/$(dir)/Rules.mk

# Standard things
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
