# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Subdirectories. Directory-specific rules are optional here. The
# order NOT matter.
#ifeq ($(JSOC_MACHINE), linux_x86_64) 
	dir	:= $(d)/apps
	-include		$(SRCDIR)/$(dir)/Rules.mk
#endif

# Standard things
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
