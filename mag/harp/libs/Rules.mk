# Standard things                                                                                             
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

dir     := $(d)/matlab
-include                $(SRCDIR)/$(dir)/Rules.mk

# Standard things                                                                                            
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
