# Standard things                                                                                                              
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables                                                                                                              
LIBDSDSMIGR     := $(d)/libdsdsmigr.a

OBJ_$(d)        := $(addprefix $(d)/, dsdsmigr.o)

DEP_$(d)        := $(OBJ_$(d):%=%.d)

CLEAN           := $(CLEAN) $(OBJ_$(d)) $(LIBDSDSMIGR) $(DEP_$(d))

TGT_LIB         := $(TGT_LIB) $(LIBDSDSMIGR)

S_$(d)          := $(notdir $(LIBDSDSMIGR))


# Local rules                                                                                                                  
$(OBJ_$(d)):    $(SRCDIR)/$(d)/Rules.mk

$(LIBDSDSMIGR):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts                                                                                                                    
.PHONY: $(S_$(d))
$(S_$(d)):      %:      $(d)/%

# Standard things
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
