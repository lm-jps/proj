# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables
LIBEGSEHMICOMP	:= $(d)/libegsehmicomp.a

OBJ_$(d)	:= $(addprefix $(d)/, egsehmicomp.o)

DEP_$(d)        := $(OBJ_$(d):%=%.d)

CLEAN           := $(CLEAN) $(OBJ_$(d)) $(LIBEGSEHMICOMP) $(DEP_$(d))

TGT_LIB         := $(TGT_LIB) $(LIBEGSEHMICOMP)

S_$(d)          := $(notdir $(LIBEGSEHMICOMP))

# Local rules
$(OBJ_$(d)):    $(SRCDIR)/$(d)/Rules.mk

$(LIBEGSEHMICOMP):	$(OBJ_$(d))
	                $(ARCHIVE)
			$(SLLIB)

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)):      %:      $(d)/%

# Standard things
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
