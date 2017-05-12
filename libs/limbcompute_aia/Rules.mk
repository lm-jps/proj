# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables
LIBCOMPUTEAIA   := $(d)/liblc_aia.a

OBJ_$(d)        := $(addprefix $(d)/, limbcompute_aia.o)

DEP_$(d)        := $(OBJ_$(d):%=%.d)

CLEAN           := $(CLEAN) $(OBJ_$(d)) $(LIBCOMPUTEAIA) $(DEP_$(d))

TGT_LIB         := $(TGT_LIB) $(LIBCOMPUTEAIA)

S_$(d)          := $(notdir $(LIBCOMPUTEAIA))

# Local rules
$(OBJ_$(d)):    $(SRCDIR)/$(d)/Rules.mk

$(LIBCOMPUTEAIA):	$(OBJ_$(d))
			$(ARCHIVE)
			$(SLLIB)

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)):      %:      $(d)/%

# Standard things
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
