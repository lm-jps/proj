# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

LIBFLATFIELD	:= $(d)/libflatfieldlib.a


OBJ_$(d)	:= $(addprefix $(d)/, module_flatfield_code.o module_flatfield_read.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBFLATFIELD) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBFLATFIELD)

S_$(d)		:= $(notdir $(LIBFLATFIELD))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk

$(LIBFLATFIELD):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
