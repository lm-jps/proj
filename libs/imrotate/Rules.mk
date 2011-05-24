# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
LIBIMROTATE	:= $(d)/libimrotate.a

OBJ_$(d)	:= $(addprefix $(d)/, imageRotate.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBIMROTATE) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBIMROTATE)

S_$(d)		:= $(notdir $(LIBIMROTATE))


# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(D)):	-I$(SRCDIR)/$(d)/

$(LIBIMROTATE):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
