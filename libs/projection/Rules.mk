# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
LIBCARTOGRAPHY	:= $(d)/libcartography.a

OBJ_$(d)	:= $(addprefix $(d)/, cartography.o imginfo.o keystuff.o mdistuff.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBCARTOGRAPHY) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBCARTOGRAPHY)

S_$(d)		:= $(notdir $(LIBCARTOGRAPHY))


# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk

$(LIBCARTOGRAPHY):	$(OBJ_$(d))
			$(ARCHIVE)
			$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
