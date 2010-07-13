# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
LIBmex2c_jsoc	:= $(d)/libmex2c_jsoc.a

OBJ_$(d)	:= $(addprefix $(d)/, mex_stubs.o mex2c_main_lib.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBmex2c_jsoc) $(DEP_$(d))

# note, the Rules.mk for .../hmi_lev1.5/libs has a TGT_LIB line...is
# this needed?  Would it clean up the build if it was there?

S_$(d)		:= $(notdir $(LIBmex2c_jsoc))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../util

$(LIBmex2c_jsoc):	$(OBJ_$(d))
		$(ARCHIVE)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
