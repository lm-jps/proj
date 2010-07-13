# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables

# two libraries
LIBmextools	:= $(d)/libmextools.a
LIBmexrng	:= $(d)/libmexrng.a
# their composition
LIBmextoolsOBJ_$(d)	:= $(addprefix $(d)/, mextool.o mexargcheck.o num_strings.o ieee_consts.o)
LIBmexrngOBJ_$(d)	:= $(addprefix $(d)/, rng.o)

OBJ_$(d)	:= $(LIBmextoolsOBJ_$(d)) $(LIBmexrngOBJ_$(d))

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBmextools) $(LIBmexrng) $(DEP_$(d))

S_$(d)		:= $(notdir $(LIBmextools) $(LIBmexrng))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../mex2c

$(LIBmextools):	$(LIBmextoolsOBJ_$(d))
		$(ARCHIVE)

$(LIBmexrng):	$(LIBmexrngOBJ_$(d))
		$(ARCHIVE)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
