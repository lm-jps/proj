# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
LIBMAGUTILS	:= $(d)/libmagutils.a

OBJ_$(d)	:= $(addprefix $(d)/, img2helioVector.o noisemaskimag4twindow.o obstime2maskid.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBMAGUTILS) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBMAGUTILS)

S_$(d)		:= $(notdir $(LIBMAGUTILS))


# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk

$(LIBMAGUTILS):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
