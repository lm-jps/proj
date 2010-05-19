# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

LIBSTATS	:= $(d)/libstats.a

OBJ_$(d)	:= $(addprefix $(d)/, fstats.o fstats2.o dstats.o dstats2.o set_statistics.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBSTATS) $(DEP_$(d))

TGT_LIB		:= $(TGT_LIB) $(LIBSTATS)

S_$(d)		:= $(notdir $(LIBSTATS))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):    CF_TGT := $(CF_TGT)
$(OBJ_$(d)):    LL_TGT := $(LL_TGT) 
$(LIBSTATS):	$(OBJ_$(d))
		$(ARCHIVE)
		$(SLLIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
