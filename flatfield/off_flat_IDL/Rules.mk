# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

CEXE_$(d)       := $(addprefix $(d)/, flatfield_iter_int)
CEXE            := $(CEXE) $(CEXE_$(d))

OBJ_$(d)	:= $(CEXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(CEXE_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(CEXE_$(d))

S_$(d)		:= $(notdir $(CEXE_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
