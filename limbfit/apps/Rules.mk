# Standard things
sp 			:= $(sp).x
dirstack_$(sp)		:= $(d)
d			:= $(dir)

# Local variables
MODEXE_USEF_$(d)	:= $(addprefix $(d)/, lfwrp)
SUPPOBJ_$(d)		:= $(addprefix $(d)/, limbfit do_one_limbfit limb expmax expfit nrutil indexx sort)

MODEXE_USEF 		:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))
#MODEXE_USEF_SOCK_$(d)	:= $(MODEXE_USEF_$(d):%=%_sock)
#MODEXE_USEF_SOCK	:= $(MODEXE_USEF_SOCK) $(MODEXE_USEF_SOCK_$(d))
SUPPOBJ_$(d)		:= $(SUPPOBJ_$(d):%=%.o)

OBJ_$(d)		:= $(MODEXE_USEF_$(d):%=%.o) $(SUPPOBJ_$(d))
OBJ_$(d) :	 	CF_TGT := $(CF_TGT) -O3 -std=c99 -Wall

DEP_$(d)		:= $(OBJ_$(d):%=%.d)
CLEAN			:= $(CLEAN) \
			   $(OBJ_$(d)) \
			   $(MODEXE_USEF_$(d)) \
#			   $(MODEXE_USEF_SOCK_$(d)) \
			   $(DEP_$(d))

#TGT_BIN	        := $(TGT_BIN) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
TGT_BIN	        := $(TGT_BIN) $(MODEXE_USEF_$(d)) 

#S_$(d)			:= $(notdir $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d)))
S_$(d)			:= $(notdir $(MODEXE_USEF_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" $(GSLH) 


$(MODEXE_USEF_$(d)):		LL_TGT := $(LL_TGT) $(GSLL) -lgsl -lgslcblas

$(MODEXE_USEF_$(d)):		$(SUPPOBJ_$(d))

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
