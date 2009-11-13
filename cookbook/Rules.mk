# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/, smpl_00 smpl_01 smpl_02 smpl_03 smpl_04 smpl_05 smpl_06 smpl_07 smpl_08)
MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)

EXE_$(d)        := $(MODEXE_$(d)) 
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

$(MODEXE_$(d)):      	LL_TGT := $(LL_TGT) $(PGLIBS) $(CFITSIOLIBS)

$(MODEXE_$(d)):	%:	%.o $(LIBDRMS_META)
			$(LINK)
			$(SLBIN)

$(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) $(PGLIBS) $(CFITSIOLIBS)

$(MODEXE_SOCK_$(d)):	%_sock:	%.o $(LIBDRMS_META_SOCK)
			$(LINK)
			$(SLBIN)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
