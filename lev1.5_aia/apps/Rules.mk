# Standard things
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)


MODEXE_$(d)	:= $(addprefix $(d)/, aia_lev1p5 aia_synoptic aia_synoptic_nrt aia_most_recent)
MODEXE	 	:= $(MODEXE) $(MODEXE_$(d))
MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock) 
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

OBJ_$(d)	:= $(MODEXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(MODEXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))



$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/imrotate/ -I$(SRCDIR)/$(d)/../../libs/astro

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	proj/libs/astro/libastro.a proj/libs/imrotate/libimrotate.a


TGT_BIN	        := $(TGT_BIN) $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) 

S_$(d)		:= $(notdir $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)))

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)):      %:      $(d)/%

# Standard things
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
