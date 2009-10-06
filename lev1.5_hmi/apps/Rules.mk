# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
GSLEXE_$(d)	:= $(addprefix $(d)/, HMI_observables)
GSLOBJ_$(d)	:= $(GSLEXE_$(d):%=%.o) 

MODEXE_$(d)	:= $(addprefix $(d)/,) 
MODEXE		:= $(MODEXE) $(MODEXE_$(d)) $(GSLEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

EXE_$(d)	:= $(MODEXE_$(d)) 
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 

DEP_$(d)	:= $(OBJ_$(d):%=%.o.d)
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(GSLOBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(GSLEXE_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)): CF_TGT := $(CF_TGT) -I $(SRCDIR)/$(d)/../libs/interp
#<<<<<<< Rules.mk
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) 
$(GSLEXE_$(d)):      LF_TGT := $(LF_TGT) 
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) $(GSLH)
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d) -I $(SRCDIR)/$(d)/../libs/interp
$(GSLEXE_$(d)):      LL_TGT := $(LL_TGT) $(GSLLIBS) -lmkl_em64t 

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBINTERP)
$(GSLEXE_$(d)): $(LIBINTERP)

#=======
#$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""
#>>>>>>> 1.4

# NOTE: Add dependent libraries with the -I compiler flag, and make the module depend
#   on that library
# $(OBJ_$(d)):				CF_TGT := -I$(SRCDIR)/$(d)/../../libs/somelib
# $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBSOMELIB)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
