# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/,  module_flatfield_combine module_flatfield write_flatfield write_offpoint write_badpix write_dark cosmic_ray_post)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

EXE_$(d)	:= $(MODEXE_$(d)) 
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(MODEXE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I $(SRCDIR)/$(d)/../libs/flatfieldlib 

# added 10/10
$(OBJ_$(d)):				CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(FFTW_INCS) 
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) $(FMATHLIBSL) -openmp -L$(FFTW_LIBS) -l$(FFTW3_LIB) $(MKL) 



$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate/
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBINTERP)
#added 10.10


# NOTE: Add dependent libraries with the -I compiler flag, and make the module depend
#   on that library
# $(OBJ_$(d)):				CF_TGT := -I$(SRCDIR)/$(d)/../../libs/somelib
# $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBSOMELIB)

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) $(FMATHLIBSL)  -lfftw3 -openmp 
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBFLATFIELD)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
