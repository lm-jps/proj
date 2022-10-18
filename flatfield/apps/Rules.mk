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
ifeq ($(JSOC_MACHINE), linux_avx2)
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) -lmkl_rt
else
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) $(FFTW3FLIBS) -lmkl_em64t
endif

$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate/

# do not use $(LIBINTERP) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	proj/libs/interpolate/libinterp.a
#added 10.10


# NOTE: Add dependent libraries with the -I compiler flag, and make the module depend
#   on that library
# $(OBJ_$(d)):				CF_TGT := -I$(SRCDIR)/$(d)/../../libs/somelib
# $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBSOMELIB)

# do not use $(LIBFLATFIELD) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	proj/flatfield/libs/flatfieldlib/libflatfieldlib.a

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
