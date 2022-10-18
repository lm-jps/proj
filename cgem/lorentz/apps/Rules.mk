# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/, lorentz fix_lorentz)
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

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

TGT_BIN         := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) $(FFTWH) $(GSLH)
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/../../../libs/interpolate -I$(SRCDIR)/$(d)/src/

ifeq ($(JSOC_MACHINE), linux_avx2)
MKL     := -lmkl_rt
else
MKL     := -lmkl_em64t
endif

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
#$(ALL_$(d)) : $(EXTRADEPS_$(d))

# do not use $(LIBASTRO) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
# do not use $(LIBSTATS) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
# do not use $(LIBINTERP) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(ALL_$(d)) : proj/libs/astro/libastro.a proj/libs/stats/libstats.a proj/libs/interpolate/libinterp.a
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(GSLLIBS) $(MKL)

# NOTE: Add dependent libraries with the -I compiler flag, and make the module depend
#   on that library
# $(OBJ_$(d)):				CF_TGT := -I$(SRCDIR)/$(d)/../../../libs/astro -DCDIR="\"$(SRCDIR)/$(d)\""
ifeq ($(JSOC_MACHINE), linux_avx2)
 $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT :=  $(LL_TGT) $(CFITSIOLIBS)
else
 $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT :=  $(LL_TGT) $(FFTW3LIBS) $(CFITSIOLIBS)
endif

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
