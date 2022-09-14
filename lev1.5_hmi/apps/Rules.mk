# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)

USE_RPC :=

ifneq ($(SUMS_USEMTSUMS),1)
	USE_RPC := yes
endif

ifneq ($(SUMS_USEMTSUMS_ALL),1)
	USE_RPC := yes
endif

DEFMAKE_$(d)	:= $(addprefix $(d)/, HMI_lookup phasemaps_voigt HMI_observables_90q HMI_observables HMI_IQUV_averaging correction_velocities ingest_dcon ingest_dcon_gen HMI_observables_dcon HMI_observables_dcon2 HMI_observables_dconS HMI_IQUV_averaging_dcon HMI_observables_135q stokes_dcon)

ifdef USE_RPC
	DEFMAKE_RPC_$(d)	:= $(addprefix $(d)/, lev1_dcon cont_dcon)
	DEFMAKE_$(d)	:= $(DEFMAKE_$(d)) $(DEFMAKE_RPC_$(d))
endif

GSLEXE_$(d)	:= $(addprefix $(d)/, HMI_observables2 undistort_lev1) $(DEFMAKE_$(d))

# THERE IS AN HMI_IQUV_averaging2.c IN THE SOURCE DIRECTORY, BUT THERE IS NO RULE TO MAKE THIS FILE.

GSLOBJ_$(d)	:= $(GSLEXE_$(d):%=%.o)

MODEXE		:= $(MODEXE) $(GSLEXE_$(d))

GSLEXE_SOCK_$(d):= $(GSLEXE_$(d):%=%_sock)
MODEXE_SOCK_$(d):= $(GSLEXE_SOCK_$(d))
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

EXE_$(d)	:= $(GSLEXE_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o)

DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(DEFMAKE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):           CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../libs/lev15 -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate -I$(SRCDIR)/$(d)/../../lev0/apps
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) $(GSLH) $(FFTWH)
$(GSLEXE_$(d)) $(GSLEXE_SOCK_$(d)):        LL_TGT := $(LL_TGT) $(GSLLIBS) $(FFTW3LIBS) -lmkl_em64t

# do not use $(LIBLEV15) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
# do not use $(LIBLIMBFITFXN) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
# do not use $(LIBINTERP) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(EXE_$(d)) $(MODEXE_SOCK_$(d)):	proj/lev1.5_hmi/libs/lev15/liblev15.a proj/lev0/apps/liblimbfitfxn.a proj/libs/interpolate/libinterp.a


# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
