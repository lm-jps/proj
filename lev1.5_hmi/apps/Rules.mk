# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
GSLEXE_$(d)	:= $(addprefix $(d)/, HMI_observables lookup phasemaps HMIfilters Leka HMI_IQUV_averaging smoothing HMI_Simulate_Doppler phasemaps_old HMI_observables2 HMI_IQUV_averaging2 phasemaps2 phasemaps_test phasemaps_FeI phasemaps_test_doublegaussian phasemaps_test_voigt HMI_observables_arta correction_velocities)
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

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):           CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../libs/lev15 -I/home/jsoc/include/ -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) -openmp -no-ipo
$(GSLEXE_$(d)) $(GSLEXE_SOCK_$(d)):        LF_TGT := $(LF_TGT) -openmp -no-ipo 
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) $(GSLH)
$(GSLEXE_$(d)) $(GSLEXE_SOCK_$(d)):        LL_TGT := $(LL_TGT) $(GSLLIBS) -lpthread -lgsl
$(GSLEXE_$(d)) $(GSLEXE_SOCK_$(d)):        LL_TGT := $(LL_TGT) $(GSLLIBS) -L$(FFTW_LIBS) -l$(FFTW3_LIB) -L/home/jsoc/lib/linux-x86_64 -lfftw3 -static-intel -lmkl_em64t

$(EXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBLEV15) $(LIBINTERP)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
