# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
GSLEXE_$(d)	:= $(addprefix $(d)/, HMI_observables lookup phasemaps HMIfilters Leka HMI_IQUV_averaging smoothing HMI_Simulate_Doppler phasemaps_old HMI_observables2 HMI_IQUV_averaging2 phasemaps2 phasemaps_test phasemaps_FeI phasemaps_test_doublegaussian phasemaps_test_voigt HMI_observables_arta correction_velocities)
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
$(OBJ_$(d)):           CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../libs/lev15 -I/home/jsoc/include/ -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) -openmp -no-ipo
$(GSLEXE_$(d)):        LF_TGT := $(LF_TGT) -openmp -no-ipo 
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) $(GSLH)
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d) -I$(SRCDIR)/$(d)/../libs/lev15 -I/home/jsoc/include/ -I$(SRCDIR)/$(d)/../../libs/interpolate -DCDIR="\"$(SRCDIR)/$(d)\""
$(GSLEXE_$(d)):        LL_TGT := $(LL_TGT) $(GSLLIBS) -lpthread -lgsl
$(GSLEXE_$(d)):        LL_TGT := $(LL_TGT) $(GSLLIBS) -L$(FFTW_LIBS) -l$(FFTW3_LIB) -L/home/jsoc/lib/linux-x86_64 -lfftw3 -L$(JSOCROOT)/_linux_x86_64/proj/libs/interpolate -linterp  -static-intel -lmkl_em64t

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBLEV15) $(LIBINTERP)
$(GSLEXE_$(d)): $(LIBLEV15)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
