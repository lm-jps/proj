# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
DEFMAKE_$(d)	:= $(addprefix $(d)/, HMI_observables HMI_IQUV_averaging correction_velocities)

GSLEXE_$(d)	:= $(addprefix $(d)/, lookup phasemaps HMIfilters Leka smoothing HMI_Simulate_Doppler phasemaps_old HMI_observables_test2 HMI_IQUV_averaging_test2 phasemaps2 phasemaps_test phasemaps_FeI phasemaps_test_doublegaussian phasemaps_test_voigt ingest_lookup ingest_core_intensity lookup_test phasemaps_test_voigt_obsmode phasemaps_test_RogerKittPeak phasemaps_test_voigt_Iripple ingest_Thomas phasemaps_test_voigt_Iripple2 phasemaps_test_voigt_Iripple3 lookup_Iripple ingest_Fleck phasemaps_test_voigt_Iripple2_all ingest_corrected_phasemaps ingest_Yang HMI_observables_test HMI_observables2 HMI_IQUV_averaging_test phasemaps_test_voigt_obsmode_Iripple_all lookup_test2 ingest_v50 undistort_lev1 HMI_observables_duvall ingest_KokLengYeo ingest_FleckS ingest_FleckN ingest_FleckSN limbfit_sc time_convert ingest_Aimee) $(DEFMAKE_$(d))

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
$(OBJ_$(d)):           CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../libs/lev15 -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate
$(GSLOBJ_$(d)):        CF_TGT := $(CF_TGT) $(GSLH) $(FFTWH)
$(GSLEXE_$(d)) $(GSLEXE_SOCK_$(d)):        LL_TGT := $(LL_TGT) $(GSLLIBS) $(FFTW3LIBS) -lmkl_em64t

$(EXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBLEV15) $(LIBINTERP)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
