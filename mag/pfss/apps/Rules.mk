# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)

# Common utilities
EXTRADEPS_$(d)		:= $(addprefix $(d)/, pfss_pkg.o fieldline_pkg.o wsa_pkg.o hccsss_pkg.o)

MODEXE_$(d)	:= $(addprefix $(d)/, jpolfil)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# Modules with external libraries
# MODEXE_USEF_$(d)	:= $(addprefix $(d)/,)
MODEXE_USEF_$(d)	:= 
# MODEXE_USEF		:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

#MODEXE_USEF_SOCK_$(d)	:= $(MODEXE_USEF_$(d):%=%_sock)
#MODEXE_USEF_SOCK	:= $(MODEXE_USEF_SOCK) $(MODEXE_USEF_SOCK_$(d))

EXE_$(d)	:= $(MODEXE_$(d)) $(MODEXE_USEF_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d)) \
		   $(EXTRADEPS_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d)) #$(MODEXE_USEF_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) $(FFTWH) $(GSLH) -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/src/ 
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

$(EXTRADEPS_$(d)):	CF_TGT := $(CF_TGT) $(FFTWH) $(GSLH) -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/src/
$(MODEXE_$(d)):			$(EXTRADEPS_$(d))
$(MODEXE_USEF_$(d)):    $(EXTRADEPS_$(d))

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
$(ALL_$(d)) : $(LIBASTRO)
ifeq ($(JSOC_MACHINE), linux_avx2)
$(ALL_$(d)) : LL_TGT :=  $(LL_TGT) $(GSLLIBS) $(CFITSIOLIBS) -lmkl_rt
else
$(ALL_$(d)) : LL_TGT :=  $(LL_TGT) $(GSLLIBS) $(FFTW3LIBS) $(CFITSIOLIBS) -lmkl_em64t
endif

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
