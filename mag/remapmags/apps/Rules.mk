# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables

# Common utilities
EXTRADEPS_$(d)		:= $(addprefix $(d)/, )

# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/, fdlos2radial resizemappingmag)
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

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/../../../libs/interpolate -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) $(FFTWH)
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

$(EXTRADEPS_$(d)):	CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d))
$(ALL_$(d)) : $(EXTRADEPS_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBSTATS) $(LIBINTERP)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(FFTW3LIBS) $(CFITSIOLIBS) -lmkl_em64t

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
ifneq ($(DEP_$(d)),)
  -include	$(DEP_$(d))
endif

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
