# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables

# Common utilities
#EXTRADEPS_$(d)		:= $(addprefix $(d)/, fresize.o)

# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/, sharp update_sharp_keys smarp)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# Modules with external libraries
MODEXE_USEF_$(d)	:= $(addprefix $(d)/, )
MODEXE_USEF		:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

MODEXE_USEF_SOCK_$(d)	:= $(MODEXE_USEF_$(d):%=%_sock)
MODEXE_USEF_SOCK	:= $(MODEXE_USEF_SOCK) $(MODEXE_USEF_SOCK_$(d))

EXE_$(d)	:= $(MODEXE_$(d)) $(MODEXE_USEF_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(MODEXE_USEF_SOCK_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/../../libs/interpolate -I$(SRCDIR)/$(d)/../../libs/stats -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) -I$(SRCDIR)/lib_third_party/include $(GSLH) $(FFTWH)
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

#$(EXTRADEPS_$(d)):	CF_TGT := $(CF_TGT) $(GSLH)

ifeq ($(JSOC_MACHINE), linux_avx2)
MKL     := -lmkl_rt
else
MKL     := -lmkl_em64t
endif

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
#$(ALL_$(d)) : $(EXTRADEPS_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBSTATS) $(LIBINTERP)
$(ALL_$(d)) : LF_TGT := $(LF_TGT) $(MKL)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(GSLLIBS) $(CFITSIOLIBS) $(MKL)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
ifneq ($(DEP_$(d)),)
  -include	$(DEP_$(d))
endif

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
