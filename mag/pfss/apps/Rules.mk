# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)

# Common utilities
EXTRADEPS_$(d)		:= $(addprefix $(d)/, pfss_pkg.o fieldline_pkg.o wsa_pkg.o hccsss_pkg.o)

MODEXE_$(d)	:= $(addprefix $(d)/, jgh2pfss jgh2wsa)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# Modules with external libraries
MODEXE_USEF_$(d)	:= $(addprefix $(d)/, jsynop2gh jgh2hccsss jpolfil)
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
		   $(DEP_$(d)) \
		   $(EXTRADEPS_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := -I$(SRCDIR)/proj/libs/astro $(FMATHLIBSH) -I/home/jsoc/include -I$(SRCDIR)/$(d)
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

$(EXTRADEPS_$(d)):	CF_TGT := $(CF_TGT) -I/home/jsoc/include -I$(SRCDIR)/proj/libs/astro -I$(SRCDIR)/$(d)
$(MODEXE_$(d)):			$(EXTRADEPS_$(d))
$(MODEXE_USEF_$(d)):    $(EXTRADEPS_$(d))


MKL     := -lmkl

ifeq ($(COMPILER), icc)
  NOIPO_$(d)    := -no-ipo
endif

ifeq ($(JSOC_MACHINE), linux_x86_64)
  MKL     := $(NOIPO_$(d)) -lmkl_em64t
endif
ifeq ($(JSOC_MACHINE), linux_ia32)
  MKL     := $(NOIPO_$(d)) -lmkl_lapack -lmkl_ia32
endif

ifeq ($(JSOC_MACHINE), linux_ia32)
  FFTW_$(d) = /home/jsoc/lib/linux-ia32
endif
ifeq ($(JSOC_MACHINE), linux_x86_64)
  FFTW_$(d) = /home/jsoc/lib/linux-x86_64
endif

SVML_$(d)       :=
GUIDE_$(d)      :=

ifeq ($(COMPILER), icc)
  SVML_$(d)     := -lsvml
  GUIDE_$(d)    := -lguide
endif

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
$(ALL_$(d)) : $(LIBASTRO)
$(ALL_$(d)) : LL_TGT :=  $(LL_TGT) $(GSLLIBS) $(CFITSIOLIBS) $(FMATHLIBS) -L$(FFTW_$(d)) -lfftw3f $(SVML_$(d)) $(MKL) $(GUIDE_$(d)) 

# Shortcuts
.PHONY:	$(S_$(d)) pfss
$(S_$(d)):	%:	$(d)/%
pfss : $(ALL_$(d))

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
