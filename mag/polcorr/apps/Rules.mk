# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

MODEXE_$(d)	:= $(addprefix $(d)/, update_poldb polcorr)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# flags for compiling and linking
MYCMPFLG	:= -openmp
MYLNKFLG	:= -lmkl_em64t -openmp -lcfitsio

EXE_$(d)	:= $(MODEXE_$(d)) 
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

TGT_BIN         := $(TGT_BIN) $(MODEXE_$(d)) $(MODEXE_SOCK_$(d))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) $(FFTWH) $(GSLH) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJ_$(d)):		CF_TGT := -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/../../../libs/interpolate -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) -I$(SRCDIR)/lib_third_party/include $(MYCMPFLG) -I/home/jsoc/include -I$(SRCDIR)/$(d)/src/ 

MKL     := -lmkl_em64t

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
#$(ALL_$(d)) : $(EXTRADEPS_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBSTATS) $(LIBINTERP)
$(ALL_$(d)) : LF_TGT := $(LF_TGT) $(MKL)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(GSLLIBS) $(CFITSIOLIBS) $(MKL)

# NOTE: Add dependent libraries with the -I compiler flag, and make the module depend
#   on that library
# $(OBJ_$(d)):				CF_TGT := -I$(SRCDIR)/$(d)/../../../libs/astro -DCDIR="\"$(SRCDIR)/$(d)\""
 $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT :=  $(LL_TGT) $(GSLLIBS) $(FFTW3LIBS) $(CFITSIOLIBS) -lmkl_em64t

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
