# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
RDVINV_$(d)	:= $(addprefix $(d)/, rdvinv)
RDFITF_$(d)	:= $(addprefix $(d)/, rdfitf)

ifeq ($(COMPILER), icc)
  MODEXE_$(d)	:= $(addprefix $(d)/, datavg maicalc maproj mtrack pspec3 rdcover rdfitc xtrackd) 
endif
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

MODEXE_USEF_$(d):= $(RDVINV_$(d)) $(RDFITF_$(d))
MODEXE_USEF	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

CEXE_$(d)	:= $(addprefix $(d)/, gentargs) 
CEXE		:= $(CEXE) $(CEXE_$(d))

EXE_$(d)        := $(MODEXE_$(d)) $(MODEXE_USEF_$(d) $(CEXE_$(d))
OBJ_F_$(d)      := $(addprefix $(d)/, cart_to_polar.o fourier_filter.o ring_pass.o polyval.o multfactor.o func1.o first_deriv_1.o second_deriv.o dogleg3.o hmifits.o ringanalysis.o ccint2.o fourierlibn32.o)

OBJ_$(d)	:= $(EXE_$(d):%=%.o) $(OBJ_F_$(d))
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
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) $(FFTWH)

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) $(FFTW3LIBS) $(FFTW3FLIBS) 

$(RDFITF_$(d)):				$(OBJ_F_$(d))
$(RDVINV_$(d)) $(RDFITF_$(d)):				LL_TGT := $(LL_TGT) -lmkl_em64t

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
