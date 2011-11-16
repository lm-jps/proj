# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
RDVINV_$(d)	:= 
ifeq ($(JSOC_MACHINE), linux_x86_64)
RDVINV_$(d)	:= $(addprefix $(d)/, rdvinv)
endif

ifeq ($(COMPILER), icc)
  MODEXE_$(d)	:= $(addprefix $(d)/, datavg maicalc mtrack pspec3 rdcover rdfitc) 
endif
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

MODEXE_USEF_$(d):= $(RDVINV_$(d))
MODEXE_USEF	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

CEXE_$(d)	:= $(addprefix $(d)/, gentargs) 
CEXE		:= $(CEXE) $(CEXE_$(d))

EXE_$(d)        := $(MODEXE_$(d)) $(MODEXE_USEF_$(d) $(CEXE_$(d))
OBJ_OLAXY_$(d)	:= $(addprefix $(d)/, ola_xy_v10.o ola_subs.o)
OBJ_$(d)	:= $(EXE_$(d):%=%.o) $(OBJ_OLAXY_$(d))
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
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I/home/jsoc/include -DCDIR="\"$(SRCDIR)/$(d)\""

ifeq ($(JSOC_MACHINE), linux_ia32)
  FFTW_$(d) = /home/jsoc/lib/linux-ia32
endif
ifeq ($(JSOC_MACHINE), linux_x86_64)
  FFTW_$(d) = /home/jsoc/lib/linux-x86_64
endif

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) -L$(FFTW_$(d)) -lfftw3 -lfftw3f

ifeq ($(JSOC_MACHINE), linux_x86_64)
$(RDVINV_$(d)):				$(OBJ_OLAXY_$(d))
$(RDVINV_$(d)):				LL_TGT := $(LL_TGT) -lmkl_em64t
endif

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
