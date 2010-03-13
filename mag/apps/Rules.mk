# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/, ) 
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

MODEXE_USEF_$(d):= $(RDVINV_$(d))
MODEXE_USEF	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

EXE_$(d)        := $(MODEXE_$(d)) $(MODEXE_USEF_$(d))
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
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I/home/jsoc/include -DCDIR="\"$(SRCDIR)/$(d)\""

ifeq ($(JSOC_MACHINE), linux_ia32)
  FFTW_$(d) = /home/jsoc/lib/linux-ia32
endif
ifeq ($(JSOC_MACHINE), linux_x86_64)
  FFTW_$(d) = /home/jsoc/lib/linux-x86_64
endif

$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	LL_TGT := $(LL_TGT) -L$(FFTW_$(d)) -lfftw3 -lfftw3f

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
