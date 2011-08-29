# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
MODEXE_$(d)	:= $(addprefix $(d)/, test_nlfff test_localpf test_preproc)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

#MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
#MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# flags for compiling and linking
MYCMPFLG	:= -openmp
MYLNKFLG	:= -lmkl_em64t -openmp

EXE_$(d)	:= $(MODEXE_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) #$(MODEXE_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d))) #$(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJ_$(d)):		CF_TGT := -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) -I$(SRCDIR)/lib_third_party/include $(MYCMPFLG)
$(EXE_$(d)):		LF_TGT := $(LF_TGT) -lmkl_em64t -openmp -fp-model precise -fp-model source

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
