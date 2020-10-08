# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
MODEXE_$(d)	:= $(addprefix $(d)/, nlfff_v1 localpf test_nlfff test_localpf test_preproc nlfff_v2 localpf_mkl decomp_mkl)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

#MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
#MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

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
$(OBJ_$(d)):		CF_TGT := -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/src/ -fp-model precise
$(EXE_$(d)):		LF_TGT := $(LF_TGT) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_rt #-lmkl_em64t 

# I removed the compiler flags "-fp-model precise" and "-fp-model source" 
# from the LINK command.  If the module really needs such precise handling 
# of floating-point operations, then the flags will have to be added to the 
# COMPILE command.  Keh-Cheng 2012.08.31    

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
