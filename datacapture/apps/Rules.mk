# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables

# Force icc compilation since this doesn't work for gcc because libhmicomp_egse was built with icc.
LOCALCC		= $(ICC_COMP)
LOCALLN		= $(ICC_LINK)

EXE_$(d)	:= $(addprefix $(d)/, ingest_tlm soc_pipe_scp)
CEXESUMS	:= $(CEXESUMS) $(EXE_$(d))

OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(EXE_$(d):%=%.o.d) 

CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(DEP_$(d)) 

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJ_$(d)):	%.o:	%.c
		$(LOCALCC)

$(EXE_$(d)):	LL_TGT :=  -L/home/production/cvs/jsoc/lib/saved/$(JSOC_MACHINE) -lhmicomp_egse -lecpg -lpq
$(EXE_$(d)):	%:	%.o $(EXELIBS)
		$(LOCALLN)
		$(SLBIN)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
