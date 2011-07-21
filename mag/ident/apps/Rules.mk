# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
MODEXE_$(d)	:= $(addprefix $(d)/, hmi_segment_module) # put module names here
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

EXE_$(d)	:= $(MODEXE_$(d))

OBJsegment_$(d) := $(addprefix $(d)/, cJSON.o segment_modelset.o)

# these files are needed by hmi_segment_module
OBJ_$(d)	:= $(EXE_$(d):%=%.o) $(OBJsegment_$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)))

# flags for compiling and linking
MYCMPFLG_$(d)  := -I$(SRCDIR)/$(d)/../libs/util -I$(SRCDIR)/$(d)/../libs/mex2c  -I$(SRCDIR)/$(d)/../libs/mexfunctions

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""	# append to $(CF_TGT)
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/src $(FMATHLIBSH) $(MYCMPFLG_$(d))

$(EXE_$(d)):		$(OBJsegment_$(d)) $(LIBsegment) $(LIBmex2c_jsoc) $(LIBmexrng) $(LIBmextools)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
