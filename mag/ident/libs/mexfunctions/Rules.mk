# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# BTW - having a file in the library directory referring to a file in the apps directory is not good practice (libraries should never depend on applications)
MYCMPFLG_$(d) := -DMEX2C_TAIL_HOOK -DUSING_MEX2LIB -DStaticP=static -I$(SRCDIR)/$(d)/../util -I$(SRCDIR)/$(d)/../mex2c -I$(SRCDIR)/$(d)/../mexfunctions -I$(SRCDIR)/$(d)/../mexfunctions/Doc -I$(SRCDIR)/$(d)/../../apps/Doc

# two libraries
LIBsegment	:= $(d)/libmxt_segment.a
LIBpatch	:= $(d)/libmxt_patch.a
# their composition
OBJsegment_$(d)	:= $(addprefix $(d)/, hmi_segment.o clean_edge_label.o makemrfdiscwts.o mixNprob2d.o mrf_segment_wts.o)
OBJpatch_$(d)	:= $(addprefix $(d)/, hmi_patch.o concomponent.o region_bb.o smoothsphere.o roi_stats_mag.o)

# globals
OBJ_$(d)	:= $(OBJsegment_$(d)) $(OBJpatch_$(d))

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBsegment) $(LIBpatch) $(DEP_$(d))

S_$(d)		:= $(notdir $(LIBsegment) $(LIBpatch))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) $(MYCMPFLG_$(d))

$(LIBsegment): $(OBJsegment_$(d))
	$(ARCHIVE)

$(LIBpatch):   $(OBJpatch_$(d))
	$(ARCHIVE)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
