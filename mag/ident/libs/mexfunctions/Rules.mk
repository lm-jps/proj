# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
MYCMPFLG_$(d) := -DMEX2C_TAIL_HOOK -DUSING_MEX2LIB -DStaticP=static -I$(SRCDIR)/$(d)/../util -I$(SRCDIR)/$(d)/../mex2c -I$(SRCDIR)/$(d) -I$(SRCDIR)/$(d)/Doc

# just one library now
LIBsegment	:= $(d)/libmxt_segment.a
# its contents
OBJsegment_$(d)	:= $(addprefix $(d)/, hmi_segment.o clean_edge_label.o makemrfdiscwts.o mixNprob2d.o mrf_segment_wts.o roi_stats_mag.o smoothsphere.o)

# globals
OBJ_$(d)	:= $(OBJsegment_$(d))

DEP_$(d)	:= $(OBJ_$(d):%=%.d)

CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(LIBsegment) $(DEP_$(d))

S_$(d)		:= $(notdir $(LIBsegment))

# Local rules
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) $(MYCMPFLG_$(d))

$(LIBsegment): $(OBJsegment_$(d))
	$(ARCHIVE)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
