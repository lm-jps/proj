# Standard things                                                                                             
sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# OUTDIR is the root directory that contains all make-generated binary files
# PDIR is the sub-directory of OUTDIR that is the parent 
OUTDIR		:= $(WORKINGDIR)/_$(MACH)
PSDIR		:= $(d)

# define here all targets built in the matlab directory
MEX2CSDIR	:= mex/src/mex2c
CSDIR		:= $(PSDIR)/$(MEX2CSDIR)
LIBMEX2C	:= $(OUTDIR)/$(PSDIR)/$(MEX2CSDIR)/libmex2c_lib.a
LIBMEX2MATL	:= $(OUTDIR)/$(PSDIR)/$(MEX2CSDIR)/libmex2matl.a

S_$(d)		:= $(notdir $(LIBMEX2C) $(LIBMEX2MATL))
S1_$(d)		:= $(addprefix $(PSDIR)/$(MEX2CSDIR)/, $(S_$(d)))

# Don't rely upon $(PWD) to tell you make's current working directory. It lies! Instead of trying 
# to come up with a path relative to where you think make is, just avoid relative paths altogether.
$(LIBMEX2C) $(LIBMEX2MATL):
	$(MAKE) -C $(WORKINGDIR)/$(PSDIR)/$(MEX2CSDIR) OUTDIR=$(OUTDIR) CSDIR=$(CSDIR) $@
#        $(MAKE) -C mfile-mex OUTDIR=$(CURDIR)/$(d) $@

.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(OUTDIR)/$(PSDIR)/$(MEX2CSDIR)/%

.PHONY:	$(S1_$(d))
$(S1_$(d)):	%:	$(OUTDIR)/%

TGT_LIB         := $(TGT_LIB) $(LIBMEX2MATL) $(LIBMEX2C)

# Standard things                                                                
-include        $(DEP_$(d))

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
