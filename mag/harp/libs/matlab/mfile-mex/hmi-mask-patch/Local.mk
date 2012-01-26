#
# Local modifications, if any
#
# These point mex and mex2c to the mexfunction libraries

# over-ride auto-selection of mex source files: exclude library mexfiles
MEX_SRC_C := hmi_patch.c

# put -lmex2matl in again because hmihelpMW needs it too
# (this seems to only be a problem with gcc 3.3?)
LIBSUSED += -lhmihelpMW -lmex2matl
#CMEXLDFLAGS += -L.
CMEXLDFLAGS += -L$(OUTDIR)/$(CSDIR)
#CMEXPPFLAGS += 
#MEX2C_LIBS += -lhmihelpCL
CPPFLAGS +=
LDFLAGS += -L.

# extra files to be removed upon "make clean"
CLEAN_FILES += libhmihelpMW.a

# the core library of mexFunctions (mathworks version)
LIBHMIHELP	:= $(OUTDIR)/$(CSDIR)/libhmihelpMW.a
LIBHMIHELP_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, concomponentMW.o region_bbMW.o smoothsphereMW.o roi_stats_magMW.o)

$(LIBHMIHELP):	$(LIBHMIHELP_OBJ)
	ar crus $@ $^

# The rule to make the .o files from the .c files is in MexCodeRules.mk

# make the mathworks library a requirement for the mexfiles
$(MEX_FILES): $(LIBHMIHELP)

