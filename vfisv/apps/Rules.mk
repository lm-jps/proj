# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# Builds on icc/x86_64 only
#ifeq ($(JSOC_MACHINE), linux_x86_64)
#ifeq ($(COMPILER), icc)
VFISV_NOT2COMP$(d)	:= $(addprefix $(d)/, vfisv vfisv_harp vfisv_fd10_old vfisv_fd10_harp_old)
VFISV_2COMP$(d)		:= $(addprefix $(d)/, vfisv_2comp)
#endif
#endif

VFISV_COBJ_$(d)	:= $(VFISV_NOT2COMP$(d):%=%.o) $(VFISV_2COMP$(d):%=%.o)

# VFISV_$(d) depends on Fortran code. 
# Here's the skinny on Fortran 90 modules. If a Fortran 90 file 
# contains a MODULE statement, then a .mod file will be produced by the compiler. To use a module file
# a Fortran 90 file uses the USE statement. However, when you COMPILE
# a Fortran 90 file with a USE statement in it, there must exist a .mod file for the module
# in the USE statement. If a Fortran 90 file has 'USE FILT_PARAM', then AT COMPILE TIME 
# filt_param.mod must exist. If you want to link a bunch of Fortran 90 files together, then this implies 
# that you must compile the Fortran 90 files in a specific order - you must compile dependent files 
# before you compile the files that require the dependent files. 
#
# This dependency chain was previous enforced with make rules between the various .o files
# generated from the .f90 files. However, this was causing the timestamps of dependent files
# to be newer than the timestamps of the files that depend on them. So make vfisv would always
# build something every time make ran, even though this was unnecessary.  Instead, we now specify the order of 
# compilation of the f90 files by listing them in the order in which they should be built.

# Compile in this order (this is ONE order that works - there are others). See CVS version
# 1.11 for the list of actual dependencies between these Fortran object files.
#   cons_param.o
#   line_param.o
#   filt_param.o
#   change_var.o
#   filt_init.o
#   inv_param.o
#   inv_init.o
#   lim_init.o
#   wfa_guess.o
#   ran_mod.o
#   line_init.o
#   svbksb.o
#   voigt_data.o
#   voigt_init.o
#   wave_init.o
#   voigt.o
#   voigt_taylor.o
#   forward_2comp.o
#   free_init.o
#   free_memory.o
#   inv_utils.o
#   invert_2comp.o

# Just list all the Fortran files since we have to build them in a specific order (see note below).
# They will be built in the order they appear in this list.
VFISV_FOBJ_$(d) := $(addprefix $(d)/, cons_param.o \
    line_param.o \
    filt_param.o \
    change_var.o \
    filt_init.o \
    inv_param.o \
    inv_init.o \
    lim_init.o \
    wfa_guess.o \
    ran_mod.o \
    line_init.o \
    svbksb.o \
    voigt_data.o \
    voigt_init.o \
    wave_init.o \
    voigt.o \
    voigt_taylor.o \
    free_init.o \
    free_memory.o \
    inv_utils.o)

VFISV_FOBJ_NOT2COMP$(d) := $(addprefix $(d)/, forward.o invert.o)
VFISV_FOBJ_2COMP$(d) := $(addprefix $(d)/, forward_2comp.o invert_2comp.o)

MODEXE_USEF_$(d):= $(VFISV_NOT2COMP$(d)) $(VFISV_2COMP$(d))
MODEXE_USEF	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

EXE_$(d)        := $(MODEXE_USEF_$(d))
OBJ_$(d)	:= $(VFISV_COBJ_$(d)) $(VFISV_FOBJ_$(d)) $(VFISV_FOBJ_NOT2COMP$(d)) $(VFISV_FOBJ_2COMP$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(MPI_INCS)
ifeq ($(JSOC_MACHINE), linux_avx2)
$(EXE_$(d)):            LL_TGT := $(LL_TGT) -lmkl_rt
else
$(EXE_$(d)):		LL_TGT := $(LL_TGT) -lmkl_em64t
endif

# Use a non-standard Fortran compiler for LINKING the modules of this make file.
$(EXE_$(d)):		FCOMPILER := $(MPIFCOMPILER)
# The line above should NOT affect the compilation of the object files, but it does!
$(VFISV_FOBJ_$(d)):		FCOMPILER := ifort
$(VFISV_FOBJ_NOT2COMP$(d)):     FCOMPILER := ifort
$(VFISV_FOBJ_2COMP$(d)):     	FCOMPILER := ifort
$(VFISV_COBJ_$(d)):		ICC_CMPLR := $(MPICOMPILER)

$(VFISV_NOT2COMP$(d)):		$(VFISV_FOBJ_$(d)) $(VFISV_FOBJ_NOT2COMP$(d))
$(VFISV_2COMP$(d)):		$(VFISV_FOBJ_$(d)) $(VFISV_FOBJ_2COMP$(d))


# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
