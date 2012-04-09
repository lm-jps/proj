# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# Builds on icc/x86_64 only
ifeq ($(JSOC_MACHINE), linux_x86_64)
ifeq ($(COMPILER), icc)
VFISV_$(d)	:= $(addprefix $(d)/, vfisv vfisv_harp)
endif
endif

VFISV_COBJ_$(d)	:= $(VFISV_$(d):%=%.o)

# VFISV_$(d) depends on Fortran code
VFISV_FOBJ_$(d)	:= $(addprefix $(d)/, $(notdir $(patsubst %.f,%.o,$(wildcard $(SRCDIR)/$(d)/*.f))))
VFISV_FOBJ_$(d)	:= $(VFISV_FOBJ_$(d)) $(addprefix $(d)/, $(notdir $(patsubst %.f90,%.o,$(wildcard $(SRCDIR)/$(d)/*.f90))))

# Actually, we don't need dmrs.f90 for now
# EXCLUDE_$(d)	:= $(d)/dmrs.o
# VFISV_FOBJ_$(d) := $(filter-out $(EXCLUDE_$(d)),$(VFISV_FOBJ_$(d)))

MODEXE_USEF_$(d):= $(VFISV_$(d))
MODEXE_USEF	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

# Not available as a sock module - uses drms_server_end_transaction()
# MODEXE_USEF_SOCK_$(d)	:= $(addsuffix _sock, $(MODEXE_USEF_$(d)))
# MODEXE_USEF_SOCK	:= $(MODEXE_USEF_SOCK) $(MODEXE_USEF_SOCK_$(d))

#EXE_$(d)        := $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
EXE_$(d)        := $(MODEXE_USEF_$(d))
OBJ_$(d)	:= $(VFISV_COBJ_$(d)) $(VFISV_FOBJ_$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I/home/jsoc/include -I/home/jsoc/mpich2/include -DCDIR="\"$(SRCDIR)/$(d)\""
$(EXE_$(d)):		LL_TGT := $(LL_TGT) -lmkl_em64t

# Use a non-standard Fortran compiler for LINKING the modules of this make file.
$(EXE_$(d)):		FCOMPILER := /home/jsoc/mpich2/bin/mpif90
# The line above should NOT affect the compilation of the object files, but it does!
$(VFISV_FOBJ_$(d)):	FCOMPILER := ifort
$(VFISV_COBJ_$(d)):	ICC_CMPLR := /home/jsoc/mpich2/bin/mpicc

# Here's the skinny on Fortran 90 modules. If a Fortran 90 file 
# contains a MODULE statement, then a .mod file will be produced by the compiler. To use a module file
# a Fortran 90 file can use the USE statement. However, when you COMPILE
# a Fortran 90 file with a USE statement in it, there must exist a .mod file for that USE statement.
# If a Fortran 90 file has USE FILT_PARAM, then AT COMPILE TIME filt_param.mod must exist. If you want to 
# link a bunch of Fortran 90 files together, then this implies that you must compile the Fortran 90
# files in a specific order - you must compile dependent files before you compile the files that require
# the dependent files. So we have to specify the order of compilation of the f90 files (by specifying
# the dependency relationships between them.
$(d)/change_var.o: $(d)/cons_param.o  $(d)/filt_param.o
$(d)/filt_init.o: $(d)/filt_param.o $(d)/cons_param.o $(d)/line_param.o
$(d)/filt_param.o: $(d)/line_param.o $(d)/cons_param.o
$(d)/forward.o: $(d)/line_param.o $(d)/cons_param.o $(d)/filt_param.o $(d)/inv_param.o $(d)/voigt.o $(d)/voigt_taylor.o $(d)/change_var.o
$(d)/free_init.o: $(d)/inv_param.o $(d)/cons_param.o $(d)/filt_param.o
$(d)/free_memory.o: $(d)/filt_param.o $(d)/line_param.o
$(d)/inv_param.o: $(d)/cons_param.o
$(d)/inv_utils.o: $(d)/cons_param.o $(d)/filt_param.o $(d)/inv_param.o $(d)/line_param.o $(d)/svdcmp.o $(d)/svbksb.o $(d)/ran_mod.o
$(d)/inv1_init.o: $(d)/inv_param.o
$(d)/wfa_guess.o: $(d)/cons_param.o $(d)/filt_param.o $(d)/inv_param.o
$(d)/ran_mod.o: $(d)/cons_param.o
$(d)/invert.o: $(d)/invert.f90 $(d)/forward.o $(d)/line_param.o $(d)/cons_param.o $(d)/filt_param.o $(d)/inv_utils.o $(d)/inv_param.o $(d)/svdcmp.o $(d)/svbksb.o $(d)/wfa_guess.o $(d)/change_var.o
$(d)/line_init.o: $(d)/line_param.o $(d)/cons_param.o
$(d)/line_param.o: $(d)/cons_param.o
$(d)/svbksb.o: $(d)/cons_param.o
$(d)/svdcmp.o: $(d)/cons_param.o
$(d)/voigt_data.o: $(d)/cons_param.o
$(d)/voigt_init.o: $(d)/voigt_data.o $(d)/cons_param.o
$(d)/voigt_taylor.o: $(d)/voigt_data.o $(d)/voigt_init.o $(d)/voigt.o $(d)/cons_param.o $(d)/line_param.o
$(d)/wave_init.o: $(d)/line_param.o $(d)/cons_param.o

$(VFISV_$(d)):		$(VFISV_FOBJ_$(d))

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
