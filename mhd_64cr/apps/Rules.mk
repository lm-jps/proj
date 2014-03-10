## Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)
## do not touch above.

## C to be compiled
MODEXE_$(d)	:= $(addprefix $(d)/, mhdtxt2jsoc_64cr)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))
## stand-alone fortran to be compiled
FEXES_$(d) 	:= $(addprefix $(d)/, hmisynofits2txt_2 mhd2equidist4_64_fits)
FEXED_$(d) 	:= $(addprefix $(d)/, mhd_64raw mhd_64pc5 hmitxt2wso mhd2equidist_64 mhd2equidist_txt4jsoc_72x64x128to72x144x80 mhdati2csv_64 mkcgsw9w_64cr slcalt15_64cr xy2uv_64)
FEXES		:= $(FEXES_$(d))
FEXED		:= $(FEXED_$(d))
FEXE_$(d) 	:= $(FEXES_$(d)) $(FEXED_$(d))
FEXE		:= $(FEXE_$(d))

## list of exec., objective, and dependency file, and clean-up list
# EXE_$(d)        := $(MODEXE_$(d))
OBJ_$(d)	:= $(MODEXE_$(d):%=%.o)
OBJFS_$(d)	:= $(FEXES_$(d):%=%.o)
OBJFD_$(d)	:= $(FEXED_$(d):%=%.o)
OBJF_$(d)       := $(OBJFS_$(d))  $(OBJFD_$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.d) $(OBJF_$(d):%=%.d) 
CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(OBJF_$(d)) $(MODEXE_$(d)) $(FEXE_$(d)) $(DEP_$(d))

## list of exec. binary to be created?
TGT_BIN	        := $(TGT_BIN) $(MODEXE_$(d)) $(FEXE_$(d))
## shortcut target list
S_$(d)		:= $(notdir $(MODEXE_$(d))) $(notdir $(FEXE_$(d)))

## Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) -I$(SRCDIR)/lib_third_party/include
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/stats
$(EXE_$(d)):		LF_TGT := $(LF_TGT)
$(OBJFS_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJFS_$(d)):		FF_TGT := $(FF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJFS_$(d)):		FF_TGT := $(FF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) -I$(SRCDIR)/lib_third_party/include
$(OBJFS_$(d)):		FF_TGT := $(FF_TGT) -I$(SRCDIR)/$(d)/../../libs/stats
$(FEXES_$(d)):		LL_TGT := $(LL_TGT) -lpq $(CFITSIOLIBS)
$(OBJFD_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJFD_$(d)):		FF_TGT := $(FF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -O3 -ip -r8
$(FEXED_$(d)):		LL_TGT := $(LL_TGT)

## Shortcuts, do not touch below
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%
# Standard things
-include	$(DEP_$(d))
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
## end of this file
