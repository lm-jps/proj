## Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)
## do not touch above.

## C to be compiled
MODEXE_$(d)	:= $(addprefix $(d)/, mhdtxt2jsoc_daily_32)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

## stand-alone fortran to be compiled
FEXES_$(d) 	:= $(addprefix $(d)/, mhd2equidist4_32_fits hmibrsynofits2bltxt_3)
FEXED_$(d) 	:= $(addprefix $(d)/, gtpotco7 mhd2equidist_32 mhd2equidist_txt4jsoc_72x32x64to36x72x80 mkcgsw9w_daily_32 whatis_latlon_now xy2uv_32 mhdati2csv_32 slcalt15_daily_32 whatistoday4_intel whatis_crlon_utnoon what_is_todaynoon what_is_todaynoon_tai mhd_32raw mhd_32pc5)
FEXE_$(d) 	:= $(FEXES_$(d)) $(FEXED_$(d))
FEXE		:= $(FEXE) $(FEXE_$(d))

## list of exec., objective, and dependency file, and clean-up list
# EXE_$(d)        := $(MODEXE_$(d))
OBJ_$(d)	:= $(MODEXE_$(d):%=%.o) $(FEXES_$(d):%=%.o)
OBJFD_$(d)	:= $(FEXED_$(d):%=%.o)
DEP_$(d)	:= $(OBJ_$(d):%=%.d) $(OBJFD_$(d):%=%.d) 
CLEAN		:= $(CLEAN) $(OBJ_$(d)) $(OBJFD_$(d)) $(MODEXE_$(d)) $(FEXE_$(d)) $(DEP_$(d))

## list of exec. binary to be created?
TGT_BIN	        := $(TGT_BIN) $(MODEXE_$(d)) $(FEXE_$(d))
## shortcut target list
S_$(d)		:= $(notdir $(MODEXE_$(d))) $(notdir $(FEXE_$(d)))

## Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/src/ $(FMATHLIBSH) -I$(SRCDIR)/$(d)/../../libs/stats
$(FEXES_$(d)):		LL_TGT := $(LL_TGT) -lpq $(CFITSIOLIBS)
$(OBJFD_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJFD_$(d)):		FF_TGT := $(FF_TGT) -O3 -ip -r8

## Shortcuts, do not touch below
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%
# Standard things
-include	$(DEP_$(d))
d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
## end of this file
