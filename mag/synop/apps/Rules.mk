# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Libs needed by others.


# Common utilities
# EXTRADEPS_$(d)          := $(addprefix $(d)/, synop-saveparm.o synop-timing.o synop-set_history.o synop-calversfunctions.o remap-img2helioVector.o remap-obs2helio.o remap-apodize.o remap-setplm2.o)

EXTRADEPS_$(d)          := $(addprefix $(d)/, synop-saveparm.o synop-timing.o synop-set_history.o synop-calversfunctions.o)

# NOTE: Add the base of the module's filename below (next to mymod)
MODEXE_$(d)	:= $(addprefix $(d)/, maproj3comperrorlonat02deg4corr hmib3compsynoptic4corr mrmlossynoptic4corr mrmlosdailysynframe4mdi mrmlossynoptic4mdi hmib3compsynoptic4cr2136 hmib3dailysynframe vectmag2helio3comp_radial vectmag2helio3comp_poten hmib3compsynoptic resizeb3comp vectmag2helio3comp_random maprojlonat02deg maproj3comperrorlonat02deg resizeb3compwitherror maproj3comperror dailysynframe dailysynframe_nrt hmisynoptic mdidailysynframe mdisynop brblossynoptic brblosdailysynframe brblosdailysynframe_nrt mrmlossynoptic mrmlosdailysynframe mrmlosdailysynframe_nrt)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# Modules with external libraries
MODEXE_USEF_$(d)	:= $(addprefix $(d)/, )
MODEXE_USEF		:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

MODEXE_USEF_SOCK_$(d)	:= $(MODEXE_USEF_$(d):%=%_sock)
MODEXE_USEF_SOCK	:= $(MODEXE_USEF_SOCK) $(MODEXE_USEF_SOCK_$(d))

CUSTOM_CARTOGRAPHY_$(d) := $(addprefix $(d)/, synop-imginfo.o synop-cartography.o)

EXE_$(d)	:= $(MODEXE_$(d)) $(MODEXE_USEF_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o) $(CUSTOM_CARTOGRAPHY_$(d)) $(EXTRADEPS_$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(MODEXE_USEF_SOCK_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := -I$(SRCDIR)/$(d) -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I$(SRCDIR)/$(d)/../../../globalhs/libs/projection -I$(SRCDIR)/$(d)/../../../libs/projection -I$(SRCDIR)/$(d)/../../libs/util -I$(SRCDIR)/$(d)/src/ $(FFTWH) -I$(SRCDIR)/$(d)/../../../globalhs/libs/pkbgn -I$(SRCDIR)/$(d)/../../../libs/interpolate
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

$(EXTRADEPS_$(d)):	CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)

ALL_$(d)	:= $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(MODEXE_USEF_SOCK_$(d))
$(ALL_$(d)) : $(EXTRADEPS_$(d)) $(CUSTOM_CARTOGRAPHY_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBCARTOGRAPHY) $(LIBMAGUTILS) $(LIBSTATS) $(LIBPROJECTION) $(LIBINTERP)
$(ALL_$(d)) : LF_TGT := $(LF_TGT) $(MKL)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) $(FFTW3LIBS) $(CFITSIOLIBS) -lmkl_em64t

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
ifneq ($(DEP_$(d)),)
  -include	$(DEP_$(d))
endif

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
