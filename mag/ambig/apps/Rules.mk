# Standard things
sp := $(sp).x
dirstack_$(sp) := $(d)
d := $(dir)

## never touch above ---------------------------------------------------------------

## edit below -----------------------------------------------------

# module name(s) to be generated and file name(s) to be compiled

## C-wrapper name (name must end with .c)
# This doesn't build on ia32
ifeq ($(JSOC_MACHINE), linux_x86_64)
  MODEXE_USEF_$(d) := $(addprefix $(d)/, test_ambig)
endif

## wrapped Fortran codes
WRAPPEDF_$(d)    := $(addprefix $(d)/, set_geometry.o sizes.o mgram_data.o pad.o point.o pot_field.o verbose.o WeightingFactor.o anneal.o ranseed.o bobs.o constant.o disk_center.o ephemeris.o pix_size.o spherical_position.o bounds.o maskvec.o deriv_coefficients.o spherical_deriv_coefficients.o energy_arrays.o recon.o ran_pix.o ran3.o trnsfrm.o transform.o buffer.o cuspl.o cuspl2d.o cuspleval.o cuspleval2d.o smooth.o mkl_dfti.o dfti_example_status_print.o potential.o pacute.o boxit.o setup_OCBP_PF_dzh_4p.o setup_spherical_PF_4p.o CalcE_OCBP_PF_dzh_4p.o CalcDE_reconfig_OCBP_PF_dzh_4p.o CalcE_spherical_PF_4p.o CalcDE_reconfig_spherical_PF_4p.o get_ran_pix.o CalcDE.o reconfig.o minimise_energy.o global.o sortrx.o nacute5.o nacute6.o nacute6p.o mollweide.o revmollweide.o tile.o colatlon.o grown.o ambig.o)

# flags for compiling and linking
MYCMPFLG := -O2 -free -nofor-main
MYLNKFLG := -lmkl_em64t -openmp

## edit above -----------------------------------------------------

## may not edit below -----------------------------------------------------

## compile and link etc.

# list .o .o.d and executables to be generated
#
MODEXE_USEF := $(MODEXE_USEF) $(MODEXE_USEF_$(d))
#
OBJ_$(d) := $(MODEXE_USEF_$(d):%=%.o) $(WRAPPEDF_$(d))
#
DEP_$(d) := $(OBJ_$(d):%=%.d)

CLEAN := $(CLEAN) $(OBJ_$(d)) $(MODEXE_USEF_$(d)) $(DEP_$(d))

# name(?) of target executable binary
TGT_BIN := $(TGT_BIN) $(MODEXE_USEF_$(d))
# what is this ?
S_$(d) := $(notdir $(MODEXE_USEF_$(d)))

# making object/executable with Rules file and extra options
$(OBJ_$(d)): $(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)): CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -DCDIR="\"$(SRCDIR)/$(d)\""

$(WRAPPEDF_$(d)): FF_TGT := $(FF_TGT) $(MYCMPFLG)
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(MYLNKFLG)
$(MODEXE_USEF_$(d)): $(WRAPPEDF_$(d))

## never touch below ---------------------------------------------------------------------------

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)): %: $(d)/%

# Standard things
-include $(DEP_$(d))

d := $(dirstack_$(sp))
sp := $(basename $(sp))

# end of this file ---------------------------------------------------------------------------------
