# Standard things
sp := $(sp).x
dirstack_$(sp) := $(d)
d := $(dir)

## never touch above ---------------------------------------------------------------

## edit below -----------------------------------------------------

# module name(s) to be generated and file name(s) to be compiled

## C-wrapper name (name must end with .c)
MODEXE_USEF_$(d) := $(addprefix $(d)/, cgem_pdfi)
MODEXE_USEF := $(MODEXE_USEF) $(MODEXE_USEF_$(d))

## wrapped Fortran codes
WRAPPEDF_OBJ_$(d)    := $(addprefix $(d)/, add_padding_ss.o ahpot_ss.o angle_be_ss.o berciktest_ss.o bhll2tp_ss.o bhpot_phot_ss.o bhpot_ss.o bhyeell2tp_ss.o bhyeetp2ll_ss.o br_voxels3d_ss.o brll2tp_ss.o brpot_ss.o bryeell2tp_ss.o bryeetp2ll_ss.o cell_ss.o cetp_ss.o coell_ss.o coetp_ss.o coll_ss.o cotp_ss.o curl_psi_rhat_ce_ss.o curl_psi_rhat_co_ss.o curlahpot_ss.o curle_ss.o curlehr_ss.o curlh_ce_ss.o curlh_co_ss.o delh2_sc.o divh_ce_ss.o divh_co_ss.o divh_sc.o e_doppler_ss.o e_flct_ss.o e_ideal_ss.o e_ptd_ss.o e_voxels3d_ss.o ehyeell2tp_ss.o ehyeetp2ll_ss.o enudge3d_ss.o enudge_ss.o eryeell2tp_ss.o eryeetp2ll_ss.o find_mask_ss.o fix_mask_ss.o gradh_ce_ss.o gradh_co_ss.o gradh_sc.o hm_ss.o hmtot_ss.o interp_data_ss.o interp_eh_ss.o interp_var_ss.o kfft.o pad_abcd_ss.o pdfi_wrapper4jsoc_ss.o pell_ss.o petp_ss.o ptdsolve_ss.o relax_psi_3d_ss.o scrbpot_ss.o sinthta_sc.o sinthta_ss.o sr_ss.o srtot_ss.o tell_ss.o tetp_ss.o)

# flags for compiling and linking
MYCMPFLG_$(d) := -O2 -nofor-main -fp-model strict
MYLNKFLG_$(d) := -lmkl_em64t

## edit above -----------------------------------------------------

## may not edit below -----------------------------------------------------

## compile and link etc.

# list .o .o.d and executables to be generated
#

OBJ_$(d) := $(MODEXE_USEF_$(d):%=%.o) $(WRAPPEDF_OBJ_$(d))
DEP_$(d) := $(OBJ_$(d):%=%.d)
CLEAN := $(CLEAN) \
         $(OBJ_$(d)) \
         $(MODEXE_USEF_$(d)) \
         $(DEP_$(d))


# name(?) of target executable binary
TGT_BIN := $(TGT_BIN) $(MODEXE_USEF_$(d))

S_$(d) := $(notdir $(MODEXE_USEF_$(d)))

# making object/executable with Rules file and extra options
$(OBJ_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):	CF_TGT := $(CF_TGT) -fp-model strict -I$(SRCDIR)/$(d)/../../../libs/astro -I$(SRCDIR)/$(d)/../../../libs/stats -I/home/jsoc/lib/$(JSOC_MACHINE)
# -I$(SRCDIR)/$(d)/../../../libs/interpolate $(FFTWH) $(FFTW3LIBS) $(FMATHLIBSH) 

$(WRAPPEDF_OBJ_$(d)): FF_TGT := $(FF_TGT) $(MYCMPFLG_$(d))
$(MODEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(MYLNKFLG_$(d)) -lfishpack_r8
$(MODEXE_USEF_$(d)): $(WRAPPEDF_OBJ_$(d))

ALL_$(d)	:= $(MODEXE_USEF_$(d))
$(ALL_$(d)) : $(LIBASTRO) $(LIBSTATS) # $(LIBINTERP)
$(ALL_$(d)) : LL_TGT := $(LL_TGT) # $(GSLLIBS)

## never touch below ---------------------------------------------------------------------------

# Shortcuts
.PHONY: $(S_$(d))
$(S_$(d)): %: $(d)/%

# Standard things
-include $(DEP_$(d))

d := $(dirstack_$(sp))
sp := $(basename $(sp))

# end of this file ---------------------------------------------------------------------------------
