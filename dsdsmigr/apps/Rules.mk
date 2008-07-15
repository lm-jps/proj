# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
#STORDSDS_$(d)  := $(addprefix $(d)/, store_dsds_2_drms)
STORDSDS_EXTRA_LIBS = $(LIBSOIJSOC)

MODEXE_$(d)	:= $(addprefix $(d)/, store_ds_fd_M_01h_lev1_8 store_ds_fd_M_96m_01d_lev1_8 store_ds_fd_V_01h_lev1_8 store_ds_fd_V_30s_01h_lev1_8 store_ds_vw_V_06h_lev1_8 ingest_dsds_a store_dsds_migrate)
MODEXESUMS_$(d)	:= $(addprefix $(d)/, sum_pe_svc sum_pe) $(STORDSDS_$(d))
MODEXE		:= $(MODEXE) $(MODEXE_$(d)) $(MODEXESUMS_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock) 
MODEXESUMS_SOCK_$(d)	:= $(MODEXESUMS_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d)) $(MODEXESUMS_SOCK_$(d))

EXE_$(d)	:= $(MODEXE_$(d)) $(MODEXESUMS_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o) 
DEP_$(d)	:= $(OBJ_$(d):%=%.o.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(MODEXE_SOCK_$(d)) \
		   $(MODEXESUMS_SOCK_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXESUMS_SOCK_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXESUMS_SOCK_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""

$(MODEXESUMS_$(d)) $(MODEXESUMS_SOCK_$(d)):	$(LIBSUMSAPI) $(LIBSUM)

#$(STORDSDS_$(d)) $(STORDSDS_$(d):%=%_sock): LL_TGT := $(STORDSDS_EXTRA_LIBS)
#$(STORDSDS_$(d)) $(STORDSDS_$(d):%=%_sock): $(STORDSDS_EXTRA_LIBS)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
