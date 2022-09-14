# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables


MODEXE_$(d)	:= $(addprefix $(d)/, store_ds_fd_M_01h_lev1_8 store_ds_fd_M_96m_01d_lev1_8 store_ds_fd_V_01h_lev1_8 store_ds_fd_V_30s_01h_lev1_8 store_ds_vw_V_06h_lev1_8 ingest_dsds_a store_dsds_migrate ingest_dsds_to_drms set_gaps_missing_v2)

MODEXESUMS_$(d)	:= $(addprefix $(d)/, sum_pe_svc sum_pe) $(STORDSDS_$(d))
MODEXE		:= $(MODEXE) $(MODEXE_$(d)) $(MODEXESUMS_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXESUMS_SOCK_$(d)	:= $(MODEXESUMS_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d)) $(MODEXESUMS_SOCK_$(d))

EXE_$(d)	:= $(MODEXE_$(d)) $(MODEXESUMS_$(d))
OBJ_$(d)	:= $(EXE_$(d):%=%.o)
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
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
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../libs/

# do not use $(LIBSUMSAPI) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
# do not use $(LIBSUM) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(MODEXESUMS_$(d)) $(MODEXESUMS_SOCK_$(d)):	base/sums/libs/api/libsumsapi.a base/sums/libs/pg/libsumspg.a
# $(LIBSUM) depends on libmisc.a, but the latter appears in the dependency list AFTER
# $(LIBSUM).  Force make to find libmisc.a
$(MODEXESUMS_$(d)) $(MODEXESUMS_SOCK_$(d)):	LL_TGT := $(LL_TGT) -Lbase/libs/misc -lmisc

# do not use $(LIBDSDSMIGR) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	 proj/dsdsmigr/libs/libdsdsmigr.a

# $(ingest_dsds_a_$(d)):     LL_TGT := $(LL_TGT) -L/usr/lib -lrt

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
