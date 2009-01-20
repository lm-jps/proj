# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
test0_$(d)		:= $(addprefix $(d)/, test0)
mainx_$(d)		:= $(addprefix $(d)/, mainx)
jpe_$(d)		:= $(addprefix $(d)/, jpe)
jpe_obj_$(d)		:= $(addprefix $(d)/, hdata.o start_pvmd.o call_dsds.o pack_key.o prod_host.o spawnds.o parse_arg.o names.o du_dir.o call_drms_in.o)

jpe2_$(d)		:= $(addprefix $(d)/, jpe2)
jpe2_obj_$(d)		:= $(addprefix $(d)/, hdata.o start_pvmd.o call_dsds.o pack_key.o prod_host.o spawnds.o parse_arg.o names.o du_dir.o call_drms_in.o)

SUMEXE_$(d)	:= $(test0_$(d)) $(jpe_$(d)) $(jpe2_$(d)) $(mainx_$(d))
CEXE_$(d)       := $(addprefix $(d))
CEXE		:= $(CEXE) $(CEXE_$(d))
MODEXESUMS	:= $(MODEXESUMS) $(SUMEXE_$(d)) $(PEEXE_$(d))

MODEXE_$(d)	:= $(addprefix $(d))
MODEXEDR_$(d)	:= $(addprefix $(d))
MODEXE_USEF_$(d)	:= $(addprefix $(d))
MODEXE_USEF 	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))
MODEXEDR	:= $(MODEXEDR) $(MODEXEDR_$(d))
MODEXE		:= $(MODEXE) $(SUMEXE_$(d)) $(MODEXE_$(d)) $(MODEXEDR_$(d)) $(PEEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock) $(MODEXEDR_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))
MODEXEDR_SOCK	:= $(MODEXEDR_SOCK) $(MODEXEDR_$(d):%=%_sock)

MODEXEDROBJ	:= $(MODEXEDROBJ) $(MODEXEDR_$(d):%=%.o)

ALLEXE_$(d)	:= $(MODEXE_$(d)) $(MODEXEDR_$(d)) $(MODEXE_USEF_$(d)) $(SUMEXE_$(d)) $(CEXE_$(d))
OBJ_$(d)	:= $(ALLEXE_$(d):%=%.o) $(jpe_obj_$(d)) $(jpe2_obj_$(d)) 
DEP_$(d)	:= $(OBJ_$(d):%=%.o.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(ALLEXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(ALLEXE_$(d)) $(MODEXE_SOCK_$(d)) 

S_$(d)		:= $(notdir $(ALLEXE_$(d)) $(MODEXE_SOCK_$(d)))

$(jpe_$(d)):     $(jpe_obj_$(d))
$(jpe2_$(d)):     $(jpe2_obj_$(d))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(SUMEXE_$(d)):		LL_TGT := -L/home/production/cvs/jsoc/lib/saved/$(JSOC_MACHINE) -lhmicomp_egse -lecpg -lpq -lpng -lsoi  -L/home/soi/CM/pvm3/lib/LINUX8664 -lpvm3 -L/home/production/cvs/JSOC/lib/$(JSOC_MACHINE) -lsoi
$(PEEXE_$(d)):		LL_TGT := -lecpg -lpq 
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I/home/soi/CM/pvm3/include -DPVM33
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)):	$(LIBASTRO)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
