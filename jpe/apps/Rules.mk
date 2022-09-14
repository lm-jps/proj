# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
mainx_$(d)		:= $(addprefix $(d)/, mainx)
jsnqry_$(d)		:= $(addprefix $(d)/, jsnqry)
jsnxqry_$(d)		:= $(addprefix $(d)/, jsnxqry)
jsnwqry_$(d)		:= $(addprefix $(d)/, jsnwqry)
jpe_$(d)		:= $(addprefix $(d)/, jpe)
jpe_obj_$(d)		:= $(addprefix $(d)/, hdata.o start_pvmd.o call_dsds.o pack_key.o prod_host.o spawnds.o parse_arg.o names.o du_dir.o call_drms_in.o)
jpeq_$(d)		:= $(addprefix $(d)/, jpeq)
jpeq_obj_$(d)		:= $(addprefix $(d)/, hdata.o start_pvmd.o call_dsds.o pack_key.o prod_host.o spawnds.o parse_arg.o names.o du_dir.o call_drms_inq.o)

shape_$(d)               := $(addprefix $(d)/, shape)
shape_obj_$(d)           := $(addprefix $(d)/, hdata.o start_pvmd.o call_dsds.o pack_key.o prod_host.o spawnds.o parse_arg.o names.o du_dir.o shape_call_drms_in.o)

SUMEXE_$(d)	:= $(jpe_$(d)) $(jpeq_$(d)) $(shape_$(d)) $(jsnqry_$(d)) $(jsnxqry_$(d)) $(jsnwqry_$(d))

MODEXE		:= $(MODEXE) $(SUMEXE_$(d)) $(MODEXE_$(d)) $(MODEXEDR_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock) $(MODEXEDR_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))


ALLEXE_$(d)	:= $(MODEXE_$(d)) $(MODEXEDR_$(d)) $(SUMEXE_$(d))
OBJ_$(d)	:= $(ALLEXE_$(d):%=%.o) $(jpe_obj_$(d)) $(jpeq_obj_$(d)) $(shape_obj_$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.o.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(ALLEXE_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(ALLEXE_$(d)) $(MODEXE_SOCK_$(d))

S_$(d)		:= $(notdir $(ALLEXE_$(d)) $(MODEXE_SOCK_$(d)))

$(jpe_$(d)):     $(jpe_obj_$(d))
$(jpeq_$(d)):     $(jpeq_obj_$(d))
$(shape_$(d)):    $(shape_obj_$(d))

# Local rules
$(MODEXEDR_$(d):%=%.o):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/proj/libs/dr

$(SUMEXE_$(d)):		$(LIBSUMSAPI) $(LIBCJSON) $(LIBSUM)

$(MODEXEDR_$(d)) $(MODEXEDR_$(d):%=%_sock):		$(LIBDR)

$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(SUMEXE_$(d)): $(LIBSOIJSOC)
$(SUMEXE_$(d)):		LL_TGT := -lecpg -lpq -lpng -L/home/soi/CM/pvm3/lib/LINUX8664 -lpvm3
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/astro -I/home/soi/CM/pvm3/include -DPVM33 -I$(SRCDIR)/$(d)/../../libs/egsehmicomp
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBASTRO)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

#endif

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
