# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
# NOTE: Add the base of the module's filename below (next to mymod)
CLEANERMODEXE_$(d)	:= $(addprefix $(d)/, query_engine)

MODEXE_$(d)	:= $(addprefix $(d)/, ingest_mharp ingest_mharp_log track_hmi_check_masks)
MODEXE		:= $(MODEXE) $(MODEXE_$(d)) $(CLEANERMODEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

CLEANERMODEXE_SOCK_$(d)	:=  $(CLEANERMODEXE_$(d):%=%_sock)

EXE_$(d)	:= $(MODEXE_$(d)) $(CLEANERMODEXE_$(d)) $(CLEANERMODEXE_SOCK_$(d)) $(MODEXE_SOCK_$(d)) 
OBJ_$(d)	:= $(MODEXE_$(d):%=%.o) $(CLEANERMODEXE_$(d):%=%.o)
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(EXE_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(EXE_$(d)) $(CLEANERMODEXE_$(d))

S_$(d)		:= $(notdir $(EXE_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" $(CFITSIOH) -I$(SRCDIR)/$(d)/../../libs/json

MODLIBS_SOCK_$(d)	:=  $(LIBJSOC_MAIN_SOCK) $(LIBDRMSCLIENT) $(LIBDEFSCLIENT) $(LIBDBCLIENT) $(LIBSUMSAPI) $(LIBTHREADUTIL) $(LIBRICECOMP) $(LIBCMDPARAMS) $(LIBTIMEIO) $(LIBFITSRW) $(LIBERRLOG) $(LIBEXPDRMS) $(LIBEXPUTL) $(LIBMISC) $(LIBDSTRUCT) $(LIBSTATS)

$(d)/query_engine_sock.o:	CF_TGT := $(CF_TGT) -DDRMS_CLIENT
$(d)/query_engine_sock.o:	$(SRCDIR)/$(d)/query_engine.c
				$(COMP)
$(d)/query_engine_sock:	LL_TGT := $(LL_TGT) $(PGLIBS) $(CFITSIOLIBS)
$(d)/query_engine_sock:	$(d)/query_engine_sock.o $(MODLIBS_SOCK_$(d))
			$(LINK)
			$(SLBIN)

# NOTE: Add dependent libraries with the -I compiler flag, and make the module depend
#   on that library
#$(OBJ_$(d)):				CF_TGT := -I$(SRCDIR)/$(d)/../../libs/somelib
#$(OBJ_$(d)):	CF_TGT := -I$(SRCDIR)/$(d)/../../libs/somelib
# $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)):	$(LIBSOMELIB)

$(EXE_$(d)):  $(LIBJSON)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
