# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables

# There are two types of modules: direct-connect modules and socket-connect modules.
# To create a direct-connect module, place the name of the executable in MODEXE.
# To create a socket-connect module, place the name of the executable (with a
# _sock suffix) in MODEXE_SOCK. In the examples below, scalesegments, opendsdsrecs, 
# and demo_td08062007 are socket-connect modules, but there are both direct- and 
# socket-connect versions of hello_world.  Modules created from .f files are always
# created as sockect-connect modules (and placed in the FMODEXE variable, ie., 
# helloworld, xinterp).  
#
# Socket-connect modules can share a single DRMS session, so they are ideal when
# transferring data among them without having to save intermediate data as series.
# Direct-connect modules cannot share data without first saving them to a series, 
# but they can also operate in "read-only" mode.  They will not start up SUMS,
# unless there is an explicit call to do so.  They are ideal for situations like
# calling from a web-interface when you do not want users making changes.
#
# For most science processing, you would want to use socket-connect modules to
# take advantage of the ability to share data in a pipeline.

MODEXE_$(d)	:= $(addprefix $(d)/, hello_world threadsigs threadalrm)
MODEXE		:= $(MODEXE) $(MODEXE_$(d))

MODEXE_SOCK_$(d):= $(addsuffix _sock, $(addprefix $(d)/, hello_world opendsdsrecs))
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

# C module that uses a third-party Fortran library
MODEXE_USEF_SOCK_$(d)	:= $(addsuffix _sock, $(addprefix $(d)/, demo_td08062007 scalesegments))
MODEXE_USEF_SOCK	:= $(MODEXE_USEF_SOCK) $(MODEXE_USEF_SOCK_$(d))

# F modules (DoIt() resides in .f file)
ifeq ($(COMPILER), icc)
  FMODEXE_$(d)	:= $(addsuffix _sock, $(addprefix $(d)/, helloworld xinterp ringfit_ssw))
  FMODEXE_GONG_$(d)	:= $(addsuffix _sock, $(addprefix $(d)/, f_ingest_gong_mrv f_dup_gong_mrv))
  GONGLIB_$(d)	:= $(d)/ingest_gong_mrv_lib.o
  FMODEXE	:= $(FMODEXE) $(FMODEXE_$(d)) $(FMODEXE_GONG_$(d))
endif

OBJ_$(d)	:= $(MODEXE_$(d):%=%.o) $(MODEXE_SOCK_$(d):%_sock=%.o) $(GONGLIB_$(d))

# .o files that use the third-party Fortran libraries
OBJUSEF_$(d)	:= $(MODEXE_USEF_SOCK_$(d):%_sock=%.o)

# .o files that depend on fdrms.mod - the interface to the Fortran versions of DRMS
OBJF_$(d)	:= $(FMODEXE_$(d):%_sock=%.o) $(FMODEXE_GONG_$(d):%_sock=%.o)


DEP_$(d)	:= $(OBJ_$(d):%=%.d) $(OBJUSEF_$(d):%=%.d) $(OBJF_$(d):%=%.d)

CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(OBJF_$(d)) \
		   $(MODEXE_$(d)) \
		   $(MODEXE_USEF_SOCK_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(FMODEXE_$(d) \
		   $(FMODEXE_GONG_$(d)) \
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) 

EXAMPLES	:=  $(MODEXE_$(d)) \
		    $(MODEXE_SOCK_$(d)) \
		    $(MODEXE_USEF_SOCK_$(d)) \
		    $(FMODEXE_$(d)) \
		    $(FMODEXE_GONG_$(d))

S_$(d)		:= $(notdir $(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_SOCK_$(d)) $(FMODEXE_$(d)) $(FMODEXE_GONG_$(d)))

# Local rules
$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJF_$(d)):		$(SRCDIR)/$(d)/Rules.mk
$(OBJF_$(d)):		CF_TGT := $(CF_TGT) $(FMATHLIBSH) -DCDIR="\"$(SRCDIR)/$(d)\""
$(OBJUSEF_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(OBJUSEF_$(d)):	CF_TGT := $(CF_TGT) $(FMATHLIBSH) -DCDIR="\"$(SRCDIR)/$(d)\""

# Don't use the make variable FDRMSMODOBJ since you really don't know when it will be 
# evaluated. Specify fdrms.o relative to root.
$(OBJF_$(d)):		base/drms/libs/api/client/fdrms.o
$(GONGLIB_$(d)):	$(SRCDIR)/$(d)/Rules.mk
$(GONGLIB_$(d)):	CF_TGT := $(CF_TGT) $(FMATHLIBSH)

$(FMODEXE_$(d)):	FF_TGT := -module base/drms/libs/api/client
$(FMODEXE_GONG_$(d)):	FF_TGT := -module base/drms/libs/api/client
$(FMODEXE_GONG_$(d)):	$(GONGLIB_$(d))

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
