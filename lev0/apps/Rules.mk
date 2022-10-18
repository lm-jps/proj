# Standard things
sp 		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)

# Local variables
LIBHKLEV0		:= $(d)/libhklev0.a
LIBLIMBFITFXN           := $(d)/liblimbfitfxn.a

test0_$(d)		:= $(addprefix $(d)/, test0)
# wtest_$(d)		:= $(addprefix $(d)/, wtest)
#ingestlev0_$(d)		:= $(addprefix $(d)/, ingest_lev0 ingest_lev0_test ingest_lev0_iris ingest_lev0_irisDUP ingest_lev0_irisDIM ingest_lev0_irisROT ingest_lev0_irisdc ingest_lev0_iris_NEW ingest_lev0_irisAmes decode_dayfile)
#ingestlev0_$(d)		:= $(addprefix $(d)/, ingest_lev0_iris ingest_lev0_irisDUP ingest_lev0_irisROT ingest_lev0_irisdc ingest_lev0_irisAmes)
ingestlev0_$(d)		:= $(addprefix $(d)/, ingest_lev0_irisdc ingest_lev0_irisdc_jim)
ingestlev0orig_$(d)	:= $(addprefix $(d)/, ingest_lev0 decode_dayfile)
#ingestlev1_$(d)		:= $(addprefix $(d)/, build_lev1)
#ingestlev1X_$(d)	:= $(addprefix $(d)/, build_lev1X)
#ingestlev1_mgr_$(d)	:= $(addprefix $(d)/, build_lev1_mgr, build_lev1_mgr_fsn, build_lev1_mgrY)
#ingestlev1_mgr_$(d)	:= $(addprefix $(d)/, build_lev1_mgr, build_lev1_mgr_jim)
ingestlev1_mgr_$(d)	:= $(addprefix $(d)/, build_lev1_mgr)
ingestlev0_obj_$(d)	:= $(addprefix $(d)/, imgdecode_iris.o decode_hk.o load_hk_config_files_iris.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms_iris.o iris_time_setting.o)
ingestlev0orig_obj_$(d) := $(addprefix $(d)/, imgdecode.o decode_hk.o load_hk_config_files.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms.o hmi_time_setting.o set_HMI_mech_values.o)

#xingestlev0_obj_$(d)	:= $(addprefix $(d)/, imgdecode_iris.o decode_hk.o load_hk_config_files.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms_iris.o iris_time_setting.o set_HMI_mech_values.o)
#yingestlev0_obj_$(d)	:= $(addprefix $(d)/, imgdecode_iris.o decode_hk.o  load_hk_config_files.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms_test.o hmi_time_setting_test.o set_HMI_mech_values.o)
#xingestlev0_$(d)	:= $(addprefix $(d)/, xingest_lev0 decode_dayfile)
#xingestlev0_$(d)	:= $(addprefix $(d)/, xingest_lev0)
#yingestlev0_$(d)	:= $(addprefix $(d)/, yingest_lev0)
#LIBHKLEV0_OBJ		:= $(addprefix $(d)/, decode_hk.o load_hk_config_files_iris.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms_iris.o )

LIBHKLEV0_OBJ		:= $(addprefix $(d)/, decode_hk.o load_hk_config_files.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms.o)
LIBLIMBFITFXN_OBJ       := $(addprefix $(d)/, limb_fit_function.o)

buildlev1iris_obj_$(d) 	:= $(addprefix $(d)/, do_flat_iris.o iris_isp2wcs.o)

printtime_$(d)		:= $(addprefix $(d)/, printtime)

#yingestlev0_$(d)	:= $(addprefix $(d)/, yingest_lev0 decode_dayfile)
#yingestlev0_obj_$(d)	:= $(addprefix $(d)/, imgdecode_iris.o decode_hk.o  load_hk_config_files.o decode_hk_vcdu.o save_packet_to_dayfile.o write_hk_to_drms.o hmi_time_setting.o set_HMI_mech_values.o)

#SUMEXE_$(d)	:= $(addprefix $(d)/, ingest_lev0)
#SUMEXE_$(d)	:= $(ingestlev0_$(d)) $(xingestlev0_$(d)) $(ingestlev1_$(d)) $(ingestlev1_mgr_$(d)) $(test0_$(d)) $(wtest_$(d))
#SUMEXE_$(d)	:= $(ingestlev0_$(d)) $(xingestlev0_$(d)) $(ingestlev1_$(d)) $(ingestlev1_mgr_$(d)) $(ingestlev1X_$(d)) $(wtest_$(d))
#SUMEXE_$(d)     := $(ingestlev0_$(d)) $(xingestlev0_$(d)) $(ingestlev1_$(d)) $(ingestlev1_mgr_$(d)) $(wtest_$(d))
SUMEXE_$(d)	:= $(ingestlev0_$(d)) $(ingestlev0orig_$(d)) $(xingestlev0_$(d)) $(yingestlev0_$(d)) $(ingestlev1_$(d)) $(ingestlev1_mgr_$(d)) $(printtime_$(d))
#SUMEXE_$(d)	:= $(ingestlev0_$(d)) $(ingestlev1_$(d)) $(ingestlev1_mgr_$(d)) $(printtime_$(d))
CEXE_$(d)       := $(addprefix $(d)/, fix_hmi_config_file_date)
CEXE		:= $(CEXE) $(CEXE_$(d))
MODEXESUMS_$(d)	:= $(SUMEXE_$(d)) $(PEEXE_$(d))

MODEXE_$(d)	:= $(addprefix $(d)/, convert_fds extract_fds_statev)


# Exclude certain apps from the ia32 build.
ifeq ($(JSOC_MACHINE), linux_ia32)
  BUILDLEV1_$(d)		:=
endif

# Remove from ia32 and gcc builds (since they don't build on ia32 and with gcc)
ifneq ($(JSOC_MACHINE), linux_ia32)
  ifeq ($(COMPILER), icc)
#    BUILDLEV1_$(d)		:=  build_lev1X build_lev1Y build_lev1 build_lev1_fsn
    BUILDLEV1IRIS_$(d)		:= $(addprefix $(d)/, build_lev1_iris build_lev1_iris_test)
#    BUILDLEV1_$(d)		:= $(addprefix $(d)/, build_lev1_aia build_lev1_hmi build_lev1_iris_NEW build_lev1_empty) $(BUILDLEV1IRIS_$(d))
    BUILDLEV1_$(d)		:= $(addprefix $(d)/, build_lev1_aia build_lev1_hmi build_lev1_empty) $(BUILDLEV1IRIS_$(d))
  endif
endif

MODEXE_USEF_$(d)	:= $(addprefix $(d)/, getorbitinfo) $(BUILDLEV1_$(d))
MODEXE_USEF 	:= $(MODEXE_USEF) $(MODEXE_USEF_$(d))

MODEXE		:= $(MODEXE) $(SUMEXE_$(d)) $(MODEXE_$(d)) $(MODEXEDR_$(d)) $(PEEXE_$(d))

MODEXE_SOCK_$(d):= $(MODEXE_$(d):%=%_sock) $(MODEXEDR_$(d):%=%_sock)
MODEXE_SOCK	:= $(MODEXE_SOCK) $(MODEXE_SOCK_$(d))

ALLEXE_$(d)	:= $(MODEXE_$(d)) $(MODEXEDR_$(d)) $(MODEXE_USEF_$(d)) $(SUMEXE_$(d)) $(CEXE_$(d))
#OBJ_$(d)	:= $(ALLEXE_$(d):%=%.o)
OBJ_$(d)	:= $(ALLEXE_$(d):%=%.o) $(TESTEXE_USEF_$(d):%=%.o) $(ingestlev0_obj_$(d)) $(ingestlev0orig_obj_$(d))
DEP_$(d)	:= $(OBJ_$(d):%=%.d)
CLEAN		:= $(CLEAN) \
		   $(OBJ_$(d)) \
		   $(ALLEXE_$(d)) \
		   $(TESTEXE_USEF_$(d)) \
		   $(MODEXE_SOCK_$(d))\
		   $(LIBHKLEV0)\
                   $(LIBLIMBFITFXN)\
		   $(DEP_$(d))

TGT_BIN	        := $(TGT_BIN) $(ALLEXE_$(d)) $(MODEXE_SOCK_$(d))
TGT_LIB		:= $(TGT_LIB) $(LIBHKLEV0) $(LIBLIMBFITFXN)

S_$(d)		:= $(notdir $(ALLEXE_$(d)) $(TESTEXE_USEF_$(d)) $(MODEXE_SOCK_$(d)))

$(ingestlev0_$(d)):	$(ingestlev0_obj_$(d))
$(ingestlev0orig_$(d)):	$(ingestlev0orig_obj_$(d))
$(xingestlev0_$(d)):	$(xingestlev0_obj_$(d))
#$(yingestlev0_$(d)):	$(yingestlev0_obj_$(d))
$(BUILDLEV1IRIS_$(d)):	$(buildlev1iris_obj_$(d))

# Local rules
$(MODEXESUMS_$(d)):	$(LIBSUMSAPI) $(LIBCJSON) $(LIBSUM)

$(MODEXEDR_$(d):%=%.o):		CF_TGT := $(CF_TGT) -I$(SRCDIR)/proj/libs/dr

$(OBJ_$(d)):		$(SRCDIR)/$(d)/Rules.mk
ifeq ($(JSOC_MACHINE), linux_avx2)
$(SUMEXE_$(d)):		LL_TGT := -L  $(POSTGRES_LIBS) -lecpg -lpq -lpng -lmkl_rt
else
$(SUMEXE_$(d)):		LL_TGT := -L  $(POSTGRES_LIBS) -lecpg -lpq -lpng $(FFTW3LIBS)
endif

$(MODEXEDR_$(d)) $(MODEXEDR_$(d):%=%_sock):		$(LIBDR)

ifeq ($(COMPILER), icc)
   ifeq ($(JSOC_MACHINE), linux_avx2)
     MKL     := -lmkl_rt
   else
     MKL     := -static-intel -lmkl_em64t
   endif
endif

#$(SUMEXE_$(d)):		LL_TGT := -L/home/production/cvs/jsoc/lib/saved/$(JSOC_MACHINE) -lhmicomp_egse -lecpg -lpq -lpng -L/SGE/lib/lx24-amd64/ -ldrmaa -Wl,-rpath,/SGE/lib/lx24-amd64

$(PEEXE_$(d)):		LL_TGT := -L $(POSTGRES_LIBS) -lecpg -lpq
$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I$(SRCDIR)/$(d)/../../libs/interpolate/ -I$(SRCDIR)/$(d)/../../libs/astro -I$(SRCDIR)/$(d)/../../libs/egsehmicomp $(FFTWH) -DLEV0SLOP

#$(OBJ_$(d)):		CF_TGT := $(CF_TGT) -DCDIR="\"$(SRCDIR)/$(d)\"" -I/home/jsoc/cvs/JSOC/proj/libs/interpolate/ -I$(SRCDIR)/$(d)/../../libs/astro -I/home/jsoc/include

ifeq ($(COMPILER), icc)
   ifneq ($(JSOC_MACHINE), linux_ia32)
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(TESTEXE_USEF_$(d)): LL_TGT := $(LL_TGT) $(MKL)
   endif
endif

# do not use $(LIBASTRO) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
# do not use $(LIBINTERP) since we can't be sure if its Rules.mk, which is where
# this variable gets set, has been read yet
$(MODEXE_$(d)) $(MODEXE_SOCK_$(d)) $(MODEXE_USEF_$(d)) $(TESTEXE_USEF_$(d)):	proj/libs/astro/libastro.a proj/libs/interpolate/libinterp.a

# decode_hk.c and load_hk_config_files.c both use egsehmicomp.h header (but not libesgehmicomp.a)
$(LIBHKLEV0_OBJ):	CF_TGT := $(CF_TGT) -I$(SRCDIR)/$(d)/../../libs/egsehmicomp
$(LIBHKLEV0):		$(LIBHKLEV0_OBJ)
			$(ARCHIVE)
			$(SLLIB)


$(LIBLIMBFITFXN_OBJ):   CF_TGT := $(CF_TGT) $(FFTWH) -I$(SRCDIR)/$(d)/../../libs/interpolate
ifeq ($(JSOC_MACHINE), linux_avx2)
$(LIBLIMBFITFXN):       LL_TGT := $(LL_TGT) -lmkl_rt
else
$(LIBLIMBFITFXN):       LL_TGT := $(LL_TGT) $(FFTW3LIBS)
endif
$(LIBLIMBFITFXN):       $(LIBLIMBFITFXN_OBJ)
			$(ARCHIVE)
			$(SLLIB)

$(ingestlev0_$(d)):	LL_TGT := $(LL_TGT) -lpng
$(ingestlev0orig_$(d)):	LL_TGT := $(LL_TGT) -lpng
$(xingestlev0_$(d)):	LL_TGT := $(LL_TGT) -lpng
#$(yingestlev0_$(d)):	LL_TGT := $(LL_TGT) -lpng
$(ingestlev1_$(d)):	LL_TGT := $(LL_TGT)
$(ingestlev1X_$(d)):	LL_TGT := $(LL_TGT)
$(ingestlev1_mgr_$(d)):	LL_TGT := $(LL_TGT)
$(test0_$(d)):	LL_TGT := $(LL_TGT)
$(printtime_$(d)): LL_TGT := $(LL_TGT)
#$(wtest_$(d)):	LL_TGT := $(LL_TGT)

# Shortcuts
.PHONY:	$(S_$(d))
$(S_$(d)):	%:	$(d)/%

# Standard things
-include	$(DEP_$(d))

d		:= $(dirstack_$(sp))
sp		:= $(basename $(sp))
