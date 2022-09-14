$(PROJOBJDIR)::
	+@[ -d $@/mag/pfss/apps ] || mkdir -p $@/mag/pfss/apps
	+@[ -d $@/mag/ambig/apps ] || mkdir -p $@/mag/ambig/apps
	+@[ -d $@/mag/harp/apps ] || mkdir -p $@/mag/harp/apps
	+@[ -d $@/mag/harp/libs/matlab/mex/src/util ] || mkdir -p $@/mag/harp/libs/matlab/mex/src/util
	+@[ -d $@/mag/harp/libs/matlab/mex/src/mex2c ] || mkdir -p $@/mag/harp/libs/matlab/mex/src/mex2c
	+@[ -d $@/mag/harp/libs/matlab/mfile-mex/standalone ] || mkdir -p $@/mag/harp/libs/matlab/mfile-mex/standalone
	+@[ -d $@/mag/harp/libs/matlab/mfile-mex/assignment ] || mkdir -p $@/mag/harp/libs/matlab/mfile-mex/assignment
	+@[ -d $@/mag/harp/libs/matlab/mfile-mex/fits ] || mkdir -p $@/mag/harp/libs/matlab/mfile-mex/fits
	+@[ -d $@/mag/harp/libs/matlab/mfile-mex/hmi-mask-patch ] || mkdir -p $@/mag/harp/libs/matlab/mfile-mex/hmi-mask-patch
	+@[ -d $@/mag/ident/apps ] || mkdir -p $@/mag/ident/apps
	+@[ -d $@/mag/ident/libs/mex2c ] || mkdir -p $@/mag/ident/libs/mex2c
	+@[ -d $@/mag/ident/libs/mexfunctions ] || mkdir -p $@/mag/ident/libs/mexfunctions
	+@[ -d $@/mag/ident/libs/util ] || mkdir -p $@/mag/ident/libs/util
	+@[ -d $@/mag/libs ] || mkdir -p $@/mag/libs
	+@[ -d $@/mag/libs/util ] || mkdir -p $@/mag/libs/util
	+@[ -d $@/mag/patch/apps ] || mkdir -p $@/mag/patch/apps
	+@[ -d $@/mag/nlfff/apps ] || mkdir -p $@/mag/nlfff/apps
	+@[ -d $@/mag/d4vm/apps ] || mkdir -p $@/mag/d4vm/apps
	+@[ -d $@/mag/remapmags/apps ] || mkdir -p $@/mag/remapmags/apps
	+@[ -d $@/mag/synop/apps ] || mkdir -p $@/mag/synop/apps
	+@[ -d $@/mag/polarfield/apps ] || mkdir -p $@/mag/polarfield/apps
	+@[ -d $@/mag/polcorr/apps ] || mkdir -p $@/mag/polcorr/apps

