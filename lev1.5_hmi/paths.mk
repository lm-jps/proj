$(PROJOBJDIR)::
	+@[ -d $@/lev1.5_hmi/libs/lev15 ] || mkdir -p $@/lev1.5_hmi/libs/lev15
	+@[ -d $@/lev1.5_hmi/apps ] || mkdir -p $@/lev1.5_hmi/apps

