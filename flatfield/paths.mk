$(PROJOBJDIR)::
	+@[ -d $@/flatfield/libs/flatfieldlib ] || mkdir -p $@/flatfield/libs/flatfieldlib
	+@[ -d $@/flatfield/apps ] || mkdir -p $@/flatfield/apps
	+@[ -d $@/flatfield/off_flat_IDL ] || mkdir -p $@/flatfield/off_flat_IDL

