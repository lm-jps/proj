$(PROJOBJDIR)::
	+@[ -d $@/cgem/apps ] || mkdir -p $@/cgem/apps
	+@[ -d $@/cgem/lorentz/apps ] || mkdir -p $@/cgem/lorentz/apps
	+@[ -d $@/cgem/pdfi/apps ] || mkdir -p $@/cgem/pdfi/apps
	+@[ -d $@/cgem/prep/apps ] || mkdir -p $@/cgem/prep/apps
