$(PROJOBJDIR)::
	+@[ -d $@/libs/astro ] || mkdir -p $@/libs/astro
	+@[ -d $@/libs/dr ] || mkdir -p $@/libs/dr
	+@[ -d $@/libs/dsputil ] || mkdir -p $@/libs/dsputil
	+@[ -d $@/libs/gapfiller ] || mkdir -p $@/libs/gapfiller
	+@[ -d $@/libs/interpolate ] || mkdir -p $@/libs/interpolate
	+@[ -d $@/libs/limbcompute_aia ] || mkdir -p $@/libs/limbcompute_aia
	+@[ -d $@/libs/stats ] || mkdir -p $@/libs/stats
	+@[ -d $@/libs/egsehmicomp ] || mkdir -p $@/libs/egsehmicomp
	+@[ -d $@/libs/imrotate ] || mkdir -p $@/libs/imrotate
	+@[ -d $@/libs/projection ] || mkdir -p $@/libs/projection

