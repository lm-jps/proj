#!/bin/csh
setenv TZ UTC
source /home/jsoc/.setJSOCenv
setenv SUMSERVER k1
echo `date`
echo $1
/home/jsoc/cvs/Development/JSOC/bin/linux_avx/aia_lev1p5  dsinp=aia.lev1_nrt2\[$1\]\[\?quality+1\>0\?\] dsout=aia_test.lev1p5 rescale=1 scale_to=0.6 do_stretchmarks=1 with_keys=1
