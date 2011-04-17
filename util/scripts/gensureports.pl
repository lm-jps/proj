#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

# All Data
`/home/jsoc/cvs/Development/JSOC/proj/util/scripts/sumstapestat.pl jsoc_sums hmidb 5434 production agg group all > /web/jsoc2/htdocs/sureports/sureportA.txt`;

# AP Data
`/home/jsoc/cvs/Development/JSOC/proj/util/scripts/sumstapestat.pl jsoc_sums hmidb 5434 production agg group ap > /web/jsoc2/htdocs/sureports/sureportB.txt`;
