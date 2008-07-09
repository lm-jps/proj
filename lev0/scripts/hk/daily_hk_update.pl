#!/usr/bin/perl
# gets latest files from cvs for hk config data for lev0 processing.
# used to update dcs0, dcs1, and dcs2.
$ENV{JSOCROOT}="/home/production/cvs/JSOC";
$ENV{CVSROOT}=":ext:sunroom.stanford.edu:/home/cvsuser/cvsroot";
system("cd /home/production/cvs/TBL_JSOC/lev0/hk_config_file; cvs update -d  . ");

