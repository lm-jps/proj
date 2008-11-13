#!/usr/bin/perl
# Gets latest files from cvs for hk config data, sdo hk config files, and map files for lev0 processing.
# Used to update production account on dcs0, dcs1, dcs2 or any machine with production TBL_JSOC directory.
# Required to be logged in as production. 
# Assumes the path to TBL_JSOC directory is $HOME/cvs/TBL_JSOC/lev0
# Required production has access to cvs. may need to log in with password during run of script

# set variables
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.";
$ENV{CVSROOT}=":ext:sunroom.stanford.edu:/home/cvsuser/cvsroot";
my $hm=$ENV{HOME};
my $lev0_dir="cvs/TBL_JSOC/lev0";

# print info to screen
print "Exiting daily_hk_update.pl script where top home directory=<$hm> and sub lev0 directory=<$lev0_dir>\n";
print "Executing cvs update for HK Configurations files in directory=<$hm/$lev0_dir>\n";

# update hmi and aia hk config files
print "..doing cvs update of hmi and aia hk config files\n";
print "..please wait this may take a few minutes\n";
system("cd $hm/$lev0_dir/hk_config_file; cvs update -d  . ");
print "..completed cvs update of hmi and aia hk config files!\n\n\n";

# update sdo hk config files
print "..doing cvs update of sdo hk config files\n";
system("cd $hm/$lev0_dir/sdo_hk_config_file; cvs update -d  . ");
print "..completed cvs update of sdo hk config files!\n\n\n";

#update jsoc version number lookup files
print "..doing cvs update of hk jsoc version number lookup files.\n";
system("cd $hm/$lev0_dir/hk_jsn_map_file/prod; cvs update -d  . ");
print "..completed cvs update of hk jsoc version number lookup files!\n";
