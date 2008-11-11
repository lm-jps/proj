#!/usr/bin/perl
# gets latest files from cvs for hk config data for lev0 processing.
# used to update dcs0, dcs1, and dcs2.
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.";
$ENV{CVSROOT}=":ext:sunroom.stanford.edu:/home/cvsuser/cvsroot";

# update hmi and aia hk config files
print "..doing cvs update of hmi and aia hk config files\n";
print "..please wait this may take a few minutes\n";
system("cd /home/production/cvs/TBL_JSOC/lev0/hk_config_file; cvs update -d  . ");
print "..completed cvs update of hmi and aia hk config files!\n\n\n";

# update sdo hk config files
print "..doing cvs update of sdo hk config files\n";
system("cd /home/production/cvs/TBL_JSOC/lev0/sdo_hk_config_file; cvs update -d  . ");
print "..completed cvs update of sdo hk config files!\n\n\n";

#update jsoc version number lookup files
print "..doing cvs update of hk jsoc version number lookup files.\n";
system("cd /home/production/cvs/TBL_JSOC/lev0/hk_jsn_map_file/prod; cvs update -d  . ");
print "..completed cvs update of hk jsoc version number lookup files!\n";
