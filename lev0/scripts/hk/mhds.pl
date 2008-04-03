#!/usr/bin/perl
#############################################################################
#NAME:   mhds.pl   - Make Housekeeping Data Series                          #
#Author: Carl                                                               # 
#Date:   March, 20, 2008 moved to JSOC                                      #
#Execution: mhds.pl stha=<stanford file version #> gtcids=<ground file ver#>#
#Example  : mhds.pl stha=1.140  gtcids=1.132                                #
#############################################################################
# main program  
##get environment variables and initialize variables.
$hm=$ENV{'HOME'};
$ENV{'HK_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file";
$ENV{'HK_JSVN_MAP_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsn_map_file";

#common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";
$ENV{'CVSROOT'}=":ext:sunroom.stanford.edu:/home/cvsuser/cvsroot";
$ENV{'CVS_BINARY_ROOT'}="/home/cvsuser/cvs_binary_root";
#$ENV{'CVS_RSH'}="/usr/bin/rsh";
$ENV{'CVS_RSH'}="ssh";

#list of apid to do use for display only-see real setting in cjds.pl
$apids_to_do=$ENV{'HK_MHDS_APID_LIST'}="0017 0019 0021 0029 0445 0475 0529 0569";

#Log machine name and environment variables
print "-->Starting running mhds.pl\n";
print "-->Environment in mhds : \n" ;
foreach $key (sort keys(%ENV)) 
{
  ##print "$key = $ENV{$key}\n";
}

#check for any arguments passed in command
&check_arguments();

#get stha and gtcids file
$new_stha=substr($ARGV[0],5) ;
$new_gtcids=substr($ARGV[1],7);

#check if have hkpdf file for file version number
$fvn= $new_stha;

#log status#
print "-->Files to process are : $ARGV[0]  and $ARGV[1] \n";
print "-->List of APIDS: $ENV{'HK_MHDS_APID_LIST'} \n";
print "-->File Version number is < $fvn >\n";
print "-->Start Run...\nIf needed, new data series will be made for these APIDS: $ENV{'HK_MHDS_APID_LIST'}\n";

# make preliminary jsd files for current version
print  "-->Executing  <make_jsd_file.pl prelim $fvn>\n";
print  "-->Executing <make_jsd_file.pl prelim $fvn> to create preliminary JSD files for file version $fvn\n";
$retlog=`perl make_jsd_file.pl prelim $fvn`;
print  "-->Execution log from make_jsd_file.pl prelim <fvn>\n$retlog\n";

# make all JSVN to PVN map files based on apids in 
# current file version number folder  
print "-->Executing  <do_jsvn_map_file.pl -g $fvn>\n";
print "-->Executing <do_jsvn_map_file.pl -g $fvn> to create latest JSOC Version Number Map Files\n";
$retlog=`perl do_jsvn_map_file.pl -g $fvn`;
print "-->Execution log from do_jsvn_map_file.pl -g <fvn>\n$retlog\n";

# make final jsd files for current version
print "-->Executing  <make_jsd_file.pl final $fvn>\n";
print "-->Executing <make_jsd_file.pl final $fvn> to create final JSD files for file version $fvn\n";
$retlog=`perl make_jsd_file.pl final $fvn`;
print "-->Execution log from make_jsd_file.pl final <$fvn>\n$retlog\n";


# check into cvs the new JSD files created when running make_jsd_file.pl final
# check new jsd files for new ones
&checkin_final_jsd();

# check new jsd files for new ones
&checkin_jsvn_maps();

# check if list of new JSD files is created in make_jsd_file.pl.
# if new JSD, then if on list of APIDS to do, create data series
# using project name and data type name specified in script
print "-->Executing <cjds.pl> to create keyword data series using new JSD file\n";
$retlog=`perl cjds.pl`;
print "-->Execution log for from cjds.pl:\n$retlog";

#close 
print "-->Finished running mhds.pl\n";
#############################################################################
# subroutine check arguments and set flags                                  #
#############################################################################
sub check_arguments()
{
  $help_flg= "0";
  $stha_flg="0";
  $gtcids_flg="0";
  
  if ($#ARGV < 0)
  {
    $help_flg="1";
  }
    
  if ($#ARGV >= 0)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    elsif (substr($ARGV[0],0,5) eq "stha=" )
    {
      #use file to get list of apids to create map files for
      $stha_flg = "1";
    }
  }
  if ($#ARGV >= 1)
  {
    if ( substr($ARGV[1],0,7) eq "gtcids=") 
    {
      #use file to get list of apids to create map files for
      $gtcids_flg = "1";
    }
    else
    {
      print "Warning:arguments entered not correct\n";
    }
  }

  if ( $help_flg eq "1")
  {
     &show_help_info;
  }
}
#############################################################################
# subroutine show_help_info: show help information                          #
#############################################################################
sub show_help_info
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Execute to make new keyword data series based on creation of new JSD file:\n
             mhds.pl stha=<STANFORD file version number> gtcids=<GROUND file version number>:\n";
  print "             Example:mhds.pl stha=1.140 gtcids=1.132\n\n";
  print "(1c)Help Information: mhds.pl -h \n"; 
  print "             Example: mhds.pl -H or mhds.pl -h \n\n";
  print "(2) General Description: Used to run a series of scripts. This includes \n";
  print "                         the tasks of creating preliminary jsd files, jsoc \n";
  print "                         version number map files, final jsd files, and based \n";
  print "                         any new jsd files create an new JSOC data series in DRMS.\n";
  print "                         This script was create to run in the clmq.pl script after\n";
  print "                         after the execution of the cicf.pl script.\n";
  print "(3) Verify before running the setting of HK_MHDS_APID_LIST in mhds.pl and cjds.pl scripts.\n";
  print "    The setting here(0019 0029 etc) tells script which APIDS to make new data series for.\n";
  print "(4) Verify before running the setting of project and data type name in cjds.pl script.\n";
  print "    The setting tells script what to use when creating data series name.\n";
  print "    For example, project name, hmi_ground and data type name, lev0 will use these.\n";
  print "    values to create series names: hmi_ground.lev_<APID VALUE>_<JVN VALUE>.\n";
 
  print "(5) Note the new JSD files created are in file, new_jsd_files.txt. \n";
  print "    The make_jsd_file script creates items in file. The cjds.pl script \n";
  print "    uses information to either create a new data series in DRMS or ignore item. \n";
  print "(6) Note the enviroment variables are set in the scripts so can run within cron job \n";
  print "(7) Note there is a ssh to n00 to create JSOC Data Series using create_series. \n";
  print "(8) Environment variable  HK_MHDS_APID_LIST \n";
  print "(8a)HK_MHDS_APID_LIST is used to define list of APIDs to data series for in cjds.pl\n";
  print "    and is used in mhds.pl to display only information in log\n";
  print "(8b)HK_MHDS_APID_LIST set to ALL will do create data series for all APIDs when a new jsd file is created.\n";
  print "(8c)HK_MHDS_APID_LIST set to 0029 0021 0019 0021 will create data series in DRMS for these APIDs only.\n";
 
  exit;
}

#################################################
# Check in Final JSD files to CVS               #
#################################################
sub checkin_final_jsd
{
  my($njsds);
  $njsds="$script_dir/new_jsd_files.txt";
  open(NJS, "<$njsds") || die "-->ERROR:Can't Open and Read $njsds: $!\n";
  while (<NJS>) { # regular expression- remove cr from file string
    $_ =~ s/\n//g;

    # add file to cvs
    $lm=`cd  $ENV{'HK_JSD_DIRECTORY'} ; cvs add $_ `;
    print "-->Adding < $_ > to $ENV{'HK_JSD_DIRECTORY'} directory. Log:$lm:\n";

    # commit file to cvs
    $lm=`cd  $ENV{'HK_JSD_DIRECTORY'} ; cvs commit -m \"auto checked into CVS new JSD file\"  $_ `;
    print "-->Commiting < $_ > file in directory $ENV{'HK_JSD_DIRECTORY'}. Log:$lm:\n";
  }
  close NJS;
}

#################################################
#  Checkin JSOC Version Number Map files to CVS #
#################################################
sub checkin_jsvn_maps()
{
  # add map file to cvs
  $lm=`cd  $ENV{'HK_JSVN_MAP_DIRECTORY'} ; cvs add  *-JSVN-TO-PVN`;
  print "-->Adding map files to $ENV{'HK_JSVN_MAP_DIRECTORY'} directory. Log:$lm\n";
                                                                                
  # commit file to cvs
  $lm=`cd  $ENV{'HK_JSVN_MAP_DIRECTORY'} ; cvs commit -m \"auto checked into CVS new map files\"  *-JSVN-TO-PVN `;
  #$lm=`cd  $ENV{'HK_JSVN_MAP_DIRECTORY'} ; pwd ; cvs status *-JSVN-TO-PVN*`;
  print  "-->Commiting map files in directory $ENV{'HK_JSVN_MAP_DIRECTORY'}. Log:$lm:\n";
}
