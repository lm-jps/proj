#!/usr/bin/perl
##############################################################################
# Name:        do_all_jsd_files.pl - do all jsd files                        #
# Description: Used to create all preliminary and final JSD file if needed.  #
#              Uses the file version numbers in HK_CONFIG_FILES to get list  #
#              of file version to do. By using argument can create prelimary #
#              or final JSDs.                                                #
# Execution:   (1)To run to create prelim jsd's: do_all_jsd_file.pl prelim   #
#              (2)To run to create final jsd's:  do_all_jsd_file.pl final    #
#              (3)For help: do_all_jsd_file.pl -h                            #
# Limitation:  The limitation is for creating preliminary and final JSD files#
#              is the HK_CONFIG_FILES variable set to place where all file   #
#              version directories are.                                      #
#              The limitation for creating final files is have JSOC version  #
#              map files created with latest preliminary files.              #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on March 27, 2008 #
##############################################################################
# set up environment variables
$hm=$ENV{'HOME'};
$ENV{'HK_CONFIG_FILES'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$script_dir=$ENV{'HK_SCRIPTS_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";

if ($#ARGV < 0 || $#ARGV > 1) 
{
  die "Usage: $0 < prelim | final >  | -h > :$!";
}

if ( $ARGV[0] eq "-h" )
{
  &show_help;
  exit;
}
elsif ( $ARGV[0] eq "prelim")
{
  $file_version_number = $ARGV[1];
  $ffp= 0; #flag for final or prelim value;
}
elsif ( $ARGV[0] eq "final")
{
  $ffp= 1;
}
else
{
  print "Error: use either prelim or final for first argument(i.e.,prelim)\n";
  &show_help;
  exit;
}

# create list of files to process
@file_ver_num_list1=`ls -1 $ENV{'HK_CONFIG_FILES'} | egrep -x "^1\.[0-9][0-9]" ` ;
@file_ver_num_list2=`ls -1 $ENV{'HK_CONFIG_FILES'} | egrep -x "^1\.[0-9][0-9][0-9]" ` ;
@file_sort_list1 =sort { $a <=> $b } @file_ver_num_list1;
push (@file_ver_num_list, @file_ver_num_list1);
@file_sort_list2 =sort { $a <=> $b } @file_ver_num_list2;
push (@file_ver_num_list, @file_ver_num_list2);
push (@file_ver_num_list, "");

# process list of all hk config versions 1.32 to 1.161
foreach $fvn (@file_ver_num_list)
{
  if ($fvn eq "")
  {
    print "breaking. final item in list \n";
    exit;
  }
  if ( $ARGV[0] eq "prelim")
  {
    # print file version number doing
    $prt=sprintf("creating preliminary jsd's for file version number:%6.6s\n", substr($fvn,0,6));
    print " $prt";
    $lm=`cd $script_dir; ./make_jsd_file.pl prelim $fvn `;
    print "$lm\n";
  }
  elsif ( $ARGV[0] eq "final")
  {
    # print file version number doing
    $prt=sprintf("creating final jsd's for file version number:%6.6s\n", substr($fvn,0,6));
    print " $prt";
    $lm=`cd $script_dir; ./make_jsd_file.pl final $fvn `;
    print "$lm\n";
  }
  else
  {
     print "Error could not create preliminary or final JSD - existing \n";
  }

}

######## show_help ##########
sub show_help
{

  print "Help Listing\n";
  print "(1)  Ways to Execute Perl Script: \n";
  print "(1a) Create All Preliminary JSD Files: do_all_jsd_file.pl prelim \n";
  print "(1b) Create All Final JSD Files:       do_all_jsd_file.pl final\n";
  print "(1c) Get Help Information:             do_all_jsd_file.pl -h \n";
  print "(2)  Environment variable HK_CONFIG_DIRECTORY\n";
  print "(2a) Used to get list of file versions number to do to process create all JSD again\n";
}
