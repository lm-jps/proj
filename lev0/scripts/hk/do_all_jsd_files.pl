#!/usr/bin/perl
##############################################################################
# Name:        do_all_jsd_files.pl - do all jsd files                        #
# Description: Used to create all preliminary and final JSD file if needed.  #
#              Uses the file version numbers in HK_CONFIG_DIRECTORY to get   #
#              list of file version to do. By using argument can create      #
#              preliminary or final JSDs.                                    #
# Execution:   (1)To run to create prelim jsd's: do_all_jsd_file.pl prelim   #
#              (2)To run to create final jsd's:  do_all_jsd_file.pl final    #
#              (3)For help: do_all_jsd_file.pl -h                            #
# Limitation:  The limitation is for creating preliminary and final JSD files#
#              is the HK_CONFIG_DIRECTORY variable set to place where all    #
#              file version directories are.                                 #
#              The limitation for creating final files is have JSOC version  #
#              map files created with latest preliminary files.              #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on March 27, 2008 #
##############################################################################
# set up environment variables
$hm=$ENV{'HOME'};

#production settings
$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$ENV{'HK_SH_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/sdo_hk_config_file";
$script_dir=$ENV{'HK_SCRIPTS_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";

#for testing
#$ENV{'HK_SH_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/jsoc-play/sdo_hk_config_file";

#setup contants process type flag
use constant NO_PROCESSING_TYPE =>  0;
use constant HK_PROCESSING_TYPE =>  1;
use constant SDO_PROCESSING_TYPE => 2;

if ($#ARGV < 0 || $#ARGV > 2) 
{
  die "Usage: $0 < prelim | final >  | -h > < SDO | HK >:$!";
}

if ( $ARGV[0] eq "-h" )
{
  &show_help;
  exit;
}
elsif ( $ARGV[0] eq "prelim")
{
  $ffp= 0; #flag for final or prelim value;
  if($ARGV[1] eq "HK")
  {
    $process_type= HK_PROCESSING_TYPE;
  }
  elsif ($ARGV[1] eq "SDO")
  {
    $process_type= SDO_PROCESSING_TYPE;
  }
  else
  {
     print "ERROR:jsoc_do_all_jsd_file.pl: no such processing type\n";
     &show_help;
  }
}
elsif ( $ARGV[0] eq "final")
{
  $ffp= 1;
  if($ARGV[1] eq "HK")
  {
    $process_type=HK_PROCESSING_TYPE;
  }
  elsif ($ARGV[1] eq "SDO")
  {
    $process_type=SDO_PROCESSING_TYPE;
  }
  else
  {
     print "ERROR:jsoc_do_all_jsd_file.pl: no such processing type:$process_type\n";
     &show_help;
  }
}
else
{
  print "Error: use either prelim or final for first argument(i.e.,prelim)\n";
  &show_help;
  exit;
}

if ($process_type == HK_PROCESSING_TYPE)
{
  # create list of files to process
  @file_ver_num_list1=`ls -1 $ENV{'HK_CONFIG_DIRECTORY'} | egrep -x "^1\.[0-9][0-9]" ` ;
  @file_ver_num_list2=`ls -1 $ENV{'HK_CONFIG_DIRECTORY'} | egrep -x "^1\.[0-9][0-9][0-9]" ` ;
  @file_sort_list1 =sort { $a <=> $b } @file_ver_num_list1;
  push (@file_ver_num_list, @file_ver_num_list1);
  @file_sort_list2 =sort { $a <=> $b } @file_ver_num_list2;
  push (@file_ver_num_list, @file_ver_num_list2);
  push (@file_ver_num_list, "");
}
elsif ($process_type == SDO_PROCESSING_TYPE)
{
  # create list of files to process
  @file_ver_num_list0 =`ls -1 $ENV{'HK_SH_CONFIG_DIRECTORY'} | egrep -x "^1\.[0-9]" ` ;
  @file_ver_num_list1=`ls -1 $ENV{'HK_SH_CONFIG_DIRECTORY'} | egrep -x "^1\.[0-9][0-9]" ` ;
  @file_ver_num_list2=`ls -1 $ENV{'HK_SH_CONFIG_DIRECTORY'} | egrep -x "^1\.[0-9][0-9][0-9]" ` ;

  @file_sort_list0 =sort { $a <=> $b } @file_ver_num_list0;
  push (@file_ver_num_list, @file_ver_num_list0);

  @file_sort_list1 =sort { $a <=> $b } @file_ver_num_list1;
  push (@file_ver_num_list, @file_ver_num_list1);

  @file_sort_list2 =sort { $a <=> $b } @file_ver_num_list2;
  push (@file_ver_num_list, @file_ver_num_list2);
  push (@file_ver_num_list, "");
}


# process list of all hk config versions 1.32 to 1.161 or sdo config versions 1.1,1.2,etc.
foreach $fvn (@file_ver_num_list)
{
  $fvn =~ s/\n//g; #regular exp rm cr
  if ($fvn eq "")
  {
    print "breaking. final item in list \n";
    exit;
  }
  if ( $ARGV[0] eq "prelim")
  {
    # print file version number doing
    $prt=sprintf("-->creating preliminary jsd's for file version number:%6.6s\n", substr($fvn,0,6));
    print " $prt";
    if ($process_type == SDO_PROCESSING_TYPE)
    {

      $lm=` $script_dir/jsoc_make_jsd_file.pl prelim $fvn SDO `;
    }
    elsif  ($process_type == HK_PROCESSING_TYPE)
    {
      $lm=` $script_dir/jsoc_make_jsd_file.pl prelim $fvn HK `;
    }
  }
  elsif ( $ARGV[0] eq "final")
  {
    # print file version number doing
    $prt=sprintf("-->creating final jsd's for file version number:%6.6s\n", substr($fvn,0,6));
    print " $prt";
    if ($process_type == SDO_PROCESSING_TYPE)
    {
      $lm=` $script_dir/jsoc_make_jsd_file.pl final $fvn SDO `;
    }
    elsif  ($process_type == HK_PROCESSING_TYPE)
    {
      $lm=` $script_dir/jsoc_make_jsd_file.pl final $fvn HK `;
    }
  }
  else
  {
     print "Error could not create preliminary or final JSD - existing \n";
  }

}

######## show_help ##########
sub show_help
{

  print "\nHelp Listing\n";
  print "(1)  Ways to Execute Perl Script: \n";
  print "(1a) Create All Preliminary SDO or HK JSD Files: do_all_jsd_file.pl prelim <SDO | HK >\n";
  print "(1b) Create All Final SDO or HK JSD Files:       do_all_jsd_file.pl final < SDO | HK >\n";
  print "     Where file type flag is either prelim or final.\n";
  print "     This means to create either preliminary or final JSD files.\n";
  print "     Where processing type flag is either SDO or HK. \n";
  print "     This means to process HK data in Stanford file or SDO data from Goddard file. \n";
  print "(1c) Get Help Information:                       do_all_jsd_file.pl -h \n";
  print "(2)  Environment variable HK_CONFIG_DIRECTORY\n";
  print "(2a) Used to get list of file versions number to do to process create all JSD again HK JSDs\n";
  print "(3)  Environment variable HK_SH_CONFIG_DIRECTORY\n";
  print "(2a) Used to get list of file versions number to do to process create all JSD again for SDO JSDs\n";
}
