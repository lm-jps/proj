#!/usr/bin/perl
##############################################################################
# Name:        do_all_hkpdf.pl   - do all hkpdf files                        #
# Description: Used to create all HKPDF small apid files.                    #
#              or final JSDs.                                                #
# Execution:   (1)To run:      do_all_hkpdf.pl                               #
#              (2)To get help: do_all_hkpdf.pl   -h                          #
# limitation:  Need to move different version of STANFORD files to lmf_dir   #
#              before running.                                               #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on March 27, 2008 #
##############################################################################
# set up environment variables
$hm=$ENV{'HOME'};
$ENV{'HK_CONFIG_FILES'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$script_dir=$ENV{'HK_SCRIPTS_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$lmf_dir="$hm/cvs/TBL_JSOC/lev0/fromlmsal";

if ( $ARGV[0] eq "-h" )
{
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

# old way use push(@stha_list, "STANFORD_TLM_HMI_AIA.txt-1.131");
# create HKPDF small apid files for all file version numbers
foreach $fvn (@file_ver_num_list)
{
  if ($fvn eq "")
  {
    print "breaking. final item in list \n";
    exit;
  }
  # regular expression- remove cr from file string
  $file="STANFORD_TLM_HMI_AIA.txt-$fvn";
  $file =~ s/\n//g;

  # do HKPDF for this file version number 
  print " doing file is <$file>\n";

  # set permission to copy this file version to STANFORD file
  $lm=`cd $lmf_dir; chmod 777 ./STANFORD_TLM_HMI_AIA.txt`;
  print "log is: $lm\n";

  # copy 
  $lm=`cd $lmf_dir;cp ./$file ./STANFORD_TLM_HMI_AIA.txt`;
  print "log is: $lm\n";

  # run script to make HKPDF files 
  #$lm=`cd $script_dir; ./make_hkpdf.pl 4 `;
  # use new script
  $all_hk_value="\"ALL_HK\"";
  $lm=`cd $script_dir; ./make_hkpdf.pl sort=4 apidlist=$all_hk_value `;
  print "log is: $lm\n";
}
######## show_help ##########
sub show_help
{

  print "Help Listing\n";
  print "(1)  Ways to Execute Perl Script: \n";
  print "(1a) Create All HKPDF small APID Files: do_all_hkpdf.pl  \n";
  print "(1b) Get Help Information:              do_all_hkpdf.pl -h \n";
  print "(2)  Environment variable HK_CONFIG_DIR\n";
  print "(2a) Used to get list of file versions number to do to process create all JSD again\n";
  print "(3)  Environment variable HK_SCRIPT_DIR\n";
  print "(3a) Used to find where make_hkpdf.pl script is.\n";
  print "(4)  Variable lmf_dir \n";
  print "(4a) Used to find Stanford files.\n";

}
