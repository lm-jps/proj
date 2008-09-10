#!/usr/bin/perl
#############################################################################
# Name:        cjds.pl - Create jsd data series                             #
# Description: Creates JSD if have a new JSD version for HK by APID JSD's   #
# Execution:   (1) To run :                                                 #
# Execution:   cjds.pl  <list of apids to do | ALL >                        # 
# Example Exc: cjds.pl  "0129 0445 0529 0475 0529"                          # 
# Limitation:  The script needs these environment variables set properly:   #
#              HK_JSD_DIRECTORY                                             #
# NOTES:       Need to add help information.                                #
# Author:      Carl                                                         #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008#
#############################################################################
# main program  
# Set environment variables and variables
$hm=$ENV{'HOME'};
$exe_dir=$ENV{'DF_DRMS_EXECUTABLES'}="$hm/cvs/JSOC/bin/linux_x86_64";
$prod_jsd_dir=$ENV{'HK_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file/prod";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";

# new jsd files created in text file.
$njsds=$ENV{'HK_NEW_JSD_FILE'}="$script_dir/new_jsd_files.txt";

# Common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$exe_dir}";
$ENV{'CVS_RSH'}="ssh";
#$ENV{'CVS_RSH'}="/usr/bin/rsh";

# check if do help
&check_arguments();

# check args passed
print "cjds:APIDs list to do is in ARGV[0] is $ARGV[0]\n";
$apids_to_do=$ARGV[0];

# set environment variables for cron job executing this script based on machine
$account=`whoami`;
chop $account;
if ( $account eq "carl")
{
  print "-->Found <$account> account\n";
  $ENV{'SSH_AUTH_SOCK'}="/tmp/ssh-XXSCh2du/agent.1132";
  $ENV{'SSH_ASKPASS'}="/usr/libexec/openssh/gnome-ssh-askpass";
}
elsif ( $account eq "production")
{
   $ENV{'SSH_ASKPASS'}="/usr/libexec/openssh/gnome-ssh-askpass";
   print LF "-->working in <$account> account\n";
   print LF "-->SSH setting for SSH_ASKPASS is <$ENV{'SSH_ASKPASS'}>\n";
}
else
{
   print "ERROR: did not find this account name <$account>- Need to update cjds.pl script.\n";
   print `date`;
   die "ERROR: did not find this account name <$account>- Need to update cjsds.pl script.";
}

# log message to standard out
print "-->Executing script cjds.pl to create JSD data series if needed\n";

# backup list of new jsd's created in case want to re-run
if ( -s $njsds )
{
  my $d=`date`;
  #regular expression to add - were there are blanks
  $d =~ s/^| /-/g;
  $lm=`cp $njsds $script_dir/new_jsd_files.txt-$d`;
}

#open file containing list of new jsd files created
open(NJS, "$njsds") || die "-->ERROR:Can't Open $njsds: $!\n";

# push new jsd created in list
while (<NJS>)
{
  chop $_;
  push(@new_jsd_list, $_) ;
}

#log message to standard out
if ( -s $njsds )
{
  print "-->List new jsd files created is:\n@new_jsd_list\n";
}
else
{
  print "-->There was no new jsd files created. Therefore no new data series created.\n";
}
# check if new jsd files created are in list of apids to do
foreach $f (@new_jsd_list)
{
  #get apid value of jsd file
  @s_line=split('_',$f );

  # get list of apids to do
  @s_list=split(' ',$apids_to_do);
  print "-->apid to do list is $apids_to_do\n";

  # for each new jsd file check if in list of apids to do
  foreach $atd  (@s_list)
  {
    # if could not find  apid in jsd file name check next
    if ((( index  $s_line[1],  $atd ) == -1) and (index  $atd,  "ALL") == -1) 
    {
      print "-->did not find <$s_line[1]> when checked $atd\n";
    }
    else
    {
      # log information to standard out
      print "-->Found <$s_line[1]> in $apids_to_do\n";
      print "-->new ds name and data series to create will be < $f >\n";

      # chmod data series file 
      $lm=` cd  $prod_jsd_dir; chmod 666  $f `;
      print "-->Finish chmod of $prod_jsd_dir and $f\nlog after chmod command is $lm\n";

      # create series using create_series in  JSOC environment 
      $lm= `cd $prod_jsd_dir; $ENV{'DF_DRMS_EXECUTABLES'}/create_series  $f `;
      print "-->Created series using jsd file at:<$prod_jsd_dir/$f>\n-->log from create_series command is $lm\n";

      # get next jsd in new jsd file list
      last;
    }
  }
}
print "-->Finished running cjds.pl - create jsoc data series\n";
close NJS;
## clear file
open(NJS, ">$njsds") || die "-->ERROR:Can't Open $njsds: $!\n";
close NJS;

###################
# check_arguments #
###################
sub check_arguments
{
  if ( $ARGV[0] eq "-h" )
  {
    &show_help_info;
  }
  elsif ($#ARGV != 0)
  {
    print "ERROR:cjds.pl: Need to enter argument for apids to do.\n\n";
    &show_help_info;
  }
}
#################
# show_help_info#
#################
sub show_help_info
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Create Data Series in DRMS based on item in new jsd file:\n";
  print "    cjds.pl < list of decimal 4 digit apids | ALL > \n";
  print "    Where list of decimal 4 digit apids are in format: <dddd><space><dddd><space>... \n";
  print "    Where ALL mean create every data series found in new_jsd_file.txt file \n";
  print "    Example: cjds.pl  ALL \n";
  print "    Example: cjds.pl  \"0129 0445 0475\" \n";
  print "(1b)Get Help Information: cjds.pl -h\n";
  print "(2)Environment variable HK_JSD_DIRECTORY \n"; 
  print "(2a)Used to locate the JSD file to use to create data series in DRMS\n";
  print "(2)Environment variable HK_NEW_JSD_FILE \n"; 
  print "(2a)Used to locate the list of new JSD series files recently create by jsoc_make_jsd_file.pl script.\n";
  print "(2b)Create data series in file if the apid is passed in list of decimal 4 digit apids or if ALL is passed.\n";

  exit;
}
