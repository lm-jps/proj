#!/usr/bin/perl
#############################################################################
# Name:        cjds.pl - Create jsd data series                             #
# Description: Creates JSD if have a new JSD version for HK by APID JSD's   #
# Execution:   (1) To run :                                                 #
# Execution:   cjds.pl                                                      # 
# Limitation:  The script needs these environment variables set properly:   #
#              HK_NEW_JSD_TARGET, HK_JSD_DIRECTORY, HK_DDF_PROJECT_NAME,    #
#              HK_DDF_DATA_ID_NAME, HK_MHDS_APID_LIST.                      #
# NOTES:       Need to add help information.                                #
# Author:      Carl                                                         #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008#
#############################################################################
# main program  
# Set environment variables and variables
$hm=$ENV{'HOME'};
$ENV{'HK_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file";
$target_jsd_dir=$ENV{'HK_NEW_JSD_TARGET'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file";
$new_jsd_dir=$ENV{'HK_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";

#set up data series name parameters.
$pjn=$ENV{'HK_DDF_PROJECT_NAME'}="instrument"; 
$dtn=$ENV{'HK_DDF_DATA_ID_NAME'}="lev0testscript";

# set list of apids to do-if apid has new jsd file create data series in drms
$apids_to_do=$ENV{'HK_MHDS_APID_LIST'}="0017 0019 0021 0029 0445 0475 0529 0569";
#$apids_to_do = "ALL";

# Common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";
$ENV{'CVS_RSH'}="ssh";
#$ENV{'CVS_RSH'}="/usr/bin/rsh";

# check if do help
&show_help_info;

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
   print "-->Found <$account> account\n";
   print "ERROR: Found account name <$account>- Need to update cjds.pl script with SSH settings.\n";
   print LF `date`;
   die "ERROR: Found account name <$account>- Need to update cjds.pl script with SSH settings";
}
else
{
   print "ERROR: did not find this account name <$account>- Need to update cjds.pl script.\n";
   print `date`;
   die "ERROR: did not find this account name <$account>- Need to update cjsds.pl script.";
}

# log message to standard out
print "-->executing script cjds.pl to create JSD data series if needed\n";

# check list of new jsd files created in text file.
$njsds="$script_dir/new_jsd_files.txt";

# backup list of new jsd's created in case want to re-run
if ( -s $njsds )
{
  my $d=`date`;
  #regular expression to add - were there are blanks
  $d =~ s/^| /-/g;
  $lm=`cp $njsds $script_dir/new_jsd_files.txt-$d`;
}

#open file containing list of new jsd files created
open(NJS, "$njsds") || die "Can't Open $njsds: $!\n";

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
  print "-->new jsd file has apid is $s_line[2]\n";
  print "-->apid to do list is $apids_to_do\n";

  # for each new jsd file check if in list of apids to do
  foreach $atd  (@s_list)
  {
    # if could not find  apid in jsd file name check next
    if ((( index  $s_line[2],  $atd ) == -1) and (index  $atd,  "ALL") == -1) 
    {
      print "-->did not find <$s_line[2]> when checked $atd\n";
    }
    else
    {
      # log information to standard out
      print "-->Found <$s_line[2]> in $apids_to_do\n";
      print "-->project name is $pjn,data type name is $dtn,apid value is $s_line[2], jsvn value is $s_line[3]\n";

      # create dynamically the instrument part of project name based on apid value
      # for prelaunch using hmi_ground or aia_ground
      if ( int $s_line[2] > 0 and int $s_line[2] < 34 )
      {
        # $pjn =~ s/^instrument/hmi/g;
         $pjn =~ s/^instrument/su_carl/g;
         print "-->instrument project name used is $pjn\n";
      }
      elsif ( int $s_line[2] > 33 and int $s_line[2] < 64 ) 
      {
         #$pjn =~ s/^instrument/aia/g;
         $pjn =~ s/^instrument/su_carl/g;
         print "-->instrument project name used is $pjn\n";
      }
      elsif ( int $s_line[2] > 399 and int $s_line[2] < 500 ) 
      {
         #$pjn =~ s/^instrument/hmi/g;
         $pjn =~ s/^instrument/su_carl/g;
         print "-->instrument project name used is $pjn\n";
      }
      elsif ( int $s_line[2] > 499 and int $s_line[2] < 600 ) 
      {
         #$pjn =~ s/^instrument/aia/g;
         $pjn =~ s/^instrument/su_carl/g;
         print "-->instrument project name used is $pjn\n";
      }
      else
      {
         print "-->not hmi or aia apid value $s_line[2]- skipping creating\n"; 
         next;
      }

      # create data series 
      $data_series_name=sprintf("%s.%s_%s_%s",$pjn,$ENV{'HK_DDF_DATA_ID_NAME'},$s_line[2],$s_line[3]);

      # log information
      print "-->new ds name and data series to create will be < $data_series_name >\n";
      print "-->used project and data name to create new jsd file name from <$f>\n";

      # replace series name with new name of data series and copy to JSOC jsd directory
      $lm=`chmod 777  $new_jsd_dir/$f`;
      print "1) new jsd dir is  $new_jsd_dir\nfile is $f\nlog is $lm\n";
      $lm=`sed \"s/$f/$data_series_name/\" $new_jsd_dir/$f >   $target_jsd_dir/$data_series_name`; 
      print "2) target jsd dir is $target_jsd_dir\n data series name is $data_series_name\nlog is $lm\n";
      $lm=`chmod 777 $target_jsd_dir/$data_series_name `;
      print "3) finish chmod of $target_jsd_dir and $data_series_name\nlog is  $lm\n";

      # create series using create_series in Carl's JSOC environment for now -move to dimp eventually
      #$lm= `ssh -l carl n00.stanford.edu \"/home1/carl/jsoc/bin/linux_ia32/create_series $target_jsd_dir/$data_series_name\" `;
      $lm= `ssh -l carl n00.stanford.edu "create_series $target_jsd_dir/$data_series_name" `;
      print "4) create series using target jsd dir $target_jsd_dir and series name $data_series_name\n log is $lm\n";

      last;
    }
  }
}
print "-->Finished running cjds.pl - create jsoc data series\n";
close NJS;
## clear file
open(NJS, ">$njsds") || die "Can't Open $njsds: $!\n";
close NJS;

#################
# show_help_info#
#################
sub show_help_info
{
  if ( $ARGV[0] eq "-h" )
  {
    print "Help Listing\n";
    print "(1)Ways to Execute Perl Script: \n";
    print "(1a)Create Data Series in DRMS based on item in file set by HK_MHDS_APID_LIST: cjds.pl\n";
    print "(1b)Get Help Information: cjds.pl -h\n";
    print "(2)Environment variable HK_JSD_DIRECTORY \n"; 
    print "(2a)Used to locate the JSD file to use to create data series in DRMS\n";
    print "(3) Environment variable HK_NEW_JSD_TARGET\n";
    print "(3a) Used to copy new JSD files to in JSOC environment.\n";
    print "(4) Environment variable HK_DDF_PROJECT_NAME\n";
    print "(4a) Used to set project name of new data series to save in DRMS(i.e.,hmi_ground, su_carl, etc)\n";
    print "(5) Environment variable HK_DDF_DATA_ID_NAME\n";
    print "(5a) Used to set data type name of new data series to save in DRMS(i.e.,lev0, lev0test,etc.)\n";
    print "(6) Environment variable HK_MHDS_APID_LIST\n";
    print "(6a) Used to set a list of APIDS to create data series in DRMS(i.e, 0029 0019 0017 0021)\n";
    exit;
  }
}
