#!/usr/bin/perl
##############################################################################
# Name:        getdf.pl  - get day files, ingest into dayfile data series    #
#                          then validates files in data series and if there  #
#                          then remove day files from file system.           #
# Description: Get dayfiles from drop directory for sources hsb.             #
#              Call ingest_dayfile.pl and validate loaded dayfiles. If so    #
#              remove dayfile from file system.                              #
# Execution:   getdf.pl hsb                                                  #
#              getdf.pl rtmon (depricated, use dsdf.pl rtmon)                #
# Limitation:  The get dayfiles process works for dayfiles for the hsb       #
#              dayfiles only. This script is not used currently to move      #
#              dayfile to hk_rtmon_dayfile directory and                     #
#              call ingest_dayfiles.pl for LMSAL dayfiles. This script does  #
#              not check if dayfile is in data series before writing to data #
#              series, therefore can overwrite.  Script not used for moc prod#
#              server dayfiles(see dsdf.pl script). The help flag is not     #
#              implemented yet.                                              #
##############################################################################
# set Environment Variables
# debug flag
$dflg=$ENV{'DF_GETDF_DEBUG'}=0;

#common setting for all environments
$ENV{'SUMSERVER'}="k1";
#$hm=$ENV{'HOME'};
$hm="/home/jsoc/cvs/Development";
$ENV{'MAILTO'}="";
#$ENV{'DF_DRMS_EXECUTABLES'}="$hm/JSOC/bin/linux_x86_64";
$mach=$ENV{'MACHINE'};
$exec_dir=$ENV{'DF_EXEC_PATH'}="$hm/JSOC/bin/$mach";
$script_dir="$hm/JSOC/proj/lev0/scripts/hk";
$log_dir="/home/jsocprod/hk/logs";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$ENV{'DF_DRMS_EXECUTABLES'}";

#set common email arguments
$from_email="\"JSOC ANS\" \<jsoc_ans\@sun.Stanford.EDU\>";
$to_email="jsoc_ans\@sun.stanford.edu";

# check arguments
&check_agruments();

# set drop off directory for location of HSB or RTMON dayfiles,etc
# set log file based on source of dayfiles
# set subject mail arguments based on either HSB or RTMON args
if ($src eq "hsb")
{
  $doff_dir=$ENV{'DF_DROPOFF_HSB_FILES'}="/tmp28/jsocprod/lev0/hk_hsb_dayfile";
  $logfile="$log_dir/log-df-hsb";
  $subject_email="JSOC:WARNING:Ingesting HSB dayfiles: status:no files loaded today";
  $subject_email_old_df="JSOC:WARNING:Ingesting HSB dayfiles: status:Found Old HSB Dayfile Not Processed->Action Required";
  $subject_email_gt_currentdate="JSOC:WARNING:Ingesting HSB dayfiles: status:Found HSB Dayfile(s) Greater Than Current Date->Action Required";
  $subject_email_olddf_gtcd="JSOC:WARNING:Ingesting HSB dayfiles: status:Found Old HSB Dayfile(s) Not Processed and HSB Dayfile(s) Greater Than Current Date->Action Required";
}
elsif ($src eq "moc")
{
  $logfile="$log_dir/log-df-moc";
}
elsif ($src eq "egsefm")
{
  $logfile="$log_dir/log-df-egsefm";
}
elsif ($src eq "rtmon")
{
  $doff_dir=$ENV{'DF_DROPOFF_HSB_FILES'}="/tmp28/jsocprod/lev0/hk_rtmon_dayfile";
  $logfile="$log_dir/log-df-rtmon";
  $subject_email="JSOC:WARNING:Ingesting RTMON dayfiles: status:no files loaded today";
}
 
# set up where to put backup logs written monthly 
$logs_dir="$log_dir/old";

# open log file and append
open(LF,">>$logfile") || die "getdf.pl:1:Can't Open $logfile: $!\n";
print LF `date`;
print LF "--->Starting script getdf.pl\n";
print LF "--->Processing day files at directory:$doff_dir\n";

#ingest all dayfiles files
# for production
if  ($src eq "hsb")
{
  # get today's dayfile and avoid picking up next days dayfile
  #($today_date)=get_today_date(); #since process at 9PM Nightly use PDT
  ($today_date)=20170414;
  $enddate=$startdate=$today_date;

  #log status to logfile
  print LF "--->Date range processing are <$startdate> to <$enddate>\n";
  print LF "--->Start processing day files to data series using ingest_dayfile.pl script\n";
  close LF;

  #call ingest_dayfile.pl script 
  $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl  apidlist=$script_dir/df_apid_list_day_file_hsb start=$startdate end=$enddate dsnlist=$script_dir/df_apid_ds_list_for_hsb src=hsb merged=0`;

  #reopen log
  open(LF,">>$logfile") || die "getdf.pl:2:Can't Open $logfile: $!\n";
  print LF "--->Completed processing day files to data series using ingest_dayfile.pl script\n";

}
elsif ($src eq "rtmon")
{
  #log status to logfile
  print LF "--->Date range is to do all existing \n";
  print LF "--->Start processing day files to data series using ingest_dayfile.pl script\n";
  close LF;

  #call ingest_dayfile.pl script 
  $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl  apidlist=$script_dir/df_apid_list_day_file_rtmon dsnlist=$script_dir/df_apid_ds_list_for_rtmon src=rtmon`;

  #reopen log
  open(LF,">>$logfile") || die "getdf.pl:3:Can't Open $logfile: $!\n";
  print LF "--->Completed processing day files to data series using ingest_dayfile.pl script\n";
}
elsif ($src eq "moc")
{
  print "--->WARNING:Not used for source equal to <$src>. Note using movedf.pl for moc dayfile processing\n";
  exit();
}
elsif ($src eq "egsefm")
{
  print "--->WARNING:Not used for source equal to <$src>. Note need to update this script to handle egsefm dayfiles\n";
  exit();
}

#Check if there and then delete all dayfiles that where ingested in dayfile data series 
open(DELFILE, "$log_dir/DF_DELETE_FILE_LIST") || die "getdf.pl:4:Can't Open $log_dir/DF_DELETE_FILE_LIST file: $!\n";
@all_del_file_lines="";
$hkt_filecount =0;
while (<DELFILE>)
{
  $_ =~ s/\n//g;
  push(@all_del_file_lines, $_) ;
  if ($dflg) {print LF "validated file saved to drms.deleting file from file system directory:$_\n";}
  $log=`rm $doff_dir/$_`;
  if ($dflg) {print LF "$log\n";}
  $hkt_filecount++;
}#while
if($hkt_filecount > 0)
{
  print LF "--->Deleted dayfiles that where successfully processed into data series. Number loaded to drms and deleted from filesystem:$hkt_filecount\n";
}
else
{
  print LF "--->Skipping deleting dayfiles because no files ingested to hk_dayfile series\n";
  sendEmail("$to_email", "$from_email", "$subject_email", "Warning Message:\n-->Received count of hkt day files of <$hkt_filecount> files from directory $doff_dir\n-->When executing <$hm/JSOC/proj/lev0/scripts/hk/getdf.pl hsb > from cron job.\n");
}

close DELFILE;

#set MF file to blank work was completed
open(DELFILE, ">$log_dir/DF_DELETE_FILE_LIST") || die "getdf.pl:5:Can't Open $log_dir/DF_DELETE_FILE_LIST file: $!\n";
close DELFILE;

# check if old hsb dayfile were left in directory and send email warning if see old dayfiles or future dayfiles
if  ($src eq "hsb")
{
  &check_for_old_dayfiles($doff_dir);
}

#close logfile
print LF "--->exiting script getdf.pl\n";
print LF `date`;
close LF;


# move log to logs directory every month
&check_log();

# check if there are any old dayfiles in directory



##############################################################################
# check_arguments()                                                          #
##############################################################################
sub check_agruments()
{
  #source index in hk_dayfile data series
  $src=$ARGV[0]; 
  #check arguments
  if ($#ARGV != 0 )
  {
    print "Usage: perl getdf.pl  <dayfile source>\nwhere source is either:hsb,moc,egsefm or rtmon.\n";
    exit;
  }
  elsif ("-h" eq substr($ARGV[0],0,2) )
  {
     print "Usage: perl getdf.pl  <dayfile source>\nwhere source is either:hsb,moc,egsefm,or rtmon.\n";
     exit;
  }
  elsif("moc" eq substr($ARGV[0],0,3) or "hsb" eq substr($ARGV[0],0,3) or "egsefm" eq substr($ARGV[0],0,6)  or "rtmon" eq substr($ARGV[0],0,5))
  {
    #print "okay\n";
  }
  else
  {
     print "ERROR: Entered incorrect source name for data series. Use moc,hsb,\n";
     print "Usage: perl getdf.pl  <dayfile source>\nwhere source is either:hsb,moc,egsefm, or rtmon.\n";
     exit;
  }
}



##############################################################################
# get_today_date()                                                           #
##############################################################################
sub get_today_date
{
  ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  $year = 1900 + $yearOffset;
  return(sprintf("%-04.4d%-02.2d%-02.2d",$year,$month+1,$dayOfMonth));
}


##########################################
# check_log: used to create monthly logs #
##########################################
sub check_log()
{
  use Time::Local 'timelocal';
  use POSIX;

  # get current time
  ($sec,$min,$hour,$mday,$monoffset,$yearoffset,$wday,$yday,$isdst) = localtime();

  # set year using year offset
  $year= $yearoffset + 1900;

  # set month using month offset
  $mon=$monoffset + 1;

  #if date is the month of day is 1th then backup log
  #since will run once a day should get one log
  if ($mday == 1 )
  {
    #check if logs directory exists
    if ( -e $logs_dir)
    {
      my $d=`date`;
      #regular expression to add - were there are blanks
      $d =~ s/^| /-/g;
      $lm=`cp $log_dir/$logfile $logs_dir/$logfile-$d`;
      #set log file to blank - copy was completed
      open(LF, ">$log_dir/$logfile") || die "getdf.pl:6:Can't Open $log_dir/$logfile file: $!\n";
      close LF;
    }
    else
    {
      print "WARNING:movedf.pl:missing logs directory:<$logs_dir>. Create one at <$log_dir>\n";
    }
  }
}

#####################
# sendEmail         #
#####################
# Simple Email Function
# ($to, $from, $subject, $message)
sub sendEmail
{
my ($to, $from, $subject, $message) = @_;
my $sendmail = '/usr/lib/sendmail';
open(MAIL, "|$sendmail -oi -t");
print MAIL "From: $from\n";
print MAIL "To: $to\n";
print MAIL "Subject: $subject\n\n";
print MAIL "$message\n";
close(MAIL);
}



##########################
# check_for_old_dayfiles #
##########################
sub check_for_old_dayfiles($)
{
  #argument passed
  my $drop_off_dir=$_[0];

  #local variables
  my $count_old_dayfiles_found=0;
  my $count_gt_date_dayfiles_found=0;
  my $dfile;
  my $today_date;
  my @old_dafile_list=();
  my @old_dayfiles_found=();
  my @gt_date_dayfiles_found=();#greater than date dayfiles found
  my @pl_old_dayfiles_found=();
  my @pl_gt_date_dayfiles_found=();
  my $item;
 
  # get today's date
  ($today_date)=get_today_date_utc(); 

  # open directory with dayfiles
  opendir(DIR_DOFF, $drop_off_dir) || die "getdf.pl:7:Can't open directory:$drop_off_dir: $!\n"; #open subdirectory

  # read in all map files and put in list 
  @old_dayfile_list = readdir(DIR_DOFF); #get a list of directory contents

  # loop thru files looking for old dayfiles and files greater than current day's date.
  foreach $dfile  (@old_dayfile_list)
  {
    if($dfile eq "." or $dfile eq ".." or substr($dfile,0,4) ne "hsb_")
    {
       next;
    }
    #strip number out of filename
    $strvalue = sprintf("%s%s%s", substr($dfile,9,4),substr($dfile,14,2),substr($dfile,17,2));
    $intvalue= int $strvalue;
    if( $intvalue < $today_date)
    {
         push( @old_dayfiles_found, $dfile);
    }
    elsif( $intvalue > $today_date)
    {
         push( @gt_date_dayfiles_found, $dfile);
    }
    #else -got current dayfiles-these ok.
  }# end foreach loop

  # close directory
  close DIR_DOFF;

  # now check if need to send emails warning on extra dayfiles
  # get count of items in lists
  $count_old_dayfiles_found=@old_dayfiles_found;
  $count_gt_date_dayfiles_found=@gt_date_dayfiles_found;

  #create printable lists for email message
  foreach $item (@gt_date_dayfiles_found)
  {
    push(@pl_gt_date_dayfiles_found,"\n");
    push(@pl_gt_date_dayfiles_found, $item);
  }
  foreach $item (@old_dayfiles_found)
  {
    push(@pl_old_dayfiles_found,"\n");
    push(@pl_old_dayfiles_found, $item);
  }

  if ( $count_old_dayfiles_found > 0 && $count_gt_date_dayfiles_found == 0)
  {
    print LF "--->there are old dayfiles in directory\n";
    # get file list ready to put in email
    sendEmail("$to_email", "$from_email", "$subject_email_old_df", "Warning Message:\n-->When checked the HSB Directory found  one issue requiring action.\n\n-->ISSUE Found:\n(1)Old Dayfile(s) were found and need to be loaded in dayfile series.\n\n-->ACTIONS To Do:\n(1)Run ingest_dayfile.pl with src=hsb_r to load these old dayfiles. Then confirm loaded in hk_dayfile series and then delete each dayfile(s).\n\n-->List of <$count_old_dayfiles_found> Old Dayfile(s) are:@pl_old_dayfiles_found\n\n");
  }
  elsif ( $count_old_dayfiles_found == 0 && $count_gt_date_dayfiles_found > 0)
  {
    print LF "--->there are old dayfiles in directory\n";
    sendEmail("$to_email", "$from_email", "$subject_email_gt_currentdate", "Warning Message:\n-->When checked the HSB Directory found one issue requiring action.\n\n-->ISSUE Found:\n(1)Found dayfile(s) greater than current date and these probably should be removed.\n\n-->ACTIONS To Do:\n(1)The greater than current date dayfile(s) should probably be removed using rm.\n\n-->List of <$count_gt_date_dayfiles_found>> Greater than Current Date Dayfiles are:@pl_gt_date_dayfiles_found\n");
  }
  elsif ( $count_old_dayfiles_found > 0 && $count_gt_date_dayfiles_found > 0)
  {
    print LF "--->there are old dayfiles in directory\n";
    sendEmail("$to_email", "$from_email", "$subject_email_olddf_gtcd", "Warning Message:\n-->When checked the HSB Directory found two issues requiring action.\n\n-->ISSUES Found:\n(1)Old Dayfiles were found and need to be loaded in dayfile series\n(2)Found dayfile(s) greater than current date and these probably should be removed.\n\n-->ACTIONS To Do:\n(1)Run ingest_dayfile.pl with src=hsb_r to load these old dayfiles. Then confirm loaded in hk_dayfile series and then delete each dayfile.\n(2)The greater than current date dayfile(s) should probably be removed using rm.\n\n-->List of <$count_old_dayfiles_found> Old Dayfiles are:@pl_old_dayfiles_found\n\n-->List of <$count_gt_date_dayfiles_found> Greater than Current Date Dayfiles are:@pl_gt_date_dayfiles_found\n");

  }
  else
  {
    print LF "--->no old dayfiles in directory\n";
  }
}
##########################
# get_today_date_utc     #
##########################
sub get_today_date_utc
{
  ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
  $year = 1900 + $yearOffset;
  return(sprintf("%-04.4d%-02.2d%-02.2d",$year,$month+1,$dayOfMonth));
}
