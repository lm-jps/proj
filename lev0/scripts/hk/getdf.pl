#!/usr/bin/perl
##############################################################################
# Name:        getdf.pl  - get day files, ingest into dayfile data series    #
#                          then valid files in data series and if there the  #
#                          remove day files from file system.                #
# Description: Get dayfiles from drop directory for sources hsb and rtmon.   #
#              Call ingest_dayfile.pl and validate loaded dayfiles. If so    #
#              remove dayfile from file system.                              #
# Execution:   getdf.pl hsb                                                  #
#              getdf.pl rtmon                                                #
# Limitation:  The get dayfiles process works for dayfiles from              #
#              the hsb and rtmon dayfiles only. This script is not used      #
#              currently to move dayfile to hk_rtmon_dayfile directory and   #
#              call ingest_dayfiles.pl for LMSAL dayfiles. This script does  #
#              not check if dayfile is in data series before writing to data #
#              series, therefore can overwrite.  Script not used for moc prod#
#              server dayfiles(see movedf.pl script). The help flag is not   #
#              implemented yet.                                              #
##############################################################################
# set Environment Variables

# debug flag
$dflg=$ENV{'DF_GETDF_DEBUG'}=0;

#common setting for all environments
$ENV{'SUMSERVER'}="j1.Stanford.edu";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$ENV{'DF_DRMS_EXECUTABLES'}="$hm/cvs/JSOC/bin/linux_x86_64";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$ENV{'DF_DRMS_EXECUTABLES'}";

#set common email arguments
$from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
$to_email="jsoc_ops\@sun.stanford.edu";
#$from_email="carl\@sun.stanford.edu";
#$to_email="carl\@sun.stanford.edu";
$subject_email="JSOC:WARNING:Ingesting HSB dayfiles: status:no files loaded today";

# check arguments
&check_agruments();

# set drop off directory for location of HSB or RTMON dayfiles,etc
# set log file based on source of dayfiles
# set subject mail arguments based on either HSB or RTMON args
if ($src eq "hsb")
{
  $doff_dir=$ENV{'DF_DROPOFF_HSB_FILES'}="/tmp21/production/lev0/hk_hsb_dayfile";
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-hsb";
  $subject_email="JSOC:WARNING:Ingesting HSB dayfiles: status:no files loaded today";
}
elsif ($src eq "moc")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-moc";
}
elsif ($src eq "egsefm")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-egsefm";
}
elsif ($src eq "rtmon")
{
  $doff_dir=$ENV{'DF_DROPOFF_HSB_FILES'}="/tmp21/production/lev0/hk_rtmon_dayfile";
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-rtmon";
  $subject_email="JSOC:WARNING:Ingesting RTMON dayfiles: status:no files loaded today";
}

# set up where to put backup logs written monthly 
$logs_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk/logs";
#$logs_dir="$hm/cvs/myprod/JSOC/proj/lev0/scripts/hk/logs";

# open log file and append
open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
print LF `date`;
print LF "--->Starting script getdf.pl\n";
print LF "--->Processing day files at directory:$doff_dir\n";

#ingest all dayfiles files
# for production
if  ($src eq "hsb")
{
  # get today's dayfile and avoid picking up next days dayfile
  ($today_date)=get_today_date(); 
  $enddate=$startdate=$today_date;

  #log status to logfile
  print LF "--->Date range processing are <$startdate> to <$enddate>\n";
  print LF "--->Start processing day files to data series using ingest_dayfile.pl script\n";
  close LF;

  #call ingest_dayfile.pl script 
  $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl  apidlist=$script_dir/df_apid_list_day_file_hsb start=$startdate end=$enddate dsnlist=$script_dir/df_apid_ds_list_for_hsb src=hsb merged=0`;

  #reopen log
  open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
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
  open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
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
open(DELFILE, "$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
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
  sendEmail("$from_email", "$to_email", "$subject_email", "Warning Message:\n-->Received count of hkt day files of <$hkt_filecount> files from directory $doff_dir\n-->When executing </home/production/cvs/JSOC/proj/lev0/scripts/hk/getdf.pl hsb > from cron job.\n");
}

close DELFILE;

#set MF file to blank work was completed
open(DELFILE, ">$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
close DELFILE;
print LF "--->exiting script getdf.pl\n";
print LF `date`;
close LF;

# move log to logs directory every month
&check_log();

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
      $lm=`cp $script_dir/$logfile $logs_dir/$logfile-$d`;
      #set log file to blank - copy was completed
      open(LF, ">$script_dir/$logfile") || die "Can't Open $script_dir/$logfile file: $!\n";
      close LF;
    }
    else
    {
      print "WARNING:movedf.pl:missing logs directory:<$logs_dir>. Create one at <$script_dir>\n";
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
