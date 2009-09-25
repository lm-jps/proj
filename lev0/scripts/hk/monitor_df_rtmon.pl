#!/usr/bin/perl
# NAME: MONITOR REAL TIME DATA FLOW TO DAYFILES
# AUTHOR: carl
# DESCRIPTION: Monitors if dayfile for today is getting updated and sends email if file not getting updated.
# Also checks if file is there, if no file there it exits and prints error message to standard out.
# If send one email notice already for no-file or no-new-data case, the script will not send another
# email until the file is there or the new data comes to dayfile and then the new data flow stops.
# So we reduce amount of emails being sent.
# RUN EXAMPLE: monitor_df_rtmon.pl apid=129
# LIMITATION: Must run on machine were dayfiles can be read(i.e.,n02) 
# CREATED: 9/25/2009

#get apid argument
if (substr($ARGV[0],0,5) eq "apid=" )
{
  #use apid following -a for apid
  $apid= substr($ARGV[0],5);
  if ($apid eq "")
  {
    print "USAGE: monitor_df_rtmon.pl  apid=129\n";
    exit;
  }
}
else
{
  print "USAGE: monitor_df_rtmon.pl apid=129\n";
  exit;
}

#set common email arguments
$hm=$ENV{'HOME'};
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'MAILTO'}="";
$from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
##$to_email="jsoc_ops\@sun.stanford.edu";
$to_email="carl\@sun.stanford.edu,rock\@solarpost.stanford.edu";
#$to_email="carl\@sun.stanford.edu";
$subject_email="JSOC:WARNING:INGESTING REALTIME DAYFILE FOR APID $apid:Dayfile is not getting new packets";
$subject_email_error="JSOC:WARNING:INGESTING REALTIME DAYFILE FOR APID $apid:Dayfile does not exist";
$to_email_error="carl\@sun.stanford.edu";

#get today's date
$td=get_todays_date();

#get file to monitor today
$file_to_monitor= sprintf("%8d.0x%04.4x", $td,$apid);

#get directory to find files
$subdir=sprintf("0x%04.4x", $apid);
$directory="/hmisdp-mon/log/packets/$subdir";

#get file and directory to monitor today
$dir_file_mon=sprintf("%s/%s", $directory, $file_to_monitor);

#get file information
print "-->FILE TO MONITOR IS::::<$dir_file_mon>\n";
($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size, $atime,$mtime,$ctime,$blksize,$blocks) = stat($dir_file_mon);
 
# set to 0 to send email notification of error once.
$already_triggered=0;

#check if got error
if ($ino eq "")
{
  print "-->WARNING:1:Got error when did stat command. No file probably there\n";
  print "-->WARNING Details: stat:dev:<$dev> ino;<$ino> nlink:<$nlink> uid:<$uid>\n";
  sendEmail("$to_email", "$from_email", "$subject_email_error","Error Message:\n-->Dayfile <$dir_file_mon> is not probably there.\n-->This monitor script will continue running and resend another email notice if the data file is there and data starts flowing and then stops again.\n-->To restart monitor enter command:  $script_dir/monitor_df_rtmon.pl apid=$apid\n-->To stop monitor run command(user=production):$script_dir/stop_monitor_df_rtmon.pl apid=$apid\n\n");
  $already_triggered=1;  
}

#initial current and previous size
$curr_size=$size;
$prev_size=0;

while (1)
{
  # sleep 5 minutes 
  print "-->WAITING 5 MINUTES THEN WILL CHECK IF TODAYS DAYFILE SIZE CHANGED.\n";
  sleep 300;
  #sleep 15;
  #get file size 
  ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size, $atime,$mtime,$ctime,$blksize,$blocks) = stat($dir_file_mon);
  #check if got error
  if ($ino eq "")
  {
    print "-->WARNING:2:Got error when did stat command. No file probably there\n";
    print "-->WARNING Details: stat:dev:<$dev> ino;<$ino> nlink:<$nlink> uid:<$uid>\n";
    if ($already_triggered != 1)
    {
      sendEmail("$to_email", "$from_email", "$subject_email_error","Error Message:\n-->Dayfile <$dir_file_mon> is not probably there.\n-->This monitor script will continue running and resend another email notice if the data file is there and data starts flowing and stops again.\n-->To restart monitor enter command:  $script_dir/monitor_df_rtmon.pl apid=$apid\n-->To stop monitor run command(user=production):$script_dir/stop_monitor_df_rtmon.pl apid=$apid\n");
      $already_triggered=1;
    }
  }
  $curr_size=$size;
  print "-->IF CURRENT-FILE-SIZE:$curr_size EQUAL PREV-FILE-SIZE:$prev_size AND HAVE NOT ALREADY SENT EMAIL THEN TRIGGER SEND OF EMAIL NOTICE\n";
  if($prev_size == $curr_size )
  {
    if ($already_triggered != 1)
    {
      sendEmail("$to_email", "$from_email", "$subject_email","Warning Message:\n-->Dayfile <$dir_file_mon> is not getting updated for 5 minutes.\n-->This monitor script will continue running and resend another email notice if the data starts flowing and stops again.\n-->To restart monitor enter command:  $script_dir/monitor_df_rtmon.pl apid=$apid\n-->To stop monitor run command(user=production):$script_dir/stop_monitor_df_rtmon.pl apid=$apid\n");
      print "-->SENDING WARNING EMAIL AFTER WAITING 5 MINUTES!\n";
      $already_triggered=1;
    }
  }
  else
  {
     $already_triggered=0;
  }
  $prev_size=$curr_size;

  #check for new day and new file
  $td=get_todays_date();
  $file_to_monitor= sprintf("%8d.0x%04.4x", $td,$apid);
  $subdir=sprintf("0x%04.4x", $apid);
  $directory="/hmisdp-mon/log/packets/$subdir";
  $dir_file_mon=sprintf("%s/%s", $directory, $file_to_monitor);
  $dd=`date`;
  print "-->FILE TO MONITOR IS::::<$dir_file_mon> at $dd\n";
}
##########################################################################
# subroutine get_todays_date()                                           #
##########################################################################
sub get_todays_date()
{
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  if($dflg == 2) {print  "DEBUG:MESSAGE:get_todays_date: year this $year month is $month day is $dayOfMonth\n";}
  #create todays date format and push on list of dates to do
  $new_date= sprintf("%4d%02.2d%02.2d", $year,$month,$dayOfMonth);
  if($dflg == 2) {print "DEBUG:MESSAGE:get_todays_date: new date is ::: $new_date :::\n";}
  return ( $new_date);
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
