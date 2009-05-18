#!/usr/bin/perl
# monitors if dayfile for today is getting updated and sends email if file not getting updated.
# also checks if file is there, if no file there it exits and prints error message to standard out.
# Example execution: monitor_df_rtmon.pl apid=129

#get apid argument
if (substr($ARGV[0],0,5) eq "apid=" )
{
  #use apid following -a for apid
  $apid= substr($ARGV[0],5);
  if ($apid eq "")
  {
    print "USAGE: monitor_df_rtmon apid=129\n";
    exit;
  }
}
else
{
  print "USAGE: monitor_df_rtmon apid=129\n";
  exit;
}

#set common email arguments
$hm=$ENV{'HOME'};
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'MAILTO'}="";
$from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
##$to_email="jsoc_ops\@sun.stanford.edu";
$to_email="carl\@sun.stanford.edu,rock\@solarpost.stanford.edu";
$subject_email="JSOC:WARNING:INGESTING REALTIME DAYFILE FOR APID $apid:Dayfile is not getting new packets";
$subject_email_error="JSOC:ERROR:INGESTING REALTIME DAYFILE FOR APID $apid:Dayfile does not exist";
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
 
#check if got error
if ($ino eq "")
{
  print "-->ERROR:Got error when did stat command. No file probably there\n";
  print "-->ERROR Details: stat:dev:<$dev> ino;<$ino> nlink:<$nlink> uid:<$uid>\n";
  sendEmail("$to_email_error", "$from_email", "$subject_email_error","Error Message:\n-->Dayfile <$dir_file_mon> is not probably there.\n-->This monitor script is exiting. To restart monitor enter command:  $script_dir/monitor_df_rtmon.pl apid=129\n");
  exit;
}

#initial current and previous size
$curr_size=$size;
$prev_size=0;
$already_triggered=0;

while (1)
{
  # sleep 5 minutes 
  sleep 300;
  print "-->WAITING 5 MINUTES THEN WILL CHECK IF TODAYS DAYFILE SIZE CHANGED.\n";
  #get file size 
  ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size, $atime,$mtime,$ctime,$blksize,$blocks) = stat($dir_file_mon);
  #check if got error
  if ($? ne 0)
  {
    print "ERROR: got error when did stat command\n";
    exit;
  }
  $curr_size=$size;
  print "-->IF CURRENT-FILE-SIZE:$curr_size EQUAL PREV-FILE-SIZE:$prev_size  THEN TRIGGER EMAIL\n";
  if($prev_size == $curr_size && $already_triggered != 1)
  {
     sendEmail("$to_email", "$from_email", "$subject_email","Warning Message:\n-->Dayfile <$dir_file_mon> is not getting updated for 5 minutes.\n-->This monitor script is exiting. To restart monitor enter command:  $script_dir/monitor_df_rtmon.pl apid=129\n");
     print "-->SEND WARNING EMAIL AFTER WAITING 5 MINUTES!\n";
     $already_triggered=1;
     exit;
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
