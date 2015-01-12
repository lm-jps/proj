#!/usr/bin/perl
#lev0a2deflev1_FULL_PZT_FSN.pl 
#
#Call with the day to process, normally n-2 (changed to 4) days from 
#current new UTC day:
#  lev0a2deflev1_FULL_PZT_FSN.pl hmi 2010.298_UTC
#
#This is normally called by the cron job:
# lev0a2deflev1_FULL_PZT_FSN_cron.pl
#
#Run daily to get from hmi.lev0a to hmi.lev1.
#These must occur in this sequence:
#  Hao verifies (TBD) that all the DDS tlm files are processed to hmi.lev0a
#  This is normally so if a file like /dds/soc2pipe/hmi/xday/Xmit_All.279
#  exists (this file is created by the datacapture system to indicate that 
#  no further tlm files are expected for day 279. (TBD put in year)) and
#  no ingest_lev0 or sum_svc interruptions occured to leave lev0a gaps.
#
#Then proceed like this, for example:
#On 27Oct2010 (2010.300_UTC) after 26Oct2010 17:00 PDT:
#
#@defdb = `lev1_definitive_db.pl hmi 2010.298_UTC`;
#grep @defdb for '**NG '. If found, give error that can't proceed
#and print @defdb.
#
#Else all the pieces are in place to do the def lev1.
#
#Make the def lev1 for day n-2
#On cl1n001 (where the lev1_def_gui runs):
#Run the build of the hmi.lev1 and fill in any hmi.lev1 gaps:
#lev1_def_gui_called_PZT_FSN -x hmi 2010.298_UTC
#
#Make the rot flat for this day using the newly build hmi.lev1
#module_flatfield_daily_cron_48_PZT_FSN.pl hmi 2010.298_UTC  
#(this is now called in lev1_def_gui_called_PZT_FSN)
#
#!!OLD below:
#  In Phil's workflow run the day like so (!!NOTE: has been nooped out)
#  /home/phil/workflow/maketicket.csh gate=hmi.LOS wantlow=2010.10.24_23:54_TAI wanthigh=2010.10.25_23:54_TAI action=5
#  Or action=4 to redo only what's missing in the day
#

sub usage {
  print "Make definitive lev1 (and then higher products) from lev0a data\n";
  print "Usage: lev0a2deflev1_FULL_PZT_FSN.pl hmi 2010.298_UTC\n";
  print " args:  instrument\n";
  print "        ord date to do\n";
  close(LOG);
  exit(0);
}

$QDIR = "/surge40/jsocprod/qsub/flat"; #dir for qsub scripts
$LOGDIR = "/usr/local/logs/lev1_def";
$PID = getppid;
$host = `hostname -s`;
chomp($host);
print "host = $host\n";
#cl1n001 is a linux_x86_64

if ($host !~ 'cl1n001|solar3') {
  print "Error: This must be run on cl1n001 or solar3\n";
  exit;
}

if ($host =~ 'cl1n001') {
  $ENV{'JSOC_MACHINE'} = "linux_x86_64";
  $JSOC_MACHINE = linux_x86_64;
} else {
  $ENV{'JSOC_MACHINE'} = "linux_avx";
  $JSOC_MACHINE = "linux_avx";
}

$ENV{'JSOC_DBUSER'}="production";
$ENV{'PATH'} = "/home/jsoc/cvs/Development/JSOC/bin/$JSOC_MACHINE:/home/jsoc/cvs/Development/JSOC/scripts:/bin:/usr/bin:/SGE/bin/lx24-amd64:";

$ENV{'SGE_ROOT'} = "/SGE";
$sgeroot = $ENV{'SGE_ROOT'};
#print "SGE_ROOT = $sgeroot\n";
$path = $ENV{'PATH'};
print "path = $path\n";  #!!!TEMP
$mach = $ENV{'JSOC_MACHINE'};
#print "JSOC_MACHINE = $mach\n";  #!!!TEMP

$date = &get_date;
$logfile = "$LOGDIR/lev0a2deflev1_FULL_PZT_FSN."."$date";
open(LOG, ">$logfile")  || die "Can't open $logfile: $!\n";
print "Start lev0a2deflev1_FULL_PZT_FSN.pl on $date\n";
print LOG "Start lev0a2deflev1_FULL_PZT_FSN.pl on $date\n";
if($#ARGV != 1) { &usage; }
$instru = $ARGV[0];
$ord_date = $ARGV[1];
if($instru ne "hmi" && $instru ne "aia") {
  print "Error: The only accepted instru are 'hmi' or 'aia'\n";
  print LOG "Error: The only accepted instru are 'hmi' or 'aia'\n";
  close(LOG);
  exit;
}
$pos1 = index($ord_date, '.');
$pos2 = index($ord_date, '_');
if($pos1==-1 || $pos2==-1) {
  &usage;
}
$yr1 = substr($ord_date, 0, $pos1);
$yrday1 = substr($ord_date, $pos1+1, ($pos2-$pos1)-1);
$zone1 = substr($ord_date, $pos2+1);
if($zone1 ne "UTC") { &usage; }
print "Call is: lev0a2deflev1_FULL_PZT_FSN.pl $instru $ord_date\n\n";
print LOG "Call is: lev0a2deflev1_FULL_PZT_FSN.pl $instru $ord_date\n\n";
#NEW 04Jan2010. Now put call to lev1_definitive_db.pl in
#lev1_def_gui_called_PZT_FSN
#$cmd = "lev1_definitive_db.pl -o $instru $ord_date"; #!!TEMP use -o
#print "$cmd\n";
#print LOG "$cmd\n";
#@defdb = `$cmd`;
#if(grep(/\*\*NG /, @defdb)) {
#  print "Can't proceed. Preconditions not meet:\n@defdb\n";
#  print LOG "Can't proceed. Preconditions not meet:\n@defdb\n";
#  close(LOG);
#  exit;
#}

#All the preconditions are meet.
#Now queue up the def lev1 processing for our given date
#Also makes the higher order products by calling Phil's workflow.
$ord_dateEX = "$ord_date"."::FEX";
#$cmd = "touch /usr/local/logs/lev1_gui/$ord_dateEX";
#do the command immediately !!!NEW
 $statlog = "$LOGDIR/lev1gui_call_$date.log";
 $cmd = "/home/jsoc/cvs/Development/JSOC/base/sums/scripts/lev1_def_gui_called_PZT_FSN -x $instru $ord_date 1> $statlog 2>&1";
print "$cmd\n";
print LOG "$cmd\n";
`$cmd`;
#Now make FF for the day. 18Nov2010 Now do this in lev1_def_gui_called_PZT_FSN
#$cmd = "module_flatfield_daily_cron_48_PZT_FSN.pl $instru $ord_date";
#print "$cmd\n";
#print LOG "$cmd\n";
#@ffrun = `$cmd`;
print "**Complete lev0a2deflev1_FULL_PZT_FSN.pl for $instru $ord_date\n";
print LOG "**Complete lev0a2deflev1_FULL_PZT_FSN.pl for $instru $ord_date\n";
close(LOG);

sub get_date {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2);
  #($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $year4 = sprintf("%04d", $year+1900);
  $date = $year4.".".$mon2.".".$mday2._.$hour2.":".$min2.":".$sec2;
  return($date);
}

