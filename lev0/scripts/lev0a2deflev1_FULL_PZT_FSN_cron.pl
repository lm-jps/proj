#!/usr/bin/perl
#lev0a2deflev1_FULL_PZT_FSN_cron.pl 
#
#This is a cron job run by cl1n001 after the start of a new UT day.
#For example, after 2010.10.26_17:00:00 PDT is UTC day 2010.300_UTC.
#It will call lev0a2deflev1_FULL_PZT_FSN.pl with the date equal n-4 from
#the new current date, e.g. after  2010.10.26_17:00:00 PDT, the call
#will be:
#  lev0a2deflev1_FULL_PZT_FSN.pl hmi 2010.296_UC
#
#See lev0a2deflev1_FULL_PZT_FSN.pl for the details of what's run.
#NOTE: the lev0a tag is really a misnomer for aia. It uses aia.lev0.

sub usage {
  print "Make definitive lev1 from lev0a data\n";
  print "Called as a cron job on cl1n001 to make the n-4 day of new data.\n";
  print "Usage: lev0a2deflev1_FULL_PZT_FSN_cron.pl hmi\n";
  print " args:  instrument\n";
  exit(0);
}

$PID = getppid;
$host = `hostname -s`;
chomp($host);
print "host = $host\n";
#cl1n001 is a linux_x86_64
if($host ne "cl1n001") {
  print "Error: This must be run on cl1n001\n";
  exit;
}
$ENV{'JSOC_MACHINE'} = "linux_x86_64";
$JSOC_MACHINE = "linux_x86_64";
$ENV{'PATH'} = "/home/production/cvs/JSOC/bin/$JSOC_MACHINE:/home/production/cvs/JSOC/scripts:/bin:/usr/bin:/SGE/bin/lx24-amd64:";

$ENV{'SGE_ROOT'} = "/SGE";
$sgeroot = $ENV{'SGE_ROOT'};
#print "SGE_ROOT = $sgeroot\n";
$path = $ENV{'PATH'};
#print "path = $path\n";  #!!!TEMP
$mach = $ENV{'JSOC_MACHINE'};
#print "JSOC_MACHINE = $mach\n";  #!!!TEMP

$date = &get_date;
print "Start lev0a2deflev1_FULL_PZT_FSN_cron.pl on $date\n";
if($#ARGV != 0) { &usage; }
$instru = $ARGV[0];
if($instru ne "hmi" && $instru ne "aia") {
  print "Error: The only accepted instru are 'hmi' or 'aia'\n";
  exit;
}
$sec = `time_convert time=$date`;
#$sec = $sec - 172800;	#n-2 days
$sec = $sec - 345600;	#n-4 days
$ord_date = `time_convert s=$sec o=ord`;
print "$ord_date\n";

#New 16Feb2011 First do a pzt flat
if($instru eq "hmi") {
  `/home/production/cvs/JSOC/proj/lev0/scripts/pzt_flat_cron.pl`;
}

#Now make the def lev1 for our given date
#Also make the higher order products by calling Phil's workflow.
$cmd = "lev0a2deflev1_FULL_PZT_FSN.pl $instru $ord_date";
print "$cmd\n";
#`$cmd`;
system($cmd);

print "**Complete lev0a2deflev1_FULL_PZT_FSN_cron.pl $instru for $ord_date\n";

#Return current _UTC time
sub get_date {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2);
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
  #($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $year4 = sprintf("%04d", $year+1900);
  $date = $year4.".".$mon2.".".$mday2._.$hour2.":".$min2.":".$sec2."_UTC";
  return($date);
}

