#!/usr/bin/perl
#cvs/JSOC/proj/lev0/scripts/module_flatfield_daily_cron_48_PZT_FSN.pl
#Called after the start of a new UT day to make the hmi.cosmic_rays 
#and  rotation flatfield and cosmic_ray_post
#for the day from the hmi.lev1 data (that is typically made just before
#this call).
#Call:
#  module_flatfield_daily_cron_48_PZT_FSN.pl hmi 2010.193_UTC firsfsn lastfsn
#
#Eventually calls: module_flatfield_daily_qsub_48_PZT_FSN.pl hmi.lev1 firstfsn lastfsn 2010.09.30
#
#This is typically called at the end of
#lev1_def_gui_called_PZT_FSN after the hmi.lev1 is made

$XmitFlgDirHMI = "/dds/soc2pipe/hmi";
$XmitFlgDirAIA = "/dds/soc2pipe/aia"; 
$NoGoFileHMI = "/tmp28/jsocprod/lev0/data/NoDefLev1HMI";
$NoGoFileAIA = "/tmp28/jsocprod/lev0/data/NoDefLev1AIA";
$LOGDIR = "/usr/local/logs/lev1_def";
@allowhost = ("cl1n001", "cl1n002", "cl1n003"); #hosts w/dcs[0,1] mounts
@mstr = ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December");
@mdays = (31,28,31,30,31,30,31,31,30,31,30,31);
@mdaysl = (31,29,31,30,31,30,31,31,30,31,30,31);
@dayofyr_mo = (1,32,60,91,121,152,182,213,244,274,305,335);
@dayofyrl_mo = (1,32,61,92,122,153,183,214,245,275,306,336);
$specialtime= 0;

$host = `hostname -s`;
chomp($host);
print "host = $host\n";
if ($host =~ 'cl1n001') {
  $ENV{'JSOC_MACHINE'} = "linux_x86_64";
  $JSOC_MACHINE = linux_x86_64;
} else {
  $ENV{'JSOC_MACHINE'} = "linux_avx";
  $JSOC_MACHINE = "linux_avx";
}

$ENV{'PATH'} = "/home/jsoc/cvs/Development/JSOC/bin/$JSOC_MACHINE:/home/jsoc/cvs/Development/JSOC/scripts:/bin:/usr/bin:/SGE/bin/lx24-amd64:";

$ENV{'SGE_ROOT'} = "/SGE";
$sgeroot = $ENV{'SGE_ROOT'};
#print "SGE_ROOT = $sgeroot\n";
$path = $ENV{'PATH'};
#print "path = $path\n";  #!!!TEMP
$mach = $ENV{'JSOC_MACHINE'};
$labeld = &labeldate;		#set $orddayUTC for today
$logfile = "$LOGDIR/module_flatfield_daily_cron_48_PZT_FSN."."$labeld";
open(LOG, ">$logfile")  || die "Can't open $logfile: $!\n";
print "log = $logfile\n";

sub usage {
  print "Called after the start of a new UT day to make the hmi.cosmic_rays 
 and rotation flatfield
 for the day from the hmi.lev1 data.\n";
  print "Normally called at the end of lev1_def_gui_called_PZT_FSN.\n";
  print "Usage: module_flatfield_daily_cron_48_PZT_FSN.pl hmi 2010.193_UTC firstfsn lastfsn\n";
  print LOG "Usage: module_flatfield_daily_cron_48_PZT_FSN.pl hmi 2010.193_UTC firstfsn lastfsn\n";
  print "Performs various validity checks before it will run (see script).\n";
  print LOG "Performs various validity checks before it will run (see script).\n";
  close(LOG);
  exit(0);
}

if($#ARGV != 3) {
  &usage;
}
$instru = $ARGV[0];
if($instru ne "hmi" && $instru ne "aia") { &usage; }
if($instru eq "hmi") {
  $hmiaiaflg = 0;
  $indata = "hmi.lev1";
} else {
  $hmiaiaflg = 1;
  $indata = "aia.lev1";
}
  $ordday = $ARGV[1];
  $pos1 = index($ordday, '.');
  $pos2 = index($ordday, '_');
  if($pos1==-1 || $pos2==-1) {
    &usage;
  }
  $yr1 = substr($ordday, 0, $pos1);
  $yrday1 = substr($ordday, $pos1+1, ($pos2-$pos1)-1);
  $zone1 = substr($ordday, $pos2+1);
  if($zone1 ne "UTC") { &usage; }
  $dodayUTC = $ordday;	#the date input
  $firstfsn = $ARGV[2];
  $lastfsn = $ARGV[3];
print "Started on $orddayUTC to process day $dodayUTC\n";
print "First fsn = $firstfsn  Last fsn = $lastfsn\n";
if($firstfsn > $lastfsn) {
  print "Error: last fsn must be >= firstfsn\n";
  exit;
}

print LOG "Started on $orddayUTC to process day $dodayUTC\n";
#make time format for module_flatfield_daily_qsub_PZT_FSN.pl call
$msec = `time_convert time=$dodayUTC`;
$moddate = `time_convert s=$msec`;
$pos = index($moddate, '_');
$mdate = substr($moddate, 0, $pos);

$localhost = `hostname -s`;
chomp($localhost);
#if(!grep(/$localhost/, @allowhost)) {
#  print "Can only be run on host with dcs[0,1] mounts: @allowhost\n";
#  print LOG "Can only be run on host with dcs[0,1] mounts: @allowhost\n";
#  close(LOG);
#  exit(0);
#}
$pos1 = index($dodayUTC, '.');
$pos2 = index($dodayUTC, '_');
$yr1 = substr($dodayUTC, 0, $pos1);
$yrday1 = substr($dodayUTC, $pos1+1, ($pos2-$pos1)-1);

#Check if no go file exists.
if($hmiaiaflg) {
  $nogo = $NoGoFileAIA;
  $xmitfile1 = "$XmitFlgDirAIA/xday/Xmit_All.$dodayUTC";
  $xmitfile2 = "$XmitFlgDirAIA/xday/Xmit_All.$yrday1";
} else {
  $nogo = $NoGoFileHMI;
  $xmitfile1 = "$XmitFlgDirHMI/xday/Xmit_All.$dodayUTC";
  $xmitfile2 = "$XmitFlgDirHMI/xday/Xmit_All.$yrday1";
}
#if(-e $nogo) {
#  print "The no-go file exist: $nogo\n";
#  print LOG "The no-go file exist: $nogo\n";
#  print "Cannot proceed with making rotational flat\n";
#  print LOG "Cannot proceed with making rotational flat\n";
#  close(LOG);
#  exit(0);
#}
print "Checking for $xmitfile1\n";
print LOG "Checking for $xmitfile1\n";
print "or $xmitfile2\n";
print LOG "or $xmitfile2\n";
#noop this check out to process the old data
#if(!-e $xmitfile1 && !-e $xmitfile2)  {
#  print "Can't run yet. Rexmit may still be pending\n";
#  print LOG "Can't run yet. Rexmit may still be pending\n";
#  close(LOG);
#  exit(0);
#}

$cmd = "module_flatfield_daily_qsub_48_PZT_FSN.pl $indata $firstfsn $lastfsn $mdate";
print "Calling: $cmd\n"; 
print LOG "Calling: $cmd\n"; 
#`$cmd`;
system($cmd);
print "module_flatfield_daily_cron_48_PZT_FSN.pl complete\n";
print LOG "module_flatfield_daily_cron_48_PZT_FSN.pl complete\n";
close(LOG);

#Return date in form for a label e.g. 1998.01.07_14:42:00
#and set global $YEAR and $MONTH and $mon and $NUMDAYS.
sub labeldate {
  local($sec,$min,$hour,$mday,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2);
  if(!$specialtime) {
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
  } else {
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime($specialtime);
  }
  if($year % 4 == 0 && $year % 100 != 0 || $year % 400 == 0) { #leap year
    $leap = 1;
  }
  else {
    $leap = 0;
  }
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $dayofmo = $mday2;                    #set global value
  $year4 = sprintf("%04d", $year+1900);
  $date = $year4.".".$mon2.".".$mday2._.$hour2.":".$min2.":".$sec2."_UTC";
  $YEAR = $year4;
  $MONTH = @mstr[$mon];
  if($leap) { $NUMDAYS = @mdaysl[$mon]; }
  else { $NUMDAYS = @mdays[$mon]; }
  $todayUTCdayofyr = $yday+1;   #1-366
  $orddayUTC = "$year4.$todayUTCdayofyr"."_UTC";
  return($date);
}

