#!/usr/bin/perl
#module_flatfield_daily_qsub_PZT.pl
#
#Run daily to get a good flatfield using the hmi.lev1_nrt.
#Usually run by cvs/JSOC/proj/lev0/scripts/module_flatfield_daily_cron.pl or
#called directly  by lev1_def_gui_called_PZT_FSN.
#This updates cosmic_rays/rotation flat for the day.
#This is done after the definitive hmi.lev1.
#
#See Richard Wachter mail 9/14/10 15:13 "flatfield module runs"
#
sub usage {
  print "Run daily to get a good flatfield using the hmi.lev1.\n";
  print "Usage: module_flatfield_daily_qsub_PZT.pl hmi.lev1 2010.05.06\n";
  print " args:  input ds name\n";
  print "        date for the FF\n";
  exit(0);
}

$IN1 = "hmi.lev1";  #use to be hmi.lev1_nrt
$IN2 = "hmi.lev1";  #use to be hmi.lev1c_nrt
$QDIR = "/surge40/jsocprod/qsub/flat"; #dir for qsub scripts
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
$ENV{'OMP_NUM_THREADS'} = 8;

#$BINDIR = "/home/production/cvs/JSOC/bin/linux_x86_64";

$date = &get_date;
print "Start module_flatfield_daily_qsub_PZT.pl on $date\n\n";
if($#ARGV != 1) { &usage; }
$inds = $ARGV[0];
$datum = $ARGV[1];
if($inds ne $IN1 && $inds ne $IN2) {
  print "Error: Currently the only accepted input ds are $IN1 or $IN2\n";
  exit;
}
#verify date is like 2010.05.06. (module_flatfield fails on non-numeric)
$pos = index($datum, '.');
if($pos != 4) { &usage; };
$pos = index($datum, '.', $pos+1);
if($pos != 7) { &usage; };
$x = $datum."X";
$pos = index($x, 'X');
if($pos != 10) { &usage; };

#Here are the 48 command that we must run:
@cmds = (
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10054 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10055 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10056 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10057 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10058 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10058 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10059 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10059 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10074 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10075 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10076 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10077 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10078 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10078 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10079 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10079 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10094 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10095 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10096 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10097 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10098 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10098 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10099 camera=1  datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10099 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10114 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10115 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10116 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10117 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10118 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10118 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10119 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10119 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10134 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10135 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10136 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10137 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10138 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10138 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10139 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10139 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10154 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10155 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10156 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10157 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10158 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10158 camera=2 datum='$datum'",
"module_flatfield input_series='$inds' cadence=135 cosmic_rays=1 flatfield=1 fid=10159 camera=1 datum='$datum'",
"module_flatfield input_series='$inds' cadence=45 cosmic_rays=1 flatfield=1 fid=10159 camera=2 datum='$datum'");


$date = &get_date;
print "Start of module_flatfield 48 command run: $date\n";
$j = 0;
  @jid = ();
  for($i=0; $i<48; $i++) {
    $file = "$QDIR/qsh.$PID.$i.csh";
    open(L0, ">$file") || die "Can't open: $file $!\n";
    print L0 "#!/bin/csh\n";
    print L0 "echo \"TMPDIR = \$TMPDIR\"\n";
    print L0 "setenv OMP_NUM_THREADS 8\n";
    $cmd = @cmds[$j++];
    print L0 "$cmd\n";
    close(L0);
    #$file = "/home/production/cvs/JSOC/proj/lev0/scripts/date.str"; #!!!TEMP
    #$qsubcmd = sprintf("qsub -o %s -e %s -q p.q %s", $QDIR, $QDIR, $file);
    #$qsubcmd = sprintf("qsub -o %s -e %s -q j8.q %s", $QDIR, $QDIR, $file);
    $qsubcmd = sprintf("qsub -o %s -e %s -q a8.q %s", $QDIR, $QDIR, $file);
    print "$qsubcmd\n";
    $x = `$qsubcmd`;
    print "$x";
    ($a,$b,$jid) = split(/\s+/, $x);
    print "jid = $jid\n\n";
    push(@jid, $jid);
  }
NOMORE:
  while(1) {
    sleep(60);
    @stat = `qstat -u production`;
    #print "@stat\n";
    shift(@stat); shift(@stat);	#header lines
    $done = 1;
    while($line = shift(@stat)) {
      ($jid) = split(/\s+/, $line);
      if(grep(/$jid/, @jid)) {
        $done = 0;
        print "Found jid=$jid\n";
        last;
      }
    }
    if(!$done) { next; }
    else { last; }
  }
$date = &get_date;
print "Stop of module_flatfield 48 command run: $date\n";
print "All module_flatfield qsub done\n";

#Now run the combine program
$lfile = "$QDIR/combine.$PID.log2";
$cmd = "/home/production/cvs/JSOC/bin/linux_x86_64/module_flatfield_combine camera=2 input_series='su_production.flatfield_fid' datum='$datum' 1>$lfile 2>&1";
print "$cmd\n";
if(system($cmd)) {
  print "ERROR: on $cmd\n";
  exit(1);
}

$lfile = "$QDIR/combine.$PID.log1";
$cmd = "/home/production/cvs/JSOC/bin/linux_x86_64/module_flatfield_combine camera=1 input_series='su_production.flatfield_fid' datum='$datum' 1>$lfile 2>&1";
print "$cmd\n";
if(system($cmd)) {
  print "ERROR: on $cmd\n";
  exit(1);
}

#New: 10Jan2011 Do cosmic_ray post processing. This populates hmi.cosmic_rays.
#See mail richard@sun.stanford.edu 01/04/11 10:40 "cosmic ray series"
#Here are the 24 command that we must run:
@postcmds = (
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=00",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=02",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=04",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=06",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=08",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=10",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=12",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=14",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=16",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=18",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=20",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=1 hour=22",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=00",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=02",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=04",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=06",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=08",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=10",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=12",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=14",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=16",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=18",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=20",
"cosmic_ray_post input_series='su_production.cosmic_rays' datum='$datum' camera=2 hour=22");

$date = &get_date;
print "Start of cosmic_ray_post 24 command run: $date\n";
$j = 0;
for($l=0; $l<12; $l++) {	#12 sets of 2 processes
  @jid = ();
  for($i=0; $i<2; $i++) {
    $file = "$QDIR/qshp.$PID.$l.$i.csh";
    open(L0, ">$file") || die "Can't open: $file $!\n";
    print L0 "#!/bin/csh\n";
    print L0 "echo \"TMPDIR = \$TMPDIR\"\n";
    $cmd = @postcmds[$j++];
    print L0 "$cmd\n";
    close(L0);
    #$qsubcmd = sprintf("qsub -o %s -e %s -q p.q %s", $QDIR, $QDIR, $file);
    #$qsubcmd = sprintf("qsub -o %s -e %s -q j8.q,o8.q %s", $QDIR, $QDIR, $file);
    $qsubcmd = sprintf("qsub -o %s -e %s -q j.q %s", $QDIR, $QDIR, $file);
    print "$qsubcmd\n";
    $x = `$qsubcmd`;
    print "$x";
    ($a,$b,$jid) = split(/\s+/, $x);
    print "jid = $jid\n\n";
    push(@jid, $jid);
  }
NOMOREP:
  while(1) {
    sleep(45);
    @stat = `qstat -u production`;
    #print "@stat\n";
    shift(@stat); shift(@stat);	#header lines
    $done = 1;
    while($line = shift(@stat)) {
      ($jid) = split(/\s+/, $line);
      if(grep(/$jid/, @jid)) {
        $done = 0;
        print "Found jid=$jid\n";
        last;
      }
    }
    if(!$done) { next; }
    else { last; }
  }
}
$date = &get_date;
print "Stop of cosmic_ray_post 24 command run: $date\n";
print "All cosmic_ray_post qsub done\n";

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

