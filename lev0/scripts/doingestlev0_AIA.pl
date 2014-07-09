#!/usr/bin/perl 
#/home/production/cvs/JSOC/proj/lev0/apps/doingestlev0_AIA.pl
#Start up the four ingest_lev0 programs and periodically cause them
#to exit and then restart them.
#ingest_lev0 vc=VC01 indir=/dds/soc2pipe/aia
#ingest_lev0 vc=VC04 indir=/dds/soc2pipe/aia
#ingest_lev0 vc=VC01 indir=/dds/soc2pipe/aia/rexmit
#ingest_lev0 vc=VC04 indir=/dds/soc2pipe/aia/rexmit

#
#$user = "388";
#if($user ne "388") {
#  print "You must be user production to run\n";
#  exit;
#}
$host = `hostname`;
if(!($host =~ /cl2n01[78]/)) {
  print "This can only be run on the cl2n0017/8 pipeline machine.\n";
  exit;
}
@vcnames = ("VC01", "VC04", "VC02", "VC05"); #primary channels
$ldate = &labeldate();
@lognames = ("VC01_$ldate.log", "VC04_$ldate.log", "VC02_$ldate.log", "VC05_$ldate.log");

$rmcmd0 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[0]_stop /usr/local/logs/lev0/@vcnames[0]_stopX /usr/local/logs/lev0/@vcnames[0]_exit /usr/local/logs/lev0/@vcnames[0]_exitX";
$rmcmd1 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[1]_stop /usr/local/logs/lev0/@vcnames[1]_stopX /usr/local/logs/lev0/@vcnames[1]_exit /usr/local/logs/lev0/@vcnames[1]_exitX";
#$rmcmd2 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[2]_stop /usr/local/logs/lev0/@vcnames[2]_exit";
#$rmcmd3 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[3]_stop /usr/local/logs/lev0/@vcnames[3]_exit";

#clean up the signal files
print "$rmcmd0\n";
if(system($rmcmd0)) {
  print "Failed: $rmcmd0\n";
  exit;
}
print "$rmcmd1\n";
if(system($rmcmd1)) {
  print "Failed: $rmcmd1\n";
  exit;
}
#print "$rmcmd2\n";
#if(system($rmcmd2)) {
#  print "Failed: $rmcmd2\n";
#  exit;
#}
#print "$rmcmd3\n";
#if(system($rmcmd3)) {
#  print "Failed: $rmcmd3\n";
#  exit;
#}

#Now do the initial start of the 4 ingest_lev0
$cmd0 = "ingest_lev0 --loopconn  vc=@vcnames[0] indir=/dds/soc2pipe/aia logfile=/usr/local/logs/lev0/@lognames[0] &";
$cmd1 = "ingest_lev0 --loopconn  vc=@vcnames[1] indir=/dds/soc2pipe/aia logfile=/usr/local/logs/lev0/@lognames[1] &";
$log4 = sprintf("/usr/local/logs/lev0/%sX", @lognames[0]);
$log5 = sprintf("/usr/local/logs/lev0/%sX", @lognames[1]);
$cmd4 = "ingest_lev0 --loopconn  vc=@vcnames[0] indir=/dds/soc2pipe/aia/rexmit logfile=$log4 &";
$cmd5 = "ingest_lev0 --loopconn  vc=@vcnames[1] indir=/dds/soc2pipe/aia/rexmit logfile=$log5 &";

#$cmd2 = "ingest_lev0 --loopconn  vc=@vcnames[2] indir=/dds/soc2pipe/hmi logfile=/usr/local/logs/lev0/@lognames[2] &";
#$cmd3 = "ingest_lev0 --loopconn  vc=@vcnames[3] indir=/dds/soc2pipe/hmi logfile=/usr/local/logs/lev0/@lognames[3] &";
print "$cmd0\n";
if(system($cmd0)) {
  print "Failed: $cmd0\n";
}
print "$cmd1\n";
if(system($cmd1)) {
  print "Failed: $cmd1\n";
}
print "$cmd4\n";
if(system($cmd4)) {
  print "Failed: $cmd4\n";
}
print "$cmd5\n";
if(system($cmd5)) {
  print "Failed: $cmd5\n";
}
#print "$cmd2\n";
#if(system($cmd2)) {
#  print "Failed: $cmd2\n";
#}
#print "$cmd3\n";
#if(system($cmd3)) {
#  print "Failed: $cmd3\n";
#}

print "\nTo cleanly stop all ingest_lev0 and commit any open images, call:\n";
#print "stop_lev0.pl\n";
print "stop_lev0_AIA.pl\n";

#$SIG{INT} = \&catchc;
#$SIG{INT} = 'IGNORE';

while(1) {
  for($i=0; $i < 720; $i++) {  #run this for 3600 seconds
    sleep 5;
    &ckingest;		#see if any ingest_lev0 are still running
    if(!$ifound) {
      $ldate = &labeldate();
      print "$ldate\n";
      $stopfile = "/usr/local/logs/lev0/VC01_stop";
      #if find a stop file, assume stop_lev0_AIA.pl was called and exit
      if(-e $stopfile) {
        print "\nAll ingest_lev0 for AIA have been stopped. Exit doingestlev0_AIA.pl\n\n";
        exit(0);
      }
      print "\nAll ingest_lev0 for AIA are stopped. May have crashed!!\n\n";
      exit(0);
      #last;  #!!TEMP this didn't seem to work here ???
      #$i = 720;
    }
  }
  &ckstop;	#see if stop_lev0_AIA.pl is running
  if($ifound) {
    exit(0);
  }
  $cmd = "touch /usr/local/logs/lev0/@vcnames[0]_stop"; #tell ingest to stop
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[0]_stopX"; #tell ingest to stop
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[1]_stop";
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[1]_stopX";
  `$cmd`;
  #$cmd = "touch /usr/local/logs/lev0/@vcnames[2]_stop";
  #`$cmd`;
  #$cmd = "touch /usr/local/logs/lev0/@vcnames[3]_stop";
  #`$cmd`;
  while(1) {
    $found = 0;
    #wait until they're all stopped
    @ps_prod = `ps -ef | grep ingest_lev0`;
    while($_ = shift(@ps_prod)) {
      #if(/vc VC01/ || /vc VC02/ || /vc VC04/ || /vc VC05/) {
      if(/vc VC01/ || /vc VC04/) {
        $found = 1;
        sleep 1;
        last;
      }
    }
    if(!$found) { last; }
  }
  sleep(5);			#make sure previous commit to db is done

      `$rmcmd0`;
      $cmd = "ingest_lev0 --loopconn -r vc=@vcnames[0] indir=/dds/soc2pipe/aia logfile=/usr/local/logs/lev0/@lognames[0] &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      `$rmcmd1`;
      $cmd = "ingest_lev0 --loopconn -r vc=@vcnames[1] indir=/dds/soc2pipe/aia logfile=/usr/local/logs/lev0/@lognames[1] &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      $cmd = "ingest_lev0 --loopconn -r vc=@vcnames[0] indir=/dds/soc2pipe/aia/rexmit logfile=$log4 &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      $cmd = "ingest_lev0 --loopconn -r vc=@vcnames[1] indir=/dds/soc2pipe/aia/rexmit logfile=$log5 &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }

}

sub ckingest {
    $ifound = 0;
    #see if they're all stopped
    @ps_prod = `ps -ef | grep ingest_lev0`;
    while($_ = shift(@ps_prod)) {
      if(/vc VC01/ || /vc VC04/) {
        $ifound = 1;
        last;
      }
    }
}

sub ckstop {
    $ifound = 0;
    #see if stop_lev0_AIA.pl is running
    @ps_prod = `ps -ef | grep stop_lev0_AIA`;
    while($_ = shift(@ps_prod)) {
      if(/stop_lev0_AIA.pl/) {
        $ifound = 1;
        last;
      }
    }
}

#Return date in form for a label e.g. 1998.01.07_14:42:00
sub labeldate {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2);
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

