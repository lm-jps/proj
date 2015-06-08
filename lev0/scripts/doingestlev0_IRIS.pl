#!/usr/bin/perl 
#/home/prodtest/cvs/JSOC/proj/lev0/scripts/doingestlev0_IRIS.pl
#Start up the two ingest_lev0 programs and periodically check that
#they're still running and restart them if not and they were not
#stopped by request (i.e. VC03_stop file exists).
#
$user = $ENV{'USER'};
if($user ne "jsocprod") {
  print "You must be user jsocprod to run\n";
  exit;
}
$host = `hostname -s`;
if(!($host =~ /cl2n01[78]/)) {
  print "This can only be run on the cl2n017/8 pipeline machine.\n";
  exit;
}
#@vcnames = ("VC01", "VC04", "VC02", "VC05"); #primary channels
@vcnames = ("VC03"); #primary channel for IRIS
$ldate = &labeldate();
@lognames = ("VC03_$ldate.log");

$rmcmd0 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[0]_stop /usr/local/logs/lev0/@vcnames[0]_stopX /usr/local/logs/lev0/@vcnames[0]_exit /usr/local/logs/lev0/@vcnames[0]_exitX";

#clean up the signal files
print "$rmcmd0\n";
if(system($rmcmd0)) {
  print "Failed: $rmcmd0\n";
  exit;
}

$cmd0 = "ingest_lev0_irisdc --loopconn -c vc=@vcnames[0] indir=/sds/soc2pipe/iris logfile=/usr/local/logs/lev0/@lognames[0] &";
$log1 = sprintf("/usr/local/logs/lev0/%sX", @lognames[0]);

$cmd1 = "ingest_lev0_irisdc --loopconn -c vc=@vcnames[0] indir=/sds/soc2pipe/iris/rexmit logfile=$log1 &";

print "$cmd0\n";
if(system($cmd0)) {
  print "Failed: $cmd0\n";
}
print "$cmd1\n";
if(system($cmd1)) {
  print "Failed: $cmd1\n";
}
#$cmd = "iris_lev0_ck.pl &"; 
#print "$cmd\n";
#if(system($cmd)) {
#  print "Failed: $cmd\n";
#}
sleep(2);

print "\nTo cleanly stop all ingest_lev0_irisdc and commit any open images, call:\n";
print "stop_lev0_IRIS.pl\n";  #!!!TBD

#$SIG{INT} = \&catchc;
#$SIG{INT} = 'IGNORE';

while(1) {
  #for($i=0; $i < 4320; $i++) {  #run this for 6hrs
  for($i=0; $i < 720; $i++) {  #run this for 3600 seconds like HMI/AIA - Hao
  sleep 5;
    &ckingest;		#see if any ingest_lev0 are still running
    if(!$foundiris && !$foundirisrexmit) {  #nothing running
      $ldate = &labeldate();
      print "$ldate\n";
      $stopfile = "/usr/local/logs/lev0/VC03_stop";
      #if find a stop file, assume stop_lev0_IRIS.pl was called and exit
      if(-e $stopfile) {
        print "\nAll ingest_lev0 for IRIS have been stopped. Exit
doingestlev0_IRIS.pl\n\n";
        exit(0);
      }
      print "\nAll ingest_lev0 for IRIS are stopped. May have crashed!!\n\n";
      exit(0);
      #last;  #!!TEMP this didn't seem to work here ???
      #$i = 720;
    }
  }
  &ckstop;          #see if stop_lev0_IRIS.pl is running then do nothing
  if($ifound) {
    exit(0);
  }
  $cmd = "touch /usr/local/logs/lev0/@vcnames[0]_stop";
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[0]_stopX";
  `$cmd`;
  while(1) {
    $found = 0;
    #wait until they're all stopped
    @ps_prod = `ps -ef | grep ingest_lev0_irisdc`;
    while($_ = shift(@ps_prod)) {
      if(/vc VC03/) {
        $found = 1;
        sleep 1;
        last;
      }
    }
    if(!$found) { last; }
  }
  sleep(5);                     #make sure previous commit to db is done
      `$rmcmd0`;
      $cmd = "ingest_lev0_irisdc --loopconn -c -r vc=@vcnames[0] indir=/sds/soc2pipe/iris logfile=/usr/local/logs/lev0/VC03_$ldate.log &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      print "Restart:\n$cmd\n";
      $cmd = "ingest_lev0_irisdc --loopconn -c -r vc=@vcnames[0] indir=/sds/soc2pipe/iris/rexmit logfile=/usr/local/logs/lev0/VC03_$ldate.logX &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      print "Restart:\n$cmd\n";
}

sub ckingest {
    $foundiris = 0;
    $foundirisrexmit = 0;
    @ps_prod = `ps -ef | grep ingest_lev0_irisdc`;
    while($_ = shift(@ps_prod)) {
      if(/iris logfile/) {
        $foundiris = 1;
      }
      elsif(/iris\/rexmit logfile/) {
        $foundirisrexmit = 1;
      }
    }
}

sub ckstop {
    $ifound = 0;
    #see if stop_lev0_IRIS.pl is running
    @ps_prod = `ps -ef | grep stop_lev0_IRIS`;
    while($_ = shift(@ps_prod)) {
      if(/stop_lev0_IRIS.pl/) {
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

