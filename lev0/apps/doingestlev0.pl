#!/usr/bin/perl 
#/home/production/cvs/JSOC/proj/lev0/apps/doingestlev0.pl
#Start up the four ingest_lev0 programs and periodically cause them
#to exit and then restart them.
#
$user = "388";
if($user ne "388") {
  print "You must be user production to run\n";
  exit;
}
$host = `hostname`;
if(!($host =~ /cl1n001/)) {
  print "This can only be run on the cl1n001 pipeline machine.\n";
print "!!!TEMP in exception mode!!!\n";
  #exit;
}
@vcnames = ("VC01", "VC04", "VC02", "VC05"); #primary channels
$ldate = &labeldate();
@lognames = ("VC01_$ldate.log", "VC04_$ldate.log", "VC02_$ldate.log", "VC05_$ldate.log");

$rmcmd0 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[0]_stop /usr/local/logs/lev0/@vcnames[0]_exit";
$rmcmd1 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[1]_stop /usr/local/logs/lev0/@vcnames[1]_exit";
$rmcmd2 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[2]_stop /usr/local/logs/lev0/@vcnames[2]_exit";
$rmcmd3 = "/bin/rm -f /usr/local/logs/lev0/@vcnames[3]_stop /usr/local/logs/lev0/@vcnames[3]_exit";

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
print "$rmcmd2\n";
if(system($rmcmd2)) {
  print "Failed: $rmcmd2\n";
  exit;
}
print "$rmcmd3\n";
if(system($rmcmd3)) {
  print "Failed: $rmcmd3\n";
  exit;
}

#Now do the initial start of the 4 ingest_lev0
$cmd0 = "ingest_lev0  vc=@vcnames[0] indir=/dds/socdc/aia logfile=/usr/local/logs/lev0/@lognames[0] &";
$cmd1 = "ingest_lev0  vc=@vcnames[1] indir=/dds/socdc/aia logfile=/usr/local/logs/lev0/@lognames[1] &";
$cmd2 = "ingest_lev0  vc=@vcnames[2] indir=/dds/socdc/hmi logfile=/usr/local/logs/lev0/@lognames[2] &";
$cmd3 = "ingest_lev0  vc=@vcnames[3] indir=/dds/socdc/hmi logfile=/usr/local/logs/lev0/@lognames[3] &";
print "$cmd0\n";
if(system($cmd0)) {
  print "Failed: $cmd0\n";
}
print "$cmd1\n";
if(system($cmd1)) {
  print "Failed: $cmd1\n";
}
print "$cmd2\n";
if(system($cmd2)) {
  print "Failed: $cmd2\n";
}
print "$cmd3\n";
if(system($cmd3)) {
  print "Failed: $cmd3\n";
}

while(1) {
  sleep 3600;		#run ingest_lev0 for this long
  #sleep 1800;		#run ingest_lev0 for this long
  $cmd = "touch /usr/local/logs/lev0/@vcnames[0]_stop"; #tell ingest to stop
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[1]_stop";
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[2]_stop";
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/@vcnames[3]_stop";
  `$cmd`;
  $f0 = $f1 = $f2 = $f3 = 0;
  while(1) {
    if(-e "/usr/local/logs/lev0/@vcnames[0]_exit") {	#wait for exit signal
      `$rmcmd0`;
      sleep(10);		#make sure previous commit to db is done
      $cmd = "ingest_lev0 -r vc=@vcnames[0] indir=/dds/socdc/aia logfile=/usr/local/logs/lev0/@lognames[0] &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      $f0 = 1;
    }
    if(-e "/usr/local/logs/lev0/@vcnames[1]_exit") {	#wait for exit signal
      `$rmcmd1`;
      sleep(10);		#make sure previous commit to db is done
      $cmd = "ingest_lev0 -r vc=@vcnames[1] indir=/dds/socdc/aia logfile=/usr/local/logs/lev0/@lognames[1] &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      $f1 = 1;
    }
    if(-e "/usr/local/logs/lev0/@vcnames[2]_exit") {	#wait for exit signal
      `$rmcmd2`;
      sleep(10);		#make sure previous commit to db is done
      $cmd = "ingest_lev0 -r vc=@vcnames[2] indir=/dds/socdc/hmi logfile=/usr/local/logs/lev0/@lognames[2] &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      $f2 = 1;
    }
    if(-e "/usr/local/logs/lev0/@vcnames[3]_exit") {	#wait for exit signal
      `$rmcmd3`;
      sleep(10);		#make sure previous commit to db is done
      $cmd = "ingest_lev0 -r vc=@vcnames[3] indir=/dds/socdc/hmi logfile=/usr/local/logs/lev0/@lognames[3] &";
      if(system($cmd)) {
        print "Failed: $cmd\n";
      }
      $f3 = 1;
    }
    if($f0 && $f1 && $f2 && $f3) { last; }
    sleep 1;
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

