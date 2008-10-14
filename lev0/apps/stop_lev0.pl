#!/usr/bin/perl
#
  print "I'm stopping all the ingest_lev0 processes...\n";
  $cmd = "touch /usr/local/logs/lev0/VC01_stop"; #tell ingest to stop
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/VC02_stop";
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/VC04_stop";
  `$cmd`;
  $cmd = "touch /usr/local/logs/lev0/VC05_stop";
  `$cmd`;
  while(1) {
    $cfound = 0;
    #wait until they're all stopped
    @ps_prod = `ps -ef | grep ingest_lev0`;
    while($_ = shift(@ps_prod)) {
      if(/vc VC01/ || /vc VC02/ || /vc VC04/ || /vc VC05/) {
        $cfound = 1;
        sleep 1;
        last;
      }
    }
    if(!$cfound) { last; }
  }

