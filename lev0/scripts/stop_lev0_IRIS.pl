#!/usr/bin/perl
#
  print "I'm stopping all the ingest_lev0 IRIS processes...\n";
  $cmd = "touch /usr/local/logs/lev0/VC03_stop";
  `$cmd`;
  while(1) {
    $cfound = 0;
    #wait until they're all stopped
    @ps_prod = `ps -ef | grep ingest_lev0_iris`;
    while($_ = shift(@ps_prod)) {
      if(/vc VC03/) {
        $cfound = 1;
        sleep 1;
        last;
      }
    }
    if(!$cfound) { last; }
  }

