#!/usr/local/bin/perl
#Find the qac file times in the /usr/local/logs/soc/soc_iris_VC03* 
#file name given else the latest on.
#The output will look like:
#
#2013.07.24_18:52:56
#*Found qac file:
#* /sds/soc2soc/iris/VC03_2013_205_11_52_42_000035c175f_00.qacx
#--
#2013.07.24_20:34:22
#*Found qac file:
#* /sds/soc2soc/iris/VC03_2013_205_15_05_56_00003683366_00.qac
#--
#
#
$logdir = "/usr/local/logs/soc";
if($#ARGV != 0) {	#no file given, use the latest on
  $cmd = "ls -t $logdir/soc_iris_VC03_prodtest*";
  @ls = `$cmd`;
  $_ = shift(@ls);
}
else {
  $_ = $ARGV[0];
  if(!/^\//) { print "Must give full path name\n"; exit; }
}
print "Given input file:\n$_\n";
$cmd = "grep -B1 -C1 \"Found qac file\" $_";
print "$cmd\n";
@ans = `$cmd`;
print "@ans\n";

