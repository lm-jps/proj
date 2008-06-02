#!/usr/bin/perl
#/home/production/cvs/JSOC/proj/datacapture/scripts/soc_logfile_scan.pl
#
#Will read the given soc log file for an instrument and vc, e.g. 
#soc_aia_VC01_production_2008.05.30_09:34:34.log
#and collect some statistics.
#

$YEARSCAN = "2008";	#look for this date stamp in log file

sub usage {
  print "Get info on a datacapture soc log file, e.g.\n";
  print "soc_aia_VC01_production_2008.05.30_09:34:34.log\n";
  print "Usage: soc_logfile_scanpl in_file\n";
  exit(1);
}

if($#ARGV != 0) {
  &usage;
}
$infile = $ARGV[0];
open(ID, $infile) || die "Can't open $infile: $!\n";
print "Scan of: $infile\n";
while(<ID>) {
  if(/^$YEARSCAN/) {
    print "$_";
  }
  if(/^\*\*Processed/) {
    print "$_";
    while(<ID>) {
      print "$_";
      last;
    }
  }
}
close(ID);

