#!/usr/bin/perl
#/home/production/cvs/JSOC/proj/datacapture/scripts/soc_logfile_scan.pl
#
#Will read the given soc log file for an instrument and vc, e.g. 
#soc_aia_VC01_production_2008.05.30_09:34:34.log
#and collect some statistics.
#NOTE: this is not designed to work for all time. It is a quick hack.

$YEARSCAN = "2009.09";	#look for this date stamp in log file

sub usage {
  print "Get info on a datacapture soc log file, e.g.\n";
  print "soc_aia_VC01_production_2008.05.30_09:34:34.log\n";
  print "Usage: soc_logfile_scan.pl in_file\n";
  exit(1);
}

if($#ARGV != 0) {
  &usage;
}
$infile = $ARGV[0];
open(ID, $infile) || die "Can't open $infile: $!\n";
print "Scan of: $infile\n";
while(<ID>) {			#first get the start time
  if(/^$YEARSCAN/) {
    print "$_";
    $date = substr($_, 0, 10);
    print "First date: $date\n";
    last;
  }
}

$tlmfilecnt = 0; $impducnt = 0;
while(<ID>) {
  if(/^$YEARSCAN/) {
    print "$_";
    $newdate = substr($_, 0, 10);
    if($date ne $newdate) {
      print "tlm file count = $tlmfilecnt  IM_PDU count = $impducnt\n";
      print "New Date: $newdate\n";
      $tlmfilecnt = 0; $impducnt = 0;
      $date = $newdate;
    }
  }
  if(/^\*\*Processed/) {
     $tlmfilecnt++;
#    print "$_";
    while(<ID>) {		#this is 'compete images ' line
#      print "$_";
      $pos1 = index($_, "and ");
      $pos1 = $pos1+4;
      $pos2 = index($_, "VCDUs");
      $vcducnt = substr($_, $pos1, ($pos2-$pos1));
      $impducnt += $vcducnt;
      #print "!!TEMP vcducnt = $vcducnt impducnt = $impducnt\n";
      last;
    }
  }
}
print "tlm file count = $tlmfilecnt  IM_PDU count = $impducnt\n";

close(ID);

