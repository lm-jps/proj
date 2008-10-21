#!/usr/bin/perl
#
#!!!TBD paramatize this fro VC, dirs, FSN, etc. so can use
#as general purpose script.
#
#use DBI;
use Term::ReadKey;

$DEST = "/dds/soc2pipe/hmi";	#place to put the .tlm and .qacx files

sub usage {
  print "Redo lev0 from the list of log files in the input file.\n";
  print "Usage: redo_lev0_from_logs.pl in_file\n";
  print "       Requires hmi password to run\n";
  print "The in_file looks like:\n";
  print "/usr/local/logs/lev0/VC02_2008.09.09_15:28:49.log\n";
  print "/usr/local/logs/lev0/VC05_2008.09.09_15:28:49.log\n";
  print "/usr/local/logs/lev0/VC02_2008.09.10_15:45:37.log\n";
  print "/usr/local/logs/lev0/VC05_2008.09.10_15:45:37.log\n";
  exit(1);
}

if($#ARGV != 0 || $ARGV[0] eq "-h") {
  &usage;
}
$IDTBL = $ARGV[0];
#print "Need hmi password to run: passwd =";
#ReadMode('noecho');
#$passwd = ReadLine(0);
#chomp($passwd);
#ReadMode('normal');
#print "\n";
#if($passwd ne "hmi4sdo") {
#  print "Invalid passwd\n";
#  exit(1);
#}
$user = "jim";
$password = "jimshoom";
$hostdb = "hmidb";      #host where Postgres runs

open (ID,"<$IDTBL")  or die "can't open $IDTBL: $!";
while(<ID>) {
  if(/^#/ || /^\n/) { #ignore any comment or blank lines
    next;
  }
  print "$_";
  chomp;
  open(LOG, "<$_")  or die "can't open $_: $!";
  $ignore = 1;
  $firstinfile = 1;
  while(<LOG>) {
    if(/^\*\*Processed /) {
      if($ignore == 0) {
        ($a, $tlmfile) = split(/\s+/);
        #print "tlmfile = $tlmfile\n";
        $cmd = "cp -p $tlmfile $DEST";
        print "$cmd\n";
        if(system($cmd)) {
          print "ERROR: on $cmd\n";
          exit(1);
        }
        $pos = rindex($tlmfile, '.');
        $base = substr($tlmfile, 0, $pos);
        $from = $base.".qac*";
        $pos = rindex($base, "VC");
        $file = substr($base, $pos);
        $to = $DEST."/".$file.".qacx";
        $cmd = "cp -p $from $to\n"; 
        print "$cmd\n";
        if(system($cmd)) {
          print "ERROR: on $cmd\n";
          exit(1);
        }
      }
      $ignore = 1;
      $firstinfile = 1;
    }
    elsif(/fsn = /) {
      #print "$_";
      ($a, $fsn) = split(/= /);
      chomp($fsn);
      if($fsn != 92481 && $fsn != 92482) { #this is a good tlm file
        if($firstinfile) {
          if(/^*Closing image for fsn/) {
            if($firstinfile) {
              $firstinfile = 0;
            }
          }
          else {
            print "fsn = $fsn\n";
            $ignore = 0;
          }
        }
      }
    }
  }
  close(LOG);
  #last;		#!!!!TEMP
}
close(ID);

