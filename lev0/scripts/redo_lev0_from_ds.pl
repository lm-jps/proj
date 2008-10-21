#!/usr/bin/perl
#redo_lev0_from_ds.pl
#
#!!!NOTE: this must be run to stage all the tlm/qacx file to $DEST
#before the ingest_lev0 processing is started.
#
#!!!TBD paramatize this fro VC, dirs, FSN, etc. so can use
#as general purpose script.
#
#Taken from redo_lev0_from_logs.pl and changed to query the lev0
#for all relevant .tlm/.qac files and reprocess these.
#
use Term::ReadKey;

$DEST = "/dds/soc2pipe/hmi";	#place to put the .tlm and .qacx files
#$FSNSTART = 1760955;		#!!TBD make an input arg
#$FSNSTOP = 1760958;		#!!TBD make an input arg
$FSNSTART = 1760959;
$FSNSTOP = 1771472;

sub usage {
  print "Redo lev0 from a query of the lev0 TLMDSNAM for the raw files\n";
  print "Usage: redo_lev0_from_ds.pl\n";
  print "       Requires hmi password to run\n";
  exit(1);
}

if($#ARGV != -1 || $ARGV[0] eq "-h") {
  &usage;
}
#$IDTBL = $ARGV[0];
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

$savfile = "xxx";
$cmd = "show_info 'hmi.lev0e[$FSNSTART-$FSNSTOP]' key=FSN,TLMDSNAM";
@ans = `$cmd`;
$x = shift(@ans);
print "$x"; 			#get label line
while($x = shift(@ans)) {
  print "$x";
  ($fsn, $tlmds) = split(/\s+/, $x);
  $cmd = "show_info -P -q $tlmds";
  $dir = `$cmd`;
  #print "$dir\n";
  chomp($dir);
  $pos1 = index($tlmds, '[');
  $pos1++;
  $pos2 = index($tlmds, ']');
  $file = substr($tlmds, $pos1, $pos2-$pos1);
  $cmd = "select sunum from hmi.tlme where filename='$file'";
  $pcmd = "echo \"$cmd\" | psql -q -h hmidb jsoc";
  @sunums = `$pcmd`;
  $x = shift(@sunums);		#skip sunum label
  $x = shift(@sunums);		#skip ----------
  $insel = " ";
  while($x = shift(@sunums)) {
    if($x =~ /^\(/) {		#e.g. (4 rows) at end of query
      last;
    }
    chomp($x);
    $insel = sprintf("%s%s,", $insel, $x);
  }
  chop($insel);			#get rid of last comma
  #print "insel = $insel\n";
  $cmd = "select online_loc from sum_main where ds_index in ($insel)";
  $pcmd = "echo \"$cmd\" | psql -q -h hmidb -p 5434 jsoc_sums";
  @online = `$pcmd`;
  $x = shift(@online);          #skip online_loc label
  $x = shift(@online);          #skip ----------
  while($x = shift(@online)) {
    if($x =~ /^\(/) {		#e.g. (4 rows) at end of query
      last;
    }
    chomp($x);			#wd like /SUM12/D9066779
    $cmd = "cp -p $x/S00000/*.tlm $DEST";
    print "$cmd\n";
    if(system($cmd)) {
      print "ERROR: on $cmd\n";
      exit(1);
    }
    $qfile = "$x/S00000/*.qac*";
    $from = `ls $qfile`;
    chomp($from);
    $pos = rindex($from, '.');
    $file = substr($from, 0, $pos);
    $pos = rindex($file, "VC");
    $file = substr($file, $pos);
    $to = "$DEST/$file."."qacx";
    $cmd = "cp -p $from $to\n";
    print "$cmd\n";
    if(system($cmd)) {
      print "ERROR: on $cmd\n";
      exit(1);
    }

  }

}
