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
$FSNSTART = 1760955;		#!!TBD make an input arg
#$FSNSTOP = 1760958;		#!!TBD make an input arg
#$FSNSTART = 1760959;
$FSNSTOP = 1771472;
#$FSNSTART = 1766530;
#$FSNSTOP = 1766542;

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
#if($passwd ne <passwd>) {
#  print "Invalid passwd\n";
#  exit(1);
#}

$savfile = "xxx";
#$cmd = "show_info 'hmi.lev0e[$FSNSTART-$FSNSTOP]' key=FSN,TLMDSNAM";
#!!TMP just do these know fsn with the 1958 problem
$cmd = "show_info 'hmi.lev0e[1760978,1760984,1761093,1761164,1761169,1761187,1761379,1761382,1761385,1761561,1761574,1761591,1761847,1761866,1761909,1762218,1762224,1762229,1762306,1762310,1762319,1762426,1762430,1762439,1762723,1762919,1762930,1763272,1763567,1763926,1763929,1763940,1763942,1764382,1764875,1765702,1765800,1765808,1765809,1766076,1766211,1766470,1766508,1766511,1766548,1766558,1766727,1766924,1766938,1767139,1767328,1767338,1767476,1767477,1767482,1767726,1767740,1767811,1767938,1767942,1768089,1768276,1768288,1768393,1768800,1768812,1768825,1769136,1769148,1769167,1769282,1769284,1769309,1769317,1769342,1769358,1769359,1769360,1770380,1770673,1770736,1770740,1771062,1771064,1771075,1771256,1771267,1771268,1771440,1771444]' key=FSN,TLMDSNAM";

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
