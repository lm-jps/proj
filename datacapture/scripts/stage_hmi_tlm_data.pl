eval 'exec /home/jsoc/bin/$JSOC_MACHINE/perl -S $0 "$@"'
    if 0;
#/home/production/cvs/JSOC/proj/datacapture/scripts/stage_hmi_tlm_data.pl
#
#Takes a spool file created by a query of the jsoc db that looks like:
#psql -h hmidb jsoc
#\o /tmp/hmi.lev0e.tlmdsnam.spool
#select tlmdsnam,fsn  from hmi.lev0e where fsn >=1800000 and fsn <= 1879935 order by fsn;
#
#            tlmdsnam             |   fsn 
#---------------------------------+---------
# hmi.tlme[VC05_2008_251_06_21_45] | 1800000
# hmi.tlme[VC02_2008_251_06_21_47] | 1800001
# hmi.tlme[VC05_2008_251_06_21_45] | 1800002
# hmi.tlme[VC02_2008_251_06_21_47] | 1800003
#
#and finds the tlmdsname file name in the jsocdc db on dcs1 (for hmi).
#Will copy the .tlm and .qac file from the found storage unit, if it is
#online, into /dds/stage.
#
use DBI;

sub usage {
  print "Copy hmi tlm data to /dds/stage\n";
  print "Usage: stage_hmi_tlm_data.pl jsoc_spool_file\n";
  exit(1);
}

#Return date in form e.g. 2008-05-29 11:38:18
#where this value is $hrsprev+120sec previous to the current time.
#NOTE: this format must match that in the jsocdc DB
#Also sets $currdate to current time of form 2008_135_12_10
sub labeldate {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2,$year2);
  $secprev = ($hrsprev*3600) + 120;	#add extra 2 mins
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time-$secprev);
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $year4 = sprintf("%04d", $year+1900);
  $date = $year4."-".$mon2."-".$mday2." ".$hour2.":".$min2.":".$sec2;
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $year4 = sprintf("%04d", $year+1900);
  $currdate = $year4."_".$yday."_".$hour2."_".$min2;
  return($date);
}

if($infile = shift(@ARGV)) {
}
else {
  print "Must give a jsoc_spool_file for input\n";
  exit;
}

$HOSTDB = `hostname -s`;
chomp($HOSTDB);
$DB = "jsocdc";
$STAGE_HMI = "/dds/stage";

$ldate = &labeldate();
$hostdb = $HOSTDB;      #host where Postgres runs
$user = $ENV{'USER'};
if($user ne "production") {
  print "You must be user production to run\n";
  exit;
}

#First connect to database
  $dbh = DBI->connect("dbi:Pg:dbname=$DB;host=$hostdb", "$user", "$password");
  if ( !defined $dbh ) {
    die "Cannot do \$dbh->connect: $DBI::errstr\n";
  }

open(IN, $infile) || die "Can't open $infile: $!\n";
while(<IN>) {
  if(/^#/ || /^\n/) { #ignore any comment or blank lines
    next; 
  }
  ($a, $tlmdsname) = split(/\s/);
  print "tlmdsname=$tlmdsname\n";
  $pos = index($tlmdsname, '[');
  $file = substr($tlmdsname, $pos+1, 22);
  print "\$file=$file\n";

    $sql = "select online_loc from sum_main where online_status='Y' and owning_series like '$file%'";
    print "\$sql = $sql\n";
    $sth = $dbh->prepare($sql);
    if ( !defined $sth ) {
      print "Cannot prepare statement: $DBI::errstr\n";
      $dbh->disconnect();
      exit; 
    }
    # Execute the statement at the database level
    $sth->execute;
    # Fetch the rows back from the SELECT statement
    @row = ();
    $count = 0;
    while ( $row = $sth->fetchrow() ) {
      print "online_loc = $row\n";
      $lastrow = $row;
      $count++;
    }
    if($count) { 
      $source = $lastrow."/*";
      $cmd = "cp $source $STAGE_HMI";
      print "$cmd\n";
      `$cmd`;
    }
}
$dbh->disconnect();

