#!/usr/bin/perl
#below has problem when run in cron job
#eval 'exec /home/jsoc/bin/$JSOC_MACHINE/perl -S $0 "$@"'
#    if 0;
#
#/home/production/cvs/JSOC/proj/datacapture/scripts/build_arc_to_dds.pl
#
#This make a .arc (archive status file) to be sent from the datacapture system 
#to DDS system. The .arc file is put in /dds/soc2dds/[aia,hmi]/arc. 
#A .arc file tells the DDS what tlm files have been successfully archive on the 
#datacapture system.
#
#This is typically run on the datacapture system as a cron job after midnight. 
#It wil query the sum_main table to find all the archived info for the
#.tlm data since the 24 hour period prior to the time of this script
#running. (e.g. if run after midnight on Wed, will find all the data 
#archived since midnight of Tues).

#!!!TBD: Also check that safe_tape and Offsite_Ack are true in the 
#        sum_main before can send ack to DDS.
#
#The .arc file built looks like:
#FILE_NAME=VC02_2008_182_19_57_55_0014211db08_1c298_00.tlm
#FILE_NAME=VC05_2008_182_19_57_56_0014211db08_1c298_00.tlm
#EOF_MARKER=C5C5
#
#NOTE: aia is always assumed to be on dsc0, and hmi on dcs1
#
use DBI;

sub usage {
  print "Make a .arc file of archive info to send to the DDS.\n";
  print "Usage: build_arc_to_dds.pl [-h] [-thrs] \n";
  print "       -h = help\n";
  print "	-t = # of hours+2mins previous to current time to find data\n";
  print "            The default is 24hrs+2mins previous to find data\n";
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
  $yday++;
  $currdate = $year4."_".$yday."_".$hour2."_".$min2;
  return($date);
}

$HOSTDB = `hostname -s`;
chomp($HOSTDB);
$DB = "jsocdc";
$ARC_ROOT_HMI = "/dds/soc2dds/hmi/arc/";
$ARC_ROOT_AIA = "/dds/soc2dds/aia/arc/";
$hrsprev = 24;
while ($ARGV[0] =~ /^-/) {
  $_ = shift;
  if (/^-h(.*)/) {
    &usage;
    exit(0);
  }
  if (/^-t(.*)/) {
    $hrsprev = $1;
  }
}

$ldate = &labeldate();
$parchmi = $ARC_ROOT_HMI."HMI_$currdate".".arc";
$parcaia = $ARC_ROOT_AIA."AIA_$currdate".".arc";
print "Going to find tlm data archived since $ldate\n";

$hostdb = $HOSTDB;      #host where Postgres runs
$user = $ENV{'USER'};
#$ENV{'SUMSERVER'} = "d02.Stanford.EDU";
if($user ne "production") {
  print "You must be user production to run\n";
  exit;
}
if(!($pgport = $ENV{'SUMPGPORT'})) {
  print "You must have ENV SUMPGPORT set to the port number, e.g. 5430\n";
  exit; 
}

#First connect to database
  $dbh = DBI->connect("dbi:Pg:dbname=$DB;host=$hostdb;port=$pgport", "$user", "$password");
  if ( !defined $dbh ) {
    die "Cannot do \$dbh->connect: $DBI::errstr\n";
  }

  if($hostdb eq "dcs1") {
    open(PARC, ">$parchmi") || die "Can't open $parchmi: $!\n";
    $resultfile = $parchmi; $instru = "HMI";
    #$sql = "select owning_series from sum_main where archive_status='Y' and Arch_Tape_Date >= '$ldate' and (owning_series like 'VC02%' or owning_series like 'VC05%' or owning_series like 'VC10%' or owning_series like 'VC13%') and safe_tape!='' and offsite_ack='Y'";
    $sql = "select owning_series from sum_main where archive_status='Y' and Arch_Tape_Date >= '$ldate' and (owning_series like 'VC02%' or owning_series like 'VC05%' or owning_series like 'VC10%' or owning_series like 'VC13%')";
  } else {
    open(PARC, ">$parcaia") || die "Can't open $parcaia: $!\n";
    $resultfile = $parcaia; $instru = "AIA";
    #$sql = "select owning_series from sum_main where archive_status='Y' and Arch_Tape_Date >= '$ldate' and (owning_series like 'VC01%' or owning_series like 'VC04%' or owning_series like 'VC09%' or owning_series like 'VC12%') and safe_tape!='' and offsite_ack='Y'";
    $sql = "select owning_series from sum_main where archive_status='Y' and Arch_Tape_Date >= '$ldate' and (owning_series like 'VC01%' or owning_series like 'VC04%' or owning_series like 'VC09%' or owning_series like 'VC12%')";
  }

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
      $file = $row.".tlm";
      print PARC "FILE_NAME=$file\n";
      #print "FILE_NAME=$file\n";
      $count++;
    }
  print PARC "EOF_MARKER=C5C5\n";
  close(PARC);
  print "\nDone: Found $count tlm files for $instru\n";
  print "Results in $resultfile\n";
  #`cp $resultfile /dds/soc2dds/save`; #!!TEMP save elsewhere
$dbh->disconnect();

