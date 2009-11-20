eval 'exec /home/jsoc/bin/$JSOC_MACHINE/perl -S $0 "$@"'
    if 0;
#/home/production/cvs/JSOC/proj/datacapture/scripts/rd_full_tape_info_input.pl
#
#Takes the output file created by rd_full_tape_info.pl (typically rd_info.log)
#and gets the tape file#, md5 checksum and ds_index values for each file of
#the tape that was read ($TAPEID gives the tapeid to do).
#Updates the sum_main and sum_file jsocdc DB tables accordingly.
#This must be run on either dcs0 or dcs1 where the given tape was originally made.
#This procedure is to fix the db for the correct info on the tape. The db got
#messed up because the mt -f /dev/nstx eod gave the wrong value for the eod. 
#The driven.c code has been fixed to not give a warning but now an error if
#this happens. So far only tape 015020L4 has this problem.
#Must be run on dcs1. If tape was first created on dcs0, then change $HNAME.
#
#Here's what an input file looks like:
#drive type = Generic SCSI-2 tape
# drive status = 1174405120
# sense key error = 0
# residue count = 0
# file number = 1
# block number = 0
# Tape block size 0 bytes. Density code 0x46 (unknown).
# Soft error count since last status=0
# General status bits on (81010000):
#  EOF ONLINE IM_REP_EN
#
# D696815/
# D696815/VC02_2009_131_14_58_50_00137629c58_1c298_00.qac
# D696815/VC02_2009_131_14_58_50_00137629c58_1c298_00.tlm
# D696816/
# D696816/VC05_2009_131_14_59_19_201376246e8_1c298_00.qac
# D696816/VC05_2009_131_14_59_19_201376246e8_1c298_00.tlm
# D696817/
# D696817/VC02_2009_131_14_59_53_00137645ef0_1c298_00.qac
# D696817/VC02_2009_131_14_59_53_00137645ef0_1c298_00.tlm
#
#md5=75e1c50616c68df2cf0b2bf19f175349

use DBI;

sub usage {
  print "Update the jsocdc DB with file from rd_full_tape_info.pl output\n";
  print "Usage: rd_full_tape_info_input.pl rd_info.log\n";
  exit(1);
}

if($infile = shift(@ARGV)) {
}
else {
  &usage;
  exit;
}

#!!!NOTE: Change these according to the tape that you're fixing
$HNAME = "dcs1";
$TAPEID = "015020L4";
$DB = "jsocdc";

$HOSTDB = `hostname -s`;
chomp($HOSTDB);
if($HOSTDB ne $HNAME) {
  print "Must be run on $HNAME\n";
  exit;
}
$hostdb = $HOSTDB;      #host where Postgres runs
$user = $ENV{'USER'};
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

#first get rid of old sum_file entry for this tapeid
  $sqlcmd = "delete from sum_file where TAPEID = '$TAPEID'";
  print "$sqlcmd\n";
  $sth = $dbh->prepare($sqlcmd);
  if ( !defined $sth ) {
    die "Cannot prepare statement: $DBI::errstr\n";
  }
  # Execute the statement at the database level
  $sth->execute;


open(IN, $infile) || die "Can't open $infile: $!\n";
$i = 0;
while(<IN>) {
  if(/^#/ || /^\n/) { #ignore any comment or blank lines
    next; 
  }
  if(/file number =/) {
    ($a, $filenum) = split(/= /);
    chomp($filenum);
    print "filenum = $filenum\n";  #!!TEMP
    next;
  }
  if(/^ D/) {
    chomp();
    if(/\/$/) {
      chop();
      $dsindex[$i] = substr($_, 2);
      print "dsindex = $dsindex[$i]\n";
      $i++;
    }
    next;
  }
  if(/md5=/) {
    chomp();
    ($a, $md5) = split(/=/);
    print "md5 = $md5\n";
  } 
  else { next; }
  #now have all the info to update the DB. Continue...
  if($i == 0) {
    print "Error: No ds_index can be found\n";
    exit;
  }
  for($j = 0; $j < $i; $j++) {
    $sqlcmd = "update sum_main set Arch_Tape='$TAPEID', Arch_Tape_FN=$filenum where ds_index=$dsindex[$j]";
    print "$sqlcmd\n";
    $sth = $dbh->prepare($sqlcmd);
    if ( !defined $sth ) {
      die "Cannot prepare statement: $DBI::errstr\n";
    }
    # Execute the statement at the database level
    $sth->execute;
  }
  $sqlcmd = "insert into sum_file (tapeid,filenum,gtarblock,md5cksum) values ('$TAPEID',$filenum,256,'$md5')";
  print "$sqlcmd\n";
  $sth = $dbh->prepare($sqlcmd);
  if ( !defined $sth ) {
    die "Cannot prepare statement: $DBI::errstr\n";
  }
  # Execute the statement at the database level
  $sth->execute;
 

  $i = 0;	#reset before get next complete entry

}
  #now put last file# in sum_tape
  $filenum++;
  $sqlcmd = "update sum_tape set nxtwrtfn=$filenum where Tapeid='$TAPEID'";
  print "$sqlcmd\n";
  $sth = $dbh->prepare($sqlcmd);
  if ( !defined $sth ) {
    die "Cannot prepare statement: $DBI::errstr\n";
  }
  # Execute the statement at the database level
  $sth->execute;
$dbh->disconnect();

