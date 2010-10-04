eval 'exec /home/jsoc/bin/$JSOC_MACHINE/perl -S $0 "$@"'
    if 0;
#lev1_definitive_db.pl
#
#Determine if can now run the definitive lev1 for the given
#ordinal date in UTC (YYYY.DDD_UTC). For example for July 12, 2010
#the day of year is 193 and to determine if the lev1 definitive can
#be run for this day and put the answers in the DB table hmi.lev1_probe.
#NOTE: Does not run anything. Just inserts or updates a DB record.
#Call:
#     lev1_definitive_db.pl [-o] hmi,aia 2010.193_UTC
#
#These requirements must be met:
#
#1. For the current UTC day# (1-366), there must be a file like so:
#/dds/soc2pipe/hmi,aia/xday/Xmit_All.198
#Optionally, the -o[verride] flag can be given for older data before
#the Xmit_All.ddd was implemented.
#Sets the DB column ReXmit to 1 (or 2 if -o).
#
#2. Do a query for the day of interest like so:
#sdo.hk_dayfile[2010.07.14_00:00:00_UTC][129][moc]
#where MERGED=1.
#If this exists, then the sdo.lev0_asd_0004 for this day is complete,
#and ok to use for the definitive lev1.
#Sets the DB column ASD to 1.
#
#3. Check sdo.fds_orbit_vectors for a record for the given ord date.
#Sets the DB column FDS to 1.
#
#4. Check that the flatfield has a flatfield_version >= 1.
#Check for camera 1&2 for hmi, and for the 23 wave_str for aia.
#Sets the DB column FF to 1.
#
#5. Check that sdo.master_pointing record exists for our date.
#Sets the DB column MP to 1.
#
#6. Check that hmi,aia.temperature_summary_300s exists for our date.
#Sets the DB column TEMP to 1.
#
#7. Check that the manually set No Go file does not exist:
#/home/production/cvs/JSOC/proj/lev0/data/NoDefLev1[HMI,AIA]
#Sets the DB column GOFLG to 1 if file does not exist.
#
#When the lev1_def_gui calls build_lev1_mgr for an ordinal date, 
#the 'executed' column is set in the DB table hmi.lev1_probe.
#
#
######################################################################

use DBI;

$DB = jsoc;
$HOSTDB = "hmidb";      #host where DB runs
$PGPORT=5432;

$XmitFlgDirHMI = "/dds/soc2pipe/hmi";
$XmitFlgDirAIA = "/dds/soc2pipe/aia";
$NoGoFileHMI = "/home/production/cvs/JSOC/proj/lev0/data/NoDefLev1HMI";
$NoGoFileAIA = "/home/production/cvs/JSOC/proj/lev0/data/NoDefLev1AIA";
#$DSFFNAMEHMI = "su_production.hmi_flatfield";
$DSFFNAMEHMI = "hmi.flatfield";
$DSFFNAMEAIA = "aia.flatfield";
#$DSFFNAMEAIA = "aia_test.flatfield";
$DSINAIA = "aia.lev0";
#$DSOUTAIA = "aia.lev1e_nrt";	#!!TEMP
#$DSOUTAIA = "su_production.aia_lev1e_nrt";	#!!TEMP
$DSOUTAIA = "aia.lev1";
$DSINHMI = "hmi.lev0a";
#$DSOUTHMI = "hmi.lev1c_nrt";
$DSOUTHMI = "hmi.lev1";
#$QSUBDIR = "/scr21/production/qsub/tmp";
$EXECUTE = 0;
$OVERRIDE = 0;
$| = 1;

@allowhost = ("cl1n001", "cl1n002", "cl1n003"); #hosts w/dcs[0,1] mounts
@wavestr = ("131_OPEN", "131_THICK", "131_THIN", "1600", "1700", "171_THICK", "171_THIN", "193_OPEN", "193_THICK", "193_THIN", "211_OPEN", "211_THICK", "211_THIN", "304_OPEN", "304_THICK", "304_THIN", "335_OPEN", "335_THICK", "335_THIN", "4500", "94_OPEN", "94_THICK", "94_THIN");


sub usage {
  print "Determine if the lev1 definitive can be made for a given ordinal date(s) and updates the DB table hmi.lev1_probe.\n";
  print "Usage: lev1_definitive_db.pl [-o] hmi,aia 2010.193_UTC 2010.198_UTC\n";
  print " where  -o = overried the check for Xmit_All.ddd (used for older data)\n";
  print "             NOTE: -o s/b used for days before 2010.182_UTC (1July2010)\n";
  exit(0);
}

$user = $ENV{'USER'};
if($user ne "production") {
  print "You must be user production to run\n";
  exit;
}
$host = `hostname -s`;
if(!grep($host, @allowhost)) {
  print "Can only be run on host with dcs[0,1] mounts: @allowhost\n";
  exit(0);
}

while ($ARGV[0] =~ /^-/) {
  $_ = shift;
  if (/^-o(.*)/) {
    $OVERRIDE = 1;
  }
}

if(!($#ARGV == 1 || $#ARGV == 2)) {
  &usage;
}
$instru = $ARGV[0];
if($instru ne "hmi" && $instru ne "aia") { &usage; }
if($instru eq "hmi") {
  $hmiaiaflg = 0;
} else {
  $hmiaiaflg = 1;
}
$orddate1 = $ARGV[1];
if($#ARGV == 2) { $orddate2 = $ARGV[2]; }
else { $orddate2 = 0; }

$pos1 = index($orddate1, '.');
$pos2 = index($orddate1, '_');
if($pos1==-1 || $pos2==-1) {
  &usage;
}
$yr1 = substr($orddate1, 0, $pos1);
$yrday1 = substr($orddate1, $pos1+1, ($pos2-$pos1)-1);
$zone1 = substr($orddate1, $pos2+1);
if($zone1 ne "UTC") { &usage; }

if($orddate2) {
  $pos1 = index($orddate2, '.');
  $pos2 = index($orddate2, '_');
  if($pos1==-1 || $pos2==-1) {
    &usage;
  }
  $yr2 = substr($orddate2, 0, $pos1);
  $yrday2 = substr($orddate2, $pos1+1, ($pos2-$pos1)-1);
  $zone2 = substr($orddate2, $pos2+1);
  if($zone2 ne "UTC") { &usage; }
  if($yrday2 < $yrday1) {
    print "Second ordinal date must be >= the first\n";
    &usage;
  }
}

#connect to database
  $dbh = DBI->connect("dbi:Pg:dbname=$DB;host=$HOSTDB;port=$PGPORT", "$user", "$password");
  if ( !defined $dbh ) {
    die "Cannot do \$dbh->connect: $DBI::errstr\n";
  }

$utcdate = &inittoday();	#set $todayUTCdayofyr 1-366 & ret date
$sec = `time_convert ord=$orddate1`;
chomp($sec);
if($orddate2) {
  $sec2 = `time_convert ord=$orddate2`;
  chomp($sec2);
} else {
  $sec2 = $sec;
}

$fulldate = `time_convert s=$sec zone=UTC`;
$nextsec = $sec + 86400;
$nextdate = `time_convert s=$nextsec zone=UTC`;
chomp($fulldate);
chomp($nextdate);
print "Run on $utcdate which is day of year $todayUTCdayofyr\n"; 
print "for lev1 on ordinal date starting at $orddate1 ($fulldate)\n";

for(; $sec <= $sec2; $sec = $sec + 86400) {

$fulldate = `time_convert s=$sec zone=UTC`;
chomp($fulldate);
$orddate = `time_convert s=$sec o=ord`;
chomp($orddate);
$pos1 = index($orddate, '.');
$pos2 = index($orddate, '_');
$yr = substr($orddate, 0, $pos1);
$yrday = substr($orddate, $pos1+1, ($pos2-$pos1)-1);
$zone = substr($orddate, $pos2+1);
print "\nProcessing: $orddate ($fulldate):\n";

if($hmiaiaflg) {
  $xmitfile = "$XmitFlgDirAIA/xday/Xmit_All.$yrday";
} else {
  $xmitfile = "$XmitFlgDirHMI/xday/Xmit_All.$yrday";
}
#print "Checking for $xmitfile\n";
if(-e $xmitfile) {
  $OKfile = 1;
  print "*OK for lev1: $xmitfile exists\n";
}
else {
  $OKfile = 0;
  print "**NG for lev1: $xmitfile does not exist\n";
  if($OVERRIDE) {
    print "*OK override flag set for $xmitfile\n";
    $OKfile = 2;
  }
}

$taidate = `time_convert s=$sec zone=TAI`;
#print "TAI date = $taidate\n";
$pos = index($taidate, '_');
$querypart = substr($taidate, 0, $pos+1);
#$qtime = $querypart."00:00:00_TAI"; #originally was TAI in the DB
$qtime = $querypart."00:00:00_UTC";
#print "query time = $qtime\n";
#print "!!TEMP force query time to 2010.07.15_00:00:00_TAI\n";
#$qtime = "2010.07.15_00:00:00_TAI";
$qsec = `time_convert time=$qtime`;
chomp($qsec);
#print "$qsec\n";
$cmd = "echo \"select sunum from sdo.hk_dayfile where obs_date=$qsec and merged=1\" | psql -h hmidb jsoc";
@select = `$cmd`;
#print "select is:\n@select";
shift(@select); shift(@select);  #elim hdr info
$ans = shift(@select);
print "sdo.hk_dayfile answer = $ans";
if($ans =~ /(0 rows)/) {
  print "**NG sdo.hk_dayfile MERGED not found\n";
  $OKhkdayfile = 0;
}
else {
  print "*OK sdo.hk_dayfile MERGED found\n";
  $OKhkdayfile = 1;
}


#Now determine if a sdo.fds_oribit_vectors record exists for our ord date.
$cmd = "echo \"select obs_date from sdo.fds_orbit_vectors where obs_date=$sec\" | psql -h hmidb jsoc";
@select = `$cmd`;
shift(@select); shift(@select);  #elim hdr info
$ans = shift(@select);
print "sdo.fds_oribit_vectors answer = $ans";
if($ans =~ /(0 rows)/) {
  print "**NG fds_orbit_vector not found\n";
  $OKfdsorbit = 0;
}
else {
  print "*OK fds_orbit_vector found\n";
  $OKfdsorbit = 1;
}

$OKff = 1;
if(!$hmiaiaflg) {		#ck for hmi flat field
  #Now see if flat field has flatfield_version >= 1
  for($i=1; $i < 3; $i++) {
    $query = sprintf("%s[? t_start <= \$(%s) and t_stop > \$(%s) and CAMERA=%d and flatfield_version=(select max(flatfield_version) from %s where t_start <= \$(%s) and t_stop > \$(%s) and CAMERA=%d) ?]", $DSFFNAMEHMI, $fulldate, $fulldate, $i, $DSFFNAMEHMI, $fulldate, $fulldate, $i);
    #print "hmi query= $query\n"; #!!TEMP
    #print "Must put single quote around the above\n";
    $cmd = "show_info key=date,flatfield_version '$query'";
    #print "$cmd\n";
    @result = `$cmd`;
    #print "Result of flatfield query for $orddate:\n";
    print "@result";
    $x = shift(@result);
    if($x =~ /date\tflatfield_version/) {
      $x = shift(@result);	#looks like 2010-07-01T17:28:23Z   1
      ($a, $ffver) = split(/\s+/, $x);
      if($ffver >= 1) { 
        print "*OK flatfield found for CAMERA=$i\n";
      } else {
        print "**NG flatfield_version not >= 1  for CAMERA=$i\n";
        $OKff = 0;
      }
    }
    else {
      print "**NG flatfield not found for CAMERA=$i\n";
      $OKff = 0;
    }
  }
}
else {				#ck for aia flat field
  while($wave = shift(@wavestr)) {
    $query = sprintf("%s[? t_start <= $sec and t_stop > $sec and WAVE_STR='%s' and flatfield_version=(select max(flatfield_version) from %s where t_start <= $sec and t_stop > $sec and WAVE_STR='%s') ?]", $DSFFNAMEAIA, $wave, $DSFFNAMEAIA, $wave);
    #print "\naia query= $query\n";
    #print "Must put double quote around the above\n";
    $cmd = "show_info key=date,flatfield_version \"$query\"";
    #print "$cmd\n";
    @result = `$cmd`;
    print "@result";
    $x = shift(@result);
    if($x =~ /date\tflatfield_version/) {
      $x = shift(@result);	#looks like 2010-07-01T17:28:23Z   1
      ($a, $ffver) = split(/\s+/, $x);
      if($ffver >= 1) { 
        print "*OK flatfield found for WAVE_STR=$wave\n";
      } else {
        print "**NG flatfield_version not >= 1  for WAVE_STR=$wave\n";
        print "**OVERRIDE: temp acceptance of flatfield_version not >= 1\n";
        #$OKff = 0;
      }
    }
    else {
      print "**NG flatfield not found for WAVE_STR=$wave\n";
      $OKff = 0;
    }
  }
}

#Check that sdo.master_pointing record exists for our date
$cmd = "show_info key=t_start,t_stop 'sdo.master_pointing[? t_stop > \$($fulldate) ?]'";
@result = `$cmd`;
print "@result";
$x = shift(@result);
if($x =~ /t_start\tt_stop/) {
  print "*OK sdo.master_pointing record found\n";
  $OKmp = 1;
}
else {
  print "**NG sdo.master_pointing record not found\n";
  $OKmp = 0;
}

#ck that hmi.temperature_summary_300s records exist
$pos = index($fulldate, '_');
$fdf = substr($fulldate, 0, $pos+1);
$full_last_5 = $fdf."23:55:00.00_UTC";
$cmd = "show_info -c 'hmi.temperature_summary_300s[? t_start >= \$('$fulldate') and t_start <= \$('$full_last_5') ?]'";
$cnt = `$cmd`;
if($cnt =~ /^288 records/) {
  print "*OK hmi.temperature_summary_300s 288 records found\n";
  $OKtemp = 1; 
}
else {
  print "**NG hmi.temperature_summary_300s incomplete (s/b 288 records):\n";
  print "$cnt\n";
  $OKtemp = 0; 
}

#Check if no go file exists.
if($hmiaiaflg) {
  $nogo = $NoGoFileAIA;
} else {
  $nogo = $NoGoFileHMI;
}
if(-e $nogo) {
  print "The no-go file exist: $nogo\n";
  $OKnogofile = 0;
}
else { $OKnogofile = 1; }

$sql = "select count(*) from hmi.lev1_probe where ord_date='$orddate'";
print "$sql\n";
$sth = $dbh->prepare($sql);
  if ( !defined $sth ) {
    print "Cannot prepare statement: $DBI::errstr\n";
    $dbh->disconnect();
    exit;
  }
    # Execute the statement at the database level
    $sth->execute;
    # Fetch the rows back from the SELECT statement
    while(@row = $sth->fetchrow()) {
      $row = shift(@row);
      print "# of rows for $orddate =  $row\n";
    }
    if(defined $sth) {
      $sth->finish;
    }

    if($row == 0) {	#insert new row, else update old one
      $sql = "SELECT NEXTVAL('hmi.lev1_probe_seq')";
      $sth = $dbh->prepare($sql);
      if ( !defined $sth ) {
        print "Cannot prepare statement: $DBI::errstr\n";
        $dbh->disconnect();
        exit;
      }
      $sth->execute;
      # Fetch the rows back from the SELECT statement
      $recnum = $sth->fetchrow();
      if(defined $sth) {
        $sth->finish;
      }
      $sql = "insert into hmi.lev1_probe (recnum, ORD_DATE,CREATE_DATE,ReXmit,ASD,FDS,FF,MP,TEMP,GOFLG,EXECUTED) values ($recnum, '$orddate', '$utcdate', $OKfile, $OKhkdayfile, $OKfdsorbit, $OKff, $OKmp, $OKtemp, $OKnogofile, 0)";
      $sth = $dbh->prepare($sql);
      if ( !defined $sth ) {
        print "Cannot prepare statement: $DBI::errstr\n";
        $dbh->disconnect();
        exit;
      }
      $sth->execute;
    }
    else {
      $sql = "update hmi.lev1_probe set CREATE_DATE='$utcdate', ReXmit=$OKfile, ASD=$OKhkdayfile, FDS=$OKfdsorbit, FF=$OKff, MP=$OKmp, TEMP=$OKtemp, GOFLG=$OKnogofile, EXECUTED=0 where ord_date='$orddate'";
      #print "!!!TEMP\n$sql\n";
      $sth = $dbh->prepare($sql);
      if ( !defined $sth ) {
        print "Cannot prepare statement: $DBI::errstr\n";
        $dbh->disconnect();
        exit;
      }
      $sth->execute;
    }

}

$dbh->disconnect();
print "**END: lev1_definitive.pl\n";  #used by lev1_def_gui to know end

#Initialize $todayUTCdayofyr when first start
sub inittoday {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2);
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
  $todayUTCdayofyr = $yday+1;   #1-366
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $year4 = sprintf("%04d", $year+1900);
  $date = $year4.".".$mon2.".".$mday2._.$hour2.":".$min2.":".$sec2."_UTC";
  return($date);
}
