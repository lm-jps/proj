eval 'exec /home/jsoc/bin/$JSOC_MACHINE/perl -S $0 "$@"'
    if 0;
#
#Determine if can now run the definitive lev1 for the given
#ordinal date in UTC (YYYY.DDD_UTC). For example for July 12, 2010
#the day of year is 193 and to determine if the lev1 definitive can
#be run for this day, call:
#     lev1_definitive.pl [-x][-o] hmi,aia 2010.193_UTC
#
#These requirements must be met:
#
#1. For the current UTC day# (1-366), there must be a file like so:
#/dds/soc2pipe/hmi,aia/xday/Xmit_All.198
#Optionally, the -o[verride] flag can be given for older data before
#the Xmit_All.ddd was implemented.
#
#2. Do a query for the day of interest like so:
#sdo.hk_dayfile[2010.07.14_00:00:00_UTC][129][moc]
#where MERGED=1.
#If this exists, then the sdo.lev0_asd_0004 for this day is complete,
#and ok to use for the definitive lev1.
#
#3. Check sdo.fds_orbit_vectors for a record for the given ord date.
#
#4. Check that the flatfield has a flatfield_version >= 1.
#Check for camera 1&2 for hmi, and for the 23 wave_str for aia.
#
#5. Check that sdo.master_pointing record exists for our date.
#
#6. Check that hmi,aia.temperature_summary_300s exists for our date.
#
#7. Check that the manually set No Go file does not exist:
#/surge40/jsocprod/lev0/data/NoDefLev1[HMI,AIA]
#
#
######################################################################

$XmitFlgDirHMI = "/dds/soc2pipe/hmi";
$XmitFlgDirAIA = "/dds/soc2pipe/aia";
$NoGoFileHMI = "/surge40/jsocprod/lev0/data/NoDefLev1HMI";
$NoGoFileAIA = "/surge40/jsocprod/lev0/data/NoDefLev1AIA";
#$DSFFNAMEHMI = "su_production.hmi_flatfield";
$DSFFNAMEHMI = "hmi.flatfield";
$DSFFNAMEAIA = "aia_test.flatfield";
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
  print "Determine if the lev1 definitive can be made for a given ordinal date.\n";
  print "Usage: lev1_definitive.pl [-x][-o] hmi,aia 2010.193_UTC\n";
  print " where -x = execute the build_lev1_mgr to make data for this day\n";
  print "       -o = overried the check for Xmit_All.ddd (used for older data)\n";
  exit(0);
}

$host = `hostname -s`;
if(!grep($host, @allowhost)) {
  print "Can only be run on host with dcs[0,1] mounts: @allowhost\n";
  exit(0);
}

while ($ARGV[0] =~ /^-/) {
  $_ = shift;
  if (/^-x(.*)/) {
    $EXECUTE = 1;
  }
  if (/^-o(.*)/) {
    $OVERRIDE = 1;
  }
}

if($#ARGV != 1) {
  &usage;
}
$instru = $ARGV[0];
if($instru ne "hmi" && $instru ne "aia") { &usage; }
if($instru eq "hmi") {
  $hmiaiaflg = 0;
} else {
  $hmiaiaflg = 1;
}
$orddate = $ARGV[1];
$pos1 = index($orddate, '.');
$pos2 = index($orddate, '_');
if($pos1==-1 || $pos2==-1) {
  &usage;
}
$yr = substr($orddate, 0, $pos1);
$yrday = substr($orddate, $pos1+1, ($pos2-$pos1)-1);
$zone = substr($orddate, $pos2+1);
if($zone ne "UTC") { &usage; }

$utcdate = &inittoday();	#set $todayUTCdayofyr 1-366 & ret date
$sec = `time_convert ord=$orddate`;
chomp($sec);
$fulldate = `time_convert s=$sec zone=UTC`;
$t_start_sec = $sec + 60;       #add a min to match TAI range in FF rec
$fulldateFF = `time_convert s=$t_start_sec zone=UTC`;
$nextsec = $sec + 86400;
$nextdate = `time_convert s=$nextsec zone=UTC`;
chomp($fulldate);
chomp($fulldateFF);
chomp($nextdate);
print "Run on $utcdate which is day of year $todayUTCdayofyr\n"; 
print "for lev1 on ordinal date = $orddate ($fulldate)\n\n";
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
    $OKfile = 1;
  }
}

$taidate = `time_convert s=$sec zone=TAI`;
#print "TAI date = $taidate\n";
$pos = index($taidate, '_');
$querypart = substr($taidate, 0, $pos+1);
#$qtime = $querypart."00:00:00_TAI"; #originally TAI in the DB
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
#$sec = 1089417600; #!!!TEMP put in bad time
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
    $query = sprintf("%s[? t_start <= \$(%s) and t_stop > \$(%s) and CAMERA=%d and flatfield_version >= 1 ?]", $DSFFNAMEHMI, $fulldateFF, $fulldateFF, $i);
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
    $query = sprintf("%s[? t_start <= $sec and t_stop > $sec and WAVE_STR='%s' and flatfield_version >= 1 ?]", $DSFFNAMEAIA, $wave);
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

if($OKfile && $OKhkdayfile && $OKfdsorbit && $OKnogofile && $OKff && $OKmp && $OKtemp) {
  print "*OK to process definitive lev1 for $orddate\n";
  print "    Building command now (may take a minute)...\n";
  if($hmiaiaflg) {
    $cmd = "show_info key=fsn '$DSINAIA\[? t_obs >= \$($fulldate) and  t_obs < \$($nextdate) ?]'";
  }
  else {
    $cmd = "show_info key=fsn '$DSINHMI\[? t_obs >= \$($fulldate) and  t_obs < \$($nextdate) ?]'";
  }
  #print "$cmd\n";
  @fsn = `$cmd`;
  shift(@fsn);
  $firstfsn = shift(@fsn);
  $lastfsn = pop(@fsn);
  chomp($firstfsn); chomp($lastfsn);
  #$slog = "/tmp/levdef_$utcdate.log";
  #$cmd1 = "$cmd 1> $slog 2>&1";
  #print "$cmd1\n"; #!!TEMP
  #system($cmd1);
  #if(!open(LOGSTAT, $slog)) {
  #  print "Can't open $slog: $!\n";
  #  exit;
  #}
  #<LOGSTAT>;		#skip first line
  #$_ = <LOGSTAT>;
  #chomp();
  #$firstfsn = $_;
  #while(<LOGSTAT>) {
  #  $lastfsn = $_;
  #}
  #close(LOGSTAT);
  #chomp($lastfsn);
  #print "first = $firstfsn  last = $lastfsn\n";  #!!TEMP
  if($hmiaiaflg) {
    $cmd = "time build_lev1_mgr mode=fsn instru=aia dsin=$DSINAIA dsout=$DSOUTAIA bfsn=$firstfsn efsn=$lastfsn logfile=/usr/local/logs/lev1/build_lev1_mgr_aia.$utcdate.log";
  }
  else {
    $cmd = "time build_lev1_mgr mode=fsn instru=hmi dsin=$DSINHMI dsout=$DSOUTHMI bfsn=$firstfsn efsn=$lastfsn logfile=/usr/local/logs/lev1/build_lev1_mgr_hmi.$utcdate.log";
  }
  print "$cmd\n";
  if($EXECUTE) {
    `$cmd`;		#make the data for the day
  }
}
else {
  print "**NG can't process definitive lev1 for $orddate\n";
}
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

