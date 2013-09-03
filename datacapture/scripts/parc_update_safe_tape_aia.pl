eval 'exec /home/jsoc/bin/$JSOC_MACHINE/perl -S $0 "$@"'
    if 0;
#/home/prodtest/cvs/JSOC/proj/datacapture/scripts/parc_update_safe_tape_aia.pl
#
#This is run on the datacapture systems as a cronjob every night to 
#see if there are new .parc files in /dds/pipe2soc/[aia,hmi].
#A .parc file entry looks like:
##owning_series_name                         tape_id  fn  date
#VC05_2013_233_06_37_26_02ce0470118_1c298_00 022382L4 1231 2013-08-22 11:53:21
#VC02_2013_233_06_37_41_02c0045e8b8_1c298_00 022382L4 1231 2013-08-22 11:53:21
#
#The tape_id/fn in the .parc file will be entered in the sum_main
#as the safe_tape (i.e. the tape that has the given tlm data on
#it on the backend system).
#
#This use to be done by the do_pipe2soc() in ingest_lev0_hmiaia.c.
#But with the new datacapture systems staring 21Aug2013, the
#owning_series entry in the sum_main table no longer is the file
#name given in the .parc file, but is now the real drms owning_series,
#i.e. hmi.tlm or aia.tlm. So the only way we have of finding this
#tlm file is to query recent entries in the sum_main that do not
#have a safe_tape and match their tlm file name with one in the
#.parc file that we're looking at.
#

use DBI;

sub usage {
  print "Update the safe_tape entry in the sum_main table from .parc files.\n";
  print "Usage: parc_update_safe_tape.pl [-h] [-thrs] \n";
  print "       -h = help\n";
  print "       -t = # of hours+2mins previous to current time to find data\n";
  print "            The default is 24hrs+2mins previous to find data\n";
}

#Return date in form e.g. 2008-05-29 11:38:18
##where this value is $hrsprev+120sec previous to the current time.
##NOTE: this format must match that in the jsocdc DB
##Also sets $currdate to current time of form 2008_135_12_10
sub labeldate {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2,$year2);
  $secprev = ($hrsprev*3600) + 120;     #add extra 2 mins
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

#NOTE:Can't use hostname to figure out which db, as this can
#run on a backup host if the dcs1x is down, so give db name. 
$DB = "aiadb_sums";
$PGPORT = 5434;
$DIR_PARC = "/dds/pipe2soc/aia";
#$DIR_PARC = "/dds/pipe2soc/aia/err"; #!!TEMP for testing
$HOSTDB = `hostname -s`;
chomp($HOSTDB);

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
print "Going to find .parc files and update safe_tape info in sum_main table\n";
$user = $ENV{'USER'};
if($user ne "prodtest") {
  print "You must be user prodtest to run\n";
  exit;
}
#First connect to database
  $dbh = DBI->connect("dbi:Pg:dbname=$DB;host=$HOSTDB;port=$PGPORT", "$user", "$password");
  if ( !defined $dbh ) {
    die "Cannot do \$dbh->connect: $DBI::errstr\n";
  }


$sql = "select online_loc from sum_main where creat_date >= '$ldate' and owning_series='aia.tlm' and safe_tape is NULL";

    $sth = $dbh->prepare($sql);
    if ( !defined $sth ) {
      print "Cannot prepare statement: $DBI::errstr\n";
      $dbh->disconnect();
      exit;
    }
    # Execute the statement at the database level
    $sth->execute;
    # Fetch the rows back from the SELECT statement
    @wd = (); 
    @ix = (); 
    $count = 0;
    while ( $wd = $sth->fetchrow() ) {
      $pos = rindex($wd, '/D');
      $ix = substr($wd, $pos+2);
      push(@wd, $wd);
      push(@ix, $ix);
    }

@parcfiles = `/bin/ls $DIR_PARC/*.parc`;
@parclines = ();
while($pfile = shift(@parcfiles)) {
  open(PF, $pfile) || die "Can't open $pfile: $!\n";
  while(<PF>) {
    if(/^#/ || /^\n/) { #ignore any comment or blank lines
      next;
    }
    push(@parclines, $_);
  }
}
  $i = 0;
  while($wd = @wd[$i]) {
    #print "wd = $wd  i=$i\n";
    #print "ix = @ix[$i]\n";
    $ix = @ix[$i];
    $cmd = "/bin/ls $wd/S00000/*.tlm";
    #print "$cmd\n";
    @tlm = `$cmd`;
    #print "@tlm\n";
    if(!@tlm) {
      print "ERROR: no dir or .tlm for $wd/S00000\n";
      $i++;
      next;
    }
    $tlm = @tlm[0];
    $pos = rindex($tlm, '/VC');
    if($pos == -1) {
      print "ERROR: no VC file in $wd/S00000\n";
      $i++;
      next;
    }
    $tlm = substr($tlm, $pos+1);
    $pos = rindex($tlm, '.tlm');
    $file = substr($tlm, 0, $pos);
    chomp($file);
    #print "$file\n";	#the name that will be in the .parc
    #push(@targetfiles, $file);
    $i++;
    if(@pline = grep(/$file/, @parclines)) {
      while($pline = shift(@pline)) { 
        ($vcfile, $tapeid, $fn, $date, $hrmin) = split(/\s+/, $pline);
        $fulldate = sprintf("%s %s", $date, $hrmin);
        $sql = "update sum_main set safe_tape='$tapeid', safe_tape_fn=$fn, safe_tape_date='$fulldate' where ds_index = $ix";
        print "$vcfile\n";
        print "$sql\n";
        $sth = $dbh->prepare($sql);
        if ( !defined $sth ) {
          print "Cannot prepare statement: $DBI::errstr\n";
          $dbh->disconnect();
          exit;
        }
        # Execute the statement at the database level
        $sth->execute;
      }
    }
  }

$dbh->disconnect();
