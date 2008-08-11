#!/usr/bin/perl
#/home/production/cvs/JSOC/proj/datacapture/scripts/rsync_prod.pl
#
#Usage: rsync_prod.pl [log_file]
#where log_file is an optional file to store the output log.
#      If none is given then the log goes to
#      /tmp/rsync_prod_YYYY.MM.DD_HH:MM:SS.log
#
#rsync's dcs0:/home/production/cvs/JSOC to d00:/home/production/dcs0_backup
#
#Typically this is run as a user production cron job like so:
#15 8-18 * * * /home/production/cvs/JSOC/proj/datacapture/scripts/rsync_prod.pl
#

#Return date in form for a label e.g. 1998.01.07_14:42:04
sub labeldate {
  local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst,$date,$sec2,$min2,$hour2,$mday2,$year2);
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  $sec2 = sprintf("%02d", $sec);
  $min2 = sprintf("%02d", $min);
  $hour2 = sprintf("%02d", $hour);
  $mday2 = sprintf("%02d", $mday);
  $mon2 = sprintf("%02d", $mon+1);
  $year4 = sprintf("%04d", $year+1900);
  $date = $year4.".".$mon2.".".$mday2._.$hour2.":".$min2.":".$sec2;
  return($date);
}

$| = 1;			#flush output as we go
if($#ARGV == -1) {
  $label = &labeldate;
  $logfile = "/tmp/rsync_prod_".$label;
}
else {
  $logfile = $ARGV[0];
}
#print "logfile = $logfile\n";

#set up for ssh w/o password
#$cmd = "source /home/production/cvs/JSOC/proj/datacapture/scripts/ssh_rsync_prod.source";

$sshfile = "/tmp/ssh-agent.env";	#made by ssh-agent setup procedure
$sshfilesource = "/tmp/ssh-agent.env.source";
if(!-e $sshfile) {
  print "Error: Can't find $sshfile\n";
  print "See log: $logfile\n";
  exit(1);
}
open (EV,"$sshfile") || die "Can't Open : $sshfile $!\n";
open (EVSRC,">$sshfilesource") || die "Can't Open : $sshfilesource $!\n";
while(<EV>) {
  if(/^setenv/) {
    ($a, $b, $arg) = split(/ /);
    print EVSRC "$b=$arg\n"; 
    print EVSRC "export $b\n";
  }
}
close(EV);
close(EVSRC);
$cmd = "source $sshfilesource; /usr/bin/rsync --rsh=/usr/bin/ssh --rsync-path=/usr/bin/rsync -avz /home/production/cvs/JSOC d00:/home/production/dcs0_backup";

print "$cmd\n";
if(system "$cmd 1> $logfile 2>&1") {
  print "Error on: $cmd\n";
  print "See log: $logfile\n";
  exit(1);
}

$cmd = "source $sshfilesource; /usr/bin/rsync --rsh=/usr/bin/ssh --rsync-path=/usr/bin/rsync -avz /srv/www d00:/home/production/dcs0_backup";

print "$cmd\n";
if(system "$cmd 1>> $logfile 2>&1") {
  print "Error on: $cmd\n";
  print "See log: $logfile\n";
  exit(1);
}
