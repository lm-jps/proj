#!/usr/bin/perl
#/home/production/cvs/JSOC/proj/datacapture/scripts/dcstab.pl
#
#Display or change the file:
#/home/production/cvs/JSOC/proj/datacapture/scripts/dsctab.txt
#which assigns the AIA and HMI processing to a datacapture machine
#hostname like so:
#
#AIA=dcs0
#HMI=dcs1
#
#The offsite dcs3 machine is typically assigned like so:
#AIA=dcs3
#HMI=dcs3

sub usage {
  print "Display or change the datacapture system assignment file.\n";
  print "Usage: dcstab [-h][-l][-e]\n";
  print "       -h = print this help message\n";
  print "       -l = list the current file contents\n";
  print "       -e = edit with vi the current file contents\n";
  exit(1);
}

$DCSTABDIR = "/home/production/cvs/JSOC/proj/datacapture/scripts";
$DCSTABFILE = "$DCSTABDIR/dcstab.txt";

$list = $editflg = 0;

while ($ARGV[0] =~ /^-/) {
  $_ = shift;
  if (/^-h(.*)/) {
    &usage;
    exit(0);
  }
  if (/^-l(.*)/) {
    $list = 1;
  }
  if (/^-e(.*)/) {
    $editflg = 1;
  }
}
if($list && $editflg) {
  print "Error: Only one option (-l or -e) can be given\n";
  exit(1);
}
if(!$list && !$editflg) {
  &usage;
  exit(1);
}
if($list) {
  @table = `cat $DCSTABFILE`;
  while($x = shift(@table)) {
    print "$x";
  }
  exit(0);
}
$cmd = "vi $DCSTABFILE";
system($cmd);

$thishost = `hostname -s`;
chomp($thishost);
print "This host = $thishost\n";
if($thishost eq 'dcs3') {
  print "Assignment complete on dcs3\n";
}

if($thishost eq 'dcs0') {
  @hosts = (dcs1,dcs2);
}
elsif($thishost eq 'dcs1') {
  @hosts = (dcs0,dcs2);
}
elsif($thishost eq 'dcs2') {
  @hosts = (dcs0,dcs1);
}

print "\nNow attempt to replicate this file other datacapture machines\n";
while($x = shift(@hosts)) {
  $cmd = "scp $DCSTABFILE $x:$DCSTABDIR";
  print "$cmd\n";
  system($cmd);
}
