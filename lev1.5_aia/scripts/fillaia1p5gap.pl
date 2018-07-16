#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
exit if -e "/tmp22/jps/stopProcessing";
$ENV{SUMSERVER} = "k1";
$ENV{SGE_ROOT} = "/SGE";
my $JSOC_MACHINE = $ENV{JSOC_MACHINE};
my ($yr, $mon, $hr, $da, $month, $year, $to, $t0, @t_obs0, @t_obs1, $qt);
my (%to0, %to1, $t_beg, $t_end);
my $mn = 0; my $nhr = 1;
my $logroot = "/tmp29/jps/logs/aia1p5/gapfill";
my $pnam = "/home/jps/bin/fillaia1p5gap.csh";
my $show_info = "/home/jsoc/cvs/Development/JSOC/bin/$JSOC_MACHINE/show_info";
$t0 = time();
($hr, $da, $month, $year) = (gmtime($t0))[2, 3, 4, 5];
$yr = $year + 1900; $mon = $month + 1;
GetOptions(
           "year=i" => \$yr,
           "month=i" => \$mon,
           "day=i" => \$da,
           "hour=i" => \$hr,
           "minute=i" => \$mn,
           "nh=i" => \$nhr
          );
$qt = sprintf "$yr.%2.2d.%2.2d_%2.2d:%2.2d_UTC/${nhr}h", $mon, $da, $hr, $mn;
@t_obs0 = `$show_info -q key=t_obs ds=\'aia.lev1_nrt2[$qt][?QUALITY>-1?]\'`;
@t_obs1 = `$show_info -q key=t_obs ds=\'aia_test.lev1p5[$qt][?QUALITY>-1?]\'`;
chomp @t_obs0; chomp @t_obs1;
foreach $to (@t_obs0) { $to =~ tr/\-T/._/; $to0{$to} = $to; } # hash of src
foreach $to (@t_obs1) { $to =~ tr/\-T/._/; $to1{$to} = $to; } # dest hash
$qt = '';                  # empty list of missing destination records
my $missing = 0;           # =1 when in range of missing records
foreach $to (@t_obs0) {
  if ($missing) {
    if ($to1{$to}) {       # 1 past end of range of missing records
       $qt .= $t_beg;      # append range of missing T_OBS to query string
       $qt .= "-$t_end" unless $t_beg eq $t_end;
       $missing = 0;
    } else {               # record still missing, update end of range
      $t_end = $to;
    }
  } else {                 # no new missing record found yet
    next if $to1{$to};     # this record not missing either, skip
    $t_beg = $t_end = $to; # begin new range of missing records, end = begin
    $qt .= "," if $qt;     # need a comma if previous range exists
    $missing = 1;
  }
}
if (substr($qt, length($qt) - 1, 1) eq ',') { chop $qt; }
my $cmd = "/SGE/bin/lx24-amd64/qsub2 -q p.q -o $logroot.log -e $logroot.err $pnam $qt";
print "$qt\n";
if ($qt) { system $cmd; } else { print "No missing records.\n"; }
