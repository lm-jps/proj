#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Time::Local;
use Getopt::Long;

my ($xband_log_nam, $vcdu_nam, $cmd, $mm, $aos_dt, $frstdt, $lastdt);
my ($pass, @passes, $yr, $mo, $da, $hr, $mn, $sc, $wd, $yd, $dst, $aos);
my ( $tlm, $sums, $full, $numd);
my ($orb, $sta, $md, $hms, $nd, $ntx, $nrx, $nmis, $pctmis, $fnam, $miss);
my ($t_beg, $t_end, $fsn_b, $fsn_e, $npkts, $seq_b, $seq_e);

my $set_info = '/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/set_info -c';
my $series = 'lm_jps.iris_dl';
GetOptions(
            "xband_log_nam=s"  => \$xband_log_nam,
            "series=s" => \$series,
            "set_info=s" => \$set_info
          );
my $ds = "ds=" . $series;
($sc, $mn, $hr, $da, $mo, $yr) = gmtime(time); $yr += 1900; $mo += 1;
my $dt = "DATE=$yr.$mo.${da}_$hr:$mn:${sc}_UTC";
if ($xband_log_nam) { open STDIN, "<<$xband_log_nam" or die "Can't dup: $!"; }
(<>);(<>);(<>);
while (<>) {
  chomp;
  unshift @passes, $_;
}
foreach my $pass (@passes) {
  ($orb, $sta, $md, $hms, $nd, $ntx, $nrx, $nmis, $pctmis, $fnam) =
                                                  split /\s+/, $pass;
  $orb += 0;
  $mm = substr($md, 0, 2); $mo = $mm - 1; $da = substr($md, 3, 2);
  $hr = substr($hms, 0 , 2); $mn = substr($hms, 3, 2);
  $sc = substr($hms, 6, 2); $yr = substr $fnam, 5, 4 if $fnam;
  $aos = timegm($sc, $mn, $hr, $da, $mo, $yr);
  $aos_dt = "AOS=$yr.$mm.${da}_$hr:$mn:${sc}_UTC";
  $miss = "NUM_MISS=\'$nmis $pctmis\'";
  my $orbit = sprintf "ORBIT=%d", $orb;
  my $station = "STA=$sta";
  $numd = "ND=$nd";
  my $tx = "TX=$ntx"; my $rx = "RX=$nrx";
  if ($fnam) {
    my $file="FILE=$fnam";
    $tlm = substr $fnam, 0, 22;
    my $name = "NAME=$tlm";
    $sums = `show_info -q -P iris.tlm[$tlm]`; chomp $sums;
    if ($sums) {
      $full = "$sums/$fnam";
      my $mtime = (stat $full)[9];
      my $line = `vcdu_time_range < $full`;
      ($t_beg, $t_end, $fsn_b, $fsn_e, $npkts, $seq_b, $seq_e) =
         split /\s+/, $line;
      ($sc, $mn, $hr, $da, $mo, $yr, $wd, $yd, $dst) = gmtime($t_beg);
      $yr += 1900; $mo += 1;
      $frstdt = "T_BEG=$yr.$mo.${da}_$hr:$mn:${sc}_UTC";
      ($sc, $mn, $hr, $da, $mo, $yr, $wd, $yd, $dst) = gmtime($t_end);
      $yr += 1900; $mo += 1;
      $lastdt = "T_END=$yr.$mo.${da}_$hr:$mn:${sc}_UTC";
      my $npckts = "NPCKTS=$npkts";
      my $frst_fsn = "FRST_FSN=$fsn_b";
      my $last_fsn = "LAST_FSN=$fsn_e";
      my $fsnb_age = sprintf "FSNB_AGE=%d", $aos - $t_beg;
      my $fsne_age = sprintf "FSNE_AGE=%d", $aos - $t_end;
      $cmd = join (' ',  $set_info, $ds, $orbit, $station, $aos_dt, $frstdt,
           $lastdt, $numd, $tx, $rx, $miss, $npckts, $frst_fsn, $last_fsn,
           $name, $file, $fsnb_age, $fsne_age, $dt);
    } else {
    $cmd = join (' ',  $set_info, $ds, $orbit,  $station, $aos_dt, $numd,
         $tx, $rx, $miss, $name, $file,  $dt);
    }
  } else {
    $cmd = join (' ',  $set_info, $ds, $orbit,  $station, $aos_dt, $numd,
         $tx, $rx, $miss, $dt,  "NAME=none");
  }
  system $cmd;
  printf "%s\n", substr $cmd, 30, 72;
}
