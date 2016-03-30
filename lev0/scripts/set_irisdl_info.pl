#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Time::Local;
use Getopt::Long;

my ($xband_log_nam, $vcdu_nam, $cmd, $cmd2, $mm, $aos_dt, $frstdt, $lastdt);
my ($pass, @passes, $yr, $mo, $da, $hr, $mn, $sc, $wd, $yd, $dst, $aos);
my ( $tlm, $sums, $full, $numd, @words);
my ($orb, $sta, $md, $hms, $nd, $ntx, $nrx, $nmis, $pctmis, $fnam, $miss);
my ($t_beg, $t_end, $fsn_b, $fsn_e, $npkts, $seq_b, $seq_e);

my $set_info = '/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/set_info -c';
my $series = 'lm_jps.iris_dl'; my $series2 = 'lm_jps.iris_dl_ids_test';
my $verbose = 1;
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
  $aos_dt = "$yr.$mm.${da}_$hr:$mn:${sc}_UTC";
  my $aos_dtkw = "AOS=$yr.$mm.${da}_$hr:$mn:${sc}_UTC";
  $miss = "NUM_MISS=\'$nmis $pctmis\'";
  my $orbit = sprintf "ORBIT=%d", $orb;
  my $station = "STA=$sta";
  $numd = "ND=$nd";
  my $tx = "TX=$ntx"; my $rx = "RX=$nrx";
  if ($fnam) {
    my $file="FILE=$fnam";
    $tlm = substr $fnam, 0, 22;
    my $name = "NAME=$tlm";
    next if `show_info -q -c \"$series\[$orb\]\[$aos_dt\]\[?NAME=\'$tlm\'?\]\"` > 0;
    $sums = `show_info -q -P iris.tlm[$tlm]`; chomp $sums;
    if ($sums) {
      my ($s0, $s1, $s2, $s3);
      $full = "$sums/$fnam";
      my $mtime = (stat $full)[9];
      my @rows = `vcdu_time_range -c < $full`;
      $cmd2 = "$set_info ds=$series2 $orbit $aos_dtkw $name";
      for (my $i=0; $i<20; $i++) {
         if ($i<=$#rows) { @words = split /\s+/, $rows[$i]; }
         else { $words[0]=0; $words[1]=0; }
        $cmd2 = join (' ',  $cmd2, "ID_$i=$words[0] NUM_$i=$words[1]");
      }
      system join (' ', $cmd2 , $dt);
      my $line = `vcdu_time_range < $full`;
      ($t_beg,$t_end,$fsn_b,$fsn_e,$npkts,$seq_b,$seq_e,$s0,$s1,$s2,$s3) =
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
      my $isysns = "ISYSN0=$s0 ISYSN1=$s1 ISYSN2=$s2 ISYSN3=$s3";
      my $fsnb_age = sprintf "FSNB_AGE=%d", $aos - $t_beg;
      my $fsne_age = sprintf "FSNE_AGE=%d", $aos - $t_end;
      $cmd = join (' ',  $set_info, $ds, $orbit, $station, $aos_dtkw, $frstdt,
           $lastdt, $numd, $tx, $rx, $miss, $npckts, $frst_fsn, $last_fsn,
           $isysns, $name, $file, $fsnb_age, $fsne_age, $dt);
    } else {
    $cmd = join (' ',  $set_info, $ds, $orbit,  $station, $aos_dtkw, $numd,
         $tx, $rx, $miss, $name, $file,  $dt);
    }
  } else {
    next if `show_info -q -c $series\[$orb\]\[$aos_dt\]` > 0;
    $cmd = join (' ',  $set_info, $ds, $orbit,  $station, $aos_dtkw, $numd,
         $tx, $rx, $miss, $dt,  "NAME=NONE");
  }
  system $cmd;
  if ($verbose > 0) { print "$orb "; }
  elsif ($verbose > 1) { printf "%s\n", substr $cmd, 30, 78; }
}
