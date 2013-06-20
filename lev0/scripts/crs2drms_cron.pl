#!/usr/bin/perl
use strict;
use warnings;

#this file is: cvs/JSOC/proj/lev0/scripts/crs2drms_cron.pl
#called by cron to detect new CRS and ingest into DRMS

$ENV{SUMSERVER} = "j1";
$ENV{JSOC_DBUSER} = "production";
my $series = 'iris_ground.window_table';
my $stage = '/scr21/jps/IRIS';
my $seg = "${stage}/CRS";
my $scriptdir = '/home/jsoc/cvs/JSOC/proj/lev0/scripts';
foreach my $flag (glob "$stage/CRS/[1-9]*") {
  unlink $flag;
  my $ser_num = `basename $flag`; chomp $ser_num;
  my $crs_root = sprintf "$stage/CRS/CRS-%.5d", $ser_num;
  my $dc_nam = sprintf "$stage/CRS/crop%d", $ser_num;
  my @xmls = glob "${crs_root}*.xml";
  warn "Too many or not enough XML for DC ser #" if $#xmls;
  my $out_base = `basename $xmls[0] .xml`; chomp $out_base;
# `/home/jps/crs2fig_twig.pl -ifile=$xmls[0] > $seg/$out_base.fig`;
  `$scriptdir/crs2svg_twig.pl -ifile=$xmls[0] > $seg/$out_base.svg`;
  `$scriptdir/crs2drms_twig.pl -ifile=$xmls[0] -series=$series -seg=$seg`;
}
