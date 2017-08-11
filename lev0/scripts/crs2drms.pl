#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

$ENV{SUMSERVER} = "k1";
$ENV{JSOC_DBUSER} = "production";
$ENV{JSOC_DBUSER} = "jps";
my $ixml;
my $icropdir = '.';
my $series = 'iris.crs_table';
my $outdir = ".";
GetOptions(
           "icropdir=s" => \$icropdir,
           "ixml=s" => \$ixml,
           "series=s" => \$series
);
if (!defined $ixml) {
  warn "Mandatory argument --ixml=inputname.xml not supplied";
  exit;
}
my $ser_num = substr $ixml, 5, 4;
my $crs_root = sprintf "$outdir/CRS-%.5d", $ser_num;
my $dc_nam = sprintf "$outdir/crop%d", $ser_num;
my $out_base = `basename $ixml .xml`; chomp $out_base;
`/home/jps/crs2svg_twig.pl -ifile=$ixml > $outdir/$out_base.svg`;
`/home/jps/crs2drms_twig.pl -ifile=$ixml -series=$series -seg=$outdir`;
