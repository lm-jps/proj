#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use XML::Twig;
use Getopt::Long;
$ENV{"TZ"} = "UTC";

#this file is: cvs/JSOC/proj/lev0/scripts/crs2svg_twig.pl
#called by:    cvs/JSOC/proj/lev0/scripts/crs2drms_cron.pl
#creates SVG figure of CRS to be ingested into CRS table series

my ($inpname, $dcrname);
my ($dt, $crs_id, $crs_desc, $crs_type, $sumsptrl, $sumspat, $numregn, $srid);
my (@StartRow, @EndRow, @StartCol, @EndCol, @Height, @Width);
my (@CRS_SX, @CRS_ER, @CRS_SY, @CRS_EC, @CRS_H, @CRS_W);
my (@Tsx, @Ter, @Tsy, @Tec, @T_H, @T_W);
my ($sc, $mn, $hr, $da, $mo, $yr, $crsid, $nr, $seg); # , $tx);
my $reg_max = 0;
my $seg_base = '.';

my $twig= new XML::Twig(twig_handlers => { 'Info/Description' => \&Info_Desc,
                                           'Info/Type' => \&Info_Type,
                                           'Header/Id' => \&Header_Id,
                                           'Header/Size' => \&Header_Size,
                                           'Header/Subregions' => \&Header_sr,
                                           'Header/Spectral' => \&Header_spec,
                                           'Header/Spatial' => \&Header_spat,
                                           'Data/SubregionId' => \&Data_srid,
                                           'Data/StartRow' => \&Data_sr,
                                           'Data/EndRow' => \&Data_er,
                                           'Data/StartCol' => \&Data_sc,
                                           'Data/EndCol' => \&Data_ec
                                         } );
GetOptions(
             "ifile=s" => \$inpname,
             "dcrname=s" => \$dcrname,
             "seg_path=s" => \$seg_base
           );
if (!defined $inpname) {
    warn "Mandatory argument --ifile=inputname.xml not supplied";
    exit;
}

foreach (@ARGV) { print "Don't understand argument '$_', missing '-'?\n"; }

($sc, $mn, $hr, $da, $mo, $yr) = gmtime(time);
$dt = sprintf "%d.%2.2d.%2.2d_%2.2d:%2.2d:%2.2d",
      $yr+1900, $mo+1, $da, $hr, $mn, $sc;

print "<?xml version=\"1.0\" standalone=\"no\"?>\n";
print "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n";
print "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
print "<!-- Creator: crs2svg_twig.pl -->\n";
print "<!-- CreationDate: $dt -->\n";
print "<!-- Magnification: 1.000 -->\n";
print "<svg    xmlns=\"http://www.w3.org/2000/svg\"\n";
print "        xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n";
print "        width=\"7.0in\" height=\"4.5in\"\n";
#print "        width=\"3.5in\" height=\"1.1in\"\n";
print "        viewBox=\"-12 -12 4200 2700\">\n";
print "<g style=\"stroke-width:.025in; fill:none\">\n";

$twig->parsefile($inpname);

if ($crs_type =~ /nuv/) { $CRS_SX[0] = 2073; } else { $CRS_SX[0] = 1; }
$Tsx[0] = 1;
$CRS_SY[0] = 1; $CRS_H[0] = 1096;
$Tsy[0] = 1411; $T_H[0] = int (1096/$sumspat);
if ($crs_type =~ /fu1/) {
  $CRS_W[0] = 4144;
  $T_W[0] = int(4144/$sumsptrl);
} else {
  $CRS_W[0] = 2072;
  $T_W[0] = int(2072/$sumsptrl);
}

for (my $i=1;  $i<=$reg_max; $i++) {
  $CRS_SX[$i] = $StartRow[$i];
  $CRS_W[$i] = $Width[$i];
  $CRS_SY[$i] = 1097 - $EndCol[$i];
  $CRS_H[$i] = $Height[$i];
}

if ($crs_type =~ /fu1/) {
  for (my $i=1;  $i<=$reg_max; $i++) {
    $Tsx[$i] = int(($StartRow[$i] - 1)/$sumsptrl) + 1;
    $T_W[$i] = int($Width[$i]/$sumsptrl);
    $Tsy[$i] = int(($StartCol[$i] - 1)/$sumspat) + 1411;
    $T_H[$i] = int($Height[$i]/$sumspat);
  }
} elsif ($crs_type =~ /nuv/) {
  for (my $i=1;  $i<=$reg_max; $i++) {
    $Tsx[$i] = int((4144 - $EndRow[$i])/$sumsptrl) + 1;
    $T_W[$i] = int($Width[$i]/$sumsptrl);
    $Tsy[$i] = int((1096 - $EndCol[$i])/$sumspat) + 1411;
    $T_H[$i] = int($Height[$i]/$sumspat);
  }
} elsif ($crs_type =~ /sji/) {
  my $sjinuv = 0;
  if ($CRS_SX[1] < 1037) { $sjinuv = 1; }
  my ($sc, $ec, $sr, $er);
  for (my $i=1;  $i<=$reg_max; $i++) {
     if ($sjinuv) {
      $sc = int(($StartRow[$i] - 1)/$sumsptrl) + 1;
    } else {
      $sc = int((2072 - $EndRow[$i])/$sumsptrl) + 1;
    }
    $Tsx[$i] = $sc;
    $T_W[$i] = int($Width[$i]/$sumsptrl);
    $Tsy[$i] = int((1096 - $EndCol[$i])/$sumspat) + 1411;
    $T_H[$i] = int(($Height[$i])/$sumspat);
  }
}

for (my $i=0;  $i<=$reg_max; $i++) {
  printf "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"0\"\n",
      $CRS_SX[$i], $CRS_SY[$i], $CRS_W[$i], $CRS_H[$i];
  if ($i) { print "style=\"stroke:#000000;stroke-width:3;\n"; }
  else { print "style=\"stroke:#000000;stroke-width:7;\n"; }
  print "stroke-linejoin:miter; stroke-linecap:butt;\n";
  print "stroke-dasharray:10 30;" unless $i;
  print "\"/>\n";
  printf "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"0\"\n",
      $Tsx[$i], $Tsy[$i], $T_W[$i], $T_H[$i];
  if ($i) { print "style=\"stroke:#000000;stroke-width:3;\n"; }
  else { print "style=\"stroke:#000000;stroke-width:7;\n"; }
  print "stroke-linejoin:miter; stroke-linecap:butt;\n";
  print "stroke-dasharray:10 30;" unless $i;
  print "\"/>\n";
}
my $txt_beg = "<text xml:space=\"preserve\" x=\"0\" y=\"1332\"";
my $txt_mid = "fill=\"#000000\"  font-family=\"Times\" font-style=\"normal\" font-weight=\"normal\" font-size=\"96\" text-anchor=\"start\">";
print join ' ', $txt_beg, $txt_mid, $crs_desc, "</text>\n";
$txt_beg = "<text xml:space=\"preserve\" x=\"0\" y=\"1212\"";
my $fnam = `basename $inpname`; chomp $fnam;
my $fstr = sprintf "$fnam ( Summing: %d x %d )", $sumsptrl, $sumspat;
print join ' ', $txt_beg, $txt_mid, $fstr, "</text>\n";
print "</g>\n</svg>\n";

sub Info_Desc
{ my ( $twig_arg, $Info_Desc) = @_;
  $crs_desc = $Info_Desc->text;
}

sub Info_Type
{ my ( $twig_arg, $Info_Type) = @_;
  $crs_type = $Info_Type->text;
}

sub Header_Id
{ my ( $twig_arg, $Header_Id) = @_;
  $crsid = $Header_Id->text;
}

sub Header_Size
{ my ( $twig_arg, $Header_Size) = @_;

}

sub Header_sr
{ my ( $twig_arg, $Header_sr) = @_;
  $nr = $Header_sr->text;
  if ($nr > $reg_max) { $reg_max = $nr; }
}

sub Header_spec
{ my ( $twig_arg, $Header_spec) = @_;
  $sumsptrl = $Header_spec->text;
}

sub Header_spat
{ my ( $twig_arg, $Header_spat) = @_;
  $sumspat = $Header_spat->text;
}

sub Data_srid
{ my ( $twig_arg, $Data) = @_;
  $srid = $Data->text;
}

sub Data_sr
{ my ( $twig_arg, $Data) = @_;
  $StartRow[$srid] = 1*$Data->text;
}

sub Data_er
{ my ( $twig_arg, $Data) = @_;
  $EndRow[$srid] = 1*$Data->text;
  $Width[$srid] = $EndRow[$srid] - $StartRow[$srid] + 1;
}

sub Data_sc
{ my ( $twig_arg, $Data) = @_;
  $StartCol[$srid] = 1*$Data->text;
}

sub Data_ec
{ my ( $twig_arg, $Data) = @_;
  $EndCol[$srid] = 1*$Data->text;
  $Height[$srid] = $EndCol[$srid] - $StartCol[$srid] + 1;
}
