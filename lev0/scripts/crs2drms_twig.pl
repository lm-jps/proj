#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use XML::Twig;
use Getopt::Long;
$ENV{"TZ"} = "UTC";

#this file is: cvs/JSOC/proj/lev0/scripts/crs2drms_twig.pl
#called by:    cvs/JSOC/proj/lev0/scripts/crs2drms_cron.pl
#cereates and populates new record in CRS table series

my ($inpname, $dcrname);
my $seg_base = '.';
my $set_info = '/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/set_info -c';
my $series = 'iris_ground.crs_table';
my $reg_max = 0;
my ($dt, $crs_id, $crs_desc, $crs_type, $sumsptrl, $sumspat, $numregn, $srid);
my (@StartRow, @EndRow, @StartCol, @EndCol);
my (@CRS_SR, @CRS_ER, @CRS_SC, @CRS_EC);
my (@Tsr, @Ter, @Tsc, @Tec, $winflip, $nwin);
my ($sc, $mn, $hr, $da, $mo, $yr, $crsid, $nr, $seg, $sptrlstr, $spatstr);

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
             "seg_path=s" => \$seg_base,
             "series=s" => \$series,
             "set_info=s" => \$set_info
           );
if (!defined $inpname) {
    warn "Mandatory argument --ifile=inputname.xml not supplied";
    exit;
}

foreach (@ARGV) { print "Don't understand argument '$_', missing '-'?\n"; }

($sc, $mn, $hr, $da, $mo, $yr) = gmtime(time);
$dt = sprintf "DATE=%d.%2.2d.%2.2d_%2.2d:%2.2d:%2.2d",
      $yr+1900, $mo+1, $da, $hr, $mn, $sc;

$twig->parsefile($inpname);

my $ds = "ds=" . $series;
my $cmd = join (' ',  $set_info, $ds, $dt, $crs_id, $crs_desc, $crs_type,
                $numregn, $sptrlstr, $spatstr, $nwin);

if ($crs_type =~ /FUV/) {
  $winflip = "WIN_FLIP=1";
  for (my $i=1;  $i<=$nr;  $i++) {
    my $sc = int(($StartRow[$i] -1)/$sumsptrl) + 1;
    $Tsc[$i] = sprintf "TSC%d=%d", $i, $sc;
    my $ec = int(($EndRow[$i] -1)/$sumsptrl) + 1;
    $Tec[$i] = sprintf "TEC%d=%d", $i, $ec;
    my $sr = int((1096 - $EndCol[$i])/$sumspat) + 1;
    $Tsr[$i] = sprintf "TSR%d=%d", $i, $sr;
    my $er = int((1096 - $StartCol[$i])/$sumspat) + 1;
    $Ter[$i] = sprintf "TER%d=%d", $i, $er;
  }
} elsif ($crs_type =~ /NUV/) {
  $winflip = "WIN_FLIP=2";
  for (my $i=1;  $i<=$nr;  $i++) {
    my $sc = int((4144 - $EndRow[$i])/$sumsptrl) + 1;
    $Tsc[$i] = sprintf "TSC%d=%d", $i, $sc;
    my $ec = int((4144 - $StartRow[$i])/$sumsptrl) + 1;
    $Tec[$i] = sprintf "TEC%d=%d", $i, $ec;
    my $sr = int(($StartCol[$i] - 1)/$sumspat) + 1;
    $Tsr[$i] = sprintf "TSR%d=%d", $i, $sr;
    my $er = int(($EndCol[$i] - 1)/$sumspat) + 1;
    $Ter[$i] = sprintf "TER%d=%d", $i, $er;
  }
} elsif ($crs_type =~ /SJI/) {
  my $sjinuv = 0;
  my ($sc, $ec, $sr, $er);
  if ($StartRow[1] < 1037) { $sjinuv = 1; }
  for (my $i=1;  $i<=$nr;  $i++) {
    if ($sjinuv) {
      $winflip = "WIN_FLIP=0";
      $sc = int(($StartRow[$i] - 1)/$sumsptrl) + 1;
      $ec = int(($EndRow[$i] - 1)/$sumsptrl) + 1;
    } else {
      $winflip = "WIN_FLIP=2";
      $sc = int((2072 - $EndRow[$i])/$sumsptrl) + 1;
      $ec = int((2072 - $StartRow[$i])/$sumsptrl) + 1;
    }
    $Tsc[$i] = sprintf "TSC%d=%d", $i, $sc;
    $Tec[$i] = sprintf "TEC%d=%d", $i, $ec;
    $sr = int(($StartCol[$i] - 1)/$sumspat) + 1;
    $Tsr[$i] = sprintf "TSR%d=%d", $i, $sr;
    $er = int(($EndCol[$i] - 1)/$sumspat) + 1;
    $Ter[$i] = sprintf "TER%d=%d", $i, $er;
  }
}
for (my $i=$reg_max+1; $i<9; $i++) {
  $CRS_SR[$i] = sprintf "CRS_SR%d=0", $i;
  $CRS_ER[$i] = sprintf "CRS_ER%d=0", $i;
  $CRS_SC[$i] = sprintf "CRS_SC%d=0", $i;
  $CRS_EC[$i] = sprintf "CRS_EC%d=0", $i;
  $Tsr[$i] = sprintf "TSR%d=0", $i;
  $Ter[$i] = sprintf "TER%d=0", $i;
  $Tsc[$i] = sprintf "TSC%d=0", $i;
  $Tec[$i] = sprintf "TEC%d=0", $i;
}
for (my $i=1;  $i<9; $i++) {
  $cmd = join (' ',  $cmd, $winflip, $CRS_SR[$i], $CRS_ER[$i], $CRS_SC[$i],
               $CRS_EC[$i], $Tsr[$i], $Ter[$i], $Tsc[$i], $Tec[$i]);
}

$seg = sprintf "decrop_table=%s/crop%d", $seg_base, $crsid;
$cmd = join (' ',  $cmd, $seg);
$seg = sprintf "crs_xml=%s", $inpname;
$cmd = join (' ',  $cmd, $seg);
#my $fignam = `basename $inpname .xml`; chomp $fignam;
#$fignam = "$seg_base/$fignam" . ".fig";
#$seg = sprintf "fig_image=%s", $fignam;
#$cmd = join (' ',  $cmd, $seg);
my $svgnam = `basename $inpname .xml`; chomp $svgnam;
$svgnam = "$seg_base/$svgnam" . ".svg";
$seg = sprintf "svg_image=%s", $svgnam;
$cmd = join (' ',  $cmd, $seg);

#printf "Cmd: '%s'\n", $cmd;
system $cmd;

sub Info_Desc
{ my ( $twig_arg, $Info_Desc) = @_;
  $crs_desc = 'CRS_DESC="' . $Info_Desc->text . '"';
}

sub Info_Type
{ my ( $twig_arg, $Info_Type) = @_;
  my $ct = $Info_Type->text;
  $ct =~ s/1/v/;
  $ct =~ s/fuv|nuv|sji/\U$&/;
  $crs_type = 'CRS_TYPE=' . $ct;
}

sub Header_Id
{ my ( $twig_arg, $Header_Id) = @_;
  $crsid = $Header_Id->text;
  $crs_id =  'IICRSID=' . $crsid;
}

sub Header_Size
{ my ( $twig_arg, $Header_Size) = @_;
  
}

sub Header_sr
{ my ( $twig_arg, $Header_sr) = @_;
  $nr = $Header_sr->text;
  $numregn = 'CRS_NREG=' . $nr;
  $nwin = 'NWIN=' . $nr;
  if ($nr > $reg_max) { $reg_max = $nr; }
}

sub Header_spec
{ my ( $twig_arg, $Header_spec) = @_;
  $sumsptrl = $Header_spec->text;
  $sptrlstr = 'SUMSPTRL=' . $sumsptrl;
}

sub Header_spat
{ my ( $twig_arg, $Header_spat) = @_;
  $sumspat = $Header_spat->text;
  $spatstr = 'SUMSPAT=' . $sumspat;
}

sub Data_srid
{ my ( $twig_arg, $Data) = @_;
  $srid = $Data->text;
}

sub Data_sr
{ my ( $twig_arg, $Data) = @_;
  $StartRow[$srid] = $Data->text;
  $CRS_SR[$srid] = sprintf "CRS_SR%d=%d", $srid, $Data->text;
}

sub Data_er
{ my ( $twig_arg, $Data) = @_;
  $EndRow[$srid] = $Data->text;
  $CRS_ER[$srid] = sprintf "CRS_ER%d=%d", $srid, $Data->text;
}

sub Data_sc
{ my ( $twig_arg, $Data) = @_;
  $StartCol[$srid] = $Data->text;
  $CRS_SC[$srid] = sprintf "CRS_SC%d=%d", $srid, $Data->text;
}

sub Data_ec
{ my ( $twig_arg, $Data) = @_;
  $EndCol[$srid] = $Data->text;
  $CRS_EC[$srid] = sprintf "CRS_EC%d=%d", $srid, $Data->text;
}
