#!/usr/bin/perl -w 

use FindBin qw($Bin);

my($kDEV) = "dev";
my($kGROUND) = "ground";
my($kLIVE) = "live";
my($kCFGPREFIX) = "dlfds_conf";
my($kLOGFILE) = "log";
my($kSCRIPTPATH) = "scrpath";
my($kBINPATH) = "binpath";
my($kPERMDATAPATH) = "permdpath";
my($kDATAPATH) = "dpath";
my($kSERIESNAME) = "series";
my($kNOTLIST) = "notify";

my($jsocmach);

my($arg);
my($argc) = scalar(@ARGV);
my($force) = "";
my($cfgfile);
my($line);

my($logfile);
my($scriptPath);
my($binPath);
my($permdataPath);
my($dataPath);
my($seriesname);
my($notlist);

my($cmd);

if ($argc != 1)
{
    die "Invalid argument list.\n";
}
else
{
    $arg = shift(@ARGV);
    if ($arg eq $kDEV || $arg eq $kGROUND || $arg eq $kLIVE)
    {
	$cfgfile = "$Bin/${kCFGPREFIX}_$arg.txt";
    }
    elsif ($arg eq "-f")
    {
	$force = "-f";
    }
}

open(CFGFILE, "< $cfgfile");

while($line = <CFGFILE>)
{
    chomp($line);

    if ($line =~ /$kLOGFILE\s*=\s*(\S+)\s*/)
    {
	$logfile = $1;
    }

    if ($line =~ /$kSCRIPTPATH\s*=\s*(\S+)\s*/)
    {
	$scriptPath = $1;
    }

    if ($line =~ /$kBINPATH\s*=\s*(\S+)\s*/)
    {
	$binPath = $1;
    }

    if ($line =~ /$kPERMDATAPATH\s*=\s*(\S+)\s*/)
    {
	$permdataPath = $1;
    }

    if ($line =~ /$kDATAPATH\s*=\s*(\S+)\s*/)
    {
	$dataPath = $1;
    }

    if ($line =~ /$kSERIESNAME\s*=\s*(\S+)\s*/)
    {
	$seriesname = $1;
    }

    if ($line =~ /$kNOTLIST\s*=\s*(.+)/)
    {
	$notlist = $1;
    }
}

close(CFGFILE);

if (!defined($ENV{"JSOC_MACHINE"}))
{
    $jsocmach = `$ENV{"JSOCROOT"}/build/jsoc_machine.csh`;
}
else
{
    $jsocmach = $ENV{"JSOC_MACHINE"};
}

local $ENV{"PATH"} = "$binPath/$jsocmach:$ENV{\"PATH\"}";

#print "l=$logfile, sp=$scriptPath, bp=$binPath, pdp=$permdataPath, dp=$dataPath, s=$seriesname, n=$notlist\n";

# Download FDS files to /surge
$cmd = "$scriptPath/dlMOCDataFiles\.pl -c $scriptPath/mocDlFdsSpec.txt -s $permdataPath/mocDlFdsStatus.txt -r $dataPath -t 30 $force";

system("$cmd 1>$logfile 2>&1");

exit(1);

# Ingest FDS files into $seriesname
$cmd = "$scriptPath/fdsIngest\.pl $dataPath -r -s $seriesname";
system("$cmd 1>$logfile 2>&1");

# Notify people who care if an error has occurred
$cmd = "$scriptPath/fdsNotification\.pl -l $logfile -n $notlist";
system($cmd);
