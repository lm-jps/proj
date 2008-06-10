#!/usr/bin/perl -w 

use FindBin qw($Bin);
use Time::Local;

my($kDEV) = "dev";
my($kGROUND) = "ground";
my($kLIVE) = "live";
my($kCFGPREFIX) = "dllzp_conf";
my($kLOGFILE) = "log";
my($kSCRIPTPATH) = "scrpath";
my($kBINPATH) = "binpath";
my($kPERMDATAPATH) = "permdpath";
my($kDATAPATH) = "dpath";
my($kNAMESPACE) = "ns";
my($kDBUSER) = "dbuser";
my($kDBNAME) = "dbname";
my($kNOTLIST) = "notify";

my($kSPECPREFIX) = "filespec:";
my($kFILESUFFIX) = "[0-9][0-9][0-9]_[0-9][0-9]\\.hkt\\S*";

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
my($namespace);
my($dbuser);
my($dbname);
my($notlist);

my($cmd);

my($fullspec);
my($pspec);  # previous years' days
my($cspec);  # current years' days
my($pbegin); # dofy to begin range, previous year
my($pend);
my($cbegin); # dofy to begin range, current year
my($cend);
my($yr);     # current year
my($dofy);   # current day of year
my($yr90);   # year 90 days ago
my($dofy90); # day of year 90 days ago
my($pyr);    # previous year

my($logcontent) = "";

if ($argc == 0)
{
    die "Invalid argument list.\n";
}
else
{
    while ($arg = shift(@ARGV))
    {
	if ($arg eq "-b")
	{
	    $branch = shift(@ARGV);

	    if ($branch ne $kDEV && $branch ne $kGROUND && $branch ne $kLIVE)
	    {	
		PrintUsage();
		exit(1);
	    }
	}
	elsif ($arg eq "-c")
	{
	    $cfgfile = shift(@ARGV);
	}
	elsif ($arg eq "-f")
	{
	    $force = "-f";
	}
	else
	{
	    PrintUsage();
	    exit(1);
	}
    }
}

if (!defined($cfgfile) && defined($branch))
{
    # Use defalt configuration file
    $cfgfile = "$Bin/${kCFGPREFIX}_$branch.txt";
}

if (!(-f $cfgfile))
{
    printf STDERR "Cannot read configuration file.";
    exit(1);
}

open(CFGFILE, "< $cfgfile");

while($line = <CFGFILE>)
{
    chomp($line);

    if ($line =~ /$kLOGFILE\s*=\s*(\S+)\s*/)
    {
	$logfile = $1;
    }
    elsif ($line =~ /$kSCRIPTPATH\s*=\s*(\S+)\s*/)
    {
	$scriptPath = $1;
    }
    elsif ($line =~ /$kBINPATH\s*=\s*(\S+)\s*/)
    {
	$binPath = $1;
    }
    elsif ($line =~ /$kPERMDATAPATH\s*=\s*(\S+)\s*/)
    {
	$permdataPath = $1;
    }
    elsif ($line =~ /$kDATAPATH\s*=\s*(\S+)\s*/)
    {
	$dataPath = $1;
    }
    elsif ($line =~ /$kNAMESPACE\s*=\s*(\S+)\s*/)
    {
	$namespace = $1;
    }
    elsif ($line =~ /$kDBUSER\s*=\s*(\S+)\s*/)
    {
	$dbuser = $1;
    }
    elsif ($line =~ /$kDBNAME\s*=\s*(\S+)\s*/)
    {
	$dbname = $1;
    }
    elsif ($line =~ /$kNOTLIST\s*=\s*(.+)/)
    {
	$notlist = $1;
    }
}

close(CFGFILE);

local $ENV{"JSOCROOT"} = "$Bin/../../../..";

if (!defined($ENV{"JSOC_MACHINE"}))
{
    $jsocmach = `$ENV{"JSOCROOT"}/build/jsoc_machine.csh`;
    chomp($jsocmach); # back-ticks cause a newline to be appended.
}
else
{
    $jsocmach = $ENV{"JSOC_MACHINE"};
}

# For JSOC modules and scripts, use the paths specified in the config file
local $ENV{"PATH"} = "$scriptPath:$binPath/$jsocmach:$ENV{\"PATH\"}";
local $ENV{"JSOC_DBUSER"} = $dbuser;
local $ENV{"JSOC_DBNAME"} = $dbname;

$logcontent = $logcontent . "which: " . `which dlMOCDataFiles\.pl`;
$logcontent = $logcontent . "which: " . `which sftpScript\.exp`;

# Create specification file
my(@timearr);
my(@threemosago);
@timearr = localtime(); # today
@threemosago = localtime(time - 90 * 86400); # 90 days ago

$yr = $timearr[5] + 1900;
$dofy = $timearr[7]; 

$yr90 = $threemosago[5] + 1900;
$dofy90 = $threemosago[7];

if ($yr90 < $yr)
{
    $pbegin = sprintf("%03d", $dofy90);
    $pend = 366;
    $pyr = $yr - 1;
    $cbegin = 1;
    $cend = $dofy;
}
else
{
    $cbegin = sprintf("%03d", $dofy90);
    $cend = sprintf("%03d", $dofy);
}

$pspec = "";
$cspec = "";

if (defined($pbegin))
{
    $pspec = "${kSPECPREFIX}lzp/${pyr}_|s[$pbegin-$pend]::000[1-9]_${pyr}_$kFILESUFFIX\n${kSPECPREFIX}lzp/${pyr}_|s[$pbegin-$pend]::00[1-5][0-9]_${pyr}_$kFILESUFFIX\n${kSPECPREFIX}lzp/${pyr}_|s[$pbegin-$pend]::006[0-3]_${pyr}_$kFILESUFFIX\n${kSPECPREFIX}lzp/${pyr}_|s[$pbegin-$pend]::0129_${pyr}_$kFILESUFFIX\n";
}

if (defined($cbegin))
{
    $cspec = "${kSPECPREFIX}lzp/${yr}_|s[$cbegin-$cend]::000[1-9]_${yr}_$kFILESUFFIX\n${kSPECPREFIX}lzp/${yr}_|s[$cbegin-$cend]::00[1-5][0-9]_${yr}_$kFILESUFFIX\n${kSPECPREFIX}lzp/${yr}_|s[$cbegin-$cend]::006[0-3]_${yr}_$kFILESUFFIX\n${kSPECPREFIX}lzp/${yr}_|s[$cbegin-$cend]::0129_${yr}_$kFILESUFFIX\n";
}

$fullspec = "$pspec$cspec";

open(SPECFILE, ">$permdataPath/mocDlLzpSpec.txt");
print SPECFILE "root:moc\n\n";
print SPECFILE $fullspec;
close(SPECFILE);

# Dump log
open(LOGFILE, ">$logfile") || die "Couldn't write to logfile '$logfile'\n";
print LOGFILE $logcontent;
close(LOGFILE);

# Download LZP files to $dataPath
$cmd = "dlMOCDataFiles\.pl -c $permdataPath/mocDlLzpSpec.txt -s $permdataPath/mocDlLzpStatus.txt -r $dataPath -t 120 $force";

system("$cmd 1>>$logfile 2>&1");

# Notify people who care if an error has occurred
$cmd = "fdsNotification\.pl -l $logfile -n $notlist";
system($cmd);

sub PrintUsage
{
    print "\tdllzp.pl {-b <branch> | -c <config file>} [-f]\n"
}
