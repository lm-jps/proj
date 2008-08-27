#!/usr/bin/perl -w 

# This should be run in the background.  It will periodically examine the input series
# and if it finds new data, then it will create a new EOF synoptic chart.  Run
# like this: makeeofsynop.pl -b live & from user jsoc on any machine that can
# run both JSOC and SOI modules and access DSDS.  This script will figure out
#  the machine types
# and run executables of the correct architecture.

use FindBin qw($Bin);

# constants
my($kDEV) = "dev";
my($kLIVE) = "live";
my($kCFGPREFIX) = "eofsynop_conf";
my($kSSHAGENTCONF) = "ssh-agent_conf";
my($kLOGFILE) = "log";
my($kSCRIPTPATH) = "scrpath";
my($kBINPATH) = "binpath";
my($kSOIBINPATH) = "soibinpath";
my($kPERMDATAPATH) = "permdpath";
my($kDATAPATH) = "dpath";
my($kNAMESPACE) = "ns";
my($kSYNOPDIR) = "synopdir";
my($kDBUSER) = "dbuser";
my($kDBNAME) = "dbname";
my($kNOTLIST) = "notify";

# vars
my($arg);
my($branch);
my($cfgfile);

my($sshagentConf);
my($logfile);
my($scriptPath);
my($binPath);
my($soiBinPath);
my($permdataPath);
my($dataPath);
my($namespace);
my($dbuser);
my($dbname);
my($notlist);

my($jsocmach);
my($soimach);
my($jsocPath);
my($soiPath);

my($dstart);
my($dend);
my($indata);

my($cmd);

while ($arg = shift(@ARGV))
{
    if ($arg eq "-b")
    {
        $branch = shift(@ARGV);

        if ($branch ne $kDEV && $branch ne $kLIVE)
        {	
            PrintUsage();
            exit(1);
        }
    }
    elsif ($arg eq "-c")
    {
        $cfgfile = shift(@ARGV);
    }
    else
    {
        PrintUsage();
        exit(1);
    }
}

if (!defined($branch) && !defined($cfgfile))
{
    PrintUsage();
    exit(1);
}

# read configuration file
if (!defined($cfgfile) && defined($branch))
{
    if ($branch eq $kLIVE)
    {
        $branch = "";
    }
    else
    {
        $branch = "_" . $branch;
    }

    # Use defalt configuration file
    $cfgfile = "$Bin/${kCFGPREFIX}${branch}.txt";
}

if (!(-f $cfgfile))
{
    printf STDERR "Cannot read configuration file '$cfgfile'.";
    exit(1);
}

open(CFGFILE, "< $cfgfile");

while($line = <CFGFILE>)
{
    chomp($line);

    if ($line =~ /^\s*$kSSHAGENTCONF\s*=\s*(\S+)\s*/)
    {
        $sshagentConf = $1;
    }
    elsif ($line =~ /^\s*$kLOGFILE\s*=\s*(\S+)\s*/)
    {
	$logfile = $1;
    }
    elsif ($line =~ /^\s*$kSCRIPTPATH\s*=\s*(\S+)\s*/)
    {
	$scriptPath = $1;
    }
    elsif ($line =~ /^\s*$kBINPATH\s*=\s*(\S+)\s*/)
    {
	$binPath = $1;
    }
    elsif ($line =~ /^\s*$kSOIBINPATH\s*=\s*(\S+)\s*/)
    {
	$soiBinPath = $1;
    }
    elsif ($line =~ /^\s*$kPERMDATAPATH\s*=\s*(\S+)\s*/)
    {
	$permdataPath = $1;
    }
    elsif ($line =~ /^\s*$kDATAPATH\s*=\s*(\S+)\s*/)
    {
	$dataPath = $1;
    }
    elsif ($line =~ /^\s*$kNAMESPACE\s*=\s*(\S+)\s*/)
    {
	$namespace = $1;
    }
    elsif ($line =~ /^\s*$kSYNOPDIR\s*=\s*(\S+)\s*/)
    {
        $synopdir = $1;
    }
    elsif ($line =~ /^\s*$kDBUSER\s*=\s*(\S+)\s*/)
    {
	$dbuser = $1;
    }
    elsif ($line =~ /^\s*$kDBNAME\s*=\s*(\S+)\s*/)
    {
	$dbname = $1;
    }
    elsif ($line =~ /^\s*$kNOTLIST\s*=\s*(.+)/)
    {
	$notlist = $1;
    }
}

close(CFGFILE);

# set environment so that script can find other scripts/exes 
$ENV{"JSOCROOT"} = "$Bin/../../../..";

if (!defined($ENV{"JSOC_MACHINE"}))
{
    $jsocmach = `$ENV{"JSOCROOT"}/build/jsoc_machine.csh`;
    chomp($jsocmach); # back-ticks cause a newline to be appended.
}
else
{
    $jsocmach = $ENV{"JSOC_MACHINE"};
}

if (!defined($ENV{"MACHINE"}))
{
    $soimach =  `$ENV{"JSOCROOT"}/base/local/scripts/soi_machine.csh`;
    chomp($soimach);
}
else
{
    $soimach = $ENV{"MACHINE"};
}

# For JSOC modules and scripts, use the paths specified in the config file
$jsocPath = "$scriptPath:$binPath/$jsocmach:$ENV{\"PATH\"}";
$soiPath = "$soiBinPath/_$soimach:$ENV{\"PATH\"}";

$ENV{"JSOC_DBUSER"} = $dbuser;
$ENV{"JSOC_DBNAME"} = $dbname;

# Specify the input data
$dstart = $dlast + 1;
$dend = $dlast + 10;
$indata = "prog:mdi_eof,level:lev1.8,series:fd_M_96m_01d\[$dstart\]";

{
    local $ENV{"PATH"} = $soiPath;

    $cmd = "$fdradialexe RADIALCORR=1 in=$indata out=$SYNOPDIR/fd_M_radial_96m_01d/fd_M_radial_96m_01d.00$day";
    system($cmd);

    $cmd = "v2helio z=1 MOFFSET=1 VCORLEV=0 MAPRMAX=0.994 MAPLGMAX=89.75 MAPLGMIN=-89.75 MAPBMAX=90 LGSHIFT=2 MAPMMAX=1800 SINBDIVS=540 in=$SYNOPDIR/fd_M_radial_96m_01d/fd_M_radial_96m_01d.00$day/ out=$SYNOPDIR/fd_Mag_remap_01d/fd_Mag_remap_01d.00$day/";

"CR=$cr checkqual=1 qualmask=0x402c01f2 center=$center nsig=3.0 -l in=$SYNOPDIR/fd_M_radial_96m_01d/fd_M_radial_96m_01d.00$day out=$SYNOPDIR/fd_Mag_remap_01d/fd_Mag_remap_01d.00$day ";
    system($cmd);

    # Now run average images to make EOF synoptic chart
    $cmd = "jsoc_eofsynop";
}

#print "jsocPath $jsocPath\n";
#print "soiPath $soiPath\n";



sub PrintUsage
{
    print "\tmakeeofsynop.pl {-b <branch> | -c <config file>}\n"
}
