#!/usr/bin/perl -w 

# This script (and sub-scripts) uses ssh to connect to the MOC Product Server.  Unless
# certain actions are taken, sshd is going to ask for a password and/or key passphrase.
# If sshd asks for a password or passphrase, the script will fail.  There are a few
# options for preventing the authentication challenge from being issued (modifying the 
# authorized_keys file on the MOC Product Server for one).  However, when run by 
# production-type users, like 'jsoc' and 'production', care must be taken to not place
# a passphrase-less pub key into the server's authorized_keys file.  Instead
# the user should add a private key (that requires a passphrase) to the list of keys
# known by ssh-agent.  Then the user may run this script.

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
my($kNAMESPACE) = "ns";
my($kFDSSERIESNAME) = "fdsseries";
my($kORBVSERIESNAME) = "orbvseries";
my($kDBUSER) = "dbuser";
my($kDBNAME) = "dbname";
my($kNOTLIST) = "notify";

my($kSSHCONFFILE) = "$ENV{'HOME'}/.ssh-agent";
my($kSSH_AUTH_SOCK) = "SSH_AUTH_SOCK";
my($kSSH_AGENT_PID) = "SSH_AGENT_PID";

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
my($fdsseriesname);
my($orbvseriesname);
my($dbuser);
my($dbname);
my($notlist);

my($cmd);

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
    elsif ($line =~ /$kFDSSERIESNAME\s*=\s*(\S+)\s*/)
    {
	$fdsseriesname = $1;
    }
    elsif ($line =~ /$kORBVSERIESNAME\s*=\s*(\S+)\s*/)
    {
	$orbvseriesname = $1;
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
$ENV{"PATH"} = "$scriptPath:$binPath/$jsocmach:$ENV{\"PATH\"}";
$ENV{"JSOC_DBUSER"} = $dbuser;
$ENV{"JSOC_DBNAME"} = $dbname;

$logcontent = $logcontent . "which: " . `which dlMOCDataFiles\.pl`;
$logcontent = $logcontent . "which: " . `which sftpScript\.exp`;
$logcontent = $logcontent . "which: " . `which fdsIngest\.pl`;
$logcontent = $logcontent . "which: " . `which extract_fds_statev`;

# Download FDS files to $dataPath
$cmd = "dlMOCDataFiles\.pl -c $scriptPath/mocDlFdsSpec.txt -s $permdataPath/mocDlFdsStatus.txt -r $dataPath -t 30 $force";

# Set the environment so that ssh-agent can be accessed by sftp and scp calls
# (can't factor out the system call - it needs to be in the same scope as the 'local' statements)
if (-e $kSSHCONFFILE)
{
    $logcontent = $logcontent . "found ssh-agent.\n";

    # 'local' ensures that calls outside of this block don't see these environmental changes.
    # It creates a copy of the environment, which then gets restored upon exit of the block.
    local %ENV = %ENV;

    # parse the file

    open(SSHCONF, "<$kSSHCONFFILE") || die "Couldn't open ssh-agent configuration file '$kSSHCONFFILE'.\n";
    while($line = <SSHCONF>)
    {
        chomp($line);
        if ($line =~ /^setenv\s+(.+)\s+(.+);/)
        {
            $logcontent = $logcontent . "got ssh env var $1 = $2\n";
            $ENV{$1} = $2;
        }
    }

    # Dump log
    open(LOGFILE, ">$logfile") || die "Couldn't write to logfile '$logfile'\n";
    print LOGFILE $logcontent;
    close(LOGFILE);

    system("$cmd 1>>$logfile 2>&1");
}
else
{
    $logcontent = $logcontent . "couldn't find ssh-agent.\n";

    # Dump log
    open(LOGFILE, ">$logfile") || die "Couldn't write to logfile '$logfile'\n";
    print LOGFILE $logcontent;
    close(LOGFILE);

    system("$cmd 1>>$logfile 2>&1");
}

# Ingest FDS files into $fdsseriesname
$cmd = "fdsIngest\.pl $dataPath/fds -s $namespace\.$fdsseriesname -r";
system("$cmd 1>>$logfile 2>&1");

# Ingest orbit vectors into $orbvSeries
$cmd = "extract_fds_statev ns=$namespace seriesIn=$fdsseriesname seriesOut=$orbvseriesname";
system("$cmd 1>>$logfile 2>&1");

# Notify people who care if an error has occurred
$cmd = "fdsNotification\.pl -l $logfile -n $notlist";
system($cmd);

sub PrintUsage
{
    print "\tdlfds.pl {-b <branch> | -c <config file>} [-f]\n"
}
