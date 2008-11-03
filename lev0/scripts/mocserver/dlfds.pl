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
#
# To set up ssh-agent, from the machine on which the cron job runs, type
# 'nohup ssh-agent -c > /home/jsoc/.ssh-agent' then
# 'ssh-add /home/jsoc/.ssh/id_rsa'
#
# There -b flag is a short-cut to specify the configuration file.
#  -b dev -   use the development configuration file.  This configuration file
#             is intended to be used by running the script as user jsoc on 
#             machine maelstrom.
#  -b live -  use the 'live' (both pre- and post- launche) configuration file.  
#             This configuration file is intended to be used by running 
#             the script as user jsoc on machine j0.

use FindBin qw($Bin);

my($kDEV) = "dev";
my($kLIVE) = "live";
my($kCFGPREFIX) = "dlfds_conf";
my($kSSHAGENTCONF) = "ssh-agent_conf";
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

my($err) = 0;

my($jsocmach);

my($arg);
my($argc) = scalar(@ARGV);
my($force) = "";
my($cfgfile);
my($line);

my($sshagentConf);
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
    printf STDERR "Cannot read configuration file.";
    exit(1);
}

open(CFGFILE, "< $cfgfile");

while($line = <CFGFILE>)
{
    chomp($line);

    if ($line =~ /$kSSHAGENTCONF\s*=\s*(\S+)\s*/)
    {
        $sshagentConf = $1;
    }
    elsif ($line =~ /$kLOGFILE\s*=\s*(\S+)\s*/)
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
if (-e $sshagentConf)
{
    my($agentpid);
    my($agentcmd);
    my($agentfound) = 0;
    
    $logcontent = $logcontent . "found ssh-agent environment variables.\n";

    # 'local' ensures that calls outside of this block don't see these environmental changes.
    # It creates a copy of the environment, which then gets restored upon exit of the block.
    local %ENV = %ENV;

    # parse the file

    open(SSHCONF, "<$sshagentConf") || die "Couldn't open ssh-agent configuration file '$sshagentConf'.\n";
    while($line = <SSHCONF>)
    {
        chomp($line);
        if ($line =~ /^setenv\s+(.+)\s+(.+);/)
        {
            $logcontent = $logcontent . "\tgot ssh env var $1 = $2\n";
            $ENV{$1} = $2;
        }
    }

    # ensure that the environment variables point to a running ssh-agent
    $agentpid = $ENV{'SSH_AGENT_PID'};
    
    if (defined($agentpid))
    {
        # $agentlink = readlink("/proc/$agentpid/exe");
        #
        # can't do what you'd like to do (read the link and see that it points to
        # the ssh-agent process file) - only root can do that.
        # Not sure why.  For other processes, you can read the links in /proc, and
        # for links to ssh-agent outside of /proc, you can read the links.
        # You can't even do if (-e "/proc/$agentpid/exe").
        #
        # You can look at /proc/$agentpid/cmdline though.
        
        $agentcmd = `cat "/proc/$agentpid/cmdline"`;
        
        if (defined($agentcmd))
        {
            if ($agentcmd =~ /ssh-agent/)
            {
                $agentfound = 1;
            }
        }
    }
    
    if (!$agentfound)
    {
        $logcontent = $logcontent . "LOGALL: ssh-agent NOT RUNNING.\n";
        $err = 1;
    }
    else
    {
        $logcontent = $logcontent . "ssh-agent running (pid $agentpid).\n";
    }


    # Dump log
    open(LOGFILE, ">$logfile") || die "Couldn't write to logfile '$logfile'\n";
    print LOGFILE $logcontent;
    close(LOGFILE);

    if (!$err)
    {
        system("$cmd 1>>$logfile 2>&1");
    }
}
else
{
    $logcontent = $logcontent . "couldn't find ssh-agent (no $sshagentConf).\n";

    # Dump log
    open(LOGFILE, ">$logfile") || die "Couldn't write to logfile '$logfile'\n";
    print LOGFILE $logcontent;
    close(LOGFILE);

    system("$cmd 1>>$logfile 2>&1");
}

if (!$err)
{
    # Ingest FDS files into $fdsseriesname
    $cmd = "fdsIngest\.pl $dataPath/fds -s $namespace\.$fdsseriesname -r";
    system("$cmd 1>>$logfile 2>&1");

    # Ingest orbit vectors into $orbvSeries
    $cmd = "extract_fds_statev ns=$namespace seriesIn=$fdsseriesname seriesOut=$orbvseriesname";

    system("$cmd 1>>$logfile 2>&1");
}

# Notify people who care if an error has occurred
$cmd = "fdsNotification\.pl -l $logfile -s \"Error downloading/ingesting MOC FDS data products\" -n $notlist";
system($cmd);

sub PrintUsage
{
    print "\tdlfds.pl {-b <branch> | -c <config file>} [-f]\n"
}
