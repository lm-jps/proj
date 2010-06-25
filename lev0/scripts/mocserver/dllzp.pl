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
use Time::Local;

my($kDEV) = "dev";
my($kLIVE) = "live";
my($kCFGPREFIX) = "dllzp_conf";
my($kSSHAGENTCONF) = "ssh-agent_conf";
my($kLOGFILE) = "log";
my($kSCRIPTPATH) = "scrpath";
my($kBINPATH) = "binpath";
my($kPERMDATAPATH) = "permdpath";
my($kDATAPATH) = "dpath";
my($kNAMESPACE) = "ns";
my($kDBUSER) = "dbuser";
my($kDBNAME) = "dbname";
my($kNOTLIST) = "notify";
my($kTIMEWINB) = "timewinb";
my($kTIMEWINE) = "timewine";

my($kSPECPREFIX) = "filespec:";
my($kFILESUFFIX) = "[0-9][0-9][0-9]_[0-9][0-9]\\.hkt\\S*";

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
my($dbuser);
my($dbname);
my($notlist);
my($timewinb);
my($timewine);

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
my($msg);    # holds messages to print to log file

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
    print STDERR "Cannot read configuration file.";
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
        if (-e $logfile)
        {
           unlink $logfile;
        }
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
    elsif ($line =~ /$kTIMEWINB\s*=\s*(.+)/)
    {
       $timewinb = $1;
    }
    elsif ($line =~ /$kTIMEWINE\s*=\s*(.+)/)
    {
       $timewine = $1;
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

# Create specification file
my(@timearr);
my(@threemosago);
@timearr = localtime(); # today
@threemosago = localtime(time - 90 * 86400); # 90 days ago

# Read file containing date (GMT mday, mon, year) of files last downloaded
my($mday);   # 1..31
my($mon);    # 0..11
my($year);   # year - 1900
my($lastdl); # seconds since epoch


if (open(LASTDATE, "<$permdataPath/datefile.txt"))
{
   $mday = <LASTDATE>;
   chomp($mday);
   $mon = <LASTDATE>;
   chomp($mon);
   $year = <LASTDATE>;
   chomp($year);
   close(LASTDATE);
}
else
{
   # No date file - this is an error. You can't make a dummy with a very old date because
   # this will cause the check for missing files at the end of this script to fail if
   # in fact there are no files to download (they've already been downloaded). So before the
   # first time you run this script, create a file with a date that is 1 day before the next
   # set of files you expect to download.
   $msg = "Missing required date file '$permdataPath/datefile.txt'.\n";
   DumpLog($logfile, $msg);
   exit(1);
}

# The date of the files last downloaded
$lastdl = timegm(0, 0, 0, $mday, $mon, $year);

# Don't run if today's files have already been downloaded.
# The new files for day $lastdl + 1d arrive between $timewinb and $timewine hours 
# after $lastdl + 1d. The 1 day and the $timewinb are approximate, and
# could be off by a leapsecond, but we don't need to worry about that 
# kind of precision.

if (!$force)
{
   # Don't check for files not being ready if the user is calling with force. The user
   # might be (and most likely is) trying to fetch older files.
   # ARTXXXXXXXX
   #if (time() + 100000 - $lastdl < ($timewinb + 24) * 60 * 60)
   if (time() - $lastdl < ($timewinb + 24) * 60 * 60)
   {
      # Too soon to run
      $msg = "Not time to check for new files yet, exiting.\n";
      DumpLog($logfile, $msg);
      exit(0);
   }
}

$yr = $timearr[5] + 1900;
$dofy = $timearr[7] + 1; 

$yr90 = $threemosago[5] + 1900;
$dofy90 = $threemosago[7] + 1;

if ($yr90 < $yr)
{
    $pbegin = sprintf("%03d", $dofy90);
    $pend = 366;
    $pyr = $yr - 1;
    $cbegin = "001";
    $cend = sprintf("%03d", $dofy);
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

# Download LZP files to $dataPath
$cmd = "dlMOCDataFiles\.pl -c $permdataPath/mocDlLzpSpec.txt -s $permdataPath/mocDlLzpStatus.txt -r $dataPath -t 120 $force";

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
    DumpLog($logfile, $logcontent);

    if (!$err)
    {
        system("$cmd 1>>$logfile 2>&1");
    }
}
else
{
    $logcontent = $logcontent . "couldn't find ssh-agent (no $sshagentConf).\n";

    # Dump log
    DumpLog($logfile, $logcontent);

    system("$cmd 1>>$logfile 2>&1");
}

if (!$err)
{
   my($filesdownloaded) = 0;
   my($dlyear);
   my($dlday);
   my($lyear);
   my($lday);
   my($filesmissing) = 0;

   # Make these group writeable so that production can delete them 
   $cmd = "chmod -R g+w $dataPath";
   system("$cmd 1>>$logfile 2>&1");

   # Check for any errors

   # Get date of most recent file downloaded (if any were downloaded at all)

   # See if any files were downloaded - dlMOCDataFiles.pl will leave some bread crumbs in the log.
   open(LOGFILE, "<$logfile") || die "Couldn't open logfile '$logfile' for reading.\n";

   $dlyear = -1;
   $dlday = -1;

   while (defined($line = <LOGFILE>))
   {
      chomp($line);
      if ($line =~ /Executing\sscp.+_(\d\d\d\d)_(\d\d\d)/)
      {
         $filesdownloaded = 1;
         $lyear = $1;
         $lday = $2;

         if ($lyear > $dlyear || ($lyear == $dlyear && $lday > $dlday))
         {
            $dlyear = $lyear;
            $dlday = $lday;
         }
      }
   }

   close(LOGFILE);

   if ($filesdownloaded)
   {
      # Convert to DD\nMM\nYYYY - I don't know how to easily do this, except by using
      # time_convert (which will properly handle leap days/secs).
      my($tcCmdLine) = "time_convert ord=${dlyear}\.${dlday}_UTC o=cal zone=UTC |";
      my($timestr);
      my($newlastdl);

      if (!open(TIMECONV, $tcCmdLine))
      {
         $msg = "LOGALL: Couldn't run time_conv: $tcCmdLine\n";
         DumpLog($logfile, $msg);
         exit(1);
      }
	
      if (!defined($timestr = <TIMECONV>))
      {
         $msg = "LOGALL: Problem running time_conv: $tcCmdLine\n";
         DumpLog($logfile, $msg);
         close TIMECONV;
         exit(1);
      }

      if ($timestr =~ /\d\d\d\d\.(\d\d)\.(\d\d)/)
      {
         $mday = $2;
         $mon = $1 - 1;
         $year = $dlyear - 1900;
      }
      else
      {
         $msg = "LOGALL: Invalid time string returned by '$tcCmdLine'.\n";
         DumpLog($logfile, $msg);
         exit(1);
      }

      # Update the datefile (if the date of files downloaded is newer than the date of the files last
      # downloaded).
      $newlastdl = timegm(0, 0, 0, $mday, $mon, $year);
      if ($newlastdl > $lastdl)
      {
         open(LASTDATE, ">$permdataPath/datefile.txt") || die "Couldn't open datefile '$permdataPath/datefile.txt' for writing.\n";
         print LASTDATE "$mday\n";
         print LASTDATE "$mon\n";
         print LASTDATE "$year\n";
         close(LASTDATE);

         # Update the date of last downloaded files.
         $lastdl = $newlastdl;
      }

      # Call Carl's ingestion script
      $msg = "Calling dsdf.pl (ingestion script):\n";
      DumpLog($logfile, $msg);
      $cmd = "/usr/bin/perl $scriptPath/../proj/lev0/scripts/hk/dsdf.pl moc";
      system("$cmd 1>>$logfile 2>&1");
   }

   # Check to see if we got all files we should have gotten (should never be more than 60
   # hours after we have downloaded the last set of files)

   # The files for day X are guaranteed to be 
   # present by 12 UTC on day X + 1.  If it is past 12 UTC on day X + 1 (this script
   # should be run at 12:15 UTC on day X + 1 to check for this)
   # and the files for day X have not been downloaded, then print a failure message.
   if (time() - $lastdl > ($timewine + 24) * 60 * 60)
   {
      $filesmissing = 1;
   }
   
   if ($filesmissing)
   {
      # Get ordinal date for for first day after date of last download of files.
      my($realyr) = $year + 1900;
      my($realmo) = $mon + 1;
      my($realdy) = $mday + 1; # Add 1 day to $lastdl (ok if this is 32)
      my($tcCmdLine) = "time_convert time=${realyr}\.${realmo}\.${realdy}_UTC o=ord zone=UTC |";
      my($timestr);
      my($nextday);

      if (!open(TIMECONV, $tcCmdLine))
      {
         $msg = "LOGALL: Couldn't run time_conv: $tcCmdLine\n";
         DumpLog($logfile, $msg);
         exit(1);
      }
	
      if (!defined($timestr = <TIMECONV>))
      {
         $msg = "LOGALL: Problem running time_conv: $tcCmdLine\n";
         DumpLog($logfile, $msg);
         close TIMECONV;
         exit(1);
      }

      if ($timestr =~ /(\d\d\d\d)\.(\d\d\d)/)
      {
         $nextday = "$1_$2";
      }
      else
      {
         $msg = "LOGALL: Invalid time string returned by '$tcCmdLine'.\n";
         DumpLog($logfile, $msg);
         exit(1);
      }

      $msg = "LOGALL: ***FAILURE TO DOWNLOAD FILES FOR DAY ${nextday}***\n";
      DumpLog($logfile, $msg);
      exit(1);
   }
}

# Notify people who care if an error has occurred
$cmd = "fdsNotification\.pl -l $logfile -s \"Error downloading MOC LZP data products\" -n $notlist";
system($cmd);

sub PrintUsage
{
    print "\tdllzp.pl {-b <branch> | -c <config file>} [-f]\n"
}

sub DumpLog
{
   my($logfile) = $_[0];
   my($content) = \$_[1];

   if (open(LOGFILE, ">>$logfile"))
   {
      print LOGFILE $$content;
      $$content = "";
      close(LOGFILE);
   }
   else
   {
      print STDERR "Couldn't write to logfile '$logfile'\n";
      exit(1);
   }
}
