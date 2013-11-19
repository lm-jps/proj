#!/usr/bin/perl 

# Here's how to run this script (from scratch)
#  ssh jsoc@j0
#  cd /home/jsoc/exports
#  rm keep_running
#  exportmanage.pl -jsocdev procser=jsoc.export_procs &
# To start the web version for jsoc.stanford.edu repeat the above but do:
#  rm keep_running_web
#  exportmanage.pl -jsocweb procser=jsoc.export_procs&
# At some point, I'll add the production version of the script, which
#  would then be run as "exportmanage.pl -jsocpro &"
# To start the debug test version for jsoc2 do:
#  rm keep_running_test
#  exportmanage.pl -jsoctest procser=jsoc.export_procs &

# To test the entire export workflow:
#  1. Run export manage like so:
#     su <USER>
#        where <USER> is the user whose $PATH contains the paths to the modules/scripts that have new code to be tested.
#        The manager program makes shell scripts that call programs and scripts by calling the program/script base file
#        name (not the full path to the program/script). So the actual program/script that runs will be the one
#        that the shell resolves using the $PATH variable. For example, if the user arta makes
#        changes to the jsoc_export_as_fits in arta's home directory, and arta's $PATH contains a pointer to the binaries
#        in arta's home directory, then arta should run 'su arta' before continuing. IMPORTANT - <USER> will need 
#        permissions to write into /home/jsoc/exports/tmp and /home/jsoc/exports/logs. Generally <USER> will be a 
#        member of the group jsoc, so this requirement will be met.
#     cd /home/jsoc/exports
#     /home/jsoc/cvs/Development/JSOC/proj/util/scripts/exportmanage.pl -root <ROOT> -dbuser <DBUSER> -dbhost <DBHOST> -manager <MANAGER> -runflag <RFLAG> &
#        where <ROOT> is the CVS code tree root containing <MANAGER>
#        and <DBUSER> is the PG user who the manager connects as (defaults to "production"). IMPORTANT - the manager will
#           write records to tables that require elevated permissions. Most likely, you'll need to connect to the database 
#           as user production to do this. This means that you'll need to place the password for user production in 
#           your .pgpass file.
#        and <DBHOST> is the host of the database server (defaults to "hmidb"). For internal exports, should be "hmidb", for exports from public db, should be "hmidb2"
#        and <MANAGER> is the name of the manager program (defaults to "jsoc_export_manage"). IMPORTANT - you should include 
#           a "-t" flag. This will cause the manager program to run in test mode, which means that it will process records 
#           in jsoc.export_new that contain the special test status of 12 (instead of the regular status of 2).
#        and <RFLAG> is the file flag that keeps this script running in a loop (defaults to keep_running in cdir)
#
#        Example : /home/jsoc/cvs/Development/JSOC/proj/util/scripts/exportmanage.pl -root /home/arta/cvs/JSOC -dbuser production -dbhost hmidb2 -manager "jsoc_export_manage -t procser=jsoc.export_procs" -runflag keepruntest.txt -logflag Test &
#
#  2. Point a browser at http://jsoc.stanford.edu/ajax/exportdatatest.html and export something.

use strict;
use warnings;

use FileHandle;
use Fcntl ':flock';
use FindBin qw($RealBin);
use lib "$RealBin/../../../localization";
use drmsparams;

use constant kExportDir => "/home/jsoc/exports";
# use constant kMailList => "arta\@sun.stanford.edu jeneen\@sun.stanford.edu phil\@sun.stanford.edu";
use constant kMailList => "arta\@sun.stanford.edu";
use constant kMailMessage1 => "exportmanage.pl could not start jsoc_export_manage. This is a critical failure.\nYou should probably contact Art, who was also notified and should respond shortly.\n";
use constant kMailMessage2 => "jsoc_export_manage died in response to an unhandled signal (e.g., a segfault).\n";
use constant kMailMessage3 => "Could not open log export-daemon log file for writing.\nThis is not a critical failure, but we should\nfix this so that we can track-down future export problems more easily.\nContact Art.\n";

use constant kMsgType1 => "msgtype1";
use constant kMsgType2 => "msgtype2";
use constant kMsgType3 => "msgtype3";

use constant kMsgQInterval => 600; # 10 minutes, at least, between message inserts
use constant kMsgQSendInterval => 120; # 2 minutes between mailing of messages

use constant kLogFlagInt => "Int";
use constant kLogFlagExt => "Ext";

my($config) = new drmsparams;

my($kINTERNALFLAG) = "/home/jsoc/exports/keep_running";
my($kWEBFLAG) = "/home/jsoc/exports/keep_running_web";
my($kTESTFLAG) = "/home/jsoc/exports/keep_running_test";
#
my($kJSOCDEV_ROOT) = "/home/jsoc/cvs/Development/JSOC";
my($kJSOCDEV_DBUSER) = "production";
my($kJSOCDEV_DBNAME) = $config->get(DBNAME);
my($kJSOCDEV_DBHOST) = $config->get(SERVER);
my($kJSOCDEV_MANAGE) = "jsoc_export_manage";
use constant kProcInfoSeriesDev => "jsoc.export_procs";
#
my($kJSOCPRO_ROOT) = "/home/jsoc/cvs/JSOC";
my($kJSOCPRO_DBUSER) = "production";
my($kJSOCPRO_DBNAME) = $config->get(DBNAME);
my($kJSOCPRO_DBHOST) = $config->get(SERVER);
my($kJSOCPRO_MANAGE) = "jsoc_export_manage";
use constant kProcInfoSeriesPro => "jsoc.export_procs";
#
my($kJSOCWEB_ROOT) = "/home/jsoc/cvs/Development/JSOC";
my($kJSOCWEB_DBUSER) = "production";
my($kJSOCWEB_DBNAME) = "jsoc";
my($kJSOCWEB_DBHOST) = "hmidb2";
my($kJSOCWEB_MANAGE) = "jsoc_export_manage";
use constant kProcInfoSeriesWeb => "jsoc.export_procs";
#
my($kJSOCTEST_ROOT) = "/home/jsoc/cvs/Development/JSOC";
my($kJSOCTEST_DBUSER) = "phil";
my($kJSOCTEST_DBNAME) = $config->get(DBNAME);
my($kJSOCTEST_DBHOST) = $config->get(SERVER);
my($kJSOCTEST_MANAGE) = "jsoc_export_manage_test";
use constant kProcInfoSeriesTst => "jsoc.export_procs";

my($runningflag) = $kINTERNALFLAG;
my($arg);
my($root);
my($dbhost) = $config->get(SERVER);
my($dbname) = $config->get(DBNAME);
my($dbuser) = "production";
my($binpath);
my($manage) = "jsoc_export_manage";
my($logfile);
my($daemonlog);
my($lckfh);
my($msg);
my($logflag);
my($procser);

while ($arg = shift(@ARGV))
{
    if ($arg eq "-root")
    {
        $root = shift(@ARGV);
        $binpath = "$root/bin";
    }
    elsif ($arg eq "-dbhost")
    {
        $dbhost = shift(@ARGV);
    }
    elsif ($arg eq "-dbuser")
    {
        $dbuser = shift(@ARGV);
    }
    elsif ($arg eq "-dbname")
    {
        $dbname = shift(@ARGV);
    }
    elsif ($arg eq "-manager")
    {
        $manage = shift(@ARGV);
    }
    elsif ($arg eq "-runflag")
    {
        $runningflag = shift(@ARGV);
    }
    elsif ($arg eq "-logflag")
    {
        $logflag = shift(@ARGV);
    }
    elsif ($arg eq "-procser")
    {
        $procser = shift(@ARGV);
    }
    elsif ($arg eq "-jsocdev")
    {
        $root = $kJSOCDEV_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCDEV_DBUSER;
        $dbname = $kJSOCDEV_DBNAME;
        $dbhost = $kJSOCDEV_DBHOST;
        $manage = $kJSOCDEV_MANAGE;
        $runningflag = $kINTERNALFLAG;
        $logflag = kLogFlagInt;
        $procser = &kProcInfoSeriesDev;
    }
    elsif ($arg eq "-jsocpro")
    {
        $root = $kJSOCPRO_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCPRO_DBUSER;
        $dbname = $kJSOCPRO_DBNAME;
        $dbhost = $kJSOCPRO_DBHOST;
        $manage = $kJSOCPRO_MANAGE;
        $runningflag = $kINTERNALFLAG;
        $logflag = kLogFlagInt;
        $procser = &kProcInfoSeriesPro;
    }
    elsif ($arg eq "-jsocweb")
    {
        $root = $kJSOCWEB_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCWEB_DBUSER;
        $dbname = $kJSOCWEB_DBNAME;
        $dbhost = $kJSOCWEB_DBHOST;
        $manage = $kJSOCWEB_MANAGE;
        $runningflag = $kWEBFLAG;
        $logflag = kLogFlagExt;
        $procser = &kProcInfoSeriesWeb;
    }
    elsif ($arg eq "-jsoctest")
    {
        $root = $kJSOCTEST_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCTEST_DBUSER;
        $dbname = $kJSOCTEST_DBNAME;
        $dbhost = $kJSOCTEST_DBHOST;
        $manage = $kJSOCTEST_MANAGE;
        $runningflag = $kTESTFLAG;
        $logflag = "Test";
        $procser = &kProcInfoSeriesTst;
    }
}

# Only run on j0.Stanford.EDU
# if ($ENV{HOSTNAME} ne "j0.Stanford.EDU") {
#    die "I will only run on j0.Stanford.EDU\n";
#}

# Don't run if somebody is already managing the export
$lckfh = FileHandle->new(">$runningflag.lck");
unless (flock($lckfh, LOCK_EX|LOCK_NB)) 
{
   print "$0 is already running. Exiting.\n";
   exit(1);
}

#if (-e $runningflag)
#{
#    die "Can't manage export; another process is already managing it.\n";
#}

if (defined($binpath))
{
    $binpath = "$binpath/$ENV{\"JSOC_MACHINE\"}";
}
else
{
    $binpath = "";
}

#local $ENV{"PATH"} = "$binpath:$ENV{\"PATH\"}";
#local $ENV{"PATH"} = "$scrpath:$ENV{\"PATH\"}";
local $ENV{"JSOCROOT"} = $root;
local $ENV{"JSOC_DBUSER"} = $dbuser;
local $ENV{"JSOC_DBNAME"} = $dbname;

#`touch $runningflag`;
`echo $$ > $runningflag`;
$logfile = kExportDir . "/logs/" . `date +"%F_%R.log"`;
open(LOG, ">> $logfile") || die "Couldn't open logfile '$logfile'.\n";
$daemonlog = kExportDir . "/logs/exportlog${logflag}.txt";

my($datenow) = `date`;
chomp($datenow);
$msg = "Started by $ENV{'USER'} at $datenow on machine $ENV{'HOST'} using $dbhost.\n";
print LOG $msg;

my($rout);
my($cmd);
my($dlogfh);
my($err) = 0;

unless (GetDLogFH(\$dlogfh, $daemonlog))
{
    print $dlogfh $msg;
    CloseDLog(\$dlogfh);
}

$cmd = "$binpath/$manage JSOC_DBHOST=$dbhost procser=$procser";

my($msgq) = {lastsend => time(), msgs => {}};

while (1)
{
   # print "running $cmd.\n";
   $rout = qx($cmd 2>&1);
   
   if ($? == -1)
   {
       QueueMessage($msgq, &kMsgType1, "Export Daemon Execution Failure!!", &kMailMessage1);
   } 
   elsif ($? & 127)
   {
      # jsoc_export_manage died in response to an unhandled signal
      my($sig) = $? & 127;
       
       QueueMessage($msgq, &kMsgType1, "Export Daemon Execution Failure!!", &kMailMessage2, "DB Host: $dbhost\n", "Unhandled signal: $sig.\n");
   } 
   elsif (($? >> 8) != 0)
   {
      # jsoc_export_manage returned with an error code
      $msg = "$manage returned with a non-zero code of $? >> 8.\n";
      unless (GetDLogFH(\$dlogfh, $daemonlog))
      {
          print $dlogfh $msg;
      }

      print LOG $msg;
   }

   if (defined($rout) && length($rout) > 0)
   {
      $msg = "$rout\n";
      unless (GetDLogFH(\$dlogfh, $daemonlog))
      {
         print $dlogfh $msg;
      }

      print LOG $msg;
   }

    SendPendingMessages($msgq);
    
   CloseDLog(\$dlogfh);
      
   if (KeepRunning($runningflag))
   {
      sleep(2);
   } 
   else
   {
      last;
   }
}

$msg = "Stopped by $ENV{'USER'} at " . `date` . ".\n";
unless (GetDLogFH(\$dlogfh, $daemonlog))
{
    print $dlogfh $msg;
    CloseDLog(\$dlogfh);
}

print LOG $msg;
close(LOG);

# Don't leave junk laying about
CleanRunFlag($runningflag);

# release the exclusive file lock
flock($lckfh, LOCK_UN);
$lckfh->close;

exit($err);

# END
sub IOwnRunFlag
{
   my($file) = $_[0];
   my($fexists);
   my($iownit);
   my($line);
   
   $fexists = (-e $file);
   if ($fexists)
   {
      if (open(FLFILE, "<$file"))
      {
         $line = <FLFILE>;
         chomp($line);
         $iownit = ($line == $$);
         close(FLFILE);
      }
   }

   return $fexists && $iownit;
}

sub KeepRunning
{
   my($file) = $_[0];

   return IOwnRunFlag($file)
}

sub CleanRunFlag
{
   my($file) = $_[0];

   if (IOwnRunFlag($file))
   {
      unlink($file);
   }
}

sub GetDLogFH
{
    my($rfh) = shift; # reference to filehandle object
    my($dlog) = shift;
    my($msgq) = shift;
    my($err);

    $err = 0;

    if (!defined($$rfh))
    {
        $$rfh = FileHandle->new(">>$dlog");

        if (!defined($$rfh))
        {
            QueueMessage($msgq, &kMsgType1, "Export Daemon Log Unavailable", &kMailMessage3);   
            $err = 1;
        }
    }

    return $err;
}

sub CloseDLog
{
    my($rfh) = $_[0]; # reference to filehandle object 

    if (defined($$rfh))
    {
        $$rfh->close();
        undef($$rfh);
    }
}

sub SendPendingMessages
{
    my($msgs) = shift;
    my($imsg);
    my($msg);
    my($subj);
    
    if (time() - $msgs->{lastsend} > &kMsgQSendInterval)
    {
        # Check for pending messages
        foreach $imsg (keys(%{$msgq->{msgs}}))
        {
            $msg = $msgq->{msgs}->{$imsg}->{msg};
            $subj = $msgq->{msgs}->{$imsg}->{subj};
            open(MAILPIPE, "| /bin/mail -s \"$subj\" " . &kMailList) || die "Couldn't open 'mail' pipe.\n";
            print MAILPIPE $msg;
            close(MAILPIPE);
            
            if ($msgq->{msgs}->{$imsg}->{ntimes} > 1)
            {
                $msgq->{msgs}->{$imsg}->{ntimes} = $msgq->{msgs}->{$imsg}->{ntimes} - 1;
            }
            else
            {
                delete($msgq->{msgs}->{$imsg});
            }
        }
        
        $msgs->{lastsend} = time();
    }
}

# Message queue:
#   key - id
#   val - hash : {instime => 10292392, msg => "export failure", ntimes => 5}
#     where instime is unix seconds identifying time message was inserted into queue.
#           subj is the mail subject
#           msg is the message to mail
#           ntimes is the number of times to send message out.
sub QueueMessage
{
    my($msgq) = shift;
    my($type) = shift;
    my($subj) = shift;
    my(@msg) = @_;
    my($oktoins);
    
    if (exists($msgq->{msgs}->{$type}))
    {
        $oktoins = time() - $msgq->{msgs}->{$type}->{instime} > &kMsgQInterval;
        
        # Message already exists. Don't add to queue until some time elapses.
        if ($oktoins)
        {
            $msgq->{msgs}->{$type}->{ntimes} = $msgq->{msgs}->{$type}->{ntimes} + 1;
        }
    }
    else
    {
        my($msgstr) = join('', @msg);
        $msgq->{msgs}->{$type} = {instime => time(), subj => $subj, msg => $msgstr, ntimes => 1};
    }
}

__DATA__
