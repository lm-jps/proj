#!/usr/bin/perl -w 

# Here's how to run this script (from scratch)
#  ssh jsoc@j0
#  cd /home/jsoc/exports
#  rm keep_running
#  exportmanage.pl -jsocdev &
# To start the web version for jsoc.stanford.edu repeat the above but do:
#  rm keep_running_web
#  exportmanage.pl -jsocweb &
# At some point, I'll add the production version of the script, which
#  would then be run as "exportmanage.pl -jsocpro &"
# To start the debug test version for jsoc2 do:
#  rm keep_running_test
#  exportmanage.pl -jsoctest &

my($kINTERNALFLAG) = "/home/jsoc/exports/keep_running";
my($kWEBFLAG) = "/home/jsoc/exports/keep_running_web";
my($kTESTFLAG) = "/home/jsoc/exports/keep_running_test";
#
my($kJSOCDEV_ROOT) = "/home/jsoc/cvs/Development/JSOC";
my($kJSOCDEV_DBUSER) = "production";
my($kJSOCDEV_DBNAME) = "jsoc";
my($kJSOCDEV_DBHOST) = "hmidb";
my($kJSOCDEV_MANAGE) = "jsoc_export_manage";
#
my($kJSOCPRO_ROOT) = "/home/jsoc/cvs/JSOC";
my($kJSOCPRO_DBUSER) = "production";
my($kJSOCPRO_DBNAME) = "jsoc";
my($kJSOCPRO_DBHOST) = "hmidb";
my($kJSOCPRO_MANAGE) = "jsoc_export_manage";
#
my($kJSOCWEB_ROOT) = "/home/jsoc/cvs/Development/JSOC";
my($kJSOCWEB_DBUSER) = "production";
my($kJSOCWEB_DBNAME) = "jsoc";
my($kJSOCWEB_DBHOST) = "hmidb2";
my($kJSOCWEB_MANAGE) = "jsoc_export_manage";
#
my($kJSOCTEST_ROOT) = "/home/jsoc/cvs/Development/JSOC";
my($kJSOCTEST_DBUSER) = "phil";
my($kJSOCTEST_DBNAME) = "jsoc";
my($kJSOCTEST_DBHOST) = "hmidb";
my($kJSOCTEST_MANAGE) = "jsoc_export_manage_test";

$runningflag = $kINTERNALFLAG;
my($arg);
my($root);
my($dbhost);
my($dbname);
my($dbuser);
my($binpath);
my($manage);
my($logfile);

while ($arg = shift(@ARGV))
{
    if ($arg eq "-root")
    {
        $root = shift(@ARGV);
    }
    if ($arg eq "-dbhost")
    {
        $dbhost = shift(@ARGV);
    }
    if ($arg eq "-dbuser")
    {
        $dbuser = shift(@ARGV);
    }
    if ($arg eq "-dbname")
    {
        $dbname = shift(@ARGV);
    }
    if ($arg eq "-jsocdev")
    {
        $root = $kJSOCDEV_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCDEV_DBUSER;
        $dbname = $kJSOCDEV_DBNAME;
        $dbhost = $kJSOCDEV_DBHOST;
        $manage = $kJSOCDEV_MANAGE;
        $runningflag = $kINTERNALFLAG;
    }
    if ($arg eq "-jsocpro")
    {
        $root = $kJSOCPRO_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCPRO_DBUSER;
        $dbname = $kJSOCPRO_DBNAME;
        $dbhost = $kJSOCPRO_DBHOST;
        $manage = $kJSOCPRO_MANAGE;
        $runningflag = $kINTERNALFLAG;
    }
    if ($arg eq "-jsocweb")
    {
        $root = $kJSOCWEB_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCWEB_DBUSER;
        $dbname = $kJSOCWEB_DBNAME;
        $dbhost = $kJSOCWEB_DBHOST;
        $manage = $kJSOCWEB_MANAGE;
        $runningflag = $kWEBFLAG;
    }
    if ($arg eq "-jsoctest")
    {
        $root = $kJSOCTEST_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCTEST_DBUSER;
        $dbname = $kJSOCTEST_DBNAME;
        $dbhost = $kJSOCTEST_DBHOST;
        $manage = $kJSOCTEST_MANAGE;
        $runningflag = $kTESTFLAG;
    }
}

# Only run on j0.Stanford.EDU
# if ($ENV{HOSTNAME} ne "j0.Stanford.EDU") {
#    die "I will only run on j0.Stanford.EDU\n";
#}

# Don't run if somebody is already managing the export
if (-e $runningflag)
{
    die "Can't manage export; another process is already managing it.\n";
}

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
$logfile = "/home/jsoc/exports/logs/" . `date +"%F_%R.log"`;
open(LOG, ">> $logfile") || die "Couldn't open logfile '$logfile'.\n";

my($datenow) = `date`;
chomp($datenow);
print LOG "Started by $ENV{'USER'} at $datenow on machine $ENV{'HOST'} using $dbhost.\n";

while (1)
{
    print LOG `$binpath/$manage JSOC_DBHOST="$dbhost"`;
    if (-e $runningflag) 
    {
        sleep(2);
    }
    else
    {
        print LOG "Stopped by $ENV{'USER'} at " . `date` . ".\n";
        close(LOG);
        exit(0);
    }
}
