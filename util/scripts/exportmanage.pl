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
#        Example : /home/jsoc/cvs/Development/JSOC/proj/util/scripts/exportmanage.pl -root /home/arta/cvs/JSOC -dbuser production -dbhost hmidb2 -manager "jsoc_export_manage -t" -runflag keepruntest.txt &
#
#  2. Point a browser at http://jsoc.stanford.edu/ajax/exportdatatest.html and export something.

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
my($dbhost) = "hmidb";
my($dbname) = "jsoc";
my($dbuser) = "production";
my($binpath);
my($manage) = "jsoc_export_manage";
my($logfile);

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
    elsif ($arg eq "-jsocdev")
    {
        $root = $kJSOCDEV_ROOT;
        $binpath = "$root/bin";
        $dbuser = $kJSOCDEV_DBUSER;
        $dbname = $kJSOCDEV_DBNAME;
        $dbhost = $kJSOCDEV_DBHOST;
        $manage = $kJSOCDEV_MANAGE;
        $runningflag = $kINTERNALFLAG;
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
#    print "running $binpath/$manage JSOC_DBHOST=\"$dbhost\"\n";
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
