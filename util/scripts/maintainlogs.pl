#!/home/jsoc/bin/linux_x86_64/perl5.12.2

# This script can be run on many 'log-sets'. Each log-set contains all log
# files for a single project, tar'ed or otherwise. All log files in a set reside in a single 
# directory. For example, the slony subscription server log-set comprises all
# log files in hmidb2:/usr/local/pgsql/log. 
# There is one file in each log-set that is 'active'. This is the log file
# to which the project code is currently appending. There may be other non-tar'd
# files in the set - the program is free to make as many log files as desired.
# Finally, there may be tar files containing archives of older log files.

# The base log-file name is a C string that identifies all log files
# in a set. All file names contain this string. The active log-file name could 
# be equivalent to the base log-file name, although this isn't assumed. Tar file
# names are created by appending to the base name a time-stamp suffix.

# This script does two things:
#   1. It will archive (move into a tar file) all non-tar files in the log set.
#      It will only do so after a configurable time interval has elapsed - each
#      log-set can have a custom 'archive interval'.
#   2. It will move 'old' tar files to a 'trash bin'. Only files older than a configurable
#      time interval will be moved to the trash.
#   Before this script can move any files into a tar archive, it must acquire a lock
#   so that it doesn't attempt to move files that are currently being used by the program
#   producing them. But this lock is optional. If the configuration for a  log-set does not 
#   specify a lock file, then this script will not attempt to acquire a lock before proceeding.

# This script will maintain a database table of information, one row per log-set.
# This information includes:
#   1. the log-set name (text, prime key)
#   2. the base log-file name (text)
#   3. the path to the directory containing the log-set (text)
#   4. the path to the trash directory (text)
#   5. the path to the lock file - optional (text)
#   6. the arhive interval - amount of time in minutes that elapses between archiving (integer)
#   7. the retention time - interval of time (days) to keep tar'd files in the log-set, after this interval
#      tar'd files are moved to a 'trash' bin (integer)
#   8. the date of the last archive (timestamp - date + time)
#   9. a cmd to run before archiving a log - the archive process will result in the active log
#      being deleted; by running this cmd, any processes that use the log can be invoked
#      one last time on that log.

# To enter log-set-specific information into this database table, use the '-a' flag, followed
# by a string containing the seven comma-separated fields.

use warnings;
use strict;

use IO::Dir;
use FileHandle;
use File::stat;
use File::Copy;
use File::Basename;
use Fcntl ':flock';
use File::lockf;
use Time::HiRes qw(gettimeofday);
use Cwd 'abs_path';
use POSIX qw(strftime);

use DBI;
use DBD::Pg;

use constant kSuccess => 0;
use constant kInvalidArg => 1;
use constant kCantArchive => 2;
use constant kCantTrash => 3;
use constant kCantUpdate => 4;
use constant kAlreadyRunning => 5;
use constant kCantDrop => 6;
use constant kCantCreate => 7;
use constant kSQLFailure => 8;
use constant kProcessorFailure => 9;
use constant kActionFailed => 10;

use constant kTarChunk => 64;

use constant kLogset => 0;
use constant kLSBase => 1;
use constant kLSPath => 2;
use constant kLSTrash => 3;
use constant kLSLock => 4;
use constant kLSArchInt => 5;
use constant kLSRetention => 6;
use constant kLSLastArch => 7;
use constant kLSProcess => 8;


my(@headers) = qw(logset basename path trashdir lckfile archinterval retention lastarch process);

my($klogset) = $headers[kLogset];
my($klsbase) = $headers[kLSBase];
my($klspath) = $headers[kLSPath];
my($klstrash) = $headers[kLSTrash];
my($klslockpath) = $headers[kLSLock];
my($klsarchint) = $headers[kLSArchInt];
my($klsretention) = $headers[kLSRetention];
my($klslastarch) = $headers[kLSLastArch];
my($klsprocess) = $headers[kLSProcess];

my(%fpaths);

my($lckfh);
my($gotlock);
my($nolock);
my($dsn);
my($dbh);
my($sth);
my($conftable);
my($stmnt);
my($rrows);
my($row);
my(@rows);
my(%unarchLogs);
my(%archLogs);
my($tableExists);

my($dbname);
my($dbhost);
my($dbport);
my($dbuser);
my($schemaname);
my($tablename);
my($tarbin);
my($gzbin);
my($action);
my($loglevel);
my($addfile);
my(@addrow);

my($rv) = &kSuccess;

# Allow only one instance of this program 
unless (flock(DATA, LOCK_EX|LOCK_NB)) 
{
   print STDERR "$0 is already running. Exiting.\n";
   exit(kAlreadyRunning);
}

# Read cmd-line arguments
if ($#ARGV < 8)
{
   print STDERR "Improper argument list.\n";
   exit(kInvalidArg);
}

$dbname = $ARGV[0];
$dbhost = $ARGV[1];
$dbport = $ARGV[2];
$dbuser = $ARGV[3];
$schemaname = $ARGV[4];
$tablename = $ARGV[5];
$tarbin = $ARGV[6];
$gzbin = $ARGV[7];
$action = $ARGV[8]; # 'add=<filename>', 'add=A,B,C,...', 'go', 'disp', 'create', 'drop'

if ($#ARGV >= 9)
{
   $loglevel = $ARGV[9];
}
else
{
   $loglevel = 0;
}

$conftable = "$schemaname\.$tablename";

# Parse action
if ($action =~ /^add=(.+)/i)
{
   $addfile = $1;
   if (!(-e $addfile))
   {
      @addrow = split(/,/, $addfile);
      $addfile = "";
      $action = "addrow";
   }
   else
   {
      $action = "addfile";
   }
}
elsif ($action !~ /^go$/i && $action !~ /^disp$/i && $action !~ /^create$/i && $action !~ /^drop$/i)
{
   print STDERR "Invalid action argument '$action'.\n";
   exit kInvalidArg;
}

if (!($loglevel =~ /\d+/ && $loglevel >= 0 && $loglevel <=1))
{
   print STDERR "Invalid log level argument '$loglevel'.\n";
   exit kInvalidArg;
}

# connect to the database
$dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
LogLevPrint("Connection to database with '$dsn' as user '$dbuser' ... ", 1, $loglevel);

# Despite ALL documentation saying otherwise, it looks like the error codes/string
# provided by DBI are all UNDEFINED, unless there is some kind of failure. So, 
# never try to look at $dbh->err or $dbh->errstr if the call succeeded. Of course,
# you can't tell if the call succeeded unless you look at $dbh->err. Stupid.
$dbh = DBI->connect($dsn, $dbuser, ''); # password should be provided by .pgpass

if (defined($dbh))
{
   $tableExists = 0;

   LogLevPrint("success!\n", 1, $loglevel);

   # Check for table existence
   $stmnt = "SELECT * FROM information_schema.tables WHERE table_schema = '$schemaname' AND table_name = '$tablename'";
   $rrows = $dbh->selectall_arrayref($stmnt, undef);

   if (NoErr($rrows, \$dbh, $stmnt))
   {
      my(@rows) = @$rrows;
      if ($#rows == 0)
      {
         # Table exists
         $tableExists = 1;
      }
   }

   if (!$tableExists)
   {
      if ($action eq "create")
      {
         if (!CreateTable(\$dbh, $conftable))
         {
            print STDERR "Unable to create table '$conftable'.\n";
            $rv = &kCantCreate;
         }
         else
         {
            LogLevPrint("Successfully created configuration table '$conftable'.\n", 0, $loglevel);
         }
      }
      else
      {
         print STDERR "Table $schemaname\.$tablename does not exist - nothing to do.\n";
         $rv = &kInvalidArg;
      }
   }
   else
   {
      # Configuration table does exist

      if ($action eq "addrow")
      {
         if (!AddRecord(\$dbh, $conftable, $loglevel, @addrow))
         {
            print STDERR "Unable to add the following record to '$conftable':\n";
            print STDERR "@addrow\n";
            $rv = &kInvalidRow;
         }
         else
         {
            LogLevPrint("Successfully added one row to '$conftable'.\n", 0, $loglevel);
         }
      }
      elsif ($action eq "addfile")
      {
         # Open file and call AddRecord() for each row
         if (open(ADDFILE, "<$addfile"))
         {
            my(@onerow);
            my($line);
            my($atleastone) = 0;

            while (defined($line = <ADDFILE>))
            {
               chomp($line);
               @onerow = split(/,/, $line);
               if (!AddRecord(\$dbh, $conftable, $loglevel, @onerow))
               {
                  print STDERR "Unable to add the following record to '$conftable':\n";
                  print STDERR "$line\n";
                  print STDERR "Continuing with the next line in '$addfile'.\n";
               }
               else
               {
                  $atleastone = 1;
               }
            }

            if (!$atleastone)
            {
               $rv = &kInvalidRow;
            }
            else
            {
               LogLevPrint("Successfully added at least one row to '$conftable'.\n", 0, $loglevel);
            }

            close(ADDFILE);
         }
         else
         {
            print STDERR "Unable to read input file '$addfile'.\n";
            $rv = &kFileIO;
         }
      }
      elsif ($action eq "create")
      {
         print STDERR "Configuration table '$conftable' already exists.\n";
         $rv = &kInvalidArg;
      }
      elsif ($action eq "drop")
      {
         if (!DropTable(\$dbh, $conftable))
         {
            print STDERR "Unable to delete configuration table.\n";
            $rv = &kCantDrop;
         }
         else
         {
            LogLevPrint("Successfully deleted configuration table '$conftable'.\n", 0, $loglevel);
         }
      }
      elsif ($action eq "disp")
      {
         DisplayTable(\$dbh, $conftable);
      } 
      elsif ($action eq "go")
      {
         if (Go(\$dbh, $conftable, \%unarchLogs, \%archLogs, $loglevel) != kSuccess)
         {
            print STDERR "Failure performing 'go' action.\n";
            $rv = &kActionFailed;
         }
      }
      else
      {
         print STDERR "Invalid action '$action'.\n";
         $rv = &kInvalidArg;
      }
   }

   $dbh->disconnect();
}
else
{
   print STDERR "Failure connecting to database: '$DBI::errstr'\n";
}

return $rv;

# FIN

sub AcquireLock
{
   my($path) =$_[0];
   my($lckfh) = $_[1];
   my($gotlock);
   my($natt);

   if (!(-e $path))
   {
      $$lckfh = FileHandle->new(">$path");
      $fpaths{fileno($$lckfh)} = $path;
   }
   else
   {
      $$lckfh = FileHandle->new(">$path");
   }
   $gotlock = 0;

   $natt = 0;
    
   while (1)
   {
       #      if (flock($$lckfh, LOCK_EX|LOCK_NB)) 
      if (!File::lockf::tlock($$lckfh, 0))
      {
         $gotlock = 1;
         last;
      }
      else
      {
         if ($natt < 10)
         {
            LogLevPrint("Lock '$path' in use - trying again in 1 second.\n", 1, $loglevel);
            sleep 1;
         }
         else
         {
            LogLevPrint("Couldn't acquire lock after $natt times; bailing.\n", 1, $loglevel);
            last;
         }
      }

      $natt++;
   }

   return $gotlock;
}

sub ReleaseLock
{
   my($lckfh) = $_[0];
   my($lckfn);
   my($lckfpath);

   $lckfn = fileno($$lckfh);
   $lckfpath = $fpaths{$lckfn};

    #flock($$lckfh, LOCK_UN);
    File::lockf::ulock($$lckfh, 0);
   $$lckfh->close;

   if (defined($lckfpath))
   {
      chmod(0666, $lckfpath);
      delete($fpaths{$lckfn});
   }
}

sub NoErr
{
   my($rv) = $_[0];
   my($dbh) = $_[2];
   my($stmnt) = $_[2];
   my($ok) = 1;

   if (!defined($rv) || !$rv)
   {
      if (defined($$dbh) && defined($$dbh->err))
      {
         print STDERR "Error " . $$dbh->errstr . ": Statement '$stmnt' failed.\n";
      }

      $ok = 0;
   } 

   return $ok;
}

# Logs are returned in an hash array, keyed by the logset name.
# unarchRef - a reference to an hash array of array references. There is one one top-level element
#   per log-set, keyed by log-set name. Each element's value is a reference to a 2-element array. The
#   first element of this array is the path to the files in the log-set, and the second element
#   is an array of files in the path.
#
#{
#    logset1 => [path, [file1, file2, ...]],
#    logset2 => [path, [file1, file2, ...]],
#    ...
#}
#
# archRef - a reference to an hash array of hash array references. There is one top-level element
#   per log-set, keyed by log-set name. Each element's value is a reference to a 2-element array. The
#   first element of this array is the path to the files in the log-set, and the second element
#   is an hash-array reference. The key for this hash is the archive filename, and the value is
#   the timestamp of the archive file.
#
#{
#    logset1 => [path, {tarfile1 => timestamp, tarfile2 => timestamp}],
#    logset2 => [path, {tarfile1 => timestamp, tarfile2 => timestamp}],
#    ...
#}

sub FindLogs
{
    my($logset) = shift;
    my($lspath) = $_[0];
    my($lsbase) = $_[1];
    my($unarchRef) = $_[2];
    my($archRef) = $_[3];
    
    my(@allfiles);
    my(@lfiles);
    my($date);
    
    # First, check to see if we already have the log files for logset
    if (!defined($unarchRef->{$logset}) || !defined($archRef->{$logset}))
    {
        tie(my(%logfiles), "IO::Dir", $lspath);
        @allfiles = keys(%logfiles);
        @lfiles = map({$_ =~ /$lsbase/ ? $_ : ()} @allfiles);
        
        foreach my $onelog (@lfiles)
        {
            if ($onelog =~ /\.tar\.gz$/i)
            {
                # An archive file
                my($sb) = stat("$lspath/$onelog");
                $date = $sb->mtime; # seconds since epoch
                
                if (!defined($archRef->{$logset}))
                {
                    $archRef->{$logset} = []; # anonymous empty array
                    $archRef->{$logset}->[0] = $lspath;
                    $archRef->{$logset}->[1] = {}; # anonymous empty hash
                }
                
                $archRef->{$logset}->[1]->{$onelog} = $date;
            }
            else
            {
                # An unarchived file
                if (!defined($unarchRef->{$logset}))
                {
                    $unarchRef->{$logset} = []; # anonymous empty array
                    $unarchRef->{$logset}->[0] = $lspath;
                    $unarchRef->{$logset}->[1] = [];  # anonymous empty array
                }
                
                my($arrref) = $unarchRef->{$logset}->[1];
                push(@$arrref, $onelog);
            }
        }
        
        untie(%logfiles);
    }
}

sub CreateTable
{
   my($dbh, $conftable) = @_;
   my($stmnt);
   my($rv);
   my($res);

   $rv = 1;

   $stmnt = "CREATE TABLE $conftable ($headers[kLogset] text NOT NULL, $headers[kLSBase] text NOT NULL, $headers[kLSPath] text NOT NULL, $headers[kLSTrash] text, $headers[kLSLock] text, $headers[kLSArchInt] interval DEFAULT '1 day', $headers[kLSRetention] interval DEFAULT '30 day', $headers[kLSLastArch] timestamp with time zone DEFAULT '1966-12-25 00:54:00 PST', $headers[kLSProcess] text, PRIMARY KEY ($klogset))";

   $res = $$dbh->do($stmnt);
   if (!NoErr($res, $dbh, $stmnt))
   {
      $rv = 0;
   }

   return $rv;
}

sub DropTable
{
   my($dbh, $conftable) = @_;
   my($stmnt);
   my($rv);
   my($res);

   $rv = 1;

   $stmnt = "DROP TABLE $conftable";
   $res = $$dbh->do($stmnt);
    if (!NoErr($res, $dbh, $stmnt))
   {
      $rv = 0;
   }

   return $rv;
}

sub DisplayTable
{
   my($dbh, $conftable) = @_;
   my($stmnt);
   my($rrows);

   $stmnt = "SELECT $headers[kLogset], $headers[kLSBase], $headers[kLSPath], $headers[kLSTrash], $headers[kLSLock], $headers[kLSArchInt], $headers[kLSRetention], $headers[kLSLastArch], $headers[kLSProcess] from $conftable";

   $rrows = $$dbh->selectall_arrayref($stmnt, undef);

   if (NoErr($rrows, $dbh, $stmnt))
   {
      my($ptbuf);

      # Print header
      $ptbuf = sprintf("%16s%16s%48s%48s%64s%16s%16s%32s%32s", $headers[kLogset], $headers[kLSBase], $headers[kLSPath], $headers[kLSTrash], $headers[kLSLock], $headers[kLSArchInt], $headers[kLSRetention], $headers[kLSLastArch], $headers[kLSProcess]);
      print "$ptbuf\n";
      print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";

      if ($#$rrows >= 0)
      {
         foreach $row (@$rrows)
         {
            my($logset, $lsbase, $lspath, $lstrash, $lslockpath, $lsarchint, $lsretention, $lslastarch, $lsprocess) = @$row;
            $ptbuf = sprintf("%16s%16s%48s%48s%64s%16s%16s%32s%32s", $logset, $lsbase, $lspath, defined($lstrash) ? $lstrash : "", defined($lslockpath) ? $lslockpath : "", $lsarchint, $lsretention, $lslastarch, defined($lsprocess) ? $lsprocess : "");
            print "$ptbuf\n";
         }
      }
   }
}

sub AddRecord
{
   my($ok);
   my($dbh, $conftable, $loglev, @rowin) = @_;
   my(@row);
   my(@colnames);
   my($stmnt);
   my($rv);
   my($ival);
   my($rowstr);
   my($colnamesstr);

   #$logset, $lsbase, $lspath, $lstrash, $lslockpath, $lsarchint, $lsretention, $lslastarch, $lsprocess

   $ival = 0;
   foreach my $item (@rowin)
   {
      if (defined($item))
      {
         push(@row, $item);
         push(@colnames, $headers[$ival]);
      }

      $ival++;
   }

   # Perl mumbo jumbo
   $" = ',';
   $rowstr = "@row";
   $colnamesstr = "@colnames";

   $stmnt = "INSERT into $conftable ($colnamesstr) VALUES ($rowstr)";
   $rv = $$dbh->do($stmnt);

   LogLevPrint("$stmnt\n", 1, $loglev);

   $ok = NoErr($rv, $dbh, $stmnt);

   return $ok;
}

sub Go
{
   my($dbh, $conftable, $unarchLogs, $archLogs, $loglev) = @_;
   my($stmnt);
   my($rrows);
   my($lckfh);
   my($gotlock);
   my($nolock);
   my($row);
   my($atleastone);
   my($procfailed);
   my($res);
   my($rv);

   $rv = kSuccess;

   # First, archive the non-archived logs, if it is time to do so.
   $stmnt = "SELECT $headers[kLogset], $headers[kLSBase], $headers[kLSPath], $headers[kLSLock], $headers[kLSProcess] FROM $conftable WHERE $headers[kLSLastArch] + $headers[kLSArchInt] < current_timestamp";
   $rrows = $$dbh->selectall_arrayref($stmnt, undef);

   if (NoErr($rrows, $dbh, $stmnt))
   {
      # $rrows will contain rows for only those log-sets that require archiving
      $atleastone = 0;

      foreach $row (@$rrows)
      {
         my($logset, $lsbase, $lspath, $lslockpath, $lsprocess) = @$row;
          
         LogLevPrint("Archiving $logset logset.\n", 0, $loglev);

         if (defined($lslockpath) && length($lslockpath) > 0)
         {
            $gotlock = AcquireLock($lslockpath, \$lckfh);
            $nolock = 0;
         } 
         else
         {
            $gotlock = 0;
            $nolock = 1;
         }

         if (!$nolock && !$gotlock)
         {
            print STDERR "Skipping $logset - unable to acquire lock.\n";
         } 
         else
         {
            # Archive the current logs for $logset
            FindLogs($logset, $lspath, $lsbase, $unarchLogs, $archLogs);

            if (defined($unarchLogs->{$logset}))
            {
               $procfailed = 0;

               # We have a log file to archive. After the file is archived, it will be deleted, so 
               # give the processing program a last chance to run.
               if (defined($lsprocess) && length($lsprocess) > 0)
               {
                  if (!CallCmd($lsprocess))
                  {
                     print STDERR "Unable to run '$lsprocess' properly.\n";
                     $procfailed = 1;
                  }
               }

               if (!$procfailed)
               {
                  LogLevPrint("Archiving logset $logset (base filename: $lsbase).\n", 0, $loglev);
                  if (!Archive($lspath, $lsbase, $tarbin, $gzbin, $unarchLogs->{$logset}->[1], $loglev))
                  {
                     print STDERR "Unable to archive log files for logset '$logset'.\n";
                  }
                  else
                  {
                     $atleastone = 1;
                  }
               }
            }
         }
          
         # Update lastarch value
         my($timenow) = strftime("%a %b %e %H:%M:%S %Y", localtime());
          
         $stmnt = "UPDATE $conftable SET $headers[kLSLastArch] = '$timenow' WHERE $headers[kLogset] = '$logset'";
         $res = $$dbh->do($stmnt);
          
         if (!NoErr($res, $dbh, $stmnt))
         {
             print STDERR "Unable to update '$klslastarch' in '$conftable'.\n";
             $rv = kCantUpdate;
             last;
         }

         # Don't hold the lock - allow the processes generating the logs a chance to write in between the processing
         # of logsets.
         if ($gotlock)
         {
            ReleaseLock(\$lckfh);
            $gotlock = 0;
         }
      } # loop over logsets
       
      if ($rv == &kSuccess)
      {
          if (!$atleastone && $#$rrows >= 0 && (keys(%$unarchLogs) > 0))
          {
              print STDERR "Failure archiving at least one log-set.\n";
              $rv = kCantArchive;
          }
          elsif ($#$rrows < 0)
          {
              LogLevPrint("Not time to archive any logs.\n", 0, $loglev);
          }
          elsif ((keys(%$unarchLogs) < 1))
          {
              LogLevPrint("No log files available for archiving.\n", 0, $loglev);
          }
      }
   }
   else
   {
      $rv = kSQLFailure;
   }

   if ($rv == kSuccess)
   {
      my($timetotrash) = 0;

      # Second, trash the archive files that have expired.
      $stmnt = "SELECT $headers[kLogset], $headers[kLSBase], $headers[kLSPath], $headers[kLSTrash], $headers[kLSLock], $headers[kLSArchInt], EXTRACT ('epoch' FROM $headers[kLSRetention]) AS $headers[kLSRetention], $headers[kLSLastArch] from $conftable";

      $rrows = $$dbh->selectall_arrayref($stmnt, undef);

      if (NoErr($rrows, $dbh, $stmnt))
      {
         $atleastone = 0;

         foreach $row (@$rrows)
         {
            my($logset, $lsbase, $lspath, $lstrash, $lslockpath, $lsarchint, $lsretention, $lslastarch) = @$row;

            if (defined($lslockpath) && length($lslockpath) > 0)
            {
               $gotlock = AcquireLock($lslockpath, \$lckfh);
               $nolock = 0;
            } 
            else
            {
               $gotlock = 0;
               $nolock = 1;
            }

            if (!$nolock && !$gotlock)
            {
               print STDERR "Skipping $logset - unable to acquire lock.\n";
            } 
            else
            {
               # Trash archives that have expired (based on file timestamp, NOT lastarch).
               my(@totrash);
               my($onearr);
               my($cpath);
               my($cfilesH);

               FindLogs($logset, $lspath, $lsbase, $unarchLogs, $archLogs);

               # Find out which archives are old enough to trash
               if (defined($archLogs->{$logset}))
               {
                  $onearr = $archLogs->{$logset}; # returns an array reference
                  $cpath = $onearr->[0];
                  $cfilesH = $onearr->[1];
                   
                  foreach my $key (keys(%$cfilesH))
                  {
                     if ($cfilesH->{$key} + $lsretention < time())
                     {
                        push(@totrash, "$cpath/$key");
                     }
                  }
               }

               if ($#totrash >= 0)
               {
                  $timetotrash = 1;

                  if (-d $lstrash)
                  {
                     # Trashing a single log-set's old archives
                     if (!Trash($lstrash, $loglev, @totrash))
                     {
                        print STDERR "Unable to move expired archives in '$lspath' to trash bin.\n";
                     }
                     else
                     {
                        $atleastone = 1;
                     }
                  }
                  else
                  {
                     print STDERR "Can't access trash bin '$lstrash'.\n";
                  }
               }

               if ($gotlock)
               {
                  ReleaseLock(\$lckfh);
                  $gotlock = 0;
               }
            }
         } # foreach

         if (!$atleastone && $timetotrash)
         {
            print STDERR "Unable to trash any old archives.\n";
            $rv = kCantTrash;
         }
         elsif (!$timetotrash && keys(%$archLogs) > 0)
         {
            LogLevPrint("Not time to trash any archives.\n", 1, $loglev);
         }
         elsif ((keys(%$archLogs) < 1))
         {
            LogLevPrint("No archives available for trashing.\n", 1, $loglev);
         }
      }
   }

   return $rv;
}

sub CallCmd
{
   my($cmd) = $_[0];

   my($rv) = 1;
   my($callstat);

   system($cmd);
   $callstat = $?;

   if ($callstat == -1)
   {
      print STDERR "Failed to execute '$cmd'.\n";
      $rv = 0;
   } 
   elsif ($callstat & 127)
   {
      print STDERR "command '$cmd' terminated abnormally.\n";
      $rv = 0;
   } 
   elsif ($callstat >> 8 != 0)
   {
      $callstat = $callstat >> 8;
      print STDERR "command '$cmd' ran unsuccessfully, status == $callstat.\n";
      $rv = 0;
   }

   return $rv;
}

sub Archive
{
   my($lspath, $lsbase, $tarbin, $gzbin, $rtotar, $loglev) = @_;
   my(@totar);
   my($rpath);
   my($isfile);
   my($fullpath);
   my($tarmcd);
   my($res);
   my($sfilelist);
   my(@ltime);
   my($datestr);
   my($rv);

   $rv = 1;

   @totar = @$rtotar;
   $rpath = abs_path(); # current wd
   chdir($lspath);

   # Clean up anything that might have gotten left behind during previous incomplete runs
   if (-e ".tmp.tar")
   {
      unlink(".tmp.tar");
   }

   if (-e ".tmp.tar.gz")
   {
      unlink(".tmp.tar.gz");
   }

   if ($#totar >= 0)
   {
      # Tar the first log file, because tar cannot create an empty archive
      $fullpath = $totar[0];
      $sfilelist = "";
      RunTar($tarbin, "cf", ".tmp.tar", "", $fullpath, $loglev);

      # Chunk tar cmds so that we don't overrun the cmd-line with too many chars.
      $isfile = 1;

      foreach $fullpath (@totar[1..$#totar])
      {
         $sfilelist = "$sfilelist $fullpath";
         if ($isfile % kTarChunk == 0)
         {
            # Execute tar cmd.
             RunTar($tarbin, "rf", ".tmp.tar", "", $sfilelist, $loglev);
            $sfilelist = "";
         }

         $isfile++;
      }

      if (length($sfilelist) > 0)
      {
          RunTar($tarbin, "rf", ".tmp.tar", "", $sfilelist, $loglev);
      }

      # Calculate a date string to insert into the tar's filename
      @ltime = localtime();
      $datestr = sprintf("%4d%02d%02d_%02d%02d%02d", $ltime[5] + 1900, $ltime[4] + 1, $ltime[3], $ltime[2], $ltime[1], $ltime[0]);

      # Validate tar file
      LogLevPrint("Validating tar '$lspath/.tmp.tar'.\n", 0, $loglev);
      $res = Validatetar(".tmp.tar", @totar);

      if ($res == 1)
      {
         my($gzcmd);
         my($tarfilename) = sprintf("%s_%s.tar.gz", $lsbase, $datestr);

         $gzcmd = "$gzbin --best .tmp.tar";
         LogLevPrint("Gzipping '$lspath/.tmp.tar' to '$lspath/.tmp.tar.gz'\n", 0, $loglev);

         $res = system($gzcmd);
         if ($res == 0) 
         {
            LogLevPrint("Moving '$lspath/.tmp.tar.gz' to '$lspath/$tarfilename'.\n", 0, $loglev);
            if (!move(".tmp.tar.gz", $tarfilename))
            {
               print STDERR "Unable to move .tmp.tar.gz to $tarfilename.\n";
               $rv = 0;
            }
            else
            {
               # Remove log files just tarred.
               while (defined($fullpath = shift(@totar)))
               {
                  LogLevPrint("Deleting tarred log file: $fullpath.\n", 0, $loglev);
                  unlink($fullpath);
               }
            }
         }
         else
         {
            print STDERR "gzip cmd '$gzcmd' failed to run properly.\n";
            $rv = 0;
         }
      }
      else
      {
         $rv = 0;
      }
   }

   chdir($rpath);
   return $rv;
}

sub RunTar
{
   my($tarbin, $op, $file, $list, $options, $loglev) = @_;
   my($tarcmd);
   my($res);
   my($rv);

   $rv = 1;

   $tarcmd = "$tarbin $op $file $list $options";
   LogLevPrint("Running tar '$tarcmd'\n", 0, $loglev);

   $res = system($tarcmd);
   ($res == 0) || ($rv = 0);

   return $rv;
}

sub Validatetar
{
   my($tarname, @logfiles) = @_;
   my($res);
   my($error);
   my($rv);

   $rv = 1;

   # examine contents
   $res = `$tarbin tf $tarname`;
   $error = $? >> 8;

   if ($error != 0)
   {
      print STDERR "ERROR executing '$tarbin tf $tarname', error '$error'.\n";
      $rv = 0;
   } 
   else 
   {
      my(@tfiles); # tar files
      my($tcount) = 0;
      my(@sfiles); # orig files
      my($scount) = 0;
      
      @tfiles = sort(map({basename($_)} split(/\s*\n\s*/, $res)));
      @sfiles = sort(map({basename($_)} @logfiles));

      if ($#tfiles >= 0)
      {
         $tcount = $#tfiles + 1;
         $scount = $#sfiles + 1;
         
         if ($tcount != $scount)
         {
            print STDERR "Number of logfiles in tar ($tcount) does not match the expected number ($scount).\n";
            $rv = 0;
         }
         else
         {
            # compare entire arrays, line by line
            my($ielem) = 0;
            while ($ielem < $tcount)
            {
               
               if ($tfiles[$ielem] ne $sfiles[$ielem])
               {
                  print STDERR "Tar validation failed.\n";
                  $rv = 0;
                  last;
               }

               $ielem++;
            }
         }
      }
      else
      {
         print STDERR "Unexpected empty tar.\n";
         $rv = 0;
      }
   }

   return $rv;
}

sub Trash
{
   my($trashbin, $loglev, @files) = @_;
   my($rv) = 1;

   if (-d $trashbin)
   {
      foreach my $afile (@files)
      {
         LogLevPrint("Trashing archive $afile.\n", 0, $loglev);

         # Move the file to the trash bin
         if (!move($afile, $trashbin))
         {
            print STDERR "Couldn't move file '$afile' into trashbin '$trashbin'.\n";
            $rv = 0;
            last;
         }
      }
   }
   else
   {
      print STDERR "'$trashbin' is not a directory.\n";
      $rv = 0;
   }

   return $rv;
}

sub LogLevPrint
{
   my($msg) = $_[0];
   my($printlev) = $_[1];
   my($loglev) = $_[2];

   if ($loglev >= $printlev)
   {
      print $msg;
   }
}
__DATA__
