#!/home/jsoc/bin/linux_x86_64/activeperl -w

# udsumsforarch.pl jsoc_sums hmidb 5434 production many /home/arta/Projects/SUMS/UpdateSUMSArch/testsm.txt

use DBI;
use DBD::Pg;
# ART - This module may not be installed. DO NOT USE CPAN to install modules. Use the ActiveState ppm. Unfortunately,
# ActiveState does not provide this module.
use Switch;
use POSIX qw(strftime);
use File::Copy "mv";

use constant kDEBUG => 0;

# SUs that are archive pending must be set to DADP
use constant kStatDADP => "2";

use constant kTQueryOnePerConn => "conn";
use constant kTQueryOnePerTrans => "xact";
use constant kTQueryManyPerTrans => "many";

use constant kUdListTableTapeFileInfo => "tmp_ulist_tapefiles";
use constant kInTypeList => 0;
use constant kInTypeMd5 => 1;
use constant kDelim => '|';

use constant kGTarBlock => "256";

use constant kDbTableSumMain => "sum_main";
use constant kDbTableSumPartnAlloc => "sum_partn_alloc";
use constant kDbTableSumFile => "sum_file";

use constant kDbTableSumMainDEBUG => "arta_main";
use constant kDbTableSumPartnAllocDEBUG => "arta_partn_alloc";
use constant kDbTableSumFileDEBUG => "arta_file";

my($err);

my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)
my($typequery); # one SU per connection, one SU per transaction, or one big query that does all SUs at the same time.
my($dfilelist); # list of file paths - each file contains data to update sum_main with.

my($dsn);       # database connection arguments
my($dbh);       # database handle
my($stmnt);
my($row);
my($rrows);
my($line);
my($rv);

if ($#ARGV < 5)
{
   print "Improper argument list.\n"; 
   exit(1);
}

$dbname = $ARGV[0];
$dbhost = $ARGV[1];
$dbport = $ARGV[2];
$dbuser = $ARGV[3];
$typequery = $ARGV[4];
$dfilelist = $ARGV[5];

$err = 0; 

if ($typequery eq kTQueryOnePerConn)
{
   print STDERR "One query per database connection unimplemented.\n";
   $err = 1;
}
elsif ($typequery eq kTQueryOnePerTrans)
{
   print STDERR "One query per database transaction unimplemented.\n";
   $err = 1;
}
elsif ($typequery eq kTQueryManyPerTrans)
{
   # batch sunums together
   my($fpath);
   my($ftype);
   my($timenow);
   my($tapeid);
   my($ttabrow);
   my(@md5csums);
   my($values);
   my($skip);
   my($maintable);
   my($partntable);
   my($filetable);

   if (dbconnect($dbname, $dbhost, $dbport, \$dbh))
   {
      # Open input file for reading - has 2 columns, group and filepath (path to a file provided by Keh-Cheng)
      if (open(FILELIST, "<$dfilelist"))
      {
         # create temporary table to hold LIST files (series, sunum, sudir, fileid, tapeid)
         $stmnt = "CREATE TEMPORARY TABLE " . kUdListTableTapeFileInfo . "(sunum bigint, fileid integer, tapeid character varying(20))";
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create temporary table " . "'kUdListTableTapeFileInfo'" . ".\n");

         $stmnt = "CREATE INDEX " . kUdListTableTapeFileInfo . "_idx on " . kUdListTableTapeFileInfo . "(sunum)";   
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create index on temporary table " . "'kUdListTableTapeFileInfo'" . ".\n");

         CommitTrans(\$dbh);

         if (kDEBUG)
         {
            $maintable = kDbTableSumMainDEBUG;
            $partntable = kDbTableSumPartnAllocDEBUG;
            $filetable = kDbTableSumFileDEBUG;
         }
         else
         {
            $maintable = kDbTableSumMain;
            $partntable = kDbTableSumPartnAlloc;
            $filetable = kDbTableSumFile;
         }
         
         # loop through filelist
         while (defined($line = <FILELIST>))
         {
            chomp($line);

            if ($line =~ /^\#/)
            {
               next;
            }

            if ($line =~ /\s*(\S+)/)
            {
               $fpath = $1;

               if ($fpath =~ /LIST\.(\w+)/)
               {
                  $tapeid = $1;
                  $ftype = kInTypeList;
               }
               elsif ($fpath =~ /MD5SUM\.(\w+)/)
               {
                  $tapeid = $1;
                  $ftype = kInTypeMd5;
               }
               else
               {
                  print STDERR "unsupported input file '$fpath'.\n";
                  next;
               }
            }
            else
            {
               print STDERR "skipping invalid input file $line.\n";
               next;
            }

            # There are two types of files that Keh-Cheng is producing:
            #   1. LIST files - cols: series, sunum, sudir, fileid, tapeid
            #   2. MD5SUM files - cols: fileid, md5sum
            $skip = 0;
            if (defined($fpath) && open(DATAFILE, "<$fpath"))
            {
               if ($ftype == kInTypeList)
               {
                  $stmnt = "COPY " . kUdListTableTapeFileInfo . " (sunum, fileid, tapeid) FROM stdin WITH DELIMITER '" . kDelim . "'";
                  ExecStatement(\$dbh, $stmnt, 1, "Troubles copying data to temporary table '" . kUdListTableTapeFileInfo . "'.\n");
               }

               while (defined($line = <DATAFILE>))
               {
                  chomp($line);

                  # The data file might have extra fields not needed for the update - strip those out
                  if ($ftype == kInTypeList)
                  {
                     if ($line =~ /\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)/)
                     {
                        $ttabrow = "$1|$2|$3";
                     }

                     print "Putting row ($ttabrow).\n";
                     $rv = $dbh->pg_putcopydata("$ttabrow\n");
                     if (!$rv)
                     {
                        print STDERR "Failure sending line '$line' to db.\n";
                        $skip = 1;
                        last;
                     }
                  }
                  else
                  {
                     # Don't save to a temporary table - this file will be relatively small
                     if ($line =~ /(\S+)\s+(\S+)/)
                     {
                        # (tapeid, filenum, gtarblock, md5cksum)
                        $values = "('$tapeid', $1, " . kGTarBlock . ", '$2')";
                        push(@md5csums, $values);
                     }
                  }
               }
               
               if ($ftype == kInTypeList)
               {
                  $rv = $dbh->pg_putcopyend();
                  if (!$rv)
                  {
                     print STDERR "Failure ending table copy.\n";
                     $skip = 1;
                  }
                  else
                  {
                     if (kDEBUG && 0)
                     {
                        # test to see if this worked.
                        $stmnt = "SELECT * from " . kUdListTableTapeFileInfo;
                        $rrows = $dbh->selectall_arrayref($stmnt, undef);
                        $err = !(NoErr($rrows, \$dbh, $stmnt));
                  
                        foreach $row (@$rrows)
                        {
                           print "tmp table row: $row->[0], $row->[1], $row->[2]\n";
                        }
                     }
                  }
               }

               close(DATAFILE);
            }
            else
            {
               print "Unable to open '$fpath' for reading.\n";
               next;
            }

            if (!$skip)
            {
               CommitTrans(\$dbh); # successfully populated temporary table ($ftype == kInTypeList).
                                   # when $ftype == kInTypeMd5, this will result in and 'end' being
                                   # issued right after a 'begin', but this is harmless.

               if ($ftype == kInTypeList)
               {
                  # Create date string
                  $timenow = strftime("%a %b %e %H:%M:%S %Y", localtime());

                  # Do a join that will update sum_main with data from temp table
                  # SQL to update sum_main:
                  #   archive_status -> character varying(5)
                  #   arch_tape ->  character varying(20)
                  #   arch_tape_fn -> integer
                  #   arch_tape_date -> timestamp(0) without time zone 
                  $stmnt = "UPDATE $maintable m SET archive_status = 'Y', arch_tape = ul.tapeid, arch_tape_fn = ul.fileid, arch_tape_date = '$timenow' FROM " . kUdListTableTapeFileInfo . " ul WHERE (m.ds_index = ul.sunum)";
                  ExecStatement(\$dbh, $stmnt, 1, "Troubles updating $maintable.\n");

                  # SQL to update sum_partn_alloc:
                  #   status -> integer
                  #   ds_index -> bigint
                  $stmnt = "UPDATE $partntable m SET status = " . kStatDADP . " FROM " . kUdListTableTapeFileInfo . " ul WHERE (m.ds_index = ul.sunum)";
                  ExecStatement(\$dbh, $stmnt, 1, "Troubles updating $partntable.\n");

                  # Drop rows from temp table.
                  $stmnt = "DELETE from " . kUdListTableTapeFileInfo;
                  ExecStatement(\$dbh, $stmnt, 1, "Couldn't drop rows from temporary table '" . kUdListTableTapeFileInfo . "'.\n");
               }
               else
               {
                  # Insert directly into sum_file (we could also use COPY to copy into a temp table, 
                  # then use a query and constant values for tapeid and gtarblock with the INSERT INTO statement):
                  #   tapeid -> character varying(20)
                  #   filenum -> integer
                  #   gtarblock -> integer
                  #   md5cksum -> character varying(36)
                  foreach $values (@md5csums)
                  {
                     $stmnt = "INSERT INTO $filetable (tapeid, filenum, gtarblock, md5cksum) VALUES $values";
                     ExecStatement(\$dbh, $stmnt, 1, "Troubles updating $filetable.\n");
                  }

                  @md5csums = ();
               }

               # Done updating dbase for $fpath - move to ${fpath}.ingested
               CommitTrans(\$dbh);

               if (mv($fpath, "$fpath\.ingested") == 0)
               {
                  print STDERR "Failed to rename $fpath to $fpath\.ingested.\n";
               }
               else
               {
                  print "Renamed $fpath to $fpath\.ingested\n";
               }
            }
            else
            {
               RollbackTrans(\$dbh); # problem populating temporary table ($ftype == kInTypeList)
            }
         } # loop over data files

         close(FILELIST);
      }

      $dbh->disconnect();
   }
}

# AL FINAL
exit($err);

# SUBROUTINES
sub dbconnect
{
   my($dbname) = $_[0];
   my($dbhost) = $_[1];
   my($dbport) = $_[2];
   my($dbh) = $_[3]; # returned by reference

   my($dsn);

   # connect to the database
   $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
   print "Connection to database with '$dsn' as user '$dbuser' ... ";
   
   # Despite ALL documentation saying otherwise, it looks like the error codes/string
   # provided by DBI are all UNDEFINED, unless there is some kind of failure. So, 
   # never try to look at $dbh->err or $dbh->errstr if the call succeeded.
   
   $$dbh = DBI->connect($dsn, $dbuser, '', { AutoCommit => 0 }); # will need to put pass in .pg_pass

   if (defined($dbh))
   {
      print "success!\n";
   }
   else
   {
      print "failure!!!!\n";
      $err = 1;
   }
}

sub ExecStatement
{
   my($dbh, $stmnt, $doit, $msg) = @_;
   my($res);

   print "executing db statement ==> $stmnt\n";

   if ($doit)
   {
      $res = $$dbh->do($stmnt);
      NoErr($res, $dbh, $stmnt) || die $msg;
   }
}

sub BeginTrans
{
   my($dbh) = $_[0];

   print "starting transaction\n";
   ExecStatement($dbh, "BEGIN", 1, "Couldn't begin transaction.\n");
}

sub CommitTrans
{
   my($dbh) = $_[0];

   print "committing transaction\n";
   ExecStatement($dbh, "COMMIT", 1, "Couldn't begin transaction.\n");
}

sub RollbackTrans
{
   my($dbh) = $_[0];

   print "rolling-back transaction\n";
   ExecStatement($dbh, "ROLLBACK", 1, "Couldn't begin transaction.\n");
}

sub NoErr
{
   my($rv) = $_[0];
   my($dbh) = $_[1];
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
