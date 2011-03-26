#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

# udsumsforarch.pl jsoc_sums hmidb 5434 production many /home/arta/Projects/SUMS/UpdateSUMSArch/filelist.txt

use DBI;
use DBD::Pg;
use Time::localtime;
use Switch;

use constant kDEBUG => 1;

use constant kStatDADP => "2";
use constant kStatDAAP => "4"; # archive pending
use constant kSubStatDAADP => "128"; # after archive completes, mark delete pending
use constant kGig => 1073741824;

use constant kTQueryOnePerConn => "conn";
use constant kTQueryOnePerTrans => "xact";
use constant kTQueryManyPerTrans => "many";

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
   $stmnt = "SELECT archive_status, arch_tape, arch_tape_fn arch_tape_date from sum_main where ds_index=123456";

   foreach $elem (@what)
   {
      if (dbconnect($dbname, $dbhost, $dbport, \$dbh))
      {
          $rrows = $dbh->selectall_arrayref($stmnt, undef);
          $err = !(NoErr($rrows, \$dbh, $stmnt));

          if (!$err)
          {
             
          }
      }
   }

}
elsif ($typequery eq kTQueryOnePerTrans)
{
   if (dbconnect($dbname, $dbhost, $dbport, \$dbh))
   {
      foreach $elem (@what)
      {
         
         $rrows = $dbh->selectall_arrayref($stmnt, undef);
         $err = !(NoErr($rrows, \$dbh, $stmnt));

         if (!$err)
         {
             
         }

      }

   }
}
elsif ($typequery eq kTQueryManyPerTrans)
{
   # batch sunums together
   my($group);
   my($fpath);
   
   if (dbconnect($dbname, $dbhost, $dbport, \$dbh))
   {
      # Open input file for reading - has 2 columns, group and filepath (path to a file provided by Keh-Cheng)
      if (open(FILELIST, "<$dfilelist"))
      {
         $stmnt = "CREATE TEMPORARY TABLE arta_updatelist (sunum bigint, loc character varying(80), bytes bigint)";
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create temporary table 'arta_updatelist'.\n");

         $stmnt = "CREATE INDEX arta_updatelist_idx on arta_updatelist (sunum)";   
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create index on temporary table 'arta_updatelist'.\n");

         # loop through filelist
         while (defined($line = <FILELIST>))
         {
            chomp($line);

            if ($line =~ /^\#/)
            {
               next;
            }

            if ($line =~ /(\S+)\s+(\S+)/)
            {
               $group = $1;
               $fpath = $2;
            }
            else
            {
               print STDERR "skipping invalid line $line.\n";
               next;
            }

            if (open(DATAFILE, "<$fpath"))
            {
               $stmnt = "COPY arta_updatelist (sunum, loc, bytes) FROM stdin WITH DELIMITER '|'";
               ExecStatement(\$dbh, $stmnt, 1, "Troubles copying data to temporary table 'arta_updatelist'.\n");

               # Create a temporary data that holds all update rows.
               while (defined($line = <DATAFILE>))
               {
                  $rv = $dbh->pg_putcopydata($line);
                  if (!$rv)
                  {
                     print STDERR "Failure sending line '$line' to db.\n";
                     $err = 1;
                  }
               }
               
               if (!$err)
               {
                  $rv = $dbh->pg_putcopyend();
                  if (!$rv)
                  {
                     print STDERR "Failure ending table copy.\n";
                     $err = 1;
                  }
               }
               
               close(DATAFILE);
            }

            if (!$err)
            {
               if (0)
               {
                  # test to see if this worked.
                  $stmnt = "SELECT * from arta_updatelist";
                  $rrows = $dbh->selectall_arrayref($stmnt, undef);
                  $err = !(NoErr($rrows, \$dbh, $stmnt));
                  
                  foreach $row (@$rrows)
                  {
                     print "row $row->[0], $row->[1], $row->[2]\n";
                  }
               }
               
            # Do a join that will update sum_main with data from arta_updatelist

            # XXXXXXXXXXXXX - using a surrogate SELECT statement (instead of the real UPDATE one)
            # because we don't want to actually change the db.
            # TBD - currently, this uses a select statement to test out run times. The real query
            # will be an update statement.
            # REAL QUERY --> "UPDATE sum_main m SET archive_status = '2', arch_tape = ul.tapeid, arch_tape_fn = ul.fileid, arch_tape_date = 'date string' FROM arta_updatelist ul WHERE (m.ds_index = ul.sunum)"
               $stmnt = "SELECT m.archive_status, m.arch_tape, m.arch_tape_fn, m.arch_tape_date, ul.sunum, ul.loc FROM arta_updatelist ul JOIN (SELECT ds_index, archive_status, arch_tape, arch_tape_fn, arch_tape_date FROM sum_main where storage_group = $group) m ON (m.ds_index = ul.sunum)";

               $rrows = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrows, \$dbh, $stmnt));

               if (!$err)
               {
                  if (kDEBUG)
                  {
                     my($col);
                     my(@outd);
                     my($ridx);

                     foreach $row (@$rrows)
                     {
                        $ridx = 0;

                        while ($ridx < scalar(@$row))
                        {
                           $outd[$ridx] = (defined($row->[$ridx]) && length($row->[$ridx]) > 0) ? $row->[$ridx] : 0;
                           $ridx++;
                        }

                        print "row $outd[0], $outd[1], $outd[2], $outd[3], $outd[4], $outd[5]\n";
                     }
                  }
               }
            }
            
            # drop rows from temp table
            $stmnt = "DELETE from arta_updatelist";
            ExecStatement(\$dbh, $stmnt, 1, "Couldn't drop rows from temporary table 'arta_updatelist'.\n");
         } # loop over data files

         close(FILELIST);
      }

      $dbh->disconnect();
   }
}


# AL FINAL
exit($err);



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
   
   $$dbh = DBI->connect($dsn, $dbuser, ''); # will need to put pass in .pg_pass

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
