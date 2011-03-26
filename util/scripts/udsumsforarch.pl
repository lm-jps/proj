#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

# udsumsforarch.pl jsoc_sums hmidb 5434 production many

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

# In raw mode, you cannot collect data on dpshort, dpmid, dplong, ap all at the same time. 
# This will exhaust machine memory. Need to provide yet another flag to select which ONE
# of these metrics to obtain - only collect data on a single metric.

my($err);

my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)
my($typequery); 

my($dsn);       # database connection arguments
my($dbh);       # database handle
my($stmnt);
my($row);
my($rrows);


if ($#ARGV < 4)
{
   print "Improper argument list.\n"; 
   exit(1);
}

$dbname = $ARGV[0];
$dbhost = $ARGV[1];
$dbport = $ARGV[2];
$dbuser = $ARGV[3];
$typequery = $ARGV[4];


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
   
   if (dbconnect($dbname, $dbhost, $dbport, \$dbh))
   {
      $stmnt = "CREATE TEMPORARY TABLE arta_updatelist (sunum bigint, loc character varying(80), bytes bigint)";
      ExecStatement(\$dbh, $stmnt, 1, "Unable to create temporary table 'arta_updatelist'.\n");

      $stmnt = "CREATE INDEX artspg_idx on arta_updatelist (sunum)";   
      ExecStatement(\$dbh, $stmnt, 1, "Unable to create index on temporary table 'arta_updatelist'.\n");

      $stmnt = "COPY arta_updatelist (sunum, loc, bytes) FROM stdin WITH DELIMITER '|'";
      ExecStatement(\$dbh, $stmnt, 1, "Troubles copying data to temporary table 'arta_updatelist'.\n");

      push(@rowdata, "148755880|/SUM19/D148755880|17684\n");
      push(@rowdata, "148756590|/SUM4/D148756590|17684\n");
      push(@rowdata, "148757174|/SUM1/D148757174|17684\n");
      push(@rowdata, "148757826|/SUM5/D148757826|17684\n");
      push(@rowdata, "148758544|/SUM18/D148758544|17684\n");

      foreach $line (@rowdata)
      {
         $rv = $dbh->pg_putcopydata($line);
      }

      $rv = $dbh->pg_putcopyend();

      if (!$err)
      {
         if (kDEBUG)
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
      }
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
