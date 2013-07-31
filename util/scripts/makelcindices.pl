#!/home/jsoc/bin/linux_x86_64/activeperl -w

# makelcindices.pl jsoc hmidb 5432 postgres
# Must run with a user that is the owner of every table on which an index is to be created, so you
# must run as user postgres for the most part.

# This is a list of namespace tables, and the columns which require lower() indices on them.

# drms_series - seriesname 
# drms_segment - seriesname, segmentname
# drms_link - seriesname
# drms_keyword - seriesname, keywordname


use DBI;
use DBD::Pg;

use constant kIdxName => "<table>_<col>_lower";


my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)
my($dsn);
my($dbh);       # database handle

my($stmnt);
my($indexname);
my(@tables);
my($table);
my($nsp);
my($rrows);
my($row);

my(%tabandcols);

my($err);

$#ARGV == 3 || die "Improper argument list.\n";

$dbname = $ARGV[0];
$dbhost = $ARGV[1];
$dbport = $ARGV[2];
$dbuser = $ARGV[3];

$err = 0; 

# create hash array of tables and their columns that need lower-cases indices
$tabandcols{'drms_series'} = [ 'seriesname' ];
$tabandcols{'drms_segment'} = [ 'seriesname', 'segmentname' ];
$tabandcols{'drms_link'} = [ 'seriesname' ];
$tabandcols{'drms_keyword'} = [ 'seriesname', 'keywordname' ];

# connect to the database
$dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
print "Connection to database with '$dsn' as user '$dbuser' ... ";

# Despite ALL documentation saying otherwise, it looks like the error codes/string
# provided by DBI are all UNDEFINED, unless there is some kind of failure. So, 
# never try to look at $dbh->err or $dbh->errstr if the call succeeded.
$dbh = DBI->connect($dsn, $dbuser, ''); # will need to put pass in .pg_pass

if (defined($dbh))
{
   print "success!\n";
  
   # Collect all namespaces (from admin.ns)
   $stmnt = "SELECT name FROM admin.ns";
   $rrows = $dbh->selectall_arrayref($stmnt, undef);
   $err = !(NoErr($rrows, \$dbh, $stmnt));

   if (!$err)
   {
      # $rrows is a reference to an array; the array is an array of refereces to an array, so $row
      # is a reference to an array that has just one element (since the SELECT statement has just
      # one column). This element is the namespace name.
      foreach $row (@$rrows)
      {
         $nsp = $row->[0];

         # set the default namespace to the namespace in the current row
         $stmnt = "SET search_path TO $nsp";
         ExecStatement(\$dbh, $stmnt, 1, "Unable to change default namespace '$stmnt'.\n");

         # create the list of tables within the namespace on which to create indices
         @tables = keys(%tabandcols);
         foreach $table (@tables)
         {
            foreach $col (@{$tabandcols{$table}})
            {
               $indexname = kIdxName;
               $indexname =~ s/<table>/$table/;
               $indexname =~ s/<col>/$col/;

               # Create indices on tables only if the indices do not already exist
               $stmnt = "SELECT c.oid, n.nspname, c.relname FROM pg_catalog.pg_class c LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace WHERE c.relname ~ '^($indexname)\$' AND n.nspname ~ '^($nsp)\$' ORDER BY 2, 3";

               $rrows = $dbh->selectall_arrayref($stmnt, undef);
               $err = !NoErr($rrows, \$dbh, $stmnt);

               if (!$err)
               {
                  my(@rows) = @$rrows;

                  if ($#rows == 0)
                  {
                     # Index exists - continue
                     print "namespace $nsp: index $indexname exists, skipping\n";
                     next;
                  } 
                  else
                  {
                     # Index does NOT exist - create it
                     $stmnt = "CREATE INDEX $indexname on $table (lower($col))";
                     ExecStatement(\$dbh, $stmnt, 1, "Unable to create index with '$stmnt'.\n");
                  }
               }
            }
         }
      }
   }
   else
   {
      print "problem...the problem is you!\n";
   }

   $dbh->disconnect();
}
else
{
   print "failure!!!!\n";
}

# AL FINAL
exit($err);

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
