#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

use DBI;
use DBD::Pg;
use Time::localtime;

use constant kStatDAAP => "4"; # archive pending
use constant kSubStatDAADP => "128"; # after archive completes, mark delete pending
use constant kGig => 1073741824;
use constant kByteThreshold => 750; # in GB

my($err);

my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)

my($dsn);       # database connection arguments
my($dbh);       # database handle
my($stmnt);
my($rrowsa);
my($rrowsb);
my($rowa);
my($rowb);
my($rowc);

my($now);

$#ARGV == 3 || die "Improper argument list.\n";

$dbname = $ARGV[0];
$dbhost = $ARGV[1];
$dbport = $ARGV[2];
$dbuser = $ARGV[3];

$err = 0; 

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

   # Loop over storage groups
   $stmnt = "SELECT group_id FROM sum_partn_alloc GROUP BY group_id ORDER BY group_id";
   $rrowsa = $dbh->selectall_arrayref($stmnt, undef);
   $err = !(NoErr($rrowsa, \$dbh, $stmnt));

   if (!$err)
   {
      # $rrowsa is a reference to an array; the array is an array of refereces to an array, so $rowa
      # is a reference to an array that has just one element (since the SELECT statement has just
      # one column). This element is the namespace name.
      foreach $rowa (@$rrowsa)
      {
         $group = $rowa->[0];

         # archive pending
         $stmnt = "SELECT main.owning_series, partn.ds_index, partn.bytes, main.online_loc FROM (SELECT ds_index, group_id, bytes FROM sum_partn_alloc WHERE status = " . kStatDAAP . " AND archive_substatus = " . kSubStatDAADP . " AND group_id = $group) AS partn, (SELECT ds_index, owning_series, online_loc FROM sum_main WHERE storage_group = $group) AS main WHERE partn.ds_index = main.ds_index ORDER BY partn.ds_index";

         $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
         $err = !(NoErr($rrowsb, \$dbh, $stmnt));

         if (!$err)
         {
            my($tbytes);
            my($rseries);
            my($rsunum);
            my($rbytes);
            my($rpath);
            my(@rstore);
            my(@sortedrstore);
            my($sbytes);
            my($sseries);
            my($line);
            my($tapeno);
            my($fileno);

            $line = sprintf("***** Group %d *****", $group);
            print "$line\n";

            # check for no SUs marked as AP
            if ($#{$rrowsb} < 0)
            {
               print "     No storage units in this group are marked archive pending.\n";
            }

            # now, read from this table until we get 750 GB
            $tbytes = 0;
            $tapeno = 1;
            foreach $rowb (@$rrowsb)
            {
               $rseries = lc($rowb->[0]);
               $rsunum = $rowb->[1]; # this is a string
               $rbytes = $rowb->[2]; # this is a string
               $rpath = $rowb->[3];

               $rbytes = $rbytes / kGig; # convert to GB

               # store row 
               push(@rstore, [$rseries, $rsunum, $rbytes, $rpath]);

               $tbytes += $rbytes;
               
               if ($tbytes > kByteThreshold)
               {
                  $line = sprintf("     ***** Tape Number %d *****", $tapeno);
                  print "$line\n";

                  $line = sprintf("        %-16s%-48s%-16s", " file number", " series", " GB");
                  print "$line\n";

                  $line = sprintf("        %-16s%-48s%-16s", " -----------", " ------", " --");
                  print "$line\n";

                  # time to group by series and sum over bytes
                  @sortedrstore = sort {$a->[0] cmp $b->[0]} @rstore;

                  # now sum GB over series
                  $sbytes = 0;
                  $sseries = undef;
                  $fileno = 1;
                  foreach $rowc (@sortedrstore)
                  {
                     if (!defined($sseries))
                     {
                        $sseries = $rowc->[0];
                     }
                     
                     if ($sseries eq $rowc->[0])
                     {
                        $sbytes += $rowc->[2];
                     }
                     else
                     {
                        # next series - save existing row and start a new row
                        $line = sprintf("         %-16d%-48s%-16f", $fileno, $sseries, $sbytes);
                        print "$line\n";
                        $sbytes = $rowc->[2];
                        $sseries = $rowc->[0];
                        $fileno++;
                     }
                  }

                  $line = sprintf("        Total GB: %f", $tbytes);
                  print "$line\n";

                  $tbytes = 0;
                  @rstore = ();
                  $tapeno++;
               }
            }

            if ($tbytes)
            {
               # These didn't make it into a tape file because there were fewer than 750GB.
               $line = sprintf("     %f GB remaining", $tbytes);
               print "$line\n";
            }
         }
      }
   }
}
else
{
   print "failure!!!!\n";
   $err = 1;
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
