#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

use DBI;
use DBD::Pg;
use Time::localtime;

use constant kStatDADP => "2";
use constant kStatDAAP => "4"; # archive pending
use constant kSubStatDAADP => "128"; # after archive completes, mark delete pending
use constant kGig => 1073741824;

my($err);

my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)

my($dsn);       # database connection arguments
my($dbh);       # database handle
my($stmnt);
my($row);
my($rowb);
my($rrows);
my($rrowsb);

my($now);
my($nowplus100d);

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
   $rrows = $dbh->selectall_arrayref($stmnt, undef);
   $err = !(NoErr($rrows, \$dbh, $stmnt));

   if (!$err)
   {
      my($timenow) = time();
      my($timeo) = localtime($timenow);
      my($timeltro) = localtime($timenow + 100 * 24 * 60 * 60);

      $now = sprintf("%04d%02d%02d%02d%02d", $timeo->year() + 1900, $timeo->mon() + 1, $timeo->mday(), $timeo->hour(), $timeo->min());
      $nowplus100d = sprintf("%04d%02d%02d%02d%02d", $timeltro->year() + 1900, $timeltro->mon() + 1, $timeltro->mday(), $timeltro->hour(), $timeltro->min());

      # $rrows is a reference to an array; the array is an array of refereces to an array, so $row
      # is a reference to an array that has just one element (since the SELECT statement has just
      # one column). This element is the namespace name.
      foreach $row (@$rrows)
      {
         $group = $row->[0];

         # Collect data for the group

         # Delete now
         $err = !GenQuery($group, kStatDADP, "effective_date < '$now'", \$stmnt);

         if (!$err)
         {
            $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
            $err = !(NoErr($rrowsb, \$dbh, $stmnt));
         }

         if (!$err)
         {
            # save results
            foreach $rowb (@$rrowsb)
            {
               $delnow{lc($rowb->[0])} = [$group, $rowb->[1]];
            }
         }

         # Delete <= 100 days AND > now
         if (!$err)
         {
            $err = !GenQuery($group, kStatDADP, "effective_date <= '$nowplus100d' AND effective_date >= '$now'", \$stmnt);
            if (!$err)
            {
               $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrowsb, \$dbh, $stmnt));
            }

            if (!$err)
            {
               # save results
               foreach $rowb (@$rrowsb)
               {
                  $delwi100d{lc($rowb->[0])} = [$group, $rowb->[1]];
               }
            }
         }

         # Delete > 100 days;
         if (!$err)
         {
            $err = !GenQuery($group, kStatDADP, "effective_date > '$nowplus100d'", \$stmnt);
            if (!$err)
            {
               $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrowsb, \$dbh, $stmnt));
            }

            if (!$err)
            {
               # save results
               foreach $rowb (@$rrowsb)
               {
                  $dellater{lc($rowb->[0])} = [$group, $rowb->[1]];
               }
            }
         }

         # Archive Pending
         if (!$err)
         {
            $err = !GenQuery($group, kStatDAAP . " AND archive_substatus = " . kSubStatDAADP, "", \$stmnt);

            if (!$err)
            {
               $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrowsb, \$dbh, $stmnt));
            }

            if (!$err)
            {
               # save results
               foreach $rowb (@$rrowsb)
               {
                  $archivepend{lc($rowb->[0])} = [$group, $rowb->[1]];
               }
            }
         }
      } # loop over storage groups

      # print
      PrintResults(\%delnow, \%delwi100d, \%dellater, \%archivepend);
   }
}
else
{
   print "failure!!!!\n";
   $err = 1;
}

# AL FINAL
exit($err);


sub GenQuery
{
   my($group) = $_[0];
   my($status) = $_[1];
   my($datewhere) = $_[2];
   my($qout) = $_[3];

   my($ok) = 1;

   if (length($datewhere) > 0)
   {
      $datewhere = " AND $datewhere";
   }

   $$qout = "SELECT main.owning_series, sum(bytes) FROM (SELECT ds_index, group_id, bytes FROM sum_partn_alloc WHERE status = $status AND group_id = $group$datewhere) AS partn, (SELECT ds_index, owning_series FROM sum_main WHERE storage_group = $group) AS main WHERE partn.ds_index = main.ds_index GROUP BY partn.group_id, main.owning_series ORDER BY main.owning_series";   

   return $ok;
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

sub PrintResults
{
   # Each of the elements in each of these hash arrays has a series and a byte count; 
   # the hash arrays do not necessarily all have the same set of series.
   my($delnow) = $_[0];
   my($delwi100d) = $_[1];
   my($dellater) = $_[2];
   my($archivepend) = $_[3];

   my(@kdelnow) = map {[$_, $delnow->{$_}->[0]]} keys(%$delnow); # series, group
   my(@kdelwi100d) = map {[$_, $delwi100d->{$_}->[0]]} keys(%$delwi100d);
   my(@kdellater) = map {[$_, $dellater->{$_}->[0]]} keys(%$dellater);
   my(@karchivepend) = map {[$_, $archivepend->{$_}->[0]]} keys(%$archivepend);
   
   my(@superduper) = (@kdelnow, @kdelwi100d, @kdellater, @karchivepend);

   my(@sdsort) = sort {$a->[0] cmp $b->[0]} @superduper;
   my(@serieslist);
   my($elem);
   my($series);
   my($group);
   my(%seen);
   my($line);

   my($dnow);
   my($d100);
   my($dlater);
   my($ap);

   foreach $elem (@sdsort)
   {
      $series =  $elem->[0];
      push(@serieslist, $elem) unless $seen{$series}++;
   }

   $line = sprintf("%-48s%-8s%-24s%-24s%-24s%-24s", "series", "group", "DP Now (GB)", "DP < 100d (GB)", "DP > 100d (GB)", "AP (GB)");
   print "$line\n";

   foreach $elem (@serieslist)
   {
      $series = $elem->[0];
      $group = $elem->[1];
      $dnow = defined($delnow->{$series}->[1]) ? $delnow->{$series}->[1] : 0;
      $d100 = defined($delwi100d->{$series}->[1]) ? $delwi100d->{$series}->[1] : 0;
      $dlater = defined($dellater->{$series}->[1]) ? $dellater->{$series}->[1] : 0;
      $ap = defined($archivepend->{$series}->[1]) ? $archivepend->{$series}->[1] : 0;

      $line = sprintf("%-48s%-8d%-24f%-24f%-24f%-24f", $series, $group, $dnow / kGig, $d100 / kGig, $dlater / kGig, $ap / kGig);
      print "$line\n";
   }
}
