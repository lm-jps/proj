#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

use DBI;
use DBD::Pg;
use Time::localtime;
use Switch;

use constant kDEBUG => 1;

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

my(%delnow);
my(%delwi100d);
my(%dellater);
my(%archivepend);

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

         if (kDEBUG)
         {
            $group = 6;
         }

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
            $err = !SaveResults($rrowsb, $group, \%delnow);
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
               $err = !SaveResults($rrowsb, $group, \%delwi100d);
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
              $err = !SaveResults($rrowsb, $group, \%dellater);
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
              $err = !SaveResults($rrowsb, $group, \%archivepend);
            }
         }

         if (kDEBUG)
         {
            last;
         }
      } # loop over storage groups

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

sub SaveResults
{
   my($rrows) = $_[0]; # reference to array
   my($group) = $_[1]; # scalar
   my($container) = $_[2]; # reference to hash

   my($row);
   my($ok) = 1;

   foreach $row (@$rrows)
   {
      if (defined($container->{lc($row->[0])}))
      {
         print "WARNING: series $row->[0] appears in more than one tape group.\n";
         $container->{lc($row->[0])}->{$group} = $row->[1];
      } 
      else
      {
         $container->{lc($row->[0])} = {$group => $row->[1]};
      }
   }

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

use constant kTypeSortNumrcAsc => 1;
use constant kTypeSortAlphaAsc => 2;

sub CombineHashKeys
{
   my($typesort) = $_[0];
   my($out) = $_[1]; # reference
   my(@hashes) = @_[2..$#_]; # array of hash references
   my($ahash);
   my(@superduper);
   my(@sorted);
   my(%seen);
   my($elem);

   my($ok) = 1;

   foreach $ahash (@hashes)
   {
      if (defined($ahash))
      {
         push(@superduper, keys(%$ahash));
      }
   }

   # sort 
   switch ($typesort) 
   {
      case kTypeSortNumrcAsc
      {
         @sorted = sort {$a <=> $b} @superduper;
      }
      case kTypeSortAlphaAsc
      {
         @sorted = sort {$a cmp $b} @superduper;
      }
      else
      {
         print "Unsupported sort operation '$typesort'.\n";
         $ok = 0;
      }
   }

   # eliminate duplicates
   foreach $elem (@sorted)
   {
      push(@$out, $elem) unless $seen{$elem}++;
   }

   return $ok;
}

sub PrintResults
{
   # Each of the elements in each of these hash arrays is a reference to a hash array.
   # The parent hash array is keyed by series name. Each child hash array is 
   # keyed by group with byte count values. The parent hash arrays do not necessarily 
   # have the same set of series.
   my($delnow) = $_[0];
   my($delwi100d) = $_[1];
   my($dellater) = $_[2];
   my($archivepend) = $_[3];

   my(@serieslist);
   my(@grouplist);
   my($elem);
   my($series);
   my($group);
   my($line);

   my($dnow);
   my($d100);
   my($dlater);
   my($ap);

   my($tdnow);
   my($td100);
   my($tdlater);
   my($tap);

   my($ok);

   $tdnow = 0;
   $td100 = 0;
   $tdlater = 0;
   $tap = 0;

   $line = sprintf("%-48s%-8s%-24s%-24s%-24s%-24s", "series", "group", "DP Now (GB)", "DP < 100d (GB)", "DP > 100d (GB)", "AP (GB)");
   print "$line\n";

   if (CombineHashKeys(kTypeSortAlphaAsc, \@serieslist, $delnow, $delwi100d, $dellater, $archivepend))
   {
      foreach $elem (@serieslist)
      {
         $series = $elem;
         @grouplist = ();

         if (CombineHashKeys(kTypeSortNumrcAsc, \@grouplist, $delnow->{$series}, $delwi100d->{$series}, $dellater->{$series}, $archivepend->{$series}))
         {
            foreach $group (@grouplist)
            {
               $dnow = defined($delnow->{$series}->{$group}) ? $delnow->{$series}->{$group} : 0;
               $d100 = defined($delwi100d->{$series}->{$group}) ? $delwi100d->{$series}->{$group} : 0;
               $dlater = defined($dellater->{$series}->{$group}) ? $dellater->{$series}->{$group} : 0;
               $ap = defined($archivepend->{$series}->{$group}) ? $archivepend->{$series}->{$group} : 0;

               $tdnow += $dnow;
               $td100 += $d100;
               $tdlater += $dlater;
               $tap += $ap;
         
               $line = sprintf("%-48s%-8d%-24f%-24f%-24f%-24f", $series, $group, $dnow / kGig, $d100 / kGig, $dlater / kGig, $ap / kGig);
               print "$line\n";
            }
         }
         else
         {
            print "Problem creating group list - continuing.\n";
         }
      }
   }
   else
   {
      print "Problem creating series list - bailing.\n";
   }
   
   # TODO - print out totals by group
}
