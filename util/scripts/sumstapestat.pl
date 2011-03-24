#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

use DBI;
use DBD::Pg;
use Time::localtime;
use Switch;

use constant kDEBUG => 0;

use constant kStatDADP => "2";
use constant kStatDAAP => "4"; # archive pending
use constant kSubStatDAADP => "128"; # after archive completes, mark delete pending
use constant kGig => 1073741824;

use constant kTypeQueryAgg => "agg";
use constant kTypeQueryRaw => "raw";
use constant kTypeOrderSeries => "series";
use constant kTypeOrderGroup => "group";

my($err);

my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)
my($typequery); # type of query to perform (aggregate bytes over series, groups, or 
                # don't aggregate)
my($order);     # series - order by series, then group, or group, order by group, then series

my($dsn);       # database connection arguments
my($dbh);       # database handle
my($stmnt);
my($row);
my($rowb);
my($rrows);
my($rrowsb);

my($queryfunc);

my($now);
my($nowplus100d);

my(%delnow);
my(%delwi100d);
my(%dellater);
my(%archivepend);

if ($#ARGV < 3)
{
   print "Improper argument list.\n"; 
   exit(1);
}

$dbname = $ARGV[0];
$dbhost = $ARGV[1];
$dbport = $ARGV[2];
$dbuser = $ARGV[3];

if ($#ARGV >= 4)
{
   switch (lc($ARGV[4]))
   {
      case kTypeQueryAgg {$typequery = kTypeQueryAgg; $queryfunc = \&GenQueryA;}
      case kTypeQueryRaw {$typequery = kTypeQueryRaw; $queryfunc = \&GenQueryB;}
      else {print "Invalid query type $ARGV[4].\n"; exit(1);}
   }
}
else
{
   $typequery = kTypeQueryAgg;
   $queryfunc = \&GenQueryA;
   $order = kTypeOrderSeries;
}

if ($#ARGV >= 5)
{
   switch (lc($ARGV[5]))
   {
      case kTypeOrderSeries {$order = kTypeOrderSeries;}
      case kTypeOrderGroup {$order = kTypeOrderGroup;}
      else {print "Invalid order specified $ARGV[5].\n"; exit(1);}
   }
}
else
{
   $order = kTypeOrderSeries;
}

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
         $err = !$queryfunc->($group, kStatDADP, "effective_date < '$now'", \$stmnt);

         if (!$err)
         {
            $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
            $err = !(NoErr($rrowsb, \$dbh, $stmnt));
         }

         if (!$err)
         {
            # save results
            $err = !SaveResults($rrowsb, $group, $typequery, $order, \%delnow);
         }

         # Delete <= 100 days AND > now
         if (!$err)
         {
            $err = !$queryfunc->($group, kStatDADP, "effective_date <= '$nowplus100d' AND effective_date >= '$now'", \$stmnt);
            if (!$err)
            {
               $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrowsb, \$dbh, $stmnt));
            }

            if (!$err)
            {
               # save results
               $err = !SaveResults($rrowsb, $group, $typequery, $order, \%delwi100d);
            }
         }

         # Delete > 100 days;
         if (!$err)
         {
            $err = !$queryfunc->($group, kStatDADP, "effective_date > '$nowplus100d'", \$stmnt);
            if (!$err)
            {
               $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrowsb, \$dbh, $stmnt));
            }

            if (!$err)
            {
               # save results
              $err = !SaveResults($rrowsb, $group, $typequery, $order, \%dellater);
            }
         }

         # Archive Pending
         if (!$err)
         {
            $err = !$queryfunc->($group, kStatDAAP . " AND archive_substatus = " . kSubStatDAADP, "", \$stmnt);

            if (!$err)
            {
               $rrowsb = $dbh->selectall_arrayref($stmnt, undef);
               $err = !(NoErr($rrowsb, \$dbh, $stmnt));
            }

            if (!$err)
            {
               # save results 
              $err = !SaveResults($rrowsb, $group, $typequery, $order, \%archivepend);
            }
         }

         if (kDEBUG)
         {
            last;
         }
      } # loop over storage groups

      SortAndPrintResults(\%delnow, \%delwi100d, \%dellater, \%archivepend, $typequery, $order);
   }
}
else
{
   print "failure!!!!\n";
   $err = 1;
}

# AL FINAL
exit($err);

# Aggregate
sub GenQueryA
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

   $$qout = "SELECT main.owning_series, sum(bytes) FROM (SELECT ds_index, group_id, bytes FROM sum_partn_alloc WHERE status = $status AND group_id = $group$datewhere) AS partn, (SELECT ds_index, owning_series FROM sum_main WHERE storage_group = $group) AS main WHERE partn.ds_index = main.ds_index GROUP BY partn.group_id, main.owning_series ORDER BY lower(main.owning_series)";   

   if (kDEBUG)
   {
      print "Query is:\n$$qout\n";
   }

   return $ok;
}

# Don't aggregate
sub GenQueryB
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

   $$qout = "SELECT main.owning_series, main.ds_index, main.online_loc, partn.bytes FROM (SELECT ds_index, group_id, bytes FROM sum_partn_alloc WHERE status = $status AND group_id = $group$datewhere) AS partn, (SELECT ds_index, owning_series, online_loc FROM sum_main WHERE storage_group = $group) AS main WHERE partn.ds_index = main.ds_index ORDER BY lower(main.owning_series), main.ds_index";

   if (kDEBUG)
   {
      print "Query is:\n$$qout\n";
   }

   return $ok;
}

sub SaveResults
{
   my($rrows) = $_[0]; # reference to array
   my($group) = $_[1]; # scalar
   my($typequery) = $_[2]; # scalar
   my($order) = $_[3]; # scalar
   my($container) = $_[4]; # reference to hash

   my($row);
   my($ok) = 1;

   switch ($typequery)
   {
      case kTypeQueryAgg
      {
         # row is series, sum(bytes)
         switch ($order)
         {
            case kTypeOrderSeries
            {
               foreach $row (@$rrows)
               {
                  if (defined($container->{lc($row->[0])}))
                  {
                     $container->{lc($row->[0])}->{$group} = $row->[1];
                  } 
                  else
                  {
                     $container->{lc($row->[0])} = {$group => $row->[1]};
                  }
               }
            }
            case kTypeOrderGroup
            {
               $container->{$group} = [];

               foreach $row (@$rrows)
               {
                  push(@{$container->{$group}}, [lc($row->[0]), $row->[1]]);
               }
            }
            else
            {
               print "Invalid column $order by which to order.\n";
               $ok = 0;
            }
         } # switch $order
      }
      case kTypeQueryRaw
      {
         # row is series, ds_index, sudir, bytes
         switch ($order)
         {
            case kTypeOrderSeries
            {
               foreach $row (@$rrows)
               {
                  if (defined($container->{lc($row->[0])}))
                  {
                     push(@{$container->{lc($row->[0])}->{$group}}, [$row->[1], $row->[2], $row->[3]])
                  } 
                  else
                  {
                     $container->{lc($row->[0])} = {$group => [[$row->[1], $row->[2], $row->[3]]]};
                  }
               }
            }
            case kTypeOrderGroup
            {
               $container->{$group} = [];

               foreach $row (@$rrows)
               {
                  push(@{$container->{$group}}, [lc($row->[0]), $row->[1], $row->[2], $row->[3]]);
               }
            }
            else
            {
               print "Invalid column $order by which to order.\n";
               $ok = 0;
            }
         }
      }
      else
      {
         print "Invalid query type $typequery.\n";
         $ok = 0;
      }
   } # switch type query

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

# Orders by series, group first. If caller requests ordering by group, series, the
# rows are re-ordered.
sub SortAndPrintResults
{
   # Each of the elements in each of these hash arrays is a reference to a hash array.
   # The parent hash array is keyed by series name. Each child hash array is 
   # keyed by group with byte count values. The parent hash arrays do not necessarily 
   # have the same set of series.
   my($delnow) = $_[0];
   my($delwi100d) = $_[1];
   my($dellater) = $_[2];
   my($archivepend) = $_[3];
   my($typequery) = $_[4];
   my($order) = $_[5];

   my(@serieslist);
   my(@grouplist);
   my($elem);
   my($series);
   my($group);
   my($line);

   my(@sorted);
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


   switch ($typequery)
   {
      case kTypeQueryAgg
      {
         switch ($order)
         {
            case kTypeOrderSeries
            {
               # type - agg; order - series
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
            }
            case kTypeOrderGroup
            {
               # type - agg; order - group
               my(@topdataelem); # 4 elements; each points to the top data element of one container
               my(@topdataidx); # indices into the top elements of the 4 containers
               my(@stack); # index of first (in sort order) data element when each of the top data elements are sorted
               my(@noelems); # number of elements in each of the 4 containers
               my($popped); 
               my($contidx);
               my(@bytes); # sum(bytes) of the 4 containers for current series.

               $line = sprintf("%-8s%-48s%-24s%-24s%-24s%-24s", "group", "series", "DP Now (GB)", "DP < 100d (GB)", "DP > 100d (GB)", "AP (GB)");
               print "$line\n";

               if (CombineHashKeys(kTypeSortNumrcAsc, \@grouplist, $delnow, $delwi100d, $dellater, $archivepend))
               {
                  foreach $group (@grouplist)
                  {
                     # $container->{$group} is a reference to an array where each element is a reference to an array
                     # of data elements (series, sum(bytes)), ordered by series. Not every container
                     # has every series!

                     # Within each container, the arrays of data elements are sorted by series. So we 
                     # need to compare the first-row data elements from each container, finding
                     # the container that has the series name that is first in alphabetical order. We then
                     # print a record containing the sum(bytes) from each container that has a data element
                     # with that series name. If a container does not have such a data element, then 
                     # we print a value of 0 for the sum(bytes) field.
                     push(@noelems, scalar(@{$delnow->{$group}}));
                     push(@noelems, scalar(@{$delwi100d->{$group}}));
                     push(@noelems, scalar(@{$dellater->{$group}}));
                     push(@noelems, scalar(@{$archivepend->{$group}}));

                     @topdataidx = (0, 0, 0, 0); # indices into 4 containers
                     
                     # iterate through each top data element, looking for the one with a series that
                     # is first when sorted. $elem is a reference to a data element.
                     while (1)
                     {
                        @stack = (); # each elem contains index of container (0 - 3)
                        @bytes = ();

                        $topdataelem[0] = $topdataidx[0] < $noelems[0] ? $delnow->{$group}->[$topdataidx[0]] : [];
                        $topdataelem[1] = $topdataidx[1] < $noelems[1] ? $delwi100d->{$group}->[$topdataidx[1]] : [];
                        $topdataelem[2] = $topdataidx[2] < $noelems[2] ? $dellater->{$group}->[$topdataidx[2]] : [];
                        $topdataelem[3] = $topdataidx[3] < $noelems[3] ? $archivepend->{$group}->[$topdataidx[3]] : [];
                        $contidx = 0;

                        foreach $elem (@topdataelem)
                        {
                           if (kDEBUG)
                           {
                              if (!defined($elem))
                              {
                                 # for some reason, noelems was 29 for delnow, even though the db returned only 25 rows for group 6.
                                 print "elem not defined; contidx = $contidx, group = $group, topdataidx = $topdataidx[$contidx], noelems $noelems[$contidx]\n";
                                 exit;
                              }
                           }

                           if (scalar(@$elem) > 0)
                           {
                              $series = $elem->[0];

                              if ($#stack >= 0)
                              {
                                 $stackseries = $topdataelem[$stack[$#stack]]->[0];

                                 if (($series cmp $stackseries) < 0)
                                 {
                                    # $series sorts before $stack[0] - pop from stack
                                    while ($#stack >= 0)
                                    {
                                       $popped = pop(@stack);
                                       
                                       # This element came from a container that doesn't have $series.
                                       $bytes[$popped] = 0;
                                    }
                                    
                                    # Push $elem onto stack now (it is the current highest-ranked data element)
                                    push(@stack, $contidx);
                                 }
                                 elsif (($series cmp $stackseries) == 0)
                                 {
                                    # $series is equivalent to $stack[0] - push onto stack
                                    push(@stack, $contidx);
                                 } 
                                 else
                                 {
                                    # $series sorts after $stack[0] - set bytes to zero for this $series
                                    $bytes[$contidx] = 0;
                                 }
                              } 
                              else
                              {
                                 # stack is empty, push
                                 push(@stack, $contidx);
                              }
                           }
                           else
                           {
                              # current container is empty - set bytes to 0
                              $bytes[$contidx] = 0;
                           }
                           
                           $contidx++;
                        }

                        if ($#stack < 0)
                        {
                           # nothing was done in sort loop (because there were no more elems) done
                           last;
                        }

                        # Elems on the stack will have non-0 sum(byte) values.
                        # The elems on the stack all have the same series.
                        foreach $contidx (@stack)
                        {
                           $bytes[$contidx] = $topdataelem[$contidx]->[1];
                           $series = $topdataelem[$contidx]->[0];

                           # advance index into data elem array
                           $topdataidx[$contidx]++;
                        }

                        $dnow = $bytes[0];
                        $d100 = $bytes[1];
                        $dlater = $bytes[2];
                        $ap = $bytes[3];

                        $line = sprintf("%-8d%-48s%-24f%-24f%-24f%-24f", $group, $series, $dnow / kGig, $d100 / kGig, $dlater / kGig, $ap / kGig);
                        print "$line\n";
                     } # while
                  } # groups
               }
            }
            else
            {
               print "Invalid column $order by which to order.\n";
            }
         } # switch $order
      } # case agg
      case kTypeQueryRaw
      {
         # row is series, ds_index, sudir, bytes
         switch ($order)
         {
            case kTypeOrderSeries
            {
               # type - raw; order - series
               my($sunum);
               my($sudir);
               my($rowdata);
               
               # First, DP Now
               print "***DP NOW***\n";
               $line = sprintf("%-24s%-8s%-16s%-24s%-24s", "series", "group", "sunum", "sudir", "DP Now (bytes)");
               print "$line\n";

               if (CombineHashKeys(kTypeSortAlphaAsc, \@serieslist, $delnow))
               {
                  foreach $elem (@serieslist)
                  {
                     $series = $elem;
                     @grouplist = ();
                     
                     if (CombineHashKeys(kTypeSortNumrcAsc, \@grouplist, $delnow->{$series}))
                     {
                        foreach $group (@grouplist)
                        {
                           if (defined($delnow->{$series}->{$group}))
                           {
                              # $delnow->{$series}->{$group} is reference to an array with array references
                              # as elements. The child arrays have sunum, sudir, and bytes
                              # as elements.
                              foreach $rowdata (@{$delnow->{$series}->{$group}})
                              {
                                 # $rowdata is a reference to an array containing sunum, sudir, and bytes
                                 # as elements.
                                 $sunum = $rowdata->[0];
                                 $sudir = $rowdata->[1];
                                 $dnow = $rowdata->[2];

                                 $line = sprintf("%-24s%-8d%-16s%-24s%-24d", $series, $group, $sunum, $sudir, $dnow);
                                 print "$line\n";
                              }
                           }
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
            }
            case kTypeOrderGroup
            {
               
            }
            else
            {
               print "Invalid column $order by which to order.\n";
               $ok = 0;
            }
         }
      }
      else
      {
         print "Invalid query type $typequery.\n";

      }
   } # switch query type
   
   # TODO - print out totals by group
}
