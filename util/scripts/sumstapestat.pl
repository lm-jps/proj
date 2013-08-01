#!/home/jsoc/bin/linux_x86_64/activeperl -w

use DBI;
use DBD::Pg;
use Time::localtime;

use constant kDEBUG => 0;

use constant kStatDADP => "2";
use constant kStatDAAP => "4"; # archive pending
use constant kStatDAAEDDP => "32"; # archive pending too I guess
use constant kSubStatDAADP => "128"; # after archive completes, mark delete pending
use constant kGig => 1073741824;

use constant kTypeQueryAgg => "agg";
use constant kTypeQueryRaw => "raw";
use constant kTypeOrderSeries => "series";
use constant kTypeOrderGroup => "group";

use constant kMetricAll => "all";
use constant kMetricDPS => "dps"; # delete pending now (short)
use constant kMetricDPM => "dpm"; # delete pending in 100 days (medium)
use constant kMetricDPL => "dpl"; # delete pending >= 100 days (long)
use constant kMetricAP => "ap";   # archive pending

use constant kTempTable => "sumsstat";

# In raw mode, you cannot collect data on dpshort, dpmid, dplong, ap all at the same time. 
# This will exhaust machine memory. Need to provide yet another flag to select which ONE
# of these metrics to obtain - only collect data on a single metric.

my($err);

my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to connect with)
my($typequery); # type of query to perform (aggregate bytes over series, groups, or 
                # don't aggregate)
my($order);     # series - order by series, then group, or group, order by group, then series
my($metric);    # while column of data to produce (delete pending short, delete pending medium, delete pending long,
                # archive pending)
my($delim);     # separator between output columns (fixed-width columns if not supplied)

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
    if (lc($ARGV[4]) eq kTypeQueryAgg)
    {
        $typequery = kTypeQueryAgg; $queryfunc = \&GenQueryA;
    }
    elsif (lc($ARGV[4]) eq kTypeQueryRaw)
    {
        $typequery = kTypeQueryRaw; $queryfunc = \&GenQueryB;
    }
    else
    {
        print "Invalid query type $ARGV[4].\n"; exit(1);
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
    if (lc($ARGV[5]) eq kTypeOrderSeries)
    {
        $order = kTypeOrderSeries;
    }
    elsif (lc($ARGV[5]) eq kTypeOrderGroup)
    {
        $order = kTypeOrderGroup;
    }
    else
    {
        print "Invalid order specified $ARGV[5].\n"; exit(1);
    }
}
else
{
   $order = kTypeOrderSeries;
}

if ($#ARGV >= 6)
{
    if (lc($ARGV[6]) eq kMetricAll)
    {
        if ($typequery eq kTypeQueryRaw)
        {
            print "Metric all cannot be used with an un-aggregated query.\n";
            exit(1);
        }
        
        $metric = kMetricAll;
    }
    elsif (lc($ARGV[6]) eq kMetricDPS)
    {
        $metric = kMetricDPS;
    }
    elsif (lc($ARGV[6]) eq kMetricDPM)
    {
        $metric = kMetricDPM;
    }
    elsif (lc($ARGV[6]) eq kMetricDPL)
    {
        $metric = kMetricDPL;
    }
    elsif (lc($ARGV[6]) eq kMetricAP)
    {
        $metric = kMetricAP;
    }
    else
    {
        print "Invalid metric specified $ARGV[6].\n"; exit(1);
    }
}
else
{
   if ($typequery eq kTypeQueryRaw)
   {
      # For non-aggregated queries, do not attempt to perform multiple queries for multiple metrics
      $metric = kMetricAP;
   }
   else
   {
      $metric = kMetricAll;
   }
}

if ($#ARGV >= 7)
{
   $delim = substr($ARGV[7], 0, 1);
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


      if ($typequery)
      {
         # Create a temporary table to hold pre-sorted results.
         $stmnt = "CREATE TEMPORARY TABLE " . kTempTable . "(tgroup integer not null, series varchar(64) not null, metric varchar(8), aggbytes bigint default 0)";
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create temporary table " . "'kTempTable'" . ".\n");
         
         $stmnt = "CREATE INDEX " . kTempTable . "_group_idx on " . kTempTable . "(tgroup, lower(series))";   
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create index on temporary table " . "'kTempTable'" . ".\n");
         
         $stmnt = "CREATE INDEX " . kTempTable . "_series_idx on " . kTempTable . "(lower(series), tgroup)";
         ExecStatement(\$dbh, $stmnt, 1, "Unable to create index on temporary table " . "'kTempTable'" . ".\n");
      }

      # $rrows is a reference to an array; the array is an array of refereces to an array, so $row
      # is a reference to an array that has just one element (since the SELECT statement has just
      # one column). This element is the namespace name.
      foreach $row (@$rrows)
      {
         $group = $row->[0];

         if (kDEBUG)
         {
            $group = 1;
         }

         # Delete now
         if ($metric eq kMetricAll || $metric eq kMetricDPS)
         {
            $err = !$queryfunc->(\$dbh, $group, kStatDAAP . " AND archive_substatus = " . kStatDAAEDDP . " OR status = " . kStatDADP, "effective_date < '$now'", $typequery, $order, \%delnow, kMetricDPS);
         }

         # Delete <= 100 days AND > now
         if ($metric eq kMetricAll || $metric eq kMetricDPM)
         {
            if (!$err)
            {
               $err = !$queryfunc->(\$dbh, $group, kStatDAAP . " AND archive_substatus = " . kStatDAAEDDP . " OR status = " . kStatDADP, "effective_date <= '$nowplus100d' AND effective_date >= '$now'", $typequery, $order, \%delwi100d, kMetricDPM);
            }
         }

         if ($metric eq kMetricAll || $metric eq kMetricDPL)
         {
            if (!$err)
            {
               $err = !$queryfunc->(\$dbh, $group, kStatDAAP . " AND archive_substatus = " . kStatDAAEDDP . " OR status = " . kStatDADP, "effective_date > '$nowplus100d'", $typequery, $order, \%dellater, kMetricDPL);
            }
         }

         # Archive Pending
         if ($metric eq kMetricAll || $metric eq kMetricAP)
         {
            if (!$err)
            {
               $err = !$queryfunc->(\$dbh, $group, kStatDAAP . " AND archive_substatus = " . kSubStatDAADP, "", $typequery, $order, \%archivepend, kMetricAP);
            }
         }

         if (kDEBUG)
         {
            last;
         }
      } # loop over storage groups

      SortAndPrintResults(\%delnow, \%delwi100d, \%dellater, \%archivepend, $typequery, $order, $metric, $delim, $dbh);
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
   my($dbh) = $_[0]; # reference to a reference
   my($group) = $_[1];
   my($status) = $_[2];
   my($datewhere) = $_[3];
   my($typequery) = $_[4];
   my($order) = $_[5];
   my($container) = $_[6]; # reference
   my($metric) = $_[7];

   my($ok) = 1;

   my($stmnt);

   if (length($datewhere) > 0)
   {
      $datewhere = " AND $datewhere";
   }

   $stmnt = "INSERT INTO " . kTempTable . " (tgroup, series, metric, aggbytes) (SELECT '$group', main.owning_series, '$metric', sum(bytes) FROM (SELECT ds_index, group_id, bytes FROM sum_partn_alloc WHERE (status = $status) AND group_id = $group$datewhere) AS partn, (SELECT ds_index, owning_series FROM sum_main) AS main WHERE partn.ds_index = main.ds_index GROUP BY partn.group_id, main.owning_series)";   

   if (kDEBUG)
   {
      print "Query is:\n$stmnt\n";
   }

   ExecStatement($dbh, $stmnt, 1, "Unable to insert data into temporary table.\n");

   return $ok;
}

# Don't aggregate
sub GenQueryB
{
   my($dbh) = $_[0]; # reference to a reference
   my($group) = $_[1];
   my($status) = $_[2];
   my($datewhere) = $_[3];
   my($typequery) = $_[4];
   my($order) = $_[5];
   my($container) = $_[6]; # reference
   my($metric) = $_[7];

   my($ok) = 1;

   my($stmnt);
   my($rrows);

   if (length($datewhere) > 0)
   {
      $datewhere = " AND $datewhere";
   }

   $stmnt = "SELECT main.owning_series, main.ds_index, main.online_loc, partn.bytes FROM (SELECT ds_index, group_id, bytes FROM sum_partn_alloc WHERE (status = $status) AND group_id = $group$datewhere) AS partn, (SELECT ds_index, owning_series, online_loc FROM sum_main) AS main WHERE partn.ds_index = main.ds_index ORDER BY lower(main.owning_series), main.ds_index";

   if (kDEBUG)
   {
      print "Query is:\n$stmnt\n";
   }

   $rrows = $$dbh->selectall_arrayref($stmnt, undef);
   $ok = NoErr($rrows, $$dbh, $stmnt);

   if ($ok)
   {
      # save results
      $ok = SaveResults($rrows, $group, $typequery, $order, $container);
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
    
    if ($typequery eq kTypeQueryAgg)
    {
        # Changed to use a temporary table to hold results. No need to use hash arrays
        # to hold the data.
    }
    elsif ($typequery eq kTypeQueryRaw)
    {
        if ($order eq kTypeOrderSeries)
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
        elsif ($order eq kTypeOrderGroup)
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
    else
    {
        print "Invalid query type $typequery.\n";
        $ok = 0;
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
    if ($typesort eq kTypeSortNumrcAsc)
    {
        @sorted = sort {$a <=> $b} @superduper;
    }
    elsif ($typesort eq kTypeSortAlphaAsc)
    {
        @sorted = sort {$a cmp $b} @superduper;
    }
    else
    {
        print "Unsupported sort operation '$typesort'.\n";
        $ok = 0;
    }
    
    # eliminate duplicates
    foreach $elem (@sorted)
    {
        push(@$out, $elem) unless $seen{$elem}++;
    }
    
    return $ok;
}

sub PrintRow
{
   my($delim) = $_[0];
   my($firstarg) = $_[1];
   my($secondarg) = $_[2];
   my($metric) = $_[3];
   my($hbytes) = $_[4]; # reference
   my($dformat) = $_[5]; # reference to delimted string format
   my($fformat) = $_[6]; # reference to fixed-width string format

   my($line);
   my($actualformat);

   $actualformat = defined($delim) ? $dformat : $fformat;

   if ($metric eq kMetricAll)
   {
      $line = sprintf($$actualformat, 
                      $firstarg,
                      $secondarg,
                      defined($hbytes->{+kMetricDPS}) ? $hbytes->{+kMetricDPS} / kGig : 0, 
                      defined($hbytes->{+kMetricDPM}) ? $hbytes->{+kMetricDPM} / kGig : 0,
                      defined($hbytes->{+kMetricDPL}) ? $hbytes->{+kMetricDPL} / kGig : 0,
                      defined($hbytes->{+kMetricAP}) ? $hbytes->{+kMetricAP} / kGig : 0);
   }
   else
   {
      $line = sprintf($$actualformat, 
                      $firstarg,
                      $secondarg,
                      defined($hbytes->{$metric}) ? $hbytes->{$metric} / kGig : 0);
   }

   print "$line\n";
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
    my($metric) = $_[6];
    my($delim) = $_[7];
    my($dbh) = $_[8];
    
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
    my($metricval);
    
    my(%metricheaders);
    my(%containers);
    my(@contkeys);
    my($contkey);
    
    my($stmnt);
    my($rrows);
    my($row);
    
    my($ok);
    
    $ok = 1;
    
    %containers = (kMetricDPS, $delnow, kMetricDPM, $delwi100d, kMetricDPL, $dellater, kMetricAP, $archivepend);
    @contkeys = keys(%containers);
    
    print "__DATA__\n";
    
    if ($typequery eq kTypeQueryAgg)
    {
        my(%hbytes); # sum(bytes) of the 4 containers for current series.
        
        %metricheaders = (kMetricDPS, "DP Now (GB)", kMetricDPM, "DP <= 100d (GB)", kMetricDPL, "DP > 100d (GB)", kMetricAP, "AP (GB)");
        
        if ($order eq kTypeOrderSeries)
        {
            # type - agg; order - series
            if (defined($delim))
            {
                if ($metric eq kMetricAll)
                {
                    $line = sprintf("series${delim}group${delim}$metricheaders{+kMetricDPS}${delim}$metricheaders{+kMetricDPM}${delim}$metricheaders{+kMetricDPL}${delim}$metricheaders{+kMetricAP}");
                }
                else
                {
                    $line = sprintf("series${delim}group${delim}$metricheaders{$metric}");
                }
            }
            else
            {
                if ($metric eq kMetricAll)
                {
                    $line = sprintf("%-48s%-8s%-24s%-24s%-24s%-24s", "series", "group", $metricheaders{+kMetricDPS}, $metricheaders{+kMetricDPM}, $metricheaders{+kMetricDPL}, $metricheaders{+kMetricAP});
                }
                else
                {
                    $line = sprintf("%-48s%-8s%-24s", "series", "group", $metricheaders{$metric});
                }
            }
            
            print "$line\n";
            
            # Just use the db to do the sorting on the temporary table containing the data.
            $stmnt = "SELECT series, tgroup, metric, aggbytes FROM " . kTempTable . " ORDER BY lower(series), tgroup";
            
            if ($ok)
            {
                $rrows = $dbh->selectall_arrayref($stmnt, undef);
                $ok = NoErr($rrows, \$dbh, $stmnt);
            }
            
            $group = "";
            $series = "";
            
            if ($ok)
            {
                %hbytes = ();
                
                if ($metric eq kMetricAll)
                {
                    if (defined($delim))
                    {
                        $dformat = "%s${delim}%s${delim}%f${delim}%f${delim}%f${delim}%f";
                    } 
                    else
                    {
                        $fformat = "%-48s%-8d%-24f%-24f%-24f%-24f";
                    }
                } 
                else
                {
                    if (defined($delim))
                    {
                        $dformat = "%s${delim}%s${delim}%f";
                    } 
                    else
                    {
                        $fformat = "%-48s%-8d%-24f";
                    }
                }
                
                foreach $row (@$rrows)
                {
                    if ((length($series) > 0 && $series ne $row->[0]) ||
                        (length($group) > 0 && $group ne $row->[1]))
                    {
                        PrintRow($delim, $series, $group, $metric, \%hbytes, \$dformat, \$fformat);
                        %hbytes = ();
                    }
                    
                    $series = $row->[0];
                    $group = $row->[1];
                    $hbytes{$row->[2]} = $row->[3];
                }
                
                # Must print last row since only the previous row is printed in the loop above
                PrintRow($delim, $series, $group, $metric, \%hbytes, \$dformat, \$fformat);
            }
        }
        elsif ($order eq kTypeOrderGroup)
        {
            # type - agg; order - group
            if (defined($delim))
            {
                if ($metric eq kMetricAll)
                {
                    $line = sprintf("group${delim}series${delim}$metricheaders{+kMetricDPS}${delim}$metricheaders{+kMetricDPM}${delim}$metricheaders{+kMetricDPL}${delim}$metricheaders{+kMetricAP}");
                }
                else
                {
                    $line = sprintf("group${delim}series${delim}$metricheaders{$metric}");
                }
            }
            else
            {
                if ($metric eq kMetricAll)
                {
                    $line = sprintf("%-8s%-48s%-24s%-24s%-24s%-24s", "group", "series", $metricheaders{+kMetricDPS}, $metricheaders{+kMetricDPM}, $metricheaders{+kMetricDPL}, $metricheaders{+kMetricAP});
                }
                else
                {
                    $line = sprintf("%-8s%-48s%-24s", "group", "series", $metricheaders{$metric});
                }
            }
            
            print "$line\n";
            
            # Just use the db to do the sorting on the temporary table containing the data.
            $stmnt = "SELECT tgroup, series, metric, aggbytes FROM " . kTempTable . " ORDER BY tgroup, lower(series)";
            
            if ($ok)
            {
                $rrows = $dbh->selectall_arrayref($stmnt, undef);
                $ok = NoErr($rrows, \$dbh, $stmnt);
            }
            
            $group = "";
            $series = "";
            
            if ($ok)
            {
                %hbytes = ();
                
                if ($metric eq kMetricAll)
                {
                    if (defined($delim))
                    {
                        $dformat = "%s${delim}%s${delim}%f${delim}%f${delim}%f${delim}%f";
                    }
                    else
                    {
                        $fformat = "%-8d%-48s%-24f%-24f%-24f%-24f";
                    }
                }
                else
                {
                    if (defined($delim))
                    {
                        $dformat = "%s${delim}%s${delim}%f";
                    }
                    else
                    {
                        $fformat = "%-8d%-48s%-24f";
                    }
                }
                
                foreach $row (@$rrows)
                {
                    if ((length($group) > 0 && $group ne $row->[0]) ||
                        (length($series) > 0 && $series ne $row->[1]))
                    {
                        PrintRow($delim, $group, $series, $metric, \%hbytes, \$dformat, \$fformat);
                        %hbytes = ();
                    }
                    
                    $group = $row->[0];
                    $series = $row->[1];
                    $hbytes{$row->[2]} = $row->[3];
                }
                
                # Must print last row since only the previous row is printed in the loop above
                PrintRow($delim, $group, $series, $metric, \%hbytes, \$dformat, \$fformat);
            }
        }
        else
        {
            print "Invalid column $order by which to order.\n";
            
        } 
    }
    elsif ($typequery eq kTypeQueryRaw)
    {
        # **** Can only work with a single container at a time when the query type is raw!! Each container uses up a 
        # huge amount of memory, so this script must be modified to handle a single one specified on the cmd-line.
        # row is series, ds_index, sudir, bytes
        my(%topheaders);
        
        # Only one of the containers should be non-empty - if that is not the case, that is an error
        if ($metric eq kMetricAll)
        {
            print "Cannot generate non-aggregate report for more than one metric.\n";
            $ok = 0;
        }
        else
        {
            %metricheaders = (kMetricDPS, "DP Now (bytes)", kMetricDPM, "DP <= 100d (bytes)", kMetricDPL, "DP > 100d (bytes)", kMetricAP, "AP (bytes)");
            %topheaders = (kMetricDPS, "*** DP Now ***", kMetricDPM, "*** DP <= 100d ***", kMetricDPL, "*** DP > 100d ***", kMetricAP, "*** AP ***");
            
            my($sunum);
            my($sudir);
            my($rowdata);
            
            if ($order eq kTypeOrderSeries)
            {
                # type - raw; order - series
                print "$topheaders{$metric}\n";
                if (defined($delim))
                {
                    $line = "series${delim}group${delim}sunum${delim}sudir${delim}$metricheaders{$metric}";
                }
                else
                {
                    $line = sprintf("%-32s%-8s%-16s%-24s%-24s", "series", "group", "sunum", "sudir", $metricheaders{$metric});
                }
                
                print "$line\n";
                
                if (CombineHashKeys(kTypeSortAlphaAsc, \@serieslist, $containers{$metric}))
                {
                    foreach $elem (@serieslist)
                    {
                        $series = $elem;
                        @grouplist = ();
                        
                        if (CombineHashKeys(kTypeSortNumrcAsc, \@grouplist, $containers{$metric}->{$series}))
                        {
                            foreach $group (@grouplist)
                            {
                                if (defined($containers{$metric}->{$series}->{$group}))
                                {
                                    # $containers{$metric}->{$series}->{$group} is reference to an array with array references
                                    # as elements. The child arrays have sunum, sudir, and bytes
                                    # as elements.
                                    foreach $rowdata (@{$containers{$metric}->{$series}->{$group}})
                                    {
                                        # $rowdata is a reference to an array containing sunum, sudir, and bytes
                                        # as elements.
                                        $sunum = $rowdata->[0];
                                        $sudir = $rowdata->[1];
                                        $metricval = $rowdata->[2];
                                        
                                        if (defined($delim))
                                        {
                                            $line = sprintf("$series${delim}$group${delim}$sunum${delim}$sudir${delim}%d", $metricval);
                                        }
                                        else
                                        {
                                            $line = sprintf("%-32s%-8d%-16s%-24s%-24d", $series, $group, $sunum, $sudir, $metricval);
                                        }
                                        
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
            elsif ($order eq kTypeOrderGroup)
            {
                my($sunum);
                my($sudir);
                my($rowdata);
                
                print "$topheaders{$metric}\n";
                if (defined($delim))
                {
                    $line =  "group${delim}series${delim}sunum${delim}sudir${delim}$metricheaders{$metric}";
                }
                else
                {
                    $line = sprintf("%-8s%-32s%-16s%-24s%-24s", "group", "series", "sunum", "sudir", $metricheaders{$metric});
                }
                
                print "$line\n";
                
                # The data in the containers are ordered by group, series. The tuples are
                # (series, ds_index, sudir, bytes). So there is very little work to do here.
                @grouplist = ();
                
                if (CombineHashKeys(kTypeSortNumrcAsc, \@grouplist, $containers{$metric}))
                {
                    foreach $group (@grouplist)
                    {
                        foreach $rowdata (@{$containers{$metric}->{$group}})
                        {
                            $series = $rowdata->[0];
                            $sunum = $rowdata->[1];
                            $sudir = $rowdata->[2];
                            $metricval = $rowdata->[3];
                            
                            if (defined($delim))
                            {
                                $line = sprintf("$group${delim}$series${delim}$sunum${delim}$sudir${delim}%d", $metricval);
                            }
                            else
                            {
                                $line = sprintf("%-8d%-32s%-16s%-24s%-24d", $group, $series, $sunum, $sudir, $metricval);
                            }
                            
                            print "$line\n";
                        }
                    }
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
    }
    
    print "__END__\n";
    
    # TODO - print out totals by group
}
