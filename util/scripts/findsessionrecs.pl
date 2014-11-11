#!/home/jsoc/bin/linux_x86_64/activeperl

# NOTE - THIS WAS MY FIRST ATTEMPT AT GETTING THIS TO WORK IN THE FACE OF HAVING REMOVED EARLY drms_session LOGS. THE
# QUERIES IN CacheRecs() RUN TOO SLOWLY THOUGH. THERE IS NO WAY TO JOIN THE SERIES TABLES WITH THE drms_session TABLES
# IN A WAY THAT RUNS QUICK ENOUGH. FOR ONE TABLE WITH 100M RECORDS, I LET THE SECOND QUERY RUN ABOUT 14 HOURS BEFORE
# TERMINATING IT.

# For each series in the specified namespaces, the script iterates through all records in the series and determines
# whether each one was created within the crash window. For each record, the sessionid is used to find the relevant
# record in the drms_session table. Then the starttime and endtime values are compared to the bcrash and ecrash
# values. If the record was created during the crash window, then it gets reported in a list of records that were
# created during the crash window.

use warnings;
use strict;

use POSIX qw(ceil floor);
use DBI;
use DBD::Pg;
use Data::Dumper;
use DateTime::Format::Strptime;
use FindBin qw($Bin);
use lib ("$Bin/../../../base/libs/perl");
use drmsArgs;

use constant DEBUG_ON => 0;

use constant kRetSuccess => 0;
use constant kRetInvalidArgs => 1;
use constant kRetDbQuery => 2;
use constant kRetFileIO => 3;

use constant kArgDbname => "dbname";
use constant kArgDbhost => "dbhost";
use constant kArgDbport => "dbport";
use constant kArgExclude => "exclude"; # A list of namespaces that contain series to exclude from examination.
use constant kArgInclude => "include";
use constant kArgSeries => "series"; # A list of series to examine. Overrides the exclude and include parameters.
use constant kArgBcrash => "bcrash";
use constant kArgEcrash => "ecrash";
use constant kArgFileOut => "out";

# Globals
my($gRecCache) = {};

my($argsinH);
my($optsinH);
my($args);
my($opts);
my($dbname);
my($dbhost);
my($dbport);
my($excludeList);
my(@excludedNs);
my($exclude);
my($ignoreExcInc);
my($bcrash);
my($includeList);
my(@includedNs);
my($seriesList);
my(@series);
my($ecrash);
my($fileout);
my($fh);
my($dbuser);
my($dsn);
my($dbh);
my($nspace);
my($table);
my($stmnt);
my($rows);
my(@rowarr);
my($rowstab);
my(@rowarrtab);
my($firstrec);
my($lastrec);
my($err);
my($rv);

$rv = &kRetSuccess;

$argsinH =
{
    &kArgDbname => 's',
    &kArgDbhost => 's',
    &kArgDbport => 'i',
    &kArgExclude => 's', # make this required so that people don't accidentally request processing of
                         # 43K series!
    &kArgBcrash => 's',
    &kArgEcrash => 's',
    &kArgFileOut => 's'
};

$optsinH =
{
    &kArgInclude => 's', # Inlcude applies after Exclude has filtered-out series.
    &kArgSeries => 's'   # If provided, then process these series only, and ignore the exclude/include flags.
};

$args = new drmsArgs($argsinH, 1);
$opts = new drmsArgs($optsinH, 0);

if (!defined($args) || !defined($opts))
{
    $rv = &kRetInvalidArgs;
}
else
{
    $dbname = $args->Get(&kArgDbname);
    $dbhost = $args->Get(&kArgDbhost);
    $dbport = $args->Get(&kArgDbport);
    $excludeList = $args->Get(&kArgExclude);
    $includeList = $opts->Get(&kArgInclude);
    $seriesList = $opts->Get(&kArgSeries);
    $bcrash = $args->Get(&kArgBcrash);
    $ecrash = $args->Get(&kArgEcrash);
    $fileout = $args->Get(&kArgFileOut);
    $dbuser = $ENV{USER};
    
    @excludedNs = split(/,|\s+/, $excludeList);
    if (defined($includeList))
    {
        @includedNs = split(/,|\s+/, $includeList);
    }
    
    if (defined($seriesList))
    {
        @series = split(/,|\s+/, $seriesList);
    }
    
    if (open($fh, ">$fileout"))
    {
        $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
        $dbh = DBI->connect($dsn, $dbuser, ''); # will need to put pass in .pg_pass
    }
    else
    {
        print STDERR "Unable to open output file '$fileout' for writing.\n";
        $rv = &kRetFileIO;
    }
    
    if (!$rv && defined($dbh))
    {
        # Get a list of all series we want to investigate
        $err = 0;
        $ignoreExcInc = 0;
        
        if ($#series < 0)
        {
            # No series specified on command line. Search through namespaces.
            $stmnt = "SELECT seriesname FROM drms_series()";
            $rows = $dbh->selectall_arrayref($stmnt, undef);
            $err = (NoErr($rows, \$dbh, $stmnt) == 1) ? 0 : 1;

            if (!$err)
            {
                @rowarr = @$rows;
                @series = map({$_->[0]} @rowarr);
            }
        }
        else
        {
            $ignoreExcInc = 1;
        }

        if (!$err)
        {
            my($seriesRecs);
            
            foreach my $aseries (@series)
            {
                ($nspace, $table) = ($aseries =~ /(\S+)\.(\S+)/);

                $exclude = 0;
                
                if ($ignoreExcInc)
                {
                    # Ignore @excludledNs and @includedNs. Process all series in @series.
                }
                else
                {
                    if ($#excludedNs >= 0)
                    {
                        foreach my $iexc (@excludedNs)
                        {
                            if ($nspace =~ /$iexc/)
                            {
                                $exclude = 1;
                                last;
                            }
                        }
                    }
                    
                    if (!$exclude)
                    {
                        if ($#includedNs >= 0)
                        {
                            $exclude = 1;
                            foreach my $iinc (@includedNs)
                            {
                                if ($nspace =~ /$iinc/)
                                {
                                    $exclude = 0;
                                    last;
                                }
                            }
                        }
                    }
                }

                if ($exclude)
                {
                    print STDERR "skipping series $aseries - it was excluded.\n";
                    next;
                }
                
                print STDERR "Processing series $aseries...\n";
                $seriesRecs = FindRecs($dbh, $aseries, $bcrash, $ecrash);
                
                # Clear record cache so we do not use up memory.
                PurgeRecs($aseries);
                
                if (!defined($seriesRecs))
                {
                    # No records in series - onto next series.
                    print STDERR "skipping series $aseries - there were no records created during the crash window.\n";
                    next;
                }

                if (scalar(@$seriesRecs) > 0)
                {
                    # At least one output line - print a header.
                    print "seriesname\trecnum\tsunum\n";
                    print $fh "seriesname\trecnum\tsunum\n";
                }

                foreach my $rec (@$seriesRecs)
                {
                    print "$rec->[0]\t$rec->[1]\t$rec->[2]\n";
                    print $fh "$rec->[0]\t$rec->[1]\t$rec->[2]\n";
                }
            } # series loop
        }
        else
        {
            $rv = kRetDbQuery;
        }
    }

    if (defined($fh))
    {
        $fh->close();
    }
    
    exit $rv;
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

sub GetSessionInfo
{
    my($dbh, $table, $recno) = @_;
    my($stmnt);
    my($rows);
    my(@rowarr);
    my($err);
    my(@rv);
    
    # Get ns and session id from series table.
    $stmnt = "SELECT sessionns, sessionid from $table WHERE recnum = $recno";
    $rows = $dbh->selectall_arrayref($stmnt, undef);
    $err = (NoErr($rows, \$dbh, $stmnt) == 1) ? 0 : 1;
    
    if (!$err)
    {
        @rowarr = @$rows;
        
        if ($#rowarr != 0)
        {
            $err = 1;
        }
        else
        {
            @rv = ($rowarr[0]->[0], $rowarr[0]->[1]);
        }
    }
    
    return @rv;
}

sub HasSession
{
    my($dbh, $table, $recno) = @_;

    my($stmnt);
    my($rows);
    my(@rowarr);
    my($ns);
    my($sessionid);
    my($err);
    my($rv);

    $err = 0;

    # We do NOT want to cache any of these records. We are testing DRMS records identified by b-searching so we are hopping around
    # lots of records. Anything we cached would be repreatedly flushed. Plus caching takes a long time.
    ($ns, $sessionid) = GetSessionInfo($dbh, $table, $recno);
    
    if (defined($ns) && defined($sessionid))
    {
        # See if a session exists for this series' record.
        $stmnt = "SELECT sessionid FROM $ns\.drms_session WHERE sessionid = $sessionid";
        
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        $err = (NoErr($rows, \$dbh, $stmnt) == 1) ? 0 : 1;
        
        if (!$err)
        {
            @rowarr = @$rows;
            
            if ($#rowarr == 0)
            {
                $rv = 1;
            }
            else
            {
                $rv = 0;
            }
        }
    }

    return $rv;
}

sub GetTimes
{
    my($dbh, $table, $recno, $nsIn, $sessionidIn) = @_;
    
    my($stmnt);
    my($rows);
    my(@rowarr);
    my($ns);
    my($sessionid);
    my($err);
    my(@rv);
    
    $err = 0;
    
    # We do NOT want to cache any of these records. We are testing DRMS records identified by b-searching so we are hopping around
    # lots of records. Anything we cached would be repreatedly flushed. Plus caching takes a long time.
    if (defined($nsIn) && defined($sessionidIn))
    {
        $ns = $nsIn;
        $sessionid = $sessionidIn;
    }
    else
    {
        ($ns, $sessionid) = GetSessionInfo($dbh, $table, $recno);
    }
    
    if (defined($ns) && defined($sessionid))
    {
        # See if a session exists for this series' record.
        $stmnt = "SELECT starttime, endtime FROM $ns\.drms_session WHERE sessionid = $sessionid";
        
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        $err = (NoErr($rows, \$dbh, $stmnt) == 1) ? 0 : 1;
        
        if (!$err)
        {
            @rowarr = @$rows;
            
            if ($#rowarr == 0)
            {
                @rv = ($rowarr[0]->[0], $rowarr[0]->[1]);
            }
        }
    }
    
    return @rv;
}

# Given a $recno, find the next valid (recno exists in series table) in the direction specified by $down. The search will not find
# any record below/above $bound.
sub GetValidRec
{
    my($dbh,
    $series,
    $recno,
    $down,
    $bound) = @_;
    
    my($rv);
    my($rows);
    my(@rowarr);
    my($err);
    
    $err = 0;

    if ($down)
    {
        $stmnt = "SELECT recnum FROM " . $series . " WHERE recnum >= " . $bound . " AND recnum <= " . $recno . " ORDER BY recnum DESC LIMIT 1";
    }
    else
    {
        $stmnt = "SELECT recnum FROM " . $series . " WHERE recnum >= " . $recno . " AND recnum <= " . $bound . " ORDER BY recnum ASC LIMIT 1";
    }

    $rows = $dbh->selectall_arrayref($stmnt, undef);
    $err = (NoErr($rows, \$dbh, $stmnt) == 0);

    if (!$err)
    {
        @rowarr = @$rows;

        if ($#rowarr != 0)
        {
            $err = 1;
        }
        else
        {
            $rv = $rowarr[0]->[0];
        }
    }
    
    return $rv;
}

# Returns an array of records (series, recnum, sunum). The first element in this array if the first record in the
# crash window, and the last element is the last record in this crash window. The array elements are sorted by
# recnum.

# Now that we are deleting old drms_session logs, we need to skip DRMS records that have no drms_session rows. It is not
# possible to evaluate all DRMS records, in a timely manner, to determine if they have drms_session rows or not. But
# we can b-search and find the first DRMS record that has a drms_session row. Use the first such DRMS record as the
# value for min.

sub FindRecs
{
    my($dbh,
       $series,
       $bcrash,
       $ecrash) = @_;
    my($stable);
    my($stmnt);
    my($state) = "before"; # in - either first or last series record in crash window.
                           #
    my($min);       # first record in series.
    my($max);       # last record in series.
    my($frec);      # Lower bound for b-searching.
    my($lrec);      # Upper bound for b-searching.
    my($loc);
    my($locF);      # location of record with minimum recnum.
    my($locL);      # location of record with maximum recnum.
    my($rec);       # Current rec.
    my($recN);
    my($recX);
    my($fincrash);  # First record in the crash window.
    my($lincrash);  # Last record in the crash window.
    my($frecF);     # Upper bound when b-searching for first record in crash window.
    my($lrecF);     # Lower bound when b-searching for first record in crash window.
    my($frecL);     # Lower bound when b-searching for last record in crash window.
    my($lrecL);     # Upper bound when b-searching for last record in crash window.
    my($rows);
    my(@rowarr);
    my($rv);
    my($err);
    
    $stable = lc($series);
    $stmnt = "SELECT recnum FROM $stable WHERE recnum = (SELECT min(recnum) FROM $stable)";
    $rows = $dbh->selectall_arrayref($stmnt, undef);

    # Stupid perl.
    $err = (NoErr($rows, \$dbh, $stmnt) == 1) ? 0 : 1;

    if (!$err)
    {
        @rowarr = @$rows;
        
        if ($#rowarr < 0)
        {
            # No rows in the series table
            print STDERR "This script was run on a series with no records in it.\n";
            $err = 1;
        }
        else
        {
            $min = $rowarr[0]->[0];
        }
    }
    
    if (!$err)
    {
        $stmnt = "SELECT recnum FROM $stable WHERE recnum = (SELECT max(recnum) FROM $stable)";
        
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        $err = (NoErr($rows, \$dbh, $stmnt) == 1) ? 0 : 1;

        if (!$err)
        {
            @rowarr = @$rows;
            
            if ($#rowarr < 0)
            {
                # No rows in the series table
                print STDERR "This script was run on a series with no records in it.\n";
                $err = 1;
            }
            else
            {
                $max = $rowarr[0]->[0];
            }
        }

        # We assume that the last record in the series has a session record. If this is not the case, then skip this series.
        if (!HasSession($dbh, $stable, $max))
        {
            $err = 1;
        }
    }
    
    if (!$err)
    {
        # Given $min and $max, find the first record in the series that has a drms_session record. We will consider records
        # before this to not exist. B-search for the first record with a drms_session record. Start with $min.
        my($tryRec);
        my($oldTryRec); # From the previous iteration.
        my($validRec);
        my($hasSession);
        my($lbound);
        my($ubound);
        
        $frec = undef;
        
        $tryRec = $min;
        $lbound = $min;
        $ubound = $max;
        
        while (1)
        {
            # If $tryRec doesn't change between iterations, then there are no more records to try and we
            # did not find any records that have a session.
            if (defined($oldTryRec) && $tryRec == $oldTryRec)
            {
                $err = 1;
                last;
            }
            
            $hasSession = HasSession($dbh, $stable, $tryRec);
            
            if (!defined($hasSession))
            {
                $err = 1;
                last;
            }

            if ($hasSession)
            {
                # If the record immediately below $tryRec does not have a drms_session row, then $tryRec is the first DRMS record
                # with a drms_session row.
                $validRec = GetValidRec($dbh, $stable, $tryRec - 1, 1, $lbound);

                if (defined($validRec))
                {
                    $hasSession = HasSession($dbh, $stable, $validRec);

                    if (defined($hasSession))
                    {
                        if (!$hasSession)
                        {
                            $frec = $tryRec;
                            last;
                        }
                    }
                    else
                    {
                        $err = 1;
                        last;
                    }
                }
                else
                {
                    # If there is no record immediately below $tryRec, then $tryRec is the first DRMS record
                    # with a drms_session row.
                    $frec = $tryRec;
                    last;
                }

                $ubound = $tryRec;
                # GetValidRec() will find the next EARLIER record (it may not have a drms_session log).
                $validRec = GetValidRec($dbh, $stable, $tryRec - floor(($ubound - $lbound) / 2), 1, $lbound);

                if (!defined($validRec))
                {
                    $err = 1;
                    last;
                }
                else
                {
                    $oldTryRec = $tryRec;
                    $tryRec = $validRec;
                }
            }
            else
            {
                $lbound = $tryRec;
                
                # GetValidRec() will find the next LATER record (it may not have a drms_session log).
                $validRec = GetValidRec($dbh, $stable, $tryRec + floor(($ubound - $lbound) / 2), 0, $ubound);
                
                if (!defined($validRec))
                {
                    $err = 1;
                    last;
                }
                else
                {
                    $oldTryRec = $tryRec;
                    $tryRec = $validRec;
                }
            }
        }
    }
    
    if (!$err)
    {
        $lrec = $max;
        
        # The first record may actually not be in the before state, but in the in state or in the after state.
        
        # Because we remove old session records, some records of the series may have no session record. We cannot
        # use those records to gauge whether or not a range of records is in the crash window. Treat such records
        # as invalid or missing records. In general, we remove only old session records. If the first record of a series
        # has no session record, then try a later record (move toward the LAST record in the series.).
        
        # Starting at the first record in the series ($frec), find the first valid record heading in the
        # direction of the LAST record in the series ($max). We actually already know that $frec is a valid record.
        ($locF, $recN) = GetLoc($dbh, $stable, $frec, 0, $frec, $lrec, $bcrash, $ecrash);
        
        DebugWrite("First valid rec for $stable is $recN (tried $frec). loc is $locF.");
        
        if ($locF eq "E" || $recN == -1)
        {
            # No valid records in the series. This is the same as there being no records created during crash window.
            $state = "notfound";
        }
        else
        {
            # Starting at the last record in the series ($lrec), find the first valid record heading in the
            # direction of the FIRST record in the series ($min).
            ($locL, $recX) = GetLoc($dbh, $stable, $lrec, 1, $frec, $lrec, $bcrash, $ecrash);

            DebugWrite("Last valid rec for $stable is $recX (tried $lrec). loc is $locL.\n");
            
            if ($locF eq "I")
            {
                if ($locL eq "I")
                {
                    # Both the min and max of the series are in the crash window.
                    $fincrash = $frec;
                    $lincrash = $lrec;
                    $state = "found";
                }
                else
                {
                    # The first record is in the crash window. This is the first record in the crash window.
                    $state = "in";
                    $fincrash = $recN;
                    $frecL = $recN; # current rec in crash window; lower bound for b-search.
                    $lrecL = $lrec; # upper bound for b-search.
                }
            }
            else
            {
                if ($locL eq "I")
                {
                    # The last record is in the crash window. This is the last record in the crash window.
                    $state = "in";
                    $lincrash = $recX;
                    $frecF = $recX; # current rec in crash window; upper bound for b-search.
                    $lrecF = $frec; # lower bound for b-search.
                }
                else
                {
                    # Neither min nor max is in the crash window.
                    if ($locF eq "B")
                    {
                        $state = "before";
                        $rec = $frec;
                    }
                    else
                    {
                        # frec is after the crash window - NOTHING TO RETURN (the series is completely
                        # outside the crash window - it is after the crash window).
                        $state = "notfound";
                    }
                }
            }
        }
        
        if (!$err)
        {

            DebugWrite("Starting FindRecs() loop in state $state.");
            while(1)
            {
                if ($state eq "found" || $state eq "notfound")
                {
                    last;
                }
                elsif ($state eq "before")
                {
                    $frec = $rec;
                    $rec = $frec + ($lrec - $frec) / 2;
                    
                    ($loc, $rec) = GetLoc($dbh, $stable, $rec, 1, $frec, $lrec, $bcrash, $ecrash);
                    DebugWrite("rec is $rec, loc is $loc.");
                    if ($rec == $frec)
                    {
                        # Was not able to find a NEW record with a record number ge to the original
                        # recnum. Try finding a new record by searching in the opposite direction
                        # (toward $max).
                        $rec = $frec + ($lrec - $frec) / 2;
                        ($loc, $rec) = GetLoc($dbh, $stable, $rec, 0, $frec, $lrec, $bcrash, $ecrash);
                        DebugWrite("rec is $rec, loc is $loc.");
                    }
                    
                    if ($rec == $lrec)
                    {
                        # No records between $frec and $lrec, and we know the $frec is before
                        # the crash window, and $lrec is after. So, there are no records in the crash window.
                        $state = "notfound";
                    }
                    else
                    {
                        if ($loc eq "B")
                        {
                            $state = "before";
                        }
                        elsif ($loc eq "A")
                        {
                            $state = "after";
                        }
                        else
                        {
                            # $frec is before the crash window, $lrec is after the crash window, and $rec
                            # is in the crash window.
                            $state = "in";
                            $frecF = $rec;
                            $lrecF = $frec;
                            $frecL = $rec;
                            $lrecL = $lrec;
                        }
                    }
                }
                elsif ($state eq "in")
                {
                    if (!defined($fincrash))
                    {
                        # $frecF is inside crash window, $lrecF is after crash window (or last record in
                        # crash window).

                        DebugWrite("Finding first in crash window: rec $frecF is in crash, rec $lrecF is before crash.");
                        $fincrash = FindFirstRec($dbh, $stable, $frec, $lrec, $bcrash, $ecrash, $frecF, $lrecF);
                    }
                    
                    if (!defined($lincrash))
                    {                        
                        # $recL is inside crash window, $frecL is before crash window (or first record
                        # in crash window).
                        DebugWrite("Finding last in crash window: rec $frecL is in crash, rec $lrecL is before crash.");
                        $lincrash = FindLastRec($dbh, $stable, $frec, $lrec, $bcrash, $ecrash, $frecL, $lrecL);
                    }                
                    
                    if (defined($fincrash) && defined($lincrash))
                    {
                        $state = "found";
                    }
                }
                elsif ($state eq "after")
                {
                    $lrec = $rec;
                    $rec = $frec + ($lrec - $frec) / 2;
                    
                    ($loc, $rec) = GetLoc($dbh, $stable, $rec, 0, $frec, $lrec, $bcrash, $ecrash);
                    DebugWrite("rec is $rec, loc is $loc.");
                    if ($rec == $lrec)
                    {
                        # Was not able to find a NEW record with a record number le to the original
                        # recnum. Try finding a new record by searching in the opposite direction
                        # (toward $min).
                        $rec = $frec + ($lrec - $frec) / 2;
                        ($loc, $rec) = GetLoc($dbh, $stable, $rec, 1, $frec, $lrec, $bcrash, $ecrash);
                        DebugWrite("rec is $rec, loc is $loc.");
                    }
                    
                    if ($rec == $frec)
                    {
                        # No records between $frec and $lrec, and we know the $frec is before
                        # the crash window, and $lrec is after. So, there are no records in the crash window.
                        $state = "notfound";
                    }
                    else
                    {
                        if ($loc eq "B")
                        {
                            $state = "before";
                        }
                        elsif ($loc eq "A")
                        {
                            $state = "after";
                        }
                        else
                        {
                            # $frec is before the crash window, $lrec is after the crash window, and $rec
                            # is in the crash window.
                            $state = "in";
                            $frecF = $rec;
                            $lrecF = $frec;
                            $frecL = $rec;
                            $lrecL = $lrec;
                        }
                    }
                }            
            } # while loop
        }
    }
    
    if (!$err && $state eq "found")
    {
        # We have the first record in the crash window, and the last. Now we need to get all the records between those two record,
        # including those two records, from the cache.

        DebugWrite("First rec in crash $fincrash, last rec in crash $lincrash.");
        
        # First cache the records (some may not be cached).
        CacheRecs($dbh, $series, $fincrash, $lincrash);
        
        # Then fetch them.
        $rv = GetCachedRecs($series, $fincrash, $lincrash);
    }
    
    return $rv
}

# Given a record number as input, determine whether that record could have been created
# within the crash window. The record number provided as input might not be a valid
# record number since its value may be the result of an arithmetic calcuation. This function
# will either repeatedly subtract one from, and add one to, the provided record number, until 
# either a valid record number is retrieved, or the min/max record number is reached (the
# min and max record numbers are provided as input).

# This function returns an array whose first element is the determined location code, 
# whose second element is a valid record number. The location codes include: "E" - error, 
# "B" - the record specified was created before the crash window, "I" - the record
# specified was created during the crash window, "A" - the record specified was created
# after the crash window.

# The record-finding algorithm is currently not efficient - there is a linear search for
# the next valid record number. A b-search could be performed instead. The way this would work
# is that you'd do a query to see if there are ANY records between $recno and $min
# (if $down == 1), and if so, cut this window in to two equal-size subwindows. Then check the 
# window closest to $recno. If there are records there, then cut that window in half. If there are 
# no records in the top-level half window between $recno and ($recno - $min) / 2, then check
# the half window between ($recno - $min) / 2 and $min. Keep cutting in half till you get
# a window with one record in it.
sub GetLoc
{
    my($dbh, 
       $stable, 
       $recno, 
       $down, 
       $min, # The first record in the series - can't go below this.
       $max, # The last record in the series - can't go above this.
       $bcrash,
       $ecrash) = @_;
    my($err);
    my($stmnt);
    my($starttime);
    my($endtime);
    my($rows);
    my(@rowarr);
    my($wincls);
    my($winfar);
    my($winhaf);
    my(@win);
    my(@rv);
    
    # $rec might not be a valid recnum - there was an arithmetic manipulation performed on it below.
    # So, make it a valid recnum. If $state is before, then round down, then start subtracting 1
    # till we find a valid recnum. If $state is after, then round up and start adding 1 till we
    # find a valid recnum. If $state is in, then subtract 1 from $frec till we find a valid recnum, and
    # add 1 to $lrec till we find a valid recnum.
    
    # Initialize return value.
    @rv = ("E", -1);
    
    # First, get a valid record, and fetch its sessionid.    
    if ($down)
    {
        $recno = floor($recno);
    }
    else
    {
        $recno = ceil($recno);
    }
    
    if ($recno < $min)
    {
        $recno = $min;
    }
    
    if ($recno > $max)
    {
        $recno = $max;
    }
    
    $err = 0;
    
    # First, check to see if $recno is already a valid record. If so, we can avoid b-searching through a lot of
    # records.
    @win = CheckWindow($dbh, $stable, $recno, $recno);
    
    if ($#win < 0)
    {
        # Create the initial b-search windows.
        if ($down)
        {
            $winfar = $min;
            $wincls = $recno;
            $winhaf = floor($min + ($recno - $min) / 2);
        }
        else
        {
            $wincls = $recno;
            $winfar = $max;
            $winhaf = ceil($recno + ($max - $recno) / 2);
        }
        
        while (1)
        {
            # Check close window
            if ($down)
            {
                @win = CheckWindow($dbh, $stable, $winhaf, $wincls);
            }
            else
            {
                @win = CheckWindow($dbh, $stable, $wincls, $winhaf);
            }
            
            if ($#win < 0)
            {
                # Check far window
                if ($down)
                {
                    @win = CheckWindow($dbh, $stable, $winfar, $winhaf -1);
                }
                else
                {
                    my($add) = $winhaf + 1;
                    @win = CheckWindow($dbh, $stable, $winhaf + 1, $winfar);
                }
                
                # in case we need to cut this window in half
                $wincls = $winhaf;
            }
            else
            {
                $winfar = $winhaf;
            }
            
            if ($#win < 0)
            {
                # No records outside of $recno - done
                @rv = ("E", -1);
                $err = 1;
                last;
            }
            elsif ($#win == 0)
            {
                # One record outside of $recno - done
                $recno = $win[0]->[0];
                ($starttime, $endtime) = GetTimes($dbh, $stable, $win[0]->[0], $win[0]->[1], $win[0]->[2]);
                last;
            }
            else
            {
                # More than one record in the window - cut it in half and go again.
                # But if there are two adjacent records in the window, then 
                # make $winhaf = $wincls.
                #if ($#win == 1 && abs($win[0]->[0] - $win[0]->[0]) == 1)
                if (abs($wincls - $winfar) == 1)
                {
                    $winhaf = $wincls;
                }
                else
                {
                    if ($down)
                    {
                        $winhaf = floor($winfar + ($wincls - $winfar) / 2);
                    }
                    else
                    {
                        $winhaf = ceil($wincls + ($winfar - $wincls) / 2);
                    }
                }
            }
        }
    }
    else
    {
        # Found $recno in the series. Set $sessionid and $sessionns.
        ($starttime, $endtime) = GetTimes($dbh, $stable, $win[0]->[0], $win[0]->[1], $win[0]->[2]);
    }
    
    if (!$err)
    {
        # We have a recno now.
        @rv = ("E", $recno);
    }
    
    if (!$err)
    {
        if (defined($starttime) && defined($endtime))
        {
            my($strp);

            $strp = DateTime::Format::Strptime->new(pattern => '%Y-%m-%d %T');

            $starttime = $strp->parse_datetime($starttime);
            $endtime = $strp->parse_datetime($endtime);
            $bcrash = $strp->parse_datetime($bcrash);
            $ecrash = $strp->parse_datetime($ecrash);

            # The record is IN the crash window if either the endtime of the record's session is in the
            # crash window, or the starttime of the record's session is in the crash window.
            # one row
            DebugWrite("Time comparisons: recno is $recno, start is $starttime, end is $endtime.");
            if (($endtime ge $bcrash && $endtime le $ecrash) ||
                ($starttime ge $bcrash && $starttime le $ecrash) ||
                ($starttime le $bcrash && $endtime ge $ecrash))
            {
                # Rec was created IN the crash window.
                @rv = ("I", $recno);
            }
            elsif ($endtime lt $bcrash)
            {
                # Rec was created BEFORE the crash window.
                @rv = ("B", $recno);
            }
            else
            {
                if ($starttime gt $endtime)
                {
                    # This shouldn't happen!
                    print STDERR "This cannot conceivably happen: starttime $starttime, endtime $endtime\n";
                    $err = 1;
                }
                else
                {
                    # Rec was created AFTER the crash window.
                    @rv = ("A", $recno);
                }
            }
        }
        else
        {
            # The record probably wasn't created properly. An error of "E" will be returned.
        }
    }
    
    return @rv;
}

sub CacheRecs
{
    my($dbh,
       $series,
       $firstRec,
       $lastRec) = @_;
    
    if ($firstRec <= $lastRec)
    {
        my($stmnt);
        my($rows);
        my(@rowarr);
        
        $stmnt = "SELECT recnum, sunum, sessionns, sessionid from $series WHERE recnum >= $firstRec AND recnum <= $lastRec";

        $rows = $dbh->selectall_arrayref($stmnt, undef);
        
        if (NoErr($rows, \$dbh, $stmnt))
        {
            @rowarr = @$rows;
            foreach my $row (@rowarr)
            {
                $gRecCache->{$series}->{$row->[0]} = [$row->[0], $row->[1], $row->[2], $row->[3]];
            }
        }
        
        # Create records for invalid/missing records in [$firstRec, $lastRec] too.
        for (my $iRec = $firstRec; $iRec <= $lastRec; $iRec++)
        {
            if (!exists($gRecCache->{$series}->{$iRec}))
            {
                $gRecCache->{$series}->{$iRec} = [];
            }
        }
    }
}

sub PurgeRecs
{
    my($series) = @_;
    
    if (exists($gRecCache->{$series}))
    {
        delete($gRecCache->{$series});
    }
}

# Returns a reference to an array of valid series records (series, recnum, sunum) with recnums in [$firstRec, $lastRec]. If any record in [$firstRec, $lastRec]
# is not cached, then an error is returned (an undefined array reference).
sub GetCachedRecs
{
    my($series, $firstRec, $lastRec) = @_;
    
    my($rv);
    my($rec);
    
    if (exists($gRecCache->{$series}))
    {
        for (my $iRec = $firstRec; $iRec <= $lastRec; $iRec++)
        {
            # DebugWrite("Fetching $series:$iRec from cache.");
            
            $rec = $gRecCache->{$series}->{$iRec};
            if (scalar(@$rec) > 0)
            {
                if (!defined($rv))
                {
                    $rv = [];
                }
                
                # If scalar(@$rec) == 0, then there is no valid record for that recnum.
                push(@$rv, [$series, $rec->[0], $rec->[1]])
            }
        }
    }
    
    return $rv;
}

# Returns (recnum, sessionns, sessionid) tuples for all valid records specified by the range [$beg-$end]. To make space in the db,
# old session records get manually deleted from time to time. But in this function, it is assumed that the recnums provided
# are for records that all have session records.
sub CheckWindow
{
    my($dbh,
       $stable,
       $beg,
       $end) = @_;
    
    my($firstRec);
    my($lastRec);
    my($stmnt);
    my($rows);
    my(@rowarr);
    my($err);
    my($subarr);
    my(@rv);
    
    
    # Cache (download from db) all valid records [$beg, $end] that have not been cached.
    $firstRec = undef;
    $lastRec = undef;
    
    for (my $iRec = $beg; $iRec <= $end; $iRec++)
    {
        if (exists($gRecCache->{$stable}->{$iRec}))
        {
            if (defined($firstRec))
            {
                $lastRec = $iRec - 1; # if $firstRec is defined, then $iRec >= 1
                CacheRecs($dbh, $stable, $firstRec, $lastRec);
                $firstRec = undef;
                $lastRec = undef;
            }
        }
        else
        {
            if (!defined($firstRec))
            {
                $firstRec = $iRec;
            }
        }
    }
    
    # Might be one more block of records to cache.
    if (defined($firstRec))
    {
        $lastRec = $end;
        CacheRecs($dbh, $stable, $firstRec, $lastRec);
        $firstRec = undef;
        $lastRec = undef;
    }

    # All records [$beg, $end] have been cached. For all invalid or missing records I, scalar(@{$gRecCache->{$stable}->{I}}) == 0.
    for (my $iRec = $beg; $iRec <= $end; $iRec++)
    {
        if (scalar(@{$gRecCache->{$stable}->{$iRec}}) != 0)
        {
            # ACK!!! Cannot push an array into an array! If you try, push() will append each element in the subarray into the
            # parent array. Instead, make the subarray a reference, and push the reference into the parent array.
            # We are pushing [recnum, sessionns, sessionid] into the returned array.
            $subarr = [$gRecCache->{$stable}->{$iRec}->[0], $gRecCache->{$stable}->{$iRec}->[2], $gRecCache->{$stable}->{$iRec}->[3]];
            push(@rv, $subarr);
        }
    }
    
    return @rv;
}

sub FindFirstRec
{
    my($dbh,
       $stable,
       $min,
       $max,
       $bcrash,
       $ecrash,
       $recin, # record inside crash window
       $recout # record before crash window (or first record in crash window)
      ) = @_;
    
    return BSearch($dbh, $stable, $min, $max, $bcrash, $ecrash, $recin, $recout);
}

sub FindLastRec
{
    my($dbh,
       $stable,
       $min,
       $max,
       $bcrash,
       $ecrash,
       $recin, # record inside crash window
       $recout # record after crash window (or last record in crash window)
      ) = @_;

    return BSearch($dbh, $stable, $min, $max, $bcrash, $ecrash, $recin, $recout);    
}

# returns last record in crash window if $recin < $recout, the first record in the crash
# window if $recin > $recout, and -1 if $recin == $recout.
sub BSearch
{
    my($dbh,
       $stable,
       $min,
       $max,
       $bcrash,
       $ecrash,
       $recin,
       $recout
      ) = @_;
    
    my($rv);
    my($down);
    my($loc);
    my($rec);
    
    if ($recin == $recout)
    {
        $rv = -1;
    }
    elsif ($recin < $recout)
    {
        $down = 0;
        
        ($loc, $recout) = GetLoc($dbh, $stable, $recout, 0, $min, $max, $bcrash, $ecrash);
        
        if ($loc eq "I")
        {
            $rv = -1; # error
        }
        else
        {
            ($loc, $recin) = GetLoc($dbh, $stable, $recin, 0, $min, $max, $bcrash, $ecrash);
            
            if ($loc ne "I")
            {
                $rv = -1; # error
            }
        }
    }
    else
    {
        $down = 1;
        
        ($loc, $recout) = GetLoc($dbh, $stable, $recout, 1, $min, $max, $bcrash, $ecrash);
        
        if ($loc eq "I")
        {
            $rv = -1; # error
        }
        else
        {
            ($loc, $recin) = GetLoc($dbh, $stable, $recin, 1, $min, $max, $bcrash, $ecrash);
            
            if ($loc ne "I")
            {
                $rv = -1; # error
            }
        }
    }

    while (!defined($rv))
    {
        if ($down)
        {
            $rec = $recin; # save current rec
            $recin = $recout + ($recin - $recout) / 2; # new rec
            
            ($loc, $recin) = GetLoc($dbh, $stable, $recin, 1, $min, $max, $bcrash, $ecrash);
            
            if ($recin == $recout)
            {
                # We were not able to find a NEW record between $recin and $recout. Try searching up
                # from $recin toward $rec (the original $recin).
                $recin = $recout + ($rec - $recout) / 2; # new rec
                ($loc, $recin) = GetLoc($dbh, $stable, $recin, 0, $min, $max, $bcrash, $ecrash);
                
                if ($recin == $rec)
                {
                    # No records between $recin and $rec, so $rec is the last record in the crash window.
                    $rv = $rec;
                }
            }
            
            if (!defined($rv))
            {
                # We found a record between $recin and $recout, check its location.
                if ($loc eq "B")
                {
                    $recout = $recin;
                    $recin = $rec;
                }
                elsif ($loc ne "I")
                {
                    # error
                    print STDERR "Something went wrong 1.\n";
                    $rv = -1;
                }
                
                # If $recin is inside the crash window, then we stay in this loop, using the new $recin as the 
                # upper bound for the next iteration of b-search.
            }
        }
        else
        {
            $rec = $recin; # save current rec
            $recin = $recin + ($recout - $recin) / 2; # new rec
            
            ($loc, $recin) = GetLoc($dbh, $stable, $recin, 0, $min, $max, $bcrash, $ecrash);
            
            if ($recin == $recout)
            {
                # We were not able to find a NEW record between $recin and $recout. Try searching down
                # from $recin toward $rec (the original $recin).
                $recin = $rec + ($recout - $rec) / 2; # new rec
                ($loc, $recin) = GetLoc($dbh, $stable, $recin, 1, $min, $max, $bcrash, $ecrash);
                
                if ($recin == $rec)
                {
                    # No records between $recin and $rec, so $rec is the last record in the crash window.
                    $rv = $rec;
                }
            }
            
            if (!defined($rv))
            {
                if ($loc eq "A")
                {
                    $recout = $recin;
                    $recin = $rec;
                }
                elsif ($loc ne "I")
                {
                    # error
                    print STDERR "Something went wrong 2.\n";
                    $rv = -1;
                }
                
                # If $recin is inside the crash window, then we stay in this loop, using the new $recin as the 
                # lower bound for the next iteration of b-search.
            }
        }
    }
    
    
    return $rv;
}

sub DebugWrite
{
    my($msg) = @_;
    
    if (&DEBUG_ON)
    {
        print STDERR $msg . "\n";
    }
}
