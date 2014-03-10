#!/home/jsoc/bin/linux_x86_64/activeperl

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
use FindBin qw($Bin);
use lib ("$Bin/../../../base/libs/perl");
use drmsArgs;

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
            $err = !(NoErr($rows, \$dbh, $stmnt));
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
                ($firstrec, $lastrec) = FindRecs($dbh, $aseries, $bcrash, $ecrash);
                
                if ($firstrec == -1 || $lastrec == -1)
                {
                    # No records in series - onto next series.
                    print STDERR "skipping series $aseries - there were no records created during the crash window.\n";
                    next;
                }
                
                $stmnt = "SELECT recnum, sunum FROM $aseries WHERE recnum >= $firstrec AND recnum <= $lastrec";
                
                $rowstab = $dbh->selectall_arrayref($stmnt, undef);
                $err = !(NoErr($rowstab, \$dbh, $stmnt));
                    
                if (!$err)
                {
                    @rowarrtab = @$rowstab;
                    
                    if ($#rowarrtab >= 0)
                    {
                        # At least one output line - print a header.
                        print "seriesname\trecnum\tsunum\n";
                        print $fh "seriesname\trecnum\tsunum\n";
                    }
                    
                    foreach my $rec (@rowarrtab)
                    {
                        print "$aseries\t$rec->[0]\t$rec->[1]\n";
                        print $fh "$aseries\t$rec->[0]\t$rec->[1]\n";
                    }
                }
                else
                {
                    $rv = kRetDbQuery;
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

sub FindRecs
{
    my($dbh,
       $series,
       $bcrash,
       $ecrash) = @_;
    my($stable);
    my($stmnt);
    my($state) = "before";
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
    my($found);     # If 1, then we found a record created during the crash window.
    my($fincrash);  # First record in the crash window.
    my($lincrash);  # Last record in the crash window.
    my($frecF);     # Upper bound when b-searching for first record in crash window.
    my($lrecF);     # Lower bound when b-searching for first record in crash window.
    my($frecL);     # Lower bound when b-searching for last record in crash window.
    my($lrecL);     # Upper bound when b-searching for last record in crash window.
    my($rows);
    my(@rowarr);
    my(@rv);
    my($err);
    
    $stable = lc($series);
    $stmnt = "SELECT recnum FROM $stable WHERE recnum = (SELECT min(recnum) FROM $stable)";
    $rows = $dbh->selectall_arrayref($stmnt, undef);
    $err = !(NoErr($rows, \$dbh, $stmnt));
    
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
        $err = !(NoErr($rows, \$dbh, $stmnt));
        
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
    }
    
    if (!$err)
    {   
        $frec = $min;
        $lrec = $max;
        
        # The first record may actually not be in the before state, but in the in state or in the after state
        ($locF, $recN) = GetLoc($dbh, $stable, $frec, 1, $min, $max, $bcrash, $ecrash);
        ($locL, $recX) = GetLoc($dbh, $stable, $lrec, 0, $min, $max, $bcrash, $ecrash);
        
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
                    $rec = $min;
                }
                else
                {
                    # frec is after the crash window - NOTHING TO RETURN (the series is completely
                    # outside the crash window).
                    $err = 1;
                }
            }
        }
        
        if (!$err)
        {
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
                    
                    if ($rec == $frec)
                    {
                        # Was not able to find a NEW record with a record number ge to the original
                        # recnum. Try finding a new record by searching in the opposite direction
                        # (toward $max).
                        $rec = $frec + ($lrec - $frec) / 2;
                        ($loc, $rec) = GetLoc($dbh, $stable, $rec, 0, $frec, $lrec, $bcrash, $ecrash);
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
                    $found = 1;
                    
                    if (!defined($fincrash))
                    {
                        # $recL is inside crash window, $frecL is before crash window (or first record
                        # in crash window)
                        $fincrash = FindFirstRec($dbh, $stable, $frec, $lrec, $bcrash, $ecrash, $frecF, $lrecF);
                    }
                    
                    if (!defined($lincrash))
                    {                        
                        # $recF is inside crash window, $lrecF is after crash window (or last record in 
                        # (crash window)
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
                    
                    if ($rec == $lrec)
                    {
                        # Was not able to find a NEW record with a record number le to the original
                        # recnum. Try finding a new record by searching in the opposite direction
                        # (toward $min).
                        $rec = $frec + ($lrec - $frec) / 2;
                        ($loc, $rec) = GetLoc($dbh, $stable, $rec, 1, $frec, $lrec, $bcrash, $ecrash);
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
    
    if ($err || $state ne "found")
    {
        return (-1, -1);
    }
    else
    {
        return ($fincrash, $lincrash);
    }
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
    my($sessionns);
    my($sessionid);
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
    @win = CheckWindow($dbh, $stable, $recno, $recno, $down);

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
            @win = CheckWindow($dbh, $stable, $wincls, $winhaf, $down);
            
            if ($#win < 0)
            {
                # Check far window
                @win = CheckWindow($dbh, $stable, $down ? $winhaf - 1 : $winhaf + 1, $winfar, $down);
                
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
                $sessionid = $win[0]->[1];
                $sessionns = $win[0]->[2];
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
        $sessionid = $win[0]->[1];
        $sessionns = $win[0]->[2];
    }
    
    if (!$err)
    {
        # The record is IN the crash window if either the endtime of the record's session is in the 
        # crash window, or the starttime of the record's session is in the crash window. So, fetch
        # the endtime and the starttime.
        $stmnt = "SELECT starttime, endtime FROM $sessionns.drms_session WHERE sessionid = $sessionid";
        
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        $err = !(NoErr($rows, \$dbh, $stmnt));
        
        if (!$err)
        {
            @rowarr = @$rows;
            
            if ($#rowarr != 0)
            {
                print STDERR "Unexpected number of rows.\n";
                $err = 1;
            }
            else
            {
                $starttime = $rowarr[0]->[0];
                $endtime = $rowarr[0]->[1];
                
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
        }
    }
    
    return @rv;
}

# Returns (recnum, sessionid, sessionns) tuples for all records specified by the range [$beg-$end]. If $down == 1, then
# the range is [$end-$beg].
sub CheckWindow
{
    my($dbh,
       $stable,
       $beg,
       $end,
       $down) = @_;
    
    my($stmnt);
    my($rows);
    my(@rowarr);
    my($err);
    my(@rv);

    if ($down)
    {
        if ($beg < $end)
        {
            @rv = [];
        }
        else
        {
            $stmnt = "SELECT recnum, sessionid, sessionns from $stable WHERE recnum <= $beg AND recnum >= $end";
        }
    }
    else
    {
        if ($beg > $end)
        {
            @rv = [];
        }
        else
        {
            $stmnt = "SELECT recnum, sessionid, sessionns from $stable WHERE recnum >= $beg AND recnum <= $end";
        }
    }

    $rows = $dbh->selectall_arrayref($stmnt, undef);
    $err = !(NoErr($rows, \$dbh, $stmnt));
    
    if ($err)
    {
        @rv = [];
    }
    else
    {   
        @rowarr = @$rows;
        
        foreach my $row (@rowarr)
        {
            push(@rv, [$row->[0], $row->[1], $row->[2]]);
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
