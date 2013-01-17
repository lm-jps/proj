#!/home/jsoc/bin/linux_x86_64/perl5.12.2 

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

use constant kArgDbname => "dbname";
use constant kArgDbhost => "dbhost";
use constant kArgDbport => "dbport";
use constant kArgExclude => "exclude"; # A list of namespaces that contain series to exclude from examination.
use constant kArgBcrash => "bcrash";
use constant kArgEcrash => "ecrash";

my($argsinH);
my($args);
my($dbname);
my($dbhost);
my($dbport);
my($exclude);
my($bcrash);
my($ecrash);
my($dbuser);
my($dsn);
my($dbh);
my($seriesR);
my($series);
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
    &kArgEcrash => 's'
};

$args = new drmsArgs($argsinH, 1);

if (!defined($args))
{
    $rv = &kRetInvalidArgs;
}
else
{
    $dbname = $args->Get(&kArgDbname);
    $dbhost = $args->Get(&kArgDbhost);
    $dbport = $args->Get(&kArgDbport);
    $exclude = $args->Get(&kArgExclude);
    $bcrash = $args->Get(&kArgBcrash);
    $ecrash = $args->Get(&kArgEcrash);
    $dbuser = $ENV{USER};
    
    $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
    $dbh = DBI->connect($dsn, $dbuser, ''); # will need to put pass in .pg_pass
    
    if (defined($dbh))
    {
        # Get a list of all series we want to investigate
        $stmnt = "SELECT seriesname FROM drms_series()";
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        $err = !(NoErr($rows, \$dbh, $stmnt));

        if (!$err)
        {
            @rowarr = @$rows;
            
            foreach my $seriesR (@rowarr)
            {
                $series = $seriesR->[0];
                if (length($exclude) > 0)
                {
                    if ($series =~ /$exclude/)
                    {
                        print STDERR "skipping series $series - it was excluded.\n";
                        next;
                    }
                }
                
                print STDERR "Processing series $series...\n";
                ($firstrec, $lastrec) = FindRecs($dbh, $series);
                
                if ($firstrec == -1 || $lastrec == -1)
                {
                    # No records in series - onto next series.
                    print STDERR "skipping series $series - there were no records created during the crash window.\n";
                    next;
                }
                
                $stmnt = "SELECT recnum, sunum FROM $series WHERE recnum >= $firstrec AND recnum <= $lastrec";
                
                $rowstab = $dbh->selectall_arrayref($stmnt, undef);
                $err = !(NoErr($rowstab, \$dbh, $stmnt));
                    
                if (!$err)
                {
                    @rowarrtab = @$rowstab;
                    
                    if ($#rowarrtab >= 0)
                    {
                        # At least one output line - print a header.
                        print "seriesname\trecnum\tsunum\n";
                    }
                    
                    foreach my $rec (@rowarrtab)
                    {
                        print "$series\t$rec->[0]\t$rec->[1]\n";
                    }
                }
                else
                {
                    $rv = kRetDbQuery;
                }
            }
        }
        else
        {
            $rv = kRetDbQuery;
        }
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
    my($dbh) = $_[0];
    my($series) = $_[1];
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
        ($locF, $rec) = GetLoc($dbh, $stable, $frec, 1, $min, $max, $bcrash, $ecrash);
        if ($locF eq "I")
        {
            # The first record is in the crash window. This is the first record in the crash window.
            $state = "in";
            $fincrash = $rec;
            $frecL = $rec; # current rec in crash window; lower bound for b-search.
            $lrecL = $lrec; # upper bound for b-search.
        }
        else
        {
            ($locL, $rec) = GetLoc($dbh, $stable, $lrec, 0, $min, $max, $bcrash, $ecrash);
            if ($locL eq "I")
            {
                # The last record is in the crash window. This is the last record in the crash window.
                $state = "in";
                $lincrash = $rec;
                $frecF = $rec; # current rec in crash window; upper bound for b-search.
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
                    
                    ($loc, $rec) = GetLoc($dbh, $stable, $rec, 1, $min, $max, $bcrash, $ecrash);
                    
                    if ($rec == $frec)
                    {
                        # Was not able to find a NEW record with a record number ge to the original
                        # recnum. Try finding a new record by searching in the opposite direction
                        # (toward $max).
                        $rec = $frec + ($lrec - $frec) / 2;
                        ($loc, $rec) = GetLoc($dbh, $stable, $rec, 0, $min, $max, $bcrash, $ecrash);
                    }
                    
                    if ($rec == $max && $loc eq "B")
                    {
                        # There is no record in this series that was created during the crash window.
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
                        $fincrash = FindFirstRec($dbh, $stable, $min, $max, $bcrash, $ecrash, $frecF, $lrecF);
                    }
                    
                    if (!defined($lincrash))
                    {                        
                        # $recF is inside crash window, $lrecF is after crash window (or last record in 
                        # (crash window)
                        $lincrash = FindLastRec($dbh, $stable, $min, $max, $bcrash, $ecrash, $frecL, $lrecL);
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
                    
                    ($loc, $rec) = GetLoc($dbh, $stable, $rec, 0, $min, $max, $bcrash, $ecrash);
                    
                    if ($rec == $lrec)
                    {
                        # Was not able to find a NEW record with a record number le to the original
                        # recnum. Try finding a new record by searching in the opposite direction
                        # (toward $min).
                        $rec = $frec + ($lrec - $frec) / 2;
                        ($loc, $rec) = GetLoc($dbh, $stable, $rec, 1, $min, $max, $bcrash, $ecrash);
                    }
                    
                    if ($rec == $min && $loc == "A")
                    {
                        # There is no record in this series that was created during the crash window.
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
# the next valid record number. A b-search could be performed instead.

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
    
    while (!$err)
    {
        $stmnt = "SELECT sessionid, sessionns FROM $stable WHERE recnum = $recno";
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        $err = !(NoErr($rows, \$dbh, $stmnt));
        
        if ($err)
        {
            last;
        }
        
        @rowarr = @$rows;
        
        if ($#rowarr < 0)
        {
            # no such record, try again.
            if ($down)
            {
                if ($recno > $min)
                {
                    $recno--;
                }
                else
                {
                    # There was no record in the series with a recnum smaller than the invalid $recno provided.
                    @rv = ("E", -1);
                    $err = 1;
                }
            }
            else
            {
                if ($recno < $max)
                {
                    $recno++;
                }
                else
                {
                    # There was no record in the series greater than the invalid $recno provided.
                    @rv = ("E", -1);
                    $err = 1;
                }
            }
        }
        elsif ($#rowarr == 0)
        {
            $sessionid = $rowarr[0]->[0];
            $sessionns = $rowarr[0]->[1];
            last;
        }
        else
        {
            print STDERR "This cannot happen!! There is more than more record with the same recnum.\n";
            $err = 1;
            last;
        }
    } # while loop
    
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
