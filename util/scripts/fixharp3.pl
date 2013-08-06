#!/home/jsoc/bin/linux_x86_64/activeperl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use FileHandle;
use Data::Dumper;
use DBI;
use DBD::Pg;

use constant kDatafile => "/home/arta/fixharp/harps-2012.07-to-discard.recnums";
use constant kSeries   => "hmi.Mharp_720s";

my($fh);
my(@lines);
my($recnum);
my($harpno);
my($trec);
my($dsn);
my($dbh);
my($sth);
my($rv);

$rv = 1;

if (defined(open($fh, &kDatafile)))
{
    @lines = <$fh>;
    $fh->close();
    
    # Make a prepared statement that will be used repeatedly to delete records from the series table.
    # Connect to the database. Hard-code the essential information (not my best script in the world).
    $dsn = "dbi:Pg:dbname=jsoc;host=hmidb;port=5432";
    $dbh = DBI->connect($dsn, "postgres", '',  { AutoCommit => 0 }); # will need to put pass in .pg_pass

    if (defined($dbh))
    {
        $sth = $dbh->prepare("DELETE FROM " . &kSeries . " WHERE recnum = ?");
        
        if (defined($sth))
        {
            foreach my $line (@lines)
            {
                chomp($line);
                ($recnum, $harpno, $trec) = split(qr/\s+/, $line);

                if (!defined($sth->execute($recnum)))
                {
                    # error - rollback
                    $rv = 1;
                    last;
                }

                print "Executing statement " . $sth->{Statement} . ", $recnum\n";
            }
        }
        else
        {
            print STDERR "Failed to create prepared statement.\n";
            $rv = 1;
        }

        if ($rv == 0)
        {
            # commit
            $dbh->commit();
        }
        else
        {
            # rollback
            $dbh->rollback();
        }

        $dbh->disconnect();
    }
    else
    {
        print STDERR "Unable to connect to the database.\n";
        $rv = 1;
    }
}
else
{
    print STDERR "Cannot open " . &kDatafile . " for reading.\n";
    $rv = 1;
}

exit($rv);
