#!/home/jsoc/bin/linux_x86_64/activeperl                                                       

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);
use DBI;
use DBD::Pg;
use lib "$Bin/../../../base/libs/perl";
use Data::Dumper;
use drmsArgs;
use drmsRunProg;

# Arguments
use constant kArgDbname   => "dbname";
use constant kArgDbhost   => "dbhost";
use constant kArgDbport   => "dbport";
use constant kArgDbuser   => "dbuser";
use constant kArgKeyVals  => "keyvals";
use constant kArgSeries   => "series";
use constant kOptDoit     => "doit";

# Return values
use constant kRetSuccess     => 0;
use constant kRetInvalidArgs => 1;
use constant kRetDbQuery     => 2;

my($argsinH);
my($optsinH);
my(@args);
my($opts);
my($dsn);
my($dbh);
my($dbname);    # name of the db instance to connect to
my($dbhost);    # name of the db host on which the db instance resides
my($dbport);    # port on $dbhost through which connections are made
my($dbuser);    # database user name (to log-in as)
my($keyvals);   # Full path to file containing key values
my($series);    # The series to modify
my($doit);      # If set, then modify database
my($lcseries);
my($fh);
my($stmnt);
my($line);
my($harpnum);
my($noaaar);
my($noaanum);
my($noaaars);
my($rv);

$rv = &kRetSuccess;

# Collect arguments
$argsinH =
{
    &kArgDbname  =>  's',
    &kArgDbhost  =>  's',
    &kArgDbport  =>  's',
    &kArgDbuser  =>  's',
    &kArgKeyVals =>  's',
    &kArgSeries  =>  's'
};

@args = GetArgs($argsinH);

$optsinH =
{
    &kOptDoit  =>   'noval'
};

$opts = GetOpts($optsinH);

if (@args)
{
    ($dbname, $dbhost, $dbport, $dbuser, $keyvals, $series) = @args;
    
    $doit = 0;
    if (defined($opts))
    {
        $doit = $opts->Get(&kOptDoit);
        if (!defined($doit))
        {
            $doit = 0;
        }
    }

    $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
    $dbh = DBI->connect($dsn, $dbuser, '', {AutoCommit => 0}); # will need to put pass in .pg_pass  
    
    if (defined($dbh))
    {
        # Read the keyword-values file.
        if (open($fh, "<$keyvals"))
        {
            $lcseries = lc($series);

            while (defined($line = <$fh>))
            {
                chomp($line);

                # skip blank/empty/ws lines
                if ($line =~ /^\s*$/)
                {
                    next;
                }
                
                if ($line =~ /^\#/)
                {
                    next;
                }
                
                if ($line =~ /^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\,0-9]+)\s*$/)
                {
                    $harpnum = $1;
                    $noaaar = $2;
                    $noaanum = $3;
                    $noaaars = $4;

                    if ($noaanum == 0)
                    {
                        $noaaars = "";
                    }
                    
                    $stmnt = "UPDATE $lcseries SET noaa_ar = $noaaar, noaa_num = $noaanum, noaa_ars = '$noaaars' WHERE harpnum = $harpnum";
                    ExeStmnt($dbh, $stmnt, $doit, 1, "Executing DB statement: $stmnt.\n");
                }
                else
                {
                    print STDERR "skipping bad line $line.\n";
                }
            }
            
            $fh->close();
        }
        
        if ($rv == &kRetSuccess && $doit)
        {
            # Commit the db changes.
            print "Committing changes.\n";
            $dbh->commit();
        }
        else
        {
            print "Rolling back changes, error code $rv.\n";
            $dbh->rollback();
        }
    }

    exit($rv);
}

sub GetArgs
{
    my($argsinH) = @_;
    my($args);
    my($arg);
    my(@rv);

    $args = new drmsArgs($argsinH, 1);

    if (defined($args))
    {
        $arg = $args->Get(&kArgDbname);
        if (defined($arg))
        {
            push(@rv, $arg);
            $arg = $args->Get(&kArgDbhost);
            if (defined($arg))
            {
                push(@rv, $arg);
                $arg = $args->Get(&kArgDbport);
                if (defined($arg))
                {
                    push(@rv, $arg);
                    $arg = $args->Get(&kArgDbuser);
                    if (defined($arg))
                    {
                        push(@rv, $arg);
                        $arg = $args->Get(&kArgKeyVals);
                        if (defined($arg))
                        {
                            push(@rv, $arg);
                            $arg = $args->Get(&kArgSeries);
                            if (defined($arg))
                            {
                                push(@rv, $arg);
                            }
                        }
                    }
                }
            }
        }
    }

    return @rv;
}

sub GetOpts
{
    my($argsinH) = @_;
    my($rv);

    $rv = new drmsArgs($optsinH, 0);
}

sub NoErr
{
    my($rv) = $_[0];
    my($dbh) = $_[2];
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

sub ExeStmnt
{
    my($dbh, $stmnt, $doit, $printit, $diag) = @_;
    my($rsp);
    my($res);
    my($rv);

    $rv = &kRetSuccess;

    if ($doit)
    {
        if ($printit)
        {
            print $diag;
        }
        
        $res = $dbh->do($stmnt);
        if (!NoErr($res, $dbh, $stmnt))
        {
            $rv = &kRetDbQuery;
        }
    }
    else
    {
        print $diag;
    }

    return $rv;
}
