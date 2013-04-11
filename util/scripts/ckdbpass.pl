#!/home/jsoc/bin/linux_x86_64/perl5.12.2                                                       

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
use constant kArgPword    => "password";
use constant kOptDoit     => "doit";

# Return values
use constant kRetSuccess     => 0;
use constant kRetInvalidArgs => 1;
use constant kRetDbconnect   => 2;
use constant kRetDbQuery     => 3;

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
my($pass);      # the password to test
my($doit);      # If set, then modify database
my($stmnt);
my($rows);
my(@rowarr);
my($tstdbh);
my($user);
my($rv);

$rv = &kRetSuccess;

# Collect arguments
$argsinH =
{
    &kArgDbname  =>  's',
    &kArgDbhost  =>  's',
    &kArgDbport  =>  's',
    &kArgDbuser  =>  's',
    &kArgPword   =>  's'
};

@args = GetArgs($argsinH);

$optsinH =
{
    &kOptDoit  =>   'noval'
};

$opts = GetOpts($optsinH);

if (@args)
{
    ($dbname, $dbhost, $dbport, $dbuser, $pass) = @args;
    
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
        # Get a list of all db users
        $stmnt = "SELECT rolname FROM pg_catalog.pg_roles";
        
        $rows = $dbh->selectall_arrayref($stmnt, undef);
        if (NoErr($rows, \$dbh, $stmnt))
        {
            @rowarr = @$rows;
            foreach my $row (@rowarr)
            {
                $user = $row->[0];

                # Try to login as $user.
                $tstdbh = DBI->connect($dsn, $user, $pass, {AutoCommit => 0, PrintError => 0});
                if (defined($tstdbh))
                {
                    print "!!! User $user has bad password '$pass' !!!\n";
                    $tstdbh->rollback();
                }
            }
        }
        else
        {
            print STDERR "Unable to obtain user names from db.\n";
            $rv = &kRetDbquery;
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
    else
    {
        print STDERR "Unable to connect to db to get user names.\n";
        $rv = &kRetDbconnect;
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
                        $arg = $args->Get(&kArgPword);
                        if (defined($arg))
                        {
                            push(@rv, $arg);
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
