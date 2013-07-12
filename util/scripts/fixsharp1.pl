#!/home/jsoc/bin/linux_x86_64/perl5.12.2

# Run like this:
# fixsharp1.pl dbname=jsoc dbhost=hmidb dbport=5432 type1=hmi.sharp_720s:/home/mbobra/pros/swharp/cea_fix/output_ccd_definitive.txt,hmi.sharp_720s_nrt:/home/mbobra/pros/swharp/cea_fix/output_ccd_nrt.txt type2=hmi.sharp_cea_720s:/home/mbobra/pros/swharp/cea_fix/output_cea_definitive.txt,hmi.sharp_cea_720s_nrt:/home/mbobra/pros/swharp/cea_fix/output_cea_nrt.txt d=0

use warnings;
use strict;

use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use DBI;
use DBD::Pg;
use FindBin qw($Bin);
use drmsLocks;
use drmsArgs;
use drmsRunProg;

use constant kRetSuccess     => 0;
use constant kRetInvalidArgs => 1;
use constant kRetDbQuery     => 2;

# Required parameters
use constant kArgDbname      => "dbname";
use constant kArgDbhost      => "dbhost";
use constant kArgDbport      => "dbport";
# The next two arguments can optionally specify a file containing a comma-separated list of recnums. The format
# in that case would look like:
# typeX=hmi.blah:/path/to/listfile1,hmi.blahblah:/path/to/listfile2,hmi.blahblahblah:/path/to/listfile3
use constant kArgType1Series => "type1";  # comma-separted list of type-1 series (like hmi.sharp_720s)
use constant kArgType2Series => "type2";  # comma-separted list of type-2 series (like hmi.sharp_cea_720s)

# Optional parameters
use constant kOptDo          => "d";

my($argsinH);
my($args);
my($optsinH);
my($opts);
my($dbname);
my($dbhost);
my($dbport);
my($dbuser);
my($dsn);
my($dbh);
my($doit);
my(@type1);
my(@type2);
my(@type1recnums); 
my(@type2recnums);
my($listfile);
my($fh);
my($line);
my($whereclz);
my($stmnt);
my($ns);
my($serieslower);
my($rv);

$argsinH =
{
    &kArgDbname => 's',
    &kArgDbhost => 's',
    &kArgDbport => 'i',
    &kArgType1Series => 's',
    &kArgType2Series => 's'
};

$optsinH = 
{
    &kOptDo => 'noval'
};

$args = new drmsArgs($argsinH, 1);
$opts = new drmsArgs($optsinH, 0);

if (!defined($args))
{
    $rv = &kRetInvalidArgs;
}
else
{
    $dbname = $args->Get(&kArgDbname);
    $dbhost = $args->Get(&kArgDbhost);
    $dbport = $args->Get(&kArgDbport);
    $dbuser = $ENV{USER};
    
    $doit = 0;
    if (defined($opts))
    {
        $doit = $opts->Get(&kOptDo);
        if (!defined($doit))
        {
            $doit = 0;
        }
    }
        
    $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
    $dbh = DBI->connect($dsn, $dbuser, '', {AutoCommit => 0}); # will need to put pass in .pg_pass
    
    if (defined($dbh))
    {
        # Type1 series
        GetSeries($args, &kArgType1Series, \@type1, \@type1recnums);
        
        foreach my $series (@type1)
        {
            $ns = ($series =~ /^\s*(\S+)\./)[0];
            $serieslower = lc($series);
            
            $fh = undef;
            if (defined($listfile = shift(@type1recnums)))
            {
                if (!defined(open($fh, $listfile)))
                {
                    $rv = &kRetFileIO;
                    last;
                }
            }
            
            while (1)
            {
                if (defined($fh))
                {
                    $line = <$fh>;
                    if (!defined($line))
                    {
                        $fh->close();
                        last; # done with records in THIS series; onto next series.
                    }
                    
                    chomp($line);
                    $whereclz = " WHERE recnum IN ($line)";
                }
                else
                {
                    $whereclz = "";
                }
                
                # Update the vaues of variable keywords CRPIX1, CRPIX2, CRVAL1, and CRVAL2
                $stmnt = "UPDATE $series SET crpix1 = imcrpix1 - crpix1 + 1, crpix2 = imcrpix2 - crpix2 + 1, (crval1, crval2) = (0.0, 0.0)" . $whereclz;
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }
                
                # Update the descriptions of keywords CRPIX1, CRPIX2, CRVAL1, and CRVAL2
                $stmnt = "UPDATE $ns.drms_keyword SET description = 'X coordinate of disk center with respect to lower-left corner of patch (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix1'" . $whereclz;
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }
                
                $stmnt = "UPDATE $ns.drms_keyword SET description = 'Y coordinate of disk center with respect to lower-left corner of patch (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix2'" . $whereclz;
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }
                
                $stmnt = "UPDATE $ns.drms_keyword SET description = 'X origin: (0,0) at disk center' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval1'" . $whereclz;
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }
                
                $stmnt = "UPDATE $ns.drms_keyword SET description = 'Y origin: (0,0) at disk center' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval2'" . $whereclz;
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
                
                if (!defined($fh))
                {
                    last;
                }
                
                if ($rv != &kRetSuccess)
                {
                    last;
                }                
            } # while
            
            if ($rv != &kRetSuccess)
            {
                last;
            }
        } # foreach
        
        if ($rv == &kRetSuccess)
        {
            # Type 2 series
            GetSeries($args, &kArgType2Series, \@type2, \@type2recnums);
            
            foreach my $series (@type2)
            {
                $ns = ($series =~ /^\s*(\S+)\./)[0];
                $serieslower = lc($series);
                
                $fh = undef;
                if (defined($listfile = shift(@type2recnums)))
                {
                    if (!defined(open($fh, $listfile)))
                    {
                        $rv = &kRetFileIO;
                        last;
                    }
                }
                
                while (1)
                {
                    if (defined($fh))
                    {
                        $line = <$fh>;
                        if (!defined($line))
                        {
                            $fh->close();
                            last; # done with records in THIS series; onto next series.
                        }
                        
                        chomp($line);
                        $whereclz = " WHERE recnum IN ($line)";
                    }
                    else
                    {
                        $whereclz = "";
                    }
                    
                    # Update the values and descriptions of the constant keywords CTYPE1 and CTYPE2
                    $stmnt = "UPDATE $ns.drms_keyword SET defaultval = 'CRLN-CEA', description = 'CRLN-CEA' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'ctype1'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    $stmnt = "UPDATE $ns.drms_keyword SET defaultval = 'CRLT-CEA', description = 'CRLT-CEA' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'ctype2'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    # Update the values, units, and descriptions of the constant keywords CUNIT1 and CUNIT2
                    $stmnt = "UPDATE $ns.drms_keyword SET unit = 'none', defaultval = 'degree', description = 'Degree' WHERE lower(seriesname) = '$serieslower' AND (lower(keywordname) = 'cunit1' OR lower(keywordname) = 'cunit2')" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    # Update the values of variable keywords CDELT1, CDELT2, CROTA2, CRPIX1, CRPIX2, CRVAL1, CRVAL2
                    $stmnt = "UPDATE $series SET cdelt1 = 0.03, cdelt2 = 0.03, crota2 = 0.0, crpix1 = sg_000_axis000::real / 2 + 0.5, crpix2 = sg_000_axis001::real / 2 + 0.5, crval1 = crln_obs + (londtmax + londtmin) / 2, crval2 = (latdtmax + latdtmin) / 2" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    # Update the units and descriptions of the keywords CDELT1, CDELT2, CRPIX1, CRPIX2, CRVAL1, CRVAL2
                    $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Map scale in X direction' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'cdelt1'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Map scale in Y direction' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'cdelt2'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    $stmnt = "UPDATE $ns.drms_keyword SET description = 'X coordinate of patch center with respect to lower-left corner (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix1'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    $stmnt = "UPDATE $ns.drms_keyword SET description = 'Y coordinate of patch center with respect to lower-left corner (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix2'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Longitude at center of patch' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval1'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }
                    
                    $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Latitude at center of patch' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval2'" . $whereclz;
                    
                    $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                    
                    if (!defined($fh))
                    {
                        last;
                    }
                    
                    if ($rv != &kRetSuccess)
                    {
                        last;
                    }                    
                } # while
                
                if ($rv != &kRetSuccess)
                {
                    last;
                }
            } # foreach
        }
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
    
    exit($rv);
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
    my($dbh, $stmnt, $doit, $diag) = @_;
    my($rsp);
    my($res);
    my($rv);
    
    $rv = &kRetSuccess;
    
    if ($doit)
    {
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

sub GetSeries
{
    my($args, $argname, $slistR, $rlistR) = @_;
    my(@list);
    
    @list = split(qr(,), $args->Get($argname));
    foreach my $elem (@list)
    {
        if ($elem =~ /([^:]+):(.+)/)
        {
            push(@$slistR, $1);
            push(@$rlistR, $2);
        }
        else
        {
            push(@$slistR, $elem);
        }
    }
}
