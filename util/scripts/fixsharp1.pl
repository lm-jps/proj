#!/home/jsoc/bin/linux_x86_64/perl5.12.2

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

use constant kArgDbname      => "dbname";
use constant kArgDbhost      => "dbhost";
use constant kArgDbport      => "dbport";
use constant kArgType1Series => "type1";  # comma-separted list of type-1 series
use constant kArgType2Series => "type2";  # comma-separted list of type-2 series

use constant kOptDo          => "do";

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
    }
    
    $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
    $dbh = DBI->connect($dsn, $dbuser, ''); # will need to put pass in .pg_pass
    
    if (defined($dbh))
    {
        # Type1 series
        @type1 = split(qr(,), $args->Get(&kArgType1Series));

        foreach my $series (@type1)
        {
            $ns = ($series =~ /^\s*(\S+)\./)[0];
            $serieslower = lc($series);
            
            # Update the vaues of variable keywords CRPIX1, CRPIX2, CRVAL1, and CRVAL2
            $stmnt = "UPDATE $series SET crpix1 = imcrpix1 - crpix1 + 1, crpix2 = imcrpix2 - crpix2 + 1, (crval1, crval2) = (0.0, 0.0)";

            $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
            if ($rv != &kRetSuccess)
            {
                last;
            }

            # Update the descriptions of keywords CRPIX1, CRPIX2, CRVAL1, and CRVAL2
            $stmnt = "UPDATE $ns.drms_keyword SET description = 'X coordinate of disk center with respect to lower-left corner of patch (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix1'";

            $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
            if ($rv != &kRetSuccess)
            {
                last;
            }

            $stmnt = "UPDATE $ns.drms_keyword SET description = 'Y coordinate of disk center with respect to lower-left corner of patch (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix2'";

            $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
            if ($rv != &kRetSuccess)
            {
                last;
            }

            $stmnt = "UPDATE $ns.drms_keyword SET description = 'X origin: (0,0) at disk center' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval1'";

            $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
            if ($rv != &kRetSuccess)
            {
                last;
            }

            $stmnt = "UPDATE $ns.drms_keyword SET description = 'Y origin: (0,0) at disk center' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval2'";

            $rv = ExeStmnt($dbh, $stmnt, $doit, "type1 series update statement: $stmnt\n");
            if ($rv != &kRetSuccess)
            {
                last;
            }
        }

        if ($rv == &kRetSuccess)
        {
            # Type 2 series
            @type2 = split(qr(,), $args->Get(&kArgType2Series));

            foreach my $series (@type2)
            {
                $ns = ($series =~ /^\s*(\S+)\./)[0];
                $serieslower = lc($series);
                
                # Update the values and descriptions of the constant keywords CTYPE1 and CTYPE2
                $stmnt = "UPDATE $ns.drms_keyword SET defaultval = 'CRLN-CEA', description = 'CRLN-CEA' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'ctype1'";
                    
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                $stmnt = "UPDATE $ns.drms_keyword SET defaultval = 'CRLT-CEA', description = 'CRLT-CEA' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'ctype2'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }
                
                # Update the values, units, and descriptions of the constant keywords CUNIT1 and CUNIT2
                $stmnt = "UPDATE $ns.drms_keyword SET unit = 'none', defaultval = 'degree', description = 'Degree' WHERE lower(seriesname) = '$serieslower' AND (lower(keywordname) = 'cunit1' OR lower(keywordname) = 'cunit2')";
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                # Update the values of variable keywords CDELT1, CDELT2, CROTA2, CRPIX1, CRPIX2, CRVAL1, CRVAL2
                $stmnt = "UPDATE $series SET cdelt1 = 0.03, cdelt2 = 0.03, crota2 = 0.0, crpix1 = sg_000_axis000 / 2 + 0.5, crpix2 = sg_000_axis001 / 2 + 0.5, crval1 = crln_obs + (londtmax + londtmin) / 2, crval2 = (latdtmax + latdtmin) / 2";
                
                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                # Update the units and descriptions of the keywords CDELT1, CDELT2, CRPIX1, CRPIX2, CRVAL1, CRVAL2
                $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Map scale in X direction' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'cdelt1'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Map scale in Y direction' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'cdelt2'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                $stmnt = "UPDATE $ns.drms_keyword SET description = 'X coordinate of patch center with respect to lower-left corner (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix1'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                $stmnt = "UPDATE $ns.drms_keyword SET description = 'Y coordinate of patch center with respect to lower-left corner (in pixels)' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crpix2'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Longitude at center of patch' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval1'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }

                $stmnt = "UPDATE $ns.drms_keyword SET unit = 'degree', description = 'Latitude at center of patch' WHERE lower(seriesname) = '$serieslower' AND lower(keywordname) = 'crval2'";

                $rv = ExeStmnt($dbh, $stmnt, $doit, "type2 series update statement: $stmnt\n");
                if ($rv != &kRetSuccess)
                {
                    last;
                }
            }
        }
    }
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
