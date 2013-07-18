#!/home/jsoc/bin/linux_x86_64/perl5.12.2 

# Run like:
#  ingestIrisFDS.pl data=<data dir> host=<db host> series=<series ingesting into>
#  ingestIrisFDS.pl data=/home/arta/iris/testdata host=hmidb series=su_arta.irisfds
use strict;
use warnings;
use Data::Dumper;
use IO::Dir;
use POSIX qw(strftime);
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use drmsLocks;
use drmsArgs;
use drmsRunProg;

# Arguments
use constant kArgDataDir         => "data";   # Path to the data files to ingest.
use constant kArgSeries          => "series"; # Series into which the data files are to be ingested.
use constant kArgDBHost          => "host";   # The database host that contains the series to ingest into.

# Return codes
use constant kRetSuccess         => 0;
use constant kRetNoLock          => 1;
use constant kRetInvalidFilename => 2;

# Product names
use constant kProdOrbit          => "orbit";
use constant kProdHLZSAA         => "hlzsaa";

# Keyword names
use constant kKeyObsDate         => "obsdate";
use constant kKeyProduct         => "product";
use constant kKeyVersion         => "version";
use constant kKeyProcDate        => "date";
use constant kKeyProcDateStr     => "datestr";

# Segment names
use constant kSegProdFile        => "prodfile";

# Other constants
use constant kLockFile           => "/home/jsoc/locks/ingestIrisFds.txt";

my($lock);
my($args);
my(@dfiles);
my($numFiles);
my($prefix);
my($obsdate);
my($prod);
my($version);
my($procdate);
my($cmd);
my($rv);

$lock = new drmsNetLocks(&kLockFile);

if (defined($lock))
{
    # Read-in command-line arguments.
    $rv = &kRetSuccess;
    $args = GetArgs();
    
    if (defined($args))
    {
        
        
        @dfiles = GetFileList($args->Get(&kArgDataDir));
        
        # Parse file name to determine the product.
        # @contents contains full absolute paths
        $numFiles = 0;
        foreach my $filename (@dfiles)
        {
            if ($rv != &kRetSuccess)
            {
                last;
            }
            
            print STDOUT "Analyzing file $filename...\n";
            
            # Extract prefix (which tells us which product the file is a part of), observation date, and file version
            # from the filename.
            if ($filename =~ /([^.]+)_([^.]+)\.([^.]+)\.txt$/)
            {
                $prefix = $1;
                $obsdate = $2;
                $version = $3;
                
                if ($obsdate =~ /^([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])$/)
                {
                    # Convert the date to something DRMSy - per Rock, assume UTC.
                    $obsdate = "$1\.$2\.$3_00:00_UTC";
                    
                    if ($version =~ /^V(\d\d)/)
                    {
                        $version = sprintf("%d", $1);
                    }
                    else
                    {
                        $rv = &kRetInvalidFilename;
                        print STDERR "File name $filename does not conform to a recognized format.\n";
                    }
                }
                else
                {
                    $rv = &kRetInvalidFilename;
                    print STDERR "File name $filename does not conform to a recognized format.\n";
                }
            }
            else
            {
                $rv = &kRetInvalidFilename;
                print STDERR "File name $filename does not conform to a recognized format.\n";
            }
            
            if ($rv == &kRetSuccess)
            {
                if ($prefix =~ /^IRIS_stanford2/ || $prefix =~ /^IRIS_HLZ_SAA/)
                {
                    $prod = &kProdOrbit;
                }
                elsif ($prefix =~ /^IRIS_stanford/ || $prefix =~ /^IRIS_orbit/)
                {
                    $prod = &kProdHLZSAA;
                }
                
                if (!defined($prod) || $prefix !~ /IRIS/ || ($prefix !~ /orbit/ && $prefix !~ /HLZ_SAA/ && $prefix !~ /stanford/))
                {
                    $rv = &kRetInvalidFilename;
                    print STDERR "File name $filename does not conform to a recognized format.\n";
                }
            }
            
            if ($rv == &kRetSuccess)
            {
                $procdate = POSIX::strftime("%Y.%m.%d_%H:%M:%S_%Z", localtime());
            }
            
            if ($rv == &kRetSuccess)
            {
                # Call set_info to ingest the product file.
                print STDOUT "File name $filename has an acceptable format.\n";
                print STDOUT "Ingesting into series " . $args->Get(&kArgSeries) . " with obsdate $obsdate and product $prod.\n";
                
                if (!SeriesExists($args->Get(&kArgSeries), $args->Get(&kArgDBHost)))
                {
                    my(@jsd);
                    my($jsdstr);
                    
                    print STDOUT "Series " . $args->Get(&kArgSeries) . " does not exist on db host " . $args->Get(&kArgDBHost) . "; creating it.\n";
                    @jsd = <DATA>;
                    $jsdstr = join("", @jsd);
                    
                    if (CreateSeries($jsdstr, $args->Get(&kArgSeries)))
                    {
                        print STDERR "Unable to create series " . $args->Get(&kArgSeries) . ".\n";
                        $rv = &kRetCreateSeries;
                    }
                    else
                    {
                        if (!SeriesExists($args->Get(&kArgSeries), $args->Get(&kArgDBHost)))
                        {
                            print STDERR "Unable to create series " . $args->Get(&kArgSeries) . ".\n";
                            $rv = &kRetCreateSeries;
                        }
                    }
                    
                    if ($rv == &kRetSuccess)
                    {
                        print STDOUT "Successfully created " . $args->Get(&kArgSeries) . ".\n";
                    }
                }
            }
            
            if ($rv == &kRetSuccess)
            {
                $cmd = "set_info -c ds=" . $args->Get(&kArgSeries) . " " . &kKeyObsDate . "=$obsdate " . &kKeyProduct . "=$prod " . &kKeyVersion . "=$version " . &kKeyProcDate . "=$procdate " . &kKeyProcDateStr . "=$procdate " . &kSegProdFile . "=" . $args->Get(&kArgDataDir) . "/$filename JSOC_DBHOST=" . $args->Get(&kArgDBHost);
                
                # Duplicate STDERR
                open(STDERRDUP, ">&STDERR");
                # Redirect STDERR
                open(STDERR, ">/dev/null");
                if (drmsSysRun::RunCmd("$cmd 1 > /dev/null 2>&1") != 0)
                {
                    # Error calling set_info
                    print STDOUT "Unable to call set_info successfully. Cmd was $cmd.\n";
                    $rv = &kRetSetInfo;
                }
                close(STDERR);
                # Restore STDERR
                open(STDERR, ">&STDERRDUP");
                close(STDERRDUP);
            }
            
            $numFiles++;
        }
    }
    else
    {
        print STDERR "Invalid command-line arguments.\n";
    }
}
else
{
    print STDERR "The IRIS data-product ingester is already running; bailing out.\n";
    $rv = &kRetNoLock;
}

exit($rv);

sub GetArgs
{
    my($argsinH);
    
    # These are required arguments.
    $argsinH =
    {
        &kArgDataDir       => 's',
        &kArgSeries        => 's',
        &kArgDBHost        => 's'
    };
    
    return new drmsArgs($argsinH, 1);
}

sub GetFileList
{
    my($dfdir) = @_; # Data-file dir
    
    my(%files);      # Hash variable initialized by the IO::Dir constructor
    my(@fileList);   # List of the names of the files in the directory specified by $dfdir.
    my($file);       # The name of a single file.
    my(@rv);
    
    if (defined(tie(%files, "IO::Dir", $dfdir)))
    {
        @fileList = keys(%files);
        
        while (defined($file = shift(@fileList)))
        {
            if ($file =~ /^\.$/ || $file =~ /^\.\.$/)
            {
                # Skip the "." and ".." files
                next;
            }
            
            push(@rv, $file);
        }
        
        # Release the directory object.
        untie(%files);
    }
    
    return @rv;
}

sub SeriesExists
{
    my($series, $dbhost) = @_;
    my($cmd);
    my($rv);
    
    # A hack - just use show_info -j. If the series does not exist, then show_info returns 1. Otherwise, 
    # it returns 0.
    $rv = 0;
    
    # Duplicate STDERR
    open(STDERRDUP, ">&STDERR");
    # Redirect STDERR
    open(STDERR, ">/dev/null");
    $cmd = "show_info JSOC_DBHOST=$dbhost -j $series 1> /dev/null 2>&1";
    if (drmsSysRun::RunCmd($cmd) == 0)
    {
        $rv = 1;
    }
    else
    {
        # Error calling show_info, or the series does not exist.
    }
    close(STDERR);
    # Restore STDERR
    open(STDERR, ">&STDERRDUP");
    close(STDERRDUP);
 
    return $rv;
}

sub CreateSeries
{
    my($jsdstr, $sname) = @_;
    my($pipe);
    my($rsp);
    my($rv);
    
    $rv = 0;
        
    if (defined($jsdstr) && length($jsdstr) > 0 && defined($sname) && length($sname) > 0)
    {
        $jsdstr =~ s/\!\!SeRIesNAmE\!\!/$sname/g;
    }
    else
    {
        print STDERR "CreateSeries(): Missing series name.\n";
        $rv = 1;
    }
    
    if ($rv == 0)
    {
        if (length($jsdstr) > 0)
        {
            # Create a bidirectional-pipe so we can pass the jsd to create_series via stdin.
            $pipe = new drmsPipeRun("create_series -i");

            if (defined($pipe))
            {
                $pipe->EnableAutoflush();
                $pipe->WritePipe($jsdstr);
                $pipe->ClosePipe(1); # close write pipe
                $pipe->ReadPipe(\$rsp);

                # close both read and write pipes (but write was already closed)
                if ($pipe->ClosePipe())
                {
                    print STDERR "Failure writing to pipe.\n";
                    $rv = 1;
                }
            }
            else
            {
                print STDERR "Unable to call create_series.\n";
                $rv = 1;
            }
        }
        else
        {
            print STDERR "Empty JSD.\n";
            $rv = 1;
        }
    }

    return $rv;
}

# The data section
__DATA__
Seriesname:             !!SeRIesNAmE!!
Author:                 "Art Amezcua"
Owner:                  production
Unitsize:               2
Archive:                0
Retention:              10
Tapegroup:              1
PrimeKeys:              obsdate, product, version
Description:            "This series contains ingested, IRIS FDS product files."

#=====Keywords=====
Keyword:obsdate, time, ts_eq, record, DRMS_MISSING_VALUE, 0, UTC, "Date embedded in product file name"
Keyword:obsdate_epoch, time, constant, record, 1993.01.01_12:00:00_UTC, 0, UTC, "MDI epoch - adjusted by 12 hours to center slots on noon of each day"
Keyword:obsdate_step, time, constant, record, 1.000000, %f, day, "Slots are 1 day wide"
Keyword:product, string, variable, record, Unidentified, %s, NA, "FDS data product"
Keyword:version, int, variable, record, DRMS_MISSING_VALUE, %d, "NA", "Version of data file"
Keyword:date, time, variable, record, DRMS_MISSING_VALUE, 0, ISO, "Date of product-file ingestion; ISO 8601"
Keyword:datestr, string, variable, record, Unknown, %s, NA, "Human-readable date of product-file ingestion"

#=====Segments=====
Data: prodfile, variable, string, 0, NA, generic, "Product file"
