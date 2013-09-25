#!/home/jsoc/bin/linux_x86_64/activeperl

# Run like:
#  ingestIrisFDS.pl data=<data dir> host=<db host> series=<series ingesting into>
#  ingestIrisFDS.pl data=/surge40/jsocprod/iris/orbit host=hmidb port=5432 user=production series=iris.fds orbseries=iris.orbit_vectors saahlzseries=iris.saa_hlz
use strict;
use warnings;
use Data::Dumper;
use IO::Dir;
use POSIX qw(strftime);
use Fcntl qw(:mode);
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use drmsLocks;
use drmsArgs;
use drmsRunProg;

# Arguments
use constant kArgDataDir          => "data";         # Path to the data files to ingest.
use constant kArgSeries           => "series";       # Series into which the data files are to be ingested.
use constant kArgDBHost           => "host";         # The database host that contains the series to ingest into.
use constant kArgDBPort           => "port";         # The port on the database host to which the db connection is made.
use constant kArgDBUser           => "user";         # The DB user whose account will be logged into.
use constant kArgOrbSeries        => "orbseries";    # The orbit series (this script will ingest orbit files into this series).
use constant kArgSAAHLZSeries     => "saahlzseries"; # The SAA_HLZ series (this script will ingest SAA_HLA files into this series).


# Return codes
use constant kRetSuccess          => 0;
use constant kRetInvalidFilename  => 1;
use constant kRetCreateSeries     => 2;
use constant kRetShowInfo         => 3;
use constant kRetSetInfo          => 4;
use constant kRetIngestIrisOrbit  => 5;
use constant kRetIngestIrisSAAHLZ => 6;
use constant kRetTimeConvert      => 7;
use constant kRetIngestPSQL       => 8;
use constant kRetNoLock           => 9;

# Product names
use constant kProdOrbit           => "orbit";
use constant kProdHLZSAA          => "hlzsaa";

# Keyword names
use constant kKeyObsDate          => "obsdate";
use constant kKeyObsDateIndex     => "obsdate_index";
use constant kKeyProduct          => "product";
use constant kKeyVersion          => "version";
use constant kKeyProcDate         => "date";
use constant kKeyProcDateStr      => "datestr";
use constant kKeyIngested         => "ingested";

# Segment names
use constant kSegProdFile         => "prodfile";

# Other constants
use constant kLockFile            => "/home/jsoc/locks/ingestIrisFds.txt";

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
my($pipe);
my($rsp);
my($rv);

$lock = new drmsNetLocks(&kLockFile);

if (defined($lock))
{
    # make the lock file world writeable (let anybody run this script)
    chmod(S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH, &kLockFile);

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
                    $prod = &kProdHLZSAA;
                }
                elsif ($prefix =~ /^IRIS_stanford/ || $prefix =~ /^IRIS_orbit/)
                {
                    $prod = &kProdOrbit;
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
                    
                    if (CreateSeries($jsdstr, $args->Get(&kArgSeries), $args->Get(&kArgDBHost), $args->Get(&kArgDBUser)))
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
                my($ingested);
                
                # Call show_info to see if this file has been ingested already.
                $cmd = "show_info -q key=" . &kKeyIngested . " JSOC_DBHOST=" . $args->Get(&kArgDBHost) . " JSOC_DBUSER=" . $args->Get(&kArgDBUser) . " ds=" . $args->Get(&kArgSeries) . "[$obsdate][$prod][$version]";
                
                # Open read pipe.
                $pipe = new drmsPipeRun($cmd, 0);
                
                if (defined($pipe))
                {
                    $pipe->ReadPipe(\$rsp);
                    
                    # close read pipe
                    if ($pipe->ClosePipe())
                    {
                        print STDERR "Failure reading from pipe.\n";
                        $rv = &kRetShowInfo;
                    }
                    else
                    {
                        chomp($rsp);
                        $ingested = ($rsp =~ /\s*y\s*/i);
                        if ($ingested)
                        {
                            print STDOUT "Already ingested, skipping.\n";
                        }
                    }
                }
                else
                {
                    print STDERR "Unable to call show_info.\n";
                    $rv = kRetShowInfo;
                }
                
                if ($rv == &kRetSuccess && !$ingested)
                {
                    # Before ingesting the file, check the FDS series to see if it has already been ingested.
                    $cmd = "set_info -c ds=" . $args->Get(&kArgSeries) . " " . &kKeyObsDate . "=$obsdate " . &kKeyProduct . "=$prod " . &kKeyVersion . "=$version " . &kKeyProcDate . "=$procdate " . &kKeyProcDateStr . "=$procdate " . &kSegProdFile . "=" . $args->Get(&kArgDataDir) . "/$filename JSOC_DBHOST=" . $args->Get(&kArgDBHost) . " " . &kKeyIngested . "=N JSOC_DBUSER=" . $args->Get(&kArgDBUser);
                    
                    # Duplicate STDERR
                    open(STDERRDUP, ">&STDERR");
                    # Redirect STDERR
                    open(STDERR, ">/dev/null");
                    
                    $rsp = drmsSysRun::RunCmd("$cmd 1 > /dev/null 2>&1");
                    if ($rsp != 0)
                    {
                        # Error calling set_info
                        print STDOUT "Unable to call set_info successfully. Cmd was $cmd.\n";
                        $rv = &kRetSetInfo;
                        
                        close(STDERR);
                        # Restore STDERR
                        open(STDERR, ">&STDERRDUP");
                        close(STDERRDUP);
                    }
                    else
                    {
                        close(STDERR);
                        # Restore STDERR
                        open(STDERR, ">&STDERRDUP");
                        close(STDERRDUP);
                        
                        # Now ingest the file contents into their product-specific series.
                        undef($cmd);
                        
                        if ($prod eq &kProdOrbit)
                        {
                            $cmd = "$Bin/ingestIrisOrbit.pl series=" . $args->Get(&kArgOrbSeries) . " source=" . $args->Get(&kArgSeries) . "[$obsdate][$prod][$version] dfile=" . $args->Get(&kArgDataDir) . "/$filename host=" . $args->Get(&kArgDBHost) . " user=" . $args->Get(&kArgDBUser);
                        }
                        elsif ($prod eq &kProdHLZSAA)
                        {
                            $cmd = "$Bin/ingestIrisSAAHLZ.pl series=" . $args->Get(&kArgSAAHLZSeries) . " source=" . $args->Get(&kArgSeries) . "[$obsdate][$prod][$version] dfile=" . $args->Get(&kArgDataDir) . "/$filename host=" . $args->Get(&kArgDBHost) . " user=" . $args->Get(&kArgDBUser);
                        }
                        
                        if (defined($cmd))
                        {
                            # Duplicate STDERR
                            open(STDERRDUP, ">&STDERR");
                            # Redirect STDERR
                            open(STDERR, ">/dev/null");
                            
                            $rsp = drmsSysRun::RunCmd("$cmd 1 > /dev/null 2>&1");
                            if ($rsp != 0)
                            {
                                # Error calling ingestion program.
                                if ($prod eq &kProdOrbit)
                                {
                                    print STDOUT "Unable to call ingestIrisOrbit.pl successfully. Cmd was $cmd.\n";
                                    $rv = &kRetIngestIrisOrbit;
                                }
                                elsif ($prod eq &kProdHLZSAA)
                                {
                                    print STDOUT "Unable to call ingestIrisSAAHLZ.pl successfully. Cmd was $cmd.\n";
                                    $rv = &kRetIngestIrisSAAHLZ;
                                }
                                
                                close(STDERR);
                                # Restore STDERR
                                open(STDERR, ">&STDERRDUP");
                                close(STDERRDUP);
                            }
                            else
                            {
                                close(STDERR);
                                # Restore STDERR
                                open(STDERR, ">&STDERRDUP");
                                close(STDERRDUP);
                                
                                # Now update the FDS series to show that the data product for this record has been successfully ingested.
                                # $obsdate is a timestring, but we need to search on the int (obsdata_index). Gotta call something
                                # to convert a time string to a index-keyword value (can't use a time value - it must be a time-index value).
                                my($slotnum);
                                
                                # Open read pipe.
                                $pipe = new drmsPipeRun("timeslot series=" . $args->Get(&kArgSeries) . " tkey=" . &kKeyObsDate . " tval=$obsdate", 0);
                                
                                if (defined($pipe))
                                {
                                    $pipe->ReadPipe(\$rsp);
                                    
                                    # close read pipe
                                    if ($pipe->ClosePipe())
                                    {
                                        print STDERR "Failure reading from pipe.\n";
                                        $rv = kRetTimeConvert;
                                    }
                                    else
                                    {
                                        chomp($rsp);
                                        $slotnum = $rsp;
                                    }
                                }
                                else
                                {
                                    print STDERR "Unable to call timeslot.\n";
                                    $rv = kRetTimeConvert;
                                }
                                
                                if ($rv == &kRetSuccess)
                                {
                                    print "Updating " . $args->Get(&kArgSeries) . "[$obsdate][$prod][$version] to reflect successful ingestion...\n";
                                    
                                    $cmd = "UPDATE " . $args->Get(&kArgSeries) . " SET " . &kKeyIngested . "='Y' WHERE " . &kKeyObsDateIndex . "=$slotnum AND " . &kKeyProduct . "='$prod' AND " . &kKeyVersion . "=$version";
                                    $cmd = "psql -h " . $args->Get(&kArgDBHost) . " -p " . $args->Get(&kArgDBPort) . " -U " . $args->Get(&kArgDBUser) . " jsoc -c \"$cmd\"";
                                    
                                    # Duplicate STDERR
                                    if (drmsSysRun::RunCmd("$cmd 1>/dev/null 2>&1") != 0)
                                    {
                                        # Error calling psql
                                        print STDOUT "Unable to call psql successfully. Cmd was $cmd.\n";
                                        $rv = &kRetIngestPSQL;
                                    }
                                }
                            }
                        }
                    }
                } # Product file not ingested.
            }
            
            $numFiles++;
        } # foreach file
    }
    else
    {
        print STDERR "Invalid command-line arguments.\n";
    }
    
    $lock->ReleaseLock();
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
        &kArgDBHost        => 's',
        &kArgDBPort        => 's',
        &kArgDBUser        => 's',
        &kArgOrbSeries     => 's',
        &kArgSAAHLZSeries  => 's'
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
    my($jsdstr, $sname, $dbhost, $dbuser) = @_;
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
            $pipe = new drmsPipeRun("create_series -i JSOC_DBHOST=$dbhost JSOC_DBUSER=$dbuser");

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
Keyword:version, int, variable, record, DRMS_MISSING_VALUE, %d, NA, "Version of data file"
# DRMS has some bug in it where you cannot use data type char to hold a character; so be wastesful and use text.
Keyword:ingested, string, variable, record, N, %c, NA, "Y/N - file contents have/have not been ingested"
Keyword:date, time, variable, record, DRMS_MISSING_VALUE, 0, ISO, "Date of product-file ingestion; ISO 8601"
Keyword:datestr, string, variable, record, Unknown, %s, NA, "Human-readable date of product-file ingestion"

#=====Segments=====
Data: prodfile, variable, string, 0, NA, generic, "Product file"
