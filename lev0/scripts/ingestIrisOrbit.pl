#!/home/jsoc/bin/linux_x86_64/activeperl

# This script runs ingestdata on IRIS oribit-file data products. It reads the product files
# and formats the content so that ingestdata can use it. We will typically be ingesting
# one file each time this script is run, so that full path to that file is provided as an
# argument. A comma-separated list of paths can also be provided, and each file will
# be ingested one at a time.

# Run like:
#   ingestIrisOrbit.pl series=su_arta.orbitvectors source=su_arta.irisfds'[2013.07.03_UTC][orbit][1]' dfile=/home/arta/Projects/Iris/ingestorbit/data/IRIS_stanford_20130703.V01.txt host=hmidb

use strict;
use warnings;
use Data::Dumper;
use IO::Dir;
use DateTime::Format::Strptime;
use Fcntl qw(:mode);
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use drmsLocks;
use drmsArgs;
use drmsRunProg;

# Return values
use constant kRetSuccess       => 0;
use constant kRetNoLock        => 1;
use constant kRetCreateSeries  => 2;
use constant kRetFileIO        => 3;
use constant kRetBustedPipe    => 4;
use constant kRetInvalidDFile  => 5;

# Parameters
use constant kArgSeries        => "series";
use constant kArgSource        => "source"; # The record-set specification of the record in the fds series containing the orbit data file whose
                                            # content will be ingested by this script.
use constant kArgDfile         => "dfile";  # A PATH to an orbit file to ingest
use constant kArgDBHost        => "host";   # DB host
use constant kArgDBUser        => "user";   # DB user

# Data file header values
use constant kHdObsDate        => "Time  (UTCG)";
use constant kHdRSunObs        => "Solar Radius (arcsec)";
use constant kHdObsVr          => "Velocity Sun wrp IRIS (km/s)";
use constant kHdHeiXObs        => "IRIS-Sun x (ECI, km)";
use constant kHdHeiYObs        => "y (km)";
use constant kHdHeiZObs        => "z (km)";
use constant kHdGeiXObs        => "IRIS x (ECI,km)";  # Do Earth X before Helio Y and Z to disambiguate the Y and Z headers
use constant kHdGeiYObs        => "y (km)"; # Sigh - duplicate header names.
use constant kHdGeiZObs        => "z(km)";  # Sigh - why is there a space in here?
use constant kHdDSunObs        => "IRIS-Sun distance (km)";

# Keys to disambiguate header values
use constant kHdIval           => "ival";
use constant kHdName           => "name";

# Orbit series keyword names
use constant kKwObsDate        => "obsdate";
use constant kKwRsunObs        => "rsunobs";
use constant kKwObsVr          => "obsvr";
use constant kKwDsunObs        => "dsunobs";
use constant kKwHeixObs        => "heixobs";
use constant kKwHeiyObs        => "heiyobs";
use constant kKwHeizObs        => "heizobs";
use constant kKwGeixObs        => "geixobs";
use constant kKwGeiyObs        => "geiyobs";
use constant kKwGeizObs        => "geizobs";

# Other constants
use constant kLockFile         => "/home/jsoc/locks/ingestIrisOrbit.txt";


my($rv);
my($lock);
my($args);
my($series);
my($source);
my($dfile);
my($dbhost);
my($dbuser);
my($cmd);
my(@filelist);
my($pipe);
my($regexp);

$lock = new drmsNetLocks(&kLockFile);

if (defined($lock))
{
    # make the lock file world writeable (let anybody run this script)
    chmod(S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH, &kLockFile);

    # Read-in command-line arguments.
    $rv = &kRetSuccess;
    $args = GetArgs();
    
    $series = $args->Get(&kArgSeries);
    $source = $args->Get(&kArgSource);
    $dfile = $args->Get(&kArgDfile);
    $dbhost = $args->Get(&kArgDBHost);
    $dbuser = $args->Get(&kArgDBUser);
    
    if (!SeriesExists($series, $dbhost))
    {
        my(@jsd);
        my($jsdstr);
        
        print STDOUT "Series $series does not exist on db host $dbhost; creating it.\n";
        @jsd = <DATA>;
        $jsdstr = join("", @jsd);

        if (CreateSeries($jsdstr, $series, $dbhost, $dbuser))
        {
            print STDERR "Unable to create series $series.\n";
            $rv = &kRetCreateSeries;
        }
        else
        {
            if (!SeriesExists($series, $dbhost))
            {
                print STDERR "Unable to create series $series.\n";
                $rv = &kRetCreateSeries;
            }
        }
    }
    
    if ($rv == &kRetSuccess)
    {
        my($headersH);
        my($headersA);
        my($fh);
        my(@dupes);
        my($pos);
        my($ival);
        my($line);
        my($iline);
        my($kwname);
        my($keywordsH);
        my(@keywords);      # Just the keys from $keywordsH
        my(@sortedKWs);     # Sorted by position of corresponding headers in the input data file.
        my($sortedKWsStr);
        
        # Map the headers that appear in the input data file to the keyword names that appear in the DRMS series.
        $headersH = {};
        $headersA = [];
        InsertHeader($headersH, $headersA, &kHdObsDate, &kKwObsDate);
        InsertHeader($headersH, $headersA, &kHdRSunObs, &kKwRsunObs);
        InsertHeader($headersH, $headersA, &kHdObsVr, &kKwObsVr);
        InsertHeader($headersH, $headersA, &kHdHeiXObs, &kKwHeixObs);
        
        # Do Earth X before Helio Y and Z to disambiguate the Y and Z headers. For all headers that are 
        # context-specific, we must first insert the disambiguating header (parent header).
        InsertHeader($headersH, $headersA, &kHdGeiXObs, &kKwGeixObs);
        InsertHeader($headersH, $headersA, &kHdHeiYObs, &kKwHeiyObs);
        InsertHeader($headersH, $headersA, &kHdHeiZObs, &kKwHeizObs);
        InsertHeader($headersH, $headersA, &kHdGeiYObs, &kKwGeiyObs);
        InsertHeader($headersH, $headersA, &kHdGeiZObs, &kKwGeizObs);
        InsertHeader($headersH, $headersA, &kHdDSunObs, &kKwDsunObs);
        
        $cmd = "ingestdata series=$series JSOC_DBHOST=$dbhost JSOC_DBUSER=$dbuser";
        # Test this puppy out.
        # $cmd = "cat";
        
        # Loop through input files, calling ingestdata for each one.
        
        # "$cmd 1 > /dev/null 2>&1"
        # Create a write-only pipe to the ingestdata process.
        $pipe = new drmsPipeRun($cmd, 1);
        
        if (defined($pipe))
        {
            $pipe->EnableAutoflush();
            
            # Send one line at a time to ingestdata.
            if (defined(open($fh, "<$dfile")))
            {
                $iline = 0;

                while (defined($line = <$fh>))
                {
                    $iline++;
                    
                    if ($iline <= 3)
                    {
                        # Loop through the lines in this file. Skip the first three lines - they offer no useful information.
                    }
                    elsif ($iline == 4)
                    {
                        # The fourth line is the keyword header. Don't assume that the headers are in any specified order.
                        
                        # Can't parse the header line since the delimiter, white space, exists in the header values. Sigh.
                        # Instead, look for expected headers, and fail if we don't find all of them.
                        chomp($line);
                        
                        $keywordsH = {};
                        foreach my $header (@$headersA)
                        {
                            @dupes = @{$headersH->{$header}};
                            $pos = 0;
                            
                            foreach my $hdinfoH (@dupes)
                            {
                                if ($hdinfoH->{&kHdIval} == -1)
                                {
                                    last;
                                }
                                else
                                {
                                    $pos = $hdinfoH->{&kHdIval} + 1;
                                }
                            }
                            
                            # Find the index in the line of the first occurrence of the header string.
                            $ival = index($line, $header, $pos);
                            
                            if ($ival == -1)
                            {
                                # Missing header. Bail.
                                print STDERR "Expected header, $header, is missing.\n";
                                $rv = &kRetInvalidDFile;
                                last;
                            }
                            
                            # Map from header name to keyword name. There might be duplicate header names, in which case
                            # $kwname will have more than one keyword name in it (separated by '|' chars).
                            $kwname = ToKey($headersH, $header, $ival);
                            
                            # Map the keyword name to its position in the string.
                            $keywordsH->{$kwname} = $ival;
                        }
                        
                        if ($rv != &kRetSuccess)
                        {
                            last;
                        }
                        
                        # Sort the keyword array, using the position the corresponding header has in the input line.
                        # $sortedKWs[i] is the DRMS keyword that corresponds to the ith header.
                        @keywords = keys(%$keywordsH);
                        @sortedKWs = sort({$keywordsH->{$a} <=> $keywordsH->{$b}} @keywords);
                        $sortedKWsStr = join(" ", @sortedKWs);
                        
                        # Need to send the DRMS Keyword list, and the keyword name "source" (the last keyword of the orbit
                        # series is source. Each value for this keyword contains the record-set specification that identifies 
                        # the record in the fds series that contains the orbit data used to create the orbit-series record) 
                        # to ingestdata.
                        $pipe->WritePipe("$sortedKWsStr source\n");
                    }
                    elsif ($iline == 5)
                    {
                        # Skip the "-------------" separator line.
                    }
                    elsif ($line =~ /^-/)
                    {
                        # For some reason they decided to insert a last line with "--" in it. Skip!
                    }
                    else
                    {
                        # A line with keyword values only. Must make a DRMS time string out of the bizarre time values provided:
                        # 03 Jul 2013 00:16:00.00 --> 2013.07.03_00:16:00.00_UTC
                        my($raw);
                        my($hms);
                        my($therest);
                        my($fracsec);
                        my($strp);
                        my($dt);
                        my($timestr);
                        my($trueline);
                        
                        chomp($line);
                        
                        if ($line =~ /^\s*(\d\d)\s+([a-zA-Z][a-zA-Z][a-zA-Z])\s+(\d\d\d\d)\s+(\S+)\s+(\S.*)/)
                        {
                            $raw = "$3\.$2\.$1_";
                            $hms = $4;
                            $therest = $5;
                            
                            if ($hms =~ /([^\.]+)\.(\d\d)/)
                            {
                                $hms = $1;
                                $fracsec = $2;
                            }
                            
                            $raw = "${raw}${hms}_UTC";
                            
                            $strp = new DateTime::Format::Strptime(pattern => '%Y.%b.%d_%H:%M:%S_%Z', locale => 'en_US', time_zone => 'UTC');
                            
                            if (defined($strp))
                            {
                                $dt = $strp->parse_datetime($raw);
                                
                                if (!defined($dt))
                                {
                                    print STDERR "Invalid date format.\n";
                                    $rv = &kRetInvalidDFile;
                                    last;
                                }
                            }
                            else
                            {
                                print STDERR "Invalid date format.\n";
                                $rv = &kRetInvalidDFile;
                                last;
                            }
                            
                            # Finally, convert the seconds-since-an-epoch into a time string usable by DRMS.
                            $timestr = $dt->strftime("%Y.%m.%d_%H:%M:%S");
                            
                            # Must append the fractional seconds because Strptime doesn't handle fractional seconds.
                            $timestr = $timestr . "\.${fracsec}_UTC";
                            
                            # Then append the time zone (which I didn't do with Strptime so I could append the fractional seconds first).
                            # Must also append the record-set specification of the record in the fds series that contains the orbit data 
                            # used to create the rest of the orbit-series record.
                            $trueline = "$timestr $therest $source\n";
                        }
                        else
                        {
                            print STDERR "Invalid date format.\n";
                            $rv = &kRetInvalidDFile;
                            last;
                        }
                        
                        $pipe->WritePipe($trueline);
                    }
                } # processing an input line
                
                $fh->close();
            }
            else
            {
                print STDERR "Unable to open file $dfile for reading. Skipping this file.\n";
                if ($pipe->ClosePipe())
                {
                    print STDERR "Failure communicating with ingestdata.\n";
                    $rv = &kRetBustedPipe;
                    last;
                }
                next;
            }
            
            # close write pipe
            if ($pipe->ClosePipe())
            {
                print STDERR "Failure communicating with ingestdata.\n";
                $rv = &kRetBustedPipe;
                last;
            }
        }
        else
        {
            print STDERR "Unable to call ingestdata.\n";
            $rv = &kRetBustedPipe;
            last;
        }
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
        &kArgSeries        => 's',
        &kArgSource        => 's',
        &kArgDfile         => 's',
        &kArgDBHost        => 's',
        &kArgDBUser        => 's'
    };
    
    return new drmsArgs($argsinH, 1);
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

sub InsertHeader
{
    # $key is the header, $value is the keyword.
    my($headersH, $headersA, $key, $value) = @_;
    
    # There might be duplicate header names - sigh. Because of this BS, maintain an ORDERED
    # list of headers. We need the context of the head to know which header it is.
    push(@{$headersH->{$key}}, {&kHdIval => -1, &kHdName => $value});
    push(@$headersA, $key);
}

sub ToKey
{
    my($headersH, $header, $ival) = @_;
    my(@dupes);
    my($hival);
    my($gival);
    
    # First, save $ival for this header.
    @dupes = @{$headersH->{$header}};
    
    # Disambiguate - must get ival for both kHdHeiXObs and kHdGeiXObs headers. Assume that for these two headers, there
    # are no duplicates.
    if ($#dupes > 0)
    {        
        $hival = $headersH->{&kHdHeiXObs}->[0]->{&kHdIval};
        $gival = $headersH->{&kHdGeiXObs}->[0]->{&kHdIval};
        
        if ($ival > $hival && $ival < $gival)
        {
            $headersH->{$header}->[0]->{&kHdIval} = $ival;
            return $headersH->{$header}->[0]->{&kHdName};
        }
        else
        {
            $headersH->{$header}->[1]->{&kHdIval} = $ival;
            return $headersH->{$header}->[1]->{&kHdName};
        }
    }
    else
    {
        $headersH->{$header}->[0]->{&kHdIval} = $ival;
        return $headersH->{$header}->[0]->{&kHdName};
    }
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
PrimeKeys:              obsdate
Description:            "This series contains IRIS FDS orbit data."

#=====Keywords=====
Keyword:obsdate, time, ts_eq, record, DRMS_MISSING_VALUE, 0, UTC, "[OBS_DATE] Date of prediction"
Keyword:obsdate_epoch, time, constant, record, 1993.01.01_00:00:30_UTC, 0, UTC, "MDI epoch - adjusted by 30 seconds to center slots on minutes"
Keyword:obsdate_step, time, constant, record, 60.000000, %f, secs, "Slots are 60 seconds wide"
Keyword:date, time, variable, record, DRMS_MISSING_VALUE, 0, ISO, "[DATE] Date of product-file ingestion; ISO 8601"
Keyword:datestr, string, variable, record, Unknown, %s, NA, "Human-readable date of product-file ingestion"
Keyword:rsunobs, float, variable, record, DRMS_MISSING_VALUE, %f, arcsec, "[RSUN_OBS] The solar radius"
Keyword:obsvr, float, variable, record, DRMS_MISSING_VALUE, %f, km/s, "[OBS_VR] The velocity of the Sun with respect to IRIS"
Keyword:dsunobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[DSUN_OBS] IRIS-Sun distance"
Keyword:heixobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[HEIX_OBS] IRIS-Sun x position (ECI)"
Keyword:heiyobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[HEIY_OBS] IRIS-Sun y position (ECI)"
Keyword:heizobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[HEIZ_OBS] IRIS-Sun z position (ECI)"
Keyword:geixobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[GEIX_OBS] IRIS x position (ECI)"
Keyword:geiyobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[GEIY_OBS] IRIS y position (ECI)"
Keyword:geizobs, double, variable, record, DRMS_MISSING_VALUE, %f, km, "[GEIZ_OBS] IRIS z position (ECI)"
Keyword:source, string, variable, record, Unknown, %s, NA, "Record query to identify source file containing orbit data"
