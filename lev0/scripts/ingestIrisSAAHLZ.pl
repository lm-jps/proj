#!/home/jsoc/bin/linux_x86_64/activeperl

# This script runs ingestdata on IRIS event-file data products. It reads the product files
# and formats the content so that ingestdata can use it. We will typically be ingesting
# one file each time this script is run, so that full path to that file is provided as an
# argument. A comma-separated list of paths can also be provided, and each file will
# be ingested one at a time.

# Run like:
#   ingestIrisSAAHLZ.pl series=iris.saa_hlz source='iris.fds[2013.07.26_00:00:00_UTC][hlzsaa][5]' dfile=/surge40/jsocprod/iris/orbit/IRIS_stanford2_20130808.V02.txt host=hmidb user=production

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
use constant kHdEventType      => "Event Type";
use constant kHdStartTime      => "Start Time (UTCG)";
use constant kHdStopTime       => "Stop Time (UTCG)";
use constant kHdDuration       => "Duration (sec)";

# SAAHLZ series keyword names
use constant kKwEventType      => "eventtype";
use constant kKwStartTime      => "starttime";
use constant kKwStopTime       => "stoptime";
use constant kKwDuration       => "duration";

# SAAHLA series keyword types
use constant kKwTEventType     => "string";
use constant kKwTStartTime     => "time";
use constant kKwTStopTime      => "time";
use constant kKwTDuration      => "time";

# Other constants
use constant kLockFile         => "/home/jsoc/locks/ingestIrisSAAHLZ.txt";


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
        my($kwnametypeH);
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
        InsertHeader($headersH, $headersA, &kHdEventType, &kKwEventType);
        InsertHeader($headersH, $headersA, &kHdStartTime, &kKwStartTime);
        InsertHeader($headersH, $headersA, &kHdStopTime, &kKwStopTime);
        InsertHeader($headersH, $headersA, &kHdDuration, &kKwDuration);
        
        # Map the series keywords to their data types (key = keyword name, value = keyword type)
        $kwnametypeH = {};
        $kwnametypeH->{&kKwEventType} = &kKwTEventType;
        $kwnametypeH->{&kKwStartTime} = &kKwTStartTime;
        $kwnametypeH->{&kKwStopTime} = &kKwTStopTime;
        $kwnametypeH->{&kKwDuration} = &kKwTDuration;
        
        $cmd = "ingestdata series=$series JSOC_DBHOST=$dbhost JSOC_DBUSER=$dbuser";
        # Loop through input files, calling ingestdata for each one.
        
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
                    
                    if ($iline <= 6)
                    {
                        # Loop through the lines in this file. Skip the first six lines - they offer no useful information.
                    }
                    elsif ($iline == 7)
                    {
                        # The fourth line is the keyword header. Don't assume that the headers are in any specified order.
                        
                        # Can't parse the header line since the delimiter, white space, exists in the header values. Sigh.
                        # Instead, look for expected headers, and fail if we don't find all of them.
                        chomp($line);
                        
                        $keywordsH = {};
                        foreach my $header (@$headersA)
                        {
                            $pos = 0;
                            
                            # Find the index in the line of the first occurrence of the header string.
                            $ival = index($line, $header, $pos);
                            
                            if ($ival == -1)
                            {
                                # Missing header. Bail.
                                print STDERR "Expected header, $header, is missing.\n";
                                $rv = &kRetInvalidDFile;
                                last;
                            }
                            
                            # Map from header name to keyword name.
                            $kwname = ToKey($headersH, $header);
                            
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
                        
                        # Need to send the DRMS Keyword list, and the keyword name "source" (the last keyword of the SAAHLZ
                        # series is source. Each value for this keyword contains the record-set specification that identifies 
                        # the record in the fds series that contains the orbit data used to create the orbit-series record) 
                        # to ingestdata.                        
                        $pipe->WritePipe("$sortedKWsStr source\n");
                    }
                    elsif ($iline == 8)
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
                        my($dtype);
                        my($raw);
                        my($hms);
                        my($therest);
                        my($fracsec);
                        my($strp);
                        my($dt);
                        my($timestr);
                        my($trueline);
                        
                        chomp($line);
                        
                        if ($line =~ /^\s*(\S+)\s\s+((\d\d)\s+([a-zA-Z][a-zA-Z][a-zA-Z])\s+(\d\d\d\d)\s+(\S+))\s\s+((\d\d)\s+([a-zA-Z][a-zA-Z][a-zA-Z])\s+(\d\d\d\d)\s+(\S+))\s\s+(\S+)/)
                        {
                            $trueline = $1;
                            if ($kwnametypeH->{$sortedKWs[1]} eq "time")
                            {
                                $timestr = FormatTimeString($2);
                                if (!defined($timestr))
                                {
                                    last;
                                }
                                $trueline = $trueline . " " . $timestr;
                            }
                            else
                            {
                                print STDERR "Invalid data field $2.\n";
                                $rv = &kRetInvalidDFile;
                                last;
                            }
                            
                            if ($kwnametypeH->{$sortedKWs[2]} eq "time")
                            {
                                $timestr = FormatTimeString($7);
                                if (!defined($timestr))
                                {
                                    last;
                                }
                                $trueline = $trueline . " " . $timestr;
                            }
                            else
                            {
                                print STDERR "Invalid data field $7.\n";
                                $rv = &kRetInvalidDFile;
                                last;
                            }

                            # Must also append the record-set specification of the record in the fds series that contains the orbit data 
                            # used to create the rest of the orbit-series record.
                            $trueline = $trueline . " $12 $source\n";
                        }
                        else
                        {
                            print STDERR "Invalid data line $line.\n";
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
                
                print "$rsp\n";
                
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

# Map header name to keyword name.
sub InsertHeader
{
    # $key is the header, $value is the keyword.
    my($headersH, $headersA, $key, $value) = @_;
    
    $headersH->{$key} = $value;
    push(@$headersA, $key);
}

sub ToKey
{
    my($headersH, $header) = @_;
    
    return $headersH->{$header};
}

sub FormatTimeString
{
    my($timestring) = @_;
    my($raw);
    my($hms);
    my($fracsec);
    my($strp);
    my($dt);
    my($err);
    my($rv);
    
    $err = 0;
    if ($timestring =~ /(\d\d)\s+([a-zA-Z][a-zA-Z][a-zA-Z])\s+(\d\d\d\d)\s+(\S+)/)
    {
        $raw = "$3\.$2\.$1_";
        $hms = $4;
        
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
                $err = 1;
            }
        }
        else
        {
            print STDERR "Invalid date format.\n";
            $err = 1;
        }
        
        if (!$err)
        {
            # Finally, convert the seconds-since-an-epoch into a time string usable by DRMS.
            $rv = $dt->strftime("%Y.%m.%d_%H:%M:%S");
            
            # Must append the fractional seconds because Strptime doesn't handle fractional seconds.
            
            # Then append the time zone (which I didn't do with Strptime so I could append the fractional seconds first).
            $rv = $rv . "\.${fracsec}_UTC";
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
PrimeKeys:              starttime, stoptime
Description:            "This series contains IRIS FDS event data."

#=====Keywords=====
Keyword:starttime, time, variable, record, DRMS_MISSING_VALUE, 0, UTC, "[TIME_BEG] Start time of the event"
Keyword:stoptime, time, variable, record, DRMS_MISSING_VALUE, 0, UTC, "[TIME_END] Start time of the event"
Keyword:duration, time, variable, record, -1, %f, secs, "[DURATION] Duration of the event"
Keyword:eventtype, string, variable, record, Unknown, %s, NA, "[TYPE] Type of event (SHLZ, NHLZ, or SAA)"
Keyword:date, time, variable, record, DRMS_MISSING_VALUE, 0, ISO, "[DATE] Date of product-file ingestion; ISO 8601"
Keyword:datestr, string, variable, record, Unknown, %s, NA, "Human-readable date of product-file ingestion"
Keyword:source, string, variable, record, Unknown, %s, NA, "[SOURCE] Record query to identify source file containing orbit data"
