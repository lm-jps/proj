#!/home/jsoc/bin/linux_x86_64/activeperl

use strict;
use warnings;
use JSON;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use Data::Dumper;
use drmsArgs;
use drmsRunProg;

# Required cmd-line arguments
use constant kArgKeyNames    => "keynames";  # Full path to file containing key names
use constant kArgKeyVals     => "keyvals";   # Full path to file containing key values
use constant kArgBitMaps     => "bmaps";     # Path to bit maps
use constant kArgBitMapFname => "fname";     # Base file name of bit-map files
use constant kArgSeries      => "series";    # The series into which we are adding records
use constant kArgSegment     => "segment";   # The segment in series 'series' that will contain data files

# Optional cmd-line arguments
use constant kOptSetDate     => "d";         # Pass to rawingest

use constant kBatchSize   => 128;

my($argsinH);
my(@args);
my($optsinH);
my($opts);
my($kevals);
my($keynames);
my($keyvals);
my($bmaps);
my($fname);
my($series);
my($segment);
my($setdate);
my($dataH);
my($cpkeyH);
my($cpkeyParentH);
my($fh);
my($line);
my($iline);
my(@pkeys); # Prime-key keyword names.
my($countpkeys);
my($npkeys);
my(@okeys); # Non-prime-key keyword names.
my(@vals);  # All vals
my(@pkeyvals); # Prime-key keyword values.
my($ival);
my($hack);
my($fpath);
my($rv);

$rv = 0; # success

# Collect arguments
$argsinH =
{
    &kArgKeyNames  =>   's',
    &kArgKeyVals  =>    's',
    &kArgBitMaps  =>    's',
    &kArgBitMapFname => 's',
    &kArgSeries =>      's',
    &kArgSegment =>     's'
};

@args = GetArgs($argsinH);

$optsinH =
{
    &kOptSetDate  =>   'noval'
};

$opts = GetOpts($optsinH);

if (@args)
{
    ($keynames, $keyvals, $bmaps, $fname, $series, $segment) = @args;
    if (defined($opts))
    {
        $setdate = $opts->Get(&kOptSetDate);
        # $setdate will be undefined if there was no flag supplied on the cmd-line
    }
    
    # First put the keynames into an array.
    if (open($fh, "<$keynames"))
    {
        $npkeys = 0;
        while (defined($line = <$fh>))
        {
            chomp($line);
            
            # skip blank/empty/ws lines
            if ($line =~ /^\s*$/)
            {
                next;
            }
            
            # skip lines that have text enclosed by []
            if ($line =~ /^\s*\[(.+)\]/)
            {
                if ($1 =~ /prime keys/i)
                {
                    # Start $npkeys counter.
                    $countpkeys = 1;
                }
                elsif ($1 =~ /keys to change/i)
                {
                    # Stop $npkeys counter.
                    $countpkeys = 0;
                }
                next;
            }
            
            if ($line =~ /^\s*(\S+)\s*$/)
            {
                if ($countpkeys)
                {
                    push(@pkeys, $1);
                    $npkeys++;
                }
                else
                {
                    push(@okeys, $1);
                }
            }
            else
            {
                print STDERR "Unexpected format of line in $keynames: $line.\n";
                $rv = 1;
                last;
            }
        }
        
        $fh->close();
    }
    else
    {
        print STDERR "Unable to open $keynames for reading.\n";
        $rv = 1;
    }
        
    if ($rv == 0 && $npkeys > 0)
    {   
        # loop through the data file identified by $keyvals
        if (open($fh, "<$keyvals"))
        {
            $iline = 0;
            $dataH = {};
            
            while (defined($line = <$fh>))
            {
                chomp($line);
                
                # Skip comment lines.
                if ($line =~ /^\#/)
                {
                    next;
                }
                
                # Skip blank/empty/ws lines
                if ($line =~ /^\s*$/)
                {
                    next;
                }

                @vals = split(/,\s*/, $line);
                
                # The first two columns in this file contain the values for the 
                # prime-key keywords: HARPNUM and T_REC. Make a hash-array layer
                # for each of these keywords.
                $fpath = $bmaps;
                $ival = 0;
                foreach my $aval (@vals)
                {
                    if ($ival < $npkeys)
                    {
                        # These are prime-key keyword values.
                        
                        # Must hack the HARPNUM keyword - the keys file does not contain the
                        # exact HARPNUM directory name. For the HARPNUM, it contains an integer,
                        # but the directory name is a 0-padded string representation of this integer 
                        # with a total of 5 digits. %05d
                        if ($ival == 0)
                        {
                            $hack = sprintf("%05d", $aval);
                            push(@pkeyvals, $hack);
                            $fpath = "$fpath/$hack";
                        }
                        else
                        {
                            push(@pkeyvals, $aval);
                            $fpath = "$fpath/$aval";
                        }
                    }
                    else 
                    {
                        $cpkeyH = $dataH;
                        
                        foreach my $pkeyval (@pkeyvals)
                        {
                            if (!exists($cpkeyH->{$pkeyval}))
                            {
                                $cpkeyH->{$pkeyval} = {};                                
                            }

                            $cpkeyH = $cpkeyH->{$pkeyval};   
                        }

                        if (!exists($cpkeyH->{"keys"}))
                        {
                            $cpkeyH->{"keys"} = {};
                            $cpkeyParentH = $cpkeyH;
                        }
                        
                        $cpkeyH = $cpkeyH->{"keys"};
                        
                        # These are non-prime-key keyword values.
                        if (!exists($cpkeyH->{$okeys[$ival - $npkeys]}))
                        {
                            $cpkeyH->{$okeys[$ival - $npkeys]} = $aval;
                        }
                        else
                        {
                            # Error - there should not be an element for this prime-key value already.
                            $rv = 1;
                            last;
                        }
                    }
                    
                    $ival++;
                }
                
                # Set the filename of the file to ingest.
                $cpkeyParentH->{"file"} = "$fpath/$fname";

                @pkeyvals = ();
                $iline++;
                                
                # Call ingestion module with batches of records.
                if ($iline % &kBatchSize == 0)
                {
                    $rv = Ingest($series, $segment, $setdate, $dataH);
                    if ($rv != 0)
                    {
                        last;
                    }
                    
                    $dataH = {};
                    $cpkeyH = undef;
                    $cpkeyParentH = undef;
                }
            } # End while
            
            if ($rv == 0)
            {
                if (keys(%$dataH) > 0 && defined($cpkeyH) && keys(%$cpkeyH) > 0)
                {
                    $rv = Ingest($series, $segment, $setdate, $dataH);
                }
            }
            
            $fh->close();
        }
        else
        {
            print STDERR "Unable to open key-names file $keynames.\n";
            $rv = 1;
        }
    }
}
else
{
    print STDERR "Invalid arguments.\n";
}

exit($rv);

sub GetArgs
{
    my($argsinH) = @_;
    my($args);
    my($arg);
    my(@rv);
    
    $args = new drmsArgs($argsinH, 1);
    
    if (defined($args))
    {
        $arg = $args->Get(&kArgKeyNames);
        if (defined($arg))
        {
            push(@rv, $arg);
            $arg = $args->Get(&kArgKeyVals);
            if (defined($arg))
            {
                push(@rv, $arg);
                $arg = $args->Get(&kArgBitMaps);
                if (defined($arg))
                {
                    push(@rv, $arg);
                    $arg = $args->Get(&kArgBitMapFname);
                    if (defined($arg))
                    {
                        push(@rv, $arg);
                        $arg = $args->Get(&kArgSeries);
                        if (defined($arg))
                        {
                            push(@rv, $arg);
                            $arg = $args->Get(&kArgSegment);
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

sub Ingest
{
    my($series, $segment, $setdate, $dataH) = @_;
    my($cmd);
    my($json);
    my($pipe);
    my($rsp);
    my($childret);
    my($rv);
    
    $rv = 0;
    
    $json = to_json($dataH);
    # print "$json\n";
    # return $rv;
    # exit;
    
    # Call the generic ingest module
    $cmd = "rawingest series=$series segment=$segment";
    if (defined($setdate) && $setdate)
    {
        $cmd = "$cmd -d";
    }

    $pipe = new drmsPipeRun($cmd);
    
    if (defined($pipe))
    {
        $pipe->EnableAutoflush();
        $pipe->WritePipe($json);
        $pipe->ClosePipe(1); # close write pipe
        $pipe->ReadPipe(\$rsp);
        
        if (defined($rsp) && length($rsp) > 0)
        {
            print "rawingest returns:\n$rsp\n";
        }
        
        # close both read and write pipes (but write was already closed)                    
        if ($pipe->ClosePipe())
        {
            print STDERR "Failure reading from pipe.\n";
            $rv = 1;
        }
        else
        {
            $childret = $pipe->GetStatus();
            if ($childret != 0)
            {
                print STDERR "Child process ran unsuccessfully, returned code $childret.\n";
                $rv = 1;
            }
        }
    }
    else
    {
        print STDERR "Unable to call geningest.\n";
        $rv = 1;
    }
    
    return $rv;
}
