#!/home/jsoc/bin/linux_x86_64/perl5.12.2                                                       

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
use constant kArgKeyNames => "keynames";  # Full path to file containing key names
use constant kArgKeyVals  => "keyvals";   # Full path to file containing key values
use constant kArgBitMaps  => "bmaps";     # Path to bit maps
use constant kBatchSize   => 128;

my($argsinH);
my(@args);
my($kevals);
my($keynames);
my($keyvals);
my($bmaps);
my($dataH);
my($cpkeyH);
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
my($rv);

$rv = 0; # success

# Collect arguments
$argsinH =
{
    &kArgKeyNames =>   's',
    &kArgKeyVals =>    's',
    &kArgBitMaps =>    's'
};

@args = GetArgs($argsinH);

if (@args)
{
    ($keynames, $keyvals, $bmaps) = @args;
    
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
                $ival = 0;
                foreach my $aval (@vals)
                {
                    if ($ival < $npkeys)
                    {
                        # These are prime-key keyword values.
                        push(@pkeyvals, $aval);
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

                @pkeyvals = ();
                $iline++;
                
                # Call ingestion module with batches of records.
                if ($iline % &kBatchSize == 0)
                {
                    $rv = Ingest($dataH);
                    if ($rv != 0)
                    {
                        last;
                    }
                    
                    $dataH = {};
                    $cpkeyH = undef;
                }
            } # End while
            
            if ($rv == 0)
            {
                if (keys(%$dataH) > 0 && defined($cpkeyH) && keys(%$cpkeyH) > 0)
                {
                    $rv = Ingest($dataH);
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
                }
            }
        }
    }
    
    return @rv;
}

sub Ingest
{
    my($dataH) = @_;
    my($json);
    my($pipe);
    my($rsp);
    my($childret);
    my($rv);
    
    $rv = 0;
    
    $json = to_json($dataH);
    print "$json\n";
    exit;
    
    # Call the generic ingest module
    $pipe = new drmsPipeRun("geningest");
    
    if (defined($pipe))
    {
        $pipe->EnableAutoflush();
        $pipe->WritePipe($json);
        $pipe->ClosePipe(1); # close write pipe
        $pipe->ReadPipe(\$rsp);
        
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
