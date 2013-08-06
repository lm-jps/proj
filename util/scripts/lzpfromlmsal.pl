#!/home/jsoc/bin/linux_x86_64/activeperl

# This script does not support HTTP requests that require authentication. If such an attempt is 
# made, the web server will respond to the request with a 401 error code, and this script
# will abort.

use strict;
use warnings;

require HTTP::Request;
use LWP::UserAgent;
use Getopt::Long;
use URI::Escape;
use File::Copy;
use File::Path qw(make_path remove_tree);
use Term::ReadKey;

use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use drmsLocks;
use drmsArgs;
# use ParseAutoIndex;

use constant kTmpDir     => "/tmp/lzp";

# Option defaults
use constant kDefBaseURI => "https://aia.lmsal.com/moc/moc_top/sdo_moc/moc/lzp";
use constant kDefAPIDs   => "193,213,283,287,303,307";
use constant kDefBegin   => "2010_040";
use constant kDefEnd     => "2012_040";
use constant kDefDataDir => "/surge40/jsocprod/sdo/mocprods/lzp";

# Return values
use constant kRetSuccess     => 0;
use constant kRetInvalidArgs => 1;
use constant kRetFileIO      => 2;

# Required cmd-line arguments



# Optional cmd-line arguments
use constant kOptBaseuri => "baseuri";
use constant kOptApids   => "apids";
use constant kOptBegin   => "begin";
use constant kOptEnd     => "end";
use constant kOptDataDir => "ddir";

# global password "file"
my(%pwords);

my($baseuri) = &kDefBaseURI;
my($apids) = &kDefAPIDs;
my($begin) = &kDefBegin;
my($end) = &kDefEnd;
my($ddir) = &kDefDataDir;
my($apidsH);
my($reqstr);
my($req);
my($ua);
my($rsp);
my($content);
my($optsinH);
my($opts);
my($doy);
my($yr);
my($rv);


# Optional arguments
$optsinH =
{
    &kOptBaseuri    => 's',
    &kOptApids      => 's',
    &kOptBegin      => 's',
    &kOptEnd        => 's',
    &kOptDataDir    => 's'
};

$opts = new drmsArgs($optsinH, 0);

$rv = &kRetSuccess;

if (!defined($opts))
{
    $rv = &kRetInvalidArgs;
}
else
{
    my($arg);
    my($beginyr);
    my($begindoy);
    my($endyr);
    my($enddoy);
    my($yr);
    my($doy);
    my($doystr);
    my($stopdoy);
    my(@flist);
    my($afile);
    my($tmpdst);
    my($dst);
    my($pword);
    
    
    GetArg($opts, &kOptBaseuri, \$baseuri);
    GetArg($opts, &kOptApids, \$apids);
    GetArg($opts, &kOptBegin, \$begin);
    GetArg($opts, &kOptEnd, \$end);
    GetArg($opts, &kOptDataDir, \$ddir);
    
    # Put apids in hash.
    $apidsH = GetIdHash($apids);
    
    # Parse $begin and $end.
    if (defined($begin) && $begin =~ /^\s*(\d\d\d\d)_(\d\d\d)\s*$/)
    {
        $beginyr = $1;
        $begindoy = $2;
    }
    else
    {
        Usage();
        $rv = &kRetInvalidArgs;
    }
    
    if ($rv == &kRetSuccess)
    {
        if (defined($end) && $end =~ /^\s*(\d\d\d\d)_(\d\d\d)\s*$/)
        {
            $endyr = $1;
            $enddoy = $2;
        }
        else
        {
            Usage();
            $rv = &kRetInvalidArgs;
        }
    }
    
    if ($rv == &kRetSuccess)
    {
        for ($yr = $beginyr; $yr <= $endyr; $yr++)
        {   
            $doy = ($yr == $beginyr ? $begindoy : 1);
            $stopdoy = ($yr == $endyr ? $enddoy : 366);
            
            for (; $doy <= $stopdoy; $doy++)
            {
                if (!defined($ua))
                {
                    $ua = LWP::UserAgent->new();
                    $ua->agent("AgentName/0.1 " . $ua->agent);
                }
                
                $doystr = sprintf("%03d", $doy);                
                $reqstr = "$baseuri/${yr}_$doystr";
                
                # percent-encode the URI string.
                uri_escape($reqstr);
                
                # Obtain a list of files in the directory specified by $reqstr.
                @flist = GetFileList($ua, $reqstr, "${yr}_$doystr", \%pwords);
                
                # Walk through list, and select files that have relevant APIDs in their names.
                # @flist might be empty (e.g., the $yr_$doy directory doesn't exist).
                if ($#flist >= 0)
                {
                    if (!(-e kTmpDir . "/${yr}_$doystr"))
                    {
                        make_path(kTmpDir . "/${yr}_$doystr");
                    }
                    
                    if (!(-d kTmpDir . "/${yr}_$doystr"))
                    {
                        print STDERR "Unable to create path " . kTmpDir . "/${yr}_$doystr";
                        $rv = &kRetFileIO;
                    }
                    
                    if ($rv == &kRetSuccess)
                    {
                        if (!(-e "$ddir/${yr}_$doystr"))
                        {
                            make_path("$ddir/${yr}_$doystr");
                        }
                        
                        if (!(-d "$ddir/${yr}_$doystr"))
                        {
                            print STDERR "Unable to create path $ddir/${yr}_$doystr";
                            $rv = &kRetFileIO;
                        }
                    }
                }
                
                if ($rv == &kRetSuccess)
                {
                    my($apid);
                    
                    foreach $afile (@flist)
                    {
                        # Filter out not-requested files
                        $apid = sprintf("%04d", ($afile =~ /^\s*(\d\d\d\d)_/)[0]);
                        
                        if (!exists($apidsH->{$apid}))
                        {
                            next;
                        }
                        
                        $reqstr = "$baseuri/${yr}_$doystr/$afile";
                        uri_escape($reqstr);
                        $tmpdst = kTmpDir . "/${yr}_$doystr/$afile";
                        $dst = "$ddir/${yr}_$doystr/$afile";
                        
                        unless (DownloadFile($ua, $reqstr, $tmpdst, \%pwords))
                        {
                            print "  Moving $tmpdst to $dst\n";
                            unless (move($tmpdst, $dst))
                            {
                                print STDERR "Unable to save file $dst.\n";
                            }
                        }
                    }
                }
                else
                {
                    last;
                }            
            }
        }
        
        # clean up
        remove_tree(kTmpDir);
    }
    
    exit($rv);
}

sub GetArg
{
    my($optsO) = shift;
    my($argname) = shift;
    my($varout) = shift;
    my($arg);
    
    $arg = $optsO->Get($argname);
    if (defined($arg))
    {
        $$varout = $arg;
    }
}

sub GetIdHash
{
    my($apids) = shift;
    my(@apidsA);
    my($iapid);
    my($apid);
    my($rv);
    
    
    @apidsA = split(/,/, $apids);
    $rv = {};
    
    foreach $iapid (@apidsA)
    {
        $apid = sprintf("%04d", $iapid);
        $rv->{$apid} = 1;
    }
    
    return $rv;
}

sub GetPword
{
    my($reqstr) = shift;
    my($ua) = shift;
    my($pwords) = shift;
    my($server);
    my($usern);
    my($pword);
    my($blob);
    my($rsp);
    
    if ($reqstr =~ /^\s*(https:\/\/[^\/]+)\//)
    {
        $server = $1;
        
        if (exists($pwords->{$server}))
        {
            # This means that we successfully authenticated already (but there could be a time-out).
            # So, do nothing.
        }
        else
        {   
            # Get username.
            $usern = "rbush";
            
            # Get password from user.
            #ReadMode(’noecho’);
            #$pword = ReadLine(0);
            $pword = "blah";
            
            $pwords->{$server} = "$usern:'$pword'";
        }
    }
    
    # Now actually authenticate with the server.
    $rsp = $ua->post("$server/auth-service", 
                      'cmd' => 'submit-login',
                      'j_username' => $usern,
    'j_password' => $pword);
    #'callback-url' => $reqstr);
    
    if ($rsp->is_error) 
    {
        print "Failed to authenticate: " . $rsp->status_line . "\n";
    }
    else
    {
        print "content : " . $rsp->decoded_content . "\n";
    }
}

sub GetFileList
{
    my($ua) = shift;
    my($reqstr) = shift;
    my($date) = shift;
    my($pwords) = shift;
    my($req);
    my($content);
    my(@contentA);
    my($fname);
    my($rsp);
    my(@rv);
    
    #    my($infoH);
    #    my($fileinfoH);
    #    my($fname);
    
    GetPword($reqstr, $ua, $pwords);
    
    $req = HTTP::Request->new(GET => $reqstr);
    $req->content_type('application/x-www-form-urlencoded');
    
    # Send the request
    $rsp = $ua->request($req);
    
    if ($rsp->is_error)
    {
        $content = $rsp->status_line;
        print STDERR "ERROR ($reqstr) : $content\n";
    }
    else
    {
        $content = $rsp->decoded_content;
    }
    
        exit;
    
    # Parse out files from response.
    # Don't use parse_index() - it doesn't work.
    # $infoH = ParseAutoIndex::parse_index($content);
    #foreach my $file (keys %$infoH)
    #{
    #    $fileinfoH = $infoH->{$file};
    #    $fname = $fileinfoH->{name};
    #    print "pushing $fname\n";
    #    pushd(@rv, $fname);
    #
    # }
    
    @contentA = split(/\n/, $content);
    foreach my $line (@contentA)
    {
        chomp($line);
        if ($line =~ /(\d\d\d\d_${date}_\d\d\.hkt.xml)/)
        {
            $fname = $1;
        }
        elsif ($line =~ /(\d\d\d\d_${date}_\d\d\.hkt)/)
        {
            $fname = $1;
        }
        else
        {
            $fname = "";
        }
        
        if (length($fname) > 0)
        {
            push(@rv, $fname);
        }
    }
    
    return @rv;
}

# returns 1 on failure
sub DownloadFile
{
    my($ua) = shift;
    my($reqstr) = shift;
    my($dest) = shift;
    my($pwords) = shift;
    my($req);
    my($rsp);
    my($content);
    my($rv);
    
    $rv = 0;
    
    
    GetPword($reqstr, $ua, $pwords);

    $req = HTTP::Request->new(GET => $reqstr);
    $req->content_type('application/x-www-form-urlencoded');
    
    # Send the request                                                                     
    $rsp = $ua->request($req);
    
    if ($rsp->is_error)
    {
        $content = $rsp->status_line;
        print STDERR "Unable to download file $reqstr.\n";
        print STDERR "Response: $content\n";
        $rv = 1;
    }
    else
    {                                                   
        $content = $rsp->decoded_content;
        
        if (open(DST, ">$dest"))
        {
            print "Downloading $reqstr to $dest.\n";
            print DST "$content";
            unless (close(DST))
            {
                print STDERR "An error occurred while writing file $dest.\n";
                $rv = 1;
            }
        }
        else
        {
            print STDERR "Unable to open file $dest for writing\n";
            $rv = 1;
        }
    }
    
    return $rv;
}

sub Usage
{
    print "Usage:\n";
    print "lzpfromlmsal.pl [ --baseuri=<base uri> ] [ --ddir=<save directory> ] [ --apids=<apid list> ] [ --begin=<first day> ] [ --end=<last day> ]\n";
    print "  <base uri> - URI of location containing YYYY_DDD directories.\n";
    print "  <save directory> - path to location to save product files.\n";
    print "  <apid list> - comma-separated list of APIDs of products to download.\n";
    print "  <first day> - YYYY_DDD date of first day of date range of products to download.\n";
    print "  <last day> - YYYY_DDD date of last day of date range of products to download.\n";
}

