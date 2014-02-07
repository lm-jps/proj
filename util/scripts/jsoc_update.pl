#!/home/jsoc/bin/linux_x86_64/activeperl -w  


use File::stat;
use Fcntl qw(:mode);
use Getopt::Long;
use FileHandle;
use Sys::Hostname;
use Cwd qw(realpath chdir); # OMG! Need to override chdir, otherwise $ENV{PWD} is NOT
                            # updated when chdir is called.
use FindBin qw($Bin);
use lib ($Bin, "$Bin/../../../base/libs/perl");
use drmsArgs;

use constant kSuccess       => 0;
use constant kBadArgs       => 1;
use constant kCantSSH       => 2;
use constant kConflicts     => 3;
use constant kFileStatus    => 4;
use constant kCantExeOnMach => 5;
use constant kCantUpdateSrc => 6;


use constant kUpdateLog => "cvsupdate.log";
use constant kMtabFile  => "/etc/mtab";
use constant kTmpDir    => "/tmp";

my(@machines) = 
(
    "n02",
    "solar3"
);

my($optsinH);
my($opts);
my($optdir);
my($machlist);
my($cvstree);
my($jsocroot);
my($nfsf);
my($volume);
my($mount);
my($netpath);
my($mach);
my($hostnm);


# If the caller is attempting to update a CVS tree that resides on a local disk, then 
# it MIGHT not be possible to update that tree after ssh'ing to a sanctioned machine
# (e.g., n02). Or the caller might be attempting to update a CVS tree that resides
# on a network disk, but the disk might not be exported to a sanctioned machine.
# This script checks for those problems.
#
# Determine whether the cvs tree resides on a local disk, or if it is NFS-
# mounted. If it is on local disk, then create a network volume name for the disk (e.g.,
# <machine>:<path>). If the tree is NFS-mounted, then use the mtab to create a
# network volume.
# 
# If the CVS tree does not reside on a disk accessible from the sanctioned machines,
# ask the caller if they want to build the updated tree on the machine running this
# script.

$optsinH = 
{
    "dir" => 's',
    "machs" => 's'
};

$opts = new drmsArgs($optsinH, 0);

if (!defined($opts))
{
    exit(kBadArgs);
}

$optdir = $opts->Get("dir");
$machlist = $opts->Get("machs");

# If user did not specify a directory to update, use the one specified by the JSOCROOT
# env variable, if it exists.
if (!defined($optdir))
{
   $jsocroot = $ENV{'JSOCROOT'};
   (!defined($jsocroot)) ? $cvstree = $ENV{'PWD'} : $cvstree = $jsocroot;
}
else
{
    # $cvstree needs to be an absolute path. If it isn't, make it so.
    if ($optdir =~ /^\s*\//)
    {
        $cvstree = $optdir;        
    }
    else
    {
        $cvstree = "$ENV{'PWD'}/$optdir"
    }
}

if (defined($machlist))
{
    my(@mintermed) = split(/,/, $machlist);
    @machines = map({ ($_ =~ /\s*(\S+)\s*/)[0] } @mintermed);
    
}

# Determine if the cvs tree specified exists
if (!(-d $cvstree))
{
   print STDERR "Specified CVS tree does not exist.\n";
   exit(kBadArgs);
}

# Determine if the cvs tree resides on a local disk.
# 1. Obtain the device number for each record in /etc/mtab.
# 2. Walk down this list, comparing the device on which the CVS tree resides with each
#    record in this list. If there is a match, then the CVS tree resides on a mounted directory.
#    Check the matching mtab record for an NFS file system. If not, then assume
#    the CVS tree is on a local disk. Otherwise, it is NFS-mounted, and we should
#    obtain the hosting volume and mount point from the mtab record.
$nfsf = (GetNFSInfo($cvstree, \$volume, \$mount));

if (!$nfsf)
{
   # Assume the tree is on a local disk. Create a 'network name' for the tree.
   $hostnm = hostname();
   $netpath = "$hostnm:" . realpath($cvstree);
}
else
{
   $netpath = realpath($cvstree);
   $netpath =~ s/$mount/$volume/;
}

# Take the local path of the JSOC tree, identify the mount root in this path
# (the mount root is the part of the path that resides in the appropriate mtab 
# record), then substitute the mount-path part of the local path with the 
# network volume (e.g., <machine>:<path>). So if the JSOC tree's local path
# is /auto/home1/arta/jsoctrees/JSOC, and the mtab record for the network
# drive that contains this tree is "sunroom:/home1 /auto/home1 ..." 
# substitute "/auto/home1" in the JSOC-tree path with "sunroom:/home1" to
# form the network path "sunroom:/home1/arta/jsoctrees/JSOC".

print "CVS JSOC tree being updated ==> $netpath\n";

# Build the JSOC tree on the sanctioned machines, if it is mounted.

# Anonymous hash whose keys are machines. For each key, the value is the 
# path, local to the machine, to the CVS tree being updated.
my($mpaths) = GetLocalPaths($netpath, @machines);
my(@mpkeys) = keys(%$mpaths);
my($impkey);
my($line);
my($cmd);

if ($#mpkeys == -1)
{
   $hostnm = hostname();
   print "The CVS tree $netpath is not mounted on any of the requested build machines (i.e., n02).\n";
   print "Would you like to build on the machine on which this script is running, $hostnm? (y/n)\n";
   $line = <STDIN>;
   chomp($line);
   $line = lc($line);

   if ($line =~ /^y/)
   {
      $mpaths->{$hostnm} = $cvstree;
      push(@mpkeys, $hostnm);
      print "                 Local path ==> $hostnm:$cvstree.\n";
   }
}

if ($#mpkeys >= 0)
{
    my($cdir) = realpath($ENV{'PWD'});
    my(@rsp);
    my($rspstr);
    
    print "\n";
    print "#############################################\n";
    print "####### Starting CVS update #################\n";
    print "#############################################\n";
    print "##\n";
    print "## changing working directory to $cvstree.\n";
    
    chdir($cvstree);
    unless (UpdateTree($cvstree))
    {
        # Run the configure script - should really check the return code.
        print "####### Running configure script ############\n";
        @rsp = `./configure 2>&1`;
        print "## Done (output from configure follows).\n";
        print "## ";
        $rspstr = join("## ", @rsp);
        print $rspstr;
        
        # Do the build on each machine
        print "####### Building binaries on machines #####\n";
        map({ BuildOnMachine($_, $mpaths); } @mpkeys);
        print "## Done building.\n";
        print "##\n";
    }
    else
    {
        # Failure updating CVS tree.
        print "## restoring working directory to $cdir.\n";
        chdir($cdir);
        
        print "###########################################\n";
        print "####### Update of CVS binaries FAILED!! ###\n";
        print "###########################################\n";
        
        exit(kCantUpdateSrc);
    }
    
    print "## restoring working directory to $cdir.\n";
    chdir($cdir);
    
    print "###########################################\n";
    print "####### Update of CVS binaries complete. ##\n";
    print "###########################################\n";
}

exit(kSuccess);

#my($st) = stat("/tmp"); # 2050
#my($st) = stat("/SUM12"); # 33 (must refer to the physical disk device on d02)
#my($st) = stat("/auto/SUM12"); # 33
#my($st) = stat("/dev/sda2"); # 2050
#my($st) = stat("/auto/home1");
#print "isblock\n" if (S_ISBLK($st->mode));
#(S_ISBLK($st->mode)) ? print $st->rdev . "\n" :  print $st->dev . "\n";

sub Usage
{
   print "jsoc_update.pl\n";
   print "  If the environment variable JSOCROOT is set, update the CVS tree specified by\n";
   print "  \$JSOCROOT. Otherwise, if the current directory is a valid DRMS CVS tree, update\n";
   print "  the current directory. If the current directory is not a valid DRMS CVS tree, exit\n";
   print "  without updating anything.\n";

   print "jsoc_update -d <path>";
   print "  Update the DRMS CVS tree specified by <path>.\n";
}

# Returns 1 and the network volume containing the JSOC tree (if the JSOC tree is
# NFS-mounted). Returns 0 otherwise.
sub GetNFSInfo
{
   my($tree) = $_[0];
   my($volr) = $_[1];
   my($mountr) = $_[2];
   
   my($rv) = 0;
   my(@mtab);
   my($fh);

   # Get the device number for the CVS tree.
   my($st) = stat($tree);
   my($devno) = (S_ISBLK($st->mode)) ? $st->rdev : $st->dev; # number of device containing 
                                                             # JSOC tree.
   my($info);

   # Open mtab file
   $fh = FileHandle->new("<" . kMtabFile);

   if (defined($fh))
   {
      @mtab = <$fh>;
      
      # Not sure how to short-circuit the map function (we don't want to examine every line
      # in the mtab if we found a relevant one).
      map({ ProcessMtabInfo(\$info, $devno, $_); } @mtab);
      $fh->close();
   }
   else
   {
      print STDERR "Fatal error - cannot find mtab.\n";
   }

   if (defined($info))
   {
      $rv = 1;
      $$volr = $info->{'vol'};
      $$mountr = $info->{'mount'};
   }

   return $rv;
}

# Find the network drive that contains the JSOC tree.
sub ProcessMtabInfo
{
    my($infor) = $_[0];
    my($devno) = $_[1]; # number of the device containing the JSOC tree
    my($line) = $_[2];
    
    my($vol);
    my($mount);
    my($type);
    my($st);
    my($mtabdevno);
    
    chomp($line);
    
    # 1 - NFS volume (machine:volume)
    # 2 - mount point
    # 3 - partition type
    if ($line =~ /\s*(\S+)\s*(\S+)\s*(\S+)/)
    {
        $vol = $1;
        $mount = $2;
        $type = $3;
        
        # Look up the device number for this line
        $st = stat($mount);
        if (defined($st))
        {
            $mtabdevno = (S_ISBLK($st->mode)) ? $st->rdev : $st->dev;
            
            if ($devno == $mtabdevno && $type =~ /nfs/i)
            {
                # match - we have the device on which the CVS tree resides, and that device
                # is NFS-mounted.
                if (!(defined($$infor)))
                {
                    # New empty hash.
                    $$infor = {};
                }
                
                # Save the network drive that contains the JSOC tree (and where it is mounted on
                # the current machine).
                $$infor->{'vol'} = $vol;
                $$infor->{'mount'} = $mount;
            }
        }
    }
}

# $line is the remote machine's (e.g., n02) mtab line.
sub ExtractMachPath
{
    my($netpath) = $_[0];
    my($line) = $_[1];
    
    my($volume);
    my($mount);
    
    # First, extract the volume from the mtab record
    chomp($line);
    
    if ($line =~ /^\s*(\S+)\s+(\S+)/)
    {
        $volume = $1;
        $mount = $2;
        
        # $volume - like sunroom:/home0
        # $1 - like /home0
        # $netpath - like sunroomg:/home0/arta/jsoctrees/JSOC
        if ($netpath =~ /^$volume(.+)/)
        {
            return $mount . $1;
        }
        else
        {
            # There might be a trailing 'g' on the $netpath (legacy - used to mean
            # gigabit for a gigabit network). But all networks are at least Gb now.
            # One physical drive might have a network name whose symbolic name ends
            # in a 'g' on one machine, but not on another. Try to remove the trailing
            # 'g' before doing a comparison. The mtab line might also have this trailing
            # 'g'.
            #
            # $netpath2 - like sunroom:/home0/arta/jsoctrees/JSOC
            my($netpath2) = $netpath;
            
            $netpath2 =~ s/g:/:/;
            $volume =~ s/g:/:/;
            
            if ($netpath2 =~ /^$volume(.+)/)
            {
                return $mount . $1; 
            }
            else
            {
                return ();
            }
        }
    }
    else
    {
        print STDERR "Unexpected response from server: $line\n";
        return ();
    }
}

sub GetLocalPaths
{
    my($netpath) = shift; # network path to CVS tree
    my(@machines) = @_;
    
    my($rv);
    my(@machpath);
    
    foreach $mach (@machines)
    {
        # make a test connection
        if (open(STATCMD, "(ssh $mach cat /etc/mtab) 2>&1 |"))
        {
            my(@remotemtab) = <STATCMD>;
            
            close(STATCMD);
            
            # Extract record for $volume.
            @machpath = map({ ExtractMachPath($netpath, $_); } @remotemtab);
            
            my($test) = $#machpath;
            
            if ($#machpath == 0)
            {
                # CVS tree lives on $mach at $machpath[0];
                print "                 Local path ==> $mach:$machpath[0].\n";
                
                if (!defined($rv))
                {
                    $rv = {};
                }
                
                $rv->{$mach} = $machpath[0];
            }
            else
            {
                print "CVS tree $netpath is not mounted on $mach; skipping machine.\n";
                next;
            }
        }
        else
        {
            print "Cannot ssh to $mach.\n"; 
            exit(kCantSSH);
        }
    }
    
    return $rv;
}

sub UpdateTree
{
   my($cvstree) = $_[0];

   my($rv) = 0;
   my(@statchk);
   my($rsp);
    my(@splitRes);

   print "####### Downloading repository changes ######\n";

   # Update local CVS tree with CVS-repository changes.
   $cmd = "$Bin/jsoc_sync.pl -l" . kUpdateLog; # use the jsoc_sync.pl relative to this script
   print "## Updating source files in $cvstree (calling '$cmd').\n";
   $rsp = qx($cmd 2>&1);
    
    if ($? >> 8)
    {
        # Error calling jsoc_sync.pl.
        print STDERR "## Failure calling jsoc_sync.pl\n";
    }
    
    @splitRes = split(qr/\n/, $rsp);
    $rsp = "";
    foreach my $line (@splitRes)
    {
        $rsp = $rsp . "## $line\n";
    }
    
   print "$rsp";
   print "##\n";

    # Log is in parent directory.
   if (open(UDL, "<" . "../" . kUpdateLog))
   {
      my(@content) = <UDL>;
      my(@conflicts) = grep({ ($_ =~ /^\s*C/) ? $_ : () } @content);
      my($iline);

      close(UDL);

      if ($#conflicts >= 0)
      {
         print "## A list of locally modified files with changes that conflict with repository changes follows:\n";
         foreach $iline (@conflicts)
         {
            print "## $iline\n";
         }

         print "## Please resolve these conflicts, then update again.\n";
         print "## (see release notes to learn how to deal with conflicting file changes.)\n";
         exit(kConflicts);
      }

      print "####### Checking file status ################\n";
      @content = CheckStatus();
      print "## Done.\n";
      
      if ($#content >= 0)
      {
         print "## A list of local files that differ from their repository counterparts follows:\n";
         foreach $iline (@content)
         {
            chomp($iline);
            print "## $iline\n";
         }

         print "##\n";
      }

      print "## All checks succeeded. Continue by entering 'cont'.\n";
      print "## ";
      $line = <STDIN>;
      chomp($line);
      $line = lc($line);

      if ($line !~ /^\s*cont/)
      {
         print "## Update aborted by user.\n";
         $rv = 1;
      }
   }
   else
   {
      # Couldn't open update log - bail.
      print "## Could not open update log file ../" . kUpdateLog . ".\n";
      $rv = 1;
   }

   print "##\n";

   return $rv;
}

sub CheckStatus
{
   my(@res);

   my $PID = getppid;
   my $user = $ENV{'USER'};
   my $logfile = "/tmp/cvs_status_".$user."_$PID.log";
   `cvs status 1> $logfile 2>&1`;
   open(CV, $logfile) || die "Can't open $logfile: $!\n";

   while (<CV>)
   {
      #print "$_";		#!!!TEMP
      if (/^File:/)
      {
         if (/Up-to-date/)
         {
            next;
         }
         push(@res, "$_");
         <CV>;
         <CV>;
         $_ = <CV>;             #get Repository revision line
         s/\s+/ /g;             #compress multiple spaces
         push(@res, "$_\n");
      }
   }

   return @res;
}

sub BuildOnMachine
{
   my($mach) = $_[0];
   my($mpathsr) = $_[1];

   my($plat);
   my($lpath) = $mpathsr->{$mach};

   if (defined($lpath))
   {
      # Should really capture the return code from this ssh cmd.
      if (open(CMD, "(ssh $mach " . "\'" . "echo \$JSOC_MACHINE" . "\'" . ") 2>&1 |"))
      {
         $plat = <CMD>;
         chomp($plat);
         close(CMD);

         print "## Starting build on platform $plat.\n";

         # Should really capture return code from this ssh cmd.
         if (open(CMD, "(ssh $mach 'cd $lpath; /home/jsoc/make_jsoc.pl') 1>make_jsoc_$plat.log 2>&1 |"))
         {
            close(CMD);
         } 
         else
         {
            print "## Unable to run cmd on $mach [1].\n";
         }

         print "## Build on platform $plat complete.\n";
      } 
      else
      {
         print "## Unable to run cmd on $mach [2].\n";
      }
   }
}
