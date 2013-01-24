#!/home/jsoc/bin/linux_x86_64/perl5.12.2

# This script maintains the 'production' set of JSOC binaries, which are located in 
# /home/jsoc/cvs/Development/JSOC/bin. It works as follows:
# 
# 1. It runs "make MACH='waystation_<arch>' <target1> <target2> ..." This will create binary files in 
#    /home/jsoc/cvs/Development/JSOC/_waystation_<arch>. There will be links created from 
#    /home/jsoc/cvs/Development/JSOC/bin/waystation_<arch> and 
#    /home/jsoc/cvs/Development/JSOC/lib/waystation_<arch> to the just-created binaries.
# 2. It then follows all these links to identify all the newly created files. Each of these files
#    is going to replace a file in /home/jsoc/cvs/Development/JSOC/_<arch>. However, these original
#    files may be in use, so each is first moved to a "save" directory, which is named
#    /home/jsoc/cvs/Development/JSOC/_<arch>_YYYYMMDDHHMMSS. This save directory
#    must be created every time this script is run. At regular intervals, these save directories 
#    must be deleted (a one-day retention time is good). Each original file is first moved to the save
#    directory, then its replacement (in /home/jsoc/cvs/Development/JSOC/_waystation_<arch>)
#    is moved into place in the production directory, /home/jsoc/cvs/Development/JSOC/_<arch>. 
# 3. The links to the files in /home/jsoc/cvs/Development/JSOC/_waystation_<arch> will now be dead.
#    This script therefore deletes those links.
# 4. /home/jsoc/cvs/Development/JSOC/waystation_<arch> is now the repository for all the object and dependency files.
#    make will use these files when it builds binaries, so we must maintain this directory. 

use strict;
use warnings;
use Cwd qw(getcwd realpath chdir);
use File::Spec;
use File::stat;
use IO::Dir;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use drmsLocks;
use drmsArgs;
use drmsRunProg;

# Required arguments
use constant kArgMods      => "mods";   # modules to build (comma-separated list)
use constant kArgFiles     => "files";  # source files to first update (comma-separated list)

# Optional arguments
# If clean is specified, then the required arguments above are ignored. In this case, the $bindir is cleaned,
# just like 'make clean' will do. And $builddir is cleaned too.
use constant kArgClean     => "clean";  # clean $bindir (i.e., _linux_x86_64) AND $builddir (i.e., _waystation_$bindir)
use constant kArgBinDir    => "bindir"; # _waystation_<value> is the directory where the binaries are saved 
                                        # (defaults to _waystation_$JSOC_MACHINE)

# Lockfile
use constant kLockFile     => "/home/jsoc/locks/prodbuildlck.txt";
use constant kJSOCtreeRoot => "JSOC";

# Return values
use constant kRetSuccess       => 0;
use constant kRetNoLock        => 1;
use constant kRetInvalidArgs   => 2;
use constant kRetDlSource      => 3;
use constant kRetInvalidWD     => 4;
use constant kRetConfigure     => 5;
use constant kRetMake          => 6;
use constant kRetMove          => 7;
use constant kRetDelLinks      => 8;

my($rv);
my($argsinH);
my($optsinH);
my($lock);
my($args);
my($opts);
my(@mods);
my($bindir);
my($cmd);
my($basetime);
my(@srcdirs);
my($wdir); # working dir (which must be a JSOC root directory)
my(@relspecs);
my($specstr);
my($usewaystation);
my($builddir);
my($clean);

$rv = &kRetSuccess;

# Required arguments
$argsinH =
{
    &kArgMods    => 's',
    &kArgFiles   => 's',
};

$optsinH =
{
    &kArgClean   => 'noval',
    &kArgBinDir  => 's'
};

# Lock this script
$lock = new drmsNetLocks(&kLockFile);

if (defined($lock))
{
    $rv = &kRetSuccess;
    
    $basetime = time();

    $opts = new drmsArgs($optsinH, 0);
    
    $clean = 0;
    if (defined($opts))
    {
        $clean = $opts->Get(&kArgClean);
    }
    else
    {
        $clean = 0;
    }
    
    if ($clean)
    {
        $bindir = $opts->Get(&kArgBinDir);
        if (!defined($bindir))
        {
            $bindir = $ENV{JSOC_MACHINE};
        }
        
        if (!defined($bindir))
        {
            print STDERR "Binary directory not specified.\n";
            $rv = &kRetInvalidArgs;
        }
        else
        {
            $cmd = "make 'MACH=$bindir' clean";            
            print "running $cmd\n";
            
            if (drmsSysRun::RunCmd($cmd) != 0)
            {
                print STDERR "Unable to clean $bindir.\n";
                $rv = &kRetMake;
            }
            else
            {
                $builddir = "waystation_$bindir";
                $cmd = "make 'MACH=$builddir' clean";
                print "runnning $cmd\n";
                
                if (drmsSysRun::RunCmd($cmd) != 0)
                {
                    print STDERR "Unable to clean $builddir.\n";
                    $rv = &kRetMake;
                }
            }
        }
    }
    else
    {
        $args = new drmsArgs($argsinH, 1);
        
        if (!defined($args))
        {
            $rv = &kRetInvalidArgs;
        }
        else
        {
            # kArgMods must be defined (otherwise new drmsArgs would have failed).
            @mods = split(qr(,), $args->Get(&kArgMods));
            @relspecs = split(qr(,), $args->Get(&kArgFiles));
            $bindir = $opts->Get(&kArgBinDir);
            if (!defined($bindir))
            {
                $bindir = $ENV{JSOC_MACHINE};
            }
            
            if (!defined($bindir))
            {
                print STDERR "Binary directory not specified.\n";
                $rv = &kRetInvalidArgs;
            }
        }

        if ($rv == &kRetSuccess)
        {
            $wdir = getcwd();
            
            # Make sure that the current directory is a JSOC-tree root.
            # I'm probably going to regret this, but to ensure that the working directory is the root of
            # a JSOC tree, let's consider this a valid tree if it contains base/jsoc_version.h. Let's not
            # count of the current directory being "JSOC" (even though this is a requirement, I know
            # Rick renames his tree root directory to DRMS).
            if (!(-f "$wdir/base/jsoc_version.h"))
            {
                print STDERR "The current directory is the root of a valid JSOC source-code tree.\n";
                $rv = kRetInvalidWD;
            }
        }
        
        if ($rv == &kRetSuccess)
        { 

            # Update source files.
            $specstr = join(',', @relspecs);
            
            $cmd = "/home/jsoc/dlsource.pl -o update -s $specstr";
            if (drmsSysRun::RunCmd($cmd) != 0)
            {
                print STDERR "Unable to update source files.\n";
                $rv = &kRetDlSource;
            }
        }
        
        # Build modules.
        if ($rv == &kRetSuccess)
        {
            my($targetstr);
            
            $targetstr = join(' ', @mods);
            
            # Run configure, just to be sure (but this assurance, which is necessary, will trigger a 
            # rebuild of pretty much all binaries since it, in essence, updates the timestamps on all
            # headers. I need to make the configure script smarter so that it doesn't re-create
            # all header links. Instead it should simply make a link to a header if the link does 
            # not already exist).
            $cmd = "./configure";
            if (drmsSysRun::RunCmd($cmd) != 0)
            {
                print STDERR "'configure' falied to run properly.\n";
                $rv = &kRetConfigure;
            }
            
            # Run make. If the $bindir already exists, then we need to build in waystation_$bindir, and
            # then copy the resulting binaries back to $bindir. But if $bindir does NOT exist, then 
            # we simply build in $bindir.
            
            if ($rv == &kRetSuccess)
            {
                print "bindir is $bindir\n";
                if (-d "_$bindir")
                {
                    $usewaystation = 1;
                    $builddir = "waystation_$bindir";
                }
                else
                {
                    $usewaystation = 0;
                    $builddir = $bindir;
                }
            }
            
            if ($rv == &kRetSuccess)
            {
                $cmd = "make 'MACH=$builddir' $targetstr";
                print "runnning $cmd\n";
                
                if (drmsSysRun::RunCmd($cmd) != 0)
                {
                    print STDERR "Unable to update source files.\n";
                    $rv = &kRetMake;
                }
            }
        }
        
        exit;
        
        # Move the old binaries to a save location, and move the new binaries into locations vacated by
        # the old binaries.
        if ($rv == &kRetSuccess)
        {
            # $basetime is the time (number of secs since the epoch) when the script started running. 
            # For each binary that was created by this script, the binary that will be replaced
            # by this new binary needs to be MOVED into a save directory. Then, the new binary
            # needs to be MOVED into the production tree.
            
            # Binaries whose timestamps are newer than $basetime are considered binaries that 
            # were created by this script.
            push(@srcdirs, "$wdir/bin/waystation_$bindir");
            
            if (MoveFiles($basetime, \@srcdirs) != 0)
            {
                print STDERR "Unable to move newly created binaries into place.\n";
                $rv = &kRetMove;
            }
            
            if ($rv == &kRetSuccess)
            {
                # Remove the links from bin/$builddir to _$builddir.
                if (RemoveLinks() != 0)
                {
                    print STDERR "Unable to remove dead links.\n";
                    $rv = &kRetDelLinks;
                }
            }
        }
    }
}
else
{
    print STDERR "This script is already running; bailing out.\n";
    $rv = &kRetNoLock;
}

exit $rv;


sub FormSpecs
{
    my($wdir) = $_[0];
    my($spec) = $_[1];
    my($fspec);
    my(@rv) = ();
    
    $fspec = File::Spec->catfile($wdir, $spec);
    
    if (-e $fspec)
    {
        # Strip off the JSOC root prefix.
        if ($fspec =~ /.+\/JSOC\/(.+)/) # 'greedy' algorithm
        {
            push(@rv, $1);            
        }
    }
    else
    {
        print STDERR "Warning: Invalid file specification '$spec'; skipping.\n";
    }
    
    return @rv;
}

sub FilterFile
{
    my($afile) = shift;
    
    if ($afile !~ /^\.$/ && $afile !~ /^\.\.$/ && $afile !~ /\.o\.d$/ && $afile !~ /\.o$/)
    {
        return $afile;
    }
    else
    {
        return ();
    }
}

sub GetNewFiles
{
    my($basetime) = shift;
    my($dir) = shift;
    my(@files) = @_;
    
    my($timestamp);
    my($tinfo);
    my(@newfiles);
    my($realfile);
    
    @files = map(FilterFile($_), @files);
    
    foreach my $afile (@files)
    {
        # Resolve links
        $realfile = realpath("$dir/$afile");
        
        # Now skip files that were not just created.
        $tinfo = stat($realfile);
        
        if (!$tinfo)
        {
            print STDERR "Unable to stat file $realfile.\n";
            next;
        }
        
        $timestamp = $tinfo->mtime;
        if ($timestamp < $basetime)
        {
            next;
        }
        
        if (-f $realfile)
        {
            # Not a subdirectory.
            push(@newfiles, $realfile);
        }
        elsif (-d $realfile)
        {
            # A subdirectory.
            my(@dirfiles);
            
            tie(my(%dirH), "IO::Dir", "$realfile");
            @dirfiles = keys(%dirH);
            push(@newfiles, GetNewFiles($basetime, "$realfile", @dirfiles));
            untie(%dirH);
        }
        else
        {
            print STDERR "Invalid file type, file $realfile.\n";
        }
    }
    
    return @newfiles;
}

sub GetOldFile
{
    my($afile) = shift;
    my($oldfile);
    
    if ($afile =~ /\/_waystation_\S+\//)
    {
        $oldfile = $afile;
        $oldfile =~ s/_waystation//;
        return $oldfile;
    }
    else
    {
        print STDERR "Unsupported source file $afile\n";
        return "";
    }
}

# Move newly built files from the _waystation_<arch> directory to the _<arch> directory.
# First, move the to-be-replaced files in the _<arch> directory to the save directory.
sub MoveFiles
{
    my($basetime) = $_[0];
    my($srcdirsR) = $_[1];
    
    my($rv) = 1;
    my($savedir);
    my(%tree);
    my(@allfiles);
    my(@newfiles);
    my(@oldfiles);
    my($srcfile);
    my($tgtfile);
    my($timestamp);
    
    # Collect a list of binary files that were created during the build phase.
    foreach my $dir (@$srcdirsR)
    {
        if (-d $dir)
        {
            tie(%tree, "IO::Dir", $dir);
            
            # @allfiles contains base file names (not paths).
            @allfiles = keys(%tree);
            # @newfiles contains full paths.
            push(@newfiles, GetNewFiles($basetime, $dir, @allfiles));
            untie(%tree);
        }
        else
        {
            print STDERR "Binary directory $dir does not exist.\n";
            last;
        }
    }
    
    # Find the original files that these new files will replace
    foreach $srcfile (@newfiles)
    {
        $tgtfile = GetOldFile($srcfile);
        print "new file $srcfile, old file $tgtfile\n";
        
        # Move old file to save directory.
# ART
    }

    
    # for now, print new files
    exit;
    
    
    # Must create link from JSOC/bin/$bindir to each binary installed, and from JSOC/lib/$bindir to each
    # library installed.
    return $rv;
}
