#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

# This script checks for the presence of vital data-processing modules in the 'production' code trees.

use XML::Simple;
use Data::Dumper;
use Fcntl ':mode';

my($xml);
my($xmlobj) = new XML::Simple;
my($xmldata);
my($environame);
my($enviroroot);
my($projname);
my(%projowner);
my(%failedcheck);
my($filename);
my($filepath);
my($fullpath);
my($fileperms);
my($actmode);
my($desmode);
my($msg);
my(@projwerr);
my($err);

$err = 0;

# Read in the configuration file to obtain the set of project files that will reside
# in the custom checkout set.
if (!ReadCfg(\$xml) && defined($xml))
{
   $xmldata = $xmlobj->XMLin($xml, ForceArray => 1);
   # print Dumper($xmldata);
}
else
{
   print STDERR "Unable to read or parse configuration.\n";
   $err = 1;
}

if (!$err)
{
   # Traverse tree, ensuring that all required files are present
   foreach $enviro (@{$xmldata->{'enviro'}})
   {
      $environame = $enviro->{'name'}->[0];
      $enviroroot = $enviro->{'root'}->[0];

      print "Processing environment '$environame' at $enviroroot...\n";

      foreach $proj (@{$enviro->{'projects'}->[0]->{'proj'}})
      {
         $projname = $proj->{'name'}->[0];
      
         print "  processing project '$projname'\n";

         $projowner{$projname} = {};

         foreach $owner (@{$proj->{'owners'}->[0]->{'owner'}})
         {
            $projowner{$projname}->{$owner->{'content'}} = 
            {
             'login' => $owner->{'login'}, 'email' => $owner->{'email'}};
         }

         foreach $file (@{$proj->{'files'}->[0]->{'file'}})
         {
            $filename = $file->{'name'}->[0];
            $filepath = $file->{'path'}->[0];
            $fileperms = $file->{'perms'}->[0];

            $fullpath = "$enviroroot/$filepath/$filename";

            # Now, actually check for presence of file with correct permissions.
            print "    checking for existence of $fullpath with permissions $fileperms...";

            $msg = "";

            if (-e $fullpath)
            {
               $actmode = (lstat($fullpath))[2] & 07777; # decimal string
               $desmode = oct($fileperms); # decimal string

               # compare modes by anding
               if (($actmode & $desmode) != $desmode)
               {
                  $msg = $msg . sprintf("Required file '$fullpath' has wrong permissions: actual %o, required %o.\n", $actmode, $desmode);
               }
            }
            else
            {
               $msg = $msg . "Required file '$fullpath' does not exist.\n";
            }

            if (length($msg) > 0)
            {
               #print "msg $msg\n";
               my($reciplist);

               if (!(defined($failedcheck{$projname})))
               {
                  print "here1, $projname\n";
                  $failedcheck{$projname} = {};
               }

               $reciplist = join(',', keys(%{$projowner{$projname}}));
               $failedcheck{$projname}->{'msg'} = $msg;

               print "failed, notifying $reciplist.\n";
            }
            else
            {
               print "passed.\n";
            }
         }
      }
   } # loop on environments
}

if (!$err)
{
   $err = Notify(\%projowner, \%failedcheck);
}

exit($err);


sub ReadCfg
{
   my($xml) = $_[0]; # reference
   my($line);
   my($rv);

   $rv = 0;
   $$xml = "";

   while (defined($line = <DATA>))
   {
      chomp($line);

      if ($line =~ /^\#/ || $line =~ /^\s*$/)
      {
         # Skip blank lines or lines beginning with # (comment lines)
         next;
      }
    
      # suck out xml
      $$xml = $$xml . "$line\n";
   } # loop over cfg

   return $rv;
}

sub Notify
{
   my($projowner) = $_[0];
   my($failedcheck) = $_[1];

   my($rv);
   my(@owners);
   my(@addrs);
   my(@projwerr);
   my($maillist);
   my($msg);
   my($user);
   my($host);
   my($cmd);

   $rv = 0;
   @projwerr = keys(%{$failedcheck});

   if ($#projwerr >= 0)
   {
      $user = getlogin() || getpwuid($<);
      $host = `hostname`;
      chomp($host);

      # Iterate through all projects, looking for failures
      foreach $proj (@projwerr)
      {
         @owners = keys(%{$projowner->{$proj}});
         @addrs = map {$projowner->{$proj}->{$_}->{'email'}} @owners;
         $maillist = join(' ', @addrs);
         $msg = $failedcheck->{$proj}->{'msg'};
         $cmd = "/bin/mail -s \"$user\@$host - production binary issue\" $maillist";

         print "Sending notification to $maillist.\n";
         print "Running $cmd\n";
exit;
         if (open(MAILPIPE, "| $cmd"))
         {
            print MAILPIPE $msg;
            close(MAILPIPE);
         }
         else
         {
            $rv = 1;
            print "Couldn't open 'mail' pipe.\n";
         }
      }
   }

   return $rv;
}

__DATA__
# This xml defines the set of 'production' modules that must be available in the production
# environment at all times. If any of these modules is missing, then we need to send out a 
# message to the appropriate owner(s) requesting that the modules be replaced.
<?xml version='1.0'?>
<enviros>
  <enviro>
    <name>dotrelease</name>
    <root>/home/jsoc/cvs/Development/JSOC</root>
    <projects>
      <proj>
        <name>observables</name>
        <owners>
          <owner login="arta" email="arta@sun.stanford.edu">Art</owner>
          <owner login="couvidat" email="couvidat@sun.stanford.edu">Sebastien</owner>
        </owners>
        <files>
          <file>
            <name>HMI_IQUV_averaging</name>
            <path>_linux_x86_64/proj/lev1.5_hmi/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>HMI_IQUV_averaging</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>HMI_observables</name>
            <path>_linux_x86_64/proj/lev1.5_hmi/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>HMI_observables</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
        </files>
      </proj>
      <proj>
        <name>datareplication</name>
        <owners>
          <owner login="arta" email="arta@sun.stanford.edu">Art</owner>
        </owners>
        <files>
          <file>
            <name>createns</name>
            <path>_suse_x86_64/base/drms/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>createns</name>
            <path>bin/suse_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>accessreplogs</name>
            <path>_suse_x86_64/base/drms/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>accessreplogs</name>
            <path>bin/suse_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>createtabstructure</name>
            <path>_suse_x86_64/base/drms/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>createtabstructure</name>
            <path>bin/suse_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>delete_series</name>
            <path>_suse_x86_64/base/util/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>delete_series</name>
            <path>bin/suse_x86_64</path>
            <perms>0777</perms>
          </file>
        </files>
      </proj>
      <proj>
        <name>magident</name>
        <owners>
          <owner login="arta" email="arta@sun.stanford.edu">Art</owner>
          <owner login="phil" email="phil@sun.stanford.edu">Phil</owner>
          <owner login="xudong" email="xudong@sun.stanford.edu">Xudong</owner>
          <owner login="yliu" email="yliu@sun.stanford.edu">Yang</owner>
        </owners>
        <files>
        <file>
            <name>hmi_segment_module</name>
            <path>_linux_x86_64/proj/mag/ident/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>hmi_segment_module</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>hmi_patch_module</name>
            <path>_linux_x86_64/proj/mag/ident/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>hmi_patch_module</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
        </files>
      </proj>
      <proj>
        <name>flatfield</name>
        <owners>
          <owner login="arta" email="arta@sun.stanford.edu">Art</owner>
          <owner login="phil" email="phil@sun.stanford.edu">Phil</owner>
          <owner login="richard" email="richard@sun.stanford.edu">Richard</owner>
          <owner login="jim" email="jim@sun.stanford.edu">Jim</owner>
          <owner login="rock" email="rock@sun.stanford.edu">Rock</owner>
        </owners>
        <files>
          <file>
            <name>module_flatfield</name>
            <path>_linux_x86_64/proj/flatfield/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>module_flatfield</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>module_flatfield_combine</name>
            <path>_linux_x86_64/proj/flatfield/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>module_flatfield_combine</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
          <file>
            <name>cosmic_ray_post</name>
            <path>_linux_x86_64/proj/flatfield/apps</path>
            <perms>0755</perms>
          </file>
          <file>
            <name>cosmic_ray_post</name>
            <path>bin/linux_x86_64</path>
            <perms>0777</perms>
          </file>
        </files>
      </proj>
    </projects>
  </enviro>
</enviros>

