#!/home/jsoc/bin/linux_x86_64/activeperl -w


# This script takes at least one argument which specifies which of three types of checkouts to perform:
#  1. A NetDRMS checkout (-net). If present, all other arguments are ignored.
#  2. A JSOC_SDP checkout (-sdp). If present, all other arugments are ignored.
#  3. A custom checkout (-custom). To perform a custom checkout, the caller must follow this argument
#     with the full path to a configuration file specifying the project directories the caller would
#     like to check-out.

# Each type of checkout contains a different subset of files that reside in the JSOC CVS module. To specify
# that subset, we maintain one "file specification" for each checkout type. A file specification is
# a list of file and directory paths relative to the CVS code-tree root. It turns out that you
# can use the CVS checkout or export command with these relative paths. So, for example, to check-out
# the sdp files, you can run "cvs checkout <relpath1> <relpath2> ... <relpathN>. This script maintains
# the file specification for each type of checkout, then performs the desired checkout or export.

# In order to determine the complete set of file names of the files that reside in a checkout file set, 
# you must download the files in that file set from the CVS repository. There is no CVS command that
# will print out all files in the repository or a subtree of the repository. 
# And it is not desirable to use the checkout command to download the files, because
# the checkout command creates additional CVS "state" files and places them in every node of the
# downloaded code tree. These intermingled state files make it difficult to isolate all the source files
# that comprise the checkout. To cope with these CVS deficiencies, when printing the list of files in a file set, 
# this script first EXPORTS all files in the JSOC module into a 
# temporary directory. When an export is performed, CVS does not introduce these extra state files.

# It is important to use this script to update the files in your working directory after the initial 
# checkout. Using the cvs update command directly can result in the download of files outside of the file set.
# And once that happens, 'make' may not work properly.

# NOTE: You cannot check-out the entire JSOC module, then delete the files not needed. If you do that, then CVS
# will think that the checkout is incomplete - cvsstatus.pl will indicate that the deleted files are 
# missing. The CVS state files, which list all files that were downloaded, including the the deleted files, 
# indicate to the CVS server that those deleted files are expected to be present in the working directory. If
# they are missing, various cvs commands will complain.

# flags:
#  -o The operation to perform, which includes:
#       checkout - The current directory must be the parent directory of the CVS working directory root.
#       export   - The current directory must be the parent directory of the CVS working directory root.
#       update   - update the set of files with changes committed to the CVS repository since the initial
#                  checkout. The current directory may be either the CVS working directory root, or its parent
#                  directory. The -f flag is ignored during an update - the checkout type is read from the
#                  configuration file. If the checkout type is sdp or net, then the file specification 
#                  for that checkout type is obtained from this script (not from the state file), since the 
#                  file specification for that checkout type may have changed since the original checkout.
#                  For a custom checkout, the file specification is re-generated. If the caller
#                  provides a -f <config file> argument, then the file specification for the proj directories
#                  is derived from <config file>. Otherwise, the proj directories file specification
#                  is obtained from the state file. Regardless, the file specification for the "core" 
#                  directories is obtained from this script.
#       tag      - tag the set of files in the CVS respository implied by the -f flag.
#       untag    - remove the tag on the  set of files in the CVS respository implied by the -f flag.
#       print    - print the set of files in the CVS repository implied by the -f flag.
#       printrel - print the set of files in the CVS repostiroy implied by the -f flag (but stripped 
#                  of files that need to be filtered out for releases).
#  -f The type of file set to operate on, which includes:
#       sdp (all files in the repository, aka the "full JSOC" source tree).
#       net (the set of files that compose the NetDRMS release).
#       <configuration file> (the set of files is specified by <configuration file>).
#  -r For the checkout, export, and update operations, this parameter is the CVS tag that identifies 
#       the revision of each file to download. For the sdp and net file-sets, this tag is applied to 
#       all files. For the custom file-set, this tag is applied to only the NetDRMS subset of files.
#  -R Applies to the custom file-set only. For the checkout, export, and update operations, this parameter is 
#       the CVS tag that identifies the revision of each non-NetDRMS project file to download. 
#  -t For the tag and untag operations, the CVS tag to apply or delete.
#  -l A log file (for the output of CVS commands for now).
#  -F If a tagged checkout occurs as a result of other options, then tell CVS to retrieve the most recent
#       file version of any file that is part of the file-set, but is not tagged with the tag provided
#       by the -r or -R option.
#  -s This argument can only be provided if the -o flag is update, tag, or untag. This argument contains 
#     a comma-separated list of file specifications. The files specified by this list must exist in the
#     working directory.
#  -d (optional) The root dir of the workspace CVS tree, which defaults to "JSOC".
#  -D (optional) The root dir of the repository CVS tree, which defaults to "JSOC".

use XML::Simple;
use IO::Dir;
use File::Copy;
use File::Basename;
use File::Path qw(mkpath remove_tree);
use File::Spec;
use Cwd qw(chdir getcwd realpath); # need to override chdir so that $ENV{'PWD'} is changed when chdir is called.
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../../../base/libs/perl";
use drmsLocks;

use constant kLockFile     => "/home/jsoc/locks/prodbuildlck.txt";

use constant kMakeDiv => "__MAKE__";
use constant kProjDiv => "__PROJ__";
use constant kFspecDev => "__FSPEC__";
use constant kEndDiv => "__END__";
use constant kStUnk => 0;
use constant kStMake => 1;
use constant kStProj => 2;

use constant kCoUnk => "unk";
use constant kCoNetDRMS => "net";
use constant kCoSdp => "sdp";
use constant kCoCustom => "custom";

use constant kDlCheckout => "checkout";
use constant kDlExport => "export";
use constant kDlUpdate => "update";
use constant kDlTag => "tag";
use constant kDlUntag => "untag";
use constant kDlPrint => "print";
use constant kDlPrintRelease => "printrel";

use constant kStrproj => "proj";
use constant kStrname => "name";

# Assume a localization directory of "localization" right in the root of the CVS tree.
use constant kDefRootDir => "JSOC/"; # default root dir, overridable by the -d argument
use constant kDefRepRootDir => "JSOC/"; # default repository root dir, overridable by the -
use constant kLocDir => "localization/";
use constant kProjSubdir => "proj";
use constant kTmpDir => "/tmp/chkout/";
use constant kTypeFile => "dlset.txt";
use constant kSuFlagFile => "suflag.txt";

use constant kMaxFileSpecs => 50;


my($arg);
my($cotype);
my($cfgfile);
my($logfile);
my($cmd);
my($err);
my(@core);
my(@netonly);
my(@sdponly);
my(@netfilter);
my(@sdpfilter);
my(@netco);
my(@sdbco);
my($curdir);
my($xmldata); # reference to hash array
my($dltype);
my($version);
my($pversion); # version of project files (when for net and custom file-set types)
my($cvstag);
my($stfile);
my($stfileold);
my($stcotype);
my($stfspec);
my($compatmode);
my($forceco);
my(@filespec); # the complete file spec for all types of file-sets
my(@pfilespec); # for custom file-set types, the file spec of the files in the configuration file
my(@bfilespec); # for custom file-set types, the file spec of the files NOT in the configuration file
my(@cmdlspec); # can specify file specifications on the cmd-line.
my($actcvs); # if the user provided a cmdlspec, then this var holds the cvs command used to update
             # the files specified by $cmdlspec.
my(@actspec); # The actual file specification that was used to download files from the repository.
my($rdir);     # root dir (either specified by -d flag, or default root dir of JSOC)
my($reprdir); # root dir (either specified by -D flag, or default root dir of JSOC)


# Don't allow more than one version of this file to run concurrently to avoid race conditions.
$lock = new drmsNetLocks(&kLockFile);

if (!defined($lock))
{
   print "$0 is already running. Exiting.\n";
   exit(1);
}

@core = qw(base/cfortran.h base/foundation.h base/jsoc.h base/jsoc_version.h base/mypng.h base/Rules.mk base/export base/drms base/libs base/sums base/util localize.py configure doc make_basic.mk Makefile make_jsoc.pl README Rules.mk target.mk build CM);
@coreDel = qw(configproj.pl customizemake.pl moreconfigure.pl getmachtype.pl);

@netonly = qw(config.local.template config.local.map seed_sums.c netdrms_setup.pl proj/example proj/myproj proj/cookbook);
@netDel = qw(gen_init.csh getuid.c);

@sdponly = qw(base/local proj configsdp.txt);
@sdpDel = qw(customizedefs.pl config.local.sutemplate);

# My botched attempt to remove tape files from NetDRMS sites that will never use our tape code
#@netfilter = qw(base/drms/doc base/drms/libs/api/test base/sums/libs/api/perl base/sums/libs/api/tape.h base/sums/libs/pg/SUMLIB_DS_DataRequest_Tst.pgc base/sums/libs/pg/SUMLIB_NC_PaRequest_AP_60d.pgc base/sums/libs/pg/SUMLIB_TapeClose.pgc base/sums/libs/pg/SUMLIB_TapeFindGroup.pgc base/sums/libs/pg/SUMLIB_TapeUpdate.pgc base/sums/apps/main.c base/sums/apps/main2.c base/sums/apps/main3.c base/sums/apps/main4.c base/sums/apps/main5.c base/sums/apps/robotn_svc.c base/sums/apps/sum_forker.c base/sums/apps/sum_test.c base/sums/apps/sum_test.pl base/sums/apps/tapearc.c base/sums/apps/tapearc0.c base/sums/apps/tapearc1.c base/sums/apps/tapearc2.c base/sums/apps/tapearc3.c base/sums/apps/tapearc4.c base/sums/apps/tapearc5.c base/sums/apps/tapearc6.c base/sums/apps/tapearc7.c base/sums/apps/tapearc8.c base/sums/apps/tapearcinfo.c base/sums/apps/tapearcX.c base/sums/apps/tape_inventory.c base/sums/apps/tapeonoff.c base/sums/apps/tape_svc.c base/sums/apps/tape_svc_proc.c base/sums/apps/tapeutil.c base/sums/apps/xsum_svc.c base/sums/apps/xsum_svc_proc.c base/sums/apps/xtape_svc.c base/sums/scripts/build_parc_file.pl base/sums/scripts/find_dir_sum_partn_alloc_dc base/sums/scripts/fixportm.pl base/sums/scripts/get_dcs_times.csh base/sums/scripts/GRAD_BLUE_LINE.gif base/sums/scripts/lev1_def_gui base/sums/scripts/lev1_def_gui_aia base/sums/scripts/lev1_def_gui_called base/sums/scripts/lev1_def_gui_called_PZT_FSN base/sums/scripts/lev1_def_gui_hmi base/sums/scripts/rsync_scr111.pl base/sums/scripts/SDO_Badge.gif base/sums/scripts/SDO_HSB_CCSDS_Data_Structures.gif base/sums/scripts/ssh_rsync.source base/sums/scripts/sum_bad_permissions.pl base/sums/scripts/sumck base/sums/scripts/sumck_j1 base/sums/scripts/sumck_j1M base/sums/scripts/sumck_n02_jim base/sums/scripts/sumlookgroup.pl base/sums/scripts/sumlook.pl base/sums/scripts/sum_start base/sums/scripts/sum_start_d00_jim base/sums/scripts/sum_start_d02 base/sums/scripts/sum_start_d02_auto base/sums/scripts/sum_start_dc base/sums/scripts/sum_start_j1 base/sums/scripts/sum_start_j1_auto base/sums/scripts/sum_start_j1_auto.MULTI base/sums/scripts/sum_start_j1.MULTI base/sums/scripts/sum_start_n02_jim base/sums/scripts/sum_start_n02_jim_auto base/sums/scripts/sum_start_xim.MULTI base/sums/scripts/sum_stop base/sums/scripts/sum_stop_d00_jim base/sums/scripts/sum_stop_d02 base/sums/scripts/sum_stop_d02_auto base/sums/scripts/sum_stop_d02_tape base/sums/scripts/sum_stop_dc base/sums/scripts/sum_stop_j1 base/sums/scripts/sum_stop_j1_auto base/sums/scripts/sum_stop_j1_auto.MULTI base/sums/scripts/sum_stop_j1.MULTI base/sums/scripts/sum_stop_n02_jim base/sums/scripts/sum_stop_n02_jim_auto base/sums/scripts/sum_stop_xim.MULTI base/sums/scripts/sum_tape_catchup_update.pl base/sums/scripts/sum_tape_insert.pl base/sums/scripts/sum_tape_insert_t50.pl base/sums/scripts/sum_tape_insert_t950.pl base/sums/scripts/t120_reachive.pl base/sums/scripts/t120stageall.pl base/sums/scripts/t120view base/sums/scripts/t50view base/sums/scripts/t950view base/sums/scripts/tapearc_do base/sums/scripts/tapearc_do_dcs1 base/sums/scripts/tape_do_0.pl base/sums/scripts/tape_do_1.pl base/sums/scripts/tape_do_2.pl base/sums/scripts/tape_do_3.pl base/sums/scripts/tape_do_4.pl base/sums/scripts/tape_do_7.pl base/sums/scripts/tape_do_8.pl base/sums/scripts/tape_do_archive.pl base/sums/scripts/tape_do.pl base/sums/scripts/tapeid.list base/sums/scripts/tapeid_t50.list base/sums/scripts/tape_verify.pl base/sums/scripts/test base/sums/scripts/tmp.pl doc/dcs2_convert_to_0_or_1.txt doc/dcs3_name_change.txt doc/dcs_warmstandby.txt doc/dsc0_just_rebooted.txt doc/HK_Level0_Debug_Guide.odt doc/HK_Level0_Debug_Guide.pdf doc/whattodo_aia_lev1.txt doc/whattodo_dcs.txt doc/whattodolev0.txt doc/whattodo_start_stop_lev1_0_sums.txt);

@netfilter = qw(base/drms/doc base/drms/libs/api/test base/sums/libs/api/perl base/sums/libs/pg/SUMLIB_DS_DataRequest_Tst.pgc base/sums/libs/pg/SUMLIB_NC_PaRequest_AP_60d.pgc base/sums/apps/main.c base/sums/apps/main2.c base/sums/apps/main3.c base/sums/apps/main4.c base/sums/apps/main5.c base/sums/apps/sum_test.c base/sums/apps/sum_test.pl base/sums/apps/xsum_svc.c base/sums/apps/xsum_svc_proc.c base/sums/apps/xtape_svc.c base/sums/scripts/build_parc_file.pl base/sums/scripts/find_dir_sum_partn_alloc_dc base/sums/scripts/fixportm.pl base/sums/scripts/get_dcs_times.csh base/sums/scripts/GRAD_BLUE_LINE.gif base/sums/scripts/lev1_def_gui base/sums/scripts/lev1_def_gui_aia base/sums/scripts/lev1_def_gui_called base/sums/scripts/lev1_def_gui_called_PZT_FSN base/sums/scripts/lev1_def_gui_hmi base/sums/scripts/rsync_scr111.pl base/sums/scripts/SDO_Badge.gif base/sums/scripts/SDO_HSB_CCSDS_Data_Structures.gif base/sums/scripts/ssh_rsync.source base/sums/scripts/sum_bad_permissions.pl base/sums/scripts/sumck base/sums/scripts/sumck_j1 base/sums/scripts/sumck_j1M base/sums/scripts/sumck_n02_jim base/sums/scripts/sumlookgroup.pl base/sums/scripts/sumlook.pl base/sums/scripts/sum_start base/sums/scripts/sum_start_d00_jim base/sums/scripts/sum_start_d02 base/sums/scripts/sum_start_d02_auto base/sums/scripts/sum_start_dc base/sums/scripts/sum_start_j1 base/sums/scripts/sum_start_j1_auto base/sums/scripts/sum_start_j1_auto.MULTI base/sums/scripts/sum_start_j1.MULTI base/sums/scripts/sum_start_n02_jim base/sums/scripts/sum_start_n02_jim_auto base/sums/scripts/sum_start_xim.MULTI base/sums/scripts/sum_stop base/sums/scripts/sum_stop_d00_jim base/sums/scripts/sum_stop_d02 base/sums/scripts/sum_stop_d02_auto base/sums/scripts/sum_stop_d02_tape base/sums/scripts/sum_stop_dc base/sums/scripts/sum_stop_j1 base/sums/scripts/sum_stop_j1_auto base/sums/scripts/sum_stop_j1_auto.MULTI base/sums/scripts/sum_stop_j1.MULTI base/sums/scripts/sum_stop_n02_jim base/sums/scripts/sum_stop_n02_jim_auto base/sums/scripts/sum_stop_xim.MULTI base/sums/scripts/sum_tape_catchup_update.pl base/sums/scripts/sum_tape_insert.pl base/sums/scripts/sum_tape_insert_t50.pl base/sums/scripts/sum_tape_insert_t950.pl base/sums/scripts/t120_reachive.pl base/sums/scripts/t120stageall.pl base/sums/scripts/t120view base/sums/scripts/t50view base/sums/scripts/t950view base/sums/scripts/tapearc_do base/sums/scripts/tapearc_do_dcs1 base/sums/scripts/tape_do_0.pl base/sums/scripts/tape_do_1.pl base/sums/scripts/tape_do_2.pl base/sums/scripts/tape_do_3.pl base/sums/scripts/tape_do_4.pl base/sums/scripts/tape_do_7.pl base/sums/scripts/tape_do_8.pl base/sums/scripts/tape_do_archive.pl base/sums/scripts/tape_do.pl base/sums/scripts/tapeid.list base/sums/scripts/tapeid_t50.list base/sums/scripts/tape_verify.pl base/sums/scripts/test base/sums/scripts/tmp.pl doc/dcs2_convert_to_0_or_1.txt doc/dcs3_name_change.txt doc/dcs_warmstandby.txt doc/dsc0_just_rebooted.txt doc/HK_Level0_Debug_Guide.odt doc/HK_Level0_Debug_Guide.pdf doc/whattodo_aia_lev1.txt doc/whattodo_dcs.txt doc/whattodolev0.txt doc/whattodo_start_stop_lev1_0_sums.txt);

@sdpfilter = qw(base/drms/doc base/drms/libs/api/test base/sums/libs/api/perl base/sums/libs/pg/SUMLIB_DS_DataRequest_Tst.pgc base/sums/libs/pg/SUMLIB_NC_PaRequest_AP_60d.pgc base/sums/apps/main.c base/sums/apps/main2.c base/sums/apps/main3.c base/sums/apps/main4.c base/sums/apps/main5.c base/sums/apps/sum_test.c base/sums/apps/sum_test.pl base/sums/apps/tapearcX.c base/sums/apps/xsum_svc.c base/sums/apps/xsum_svc_proc.c base/sums/apps/xtape_svc.c base/sums/scripts/fixportm.pl base/sums/scripts/get_dcs_times.csh base/sums/scripts/GRAD_BLUE_LINE.gif base/sums/scripts/lev1_def_gui base/sums/scripts/lev1_def_gui_aia base/sums/scripts/lev1_def_gui_called base/sums/scripts/lev1_def_gui_called_PZT_FSN base/sums/scripts/lev1_def_gui_hmi base/sums/scripts/rsync_scr111.pl base/sums/scripts/SDO_Badge.gif base/sums/scripts/SDO_HSB_CCSDS_Data_Structures.gif base/sums/scripts/ssh_rsync.source base/sums/scripts/sum_bad_permissions.pl base/sums/scripts/sumck base/sums/scripts/sumck_j1 base/sums/scripts/sumck_j1M base/sums/scripts/sumck_n02_jim base/sums/scripts/sumlookgroup.pl base/sums/scripts/sumlook.pl base/sums/scripts/sum_start base/sums/scripts/sum_start_d00_jim base/sums/scripts/sum_start_d02 base/sums/scripts/sum_start_d02_auto base/sums/scripts/sum_start_dc base/sums/scripts/sum_start_j1 base/sums/scripts/sum_start_j1_auto base/sums/scripts/sum_start_j1_auto.MULTI base/sums/scripts/sum_start_j1.MULTI base/sums/scripts/sum_start_n02_jim base/sums/scripts/sum_start_n02_jim_auto base/sums/scripts/sum_start_xim.MULTI base/sums/scripts/sum_stop base/sums/scripts/sum_stop_d00_jim base/sums/scripts/sum_stop_d02 base/sums/scripts/sum_stop_d02_auto base/sums/scripts/sum_stop_d02_tape base/sums/scripts/sum_stop_dc base/sums/scripts/sum_stop_j1 base/sums/scripts/sum_stop_j1_auto base/sums/scripts/sum_stop_j1_auto.MULTI base/sums/scripts/sum_stop_j1.MULTI base/sums/scripts/sum_stop_n02_jim base/sums/scripts/sum_stop_n02_jim_auto base/sums/scripts/sum_stop_xim.MULTI base/sums/scripts/sum_tape_catchup_update.pl base/sums/scripts/sum_tape_insert.pl base/sums/scripts/sum_tape_insert_t50.pl base/sums/scripts/sum_tape_insert_t950.pl base/sums/scripts/t120_reachive.pl base/sums/scripts/t120stageall.pl base/sums/scripts/t120view base/sums/scripts/t50view base/sums/scripts/t950view base/sums/scripts/tapearc_do base/sums/scripts/tapearc_do_dcs1 base/sums/scripts/tape_do_0.pl base/sums/scripts/tape_do_1.pl base/sums/scripts/tape_do_2.pl base/sums/scripts/tape_do_3.pl base/sums/scripts/tape_do_4.pl base/sums/scripts/tape_do_7.pl base/sums/scripts/tape_do_8.pl base/sums/scripts/tape_do_archive.pl base/sums/scripts/tape_do.pl base/sums/scripts/tapeid.list base/sums/scripts/tapeid_t50.list base/sums/scripts/tape_verify.pl base/sums/scripts/test base/sums/scripts/tmp.pl doc/dcs2_convert_to_0_or_1.txt doc/dcs3_name_change.txt doc/dcs_warmstandby.txt doc/dsc0_just_rebooted.txt doc/HK_Level0_Debug_Guide.odt doc/HK_Level0_Debug_Guide.pdf doc/whattodo_aia_lev1.txt doc/whattodo_dcs.txt doc/whattodolev0.txt doc/whattodo_start_stop_lev1_0_sums.txt);

$err = 0;
$cotype = kCoSdp;
$dltype = kDlCheckout;
$version = "";
$cvstag = "";
$compatmode = 0;
$forceco = 0;
$rdir = &kDefRootDir;
$reprdir = &kDefRepRootDir;

while ($arg = shift(@ARGV))
{
    if ($arg eq "-o")
    {
        # download type
        $arg = shift(@ARGV);
        if ($arg eq kDlCheckout ||
            $arg eq kDlExport ||
            $arg eq kDlUpdate ||
            $arg eq kDlTag ||
            $arg eq kDlUntag ||
            $arg eq kDlPrint ||
            $arg eq kDlPrintRelease)
        {
            $dltype = $arg;
        }
        else
        {
            print STDERR "Invalid download type - please choose from 'checkout', 'export', 'update', 'tag', 'untag', 'print', or 'printrel'.\n";
            $err = 1;
            last;
        }
    }
    elsif ($arg eq "-r")
    {
        # revision (version)
        $arg = shift(@ARGV);
        $version = $arg;
    }
    elsif ($arg eq "-R")
    {
        $arg = shift(@ARGV);
        $pversion = $arg;
    }
    elsif ($arg eq "-t")
    {
        # CVS tag to set/remove
        $arg = shift(@ARGV);
        $cvstag = $arg;
    }
    elsif ($arg eq "-l")
    {
        $arg = shift(@ARGV);
        $logfile = $arg;
    }
    elsif ($arg eq "-f")
    {
        # file set
        $arg = shift(@ARGV);
        
        if ($arg eq kCoSdp)
        {
            $cotype = kCoSdp;
        }
        elsif ($arg eq kCoNetDRMS)
        {
            $cotype = kCoNetDRMS;
        }
        else
        {
            # custom - argument must be a configuration file
            if (-f $arg)
            {
                $cotype = kCoCustom;
                $cfgfile = $arg;
            }
            else
            {
                print STDERR "Invalid custom-download configuration file $arg.\n";
                $err = 1;
            }
        }
    }
    elsif ($arg eq "-F")
    {
        $forceco = 1;
    }
    elsif ($arg eq "-s")
    {
        $arg = shift(@ARGV);
        @cmdlspec = split(/,/, $arg);
    }
    elsif ($arg eq "-d")
    {
        $arg = shift(@ARGV);
        $rdir = $arg;
    }
    elsif ($arg eq "-D")
    {
        $arg = shift(@ARGV);
        $reprdir = $arg;
    }
}

if (!$err)
{
    my($inparent); # if 1, curr dir is the parent of JSOC root dir.
    my($crootdir);  # canonical root dir
    my($creprootdir); # canonical repository root dir (i.e., JSOC)
    
    $inparent = 0;
   if (defined($cfgfile))
   {
      # Custom checkout - if the user is performing an update, then use the file spec saved
      # in the TYPEFILE (no need to re-read a config file).
      my($xml);
      my($xmlobj) = new XML::Simple;

      # Read in the configuration file to obtain the set of project files that will reside
      # in the custom checkout set.
      if (!ReadCfg($cfgfile, \$xml) && defined($xml))
      {
         $xmldata = $xmlobj->XMLin($xml, ForceArray => 1);
      }
      else
      {
         print STDERR "Unable to read or parse configuration file $cfgfile.\n";
         $err = 1;
      }
   }

   if ($dltype eq kDlCheckout || $dltype eq kDlExport || $dltype eq kDlUpdate)
   {
       # Set the state file path.
       $crootdir = File::Spec->catdir($rdir);
       $creprootdir = File::Spec->catdir($reprdir);
       #my($cdir) = File::Spec->catdir($ENV{'PWD'});
       my($cdir) = realpath($ENV{'PWD'});
       
       $stfile = kLocDir . kTypeFile;
       $stfileold = kSuFlagFile;
       
       if ($cdir !~ /$crootdir\s*$/)
       {
           # Assume that the current directory is the parent of the JSOC code tree.
           $stfile = $rdir . $stfile;
           $stfileold = $rdir . $stfileold;
           $inparent = 1;
       }
   }

   if ($dltype eq kDlUpdate)
   {
      # If this is an update, obtain the checkout type from the statefile. The 
      # state file will exist at this point only for kDlUpdate.
      if (open(STFILE, "<$stfile"))
      {
         # Modify $cotype - should be determined by the first line of the state file.
         my($line);

         $line = <STFILE>;
         chomp($line);
         $stcotype = $line;
         $cotype = $stcotype;
         $line = <STFILE>;
         chomp($line);
         $stfspec = $line;

         close(STFILE);
      }
      elsif (!(-e $stfile))
      {
         # Backward compatibility for previous versions of the cvs tree.
         if (-e $stfileold)
         {
            # Assume old sdp tree.
            $cotype = kCoSdp;
         }
         else
         {
            # Assume old net tree.
            $cotype = kCoNetDRMS;
         }

         $compatmode = 1;
      }
      else
      {
         print STDERR "Unable to open state file '$stfile' for reading.\n";
         $err = 1;
      }
   }

   if ($err)
   {
      # Do nothing - this will essentially cause this script to exit.
   }
   elsif ($cotype ne kCoSdp && $cotype ne kCoNetDRMS && $cotype ne kCoCustom)
   {
      print STDERR "Invalid file set identifier '$cotype'.\n";
      $err = 1;
   }
   elsif ($dltype ne kDlCheckout && $dltype ne kDlExport && $dltype ne kDlUpdate && 
          $dltype ne kDlTag && $dltype ne kDlUntag && $dltype ne kDlPrint && $dltype ne kDlPrintRelease)
   {
      print STDERR "Invalid operation '$dltype'.\n";
      $err = 1;
   }
   else
   { 
       my(@cmdlspecrel);
       
       if ($#cmdlspec >= 0)
       {
           if ($inparent)
           {
               # Prepend each spec with the $rootdir
               @cmdlspecrel = map({($creprootdir . $_)} @cmdlspec);
           }
           else
           {
               @cmdlspecrel = @cmdlspec;
           }
       }
       
      if (BuildFilespec($cotype, $dltype, $stfspec, $xmldata, \@core, \@netonly, \@sdponly, \@filespec, \@bfilespec, \@pfilespec, \@cmdlspecrel, \@coreDel, \@netDel, \@sdpDel))
      {
         print STDERR "Unable to build filespec.\n";
         $err = 1;
      }
      else
      {
         undef($curdir);

         if ($dltype eq kDlTag || $dltype eq kDlUntag || $dltype eq kDlPrint || $dltype eq kDlPrintRelease)
         {
            # cd to tmp directory for these commands
            if (!(-d kTmpDir))
            {
               # no need to check this call, because the chdir() cmd is being checked.
               mkpath(kTmpDir);
            }
           
            $curdir = $ENV{'PWD'};
            $err = (chdir(kTmpDir) == 0);
            if ($err)
            {
               print STDERR "Unable to cd to " . kTmpDir . ".\n";
            }
         }
         elsif ($dltype eq kDlUpdate)
          {
              # If the current directory is the root directory of the CVS working directory, then
              # cd up to the parent directory (DownloadTree assumes the current directory is the
              # parent of the CVS working directory).
              my($crootdir) = File::Spec->catdir($rdir);
              my($creprootdir) = File::Spec->catdir($reprdir);
              #my($cdir) = File::Spec->catdir($ENV{'PWD'});
              my($cdir) = realpath($ENV{'PWD'});
              
              print "dlsource.pl: rootdir == $crootdir, cdir == $cdir.\n";
              if ($cdir =~ /$crootdir\s*$/)
              {
                  # The current directory is the CVS working directory.
                  $curdir = $ENV{'PWD'};
                  $err = (chdir('..') == 0) ? 1 : 0;
              }
          }

         # Do a cvs checkout, export, or update into the current directory
         if (!$err)
         {
             # @filespec - the complete set of files that reside on the server for the current check-out type.
             # @bfilespec - for non-custom check-out, the set of files to update or checkout. For 
             #              custom check-out, the set of files in the base dir to update or checkout.
             # @pfilespec - for custom check-out, the set of files in the proj dir to update or checkout.
             # @actspec - the spec actually used with the cvs command. Might differ from the others if
             #            an export was done and undesirable files were filtered out.
            $err = DownloadTree($cotype, $dltype, $version, $pversion, \@filespec, \@bfilespec, \@pfilespec, \@cmdlspec, $logfile, $forceco, \@netfilter, \@sdpfilter, \$actcvs, \@actspec);

            if ($err)
            {
               print STDERR "Unable to $dltype CVS tree.\n";
            }
            elsif ($compatmode)
            {
               # Must remove kSuFlagFile if it was present.
               # compatmode is used only during an update, so the current
               # directory is the parent of the code root directory.
               if (-e $rdir . kSuFlagFile)
               {
                  my($cvscmd) = "cvs update -A " . $rdir . kSuFlagFile . " " . $rdir . "jsoc_sync.pl " . $rdir . "jsoc_update.pl";

                  if (CallCVS($cvscmd, $logfile, undef, 0))
                  {
                     print STDERR "Unable to run $cvscmd.\n";
                     $err = 1;
                  }
                  elsif (-e $rdir . kSuFlagFile || -e $rdir . "jsoc_sync.pl" || -e $rdir . "jsoc_update.pl")
                  {
                     print STDERR "Unable to delete old suflag.txt, jsoc_sync.pl, or jsoc_update.pl\n";
                     $err = 1;
                  }
               }
            }
         }

         if (!$err)
         {
            # Assumes that the cdir is the one containing kRootDir
            if ($dltype eq kDlCheckout || $dltype eq kDlExport || $dltype eq kDlUpdate)
            {
                my($typefile);
                my($tmptypefile);
                my($ierr) = 0;
                
                $typefile = $rdir . &kLocDir . &kTypeFile;
                $tmptypefile = $rdir . &kLocDir . "." . &kTypeFile . ".tmp";
                if (!(-d $rdir . kLocDir))
                {
                    mkpath($rdir . kLocDir);
                }
                
                # Copy original version of TYPEFILE
                if (-e $typefile)
                {
                    copy($typefile, $tmptypefile);
                }

                # save state file
                if (open(TYPEFILE, ">" . $typefile))
                {
                    # save check-out type                 
                    print TYPEFILE "$cotype\n";
                    
                    # print file spec used during checkout
                    my($fs) = join(' ', @filespec);
                    print TYPEFILE "$fs\n";                  
                    
                    # Now print list of files that compose the file set.
                    if ($dltype eq kDlUpdate)
                    {
                        # Must print each tree in file specification, since
                        # the tree rooted at kRootDir may contain
                        # files other than the files originally downloaded
                        # from CVS (e.g., running make will create new 
                        # files).
                        
                        # The update may have downloaded additional files that were not
                        # in the original check-out set. And, only a subset of files 
                        # may have been updated. So, read in the previous set of 
                        # files in the TYPEFILE, then add to this list files new files
                        # that were downloaded by the update.
                        my(@combined);
                        my(@sorted);
                        my(@listf);
                        my($sdir);
                        my(@oldflist);
                        my(@dlist);
                        my(%seen);
                        my($lastindx);
                        
                        # Read-in the list of files in the old TYPEFILE.
                        if (open(OLDTYPEFILE, "<" . $tmptypefile))
                        {
                            @oldflist = <OLDTYPEFILE>;
                            $lastindx = $#oldflist;          
                            push(@combined, @oldflist[2..$lastindx]);
                            close(OLDTYPEFILE);

                            # Export-download the update filespec into the tmp directory. The 
                            # user may have provided a cmd-line spec, in which case it should
                            # be used for the download, and not the full file spec. The previous
                            # DownloadTree() did perform this logic, and saved the cvs command
                            # used in $actcvs.
                            if (!(-d kTmpDir))
                            {
                                # no need to check this call, because the chdir() cmd is being checked.
                                mkpath(kTmpDir);
                            }

                            $sdir = $ENV{'PWD'};
                            if (chdir(kTmpDir) == 0)
                            {
                                print STDERR "Unable to cd to " . kTmpDir . ".\n";
                                $ierr = 1;
                            }
                            else
                            {
                                # The $actcvs command won't be quite right. We need to change the
                                # word 'update -APd' to export, and if there is no '-r' flag, we need to 
                                # add '-r HEAD' after update.
                                $actcvs =~ s/update\s+-APd\s+/export /;
                                if ($actcvs !~ /-r\s/)
                                {                                    
                                    $actcvs =~ s/export\s/export -r HEAD /;
                                }

                                unless (CallCVSInChunks($actcvs, \@actspec, undef, undef, 1))
                                {
                                    unless (GetFileList(&kTmpDir . $rdir, "", \@dlist))
                                    {
                                        # Combine the list of files in the downloaded tree with the list of
                                        # files in @combined.
                                        my($tmpdiri) = &kTmpDir;
                                        foreach my $afile (@dlist)
                                        {
                                            # remove &kTmpDir
                                            $afile =~ s/$tmpdiri//;
                                            push(@combined, "$afile\n");
                                        }
                                        
                                        # Sort and extract list of unique file names.
                                        @sorted = sort(@combined);
                                        @listf = map({ unless ($seen{$_}++){($_)}else{()} } @sorted);
                                        
                                        # Finally, print list.
                                        foreach my $afile (@listf)
                                        {
                                            print TYPEFILE $afile;
                                        }
                                    }
                                    else
                                    {
                                        $ierr = 1;
                                    }
                                                   
                                    remove_tree(&kTmpDir);
                                }
                                else
                                {
                                    $ierr = 1;
                                }
                                
                                chdir($sdir);
                            }
                        }
                        else
                        {
                            print STDERR "Unable to open $tmptypefile for reading.\n";
                            $ierr = 1;
                        }
                    }
                    else
                    {                        
                        if (PrintFilenames(*TYPEFILE, $rdir, 1, qr(\/CVS\/)))
                        {
                            print STDERR "Unable to print file-set file names.\n";
                            $ierr = 1;
                        }
                    }
                    
                    close(TYPEFILE);
                }
                else
                {
                    print STDERR "Unable to open file " . $rdir . kLocDir . kTypeFile . " for writing.\n";
                    $ierr = 1;
                }
                
                if (!$ierr && -e $tmptypefile)
                {
                    unlink($tmptypefile);
                }
                
                # Copy the cvs update log, if it exists, back to the kRootDir (programs calling this script)
                # expect it in kRootDir.
                if (defined($logfile) && -e $logfile)
                {
                    # cp $logfile kRootDir
                    if (!copy($logfile, $rdir))
                    {
                        # copy failure
                        print STDERR "Unable to copy log file $logfile to " . $rdir . ".\n";
                        $err = 1;
                    }
                }
            }
            elsif ($dltype eq kDlTag || $dltype eq kDlUntag)
            {
               # Can only use the tag/untag command to tag/untag releases, 
               # which can only be done on either the sdp or netdrms checkouts types.
               if ($cotype eq kCoSdp || $cotype eq kCoNetDRMS)
               {
                  my($curdirin) = $ENV{'PWD'};
                  $err = (chdir(kTmpDir . $rdir) == 0);
                  if ($err)
                  {
                     print STDERR "Unable to cd to " . kTmpDir . $rdir . ".\n";
                  }
                  else
                  {
                     if (TagFiles($cvstag, $dltype, $logfile))
                     {
                        print STDERR "Unable to tag/untag files in file specification.\n";
                        $err = 1;
                     }
                  }

                  if (chdir($curdirin) == 0)
                  {
                     print STDERR "Unable to cd to $curdir.\n";
                     if (!$err)
                     {
                        $err = 1;
                     }
                  }
               }
               else
               {
                  print STDERR "Checkout type $cotype not compatible with tag command.\n";
                  $err = 1;
               }
            }
            elsif ($dltype eq kDlPrint || $dltype eq kDlPrintRelease)
            {
               if (PrintFilenames(*STDOUT, kTmpDir . $rdir, 0))
               {
                  print STDERR "Unable to print file set file names.\n";
                  $err = 1;
               }
            }
         }

         if (defined($curdir))
         {
            # This implies that a successful chdir was done previously.
            if (chdir($curdir) == 0)
            {
               print STDERR "Unable to cd to $curdir.\n";
               $err = 1;
            }
         }

         # Delete all files from temporary directory, but only if there were no errors
         if (!$err)
         {
            if (-d kTmpDir)
            {
               remove_tree(kTmpDir, {error => \my $errlist});
               if (@{$errlist})
               {
                  print STDERR "Unable to properly remove temporary subdirectory $tmpdir.\n";
                  $err = 1;
               }
            }
         }
      }
   }
}

exit($err);

sub ReadCfg
{
   my($cfgfile) = $_[0];
   my($xml) = $_[1]; # reference
   my($rv);

   $rv = 0;

   # initialize xml string variable
   $$xml = "";
   
   if (defined($cfgfile) && -e $cfgfile)
   {
      if (open(CFGFILE, "<$cfgfile"))
      {
         my($projdone);
         my($line);

         $st = kStUnk;
         $pdiv = kProjDiv;
         $ediv = kEndDiv;

         $projdone = 0;

         while (defined($line = <CFGFILE>))
         {
            chomp($line);

            if ($line =~ /^\#/ || $line =~ /^\s*$/)
            {
               # Skip blank lines or lines beginning with # (comment lines)
               next;
            }
            elsif ($line =~ /^$pdiv/)
            {
               $st = kStProj;

              
               next;
            }
            elsif ($line =~ /^$ediv/)
            {
               if ($st == kStProj)
               {
                  $projdone = 1;
               }

               $st = kStUnk;

               if ($projdone)
               {
                  last;
               }

               next;
            }

            if ($st == kStProj)
            {
               # suck out xml
               $$xml = $$xml . "$line\n";
            }
         } # loop over cfg file

         close(CFGFILE);
      }
      else
      {
         print STDERR "Unable to open $cfgfile for reading.\n";
         $rv = 1;
      }
   }
   else
   {
      print STDERR "Unable to open $cfgfile for reading.\n";
      $rv = 1;
   }

   return $rv;
}

# fsout - The full file spec for the given type of checkout.
# bfsout - Same as fsout. Used only for custom checkouts.
# pfsout - For custom checkouts only. This is the project-directory file spec.
sub BuildFilespec
{
   my($cotype) = $_[0];
   my($dltype) = $_[1];
   my($stfspec) = $_[2];
   my($xmldata) = $_[3]; # reference to hash
   my($core) = $_[4];
   my($netonly) = $_[5];
   my($sdponly) = $_[6];
   my($fsout) = $_[7];
   my($bfsout) = $_[8];
   my($pfsout) = $_[9];
   my($cmdlspec) = $_[10]; # Altenate file spec - use in place of $core, $netonly, and $sdponly.
                            # Can be used only for an update, tag, or untag download type.
   my($coreDel) = $_[11];
   my($netDel) = $_[12];
   my($sdpDel) = $_[13];
    
   my($rv);
   my($strproj) = &kStrproj;
   my($strname) = &kStrname;

   $rv = 0;
    
    if ($#cmdlspec >= 0)
    {
        if ($cotype ne kCoCustom && ($dltype eq kDlUpdate || $dltype eq kDlTag || $dltype eq kDlTag))
        {
            my($compflistH);
            my(@fset);
            
            # Ensure that files specified on the cmd-line are actually part of the current 
            # check-out type (net, sdp, or custom). To obtain the list of the current check-out set, 
            # we need to check out everything in the appropriate definitive file spec (as an export).
            
            # This function returns a hash whose element keys are file names. If the filename-key
            # represents a regualar file, then the value for the key is 'f'. If the filename-key
            # represents a directory, then the value for this key is a hash with one element
            # per file in this directory. 
            # 
            # If there is a regular file in $cmdlspec, then traverse $compflistH until the 
            # parent directory is found. Then verify that the regular file is in the hash
            # of this parent directory. If it is found, then add the regular file
            # to the array of file-specs returned to the caller. If there is a directory in $cmdlspec, 
            # traverse $compflistH until the directory is found. Then collect ALL regular files rooted
            # at this directory and add those to the array of file-specs returned to the caller.
            if (!$rv)
            {
                $compflistH = GetDefinitiveFileSet($cotype, $core, $netonly, $sdponly, $xmldata, $stfspec, $version, $pversion, $reprdir, $rdir);
                
                foreach my $spec (@$cmdlspec)
                {
                    @fset = CollectFiles($spec, $compflistH);
                    if ($#fset >= 0)
                    {
                        # ART - If this is a custom check-out, then this is wrong. We'd need to 
                        # sort the files into those in the proj dir, and those in the base dir,
                        # then put the proj-dir ones into $fpsout. But there is too much to do,
                        # so for now ignore the custom check-out case.
                        push(@$bfsout, @fset);
                    }
                    else
                    {
                        print STDERR "Warning: Invalid cmd-line file specification $spec.\n";
                    }
                }
                
                # fill in the complete set of files
                if ($cotype eq kCoNetDRMS)
                {
                    push(@$fsout, @$core);
                    push(@$fsout, @$netonly);
                }
                elsif ($cotype eq kCoSdp)
                {
                    push(@$fsout, @$core);
                    push(@$fsout, @$sdponly);
                }
            }
        }
        else
        {
            print STDERR "No support for cmd-line file specification with this check-out type/operation.\n";
            $rv = 1;
        }
    }    
    else
    {
        if ($cotype eq kCoCustom && $dltype eq kDlUpdate && (!defined($xmldata) || $#{$xmldata->{$strproj}} < 0))
        {
            # If this is an update of a custom checkout, and the user did not provide a config-file argument,
            # then use the file spec stored in the existing typefile
            if (defined($stfspec) && length($stfspec) > 0)
            {
                @$fsout = qw($stfspec);
            }
            else
            {
                print STDERR "Invalid file specification in state file.\n";
                $rv = 1;
            }
        }
        else
        {
            # If doing an update, add to the file set the list of deleted files. This will
            # result in removal of files in the sandbox that have already been removed from
            # the repository. This only matters for files directory in the root directory.
            # Files is subdirectories will be processed via a cvs update of the subdirectory.
            if ($cotype eq kCoNetDRMS)
            {
                push(@$fsout, @$core);
                if ($dltype eq kDlUpdate)
                {
                    push(@$fsout, @$coreDel);
                }
                push(@$fsout, @$netonly);
                if ($dltype eq kDlUpdate)
                {
                    push(@$fsout, @$netDel);
                }
                push(@$bfsout, @$fsout);
            }
            elsif ($cotype eq kCoSdp)
            {
                push(@$fsout, @$core);
                if ($dltype eq kDlUpdate)
                {
                    push(@$fsout, @$coreDel);
                }
                push(@$fsout, @$sdponly);
                if ($dltype eq kDlUpdate)
                {
                    push(@$fsout, @$sdpDel);
                }
                push(@$bfsout, @$fsout);
            }
            elsif ($cotype eq kCoCustom)
            {
                push(@$fsout, @$core);
                push(@$fsout, @$netonly);
                push(@$bfsout, @$fsout);
                
                # If there is a -R flag, then it applies only to the checkout of project directories specified in 
                # the config file - so keep these files separate from the rest. @$bfsout contains the files in base, 
                # @$pfsout contains the files in proj.
                
                # Use $xmldata to populate @pfsout and @fsout
                foreach $proj (@{$xmldata->{$strproj}})
                {
                    push(@$fsout, kProjSubdir . "/" . $proj->{$strname}->[0]);
                    push(@$pfsout, kProjSubdir . "/" . $proj->{$strname}->[0]);
                }
            }
            else
            {
                print STDERR "Unsupported checkout type $cotype.\n";
                $rv = 1;
            }
        }
    }

   return $rv;
}

sub IsBadGuy
{
   my($file) = $_[0]; 
   my($badguys) = $_[1];

   my($dir);

   if (defined($badguys->{$file}))
   {
      return 1;
   }
   else
   {
      if (substr($file, -1, 1) eq '/')
      {
         $file = substr($file, 0, length($file) - 1);
      }

      while (1)
      {        
         my($fn, $dir, $sfx) = fileparse($file);

         if (substr($dir, -1, 1) eq '/')
         {
            $dir = substr($dir, 0, length($dir) - 1);
         }

         if (!defined($dir) || length($dir) == 0 || $dir eq ".")
         {
            return 0;
         }

         if (defined($badguys->{$dir}))
         {
            return 1;
         }

         if ($file ne $dir)
         {
            $file = $dir;
         }
         else
         {
            return 0;
         }
      }
   }

   return 0;
}

sub DownloadTree
{
   my($cotype) = $_[0];
   my($dltype) = $_[1];
   my($version) = $_[2];
   my($pversion) = $_[3];
   my($fspec) = $_[4];
   my($bfspec) = $_[5];
   my($pfspec) = $_[6];
    my($cmdlspec) = $_[7];
   my($logfile) = $_[8];
   my($forceco) = $_[9];
   my($netfilter) = $_[10];
   my($sdpfilter) = $_[11];
    my($actualcvs) = $_[12];
    my($actualspec) = $_[13];

   my($rv) = 0;
   my($curdir);
   my($callstat);
   my($cvscmd);
   my($cmd);
   my($rev);
   my($prev);
   my(@relpaths);
   my($forcecostr) = $forceco ? "-f" : "";
   my($filterdir);
   my(@specfiles);
   my($filter);
   my($log);

   if ($dltype eq kDlPrint || $dltype eq kDlPrintRelease)
   {
      # Output should be a list of files only.
      $log = "/dev/null";
   }
   else
   {
      $log = $logfile;
   }

   if (length($dltype) > 0)
   {
      if ($dltype eq kDlCheckout || $dltype eq kDlTag || $dltype eq kDlUntag)
      {
         $cvscmd = "checkout -AP";
      }
      elsif ($dltype eq kDlExport || $dltype eq kDlPrint || $dltype eq kDlPrintRelease)
      {
         $cvscmd = "export";
      }
      elsif ($dltype eq kDlUpdate)
      {
         # If a new directory is added to the repository AND the new directory is added to the
         # file specifications above, then the -d flag to cvs update will cause the new directory
         # to be downloaded to the client.
         $cvscmd = "update -APd";
      }
      else
      {
         print STDERR "Unsupported download type $dltype.\n";
         $rv = 1;
      }
   }
   else
   {
      # Default to a checkout.
      $cvscmd = "checkout -AP";
   }

   if (!$rv)
   {
      if (defined($version) && length($version) > 0)
      {
         $rev = "-r $version";
      } 
      elsif ($dltype eq kDlExport || $dltype eq kDlPrint || $dltype eq kDlPrintRelease)
      {
         # Only export requires a revision argument
         $rev = "-r HEAD";
      } 
      else
      {
         $rev = "";
      }

      if (defined($pversion) && length($pversion) > 0)
      {
         $prev = "-r $pversion";
      } 
      elsif ($dltype eq kDlExport || $dltype eq kDlPrint || $dltype eq kDlPrintRelease)
      {
         # Only export requires a revision argument
         $prev = "-r HEAD";
      } 
      else
      {
         $prev = "";
      }
   }

   if (!$rv)
   {
      # Check the existence of the proper directories required for each operation.
      if ($dltype eq kDlCheckout || $dltype eq kDlExport)
      {
         # If $dltype is kDlCheckout, kDlExport,then there must not be a JSOC subdirectory in the current directory. 
         if (-e $rdir)
         {
            print STDERR "Root directory " . $rdir . " already exists.\n";
            $rv = 1;
         }
      }
      elsif ($dltype eq kDlTag || $dltype eq kDlUntag || $dltype eq kDlPrint || $dltype eq kDlPrintRelease)
      {
         # If $dltype is kDlTag,  kDlUntag, kDlPrint, or kDlPrintRelease then it is okay to delete the JSOC subdirectory 
         # because these three operations create a temporary JSOC directory.
         if (-e kTmpDir . $rdir)
         {
            remove_tree(kTmpDir . $rdir, {error => \my $errlist});
            if (@{$errlist})
            {
               print STDERR "Unable to properly remove temporary subdirectory $tmpdir.\n";
               $rv = 1;
            }
         }
      }
      elsif ($dltype eq kDlUpdate)
      {
          # If $dltype is kDlUpdate, then there MUST be a working directory root.
          # The current directory is the parent of the root directory (if it exists).
          my($rootdir) = File::Spec->catdir($rdir);
          
          if (!(-d $rootdir))
          {
              print STDERR "Working dir is $ENV{PWD}; expected DRMS root dir is $rootdir.\n";
              print STDERR "No CVS working directory exists in $ENV{'PWD'}.\n";
              $rv = 1;
          }
      }
      else
      {
         print STDERR "Unsupported download type $dltype.\n";
         $rv = 1;
      }
   }

   if (!$rv)
   {
      # If this is an export or print command, download to a temporary directory, 
      # obtain the full list of files from the download directory, 
      # apply the filters to remove files from this list, then call CVS again
      # providing this list of file paths. After obtaining this list, 
      # delete the files in the temporary directory.
      #
      # Must cd to the temporary directory because if this is an export, 
      # then the current directory is not the temporary directory.
      if ($dltype eq kDlExport || $dltype eq kDlPrint || $dltype eq kDlPrintRelease)
      {
         $filterdir = $ENV{'PWD'};
         $filter = ($cotype eq kCoSdp) ? $sdpfilter : $netfilter;

         # The temp directory may not exist if this is an export.
         if (!(-d kTmpDir))
         {
            # no need to check this call, because the chdir() cmd is being checked.
            mkpath(kTmpDir);
         }
      }
   }

   if (!$rv)
   {
      # Checkout the project files with a separate cvs command - using the version tag specific to the project
      # files.

      # THIS SECTION CHECKS-OUT THE BASE FILES.
      if (length($rev) > 0)
      {
         $cvscmd = "$cvscmd $forcecostr";
      }

       @relpaths = map({$reprdir . "$_"} @{$bfspec});           
       
       # Filter out undesired files. $filter contains a black list. DOES NOT APPLY TO update.
      if (defined($filterdir))
      {
         # We need to filter out undesirables - checkout into temp dir now (which we cded to above).
         my(%badguys);
         my($tcmd) = join(' ', "cvs", $cvscmd, $rev, @relpaths);

         $rv = (chdir(kTmpDir) == 0);
         if ($rv)
         {
            print STDERR "Unable to cd to " . kTmpDir . ".\n";
         } 
         else
         {
            # Check-out all files from repository (will filter out the non-release files below).
            if (CallCVS($tcmd, $log, undef, 0))
            {
               print STDERR "Unable to $dltype repository files.\n";
               $rv = 1;
            }
            else
            {
               # Get a list of all files from temp dir.
               @specfiles = ();
               if (GetFileList($rdir, "", \@specfiles))
               {
                  print STDERR "Unable to retrieve list of files rooted at " . kTmpDir . $rdir . ".\n";
                  $rv = 1;
               }
               else
               {
                  # At long last we can filter out the bad guys.
                  foreach $guy (@{$filter})
                  {
                     $badguys{$rdir . $guy} = 1;
                  }
                  
                  # Since badguys might contain directories, not just plain files, we have to 
                  # check ALL parent directories for existence in badguys. Blah.
                  @relpaths = map({IsBadGuy($_, \%badguys) ? () : $_;} @specfiles);
               }
            }
         }

         # Back to the original working directory
         remove_tree(kTmpDir . $rdir, {error => \my $errlist});
         if (@{$errlist})
         {
            print STDERR "Unable to properly remove temporary subdirectory $tmpdir.\n";
            $rv = 1;
         }

         chdir($filterdir);
      }

      $cmd = join(' ', $cvscmd, $rev);
       
       if (defined($actualcvs))
       {
           $$actualcvs = $cmd;
       }
       
       if (defined($actualspec))
       {
           @$actualspec = @relpaths;
       }

       if (CallCVSInChunks($cmd, \@relpaths, $log, undef, 0))
      {
         print STDERR "Unable to $dltype repository files.\n";
         $rv = 1;
      }
      else
      {
         # If this is a sdp or net checkout, then $pfspec will be empty, so we can skip this second
         # checkout.

         if (defined($pfspec) && $#{$pfspec} >= 0)
         {
            # THIS SECTION CHECKS-OUT THE PROJECT FILES (if this is a custom checkout).
            if (length($prev) > 0)
            {
               $cvscmd = "$cvscmd $forcecostr";
            }

            @relpaths = map({$reprdir . "$_"} @{$pfspec});

            if (defined($filterdir))
            {
               # We need to filter out undesirables - checkout into temp dir now (which we cded to above).
                # DOES NOT APPLY TO update.
               my(%badguys);
               my($tcmd) = join(' ', "cvs", $cvscmd, $prev, @relpaths);

               $rv = (chdir(kTmpDir) == 0);
               if ($rv)
               {
                  print STDERR "Unable to cd to " . kTmpDir . ".\n";
               } 
               else
               {
                  # Check-out all files from repository (will filter out the non-release files below).
                  if (CallCVS($tcmd, $log, undef, 0))
                  {
                     print STDERR "Unable to $dltype repository files.\n";
                     $rv = 1;
                  }
                  else
                  {
                     # Get a list of all files from temp dir.
                     @specfiles = ();
                     if (GetFileList($rdir, "", \@specfiles))
                     {
                        print STDERR "Unable to retrieve list of files rooted at " . kTmpDir . $rdir . ".\n";
                        $rv = 1;
                     }
                     else
                     {
                        # At long last we can filter out the bad guys.
                        foreach $guy (@{$filter})
                        {
                           $badguys{$rdir . $guy} = 1;
                        }
                        
                        @relpaths = map({IsBadGuy($_, \%badguys) ? () : $_;} @specfiles);
                     }
                  }
               }

               remove_tree(kTmpDir . $rdir, {error => \my $errlist});
               if (@{$errlist})
               {
                  print STDERR "Unable to properly remove temporary subdirectory $tmpdir.\n";
                  $rv = 1;
               }

               chdir($filterdir);
            }

            $cmd = join(' ', "cvs", $cvscmd, $prev, @relpaths);
             
            if (CallCVS($cmd, $log, undef, 0))
            {
               print STDERR "Unable to $dltype repository files.\n";
               $rv = 1;
            }
         }
      }
   }
    
    if (!$rv && $dltype eq kDlUpdate)
    {
        # An update will not remove files that were deleted from the repository from the sandbox.
        # So, we must do a cvs update on every file that was deleted from the server. Shit. 
        # How do we figure out which files were deleted from the server? It's not like
        # we can ask cvs that question? But we can run any old cvs command using the global
        # -n flag, which will forestall any sandbox changes. Call cvs -n update on the root dir.
        # Then parse STDERR results (for some crazy reason cvs prints this output to stderr), 
        # saving all lines with "is no longer in the repository". The text immediately before this
        # string is the name of the file. For all these files, call cvs update.
        my(@resp);
        my(@orig);
        my(@delfiles);
        
        # The working directory is the root directory.
        
        # Unfortunately, we have to call CVS twice - once for base specs and once for proj 
        # specs. Each type of spec will use a different -r flag.
        $cmd = "-n update $rev";

        if (length($rev) > 0)
        {
            $cmd = "$cmd $forcecostr";
        }
        
        # Original file specs.
        if ($#cmdlspec >= 0)
        {
            @orig = map({$reprdir . "$_"} @{$cmdlspec});
        }
        else
        {
            @orig = map({$reprdir . "$_"} @{$bfspec});
        }
        
        if (CallCVSInChunks($cmd, \@orig, undef, \@resp, 1))
        {
            print STDERR "Unable to locate sandbox files that were deleted from the repository.\n";
            $rv = 1;
        }
        else
        {
            push(@delfiles, @resp);

            if ($#cmdlspec < 0)
            {
                # If a cmdlspec was given, then all specs are in cmdlspec.
                $cmd = "-n update $prev";
                if (length($prev) > 0)
                {
                    $cmd = "$cmd $forcecostr";
                }
                
                if (CallCVSInChunks($cmd, \@$pfspec, undef, \@resp, 1))
                {
                    print STDERR "Unable to locate sandbox files that were deleted from the repository.\n";
                    $rv = 1;
                }
                else
                {
                    push(@delfiles, @resp);
                }
            }
        }

        if (!$rv)
        {
            my(@filestoud);
            
            foreach my $line (@delfiles)
            {
                chomp($line);
                if ($line =~ /(\S+)\s+is no longer in the repository/)
                {
                    push(@filestoud, $1);
                }
            }
            
            # Finally, call CVS yet again, but without the "-n" flag. This will remove from the sandbox
            # files deleted from the repository.
            
            my($removecmd) = $cmd;
            $removecmd =~ s/-n//;

            if (CallCVSInChunks($removecmd, \@filestoud, $log, undef, 0))
            {
                print STDERR "Unable to locate sandbox files that were deleted from the repository.\n";
                $rv = 1;
            }
        }
    }
    
   return $rv;
}

# tag all files in cwd with tag $tag
sub TagFiles
{
   my($tag) = $_[0];
   my($dltype) = $_[1];
   my($logfile) = $_[2];

   my($rv) = 0;
   my(@allfiles);
   my($curdir) = $ENV{'PWD'};
   my($cmd);

   if ($dltype eq kDlTag)
   {
      # call cvs tag -c <tag> 
      $cmd = "cvs tag -c $tag ."
   } else
   {
      # call cvs tag -d <tag> 
      $cmd = "cvs tag -d $tag ."
   }
   
   if (CallCVS($cmd, $logfile, undef, 0))
   {
      print STDERR "Unable to tag repository files.\n";
      $rv = 1;
   }

   return $rv;
}

sub PrintFilenames
{
   my($fh) = $_[0];
   my($froot) = $_[1];
    my($sort) = $_[2];
    my($regexp) = $_[3];

   my($rv) = 0;
   my(@allfiles);

   if (GetFileList($froot, "", \@allfiles))
   {
      print STDERR "Unable to retrieve list of files rooted at $froot.\n";
      $rv = 1;
   }
   else
   {
       my(@final);
       
       if ($sort)
       {
           @final = sort(@allfiles);
       }
       else
       {
           @final = @allfiles;
       }
       
       foreach $afile (@final)
       {
           if (defined($regexp))
           {
               if ($afile !~ $regexp)
               {
                   print $fh "$afile\n";
               }
           }
           else
           {
               print $fh "$afile\n";
           }
       }
   }

   return $rv;
}

sub SPrintFilenames
{
    my($listA) = $_[0];
    my($froot) = $_[1];
    
    my($rv) = 0;
    my(@allfiles);
    
    if (GetFileList($froot, "", \@allfiles))
    {
        print STDERR "Unable to retrieve list of files rooted at $froot.\n";
        $rv = 1;
    }
    else
    {
        foreach $afile (@allfiles)
        {
            push(@$listA, "$afile");
        }
    }
    
    return $rv;
}

# Returns a list of all files in the code tree rooted at $spec. The names of each file
# will be prefixed by $dir, unless $dir is the empty string. $listout points to the
# returned list of files.
sub GetFileList
{
   my($spec) = $_[0];
   my($dir) = $_[1]; # the directory to prepend to files when inserting into $listout
   my($listout) = $_[2];

   my($rv) = 0;
   my($prefix);

   if (length($dir) > 0)
   {
      if (substr($dir, length($dir) - 1, 1) eq "/")
      {
         $prefix = $dir;
      }
      else
      {
         $prefix = "$dir/";
      }
   } 
   else
   {
      $prefix = "";
   }

   if (!(-e $spec))
   {
      print STDERR "File $spec does not exist in $ENV{'PWD'}.\n";
      $rv = 1;
   }
   else
   {

      if (-f "$spec")
      {
         # $spec is a file - just push onto output list
         push(@{$listout}, "$prefix$spec");
      }
      elsif (-d "$spec")
      {
         my(@alltreefiles);
         my(@treefiles);
         my($curdir);

         # This is a directory, collect all files (excluding "." and "..") in the tree rooted at
         # this directory.
         $curdir = "$ENV{'PWD'}";
         chdir($spec);
         tie(my(%tree), "IO::Dir", ".");
         @alltreefiles = keys(%tree);
         @treefiles = map({$_ !~ /^\.$/ && $_ !~ /^\.\.$/ ? $_ : ()} @alltreefiles);


         # Now recursively call GetFileList() for each item in @treefiles
         foreach $childspec (@treefiles)
         {
            GetFileList($childspec, "$prefix$spec", $listout);
         }

         chdir($curdir);
      }
      else
      {
         print STDERR "File '$spec' is not a supported file type.\n";
         $rv = 1;
      }
   }

   return $rv;
}

sub CallCVS
{
    my($cmd) = $_[0];
    my($log) = $_[1];
    my($rsp) = $_[2];
    my($silent) = $_[3];
    
    my($rv) = 0;
    my($callstat);
    
    if (defined($log) && length($log) > 0)
    {
        system("$cmd 1>$log 2>&1");
    }
    elsif (defined($rsp))
    {
        if (open(PIPE, "$cmd 2>&1 |"))
        {
            @$rsp = <PIPE>;
            close(PIPE);
        }
        else
        {
            print STDERR "Could not run '$cmd'\n";
            $rv = 1;
        }
        
    }
    elsif ($silent)
    {
        system("$cmd 1>/dev/null 2>&1");
    }
    else
    {
        system($cmd);
    }
    
    $callstat = $?;
    
    if ($callstat == -1)
    {
        print STDERR "Failed to execute '$cmd'.\n";
        $rv = 1;
    }
    elsif ($callstat & 127)
    {
        print STDERR "cvs command terminated abnormally.\n";
        $rv = 1;
    }
    elsif ($callstat >> 8 != 0)
    {
        $callstat = $callstat >> 8;
        print STDERR "cvs command ran unsuccessfully, status == $callstat.\n";
        $rv = 1;
    }
    
    return $rv;
}

sub CallCVSInChunks
{
    my($cmd) = $_[0];
    my($specs) = $_[1];
    my($log) = $_[2];
    my($rsp) = $_[3];
    my($silent) = $_[4];
    
    my(@chunk);
    my(@checkout);
    my($fullcmd);
    my($rv);
    
    $rv = 0;

    # Submit with kMaxFileSpecs file specs at most.
    foreach my $aspec (@$specs)
    {
        # %@$&*(~ CVS! If $aspec is part of a path that does not exist in the current working directory, 
        # CVS will choke if the CVS command is not a checkout or export command. The work-around in this case is to run
        # cvs checkout -A first. But you cannot run this command on a directory, no that would be too easy.
        # If you do that CVS will checkout all files in the directory. We have to instead checkout a file in 
        # the directory that we want to really have in the working directory after CallCVSInChunks completes.
        # That will result in the checkout of just the file we want to download and the creation of the path 
        # leading to the file. 
        #
        # So, if $aspec is a file that doesn't exist, first run cvs checkout -A on that file. But only do this
        # if $cmd is not checkout or export!

        if ($cmd !~ /^\s*checkout/ && $cmd !~ /^\s*co/ && $cmd !~ /^export/)
        {
            if (!(-e $aspec))
            {
                push(@checkout, $aspec);
            }
        }

        push(@chunk, $aspec);
        
        if ($#chunk == &kMaxFileSpecs - 1)
        {
            # Checkout files that do not exist first.
            if ($#checkout >= 0)
            {
                $fullcmd = join(' ', "cvs checkout -A", @checkout);
                if (CallCVS($fullcmd, $log, $rsp, $silent))
                {
                    $rv = 1;
                    last;
                }

                # Must remove the files just checked out or else the next call to cvs might fail. But
                # do not remove the temp dir we just downloaded files into.
                remove_tree(&kTmpDir, {keep_root => 1});
                
                @checkout = ();
            }

            $fullcmd = join(' ', "cvs", $cmd, @chunk);
            if (CallCVS($fullcmd, $log, $rsp, $silent))
            {
                $rv = 1;
                last;
            }
            
            @chunk = ();
        }
    }
    
    if (!$rv && $#chunk >= 0)
    {
        if ($#checkout >= 0)
        {
            $fullcmd = join(' ', "cvs checkout -A", @checkout);
            if (CallCVS($fullcmd, $log, $rsp, $silent))
            {
                $rv = 1;
            }

            # Must remove the files just checked out or else the next call to cvs might fail. But
            # do not remove the temp dir we just downloaded files into.
            remove_tree(&kTmpDir, {keep_root => 1});
            
            @checkout = ();
        }

        $fullcmd = join(' ', "cvs", $cmd, @chunk);
        if (CallCVS($fullcmd, $log, $rsp, $silent))
        {
            $rv = 1;
        }
     
        @chunk = ();
    }
 
    return $rv;
}

sub MoveFiles
{
   my($files) = $_[0];
   my($srcroot) = $_[1];
   my($tgtroot) = $_[2];

   my($rv) = 0;
   my($fullsrc);
   my($fulltgt);
   my($tgtdir);

   foreach $onefile (@{$files})
   {
      $fullsrc = "$srcroot/$onefile";
      $fulltgt = "$tgtroot/$onefile";
      $tgtdir = dirname($fulltgt);

      if (-e $fullsrc)
      {
         # Make sure the tgt directory exists.
         if (!(-d $tgtdir))
         {
            mkpath($tgtdir);
         }

         if (-d $tgtdir)
         {
            if (!move($fullsrc, $fulltgt))
            {
               print STDERR "Unable to move $fullsrc to $fulltgt; skipping.\n";
            }
         }
         else
         {
            print STDERR "Unable to make target directory $tgtdir; skipping.\n";
         }
      }
      else
      {
         print STDERR "File $fullsrc does not exist; skipping.\n";
      }
   }

   return $rv;
}

sub InsertFile
{
    my($root) = shift;
    my($afile) = shift; 
    my($hash) = $_[0];
    
    my($top);
    my($bot);
    
    if ($afile =~ /^([^\/]+)\/(.+)$/)
    {
        $top = $1;
        $bot = $2;
        
        if (!exists($hash->{$top}))
        {
            $hash->{$top} = {};
        }
        
        if (length($root) > 0)
        {
            InsertFile("$root/$top", $bot, $hash->{$top});
        }
        else
        {
            InsertFile($top, $bot, $hash->{$top});
        }
    }
    else
    {
        $hash->{$afile} = $root;
    }
}

sub GetDefinitiveFileSet
{
    my($cotype) = $_[0];
    my($core) = $_[1];
    my($netonly) = $_[2];
    my($sdponly) = $_[3];
    my($xmldata) = $_[4];
    my($stfspec) = $_[5];
    my($version) = $_[6];
    my($pversion) = $_[7];
    my($reprdir) = $_[8];
    my($rdir) = $_[9];
    
    my($sdir);
    my(@dlist);
    my(@flist);
    my(@sorted);
    my($listf);
    my($cmdlcvs);
    my(%seen);
    my($ierr);
    my($rv);
    
    $rv = {};
    
    if (!(-d &kTmpDir))
    {
        # no need to check this call, because the chdir() cmd is being checked.
        mkpath(kTmpDir);
    }
    
    $sdir = $ENV{'PWD'};
    if (chdir(kTmpDir) == 0)
    {
        print STDERR "Unable to cd to " . kTmpDir . ".\n";
        $ierr = 1;
    }
    else
    {
        my(@allspecs);
        my(@wroot);

        $cmdlcvs = "cvs export ";
        
        if ($cotype eq &kCoSdp || $cotype eq &kCoNetDRMS)
        {
            if (length($version) > 0)
            {
                $cmdlcvs = $cmdlcvs . "-r $version ";
            }
            else
            {
                $cmdlcvs = $cmdlcvs . "-r HEAD ";
            }
            
            if ($cotype eq &kCoSdp)
            {
                @wroot = map({$reprdir . "$_"} @{$core});
                push(@allspecs, @wroot);
                @wroot = map({$reprdir . "$_"} @{$sdponly});
                push(@allspecs, @wroot);

                $cmdlcvs = join(' ', $cmdlcvs, @allspecs);

            }
            else
            {
                @wroot = map({$reprdir . "$_"} @{$core});
                push(@allspecs, @wroot);
                @wroot = map({$reprdir . "$_"} @{$sdponly});
                push(@allspecs, @wroot);

                $cmdlcvs = join(' ', $cmdlcvs, @allspecs);
            }
        }
        else
        {
            # ??
            print STDERR "Can't specify cmd-line file spec for custom check-out type.\n";
            $ierr = 1;
        }
        
        if (!$ierr)
        {
            unless (CallCVS($cmdlcvs, undef, undef, 1))
            {
                unless (GetFileList(&kTmpDir . $rdir, "", \@dlist))
                {
                    my($prefdiri) = &kTmpDir . $rdir;

                    foreach my $afile (@dlist)
                    {
                        # remove &kTmpDir . &kRootDir
                        $afile =~ s/$prefdiri//;
                        push(@flist, "$afile");
                    }
                    
                    # Sort and extract list of unique file names.
                    @sorted = sort(@flist);
                    @listf = map({ unless ($seen{$_}++){($_)}else{()} } @sorted);
                    
                    foreach my $afile (@listf)
                    {
                        # Insert into $rv.
                        InsertFile("", $afile, $rv);
                    }
                }
                else
                {
                    $ierr = 1;
                }

                remove_tree(&kTmpDir);
            }
            else
            {
                $ierr = 1;
            }
            
            chdir($sdir);
        }
    }

    return $rv;
}

sub CollectFiles
{
    my($spec) = $_[0];
    my($deflistH) = $_[1];
    
    my($ierr);
    my(@rv);
    
    $ierr = 0;

    # Peel off the top directory name (must be relative to jsocroot).
    if ($spec =~ /^([^\/]+)\/(.+)$/)
    {
        my($top) = $1;
        my($bottom) = $2;
        my($listH) = $deflistH->{$top};
        
        if (ref($listH))
        {
            @rv = CollectFiles($bottom, $listH);
        }
        else
        {
            print STDERR "The directory $top should have a hash array associated with it.\n";
            $ierr = 1;
        }
    }
    else
    {
        my($listH);
        
        # $spec could have a trailing '/'.
        if ($spec =~ /\s*(.+)\/\s*$/)
        {
            $spec = $1;
        }
        
        # This is a node - either a directory or a regular file.
        $listH = $deflistH->{$spec};
        
        if (!defined($listH))
        {
            print STDERR "Unknown file $spec.\n";
            $ierr = 1;
        }
        else
        {
            if (ref($listH))
            {
                # Directory.
                my(@dirfiles);
                
                foreach my $elem (keys(%$listH))
                {
                    @dirfiles = CollectFiles($elem, $listH);
                    push(@rv, @dirfiles);
                }
            }
            else
            {
                # Regular file.
                @rv = $listH . "/" . $spec;
            }
        }
    }

    return @rv;
}


__DATA__
I think we need data here in avx.
__END__
