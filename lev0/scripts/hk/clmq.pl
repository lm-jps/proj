#!/usr/bin/perl
##############################################################################
# Name:        clmq.pl - Check Lockheed Martin Queue file                    #
# Description: Used in cron job on production to check queue file for new    #
#              STANFORD and GROUND files from LMSAL. If new files run        #
#              scripts to bring over STANFORD and GROUND files from solserv  #
#              server. Build HKPDF and JSD files. Checkin to cvs all files.  #
#              If needed create new JSD in DRMS based on list of needed apids#
#              Send email status to jsoc_ops team.                           #
# Execution:   (1)To run standalone:clmq.pl                                  #
#              (2)To run in cron is normal way to run. See example below.    #
#              (3)For help: clmq.pl -h                                       #
# Example      Cron:0,30 9-18 * * 1-5 /usr/bin/perl                          #
#               /home1/carl/cvs/JSOC/proj/lev0/scripts/hk/clmq.pl >          #
#               /dev/null 2>&1                                               #
# Example:     % clmuq.pl                                                    #
# Limitation:  Setup required to ssh without passwords to environment running#
#              DRMS libraies                                                 #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008 #
##############################################################################
##get environment variables and initialize variables.
$hm=$ENV{'HOME'};
$mailfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/dailyStatusEmail.txt";
$tar_src_dir="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$lmq_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$lmq_file="UpdateLMQueue.txt";
$lmf_dir="$hm/cvs/TBL_JSOC/lev0/fromlmsal";
$solserv_dir="/home/jim/stanford-files";

#common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";
$logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-clmuq";
$ENV{'CVS_RSH'}="ssh";

# open log file
open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
print LF `date`;

# set environment variables for cron job executing this script based on machine
$account=`whoami`;
chop $account;
if ( $account eq "carl")
{
  print LF "-->working in <$account> account\n";
  $ENV{'SSH_AUTH_SOCK'}="/tmp/ssh-XXSCh2du/agent.1132";
  $ENV{'SSH_ASKPASS'}="/usr/libexec/openssh/gnome-ssh-askpass";
}
elsif ( $account eq "production")
{
   print LF "-->working in <$account> account\n";
   print LF "ERROR: Found account name <$account>- Need to update clmq.pl script with SSH settings.\n";
   print LF `date`;
   die "ERROR: Found account name <$account>- Need to update clmq.pl script with SSH settings";
}
else
{
   print LF "ERROR: did not find this account name <$account>- Need to update clmq.pl script.\n";
   print LF `date`;
   die "ERROR: did not find this account name <$account>- Need to update script.";
}

# open mail text file to send to users of hk config data
open(MF,">$mailfile") || die "Can't Open $mailfile: $!\n";

print LF "-->accout home directory is <$hm>\n";
print LF "-->check if time to do update\n";

##get directory and filename for queue file
$dn=$lmq_dir;
$fn=$lmq_file;

##check Queue file has work to do by getting file on solserv
print LF "-->check for work to do in Queue file on solserv machine\n";
$lm=`scp -o user=jim solserv.lmsal.com:$solserv_dir\/$fn $lmq_dir`;

##clear Q file on solserv to confirm have work item.
$lm=`ssh -l jim solserv.lmsal.com \"echo \"\" > $solserv_dir/$fn\" `;

#open file and load all lines
open(QFILE, "$dn/$fn") || die "Can't Open $dn/$fn file: $!\n";

##loop in Q file and process each message
$uflag="n";
while (<QFILE>)
{
  ## set update flag. This is use to print correct message to log
  $uflag="y";

  ##push lines in Qfile
  push(@all_lines, $_) ;

  ##log
  print LF "-->Work to do, line in Q file is  $_";
  print LF "-->Doing update for <$account> account \n";

  ##get latest GTCID and STHA files based on file ver# from solserv server
  $pos = index $_, " ";
  $sfline=substr($_,0, $pos); 
  $gfline=substr($_, $pos + 1); 
  $pos= index  $sfline, "=";
  $sfname=substr($sfline,$pos + 1 ); 
  $pos= index  $gfline, "=";
  $gfname=substr($gfline,$pos + 1); 
  chomp $gfname; #knock off cr at end of line
  print LF "-->Filenames are $gfname and  $sfname\n";

  # scp STHA file and GTCIDS file using filenames from solserv to dimp
  $lm=`scp -o user=jim solserv.lmsal.com:$solserv_dir/$gfname $lmf_dir`;
  print LF "-->Secure copy $gfname finished. $lm\n";
  $lm=`scp -o user=jim solserv.lmsal.com:$solserv_dir/$sfname $lmf_dir`;
  print LF "-->Secure copy $sfname finished. $lm\n";

  # chmod on STHA and GTCID files
  $lm=`chmod 777 $lmf_dir/$sfname`;
  $lm=`chmod 777 $lmf_dir/$gfname`;
  print LF "-->completed chmod of stha and gtcid files. $lm\n";

  ##log email text
  print MF `date`;
  print MF "-->Doing update for Stanford <$account> account.\n";
  #Pprint MF "-->Doing update for Stanford production account on ?yeti?\n";
  print MF "-->Doing updates based on $_";

  ##execute script to build hkpdf files, checkin files, update hmi0
  print LF "-->Running cicf.pl - check in configuration files for $_\n";
  $lm=`perl $script_dir/cicf.pl $_ `;
  print LF "-->Log for cicf: $lm\n";

  ##log email text
  print MF "-->Finished update for Stanford <$account> account\n";
 
  ## get stanford file version number to name of tar file
  $pos= index $sfname, "-";
  $sfvn=substr($sfname,$pos + 1);
  print LF "-->Stanford file version number is $sfvn\n";

  ## tar up HKPDF files and gtcids file
  $lm=`cd $tar_src_dir ;tar -cvf TAR.$sfvn $sfvn/apid* gtcids.txt`;
  print LF "-->tar up files completed. Log:\n$lm\n";

  ## secure copy files to solserv for distribution via hmifsw1 to EGSE machines
  #$lm=`scp -o user=jim  $tar_src_dir/TAR.$sfvn solserv.lmsal.com:$solserv_dir/.`;
  #print LF "-->Secure copied $tar_src_dir/TAR.$sfvn file to solserv. $lm\n";
  #print MF "-->Distributed data files to solserv server. Ready to move tar file to outside machine if need\n";
  #print MF `date`;
  # don't need to send to solserv,see above, since EGSE scripts aleady do this.
  print MF "-->Ready to move tar file to outside machine if needed. TAR file located at $tar_src_dir/TAR.$sfvn\n";
  print MF "-->Or can do a cvs update to get hk configuration files checked into JSOC CVS.\n";
  print MF `date`;

  ## send mail status on update
  $lm=`mail -s "For carl\@yeti only: update completed for hk config files" carl\@sun.stanford.edu < $mailfile`;
 ##P$lm=`mail -s "For production\@yeti only: update completed for hk config files" jsoc_ops\@sun.stanford.edu < $mailfile`;

  ##execute script to build of jsd, jsvn map files, create series if needed 
  ## get ground file version number to name of tar file
  $pos= index $gfname, "-";
  $gfvn=substr($gfname,$pos + 1);
  print LF "-->Ground file version number is $gfvn\n";
  print LF "-->STANFORD file version number is $sfvn\n";

  #set gfname to pass to mhds.pl 
  print LF "-->Running make data series script... \n";
  $lm=`cd $script_dir;perl mhds.pl stha=$sfvn gtcids=$gfvn`;
  print LF "-->Log from mhds:\n$lm\n";

}
if ($uflag eq "n")
{
  print LF "-->No work to do in Q file.\n";
}
close QFILE;
#set Q file to blank work was completed
open(QFILE, ">$dn/$fn") || die "Can't Open $dn/$fn file: $!\n";
close QFILE;
print LF `date`;
close LF;
close MF;
