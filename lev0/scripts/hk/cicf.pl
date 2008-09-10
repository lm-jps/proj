#!/usr/bin/perl
#############################################################################
#Name: cicf.pl - Check in hk configuration files                            #
#Description: Checks into cvs Stanford and Ground files from LMSAL. Builds  #
#             so called housekeeping packet data format(HKPDF) files. Checks#
#             into CVS automatically in folder based on file version name,  #
#             the HKPDF files. Checks in the Ground file as gtcids.txt.     #    
#             These files are used by Level 0 C code and scripts to decode  #
#             housekeeping(HK) keyword values.                              #
#Execution:   To run :                                                      #
#             cicf.pl stha=<file with version #>                            #
#                     gtcids=<file with version number>                     #
#             For help: cicf.pl -h                                          #
#Example:     cicf.pl stha=STANFORD_TLM_HMI_AIA.txt-1.111                   #
#                 gtcids=GROUND_to_CODE_ids.txt-1.103                       #
#Limitationr: User should be careful entering arguments since these are the #
#             the files that will get checked in and updated on production  #
#Author:      Carl                                                          #
#Date:        Move from EGSE to JSOC software environment on March 21, 2008 #
#############################################################################
# main program                                                              #
# set up environment variables for cron job executing this script via clmq.pl
$hm=$ENV{'HOME'};
$ENV{'CVSROOT'}=":ext:sunroom.stanford.edu:/home/cvsuser/cvsroot";
$ENV{'CVS_BINARY_ROOT'}="/home/cvsuser/cvs_binary_root";
$lmf_dir="$hm/cvs/TBL_JSOC/lev0/fromlmsal";
$tar_src_dir="$hm/cvs/TBL_JSOC/lev0/hk_config_file";

#common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";
$logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-clmuq";
$ENV{'CVS_RSH'}="ssh";
#$ENV{'CVS_RSH'}="/usr/bin/rsh";

# set environment variables for cron job executing this script based on machine
$account=`whoami`;
chop $account;
if ( $account eq "carl")
{
  print LF "-->Found <$account> account\n";
  $ENV{'SSH_AUTH_SOCK'}="/tmp/ssh-XXSCh2du/agent.1132";
  $ENV{'SSH_ASKPASS'}="/usr/libexec/openssh/gnome-ssh-askpass";
}
elsif ( $account eq "production")
{
   print LF "-->Found <$account> account\n";
   $ENV{'SSH_ASKPASS'}="/usr/libexec/openssh/gnome-ssh-askpass";
   print LF "-->working in <$account> account\n";
   print LF "-->SSH setting for SSH_ASKPASS is <$ENV{'SSH_ASKPASS'}>\n";
}
else
{
   print LF "ERROR: did not find this account name <$account>- Need to update cicf.pl script.\n";
   print LF `date`;
   die "ERROR: did not find this account name <$account>- Need to update cicf.pl script.";
}


# show environment variable settings
print "-->Environment for cicf : \n";
print `env`;

#check for any arguments passed in command
&check_arguments();
print "\n-->Files to process are : $ARGV[0]  and $ARGV[1]\n" ;

#get stha and gtcids file
$new_stha=substr($ARGV[0],5) ;
$new_gtcids=substr($ARGV[1],7);

#copy STANFORD file
print "-->Copying  $lmf_dir/$new_stha to $lmf_dir/STANFORD_TLM_HMI_AIA.txt file\n" ;
$lm=`rm  -f $lmf_dir/STANFORD/STANFORD_TLM_HMI_AIA.txt`;
$lm=`cp  -f $lmf_dir/$new_stha $lmf_dir/STANFORD_TLM_HMI_AIA.txt`;
$lm=`chmod 777 $lmf_dir/STANFORD_TLM_HMI_AIA.txt`;

#copy GROUND file
print "-->Copying $lmf_dir/$new_gtcids to  $lmf_dir/GROUND_to_CODE_ids.txt file\n";
$lm=`rm -f $lmf_dir/GSE/GROUND_to_CODE_ids.txt`;
$lm=`cp  -f $lmf_dir/$new_gtcids    $lmf_dir/GROUND_to_CODE_ids.txt`;
$lm=`chmod 777 $lmf_dir/GROUND_to_CODE_ids.txt`;

#copy gtcidsi.txt file to lev0 hk configuration files directory
print "-->Copying $lmf_dir/$new_gtcids to $tar_src_dir/gtcids.txt\n";
$lm=`rm -f $tar_src_dir/gtcids.txt`;
$lm=`cp  -f  $lmf_dir/$new_gtcids  $tar_src_dir/gtcids.txt`;
$lm=`chmod 777 $tar_src_dir/gtcids.txt`;

#run script to build Housekeeping Packet Data format(HKPDF) files
print "-->Building HKPDF files:\n";
$arg1= "/usr/bin/perl" ;
$arg2= " make_hkpdf.pl ";
#$arg3= " 4 ";
#@allarg= ("$arg1","$arg2","$arg3");
# new make_hkpdf.pl sort=4 apidlist="ALL_HK"
$arg3= " sort=4 ";
$arg4= " apidlist=\"ALL_HK\" ";
@allarg= ("$arg1","$arg2","$arg3","$arg4");
print "@allarg \n";
$lm=`cd $script_dir;  @allarg`; 
print "$lm\n";
   
#get gtcids's fvn and stha's fvn
$pos= index $ARGV[0], "-";
$sfvn=substr($ARGV[0],$pos + 1);
$pos= index $ARGV[1], "-";
$gfvn=substr($ARGV[1],$pos + 1);

#check in files
$lm=`cd $tar_src_dir/. ;  cvs commit  -m \"auto checked into CVS with version $gfvn\"  gtcids.txt`;
print "-->Committing into CVS gtcids.txt file version $gfvn \n$lm\n";
$lm=`cd   $lmf_dir/. ; cvs commit -m \"auto checked into CVS with version $gfvn\"  -r$gfvn  GROUND_to_CODE_ids.txt`;
print "-->Committing into CVS GROUND_to_CODE_ids.txt file version $gfvn \n$lm\n";
$lm=`cd  $lmf_dir/.  ; cvs commit -m \"auto checked into CVS with version $sfvn\" -r$sfvn  STANFORD_TLM_HMI_AIA.txt`;
print "-->Commiting into CVS STANFORD_TLM_HMI_AIA.txt version $sfvn file \n$lm\n";
$lm=`cd  $tar_src_dir/. ; cvs add $sfvn`;
print "-->Adding $sfvn directory to $tar_src_dir directory \n$lm\n";
$lm=`cd  $tar_src_dir/$sfvn/. ; cvs add apid* `;
print "-->Adding HKPDF files to directory to $tar_src_dir/$sfvn \n$lm\n";
$lm=`cd  $tar_src_dir/$sfvn/. ; cvs commit -m \"auto checked into CVS with build files using STHA file version $sfvn\"  apid* `;
print "-->Commiting HKPDF files in directory $tar_src_dir/$sfvn \n$lm\n";

#############################################################################
# subroutine check arguments and set flags                                  #
#############################################################################
sub check_arguments()
{
  $help_flg= "0";
  $stha_flg="0";
  $gtcids_flg="0";
  
  if ($#ARGV < 0)
  {
    $help_flg="1";
  }
    
  if ($#ARGV >= 0)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    elsif (substr($ARGV[0],0,5) eq "stha=" )
    {
      #use file to get list of apids to create map files for
      $stha_flg = "1";
    }
  }
  if ($#ARGV >= 1)
  {
    if ( substr($ARGV[1],0,7) eq "gtcids=") 
    {
      #use file to get list of apids to create map files for
      $gtcids_flg = "1";
    }
    else
    {
      print "Warning:arguments entered not correct\n";
    }
  }

  if ( $help_flg eq "1")
  {
     &show_help_info;
  }
}
#############################################################################
# subroutine show_help_info: show help information                          #
#############################################################################
sub show_help_info
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Execute to build new HKPDF files and check into CVS the HKPDF, STHA ,GTCIDS, gtcids.txt files:\n
             cicf.pl stha=<STANFORD file to process> gtcids=<GROUND file to process>:\n";
  print "             Example:cicf.pl stha=STANFORD_TLM_HMI_AIA.txt-1.109 gtcids=GROUND_to_CODE_ids.txt-1.101\n\n";
  print "(1b)Help Information: cicf.pl -H \n"; 
  print "             Example: cicf.pl -H or cicf.pl -h\n\n";
  print "(2) Requires using send_GTCIDS_File.pl and send_STHA_File.pl scripts on LM CVS machine. \n";
  print "(3) The send scripts will save files in correct directories with file version number tag at end. \n";
  print "(4) Requires using file version numbers like STANFORD_TLM_HMI_AIA.txt-1.109 & GROUND_to_CODE_ids.txt-1.101 as input\n";
  print "(5) Review Log or printed output for errors \n";
  exit;
}
