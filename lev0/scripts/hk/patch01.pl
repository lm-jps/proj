#!/usr/bin/perl
################################################################################
# NAME: patch01.pl
# AUTHOR: carl
# DESCRIPTION: Provides patch for HK Configuration files  AIA-02C-version-1.163 and
#              HMI-010-version-1.163 in 1.163 directory, HMI-010-version-1.162
#              in directory 1.162, and HMI-010-version-1.161 in directory 1.161.
#              Fixes  the following items for bug #261 in TRAC bug reporting tool:
#              --->"IS1 A" in config files to "IU1 A"
#              --->"IL1 A" in config files to "UL1 A"
#              --->" 0 6  -1" to " 0 6  1"
# Execution: patch01.pl
# LIMITATION:
# --Hard coded filename shown below(AIA-02C-version-1.163,HMI-010-version-1.163).
# --Set correct directory to file HK configuration files to update
#   (i.e., $hm/cvs/TBL_JSOC/lev0/hk_config_file)
# --To use for other version like 1.164. Just update filename(i.e., 
#   AIA-02C-version-1.164,HMI-010-version-1.164) in push command below. May want
#   to rename to patch02.pl and update comment on patch in script.
# CREATED:03/05/2010
################################################################################
 
#set environment variables and initialize variables.
$hm=$ENV{'HOME'};

#check for any arguments passed in command
&check_arguments();

# HK environment variable for processing STANFORD file
$cf_dir=$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";

# HK Config files to Patch-if need to do other then replace or add file names below
push (@filelist,"HMI-010-version-1.161");
push (@filelist,"HMI-010-version-1.162");
push (@filelist,"AIA-02C-version-1.163");
push (@filelist,"HMI-010-version-1.163");

#loop through each file to patch and update
foreach $f (@filelist)
{

  #set version directory using filename's file version number
  @s_filename= split("-",$f);
  $cf_ver= $s_filename[3];

  #print to standard out message
  print "\n-->Patching file <$cf_dir/$cf_ver/$f> with patch01.pl script\n";

  #copy orig file to BACKUP_PATCH01_<filename>
  $log=`cp $cf_dir/$cf_ver/$f  $cf_dir/$cf_ver/BACKUP-PATCH01-$f`;

  # open copied file for reading and search for items to patch in while loop
  open(INFILE, "+<$cf_dir/$cf_ver/BACKUP-PATCH01-$f") || die "Can't Open: <$cf_dir/$cf_ver/BACKUP-PATCH01-$f> file: $!\n";

  # open actual filename to write updated patch lines  and unchanged line to.(note if need orig file content get from cvs)
  open(OUTFILE, ">$cf_dir/$cf_ver/$f") || die "Can't Open: <$cf_dir/$cf_ver/$f> file: $!\n";

  #set count of number of lines patches- for patch01 should have 32 lines to AIA-02C-version-1.163 and 16 lines in HMI-010-version-1.163
  $count=0;
  while (<INFILE> ) 
  {
    #remove cr from printed line in file only not actual line
    $prtline=$_;
    $prtline=~s/\n//;

    #check if got three patterns to patch or replace and print out to screen
    if ( /IS1 A  20070904/ ) 
    {     
        print "-->Found line with substring <IS1 A 20070904>.\n-->Patching line <$prtline>.\n-->Replacing with substring <IU1 A  20100303>\n";
        $count++;
    }
    if ( /IL1 A  20070904/ ) 
    {     
        print "-->Found line with substring <IL1 A 20070904>.\n-->Patching line <$prtline>.\n-->Replacing with substring <UL1 A  20100303>\n";
        $count++;
    }
    if ( / 6 0  -1 0  0  0  0  20070904/ ) 
    {     
        print "-->Found line with substring < 6 0  -1 0  0  0  0  20070904>.\n-->Patching line <$prtline>.\n-->Replacing with substring < 6 0  1 0  0  0  0  20100303>\n";
        $count++;
    }
    if ( /script version/ )
    {
       if ( /A PATCH01/ )
       {
         ;#skip-write patch note once
       }
       else
       {
         #add comment on patch and date- use prtline variable since has no cr to rewrite version line in file with patch notes
         $datepatch=`date`;
         $datepatch=~s/\n//;
         $commentline= sprintf("%s (%s%s)%s", $prtline,"NOTE: A PATCH01 was applied using patch01.pl script on:",$datepatch,"\n");
         $_= $commentline;
       }
    }

     #replace here the pattern with patch value
     s/IS1 A  20070904/IU1 A  20100303/; 
     s/IL1 A  20070904/UL1 A  20070303/; 
     s/ 6 0  -1 0  0  0  0  20070904/ 6 0  1 0  0  0  0  20100303/;

     # print updated line to updated file or just print unchanged lines to updated file
     print OUTFILE "$_";
  }
  #print number of lines patched!
  print "Patched <$count> lines in file <$cf_dir/$cf_ver/$f>\n";

  #close in and out files and do next file in foreach list
  close INFILE;
  close OUTFILE;

}
#print completed line
print "\n-->COMPLETE run. View above data for patches completed.\n";
#end main functure

#############################################################################
# subroutine check arguments and set flags                                  #
#############################################################################
sub check_arguments()
{
  $help_flg = "0";
  if ($#ARGV >= 0 && $#ARGV <= 1)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    else
    {
       print "-->ERROR:Too many arguments or bad arguments. View usage below.\n";
       $help_flg = "1";
    }
  }
  if ( $help_flg eq "1")
  {
    show_help_info();
  }
}
#############################################################################
# subroutine show_help_info: show help information                          #
#############################################################################
sub show_help_info
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Run patch01.pl perl script and direct change messages to screen:\n\n";
  print "         patch01.pl\n\n";
  print "(1b)Run patch01.pl perl scripit and redirect change messages to log file:\n\n";
  print "         patch01.pl > LOG-MESSAGE-PATCH01\n\n";
  print "(1c)Get Help Information:\n\n";
  print "         patch01.pl -h  or  patch01.pl -help\n\n";
  print "(2)****REQUIRED ITEMS TO SET****\n";
  print "(2a) Required that the hk config directory environment variable are setup correctly at top of script(HK_CONFIG_DIRECTORY).\n";
  print "(2c) Required that we add push command for HK Config files on list of files to process. Currently we process only \n";
  print "     AIA-02C-version-1.163 and HMI-010-version-1.163 files for this patch01.\n";
  print "(2d) Required that have access to directories and files used by HK_CONFIG_DIRECTORY set in script.\n";
  print "(3)****Limitation****:\n";
  print "(3a)Need to Hard code filenames to process(i.e.,AIA-02C-version-1.163,HMI-010-version-1.163).\n";
  print "(3b)Need to set correct directory to HK configuration files to update(HK_CONFIG_DIRECTORY).\n";
  print "(3c)To use for other file versions like 1.164. Probably copy patch01.pl to patch02.pl. Then change filenames\n";
  print "    in push commands in script(i.e., push (filelist,AIA-02C-version-1.164),push(filename,HMI-010-version-1.164).\n";
  print "    Finally update patch note in script.\n";
  print "(3d)Patterns to patch are hardcoded in script. If pattern changes will not work(see \"IS1 A\" and \"IL1 A\" in this script).\n";
  print "    May need to update patten if doing another file version or another patch fix.\n";
  exit;
}
