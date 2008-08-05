#!/usr/bin/perl
##############################################################################
# Name:        getdf.pl  - get day files, ingest into dayfile data series    #
#                          then remove file                                  #
# Description: Get dayfiles from drop directory for hsb, egsefm(in future)   #
#              Call ingest_dayfile.pl and valid loaded dayfile . If so       #
#              remove dayfile from file system.                              #
# Execution:   getdf.pl hsb                                                  #
# Limitation:   The move-and-ingest-dayfiles process works for dayfiles from #
#               the hsb dayfiles only. This script is not used               #
#               currently to move and  call ingest_dayfiles.pl for           #
#               LMSAL dayfiles.                                              #
##############################################################################
# set Environment Variables

# drop off directory for HSB dayfiles
$doff_dir=$ENV{'DF_DROPOFF_HSB_FILES'}="/tmp21/production/lev0/hk_hsb_dayfile";

# debug flag
$dflg=$ENV{'DF_GETDF_DEBUG'}="0";

#common setting for all environments
$ENV{'SUMSERVER'}="d02.Stanford.EDU";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$ENV{'DF_DRMS_EXECUTABLES'}="$hm/cvs/JSOC/bin/linux_x86_64";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$ENV{'DF_DRMS_EXECUTABLES'}";

# check arguments
&check_agruments();

# set log file based on sdo, moc, or egsefm
if ($src eq "hsb")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-hsb";
}
elsif ($src eq "moc")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-moc";
}
elsif ($src eq "egsefm")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-egsefm";
}

# open log file
open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
print LF `date`;
print LF "--->starting script getdf.pl\n";
print LF "--->processing day files at directory:$doff_dir\n";
close LF;

#ingest all dayfiles files
# for production
if  ($src eq "hsb")
{
  $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl  apidlist=$script_dir/df_apid_list_day_file_hsb start=20080801 end=20080830 dsnlist=$script_dir/df_apid_ds_list_for_hsb src=hsb`;
  print LF "--->Processed day files to data series using ingest_dayfile.pl script\n";
}
elsif ($src eq "moc")
{
  print "--->WARNING:Not used for source equal to <$src>. Note using movedf.pl for moc dayfile processing\n";
  exit();
}
elsif ($src eq "egsefm")
{
  print "--->WARNING:Not used for source equal to <$src>. Note need to update this script to handle egsefm dayfiles\n";
  exit();
}

#reopen log
open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
#Check if there and then delete all dayfiles that where ingested in dayfile data series 
open(DELFILE, "$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
@all_del_file_lines="";
while (<DELFILE>)
{
  $_ =~ s/\n//g;
  push(@all_del_file_lines, $_) ;
  if ($dflg eq "1") {print LF "validated file saved to drms.deleting file from file system directory:$_\n";}
  $log=`rm $doff_dir/$_`;
  if ($dflg eq "1") {print LF "$log\n";}
  $hkt_filecount =1;
}#while
if($hkt_filecount > 0)
{
  print LF "--->Deleting dayfiles that where successfully processed into data series\n";
}

close DELFILE;
#set MF file to blank work was completed
open(DELFILE, ">$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_HSB_DELETE_FILE_LIST file: $!\n";
close DELFILE;
print LF "--->exiting script getdf.pl\n";
print LF `date`;
close LF;
##############################################################################
# check_arguments()                                                          #
##############################################################################
sub check_agruments()
{
  #source index in hk_dayfile data series
  $src=$ARGV[0]; 
  #check arguments
  if ($#ARGV != 0 )
  {
    print "Usage: perl getdf.pl  <dayfile source>\nwhere source is either:hsb,moc,egsefm.\n";
    exit;
  }
  elsif ("-h" eq substr($ARGV[0],0,2) )
  {
     print "Usage: perl getdf.pl  <dayfile source>\nwhere source is either:hsb,moc,egsefm.\n";
     exit;
  }
  elsif("moc" eq substr($ARGV[0],0,3) or "hsb" eq substr($ARGV[0],0,3) or "egsefm" eq substr($ARGV[0],0,6))
  {
    #print "okay\n";
  }
  else
  {
     print "ERROR: Entered incorrect source name for data series. Use moc,hsb,\n";
     print "Usage: perl getdf.pl  <dayfile source>\nwhere source is either:hsb,moc,egsefm.\n";
     exit;
  }
}
