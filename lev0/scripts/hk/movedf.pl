#!/usr/bin/perl
##############################################################################
# Name:        movedf.pl  - Move or get day files for processing then remove #
# Description: Get dayfiles from moc server drop directory. Move files over  #
#              to /tmp21/production/lev0/. Then call ingest_dayfile.pl       #
#              Check if successfully loaded dayfile. If so remove dayfile    #
#              from file system.                                             #
# Execution:   movedf.pl                                                     #
# Setup Process: (1)Setup path to pickup directory where Art's scripts puts  #
#                   dayfile from moc product server.                         #
#                (2)Setup path to dropoff directory where this script drops  #
#                   off dayfile using move for ingest_dayfile.pl to get files#
#                   from.                                                    #
#                (3)Setup map file which tells ingest_dayfile.pl script      #
#                   where to put each apid to which dayfile series.          #
#                   (i.e., apid 1 goes to hmi.hk_dayfile and apid 40 goes to #
#                    aia.hk_dayfile series)                                  #
#                (4)Setup apid list file which tells ingest_dayfile.pl       #
#                   script which apids to process.                           #
#                (5)Setup DF_DAYFILE_DIRECTORY environment variable use in   #
#                   ingest_dayfile.pl script to pickup files to ingest to    #
#                   value of the dropoff directory in this script($doff_dir).#
# Limitation:   The move-and-ingest-dayfiles process works for dayfiles from #
#               the sdo moc product server only. This script is not used     #
#               currently to move and  call ingest_dayfiles.pl for the HSB   #
#               dayfiles and LMSAL dayfiles.                                 #
##############################################################################
# set Environment Variables
# Process setup (1):
# for test pickup dayfile location:
# $pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/home1/jsoc/sdo/lzp/MOCFiles/moc/lzp/2008_111";
# for test art's test directory
$pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/tmp21/jsoc/sdo_dev/mocprods/lzp";
# for production use:
#$pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/home1/jsoc/sdo/lzp/MOCFiles/moc/lzp/";

# Process setup (2)
# for production - tested with this initially
my $doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp21/production/lev0/hk_moc_dayfile";
# for lev0 cpt dropoff
#my $doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp21/production/lev0/hk_moc_dayfile";

# debug flag
my $dflg=$ENV{'DF_MOVEDF_DEBUG'}="0";

#common setting for all environments
$ENV{'SUMSERVER'}="d02.Stanford.EDU";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$ENV{'DF_DRMS_EXECUTABLES'}="$hm/cvs/JSOC/bin/linux_x86_64";
my $script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
my $source="moc";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$ENV{'DF_DRMS_EXECUTABLES'}";

# pick up dayfiles and xml files there
#carl test with files from 2008_111
#@list_hkt_files=`find $pup_dir | grep \.hkt\$`;
#@list_xml_files=`find $pup_dir | grep \.hkt.xml\$`;
#for production use - gather day and xml files there starting at june 1, 2008 to 2029
@list_hkt_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\$`;
@list_xml_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\.xml\$`;

# set log file based on sdo, moc, or egsefm
if ($source eq "hsb")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-hsb";
}
elsif ($source eq "moc")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-moc";
}
elsif ($source eq "egsefm")
{
  $logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-df-egsefm";
}

# open log file
open(LF,">>$logfile") || die "Can't Open $logfile: $!\n";
print LF `date`;

# move files over to /tmp21/production/lev0/hk_moc_dayfile
foreach $hkt (@list_hkt_files)
{
 $hkt =~ s/\n//g;
 if ($dflg eq "1") {print LF " Move HKT is <$hkt>\n";}
 $log=`cp  $hkt  $doff_dir`;
 # to use during production
 #$log=`mv  $hkt  $doff_dir`;
}
foreach $xml (@list_xml_files)
{
 $xml =~ s/\n//g;
 if ($dflg eq "1") {print LF " Move XML is <$xml>\n";}
 $log=`cp  $xml  $doff_dir`;
 # to use during production
 #$log=`mv  $xml  $doff_dir`;
}
print LF "--->Moved df and xml files to :$doff_dir\n";
$log=`chmod 777 $doff_dir/*`;
print LF "--->chmod to 777 for df and xml files in :$doff_dir\n";
close LF;

#Process setup (3) and (4)
#ingest all dayfiles and xml files
# for production
print LF "--->executing </usr/bin/perl  $script_dir/ingest_dayfile.pl apidlist=$script_dir/df_apid_list_day_file_moc dsnlist=$script_dir/df_apid_ds_list_for_moc src=moc>\n";
$log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl apidlist=$script_dir/df_apid_list_day_file_moc dsnlist=$script_dir/df_apid_ds_list_for_moc src=moc`;
# for testing
#$log=`perl ingest_dayfile.pl apidlist=./df_apid_list_day_file_moc-CARLS dsnlist=./df_apid_ds_list_for_moc-CARLS src=moc`;
print LF "--->Processed df and xml files to data series\n";

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
}#while
print LF "--->Deleting df and xml files that where successfully processed into data series\n";

close DELFILE;
#set MF file to blank work was completed
open(DELFILE, ">$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
close DELFILE;
print LF `date`;
close LF;
