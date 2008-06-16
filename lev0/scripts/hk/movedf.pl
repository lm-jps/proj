#!/bin/perl
##############################################################################
# Name:        movedf.pl  - Move or get day files for processing then remove #
# Description: Get dayfiles from moc server drop directory. Move files over  #
#              to /tmp21/production/lev0/. Then call ingest_dayfile.pl       #
#              Check if successfully loaded dayfile. If so remove dayfile    #
#              from file system.                                             #
# Execution:   movedf.pl                                                     #
##############################################################################
# set Environment Variables
# for test pickup dayfile location:
$pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/home1/jsoc/sdo/lzp/MOCFiles/moc/lzp/2008_111";
# for production use:
#$pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/home1/jsoc/sdo/lzp/MOCFiles/moc/lzp/";

# for 60d test dropoff - tested with this initially
my $doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp21/production/lev0_60d/hk_moc_dayfile";
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

# pick up dayfiles and xml files there
@list_hkt_files=`find $pup_dir | grep \.hkt\$`;
@list_xml_files=`find $pup_dir | grep \.hkt.xml\$`;
#for production use - gather files there starting at june 1, 2008 to 2029
#@list_hkt_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\$`;
#@list_xml_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\.xml\$`;

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
 $log=`cp  $hkt /tmp21/production/lev0_60d/hk_moc_dayfile`;
}
foreach $xml (@list_xml_files)
{
 $xml =~ s/\n//g;
 if ($dflg eq "1") {print LF " Move XML is <$xml>\n";}
 $log=`cp  $xml /tmp21/production/lev0_60d/hk_moc_dayfile`;
}
print LF "--->Moved df and xml files to :$doff_dir\n";
close LF;

#ingest all dayfiles and xml files
#$log=`perl ingest_dayfile.pl apidlist=./df_apid_list_day_file_moc dsnlist=./df_apid_ds_list_for_moc src=moc`
$log=`perl ingest_dayfile.pl apidlist=./df_apid_list_day_file_moc-CARLS dsnlist=./df_apid_ds_list_for_moc-CARLS src=moc`;
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

