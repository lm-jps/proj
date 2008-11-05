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
# $pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/tmp21/jsoc/sdo_dev/mocprods/lzp";
# for production use:
$pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/tmp21/jsoc/sdo/mocprods/lzp";

# Process setup (2)
# for production - tested with this initially
$doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp21/production/lev0/hk_moc_dayfile";
# for lev0 cpt dropoff
#my $doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp21/production/lev0/hk_moc_dayfile";

# debug flag
my $dflg=$ENV{'DF_MOVEDF_DEBUG'}="0";

#set common email arguments
$from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
$to_email="jsoc_ops\@sun.stanford.edu";
$subject_email_no_files="JSOC:WARNING:Ingesting MOC dayfiles: status:no files loaded today";
$subject_email_not_sure="JSOC:WARNING:Ingesting MOC dayfiles: status:not sure of status";

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
@list_hkt_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]|2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\$`;
@list_xml_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]|2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\.xml\$`;

# set log file based on sdo, moc, or egsefm
if ($source eq "hsb")
{
  $logfile="log-df-hsb";
}
elsif ($source eq "moc")
{
  $logfile="log-df-moc";
}
elsif ($source eq "egsefm")
{
  $logfile="log-df-egsefm";
}

# set up where to put backup logs written monthly 
$logs_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk/logs";

# open log file
open(LF,">>$script_dir/$logfile") || die "Can't Open $script_dir/$logfile: $!\n";
print LF `date`;
print LF "--->starting script movedf.pl\n";

# move files over to /tmp21/production/lev0/hk_moc_dayfile
foreach $hkt (@list_hkt_files)
{
 $hkt =~ s/\n//g;
 if ($dflg eq "1") {print LF " Move HKT is <$hkt>\n";}
 #$log=`cp  $hkt  $doff_dir`;
 # to use during production
 $log=`/bin/mv  $hkt  $doff_dir`;
 if ($dflg eq "1") {print LF "log after mv: $log"};
}
foreach $xml (@list_xml_files)
{
 $xml =~ s/\n//g;
 if ($dflg eq "1") {print LF " Move XML is <$xml>\n";}
 #$log=`cp  $xml  $doff_dir`;
 # to use during production
 $log=`/bin/mv  $xml  $doff_dir`;
 if ($dflg eq "1") {print LF "log after mv: $log"};
}
print LF "--->Moved df and xml files to :$doff_dir\n";
$hkt_filecount=@list_hkt_files;   
$xml_filecount=@list_xml_files;    
print LF "--->hkt file count is:$hkt_filecount xml file count is:$xml_filecount\n";

# check if files are there and chmod on files
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  ##$log=`chmod 777 $doff_dir/*`;
  print LF "--->skip chmod to 777 for df and xml files in :<$doff_dir>\n";
}
close LF;
# check if dayfile there before calling ingest_dayfile.pl script
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  #Process setup (3) and (4)
  #ingest all dayfiles and xml files
  # for production
  $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl apidlist=$script_dir/df_apid_list_day_file_moc dsnlist=$script_dir/df_apid_ds_list_for_moc src=moc`;
  print LF "--->Processed df and xml files to data series. Log results:$log\n";
}

#reopen log
open(LF,">>$script_dir/$logfile") || die "Can't Open $script_dir/$logfile: $!\n";
#Check if there and then delete all dayfiles that where ingested in dayfile data series 
open(DELFILE, "$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
@all_del_file_lines="";
$delcount=0;
while (<DELFILE>)
{
   $_ =~ s/\n//g;
   push(@all_del_file_lines, $_) ;
   if ($dflg eq "1") {print LF "validated file saved to drms.deleting file from file system directory:$_\n";}
   $log=`rm $doff_dir/$_`;
   if ($dflg eq "1") {print LF "$log\n";}
   $delcount++;
}#while
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  print LF "--->Deleting df and xml files that where successfully processed into data series\n";
}
print LF "--->movedf.pl deleted <$delcount> xml and hkt files\n";

#check if number of file got equal number loaded and deleted
if ($delcount == ($hkt_filecount +   $xml_filecount) and $delcount != 0)
{
  print LF "--->status:successfully got all files loaded into DRMS\n";
}
elsif (($hkt_filecount +   $xml_filecount) == 0)
{
  print LF "--->status:warning got no files to loaded into DRMS\n";
  ##note sent email here if occurs
  sendEmail("$to_email", "$from_email", "$subject_email_no_files","Warning Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script </home/production/cvs/JSOC/proj/lev0/scripts/hk/movedf.pl> from cron job.\n");
}
else
{
  print LF "--->status:not sure of status but got delcount:$delcount hkt_filecount:$hkt_filecount xml_filecount:$xml_filecount\n";
  sendEmail("carl\@sun.stanford.edu", "$from_email", "$subject_email_not_sure", "Warning Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script </home/production/cvs/JSOC/proj/lev0/scripts/hk/movedf.pl> from cron job.\n");
}

close DELFILE;
#set MF file to blank work was completed
open(DELFILE, ">$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
close DELFILE;
print LF "--->exiting script movedf.pl\n";
print LF `date`;
close LF;

# move log to logs directory every month
&check_log();


##########################################
# check_log: used to create monthly logs #
##########################################
sub check_log()
{
  use Time::Local 'timelocal';
  use POSIX;

  # get current time
  ($sec,$min,$hour,$mday,$monoffset,$yearoffset,$wday,$yday,$isdst) = localtime();

  # set year using year offset
  $year= $yearoffset + 1900;

  # set month using month offset
  $mon=$monoffset + 1;

  #if date is the end of month  then backup log
  #since will run once a day should get one log
  if (($mday == 31 && ($mon == 1 || $mon == 3 || $mon == 5 || $mon == 7 || $mon == 8 || $mon == 10 || $mon == 12)) || ($mday == 30 && ($mon == 4 || $mon == 6 || $mon == 9 ||  $mon == 11)) || ($mon == 2 && mday == 28))
  {
    #check if logs directory exists
    if ( -e $logs_dir)
    {
      $lm=`cp $script_dir/$logfile $logs_dir/$logfile-$mon-$year`;
      #set log file to blank - copy was completed
      open(LF, ">$script_dir/$logfile") || die "Can't Open $script_dir/$logfile file: $!\n";
      close LF;
    }
    else
    {
      print "WARNING:movedf.pl:missing logs directory:<$logs_dir>. Create one at <$script_dir>\n";
    }
  }
}

#####################
# sendEmail         #
#####################
# Simple Email Function
# ($to, $from, $subject, $message)
sub sendEmail
{
my ($to, $from, $subject, $message) = @_;
my $sendmail = '/usr/lib/sendmail';
open(MAIL, "|$sendmail -oi -t");
print MAIL "From: $from\n";
print MAIL "To: $to\n";
print MAIL "Subject: $subject\n\n";
print MAIL "$message\n";
close(MAIL);
}
