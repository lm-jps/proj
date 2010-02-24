#!/usr/bin/perl
##############################################################################
# Name:        savedf.pl  - Move or copy day files for processing & remove  #
# Description: Get dayfiles from moc server drop directory. Move files over  #
#              to /tmp21/production/lev0/. Then call ingest_dayfile.pl       #
#              Check if successfully loaded dayfile. If so remove dayfile    #
#              from file system.                                             #
# Execution:   savedf.pl < rtmon | moc >                                     #
# Setup Process: (1)Setup path to pickup directory where Art's scripts puts  #
#                   dayfiles from moc product server. The path is also for   #
#                   the location of rtmon dayfiles and xml files on the      #
#                   hmisdp-mon server.                                       #
#                (2)Setup path to dropoff directory where this script drops  #
#                   off dayfile using move or copy for ingest_dayfile.pl to  #
#                   get files from(i.e.,DF_DROPOFF_MOC_FILES and             #
#                   DF_PICKUP_RTMON_FILES )                                  #
#                (3)Setup map file which tells ingest_dayfile.pl script      #
#                   where to put each apid to which dayfile series.          #
#                   (i.e., apid 1 goes to hmi.hk_dayfile and apid 40 goes to #
#                    aia.hk_dayfile series)                                  #
#                (4)Setup apid list file which tells ingest_dayfile.pl       #
#                   script which apids to process.                           #
#                (5)Setup DF_DAYFILE_DIRECTORY environment variable use in   #
#                   ingest_dayfile.pl script to pickup files to ingest to    #
#                   value of the dropoff directory in this script($doff_dir).#
# Limitation:   The move-and-ingest-dayfiles or copy-and-ingest-dayfiles     #
#               process works for dayfiles from the sdo moc product server   #
#               or rtmon server only. This script is not used currently to   #
#               move and  call ingest_dayfiles.pl for the HSB dayfiles.      #
##############################################################################
#
# Process setup (1):
# check argument passed
$source=check_agruments($ARGV[0]);

# set Environment Variables
if ($source eq "moc")
{
  # Process setup (2):
  $pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/tmp21/jsoc/sdo/mocprods/lzp";

  # Process setup (3)
  $doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp21/production/lev0/hk_moc_dayfile";

  #set common email arguments(4)
  $from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
  $to_email="jsoc_ops\@sun.stanford.edu";
  $subject_email_no_files="JSOC:WARNING:Ingesting MOC dayfiles: status:No files loaded today";
  $subject_email_not_sure="JSOC:WARNING:Ingesting MOC dayfiles: status:Possible error ingesting dayfile.";

}
elsif ($source eq "rtmon")
{
  # Process setup (2):
  $pup_dir=$ENV{'DF_PICKUP_RTMON_FILES'}="/hmisdp-mon/log/packets";

  # Process setup (3)
  $doff_dir=$ENV{'DF_DROPOFF_RTMON_FILES'}="/tmp02/production/lev0/hk_rtmon_dayfile";

  #set common email arguments(4)
  $from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
  $to_email="jsoc_ops\@sun.stanford.edu";
  $subject_email_no_files="JSOC:WARNING:Ingesting RTMON dayfiles: status:No files loaded today";
  $subject_email_not_sure="JSOC:WARNING:Ingesting RTMON dayfiles: status:Possible error ingesting dayfile.";
}
else 
{
  print "ERROR: unknown source value. Exting script savedf.pl!\n";
  exit;
}

# debug flag (5)
my $dflg=$ENV{'SAVE_DF_DEBUG'}="0";

# common setting for all environments(6)
$ENV{'SUMSERVER'}="j1.Stanford.edu";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$ENV{'DF_DRMS_EXECUTABLES'}="$hm/cvs/JSOC/bin/linux_x86_64";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$ENV{'DF_DRMS_EXECUTABLES'}";

# pick up dayfiles and xml files there(7)
 
if ($source eq "moc")
{
  #for production use - gather day and xml files there starting at june 1, 2008 to 2029
  @list_hkt_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]|2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\$`;
  @list_xml_files=`find $pup_dir  | egrep '(2008_[1][5][2-9]|2008_[1][6-9][0-9]|2008_[2-3][0-9][0-9]|2009_[0-3][0-9][0-9]|20[1-2][0-9]_[0-3][0-9][0-9])' | grep \.hkt\.xml\$`;
}
elsif ($source eq "rtmon")
{
  #get current date or today's date
  $new_date=get_current_time_rtmon();

  #for production rtmon dayfiles- gather day and xml files based on 
  #today's date and the apids in apidlist file.
  $apidlist="$script_dir/df_apid_list_day_file_rtmon";

  # get list of apids in file and create search arg for dayfiles
  $search_str=create_apid_searchlist("df",$apidlist);

  # get list of dayfiles to save
  @list_hkt_files=`find $pup_dir  | egrep '($new_date)' | egrep '$search_str'`;

  # get list of apids in file and create search arg for dayfiles
  $search_str=create_apid_searchlist("xml",$apidlist);

  # get list of dayfiles to save
  @list_xml_files=`find $pup_dir  | egrep '($new_date)' | egrep '$search_str'`;
  #print  "savedf.pl:LIST HKT: @list_hkt_files\n";
  #print  "savedf.pl:LIST XML: @list_xml_files\n";
}
else
{
  print "ERROR: Unknown source value:<$source>. Exiting script savedf.pl!\n";
  exit;
}


# (8)set log file based on rtmon,hsb, moc, or egsefm
if ($source eq "rtmon")
{
  $logfile="log-df-rtmon";
}
elsif ($source eq "hsb")
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
else
{
  $logfile="log-df-default";
}

# (9)set up where to put backup logs written monthly 
$logs_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk/logs";

# (10)open log file
open(LF,">>$script_dir/$logfile") || die "Can't Open $script_dir/$logfile: $!\n";
print LF `/bin/date`;
print LF "--->starting script savedf.pl\n";

# (11)"move" files over to /tmp21/production/lev0/hk_moc_dayfile for moc
# or "copy" files to /tmp02/production/lev0/hk_rtmon_dayfile for rtmon
if ($source eq "moc")
{
  foreach $hkt (@list_hkt_files)
  {
    $hkt =~ s/\n//g;
    if ($dflg eq "1") {print LF " Move HKT is <$hkt>\n";}
    $log=`/bin/mv  $hkt  $doff_dir`;
    if ($dflg eq "1") {print LF "log after mv: $log"};
  }
  foreach $xml (@list_xml_files)
  {
    $xml =~ s/\n//g;
    if ($dflg eq "1") {print LF " Move XML is <$xml>\n";}
    $log=`/bin/mv  $xml  $doff_dir`;
    if ($dflg eq "1") {print LF "log after mv: $log"};
  }
  print LF "--->Moved df and xml files to :$doff_dir\n";
}
elsif ($source eq "rtmon")
{
  foreach $hkt (@list_hkt_files)
  {
    $hkt =~ s/\n//g;
    if ($dflg eq "1") {print LF " Copy HKT is <$hkt>\n";}
    $log=`/bin/cp  $hkt  $doff_dir`;
    if ($dflg eq "1") {print LF "log after cp: $log"};
  }
  foreach $xml (@list_xml_files)
  {
    $xml =~ s/\n//g;
    if ($dflg eq "1") {print LF " Copy XML is <$xml>\n";}
    $log=`/bin/cp  $xml  $doff_dir`;
    if ($dflg eq "1") {print LF "log after cp: $log"};
  }
  print LF "--->Copied df and xml files to :$doff_dir\n";
}

# (12)get file counts for xml and dayfiles
$hkt_filecount=@list_hkt_files;   
$xml_filecount=@list_xml_files;    
print LF "--->hkt file count is:$hkt_filecount xml file count is:$xml_filecount\n";

# (13)check if files are there and chmod on files
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  #$log=`chmod 777 $doff_dir/*`;
  print LF "--->skip chmod to 777 for df and xml files in :<$doff_dir>\n";
}
close LF;#because ingest_dayfile.pl is writing to log now. after complete reopen log below

# (14)check if dayfile there before calling ingest_dayfile.pl script
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  #ingest all dayfiles and xml files
  # for production
  if ($source eq "moc")
  {
    $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl apidlist=$script_dir/df_apid_list_day_file_moc dsnlist=$script_dir/df_apid_ds_list_for_moc src=moc`;

  }
  elsif ($source eq "rtmon")
  {
    $log=`/usr/bin/perl  $script_dir/ingest_dayfile.pl apidlist=$script_dir/df_apid_list_day_file_rtmon dsnlist=$script_dir/df_apid_ds_list_for_rtmon src=rtmon`;
  }
}

# (15)reopen log
open(LF,">>$script_dir/$logfile") || die "Can't Open $script_dir/$logfile: $!\n";
print LF "--->Processed df and xml files to data series. Log results:$log\n";

# (16)Check if there and then delete all dayfiles that where ingested in dayfile data series 
#     note:These files gets loaded with files by ingest_dayfile.pl, if ingest was successful.
if ($source eq "moc")
{
  open(DELFILE, "$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
}
elsif ($source eq "rtmon")
{
  open(DELFILE, "$script_dir/DF_DELETE_FILE_LIST_RTMON") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST_RTMON file: $!\n";
}
else
{
  print "ERROR: unknown source value when trying to open DF_DELETE_FILES_LIST!\n";
}

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


#(16)Check status of run and send emails if have possible errors
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  print LF "--->Deleting df and xml files that where successfully processed into data series\n";
}
print LF "--->savedf.pl deleted <$delcount> xml and hkt files\n";

#check if number of file got equal number loaded and deleted
##if ($delcount == ($hkt_filecount +   $xml_filecount) and $delcount != 0)
if ($delcount == ($hkt_filecount +   $xml_filecount ) and $delcount != 0)
{
  print LF "--->status:successfully got all files loaded into DRMS\n";
}
elsif (($hkt_filecount +   $xml_filecount) == 0)
{
  print LF "--->status:warning got no files to loaded into DRMS\n";
  ##note sent email here if occurs
  sendEmail("$to_email", "$from_email", "$subject_email_no_files","Warning Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script </home/production/cvs/JSOC/proj/lev0/scripts/hk/savedf.pl> from cron job.\n");
}
else
{
  print LF "--->status:not sure of status but got delcount:$delcount hkt_filecount:$hkt_filecount xml_filecount:$xml_filecount\n";
  print LF "--->Check if there is problem. The dayfiles to ingest into data series and delete from directory did not match the count of the dayfiles received.\n--->Possible problem ingesting dayfiles in series because of bad setting of SUMSERVER parameter or SUMS could be not available.\n--->Possibly not an issue which was caused by the dayfiles not being ingested on previous day(s) therefore the file count received today does not match files ingested.\n";
  sendEmail("$to_email", "$from_email", "$subject_email_not_sure", "Warning Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script </home/production/cvs/JSOC/proj/lev0/scripts/hk/savedf.pl> from cron job.\n-->Check if there is problem. The dayfiles to ingest into data series and delete from directory did not match the count of the dayfiles received.\n-->Possible problem ingesting dayfiles in series because of bad setting of SUMSERVER parameter or SUMS could be not available.\n-->Possibly not an issue which was caused by the dayfiles not being ingested on previous day(s) therefore the file count received today does not match files ingested.\n");
}

close DELFILE;

#(17)set MF file to blank work was completed
if ($source eq "moc")
{
open(DELFILE, ">$script_dir/DF_DELETE_FILE_LIST") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST file: $!\n";
}
elsif ($source eq "rtmon")
{
open(DELFILE, ">$script_dir/DF_DELETE_FILE_LIST_RTMON") || die "(6)Can't Open $script_dir/DF_DELETE_FILE_LIST_RTMON file: $!\n";
}

#(18)close log and delete list files
close DELFILE;
print LF "--->Exiting script savedf.pl\n";
print LF `date`;
close LF;

#(19) move log to logs directory every month
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

##############################################################################
# check_arguments()                                                          #
##############################################################################
sub check_agruments($)
{
  #source index in hk_dayfile data series
  $src=$_[0];
  #check arguments
  if ("-h" eq substr($src,0,2) )
  {
     print "Usage: perl savedf.pl  <dayfile source>\nwhere source is either:hsb, moc, egsefm or rtmon.\n->Note currently script handles rtmon or moc as source values.\n->Note that the setup for pickup and dropoff directories is for production user to run on n02.\n->Note that you can read header of script for more help information.\n";
     exit;
  }
  elsif("moc" eq substr($src,0,3) or "hsb" eq substr($src,0,3) or "egsefm" eq substr($src,0,6)  or "rtmon" eq substr($src,0,5))
  {
    #print "okay\n";
  }
  else
  {
     print "ERROR: Entered incorrect source name for data series. Use moc,hsb,\n";
     print "Usage: perl savedf.pl  <dayfile source>\nwhere source is either:hsb, moc, egsefm or rtmon.\nNote currently script handles rtmon or moc as source values.\n";
     exit;
  }
  return ($src);
}


##########################################################################
# subroutine get_current_time()                                          #
##########################################################################
sub get_current_time_rtmon()
{
  my($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings, $year, $month, $new_date);
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  #($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  #create todays date and time format 
  $new_date= sprintf("%4d%02.2d%02.2d",$year,$month,$dayOfMonth);
  return ( $new_date);
}


#########################################################################
# subroutine create_apid_searchlist()                                   #
#########################################################################
sub  create_apid_searchlist($,$)  
{
  my($apids,$dir_filename,@all_apids);
  $filetype=$_[0];
  $dir_filename=$_[1];
  open(FILE_APID_LIST, "$dir_filename") || die "Can't Open: <$dir_filename> file: $!\n";
  while (<FILE_APID_LIST>)
  {
    if( $filetype eq "df")
    {
      $_=~s/\n/\$|/g; #substitute in parameters needed to search list of dayfiles for today
    }
    elsif(  $filetype eq "xml")
    {
      $_=~s/\n/x\$|/g; #substitute in parameters needed to search list of dayfiles for today
    }
    push(@all_apids, $_) ;
  }
  # join array of items 
  $apids=join("", @all_apids);
  $apids=~ s/\|\z//;#remove last |
  close FILE_APID_LIST ;
  return ( $apids);
}
