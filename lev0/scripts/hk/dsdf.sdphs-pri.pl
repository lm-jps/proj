#!/usr/bin/perl
#Hacked from dsdf.pl to get packets from sdphs-pri instead of hmisdp-mon for rtmon flag
##############################################################################
# Name:        dsdf.sdphs-pri.pl  - Decode and Save Day Files                          #
# Description: Get dayfiles from moc server dropoff directory or rtmon       #
#              directory. Move files over to /tmp22/production/lev0/.        #
#              Then if dayfiles apid is in list of apids to decode then      #
#              execute decode_dayfile on dayfile. Then if dayfile apid is in #
#              list to process min,max,mean and standard deviation keywords  #
#              then execute load_m3sd on dayfile. If successful decode, then #
#              execute ingest_dayfile.pl with merged value equal to 1. If    #
#              not on list to decode, then execute ingest_dayfile.pl with    #
#              merged value equal to 0. If on decode list but failed to      #
#              decode the dayfile, then execute ingest_dayfile.pl with       #
#              merged value equal to -1. Always check if successfully loaded #
#              dayfile using ingest_dayfile. If so remove dayfile from file  #
#              system.                                                       #
# Execution:   dsdf.sdphs-pri.pl < rtmon | moc >                                       #
#              dsdf.sdphs-pri.pl -h                                                    #
# Examples:    dsdf.sdphs-pri.pl moc                                                   #
# Examples:    dsdf.sdphs-pri.pl rtmon                                                 #
# Setup Process: (1)Setup path to pickup directory where Art's scripts puts  #
#                   moc dayfiles from moc product server. The path is also   #
#                   for the location of rtmon dayfiles and xml files on the  #
#                   sdphs-pri server.                                       #
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
#                (6)Setup file with apids to decode(dsdf_apid_list_decode_moc#
#                   or dsdf_apid_list_decode_rtmon).                         #
#                (7)Setup file with apids to process min,max,mean, and       #
#                   standard deviation values for(dsdf_apid_list_m3sd_moc    #
#                   or dsdf_apid_list_m3sd_rtmon).                           #
#                (8)Run dsdf.sdphs-pri.pl with source value so that script knows which #
#                   place to pickup and dropoff dayfiles, which file to use  #
#                   for the decode and m3sd list of apids, how to parse      #
#                   differently formatted files and were to put log          #
#                   information.                                             #
# Limitation:    (1)The move-and-ingest-dayfiles or copy-and-ingest-dayfiles #
#                   process works for dayfiles from the sdo moc product      #
#                   server or rtmon server only. This script is not used     #
#                   currently to move and  call ingest_dayfiles.pl for the   #
#                   HSB dayfiles. But can be enhanced to do hsb dayfiles.    #
#                (2)Setup in script is for production environment.           #
#                (3)Setup parameters at top of script in main of script      #
##############################################################################

#
# (1)check input argument
# check argument passed
$source=check_agruments($ARGV[0]);

# (2)setup common  variables
$ENV{'SUMSERVER'}="k1";
$ENV{'JSOC_DBUSER'}="production";
$hm="/home/jsoc/cvs/Development";
$exec_dir=$ENV{'DF_EXEC_PATH'}="$hm/JSOC/bin/linux_x86_64";
$ENV{'MAILTO'}="";
$script_dir="$hm/JSOC/proj/lev0/scripts/hk";
$log_dir="/home/jsocprod/hk/logs";
$script_lev1_dir="$hm/JSOC/proj/lev1/scripts";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$script_lev1_dir:$ENV{'DF_EXEC_PATH'}";

# (3)setup instruction filename for each apid to use when running load_m3sd
$inst_apid_17="$hm/TBL_JSOC/lev1/instruction_file/prod/hmi.leg_status.txt";
$inst_apid_19="$hm/TBL_JSOC/lev1/instruction_file/prod/hmi_thermal_300s_template.txt";
$inst_apid_21="$hm/TBL_JSOC/lev1/instruction_file/prod/hmi.iss_status.txt";
$inst_apid_38="$hm/TBL_JSOC/lev1/instruction_file/prod/aia.iss_3_4_status.txt";
$inst_apid_40="$hm/TBL_JSOC/lev1/instruction_file/prod/aia.iss_1_2_status.txt";
$inst_apid_44="$hm/TBL_JSOC/lev1/instruction_file/prod/aia_thermal_300s_template.txt";

# (4)setup pickup and dropoff directories and email variables based on source value
if ($source eq "moc")
{
  # moc pickup directory setup:
  $pup_dir=$ENV{'DF_PICKUP_MOC_FILES'}="/tmp28/jsocprod/sdo/mocprods/lzp";

  # moc dropoff directory setup
  $doff_dir=$ENV{'DF_DROPOFF_MOC_FILES'}="/tmp28/jsocprod/lev0/hk_moc_dayfile";

  # moc common email arguments
  $from_email="\"JSOC AND\" \<jsoc_ans\@sun.Stanford.EDU\>";
  $to_email="jsoc_ans\@sun.stanford.edu";
  $subject_email_no_files="JSOC:WARNING:Ingesting MOC dayfiles: status:No files loaded today";
  $subject_email_not_sure="JSOC:WARNING:Ingesting MOC dayfiles: status:Possible error ingesting dayfile.";
  $subject_email_no_decode="JSOC:ERROR:Decoding MOC dayfiles into housekeeping data series.";
  $subject_email_no_m3sd="JSOC:ERROR:M3SD Processing of MOC dayfiles into housekeeping data series.";
  $subject_email_no_decode_m3sd="JSOC:ERROR:Decoding and M3SD Processing of MOC dayfiles into housekeeping data series.";
  $subject_email_not_equal_count="JSOC:ERROR:Ingesting MOC dayfiles: Did not ingest dayfiles because some hkt files don't have corresponding xml files";
  $subject_email_warn_m3sd="JSOC:WARNING:M3SD Processing of all six MOC dayfiles into housekeeping data series was NOT done.";
  $subject_email_warn_decode="JSOC:WARNING:Decode Dayfile Processing of MOC dayfiles into asd housekeeping data series was NOT done.";
  $subject_email_warn_m3sd_decode="JSOC:WARNING:M3SD Processing and Decode Dayfile Processing of all MOC dayfiles NOT done.";

}
elsif ($source eq "rtmon")
{
  # rtmon pickup directory:
  $pup_dir=$ENV{'DF_PICKUP_RTMON_FILES'}="/sdphs-pri/log/packets";

  # rtmon dropoff directory setup 
  $doff_dir=$ENV{'DF_DROPOFF_RTMON_FILES'}="/tmp28/jsocprod/lev0/hk_rtmon_dayfile";

  # rtmon common email variables
  $from_email="\"JSOC ANS\" \<jsoc_ans\@sun.Stanford.EDU\>";
  $to_email="jsoc_ans\@sun.stanford.edu";
  $subject_email_no_files="JSOC:WARNING:Ingesting RTMON dayfiles: status:No files loaded today";
  $subject_email_not_sure="JSOC:WARNING:Ingesting RTMON dayfiles: status:Possible error ingesting dayfile.";
  $subject_email_no_decode="JSOC:ERROR:Decoding RTMON dayfiles into housekeeping data series.";
  $subject_email_not_equal_count="JSOC:ERROR:Ingesting RTMON dayfiles: Did not ingest dayfiles since some hkt files don't have corresponding xml files.";
  $subject_email_warn_m3sd="JSOC:WARNING:M3SD Processing of six RTMOC dayfiles into housekeeping data series was NOT done.";
  $subject_email_warn_decode="JSOC:WARNING:Decode Dayfile Processing of RTMON dayfiles into asd housekeeping data series was NOT done.";
  $subject_email_warn_m3sd_decode="JSOC:WARNING:M3SD Processing and Decode Dayfile Processing of all RTMON dayfiles NOT done.";
}
elsif ( $source == "hsb" ) 
{
  print "ERROR: unknown source value. Exting script dsdf.sdphs-pri.pl!\n";
  print "ERROR: note hsb is not implemented in dsdf.sdphs-pri.pl yet!\n";
  exit;
}
else 
{
   print "ERROR: unknown source value. Exting script dsdf.sdphs-pri.pl!\n";
   exit;
}

# (5)debug flag where 0=minimum log message 1=debug mode or max log messages
$dflg=$ENV{'SAVE_DF_DEBUG'}=0;


# (6)local variables for capturing processing errors
# failed list of dayfiles processed using do_decode_dayfile()
my @failed_to_decode_list=();#initialize to empty, this is global variable
# failed list of dayfiles processed using do_m3sd_dayfile()
my @failed_to_m3sd_list=();#initialize to empty, this is global variable
# reference to lists of error so can pass in functions and set
my $ref_failed_to_decode_list = \@failed_to_decode_list; #ref setup
my $ref_failed_to_m3sd_list = \@failed_to_m3sd_list;

# (7)pick up dayfiles and xml files from the pickup directory
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
  #print  "dsdf.sdphs-pri.pl:LIST HKT: @list_hkt_files\n";
  #print  "dsdf.sdphs-pri.pl:LIST XML: @list_xml_files\n";
}
else
{
  print "ERROR: Unknown source value:<$source>. Exiting script dsdf.sdphs-pri.pl!\n";
  exit;
}


# (8)set log file based on rtmon,hsb, moc, or egsefm
if ($source eq "rtmon")
{
  $logfile="log-df-rtmon";
}
elsif ($source eq "moc")
{
  $logfile="log-df-moc";
}
else
{
  $logfile="log-df-default";
}

# (9)set up where to put backup logs written monthly
$logs_dir="$log_dir/old";

# (10)open log file
open(LF,">>$log_dir/$logfile") || die "Can't Open $log_dir/$logfile: $!\n";
print LF `/bin/date -u`;
print LF "--->Starting script dsdf.sdphs-pri.pl\n";
# send message to users running script at command line so they know something is working
print "...running dsdf.sdphs-pri.pl script\n...please wait to finish\n...view status on run by doing: tail -f $log_dir\/$logfile\n";

# (11)"move" files over to /tmp28/jsocprod/lev0/hk_moc_dayfile for moc
# or "copy" files to /tmp28/jsocprod/lev0/hk_rtmon_dayfile for rtmon
if ($source eq "moc")
{
  foreach $hkt (@list_hkt_files)
  {
    $hkt =~ s/\n//g;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE: Move HKT is <$hkt>\n";}
    $log=`/bin/mv  $hkt  $doff_dir`;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:log after mv: $log\n"};
  }
  foreach $xml (@list_xml_files)
  {
    $xml =~ s/\n//g;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:Move XML is <$xml>\n";}
    $log=`/bin/mv  $xml  $doff_dir`;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:log after mv: $log\n"};
  }
  print LF "--->Moved df and xml files to :$doff_dir\n";
}
elsif ($source eq "rtmon")
{
  foreach $hkt (@list_hkt_files)
  {
    $hkt =~ s/\n//g;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:Copy HKT is <$hkt>\n";}
    $log=`/bin/cp  $hkt  $doff_dir`;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:log after cp: $log\n"};
  }
  foreach $xml (@list_xml_files)
  {
    $xml =~ s/\n//g;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:Copy XML is <$xml>\n";}
    $log=`/bin/cp  $xml  $doff_dir`;
    if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:log after cp: $log\n"};
  }
  print LF "--->Copied df and xml files to :$doff_dir\n";
}

# (12)get file counts for xml and dayfiles
$hkt_filecount=@list_hkt_files;   
$xml_filecount=@list_xml_files;    
print LF "--->hkt file count is:$hkt_filecount xml file count is:$xml_filecount\n";

# (14)check if files are there and if needed chmod on files
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  #$log=`chmod 777 $doff_dir/*`;
  print LF "--->skip chmod to 777 for df and xml files in :<$doff_dir>\n";
}


# (15)check if dayfiles and xml files are there before starting processing
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  #ingest all dayfiles and xml files
  # for production
  if ($source eq "moc")
  {
    # (16.1)do decode dayfile and ingest dayfile for moc dayfiles
    $lu_map_file="df_apid_ds_list_for_moc";
    ($success_decode_count,$success_m3sd_count)=do_process_dayfile($doff_dir, $lu_map_file, $source, $ref_failed_to_decode_list, $ref_failed_to_m3sd_list );
  }
  elsif ($source eq "rtmon")
  {
    # (16.2)do decode dayfile and ingest dayfile for rtmon dayfiles
    $lu_map_file="df_apid_ds_list_for_rtmon";
     ($success_decode_count,$success_m3sd_count)=do_process_dayfile($doff_dir, $lu_map_file, $source, $ref_failed_to_decode_list, $ref_failed_to_m3sd_list);
  }
}

# (17)reopen log
open(LF,">>$log_dir/$logfile") || die "Can't Open $log_dir/$logfile: $!\n";
print LF "--->Processed df and xml files to data series. Log results:$log\n";

# (18)Check if there and then delete all dayfiles that where ingested in dayfile data series 
#     note:These files gets loaded with files by ingest_dayfile.pl, if ingest was successful.
if ($source eq "moc")
{
  open(DELFILE, "$log_dir/DF_DELETE_FILE_LIST_MOC") || die "(6)Can't Open $log_dir/DF_DELETE_FILE_LIST_MOC file: $!\n";
}
elsif ($source eq "rtmon")
{
  open(DELFILE, "$log_dir/DF_DELETE_FILE_LIST_RTMON") || die "(6)Can't Open $log_dir/DF_DELETE_FILE_LIST_RTMON file: $!\n";
}
else
{
  print LF "ERROR: unknown source value when trying to open DF_DELETE_FILES_LIST!\n";
}

@all_del_file_lines=();
$delcount=0;
while (<DELFILE>)
{
   $_ =~ s/\n//g;
   push(@all_del_file_lines, $_) ;
   if ($dflg == 1) {print LF "--->DEBUG:MESSAGE:validated file saved to drms.deleting file from file system directory:$_\n";}
   $log=`rm $doff_dir/$_`;
   if ($dflg == 1) {print LF "--->DEBUG:MESSAGE: log after removed files:<$log>\n";}
   $delcount++;
}#while

# (19)Check status of run and send emails if have possible errors
if($hkt_filecount > 0 or $xml_filecount > 0)
{
  print LF "--->Deleting df and xml files that where successfully processed into data series\n";
}
print LF "--->Deleted <$delcount> xml and hkt files\n";

# (20)check if number of file got equal number loaded and deleted
if ($delcount == ($hkt_filecount +   $xml_filecount ) and $delcount != 0)
{
  print LF "--->status:successfully got all day files loaded into DRMS\n";
}
elsif (($hkt_filecount +   $xml_filecount) == 0)
{
  print LF "--->status:warning got no files to loaded into DRMS\n";
  ##note sent email here if occurs
  sendEmail("$to_email", "$from_email", "$subject_email_no_files","Warning Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script <$hm/JSOC/proj/lev0/scripts/hk/dsdf.sdphs-pri.pl> from cron job.\n");
}
elsif($hkt_filecount !=  $xml_filecount)
{
  print LF "--->status:failed to load hkt and xml file in dayfile series.\n";;
  print LF "--->delcount:$delcount hkt_filecount:$hkt_filecount xml_filecount:$xml_filecount\n";
  print LF "--->hkt and xml file counts should match therefore could not load dayfile and xml to dayfile series.\n";
  print LF "--->Check why each hkt file does not have a corresponding xml.\n";
  sendEmail("$to_email", "$from_email", "$subject_email_not_equal_count", "Error Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script <$hm/JSOC/proj/lev0/scripts/hk/dsdf.sdphs-pri.pl> from cron job.\n-->All of the  hkt files don't have corresponding xml files therefore can not ingest dayfile.\n-->Check why each hkt file does not have a corresponding xml.\n");

}
else
{
  print LF "--->status:not sure of status but got delcount:$delcount hkt_filecount:$hkt_filecount xml_filecount:$xml_filecount\n";
  print LF "--->Check if there is problem. The dayfiles to ingest into data series and delete from directory did not match the count of the dayfiles received.\n--->Possibly a problem ingesting dayfiles in series because of bad setting of SUMSERVER parameter or SUMS could be not available.\n--->Possibly not an issue which was caused by the dayfiles not being ingested on previous day(s) therefore the file count received today does not match files ingested.\n";
  sendEmail("$to_email", "$from_email", "$subject_email_not_sure", "Warning Message:\n-->Received count of <$hkt_filecount> hkt files and count of <$xml_filecount> xml files from directory <$pup_dir>.\n-->When executing script <$hm/JSOC/proj/lev0/scripts/hk/dsdf.sdphs-pri.pl> from cron job.\n-->Check if there is problem. The dayfiles to ingest into data series and delete from directory did not match the count of the dayfiles received.\n-->Possibly a problem ingesting dayfiles in series because of bad setting of SUMSERVER parameter or SUMS could be not available.\n-->Possibly not an issue which was caused by the dayfiles not being ingested on previous day(s) therefore the file count received today does not match files ingested.\n");
}

close DELFILE;

# (21)set DELFILE file to blank since work was completed
if ($source eq "moc")
{
open(DELFILE, ">$log_dir/DF_DELETE_FILE_LIST_MOC") || die "(6)Can't Open $log_dir/DF_DELETE_FILE_LIST_MOC file: $!\n";
}
elsif ($source eq "rtmon")
{
open(DELFILE, ">$log_dir/DF_DELETE_FILE_LIST_RTMON") || die "(6)Can't Open $log_dir/DF_DELETE_FILE_LIST_RTMON file: $!\n";
}


# (22)close log and delete list files
close DELFILE;

# (23)check status on decode of dayfile keywords and send separate email if get decode failure
# get number of files failed to decode
$failed_ddf_count= @$ref_failed_to_decode_list;
$failed_m3sd_count= @$ref_failed_to_m3sd_list;

# (23)check if failed to decode dayfile or/and load m3sd using decode_dayfile and load_m3sd executables
if($failed_ddf_count == 0 and $failed_m3sd_count == 0)
{
  print  LF "--->status:successfully decoded dayfiles with no errors. Failed count <$failed_ddf_count>\n";
  print  LF "--->status:successfully loaded m3sd for dayfiles with no errors. Failed count <$failed_m3sd_count>\n";
}
elsif($failed_ddf_count == 0 and $failed_m3sd_count > 0)
{
  print  LF "--->status:successfully decoded dayfiles with no errors. Failed count <$failed_ddf_count>\n";
  print  LF "--->status:error:failed to load m3sd for all required dayfiles. Number of files failed to decode is <$failed_m3sd_count>\n";
  sendEmail("$to_email", "$from_email", "$subject_email_no_m3sd","Error Message:\n-->Failed to load m3sd data into keywords into hk data series.\n-->Number of dayfiles the could not be processed for m3sd data  was <$failed_m3sd_count>.\n-->Day Files not processed are:<@$ref_failed_to_m3sd_list>");
}
elsif($failed_ddf_count > 0 and $failed_m3sd_count == 0)
{
  print  LF "--->status:successfully loaded m3sd for dayfiles with no errors. Failed count <$failed_m3sd_count>\n";
  print  LF "--->status:error:failed to decode for all required dayfiles. Number of files failed to decode is <$failed_ddf_count>\n";
  sendEmail("$to_email", "$from_email", "$subject_email_no_decode","Error Message:\n-->Failed to decode dayfiles into keywords into hk data series.\n-->Number of dayfiles the could not be decoded was <$failed_ddf_count>.\n-->Day Files not decoded are:<@$ref_failed_to_decode_list>");
}
else
{
  print  LF "--->status:error:failed to decode all required dayfiles. Number of files failed to decode is <$failed_ddf_count>\n";
  print  LF "--->status:error:failed to load m3sd for all required dayfiles. Number of files failed to decode is <$failed_m3sd_count>\n";
  sendEmail("$to_email", "$from_email", "$subject_email_no_decode_m3sd","Error Message #1:\n-->Failed to decode dayfiles into keywords into hk data series.\n-->Number of dayfiles the could not be decoded was <$failed_ddf_count>.\n-->Day Files not decoded are:<@$ref_failed_to_decode_list>\nError Message #2:\n-->Failed process dayfiles into min,max,mean and standard deviation(m3sd) keywords into hk data series.\n-->Number of dayfiles the could not process m3sd keywords was <$failed_m3sd_count>.\n-->Day Files not decoded are:<@$ref_failed_to_m3sd_list>\n");
}

# (24)check if processed all six files dayfiles using load_m3sd and check if processed 1 file using decode_dayfile.
if($source eq "moc")
{
  #setup expected count of files to process and then check if processed successfully
  $expected_count_m3sd=6; #later can automate better by getting count via file
  $expected_count_decode=1; #later can automate better by getting count via file

  # check count not expected and check did not send other emails above already
  if (($success_m3sd_count ne $expected_count_m3sd) and ($failed_m3sd_count == 0) and ($success_decode_count ne $expected_count_decode) and ($failed_ddf_count == 0))
  {
     sendEmail("$to_email", "$from_email", "$subject_email_warn_m3sd_decode","Warning Message(1):\n-->Did not process the expected  <$expected_count_m3sd> dayfiles using load_m3sd executable.\n-->Only processed <$success_m3sd_count> dayfile(s)\n-->Check if the dayfiles for apid 17, 19, 21, 38, 40 and 44 have been retrieved and ingested in < hmi | aia >.hk_dayfile series .\n-->If dayfiles exist, get dayfiles and rerun with load_m3sd executable.\n\nWarning Message(2):\n-->Did not process the expected  <$expected_count_decode> dayfiles using decode_dayfile executable.\n-->Note sdo.lev0_asd_0004 was not updated today.\n-->Only processed <$success_decode_count> dayfile(s).\n-->Check if the dayfiles for apid 129 has been retrieved and ingested in sdo.hk_dayfile series .\n-->If dayfile exists, get dayfiles and rerun with decode_dayfile executable.");
    print LF "--->status:WARNING:did not process expected <$expected_count_m3sd> dayfiles using load_m3sd executable.\n";
    print LF "--->status:WARNING:did not process expected <$expected_count_decode> dayfiles using decode_dayfile executable. sdo.lev0_asd_0004 was not updated today.\n";
  } 
  elsif (($success_m3sd_count ne $expected_count_m3sd) and ($failed_m3sd_count == 0))
  {
     sendEmail("$to_email", "$from_email", "$subject_email_warn_m3sd","Warning Message:\n-->Did not process the expected  <$expected_count_m3sd> dayfiles using load_m3sd executable.\n-->Only processed <$success_m3sd_count> dayfile(s)\n-->Check if the dayfiles for apid 17, 19, 21, 38, 40 and 44 have been retrieved and ingested in < hmi | aia >.hk_dayfile series .\n-->If dayfiles exist, get dayfiles and rerun with load_m3sd executable.");
    print LF "--->status:WARNING:did not process expected <$expected_count_m3sd> dayfiles using load_m3sd executable.\n";
  }
  # check count not expected and check did not send other emails above already
  elsif (($success_decode_count ne $expected_count_decode) and ($failed_ddf_count == 0))
  {
     sendEmail("$to_email", "$from_email", "$subject_email_warn_decode","Warning Message:\n-->Did not process the expected  <$expected_count_decode> dayfile(s) using decode_dayfile executable.\n-->Note sdo.lev0_asd_0004 was not updated today.\n-->Only processed <$success_decode_count> dayfile(s).\n-->Check if the dayfiles for apid 129 has been retrieved and ingested in sdo.hk_dayfile series .\n-->If dayfile exists, get dayfiles and rerun with decode_dayfile executable.");
    print LF "--->status:WARNING:did not process expected <$expected_count_decode> dayfiles using decode_dayfile executable. sdo.lev0_asd_0004 was not updated today.\n";
  }
  
}#end if source is moc
elsif ($source eq "rtmon") 
{
  #setup expected count of files to process and then check if processed successfully
  $expected_count_m3sd=6; #later can automate better by getting count via file
  $expected_count_decode=0; #later can automate better by getting count via file

  # check count not expected and check did not send other emails above already
  if (($success_m3sd_count ne $expected_count_m3sd) and ($failed_m3sd_count == 0) and ($success_decode_count ne $expected_count_decode) and ($failed_ddf_count == 0))
  {
     sendEmail("$to_email", "$from_email", "$subject_email_warn_m3sd_decode","Warning Message(1):\n-->Did not process the expected  <$expected_count_m3sd> dayfiles using load_m3sd executable.\n-->Only processed <$success_m3sd_count> dayfile(s)\n-->Check if the dayfiles for apid 17, 19, 21, 38, 40 and 44 have been retrieved and ingested in < hmi | aia >.hk_dayfile series .\n-->If dayfiles exist, get dayfiles and rerun with load_m3sd executable.\n\nWarning Message(2):\n-->Did not process the expected  <$expected_count_decode> dayfiles using decode_dayfile executable.\n-->Note sdo.lev0_asd_0004 was not updated today.\n-->Only processed <$success_decode_count> dayfile(s).\n-->Check if the dayfiles for apid 129 has been retrieved and ingested in sdo.hk_dayfile series .\n-->If dayfile exists, get dayfiles and rerun with decode_dayfile executable.");
    print LF "--->status:WARNING:did not process expected <$expected_count_m3sd> dayfiles using load_m3sd executable.\n";
    print LF "--->status:WARNING:did not process expected <$expected_count_decode> dayfiles using decode_dayfile executable. sdo.lev0_asd_0004 was not updated today.\n";
  } 
  elsif (($success_m3sd_count ne $expected_count_m3sd) and ($failed_m3sd_count == 0))
  {
     sendEmail("$to_email", "$from_email", "$subject_email_warn_m3sd","Warning Message:\n-->Did not process the expected  <$expected_count_m3sd> dayfiles using load_m3sd executable.\n-->Only processed <$success_m3sd_count> dayfile(s)\n-->Check if the dayfiles for apid 17, 19, 21, 38, 40 and 44 have been retrieved and ingested in < hmi | aia >.hk_dayfile series .\n-->If dayfiles exist, get dayfiles and rerun with load_m3sd executable.");
    print LF "--->status:WARNING:did not process expected <$expected_count_m3sd> dayfiles using load_m3sd executable.\n";
  }
  # check count not expected and check did not send other emails above already
  elsif (($success_decode_count ne $expected_count_decode) and ($failed_ddf_count == 0))
  {
     sendEmail("$to_email", "$from_email", "$subject_email_warn_decode","Warning Message:\n-->Did not process the expected  <$expected_count_decode> dayfile(s) using decode_dayfile executable.\n-->Note sdo.lev0_asd_0004 was not updated today.\n-->Only processed <$success_decode_count> dayfile(s).\n-->Check if the dayfiles for apid 129 has been retrieved and ingested in sdo.hk_dayfile series .\n-->If dayfile exists, get dayfiles and rerun with decode_dayfile executable.");
    print LF "--->status:WARNING:did not process expected <$expected_count_decode> dayfiles using decode_dayfile executable. sdo.lev0_asd_0004 was not updated today.\n";
  }
}

# (25)close log 
print LF "--->Exiting script dsdf.sdphs-pri.pl\n";
print LF `date -u`;
close LF;

# (26) check if need to move log to logs directory. this is done every month
&check_log();

# (27)display message to users running script at command line instead of as cronjob
print "...completed running dsdf.sdphs-pri.pl.\n";
# End of main



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
      $lm=`cp $log_dir/$logfile $logs_dir/$logfile-$mon-$year`;
      #set log file to blank - copy was completed
      open(LF, ">$log_dir/$logfile") || die "Can't Open $log_dir/$logfile file: $!\n";
      close LF;
    }
    else
    {
      print LF "WARNING:movedf.pl:missing logs directory:<$logs_dir>. Create one at <$log_dir>\n";
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
     print "Usage: perl dsdf.sdphs-pri.pl  <dayfiles from moc or rtmon directory on filesystem >\nwhere can be  moc or rtmon.\n->Note currently script handles rtmon or moc as source values.\n->Note that the setup for pickup and dropoff directories is for production user to run on j0.\n->Note that you can read header of script for more help information.\n->Note create list of apids to decode using files dsdf_apid_list_decode_rtmon(not implemented yet) and/or dsdf_apid_list_decode_moc.\n->Note create mapping file to identify series to write dayfiles to using df_apid_ds_list_for_rtmon or df_apid_ds_list_for_moc files.\n->Note that this routine uses decode_dayfile executable and ingest_dayfile.pl script therefore for more information view help information on these routines by using -h flag\n->Note create list of apids to run load_m3sd executable on using file dsdf_apid_list_m3sd_moc.\n";
     print "Limitations:\n";
     print "(1) Only used to process and save dayfile(s) with source equal to moc or rtmon.\n";
     print "(2) Only implemented to do decode_dayfile and load_m3sd on MOC dayfiles.\n";
     print "(3) Only implemented to do decode_dayfile and load_m3sd on RTMON dayfiles.\n";
     print "(4) Only setup to run on production-at-j0.\n";
     print "(5) Only setup for load_m3sd to use instruction files on production and for apid 17, 19, 21, 38, 40 and 44 .\n";
     print "(6) Only setup to check processing status is done on files for one MOC day file for decode_dayfile and six MOC day files for load_m3sd\n";
     print "(7) Only setup to check processing status is done on files for zero RTMON day files for decode_dayfile and six RTMON day files for load_m3sd\n";
     print "(8) Need to update expected variable setting when add apid to do processing decode_dayfile and load_m3sd to config files(dsdf_apid_list_m3sd_moc, df_apid_ds_list_for_rtmon, etc).\n";
     print "(9) To run in another environment check and redo variable setting in steps(1)-(9) settings\n";

     #exist
     exit;
  }
  elsif("moc" eq substr($src,0,3) or  "rtmon" eq substr($src,0,5))
  {
    #print "okay\n";
  }
  else
  {
     print "ERROR: Entered incorrect source name for data series. Use moc or rtmon only!\n";
     print "Usage: perl dsdf.sdphs-pri.pl  <dayfile source from moc or rtmon directory on filesystem >\nwhere source is either: moc or rtmon.\nNote currently script handles rtmon or moc as source dayfile values but script can be updated for handling src=hsb.\n";
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

#########################################################################
#@apid_decode_list = get_apid_decode_list();
#########################################################################
sub get_apid_decode_list($)
{
  my $input=$_[0];
  my @dlist="";
  my $fn_moc="$script_dir/dsdf_apid_list_decode_moc";
  my $fn_rtmoc="$script_dir/dsdf_apid_list_decode_rtmon";
  my $fn_hsb="$script_dir/dsdf_apid_list_decode_hsb";
  if($input eq "moc")
  {
    open(FILE_DECODE_APID_LIST, "$fn_moc") || die "Can't Open: <$fn_moc> file: $!\n";
    while (<FILE_DECODE_APID_LIST>)
    {
        push(@dlist, $_); 
    }
    close(FILE_APID_LIST);
  }
  else
  {
    if($dflg == 1) {print LF "--->DEBUG:MESSAGE:get_apid_decode_list:got not valid value <$input>. Exiting.\n";}
  }
  if($dflg == 1) {print LF "--->DEBUG:MESSAGE:get_apid_decode_list:returning decode list:@dlist\n";}

  return(@dlist);
}



#########################################################################
#@apid_m3sd_list = get_apid_m3sd_list();
#########################################################################
sub get_apid_m3sd_list($)
{
  my $input=$_[0];
  my @mlist=();
  my $fn_moc="$script_dir/dsdf_apid_list_m3sd_moc";
  my $fn_rtmon="$script_dir/dsdf_apid_list_m3sd_rtmon";
  my $fn_hsb="$script_dir/dsdf_apid_list_m3sd_hsb";

  # current only use moc option for src value
  if($input eq "moc")
  {
    open(FILE_MEAN_APID_LIST, "$fn_moc") || die "Can't Open: <$fn_moc> file: $!\n";
    while (<FILE_MEAN_APID_LIST>)
    {
        push(@mlist, $_); 
    }
    close(FILE_MEAN_LIST);
  }
  elsif($input eq "rtmon")
  {
    open(FILE_MEAN_APID_LIST, "$fn_rtmon") || die "Can't Open: <$fn_rtmon> file: $!\n";
    while (<FILE_MEAN_APID_LIST>)
    {
        push(@mlist, $_); 
    }
    close(FILE_MEAN_LIST);
  }
  else
  {
    if($dflg == 1) {print LF "--->DEBUG:MESSAGE:get_apid_mean_list:got not valid value <$input>. Exiting.\n";}
  }
  if($dflg == 1) {print LF "--->DEBUG:MESSAGE:get_apid_mean_list:returning decode list:@mlist\n";}

  return(@mlist);
}



#########################################################################
#%decode_status_hk = do_process_dayfile(@apid_decode_list);
#########################################################################
sub do_process_dayfile($,$,$,$,$)
{
  #local variables
  my $found_decode_apid, $found_m3sd_apid;
  my(@filelist,$file,$lu,$decode_apid,$src,$arg1,$arg2,$arg3,$arg4,$arg5,$arg6);
  $successful_m3sd_files_processed=0;
  $successful_decode_files_processed=0;

  #passed arguments dropoff directory for dayfiles, source and decode apid list 
  my ($indir,$lu,$src,$ref_failed_decode_list,$ref_failed_m3sd_list)=@_;

  #get dayfiles in doff_dir
  #open dir and read all files there
  opendir(DIR_INFILE, $indir) || die ":(1)Can't open:$!\n"; #open subdirectory

  # read in all map files and put in list 
  @filelist = readdir(DIR_INFILE); #get a list of directory contents

  # LOG #
  if ($dflg == 1)
  {
    print LF "-->DEBUG:MESSAGE:do_process_dayfile:0:indir is <$indir>\n";
    print LF "-->DEBUG:MESSGE:do_process_dayfile:0:$lu is <$lu>\n";
    print LF "-->DEBUG:MESSAGE:do_process_dayfile:0:src is <$src>\n";
    print LF "-->DEBUG:MESSAGE:do_process_dayfile:0:decode apid list: <@dlist>\n";
    print LF "-->DEBUG:MESSAGE:do_process_dayfile:0:filelist: <@filelist>\n";
  }

  # loop through each file and skip doing decode dayfile on xml files
  foreach $file  (@filelist)
  {
    if( $src eq "moc" && $file =~ m/(xml)/ )
    {
      if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:1:got xml file <$file>-skip!\n";}
      next;
    }
    if ( $src eq "rtmon" && $file =~ m/(x$)/ )
    {
      if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:1.5:got x file <$file>-skip!\n";}
      next;
    }

    if( $file eq "." | $file eq "..")
    {
     if ($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:2:got non-dayfile <$file>-skip!\n";}
      next;
    }
   
    # call do_decode_dayfile to process dayfiles for decode-list of apids
    # set found flag to 0 meaning meaning set merged flag to 0 because apid not on decode list
    $found_decode_apid=0;
    $found_decode_apid=do_decode_dayfile($indir,$lu,$src,$file,$ref_failed_decode_list);
    if($found_decode_apid eq 1)
    {
      $successful_decode_files_processed++;
    }
    if ($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:3.1:got MERGE value <$found_decode_apid>\n";}

    # call do_m3sd_dayfile to process dayfiles for m3sd-list of apids
    $found_m3sd_apid=1;#not used to set merged flag, since decode_dayfile results are used only
    $found_decode_m3sd=do_m3sd_dayfile($indir,$lu,$src,$file, $ref_failed_m3sd_list);
    if($found_decode_m3sd eq 1)
    {
      $successful_m3sd_files_processed++;
    }
    if ($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:3.2:got ret value when ran load_m3sd <$found_m3sd_apid>\n"};

    # ingest dayfile with merged value based on found_decode_apid flag
    if ($found_decode_apid == 0)
    {
       $merged=0;
       &do_ingest_dayfile($file,$src,$merged,$lu);
    }
    elsif($found_decode_apid == 1)
    {
       $merged=1;
       &do_ingest_dayfile($file,$src,$merged,$lu);
    }
    elsif($found_decode_apid == -1)
    {
       $merged=-1;
       &do_ingest_dayfile($file,$src,$merged,$lu);
    }
  }#foreach file
 
 
  # close directory
  closedir (DIR_INFILE);
  return ($successful_decode_files_processed,$successful_m3sd_files_processed);
}



##########################################################################
# subroutine get_current_time()                                          #
##########################################################################
sub get_current_time()
{
  my($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings, $year, $month, $new_date);
  #($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  #create todays date and time format 
  $new_date= sprintf("%4d.%02.2d.%02.2d_%02.2d:%02.2d:%02.2d UTC",$year,$month,$dayOfMonth,$hour,$minute,$second);
  return ( $new_date);
}

#############################################################################
# subroutine do_ingest_dayfile                                              #
#############################################################################
sub do_ingest_dayfile($,$,$,$)
{
  my $f = $_[0];
  my $src = $_[1];
  my $merge = $_[2];
  my $apid_ds_lu = $_[3];

  # get apid and start/end date from dayfile filename
  if ($src eq "moc")
  {
    #get apid value
    $arg1=substr($f,0,4);
    if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:1:apid is <$arg1)>\n";}

    #translate 2010_044(year and doy) to yyyymmdd for start and end times--OPEN
    $arg2 = get_df_date(substr($f,5,4), substr($f,10,3));
    $arg2=~ s/\.//g;
  }
  elsif ($src eq "rtmon")
  {
    #get apid value
    $arg1= substr($f,11,4);
    if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:1:apid is <$arg1)>\n";}

    #parse file 20100329.0x0081 date
    $arg2 =  substr($f,0,8);
  }

  # set end time to start
  $arg3=$arg2;
  if($dflg == 1) 
  {
    print LF "-->DEBUG:MESSAGE:do_ingest_dayfile;2:startdate is <$arg2>\n";
    print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:3:enddate is <$arg3>\n";
  }

  #use series lookup file - use the full path to file!
  $arg4="$script_dir/$apid_ds_lu"; 
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:4.1u_map:<$arg4>\n";}

  #set source
  $arg5=$src;
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:5:source:<$arg5>\n";}

  #set merged value
  $arg6=$merge;
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:6:merged::<$arg6>\n";}

  #example exec
   if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_ingest_dayfile:7:Running ingest_dayfile.pl apid=$arg1 start=$arg2 end=$arg3 dsnlist=$arg4 src=$arg5 merged=$arg6\n";}  

  #close to log message for ingest_dayfile
  close(LF);

  $log=`$script_dir/ingest_dayfile.pl apid=$arg1 start=$arg2 end=$arg3 dsnlist=$arg4 src=$arg5 merged=$arg6`;

  #reopen log appending more messages
  open(LF,">>$log_dir/$logfile") || die "Can't Open $log_dir/$logfile: $!\n";

}



#############################################################################
# subroutine get_df_date                                                    #
#############################################################################
sub get_df_date($,$)
{
  use Time::Local 'timelocal';
  use POSIX;

  # set day in year to day1
  #$day1=$dyear ;
  $day1=$_[1] ;
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:get_df_date:day is <$day1>\n"};

  # set month to 0 to use day of year to get time
  $mon1=0;

  # set year from filename
  #$year1=$year - 1900;
  $year1=$_[0] - 1900;
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:get_df_date:year - 1900 is <$year1>\n"};

  # call mktime
  $unixtime=mktime('','','',$day1,$mon1,$year1,'','');
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:get_df_date:unixtime is <$unixtime>\n"};

  #get broke up time
  ($sec,$min,$hour,$mday,$monoffset,$yearoffset,$wday,$yday,$isdst) = localtime($unixtime);

  # set year using year offset
  $year= $yearoffset + 1900;

  # set month using month offset
  $mon=$monoffset + 1;

  # create DATE for index for dayfile
  $date_for_dayfile=sprintf("%-4.4d.%-02.2d.%-02.2d",$year,$mon,$mday);
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:get_df_date:date for dayfile is <$date_for_dayfile>\n"};

  # get time in seconds
  $tm_seconds=timelocal($sec,$min,$hour,$mday,$monoffset,$year);
  if($dflg == 1) {print LF "-->DEBUG:MESSAGE:get_df_date:timelocal returned seconds are <$tm_seconds>\n"};

  #return value yyyy.mm.dd
  return ($date_for_dayfile);

}



#############################################################################
# subroutine do_decode_dayfile -call decode_dayfile executable              #
#############################################################################
sub do_decode_dayfile($,$,$,$,$)
{

  # set input arguments
  my($indir,$lu,$src,$file,$ref_failed_list)=@_;

  # local variable
  my @decode_list;
  my $found_decode_apid;

  # set found flag to 0 meaning meaning set merged flag to 0 because apid not on decode list
  $found_decode_apid=0;

  # get apidlist_to_decode_moc or apidlist_to_decode_rtmon based on $source value
  @decode_list = get_apid_decode_list($source);
  if($dflg == 1) {print LF "--->DEBUG:MESSAGE:dsdf.sdphs-pri.pl:main:apid to decode::@decode_list\n";}

  # check if file is on list to decode
  foreach $decode_apid  (@decode_list)
  {
    $decode_apid =~ s/\n//;
    if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:4:check if have match for <$decode_apid> value for file<$file>\n";}

    # if file is on list to decode then decode
    if(substr($file,0,4) eq  $decode_apid )
    {
      #pick out files in @dlist for apids to decode
      #got one to decode
      if ($dflg == 1)
      {
        print LF "-->DEBUG:MESSAGE:do_process_dayfile:5:FOUND apid to process:<$decode_apid>\n";
        print LF "-->DEBUG:MESSAGE:do_process_dayfile:6:FOUND file to process <$file>\n";
      }

      # Execute decode dayfile on day file
      # decode and write to drms only - create NO report!
      # LOG #
      if ($dflg == 1)
      {
        print LF "-->DEBUG:MESSAGE:do_process_dayfile:7:exec command <$exec_dir/decode_dayfile>\n";
        print LF "-->DEBUG:MESSAGE:do_process_dayfile:8:source:<$src>\n";
        print LF "-->DEBUG:MESSAGE:do_process_dayfile:9:in file is <$indir/$file>\n";
      }
      printf(LF "-->Start running decode_dayfile for src=<%s> and dayfile=<%s> at %s\n",$src,$file,get_current_time());

      # decode dayfile
      $log=`$exec_dir/decode_dayfile  src=$src in=$indir/$file   2>&1`; 

      # check status returned when executing decode_dayfile
      $log =~ /(ERROR)/g ; #regular exp - look for field
      if( $1 eq "ERROR" )
      {
        # if decode_dayfile exec failed set merged=-1 for ingest_dayfile.pl
        print LF "ERROR:Status failed when executing decode_dayfile. Exiting script dsdf.sdphs-pri.pl.\nLOG returned:$log:\n";
        # set found flag to -1 meaning set merged flag to -1 because on decode list but unsuccessfully decoded dayfile
        $found_decode_apid = -1;

        # push on list dayfiles failed to decode 
        push(@$ref_failed_list,$file); 

      }#end if error found in log results of decode_dayfile
      else
      {
        #check for warning
        $log =~ /(WARNING|Warning)/g;
        if ($1 eq "WARNING" or $1 eq "Warning")
        {
          # Log #
          printf(LF "-->Completed executing decode_dayfile with status PASSED at %s\n-->But got warning messages:<%s>\n",get_current_time(),$log);
          # set found flag to 1 meaning set merged flag to 1 because on decode list and successfully decoded dayfile
          $found_decode_apid = 1;
        }
        else
        {
          # Log #
          printf(LF "-->Completed executing decode_dayfile with status PASSED at %s\n",get_current_time());

          # set found flag to 1 meaning set merged flag to 1 because on decode list and successfully decoded dayfile
          $found_decode_apid = 1;
        }

      }#end of else decode_dayfile had no errors in log results

      # LOG #
      if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_process_dayfile:14:file is <$file>\n\n";}

      #found file's apid on decode-list and processed using decode_dayfile so exit for foreach thru decode list and get next file to do
      last;

    }#if in list to decode
  }#foreach file check if should run decode_dayfile on file

  # Return MERGE value setting to use to ingest dayfile
  return($found_decode_apid);

}# end do_decode_dayfile function



#############################################################################
# subroutine do_m3sd_dayfile -call load_m3sd executable                     #
#############################################################################
sub do_m3sd_dayfile($,$,$,$,$)
{

  # set input arguments
  my($indir,$lu,$src,$file,$ref_failed_list)=@_;

  # local variable
  my @m3sd_list;
  my $found_m3sd_apid;

  # set found flag to 0 meaning meaning set merged flag to 0 because apid not on m3sd list
  $found_m3sd_apid=0;

  # get apidlist_to_m3sd_moc or apidlist_to_m3sd_rtmon based on $source value
  @m3sd_list = get_apid_m3sd_list($src);
  if($dflg == 1) {print LF "--->DEBUG:MESSAGE:do_m3sd_dayfile:apid to do load_m3sd processing::@m3sd_list\n";}

  # check if file is on list to decode
  foreach $m3sd_apid  (@m3sd_list)
  {
    $m3sd_apid =~ s/\n//;
    if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:4:check if have match for <$m3sd_apid> value for file<$file>\n";}

    # if file is on list to m3sd then do m3sd processing
    if((substr($file,2,2) eq  $m3sd_apid && $src eq "moc") or (substr($file,11,4) eq  $m3sd_apid && $src eq "rtmon"))
    {
      #pick out files in @m3sd_list for apids to process using load_m3sd
      #got one to decode
      if ($dflg == 1)
      {
        print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:5:FOUND apid to process:<$m3sd_apid>\n";
        print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:6:FOUND file to process <$file>\n";
        print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:7:FOUND src <$src>\n";
      }

      # get instruction file based on apid
      $instruct_file=get_instruction_file( $src,$m3sd_apid);
      if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:-:Check FOUND instruction file to process <$instruct_file>\n"};
      if ($instruct_file eq "")
      {
         print LF "-->ERROR:do_m3sd_dayfile:-: Could not find instruction file to process <$instruct_file>. Skipping processing of file<$file>\n";
         push(@$ref_failed_list,$file); 
         $found_m3sd_apid = -1;
         last;
      }

      # Execute decode dayfile on day file
      # decode and write to drms only - create NO report!
      # LOG #
      if ($dflg == 1)
      {
        print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:7:exec command <$exec_dir/load_m3sd>\n";
        print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:8:source:<$src>\n";
        print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:9:in file is <$indir/$file>\n";
      }
      printf(LF "-->Start running load_m3sd for src=<%s> and dayfile=<%s> and instuction file=<%s> at %s\n",$src,$file,$instruct_file, get_current_time());

      # decode dayfile
      $log=`$exec_dir/load_m3sd  src=$src in=$indir/$file isf=$instruct_file  2>&1`; 
      $log =~ s/\n// ; #regular exp - look for field
      print LF "-->Results of running load_m3sd:<$log>\n";

      # check status returned when executing decode_dayfile
      $log =~ /(ERROR)/g ; #regular exp - look for field
      if( $1 eq "ERROR" )
      {
        # if decode_dayfile exec failed set merged=-1 for ingest_dayfile.pl
        print LF "ERROR:Status failed when executing load_m3sd. Exiting script dsdf.sdphs-pri.pl.\nLOG returned:$log:\n";
        # set found flag to -1 meaning set merged flag to -1 because on decode list but unsuccessfully decoded dayfile
        $found_m3sd_apid = -1;

        # push on list dayfiles failed to process using load_m3sd
        #push(@failed_to_m3sd_list,$file); 
        push(@$ref_failed_list,$file); 

      }#end if error found in log results of decode_dayfile
      else
      {
        #check for warning
        $log =~ /(WARNING|Warning)/g;
        if ($1 eq "WARNING" or $1 eq "Warning")
        {
          # Log #
          printf(LF "-->Completed executing load_m3sd with status PASSED at %s\n-->But got warning messages:<%s>\n",get_current_time(),$log);
          # set found flag to 1 meaning set merged flag to 1 because on decode list and successfully decoded dayfile
          $found_m3sd_apid = 1;
        }
        else #got no errors or warning- passed
        { 
          # Log #
          printf(LF "-->Completed executing load_m3sd with status PASSED at %s\n",get_current_time());

          # set found flag to 1 meaning set merged flag to 1 because on decode list and successfully decoded dayfile
          $found_m3sd_apid = 1;
        }

      }#end of else load_m3sd had no errors in log results

      # LOG #
      if($dflg == 1) {print LF "-->DEBUG:MESSAGE:do_m3sd_dayfile:14:file is <$file>\n\n";}

      #found file's apid on m3sd-list and processed using load_m3sd so exit for foreach thru decode list and get next file to do
      last;

    }#if in list to decode
  }#foreach file check if should run load_m3sd on file

  # Return MERGE value setting to use to ingest dayfile
  return($found_m3sd_apid);

}




#############################################################################
# subroutine get_instruction_file -gets file for 17,19,21,38,40 and 44 only!#
#############################################################################
sub get_instruction_file($,$)
{
  #local variables
  my $inst;
  $inst="";
  #arguments passed in function
  my $src=$_[0];
  my $apid=$_[1];
  
  #set instruction files for moc and rtmon dayfiles only
  if(($apid == 17 and $src eq "moc") or (hex $apid == 0x11 and $src eq "rtmon"))
  {
    $inst=$inst_apid_17;
  }
  elsif(($apid == 19 and $src eq "moc") or (hex $apid == 0x13 and $src eq "rtmon"))
  {
    $inst=$inst_apid_19;
  }
  elsif(($apid == 21 and $src eq "moc") or (hex $apid == 0x15 and $src eq "rtmon"))
  {
  $inst=$inst_apid_21;
  }
  elsif(($apid == 38 and $src eq "moc") or (hex $apid == 0x26 and $src eq "rtmon"))
  {
    $inst=$inst_apid_38;
  }
  elsif(($apid == 40 and $src eq "moc") or (hex $apid == 0x28 and $src eq "rtmon"))
  {
      $inst=$inst_apid_40;
  }
  elsif(($apid == 44 and $src eq "moc") or (hex $apid == 0x2c and $src eq "rtmon"))
  {
    $inst=$inst_apid_44;
  }
  else
  {
     print LF "-->ERROR:do_m3sd_dayfile:Cannot find instruction file for this apid <$m3sd_apid>. ";
     print LF "Add to script code for instruction file for this apid <$m3sd_apid>\n";
     #return empty string
  }

  return ($inst);

}
