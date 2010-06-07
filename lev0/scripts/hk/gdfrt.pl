#!/usr/bin/perl
##############################################################################
# Name:        gdfrt.pl - Get  dayfiles from real time files, break up into  #
#                         minute files and send to decode_dayfile for writing#
#                         keyword values to DRMS data series.                #
# Description: Used to gather Level 0 SDO formatted dayfiles from low speed  #
#              bus using today date. Then using dd unix command breaks a     #
#              minute of data and creates a minute dayfile. Then the minute  #
#              dayfile is sent to decode_dayfile executable. decode_dayfile  #
#              executable decode packets in dayfile and writes values to     #
#              DRMS data series. Then the minute dayfile is deleted.         #
#              Currently reads in all data in dayfile for first minute. Then #
#              after sleeping 60 seconds reads next new set of packets for   #
#              next minute, etc. Set the debug flag to 0,1 or 2 to control   #
#              the amount of debug to see($dflg).                            #
#              To be done in future maybe:                                   #
#              Add flexibility to write report of processing of              #
#              keywords when running decode_dayfile with -p flag.            #
# Execution:   (1)The run options are show in help listing:                  #
#                   gdfdrms  -h                                              #
#              (2)The execution options are shown using "gdfdrms  -h" or see #
#                 below the execution format and an examle.                  #
#                   gdfrt.pl apid=<decimal apid value> pkt_size=<value>      #
#                             Where value is packet size in packet + 7       #
#                   Example:  gdfrt.pl apid=129  pkt_size=128                #
#                   Example:  gdfrt.pl apid=29  pkt_size=136                 #
# Limitation:  Setup required environment variables at top of file based on  #
#              where running. decode_dayfile setup is needed since this      #
#              script executes decode_dayfile. Script works only using real  #
#              time dayfile file with LMSAL file formats for source equal to #
#              rtmon. Here are examples input dayfile file formats.          #
#              -For src=rtmon for HMI ISP dayfile: 200804019.0x01d           #
#              -For src=rtmon for AIA ISP dayfile: 200804019.0x027           #
#              -For src=rtmon for SDO ASD dayfile: 200804019.0x081           #
#              Make sure to get pkt_size value correct or this script will   #
#              fail to pass minute files of data correctly to decode_dayfile #
#              and cause missing records in HK By APID data series           #
# Author:      Carl                                                          #
# Date:        Created April, 7, 2009                                        #
##############################################################################
# main program                                                               #
##############################################################################

#(1)common setting for all environmnents
$ENV{'SUMSERVER'}="j1.Stanford.edu";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$exec_dir=$ENV{'DF_EXEC_PATH'}="$hm/cvs/JSOC/bin/linux_x86_64";
$script_dir=$ENV{'HK_SCRIPT_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'HK_DF_RT_LOGFILE'}="log-gdfrt-apid-";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$exec_dir";
#username and user id used to check if need to restart process
$user=$ENV{'USER'};
$uid =`id -u`;
$uid =~ s/\n//g; #regular exp rm cr
# set up where to put backup logs written monthly 
$logs_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk/logs";

#(2)set common email arguments
$from_email="\"JSOC OPS\" \<jsoc_ops\@sun.Stanford.EDU\>";
$to_email="jsoc_ops\@sun.stanford.edu";
$subject_email_failed_decode_dayfile="JSOC:ERROR:HK Level 0 Real Time Processing failed to decode minute dayfiles.";

#(3) set choices for input dayfile directory
# input dayfile directory location with subdirectories with apid names(0x081,0x029,etc).
$ENV{'HK_DDF_RT_DIRECTORY'}="/hmisdp-mon/log/packets";

#(4) set choices for minute file directory
# minute dayfile directory location - this is temp location- since all files are deleted  
$minfiles=$ENV{'HK_DDF_MINUTE_FILE_DIRECTORY'}="/tmp22/production/lev0/hk_minute_dayfiles";

#(5) set debug flag 
$dflg=$ENV{'DF_GDFRT_DEBUG'}=2;#use 0 to display min debug to standard out, 1 to send log to file or use 2 to send max debug log to display

#(6) check for any arguments passed in command
&check_arguments();

#(7) setup log file 
$apid= substr($ARGV[0],5);#get apid for log file
$logfile=`echo $ENV{'HK_DF_RT_LOGFILE'}$apid`;
$dir_logfile=`echo $script_dir/$ENV{'HK_DF_RT_LOGFILE'}$apid`;
$logfile =~ s/\n//g; #regular exp rm cr
$dir_logfile =~ s/\n//g; #regular exp rm cr

if($dir_logfile eq "")
{
   #use default name
   $logfile="$script_dir/log-gdfrt-default";#send log here when logfile not set
}

#(8) open log file
open(LF,">>$dir_logfile") || die  "ERROR in gdfrt.pl:(0):Can't Open log file <$dir_logfile>: $!\n; exit;";
if($dflg == 1) {print LF `date`;}

#(9)check if need to restart
&check_if_need_to_restart($script_dir,"gdfrt.pl",$ARGV[0],$ARGV[1], $user,$uid);

#(10) get list of apids to do 
$av=&get_apid_to_do();
if($dflg == 2) {print LF "DEBUG:MESSAGE:gdfrt: apid value is $av\n";}


#(11)loop forever thru today and next days 
while (1)
{
  #(12) check if need to move log to logs directory every month when the 1st-UTC is detected. gzip backup log too.
  &check_log($logs_dir,$logfile, $script_dir);

  #(14)get list of dates to do 
  $date_to_do=get_date_to_do();
 
  ### LOG ###
  if($dflg == 2) {print LF "DEBUG:MESSAGE:gdfrt:calling get-rt-file function with apid=$av date=$date_to_do\n";}

  #(15)go through day files in directory looking for all files for given apid and date range
  $rtfile=&get_rt_files( $date_to_do, $av);

  ### LOG ###
  if($dflg == 0) {print "-->Listing of rt file to do: <$rtfile> based on date <$date_to_do>.\n";}
  if($dflg == 1) {print LF "-->(3)Listing of rt file to do: <$rtfile> based on date  <$date_to_do>.\n";}
  if($dflg == 2) {print LF "DEBUG:MESSAGE:gdfrt:Listing of rt file to do: $rtfile\n";}

  #(16) get dayfile to read packets and write to minute files
  if ($rtfile eq "")
  {
    ### LOG ###
    if($dflg == 2) {print LF "DEBUG:MESSAGE:gdfrt:got no files to process for apid:$av . Sleep 1 minute and check again\n";}

    # sleep 1 minute and go to top of while loop and check again
    sleep 60;
    next;
  }

  #(17) get today date to handle reading next day's dayfile
  $c_date=&get_todays_date();

  #(18)get packet size which is done by just getting argument passed for pkt_size=<value>
  $packet_size= &get_packet_size();
  
  #(19)read as many packets that are in dayfile and then read next new set of packets every 1 minute for rest of day
  &read_n_packets($rtfile,$c_date, $packet_size);

  #(20) continue while loop forever reading one dayfile at a time for each of the next days to come.
}
 
#(21)close logfile
print LF `date`;
close LF;



#############################################################################
# subroutine check arguments and set flags                                  #
#############################################################################
sub check_arguments()
{
  my($temp);
  $help_flg= "0";
  $date_range_flg="0";
  
  if ($#ARGV < 0)
  {
    $help_flg="1";
  }
  if ($#ARGV >= 0 && $#ARGV <= 1)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    elsif (substr($ARGV[0],0,5) eq "apid=" )
    {
      #use apid following -a for apid
      $temp= substr($ARGV[0],5);
      if($dflg == 2) { print LF "DEBUG:MESSAGE:check_arguments:apid is $temp\n";}

      # use packet size to determine slices of data to use for dd command later
      if (substr($ARGV[1],0,9) eq "pkt_size=")
      {
        $temp= substr($ARGV[1],9);
        if($dflg == 2) { print LF "DEBUG:MESSAGE:check_arguments:pkt size is $temp\n";}
      }
      else
      {
        print  LF "gdfrt.pl:ERROR: Did not use pkt_size argument. Exiting script\n";
        print  "gdfrtf.pl:1:ERROR: Did not use pkt_size argument. Exiting script\n";
        exit;
      }
    }
  }
  else
  {
    print  LF "gdfrt.pl:ERROR: Entered incorrect arguments. Exiting script\n";
    print  "\ngdfrtf.pl:3:ERROR: Entered incorrect arguments. Exiting script\n";
    print "NOTE:Usage:   gdfrt.pl  apid=< APID value in decimal number format>  pkt_size=<packet size in packet + 7>\n";
    print "NOTE:Example: gdfrt.pl apid=129  pkt_size=128\n\n";
    exit;
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
  print "(1a)Get dayfile real time data using todays date and apid value to determine which dayfile to process.\n\n";
  print "         gdfrt.pl  apid=< APID value in decimal number format>  pkt_size=<packet size in packet + 7>\n\n";
  print "         Example: gdfrt.pl apid=129  pkt_size=128\n";
  print "         Example: gdfrt.pl apid=29   pkt_size=136 \n\n";
  print "(1b)Get Help Information:\n\n";
  print "         gdfrt.pl -h  or  gdfrt.pl -help\n\n";
  print "(2) Requires that hk level0  data series(i.e.sdo.lev0_asd_0003,hmi.lev0_isp_0021,etc.) be already created.\n";
  print "(3) Required that the environment variable are setup correctly at top of script(HK_DDF_RT_DIRECTORY, etc).\n";
  print "(4) Required that the packet_size argument be correct value for apids packet.\n";
  print "(5)****Limitation****:\n";
  print "(5a)Works only on rtmon dayfile formats only(20090407.0x081)\n";
  print "(5b)Enter arguments in specified order when running script is required.\n";
  print "(5c)It is required that a correct setup for decode_dayfile execution is done since this scripts uses decode_dayfile.\n";
  print "(5d)It is required that a correct value be use for packet size argument. Get packet size in data file or hk config\n";
  print "    file then add 7 to get value to use(for sdo asd packet, 121 + 7 = 128).\n\n\n";
  exit;
}

#############################################################################
# subroutine get_apid_to_do: gets list of apid to create maps files for     #
#############################################################################
sub get_apid_to_do
{
  # create list of files to process using apid values in file
  #push decimal character value in this format dddd.Example 0001
  $apid= substr($ARGV[0],5);
  if($dflg == 1) {print  LF "-->(1)Doing decimal formatted apids: $apid\n";}
  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_apid_to_do:apid doing = $apid\n";}
  return($apid);
}




##########################################################################
# subroutine get_todays_date()                                           #
##########################################################################
sub get_todays_date()
{
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  if($dflg == 2) {print  LF "DEBUG:MESSAGE:get_todays_date: year:<$year> month:<$month> day:<$dayOfMonth>\n";}
  #create todays date format and push on list of dates to do
  $new_date= sprintf("%4d%02.2d%02.2d", $year,$month,$dayOfMonth);
  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_todays_date: new date is ::: $new_date :::\n";}
  return ( $new_date);
}


#############################################################################
# subroutine get_date_to_do: gets date to do                                #
#############################################################################
sub get_date_to_do()
{
  my($new_date);
  # initialize new_date to null
  $new_date=();

  # get UTC today
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();

  ### LOG ###
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  if($dflg == 0) {print  "-->Today's Date is $month-$dayOfMonth-$year\n";}
  if($dflg == 1) {print LF "-->(2)Today's Date is $month-$dayOfMonth-$year\n";}
  if($dflg == 2) {print LF "DEBUG:MESSAGE:gdfrt: year this $year month is $month day is $dayOfMonth\n";}

  #create todays date format and set new_date
  $new_date= sprintf("%4d%02.2d%02.2d", $year,$month,$dayOfMonth);

  ### LOG ###
  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_date_to_do: Date is $new_date\n";}

  return ($new_date);
}



###############################################################################
# subroutine get_list_from_drms: get list of filenames to decode keywords for #
###############################################################################
sub check_lastday_in_month($$$)
{
     my($loop_yr, $loop_mo, $loop_day);
     $loop_yr= $_[0];
     $loop_mo= $_[1];
     $loop_day= $_[2];

     if( (($loop_mo == 4 || $loop_mo == 6 || $loop_mo == 9 || $loop_mo == 11 ) && $loop_day == 30) || (($loop_mo == 1 || $loop_mo == 3 || $loop_mo == 5 || $loop_mo == 7 || $loop_mo == 8|| $loop_mo == 10 || $loop_mo == 12) && $loop_day == 31))
     {
        return 1;
     }
     elsif ( $loop_mo == 2 && (($loop_day == 29 && $loop_yr == 2008) || ($loop_day == 29 && $loop_yr == 2012) || ($loop_day == 29 && $loop_yr == 2016) || ($loop_day == 29 && $loop_yr == 2020) || ($loop_day == 29 && $loop_yr == 2024)))
     {
        return 1;
     }
     elsif ( $loop_mo == 2 && $loop_day == 28 && $loop_yr != 2008 &&  $loop_yr != 2012 && $loop_yr != 2016  && $loop_yr != 2020 && $loop_yr != 2024 )
     {
        return 1;
     }
     else
     {
        return 0;
     }
}



######################################################################################
# subroutine get_rt_files(): get real time dayfiles and send back file               #
######################################################################################
sub get_rt_files($,$)
{
  # setup my variables
  my($dr, $dn_fn, @apiddir_to_do_list, @dir_list, $apid_dir,$apid_file,@apid_filelist,$fullpathname,$file_in_date_range);
  $dr=$_[0];#date 
  $apid_value=$_[1];#apid value

  #initialize file in date range to null
  @apiddir_to_do_list=(); #added 5-1-2009
  @dir_list=(); #added 5-1-2009
  @apid_filelist=(); #added 5-1-2009
  $apid_dir=""; #added 5-1-2009
  $apid_file=""; #added 5-1-2009
  $fullpathname=""; #added 5-1-2009
  $dn_fn=""; #added 5-1-2009
  $file_in_date_range="";#added 5-1-2009

  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files: Passed argument values:$dr:$apid_value:\n";}

  # get directory name to find apid directories 
  $dn = sprintf("%s",  $ENV{'HK_DDF_RT_DIRECTORY'});
  opendir(DIR_RT, $dn) || die "ERROR in gdfrt.pl:(3):Can't open rt directory <$dn>. Check setting of HK_DDF_RT_DIRECTORY environment variable:$!\n"; #open subdirectory
  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files: Directory to find dayfile directories:$dn:\n";}

  # read in all directories  and put in list 
  @dir_list = readdir(DIR_RT); #get a list of directory contents

  ### LOG ###
  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files: Looking at these directory for apid: $apid_value: Got directory:@dir_list:\n";}

  #check which apid directories to process
  foreach $dir  (@dir_list)
  {
    if( hex substr($dir,2,5) ne  $apid_value)
    {
       next;
    }
    else
    {
      if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files: Found apid directory:$dir:\n";}
      push(@apiddir_to_do_list, $dir);
      last;
    }
  }

  # check each apid directory for dayfiles
  foreach $apid_dir  (@apiddir_to_do_list)
  {
    #found directory for this APID
    #check day files to do in each directory
    $dn_fn = sprintf("%s/%s",  $ENV{'HK_DDF_RT_DIRECTORY'}, $apid_dir);
    if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files:Checking this directory to find this apid's dayfiles:$dn_fn:\n";}

    #open directory looking for this apids dayfiles
    opendir(DIR_RT_APID, $dn_fn) || die "ERROR in gdfrt.pl:(4):Can't open apid directory <$dn_fn>:$!\n"; #open subdirectory

    #read all files for apid directory
    @apid_filelist = readdir(DIR_RT_APID); #get a list of directory contents

    #loop thru each file looking for file with today's date
    foreach $apid_file  (@apid_filelist)
    {
      # Filter out 20090327.0x0081x 20090320.0x0081x and . and .. files
      if (substr($apid_file,15,1) eq "x" || $apid_file eq ".." || $apid_file eq ".")
      { 
         next;
      }#outer if
      else
      {
        # check first 8 characters of date (20090320.0x0081) is within range
        if (int substr($apid_file,0,8) == int $dr )
        {
          if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files:Found dayfile with good date and apid value :$dr:$apid_file\n";}
          $fullpathname = sprintf("%s/%s",$dn_fn, $apid_file);
          $file_in_date_range= $fullpathname;
          last;
        }#inner if
        else
        {
          ;#skip not in date range
          #if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files:skip no files found in date range\n";}
        }#inner else
      }#outer else
    }# inner foreach file in apid directory
    last if(int substr($apid_file,0,8) == int $dr );#if found file break out of outer loop too
  }# outer foreach
  if($dflg == 2) {print LF "DEBUG:MESSAGE:get_rt_files:Returned this file for processing :$file_in_date_range\n";}
  return ($file_in_date_range );
}



#####################################################################################################
# subroutine read_n_packets(): read n number of packet per minute and send to decode dayfile execute#
#####################################################################################################
sub read_n_packets($,$,$)
{
   # set my variables
   my($file,$apid, $f, $d, $source,$curr_date, $x,$y,$skip,$psize);
   # passed arguement
   $file=$_[0];
   $curr_date=$_[1];
   $psize=$_[2];#packet_size;

   #local variable initialized
   $finished_curr_df_flg=0;
   my $lasted_trigger_mail=0; #inital start up value to at least trigger email at start up if needed
   my $trigger_timer=60;#minutes value, trigger every hour

   ### LOG ###
   if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets:Called read_n_packets with arguments file:$file and todays date:$curr_date\n";}

   # set up delay to exiting loop when today's date changes 
   $wait_clock_trigger=10;
   $wait_clock=0;

   #get directory and filename and apid value
   ($d, $f) = $file=~ m/(.*\/)(.*)$/;

   #get apid from filename
   $apid= hex substr($f,11,4);

   ### LOG ###
   if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:Apid used for display is $apid\n";}
   if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:Packet size used for dd command is $psize\n";}

   # open dayfile to process to get current size of file
   open(DAYFILE, "$d/$f") || die "ERROR in gdfrt.pl:(5):Can't Open dayfile to read:$d/$f $!\n";

   #get file size
   $size = (stat(DAYFILE))[7];
   if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:File:$d/$f: has a file size :$size:\n";}

   #pkt count
   $current_pkts_in_file= $size/$psize;
   $pkt_count= $current_pkts_in_file;

   ### LOG ###
   if($dflg == 1) {print LF "-->(5)Packet count to be used in dd command is $pkt_count\n";}
   if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:Packet count to be used in dd command is $pkt_count\n";}

   #skip initialize to 0
   $skip=0;

   # set up parameter to do if compare in while loop
   # max is total number of packets in current file
   $max= $current_pkts_in_file;
   # x is the number of packet did so far
   $x=0;
   # use y to increment every minute if execute dd command
   $y=0; 

   while ( 1 )
   {
      # check if time to process in packets in file
      if ( $x < $max)
      {
        $wait_clock=0;

        # set pkt count base on how many pkts are left to do in file.
        $pkt_count = $current_pkts_in_file;
  
        ### LOG ###
        if($dflg == 0) {printf( "-->Executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s UTC\n",get_current_time());}
        if($dflg == 1) {printf( LF "-->(4)Executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s UTC\n",get_current_time());}
        if($dflg == 2) {printf( LF "DEBUG:MESSAGE:read_n_packets:Executing dd command at %s UTC\n", get_current_time());}

        # Execute dd command to create minue files.
        $log=`dd if=$d/$f of=$minfiles/$f-$y-minute bs=$psize count=$pkt_count skip=$skip`;

        ### LOG ###
        # check for errors when running dd command
        if ($? ne 0)
        {
          #print " $? \n";
          print "ERROR: Got error executing dd command. Check packet size passed as argument is correct. Packet size used is $psize. ";
          print "When packet size is divided by file size to process. Should get even packet count number. Getting packet count to do as $pkt_count . Exiting script.\n";

          print LF "ERROR: Got error executing dd command. Check packet size passed as argument is correct. Packet size used is $psize. ";
          print LF "When packet size is divided by file size to process. Should get even packet count number. Getting packet count to do as $pkt_count . Exiting script.\n";
          exit;
        }
        else
        {
          ### LOG ###
          if($dflg == 0) {printf("-->Completed executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s UTC\n",get_current_time());}
          if($dflg == 1) {printf( LF "-->(5)Completed executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s UTC\n",get_current_time());}
          if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:Complete executing dd command at %s UTC\n",get_current_time());}
        }


        ### LOG ###
        if($dflg == 0) {printf("-->Executing decode_dayfile command writing <$pkt_count> packets to drms series at %s UTC\n",get_current_time());}
        if($dflg == 1) {printf(LF "-->(6)Executing decode_dayfile command writing <$pkt_count> packets to drms series at %s UTC\n",get_current_time());}
        if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:Executing decode_dayfile command writing <$pkt_count> packets to drms series at %s UTC\n",get_current_time());}

         
        # Execute decode dayfile on minute file
        $source="rtmon"; #hardcoded since this script only used for real time dayfiles
        # decode and write to drms and create minute report in minfiles directory- COMMENT OUT and used instead command below
        #$log=`$exec_dir/decode_dayfile -p src=$source in=$minfiles/$f-$y-minute > $minfiles/Report-$f-$y-minute`;
        # decode and write to drms only - create NO report!
        $log=`$exec_dir/decode_dayfile src=$source in=$minfiles/$f-$y-minute   2>&1`;

        ### LOG ###
        # check status returned when executing decode_dayfile
        $log =~ m/(ERROR)/ ; #regular exp - look for field
        if($log eq "" )
        {
          $failed_processing_flag=0;#set to 0 if had no errors
          if($dflg == 0) {printf("-->Completed executing decode_dayfile with status PASSED at %s UTC\n",get_current_time());}
          if($dflg == 1) {printf(LF "-->(7)Completed executing decode_dayfile with status PASSED at %s UTC\n",get_current_time());}
          if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:Completed executing decode_dayfile with status PASSED at %s UTC\n",get_current_time());}
        }
        else
        {
          $failed_processing_flag=1;#set to 1 if had errors returned from decode_dayfile
          if($dflg == 0) {print  "ERROR:Status FAILED when executing decode_dayfile.DRMS may be down.\nLOG returned:\n$log:\n";}
          if($dflg == 1) {print LF "ERROR:Status FAILED when executing decode_dayfile.DRMS may be down\nLOG returned:\n$log:\n";}
          if($dflg == 2) {print LF "ERROR:read_n_packets:Status FAILED when executing decode_dayfile. DRMS Database may be down.\n";}
          if($dflg == 2) {print LF "LOG returned:$log:\n";}
          if(trigger_send_email($lasted_trigger_mail,$trigger_timer))
          {
            sendEmail("$to_email", "$from_email", "$subject_email_failed_decode_dayfile","ERROR Message:\n-->Failed to decode dayfiles into keywords and write keywords to hk data series for apid <$apid>.\n-->Cause could be that DRMS is down or the drms series does not exist.\n-->When DRMS comes back up, this script will start to process dayfile were last processed data successfully.\n-->Minute Day File not decoded is:<$minfiles/$f-$y-minute>\n-->Next Email Alert will be sent in 1 hour if problem persists.\n-->Here is log from decode_dayfile executable showing error messages:\n$log");
            $lasted_trigger_mail=time();#set to time in seconds to use in trigger_send_email function to test if 60 minutes passed
          }
        }

        ### LOG ###
        if($dflg == 0) {printf("-->Executing  delete of  minute file <$f-$y-minute> at %s UTC\n",get_current_time());}
        if($dflg == 1) {printf(LF  "-->(8)Executing  delete of minute file <$f-$y-minute> at %s UTC\n",get_current_time());}
        if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:Executing  delete  of minute file <$f-$y-minute> at %s UTC\n",get_current_time());}

        # delete minute files that was already loaded into data series
        unlink("$minfiles/$f-$y-minute");

        ### LOG ###
        if($dflg == 0) {printf("-->Completed executing  delete of  minute file <$f-$y-minute> at %s UTC\n",get_current_time());}
        if($dflg == 1) {printf(LF  "-->(9)Completed executing  delete of minute file <$f-$y-minute> at %s UTC\n",get_current_time());}
        if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:Completed executing  delete  of minute file <$f-$y-minute> at %s UTC\n",get_current_time());}
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets:sleeping 60\n\n";}

        # sleep 1 minute and then repeat above
        sleep 60;

        #Reopen and get current file size 
        close (DAYFILE);

        # reopen file and get latest file size
        open(DAYFILE, "$d/$f") || die "ERROR in gdfdrms.pl:(6):Can't Open day file to read:$d/$f  $!\n";
        $size = (stat(DAYFILE))[7];

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:File:$d/$f: has a file size :$size:\n";}

        # get total number of packets in file
        $total_pkts_in_file= $size/$psize;

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets:(if)After reopen, total packets in file: $total_pkts_in_file\n";}

        # Total Number of packets processed from  minute files so far are the  current_pkts_in_file amount just did and x amount did before
        if($failed_processing_flag == 0)
        {
          # Total Number of packets processed from  minute files so far are the  current_pkts_in_file amount just did and x amount did before
          $x=$current_pkts_in_file + $x ; #number of packets did incrementally so far

          # Update skip value to amount skipped last time with the amount did just now
          $skip= $pkt_count + $skip ; #watch pkt-count value and current-pkts-in-file are same most the time but not all the time.
        }
        else
        {
          #if decode_dayfile failed processing minute of dayfile data then skip incrementing x and skip values
	  if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets:Since decode_dayfile failed. Will retry number of packets failed to do.\n";}

	  if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets:Preparing to rerun with decode_dayfile when is working starting at packet $x.\n";}
        }
 
        # Reset max value in loop to total pkts in file
        $max=$total_pkts_in_file;

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:Next dd command will skip  $skip\n";}

        # calculate number of packet left to do now
        $current_pkts_in_file = $total_pkts_in_file - $x;

        # increment y
        $y++;
      
      }# if (x < max) 
      else
      {
        # when there is no new packets to do in current file, then sleep 1 minute or check for new days dayfile
        # get today date and time
        $td_date=get_todays_date();

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet:(else)If today's date:$td_date > $curr_date than previous day and next dayfile there and did old day processing one last time then break from loop.\n";}


        # check if next days dayfile there do old dayfile one last time then start doing next days dayfile
        if( $td_date > $curr_date && $finished_curr_df_flg == 0 && check_for_nextdays_df($td_date,$av))
        {
          #do current dayfile one last time
          $finished_curr_df_flg= 1;
          if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets: Next day detected and next dayfile detected. Do current dayfile one last time.\n";}
        }
        elsif($td_date > $curr_date && $finished_curr_df_flg == 1 && check_for_nextdays_df($td_date,$av))
        {
          ### LOG ###
          if($dflg == 0) {print  "-->Next day detected and next days dayfile detected. Today date is now <$td_date> and previous day was <$curr_date>. Try getting new dayfile to process.\n\n";}
          if($dflg == 1) {print LF "-->(11)Next day detected and next days dayfile detected. Today date is now <$td_date> and previous day was <$curr_date>. Try getting new dayfile to process.\n\n";}
          if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets: Next day detected and next days dayfile detected. Break from loop and reread new dayfile for next day.\n\n";}

          #do next days dayfile 
          $finished_curr_df_flg= 0;

          # reset current date
          $curr_date=$td_date;

          # break from this loop to check for dayfiles for next day
          last;
           
        }

        #if at end of day and at end of current day-check finish-current-dayfile-flag if need to not sleep and just do old dayfile one last time
        if ($finished_curr_df_flg == 1)
        {
          #### LOG ###
          if($dflg == 0) {printf("-->No new packet in dayfile, don't write to minute file, skip sleeping 1 minute and check old dayfile one last time for new packets at %s UTC\n\n",get_current_time());}
          if($dflg == 1) {printf(LF  "-->(10)No new packet in dayfile, don't write minute file, skip sleeping 1 minute and will check old dayfile one last last time for new packets at %s UTC\n\n",get_current_time());}
          if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:No new packets in dayfile, don't write minute file, skip sleep 1 minute and check old dayfile one last time for new packets at %s UTC\n\n",get_current_time());}

        }
        else
        {
          #### LOG ###
          if($dflg == 0) {printf("-->No new packet in dayfile sleeping 1 minute and will check again for new packets in dayfile at %s UTC\n\n",get_current_time());}
          if($dflg == 1) {printf(LF  "-->(10)No new packet in dayfile sleeping 1 minute and will check again for new packets in dayfile at %s UTC\n\n",get_current_time());}
          if($dflg == 2) {printf(LF "DEBUG:MESSAGE:read_n_packets:No new packets in dayfile, skip writing minute file, sleep 60 seconds and wait for more at %s UTC\n\n",get_current_time());}

          # sleep 1 minute and then repeat above
          sleep 60;
        }

        #Reopen file and get current file size
        close (DAYFILE);
        open(DAYFILE, "$d/$f") || die "ERROR in gdfdrms.pl:(7):Can't Open dayfile to read file:$d/$f  $!\n";
        $size = (stat(DAYFILE))[7];

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet: (else)File:$d/$f: has a file size :$size:\n";}

        # get total number of packets in file
        $total_pkts_in_file= $size/$psize;

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packets: (else)After reopen, total packets in file: $total_pkts_in_file\n";}

        # Number of packets created minute files for so far is current_pkts_in_file amount just did and x amount did before
        $x=$current_pkts_in_file + $x ; #number of packets did incrementally

        # Reset max value in loop to total pkts in file
        $max=$total_pkts_in_file;

        ### LOG ###
        if($dflg == 2) {print LF "DEBUG:MESSAGE:read_n_packet: (else)during next dd command will skip  $skip\n";}

        # calculate number of packet left to do now
        $current_pkts_in_file = $total_pkts_in_file - $x;
         
      }#else if x = max
   } # end while
   close (DAYFILE);
}


#####################################################################################################
# subroutine get_packet_size(): gets packet size argument passed                                    #
#####################################################################################################
sub get_packet_size($)
{
  my($p_size);
  $p_size= substr($ARGV[1],9);#get pkt size
  return ($p_size);

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
  $new_date= sprintf("%4d.%02.2d.%02.2d_%02.2d:%02.2d:%02.2d",$year,$month,$dayOfMonth,$hour,$minute,$second);
  return ( $new_date);
}
##########################################################################
# subroutine check_for_nextdays_df():check for next days dayfile         #
##########################################################################
sub check_for_nextdays_df($$)
{
  my($tdate,$apid,$f);
  # passed arguement
  $tdate=$_[0];#today's date using
  $apid=$_[1];#application id

  ### LOG ###
  print LF "DEBUG:MESSAGE:check_for_nextdays_df:year-month-day passed is $tdate apid passed is $apid\n";

  $f=&get_rt_files( $tdate, $apid);
  if($f eq "")
  {
     print LF "DEBUG:MESSAGE:check_for_nextdays_df:no file found :<$f>\n";
     return (0);
  }
  else
  {
     print LF "DEBUG:MESSAGE:check_for_nextdays_df:file found :<$f>\n";
     return (1);
  }
      
}



########################################################################################
# check_if_need_to_restart(): checks if need to restart. Use with cronjob required     #
########################################################################################
sub check_if_need_to_restart($,$,$,$,$,$)
{
  my $script_dir=$_[0];
  my $scriptname= $_[1];
  my $arg0_apid=$_[2];
  my $arg1_pktsize=$_[3];
  my $username=$_[4];#username for his current process
  my $userid=$_[5];#userid for this current process-use if run in background for username

  # get PID for this process
  my $pid=`echo $$`;
  $pid =~ s/\n//g;

  #Log#
  if($dflg == 2) {print LF "\n------------Start Check for Restart-----------------------------\n";}
  if($dflg == 2) {my $d_date=`date`;printf(LF "DEBUG:MESSAGE:check_if_need_to_restart: Date for check to start of script is %s",$d_date);}
  if($dflg == 2) {print LF "DEBUG:MESSAGE:check_if_need_to_restart: New script or cronjob script execution of gdfrt.pl where pid<$pid>, user<$username>, userid:<$userid> and apid<$arg0_apid>\n";}

  #search for process using script name
  my $cmd = "ps -ef | grep $scriptname | grep -v grep | grep -v sh";
  my $result = `$cmd`;
  my @lines = split("\\n", $result);

  # go through lines with name and check for this script name with apid and pkt size arguments
  foreach(@lines) 
  {
    # split up line to get pid and username data
    my  @s_line= split / \s*/,  $_;

    #LOG#
    if($dflg == 2) {print LF "DEBUG:MESSAGE:check_if_need_to_restart: Currently running gdfrt.pl script where pid<$s_line[1]>, user/userid<$s_line[0]> and apid<$s_line[9]>\n";}

    # Check if this scriptname and args is already running for same username
    if( ($_  =~ m/$scriptname $arg0_apid $arg1_pktsize/) and ($s_line[1] ne $pid) and ($s_line[0] eq $username or $s_line[0] eq $userid))
    {
      # this script is running, so exit and don't start running another process when one is running already
      if($dflg == 0) {print  "-->Check if process running status is <running already>. Exiting from this run <$scriptname $arg0_apid $arg1_pktsize>.\n";}
      if($dflg == 1) {print LF "-->(0)Check if process running status is <running already>. Exiting from this run <$scriptname $arg0_apid $arg1_pktsize>.\n";}
      if($dflg == 2) {print LF "DEBUG:MESSAGE:check_if_need_to_restart: FOUND ANOTHER PROCESS RUNNING FOR THIS USER & APID. Exiting the restart of <$scriptname $arg0_apid $arg1_pktsize> since is already running\n";}
      if($dflg == 2) {print LF "------------End Check for Restart-----------------------------\n\n";}
      print "Exiting script because script is already running on this system for this user:<$username> or userid:<$userid> with these apid arguments:<$arg0_apid>\n";
      close(LF);
      exit(0);
    }
    #else 
    if($dflg == 2) {print LF "DEBUG:MESSAGE:check_if_need_to_restart: DID NOT FIND PROCESS RUNNING FOR THIS USER & APID. The only process running is this process.\n";}
  }
  if($dflg == 0) {print  "-->Check if process running status is <not running>. Start running this process <$scriptname $arg0_apid $arg1_pktsize>.\n";}
  if($dflg == 1) {print LF "-->(0)Check if process running status is <not running>. Start running this process <$scriptname $arg0_apid $arg1_pktsize>.\n";}
  if($dflg == 2) {print LF "DEBUG:MESSAGE:check_if_need_to_restart: PROCESS NOT RUNNING AS THIS USER & APID -so run this script or cronjob-script to restart  <$scriptname $arg0_apid $arg1_pktsize> process.\n";}
  if($dflg == 2) {print LF "------------End Check for Restart-----------------------------\n\n";}
  return;# Return and start running the script in forever loop
}



##########################################
# check_log: used to create monthly logs #
##########################################
sub check_log($,$,$)
{
  use Time::Local 'timelocal';
  use POSIX;

  # args passed
  my $logs_directory= $_[0];
  my $log_filename= $_[1];
  my $script_directory= $_[2];

  #local params
  my $trigger_day=1;#day to trigger save of log file to backup logs directory
  my $backup_fn;
  my $log_cp;
  my $log_gzip;

  # get current time
  #($sec,$min,$hour,$mday,$monoffset,$yearoffset,$wday,$yday,$isdst) = localtime();
  # use UTC time
  ($sec,$min,$hour,$mday,$monoffset,$yearoffset,$wday,$yday,$isdst) = gmtime();
  if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log for Today's day(UTC) =<$mday> and trigger day is =<$trigger_day>\n"};

  if($mday == $trigger_day)
  {
    #LOG#
    if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log:Found trigger day or first day of month\n"};

    # set year using year offset
    $year= $yearoffset + 1900;

    # set month using last month by not considering offset 
    $mon=$monoffset ;#check if last months log is there

    # get directory and filename for backup log using month and year values
    $backup_fn= sprintf("%s/%s-%4d%02.2d", $logs_directory,$log_filename,$year,$mon);
    $backup_fn =~ s/\n//g; #regular exp rm cr
    $backup_gz_fn= sprintf("%s.%s", $backup_fn,"gz");
    $backup_gz_fn =~ s/\n//g; #regular exp rm cr

    #LOG#
    if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log:File to copy backup log to is <$backup_fn>\n"};
    if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log:gzip file is <$backup_gz_fn>\n"};

    #LOG#
    if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log:Checking if backup file already exists. If does not, create backup file.\n"};

    #check if last months log there.
    if (-e $backup_fn or -e $backup_gz_fn)
    {
      #if there then already copied log over so just return
      #since we are into next day we don't worry about restarts since log is concatenated during last month 
      if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log:Log file from last month already copied. skip copying again.\n"};
      return;
    }
    else
    {
   
      #check if logs directory exists
      if ( -e $logs_directory)
      {
        # copy over log to logs directory
        if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log:copy <$script_directory/$log_filename> to <$backup_fn>\n"};
        $log_cp=`cp $script_directory/$log_filename  $backup_fn`;
        $log_gzip=`gzip $backup_fn`;
        

        #set log file to blank- copy was completed
        local $|=1;
        truncate( LF, 0 );
        if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log: Completed copying logfile. Log:<$log_cp>\n"};
        if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log: Completed gzip of logfile in logs directory. Log:<$log_gzip>\n"};
        if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log: Completed copying log to logs directory & setting logfile for log to blank\n\n"};
        #close LF;#keep open for new never ending log messages
      }
      else
      {
        if($dflg == 0) {print  "-->gdfrt.pl:check_log:WARNING:missing logs directory:<$logs_directory>. Create one at <$script_directory>\n"};
        if($dflg == 2) {print LF "WARNING:check_log:missing logs directory:<$logs_directory>. Create one at <$script_directory>\n"};
      }
    }
  }
  else
  {
    #else return and do nothing since not the 1th
   if($dflg == 2) {print LF "DEBUG:MESSAGE:check_log: Day is not the first or trigger day, so no need to copy log.\n"};
  }
  return;
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

######################
# trigger send Email #
######################
sub trigger_send_email($,$)
{
  my $last_mail = $_[0];
  my $trigger_minutes = $_[1];
  my $trigger_secs;
  my $current_time;

  my $trigger_secs = $trigger_minutes * 60;#get minutes into seconds
  my $current_time=time();
  if (($last_mail + $trigger_secs) > $current_time)
  {
     #don't send email
     return 0;
  }
  else
  {
     #send email
     return 1;
  }
}
